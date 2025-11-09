#include <cmath>
#include <algorithm>
#include <iostream>
#include <pthread.h>
#include "solver.h"

// ===================================================
// Структуры для многопоточности
// ===================================================
struct ThreadArgs {
    int tid;
    int n;
    int k;
    double* A;
    double* b;
    int* col_perm;
    double pivot;
    double* row_k;
    int num_threads;
    pthread_barrier_t* barrier;
};

void* elimination_worker(void* arg) {
    ThreadArgs* args = (ThreadArgs*)arg;
    int n = args->n;
    int k = args->k;
    int tid = args->tid;
    int num_threads = args->num_threads;
    double* A = args->A;
    double* b = args->b;
    int* col_perm = args->col_perm;
    double pivot = args->pivot;
    double* row_k = args->row_k;

    // Каждому потоку — своя порция строк
    for (int i = k + 1 + tid; i < n; i += num_threads) {
        double* row_i = A + i * n;
        double factor = row_i[col_perm[k]] / pivot;
        for (int j = k; j < n; ++j)
            row_i[col_perm[j]] -= factor * row_k[col_perm[j]];
        b[i] -= factor * b[k];
    }

    pthread_barrier_wait(args->barrier);
    return nullptr;
}

// ===================================================
// Метод Гаусса с выбором главного по всей матрице
// ===================================================
bool gauss_solve(int n, double* A, double* b, double* x) {
    const int num_threads = std::min(4, n);
    pthread_t* threads = new pthread_t[num_threads];
    ThreadArgs* targs = new ThreadArgs[num_threads];
    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier, nullptr, num_threads);

    int* col_perm = new int[n];
    for (int j = 0; j < n; ++j) col_perm[j] = j;

    for (int k = 0; k < n; ++k) {
        // --- Поиск главного элемента ---
        double max_val = 0.0;
        int i_max = k, j_max = k;
        for (int i = k; i < n; ++i) {
            double* row_i = A + i * n;
            for (int j = k; j < n; ++j) {
                double val = std::abs(row_i[col_perm[j]]);
                if (val > max_val) {
                    max_val = val;
                    i_max = i;
                    j_max = j;
                }
            }
        }

        if (max_val == 0.0) {
            delete[] col_perm;
            delete[] threads;
            delete[] targs;
            pthread_barrier_destroy(&barrier);
            return false;
        }

        // --- Перестановка строк ---
        if (i_max != k) {
            double* row_k = A + k * n;
            double* row_im = A + i_max * n;
            for (int j = 0; j < n; ++j)
                std::swap(row_k[j], row_im[j]);
            std::swap(b[k], b[i_max]);
        }

        // --- Перестановка столбцов ---
        if (j_max != k)
            std::swap(col_perm[k], col_perm[j_max]);

        // --- Параллельное зануление ---
        double* row_k = A + k * n;
        double pivot = row_k[col_perm[k]];

        for (int t = 0; t < num_threads; ++t) {
            targs[t] = {t, n, k, A, b, col_perm, pivot, row_k, num_threads, &barrier};
            pthread_create(&threads[t], nullptr, elimination_worker, &targs[t]);
        }
        for (int t = 0; t < num_threads; ++t)
            pthread_join(threads[t], nullptr);
    }

    // Обратный ход (последовательный)
    for (int i = n - 1; i >= 0; --i) {
        double* row_i = A + i * n;
        double sum = b[i];
        for (int j = i + 1; j < n; ++j)
            sum -= row_i[col_perm[j]] * x[j];
        x[i] = sum / row_i[col_perm[i]];
    }

    pthread_barrier_destroy(&barrier);
    delete[] col_perm;
    delete[] threads;
    delete[] targs;
    return true;
}

// ===================================================
// compute_b — pthread-параллельная версия
// ===================================================
struct BArgs {
    int tid;
    int n;
    double* A;
    double* b;
    int num_threads;
};

void* compute_b_worker(void* arg) {
    BArgs* args = (BArgs*)arg;
    int n = args->n;
    int tid = args->tid;
    int num_threads = args->num_threads;
    double* A = args->A;
    double* b = args->b;

    for (int i = tid; i < n; i += num_threads) {
        double sum = 0.0;
        for (int k = 0; k <= (n - 1) / 2; ++k) {
            int j = 2 * k;
            if (j < n) sum += A[i * n + j];
        }
        b[i] = sum;
    }
    return nullptr;
}

void compute_b(int n, double* A, double* b) {
    const int num_threads = std::min(4, n);
    pthread_t threads[num_threads];
    BArgs args[num_threads];

    for (int t = 0; t < num_threads; ++t) {
        args[t] = {t, n, A, b, num_threads};
        pthread_create(&threads[t], nullptr, compute_b_worker, &args[t]);
    }
    for (int t = 0; t < num_threads; ++t)
        pthread_join(threads[t], nullptr);
}

// ===================================================
// compute_residual_norm — pthread-параллельная версия
// ===================================================
struct NormArgs {
    int tid;
    int n;
    const double* A;
    const double* x;
    const double* b;
    int num_threads;
    double local_num;
    double local_den;
};

void* residual_worker(void* arg) {
    NormArgs* args = (NormArgs*)arg;
    int n = args->n;
    int tid = args->tid;
    int num_threads = args->num_threads;
    const double* A = args->A;
    const double* x = args->x;
    const double* b = args->b;
    double local_num = 0.0;
    double local_den = 0.0;

    for (int i = tid; i < n; i += num_threads) {
        double sum = 0.0;
        for (int j = 0; j < n; ++j)
            sum += A[i * n + j] * x[j];
        local_num += (sum - b[i]) * (sum - b[i]);
        local_den += b[i] * b[i];
    }
    args->local_num = local_num;
    args->local_den = local_den;
    return nullptr;
}

double compute_residual_norm(int n, const double* A, const double* x, const double* b) {
    const int num_threads = std::min(4, n);
    pthread_t threads[num_threads];
    NormArgs args[num_threads];

    for (int t = 0; t < num_threads; ++t) {
        args[t] = {t, n, A, x, b, num_threads, 0.0, 0.0};
        pthread_create(&threads[t], nullptr, residual_worker, &args[t]);
    }
    double norm_num = 0.0, norm_den = 0.0;
    for (int t = 0; t < num_threads; ++t) {
        pthread_join(threads[t], nullptr);
        norm_num += args[t].local_num;
        norm_den += args[t].local_den;
    }
    return (norm_den != 0.0) ? sqrt(norm_num / norm_den) : 0.0;
}
