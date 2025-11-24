#include <cstdio>
#include <cmath>
#include <algorithm>
#include <pthread.h>
#include "print.h"

void print_matrix(int rows, int cols, const double* A, int m) {
    int r = std::min(rows, m);
    int c = std::min(cols, m);
    for (int i = 0; i < r; ++i) {
        for (int j = 0; j < c; ++j) {
            printf("%10.3e ", A[i * cols + j]);
        }
        printf("\n");
    }
}

void print_vector(const double* b, int n, int m) {
    int r = std::min(n, m);
    for (int i = 0; i < r; ++i) {
        printf("%10.3e ", b[i]);
    }
    printf("\n");
}

// Функция потока для вычисления невязки
void* residual_thread(void* arg) {
    residual_data_t* data = (residual_data_t*)arg;
    
    data->local_norm_num = 0.0;
    data->local_norm_den = 0.0;
    
    for (int i = data->start_row; i < data->end_row; ++i) {
        double sum = 0.0;
        for (int j = 0; j < data->n; ++j) {
            // Используем исходную матрицу A и решение x
            sum += data->A[i * data->n + j] * data->x[j];
        }
        // Используем исходный вектор b
        data->local_norm_num += (sum - data->b[i]) * (sum - data->b[i]);
        data->local_norm_den += data->b[i] * data->b[i];
    }
    
    pthread_exit(NULL);
}

// Параллельная версия вычисления невязки
double compute_residual_norm(int n, const double* A, const double* x, const double* b, int num_threads) {
    if (num_threads == 1) {
        // Для 1 потока - последовательная версия в основном потоке
        double norm_num = 0.0;
        double norm_den = 0.0;
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < n; ++j) {
                sum += A[i * n + j] * x[j];
            }
            norm_num += (sum - b[i]) * (sum - b[i]);
            norm_den += b[i] * b[i];
        }
        return (norm_den != 0.0) ? sqrt(norm_num / norm_den) : 0.0;
    }
    
    pthread_t threads[num_threads];
    residual_data_t thread_data[num_threads];
    
    int rows_per_thread = n / num_threads;
    int extra_rows = n % num_threads;
    
    int current_row = 0;
    for (int i = 0; i < num_threads; ++i) {
        thread_data[i].n = n;
        thread_data[i].A = A;
        thread_data[i].x = x;
        thread_data[i].b = b;
        thread_data[i].start_row = current_row;
        
        int rows_for_this_thread = rows_per_thread + (i < extra_rows ? 1 : 0);
        thread_data[i].end_row = current_row + rows_for_this_thread;
        
        current_row = thread_data[i].end_row;
        
        pthread_create(&threads[i], NULL, residual_thread, &thread_data[i]);
    }
    
    double total_norm_num = 0.0;
    double total_norm_den = 0.0;
    
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(threads[i], NULL);
        total_norm_num += thread_data[i].local_norm_num;
        total_norm_den += thread_data[i].local_norm_den;
    }
    
    return (total_norm_den != 0.0) ? sqrt(total_norm_num / total_norm_den) : 0.0;
}