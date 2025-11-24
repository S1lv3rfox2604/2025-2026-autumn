#include <iostream>
#include <fstream>
#include <cmath>
#include <pthread.h>
#include "init.h"

// Чтение матрицы из файла
bool read_matrix_from_file(int n, double* A, const char* filename) {
    std::ifstream in(filename);
    if (!in) return false;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (!(in >> A[i * n + j])) return false;
        }
    }
    return true;
}

// Заполнение матрицы по формуле
bool init_matrix_formula(int k, int n, double* A) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            switch(k) {
                case 1:
                    A[i*n + j] = n - std::max(i+1,j+1) + 1;
                    break;
                case 2:
                    A[i*n + j] = std::max(i+1,j+1);
                    break;
                case 3:
                    A[i*n + j] = std::abs(i-j);
                    break;
                case 4:
                    A[i*n + j] = 1.0 / (i + j + 1);
                    break;
                default:
                    return false;
            }
        }
    }
    return true;
}

// Функция для потока вычисления вектора b
void* compute_b_thread(void* arg) {
    thread_data_t* data = (thread_data_t*)arg;
    
    for (int i = data->start_row; i < data->end_row; ++i) {
        double sum = 0.0;
        for (int k = 0; k <= (data->n-1)/2; ++k) {
            int j = 2*k;
            if (j < data->n) sum += data->A[i * data->n + j];
        }
        data->b[i] = sum;
    }
    
    pthread_exit(NULL);
}

// Параллельная версия вычисления вектора b
void compute_b_parallel(int n, double* A, double* b, int num_threads) {
   
    pthread_t threads[num_threads];
    thread_data_t thread_data[num_threads];
    
    int rows_per_thread = n / num_threads;
    int extra_rows = n % num_threads;
    
    int current_row = 0;
    for (int i = 0; i < num_threads; ++i) {
        thread_data[i].n = n;
        thread_data[i].A = A;
        thread_data[i].b = b;
        thread_data[i].start_row = current_row;
        
        int rows_for_this_thread = rows_per_thread + (i < extra_rows ? 1 : 0);
        thread_data[i].end_row = current_row + rows_for_this_thread;
        
        current_row = thread_data[i].end_row;
        
        pthread_create(&threads[i], NULL, compute_b_thread, &thread_data[i]);
    }
    
    // Ожидание завершения всех потоков
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(threads[i], NULL);
    }
}