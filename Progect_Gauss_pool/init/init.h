#ifndef INIT_H
#define INIT_H

#include <pthread.h>

bool read_matrix_from_file(int n, double* A, const char* filename);
bool init_matrix_formula(int k, int n, double* A);

// Структура для передачи данных в потоки
typedef struct {
    int n;
    double* A;
    double* b;
    int start_row;
    int end_row;
} thread_data_t;

void compute_b_parallel(int n, double* A, double* b, int num_threads);

#endif