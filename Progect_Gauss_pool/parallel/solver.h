#ifndef PARALLEL_SOLVER_H
#define PARALLEL_SOLVER_H

#include <pthread.h>

bool gauss_solve_parallel(int n, double* A, double* b, double* x, int num_threads);

// Структура для статических данных потоков
typedef struct {
    int n;
    double* A;
    double* b;
    int* col_perm;
    int num_threads;
    int thread_id;
} static_thread_data_t;

#endif