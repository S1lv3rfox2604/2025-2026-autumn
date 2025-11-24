#ifndef UTIL_H
#define UTIL_H

#include <pthread.h>

void print_matrix(int rows, int cols, const double* A, int m);
void print_vector(const double* b, int n, int m);
double compute_residual_norm(int n, const double* A, const double* x, const double* b, int num_threads);

// Структура для параллельного вычисления невязки
typedef struct {
    int n;
    const double* A;
    const double* x;
    const double* b;
    int start_row;
    int end_row;
    double local_norm_num;
    double local_norm_den;
} residual_data_t;

#endif