#include <cstdio>
#include <cmath>
#include <algorithm>
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

double compute_residual_norm(int n, const double* A, const double* x, const double* b) {
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
