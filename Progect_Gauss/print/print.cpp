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
