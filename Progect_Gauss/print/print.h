#ifndef UTIL_H
#define UTIL_H

void print_matrix(int rows, int cols, const double* A, int m);
void print_vector(const double* b, int n, int m);
double compute_residual_norm(int n, const double* A, const double* x, const double* b);

#endif
