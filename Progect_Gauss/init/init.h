#ifndef INIT_H
#define INIT_H

bool read_matrix_from_file(int n, double* A, const char* filename);
void compute_b(int n, double* A, double* b);
bool init_matrix_formula(int k, int n, double* A);

#endif
