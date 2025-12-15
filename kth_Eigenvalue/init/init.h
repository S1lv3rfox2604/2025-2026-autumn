#ifndef INIT_H
#define INIT_H

#include <cstddef>

// read matrix from file: A is stored row-major length n*n
bool read_matrix_from_file(int n, double* A, const char* filename);

// initialize by formula number k (1..4)
bool init_matrix_formula(int k, int n, double* A);

#endif // INIT_H

