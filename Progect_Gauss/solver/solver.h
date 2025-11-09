#ifndef SOLVER_H
#define SOLVER_H

bool gauss_solve(int n, double* A, double* b, double* x);
void compute_b(int n, double* A, double* b);
double compute_residual_norm(int n, const double* A, const double* x, const double* b);

#endif
