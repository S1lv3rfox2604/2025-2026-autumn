#include <iostream>
#include <fstream>
#include <cmath>
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

// Вычисление вектора b
/*void compute_b(int n, double* A, double* b) {
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int k = 0; k <= (n-1)/2; ++k) {
            int j = 2*k;
            if (j < n) sum += A[i * n + j];
        }
        b[i] = sum;
    }
}*/

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
