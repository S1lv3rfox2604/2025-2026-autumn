#include <iostream>
#include <fstream>
#include <cmath>
#include "init.h"

// ------------------------------------------------------------
// ЧТЕНИЕ МАТРИЦЫ ИЗ ФАЙЛА
// ------------------------------------------------------------
bool read_matrix_from_file(int n, double* A, const char* filename) {
    std::ifstream in(filename);
    if (!in) return false;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (!(in >> A[i*n + j])) return false;
        }
    }
    return true;
}


// ------------------------------------------------------------
// ИНИЦИАЛИЗАЦИЯ МАТРИЦЫ ПО ФОРМУЛАМ 1–4
// ------------------------------------------------------------
bool init_matrix_formula(int k, int n, double* A) {

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {

            switch(k) {

                case 1:
                    A[i*n + j] = n - std::max(i+1, j+1) + 1;
                    break;

                case 2:
                    if (i == j)
                        A[i*n + j] = 2.0;
                    else if (std::abs(i - j) == 1)
                        A[i*n + j] = -1.0;
                    else
                        A[i*n + j] = 0.0;
                    break;

                case 3:
                    if (i == j && i != n-1)
                        A[i*n + j] = 1.0;

                    else if (j == n-1)    // последний столбец
                        A[i*n + j] = (double)(i + 1);

                    else if (i == n-1)   // последняя строка
                        A[i*n + j] = (double)(j + 1);

                    else
                        A[i*n + j] = 0.0;

                    break;

                case 4:
                    A[i*n + j] = 1.0 / (double)(i + j + 1);
                    break;

                default:
                    return false;
            }
        }
    }
    return true;
}
