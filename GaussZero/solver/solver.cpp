#include <cmath>
#include <algorithm>
#include <iostream>
#include <limits>
#include "solver.h"

bool gauss_solve(int n, double* A, double* b, double* x) {
    int* col_perm = new int[n];
    for (int j = 0; j < n; ++j) col_perm[j] = j;

    // чтобы избежать деления на почти ноль.
    const double EPS = 1e-12;

    // Прямой ход
    for (int k = 0; k < n; ++k) {
        // Поиск главного элемента
        double max_val = 0.0;
        int i_max = k, j_max = k;
        for (int i = k; i < n; ++i) {
            double* row_i = A + i*n;
            for (int j = k; j < n; ++j) {
                double val = std::abs(row_i[col_perm[j]]);
                if (val > max_val) {
                    max_val = val;
                    i_max = i;
                    j_max = j;
                }
            }
        }

        // Проверка на "нулевую" подматрицу
        if (max_val < EPS) {
            delete[] col_perm;
            return false; 
        }

        // Перестановка строк
        if (i_max != k) {
            double* row_k = A + k*n;
            double* row_im = A + i_max*n;
            for (int j = 0; j < n; ++j) {
                std::swap(row_k[j], row_im[j]);
            }
            std::swap(b[k], b[i_max]);
        }

        // Перестановка столбцов
        if (j_max != k) {
            std::swap(col_perm[k], col_perm[j_max]);
        }

        // Зануление чисел под главной диагональю
        double* row_k = A + k*n;
        double pivot = row_k[col_perm[k]];
        for (int i = k+1; i < n; ++i) {
            double* row_i = A + i*n;
            double factor = row_i[col_perm[k]] / pivot;
            for (int j = k; j < n; ++j) {
                row_i[col_perm[j]] -= factor * row_k[col_perm[j]];
            }
            b[i] -= factor * b[k];
        }
    }

    // Обратный ход
    for (int i = n-1; i >= 0; --i) {
        double* row_i = A + i*n;
        double sum = b[i];
        for (int j = i+1; j < n; ++j) {
            sum -= row_i[col_perm[j]] * x[j];
        }
        x[i] = sum / row_i[col_perm[i]];
    }

    delete[] col_perm;
    return true;
}
