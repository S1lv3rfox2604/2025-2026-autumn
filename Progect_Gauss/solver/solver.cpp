#include <cmath>
#include <algorithm>
#include <iostream>
#include "solver.h"

// Метод Гаусса с выбором главного элемента (работаем с A и b отдельно)
bool gauss_solve(int n, double* A, double* b, double* x) {
    // Прямой ход
    for (int k = 0; k < n; ++k) {
        // Поиск главного элемента
        double max_val = 0.0;
        int i_max = k, j_max = k;
        for (int i = k; i < n; ++i) {
            for (int j = k; j < n; ++j) {
                if (std::abs(A[i*n + j]) > max_val) {
                    max_val = std::abs(A[i*n + j]);
                    i_max = i;
                    j_max = j;
                }
            }
        }

        if (max_val == 0.0) {
            return false; // Полностью нулевая подматрица
        }

        // Перестановка строк
        if (i_max != k) {
            for (int j = 0; j < n; ++j) {
                std::swap(A[k*n + j], A[i_max*n + j]);
            }
            std::swap(b[k], b[i_max]);
        }

        // Перестановка столбцов
        if (j_max != k) {
            for (int i = 0; i < n; ++i) {
                std::swap(A[i*n + k], A[i*n + j_max]);
            }
        }

        // Зануление чисел под главной диагональю
        double pivot = A[k*n + k];
        for (int i = k+1; i < n; ++i) {
            double factor = A[i*n + k] / pivot;
            for (int j = k; j < n; ++j) {
                A[i*n + j] -= factor * A[k*n + j];
            }
            b[i] -= factor * b[k];
        }
    }

    // Обратный ход
    for (int i = n-1; i >= 0; --i) {
        double sum = b[i];
        for (int j = i+1; j < n; ++j)
            sum -= A[i*n + j] * x[j];
        x[i] = sum / A[i*n + i];
    }

    return true;
}

