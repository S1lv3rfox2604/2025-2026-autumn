#include "solver.h"
#include <cmath>
#include <algorithm>
#include <cstring>

// Вспомогательная функция для вычисления нормы вектора
static double vector_norm(int n, double* v) {
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        sum += v[i] * v[i];
    }
    return std::sqrt(sum);
}

// Вспомогательная функция для вычисления скалярного произведения
static double dot_product(int n, double* a, double* b) {
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        sum += a[i] * b[i];
    }
    return sum;
}

// Приведение симметричной матрицы к трёхдиагональному виду методом отражений
static void tridiagonal(int n, double* A, double* d, double* e) {
    if (n <= 1) {
        if (n == 1) d[0] = A[0];
        return;
    }
    
    // Рабочая копия матрицы
    double* Q = new double[n*n];
    std::memcpy(Q, A, n*n*sizeof(double));
    
    for (int k = 0; k < n-2; ++k) {
        // Вычисляем вектор x (столбец под диагональю)
        double* x = new double[n-k-1];//тут с паматью беда
        for (int i = 0; i < n-k-1; ++i) {
            x[i] = Q[(k+i+1)*n + k];
        }
        
        double norm_x = vector_norm(n-k-1, x);
        if (norm_x < 1e-15) {
            e[k] = 0.0;
            delete[] x;
            continue;
        }
        
        // Вектор отражения v
        double* v = new double[n-k-1];
        std::memcpy(v, x, (n-k-1)*sizeof(double));
        v[0] += (x[0] >= 0 ? norm_x : -norm_x);
        
        double norm_v = vector_norm(n-k-1, v);
        if (norm_v < 1e-15) {
            delete[] x;
            delete[] v;
            continue;
        }
        
        // Нормализуем v
        for (int i = 0; i < n-k-1; ++i) {
            v[i] /= norm_v;
        }
        
        // Применяем преобразование к подматрице
        for (int i = k; i < n; ++i) {
            double* col = new double[n-k-1];
            for (int j = 0; j < n-k-1; ++j) {
                col[j] = Q[(k+j+1)*n + i];
            }
            
            double dot = dot_product(n-k-1, v, col);
            
            for (int j = 0; j < n-k-1; ++j) {
                Q[(k+j+1)*n + i] -= 2.0 * v[j] * dot;
            }
            
            delete[] col;
        }
        
        // Аналогично для строк
        for (int i = k; i < n; ++i) {
            double* row = new double[n-k-1];
            for (int j = 0; j < n-k-1; ++j) {
                row[j] = Q[i*n + (k+j+1)];
            }
            
            double dot = dot_product(n-k-1, v, row);
            
            for (int j = 0; j < n-k-1; ++j) {
                Q[i*n + (k+j+1)] -= 2.0 * v[j] * dot;
            }
            
            delete[] row;
        }
        
        // Сохраняем поддиагональный элемент
        e[k] = (x[0] >= 0 ? -norm_x : norm_x);
        
        delete[] x;
        delete[] v;
    }
    
    // Извлекаем диагональ и поддиагональ
    for (int i = 0; i < n; ++i) {
        d[i] = Q[i*n + i];
    }
    
    // Для последнего поддиагонального элемента
    if (n > 1) {
        e[n-2] = Q[(n-1)*n + (n-2)];
    }
    
    delete[] Q;
}

// LU-разложение трёхдиагональной матрицы T - lambda*I
// Возвращает количество отрицательных собственных значений < lambda
static int count_negative_pivots(int n, double* d, double* e, double lambda) {
    int count = 0;
    double u = d[0] - lambda;
    
    // Первый пивот
    if (u < 0.0) count++;
    if (std::abs(u) < 1e-15) u = 1e-15; // для устойчивости
    
    for (int i = 1; i < n; ++i) {
        // L[i] = e[i-1] / u_prev
        // U[i] = d[i] - lambda - L[i] * e[i-1]
        double l = e[i-1] / u;
        u = d[i] - lambda - l * e[i-1];
        
        if (u < 0.0) count++;
        if (std::abs(u) < 1e-15) u = 1e-15; // для устойчивости
    }
    
    return count;
}

// Метод бисекции для нахождения k-го собственного значения
double find_kth_eigenvalue(int n, double* A, int k, double eps) {
    if (n <= 0 || k <= 0 || k > n) {
        return 0.0;
    }
    
    // Для матрицы 1x1 сразу возвращаем элемент
    if (n == 1) {
        return A[0];
    }
    
    double* d = new double[n];
    double* e = new double[n-1];
    
    // Приводим матрицу к трёхдиагональному виду
    tridiagonal(n, A, d, e);
    
    // Определяем границы интервала для бисекции
    double left = d[0];
    double right = d[0];
    
    for (int i = 0; i < n; ++i) {
        // Оценка по теореме Гершгорина
        double radius = 0.0;
        if (i > 0) radius += std::abs(e[i-1]);
        if (i < n-1) radius += std::abs(e[i]);
        
        left = std::min(left, d[i] - radius);
        right = std::max(right, d[i] + radius);
    }
    
    // Добавляем запас для надёжности
    left -= 1.0;
    right += 1.0;
    
    // Метод бисекции
    int iterations = 0;
    const int max_iterations = 100;
    
    while (right - left > eps && iterations < max_iterations) {
        double mid = 0.5 * (left + right);
        int count_less = count_negative_pivots(n, d, e, mid);
        
        if (count_less < k) {
            left = mid;
        } else {
            right = mid;
        }
        iterations++;
    }
    
    double result = 0.5 * (left + right);
    
    delete[] d;
    delete[] e;
    
    return result;
}