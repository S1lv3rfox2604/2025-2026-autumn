#include <iostream>
#include <cstdlib>
#include <ctime>
#include <pthread.h>
#include <sys/sysinfo.h>
#include "init.h"
#include "print.h"

// Подключаем обе версии solver
#include "sequential/solver.h"
#include "parallel/solver.h"

int main(int argc, char* argv[]) {
    if (argc != 4 && argc != 5 && argc != 6) {
        std::cerr << "Usage: " << argv[0] << " n m k [num_threads] [filename]\n";
        std::cerr << "Examples:\n";
        std::cerr << "  " << argv[0] << " 2000 6 1           # sequential\n";
        std::cerr << "  " << argv[0] << " 2000 6 1 4         # parallel, 4 threads\n";
        std::cerr << "  " << argv[0] << " 2000 6 0 matrix.txt # from file\n";
        return 1;
    }

    int n = std::atoi(argv[1]);
    int m = std::atoi(argv[2]);
    int k = std::atoi(argv[3]);
    int num_threads = 1;
    const char* filename = nullptr;

    // Разбор аргументов
    if (argc == 5) {
        if (k == 0) {
            filename = argv[4];
        } else {
            num_threads = std::atoi(argv[4]);
        }
    } else if (argc == 6) {
        num_threads = std::atoi(argv[4]);
        filename = argv[5];
    }

    // Корректировка количества потоков
    if (num_threads <= 0) {
        num_threads = 1;
    }
    if (num_threads > n) {
        num_threads = n;
    }

    // Выделение памяти
    double* A = nullptr;
    double* b = nullptr;
    double* x = nullptr;
    
    try {
        A = new double[n * n];
        b = new double[n];
        x = new double[n];
    } catch (const std::bad_alloc& e) {
        std::cerr << "Memory allocation failed: " << e.what() << "\n";
        return 1;
    }

    bool ok = false;

    // Инициализация матрицы
    if (k == 0) {
        if (filename == nullptr) {
            std::cerr << "Filename required when k = 0\n";
            delete[] A; delete[] b; delete[] x;
            return 1;
        }
        ok = read_matrix_from_file(n, A, filename);
        if (!ok) {
            std::cerr << "Cannot open or read file: " << filename << "\n";
            delete[] A; delete[] b; delete[] x;
            return 1;
        }
    } else if (k >= 1 && k <= 4) {
        ok = init_matrix_formula(k, n, A);
        if (!ok) {
            std::cerr << "Error initializing matrix with formula " << k << "\n";
            delete[] A; delete[] b; delete[] x;
            return 1;
        }
    } else {
        std::cerr << "Invalid formula k: " << k << "\n";
        delete[] A; delete[] b; delete[] x;
        return 1;
    }

    std::cout << "=== Gauss Solver ===\n";
    std::cout << "Matrix size: " << n << "x" << n << "\n";
    std::cout << "Threads: " << num_threads << "\n";
    std::cout << "Matrix formula: " << k << "\n";
    if (filename) std::cout << "Input file: " << filename << "\n";
    
    // Выбор версии решателя
    std::string solver_version;
    if (num_threads == 1) {
        solver_version = "SEQUENTIAL";
        std::cout << "Solver: Sequential version\n";
    } else {
        solver_version = "PARALLEL";
        std::cout << "Solver: Parallel version (" << num_threads << " threads)\n";
    }
    std::cout << "=============================\n";

    std::cout << "Matrix A:\n";
    print_matrix(n, n, A, m);

    // Начало измерения общего времени
    clock_t total_start = clock();

    //Вычисление вектора b
    clock_t b_start = clock();
    compute_b_parallel(n, A, b, num_threads);
    clock_t b_end = clock();
    double b_time = double(b_end - b_start) / CLOCKS_PER_SEC;

    std::cout << "Vector b:\n";
    print_vector(b, n, m);

    //Решение СЛАУ
    clock_t solve_start = clock();
    
    if (num_threads == 1) {
        ok = gauss_solve_sequential(n, A, b, x);
    } else {
        ok = gauss_solve_parallel(n, A, b, x, num_threads);
    }
    
    clock_t solve_end = clock();
    double solve_time = double(solve_end - solve_start) / CLOCKS_PER_SEC;

    double residual_norm = 0.0;
    double residual_time = 0.0;

    if (!ok) {
        std::cerr << "Solver failed - matrix is singular\n";
    } else {
        std::cout << "Solution x:\n";
        print_vector(x, n, m);

        //Вычисление невязки
        clock_t residual_start = clock();
        residual_norm = compute_residual_norm(n, A, x, b, num_threads);
        clock_t residual_end = clock();
        residual_time = double(residual_end - residual_start) / CLOCKS_PER_SEC;

        std::cout << "Residual norm ||Ax-b||/||b|| = " << residual_norm << "\n";
    }

    // Конец измерения общего времени
    clock_t total_end = clock();
    double total_time = double(total_end - total_start) / CLOCKS_PER_SEC;

    // Создавалось для того чтобы понять где потеря по времени
    std::cout << "\n=== PERFORMANCE RESULTS ===\n";
    std::cout << "Solver version:       " << solver_version << "\n";
    std::cout << "Compute b time:        " << b_time << " sec\n";
    std::cout << "Gauss solve time:      " << solve_time << " sec\n";
    if (ok) {
        std::cout << "Residual compute time: " << residual_time << " sec\n";
    }
    std::cout << "Total time:            " << total_time << " sec\n";
    std::cout << "Threads used:          " << num_threads << "\n";

    delete[] A;
    delete[] b;
    delete[] x;
    return 0;
}