#include <iostream>
#include <cstdlib>
#include <ctime>
#include "init.h"
#include "solver.h"
#include "print.h"

{
    if (argc < 6) {
        std::cerr << "Использование:\n"
                  << "  Инициализация по формуле: ./main n m eps k eig_index   (k = 1..4)\n"
                  << "  Чтение из файла:         ./main n m eps 0 filename eig_index\n";
        return 1;
    }

    int n = std::atoi(argv[1]);
    int m = std::atoi(argv[2]);
    double eps = std::atof(argv[3]);
    int k = std::atoi(argv[4]);

    if (n <= 0 || m <= 0 || eps <= 0.0) {
        std::cerr << "Ошибка: неверные параметры n, m или eps\n";
        return 1;
    }

    const char* filename = nullptr;
    int eig_index = 1;

    if (k == 0) {
        // ожидаем argv[5] = filename, argv[6] = eig_index
        if (argc < 7) {
            std::cerr << "Ошибка: при k=0 нужно указать filename и eig_index\n";
            return 1;
        }
        filename = argv[5];
        eig_index = std::atoi(argv[6]);
    } else {
        // k != 0: argv[5] = eig_index
        if (argc < 6) {
            std::cerr << "Ошибка: недостаточно аргументов для инициализации по формуле\n";
            return 1;
        }
        eig_index = std::atoi(argv[5]);
    }

    if (k < 0 || k > 4) {
        std::cerr << "Ошибка: параметр k должен быть 0..4 (0 — чтение из файла, 1..4 — формулы)\n";
        return 1;
    }

    if (eig_index < 1 || eig_index > n) {
        std::cerr << "Ошибка: eig_index должен быть в диапазоне 1.." << n << "\n";
        return 1;
    }

    // выделяем матрицу
    double* A = new double[n * n];
    bool ok = false;

    if (k == 0) {
        ok = read_matrix_from_file(n, A, filename);
        if (!ok) {
            std::cerr << "Ошибка чтения файла: " << filename << "\n";
            delete[] A;
            return 1;
        }
    } else {
        ok = init_matrix_formula(k, n, A);
        if (!ok) {
            std::cerr << "Ошибка инициализации матрицы формулой k=" << k << "\n";
            delete[] A;
            return 1;
        }
    }

    // вывод фрагмента матрицы
    std::cout << "Фрагмент матрицы A (" << std::min(n,m) << "x" << std::min(n,m) << "):\n";
    print_matrix(n, n, A, m);

    // поиск k-го собственного (eig_index) — 1-based, минимальное = 1
    std::clock_t t0 = std::clock();
    double eig = find_kth_eigenvalue(n, A, eig_index, eps);
    std::clock_t t1 = std::clock();

    std::cout << "\nСобственное значение (индекс " << eig_index << "): " << eig << "\n";
    std::cout << "Время поиска: " << double(t1 - t0) / CLOCKS_PER_SEC << " сек\n";

    delete[] A;
    return 0;
}