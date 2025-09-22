#include <iostream>
#include <cstdlib>
#include <ctime>
#include "init.h"
#include "solver.h"
#include "print.h"

int main(int argc, char* argv[]) {
    if (argc < 4 || argc > 5) {
        std::cerr << "Usage: " << argv[0] << " n m k [filename]\n";
        return 1;
    }

    int n = std::atoi(argv[1]);
    int m = std::atoi(argv[2]);
    int k = std::atoi(argv[3]);

    double* A = new double[n * n];
    double* b = new double[n];
    double* x = new double[n];

    bool ok = false;

    if (k == 0) {
        if (argc != 5) {
            std::cerr << "Filename required when k = 0\n";
            delete[] A; delete[] b; delete[] x;
            return 1;
        }
        const char* filename = argv[4];
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

    compute_b(n, A, b);

    std::cout << "Matrix A:\n";
    print_matrix(n, n, A, m);
    std::cout << "Vector b:\n";
    print_vector(b, n, m);

    clock_t start = clock();
    ok = gauss_solve(n, A, b, x);
    clock_t end = clock();

    if (!ok) {
        std::cerr << "Solver failed\n";
    } else {
        std::cout << "Solution x:\n";
        print_vector(x, n, m);

        double res_norm = compute_residual_norm(n, A, x, b);
        std::cout << "Residual norm ||Ax-b||/||b|| = " << res_norm << "\n";
    }

    std::cout << "Time elapsed: " << double(end - start) / CLOCKS_PER_SEC << " sec\n";

    delete[] A;
    delete[] b;
    delete[] x;
    return 0;
}
