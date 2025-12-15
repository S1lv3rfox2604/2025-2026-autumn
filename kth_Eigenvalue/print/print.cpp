#include <iostream>
#include <iomanip>
#include <algorithm>
#include "print.h"

void print_matrix(int n, int ncols, const double* A, int m)
{
    int mm = std::min(std::min(n, ncols), m);

    std::cout << std::fixed << std::setprecision(6);

    for (int i = 0; i < mm; ++i) {
        for (int j = 0; j < mm; ++j) {
            std::cout << std::setw(12) << A[i * ncols + j] << " ";
        }
        std::cout << "\n";
    }
}
