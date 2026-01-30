#include "LUBand.h"

int main(){
    int n = 6, h = 3, h1 = h - 1;
    
    int maxIdx = (n - 1) + h1 * n;
    std::vector<double> a(maxIdx + 1, 0.0);

    // set U (upper triangular matrix)
    auto set_matrix = [&](int row, int col, double val){
        if (row > col) return;
        int idx = row + h1 * (col + 1);
        a[idx] = val; 
    };

    // manually set the matrix
    set_matrix(0, 0, 3.0);  set_matrix(0, 1, -1.0); set_matrix(1, 1, 2.0);
    set_matrix(0, 2, -1.0); set_matrix(1, 2, 0.0);  set_matrix(2, 2, 3.0);
    set_matrix(1, 3, -2.0); set_matrix(2, 3, -1.0); set_matrix(3, 3, 2.0);
    set_matrix(2, 4, 0.0);  set_matrix(3, 4, 0.0);  set_matrix(4, 4, 1.0);
    set_matrix(3, 5, 0.0);  set_matrix(4, 5, -2.0); set_matrix(5, 5, 1.0);

    std::vector<double> b = {-2.0, -5.0, 4.0, 1.0, -7.0, -4.0}; // RHS 

    std::cout << "h = " << h << ", n = " << n << std::endl;

    //solve
    LUband solver;
    solver.decompose(a, n, h);
    solver.solve(a, n, h, b);

    std::cout << "Solution: " << std::endl;
    std::cout << std::fixed << std::setprecision(4);
    for(int i = 0; i < n; ++i){
        std::cout << "x" << i + 1 << " = " << b[i] << std::endl;
    }

    return 0;
}