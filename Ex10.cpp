#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>

constexpr double EPSILON = 1.0e-10;

enum Method {
    LIEBMAN,
    SOR
};

void solve_poisson(int n, double h, Method methodType, double omega) {
    // Grid size is (n+1) x (n+1) to include boundaries from 0 to 1
    int size = n + 1;
    std::vector<std::vector<double>> u(size, std::vector<double>(size, 0.0));

    int iter = 0;
    double max_error;
    
    // Constant term from the equation: 4*u_ij = u_neighbors + 2*h^2
    // So u_ij_GS = 0.25 * (u_neighbors + 2*h^2)
    double constant_term = 2.0 * h * h;
    do {
        max_error = 0.0;
        iter++;

        for (int i = 1; i < n; ++i) {     // x index
            for (int j = 1; j < n; ++j) { // y index
                
                // u_new = (u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] + 2h^2) / 4
                double u_gs = (u[i + 1][j] + u[i - 1][j] + u[i][j + 1] + u[i][j - 1] + constant_term) / 4.0;

                double u_new;
                if (methodType == LIEBMAN) {
                    u_new = u_gs;
                } else {
                    //u_new = (1-w)*u_old + w*u_gs
                    u_new = (1.0 - omega) * u[i][j] + omega * u_gs;
                }

                double current_error = std::abs(u_new - u[i][j]);
                if (current_error > max_error) {
                    max_error = current_error;
                }

                u[i][j] = u_new;
            }
        }
        if (iter > 100000) {
            std::cout << "Max iterations reached!" << std::endl;
            break;
        }

    } while (max_error > EPSILON);

    std::string methodName = (methodType == LIEBMAN) ? "Liebman method" : "SOR method, overrelaxation factor = 1.6";
    
    std::cout << methodName << std::endl;
    std::cout << "n = " << n << "   h = " << std::fixed << std::setprecision(3) << h << std::endl;
    std::cout << "iter = " << iter << std::endl;

    // Print values along the line y = 0.5
    // The index for y = 0.5 is n / 2
    int y_index = n / 2;
    
    std::cout << std::fixed << std::setprecision(2);
    
    for (int i = 0; i <= n; ++i) {
        double x = i * h;
        double val = u[i][y_index];
        
        std::cout << "x = " << x << "   u = " 
                  << std::fixed << std::setprecision(5) << val << std::endl;
        
        // Reset precision for x in next loop iteration
        std::cout << std::fixed << std::setprecision(2);
    }
    std::cout << std::endl; // Spacing between methods
}

int main() {
    std::cout << "d^2u/dx^2 + d^2u/dy^2 + 2 = 0, edge length=1" << std::endl << std::endl;

    // --- Case 1: h = 0.1 (n = 10) ---
    int n1 = 10;
    double h1 = 0.1;
    
    // a) Liebman
    solve_poisson(n1, h1, LIEBMAN, 1.0);
    
    // b) SOR (omega = 1.6)
    solve_poisson(n1, h1, SOR, 1.6);

    // --- Case 2: h = 0.05 (n = 20) ---
    int n2 = 20;
    double h2 = 0.05;
    
    // a) Liebman
    solve_poisson(n2, h2, LIEBMAN, 1.0);
    
    // b) SOR (omega = 1.6)
    solve_poisson(n2, h2, SOR, 1.6);

    return 0;
}