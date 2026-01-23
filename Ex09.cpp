#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <string>

double dydx(double x, double y) {
    return x * y;
}

double exact_solution(double x) {
    return std::exp(0.5 * x * x);
}

// y_next = y_current + h * dydx(x_current, y_current)
double solve_euler(double h, double x0, double y0, double x_target) {
    double x = x0;
    double y = y0;
    int steps = std::round((x_target - x0) / h);
    for (int i = 0; i < steps; ++i) {
        y = y + h * dydx(x, y);
        x += h;
    }
    return y;
}

// weighted average of 4 slopes
double solve_rk4(double h, double x0, double y0, double x_target) {
    double x = x0;
    double y = y0;
    int steps = std::round((x_target - x0) / h);

    for (int i = 0; i < steps; ++i) {
        double k1 = h * dydx(x, y);
        double k2 = h * dydx(x + h / 2.0, y + k1 / 2.0);
        double k3 = h * dydx(x + h / 2.0, y + k2 / 2.0);
        double k4 = h * dydx(x + h, y + k3);

        y = y + (1.0 / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        x += h;
    }
    return y;
}

int main() {
    const double x0 = 0.0;
    const double y0 = 1.0;
    const double x_target = 1.0;
    
    double y_theory = exact_solution(x_target);
    
    std::cout << "dy/dx = xy, a = 0, b = 1, y(a) = 1" << std::endl;
    std::vector<double> h_values = {0.1, 0.01};

    for (double h : h_values) {
        // Euler
        double y_euler = solve_euler(h, x0, y0, x_target);
        double err_euler = y_euler - y_theory;
        std::cout << std::left << std::setw(13) << "Euler" << ": "
              << std::fixed << std::setprecision(3) 
              << "h = " << h << "  "
              << "x = " << x_target << "  "
              << std::setprecision(8) 
              << "y = " << y_euler << "  "
              << std::scientific << std::setprecision(4)
              << "Error = " << err_euler << std::endl;

        // Runge-Kutta
        double y_rk4 = solve_rk4(h, x0, y0, x_target);
        double err_rk4 = y_rk4 - y_theory;
        std::cout << std::left << std::setw(13) << "Runge-Kutta" << ": "
              << std::fixed << std::setprecision(3) 
              << "h = " << h << "  "
              << "x = " << x_target << "  "
              << std::setprecision(8) 
              << "y = " << y_rk4 << "  "
              << std::scientific << std::setprecision(4)
              << "Error = " << err_rk4 << std::endl;
    }

    return 0;
}