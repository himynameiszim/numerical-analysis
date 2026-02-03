#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

class Ex12 {
private:
    long long r;
    long long a, b;
    long long m;

public:
    Ex12() : r(0), a(0), b(0), m(0) {}

    void InitGenerator(long long seed, long long c1, long long c2, long long mod) {
        r = seed;
        a = c1;
        b = c2;
        m = mod;
    }

    double Rand() {
        r = (a * r + b) % m;
        return (double)r / (double)m;
    }

    static double f(double x) { 
        return 1.0 / (1.0 + x * x); 
    }
    
    static double exact() { 
        return 0.25 * M_PI; 
    }

    // Function for 2D Ellipse
    static double f2(double x, double y) { 
        return 13.0 * x * x + 34.0 * x * y + 25.0 * y * y - 1.0; 
    }
    
    static double exact2() { 
        return M_PI / 6.0; 
    }

    // Function for 3D Sphere
    static double f3(double x, double y, double z) {
        return (x - 1.0) * (x - 1.0) + (y - 0.5) * (y - 0.5) + z * z - 1.0;
    }
    
    static double exact3() { 
        return 4.0 * M_PI / 3.0; 
    }

    // Type 1: Mean Value Method
    void MonteCarlo1D(int n) {
        double sum = 0.0;

        for (int i = 0; i < n; i++) {
            double x = Rand(); // x is in [0, 1]
            sum += f(x);
        }

        // Integral = (b-a) * average. Here (b-a) = (1-0) = 1.
        double d = sum / n; 

        std::cout << "Num: " << n 
                  << " Result = " << d 
                  << " Error = " << (d - exact()) << std::endl;
    }

    // Type 2: Hit-or-Miss Method (2D)
    void MonteCarlo2D(int n) {
        int hits = 0;

        for (int i = 0; i < n; i++) {
            // Generate x in [-1, 1]
            double x = -1.0 + 2.0 * Rand();
            // Generate y in [-1, 1]
            double y = -1.0 + 2.0 * Rand();

            if (f2(x, y) < 0) {
                hits++;
            }
        }

        // Bounding Box Area = 2 * 2 = 4.0
        double d = 4.0 * (double)hits / (double)n;

        std::cout << "Num: " << n 
                  << " Result = " << d 
                  << " Error = " << (d - exact2()) << std::endl;
    }

    // Type 2: Hit-or-Miss Method (3D)
    void MonteCarlo3D(int n) {
        int hits = 0;

        for (int i = 0; i < n; i++) {
            // Generate x in [0, 2]
            double x = 0.0 + 2.0 * Rand();
            // Generate y in [-0.5, 1.5]
            double y = -0.5 + 2.0 * Rand();
            // Generate z in [-1, 1]
            double z = -1.0 + 2.0 * Rand();

            if (f3(x, y, z) < 0) {
                hits++;
            }
        }

        // Bounding Box Volume = 2 * 2 * 2 = 8.0
        double d = 8.0 * (double)hits / (double)n;

        std::cout << "Num: " << n 
                  << " Result = " << d 
                  << " Error = " << (d - exact3()) << std::endl;
    }
};

int main() {
    std::cout << std::fixed << std::setprecision(16);

    Ex12 e;

    std::cout << std::endl;
    std::cout << "Monte Carlo with MSRNG" << std::endl;
    std::cout << std::endl;

    std::cout << "type1: 1D" << std::endl;
    for (int i = 1000; i <= 1000000; i *= 10) {
        // Init MSRNG: seed=1, a=16807, b=0, m=2^31-1 (0x7fffffff)
        e.InitGenerator(1, 16807, 0, 0x7fffffff);
        e.MonteCarlo1D(i);
    }

    std::cout << std::endl;
    std::cout << "type2: 2D" << std::endl;
    for (int i = 1000; i <= 1000000; i *= 10) {
        e.InitGenerator(1, 16807, 0, 0x7fffffff);
        e.MonteCarlo2D(i);
    }

    std::cout << std::endl;
    std::cout << "type2: 3D" << std::endl;
    for (int i = 1000; i <= 1000000; i *= 10) {
        e.InitGenerator(1, 16807, 0, 0x7fffffff);
        e.MonteCarlo3D(i);
    }

    return 0;
}