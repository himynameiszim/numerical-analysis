#include<iostream>
#include<cmath>
#include<vector>
#include<iomanip>

template<typename T>
T f(T x){
    return std::sin(x);
}

template<typename T>
T df_fwdDiff(T x, T h){
    return (f(x + h) - f(x)) / h;
}

template<typename T>
T df_centDiff(T x, T h){
    return (f(x + h) - f(x - h)) / (2.0 * h);
}

template<typename T>
T df_exact(T x){
    return std::cos(x);
}

void Task1(){
    double h = .025 * M_PI;

    std::vector<double> p;
    for(int i = 0; i <= 5; ++i)
        p.push_back(i * .1 * M_PI);
    
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Numerical differentiation\n" << std::endl;

    std::cout << "Forward difference" << std::endl;

    for(auto x : p){
        double dfx_fwdDiff = df_fwdDiff<double>(x, h);
        double dfx_exact = df_exact<double>(x);
        double error = dfx_fwdDiff - dfx_exact;

        std::cout << "x = " << std::setprecision(2) <<  x/M_PI 
                  << " PI  Deriv = " << std::setprecision(6)
                  << dfx_fwdDiff << " Error = " << error
                  << std::endl;
    }

    std::cout << "\nCentral difference" << std::endl;

    for(auto x : p){
        double dfx_centDiff = df_centDiff<double>(x, h);
        double dfx_exact = df_exact<double>(x);
        double error = dfx_centDiff - dfx_exact;

        std::cout << "x = " << std::setprecision(2) <<  x/M_PI 
                  << " PI  Deriv = " << std::setprecision(6)
                  << dfx_centDiff << " Error = " << error
                  << std::endl;
    }
}

template<typename T>
T g(T x){
    return 1.0 / (x * x  + 1.0);
}

void Task2(){
    double a = 0.0;
    double b = 1.0;
    int n = 24;
    double h = (b - a) / n;
    double int_exact = M_PI / 4.0;

    std::cout << "Numerical integration\n" << "I = Integral(1/(x^2+1)), limits: [0,1]\n" << std::endl;

    std::cout << std::fixed << std::setprecision(9);

    std::cout << "Trapezoidal rule, nsub = 24" << std::endl;
    double trap_sum = g<double>(a) + g<double>(b);

    for(int i = 1; i < n; ++i){
        double x = a + i * h;
        trap_sum += 2.0 * g<double>(x);
    }
    double trap_int = (h / 2.0) * trap_sum;
    double trap_error = trap_int - int_exact;
    double trap_relError = (trap_error) / int_exact;

    std::cout << "I = " << trap_int
                << " Error = " << std::scientific << std::setprecision(4) 
                << trap_error << " R = " << trap_relError * 100 << " %"
                << std::endl;

    std::cout << "\nSimpson's rule, nsub = 24" << std::endl;
    double simp_sum = g<double>(a) + g<double>(b);
    for(int i = 1; i < n; ++i){
        double x = a + i * h;
        if(i % 2 == 0)
            simp_sum += 2.0 * g<double>(x);
        else
            simp_sum += 4.0 * g<double>(x);
    }
    double simp_int = (h / 3.0) * simp_sum;
    double simp_error = simp_int - int_exact;
    double simp_relError = (simp_error) / int_exact;

    std::cout << "I = " << std::fixed << std::setprecision(9) << simp_int
                << " Error = " << std::scientific << std::setprecision(4) 
                << simp_error << " R = " << simp_relError * 100 << " %"
                << std::endl;

}


int main(){
    Task1();
    Task2();
    return EXIT_SUCCESS;

}