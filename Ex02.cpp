// FU11 - Exercise 02 - s1312004

#include<bits/stdc++.h>
#include<iomanip>

using namespace std;

template <typename T>
pair<T, T> Ex01(T a, T b, T c){
    T delta = b * b - 4 * (a * c);
    T x1 = ( (-1 * b) + sqrt(delta) ) / (2 * a);
    T x2 = ( (-1 * b) - sqrt(delta) ) / (2 * a);

    return make_pair(x1, x2);
}

template <typename T>
pair<T, T> Ex01Improved(T a, T b, T c){
    /*
    Since we can pre-calculate delta, we know that the process of calculating x2 involves
        the subtraction of two nearly equal positive numbers in the numerator which yields huge relative error, causing a bad approximation.
    Thus, we apply Vieta's. Since x1.x2 = c/a -> x2 = c / (a.x1)
    */
    T delta = b * b - 4 * (a * c);

    T x1 = ( (-1 * b) + sqrt(delta) ) / (2 * a);
    T x2 = c / ( a * x1 );

    return make_pair(x1, x2);
}

const double PI = M_PI;
pair<double, double> Ex02(double x, char* taskname){
    /*
    a. Each log calculation will truncate the solution, causing precision problem.
        Reducing the formula down to one log calculation will reduce loss of significance.
    b. For extremely large x, numerator will be truncated to 0, losing numerical precision.
        By rationalizing the process of computation, we avoid round-off errors.
    c. Reduce computation complexity and precision 
        by reducing two truncated square operations, one subtraction down to only evaluation.
    */
    if(taskname == "a")   return {log(x + 1.0) - log(x), log(1.0 + 1.0 / x)};
    if(taskname == "b")   return {sqrt((x * x) + 1.0) - x, 1.0 / (sqrt((x * x) + 1.0) + 1.0)};
    if(taskname == "c")   return {(cos(x) * cos(x)) - (sin(x) * sin(x)), cos(2.0 * x)}; // one need to convert x to rad first.
}

int main(){
    cout << "Task 1. Quadratic equation ax^2+bx+c=0" << endl;
    cout << "--------------------" << endl;
    cout << "a = 1.0 \t b = -8000.0 \t c = 1.0" << endl;

    cout << "----------Traditional technique----------" << endl;
    cout << setprecision(16);
    cout << "Sinlge precision x1 = " << scientific << Ex01<float>(1.0f, -8000.0f, 1.0f).first << "; x2 = " << Ex01<float>(1.0f, -8000.0f, 1.0f).second << endl;
    cout << "Double precision x1 = " << scientific << Ex01<double>(1.0, -8000.0, 1.0).first << "; x2 = " << Ex01<double>(1.0, -8000.0, 1.0).second << endl;
    
    cout << "----------Improved for single precision----------" << endl;
    cout << "Sinlge precision x1 = " << scientific << Ex01Improved<float>(1.0f, -8000.0f, 1.0f).first << "; x2 = " << Ex01Improved<float>(1.0f, -8000.0f, 1.0f).second << endl;
    
    cout << "\nTask 2. Cancellation" << endl;
    cout << "--------------------" << endl;

    cout << "ln(x+1) - ln(x), \t x = 1.0e12 : " << Ex02(1e12, "a").first << endl;
    cout << "ln(1 + 1/x), \t x = 1.0e12 : " << Ex02(1e12, "a").second << endl;

    cout << "--------------------" << endl;

    cout << "sqrt(x^2 + 1) - x, \t x = 1.0e9 : " << Ex02(1e9, "b").first << endl;
    cout << "1 / (sqrt(x^2 + 1) + x, \t x = 1.0e9 : " << Ex02(1e9, "b").second << endl;

    cout << "--------------------" << endl;
    double x_rad = 44.999999999999 * PI / 180.0;
    cout << "cos^2(x) - sin^2(x), \t x = 44.999999999999 : " << Ex02(x_rad, "c").first << endl;
    cout << "cos(2x), \t x = 44.999999999999 : " << Ex02(x_rad, "c").second << endl;

    return 0;
}