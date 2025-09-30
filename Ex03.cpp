// FU11 - Exercise 03 - s1312004

#include<bits/stdc++.h>
#include<iomanip>

using namespace std;

const double tol = 1e-7;

template<typename T>
T f(T x){ return (x * x * x) - (3.0 * x * x) + (2.0 * x) - 0.2; }

template<typename T>
T df(T x){ return (3.0 * x * x) - (6.0 * x) + 2.0; }

template <typename T>
vector<pair<int, T>> solveBisection(T xPos, T xNeg, double tol){
    /*
    Traditional bisection method

    Convergence rate: Linear
    */
    if(f<double>(xPos) * f<double>(xNeg) > 0){
        cerr << "Root not bracketted" << endl;
        return {};
    }
    T x0;
    int iterCount = 0;
    vector<pair<int, T>> res;

    while(fabs(xPos - xNeg) >= tol){
        x0 = xNeg + (xPos - xNeg) / 2.0;

        if(f<double>(x0) == 0.0){
            iterCount++;
            res.push_back(make_pair(iterCount, x0));
            break;
        }

        (f<double>(x0) > 0) ? xPos = x0 : xNeg = x0;
        
        iterCount++;
        res.push_back(make_pair(iterCount, x0));
    }

    return res;
}

template<typename T>
vector<pair<int, T>> solveSecant(T x0, T x1, double tol){
    /*
    Traditional Secant method
    A slight improvement from bisection method by obtaining next iterate point by intersection with Ox instead of midpoint.

    Convergence rate: Superlinear.
    */
    T x2;
    int iterCount = 0;
    vector<pair<int, T>> res;
    (fabs(f<double>(x0)) < fabs(f<double>(x1))) ? swap(x0, x1) : (void)0;

    while(fabs(x1 - x0) >= tol){
        x2 = x0 - ( (f<double>(x0) * (x1 - x0)) / (f<double>(x1) - f<double>(x0)) );
        x0 = x1;
        x1 = x2;
        iterCount++;
        res.push_back(make_pair(iterCount, x1));
    }

    return res;
}

template<typename T>
vector<pair<int, T>> solveFalsePosition(T x0, T x1, double tol){
    /*
    Linear interpolation (False position method)
        Fixing the potential stucking problem of Secant method by keeping the root always bracketted in between the positive-negative interval.
    
    Convergence rate: Linear
    */
   T xNext, x = x0;
   int iterCount = 0;
   vector<pair<int, T>> res;
   (fabs(f<double>(x0)) < fabs(f<double>(x1))) ? swap(x0, x1) : (void)0;

   while(fabs(xNext - x0) >= tol){
        xNext = x;
        x = x0 - (f<double>(x0) * (x1 - x0)) / (f<double>(x1) - f<double>(x0));
        
        iterCount++;
        res.push_back(make_pair(iterCount, x));
        (f<double>(x0) < 0 || f<double>(x) < 0) ? x1 = x : x0 = x;
   }

    return res;
}

template<typename T>
vector<pair<int, T>> solveNewton(T x0, double tol){
    /*
    Traditional Newton's method

    Convergence rate: Quadratic
    */
    if(f<double>(x0) <= 1e-9 || df<double>(x0) <= 1e-9){
        cerr << "Initial guess or derivative too close to zero" << endl;
        return {};
    }
    T x1 = x0;
    int iterCount = 0;
    vector<pair<int, T>> res;

    while(true){
        x0 = x1;
        x1 = x0 - (f<double>(x0) / df<double>(x0));
        iterCount++;
        res.push_back(make_pair(iterCount, x1));

        if(fabs(x0 - x1) < tol)  break;        
    }

    return res;
}

int main(){
    cout << "f(x) = x^3 -3*x^2 + 2*x - 0.2 = 0" << endl;

    cout << "Bisection method" << endl;
    cout << "------------------------------" << endl;

    cout << setprecision(8);

    vector<pair<int, double>> taskA = solveBisection<double>(3.0, 1.5, tol);
    for(auto i : taskA)
        cout << "Iteration\t" << i.first << "\tx = \t" << scientific <<  i.second << endl;

    cout << "\nNewton's method" << endl;
    cout << "------------------------------" << endl;

    vector<pair<int, double>> taskB = solveNewton<double>(3.0, tol);
    for(auto i : taskB)
        cout << "Iteration\t" << i.first << "\tx = \t" << scientific << i.second << endl;

    cout << "\nSecant method" << endl;
    cout << "------------------------------" << endl;

    vector<pair<int, double>> taskC = solveSecant<double>(3.0, 1.5, tol);
    for(auto i : taskC)
        cout << "Iteration\t" << i.first << "\tx = \t" << scientific << i.second << endl;

    cout << "\nFalse position method" << endl;
    cout << "------------------------------" << endl;

    vector<pair<int, double>> taskD = solveFalsePosition<double>(3.0, 1.5, tol);
    for(auto i : taskC)
        cout << "Iteration\t" << i.first << "\tx = \t" << scientific << i.second << endl;

    return 0;
}