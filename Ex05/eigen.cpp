#include<bits/stdc++.h>

#include "eigen.hpp"

using namespace std;

template<typename T>
pair<vector<T>, vector<vector<T>>> EigenJacobi(const vector<vector<T>>& A, double tol){
    /*
    Jacobi diagonalization algorithm. Linear convergence.

    :return
        a pair of eigenvalues and a matrix with eigenvectors embedded in it columns.
    */

    int n = A.size();
    (n == 0 || A[0].size() != n) ? throw invalid_argument("Matrix must be squared. \n") : true;

    vector<vector<T>> currentA = A;
    vector<T> eigenvals(n);
    vector<vector<T>> eigenvecs(n, vector<T>(n, 0.0));

    for(int i = 0; i < n; i++){
        eigenvecs[i][i] = 1.0;
    }

    int iterCount = 0, maxIter = 100;
    while(iterCount++ < maxIter){
        T maxVal = 0.0;
        int p = 0, q = 1;
        for(int i = 0; i < n; i++){
            for(int j = i + 1; j < n; j++){
                if(abs(currentA[i][j]) > maxVal){
                    maxVal = abs(currentA[i][j]);
                    p = i;
                    q = j;
                }
            }
        }

        if(maxVal < tol){
            break;
        }

        // Build the hideous Jacobi transformation matrix
        T a_pq = currentA[p][q];
        T a_pp = currentA[p][p];
        T a_qq = currentA[q][q];

        T thetaVal = 0.5 * atan2(2 * a_pq, a_qq - a_pp);
        T cosTheta = cos(thetaVal);
        T sinTheta = sin(thetaVal);

        T a_ppNew = (cosTheta * cosTheta * a_pp) + (sinTheta * sinTheta * a_qq) - (2 * cosTheta * sinTheta * a_pq);
        T a_qqNew = (sinTheta * sinTheta * a_pp) + (cosTheta * cosTheta * a_qq) + (2 * cosTheta * sinTheta * a_pq);

        currentA[p][p] = a_ppNew;
        currentA[q][q] = a_qqNew;
        currentA[p][q] = currentA[q][p] = 0.0;

        for(int i = 0; i < n; i++){
            if(i != p && i != q){
                T a_ip = currentA[i][p];
                T a_iq = currentA[i][q];

                currentA[i][p] = currentA[p][i] = (cosTheta * a_ip) - (sinTheta * a_iq);
                currentA[i][q] = currentA[q][i] = (sinTheta * a_ip) + (cosTheta * a_iq);
            }
        }

        for(int i = 0; i < n; i++){
            T v_ip = eigenvecs[i][p];
            T v_iq = eigenvecs[i][q];

            eigenvecs[i][p] = (cosTheta * v_ip) - (sinTheta * v_iq);
            eigenvecs[i][q] = (sinTheta * v_ip) + (cosTheta * v_iq);
        }
    }

    for(int i = 0; i < n; i++){
        eigenvals[i] = currentA[i][i];
    }

    return make_pair(eigenvals, eigenvecs);
}

template<typename T>
pair<T, vector<T>> EigenPowerIteration(const vector<vector<T>>& A, const vector<T>& v0){
    /*
    Find a single eigenvalue using the Power Iteration method.

    :return
        one eigenvalue of the matrix A.
    */

    //normalizing init vector
    int iterCount = 0, maxIter = 100, n = A.size();
    T norm = 0.0, sum = 0.0, eigenvalue1 = 0.0;

    for(T val : v0){
        sum += (val * val);
    }
    vector<T> v = v0;
    for(int i = 0; i < n; i++){
        v[i] /= sqrt(sum);
    }

    (n == 0 || A[0].size() != n || v0.size() != n) ? throw invalid_argument("Matrix/vector in wrong dimension. \n") : true;

    while(iterCount++ < maxIter){
        vector<T> Av(n, 0.0);
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                Av[i] += A[i][j] * v[j];
            }
        }

        // normalize vector by L2 norm
        norm = 0.0;
        for(double val : Av)
            norm += (val * val);
        norm = sqrt(norm);
        
        for(int i = 0; i < n; i++){
            v[i] = Av[i] / norm;
        }

        if(abs(norm - eigenvalue1) < tol){
            eigenvalue1 = norm;
            break;
        }

        eigenvalue1 = norm;
    }   
    return make_pair(eigenvalue1, v);
}

template pair<vector<double>, vector<vector<double>>> EigenJacobi(const vector<vector<double>>& A, double tol);
template pair<vector<float>, vector<vector<float>>> EigenJacobi(const vector<vector<float>>& A, double tol);
template pair<double, vector<double>> EigenPowerIteration(const vector<vector<double>>& A, const vector<double>& v0);
template pair<float, vector<float>> EigenPowerIteration(const vector<vector<float>>& A, const vector<float>& v0);