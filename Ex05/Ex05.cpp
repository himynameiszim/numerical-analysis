#include<bits/stdc++.h>
#include<iomanip>

#include "helper.hpp"
#include "eigen.hpp"

using namespace std;

int main(){
    // Problem 1: Calculating the inverse of a matrix
        vector<vector<double>> A = {
            {4, 2, 2, 1},
            {2, -3, 1, 1},
            {2, 1, 3, 1},
            {1, 1, 1, 2}
        };
        vector<vector<double>> invA = invertMatrix<double>(A);

        printMatrix<double>(invA);

        vector<vector<double>> I = matrixMultiply<double>(A, invA);

        printMatrix<double>(I);

    cout << "------------------" << endl;

    // Problem 2: Finding eigenvalues, eigenvectors

        // With Jacobi diagonalization method
        pair<vector<double>, vector<vector<double>>> Jacobi = EigenJacobi(A, tol);
        
        printVector(Jacobi.first);
        printMatrix(Jacobi.second);

        // With power iteration
        vector<double> v0 = {1, 2, 4, -1};
        pair<double, vector<double>> PowerIteration = EigenPowerIteration(A, v0);

        cout << "Lambda1 : " << PowerIteration.first << endl;
        printVector(PowerIteration.second);

        // With Deflated power iteration
        pair<vector<double>, vector<vector<double>>> DeflatedPowerIteration = EigenDeflatedPowerIteration(A, v0, tol, 3);

        printVector(DeflatedPowerIteration.first);
        printMatrix(DeflatedPowerIteration.second);

    return 0;
}