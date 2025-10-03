#include<bits/stdc++.h>

#include "helper.hpp"

using namespace std;

template<typename T>
vector<vector<T>> invertMatrix(const vector<vector<T>>& A){
    /*
    Finding the inverse of a matrix using Gaussian Elimination of the augmented matrix [A | I]

    :return
        the inverse of the input matrix
    */
   int n = A.size();
   (n == 0 || A[0].size() != n) ? throw invalid_argument("Matrix not square. \n") : false;

    // create the augemented matrix [A | I]
   vector<vector<T>> augmented(n, vector<T>(2 * n));
   for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            augmented[i][j] = A[i][j];
        }
        augmented[i][i+n] = T(1);
   }

   // GE here
   for(int i = 0; i < n; i++){
        // pivotting
        int maxIdxRow = i;  
        for(int j = i; j < n; j++)
            (abs(augmented[j][i]) > abs(augmented[maxIdxRow][i])) ? maxIdxRow = j : true;

        swap(augmented[i], augmented[maxIdxRow]);
        T pivot = augmented[i][i];
        if(abs(pivot) < 1e-9){
            throw runtime_error("Matrix is near singular, very ill-conditioned.");
        }

        for(int j = i; j < 2 * n; j++)
            augmented[i][j] /= pivot; 
        
        for(int j = 0; j < n; j++){
            if(i != j){
                T factor = augmented[j][i];
                for(int k = i; k < 2 * n; k++){
                    augmented[j][k] -= factor * augmented[i][k];
                }
            }
        }
   }

   vector<vector<T>> inverseA(n, vector<T>(n));
   for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            inverseA[i][j] = augmented[i][j+n];
        }
   }

   return inverseA;
}

template<typename T>
vector<vector<T>> matrixMultiply(const vector<vector<T>>& A, const vector<vector<T>>& B){
    /*
    Multiply two matrices together A \in R^nxm; B \in B \in R^mxn

    :return
        The product of the two input matrices.
    */

   (A[0].size() != B[0].size()) ? throw invalid_argument("Invalid dimension") : true;
   int n = A.size(), m = A[0].size();
   vector<vector<T>> result(n, vector<T>(n));
   int k = 0;
   for(int i = 0; i < n; i++){
        T sum = T(0);
        for(int j = 0; j < m; j++){
            sum += A[i][j] * B[j][i];
        }
        result[i][k] = sum;
        k++;
   }

   return result;
} 

template<typename T>
T dotProduct(const vector<T>& a, const vector<T>& b){
    /*
    Dot product between 2 vectors.
    */
    int n = a.size();
    T sum = T(0);
    for(int i = 0; i < n; i++){
        sum += a[i] * b[i];
    }
    return sum;
}

template<typename T>
vector<T> vectorSubtract(const vector<T>& a, const vector<T>& b){
    /*
    Subtract 2 vectors.
    */
    int n = a.size();
    vector<T> result(n, T(0));
    for(int i = 0; i < n; i++){
            result[i] = a[i] - b[i];
    }
    return result;
}

template<typename T>
vector<T> scalerMultiply(const vector<T>& a, double s){
    /*
    Multiply scaler by vector.
    */
    int n = a.size();
    vector<T> result(n, T(0));
    for(int i = 0; i < n; i++){
        result[i] = s * a[i];
    }
    return result;
}

template<typename T>
vector<vector<T>> transposeMatrix(const vector<vector<T>>& A){
    /*
    Invert a matrix
    */
    int n = A.size();
    int m = A[0].size();
    vector<vector<T>> B = A;

    for(int i = 0; i < n; ++i){
        for(int j = 0; j < m; ++j)
            B[i][j] = A[j][i];
    }
    return B;
}

template<typename T>
void printMatrix(const vector<vector<T>>& A){
    /*
    Print a mtrix
    */
    for(const auto& row : A){
        for(const auto& i : row){
            cout << setprecision(5) << i << "  ";
        }
        cout << endl;
    }
    cout << endl;
}

template<typename T>
void printVector(const vector<T>& b){
    /*
    Print a vector
    */
   for(auto& i : b){
        cout << setprecision(5) << "( " << i << " )\n";
   }
   cout << endl;
}

template vector<vector<double>> invertMatrix(const vector<vector<double>>& A);
template vector<vector<float>> invertMatrix(const vector<vector<float>>& A);
template vector<vector<double>> matrixMultiply(const vector<vector<double>>& A, const vector<vector<double>>& B);
template vector<vector<float>> matrixMultiply(const vector<vector<float>>& A, const vector<vector<float>>& B);
template void printMatrix(const vector<vector<double>>& A);
template void printMatrix(const vector<vector<float>>& A);
template void printVector(const vector<double>& b);
template void printVector(const vector<float>& b);
template double dotProduct(const vector<double>& a, const vector<double>& b);
template float dotProduct(const vector<float>& a, const vector<float>& b);
template vector<double> vectorSubtract(const vector<double>& a, const vector<double>& b);
template vector<float> vectorSubtract(const vector<float>& a, const vector<float>& b);
template vector<double> scalerMultiply(const vector<double>& a, double s);
template vector<float> scalerMultiply(const vector<float>& a, double s);
template vector<vector<double>> transposeMatrix(const vector<vector<double>>& A);
template vector<vector<float>> transposeMatrix(const vector<vector<float>>& A);

