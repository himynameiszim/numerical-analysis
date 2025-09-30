#include<bits/stdc++.h>
#include<iomanip>

using namespace std;

const double tol = 1e-10;

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
void printMatrix(const vector<vector<T>>& A){
    for(const auto& row : A){
        for(const auto& i : row){
            cout << setprecision(5) << i << "  ";
        }
        cout << endl;
    }
    cout << endl;
}

pair<vector<double>, vector<vector<double>>> EigenJacobi(const vector<vector<double>>& A, double tol){
    /*
    Jacobi diagonalization algorithm. Linear convergence.

    :return
        a pair of eigenvalues and a matrix with eigenvectors embedded in it columns.
    */

    int n = A.size();
    (n == 0 || A[0].size() != n) ? throw invalid_argument("Matrix must be squared. \n") : true;

    vector<vector<double>> currentA = A;
    vector<double> eigenvals(n);
    vector<vector<double>> eigenvecs(n, vector<double>(n, 0.0));

    for(int i = 0; i < n; i++){
        eigenvecs[i][i] = 1.0;
    }

    int iterCount = 0, maxIter = 100;
    while(iterCount++ < maxIter){
        double maxVal = 0.0;
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
        double a_pq = currentA[p][q];
        double a_pp = currentA[p][p];
        double a_qq = currentA[q][q];

        double thetaVal = 0.5 * atan2(2 * a_pq, a_qq - a_pp);
        double cosTheta = cos(thetaVal);
        double sinTheta = sin(thetaVal);

        double a_ppNew = (cosTheta * cosTheta * a_pp) + (sinTheta * sinTheta * a_qq) - (2 * cosTheta * sinTheta * a_pq);
        double a_qqNew = (sinTheta * sinTheta * a_pp) + (cosTheta * cosTheta * a_qq) + (2 * cosTheta * sinTheta * a_pq);

        currentA[p][p] = a_ppNew;
        currentA[q][q] = a_qqNew;
        currentA[p][q] = currentA[q][p] = 0.0;

        for(int i = 0; i < n; i++){
            if(i != p && i != q){
                double a_ip = currentA[i][p];
                double a_iq = currentA[i][q];

                currentA[i][p] = currentA[p][i] = (cosTheta * a_ip) - (sinTheta * a_iq);
                currentA[i][q] = currentA[q][i] = (sinTheta * a_ip) + (cosTheta * a_iq);
            }
        }

        for(int i = 0; i < n; i++){
            double v_ip = eigenvecs[i][p];
            double v_iq = eigenvecs[i][q];

            eigenvecs[i][p] = (cosTheta * v_ip) - (sinTheta * v_iq);
            eigenvecs[i][q] = (sinTheta * v_ip) + (cosTheta * v_iq);
        }
    }

    for(int i = 0; i < n; i++){
        eigenvals[i] = currentA[i][i];
    }

    return make_pair(eigenvals, eigenvecs);
}

pair<double, vector<double>> EigenPowerIteration(const vector<vector<double>>& A, const vector<double>& v0){
    /*
    Find a single eigenvalue using the Power Iteration method.

    :return
        one eigenvalue of the matrix A.
    */

    //normalizing init vector
    int iterCount = 0, maxIter = 100, n = A.size();
    double norm = 0.0, sum = 0.0, eigenvalue1 = 0.0;

    for(double val : v0){
        sum += (val * val);
    }
    vector<double> v = v0;
    for(int i = 0; i < n; i++){
        v[i] /= sqrt(sum);
    }

    (n == 0 || A[0].size() != n || v0.size() != n) ? throw invalid_argument("Matrix/vector in wrong dimension. \n") : true;

    while(iterCount++ < maxIter){
        vector<double> Av(n, 0.0);
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

void printVector(const vector<double>& b){
    /*
    Print a vector
    */
   for(auto& i : b){
        cout << setprecision(5) << "( " << i << " )\n";
   }
   cout << endl;
}

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
    pair<vector<double>, vector<vector<double>>> Jacobi = EigenJacobi(A, tol);
    
    printVector(Jacobi.first);
    printMatrix(Jacobi.second);

    // With power iteration
    vector<double> v0 = {1, 2, 4, -1};
    pair<double, vector<double>> PowerIteration = EigenPowerIteration(A, v0);

    cout << "Lambda1 : " << PowerIteration.first << endl;
    printVector(PowerIteration.second);

    return 0;
}