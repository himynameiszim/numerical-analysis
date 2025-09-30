#include<bits/stdc++.h>
#include<iomanip>

using namespace std;

// #define double long double

vector<double> gaussianEliminationSolve(vector<vector<double>> A, vector<double> b){
    /*
    Solving Ax = b where A is a square matrix, b is a column vector.

    :return
        the vector b containing the solutions of x
    */
   int n = A.size();

   // Forward substitution
   for(int i = 0; i < n; i++){

        // Partial pivoting: pick the pivotal row by choosing the row with the maximum absolute value element in the current column
        int maxRowIdx = i;
        for(int j = i + 1; j < n; j++)
            (abs(A[j][i])) > abs(A[maxRowIdx][j]) ? maxRowIdx = j : true;

        // Swapping rows, pushing the pivotal row up
        swap(A[i], A[maxRowIdx]);   swap(b[i], b[maxRowIdx]);

        if(abs(A[i][i]) < 1e-9)
            cout << "A is near singular, very ill-conditioned; continue anyway" << endl;;

        for(int j = i + 1; j < n; j++){
            int r = A[i][j] / A[i][i];
            for(int k = j; k < n; k++)
                A[j][k] -= r * A[i][k];

            b[j] -= r * b[i];
        }

        // A is now an upper-triangular matrix
   }

   // Backward substitution
   vector<double> x(n);
   for(int i = n - 1; i >= 0; i--){
        double sum = 0.0;
        for(int j = i + 1; j < n; j++){
            sum += A[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / A[i][i];
   }

   return x;
}

pair<vector<vector<double>>, vector<vector<double>>> factorizationLU(vector<vector<double>> A){
    /*
    LU Factorization pipeline. A = LU; L is lower-triangular matrix, U is a upper-triangular matrix, diag(U) = 1.

    :return
        a pair (L, U).
    */
   int n = A.size();
   vector<vector<double>> L(n, vector<double>(n, 0.0));
   vector<vector<double>> U(n, vector<double>(n, 0.0));

   for(int j = 0; j < n; j++){
        U[j][j] = 1.0;

        for(int i = j; i < n; i++){
            double sum = 0.0;
            for(int k = 0; k < j; k++){
                sum += L[i][k] * U[k][j];
            }
            L[i][j] = A[i][j] - sum;
        }

        for(int i = j + 1; i < n; i++){
            double sum = 0.0;
            for(int k = 0; k < j; k++){
                sum += L[j][k] * U[k][i];
            }
            U[j][i] = (A[j][i] - sum) / L[j][j];
        }
   }

    return make_pair(L, U);
}

void printSystem(vector<vector<double>> A, vector<double> b){
    cout << "A =    ";
    for(const vector<double>& i : A){
        for(const double& j : i){
            cout << "\t" << j << " ";
        }
        cout << endl;
    }
    cout << "\nb =   ";
    for(const double& j : b){
        cout << "\t" << j << endl;
    }
    cout << endl;
}

int main(){
    vector<vector<double>> A = {
        {2, -1, -2},
        {-4, 6, 3},
        {-4, -2, 8}
    };

    vector<double> b = {1, 5, 6};

    cout << "The system Ax = b :"  << endl;
    printSystem(A, b);

    vector<double> x = gaussianEliminationSolve(A, b);
    cout << "x =    ";
    for(auto i : x){
        cout << setprecision(8) << "\t" << i << endl;
    }
    
    cout << "A = LU" << endl;

    pair<vector<vector<double>>, vector<vector<double>>> LU = factorizationLU(A);
    vector<vector<double>> L = LU.first;
    vector<vector<double>> U = LU.second;

    cout << "L =    ";
    for(vector<double> i : factorizationLU(A).first){
        for(double j : i){
            cout << "\t" << fixed << setprecision(4) << j;
        }
        cout << endl;
    }
    cout << endl;
    cout << "U =    ";
    for(vector<double> k : factorizationLU(A).second){
        for(double t : k){
            cout << "\t" << fixed << setprecision(4) <<  t;
        }
        cout << endl;
    }

    return 0;
}