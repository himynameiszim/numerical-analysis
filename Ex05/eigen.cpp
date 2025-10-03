#include<bits/stdc++.h>

#include "eigen.hpp"
#include "helper.hpp"

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

template<typename T>
pair<vector<T>, vector<vector<T>>> EigenDeflatedPowerIteration(const vector<vector<T>>& A, const vector<T>& v0, double tol, int k){
    /*
    Deflated Power Iteration to find top k-largest, eigenvalues and their corresponding eigenvectors.
        After each convergence of Power Iteration, there is a small chance that the initial guess vector has no component parallel to the rest eigenvectors.
        Thus, keep multiplying does not do anything. So after each k-th convergence, we project xk out of vk, since eigenvectors of A are orthogonal,
        we obtain the (k+1)-th largest eigenvalue.
    (This requires A to be a symmetric matrix). One can try other deflation techniques (Householder transformation, ...) still to achieve the same result.
    */
    vector<T> eigenvalues;
    eigenvalues.push_back(EigenPowerIteration(A, v0).first);
    vector<vector<T>> eigenvectors;
    eigenvectors.push_back(EigenPowerIteration(A, v0).second);
    vector<vector<T>> B = A;
    int n = A.size();
    
    for(int l = 0; l < k; l++){
        T lambda = T(0);
        vector<T> v_l(n);
        random_device rd;
        std::mt19937 gen(rd());
        uniform_real_distribution<> dis(0.0, 1.0);
        for(int i = 0; i < n; ++i) v_l[i] = dis(gen);
        int maxIter = 50;
        while(maxIter--){
            vector<T> u_l = v_l;
            for(const auto& prev_vec : eigenvectors)
                u_l = vectorSubtract(u_l, scalerMultiply(prev_vec, dotProduct(v_l, prev_vec)));
            
            vector<T> w_l(n, T(0));
            for(int i = 0; i < n; i++){
                for(int j = 0; j < n; j++){
                    w_l[i] += B[i][j] * u_l[j];
                }
            }

            double norm = sqrt(dotProduct(w_l, w_l));
            for(int i = 0; i < n; i++){
                v_l[i] = w_l[i] / norm;
            }
            if(abs(norm - lambda) < tol) break;
            lambda = norm;
        }
        eigenvalues.push_back(lambda);
        eigenvectors.push_back(v_l);
    }  
    return make_pair(eigenvalues, transposeMatrix(eigenvectors));

}

template<typename T>
vector<T> EigenQRIteration(const vector<vector<T>>& A){
    /*
    Finding all eigenvalues from a matrix A, A does not have to be symmetric.
        Since the iteration takes action through orthogonalization (Gram-Schmidt), it is numerically stable.
        QR iteration is not guranteed to converge. Gershgorin circle theorem provides an error bound (thus, it gurantees convergence for diagonally dominant matrices).

    :return
        a vector, with all eigenvalues of A.
    */
    int n = A.size(), maxIter = 50;
    vector<vector<T>> Ak = A;
    for(int k = 0; k < maxIter; k++){
        auto [Q, R] = factorizeQR(Ak);
        Ak = matrixMultiply(R, Q); // note: this is multiplied in reverse order
    }
    vector<T> eigenvalues(n, T(0));
    for(int i = 0; i < n; i++){
        eigenvalues[i] = Ak[i][i];
    }
    return eigenvalues;
}

template pair<vector<double>, vector<vector<double>>> EigenJacobi(const vector<vector<double>>& A, double tol);
template pair<vector<float>, vector<vector<float>>> EigenJacobi(const vector<vector<float>>& A, double tol);
template pair<double, vector<double>> EigenPowerIteration(const vector<vector<double>>& A, const vector<double>& v0);
template pair<float, vector<float>> EigenPowerIteration(const vector<vector<float>>& A, const vector<float>& v0);
template pair<vector<double>, vector<vector<double>>> EigenDeflatedPowerIteration(const vector<vector<double>>& A, const vector<double>& v0, double tol, int k);
template pair<vector<float>, vector<vector<float>>> EigenDeflatedPowerIteration(const vector<vector<float>>& A, const vector<float>& v0, double tol, int k);
template vector<double> EigenQRIteration(const vector<vector<double>>& A);
template vector<float> EigenQRIteration(const vector<vector<float>>& A);