#ifndef EIGEN_HPP
#define EIGEN_HPP

#include<bits/stdc++.h>

using namespace std;

const double tol = 1e-10;

template<typename T>
pair<vector<T>, vector<vector<T>>> EigenJacobi(const vector<vector<T>>& A, double tol);
template<typename T>
pair<T, vector<T>> EigenPowerIteration(const vector<vector<T>>& A, const vector<T>& v0);
template<typename T>
pair<vector<T>, vector<vector<T>>> EigenDeflatedPowerIteration(const vector<vector<T>>& A, const vector<T>& v0, double tol, int k);
template<typename T>
vector<T> EigenQRIteration(const vector<vector<T>>& A);

#endif