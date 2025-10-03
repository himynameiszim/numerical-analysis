#ifndef HELPER_HPP
#define HELPER_HPP

#include<bits/stdc++.h>

using namespace std;

template<typename T>
vector<vector<T>> invertMatrix(const vector<vector<T>>& A);
template<typename T>
vector<vector<T>> matrixMultiply(const vector<vector<T>>& A, const vector<vector<T>>& B);
template<typename T>
void printMatrix(const vector<vector<T>>& A);
template<typename T>
void printVector(const vector<T>& b);
template<typename T>
T dotProduct(const vector<T>& a, const vector<T>& b);
template<typename T>
vector<T> vectorSubtract(const vector<T>& a, const vector<T>& b);
template<typename T>
vector<T> scalerMultiply(const vector<T>& a, double s);
template<typename T>
vector<vector<T>> transposeMatrix(const vector<vector<T>>& A);
template<typename T>
pair<vector<vector<T>>, vector<vector<T>>> factorizeQR(const vector<vector<T>>& A);

#endif