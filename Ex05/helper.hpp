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

#endif