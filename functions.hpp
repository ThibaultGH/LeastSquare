#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include "Mat.hpp"

#include <vector>
#include <iostream>
#include <assert.h>
#include <string>

using namespace std;

void decomp_LU(Mat& M);

vector<double> inv_triang_inf_LU(const Mat& M, const vector<double>& V);

vector<double> inv_triang_sup(const Mat& M, const vector<double>& v);

Mat transpose(const Mat& M);

vector<string> split(const string& str, const char* separator);

vector<double> VscaProd(const vector<double>& v, double l);

vector<double> VVSum(const vector<double>& v, const vector<double>& w);

double ProdL2(const vector<double>& v, const vector<double>& w);

double NormL2(const vector<double>& v);

vector<double> Col(const Mat& A, int i);

vector<double> Row(const Mat& A, int i);

Mat GramSchmidt(const Mat& A);

Mat Id(const int n);

Mat H(const vector<double>& v);





#endif
