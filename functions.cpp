#include "functions.hpp"
#include "Mat.hpp"

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

using namespace std;

void decomp_LU(Mat& M)
{
  assert(M.Col() == M.Row());
  int n = M.Row();
  for(int k = 0; k < n; ++k){
    for(int i = k+1; i < n; ++i){
      M(i,k) = M(i,k)/M(k,k);
      for(int j = k+1; j < n; ++j){
	M(i,j) -= M(i,k)*M(k,j);
      }
    }
  }
}

vector<double> inv_triang_inf_LU(const Mat& M, const vector<double>& v)
{
  assert(M.Row() == v.size());
  int n = v.size();
  vector<double> x(n);
  for (int i = 0; i < n; ++i){
    for (int j = 0; j < i; ++j){
      x[i] += M(i,j)*x[j];
    }
    x[i] = (v[i]-x[i]);
  }
  return x;
}

vector<double> inv_triang_sup(const Mat& M, const vector<double>& v)
{
  assert(M.Row() == v.size());
  int n = v.size();
  vector<double> x(n);
  for (int i = n-1; i >= 0; --i){
    for (int j = n-1; j > i; --j){
      x[i] += M(i,j)*x[j];
    }
    x[i] = (1/M(i,i))*(v[i]-x[i]);
  }
  return x;
}

Mat transpose(const Mat& M)
{
  int n = M.Row();
  int m = M.Col();
  Mat result(m,n);
  for(int i = 0; i < m; ++i){
    for(int j = 0; j < n; ++j){
      result(i,j) = M(j,i);
    }
  }
  return result;
}


vector<string> split(const string& str, const char* separator)
{
  string stri = str;
  int n = stri.length();
  string word = "";
  vector<string> splited;
  if(stri[n-1] != *separator){
    stri += *separator;}
  n = stri.length();
  for(int i = 0; i < n; ++i){      
    if(stri[i] != *separator){
      word += stri[i];
    }
    else{
      if(word != ""){
	splited.push_back(word);
	word = "";}
    }
  }
  return splited;
}

vector<double> VscaProd(const vector<double>& v, double l)
{
  vector<double> result(v.size());
  for(int i = 0; i < v.size(); ++i){
    result[i] = v[i]*l;}
  return result;
}

vector<double> VVSum(const vector<double>& v, const vector<double>& w)
{
  assert(v.size() == w.size());
  vector<double> result(v.size());
  for(int i = 0; i < v.size(); ++i){
    result[i] = v[i] + w[i];}
  return result;
}

double ProdL2(const vector<double>& v, const vector<double>& w)
{
  assert(v.size() == w.size());
  double result = 0;
  for(int i = 0; i < v.size(); ++i){
    result += v[i]*w[i];}
  return result;
}

double NormL2(const vector<double>& v)
{
  double result = sqrt(ProdL2(v,v));
  return result;
}

vector<double> Col(const Mat& A, int i)
{
  vector<double> col(A.Row());
  for(int j = 0; j < A.Row(); ++j){
    col[j] = A(j,i);}
  return col;
}

vector<double> Row(const Mat& A, int i)
{
  vector<double> row(A.Col());
  for(int j = 0; j < A.Col(); ++j){
    row[j] = A(i,j);}
  return row;
}

Mat GramSchmidt(const Mat& A)
{
  Mat result(A.Row(),A.Col());
  vector<double> w = Col(A,0);
  vector<double> q = VscaProd(w,1./NormL2(w));
  result.Col(q,0);  
  for(int k = 1; k < A.Col(); ++k){
    w = Col(A,k);
    for(int j = k-1; j >= 0; --j){
      /*double x1 = ProdL2(Col(result,j),Col(A,k));
      vector<double> vec = VscaProd(Col(result,j),x1);
      vec = VscaProd(vec,-1);
      w = VVSum(w,vec);
      */
      w = VVSum(w,VscaProd(VscaProd(Col(result,j),ProdL2(Col(result,j),Col(A,k))),-1));
    }
    q = VscaProd(w,1./NormL2(w));
    result.Col(q,k);
  }
  return result;
}


Mat Id(const int n)
{
  Mat result(n,n);
  for(int i = 0; i < n; ++i){
    result(i,i) = 1;}
  return result;
}

Mat H(const vector<double>& v)
{
  int n = v.size();
  Mat H = Id(n);
  double N = NormL2(v);
  for(int i = 0; i < n; ++i){
    for(int j = 0; j < n; ++j){
      H(i,j) -= 2*v[i]*v[j]/(N*N);}
  }
  return H;
}
