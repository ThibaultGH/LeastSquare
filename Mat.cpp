#include <vector>
#include <iostream>
#include <assert.h>
#include <math.h>
#include <fstream>
#include <string>
#include <typeinfo>

#include "Mat.hpp"
#include "functions.hpp"

using namespace std;

Mat::Mat(const int& nr, const int& nc)
{
  m_NbRow = nr;
  m_NbCol = nc;
  m_val = vector<double> (nr*nc);
}

Mat::Mat(const vector<double>& val, const int nr, const int nc)
{
  assert(val.size() == nr*nc);
  m_NbRow = nr;
  m_NbCol = nc;
  m_val = val;
}

Mat::Mat(const Mat& M) : m_val(M.m_val), m_NbRow(M.m_NbRow), m_NbCol(M.m_NbCol) {}

int Mat::Col() const {return m_NbCol;}

int Mat::Row() const {return m_NbRow;}

void Mat::resize(int nr, int nc)
{
  m_NbRow = nr;
  m_NbCol = nc;
  m_val.resize(nr*nc);
}

double Mat::operator() (int i, int j) const {return m_val[j+i*m_NbCol];}

double& Mat::operator() (int i, int j) {return m_val[j+i*m_NbCol];}

Mat & Mat::operator= (const Mat& M)
{
  if (this != &M){
    m_NbRow = M.Row();
    m_NbCol = M.Col();
    for(int i = 0; i < M.Row(); ++i){
      for(int j = 0; j < M.Col(); ++j){
        (*this)(i,j) = M(i,j);
      }
    }
  return *this;
  }
}


void Mat::Load (char* const filename)
{
  ifstream fichier(filename, ios::in);
  
  if (fichier){
    
    string line;

    vector<string> contenu;
    
    while(getline(fichier,line)){
      contenu.push_back(line);
    }
    
    vector<string> dim = split(contenu[1],"\t");
    
    int nr;
    int nc;
    
    sscanf(dim[0].c_str(), "%d", &nr);
    sscanf(dim[1].c_str(), "%d", &nc);
    
    vector<string> val = split(contenu[3],"   ");

    int k = val.size();
    
    vector<double> valu(k);

    for(int i = 0; i < k; ++i){
      sscanf(val[i].c_str(), "%lf", &valu[i]);
    }
    
    m_NbRow = nr;
    m_NbCol = nc;
    m_val.resize(k);
    
    for(int i = 0; i < k; ++i){
      m_val[i] = valu[i];
    }
    
    fichier.close();
  }
  else{
    cout << "error : cannot open file" << endl;
  }
}

ostream& operator<< (ostream& os, const Mat& M)
{
  for(int i = 0; i < M.m_NbRow; ++i){
    for(int j = 0; j < M.m_NbCol; ++j){
      os << M.m_val[j+i*M.m_NbCol] << "\t";
    }
     os << endl;
  }
  return os;
}

void Mat::Col(vector<double>& v, int i)
{
  assert(v.size() == m_NbRow);
  for(int j = 0; j < m_NbRow; ++j){
    (*this)(j,i) = v[j];}
  return ;
}

void Mat::Row(vector<double>& v, int i)
{
  assert(v.size() == m_NbCol);
  for(int j = 0; j < m_NbCol; ++j){
    (*this)(i,j) = v[j];}
  return ;
}

void MvProd(const Mat& A, const vector<double>& x, vector<double>& b)
{
  assert(A.m_NbCol == x.size());
  int n = A.m_NbRow;
  int m = A.m_NbCol;
  b.resize(n);
  for(int k = 0; k < n; ++k){b[k] = 0;}
  for(int i = 0; i < n; ++i){
    for(int j = 0; j < m; ++j){
      b[i] += A(i,j)*x[j];
    }
  }
}

void MMProd(const Mat& A, const Mat& B, Mat& C)
{
  assert(A.m_NbCol == B.m_NbRow);
  C.resize(A.m_NbRow,B.m_NbCol);
  for(int i = 0; i < C.m_NbRow; ++i){
    for(int j = 0; j < C.m_NbCol; ++j){
      C(i,j) = 0;
      for(int k = 0; k < A.m_NbCol; ++k){
	C(i,j) += A(i,k)*B(k,j);
      }
    }
  }
}

void Solve(const Mat& A,vector<double>& x,const vector<double>& b)
{
  assert(A.m_NbRow == b.size());
  if(A.m_NbRow == A.m_NbCol){
    int n = A.m_NbRow;
    Mat LU(A);
    decomp_LU(LU);
    vector<double> y = inv_triang_inf_LU(LU,b);
    x.resize(n);
    x = inv_triang_sup(LU,y);
  }
  else{
    cout << "error : matrix given is not square" << endl;
  }
}

void NormalSolve(const Mat& A,vector<double>& x,const vector<double>& b)
{
  assert(A.m_NbRow == b.size());
  int n = A.m_NbRow;
  int m = A.m_NbCol;
  Mat G(m,m);
  vector<double> D(m);
  MMProd(transpose(A),A,G);
  MvProd(transpose(A),b,D);
  Solve(G,x,D);
}

void GramSchmidtSolve(const Mat& A, vector<double>& x, const vector<double>& b)
{
  Mat R(A.Col(),A.Col());
  Mat GS(A.Row(),A.Col());
  vector<double> w = Col(A,0);
  vector<double> q = VscaProd(w,1./NormL2(w));
  GS.Col(q,0);
  R(0,0) = NormL2(w);
  for(int k = 1; k < A.Col(); ++k){
    w = Col(A,k);
    for(int j = k-1; j >= 0; --j){
      w = VVSum(w,VscaProd(VscaProd(Col(GS,j),ProdL2(Col(GS,j),Col(A,k))),-1));
    }
    R(k,k) = NormL2(w);
    q = VscaProd(w,1./NormL2(w));
    GS.Col(q,k);
  }
  
  for(int k = 0; k < A.Col(); ++k){
    for(int j = 0; j < k; ++j){
      R(j,k) = ProdL2(Col(GS,j),Col(A,k));
    }
  }
  
  Mat invR(A.Col(),A.Col());
  for(int i = 0; i < A.Col(); ++i){
    vector<double> e(A.Col());
    e[i] = 1;
    vector<double> c = inv_triang_sup(R,e);
    invR.Col(c,i);
  }
  Mat Q(A.Row(),A.Col());
  MMProd(A,invR,Q);

  vector<double> b1(A.Col());
  MvProd(transpose(Q),b,b1);

  MvProd(invR,b1,x);
  
  //x = inv_triang_sup(R,b1);
}

void HouseholderSolve(const Mat& A, vector<double>& x, const vector<double>& b)
{
  int n = A.Row();
  int m = A.Col();
  Mat I = Id(n);
  vector<double> u = VVSum(Col(A,0),VscaProd(Col(I,0),(-1)*NormL2(Col(A,0))));
  Mat S = H(u);
  Mat R(n,m);
  MMProd(S,A,R);
  Mat Q = S;
  for(int p = 0; p < m-1; ++p){
    vector<double> w(n);
    for(int k = p+1; k < n; ++k){w[k] = Col(R,p+1)[k];}
    u = VVSum(w,VscaProd(Col(I,p+1),(-1)*NormL2(w)));    
    S = H(u);
    Mat temp1(n,n);
    MMProd(Q,S,temp1);
    Q = temp1;
    Mat temp2(n,m);
    MMProd(S,R,temp2);
    R = temp2;
  }

  Mat r(m,m);
  Mat q(n,m);
  for(int k = 0; k < m; ++k){
    vector<double> t1 = Row(R,k);
    r.Row(t1,k);
  }
  
  for(int k = 0; k < m; ++k){
    vector<double> t2 = Col(Q,k);
    q.Col(t2,k);
  }
  
  vector<double> b1(A.Col());

  MvProd(transpose(q),b,b1);

  x = inv_triang_sup(r,b1);
}
