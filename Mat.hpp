#ifndef MAT_HPP
#define MAT_HPP

#include <vector>
#include <iostream>


using namespace std;


class Mat {
  
  friend void MvProd(const Mat& A, const vector<double>& x, vector<double>& b);

  friend void MMProd(const Mat& A, const Mat& B, Mat& C);

  friend ostream& operator<< (ostream& os, const Mat& M);

  friend void Solve(const Mat& A, vector<double>& x, const vector<double>& b);

  friend void NormalSolve(const Mat& A, vector<double>& x, const vector<double>& b);

  friend void GramSchmidtSolve(const Mat& A, vector<double>& x, const vector<double>& b);

  friend void HouseholderSolve(const Mat& A, vector<double>& x, const vector<double>& b);

private :

  int m_NbRow;

  int m_NbCol;

  vector<double> m_val;

public :
  
  Mat(const int& nr = 0, const int& nc = 0);

  Mat(const vector<double>& val, const int nr, const int nc);

  Mat(const Mat& M);

  int Col() const;

  int Row() const;
  
  void Col(vector<double>& v, int i);

  void Row(vector<double>& v, int i);

  void resize(int nr, int nc);

  double operator() (int i, int j) const;

  double& operator() (int i, int j);
  
  Mat & operator= (const Mat& M);

  void Load (char* const filename);


};


#endif
