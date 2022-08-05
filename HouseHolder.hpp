#ifndef HOUSEHOLDER_HPP
#define HOUSEHOLDER_HPP

#include <vector>
#include <iostream>
#include <assert.h>

using namespace std;

class HH{

private :

  vector<double> m_vec;

public :

  HH(int n = 0);

  HH(const vector<double>& vec);

  HH(const HH& H);

  int size ();

  double& operator() (int i, int j);
  
  double operator() (int i, int j) const;

  vector<double> operator* (const vector<double>& vec);

};


#endif
