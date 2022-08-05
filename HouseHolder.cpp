#include <vector>
#include <assert.h>

#include "HouseHolder.hpp"
#include "functions.hpp"

using namespace std;

HH::HH(int n){
  m_vec = vector<double> (n,0);
}

HH::HH(const vector<double>& vec){
  int n = vec.size();
  m_vec.resize(n);
  for(int i = 0; i < n; ++i){
    m_vec[i] = vec[i];
  }
}

HH::HH(const HH& H) : m_vec(H.m_vec) {}

int HH::size () {return m_vec.size();}

double& HH::operator() (int i, int j)
{
  double result;
  if(i == j){
    result  = 1.0-2*m_vec[i]*m_vec[i]/ProdL2(m_vec,m_vec);
  }
  else{
    result  = -2*m_vec[i]*m_vec[j]/ProdL2(m_vec,m_vec);
  }
  return result;
}

double HH::operator() (int i, int j) const
{
  double result;
  if(i == j){
    result  = 1-2*m_vec[i]*m_vec[i]/ProdL2(m_vec,m_vec);
  }
  else{
    result  = -2*m_vec[i]*m_vec[j]/ProdL2(m_vec,m_vec);
  }
  return result;
}

vector<double> HH::operator* (const vector<double>& vec)
{
  assert( (*this).size() == vec.size() );
  int n = vec.size();
  vector<double> result = vec;
  double sum = 0;
  for(int i = 0; i < n; ++i){
    sum += m_vec[i]*vec[i];
  }
  for(int i = 0; i < n; ++i){
    result[i] -= 2*vec[i]*sum/(ProdL2(m_vec,m_vec));
  }
  return result;
}


