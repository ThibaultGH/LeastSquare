#include <iostream>
#include <vector>
#include <string>
#include <math.h>



#include "Mat.hpp"
#include "HouseHolder.hpp"
#include "functions.hpp"



int main(){

  int l = 5;

  vector<double> vec(l,1);

  HH H(vec);

  HH G(H);

  cout << H(0,0) << endl;

  

  

  Mat A;

  char filename[] = "A1.txt";

  A.Load(filename);

  //cout << A;

  int n = A.Row();
  int m = A.Col();

  vector<double> xN(m);
  vector<double> xGS(m);
  vector<double> xH(m);

  vector<double> b(n,1);

  
  NormalSolve(A,xN,b);
  
  for(int i = 0; i < xN.size(); ++i){cout << xN[i] << endl;}

  cout << endl << endl;

  GramSchmidtSolve(A,xGS,b);
  
  for(int i = 0; i < xGS.size(); ++i){cout << xGS[i] << endl;}

  cout << endl << endl;

  HouseholderSolve(A,xH,b);

  for(int i = 0; i < xH.size(); ++i){cout << xH[i] << endl;}

}
