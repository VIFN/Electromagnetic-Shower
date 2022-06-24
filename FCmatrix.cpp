#include "FCmatrix.h"
#include "Vect.h"
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
using namespace std;

FCmatrix::FCmatrix() { classname="FCmatrix";}

FCmatrix::FCmatrix(double** fM, int fm , int fn) {
  classname="FCmatrix";
      #ifdef DEBUG
      cout << __PRETTY_FUNCTION__  <<endl;
      #endif
  for (int i=0;i<fm;i++)
    {
      Vec aux(fn,fM[i]);
      M.push_back(aux);
    }
   #ifdef DEBUG
      cout << __PRETTY_FUNCTION__  <<endl;
      #endif
}


FCmatrix::FCmatrix(double* fM, int fm , int fn) {
  #ifdef DEBUG
cout << __PRETTY_FUNCTION__  <<endl;
#endif
    classname="FCmatrix";  
 double** M1= new double*[fm];
  for (int i=0;i<fm;i++){
    M1[i]=new double [fn];
    for (int j=0;j<fn;j++){
      M1[i][j]=fM[i*fn+j];
    }
  }
for (int i=0;i<fm;i++){
      Vec a(fn,M1[i]);
      M.push_back(a);
 }

  for (int i=0;i<fm;i++)
    delete []M1[i];
  delete []M1;
   #ifdef DEBUG
      cout << __PRETTY_FUNCTION__  <<endl;
      #endif
}


FCmatrix::FCmatrix(vector <Vec> copy_from){
  for(int i=0;i<copy_from.size();++i)
    M.push_back(copy_from[i]);
    #ifdef DEBUG
      cout << __PRETTY_FUNCTION__  <<endl;
      #endif
}







