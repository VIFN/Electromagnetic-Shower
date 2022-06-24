#ifndef __FCmatrix__
#define __FCmatrix__

//author Botão
//Fis computacional 2016-17
//history
//15/11/2016: acabou FCmatrixFull
//24/11/2016: acabou FCmatrixBanded


#include "Vect.h"
#include <vector>
#include <string>
using namespace std;
class FCmatrix{
 public:
  FCmatrix();
  FCmatrix(double** fM, int fm, int fn);
  FCmatrix(double* fM, int fm, int fn);
  FCmatrix(vector<Vec>);
  ~FCmatrix(){;};

// operators  
virtual Vec& operator[] (int) =0;
virtual Vec operator[] (int) const=0;
// methods
 virtual Vec GetRow(int i)=0 ; // retrieve row i
  virtual Vec GetRow(int i) const=0 ; // retrieve row i
 virtual  Vec GetCol(int i)=0 ; // retrieve column i
  virtual  Vec GetCol(int i) const=0 ; // retrieve column i
 virtual double Determinant()=0 ;
virtual void Print()=0;

// row max element index (scaled by s, from i on)
virtual int GetRowMax(int i=0, int s=0)=0 ;
// col max element index  (scaled by s, from j on)
virtual int GetColMax(int j=0, int s=0)=0;
virtual void swapRows(int,int)=0;

//metodos da banded para o eqsolver
virtual double* GetRaw(int)=0; 
virtual int Getn() const=0;

string getname() const {return classname;}
string& getname() {return classname;}
int nrows() const{return M.size();}
  int ncols() const {return M[0].size();}

protected:
vector<Vec> M;
string classname;
};
#endif
