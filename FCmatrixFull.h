//author Botão e Vânia
//Fis computacional 2016-17
#ifndef __FCmatrixFull__
#define __FCmatrixFull__

#include "Vect.h"
#include "FCmatrix.h"
#include <vector>
#include <string>

class FCmatrixFull : public FCmatrix {
public:
// constructors
FCmatrixFull(double** fM, int fm, int fn); //matrix fm x fn

FCmatrixFull(int fm=0,int fn=0);

FCmatrixFull(double* fM, int fm, int fn);

FCmatrixFull(vector<Vec>);

~FCmatrixFull();

// copy constructor
FCmatrixFull(const FCmatrixFull&);


// operators
Vec& operator[] (int) ;
Vec operator[] (int) const;

FCmatrixFull operator=(const FCmatrix&);// copy constructor

FCmatrixFull operator+(const FCmatrix&); // add 2 matrices of any kind

FCmatrixFull operator-(const FCmatrix&); // sub 2 matrices of any kind

FCmatrixFull operator*(const FCmatrix&); // mul 2 matrices of any kind

FCmatrixFull operator*(double lambda); // mul matrix of any kind by scalar

Vec operator*(const Vec&); // mul matrix by Vec

friend FCmatrixFull operator* (double, const FCmatrixFull&);

// virtual inherited
Vec GetRow(int i); // retrieve row i
Vec GetRow(int i) const; // retrieve row i from a constant matrix class
Vec GetCol(int i); // retrieve column i
Vec GetCol(int i) const; // retrieve column i from a constant matrix class
double Determinant();
void Print();
void swapRows(int,int);

//metodos herdados mas para a banded e sparse
double* GetRaw(int i=0);
int Getn() const {return M.size();};

// row max element index (scaled by s, from i on)
int GetRowMax(int i=0, int s=0);
// col max element index  (scaled by s, from j on)
int GetColMax(int j=0, int s=0);


private:
int *rowindices; // row indices (0,1,2,...)
int *colindices; // column indices (0,1,2,..)
};
#endif
