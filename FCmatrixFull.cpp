//author Boto e Vânia


#include "FCmatrix.h"
#include "FCmatrixFull.h"
#include "Vect.h"
#include <cmath>
#include <stdlib.h>
#include <vector>
#include <iostream>
using namespace std;

FCmatrixFull::FCmatrixFull(double** fM, int fm, int fn): FCmatrix(fM,fm,fn){

  classname="FCmatrixFull";
   #ifdef DEBUG
      cout << __PRETTY_FUNCTION__  <<endl;
      #endif
  rowindices= new int[M.size()];
  colindices= new int[M[0].size()];
}

FCmatrixFull::FCmatrixFull(int fm,int fn):FCmatrix(){
  classname="FCmatrixFull";
  Vec V(fn,0.);
  for (int i = 0; i < fm; ++i)
  {
    M.push_back(V);
  }
  rowindices= new int[M.size()];
  colindices= new int[M[0].size()];
}

FCmatrixFull::FCmatrixFull(double* fM,int fm,int fn): FCmatrix(fM,fm,fn) {
  #ifdef DEBUG
cout << __PRETTY_FUNCTION__  <<endl;
#endif
  classname="FCmatrixFull";
  rowindices= new int[M.size()];
  colindices= new int[M[0].size()];
    #ifdef DEBUG
      cout << __PRETTY_FUNCTION__  <<endl;
      #endif
}

FCmatrixFull::FCmatrixFull(vector <Vec> copy_from): FCmatrix(copy_from) {
  classname="FCmatrixFull";
   #ifdef DEBUG
      cout << __PRETTY_FUNCTION__  <<endl;
      #endif
  rowindices=new int[M.size()];
  colindices=new int[M[0].size()];
  }

FCmatrixFull::~FCmatrixFull() {
     #ifdef DEBUG
      cout << __PRETTY_FUNCTION__  <<endl;
      #endif
  delete rowindices;
  delete colindices;
}


FCmatrixFull::FCmatrixFull(const FCmatrixFull& copy_from)
{
  classname="FCmatrixFull";
  rowindices=new int[copy_from.M.size()];
  colindices=new int[copy_from[0].size()];
   #ifdef DEBUG
      cout << __PRETTY_FUNCTION__  <<endl;
      #endif
  M.clear();
  Vec auxiliar(copy_from.ncols());
  for(int j=0;j<copy_from.nrows();++j)
    {
      for(int i=0;i<copy_from.ncols();++i)
	{
	  auxiliar[i]=copy_from[j][i];
	}
      M.push_back(auxiliar);
    }
  
}


//DETERMINANTE FUNCIONAL MAS SEM TESTES DE PIVOT (PODE DIVIDIR POR 0 ASSIM)

double FCmatrixFull::Determinant(){
    #ifdef DEBUG
      cout << __PRETTY_FUNCTION__  <<endl;
      #endif
        double determinante=0;
//Método de laplace
 if(M.size()<3)
    {     
       cout << "laplace" <<endl;
       if (M[0].size()!=M.size())
          {
            cout <<"Inválido! Matriz não quadrada"<<endl;
            return 0;
          }
        
        if(M.size()==2)
          return (M[0][0]*M[1][1]-M[0][1]*M[1][0]);
        else
          for(int i=0;i<M.size();++i)
            {
              vector<Vec> vv;
              for(int j=0;j<M[0].size();j++)
              	{
                  if(j!=i)
                  {
              	    double *aux=new double[M.size()-1];
              	    for(int k=1;k<M.size();k++)
              	      aux[k-1]=M[j][k];
              	    Vec aux2 (M.size()-1,aux);
              	    delete []aux;
              	    vv.push_back(aux2);
                  }
              	   
              	}
              FCmatrixFull aux3(vv);
              vv.clear();
              determinante+=pow(-1,i)*M[i][0]*aux3.Determinant();
            }
    }

//Método de Gauss
  else
    {
      cout << "determinante por gauss" << endl;
      determinante=1.0;
      vector <Vec> A = M;
      FCmatrixFull B(A);
      Vec auxiliar(A.size(),0.0);
      for(int k=0 ; k<(B.nrows()-1);++k)
        {
                int pivot= B.GetColMax(k,k);
                if(pivot!=k)
                {
                  B.swapRows(k,pivot);
                  determinante*=-1;
                  cout << "troca de linhas " << k << " e " << pivot << endl;
                }
                if(B[k][k]<(0.0001)) 
                {
                  cout << "matriz singular"<<endl;
                  return 0;
                }  
              
          for (int i = k+1; i < B.nrows() ; ++i)
          {
            double lambda=B[i][k]/B[k][k];
            if(lambda!=0)
              B[i]-=lambda*B[k];
          }
        }

      for (int i = 0; i < B.nrows(); ++i)
      {
        determinante*=B[i][i];
      }
    }

  return determinante;
}



Vec& FCmatrixFull::operator[](int i){
  return M[i];
   #ifdef DEBUG
      cout << __PRETTY_FUNCTION__  <<endl;
      #endif
}


Vec FCmatrixFull::operator[](int i) const{
  
  return M[i];
   #ifdef DEBUG
      cout << __PRETTY_FUNCTION__  <<endl;
      #endif
}


FCmatrixFull FCmatrixFull::operator=(const FCmatrix& build)
{
    #ifdef DEBUG
      cout << __PRETTY_FUNCTION__  <<endl;
      #endif
  M.clear();
  for(int i=0; i<build.ncols(); ++i)
  {
    M.push_back(build.GetRow(i));
  }
    #ifdef DEBUG
      cout << __PRETTY_FUNCTION__  <<endl;
      #endif
}

FCmatrixFull FCmatrixFull::operator+(const FCmatrix& A) {
  if(M.size()==A.nrows() && M[0].size()==A.ncols())
    {
      for(int i=0; i<M.size(); ++i)
      	for(int j=0;j<M[0].size();++j)
    	    M[i][j]+=A[i][j];
      return M;
    }
else
  cout<<"As matrizes não têm o mesmo tamanho. Nada foi feito à matriz armazenada" <<endl;
  return M;
}

FCmatrixFull FCmatrixFull::operator-(const FCmatrix& A) {
  if(M.size()==A.nrows() && M[0].size()==A.ncols())
    {
      for(int i=0; i<M.size(); ++i)
	for(int j=0;j<M[0].size();++j)
	  M[i][j]-=A[i][j];
      return M;		    
    }
else
  cout<<"As matrizes não têm o mesmo tamanho. Nada foi feito à matriz armazenada" <<endl;
  return M;
}


FCmatrixFull FCmatrixFull::operator*(double lambda) {
  for(int i=0; i<M.size(); ++i)
    for(int j=0;j<M[0].size();++j)
      M[i][j]=M[i][j]*lambda;
  return M;
}

FCmatrixFull operator*(double lambda, const FCmatrixFull& A) {
    #ifdef DEBUG
  cout << __PRETTY_FUNCTION__  <<endl;
    #endif
  FCmatrixFull K(A);
   for(int i=0; i<K.M.size(); ++i)
    for(int j=0;j<K[0].size();++j)
      K.M[i][j]=K.M[i][j]*lambda;
  return K;
}



Vec FCmatrixFull::operator*(const Vec& v) {
  Vec b(M.size(),0.0);
  for(int i=0; i<M.size(); ++i)
    for(int j=0;j<M[0].size();++j)
    {
      M[i][j]=M[i][j]*v[j];
      b[i]+=M[i][j];
    }
  return b;
}


FCmatrixFull FCmatrixFull::operator*(const FCmatrix& copyfrom) {
 

  #ifdef DEBUG
  cout << __PRETTY_FUNCTION__  <<endl;
    #endif
  if (M[0].size()!=copyfrom.nrows())
    {
      cout <<"Dimensões inválidas para multiplicação!"<<endl;
      return *this;
    }

  vector<Vec> vv;
  for (int i=0;i<M.size();i++) 
  {
    Vec ve(copyfrom.ncols(),0.0);

        for (int j=0;j<copyfrom.ncols();j++)
        {

          for (int v = 0; v < M[0].size(); ++v)
          {
                ve[j]+=M[i][v]*copyfrom[v][j];
          }
        }
      
    vv.push_back(ve);
  }

  FCmatrixFull res(vv);
  vv.clear();
  return res;
}


Vec FCmatrixFull::GetRow(int i){
  return M[i];
}

Vec FCmatrixFull::GetRow(int i) const{
  return M[i];
}

Vec FCmatrixFull::GetCol(int i){
  Vec auxiliar(M.size());
  for(int j=0;j<M.size();++j)
    {
      auxiliar[j]=M[j][i];
    }
  return auxiliar;
}

Vec FCmatrixFull::GetCol(int i) const{
  Vec auxiliar(M.size());
  for(int j=0;j<M.size();++j)
    {
      auxiliar[j]=M[j][i];
    }
  return auxiliar;
}

void FCmatrixFull::swapRows( int i, int j)
{
  Vec auxiliar(M[0].size());
  auxiliar=M[i];
  M[i]=M[j];
  M[j]=auxiliar;
}

void FCmatrixFull::Print() {
  cout << "matrix printing " << this << endl;
    #ifdef DEBUG
      cout << __PRETTY_FUNCTION__  <<endl;
      #endif
      cout << "_________" << endl;
  for(int i=0;i<(M.size());++i)
    {
      cout << "|vector"<< i << "=";
      M[i].Print();
    }
       cout << "_________" << endl;
}



int FCmatrixFull::GetRowMax (int i, int s)
{
  double teste=0;
  int max=0;
  if(i<=M.size())
    {
      for(int j=s;j<M[i].size();++j)
      {
        if(fabs(teste)<fabs(M[i][j]))
        {
          teste=M[i][j];
          max=j;
        }
      }
    }
  return max;
    #ifdef DEBUG
      cout << __PRETTY_FUNCTION__  <<endl;
      #endif
}

int FCmatrixFull::GetColMax (int i, int s)
{
  double teste=0;
  int max=0;
  if(i<=M[0].size())
    {
      for(int j=s;j<M.size();++j)
      {
        if(fabs(teste)<fabs(M[j][i]))
        {
          teste=M[j][i];
          max=j;
        }
      }
    }
  return max;
    #ifdef DEBUG
      cout << __PRETTY_FUNCTION__  <<endl;
      #endif
}


//useless one that isnt ever called

double* FCmatrixFull::GetRaw(int i){
  double *a=new double[1];
  a[0]=9999;
  return a;
}
