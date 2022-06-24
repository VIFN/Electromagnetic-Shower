
#include "Vect.h"
#include <iostream>
using namespace std;

Vec :: Vec(int i) : N(i) {entries = new double[N];}
Vec :: Vec(int i, double a) : N(i) {entries = new double[N];
  for(int i = 0; i < N; ++i)
    entries[i] = a;
}

Vec::Vec(Vec &&a) : N(a.N) {
  entries = new double[N];
  for (int i = 0; i < N; ++i)
    entries[i] = a.entries[i];
  a.entries = NULL;
  a.N = 0;
   #ifdef DEBUG
      cout << __PRETTY_FUNCTION__ << "move" <<endl;
      #endif
  }
Vec :: Vec(int i,const double* a) : N(i) {
  entries = new double[N];
  for (int j = 0; j < N; ++j)
    entries[j] = a[j];
}
Vec :: Vec(const Vec& a) {
  N = a.N;
  entries = new double[N];
   #ifdef DEBUG
      cout << __PRETTY_FUNCTION__  <<endl;
      #endif
  for (int i = 0; i < N; ++i)
    entries[i] = a.entries[i];
}
Vec :: ~Vec()
{
     #ifdef DEBUG
      cout << __PRETTY_FUNCTION__  <<endl;
  cout << "Destructor " << this << endl;
   #endif
  delete[] entries;
  #ifdef DEBUG
  cout << "Destroyed" << endl;
#endif
}
void Vec :: SetEntries(int i, double* a)
{
  delete [] entries;
  entries = new double[i];
  N = i;
  for (int j = 0; j < N;++j)
    entries[j] = a[j];
}
int Vec :: size() const {
  #ifdef DEBUG
      cout << __PRETTY_FUNCTION__  <<endl;
      #endif
	return N;}
	
double Vec :: dot(const Vec& A) {
	 cout << __PRETTY_FUNCTION__ << endl;
  if (A.N <= N)
    {
      double sum = 0;
      for (int i = 0; i < A.N;++i)
	sum += A.entries[i]*entries[i];
      return sum;
    }
  else
    {
      double sum = 0;
      for (int i = 0;i < N; ++i)
	sum += A.entries[i]*entries[i];
      return sum;
    }
    
}
void Vec :: swap(int c, int x) {
  if(c < N && x < N)
    {
      double dum;
      dum = entries[c];
      entries[c] = entries[x];
      entries[x] = dum;
    }

}
void Vec :: Print() const {
  cout << "Vector " << this << " : (";
  for (int j = 0; j < N; ++j)
    cout << " " << entries[j] << ",";
  cout << ")" << endl;
} 

double& Vec :: operator[] (int i){return entries[i];}
double Vec :: operator[] (int i) const {return entries[i];}
void Vec :: operator=(const Vec& A){
  if (A.N == N)
    {
      for (int i = 0; i < N; ++i)
	entries[i] = A.entries[i];
    }
  else
    {
      delete[] entries;
      entries = new double[A.N];
      for(int i = 0 ; i < A.N; ++i)
	entries[i] = A.entries[i];
      N = A.N;
    }
}
const Vec& Vec :: operator+=(const Vec& A) {
  *this = *this + A;
  #ifdef DEBUG
  cout << __PRETTY_FUNCTION__  <<endl;
    #endif
  return *this;
 }

Vec Vec :: operator+(const Vec& A) const
{
  if (N == A.N)
    {
      Vec C(N);
      cout << "soma" << &C << endl;
      for (int i = 0; i < N; ++i)
	C[i] = entries[i] + A.entries[i];
      return C;
    }
  if (N < A.N)
    {
      Vec C(A.N);
      cout << "soma1 " << &C <<endl;
      for (int i = 0; i < N; ++i)
	C[i] = entries[i] + A.entries[i];
      for (int j = N; j < A.N; ++j)
	C[j] = A.entries[j];
      return C;
    }
  if (N > A.N)
    {
      Vec C(N);
      cout << "soma2 " << &C <<endl;
      for (int i = 0; i < A.N; ++i)
	C[i] = entries[i] + A.entries[i];
      for (int j = A.N; j < N; ++j)
	C[j] = entries[j];
      return C;
    }
}

Vec Vec :: operator-=(const Vec& A) {
  #ifdef DEBUG
  cout << "-=" << endl;
  #endif
  *this = *this - A;
  #ifdef DEBUG 
  cout << " Outro " << &A << endl;
  #endif
  return *this;
 }


Vec Vec :: operator-(const Vec& A) const
{
  if (N == A.N)
    {
      Vec C(N);
      for (int i = 0; i < N; ++i)
	C[i] = entries[i] - A.entries[i];
      return C;
    }
  if (N < A.N)
    {
      Vec C(A.N);
      for (int i = 0; i < N; ++i)
	C[i] = entries[i] - A.entries[i];
      for (int j = N; j < A.N; ++j)
	C[j] = A.entries[j];
      return C;
    }
  if (N > A.N)
    {
      Vec C(N);
      for (int i = 0; i < A.N; ++i)
	C[i] = entries[i] - A.entries[i];
      for (int j = A.N; j < N; ++j)
	C[j] = entries[j];
      return C;
    }
}

Vec Vec::operator*(const Vec& A) const
{
  if (N == A.N)
    {
      Vec C(N);
      cout << "mult " << &C << endl;
      for (int i = 0; i < N; ++i)
	C[i] = entries[i] * A.entries[i];
      return C;
    }
  else
    {
	     cout << "Não têm as mesmas dimensões, não se alterou o vetor armazenado" << endl;
	     Vec C(N);
	     for (int i = 0; i < N; ++i)
		 C[i] = entries[i];
	     return C;
	}

}

Vec Vec::operator*(double a) const
{
      Vec C(N);
      cout << "mult por constante " << a << " para " << &C << endl;
      for (int i = 0; i < N; ++i)
	C[i] = a*entries[i];
      return C;
}


Vec Vec::operator-()
{
      Vec C(N);
      cout << "simetrico de vetor " << &C << endl;
      for (int i = 0; i < N; ++i)
	C[i] = (-1)*entries[i];
      return C;
}

Vec operator*(double a, const Vec& A)
{
 /*multiplicar a* vec em vez de fazer vec*a (feito na friend fuction)*/
  #ifdef DEBUG
  cout << __PRETTY_FUNCTION__  <<endl;
    #endif


Vec C(A.N);


	for(int i=0; i<A.N;++i)
	{
	C[i]=A[i]*a;
	}
return C;
}

