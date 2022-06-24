//author Botão e Vânia
//Fis computacional 2016-17
#ifndef __Vec__
#define __Vec__
class Vec {
 public:
  Vec(int i=1); //default constructor
  Vec(int, double); //set N elements equal to value
  Vec(Vec&&) ; // move constructor
  Vec(int, const double* ); //set N elements from array
  Vec(const Vec&); // copy constructor
  ~Vec(); //destructor
  void SetEntries (int, double* );
  int size() const; //Vec size

  //ODEPOINT METHODS (same as others in vec but renamed)
  double T() const {return entries[0];};
  double X(int i) const {return entries[i+1];};
  double* GetArray() { return entries;}


  double dot(const Vec&); //produto interno
  void swap(int, int); //swap Vec elements
  void Print() const; //class dump
  double& operator[] (int);
  double operator[] (int) const; //Vec is declared as const
  void operator=(const Vec&);
  const Vec& operator+=(const Vec&);
  Vec operator+(const Vec&) const;
  Vec operator-=(const Vec&);
  Vec operator-(const Vec& ) const;
  Vec operator*(const Vec& ) const; //x1x2,y1y2,z1z2
  Vec operator*(double ) const; //Vec.operator*(k) = Vec*scalar
  Vec operator-();
  friend Vec operator* (double, const Vec&);
 private:
  int N=0;
  double* entries;
};
#endif
