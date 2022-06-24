#include "Particle.h"
#include "Vect.h"
#include <cmath>
#include "TMath.h"
#include <iostream>
using namespace std;

int main()
{
	Vec A(3.,1.);
	//A.Print();
	Particle X(0.,1000,A,0.);
	X.SetTheta(TMath::Pi()/4.);
	X.SetPhi(TMath::Pi()/4.);
	Vec B(3.,10.);
	X.SetDeath(B, 100.,1000.);
 	 X.Print();

return 0;
}