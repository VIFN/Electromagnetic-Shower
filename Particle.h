#ifndef __Particle__
#define __Particle__
#include "Vect.h"
#include "TVector3.h"
#include "TRotation.h"
class Particle {
public:
Particle(int type,double Energy, Vec Positions, double Time); //energia convertida em momento segundo Z

//getting a copy of the raw datamembers
TVector3 GetMoment();
TVector3 GetMomentf();
int GetType();
Vec GetPositions(); //returns the initial positions
Vec GetFPositions(); //returns the final positions


//alterar o vetor momento
void SetXYZ(double x,double y,double z);
void SetTheta(double);
void SetPhi(double);
void SetMag(double); // alterar módulo do momento linear diretamente
void SetDeath(Vec Positions, double Timef, double Energyf);
void SetEnergy(double Energy);
void SetEnergyLoss(double Energy);
void SetInitialEnergy(double Energy);

//receber caraterísticas rapidamente
double GetMass();
double GetEnergy();
double GetInitialEnergy();
double GetMag(); //receber módulo do momento linear diretamente
double GetTheta();
double GetPhi();
double GetEnergyLoss();
double GetInitTime();
double GetTime();

//rodar vector com uma matriz de rotação (não foi usado)
void Rotate(TRotation);

//orientados para debugging
void Print(int modo=0);


private:
int type; //type 0=fotão, 1=eletrão, 2=positrão, 
TVector3 moment; // momento linear da partícula nos eixos x,y,z inicial
//Energia pode ser obtida diretamente com E2= p2 c2 + mo2 c4 (usamos unidades naturais com c=1)
//massa obtida a partir do tipo diretamente (precisa de método)
TVector3 momentf; //momento linear da particula nos eixos x,y,z que é constantemente alterado nos eletroes e positroes

double time; // tempo em que se formou
double timef; //tempo em que desapareceu ou ocorreu interação
double energyloss; //dE/dx nesse intervalo

Vec positions; // Vector com 3 entradas da seguinte form: 
//coordenadas de origem :  x0, y0, z0

Vec positionsf; 
//cooordenadas de desaparecimento : xf , yf, zf
};
#endif