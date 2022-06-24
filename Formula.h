#ifndef __Formula__
#define __Formula__
#include "TFormula.h"
#include "TVector3.h"
#include "Particle.h"
#include "FCmatrixFull.h"
#include "TF1.h"
class Formula {
public:
	Formula(); 
~Formula(){;}
double PhotonEnergyDivider(double energy);
void PhotonAngles(double Theta,double Phi,Particle& Pe);


double BremRadiation(double energy);
double EnergyBrem(double energy);
void BremAngles(double Theta, double Phi, Particle& Pe);

double AniqPositron(double energy);
double AniqEnergy(double energy);
void AniqAngles(double Theta,double Phi,Particle& F1,Particle& F2,double u, int mode);

double EnergyLoss(int type,double energy);


private:
//relacionado com divisão de energia nos fotões
	double Emin;
	double Emax;
	double a;
	double g0;
	double PMAX;
	TFormula b;
	TF1 g1;
	TF1 g2;
	TF1 phi;
	TF1 f;

//relacionado com angulos
	FCmatrixFull mat; //matriz de transformação
	TVector3 pos_aux;
	Vec pos;
	TF1 auxiliarangles;	
	TF1 anglesbrem;
	double BremanglesAux(double Energy);

//relacionado com efeito brem
	double xr; //geracao de aleatorios
	double y;
	double ur;
	double auxA;
	//ate aqui relacionado com aleatorios
	TF1 aux2;
	TF1 auxiliar;
	TF1 dsig;
	double adaptiveSimpsonsAux(double a, double b, double epsilon, double S, double fa, double fb, double fc, int bottom);
	double adaptiveSimpsons(double a, double b,  double epsilon, int maxRecursionDepth); 

//relacionado com aniquilacao de positroes
	double Umax;
	double Umin;
	TF1 Aniq;


//relacionado com perda energia
	double betha;
	double I;
	double Gamma;
	double bethagamma;
	double B;
	double aa;
	double delta;
	double C0;
	double m;
	double U1;
	double U0;
	double C;
	double fm;
	TFormula de;
};
#endif