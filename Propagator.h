#ifndef __Propagator__
#define __Propagator__
#include "Particle.h"
#include "Formula.h"
#include <vector>
using namespace std;

class Propagator {
public:
Propagator(int N, double Energy); //gerar N fotões com energia igual e vetor momento segundo z e analisar a sua cascata
~Propagator(){delete F;}
vector <Particle> GetResults(int& i){
	i=P.size();
	return P;}

private:
vector<Particle> P;
Formula *F;
void PhotonPropagator(int i); //propagar fotao i
double Generator_r();//gerar distancia ate fotao se dividir de acordo com distribuicao PDF


void EletronPropagator(int i); //propagar eletrao i
void PositronPropagator(int i); //propagar positrao i


};
#endif