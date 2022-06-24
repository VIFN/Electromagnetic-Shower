#include "Particle.h"
#include "TMath.h"
#include <iostream>
using namespace std;

Particle::Particle(int Type, double Energy, Vec Positions, double Time): type(Type), positions(Positions), time(Time) , positionsf(Positions), timef(Time), energyloss(0)//energia recebida em MeV
{
	#ifdef DEBUG
	cout << __PRETTY_FUNCTION__  <<endl;
	#endif	

	double mass=0.;
	double momentaux=0.;
	switch(type)
	{
		case 0: //tipo corresponde a fotão
				mass=0.;
				moment.SetZ(Energy);
				momentf.SetZ(Energy);
				#ifdef DEBUG
				cout << __PRETTY_FUNCTION__  <<endl;
				moment.Print();
				#endif	
				break;

		case 1: //corresponde a eletrão
				mass= 0.5109989461; //massa em MeV c^-2
				momentaux=Energy*Energy-mass*mass;
				moment.SetZ(sqrt(momentaux));
				momentf.SetZ(sqrt(momentaux));
				#ifdef DEBUG
				cout << __PRETTY_FUNCTION__  <<endl;
				moment.Print();
				#endif	
				break;

		case 2: //corresponde a positrao
				mass= 0.5109989461; //massa em MeV c^-2
				momentaux=Energy*Energy-mass*mass;
				moment.SetZ(sqrt(momentaux));
				momentf.SetZ(sqrt(momentaux));
				#ifdef DEBUG
				cout << __PRETTY_FUNCTION__  <<endl;
				moment.Print();
				#endif	
				break;

		default:
				cout << __PRETTY_FUNCTION__ << "Non existing type" << endl;
				exit(1);
				break;
	}

	#ifdef DEBUG
	cout << __PRETTY_FUNCTION__  <<endl;
	#endif	
}

TVector3 Particle::GetMoment()
{
	#ifdef DEBUG
	cout << __PRETTY_FUNCTION__  <<endl;
	#endif		
	return moment;
}

TVector3 Particle::GetMomentf()
{
	#ifdef DEBUG
	cout << __PRETTY_FUNCTION__  <<endl;
	#endif		
	return momentf;
}


int Particle::GetType()
{
	#ifdef DEBUG
	cout << __PRETTY_FUNCTION__  <<endl;
	#endif		
	return type;
}

double Particle::GetTime()
{
	#ifdef DEBUG
	cout << __PRETTY_FUNCTION__  <<endl;
	#endif		
	return timef;
}


Vec Particle::GetPositions()
{
	#ifdef DEBUG
	cout << __PRETTY_FUNCTION__  <<endl;
	#endif	


	return positions;
}

Vec Particle::GetFPositions()
{
	#ifdef DEBUG
	cout << __PRETTY_FUNCTION__  <<endl;
	#endif	


	return positionsf;	
}

void Particle::SetXYZ(double x, double y, double z)
{
	#ifdef DEBUG
	cout << __PRETTY_FUNCTION__  <<endl;
	#endif	
		moment.SetXYZ(x,y,z);
	momentf.SetXYZ(x,y,z);
}


void Particle::SetTheta(double theta)
{
	#ifdef DEBUG
	cout << __PRETTY_FUNCTION__  <<endl;
	#endif	
	momentf.SetTheta(theta);
}

void Particle::SetPhi(double phi)
{
	#ifdef DEBUG
	cout << __PRETTY_FUNCTION__  <<endl;
	#endif	
	momentf.SetPhi(phi);
}

void Particle::SetMag(double mag)
{
	#ifdef DEBUG
	cout << __PRETTY_FUNCTION__  <<endl;
	#endif	
	momentf.SetMag(mag);
}

void Particle::SetDeath(Vec Positions, double Timef,double Energyf)
{
	#ifdef DEBUG
	cout << __PRETTY_FUNCTION__  <<endl;
	#endif //se quiser atualizar entre cada iteração, alterar aqui. Caso queiramos fazer perfil temporal (não pedido)
	/*time=timef;
	moment=momentf;
	positions=positionsf;*/

    timef=Timef;
	Particle::SetEnergy(Energyf);
	positionsf=Positions;

}

void Particle::SetEnergyLoss(double Energyloss)
{
	energyloss=Energyloss;
}

void Particle::SetEnergy(double Energy)
{
	#ifdef DEBUG
	cout << __PRETTY_FUNCTION__  <<endl;
	#endif
	switch(type)
	{
		case 0:
			momentf.SetMag(Energy);
			break;
		case 1:
			momentf.SetMag(sqrt(Energy*Energy- 0.5109989461* 0.5109989461));
			break;
		case 2:
			momentf.SetMag(sqrt(Energy*Energy- 0.5109989461* 0.5109989461));
			break;		
		default:
			cout << __PRETTY_FUNCTION__ << "Non existing type" << endl;
			exit(1);
			break;	
	}
}

void Particle::SetInitialEnergy(double Energy)
{
	#ifdef DEBUG
	cout << __PRETTY_FUNCTION__  <<endl;
	#endif
	switch(type)
	{
		case 0:
			moment.SetMag(Energy);
			break;
		case 1:
			moment.SetMag(sqrt(Energy*Energy- 0.5109989461* 0.5109989461));
			break;
		case 2:
			moment.SetMag(sqrt(Energy*Energy- 0.5109989461* 0.5109989461));
			break;		
		default:
			cout << __PRETTY_FUNCTION__ << "Non existing type" << endl;
			exit(1);
			break;	
	}
}

double Particle::GetMass()
{
	double mass=0.;
	switch(type)
	{
		case 0: //tipo corresponde a fotão
				return mass;
				break;

		case 1: //corresponde a eletrão
				mass= 0.5109989461; //massa em MeV c^-2
				return mass;
				break;

		case 2:
				mass= 0.5109989461; //massa em MeV c^-2
				return mass;
				break;

		default:
				cout << __PRETTY_FUNCTION__ << "Non existing type" << endl;
				exit(1);
				break;
	}
}

double Particle::GetEnergy()
{
	double mass=0.;
	double Energy=0.;
	switch(type)
	{
		case 0: //tipo corresponde a fotão
				Energy= sqrt(momentf.X()*momentf.X()+momentf.Y()*momentf.Y()+momentf.Z()*momentf.Z());
				return Energy;
				break;

		case 1: //corresponde a eletrão
				mass= 0.5109989461; //massa em MeV c^-2
				Energy= sqrt(momentf.X()*momentf.X()+momentf.Y()*momentf.Y()+momentf.Z()*momentf.Z())*sqrt(momentf.X()*momentf.X()+momentf.Y()*momentf.Y()+momentf.Z()*momentf.Z())+mass*mass;
				return sqrt(Energy);
				break;

		case 2:
				mass= 0.5109989461; //massa em MeV c^-2
				Energy=  sqrt(momentf.X()*momentf.X()+momentf.Y()*momentf.Y()+momentf.Z()*momentf.Z())* sqrt(momentf.X()*momentf.X()+momentf.Y()*momentf.Y()+momentf.Z()*momentf.Z())+mass*mass;
				return sqrt(Energy);
				break;

		default:
				cout << __PRETTY_FUNCTION__ << "Non existing type" << endl;
				exit(1);
				break;
	}
}

double Particle::GetInitialEnergy()
{
	double mass=0.;
	double Energy=0.;
	switch(type)
	{
		case 0: //tipo corresponde a fotão
				Energy= sqrt(moment.X()*moment.X()+moment.Y()*moment.Y()+moment.Z()*moment.Z());
				return Energy;
				break;

		case 1: //corresponde a eletrão
				mass= 0.5109989461; //massa em MeV c^-2
				Energy= sqrt(moment.X()*moment.X()+moment.Y()*moment.Y()+moment.Z()*moment.Z())*sqrt(moment.X()*moment.X()+moment.Y()*moment.Y()+moment.Z()*moment.Z())+mass*mass;
				return sqrt(Energy);
				break;

		case 2:
				mass= 0.5109989461; //massa em MeV c^-2
				Energy=  sqrt(moment.X()*moment.X()+moment.Y()*moment.Y()+moment.Z()*moment.Z())* sqrt(moment.X()*moment.X()+moment.Y()*moment.Y()+moment.Z()*moment.Z())+mass*mass;
				return sqrt(Energy);
				break;

		default:
				cout << __PRETTY_FUNCTION__ << "Non existing type" << endl;
				exit(1);
				break;
	}
}

double Particle::GetMag()
{
	#ifdef DEBUG
	cout << __PRETTY_FUNCTION__  <<endl;
	#endif	
	return momentf.Mag();
}

double Particle::GetTheta()
{
	#ifdef DEBUG
	cout << __PRETTY_FUNCTION__  <<endl;
	#endif	
	if(momentf.X()<0)
		return atan(momentf.Y()/momentf.X())+TMath::Pi();
	if (momentf.X()>0)
		return atan(momentf.Y()/momentf.X());
	if(momentf.X()==0)
		{
			if (momentf.Y()<0)
			{
				return TMath::Pi()/2.;
			}
			if (momentf.Y()>0)
			{
				return TMath::Pi()/2.;
			}
			if (momentf.Y()==0)
			{
				return 0;
			}
		}
}

double Particle::GetPhi()
{
	#ifdef DEBUG
	cout << __PRETTY_FUNCTION__  <<endl;
	#endif	

	return atan(sqrt(momentf.X()*momentf.X()+momentf.Y()*momentf.Y())/momentf.Z());	
}

double Particle::GetEnergyLoss()
{
	return energyloss;
}

double Particle::GetInitTime()
{
	return time;
}

void Particle::Rotate(TRotation m)
{
	#ifdef DEBUG
	cout << __PRETTY_FUNCTION__  <<endl;
	#endif	
	momentf.Transform(m);
}

void Particle::Print(int modo)
{
	if(modo==0)
	{
	cout << "A partícula é do tipo " << type << " (0= fotao, 1= eletrao, 2=positrao) " << endl;
	cout << "Tem massa= " << Particle::GetMass() << " MeV c^{-2} e Energia= "<< Particle::GetEnergy() << " MeV " << endl;
	cout << "Esta partícula formou-se no instante de tempo " << time << " nano s e na posição ("<< positions[0] << ", " << positions[1] << ", " << positions[2] << ") cm" <<endl;
	cout  << "o Vetor momento é da forma : ";
		moment.Print();
	}

	momentf.Print();

	if (timef!=time && modo!=2)
	{
		cout  <<"Esta partícula desapareceu em " ;
		printf("%.10f nano s",timef);
		cout  << " na posição ("<< positionsf[0] << ", " << positionsf[1] << ", " << positionsf[2] << ") cm" <<endl ;
		cout << endl;
	}
}