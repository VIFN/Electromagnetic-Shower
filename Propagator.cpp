#include "Propagator.h"
#include <ctime>
#include <cmath>
#include <vector>
#include "TRandom3.h"
#include "TF1.h"
#include "TFormula.h"
#include "Vect.h"
#include "TMath.h"
#include "cFCgraphics.h"
#include "TAxis.h"

#include <iostream>
using namespace std;



const double c=299792458; //velocidade da luz
const double masseletron=0.5109989461; //massa eletrao/positrao



Propagator::Propagator(int N, double Energy)
{
	#ifdef DEBUG
	cout << __PRETTY_FUNCTION__  <<endl;
	#endif	
	F= new Formula();
	Vec position(3.,0.);
	for (int i = 0; i < N; ++i)
	{
		Particle newp(0,Energy,position,0.);
		P.push_back(newp);
	}

	gRandom= new TRandom3(0.);
	double test=0;
	int i=0;
	for ( i = 0; i < P.size(); ++i)
	{
		switch(P[i].GetType())
		{
			case 0: //é do tipo fotao 
				Propagator::PhotonPropagator(i);
				//cout << "one photon died"<<endl<<endl;
				
				break;

			case 1: //e do tipo eletrao
				Propagator::EletronPropagator(i);
				//cout << "one eletron died"<<endl<<endl;

				break;

			case 2: // e do tipo positrao
				Propagator::PositronPropagator(i);

				//cout << "One positron died  "<<  endl<<endl;

				break;

			default:
				cout << __PRETTY_FUNCTION__ << "Non existing type" << endl;
				exit(1);
				break;

		}
	}


	double time= P[i-1].GetTime();

	cout << i-1 << " particulas foram detetadas em "<< time <<" nano segundos"<<endl;
}




///////// RELACIONADO COM FOTOES //////////////////////////////////

void Propagator::PhotonPropagator(int i)
{
	#ifdef DEBUG
	cout << __PRETTY_FUNCTION__  <<endl;
	#endif	
	double Energy= P[i].GetEnergy();
	//cout << "Propagando fotão com energia total " << Energy << " MeV."<<endl;

	//descobrir a distancia percorrida por este fotão antes de decair
	double Rint= Propagator::Generator_r();
	//cout << "Percorreu " << Rint << "cm até decair.";

	double timef= pow(10,7)/c;
	timef= P[i].GetTime()+  Rint*timef;
	//printf("%.10f nano s\n",timef);

	//descobrir posição final
	Vec V0(P[i].GetPositions()); 
	double Theta= P[i].GetTheta(); 
	double Phi= P[i].GetPhi();
	
	//cout << "antes com raio int "<< Rint<< " theta " << Theta << " phi "<<Phi<<endl;
	//V0.Print();
	//coordinates in the lab referencial
	V0[0]+=Rint*sin(Phi)*cos(Theta);
	V0[1]+=Rint*sin(Phi)*sin(Theta);
	V0[2]+=Rint*cos(Phi);
	//cout <<"after party "<<endl;
	//V0.Print();
	P[i].SetDeath(V0,timef,Energy);
	//P[i].Print();

	//dividir a energia

	double RedEnergy= F->PhotonEnergyDivider(Energy);
	//cout << RedEnergy << endl;
	double EnergyElectron= RedEnergy*Energy; //energia cinetica
	double EnergyPositron= Energy- EnergyElectron; //energia cinetica
	//cout << " energia total eletrao " <<  EnergyElectron+masseletron << " " << " energia total positrao " << EnergyPositron+masseletron << endl;
 	
	Particle Pe(1,EnergyElectron,V0,timef);
	Particle Pp(2,EnergyPositron,V0,timef);
	F->PhotonAngles(Theta,Phi,Pe);
	F->PhotonAngles(Theta,Phi,Pp);


	P.push_back(Pe);
	P.push_back(Pp);

	#ifdef DEBUG
	cout << __PRETTY_FUNCTION__  <<endl;
	#endif	
}



#ifdef DEBUG
int i=0.;
#endif
//divisor energia



double Propagator::Generator_r(){

	#ifdef DEBUG
	cout << __PRETTY_FUNCTION__  <<endl;
	#endif
	 double Densidade=2.699; //g/cm³ densidade do aluminio
	 double x0= 24.01;  //comprimento de radiação associado ao material (g/cm²)
	#ifdef IRON
 	Densidade=7.874;
 	x0=13.84;
 	#endif
 	#ifdef CHUMBO
 	x0=6.37;
 	Densidade=11.35;
 	#endif
 	#ifdef GRAFITE
 	x0=42.70;
 	Densidade=2; //segundo tabela.Leo em anexo (diferente do valor da tabela no site pdg.lbl.gov)
 	#endif

	double y=gRandom->Uniform(0.,1.);
	double ppar=(7*Densidade)/(9*x0);
	double rint;
	//cout << ppar << endl;
	rint=log(y)/(-ppar);

	
	return rint;


}




////////////////////////////////////////////////////////////////// 
                 // ELETRON RELATED //
///////////////////////////////////////////////////////////////


void Propagator::EletronPropagator(int i)
{
	#ifdef DEBUG
	cout << __PRETTY_FUNCTION__  <<endl;
	#endif	
	double Energy= P[i].GetEnergy();
	double rint=0.;
	double probability;
	double energyloss;
	double dedx;
	double y;
	double Finalenergy;

	while((Energy-masseletron)>5e-5)
	{
		//cout << "Propagando eletrão com energia " << Energy << " MeV."<<endl;

		//descobrir a distancia percorrida por este eltrao antes de decair
		
		probability= F->BremRadiation(Energy-masseletron); //integral de 50 eV até à energia da particula

		//cout <<"probabilidade"<< probability<<endl;

		y=gRandom->Uniform(0.,1.); //gerando R interacao de acordo log(y)/-pint
		rint=log(y)/(-probability);
		//cout <<"raios de interacao " << rint <<endl;
		dedx= F->EnergyLoss(1,Energy-masseletron); 
		energyloss=rint*dedx; //energia cinetica perdida
		//cout << "energy lost ="<< energyloss << " MeV and walked " << rint << "cm" << endl;
		P[i].SetEnergyLoss(energyloss);

		Finalenergy= Energy-energyloss;
		//cout<<"Finalenergy="<<Finalenergy<<endl;
		//cout << "energia final "<<Finalenergy<<endl;
		double Theta= P[i].GetTheta(); 
		double Phi= P[i].GetPhi();
			//divisao de energia
		double EnergyFotao=0.;
		if (Finalenergy-masseletron > 5e-5)
		{
			EnergyFotao= F->EnergyBrem(Finalenergy-masseletron);
			//cout << "energy fotao " <<endl;
			//cout << EnergyFotao <<endl;
			Finalenergy-=EnergyFotao;
	
		}


		if(Finalenergy-masseletron < 5e-5)//correcao no caso de ter perdido demasiada energia e deixar de ser detetado
		{
			//cout << "raio "<<rint<<endl;
			Finalenergy=5e-5+masseletron;
			energyloss= Energy-Finalenergy;
			//cout << energyloss << endl;
			rint=energyloss/dedx;
			//cout <<"raios "<<rint<<endl;
		}
		//printf("%.50f c initial velocity\n",initialvelocity);

		//printf("%.50f c final velocity\n",finalvelocity);
		double finalvelocity=sqrt((Finalenergy-masseletron)*(Finalenergy+masseletron))/Finalenergy;
		double initialvelocity=sqrt((Energy-masseletron)*(Energy+masseletron))/Energy;
		//printf("%.50f c initial velocity\n",initialvelocity);
//cout<<"finalvelocity "<<finalvelocity<<endl;
		//printf("%.50f c final velocity\n",finalvelocity);
		double timef=fabs(P[i].GetTime() + rint*1e7*log(initialvelocity/finalvelocity));
		//cout << timef << " nano s percorridos " <<endl;

		//descobrir posição final
		Vec V0(P[i].GetFPositions()); 
		


		//coordinates in the lab referencial
		V0[0]=V0[0]+ rint*sin(Phi)*cos(Theta);
		V0[1]=V0[1]+ rint*sin(Phi)*sin(Theta);
		V0[2]=V0[2]+ rint*cos(Phi);

		if (EnergyFotao>5)
		{
		//	cout << "energia a entrar "<< EnergyFotao << endl;
			Particle Pf(0,EnergyFotao,V0,timef);
			F->BremAngles(Theta,Phi,Pf);
		//	cout << "energia a sair "<< Pf.GetEnergy() << endl;
			P.push_back(Pf);
		}

	P[i].SetDeath(V0,timef,Finalenergy); //não é realmente morte *.*
			//P[i].Print(1);

	Energy=Finalenergy;
	double velocityelectron=finalvelocity;

	}

//cout << P[i-1].GetTime()<<endl;
	//P[i].Print();
	#ifdef DEBUG
	cout << __PRETTY_FUNCTION__  <<endl;
	#endif	
}



void Propagator::PositronPropagator(int i)
{
	#ifdef DEBUG
	cout << __PRETTY_FUNCTION__  <<endl;
	#endif	
	double Energy= P[i].GetEnergy();
	double rintbrem=0.;
	double rintaniq=0;
	double rint=0;
	double probability;
	double energyloss;
	double dedx;
	double y;
	double Finalenergy;
	double u=0;
	double EnergyFotao=0.;
	double EnergyFotao2=0.;

	while((Energy-masseletron)>5e-5)
	{
		//cout << "Propagando positrão com energia " << Energy << " MeV."<<endl;

		//descobrir a distancia percorrida por este fotão antes de decair
		
		probability= F->BremRadiation(Energy-masseletron); //integral de 50 eV até à energia da particula

//		cout <<"probabilidade de brem "<< probability<<endl;

		y=gRandom->Uniform(0.,1.); //gerando R interacao de acordo log(y)/-pint
		rintbrem=log(y)/(-probability);


		probability= F->AniqPositron(Energy-masseletron); //integral de 50 eV até à energia da particula
		//cout << "probabilidade de aniquilar "<<  probability <<endl;
		y=gRandom->Uniform(0.,1.);
		rintaniq=log(y)/(-probability);
		//cout << "aniq "<<rintaniq<<" brem "<< rintbrem<<endl;
		if(rintaniq<=rintbrem)
		{
		//	cout << "aniq"<<endl;
			rint=rintaniq;
		}
		if(rintbrem<rintaniq)
		{
			//cout << "brem"<<endl;
			rint=rintbrem;
		}

		//cout <<"raios de interacao " << rint <<endl;
		dedx= F->EnergyLoss(2,Energy-masseletron); 
		energyloss=rint*dedx; //energia cinetica perdida
		//cout << "energy lost ="<< energyloss << " MeV and walked " << rint << "cm" << endl;
		P[i].SetEnergyLoss(energyloss);

		Finalenergy= Energy-energyloss;
		//cout<<"Finalenergy="<<Finalenergy<<endl;
		//cout << "energia final "<<Finalenergy<<endl;
		double Theta= P[i].GetTheta(); 
		double Phi= P[i].GetPhi();

	
		//divisao de energia

		if (Finalenergy-masseletron > 5e-5)
		{
		
			if(rintaniq<=rintbrem)
			{
				//cout << "aniquilação do positrao" <<endl;
				u= F->AniqEnergy(Finalenergy-masseletron); //retorna u = energia fotao 1 / (E+m)
				EnergyFotao=u*(Finalenergy+masseletron);
				EnergyFotao2=Finalenergy+masseletron-EnergyFotao;
				//cout << "Energia positrao "<< Finalenergy-masseletron<<endl;
				//cout << "u "<< u << " Energia do fotao 1 "<< EnergyFotao  << "  e do 2 "<<EnergyFotao2<<endl;

				Finalenergy=5e-5;
			}
			if(rintbrem<rintaniq)
			{
				u=-1;
				EnergyFotao= F->EnergyBrem(Finalenergy-masseletron);

				Finalenergy-=EnergyFotao;
				//cout << "energy fotao " ;
			    //cout << EnergyFotao <<endl;
			}

		}

		if(Finalenergy-masseletron < 5e-5)//correcao no caso de ter perdido demasiada energia e deixar de ser detetado
		{
			//cout << "raio "<<rint<<endl;
			Finalenergy=5e-5+masseletron;
			energyloss= Energy-Finalenergy;
			//cout << energyloss << endl;
			rint=energyloss/dedx;
			//cout <<"raios "<<rint<<endl;
		}


		//printf("%.50f c initial velocity\n",initialvelocity);

		//printf("%.50f c final velocity\n",finalvelocity);
		double finalvelocity=sqrt((Finalenergy-masseletron)*(Finalenergy+masseletron))/Finalenergy;
		double initialvelocity=sqrt((Energy-masseletron)*(Energy+masseletron))/Energy;
		//printf("%.50f c initial velocity\n",initialvelocity);
		//printf("%.50f c final velocity\n",finalvelocity);
		double timef=fabs(P[i].GetTime() + rint*1e7*log(initialvelocity/finalvelocity));
		//cout << timef << " nano s percorridos " <<endl;

		//descobrir posição final
		Vec V0(P[i].GetFPositions()); 
	

		//coordinates in the lab referencial
		V0[0]=V0[0]+ rint*sin(Phi)*cos(Theta);
		V0[1]=V0[1]+ rint*sin(Phi)*sin(Theta);
		V0[2]=V0[2]+ rint*cos(Phi);
		if(u!=-1)
		{
			if(EnergyFotao>5 && EnergyFotao2>5)
			{
				Particle Pf1(0,EnergyFotao,V0,timef);
				Particle Pf2(0,EnergyFotao2,V0,timef);
				F->AniqAngles(Theta,Phi,Pf1,Pf2,u,0);
				P.push_back(Pf1);
				P.push_back(Pf2);
				EnergyFotao==0;
			}
			if (EnergyFotao2>5)
			{
				Particle Pf1(0,EnergyFotao2,V0,timef);
				F->AniqAngles(Theta,Phi,Pf1,Pf1,u,1);
				P.push_back(Pf1);
				EnergyFotao==0;
			}
				
		}
		if (EnergyFotao>5 )
		{
		//	cout << "Energy Fotao "<<EnergyFotao<<endl;
			Particle Pf(0,EnergyFotao,V0,timef);
		//	cout << "energia "<<Pf.GetEnergy() << " a entrar em bremangles "<<endl;
			F->BremAngles(Theta,Phi,Pf);
			//cout << Pf.GetEnergy() <<endl;
			P.push_back(Pf);
		}
		//cout << Finalenergy <<endl;
		//cout << timef<<endl;
		
	

	P[i].SetDeath(V0,timef,Finalenergy); //não é realmente morte *.*
		//	P[i].Print(1);

	Energy=Finalenergy;

	}


	//P[i].Print();
	//cout << P[i].GetTime()<< " positrao"<<endl;
	#ifdef DEBUG
	cout << __PRETTY_FUNCTION__  <<endl;
	#endif	
}