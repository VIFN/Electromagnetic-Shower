#include "Formula.h"
#include <ctime>
#include <cmath>
#include "TRandom3.h"
#include "TF1.h"
#include "Vect.h"
#include "FCmatrixFull.h"
#include "TFormula.h"
#include "Vect.h"
#include "TMath.h"

#include "cFCgraphics.h"
#include "TAxis.h"
#include <iostream>
using namespace std;

 double X0= 24.01;  //comprimento de radiação associado ao material (g/cm²)
 double Z=13.;
const double c=299792458; //velocidade da luz
const double hr=pow(6582119514,-22); //planck constante, reduzida
const double e=pow(1.6021766208,-19); //carga eletrao

const double pire=0.1535; //constante 2Pi Na re*re me c*c MeV.cm^2 /g
const double ALPHA=137.035999139; //inverso da constante de estruturas finas
const double re= 2.81794*pow(10,-13); //raio classico do eletrao em cm
const double masseletron= 0.5109989461;
 double densidade=2.699; //g/cm³ densidade do aluminio
 double A=26.9815385; //massa atomica g mol-1 
const double Na=6.022140857e23; //constante avogadro mol-1
 double n=Na*densidade/A;
const double Wc=5e-5;
double Rmech=44.503;


Formula::Formula(): Emin(0.),Emax(1.),a(Z/ALPHA),g0(0.),I(0.),aa(0.0802),C0(-4.24),m(3.63),U1(3.01),U0(0.1708),delta(0.),Gamma(0.),bethagamma(0.),B(0.),C(0.),fm(0.),PMAX(0.),y(0.),xr(0.),auxA(0.),ur(0.),pos_aux(0.,0,0.),pos(3,0.),mat(3,3),   Umax(0.5),Umin(0.),
	b("b", "[0]/(x*(1-x))") , 
    g1("g1", "7/3. - 2*log(1+b*b) -6*b*atan(1/b) -b*b*(4-4*b*atan(1/b) - 3*log(1+1/(b*b)) )"), g2("g2", "11/6. - 2*log(1+b*b) -3*b*atan(1/b) +0.5*b*b*(4-4*b*atan(1/b) - 3*log(1+1/(b*b)) )"),
	phi("phi", "2*( (0.5-x)*(0.5-x)  )*(g1+[0])+g2+[0]",Emin,Emax),
	f("f","phi/[0]",Emin,Emax),
	//efeito brem
	dsig("dsig","([0]/x)*((x/[1])*(x/[1])*(4*log([3]) + 2 - 2*log(1+([3]*0.5*[2] *(x/[1]) / (1-(x/[1])))*([3]*0.5*[2] *(x/[1]) / (1-(x/[1]))))- 4*([3]*0.5*[2] *(x/[1]) / (1-(x/[1])))*atan(1/([3]*0.5*[2] *(x/[1]) / (1-(x/[1])))))+(4/3.)*(1-(x/[1]))*(4*log([3]) + 7/3 - 2*log(1+([3]*0.5*[2] *(x/[1]) / (1-(x/[1])))*([3]*0.5*[2] *(x/[1]) / (1-(x/[1]))))-6*([3]*0.5*[2] *(x/[1]) / (1-(x/[1])))*atan(1/([3]*0.5*[2] *(x/[1]) / (1-(x/[1]))))-([3]*0.5*[2] *(x/[1]) / (1-(x/[1])))*([3]*0.5*[2] *(x/[1]) / (1-(x/[1])))*( 4 -4*([3]*0.5*[2] *(x/[1]) / (1-(x/[1])))*atan(1/([3]*0.5*[2] *(x/[1]) / (1-(x/[1]))))-3*log(1+1/(([3]*0.5*[2] *(x/[1]) / (1-(x/[1])))*([3]*0.5*[2] *(x/[1]) / (1-(x/[1])))) )  )) )",0,5000),

	Aniq("Aniq", "[0]*[0]/([1]+1)/([1]*[1]-1)*(-([1]+1)*([1]+1)+([1]*[1]+4*[1]+1)/x-1/x/x-([1]+1)*([1]+1)+([1]*[1]+4*[1]+1)/x-1/(1-x)/(1-x))", 0, 0.5),
	auxiliar("aux", "[0]/x",0,5000),
	anglesbrem("1", "[0]*(1/pow(1-[1]*x,2)+pow(x-[1],2)/pow(1-[1]*x,4))", -1, 1),
	auxiliarangles("auxang","2*[0]/(1-[1]*x)/(1-[1]*x)",-1,1),
	//perda de energia
	de("dE", "([0]/[1])*( log( [2]/[3]*( ([4]+1)/2 ) )+ [8] - [5] -2*[6]/[7] )")
{
	gRandom= new TRandom3(0.);  

	//outros elementos alterar X0,Z,densidade,A,n,Rmech,C0,aa,m,U1,U0
	#ifdef IRON //Z=26
	Z=26;
 	X0=13.84;
 	densidade=7.874;
	A=55.845;
	n=Na*densidade/A;
	Rmech=36.663;
	C0=-4.29;
	aa=0.1468;
	m=2.96;
	U1=3.15;
	U0=-0.0012;
	#endif

	#ifdef CHUMBO //Z=82
	Z=82;
 	X0=6.37;
 	densidade=11.35;
	A=207.2;
	n=Na*densidade/A;
	Rmech=23.922;
	C0=-6.20;
	aa=0.0936;
	m=3.16;
	U1=3.81;
	U0=0.3776;
	#endif

	#ifdef GRAFITE //Z=6
	Z=6;
 	X0=42.70;
 	densidade=2; //segundo tabela.Leo em anexo (diferente do valor da tabela no site pdg.lbl.gov)
	A=12.0107;
	n=Na*densidade/A;
	Rmech=61.228;
	C0=-2.99;
	aa=0.2024;
	m=3.00;
	U1=2.486;
	U0=-0.0351;
	#endif
	dsig.SetParameter(0,re*re*13*13/(ALPHA));
	dsig.SetParameter(3,Rmech);
	Aniq.SetParameter(0, TMath::Pi()*re);
	if (Z<13)
			I=(12+7/Z)/Z;
	else
			I= (9.76 + 58.8 * pow(Z,-1.19))*Z; //em eV (CONVERSAO NECESSARIA)

}


double Formula::PhotonEnergyDivider(double Energy)
{
	Emin=masseletron/Energy;
	Emax=1-masseletron/Energy;

	g0= 4*log(Rmech) - 4*a*a*( (1/(1+a*a)) + 0.202059 -0.03693*pow(a,2) + 0.00835*pow(a,4)*-0.00201*pow(a,6)+0.00049*pow(a,8)-0.00012*pow(a,10)+0.00003*pow(a,12) )  + (-0.1774- 12.10*a + 11.18*a*a)*( sqrt(2*Emin)) +  (8.523+ 73.26*a - 44.41*a*a)*( 2*Emin)  -  (13.52+ 121.1*a - 96.41*a*a)*( pow(2*Emin,1.5) ) + (8.946+ 62.05*a - 63.41*a*a)*( pow(2*Emin,2.) ) ;
#ifdef DEBUG
cout << g0<<endl;
#endif
	b.SetParameter(0,Rmech*0.5*Emin); //mass eletron

	#ifdef DEBUG
cout << b.Eval(0.2)<<endl;
#endif
	
	phi.SetParameter(0.,g0);
	PMAX= phi.GetMaximum();
	//cout << PMAX << endl;
	f.SetParameter(0,PMAX);
	#ifdef DEBUG
cout << "hi"<<endl;
#endif


	double R2=1.;
	double ER=0.5;
	while(R2>=f.Eval(ER))
	{
		ER=gRandom->Uniform(Emin,Emax);
		R2=gRandom->Uniform(0,1);
	}
	#ifdef DEBUG
	cout << ER<<endl;
	#endif
	return ER;

}

void Formula::PhotonAngles(double Theta,double Phi,Particle& Pe)
{
	/*double mat[3][3]={{cos(Theta)*cos(Phi), -sin(Phi),  sin(Theta)*cos(Phi)}, 
						{cos(Theta)*sin(Phi), cos(Phi), sin(Theta)*sin(Phi)}, 
							{-sin(Theta), 0, cos(Theta)}};*/
	//cout << "e ai "<<Theta << " "<<Phi <<endl;
	mat[0][0]=cos(Phi)*cos(Theta);
	mat[0][1]=-sin(Theta);
	mat[0][2]=sin(Phi)*cos(Theta);
	mat[1][0]=cos(Phi)*sin(Theta);
	mat[1][1]=cos(Theta);
	mat[1][2]=sin(Phi)*sin(Theta);
	mat[2][0]=-sin(Phi);
	mat[2][2]=cos(Phi);
	y=gRandom->Uniform(0, 1);
	 //Energia cinetica do eletrao
	Emin=Pe.GetEnergy()-masseletron;//energia 

    //cout << "energy "<< Emin <<endl;
	betha=sqrt(Emin*Emin+2*Emin*masseletron)/(Emin+masseletron); 
	//cout << betha <<endl;
    delta= (-1+betha*betha)/2; //constante de normalizacao de 1/(1-bx)^2
    //delta corresponde ao desvio da particula em relacao a direcao do fotao
    delta= (delta-betha*y+y)/(delta-betha*y*(betha-1));
    //cout << delta<<endl;
    if(isnan(delta)==1)
    	delta=1;
    C=gRandom->Uniform(0, 2*(TMath::Pi()));
    //cout<<"Azimutal angle"<<C<<"      delta="<<delta<<endl;
	pos[0]=(Emin+masseletron)*sqrt(1-delta*delta)*cos(C);
	pos[1]=(Emin+masseletron)*sqrt(1-delta*delta)*sin(C);
	pos[2]=(Emin+masseletron)*delta;
	//pos.Print();
	//	Pe.Print(1);

	//Pe.SetTheta(C);
	//Pe.SetPhi(delta);
	//cout << "cheguei"<<endl;

	
	//pos.Print();
	pos=(mat)*pos;
	//pos.Print();

	Pe.SetXYZ(pos[0],pos[1],pos[2]);
 	//Pe.Print(1);
 	//cout << "theta "<<Pe.GetTheta()<< " phi "<< Pe.GetPhi()<<endl;
}



double Formula::BremRadiation(double Energy)
{ //emin representa probabilidade de aniq, apenas para nao declarar nova



	dsig.SetParameter(1,n*(Energy+masseletron));
	dsig.SetParameter(2,masseletron/(masseletron+Energy ) );

	Emin= adaptiveSimpsons(Wc,Energy,1e-27,50); //integral de 50 eV at
	//cout << Emin*n << "probs"<<endl;
	/*double test5= dsig.Eval(1);
	printf("%.50f dsig eval em 1 para energia 1000\n",test5);
	cout <<"integral "<< integral <<endl;*/

	return Emin*n;
}

double Formula::AniqPositron(double Energy)
{ //emin representa probabilidade de aniq, apenas para nao declarar nova
	Gamma= 1 + Energy/masseletron;
	Emin=TMath::Pi()*re*re/((Gamma+1)*(Gamma*Gamma-1));
	Emin= Emin*( (Gamma*Gamma+4*Gamma+1)*log(Gamma+sqrt(Gamma*Gamma-1)) - (3+Gamma)*sqrt(Gamma*Gamma-1) );
	return Emin*n*Z;
}

double Formula::EnergyBrem(double Energy)
{
//utilizou-se uma função auxiliar da forma A/x com eficiencia de 93 a 96 %

 	dsig.SetParameter(1,Energy+masseletron);
	dsig.SetParameter(2,masseletron/(masseletron+Energy ) );

auxiliar.SetParameter(0,Wc*dsig.Eval(Wc));

 while(1)
 {
	 y=gRandom->Uniform(0.,1.);
	 xr= Wc * pow(Energy/Wc,y);
	 ur=gRandom->Uniform();
	 if(ur<=dsig.Eval(xr)/auxiliar.Eval(xr))
	 {
	 	
	 	return xr;
	 }
	
 }

}

double Formula::AniqEnergy(double Energy)
{
	Gamma= 1 + Energy/masseletron;
	Umin=1/(1+Gamma+(sqrt(Gamma*Gamma-1)));
	Aniq.SetParameter(1, Gamma);
	auxiliar.SetParameter(0,Umin*Aniq.Eval(Umin));

	 while(1)
	 {
		 y=gRandom->Uniform(0.,1.);
		 xr= Umin * pow(Umax/Umin,y);
		 ur=gRandom->Uniform();
		 if(ur<=Aniq.Eval(xr)/auxiliar.Eval(xr))
		 {
		 	
		 	return xr;
		 }
		
	 }

}


void Formula::BremAngles(double Theta,double Phi, Particle& Pe)
{
	/*double mat[3][3]={{cos(Theta)*cos(Phi), -sin(Phi),  sin(Theta)*cos(Phi)}, 
						{cos(Theta)*sin(Phi), cos(Phi), sin(Theta)*sin(Phi)}, 
							{-sin(Theta), 0, cos(Theta)}};*/
	//cout << "e ai "<<Theta << " "<<Phi <<endl;
	mat[0][0]=cos(Phi)*cos(Theta);
	mat[0][1]=-sin(Theta);
	mat[0][2]=sin(Phi)*cos(Theta);
	mat[1][0]=cos(Phi)*sin(Theta);
	mat[1][1]=cos(Theta);
	mat[1][2]=sin(Phi)*sin(Theta);
	mat[2][0]=-sin(Phi);
	mat[2][2]=cos(Phi);
	y=gRandom->Uniform(0, 1);
	 //Energia cinetica do eletrao
	Emin=Pe.GetEnergy();
	//cout << Emin << endl;
	//cout << Theta << " theta "<< Phi << " Phi "<<endl;
	delta= Formula::BremanglesAux(Emin); //cos delta obtido
	if(isnan(delta)==1)
    	delta=1;
    C=gRandom->Uniform(0, 2*(TMath::Pi()));
  // cout<<"Brem:Azimutal angle "<<C<<"      delta="<<delta<<endl;
	pos[0]=Emin*sqrt(1-delta*delta)*cos(C);
	pos[1]=Emin*sqrt(1-delta*delta)*sin(C);
	pos[2]=Emin*delta;
	//pos.Print();
	//	Pe.Print(1);

	//Pe.SetTheta(C);
	//Pe.SetPhi(delta);
	//cout << "cheguei"<<endl;

	
	//pos.Print();
	pos=(mat)*pos;
	//pos.Print();

	Pe.SetXYZ(pos[0],pos[1],pos[2]);
	//cout << Pe.GetEnergy()<<endl;
 	//Pe.Print(1);
 	//cout << "theta "<<Pe.GetTheta()<< " phi "<< Pe.GetPhi()<<endl;
}

double Formula::BremanglesAux(double Energy)
{
betha= sqrt(Energy*(Energy+2*masseletron))/(Energy+ masseletron);


anglesbrem.SetParameter(0,3./(16.*TMath::Pi())*(masseletron/(Energy+masseletron))*(masseletron/(Energy+masseletron)) );

anglesbrem.SetParameter(1,betha);
auxiliarangles.SetParameter(0,3./(16.*TMath::Pi())*(masseletron/(Energy+masseletron))*(masseletron/(Energy+masseletron)) );

auxiliarangles.SetParameter(1,betha);

 while(1)
 {
	 y=gRandom->Uniform(0.,1.);
	 xr= (betha+2*y-1)/(betha*(2*y-1)+1);
	 ur=gRandom->Uniform();
	 if(ur<=anglesbrem.Eval(xr)/auxiliarangles.Eval(xr))
	 {
	 
	 	//cout << xr<< endl;
	 	return xr;
	 }
	
 }
}

void Formula::AniqAngles(double Theta,double Phi,Particle& F1,Particle& F2,double u, int mode)
{
	mat[0][0]=cos(Phi)*cos(Theta);
	mat[0][1]=-sin(Theta);
	mat[0][2]=sin(Phi)*cos(Theta);
	mat[1][0]=cos(Phi)*sin(Theta);
	mat[1][1]=cos(Theta);
	mat[1][2]=sin(Phi)*sin(Theta);
	mat[2][0]=-sin(Phi);
	mat[2][2]=cos(Phi);
	FCmatrixFull matrixsafeguard(mat);
	y=gRandom->Uniform(0, 1);
	 //Energia cinetica do eletrao
	if(mode==0)
	{	 //indica que F1 é o fotao de energia mais baixa
		Emin=F1.GetEnergy();
		Emin=Emin/u - 2*masseletron;
		//cout << "energia cinetica do positrao "<< Emin <<endl;
	}

	if(mode==1)
	{ //indica que F1=F2 e apenas o fotao mais energetico foi aceite
		Emin=F2.GetEnergy();
		Emin=Emin/(1-u)-2*masseletron;
		//cout << "energia cinetica do positrao "<< Emin <<endl;
	}
	Gamma=1+ Emin/masseletron;

	//angulos

	if(mode==0)
	{
		delta= (1+Gamma-1/(1-u))/sqrt((Gamma)*(Gamma)-1);
		ur= (1+Gamma-1/u)/sqrt((Gamma)*(Gamma)-1);
		//cout << "angulo do fotao 1 "<< ur <<" angulo do fotao 2 "<< delta <<endl;
	}
	if (mode==1)
	{
		delta= (1+Gamma-1/(1-u))/sqrt((Gamma)*(Gamma)-1);
		//cout <<" angulo do fotao mais energetico "<< delta <<endl;
	}
	if(isnan(delta)==1)
    	delta=1;
    if(isnan(ur)==1)
    	ur=1;
    C=gRandom->Uniform(0, 2*(TMath::Pi()));

	//cout<<"angulo azimutal fotao menos energetico= "<<azimutal1<< " se mode " <<mode<< "for 1, este é o único"<<endl;
	if (mode==0)
	{
		fm=C+TMath::Pi();
		if(fm>2*TMath::Pi())
			fm-=2*TMath::Pi();
		//cout<<"angulo azimutal fotao menos energetico="<<fm << " e do mais "<<C <<endl;
	}
	Emin=F2.GetEnergy();
	pos[0]=Emin*sqrt(1-delta*delta)*cos(C);
	pos[1]=Emin*sqrt(1-delta*delta)*sin(C);
	pos[2]=Emin*delta;
	//pos.Print();
	//	Pe.Print(1);
	//pos.Print();
	pos=(matrixsafeguard)*pos;
	//pos.Print();
	F2.SetXYZ(pos[0],pos[1],pos[2]);

	if (mode==0)
	{
		Emin=F1.GetEnergy();
		pos[0]=Emin*sqrt(1-ur*ur)*cos(fm);
		pos[1]=Emin*sqrt(1-ur*ur)*sin(fm);
		pos[2]=Emin*ur;
		pos=(mat)*pos;
		//pos.Print();
		F1.SetXYZ(pos[0],pos[1],pos[2]);
		//F1.Print();
		//F2.Print();
	}
// 	F2.Print();

 	//cout << "theta "<<Pe.GetTheta()<< " phi "<< Pe.GetPhi()<<endl;
}




double Formula::EnergyLoss(int type, double Energy)
{
	#ifdef DEBUG
	cout << __PRETTY_FUNCTION__  <<endl;
	#endif
	betha= sqrt(Energy*(Energy+2*masseletron))/(Energy+ masseletron);


	Gamma=1+ Energy/masseletron;
	

	bethagamma=betha*Gamma;
	B= log10(bethagamma);


	if (B<U0)
			delta=0.;
	else if(U0<=B && B<U1)
			delta=4.6052*B+C0+aa*pow(U1-B,m);
	else if(B>=U1)
			delta= 4.6052*B+C0;


	C= 1e-6 * I*I*( (0.422377/(bethagamma*bethagamma) +0.0304043/(bethagamma*bethagamma*bethagamma*bethagamma) - 0.00038106/(bethagamma*bethagamma*bethagamma*bethagamma*bethagamma*bethagamma) ) + ( 3.850190/(bethagamma*bethagamma) +0.1667989/(bethagamma*bethagamma*bethagamma*bethagamma) - 0.00157955/(bethagamma*bethagamma*bethagamma*bethagamma*bethagamma*bethagamma) )*1e-3*I    );

	double answer=0;
	
	fm=0;

	switch(type)
	{
		case 1:			
		fm= 1- betha*betha - ( (2*Gamma-1)/(Gamma*Gamma) )*log(2) + (1/8.) * (Gamma-1)*(Gamma-1)/(Gamma*Gamma);
		break;
		case 2:
		fm= 2*log(2)- (betha*betha/12) * ( 23 + 14/(Gamma+1) + 10/pow(Gamma+1,2) + 4/pow(Gamma+1,3) ) ;
		break;
	}


			de.SetParameter(0,pire*densidade*Z);
			de.SetParameter(1,A*betha*betha);
			de.SetParameter(2,Energy*Energy);
			de.SetParameter(3,I*I*1e-12); //em ev para passar a mev
			de.SetParameter(4,Gamma);
			de.SetParameter(5,delta);
			de.SetParameter(6,C); // C era em ev
			de.SetParameter(7,Z);
			de.SetParameter(8,fm);
	answer= de.Eval(0.);

	#ifdef DEBUG
	cout << __PRETTY_FUNCTION__  <<endl;
	cout << answer<<endl;
	#endif
	return answer;
		
}

double Formula::adaptiveSimpsonsAux(double a, double b, double epsilon,                 
                         double S, double fa, double fb, double fc, int bottom) {                 
  double c = (a + b)/2;
  double h = b - a;                                                                  
  double d = (a + c)/2;
  double e = (c + b)/2;                                                              
  double fd = dsig.Eval(d);
  double fe = dsig.Eval(e);                                                     
  double Sleft = (h/12)*(fa + 4*fd + fc);                                                           
  double Sright = (h/12)*(fc + 4*fe + fb);                                                          
  double S2 = Sleft + Sright;                                                                       
  if (bottom <= 0 || fabs(S2 - S) <= 15*epsilon)   // magic 15 comes from error analysis                                       
    return S2 + (S2 - S)/15;                                                                        
  return adaptiveSimpsonsAux(a, c, epsilon/2, Sleft,  fa, fc, fd, bottom-1) +                    
         adaptiveSimpsonsAux(c, b, epsilon/2, Sright, fc, fb, fe, bottom-1);                     
}

double Formula::adaptiveSimpsons(   // ptr to function
                           double a, double b,  // interval [a,b]
                           double epsilon,  // error tolerance
                           int maxRecursionDepth) {   // recursion cap        
  double c = (a + b)/2;
  double h = b - a;                                                                  
  double fa = dsig.Eval(a); 
  double fb = dsig.Eval(b);
  double fc = dsig.Eval(c);                                                           
  double S = (h/6)*(fa + 4*fc + fb);                                                                
  return adaptiveSimpsonsAux( a, b, epsilon, S, fa, fb, fc, maxRecursionDepth);                   
}  