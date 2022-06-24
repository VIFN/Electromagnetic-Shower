#include "math.h"
#include <fstream>
#include "TH1F.h"
#include "cFCgraphics.h"
#include <iostream>
#include "TROOT.h"
#include "TStyle.h"
#include <string>
#include "TPaveText.h"
using namespace std;
int main()
{
		double step=1.; // se quiser alterar numero de divisoes
	cFCgraphics G;
	TPad *pad1=G.CreatePad("pad1");
	TPad *pad2=G.CreatePad("pad2");
	TPad *pad3=G.CreatePad("pad3");
	TPad *pad4=G.CreatePad("pad4");
	double x1;
	double y1;
	double z1;
	double x2;
	double y2;
	double z2;
	double Energy=1000;
	#ifdef E3
	Energy=3000;
	#endif
	#ifdef E5
	Energy=5000;
	#endif
	int type;

	//////////////////////////////////READ FILE ////////////////////////////
	ifstream fs;
	fs.open("results.data");
	#ifdef IRON 
	fs.close();
	fs.open("iron.data");
	#endif
	#ifdef CHUMBO 
	fs.close();
	fs.open("chumbo.data");
	#endif
	#ifdef GRAFITE 
	fs.close();
	fs.open("grafite.data");
	#endif
	if (!fs.is_open()){
			cout << "FILE NOT FOUND!!! "<<endl;
		return 1;
	}
	int count=0;
	int countpos=0;
	int countfot=0;
	double limitMolier=0;
	double lim_max=0, lim_min=0;
	while(fs>>type>>x1>>y1>>z1>>x2>>y2>>z2)
	{ 
		count++;
		if (type==0)
			++countfot;
		if(type!=0)
			++countpos;
		if(z2>lim_max)
			lim_max=z2;
		if(x2>limitMolier)
			limitMolier=x2;
		if(y2>limitMolier)
			limitMolier=y2;
	
	//	cout <<type<<"    "<<x1<<"    "<<y1<<"    "<<z1<<"    "<<x2<<"     "<<y2<<"     "<<z2<<endl;
	}
	fs.close();
	////////////////////////////////// RAIO MÒLIER ////////////////////
	//Molière radius eq33.38 in PASSAGE OF PARTICLES TROUGH MATTER (Rm=X0 Es / Ec) Es=21.2052
	//utilizámos os valores experimentais tabelados
	double RMolier= 4.419; //esta analise requer Energia Elevada para ser corretamente analisada
	//double RMolier=24.01*21.2052/43.;
	#ifdef IRON 
	RMolier= 1.719;
	//RMolier=13.84*21.2052/21.68;
	#endif

	#ifdef CHUMBO
	RMolier=1.602;
	//RMolier=6.37*21.2052/7.43;
	#endif

	#ifdef GRAFITE 
	RMolier=5.012;
	//RMolier=42.70*21.2052/81.74;
	#endif
	cout <<"Largura máxima do chuveiro (Mòlier) deveria ser "<< RMolier <<" cm e foi "<<limitMolier<<" cm"<< endl;
	cout << "É expectável o valor obtido não ser concordante com o teórico já que nao foram considerados desvios dos eletroes e positroes após uma interacao de bremsstrahlung"<<endl;
	//////////////////////////////////////Xmax medio ////////////////////////
	// <Xmax>= X0 ln(E0/Ec) +0.5 (artigo na bibliografia) in Air showers and Hadronic Interactions pag.5
	double Xmax_teorico=  24.01*log(Energy/43)/log(2);
	#ifdef IRON 
	Xmax_teorico= 13.84*log(Energy/21.68)/log(2);
	#endif
	#ifdef CHUMBO 
	Xmax_teorico=6.37*log(Energy/7.43)/log(2);
	#endif

	#ifdef GRAFITE 
	Xmax_teorico=42.7*log(Energy/81.74)/log(2);
	#endif
	cout <<"Profundidade máxima do chuveiro (Rossi) deveria ser "<< Xmax_teorico <<" cm e foi "<<lim_max<<" cm"<< endl;

//////////	//////////////////////HISTOGRAMAS ///////////////////////////////////////
	double max;
	modf(lim_max, &max);

	//cout<<"lim_max="<<max<<endl;
	int m=(max+1)/step;
	cout<<"Total de particulas "<<count<< " . Total de fotoes "<< countfot << " e total de positroes/eletroes "<<countpos<<endl;
cout<<"m="<<m<<endl;
	string namefotoes("Numero de fot#tilde{o}es (");
	namefotoes += std::to_string((int) countfot);
	namefotoes += ") por unidade de comprimento";
	string namepos("Numero de eletr#tilde{o}es/positr#tilde{o}es (");
	namepos += std::to_string((int) countpos);
	namepos += ") por unidade de comprimento";
	TH1F*h = new TH1F("part#acute{i}culas", "Numero total de particulas por unidade de comprimento",  (max+1)/step, 0., max+1);
	TH1F*f = new TH1F("fot#tilde{o}es", namefotoes.c_str(),  (max+1)/step, 0., max+1);
	TH1F*g = new TH1F("eletr#tilde{o}es+positr#tilde{o}es", namepos.c_str(),  (max+1)/step, 0., max+1);

	fs.open("results.data");
	#ifdef IRON 
	fs.close();
	fs.open("iron.data");
	#endif
	#ifdef CHUMBO 
	fs.close();
	fs.open("chumbo.data");
	#endif
	#ifdef GRAFITE 
	fs.close();
	fs.open("grafite.data");
	#endif
	while (fs>>type>>x1>>y1>>z1>>x2>>y2>>z2) 
	{
		double i= z1;
		double z2max=(int)z2 +step;
		while(i< z2max)
			{
				h->Fill(i);
			    if(type==0)
					f->Fill(i);
				if(type!=0)
					g->Fill(i);
				i+=step;
			}
		
	}
	fs.close();
	int n_max1=h->GetMaximum();
	double z_max1=h->GetMaximumBin()/step;
	
	cout<<endl<<endl<<"Distância z para a qual o perfil longitudinal é máximo="<<z_max1<<endl;
	cout<<"Numero de particulas correspondentes a esta distância="	<<n_max1<<endl<<endl;


//////////////////////	//Modelo de Heitler //////////////////////7
	  
	double N_Heitler=round(Energy/42.7); //em Mev 
	#ifdef IRON 
	N_Heitler=round(Energy/21.68); //Ec tabelado
	#endif

	#ifdef CHUMBO 
	N_Heitler=round(Energy/7.43);
	#endif

	#ifdef GRAFITE 
	N_Heitler=round(Energy/81.74);
	#endif
	cout<<"Energy="<<Energy<<endl;
	if(N_Heitler ==n_max1){cout<<"Valores coincidentes com o metodo de Heitler para o numero de particulas no maximo "<<endl;}
	else
	{cout<<"Valor nao coincidente com o expectavel atraves do metodo de Heitler. Diferenca de "<<fabs(N_Heitler-n_max1)<<"  particulas"<<endl;}



/////////////////////////////// PARTE GRAFICA ///////////////////////////7
	string namemax("N_{max} experimental de "); //making strings
	namemax += std::to_string(n_max1);
	string namemaxHeitler("N_{max} de Heitler ");
	namemaxHeitler += std::to_string((int) N_Heitler);
	TPaveText *text= new TPaveText(0.01,0.75,0.99,0.97);
	string nameEnergy("Gerado por 1 fot#tilde{a}o com Energia= ");
	nameEnergy += std::to_string((int)Energy);
	nameEnergy += " MeV";
	string nameXmaxteo("<X_{max}>= ");
	nameXmaxteo += std::to_string((int)Xmax_teorico);
	nameXmaxteo += " cm ";
	nameXmaxteo += "X_{max_{exp}}= ";
	nameXmaxteo += std::to_string((int)lim_max);
	nameXmaxteo += " cm";
	string nameRmteo("<R_{M}>= ");
	nameRmteo += std::to_string((int)RMolier);
	nameRmteo += " cm ";
	nameRmteo += "R_{M_{exp}}= ";
	nameRmteo += std::to_string((int)limitMolier);
	nameRmteo += " cm";
	string nameNum("");
	nameNum += std::to_string((int)countfot);
	nameNum += " Fot#tilde{o}es e ";
	nameNum += std::to_string((int)countpos);
	nameNum += " Eletr#tilde{o}es/Positr#tilde{o}es";
	text->AddText("Cascata eletromagn#acute{e}tica numa placa de Alum#acute{i}nio");
	#ifdef IRON 
	text->Clear();
	text->AddText("Cascata eletromagn#acute{e}tica numa placa de Ferro");
	#endif

	#ifdef CHUMBO 
	text->Clear();
	text->AddText("Cascata eletromagn#acute{e}tica numa placa de Chumbo");
	#endif

	#ifdef GRAFITE 
	text->Clear();
	text->AddText("Cascata eletromagn#acute{e}tica numa placa de Grafite");
	#endif
	
	text->AddText(nameEnergy.c_str());
	text->SetFillColor(10);
	TPaveText *text2= new TPaveText(0.1,0.51,0.9,0.67);
	TPaveText *text3= new TPaveText(0.1,0.40,0.9,0.47);
	TPaveText *text4= new TPaveText(0.1,0.30,0.9,0.37);
	TPaveText *text5= new TPaveText(0.1,0.11,0.9,0.27);
	text2->AddText(namemax.c_str());
	text2->AddText(namemaxHeitler.c_str());
	text3->AddText(nameXmaxteo.c_str());
	text4->AddText(nameRmteo.c_str());
	text5->AddText(nameNum.c_str());
	text2->SetFillColor(10);
	text3->SetFillColor(10);
	text4->SetFillColor(10);
	text5->SetFillColor(10);
	G.AddObject(text,"pad4");
	G.AddObject(text2,"pad4");
	G.AddObject(text3,"pad4");
	G.AddObject(text4,"pad4");
	G.AddObject(text5,"pad4");
	h->GetXaxis()->SetRangeUser(0, max+1);
	h->GetXaxis()->SetTitle("cm");
	h->GetYaxis()->SetTitle("Numero total de particulas");
	h->SetLineWidth(3);
	h->SetLineColor(6);
	f->GetXaxis()->SetRangeUser(0, max+1);
	f->GetXaxis()->SetTitle("cm");
	f->GetYaxis()->SetTitle("Numero de fot#tilde{o}es");
	f->SetLineWidth(3);
	f->SetLineColor(1);
	g->GetXaxis()->SetRangeUser(0, max+1);
	g->GetXaxis()->SetTitle("cm");
	g->GetYaxis()->SetTitle("Numero de positr#tilde{o}es e eletr#tilde{o}es");
	g->SetLineWidth(3);
	g->SetLineColor(9);
	G.AddObject(h, "pad1");
	G.AddObject(f, "pad2");
	G.AddObject(g, "pad3");
	G.AddObject(pad1); //para ter total
	//G.AddObject(pad2); //para ter fotoes apenas
	//G.AddObject(pad3); //para ter positroes e eletroes
	G.AddObject(pad4); //para ter caixas de texto
	G.Draw();

	///////////////////////PRINTS////////////////////
	G.Print("histogramaAl.png");
	#ifdef IRON 
	G.Print("histogramaFe.png");
	#endif

	#ifdef CHUMBO

	G.Print("histogramaPb.png");
	#endif

	#ifdef GRAFITE 
	G.Print("histogramaC.png");
	#endif
	return 0;
}