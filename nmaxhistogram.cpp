#include "math.h"
#include <fstream>
#include "TH1F.h"
#include "cFCgraphics.h"
#include <iostream>
#include <string>
using namespace std;

int main() 
{
	ifstream fs;
	cFCgraphics G;
	TPad *pad1=G.CreatePad("pad1");

	double Energy=3000;
	int count=0;
	int limitlow=1000;
	int limithigh=0;
	int nmax=0;
	fs.open("nmax.data");
	while(fs>>nmax)
	{ 
		count++;
		if (limitlow>nmax)
			limitlow=nmax;
		if(limithigh<nmax)
			limithigh=nmax;
	}
	fs.close();
	cout << "Foram propagadas "<< count << " particulas"<<endl;
	double N_Heitler=round(Energy/42.7); //em Mev 
	string namemaxHeitler("N_{max} de part#acute{i}culas para uma cascata eletromagn#acute{e}tica, E=3 GeV, <N_{max}>=");
	namemaxHeitler += std::to_string((int) N_Heitler);
	fs.open("nmax.data");
	TH1F*h = new TH1F("N_{max}", namemaxHeitler.c_str(),  (limithigh-limitlow), limitlow, limithigh);
	while(fs>>nmax)
	{ 
		h->Fill(nmax);
	}
	h->GetXaxis()->SetTitle("N_{max} (cm)");
	h->GetYaxis()->SetTitle("events");
	fs.close();

	
	G.AddObject(h,"pad1");
	G.AddObject(pad1);
	G.Draw();
	G.Print("nmax.png");
}