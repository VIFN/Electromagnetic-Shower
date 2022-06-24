#include "cFCgraphics.h"
#include "TLine.h"
#include "iostream"
#include "TFrame.h"
#include "TH1F.h"
#include <fstream>
#include "TStyle.h"
using namespace std;

int main()
{


	double Energy=1000;
	#ifdef E3
	Energy=3000;
	#endif
	#ifdef E5
	Energy=5000;
	#endif
	double x1;
	double z1;
	double y1;
	double x2;
	double y2;
	double z2;
	int type;
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
	if (!fs.is_open())
	{
		cout << "FILE NOT FOUND!!! "<<endl;
		return 1;
	}
	double limitx=0;
	double limitz=0;
	int count=0;
	while(fs>>type>>x1>>y1>>z1>>x2>>y2>>z2)
	{
		count++;
		if(limitx<fabs(x2))
			limitx=fabs(x2);
		if(limitz<fabs(z2))
			limitz=fabs(z2);
		if(limitx<fabs(x1))
			limitx=fabs(x1);
		if(limitz<fabs(z1))
			limitz=fabs(z1);
	}
	cout << count <<endl;
	fs.close();


	cFCgraphics G;
	TPad *pad1=G.CreatePad("pad");
	cout << limitx << " z "<< limitz <<endl;
	TH1F *frame= pad1->DrawFrame(0,-limitx,limitz,limitx);
	frame->SetTitle("Eletromagnetic Shower");
	frame->GetXaxis()->SetTitle("z (cm)");
	frame->GetYaxis()->SetTitle("x (cm)");
	frame->SetLineColor(0);
	
	TLine **lines= new TLine*[count];
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
	int i=0;
	int firstfotao=0;
	int firsteletrao=0;
	int firstpositrao=0;

	
	while(fs>>type>>x1>>y1>>z1>>x2>>y2>>z2)
	{
		lines[i]= new TLine(z1,x1,z2,x2);
		if (type==0)
		{
			lines[i]->SetLineColor(kRed);
			lines[i]->SetLineStyle(2);
			lines[i]->SetLineWidth(4);
		}
		if (type==1)
		{
			lines[i]->SetLineColor(kBlue);
			lines[i]->SetLineStyle(1);
			lines[i]->SetLineWidth(4);
		}
		if (type==2)
		{
			lines[i]->SetLineColor(kGreen+2);
			lines[i]->SetLineStyle(1);
			lines[i]->SetLineWidth(4);
		}
		cout << "z1 "<<z1<< " z2 "<< z2<< " x1 "<< x1<< " x2 "<< x2<<endl;
		++i;
	}
	/*TLine *photon= new TLine(10,0,15,0.2);
	photon->SetLineWidth(2);
	photon->SetLineStyle(4);
	G.AddObject(photon,"pad1");
	G.DrawPadFlush("pad1");*/
	for(int j=0;j<i;++j)
		G.AddObject(lines[j],"pad1");
	fs.close();

	G.DrawPad("pad1");
	G.Print("DrawFrame.png");
	return 0;
}