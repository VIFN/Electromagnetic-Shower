#include "TPolyLine3D.h"
#include "TView.h"
#include "TAxis.h"
#include <fstream>
#include "TGraph2D.h"
#include "cFCgraphics.h"
#include <iostream>
#include "TFrame.h"
#include "TLegendEntry.h"
#include "TLegend.h"
#include "TAxis3D.h"
#include "TView.h"
#include "TView3D.h"
#include "TStyle.h"
#include <cmath>
#include "TPaveText.h"
using namespace std;
int main()
{
	cFCgraphics G;
	TPad *pad1= G.CreatePad("pad1");
	double Energy=1000;
	#ifdef E3
	Energy=3000;
	#endif
	#ifdef E5
	Energy=5000;
	#endif
	double x1;
	double y1;
	double z1;
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
	int count=0;
	while(fs>>type>>x1>>y1>>z1>>x2>>y2>>z2)
		count++;
	cout << count <<endl;
	fs.close();

	TPolyLine3D **lines= new TPolyLine3D*[count+3];
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
	double limitx=0;
	double limity=0;
	while(fs>>type>>x1>>y1>>z1>>x2>>y2>>z2)
	{
		lines[i]= new TPolyLine3D(1);
		lines[i]->SetPoint(0,x1,y1,z1);
		lines[i]->SetPoint(1,x2,y2,z2);
		if(limitx<fabs(x2))
			limitx=fabs(x2);
		if(limity<fabs(y2))
			limity=fabs(y2);
		if(limitx<fabs(x1))
			limitx=fabs(x1);
		if(limity<fabs(y1))
			limity=fabs(y1);
		if (firstfotao==0 && type==0)
			firstfotao=i;
		if (firsteletrao==0 && type==1)
			firsteletrao=i;
		if (firstpositrao==0 && type==2)
			firstpositrao=i;
		
		if (type==0)
		{
			lines[i]->SetLineColor(kRed);
			lines[i]->SetLineWidth(1);
			
		}
		if (type==1)
		{
			lines[i]->SetLineColor(kBlue);
			lines[i]->SetLineStyle(1);
			lines[i]->SetLineWidth(2);
			lines[i]->SetLineStyle(3);
		}
		if (type==2)
		{
			lines[i]->SetLineColor(kGreen+1);
			lines[i]->SetLineStyle(1);
			lines[i]->SetLineWidth(1);
		}
		G.AddObject(lines[i],"pad1");
		++i;
	}

	lines[i+1]= new TPolyLine3D(4);//limites laterais
	lines[i+1]->SetPoint(0,limitx,-limity,0);
	lines[i+1]->SetPoint(1,-limitx,-limity,0);
	lines[i+1]->SetPoint(2,-limitx,limity,0);
	lines[i+1]->SetPoint(3,limitx,limity,0);
	lines[i+1]->SetPoint(4,limitx,-limity,0);
	lines[i+1]->SetLineColor(kOrange);
	lines[i+1]->SetLineStyle(9);
	G.AddObject(lines[i+1],"pad1");//linhas transversais e diagonais
	lines[i+2]= new TPolyLine3D(8);
	lines[i+2]->SetPoint(0,0,-limity,0);
	lines[i+2]->SetPoint(1,0,limity,0);
	lines[i+2]->SetPoint(2,limitx,limity,0);
	lines[i+2]->SetPoint(3,limitx,0,0);
	lines[i+2]->SetPoint(4,-limitx,0,0);
	lines[i+2]->SetPoint(5,-limitx,limity,0);
	lines[i+2]->SetPoint(6,limitx,-limity,0);
	lines[i+2]->SetPoint(7,limitx,limity,0);
	lines[i+2]->SetPoint(8,-limitx,-limity,0);
	lines[i+2]->SetLineColor(kOrange);
	lines[i+2]->SetLineStyle(9);
	G.AddObject(lines[i+2],"pad1");


	fs.close();
	TPaveText *text= new TPaveText(0.1,0.81,0.9,0.97);
	text->AddText("Cascata eletromagnetica no aluminio");
	#ifdef IRON 
	text->Clear();
	text->AddText("Cascata eletromagnetica numa placa de Ferro (Z=26)");
	#endif

	#ifdef CHUMBO 
	text->Clear();
	text->AddText("Cascata eletromagnetica numa placa de Chumbo (Z=82)");
	#endif

	#ifdef GRAFITE 
	text->Clear();
	text->AddText("Cascata eletromagnetica numa placa de Grafite (Z=6)");
	#endif
	text->SetFillColor(10);
	string nameEnergy("Energia= ");
	nameEnergy += std::to_string((int)Energy);
	nameEnergy += " MeV";
	text->AddText(nameEnergy.c_str());
	G.AddObject(text,"pad1");
	//legend
	TLegend *legend=new TLegend(0.12,0.75,0.4,0.98,"legenda","brNDC");
	legend->SetTextSize(0.03);
	legend->SetBorderSize(0);
	legend->SetTextAlign(11);
	legend->SetTextFont(1);
	legend->SetHeader("Legenda");
	legend->AddEntry(lines[firstfotao],"Fotao","l");
	legend->AddEntry(lines[firsteletrao],"Eletrao","l");
	legend->AddEntry(lines[firstpositrao],"Positrao","l");
	G.AddObject(legend,"pad1");

 //axis
	TAxis3D *axis= new TAxis3D();
	//axis->SetAxisRange(-limit,limit,"Y");
	//axis->SetAxisRange(-limit,limit,"X");
	axis->GetXaxis()->SetTitle("x (cm)");
	axis->GetYaxis()->SetTitle("y (cm)");
	axis->GetYaxis()->CenterTitle();
	axis->GetXaxis()->CenterTitle();
	axis->GetZaxis()->SetTitle("z (cm)");
	axis->GetXaxis()->SetTitleFont(1);
	axis->GetYaxis()->SetTitleFont(1);
	axis->GetZaxis()->SetTitleFont(1);
	axis->GetXaxis()->SetTitleOffset(1.3);
	axis->GetYaxis()->SetTitleOffset(1.3);
	axis->GetZaxis()->SetTitleOffset(1.3);
	axis->GetXaxis()->SetAxisColor(kBlue);
	axis->GetYaxis()->SetAxisColor(kBlue);
	axis->GetZaxis()->SetAxisColor(kBlue);
	axis->GetXaxis()->SetLabelColor(kBlack);
	axis->GetYaxis()->SetLabelColor(kBlack);
	axis->GetZaxis()->SetLabelColor(kBlack);

	//axis->SetDrawOption("AP");
	G.AddObject(axis,"pad1");
	G.AddObject(pad1);

	G.Draw();

	G.Print("result.png");

	#ifdef IRON 
	G.Print("iron.png");
	#endif
	#ifdef CHUMBO 
	G.Print("lead.png");
	#endif
	#ifdef GRAFITE 
	G.Print("graphite.png");
	#endif
	return 0;
}