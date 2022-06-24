#include "Propagator.h"
#include "Particle.h"
#include <vector>
#include "Vect.h"
#include <iostream>
#include "TH1F.h"
#include "cFCgraphics.h"
#include <fstream>
using namespace std;
int main()
{	
	double Energy=3000;



	ofstream fs;
	fs.open("nmax.data");
	double step=1;
		Vec posinicial(3,0.);
		Vec posfinal(3,0.);
	int *n_max1= new int[2000];
	vector<Particle> Results;
	for (int i = 0; i < 2000 ; ++i)
	{
		Propagator P(1,Energy);
		
		int k;
		Results= P.GetResults(k);
		TH1F*h = new TH1F("part#acute{i}culas", "Numero total de particulas por unidade de comprimento",  (200)/step, 0., 200);

	
		for(int j=0;j<k;++j)
		{
			posinicial=Results[j].GetPositions();
			posfinal=Results[j].GetFPositions();
			double jj= posinicial[2];
			double z2max=(int)posfinal[2] +step;
			while(jj< z2max)
			{
				h->Fill(jj);
				jj+=step;
			}
		}
		n_max1[i]= h->GetMaximum();
		fs << n_max1[i]<<endl;
		cout << "done one "<< n_max1[i]<< " Energy" << Energy<<endl<<endl<<endl <<endl;
		Results.clear();
		delete h;
	}
	fs.close();

return 0;
}