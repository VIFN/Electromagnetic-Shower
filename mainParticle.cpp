#include "Propagator.h"
#include "Particle.h"
#include <vector>
#include "Vect.h"
#include <iostream>
#include <fstream>
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
	Propagator P(1,Energy);
cout << "Done propagating"<<endl;
int i;
	vector<Particle> Results= P.GetResults(i);
	cout << Results[0].GetType() << " size "<< i<<endl;
	ofstream fs;

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
	if(!fs  )
	{
		cout << "cannot open the output file!"<<endl;
		return 1;
	}
	Vec posinicial(3,0.);
	Vec posfinal(3,0.);
	for(int j=0;j<i;++j)
	{
		posinicial=Results[j].GetPositions();
		posfinal=Results[j].GetFPositions();
		fs<<Results[j].GetType()<<" "<< posinicial[0]<<" "<< posinicial[1]<< " "<<posinicial[2]<<" "<<posfinal[0]<<" "<<posfinal[1]<<" "<<posfinal[2]<<endl;
	}
	fs.close();
Results.clear();
return 0;
}