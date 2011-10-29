#include <iostream>
#include "strain.h"


int stotal=0;
vector<strain> strains;
strain s(stotal);

void Initial_Conditions(){
	stotal=0;
	stotal++;
	s.N=10;
	//strains.push_back(s);
}

void Update(){
	s.fitness=f0*(1-beta0*s.sumM());
	s.N=s.N*(1+s.fitness*dt);
	//for(int i=0; i<rmax; i++){
		s.M[0]+=nu*s.N*dt;
		//cerr<< s.N <<endl;
	//}
}

double t;
void Run(){
	for(t=0; t<=tMax; t+=dt){
		Update();
		cout<< t<<"    "<< s.N <<"    "<< s.sumM() <<endl;
		}
}

int main(){
	Initial_Conditions();
	Run();
return 0;
}
