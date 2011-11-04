#ifndef STRAIN_H
#define STRAIN_H
#include <vector>
#include"parameters.h"

using namespace std;

class CStrain{
	public:
	CStrain(int i=-1, CStrain *f=NULL);
	~CStrain();
	double sumM();
	double WeightedSumM(double *chi);
	int count_neigh();
	void die();
	CStrain *father(){return neighbours.at(0);}
	void add_neighbour(CStrain *ps){return neighbours.push_back(ps);};
	void get_infected(double *sumN, int distance=0, CStrain *exclude=NULL);

	double fitness;
	double N;
	int ID;
	bool dead;
	double *M;
	vector<CStrain*> neighbours;
};

void CStrain::get_infected(double *sumN, int distance, CStrain *exclude){
	sumN[distance]+=N;
	if(distance==rmax) return;

	for(int i=0; i<neighbours.size(); i++){
		if(neighbours.at(i)==exclude) continue;
		neighbours.at(i)->get_infected(sumN, distance+1, this);
	}
	return;
}


int CStrain::count_neigh(){
	int count=0;
	for (int i=0; i<neighbours.size();i++){
		if(! neighbours.at(i)->dead) count++;
	}
	return count;
}



CStrain::~CStrain(){
	delete[] M;
}

CStrain::CStrain(int i, CStrain *f){
	M=new double[rmax+1];
	ID=i; 
	N=0.0; 
	if(f!=NULL){
		neighbours.push_back(f);
		f->add_neighbour(this);
	}

	fitness=0.0;
	for(int i=0;i<=rmax;i++){
		M[i]=0.0;
	}
	dead=false;
}
double CStrain::sumM(){
	double sum=0.0;
	for(int i=0; i<=rmax; i++){
		sum+=M[i];
	}
	return sum;
}

double CStrain::WeightedSumM(double *chi){
	double sum=0.0;
	for(int i=0; i<=rmax; i++){
		sum+=M[i]*chi[i];
	}
	return sum;
}

void CStrain::die(){
	N=0.0;
	dead=true;
	delete[] M;
}



#endif
