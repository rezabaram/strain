#ifndef STRAIN_H
#define STRAIN_H
#include <vector>
#include"parameters.h"

using namespace std;

class strain{
	public:
	strain(int i=-1, strain *f=NULL){
		ID=i; 
		N=0.0; 
		father=f;
		fitness=0.0;
		for(int i=0;i<rmax;i++){
			M[i]=0.0;
		}
	}
	double sumM(){
		double sum=0;
		for(int i=0; i<rmax; i++){
			sum+=M[i];
		}
		return sum;
	}
	double fitness;
	double N;
	double M[rmax];
	int ID;
	strain* father;
	vector<strain*> daughters;
};

#endif
