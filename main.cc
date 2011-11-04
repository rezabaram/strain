#include <iostream>
#include <list>
#include <tr1/random>
#include "strain.h"
#include <ctime>
std::tr1::ranlux64_base_01 eng;
std::tr1::uniform_real<double> unif(0, 1);




int stotal=0;
list<CStrain*> strains;
CStrain *top=NULL;
double t;
double chi[rmax+1];

void define_cross_im(){


	double a=1;
	for(int i=0; i<=rmax; i++){
	chi[i]=a;
	a-=1.0/(rmax+1);
	cout << chi[i] << "    ";
	}
	cout << endl;

/*
	chi[0]=0.1;
	chi[1]=0.9;
	chi[2]=0.85;
	chi[3]=0.3;
	chi[4]=0.1;
	chi[5]=0.01;
	chi[5]=0.5;
	chi[6]=0.4;
	chi[7]=0.0;
	chi[8]=0.0;
	chi[9]=0.0;
	chi[10]=0.0;
*/

}


void Initial_Conditions(){
	eng.seed(time(0));
	stotal=0;
	top=new CStrain(stotal);
	stotal++;
	top->N=10;
	strains.push_back(top);
	define_cross_im();
}

void Mutate(CStrain *pfather){
	CStrain *ps =new CStrain(stotal,pfather);
	ps->N=1;

	strains.push_back(ps);

	stotal++;
}

void Immune_Selection(){
	list<CStrain*>::iterator it;
	for(it=strains.begin(); it!=strains.end(); it++){
	//applying to all strains
		CStrain *s=(*it);
		s->fitness=f0*(1-beta0*s->WeightedSumM(chi) );
		s->N=s->N*(1+s->fitness*dt);//make sure about the order of update N and M
	}
}

void Genetic_Drift(){
	//applying to all strains
	list<CStrain*>::iterator it=strains.begin();

	while(it!=strains.end()) {

		std::tr1::poisson_distribution<double> poisson((*it)->N);
		double rnd = poisson(eng);

		if(rnd<1) {
			(*it)->die();
		 	 it=strains.erase(it);
		}
		else
			(*it)->N=rnd;
	++it;
	}

}

void Mutations(){
	list<CStrain*>::iterator it;
	for(it=strains.begin(); it!=strains.end(); it++){
		if(unif(eng)<=mut_rate)
			Mutate(*it);
	}
}

void Update_Immunes(){
	double sumN[rmax+1];
	list<CStrain*>::iterator it;
	int ii=0;
	for(it=strains.begin(); it!=strains.end(); it++){
		ii++;
		for(int i=0; i<rmax+1; i++) sumN[i]=0; 
		(*it)->get_infected(sumN);
		for(int i=0; i<rmax+1; i++) (*it)->M[i]+=nu*sumN[i]*dt; 
	}
	cerr<< "Alive: "<< ii <<endl;
	
	//for(int i=0; i<=rmax; i++){
		//s.M[0]+=nu*s.N*dt;
		//cerr<< s.N <<endl;
	//}

}

void Update(){
	Immune_Selection();
	Genetic_Drift();
	Mutations();
	Update_Immunes();
}

void Run(){
	double sumAllI;

	for(t=dt; t<=tMax; t+=dt){
		Update();
		if(strains.size()==0) break;

		cout << t <<"    "<< stotal <<"    "<< strains.size() <<"    ";

		list<CStrain*>::iterator it;

		sumAllI=0.;

		for(it=strains.begin(); it!=strains.end(); it++){
			sumAllI+=(*it)->N;
			//cout << sumAllI <<"    "<<"    "<< (*it)->N;
		}
		cout << sumAllI << endl;
	}
}

void finish(){

	delete top;
}

int main(){
	Initial_Conditions();
	Run();
return 0;
}


//cout << f0 <<"    " << beta0 <<"    "<< mut_rate <<"     ";