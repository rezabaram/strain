#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <tr1/random>
#include "strain.h"
#include "tools.h"
#include"parameters.h"
#include <ctime>
#include <math.h>
std::tr1::ranlux64_base_01 eng;
std::tr1::uniform_real<double> unif(0, 1);

using namespace std;

//To to be used for assigning ID's
int stotal=0;
//List of all alive strains
list<CStrain*> strains;
CStrain *top=NULL;
//time
double t;
unsigned int iTime=0;
//Cross immunity matrix
double chi[rmax+1];
int Nfiles=1021;

void define_cross_im(){

/*
	double a=1;
	for(int i=0; i<=rmax; i++){
		chi[i]=a;
		a-=1.0/(rmax+1);
		//cout << chi[i] << "    ";
	}
	//cout << endl;
*/


	chi[0]=1.;
	chi[1]=0.99;
	chi[2]=0.95;
	chi[3]=0.9;
	chi[4]=0.85;
	chi[5]=0.8;
	chi[5]=0.7;
	chi[6]=0.65;
	chi[7]=0.5;
	chi[8]=0.2;
	chi[9]=0.2;
	chi[10]=0.2;
	chi[11]=0.2;
	chi[12]=0.2;
	chi[13]=0.2;
	chi[14]=0.2;
	chi[15]=0.2;

}


void Initial_Conditions(){
	eng.seed(time(0));
	//eng.seed(1);
	stotal=0;
	//creating the root node
	top=new CStrain(stotal,NULL);
	stotal++;
	top->N=N0;
	strains.push_back(top);
	define_cross_im();
}

//create a new strain through mutation
//and add to the list
void Mutate(CStrain *pfather){
	CStrain *ps = new CStrain(stotal,pfather);
	ps->N=1;
	pfather->N--;
	if(pfather->N<1) pfather->die();

	strains.push_back(ps);

	stotal++;
}

void Immune_Selection(){
	list<CStrain*>::iterator it;
	//applying to all strains
	for(it=strains.begin(); it!=strains.end(); it++){
		CStrain *s=(*it);
		s->fitness=f0*(1-beta0*s->WeightedSumM(chi) );
		s->N=s->N*(1+s->fitness*dt);//make sure about the order of update N and M
	}
}

//removes an element from the list and returns 
//the iterator to the PREVIOUS element
list<CStrain*>::iterator remove(list<CStrain*>::iterator &it){
		list<CStrain*>::iterator it0;
		//make sure we dont jump over a strain while erasing 
		(*it)->die();
		it0=it;
		it--;
		strains.erase(it0);
		return it;
}

void Genetic_Drift(){
	//applying to all strains
	list<CStrain*>::iterator it=strains.begin();

	while(it!=strains.end()) {

		std::tr1::poisson_distribution<double> poisson((*it)->N);
		double rnd = poisson(eng);

		if(rnd<1) {
			it=remove(it);
		}
		else{
			(*it)->N=rnd;
		}
	it++;
	}

}

void Mutations(){
	list<CStrain*>::iterator it;
	for(it=strains.begin(); it!=strains.end(); it++){
		if(unif(eng)<=mut_rate){
			Mutate(*it);
			if((*it)->N<1) remove(it);
			}
	}
}

void Update_Immunes(){
	double sumN[rmax+1];
	list<CStrain*>::iterator it;
	for(it=strains.begin(); it!=strains.end(); it++){
		for(int i=0; i<rmax+1; i++) sumN[i]=0; 

		(*it)->get_infected(sumN);

		for(int i=0; i<rmax+1; i++) (*it)->M[i]+=nu*sumN[i]*dt; 
	}

}

void Update(){
	Immune_Selection();
	if(iTime%inf_period==0) Genetic_Drift();
	Mutations();
	Update_Immunes();
	//trims the dead leaves
	//if(iTime%100==0) top->trim();
}

void Run(){
	double sumAllI;
	int s=0;
	ofstream singleouts[Nfiles];	
	for(int i=0; i<Nfiles; i++){
		string name ="single"+stringify(i,5,'0');
		singleouts[i].open(name.c_str());
		}

	unsigned int iTimeMax=tMax/dt;
	for(iTime=1; iTime<=iTimeMax; iTime++){
		t=iTime*dt;
		Update();
		if(strains.size()==0) break;
	
		//if(strains.size()>1000 and iTime%(2*inf_period)==0) {mut_rate-=0.0001; s=1;}

		//if(strains.size()<1000 and iTime%(2*inf_period)==0 and s==1) {mut_rate+=0.0001;}
		//if(strains.size()<500 and iTime%inf_period==0) mut_rate+=0.0001;
	
		cout << t <<"    "<< CStrain::stotal <<"    "<< strains.size() <<"    "<< mut_rate <<"    ";

		sumAllI=0.;

		list<CStrain*>::iterator it;
		for(it=strains.begin(); it!=strains.end(); it++){
			sumAllI+=(*it)->N;
			//cout << sumAllI <<"    "<<"    "<< (*it)->N;
			if((*it)->ID<Nfiles and (*it)->ID>=0){
				singleouts[(*it)->ID]<< t<<"  "<<(*it)->N <<endl;
			}
		}
		cout << CStrain::max_dist<<"  ";
		cout << sumAllI << endl;
		if(iTime==500){
			ofstream out("tree");
			top->print(out);
			out.close();
		}
	}
	
}

void finish(){
	
	delete top;
}

int main(){
	Initial_Conditions();
	Run();
	finish();
return 0;
}


//cout << f0 <<"    " << beta0 <<"    "<< mut_rate <<"     ";
