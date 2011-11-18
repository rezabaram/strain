#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <ctime>
#include <math.h>
#include "timer.h"
#include <tr1/random>
#include "strain.h"
#include "tools.h"
#include"parameters.h"
#include"io.h"
#include"signal.h"
std::tr1::ranlux64_base_01 eng;
std::tr1::uniform_real<double> unif(0, 1);

ofstream logtime("logtime");
CTimer timer;

using namespace std;

//To to be used for assigning ID's
int stotal=0;
//List of all alive strains
vector<CStrain *> allstrains;
list<CStrain*> strains;
CStrain *top=NULL;
//time
double t;
unsigned int iTime=0;
//Cross immunity matrix
double *chi;
const unsigned int Nfiles=1000;

void define_cross_im(){

	chi=new double[rmax+1];
/*
	double a=1;
	for(int i=0; i<=rmax; i++){
		chi[i]=a;
		a-=1.0/(rmax+1);
		//cout << chi[i] << "    ";
	}
	//cout << endl;
*/

// step function with rmax=10;
/*
	double a=1.;

	for(int d=0; d<=rmax; d++){
		chi[d]=a;
	}	
*/

/*
	double A=0.; //lower asymptote
	double K=1.; //upper asymptote
	double B=2.5; 
	double Q=1.;
	double d0=10.;
*/

/*
	double A=0.; //lower asymptote
	double K=1.; //upper asymptote
	double B=0.5; 
	double Q=1.;
	double d0=10.;
*/


	double A=0.3; //lower asymptote
	double K=1.; //upper asymptote
	double B=0.5; 
	double Q=1.;
	double d0=10.;

	for(int d=0; d<=rmax; d++){
		chi[d] = A + (K-A)/(1.+Q*exp(B*(d-d0)));
		//cout << chi[d] << "    ";	
	} 
	//cout << endl;

}


void Initial_Conditions(){
	int seed=time(0);
	cerr<< "Seed: "<<seed <<endl;
	eng.seed(seed);
	//eng.seed(10);
	stotal=0;
	//creating the root node
	top=new CStrain(stotal,NULL);
	stotal++;
	top->N=N0;
	top->M[0]=100000.;
	strains.push_back(top);
	allstrains.push_back(top);
	define_cross_im();
}

//create a new strain through mutation
//and add to the list
void Mutate(CStrain *pfather){
	CStrain *ps = new CStrain(stotal,pfather);
	pfather->N--;
	if(pfather->N<1) pfather->die();

	strains.push_back(ps);
	allstrains.push_back(ps);

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

double Diversity(){
	list<CStrain*>::iterator it;
	double diversity=0;
	for(it=strains.begin(); it!=strains.end(); it++){
		double div=0;
		(*it)->get_diversity(div,0,(*it));
		//(*it)->get_diversity(div);
		diversity+=(*it)->N*div/Npop;
	}
	diversity=diversity/2/Npop;
	return diversity;
}

void Update(){
	Immune_Selection();
	if(iTime%inf_period==0) Genetic_Drift();
	Mutations();
	Update_Immunes();
	//trims the dead leaves
	if(iTime%100==0) top->trim();
}


void output(ostream &out){
	double sumAllI=0.;

	list<CStrain*>::iterator it;
	for(it=strains.begin(); it!=strains.end(); it++){
		sumAllI+=(*it)->N;
	}
	out << t <<"    "<< CStrain::stotal <<"    "<< strains.size() <<"    "<< mut_rate <<"    ";
	out << CStrain::max_dist<<"  ";
	//out << Diversity() <<"  ";
	out << sumAllI <<"   ";
	out<< endl;
}

ofstream singleouts[Nfiles];	

void PrintSingleInfected(){
	list<CStrain*>::iterator it;
	for(it=strains.begin(); it!=strains.end(); it++){
		if((*it)->ID<Nfiles and (*it)->ID>=0){
			singleouts[(*it)->ID]<< t<<"  "<<(*it)->N <<endl;
		}
	}
}

void Run(){
	ofstream out("out");
	int s=0;
	for(int i=0; i<Nfiles; i++){
		string name ="single"+stringify(i,5,'0');
		singleouts[i].open(name.c_str());
	}

	unsigned int iTimeMax=tMax/dt;
	for(iTime=1; iTime<=iTimeMax; iTime++){

		if(bSignal==0) {
        		mut_rate*=0.95;
			//mut_rate-=0.0001;
                  	bSignal=-1;
                       	cerr<< "New Mutation rate:" << mut_rate <<"   "<< "Time:" << iTime*dt << endl;
                }
               	if(bSignal==1) {
			mut_rate*=1.15;
           		//mut_rate+=0.0001;
                	bSignal=-1;
                       	cerr<< "New Mutation rate:" << mut_rate <<"   "<< "Time:" << iTime*dt << endl;
                }


		//logtime<<timer.read()/strains.size()<<"   "<<t<<endl;
		logtime<<timer.read()<<"   "<<t<<endl;
		t=iTime*dt;
		Update();
		if(strains.size()==0) break;
	
		PrintSingleInfected();
		//if(iTime%200==0) 
		output(out);
		if(iTime==2500){
			//prints two versions of the tree
			//ofstream tree1("tree1");
			//top->print(tree1);
			//tree1.close();

			//ofstream tree("tree");
			//SaveState(tree, allstrains);
			//tree.close();
			
			cerr<< "tree printed " <<endl;
			ofstream tree("tree");
			top->print(tree, 0, 0.5);
			tree.close();
			//testing if the file was read correctly
			//vector <CStrain*> temp;
			//ifstream treein("tree");
			//ReadState(treein, temp);
			//treein.close();

			//ofstream tree2("tree2");
			//SaveState(tree2, temp);
			//tree2.close();
		}
	}
	
}

void finish(){
	
	delete top;
}

int main(){

	signal(30,SignalControlReport);
	signal(31,SignalControlReport);
	signal(32,SignalControlReport);

	Initial_Conditions();

	Run();
	finish();

return 0;

}


//cout << f0 <<"    " << beta0 <<"    "<< mut_rate <<"     ";
