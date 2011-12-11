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
#include"config.h"
#include"chi.h"
int seed=time(0);
std::tr1::ranlux64_base_01 eng;
std::tr1::uniform_real<double> unif(0, 1);

ofstream logtime("logtime");
ofstream histout("pair_distances");
bool print_hist=false;
bool print_hist_on=true;

CTimer timer;

CConfig config("config");


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
const unsigned int Nfiles=1;//1000;
const int IDstart=1;//62300;
double mtotal=0.;
double mutants=0.;

void Initial_Conditions(){
	//seed=1321702805;
	cerr<< "Seed: "<<seed <<endl;
	eng.seed(seed);
	//eng.seed(10);
	rmax=config.get_param<size_t>("rmax");
	stotal=0;
	//creating the root node
	top=new CStrain(stotal,NULL);
	stotal++;
	top->N=N0;
	top->M0=0.;
	strains.push_back(top);
	allstrains.push_back(top);
	//define_cross_im();
}

//create a new strain through mutation
//and add to the list

void Immune_Selection(){
	list<CStrain*>::iterator it;
	//applying to all strains
	for(it=strains.begin(); it!=strains.end(); it++){
		CStrain *s=(*it);
		s->M0+=s->WeightedSumM(chi_at_d)*nu*dt;
		s->fitness=f0*(1-beta0*s->M0);
		s->N=s->N*(1+s->fitness*dt);//make sure about the order of update N and M
		assert(s->N>0);
	}
}


void Genetic_Drift(){
	//applying to all strains
	list<CStrain*>::iterator it=strains.begin();

	while(it!=strains.end()) {

		std::tr1::poisson_distribution<double> poisson((*it)->N);
		double rnd = poisson(eng);

		//cerr << "random number" << "    " << rnd << "    " << endl;		

		if(rnd<1) {
			//cerr << "random number smaller than 1" << "    " << rnd << "    " << endl;
			(*it)->die();
			//this also sets "it" to next value
			it=strains.erase(it);
		}
		else{
			(*it)->N=rnd;
			//cerr << (*it)->N << endl;
			it++;
		}
	}

	
	//for(it=strains.begin(); it!=strains.end(); it++){
	//	cerr << (*it)->N << endl;
	//}

}

void Mutate(CStrain *pfather){
	CStrain *ps = new CStrain(stotal,pfather);

	ps->M0=ps->WeightedSumM0(chi_at_d);

	/*
	for(size_t i=1;i<=rmax;i++){
		ps->M[i]=pfather->M[i-1];
	}
	*/

	pfather->N--;

	strains.push_back(ps);
	allstrains.push_back(ps);

	stotal++;
}

void Mutations(){
	list<CStrain*>::iterator it=strains.begin();
	mtotal=0.;
	while(it!=strains.end()) {
		int nn=(*it)->N;
		for(int i=0; i<nn; i++){
			if(unif(eng)<=mut_rate){
                		Mutate(*it);
				mtotal++;
        		}
        		if((*it)->N<1) {
                		(*it)->die();
                		//this also sets "it" to next value
                		it=strains.erase(it);
                		continue;
            		}
        	}
        	it++;
    	}
}


void Mutations2(){
	list<CStrain*>::iterator it=strains.begin();
	mtotal=0.;
	while(it!=strains.end()) {
		double nn=(*it)->N;
		std::tr1::poisson_distribution<double> poisson( mut_rate*nn );
		double rnd = poisson(eng);
		int num_mutants = min(rnd,nn);

		mtotal+=num_mutants;
		//cerr << "number of mutants" << "    " << rnd << "    " << endl;
		for(int i=1; i<=num_mutants; i++){
                	Mutate(*it);
        		if((*it)->N<1) {
                		(*it)->die();
                		//this also sets "it" to next value
                		it=strains.erase(it);
                		continue;
            		}
        	}
        	it++;
    	}
}

void Update_Immunes(){
	list<CStrain*>::iterator it;
	for(it=strains.begin(); it!=strains.end(); it++){
		(*it)->accN+=nu*dt*(*it)->N;
	}
}

double Diversity(){
	list<CStrain*>::iterator it;
	double diversity=0.;
	double sumAllN=0.;
	for(it=strains.begin(); it!=strains.end(); it++){
		double div=0;
		(*it)->get_diversity2(div,0,(*it));
		//(*it)->get_diversity(div);
		diversity+=(*it)->N*div;
		sumAllN+=(*it)->N;
	}
	diversity=diversity/2/(sumAllN*sumAllN);
	return diversity;
}

void Update(){
	Immune_Selection();
	if(iTime%inf_period==0) Genetic_Drift();
	Mutations2();
	Update_Immunes();
	//trims the dead leaves
	//if(iTime%10==0) top->trim(); not necessary!
	if(iTime%inf_period==0) top->make_bridges();
	if(print_hist) histout<<endl;
}

void output(ostream &out){

	double sumAllN=0.;

	list<CStrain*>::iterator it;
	for(it=strains.begin(); it!=strains.end(); it++){
		//cerr << (*it)->N << endl;
		
		assert((*it)->N>0);
		sumAllN+=(*it)->N;	
	}
	out << t <<"    "<< CStrain::stotal <<"    "<< strains.size() <<"    "<< mut_rate <<"    ";
	out << CStrain::max_dist<<"    ";
	//out << Diversity() <<"   ";
	out << sumAllN <<"    ";
	out << mtotal << "    ";
	out << endl;
}

void print_diversity(ostream &out){
	out << t <<"    ";
	out << Diversity() <<"    ";
	out << endl;
}

void print_fitness(ostream &out){

	out << t <<"    ";

	list<CStrain*>::iterator it;

	for(it=strains.begin(); it!=strains.end(); it++){
		out << (*it)->fitness <<"    ";
	}
	out << endl;
}

void print_N(ostream &out){

	out << t <<"    ";

	list<CStrain*>::iterator it;

	for(it=strains.begin(); it!=strains.end(); it++){
		out << (*it)->N <<"    ";
	}
	out << endl;
}

ofstream singleouts[Nfiles];

void PrintSingleInfected(){
	list<CStrain*>::iterator it;
	for(it=strains.begin(); it!=strains.end(); it++){
		if((*it)->ID<(int)Nfiles+IDstart and (*it)->ID>=IDstart){
			singleouts[(*it)->ID-IDstart]<< t<< "    " <<(*it)->N << "    " << (*it)->fitness << endl;
		}
	}
}

void output_graphic_tree(){
	static int outcount=0;
	outcount++;
	//prints two versions of the tree
	//ofstream tree1("tree1");
	//top->print(tree1);
	//tree1.close();

	//ofstream tree("tree");
	//SaveState(tree, allstrains);
	//tree.close();
	
	//cerr<< "tree printed " <<endl;
	ofstream tree(("tree"+stringify(outcount,5, '0')).c_str());
//	ofstream tree2("tree2");
	//top->cal_print_spaces();
	COffset offset=top->cal_offsets();
	tree<<offset.width()<<"  "<<offset.bottom<<"   "<<offset.top<<endl;
	top->print(tree, -0, 0.0);
	top->print_bridges(tree);
	//top->print2(tree2, -0, 0.0);
	tree.close();
	//tree2.close();
	//exit(0); to go out
	//testing if the file was read correctly
	//vector <CStrain*> temp;
	//ifstream treein("tree");
	//ReadState(treein, temp);
	//treein.close();

	//ofstream tree2("tree2");
	//SaveState(tree2, temp);
	//tree2.close();
}

void Run(){
	ofstream out("out");
	ofstream outdiv("diversity");
	ofstream outfit("fitness");
	ofstream outN("N");
	int s=0;
	for(size_t i=0; i<Nfiles; i++){
		string name ="single"+stringify(i,5,'0');
		singleouts[i].open(name.c_str());
	}

	unsigned int iTimeMax=tMax/dt;
	for(iTime=1; iTime<=iTimeMax; iTime++){
		t=iTime*dt;

		print_hist=false;
		if(print_hist_on and iTime%(inf_period*90)==0 and t>=15.){
			print_hist=true;
			histout<<t<<"   ";
		}

		CStrain::max_dist=0;
		if(bSignal==0) {
        		//mut_rate*=0.95;
			mut_rate-=0.0001;
                  	bSignal=-1;
                       	cerr<< "New Mutation rate:" << mut_rate <<"   "<< "Time:" << iTime*dt << endl;
                }
               	if(bSignal==1) {
			//mut_rate*=1.15;
           		mut_rate+=0.0001;
                	bSignal=-1;
                       	cerr<< "New Mutation rate:" << mut_rate <<"   "<< "Time:" << iTime*dt << endl;
                }


		//logtime<<timer.read()/strains.size()<<"   "<<t<<endl;
		logtime<<timer.read()<<"   "<<t<<endl;

		Update();
		if(strains.size()==0) {
			output_graphic_tree();
			break;
		}
	
		PrintSingleInfected();
		//if(iTime%200==0)
		output(out);

		//if(iTime%(inf_period*7)==0 and t >= 10.) 
		//print_diversity(outdiv);

		//if(iTime%(inf_period*7)==0) 
		print_fitness(outfit);

		print_N(outN);

		//if(iTime%200==0 and t <= 4.) output_graphic_tree();
	}
}

void finish(){
	
	delete top;
}

int main(int argc, char **argv){

	if(argc>1)seed=atoi(argv[1]);
	signal(30,SignalControlReport);
	signal(31,SignalControlReport);
	signal(32,SignalControlReport);

	Initial_Conditions();

	Run();
	finish();

return 0;

}


//cout << f0 <<"    " << beta0 <<"    "<< mut_rate <<"     ";
