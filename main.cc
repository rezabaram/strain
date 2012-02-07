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
int seed=time(0);
std::tr1::ranlux64_base_01 eng;
std::tr1::uniform_real<double> unif(0, 1);

ofstream logtime("logtime");
CTimer timer;

CConfig config("config");


using namespace std;

//To to be used for assigning ID's
int stotal=0;
//List of all alive strains
vector<CStrain*> allstrains;
list<CStrain*> strains;
CStrain *top=NULL;
//time
double t;
unsigned int iTime=0;
//Cross immunity matrix
double *chi;
const unsigned int Nfiles=1;//1000;
const int IDstart=1;//62300;
double mtotal=0.;
double mutants=0.;

double chi_at_d(double d){
	double chi;
	chi=4./((double)d+4.);
	//chi=1.;

	return chi;
}

void define_cross_im(){

	chi=new double[rmax+1];


// constant function

/*
	double a=0.5;

	for(size_t d=0; d<=rmax; d++){
		chi[d]=a;
	}
*/

// step function
/*
	double a=1.;
	size_t cutoff=3;

	for(size_t d=0; d<=cutoff; d++){
		chi[d]=a;
		//cout << chi[d] << "    ";
	}

	for(size_t d=cutoff+1; d<=rmax; d++){
		chi[d]=0.;
		//cout << chi[d] << "    ";
	}
*/
// double step

/*
	double a=1.;
	double b=0.5;
	size_t cutoff=3;

	for(size_t d=0; d<=cutoff; d++){
		chi[d]=a;
		//cout << chi[d] << "    " << endl;
	}

	for(size_t d=cutoff+1; d<=2*cutoff+1; d++){
		chi[d]=b;
		//cout << chi[d] << "    " << endl;
	}

	for(size_t d=2*cutoff+1; d<=rmax; d++){
		chi[d]=0.;
		//cout << chi[d] << "    " << endl;
	}
*/


// generalized logistic function

/*
	double A=2./10.; //lower asymptote
	double K=1.; //upper asymptote
	double B=9./2.;
	double Q=5./3.;
	double d0=4.;

	for(size_t d=0; d<=rmax; d++){
		chi[d] = A + (K-A)/(1.+Q*exp(B*(d-d0)));
		//cout << chi[d] << "    " << endl;
	}
*/

// hyperbola

/*
	double m=4.;
	double y0=0.; // asymptote
	double x0=-4.;

	for(size_t d=0; d<=rmax; d++){
		chi[d] = m/(d-x0) + y0;
		//cout << chi[d] << "    " << endl;
	}
*/

// inverse

/*
	chi[0]=1;

	for(size_t d=1; d<=rmax; d++){
		chi[d] = 1./(double)d;
		//cout << chi[d] << "    " << endl;
	}
*/

// shifted inverse
	
	for(size_t d=0; d<=rmax; d++){
		chi[d] = 1./((double)d+1.);
		//cout << chi[d] << "    " << endl;
	}

}


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
double Nall(){
	double sumAllN=0.;

	list<CStrain*>::iterator it;
	for(it=strains.begin(); it!=strains.end(); it++){
		//cerr << (*it)->N << endl;
		
		assert((*it)->N>0);
		sumAllN+=(*it)->N;	
	}
	return sumAllN;	
}

void Immune_Selection(){
	list<CStrain*>::iterator it;
	//applying to all strains
	for(it=strains.begin(); it!=strains.end(); it++){
		CStrain *s=(*it);
		s->M0+=s->WeightedSumM(chi_at_d)*nu*dt;
		s->fitness=f0*(1-beta0*s->M0*Nall()*cp) - s->red_m*Cf;
		s->N=s->N*(1+s->fitness*dt);//make sure about the order of update N and M
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

	double dist=0.;
	bool red=true;
	if (unif(eng) <= 1.) {dist=1.; red=false;}

	CStrain *ps = new CStrain(stotal,pfather,dist);
	if(red) ps->red_m++;

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
}


void output(ostream &out){

	vector<CStrain*>::iterator it; 

	double totalN=top->calSubN();

	for(it=allstrains.begin(); it!=allstrains.end(); it++){
		(*it)->setFreq((*it)->SubN/totalN);
	}			

	out << t <<"    "<< CStrain::stotal <<"    "<< strains.size() <<"    "<< mut_rate <<"    ";
	out << CStrain::max_dist<<"    ";
	//out << Diversity() <<"   ";
	out << Nall() <<"    ";
	//out << mtotal << "    ";
	out << totalN;
	out << endl;

}

void FreqDist(){
	double bins[20];
	int max_i;
	vector<CStrain*>::iterator it;
	for(it=allstrains.begin(); it!=allstrains.end(); it++){
		cout<<(*it)->SubN<<"    "<<endl;
		max_i=19*(*it)->maxFreq;	
		for(int i=0; i<=max_i;i++){
			bins[i]++;
		}
	}
	
	for(int i=0; i<20;i++){
		cout<<i/20.0<<"   "<<bins[i]<<endl;
	}

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
		t=iTime*dt;
		Update();
		if(strains.size()==0) {
			output_graphic_tree();
			break;
		}
	
		PrintSingleInfected();
		//if(iTime%200==0)

		output(out);

		//if(iTime%(inf_period*7)==0 and t >= 20.) 
		//print_diversity(outdiv);

		//if(iTime%(inf_period*7)==0) 
		print_fitness(outfit);

		print_N(outN);

		//if(iTime%200==0 and t <= 4.) output_graphic_tree();
	}
}

void finish(){
	FreqDist();
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
