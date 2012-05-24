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

#include<fstream>
bool secondrun=false;
vector<unsigned int> TrackedIDs;
ofstream * trackedouts;
unsigned int NTrackedFiles=0;

int seed=time(0);
std::tr1::ranlux64_base_01 eng;
std::tr1::uniform_real<double> unif(0, 1);

ofstream logtime("logtime");
ofstream outfixtimes("FixTimes");
ofstream outfixtimesgreen("FixTimesGreen");
ofstream outselcoef("SelectionCoefficient");
ofstream outfixperyear("FixPerYear");
ofstream outfixdist("FixDist");
ofstream outmeanfixdist("MeanFixDist");
ofstream outlifetime("LifeTime");
ofstream outfractionspropagators("FractionsPropagators");
ofstream outpropagators("Propagators");
ofstream outIDsfixed("IDsFixed");
CTimer timer;

//CConfig config("config");


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
const unsigned int Nfiles=1000;
const int IDstart=40000;
double mtotal=0.;
double mutants=0.;

double A=2./10.; //lower asymptote
double K=1.; //upper asymptote
double B=3./2.;
double Q=5./3.;
double d0=4.;

double chi_at_d(double d){
	double chi;
	//chi=4./((double)d+4.);
	//chi=1.;

	chi = A + (K-A)/(1.+Q*exp(B*(d-d0)));
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

	if(secondrun){
		ifstream id_file("IDsFixedInput");
		unsigned int temp_id;
		while(!id_file.eof()){ //Reading the IDs to be tracked
			id_file>>temp_id;
			if(id_file.eof())break;
			TrackedIDs.push_back(temp_id);
			}
		
		NTrackedFiles=TrackedIDs.size();
		trackedouts=new ofstream [NTrackedFiles];
		for(int i=0; i<TrackedIDs.size(); i++){
		        string name ="tracked"+stringify(i,5,'0');
                	trackedouts[i].open(name.c_str());

			}
		}

	//seed=1321702805;
	t=0.;
	cerr<< "Seed: "<<seed <<endl;
	eng.seed(seed);
	//eng.seed(10);
	//rmax=10;
	//rmax=config.get_param<size_t>("rmax");
	stotal=0;
	//creating the root node
	top=new CStrain(stotal,NULL);
	stotal++;
	top->N=N0;
	top->M0=0.;
	top->crtime=t;
	top->mut_type=1;
	top->cost=0.;
	strains.push_back(top);
	allstrains.push_back(top);
	//define_cross_im();
}

//create a new strain through mutation
//and add to the list
double Nall(){
	double Ntot=0.;

	list<CStrain*>::iterator it;
	for(it=strains.begin(); it!=strains.end(); it++){
		//cerr << (*it)->N << endl;
		assert((*it)->N>0);
		Ntot+=(*it)->N;	
	}
	return Ntot;	
}

double average_red_m(){
	double av_red_m=0.;
	double Ntot=Nall();

	list<CStrain*>::iterator it;
	for(it=strains.begin(); it!=strains.end(); it++){
		av_red_m+=(*it)->red_m*(*it)->N/Ntot;
	}

	return av_red_m;	
}

double average_cost(){
	double av_cost=0.;
	double Ntot=Nall();	

	list<CStrain*>::iterator it;
	for(it=strains.begin(); it!=strains.end(); it++){
		av_cost+=(*it)->cost*(*it)->N/Ntot;
	}
	return av_cost;	
}


void Immune_Selection(){
	list<CStrain*>::iterator it;
	//applying to all strains

	double Ntot=Nall();

	for(it=strains.begin(); it!=strains.end(); it++){
		CStrain *s=(*it); 
		s->M0+=s->WeightedSumM(chi_at_d)*nu*dt;
		//s->fitness=f0*(1-beta0*s->M0*Nall()*cp) - s->red_m*Cf;
		s->fitness=f0*(1-beta0*s->M0*Ntot*cp) - s->cost;
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
			//cout<<(*it)->fitness<<endl;
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

//ofstream outdist("dist.dat");

void Mutate(CStrain *pfather){

	double dist;
	int ii=-1;

	double prob=unif(eng);

	if(prob<=(double)Lep/(double)L){ ii=1; dist=1.; } //green
	else if(prob<=((double)Lep+(double)Lnep)/(double)L){ ii=0; dist=0.; } //red
	else { ii=2; dist=0.; } //blue

	CStrain *ps = new CStrain(stotal,pfather,dist);

	if(ii==0) { 
		//ps->cost+=Cf;
		ps->mut_type = 0;
		std::tr1::exponential_distribution<double> exponential(1./Cf);
		double rnd = exponential(eng);

		if(unif(eng)<((double)Lnep-ps->red_m)/(double)Lnep) {
			ps->red_m++;
			ps->cost+=rnd;
		}
		else {//if(ps->red_m>=1 && ps->cost>0.){
			ps->red_m=ps->red_m-1;
			if(ps->cost>=rnd) ps->cost=ps->cost-rnd;
			else ps->cost=0.;
			//}
		}
	} //red
	if(ii==1) {ps->mut_type = 1;} //green
	if(ii==2) {ps->mut_type = 2;} //blue
	if(ii<0) {cerr<<"Error in Mutate function"<<endl;}

	ps->M0=ps->WeightedSumM0(chi_at_d);
	ps->crtime=t;

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
				//cout<<(*it)->fitness<<endl;
                		(*it)->die();
                		//this also sets "it" to next value
                		it=strains.erase(it);
                		continue;
            		}
        	}
        	it++;
    	}
}

double infperyear=0.;
//int iin=0;

void Update_Immunes(){
	list<CStrain*>::iterator it;
	for(it=strains.begin(); it!=strains.end(); it++){
		(*it)->accN+=nu*dt*(*it)->N;
		infperyear+=nu*dt*(*it)->N;
	}

	//if (t>=40. && t<=41.) {iin++;}
}

double Diversity(){
	list<CStrain*>::iterator it;
	double diversity=0.;
	double Ntot=Nall();
	for(it=strains.begin(); it!=strains.end(); it++){
		double div=0;
		(*it)->get_diversity2(div,0,(*it));
		//(*it)->get_diversity(div);
		diversity+=(*it)->N*div;
	}
	diversity=diversity/2./(Ntot*Ntot);
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
		(*it)->setFreq((*it)->SubN/totalN,t);
	}

	out << t << "    " << CStrain::stotal << "    " << strains.size() << "    " << average_red_m() << "    " << average_cost() << "    " << CStrain::max_dist << "    " << Nall() << endl;
}


void FreqDist(){
	double binsGreen[nbins], binsRed[nbins], binsBlue[nbins];
	
	int fixperyear[(int)tMax+1];
	int fixperyearred[(int)tMax+1];
	int fixperyeargreen[(int)tMax+1];
	int fixperyearblue[(int)tMax+1];
	int fixdist[maxmut];
	int fixdistred[maxmut];
	int fixdistgreen[maxmut];
	int fixdistblue[maxmut];
	
	for(int i=0; i<=(int)tMax; i++){
		fixperyear[i]=0;
		fixperyearred[i]=0;
		fixperyeargreen[i]=0;
		fixperyearblue[i]=0;
	}

	for(int i=0; i<maxmut; i++){
		fixdist[i]=0;
		fixdistred[i]=0;
		fixdistgreen[i]=0;
		fixdistblue[i]=0;
	}

	for(int i=0; i<=(nbins-1); i++){
		binsGreen[i]=0;
		binsRed[i]=0;
		binsBlue[i]=0;
	}

	int max_i;
	vector<CStrain*>::iterator it;
	int ynum;
	double timetofix=0.;

	for(it=allstrains.begin(); it!=allstrains.end(); it++){
		//cout<<(*it)->fitness<<endl;
		//cout<<(*it)->maxFreq<<"    ";

		//if( (double)timestart<=(*it)->crtime && (*it)->crtime<=(double)timeend /*&& (*it)->maxFreq > 0.1*/ ){
			//cerr<<(*it)->red_m<<"   ";
			max_i=(nbins-1)*(*it)->maxFreq;

			if(max_i==(nbins-1)){
				outfixtimes << (*it)->mut_type << "    " << (*it)->fixtime << "    " << (*it)->fixtime-(*it)->crtime << "    " << (*it)->cost << endl;
				
				timetofix+=(*it)->fixtime-(*it)->crtime;
				ynum=(*it)->fixtime;
				fixperyear[ynum]++;
				if( (*it)->mut_type==0 ){ fixperyearred[ynum]++; }
				else if ( (*it)->mut_type==1 ){ fixperyeargreen[ynum]++;
					outfixtimesgreen << (*it)->ID << "    " << (*it)->crtime << "    " << (*it)->fixtime << "    " << (*it)->fixtime-(*it)->crtime << "    " << (*it)->cost << endl;
					outIDsfixed << (*it)->ID << endl;}
				else if ( (*it)->mut_type==2 ){ fixperyearblue[ynum]++; }
			}

			for(int i=0; i<=max_i; i++){
				if ((*it)->mut_type==1){binsGreen[i]++;}
				else if ((*it)->mut_type==0){binsRed[i]++;}
				else if ((*it)->mut_type==2){binsBlue[i]++;}
				//else cerr<<"Error in the FreqDist function"<<endl;
				//else cout << "ERROR"<<endl;
			}
		//}
	}
	cerr << endl;

	int mean=0;

	for(int i=0; i<=(int)tMax; i++){
		//cerr<<fixperyear[i]<<"   "<<endl;
		fixdist[fixperyear[i]]++;
		fixdistred[fixperyearred[i]]++;
		fixdistgreen[fixperyeargreen[i]]++;
		fixdistblue[fixperyearblue[i]]++;
	}

	for(int i=0; i<=(int)tMax; i++){
		outfixperyear << i << "    " << fixperyear[i] << "    " << fixperyearred[i] << "    " << fixperyeargreen[i] << "    " << fixperyearblue[i] << "    " << endl;
	}

	for(int i=0; i<maxmut; i++){
		mean+=i*fixdist[i];
		//if(fixdist[i]!=0){
		outfixdist << i << "    " << fixdist[i] << "    " << fixdistred[i] << "    " << fixdistgreen[i] << "    " << fixdistblue[i] << "    " << endl;
		//}
	}	

	outmeanfixdist << (double)mean/tMax << endl;
	outlifetime << timetofix/(binsBlue[nbins-1]+binsGreen[nbins-1]+binsRed[nbins-1]) << endl;
	
	for(int i=0; i<=(nbins-1); i++){
		//cout<<binsRed[i]<<endl;
		outpropagators<<i/(nbins-1.)<<"   "<<log((double)binsGreen[i]/(double)binsGreen[0])/log(10.)<<"   "<<(double)binsGreen[i]/(double)binsGreen[0]<<"    "<<binsGreen[i]<<"    "<<log((double)binsRed[i]/(double)binsRed[0])/log(10.)<<"    "<<(double)binsRed[i]/(double)binsRed[0]<<"    "<<binsRed[i]<<"   "<<log((double)binsBlue[i]/(double)binsBlue[0])/log(10.)<<"   "<<(double)binsBlue[i]/(double)binsBlue[0]<<"    "<<binsBlue[i]<<endl;
	}
	//cout<<endl;

	for(int i=0; i<=(nbins-1); i++){
		outfractionspropagators<<i/(nbins-1.)<<"   "<< log((double)binsGreen[i]/(double)binsGreen[0]/((double)binsBlue[i]/(double)binsBlue[0]))/log(10.)<<"   "<<log((double)binsRed[i]/(double)binsRed[0]/((double)binsBlue[i]/(double)binsBlue[0]))/log(10.)<<endl;
	}
	//cout<<endl;

	for(int i=0; i<=(nbins-1); i++){
		outfractionspropagators<<i/(nbins-1.)<<"   "<< (double)binsGreen[i]/(double)binsGreen[0]/((double)binsBlue[i]/(double)binsBlue[0])<<"   "<<(double)binsRed[i]/(double)binsRed[0]/((double)binsBlue[i]/(double)binsBlue[0])<<endl;
	}


}

void print_diversity(ostream &out){
	out << t << "    " << Diversity() << "    " << endl;
}

void print_fitness(ostream &out){

	out << t << "    ";

	list<CStrain*>::iterator it;
	double Ntot=Nall();
	double mean_fitness=0.;

	for(it=strains.begin(); it!=strains.end(); it++){
		mean_fitness+=(*it)->fitness*(*it)->N/Ntot;
	}

	out << mean_fitness << "    ";

	for(it=strains.begin(); it!=strains.end(); it++){
		out << (*it)->fitness << "    ";
	}

	out << endl;
}

void print_N(ostream &out){

	out << t <<"    ";

	list<CStrain*>::iterator it;

	for(it=strains.begin(); it!=strains.end(); it++){
		out << (*it)->N << "    ";
	}
	out << endl;
}

void print_recovered(ostream &out){

	out << t << "    " << infperyear << endl;
	infperyear=0.;
}

//checkes whether an id is in the tracked vector
//if so returns the its location in the vector
//otherwise -1
int tracked_index(unsigned int id){
	
	for(int i=0; i<NTrackedFiles; i++){
	if(TrackedIDs.at(i)==id) return i;
	}
	return -1;
}

void PrintTracked(){

	list<CStrain*>::iterator it;
	double Ntot=Nall();
	double mean_fitness=0.;

	for(it=strains.begin(); it!=strains.end(); it++){
		mean_fitness+=(*it)->fitness*(*it)->N/Ntot;
	}

	for(it=strains.begin(); it!=strains.end(); it++){
		int ind=tracked_index((*it)->ID);
		if(ind>=0){

			double weighsumf=0., subtrN=0.;
			(*it)->calSubMeanFitness(weighsumf, subtrN);

			double mean_fitness_subtree=weighsumf/subtrN;
			double mean_fitness_rest=(mean_fitness-mean_fitness_subtree*subtrN/Ntot)*Ntot/(Ntot-subtrN);

			trackedouts[ind] << t << "    " << strains.size() << "    " << subtrN << "    " << (*it)->N << "    " << (*it)->fitness << "    " << mean_fitness << "    " << mean_fitness_subtree << "    " << mean_fitness_rest << endl;
		}
	}
}

ofstream singleouts[Nfiles];

void PrintSingleInfected(){

	list<CStrain*>::iterator it;
	double Ntot=Nall();
	double mean_fitness=0.;

	for(it=strains.begin(); it!=strains.end(); it++){
		mean_fitness+=(*it)->fitness*(*it)->N/Ntot;
	}

	for(it=strains.begin(); it!=strains.end(); it++){
		if((*it)->ID<(int)Nfiles+IDstart and (*it)->ID>=IDstart){

			double weighsumf=0., subtrN=0.;
			(*it)->calSubMeanFitness(weighsumf, subtrN);

			double mean_fitness_subtree=weighsumf/subtrN;
			double mean_fitness_rest=(mean_fitness-mean_fitness_subtree*subtrN/Ntot)*Ntot/(Ntot-subtrN);

			singleouts[(*it)->ID-IDstart] << t << "    " << strains.size() << "    " << subtrN << "    " << (*it)->N << "    " << (*it)->fitness << "    " << mean_fitness << "    " << mean_fitness_subtree << "    " << mean_fitness_rest << endl;
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
	ofstream outrecovered("recovered");
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

		if( iTime%365==0 ) print_recovered(outrecovered);

		
		if(strains.size()==0) {
			//output_graphic_tree();
			break;
		}
		
		//if(iTime%(inf_period*7)==0 and t >= 20.)

		//PrintSingleInfected();
		PrintTracked();
		output(out);
		//print_diversity(outdiv);
		//print_fitness(outfit);
		//print_N(outN);

		//if(iTime%200==0 and t <= 4.) output_graphic_tree();
	}
}

void finish(){
	//cerr<<infperyear<<"   "/*<<iin*/<<endl;
	FreqDist();
	delete top;
}

int main(int argc, char **argv){

	if(argc>1)seed=atoi(argv[1]);
	if(argc>2)secondrun=atoi(argv[2]);
	signal(30,SignalControlReport);
	signal(31,SignalControlReport);
	signal(32,SignalControlReport);

	Initial_Conditions();

	Run();
	finish();

return 0;

}


//cout << f0 <<"    " << beta0 <<"    "<< mut_rate <<"     ";
