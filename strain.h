#ifndef STRAIN_H
#define STRAIN_H
#include <vector>
#include<assert.h>
#include"parameters.h"

using namespace std;

class CStrain{
	public:
	CStrain(int i, CStrain *f);
	~CStrain();
	void SetAlive();
	double sumM();
	double WeightedSumM(double *chi);
	int count_neigh();
	void die();
	void trim();
	CStrain *father();
	void add_neighbour(CStrain *ps){neighbours.push_back(ps);is_leaf=false;};
	void get_infected(double *sumN, int distance=0, CStrain *exclude=NULL);
	void get_diversity(double &div, int distance=0, CStrain *exclude=NULL);
	void print(ostream &out);
	void print_node(ostream &out)const;
	void print(ostream &out, double x, double y);
	double cal_print_spaces();

	double fitness;
	double N;
	int ID;
	bool dead;
	bool is_leaf;
	double *M;
	static int max_dist;
	std::vector<CStrain*> neighbours;
	static unsigned int stotal;
	private:
	double print_space;
};
unsigned int CStrain::stotal=0;

//this is just used in function get_infected()
//to find out what was the maximum distance
//ever reached between the alive nodes
int CStrain::max_dist=0;

//Constructor take an int for ID and the point of the
//father node
CStrain::CStrain(int i, CStrain *f){
	print_space=0.05;
	stotal++;
	ID=i; 
	N=0.0; 
	fitness=0.0;
	if(f!=NULL){
		f->add_neighbour(this);
	}
	neighbours.push_back(f);

	dead=true;
	is_leaf=true;
	if(ID>=0) {
		SetAlive();
		N=1;
	}
}

void CStrain::SetAlive(){
	dead=false;
	M=new double[rmax+1];
	for(int i=0;i<=rmax;i++){
		M[i]=0.0;
	}
}

//Cleans the allocated memory if not yet cleaned 
CStrain::~CStrain(){
	if( M!=NULL) delete[] M;
	stotal--;
}

//Returns the pointer to the father if not root node
CStrain* CStrain::father(){
	return neighbours.at(0);
}

//Returns the sum of N of all strains at distances up to rmax
//sumN[0] is the N of the strain itself
//sumN[1] is the sum of N's of all strains at distance 1
//....
//sumN[rmax] is the sum of N's of all strains at distance rmax
void CStrain::get_infected(double *sumN, int distance, CStrain *exclude){
	sumN[distance]+=N;
	if(distance>max_dist)max_dist=distance;
	if(distance==rmax) return;

	
	if(father()!=NULL and father()!=exclude) 
		father()->get_infected(sumN, distance+1, this);

	if(is_leaf) return;
	for(int i=1; i<neighbours.size(); i++){
		if(neighbours.at(i)==exclude) continue;
		neighbours.at(i)->get_infected(sumN, distance+1, this);
	}
	return;
}

void CStrain::get_diversity(double &diversity, int distance, CStrain *exclude){
	diversity+=N*distance;

	for(int i=0; i<neighbours.size(); i++){
		if(neighbours.at(i)==NULL) continue;
		if(neighbours.at(i)==exclude) continue;
		neighbours.at(i)->get_diversity(diversity, distance+1, this);
	}
	return;
}


//Returns the number of alive neighbours of the strain
int CStrain::count_neigh(){
	int count=0;
	for (int i=0; i<neighbours.size();i++){
		if(neighbours.at(i)==NULL) continue;
		if(! neighbours.at(i)->dead) count++;
	}
	return count;
}


double CStrain::sumM(){
	double sum=0.0;
	for(int i=0; i<=rmax; i++){
		sum+=M[i];
	}
	return sum;
}


//calculates the weighted sum of M, the weight is passed 
//to the function through the pointer of an array (chi) 
//with same lenght as M
double CStrain::WeightedSumM(double *chi){
	double sum=0.0;
	for(int i=0; i<=rmax; i++){
		sum+=M[i]*chi[i];
	}
	return sum;
}

// To same some memory we clean M of 
// dead strains
void CStrain::die(){
	N=0.0;
	dead=true;
	delete[] M;
	M=NULL;
}

void CStrain::trim(){
	if(is_leaf)return;
	std::vector<CStrain*>::iterator it, it0;
	int alive_branches=0;
	for(it=neighbours.begin(), it++; it!=neighbours.end(); it++){
		if(!(*it)->is_leaf or !(*it)->dead){
			(*it)->trim();
			alive_branches++;
		}
	}
	if(alive_branches==0) is_leaf=true;

}

//prints the whole tree starting from this node
void CStrain::print(ostream &out){
	out<<ID<< "   "<<N<<"   "<<neighbours.size()-1<<"  ";
	std::vector<CStrain*>::iterator it;
	for(it=neighbours.begin(); it!=neighbours.end(); it++){
		if((*it)==NULL)continue;
		out<<(*it)->ID<< "   ";
	}
	out<<endl;
	for(it=neighbours.begin(), it++; it!=neighbours.end(); it++){
		(*it)->print(out);
	//	out<<endl;
	}
}

//prints N, number of neighbours and their IDs
void CStrain::print_node(ostream &out)const{
	
	out<<ID<< "   "<<N<<"   ";
	out<<fitness<<"  ";
	if(!dead) {
		assert(N>0);
		for(int i=0; i<rmax+1; i++){
	 		out<<M[i]<<"  ";
		}
	}

	if(is_leaf) {
		out<<0<<endl;	
		return;
	}

	std::vector<CStrain*>::const_iterator it=neighbours.begin();
	int neigh=neighbours.size();
	if((*it)==NULL)neigh--;
	out<<neigh<<"  ";
	for(it=neighbours.begin(); it!=neighbours.end(); it++){
		if((*it)==NULL)continue;
		out<<(*it)->ID<< "   ";
	}
	out<<endl;
}
/*
//FIXME
CStrain::CStrain(stringstream &ss){
	stotal++;
	input>>ID>>N;
	fitness=0.0;
	

	dead=false;
	is_leaf=true;

	M=new double[rmax+1];
	for(int i=0;i<=rmax;i++){
		M[i]=0.0;
	}

	int fID;
	input>>fID
	if(f!=NULL){
		f->add_neighbour(this);
	}
	neighbours.push_back(f);
}

void CStrain::read_node(stringsream &ss){

	//Print M's if not dead
	if(!dead) {
		input<<rmax+1<<"  ";
		for(int i=0; i<rmax+1; i++){
	 		input<<M[i]<<"  ";
		}
	}
	else input<<0<<"  ";
	
	std::vector<CStrain*>::iterator it=neighbours.begin();
	int neigh=neighbours.size();
	if((*it)==NULL)neigh--;
	input<<neigh<<"  ";
	for(it=neighbours.begin(); it!=neighbours.end(); it++){
		if((*it)==NULL)continue;
		input<<(*it)->ID<< "   ";
	}
}
*/
void CStrain::print(ostream &out, double x, double y){
	if(neighbours.size()<2)
		out<<"c "<<x<<"  "<<y<<" 0.02"<<endl;
	static double dL=0.05;
	cal_print_spaces();
	double L=print_space;
	out<<"l "<<x<<"  "<<y<<"  ";
	x+=dL;
	out<<"l "<<x<<"   "<<y<<endl;
	y+=L/2.;
	double nn=neighbours.size()-1.;
	for(int i=1; i<neighbours.size(); i++){
		y-=(i-1.)*L/(nn);
		out<<"l "<<x<<"  "<<y<<"  "<<x+dL<<"   "<<y<<endl;
		neighbours.at(i)->print(out,x+dL, y);
	}
}

double CStrain::cal_print_spaces(){
	if(neighbours.size()==1)return print_space;
	
	cerr<< print_space <<endl;
	for(int i=1; i<neighbours.size(); i++)
		print_space+=neighbours.at(i)->cal_print_spaces();

}

#endif
