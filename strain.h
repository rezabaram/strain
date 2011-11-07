#ifndef STRAIN_H
#define STRAIN_H
#include <vector>
#include"parameters.h"

using namespace std;

class CStrain{
	public:
	CStrain(int i, CStrain *f);
	~CStrain();
	double sumM();
	double WeightedSumM(double *chi);
	int count_neigh();
	void die();
	void trim();
	CStrain *father();
	void add_neighbour(CStrain *ps){neighbours.push_back(ps);is_leaf=false;};
	void get_infected(double *sumN, int distance=0, CStrain *exclude=NULL);
	void print(ostream &out);

	double fitness;
	double N;
	int ID;
	bool dead;
	bool is_leaf;
	double *M;
	static int max_dist;
	std::vector<CStrain*> neighbours;
	static unsigned int stotal;
};
unsigned int CStrain::stotal=0;

//this is just used in function get_infected()
//to find out what was the maximum distance
//ever reached between the alive nodes
int CStrain::max_dist=0;

//Constructor take an int for ID and the point of the
//father node
CStrain::CStrain(int i, CStrain *f){
	stotal++;
	ID=i; 
	N=0.0; 
	fitness=0.0;
	if(f!=NULL){
		f->add_neighbour(this);
	}
	neighbours.push_back(f);

	M=new double[rmax+1];
	for(int i=0;i<=rmax;i++){
		M[i]=0.0;
	}
	dead=false;
	is_leaf=true;
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

	for(int i=0; i<neighbours.size(); i++){
		if(neighbours.at(i)==NULL) continue;
		if(neighbours.at(i)==exclude) continue;
		neighbours.at(i)->get_infected(sumN, distance+1, this);
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
	std::vector<CStrain*>::iterator it, it0;
	for(it=neighbours.begin(), it++; it!=neighbours.end(); it++){
		if((*it)->is_leaf and (*it)->dead) {
			//cut_branch=true;
			it0=it;
			it--;
			delete (*it0);
			neighbours.erase(it0);
		}
		else{
			(*it)->trim();
		}
	}
	if(neighbours.size()==1) is_leaf=true;

}

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

/*
void print(ostream &out, double x, double y, double &dx){
	static double dx=0.05;
	out<<x<<"  "<<y<<"  "<<x+dx<<"   "<<y<<endl;
	double l=(neighbours.size()-1)*dx;
	out<<x<<"  "<<y<<"  "<<x+dx<<"   "<<y<<endl;
	for(int i=0; i<neighbours.size(); i++){
		out<<x<<"  "<<y-l/2.0+i*dx<<"  "<<x+dx<<"   "<<y-l/2.0+i*dx<endl;
	}
}
*/

#endif
