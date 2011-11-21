#ifndef STRAIN_H
#define STRAIN_H
#include <vector>
#include<assert.h>
#include"parameters.h"
#include"tools.h"

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
	void trim_links();
	CStrain *father();
	void add_neighbour(CStrain *ps){neighbours.push_back(ps);is_leaf=false;};
	void add_link(CStrain *ps, int d){links.push_back(CLink<CStrain>(ps, d));};
	void get_infected(double *sumN, size_t distance=0, CStrain *exclude=NULL);
	void get_infected2(double *sumN, size_t distance=0, CStrain *exclude=NULL);
	CLink<CStrain> to_be_bridged(int distance);
	void make_bridges();
	void get_diversity(double &div, size_t distance=0, CStrain *exclude=NULL);
	void get_diversity2(double &div, size_t distance=0, CStrain *exclude=NULL);
	void print(ostream &out);
	void print_node(ostream &out)const;
	void print_bridges(ostream &out);
	void print(ostream &out, double x, double y);
	void print2(ostream &out, double x, double y);
	double cal_print_widths();
	COffset &cal_offsets();

	double fitness;
	double N;
	int ID;
	bool dead;
	bool is_leaf;
	double *M;
	static unsigned int max_dist;
	std::vector<CStrain*> neighbours;
	std::vector<CLink<CStrain> > links;
	static unsigned stotal;
	COffset offset;
	double print_width;
	int color;
	double x, y;
	private:
	static const double base_print_width=0.01;
};
unsigned int CStrain::stotal=0;

//this is just used in function get_infected()
//to find out what was the maximum distance
//ever reached between the alive nodes
unsigned int CStrain::max_dist=0;

//Constructor take an int for ID and the point of the
//father node
CStrain::CStrain(int i, CStrain *f){
	print_width= base_print_width;
	stotal++;
	ID=i; 
	N=0.0; 
	fitness=0.0;
	if(f!=NULL){
		f->add_neighbour(this);
		f->add_link(this,1);
	}
	neighbours.push_back(f);
	add_link(f,1);

	dead=true;
	is_leaf=true;
	color=1;
	if(ID>=0) {
		SetAlive();
		N=1;
	}
}

void CStrain::SetAlive(){
	dead=false;
	color=0;
	M=new double[rmax+1];
	for(size_t i=0;i<=rmax;i++){
		M[i]=0.0;
	}
}

//Cleans the allocated memory if not yet cleaned 
CStrain::~CStrain(){
	if(M!=NULL) delete[] M; 
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
void CStrain::get_infected(double *sumN, size_t distance, CStrain *exclude){
	if(distance>=rmax) return;
	if(distance>max_dist)max_dist=distance;
	sumN[distance]+=N;

	
	if(father()!=NULL and father()!=exclude) 
		father()->get_infected(sumN, distance+1, this);

	if(is_leaf) return;
	for(size_t i=1; i<neighbours.size(); i++){
		if(neighbours.at(i)==exclude) continue;
		neighbours.at(i)->get_infected(sumN, distance+1, this);
	}
	return;
}

void CStrain::get_infected2(double *sumN, size_t distance, CStrain *exclude){
	sumN[distance]+=N;
	if(distance>=rmax) return;
	if(distance>max_dist)max_dist=distance;

	
	if(links.at(0).head!=NULL and links.at(0).head!=exclude) 
		links.at(0).head->get_infected2(sumN, distance+links.at(0).length, this);

	if(is_leaf) return;
	for(size_t i=1; i<links.size(); i++){
		if(links.at(i).head==exclude) continue;
		links.at(i).head->get_infected2(sumN, distance+links.at(i).length, this);
	}
	return;
}

CLink<CStrain> CStrain::to_be_bridged(int distance ){
	
	int alive_branches=0, branch=-1;;
	for(size_t i=1; i<links.size(); i++){
		if(!(links.at(i).head->is_leaf) or !(links.at(i).head->dead)){
			alive_branches++;
			branch=i;
		}
	}
	if(alive_branches==1 and this->dead){
		return links.at(branch).head->to_be_bridged(distance+links.at(branch).length );
	}
	if(alive_branches==0 and dead){
		 is_leaf=true;
		  if(links.size()>1) color=2;//color for dead branch
	}
	return CLink<CStrain>(this,distance);
	
}

//stopped here, wrong bridges 
void CStrain::make_bridges(){

	int distance=0;
	for(size_t i=1; i<links.size(); i++){
		CLink<CStrain> l=links.at(i).head->to_be_bridged(distance+links.at(i).length);
		if(l.head!=NULL and l.head!=links.at(i).head) {
			links.at(i).head=l.head;
			links.at(i).length=l.length;
			l.head->links.at(0).head=this;
			l.head->links.at(0).length=l.length;
		}
		l.head->make_bridges();
	}
}

void CStrain::get_diversity(double &diversity, size_t distance, CStrain *exclude){
	diversity+=N*distance;

	for(size_t i=0; i<neighbours.size(); i++){
		if(neighbours.at(i)==NULL) continue;
		if(neighbours.at(i)==exclude) continue;
		neighbours.at(i)->get_diversity(diversity, distance+1, this);
	}
	return;
}

void CStrain::get_diversity2(double &diversity, size_t distance, CStrain *exclude){
	diversity+=N*distance;

	for(size_t i=0; i<links.size(); i++){
		if(links.at(i).head==NULL) continue;
		if(links.at(i).head==exclude) continue;
		links.at(i).head->get_diversity2(diversity, distance+links.at(i).length, this);
	}
	return;
}

//Returns the number of alive neighbours of the strain
int CStrain::count_neigh(){
	int count=0;
	for (size_t i=0; i<neighbours.size();i++){
		if(neighbours.at(i)==NULL) continue;
		if(! neighbours.at(i)->dead) count++;
	}
	return count;
}


double CStrain::sumM(){
	double sum=0.0;
	for(size_t i=0; i<=rmax; i++){
		sum+=M[i];
	}
	return sum;
}


//calculates the weighted sum of M, the weight is passed 
//to the function through the pointer of an array (chi) 
//with same lenght as M
double CStrain::WeightedSumM(double *chi){
	double sum=0.0;
	for(size_t i=0; i<=rmax; i++){
		sum+=M[i]*chi[i];
	}
	return sum;
}

// To same some memory we clean M of 
// dead strains
void CStrain::die(){
	assert(!dead);
	N=0.0;
	dead=true;
	color=1;
	delete[] M;
	M=NULL; 
}

void CStrain::trim_links(){
       
       int alive_branches=0;
       for(size_t i=1; i<links.size(); i++){
               if(!(links.at(i).head->is_leaf) or !(links.at(i).head->dead)){
                       links.at(i).head->trim_links();
                       alive_branches++;
               }
       }
       if(alive_branches==0 and dead){
                is_leaf=true;
                 if(links.size()>1) color=2;//color for dead branch
       }
}

void CStrain::trim(){
	if(is_leaf)return;
	std::vector<CStrain*>::iterator it, it0;
	int alive_branches=0;
	assert(neighbours.size()>1);
	for(it=neighbours.begin(), it++; it!=neighbours.end(); it++){
		if(!(*it)->is_leaf or !(*it)->dead){
			(*it)->trim();
			alive_branches++;
		}
	}
	if(alive_branches==0 and dead){
		is_leaf=true;
		if(neighbours.size()>1) color=2;//color for dead branch
	}

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
		for(size_t i=0; i<rmax+1; i++){
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
void CStrain::print2(ostream &out, double x, double y){
	static double dL=0.02;
	out<<"c "<<x<<"  "<<y<<"  "<<base_print_width/4<<"  "<<color<<endl;
	size_t nn=neighbours.size()-1.;
	if(nn<1) return;
	out<<"l "<<x<<"  "<<y<<"  ";
	x+=dL;
	out<<x<<"   "<<y<<endl;
	if(nn==1){
		neighbours.at(1)->print2(out,x, y);
		return;
	}
	double L=offset.width();
	L-=neighbours.at(1)->offset.top;
	L-=neighbours.at(nn)->offset.bottom;
	y+=L/2.;
	out<<"l "<<x<<"  "<<y<<"  ";
	out<<x<<"   "<<y-L<<endl;
	
	for(size_t i=1; i<neighbours.size(); i++){
		//out<<"l "<<x<<"  "<<y<<"  "<<x+dL<<"   "<<y<<endl;
		neighbours.at(i)->print2(out,x, y);
		y-=neighbours.at(i)->offset.bottom;
		if(i<nn)y-=neighbours.at(i+1)->offset.top;
	}
}
void CStrain::print(ostream &out, double x, double y){
	static double dL=0.02;
	out<<"c "<<x<<"  "<<y<<"  "<<base_print_width/4<<"  "<<color<<endl;
	this->x=x;this->y=y;

	if(neighbours.size()<2) return;
	size_t nn=neighbours.size()-1.;
	size_t ind[neighbours.size()];
	mysort(neighbours,ind, neighbours.size());
	if(nn==1){
		out<<"l "<<x<<"  "<<y<<"  "<<x+dL<<"   "<<y<<endl;
		neighbours.at(ind[1])->print(out,x+dL, y);
		return;
	}
	double L=offset.width();
	L-=neighbours.at(ind[1])->offset.top;
	L-=neighbours.at(ind[nn])->offset.bottom;
	y+=L/2.;
	out<<"l "<<x<<"  "<<y<<"  ";
	out<<x<<"   "<<y-L<<endl;
	
	for(size_t i=1; i<=nn; i++){
		out<<"l "<<x<<"  "<<y<<"  "<<x+dL<<"   "<<y<<endl;
		neighbours.at(ind[i])->print(out,x+dL, y);
		y-=neighbours.at(ind[i])->offset.bottom;
		if(i<nn)y-=neighbours.at(ind[i+1])->offset.top;
	}
}
void CStrain::print_bridges(ostream &out){

	for(size_t i=1; i<links.size(); i++){
		CLink<CStrain> &l=links.at(i);
		if(l.length>1)out<<"a "<<x<<"  "<<y<<"  "<<l.head->x<<"   "<<l.head->y<<endl;
		l.head->print_bridges(out);
		}
}
COffset &CStrain::cal_offsets(){
	offset.top=0.0;
	offset.bottom=0.0;
	size_t n=neighbours.size()-1;//no. of children
	cal_print_widths();
	if(n==0){
		offset.top=print_width/2;
		offset.bottom=print_width/2;
		return offset;
		}
	
	size_t ind[neighbours.size()];
	mysort(neighbours,ind, neighbours.size());

	COffset &off1=neighbours.at(ind[1])->cal_offsets();
	if(n==1){
		offset=off1;
		return offset;
		}
	
	double size=0;
	size+=off1.bottom;
	for(size_t i=2; i<n; i++){
		COffset &off=neighbours.at(ind[i])->cal_offsets();
		size+=off.width();
		}
	COffset &off2=neighbours.at(ind[n])->cal_offsets();
	size+=off2.top;

	offset.top=size/2+off1.top;
	offset.bottom=size/2+off2.bottom;

	return offset;
}

double CStrain::cal_print_widths(){
	if(neighbours.size()==1)return base_print_width;
	
	print_width=0;
	for(size_t i=1; i<neighbours.size(); i++)
		print_width+=neighbours.at(i)->cal_print_widths();

	return print_width;
}

#endif
