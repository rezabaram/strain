#ifndef STRAIN_H
#define STRAIN_H
#include <vector>
#include<assert.h>
#include"parameters.h"
#include"tools.h"

using namespace std;


class CStrain{
	public:
	explicit CStrain(int i, CStrain *f, double im_d=1);
	~CStrain();
	void SetAlive();
	double sumM();
	//double WeightedSumM(double chi(double) );
	int count_neigh();
	void die();
	void trim();
	void trim_links();
	CStrain *father();
	void add_neighbour(CStrain *ps, int d, double im_d=1){neighbours.push_back(CLink<CStrain>(ps, d, im_d));is_leaf=false;};
	void add_link(CStrain *ps, int d, double im_d=1){links.push_back(CLink<CStrain>(ps, d, im_d));};
	CLink<CStrain> to_be_bridged(int distance, double imm_dist);
	void make_bridges();
	void get_diversity(double &div, size_t distance=0, CStrain *exclude=NULL);
	void get_diversity2(double &div, size_t distance=0, CStrain *exclude=NULL);
	double WeightedSumM0(double chi(double), double distance=0, CStrain *exclude=NULL);
	double WeightedSumM(double chi(double), double distance=0, CStrain *exclude=NULL);
	double calSubN();
	void setFreq(double f,double tt);
	void print(ostream &out);
	void print_node(ostream &out)const;
	void print_bridges(ostream &out);
	void print(ostream &out, double x, double y);
	void print2(ostream &out, double x, double y);
	double cal_print_widths();
	COffset &cal_offsets();
	double Freq, maxFreq;
	double SubN;
	double accN;
	double fitness;
	unsigned int red_m;
	double N;
	int ID;
	double crtime;
	double fixtime;
	bool notfixed;
	bool dead;
	bool is_leaf;
	int mut_type;
	static unsigned int max_dist;
	std::vector<CLink<CStrain> > neighbours;
	std::vector<CLink<CStrain> > links;
	static unsigned stotal;
	COffset offset;
	double print_width;
	int color;
	double x, y;
	double M0;
	private:
	static const double base_print_width=0.005;
};
unsigned int CStrain::stotal=0;

//this is just used in function get_infected()
//to find out what was the maximum distance
//ever reached between the alive nodes
unsigned int CStrain::max_dist=0;

//Constructor take an int for ID and the point of the
//father node
CStrain::CStrain(int i, CStrain *f, double im_d){
	print_width = base_print_width;
	stotal++;
	ID=i;
	N=0.0;
	fitness=0.0;
	accN=0.0;
	M0=0.0;
	SubN=0.;
	Freq=0.;
	maxFreq=0.;
	notfixed=true;
	red_m=0;
	if(f!=NULL){
		f->add_neighbour(this, 1, im_d);
		f->add_link(this,1, im_d);
		red_m=f->red_m;
	}
	add_neighbour(f,1, im_d);
	add_link(f,1, im_d);
	

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
}

//Cleans the allocated memory if not yet cleaned 
CStrain::~CStrain(){
	stotal--;
}

//Returns the pointer to the father if not root node
CStrain* CStrain::father(){
	return neighbours.at(0).head;
}

void CStrain::setFreq(double f,double tt){
	  int max_bin=(nbins-1)*maxFreq;
	  if( max_bin==(nbins-1) && notfixed){ fixtime=tt; notfixed=false;}
          Freq=f;
          if(Freq>maxFreq)maxFreq=Freq;
}

//Returns the sum of N of all strains at distances up to rmax
//sumN[0] is the N of the strain itself
//sumN[1] is the sum of N's of all strains at distance 1
//....
//sumN[rmax] is the sum of N's of all strains at distance rmax
double CStrain::WeightedSumM0(double chi(double), double distance, CStrain *exclude){
	if(distance>=rmax) return 0;
	if(distance>max_dist)max_dist=distance;
	double weightedsum=chi(distance)*accN;

	
	if(father()!=NULL and father()!=exclude) 
		weightedsum+=neighbours.at(0).head->WeightedSumM0(chi, distance+neighbours.at(0).immune_distance, this);

	//if(is_leaf) return weightedsum;
	for(size_t i=1; i<neighbours.size(); i++){
		if(neighbours.at(i).head==exclude) continue;
		weightedsum+=neighbours.at(i).head->WeightedSumM0(chi, distance+neighbours.at(i).immune_distance, this);
	}
	return weightedsum;
}

double CStrain::WeightedSumM(double chi(double), double distance, CStrain *exclude){
	if(distance>=rmax) return 0;
	if(distance>max_dist)max_dist=distance;
	double weightedsum=chi(distance)*N;

	
	if(links.at(0).head!=NULL and links.at(0).head!=exclude)
		weightedsum+=links.at(0).head->WeightedSumM(chi, distance+links.at(0).immune_distance, this);

	if(is_leaf) return weightedsum;

	for(size_t i=1; i<links.size(); i++){
		if(links.at(i).head==exclude) continue;
		weightedsum+=links.at(i).head->WeightedSumM(chi, distance+links.at(i).immune_distance, this);
	}
	return weightedsum;
}



CLink<CStrain> CStrain::to_be_bridged(int distance, double imm_dist ){
	
	int alive_branches=0, branch=-1;;
	for(size_t i=1; i<links.size(); i++){
		if(!(links.at(i).head->is_leaf) or !(links.at(i).head->dead)){
			alive_branches++;
			branch=i;
		}
	}
	if(alive_branches==1 and this->dead){
		return links.at(branch).head->to_be_bridged(distance+links.at(branch).length, imm_dist+links.at(branch).immune_distance );
	}
	if(alive_branches==0 and dead){
		 is_leaf=true;
		  if(links.size()>1) color=2;//color for dead branch
	}
	return CLink<CStrain>(this,distance, imm_dist);
	
}
 
void CStrain::make_bridges(){
	assert(links.size()==neighbours.size());
	int distance=0;
	double imm_dist=0;
	for(size_t i=1; i<links.size(); i++){
		CLink<CStrain> l=links.at(i).head->to_be_bridged(distance+links.at(i).length, imm_dist+links.at(i).immune_distance);
		if(l.head!=NULL and l.head!=links.at(i).head) {
			links.at(i).head=l.head;
			links.at(i).length=l.length;
			links.at(i).immune_distance=l.immune_distance;
			l.head->links.at(0).head=this;
			l.head->links.at(0).length=l.length;
			l.head->links.at(0).immune_distance=l.immune_distance;
		}
		l.head->make_bridges();
	}
}

void CStrain::get_diversity(double &diversity, size_t distance, CStrain *exclude){
	diversity+=N*distance;

	for(size_t i=0; i<neighbours.size(); i++){
		if(neighbours.at(i).head==NULL) continue;
		if(neighbours.at(i).head==exclude) continue;
		neighbours.at(i).head->get_diversity(diversity, distance+1, this);
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
		if(neighbours.at(i).head==NULL) continue;
		if(! neighbours.at(i).head->dead) count++;
	}
	return count;
}

//calculates the weighted sum of M, the weight is passed
//to the function through the pointer of an array (chi)
//with same length as M
/*
double CStrain::WeightedSumM(double chi(double) ){
	double sum=0.0;
	for(size_t i=0; i<=rmax; i++){
		sum+=M[i]*chi(i);
	}
	return sum;
}
*/
// To save some memory we clean M of 
// dead strains
void CStrain::die(){
	assert(!dead);
	N=0.0;
	dead=true;
	color=1;
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

double CStrain::calSubN(){
	SubN=N;
	std::vector<CLink<CStrain> >::iterator it;
	for(it=neighbours.begin(), it++; it!=neighbours.end(); it++){
		SubN+=(*it).head->calSubN();
	}
	
	return SubN;
}

void CStrain::trim(){
	if(is_leaf)return;
	std::vector<CLink<CStrain> >::iterator it, it0;
	int alive_branches=0;
	assert(neighbours.size()>1);
	for(it=neighbours.begin(), it++; it!=neighbours.end(); it++){
		if(!(*it).head->is_leaf or !(*it).head->dead){
			(*it).head->trim();
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
	std::vector<CLink<CStrain> >::iterator it;
	for(it=neighbours.begin(); it!=neighbours.end(); it++){
		if((*it).head==NULL)continue;
		out<<(*it).head->ID<< "   ";
	}
	out<<endl;
	for(it=neighbours.begin(), it++; it!=neighbours.end(); it++){
		(*it).head->print(out);
	//	out<<endl;
	}
}

//prints N, number of neighbours and their IDs
/*
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
*/
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
		neighbours.at(1).head->print2(out,x, y);
		return;
	}
	double L=offset.width();
	L-=neighbours.at(1).head->offset.top;
	L-=neighbours.at(nn).head->offset.bottom;
	y+=L/2.;
	out<<"l "<<x<<"  "<<y<<"  ";
	out<<x<<"   "<<y-L<<endl;
	
	for(size_t i=1; i<neighbours.size(); i++){
		//out<<"l "<<x<<"  "<<y<<"  "<<x+dL<<"   "<<y<<endl;
		neighbours.at(i).head->print2(out,x, y);
		y-=neighbours.at(i).head->offset.bottom;
		if(i<nn)y-=neighbours.at(i+1).head->offset.top;
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
		neighbours.at(ind[1]).head->print(out,x+dL, y);
		return;
	}
	double L=offset.width();
	L-=neighbours.at(ind[1]).head->offset.top;
	L-=neighbours.at(ind[nn]).head->offset.bottom;
	y+=L/2.;
	out<<"l "<<x<<"  "<<y<<"  ";
	out<<x<<"   "<<y-L<<endl;
	
	for(size_t i=1; i<=nn; i++){
		out<<"l "<<x<<"  "<<y<<"  "<<x+dL<<"   "<<y<<endl;
		neighbours.at(ind[i]).head->print(out,x+dL, y);
		y-=neighbours.at(ind[i]).head->offset.bottom;
		if(i<nn)y-=neighbours.at(ind[i+1]).head->offset.top;
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

	COffset &off1=neighbours.at(ind[1]).head->cal_offsets();
	if(n==1){
		offset=off1;
		return offset;
		}
	
	double size=0;
	size+=off1.bottom;
	for(size_t i=2; i<n; i++){
		COffset &off=neighbours.at(ind[i]).head->cal_offsets();
		size+=off.width();
		}
	COffset &off2=neighbours.at(ind[n]).head->cal_offsets();
	size+=off2.top;

	offset.top=size/2+off1.top;
	offset.bottom=size/2+off2.bottom;

	return offset;
}

double CStrain::cal_print_widths(){
	if(neighbours.size()==1)return base_print_width;
	
	print_width=0;
	for(size_t i=1; i<neighbours.size(); i++)
		print_width+=neighbours.at(i).head->cal_print_widths();

	return print_width;
}

#endif
