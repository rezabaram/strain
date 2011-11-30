#ifndef IO_H
#define IO_H 
#include<iostream>
#include<string>
#include<assert.h>
#include"strain.h"
using namespace std;

/*
void SaveState(ostream &out, const vector<CStrain*> &strains){
	out<<strains.size()<<"   "<<rmax<<endl;
	vector<CStrain*>::const_iterator it;
	for(it=strains.begin(); it!=strains.end(); it++){
		(*it)->print_node(out);
	}
}
*/

/*
void ReadState(istream &input, vector<CStrain*> &strains){
	strains.clear();
	string line;

	//read the header of the file
	int ntotal_read;
	if(getline(input,line)){
		stringstream ss(line);
		ss>>ntotal_read;
		ss>>rmax;
		}

	//first reading all the lines and
	//create backbone of the strain
	vector<string> lines;
	while(getline(input,line))
	{
		lines.push_back(line);
		CStrain *s=new CStrain(-1,NULL);
		strains.push_back(s);
	}

	//Insert the line string into a stream
	for(size_t i=0; i<strains.size(); i++){
		stringstream ss(lines.at(i));
		

		//alias for current strain (for conciseness)
		CStrain &s=*strains.at(i);
		
		//readin ID, N, etc

		ss>> s.ID; 
		ss>> s.N; 
		ss>> s.fitness; 

		if(s.N>0) s.SetAlive();

		if(! s.dead){
			//read M's
			for(size_t i=0; i<=rmax; i++){
				ss>>s.M[i];
			}
		}
		size_t nneigh;
		ss>>nneigh;
		for(size_t i=0; i<nneigh; i++){
			int neigh_id;
			ss>>neigh_id;
			s.add_neighbour(strains.at(neigh_id));
		}
	}

}
*/

#endif /* IO_H */
