#include<iostream>
#include"strain.h"

using namespace std;

int main(){

double *a;
a=new double[10];
cerr<< a[4] <<endl;

CStrain *n[10];
n[0]=new CStrain(0);
n[1]=new CStrain(1,n[0]);
n[2]=new CStrain(2,n[0]);
n[3]=new CStrain(3,n[0]);
n[4]=new CStrain(4,n[1]);
n[5]=new CStrain(5,n[1]);
n[6]=new CStrain(6,n[1]);
n[7]=new CStrain(7,n[2]);
n[8]=new CStrain(8,n[2]);
n[9]=new CStrain(9,n[3]);

double m[5+1];
for(int i=0; i<10; i++){
	continue;
	m[0]=0;
	m[1]=0;
	m[2]=0;
	m[3]=0;
	m[4]=0;
	m[5]=0;
	n[i]->get_infected(m);
	for(int j=0; j<rmax+1; j++){
		cout<< m[j] <<"  ";
	}
	cout<<endl;
	}
return 0;
}
