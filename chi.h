#ifndef CHI_H
#define CHI_H 

double *chi;

double chi_at_d(double d){
	double chi;
	chi=1./(d*d*d+1);

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
#endif /* CHI_H */
