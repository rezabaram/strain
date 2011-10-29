#ifndef PARAMETERS_H
#define PARAMETERS_H 

const int rmax=4;
const double N=1e+5;
const double R0=1.5;
const double nu=365.0/4.0;
const double beta=R0*nu;
const double f0=beta-nu;
const double beta0=beta/(N*(beta-nu));
const double dt=1./365;
const double tMax=1;

#endif /* PARAMETERS_H */
