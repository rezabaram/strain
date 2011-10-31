#ifndef PARAMETERS_H
#define PARAMETERS_H 

const int rmax=50;
const double Npop=1e+5;
const double R0=1.5;
const double nu=365.0/4.0;
const double mut_rate=100.0*1e-3;
const double beta=R0*nu;
const double f0=beta-nu;
const double beta0=beta/(Npop*(beta-nu));
const double dt=1./365;
const double tMax=100;

#endif /* PARAMETERS_H */
