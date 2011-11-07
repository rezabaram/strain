#ifndef PARAMETERS_H
#define PARAMETERS_H

const int rmax=500;
const int N0=10;
const double Npop=1e+5;
const double R0=1.5;
const unsigned int inf_period=4;
const double nu=365.0/(double)inf_period;
const int L=120; // 60 epitope codons in the HA1 domain
const double mu=5.8*1e-3; // "ideal" mut_rate = L*mu/365 = 0.001906849;
const double mut_rate=500.*1e-4;
const double beta=R0*nu;
const double f0=nu*(R0-1.); // beta-nu
const double beta0=beta/(Npop*(beta-nu));
const double dt=1./365.;
const double tMax=1000.;

#endif /* PARAMETERS_H */
