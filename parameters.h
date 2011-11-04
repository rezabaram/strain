#ifndef PARAMETERS_H
#define PARAMETERS_H

const int rmax=50;
const double Npop=1e+7;
const double R0=1.5;
const double nu=365.0/4.0;
const int L=120; // 60 epitope codons in the HA1 domain
const double mu=5.8*1e-3; // "ideal" mut_rate = L*mu/365 = 0.001906849;
const double mut_rate=300.*1e-4;
const double beta=R0*nu;
const double f0=nu*(R0-1.); // beta-nu
const double beta0=beta/(Npop*(beta-nu));
const double dt=1./365.;
const double tMax=1000.;

#endif /* PARAMETERS_H */