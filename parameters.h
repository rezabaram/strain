#ifndef PARAMETERS_H
#define PARAMETERS_H

size_t rmax=10000;
const int N0=10;
const double Npop=1e+6;
const double R0=1.2;
const unsigned int inf_period=5;
const double nu=365.0/(double)inf_period;
const int L=120; // 60 epitope codons in the HA1 domain
const double mu=5.8*1e-3; // "ideal" mut_rate = L*mu/365 = 0.001906849;
//const 
double mut_rate=0.009;
const double beta=R0*nu;
const double f0=nu*(R0-1.); // beta-nu
const double beta0=beta/(Npop*(beta-nu));
const double dt=1./365.;
const double tMax=1000.;

// ./plot.sh [0:10] 1:2 single00*

// git stash -- to delete current
// git pull reza-github ganna:ganna -- to download
// git log -- to check the version

// git commit -a -m "your comment"
// git push origin ganna:ganna

// git config --global color.log
// git config --global color.diff always

// git add signal.h

// ./strain &
// kill -30 27241 (the latter is the process number)
// ps
// fg

#endif /* PARAMETERS_H */
