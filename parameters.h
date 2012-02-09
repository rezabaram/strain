#ifndef PARAMETERS_H
#define PARAMETERS_H

const double Cf=0.;
const double cp=0.01;
size_t rmax=10;
const int N0=10;
const double Npop=1e+6;
const double R0=1.1;
const unsigned int inf_period=4;
const double nu=365.0/(double)inf_period; // units 1/year
const int L=120;//120+160; // 60 epitope codons in the HA1 domain? 80 non-epitope codons
const double mu=5.8*1e-3; // units 1/year
double mut_rate=mu*L/365.;
const double beta=R0*nu;
const double f0=nu*(R0-1.); // beta-nu 1.7
const double beta0=beta/(Npop*(beta-nu)); // 0.00001/f0
const double dt=1./365.;
const double tMax=100.;
const int nbins=100;

// ./plot.sh [0:10] 1:2 single00*

// git stash -- to delete current
// git pull reza-github ganna:ganna -- to download
// git log -- to check the version

// git commit -a -m "your comment"
// git push origin ganna:ganna

// git config --global color.log
// git config --global color.diff always

// git add signal.h

// git checkout branch

// ./strain &
// kill -30 27241 (the latter is the process number)
// ps
// fg

//gawk '{ print $0 > "line"NR}' fitness

//:wq

//mplayer N_distribution.avi -loop 5
//mplayer N_distribution.avi -speed 0.5

#endif /* PARAMETERS_H */
