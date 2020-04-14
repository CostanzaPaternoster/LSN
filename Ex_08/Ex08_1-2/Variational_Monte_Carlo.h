/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __fluid_
#define __fluid_

using namespace std;

//Random numbers
#include "random.h"
int seed[4];
Random rnd("seed.in");

//parameters, observables
const int m_props=1000;
int n_props, iv, ik, ie;
double walker[m_props];

// averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props], err[m_props];
double stima[m_props];

//configuration
double x;

//parameters
double mean, sigma;

// simulation
int nstep, nblk, eqstep;
double delta;
bool print_xyz, print_instant, print_blocks;

//variation
int nparam=1;
bool mean_set = false;
bool sigma_set = false;
double energy_new, energy_old, energy_opt, err_opt;
int nvar;
double deltamean_var,deltasigma_var, acceptedmean_var,acceptedsigma_var, attemptedmean_var,attemptedsigma_var;
double mean_old, sigma_old, mean_opt, sigma_opt;

//histogram
int  nbins;
double histo_start, histo_end, bin;


//functions
void Input(void);
void SetParam(int);
void Equilibrate();
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
void ConfXYZ();
void Measure(void);
double Potential(double);
double PsiT(double);
void Final(int,int);
double Error(double,double,int);
void Variation();
void Histogram();
void FinalizeHistogram(int);

#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
