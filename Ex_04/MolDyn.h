/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

int n_eq; // number of equilibration steps
bool mode1; // 1: equilibration - 0: simulation
bool mode2; // 1: start from previous config. - 0: start from fcc

//parameters, observables
const int m_props=500;
int n_props;
int iv,ik,it,ie,ip,igofr;
double insta_pot, insta_kin, insta_etot, insta_temp, insta_pres;
double stima_pot, stima_kin, stima_etot, stima_temp, stima_pres, stima_gdir;
double err_pot, err_kin, err_temp, err_etot, err_pres, err_gdir;
double walker[m_props];

// averages
double acc,att;
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, nblk, block_steps, iprint, seed, nconf;
double nbins, bin_size, delta_V; // Histogram g(r)
double delta;

//pigreco
const double pi=3.1415927;

//functions
void Input(void);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);

void Reset(int);
void Accumulate(void);
void Averages(int);

double Force(int, int);
double Pbc(double);
void Equilibrate(void); // equilibrate the system before starting the simulation
void Rescale(void); // rescale the velocities according to the desired temperature
void StartingConf(void); // set the starting configuration: fcc or from previous simulation 
void Setting(void); // choose: equilibration/simulation, starting configuration
double GetTemp(void); // get the actual temperature
double Error(double,double,int);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
