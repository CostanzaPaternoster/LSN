/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "funz.h"

using namespace std;

int main()
{ 
  for( int T_step=0; T_step <=nbins; T_step++) // cycle over temperatures
  {
        Input(T_step); //Inizialization
        //Equilibrate(); //Equilibration
      
        for(int iblk=1; iblk <= nblk; ++iblk) // MC simulation
        {
            Reset(iblk);   //Reset block averages
            for(int istep=1; istep <= nstep; ++istep)
            {
                Move(metro);
                Measure();
                Accumulate(); //Update block averages
            }
            Averages(iblk);   //Print results for current block
        }
        ConfFinal(); //Write final configuration
  }

  return 0;
}


void Input(int T_step)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp_max;
  ReadInput >> temp_min;
  ReadInput >> nbins;

  temp = temp_max - (temp_max-temp_min)*T_step/(double)nbins;
    
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;
  
  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;
    
  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
    
  ReadInput.close();

//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables
  
  if(T_step==0) Setting(); // choose whether to start from default (inf) or previous config
  else mode2 = 1;
  
//Prepare initial config. starting with the previous one
  if (mode2){
    
    ifstream ReadConfig("config.final");  
    if (ReadConfig.is_open()){
        cout << "Preparing initial spin configuration from file 'config.final' ." << endl << endl;
        for (int i=0; i<nspin; ++i){
            ReadConfig >> s[i];
        }        
    } else{
        cerr << "PROBLEM: Unable to open 'config.final' " << endl << endl;
        //Prepare default initial configuration (T = inf)
        cout << "Preparing default initial spin configuration (T = inf) " << endl << endl;
        for (int i=0; i<nspin; ++i)
            {
                if(rnd.Rannyu() >= 0.5) s[i] = 1;
                else s[i] = -1;
            }
        }
        ReadConfig.close();    
 }
 else{
    
//Prepare default initial configuration (T = inf)
     cout << "Preparing default initial spin configuration (T = inf) " << endl << endl;
     for (int i=0; i<nspin; ++i)
    {
        if(rnd.Rannyu() >= 0.5) s[i] = 1;
        else s[i] = -1;
    }
 }
    
//Evaluate energy etc. of the initial configuration (aggiungere C,M,X)
Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}


void Move(int metro)
{
  int o;
  //double p, energy_old, energy_new, sm;
  //double energy_up, energy_down;

  for(int i=0; i<nspin; ++i)
  {
    //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1) Metropolis(o);  //Metropolis
    else Gibbs(o); //Gibbs sampling
      
    attempted += 1;
  }    
}

double Boltzmann(int sm, int ip) // not used
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Metropolis(int o)
{
    //Trial: flip the randomly selected spin
    double s_move = - s[o];
    
    //Evaluate the energy difference between the actual and the trial spin config.
    double deltaE = -2 * J * s_move * ( s[Pbc(o-1)] + s[Pbc(o+1)] ) -2 * h * s_move;
    
    //Accept the move with p = min[1,exp{-deltaE * beta}]
    double p = exp(- beta * deltaE);
    if(p > 1){
        s[o] = s_move; //Accept move: flip the spin s[o] 
        accepted += 1;
    }
    else{
        double r = rnd.Rannyu();
        if(r < p){ 
            s[o] = s_move; //Accept move: flip the spin s[o]
            accepted += 1;
        }
    }  
}

void Gibbs(int o)
{
    // Compute the probability p to assign value +1 at spin s[o]
    double p = 1. / ( 1 + exp( -2 * beta * (J * (s[Pbc(o-1)] + s[Pbc(o+1)]) + h )) );
    
    // Assign value +1 at spin s[o] with probability p
    double r = rnd.Rannyu();
    
    if( r < p ) s[o] = 1;
    else s[o] = -1;   
}

void Measure(void)
{

  double u = 0.0, m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     m += s[i];

  }
  walker[iu] = u;
  walker[ic] = u*u;
  walker[im] = m;
  walker[ix] = m*m;

}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    const int wd=12;
    ofstream Ene,Heat, Mag, Chi;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    //Ene.open("output.ene",ios::app);    // stampa i valori medi per ogni blocco (distattivato)
    //Heat.open("output.heat",ios::app);
    //Mag.open("output.mag",ios::app);
    //Chi.open("output.chi",ios::app);
    
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    
    stima_c = beta*beta*(blk_av[ic]/blk_norm - pow(blk_av[iu]/blk_norm, 2))/(double)nspin; //Energy
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c=Error(glob_av[ic],glob_av2[ic],iblk); // Heat capacity
    
    stima_m = (blk_av[im]/blk_norm)/(double)nspin; //Magnetisation
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im],glob_av2[im],iblk);
    
    stima_x = beta*(blk_av[ix]/blk_norm - pow(blk_av[im]/blk_norm, 2))/(double)nspin; //Magnetic susceptibility
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);
    
    //Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    //Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
    //Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
    //Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
    
    //Ene.close();
    //Heat.close();
    //Mag.close();
    //Chi.close();
        
    if(iblk == nblk){ // stampa i valori medi all'ultimo blocco
        
        Ene.open("Temp/"+fold_name[metro]+"/Energy.txt",ios::app);
        Heat.open("Temp/"+fold_name[metro]+"/Heat.txt",ios::app);
        Mag.open("Temp/"+fold_name[metro]+"/Mag.txt",ios::app);
        Chi.open("Temp/"+fold_name[metro]+"/Chi.txt",ios::app);
        
        Ene << setw(wd) << temp << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
        Heat << setw(wd) << temp << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
        Mag << setw(wd) << temp << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
        Chi << setw(wd) << temp << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
        
        Ene.close();
        Heat.close();
        Mag.close();
        Chi.close();
    }
    

// INCLUDE YOUR CODE HERE

    cout << "----------------------------" << endl << endl;
}

void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==0) return 0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

// Function to choose equilibration/simulation and starting config
void Setting(void){
   
  char b;
//Choose whether to start from a new configuration (T = inf) or from a previous one
  cout << "Do you want to start from the previous configuration [y/n]" << endl << endl;
  cin >> b;

  while( b != 'y' && b != 'n'){
     cout << "Input not valid. Choose y to start from previous configuration or n to start from the default configuration (T = inf). " << endl << endl;
     cin >> b;
 }
    
 if(b == 'y') mode2 = 1;
 else mode2 = 0;
    
 return;
}

// Equilibration function
void Equilibrate(void){
    
    char a,b;
    double var[4];
    double u,m;
    ofstream Ene;
    
    // Choose whether to equilibrate the system or to start the entire simulation
    cout << "Do you want to equilibrate the system? [y/n]" << endl << endl;
    cin >> a;
    
    while( a != 'y' && a != 'n'){
        cout << "Input not valid. Choose y to equilibrate the system or n to start the simulation. " << endl << endl;
        cin >> a;
    }
    
    if( a == 'y' ){
        cout << "Choose the number of MC steps for equilibration: " << endl << endl;
        cin >> nstep_eq;
        cout << "****** EQUILIBRATION  ****** "<< endl << endl;
        mode = 1;
    }else{
        cout << "****** SIMULATION  ****** "<< endl << endl;
        mode = 0;
    }
    
    while(mode){
        
        Ene.open("Temp/"+fold_name[metro]+"Equilibrazione/instant_ene.txt");
        
        for (int j = 0; j < 4; j++ ) var[j] = 0;
        
        for(int ieq = 0; ieq < nstep_eq; ieq++){
	  
	       u = 0;
           m = 0;
           for(int i = 0; i < nspin; i++) {
		      
              u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]); // energy
              m += s[i];
	       }
            
           var[0] += u; // energy
		   var[1] += m; // magnetization
		   var[2] += u*u; // heat capacity
		   var[3] += m*m; // magnetic susceptibility
            
           Ene << u << endl;
        
           Move(metro);
	       
	   }
        
       Ene.close();
        
       var[2] = beta*beta*(var[2]/(double)nstep_eq - var[0]*var[0]/(double)(nstep_eq*nstep_eq));
       var[3] = beta*(var[3]/(double)nstep_eq - var[1]*var[1]/(double)(nstep_eq*nstep_eq));
       
    // Printing to video
       cout << "The mean energy per particle is: " << (var[0]/(double)nstep_eq)/(double)nspin << endl;
       cout << "The mean magnetisation per particle is: " << (var[1]/(double)nstep_eq)/(double)nspin  << endl;
       cout << "The mean heat capacity per particle is: " << (var[2])/(double)nspin  << endl;
       cout << "The mean magnetic susceptibility is: " << (var[3])/(double)nspin  << endl;
        
        cout << "Do you want to start the simulation? [y/n]. If n the program will go on with the equilibration. " << endl << endl;
        
        cin >> b;
        
        while( b != 'y' && b != 'n'){
            cout << endl << "Invalid input. Choose y or n" << endl;
            cin >> b;
        }
        if (b == 'y') {
            cout << "****** SIMULATION  ****** "<< endl << endl; 
            mode = 0;
        }
        else mode= 1;
        
        if(mode){
            ConfFinal();     
        }
    }
    return;    
    
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
