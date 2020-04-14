#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <string>
#include "Variational_Monte_Carlo.h"

using namespace std;

int main()
{
  Input();  // Inizialization
  
  for(int i=0; i<nparam;i++){ // Cycle over starting parameters mu and sigma
      
      SetParam(i);  // Set starting mu and sigma
      if(nvar==0) Equilibrate();    // Equilibrate the system (not during variation)

      for(int ivar = 0; ivar <= nvar; ++ivar){
          if(ivar>0) Variation();   // Make a move in paramters' space

          for(int iblk=1; iblk <= nblk; ++iblk){ // Simulation with new parameters
              Reset(iblk);
              for(int istep=1; istep <= nstep; ++istep){
                  Move();
                  Measure();
                  Accumulate();
                  if((istep%10 == 0) && (print_xyz))  ConfXYZ();//Write actual configuration in XYZ format
              }
              Averages(iblk);
          }
          if (nvar>0) Final(ivar,i); // Accept/reject new parameters
      }
    }
  return 0;
}


void Input(void)
{
  ifstream ReadInput,ReadConf;
  cout << "Variational Monte Carlo simulation             " << endl << endl;
  cout << "Ground state of 1D quantum particle in the external potential:" << endl;
  cout << "   V(x) = x^4 - 2.5x^2" << endl << endl;
  cout << "Test wave function:" << endl;
  cout << "   Psi_T(x) = exp[ -(x-mean)^2/2sigma^2 ] + exp[ -(x+mean)^2/2sigma^2 ]" << endl << endl;


//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> x;

  ReadInput >> mean;

  ReadInput >> sigma;

  ReadInput >> delta;

  ReadInput >> nblk;

  ReadInput >> nstep;
    
  ReadInput >> eqstep;

  ReadInput >> deltamean_var;
    
  ReadInput >> deltasigma_var;

  ReadInput >> nvar;

  ReadInput >> print_xyz;

  ReadInput >> print_instant;

  ReadInput >> print_blocks;
    
  ReadInput >> nbins;
      
  ReadInput >> histo_start;				
  
  ReadInput >> histo_end;
    
  bin = (histo_end-histo_start)/nbins;

  cout << "The program perform Metropolis moves with uniform translations" << endl;
  cout << "Moves parameter = " << delta << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl;
  
  if(nvar>0){
      ifstream  file;
      string line;
      file.open("start_param.txt");
      
      while(!file.eof()){
          getline(file,line);
          nparam += 1;
      }
      nparam -= 1;
      
      file.close();
      
      cout << "Starting " << nparam << " Variational Monte Carlo Simulations" << endl << endl;
  }
    
  ReadInput.close();

//Prepare arrays for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic
  ie = 2; //Total

  if(nvar != 0) nbins = 0;
  n_props = 3 + nbins; //Number of observables + number of histogram bins

}

void Equilibrate(){
    
    cout << "Starting " << eqstep << " equilibration steps." << endl << endl;
    
    for(int i=0; i<eqstep; i++) Move();
   
    cout << "Starting simulation" << endl << endl;
}


void Move(void)
{
  double p, p_old, p_new;
  double xold,xnew;

  //Old
    xold = x;
    p_old = PsiT(xold);

  //New
    xnew = x + delta*rnd.Rannyu(-1,1);
    p_new = PsiT(xnew);

  //Metropolis test
    p = p_new/p_old;
    if(p >= rnd.Rannyu())
    {
       x = xnew;
       accepted += 1;
    }
    attempted += 1;
}

double PsiT(double x){
  
  return exp( -(x-mean)*(x-mean)/(sigma*sigma) ) + exp( -(x+mean)*(x+mean)/(sigma*sigma) ) + 2*exp( -(x*x+mean*mean)/(sigma*sigma) );

}

double Potential(double x){
  return x*x*x*x - 2.5*x*x;
}

double Kinetic(double x){
  double m = (x-mean)*(x-mean)/(sigma*sigma);
  double p = (x+mean)*(x+mean)/(sigma*sigma);
  double psi = exp(-m/2) + exp(-p/2);
  return -0.5 * ( (m-1)*exp(-m/2)  + (p-1)*exp(-p/2) ) / (sigma*sigma) / psi;
}

void Measure()
{
  double v, k;

// contribution to energies
  v = Potential(x);
  k = Kinetic(x);

//update walker
  walker[iv] = v;
  walker[ik] = k;
  walker[ie] = k+v;

  //Print instant values
  if (print_instant==1){
    ofstream Ekin, Epot, Etot;

    Epot.open("Temp/Instant/insta_epot.txt",ios::app);
    Ekin.open("Temp/Instant/insta_ekin.txt",ios::app);
    Etot.open("Temp/Instant/insta_etot.txt",ios::app);

    Epot << walker[iv] << endl; //Potential energy
    Ekin << walker[ik] << endl; //Kinetic energy
    Etot << walker[ie] << endl; //Total energy

    Epot.close();
    Ekin.close();
    Etot.close();
  }
  
  // update histogram (only with optimized parameters)
  if(nvar == 0) Histogram();
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

void Histogram(){
    
  // histogram filling
  
  int k=1;
  while((histo_start + k*bin) < x){
    k++;
  }
  walker[2+k]+=1;

    
}

void FinalizeHistogram(int iblk){
    
    for(int i=0; i<nbins;i++){
        
        stima[3+i] = blk_av[3+i]/blk_norm;
        glob_av[3+i] += stima[3+i];
        glob_av2[3+i] += stima[3+i]*stima[3+i];

    // stampa su file (solo ultimo blocco)
    
        if(iblk == nblk){
    
            ofstream Histo;
            Histo.open("Temp/Average/histo.txt",ios::app);
            const int wd=12;

            err[3+i]=Error(glob_av[3+i],glob_av2[3+i],nblk);
            
            Histo << setw(wd) << glob_av[3+i] << setw(wd) << err[3+i] << endl;
            
            Histo.close();
        }
    }
}

void Averages(int iblk) //Print results for current block
{
   ofstream Epot, Ekin, Etot;
   const int wd=12;


    stima[iv] = blk_av[iv]/blk_norm; //Potential energy
    glob_av[iv] += stima[iv];
    glob_av2[iv] += stima[iv]*stima[iv];
    err[iv]=Error(glob_av[iv],glob_av2[iv],iblk);

    stima[ik]  = blk_av[ik]/blk_norm; //Potential energy
    glob_av[ik] += stima[ik];
    glob_av2[ik] += stima[ik]*stima[ik];
    err[ik]=Error(glob_av[ik],glob_av2[ik],iblk);

    stima[ie] = blk_av[ie]/blk_norm; //Potential energy
    glob_av[ie] += stima[ie];
    glob_av2[ie] += stima[ie]*stima[ie];
    err[ie]=Error(glob_av[ie],glob_av2[ie],iblk);


    if (print_blocks) {
      Epot.open("Temp/Average/epot.txt", ios::app);
      Ekin.open("Temp/Average/ekin.txt", ios::app);
      Etot.open("Temp/Average/etot.txt", ios::app);

//Potential energy per particle
      Epot << setw(wd) << stima[iv] << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err[iv] << endl;
//Kinetic
      Ekin << setw(wd) << stima[ik] << setw(wd) << glob_av[ik]/(double)iblk << setw(wd) << err[ik] << endl;
//Total
      Etot << setw(wd) << stima[ie] << setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err[ie] << endl;
 
//Histogram        
      if(nvar == 0) FinalizeHistogram(iblk);

      cout << "Block number " << iblk << endl;
      cout << "Acceptance rate " << double(accepted)/double(attempted) << endl << endl;
      cout << "----------------------------" << endl << endl;

      Epot.close();
      Ekin.close();
      Etot.close();
    }

}

void ConfXYZ(){ //Write configuration in .xyz format
  ofstream WriteX;

  WriteX.open("config.final", ios::app);
  WriteX << x << endl;

  WriteX.close();
}

void Final(int ivar,int iparam)
{
  const int wd = 15;
  energy_new = glob_av[ie]/(double)nblk;

  if (ivar >0){
    cout << "Step number: " << ivar << endl;
    cout << "Acceptance rate mean: " << double(acceptedmean_var)/double(attemptedmean_var) << endl;
    cout << "Acceptance rate sigma: " << double(acceptedsigma_var)/double(attemptedsigma_var) << endl;
    cout << "Delta_mean: " << deltamean_var << endl;
    cout << "Delta_sigma: " << deltasigma_var << endl << endl;

    cout << "Final values: " << endl;
    cout << setw(wd) << "mean = " << setw(wd) << mean << endl;
    cout << setw(wd) << "sigma = " << setw(wd) << sigma << endl;
    cout << setw(wd) << "<H> = " << setw(wd) << energy_new  <<" +/- " << err[ie] << endl;
  }
  ofstream Var;

    
  Var.open("Temp/Paths/path_"+to_string(iparam+1)+".txt", ios::app);

  if (ivar == 0){
    energy_old = energy_new;
    energy_opt = energy_new;
    mean_opt = mean;
    sigma_opt = sigma;
    err_opt = err[ie];
  }

    if ((energy_new < energy_old)) {
        cout << "<H> decreased: new parameters accepted!" << endl;
        energy_old = energy_new;
        
        if(mean_set == true){
            //deltamean_var /= 1.15;
            acceptedmean_var += 1.0;
            
        } else {
            //deltasigma_var /= 1.25;
            acceptedsigma_var += 1.0;
        }
        
        energy_opt = energy_new;
        mean_opt = mean;
        sigma_opt = sigma;
        err_opt = err[ie];
    }
    if (energy_new > energy_old){
        cout << "<H> increased: new parameters refused!" << endl;
        mean = mean_old;
        sigma = sigma_old;
    }

    
  Var << setw(wd) << acceptedmean_var << setw(wd) << mean_opt <<  setw(wd)<< acceptedsigma_var << setw(wd) << sigma_opt << setw(wd) << energy_opt << setw(wd) << err_opt << endl;
  Var.close();    
    
  cout << endl;
  cout << "Optimal values: " << endl;
  cout << setw(wd) << "mean = " << setw(wd) << mean_opt << endl;
  cout << setw(wd) <<"sigma = " << setw(wd) << sigma_opt << endl;
  cout << setw(wd) << "<H> = " << setw(wd) << energy_opt << " +/- " << err_opt << endl;

  cout << "----------------------------" << endl << endl;
   
  if(ivar==nvar){
      acceptedmean_var = 0;
      attemptedmean_var = 0;
      acceptedsigma_var = 0;
      attemptedsigma_var = 0;
  }

}

void Variation(){

  mean_old = mean;
  sigma_old = sigma;
  //int sign = 1;
  //if (rnd.Rannyu() < 0.5) sign = -1;
  if (rnd.Rannyu() < 0.5) {
      mean += deltamean_var * rnd.Rannyu(-1,1);//*sign;
      mean_set = true;
      sigma_set = false;
  }
  else {
      sigma += deltasigma_var * rnd.Rannyu(-1,1);//*sign;
      sigma_set = true;
      mean_set = false;
  } 
      

  attemptedmean_var += 1;
  attemptedsigma_var += 1;

}
    
void SetParam(int i_param){
    
    if(nvar>0){
        ifstream file;
        string line;
        file.open("start_param.txt");
        for(int i=0;i<=i_param;i++){
            getline(file,line);
            mean = stod(line.substr(0,4));
            sigma = stod(line.substr(5,4));
        }
        file.close();
    }
    
    //Evaluate potential energy and virial of the initial configuration
    Measure();

    //Print initial values for the potential energy and virial
    cout << "Starting position = " << x << endl;
    cout << "Initial total energy  = " << walker[ie] << endl << endl;
    cout << "Wave function parameters: " << endl << "mean = " << mean << " , sigma= " << sigma << endl << endl;
    
    acceptedmean_var = 0;
    attemptedmean_var = 0;
    acceptedsigma_var = 0;
    attemptedsigma_var = 0;
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}




/****************************************************************
*****************************************************************
*****************************************************************
*****************************************************************/
