/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <iomanip>
#include "MolDyn.h"
#include "random.h"

using namespace std;

int main(){ 

  Input(); 
  Equilibrate();    
   
  nconf = 1;    
  int frame = 1;
  for(int iblk=1; iblk <= nblk; ++iblk){
     
    Reset(iblk);
    for(int istep=1; istep <= block_steps; ++istep){
        Move();           //Move particles with Verlet algorithm
        Measure();
        Accumulate();
        if(((iblk-1) * block_steps + istep)%iprint == 0) {
            cout << "Number of time-steps: " << (iblk-1) * block_steps + istep << endl;
            cout << "The instant value of temperature is " << insta_temp << endl;
        }
        if(istep%10 == 0){
                 //Properties measurement
            ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
            frame += 1;
        }
         nconf += 1;
    }
    Averages(iblk);
  } 
  ConfFinal();         //Write final configuration to restart

  return 0;
}

void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> nblk;
  ReadInput >> iprint;
  ReadInput >> nbins;

  block_steps = nstep/nblk; // number of steps in each block  

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in each block = " << block_steps << endl;
  cout << "Number of bins for the g(r) = " << nbins << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  ip = 4;
  n_props = 5; //Number of observables
    
  //measurement of g(r)
  igofr = 5;
  n_props = n_props + nbins;
  bin_size = (box/2.0)/(double)nbins;
    
  Setting(); // ask 1)whether to equilibrate the system or to start the simulation 2)whether to start from a previous configuration of from fcc configuration

  if(mode2) return;
  else {
    
//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

//Prepare initial velocities
   Random rnd;
      
   cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   double sumv[3] = {0.0, 0.0, 0.0};
   for (int i=0; i<npart; ++i){
     vx[i] = rnd.Gauss(0,sqrt(temp));
     vy[i] = rnd.Gauss(0,sqrt(temp));
     vz[i] = rnd.Gauss(0,sqrt(temp));

     sumv[0] += vx[i];
     sumv[1] += vy[i];
     sumv[2] += vz[i];
   }
   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;

   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = Pbc(x[i] - vx[i] * delta);
     yold[i] = Pbc(y[i] - vy[i] * delta);
     zold[i] = Pbc(z[i] - vz[i] * delta);
   }
  }
    
  return;
}

void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
  
  double v, t, p, vij, pij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp, Pres;
    
  //variables for g(r)
  int step;
  bool control;

  //reset the walker of g(r)
  for (int i=igofr; i<igofr+nbins; ++i) walker[i]=0.0;

  Epot.open("Temp/Solid/Instant/instant_epot.txt",ios::app);
  Ekin.open("Temp/Solid/Instant/instant_ekin.txt",ios::app);
  Temp.open("Temp/Solid/Instant/instant_temp.txt",ios::app);
  Etot.open("Temp/Solid/Instant/instant_etot.txt",ios::app);
  Pres.open("Temp/Solid/Instant/instant_pres.txt",ios::app);

  v = 0.0; //reset observables
  t = 0.0;
  p = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( x[i] - x[j] );
     dy = Pbc( y[i] - y[j] );
     dz = Pbc( z[i] - z[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
       pij = 16.0/pow(dr,12) - 8.0/pow(dr,6);
         
       v += vij; // Potential energy
       p += pij/vol; // Pressure
     }
        
//update g(r) histogram    
    step = 0;
    control = false;
        
    while(control == false && step < nbins){
        if(dr >= step*bin_size && dr < (step+1)*bin_size){
             walker[igofr+step] += 2;
             control = true;
        }else step += 1;   
     }    
        
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    insta_pot = v/(double)npart; //Potential energy per particle
    insta_kin = t/(double)npart; //Kinetic energy per particle
    insta_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    insta_etot = (t+v)/(double)npart; //Total energy per particle
    insta_pres = p + rho * insta_temp; // Pressure
    
    walker[iv] = insta_pot;
    walker[ik] = insta_kin;
    walker[it] = insta_temp;
    walker[ie] = insta_etot;
    walker[ip] = insta_pres;
    
    if(nconf%10 == 0){ // Print instant values every 10 steps

    Epot << insta_pot  << endl;
    Ekin << insta_kin  << endl;
    Temp << insta_temp << endl;
    Etot << insta_etot << endl;
    Pres << insta_pres << endl;

    }
    
    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    Pres.close();

    return;
}

void Reset(int iblk){ // Reset block averages
    
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
    
}

void Accumulate(void){ // Update block averages
    
    for(int i=0; i<n_props; ++i){
        blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}

void Averages(int iblk){ // Print block averages 
    
    ofstream Epot, Ekin, Etot, Temp, Pres, Gofr, Gave;
    const int wd=12;
    
    cout << "Block number " << iblk << endl << endl;
    
    Epot.open("Temp/Solid/Average/epot.txt",ios::app);
    Ekin.open("Temp/Solid/Average/ekin.txt",ios::app);
    Temp.open("Temp/Solid/Average/temp.txt",ios::app);
    Etot.open("Temp/Solid/Average/etot.txt",ios::app);
    Pres.open("Temp/Solid/Average/pres.txt",ios::app);
    Gofr.open("Temp/Solid/gdr/gofr.txt",ios::app);
    Gave.open("Temp/Solid/gdr/gave.txt",ios::app);
    
    stima_pot = blk_av[iv]/blk_norm; //Potential energy per particle
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot*stima_pot;
    err_pot=Error(glob_av[iv],glob_av2[iv],iblk);
    
    stima_kin = blk_av[ik]/blk_norm; //Kinetic energy per particle
    glob_av[ik] += stima_kin;
    glob_av2[ik] += stima_kin*stima_kin;
    err_kin=Error(glob_av[ik],glob_av2[ik],iblk);
    
    stima_temp = blk_av[it]/blk_norm; //Temperature
    glob_av[it] += stima_temp;
    glob_av2[it] += stima_temp*stima_temp;
    err_temp=Error(glob_av[it],glob_av2[it],iblk);
    
    stima_etot = blk_av[ie]/blk_norm; //Total energy per particle
    glob_av[ie] += stima_etot;
    glob_av2[ie] += stima_etot*stima_etot;
    err_etot=Error(glob_av[ie],glob_av2[ie],iblk);
    
    stima_pres = blk_av[ip]/blk_norm; // Pressure
    glob_av[ip] += stima_pres;
    glob_av2[ip] += stima_pres*stima_pres;
    err_pres=Error(glob_av[ip],glob_av2[ip],iblk);

//Potential energy per particle
    Epot << iblk <<  setw(wd) << stima_pot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_pot << endl;
//Kinetic energy per particle
    Ekin << iblk <<  setw(wd) << stima_kin << setw(wd) << glob_av[ik]/(double)iblk << setw(wd) << err_kin << endl;
//Temperature
    Temp << iblk <<  setw(wd) << stima_temp << setw(wd) << glob_av[it]/(double)iblk << setw(wd) << err_temp << endl;
//Total energy per particle
    Etot << iblk <<  setw(wd) << stima_etot << setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err_etot << endl;
//Pressure
    Pres << iblk <<  setw(wd) << stima_pres << setw(wd) << glob_av[ip]/(double)iblk << setw(wd) << err_pres << endl;

//g(r)
    Gofr << setw(wd) << iblk << endl;
    
    for(int i=0; i < nbins; i++){
        delta_V = 4 * pi * ( pow((i * bin_size + bin_size),3) - pow((i * bin_size),3) )/3.;
        stima_gdir = blk_av[igofr + i]/blk_norm/(double)npart/rho/delta_V;
        glob_av[igofr + i] += stima_gdir;
        glob_av2[igofr + i] += stima_gdir*stima_gdir;
        err_gdir = Error(glob_av[igofr + i],glob_av2[igofr + i],iblk);
        
        Gofr << setw(wd) << i << setw(wd) << glob_av[igofr + i]/(double)iblk << setw(wd) << err_gdir << endl;
        
        if(iblk == nblk){
            Gave << setw(wd) << i << setw(wd) << glob_av[igofr + i]/(double)iblk << setw(wd) << err_gdir << endl;    
        }    
    }
    

    cout << "----------------------------" << endl << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    Pres.close();
    Gofr.close();
    Gave.close();
}

void ConfFinal(void){ //Write final and old configuration
  ofstream WriteConf_final, WriteConf_old;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf_final.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf_final << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
    
  cout << "Print old configuration to file config.old " << endl << endl;
  WriteConf_old.open("config.old");

  for (int i=0; i<npart; ++i){
    WriteConf_old << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }  
    
  WriteConf_final.close();  
  WriteConf_old.close();
  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("Temp/Solid/frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

void Equilibrate(void){

    double t;
    char a;
    //ofstream Temp;
    
    //Temp.open("Solid/Instant/instant_temp.txt");
    
    while(mode1){
        
        for(int i = 1; i <= n_eq; i++){
            Move();
            t = GetTemp();
            //Temp << t << endl;
        }
        
        cout << "The temperature, after " << n_eq << " steps is " << t << endl << endl;
        cout << "Do you want to start the simulation with this temperature? [Y/N]. If N the program will go on with the equilibration rescaling velocities according to the desired temperature. " << endl << endl;
        
        cin >> a;
        
        while( a != 'Y' && a != 'N'){
            cout << endl << "Invalid input. Choose Y or N" << endl;
            cin >> a;
        }
        if (a == 'Y') {
            cout << "****** SIMULATION (temperature " << t << ") ****** "<< endl << endl; 
            mode1 = 0;
        }
        else mode1= 1;
        
        if(mode1){
            ConfFinal();
            Rescale();     
        } else{
            //Temp.close();
            return;
        }
    }
    //Temp.close();
    return;
}

void Rescale(void){
  
  double fs;
  double v[3] = {0.0, 0.0, 0.0};
  double sumv[3] = {0.0, 0.0, 0.0};
  double sumv2 = 0.0;  
    
  Move(); // move with Verlet algorithm

  // compute the scale factor fs
  for (int i=0; i<npart; ++i){
  sumv[0] += vx[i];
  sumv[1] += vy[i];
  sumv[2] += vz[i];
  }

  for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   
  for (int i=0; i<npart; ++i){
     v[0] = vx[i] - sumv[0];
     v[1] = vy[i] - sumv[1];
     v[2] = vz[i] - sumv[2];

     sumv2 += v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
   }
   sumv2 /= (double)npart;

   fs = sqrt(3 * temp / sumv2); 

   for (int i=0; i<npart; ++i){  
     vx[i] *= fs; // rescale the velocities
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = Pbc(x[i] - vx[i] * delta); // compute "new" old config. according to the rescaled velocities
     yold[i] = Pbc(y[i] - vy[i] * delta);
     zold[i] = Pbc(z[i] - vz[i] * delta);
   } 

    return;
}

void StartingConf(void){
    
  char a;
  ifstream input; 
    
  cout << "Do you want to start from a previous configuration? [Y/N]. If N (or if a previous configuration doesn't exist) the program will use a fcc configuration" << endl << endl;
    
  cin >> a;
    
  while( a != 'Y' && a != 'N'){
  	cout << endl << "Invalid input. Choose Y or N" << endl;
  	cin >> a;
    cout << endl << endl;
  }
    
  if (a == 'Y') {
  	input.open("config.final");
	if(input){
		cout << "Loading configurations..." << endl;
		for (int i=0; i<npart; ++i){
    			input >> x[i] >> y[i] >> z[i];
			x[i] *= box;
			y[i] *= box;
			z[i] *= box;
  		}
        mode2 = 1;
	}
	
	else {
		
		cout << "ERROR: file config.final doesn't exist! The program will use fcc configuration" << endl;
        mode2 = 0;
	}
      
    input.close();
    
    if(mode2){
	   input.open("config.old");
	   if(input){
		  cout << "Loading old configurations..."<< endl;
          for (int i=0; i<npart; ++i){
            input >> xold[i] >> yold[i] >> zold[i];
			xold[i] *= box;
			yold[i] *= box;
			zold[i] *= box;
  		    }
	   }
	
	   else {
		  cout << "ERROR: file config.old doesn't exist! The program will use fcc configuration" << endl;
		  mode2 = 0;	
	   }
        
       cout << "Do you want to rescale the velocities in order for them to match the desired temperature?" << endl;
   	   cin >> a;
   	   while( a != 'Y' && a != 'N'){
           cout << endl << "Invalid input. Choose Y or N" << endl;
  		   cin >> a;
   	   }
   	
	   if(a == 'Y') Rescale();
            
    } 
	    
} else { mode2 = 0; }
    
return;    
    
}

void Setting(void){
  
  char a;
    
  cout << "Do you want to equilibrate the system? [Y/N]" << endl;
  cin >> a;
  cout << endl;
    
  while(a != 'Y' && a != 'N' ){
      cout << "Invalid input. Choose Y or N" << endl;
      cin >> a;
      cout << endl << endl;
  }  
    
  if(a == 'Y'){
    mode1 = 1; // equilibration
    cout << "Choose the number of steps to equilibrate the system: " << endl;
    cin >> n_eq;
    cout << endl;
    cout << "****** EQUILIBRATION ******" << endl << endl;
    
  } else{
      mode1 = 0; // simulation
      cout << "****** SIMULATION ******" << endl << endl;
  } 
    
  StartingConf();
 
    return;
}

double GetTemp(void){
    
  double t = 0;
  
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
 
  return (2.0 / 3.0) * t/(double)npart;
      
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
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
