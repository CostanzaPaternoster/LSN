// Lezione 2 - Es. 2.2

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <math.h>
#include <cmath>
#include <iomanip>

using namespace std;

//  Funzione per calcolo errore
double Error(double AV,double AV2,int n);


int main (){

Random rnd;
    
// RANDOM WALK 3D DISCRETA
    
// Costanti
int a=1; //lunghezza passo
int N=100; //numero blocchi
int M=10000; //RW totali
int L=M/N; //RW in un blocco
int Nstep=100; //passi RW

// Variabili
double x[3]; //vettore posizione
double r,r_mean[Nstep+1],r_ave[Nstep+1],r2_ave[Nstep+1];
double err;

int c,p;
double theta,phi,y;
ofstream file1,file2;
    
file1.open("Temp/Discrete_ave.txt");
file2.open("Temp/Discrete_instant.txt");

for(int i=0; i<Nstep+1; i++){ //set a zero delle distanze dall'origine
    r_ave[i]=0;
    r2_ave[i]=0;
}    

for(int i=0; i<N; i++){ //ciclo sui blocchi
    for(int k=0; k<Nstep+1; k++) r_mean[k]=0; //reset media del blocco
    
    for(int j=0; j<L; j++){ //ciclo nel blocco
        for(int k=0; k<3; k++) x[k]=0;//reset del vettore posizione
        
        for(int k=1; k<Nstep+1; k++){ // generazione passo k
            c = int(rnd.Rannyu(0,3)); //scelta coordinata x,y,z con p=1/3
            p = int(rnd.Rannyu(0,2)); //scelta passo avanti o indietro con p=1/2
            x[c] += a*p - a*(1-p); //passo della j-esima RW del blocco
            r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]); // aggiornamento distanza ad ogni passo
            r_mean[k] += r/L; // aggiornamento media distanza nel blocco
            if(i==0 && j<10) file2 << r << endl; // stampa 10 cammini sul file 'Discrete_instant.txt'
        }
    }
    for(int k=1;k<Nstep+1;k++){
        r_ave[k] += r_mean[k]/N; // aggiornamento medie dei blocchi
        r2_ave[k] += r_mean[k]*r_mean[k]/N;
    }
}

// Calcolo errore su distanza media e stampa su file file 'Discrete_r.txt': distanza + errore (blocco N)

for(int k=0;k<Nstep+1;k++){
    err = Error(r_ave[k],r2_ave[k],N-1);
    file1 << fixed << setprecision(7) << r_ave[k] << " " << err << " " << err/r_ave[k] << endl;
}

file1.close();
file2.close();


// RANDOM WALK 3D CONTINUO

// Calcolo errore su r e stampa su file file 'Continuum_r.txt': r + errore (blocco N)
file1.open("Temp/Continuum_ave.txt");
file2.open("Temp/Continuum_instant.txt");

for(int i=0; i<Nstep+1; i++){ //set a zero delle distanze dall'origine
    r_ave[i]=0;
    r2_ave[i]=0;
}

for(int i=0; i<N; i++){ //ciclo sui blocchi
    for(int k=0;k<Nstep+1;k++) r_mean[k]=0;
    
    for(int j=0; j<L; j++){ //reset media del blocco
        for(int k=0;k<3;k++) x[k]=0; //reset del vettore posizione
        
        for(int k=1; k<Nstep+1; k++){ // generazione passo k
            phi = rnd.Rannyu(0,2*M_PI); //scelta angolo polare con p uniforme in [0,2pi)
            y = rnd.Rannyu();
            theta = acos(1-2*y); //scelta angolo azimutale theta in [0,pi) con p=sin(theta)/2
            
            x[0]+=a*sin(theta)*cos(phi); //aggiornamento coordinate al passo k
            x[1]+=a*sin(theta)*sin(phi);
            x[2]+=a*cos(theta);
            
            r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
            r_mean[k] += r/L; //somma delle distanze al passo k
            if(i==0 && j<10) file2 << r << endl; // stampa 10 cammini sul file 'Discrete_instant.txt'
        }
    }
    for(int k=1;k<Nstep+1;k++){
        r_ave[k] += r_mean[k]/N; //distanza dall'origine al passo k mediata sugli N blocchi
        r2_ave[k] += r_mean[k]*r_mean[k]/N;//quadrato della distanza dall'origine al passo k mediata sugli N blocchi
    }
}
    
    
for(int k=0;k<Nstep+1;k++){
    err = Error(r_ave[k],r2_ave[k],N-1);
    file1 << fixed << setprecision(7) << r_ave[k] << " " << err << " " << err/r_ave[k] << endl;
}
    
file1.close();
file2.close();
    
return 0;
}

// funzione per calcolo errore
double Error(double AV,double AV2,int n){
    
    if(n==0) return 0;
    else return sqrt(abs(AV2-AV*AV)/n);
}




