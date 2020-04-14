// Es. 1.3 - SIMULAZIONE ESPERIMENTO DELL'AGO DI BUFFON

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <math.h>
#include <iomanip>
#include <cmath>

using namespace std;

// funzione per calcolo errore
double Error(double AV,double AV2,int n);

int main (){

Random rnd;

// costanti
int M=10000000;
int N=100; // numero blocchi
int Nth=M/N; //numero lanci in ciascun blocco
int Nhit; // numero lanci con intersezione
double d=3; //distanza tra le linee
double l=2; //lunghezza ago
  
// variabili
double x,x1,x2;
double inst,ave1 = 0,ave2 = 0;
double error;

ofstream file;
file.open("Temp/Pi.txt");

// calcolo stima pi ed errore statistico
for(int i=0; i<N; i++){
    Nhit=0;
    for(int j=0; j<Nth; j++){
        x=rnd.Rannyu(0,d/2);
        if(x <= l/2){
            x1=rnd.Rannyu(0,1);
            x2=rnd.Rannyu(0,1);
            while((x1*x1 + x2*x2)>1){
                x1=rnd.Rannyu(0,1);
                x2=rnd.Rannyu(0,1);
            }
            double theta = atan2(x2,x1);
            if(sin(theta)*l/2 >= x) Nhit += 1;
        }
    }
    inst = Nth*2*l/(Nhit*d);
    ave1 += inst;
    ave2 += inst*inst;
    error = Error(ave1/(i+1),ave2/(i+1),i);
    
    // stampa su file: stima di pi-greco del blocco i-esimo, stima errore del blocco i-esimo
    file << fixed << setprecision(7) << ave1/(i+1) << " " << error << endl;
    
}

file.close();
           
return 0;}

// funzione per calcolo errore
double Error(double AV,double AV2,int n){
	
	if(n==0)
		return 0;
	else
		return sqrt((AV2-AV*AV)/n);
}
