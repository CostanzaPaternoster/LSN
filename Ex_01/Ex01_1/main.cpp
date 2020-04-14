
#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <math.h>
#include <iomanip>

using namespace std;

// Funzione per il calcolo dell'errore
double Error(double AV,double AV2,int n);

// Funzione Chi-quadro
double Chi(int counts[],int, int);

// Costanti
int M=100000;	// # punti generati
int N=100;	// # blocchi
int L=int(M/N); // # punti per blocco

int S = 100; // # sottointervalli di [0,1)
int I = 10000; // # punti generati per test Chi-quadro

int main (){

Random rnd;

// Es. 1.1.1 e 1.1.2 - Calcolo integrali con metodo MC (campionamento distribuzione uniforme in [0,1) con generatore di numeri pseudocasuali della classe Random)    
    
// Variabili     
double x; // variabile per numero casuale
double mean, ave1_mean = 0, ave2_mean = 0; // variabili per stima media <r>
double var, ave1_var = 0, ave2_var = 0; // variabili per stima varianza <(r-0.5)^2>
double error_mean, error_var; // variabili per errori finali
ofstream out_mean,out_var;
    
out_mean.open("Temp/mean.out");
out_var.open("Temp/var.out");

// Calcolo stima <r> ed errore su stima <r> in funzione del numero di blocchi
// Calcolo stima <(r-0.5)^2> ed errore su stima <(r-0.5)^2> in funzione del numero di blocchi
    
for(int i=0; i<N; i++){
	mean=0;
    var=0;
	
	for(int j=0; j<L; j++){
        x = rnd.Rannyu();
        mean += x/L;
        var += (x-0.5)*(x-0.5)/L;
	}
	ave1_mean += mean;
	ave2_mean += mean*mean;
    error_mean = Error(ave1_mean/(i+1),ave2_mean/(i+1),i);
    
    ave1_var += var;
    ave2_var += var*var;
    error_var = Error(ave1_var/(i+1),ave2_var/(i+1),i);
    
    // Stampa su file 'mean.out' e 'var.out': ave1/(i+1) (media cumulata), error (errore su media cumulata)
    out_mean << fixed << setprecision(6) << ave1_mean/(i+1) << " " << error_mean << endl;
    out_var << fixed << setprecision(6) << ave1_var/(i+1) << " " << error_var << endl;
    
}
    
out_mean.close();
out_var.close();


// Es. 1.1.3 - Test statistico di Pearson
    
// Variabili     
int counts[S];
double chi, ave1_chi = 0, ave2_chi = 0;
double error_chi;
ofstream out_chi;
    
out_chi.open("Temp/chi.out");    
 
// Calcolo stima Chi-quadro ed errore su stima Chi-quadro in funzione del numero di blocchi    
for(int i=0;i<N;i++){
    for(int l=0;l<S;l++) counts[l] = 0;
    for(int j=0;j<I;j++){
        x = rnd.Rannyu();
        for(int k=0;k<S;k++){
            if(x>=(double)k/(double)S && x<(double)(k+1)/(double)S) counts[k]++;
        }
    }
    chi = Chi(counts,S,I);
    ave1_chi += chi;
    ave2_chi += chi*chi;
    error_chi = Error(ave1_chi/(i+1),ave2_chi/(i+1),i);
    
    // Stampa su file 'chi.out': ave1/(i+1) (media cumulata), error (errore su media cumulata)
    out_chi << fixed << setprecision(6) << ave1_chi/(i+1) << " " << error_chi << endl;    
}

out_chi.close();

return 0;
}


// funzione per calcolo errore
double Error(double AV,double AV2,int n){
	
	if(n==0)
		return 0;
	else
		return sqrt((AV2-AV*AV)/n);
}

// funzione Chi-quadro
double Chi(int counts[],int intervals, int throws){

    double chi2=0;
    for(int i=0; i<intervals; i++){
        chi2 += double(counts[i] - throws/intervals)*double(counts[i] - throws/intervals)/(throws/intervals);
    }
    return chi2;
}







