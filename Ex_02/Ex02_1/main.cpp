// Ex. 2.1
// INTEGRALE 1D - SAMPLING CON DISTRIBUZIONE UNIFORME vs IMPORTANCE SAMPLING

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <math.h>
#include <iomanip>

using namespace std;

//  Funzione per calcolo errore
double Error(double AV,double AV2,int n);

int main (){

Random rnd;

// Costanti
int M=100000;	// # punti generati
int N=100;	// # blocchi
int L=int(M/N); // # punti per blocco
    
// Variabili
double x,r; // variabile per numero casuale
double mean_uni, ave1_uni = 0, ave2_uni = 0; // variabili per stima integrale con p uniforme
double mean_imp, ave1_imp = 0, ave2_imp = 0; // variabili per stima integrale con importance sampling (p(x) = 2(1-x)
double error_uni, error_imp; // variabile errore statistico
ofstream file1,file2;

file1.open("Temp/Uniform.txt");
file2.open("Temp/Importance.txt");
    
// Calcolo stima integrale ed errore statistico in funzione del numero di blocchi (p(x)=1 e p(x)=2(1-x))
for(int i=0; i<N; i++){
	mean_uni=0;
    mean_imp=0;
	for(int j=0; j<L; j++){
		x = rnd.Rannyu();
        mean_uni += M_PI*0.5*cos(M_PI*0.5*x)/L;
        r = 1-sqrt(1-x);
        mean_imp += M_PI*0.5*0.5*cos(M_PI*0.5*r)/(1-r)/L;
	}
    ave1_uni += mean_uni;
    ave2_uni += mean_uni*mean_uni;
    
    ave1_imp += mean_imp;
    ave2_imp += mean_imp*mean_imp;
    
    error_uni = Error(ave1_uni/(i+1),ave2_uni/(i+1),i);
    error_imp = Error(ave1_imp/(i+1),ave2_imp/(i+1),i);
    
    // Stampa su file 'Uniform.txt' e 'Importance.txt': stima integrale e stima errore in funzione del n d blocchi
    file1 << fixed << setprecision(7) << ave1_uni/(i+1) << " " << error_uni << endl;
    file2 << fixed << setprecision(7) << ave1_imp/(i+1) << " " << error_imp << endl;
}

file1.close();
file2.close();
         
return 0;
}

// funzione per calcolo errore
double Error(double AV,double AV2,int n){
    
    if(n==0)
        return 0;
    else
        return sqrt((AV2-AV*AV)/n);
}






