// EX. 05 - METROPOLIS - SAMPLING ORBITALI ATOMICI

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <math.h>
#include <cmath>
#include <iomanip>

using namespace std;

void Input();
void Equilibrate();
void Trial(int trial_mode);
void Metropolis(int orbital_mode);
double S(double x,double y, double z); // distribuzione di prob. orbitale 1S
double P(double x, double y, double z); // distribuzione di prob. orbitale  2P
bool Test(double p_new, double p_old); //funzione accettazione/rigetto mossa
double Error(double AV,double AV2,int n); // funzione d'errore


Random rnd;



// Costanti
int N=100; //numero blocchi
int M=100000; //passi totali
int L = M/N;
int mode_mossa = 0;
int mode_orbitale = 0;
string filename_mossa;
string filename_orbital;
double delta;
int eq_steps;

// Variabili
double phi,c,p_mossa,p;
bool accept;
double X[2],Y[2],Z[2];
double r, mean_r, ave1_r = 0, ave2_r = 0;
double err_r;
int count_accept = 0;
ofstream file_xyz, file_r;

int main(int argc, char *argv[]){
 
// posizione inziale
X[0]=0;
Y[0]=0;
Z[0]=1;
    
Input();
Equilibrate();
    
file_xyz.open("Temp/"+filename_mossa+"/"+filename_orbital+"_xyz.txt");
file_r.open("Temp/"+filename_mossa+"/"+filename_orbital+"_r.txt");

// stampa su file "Uniform/1S_xyz.txt" la posizione iniziale
file_xyz << fixed << X[0] << " " << Y[0] << " " << Z[0] << endl;
    
for(int i=0;i<N;i++){
    mean_r = 0;
    for(int j=0; j<L; j++){
        
        Trial(mode_mossa);
        
        Metropolis(mode_orbitale);
        
        // stampa su file "Uniform/1S_xyz.txt" le coordinate aggiornate della posizione
        file_xyz << fixed << setprecision(7) << X[0] << " " << Y[0] << " " << Z[0] << endl;
        
        r = sqrt(X[0]*X[0]+Y[0]*Y[0]+Z[0]*Z[0]); // distanza dall'origine
        mean_r += r/L;
    }
    
    // medie cumulative
    ave1_r += mean_r;
    ave2_r += mean_r*mean_r;
        
    // errore nel blocco i-esimo
    err_r = Error(ave1_r/(i+1),ave2_r/(i+1),i);
        
    // stampa su file "Uniform/1S_r.txt": distanza media dall'origine e errore statistico in funzione del blocco
    file_r << fixed << setprecision(7) << ave1_r/(i+1) << " " << err_r << endl;
}
 
file_xyz.close();
file_r.close();
    
cout << "Acceptance rate: " <<  int(100*double(count_accept)/double(M)) << " % " << endl;
    
return 0;

}
    
/************************************  FUNZIONI ************************************/

void Input(){
    
    while(mode_orbitale!=1 && mode_orbitale!=2){
        cout << "Choose orbital to sample: type 1 for 1S or type 2 for 2P" << endl;
        cin >> mode_orbitale;
    }
    
    while(mode_mossa!=1 && mode_mossa!=2){
        cout << "Choose how to make a move: type 1 for uniform movements or type 2 for gaussian movements" << endl;
        cin >> mode_mossa;
    }
    
    cout << "Choose step length. Suggestion: 1.1 for uniform sampling of 1S, 0.8 for gaussian sampling of 1S, 2.5 for uniform sampling of 2P, 1.8 for gaussian sampling of 2P" << endl << endl;
    
    cin >> delta;
    
    if(mode_orbitale==1) filename_orbital = "1S";
    else filename_orbital = "2P";
    
    if(mode_mossa==1) filename_mossa = "Uniform";
    else filename_mossa = "Gaussian";
    
    cout << endl << "Sampling orbital "+filename_orbital+" with Metropolis algorithm having "+filename_mossa+" movements of length " << delta << endl << "Starting position: X = "  << X[0] << " Y = " << Y[0] << " Z = " << Z[0] << endl << endl;
}

// funzione per l'equilibrazione del sistema
void Equilibrate(){
    
    string mode;
    while(mode!="y" && mode!="n"){
        cout << "Do you want to equilibrate the system? [y/n]" << endl;
        cin >> mode;
    }
    if(mode=="y"){
        cout << "Choose the number of steps for equilibration. " << endl;
        cin >> eq_steps;
        ofstream file;
        file.open("Temp/"+filename_mossa+"/"+filename_orbital+"_instant_r.txt");
        cout << endl << "Starting " << eq_steps << " steps of equilibration" << endl << endl;
        for(int i=0; i<eq_steps; i++){
            Trial(mode_mossa);
            Metropolis(mode_orbitale);
            r = sqrt(X[0]*X[0]+Y[0]*Y[0]+Z[0]*Z[0]);
            file << r << endl; // stampa su file distanza istantanea in fase di equilibrazione
        }
    cout << "Acceptance rate: " <<  int(100*double(count_accept)/double(eq_steps)) << " % " << endl;
        cout << endl << "Starting " << M << " step simulation " << endl << endl;
    count_accept = 0;
    file.close();
    } else cout << "Starting "  << M << " step simulation without equilibration" << endl << endl;
}
        
// funzione mossa
void Trial(int trial_mode){
    if(trial_mode==1){
        // tentativo di mossa con campionamento uniforme
        phi=rnd.Rannyu(0,2*M_PI);
        c=rnd.Rannyu(-1,1);
        
        // calcolo nuove posizioni
        X[1]=X[0]+delta*sqrt(1-c*c)*cos(phi);
        Y[1]=Y[0]+delta*sqrt(1-c*c)*sin(phi);
        Z[1]=Z[0]+delta*c;
    } else {
        // tentativo di mossa con campionamento gaussiano
        X[1]=X[0]+rnd.Gauss(0,delta);
        Y[1]=Y[0]+rnd.Gauss(0,delta);
        Z[1]=Z[0]+rnd.Gauss(0,delta);
    }
}
        
// algoritmo Metropolis
void Metropolis(int orbital_mode){
    if(orbital_mode==1){
        // campionamento orbitale 1S
        p_mossa = S(X[1],Y[1],Z[1]); // densità di prob. nel nuovo punto generato
        p = S(X[0],Y[0],Z[0]); // densità di prob. nel punto attuale
    } else {
        // campionamento orbitale 2P
        p_mossa = P(X[1],Y[1],Z[1]); //valutazione densità di prob. nel punto generato
        p = P(X[0],Y[0],Z[0]); // valutazione densità di prob. nel punto attuale
    }
    
    accept = Test(p_mossa,p); //accettazione/respingimento mossa
    
    if(accept){
        X[0]=X[1];
        Y[0]=Y[1];
        Z[0]=Z[1];
        count_accept +=1;
    }
}

// funzione accettazione/rigetto mossa
bool Test(double p_new, double p_old){
    
    double rand = rnd.Rannyu();
    if(p_new < p_old){
        if(p_new/p_old >= rand)return true;
        else return false;
    }else return true;    
}

// distribuzione di probabilità 1S
double S(double x,double y, double z){
    
    double R= sqrt(x*x+y*y+z*z);
    double P = exp(-2*R); // r in unità del raggio di Bohr
    
    return P;
}

// distribuzione di probabilità 2P
double P(double x, double y, double z){
    
    double R= sqrt(x*x+y*y+z*z);
    double cos_theta=z/R;
    double P= R*R*exp(-R)*cos_theta*cos_theta; // r in unità del raggio di Bohr P= R*R*exp(-R)*cos_theta*cos_theta/(32*M_PI);
    
    return P;
    
}

// funzione d'errore
double Error(double AV,double AV2,int n){
	
	if(n==0)
		return 0;
	else
		return sqrt((AV2-AV*AV)/n);
}




