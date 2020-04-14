#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <stdlib.h> 
//#include "vector"
#include "funzioni.h"
#include "costanti.h"

using namespace std;

int main(){
 
 // Carica parametri da file di input parameters.dat
 Input();
 
 // Crea oggetto SimulatedAnnealing
 SimulatedAnnealing SA(init_shuffle,pmutation,beta0,deltabeta);    
 
 // Crea un individuo (configurazione iniziale) a partire dalla mappa delle città (se non esiste, la crea)
 Individual individuo(ncities,geometry);
 SA.Shuffle(individuo);    
 // Stampa a video la configurazione iniziale e la sua fitness
 individuo.ShowInd();
 individuo.ShowFitness();
 
 // Simulated annealing per 1000 "temperature" diverse da T=inf a T<<1
 for(int i=0; i<ntemp; i++){
     
     // Imposta la temperatura
     SA.SetBeta(i);
     cout << "----------------------------------------------------" << endl << endl;
     cout << "Starting Simulated Annealing with inverse temperature beta: " << SA.GetBeta() << " (T = " << 1/SA.GetBeta() << ") " << endl;
     
     // Esegui nstep passi di Metropolis
     for(int j=0; j<nstep; j++){ SA.Metropolis(individuo); }
     
     // Stampa a video la configurazione finale e la sua fitness alla temperatura i-esima
     individuo.ShowInd();
     individuo.ShowFitness();
     
     // Salva su file Circle(Square)/Circle(Square)Path.i.final la configurazione finale alla temperatura i-esima
     // Salva su file Circle(Square)/fitness.dat la temperatura i-esima e la fitness della conf finale
     SA.Print(i,individuo);
     
     // Stampa a video i rate di accettazione del Metropolis
     cout << "Acceptance rate: " << 100 * SA.GetAccepted()/(double)nstep << " % " << endl << endl;
     
     // Resetta i rate di accettazione
     SA.ResetAccepted();    
 }  
    
 return 0;
}

void Input(){
  
    ifstream input;
    input.open("parameters.dat");
    if(input){
 
        input >> geometry;
        input >> ncities;
        input >> init_shuffle;
        input >> pmutation;
        input >> nstep;
        input >> beta0;
        input >> deltabeta;
        input >> ntemp;
        
        input.close();
    }else cout << "Unable to open file parameters.dat ! ";
    
    input.close();
    
    return;
}

