#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <stdlib.h>
#include "funzioni.h"
#include "costanti.h"


using namespace std;

int main(){
    
    // Carica parametri da file di input parameters.dat
    Input();
    
    // Crea la mappa delle città (geometry = 1 -> città su circonferenza, geometry = 0 -> città nel quadrato [0,1]x[0,1])
    /*CityMap mappa(ncities,geometry);
    mappa.GenerateCityMap(rnd); // crea la mappa
    mappa.SaveCityMap();    // Salva la mappa su file
    mappa.PrintCityMap();   // Stampa a video la mappa
    */
     
    // Dichiara oggetto GAmethods che implementa gli operatori di mutazione genetica e crossover
    GAmethods GA(init_shuffle,selection_exp,pmutation,pcrossover);
    
    // Genera la popolazione a partire dalla mappa (se non esiste, la crea)
    Population popolazione(dim_population,ncities,geometry);
    
    popolazione.GeneratePopulation(GA);   // Genera la popolazione
    popolazione.CheckPopulation();  // Verifica che ciascun individuo soddisifi i vincoli
    popolazione.SortPopulation();   // Riordina la popolazione in base alla fitness di ciascun individuo (per selezione) 
    
    popolazione.AveFitness(); // Calcola la fitness mediata sulla metà migliore della popolazione
    popolazione.ShowBest();     //  Stampa a video l'individuo migliore
    popolazione.ShowAveFitness(); // Stampa a video la fitness media
    popolazione.Print(0);   // Stampa su file l'individuo migliore e la fitness media
    
    // Fai evolvere la popolazione nstep volte
    for(int i=1; i<=nstep; i++){
    
        popolazione.EvolvePopulation(GA); // Fai evolvere la popolazione favorendo la selezione di individui con fitness più alta (funzione costo più bassa)
        popolazione.CheckPopulation();
        popolazione.SortPopulation();
        
        cout << "-----------------------" << endl << "Generation n° " << i << endl;
        popolazione.ShowBest();
        popolazione.ShowAveFitness();
        popolazione.Print(i);   // Stampa su file "Circle(Square)/Circle(Square)Path.i.final" l'individuo migliore della generazione i-esima, stampa su file "Circle(Square)/fitness.dat" la fitness media
    }
    
    //GA.SaveSeed();

    
    return 0;
}


void Input(){
  
    ifstream input;
    input.open("parameters.dat");
    if(input){
 
        input >> geometry;
        input >> ncities;
        input >> init_shuffle;
        input >> selection_exp;
        input >> pmutation;
        input >> pcrossover;
        input >> nstep;
        input >> dim_population;
        
        input.close();
    }else cout << "Unable to open file parameters.dat ! ";
    
    input.close();
    
    return;
}

