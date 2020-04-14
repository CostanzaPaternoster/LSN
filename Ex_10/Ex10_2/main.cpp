#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <stdlib.h>
#include "funzioni.h"
#include "costanti.h"

#include "mpi.h"

using namespace std;

int main(int argc, char* argv[]){
 
 MPI::Init(argc,argv);

 int size = MPI::COMM_WORLD.Get_size();
 int myrank = MPI::COMM_WORLD.Get_rank();
    
 Input();
 
 SimulatedAnnealing SA(init_shuffle,pmutation,beta0,deltabeta,myrank);    
    
 Individual individuo(ncities,geometry);
 SA.Shuffle(individuo);    
 //individuo.ShowInd();
 SA.Print(0,myrank,individuo);
 
 for(int i=1; i<ntemp; i++){
     
     SA.SetBeta(i);
     
     /*cout << "----------------------------------------------------" << endl << endl;
     cout << "Starting Simulated Annealing with inverse temperature beta: " << SA.GetBeta() << " (T = " << 1/SA.GetBeta() << ") " << endl;*/
     
     for(int j=0; j<nstep; j++){ SA.Metropolis(individuo); }
     
     //individuo.ShowInd();
     SA.Print(0,myrank,individuo);
     
     //cout << "Acceptance rate: " << 100 * SA.GetAccepted()/(double)nstep << " % " << endl << endl;
     
     SA.ResetAccepted();    
 }
 
 individuo.PrintInd(1,myrank); //stampa solo l'ultimo individuo
 
 // CONFRONTO PROCESSI PARALLELI    
 
 int bestrank;
 struct{
     double fit;
     int rank;
 } in,out;
    
 in.fit = individuo.SelfFitness();
 in.rank = myrank;
 cout << "Path's fitness of process " << myrank << " : " << in.fit << endl;
 MPI_Barrier(MPI_COMM_WORLD);
 MPI_Reduce(&in,&out,1,MPI_DOUBLE_INT,MPI_MINLOC,0,MPI::COMM_WORLD);
 
 if(myrank==0) bestrank = out.rank;
 
 MPI_Bcast(&bestrank,1,MPI_INTEGER,0,MPI::COMM_WORLD);
    
 if(myrank == bestrank){
     cout << "The best path was found in process " << myrank << endl;
     individuo.ShowInd();
     SA.PrintBest(individuo);
 }
    
 MPI::Finalize();
 
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

