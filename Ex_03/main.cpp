#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <math.h>
#include <cmath>
#include <iomanip>

using namespace std;

// funzione per calcolo errore
double Error(double AV,double AV2,int n);

int main (){

// Es. 03 - PLAIN VANILLA OPTION PRICING
    
int N=100; // number of blocks
int M=100000; // number of time intervals
int L = M/N;
int Nstep=100; // number of time steps in interval [0,T]
    
double s_0=100; // asset price at t=0
double T=1; // delivery time
double K=100; // strike price
double r=0.1; // risk-free interest rate
double sigma=0.25; // volatility

double s_T; // asset price at time T
double C,P; // call(put)-option price and squared call(put) option price
double aveC=0, aveP=0, aveC2=0, aveP2=0;
double errC, errP;
    
double z;
ofstream file;
    
Random rnd;
  
// Computation via MC of the European call/put option price at time t=0 by SAMPLING THE FINAL ASSET PRICE S(T) for a GBM(r,σ^2)

file.open("Temp/Direct.txt");
    
for(int i=0; i<N; i++){ // cycle on number of blocks
    C=0,P=0;
    for(int j=0; j<L; j++){ //cyle in each block
        z = rnd.Gauss(0,1);
        s_T = s_0*exp((r-sigma*sigma*0.5)*T+sigma*z*sqrt(T)); // computation of S(T) which exhibits a GBM
    
        if(s_T-K >0) C += exp(-r*T)*(s_T-K)/L; // computation of the profit - discounted by a factor exp(-rT)
        else P += exp(-r*T)*(K-s_T)/L;
    }
    
    aveC += C;
    aveC2 += C*C;
    aveP += P;
    aveP2 += P*P;
    
    errC = Error(aveC/(i+1),aveC2/(i+1),i);
    errP = Error(aveP/(i+1),aveP2/(i+1),i);
    
    // streaming on file "Direct.txt" (aveC, err_C, aveP, err_P as a function of number of blocks)
    file << fixed << setprecision(3) <<  aveC/(i+1) << " " << errC << " " << aveP/(i+1) << " " << errP << endl;

}
    
file.close();
    
// Computation via MC of the European call/put option price at time t=0 SAMPLING THE DISCRETIZED GBM(r,σ2) PATH OF THE ASSET PRICE dividing [0,T] in 100 time intervals

aveC=0, aveP=0, aveC2=0, aveP2=0;
    
file.open("Temp/Discrete.txt");

for(int i=0; i<N; i++){ // cycle on number of blocks
    C=0,P=0;
    for(int j=0; j<L; j++){ //cyle in each block
        s_T = s_0;
        
        for(int k=0; k<Nstep; k++){
            z = rnd.Gauss(0,1);
            s_T = s_T*exp((r-sigma*sigma*0.5)*(T/Nstep) + z*sigma*sqrt(T/Nstep)); // computation of discretized value of the asset price
        }
        
        if(s_T-K >0) C += exp(-r*T)*(s_T-K)/L;  // computation of the profit - discounted by a factor exp(-rT)
        else P += exp(-r*T)*(K-s_T)/L;
    }
    
    aveC += C;
    aveC2 += C*C;
    aveP += P;
    aveP2 += P*P;
    
    errC = Error(aveC/(i+1),aveC2/(i+1),i);
    errP = Error(aveP/(i+1),aveP2/(i+1),i);
    
    // streaming on file "Discrete.txt" (aveC, err_C, aveP, err_P as a function of number of blocks)
    file << fixed << setprecision(3) <<  aveC/(i+1) << " " << errC << " " << aveP/(i+1) << " " << errP << endl;
}
    
file.close();
    
return 0;
}

// funzione per calcolo errore
double Error(double AV,double AV2,int n){
	
	if(n==0) return 0;
	else return sqrt((AV2-AV*AV)/n);
}


