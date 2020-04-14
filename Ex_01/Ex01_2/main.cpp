#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <math.h>
#include <iomanip>

using namespace std;

int main (){

    Random rnd;

    // ESERCIZIO 01.2 - VERIFICA TEOREMA DEL LIMITE CENTRALE

    int M=10000;
    int N[5]={0,1,2,10,100};
    string Data[3]={"Temp/Uniform.txt","Temp/Exp.txt","Temp/Lorentz.txt"};
    double sum[3];

    ofstream file[3];

    for(int i=0; i<3; i++) file[i].open(Data[i]);
    
    for(int i=0; i<M; i++){
        for(int l=0; l<3; l++) sum[l]=0;
        
        for(int k=1; k<5; k++){
            for(int j=0; j<N[k]-N[k-1]; j++){
                sum[0] += rnd.Rannyu();
                sum[1] += rnd.Exp(1);
                sum[2] += rnd.Lorentz(1,0);
            }
            for(int l=0; l<3; l++) file[l] << fixed << sum[l]/N[k] << " ";
        }
        for(int l=0; l<3; l++) file[l] << endl;
    }

    for(int i=0; i<3; i++) file[i].close();
    
    return 0;
}