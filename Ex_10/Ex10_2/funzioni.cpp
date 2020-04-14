#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib> 
#include <iomanip>
#include "vector"
#include <algorithm>
#include "funzioni.h"
#include <string>

using namespace std;
 

  ////////////////////////////////////////////
 /////////// Metodi classe Random ///////////
////////////////////////////////////////////

Random :: Random(){
 int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();
 
 ifstream input("seed.in");
   string property;
   if (input.is_open()){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            SetRandom(seed,p1,p2);

	input.close();
      }
      
    else cerr << "PROBLEM: Unable to open seed.in" << endl;
return;
}

// costruttore per processi in parallelo: utilizza il rank del processo per ottenere il seme
Random :: Random(int myrank){
 
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();
    
   for(int i=0;i<4;i++) seed[i] = myrank + 2*i;
   SetRandom(seed,p1,p2);
}


Random :: ~Random(){}

void Random :: SaveSeed(){
   ofstream WriteSeed;
   WriteSeed.open("seed.in");
   if (WriteSeed.is_open()){
      WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}

double Random :: Gauss(double mean, double sigma) {
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random :: Rannyu(double min, double max){

   return min+(max-min)*Rannyu();
}

double Random :: Rannyu(void){
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));
//cout << "Numero random: " <<r << endl;
  return r;
}

void Random :: SetRandom(int * s, int p1, int p2){
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0]%4096;
  l2 = s[1]%4096;
  l3 = s[2]%4096;
  l4 = s[3]%4096;
  l4 = 2*(l4/2)+1;
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}

void Random :: SetRank(int myrank){
    
    
}


  ////////////////////////////////////////////
 /////////// Metodi classe CityMap //////////
////////////////////////////////////////////

CityMap :: CityMap(int ncity, bool geom){
    citynum = ncity;
    geometry = geom;
}

CityMap :: ~CityMap(){} 

// Create the city map
void CityMap :: GenerateCityMap(Random rand){
    
// Generate cities on circle
    if(geometry){
	   cout << "Creating map with cities on a unitary circumference" << endl;
       vector<double> theta;

	   for (int i = 0; i < citynum; i++) theta.push_back(rand.Rannyu(0, 2*PI));
	   sort(theta.begin(), theta.end());
	   for (int i = 0; i < citynum; i++) {
		  x.push_back(cos(theta.at(i)));
		  y.push_back(sin(theta.at(i)));
	   }
    }
 // Generate cities inside a square
    else{
       cout << "Creating map with cities inside a square [0,1]x[0,1]" << endl;    
	   for (int i = 0; i < citynum; i++) {
		  x.push_back( rand.Rannyu());
		  y.push_back( rand.Rannyu());
	   }
    }

    return;    
}

// Load cities from file
void CityMap :: LoadCityMap(){
    double a, b;
    ifstream input;
    if(geometry) input.open("circle.start");
    else input.open("square.start");
 
    if(input){
 	  for (int i = 0; i < citynum; i++) {
	   input >> a >> b;
	   x.push_back(a);
	   y.push_back(b);
	   }
    }
    else cout << "Map of cities not found " << endl;
    input.close();

    return;
}

// Save map on file
void CityMap :: SaveCityMap(){

    ofstream print;
    if(geometry) print.open("circle.start");
    else print.open("square.start");
    for (int i = 0; i < citynum; i++) print <<  fixed << setprecision(7) << x[i] << " " << y[i] << endl;
 
    print.close();
}

// Print map on video
void CityMap :: PrintCityMap(){
 
    if(geometry) cout << "Coordinate delle città sulla mappa (città sulla circonferenza): " << endl;
    else cout << "Coordinate delle città sulla mappa (città nel quadrato): " << endl;
 
    cout << "    x    " << "    y    " << endl;       
    for (int i = 0; i < citynum; i++) cout <<  fixed << setprecision(7) << x[i] << " " << y[i] << endl;
}

int CityMap :: GetCityNum(){

    return citynum;
}



  ////////////////////////////////////////////
 ///////// Metodi classe Individual /////////
////////////////////////////////////////////

Individual :: Individual(int ncities, bool geometry) : CityMap(ncities, geometry){

    // caricamento vettori x,y con mappa città da file esistente. Se la mappa non esiste già, crea la mappa 
    Random rand;
    
    ifstream input1,input0;
    input1.open("circle.start");
    input0.open("square.start");
    if((geometry && input1) || (!geometry && input0) ) {
        cout << "the right map already exists!" << endl;
        LoadCityMap();   
    }
    else{
        cout << "The right map doesn't exist!" << endl;
        cout << "Let's create a map..." << endl;
        GenerateCityMap(rand);
        SaveCityMap();
    }
        
    input1.close();
    input0.close();
    
    // creazione individuo (default: [0,1,2,3,...])
    dim_ID = GetCityNum();
    
    for(int i=0; i<dim_ID; i++){
        ID.push_back(i);
    }  
}

Individual :: ~Individual(){}

int Individual :: GetElement(int i){
 
 while(i >= dim_ID || i < 0) {
	if(i >= dim_ID) i = i - dim_ID;
	if(i < 0) i = i + dim_ID;
 }

 return ID.at(i);
}

void Individual :: SetElement(int newid,int i){
 
    while(i >= dim_ID || i < 0) {
	if(i >= dim_ID) i = i - dim_ID;
	if(i < 0) i = i + dim_ID;
    }
    ID.at(i) = newid;

 return;       
}

bool Individual :: Check(vector<int> Seq){
    
 double x;
 
 if (Seq.size() != 30) {
    cout << "Check Function Error : incorrect individual length " << endl;
    return 0;
 }

 for (int i = 0; i < dim_ID-1; i++){
    x = Seq.at(i);
	for (int j = i+1; j < dim_ID; j++){
 	  if(x==Seq.at(j)) {
        cout << "Check Function Error : individual with equal entries! " << endl;
		return 0;
      }
	} 
 }

 return 1;    
}

bool Individual :: SelfCheck(){
    
    return Check(ID);   
}

void Individual :: SetNewID(vector<int> NewID){
    
     if(!Check(NewID)){
        cout << "ERROR: Trying to set corrupted ID in sane Individual." << endl;
     return;
    }
 
    for(int i = 0; i < dim_ID; i++){
        ID.at(i) = NewID.at(i);
    }

    return;    
}

void Individual :: ShowInd(){
    
//Stampa a video il vettore
cout << "Sequenza città: " ;
 for(int i = 0; i < dim_ID; i++) {
	 cout << ID.at(i) << "-";
 }
 cout << endl;
 return;
}

void Individual :: ShowFitness(){
    
    cout << "fitness: " << SelfFitness() << endl;
}

void Individual :: PrintInd(int file,int rank){
 
 ofstream print;
 if(geometry) print.open("Temp/Circle_"+to_string(rank)+"/CirclePath.final");
 else print.open("Temp/Square_"+to_string(rank)+"/SquarePath.final");
 
 for(int i = 0; i < dim_ID; i++) {	
     	   print << ID.at(i) <<  endl;
 	}
 
 print.close();

return;
}

double Individual :: Fitness(vector<int> v){
    
    double L = 0;
    for(int i = 0; i < v.size()-1; i++) {

	   L += pow(x.at(v.at(i)) - x.at(v.at(i+1)), 2) + pow(y.at(v.at(i)) - y.at(v.at(i+1)),2);
    }

    L += pow(x.at(v.at(0)) - x.at(v.at(v.size()-1)), 2) + pow(y.at(v.at(0)) - y.at(v.at(v.size()-1)),2);
    return L;
}

double Individual :: SelfFitness(){

   return Fitness(ID);
}

int Individual :: GetDim(){
    return dim_ID;
}


  ////////////////////////////////////////////
 ////// Metodi classe SimulatedAnnealing/////
////////////////////////////////////////////

SimulatedAnnealing :: SimulatedAnnealing(int init_shuffle, double pmutation, double beta0, double deltabeta,int myrank) : Random(myrank){

n_shuffle = init_shuffle;
p_mutation = pmutation;
initial_beta = beta0;
delta_beta = deltabeta;

}

void SimulatedAnnealing :: PairPermutation(Individual & cavy){
    
 int rand1, rand2, temp;
 int range = cavy.GetCityNum();
 rand1 = (int) Rannyu(0, range);
 rand2 = (int) Rannyu(0, range);
 while(rand1 == rand2) rand2 = (int) Rannyu(0, range);

 temp = cavy.GetElement(rand1);
 cavy.SetElement(cavy.GetElement(rand2), rand1);
 cavy.SetElement(temp, rand2); 

return;
}

void SimulatedAnnealing :: Shuffle(Individual & cavy){

 for(int i = 0; i < n_shuffle; i++){
 PairPermutation(cavy);
 }

return;
}

void SimulatedAnnealing :: Shift(Individual & cavy){

 int dim = cavy.GetCityNum();
 int K_shift = (int)Rannyu(1, dim);
 vector<int> temp(dim);
 
 for (int i = 0; i < dim;i++){
 temp[i] = cavy.GetElement(i - K_shift);
 }
 
 for (int i = 0; i < dim;i++){
 cavy.SetElement(temp[i], i);
 }
 
return;
}

void SimulatedAnnealing :: GroupPermutation(Individual & cavy){

 int dim = cavy.GetCityNum();
 int M_perm = (int) Rannyu(2, (int) dim/2.0);
 int max_ind = (int) dim/(double)M_perm;
 int first_blk = (int) Rannyu(0, max_ind);
 int sec_blk = (int) Rannyu(0, max_ind);
 while(sec_blk == first_blk) sec_blk = (int) Rannyu(0, max_ind);
 
 vector<int> temp1(M_perm);
 vector<int> temp2(M_perm);

 for (int i = 0; i < M_perm; i++){
 temp1[i] = cavy.GetElement(i + first_blk*M_perm);
 temp2[i] = cavy.GetElement(i + sec_blk*M_perm);
 }
 
 for (int i = 0; i < M_perm; i++){
 cavy.SetElement(temp2[i], i+first_blk*M_perm) ;
 cavy.SetElement(temp1[i], i+sec_blk*M_perm) ;
 }
 
return;
}

void SimulatedAnnealing :: Inversion(Individual & cavy){

 int dim = cavy.GetCityNum();
 int M_inv = (int) Rannyu(2, dim-1);
 int start_ind = (int) Rannyu(0, dim);
 vector<int> temp(M_inv);

 for (int i = 0; i < M_inv; i++){
 temp[i] = cavy.GetElement(i + start_ind);
 }
 
 for (int i = 0; i < M_inv;i++){
 cavy.SetElement(temp[M_inv - i - 1], i+start_ind);
 }
 
return;
}

void SimulatedAnnealing :: RandomSearch(Individual & cavy){

 double r = Rannyu();

 if(r < p_mutation) PairPermutation(cavy), PP += 1;
 	
 if(r > p_mutation && r < 2*p_mutation) Shift(cavy) , S += 1;
	
 if(r > 2*p_mutation && r < 3*p_mutation) GroupPermutation(cavy),GP += 1;
	
 if(r > 3*p_mutation && r < 4*p_mutation) Inversion(cavy),I += 1;
 
 return;
}

void SimulatedAnnealing :: SetBeta(int step){
    
    beta = initial_beta + step*delta_beta;
}

double SimulatedAnnealing :: GetBeta(){
    
    return beta;
}

double SimulatedAnnealing :: Boltzmann(Individual & new_ind,Individual & ind){
    
    double exponent = (ind.SelfFitness() - new_ind.SelfFitness());  
    //cout << "(ind.SelfFitness() - new_ind.SelfFitness()): " << exponent << endl;
    double weight = exp(beta*exponent);
    //cout << "exp(exponent): " << weight << endl;
    return weight;
}

void SimulatedAnnealing :: Metropolis(Individual & ind){
    
    Individual new_ind = ind;
    RandomSearch(new_ind);
    
    if(!new_ind.SelfCheck()) return;
    
    double p = Boltzmann(new_ind,ind);
    //cout << "p: " << p << endl;
    double rand = Rannyu();
    //cout << "rand: " << rand << endl;
    if(rand < p ){
        accepted += 1;
        ind = new_ind;
    } 
}

void SimulatedAnnealing :: ResetAccepted(){
    
    accepted = 0;
    PP = 0;
    GP = 0;
    S = 0;
    I = 0;
}

int SimulatedAnnealing :: GetAccepted(){
    
    return accepted;
}

void SimulatedAnnealing :: Print(int mode,int rank,Individual ind){
    
    if(mode==0){ 
        
        ofstream print;
        if(ind.geometry) print.open("Temp/Circle_"+to_string(rank)+"/fitness.dat",ios::app);
        else print.open("Temp/Square_"+to_string(rank)+"/fitness.dat",ios::app);
    
        print.precision(1);
        print << fixed << beta;
        print.precision(5);
        print << setw(10) << ind.SelfFitness() << endl;
    
        print.close();
    }else{
        
        ind.PrintInd(mode,rank);
    }
    
}

void SimulatedAnnealing :: PrintBest(Individual & ind){
    
    ofstream print;
    if(ind.geometry) print.open("Temp/Circle/BestPath.dat");
    else print.open("Temp/Square/BestPath.dat");
    
    for(int i=0; i<ind.GetDim();i++) print << ind.GetElement(i) << endl;
    
    print.close();
    
    if(ind.geometry) print.open("Temp/Circle/BestFitness.dat");
    else print.open("Temp/Square/BestFitness.dat");
    
    print.precision(5);
    print << ind.SelfFitness();
    
    print.close();

}




