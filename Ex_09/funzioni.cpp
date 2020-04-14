#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib> 
#include <iomanip>
#include "vector"
#include "funzioni.h"

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

Random :: ~Random(){}

void Random :: SaveSeed(){
   ofstream WriteSeed;
   WriteSeed.open("seed.in");
   if (WriteSeed.is_open()){
       WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
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


  ////////////////////////////////////////////
 /////////// Metodi classe CityMap //////////
////////////////////////////////////////////

CityMap :: CityMap(int ncity, bool geom){
    citynum = ncity;
    geometry = geom;
}

CityMap :: ~CityMap(){} 

// crea mappa della città
void CityMap :: GenerateCityMap(Random & rand){
    
// genera città sul cerchio
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
 // genera città nel quadrato
    else{
       cout << "Creating map with cities inside a square [0,1]x[0,1]" << endl;    
	   for (int i = 0; i < citynum; i++) {
		  x.push_back( rand.Rannyu());
		  y.push_back( rand.Rannyu());
	   }
    }

    return;    
}

// carica città da file
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
    cout << "Check Function Error : incorrect individual length!" << endl;
    return 0;
 }

 for (int i = 0; i < dim_ID-1; i++){
    x = Seq.at(i);
	for (int j = i+1; j < dim_ID; j++){
 	  if(x==Seq.at(j)) {
        cout << "Check Function Error : individual with equal entries!" << endl;
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
    
cout << "ID (city sequence): " ;
 for(int i = 0; i < dim_ID; i++) {
	 cout << ID.at(i) << "-";
 }
 cout << endl;
 return;
}

void Individual :: ShowFitness(){
    
    cout << "fitness: " << SelfFitness() << endl;
}

void Individual :: PrintInd(int file){
    
 ofstream print;
 if(geometry) print.open("Temp/Circle/CirclePath."+to_string(file)+".final");
 else print.open("Temp/Square/SquarePath."+to_string(file)+".final");
 
 for(int i = 0; i < dim_ID; i++) {	
     	   print << ID.at(i) <<  endl;
 	}
 
 print.close();

return;
}

double Individual :: Fitness(vector<int> v){
    
    double L = 0;
    for(int i = 0; i < v.size()-1; i++) L += pow(x.at(v.at(i)) - x.at(v.at(i+1)), 2) + pow(y.at(v.at(i)) - y.at(v.at(i+1)),2);

    return L + pow(x.at(v.at(0)) - x.at(v.at(v.size()-1)), 2) + pow(y.at(v.at(0)) - y.at(v.at(v.size()-1)),2);
}

double Individual :: SelfFitness(){

   return Fitness(ID);
}


  ////////////////////////////////////////////
 ////////// Metodi classe GAmethods//////////
////////////////////////////////////////////

GAmethods :: GAmethods(int init_shuffle, double pselection, double pmutation, double pcrossover){

n_shuffle = init_shuffle;
sel_exponent = pselection;
p_mutation = pmutation;
p_crossover = pcrossover;

}

GAmethods :: ~GAmethods(){};

void GAmethods :: PairPermutation(Individual & ind){
    
 int rand1, rand2, temp;
 int range = ind.GetCityNum();
 rand1 = (int) Rannyu(0, range);
 rand2 = (int) Rannyu(0, range);
 while(rand1 == rand2) rand2 = (int) Rannyu(0, range);

 temp = ind.GetElement(rand1);
 ind.SetElement(ind.GetElement(rand2), rand1);
 ind.SetElement(temp, rand2); 

return;
}

void GAmethods :: Shuffle(Individual & ind){

 for(int i = 0; i < n_shuffle; i++){
 PairPermutation(ind);
 }

return;
}

void GAmethods :: Shift(Individual & ind){

 int dim = ind.GetCityNum();
 int K_shift = (int)Rannyu(1, dim);
 vector<int> temp(dim);
 
 for (int i = 0; i < dim;i++){
 temp[i] = ind.GetElement(i - K_shift);
 }
 
 for (int i = 0; i < dim;i++){
 ind.SetElement(temp[i], i);
 }
 
return;
}

void GAmethods :: GroupPermutation(Individual & ind){

 int dim = ind.GetCityNum();
 int M_perm = (int) Rannyu(2, (int) dim/2.0);
 int max_ind = (int) dim/(double)M_perm;
 int first_blk = (int) Rannyu(0, max_ind);
 int sec_blk = (int) Rannyu(0, max_ind);
 while(sec_blk == first_blk) sec_blk = (int) Rannyu(0, max_ind);
 
 vector<int> temp1(M_perm);
 vector<int> temp2(M_perm);

 for (int i = 0; i < M_perm; i++){
 temp1[i] = ind.GetElement(i + first_blk*M_perm);
 temp2[i] = ind.GetElement(i + sec_blk*M_perm);
 }
 
 for (int i = 0; i < M_perm; i++){
 ind.SetElement(temp2[i], i+first_blk*M_perm) ;
 ind.SetElement(temp1[i], i+sec_blk*M_perm) ;
 }
 
return;
}

void GAmethods :: Inversion(Individual & ind){

 int dim = ind.GetCityNum();
 int M_inv = (int) Rannyu(2, dim-1);
 int start_ind = (int) Rannyu(0, dim);
 vector<int> temp(M_inv);

 for (int i = 0; i < M_inv; i++){
 temp[i] = ind.GetElement(i + start_ind);
 }
 
 for (int i = 0; i < M_inv;i++){
 ind.SetElement(temp[M_inv - i - 1], i+start_ind);
 }
 
return;
}

void GAmethods :: Mutation(Individual & ind){

 double r = Rannyu();

 if(r < p_mutation) PairPermutation(ind);
 	
 if(r > p_mutation && r < 2*p_mutation) Shift(ind);
	
 if(r > 2*p_mutation && r < 3*p_mutation) GroupPermutation(ind);
	
 if(r > 3*p_mutation && r < 4*p_mutation) Inversion(ind);
 
 return;
}

void GAmethods :: CrossOver(Individual & Father, Individual & Mother){

vector<int> Son1, Son2;
int dim = Mother.GetCityNum();
int Cut_ind = (int) Rannyu(1, dim-2);		// punto in cui effettuare il taglio
for(int i = 0; i < Cut_ind; i++){
	Son1.push_back(Mother.GetElement(i));	// copia nei figli la prima parte del DNA 
	Son2.push_back(Father.GetElement(i));
 }   
 vector<int> FCutted(dim - Cut_ind);
 vector<int> FPosition(dim - Cut_ind);
 vector<int> MCutted(dim - Cut_ind);
 vector<int> MPosition(dim - Cut_ind);
for(int i= 0; i< dim - Cut_ind; i++){
	FCutted[i] = Father.GetElement(i+Cut_ind);	// Salva i pezzi tagliati in nuovi vettori
	MCutted[i] = Mother.GetElement(i+Cut_ind);
 }
for (int i = 0; i < dim - Cut_ind; i++){				// "Se l'elemento j-esimo del genitore1 è
for (int j = 0; j < dim ; j++){					//  uguale all'elemento i-esimo della parte tagliata
if(Mother.GetElement(j) == FCutted[i]) MPosition[i] = j;//  al genitore2, allora salva la posizione j all'i-esimo posto
if(Father.GetElement(j) == MCutted[i]) FPosition[i] = j;//  di un nuovo vettore *Position[]
	}
 }
sort(FPosition.begin(),FPosition.end());	
sort(MPosition.begin(),MPosition.end());
for(int i = 0; i < dim - Cut_ind; i++){
	Son1.push_back(Father.GetElement(FPosition[i]));
	Son2.push_back(Mother.GetElement(MPosition[i]));
 }   
 Mother.SetNewID(Son1);
 Father.SetNewID(Son2);

return ; 
}


  ////////////////////////////////////////////
 ////////// Metodi classe Population/////////
////////////////////////////////////////////

Population :: Population(int dim_pop,int n_città,bool shape){
    
    n_individui = dim_pop;
    geometry = shape;
    
    Individual ind(n_città,shape);
    ind.ShowInd();
    
    for(int i=0;i<n_individui;i++){
        
        generazione.push_back(ind);
    
    }

}

void Population :: GeneratePopulation(GAmethods & GA){
    
    for(int i=0; i< n_individui; i++){
        
        GA.Shuffle(generazione.at(i));
        if(!generazione.at(i).SelfCheck()) cout << "problem during population generation!" << endl;
    }  
}

void Population :: SortPopulation(){    
    
 if(generazione.size() == 1) {
    cout << "ERROR: Cannot order a size=1 array" << endl;
    return;
 }
 vector<double> OriginalFitness, OrderedFitness;
 for(int i = 0; i< generazione.size(); i++) {
 	OriginalFitness.push_back(generazione.at(i).SelfFitness());
	OrderedFitness.push_back(OriginalFitness.at(i));
 }

sort(OrderedFitness.begin(), OrderedFitness.end());
vector<Individual> CopyPopulation;

for(int i = 0; i < OriginalFitness.size(); i++){
    for(int k = 0; k < OriginalFitness.size(); k++){
        if(OriginalFitness.at(k) == OrderedFitness.at(i)){
			CopyPopulation.push_back(generazione.at(k));
            break;
        }		
    }
}
// Set the old disordered population equal to the sorted one
generazione = CopyPopulation;
return;   
}

void Population :: ShowInd(int i){
    
    cout << "Individuo " << i << " della popolazione: " << endl;
    generazione.at(i).ShowInd();
    
}

double Population :: GetFitnessInd(int i){
    
    return generazione.at(i).SelfFitness();
}

void Population :: CheckPopulation(){
    
    for(int i=0;i<n_individui;i++){
        if(!generazione.at(i).SelfCheck()) cout << "Individual " << i << " : CORRUPTED! ";
    }
}

// da usare dopo la funzione SortPopulation() !!
void Population :: EvolvePopulation(GAmethods & GA){
    
    int dim = generazione.size(), count = 0;
    vector<Individual> New_Population;	// Nuova popolazione
    int pos1, pos2;
    
    while(New_Population.size() < dim){ 	// ripetere finchè la dimensione della nuova popolazione è uguale a quella vecchia
        pos1 = int(pow(GA.Rannyu(), GA.GetExp())*dim);	// Selezione a roulette truccata che favorisce
        pos2 = int(pow(GA.Rannyu(), GA.GetExp())*dim); 	// i primi elementi della popolazione (già ordinata in questa fase)
        
        if(GA.Rannyu()< GA.GetP_cross()){
            New_Population.push_back(generazione.at(dim - pos1 - 1)); // Creazione a 2 a 2 della nuova popolazione
            New_Population.push_back(generazione.at(dim - pos2 - 1));
            for (int i = count; i < count + 2; i++) GA.Mutation(New_Population.at(i)); // Mutazione
            GA.CrossOver(New_Population.at(count), New_Population.at(count+1)); // Crossover
            
            count += 2;
        }
    }
    generazione = New_Population;	// Aggiorna la popolazione
} 

// da usare dopo la funzione SortPopulation() !!
void Population :: Print(int file){
    
    generazione.at(0).PrintInd(file);
    
    AveFitness();

    ofstream output;
    if(geometry) output.open("Temp/Circle/fitness.dat",ios::app);
    else output.open("Temp/Square/fitness.dat",ios::app);
    // stampa la fitness mediata sulla prima metà della popolazione e la fitness del miglior individuo
    output << ave_fitness << " " << generazione.at(0).SelfFitness() << endl;
    output.close();
    
}

// da usare dopo la funzione SortPopulation() !!
void Population :: AveFitness(){
    
    ave_fitness = 0;
    
    
    for(int i=0;i<(n_individui/2);i++){
         ave_fitness += generazione.at(i).SelfFitness();
    }
    ave_fitness /= (double)(n_individui/2);
}

// da usare dopo la funzione SortPopulation() !!
void Population :: ShowAveFitness(){
    
    cout << "Average fitness: " << ave_fitness << endl;
}

// da usare dopo la funzione SortPopulation() !!
void Population :: ShowBest(){
    
    cout << "Best individual in the population: " << endl;
    generazione.at(0).ShowInd();
    generazione.at(0).ShowFitness();
    
}


