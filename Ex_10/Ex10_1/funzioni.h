
#include <vector>
using namespace std;
const double PI = 3.14159265358979323846;


//carica parametri da file di input parameters.dat
void Input(); 

//Generatore di numeri casuali
class Random {

private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

protected:

public:
  // constructors
  Random();
  // destructor
  ~Random();
  // methods

  void SetRandom(int * , int, int);
  void SaveSeed();
  double Rannyu(void);
  double Rannyu(double min, double max);
  double Gauss(double mean, double sigma);
};


// creazione mappa della città
class CityMap {

private:

 int citynum;

public:
//costruttore
CityMap(int, bool);
// distruttore
~CityMap();
// se geometry = 1 -> città sul cerchio, altrimenti città nel quadrato
bool geometry;
// Coordinate delle città
vector<double> x;
vector<double> y;

void GenerateCityMap(Random);   // genera la mappa delle città
void LoadCityMap();             // carica da file la mappa delle città
void PrintCityMap();            // stampa a video la mappa delle città
void SaveCityMap();             // salva su file la mappa delle città
int GetCityNum();               // restituisce il numero di città

};

// Individuo
class Individual : public CityMap {
    
private:
    vector<int> ID;
    int dim_ID;
    double fitness;
    
public:
  //Costruttore		
  Individual(int, bool);
  //Distruttore
  ~Individual();
 
  void SetElement(int, int);            // cambia elemento dell'ID
  int GetElement(int);                  // restituisce elemente dell'ID
  void SetNewID(vector<int>);           // cambia ID
  double GetFitness(){return fitness;}; // restituisce la fitness
  void ShowInd();                       // stampa l'ID a video
  void ShowFitness();                   // stampa la fitness a video
  void PrintInd(int);                   // salva l'ID su file
  bool SelfCheck();                     // controlla la correttezza dell'ID
  bool Check(vector<int>);
  double SelfFitness();
  double Fitness(vector<int>);          // calcola la fitness
    
};

// metodi per simulated annealing
class SimulatedAnnealing : public Random {
    
private:
        
    double p_mutation;
    int n_shuffle;
    int accepted = 0;
    double beta;
    double initial_beta;
    double delta_beta;
   

public:
    // contatori PairPermutation,GroupPermutation,Shift,Inversion
    
    int PP = 0;
    int GP = 0;
    int S = 0;
    int I = 0;
        
    SimulatedAnnealing(int,double,double,double);
    ~SimulatedAnnealing(){}
   
    // operatori di mutazione per effettuare mosse
    void PairPermutation(Individual &);
    void Shuffle(Individual &);
    void Shift(Individual &);
    void GroupPermutation(Individual &);
    void Inversion(Individual &);
    void RandomSearch(Individual &);            // chiama una sola delle 4 mutazioni genetiche con una certa probabilità
    
    void SetBeta(int);                          // imposta la temperatura T
    double GetBeta();                           // restituisce l'inverso di T
    double Boltzmann(Individual &,Individual &);// calcola la probabilità di un passo del Metropolis
    void Metropolis(Individual &);              // effettua un passo del Metropolis
    void ResetAccepted();                       // reset del contatore di mosse accettate
    int GetAccepted();                          // restituisce il numero di mosse accettate
    void Print(int,Individual);                 // stampa su file individuo e fitness
    


};