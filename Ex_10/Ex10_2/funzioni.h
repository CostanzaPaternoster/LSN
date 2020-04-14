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
  Random(int); //costruttore per processi in parallelo: utilizza il rank del processo per ottenere il seme 
  // destructor
  ~Random();
  // methods

  void SetRandom(int * , int, int);
  void SaveSeed();
  void SetRank(int);
  double Rannyu(void);
  double Rannyu(double min, double max);
  double Gauss(double mean, double sigma);
};


//Mappa delle città
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

void GenerateCityMap(Random);   // genera mapa delle città
void LoadCityMap();             // carica mappa delle città  nei vettori x,y
void PrintCityMap();            // stampa a video la mappa delle città
void SaveCityMap();             // salva su file la mappa delle città

int GetCityNum();               // restituisce il numero di città

};

//Individuo
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
  void SetNewID(vector<int>);			// cambia ID
  double GetFitness(){return fitness;}; // restituisce la fitness
  void ShowInd();                       // stampa l'ID a video
  void ShowFitness();                   // stampa la fitness a video
  void PrintInd(int,int);				// salva l'ID su file, in cartella diversa diversa a seconda del rank del processo
  bool SelfCheck();                     // controlla la correttezza dell'ID
  bool Check(vector<int>);
  double SelfFitness();                 // calcola il costo del cammino (fitness dell'individuo)
  double Fitness(vector<int>);
  int GetDim();                         // restituisce la dimensione dell'individuo
    
};


class SimulatedAnnealing : public Random {
    
private:
        
    double p_mutation;
    int n_shuffle;
    int accepted = 0;
    double beta;
    double initial_beta;
    double delta_beta;
   

public:
    
    int PP = 0;
    int GP = 0;
    int S = 0;
    int I = 0;
        
    SimulatedAnnealing(int,double,double,double,int);
    ~SimulatedAnnealing(){SaveSeed();}
    
    // operatori di mutazione per la generazione delle mosse
    void PairPermutation(Individual &);
    void Shuffle(Individual &);
    void Shift(Individual &);
    void GroupPermutation(Individual &);
    void Inversion(Individual &);
    void RandomSearch(Individual &); // seleziona, con uguale probabilità, una tra le 4 mosse possibili
    
    void SetBeta(int);                              // imposta la temperatura T
    double GetBeta();                               // restituisce l'inverso di T
    double Boltzmann(Individual &,Individual &);    // calcola la probabilità di un passo del Metropolis
    void Metropolis(Individual &);                  // effettua un passo del Metropolis
    void ResetAccepted();                           // reset del contatore di mosse accettate
    int GetAccepted();                              // restituisce il numero di mosse accettate
    void Print(int mode,int rank ,Individual);      // stampa su file individuo e sua fitness, in cartelle diverse a seconda del rank del processo. mode = 0 -> stampa solo fitness, mode ≠ 0 -> stampa fitness e individuo
    void PrintBest(Individual &);                   //stampa individuo e fitness: viene utilizzata dopo aver confrontato i risultati dei diversi processi
    


};
