#include <vector>

using namespace std;

const double PI = 3.14159265358979323846;


//carica parametri da file di input parameters.dat
void Input(); 

//Generatore di numeri casuali
class Random {

private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

public:

  // constructors
  Random();
  //Random(int *);
  // destructor
  ~Random();
  // methods

  void SetRandom(int * , int, int);
  void SaveSeed();
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

void GenerateCityMap(Random &);     // genera la mappa delle città
void LoadCityMap();                 // carica da file la mappa delle città
void PrintCityMap();                // stampa a video la mappa delle città
void SaveCityMap();                 // salva su file la mappa delle città

int GetCityNum();                   // restituisce il numero di città

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
  void SetNewID(vector<int>);           // cambia ID
  double GetFitness(){return fitness;}; // restituisce la fitness
  void ShowInd();                       // stampa l'ID a video
  void ShowFitness();                   // stampa la fitness a video
  void PrintInd(int);                   // salva l'ID su file
  bool SelfCheck();                     // controlla la correttezza dell'ID
  bool Check(vector<int>);
  double SelfFitness();                 // calcola la fitness
  double Fitness(vector<int>);
    
};

//metodi AG
class GAmethods : public Random{
    
private:
    
    double p_crossover;
    double p_mutation;
    double sel_exponent;
    int n_shuffle;

public:
    
    GAmethods(int,double,double,double);
    ~GAmethods();
    
    void PairPermutation(Individual &);         // operatori di mutazione genetica
    void Shuffle(Individual &);
    void Shift(Individual &);
    void GroupPermutation(Individual &);
    void Inversion(Individual &);
    void Mutation(Individual &);                // chiama una sola delle 4 mutazioni genetiche con una certa probabilità
    void CrossOver(Individual &,Individual &);  // effettua il cross-over tra due genitori
    double GetP_cross(){return p_crossover;};   // restituisce la probabilità di cross-over
    double GetExp(){return sel_exponent;};      // restituisce l'esponente per la selezione truccata dei genitori
    double GetP_mutation(){return p_mutation;}; // restituisce la probabilità di subire una delle 4 possibii mutazioni
 };


//Popolazione
class Population {
    
private:
    int n_individui;
    bool geometry;
    double ave_fitness;
    vector<Individual> generazione;
    //vector<Individual> nuova_generazione;

public:
    Population(int,int,bool);
    ~Population(){};
    
    void GeneratePopulation(GAmethods &);     // genera la popolazione
    void SortPopulation();                  // riordina la popolazione in base alla fitness dell'individuo
    void AveFitness();                      // calcola la fitness mediata sulla metà migliore della popolazione
    void ShowBest();                        // stampa a video l'individuo migliore e la sua fitness
    void ShowAveFitness();                  // stampa a video la fitness mediata sulla metà migliore della popolazione
    void Print(int);                        // stampa su file l'individuo migliore e la fitness media
    void EvolvePopulation(GAmethods &);       // fai evolvere la popolazione con i metodi in GAmethods
    void CheckPopulation();                 // controlla che l'ID di ciascun individuo sia corretto
    void ShowInd(int);                      // stampa a video l'individuo selezionato
    double GetFitnessInd(int);              // restituisce la fitness dell'individuo selezionato
    
    
};









