/****************************************************************
*****************************************************************
Laboratorio Simulazione Numerica 2019
Random.h

Francesco Ariele Piziali
*****************************************************************
*****************************************************************/

#ifndef __Random__
#define __Random__

class Random {

private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

protected:
  int seed[4], p1, p2;

public:
  // constructors
  Random();
  Random(std::string seed_file);
  // destructor
  ~Random();
  // methods
  void SetRandom(int * , int, int);
  void SaveSeed();
  double Rannyu(void);
  double Rannyu(double min, double max);
  double Gauss(double mean, double sigma);
  double Expo(double rate);
  double CauchyLorentz(double mean, double gamma);
  int RanInt(int min, int max);
};

#endif // __Random__

/****************************************************************
*****************************************************************
*****************************************************************
*****************************************************************/
