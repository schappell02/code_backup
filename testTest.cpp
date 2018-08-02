#include "integrator.hpp"
#include "gauss_legendre.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <string>
#include <float.h>
#include <iostream>
#include <iomanip>
//#include <boost/numeric/quadrature/adaptive.hpp>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
//using namespace std;
using std::cin;
using std::cout;
using std::vector;
using std::ofstream;
using std::ifstream;
using std::ios;
using std::string;
using std::stoi;
using std::endl;
//using namespace boost::numeric;


double func(double Xp, void* data)
{
  return Xp*Xp;
}

int main()
{
  double tmp = gauss_legendre(100,func,NULL,0.0,1.0);
  cout << "Gauss " << tmp << endl;
  integrator::QuadIntegral inti(1.0e-6,1000);
  auto funcP = [&](double x){return x*x;};
  double tmpP = inti.RombInt(funcP,0.0,1.0);
  cout << "Romb " << tmpP << endl;
}
