/* 			Run MultiNest on late type pop


*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <string>
#include <float.h>
#include "multinest.h"
#include "gauss_legendre.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

const double PI = 3.14159265358979;
const double mass = 4.07e6; //Ghez 2008
const double dist = 7960.0; //Ghez 2008
const double G = 6.6726e-8;
const double msun = 1.99e33;
const double sec_in_yr = 3.1557e7;
const double cm_in_au = 1.496e13;
const double cm_in_pc = 3.086e18;
const double km_in_pc = 3.086e13;
const double au_in_pc = 206265.0;
const double asy_to_kms = dist * cm_in_au / (1e5 * sec_in_yr);
const double as_to_km = dist * cm_in_au / (1e5);
const double GM = G * mass* msun;

//char root[100] = "/u/schappell/pmnOld/c_accel_test_pert_more_maxr0.4_-0.2_4_4_0.5_";// root for output files
//define limits of priors
double min_g = -3.0;
double max_g = 2.0;
double min_d = 0.01;
double max_d = 10.0;
double min_a = 0.0;
double max_a = 10.0;
double max_b = 1.0 * cm_in_pc;
//double max_Z = 10.0 * cm_in_pc;
//double max_R = 1.8 * dist * cm_in_au;
//double maxr = sqrt(max_Z*max_Z + max_R*max_R);

int situation;

vector<double> r2d, r2dm;
vector<double> ar;
vector<double> are;
vector<double> pOld, pOldm;
vector<double> min_accel;
vector<double> max_accel;
vector<double> maxz_star;

double gmod, almod, demod, brmod;
double max_Z, max_R, maxr;
//double almod = 4.0;
//double demod = 4.0;
//double brmod = 0.5 * cm_in_pc; //hold alpha, delta, and r_break constant
double use_maxZ, density_norm; //for certain models the max Z needs to be smaller for integration to work
double density_normM, density_box, density_cyl;
int stardex, maserdex;
int num_stars, num_maser;
double max_X,max_Y,min_X,min_Y,innerCut,outerCut,ZprimeM;

double density_int(double Rprime, double Zprime, void* data)
{
  double tmpvalue = pow((sqrt(Rprime*Rprime + Zprime*Zprime)/brmod),demod);
  return fabs(2.0*PI*Rprime*pow((Rprime*Rprime + Zprime*Zprime),(gmod/-2.0))*pow((1.0+tmpvalue),((gmod-almod)/demod)));
}

double density_xyz(double Xprime, double Yprime, void* data)
{
  double rprime = sqrt(pow(Xprime,2.0) + pow(Yprime,2.0) + pow(ZprimeM,2.0));
  double tmpvalue = pow(rprime/brmod,demod);
  return fabs(pow(rprime,(-1.0*gmod))*pow(1.0+tmpvalue,(gmod-almod)/demod));
}

double density_int_box(double ZprimeM, void* data)
{
  return gauss_legendre_2D_cube(40,density_xyz,NULL,min_X,max_X,min_Y,max_Y);
}

int main(int argc, char *argv[])
{
  //Model info
  string tmp = argv[1];
  char root[100];
  strcpy(root,tmp.c_str());
  max_Z = atof(argv[5]);
  max_Z *= cm_in_pc;
  max_R = atof(argv[6]);
  max_R *= cm_in_pc;
  maxr = sqrt(max_Z*max_Z + max_R*max_R);
  situation = stoi(argv[7]);

  almod = atof(argv[2]);
  demod = atof(argv[3]);
  brmod = atof(argv[4]);
  brmod *= cm_in_pc;

  innerCut = atof(argv[8]) * cm_in_pc;
  outerCut = atof(argv[9]) * cm_in_pc;
  min_X = -1.0 * outerCut;
  min_Y = -1.0 * outerCut;
  max_X = 1.0 * outerCut;
  max_Y = 1.0 * outerCut;

  double tmpval = gauss_legendre(40,density_int_box,NULL,0.0,max_Z);
  cout << tmpval << "\n";
  tmpval = gauss_legendre_2D_cube(40,density_int,NULL,0.0,innerCut,0.0,max_Z);
  cout << tmpval << "\n";

}
