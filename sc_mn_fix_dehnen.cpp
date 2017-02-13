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

//char root[100] = "/u/schappell/pmnOld/c_accel_test_dehnen_";// root for output files
//define limits of priors
double min_g = -3.0;
double max_g = 2.0;
//double max_Z = 10.0 * cm_in_pc;
//double max_R = 1.8 * dist * cm_in_au;
//double maxr = sqrt(max_Z*max_Z + max_R*max_R);

vector<double> r2d;
vector<double> ar;
vector<double> are;
vector<double> pOld;
vector<double> min_accel;
vector<double> max_accel;
vector<double> maxz_star;

double gmod, scamod, max_Z, max_R, maxr;
//double amod = 0.5 * cm_in_pc; //hold alpha, delta, and r_break constant
double use_maxZ, density_norm; //for certain models the max Z needs to be smaller for integration to work
int stardex;
int num_stars;

double density_int(double Rprime, double Zprime, void* data)
{
  return Rprime*(3.0-gmod)*scamod*pow((Rprime*Rprime+Zprime*Zprime),(gmod/-2.0))*pow((sqrt(Rprime*Rprime+Zprime*Zprime)+scamod),
										   (-4.0+gmod));
}

double pos_int(double ax, void* data){return exp(-1.0*(ar[stardex]-ax)*(ar[stardex]-ax)/(2.0*are[stardex]*are[stardex]));}

double star_likeZ(double z0mod, void* data)
{
  //cout << z0mod << "\n";
  double like_pos, norm_pos;
  if (ar[stardex] < 0.0)
    {
      double amod = -1.0*GM*r2d[stardex] / pow((sqrt(r2d[stardex]*r2d[stardex] + z0mod*z0mod)),3.0);
      like_pos = exp(-1.0*(ar[stardex]-amod)*(ar[stardex]-amod)/(2.0*are[stardex]*are[stardex]));
      norm_pos = are[stardex]*sqrt(PI/2.0)*(erf((ar[stardex]-min_accel[stardex])/(sqrt(2.0)*are[stardex]))
					    -erf((ar[stardex]-max_accel[stardex])/(sqrt(2.0)*are[stardex])));
      if(norm_pos==0.0){norm_pos=gauss_legendre(100,pos_int,NULL,min_accel[stardex],max_accel[stardex]);}
      norm_pos=fabs(norm_pos);
      if(norm_pos==0.0)
	{
	  cout << "Norm pos is zero";
	  cout <<"\n";
	}
    }
  else
    {
      like_pos = pow((r2d[stardex]*r2d[stardex] + z0mod*z0mod),2.5) / (3.0*z0mod*GM*r2d[stardex]);
      norm_pos = maxz_star[stardex];
    }

  double like_den = pow((sqrt(r2d[stardex]*r2d[stardex]+z0mod*z0mod)+scamod),(-4.0+gmod));
  if(isinf(like_den)==1)
    {
      cout << "It is power of 4 - gamma that's INF";
      cout << "\n";
      cout << "Gamma is " << gmod << " and a is " << scamod;
      cout << "\n";
    }
  like_den *= pow((r2d[stardex]*r2d[stardex]+z0mod*z0mod),(gmod/-2.0))*(3.0-gmod)*scamod;
  double like_return = like_pos*like_den/ (norm_pos * density_norm);
  if(like_return==0.0){like_return=1e-323;}
  if(like_return != like_return){like_return=1e-323;}
  if(isinf(like_return)==1)
    {
      cout << "Got oo in summation";
      cout << "\n";
    }
  //cout << "z mod = " << z0mod << "\n";
  //cout << "like pos = " << like_pos << "\n";
  //cout << "like den = " << like_den << "\n";
  //cout << "den norm = " << density_norm << "\n";
  //cout << "pos norm = " << norm_pos << "\n";
  //cout << "tmp value = " << tmp_value << "\n";
  return like_return;
}

void LogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{
  double total_lnL = 0.0;
  double starlike;
  gmod = Cube[0] * (max_g - min_g) + min_g;
  //demod = Cube[1] * (max_d - min_d) + min_d;
  //cout << "Density norm start" << "\n";
  density_norm = gauss_legendre_2D_cube(40,density_int,NULL,0,max_R,0,max_Z);
  if (density_norm <= 0.0)
    {
      cout << "Still have a problem with density integration";
      cout <<"\n";
      cout <<"Gamma is "<<gmod<<" and a is "<<scamod;
      cout <<"\n";
    }
  //cout <<"Gamma is "<<gmod<<" and a is "<<scamod<<"\n";
  //cout <<"Density norm = " << density_norm << "\n";
  //cout <<"Max R = " << max_R << " and Max Z = " << max_Z << "\n";
  //cout << "Density norm end" << "\n";
  for(int i=0; i<num_stars; i++)
    {
      stardex = i;
      //cout << ar[i] << " " << are[i] << " " << pOld[i] << "\n";
      //cout << "Begin starlike" << "\n";
      starlike = gauss_legendre(100,star_likeZ,NULL,0,max_Z);
      //cout << starlike << "\n";
      if (starlike == 0.0){starlike=1e-323;}
      total_lnL += pOld[i]*log(starlike);
      //cout << starlike << "\n";
      //cout << total_lnL << "\n";
      //cout << "End starlike" << "\n";
    }
  lnew = total_lnL * 1.0;
  //cout << "Total ln L = " << lnew << "\n";
}

void dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &logZerr, void *context)
{
	// convert the 2D Fortran arrays to C++ arrays
	
	
	// the posterior distribution
	// postdist will have nPar parameters in the first nPar columns & loglike value & the posterior probability in the last two columns
	
	int i, j;
	
	double postdist[nSamples][nPar + 2];
	for( i = 0; i < nPar + 2; i++ )
		for( j = 0; j < nSamples; j++ )
			postdist[j][i] = posterior[0][i * nSamples + j];
	
	
	
	// last set of live points
	// pLivePts will have nPar parameters in the first nPar columns & loglike value in the last column
	
	double pLivePts[nlive][nPar + 1];
	for( i = 0; i < nPar + 1; i++ )
		for( j = 0; j < nlive; j++ )
			pLivePts[j][i] = physLive[0][i * nlive + j];
}


int main(int argc, char *argv[])
{
  //Model info
  string tmp = argv[1];
  char root[100];
  strcpy(root,tmp.c_str());
  scamod = atof(argv[2]);
  max_Z = atof(argv[3]);
  max_R = atof(argv[4]);
  max_Z *= cm_in_pc;
  max_R *= cm_in_pc;
  maxr = sqrt(max_Z*max_Z + max_R*max_R);

  //Read in values for stars
  ifstream in("stars_mn.dat",ios::in);
  double tmp1,tmp2,tmp3,tmp4,tmpvalue;
  if (in.is_open())
    {
      while (!in.eof())
	{
	  in >> tmp1 >> tmp2 >> tmp3 >> tmp4;
	  r2d.push_back(tmp1);
	  ar.push_back(tmp2);
	  are.push_back(tmp3);
	  pOld.push_back(tmp4);
	  tmpvalue = -1.0*GM / (tmp1*tmp1);
	  min_accel.push_back(tmpvalue);
	  tmpvalue = sqrt(tmp1*tmp1 + max_Z*max_Z);
	  tmpvalue = -1.0*GM*tmp1 / (tmpvalue*tmpvalue*tmpvalue);
	  max_accel.push_back(tmpvalue);
	  tmpvalue = sqrt(maxr*maxr - tmp1*tmp1);
	  maxz_star.push_back(tmpvalue);
	}
    }
  in.close();
  num_stars = r2d.size();

  std::string pfile;
  pfile.append(root);
  pfile.append("priors.txt");

  //Writing priors to file
  ofstream priors;
  priors.open(pfile);
  priors << "#Gamma priors:\n" << min_g << " " << max_g << "\n";
  //priors << "#Delta priors:\n" << 0.0 << " " << max_d << "\n";
  priors.close();

	
  // set the MultiNest sampling parameters
	
  int mmodal = 0;// do mode separation?
	
  int ceff = 0;// run in constant efficiency mode?
	
  int nlive = 1000;// number of live points
	
  double efr = 0.8;// set the required efficiency
	
  double tol = 0.5;// tol, defines the stopping criteria
	
  int ndims = 1;// dimensionality (no. of free parameters)
	
  int nPar = 1;// total no. of parameters including free & derived parameters
	
  int nClsPar = 1;// no. of parameters to do mode separation on
	
  int updInt = 100;// after how many iterations feedback is required & the output files should be updated
                    // note: posterior files are updated & dumper routine is called after every updInt*10 iterations
	
  double Ztol = -1E90;// all the modes with logZ < Ztol are ignored
	
  int maxModes = 100;// expected max no. of modes (used only for memory allocation)
	
  int pWrap[ndims];// which parameters to have periodic boundary conditions?
  for(int i = 0; i < ndims; i++) pWrap[i] = 0;
	
  int seed = -1;// random no. generator seed, if < 0 then take the seed from system clock
	
  int fb = 1;// need feedback on standard output?
	
  int resume = 1;// resume from a previous job?
	
  int outfile = 1;// write output files?
	
  int initMPI = 1;// initialize MPI routines?, relevant only if compiling with MPI
		  // set it to F if you want your main program to handle MPI initialization
	
  double logZero = -1E20;// points with loglike < logZero will be ignored by MultiNest
	
  int maxiter = 0;// max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it 
		  // has done max no. of iterations or convergence criterion (defined through tol) has been satisfied
	
  void *context = 0;// not required by MultiNest, any additional information user wants to pass	
	
  // calling MultiNest
	
  nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI,logZero, maxiter, LogLike, dumper, context);
}
