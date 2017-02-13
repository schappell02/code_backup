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
double max_b = 15.0 * dist * cm_in_au;
//double max_Z = 10.0 * cm_in_pc;
//double max_R = 1.8 * dist * cm_in_au;
//double maxr = sqrt(max_Z*max_Z + max_R*max_R);

int situation;

vector<double> r2dm;
vector<double> pOldm;

double gmod, almod, demod, brmod;
double max_Z, max_R, maxr, min_R;
//double almod = 4.0;
//double demod = 4.0;
//double brmod = 0.5 * cm_in_pc; //hold alpha, delta, and r_break constant
double use_maxZ, density_norm; //for certain models the max Z needs to be smaller for integration to work
//double density_normM, density_box, density_cyl;
int maserdex;
int num_maser;
double innerCut,outerCut;

double density_int(double Rprime, double Zprime, void* data)
{
  double tmpvalue = pow((sqrt(Rprime*Rprime + Zprime*Zprime)/brmod),demod);
  return fabs(2.0*PI*Rprime*pow((Rprime*Rprime + Zprime*Zprime),(gmod/-2.0))*pow((1.0+tmpvalue),((gmod-almod)/demod)));
}

double star_likeZmaser(double z0mod, void* data)
{
  //cout << z0mod << "\n";
  double tmp_value = pow((sqrt(r2dm[maserdex]*r2dm[maserdex]+z0mod*z0mod)/brmod),demod);
  double like_den = pow((1.0+tmp_value),((gmod-almod)/demod));
  if(isinf(like_den)==1)
    {
      cout << "It is power of gamma - alpha / delta that's INF";
      cout << "\n";
      cout << "Gamma is " << gmod << ", alpha is " << almod << ", and delta is " << demod;
      cout << "\n";
    }
  like_den *= pow((r2dm[maserdex]*r2dm[maserdex]+z0mod*z0mod),(gmod/-2.0));
  double like_return = like_den / density_norm;
  if(like_return==0.0){like_return=1e-323;}
  if(like_return != like_return){like_return=1e-323;}
  if(isinf(like_return)==1)
    {
      cout << "Got oo in summation, maser";
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
  if (situation < 2.0)
    {
      almod = Cube[1] * (max_a - min_a) + min_a;
      demod = Cube[2] * (max_d - min_d) + min_d;
      brmod = Cube[3] * cm_in_pc;
    }
  //demod = Cube[1] * (max_d - min_d) + min_d;
  //cout << "Density norm start" << "\n";
  //if (situation > 2)
  //{
      //density_box = gauss_legendre(40,density_int_box,NULL,0.0,max_Z);
      //if (density_box <= 0.0)
      //{
      //use_maxZ = 5.0 * brmod * pow(fabs(gmod/almod),(1.0/demod));
      //density_box = gauss_legendre(40,density_int_box,NULL,0.0,use_maxZ);
      //cout << "Using 5 times Z_rho,max as int limits MASER" << "\n";
      //if (density_box <= 0.0)
      //{
      //use_maxZ = 3.0 * brmod * pow(fabs(gmod/almod),(1.0/demod));
      //density_box = gauss_legendre(40,density_int_box,NULL,0.0,use_maxZ);
      //cout << "Using 3 times Z_rho,max as int limits MASER" << "\n";
      //if (density_box <= 0.0)
      //{
      //use_maxZ = 1.5 * brmod * pow(fabs(gmod/almod),(1.0/demod));
      //density_box = gauss_legendre(40,density_int_box,NULL,0.0,use_maxZ);
      //cout << "Using 1.5 times Z_rho,max as int limits" << "\n";
      //if (density_box <= 0.0)
      //{
      //cout << "Still have a problem with density integration MASER outer";
      //cout <<"\n";
      //cout <<"Gamma is "<<gmod<<", alpha is "<<almod<<", delta is "<<demod<<", and break r is "<<brmod;
      //cout <<"\n";
      //}
      //}
      //}
      //}
      //density_cyl = gauss_legendre_2D_cube(40,density_int,NULL,0.0,innerCut,0.0,max_Z);
      //if (density_cyl <= 0.0)
      //{
      //use_maxZ = 5.0 * brmod * pow(fabs(gmod/almod),(1.0/demod));
      //density_cyl = gauss_legendre_2D_cube(100,density_int,NULL,0,max_R,0,use_maxZ);
      //cout << "Using 5 times Z_rho,max as int limits MASER" << "\n";
      //if (density_cyl <= 0.0)
      //{
      //use_maxZ = 3.0 * brmod * pow(fabs(gmod/almod),(1.0/demod));
      //density_cyl = gauss_legendre_2D_cube(100,density_int,NULL,0,max_R,0,use_maxZ);
      //cout << "Using 3 times Z_rho,max as int limits MASER" << "\n";
      //if (density_cyl <= 0.0)
      //{
      //use_maxZ = 1.5 * brmod * pow(fabs(gmod/almod),(1.0/demod));
      //density_cyl = gauss_legendre_2D_cube(100,density_int,NULL,0,max_R,0,use_maxZ);
      //cout << "Using 1.5 times Z_rho,max as int limits" << "\n";
      //if (density_cyl <= 0.0)
      //{
      //cout << "Still have a problem with density integration MASER inner";
      //cout <<"\n";
      //cout <<"Gamma is "<<gmod<<", alpha is "<<almod<<", delta is "<<demod<<", and break r is "<<brmod;
      //cout <<"\n";
      //}
      //}
      //}
      //}
      //density_normM = density_box - density_cyl;
      //if (density_normM <= 0.0){cout<<"Density norm is less than zero for maser stars"<<"\n";}
      //}
  density_norm = gauss_legendre_2D_cube(40,density_int,NULL,min_R,max_R,0,max_Z);
  if (density_norm <= 0.0)
    {
      use_maxZ = 5.0 * brmod * pow(fabs(gmod/almod),(1.0/demod));
      density_norm = gauss_legendre_2D_cube(100,density_int,NULL,min_R,max_R,0,use_maxZ);
      cout << "Using 5 times Z_rho,max as int limits" << "\n";
      if (density_norm <= 0.0)
	{
	  use_maxZ = 3.0 * brmod * pow(fabs(gmod/almod),(1.0/demod));
	  density_norm = gauss_legendre_2D_cube(100,density_int,NULL,min_R,max_R,0,use_maxZ);
	  cout << "Using 3 times Z_rho,max as int limits" << "\n";
	  if (density_norm <= 0.0)
	    {
	      use_maxZ = 1.5 * brmod * pow(fabs(gmod/almod),(1.0/demod));
	      density_norm = gauss_legendre_2D_cube(100,density_int,NULL,min_R,max_R,0,use_maxZ);
	      cout << "Using 1.5 times Z_rho,max as int limits" << "\n";
	      if (density_norm <= 0.0)
	      {
		cout << "Still have a problem with density integration";
		cout <<"\n";
		cout <<"Gamma is "<<gmod<<", alpha is "<<almod<<", delta is "<<demod<<", and break r is "<<brmod;
		cout <<"\n";
	      }
	    }
	}
    }
  //cout <<"Gamma is "<<gmod<<", alpha is "<<almod<<", delta is "<<demod<<", and break r is "<<brmod<<"\n";
  //cout <<"Density norm = " << density_norm << "\n";
  //cout <<"Max R = " << max_R << " and Max Z = " << max_Z << "\n";
  //cout << "Density norm end" << "\n";

  for(int m=0; m<num_maser; m++)
    {
      maserdex = m;
	  //cout << ar[i] << " " << are[i] << " " << pOld[i] << "\n";
	  //cout << "Begin starlike" << "\n";
      starlike = gauss_legendre(100,star_likeZmaser,NULL,0,max_Z);
	  //cout << starlike << "\n";
      if (starlike == 0.0){starlike=1e-323;}
      total_lnL += pOldm[m]*log(starlike);
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
  max_Z = atof(argv[5]);
  max_Z *= cm_in_pc;
  min_R = atof(argv[7]);
  max_R = atof(argv[8]);
  min_R *= dist * cm_in_au;
  max_R *= dist * cm_in_au;
  maxr = sqrt(max_Z*max_Z + max_R*max_R);
  situation = stoi(argv[6]);
  if (situation > 1.0)
    {
      cout << "Fixing alpha, delta, and r_break"<<"\n";
      almod = atof(argv[2]);
      demod = atof(argv[3]);
      brmod = atof(argv[4]);
      brmod *= cm_in_pc;
    }

      //Read in values for stars in maser field
  cout <<"Reading in stars in maser fields"<<"\n";
  ifstream in("onlySchodel_mn.dat",ios::in);
  double tmp1,tmp2,tmp3,tmp4,tmpvalue;
  if (in.is_open())
    {
      while (!in.eof())
	{
	  in >> tmp1 >> tmp2;
	  r2dm.push_back(tmp1);
	  pOldm.push_back(tmp2);
	}
    }
  in.close();
  num_maser = r2dm.size();

  std::string pfile;
  pfile.append(root);
  pfile.append("priors.txt");

  //Writing priors to file
  ofstream priors;
  priors.open(pfile);
  priors << "#Gamma priors:\n" << min_g << " " << max_g << "\n";
  //cout << "situation equals " << situation << "\n";
  //cout << "Remainder of situation " << fabs(remainder(situation,2)) << "\n";
  if (situation < 2.0)
    {
      priors << "#Alpha priors:\n" << min_a << " " << max_a << "\n";
      priors << "#Delta priors:\n" << 0.0 << " " << max_d << "\n";
      priors << "#Break r priors (pc):\n" << 0.0 << " " << max_b/cm_in_pc << "\n";
    }
  priors.close();

	
  // set the MultiNest sampling parameters
	
  int mmodal = 0;// do mode separation?
	
  int ceff = 0;// run in constant efficiency mode?
	
  int nlive = 1000;// number of live points
	
  double efr = 0.8;// set the required efficiency
	
  double tol = 0.5;// tol, defines the stopping criteria
	
  int ndims = 4;// dimensionality (no. of free parameters)
	
  int nPar = 4;// total no. of parameters including free & derived parameters
	
  int nClsPar = 4;// no. of parameters to do mode separation on

  if (situation > 1.0)
    {
      cout <<"Fixing alpha, delta, and r_break"<< "\n";
      ndims = 1;
      nPar = 1;
      nClsPar = 1;
    }
	
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
