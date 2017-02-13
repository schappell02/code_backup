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
//#include "gsl/gsl/gsl_integration.h"
#include <iostream>
#include <iomanip>
#include "boost/config/user.hpp"
//#include <boost/numeric/quadrature/adaptive.hpp>
#include <fstream>
#include <string>
#include <vector>
using namespace std;
//namespace quadrature=boost::numeric::quadrature;

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

//define limits of priors
double min_g = -3.0;
double max_g = 2.0;
double plate_scale = 0.00995;
double dxy = plate_scale * dist * cm_in_au; //dx or dy of a pixel

int nonRadial;

vector<double> r2d;
vector<double> ar;
vector<double> are;
vector<double> pOld;
vector<double> min_accel;
vector<double> max_accel;
vector<double> maxz_star;
vector<int> gcows_rows, gcows_cols;

double gmod;
double max_r, max_R;
double use_maxZ, density_norm;
//for certain models the max Z needs to be smaller for integration to work
double density_box, density_cyl, use_min_accel;
int stardex;
int num_stars, num_gcows;
double max_X,max_Y,min_X,min_Y,Rcut,Xprime,Yprime,Rprime,Rcyl,norm_pos,max_Z;
  

double density_int(double Rprime, double Zprime, void* data)
{
  //cout << "Density integration over R and z" << "\n";
  return fabs(2.0*PI*Rprime*pow((Rprime*Rprime + Zprime*Zprime),(gmod/-2.0)));
}

double density_intZ(double Zprime, void* data)
{
  return fabs(dxy*dxy*pow((Xprime*Xprime + Yprime*Yprime + Zprime*Zprime),(gmod/-2.0)));
}

double density_intZcyl(double Zprime, void* data)
{
  return fabs(2.0*PI*Rcyl*pow((Rcyl*Rcyl + Zprime*Zprime),(gmod/-2.0)));
}


double density_intR(double Rprime, void* data)
{
  double max_Z = sqrt(max_r*max_r - Rprime * Rprime);
  //double min_Z;
  Rcyl = Rprime * 1.0;
  return gauss_legendre(100,density_intZcyl,NULL,0.0,max_Z);
}


double density_sphere(double r3dprime, void* data)
{
  return fabs(4.0*PI*r3dprime*r3dprime*pow(r3dprime,gmod));
}

double pos_int(double ax, void* data){return exp(-1.0*(ar[stardex]-ax)*(ar[stardex]-ax)/(2.0*are[stardex]*are[stardex]));}

double star_likeZ(double z0mod, void* data)
{
  //cout << z0mod << "\n";
  double like_pos;
  if ((ar[stardex]) < 0.0)
    {
      double amod = -1.0*GM*r2d[stardex] / pow((sqrt(r2d[stardex]*r2d[stardex] + z0mod*z0mod)),3.0);
      like_pos = exp(-1.0*(ar[stardex]-amod)*(ar[stardex]-amod)/(2.0*are[stardex]*are[stardex]));
    }
  else
    {
      like_pos = 1.0;
    }


  double like_den = pow((r2d[stardex]*r2d[stardex]+z0mod*z0mod),(gmod/-2.0));
  double like_return = like_pos*like_den;

  if(like_return==0.0){like_return=1e-323;}
  if(like_return != like_return){like_return=1e-323;}
  if(isinf(like_return)==1)
    {
      cout << "Got oo in summation";
      cout << "\n";
    }

  return like_return;
}


void LogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{
  double total_lnL = 0.0;
  double starlike;
  gmod = Cube[0] * (max_g - min_g) + min_g;

  if (nonRadial > 0)
    {
      //cout << "Density norm for GCOWS" << "\n";
      density_norm = 0.0;
      //density_norm_ad = 0.0;
      for(int i=0; i<num_gcows; i++)
	{
	  //cout << "In row loop" << "\n";
	  //for(int col=0; col<gcowsCols; col++)
	  //{
	      //cout << "In column loop" << "\n";
	  Rprime = sqrt((gcows_rows[i]-1500.0)*(gcows_rows[i]-1500.0)+(gcows_cols[i]-1500.0)*(gcows_cols[i]-1500.0))*dxy;
	      //cout << "R of star " << Rprime << " and in gcows " << gcowsField[row][col] << "\n";
	  max_Z = sqrt(max_r*max_r - Rprime*Rprime);
	  //double min_Z;
	  if ((Rprime < Rcut) & (Rprime > 0.0))
	    {
	      //Rcyl = Rprime * 1.0;
	      Xprime = sqrt((gcows_cols[i]-1500.0)*(gcows_cols[i]-1500.0))*dxy;
	      Yprime = sqrt((gcows_rows[i]-1500.0)*(gcows_rows[i]-1500.0))*dxy;
	      
	      double tmpDen = gauss_legendre(100,density_intZ,NULL,0.0,max_Z);

	      //Test other integration

	      //double result,result_error;
	      //quadrature::adaptive().relative_accuracy(1e-5).absolute_accuracy(1e-7)(density_intZ_ad,min_Z,max_Z,result,result_error);
	      //double tmpDen_ad = result * 1.0;

	      //cout << " " << "\n";
	      //cout << "Gauss Legendre: " << tmpDen << "\n";
	      //cout << "Adaptive Int: "<< result << "\n";

	      //cout << tmpDen << "\n";
	      if (tmpDen <= 0.0)
		{
		  cout << "Still have a problem with density integration, i " << i << "\n";
		  cout <<"\n";
		  cout <<"Gamma is "<<gmod;
		  cout <<"\n";
		  tmpDen = 0.0;
		}
	      density_norm += tmpDen;
	      //density_norm_ad += tmpDen_ad;
	    }
	}
    }
  else
    {
      //cout << "Density integration for accelerating sample with R2d" << "\n";
      density_norm = gauss_legendre(100,density_intR,NULL,0,Rcut);
      //double density_hole = gauss_legendre(100,density_sphere,NULL,0,rhmod);
      //density_norm -= density_hole;

      if (density_norm <= 0.0)
	{
	  cout <<"GCOWS part of density norm is 0, make max r larger"<<"\n";
	}
    }

  //cout << "Density norm " << density_norm << "\n";
  double for_print = gauss_legendre(100,density_intR,NULL,0,Rcut);
  //cout << "Test norm " << for_print << "\n";
  double norm_pos,min_Z;
  for(int i=0; i<num_stars; i++)
    {
      stardex = i;

      use_min_accel = min_accel[stardex] * 1.0;
      if ((ar[stardex]) < 0.0)
	{
	  norm_pos = are[stardex]*sqrt(PI/2.0)*(erf((ar[stardex]-use_min_accel)/(sqrt(2.0)*are[stardex]))
						-erf((ar[stardex]-max_accel[stardex])/(sqrt(2.0)*are[stardex])));
	  if(norm_pos==0.0){norm_pos=gauss_legendre(40,pos_int,NULL,use_min_accel,max_accel[stardex]);}
	  norm_pos=fabs(norm_pos);
	  if(norm_pos==0.0)
	    {
	      cout << "Norm pos is zero";
	      cout <<"\n";
	    }
	}
      else
	{
	  norm_pos = fabs(max_accel[stardex] - use_min_accel);
	}
      double max_Z = sqrt(max_r*max_r - r2d[stardex]*r2d[stardex]);
      starlike = gauss_legendre(100,star_likeZ,NULL,0.0,max_Z);
      if (starlike == 0.0){starlike=1e-323;}
      //cout << "pos and density like " << starlike << "\n";
      starlike *= 1.0 / (norm_pos * density_norm);
      //cout << "norm pos " << norm_pos << "\n";
      //cout << "density norm " << density_norm << "\n";
      total_lnL += pOld[i]*log(starlike);
    }

  lnew = total_lnL * 1.0;
  //cout << "Gamma = " << gmod << "\n";
  //cout << "Alpha = " << almod << "\n";
  //cout << "Delta = " << demod << "\n";
  //cout << "R_break = " << brmod << "\n";
  //cout << "Constant = " << cmod << "\n";
  //cout << "R_hole = " << rhmod << "\n";
  //  cout << "Total ln L = " << lnew << "\n";
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
  tmp = argv[5];
  char ending[100];
  strcpy(ending,tmp.c_str());
  max_r = atof(argv[2]);
  max_r *= cm_in_pc;
  Rcut = atof(argv[3]);
  Rcut *= cm_in_pc;
  nonRadial = stoi(argv[4]);


  //Read in values for stars
  std::string gcowsfile;
  gcowsfile.append("stars_mn");
  gcowsfile.append(ending);
  gcowsfile.append(".dat");
  ifstream in(gcowsfile,ios::in);
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
	  tmpvalue = -1.0*GM*tmp1 / (max_r*max_r*max_r);
	  max_accel.push_back(tmpvalue);
	}
    }
  in.close();
  num_stars = r2d.size();

  if (nonRadial > 0)
    {
      //Read in gcows field
      cout << "Reading in GCOWS field info"<<"\n";
      ifstream in("gcows_field.dat",ios::in);
      double tmp1,tmp2;
      if (in.is_open())
	{
	  while(!in.eof())
	    {
	      in >> tmp1 >> tmp2;
	      gcows_rows.push_back((int) tmp1);
	      gcows_cols.push_back((int) tmp2);
	    }
	}

      in.close();
      num_gcows = gcows_rows.size();
    }
  
  std::string pfile;
  pfile.append(root);
  pfile.append("priors.txt");

  //Writing priors to file
  ofstream priors;
  priors.open(pfile);
  priors << "#Gamma priors:\n" << min_g << " " << max_g << "\n";
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
