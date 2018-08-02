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
//#include "sc_mn_lib.h"
#include "ez_thread.hpp"
//#include "gsl/gsl/gsl_integration.h"
#include <iostream>
#include <iomanip>
//#include <boost/config/user.hpp>
//#include <boost/math/special_functions/erf.hpp>
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

std::mutex mutex2;
int NThreads = 14;

//define limits of priors
double min_g = -3.0;
double max_g = 2.0;
double min_d = 2.0;
double max_d = 10.0;
double min_a = 0.0;
double max_a = 10.0;
double min_b = 5.0*dist*cm_in_au;
double max_b = 2.0*cm_in_pc;
double meu_b = 10.25 * (4.0/3.0) * dist * cm_in_au;
double sigma_b = 5.75 * dist * cm_in_au;
double plate_scale = 0.00995;
double dxy = plate_scale * dist * cm_in_au; //dx or dy of a pixel

int situation,nonRadial;

vector<double> r2dv, r2dvm, arv,arve, amodv, like_returnv, starlikev, starlikevm, like_returnvm;
vector<double> Xgcows,Ygcows,Rgcows,Zmax_gcows,rho_gcows;
vector<double> pOldv, pOldvm, Zmax_star, Zmax_starm;
vector<double> min_accelv,max_accelv,norm_posv;
vector<int> gcows_vrows, gcows_vcols;

double density_normv, max_r, Rcut, innerCut, outerCut;
double gmodv, almodv, demodv, brmodv, cmodv;
int num_stars, num_gcows, num_maser;



//vector<double> amodv, like_returnv, like_returnvm, Xgcows, Ygcows, rho_gcows, Rgcows, Zmax_gcows;
//vector<double> Zmax_star, Zmax_starm, starlikev, starlikevm, norm_posv, min_accelv, max_accelv;

double density_intZcyl(double Zprime, void* data)
{
  double Rcyl = *(double *)data;
  double tmpvalue = pow((sqrt(Rcyl*Rcyl + Zprime*Zprime)/brmodv),demodv);
  return fabs(2.0*PI*Rcyl*pow((Rcyl*Rcyl + Zprime*Zprime)/(brmodv*brmodv),(gmodv/-2.0))*
	      pow((1.0+tmpvalue),((gmodv-almodv)/demodv)));
}


double density_intR(double Rprime, void* data)
{
  double max_Z = sqrt(max_r*max_r - Rprime * Rprime);
  double Rcyl = Rprime * 1.0; 
  return gauss_legendre(100,density_intZcyl,&Rcyl,0.0,max_Z);
}



double density_intZ(double Zprime, void* data)
{
  int iii = *(int *)data;
  return fabs(dxy*dxy*pow((Xgcows[iii]*Xgcows[iii] + Ygcows[iii]*Ygcows[iii] + Zprime*Zprime)/(brmodv*brmodv),(gmodv/-2.0))*
	      pow((1.0+pow((sqrt(Xgcows[iii]*Xgcows[iii] + Ygcows[iii]*Ygcows[iii] + Zprime*Zprime)/brmodv),demodv)),
		  ((gmodv-almodv)/demodv)));
}

//double density_intZBLAH(double Zprime, void *data)
//{
//  double tmpvalue = pow((sqrt(Xprime*Xprime + Yprime*Yprime + Zprime*Zprime)/brmodv),demodv);
  //cout << "X and Y " << Xprime << " " << Yprime << endl;
  //cout << "tmpvalue " << tmpvalue << endl;
//  return fabs(dxy*dxy*pow((Xprime*Xprime + Yprime*Yprime + Zprime*Zprime)/(brmodv*brmodv),(gmodv/-2.0))*
//	      pow((1.0+tmpvalue),((gmodv-almodv)/demodv)));
//}


double pos_int(double ax, void* data)
{
  int iii = *(int *)data;
  return exp(-1.0*(arv[iii]-ax)*(arv[iii]-ax)/(2.0*arve[iii]*arve[iii]));
}



double star_likeZ(double z0modv, void* data)
{
  int iii = *(int *)data;
  if (arv[iii] < 0.0)
  {
    amodv[iii] = -1.0*GM*r2dv[iii] / pow((sqrt(r2dv[iii]*r2dv[iii] + z0modv*z0modv)),3.0);
    like_returnv[iii] = exp(-1.0*(arv[iii]-amodv[iii])*(arv[iii]-amodv[iii])/(2.0*arve[iii]*arve[iii]));
  }
  else
  {
      like_returnv[iii] = 1.0;
    }

  like_returnv[iii] *= pow((1.0+pow((sqrt(r2dv[iii]*r2dv[iii]+z0modv*z0modv)/brmodv),demodv)),((gmodv-almodv)/demodv));
  like_returnv[iii] *= pow((r2dv[iii]*r2dv[iii]+z0modv*z0modv)/(brmodv*brmodv),(gmodv/-2.0));
  if(like_returnv[iii]==0.0){like_returnv[iii]=1e-323;}
  if(like_returnv[iii] != like_returnv[iii]){like_returnv[iii]=1e-323;}
  if(isinf(like_returnv[iii])==1)
    {
      cout << "Got oo in summation for gcows stars";
      cout << endl;
    }
  return like_returnv[iii];
}


double star_likeZmaser(double z0modv, void* data)
{
  int iii = *(int *)data;
  like_returnvm[iii] = pow((1.0+pow((sqrt(r2dvm[iii]*r2dvm[iii]+z0modv*z0modv)/brmodv),demodv)),((gmodv-almodv)/demodv));
  like_returnvm[iii] *= pow((r2dvm[iii]*r2dvm[iii]+z0modv*z0modv)/(brmodv*brmodv),(gmodv/-2.0));
  if(like_returnvm[iii]==0.0){like_returnvm[iii]=1e-323;}
  if(like_returnvm[iii] != like_returnvm[iii]){like_returnvm[iii]=1e-323;}
  if(isinf(like_returnvm[iii])==1)
    {
      cout << "Got oo in summation for schodel stars";
      cout << endl;
    }
  return like_returnvm[iii];
}  



void LogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{
  lnew = 0.0;
  gmodv = Cube[0] * (max_g - min_g) + min_g;
  if (fabs(remainder(situation,2)) > 0.0)
    {
      almodv = Cube[1] * (max_a - min_a) + min_a;
      demodv = Cube[2] * (max_d - min_d) + min_d;
      brmodv = Cube[3] * (max_b - min_b) + min_b; //flat prior on r_break
      //brmodv = exp((log(max_b)-log(min_b))*Cube[3] + log(min_b)); //log prior on r_break
      //brmodv = sqrt(2.0) * sigma_b * boost::math::erf_inv(2.0*Cube[3] - 1.0) + meu_b; //gaussian prior on r_break


      if (situation == 3){cmodv = Cube[4];}
    }
  density_normv= 0.0;
  //double density_gcowsv = 0.0;
  //double density_schodelv = 0.0;
  double total_lnLv = 0.0;
  //double blahDex = 0.0;

  if (nonRadial > 0)
    {
      ez_thread(threadnum, NThreads)
      {
	for(int i=threadnum; i<num_gcows; i+=NThreads)
	  {
	    int iii = i;
	    //Rprime = sqrt((gcows_vrows[iii]-1500.0)*(gcows_vrows[iii]-1500.0)+(gcows_vcols[iii]-1500.0)*(gcows_vcols[iii]-1500.0))*dxy;
	    //double max_Zv = sqrt(max_r*max_r - Rprime*Rprime);
	    if ((Rgcows[iii] < Rcut) & (Rgcows[iii] >= 0.0))
	      {
		//Xprime = sqrt((gcows_vcols[iii]-1500.0)*(gcows_vcols[iii]-1500.0))*dxy;
		//Yprime = sqrt((gcows_vrows[iii]-1500.0)*(gcows_vrows[iii]-1500.0))*dxy;
		//double tmpDen = gauss_legendre(100,density_intZBLAH,NULL,0.0,max_Zv);
		//double tmpDen, result_error;
		//quadrature::adaptive().relative_accuracy(1e-5).absolute_accuracy(1e-7)(density_intZBLAH,0.0,max_Zv,tmpDen,result_error);
		rho_gcows[iii] = gauss_legendre(100,density_intZ,&iii,0.0,Zmax_gcows[iii]);
		if (rho_gcows[iii] <= 0.0)
		  {
		    cout << "Integration problem with GCOWS density norm" << endl;
		  }
	      //cout << iii << " max Z " << std::setprecision(20) << max_Zv << " and density " << std::setprecision() << tmpDen << endl;
	      //tmpDen = 1e47;
		//density_normv += tmpDen;
		//density_gcowsv += ;
		//mutex2.lock();
		//blahDex += 1.0;
		//density_gcowsv += rho_gcows[iii];
		//mutex2.unlock();
          }
	  }
      };
      for (auto &&elem : rho_gcows)
	density_normv += elem;
    }

  else
    {
      density_normv += gauss_legendre(100,density_intR,NULL,0.0,Rcut);
      if (density_normv <= 0.0){cout << "GCOWS part of density, radial integration, is 0"<<endl;}
    }


  density_normv *= cmodv;
  //density_gcowsv *= cmodv;
  if (situation > 2)
    {
      double tmp = gauss_legendre(100,density_intR,NULL,innerCut,outerCut);
      if (tmp <= 0.0)
	{
	  cout << "Schodel density norm integration is 0, something is wrong" << endl;
	}
      density_normv += (1.0 - cmodv) * tmp;
      //density_schodelv += (1.0 - cmodv) * tmp;
    }

  ez_thread(threadnum,NThreads)
    {
      for(int i=threadnum; i<num_stars; i+=NThreads)
	{
	  int iii = i;
	  starlikev[iii] = gauss_legendre(100,star_likeZ,&iii,0.0,Zmax_star[iii]);
	  //cout << "INDEX " << iii << endl;
	  //cout << "starlikev  " << starlikev[iii] << endl;
	  //cout << "norm pos " << norm_posv[iii] << endl;
	  //cout << "Density norm" << density_normv << endl;
	  //cout << "prob old " << pOldv[iii] << endl;
	  starlikev[iii] *= cmodv / (density_normv * norm_posv[iii]);
	  //mutex.lock();
	  //double tmpvalue = pOldv[iii] * log(starlikev[iii]);
	  //if(tmpvalue < -1e20){tmpvalue=0.0;}
	  //if(tmpvalue != tmpvalue){tmpvalue = 0.0;}
	  //if(isinf(tmpvalue)==1){tmpvalue = 0.0;}
	  starlikev[iii] = pOldv[iii] * log(starlikev[iii]);
	  if(starlikev[iii] < -1e20){starlikev[iii]=0.0;}
	  if(starlikev[iii] != starlikev[iii]){starlikev[iii] = 0.0;}
	  if(isinf(starlikev[iii])==1){starlikev[iii] = 0.0;}
	  //cout << "Final " << tmpvalue << endl;
	  //cout << "Current total " << total_lnLv+tmpvalue << endl;
	  //cout << " " << endl;
	  //total_lnLv += tmpvalue;
	  //mutex.unlock();
	}
    };
      for (auto &&elem : starlikev)
	total_lnLv += elem;
  if (situation > 2)
    {
      ez_thread(threadnum,NThreads)
	{
	  for(int i=threadnum; i<num_maser; i+=NThreads)
	    {
	      int iii = i;
	      starlikevm[iii] = gauss_legendre(100,star_likeZmaser,&iii,0.0,Zmax_starm[iii]);
	      //cout << "INDEX " << iii << endl;
	      //cout << "starlikev  " << starlikevm[iii] << endl;
	      //cout << "Density norm " << density_normv << endl;
	      //cout << "prob old " << pOldvm[iii] << endl;
	      starlikevm[iii] *= (1.0 - cmodv) / density_normv;
	      //mutex.lock();
	      //double tmpvalue = pOldvm[iii] * log(starlikevm[iii]);
	      //if(tmpvalue < -1e20){tmpvalue=0.0;}
	      //if(tmpvalue != tmpvalue){tmpvalue = 0.0;}
	      //if(isinf(tmpvalue)==1){tmpvalue = 0.0;}
	      starlikevm[iii] = pOldvm[iii] * log(starlikevm[iii]);
	      if(starlikevm[iii] < -1e20){starlikevm[iii]=0.0;}
	      if(starlikevm[iii] != starlikevm[iii]){starlikevm[iii] = 0.0;}
	      if(isinf(starlikevm[iii])==1){starlikevm[iii] = 0.0;}
	      //cout << "Final " << tmpvalue << endl;
	      //cout << "Current total " << total_lnLv+tmpvalue << endl;
	      //cout << " " << endl;
	      //total_lnLv += tmpvalue;
	      //mutex.unlock(); 
	    }
	};
      for (auto &&elem : starlikevm)
	total_lnLv += elem;
    }

  //vector<double> testVal = LogLike2(Cube);
  //cout << "Total ln L " << testVal[4] << endl;
  //cout << "Total ln L VVVV " << total_lnLv << endl;
  //cout << "Diff density " << (testVal[0] - density_normv) << endl;
  //cout << "GCOWS diff " << (testVal[1] - density_gcowsv) << endl;
  //cout << "Schodel diff " << (testVal[2] - density_schodel) << endl; 
  //cout << "Diff index " << (testVal[3] - blahDex) << endl;
  //cout << "Diff ln L " << (testVal[4] - total_lnLv) << endl;
  //if (abs(testVal - total_lnLv) > 0.0)
  // {
  //   cout << "NOT THE SAME: " << endl;
  //   cout << "Difference, abs " << (testVal - total_lnLv) << endl;
  //   cout << gmodv << " " << almodv << " " << demodv << " " << brmodv << " " << cmodv << endl;
  // }

  //cout << gmodv << " " << almodv << " " << demodv << " " << brmodv << " " << cmodv << endl;
  lnew = total_lnLv * 1.0;

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
  tmp = argv[11];
  char ending[100];
  strcpy(ending,tmp.c_str());
  max_r = atof(argv[5]);
  max_r *= cm_in_pc;
  Rcut = atof(argv[6]);
  Rcut *= cm_in_pc;
  situation = stoi(argv[7]);
  nonRadial = stoi(argv[10]);
  if (fabs(remainder(situation,2)) < 1.0)
    {
      cout << "Fixing alpha, delta, and r_break"<<endl;
      almodv = atof(argv[2]);
      demodv = atof(argv[3]);
      brmodv = atof(argv[4]);
      brmodv *= cm_in_pc;
    }
  if (situation != 3){cmodv = 0.5;}
  innerCut = atof(argv[8]) * cm_in_pc;
  outerCut = atof(argv[9]) * cm_in_pc;


  //Read in values for stars
  std::string gcowsfile;
  gcowsfile.append("stars_mn");
  gcowsfile.append(ending);
  gcowsfile.append(".dat");
  ifstream in(gcowsfile,ios::in);
  double tmp1,tmp2,tmp3,tmp4,tmpmin,tmpmax,tmpval;
  if (in.is_open())
    {
      while (!in.eof())
	{
	  in >> tmp1 >> tmp2 >> tmp3 >> tmp4;
	  r2dv.push_back(tmp1);
	  arv.push_back(tmp2);
	  arve.push_back(tmp3);
	  pOldv.push_back(tmp4);
	  tmpmin = -1.0*GM / (tmp1*tmp1);
	  min_accelv.push_back(tmpmin);
	  tmpmax = -1.0*GM*tmp1 / (max_r*max_r*max_r);
	  max_accelv.push_back(tmpmax);
	  if (tmp2 < 0.0)
	    {
	      tmpval = tmp3*sqrt(PI/2.0)*(erf((tmp2-tmpmin)/(sqrt(2.0)*tmp3))-erf((tmp2-tmpmax)/(sqrt(2.0)*tmp3)));
	      if (tmpval==0.0)
		{
		  int iii = r2dv.size() - 1;
		  tmpval=gauss_legendre(40,pos_int,&iii,tmpmin,tmpmax);
		}
	      if (tmpval==0.0){cout << "Error in norm pos calculation" << endl;}
	    }
	  else{tmpval = fabs(tmpmax - tmpmin);}
	  norm_posv.push_back(tmpval);
	  tmpval = sqrt(max_r*max_r - tmp1*tmp1);
	  Zmax_star.push_back(tmpval);
	  starlikev.push_back(0.0);
	  amodv.push_back(0.0);
	  like_returnv.push_back(0.0);
	}
    }
  in.close();
  num_stars = r2dv.size();

  if (situation > 2)
    {
      //Read in values for stars in maser field
      cout <<"Reading in stars in maser fields"<< endl;
      std::string maserfile;
      maserfile.append("maser_mn");
      maserfile.append(ending);
      maserfile.append(".dat");
      ifstream in(maserfile,ios::in);
      double tmp1,tmp2,tmp3,tmp4,tmpvalue;
      if (in.is_open())
	{
	  while (!in.eof())
	    {
	      in >> tmp1 >> tmp2;
	      r2dvm.push_back(tmp1);
	      pOldvm.push_back(tmp2);
	      tmpvalue = sqrt(max_r*max_r - tmp1*tmp1);
	      Zmax_starm.push_back(tmpvalue);
	      starlikevm.push_back(0.0);
	      like_returnvm.push_back(0.0);
	    }
	}
      in.close();
      num_maser = r2dvm.size();
    }



  if (nonRadial > 0)
    {
      //Read in gcows field
      cout << "Reading in GCOWS field info"<< endl;
      ifstream in("gcows_field.dat",ios::in);
      double tmp1,tmp2, tmpx, tmpy, tmpR, tmpz;
      if (in.is_open())
	{
	  while(!in.eof())
	    {
	      in >> tmp1 >> tmp2;
	      gcows_vrows.push_back((int) tmp1);
	      gcows_vcols.push_back((int) tmp2);
	      tmpx = sqrt((tmp2-1500.0)*(tmp2-1500.0))*dxy;
	      tmpy = sqrt((tmp1-1500.0)*(tmp1-1500.0))*dxy;
	      tmpR = sqrt(tmpx*tmpx + tmpy*tmpy);
	      Xgcows.push_back(tmpx);
	      Ygcows.push_back(tmpy);
	      Rgcows.push_back(tmpR);
	      tmpz = sqrt(max_r*max_r - tmpR*tmpR);
	      Zmax_gcows.push_back(tmpz);
	      rho_gcows.push_back(0.0);
	    }
	}
      in.close();
      num_gcows = gcows_vrows.size();
    }
  
  std::string pfile;
  pfile.append(root);
  pfile.append("priors.txt");

  //Writing priors to file
  ofstream priors;
  priors.open(pfile);
  priors << "#Gamma priors:\n" << min_g << " " << max_g << "\n";
  if (fabs(remainder(situation,2)) > 0.0)
    {
      priors << "#Alpha priors:\n" << min_a << " " << max_a << "\n";
      priors << "#Delta priors:\n" << min_d << " " << max_d << "\n";
      priors << "#Break r priors (pc):\n" << min_b/cm_in_pc << " " << max_b/cm_in_pc << "\n";
      priors << "#Break r meu and sigma (pc):\n" << meu_b/cm_in_pc << " " << sigma_b/cm_in_pc << "\n";
    }
  priors.close();

	
  // set the MultiNest sampling parameters
	
  int mmodal = 0;// do mode separation?
	
  int ceff = 0;// run in constant efficiency mode?
	
  int nlive = 1000;// number of live points
	
  double efr = 0.8;// set the required efficiency
	
  double tol = 0.5;// tol, defines the stopping criteria
	
  int ndims = 5;// dimensionality (no. of free parameters)
	
  int nPar = 5;// total no. of parameters including free & derived parameters
	
  int nClsPar = 5;// no. of parameters to do mode separation on

  if (fabs(remainder(situation,2)) < 1.0)
    {
      cout <<"Fixing alpha, delta, and r_break"<< "\n";
      ndims = 1;
      nPar = 1;
      nClsPar = 1;
    }
  if (situation == 1.0)
    {
      cout <<"Fixing C factor, not using Schodel/maser"<< "\n";
      ndims = 4;
      nPar = 4;
      nClsPar = 4;
    }
	
  int updInt = 100;// after how many iterations feedback is required & the output files should be updated
                    // note: posterior files are updated & dumper routine is called after every updInt*10 iterations
	
  double Ztol = -1E90;// all the modes with logZ < Ztol are ignored
	
  int maxModes = 100;// expected max no. of modes (used only for memory allocation)
	
  int pWrap[ndims];// which parameters to have periodic boundary conditions?
  for(int i = 0; i < ndims; i++) pWrap[i] = 0;
	
  int seed = -1;// random no. generator seed, if < 0 then take the seed from system clock
	
  int fb = 1;// need feedback on standard output?
	
  int resume = 0;// resume from a previous job?
	
  int outfile = 1;// write output files?
	
  int initMPI = 1;// initialize MPI routines?, relevant only if compiling with MPI
		  // set it to F if you want your main program to handle MPI initialization
	
  double logZero = -1E100;// points with loglike < logZero will be ignored by MultiNest
	
  int maxiter = 0;// max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it 
		  // has done max no. of iterations or convergence criterion (defined through tol) has been satisfied
	
  void *context = 0;// not required by MultiNest, any additional information user wants to pass	
	
  // calling MultiNest
	
  nested::run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI,logZero, maxiter, LogLike, dumper, context);
}
