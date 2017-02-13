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

//char root[100] = "/u/schappell/pmnOld/c_accel_test_pert_more_maxr0.4_-0.2_4_4_0.5_";// root for output files
//define limits of priors
double min_g = -3.0;
double max_g = 2.0;
double min_d = 2.0;
double max_d = 10.0;
double min_a = 0.0;
double max_a = 10.0;
double min_b = 5.0*dist*cm_in_au;
double max_b = 2.0*cm_in_pc;
//double max_rh = 5.0*dist*cm_in_au;
double plate_scale = 0.00995;
double dxy = plate_scale * dist * cm_in_au; //dx or dy of a pixel
//double max_Z = 10.0 * cm_in_pc;
//double max_R = 1.8 * dist * cm_in_au;
//double maxr = sqrt(max_Z*max_Z + max_R*max_R);

int situation,nonRadial;

vector<double> r2d, r2dm;
vector<double> ar;
vector<double> are;
vector<double> pOld, pOldm;
vector<double> min_accel;
vector<double> max_accel;
vector<double> maxz_star;
vector<int> gcows_rows, gcows_cols;
//int gcowsRows = 3000;
//int gcowsCols = 3000;
//int gcowsField[3000][3000];
//vector< vector<int> > gcowsField;

double gmod, almod, demod, brmod, cmod, rhmod,rho0mod;
double max_r, max_R, max_rho0,max_Z,min_Z;
double use_maxZ, density_norm, density_norm_maser;
//for certain models the max Z needs to be smaller for integration to work
double density_normM, density_box, density_cyl, use_min_accel;
int stardex, maserdex;
int num_stars, num_maser, num_gcows;
double max_X,max_Y,min_X,min_Y,innerCut,outerCut,ZprimeM,Rcut,Xprime,Yprime,Rprime,Rcyl,norm_pos;
  

double density_int(double Rprime, double Zprime, void* data)
{
  //cout << "Density integration over R and z" << "\n";
  double tmpvalue = pow((sqrt(Rprime*Rprime + Zprime*Zprime)/brmod),demod);
  return fabs(2.0*PI*Rprime*pow((Rprime*Rprime + Zprime*Zprime)/(brmod*brmod),(gmod/-2.0))*
	      pow((1.0+tmpvalue),((gmod-almod)/demod)));
}

double density_intZ(double Zprime, void* data)
{
  double tmpvalue = pow((sqrt(Xprime*Xprime + Yprime*Yprime + Zprime*Zprime)/brmod),demod);
  return fabs(dxy*dxy*pow((Xprime*Xprime + Yprime*Yprime + Zprime*Zprime)/(brmod*brmod),(gmod/-2.0))*
	      pow((1.0+tmpvalue),((gmod-almod)/demod)));
}

double density_intZcyl(double Zprime, void* data)
{
  double tmpvalue = pow((sqrt(Rcyl*Rcyl + Zprime*Zprime)/brmod),demod);
  return fabs(2.0*PI*Rcyl*pow((Rcyl*Rcyl + Zprime*Zprime)/(brmod*brmod),(gmod/-2.0))*
	      pow((1.0+tmpvalue),((gmod-almod)/demod)));
}


double density_intR(double Rprime, void* data)
{
  max_Z = sqrt(max_r*max_r - Rprime * Rprime);
  if (Rprime >= rhmod)
    {
      min_Z = 0.0;
    }
  else
    {
      min_Z = sqrt(rhmod*rhmod - Rprime*Rprime);
    }
  Rcyl = Rprime * 1.0;
  return gauss_legendre(100,density_intZcyl,NULL,min_Z,max_Z);
}

//double density_intZ_ad(double Zprime)
//{
//  double tmpvalue = pow((sqrt(Xprime*Xprime + Yprime*Yprime + Zprime*Zprime)/brmod),demod);
//  return fabs(dxy*dxy*pow((Xprime*Xprime + Yprime*Yprime + Zprime*Zprime)/(brmod*brmod),(gmod/-2.0))*
//	      pow((1.0+tmpvalue),((gmod-almod)/demod)));
//}

double density_sphere(double r3dprime, void* data)
{
  //cout << "Density integration over R and z" << "\n";
  double tmpvalue = pow((r3dprime/brmod),demod);
  return fabs(4.0*PI*r3dprime*r3dprime*pow((r3dprime/brmod),gmod)*
	      pow((1.0+tmpvalue),((gmod-almod)/demod)));
}

//double density_xyz(double Xprime, double Yprime, void* data)
//{
//  double rprime = sqrt(pow(Xprime,2.0) + pow(Yprime,2.0) + pow(ZprimeM,2.0));
//  double tmpvalue = pow(rprime/brmod,demod);
//  return fabs(pow(rprime,(-1.0*gmod))*pow(1.0+tmpvalue,(gmod-almod)/demod));
//}

//double density_int_box(double ZprimeM, void* data)
//{
//  return gauss_legendre_2D_cube(40,density_xyz,NULL,min_X,max_X,min_Y,max_Y);
//}

double pos_int(double ax, void* data){return exp(-1.0*(ar[stardex]-ax)*(ar[stardex]-ax)/(2.0*are[stardex]*are[stardex]));}

double star_likeZ(double z0mod, void* data)
{
  //cout << z0mod << "\n";
  //double like_pos;
  //if ((ar[stardex]) < 0.0)
  //{
  //  double amod = -1.0*GM*r2d[stardex] / pow((sqrt(r2d[stardex]*r2d[stardex] + z0mod*z0mod)),3.0);
  //  like_pos = exp(-1.0*(ar[stardex]-amod)*(ar[stardex]-amod)/(2.0*are[stardex]*are[stardex]));
  //}
  //else
  //{
  //  like_pos = 1.0;
      //like_pos = pow((r2d[stardex]*r2d[stardex] + z0mod*z0mod),2.5) / (3.0*z0mod*GM*r2d[stardex]);
      //norm_pos = maxz_star[stardex];
  //}

  double tmp_value = pow((sqrt(r2d[stardex]*r2d[stardex]+z0mod*z0mod)/brmod),demod);
  double like_den = pow((1.0+tmp_value),((gmod-almod)/demod));
  if(isinf(like_den)==1)
    {
      cout << "It is power of gamma - alpha / delta that's INF";
      cout << "\n";
      cout << "Gamma is " << gmod << ", alpha is " << almod << ", and delta is " << demod;
      cout << "\n";
    }
  like_den *= rho0mod * pow((r2d[stardex]*r2d[stardex]+z0mod*z0mod)/(brmod*brmod),(gmod/-2.0));
  //like_den *= pow((r2d[stardex]*r2d[stardex]+z0mod*z0mod)/(brmod*brmod),(gmod/-2.0));
  //double like_return = like_pos*like_den;// / (norm_pos * density_norm);
  double like_return = like_den * 1.0;
  //like_return *= cmod;
  if(like_return==0.0){like_return=1e-323;}
  if(like_return != like_return){like_return=1e-323;}
  if(isinf(like_return)==1)
    {
      cout << "Got oo in summation";
      cout << "\n";
    }
  //cout << "index " << stardex << "\n";
  //cout << "Gamma is " << gmod << ", alpha is " << almod << ", and delta is " << demod;
  //cout << "\n";
  //cout << "z mod = " << z0mod << "\n";
  //cout << "like pos = " << like_pos << "\n";
  //cout << "like den = " << like_den << "\n";
  //cout << "den norm = " << density_norm << "\n";
  //cout << "pos norm = " << norm_pos << "\n";
  //cout << "tmp value = " << tmp_value << "\n";
  return like_return;
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
  like_den *= pow((r2dm[maserdex]*r2dm[maserdex]+z0mod*z0mod)/(brmod*brmod),(gmod/-2.0));
  //like_den *= rho0mod * cmod;
  double like_return = like_den * 1.0; // / density_norm;//_maser;
  //like_return *= (1.0 - cmod);
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
  if (fabs(remainder(situation,2)) > 0.0)
    {
      almod = Cube[1] * (max_a - min_a) + min_a;
      demod = Cube[2] * (max_d - min_d) + min_d;
      brmod = Cube[3] * (max_b - min_b) + min_b;
      if (situation == 3){cmod = Cube[4];}
      //rho0mod = Cube[5] * max_rho0;
      //rhmod = Cube[5] * max_rh;
    }

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
	      if (Rprime >= rhmod)
		{
		  min_Z = 0.0;
		}
	      else
		{
		  min_Z = sqrt(rhmod*rhmod - Rprime*Rprime);
		}
	      //Rcyl = Rprime * 1.0;
	      Xprime = sqrt((gcows_cols[i]-1500.0)*(gcows_cols[i]-1500.0))*dxy;
	      Yprime = sqrt((gcows_rows[i]-1500.0)*(gcows_rows[i]-1500.0))*dxy;
	      
	      double tmpDen = gauss_legendre(100,density_intZ,NULL,min_Z,max_Z);

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
		  use_maxZ = 5.0 * brmod * pow(fabs(gmod/almod),(1.0/demod));
		  if (use_maxZ < max_Z)
		    {
		      tmpDen = gauss_legendre(100,density_intZ,NULL,min_Z,use_maxZ);
		      cout << "Using 5 times Z_rho,max as int limits" << "\n";
		      if (tmpDen <= 0.0)
			{
			  use_maxZ = 3.0 * brmod * pow(fabs(gmod/almod),(1.0/demod));
			  tmpDen = gauss_legendre(100,density_intZ,NULL,min_Z,use_maxZ);
			  cout << "Using 3 times Z_rho,max as int limits" << "\n";
			  if (tmpDen <= 0.0)
			    {
			      use_maxZ = 1.5 * brmod * pow(fabs(gmod/almod),(1.0/demod));
			      tmpDen = gauss_legendre(100,density_intZ,NULL,min_Z,use_maxZ);
			      cout << "Using 1.5 times Z_rho,max as int limits" << "\n";
			      if (tmpDen <= 0.0)
				{
				  Rcyl = Rprime * 1.0;
				  cout << "Still have a problem with density integration, i " << i << "\n";
				  cout <<"\n";
				  cout <<"Gamma is "<<gmod<<", alpha is "<<almod<<", delta is "<<demod<<", and break r is "<<brmod;
				  cout <<"\n";
				  tmpDen = 0.0;
				}
			    }
			}
		    }
		  else{cout <<"GCOWS part of density norm is 0, make max r larger"<<"\n";}
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
      double density_hole = gauss_legendre(100,density_sphere,NULL,0,rhmod);
      density_norm -= density_hole;

      if (density_norm <= 0.0)
	{
	  cout <<"GCOWS part of density norm is 0, make max r larger"<<"\n";
	}
    }
  density_norm *= cmod;
  //density_norm_ad *= cmod;
  //cout << " " << "\n";
  //cout << "GCOWS norm: " << density_norm/cmod << "\n";
  //cout << "New norm: " << density_norm_ad << "\n";
  //cout << "Acceleration density norm = " << density_norm << "\n";
  if(situation > 2)
    {
      //density_norm_maser = 0.0;
      //cout << "Density integration for Schodel/maser" << "\n";
      double density_tmp = gauss_legendre(100,density_intR,NULL,innerCut,outerCut);
      if (density_tmp <= 0.0)
	{
	  cout << "Schodel part of density norm is 0, make max r larger" << "\n";
	}
      //density_tmp *= (1.0 - cmod);
      density_norm += (1.0 - cmod) * density_tmp;
      //density_norm_maser += density_tmp;
      //cout << "Schodel density norm = " << density_tmp << "\n";
      //double for_print = gauss_legendre(100,density_intR,NULL,0.0,innerCut);
      //cout << "Cyl integral between 0 and 5 as " << for_print << "\n";
    }
  //cout <<"Gamma is "<<gmod<<", alpha is "<<almod<<", delta is "<<demod<<", and break r is "<<brmod<<"\n";
  //cout <<"Density norm = " << density_norm << "\n";
  //cout <<"Max R = " << max_R << " and Max Z = " << max_Z << "\n";
  //cout << "Density norm end" << "\n";

  //density_norm *= rho0mod;
  //density_norm_maser *= rho0mod;

  double norm_pos,min_Z;
  for(int i=0; i<num_stars; i++)
    {
      stardex = i;
      //cout << ar[i] << " " << are[i] << " " << pOld[i] << "\n";
      //cout << "Begin starlike" << "\n";
      if (r2d[stardex] >= rhmod)
	{
	  use_min_accel = min_accel[stardex] * 1.0;
	  min_Z = 0.0;
	}
      else
	{
	  use_min_accel = -1.0 * GM / (rhmod*rhmod);
	  min_Z = sqrt(rhmod*rhmod - r2d[stardex]*r2d[stardex]);
	}
      //if ((ar[stardex]) < 0.0)
      //{
      //  norm_pos = are[stardex]*sqrt(PI/2.0)*(erf((ar[stardex]-use_min_accel)/(sqrt(2.0)*are[stardex]))
      //					-erf((ar[stardex]-max_accel[stardex])/(sqrt(2.0)*are[stardex])));
      //  if(norm_pos==0.0){norm_pos=gauss_legendre(40,pos_int,NULL,use_min_accel,max_accel[stardex]);}
      //  norm_pos=fabs(norm_pos);
      //  if(norm_pos==0.0)
      //    {
      //      cout << "Norm pos is zero";
      //      cout <<"\n";
      //    }
      //}
      //else
      //{
      //  norm_pos = fabs(max_accel[stardex] - use_min_accel);
      //}
      double max_Z = sqrt(max_r*max_r - r2d[stardex]*r2d[stardex]);
      starlike = gauss_legendre(100,star_likeZ,NULL,min_Z,max_Z);
      //cout << starlike << "\n";
      if (starlike == 0.0){starlike=1e-323;}
      //cout << "pos and density like " << starlike << "\n";
      //starlike *= cmod / (norm_pos * density_norm);
      starlike *= cmod / (density_norm);
      total_lnL += pOld[i]*log(starlike);
      //cout << starlike << "\n";
      //cout << "norm pos " << norm_pos << "\n";
      //cout << "norm density " << density_norm << "\n";
      //cout << "c mod " << cmod << "\n";
      //cout << total_lnL << "\n";
      //cout << "End starlike" << "\n";
    }

  if (situation > 2)
    {
      for(int m=0; m<num_maser; m++)
	{
	  //cout << "Cycle through schodel/maser stars" << "\n";
	  maserdex = m;
	  //cout << ar[i] << " " << are[i] << " " << pOld[i] << "\n";
	  //cout << "Begin starlike" << "\n";
	  double max_Z = sqrt(max_r*max_r - r2dm[maserdex]*r2dm[maserdex]);
	  starlike = gauss_legendre(100,star_likeZmaser,NULL,0,max_Z);
	  //cout << starlike << "\n";
	  if (starlike == 0.0){starlike=1e-323;}
	  //cout << "density like " << starlike << "\n";
	  starlike *= (1.0 - cmod) / density_norm;
	  total_lnL += pOldm[m]*log(starlike);
	  //cout << starlike << "\n";
	  //cout << "\n";
	  //cout << "R2d " << r2dm[maserdex] << "\n";
	  //cout << "norm density " << density_norm << "\n";
	  //cout << "c mod " << (1.0 - cmod) << "\n";
	  //cout << starlike << "\n";
	  //cout << total_lnL << "\n";
	  //cout << "End starlike" << "\n";
	}
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
  tmp = argv[11];
  char ending[100];
  strcpy(ending,tmp.c_str());
  max_r = atof(argv[5]);
  max_r *= cm_in_pc;
  Rcut = atof(argv[6]);
  Rcut *= cm_in_pc;
  //maxr = sqrt(max_Z*max_Z + max_R*max_R);
  situation = stoi(argv[7]);
  nonRadial = stoi(argv[10]);
  if (fabs(remainder(situation,2)) < 1.0)
    {
      cout << "Fixing alpha, delta, and r_break"<<"\n";
      almod = atof(argv[2]);
      demod = atof(argv[3]);
      brmod = atof(argv[4]);
      brmod *= cm_in_pc;
    }
  if (situation != 3){cmod = 0.5;}
  //cmod = 0.31;
  rhmod = 0.0;
  rho0mod = 1.0;
  //max_R = Rcut * 1.0;
  innerCut = atof(argv[8]) * cm_in_pc;
  outerCut = atof(argv[9]) * cm_in_pc;
      //min_X = -1.0 * outerCut;
      //min_Y = -1.0 * outerCut;
      //max_X = 1.0 * outerCut;
      //max_Y = 1.0 * outerCut;

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
	  //tmpvalue = sqrt(tmp1*tmp1 + max_Z*max_Z);
	  tmpvalue = -1.0*GM*tmp1 / (max_r*max_r*max_r);
	  max_accel.push_back(tmpvalue);
	  // tmpvalue = sqrt(maxr*maxr - tmp1*tmp1);
	  //maxz_star.push_back(tmpvalue);
	}
    }
  in.close();
  num_stars = r2d.size();

  if (situation > 2)
    {
      //Read in values for stars in maser field
      cout <<"Reading in stars in maser fields"<<"\n";
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
	      r2dm.push_back(tmp1);
	      pOldm.push_back(tmp2);
	    }
	}
      in.close();
      num_maser = r2dm.size();
    }

  if (nonRadial > 0)
    {
      //Read in gcows field
      cout << "Reading in GCOWS field info"<<"\n";
      ifstream in("gcows_field.dat",ios::in);
      //cout << "opened file" << "\n";
      //int gcowsField[gcowsRows][gcowsCols];
      //cout << "defined array" << "\n";
      double tmp1,tmp2;
      if (in.is_open())
	{
	  //cout << "In if statement" << "\n";
	  while(!in.eof())
	    {
	      in >> tmp1 >> tmp2;
	      //cout << tmp1 << " " << tmp2 << "\n";
	      gcows_rows.push_back((int) tmp1);
	      gcows_cols.push_back((int) tmp2);
	    }
	}
	  //cout << "In while statement" << "\n";
	  //for(int row=0; row<gcowsRows; row++)
	      //{
	      //cout << "In for statement" << "\n";
	      //cout << row << "\n";
	      //for(int col=0; col<gcowsCols; col++)
	      //{
		  //cout << "In for for statement" << "\n";
		  //gcowsFile >> tmp;
		  //if (tmp > 0){
	      //cout << "in GCOWS field" << "\n";
	      // }
	      //  gcowsField[row][col] = tmp * 1;
	      //}
	      //}
	  //}
	      //}
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
  //cout << "situation equals " << situation << "\n";
  //cout << "Remainder of situation " << fabs(remainder(situation,2)) << "\n";
  if (fabs(remainder(situation,2)) > 0.0)
    {
      priors << "#Alpha priors:\n" << min_a << " " << max_a << "\n";
      priors << "#Delta priors:\n" << min_d << " " << max_d << "\n";
      priors << "#Break r priors (pc):\n" << min_b/cm_in_pc << " " << max_b/cm_in_pc << "\n";
      //priors << "#Hole r priors (pc):\n" << 0.0 << " " << max_rh/cm_in_pc << "\n";
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
