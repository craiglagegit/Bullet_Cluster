/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, NYU CCPP
  Date: Feb 14, 2014

  Triaxial halo generation program

*/

//****************** halo.h **************************

#include <stdio.h>       
#include <stdlib.h>      
#include <math.h>        
#include <string.h>
#include <vector>
#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/iter_find.hpp>
#include <sstream>
#include <iostream>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include <globals.h>
#include <fileio.h>
#include <potential.h>

using namespace std;

struct Particle
{
  int ID;             // ID
  double mass;        // Mass
  double x[3];        // Position
  double v[3];        // Velocity
  int grid;           // Identifies in which grid it falls
};

class Halo
{
 public:
  double	G;		// Gravitational Constant
  double	H;		// Hubble Constant
  double	P;		// Ratio of Z to X axes - P<1, Q=1 => oblate
  double	Q;		// Ratio of Y to X axes - P=1, Q<1 +> prolate

  int	        NPot;		// Number of steps in Potential array
  int           NLook;          // Number of steps in look-up table
  double	M200;		// Total mass of halo
  double	C;		// Halo concentration parameter
  double	R200;		// Radius at M200
  double	Rs;		// Scale Radius
  double        MaxR;           // Maximum DM particle radius
  double 	rmin;		// min R used in U(R) lookup
  double 	mmin;		// min m used in R(m) llokup
  double	GF;		// Gas Fraction

  double	RhoNorm;
  double	RhoDMNorm;
  double	GasAlpha;	// Parameters used in the gas density model
  double        GasBeta;
  double	GasRcool;
  double 	GasEpsilon;
  double 	GasRs;
  double	GasBeta2;
  double	GasRcool2;
  double	GasN2;

  int	        NTot;		// Total number of particles
  int	        NDM;		// Number of DM particles
  int	        NGas;		// Number of gas particles
  int	        NDM200;		// Number of DM particles within R200
  int	        NGas200;	// Number of gas particles within R200

  double* Xoffset;            // Shift of halo origin relative to collision origin
  double* Voffset;            // Velocity shift of halo origin relative to collision origin
  double Phi;                 // Halo rotation Euler angle
  double Theta;               // Halo rotation Euler angle
  double Psi;                 // Halo rotation Euler angle
  double** EulerPlus;         // Rotation matrix from collision system to halo system
  double** EulerMinus;        // Rotation matrix from halo system to collision system

  double Ufactor;               // Internal energy scale factor
  double Rhofactor;             // Density scale factor

  double Fntp0;
  double Fntpexp;
  string outfile;		// output file name
  string profile;		// profiles output file name
  string potentialfile;	        // file with potential data
  string nbodyfile;	        // file with nbody data

  gsl_spline* uinterp;   // For U(R)
  gsl_interp_accel* uacc;

  Particle* Part;// This contains the DM particle info

  Potential* Pot; // This contains the potential information

  Halo() {};
  Halo(string);
  ~Halo();

  double Rho(double);
  double RhoDM(double);
  double MassDM(double);
  double Mass(double);
  void SetRhoNorm();
  double UCalc(double);
  double U(double);
  double NormR(double, double);
  double Rxyz(double*);
  double Fntp(double);
  void ReadNbodyData();
  void EulerRotateShift(int);
  void EulerInitialize();
};

double MassShellDM(double, void*); // These are at top level to enable simple passing to gsl_integrate
double MassShell(double, void*); // These are at top level to enable simple passing to gsl_integrate
double UIntegrand(double, void*); // These are at top level to enable simple passing to gsl_integrate
