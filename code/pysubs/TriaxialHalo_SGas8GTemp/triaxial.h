/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, NYU CCPP
  Date: May 25, 2012

  Triaxial halo generation program

*/

#include <stdio.h>              // for printf, fprintf, fscanf
#include <stdlib.h>             // for calloc, free
#include <math.h>               // for sqrt, exp
#include <string.h>             // strcpy, strtok...
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#define pi 3.14159265358979     // Pi
#define STRING_LENGTH 128

#define GRAVITY    6.672e-8
#define MSOL       1.989e33
#define CM_PER_MPC 3.085678e24
#define HUBBLE     0.72		// Hubble Constant in 100km/sec / Mpc

#define UMASS      ( MSOL       * 1e10 )
#define ULENGTH    (CM_PER_MPC / 1000.0 )
#define UVELOCITY  (             1e5   )
#define UTIME      ( ULENGTH / UVELOCITY )
#define UDENSITY   ( UMASS / ( ULENGTH * ULENGTH * ULENGTH ) )
#define UPRESSURE  ( UMASS / ( ULENGTH * UTIME * UTIME ) )
#define UENERGY    ( UMASS * UVELOCITY * UVELOCITY )

#define CM         ( 1.0 / ULENGTH )
#define GRAM       ( 1.0 / UMASS )
#define SEC        ( 1.0 / UTIME )
#define ERGS       ( GRAM * ( CM * CM ) / ( SEC * SEC ) )

#define keV        ( 1.602e-9 * ERGS )

#define GAMMA_MINUS1 0.666
#define PROTONMASS   1.6726e-24
#define BOLTZMANN    1.3806e-16
#define HYDROGEN_MASSFRAC 0.76


//************DATA STRUCTURES*************

struct Halo
// This data structure holds the various values which determine the halo structure and the various constants. 
{
double	G;		// Gravitational Constant
double	H;		// Hubble Constant
double	P;		// Ratio of Z to X axes - P<1, Q=1 => oblate
double	Q;		// Ratio of Y to X axes - P=1, Q<1 +> prolate

int	NPot;		// Number of steps in Potential array
int	NLook;		// Number of steps in Look-up table for R(m) and U(R)

double	M200;		// Total mass of halo
double	C;		// Halo concentration parameter
double	R200;		// Radius at M200
double	Rs;		// Scale Radius
double  MaxR;           // Maximum DM particle radius
double 	rmin;		// min R used in U(R) lookup
double 	mmin;		// min m used in R(m) llokup
double	GF;		// Gas Fraction

double	RhoNorm;
double	RhoDMNorm;
double	GasAlpha;	// Parameters used in the gas density model
double  GasBeta;
double	GasRcool;
double 	GasEpsilon;
double 	GasRs;
double	GasBeta2;
double	GasRcool2;
double	GasN2;

int	NTot;		// Total number of particles
int	NDM;		// Number of DM particles
int	NGas;		// Number of gas particles
int	NDM200;		// Number of DM particles within R200
int	NGas200;	// Number of gas particles within R200

char	outfile[STRING_LENGTH];		// output file name
char	profile[STRING_LENGTH];		// profiles output file name
char	potentialfile[STRING_LENGTH];	// file with potential data
char	nbodyfile[STRING_LENGTH];	// file with nbody data
};

struct Potential
// This data structure holds the potential data. 
{
double*	xa[3];
double* phia[3];
double* fa[3];	
gsl_spline* phi[3];
gsl_spline* f[3];
gsl_interp_accel* phiacc[3];
gsl_interp_accel* facc[3];
};

struct UMass
// This data structure holds look-up arrays for the gas mass, and U(R). 
{
double*	rx;
double* ux;
double*	r;
double* u;
gsl_spline* rinterp;
gsl_spline* uinterp;
gsl_interp_accel* racc;
gsl_interp_accel* uacc;
};	

struct PQPhi
// This data structure holds look-up arrays for Pphi, and Qphi. 
{
double*	rx;
double* pphi;
double* qphi;
gsl_spline* pphiinterp;
gsl_spline* qphiinterp;
gsl_interp_accel* pphiacc;
gsl_interp_accel* qphiacc;
};	
	

struct UParams
{
struct Halo *H;
struct Potential *Pot;
struct PQPhi *PQ;
};

// Definition of Data Structures in GADGET files

struct GadgetHeader
{ 
  int   iNumberOfParticles[6];
  double  rParticleMass[6];            // mass of each particle type in units of 10^10 Msol/h
  double  rExpansionFactor;            // ( Time / a )
  double  rRedShift;                   // Set only for cosmological simulations
  int   iStarFormationFlag;
  int   iFeedbackFlag;
  int   iTotalNumberOfParticles[6];  // Only differs from NumberOfParticles if file is split into chunks.
  int   iCoolingFlag;
  int   iNumberOfFiles;
  double  rBoxSize;                    // Relevant only for periodic boundaries
  double  rOmega0;                     // matter density at z=0 (cos)
  double  rOmegaLambda;                // Vacuum energy density at z=0 (cos)
  double  rHubbleConstant;             // Value of the hubble constant in units of 100 km/Mpc*sec
  int FlagAge;
  int FlagMetals;
  int NallHW[6];
  int flag_entr_ics;
  char   bUnused[60];
};


struct GadgetParticles
{
  int   iParticleID;           // unique particle identification number
  int   iParticleType;	       // 0 - Gas; 1 - Dark Matter
  float  rPosition[3];         // Comoving coordinates in internal length units - kpc
  float  rVelocity[3];         // Particle velocity in internal velocity units = km/sec
  float  rMass;                // Mass in internal units 10^10 Msol
  float  rU;                   // Internal energy per mass unit for gas particles (unit: (km/sec)^2)
  float  rRho;                 // Comoving density of SPH particles in units of:  10^10 (Msol/h) / (kpc/h)^3 
  float  rNe;                  // Electron Density in units of the number of hydrogen nuclei.
  float  rAxisDistance;        // UNOFFICIAL! Used for z-Densities, etc.
  float  rOriginDistance;      // UNOFFICIAL!
  float  rLongitude;
  float  rLatitude;  
  float  e;
  float  i2;
  float  i3;
};

enum {
  OK,
  INVALID_PARAMETERS,
  COULD_NOT_OPEN_FILE,
  READ_ERROR,
  INSUFFICIENT_MEMORY,
  WRITE_ERROR
};

// *******************FUNCTION DECLARATIONS**************

// ****************** fileio.c **************************

struct Halo ReadConfig(const char *inname);
struct Potential ReadPotentialData(const char *inname, struct Halo *H);
void ReadNbodyData(const char *inname, struct GadgetParticles *P, struct Halo *H);
int ReadPar(const char *fname, const char *parnam, char *parval);
double GetDoubleParam(const char *fname, const char *parnam, double defval);
int GetIntParam(const char *fname, const char *parnam, int defval);
void GetStringParam(const char *fname, const char *parnam, char *defval, char* result);
struct GadgetHeader CreateGadgetHeader(struct Halo *H);
int WriteGadgetFile(char *filename, struct GadgetHeader Header, struct GadgetParticles *Particles);
void WriteProfiles(char *filename,struct Halo *H, struct UMass *UM, struct PQPhi *PQ, struct Potential *Pot);

// ****************** utils.c **************************

double Pphi(double R, struct PQPhi *PQ, struct Halo *H);
double Qphi(double R, struct PQPhi *PQ, struct Halo *H);
double AEllipse(double R, struct PQPhi *PQ, struct Halo *H);
double AEllipseDM(double R, struct Halo *H);
double GaussianRandom(double mean, double sigma);
void RandomDirection(double R, double *theta, double *phi, double *x, double *y, double *z, double P, double Q);
//void CharacterizeDistribution(struct GadgetParticles *P, double *G, struct Halo *H);
//void CleanUp(struct GadgetParticles *P, struct Halo *H);
//void Symmetrize(struct GadgetParticles *P, struct Halo *H);
//void Randomize(int axis, struct GadgetParticles *P, struct Halo *H);

// ****************** models.c **************************

double Rho(double R, struct Halo *H);
double RhoDM(double R, struct Halo *H);
double MassShell(double R, void *p);
double MassShellDM(double R, void *p);
double Mass(double R, struct PQPhi *PQ, struct Halo *H);
double MassDM(double R, struct Halo *H);
void SetRhoNorm(struct PQPhi *PQ, struct Halo *H);
double UIntegrand(double R, void *params);
double U(double R, struct Potential *Pot, struct Halo *H);
double NormR(double r, double m, struct PQPhi *PQ, struct Halo *H);
double Radius(double m, struct PQPhi *PQ, struct Halo *H);
struct UMass UMLookup(struct Potential *Pot, struct PQPhi *PQ, struct Halo *H);
struct PQPhi PQLookup(struct Potential *Pot, struct Halo *H);
double RPhi(int axis ,double phi, double rtry, struct Potential *Pot, struct Halo *H);





