/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, NYU CCPP
  Date: Feb 14, 2014

  Triaxial halo generation program

*/

//****************** halo.cpp **************************

#include "halo.h"

Halo::Halo(string inname) //Constructor                                                                                            
{
  // This reads in the data for the halo from the configuration file
  // and initializes the various halo functions
  int i;
  G = GRAVITY * ((UMASS * gsl_pow_2(UTIME)) / (ULENGTH * ULENGTH * ULENGTH));
  H = HUBBLE * ((100 * 1e5 * ULENGTH) / (CM_PER_MPC * UVELOCITY)) * (0.315 * pow(1.296,3.0) + 0.685);

  P = GetDoubleParam(inname, "P", 0.6);
  Q = GetDoubleParam(inname, "Q", 0.8);
  P = GSL_MIN(P, 1.00);
  Q = GSL_MIN(Q, 1.00);
  NPot = GetIntParam(inname, "NPot", 100);
  NLook = GetIntParam(inname, "NLook", 100);
  M200 = GetDoubleParam(inname, "M200", 10000.0);
  C = GetDoubleParam(inname, "C", 3.0);
  R200 = pow((M200 * G) / (100.0 * gsl_pow_2(H)), 1.0 / 3.0);
  Rs = R200 / C;

  mmin = 1.0E-6; //Min value in lookup tables

  RhoNorm = 1.0;
  RhoDMNorm = 1.0;
  GasAlpha = GetDoubleParam(inname, "GasAlpha", 0.0);
  GasBeta = GetDoubleParam(inname, "GasBeta", 0.5);
  GasRcool = GetDoubleParam(inname, "GasRcool", 30.0);
  GasEpsilon = GetDoubleParam(inname, "GasEpsilon", 0.0);
  GasRs = GetDoubleParam(inname, "GasRs", 1.0);
  GasBeta2 = GetDoubleParam(inname, "GasBeta2", 2.0);
  GasRcool2 = GetDoubleParam(inname, "GasRcool2", 10.0);
  GasN2 = GetDoubleParam(inname, "GasN2", 100.0);
  GasRcool = GasRcool / Rs;
  GasRs = GasRs / Rs;
  GasRcool2 = GasRcool2 / Rs;

  NGas = GetIntParam(inname, "NGas", 10000);
  NDM = GetIntParam(inname, "NDM", 10000);
  NTot = NDM + NGas;
  GF = GetDoubleParam(inname, "GF", 0.16);

  Ufactor = G * M200 / Rs;
  Rhofactor = M200 * GF / pow(Rs, 3.0);

  Fntp0 = GetDoubleParam(inname, "Fntp0", 0.5);
  Fntpexp = GetDoubleParam(inname, "Fntpexp", 0.8);

  Xoffset = new double[3];
  Voffset = new double[3];
  outfile = GetStringParam(inname, "outfile", "output.dat");
  profile = GetStringParam(inname, "profile", "profile.dat");
  potentialfile = GetStringParam(inname, "potentialfile", "potential.dat");
  nbodyfile = GetStringParam(inname, "nbodyfile", "nbody.dat");

  Pot = new Potential(potentialfile, NPot); //This reads in and initializes the potential

  // Now we build a look-up table for U(R).  R is on a logarithmic scale
  double rmin = Pot->rmin, rmax = 30.0 * C, dr;
  dr = log10(rmax / rmin) / ((double)NLook - 1);
  double* ux = new double[NLook];
  double* u = new double[NLook];

  for (i=0; i<NLook; i++)
  {
	ux[i] = rmin * pow(10.0,(double)i * dr);
	u[i] = UCalc(ux[i]);
  }

  uinterp = new gsl_spline;
  uacc = new gsl_interp_accel;

  uacc = gsl_interp_accel_alloc();
  uinterp = gsl_spline_alloc (gsl_interp_cspline, NLook);  
  gsl_spline_init (uinterp,  ux, u, NLook);
  delete[] ux;
  delete[] u;
  
  SetRhoNorm();
  ReadNbodyData();
  return;
}

Halo::~Halo() //Destructor                                                                                            
{
  int i;
  delete Pot;
  delete uacc;
  delete uinterp;
  delete[] Part;
  delete[] Xoffset;
  delete[] Voffset;
  for (i=0; i<3; i++)
    {
      delete[] EulerPlus[i];
      delete[] EulerMinus[i];
    }
  delete[] EulerPlus;
  delete[] EulerMinus;
}

double Halo::RhoDM(double R)
{
  // This is the DM density at radius R
  return RhoDMNorm /((R / Rs) * gsl_pow_2(1.0 + (R / Rs)));
}

double MassShellDM(double R, void* p)
{
  // This is the mass in a shell of radius R
  // It is at top level to enable passing to gsl_integration
  Halo* Haloptr = (Halo*)p;
  return Haloptr->RhoDM(R) * 4.0 * M_PI * gsl_pow_2(R) * Haloptr->P * Haloptr->Q;
}

double Halo::MassDM(double R)
{
    // This is the DM mass inside radius R
    gsl_function F;
    F.function=&MassShellDM;
    F.params=this;
    double result, error;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
    gsl_integration_qag(&F, 0.0, R, 0.0, 1.0E-3, 1000, GSL_INTEG_GAUSS51, w, &result, &error);
    gsl_integration_workspace_free(w);
    return result;
}

double Halo::Rho(double R)
{
  // This is the gas density at radius R
  return RhoNorm * (pow(1.0 + gsl_pow_2(R / GasRs), - GasAlpha) * pow((1 + gsl_pow_2(R / GasRcool)), - (GasBeta - GasAlpha)) * pow((1 + gsl_pow_2(R / GasRcool2)), - (GasBeta2 - GasBeta))) ;
}

double MassShell(double R, void* p)
{
  // This is the mass in a shell of radius R
  // It is at top level to enable passing to gsl_integration
  Halo* Haloptr = (Halo*)p;
  return Haloptr->Rho(R) * 4.0 * M_PI * gsl_pow_2(R) * Haloptr->Pot->PPhi(R) * Haloptr->Pot->QPhi(R);
}

double Halo::Mass(double R)
{
    // This is the gas mass inside radius R
    gsl_function F;
    F.function=&MassShell;
    F.params=this;
    double result, error;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
    gsl_integration_qag(&F, 0.0, R, 0.0, 1.0E-3, 1000, GSL_INTEG_GAUSS51, w, &result, &error);
    gsl_integration_workspace_free(w);
    return result;
}

void Halo::SetRhoNorm()
{
  RhoNorm = 1.0 / Mass(C);
  RhoDMNorm = 1.0 / MassDM(C);
  return;
}

double UIntegrand(double R, void *p)
{
  // This is the integrand for calculating the internal energy
  // It is at top level to enable passing to gsl_integration
  Halo* Haloptr = (Halo*)p;
  return Haloptr->Rho(R) * gsl_spline_eval(Haloptr->Pot->f[0],R,Haloptr->Pot->facc[0]);
}

double Halo::UCalc(double R)
{
    // This calculates the internal energy at radius R
    gsl_function F;
    F.function=&UIntegrand;
    F.params=this;
    double result, error;
    double Rmax = 30.0 * C;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(100000);
    gsl_integration_qag(&F, R, Rmax, 0.0, 1.0E-5, 100000, GSL_INTEG_GAUSS51, w, &result, &error);
    gsl_integration_workspace_free(w);

    return 3.0 / 2.0 * result / Rho(R);
}

double Halo::U(double R)
{
  // This looks up U(R) from the look-up table
  if (R > 30.0 * C)
    return gsl_spline_eval(uinterp, 30.0 * C, uacc);
  else if (R < Pot->rmin)
    return gsl_spline_eval(uinterp, Pot->rmin, uacc);
  else
    return gsl_spline_eval(uinterp, R, uacc);
}

double Halo::NormR(double r, double m)
{
  return gsl_pow_2(Mass(r) - m);
}

double Halo::Rxyz(double* x)
// Function for determining R, given x,y,z, through iteration
{
  // First we translate the given coordinates back into the halo frame
  int i, j;
  double* xp = new double[3];

  for (i=0; i<3; i++)
    {
      xp[i] = 0.0;
      for (j=0; j<3; j++)
	{
          xp[i] += ((x[j] - Xoffset[j]) * EulerMinus[i][j]) / Rs;
        }
    }
  //Now we calculate R in the Halo frame
  double normf, r, newr, P, Q; 
  double tol = 1E-12; // Iteration tolerance
  int counter = 0;
  r = sqrt(xp[0]*xp[0] + xp[1]*xp[1] + xp[2]*xp[2]);
  normf = r * r;
  while ((normf > tol) && counter < 100)
	{
	  counter = counter + 1;
	  P = Pot->PPhi(r);
	  Q = Pot->QPhi(r);
	  newr = sqrt(xp[0]*xp[0] + xp[1]*xp[1] / (Q*Q) + xp[2]*xp[2] / (P*P));
	  normf = (newr - r) * (newr - r);
	  //printf("R = %f, NewR = %f, P = %f, Q = %f\n",r,newr,P,Q);
	  r = newr;
	}
  if (counter > 99)  printf("Number of iterations = %d, normf = %g\n", counter,normf);
  delete[] xp;
  return r; 
}

double Halo::Fntp(double R)
// This determines the non-thermal pressure fraction as a function of R
// Currently, this is just an estimate of the Battaglia function
{
  double FluidLossFactor = 0.6;// This is energy of fluid motion converted to thermal energy during the initial equilibration
  return Fntp0 * pow(R / C, Fntpexp) * FluidLossFactor; 
}

void Halo::ReadNbodyData()
{
  // This reads in the nbody data
  int i, j, counter = 0, test, mcounter = 0;
  double  rfactor, mfactor, vfactor, x[3], v[3], m, R, MaxR = 0.0;
  rfactor = Rs;
  mfactor =  M200 * (1.0 - GF);
  vfactor = sqrt(G * mfactor / rfactor);
  Part = new Particle[NDM]; // Particle data
  string line;
  vector<string> values;
  ifstream infile(nbodyfile.c_str());

  // Now we read the data from the file
  getline (infile,line); // Discard the first line, which is labels
  for (i=0; i<NDM; i++)
    {
      getline (infile,line);
      boost::split(values, line, boost::is_any_of(" ,\t"));
      test = (int) values.size();
      if (test == 7)
	{
	  for (j=0; j<3; j++)
	    {
	      x[j] = atof(values[j].c_str());
	      v[j] = atof(values[j + 3].c_str());
	      Part[counter].x[j] = (double)x[j] * rfactor;
	      Part[counter].v[j] = (double)v[j] * vfactor;
	    }
	  m = atof(values[6].c_str());
	  R = sqrt(gsl_pow_2(x[0]) + gsl_pow_2(x[1]/Q) + gsl_pow_2(x[2]/P));
	  if (R > MaxR) MaxR = R; // maximum radius
	  Part[counter].ID = counter + 1;
	  Part[counter].mass = m * mfactor;
	  Part[counter].grid = -1; // Used to identify which grid its in: -1 means unassigned
	  counter++;
	  if (R < C)
	    {
	      mcounter++;
	    }
	}
      else 
	{
	  printf("Some bad lines!! Exiting!! Line size = %d, i = %d\n",test,i);
	  exit(0);
	  continue;
	}
    }
  NDM200 = mcounter;
  MaxR = MaxR;
  printf("Successfully read in the nbody data from %s.\n",nbodyfile.c_str());
  printf("A total of %d DM particles were within R200.\n",mcounter);
  infile.close();
  return;
}

void Halo::EulerInitialize()
{
  int i;
  EulerPlus = new double*[3];
  EulerMinus = new double*[3];
  for (i=0; i<3; i++)
    {
      EulerPlus[i] = new double[3];
      EulerMinus[i] = new double[3];
    }

  EulerPlus[0][0] = -sin(Psi)*sin(Phi)+cos(Theta)*cos(Phi)*cos(Psi);
  EulerPlus[0][1] = -sin(Psi)*cos(Phi)-cos(Theta)*sin(Phi)*cos(Psi);
  EulerPlus[0][2] = cos(Psi)*sin(Theta);
  EulerPlus[1][0] = cos(Psi)*sin(Phi)+cos(Theta)*cos(Phi)*sin(Psi);
  EulerPlus[1][1] = cos(Psi)*cos(Phi)-cos(Theta)*sin(Phi)*sin(Psi);
  EulerPlus[1][2] = sin(Psi)*sin(Theta);
  EulerPlus[2][0] = -sin(Theta)*cos(Phi);
  EulerPlus[2][1] = sin(Theta)*sin(Phi);
  EulerPlus[2][2] = cos(Theta);

  EulerMinus[0][0] = -sin(-Phi)*sin(-Psi)+cos(-Theta)*cos(-Psi)*cos(-Phi);
  EulerMinus[0][1] = -sin(-Phi)*cos(-Psi)-cos(-Theta)*sin(-Psi)*cos(-Phi);
  EulerMinus[0][2] = cos(-Phi)*sin(-Theta);
  EulerMinus[1][0] = cos(-Phi)*sin(-Psi)+cos(-Theta)*cos(-Psi)*sin(-Phi);
  EulerMinus[1][1] = cos(-Phi)*cos(-Psi)-cos(-Theta)*sin(-Psi)*sin(-Phi);
  EulerMinus[1][2] = sin(-Phi)*sin(-Theta);
  EulerMinus[2][0] = -sin(-Theta)*cos(-Psi);
  EulerMinus[2][1] = sin(-Theta)*sin(-Psi);
  EulerMinus[2][2] = cos(-Theta);

  return;
}

void Halo::EulerRotateShift(int n)
{
  //This rotates and shifts the dark matter positions and velocities
  int i, j ;
  double* xp = new double[3];
  double* vp = new double[3];

  for (i=0; i<3; i++)
    {
      xp[i] = 0.0;
      for (j=0; j<3; j++)
	{
          xp[i] += Part[n].x[j] * EulerPlus[i][j];
          vp[i] += Part[n].v[j] * EulerPlus[i][j];
        }
    }
  for (i=0; i<3; i++)
    {
      Part[n].x[i] = xp[i] + Xoffset[i];
      Part[n].v[i] = vp[i] + Voffset[i];
    }
  delete[] xp;
  delete[] vp;
  return;
}



