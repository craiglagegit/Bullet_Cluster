//****************** fileio.c **************************

#include "triaxial.h"

struct Halo ReadConfig(const char *inname)
// Reads in the configuration file and sets all of the various constants
{
  struct Halo H;
  H.G = GRAVITY * ((UMASS * gsl_pow_2(UTIME)) / (ULENGTH * ULENGTH * ULENGTH));
  H.H = HUBBLE * ((100 * 1e5 * ULENGTH) / (CM_PER_MPC * UVELOCITY));

  H.P = GetDoubleParam(inname, "P", 0.6);
  H.Q = GetDoubleParam(inname, "Q", 0.8);
  H.P = GSL_MIN(H.P, 1.00);
  H.Q = GSL_MIN(H.Q, 1.00);
  H.alpha = - 1.0;
  H.beta =  - gsl_pow_2(H.Q);
  H.gamma = - gsl_pow_2(H.P);
  H.NPot = GetIntParam(inname, "NPot", 100);
  H.NLook = GetIntParam(inname, "NLook", 100);
  H.MaxProb = MaxProb(H.alpha, H.beta, H.gamma);

  H.M200 = GetDoubleParam(inname, "M200", 10000.0);
  H.C = GetDoubleParam(inname, "C", 3.0);
  //H.rmin = GetDoubleParam(inname, "rmin", 0.001);
  H.R200 = pow((H.M200 * H.G) / (100.0 * gsl_pow_2(H.H)), 1.0 / 3.0);
  H.Rs = H.R200 / H.C;
  H.rmin = 1.0E-6; //Min value in lookup tables
  H.mmin = 1.0E-6; //Min value in lookup tables

  H.RhoNorm = 1.0;
  H.RhoDMNorm = 1.0;
  H.GasAlpha = GetDoubleParam(inname, "GasAlpha", 0.0);
  H.GasBeta = GetDoubleParam(inname, "GasBeta", 0.5);
  H.GasRcool = GetDoubleParam(inname, "GasRcool", 30.0);
  H.GasEpsilon = GetDoubleParam(inname, "GasEpsilon", 0.0);
  H.GasRs = GetDoubleParam(inname, "GasRs", 1.0);
  H.GasBeta2 = GetDoubleParam(inname, "GasBeta2", 2.0);
  H.GasRcool2 = GetDoubleParam(inname, "GasRcool2", 10.0);
  H.GasN2 = GetDoubleParam(inname, "GasN2", 100.0);
  H.GasRcool = H.GasRcool / H.Rs;
  H.GasRs = H.GasRs / H.Rs;
  H.GasRcool2 = H.GasRcool2 / H.Rs;

  H.NGas = GetIntParam(inname, "NGas", 10000);
  H.NDM = GetIntParam(inname, "NDM", 10000);
  H.NTot = H.NDM + H.NGas;
  H.GF = GetDoubleParam(inname, "GF", 0.16);
  GetStringParam(inname, "outfile", "output.dat", H.outfile);
  GetStringParam(inname, "profile", "profiles.dat", H.profile);
  GetStringParam(inname, "potentialfile", "potential", H.potentialfile);
  GetStringParam(inname, "nbodyfile", "nbody", H.nbodyfile);
  
  return H;
}


struct Potential ReadPotentialData(const char *inname, struct Halo *H)
{
  // This reads in the potential data and initializes interpolation objects
  struct Potential Pot;
  char dummy[128];
  FILE *infile = fopen(inname, "r");// This is the input file
  int i, j;
  float x[3], phi, fx[3], rho;

  for (i=0; i<3; i++)
  {
	Pot.xa[i]= calloc(H->NPot, sizeof(double));
	Pot.phia[i]= calloc(H->NPot, sizeof(double));
	Pot.fa[i]= calloc(H->NPot, sizeof(double));
	//printf("In ReadPotential\n");
	//fflush(stdout);
  	fgets(dummy,128,infile);// Read a dummy line
	for (j=0; j<H->NPot; j++)
	{
		if (fscanf(infile, "%f %f %f %f %f %f %f %f",  &x[0], &x[1], &x[2], &phi, &fx[0], &fx[1], &fx[2], &rho)==8)
		{
			Pot.xa[i][j] = (double)x[i];
			Pot.phia[i][j] = (double)phi;
			Pot.fa[i][j] = (double)fx[i];
			//printf("i = %d, j = %d, x = %f, phi = %f, fx = %f, rho = %f\n",i,j,Pot.xa[i][j],Pot.phia[i][j],Pot.fa[i][j],rho);
		}
		else continue;
	}
	//fflush(stdout);
        Pot.phiacc[i] = gsl_interp_accel_alloc();
        Pot.phi[i] = gsl_spline_alloc (gsl_interp_cspline, H->NPot);  
        gsl_spline_init (Pot.phi[i],  Pot.xa[i], Pot.phia[i], H->NPot);
        Pot.facc[i] = gsl_interp_accel_alloc();
        Pot.f[i] = gsl_spline_alloc (gsl_interp_cspline, H->NPot);  
        gsl_spline_init (Pot.f[i],  Pot.xa[i], Pot.fa[i], H->NPot);

  }
  fclose(infile);
  if (i != 3 || j != H->NPot)
  {
	printf("Error reading Potential file! Exiting.\n");
	exit(0);
  }
  return Pot;
}

void ReadNbodyData(const char *inname, struct GadgetParticles *P, struct Halo *H)
{
  // This reads in the nbody data
  char dummy[128];
  FILE *infile = fopen(inname, "r");// This is the input file
  fgets(dummy,128,infile);// Read a dummy line
  int i, j, test, counter = 0, mcounter = 0;
  float x[3], v[3], m, R;
  double  rfactor, mfactor, vfactor, MaxR = 0.0;
  rfactor = H->Rs;
  mfactor =  H->M200 * (1.0 - H->GF);
  vfactor = sqrt(H->G * mfactor / rfactor);
  for (i=H->NGas; i<H->NTot; i++)
  {
	test = fscanf(infile, "%f %f %f %f %f %f %f\n",  &x[0], &x[1], &x[2], &v[0], &v[1], &v[2], &m);
	if (test == 7)
	{
		R = sqrt(gsl_pow_2(x[0]) + gsl_pow_2(x[1]/H->Q) + gsl_pow_2(x[2]/H->P));
		if (R > MaxR) MaxR = R; // maximum radius
		for (j=0; j<3; j++)
		  {
		    P[counter].rPosition[j] = (double)x[j] * rfactor;
		    P[counter].rVelocity[j] = (double)v[j] * vfactor;
		  }
		P[counter].iParticleID = counter + 1;
		P[counter].iParticleType = 1;
		P[counter].rMass = m * mfactor;
		counter = counter + 1;
		if (R < H->C)
			{
			  mcounter = mcounter + 1;
			}
	}
	else 
	{
		printf("Some bad lines!! Exiting!! Test = %d, i = %d\n",test,i);
		exit(0);
		continue;
	}
  }
  H->NDM200 = mcounter;
  H->MaxR = MaxR;
  printf("A total of %d DM particles were within R200.\n",mcounter);
  fclose(infile);
  return;
}


double GetDoubleParam(const char *fname, const char *parnam, double defval)
{
  char parval[STRING_LENGTH];
  if (ReadPar(fname, parnam, parval))
    {
      printf("Parameter %s not found in %s. Using default: %f.\n", parnam, fname, defval);
      return defval;
    }
  else
    {
      return atof(parval);
    }
}

int GetIntParam(const char *fname, const char *parnam, int defval)
{
  char parval[STRING_LENGTH];
  if (ReadPar(fname, parnam, parval))
    {
      printf("Parameter %s not found in %s. Using default: %d.\n", parnam, fname, defval);
      return defval;
    }
  else
    {
      return atoi(parval);
    }
}

void GetStringParam(const char *fname, const char *parnam, char *defval, char* result)
{
  char parval[STRING_LENGTH];
  if (ReadPar(fname, parnam, parval))
    {
      printf("Parameter %s not found in %s. Using default: %s.\n", parnam, fname, defval);
      strcpy(result, defval);
    }
  else
    {
      strcpy(result, parval);
    }
  return;
}


int ReadPar(const char *fname, const char *parnam, char *parval)
{
  FILE *pfile = fopen(fname, "r");
  char *charp;
  char  line[STRING_LENGTH];
  char  lhs [STRING_LENGTH];

  while ((charp = fgets(line, STRING_LENGTH, pfile)))
    {

      // skip lines beginning with #
      if ((charp = strchr(line, '#')))
	*charp = '\0';

      // break the line into variable name and value 
      if (strchr(line, '='))
	{
	  int i;
	  char *cp0, *cp1;
	  cp0 = strtok(line, "=");
	  cp1 = strtok(NULL, "=");

	  for (i=0; i<STRING_LENGTH; ++i)
	    {
	      if ((charp=strchr(cp0, '\t'))) *charp = ' ';
	      else break;
	    }

	  cp0 = strtok(cp0, " ");

	  for (i=0; i<STRING_LENGTH; ++i)
	    {
	      if ((charp = strchr(cp1, '\t'))) *charp = ' ';
	      else break;
	    }
	  for (i=0; i<STRING_LENGTH; ++i)
	    {
	      if ((charp = strchr(cp1, '"'))) *charp = ' ';
	      else break;
	    }
	  for (i=0; i<STRING_LENGTH; ++i)
	    {
	      if ((charp = strchr(cp1, '\n'))) *charp = ' ';
	      else break;
	    }

	  cp1 = strtok(cp1, " ");
	  strcpy(lhs, cp0);

	  if (strcmp(lhs, parnam) == 0)
	    {
	      strcpy(parval, cp1);
	      fclose(pfile);
	      return 0;
	    }
	}
    }

  fclose(pfile);
  return 1;
}

struct GadgetHeader CreateGadgetHeader(struct Halo *H)
{
  struct GadgetHeader Header;
  int i;
  // First, zero things out
  for (i=0; i<6; i++)
  {
	Header.iNumberOfParticles[i] = 0;
	Header.iTotalNumberOfParticles[i] = 0;
	Header.rParticleMass[i] = 0;
	Header.NallHW[i] = 0;

  }
  Header.iNumberOfParticles[0] = H->NGas;
  Header.iNumberOfParticles[1] = H->NDM;
  Header.iTotalNumberOfParticles[0] = H->NGas;
  Header.iTotalNumberOfParticles[1] = H->NDM;
// Most of the parameters below are not used, and are over-written by the .param file when Gadget is run.
  Header.rExpansionFactor = 1.0;       // ( Time / a )
  Header.rRedShift = 0.0;              // Set only for cosmological simulations

  Header.iStarFormationFlag = 0;
  Header.iFeedbackFlag = 0;
  Header.iCoolingFlag = 0;

  Header.iNumberOfFiles = 1;
  Header.rBoxSize = 0.0;                    // Relevant only for periodic boundaries
  Header.rOmega0 = 0.0;                     // matter density at z=0 (cos)
  Header.rOmegaLambda = 0.0;                // Vacuum energy density at z=0 (cos)
  Header.rHubbleConstant = 0.0;             // Value of the hubble constant in units of 100 km/Mpc*sec
  Header.FlagAge = 0;
  Header.FlagMetals = 0;
  Header.flag_entr_ics = 0;
  for (i=0; i<60; i++) Header.bUnused[i]=0;

  return Header;
}

int WriteGadgetFile(char *filename, struct GadgetHeader Header, struct GadgetParticles *Particles)
{
  //Borrowed from SphericalHalo
  //Written 2001 by Martin Jubelgas,
  //Max-Planck-Institute for Astronomy, Garching

  FILE *GadgetFile;		// Handle for Data File

  int iCounter;			// used for loops through Particle Field
  int iTypeCounter;		// used for loops through Particle Types
  int iCounterStart;		// Offset for first particle of current type
  int iTotalParticles;		// Total number of Particles
  int iWriteNumber;		// Used for verification of Write operations
  int iInvert;			// Boolean, Boolean Conversion required == 1
  int iFirstInt;
  float *rBuffer;		// Buffer for read operations
  int *iBuffer;
  void *ReallocAdress;
  int iSize;

  // Open the Gadget File
  GadgetFile = fopen(filename, "wb");
  if(GadgetFile == NULL)
    {
      // If File could not be opened, return with error
      return COULD_NOT_OPEN_FILE;
    }
  // Determine number of particles
  iTotalParticles = 0;
  for(iCounter = 0; iCounter < 6; iCounter++)
    {
      iTotalParticles += Header.iNumberOfParticles[iCounter];
    }
  // Reserve Buffer for float - Values, like Coordinates, Velocities, ...
  rBuffer = calloc(iTotalParticles, 3 * sizeof(float));
  if(rBuffer == NULL)
    {
      // If not enough memory available, dealloc memory for Particles
      //   and return with an error
      fclose(GadgetFile);
      return INSUFFICIENT_MEMORY;
    }
  iBuffer = (int *) rBuffer;
  iSize = 256;
  fwrite(&iSize, sizeof(int), 1, GadgetFile);
  iWriteNumber = fwrite(&Header, sizeof(struct GadgetHeader), 1, GadgetFile);
  if(iWriteNumber != 1)
    {
      // If header could not be read properly, return with error
      fclose(GadgetFile);
      return WRITE_ERROR;
    }
  fwrite(&iSize, sizeof(int), 1, GadgetFile);
  // START WRITING PARTICLE DATA
  // +--> Positions (all Particle Types)
  iSize = 3 * sizeof(int) * iTotalParticles;
  fwrite(&iSize, sizeof(int), 1, GadgetFile);
  printf("Storing Positions...\n");
  fflush(stdout);
  for(iCounter = 0; iCounter < iTotalParticles; iCounter++)
    {
      rBuffer[(3 * iCounter)] = Particles[iCounter].rPosition[0];
      rBuffer[(3 * iCounter) + 1] = Particles[iCounter].rPosition[1];
      rBuffer[(3 * iCounter) + 2] = Particles[iCounter].rPosition[2];
    }
  printf("Done.\n");
  fflush(stdout);
  fwrite(rBuffer, 3 * sizeof(float), iTotalParticles, GadgetFile);
  fwrite(&iSize, sizeof(int), 1, GadgetFile);
  // +--> Velocities (all Particle Types)
  fwrite(&iSize, sizeof(int), 1, GadgetFile);
  for(iCounter = 0; iCounter < iTotalParticles; iCounter++)
    {
      rBuffer[(3 * iCounter)] = Particles[iCounter].rVelocity[0];
      rBuffer[(3 * iCounter) + 1] = Particles[iCounter].rVelocity[1];
      rBuffer[(3 * iCounter) + 2] = Particles[iCounter].rVelocity[2];
    }

  fwrite(rBuffer, 3 * sizeof(float), iTotalParticles, GadgetFile);
  fwrite(&iSize, sizeof(int), 1, GadgetFile);
  // +--> Particle ID's
  iSize = iTotalParticles * sizeof(int);
  fwrite(&iSize, sizeof(int), 1, GadgetFile);
  for(iCounter = 0; iCounter < iTotalParticles; iCounter++)
    {
      iBuffer[iCounter] = Particles[iCounter].iParticleID;
    }
  fwrite(iBuffer, sizeof(int), iTotalParticles, GadgetFile);
  fwrite(&iSize, sizeof(int), 1, GadgetFile);
  // +--> Masses (ONLY for Particle Types with variable Mass,
  //               defined in header with Mass = 0)
  // RIGHT NOW: Only Mass = 0 in header allowed!!!
  iCounterStart = 0;
  fwrite(&iSize, sizeof(int), 1, GadgetFile);
  for(iCounter = 0; iCounter < iTotalParticles; iCounter++)
    {
      rBuffer[iCounter] = Particles[iCounter].rMass;
    }
  fwrite(rBuffer, sizeof(float), iTotalParticles, GadgetFile);
  fwrite(&iSize, sizeof(int), 1, GadgetFile);
  // +--> Internal Energy (GAS ONLY= TYPE 0)
  iTotalParticles = Header.iNumberOfParticles[0];	// We're interested in
  // gas particles only from
  // now on
  iSize = sizeof(int) * iTotalParticles;
  if(iTotalParticles != 0)
    {
      for(iCounter = 0; iCounter < iTotalParticles; iCounter++)
	{
	  rBuffer[iCounter] = Particles[iCounter].rU;
	}
      fwrite(&iSize, sizeof(int), 1, GadgetFile);
      fwrite(rBuffer, sizeof(float), iTotalParticles, GadgetFile);
      fwrite(&iSize, sizeof(int), 1, GadgetFile);
      // WE'RE FINISHED! RETURN SUCCESS !!
    }
  free(rBuffer);
  fclose(GadgetFile);
  return OK;
}

void WriteProfiles(char *filename,struct Halo *H, struct UMass *UM, struct Potential *Pot)
// This writes out the density and temperature profiles to a file
// Units are g/cm^3 and keV
// Logarithmic steps in radius are used
{
  int i, npoints;
  double rhodm, rhogas, mdm, mgas, t, fgrav, r, rmin, rmax, dr, tfactor, mstar;
  mstar = 2.15E-24; // Average ion mass in CGS
  npoints = 256;
  rmin = .0001; // min r/Rs
  rmax = 3.0 * H->C; // max r/Rs
  dr = pow(rmax / rmin, 1.0 / (npoints - 1));
  //tfactor = 2.0 * (mstar / UMASS) * H->G * H->M200 * (1.0 - H->GF) / (3.0 * H->Rs) / keV;
  tfactor = 2.0 * (mstar / UMASS) * H->G * H->M200 / (3.0 * H->Rs) / keV;

  FILE *outfile = fopen(filename, "w");//This is the output file
  fprintf(outfile, "# R     RhoDM      RhoGas     MDM      MGas      TGas      FGrav\n");
  for (i=0; i<npoints; i++)
    {
      r = rmin * pow(dr, i);
      rhodm = H->M200 * (1.0 - H->GF) * RhoDM(r, H) / gsl_pow_3(H->Rs) * UDENSITY;
      rhogas = H->M200 * H->GF * Rho(r,H) / gsl_pow_3(H->Rs) * UDENSITY;
      mdm = H->M200 * (1.0 - H->GF) * MassDM(r,H) * UMASS;
      mgas = H->M200 * H->GF * Mass(r,H) * UMASS;
      t = gsl_spline_eval(UM->uinterp,r,UM->uacc) * tfactor;
      fgrav = gsl_spline_eval(Pot->f[0],r,Pot->facc[0]);
      fprintf(outfile, "%.4g %.4g %.4g %.4g %.4g %.4g %.4g\n",r,rhodm,rhogas,mdm,mgas,t,fgrav);
      fflush(outfile);
    }
  fclose(outfile);
  return;
}
