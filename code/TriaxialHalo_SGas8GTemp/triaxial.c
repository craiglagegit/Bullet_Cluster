/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, NYU CCPP
  Date: Aug 21, 2012

  Triaxial halo generation program - takes an input of DM halo data from SMILE and adds gas.

*/

#include "triaxial.h"

//***************MAIN PROGRAM***********************

int main(int argc, char *argv[])
{
  if(argc != 2)
    {
      printf("\nwrong number of arguments\n");
      printf("Only argument should be name of .cfg file");
      exit(0);
    }
  int i,io,j,seed;
  seed = (int) time(NULL);
  srand48(seed);
  struct Halo H = ReadConfig(argv[1]);  //Read configuration file and initialize parameters
  printf("P=%f, Q=%f, G = %f, R200 = %f, Rs = %f, M200 = %f NDM = %d, Ngas = %d, NTot = %d\n",H.P,H.Q,H.G,H.R200,H.Rs,H.M200,H.NDM,H.NGas,H.NTot);
  struct Potential Pot = ReadPotentialData(H.potentialfile, &H);  //Read Potential data and create interpolation objects.
  struct PQPhi PQ = PQLookup(&Pot, &H);  //Create lookup table for Potential shape parameters Pphi(R), and Qphi(R).
  SetRhoNorm(&PQ,&H); // Normalize density profiles
  struct UMass UM = UMLookup(&Pot,&PQ,&H);  //Create lookup table for R(m), and U(R).
  /*
  double xx;
  for (i = 0; i<100; i++)
  {
    xx = 1.0E-6 * pow(pow(10.0,.10),(double)i);
    printf("R = %f, Pphi(R) = %f, Qphi(R) = %f AEllipse = %f, Rho(R) = %g, Mass(R) = %f\n",xx,Pphi(xx,&PQ,&H), Qphi(xx,&PQ,&H),AEllipse(xx,&PQ,&H),Rho(xx,&H),Mass(xx,&PQ,&H));
  }
  exit(0);
  */
  double R,m,x,y,z,r,theta,phi,vx,vy,vz, ufactor, mtot, Pgas, Qgas;
  struct GadgetParticles *PTemp = (struct GadgetParticles*) calloc(H.NTot, sizeof(struct GadgetParticles));  //Particle data.
  ReadNbodyData(H.nbodyfile, PTemp, &H); // This reads in the halo data generated by SMILE
  struct GadgetParticles *P = (struct GadgetParticles*) calloc(H.NTot, sizeof(struct GadgetParticles));  //Particle data.
  struct GadgetHeader Header = CreateGadgetHeader(&H); 
  mtot = (double) H.NDM / (double) H.NDM200;
  printf("Maximum DM Particle radius is %f * R200. Total mass is %f * M200\n",H.MaxR,mtot);
  fflush(stdout);

  //ufactor = H.G * H.M200 * (1.0 - H.GF) / H.Rs;
  ufactor = H.G * H.M200 / H.Rs;
  int mcounter = 0; // Number of particles within R200
  // Now we add the gas particle data.
  double R500, R500_15, Ave_Temp, mstar, tfactor, ucalc;
  int nTemp = 0;
  mstar = 2.15E-24; // Average ion mass in CGS
  tfactor = 2.0 * (mstar / UMASS) * H.G * H.M200 / (3.0 * H.Rs) / keV;

  R500 = gsl_spline_eval(UM.rinterp,0.4,UM.racc);  // R200 is m = 1, R500 is m=2/5
  R500_15 = R500*0.15;
  for (i=0; i<H.NGas; i++) 
  {
	P[i].iParticleID = i + 1;	
	m = 1.5 * mtot * drand48(); // Larger mass ratio than DM, since gas falls off more slowly
	R = gsl_spline_eval(UM.rinterp,m,UM.racc);  // R specifies which ellipsoid we're on
	if (R < H.C)  mcounter = mcounter +1;
  	Pgas = Pphi(R,&PQ,&H); // These are the shape parameters of the potential ellipsoids 
  	Qgas = Qphi(R,&PQ,&H); // These are the shape parameters of the potential ellipsoids 
	RandomDirection(R,&theta,&phi,&x,&y,&z,Pgas,Qgas); // Choose a random (theta, phi)

	if (i%100000 == 0)
	  {
	    printf("m = %f, R = %f, Pgas = %f, Qgas = %f, PDM = %f, QDM = %f\n",m, R, Pgas, Qgas, H.P, H.Q);
	  }

	r = sqrt(gsl_pow_2(x) + gsl_pow_2(y) + gsl_pow_2(z));
	P[i].rOriginDistance = H.Rs * r; 
	P[i].rAxisDistance = H.Rs * sqrt(gsl_pow_2(x) + gsl_pow_2(y));
	P[i].rLongitude = theta;
	P[i].rLatitude = phi;
	P[i].rPosition[0] = H.Rs * x;
	P[i].rPosition[1] = H.Rs * y;
	P[i].rPosition[2] = H.Rs * z;
	P[i].rVelocity[0] = 0.0;
	P[i].rVelocity[1] = 0.0;
	P[i].rVelocity[2] = 0.0;
	P[i].iParticleType = 0;
	ucalc = gsl_spline_eval(UM.uinterp,R,UM.uacc);
	P[i].rU = ucalc * ufactor;
	if (R > R500_15 && R < R500)
	  {
	    Ave_Temp = Ave_Temp + ucalc;
	    nTemp = nTemp + 1;
	  }
	P[i].rRho = Rho(R,&H);// Is this right?
	P[i].rNe = 0.0;
	if (i % 100000 == 0)   printf("Particle %d, ID = %d, x = %f, y = %f, z = %f, vx = %f, vy = %f, vz = %f, u = %f, mass = %f\n",i,P[i].iParticleID,P[i].rPosition[0],P[i].rPosition[1],P[i].rPosition[2],P[i].rVelocity[0],P[i].rVelocity[1],P[i].rVelocity[2],P[i].rU,P[i].rMass);
	fflush(stdout);
  }
  Ave_Temp = Ave_Temp / (double) nTemp * tfactor;
  printf("A total of %d gas particles were between 0.15 R500 and R500, with an Average Temp of %f keV.\n",nTemp,Ave_Temp);
  // Now set the gas particle mass, based on the number of particles within R200
  double gmass = (H.M200 * H.GF) / (double) mcounter;
  for (i=0; i<H.NGas; i++) 
  {
    P[i].rMass = gmass;
  }
  printf("A total of %d gas particles were within R200.\n",mcounter);
  // Next we transfer the DM data from PTemp to P
  for (i=H.NGas; i<H.NTot; i++) 
  {
	io = i - H.NGas;
	P[i].iParticleID = i+1;	
	P[i].iParticleType = 1;	
	P[i].rMass = PTemp[io].rMass;
			for (j=0; j<3; j++)
			{
				P[i].rPosition[j] = PTemp[io].rPosition[j];
				P[i].rVelocity[j] = PTemp[io].rVelocity[j];
			}
	if (i % 100000 == 0)   printf("Particle %d, ID = %d, x = %f, y = %f, z = %f, vx = %f, vy = %f, vz = %f, u = %f, mass = %f\n",i,P[i].iParticleID,P[i].rPosition[0],P[i].rPosition[1],P[i].rPosition[2],P[i].rVelocity[0],P[i].rVelocity[1],P[i].rVelocity[2],P[i].rU,P[i].rMass);
	fflush(stdout);
	
  }  
  //CleanUp(P,&H);
  WriteProfiles(H.profile, &H, &UM, &PQ, &Pot);  
  WriteGadgetFile(H.outfile, Header, P);
  //CharacterizeDistribution(P,G,&H);*/

  free(P);
  free(PTemp);

  return 0;
}
//***************END MAIN PROGRAM***********************
