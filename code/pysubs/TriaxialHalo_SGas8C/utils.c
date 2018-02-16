//****************** utils.c **************************

#include "triaxial.h"


double ERatio(double R)
     // This is an approximation of the ratio of the potential ellipsoid to the halo ellipsoid
     // See Lee and Suto - 2002
{
  if (R < 1.0E-3) return 0.5;
  else return (-6.0 * R - 3.0 * gsl_pow_2(R) + gsl_pow_3(R) + 6.0 * (1.0  + R) * log(1.0 + R)) / (2.0 * gsl_pow_2(R) * (-R + (1.0 + R) * log(1.0 + R)));
}

double Pphi(double R, struct Halo *H)
     // This gives p of the potential given p of the halo
{

  return sqrt(1.0 + (gsl_pow_2(H->P) - 1.0) * ERatio(R));
}

double Qphi(double R, struct Halo *H)
     // This gives q of the potential given q of the halo
{

  return sqrt(1.0 + (gsl_pow_2(H->Q) - 1.0) * ERatio(R));
}

double AEllipse(double R, struct Halo *H)
     // This gives the surface area of the potential ellipsoid
     // It is an approximation good to better than 1%
{
  double px, qx, x = 1.6075;
  px = pow(Pphi(R,H),x);
  qx = pow(Qphi(R,H),x);
 
  return pow((px + qx + px * qx) / 3.0, 1 / x);
}

double AEllipseDM(double R, struct Halo *H)
     // This gives the surface area of the potential ellipsoid
     // It is an approximation good to better than 1%
{
  double px, qx, x = 1.6075;
  px = pow(H->P,x);
  qx = pow(H->Q,x);
 
  return pow((px + qx + px * qx) / 3.0, 1 / x);
}

int RandomSign() 
{
    return (rand() % 2) * 2 - 1;
}

double GaussianRandom(double mean, double sigma)
{
  int i;
  double sum = 0.0;
  for(i=0; i<12; i++)
    {
      sum = sum + drand48();
    }
  sum = sum - 6.0;
  return (mean + (sigma * sum));
}

double RetainProb(double x, double phi, double alpha, double beta, double gamma)
   // This calculates the probability of retaining an (x,phi) pair
   // based on the (theta, phi) Area element on the ellipsoid.  Note x = cos(theta).
{
  double prob;

  prob = sqrt((gsl_pow_2(alpha)*gsl_pow_2(beta)*gsl_pow_2(gamma)*(gsl_pow_2(alpha)*gsl_pow_2(beta)*gsl_pow_2(x) - 
        (gsl_pow_2(gamma)*(-1 + gsl_pow_2(x))*gsl_pow_2(alpha + beta + (-alpha + beta)*cos(2*phi)))/4.)*
      (1 - (gsl_pow_2(alpha - beta)*gsl_pow_2(gamma)*(-1 + gsl_pow_2(x))*gsl_pow_2(cos(phi))*gsl_pow_2(sin(phi)))/
         gsl_pow_2(alpha*beta*gsl_pow_2(x) - gamma*(-1 + gsl_pow_2(x))*(beta*gsl_pow_2(cos(phi)) + alpha*gsl_pow_2(sin(phi))))))/
    gsl_pow_2(gsl_pow_2(alpha*beta*gsl_pow_2(x) - gamma*(-1 + gsl_pow_2(x))*(beta*gsl_pow_2(cos(phi)) + alpha*gsl_pow_2(sin(phi))))));

  return prob;
}

void RandomDirection(double R, double *theta, double *phi, double *x, double *y, double *z, double alpha, double beta, double gamma, struct Halo *H)
   // Function for determining (x,y,z) point on an ellipsoidal surface,
   // randomly distributed on the surface, given a value for the ellipsoid R.
   // It is a rejection sampling algorithm where an (X, phi) pair is generated and then 
   // samples are discarded until they meet the desired distribution.  Then (theta=arccos(X),phi)
   // is transformed to (x,y,z)
{
  double X, r, fthetaphi;
  int i;
  //maxprob = MaxProb(alpha, beta, gamma);
  for(i=0; i<1000; i++)
  {
	X = -1.0 + 2.0 * drand48();
	*phi = 2.0 * pi * drand48();
	if (RetainProb(X, *phi, alpha, beta, gamma) / H->MaxProb  > drand48()) break;
  }
 if (i > 999) 
    {
	printf("Count iteration exceeded in Direction generation.  i = %d.  Something's wrong!\n",i);
    }

  *theta = acos(X);

  fthetaphi = -(gsl_pow_2(cos(*phi))*gsl_pow_2(sin(*theta)))/alpha - (gsl_pow_2(sin(*phi))*gsl_pow_2(sin(*theta)))/beta - gsl_pow_2(cos(*theta))/gamma;
  r =  R / sqrt(fthetaphi);
  *x = r * sin(*theta) * cos(*phi);
  *y = r * sin(*theta) * sin(*phi);
  *z = r * cos(*theta);

  return; 
}

double MaxProb(double alpha, double beta, double gamma)
   // Function for determining the maximum probability
   // of the (x, phi) function to use in normalizing RetainProb.
   // It is only a function of ellipsoid shape.
{
  double x, phi, prob, maxprob = 0.0;
  phi = pi;  // The function max is near here.
  int i;
  for(i=0; i<100; i++)
  {
	x = 0.25 * drand48(); // The function max is in this region
	prob = RetainProb(x, phi, alpha, beta, gamma);
	if (prob  > maxprob) maxprob = prob;
  }
  return 1.2 * maxprob; // 1.2 is just an increase for good measure in case we missed the maximum slightly
}

/*
void CharacterizeDistribution(struct GadgetParticles *P, double *G, struct Halo *H)
   // Function for calculating sigma v's,
   // used for debug.
{
  int i;
  double vx, vy, vz, variancevx=0.0, variancevy=0.0, variancevz=0.0;
  double meanvx = 0.0, meanvy = 0.0, meanvz = 0.0;
  double sigmavx, sigmavy, sigmavz;
  double x, y, z, lambda, mu, nu, lx, ly, lz, lmag;
  double KE, PE, PElam, PEmu, PEnu, TotalKE = 0.0, TotalPE = 0.0;
  double vfactor, Eratio;
  vfactor = sqrt(H->M200 / H->Rs);

  for (i=0; i<H->NTot; i++) 
  {
	  vx = P[i].rVelocity[0];
	  vy = P[i].rVelocity[1];
	  vz = P[i].rVelocity[2];

	  x = P[i].rPosition[0] / H->Rs;
	  y = P[i].rPosition[1] / H->Rs;
	  z = P[i].rPosition[2] / H->Rs;

	  lx = y * vz - z * vy;
	  ly = z * vx - x * vz;
	  lz = x * vy - y * vx;
	  lmag = sqrt (lx*lx + ly*ly + lz*lz);
	 
	  meanvx = meanvx + vx;
	  meanvy = meanvy + vy;
	  meanvz = meanvz + vz;

	  variancevx = variancevx + gsl_pow_2(vx);
	  variancevy = variancevy + gsl_pow_2(vy);
	  variancevz = variancevz + gsl_pow_2(vz);

	  KE = (gsl_pow_2(vx) + gsl_pow_2(vy) + gsl_pow_2(vz)) / (2.0 * gsl_pow_2(vfactor));
	  LMN(x,y,z,&lambda,&mu,&nu,H);
	  PElam = -(lambda + H->alpha) * (lambda + H->gamma) / ((lambda - mu) * (lambda - nu)) * GLookup(lambda,G,H);
	  PEmu =  - (mu + H->alpha) * (mu + H->gamma) / ((mu - nu) * (mu - lambda)) * GLookup(mu,G,H);
	  PEnu =  - (nu + H->alpha) * (nu + H->gamma) / ((nu - lambda) * (nu - mu)) * GLookup(nu,G,H);
	  PE = PElam + PEmu + PEnu;
	  TotalKE = TotalKE + KE;
	  TotalPE = TotalPE + PE;
	  //if (i%10000 == 0)printf(" x = %f, y = %f, z = %f, lambda = %f, mu = %f, nu = %f\n",x,y,z,lambda,mu,nu);
	  //if (i%10000 == 0)printf(" KE = %f, PE = %f\n",KE, PE);
  }
  meanvx = meanvx / H->NTot;
  meanvy = meanvy / H->NTot;
  meanvz = meanvz / H->NTot;

  sigmavx = sqrt(variancevx / H->NTot - gsl_pow_2(meanvx));
  sigmavy = sqrt(variancevy / H->NTot - gsl_pow_2(meanvy));
  sigmavz = sqrt(variancevz / H->NTot - gsl_pow_2(meanvz));
  Eratio = fabs(TotalPE / TotalKE);
  printf("Meanvx = %f, Meanvy = %f, Meanvz = %f, Sigmavx = %f, Sigmavy = %f, Sigmavz = %f\n",meanvx,meanvy,meanvz,sigmavx,sigmavy,sigmavz);
  printf( "TotalKE = %g, TotalPE = %g, ERatio = %f\n",TotalKE,TotalPE,Eratio);
  printf( "Angular Momentum: Lx = %g, Ly = %g, Lz = %g, Lmag = %g\n",lx,ly,lz,lmag);
  fflush(stdout);

  return;
}


void Symmetrize(struct GadgetParticles *P, struct Halo *H)
   // Randomizes the particles about any axis of symmetry
   // This improves the uniformity in symmetrical situations
{
  if (H->P >= 0.989999) // Oblate case - symmetrical about z-axis
  {
	Randomize(2,P,H);
  }
  if (H->Q >= 0.989999) // Prolate case - symmetrical about x-axis
  {
	Randomize(0,P,H);
  }

  if (H->P >= 0.989999 && H->Q >= 0.989999)  // Spherical case, symmetrical about all axes.
  {
	Randomize(1,P,H);
	Randomize(0,P,H);
	Randomize(2,P,H);
  }
  return;
}

void Randomize(int axis, struct GadgetParticles *P, struct Halo *H)
   // Randomizes the particles about an axis
   // This improves the uniformity in symmetrical situations
{
  int i, axis1, axis2;
  double theta, costheta, sintheta;
  axis1 = (axis + 1) % 3; // One axis perpendicular to the given axis
  axis2 = (axis + 2) % 3; // Second axis perpendicular to the given axis
  for (i=0; i<H->NTot; i++) 
  {
	  theta = 2.0 * pi * drand48();
	  costheta = cos(theta);
	  sintheta = sin(theta);
	  P[i].rPosition[axis1] = P[i].rPosition[axis1] * costheta + P[i].rPosition[axis2] * sintheta;
	  P[i].rPosition[axis2] = -P[i].rPosition[axis1] * sintheta + P[i].rPosition[axis2] * costheta;
	  P[i].rVelocity[axis1] = P[i].rVelocity[axis1] * costheta + P[i].rVelocity[axis2] * sintheta;
	  P[i].rVelocity[axis2] = -P[i].rVelocity[axis1] * sintheta + P[i].rVelocity[axis2] * costheta;
  }
  return;
}*/

