//****************** lmnarray.c **************************

#include "triaxial.h"

double RhoDM(double R, struct Halo *H)
{
  // This is the DM density at radius R
  return H->RhoDMNorm /((R / H->Rs) * gsl_pow_2(1.0 + (R / H->Rs)));
}

double MassShellDM(double R, void *p)
{
   // This is the mass in a shell of radius R - AEllipse is the surface area of the ellipse - it is ~ 1.0
   struct UParams * params = (struct UParams *)p;
   return RhoDM(R,params->H) * 4.0 * M_PI * gsl_pow_2(R) * params->H->P * params->H->Q;
}

double MassDM(double R, struct Halo *H)
{
    // This is the DM mass inside radius R
    struct UParams params;
    params.H = H;
    gsl_function F;
    F.function=&MassShellDM;
    F.params=&params;
    double result, error;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
    //printf("R = %f\n",R);
    //fflush(stdout);
    gsl_integration_qag(&F, 0.0, R, 0.0, 1.0E-3, 1000, GSL_INTEG_GAUSS51, w, &result, &error);
    //printf("R = %f, result = %g, error = %g\n",R,result,error);
    //fflush(stdout);
    gsl_integration_workspace_free(w);
    return result;
}

double Rho(double R, struct Halo *H)
{
  // This is the gas density at radius R
  return H->RhoNorm * (pow(1.0 + gsl_pow_2(R / H->GasRs), - H->GasAlpha) * pow((1 + gsl_pow_2(R / H->GasRcool)), - (H->GasBeta - H->GasAlpha)) * pow((1 + gsl_pow_2(R / H->GasRcool2)), - (H->GasBeta2 - H->GasBeta))) ;
}

double MassShell(double R, void *p)
{
   // This is the mass in a shell of radius R - AEllipse is the surface area of the ellipse - it is ~ 1.0
   struct UParams * params = (struct UParams *)p;
   return Rho(R,params->H) * 4.0 * M_PI * gsl_pow_2(R) * Pphi(R,params->PQ,params->H) * Qphi(R,params->PQ,params->H);
}

double Mass(double R, struct PQPhi *PQ, struct Halo *H)
{
    // This is the gas mass inside radius R
    struct UParams params;
    params.H = H; params.PQ = PQ;
    gsl_function F;
    F.function=&MassShell;
    F.params=&params;
    double result, error;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
    gsl_integration_qag(&F, 0.0, R, 0.0, 1.0E-3, 1000, GSL_INTEG_GAUSS51, w, &result, &error);
    gsl_integration_workspace_free(w);
    return result;
}

void SetRhoNorm(struct PQPhi *PQ, struct Halo *H)
{
  H->RhoNorm = 1.0 / Mass(H->C, PQ, H);
  H->RhoDMNorm = 1.0 / MassDM(H->C, H);
  return;
}

double UIntegrand(double R, void *p)
{
   // This is the integrand for calculating the internal energy
   struct UParams * params = (struct UParams *)p;
   return Rho(R,params->H) * gsl_spline_eval(params->Pot->f[0],R,params->Pot->facc[0]);
}

double U(double R, struct Potential *Pot, struct Halo *H)
{
    // This is the internal energy at radius R
    struct UParams params;
    params.H = H; params.Pot = Pot;
    gsl_function F;
    F.function=&UIntegrand;
    F.params=&params;
    double result, error;
    double Rmax = 3.0 * H->C;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc(100000);
    gsl_integration_qag(&F, R, Rmax, 0.0, 1.0E-5, 100000, GSL_INTEG_GAUSS51, w, &result, &error);
    gsl_integration_workspace_free(w);

    return 3.0 / 2.0 * result / Rho(R,H);
}

double Radius(double m, struct PQPhi *PQ, struct Halo *H)
   // Function for determining radius, given m
   // inverts M(R) using Newton's method

{
  double normf, oldnormf, divisor, newr, deltar; 
  double tol = 1E-18; // Newton's method conversion tolerance
  double rtry = m; // initial guess
  int counter = 0;
  double NormR(double r, double m, struct PQPhi *PQ, struct Halo *H)
	{
	  return gsl_pow_2(Mass(r,PQ,H) - m);
	}

  normf = NormR(rtry,m,PQ,H);
  while ((normf > tol) && counter < 100)
	{
  	  //printf("In counter loop, counter = %d, normf = %f, rtry = %f\n", counter, normf, rtry);
	  counter = counter + 1;
	  // This is the full Newton's step delta
	  deltar = (Mass(rtry,PQ,H) - m) / (4 * pi * rtry * rtry * Rho(rtry,H));
	  //printf("Deltar = %f\n", deltar);
	  divisor = 0.5;
	  oldnormf = normf;
	  while(normf >= oldnormf)
	    {
		divisor = divisor * 2.0;
		if (divisor > 256.0) 
		  {
			//printf("Binary search count exceeded. Trying new guess. normf = %f, # of iterations = %d\n",normf,counter);
			newr = 10.0 * drand48();// new guess
			normf = NormR(newr,m,PQ,H);
			break;
		  }
		newr = rtry - deltar / divisor;
  	  	normf = NormR(newr,m,PQ,H);
		//printf("In divisor loop. divisor = %f, newr = %f, normf = %f\n",divisor,newr, normf);
	    }
	  rtry = newr;
	}
  if (counter > 99)  printf("Number of iterations = %d, normf = %g\n", counter,normf);
  return rtry; 
}

struct UMass UMLookup(struct Potential *Pot, struct PQPhi *PQ, struct Halo *H)
{
  // This builds look-up tables for R(m) and U(R).  R is on a logarithmic scale
  struct UMass umass;
  int i;
  double mmin = H->mmin, mmax = 2.0, rmin = H->rmin, rmax = 30.0 * H->C,  dr, dm;
  dr = log10(rmax / rmin) / ((double)H->NLook - 1);
  dm = log10(mmax / mmin) / ((double)H->NLook - 1);
  umass.rx = calloc(H->NLook, sizeof(double));
  umass.ux = calloc(H->NLook, sizeof(double));
  umass.r = calloc(H->NLook, sizeof(double));
  umass.u = calloc(H->NLook, sizeof(double));

  for (i=0; i<H->NLook; i++)
  {
	umass.rx[i] = mmin * pow(10.0,(double)i * dm);
	umass.ux[i] = rmin * pow(10.0,(double)i * dr);
	umass.r[i] = Radius(umass.rx[i], PQ, H);
	umass.u[i] = U(umass.ux[i], Pot, H);
	//if (i % 5 ==0) printf("i = %d, rx = %f, r = %f, ux = %f, u = %f\n",i,umass.rx[i],umass.r[i],umass.ux[i],umass.u[i]);
  }
  umass.racc = gsl_interp_accel_alloc();
  umass.uacc = gsl_interp_accel_alloc();
  umass.rinterp = gsl_spline_alloc (gsl_interp_cspline, H->NLook);  
  gsl_spline_init (umass.rinterp,  umass.rx, umass.r, H->NLook);
  umass.uinterp = gsl_spline_alloc (gsl_interp_cspline, H->NLook);  
  gsl_spline_init (umass.uinterp,  umass.ux, umass.u, H->NLook);
  return umass;
}

struct PQPhi PQLookup(struct Potential *Pot, struct Halo *H)
{
  // This builds look-up tables for Pphi(R), and Qphi(R).  R is on a logarithmic scale
  struct PQPhi pqphi;
  int i;
  double rmin = H->rmin, rmax = 30.0 * H->C,  dr, phix;
  dr = log10(rmax / rmin) / ((double)H->NLook - 1);
  pqphi.rx = calloc(H->NLook, sizeof(double));
  pqphi.qphi = calloc(H->NLook, sizeof(double));
  pqphi.pphi = calloc(H->NLook, sizeof(double));

  for (i=0; i<H->NLook; i++)
  {
	pqphi.rx[i] = rmin * pow(10.0,(double)i * dr);
	phix = gsl_spline_eval(Pot->phi[0],pqphi.rx[i],Pot->phiacc[0]);
	pqphi.pphi[i] = RPhi(2, phix, pqphi.rx[i], Pot, H) / pqphi.rx[i];
	pqphi.qphi[i] = RPhi(1, phix, pqphi.rx[i], Pot, H) / pqphi.rx[i];
  }
  pqphi.pphiacc = gsl_interp_accel_alloc();
  pqphi.qphiacc = gsl_interp_accel_alloc();
  pqphi.pphiinterp = gsl_spline_alloc (gsl_interp_cspline, H->NLook);  
  gsl_spline_init (pqphi.pphiinterp,  pqphi.rx, pqphi.pphi, H->NLook);
  pqphi.qphiinterp = gsl_spline_alloc (gsl_interp_cspline, H->NLook);  
  gsl_spline_init (pqphi.qphiinterp,  pqphi.rx, pqphi.qphi, H->NLook);
  return pqphi;
}

double RPhi(int axis ,double phi, double rtry, struct Potential *Pot, struct Halo *H)
// Function for determining R(Phi) for determining the shape parameters of the
// gas ellipsoids. Inverts Phi(R) using Newton's method
// axis=0:x, axis=1:y, axis=2:z
{
  double normf, oldnormf, divisor, newr, deltar; 
  double tol = 1E-12; // Newton's method conversion tolerance
  int counter = 0;

  normf = gsl_pow_2(gsl_spline_eval(Pot->phi[axis],rtry,Pot->phiacc[axis]) - phi);
  while ((normf > tol) && counter < 100)
	{
	  counter = counter + 1;
	  // This is the full Newton's step delta
	  deltar = (gsl_spline_eval(Pot->phi[axis],rtry,Pot->phiacc[axis]) - phi) / gsl_spline_eval(Pot->f[axis],rtry,Pot->facc[axis]); 
	  divisor = 0.5;
	  oldnormf = normf;
	  while(normf >= oldnormf)
	    {
		divisor = divisor * 2.0;
		if (divisor > 256.0) 
		  {
			newr = 10.0 * drand48();// new guess
			normf = gsl_pow_2(gsl_spline_eval(Pot->phi[axis],newr,Pot->phiacc[axis]) - phi);
			break;
		  }
		newr = rtry - deltar / divisor;
  	  	normf = gsl_pow_2(gsl_spline_eval(Pot->phi[axis],newr,Pot->phiacc[axis]) - phi);
	    }
	  rtry = newr;
	}
  if (counter > 99)  printf("Number of iterations = %d, normf = %g\n", counter,normf);
  return rtry; 
}






