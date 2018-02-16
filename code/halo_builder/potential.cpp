/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, NYU CCPP
  Date: Feb 14, 2014

  Triaxial halo generation program

*/

//****************** potential.cpp **************************

#include "potential.h"

Potential::Potential(string inname, int NPot) //Constructor                                                                                            
{
  // This reads in the potential data and initializes interpolation objects
  // for the potential and for PPhi(R) and QPhi(R), which determine the 
  // shape factors for the equipotential surfaces.

  string line;
  vector<string> values;
  int i, j;
  double rx, phix;
  double** xa = new double*[3];
  double ** phia = new double*[3];
  double ** fa = new double*[3];
  phi = new gsl_spline*[3];
  f = new gsl_spline*[3];
  phiacc = new gsl_interp_accel*[3];
  facc = new gsl_interp_accel*[3];

  ifstream infile(inname.c_str());

  // First we initialize the potential interpolation objects
  for (i=0; i<3; i++)
  {
    xa[i]= new double[NPot];
    phia[i]= new double[NPot];
    fa[i]= new double[NPot];

    phi[i] = new gsl_spline;
    f[i] = new gsl_spline;
    phiacc[i] = new gsl_interp_accel;
    facc[i] = new gsl_interp_accel;

  // Now we read the data from the file
    getline (infile,line); // Discard the first line, which is labels
    for (j=0; j<NPot; j++)
      {
	getline (infile,line);
	boost::split(values, line, boost::is_any_of(" ,\t"));
	if (values.size() == 8)
	  {
	    xa[i][j] = atof(values[i].c_str());
	    phia[i][j] = atof(values[3].c_str());
	    fa[i][j] = atof(values[i + 4].c_str());
	  }
	else continue;
      }
    phiacc[i] = gsl_interp_accel_alloc();
    phi[i] = gsl_spline_alloc (gsl_interp_cspline, NPot);  
    gsl_spline_init (phi[i],  xa[i], phia[i], NPot);
    facc[i] = gsl_interp_accel_alloc();
    f[i] = gsl_spline_alloc (gsl_interp_cspline, NPot);  
    gsl_spline_init (f[i],  xa[i], fa[i], NPot);

  }
  infile.close();

  rmin = xa[0][0];
  rmax = xa[0][NPot - 1];

  // Now we initialize the P(R), Q(R) interpolation objects
  double* qphia = new double[NPot];
  double* pphia = new double[NPot];
  qphi = new gsl_spline;
  pphi = new gsl_spline;
  qphiacc = new gsl_interp_accel;
  pphiacc = new gsl_interp_accel;

  for (j=0; j<NPot; j++)
    {
      rx = xa[0][j];
      phix = phia[0][j];
      pphia[j] = RPhi(2, phix, rx) / rx;
      qphia[j] = RPhi(1, phix, rx) / rx;
    }
  pphiacc = gsl_interp_accel_alloc();
  pphi = gsl_spline_alloc (gsl_interp_cspline, NPot);  
  gsl_spline_init (pphi,  xa[0], pphia, NPot);
  qphiacc = gsl_interp_accel_alloc();
  qphi = gsl_spline_alloc (gsl_interp_cspline, NPot);  
  gsl_spline_init (qphi,  xa[0], qphia, NPot);

  if (i != 3 || j != NPot)
  {
	printf("Error reading Potential file! Exiting.\n");
	exit(0);
  }

  for (i=0; i<3; i++)
  {
    delete[] xa[i];
    delete[] phia[i];
    delete[] fa[i];
  }
  delete[] xa;
  delete[] phia;
  delete[] fa;
  delete[] qphia;
  delete[] pphia;
  return;
}

Potential::~Potential() //Destructor                                                                                            
{
  int i;

  for (i=0; i<3; i++)
  {
    delete phi[i];
    delete f[i];
    delete phiacc[i];
    delete facc[i];
  }
  delete[] phi;
  delete[] f;
  delete[] phiacc;
  delete[] facc;
  delete qphi;
  delete pphi;
  delete qphiacc;
  delete pphiacc;
}

double Potential::RPhi(int axis ,double Phi, double rtry)
// Function for determining R(Phi) for determining the shape parameters of the
// gas ellipsoids. Inverts Phi(R) using Newton's method
// axis=0:x, axis=1:y, axis=2:z
{
  double normf, oldnormf, divisor, newr = rtry, deltar; 
  double tol = 1E-12; // Newton's method conversion tolerance
  int counter = 0;

  normf = gsl_pow_2(gsl_spline_eval(phi[axis],rtry,phiacc[axis]) - Phi);
  while ((normf > tol) && counter < 100)
	{
	  counter = counter + 1;
	  // This is the full Newton's step delta
	  deltar = (gsl_spline_eval(phi[axis],rtry,phiacc[axis]) - Phi) / gsl_spline_eval(f[axis],rtry,facc[axis]); 
	  divisor = 0.5;
	  oldnormf = normf;
	  while(normf >= oldnormf)
	    {
		divisor = divisor * 2.0;
		if (divisor > 256.0) 
		  {
			newr = 10.0 * drand48();// new guess
			normf = gsl_pow_2(gsl_spline_eval(phi[axis],newr,phiacc[axis]) - Phi);
			break;
		  }
		newr = rtry - deltar / divisor;
  	  	normf = gsl_pow_2(gsl_spline_eval(phi[axis],newr,phiacc[axis]) - Phi);
	    }
	  rtry = newr;
	}
  if (counter > 99)  printf("Number of iterations = %d, normf = %g\n", counter,normf);
  return rtry; 
}

double Potential::PPhi(double R)
     // This gives p of the potential from the look-up table
{
  if (R > rmax)
    return gsl_spline_eval(pphi, rmax, pphiacc);
  else if (R < rmin)
    return gsl_spline_eval(pphi, rmin, pphiacc);
  else
    return gsl_spline_eval(pphi, R, pphiacc);
}

double Potential::QPhi(double R)
     // This gives q of the potential from the look-up table
{
  if (R > rmax)
    return gsl_spline_eval(qphi, rmax, qphiacc);
  else if (R < rmin)
    return gsl_spline_eval(qphi, rmin, qphiacc);
  else
    return gsl_spline_eval(qphi, R, qphiacc);
}








