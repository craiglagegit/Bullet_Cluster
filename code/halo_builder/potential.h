/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, NYU CCPP
  Date: Feb 14, 2014

  Triaxial halo generation program

*/

//****************** potential.h **************************

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
using namespace std;

class Potential
{
 public:
  double rmin, rmax;  // Min and max of the interpolation routines
  gsl_spline** phi; // For the potential Phi(R)
  gsl_spline** f;   // For the field f(R)
  gsl_spline* pphi;   // For the shape parameter PPhi(R)
  gsl_spline* qphi;   // For the shape parameter QPhi(R)
  gsl_interp_accel** phiacc;
  gsl_interp_accel** facc;
  gsl_interp_accel* pphiacc;
  gsl_interp_accel* qphiacc;

  Potential() {};
  Potential(string, int);
  ~Potential();
  double RPhi(int, double, double); // Inverts Phi(R) to give R(Phi)
  double PPhi(double);  // Shape parameter PPhi(R)
  double QPhi(double);  // Shape parameter QPhi(R)
};
