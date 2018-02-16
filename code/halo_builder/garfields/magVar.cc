/*
The "Garfields" code is a Gaussian random field generator, which can also
generate divergenceless and helical fields.
Copyright (C) 2007  Francisco S. Kitaura

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Francisco Kitaura can be reached at kitaura@mpa-garching.mpg.de
*/

/*
F.S. Kitaura would like to thank Torsten Ensslin and Jens Jasche for useful
contributions for the development of this code.
*/
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>         
#include <cassert>
#include <cfloat>
#include <math.h>
#include <gsl/gsl_randist.h>

#include "fftOperators.h" 
#include "functions.h"
#include "garFields.h"
#include "fftw_array.h"
#include "variables.h"
#include "gaussRandom.h"
#include "kvector.h"
#include "divclean.h"
#include "helicity.h"

using namespace std;

void magVar(int *N, double Bvar, fftw_real *Bx, fftw_real *By, fftw_real *Bz)
{
  
  int  Ntot = N[0]*N[1]*N[2];

  if (Bvar == 0)
    Bvar = 1.0;

  double Bmod2=0.0;
 
  for (int i=0;i<N[0];i++)
    for (int j=0;j<N[1];j++)
      for (int k=0;k<N[2];k++)
	{
	  double Bxx=Bx[k+N[2]*(j+N[1]*i)];
	  double Byy=By[k+N[2]*(j+N[1]*i)];
	  double Bzz=Bz[k+N[2]*(j+N[1]*i)];
	  
	  Bmod2+=Bxx*Bxx+Byy*Byy+Bzz*Bzz;
	}
  
  Bmod2/=Ntot;

  
  double Bmod=sqrt(Bmod2);
  double Bvarsq=sqrt(Bvar);
  
  for (int i=0;i<N[0];i++)
    for (int j=0;j<N[1];j++)
      for (int k=0;k<N[2];k++)
	{
	  Bx[k+N[2]*(j+N[1]*i)]*=Bvarsq/Bmod;
	  By[k+N[2]*(j+N[1]*i)]*=Bvarsq/Bmod;
	  Bz[k+N[2]*(j+N[1]*i)]*=Bvarsq/Bmod;
	}
  
  Bmod2=0.0;
 
  for (int i=0;i<N[0];i++)
    for (int j=0;j<N[1];j++)
      for (int k=0;k<N[2];k++)
	{
	  double Bxx=Bx[k+N[2]*(j+N[1]*i)];
	  double Byy=By[k+N[2]*(j+N[1]*i)];
	  double Bzz=Bz[k+N[2]*(j+N[1]*i)];
	  
	  Bmod2+=Bxx*Bxx+Byy*Byy+Bzz*Bzz;
	}
  
  Bmod2/=Ntot;
 
  cout<<endl;
 
  cout<<" >>> Variance of the field realization: "<<Bmod2<<endl;
}
