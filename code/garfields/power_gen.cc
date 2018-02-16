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
#include <math.h>
#include <iostream>
#include "fftOperators.h" 
#include "functions.h"
#include "variables.h"
#include "fftw_array.h"
#include "kvector.h"
#include "power_gen.h"

using namespace std;


void power_gen(int *N, double Vol, double PowIndex, double kmin, double kmax, double *L, double *Power)
{

  /* Definition of parameters for the Power Spectrum !!! */

  
  if (PowIndex==0)
    PowIndex = 5./3.;  

  double Ndb[3];
  Ndb[0]=N[0];
  Ndb[1]=N[1];
  Ndb[2]=N[2];

  double kmaxBox=2.*cgs_pi/min(3,L)*min(3,Ndb)/2.;
  if (kmax == 0 || kmax>kmaxBox)
    kmax = kmaxBox;

  if (kmin  < 0 ) 
    kmin = kmax / (-kmin);
  // This allows kmin to be a ratio of kmax - Lage - 30-Sep13

  double kminBox=2.*cgs_pi/min(3,L);
  if (kmin  == 0 || kmin<kminBox) 
    kmin = kminBox;

  printf("kmin = %.4f, kmax = %.4f\n",kmin,kmax);

  double norm = 1.0/(4.*cgs_pi);

  fftw_array<fftw_real> kx(N[0]), ky(N[1]), kz(N[2]);

  kvector(N, L, kx, ky, kz);  
     
  for (int i=0 ; i<N[0];i++)
    for (int j=0 ; j<N[1];j++)
      for (int k=0 ; k<N[2];k++)
	{
	  double kvec=sqrt(kx[i]*kx[i]+ky[j]*ky[j]+kz[k]*kz[k]);

	  if ((i==0 && j==0 && k==0) || kvec<kmin || kvec>kmax)
	    {		     
	      Power[k+N[2]*(j+N[1]*i)] = 0.0;
	    }
	  else	
	    {	
	      Power[k+N[2]*(j+N[1]*i)]=norm*pow(kvec/kmin,-PowIndex-2);
	    }
	}
  
  Power[0]=Power[N[0]-1];
  return;
}


