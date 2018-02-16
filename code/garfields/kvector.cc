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
#include "kvector.h"

void kvector(int *N, double *L, double *kx, double *ky, double *kz)
{
  double  kfac=2*cgs_pi;

	   for (int i=0 ; i<N[0];i++)
	     {	     
	       if (i<=N[0]/2)
		 kx[i] = kfac/L[0]*i;
	       else
		 kx[i] = -kfac/L[0]*(N[0]-i);		
	     }	  		 

	   for (int i=0 ; i<N[1];i++)
	     {
	       if (i<=N[1]/2)
		 ky[i] = kfac/L[1]*i;		
	       else
		 ky[i] = -kfac/L[1]*(N[1]-i);		
	     }	     

	   for (int i=0 ; i<N[2];i++)
	     {
	       if (i<=N[2]/2)
		 kz[i] = kfac/L[2]*i;		
	       else
		 kz[i] = -kfac/L[2]*(N[2]-i);		
	     }
}
