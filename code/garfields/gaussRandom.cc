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
#include <gsl/gsl_randist.h>
#include "fftOperators.h" 
#include "functions.h" 

void gaussRandom(int *N, fftw_real *Power, fftw_complex *gauss_random, gsl_rng *seed)
{
 for(int i=0;i<N[0];i++)
     for(int j=0;j<N[1];j++)
       for(int k=0;k<N[2];k++)
	 {			
	   double sigma = sqrt(Power[k+N[2]*(j+N[1]*i)]/ 2.0);

	   if((i == 0) && (j == 0) && (k == 0))
	     {
	       re(gauss_random[k+N[2]*(j+N[1]*i)])=0.0;
	       im(gauss_random[k+N[2]*(j+N[1]*i)])=0.0;
	     }
	   else if((i == N[0]/2) && (j == N[1]/2) && (k == N[2]/2))
	     {					
	       re(gauss_random[k+N[2]*(j+N[1]*i)])=sqrt(2.)*gsl_ran_gaussian(seed,sigma);;
	       im(gauss_random[k+N[2]*(j+N[1]*i)])=0.0;
	     }
	   else if((i == 0) && (j == N[1]/2) && (k == N[2]/2))
	     {					
	       re(gauss_random[k+N[2]*(j+N[1]*i)])=sqrt(2.)*gsl_ran_gaussian(seed,sigma);;
	       im(gauss_random[k+N[2]*(j+N[1]*i)])=0.0;
	     }
	   else if((i == N[0]/2) && (j == 0) && (k == N[2]/2))
	     {					
	       re(gauss_random[k+N[2]*(j+N[1]*i)])=sqrt(2.)*gsl_ran_gaussian(seed,sigma);;
	       im(gauss_random[k+N[2]*(j+N[1]*i)])=0.0;
	     }
	   else if((i == N[0]/2) && (j == N[1]/2) && (k == 0))
	     {					
	       re(gauss_random[k+N[2]*(j+N[1]*i)])=sqrt(2.)*gsl_ran_gaussian(seed,sigma);;
	       im(gauss_random[k+N[2]*(j+N[1]*i)])=0.0;
	     }
	   else if((i == 0) && (j == 0) && (k == N[2]/2))
	     {					
	       re(gauss_random[k+N[2]*(j+N[1]*i)])=sqrt(2.)*gsl_ran_gaussian(seed,sigma);;
	       im(gauss_random[k+N[2]*(j+N[1]*i)])=0.0;
	     }
	   else if((i == 0) && (j == N[1]/2) && (k == 0))
	     {					
	       re(gauss_random[k+N[2]*(j+N[1]*i)])=sqrt(2.)*gsl_ran_gaussian(seed,sigma);;
	       im(gauss_random[k+N[2]*(j+N[1]*i)])=0.0;
	     }
	   else if((i == N[0]/2) && (j == 0) && (k == 0))
	     {					
	       re(gauss_random[k+N[2]*(j+N[1]*i)])=sqrt(2.)*gsl_ran_gaussian(seed,sigma);;
	       im(gauss_random[k+N[2]*(j+N[1]*i)])=0.0;
	     }
	   else
	     {
	       re(gauss_random[k+N[2]*(j+N[1]*i)])=gsl_ran_gaussian(seed,sigma);
	       im(gauss_random[k+N[2]*(j+N[1]*i)])=gsl_ran_gaussian(seed,sigma);
	       
	     }
	 }

 HERMITICITY(gauss_random,N[0],N[1],N[2]);
}

