
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
#include "fftOperators.h" 

double cross_prod(int index, double Ax, double Ay, double Az, double Bx, double By, double Bz)
{
  double out;

  switch(index)
    {
    case 1:
      {
	out=Ay*Bz-Az*By;
	break;
      }
    case 2:
      {
	out=Az*Bx-Ax*Bz;
	break;
      }
    case 3:
      {
	out=Ax*By-Ay*Bx;
	break;
      }
    }
  
  return out;
}


void helicity(double alpha, int *N, fftw_real *kx, fftw_real *ky, fftw_real *kz, fftw_complex *Bx, fftw_complex *By, fftw_complex *Bz, fftw_complex *Bxhel, fftw_complex *Byhel, fftw_complex *Bzhel)
{
  double eps=1e-30;
	     
  for (int i=0;i<N[0];i++)
    for (int j=0;j<N[1];j++)
      for (int k=0;k<N[2];k++)
	{	
	  double norm=1./sqrt(1.+alpha*alpha);
	  double kmod=sqrt(kx[i]*kx[i]+ky[j]*ky[j]+kz[k]*kz[k]);

	  double re_Bx=re(Bx[k+N[2]*(j+N[1]*i)]);
	  double im_Bx=im(Bx[k+N[2]*(j+N[1]*i)]);

	  double re_By=re(By[k+N[2]*(j+N[1]*i)]);
	  double im_By=im(By[k+N[2]*(j+N[1]*i)]);
	  
	  double re_Bz=re(Bz[k+N[2]*(j+N[1]*i)]);
	  double im_Bz=im(Bz[k+N[2]*(j+N[1]*i)]);
	  
	  
	  double re_kxB_x=cross_prod(1, kx[i], ky[j], kz[k], re_Bx, re_By, re_Bz);
	  double im_kxB_x=cross_prod(1, kx[i], ky[j], kz[k], im_Bx, im_By, im_Bz);	 

	  double re_kxB_y=cross_prod(2, kx[i], ky[j], kz[k], re_Bx, re_By, re_Bz);
	  double im_kxB_y=cross_prod(2, kx[i], ky[j], kz[k], im_Bx, im_By, im_Bz);	 

	  double re_kxB_z=cross_prod(3, kx[i], ky[j], kz[k], re_Bx, re_By, re_Bz);
	  double im_kxB_z=cross_prod(3, kx[i], ky[j], kz[k], im_Bx, im_By, im_Bz);

	  if (kmod>eps)
	    { 
	      re(Bxhel[k+N[2]*(j+N[1]*i)])=norm*(re_Bx-alpha*im_kxB_x/kmod);
	      im(Bxhel[k+N[2]*(j+N[1]*i)])=norm*(im_Bx+alpha*re_kxB_x/kmod);
	 
	      re(Byhel[k+N[2]*(j+N[1]*i)])=norm*(re_By-alpha*im_kxB_y/kmod);
	      im(Byhel[k+N[2]*(j+N[1]*i)])=norm*(im_By+alpha*re_kxB_y/kmod);
	      
	      re(Bzhel[k+N[2]*(j+N[1]*i)])=norm*(re_Bz-alpha*im_kxB_z/kmod);
	      im(Bzhel[k+N[2]*(j+N[1]*i)])=norm*(im_Bz+alpha*re_kxB_z/kmod);
	    }
	  else
	    {
	      re(Bxhel[k+N[2]*(j+N[1]*i)])=0.0;
	      im(Bxhel[k+N[2]*(j+N[1]*i)])=0.0;
	 				   
	      re(Byhel[k+N[2]*(j+N[1]*i)])=0.0;
	      im(Byhel[k+N[2]*(j+N[1]*i)])=0.0;
	      				   
	      re(Bzhel[k+N[2]*(j+N[1]*i)])=0.0;
	      im(Bzhel[k+N[2]*(j+N[1]*i)])=0.0;
	    }	  
	}
}

