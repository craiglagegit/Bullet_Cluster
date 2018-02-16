
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

double gramschmidt(int index, double Ax, double Ay, double Az, double Bx, double By, double Bz)
{
  double out, Bl = 0.0, Al = 0.0, Amod2;
  double eps=1e-30;

  switch(index)
    {
    case 1:
      {
	Bl=Bx;
	Al=Ax;
	break;
      }
    case 2:
      {
	Bl=By;
	Al=Ay;
	break;
      }
    case 3:
      {
	Bl=Bz;
	Al=Az;
	break;
      }
    }

  Amod2=Ax*Ax+Ay*Ay+Az*Az;

  if (Amod2>eps)
    out = Bl-((Bx*Ax+By*Ay+Bz*Az)/Amod2)*Al;
  else
    out = 0.0;
  
  return out;
}


void divclean(int *N, fftw_real *kx, fftw_real *ky, fftw_real *kz, fftw_complex *Bx, fftw_complex *By, fftw_complex *Bz, fftw_complex *BxDivfree, fftw_complex *ByDivfree, fftw_complex *BzDivfree)
{
 	     
  double norm=sqrt(3./2.);

  for (int i=0;i<N[0];i++)
    for (int j=0;j<N[1];j++)
      for (int k=0;k<N[2];k++)
	{
	  double re_Bx=re(Bx[k+N[2]*(j+N[1]*i)]);
	  double im_Bx=im(Bx[k+N[2]*(j+N[1]*i)]);

	  double re_By=re(By[k+N[2]*(j+N[1]*i)]);
	  double im_By=im(By[k+N[2]*(j+N[1]*i)]);
	  
	  double re_Bz=re(Bz[k+N[2]*(j+N[1]*i)]);
	  double im_Bz=im(Bz[k+N[2]*(j+N[1]*i)]);
	  
	  re(BxDivfree[k+N[2]*(j+N[1]*i)])=norm*gramschmidt(1, kx[i], ky[j], kz[k], re_Bx, re_By, re_Bz);	
	  im(BxDivfree[k+N[2]*(j+N[1]*i)])=norm*gramschmidt(1, kx[i], ky[j], kz[k], im_Bx, im_By, im_Bz);

	  re(ByDivfree[k+N[2]*(j+N[1]*i)])=norm*gramschmidt(2, kx[i], ky[j], kz[k], re_Bx, re_By, re_Bz);
	  im(ByDivfree[k+N[2]*(j+N[1]*i)])=norm*gramschmidt(2, kx[i], ky[j], kz[k], im_Bx, im_By, im_Bz);

	  re(BzDivfree[k+N[2]*(j+N[1]*i)])=norm*gramschmidt(3, kx[i], ky[j], kz[k], re_Bx, re_By, re_Bz);
	  im(BzDivfree[k+N[2]*(j+N[1]*i)])=norm*gramschmidt(3, kx[i], ky[j], kz[k], im_Bx, im_By, im_Bz);
	}
}

