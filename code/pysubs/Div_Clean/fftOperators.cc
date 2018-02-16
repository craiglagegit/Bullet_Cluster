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

/*
1)  FOURIER_DEF_1  -----> FT[x]=1/N * Sigma w x ; IFT[x]=Sigma w x
     
2)  FOURIER_DEF_2  -----> FT[x]=Sigma w x; IFT[x]=1/N * Sigma w x
*/
#define ZERO_PAD
#define FOURIER_DEF_1
#undef FOURIER_DEF_3
#define CONV_DEF_1

#include "fftOperators.h"
#include "functions.h"
#include "fftw_array.h"

#include <iostream>

void FFT ( int dim, int *sz, int factor, bool direction, fftw_complex *in, fftw_complex *out )
{
  fftw_plan fftp;

      if (direction == true)
	{
	  fftp = fftw_plan_dft(dim, sz, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	  fftw_execute(fftp);
#ifdef	  FOURIER_DEF_1 
	  complexfactor_mult(factor, 1./factor, out, out);
#endif
	}
      else
	{
	  fftp = fftw_plan_dft(dim, sz, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
	  fftw_execute(fftp);
#ifdef	  FOURIER_DEF_2 
	  complexfactor_mult(factor, 1./factor, out, out);
#endif
	}
  fftw_destroy_plan(fftp);
}

void FFTR2CC ( int dim, int *sz, int factor, fftw_real *in, fftw_complex *out )
{
  fftw_plan fftp;

  fftw_array<fftw_complex> in2(factor);

  realcomplex_eq (factor, in, in2);

	  fftp = fftw_plan_dft(dim, sz, in2, out, FFTW_FORWARD, FFTW_ESTIMATE);

	  fftw_execute(fftp);

#ifdef	  FOURIER_DEF_1 
	  complexfactor_mult(factor, 1./factor, out, out);
#endif

	  fftw_destroy_plan(fftp);
}

void FFTC2RC ( int dim, int *sz, int factor, fftw_complex *in, fftw_real *out )
{
  fftw_plan fftp;
  fftw_array<fftw_complex> out2(factor);
	
  fftp = fftw_plan_dft(dim, sz, in, out2, FFTW_BACKWARD, FFTW_ESTIMATE);

  fftw_execute(fftp);

  complexreal_eq(factor, out2, out);

#ifdef	  FOURIER_DEF_2 
	  realfactor_mult(factor, 1./factor, out, out);
#endif
	  fftw_destroy_plan(fftp);
}

void FFTR2C ( int dim, int *sz, int factor, fftw_real *in, fftw_complex *out )
{
  fftw_plan fftp;
  fftp = fftw_plan_dft_r2c(dim, sz, in, out, FFTW_ESTIMATE);

  fftw_execute(fftp);

#ifdef	  FOURIER_DEF_1 
	  complexfactor_mult(factor, 1./factor, out, out);
#endif
	  fftw_destroy_plan(fftp);
}


void FFTC2R ( int dim, int *sz, int factor, fftw_complex *in, fftw_real *out )
{
  fftw_plan fftp;
  fftp = fftw_plan_dft_c2r(dim, sz, in, out, FFTW_ESTIMATE);

  fftw_execute(fftp);

#ifdef	  FOURIER_DEF_2 
	  realfactor_mult(factor, 1./factor, out, out);
#endif
	  fftw_destroy_plan(fftp);
}

void CONVPREPR2C ( int dim, int *sz, int factor, fftw_real *in, fftw_complex *out )
{
  fftw_array<fftw_complex> dummy(factor);

  int Nx, Ny, Nz, factorpad, szpad[3];

  int fac_pad=2, fac_pad2=1;

  Nx = sz[0];
  Ny = sz[1];
  Nz = sz[2]; 

#ifdef ZERO_PAD 
  szpad[0] = fac_pad*Nx-fac_pad2;
  szpad[1] = fac_pad*Ny-fac_pad2;
  szpad[2] = fac_pad*Nz-fac_pad2;
#else
  szpad[0] = Nx;
  szpad[1] = Ny;
  szpad[2] = Nz;
#endif
 
  factorpad=factorize(dim,szpad);
 
  realcomplex_eq(factor, in, dummy);

#ifdef ZERO_PAD 
  complex_zeropad(sz, factorpad, dummy, out);
#else
  complex_eq(factorpad, dummy, out);  
#endif
  
#ifdef	  CONV_DEF_1 
  FFT (dim, szpad, factorpad, inv, out, out);
#endif
#ifdef	  CONV_DEF_2 
  FFT (dim, szpad, factorpad, fwd, out, out);
#endif
}

void CONVPREPC2C ( int dim, int *sz, int factor, fftw_complex *in, fftw_complex *out )
{
  int Nx, Ny, Nz, factorpad, szpad[3];

  int fac_pad=2, fac_pad2=1;

  Nx = sz[0];
  Ny = sz[1];
  Nz = sz[2];


#ifdef ZERO_PAD 
  szpad[0] = fac_pad*Nx-fac_pad2;
  szpad[1] = fac_pad*Ny-fac_pad2;
  szpad[2] = fac_pad*Nz-fac_pad2;
#else
  szpad[0] = Nx;
  szpad[1] = Ny;
  szpad[2] = Nz;
#endif

  factorpad=factorize(dim,szpad);
 

#ifdef ZERO_PAD 
  complex_zeropad(sz, factorpad, in, out);
#else
  complex_eq(factorpad, in, out);  
#endif

#ifdef	  CONV_DEF_1 
  FFT (dim, szpad, factorpad, inv, out, out);
#endif
#ifdef	  CONV_DEF_2 
  FFT (dim, szpad, factorpad, fwd, out, out);
#endif
}

void convolution ( int dim, int *sz, int factor, bool space_ina, bool space_inb, bool space_out, fftw_complex *in_a, fftw_complex *in_b, fftw_complex *out )
{
  int Nx, Ny, Nz, factorpad, szpad[3];

  int fac_pad=2, fac_pad2=1;

  Nx = sz[0];
  Ny = sz[1];
  Nz = sz[2];

 
#ifdef ZERO_PAD 
  szpad[0] = fac_pad*Nx-fac_pad2;
  szpad[1] = fac_pad*Ny-fac_pad2;
  szpad[2] = fac_pad*Nz-fac_pad2;
#else
  szpad[0] = Nx;
  szpad[1] = Ny;
  szpad[2] = Nz;
#endif

  factorpad=factorize(dim,szpad);

  fftw_array<fftw_complex> X(factorpad), Y(factorpad), Z(factorpad);
 
 if (space_ina == true)
   {
#ifdef ZERO_PAD 
     complex_zeropad(sz, factorpad, in_a, X);
#else
     complex_eq(factorpad, in_a, X);  
#endif

#ifdef	  CONV_DEF_1 
     FFT (dim, szpad, factorpad, inv, X, X);
#endif
#ifdef	  CONV_DEF_2 
     FFT (dim, szpad, factorpad, fwd, X, X);
#endif
   }
 else
   complex_eq (factorpad, in_a, X);
 
  if (space_inb == true)
    {
#ifdef ZERO_PAD 
      complex_zeropad(sz, factorpad, in_b, Y);
#else
      complex_eq(factorpad, in_b, Y);  
#endif

#ifdef	  CONV_DEF_1 
      FFT (dim, szpad, factorpad, inv, Y, Y);
#endif
#ifdef	  CONV_DEF_2 
      FFT (dim, szpad, factorpad, fwd, Y, Y);
#endif
    }
  else
   complex_eq (factorpad, in_b, Y);

  complex_mult (factorpad, X, Y, Z);
  
 if (space_out == true)
   {
#ifdef	  CONV_DEF_1 
     FFT (dim, szpad, factorpad, fwd, Z, Z);
#endif
#ifdef	  CONV_DEF_2 
     FFT (dim, szpad, factorpad, inv, Z, Z);
#endif
#ifdef	  FOURIER_DEF_3 
     double norm=factorpad;// Has to be checked for each definition !!!
     
     complexfactor_mult(factorpad, norm, Z, Z);
#endif
     
#ifdef ZERO_PAD 
    switch(dim)
      {
      case 1:

	for (int i=0 ; i<Nx;i++)
	  { 
	    re(out[i])=re(Z[i]);
	    im(out[i])=im(Z[i]);
	  }
	
	    break;
      case 2:

	for (int i=0 ; i<Nx;i++)
	  { 
	    for (int j=0 ; j<Ny;j++)
	      {		
		re(out[j+Ny*i])=re(Z[j+(fac_pad*Ny-fac_pad2)*i]);
		im(out[j+Ny*i])=im(Z[j+(fac_pad*Ny-fac_pad2)*i]);
	      }
	  }
	
	break;

      case 3:
	
	for (int i=0 ; i<Nx;i++)
	  { 
	    for (int j=0 ; j<Ny;j++)
	      { 
		for (int k=0 ; k<Nz;k++)
		  {
		    re(out[k+Nz*(j+Ny*i)])=re(Z[k+(fac_pad*Nz-fac_pad2)*(j+(fac_pad*Ny-fac_pad2)*i)]);
		    im(out[k+Nz*(j+Ny*i)])=im(Z[k+(fac_pad*Nz-fac_pad2)*(j+(fac_pad*Ny-fac_pad2)*i)]);
		  }
	      }
	  }

	break;
      }
#else
     complex_eq(factorpad, Z, out);     
#endif

   }
 else
   {
#ifdef	  CONV_DEF_1 
     FFT (dim, szpad, factorpad, fwd, Z, Z);
#endif
#ifdef	  CONV_DEF_2 
     FFT (dim, szpad, factorpad, inv, Z, Z);
#endif
#ifdef	  FOURIER_DEF_3 
     double norm=factorpad;// Has to be checked for each definition !!!
     
     complexfactor_mult(factorpad, norm, Z, Z);
#endif
     complex_eq (factorpad, Z, out);
   }
}

void convolutioninv ( int dim, int *sz, int factor, bool space_ina, bool space_inb, bool space_out, fftw_complex *in_a, fftw_complex *in_b, fftw_complex *out )
{
  int Nx, Ny, Nz, factorpad, szpad[3];

  int fac_pad=2, fac_pad2=1;

  Nx = sz[0];
  Ny = sz[1];
  Nz = sz[2];

 
#ifdef ZERO_PAD 
  szpad[0] = fac_pad*Nx-fac_pad2;
  szpad[1] = fac_pad*Ny-fac_pad2;
  szpad[2] = fac_pad*Nz-fac_pad2;
#else
  szpad[0] = Nx;
  szpad[1] = Ny;
  szpad[2] = Nz;
#endif

  factorpad=factorize(dim,szpad);

  fftw_array<fftw_complex> X(factorpad), Y(factorpad), Z(factorpad);
 
 if (space_ina == true)
   {
#ifdef ZERO_PAD 
     complex_zeropad(sz, factorpad, in_a, X);
#else
     complex_eq(factorpad, in_a, X);  
#endif

#ifdef	  CONV_DEF_2 
     FFT (dim, szpad, factorpad, inv, X, X);
#endif
#ifdef	  CONV_DEF_1 
     FFT (dim, szpad, factorpad, fwd, X, X);
#endif
   }
 else
   complex_eq (factorpad, in_a, X);
 
  if (space_inb == true)
    {
#ifdef ZERO_PAD 
      complex_zeropad(sz, factorpad, in_b, Y);
#else
      complex_eq(factorpad, in_b, Y);  
#endif

#ifdef	  CONV_DEF_2 
      FFT (dim, szpad, factorpad, inv, Y, Y);
#endif
#ifdef	  CONV_DEF_1 
      FFT (dim, szpad, factorpad, fwd, Y, Y);
#endif
    }
  else
   complex_eq (factorpad, in_b, Y);

  complex_mult (factorpad, X, Y, Z);
  
 if (space_out == true)
   {
#ifdef	  CONV_DEF_2 
     FFT (dim, szpad, factorpad, fwd, Z, Z);
#endif
#ifdef	  CONV_DEF_1 
     FFT (dim, szpad, factorpad, inv, Z, Z);
#endif
#ifdef	  FOURIER_DEF_3 
     double norm=factorpad;// Has to be checked for each definition !!!
     
     complexfactor_mult(factorpad, norm, Z, Z);
#endif
     
#ifdef ZERO_PAD 
    switch(dim)
      {
      case 1:

	for (int i=0 ; i<Nx;i++)
	  { 
	    re(out[i])=re(Z[i]);
	    im(out[i])=im(Z[i]);
	  }
	
	    break;
      case 2:

	for (int i=0 ; i<Nx;i++)
	  { 
	    for (int j=0 ; j<Ny;j++)
	      {		
		re(out[j+Ny*i])=re(Z[j+(fac_pad*Ny-fac_pad2)*i]);
		im(out[j+Ny*i])=im(Z[j+(fac_pad*Ny-fac_pad2)*i]);
	      }
	  }
	
	break;

      case 3:
	
	for (int i=0 ; i<Nx;i++)
	  { 
	    for (int j=0 ; j<Ny;j++)
	      { 
		for (int k=0 ; k<Nz;k++)
		  {
		    re(out[k+Nz*(j+Ny*i)])=re(Z[k+(fac_pad*Nz-fac_pad2)*(j+(fac_pad*Ny-fac_pad2)*i)]);
		    im(out[k+Nz*(j+Ny*i)])=im(Z[k+(fac_pad*Nz-fac_pad2)*(j+(fac_pad*Ny-fac_pad2)*i)]);
		  }
	      }
	  }

	break;
      }
#else
     complex_eq(factorpad, Z, out);     
#endif

   }
 else
   {
#ifdef	  CONV_DEF_2 
     FFT (dim, szpad, factorpad, fwd, Z, Z);
#endif
#ifdef	  CONV_DEF_1 
     FFT (dim, szpad, factorpad, inv, Z, Z);
#endif
#ifdef	  FOURIER_DEF_3 
     double norm=factorpad;// Has to be checked for each definition !!!
     
     complexfactor_mult(factorpad, norm, Z, Z);
#endif
     complex_eq (factorpad, Z, out);
   }
}
