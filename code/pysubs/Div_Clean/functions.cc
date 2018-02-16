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

#define FOURIER_DEF_1

#define max_element(x,y) (x>y?x:y)      // max macro definition
#define min_element(x,y) (x>y?y:x)      // min macro definition

#include <math.h>
#include "fftOperators.h"
#include "functions.h"
#include "fftw_array.h"

#include <iostream>
using namespace std;
 

double min ( int factor, double *in )
{
  double *firstn = in;
  double *lastn  = in + factor;
  double minn;

       //get min
      double *min = min_element(firstn,lastn);
      minn = *min;
 
  return minn;
}
    
double max ( int factor, double *in )
{
  double *firstn = in;
  double *lastn  = in + factor;
  double max;

       //get max
      double *maxn = max_element(firstn,lastn);
      
      max = *maxn;

  return max;
}

int factorize ( int dim, int *sz )
{
  int factor = 1;
  for (int i = 0; i < dim; i++)
    {
      factor *=  sz[i];
    }
  return factor;
}


fftw_real factorial(int term)
{
  double out=1.0;

  for (int i=1;i<=term;i++)    
    out*=i*1.0;

  return out;
}

void complex_mult ( int factor, fftw_complex *in_a, fftw_complex *in_b, fftw_complex *out )
{
  for (int i = 0; i < factor; i++) 
    {
      re(out[i]) = re(in_a[i])*re(in_b[i]) - im(in_a[i])*im(in_b[i]);
      im(out[i]) = re(in_a[i])*im(in_b[i]) + im(in_a[i])*re(in_b[i]);
    }
}

void complex_cfactor_mult ( int factor, fftw_complex in_a, fftw_complex *in_b, fftw_complex *out )
{
  for (int i = 0; i < factor; i++) 
    {
      re(out[i]) = re(in_a)*re(in_b[i]) - im(in_a)*im(in_b[i]);
      im(out[i]) = re(in_a)*im(in_b[i]) + im(in_a)*re(in_b[i]);
    }
}

void complexfactor_mult ( int factor, fftw_real in_a, fftw_complex *in_b, fftw_complex *out )
{
  for (int i = 0; i < factor; i++) 
    {
      re(out[i]) = in_a*re(in_b[i]);
      im(out[i]) = in_a*im(in_b[i]);
    }
}

void complex_mult1 (  fftw_complex in_a, fftw_complex in_b, fftw_complex out )
{
      re(out) = re(in_a)*re(in_b) - im(in_a)*im(in_b);
      im(out) = re(in_a)*im(in_b) + im(in_a)*re(in_b);
}

void complex_vec_mult ( int factor, fftw_complex *in_a, fftw_complex *in_b, fftw_complex out )
{
  fftw_complex dummy;

  /*  dummy = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) );
   */
      re(out) = 0.0;
      im(out) = 0.0;

   for (int i = 0; i < factor; i++) 
    {
      re(dummy) = re(in_a[i])*re(in_b[i]) - im(in_a[i])*im(in_b[i]);
      im(dummy) = re(in_a[i])*im(in_b[i]) + im(in_a[i])*re(in_b[i]);

      re(out) = re(out) + re(dummy);
      im(out) = im(out) + im(dummy);
    }

}

void complexreal_mult ( int factor, fftw_complex *in_a, fftw_complex *in_b, fftw_complex *out )
{
  for (int i = 0; i < factor; i++) 
    {
      re(out[i]) = re(in_a[i])*re(in_b[i]) - im(in_a[i])*im(in_b[i]);
      im(out[i]) = 0.0;
    }
}

void realcomplex_mult ( int factor, fftw_real *in_a, fftw_complex *in_b, fftw_complex *out )
{
  for (int i = 0; i < factor; i++) 
    {
      re(out[i]) = in_a[i]*re(in_b[i]);
      im(out[i]) = in_a[i]*im(in_b[i]);
    }
}

void realfactor_mult ( int factor, fftw_real in_a, fftw_real *in_b, fftw_real *out )
{
  for (int i = 0; i < factor; i++) 
    {
      out[i] = in_a*in_b[i];
    }
}

void real_mult ( int factor, fftw_real *in_a, fftw_real *in_b, fftw_real *out )
{
  for (int i = 0; i < factor; i++) 
    {
      out[i] = in_a[i]*in_b[i];
    }
}


fftw_real real_average ( int factor, fftw_real *in )
{
  fftw_real out=0.0;
  for (int i = 0; i < factor; i++) 
    {
      out += in[i];
    }
  out /= factor;

  return out;
}


void real_delta ( int factor, fftw_real *in, fftw_real *out )
{
  fftw_real med;
  med = real_average ( factor, in);

  for (int i = 0; i < factor; i++) 
    {
      out[i] = (in[i]-med)/med;
    }
}

void real_subtract ( int factor, fftw_real *in_a, fftw_real *in_b, fftw_real *out )
{
  for (int i = 0; i < factor; i++) 
    {
      out[i] = in_a[i]-in_b[i];
    }
}

void real_eq ( int factor, fftw_real *in, fftw_real *out )
{
  for (int i = 0; i < factor; i++) 
    {
      out[i] = in[i];
    }
}

void real_identity ( int factor, fftw_real *out )
{
  for (int i = 0; i < factor; i++) 
    {
      out[i] = 1.0;
    }
}

void complex_identity ( int factor, fftw_complex in, fftw_complex *out )
{
  for (int i = 0; i < factor; i++) 
    {
      re(out[i]) = re(in);
      im(out[i]) = im(in);
    }
}


void realcomplex_eq ( int factor, fftw_real *in, fftw_complex *out )
{
  for (int i = 0; i < factor; i++) 
    {
      re(out[i]) = in[i];
      im(out[i]) = 0.0;
    }
}

void realcomplex1_eq ( fftw_real in, fftw_complex out )
{
      re(out) = in;
      im(out) = 0.0;
 }

void complexreal_eq ( int factor, fftw_complex *in, fftw_real *out  )
{
  for (int i = 0; i < factor; i++) 
    {
      out[i] = re(in[i]);
    }
}

void complexreal1_eq (fftw_complex in, fftw_real out  )
{
      out = re(in);
}

void complex_eq ( int factor, fftw_complex *in, fftw_complex *out )
{
  for (int i = 0; i < factor; i++) 
    {
      re(out[i]) = re(in[i]);
      im(out[i]) = im(in[i]);
    }
}

void complex_zero ( int factor, fftw_complex *out )
{
  for (int i = 0; i < factor; i++) 
    {
      re(out[i]) = 0.0;
      im(out[i]) = 0.0;
    }
}

void real_zero ( int factor, fftw_real *out )
{
  for (int i = 0; i < factor; i++) 
    {
      out[i] = 0.0;
    }
}

fftw_real variance ( int factor, fftw_real * in, fftw_real mu)
{
  double out=0.0;
  for (int i = 0; i < factor; i++) 
    {
      out += (in[i]-mu)*(in[i]-mu);
    }
  out /= factor;

  return out;
}

void complex_zeropad ( int *sz, int factor, fftw_complex *in, fftw_complex *out )
{
    int Nx, Ny, Nz;

    Nx = sz[0];
    Ny = sz[1];
    Nz = sz[2];
    
    int fac_pad=2, fac_pad2=1;    

    complex_zero (factor, out);

     for (int i=0 ; i<Nx;i++)
       { 
	 for (int j=0 ; j<Ny;j++)
	   { 
	     for (int k=0 ; k<Nz;k++)
	       {
		 re(out[k+(fac_pad*Nz-fac_pad2)*(j+(fac_pad*Ny-fac_pad2)*i)])=re(in[k+Nz*(j+Ny*i)]);
		 im(out[k+(fac_pad*Nz-fac_pad2)*(j+(fac_pad*Ny-fac_pad2)*i)])=im(in[k+Nz*(j+Ny*i)]);
	       }
	   }
       }
}

void real_zeropad ( int *sz, int factor, fftw_real *in, fftw_real *out )
{
    int Nx, Ny, Nz;

    Nx = sz[0];
    Ny = sz[1];
    Nz = sz[2];
    
    int fac_pad=2, fac_pad2=1;    

    real_zero (factor, out);

     for (int i=0 ; i<Nx;i++)
       { 
	 for (int j=0 ; j<Ny;j++)
	   { 
	     for (int k=0 ; k<Nz;k++)
	       {
		 out[k+(fac_pad*Nz-fac_pad2)*(j+(fac_pad*Ny-fac_pad2)*i)]=in[k+Nz*(j+Ny*i)];
	       }
	   }
       }
}

void complex_eq1 ( fftw_complex in, fftw_complex out )
{
      re(out) = re(in);
      im(out) = im(in);
}

void complex_sum ( int factor, fftw_complex *in_a, fftw_complex *in_b, fftw_complex *out )
{
  for (int i = 0; i < factor; i++) 
    {
      re(out[i]) = re(in_a[i]) + re(in_b[i]);
      im(out[i]) = im(in_a[i]) + im(in_b[i]);
    }
}

void complex_subtract ( int factor, fftw_complex *in_a, fftw_complex *in_b, fftw_complex *out )
{
  for (int i = 0; i < factor; i++) 
    {
      re(out[i]) = re(in_a[i]) - re(in_b[i]);
      im(out[i]) = im(in_a[i]) - im(in_b[i]);
    }
}

void complexreal_sum ( int factor, fftw_complex *in_a, fftw_complex *in_b, fftw_complex *out )
{
  for (int i = 0; i < factor; i++) 
    {
      re(out[i]) = re(in_a[i]) + re(in_b[i]);
      im(out[i]) = 0.0;
    }
}


void real_sum ( int factor, fftw_real *in_a, fftw_real *in_b, fftw_real *out )
{
  for (int i = 0; i < factor; i++) 
    {
      out[i] = in_a[i] + in_b[i];
    }
}

fftw_real action (int factor, fftw_real *in_a,  fftw_real *in_b)
{
  fftw_real rdummy, out;

  out=0.0;
  for (int i = 0; i < factor; i++) 
    {
      rdummy = fabs(in_a[i]-in_b[i]);
      //if (fabs(rdummy) > 1.0e-7) !!! why???
	out +=  rdummy*rdummy;      
    }
  out=fabs(out)/(factor*1.0);  
  return out;
}

void power_spectrum ( int factor, fftw_complex *in, fftw_real *out )
{
  /*fftw_complex *cin;
    cin = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * factor);
    conjugate(factor, in, cin);
    complex_mult( factor, in, cin, out);
    fftw_free (cin);*/
    
  for (int i = 0; i < factor; i++)
    {
      out[i] = re(in[i])*re(in[i]) + im(in[i])*im(in[i]);

#ifdef	FOURIER_DEF_2 
      out[i]/=(factor*factor);
#endif
    }
}

void power_fourier ( int dim, int *sz, int factor, fftw_real *in_s, fftw_real *out )
{
  fftw_array<fftw_complex> in(factor);
 
  FFTR2C (dim, sz, factor, in_s, in);

  power_spectrum (factor, in, out);
}

void complex_inverse ( int factor, fftw_complex *in, fftw_complex *out )
{
  for (int i = 0; i < factor; i++)
    {
 	  re(out[i]) = re(in[i])/(re(in[i])*re(in[i])+im(in[i])*im(in[i]));
	  im(out[i]) = -im(in[i])/(re(in[i])*re(in[i])+im(in[i])*im(in[i]));

    }
}

void complex_inverse1 ( fftw_complex in, fftw_complex out )
{
 	  re(out) = re(in)/(re(in)*re(in)+im(in)*im(in));
	  im(out) = -im(in)/(re(in)*re(in)+im(in)*im(in));
}

void inverse ( int factor,  fftw_real *in, fftw_real *out )
{
  double prec, prec2, prec3, prec4;
  
  /* delicate precision parameters !!! */
  prec=1.0e-30;
  prec2=1.0e30;
  prec3=prec2;
  prec4=prec;

  for (int i = 0; i < factor; i++)
    {
      if (fabs(in[i]) > prec && fabs(in[i]) < prec3)
	  out[i] = 1.0/in[i];
      if (fabs(in[i]) < prec)
	  out[i] = prec2;
      if (fabs(in[i]) > prec3)
	  out[i] = prec4;
    }
}

void rootdiag ( int factor,  fftw_real *in, fftw_real *out )
{
  for (int i = 0; i < factor; i++)
    {
      if (in[i] > 0.0) 
	out[i] = sqrt(in[i]);
      else
	cout<<"imaginary number in rootdiag"<<endl;
    }
}

void conjugate ( int factor, fftw_complex *in, fftw_complex *out )
{
  for (int i=0; i < factor; i++)
    {
      re(out[i]) = re(in[i]) ;
      im(out[i]) = -im(in[i]);
    }  
}

void conjugate1 (fftw_complex in, fftw_complex out )
{
      re(out) = re(in) ;
      im(out) = -im(in);
}

void hermRedun ( int dim, int *sz, int factor, fftw_complex *in, fftw_complex *out)
{
 
  int Nx, Ny, Nz;
  
  Nx = sz[0];
  Ny = sz[1];
  Nz = sz[2];


  switch(dim)
    {
    case 1:
      {	  
	for (int i=0 ; i<Nx/2+1;i++)
	  {
	    re(out[i]) = re(in[i]); 
	    im(out[i]) = im(in[i]); 
	  }
	
      } 
      break;
    case 2:
      {	   
	for (int i=0 ; i<Nx;i++)
	  {
	    for (int j=0 ; j<Ny/2+1;j++)
	      {
		re(out[(Ny/2+1)*i+j]) = re(in[Ny*i+j]); 
		im(out[(Ny/2+1)*i+j]) = im(in[Ny*i+j]); 
	      }
	  }
      }
      break;
    case 3:
      {	
	for (int i=0 ; i<Nx;i++)
	  {
	    for (int j=0 ; j<Ny;j++)
	      {
		for (int k=0 ; k<Nz/2+1;k++)
		  {
		    re(out[k+(Nz/2+1)*(Ny*i+j)]) = re(in[k+Nz*(Ny*i+j)]); 
		     im(out[k+(Nz/2+1)*(Ny*i+j)]) = im(in[k+Nz*(Ny*i+j)]); 
		  }
	      }
	  }
      }
      break;
    }
}

void HERMITICITY(fftw_complex *A_C,int N1, int N2, int N3)
{
// (0,0,0) must be real
{
	int i=0;
	int j=0;
	int k=0;

	im(A_C[k+N3*(j+N2*i)]) = 0.0;
}

//(N1/2,N2/2,N3/2) Nyquist frequency must be real
{
	int i=N1/2;
	int j=N2/2;
	int k=N3/2;

	im(A_C[k+N3*(j+N2*i)]) = 0.0;
}

// frequency (N1/2,N2/2,0)  must be real
{
	int i=N1/2;
	int j=N2/2;
	int k=0;

	im(A_C[k+N3*(j+N2*i)]) = 0.0;
	
}

// frequency (N1/2,0,N3/2) must be real
{
	int i=N1/2;
	int j=0;
	int k=N3/2;
	
	im(A_C[k+N3*(j+N2*i)]) = 0.0;
	
}

// frequency (0,N2/2,N3/2) must be real
{
	int i=0;
	int j=N2/2;
	int k=N3/2;
	
	im(A_C[k+N3*(j+N2*i)]) = 0.0;
}

// frequency (N1/2,0,0) must be real
{
	int i=N1/2;
	int j=0;
	int k=0;
	
	im(A_C[k+N3*(j+N2*i)]) = 0.0;
}

// frequency (0,N2/2,0) must be real
{
	int j=N2/2;
	int i=0;
	int k=0;
	
	im(A_C[k+N3*(j+N2*i)]) = 0.0;
}

// frequency (0,N2/2,0) must be real
{
	int i=0;
	int j=0;
	int k=N3/2;

	im(A_C[k+N3*(j+N2*i)]) = 0.0;
}

/* hermiticity: complete volume for k < N3/2 */
for(int i=1;i<N1;i++)
{
	for(int j=1;j<N2;j++)
	{
		for(int k=1;k<N3/2;k++)
		{
			int ii = N1 - i;
			int jj = N2 - j;
			int kk = N3 - k;

			re(A_C[kk+N3*(jj+N2*ii)]) = re(A_C[k+N3*(j+N2*i)]);
	   		im(A_C[kk+N3*(jj+N2*ii)]) = -im(A_C[k+N3*(j+N2*i)]);
			
		}
	}
}

/* hermiticity: plane for k = N3/2, j<N2/2 */

for(int i=1;i<N1;i++)
{
	for(int j=1;j<N2/2;j++)
	{
	int k=N3/2;
	int ii = N1 - i;
	int jj = N2 - j;
	int kk = N3 - k;

	re(A_C[kk+N3*(jj+N2*ii)]) = re(A_C[k+N3*(j+N2*i)]);
	im(A_C[kk+N3*(jj+N2*ii)]) = -im(A_C[k+N3*(j+N2*i)]);
	}
}

/* hermiticity: plane for k = 0, j<N2/2 */
for(int i=1;i<N1;i++)
{
	for(int j=1;j<N2/2;j++)
	{
	int k=0;
	int ii = N1 - i;
	int jj = N2 - j;
	
	re(A_C[k+N3*(jj+N2*ii)]) = re(A_C[k+N3*(j+N2*i)]);
	im(A_C[k+N3*(jj+N2*ii)]) = -im(A_C[k+N3*(j+N2*i)]);
	}
}

/* hermiticity: plane for j = 0, k<N3/2 */
for(int i=1;i<N1;i++)
{
	for(int k=1;k<N3/2;k++)
	{
	int j=0;
	int ii = N1 - i;
	int kk = N3 - k;
	
	re(A_C[kk+N3*(j+N2*ii)]) = re(A_C[k+N3*(j+N2*i)]);
	im(A_C[kk+N3*(j+N2*ii)]) = -im(A_C[k+N3*(j+N2*i)]);
	}
}

/* hermiticity: plane for i = 0, k<N3/2 */
for(int j=1;j<N2;j++)
{
	for(int k=1;k<N3/2;k++)
	{
	int i=0;
	int jj = N2 - j;
	int kk = N3 - k;
		
	re(A_C[kk+N3*(jj+N2*i)]) = re(A_C[k+N3*(j+N2*i)]);
	im(A_C[kk+N3*(jj+N2*i)]) = -im(A_C[k+N3*(j+N2*i)]);
	}
}

/* hermiticity: line for k = N3/2, j=N2/2 */
for(int i=1;i<N1/2;i++)
{
	int j=N2/2;
	int k=N3/2;
	int ii = N1 - i;
	int jj = N2 - j;
	int kk = N3 - k;

	re(A_C[kk+N3*(jj+N2*ii)]) = re(A_C[k+N3*(j+N2*i)]);
	im(A_C[kk+N3*(jj+N2*ii)]) = -im(A_C[k+N3*(j+N2*i)]);
	
}

/* hermiticity: line for k = 0, j=N2/2 */
for(int i=1;i<N1/2;i++)
{
	int j=N2/2;
	int k=0;
	int ii = N1 - i;
	
	re(A_C[k+N3*(j+N2*ii)]) = re(A_C[k+N3*(j+N2*i)]);
	im(A_C[k+N3*(j+N2*ii)]) = -im(A_C[k+N3*(j+N2*i)]);
	
}

/* hermiticity: line for k = N3/2, j=0 */
for(int i=1;i<N1/2;i++)
{
	int k=N3/2;
	int j=0;
	int ii = N1 - i;
	int kk =N3/2;

	re(A_C[kk+N3*(j+N2*ii)]) = re(A_C[k+N3*(j+N2*i)]);
	im(A_C[kk+N3*(j+N2*ii)]) = -im(A_C[k+N3*(j+N2*i)]);
}


/* hermiticity: line for k = N3/2, i=0 */
for(int j=1;j<N2/2;j++)
{
	int k=N3/2;
	int i=0;
	int jj = N2 - j;
	int kk =N3/2;

	re(A_C[kk+N3*(jj+N2*i)]) = re(A_C[k+N3*(j+N2*i)]);
	im(A_C[kk+N3*(jj+N2*i)]) = -im(A_C[k+N3*(j+N2*i)]);
}


/* hermiticity: line for k = 0, j=0 */
for(int i=1;i<N1/2;i++)
{
	int k=0;
	int j=0;
	int ii = N1 - i;
	
	re(A_C[k+N3*(j+N2*ii)]) = re(A_C[k+N3*(j+N2*i)]);
	im(A_C[k+N3*(j+N2*ii)]) = -im(A_C[k+N3*(j+N2*i)]);
}



/* hermiticity: line for k = 0, i=0 */
for(int j=1;j<N2/2;j++)
{
	int k=0;
	int i=0;
	int jj = N2 - j;
	
	re(A_C[k+N3*(jj+N2*i)]) = re(A_C[k+N3*(j+N2*i)]);
	im(A_C[k+N3*(jj+N2*i)]) = -im(A_C[k+N3*(j+N2*i)]);
}



/* hermiticity: line for i = 0, j=0 */
for(int k=1;k<N3/2;k++)
{
	int i=0;
	int j=0;
	int kk = N3 - k;
	
	re(A_C[kk+N3*(j+N2*i)]) = re(A_C[k+N3*(j+N2*i)]);
	im(A_C[kk+N3*(j+N2*i)]) = -im(A_C[k+N3*(j+N2*i)]);
}
}
