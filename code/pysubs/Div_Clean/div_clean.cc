/*
Craig Lage - 22-Jan-14

This cleans the divergence from a B-field array
It basically packages my 3D arrays so I can use code from the GarFields software package.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <time.h>
#define pi 3.14159265358979     // Pi

#include "fftOperators.h" 
#include "functions.h"
#include "fftw_array.h"
#include "kvector.h"
#include "divclean.h"

//  DATA STRUCTURES:  

struct Arr //This packages the 3D data sets.
{
  int nx, ny, nz;
  double xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, *x, *y, *z, *data;
};

// FUNCTION PROTOTYPES

extern "C" int DivClean(Arr* Bxptr, Arr* Byptr, Arr* Bzptr);

//THIS IS THE TOP LEVEL METHOD

extern "C" int DivClean(Arr* Bxptr, Arr* Byptr, Arr* Bzptr)
// 
{
  int i, N[3];

  N[0] = Bxptr->nx;
  N[1] = Bxptr->ny;
  N[2] = Bxptr->nz;

  int  Ntot = N[0]*N[1]*N[2];

  fftw_real L[3];
  L[0] = Bxptr->xmax - Bxptr->xmin;
  L[1] = Bxptr->ymax - Bxptr->ymin;
  L[2] = Bxptr->zmax - Bxptr->zmin;

/*
    ****************************************************** 
    >>> Put data in FFT arrays
    ******************************************************
*/

  fftw_array<fftw_real> Bx(Ntot), By(Ntot), Bz(Ntot);
  for (i=0; i<Ntot; i++)
    {
      Bx[i] = Bxptr->data[i];
      By[i] = Byptr->data[i];
      Bz[i] = Bzptr->data[i];
    }
  
/*
    ****************************************************** 
    >>> define the wave vector k <<< 
    ******************************************************
*/  

  fftw_array<fftw_real> kx(N[0]), ky(N[1]), kz(N[2]);
  kvector(N, L, kx, ky, kz);

/*
    ****************************************************** 
    >>> Transform to Fourier-space <<< 
    ******************************************************
*/
  fftw_array<fftw_complex> BxF(Ntot), ByF(Ntot), BzF(Ntot);
  FFTR2CC(3, N, Ntot, Bx, BxF);
  FFTR2CC(3, N, Ntot, By, ByF);
  FFTR2CC(3, N, Ntot, Bz, BzF);

/*
    ****************************************************** 
    >>> Clean the vector field from divergence <<< 
    ******************************************************
*/  
     
  divclean(N, kx, ky, kz, BxF, ByF, BzF, BxF, ByF, BzF);

/*
    ****************************************************** 
    >>> Transform to real-space <<< 
    ******************************************************
*/
 
  FFTC2RC(3, N, Ntot, BxF, Bx);
  FFTC2RC(3, N, Ntot, ByF, By);
  FFTC2RC(3, N, Ntot, BzF, Bz);

/*
    ****************************************************** 
    >>> Put data back in incoming arrays
    ******************************************************
*/

  for (i=0; i<Ntot; i++)
    {
      Bxptr->data[i] = Bx[i];
      Byptr->data[i] = By[i];
      Bzptr->data[i] = Bz[i];
    }
 return 0;
}


