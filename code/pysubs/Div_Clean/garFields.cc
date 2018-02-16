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
   **************************************************************************
   GARFIELDS: GAussian Random FIELDS generator

   Author:  Francisco Kitaura

            Max Planck Institute for Astrophysics

   Date:    02/05/2006

   Revision: 05/07/2008

   Version: 2.1

   compile: make -f garFields.make

   parameter-file: garFields_par

   generates: gaussian random + vector fields
                              + divergence free vector fields
			      + helical vector fields
   **************************************************************************
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
#include "magVar.h"

using namespace std;
/*#include "bstream.h"*/

int main ()
{
/*
    ****************************************************** 
    >>> Starting Logo <<< 
    ******************************************************
*/  
  banner_garf(1);

/*
    ****************************************************** 
    >>> Read parameter file <<< 
    ******************************************************
*/  

  read_parfil("garFields_par");

  cout<<"\n >>> Reading parameter file: "<< "garFields_par" <<"...\n";
  

/*
    ****************************************************** 
    >>> declare variables <<< 
    ******************************************************
*/  

  int N[3];
  N[0] = Nx;
  N[1] = Ny;
  N[2] = Nz;

  int  Ntot = N[0]*N[1]*N[2];

  Lx=6000;
  Ly=3000;
  Lz=3000;

  fftw_real L[3];
  L[0] = Lx;
  L[1] = Ly;
  L[2] = Lz;

  Vol = L[0]*L[1]*L[2];


  fftw_array<fftw_real> Power(Ntot);

  fftw_array<fftw_real> Bx(Ntot), By(Ntot), Bz(Ntot);
  
  fftw_array<fftw_complex> BxF(Ntot), ByF(Ntot), BzF(Ntot);

  fftw_array<fftw_real> kx(Nx), ky(Ny), kz(Nz);


/*
    ****************************************************** 
    >>> Generate the power spectrum <<< 
    ******************************************************
*/  

  cout<<"\n >>> Generating the power spectrum ...\n";

  power_gen(N, Vol, L, Power);

  
/*
    ****************************************************** 
    >>> Generate Gaussian vector field 
    given the power spectrum defined above <<< 
    ******************************************************
*/  

 cout<<"\n >>> Generating the Gaussian vector field ...\n";

 // Initialize random seed 

 gsl_rng* gseed;
 const gsl_rng_type *T;
 srand(time(NULL));        
 unsigned long randSeed;	
 randSeed = seed;
 gseed = gsl_rng_alloc(gsl_rng_mt19937);
 gsl_rng_set(gseed,randSeed);
  
 gaussRandom(N, Power, BxF, gseed);
 gaussRandom(N, Power, ByF, gseed);
 gaussRandom(N, Power, BzF, gseed);

 
/*
    ****************************************************** 
    >>> define the wave vector k <<< 
    ******************************************************
*/  

 kvector(N, L, kx, ky, kz);
 

/*
    ****************************************************** 
    >>> Clean the vector field from divergence <<< 
    ******************************************************
*/  
     
 if (clean_div == true)      
   {
     cout<<"\n >>> Cleaning the divergence ...\n";
     divclean(N, kx, ky, kz, BxF, ByF, BzF, BxF, ByF, BzF);
   }
 

 if (alpha!=0)
   {
     cout<<"\n >>> Adding the helical field ...\n";
     helicity(alpha, N, kx, ky, kz, BxF, ByF, BzF, BxF, ByF, BzF);
   }
 

/*
    ****************************************************** 
    >>> Transform to real-space <<< 
    ******************************************************
*/
 
 
 FFTC2RC(3, N, Ntot, BxF, Bx);
 FFTC2RC(3, N, Ntot, ByF, By);
 FFTC2RC(3, N, Ntot, BzF, Bz);

 magVar(N, Bx, By, Bz);
 

/*
    ****************************************************** 
    >>> Regular field  write output file <<< 
    ******************************************************
*/


 
  string outputFileName=out_file;
  cout<<"\n >>> Writing the magnetic vector field on: "<<outputFileName<<" ...\n";

  
  ofstream outStream;
  outStream.open(outputFileName.data());

  for (int i=0;i<Ntot;i++)
    {   
      outStream<<Bx[i]<<" "<<By[i]<<" "<<Bz[i]<<endl;
    }

    /*
  bofstream outStream(outputFileName.data(),file_is_natural);
  assert(outStream.is_open());
  
  outStream.put(Bx.data,Ntot);
  outStream.put(By.data,Ntot);
  outStream.put(Bz.data,Ntot);
    */

  outStream.close();

 
  cout<<"\n >>> GARFIELDS finished correctly !\n\n";
 
}


