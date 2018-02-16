/*
Craig Lage - 2-May-12

This uses an MCMC to optimize pattern alignment

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <time.h>
#define pi 3.14159265358979     // Pi

//  DATA STRUCTURES:  

struct Align  
/*  This is the x, y, and theta that optimally aligns the two data sets
    dx=d[0], dy=d[1], theta=d[2].  This also contains the allowable limits
    on dx, dy, and theta.
    This is all passed from Python.
 */
{
  double d[5], dmin[5], dmax[5];
};

struct Arr //This packages the 2D data sets.
{
  int nx, ny;
  double xmin, xmax, ymin, ymax, dx, dy, *x, *y, *data;
};

struct ArraySet 
/*This packages a set of arrays.
    data1 is the larger, unshifted simulation.  
    data2 is the smaller, measured data. 
    shifteddata1 is the simulation aligned to the data
    sigma is the sigma of the measurement. 
    numarrays is the number of arrays to be aligned.
    This is all passed from Python.
*/
{
  int numarrays;
  Arr *data1, *data2, *shifteddata1, *sigma, *mask;
};


struct Func //This is the function to be minimized
{
  ArraySet *arrayset;
  Align *align;
  double EPS;
  double f; //Holds the most recent value of the function
  int evals; // Counts the number of function evaluations
  double shift();
  double operator()(double*);
};

// FUNCTION PROTOTYPES
double kernel(double, double);
  
double interp(Arr*, double, double);

void newpoint(double*, double*, double*, double*, double*);

double prob(double, double, double);

extern "C" int optimize(ArraySet*, Align*, double*, const double);

//THIS IS THE TOP LEVEL METHOD

extern "C" int optimize(ArraySet* setptr, Align* aptr, double* fptr, const double gtol)
// 

{
  srand (time(NULL));
  Func func;

  func.evals=0;
  func.align=aptr;
  func.arrayset=setptr;

  double p[5], pmax[5], pmin[5], newp[5], bestp[5], r[5];
  double f, newf, bestf, randy, pr, T, T0 = 1.0;
  int i, kmax=10000; // 10000 iterations total

  for (i=0; i<5; i++)
    {
      p[i]=bestp[i]=func.align->d[i];
      pmin[i]=func.align->dmin[i];
      pmax[i]=func.align->dmax[i];
    }

  f = bestf = func(p);

  while (func.evals < kmax)
    {
      T = T0 * (kmax - func.evals) / kmax; // Decrease the temperature as evaluation proceeds
      for (i=0; i<5; i++)
	    {
	      r[i]=(pmax[i] - pmin[i]) * (kmax - func.evals) / kmax; // Decrease the radius of search as the evaluation proceeds
	    }
      newpoint(p, pmin, pmax, newp, r); // Find a new randomly chosen point within the radius of search
      newf = func(newp); // Evaluate the function at the new point
     
      randy = (double) rand()/RAND_MAX;
      pr = prob(f,newf,T);
      if (pr > randy)
	{
	  f = newf;
	  for (i=0; i<5; i++)
	    {
	      p[i] = newp[i];
	    }
	}
      if (f < bestf)
	{
	  bestf = f;
	  for (i=0; i<5; i++)
	    {
	      bestp[i] = p[i];
	    }
	}
    }
  for (i=0; i<5; i++)
    {
      func.align->d[i]=bestp[i];
    }
  *fptr = bestf;
  printf("Best point found,  dx=%.2f  dy=%.2f  shiftx = %.2f, shifty = %.2f, theta=%.2f  fom=%.6f\n", bestp[0],bestp[1],bestp[2],bestp[3],bestp[4],*fptr);

  return 0;
}

//  SUBROUTINES:

double Func::shift()
//This shifts and rotates one data set to align with another. 
//Data1 is the larger, unshifted simulation.  Data2 is the smaller, 
//measured data. Shiftx and Shifty shift dataset 0 relative to the others.

{
  int i,j,k,pixelsum=0;
  double fomsum=0.0, fom[arrayset->numarrays], xprime, yprime, x, y;
  for (k=0; k<arrayset->numarrays; k++)
	  {
	  fom[k] = 0.0;
	  for (i=0; i<arrayset->data2[k].nx; i++) 
	    {
	      x=arrayset->data2[k].x[i];
	      for (j=0; j<arrayset->data2[k].ny; j++)
		{
		  y=arrayset->data2[k].y[j];
		  xprime=x*cos(align->d[4])+y*sin(align->d[4])+align->d[0];
		  yprime=-x*sin(align->d[4])+y*cos(align->d[4])+align->d[1];
		  if (k == 0)
		    {
		      xprime=xprime+align->d[2];
		      yprime=yprime+align->d[3];
		    }
		  arrayset->shifteddata1[k].data[i+j*arrayset->data2[k].nx]=interp(&arrayset->data1[k],xprime,yprime);
		  fom[k] = fom[k] + pow((arrayset->shifteddata1[k].data[i+j*arrayset->data2[k].nx] - arrayset->data2[k].data[i+j*arrayset->data2[k].nx]) / arrayset->sigma[k].data[i+j*arrayset->data2[k].nx],2)  * arrayset->mask[k].data[i+j*arrayset->data2[k].nx];
		  pixelsum = pixelsum + arrayset->mask[k].data[i+j*arrayset->data2[k].nx];
		}
	    }
	//printf("fom[%d] = %.3f\n",k,fom[k]);
	fomsum = fomsum + fom[k];
	  }
  return fomsum / pixelsum;
}

double Func::operator()(double *x)
  // This is the function being minimized.
{
  int j;
  for (j=0; j<5; j++) align->d[j]=x[j];
  evals += 1;
  f=shift();
  return f;
}
 
double kernel(double deltax, double deltay)
//This is a pyramidal interpolation kernel
{
  if (deltax>=1.0 || deltay>=1.0)
    { return 0.0; }
  else
    { return (1.0-deltax)*(1.0-deltay);  }
}

double interp(Arr *data, double xprime, double yprime)
//This interpolates using the kernel
{
  int i,j,m,n;
  double deltax, deltay, d=0.0;
  i=(int)floor((xprime-data->xmin)/data->dx);
  j=(int)floor((yprime-data->ymin)/data->dy);
  for (m=i-1; m<i+2; m++) 
    {
      deltax = fabs((xprime-data->x[m])/data->dx);
      for (n=j-1; n<j+2; n++)
	{
	  deltay = fabs((yprime-data->y[n])/data->dy);
	  d = d + kernel(deltax, deltay) * data->data[m+n*data->nx];
	}
    }
  return d;
}

void newpoint(double x[], double xmin[], double xmax[], double newx[], double r[])
{
  int i;
  double randy;

  for (i=0; i<5; i++)
    {
	newx[i] = xmin[i] - 0.01;
	while (newx[i] < xmin[i] || newx[i] > xmax[i])
	  {
		randy = (double) rand() / RAND_MAX;
                //printf("randy=%.6f\n",randy);
		newx[i] = (x[i] - r[i]) + 2 * r[i] * randy;
	  } 
    }
}

double prob(double f, double newf, double T)
{
    return exp((f - newf) / T);
}

