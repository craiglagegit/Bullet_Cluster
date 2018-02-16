/*
Craig Lage - 2-May-12

This uses a simulated annealing method to determine a function optimum

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
  double d[3], dmin[3], dmax[3];
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
  void df(double*, double*);
};

// FUNCTION PROTOTYPES

void dfpmin(double*, double*, double*, const double, int &, double &, Func &);

void lnsrch(double*, const double, double*, double*, double*, double &, double*, double*, int &, Func &);

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

  printf("rand1 = %d, rand2 = %d, RAND_MAX = %d\n",rand(), rand(), RAND_MAX);

  func.evals=0;
  func.align=aptr;
  func.arrayset=setptr;

  double p[3], pmax[3], pmin[3], newp[3], bestp[3], r[3];
  double f, newf, bestf, randy, pr, T, T0 = 1.0;
  int i, kmax=5000;

  for (i=0; i<3; i++)
    {
      p[i]=bestp[i]=func.align->d[i];
      pmin[i]=func.align->dmin[i];
      pmax[i]=func.align->dmax[i];
    }

  f = bestf = func(p);

  printf("Getting started\n");

  while (func.evals < kmax)
    {
      printf("In while loop, func.evals = %d\n",func.evals);
      T = T0 * (kmax - func.evals) / kmax; // Decrease the temperature as evaluation proceeds
      for (i=0; i<3; i++)
	    {
	      r[i]=(pmax[i] - pmin[i]) * (kmax - func.evals) / kmax; // Decrease the radius of search as the evaluation proceeds
	    }
      printf("r[0] = %.6f, r[1] = %.6f, r[2] = %.6f\n",r[0],r[1],r[2]);
      newpoint(p, pmin, pmax, newp, r); // Find a new randomly chosen point within the radius of search
      printf("newp[0] = %.6f, newp[1] = %.6f, newp[2] = %.6f\n",newp[0],newp[1],newp[2]);
      newf = func(newp); // Evaluate the function at the new point
      printf("newf = %.6f\n",newf);
     
      randy = (double) rand()/RAND_MAX;
      pr = prob(f,newf,T);
      printf("Prob of new point = %.6f, randy = %.6f\n",pr,randy);
      if (pr > randy)
	{
	  f = newf;
	  for (i=0; i<3; i++)
	    {
	      p[i] = newp[i];
	    }
	}
      if (f < bestf)
	{
	  bestf = f;
	  for (i=0; i<3; i++)
	    {
	      bestp[i] = p[i];
	    }
	}
    printf("Current point,  dx=%.2f  dy=%.2f  theta=%.2f  fom=%.6f, bestfom=%.6f\n", p[0],p[1],p[2],f,bestf);

    }
  fptr = &bestf;
  printf("Best point found,  dx=%.2f  dy=%.2f  theta=%.2f  fom=%.6f\n", p[0],p[1],p[2],*fptr);
  printf("# of Evaluations=%d\n",func.evals);

  return 0;
}

//  SUBROUTINES:

double Func::shift()
//This shifts and rotates one data set to align with another. 
//Data1 is the larger, unshifted simulation.  Data2 is the smaller, 
//measured data. 

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
		  xprime=x*cos(align->d[2])+y*sin(align->d[2])+align->d[0];
		  yprime=-x*sin(align->d[2])+y*cos(align->d[2])+align->d[1];
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
  for (j=0; j<3; j++) align->d[j]=x[j];
  evals += 1;
  f=shift();
  return f;
}

void Func::df(double x[], double df[])
{
  int j;
  double xh[3];
  for (j=0; j<3; j++) xh[j]=x[j];
  double fold=f;
  for (j=0; j<3; j++)
    {
      double temp=x[j];
      double h=EPS*fabs(temp);
      if (h==0.0) h=EPS;
      xh[j]=temp+h;
      h=xh[j]-temp;
      double fh=operator()(xh);
      xh[j]=temp;
      df[j]=(fh-fold)/h;
    }
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

  for (i=0; i<3; i++)
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

void lnsrch(double xold[], const double fold, double g[], double p[], double x[], double &f, double xmin[], double xmax[], int &check, Func &func)
//This does a linear minimum search. It starts with the full Newton step, 
//then backtracks to find a minimum.  This is adapted from Numerical Recipes,
//but modified to include max and min tests on the input vector
//instead of a maximum step size.
{
  const double ALF=1.0E-6, TOLX=1.0E-9;
  double a,alam,alam2=0.0,alamin,b,disc,f2=0.0;
  double rhs1,rhs2,slope=0.0,temp,test,tmplam;
  int i,limflag,n=3;
  check=0;
  for (i=0; i<n; i++)
    slope += g[i]*p[i];
  if (slope>0.0) printf("Roundoff problem in lnsrch.");
  test=0.0;
  for (i=0; i<n; i++)
    {
      temp=fabs(p[i])/fmax(fabs(xold[i]),1.0);
      if (temp>test) test=temp;
    }
  alamin=TOLX/test;
  alam=1.0;
  for (;;)
    {
      limflag=0;
      printf("alam=%.9f, alamin=%.9f\n",alam,alamin);
      for (i=0; i<n; i++)
	{
	  x[i]=xold[i]+alam*p[i];
	  if (x[i]>=xmax[i] || x[i]<=xmin[i]) limflag=1;
	}
      if (limflag==1)
	{
	  for (i=0; i<n; i++) x[i]=xold[i];
	  f=fold;
	}
      else f=func(x);
      if (alam<alamin)
	{
	  for (i=0; i<n; i++) x[i]=xold[i];
	  check=1;
	  return;
	}
      else if (f<=fold+ALF*alam*slope) return;
      else
	{
	  if (alam==1.0)
	    tmplam=-slope/(2.0*(f-fold-slope));
	  else
	    {
	      rhs1=f-fold-alam*slope;
	      rhs2=f2-fold-alam2*slope;
	      a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	      b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
	      if (a==0.0) tmplam=-slope/(2.0*b);
	      else
		{
		  disc=b*b-3.0*a*slope;
		  if (disc<0.0) tmplam=0.5*alam;
		  else if (b<=0.0) tmplam=(-b+sqrt(disc))/(3.0*a);
		  else tmplam=-slope/(b+sqrt(disc));
		}
	      if (tmplam>0.5*alam)
		tmplam=0.5*alam;
	    }
	}
      alam2=alam;
      f2=f;
      alam=fmax(tmplam,0.1*alam);
    }
}

void dfpmin(double p[], double pmin[] ,double pmax[], const double gtol, int &iter, double &fret, Func &func)
//This finds the minimum of a multidimensional function.
//This is adapted from Numerical Recipes,
//but modified to include max and min tests on the input vector
//instead of a maximum step.
{
  const int ITMAX=200;
  const double EPS=1.0E-1;
  const double TOLX=1.0E-9;//4*EPS;
  int check;
  double den,fac,fad,fae,fp,sum=0.0,sumdg,sumxi,temp,test;
  int i,j,its,n=3;
  double dg[3], g[3], hdg[3], pnew[3], xi[3];
  double hessin[n*n];
  fp=func(p);
  func.df(p,g);
  for (i=0; i<n; i++)
    {
      for (j=0; j<n; j++) hessin[i+j*n]=0.0;
      hessin[i+i*n]=1.0;
      xi[i]=-g[i];
      sum+=p[i]*p[i];
    }
      for (its=0; its<ITMAX; its++)
	{
	  iter=its;
	  lnsrch(p,fp,g,xi,pnew,fret,pmin,pmax,check,func);
	  if (check==1) printf("LNSRCH FAILURE!\n");
	  fp=fret;
	  for (i=0; i<n; i++)
	    {
	      printf("i=%d, p[%d]=%.6f, pnew[%d]=%.6f, g[%d]=%.9f\n",i,i,p[i],i,pnew[i],i,g[i]);
	      xi[i]=pnew[i]-p[i];
	      p[i]=pnew[i];
	    }
	  printf("func(p) = %.6f\n",func(p));
	  test=0.0;
	  for (i=0; i<n; i++)
	    {
	      temp=fabs(xi[i])/fmax(fabs(p[i]),1.0);
	      if (temp>test) test=temp;
	    }
	  if (test<TOLX) return;
	  for (i=0; i<n; i++) dg[i]=g[i];
	  func.df(p,g);
	  test=0.0;
	  den=fmax(fabs(fret),1.0);
	  for (i=0; i<n; i++)
	    {
	      temp=fabs(g[i])*fmax(fabs(p[i]),1.0)/den;
	      if (temp>test) test=temp;
	    }
	  if (test<gtol) return;
	  for (i=0; i<n; i++) dg[i]=g[i]-dg[i];
	  for (i=0; i<n; i++)
	    {
	      hdg[i]=0.0;
	      for (j=0; j<n; j++) hdg[i]+=hessin[i+j*n]*dg[j];
	    }
	  fac=fae=sumdg=sumxi=0.0;
	  for (i=0; i<n; i++)
	    {
	      fac += dg[i]*xi[i];
	      fae += dg[i]*hdg[i];
	      sumdg += dg[i]*dg[i];
	      sumxi += xi[i]*xi[i];
	    }
	  if (fac>sqrt(EPS*sumdg*sumxi))
	    {
	      fac=1.0/fac;
	      fad=1.0/fae;
	      for (i=0; i<n; i++) dg[i]=fac*xi[i]-fad*hdg[i];
	      for (i=0; i<n; i++)
		{
		  for (j=i; j<n; j++)
		    {
		      hessin[i+j*n] += fac*xi[i]*xi[j] - fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
		      hessin[j+i*n]=hessin[i+j*n];
		    }
		}
	    }
	  for (i=0; i<n; i++)
	    {
	      xi[i]=0.0;
	      for (j=0; j<n; j++) xi[i] -= hessin[i+j*n]*g[j];
	    }
	}
      printf("too many iterations in dfpmin");
      return;
}

