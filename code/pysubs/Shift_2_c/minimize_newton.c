/*
Craig Lage - 25-Jul-10

This uses a Newton's method function minimization routine from Numerical Recipes to align two 2d functions

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#define pi 3.14159265358979     // Pi

//  DATA STRUCTURES:  

struct Align  
/*  This is the x, y, and theta that optimally aligns the two data sets
    dx=d[0], dy=d[1], theta=d[2].  This also contains the allowable limits
    on dx, dy, and theta.
 */
{
  double d[3], dmin[3], dmax[3];
};

struct Arr //This packages the 2D data sets.
{
  Arr(double, double, int, double, double, int);
  ~Arr();
  int nx, ny;
  double xmin, xmax, ymin, ymax, dx, dy, *x, *y, *data;
};

struct Func2 //This is the function to be minimized
{
  Arr *data1A, *data2A, *shifteddata1A, *data1B, *data2B, *shifteddata1B;
  Align *align;
  double EPS, weightA, weightB;
  double f; //Holds the most recent value of the function
  int evals; // Counts the number of function evaluations
  Func2(Arr *, Arr *, Arr *, Arr *, Arr *, Arr *, Align *); //Constructor
  double shift();
  double operator()(std::vector<double>&);
  void df(std::vector<double>&, std::vector<double>&);
};

// FUNCTION PROTOTYPES

void dfpmin(std::vector<double> &, std::vector<double> &, std::vector<double> &, const double, int &, double &, Func2 &);

void lnsrch(std::vector<double> &, const double, std::vector<double> &, std::vector<double> &, std::vector<double> &, double &, std::vector<double> &, std::vector<double> &, int &, Func2 &);

double kernel(double, double);
  
double interp(Arr*, double, double);

extern "C" int optimize2(Arr*, Arr*, Arr*, Arr*, Arr*, Arr*, Align*, double*,                            double, double, const double);

//THIS IS THE TOP LEVEL METHOD

extern "C" int optimize2(Arr* d1Aptr, Arr* sd1Aptr, Arr* d2Aptr,                                         Arr* d1Bptr, Arr* sd1Bptr, Arr* d2Bptr, Align* aptr,                            double* fptr, double weightA, double weightB,                                   const double gtol)
// This is a wrapper function to translate the Python call
// to use the appropriate c++ optimization routines.

{
  Func2 func2(d1Aptr,sd1Aptr,d2Aptr,d1Bptr,sd1Bptr,d2Bptr,aptr);

  func2.weightA=weightA;
  func2.weightB=weightB;

  std::vector<double> p(3), pmax(3), pmin(3);
  int i, iter;

  for (i=0; i<3; i++)
    {
      p[i]=func2.align->d[i];
      pmin[i]=func2.align->dmin[i];
      pmax[i]=func2.align->dmax[i];
    }

  dfpmin(p, pmin, pmax, gtol, iter, *fptr, func2);
  printf("Back from dfpmin, iter=%d, dx=%.2f  dy=%.2f  theta=%.2f  diff=%.6f\n", iter, p[0],p[1],p[2],*fptr);
  printf("Evaluations=%d\n",func2.evals);

  return 0;
}

//  SUBROUTINES:

Arr::Arr(double Xmin, double Xmax, int Nx, double Ymin, double Ymax, int Ny)//Array constructor
{
  xmin=Xmin; xmax=Xmax; ymin=Ymin; ymax=Ymax; nx=Nx; ny=Ny;
  data=new double[nx*ny];
  x = new double[nx];
  y = new double[ny];
  dx=(xmax-xmin)/nx;
  dy=(ymax-ymin)/ny;
  int i, j;
  for (i=0; i<nx; i++)
    {
      x[i] = xmin + dx/2 + i*dx;
      for(j=0; j<ny; j++)
	{
	  y[j] = ymin + dy/2 + j*dy;
	  data[i+j*nx]=0.0;
	}
    }
}

Arr::~Arr()//Array destructor
{
  delete[] data;
  delete[] x;
  delete[] y;
}

Func2::Func2(Arr* d1Aptr, Arr* sd1Aptr, Arr* d2Aptr,                                         Arr* d1Bptr, Arr* sd1Bptr, Arr* d2Bptr, Align* aptr)
// Func2 constructor
{
  evals=0;
  EPS=1.0E-4;
  data1A=d1Aptr;
  shifteddata1A=sd1Aptr;
  data2A=d2Aptr;
  data1B=d1Bptr;
  shifteddata1B=sd1Bptr;
  data2B=d2Bptr;
  align=aptr;
  f=0;
}

double Func2::shift()
//This shifts and rotates one data set to align with another. 
//Data1 is the larger, unshifted simulation.  Data2 is the smaller, 
//measured data.

{
  int i,j;
  double fom=0.0, xprime, yprime, x, y;
  for (i=0; i<data2A->nx; i++) 
    {
      x=data2A->x[i];
      for (j=0; j<data2A->ny; j++)
	{
	  y=data2A->y[j];
	  xprime=x*cos(align->d[2])+y*sin(align->d[2])+align->d[0];
	  yprime=-x*sin(align->d[2])+y*cos(align->d[2])+align->d[1];
	  shifteddata1A->data[i+j*data2A->nx]=interp(data1A,xprime,yprime);
	  fom = fom + weightA*pow((shifteddata1A->data[i+j*data2A->nx] -                        data2A->data[i+j*data2A->nx]),2);
	}
    }
  for (i=0; i<data2B->nx; i++) 
    {
      x=data2B->x[i];
      for (j=0; j<data2B->ny; j++)
	{
	  y=data2B->y[j];
	  xprime=x*cos(align->d[2])+y*sin(align->d[2])+align->d[0];
	  yprime=-x*sin(align->d[2])+y*cos(align->d[2])+align->d[1];
	  shifteddata1B->data[i+j*data2B->nx]=interp(data1B,xprime,yprime);
	  fom = fom + weightB*pow((shifteddata1B->data[i+j*data2B->nx] -                        data2B->data[i+j*data2B->nx]),2);
	}
    }
  return fom;
}

double Func2::operator()(std::vector<double> &x)
  // This is the function being minimized.
{
  int j, n=x.size();
  for (j=0; j<n; j++) align->d[j]=x[j];
  evals += 1;
  f=shift();
  //printf("f=%.10f\n",f);
  return f;
}

void Func2::df(std::vector<double> &x, std::vector<double> &df)
{
  int j, n=x.size();
  std::vector<double> xh=x;
  double fold=f;
  for (j=0; j<n; j++)
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

void lnsrch(std::vector<double> &xold, const double fold, std::vector<double> &g, std::vector<double> &p, std::vector<double> &x, double &f, std::vector<double> &xmin, std::vector<double> &xmax, int &check, Func2 &func)
//This does a linear minimum search. It starts with the full Newton step, 
//then backtracks to find a minimum.  This is adapted from Numerical Recipes,
//but modified to include max and min tests on the input vector
//instead of a maximum step size.
{
  const double ALF=1.0E-4, TOLX=1.0E-6;
  double a,alam,alam2=0.0,alamin,b,disc,f2=0.0;
  double rhs1,rhs2,slope=0.0,temp,test,tmplam;
  int i,limflag,n=xold.size();
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
      //printf("alam=%.9f\n",alam);
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

void dfpmin(std::vector<double> &p, std::vector<double> &pmin ,std::vector<double> &pmax, const double gtol, int &iter, double &fret, Func2 &func)
//This finds the minimum of a multidimensional function.
//This is adapted from Numerical Recipes,
//but modified to include max and min tests on the input vector
//instead of a maximum step.
{
  const int ITMAX=200;
  const double EPS=1.0E-4;
  const double TOLX=4*EPS;
  int check;
  double den,fac,fad,fae,fp,sum=0.0,sumdg,sumxi,temp,test;
  int i,j,its,n=p.size();
  std::vector<double> dg(n), g(n), hdg(n), pnew(n), xi(n);
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
	      xi[i]=pnew[i]-p[i];
	      p[i]=pnew[i];
	    }
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
		      hessin[i+j*n] += fac*xi[i]*xi[j]                                                               - fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
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

