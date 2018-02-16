#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

void findHsmlAndProject(int numpart, struct particle_data *p, float *hsml, float *mass, float *quantity, float xmin, float xmax, float ymin, float ymax, float zmin, float zmax, int xpixels, int ypixels, int desdensngb, int axis1, int axis2, int axis3, float hmax, double boxsize, float *value, float *valuequantity)
{
/* These statements set the Global Variables (with initial capitals) to equal the variables which are passed into the function statement.  A little clumsy, but the easiest way to adapt the softawre that is already written to run in my python code.  Craig Lage - 17-Jun-10. */

  NumPart = numpart;
  P = p;
  Hsml = hsml;
  Mass = mass;
  Quantity = quantity;

  Xmin = xmin;
  Xmax = xmax;
  Ymin = ymin;
  Ymax = ymax;
  Zmin = zmin;
  Zmax = zmax;

  Xpixels = xpixels;
  Ypixels = ypixels;

  DesDensNgb =  desdensngb;

  Axis1 = axis1;
  Axis2 = axis2;
  Axis3 = axis3;

  Hmax = hmax;

  BoxSize = boxsize;

  Value =         value;
  ValueQuantity = valuequantity;

  printf("N=%d\n", NumPart);

  printf("peano-hilbert order...\n");
  peano_hilbert_order();
  printf("done\n");

  tree_treeallocate(2.0 * NumPart, NumPart);
  Softening = (Xmax-Xmin)/Xpixels/100;

  printf("build tree...\n");
  tree_treebuild();
  printf("done.\n"); fflush(stdout);

  printf("finding neighbours...\n"); fflush(stdout);
  determine_hsml();
  printf("done.\n");

  tree_treefree();

  printf("projecting\n");  fflush(stdout);
  make_map();
  printf("done\n");
}

void determine_hsml(void)
{
  int i, signal;
  double h;

  for(i = 0, signal = 0, h = 0; i < NumPart; i++)
    {
      if(i > (signal / 100.0) * NumPart)
        {
          printf("x");
          fflush(stdout);
          signal++;
        }

      if(Hsml[i] == 0)
        Hsml[i] = h  = ngb_treefind(P[i].Pos, DesDensNgb, h * 1.1);
    }

  printf("\n");
}


#ifdef PERIODIC
#define NGB_PERIODIC(x) (((x)>BoxHalf)?((x)-BoxSize):(((x)<-BoxHalf)?((x)+BoxSize):(x)))
#else
#define NGB_PERIODIC(x) (x)
#endif



void make_map(void)
{
  int i, j, n;
  int dx, dy, nx, ny;
  double h, r, u, wk;
  double pos[2];
  double LengthX;
  double LengthY;
  double r2, h2;
  double sum, hmin, hmax, x, y, xx, yy, xxx, yyy;
  double pixelsizeX, pixelsizeY;

  BoxHalf = 0.5 * BoxSize;

  for(i = 0; i < Xpixels; i++)
    for(j = 0; j < Ypixels; j++)
      {
        Value[i + Xpixels * j] = 0;
        ValueQuantity[i + Xpixels * j] = 0;
      }
 
  LengthX = Xmax-Xmin;
  LengthY = Ymax-Ymin;

  Xc = 0.5 * (Xmax + Xmin);
  Yc = 0.5 * (Ymax + Ymin);

  pixelsizeX = LengthX / Xpixels;
  pixelsizeY = LengthY / Ypixels;

  if(pixelsizeX < pixelsizeY)
    hmin = 1.001 * pixelsizeX / 2;
  else
    hmin = 1.001 * pixelsizeY / 2;


  hmax = Hmax;


  for(n = 0; n < NumPart; n++)
    {
      if((n % (NumPart / 100)) == 0)
	{
	  printf(".");
	  fflush(stdout);
	}

      if(P[n].Pos[Axis3]< Zmin || P[n].Pos[Axis3] > Zmax)
        continue;

      pos[0]= P[n].Pos[Axis1]-Xmin;
      pos[1]= P[n].Pos[Axis2]-Ymin;


      h = Hsml[n];
      
      if(h < hmin)
        h = hmin;

      if(h > hmax)
        h = hmax;

#ifdef PERIODIC
      if((NGB_PERIODIC(Xc - P[n].Pos[Axis1]) + 0.5 * (Xmax - Xmin)) < -Hsml[n])
        continue;
      if((NGB_PERIODIC(Xc - P[n].Pos[Axis1]) - 0.5 * (Xmax - Xmin)) > Hsml[n])
        continue;
      
      if((NGB_PERIODIC(Yc - P[n].Pos[Axis2]) + 0.5 * (Ymax - Ymin)) < -Hsml[n])
        continue;
      if((NGB_PERIODIC(Yc - P[n].Pos[Axis2]) - 0.5 * (Ymax - Ymin)) > Hsml[n])
        continue;
#else
      if(pos[0] + h < 0 || pos[0] - h >  LengthX
         || pos[1] + h < 0 || pos[1] - h > LengthY)
        continue;
#endif
 
      h2 = h * h;

      nx = h / pixelsizeX + 1;
      ny = h / pixelsizeY + 1;

      /* x,y central pixel of region covered by the particle on the mesh */
      
      x = (floor(pos[0] / pixelsizeX) + 0.5) * pixelsizeX;
      y = (floor(pos[1] / pixelsizeY) + 0.5) * pixelsizeY;

      /* determine kernel normalization */

      
      sum = 0;

      
      for(dx = -nx; dx <= nx; dx++)
        for(dy = -ny; dy <= ny; dy++)
          {
            xx = x + dx * pixelsizeX - pos[0];
            yy = y + dy * pixelsizeY - pos[1];
            r2 = xx * xx + yy * yy;
            
            if(r2 < h2)
              {
                r = sqrt(r2);
                u = r / h;
                
                if(u < 0.5)
                  wk = (2.546479089470 + 15.278874536822 * (u - 1) * u * u);
                else
                  wk = 5.092958178941 * (1.0 - u) * (1.0 - u) * (1.0 - u);
                
                sum += wk;
              }
          }
      
      if(sum < 1.0e-10)
        continue;

      for(dx = -nx; dx <= nx; dx++)
        for(dy = -ny; dy <= ny; dy++)
          {
            xxx = x + dx * pixelsizeX;
            yyy = y + dy * pixelsizeY;

#ifdef PERIODIC
            xxx = NGB_PERIODIC(xxx + Xmin - Xc) + Xc - Xmin;
            yyy = NGB_PERIODIC(yyy + Ymin - Yc) + Yc - Ymin;
#endif
           
            if(xxx >= 0 && yyy >= 0)
              {
                i = xxx / pixelsizeX;
                j = yyy / pixelsizeY;
                
                if(i >= 0 && i < Xpixels)
                  if(j >= 0 && j < Ypixels)
                    {
                      xx = x + dx * pixelsizeX - pos[0];
                      yy = y + dy * pixelsizeY - pos[1];
                      r2 = xx * xx + yy * yy;
                      
                      if(r2 < h2)
                        {
                          r = sqrt(r2);
                          u = r / h;
                          
                          if(u < 0.5)
                            wk = (2.546479089470 + 15.278874536822 * (u - 1) * u * u);
                          else
                            wk = 5.092958178941 * (1.0 - u) * (1.0 - u) * (1.0 - u);
                          
                          Value[i + Xpixels * j] += Mass[n] * wk / sum;
                          ValueQuantity[i + Xpixels * j] += Mass[n]*Quantity[n]*wk / sum;
                        }
                    }
              }
          }
    }


  for(i = 0; i < Xpixels; i++)
    for(j = 0; j < Ypixels; j++)
      if(Value[i + Xpixels * j]>0)
        ValueQuantity[i + Xpixels * j] /= Value[i + Xpixels * j]; 
}

