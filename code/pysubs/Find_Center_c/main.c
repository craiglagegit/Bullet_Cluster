#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define max(x,y) (x>y?x:y)      // max macro definition
#define min(x,y) (x>y?y:x)      // min macro definition

#include "allvars.h"
#include "proto.h"

void Find_Center(int numpart, struct particle_data *p, float *hsml, float *mass, float xcen, float ycen, float zcen, float rmax, int desdensngb, float hmax, double boxsize, float *center)
{
/* These statements set the Global Variables (with initial capitals) to equal the variables which are passed into the function statement.  A little clumsy, but the easiest way to adapt the software that is already written to run in my python code.  Craig Lage - 17-Jun-10. */
/* Reworking this to create a 3D grid instead of projecting out a 2D grid.  Craig Lage - 26-Oct-10. */
/* Reworking to calculate a radial profile. Craig Lage - 26-May-11. */

  NumPart = numpart;
  P = p;
  Hsml = hsml;
  Mass = mass;

  DesDensNgb =  desdensngb;

  Hmax = hmax;
  Rmax = rmax;
  Xcen=xcen;
  Ycen = ycen;
  Zcen = zcen;
  BoxSize = boxsize;

  Center =         center;

  printf("N=%d\n", NumPart);

  printf("peano-hilbert order...\n");
  peano_hilbert_order();
  printf("done\n");

  tree_treeallocate(2.0 * NumPart, NumPart);
  Softening = (Rmax)/Rpixels/100;

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
  int i, n;
  int dR, nR;
  double h, r, u, wk;
  double pos[3];
  double sum, hmin, hmax, Rpart, Rcen, rr;
  double pixelsizeR;
 
  BoxHalf = 0.5 * BoxSize;

  for(i = 0; i < Rpixels; i++)
        {
          Value[i] = 0;
          ValueQuantity[i] = 0;
        }

  pixelsizeR = Rmax / Rpixels;

  hmin = 1.001 * pixelsizeR / 2;

  hmax = Hmax;

  for(n = 0; n < NumPart; n++)
    {
      if((n % (NumPart / 100)) == 0)
	{
	  printf(".");
	  fflush(stdout);
	}

      pos[0]= P[n].Pos[0]-Xcen;
      pos[1]= P[n].Pos[1]-Ycen;
      pos[2]= P[n].Pos[1]-Zcen;

      Rpart = sqrt(pow(pos[0],2)+pow(pos[1],2)+pow(pos[2],2)); /* Particle center */

      h = Hsml[n];
      
      if(h < hmin)
        h = hmin;

      if(h > hmax)
        h = hmax;

#ifdef PERIODIC /* PERIODIC IS NOT IMPLEMENTED - Lage - 26-May-11 */
      if((NGB_PERIODIC(Xc - P[n].Pos[Axis1]) + 0.5 * (Xmax - Xmin)) < -Hsml[n])
        continue;
      if((NGB_PERIODIC(Xc - P[n].Pos[Axis1]) - 0.5 * (Xmax - Xmin)) > Hsml[n])
        continue;
      
      if((NGB_PERIODIC(Yc - P[n].Pos[Axis2]) + 0.5 * (Ymax - Ymin)) < -Hsml[n])
        continue;
      if((NGB_PERIODIC(Yc - P[n].Pos[Axis2]) - 0.5 * (Ymax - Ymin)) > Hsml[n])
        continue;

      if((NGB_PERIODIC(Zc - P[n].Pos[Axis3]) + 0.5 * (Zmax - Zmin)) < -Hsml[n])
        continue;
      if((NGB_PERIODIC(Zc - P[n].Pos[Axis3]) - 0.5 * (Zmax - Zmin)) > Hsml[n])
        continue;
#else
      if(Rpart + h > Rmax)
        continue;
#endif
      nR = h / pixelsizeR + 1;

      /* Rcen is the center of the shell covered by the particle on the mesh */
      
      Rcen = floor(Rpart / pixelsizeR) * pixelsizeR;

      /* determine kernel normalization */
  
      sum = 0;      
      for(dR = -nR; dR <= nR; dR++)
            {
              r = abs(Rcen + dR * pixelsizeR - Rpart);
            
              if(r < h)
                {
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
      for(dR = -nR; dR <= nR; dR++)
              {
                rr = abs(Rcen + dR * pixelsizeR);

#ifdef PERIODIC /* PERIODIC IS NOT IMPLEMENTED - Lage - 26-May-11 */
                xxx = NGB_PERIODIC(xxx + Xmin - Xc) + Xc - Xmin;
                yyy = NGB_PERIODIC(yyy + Ymin - Yc) + Yc - Ymin;
                zzz = NGB_PERIODIC(zzz + Zmin - Zc) + Zc - Zmin;
#endif         
                    i = rr / pixelsizeR;                
                    if(i < Rpixels)
                          {
              		    r = abs(Rcen + dR * pixelsizeR - Rpart);            
              		    if(r < h)
                              {
                                u = r / h;                          
                                if(u < 0.5)
                                  wk = (2.546479089470 + 15.278874536822 * (u - 1) * u * u);
                                else
                                  wk = 5.092958178941 * (1.0 - u) * (1.0 - u) * (1.0 - u);
                                Value[i] += Mass[n] * wk / sum;
                                ValueQuantity[i] += Mass[n]*Quantity[n]*wk / sum;
                              }
                          }
              }
    }


  for(i = 0; i < Rpixels; i++)
        if(Value[i]>0)
          ValueQuantity[i] /= Value[i];
}

