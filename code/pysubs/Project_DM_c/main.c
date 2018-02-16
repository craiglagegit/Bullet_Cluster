#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

void DMProject(int numpart, struct particle_data *p, float *mass, float xmin, float xmax, float ymin, float ymax, float zmin, float zmax, int xpixels, int ypixels, float theta, float phi, float psi, float *value)
{
/* This code projects out the mass of dark matter in a given region, rotated by Euler angles  Craig Lage - 15-Jun-11. */
  NumPart = numpart;
  P = p;
  Mass = mass;
  Value = value;
  Xmin = xmin;
  Xmax = xmax;
  Ymin = ymin;
  Ymax = ymax;
  Zmin = zmin;
  Zmax = zmax;

  Xpixels = xpixels;
  Ypixels = ypixels;

  Theta = theta;
  Phi = phi;
  Psi = psi;

  printf("N=%d\n", NumPart);
  printf("projecting\n");  fflush(stdout);
  make_map_dm();
  printf("done\n");
}

void make_map_dm(void)
{
  float R[3][3]; /* These are the Euler angles used for rotating the particle coordinates. */
  R[0][0]=cos(Psi)*cos(Phi)-cos(Theta)*sin(Phi)*sin(Psi);
  R[0][1]=cos(Psi)*sin(Phi)+cos(Theta)*cos(Phi)*sin(Psi);
  R[0][2]=sin(Psi)*sin(Theta);
  R[1][0]=-sin(Psi)*cos(Phi)-cos(Theta)*sin(Phi)*cos(Psi);
  R[1][1]=-sin(Psi)*sin(Phi)+cos(Theta)*cos(Phi)*cos(Psi);
  R[1][2]=cos(Psi)*sin(Theta);
  R[2][0]=sin(Theta)*sin(Phi);
  R[2][1]=-sin(Theta)*cos(Phi);
  R[2][2]=cos(Theta);

  int i, j, n;
  double pixelsizeX, pixelsizeY, z;

  for(i = 0; i < Xpixels; i++)
    for(j = 0; j < Ypixels; j++)
      {
        Value[i + Xpixels * j] = 0;
      }
  pixelsizeX = (Xmax - Xmin) / Xpixels;
  pixelsizeY = (Ymax - Ymin) / Ypixels;

  for(n = 0; n < NumPart; n++)
    {
      if((n % (NumPart / 100)) == 0)
	{
	  printf(".");
	  fflush(stdout);
	}
      z = R[2][0] * P[n].Pos[0] + R[2][1] * P[n].Pos[1] + R[2][2] * P[n].Pos[2];
      if(z < Zmin || z > Zmax)
        continue;

      i = (int) (R[0][0] * P[n].Pos[0] + R[0][1] * P[n].Pos[1] + R[0][2] * P[n].Pos[2]-Xmin) / pixelsizeX;
      j = (int) (R[1][0] * P[n].Pos[0] + R[1][1] * P[n].Pos[1] + R[1][2] * P[n].Pos[2]-Ymin) / pixelsizeY;

      if (i >= 0 && i < Xpixels && j >= 0 && j < Ypixels)
          Value[i + Xpixels * j] = Value[i + Xpixels * j] + Mass[n];
    }

}

