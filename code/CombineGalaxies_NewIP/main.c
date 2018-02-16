/* code to combine two galaxies such that they collide
 * on a parabolic encounter 
 *
 * written by V. Springel, MPA
 * Modified by Craig Lage - NYU 27-Nov-12
 * Now uses an Impact Parameter (ip) rotated in the y-z plane
 * by an amount theta
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "globvars.h"




int main(int argc, char *argv[])
{
  double phi1, theta1, psi1, phi2, theta2, psi2;
  double ip, theta, rstart, velincrease;
  char gal_fname1[100], gal_fname2[100], gal_output[100];

  if(argc != 14)
    {
      printf("\nwrong number of arguments\n");
      printf
	("call with:\n\n<fname_gal1> <phi1> <theta1> <psi1>\n<fname_gal2> <phi2> <theta2> <psi2>\n<ip> <theta> <rstart>\n<fname_galout>\n\n");
      printf("(angles in degrees.)\n\n");
      exit(0);
    }

  strcpy(gal_fname1, argv[1]);
  phi1 = atof(argv[2]);
  theta1 = atof(argv[3]);
  psi1 = atof(argv[4]);

  strcpy(gal_fname2, argv[5]);
  phi2 = atof(argv[6]);
  theta2 = atof(argv[7]);
  psi2 = atof(argv[8]);

  ip = atof(argv[9]);
  theta = atof(argv[10]);
  rstart = atof(argv[11]);
  velincrease = atof(argv[12]);

  strcpy(gal_output, argv[13]);


  load_particles(gal_fname1, &gal[0], &header[0]);
  load_particles(gal_fname2, &gal[1], &header[1]);

  turn_galaxy(&gal[0], phi1, theta1, psi1);
  turn_galaxy(&gal[1], phi2, theta2, psi2);


  move_galaxies(&gal[0], &gal[1], ip, theta, rstart, velincrease);

  save_combined(gal_output, &gal[0], &gal[1]);

  return 0;
}
