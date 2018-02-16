

extern int  NumPart;

extern float *Mass;

extern float Xmin, Ymin, Xmax, Ymax, Zmin, Zmax, Theta, Phi, Psi;

extern int Xpixels, Ypixels;

extern float *Value;

extern struct particle_data
{
  float Pos[3];			/*!< particle position at its current time */
}
*P;                            /*!< points to particles on this processor */







