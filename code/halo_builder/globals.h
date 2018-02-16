/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, NYU CCPP
  Date: Feb 14, 2014

  Triaxial halo generation program

*/

//****************** globals.h **************************

#define pi 3.14159265358979     // Pi

#define GRAVITY    6.672e-8
#define MSOL       1.989e33
#define CM_PER_MPC 3.085678e24
#define HUBBLE     0.72		// Hubble Constant in 100km/sec / Mpc

#define UMASS      ( MSOL       * 1e10 )
#define ULENGTH    (CM_PER_MPC / 1000.0 )
#define UVELOCITY  (             1e5   )
#define UTIME      ( ULENGTH / UVELOCITY )
#define UDENSITY   ( UMASS / ( ULENGTH * ULENGTH * ULENGTH ) )
#define UPRESSURE  ( UMASS / ( ULENGTH * UTIME * UTIME ) )
#define UENERGY    ( UMASS * UVELOCITY * UVELOCITY )

#define CM         ( 1.0 / ULENGTH )
#define GRAM       ( 1.0 / UMASS )
#define SEC        ( 1.0 / UTIME )
#define ERGS       ( GRAM * ( CM * CM ) / ( SEC * SEC ) )

#define keV        ( 1.602e-9 * ERGS )

#define GAMMA_MINUS1 0.666
#define PROTONMASS   1.6726e-24
#define BOLTZMANN    1.3806e-16
#define HYDROGEN_MASSFRAC 0.76
