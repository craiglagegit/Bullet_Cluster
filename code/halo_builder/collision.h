/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, NYU CCPP
  Date: Feb 14, 2014

  Triaxial halo generation program

*/

//****************** collision.h **************************

#include <stdio.h>            
#include <stdlib.h>           
#include <math.h>             
#include <string.h>
#include <vector>
#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/iter_find.hpp>
#include <sstream>
#include <iostream>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

using namespace std;

#include <globals.h>
#include <fileio.h>
#include <hdf5write.h>
#include <halo.h>
#include <garfields/garFields.h>

class Collision
{
 public:
  int NumHalos;              // Number of Halos being collided
  int NumGrids;              // Number of grids in the initial enzo file
  int* NDM;                  // Number of DM particles in each grid
  double Umin;               // Minimum internal energy
  double Rhomin;             // Minimum density
  Halo** Halolist;           // Array of Halos

  int** GridDimension;       // Array of grid dimensions
  int*  GridLevel;           // Level of Grid Nesting
  int HighestGridLevel;     // Deepest level of nesting
  int HighestGridIndex;     // Grid with deepest level of nesting
  double** GridLeftEdge;     // Array of grid left edge
  double** GridRightEdge;    // Array of grid right edge
  double** GridSpan;         // Array of grid width
  double** GridDelta;        // Array of grid dx
  double*  PixelVolume;      // Pixel volume of grid

  int VelocityFieldGrid;     // Which grid is used for velocity field
  double VelocityFieldCoherenceLength;
  double VelocityFieldSpectralIndex;
  int VelocityFieldSeed;

  int MagneticFieldGrid;     // Which grid is used for velocity field
  double MagneticFieldCoherenceLength;
  double MagneticFieldSpectralIndex;
  int MagneticFieldSeed;
  double MagneticFieldStrength;

  int NumHDFFiles;
  string* HDFNames;           // HDF file names
  string enzofile;	      // .enzo file name

  Collision() {};
  Collision(string);
  ~Collision();
  void ReadHaloInfo(string);
  void ReadCollisionInfo(string);
  void ReadEnzoGridInfo();
  void ReadHDFfilenames();
  void WriteHydroFiles();
  void WriteParticleFiles();

  class RandomField
  {
  public:
    int FieldGrid;
    int GridSizeRatio;
    int NTot;
    double** Field;
    double kmin;
    double kmax;
    double SpectralIndex;
    bool Clean_Divergence;
    double alpha;
    int Seed;
    int FieldDimension[3];
    double FieldSpan[3];
    double FieldLeftEdge[3];
    double FieldDelta[3];
    int HighestGridLevel;

    RandomField(Collision*, int Type);
    ~RandomField();
    void GetFieldValue(double*, int, int, double*);
  };

  RandomField* VField;
  RandomField* BField;

};

