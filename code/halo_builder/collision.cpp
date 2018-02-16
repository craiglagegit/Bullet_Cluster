/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, NYU CCPP
  Date: Feb 14, 2014

  Triaxial halo generation program

*/

//****************** collision.cpp **************************

#include "collision.h"

Collision::Collision(string inname) //Constructor                                                                                            
{
  // This reads in the data from the collision.cfg
  // file and sets up the various halos

  // First we read in the Halo information and initialize the halos
  Rhomin = 1.0E-7;
  Umin = 6.0E5;
  ReadHaloInfo(inname);

  // Then, we read in the grid parameters from the .enzo file
  enzofile = GetStringParam(inname, "enzofile", "AMRTest.enzo"); // Get enzo file name
  ReadEnzoGridInfo();

  // Then, we read in the collision parameters
  ReadCollisionInfo(inname);
  // Then, we read in the hdf filenames from the .enzo file
  NumHDFFiles = 13;
  HDFNames = new string[NumHDFFiles];
  ReadHDFfilenames();

  VelocityFieldGrid = GetIntParam(inname, "VelocityFieldGrid",0);
  VelocityFieldCoherenceLength = GetDoubleParam(inname, "VelocityFieldCoherenceLength",100.0);
  VelocityFieldSpectralIndex = GetDoubleParam(inname, "VelocityFieldSpectralIndex",1.67);
  VelocityFieldSeed = GetIntParam(inname, "VelocityFieldSeed",13);
  MagneticFieldGrid = GetIntParam(inname, "MagneticFieldGrid",0);
  MagneticFieldCoherenceLength = GetDoubleParam(inname, "MagneticFieldCoherenceLength",100.0);
  MagneticFieldSpectralIndex = GetDoubleParam(inname, "MagneticFieldSpectralIndex",1.67);
  MagneticFieldSeed = GetIntParam(inname, "MagneticFieldSeed",13);
  MagneticFieldStrength = GetDoubleParam(inname, "MagneticFieldStrength",10.0E-6);
  // Now we write out the HDF hydro fields and particle data

  VField = new RandomField(this, 0);
  BField = new RandomField(this, 1);
  WriteHydroFiles();
  WriteParticleFiles();
  return;
}

Collision::~Collision() //Destructor                                                                                            
{

  int i;
  for (i=0; i<NumHalos; i++)
    {
      delete Halolist[i];
    }
  for (i=0; i<NumGrids; i++)
    {
      delete[] GridDimension[i];
      delete[] GridLeftEdge[i];
      delete[] GridRightEdge[i];
      delete[] GridDelta[i];
      delete[] GridSpan[i];
    }
  delete[] NDM;
  delete[] GridLevel;
  delete[] GridDimension;
  delete[] GridLeftEdge;
  delete[] GridRightEdge;
  delete[] GridDelta;
  delete[] GridSpan;
  delete[] PixelVolume;
  delete[] Halolist;
  delete[] HDFNames;
  delete VField;
  delete BField;
}

void Collision::ReadHaloInfo(string inname)
{
  int i;
  string filename, halonum;
  // First, we read in the halo files and initialize the halos
  NumHalos = GetIntParam(inname, "NumHalos", 1);
  Halolist = new Halo*[NumHalos];
  printf("\n############################\n\n");
  printf("NumHalos = %d\n",NumHalos);
  printf("\n############################\n\n");
  for (i=0; i<NumHalos; i++)
    {
      halonum = boost::lexical_cast<std::string>(i);
      filename = GetStringParam(inname,"HaloConfigFile["+halonum+"]", "triaxial.cfg"); //Halo config file name
      printf("Halo configuration %d filename is %s\n",i, filename.c_str());
      Halolist[i]= new Halo(filename); //Read configuration file and initialize parameters
      printf("Halo %d, P=%f, Q=%f, G = %f, R200 = %f, Rs = %f, M200 = %f NDM = %d, Ngas = %d, NTot = %d\n",i, Halolist[i]->P,Halolist[i]->Q,Halolist[i]->G,Halolist[i]->R200,Halolist[i]->Rs,Halolist[i]->M200,Halolist[i]->NDM,Halolist[i]->NGas,Halolist[i]->NTot);
      printf("\n############################\n\n");
    }
  return;
}

void Collision::ReadCollisionInfo(string inname)
{
  // Reads in the shifts and rotations for the halos
  int i;
  string halonum;
  for (i=0; i<NumHalos; i++)
    {
      halonum = boost::lexical_cast<std::string>(i);
      Halolist[i]->Xoffset[0] = GetDoubleParam(inname, "Deltax["+halonum+"]", 0.0);
      Halolist[i]->Xoffset[1] = GetDoubleParam(inname, "Deltay["+halonum+"]", 0.0);
      Halolist[i]->Xoffset[2] = GetDoubleParam(inname, "Deltaz["+halonum+"]", 0.0);
      Halolist[i]->Voffset[0] = GetDoubleParam(inname, "Deltavx["+halonum+"]", 0.0);
      Halolist[i]->Voffset[1] = GetDoubleParam(inname, "Deltavy["+halonum+"]", 0.0);
      Halolist[i]->Voffset[2] = GetDoubleParam(inname, "Deltavz["+halonum+"]", 0.0);
      Halolist[i]->Phi = GetDoubleParam(inname, "Phi["+halonum+"]", 0.0) * pi / 180.0;
      Halolist[i]->Theta = GetDoubleParam(inname, "Theta["+halonum+"]", 0.0) * pi / 180.0;
      Halolist[i]->Psi = GetDoubleParam(inname, "Psi["+halonum+"]", 0.0) * pi / 180.0;
      Halolist[i]->EulerInitialize();

      printf("Halo %d, Xoffset = (%f, %f, %f), Voffset = (%f, %f, %f)\n",i, Halolist[i]->Xoffset[0],Halolist[i]->Xoffset[1],Halolist[i]->Xoffset[2],Halolist[i]->Voffset[0],Halolist[i]->Voffset[1],Halolist[i]->Voffset[2]);
      printf("Halo %d, Phi = %f, Theta = %f, Psi = %f\n",i,Halolist[i]->Phi*180.0/pi,Halolist[i]->Theta*180.0/pi,Halolist[i]->Psi*180.0/pi);
    }
  printf("\n############################\n\n");
  return;
}


void Collision::ReadEnzoGridInfo()
{
  int i, j;
  string gridnum;
  NumGrids = GetIntParam(enzofile, "CosmologySimulationNumberOfInitialGrids", 1);
  printf("Number of Enzo Grids is %d\n",NumGrids);

  GridDimension = new int*[NumGrids];
  GridLeftEdge = new double*[NumGrids];
  GridRightEdge = new double*[NumGrids];
  GridSpan = new double*[NumGrids];
  GridDelta = new double*[NumGrids];
  GridLevel = new int[NumGrids];
  PixelVolume = new double[NumGrids];
  NDM = new int[NumGrids];
  HighestGridLevel = 0;

  for (i=0; i<NumGrids; i++)
    {
      gridnum = boost::lexical_cast<std::string>(i);
      NDM[i] = 0;
      GridDimension[i] = new int[3];
      GridLeftEdge[i] = new double[3];
      GridRightEdge[i] = new double[3];
      GridSpan[i] = new double[3];
      GridDelta[i] = new double[3];
      for (j=0; j<3; j++)
	{
	  GridDimension[i][j] = 1;
	  GridLeftEdge[i][j] = -1.0;
	  GridRightEdge[i][j] = 1.0;
	}
      if (i == 0)
	{
	  GridDimension[i] = GetIntList(enzofile, "TopGridDimensions", 3, GridDimension[i]);
	  GridLeftEdge[i] = GetDoubleList(enzofile, "DomainLeftEdge", 3, GridLeftEdge[i]);
	  GridRightEdge[i] = GetDoubleList(enzofile, "DomainRightEdge", 3, GridRightEdge[i]);
	  GridLevel[i] = 0;
	}
      else
	{
	  GridDimension[i] = GetIntList(enzofile, "CosmologySimulationGridDimension["+gridnum+"]", 3, GridDimension[i]);
	  GridLeftEdge[i] = GetDoubleList(enzofile, "CosmologySimulationGridLeftEdge["+gridnum+"]", 3, GridLeftEdge[i]);
	  GridRightEdge[i] = GetDoubleList(enzofile, "CosmologySimulationGridRightEdge["+gridnum+"]", 3, GridRightEdge[i]);
	  GridLevel[i] = GetIntParam(enzofile, "CosmologySimulationGridLevel["+gridnum+"]", 1);
	  if (GridLevel[i] > HighestGridLevel) 
	    {
	      HighestGridLevel = GridLevel[i];
	      HighestGridIndex = i;
	    }
	}
      PixelVolume[i] = 1.0;
      for (j=0; j<3; j++)
	{
	  GridDelta[i][j] = (GridRightEdge[i][j] - GridLeftEdge[i][j]) / (double)GridDimension[i][j];
	  GridSpan[i][j] = GridRightEdge[i][j] - GridLeftEdge[i][j];
	  PixelVolume[i] *= GridDelta[i][j];
	}
      printf("Grid Dimension %d = %d, %d, %d\n",i,GridDimension[i][0],GridDimension[i][1],GridDimension[i][2]);
      printf("Grid Left Edge %d = %f, %f, %f\n",i,GridLeftEdge[i][0],GridLeftEdge[i][1],GridLeftEdge[i][2]);
      printf("Grid Right Edge %d = %f, %f, %f\n",i,GridRightEdge[i][0],GridRightEdge[i][1],GridRightEdge[i][2]);
      printf("PixelVolume %d = %f\n",i,PixelVolume[i]);
      printf("\n############################\n\n");
    }
  return;
}

void Collision::ReadHDFfilenames()
{
  int i;
  string names[13] = {"Density", "Velocity1", "Velocity2", "Velocity3", "Bfield1", "Bfield2", "Bfield3", "PhiField", "TotalEnergy", "GasEnergy", "ParticlePosition", "ParticleVelocity", "ParticleMass"};
  for (i=0; i<NumHDFFiles; i++)
    {
       HDFNames[i] = GetStringParam(enzofile, "CosmologySimulation"+names[i]+"Name", "CosmologySimulation"+names[i]+"Name");
    }
  return;
}

void Collision::WriteHydroFiles()
{
  // This writes the hydro data to the HDF files
  // This is the meat of where the physical quantities get calculated
  int i, j, k, jj, n, h, f, NX, NY, NZ, NTot;
  int NumFields = 10; // Number of hydro fields
  double x[3], v[3], cmv[3], b[3], R, rho, rhoh, rhomax = 0.0, tote, inte, phi;
  double rvfactor, fntp, cmvfactor;
  string gridnum,  dot = ".", hdfname;
  int* attr_data1  = new int[1];
  int* attr_data3  = new int[3];
  int attr_dims = 1;

  // Find peak gas density - used for magnetic field scaling
  for (h=0; h<NumHalos; h++)
    {
      if (Halolist[h]->Rho(0.0) * Halolist[h]->Rhofactor > rhomax)
	{
	  rhomax = Halolist[h]->Rho(0) * Halolist[h]->Rhofactor;
	}
    }

  for (n=0; n<NumGrids; n++)
    {
      double** data = new double*[NumFields];
      gridnum = boost::lexical_cast<std::string>(n);
      NX = GridDimension[n][0]; NY = GridDimension[n][1]; NZ = GridDimension[n][2];
      NTot = NX * NY * NZ;
      for (f=0; f<NumFields; f++)
	{
	  data[f] = new double[NTot]; 
	}
      printf("Generating Hydro files for Grid %d",n);
      for (i=0; i<NX; i++)
	{
	  printf(".");
	  fflush(stdout);
	  x[0] = GridLeftEdge[n][0] + GridDelta[n][0] * ((double)i + 0.5);
	  for (j=0; j<NY; j++)
	    {
	      x[1] = GridLeftEdge[n][1] + GridDelta[n][1] * ((double)j + 0.5);
	      for (k=0; k<NZ; k++)
		{
		  x[2] = GridLeftEdge[n][2] + GridDelta[n][2] * ((double)k + 0.5);
		  rho = 0.0;
		  for (jj=0; jj<3; jj++)
		    {
		      cmv[jj] = 0.0;
		    }
		  rvfactor = 0.0;
		  phi = 0.0;
		  tote = 0.0;
		  inte = 0.0;
		  
		  for (h=0; h<NumHalos; h++)
		    {
		      R = Halolist[h]->Rxyz(x);
		      if (R < Halolist[h]->C) cmvfactor = 1.0;
		      else cmvfactor = exp(-R + Halolist[h]->C);
		      fntp = Halolist[h]->Fntp(R);
		      for (jj=0; jj<3; jj++)
			{
			  cmv[jj] += Halolist[h]->Voffset[0] * cmvfactor;
			}
		      rhoh = Halolist[h]->Rho(R) * Halolist[h]->Rhofactor;
		      rvfactor += sqrt(Halolist[h]->U(R) * fntp * Halolist[h]->Ufactor) * rhoh;
		      inte += Halolist[h]->U(R) * (1.0 - fntp) * Halolist[h]->Ufactor * rhoh;
		      rho += rhoh;
		    }
		  if (rho < Rhomin)
		    {
		      rho = Rhomin;
		      inte = Umin;
		      rvfactor = 0.0;
		    }
		  else
		    {
		      inte /= rho;
		      rvfactor /= rho;
		    }
		  VField->GetFieldValue(x, n, GridLevel[n], v);	      
		  BField->GetFieldValue(x, n, GridLevel[n], b);	      
		  for (jj=0; jj<3; jj++)
		    {
		      v[jj] = v[jj] * rvfactor + cmv[jj];
		      b[jj] = b[jj] * MagneticFieldStrength * pow(rho / rhomax, 2.0/3.0);
		    }
		  tote = inte + (v[0]*v[0] + v[1]*v[1] + v[2]*v[2])/2.0 + (b[0]*b[0] + b[1]*b[1] + b[2]*b[2]) / (2.0 * rho);
		  data[0][k*NX*NY + j*NX + i] = rho;
		  data[1][k*NX*NY + j*NX + i] = v[0];
		  data[2][k*NX*NY + j*NX + i] = v[1];
		  data[3][k*NX*NY + j*NX + i] = v[2];
		  data[4][k*NX*NY + j*NX + i] = b[0];
		  data[5][k*NX*NY + j*NX + i] = b[1];
		  data[6][k*NX*NY + j*NX + i] = b[2];
		  data[7][k*NX*NY + j*NX + i] = phi;
		  data[8][k*NX*NY + j*NX + i] = tote;
		  data[9][k*NX*NY + j*NX + i] = inte;
		}
	    }
	}
      printf("\n");
      //printf("n = %d, maxv = %f, minv = %f\n",n,maxv,minv);
      for (f=0; f<NumFields; f++)
	{
	  // First we write the data
	  if (NumGrids == 1)  hdfname = HDFNames[f];
	  else  hdfname = HDFNames[f]+dot+gridnum;
	  WriteHDF5File2(hdfname, hdfname, 1, NX*NY*NZ, data[f]);
	  delete[] data[f];
	  // Now we write the attributes
	  attr_data1[0] = 1;
	  attr_dims = 1;
	  WriteHDF5Attribute(hdfname, hdfname, "Component_Rank", attr_dims, attr_data1);
	  attr_data1[0] = NX*NY*NZ;
	  WriteHDF5Attribute(hdfname, hdfname, "Component_Size", attr_dims, attr_data1);
	  attr_data3[0] = NX; attr_data3[1] = NY; attr_data3[2] = NZ; 
	  attr_dims = 3;
	  WriteHDF5Attribute(hdfname, hdfname, "Dimensions", attr_dims, attr_data3);
	  attr_data1[0] = 3;
	  attr_dims = 1;
	  WriteHDF5Attribute(hdfname, hdfname, "Rank", attr_dims, attr_data1);
	  attr_data3[0] = GridDimension[0][0]; attr_data3[1] = GridDimension[0][1]; attr_data3[2] = GridDimension[0][2];
	  attr_dims = 3;
	  WriteHDF5Attribute(hdfname, hdfname, "TopGridDims", attr_dims, attr_data3);
	  for (jj=0; jj<3; jj++)
	    {
	      attr_data3[jj] = (int)round((GridLeftEdge[n][jj] - GridLeftEdge[0][jj]) / GridDelta[0][jj]);
	    }
	  WriteHDF5Attribute(hdfname, hdfname, "TopGridStart", attr_dims, attr_data3);
	  for (jj=0; jj<3; jj++)
	    {
	      attr_data3[jj] = (int)round((GridRightEdge[n][jj] - GridLeftEdge[0][jj]) / GridDelta[0][jj]) - 1;
	    }
	  WriteHDF5Attribute(hdfname, hdfname, "TopGridEnd", attr_dims, attr_data3);
	}
      delete[] data;
    }
  delete[] attr_data1;
  delete[] attr_data3;
  return;
}

void Collision::WriteParticleFiles()
{
  // This writes the particle data to the HDF files
  int j, k, n, h, NX, NY, NZ, counter;
  string gridnum,  dot = ".", hdfname;
  int* attr_data1  = new int[1];
  int* attr_data3  = new int[3];
  int attr_dims = 1;

  // First we put the particles into the appropriate grids, counting how many in each grid
  for (h=0; h<NumHalos; h++) // h is the halo counter
    {
      for (k=0; k<Halolist[h]->NDM; k++) // k is the particle counter
	{
	  // First we rotate and shift the position and velocity data for this particle
	  Halolist[h]->EulerRotateShift(k);
	  for (n=NumGrids-1; n>-1; n--) //n is the grid counter.  
	    {
	      //We step from the smallest grid outwards checking if the particle is in this grid
	      Halolist[h]->Part[k].grid = n;  // Assume it's in this grid until we find otherwise
	      for (j=0; j<3; j++)
		{
		  if (Halolist[h]->Part[k].x[j] < GridLeftEdge[n][j] || Halolist[h]->Part[k].x[j] > GridRightEdge[n][j]) 
		    {
		      Halolist[h]->Part[k].grid = -1; // It's not in this grid
		      break; // If it's outside one coordinate, no need to test the others
		    }
		}
	      if (Halolist[h]->Part[k].grid == n) // If it's in this grid
		{
		  NDM[n] += 1;
		  Halolist[h]->Part[k].mass /= PixelVolume[n]; // Scale mass by pixel volume
		  break; //Leave the loop once you've found the smallest grid that contains the particle
		}
	    }
	}
    }
  for (n=0; n<NumGrids; n++)
    {
      printf("Number of particles in Grid  %d is %d\n",n, NDM[n]);
    }

  // Now we collect the data and write the HDF files
  for (n=0; n<NumGrids; n++)
    {
      NX = GridDimension[n][0]; NY = GridDimension[n][1]; NZ = GridDimension[n][2];
      gridnum = boost::lexical_cast<std::string>(n);
      double* massdata = new double[NDM[n]];
      double* posdata = new double[3*NDM[n]];
      double* veldata = new double[3*NDM[n]];
      counter = 0;
      for (h=0; h<NumHalos; h++) // h is the halo counter
	{
	  for (k=0; k<Halolist[h]->NDM; k++) // k is the particle counter
	    {
	      if (Halolist[h]->Part[k].grid == n) // If it's in this grid
		{
		  massdata[counter] = Halolist[h]->Part[k].mass;
		  for (j=0; j<3; j++)
		    {
		      posdata[counter + NDM[n] * j] = Halolist[h]->Part[k].x[j];
		      veldata[counter + NDM[n] * j] = Halolist[h]->Part[k].v[j];
		    }
		  counter++;
		}
	    }
	}
      // First the Position
      if (NumGrids == 1)  hdfname = HDFNames[10];
      else  hdfname = HDFNames[10]+dot+gridnum;
      WriteHDF5File2(hdfname, hdfname, 3, NDM[n], posdata);
      attr_dims = 1;
      attr_data1[0] = 3;
      WriteHDF5Attribute(hdfname, hdfname, "Component_Rank", attr_dims, attr_data1);
      attr_data1[0] = NDM[n];
      WriteHDF5Attribute(hdfname, hdfname, "Component_Size", attr_dims, attr_data1);
      attr_data1[0] = NDM[n];
      WriteHDF5Attribute(hdfname, hdfname, "Dimensions", attr_dims, attr_data1);
      attr_data1[0] = 1;
      WriteHDF5Attribute(hdfname, hdfname, "Rank", attr_dims, attr_data1);
      attr_data3[0] = GridDimension[0][0]; attr_data3[1] = GridDimension[0][1]; attr_data3[2] = GridDimension[0][2];
      attr_dims = 3;
      WriteHDF5Attribute(hdfname, hdfname, "TopGridDims", attr_dims, attr_data3);
      attr_data3[0] = 0; attr_data3[1] = 0; attr_data3[2] = 0; 
      WriteHDF5Attribute(hdfname, hdfname, "TopGridStart", attr_dims, attr_data3);
      attr_data3[0] = NX - 1; attr_data3[1] = NY - 1; attr_data3[2] = NZ - 1; 
      WriteHDF5Attribute(hdfname, hdfname, "TopGridEnd", attr_dims, attr_data3);

      // Then the Velocity
      if (NumGrids == 1)  hdfname = HDFNames[11];
      else  hdfname = HDFNames[11]+dot+gridnum;
      WriteHDF5File2(hdfname, hdfname, 3, NDM[n], veldata);
      attr_dims = 1;
      attr_data1[0] = 3;
      WriteHDF5Attribute(hdfname, hdfname, "Component_Rank", attr_dims, attr_data1);
      attr_data1[0] = NDM[n];
      WriteHDF5Attribute(hdfname, hdfname, "Component_Size", attr_dims, attr_data1);
      attr_data1[0] = NDM[n];
      WriteHDF5Attribute(hdfname, hdfname, "Dimensions", attr_dims, attr_data1);
      attr_data1[0] = 1;
      WriteHDF5Attribute(hdfname, hdfname, "Rank", attr_dims, attr_data1);
      attr_data3[0] = GridDimension[0][0]; attr_data3[1] = GridDimension[0][1]; attr_data3[2] = GridDimension[0][2];
      attr_dims = 3;
      WriteHDF5Attribute(hdfname, hdfname, "TopGridDims", attr_dims, attr_data3);
      attr_data3[0] = 0; attr_data3[1] = 0; attr_data3[2] = 0; 
      WriteHDF5Attribute(hdfname, hdfname, "TopGridStart", attr_dims, attr_data3);
      attr_data3[0] = NX - 1; attr_data3[1] = NY - 1; attr_data3[2] = NZ - 1; 
      WriteHDF5Attribute(hdfname, hdfname, "TopGridEnd", attr_dims, attr_data3);

      // Finally the Masses
      if (NumGrids == 1)  hdfname = HDFNames[12];
      else  hdfname = HDFNames[12]+dot+gridnum;
      WriteHDF5File2(hdfname, hdfname, 1, NDM[n], massdata);
      attr_dims = 1;
      attr_data1[0] = 1;
      WriteHDF5Attribute(hdfname, hdfname, "Component_Rank", attr_dims, attr_data1);
      attr_data1[0] = NDM[n];
      WriteHDF5Attribute(hdfname, hdfname, "Component_Size", attr_dims, attr_data1);
      attr_data1[0] = NDM[n];
      WriteHDF5Attribute(hdfname, hdfname, "Dimensions", attr_dims, attr_data1);
      attr_data1[0] = 1;
      WriteHDF5Attribute(hdfname, hdfname, "Rank", attr_dims, attr_data1);
      attr_data3[0] = GridDimension[0][0]; attr_data3[1] = GridDimension[0][1]; attr_data3[2] = GridDimension[0][2];
      attr_dims = 3;
      WriteHDF5Attribute(hdfname, hdfname, "TopGridDims", attr_dims, attr_data3);
      attr_data3[0] = 0; attr_data3[1] = 0; attr_data3[2] = 0; 
      WriteHDF5Attribute(hdfname, hdfname, "TopGridStart", attr_dims, attr_data3);
      attr_data3[0] = NX - 1; attr_data3[1] = NY - 1; attr_data3[2] = NZ - 1; 
      WriteHDF5Attribute(hdfname, hdfname, "TopGridEnd", attr_dims, attr_data3);

      delete[] massdata;
      delete[] posdata;
      delete[] veldata;
    }
  delete[] attr_data1;
  delete[] attr_data3;
  return;
}


Collision::RandomField::RandomField(Collision* Coll, int Type)
{
  // This generates a Gaussian random field with Kolmogorov spectrum
  // Type = 0 => Velocity Field, Type = 1 => Magnetic Field
  int i;
  double LMin = 1.0E12;
  int NMin = 1E6;
  Field = new double*[3];
  if (Type ==0)
    {
      FieldGrid = Coll->VelocityFieldGrid;
      SpectralIndex = Coll->VelocityFieldSpectralIndex;
      Seed = Coll->VelocityFieldSeed;
      Clean_Divergence = false;
      alpha = 0.0;
      kmin = 2.0 * pi / Coll->VelocityFieldCoherenceLength;
    }
  else if (Type == 1)
    {
      FieldGrid = Coll->MagneticFieldGrid;
      SpectralIndex = Coll->MagneticFieldSpectralIndex;
      Seed = Coll->MagneticFieldSeed;
      Clean_Divergence = true;
      alpha = 0.0;
      kmin = 2.0 * pi / Coll->MagneticFieldCoherenceLength;
    }
  else
    {
      printf(" Invalid Type specification in Random Field generation. Quitting\n");
      exit(0);
    }
  HighestGridLevel = Coll->HighestGridLevel;
  GridSizeRatio = (int) pow(2.0, Coll->HighestGridLevel - Coll->GridLevel[FieldGrid]);
  NTot = 1;
  for (i=0; i<3; i++)
    {
      FieldSpan[i] = Coll->GridSpan[FieldGrid][i];
      if (FieldSpan[i] < LMin) LMin = FieldSpan[i];
      FieldDimension[i] = Coll->GridDimension[FieldGrid][i] * GridSizeRatio;
      if (FieldDimension[i] < NMin) NMin = FieldDimension[i];
      NTot *= FieldDimension[i];
      FieldLeftEdge[i] = Coll->GridLeftEdge[FieldGrid][i];
      FieldDelta[i] = Coll->GridDelta[Coll->HighestGridIndex][i];
    }
  printf("FieldDimension = (%d, %d, %d)\n",FieldDimension[0],FieldDimension[1],FieldDimension[2]);
  if (NTot > 1024 * 512 * 512)
    {
      printf(" Initial grid exceeds memory allocation. Quitting\n");
      exit(0);
    }
  for (i=0; i<3; i++)
    {
      Field[i] = new double[NTot];
    }
  kmax = 2.0 * pi / LMin * NMin / 2.0;
  garFields(FieldDimension, FieldSpan, 1.0, SpectralIndex, kmin, kmax, Seed, Clean_Divergence, alpha, Field[0], Field[1], Field[2]);
  return;
}

Collision::RandomField::~RandomField()
{
  int i;
  for (i=0; i<3; i++)
    {
      delete[] Field[i];
    }
  delete[] Field;
}

void Collision::RandomField::GetFieldValue(double* x, int Grid, int GridLevel, double* f)
{
  int i, j, k, ii, jj, kk, jjj, gridcounter, index;
  int GridRatio = (int) pow(2.0, HighestGridLevel - GridLevel);
  double* xp = new double[3];
  
  for (j=0; j<3; j++)
    {
      xp[j] = x[j];
      f[j] = 0;
    }
  if (Grid < FieldGrid)
    {
      return;
    }
  else
    {
      gridcounter = 0;
      for (i=0; i<HighestGridLevel - GridLevel; i++)
	{
	  for (j=0; j<3; j++)
	    {
	      xp[j] -= (pow(2.0, i) * FieldDelta[j]) / 2.0;
	    }
	}
      i = (int) (xp[0] - FieldLeftEdge[0]) / FieldDelta[0];
      j = (int) (xp[1] - FieldLeftEdge[1]) / FieldDelta[1];
      k = (int) (xp[2] - FieldLeftEdge[2]) / FieldDelta[2];
      for (ii=i; ii<i+GridRatio; ii++)
	{
	  for (jj=j; jj<j+GridRatio; jj++)
	    {
	      for (kk=k; kk<k+GridRatio; kk++)
		{
		  index = kk + FieldDimension[2] * jj + FieldDimension[1] * FieldDimension[2] * ii;
		  for (jjj=0; jjj<3; jjj++)
		    {
		      f[jjj] += Field[jjj][index];
		    }
		  gridcounter++;
		}
	    }
	}
      for (j=0; j<3; j++)
	{
	  f[j] /= (double) gridcounter;
	}
    }
  delete[] xp;
  return;
}

