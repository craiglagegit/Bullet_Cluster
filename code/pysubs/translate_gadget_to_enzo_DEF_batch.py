#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 8-Aug-11


#This program translates a Gadget snapshot file into enzo init files.
#Added a random magnetic field read-in

import lageconfig # system specific path information
import sys
sys.path.append(lageconfig.bulletpath)
import pysubs_nas_18Sep12 as pysubs
import BulletConstants
from pylab import *
from scipy.ndimage import gaussian_filter
import h5py

#****************MAIN PROGRAM*****************
# There are three inputs required:
#	(1) The Gadget snapshot file to be translated
#	(2) The Gadget parameter file used
#	(3) The Enzo parameter file to be run
#
# If using Cosmology (If ComovingCoordinates==1) then the program assumes a bounding box to be translated in 
# Gadget coordinates which is a cube with edge length equal to the Enzo parameter CosmologyComovingBoxSize, 
# with origin in the center of the cube. In this case, the data in Enzo will be output in a bounding box 
# with coordinates (0.0,0.0,0.0) to (1.0,1.0,1.0)
#
# If not using Cosmology (If ComovingCoordinates==0) then the bounding box will be assumed to be specified
# by the parameters DomainLeftEdge and DomainRightEdge in the Enzo parameter file.
#

cmd=sys.argv
MagScaleFactor = float(cmd[1]) # Generated random fields with peak value of 1.0. This multiplier sets the actual value in uG.
FieldMode = cmd[2] # Generated random fields with peak value of 1.0. This multiplier sets the actual value in uG.
Tbg = float(cmd[3]) # Generated random fields with peak value of 1.0. This multiplier sets the actual value in uG.
LogRhoMin = float(cmd[4]) # Minimum density outside the initial clusters

GadgetICFile  = 'collision.dat'
GadgetParameterFile  = '../../../../collision.param'
EnzoParameterFile   = '../../../../AMRTest.enzo'

snapfile = pysubs.ReadGadgetICFile(GadgetICFile)
Ngas = snapfile.header.Npart[0] # Number of SPH particles in Gadget input
Ndm = snapfile.header.Npart[1]  # Number of DM particles in Gadget input
gadparam = pysubs.ReadGadgetParameterFile(GadgetParameterFile)
enzparam = pysubs.ReadEnzoParameterFile(EnzoParameterFile)

# Calculate the conversion factors between Gadget and Enzo parameters
RhoCrit        = 1.8788e-29 # Critical Density * h^2 in g/cm^3
GNewton        = 6.67428e-8 # cm^3 g^-1 s^-2
cm_per_mpc     = 3.0857e24
MSun           = 1.9891e33 # Solar mass in grams
GadgetLength   = gadparam.UnitLength_in_cm
GadgetMass     = gadparam.UnitMass_in_g
GadgetVelocity = gadparam.UnitVelocity_in_cm_per_s
GadgetTime     = GadgetLength / GadgetVelocity
GadgetDensity  = GadgetMass / GadgetLength**3.0

if enzparam.ComovingCoordinates==1:
	EnzoLength     = enzparam.CosmologyComovingBoxSize * cm_per_mpc / enzparam.CosmologyHubbleConstantNow / (1 + enzparam.CosmologyInitialRedshift)
	EnzoDensity    = RhoCrit * enzparam.CosmologyOmegaMatterNow*(1 + enzparam.CosmologyInitialRedshift)**3
	EnzoMass       = EnzoDensity * EnzoLength**3.0
	EnzoTime       = 1 / sqrt(4 * pi * GNewton * RhoCrit * (1 + enzparam.CosmologyInitialRedshift)**3)
	EnzoVelocity   = EnzoLength / EnzoTime
	EnzoBField     = sqrt(4 * pi * EnzoDensity) * EnzoVelocity
	enzparam.CosmologySimulationGridLeftEdge[0]  = [0.0,0.0,0.0] # Default
	enzparam.CosmologySimulationGridRightEdge[0] = [1.0,1.0,1.0] # Default
	LengthConversion   = GadgetLength / EnzoLength
	GadgetBoxLeftEdge   = array([-0.5, -0.5, -0.5]) / LengthConversion  # Box Edge is 1.0 in Enzo coordinates - this puts it back into Gadget coordinates
	GadgetBoxRightEdge  = array([0.5, 0.5, 0.5]) / LengthConversion  # Box Edge is 1.0 in Enzo coordinates - this puts it back into Gadget coordinates
else:
	EnzoLength     = enzparam.LengthUnits
	EnzoDensity    = enzparam.DensityUnits
	EnzoMass       = EnzoDensity * EnzoLength**3.0
	EnzoTime       = enzparam.TimeUnits
	EnzoVelocity   = EnzoLength / EnzoTime
	EnzoBField     = sqrt(4 * pi * EnzoDensity) * EnzoVelocity
	enzparam.CosmologySimulationGridLeftEdge[0]  = enzparam.DomainLeftEdge
	enzparam.CosmologySimulationGridRightEdge[0] = enzparam.DomainRightEdge
	GadgetBoxLeftEdge   = array(enzparam.DomainLeftEdge)
	GadgetBoxRightEdge  = array(enzparam.DomainRightEdge)

LengthConversion   = GadgetLength / EnzoLength
TimeConversion     = GadgetTime / EnzoTime
VelocityConversion = GadgetVelocity / EnzoVelocity 
MassConversion     = GadgetMass / EnzoMass
DensityConversion  = GadgetDensity / EnzoDensity

print "Length Conversion Factor=%e\nTime Conversion Factor=%e\nVelocity Conversion Factor=%e\nMass Conversion Factor=%e\nDensity Conversion Factor=%e"%(LengthConversion, TimeConversion, VelocityConversion, MassConversion, DensityConversion)

GadgetBoxCenter = (GadgetBoxLeftEdge + GadgetBoxRightEdge) / 2.0
EnzoBoxCenter   = (array(enzparam.CosmologySimulationGridLeftEdge[0]) + array(enzparam.CosmologySimulationGridRightEdge[0])) / 2.0

#print GadgetBoxCenter, EnzoBoxCenter

NumGrids=enzparam.CosmologySimulationNumberOfInitialGrids
enzparam.CosmologySimulationGridDimension[0] = enzparam.TopGridDimensions
enzparam.CosmologySimulationGridLevel[0]     = 0
PixelVolume=zeros([NumGrids])

# Now to grid the Gadget data, once for each subgrid
for i in range(NumGrids): 

	if enzparam.ComovingCoordinates==1:
		xmin = GadgetBoxLeftEdge[0] + enzparam.CosmologySimulationGridLeftEdge[i][0] * (GadgetBoxRightEdge[0]-GadgetBoxLeftEdge[0])
		xmax = GadgetBoxLeftEdge[0] + enzparam.CosmologySimulationGridRightEdge[i][0] * (GadgetBoxRightEdge[0]-GadgetBoxLeftEdge[0])
		ymin = GadgetBoxLeftEdge[1] + enzparam.CosmologySimulationGridLeftEdge[i][1] * (GadgetBoxRightEdge[1]-GadgetBoxLeftEdge[1])
		ymax = GadgetBoxLeftEdge[1] + enzparam.CosmologySimulationGridRightEdge[i][1] * (GadgetBoxRightEdge[1]-GadgetBoxLeftEdge[1])
		zmin = GadgetBoxLeftEdge[2] + enzparam.CosmologySimulationGridLeftEdge[i][2] * (GadgetBoxRightEdge[2]-GadgetBoxLeftEdge[2])
		zmax = GadgetBoxLeftEdge[2] + enzparam.CosmologySimulationGridRightEdge[i][2] * (GadgetBoxRightEdge[2]-GadgetBoxLeftEdge[2])
	else:
		xmin = enzparam.CosmologySimulationGridLeftEdge[i][0] 
		xmax = enzparam.CosmologySimulationGridRightEdge[i][0] 
		ymin = enzparam.CosmologySimulationGridLeftEdge[i][1] 
		ymax = enzparam.CosmologySimulationGridRightEdge[i][1] 
		zmin = enzparam.CosmologySimulationGridLeftEdge[i][2] 
		zmax = enzparam.CosmologySimulationGridRightEdge[i][2] 

	nx   = enzparam.CosmologySimulationGridDimension[i][0]
	ny   = enzparam.CosmologySimulationGridDimension[i][1]
	nz   = enzparam.CosmologySimulationGridDimension[i][2]

	mass = pysubs.Array3d(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz) # This is the grid file used to grid the Gadget input.
	TotalEnergy = pysubs.Array3d(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz) 
	PixelVolume[i] = mass.dx * mass.dy * mass.dz # This is the pixel volume in Gadget coordinates

	[mass,Density,InternalEnergy,Vx,Vy,Vz]  = pysubs.GridGadgetGasInit(snapfile,mass) # This grids the Gadget SPH particles
	masssum = mass.data.sum() * GadgetMass / MSun # This is just for the users benefit
	print "Total Gas Mass for grid %d is %e MSun\n"%(i,masssum)
	Tmin = Tbg * 1000.0 # Temp of background gas in eV
	Umin = Tmin / BulletConstants.Tconv
	#print 'Max internal energy = %.2g, Min Internal Energy = %.2g, Umin = %.2g \n'%(InternalEnergy.data.max(),InternalEnergy.data.min(),Umin)
	Rhomin = 10**LogRhoMin

	# Read in the magnetic field
	bx = pysubs.Array3d(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz)
	by = pysubs.Array3d(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz)
	bz = pysubs.Array3d(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz)	
	phi = pysubs.Array3d(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz)
	if i==0:
		bx1 = pysubs.Array3d(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz)
		by1 = pysubs.Array3d(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz)
		bz1 = pysubs.Array3d(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz)	
		for ii in range(nx):
			for jj in range(ny):
				for kk in range(nz):
					if Density.data[ii,jj,kk] < Rhomin:
						Density.data[ii,jj,kk] = Rhomin
						InternalEnergy.data[ii,jj,kk] = Umin
		if FieldMode == 'Random':
			infile = open('../../../../magField.dat','r')
			lines = infile.readlines()
			infile.close()
			for nn in range(len(lines)):
				kk = (nn % (nz*ny)) % nz
				jj = (nn % (nz*ny) - kk) / nz
				ii = (nn - jj*nz - kk) / (nz * ny)
				bx1.data[ii,jj,kk] = float(lines[nn].split()[0]) * MagScaleFactor * 1.0E-6 / EnzoBField
				by1.data[ii,jj,kk] = float(lines[nn].split()[1]) * MagScaleFactor * 1.0E-6 / EnzoBField
				bz1.data[ii,jj,kk] = float(lines[nn].split()[2]) * MagScaleFactor * 1.0E-6 / EnzoBField
		elif FieldMode == 'RandomRadial':
			infile = open('../../../../magField.dat','r')
			lines = infile.readlines()
			infile.close()
			SmoothedDensity=gaussian_filter(Density.data,3.0)
			RhoMax = SmoothedDensity.max()
			for nn in range(len(lines)):
				kk = (nn % (nz*ny)) % nz
				jj = (nn % (nz*ny) - kk) / nz
				ii = (nn - jj*nz - kk) / (nz * ny)
				bx1.data[ii,jj,kk] = float(lines[nn].split()[0]) * MagScaleFactor * SmoothedDensity[ii,jj,kk] / RhoMax * 1.0E-6 / EnzoBField
				by1.data[ii,jj,kk] = float(lines[nn].split()[1]) * MagScaleFactor * SmoothedDensity[ii,jj,kk] / RhoMax * 1.0E-6 / EnzoBField
				bz1.data[ii,jj,kk] = float(lines[nn].split()[2]) * MagScaleFactor * SmoothedDensity[ii,jj,kk] / RhoMax * 1.0E-6 / EnzoBField
		elif FieldMode == 'Constant':
			for ii in range(nx):
				for jj in range(ny):
					for kk in range(nz):
						bx1.data[ii,jj,kk] = 0.0
						by1.data[ii,jj,kk] = MagScaleFactor * 1.0E-6 / EnzoBField
						bz1.data[ii,jj,kk] = 0.0
		bx.data = bx1.data
		by.data = by1.data
		bz.data = bz1.data
	else:
		for ii in range(nx):
			for jj in range(ny):
				for kk in range(nz):
					if Density.data[ii,jj,kk] < Rhomin:
						Density.data[ii,jj,kk] = Rhomin
						InternalEnergy.data[ii,jj,kk] = Umin
					bx.data[ii,jj,kk] = pysubs.DataInterpolate3d(bx1,bx.x[ii],bx.y[jj],bx.z[kk])
					by.data[ii,jj,kk] = pysubs.DataInterpolate3d(by1,by.x[ii],by.y[jj],by.z[kk])
					bz.data[ii,jj,kk] = pysubs.DataInterpolate3d(bz1,bz.x[ii],bz.y[jj],bz.z[kk])
	for ii in range(nx):
		for jj in range(ny):
			for kk in range(nz):
				TotalEnergy.data[ii,jj,kk] = InternalEnergy.data[ii,jj,kk] + (Vx.data[ii,jj,kk]**2.0 + Vy.data[ii,jj,kk]**2.0 + Vz.data[ii,jj,kk]**2.0)/2.0 + (bx.data[ii,jj,kk]**2.0 + by.data[ii,jj,kk]**2.0 + bz.data[ii,jj,kk]**2.0) / (2.0 * Density.data[ii,jj,kk]) 


	print 'Max internal energy = %.2g, Min Internal Energy = %.2g, Umin = %.2g \n'%(InternalEnergy.data.max(),InternalEnergy.data.min(),Umin)
	print 'Max total energy = %.2g, Min total Energy = %.2g \n'%(TotalEnergy.data.max(),TotalEnergy.data.min())


	mass.data                = MassConversion * mass.data.reshape((1,)+mass.data.shape)
	Density.data             = DensityConversion * Density.data.reshape((1,)+Density.data.shape) 
	TotalEnergy.data         = (VelocityConversion**2) * TotalEnergy.data.reshape((1,)+TotalEnergy.data.shape)
	InternalEnergy.data      = (VelocityConversion**2) * InternalEnergy.data.reshape((1,)+InternalEnergy.data.shape)
	Vx.data                  = VelocityConversion * Vx.data.reshape((1,)+Vx.data.shape)
	Vy.data                  = VelocityConversion * Vy.data.reshape((1,)+Vy.data.shape)
	Vz.data                  = VelocityConversion * Vz.data.reshape((1,)+Vz.data.shape)
	bx.data = bx.data.reshape((1,)+bx.data.shape)
	by.data = by.data.reshape((1,)+by.data.shape)
	bz.data = bz.data.reshape((1,)+bz.data.shape)
	phi.data = phi.data.reshape((1,)+phi.data.shape)	
					
	if NumGrids==1:
		pysubs.WriteEnzoInitFile(Density.data, enzparam.CosmologySimulationDensityName,enzparam.CosmologySimulationDensityName)
		pysubs.WriteEnzoInitFile(Vx.data, enzparam.CosmologySimulationVelocity1Name, enzparam.CosmologySimulationVelocity1Name)
		pysubs.WriteEnzoInitFile(Vy.data, enzparam.CosmologySimulationVelocity2Name, enzparam.CosmologySimulationVelocity2Name)
		pysubs.WriteEnzoInitFile(Vz.data, enzparam.CosmologySimulationVelocity3Name, enzparam.CosmologySimulationVelocity3Name)
		pysubs.WriteEnzoInitFile(InternalEnergy.data, enzparam.CosmologySimulationGasEnergyName, enzparam.CosmologySimulationGasEnergyName)
		pysubs.WriteEnzoInitFile(TotalEnergy.data, enzparam.CosmologySimulationTotalEnergyName, enzparam.CosmologySimulationTotalEnergyName)
		pysubs.WriteEnzoInitFile(bx.data, enzparam.CosmologySimulationBfield1Name, enzparam.CosmologySimulationBfield1Name)
		pysubs.WriteEnzoInitFile(by.data, enzparam.CosmologySimulationBfield2Name, enzparam.CosmologySimulationBfield2Name)
		pysubs.WriteEnzoInitFile(bz.data, enzparam.CosmologySimulationBfield3Name, enzparam.CosmologySimulationBfield3Name)
		pysubs.WriteEnzoInitFile(phi.data, enzparam.CosmologySimulationPhiFieldName, enzparam.CosmologySimulationPhiFieldName)
	else:
		pysubs.WriteEnzoInitFile(Density.data, enzparam.CosmologySimulationDensityName+"."+str(i),enzparam.CosmologySimulationDensityName+"."+str(i))
		pysubs.WriteEnzoInitFile(Vx.data, enzparam.CosmologySimulationVelocity1Name+"."+str(i), enzparam.CosmologySimulationVelocity1Name+"."+str(i))
		pysubs.WriteEnzoInitFile(Vy.data, enzparam.CosmologySimulationVelocity2Name+"."+str(i), enzparam.CosmologySimulationVelocity2Name+"."+str(i))
		pysubs.WriteEnzoInitFile(Vz.data, enzparam.CosmologySimulationVelocity3Name+"."+str(i), enzparam.CosmologySimulationVelocity3Name+"."+str(i))
		pysubs.WriteEnzoInitFile(InternalEnergy.data, enzparam.CosmologySimulationGasEnergyName+"."+str(i), enzparam.CosmologySimulationGasEnergyName+"."+str(i))
		pysubs.WriteEnzoInitFile(TotalEnergy.data, enzparam.CosmologySimulationTotalEnergyName+"."+str(i), enzparam.CosmologySimulationTotalEnergyName+"."+str(i))
		pysubs.WriteEnzoInitFile(bx.data, enzparam.CosmologySimulationBfield1Name+"."+str(i), enzparam.CosmologySimulationBfield1Name+"."+str(i))
		pysubs.WriteEnzoInitFile(by.data, enzparam.CosmologySimulationBfield2Name+"."+str(i), enzparam.CosmologySimulationBfield2Name+"."+str(i))
		pysubs.WriteEnzoInitFile(bz.data, enzparam.CosmologySimulationBfield3Name+"."+str(i), enzparam.CosmologySimulationBfield3Name+"."+str(i))
		pysubs.WriteEnzoInitFile(phi.data, enzparam.CosmologySimulationPhiFieldName+"."+str(i), enzparam.CosmologySimulationPhiFieldName+"."+str(i))

# Now the DM particles - first we need to identify which grid they fall into
DMPos=zeros([3,Ndm,NumGrids])
DMVel=zeros([3,Ndm,NumGrids])
DMMass=zeros([1,Ndm,NumGrids])
Pos=zeros([3])
Vel=zeros([3])
NumParticles=zeros([NumGrids],'int')
print 'GadgetBoxCenter=',GadgetBoxCenter,' EnzoBoxCenter=',EnzoBoxCenter
for n in range(Ndm):
	no=Ngas+n # offset for gas particles
	for m in range(3):
	    Pos[m] = LengthConversion * (snapfile.data.pos[3*no+m] - GadgetBoxCenter[m]) + EnzoBoxCenter[m]
	    Vel[m] = VelocityConversion * snapfile.data.vel[3*no+m]
	if snapfile.data.masses[no] > 10.0:
		print 'ID = %d, OldX = %.3f, OldY = %.3f, OldZ = %.3f,\n'%(snapfile.data.id[no], snapfile.data.pos[3*no+0],snapfile.data.pos[3*no+1],snapfile.data.pos[3*no+2])
		print 'NewX = %.3f, NewY = %.3f, NewZ = %.3f,\n'%(Pos[0],Pos[1],Pos[2])
	for i in range(NumGrids-1, -1, -1): # Step from the smallest grid outwards checking if the particle is in this grid
		InThisGrid = True
		for m in range(3):
			if Pos[m] < enzparam.CosmologySimulationGridLeftEdge[i][m] or Pos[m] > enzparam.CosmologySimulationGridRightEdge[i][m]:
				InThisGrid=False
				break # If it's outside one coordinate, no need to test the others
		if InThisGrid:
			for m in range(3):
				DMPos[m,NumParticles[i],i] = Pos[m]
				DMVel[m,NumParticles[i],i] = Vel[m]
			if snapfile.header.Massarr[1]==0:
				DMMass[0,NumParticles[i],i] = DensityConversion * snapfile.data.masses[no] / PixelVolume[i]
			else:
				DMMass[0,NumParticles[i],i] = DensityConversion * snapfile.header.Massarr[1] / PixelVolume[i]
			NumParticles[i]=NumParticles[i]+1
			break # Leave the loop once you've found the smallest grid that contains the particle

# Now we've placed the particles in the appropriate grids, and can write the files
totsum=0.0
totpart=0
for i in range(NumGrids): 
	NPart=NumParticles[i]
	masssum=DMMass[0,0:Ndm,i].sum() * PixelVolume[i]  / DensityConversion * GadgetMass / MSun # Again, just for the users benefit
	totsum=totsum+masssum
	totpart=totpart+NPart
	print "DM Mass for %d particles in grid %d is %e MSun\n"%(NPart,i,masssum)
	PosData=zeros([3,NPart])
	VelData=zeros([3,NPart])
	MassData=zeros([1,NPart])
	for n in range(NPart):
		MassData[0,n]=DMMass[0,n,i] # Need to place the data in a contiguous array for the hdf writing to work
		for m in range(3):
			PosData[m,n]=DMPos[m,n,i] # Need to place the data in a contiguous array for the hdf writing to work
			VelData[m,n]=DMVel[m,n,i] # Need to place the data in a contiguous array for the hdf writing to work
	#if NPart != 0:
	if NumGrids==1:
		pysubs.WriteEnzoInitFile(PosData, enzparam.CosmologySimulationParticlePositionName, enzparam.CosmologySimulationParticlePositionName)
		pysubs.WriteEnzoInitFile(VelData, enzparam.CosmologySimulationParticleVelocityName, enzparam.CosmologySimulationParticleVelocityName)
		pysubs.WriteEnzoInitFile(MassData, enzparam.CosmologySimulationParticleMassName, enzparam.CosmologySimulationParticleMassName)
	else:
		pysubs.WriteEnzoInitFile(PosData, enzparam.CosmologySimulationParticlePositionName+"."+str(i), enzparam.CosmologySimulationParticlePositionName+"."+str(i))
		pysubs.WriteEnzoInitFile(VelData, enzparam.CosmologySimulationParticleVelocityName+"."+str(i), enzparam.CosmologySimulationParticleVelocityName+"."+str(i))
		pysubs.WriteEnzoInitFile(MassData, enzparam.CosmologySimulationParticleMassName+"."+str(i), enzparam.CosmologySimulationParticleMassName+"."+str(i))

print "Total DM Mass for %d particles is %e MSun\n"%(totpart,totsum)

#************END MAIN PROGRAM*************************

