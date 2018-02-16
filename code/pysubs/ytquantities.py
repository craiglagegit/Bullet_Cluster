#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 5-Mar-11


#This program makes a series of plots using yt.

import lageconfig # system specific path information
import sys
sys.path.append(lageconfig.bulletpath)
#import pysubs_nas_22Jan14 as pysubs
import BulletConstants
import time
from pylab import *
from yt.mods import *

def GetPF(snap, dir):
	# Retrieves a yt - pf file given a snapshot number
	if snap<10:
	    filename = dir+"DD000"+str(snap)+"/output_000"+str(snap) 
	elif snap<100:
	    filename = dir+"DD00"+str(snap)+"/output_00"+str(snap) 
	elif snap<1000:
	    filename = dir+"DD0"+str(snap)+"/output_0"+str(snap) 
	else:
	    filename = dir+"DD"+str(snap)+"/output_"+str(snap) 
	
	pf = load(filename)
	return pf

def _MyGasEnergy(field, data):
	return data["Density"] * data["Temperature"] * data["CellVolume"] * BulletConstants.kBoltzmann / BulletConstants.mstar * (1.0 + BulletConstants.NeNi) * 3.0 / 2.0

def _MyTotalEnergy(field, data):
	EV = data["Density"] * (data["x-velocity"]*data["x-velocity"]+data["y-velocity"]*data["y-velocity"]+data["z-velocity"]*data["z-velocity"])/2.0
	BV = (data["Bx"]*data["Bx"]+data["By"]*data["By"]+data["Bz"]*data["Bz"])/(8.0*pi)
	GV = data["Density"] * data["Temperature"] * BulletConstants.kBoltzmann / BulletConstants.mstar * (1.0 + BulletConstants.NeNi) * 3.0 / 2.0
	return data["CellVolume"] * (GV + EV + BV)

def _MyMass(field, data):
	return data["Density"] * data["CellVolume"]

def _TotalXRay(field, data):
	ApecData = data.get_field_parameter('ApecData')
	TFudge = data.get_field_parameter('TFudge')
	T = data['Temperature'] * TFudge
	logT = np.log10(T)*100
	np.putmask(logT,logT>899.0,899.0)
	np.putmask(logT,logT<400.0,400.0)
	minT = logT.astype(int)
	flux = (minT+1.0-logT) * ApecData[minT-399] + (logT-minT) * ApecData[minT-398]
	return 10**flux * data['Density']**2 / BulletConstants.mp**2 * data["CellVolume"]

def ReadLookups(Z):
	# This subroutine reads in the data from the APEC cooling tables and places it in an array
	# The data is interpolated to get the data for the required Z (metallicity)
	ApecData = zeros([500]) # Array to hold the data
	infile = open(lageconfig.toppath+'bullet/code/pysubs/apec_cooling.in','r')
	dump = infile.readline() # Skips header line
	dump = infile.readline() # Skips header line
	lines = infile.readlines()
	infile.close()
	counter = 0
	for line in lines:
		minZ = max(0, int(round((Z * 10))))
		maxZ = minZ + 1
		if maxZ > 10:
			maxZ = 11
			minZ = 10
		f = (maxZ-10.0*Z) * float(line.strip().split()[minZ+1]) + (10.0*Z-minZ) * float(line.strip().split()[maxZ+1])
		ApecData[counter] = f
		counter=counter+1
	return ApecData

#****************MAIN PROGRAM*****************

cmd =sys.argv

dir = 'ddfiles/run'+cmd[1]+'/'
snapmin = int(cmd[2])
snapmax = int(cmd[3])


for snap in range(snapmin,snapmax):

	pf = GetPF(snap, dir)
	simtime = pf.h.parameters['InitialTime'] * BulletConstants.TimeConversion # Put time in Gyears

	ApecData = ReadLookups(0.7) # Reads the APEC lookup tables.
	for grid in pf.h.grids:
		grid.set_field_parameter('ApecData',ApecData)
                grid.set_field_parameter('TFudge',1.0)
                grid.set_field_parameter('SpectralIndex',3.5)


		#add_field('XRay', function = pysubs._EnzoXRay, units=r"\rm{cm}^{-3}\rm{s}^{-1}", projected_units=r"\rm{cm}^{-2}\rm{s}^{-1}", validators=[ValidateParameter('ApecData')], take_log=False)

        add_field('TotalXRay', function = _TotalXRay, units=r"\rm{s}^{-1}", projected_units=r"\rm{s}^{-1}", validators=[ValidateParameter('ApecData')], take_log=False)

        add_field('MyGasEnergy', function = _MyGasEnergy, units=r"\rm{erg}{cm}^{-3}\rm{s}^{-1}", projected_units=r"\rm{erg}{cm}^{-2}\rm{s}^{-1}", take_log=False)

        add_field('MyTotalEnergy', function = _MyTotalEnergy, units=r"\rm{erg}{cm}^{-3}\rm{s}^{-1}", projected_units=r"\rm{erg}{cm}^{-2}\rm{s}^{-1}", take_log=False)
        add_field('MyMass', function = _MyMass, units=r"\rm{g}{cm}^{-3}\rm{s}^{-1}", projected_units=r"\rm{g}{cm}^{-2}\rm{s}^{-1}", take_log=False)
    
        sp = pf.h.sphere("max",(2.0, "mpc"));
	sp.set_field_parameter('ApecData',ApecData)
	sp.set_field_parameter('TFudge',1.0)
	sp.set_field_parameter('SpectralIndex',3.5)

        [ParticleMass, GasMass] = sp.quantities["TotalQuantity"](["ParticleMass","CellMass"])
        print "DMMass = %g, GasMass = %g\n"%(ParticleMass,GasMass)
        MyMass = sp.quantities["TotalQuantity"]("MyMass")
        print "MyMass = %g\n"%MyMass[0]

        [ParticleMassMsun, GasMassMsun] = sp.quantities["TotalQuantity"](["ParticleMassMsun","CellMassMsun"])
        print "DMMassMsun = %g, GasMassMsun = %g, GasFraction = %f\n"%(ParticleMassMsun,GasMassMsun,GasMassMsun/ParticleMassMsun)

        MyGasEnergy = sp.quantities["TotalQuantity"]("MyGasEnergy")[0]
        MyTotalEnergy = sp.quantities["TotalQuantity"]("MyTotalEnergy")[0]
        print "TotalEnergy = %g, GasEnergy = %g, GasFraction = %f\n"%(MyTotalEnergy,MyGasEnergy,MyGasEnergy/MyTotalEnergy)

        MyTotalXRay = sp.quantities["TotalQuantity"]("TotalXRay")[0]
        print "TotalXRay = %g\n"%MyTotalXRay

	print "Decay Time = %g Gy \n"%(MyGasEnergy/MyTotalXRay/BulletConstants.sec_per_Gy)


#************END MAIN PROGRAM*************************

