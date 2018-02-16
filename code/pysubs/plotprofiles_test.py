#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 1-Jun-11


#This program parses the snapshot files from Gadget2 and plots out a 3d plot

import lageconfig # system specific path information
import sys
import matplotlib
matplotlib.use("PDF")
sys.path.append(lageconfig.bulletpath)
from matplotlib.backends.backend_pdf import PdfPages
import  pysubs_nas_11Feb13_test as pysubs
from scipy.ndimage import gaussian_filter, convolve
from subprocess import *
from pylab import *
import BulletConstants
import time
#****************MAIN PROGRAM*****************

xmin = -1500.0
xmax = 1500.0
ymin = zmin = -1500.0
ymax = zmax = 1500.0
nx = ny = 128 
mass = pysubs.Array2d(xmin,xmax,nx,ymin,ymax,ny)
figure()
subplots_adjust(hspace=0.8, wspace=0.4)
for snap in [0,4,51,99]:
	pf = pysubs.GetPF(snap)
	simtime = pf.h.parameters['InitialTime'] * BulletConstants.TimeConversion # Put time in Gyears
	ApecData = pysubs.ReadLookups(1.0) # Reads the APEC lookup tables.
	for grid in pf.h.grids:
		grid.set_field_parameter('ApecData',ApecData)
                grid.set_field_parameter('TFudge',1.0)

	[DM, mass,Xray1,Xray2,Xray3,SZ] = pysubs.ProjectEnzoData(pf,mass,phi=0.0,theta=0.0,psi=0.0,zmin=zmin,zmax=zmax,DMProject=False)


	print "Snap = %d, Total DM mass = %f, Total Gas mass = %f"%(snap,DM.data.sum(),mass.data.sum())

#****************END MAIN PROGRAM*****************
