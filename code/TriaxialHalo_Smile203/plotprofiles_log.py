#!/usr/bin/env python

#Author: Craig Lage1 NYU; 
#Date: 1-Jun-11


#This program parses the snapshot files from Gadget0 and plots out a 3d plot

import lageconfig # system specific path information
import sys
#import matplotlib
#matplotlib.use("PDF")
sys.path.append(lageconfig.bulletpath)
import  pysubs_mhd_apec_mass_only_1May12 as pysubs
from subprocess import *
from pylab import *
import BulletConstants
import time





#****************MAIN PROGRAM*****************

filenamebase = 'test2/snapshots/test2'
xmin = ymin = zmin = 0.0
xmax = ymax = zmax = 2000.0
nx = ny = nz =128 
mass = pysubs.Array3d(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz)
figure()
for snap in  [0,4,9]:#[0,10,20,30,40,50]:
	if snap<10:
	    filename = filenamebase+"_00"+str(snap) 
	elif snap<100:
	    filename = filenamebase+"_0"+str(snap) 
	else:
	    filename = filenamebase+"_"+str(snap) 

	snapfile = pysubs.ReadGadgetSnapshot(filename)

	[NumPart, Masses, Centroids] = pysubs.FindGadgetCentroids(snapfile)

	print "Snap = %d, Centroids = %f, %f, %f"%(snap,Centroids[1,1,0],Centroids[1,1,1],Centroids[1,1,2])

	xmin = Centroids[1,1,0]
	xmax = 2000.0 + Centroids[1,1,0]
	ymin = Centroids[1,1,1]
	ymax = 2000.0 + Centroids[1,1,1]
	zmin = Centroids[1,1,2]
	zmax = 2000.0 + Centroids[1,1,2]
	mass = pysubs.Array3d(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz)
	rho = pysubs.Array3d(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz)


	[mass,rho] = pysubs.GridGadgetDM(snapfile,mass,1)
	simtime=snapfile.header.Time[0]*BulletConstants.TimeConversion # Put time in Gyears
	subplot(1,3,1)
	title('X-axis Density profile')	
	plot(log10(abs(mass.x-Centroids[1,1,0])), log10(rho.data[:,0,0]), label = 'Time = %.2f Gy'%simtime)
	xlim(0.0,4.0)
	ylim(-5.0,-1.0)
	subplot(1,3,2)
	title('Y-axis Density profile')
	plot(log10(abs(mass.y-Centroids[1,1,1])), log10(rho.data[0,:,0]), label = 'Time = %.2f Gy'%simtime)
	xlim(0.0,4.0)
	ylim(-5.0,-1.0)
	subplot(1,3,3)
	title('Z-axis Density profile')
	plot(log10(abs(mass.z-Centroids[1,1,2])), log10(rho.data[0,0,:]), label = 'Time = %.2f Gy'%simtime)
	xlim(0.0,4.0)
	ylim(-5.0,-1.0)

subplot(1,3,1)
plot([1.0,2.0], [-1.5,-2.5], label = 'r^-1')
plot([2.0,2.5], [-2.5,-4.0], label = 'r^-3')


legend(loc=1,bbox_to_anchor=[1.00, 0.97])
show()



#****************END MAIN PROGRAM*****************
