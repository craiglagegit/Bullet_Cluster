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

filenamebase = 'test2/snapshots/fulltest3'
xmin = ymin = zmin = -2000.0
xmax = ymax = zmax = 2000.0
nx = ny = nz = 128 
plotscale = 1000.0
mass = pysubs.Array3d(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz)
figure()
for snap in  [0,9,18]:#[0,10,20,30,40,50]:
	if snap<10:
	    filename = filenamebase+"_00"+str(snap) 
	elif snap<100:
	    filename = filenamebase+"_0"+str(snap) 
	else:
	    filename = filenamebase+"_"+str(snap) 

	snapfile = pysubs.ReadGadgetSnapshot(filename)

	[NumPart, Masses, Centroids] = pysubs.FindGadgetCentroids(snapfile)

	print "Snap = %d, Centroids = %f, %f, %f"%(snap,Centroids[1,1,0],Centroids[1,1,1],Centroids[1,1,2])

	xmin = xmin + Centroids[1,1,0]
	xmax = xmax + Centroids[1,1,0]
	ymin = ymin + Centroids[1,1,1]
	ymax = ymax + Centroids[1,1,1]
	zmin = zmin + Centroids[1,1,2]
	zmax = zmax + Centroids[1,1,2]
	mass = pysubs.Array3d(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz)
	rho = pysubs.Array3d(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz)


	[mass,rho] = pysubs.GridGadgetDM(snapfile,mass,1)
	print "Total DM mass for snap %d = %f\n"%(snap,mass.data.sum())
	simtime=snapfile.header.Time[0]*BulletConstants.TimeConversion # Put time in Gyears

	subplot(3,3,1)
	title('X-axis DM Density profile')	
	plot(mass.x-Centroids[1,1,0], log10((rho.data[:,ny/2,nz/2]+rho.data[:,ny/2-1,nz/2]+rho.data[:,ny/2,nz/2-1]+rho.data[:,ny/2-1,nz/2-1])/4.0), label = 'Time = %.2f Gy'%simtime)
	xlim(-plotscale,plotscale)
	ylim(-6.0,-2.0)

	subplot(3,3,2)
	title('Y-axis DM Density profile')
	plot(mass.y-Centroids[1,1,1], log10((rho.data[nx/2,:,nz/2]+rho.data[nx/2-1,:,nz/2]+rho.data[nx/2,:,nz/2-1]+rho.data[nx/2-1,:,nz/2-1])/4.0), label = 'Time = %.2f Gy'%simtime)
	xlim(-plotscale,plotscale)
	ylim(-6.0,-2.0)

	subplot(3,3,3)
	title('Z-axis DM Density profile')
	plot(mass.z-Centroids[1,1,2], log10((rho.data[nx/2,ny/2,:]+rho.data[nx/2-1,ny/2,:]+rho.data[nx/2,ny/2-1,:]+rho.data[nx/2-1,ny/2-1,:])/4.0), label = 'Time = %.2f Gy'%simtime)
	xlim(-plotscale,plotscale)
	ylim(-6.0,-2.0)

	[mass,temp,rho,xray,pressure,SZ] = pysubs.GridGadgetGas(snapfile,mass)
	print "Total Gas mass for snap %d = %f\n"%(snap,mass.data.sum())

	subplot(3,3,4)
	title('X-axis Gas Density profile')
	plot(mass.x-Centroids[1,1,0], log10((rho.data[:,ny/2,nz/2]+rho.data[:,ny/2-1,nz/2]+rho.data[:,ny/2,nz/2-1]+rho.data[:,ny/2-1,nz/2-1])/4.0), label = 'Time = %.2f Gy'%simtime)
	xlim(-plotscale,plotscale)
	ylim(-7.0,-3.0)

	subplot(3,3,5)
	title('Y-axis Gas Density profile')
	plot(mass.y-Centroids[1,1,1], log10((rho.data[nx/2,:,nz/2]+rho.data[nx/2-1,:,nz/2]+rho.data[nx/2,:,nz/2-1]+rho.data[nx/2-1,:,nz/2-1])/4.0), label = 'Time = %.2f Gy'%simtime)
	xlim(-plotscale,plotscale)
	ylim(-7.0,-3.0)

	subplot(3,3,6)
	title('Z-axis Gas Density profile')
	plot(mass.z-Centroids[1,1,2], log10((rho.data[nx/2,ny/2,:]+rho.data[nx/2-1,ny/2,:]+rho.data[nx/2,ny/2-1,:]+rho.data[nx/2-1,ny/2-1,:])/4.0), label = 'Time = %.2f Gy'%simtime)
	xlim(-plotscale,plotscale)
	ylim(-7.0,-3.0)

	subplot(3,3,7)
	title('X-axis Gas Temp profile')
	plot(mass.x-Centroids[1,1,0], log10((temp.data[:,ny/2,nz/2]+temp.data[:,ny/2-1,nz/2]+temp.data[:,ny/2,nz/2-1]+temp.data[:,ny/2-1,nz/2-1])/4.0), label = 'Time = %.2f Gy'%simtime)
	xlim(-plotscale,plotscale)
	ylim(0.0,5.0)

	subplot(3,3,8)
	title('Y-axis Gas Temp profile')
	plot(mass.y-Centroids[1,1,1], log10((temp.data[nx/2,:,nz/2]+temp.data[nx/2-1,:,nz/2]+temp.data[nx/2,:,nz/2-1]+temp.data[nx/2-1,:,nz/2-1])/4.0), label = 'Time = %.2f Gy'%simtime)
	xlim(-plotscale,plotscale)
	ylim(0.0,5.0)

	subplot(3,3,9)
	title('Z-axis Gas Temp profile')
	plot(mass.z-Centroids[1,1,2], log10((temp.data[nx/2,ny/2,:]+temp.data[nx/2-1,ny/2,:]+temp.data[nx/2,ny/2-1,:]+temp.data[nx/2-1,ny/2-1,:])/4.0), label = 'Time = %.2f Gy'%simtime)
	xlim(-plotscale,plotscale)
	ylim(0.0,5.0)



legend(loc=1,bbox_to_anchor=[1.25, 0.25])
show()



#****************END MAIN PROGRAM*****************
