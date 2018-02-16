#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 1-Jun-11


#This program parses the snapshot files from Gadget2 and plots out a 3d plot

import lageconfig # system specific path information
import sys
#import matplotlib
#matplotlib.use("PDF")
sys.path.append(lageconfig.bulletpath)
import  pysubs_mhd_apec_mass_only_1May12 as pysubs
from scipy.ndimage import gaussian_filter, convolve
from subprocess import *
from pylab import *
import BulletConstants
import time
#****************MAIN PROGRAM*****************

filenamebase = 'test2/snapshots/test2'
xlimit = 2000.0
zlimit = 100.0
xmin = ymin = -xlimit
xmax = ymax = xlimit
nx = ny = 128 
zmin = -zlimit
zmax = zlimit
figure()
subplots_adjust(hspace=0.4, wspace=0.4)
counter = 0
levels = [1,3.16,10,31.6,100.0,316.0]#[10**-4.0,10**-3.5,10**-3.0,10**-2.5,-2.0]
for snap in [0,25,50]:
	if snap<10:
	    filename = filenamebase+"_00"+str(snap) 
	elif snap<100:
	    filename = filenamebase+"_0"+str(snap) 
	else:
	    filename = filenamebase+"_"+str(snap) 

	snapfile = pysubs.ReadGadgetSnapshot(filename)
	#print snapfile.data.masses[3],snapfile.header.Npart[1],data.data.shape, data.data.sum(),data.data.max()
	simtime=snapfile.header.Time[0]*BulletConstants.TimeConversion # Put time in Gyears


	[NumPart, Masses, Centroids] = pysubs.FindGadgetCentroids(snapfile)

	print "Snap = %d, Centroids = %f, %f, %f"%(snap,Centroids[1,1,0],Centroids[1,1,1],Centroids[1,1,2])

	xmin = -xlimit + Centroids[1,1,1]
	xmax = xlimit + Centroids[1,1,1]
	ymin = -xlimit + Centroids[1,1,2]
	ymax = xlimit + Centroids[1,1,2]
	zmin = -zlimit + Centroids[1,1,0]
	zmax = zlimit + Centroids[1,1,0]
	data = pysubs.Array2d(xmin,xmax,nx,ymin,ymax,ny)
	[xx,yy]=meshgrid(data.x,data.y)


	data = pysubs.ProjectGadgetDM(snapfile,data,1,zmin,zmax,theta = pi/2.0, phi = pi/2.0)
	data.data=gaussian_filter(data.data,3.0)
	subplot(3,3,1+3*counter, aspect = 1.0)
	title('X Density, T=%.2f Gy'%simtime)	
	contour(xx,yy,data.data, levels,colors='b', linewidths=2.0)
	xlim(xmin,xmax)
	ylim(ymin,ymax)
	gca().set_xticks([-xlimit/2,0.0,xlimit/2])
	gca().set_yticks([-xlimit/2,0.0,xlimit/2])

	xmin = -xlimit + Centroids[1,1,0]
	xmax = xlimit + Centroids[1,1,0]
	ymin = -xlimit + Centroids[1,1,2]
	ymax = xlimit + Centroids[1,1,2]
	zmin = -zlimit - Centroids[1,1,1]
	zmax = zlimit - Centroids[1,1,1]
	data = pysubs.Array2d(xmin,xmax,nx,ymin,ymax,ny)
	[xx,yy]=meshgrid(data.x,data.y)

	data = pysubs.ProjectGadgetDM(snapfile,data,1,zmin,zmax,theta = pi/2.0, phi = 0.0)
	data.data=gaussian_filter(data.data,3.0)
	subplot(3,3,2+3*counter, aspect = 1.0)
	title('Y Density, T=%.2f Gy'%simtime)	
	contour(xx,yy,data.data, levels,colors='g', linewidths=2.0)
	xlim(xmin,xmax)
	ylim(ymin,ymax)
	gca().set_xticks([-xlimit/2,0.0,xlimit/2])
	gca().set_yticks([-xlimit/2,0.0,xlimit/2])


	xmin = -xlimit + Centroids[1,1,0]
	xmax = xlimit + Centroids[1,1,0]
	ymin = -xlimit + Centroids[1,1,1]
	ymax = xlimit + Centroids[1,1,1]
	zmin = -zlimit + Centroids[1,1,2]
	zmax = zlimit + Centroids[1,1,2]
	data = pysubs.Array2d(xmin,xmax,nx,ymin,ymax,ny)
	[xx,yy]=meshgrid(data.x,data.y)

	data = pysubs.ProjectGadgetDM(snapfile,data,1,zmin,zmax,theta = 0.0, phi = 0.0)
	data.data=gaussian_filter(data.data,3.0)
	subplot(3,3,3+3*counter, aspect = 1.0)
	title('Z Density, T=%.2f Gy'%simtime)	
	contour(xx,yy,data.data, levels,colors='r', linewidths=2.0)
	xlim(xmin,xmax)
	ylim(ymin,ymax)
	gca().set_xticks([-xlimit/2,0.0,xlimit/2])
	gca().set_yticks([-xlimit/2,0.0,xlimit/2])





	counter = counter + 1
show()



#****************END MAIN PROGRAM*****************
