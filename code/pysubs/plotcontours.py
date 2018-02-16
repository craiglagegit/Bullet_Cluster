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
import  pysubs_nas_11Feb13 as pysubs
from scipy.ndimage import gaussian_filter, convolve
from subprocess import *
from pylab import *
import BulletConstants
import time
#****************MAIN PROGRAM*****************

filenamebase = 'plots/'
pp=PdfPages('plots/Graph_Contours.pdf')
xlimit = 800.0
zlimit = 100.0
xmin = ymin = -xlimit
xmax = ymax = xlimit
nx = ny = 128 
zmin = -zlimit
zmax = zlimit
counter = 0
levels = [1.0,1.78,3.16,5.62,10.0,17.8,31.6,56.2,100.0,178.0,316.0,562.0,1000.0]#[10**-4.0,10**-3.5,10**-3.0,10**-2.5,-2.0]
figure()
figure()
subplots_adjust(hspace=0.8, wspace=0.4)
for snap in [10,51,99]:
	pf = pysubs.GetPF(snap)
	simtime = pf.h.parameters['InitialTime'] * BulletConstants.TimeConversion # Put time in Gyears
	ApecData = pysubs.ReadLookups(1.0) # Reads the APEC lookup tables.
	for grid in pf.h.grids:
		grid.set_field_parameter('ApecData',ApecData)
                grid.set_field_parameter('TFudge',0.5)

	#[NumPart, Masses, Centroids] = pysubs.FindEnzoCentroids(pf)
	Centroids=zeros([2,2,3])
	print "Snap = %d, Centroids = %f, %f, %f"%(snap,Centroids[1,1,0],Centroids[1,1,1],Centroids[1,1,2])

	xmin = -xlimit + Centroids[1,1,1]
	xmax = xlimit + Centroids[1,1,1]
	ymin = -xlimit + Centroids[1,1,2]
	ymax = xlimit + Centroids[1,1,2]
	zmin = -zlimit + Centroids[1,1,0]
	zmax = zlimit + Centroids[1,1,0]
	mass = pysubs.Array2d(xmin,xmax,nx,ymin,ymax,ny)
	[xx,yy]=meshgrid(mass.x,mass.y)

	[DM, mass,Xray1,Xray2,Xray3,SZ] = pysubs.ProjectEnzoData(pf,mass,phi=0.0,theta=0.0,psi=0.0,zmin=zmin,zmax=zmax,DMProject=False)

	DM.data=gaussian_filter(DM.data,3.0)
	subplot(3,3,1+3*counter, aspect = 1.0)
	title('DM Density, T=%.2f Gy'%simtime,fontsize=10)	
	contour(xx,yy,DM.data, levels,colors='b', linewidths=2.0)
	xlim(xmin,xmax)
	ylim(ymin,ymax)
	gca().set_xticks([-xlimit/2,xlimit/2])
	gca().set_yticks([-xlimit/2,0.0,xlimit/2])
	xlabel('Y')
	ylabel('X',rotation='horizontal')
	xmin = -xlimit + Centroids[1,1,0]
	xmax = xlimit + Centroids[1,1,0]
	ymin = -xlimit + Centroids[1,1,2]
	ymax = xlimit + Centroids[1,1,2]
	zmin = -zlimit - Centroids[1,1,1]
	zmax = zlimit - Centroids[1,1,1]
	mass = pysubs.Array2d(xmin,xmax,nx,ymin,ymax,ny)
	[xx,yy]=meshgrid(mass.x,mass.y)

	[DM, mass,Xray1,Xray2,Xray3,SZ] = pysubs.ProjectEnzoData(pf,mass,phi=0.0,theta=pi/2.0,psi=0.0,zmin=zmin,zmax=zmax,DMProject=False)

	DM.data=gaussian_filter(DM.data,3.0)
	subplot(3,3,2+3*counter, aspect = 1.0)
	title('DM Density, T=%.2f Gy'%simtime,fontsize=10)	
	contour(xx,yy,DM.data, levels,colors='g', linewidths=2.0)
	xlim(xmin,xmax)
	ylim(ymin,ymax)
	gca().set_xticks([-xlimit/2,xlimit/2])
	gca().set_yticks([-xlimit/2,0.0,xlimit/2])
	xlabel('Z')
	ylabel('X',rotation='horizontal')


	xmin = -xlimit + Centroids[1,1,0]
	xmax = xlimit + Centroids[1,1,0]
	ymin = -xlimit + Centroids[1,1,1]
	ymax = xlimit + Centroids[1,1,1]
	zmin = -zlimit + Centroids[1,1,2]
	zmax = zlimit + Centroids[1,1,2]
	mass = pysubs.Array2d(xmin,xmax,nx,ymin,ymax,ny)
	[xx,yy]=meshgrid(mass.x,mass.y)

	[DM, mass,Xray1,Xray2,Xray3,SZ] = pysubs.ProjectEnzoData(pf,mass,phi=0.0,theta=pi/2.0,psi=pi/2.0,zmin=zmin,zmax=zmax,DMProject=False)

	DM.data=gaussian_filter(DM.data,3.0)
	subplot(3,3,3+3*counter, aspect = 1.0)
	title('DM Density, T=%.2f Gy'%simtime,fontsize=10)	
	contour(xx,yy,DM.data, levels,colors='r', linewidths=2.0)
	xlim(xmin,xmax)
	ylim(ymin,ymax)
	gca().set_xticks([-xlimit/2,xlimit/2])
	gca().set_yticks([-xlimit/2,0.0,xlimit/2])
	xlabel('Z')
	ylabel('Y',rotation='horizontal')

	


	counter = counter + 1
pp.savefig()
plt.clf()
pp.close()
#show()



#****************END MAIN PROGRAM*****************
