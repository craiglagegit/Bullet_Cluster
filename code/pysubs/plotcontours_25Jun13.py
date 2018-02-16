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
import  pysubs_nas_20May13 as pysubs
from scipy.ndimage import gaussian_filter, convolve
from subprocess import *
from pylab import *
import BulletConstants
import time
#****************MAIN PROGRAM*****************

filenamebase = 'plots/'
pp=PdfPages('plots/Graph_Contours_New.pdf')
xlimit = 800.0
zlimit = 100.0
plotlimit = 1000.0
xmin = ymin = -xlimit
xmax = ymax = xlimit
nx = ny = 128 
zmin = -zlimit
zmax = zlimit
counter = 0
levels = [1.0,1.78,3.16,5.62,10.0,17.8,31.6,56.2,100.0,178.0,316.0,562.0,1000.0]#[10**-4.0,10**-3.5,10**-3.0,10**-2.5,-2.0]
fig = figure()
Ax = []
snaps = [10,51,99]
for i in range(3):
	if i == 2: lp = -25
	else: lp = -10
	yplot = 0.65 - 0.30 * i
	pf = pysubs.GetPF(snaps[i])
	simtime = pf.h.parameters['InitialTime'] * BulletConstants.TimeConversion # Put time in Gyears
	ApecData = pysubs.ReadLookups(1.0) # Reads the APEC lookup tables.
	for grid in pf.h.grids:
		grid.set_field_parameter('ApecData',ApecData)
                grid.set_field_parameter('TFudge',0.5)

	mass = pysubs.Array2d(xmin,xmax,nx,ymin,ymax,ny)
	[xx,yy]=meshgrid(mass.x,mass.y)
	[DM, mass,Xray1,Xray2,Xray3,SZ] = pysubs.ProjectEnzoData(pf,mass,phi=0.0,theta=0.0,psi=0.0,zmin=zmin,zmax=zmax,DMProject=False)
	DM.data=gaussian_filter(DM.data,3.0)
	xplot = 0.17
	Ax.append(axes([xplot,yplot,0.30,0.30],aspect = 1.0))
	Ax[-1].set_title('DM Density, T=%.2f Gy'%simtime,fontsize=10,x=0.5,y=0.88)	
	contour(xx,yy,DM.data, levels,colors='b', linewidths=2.0)
	Ax[-1].set_xlim(-plotlimit,plotlimit)
	Ax[-1].set_ylim(-plotlimit,plotlimit)
	if i ==2:
		Ax[-1].set_xticks([-xlimit/2,xlimit/2])
	else:
		Ax[-1].set_xticks([])
	Ax[-1].set_xlabel('Y',labelpad=lp)
	Ax[-1].set_ylabel('X',rotation='horizontal',labelpad=-50)
	Ax[-1].set_yticks([-xlimit/2,0.0,xlimit/2])
	Ax[-1].yaxis.tick_left()

	[DM, mass,Xray1,Xray2,Xray3,SZ] = pysubs.ProjectEnzoData(pf,mass,phi=0.0,theta=pi/2.0,psi=0.0,zmin=zmin,zmax=zmax,DMProject=False)
	DM.data=gaussian_filter(DM.data,3.0)
	xplot = 0.39
	Ax.append(axes([xplot,yplot,0.30,0.30],aspect = 1.0))
	Ax[-1].set_title('DM Density, T=%.2f Gy'%simtime,fontsize=10,x=0.5,y=0.88)	
	contour(xx,yy,DM.data, levels,colors='g', linewidths=2.0)
	Ax[-1].set_xlim(-plotlimit,plotlimit)
	Ax[-1].set_ylim(-plotlimit,plotlimit)
	if i ==2:
		Ax[-1].set_xticks([-xlimit/2,xlimit/2])
	else:
		Ax[-1].set_xticks([])
	Ax[-1].set_yticks([])
	Ax[-1].set_xlabel('Z',labelpad=lp)
	Ax[-1].set_ylabel('X',rotation='horizontal',labelpad=-10)

	[DM, mass,Xray1,Xray2,Xray3,SZ] = pysubs.ProjectEnzoData(pf,mass,phi=0.0,theta=pi/2.0,psi=pi/2.0,zmin=zmin,zmax=zmax,DMProject=False)
	DM.data=gaussian_filter(DM.data,3.0)
	xplot = 0.61
	Ax.append(axes([xplot,yplot,0.30,0.30],aspect = 1.0))
	Ax[-1].set_title('DM Density, T=%.2f Gy'%simtime,fontsize=10,x=0.5,y=0.88)	
	contour(xx,yy,DM.data, levels,colors='r', linewidths=2.0)
	Ax[-1].set_xlim(-plotlimit,plotlimit)
	Ax[-1].set_ylim(-plotlimit,plotlimit)
	if i ==2:
		Ax[-1].set_xticks([-xlimit/2,xlimit/2])
	else:
		Ax[-1].set_xticks([])
	Ax[-1].set_yticks([])
	Ax[-1].set_xlabel('Z',x=0.5,y=0.1,labelpad=lp)
	Ax[-1].set_ylabel('Y',rotation='horizontal',labelpad=-10)
pp.savefig()
plt.clf()
pp.close()
#show()



#****************END MAIN PROGRAM*****************
