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
import  pysubs_nas_3Apr14 as pysubs
from scipy.ndimage import gaussian_filter, convolve
from subprocess import *
from pylab import *
import BulletConstants
import time
from yt.pmods import *
from mpi4py import MPI

#****************MAIN PROGRAM*****************                                                                                                               
comm = MPI.COMM_WORLD
rank = comm.Get_rank()


filenamebase = 'plots/'

if rank == 0:
	pp=PdfPages('plots/Graph_Profiles_2.pdf')

xmin = -2000.0
xmax = 2000.0
ymin = zmin = -20.0
ymax = zmax = 20.0
nx = 128 
ny = 3
plotscale = 1000.0
mass = pysubs.Array2d(xmin,xmax,nx,ymin,ymax,ny)
if rank == 0:
	figure()
	subplots_adjust(hspace=0.8, wspace=0.4)
for snap in [65]:#[0,10,39,59]:#[1,40,80]:#,97]:
	if snap == 65:#0:
		ls = '-'
	elif snap == 10:
		ls = '--'
	elif snap == 39:
		ls = '--'
	elif snap == 59:
		ls = '-'
	pf = pysubs.GetPF(snap)
	simtime = pf.h.parameters['InitialTime'] * BulletConstants.TimeConversion # Put time in Gyears
	ApecData = pysubs.ReadLookups(1.0) # Reads the APEC lookup tables.
	for grid in pf.h.grids:
		grid.set_field_parameter('ApecData',ApecData)
                grid.set_field_parameter('TFudge',1.0)
                grid.set_field_parameter('SpectralIndex',3.5)

	#[NumPart, Masses, Centroids] = pysubs.FindEnzoCentroids(pf)
	Centroids=zeros([2,2,3])

	#print "Snap = %d, Centroids = %f, %f, %f"%(snap,Centroids[1,1,0],Centroids[1,1,1],Centroids[1,1,2])

	xmin = xmin + Centroids[1,1,0]
	xmax = xmax + Centroids[1,1,0]
	ymin = ymin + Centroids[1,1,1]
	ymax = ymax + Centroids[1,1,1]
	zmin = zmin + Centroids[1,1,2]
	zmax = zmax + Centroids[1,1,2]

	[DM, mass,Xray1,Xray2,Xray3,SZ] = pysubs.ProjectEnzoData(pf,mass,phi=0.0,theta=0.0,psi=0.0,zmin=zmin,zmax=zmax,DMProject=False)
	[Temp] = pysubs.ProjectEnzoTemp(pf,mass,phi=0.0,theta=0.0,psi=0.0,zmin=zmin,zmax=zmax,MHD=False)	
	if rank == 0:
		#print "X-axis, TrueTemp = ",TrueTemp.data, "\nTemp = ",Temp.data
		subplot(4,3,1)
		title('X-axis Log DM Density',fontsize=10)	
		plot(mass.x-Centroids[1,1,0], (log10(DM.data[:,0])+log10(DM.data[:,1])+log10(DM.data[:,2]))/3.0, label = 'Time = %.2f Gy'%simtime, ls = ls)
		xlim(-plotscale,plotscale)
		ylim(-1.0,2.0)
		gca().set_xticks([-plotscale,0.0,plotscale])
		gca().set_yticks([-1.0,0.0,1.0,2.0])
		subplot(4,3,4)
		title('X-axis Log Gas Density',fontsize=10)	
		plot(mass.x-Centroids[1,1,0], (log10(mass.data[:,0])+log10(mass.data[:,1])+log10(mass.data[:,2]))/3.0, label = 'Time = %.2f Gy'%simtime, ls = ls)
		xlim(-plotscale,plotscale)
		ylim(-2.0,1.0)
		gca().set_xticks([-plotscale,0.0,plotscale])
		subplot(4,3,7)
		title('X-axis Gas Temp(keV)',fontsize=10)	
		plot(mass.x-Centroids[1,1,0], (Temp.data[:,0]+Temp.data[:,1]+Temp.data[:,2])/3.0 , label = 'Time = %.2f Gy'%simtime, ls = ls)
		xlim(-plotscale,plotscale)
		ylim(0.0,60.0)
		gca().set_xticks([-plotscale,0.0,plotscale])
		gca().set_yticks([0.0,10.0,20.0,30.0,40.0, 50.0, 60.0])
		xlabel('kpc')

	[DM, mass,Xray1,Xray2,Xray3,SZ] = pysubs.ProjectEnzoData(pf,mass,phi=pi/2.0,theta=0.0,psi=0.0,zmin=zmin,zmax=zmax,DMProject=False)
	[Temp] = pysubs.ProjectEnzoTemp(pf,mass,phi=pi/2.0,theta=0.0,psi=0.0,zmin=zmin,zmax=zmax,MHD=False)	
	if rank == 0:
		#print "Y-axis, TrueTemp = ",TrueTemp.data, "\nTemp = ",Temp.data
		subplot(4,3,2)
		title('Y-axis Log DM Density',fontsize=10)	
		plot(mass.x-Centroids[1,1,0], (log10(DM.data[:,0])+log10(DM.data[:,1])+log10(DM.data[:,2]))/3.0, label = 'Time = %.2f Gy'%simtime, ls = ls)
		xlim(-plotscale,plotscale)
		ylim(-1.0,2.0)
		gca().set_xticks([-plotscale,0.0,plotscale])
		gca().set_yticks([-1.0,0.0,1.0,2.0])
		subplot(4,3,5)
		title('Y-axis Log Gas Density',fontsize=10)	
		plot(mass.x-Centroids[1,1,0], (log10(mass.data[:,0])+log10(mass.data[:,1])+log10(mass.data[:,2]))/3.0, label = 'Time = %.2f Gy'%simtime, ls = ls)
		xlim(-plotscale,plotscale)
		ylim(-2.0,1.0)
		gca().set_xticks([-plotscale,0.0,plotscale])
		subplot(4,3,8)
		title('Y-axis Gas Temp(keV)',fontsize=10)	
		plot(mass.x-Centroids[1,1,0], (Temp.data[:,0]+Temp.data[:,1]+Temp.data[:,2])/3.0 , label = 'Time = %.2f Gy'%simtime, ls = ls)
		xlim(-plotscale,plotscale)
		ylim(0.0,60.0)
		gca().set_xticks([-plotscale,0.0,plotscale])
		gca().set_yticks([0.0,10.0,20.0,30.0,40.0, 50.0, 60.0])
		xlabel('kpc')

	[DM, mass,Xray1,Xray2,Xray3,SZ] = pysubs.ProjectEnzoData(pf,mass,phi=pi/2.0,theta=pi/2.0,psi=0.0,zmin=zmin,zmax=zmax,DMProject=False)
	[Temp] = pysubs.ProjectEnzoTemp(pf,mass,phi=pi/2.0,theta=pi/2.0,psi=0.0,zmin=zmin,zmax=zmax,MHD=False)	
	if rank == 0:
		#print "Z-axis, TrueTemp = ",TrueTemp.data, "\nTemp = ",Temp.data
		subplot(4,3,3)
		title('Z-axis Log DM Density',fontsize=10)	
		plot(mass.x-Centroids[1,1,0], (log10(DM.data[:,0])+log10(DM.data[:,1])+log10(DM.data[:,2]))/3.0, label = 'Time = %.2f Gy'%simtime, ls = ls)
		xlim(-plotscale,plotscale)
		ylim(-1.0,2.0)
		gca().set_xticks([-plotscale,0.0,plotscale])
		gca().set_yticks([-1.0,0.0,1.0,2.0])
		subplot(4,3,6)
		title('Z-axis Log Gas Density',fontsize=10)	
		plot(mass.x-Centroids[1,1,0], (log10(mass.data[:,0])+log10(mass.data[:,1])+log10(mass.data[:,2]))/3.0, label = 'Time = %.2f Gy'%simtime, ls = ls)
		xlim(-plotscale,plotscale)
		ylim(-2.0,1.0)
		gca().set_xticks([-plotscale,0.0,plotscale])
		subplot(4,3,9)
		title('Z-axis Gas Temp(keV)',fontsize=10)	
		plot(mass.x-Centroids[1,1,0], (Temp.data[:,0]+Temp.data[:,1]+Temp.data[:,2])/3.0 , label = 'Time = %.2f Gy'%simtime, ls = ls)
		xlim(-plotscale,plotscale)
		ylim(0.0,60.0)
		gca().set_xticks([-plotscale,0.0,plotscale])
		gca().set_yticks([0.0,10.0,20.0,30.0,40.0, 50.0, 60.0])
		xlabel('kpc')
if rank == 0:
	legend(loc=1,bbox_to_anchor=[-0.2, -0.41])
	ltext = gca().get_legend().get_texts()
	setp(ltext, fontsize = 10)
	pp.savefig()
	pp.close()

#****************END MAIN PROGRAM*****************
