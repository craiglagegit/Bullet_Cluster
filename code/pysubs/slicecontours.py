#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 28-Mar-13

#This program makes a series of plots using yt.
import lageconfig # system specific path information
import sys
import matplotlib
matplotlib.use("PDF")
sys.path.append(lageconfig.bulletpath)
from matplotlib.backends.backend_pdf import PdfPages
import  pysubs_nas_21Mar13 as pysubs
from scipy.ndimage import gaussian_filter, convolve
from subprocess import *
from pylab import *
import BulletConstants
import time

def GetFomParams(filename):
    file = open(filename,'r')
    lines = file.readlines()
    file.close()

    Z = float(lines[0].strip().split()[1])
    TFudge = float(lines[8].strip().split()[1])
    return [Z,TFudge]

def GetAngles(filename):
    file = open(filename,'r')
    line = file.readline()
    file.close()
    phi = float(line.split()[9].strip(","))
    theta = float(line.split()[12].strip(","))
    psi = float(line.split()[15].strip(","))
    return [phi,theta,psi]



#****************MAIN PROGRAM*****************
cmd =sys.argv
snapmin = int(cmd[1])
snapmax = int(cmd[2])
zmax = 100.0
zmin = -zmax
xmax = ymax = 1000.0
ymin = -ymax
xmin = -xmax
nx = ny = 512
[Z,TFudge] = GetFomParams('fominput')
[phi,theta,psi] = GetAngles('newfom_massx1.out')
phi=theta=psi=0.0
theta= 0.78
mass = pysubs.Array2d(xmin,xmax,nx,ymin,ymax,ny)
[dyy,dxx] = meshgrid(mass.y,mass.x)# Grid for plots
pp=PdfPages('Graph_MassContours_NoRot4.pdf')
gaslevels = [0.01,0.0178,0.0316,0.0562,0.1,0.178,0.316,0.562,1.0,1.78,3.16,5.62,10.0,17.8,31.6,56.2,100.0]#[10**-4.0,10**-3.5,10**-3.0,10**-2.5,-2.0]
dmlevels = [0.1,0.178,0.316,0.562,1.0,1.78,3.16,5.62,10.0,17.8,31.6,56.2,100.0,178.0,316.0,562.0,1000.0]#[10**-4.0,10**-3.5,10**-3.0,10**-2.5,-2.0]
for snap in range(snapmin,snapmax):
    pf = pysubs.GetPF(snap)
    simtime = pf.h.parameters['InitialTime'] * BulletConstants.TimeConversion # Put time in Gyears
    figure()
    subplots_adjust(hspace=0.8, wspace=0.4)
    suptitle('Bullet Cluster Collision, T = %.2f Gy'%simtime, fontsize=16)
    ApecData = pysubs.ReadLookups(Z) # Reads the APEC lookup tables.
    for grid in pf.h.grids:
        grid.set_field_parameter('ApecData',ApecData)
        grid.set_field_parameter('TFudge',TFudge)
        
    [DM, mass,Xray1,Xray2,Xray3,SZ] = pysubs.ProjectEnzoData(pf,mass,phi=phi,theta=theta,psi=psi,zmin=zmin,zmax=zmax,DMProject=False)
    DM.data=gaussian_filter(DM.data,3.0)
    mass.data=gaussian_filter(mass.data,3.0)
    subplot(1,2,1)
    title('DM Contours',fontsize=10)	
    contour(dxx,dyy,DM.data, dmlevels,colors='b', linewidths=2.0)
    xlim(xmin,xmax)
    ylim(ymin,ymax)
    gca().set_xticks([-500.0,0.0,500.0])

    subplot(1,2,2)
    title('Gas Contours',fontsize=10)	
    contour(dxx,dyy,mass.data, gaslevels,colors='b', linewidths=2.0)
    xlim(xmin,xmax)
    ylim(ymin,ymax)
    gca().set_xticks([-500.0,0.0,500.0])

    pp.savefig()
plt.clf()
pp.close()


#************END MAIN PROGRAM*************************

