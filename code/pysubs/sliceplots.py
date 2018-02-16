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
zmax = 500.0
zmin = -zmax
xmax = 3000.0
ymax = 1500.0
ymin = -ymax
xmin = -xmax
nx = 512
ny = 256
ploticks = 2000.0
[Z,TFudge] = GetFomParams('fominput')
[phi,theta,psi] = GetAngles('newfom_massx1.out')
mass = pysubs.Array2d(xmin,xmax,nx,ymin,ymax,ny)
[dyy,dxx] = meshgrid(mass.y,mass.x)# Grid for plots
pp=PdfPages('Graph_Slices.pdf')

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
    [Temp,BMag] = pysubs.ProjectEnzoTemp(pf,mass,phi=phi,theta=theta,psi=psi,zmin=zmin,zmax=zmax)	

    subplot(2,2,1,aspect=2)
    title('Log DM Density',fontsize=10)	
    contourf(dxx,dyy,log10(DM.data))
    for c in contourf(dxx,dyy,log10(DM.data)).collections:
        c.set_linewidth(0.1)
        c.set_alpha(1.0)
    #colorbar()
    gca().set_xticks([-ploticks,0.0,ploticks])

    subplot(2,2,2,aspect=2)
    title('Log Gas Density',fontsize=10)	
    contourf(dxx,dyy,log10(mass.data))
    for c in contourf(dxx,dyy,log10(mass.data)).collections:
        c.set_linewidth(0.1)
        c.set_alpha(1.0)
    #colorbar()
    gca().set_xticks([-ploticks,0.0,ploticks])

    subplot(2,2,3,aspect=2)
    title('Gas Temperature (keV)',fontsize=10)	
    contourf(dxx,dyy,Temp.data)
    for c in contourf(dxx,dyy,Temp.data).collections:
        c.set_linewidth(0.1)
        c.set_alpha(1.0)
    colorbar()
    gca().set_xticks([-ploticks,0.0,ploticks])

    subplot(2,2,4,aspect=2)
    title('Magnetic Field Magnitude ($\mu G$)',fontsize=10)	
    contourf(dxx,dyy,BMag.data*1.0E6)
    for c in contourf(dxx,dyy,BMag.data*1.0E6).collections:
        c.set_linewidth(0.1)
        c.set_alpha(1.0)
    colorbar()
    gca().set_xticks([-ploticks,0.0,ploticks])
    pp.savefig()
plt.clf()
pp.close()


#************END MAIN PROGRAM*************************

