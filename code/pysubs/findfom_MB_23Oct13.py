
#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 26-Nov-12

# Calculates best fom given a limited input range
# This version incorporates a limit on Vz


import lageconfig # system specific path information
import sys
sys.path.append(lageconfig.bulletpath)
import BulletConstants # Constants used in Bullet cluster simulations
import matplotlib
matplotlib.use("PDF")
import pysubs_nas_MB_23Oct13 as pysubs
from pylab import *
import time
from subprocess import *
import os
import copy

#****************MAIN PROGRAM*****************
filename = 'fominput'
file = open(filename,'r')
lines = file.readlines()
file.close()

Z = float(lines[0].strip().split()[1])
snapmin = int(lines[1].strip().split()[1])
snapmax = int(lines[2].strip().split()[1])
theta = float(lines[3].strip().split()[1])
psi = float(lines[4].strip().split()[1])
ConstrainPhi = True
phi=0.0
TFudge = float(lines[8].strip().split()[1])
MaxShift = float(lines[9].strip().split()[1])
SpectralIndex = float(lines[10].strip().split()[1])

Mask=(1,1,0,0,0,0)
(bestmassx1fom,xfom,bestmassx1snap,bestmassx1phi,bestmassx1theta,bestmassx1psi,bestmassx1time,counter) = pysubs.FindFom(snapmin,snapmax,Z,phi,theta,psi,ConstrainPhi=ConstrainPhi, TFudge=TFudge,Mask=Mask, MaxShift=MaxShift, SpectralIndex=SpectralIndex)
PlotSuffix = 'MassX1_NewRelVel2'
efom = pysubs.BestPlot(bestmassx1snap,Z,phi,bestmassx1theta,bestmassx1psi,PlotSuffix=PlotSuffix,TFudge=TFudge,FomLocate=True,ConstrainPhi=ConstrainPhi,Mask=Mask, MaxShift=MaxShift, SpectralIndex=SpectralIndex)
result="FOM=%.6f, XFOM=%.6f, EFOM=%.6f, BestSnap = %d, BestTime =  %.4f, BestPhi = %.3f, BestTheta = %.3f, BestPsi = %.3f\n"%(bestmassx1fom,xfom,efom,bestmassx1snap,bestmassx1time,bestmassx1phi,bestmassx1theta,bestmassx1psi)
fomfile=open('newfom_massx1.out','w')
fomfile.write(result)
fomfile.close

fomfile=open('FomFinished','w')
fomfile.write('All fom runs finished')
fomfile.close

#************END MAIN PROGRAM*************************
