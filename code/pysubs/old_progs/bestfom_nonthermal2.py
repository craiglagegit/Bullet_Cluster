
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
import pysubs_nas_11Feb13_newshock as pysubs
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
snap = int(lines[1].strip().split()[1])
phi = float(lines[2].strip().split()[1])
theta = float(lines[3].strip().split()[1])
psi = float(lines[4].strip().split()[1])
RelVel0 = float(lines[5].strip().split()[1])
RelVel1 = float(lines[6].strip().split()[1])
RelVel2 = float(lines[7].strip().split()[1])
RelVel=[RelVel0,RelVel1,RelVel2]
ConstrainPhi = True

TFudge = float(lines[8].strip().split()[1])

#Mask=(1,0,0,0,0,0)

#PlotSuffix = 'Mass_NonThermal2'
#pysubs.BestPlot(snap,Z,phi,theta,psi,PlotSuffix=PlotSuffix,TFudge=TFudge,FomLocate=True,ConstrainPhi=ConstrainPhi,Mask=Mask)


Mask=(1,1,0,0,0,0)

PlotSuffix = 'MassX1_NewShock'
pysubs.BestPlot(snap,Z,phi,theta,psi,PlotSuffix=PlotSuffix,TFudge=TFudge,FomLocate=True,ConstrainPhi=ConstrainPhi,Mask=Mask)



#************END MAIN PROGRAM*************************

