#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 27-Mar-12


#This is a standalone python program to run findfom in batch mode

import lageconfig # system specific path information
import sys
import matplotlib
matplotlib.use("PDF")
sys.path.append(lageconfig.bulletpath)
import pysubs_nas_21Mar13 as pysubs
import time


#****************MAIN PROGRAM*****************

cmd=sys.argv

snap=int(cmd[1])
Z=float(cmd[2])
phi=float(cmd[3])
theta=float(cmd[4])
psi=float(cmd[5])
TFudge=float(cmd[6])

pysubs.BestPlot(snap,Z,phi=phi,theta=theta,psi=psi,PlotSuffix='MassX1_Test3',FomLocate=True,TFudge=TFudge,ConstrainPhi=True, Mask=(1,1,0,0,0,0))


#************END MAIN PROGRAM*************************

