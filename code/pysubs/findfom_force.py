#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 2-Sep-11


#This is a standalone python program to run findfom in batch mode

import lageconfig # system specific path information
import sys
import matplotlib
matplotlib.use("PDF")
sys.path.append(lageconfig.bulletpath)
import pysubs_mhd_force as pysubs


#****************MAIN PROGRAM*****************
cmd=sys.argv

toppath=lageconfig.toppath
dir='./'
snapmin=int(cmd[1])
snapmax=int(cmd[2])
counter=1
phi=float(cmd[3])
theta=float(cmd[4])
psi=float(cmd[5])

[besttime,bestfom]=pysubs.FindFom(toppath,dir,snapmin,snapmax,counter,phi,theta,psi,simulator='Enzo',PlotSuffix='Test',FomLocate=True)

#[besttime,bestfom]=pysubs.FindFom(toppath,dir,snapmin,snapmax,counter,phi,theta,psi,simulator='Enzo',PlotSuffix='Test_%d'%snapmin,FomLocate=True)

#************END MAIN PROGRAM*************************

