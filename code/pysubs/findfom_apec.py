#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 2-Sep-11


#This is a standalone python program to run findfom in batch mode

import lageconfig # system specific path information
import sys
import matplotlib
matplotlib.use("PDF")
sys.path.append(lageconfig.bulletpath)
import pysubs_mhd_new_apec as pysubs


#****************MAIN PROGRAM*****************
cmd=sys.argv

dir='./'
snapmin=int(cmd[1])
snapmax=int(cmd[2])
counter=1
Z=float(cmd[3])
phi=float(cmd[4])
theta=float(cmd[5])
psi=float(cmd[6])
PlotSuffix = cmd[7]
if PlotSuffix == 'None':
	PlotSuffix = None
if cmd[7] == 'IncludeTemp':
	IncludeTemp = True
else:
	IncludeTemp = False

[besttime,bestfom]=pysubs.FindFom(dir,snapmin,snapmax,counter,Z,phi,theta,psi,simulator='Enzo',PlotSuffix=PlotSuffix,IncludeTemp=IncludeTemp,FomLocate=True)

if PlotSuffix == None:
	if IncludeTemp:
		filename='newfom_apec.out'
	else:
		filename='newfom_apec_noT.out'
	result="FOM=%.6f at time =  %.4f"%(bestfom,besttime)
	fomfile=open(filename,'w')
	fomfile.write(result)
	fomfile.close

#************END MAIN PROGRAM*************************

