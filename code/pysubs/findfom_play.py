#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 27-Mar-12


#This is a standalone python program to run findfom in batch mode

import lageconfig # system specific path information
import sys
sys.path.append(lageconfig.bulletpath)
import pysubs_mhd_new_apec_play as pysubs
import time



#****************MAIN PROGRAM*****************

cmd=sys.argv

snapmin=int(cmd[1])
snapmax=int(cmd[2])
Z=float(cmd[4])
phi=float(cmd[4])
theta=float(cmd[5])
psi=float(cmd[6])
PlotSuffix = cmd[7]
if cmd[8] == 'IncludeTemp':
	IncludeTemp = True
else:
	IncludeTemp = False

"""
if PlotSuffix == None:
	if IncludeTemp:
		filename='newfom_apec.out'
	else:
		filename='newfom_apec_noT.out'
	result="FOM=%.6f at time =  %.4f"%(bestfom,besttime)
	fomfile=open(filename,'w')
	fomfile.write(result)
	fomfile.close
"""
#start = time.time()
#(fom,snap,theta,phi,simtime,counter) = pysubs.FindFom2(snapmin,snapmax,Z,phi,theta,psi)

#print fom, snap, theta, phi
#elapsed = time.time()-start
#print "Elapsed time = %.2f"%elapsed
#pysubs.BestPlot(snap,Z,phi,theta,psi=0.0,PlotSuffix=PlotSuffix,FomLocate=True,IncludeTemp=IncludeTemp)
pysubs.BestPlot(snapmin,Z,phi,theta,psi=0.0,PlotSuffix=PlotSuffix,FomLocate=True,IncludeTemp=IncludeTemp)

#************END MAIN PROGRAM*************************

