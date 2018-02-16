#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 27-Mar-12


#This is a standalone python program to run findfom in batch mode

import lageconfig # system specific path information
import sys
import matplotlib
matplotlib.use("PDF")
sys.path.append(lageconfig.bulletpath)
import pysubs_mhd_apec_27Mar12 as pysubs
import time



#****************MAIN PROGRAM*****************

cmd=sys.argv

snapmin=int(cmd[1])
snapmax=int(cmd[2])
Z=float(cmd[3])
phi=float(cmd[4])
theta=float(cmd[5])
psi=float(cmd[6])
PlotSuffix = 'APEC'
IncludeTemp = False

start = time.time()
(fom,snap,phi,theta,simtime,counter) = pysubs.FindFom2(snapmin,snapmax,Z,phi,theta,psi)
	
result="FOM=%.6f, BestSnap = %d, BestTime =  %.4f, BestPhi = %.3f, BestTheta = %.3f\n"%(fom,snap,simtime,phi,theta)

elapsed = time.time()-start
print result
print "Elapsed time = %.2f"%elapsed

pysubs.BestPlot(snap,Z,phi,theta,psi=0.0,PlotSuffix=PlotSuffix,FomLocate=True,IncludeTemp=IncludeTemp)

fomfile=open('newfom_apec_noT.out','w')
fomfile.write(result)
fomfile.close

#************END MAIN PROGRAM*************************

