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
import pysubs_nas_10Dec12 as pysubs
from pylab import *
import time
from subprocess import *
import os
import copy

def EnzoRunComplete(snap):
	# This subroutine checks to see if the Enzo file exists and is complete
	NumTries = 0
	LastFileSize = 0
	if snap<10:
		EnzoWaitFile = "DD000"+str(snap)
	elif snap<100:
		EnzoWaitFile = "DD00"+str(snap)
	elif snap<1000:
		EnzoWaitFile = "DD0"+str(snap)
	else:
		EnzoWaitFile = "DD"+str(snap)

	while NumTries < 4:
		try:
			NumTries = NumTries + 1
			FileSize = int(Popen('du -cs '+EnzoWaitFile, shell=True, stdout = PIPE).communicate()[0].split()[0])
			#print "WaitFile = %s, FileSize = %d, LastFileSize = %d\n"%(self.waitfile,FileSize,LastFileSize)
			if FileSize > 0 and FileSize == LastFileSize :
                            time.sleep(5.0)
			    return True
			else:
			    LastFileSize = FileSize
			    time.sleep(1.0)
			    continue
		except OSError:
			time.sleep(1.0)
	return False			


#****************MAIN PROGRAM*****************
cmd=sys.argv 
Z = float(cmd[1])
snapmin = int(cmd[2])
snapmax = int(cmd[3])
theta = float(cmd[4])
psi = float(cmd[5])
RelVel0 = float(cmd[6])
RelVel1 = float(cmd[7])
RelVel2 = float(cmd[8])
RelVel=[RelVel0,RelVel1,RelVel2]
ConstrainPhi = True
phi=0.0

Mask=(1,0,0,0,0,0)
(bestmassfom,bestmasssnap,bestmassphi,bestmasstheta,bestmasspsi,bestmasstime,counter) = pysubs.FindFom(snapmin,snapmax,Z,phi,theta,psi,RelVel,ConstrainPhi=ConstrainPhi, Mask=Mask)
result="FOM=%.6f, BestSnap = %d, BestTime =  %.4f, BestPhi = %.3f, BestTheta = %.3f, BestPsi = %.3f\n"%(bestmassfom,bestmasssnap,bestmasstime,bestmassphi,bestmasstheta,bestmasspsi)
PlotSuffix = 'Mass'
pysubs.BestPlot(bestmasssnap,Z,phi,bestmasstheta,bestmasspsi,PlotSuffix=PlotSuffix,FomLocate=True,ConstrainPhi=ConstrainPhi,Mask=Mask)
fomfile=open('newfom_mass.out','w')
fomfile.write(result)
fomfile.close

Mask=(1,1,0,0,0,0)
(bestmassx1fom,bestmassx1snap,bestmassx1phi,bestmassx1theta,bestmassx1psi,bestmassx1time,counter) = pysubs.FindFom(snapmin,snapmax,Z,phi,theta,psi,RelVel,ConstrainPhi=ConstrainPhi, Mask=Mask)
result="FOM=%.6f, BestSnap = %d, BestTime =  %.4f, BestPhi = %.3f, BestTheta = %.3f, BestPsi = %.3f\n"%(bestmassx1fom,bestmassx1snap,bestmassx1time,bestmassx1phi,bestmassx1theta,bestmassx1psi)
PlotSuffix = 'MassX1'
pysubs.BestPlot(bestmassx1snap,Z,phi,bestmassx1theta,bestmassx1psi,PlotSuffix=PlotSuffix,FomLocate=True,ConstrainPhi=ConstrainPhi,Mask=Mask)
fomfile=open('newfom_massx1.out','w')
fomfile.write(result)
fomfile.close

fomfile=open('FomFinished','w')
fomfile.write('All fom runs finished')
fomfile.close


#************END MAIN PROGRAM*************************

