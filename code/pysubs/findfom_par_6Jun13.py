
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
import pysubs_nas_par_6Jun13 as pysubs
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
filename = 'fominput'
file = open(filename,'r')
lines = file.readlines()
file.close()

Z = float(lines[0].strip().split()[1])
snapmin = int(lines[1].strip().split()[1])
snapmax = int(lines[2].strip().split()[1])
theta = float(lines[3].strip().split()[1])
psi = float(lines[4].strip().split()[1])
RelVel0 = float(lines[5].strip().split()[1])
RelVel1 = float(lines[6].strip().split()[1])
RelVel2 = float(lines[7].strip().split()[1])
RelVel=[RelVel0,RelVel1,RelVel2]
ConstrainPhi = True
phi=0.0
TFudge = float(lines[8].strip().split()[1])
MaxShift = float(lines[9].strip().split()[1])
SpectralIndex = float(lines[10].strip().split()[1])

Mask=(1,1,0,0,0,0)
(bestmassx1fom,xfom,bestmassx1snap,bestmassx1phi,bestmassx1theta,bestmassx1psi,bestmassx1time,counter) = pysubs.FindFom(snapmin,snapmax,Z,phi,theta,psi,RelVel,ConstrainPhi=ConstrainPhi, TFudge=TFudge,Mask=Mask, MaxShift=MaxShift, SpectralIndex=SpectralIndex)
result="FOM=%.6f, XFOM=%.6f, BestSnap = %d, BestTime =  %.4f, BestPhi = %.3f, BestTheta = %.3f, BestPsi = %.3f\n"%(bestmassx1fom,xfom,bestmassx1snap,bestmassx1time,bestmassx1phi,bestmassx1theta,bestmassx1psi)
PlotSuffix = 'MassX1'
pysubs.BestPlot(bestmassx1snap,Z,phi,bestmassx1theta,bestmassx1psi,PlotSuffix=PlotSuffix,TFudge=TFudge,FomLocate=True,ConstrainPhi=ConstrainPhi,Mask=Mask, MaxShift=MaxShift, SpectralIndex=SpectralIndex)
fomfile=open('newfom_massx1.out','w')
fomfile.write(result)
fomfile.close

fomfile=open('FomFinished','w')
fomfile.write('All fom runs finished')
fomfile.close

#************END MAIN PROGRAM*************************
