
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
import pysubs_nas_test_17Jun13 as pysubs
from pylab import *
import time
from subprocess import *
import os
import copy



#****************MAIN PROGRAM*****************
cenfilename = 'centroids.out'
cmd=sys.argv
snapmin=int(cmd[1])
snapmax=int(cmd[2])
lastsimtime = 1000.0 # Nonsense value for first time through
LastBulletDMPos = [0.0,0.0,0.0]
LastMainDMPos = [0.0,0.0,0.0]
for snap in range(snapmin,snapmax):
	pf = pysubs.GetPF(snap)
	simtime = pf.h.parameters['InitialTime'] * BulletConstants.TimeConversion # Put time in Gyears
	CenOut = pysubs.FindEnzoCentroids(pf) # CenOut = [NumPart, Masses, Centroids, MassWithin250K]
	BulletDMPos = CenOut[2][0]
	MainDMPos = CenOut[2][1]
	dt = simtime - lastsimtime
	lastsimtime = simtime
	BulletDMVel = (BulletDMPos-LastBulletDMPos)/dt
	MainDMVel = (MainDMPos-LastMainDMPos)/dt
	LastBulletDMPos = BulletDMPos 
	LastMainDMPos = MainDMPos 
	RelVel = BulletDMVel - MainDMVel
	if snap > snapmin: # RelVel not valid for first snap
		cenfile = open(cenfilename, 'a')
		result = "Time = %.3f BulletDMPos = %.2f %.2f %.2f MainDMPos = %.2f %.2f %.2f RelVel = %.2f %.2f %.2f\n"%(simtime, BulletDMPos[0], BulletDMPos[1], BulletDMPos[2], MainDMPos[0], MainDMPos[1], MainDMPos[2], RelVel[0], RelVel[1], RelVel[2])
		cenfile.write(result)
		cenfile.close()

#************END MAIN PROGRAM*************************
