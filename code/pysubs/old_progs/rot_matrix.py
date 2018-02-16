#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 9-Sep-11


#This program looks at a matrix of rotation angles on a databse with an IP of 4 kpc

import lageconfig # system specific path information
import sys
sys.path.append(lageconfig.bulletpath)
import pysubs_nas_18Sep12 as pysubs
from pylab import *


#****************MAIN PROGRAM*****************

cmd=sys.argv

Z = float(cmd[1])
snapmin=int(cmd[2])
snapmax=int(cmd[3])
phi=0.0
filenamemassx1='fomrotmassx1_new.out'
filenamemass='fomrotmass_new.out'

for thetamult in range(16):
	theta = pi/16.0 * thetamult
	for psimult in range(32):
		for snap in range(snapmin,snapmax):
			psi = pi/16.0 * psimult

			(fommassx1,fommass)=pysubs.SimpleFom(snap,Z,phi=phi,theta=theta,psi=psi,ConstrainPhi=True)

			result='Psi = %.4f, Theta = %.4f, Snap = %d, FOM=%.4f\n'%(psi,theta,snap,fommassx1)
			fomfile=open(filenamemassx1,'a')
			fomfile.write(result)
			fomfile.close

			result='Psi = %.4f, Theta = %.4f, Snap = %d, FOM=%.4f\n'%(psi,theta,snap,fommass)
			fomfile=open(filenamemass,'a')
			fomfile.write(result)
			fomfile.close


#************END MAIN PROGRAM*************************

