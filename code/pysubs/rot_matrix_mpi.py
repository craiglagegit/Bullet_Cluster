#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 9-Sep-11


#This program looks at a matrix of rotation angles on a databse with an IP of 4 kpc

import lageconfig # system specific path information
import sys
sys.path.append(lageconfig.bulletpath)
import pysubs_nas_18Sep12 as pysubs
from pylab import *
from mpi4py import MPI

#****************MAIN PROGRAM*****************

cmd=sys.argv

Z = float(cmd[1])
snapmin=int(cmd[2])
snapmax=int(cmd[3])
phi=0.0
filenamemassx1='fomrotmassx1_mpi.out'
filenamemass='fomrotmass_mpi.out'

numtheta = 16
numpsi = 32
numsnaps = snapmax - snapmin
numtot = numtheta * numpsi * numsnaps

comm = MPI.COMM_WORLD
rank = comm.Get_rank() # My Rank
size = comm.Get_size() # Total number of processors

step = numtot / size # Number done by each processor
mymin = step * rank
if rank != size - 1:
	mymax = mymin + step
else:
	mymax = numtot

for i in range(mymin, mymax):
	snap = int(i / (numtheta * numpsi))
	psimult = int((i - snap * (numtheta * numpsi)) / numtheta)
	thetamult = (i - snap * (numtheta * numpsi) - psimult * numtheta)
	snap = snap + snapmin
	#print "Size = %d, Rank = %d, Mymin = %d, Mymax = %d, Snap = %d, Thetamult = %d, Psimult = %d\n"%(size,rank,mymin,mymax,snap,thetamult,psimult)
	theta = pi/numtheta * thetamult
	psi = 2.0 * pi/numpsi * psimult
	
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

