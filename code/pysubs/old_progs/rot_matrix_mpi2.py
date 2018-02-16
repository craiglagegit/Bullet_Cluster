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
import time
from subprocess import *
#****************MAIN PROGRAM*****************
phi=0.0
numtheta  = 16
numpsi = 32

comm = MPI.COMM_WORLD
rank = comm.Get_rank() # My Rank
size = comm.Get_size() # Total number of processors
if rank == 0:
	done = ones([size],dtype=int) # Array to check if we're done
	cmd=sys.argv 
	Z = float(cmd[1])
	snapmin=int(cmd[2])
	snapmax=int(cmd[3])
	for i in range(1,size):
		comm.send(Z,dest=i,tag=1000+i)
		comm.send(snapmin,dest=i,tag=2000+i)
		comm.send(snapmax,dest=i,tag=3000+i)
else:
	Z = comm.recv(source=0,tag=1000+rank)
	snapmin = comm.recv(source=0,tag=2000+rank)
	snapmax = comm.recv(source=0,tag=3000+rank)
numsnaps = snapmax - snapmin
numtot = numtheta * numpsi * numsnaps

filenamemassx1='fomrot/massx1_mpi_%d.out'%rank
filenamemass='fomrot/mass_mpi_%d.out'%rank

step = numtot / size # Number done by each processor
mymin = step * rank
if rank != size - 1:
	mymax = mymin + step
else:
	mymax = numtot

for i in range(mymin, mymax):
	
	#memusage = Popen("ps -u clage -o rss | awk '{sum+=$1} END {print sum}'",shell=True,stdout=PIPE)
        #print "Rank = %d, i = %d, Memusage="%(rank,i), int(memusage.communicate()[0].split('\n')[0])
	snap = int(i / (numtheta * numpsi))
	psimult = int((i - snap * (numtheta * numpsi)) / numtheta)
	thetamult = (i - snap * (numtheta * numpsi) - psimult * numtheta)
	snap = snap + snapmin
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
if rank != 0:
	comm.send(0,dest=0,tag=rank)
else:
	while done.sum() != 0:
		done[0] = 0
		for i in range(1,size):
			done[i] = comm.recv(source=i, tag=i)
			time.sleep(0.5)
	finfile=open('fomrot/Rot_Finished','w')
	finfile.write('Finished')
	finfile.close


#************END MAIN PROGRAM*************************

