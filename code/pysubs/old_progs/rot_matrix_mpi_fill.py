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
cmd=sys.argv 
Z = float(cmd[1])
snapmin = int(cmd[2])
snapmax = int(cmd[3])
phi=0.0
numtheta  = 16
numpsi = 32

comm = MPI.COMM_WORLD
rank = comm.Get_rank() # My Rank
size = comm.Get_size() # Total number of processors

if rank == 0:

	file = open('newfom_mass.out','r')
	line = file.readline()
	file.close()

	snapmin = int(line.strip().split()[3].split(',')[0]) - 2
	snapmax = min(snapmax, snapmin + 5)

	print "Snap finding done - Global snapmin = %d, Global snapmax = %d\n"%(snapmin,snapmax)
# Now we send the beginning and end snaps to everybody
snapmin = comm.bcast(snapmin,root=0)
snapmax = comm.bcast(snapmax,root=0)


# Now we start the angle matrices
numsnaps = snapmax - snapmin
numtot = numtheta * numpsi * numsnaps

filenamemassx1='fomrot/massx1_mpi_%d.out'%rank
filenamemass='fomrot/mass_mpi_%d.out'%rank

step = float(numtot) / float(size) # Number done by each processor
mymin = int(step * rank)
if rank != size - 1:
	mymax = int(step * (rank + 1))
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
	alpha = abs(arcsin( sin(theta) * sin(psi)))
	if alpha > pi / 4.0 or alpha <= pi / 6.0: # Don't bother if angle to LOS is > 45 degrees or < 30 degrees.
		continue
	else:
		(fommassx1,fommass)=pysubs.SimpleFom(snap,Z,phi=phi,theta=theta,psi=psi,ConstrainPhi=True)

	result='Psi = %.4f, Theta = %.4f, Snap = %d, FOM=%.4f\n'%(psi,theta,snap,fommassx1)
	fomfile=open(filenamemassx1,'a')
	fomfile.write(result)
	fomfile.close

	result='Psi = %.4f, Theta = %.4f, Snap = %d, FOM=%.4f\n'%(psi,theta,snap,fommass)
	fomfile=open(filenamemass,'a')
	fomfile.write(result)
	fomfile.close

# This writes out a file when we're all done
if rank != 0:
	comm.send(0,dest=0,tag=rank)
else:
	while done.sum() != 0:
		done[0] = 0
		for i in range(1,size):
			done[i] = comm.recv(source=i, tag=i)
			time.sleep(1.0)
	finfile=open('fomrot/Rot_Finished','w')
	finfile.write('Finished')
	finfile.close


#************END MAIN PROGRAM*************************

