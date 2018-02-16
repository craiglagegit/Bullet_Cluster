#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 16-Oct-12


import lageconfig # system specific path information
import sys
sys.path.append(lageconfig.bulletpath)
import BulletConstants # Constants used in Bullet cluster simulations
import matplotlib
matplotlib.use("PDF")
import pysubs_nas_23Oct12 as pysubs
from pylab import *
from mpi4py import MPI
import time
from subprocess import *
import os

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
thetamin = float(cmd[4])
thetamax = float(cmd[5])
psimin = float(cmd[6])
psimax = float(cmd[7])
ConstrainPhi = True
phi=0.0
MaxLOSAngle = 15.0 # Maximum angle to line-of-sight in degrees.

comm = MPI.COMM_WORLD
rank = comm.Get_rank() # My Rank
size = comm.Get_size() # Total number of processors
filenamex1='fomrot/x1_mpi_%d.out'%rank
filenamemass='fomrot/mass_mpi_%d.out'%rank

# First, processor zero goes through the angles and counts how many valid angles there are
angles = list()
if rank == 0:
        V = [0.994894, 0.099724, 0.0155452] # Approximate bullet velocity vector at time of best fit
	for Theta in range(int(100*thetamin), int(100*thetamax), 20):
		theta = Theta / 100.0
		for Psi in range(int(100*psimin), int(100*psimax), 20):
			psi = Psi / 100.0
			R = pysubs.EulerAngles(-psi,-theta,0.0)
			Vp = dot(R,V)
			alpha = arcsin(Vp[2])*180.0/pi
			if abs(alpha) < MaxLOSAngle:
				angles.append([theta,psi,100.0,100.0])
			else:
				continue

angles = comm.bcast(angles,root=0) # Next, send everybody the angles to be run
numangles = len(angles)
step = float(numangles) / float(size) # Number done by each processor
mymin = int(step * rank)
if rank != size - 1:
	mymax = int(step * (rank + 1))
else:
	mymax = numangles
print "Processor #%d, numangles = %d, mymin = %d, mymax = %d\n"%(rank,numangles,mymin,mymax)			
if rank == 0:
	bestmassfom = 100.0
	bestmasssnap = 1
	bestmasstime = 0.0
	bestmasstheta = 0.0
	bestmasspsi = 0.0
	bestx1fom = 100.0
	bestx1snap = 1
	bestx1time = 0.0
	bestx1theta = 0.0
	bestx1psi = 0.0

# Now, we cycle through the snaps, calculating the FOM for all valid angles

for snap in range(snapmin, snapmax):
	if rank == 0:
		while True:
			if EnzoRunComplete(snap):
				break
		print "Processor %d, Enzo run %d is complete\n"%(rank,snap)
		sys.stdout.flush()
		try:
			pf = pysubs.GetPF(snap)
			simtime = pf.h.parameters['InitialTime'] * BulletConstants.TimeConversion # Put time in Gyears
			[NumPart, Masses, Centroids] = pysubs.FindEnzoCentroids(pf)
			#print "Processor #%d, snap = %d"%(rank,snap)," Centroids=",Centroids 
			sys.stdout.flush()
			CenterSep = 0.0
			for i in range(3):
				CenterSep = CenterSep + (Centroids[0][i] - Centroids[1][i])**2
			CenterSep = sqrt(CenterSep)
		except:
			continue
		print "Processor #%d, snap = %d, Sep = %f\n"%(rank,snap,CenterSep)
		sys.stdout.flush()
		if CenterSep < 500.0 or Centroids[0][0] < Centroids[1][0]: #Bullet hasn't yet reached 500 kpc past main
			SnapFlag = 'Wait'
		elif CenterSep > 750.0 and Centroids[0][0] > Centroids[1][0]: #Bullet has passed main and sep > 750 kpc means we're done
			SnapFlag = 'Break'
		else:
			SnapFlag = 'Run'
	else:
		SnapFlag = None
	SnapFlag = comm.bcast(SnapFlag,root=0) # Communicate CenterSep results to everyone
	if SnapFlag == 'Break':
		break
	elif SnapFlag == 'Wait':
		continue
	# If SnapFlag == 'Run', then we start the angle runs, logging them in case we want to look later
	for i in range(mymin, mymax):
		theta = angles[i][0]
		psi = angles[i][1]
		(fomx1,fommass)=pysubs.SimpleFom(snap,Z,phi=phi,theta=theta,psi=psi,ConstrainPhi=ConstrainPhi)
		angles[i][2] = fommass
		angles[i][3] = fomx1
		result='Psi = %.4f, Theta = %.4f, Snap = %d, FOM=%.4f\n'%(psi,theta,snap,fomx1)
		fomfile=open(filenamex1,'a')
		fomfile.write(result)
		fomfile.close

		result='Psi = %.4f, Theta = %.4f, Snap = %d, FOM=%.4f\n'%(psi,theta,snap,fommass)
		fomfile=open(filenamemass,'a')
		fomfile.write(result)
		fomfile.close
	# Now we collect the results
	if rank != 0:
		comm.send(angles, dest=0, tag=100*snap+rank)
	# Then find the best results 
	if rank == 0:
		for proc in range(size):
			if proc == 0:
				procangles = copy(angles)
			else:
				procangles = comm.recv(source=proc, tag=100*snap+proc)
			for angle in procangles:
				if angle[2] < bestmassfom:
					bestmassfom = angle[2]
					bestmasssnap = snap
					bestmasstime = simtime
					bestmasstheta = angle[0]
					bestmasspsi = angle[1]
				if angle[3] < bestx1fom:
					bestx1fom = angle[3]
					bestx1snap = snap
					bestx1time = simtime
					bestx1theta = angle[0]
					bestx1psi = angle[1]
		print "Processor #%d, snap = %d, bestmassfom = %f, bestmasssnap = %d, bestmasstheta = %f, bestmasspsi = %f\n"%(rank,snap,bestmassfom, bestmasssnap, bestmasstheta, bestmasspsi)
		print "Processor #%d, snap = %d, bestx1fom = %f, bestx1snap = %d, bestx1theta = %f, bestx1psi = %f\n"%(rank,snap,bestx1fom, bestx1snap, bestx1theta, bestx1psi)
		sys.stdout.flush()

# This takes the best result and then does a finer angle matrix around it, then writes out a file and makes plots when we're all done
if rank == 0:
	comm.send([bestx1fom,bestx1snap,bestx1time,bestx1theta,bestx1psi], dest=1, tag=2000)
	Mask=(1,0,0,0,0,0)
	(bestmassfom,bestmasssnap,bestmasspsi,bestmasstheta,bestmasstime,counter) = pysubs.FindFom(max(snapmin,bestmasssnap-1),min(snapmax,bestmasssnap+2),Z,phi,bestmasstheta,bestmasspsi,ConstrainPhi=ConstrainPhi, Mask=Mask)
	result="FOM=%.6f, BestSnap = %d, BestTime =  %.4f, BestPsi = %.3f, BestTheta = %.3f\n"%(bestmassfom,bestmasssnap,bestmasstime,bestmasspsi,bestmasstheta)
	PlotSuffix = 'Mass'
	pysubs.BestPlot(bestmasssnap,Z,phi,bestmasstheta,bestmasspsi,PlotSuffix=PlotSuffix,FomLocate=True,ConstrainPhi=ConstrainPhi,Mask=Mask)
	fomfile=open('newfom_mass.out','w')
	fomfile.write(result)
	fomfile.close
	print "Processor %d - Normal Exit\n"%rank

elif rank == 1:
	[bestx1fom,bestx1snap,bestx1time,bestx1theta,bestx1psi] = comm.recv(source=0, tag=2000)
	Mask=(0,1,0,0,0,0)
	(bestx1fom,bestx1snap,bestx1psi,bestx1theta,bestx1time,counter) = pysubs.FindFom(max(snapmin,bestx1snap-1),min(snapmax,bestx1snap+2),Z,phi,bestx1theta,bestx1psi,ConstrainPhi=ConstrainPhi, Mask=Mask)
	result="FOM=%.6f, BestSnap = %d, BestTime =  %.4f, BestPsi = %.3f, BestTheta = %.3f\n"%(bestx1fom,bestx1snap,bestx1time,bestx1psi,bestx1theta)
	PlotSuffix = 'X1'
	pysubs.BestPlot(bestx1snap,Z,phi,bestx1theta,bestx1psi,PlotSuffix=PlotSuffix,FomLocate=True,ConstrainPhi=ConstrainPhi,Mask=Mask)
	fomfile=open('newfom_x1.out','w')
	fomfile.write(result)
	fomfile.close
	print "Processor %d - Normal Exit\n"%rank

else:
	print "Processor %d - Normal Exit\n"%rank

#************END MAIN PROGRAM*************************

