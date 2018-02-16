#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 25-Sep-12


#This program plots the results of the rotattion matrix
import lageconfig # system specific path information
import sys
from subprocess import *
import matplotlib
matplotlib.use("PDF")
sys.path.append(lageconfig.bulletpath)
from pylab import *
import pysubs_nas_18Sep12 as pysubs
import time

#****************MAIN PROGRAM*****************
cmd = sys.argv
numtheta = 16
numpsi = 32
ncpus = int(cmd[1])
Z = float(cmd[2])
V = [3393,341,53]
delfile=Popen('rm -f fomrot/rot.out',shell=True)# Delete the large output file
Popen.wait(delfile)

for fname in ['mass', 'massx1']:
	psi=zeros([numpsi])
	theta=zeros([numtheta])
	fom=zeros([numpsi,numtheta])+10.0
	minfom = 100.0
	for i in range(ncpus):
		filename = 'fomrot/%s_mpi_%d.out'%(fname,i)
		file=open(filename,'r')
		lines = file.readlines()
		file.close
		for line in lines:
			thispsi = float(line.split()[2].strip(','))
			i = int(round(numpsi * thispsi / (2.0 * pi)))
			psi[i] = thispsi
			thistheta = float(line.split()[5].strip(','))
			j = int(round(numtheta * thistheta / pi))
			theta[j] = thistheta
			snap = int(line.split()[8].strip(','))
			thisfom = float(line.split()[9].strip('FOM='))
			R = pysubs.EulerAngles(-psi[i],-theta[j],0.0)
			Vp = dot(R,V)
			alpha = arcsin(Vp[2]/sqrt(dot(Vp,Vp)))*180.0/pi
			if abs(alpha) > 20.0:
				thisfom = 10.0

			if thisfom < fom[i,j]:
				fom[i,j] = thisfom
			if fom[i,j] < minfom:
				minfom = fom[i,j]
				minpsi = psi[i]
				mintheta = theta[j]
				minsnap = snap


	#alpha = abs(arcsin( sin(mintheta) * sin(minpsi))) * 180/pi
	plotfile = 'fomrot/Rot_results_new_%s.pdf'%fname
	figure()		

	title('IP = 25, Min FOM1 = %.3f, at Psi = %.3f, Theta = %.3f, Snap = %d, Angle to LOS =%.2f degrees\n'%(minfom,minpsi,mintheta,minsnap, alpha),fontsize = 12)
	xx,yy = meshgrid(theta,psi)
	levels = linspace(1.0,7.0,21)	
	contourf(xx,yy,fom,levels)
	xlim(0,3.14)
	ylim(0,6.28)
	colorbar()
	plt.savefig(plotfile)
	print 'Wrote file', plotfile
	plt.clf()
	sys.exit()
	snap = minsnap
	Phi=0.0
	Theta=mintheta
	Psi=minpsi
	ConstrainPhi=True
	simtime = minsnap / 100.0 * 0.9777
	result="FOM=%.6f, BestSnap = %d, BestTime =  %.4f, BestPsi = %.3f, BestTheta = %.3f\n"%(minfom,minsnap,simtime,minpsi,mintheta)
	if fname == 'mass':
		PlotSuffix = 'Mass'
		Mask=(1,0,0,0,0,0)
		fomfile=open('newfom_mass.out','w')
		fomfile.write(result)
		fomfile.close

	else:
		PlotSuffix = 'MassX1'
		Mask=(1,1,0,0,0,0)
		fomfile=open('newfom_massx1.out','w')
		fomfile.write(result)
		fomfile.close

	pysubs.BestPlot(snap,Z,Phi,Theta,Psi,PlotSuffix=PlotSuffix,FomLocate=True,ConstrainPhi=ConstrainPhi,Mask=Mask)


#************END MAIN PROGRAM*************************

