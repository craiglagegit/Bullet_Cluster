#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 9-Sep-11


#This program looks at a matrix of rotation angles on a databse with an IP of 4 kpc

from pylab import *
import sys


#****************MAIN PROGRAM*****************
cmd = sys.argv
numtheta = 16
numpsi = 32
psi=zeros([numpsi])
theta=zeros([numtheta])
fom=zeros([numpsi,numtheta])+10.0
run = int(cmd[1])
ncpus = int(cmd[2])
fname = cmd[3]
minfom = 100.0
for i in range(ncpus):
	filename = 'ddfiles/run%d/fomrot/%s_mpi_%d.out'%(run,fname,i)
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
		if thisfom < fom[i,j]:
			fom[i,j] = thisfom
		if fom[i,j] < minfom:
			minfom = fom[i,j]
			minpsi = psi[i]
			mintheta = theta[j]
			minsnap = snap

alpha = abs(arcsin( sin(mintheta) * sin(minpsi))) * 180/pi	
title('IP = 25, Min FOM1 = %.3f, at Psi = %.3f, Theta = %.3f, Snap = %d, Angle to LOS =%.2f degrees\n'%(minfom,minpsi,mintheta,minsnap, alpha),fontsize = 12)
xx,yy = meshgrid(theta,psi)
levels = linspace(1.0,7.0,21)	
contourf(xx,yy,fom,levels)
xlim(0,3.14)
ylim(0,6.28)
colorbar()
show()

#************END MAIN PROGRAM*************************

