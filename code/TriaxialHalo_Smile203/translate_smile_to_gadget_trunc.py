#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 10-Aug-12


#This program translates a smaile nbody file to a Gadget input file
#Added a random magnetic field read-in

import lageconfig # system specific path information
import sys
sys.path.append(lageconfig.bulletpath)
import pysubs_mhd_apec_mass_x1_3May12 as pysubs
import BulletConstants
from pylab import *

#****************MAIN PROGRAM*****************

cmd=sys.argv
inputfile = cmd[1]
outputfile = cmd[2]
M200 = float(cmd[3]) # Mass
R200 = float(cmd[4]) # Radius
mfactor = M200 / .988287
rfactor = R200
G = 43007.1
q = 0.8333
p=0.4
c=5.0
vfactor = sqrt (G * mfactor/rfactor)

file=open(inputfile)
lines=file.readlines()
file.close()

Npart = len(lines) - 1

P = pysubs.ParticleData(Npart,0)

n = 0
for line in lines:
	if line.strip().split()[0] == 'x': continue # skip first line
	X=float(line.strip().split()[0])
	Y=float(line.strip().split()[1])
	Z=float(line.strip().split()[2])
	m = sqrt(X**2 + (Y/q)**2 + (Z/p)**2)
	if m > c: continue

	P.x[n]=X * rfactor
	P.y[n]=Y * rfactor
	P.z[n]=Z * rfactor
	P.vx[n]=float(line.strip().split()[3]) * vfactor
	P.vy[n]=float(line.strip().split()[4]) * vfactor
	P.vz[n]=float(line.strip().split()[5]) * vfactor
	P.m[n]=float(line.strip().split()[6]) * mfactor
	P.id[n] = n
	n = n + 1

header = pysubs.CreateGadgetHeader(n,0)
data = pysubs.CreateGadgetData(P,header,1)
class ICfile:
    pass

icfile = ICfile()
icfile.header = header
icfile.data = data
pysubs.WriteGadgetICFile(icfile,outputfile)

#************END MAIN PROGRAM*************************

