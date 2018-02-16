
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
import pysubs_nas_19Mar13 as pysubs
from pylab import *
import time
from subprocess import *
import os
import copy

#****************MAIN PROGRAM*****************
filename = 'fominput'
file = open(filename,'r')
lines = file.readlines()
file.close()

Z = float(lines[0].strip().split()[1])
snap = 5
theta = 0.85
psi = 3.00
TFudge = float(lines[8].strip().split()[1])

for i in range(10):
	fommassx1 = pysubs.SimpleFom(snap,Z,phi=0.0,theta=theta,psi=psi,ConstrainPhi=True,TFudge=TFudge)
	print "Run %d, fom = %.4f"%(i,fommassx1)
	sys.stdout.flush()
#************END MAIN PROGRAM*************************

