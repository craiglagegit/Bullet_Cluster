#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 27-Mar-12


#This is a standalone python program to run findfom in batch mode

import lageconfig # system specific path information
import sys
import matplotlib
matplotlib.use("PDF")
sys.path.append(lageconfig.bulletpath)
import pysubs_nas_18Sep12 as pysubs
import time


#****************MAIN PROGRAM*****************

cmd=sys.argv

snap=int(cmd[1])
Z=float(cmd[2])
Phi=float(cmd[3])
Theta=float(cmd[4])
Psi=float(cmd[5])
ConstrainPhi=cmd[6]

start = time.time()
if ConstrainPhi in ['True','true','T','t','1']:
    ConstrainPhi = True
else:
    ConstrainPhi = False

(fommassx1,fommass)=pysubs.SimpleFom(snap,Z,phi=Phi,theta=Theta,psi=Psi,ConstrainPhi=ConstrainPhi)

print "Mass FOM = %f, Mass_X1 FOM = %f\n"%(fommass,fommassx1)



#************END MAIN PROGRAM*************************

