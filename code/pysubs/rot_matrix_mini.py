#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 9-Sep-11


#This program looks at a matrix of rotation angles on a databse with an IP of 4 kpc

import lageconfig # system specific path information
import sys
sys.path.append(lageconfig.bulletpath)
import pysubs_mhd_new as pysubs
from pylab import *


#****************MAIN PROGRAM*****************
cmd=sys.argv

toppath=lageconfig.toppath
dir='./'
snapmin=int(cmd[1])
snapmax=int(cmd[2])
counter=1
psi=0.0
filename='fomrot.out'

for theta in [2.55,2.65,2.75]:
	for phi in [3.09,3.14,3.19]:

		[besttime,bestfom]=pysubs.FindFom(toppath,dir,snapmin,snapmax,counter,phi,theta,psi,simulator='Enzo',PlotSuffix=None,IncludeTemp=False)


		result='Phi = %.4f, Theta = %.4f, FOM=%.4f at time =  %.4f\n'%(phi,theta,bestfom,besttime)
		fomfile=open(filename,'a')
		fomfile.write(result)
		fomfile.close

#************END MAIN PROGRAM*************************

