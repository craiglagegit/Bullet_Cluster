#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 1-Jun-11


#This program parses the snapshot files from Gadget2 and plots out a 3d plot

import lageconfig # system specific path information
import sys
#import matplotlib
#matplotlib.use("PDF")
sys.path.append(lageconfig.bulletpath)
import  pysubs_mhd_apec_mass_x1_3May12 as pysubs
from subprocess import *
from pylab import *
import BulletConstants
from mpl_toolkits.mplot3d.axes3d import Axes3D
import time


def PrintIDData(snapfile,parttype,id,filetype): 
    # Prints selected data
    # filetype=0->snapshot file, filetype=1->IC file

    offset=0
    for k in range(parttype):
        offset=offset+snapfile.header.Npart[k]

    for i in range(snapfile.header.Npart[parttype]):
        io = i+offset
	if id == snapfile.data.id[io]:
		for j in range(3):
		    print "PartType ",parttype,"Particle ",i," X",j," = ",snapfile.data.pos[3*io+j]
		    print "PartType ",parttype,"Particle ",i," V",j," = ",snapfile.data.vel[3*io+j]
		print "PartType ",parttype,"Particle ",i," id = ",snapfile.data.id[io]
		print "PartType ",parttype,"Particle ",i," mass = ",snapfile.data.masses[io]
		if parttype==0:
		    print "PartType ",parttype,"Particle ",i," u = ",snapfile.data.u[io]
		    if filetype==0:
		    	print "PartType ",parttype,"Particle ",i," rho = ",snapfile.data.rho[io]
		    	print "PartType ",parttype,"Particle ",i," ne = ",snapfile.data.ne[io]
		    	print "PartType ",parttype,"Particle ",i," np = ",snapfile.data.np[io]
		    	print "PartType ",parttype,"Particle ",i," hsml = ",snapfile.data.hsml[io]
		#break
    return




#****************MAIN PROGRAM*****************


filename = 'test2/datfiles/fulltest1.dat'
id = 0

snapfile = pysubs.ReadGadgetICFile(filename)
pysubs.PrintGadgetHeader(snapfile)
pysubs.PrintGadgetData(snapfile,0,0,3,1)
pysubs.PrintGadgetData(snapfile,1,0,3,1)
for i in [0,1,2,3,99998,99999,100000,100001,100002,999999,1000000,1000001,1099998,1099999]:
	print "ID = %d\n"%snapfile.data.id[i]
PrintIDData(snapfile,1,id,1)
PrintIDData(snapfile,0,id,1)


#****************END MAIN PROGRAM*****************
