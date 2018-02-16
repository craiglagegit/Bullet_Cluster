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


def FindIDData(snapfile,parttype,id,filetype): 
    # Prints selected data
    # filetype=0->snapshot file, filetype=1->IC file

    x = zeros([1])
    y = zeros([1])
    z = zeros([1])
    offset=0
    for k in range(parttype):
        offset=offset+snapfile.header.Npart[k]

    for i in range(snapfile.header.Npart[parttype]):
        io = i+offset
	if id == snapfile.data.id[io]:
		x[0] = snapfile.data.pos[3*io]
		y[0] = snapfile.data.pos[3*io+1]
		z[0] = snapfile.data.pos[3*io+2]
		break
    return [x,y,z]




#****************MAIN PROGRAM*****************


filenamebase = 'test2/snapshots/test2'
id = 334567
parttype = 1
scale = 500.0

fig=figure()
ax=Axes3D(fig,aspect=1)


for snap in  range(50):#[0,10,20,30,40,50]:
	if snap<10:
	    filename = filenamebase+"_00"+str(snap) 
	elif snap<100:
	    filename = filenamebase+"_0"+str(snap) 
	else:
	    filename = filenamebase+"_"+str(snap) 

	snapfile = pysubs.ReadGadgetSnapshot(filename)
	[x,y,z] = FindIDData(snapfile,1,id,0)
    	ax.plot(x,y,z,marker='.', markersize=10.0, linewidth=1.0, c='r')

ax.set_xlim3d(-scale,scale)
ax.set_ylim3d(-scale,scale)
ax.set_zlim3d(-scale,scale)
ax.view_init(elev=0, azim=0)

show()


#****************END MAIN PROGRAM*****************
