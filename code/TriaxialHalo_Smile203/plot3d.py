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
#****************MAIN PROGRAM*****************

scale = 5000.0
#filenames = ['snapshots/SHspherical_000','snapshots/SHspherical_005','snapshots/SHspherical_010']
#filenames = ['snapshots/spherical1_000','snapshots/spherical1_009','snapshots/spherical1_018']
#filenames = ['test2/snapshots/test2_050']#,'test1/snapshots/test1_002', 'test1/snapshots/test1_004']
#filenames = ['datfiles/triaxial11.dat']#,'datfiles/triaxial10.dat']
filenames=['test2/datfiles/test1.dat']
for filename in filenames:
	snapfile = pysubs.ReadGadgetICFile(filename)

	fig=figure()
	ax=Axes3D(fig,aspect=1)

	pysubs.ScatterPlot3d(snapfile,ax,1,0.4,'r')
	pysubs.ScatterPlot3d(snapfile,ax,0,0.4,'g')

	ax.set_xlim3d(-scale,scale)
	ax.set_ylim3d(-scale,scale)
	ax.set_zlim3d(-scale,scale)
	ax.view_init(elev=0, azim=0)

	fig=figure()
	ax=Axes3D(fig,aspect=1)

	pysubs.ScatterPlot3d(snapfile,ax,1,0.4,'r')
	pysubs.ScatterPlot3d(snapfile,ax,0,0.4,'g')

	ax.set_xlim3d(-scale,scale)
	ax.set_ylim3d(-scale,scale)
	ax.set_zlim3d(-scale,scale)
	ax.view_init(elev=0, azim=90)

	fig=figure()
	ax=Axes3D(fig,aspect=1)

	pysubs.ScatterPlot3d(snapfile,ax,1,0.4,'r')
	pysubs.ScatterPlot3d(snapfile,ax,0,0.4,'g')

	ax.set_xlim3d(-scale,scale)
	ax.set_ylim3d(-scale,scale)
	ax.set_zlim3d(-scale,scale)
	ax.view_init(elev=90, azim=90)


show()



#****************END MAIN PROGRAM*****************
