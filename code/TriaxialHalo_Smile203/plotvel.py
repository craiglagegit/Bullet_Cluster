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

def ScatterVel3d(snapfile,ax,parttype,size,col): 
    # Makes a 3D scatter plot
    # Need to clean this up to use Array classes
    Npart=snapfile.header.Npart
    num=Npart[parttype]

    offset=0
    for k in range(parttype):
        offset=offset+Npart[k]

    x=zeros([num])
    y=zeros([num])
    z=zeros([num])

    for j in range(num):

        jo=offset+j
        x[j] = snapfile.data.vel[3*jo]
        y[j] = snapfile.data.vel[3*jo+1]
        z[j] = snapfile.data.vel[3*jo+2]

    ax.plot(x,y,z,marker='.', markersize=size, linewidth=0, c=col)


#****************MAIN PROGRAM*****************


#filenames = ['datfiles/SHspherical.dat', 'datfiles/spherical1.dat','datfiles/spherical2.dat']#['snapshots/SHspherical1_000','snapshots/spherical7_000']#,'snapshots/spherical3_011']
#filenames = ['test1/datfiles/nbody2_100K.dat','datfiles/triaxial11.dat','datfiles/triaxial12.dat','datfiles/triaxial18.dat','datfiles/triaxial_log1.dat','datfiles/triaxial_log4.dat']
filenames=['test2/datfiles/test1_trunc.dat']#,'datfiles/tri7F.dat','datfiles/tri8F.dat','datfiles/tri9F.dat']
for filename in filenames:
	snapfile = pysubs.ReadGadgetICFile(filename)

	fig=figure()
	ax=Axes3D(fig,aspect=1)

	ScatterVel3d(snapfile,ax,1,0.4,'r')
	#ScatterVel3d(snapfile,ax,0,0.4,'g')

	ax.set_xlim3d(-5000.0,5000.0)
	ax.set_ylim3d(-5000.0,5000.0)
	ax.set_zlim3d(-5000.0,5000.0)
	ax.view_init(elev=0, azim=0)

	fig=figure()
	ax=Axes3D(fig,aspect=1)

	ScatterVel3d(snapfile,ax,1,0.4,'r')
	#ScatterVel3d(snapfile,ax,0,0.4,'g')

	ax.set_xlim3d(-5000.0,5000.0)
	ax.set_ylim3d(-5000.0,5000.0)
	ax.set_zlim3d(-5000.0,5000.0)
	ax.view_init(elev=90, azim=0)

	fig=figure()
	ax=Axes3D(fig,aspect=1)

	ScatterVel3d(snapfile,ax,1,0.4,'r')
	#ScatterVel3d(snapfile,ax,0,0.4,'g')

	ax.set_xlim3d(-5000.0,5000.0)
	ax.set_ylim3d(-5000.0,5000.0)
	ax.set_zlim3d(-5000.0,5000.0)
	ax.view_init(elev=0, azim=90)


show()



#****************END MAIN PROGRAM*****************
