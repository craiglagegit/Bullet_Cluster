#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 1-Jun-11


#This program parses the snapshot files from Gadget2 and plots out a 3d plot

import lageconfig # system specific path information
import sys
#import matplotlib
#matplotlib.use("PDF")
sys.path.append(lageconfig.bulletpath)
import  pysubs_mhd_apec_mass_only_1May12 as pysubs
from subprocess import *
from pylab import *
import BulletConstants
from mpl_toolkits.mplot3d.axes3d import Axes3D
import time

def Scatter3d(filename): 
    # Makes a 3D scatter plot
    file = open(filename,'r')
    lines = file.readlines()
    file.close
    Npart=len(lines)

    x=zeros([Npart])
    y=zeros([Npart])
    z=zeros([Npart])
    j = 0

    for line in lines:

        x[j] = float(line.strip().split()[0])
        y[j] = float(line.strip().split()[1])
        z[j] = float(line.strip().split()[2])
	j = j + 1

    return[x,y,z]

#****************MAIN PROGRAM*****************


#filenames = ['snapshots/SHspherical_000','snapshots/SHspherical_005','snapshots/SHspherical_010']
#filenames = ['snapshots/spherical1_000','snapshots/spherical1_009','snapshots/spherical1_018']
#filenames = ['snapshots/tri6NC_000','snapshots/tri6NC_002', 'snapshots/tri6NC_004']
#filenames = ['datfiles/triaxial11.dat']#,'datfiles/triaxial10.dat']
filename='test6/datfiles/density_NFW.dat'

scale = 2.0

[x,y,z] = Scatter3d(filename)
fig=figure()
ax=Axes3D(fig,aspect=1)
ax.plot(x,y,z,marker='.', markersize=0.1, linewidth=0, c='r')

ax.set_xlim3d(-scale,scale)
ax.set_ylim3d(-scale,scale)
ax.set_zlim3d(-scale,scale)
ax.view_init(elev=0, azim=0)

show()



#****************END MAIN PROGRAM*****************
