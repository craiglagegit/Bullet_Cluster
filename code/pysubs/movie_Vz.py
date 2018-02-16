#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 18-May-11


#This program generates a 'publication quality' movie.

import lageconfig # system specific path information
import matplotlib
matplotlib.use("Agg")
import sys
sys.path.append(lageconfig.bulletpath)
import pysubs_nas_11Nov13 as pysubs

def GetFomParams(filename):
    file = open(filename,'r')
    lines = file.readlines()
    file.close()

    Z = float(lines[0].strip().split()[1])
    TFudge = float(lines[8].strip().split()[1])
    return [Z,TFudge]

def GetAngles(filename):
    file = open(filename,'r')
    line = file.readline()
    file.close()
    phi = float(line.split()[11].strip(","))
    theta = float(line.split()[14].strip(","))
    psi = float(line.split()[17].strip(","))
    return [phi,theta,psi]


#****************MAIN PROGRAM*****************
cmd=sys.argv

snapmin=int(cmd[1])
snapmax=int(cmd[2])
#Z=float(cmd[3])
#TFudge=float(cmd[4])
#phi=float(cmd[5])
#theta=float(cmd[6])
#psi=float(cmd[7])

zmax = 3000.0
zmin = -zmax
#xmax = 3000.0
#ymax = 1500.0
xmax = 1500.0
ymax = 750.0
ymin = -ymax
xmin = -xmax
nx = 256
ny = 128
[Z,TFudge] = GetFomParams('fominput')
[phi,theta,psi] = GetAngles('newfom_massx1.out')
xyarray = pysubs.Array2d(xmin,xmax,nx,ymin,ymax,ny)

pysubs.MoviePlotsRMVz(snapmin,snapmax,xyarray,zmin,zmax,Z,phi=phi,theta=theta,psi=psi,TFudge=TFudge)
#************END MAIN PROGRAM*************************

