
#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: May 22, 2013

# Testing the data interpolation for shifts and rotations
# Errors always less than 0.5%

import lageconfig # system specific path information
import sys
sys.path.append(lageconfig.bulletpath)
import BulletConstants # Constants used in Bullet cluster simulations
import pysubs_nas_20May13 as pysubs
from pylab import *
import time
from subprocess import *

#****************MAIN PROGRAM*****************
(dataA,sigmaA,maskA,maskANull,dataB1,sigmaB1,dataB2,sigmaB2,dataB3,sigmaB3,dataC,sigmaC,dataD,sigmaD,maskD,maskDNull,dataE,sigmaE) = pysubs.GetData()
dataA_sum = dataA.data.sum()
[dyyA,dxxA] = meshgrid(dataA.y,dataA.x)# Data grid for plots

# New grid smaller than existing grid
grid_ratio = 1.0/12.0
testarray=pysubs.Array2d(2.0*dataA.xmin,2.0*dataA.xmax,int(round(2*grid_ratio*dataA.nx)),2.0*dataA.ymin,2.0*dataA.ymax,int(round(2*grid_ratio*dataA.ny)))
[syyA,sxxA] = meshgrid(testarray.y,testarray.x) # Sim grid for plots
for (theta,shiftx,shifty) in [(0.0,0.0,0.0),(0.0,3.8,-19.7),(0.78,3.8,-19.7),(1.33,38.7,-21.2)]:
	figure()
	for i in range(testarray.nx):
		for j in range(testarray.ny):
			x = testarray.x[i]
			y = testarray.y[j]
			xp = x*cos(theta) + y*sin(theta) + shiftx
			yp = -x*sin(theta) + y*cos(theta) + shifty
			if xp < dataA.xmin or xp > dataA.xmax or yp < dataA.ymin or yp >dataA.ymax:
				testarray.data[i,j] = 0.0
			else:
				testarray.data[i,j] = pysubs.DataInterpolate(dataA,xp,yp)
	testsum = testarray.data.sum() / grid_ratio**2
	print "Data shape = ",dataA.data.shape, " Test shape = ",testarray.data.shape
	print "Theta = %f, Shiftx = %f, Shifty = %f, Data sum = %f, Test sum = %f, Percent loss = %f percent"%(theta,shiftx,shifty,dataA_sum,testsum, abs((testsum-dataA_sum)/dataA_sum*100))
	subplot(1,2,1,aspect=1)
	contourf(dxxA,dyyA,dataA.data)
	subplot(1,2,2,aspect=1)
	contourf(sxxA,syyA,testarray.data)
show()
"""
# Grids same size - just shifted and rotated.

testarray=pysubs.Array2d(2.0*dataA.xmin,2.0*dataA.xmax,2*dataA.nx,2.0*dataA.ymin,2.0*dataA.ymax,2*dataA.ny)
[syyA,sxxA] = meshgrid(testarray.y,testarray.x) # Sim grid for plots
for (theta,shiftx,shifty) in [(0.0,0.0,0.0),(0.0,3.8,-19.7),(0.78,3.8,-19.7),(1.33,38.7,-21.2)]:
	figure()
	for i in range(testarray.nx):
		for j in range(testarray.ny):
			x = testarray.x[i]
			y = testarray.y[j]
			xp = x*cos(theta) + y*sin(theta) + shiftx
			yp = -x*sin(theta) + y*cos(theta) + shifty
			if xp < dataA.xmin or xp > dataA.xmax or yp < dataA.ymin or yp >dataA.ymax:
				testarray.data[i,j] = 0.0
			else:
				testarray.data[i,j] = pysubs.DataInterpolate(dataA,xp,yp)
	testsum = testarray.data.sum()
	print "Theta = %f, Shiftx = %f, Shifty = %f, Data sum = %f, Test sum = %f, Percent loss = %f percent"%(theta,shiftx,shifty,dataA_sum,testsum, abs((testsum-dataA_sum)/dataA_sum*100))
	subplot(1,2,1,aspect=1)
	contourf(dxxA,dyyA,dataA.data)
	subplot(1,2,2,aspect=1)
	contourf(sxxA,syyA,testarray.data)
show()
"""
#************END MAIN PROGRAM*************************

