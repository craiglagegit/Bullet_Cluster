#Author: Craig Lage, NYU; 
#Date: 20-Mar-13


#This program centers the densest grid on the cluster centers

import sys
from subprocess import *
from pylab import *

def CalcGridEdges(OuterGridLeftEdge, OuterGridRightEdge, OuterGridN, InnerGridN, InnerGridCenter):
	OuterGridD=[]
	InnerGridLeftEdge=[]
	InnerGridRightEdge=[]
	InnerGridD=[]
	InnerGridSize=[]
	for i in range(3):
		Step = (OuterGridRightEdge[i] - OuterGridLeftEdge[i]) / OuterGridN[i]
		OuterGridD.append(Step)
		InnerGridD.append(Step / 2.0)
		InnerGridSize.append(Step/2.0 * InnerGridN[i])
	for i in range(3):
		InnerGridLeftEdge.append(OuterGridLeftEdge[i] + OuterGridD[i] * round(((InnerGridCenter[i]-InnerGridSize[i]/2.0) - OuterGridLeftEdge[i]) / OuterGridD[i]))
		InnerGridRightEdge.append(InnerGridLeftEdge[i] + InnerGridSize[i])

	return [InnerGridLeftEdge, InnerGridRightEdge]

def ReadClusterCenters(filename):
	file = open(filename,'r')
	lines = file.readlines()
	file.close()
	line1 = lines[0].split()
	line2 = lines[1].split()
	BulletCenter=[float(line1[5]),float(line1[7]),float(line1[9])]
	MainCenter=[float(line2[5]),float(line2[7]),float(line2[9])]
	return [MainCenter,BulletCenter]

def ModAMRTest(filename,BulletLeftEdge,BulletRightEdge,MainLeftEdge,MainRightEdge):
	file = open(filename,'r')
	lines = file.readlines()
	file.close()

	lines[46] = "CosmologySimulationGridLeftEdge[4]      = %.4f %.4f %.4f\n"%(BulletLeftEdge[0],BulletLeftEdge[1],BulletLeftEdge[2])
	lines[47] = "CosmologySimulationGridRightEdge[4]     = %.4f %.4f %.4f\n"%(BulletRightEdge[0],BulletRightEdge[1],BulletRightEdge[2])
	lines[50] = "CosmologySimulationGridLeftEdge[5]      = %.4f %.4f %.4f\n"%(MainLeftEdge[0],MainLeftEdge[1],MainLeftEdge[2])
	lines[51] = "CosmologySimulationGridRightEdge[5]     = %.4f %.4f %.4f\n"%(MainRightEdge[0],MainRightEdge[1],MainRightEdge[2])

	file = open(filename,'w')
	for line in lines:
		file.write(line)
	file.close()
		
	return


#*************************************MAIN PROGRAM*******************************************
collisionfilename = 'collision.out'
AMRfilename = 'AMRTest.enzo'
OuterGridLeftEdge = [-2906.25,-750.0,-750.0]
OuterGridRightEdge = [843.75,750.0,750.0]
OuterGridN = [160,64,64]
InnerGridN = [32,32,32]

[MainCenter,BulletCenter] = ReadClusterCenters(collisionfilename)
[BulletLeftEdge,BulletRightEdge] = CalcGridEdges(OuterGridLeftEdge, OuterGridRightEdge, OuterGridN, InnerGridN, BulletCenter)
[MainLeftEdge,MainRightEdge] = CalcGridEdges(OuterGridLeftEdge, OuterGridRightEdge, OuterGridN, InnerGridN, MainCenter)

#print BulletLeftEdge, BulletRightEdge
#print MainLeftEdge, MainRightEdge

ModAMRTest(AMRfilename,BulletLeftEdge,BulletRightEdge,MainLeftEdge,MainRightEdge)

#*************************************END MAIN PROGRAM*******************************************

		  
