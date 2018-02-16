#Author: Craig Lage, NYU; 
#Date: 27-Jan-12

#This code generates maximally spaced points to fill out the remaining space of a simulation region

import lageconfig # system specific path information
from pylab import *
import sys
import cython
from fill import WeedPoints

def ReadStartingFile(startingfile,startingcounter):
	file = open(startingfile,'r')
	lines = file.readlines()
	file.close()
	entries = lines[startingcounter-1].split()
	scale = []
	names = []
	factors = []
	factormin = []
	factormax = []
	for entry in entries:
		name = entry.split('=')[0]
		value = entry.split('=')[1]
		if name == 'Counter' or name == 'FieldMode':
			continue
		else:
			names.append(name)
			value = float(value)
			if abs(value) > 1.0E-6:
				scale.append(float(value))
				factors.append(1.0)
				factormin.append(1.0)
				factormax.append(1.0)
			else:
				scale.append(1.0)
				factors.append(0.0)
				factormin.append(0.0)
				factormax.append(0.0)
	return [scale,names,factors,factormin,factormax]



def ReadFactorFile(factorfile,names,factors,factormin,factormax):
	file = open(factorfile,'r')
	lines = file.readlines()
	file.close()
	for line in lines:
		fname = line.split()[0]
		factor = float(line.split()[1])
		fmin = float(line.split()[2])
		fmax = float(line.split()[3])
		for i, name in enumerate(names):
			if name == fname:
				factors[i] = factor
				factormin[i] = fmin
				factormax[i] = fmax
	return [factors,factormin,factormax]

#****************MAIN PROGRAM*****************
cmd=sys.argv 
startingfile = cmd[1]
startingcounter = int(cmd[2])
NumPoints = int(cmd[3])
Kappa = int(cmd[4])
factorfile = cmd[5]
[scale,names,factors,factormin,factormax] = ReadStartingFile(startingfile,startingcounter)
[factors,factormin,factormax] = ReadFactorFile(factorfile,names,factors,factormin,factormax)

filename = 'points.out'
NumParam	= len(scale)	# Number of Parameters

Betas 	= zeros([NumPoints*Kappa,NumParam])
NewBetas 	= zeros([NumPoints,NumParam])

print "Got started"

# Create a list of random points, Kappa times larger than needed
for i in range(NumPoints*Kappa):
	for j in range(NumParam):
		Betas[i,j] = factormin[j] + random() * (factormax[j] - factormin[j])
	for PQPairs in [(29,30),(31,32)]: # This guarantees that P < Q, as Smile requires
		if Betas[i,PQPairs[0]] * scale[PQPairs[0]] > Betas[i,PQPairs[1]] * scale[PQPairs[1]]:
			Betas[i,PQPairs[0]], Betas[i,PQPairs[1]] = Betas[i,PQPairs[1]], Betas[i,PQPairs[0]] 

	for MGFPairs in [(2,20),(3,21)]: # This adjusts GF as 1/M
		Betas[i,MGFPairs[1]] = 1.0 / Betas[i,MGFPairs[0]]
print "Points created"
NewBetas = WeedPoints(Betas, NumPoints, Kappa, NumParam)
# WeedPoints removes the point closest to any of the others, iterating until we have the desired number of points.

file = open(filename,'w')
for i in range(NumPoints):
	outstring = 'Counter='+str(i+1)
	for j in range(NumParam):
		value=NewBetas[i,j] * scale[j]
		outstring = outstring+' '+names[j]+'=%.4f'%value
		if j == 23:
			outstring = outstring+' FieldMode=RandomRadial '
	outstring = outstring + '\n'
	file.write(outstring)
file.close()

		
#************END MAIN PROGRAM*************************

