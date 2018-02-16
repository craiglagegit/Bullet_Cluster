#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 29-Feb-12

'''
This code reads the xray lookup tables and puts the data in a dictionary.
Then it interpolates xray flux as a function of T and Z.
'''
from pylab import *
import lageconfig
import time
from scipy.special import exp1
from yt.mods import *
#****************************SUBROUTINES**********************************

def ReadLookups(Z):
	# This subroutine reads in the data from the APEC lookup tables and places it in a nested dictionary
	# The data is interpolated to get the data for the required Z (metallicity)
	global ApecData
	ApecData = {} # Dictionary for holding the data in the format:
		  # {(100*logT,(0.5-2.0)):Flux1,(100*logT,(2.0-5.0)):Flux2,(100*logT,(5.0-8.0)):Flux3,(100*logT,(0.5-8.0)):Flux4}
	for bin in [(0.5,2.0),(2.0,5.0),(5.0,8.0)]:
		infile = open(lageconfig.toppath+'bullet/data/apec/apec_xray_'+str(bin[0])+'_'+str(bin[1])+'.txt','r')
		lines = infile.readlines()
		infile.close()
		for line in lines:
			if line.strip().split()[0] == 'LogT': # Skips header line
				continue
			logT = int(round(100*float(line.strip().split()[0])))
			print float(line.strip().split()[0]),logT
			minZ = max(0, int(round((Z * 10))))
			maxZ = minZ + 1
			if maxZ > 10:
				maxZ = 10
				minZ = 9
			f = (maxZ-10.0*Z) * float(line.strip().split()[minZ+1]) + (10.0*Z-minZ) * float(line.strip().split()[maxZ+1])
			ApecData[(logT,bin)] = f
			if bin == (0.5,2.0):
				ApecData[(logT,(0.5,8.0))] = f # This bin is just the sum of the other three
			else:
				ApecData[(logT,(0.5,8.0))] = ApecData[(logT,(0.5,8.0))] + f	
	return 

def ApecXRay(T, bin):
	#This interpolates the Xray flux in terms of T (in kev) and Z for a given bin
	global ApecData
	logT = log10(T)
	minT = max(-100, int(round((logT * 100))))
	maxT = minT + 1
	if maxT > 180:
		maxT = 180
		minT = 179
	flux = (maxT-100.0*logT) * ApecData[(minT,bin)] + (100.0*logT-minT) * ApecData[(maxT,bin)]
	return flux

def OldXRay(T, bin):
	flux = 100*2000/3.745 * (exp1(.37) - exp1(2.7))
	return flux

#****************************END SUBROUTINES**********************************

#****************************MAIN PROGRAM**********************************
global ApecData
Z=0.47
ReadLookups(Z)
T=zeros([20,20,20])

for Tel in np.nditer(T,op_flags=['readwrite']):
	Tel[...] = 1.0+50.0*rand()

start = time.time()
for Tel in np.nditer(T,op_flags=['readwrite']):
	Tel[...] = ApecXRay(Tel,(0.5,8.0))

stop = time.time()
print 'Elapsed time = %.3f\n'%(stop-start)
start = time.time()
for i in range(20):
	for j in range(20):
		for k in range(20):
			T[i,j,k] =  ApecXRay(T[i,j,k],(0.5,8.0))
stop = time.time()
print 'Elapsed time = %.3f\n'%(stop-start)
start = time.time()
for i in range(8000):
	T = 10.0#1.0+50.0*rand()
	flux = OldXRay(T,(0.5,8.0))
	#print 'T=%.2f, Z=%.2f, Flux=%.4f\n'%(T,Z,flux)
stop = time.time()
print 'Elapsed time = %.3f\n'%(stop-start)

#****************************END MAIN PROGRAM**********************************


