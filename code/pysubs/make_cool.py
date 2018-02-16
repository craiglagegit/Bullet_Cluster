#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 2-Sep-11

#This program generates an interpolated cooling file, given an input metallicity

import lageconfig # system specific path information
import sys
from pylab import *

#****************MAIN PROGRAM*****************

cmd = sys.argv
Z = float(cmd[1])

if Z < 0.0 or Z > 1.0:
	print 'Z out of range, Z = %.3f\n'%Z
	sys.exit()

file = open('/home1/clage/Research/bullet/code/pysubs/apec_cooling.in','r')
lines=file.readlines()
file.close()
NumTemps = len(lines)-2

EF=zeros([11])

LowerColumn = int(Z/0.1)
UpperColumn = LowerColumn + 1
InterpolatingFactor =  (Z / 0.1 - LowerColumn )

file=open('cool_rates.in','w')
file.write('#  ################################Cooling Rate in erg/cm^3/sec, calculated using APEC##################################\n')
header = '#  logT(K)   A = %.3f   A = %.3f   \n'%(Z,Z)
file.write(header)

for i in range(NumTemps):
	try:
		if lines[i].split()[0]=='#': continue
	except: continue
	logT=float(lines[i].split()[0])
	line = '   %.5f  '%logT
	for j in range(11):
		Z = j/10.0
		EF[j] = float(lines[i].split()[j+1])

	EFout = EF[LowerColumn] + InterpolatingFactor * (EF[UpperColumn] - EF[LowerColumn])
	line = line + '%.5f  '%EFout
	line = line + '%.5f  '%EFout
	line = line + '\n'
	file.write(line)
file.close()

#************END MAIN PROGRAM*************************

