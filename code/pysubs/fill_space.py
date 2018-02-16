#Author: Craig Lage, NYU; 
#Date: 27-Jan-12

#This code generates maximally spaced points to fill out the remaining space of a simulation region

import lageconfig # system specific path information
from pylab import *
import cython
from fill import WeedPoints

#****************MAIN PROGRAM*****************
# mhd_lowresn29 run32
# Counter=32 C1=2.2700 C2=8.1200 M1=171110.0000 M2=25085.0000 RC1=67.8400 RC2=47.1900 Beta1=0.5249 Beta2=0.5000 Alpha1=0.0496 Alpha2=0.2757 Epsilon1=0.0000 Epsilon2=0.0000 Rs1=1.0000 Rs2=1.0000 Beta21=1.0000 Beta22=1.0000 RC21=1.0000 RC22=1.0000 N21=0.0000 N22=0.0000 GF1=0.1200 GF2=0.1175 Z=0.9000 Mag=6.9500 FieldMode=RandomRadial  Tbg=4.0000 LogRhoMin=-8.0000 IP=158.6700 IPTheta=81.2900 Vinc=0.9985 P1=0.4194 Q1=0.7342 P2=1.0000 Q2=1.0000 Phi1=160.9420 Theta1=40.4725 Psi1=200.7032 Phi2=0.0000 Theta2=0.0000 Psi2=0.0000 Visc=0.1400


scale=[2.27, 8.12, 171110.0, 25085.0, 67.84, 47.19, 0.5249, 0.50, 0.0496, 0.2757, 0.00, 0.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.00, 0.00, 0.12, 0.1175, 0.90, 6.95, 4.0, -8.0, 158.67, 81.29, 0.9985, 0.419, 0.734, 1.0, 1.0, 160.9,40.4,200.7,180.0,180.0,180.0,0.14]

names=['C1','C2', 'M1', 'M2', 'RC1','RC2','Beta1', 'Beta2', 'Alpha1', 'Alpha2', 'Epsilon1', 'Epsilon2', 'Rs1', 'Rs2', 'Beta21', 'Beta22', 'RC21', 'RC22', 'N21', 'N22', 'GF1','GF2','Z','Mag',  'Tbg', 'LogRhoMin', 'IP', 'IPTheta', 'Vinc', 'P1', 'Q1', 'P2', 'Q2', 'Phi1', 'Theta1', 'Psi1', 'Phi2', 'Theta2', 'Psi2', 'Visc']

factors=[1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00]

factormin=[1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,0.90,0.90,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,0.00,0.00,0.00,1.00]
factormax=[1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.20,1.00,1.20,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.10,1.10,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,0.00,0.00,0.00,1.00]


filename = 'points.out'
NumPoints  	= 96	# Number of desired points
NumParam	= len(scale)	# Number of Parameters
Kappa	 	= 5	# Multiplier

Betas 	= zeros([NumPoints*Kappa,NumParam])
NewBetas 	= zeros([NumPoints,NumParam])

print "Got started"

# Create a list of random points, Kappa times larger than needed
for i in range(NumPoints*Kappa):
	for j in range(NumParam):
		Betas[i,j] = factormin[j] + random() * (factormax[j] - factormin[j])

print "Points created"
NewBetas = WeedPoints(Betas, NumPoints, Kappa, NumParam)
# WeedPoints removes the point closest to any of the others, iterating until we have the desired number of points.

file = open(filename,'w')
for i in range(NumPoints):
	outstring = 'Counter='+str(i+1)
	for j in range(NumParam):
		value=Betas[i,j] * scale[j]
		outstring = outstring+' '+names[j]+'=%.4f'%value
		if j == 23:
			outstring = outstring+' FieldMode=RandomRadial '
	outstring = outstring + '\n'
	file.write(outstring)
file.close()

		
#************END MAIN PROGRAM*************************

