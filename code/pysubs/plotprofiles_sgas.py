#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 1-Jun-11


#This program parses the snapshot files from Gadget2 and plots out a 3d plot

import sys
import matplotlib
matplotlib.use("PDF")
from matplotlib.backends.backend_pdf import PdfPages
from pylab import *

def ReadProfiles(filename):
	file = open(filename,'r')
	lines = file.readlines()
	file.close
	R = list()
	RhoDM = list()
	RhoGas = list()
	MDM = list()
	MGas = list()
	TGas = list()
	FG = list()
	for line in lines:
		if line.split()[0] == "#":
			continue
		R.append(float(line.split()[0]))
		RhoDM.append(float(line.split()[1]))
		RhoGas.append(float(line.split()[2]))
		MDM.append(float(line.split()[3]))
		MGas.append(float(line.split()[4]))
		TGas.append(float(line.split()[5]))
		FG.append(float(line.split()[6]))

	return [R,RhoDM,RhoGas,MDM,MGas,TGas,FG]

def ReadParameters(filename):
	file = open(filename,'r')
	lines = file.readlines()
	file.close
	P = float(lines[1].split()[2])
	Q = float(lines[2].split()[2])
	M200 = float(lines[4].split()[2])
	C = float(lines[5].split()[2])
	GF = float(lines[11].split()[2])
	Alpha = float(lines[13].split()[2])
	Rs = float(lines[17].split()[2])
	Beta = float(lines[14].split()[2])
	RC = float(lines[15].split()[2])
	Beta2 = float(lines[18].split()[2])
	RC2 = float(lines[19].split()[2])

	return [P,Q,M200,C,GF,Alpha,Rs,Beta,RC,Beta2,RC2]

#****************MAIN PROGRAM*****************

pp=PdfPages('Graph_Profiles.pdf')

for cluster in ["main","bullet"]:
	figure()
	subplots_adjust(hspace=0.8, wspace=0.4)
	[R,RhoDM,RhoGas,MDM,MGas,TGas,FG] = ReadProfiles(cluster+"profile.dat")
	[P,Q,M200,C,GF,Alpha,Rs,Beta,RC,Beta2,RC2] = ReadParameters("triaxial"+list(cluster)[0]+".cfg")
	subplot(3,1,1)
	title('Log Density (g/cm^3)',fontsize=10)	
	RNorm = array(R) / C # This is R / R200
	plot(log10(RNorm),log10(RhoDM), label = "RhoDM")
	plot(log10(RNorm),log10(RhoGas), label = "RhoGas")
	xlim(-3.0,0.5)
	ylim(-29.0, -23.0)
	xlabel ("log(R / R200)")
	legend(loc=1,bbox_to_anchor=[1.0,1.2])
	subplot(3,1,2)
	title('Log Density (g/cm^3)',fontsize=10)	
	plot(RNorm,log10(RhoDM), label = "RhoDM")
	plot(RNorm,log10(RhoGas), label = "RhoGas")
	xlim(0.0,3.0)
	ylim(-29.0, -23.0)
	xlabel ("R / R200")
	legend(loc=1,bbox_to_anchor=[1.0,1.2])
	#gca().set_xticks([-plotscale,0.0,plotscale])
	#gca().set_yticks([-1.0,0.0,1.0,2.0])

	subplot(3,1,3)
	title('Gas Temp(keV)',fontsize=10)	
	plot(RNorm,TGas)
	xlim(0.0,3.0)
	ylim(0.0,60.0)
	#gca().set_xticks([-plotscale,0.0,plotscale])
	#gca().set_yticks([0.0,10.0,20.0,30.0,40.0])
	xlabel ("R / R200")
	pp.savefig()
figure()
title('Mass Inside R',fontsize=10)	
for cluster in ["main","bullet"]:
	[R,RhoDM,RhoGas,MDM,MGas,TGas,FG] = ReadProfiles(cluster+"profile.dat")
	[P,Q,M200,C,GF,Alpha,Rs,Beta,RC,Beta2,RC2] = ReadParameters("triaxial"+list(cluster)[0]+".cfg")
	RNorm = array(R) / C # This is R / R200
	MDMNorm = array(MDM) / (M200 * 1.93E43)
	MGasNorm = array(MGas) / (M200 * 1.93E43)
	plot(RNorm,MDMNorm, label = "MDM_"+cluster)
	plot(RNorm,MGasNorm, label = "MGas_"+cluster)
xlim(0.0,3.0)
ylim(0.0,4.0)
xlabel ("R / R200")
legend(loc=1,bbox_to_anchor=[1.0,0.8])
pp.savefig()
figure()
title('Log Mass Inside R',fontsize=10)	
for cluster in ["main","bullet"]:
	[R,RhoDM,RhoGas,MDM,MGas,TGas,FG] = ReadProfiles(cluster+"profile.dat")
	[P,Q,M200,C,GF,Alpha,Rs,Beta,RC,Beta2,RC2] = ReadParameters("triaxial"+list(cluster)[0]+".cfg")
	RNorm = array(R) / C # This is R / R200
	MDMNorm = array(MDM) / (M200 * 1.93E43)
	MGasNorm = array(MGas) / (M200 * 1.93E43)
	plot(RNorm,log10(MDMNorm), label = "MDM_"+cluster)
	plot(RNorm,log10(MGasNorm), label = "MGas_"+cluster)
xlim(0.0,3.0)
ylim(-3.0,1.0)
xlabel ("R / R200")
legend(loc=1,bbox_to_anchor=[1.0,0.6])
pp.savefig()

figure()
title('Mass Ratio Inside R',fontsize=10)	
for cluster in ["main","bullet"]:
	[R,RhoDM,RhoGas,MDM,MGas,TGas,FG] = ReadProfiles(cluster+"profile.dat")
	[P,Q,M200,C,GF,Alpha,Rs,Beta,RC,Beta2,RC2] = ReadParameters("triaxial"+list(cluster)[0]+".cfg")
	RNorm = array(R) / C # This is R / R200
	MDMNorm = array(MDM) / (M200 * 1.93E43)
	MGasNorm = array(MGas) / (M200 * 1.93E43)
	plot(RNorm,MGasNorm / MDMNorm, label = "Ratio_"+cluster)
xlim(0.0,3.0)
ylim(0.0,1.0)
xlabel ("R / R200")
legend(loc=1,bbox_to_anchor=[1.0,0.8])
pp.savefig()

figure()
title('Gravitational Force',fontsize=10)	
for cluster in ["main","bullet"]:
	[R,RhoDM,RhoGas,MDM,MGas,TGas,FG] = ReadProfiles(cluster+"profile.dat")
	[P,Q,M200,C,GF,Alpha,Rs,Beta,RC,Beta2,RC2] = ReadParameters("triaxial"+list(cluster)[0]+".cfg")
	RNorm = array(R) / C # This is R / R200
	plot(RNorm,FG, label = "FGrav_"+cluster)
xlim(0.0,3.0)
#ylim(0.0,1.0)
xlabel ("R / R200")
legend(loc=1,bbox_to_anchor=[1.0,0.8])
pp.savefig()

#legend(loc=1,bbox_to_anchor=[-0.2, -0.4])
#ltext = gca().get_legend().get_texts()
#setp(ltext, fontsize = 10)
plt.clf()
pp.close()

#****************END MAIN PROGRAM*****************
