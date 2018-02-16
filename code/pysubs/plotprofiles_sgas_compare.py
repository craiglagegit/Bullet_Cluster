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

cmd =sys.argv
dirs = []
names=[]
dirs.append(cmd[1]+'/')
names.append(cmd[2])
dirs.append(cmd[3]+'/')
names.append(cmd[4])


pp=PdfPages('profiles/Graph_Profiles_'+names[0]+'_'+names[1]+'.pdf')

for cluster in ["main","bullet"]:
	figure()
	suptitle(cluster+" profiles")
	subplots_adjust(hspace=0.8, wspace=0.4)
	for ii in range(2):
		[R,RhoDM,RhoGas,MDM,MGas,TGas,FG] = ReadProfiles(dirs[ii]+cluster+"profile.dat")
		[P,Q,M200,C,GF,Alpha,Rs,Beta,RC,Beta2,RC2] = ReadParameters(dirs[ii]+"triaxial"+list(cluster)[0]+".cfg")

		ax1 = axes([0.1,0.7,0.6,0.2])
		ax1.set_title('Log Density (g/cm^3)',fontsize=10)	
		if (C>1.35 and C<2.5) or C>6.5: # Old Masses
			RFactor = 1.0
			RhoFactor = 1.0
		else: # New Masses
			RFactor = 1.0 * 1.37
			RhoFactor = 1.0 / (1.37**3.0)
		RNorm = array(R) * RFactor / C # This is R / R200
		RhoGas = array(RhoGas) * RhoFactor
		RhoDM = array(RhoDM) * RhoFactor
		ax1.plot(log10(RNorm),log10(RhoDM), label = "RhoDM"+names[ii])
		ax1.plot(log10(RNorm),log10(RhoGas), label = "RhoGas"+names[ii])
		ax1.set_xlim(-3.0,0.5)
		ax1.set_ylim(-29.0, -23.0)
		ax1.set_xlabel ("log(R / R200)")
		ax1.legend(loc=1,bbox_to_anchor=[1.5,1.2])
		ltext1 = ax1.get_legend().get_texts()
		setp(ltext1, fontsize = 12)

		ax2 = axes([0.1,0.4,0.6,0.2])
		ax2.set_title('Log Density (g/cm^3)',fontsize=10)	
		ax2.plot(RNorm,log10(RhoDM), label = "RhoDM"+names[ii])
		ax2.plot(RNorm,log10(RhoGas), label = "RhoGas"+names[ii])
		ax2.set_xlim(0.0,3.0)
		ax2.set_ylim(-29.0, -23.0)
		ax2.set_xlabel ("R / R200")
		ax2.legend(loc=1,bbox_to_anchor=[1.5,1.2])
		ltext2 = ax2.get_legend().get_texts()
		setp(ltext2, fontsize = 12)
		#gca().set_xticks([-plotscale,0.0,plotscale])
		#gca().set_yticks([-1.0,0.0,1.0,2.0])

		ax3 = axes([0.1,0.1,0.6,0.2])
		ax3.set_title('Gas Temp(keV)',fontsize=10)	
		ax3.plot(RNorm,TGas, label = names[ii])
		ax3.set_xlim(0.0,3.0)
		ax3.set_ylim(0.0,60.0)
		#gca().set_xticks([-plotscale,0.0,plotscale])
		#gca().set_yticks([0.0,10.0,20.0,30.0,40.0])
		ax3.set_xlabel ("R / R200")
		ax3.legend(loc=1,bbox_to_anchor=[1.5,1.2])
		ltext3 = ax3.get_legend().get_texts()
		setp(ltext3, fontsize = 12)

	pp.savefig()
plt.clf()
pp.close()

#****************END MAIN PROGRAM*****************
