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
fig = figure()

ax1 = axes([0.1,0.5,0.4,0.4])
ax1.set_title('Log Density (g/cm^3)',fontsize=12)	
ax1.set_xlim(0.0,3.0)
ax1.set_ylim(-29.0,-23.0)
ax1.set_yticks([-28.0,-27.0,-26.0,-25.0,-24.0])
ax1.set_xticks([])
ax1.yaxis.tick_left()
ax2 = axes([0.1,0.1,0.4,0.4])
ax2.set_title('Gas Temp (keV)',fontsize=12, x=0.5, y=0.85)	
ax2.set_xlim(0.0,3.0)
ax2.set_ylim(0.0,60.0)
ax2.set_xticks([0.0,0.5,1.0,1.5,2.0,2.5])
ax2.set_xlabel("R/R200")
ax2.set_yticks([0.0,20.0,40.0])
ax2.yaxis.tick_left()

ax3 = axes([0.5,0.5,0.4,0.4])
ax3.set_title('Log Density (g/cm^3)',fontsize=12)	
ax3.set_xlim(0.0,3.0)
ax3.set_ylim(-29.0,-23.0)
ax3.set_xticks([])
ax3.set_yticks([-28.0,-27.0,-26.0,-25.0,-24.0])
ax3.yaxis.tick_right()
ax4 = axes([0.5,0.1,0.4,0.4])
ax4.set_title('Gas Temp (keV)',fontsize=12, x=0.5, y=0.85)	
ax4.set_xlim(0.0,3.0)
ax4.set_ylim(0.0,30.0)
ax4.set_xticks([0.0,0.5,1.0,1.5,2.0,2.5])
ax4.set_xlabel("R/R200")
ax4.set_yticks([0.0,10.0,20.0])
ax4.yaxis.tick_right()

for ii in range(2):
	[R,RhoDM,RhoGas,MDM,MGas,TGas,FG] = ReadProfiles(dirs[ii]+"mainprofile.dat")
	[P,Q,M200,C,GF,Alpha,Rs,Beta,RC,Beta2,RC2] = ReadParameters(dirs[ii]+"triaxialm.cfg")
	if (C>1.35 and C<2.5) or C>6.5: # Old Masses
		RFactor = 1.0
		RhoFactor = 1.0
		print "ii = %d, Old Masses, M200 =  %.3f, RFactor = %.3f"%(ii,M200,RFactor)
	else: # New Masses
		RFactor = 1.0 * 1.37
		RhoFactor = 1.0 / (1.37**3.0)
		print "ii = %d, New Masses, M200 = %.3f, RFactor = %.3f"%(ii,M200,RFactor)
	RNorm = array(R) * RFactor / C # This is R / R200
	RhoGas = array(RhoGas) * RhoFactor
	RhoDM = array(RhoDM) * RhoFactor
	if ii == 0:
		linestyle = '-'
	else:
		linestyle='--'
	ax1.plot(RNorm,log10(RhoDM), ls = linestyle, label = "RhoDM_"+names[ii], lw=2, color = 'b', marker = 's', mec = 'b', ms = 4, markevery=10)
	ax1.plot(RNorm,log10(RhoGas), ls = linestyle, label = "RhoGas_"+names[ii], lw=2, color = 'g', marker = 'o', mec = 'g', ms = 4, markevery=10)
	ax2.plot(RNorm,TGas, ls = linestyle,label = "Temp_"+names[ii], lw=2, color = 'r')
	[R,RhoDM,RhoGas,MDM,MGas,TGas,FG] = ReadProfiles(dirs[ii]+"bulletprofile.dat")
	[P,Q,M200,C,GF,Alpha,Rs,Beta,RC,Beta2,RC2] = ReadParameters(dirs[ii]+"triaxialb.cfg")
	RNorm = array(R) / C # This is R / R200
   	ax3.plot(RNorm,log10(RhoDM), ls = linestyle, label = "RhoDM_"+names[ii], lw=2, color = 'b', marker = 's', mec = 'b', ms = 4, markevery=10)
	ax3.plot(RNorm,log10(RhoGas), ls = linestyle, label = "RhoGas_"+names[ii], lw=2, color = 'g', marker = 'o', mec = 'g', ms = 4, markevery=10)
	ax4.plot(RNorm,TGas, ls = linestyle,label = "Temp_"+names[ii], lw=2, color = 'r')

ax1.legend(loc=1,bbox_to_anchor=[1.0,1.0])
l1text = ax1.get_legend().get_texts()
setp(l1text, fontsize = 12, color = 'k')
ax2.legend(loc=1,bbox_to_anchor=[1.0,0.8])
l2text = ax2.get_legend().get_texts()
setp(l2text, fontsize = 12, color = 'k')
ax3.legend(loc=1,bbox_to_anchor=[1.0,1.0])
l3text = ax3.get_legend().get_texts()
setp(l3text, fontsize = 12, color = 'k')
ax4.legend(loc=1,bbox_to_anchor=[1.0,0.8])
l4text = ax4.get_legend().get_texts()
setp(l4text, fontsize = 12, color = 'k')

pp.savefig()
plt.clf()
pp.close()

#****************END MAIN PROGRAM*****************
