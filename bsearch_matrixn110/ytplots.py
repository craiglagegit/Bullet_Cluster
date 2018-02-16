#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 5-Mar-11


#This program makes a series of plots using yt.

from yt.mods import * # set up our namespace
import numpy as np
import sys

def _BMag(field, data):
    return np.sqrt(data["Bx"] * data["Bx"] + data["By"] * data["By"] + data["Bz"] * data["Bz"])

def _Beta(field, data):
    mu = 1.673E-24 / 2.0
    kb = 1.38065E-16
    PGas = data["Density"] * data["Temperature"] * kb / mu
    PB = data["BMag"] * data["BMag"] / (8 * np.pi)
    return PB / PGas

def _EnergyRatio(field, data):
    ER = data["GasEnergy"] / data["TotalEnergy"]
    return ER


def _TemperatureRatio(field, data):
    ER = data["TotalEnergy"] / data["GasEnergy"]
    eta = 0.2
    TR = 1.0 - eta * (ER - 1.0)
    print TR.min()
    for i,x in enumerate(TR):
        if x < 0.4:
            TR[i] = 0.4
    print TR.min()
    #sys.exit()
    return TR

#****************MAIN PROGRAM*****************

cmd =sys.argv



dir = 'batchfiles/batch'+cmd[1]+'/ddfiles/run'+cmd[2]+'/'
snapmin = int(cmd[3])
snapmax = int(cmd[4])

for i in range(snapmin,snapmax):
        if i<10:
            filename = dir+"DD000"+str(i)+"/output_000"+str(i) # parameter file to load
        elif i<100:
            filename = dir+"DD00"+str(i)+"/output_00"+str(i) # parameter file to load
        elif i<1000:
            filename = dir+"DD0"+str(i)+"/output_0"+str(i) # parameter file to load
        else:
            filename = dir+"DD"+str(i)+"/output_"+str(i) # parameter file to load

	pf = load(filename) # load data
	pc = PlotCollection(pf,center=[0.0,0.0,0.0])
	add_field("BMag", function=_BMag)
	add_field("Beta", function=_Beta)
        #add_field("TemperatureRatio", function=_TemperatureRatio, take_log=False)

        minrho = pf.h.find_min('Density')
        maxrho = pf.h.find_max('Density')
        print "minimum Density = %.3g at "%minrho[0], minrho[1],"\n"
        print "maximum Density = %.3g at "%maxrho[0], maxrho[1],"\n"

        minbeta = pf.h.find_min('Beta')
        maxbeta = pf.h.find_max('Beta')
        print "minimum Beta = %.3g at "%minbeta[0], minbeta[1],"\n"
        print "maximum Beta = %.3g at "%maxbeta[0], maxbeta[1],"\n"

	x=pc.add_slice("Density", 0) # 0 = x-axis
	y=pc.add_slice("Density", 1) # 1 = y-axis
	y.modify["grids"]()
	z=pc.add_slice("Density", 2) # 2 = z-axis
	pc.add_slice("Dark_Matter_Density", 0) # 0 = x-axis
        pc.add_slice("Dark_Matter_Density", 1) # 0 = y-axis
        pc.add_slice("Dark_Matter_Density", 2) # 0 = z-axis
	#pc.add_slice("B_Div", 0) # 0 = x-axis - By
	#pc.add_slice("B_Div", 1) # 1 = y-axis - By
	#pc.add_slice("B_Div", 2) # 2 = z-axis - By
	pc.add_slice("Beta", 0) # 0 = x-axis - By
	pc.add_slice("Beta", 1) # 1 = y-axis - By
	pc.add_slice("Beta", 2) # 2 = z-axis - By
	pc.add_slice("Temperature", 0) # 0 = x-axis
	pc.add_slice("Temperature", 1) # 1 = y-axis
	pc.add_slice("Temperature", 2) # 2 = z-axis
	#pc.add_slice("TemperatureRatio", 0) # 0 = x-axis
	#pc.add_slice("TemperatureRatio", 1) # 1 = y-axis
	#pc.add_slice("TemperatureRatio", 2) # 2 = z-axis
	pc.save(filename) # save all plots


#************END MAIN PROGRAM*************************

