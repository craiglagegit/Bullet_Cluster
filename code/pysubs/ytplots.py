#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 5-Mar-11


#This program makes a series of plots using yt.

from yt.mods import * # set up our namespace
import numpy as np
import sys
from subprocess import *

def _BMag(field, data):
    return np.sqrt(data["Bx"] * data["Bx"] + data["By"] * data["By"] + data["Bz"] * data["Bz"])

def _Beta(field, data):
    mu = 1.673E-24 / 2.0
    kb = 1.38065E-16
    PGas = data["Density"] * data["Temperature"] * kb / mu
    PB = data["BMag"] * data["BMag"] / (8 * np.pi)
    return PB / PGas

#****************MAIN PROGRAM*****************

cmd =sys.argv


snapmin = int(cmd[1])
snapmax = int(cmd[2])

for i in range(snapmin,snapmax):
        if i<10:
            filename = "DD000"+str(i)+"/output_000"+str(i) # parameter file to load
        elif i<100:
            filename = "DD00"+str(i)+"/output_00"+str(i) # parameter file to load
        elif i<1000:
            filename = "DD0"+str(i)+"/output_0"+str(i) # parameter file to load
        else:
            filename = "DD"+str(i)+"/output_"+str(i) # parameter file to load

	pf = load(filename) # load data
	pc = PlotCollection(pf,center=[0.0,0.0,0.0])
	add_field("BMag", function=_BMag)
	add_field("Beta", function=_Beta)

	x=pc.add_slice("Density", 0) # 0 = x-axis
	y=pc.add_slice("Density", 1) # 1 = y-axis
	y.modify["grids"]()
	z=pc.add_slice("Density", 2) # 2 = z-axis
	pc.add_slice("BMag", 0) # 0 = x-axis - By
	pc.add_slice("BMag", 1) # 1 = y-axis - By
	pc.add_slice("BMag", 2) # 2 = z-axis - By
	pc.add_slice("B_Div", 0) # 0 = x-axis - By
	pc.add_slice("B_Div", 1) # 1 = y-axis - By
	pc.add_slice("B_Div", 2) # 2 = z-axis - By
	pc.add_slice("Beta", 0) # 0 = x-axis - By
	pc.add_slice("Beta", 1) # 1 = y-axis - By
	pc.add_slice("Beta", 2) # 2 = z-axis - By
	pc.add_slice("Temperature", 0) # 0 = x-axis
	pc.add_slice("Temperature", 1) # 1 = y-axis
	pc.add_slice("Temperature", 2) # 2 = z-axis
	pc.save(filename) # save all plots
	pc = PlotCollection(pf,center=[3000.0,0.0,0.0])
	pc.add_projection("Dark_Matter_Density", 0) # 0 = x-axis
	pc.add_projection("Dark_Matter_Density", 1) # 1 = y-axis
	pc.add_projection("Dark_Matter_Density", 2) # 2 = z-axis
	pc.save(filename) # save all plots

Popen('touch PlotFinished',shell=True)


#************END MAIN PROGRAM*************************

