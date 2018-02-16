#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 2-Sep-11


#This is a standalone python program to run findfom in batch mode

import lageconfig # system specific path information
import sys
import matplotlib
matplotlib.use("PDF")
sys.path.append(lageconfig.bulletpath)
import pysubs_mhd_play as pysubs


#****************MAIN PROGRAM*****************
pysubs.CreateGadgetMarkerParticles('collision.dat','collision_mod.dat',5,400.0,5,400.0)


#************END MAIN PROGRAM*************************

