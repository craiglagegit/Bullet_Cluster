
#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 26-Nov-12

# Finds DM Histograms

import lageconfig # system specific path information
import matplotlib
matplotlib.use("Agg")
import sys
sys.path.append(lageconfig.bulletpath)
import pysubs_nas_11Nov13 as pysubs

def GetFomParams(filename):
    file = open(filename,'r')
    lines = file.readlines()
    file.close()

    Z = float(lines[0].strip().split()[1])
    TFudge = float(lines[8].strip().split()[1])
    return [Z,TFudge]

def GetAngles(filename):
    file = open(filename,'r')
    line = file.readline()
    file.close()
    phi = float(line.split()[11].strip(","))
    theta = float(line.split()[14].strip(","))
    psi = float(line.split()[17].strip(","))
    return [phi,theta,psi]

#****************MAIN PROGRAM*****************
cmd=sys.argv

snapmin=int(cmd[1])
snapmax=int(cmd[2])

[Z,TFudge] = GetFomParams('fominput')
[phi,theta,psi] = GetAngles('newfom_massx1.out')
pysubs.DMHist(snapmin, snapmax, phi, theta, psi)
  
#************END MAIN PROGRAM*************************
