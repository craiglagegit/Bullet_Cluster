#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 17-Oct-12

#These are Cython routines to speed up various functions

import time
import numpy as np
cimport numpy as np

from yt.mods import * # set up our namespace

def GetPF(snap):
    # Retrieves a yt - pf file given a snapshot number
    if snap<10:
        filename = "DD000"+str(snap)+"/output_000"+str(snap) 
    elif snap<100:
        filename = "DD00"+str(snap)+"/output_00"+str(snap) 
    elif snap<1000:
        filename = "DD0"+str(snap)+"/output_0"+str(snap) 
    else:
        filename = "DD"+str(snap)+"/output_"+str(snap) 
	
    pf = load(filename)
    return pf

def FindEnzoCentroids(pf): 
    start = time.time()
    cdef double PixelVolume, TempMass, TempCentroid
    cdef int i, n, ParticleCounter, NGrids, NPart, NGridPart, TempNum
    NPart = pf.h.grid_particle_count.sum()
    NGrids = pf.h.num_grids
    cdef np.ndarray[double,ndim=1] mass, x, y, z, thismass, thisx, thisy, thisz

    mass=np.zeros([NPart])
    x=np.zeros([NPart])
    y=np.zeros([NPart])
    z=np.zeros([NPart])
    ParticleCounter = 0
    #elapsed = time.time()-start
    #print "Finished Initializing. Time = "+str(elapsed)

    for i in range(NGrids): # Read in all of the particle masses and positions
        PixelVolume = pf.h.grids[i]['dx'] * pf.h.grids[i]['dy'] * pf.h.grids[i]['dz']
        NGridPart = int(pf.h.grid_particle_count[i])
        thismass = pf.h.grids[i]['particle_mass']
        thisx = pf.h.grids[i]['particle_position_x']
        thisy = pf.h.grids[i]['particle_position_y']
        thisz = pf.h.grids[i]['particle_position_z']
        for j in range(NGridPart):
            mass[ParticleCounter] = thismass[j] * PixelVolume
            x[ParticleCounter] = thisx[j]
            y[ParticleCounter] = thisy[j]
            z[ParticleCounter] = thisz[j]
            ParticleCounter = ParticleCounter + 1
    cdef np.ndarray[double,ndim=1] Masses = np.zeros([2])
    cdef np.ndarray[long,ndim=1] NumPart = np.zeros([2], dtype = long)
    cdef np.ndarray[double,ndim=2] Centroids=np.zeros([2,3])
    #elapsed = time.time()-start
    #print "Finished Loading data in arrays. Time = "+str(elapsed)
    
    # Finding the 2 different masses
    Masses[0] = mass[0]
    while Masses[1] < 1.0E-8:
        ParticleCounter = int(NPart * np.random.rand()) # Keep randomly choosing particles until I have found 2 different masses
        if abs(mass[ParticleCounter] - Masses[0]) > 1E-8:
            Masses[1] = mass[ParticleCounter]
            break

    for n in range(NPart): # Cycle through the particles, finding the two centroids
        if abs(mass[n] - Masses[0]) < 1E-8:	
            NumPart[0] = NumPart[0] + 1
            Centroids[0,0] = Centroids[0,0] + x[n]
            Centroids[0,1] = Centroids[0,1] + y[n]
            Centroids[0,2] = Centroids[0,2] + z[n]
        else:	
            NumPart[1] = NumPart[1] + 1
            Centroids[1,0] = Centroids[1,0] + x[n]
            Centroids[1,1] = Centroids[1,1] + y[n]
            Centroids[1,2] = Centroids[1,2] + z[n]
    for k in range(2): # k denotes the bullet or main particles
        if NumPart[0] > NumPart[1]: # Swap everything to put the bullet particles first
            TempNum = NumPart[0]
            NumPart[0] = NumPart[1]
            NumPart[1] = TempNum
            TempMass = Masses[0]
            Masses[0] = Masses[1]
            Masses[1] = TempMass
            for m in range(3):
                TempCentroid = Centroids[0,m]
                Centroids[0,m] = Centroids[1,m]
                Centroids[1,m] = TempCentroid

        for m in range(3):
                Centroids[k,m] = Centroids[k,m] / float(NumPart[k])
    elapsed = time.time()-start
    print "Elapsed time to locate centroids = "+str(elapsed)
    # Returns the number of particles of each mass, the two masses, and the centroids.
    return [NumPart, Masses, Centroids]


