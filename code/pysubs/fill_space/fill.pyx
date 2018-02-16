#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 8-Dec-12


#Cython subroutines for fill_space.py

import time
import numpy as np
cimport numpy as np

cdef double phi(np.ndarray[double, ndim=2] Betas, int i, int j, int NumParam):
	cdef double rsquared = 0.0
	cdef int k
	for k in range(NumParam):
		rsquared = rsquared + (Betas[i,k] - Betas[j,k]) * (Betas[i,k] - Betas[j,k])
	return rsquared

def WeedPoints(np.ndarray[double, ndim=2] Betas, int NumPoints, int Kappa, int NumParam):
	cdef int i, j, k, NBetas, WorstPoint, NBetasInit
	cdef double Distij, MinDistance
	NBetasInit = NumPoints * Kappa
	NBetas = NBetasInit
	while NBetas > NumPoints:
		WorstPoint = 0
		MinDistance = 1.0E12
		for i in range(NBetasInit):
			if Betas[i,0] < 0.0:
				continue
			for j in range(i+1,NBetasInit):
				Distij = phi(Betas,i,j,NumParam)
				if Distij < MinDistance:
					MinDistance = Distij
					WorstPoint = i
		Betas[WorstPoint,0] = -1.0
		NBetas = NBetas - 1
		print "NumPoints = %d, MinDistance = %.9f"%(NBetas,MinDistance)
	cdef np.ndarray[double, ndim = 2] NewBetas
	NewBetas = np.zeros([NumPoints,NumParam])
	j = 0
	for i in range(NBetasInit):
		if Betas[i,0] < 0.0:
			continue
		else:
			for k in range(NumParam):
				NewBetas[j,k] = Betas[i,k]
			j = j + 1
	return NewBetas
