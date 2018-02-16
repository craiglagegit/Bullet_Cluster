#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 26-Jan-11


#This program tests convolution of images.

import lageconfig # system specific path information
import sys
sys.path.append(lageconfig.bulletpath)
import pysubs

from pylab import *
import pyfits
from scipy import ndimage

#****************MAIN PROGRAM*****************
k=pyfits.open("bullet_map.fits")
kernel=k[0].data
oldpixels=39
OriginalData=pysubs.Array2d(0.0,oldpixels*10.0,oldpixels,0.0,oldpixels*10.0,oldpixels)

print OriginalData.x[0],OriginalData.x[1],OriginalData.x[2],OriginalData.x[oldpixels-1]

newpixels = int(round(oldpixels*10.0/3.55179888))

for i in range(oldpixels):
	for j in range(oldpixels):
		offset = 221
		OriginalData.data[i,j]=kernel[offset+j,offset+oldpixels-i]

NewData=pysubs.Array2d(0.0,oldpixels*10.0,newpixels,0.0,oldpixels*10.0,newpixels)
print NewData.x[0],NewData.x[1],NewData.x[2],NewData.x[newpixels-1]
for i in range(3,newpixels-3):
	for j in range(3,newpixels-3):
		NewData.data[i,j]=pysubs.DataInterpolate(OriginalData,NewData.x[i],NewData.y[j])

print oldpixels, newpixels, OriginalData.data.sum(), NewData.data.sum()

hdu=pyfits.PrimaryHDU(NewData.data)
hdulist=pyfits.HDUList([hdu])
hdulist.writeto('bullet_map_rescaled.fits')

subplot("221")
contourf(kernel)
subplot("222")
contourf(OriginalData.data)
subplot("223")
contourf(NewData.data)
show()


#************END MAIN PROGRAM*************************

