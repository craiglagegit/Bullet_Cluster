#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 26-Jan-11


#This program tests convolution of images.


from pylab import *
import pyfits
from scipy import ndimage

#****************MAIN PROGRAM*****************
sze=pyfits.open("bullet_center25.fits")
data=sze[0].data
szemax=data.min()
szesum=data.sum()
print "SZE max = %f, SZE sum = %f"%(szemax,szesum)

k=pyfits.open("bullet_transfer_rescaled.fits")
kernel=k[0].data

newdata=ndimage.convolve(data,kernel,mode='constant',cval=0.0)
subplot("221")
contourf(data)
subplot("222")
contourf(kernel)
subplot("223")
contourf(newdata)
show()


#************END MAIN PROGRAM*************************

