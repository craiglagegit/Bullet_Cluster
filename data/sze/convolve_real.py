#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 26-Jan-11


#This program tests convolution of images.


from pylab import *
import pyfits
from scipy import ndimage

#****************MAIN PROGRAM*****************
sze=pyfits.open("bullet_center30.fits")
data=sze[0].data
szemax=data.min()
szesum=data.sum()
print "SZE max = %f, SZE sum = %f"%(szemax,szesum)

k=pyfits.open("bullet_transfer.fits")
kernel=k[0].data

kernelsum=kernel.sum()
print "Kernel sum = %f"%kernelsum

SZlevels=linspace(0,-800E-6,9)

newdata=ndimage.convolve(data,kernel,mode='constant',cval=0.0)

newmax=newdata.min()
newsum=newdata.sum()
#newdata=newdata*szesum/newsum

newmax=newdata.min()
newsum=newdata.sum()

print "New SZE max = %f, New SZE sum = %f"%(newmax,newsum)

#hdu=pyfits.PrimaryHDU(newdata)
#hdulist=pyfits.HDUList([hdu])
#hdulist.writeto('newdata.fits')

subplot("221")
contourf(data,SZlevels)
subplot("222")
contourf(kernel)
subplot("223")
contourf(newdata,SZlevels)
show()


#************END MAIN PROGRAM*************************

