#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 20-Jan-11


#This program is for reading and plotting the FITS files with bullet cluster SZE data.


from pylab import *
import pyfits

#****************MAIN PROGRAM*****************
sze=pyfits.open("bullet_transfer.fits")
data=sze[0].data
szemax=data.min()
szesum=data.sum()
print "SZE max = %f, SZE sum = %f"%(szemax,szesum)

subplot("221")
contourf(data)
subplot("222")
plot(data[:,25])
show()


#************END MAIN PROGRAM*************************

