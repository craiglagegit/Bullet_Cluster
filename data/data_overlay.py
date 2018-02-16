#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 31-May-11


#This program overlays the lensing and xray data from the bullet cluster
#to check the offsets. Adding XRay Temp data.
#


import lageconfig # system specific path information
import sys
sys.path.append(lageconfig.bulletpath)
import pysubs
from pylab import *
from subprocess import *
import array
from scipy.ndimage import gaussian_filter
import BulletConstants
	
#****************MAIN PROGRAM*****************

toppath = lageconfig.toppath

datapathA=toppath+'bullet/data/kappa_mod_2.dat'
datapathB=toppath+'bullet/data/xray_500_6000_No_PS.dat'
datapathC=toppath+'bullet/data/sze_data.dat'
datapathD=toppath+'bullet/data/xray_temp.dat'

datapathsA=toppath+'bullet/data/mass_sigma.dat'
datapathsB=toppath+'bullet/data/xray_sigma.dat'
datapathsC=toppath+'bullet/data/sze_sigma.dat'
datapathszekernel=toppath+'bullet/data/sze/bullet_transfer_rescaled_small.fits'
datapathkernel60=toppath+'bullet/data/kernel60.dat'

# First, the Mass data - designated as A
xscaleA=yscaleA=BulletConstants.AngularScale*3600*BulletConstants.MassDataScale 
# Pixel size in kpc. 4.413 kpc/" is the angular scale at the bullet cluster
# 3600 "/degree, 9.86E-4 degrees/pixel is the data scale

xpixelsA=ypixelsA=110
sxpixelsA=sypixelsA=220
dxmaxA = xscaleA*xpixelsA/2
dymaxA = yscaleA*ypixelsA/2
dxminA = -xscaleA*xpixelsA/2
dyminA = -yscaleA*ypixelsA/2
sxmaxA=xscaleA*xpixelsA
symaxA=yscaleA*ypixelsA

dataA=pysubs.Array2d(dxminA,dxmaxA,xpixelsA,dyminA,dymaxA,ypixelsA)
dataA=pysubs.GetBulletData(datapathA,dataA)# Mass density measured data
sigmaA=pysubs.Array2d(dxminA,dxmaxA,xpixelsA,dyminA,dymaxA,ypixelsA)
sigmaA=pysubs.GetBulletData(datapathsA,sigmaA)#  Mass sigma
[dyyA,dxxA] = meshgrid(dataA.y,dataA.x)# Data grid for plots

xscaleB=yscaleB=BulletConstants.AngularScale*3600*BulletConstants.XRayDataScale 
# Pixel size in kpc. 4.413 kpc/" is the angular scale at the bullet cluster
# 3600 "/degree, 9.86E-4 degrees/pixel is the data scale

xpixelsB=ypixelsB=110
sxpixelsB=sypixelsB=220
dxmaxB = xscaleB*xpixelsB/2
dymaxB = yscaleB*ypixelsB/2
sxmaxB = dxmaxB*2
symaxB = dymaxB*2
dataB=pysubs.Array2d(-dxmaxB,dxmaxB,xpixelsB,-dymaxB,dymaxB,ypixelsB)
dataB=pysubs.GetBulletData(datapathB,dataB)# XRay measured data
sigmaB=pysubs.Array2d(-dxmaxB,dxmaxB,xpixelsB,-dymaxB,dymaxB,ypixelsB)
sigmaB=pysubs.GetBulletData(datapathsB,sigmaB)# XRay sigma
[dyyB,dxxB] = meshgrid(dataB.y,dataB.x)# Data grid for plots
dataB.data=gaussian_filter(dataB.data,0.5) 
# 1-sigma smoothing of X-Ray data

dataC=pysubs.Array2d(dxminA,dxmaxA,xpixelsA,dyminA,dymaxA,ypixelsA)
dataC=pysubs.GetBulletData(datapathC,dataC)# Measured SZE data
sigmaC=pysubs.Array2d(dxminA,dxmaxA,xpixelsA,dyminA,dymaxA,ypixelsA)
sigmaC=pysubs.GetBulletData(datapathsC,sigmaC)# SZE Sigma

dataC.data = abs(dataC.data*1E6) #convert to micro Kelvin
sigmaC.data = sigmaC.data*1E6 #convert to micro Kelvin

szekernel=pysubs.Array2d(0.0,400.0,113,0.0,400.0,113)
szekernel=pysubs.GetBulletFits(datapathszekernel,szekernel)

kernel60=pysubs.Array2d(0.0,1.0,110,0.0,1.0,110)
kernel60=pysubs.GetBulletData(datapathkernel60,kernel60)

xscaleD=yscaleD=BulletConstants.AngularScale*3600*BulletConstants.XRayTempDataScale 
# Pixel size in kpc. 4.413 kpc/" is the angular scale at the bullet cluster
# 3600 "/degree, 1.09333 degrees/pixel is the data scale

xpixelsD=ypixelsD=80
sxpixelsD=sypixelsD=160
dxmaxD = xscaleD*xpixelsD/2
dymaxD = yscaleD*ypixelsD/2
sxmaxD = dxmaxD*2
symaxD = dymaxD*2
dataD=pysubs.Array2d(-dxmaxD,dxmaxD,xpixelsD,-dymaxD,dymaxD,ypixelsD)
dataD=pysubs.GetBulletData(datapathD,dataD)# XRay measured data
#sigmaD=pysubs.Array2d(-dxmaxD,dxmaxD,xpixelsD,-dymaxD,dymaxD,ypixelsD)
#sigmaD=pysubs.GetBulletData(datapathsD,sigmaD)# XRay sigma
[dyyD,dxxD] = meshgrid(dataD.y,dataD.x)# Data grid for plots

dmlevels=linspace(0,0.4,40) 
for i in range(40): dmlevels[i]=int(dmlevels[i]*100)/100.0
xlevels=linspace(0.0,5E-6,25) 
#for i in range(25): xlevels[i]=int(xlevels[i]*10)/10.0

templevels=linspace(0.0,20,25) 
for i in range(25): xlevels[i]=int(xlevels[i]*10)/10.0
figure()
suptitle('Bullet Cluster Data Overlay', fontsize=16)
subplots_adjust(hspace=0.4, wspace=0.4)

subplot("221")
title("Lensing data")
contourf(dxxA,dyyA,dataA.data)
for c in contourf(dxxA,dyyA,dataA.data).collections:
    c.set_linewidth(0.1)
    c.set_alpha(1.0)
colorbar()
xlabel('kpc')
ylabel('kpc')
gca().set_xticks([-800.0,-400.0,0.0,400.0,800.0])

subplot("222")
title("X-Ray data")
contourf(dxxB,dyyB,dataB.data,cmap=cm.spectral)
for c in contourf(dxxB,dyyB,dataB.data,cmap=cm.spectral).collections:
    c.set_linewidth(0.1)
    c.set_alpha(1.0)
colorbar()
xlabel('kpc')
ylabel('kpc')
gca().set_xticks([-800.0,-400.0,0.0,400.0,800.0])

subplot("223")

title('Both Data Sets')

contourf(dxxB,dyyB,dataB.data,cmap=cm.spectral)
for c in contourf(dxxB,dyyB,dataB.data,cmap=cm.spectral).collections:
    c.set_linewidth(0.1)
    c.set_alpha(1.0)
colorbar()

contour(dxxA,dyyA,dataA.data, colors='y', linewidth=0.1)
xlabel('kpc')
ylabel('kpc')
gca().set_xticks([-800.0,-400.0,0.0,400.0,800.0])
savefig("XRay_Mass_Overlay.png")

figure()
suptitle('Bullet Cluster Data Overlay', fontsize=16)
subplots_adjust(hspace=0.4, wspace=0.4)

subplot("221")
title("XRay Temp data")
contourf(dxxD,dyyD,dataD.data)
for c in contourf(dxxD,dyyD,dataD.data).collections:
    c.set_linewidth(0.1)
    c.set_alpha(1.0)
colorbar()
xlabel('kpc')
ylabel('kpc')
gca().set_xticks([-800.0,-400.0,0.0,400.0,800.0])

subplot("222")
title("X-Ray data")
contourf(dxxB,dyyB,dataB.data,cmap=cm.spectral)
for c in contourf(dxxB,dyyB,dataB.data,cmap=cm.spectral).collections:
    c.set_linewidth(0.1)
    c.set_alpha(1.0)
colorbar()
xlabel('kpc')
ylabel('kpc')
gca().set_xticks([-800.0,-400.0,0.0,400.0,800.0])

subplot("223")

title('Both Data Sets')

contourf(dxxD,dyyD,dataD.data,cmap=cm.spectral)
for c in contourf(dxxD,dyyD,dataD.data,cmap=cm.spectral).collections:
    c.set_linewidth(0.1)
    c.set_alpha(1.0)
colorbar()

contour(dxxB,dyyB,dataB.data, colors='k', linewidth=0.2)
xlabel('kpc')
ylabel('kpc')
gca().set_xticks([-800.0,-400.0,0.0,400.0,800.0])
savefig("XRay_Temp_Overlay.png")


#************END MAIN PROGRAM*************************


