#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 27-Mar-12


#This is a standalone python program to run findfom in batch mode

import matplotlib
matplotlib.use("PDF")
from pylab import *
from scipy.ndimage import gaussian_filter, convolve
from scipy.ndimage.filters import sobel
from matplotlib.backends.backend_pdf import PdfPages
import pickle

class Array1d:
    def __init__(self,xmin,xmax,nx):
        self.nx=nx
        self.xmin=xmin
        self.xmax=xmax
        self.dx=(xmax-xmin)/nx
        self.x=linspace(xmin+self.dx/2,xmax-self.dx/2,nx)
        self.y=zeros([nx])# Used for 1D slices
        self.data=zeros([nx])

class Array2d:
    def __init__(self,xmin,xmax,nx,ymin,ymax,ny):
        self.nx=nx
        self.ny=ny

        self.xmin=xmin
        self.ymin=ymin
        
        self.xmax=xmax
        self.ymax=ymax
        
        self.dx=(xmax-xmin)/nx
        self.dy=(ymax-ymin)/ny
        
        self.x=linspace(xmin+self.dx/2,xmax-self.dx/2,nx)
        self.y=linspace(ymin+self.dy/2,ymax-self.dy/2,ny)

        self.data=zeros([nx,ny])

def OneDSlice2Fixed(dataA,dataB,theta,imax,jmax):
	# This version allows you to specify the point it goes through.
	nx = dataA.nx
	xminline = min( dataA.x[imax] - (dataA.y[jmax] - dataA.ymin) / tan(theta) ,dataA.x[imax] + (-dataA.y[jmax] + dataA.ymax) / tan(theta))
	xmaxline = max( dataA.x[imax] - (dataA.y[jmax] - dataA.ymin) / tan(theta) ,dataA.x[imax] + (-dataA.y[jmax] + dataA.ymax) / tan(theta))
	xmin = max(dataA.xmin, xminline)
	xmax = min(dataA.xmax, xmaxline)
        sliceA=Array1d(xmin,xmax,nx)
        sliceB=Array1d(xmin,xmax,nx)
	for i in range(nx):
            xprime=sliceA.x[i]
            yprime=dataA.y[jmax] - (dataA.x[imax] - xprime) * tan(theta)
            sliceA.y[i]=sliceB.y[i]=yprime
            sliceA.data[i]=DataInterpolate(dataA,xprime,yprime)
            sliceB.data[i]=DataInterpolate(dataB,xprime,yprime)
        return [sliceA,sliceB]

def PyramidalKernel(deltax,deltay):
	if deltax>=1.0 or deltay>=1.0:
		return 0.0
	else:
		return (1.0-deltax)*(1.0-deltay)

def DataInterpolate(data,xprime,yprime):
	i=int((xprime-data.xmin)/data.dx)
	j=int((yprime-data.ymin)/data.dy)
	d=0.0
	for m in [i-1,i,i+1]:
	    ml = max(0,min(data.nx-1,m))
            deltax=abs((xprime-data.x[ml])/data.dx)
            for n in [j-1,j,j+1]:
	        nl = max(0,min(data.ny-1,n))
                deltay=abs((yprime-data.y[nl])/data.dy)
                d=d+PyramidalKernel(deltax,deltay)*data.data[ml,nl]
	return d


#****************MAIN PROGRAM*****************

# Unpickle the data
pickled_data = open('plotdata.pkl', 'rb')

data = pickle.load(pickled_data)
sigmas = pickle.load(pickled_data)
masks = pickle.load(pickled_data)
sims = pickle.load(pickled_data)
pickled_data.close()

dataA=data['dataA'];dataB1=data['dataB1'];dataB2=data['dataB2'];dataB3=data['dataB3'];dataC=data['dataC'];dataD=data['dataD'];dataE=data['dataE'];
sigmaA=sigmas['sigmaA'];sigmaB1=sigmas['sigmaB1'];sigmaB2=sigmas['sigmaB2'];sigmaB3=sigmas['sigmaB3'];sigmaC=sigmas['sigmaC'];sigmaD=sigmas['sigmaD'];sigmaE=sigmas['sigmaE'];
maskA=masks['maskA'];maskANull=masks['maskANull'];maskD=masks['maskD'];maskDNull=masks['maskDNull']

shiftedmsim=sims['shiftedmsim'];shiftedxsim1=sims['shiftedxsim1'];shiftedxsim2=sims['shiftedxsim2'];shiftedxsim3=sims['shiftedxsim3'];shiftedszsim=sims['shiftedszsim'];shiftedtempsim=sims['shiftedtempsim'];shiftedsynchsim=sims['shiftedsynchsim'];shiftedbmagsim=sims['shiftedbmagsim'];

pp=PdfPages('Graph_Test.pdf')
simtime=0.79
fom = 4.60
[dyyA,dxxA] = meshgrid(dataA.y,dataA.x)# Data grid for plots
[dyyB,dxxB] = meshgrid(dataB1.y,dataB1.x)# Data grid for plots
[dyyD,dxxD] = meshgrid(dataD.y,dataD.x)# Data grid for plots


dmlevels=linspace(0,5.0,25) 
for i in range(25): dmlevels[i]=int(dmlevels[i]*10)/10.0

xlevels=linspace(0.0,10.0,25) 
for i in range(25): xlevels[i]=int(xlevels[i]*10)/10.0

# First the mass plots
figure()
suptitle('Mass Lensing, T = %.2f Gy'%simtime, fontsize=16)
subplots_adjust(hspace=0.2, wspace=0.0)

subplot("221",aspect=1)
title("Data")
contourf(dxxA,dyyA,dataA.data)
for c in contourf(dxxA,dyyA,dataA.data).collections:
    c.set_linewidth(0.1)
    c.set_alpha(1.0)
gca().set_xticks([])

subplot("222",aspect=1)
title("Simulation")
contourf(dxxA,dyyA,shiftedmsim.data)
for c in contourf(dxxA,dyyA,shiftedmsim.data).collections:
    c.set_linewidth(0.1)
    c.set_alpha(1.0)
colorbar()
gca().set_xticks([])
gca().set_yticks([])

[dataslice,simslice] = OneDSlice2Fixed(dataA,shiftedmsim,0.26,80,53)

subplot("223",aspect=1)

title('Overlay: $\chi^2$ = %.2f '%fom)

contourf(dxxA,dyyA,dataA.data)
for c in contourf(dxxA,dyyA,dataA.data).collections:
    c.set_linewidth(0.1)
    c.set_alpha(1.0)
contour(dxxA,dyyA,shiftedmsim.data, colors='k')
plot(dataslice.x,dataslice.y,'w--',linewidth=2.0)# Plot the 1D slice
gca().set_xticks([-800.0,-400.0,0.0,800.0,400.0])
xlabel('kpc')

subplot("224")
title('1D Slice through data max')
plot(dataslice.x,dataslice.data,label="Data",linewidth=3.0)
plot(simslice.x,simslice.data,label="Sim",linewidth=3.0)

xlabel('kpc')
legend(loc=1,bbox_to_anchor=[0.7, 0.3])
ltext = gca().get_legend().get_texts()
setp(ltext, fontsize = 9, color = 'k')
pp.savefig()
pp.close()
#************END MAIN PROGRAM*************************

