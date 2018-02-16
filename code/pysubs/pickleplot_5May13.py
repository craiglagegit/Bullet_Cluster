#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 27-Mar-12


#This is a standalone python program to run findfom in batch mode

import lageconfig # system specific path information
import sys
import matplotlib
matplotlib.use("PDF")
sys.path.append(lageconfig.bulletpath)
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

def SetContourLevels(min,max):
    filled_levels=linspace(min,max,25) 
    for i in range(25): filled_levels[i]=int(filled_levels[i]*10)/10.0
    line_levels=linspace(min,max,10) 
    return [filled_levels,line_levels]

def ComparisonPlot(data,sim,title,filled_levels=None,line_levels=None,cmap=None,line_plot_multiplier=1.0,line_plot_yticks=[0.0,1.0],simtime=0.80,fom=4.0,line=[80,53,0.26],legend_location=[0.50,1.0],take_log=False):
    fig = figure()
    [dyy,dxx] = meshgrid(data.y,data.x)# Data grid for plots
    [dataslice,simslice] = OneDSlice2Fixed(data,sim,line[2],line[0],line[1])
    plotdata=Array2d(data.xmin,data.xmax,data.nx,data.ymin,data.ymax,data.ny)
    plotsim=Array2d(sim.xmin,sim.xmax,sim.nx,sim.ymin,sim.ymax,sim.ny)
    if take_log:
        plotdata.data = log10(data.data)
        plotsim.data = log10(sim.data)
    else:
        plotdata.data = data.data
        plotsim.data = sim.data

    suptitle(title+', T = %.2f Gy'%simtime, fontsize=16)
    ax1=axes([0.15,0.5,0.4,0.4],aspect=1)
    ax1.set_title("Data")
    cont1 = ax1.contourf(dxx,dyy,plotdata.data,filled_levels,cmap=cmap)
    for c in cont1.collections:
        c.set_linewidth(0.1)
        c.set_alpha(1.0)
    ax1.set_xticks([])
    ax1.set_yticks([-500.0,0.0,500.0])
    ax1.set_ylabel('kpc',labelpad = -15)

    ax2=axes([0.475,0.5,0.4,0.4],aspect=1)
    ax2.set_title("Simulation")
    cont2 = ax2.contourf(dxx,dyy,plotsim.data,filled_levels,cmap=cmap)
    for c in cont2.collections:
        c.set_linewidth(0.1)
        c.set_alpha(1.0)
    cb = colorbar(cont2, pad = 0.02)
    ax2.set_xticks([])
    ax2.set_yticks([])

    ax3=axes([0.15,0.09,0.4,0.4],aspect=1)
    cont3A = ax3.contourf(dxx,dyy,plotdata.data,filled_levels,cmap=cmap)
    for c in cont3A.collections:
        c.set_linewidth(0.1)
        c.set_alpha(1.0)

    cont3B = ax3.contour(dxx,dyy,plotsim.data,line_levels,linestyles='solid',colors='k')
    ax3.plot(dataslice.x,dataslice.y,'w--',linewidth=2.0)# Plot the 1D slice
    ax3.set_xticks([-500.0,0.0,500.0])
    ax3.set_yticks([-500.0,0.0,500.0])
    ax3.set_xlabel('kpc')
    ax3.set_ylabel('kpc',labelpad = -15)

    ax4=axes([0.507,0.09,0.30,0.4])
    ax4.plot(dataslice.x,dataslice.data*line_plot_multiplier,label="Data",linewidth=3.0)
    ax4.plot(simslice.x,simslice.data*line_plot_multiplier,label="Sim",linewidth=3.0)
    ax4.set_xticks([-500.0,0.0,500.0])
    ax4.set_xlabel('kpc')
    ax4.yaxis.tick_right()
    ax4.set_yticks(line_plot_yticks)
    ax4.legend(loc=1,bbox_to_anchor=legend_location)
    ltext = ax4.get_legend().get_texts()
    setp(ltext, fontsize = 9, color = 'k')
    return fig


def MultiLinePlot(data,sim,title,filled_levels=None,line_levels=None,cmap=None,line_plot_yticks=[0.0,1.0],simtime=0.80,fom=4.0,line1=[78,56,-0.24],line2=[78,56,0.24],line3=[78,56,1.02],legend_location1=[0.50,1.0],legend_location2=[0.50,1.0],legend_location3=[0.50,1.0],take_log=False):
    fig = figure()
    [dyy,dxx] = meshgrid(data.y,data.x)# Data grid for plots
    [dataslice1,simslice1] = OneDSlice2Fixed(data,sim,line1[2],line1[0],line1[1])
    [dataslice2,simslice2] = OneDSlice2Fixed(data,sim,line2[2],line2[0],line2[1])
    [dataslice3,simslice3] = OneDSlice2Fixed(data,sim,line3[2],line3[0],line3[1])
    plotdata=Array2d(data.xmin,data.xmax,data.nx,data.ymin,data.ymax,data.ny)
    plotsim=Array2d(sim.xmin,sim.xmax,sim.nx,sim.ymin,sim.ymax,sim.ny)
    if take_log:
        plotdata.data = log10(data.data)
        plotsim.data = log10(sim.data)
    else:
        plotdata.data = data.data
        plotsim.data = sim.data
    suptitle(title+', T = %.2f Gy'%simtime, fontsize=16)
    ax1=axes([0.15,0.5,0.4,0.4],aspect=1)
    ax1.set_title("Data")
    cont1 = ax1.contourf(dxx,dyy,plotdata.data,filled_levels,cmap=cmap)
    for c in cont1.collections:
        c.set_linewidth(0.1)
        c.set_alpha(1.0)
    ax1.plot(dataslice1.x,dataslice1.y,'w--',linewidth=2.0)# Plot the 1D slice
    ax1.plot(dataslice2.x,dataslice2.y,'y--',linewidth=2.0)# Plot the 1D slice
    ax1.plot(dataslice3.x,dataslice3.y,'k--',linewidth=2.0)# Plot the 1D slice
    ax1.set_xticks([])
    ax1.set_yticks([-500.0,0.0,500.0])
    ax1.set_ylabel('kpc',labelpad = -15)

    ax2=axes([0.507,0.50,0.30,0.4])
    ax2.set_title('1D Slice (Log - White)',x=0.50,y=0.90)
    ax2.plot(dataslice1.x,log10(dataslice1.data),label="Data",linewidth=3.0)
    ax2.plot(simslice1.x,log10(simslice1.data),label="Sim",linewidth=3.0)
    ax2.set_xticks([-500.0,0.0,500.0])
    ax2.yaxis.tick_right()
    ax2.set_yticks(line_plot_yticks)
    ax2.legend(loc=1,bbox_to_anchor=legend_location1)
    ltext = ax2.get_legend().get_texts()
    setp(ltext, fontsize = 9, color = 'k')

    ax3=axes([0.20,0.09,0.30,0.4])
    ax3.set_title('1D Slice (Log - Yellow)',x=0.50,y=0.90)
    ax3.plot(dataslice2.x,log10(dataslice2.data),label="Data",linewidth=3.0)
    ax3.plot(simslice2.x,log10(simslice2.data),label="Sim",linewidth=3.0)
    ax3.set_xticks([-500.0,0.0,500.0])
    ax3.set_xlabel('kpc')
    ax3.set_yticks(line_plot_yticks)
    ax3.legend(loc=1,bbox_to_anchor=legend_location2)
    ltext = ax3.get_legend().get_texts()
    setp(ltext, fontsize = 9, color = 'k')

    ax4=axes([0.507,0.09,0.30,0.4])
    ax4.set_title('1D Slice (Log - Black)',x=0.50,y=0.90)
    ax4.plot(dataslice3.x,log10(dataslice3.data),label="Data",linewidth=3.0)
    ax4.plot(simslice3.x,log10(simslice3.data),label="Sim",linewidth=3.0)
    ax4.set_xticks([0.0,500.0])
    ax4.set_xlabel('kpc')
    ax4.yaxis.tick_right()
    ax4.set_yticks(line_plot_yticks)
    ax4.legend(loc=1,bbox_to_anchor=legend_location3)
    ltext = ax4.get_legend().get_texts()
    setp(ltext, fontsize = 9, color = 'k')
    return fig

def EdgeDetectionPlot(data,sim,title,filled_levels=None,ED_filled_levels=None,ED_line_levels=None,ED_multiplier=1.0,cmap=None,simtime=0.80,fom=4.0,take_log=False):
    fig = figure()
    [dyy,dxx] = meshgrid(data.y,data.x)# Data grid for plots
    plotdata=Array2d(data.xmin,data.xmax,data.nx,data.ymin,data.ymax,data.ny)
    dataedge=Array2d(data.xmin,data.xmax,data.nx,data.ymin,data.ymax,data.ny)
    datatrunc=Array2d(data.xmin,data.xmax,data.nx,data.ymin,data.ymax,data.ny)
    plotsim=Array2d(sim.xmin,sim.xmax,sim.nx,sim.ymin,sim.ymax,sim.ny)
    simedge=Array2d(sim.xmin,sim.xmax,sim.nx,sim.ymin,sim.ymax,sim.ny)
    simtrunc=Array2d(sim.xmin,sim.xmax,sim.nx,sim.ymin,sim.ymax,sim.ny)

    xraymax = pow(10, -6.8)
    xraymin = pow(10, -8.0)
    for ii in range(data.nx):
        for jj in range(data.ny):
            datatrunc.data[ii,jj] = max(min(xraymax, data.data[ii,jj]),xraymin)
            simtrunc.data[ii,jj] = max(min(xraymax, sim.data[ii,jj]),xraymin)

    dataedge.data = hypot(sobel(datatrunc.data, 0),sobel(datatrunc.data, 1))
    simedge.data = hypot(sobel(simtrunc.data, 0),sobel(simtrunc.data, 1))

    if take_log:
        plotdata.data = log10(data.data)
        plotsim.data = log10(sim.data)
    else:
        plotdata.data = data.data
        plotsim.data = sim.data

    suptitle(title+', T = %.2f Gy'%simtime, fontsize=16)
    ax1=axes([0.165,0.5,0.4,0.4],aspect=1)
    ax1.set_title("Data")
    cont1 = ax1.contourf(dxx,dyy,plotdata.data,filled_levels,cmap=cmap)
    for c in cont1.collections:
        c.set_linewidth(0.1)
        c.set_alpha(1.0)
    ax1.set_xticks([])
    ax1.set_yticks([-500.0,0.0,500.0])
    ax1.set_ylabel('kpc',labelpad = -15)

    ax2=axes([0.475,0.5,0.4,0.4],aspect=1)
    ax2.set_title("Edge Detection")
    ax2.text(0.0,0.0,"Data",color="w",ha='center',va='center')
    cont2 = ax2.contourf(dxx,dyy,dataedge.data*ED_multiplier,ED_filled_levels,cmap=cmap)
    for c in cont2.collections:
        c.set_linewidth(0.1)
        c.set_alpha(1.0)
    ax2.set_xticks([])
    ax2.yaxis.tick_right()
    ax2.set_yticks([])

    ax3=axes([0.165,0.09,0.4,0.4],aspect=1)
    ax3.text(0.0,0.0,"Simulation",color="w",ha='center',va='center')
    cont3 = ax3.contourf(dxx,dyy,simedge.data*ED_multiplier,ED_filled_levels,cmap=cmap)
    for c in cont3.collections:
        c.set_linewidth(0.1)
        c.set_alpha(1.0)
    ax3.set_xticks([-500.0,0.0,500.0])
    ax3.set_yticks([-500.0,0.0,500.0])
    ax3.set_xlabel('kpc')
    ax3.set_ylabel('kpc',labelpad = -15)

    ax4=axes([0.475,0.09,0.4,0.4],aspect=1)
    ax4.text(0.0,0.0,"Overlay",color="w",ha='center',va='center')
    cont4A = ax4.contourf(dxx,dyy,dataedge.data*ED_multiplier,ED_filled_levels,cmap=cmap)
    for c in cont4A.collections:
        c.set_linewidth(0.1)
        c.set_alpha(1.0)

    cont4B = ax4.contour(dxx,dyy,simedge.data*ED_multiplier,ED_line_levels,linestyles='solid',colors='k')
    ax4.set_xticks([-500.0,0.0,500.0])
    ax4.set_xlabel('kpc')
    ax4.yaxis.tick_right()
    ax4.set_yticks([])
    return fig

def FomLocationPlot(datasets,subtitles,plotsets,mask_sum,cmap=None,simtime=0.80,fom=4.0):
    fig = figure()
    suptitle("$\chi^2$ by dataset, T = %.2f Gy"%simtime, fontsize=16)
    counter=0
    for j in range(2):
        for i in range(3):
            [dyy,dxx] = meshgrid(plotsets[counter].y,plotsets[counter].x)# Data grid for plots
            chisquared = datasets[counter].sum() / mask_sum
            plotx = 0.08 + 0.3 * i
            ploty = 0.60 - 0.3 * j
            ax = axes([plotx,ploty,0.3,0.3])
            ax.set_title(subtitles[counter]+", $\chi^2$ = %.2f"%chisquared,x=0.5,y=0.85,color='w',fontsize=14)
            cont = ax.contourf(dxx,dyy,datasets[counter],cmap=cmap)
            for c in cont.collections:
                c.set_linewidth(0.1)
                c.set_alpha(1.0)
            if i == 0:
                ax.set_yticks([-500.0,0.0,500.0])
                ax.set_ylabel('kpc',labelpad = -15)
            else:
                ax.set_yticks([])
            if j == 1:
                ax.set_xticks([-500.0,0.0,500.0])
                ax.set_xlabel('kpc')
            else:
                ax.set_xticks([])
            counter = counter + 1
    return fig

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
FomLocate=True

[filled_levels,line_levels]=SetContourLevels(0.0,70.0)
fig = ComparisonPlot(dataA,shiftedmsim,"Mass Lensing",filled_levels=filled_levels,line_levels=line_levels,simtime=simtime,fom=fom,line=[80,53,0.26],legend_location=[0.80, 0.3],line_plot_yticks=[0.0,20.0,40.0,60.0])
pp.savefig(fig)

[filled_levels,line_levels]=SetContourLevels(-10.0,-4.0)
fig = ComparisonPlot(dataB1,shiftedxsim1,"X-ray Flux - 500-2000eV",filled_levels=filled_levels,line_levels=line_levels,cmap=cm.spectral,simtime=simtime,fom=fom,line=[70,55,0.24],legend_location=[0.50, 1.0],line_plot_yticks=[0.0,1.0,2.0,3.0],line_plot_multiplier=1.0E6,take_log=True)
pp.savefig(fig)

fig = ComparisonPlot(dataB2,shiftedxsim2,"X-ray Flux - 2000-5000eV",filled_levels=filled_levels,line_levels=line_levels,cmap=cm.spectral,simtime=simtime,fom=fom,line=[70,55,0.24],legend_location=[0.50, 1.0],line_plot_yticks=[0.0,0.5,1.0],line_plot_multiplier=1.0E6,take_log=True)
pp.savefig(fig)

fig = ComparisonPlot(dataB3,shiftedxsim3,"X-ray Flux - 5000-8000eV",filled_levels=filled_levels,line_levels=line_levels,cmap=cm.spectral,simtime=simtime,fom=fom,line=[70,55,0.24],legend_location=[0.50, 1.0],line_plot_yticks=[0.0,1.0,2.0,3.0],line_plot_multiplier=1.0E7,take_log=True)
pp.savefig(fig)

[filled_levels,line_levels]=SetContourLevels(-10.0,-4.0)
fig = MultiLinePlot(dataB1,shiftedxsim1,"X-ray Flux - 500-2000eV",filled_levels=filled_levels,line_levels=line_levels,cmap=cm.spectral,simtime=simtime,fom=fom,line1=[78,56,-0.24],line2=[78,56,0.24],line3=[78,56,1.02],legend_location1=[0.70,0.30],legend_location2=[0.70,0.30],legend_location3=[0.70,0.30],line_plot_yticks=[-8.0,-7.0,-6.0],take_log=True)

[filled_levels,line_levels]=SetContourLevels(-10.0,-4.0)
ED_filled_levels=linspace(0,10.0,11) 
ED_line_levels=linspace(1.0,5.0,5)
fig = EdgeDetectionPlot(dataB1,shiftedxsim1,"X-ray Shock Edge Detection",filled_levels=filled_levels,ED_filled_levels=ED_filled_levels,ED_line_levels=ED_line_levels,ED_multiplier=1.0E7,cmap=cm.spectral,simtime=simtime,fom=fom,take_log=True)
pp.savefig(fig)

[filled_levels,line_levels]=SetContourLevels(0.0,20.0)
fig = ComparisonPlot(dataD,shiftedtempsim,"Gas Temperature (keV)",filled_levels=filled_levels,line_levels=line_levels,cmap=cm.spectral,simtime=simtime,fom=fom,line=[70,52,0.24],legend_location=[0.60, 1.0],line_plot_yticks=[5.0,10.0,15.0])
pp.savefig(fig)

[filled_levels,line_levels]=SetContourLevels(0.0,400.0)
fig = ComparisonPlot(dataC,shiftedszsim,"SZ Temperature Decrement ($\mu K$)",filled_levels=filled_levels,line_levels=line_levels,simtime=simtime,fom=fom,line=[70,59,0.24],legend_location=[0.70, 0.32],line_plot_yticks=[0.0,100.0,200.0,300.0])
pp.savefig(fig)

filled_levels = [0.0,6.0,12.0,24.0,36.0,48.0,96.0]
line_levels = [0.0,6.0,12.0,24.0,36.0,48.0,96.0]
fig = ComparisonPlot(dataE,shiftedsynchsim,"Radio Flux - 1.3 GHz ($\mu$ Jy / pixel)",filled_levels=filled_levels,line_levels=line_levels,cmap=cm.spectral,simtime=simtime,fom=fom,line=[70,59,0.24],legend_location=[1.02, 1.0],line_plot_yticks=[0.0,25.0,50.0,75.0])
pp.savefig(fig)

if FomLocate:
    numpixels = dataA.nx * dataA.ny
    masschisquared = pow((dataA.data-shiftedmsim.data)/sigmaA.data,2) * maskA.data
    xray1chisquared = pow((dataB1.data-shiftedxsim1.data)/sigmaB1.data,2) * maskA.data
    xray2chisquared = pow((dataB2.data-shiftedxsim2.data)/sigmaB2.data,2) * maskA.data
    xray3chisquared = pow((dataB3.data-shiftedxsim3.data)/sigmaB3.data,2) * maskA.data
    szechisquared = pow((dataC.data-shiftedszsim.data)/sigmaC.data,2) * maskA.data
    tempchisquared = pow((dataD.data-shiftedtempsim.data)/sigmaD.data,2) * maskD.data
    mask_sum = maskA.data.sum()
    datasets = [masschisquared, xray1chisquared, xray2chisquared, xray3chisquared, tempchisquared, szechisquared]
    plotsets = [dataA,dataB1,dataB2,dataB3,dataD,dataC]
    subtitles = ["Mass","X-ray 1","X-ray 2","X-ray 3","Temp","SZE"]
    fig = FomLocationPlot(datasets,subtitles,plotsets,mask_sum,cmap=cm.spectral,simtime=simtime,fom=fom)
    pp.savefig(fig)
pp.close()
#************END MAIN PROGRAM*************************

