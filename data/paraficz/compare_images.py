
#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: May 22, 2013

# Testing the data interpolation for shifts and rotations
# Errors always less than 0.5%

import lageconfig # system specific path information
import sys
sys.path.append(lageconfig.bulletpath)
import BulletConstants # Constants used in Bullet cluster simulations
from pylab import *
import time
from subprocess import *
import pysubs_nas_20May13 as pysubs
import pyfits

#***************SUBROUTINES******************************

def ComparisonPlot(data,sim,title,filled_levels=None,line_levels=None,cmap=None,line_plot_multiplier=1.0,line_plot_yticks=[0.0,1.0],simtime=0.80,fom=4.0,line=[80,53,0.26],legend_location=[0.50,1.0],take_log=False):
    fig = figure()
    [dyy,dxx] = meshgrid(data.y,data.x)# Data grid for plots
    [dataslice,simslice] = pysubs.OneDSlice2Fixed(data,sim,line[2],line[0],line[1])
    plotdata=pysubs.Array2d(data.xmin,data.xmax,data.nx,data.ymin,data.ymax,data.ny)
    plotsim=pysubs.Array2d(sim.xmin,sim.xmax,sim.nx,sim.ymin,sim.ymax,sim.ny)
    if take_log:
        plotdata.data = log10(data.data)
        plotsim.data = log10(sim.data)
    else:
        plotdata.data = data.data
        plotsim.data = sim.data

    suptitle('Mass Lensing Map Comparison', fontsize=16)
    ax1=axes([0.15,0.5,0.4,0.4],aspect=1)
    ax1.set_title("Bradac Map")
    cont1 = ax1.contourf(dxx,dyy,plotdata.data,filled_levels,cmap=cmap)
    for c in cont1.collections:
        c.set_linewidth(0.1)
        c.set_alpha(1.0)
    ax1.set_xticks([])
    ax1.set_yticks([-500.0,0.0,500.0])
    ax1.set_ylabel('kpc',labelpad = -15)

    ax2=axes([0.475,0.5,0.4,0.4],aspect=1)
    ax2.set_title("Paraficz Map")
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
    ax4.plot(dataslice.x,dataslice.data*line_plot_multiplier,label="B Data",linewidth=3.0)
    ax4.plot(simslice.x,simslice.data*line_plot_multiplier,label="P Data",linewidth=3.0)
    ax4.set_xticks([-500.0,0.0,500.0])
    ax4.set_xlabel('kpc')
    ax4.set_yticks([0.0,20.0,40.0,60.0,80.0])
    ax4.set_ylim(0.0,100.0)
    ax4.yaxis.tick_right()
    ax4.legend(loc=1,bbox_to_anchor=(1.04,1.0),ncol=2)
    ltext = ax4.get_legend().get_texts()
    setp(ltext, fontsize = 8, color = 'k')
    return fig



#**************************MAIN PROGRAM********************************

toppath = lageconfig.toppath
datapathA=toppath+'bullet/data/kappa_25Apr12.dat'

xscaleA=yscaleA=BulletConstants.AngularScale*3600*BulletConstants.MassDataScale 
# Pixel size in kpc. 4.413 kpc/" is the angular scale at the bullet cluster
# 3600 "/degree, 9.86E-4 degrees/pixel is the data scale
start=time.time()
xpixelsA=ypixelsA=110
dxmaxA = xscaleA*xpixelsA/2
dymaxA = yscaleA*ypixelsA/2
dxminA = -xscaleA*xpixelsA/2
dyminA = -yscaleA*ypixelsA/2
dataA=pysubs.Array2d(dxminA,dxmaxA,xpixelsA,dyminA,dymaxA,ypixelsA)

dataA=pysubs.GetBulletData(datapathA,dataA)# Mass density measured data

bfitsfile=toppath+'bullet/data/kappa_23-Apr-12/1e0657.files_kappa_rs_0.9_rescaled.fits'
bfits=pysubs.Array2d(dxminA,dxmaxA,xpixelsA,dyminA,dymaxA,ypixelsA)
NormFactor = 67.968

fits = pyfits.open(bfitsfile)
bfits.data = fits[0].data * NormFactor

pfitsfile=toppath+'bullet/data/paraficz/mass_rescaled.fits'
pfits=pysubs.Array2d(dxminA,dxmaxA,xpixelsA,dyminA,dymaxA,ypixelsA)
NormFactor = 100.0

fits = pyfits.open(pfitsfile)
pfits.data = fits[0].data * NormFactor

elapsed = time.time()-start
print "Elapsed time for data retrieval = %f"%elapsed
start=time.time()

btemp = copy(bfits.data)
ptemp = copy(pfits.data)
for i in range(dataA.nx):
    for j in range(dataA.ny):
        bfits.data[i,j] = btemp[j,i]
        pfits.data[i,j] = ptemp[j,i]

elapsed = time.time()-start
print "Elapsed time for rotation = %f"%elapsed
start=time.time()
print bfits.data.shape, pfits.data.shape

[filled_levels,line_levels]=pysubs.SetContourLevels(0.0,70.0)
fig = ComparisonPlot(bfits,pfits,"Mass Lensing",filled_levels=filled_levels,line_levels=line_levels,simtime=0.0,fom=0.0,line=[80,53,0.26],legend_location=[0.70, 0.25],line_plot_yticks=[0.0,20.0,40.0,60.0])
elapsed = time.time()-start
print "Elapsed time for plotting = %f"%elapsed

show()
