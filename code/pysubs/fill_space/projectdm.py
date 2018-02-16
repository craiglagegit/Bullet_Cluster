#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 5-Sep-11

'''
Playing with Ctypes vs Python vs Cython

'''
import lageconfig # system specific path information
from pylab import *
import sys
sys.path.append(lageconfig.bulletpath)
import BulletConstants # Constants used in Bullet cluster simulations
##from subprocess import *
import time
import array

import cython

#from projectdm3 import ProjectGadgetDMCython
from projectdm4 import ProjectGadgetDMCython, ProjectGadgetDMCython2
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



class Array3d:
    def __init__(self,xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz):
        self.nx=nx
        self.ny=ny
        self.nz=nz

        self.xmin=xmin
        self.ymin=ymin
        self.zmin=zmin
        
        self.xmax=xmax
        self.ymax=ymax
        self.zmax=zmax
        
        self.dx=(xmax-xmin)/nx
        self.dy=(ymax-ymin)/ny
        self.dz=(zmax-zmin)/nz
        
        self.x=linspace(xmin+self.dx/2,xmax-self.dx/2,nx)
        self.y=linspace(ymin+self.dy/2,ymax-self.dy/2,ny)
        self.z=linspace(zmin+self.dz/2,zmax-self.dz/2,nz)

        self.data=zeros([nx,ny,nz])

class Grid3d:
    def __init__(self,xmin,xmax,ymin,ymax,zmin,zmax):
       
        self.i=range(xmin,xmax+1)
        self.j=range(ymin,ymax+1)
        self.k=range(zmin,zmax+1)

class Align:
    def __init__(self):
        self.d=zeros([3])
        self.dmin=zeros([3])
        self.dmax=zeros([3])


def ReadGadgetSnapshot(filename):# Reads a snapshot file

    class Snapfile:
        pass

    snap = Snapfile()
    file = open(filename,"rb")
    snap.header = ReadGadgetHeader(file) # Header data
    snap.data = ReadGadgetData(file, snap.header,0) # Particle data, 0 = snapshot file, 1-IC file
    file.close()

    return snap


def ReadGadgetHeader(file): # Reads the header of the snapshot or IC file

    class HeadData:
        pass
    header=HeadData()

    header.blocksize = array.array("i")
    header.Npart = array.array("I")
    header.Massarr = array.array("d")
    header.Time = array.array("d")
    header.Redshift = array.array("d")
    header.FlagSfr = array.array("i")
    header.FlagFeedback = array.array("i")
    header.Nall = array.array("i")
    header.FlagCooling = array.array("i")
    header.NumFiles = array.array("i")
    header.BoxSize = array.array("d")
    header.Omega0 = array.array("d")
    header.OmegaLambda = array.array("d")
    header.HubbleParam = array.array("d")
    header.FlagAge = array.array("i")
    header.FlagMetals = array.array("i")
    header.NallHW = array.array("i")
    header.flag_entr_ics = array.array("i")
    header.unused = array.array("i")

    header.blocksize.fromfile(file,1)
    header.Npart.fromfile(file,6)
    header.Massarr.fromfile(file,6)
    header.Time.fromfile(file,1)
    header.Redshift.fromfile(file,1)
    header.FlagSfr.fromfile(file,1)
    header.FlagFeedback.fromfile(file,1)
    header.Nall.fromfile(file,6)
    header.FlagCooling.fromfile(file,1)
    header.NumFiles.fromfile(file,1)
    header.BoxSize.fromfile(file,1)
    header.Omega0.fromfile(file,1)
    header.OmegaLambda.fromfile(file,1)
    header.HubbleParam.fromfile(file,1)
    header.FlagAge.fromfile(file,1)
    header.FlagMetals.fromfile(file,1)
    header.NallHW.fromfile(file,6)
    header.flag_entr_ics.fromfile(file,1)
    header.unused.fromfile(file,15)
    header.blocksize.fromfile(file,1)

    return header


def ReadGadgetData(file,header,filetype): 
    # Reads the data (position, velocity, mass, ...) 
    # filetype=0->snapshot file, filetype=1->IC file

    class Data:
        pass
    data=Data()

    Npart=header.Npart
    Massarr=header.Massarr
    Ngas=header.Npart[0]
    N=0
    for i in range(6):
        N = N+Npart[i]
    Nm=0
    for i in range(6):
        if Massarr[i]==0:
            Nm = Nm+Npart[i]
    
    data.blocksize = array.array("i")
    data.pos = array.array("f")
    data.vel = array.array("f")
    data.id = array.array("I")
    data.masses = array.array("f")
    data.u = array.array("f")
    data.rho = array.array("f")
    if header.FlagCooling[0]==1:
	    data.ne = array.array("f")
	    data.np = array.array("f")
    data.hsml = array.array("f")

    data.blocksize.fromfile(file,1)
    data.pos.fromfile(file,3*N)
    data.blocksize.fromfile(file,1)

    data.blocksize.fromfile(file,1)
    data.vel.fromfile(file,3*N)
    data.blocksize.fromfile(file,1)

    data.blocksize.fromfile(file,1)
    data.id.fromfile(file,N)
    data.blocksize.fromfile(file,1)

    if Nm!=0:
        data.blocksize.fromfile(file,1)
        data.masses.fromfile(file,Nm)
        data.blocksize.fromfile(file,1)
    if Ngas!=0:      
        data.blocksize.fromfile(file,1)
        data.u.fromfile(file,Ngas)
        data.blocksize.fromfile(file,1)

   	if filetype==0:

	        data.blocksize.fromfile(file,1)
	        data.rho.fromfile(file,Ngas)
	        data.blocksize.fromfile(file,1)

		if header.FlagCooling[0]==1:

	        	data.blocksize.fromfile(file,1)
	        	data.ne.fromfile(file,Ngas)
	        	data.blocksize.fromfile(file,1)

	        	data.blocksize.fromfile(file,1)
	        	data.np.fromfile(file,Ngas)
	        	data.blocksize.fromfile(file,1)

        	data.blocksize.fromfile(file,1)
        	data.hsml.fromfile(file,Ngas)
        	data.blocksize.fromfile(file,1)

    return data


def PrintGadgetHeader(snapfile): # Prints the header data

    for i in range(6):
        print "Npart[",i,"] = ",snapfile.header.Npart[i]
        print "Massarr[",i,"] = ",snapfile.header.Massarr[i]
        print "Nall[",i,"] = ",snapfile.header.Nall[i]
        print "NallHW[",i,"] = ",snapfile.header.NallHW[i]

    print "Time = ", snapfile.header.Time[0]
    print "Redshift = ", snapfile.header.Redshift[0]
    print "FlagSfr = ", snapfile.header.FlagSfr[0]
    print "FlagFeedback = ", snapfile.header.FlagFeedback[0]
    print "FlagCooling = ", snapfile.header.FlagCooling[0]
    print "NumFiles = ", snapfile.header.NumFiles[0]
    print "BoxSize = ", snapfile.header.BoxSize[0]
    print "Omega0 = ", snapfile.header.Omega0[0]
    print "OmegaLambda = ", snapfile.header.OmegaLambda[0]
    print "HubbleParam = ", snapfile.header.HubbleParam[0]
    print "FlagAge = ", snapfile.header.FlagAge[0]
    print "FlagMetals = ", snapfile.header.FlagMetals[0]
    print "flag_entr_ics = ", snapfile.header.flag_entr_ics[0]

def PrintGadgetData(snapfile,parttype,nmin,nmax,filetype=0): 
    # Prints selected data
    # filetype=0->snapshot file, filetype=1->IC file

    offset=0
    for k in range(parttype):
        offset=offset+snapfile.header.Npart[k]

    for i in range(nmin,nmax):
        io = i+offset
        for j in range(3):
            print "PartType ",parttype,"Particle ",i," X",j," = ",snapfile.data.pos[3*io+j]
            print "PartType ",parttype,"Particle ",i," V",j," = ",snapfile.data.vel[3*io+j]
        print "PartType ",parttype,"Particle ",i," id = ",snapfile.data.id[io]
        print "PartType ",parttype,"Particle ",i," mass = ",snapfile.data.masses[io]
        if parttype==0:
            print "PartType ",parttype,"Particle ",i," u = ",snapfile.data.u[io]
	    if filetype==0:
            	print "PartType ",parttype,"Particle ",i," rho = ",snapfile.data.rho[io]
            	print "PartType ",parttype,"Particle ",i," ne = ",snapfile.data.ne[io]
            	print "PartType ",parttype,"Particle ",i," np = ",snapfile.data.np[io]
            	print "PartType ",parttype,"Particle ",i," hsml = ",snapfile.data.hsml[io]

def ReadGadgetParameterFile(filename):
    # This reads a Gadget parameter file and returns a set of lists with the parameters
    class ParamData:
        pass
    param=ParamData()

    try:
        file = open(filename,'r')
        lines=file.readlines()
        file.close()
        for line in lines:
		ThisLine=line.strip().split()
		ThisLineLength=len(ThisLine)
		if line[0]=='%' or line[0]=='\n' or ThisLineLength<2:
			continue
		#print ThisLine
		ParamName=ThisLine[0]
		#print ParamName
		try:
			ThisParam=int(ThisLine[1])
		except ValueError:
			try:
				ThisParam=float(ThisLine[1])
			except ValueError:
				ThisParam=ThisLine[1]
		#print ThisParam
		exec("param.%s=ThisParam" %ParamName)
    except IOError:
        print "Error reading Gadget parameter file"

    return param


def EulerAngles(phi=0.0,theta=0.0,psi=0.0):
    R=zeros([3,3])
    R[0,0]=cos(psi)*cos(phi)-cos(theta)*sin(phi)*sin(psi)
    R[0,1]=cos(psi)*sin(phi)+cos(theta)*cos(phi)*sin(psi)
    R[0,2]=sin(psi)*sin(theta)
    R[1,0]=-sin(psi)*cos(phi)-cos(theta)*sin(phi)*cos(psi)
    R[1,1]=-sin(psi)*sin(phi)+cos(theta)*cos(phi)*cos(psi)
    R[1,2]=cos(psi)*sin(theta)
    R[2,0]=sin(theta)*sin(phi)
    R[2,1]=-sin(theta)*cos(phi)
    R[2,2]=cos(theta)
    return R



def ProjectGadgetDMPython(snapfile,data,parttype,zmin,zmax,phi=0.0,theta=0.0,psi=0.0): 
    # Projects 3D mass density data onto a 2D grid
    # Euler angles are used to rotate the data if desired
    xpixels=data.nx
    ypixels=data.ny
    xy=xpixels*ypixels
    R=EulerAngles(phi,theta,psi)
    x=zeros([3]) # Position vector
    xp=zeros([3])# Rotated position vector
    #data.data=1E-12# Mass density

    Npart=snapfile.header.Npart
    num=Npart[parttype] # Chosen particle type

    start = time.time()

    offset=0
    for k in range(parttype):
	offset=offset+Npart[k]

    for n in range(num):

	no=offset+n
        for m in range(3):
            x[m]=snapfile.data.pos[3*no+m]
        xp=dot(R,x) # Rotate the position by the Euler angles


	i = int((xp[0]-data.xmin)/data.dx)
        j = int((xp[1]-data.ymin)/data.dy)
 	z = xp[2]

        if snapfile.header.Massarr[parttype]==0:
        	mass=snapfile.data.masses[no]
        else:
        	mass=snapfile.header.Massarr[parttype]

        if i>=0 and i<data.nx and j>=0 and j<data.ny and z>zmin and z<zmax:
            data.data[i,j]=data.data[i,j]+mass

    elapsed = (time.time()-start)
    print "Elapsed time to run dm projection (Python) = "+str(elapsed)

    return data


#************************ MAIN PROGRAM ********************************
filename = 'snapshots/snapshot19_000'
zmin = -3000.0
zmax =  3000.0
parttype = 1
data1 = Array2d(-3000.0,3000.0,128,-3000.0,3000.0,128)
data2 = Array2d(-3000.0,3000.0,128,-3000.0,3000.0,128)
data3 = Array2d(-3000.0,3000.0,128,-3000.0,3000.0,128)
data4 = Array2d(-3000.0,3000.0,128,-3000.0,3000.0,128)

snapfile=ReadGadgetSnapshot(filename)

data1 = ProjectGadgetDMPython(snapfile,data1,parttype,zmin,zmax,phi=0.0,theta=0.0,psi=0.0)
data3 = ProjectGadgetDMCython(snapfile,data3,parttype,zmin,zmax,phi=0.0,theta=0.0,psi=0.0)
data4 = ProjectGadgetDMCython2(snapfile,data4,parttype,zmin,zmax,phi=0.0,theta=0.0,psi=0.0)


print data1.data.max(), data2.data.max(), data3.data.max()
print data1.data.sum(), data2.data.sum(), data3.data.sum()
"""
for i in range(128):
	for j in range(128):
		if rand() < .01 and data1.data[i,j]>1.0:
			print 'Python = %.3g, CTypes = %.3g, Cython = %.3g\n'%(data1.data[i,j],data2.data[i,j],data3.data[i,j])
"""
