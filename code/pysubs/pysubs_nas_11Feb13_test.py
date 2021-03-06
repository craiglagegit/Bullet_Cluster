#!/usr/bin/env python
#Author: Craig Lage, NYU; 
#Date: 14-Sep-12

'''
These subroutines perform various reading, writing, and plotting functions for Gadget 3 and Enzo files for simulating the Bullet cluster.

'''
import lageconfig # system specific path information
from pylab import *
from subprocess import *
from ctypes import *
import time
import array
from scipy.ndimage import gaussian_filter, convolve
from scipy.ndimage.filters import sobel
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d.axes3d import Axes3D
from scipy.special import exp1 # Exponential integral
import h5py
import sys
import pyfits

import BulletConstants # Constants used in Bullet cluster simulations
from yt.mods import * # set up our namespace
from yt.visualization.volume_rendering.api import ProjectionCamera

# C libraries:
projectlib=CDLL(lageconfig.bulletpath+'HsmlAndProject.so')
projectnlib=CDLL(lageconfig.bulletpath+'HsmlAndProjectN.so')
projectdmlib=CDLL(lageconfig.bulletpath+'Project_DM.so')
gridlib=CDLL(lageconfig.bulletpath+'HsmlAndGrid.so')
gridinitlib=CDLL(lageconfig.bulletpath+'Grid_Init.so')
shiftnlib=CDLL(lageconfig.bulletpath+'Shift_n_Mask.so') 
profilelib=CDLL(lageconfig.bulletpath+'Profile.so')

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

class CArray2d(Structure):
    # This is for calling c subroutines
    _fields_ = [("nx", c_int),("ny", c_int),("xmin",c_double),("xmax",c_double),("ymin",c_double),("ymax",c_double),("dx",c_double),("dy",c_double),("x",POINTER(c_double)),("y",POINTER(c_double)),("data",POINTER(c_double))]

    def __init__(self,array):
	# Copies from Python to c
	self.nx=array.nx
	self.ny=array.ny
	self.dx=array.dx
	self.dy=array.dy
	self.xmin=array.xmin
	self.xmax=array.xmax
	self.ymin=array.ymin
	self.ymax=array.ymax
	size=array.nx*array.ny  
	self.data=(c_double*size)()
	self.x=(c_double*array.nx)()
	self.y=(c_double*array.ny)()

	for j in range(array.ny):
		self.y[j]=array.y[j]

	for i in range(array.nx):
		self.x[i]=array.x[i]
		for j in range(array.ny):
			self.data[i+j*array.nx]=array.data[i,j]

class ArraySet(Structure):
    # This carries a set of c arrays for shipping to the alignment subroutines
    _fields_ = [("numarrays", c_int), ("data1",POINTER(CArray2d)), ("data2",POINTER(CArray2d)), ("shifteddata1",POINTER(CArray2d)), ("sigma",POINTER(CArray2d)), ("mask",POINTER(CArray2d))]
    def __init__(self,numarrays,data1list, data2list, shifteddata1list, sigmalist, masklist):
	self.numarrays=numarrays
	self.data1=(CArray2d*numarrays)()
	self.data2=(CArray2d*numarrays)()
	self.shifteddata1=(CArray2d*numarrays)()
	self.sigma=(CArray2d*numarrays)()
	self.mask=(CArray2d*numarrays)()
	for i in range(numarrays):
		self.data1[i]=CArray2d(data1list[i])
		self.data2[i]=CArray2d(data2list[i])
		self.shifteddata1[i]=CArray2d(shifteddata1list[i])
		self.sigma[i]=CArray2d(sigmalist[i])
		self.mask[i]=CArray2d(masklist[i])

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

class ParticleData:
    def __init__(self,Ndm,Ngas):
	Npart = Ndm + Ngas
	self.id=zeros([Npart],dtype=int)
        self.x=zeros([Npart])
        self.y=zeros([Npart])
        self.z=zeros([Npart])
        self.vx=zeros([Npart])
        self.vy=zeros([Npart])
        self.vz=zeros([Npart])
        self.m=zeros([Npart])
        self.rho=zeros([Ngas])
        self.u=zeros([Ngas])
        self.ne=zeros([Ngas])
        self.np=zeros([Ngas])
        self.hsml=zeros([Ngas])

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

class CAlign(Structure):
    # This is for calling c subroutines
    _fields_ = [("d", c_double*3),("dmin", c_double*3),("dmax", c_double*3)]
    def __init__(self,align):
	# Copies from Python to c
	for i in range(3):
		self.d[i]=align.d[i]
		self.dmin[i]=align.dmin[i]
		self.dmax[i]=align.dmax[i]

def ReadGadgetSnapshot(filename):# Reads a snapshot file

    class Snapfile:
        pass

    snap = Snapfile()
    file = open(filename,"rb")
    snap.header = ReadGadgetHeader(file) # Header data
    snap.data = ReadGadgetData(file, snap.header,0) # Particle data, 0 = snapshot file, 1-IC file
    file.close()

    return snap

def ReadGadgetICFile(filename):# Reads an initial conditions file

    class ICfile:
        pass

    icfile = ICfile()
    file = open(filename,"rb")
    icfile.header = ReadGadgetHeader(file) # Header data
    icfile.data = ReadGadgetData(file, icfile.header,1) 
    file.close()

    return icfile

def WriteGadgetICFile(icfile,filename):# Writes an initial conditions file

    file = open(filename,"wb")
    WriteGadgetHeader(file, icfile) # Header data
    WriteGadgetData(file, icfile) # Particle data
    file.close()

def ShiftGadgetICFile(filename,DeltaV):# Shifts a Gadget ic file by a velocity increment

    icfile = ReadGadgetICFile(filename)
    Npart=icfile.header.Npart
    N=0
    for i in range(6):
        N = N+Npart[i]
    for n in range(N):
	for m in range(3):
		icfile.data.vel[3*n+m] = icfile.data.vel[3*n+m] + DeltaV[m]

    WriteGadgetICFile(icfile,filename)	

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

def CreateGadgetHeader(Ndm, Ngas): 

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

    header.blocksize.append(256)
    header.Npart.append(Ngas)
    header.Npart.append(Ndm)
    header.Npart.append(0)
    header.Npart.append(0)
    header.Npart.append(0)
    header.Npart.append(0)
    header.Massarr.append(0)
    header.Massarr.append(0)
    header.Massarr.append(0)
    header.Massarr.append(0)
    header.Massarr.append(0)
    header.Massarr.append(0)
    header.Time.append(0.0)
    header.Redshift.append(0.0)
    header.FlagSfr.append(0)
    header.FlagFeedback.append(0)
    header.Nall.append(Ngas)
    header.Nall.append(Ndm)
    header.Nall.append(0)
    header.Nall.append(0)
    header.Nall.append(0)
    header.Nall.append(0)
    header.FlagCooling.append(0)
    header.NumFiles.append(1)
    header.BoxSize.append(0.0)
    header.Omega0.append(0.0)
    header.OmegaLambda.append(0.0)
    header.HubbleParam.append(0.0)
    header.FlagAge.append(0)
    header.FlagMetals.append(0)
    header.NallHW.append(0)
    header.NallHW.append(0)
    header.NallHW.append(0)
    header.NallHW.append(0)
    header.NallHW.append(0)
    header.NallHW.append(0)
    header.flag_entr_ics.append(0)
    for i in range(15):
    	header.unused.append(0)
    header.blocksize.append(256)

    return header

def CreateGadgetData(P,header,filetype): 
    # Creates Gadget data, given a header and an instance of ParticleData

    class Data:
        pass
    data=Data()

    sizeofint = sizeoffloat = 4
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
    data.id = array.array("i")
    data.masses = array.array("f")
    data.u = array.array("f")
    data.rho = array.array("f")
    if header.FlagCooling[0]==1:
	    data.ne = array.array("f")
	    data.np = array.array("f")
    data.hsml = array.array("f")

    data.blocksize.append(N * 3 * sizeoffloat)
    for i in range(N):
    	data.pos.append(P.x[i])
    	data.pos.append(P.y[i])
    	data.pos.append(P.z[i])
    data.blocksize.append(N * 3 * sizeoffloat)

    data.blocksize.append(N * 3 * sizeoffloat)
    for i in range(N):
    	data.vel.append(P.vx[i])
    	data.vel.append(P.vy[i])
    	data.vel.append(P.vz[i])
    data.blocksize.append(N * 3 * sizeoffloat)

    data.blocksize.append(N * sizeofint)
    for i in range(N):
    	data.id.append(P.id[i])
    data.blocksize.append(N * sizeofint)

    if Nm!=0:
    	data.blocksize.append(Nm * sizeoffloat)
    	for i in range(Nm):
    		data.masses.append(P.m[i])
    	data.blocksize.append(Nm * sizeoffloat)
    if Ngas!=0:      
    	data.blocksize.append(Ngas * sizeoffloat)
    	for i in range(Ngas):
    		data.u.append(P.u[i])
    	data.blocksize.append(Ngas * sizeoffloat)

   	if filetype==0:

	        data.blocksize.append(Ngas * sizeoffloat)
	        for i in range(Ngas):
    			data.rho.append(P.rho[i])
	        data.blocksize.append(Ngas * sizeoffloat)

		if header.FlagCooling[0]==1:

	        	data.blocksize.append(Ngas * sizeoffloat)
	        	for i in range(Ngas):
    				data.ne.append(P.ne[i])
	        	data.blocksize.append(Ngas * sizeoffloat)

	        	data.blocksize.append(Ngas * 3 * sizeoffloat)
	        	for i in range(Ngas):
    				data.np.append(P.np[i])
	        	data.blocksize.append(Ngas * 3 * sizeoffloat)

        	data.blocksize.append(Ngas * 3 * sizeoffloat)
	        for i in range(Ngas):
    			data.hsml.append(P.hsml[i])
        	data.blocksize.append(Ngas * 3 * sizeoffloat)

    return data


def WriteGadgetHeader(file, icfile): # Writes the header to the IC file

    blocksize = array.array("i")
    blocksize.append(icfile.header.blocksize[0])
    blocksize.tofile(file)
    icfile.header.Npart.tofile(file)
    icfile.header.Massarr.tofile(file)
    icfile.header.Time.tofile(file)
    icfile.header.Redshift.tofile(file)
    icfile.header.FlagSfr.tofile(file)
    icfile.header.FlagFeedback.tofile(file)
    icfile.header.Nall.tofile(file)
    icfile.header.FlagCooling.tofile(file)
    icfile.header.NumFiles.tofile(file)
    icfile.header.BoxSize.tofile(file)
    icfile.header.Omega0.tofile(file)
    icfile.header.OmegaLambda.tofile(file)
    icfile.header.HubbleParam.tofile(file)
    icfile.header.FlagAge.tofile(file)
    icfile.header.FlagMetals.tofile(file)
    icfile.header.NallHW.tofile(file)
    icfile.header.flag_entr_ics.tofile(file)
    icfile.header.unused.tofile(file)
    blocksize.tofile(file)

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
    data.id = array.array("i")
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

def WriteGadgetData(file, icfile): # Writes the data to the IC file

    Npart=icfile.header.Npart
    Massarr=icfile.header.Massarr
    Ngas=Npart[0]
    Nm=0
    for i in range(6):
        if Massarr[i]==0:
            Nm = Nm+Npart[i]

    blocksize = array.array("i")
    blocksize.append(icfile.data.blocksize[0])
    blocksize.tofile(file)
    icfile.data.pos.tofile(file)
    blocksize.tofile(file)

    blocksize.pop()
    blocksize.append(icfile.data.blocksize[2])
    blocksize.tofile(file)
    icfile.data.vel.tofile(file)
    blocksize.tofile(file)

    blocksize.pop()
    blocksize.append(icfile.data.blocksize[4])
    blocksize.tofile(file)
    icfile.data.id.tofile(file)
    blocksize.tofile(file)

    if Nm!=0:
    	blocksize.pop()
    	blocksize.append(icfile.data.blocksize[6])
    	blocksize.tofile(file)
	print "Number of masses = %d\n"%len(icfile.data.masses)
    	icfile.data.masses.tofile(file)
    	blocksize.tofile(file)

    if Ngas!=0:
    	blocksize.pop()
    	blocksize.append(icfile.data.blocksize[8])
    	blocksize.tofile(file)
    	icfile.data.u.tofile(file)
    	blocksize.tofile(file)

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

def PrintGadgetData(snapfile,parttype,nmin,nmax,filetype): 
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

def ReadEnzoParameterFile(filename):
    # This reads an Enzo parameter file and returns a set of lists with the parameters
    class ParamData:
        pass
    param=ParamData()

    try:
        file = open(filename,'r')
        lines=file.readlines()
        file.close()
        for line in lines:
		if line[0]=='#' or line[0]=='\n':
			continue
		ThisLine=line.strip().split()
		ThisLineLength=len(ThisLine)
		#print ThisLine
		ParamName=ThisLine[0]
		HasIndex=False
		if ParamName.find('[')!=-1:
			ParamIndex=int(ParamName.split('[')[1].rstrip(']'))
			ParamName=ParamName.split('[')[0]
			HasIndex=True
		#print ParamName
		ThisParam=list()
		for j in range(2,ThisLineLength):
			try:
				value=int(ThisLine[j])
			except ValueError:
				try:
					value=float(ThisLine[j])
				except ValueError:
					value=ThisLine[j]
			ThisParam.append(value)
		if len(ThisParam)==1:
			ThisParam = ThisParam[0]
		#print ThisParam
		if HasIndex:
			if ParamIndex==1:
				exec("param.%s=list()" %ParamName)
				exec("param.%s.append('0')"%ParamName)
				exec("param.%s.append(ThisParam)" %ParamName)
			else:				
				exec("param.%s.append(ThisParam)" %ParamName)
		else:
			exec("param.%s=ThisParam" %ParamName)
    except IOError:
        print "Error reading Enzo parameter file"

    return param

def WriteEnzoInitFile(data,name,filename):
    # Writes an initial conditions file
    outfile = h5py.File(filename,"w")
    if len(data.shape)==4: # This loop puts data in fortran order before writing to the HDF5 file
	newdata=zeros([1,data.size])
	nx = data.shape[1]
	ny = data.shape[2]
	nz = data.shape[3]
	for i in range(nx):
		for j in range(ny):
			for k in range(nz):
				newdata[0,i + nx*j + nx*ny*k] = data[0,i,j,k]
    	dataset = outfile.create_dataset(name,newdata.shape,'>f8')
    	dataset.write_direct(newdata)
    else: 
	dataset = outfile.create_dataset(name,data.shape,'>f8')
	dataset.write_direct(data)
    dataset.attrs["Component_Rank"]=data.shape[0]
    dataset.attrs["Component_Size"]=data.size/data.shape[0]
    dataset.attrs["Rank"]=len(data.shape)-1 
    shapelist=list(data.shape)
    shapelist.remove(shapelist[0])
    endlist=list(shapelist)
    for i in range(len(endlist)):
	endlist[i]=endlist[i]-1
    dataset.attrs["Dimensions"]=tuple(shapelist)
    dataset.attrs["TopGridStart"]=(0,0,0)	        # Don't think these are used at present
    dataset.attrs["TopGridEnd"]=tuple(endlist)	        # Don't think these are used at present
    dataset.attrs["TopGridDims"]=tuple(shapelist)	# Don't think these are used at present
    outfile.close()

def ReadEnzoData(filename): 
    # Reads the data from a single Enzo output file and returns a Python object with the data in NumPy arrays.

    class Data:
        pass
    data=Data()

    names=list()
    outfile=h5py.File(filename,'r')
    numparams=len(outfile.items()[0][1])
    
    for i in range(numparams):
	name=outfile.items()[0][1].items()[i][0]
	name=name.replace("-","_")	
	#print name
    	dataset=outfile.items()[0][1].items()[i][1]
    	temparray=zeros(dataset.shape)
	dataset.read_direct(temparray)
    	exec("data.%s=temparray" %name)

    outfile.close()
    return data

def PrintEnzoData(data,nmin,nmax,grid): 
    # Prints selected output file data.  nmin and nmax are the particle indices.  grid is a 3D array telling what part of the grid should be printed.

    for i in range(nmin,nmax):
        print "Particle ",i," index = ",data.particle_index[i]
        print "Particle ",i," mass = ",data.particle_mass[i]
        print "Particle ",i," type = ",data.particle_type[i]
        print "Particle ",i," X = ",data.particle_position_x[i]
        print "Particle ",i," Y = ",data.particle_position_y[i]
        print "Particle ",i," Z = ",data.particle_position_z[i]
        print "Particle ",i," Vx = ",data.particle_velocity_x[i]
        print "Particle ",i," Vy = ",data.particle_velocity_y[i]
        print "Particle ",i," Vz = ",data.particle_velocity_z[i]

    for i in range(grid.nx):
	for j in range(grid.ny):
	    for k in range(grid.nz):
        	print "Grid point ",i,j,k,"Dark_Matter_Density = ",data.Dark_Matter_Density[i,j,k]
        	print "Grid point ",i,j,k,"Density = ",data.Density[i,j,k]
        	print "Grid point ",i,j,k,"GasEnergy = ",data.GasEnergy[i,j,k]
        	print "Grid point ",i,j,k,"Temperature = ",data.Temperature[i,j,k]
        	print "Grid point ",i,j,k,"TotalEnergy = ",data.TotalEnergy[i,j,k]
        	print "Grid point ",i,j,k,"X-Velocity = ",data.x_velocity[i,j,k]
        	print "Grid point ",i,j,k,"Y-Velocity = ",data.y_velocity[i,j,k]
        	print "Grid point ",i,j,k,"Z-Velocity = ",data.z_velocity[i,j,k]

def GadgetKernel(r,h): 
    # This is the spline kernel used in Gadget2.  
    # Currently not used, but I might need it.
    m=8/(pi*h**3)
    rh = r/h
    if rh>=0 and rh<=1/2:
        return m(1-6*rh**2+6*rh**3)
    elif rh>1/2 and rh<=1.0:
        return 2*m*(1-rh)**3
    else:
        return 0.0

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

def AddGadgetGalaxies(infile,outfile,N,Mass):
    # Adds a number of galaxies to a gadget IC file
    # Converts DM particles (type 1) to galaxies (type 4)

    icfile = ReadGadgetICFile(infile)
    icfile.header.Npart[1] = icfile.header.Npart[1] - N
    icfile.header.Npart[4] = N
    icfile.header.Nall[1] = icfile.header.Nall[1] - N
    icfile.header.Nall[4] = N
    Nstart = icfile.header.Npart[0] + icfile.header.Npart[1]
    for i in range(Nstart,Nstart+N):
	icfile.data.masses[i] = Mass

    WriteGadgetICFile(icfile,outfile)	

def CreateGadgetMarkerParticles(infile,outfile,Ng,MassGas,Nd,MassDM):
    # Increases the mass of a number of DM particles
    # to serve as markers for testing the projection algorithms

    icfile = ReadGadgetICFile(infile)
    Ngas = icfile.header.Npart[0] 
    NDM = icfile.header.Npart[1]
    for i in range(Ng):
	part = int(Ngas*rand())
	icfile.data.masses[part] = icfile.data.masses[part] + MassGas
	print 'XGas = %.3f, YGas = %.3f, ZGas = %.3f\n'%(icfile.data.pos[3*part],icfile.data.pos[3*part+1],icfile.data.pos[3*part+2])
    for i in range(Nd):
	part = Ngas + int(NDM*rand())
	icfile.data.masses[part] = icfile.data.masses[part] + MassDM
	print 'XDM = %.3f, YDM = %.3f, ZDM = %.3f\n'%(icfile.data.pos[3*part],icfile.data.pos[3*part+1],icfile.data.pos[3*part+2])
    WriteGadgetICFile(icfile,outfile)	


def AddGadgetCDGalaxy(infile,outfile,Mass,r):
    # Adds a single CD galaxy to a gadget IC file
    # Places one galaxy (type 4) at rest near the center of the DM cloud
    # Will be within radius r

    icfile = ReadGadgetICFile(infile)
    icfile.header.Npart[4] = 1
    icfile.header.Nall[4] = 1
    NCD = icfile.header.Npart[0] + icfile.header.Npart[1]
    icfile.data.masses.append(Mass)
    icfile.data.id.append(NCD + 1)
    for i in range(3):
	icfile.data.vel.append(0.0)
	icfile.data.pos.append(-r + 2 * r * rand())

    WriteGadgetICFile(icfile,outfile)	


def ScatterPlot3d(snapfile,ax,parttype,size,col): 
    # Makes a 3D scatter plot
    # Need to clean this up to use Array classes
    Npart=snapfile.header.Npart
    num=Npart[parttype]

    offset=0
    for k in range(parttype):
        offset=offset+Npart[k]

    x=zeros([num])
    y=zeros([num])
    z=zeros([num])

    for j in range(num):

        jo=offset+j
        x[j] = snapfile.data.pos[3*jo]
        y[j] = snapfile.data.pos[3*jo+1]
        z[j] = snapfile.data.pos[3*jo+2]

    ax.plot(x,y,z,marker='.', markersize=size, linewidth=0, c=col)



def GridGadgetGas(snapfile,mass,phi=0.0,theta=0.0,psi=0.0): 
    # Projects 3D gas particle data onto a 3D grid
    # Returns mass data, Temp, Pressure, XRay intensity 
    # Euler angles are used to rotate the data if desired
    xpixels=mass.nx
    ypixels=mass.ny
    zpixels=mass.nz
    xyz=xpixels*ypixels*zpixels
    PixelVolume=mass.dx * mass.dy * mass.dz
    PixelArea=mass.dx * mass.dy
    R=EulerAngles(phi,theta,psi)
    x=zeros([3]) # Position vector
    xp=zeros([3])# Rotated position vector
    # Mass will go in the input 3d array. 
    # Will create additional 3d arrays for the Temp, Rho, Xray, Pressure, and Velocity data
    Temp=Array3d(mass.xmin,mass.xmax,mass.nx,mass.ymin,mass.ymax,mass.ny,mass.zmin,mass.zmax,mass.nz)
    Rho=Array3d(mass.xmin,mass.xmax,mass.nx,mass.ymin,mass.ymax,mass.ny,mass.zmin,mass.zmax,mass.nz)
    Xray=Array3d(mass.xmin,mass.xmax,mass.nx,mass.ymin,mass.ymax,mass.ny,mass.zmin,mass.zmax,mass.nz)
    Pressure=Array3d(mass.xmin,mass.xmax,mass.nx,mass.ymin,mass.ymax,mass.ny,mass.zmin,mass.zmax,mass.nz)
    SZ=Array3d(mass.xmin,mass.xmax,mass.nx,mass.ymin,mass.ymax,mass.ny,mass.zmin,mass.zmax,mass.nz)

    Npart=snapfile.header.Npart
    Ngas=Npart[0] # Gas particles

    class P(Structure):
        _fields_ = [("x", c_float),("y", c_float),("z", c_float)]

    Pptr = (P*Ngas)()  # Pointer to an array of Pos values
    massptr = (c_float*Ngas)()  # Pointer to an array of input mass values
    hsmlptr = (c_float*Ngas)()  # Pointer to an array of hsml vlues
    quantityptr = (c_float*Ngas)()  # Pointer to an array of input values
    valueptr = (c_float*xyz)() # Pointer to an array of output mass values
    valuequantityptr = (c_float*xyz)() # Pointer to an array of output values

    for j in range(Ngas):
        massptr[j] = snapfile.data.masses[j]
        hsmlptr[j] = snapfile.data.hsml[j]
	quantityptr[j] = snapfile.data.u[j] * BulletConstants.Tconv # Tconv converts u into T in eV.

        for m in range(3):
            x[m]=snapfile.data.pos[3*j+m]

        xp=dot(R,x) # Rotate the position by the Euler angles
     	Pptr[j].x = xp[0]
     	Pptr[j].y = xp[1]
   	Pptr[j].z = xp[2]

    xmin=c_float(mass.xmin)
    xmax=c_float(mass.xmax)
    ymin=c_float(mass.ymin)
    ymax=c_float(mass.ymax)
    zmin=c_float(mass.zmin)
    zmax=c_float(mass.zmax)

    desdensngb = BulletConstants.Nsph # SPH number of neighbors

    axis1=0
    axis2=1
    axis3=2 # This is the axis which is projected out
    # I always project out the Z-axis and use Euler angles to rotate

    hmax= c_float(BulletConstants.HsmlMax) # Max Hsml size
    boxsize = c_double(snapfile.header.BoxSize[0]) # used for Periodic boundary conditions

    start = time.time()

    gridlib.findHsmlAndGrid(Ngas, Pptr, hsmlptr, massptr, quantityptr, xmin, xmax, ymin, ymax, zmin, zmax, xpixels, ypixels, zpixels, desdensngb, axis1, axis2, axis3, hmax, boxsize, valueptr, valuequantityptr)

    XFactor = BulletConstants.PreFactor * BulletConstants.TimeFactor * BulletConstants.AreaFactor  * BulletConstants.nconv**2 / PixelVolume
    RhoMin = 1.02 * 10**(BulletConstants.LogRhoMin)
    XrayMin = 1.02 * 10**(BulletConstants.LogXrayMin)
    PressureMin = 1.02 * 10**(BulletConstants.LogPressureMin)
    for i in range(xpixels):
    	for j in range(ypixels):
	    for k in range(zpixels):
                mass.data[i,j,k] = m = max(1E-24,valueptr[i+j*xpixels+k*xpixels*ypixels])
                Temp.data[i,j,k] = T = max(1E-12,valuequantityptr[i+j*xpixels+k*xpixels*ypixels])
                Rho.data[i,j,k] = rho = max(RhoMin, m / PixelVolume)
                Xray.data[i,j,k] = max(XrayMin, XFactor * m**2 * sqrt(BulletConstants.me/T) * (exp1(BulletConstants.Emin/T)-exp1(BulletConstants.Emax/T)))
                Pressure.data[i,j,k] = max(PressureMin, rho * T) # Not absolutely calibrated
		SZ.data[i,j,k] = BulletConstants.SZFactor * m * T / PixelArea

    elapsed = (time.time()-start)
    print "Elapsed time to run Gadget gas gridding = "+str(elapsed)

    return [mass,Temp,Rho,Xray,Pressure,SZ] 

def FindGadgetCentroids(snapfile): 
    # Finds centroids of dark matter and gas particles
    Npart=snapfile.header.Npart
    Masses=zeros([2,2])
    NumPart=zeros([2,2])
    Centroids=zeros([2,2,3])

    start = time.time()
    for parttype in range(2): # Finding the 4 different masses
	if parttype == 0:
		offset = 0
	else:
		offset = Npart[0]
	Masses[parttype,0] = snapfile.data.masses[offset]
	print 'Mass of particle type %d ; cluster 0 = %.4f'%(parttype,Masses[parttype,0])
	while (Masses[parttype,1] == 0):
		pc = offset + int(Npart[parttype] * rand()) # Keep randomly choosing particles until I have found 4 different masses
		print 'Trying another mass = %.3f'%snapfile.data.masses[pc]
		if abs(snapfile.data.masses[pc] - Masses[parttype,0]) > 1E-8:
			Masses[parttype,1] = snapfile.data.masses[pc]
			print 'Mass of particle type %d ; cluster 1 = %.4f'%(parttype,Masses[parttype,1])
    for parttype in range(2): 
	if parttype == 0:
		offset = 0
	else:
		offset = Npart[0]
	for n in range(offset , offset + Npart[parttype]): # n cycles through the number of particles of this type
		#if n%100000 == 0:
			#print '*'
		if abs(snapfile.data.masses[n] - Masses[parttype,0]) < 1E-8:	
			NumPart[parttype,0] = NumPart[parttype,0] + 1
			for m in range(3):
				Centroids[parttype,0,m] = Centroids[parttype,0,m] + snapfile.data.pos[3*n+m]
		else:	
			NumPart[parttype,1] = NumPart[parttype,1] + 1
			for m in range(3):
				Centroids[parttype,1,m] = Centroids[parttype,1,m] + snapfile.data.pos[3*n+m]

    for parttype in range(2):
	    for k in range(2): # k denotes the bullet or main particles
			if NumPart[parttype,0] > NumPart[parttype,1]: # Swap everything to put the bullet particles first
				TempNum = NumPart[parttype,0]
				NumPart[parttype,0] = NumPart[parttype,1]
				NumPart[parttype,1] = TempNum
				TempMass = Masses[parttype,0]
				Masses[parttype,0] = Masses[parttype,1]
				Masses[parttype,1] = TempMass
				for m in range(3):
					TempCentroid = Centroids[parttype,0,m]
					Centroids[parttype,0,m] = Centroids[parttype,1,m]
					Centroids[parttype,1,m] = TempCentroid

			for m in range(3):
				Centroids[parttype,k,m] = Centroids[parttype,k,m] / NumPart[parttype,k]

    elapsed = (time.time()-start)
    print "Elapsed time to locate centroids = "+str(elapsed)

    return [NumPart, Masses, Centroids]

def ProjectEnzoData(pf,mass,phi=0.0,theta=0.0,psi=0.0,zmin=-3000.0,zmax=3000.0,DMProject=False):
    # Projects 3D Enzo data onto a 2D grid
    # Returns DM mass, Baryon mass, Three Xray intensities, and SZ data.
    # DMProject = True runs a separate DM Projection.
    # DMProject = False uses the yt raycasting.

    #print 'PATH = ',os.environ['PATH']
    #print 'PYTHONPATH = ',os.environ['PYTHONPATH']
    #print 'LD_LIBRARY_PATH = ',os.environ['LD_LIBRARY_PATH']
    #envir = Popen('env',shell=True, stdout = PIPE).communicate()[0].split('\n')
    #print envir

    DM = Array2d(mass.xmin,mass.xmax,mass.nx,mass.ymin,mass.ymax,mass.ny) 
    if DMProject:
	print "Doing Separate DM Projection"
    	DM = ProjectEnzoDM(pf,mass,1,zmin,zmax,-psi,-theta,-phi)
    else:
	print "Skipping Separate DM Projection"
    start = time.time()
    xpixels=mass.nx
    ypixels=mass.ny
    PixelArea = mass.dx * mass.dy
    # Mass density will go in the input 2d array. 
    # Will create additional 2d arrays for the Xray and SZ data
    Xray1=Array2d(mass.xmin,mass.xmax,mass.nx,mass.ymin,mass.ymax,mass.ny)
    Xray2=Array2d(mass.xmin,mass.xmax,mass.nx,mass.ymin,mass.ymax,mass.ny)
    Xray3=Array2d(mass.xmin,mass.xmax,mass.nx,mass.ymin,mass.ymax,mass.ny)
    SZ=Array2d(mass.xmin,mass.xmax,mass.nx,mass.ymin,mass.ymax,mass.ny)
    add_field('XRay1', function = _EnzoXRay1, units=r"\rm{cm}^{-3}\rm{s}^{-1}", projected_units=r"\rm{cm}^{-2}\rm{s}^{-1}", validators=[ValidateParameter('ApecData')], take_log=False)
    add_field('XRay2', function = _EnzoXRay2, units=r"\rm{cm}^{-3}\rm{s}^{-1}", projected_units=r"\rm{cm}^{-2}\rm{s}^{-1}", validators=[ValidateParameter('ApecData')], take_log=False)
    add_field('XRay3', function = _EnzoXRay3, units=r"\rm{cm}^{-3}\rm{s}^{-1}", projected_units=r"\rm{cm}^{-2}\rm{s}^{-1}", validators=[ValidateParameter('ApecData')], take_log=False)
    add_field('SZ', function =   _EnzoSZ, units =r"\rm{K}\rm{cm}^{-1}" , projected_units=r"\rm{K}",take_log=False)

    pf.field_info['Density'].take_log=False
    pf.field_info['XRay1'].take_log=False
    pf.field_info['XRay2'].take_log=False
    pf.field_info['XRay3'].take_log=False
    pf.field_info['SZ'].take_log=False
    pf.field_info['Dark_Matter_Density'].take_log=False

    center = [(mass.xmin+mass.xmax)/2.0,(mass.ymin+mass.ymax)/2.0,(zmin+zmax)/2.0] # Data Center
    normal_vector=(0.0,0.0,1.0)
    north_vector = (0.0,1.0,0.0)

    R = EulerAngles(phi,theta,psi)
    normal_vector = dot(R,normal_vector)
    north_vector = dot(R,north_vector)
    width = (mass.xmax - mass.xmin, mass.ymax - mass.ymin, zmax - zmin)
    resolution = (mass.nx,mass.ny)

    MassFactor = BulletConstants.cm_per_kpc**2 * PixelArea / (BulletConstants.g_per_Msun * 1E10)
    XFactor    = BulletConstants.cm_per_kpc**2 * PixelArea * BulletConstants.AreaFactor

    if not DMProject:
	print "Using yt Ray casting for DM"
        projcam=ProjectionCamera(center,normal_vector,width,resolution,"Dark_Matter_Density",north_vector=north_vector,pf=pf, interpolated=True)
        DM.data =  projcam.snapshot()[:,:,0] * MassFactor 

    projcam=ProjectionCamera(center,normal_vector,width,resolution,"Density",north_vector=north_vector,pf=pf, interpolated=True)
    mass.data =  projcam.snapshot()[:,:,0] * MassFactor
    projcam=ProjectionCamera(center,normal_vector,width,resolution,"XRay1",north_vector=north_vector,pf=pf, interpolated=True)
    Xray1.data = projcam.snapshot()[:,:,0] * XFactor
    projcam=ProjectionCamera(center,normal_vector,width,resolution,"XRay2",north_vector=north_vector,pf=pf, interpolated=True)
    Xray2.data = projcam.snapshot()[:,:,0] * XFactor
    projcam=ProjectionCamera(center,normal_vector,width,resolution,"XRay3",north_vector=north_vector,pf=pf, interpolated=True)
    Xray3.data = projcam.snapshot()[:,:,0] * XFactor
    projcam=ProjectionCamera(center,normal_vector,width,resolution,"SZ",north_vector=north_vector,pf=pf, interpolated=True)
    SZ.data   = projcam.snapshot()[:,:,0]
    elapsed = (time.time()-start)
    print "Elapsed time to run Enzo gas projection = "+str(elapsed)

    return [DM, mass,Xray1,Xray2,Xray3,SZ]

def ProjectEnzoTemp(pf,mass,phi=0.0,theta=0.0,psi=0.0,zmin=-3000.0,zmax=3000.0):
    # Projects 3D Enzo data onto a 2D grid
    # Returns mass data, Three Xray intensities, and SZ data.

    start = time.time()
    xpixels=mass.nx
    ypixels=mass.ny
    PixelArea = mass.dx * mass.dy
    Temp=Array2d(mass.xmin,mass.xmax,mass.nx,mass.ymin,mass.ymax,mass.ny)
    SZ=Array2d(mass.xmin,mass.xmax,mass.nx,mass.ymin,mass.ymax,mass.ny)
    BMag=Array2d(mass.xmin,mass.xmax,mass.nx,mass.ymin,mass.ymax,mass.ny)
    add_field('TXRay', function = _EnzoTXRay, units=r"\rm{cm}^{-3}\rm{s}^{-1}", projected_units=r"\rm{cm}^{-2}\rm{s}^{-1}", validators=[ValidateParameter('ApecData')], take_log=False)
    add_field('XRay', function = _EnzoXRay, units=r"\rm{cm}^{-3}\rm{s}^{-1}", projected_units=r"\rm{cm}^{-2}\rm{s}^{-1}", validators=[ValidateParameter('ApecData')], take_log=False)
    add_field('BMag', function = _EnzoBMag,take_log=False)

    center = [(mass.xmin+mass.xmax)/2.0,(mass.ymin+mass.ymax)/2.0,(zmin+zmax)/2.0] # Data Center
    normal_vector=(0.0,0.0,1.0)
    north_vector = (0.0,1.0,0.0)
    R = EulerAngles(phi,theta,psi)
    normal_vector = dot(R,normal_vector)
    north_vector = dot(R,north_vector)
    width = (mass.xmax - mass.xmin, mass.ymax - mass.ymin, zmax - zmin)
    resolution = (mass.nx,mass.ny)

    ave_weight = "CellMassMsun" # The weight for the average
    sp = pf.h.sphere( (0.0,0.0,0.0), (10.0,"kpc"))
    average_temp = sp.quantities["WeightedAverageQuantity"]("Temperature", ave_weight)
    #average_xray = sp.quantities["WeightedAverageQuantity"]("XRay", ave_weight)
    #average_txray = sp.quantities["WeightedAverageQuantity"]("TXRay", ave_weight)

    print "Average Temp = %f keV"%average_temp
    #print "Average XRay = %g "%average_xray
    #print "Average TXRay = %g "%average_txray




    MassFactor = BulletConstants.cm_per_kpc**2 * PixelArea / (BulletConstants.g_per_Msun * 1E10)
    XFactor    = BulletConstants.cm_per_kpc**2 * PixelArea * BulletConstants.AreaFactor
    BFactor    = 1.0 / (BulletConstants.cm_per_kpc * width[2]) 

    projcam=ProjectionCamera(center,normal_vector,width,resolution,"TXRay",north_vector=north_vector,pf=pf, interpolated=True)
    Temp.data = projcam.snapshot()[:,:,0]
    projcam=ProjectionCamera(center,normal_vector,width,resolution,"XRay",north_vector=north_vector,pf=pf, interpolated=True)
    Temp.data = Temp.data / projcam.snapshot()[:,:,0]
    projcam=ProjectionCamera(center,normal_vector,width,resolution,"BMag",north_vector=north_vector,pf=pf, interpolated=True)
    BMag.data = projcam.snapshot()[:,:,0] * BFactor
    TrueTemp=Array2d(mass.xmin,mass.xmax,mass.nx,mass.ymin,mass.ymax,mass.ny)
    projcam=ProjectionCamera(center,normal_vector,width,resolution,"Temperature",north_vector=north_vector,pf=pf, interpolated=True)
    TrueTemp.data = projcam.snapshot()[:,:,0] * BulletConstants.DegreesK_to_keV / (BulletConstants.cm_per_kpc * width[2])
    elapsed = (time.time()-start)
    print "Elapsed time to run Enzo temp projection = "+str(elapsed)
    return [TrueTemp,Temp,BMag]


def _EnzoTXRay(field,data):
	ApecData = data.get_field_parameter('ApecData')
	TFudge = data.get_field_parameter('TFudge')
	T = data['Temperature'] * BulletConstants.DegreesK_to_keV * TFudge
	logT = np.log10(T) * 100.0
	np.putmask(logT,logT>179.0,179.0)
	np.putmask(logT,logT<-100.0,-100.0)
	minT = logT.astype(int)
	flux = (minT+1.0-logT) * ApecData[minT+100,3] + (logT-minT) * ApecData[minT+101,3]
	return T * flux * data['Density']**2 / BulletConstants.mp**2 * 1E-14

def _EnzoXRay(field,data):
	ApecData = data.get_field_parameter('ApecData')
	TFudge = data.get_field_parameter('TFudge')
	logT = np.log10(data['Temperature'] * BulletConstants.DegreesK_to_keV * TFudge) * 100.0
	np.putmask(logT,logT>179.0,179.0)
	np.putmask(logT,logT<-100.0,-100.0)
	minT = logT.astype(int)
	flux = (minT+1.0-logT) * ApecData[minT+100,3] + (logT-minT) * ApecData[minT+101,3]
	return flux * data['Density']**2 / BulletConstants.mp**2 * 1E-14

def _EnzoXRay1(field,data):
	ApecData = data.get_field_parameter('ApecData')
	TFudge = data.get_field_parameter('TFudge')
	logT = np.log10(data['Temperature'] * BulletConstants.DegreesK_to_keV * TFudge) * 100.0
	np.putmask(logT,logT>179.0,179.0)
	np.putmask(logT,logT<-100.0,-100.0)
	minT = logT.astype(int)
	flux = (minT+1.0-logT) * ApecData[minT+100,0] + (logT-minT) * ApecData[minT+101,0]
	return flux * data['Density']**2 / BulletConstants.mp**2 * 1E-14 

def _EnzoXRay2(field,data):
	ApecData = data.get_field_parameter('ApecData')
	TFudge = data.get_field_parameter('TFudge')
	logT = np.log10(data['Temperature'] * BulletConstants.DegreesK_to_keV * TFudge) * 100.0
	np.putmask(logT,logT>179.0,179.0)
	np.putmask(logT,logT<-100.0,-100.0)
	minT = logT.astype(int)
	flux = (minT+1.0-logT) * ApecData[minT+100,1] + (logT-minT) * ApecData[minT+101,1]
	return flux * data['Density']**2 / BulletConstants.mp**2 * 1E-14

def _EnzoXRay3(field,data):
	ApecData = data.get_field_parameter('ApecData')
	TFudge = data.get_field_parameter('TFudge')
	logT = np.log10(data['Temperature'] * BulletConstants.DegreesK_to_keV * TFudge) * 100.0
	np.putmask(logT,logT>179.0,179.0)
	np.putmask(logT,logT<-100.0,-100.0)
	minT = logT.astype(int)
	flux = (minT+1.0-logT) * ApecData[minT+100,2] + (logT-minT) * ApecData[minT+101,2]
	return flux * data['Density']**2 / BulletConstants.mp**2 * 1E-14

def _EnzoSZ(field,data):
	TFudge = data.get_field_parameter('TFudge')
	T = data['Temperature'] * BulletConstants.DegreesK_to_eV * TFudge
	SZ = BulletConstants.EnzoSZFactor * data['Density'] * T 
	return SZ

def _EnzoBMag(field, data):
    return np.sqrt(data["Bx"] * data["Bx"] + data["By"] * data["By"] + data["Bz"] * data["Bz"])

def ReadLookups(Z):
	# This subroutine reads in the data from the APEC lookup tables and places it in an array
	# The data is interpolated to get the data for the required Z (metallicity)
	ApecData = zeros([281,4]) # Array to hold the data
	for bin in range(3): #bin 0 is 0.5-2kev, bin 1 is 2-5kev, bin 2 is 5-8kev, bin 3 is 0.5-8kev
		if bin == 0:
			infile = open(lageconfig.toppath+'bullet/data/apec/apec_xray_0.5_2.0.txt','r')
		elif bin == 1:
			infile = open(lageconfig.toppath+'bullet/data/apec/apec_xray_2.0_5.0.txt','r')
		elif bin == 2:
			infile = open(lageconfig.toppath+'bullet/data/apec/apec_xray_5.0_8.0.txt','r')
		lines = infile.readlines()
		infile.close()
		counter = 0
		for line in lines:
			if line.strip().split()[0] == 'LogT': # Skips header line
				continue
			minZ = max(0, int(round((Z * 10))))
			maxZ = minZ + 1
			if maxZ > 10:
				maxZ = 10
				minZ = 9
			f = (maxZ-10.0*Z) * float(line.strip().split()[minZ+1]) + (10.0*Z-minZ) * float(line.strip().split()[maxZ+1])
			ApecData[counter,bin] = f
			if bin == 0:
				ApecData[counter,3] = f # This bin is just the sum of the other three
			else:
				ApecData[counter,3] = ApecData[counter,3] + f	
			counter=counter+1
	return ApecData

def GridGadgetGasInit(snapfile,mass,phi=0.0,theta=0.0,psi=0.0): 
    # Projects 3D gas particle data onto a 3D grid
    # Returns mass data, internal energy and velocity
    # This is used to translate Gadget initial data files into Enzo
    # Euler angles are used to rotate the data if desired

    xpixels=mass.nx
    ypixels=mass.ny
    zpixels=mass.nz
    xyz=xpixels*ypixels*zpixels
    PixelVolume=mass.dx * mass.dy * mass.dz
    R=EulerAngles(phi,theta,psi)
    x=zeros([3]) # Position vector
    xp=zeros([3])# Rotated position vector
    v=zeros([3]) # Position vector
    vp=zeros([3])# Rotated position vector

    # Mass will go in the input 3d array. 
    # Will create additional 3d arrays for the Density, Internal Energy, and Velocity data
    Density=Array3d(mass.xmin,mass.xmax,mass.nx,mass.ymin,mass.ymax,mass.ny,mass.zmin,mass.zmax,mass.nz)
    InternalEnergy=Array3d(mass.xmin,mass.xmax,mass.nx,mass.ymin,mass.ymax,mass.ny,mass.zmin,mass.zmax,mass.nz)
    Vx=Array3d(mass.xmin,mass.xmax,mass.nx,mass.ymin,mass.ymax,mass.ny,mass.zmin,mass.zmax,mass.nz)
    Vy=Array3d(mass.xmin,mass.xmax,mass.nx,mass.ymin,mass.ymax,mass.ny,mass.zmin,mass.zmax,mass.nz)
    Vz=Array3d(mass.xmin,mass.xmax,mass.nx,mass.ymin,mass.ymax,mass.ny,mass.zmin,mass.zmax,mass.nz)

    Npart=snapfile.header.Npart
    Ngas=Npart[0] # Gas particles

    class P(Structure):
        _fields_ = [("x", c_float),("y", c_float),("z", c_float)]

    Pptr = (P*Ngas)()  # Pointer to an array of Pos values
    massptr = (c_float*Ngas)()  # Pointer to an array of input mass values
    hsmlptr = (c_float*Ngas)()  # Pointer to an array of hsml vlues
    quantity1ptr = (c_float*Ngas)()  # Pointer to an array of input values - Vx
    quantity2ptr = (c_float*Ngas)()  # Pointer to an array of input values - Vy
    quantity3ptr = (c_float*Ngas)()  # Pointer to an array of input values - Vz
    quantity4ptr = (c_float*Ngas)()  # Pointer to an array of input values - u
    valueptr = (c_float*xyz)() # Pointer to an array of output mass values
    valuequantity1ptr = (c_float*xyz)() # Pointer to an array of output values - Vx
    valuequantity2ptr = (c_float*xyz)() # Pointer to an array of output values - Vy
    valuequantity3ptr = (c_float*xyz)() # Pointer to an array of output values - Vz
    valuequantity4ptr = (c_float*xyz)() # Pointer to an array of output values - u
    
    for j in range(Ngas):
        massptr[j] = snapfile.data.masses[j]
        hsmlptr[j] = 50.0#snapfile.data.hsml[j]
	quantity4ptr[j] = snapfile.data.u[j] #Internal Energy

        for m in range(3):
            x[m]=snapfile.data.pos[3*j+m]
            v[m]=snapfile.data.vel[3*j+m]

        xp=dot(R,x) # Rotate the position by the Euler angles
        vp=dot(R,v) # Rotate the velocities by the Euler angles
     	Pptr[j].x = xp[0]
     	Pptr[j].y = xp[1]
   	Pptr[j].z = xp[2]
        quantity1ptr[j] = vp[0] # Vx
        quantity2ptr[j] = vp[1] # Vy
        quantity3ptr[j] = vp[2] # Vz
    
    xmin=c_float(mass.xmin)
    xmax=c_float(mass.xmax)
    ymin=c_float(mass.ymin)
    ymax=c_float(mass.ymax)
    zmin=c_float(mass.zmin)
    zmax=c_float(mass.zmax)

    desdensngb = BulletConstants.Nsph # SPH number of neighbors

    axis1=0
    axis2=1
    axis3=2 # This is the axis which is projected out
    # I always project out the Z-axis and use Euler angles to rotate

    hmax= c_float(BulletConstants.HsmlMax) # Max Hsml size
    boxsize = c_double(snapfile.header.BoxSize[0]) # used for Periodic boundary conditions

    start = time.time()
    gridinitlib.findGridInit(Ngas, Pptr, hsmlptr, massptr, quantity1ptr, quantity2ptr, quantity3ptr, quantity4ptr, xmin, xmax, ymin, ymax, zmin, zmax, xpixels, ypixels, zpixels, desdensngb, axis1, axis2, axis3, hmax, boxsize, valueptr, valuequantity1ptr, valuequantity2ptr, valuequantity3ptr, valuequantity4ptr)

    for i in range(xpixels):  
    	for j in range(ypixels):
	    for k in range(zpixels):
                mass.data[i,j,k] = valueptr[i+j*xpixels+k*xpixels*ypixels]
		Density.data[i,j,k] = mass.data[i,j,k] / PixelVolume
                InternalEnergy.data[i,j,k] = valuequantity4ptr[i+j*xpixels+k*xpixels*ypixels]
                Vx.data[i,j,k] = valuequantity1ptr[i+j*xpixels+k*xpixels*ypixels]
                Vy.data[i,j,k] = valuequantity2ptr[i+j*xpixels+k*xpixels*ypixels]
                Vz.data[i,j,k] = valuequantity3ptr[i+j*xpixels+k*xpixels*ypixels]

    elapsed = (time.time()-start)
    print "Elapsed time to run Gadget gas gridding = "+str(elapsed)

    return [mass,Density,InternalEnergy,Vx,Vy,Vz] 
 
def GridGadgetDM(snapfile,mass,parttype,phi=0.0,theta=0.0,psi=0.0): 
    # Projects 3D DM particle data onto a 3D grid
    # Euler angles are used to rotate the data if desired

    R=EulerAngles(phi,theta,psi)
    x=zeros([3]) # Position vector
    xp=zeros([3])# Rotated position vector
    Rho=Array3d(mass.xmin,mass.xmax,mass.nx,mass.ymin,mass.ymax,mass.ny,mass.zmin,mass.zmax,mass.nz)
    PixelVolume=mass.dx * mass.dy * mass.dz
    Npart=snapfile.header.Npart
    num=Npart[parttype] # Chosen particle type

    start = time.time()

    offset=0
    for k in range(parttype):
	offset=offset+Npart[k]

    mass.data=zeros([mass.nx,mass.ny,mass.nz])+1E-12
    for n in range(num):
	no=offset+n
        for m in range(3):
            x[m]=snapfile.data.pos[3*no+m]
        xp=dot(R,x) # Rotate the position by the Euler angles

	i = int((xp[0]-mass.xmin)/mass.dx)
        j = int((xp[1]-mass.ymin)/mass.dy)
        k = int((xp[2]-mass.zmin)/mass.dz)

        if snapfile.header.Massarr[parttype]==0:
        	mp=snapfile.data.masses[no]
        else:
        	mp=snapfile.header.Massarr[parttype]

        if i>=0 and i<mass.nx and j>=0 and j<mass.ny and k>=0 and k<mass.nz:
            mass.data[i,j,k]=mass.data[i,j,k]+mp
    Rho.data=mass.data/PixelVolume

    elapsed = (time.time()-start)
    print "Elapsed time to run Gadget DM projection = "+str(elapsed)

    return [mass,Rho]


def ProjectEnzoDM(pf,data,parttype,zmin,zmax,phi=0.0,theta=0.0,psi=0.0): 
    # Projects 3D mass data onto a 2D grid
    # Euler angles are used to rotate the data if desired
    xpixels=data.nx
    ypixels=data.ny
    xy=xpixels*ypixels
    DM = Array2d(data.xmin,data.xmax,data.nx,data.ymin,data.ymax,data.ny) 

    num = int(pf.h.grid_particle_count.sum()) # of particles

    class P(Structure):
        _fields_ = [("x", c_float),("y", c_float),("z", c_float)]

    Pptr = (P*num)()  # Pointer to an array of Pos values
    massptr = (c_float*num)()  # Pointer to an array of input mass values
    valueptr = (c_float*xy)() # Pointer to an array of output mass values

    start = time.time()
    n = 0 # Particle counter	
    for i in range(pf.h.num_grids):
	PixelVolume = pf.h.grids[i]['dx'] * pf.h.grids[i]['dy'] * pf.h.grids[i]['dz']
	for j in range(int(pf.h.grid_particle_count[i])):
		massptr[n] = pf.h.grids[i]['particle_mass'][j] * PixelVolume	
     		Pptr[n].x = pf.h.grids[i]['particle_position_x'][j]
     		Pptr[n].y = pf.h.grids[i]['particle_position_y'][j]
   		Pptr[n].z = pf.h.grids[i]['particle_position_z'][j]
		if massptr[n]>10.0:
			#massptr[n] = massptr[n]*100.0
			print 'X = %.3f, Y=%.3f, Z=%.3f\n'%(Pptr[n].x, Pptr[n].y, Pptr[n].z)		
		n = n + 1

    xmin=c_float(data.xmin)
    xmax=c_float(data.xmax)
    ymin=c_float(data.ymin)
    ymax=c_float(data.ymax)
    zmin=c_float(zmin)
    zmax=c_float(zmax)
    phi=c_float(phi)
    theta=c_float(theta)
    psi=c_float(psi)
    preptime = time.time()-start

    projectdmlib.DMProject(num, Pptr, massptr, xmin, xmax, ymin, ymax, zmin, zmax, xpixels, ypixels, theta, phi, psi, valueptr)

    for i in range(xpixels):
	for j in range(ypixels):
            DM.data[i,j] = valueptr[i+j*xpixels] #valueptr[i+j*xpixels] - Why did I need to flip these? valueptr[j+i*ypixels]

    elapsed = (time.time()-start)
    print "Elapsed time to run Enzo DM projection = "+str(elapsed), "Prep time = "+str(preptime)

    return DM

def FindEnzoCentroids(pf): 

    start = time.time()
    #print 'Starting Centroid finding\n'
    mass=list()
    x=list()
    y=list()
    z=list()
    NPart = 0 # Particle counter	
    for i in range(pf.h.num_grids): # Read in all of the particle masses and positions
	PixelVolume = pf.h.grids[i]['dx'] * pf.h.grids[i]['dy'] * pf.h.grids[i]['dz']
	for j in range(int(pf.h.grid_particle_count[i])):
		#if NPart%100000 == 0:
			#print '*'
		mass.append(pf.h.grids[i]['particle_mass'][j] * PixelVolume)	
     		x.append(pf.h.grids[i]['particle_position_x'][j])
     		y.append(pf.h.grids[i]['particle_position_y'][j])
   		z.append(pf.h.grids[i]['particle_position_z'][j])
		NPart = NPart + 1

    Masses=zeros([2])
    NumPart=zeros([2])
    Centroids=zeros([2,3])
    # Finding the 2 different masses
    Masses[0] = mass[0]
    #print 'Mass of particle type 0  = %.4f'%Masses[0]
    while (Masses[1] == 0):
		pc = int(NPart * rand()) # Keep randomly choosing particles until I have found 2 different masses
		#print 'Trying another mass = %.3f'%mass[pc]
		if abs(mass[pc] - Masses[0]) > 1E-8:
			Masses[1] = mass[pc]
			#print 'Mass of particle type 1  = %.4f'%Masses[1]
    for n in range(NPart): # n cycles through the number of particles of this type
		#if n%100000 == 0:
			#print '*'
		if abs(mass[n] - Masses[0]) < 1E-8:	
			NumPart[0] = NumPart[0] + 1
			Centroids[0,0] = Centroids[0,0] + x[n]
			Centroids[0,1] = Centroids[0,1] + y[n]
			Centroids[0,2] = Centroids[0,2] + z[n]
		else:	
			NumPart[1] = NumPart[1] + 1
			Centroids[1,0] = Centroids[1,0] + x[n]
			Centroids[1,1] = Centroids[1,1] + y[n]
			Centroids[1,2] = Centroids[1,2] + z[n]

    for k in range(2): # k denotes the bullet or main particles
		if NumPart[0] > NumPart[1]: # Swap everything to put the bullet particles first
			TempNum = NumPart[0]
			NumPart[0] = NumPart[1]
			NumPart[1] = TempNum
			TempMass = Masses[0]
			Masses[0] = Masses[1]
			Masses[1] = TempMass
			for m in range(3):
				TempCentroid = Centroids[0,m]
				Centroids[0,m] = Centroids[1,m]
				Centroids[1,m] = TempCentroid

		for m in range(3):
			Centroids[k,m] = Centroids[k,m] / NumPart[k]

    elapsed = (time.time()-start)
    print "Elapsed time to locate centroids = "+str(elapsed)

    return [NumPart, Masses, Centroids]

def ProjectGadgetGas(snapfile,mass,zmin,zmax,Emin,Emax,phi=0.0,theta=0.0,psi=0.0): 
    # Projects 3D gas particle data onto a 2D grid
    # Returns both mass data and Xray intensity
    # Euler angles are used to rotate the data if desired
    xpixels=mass.nx
    ypixels=mass.ny
    xy=xpixels*ypixels
    PixelVolume=mass.dx * mass.dy * (zmax-zmin)
    R=EulerAngles(phi,theta,psi)
    x=zeros([3]) # Position vector
    xp=zeros([3])# Rotated position vector
    # Mass density will go in the input 2d array. 
    # Will create a second 2d array for the Xray data
    Xray=Array2d(mass.xmin,mass.xmax,mass.nx,mass.ymin,mass.ymax,mass.ny)

    Npart=snapfile.header.Npart
    Ngas=Npart[0] # Gas particles

    class P(Structure):
        _fields_ = [("x", c_float),("y", c_float),("z", c_float)]

    Pptr = (P*Ngas)()  # Pointer to an array of Pos values
    massptr = (c_float*Ngas)()  # Pointer to an array of input mass values
    hsmlptr = (c_float*Ngas)()  # Pointer to an array of hsml vlues
    quantityptr = (c_float*Ngas)()  # Pointer to an array of input values
    valueptr = (c_float*xy)() # Pointer to an array of output mass values
    valuequantityptr = (c_float*xy)() # Pointer to an array of output values
    desdensngb = BulletConstants.Nsph # SPH number of neighbors

    for j in range(Ngas):
        massptr[j] = snapfile.data.masses[j]
        hsmlptr[j] = snapfile.data.hsml[j]
	T = snapfile.data.u[j] * BulletConstants.Tconv
	quantityptr[j] =  BulletConstants.PreFactor * BulletConstants.TimeFactor * BulletConstants.AreaFactor * (snapfile.data.rho[j]*BulletConstants.nconv*BulletConstants.nconv) * sqrt(BulletConstants.me/T) * (exp1(Emin/T)-exp1(Emax/T))
        # Estimate of X-Ray emissions. T is in eV. Note that another factor of mass is in the subroutine.

        for m in range(3):
            x[m]=snapfile.data.pos[3*j+m]

        xp=dot(R,x) # Rotate the position by the Euler angles
     	Pptr[j].x = xp[0]
     	Pptr[j].y = xp[1]
   	Pptr[j].z = xp[2]

    xmin=c_float(mass.xmin)
    xmax=c_float(mass.xmax)
    ymin=c_float(mass.ymin)
    ymax=c_float(mass.ymax)
    zmin=c_float(zmin)
    zmax=c_float(zmax)

    axis1=0
    axis2=1
    axis3=2 # This is the axis which is projected out
    # I always project out the Z-axis and use Euler angles to rotate

    hmax = c_float(BulletConstants.HsmlMax) # Max Hsml size
    boxsize = c_double(snapfile.header.BoxSize[0]) # Used in Periodic boundary conditions

    start = time.time()

    projectlib.findHsmlAndProject(Ngas, Pptr, hsmlptr, massptr, quantityptr, xmin, xmax, ymin, ymax, zmin, zmax, xpixels, ypixels,                                 desdensngb, axis1, axis2, axis3, hmax, boxsize, valueptr, valuequantityptr)

    for i in range(xpixels):
	for j in range(ypixels):
            mass.data[i,j] = max(1E-12,valueptr[i+j*xpixels])
            Xray.data[i,j] = max(1E-12,valuequantityptr[i+j*xpixels] * valueptr[i+j*xpixels])

    elapsed = (time.time()-start)
    print "Elapsed time to run Gadget gas projection = "+str(elapsed)

    return [mass,Xray] #mass contains the mass, Xray the Xray intensity

def ProfileGadgetGas(snapfile,center,data): 
    # Projects 3D gas onto a spherical profile
    # center is the coordinates of the center
    # data is a 1d array with data extents
    rpixels=data.nx
    x=zeros([3]) # Position vector
    # Mass will go in the input 1d array. 
    # Will create additional 1d arrays for density, temperature, and pressure
    Density=Array1d(data.xmin,data.xmax,data.nx)
    Temp=Array1d(data.xmin,data.xmax,data.nx)
    Pressure=Array1d(data.xmin,data.xmax,data.nx)

    Npart=snapfile.header.Npart
    Ngas=Npart[0] # Gas particles

    class P(Structure):
        _fields_ = [("x", c_float),("y", c_float),("z", c_float)]

    Pptr = (P*Ngas)()  # Pointer to an array of Pos values
    massptr = (c_float*Ngas)()  # Pointer to an array of input mass values
    hsmlptr = (c_float*Ngas)()  # Pointer to an array of hsml vlues
    quantityptr = (c_float*Ngas)()  # Pointer to an array of input values
    valueptr = (c_float*rpixels)() # Pointer to an array of output mass values
    valuequantityptr = (c_float*rpixels)() # Pointer to an array of output values
    desdensngb = BulletConstants.Nsph # SPH number of neighbors

    for j in range(Ngas):
        massptr[j] = snapfile.data.masses[j]
        hsmlptr[j] = snapfile.data.hsml[j]
	quantityptr[j] =  snapfile.data.u[j] * BulletConstants.Tconv

        for m in range(3):
            x[m]=snapfile.data.pos[3*j+m]

     	Pptr[j].x = x[0]
     	Pptr[j].y = x[1]
   	Pptr[j].z = x[2]

    rmax=c_float(data.xmax)
    xcen=c_float(center[0])
    ycen=c_float(center[1])
    zcen=c_float(center[2])

    hmax = c_float(BulletConstants.HsmlMax) # Max Hsml size
    boxsize = c_double(snapfile.header.BoxSize[0]) # Used in Periodic boundary conditions

    start = time.time()
    profilelib.Profile(Ngas, Pptr, hsmlptr, massptr, quantityptr, xcen, ycen, zcen, rmax, rpixels, desdensngb, hmax, boxsize, valueptr, valuequantityptr)

    for i in range(rpixels):
	    if i == 0:
		ShellVolume = 4.0 * pi * data.dx**3.0 / 3.0
	    else:
		ShellVolume = 4.0 * pi * data.x[i]**2.0 * data.dx
            data.data[i] = max(1E-12,valueptr[i])
            Density.data[i] = data.data[i] / ShellVolume
	    Temp.data[i] = max(1E-12,valuequantityptr[i])
	    Pressure.data[i] = Density.data[i] * Temp.data[i] / BulletConstants.Tconv * BulletConstants.GammaMinus1

    elapsed = (time.time()-start)
    print "Elapsed time to run Gadget gas profile = "+str(elapsed)

    return [data,Density,Temp,Pressure] #data contains the mass

def ProfileGadgetGasInit(icfile,center,data): 
    # Projects 3D gas onto a spherical profile
    # center is the coordinates of the center
    # data is a 1d array with data extents
    rpixels=data.nx
    x=zeros([3]) # Position vector
    # Mass will go in the input 1d array. 
    # Will create additional 1d arrays for density, temperature, and pressure
    Density=Array1d(data.xmin,data.xmax,data.nx)
    Temp=Array1d(data.xmin,data.xmax,data.nx)
    Pressure=Array1d(data.xmin,data.xmax,data.nx)

    Npart=icfile.header.Npart
    Ngas=Npart[0] # Gas particles

    class P(Structure):
        _fields_ = [("x", c_float),("y", c_float),("z", c_float)]

    Pptr = (P*Ngas)()  # Pointer to an array of Pos values
    massptr = (c_float*Ngas)()  # Pointer to an array of input mass values
    hsmlptr = (c_float*Ngas)()  # Pointer to an array of hsml vlues
    quantityptr = (c_float*Ngas)()  # Pointer to an array of input values
    valueptr = (c_float*rpixels)() # Pointer to an array of output mass values
    valuequantityptr = (c_float*rpixels)() # Pointer to an array of output values
    desdensngb = BulletConstants.Nsph # SPH number of neighbors

    for j in range(Ngas):
        massptr[j] = icfile.data.masses[j]
        hsmlptr[j] = 10.0 # Temp hack icfile.data.hsml[j]
	quantityptr[j] =  icfile.data.u[j] * BulletConstants.Tconv

        for m in range(3):
            x[m]=icfile.data.pos[3*j+m]

     	Pptr[j].x = x[0]
     	Pptr[j].y = x[1]
   	Pptr[j].z = x[2]

    rmax=c_float(data.xmax)
    xcen=c_float(center[0])
    ycen=c_float(center[1])
    zcen=c_float(center[2])

    hmax = c_float(BulletConstants.HsmlMax) # Max Hsml size
    boxsize = c_double(icfile.header.BoxSize[0]) # Used in Periodic boundary conditions

    start = time.time()
    profilelib.Profile(Ngas, Pptr, hsmlptr, massptr, quantityptr, xcen, ycen, zcen, rmax, rpixels, desdensngb, hmax, boxsize, valueptr, valuequantityptr)

    for i in range(rpixels):
	    if i == 0:
		ShellVolume = 4.0 * pi * data.dx**3.0 / 3.0
	    else:
		ShellVolume = 4.0 * pi * data.x[i]**2.0 * data.dx
            data.data[i] = max(1E-12,valueptr[i])
	    #print "Radius of shell %d = %f, Mass in shell %d = %f\n"%(i,data.x[i],i,valueptr[i])
            Density.data[i] = data.data[i] / ShellVolume
	    Temp.data[i] = max(1E-12,valuequantityptr[i])
	    Pressure.data[i] = Density.data[i] * Temp.data[i] / BulletConstants.Tconv * BulletConstants.GammaMinus1

    elapsed = (time.time()-start)
    print "Elapsed time to run Gadget gas profile = "+str(elapsed)

    return [data,Density,Temp,Pressure] #data contains the mass


def ProjectGadgetGasSlow(snapfile,mass,zmin,zmax,nz,phi=0.0,theta=0.0,psi=0.0): 
    # Projects 3D gas particle data onto a 2D grid, which is input. nz is the number of zsteps.
    # Returns both mass data and Xray intensity
    # Euler angles are used to rotate the data if desired
    # This works by gridding the gas, and then adding up the grid elements, and is much slower than projectgas
    xpixels=mass.nx
    ypixels=mass.ny

    mass3d=Array3d(mass.xmin,mass.xmax,mass.nx,mass.ymin,mass.ymax,mass.ny,zmin,zmax,nz)

    [mass3d,temp3d,rho3d,xray3d,press3d,SZ3d]=GridGadgetGas(snapfile,mass3d,phi,theta,psi)

    Xray=Array2d(mass.xmin,mass.xmax,mass.nx,mass.ymin,mass.ymax,mass.ny)
    SZ=Array2d(mass.xmin,mass.xmax,mass.nx,mass.ymin,mass.ymax,mass.ny)
    for i in range(xpixels):
	for j in range(ypixels):
		mass.data[i,j]=0.0
		Xray.data[i,j]=0.0
		for k in range(nz):
            		mass.data[i,j] = mass.data[i,j] + mass3d.data[i,j,k] # Sum up along a column
      	      		Xray.data[i,j] = Xray.data[i,j] + xray3d.data[i,j,k] # Sum up along a column
			SZ.data[i,j]   = SZ.data[i,j]   + SZ3d.data[i,j,k]   # Sum up along a column
    return [mass,Xray,SZ] # mass contains the mass, Xray the Xray intensity, SZ the SZ effect data

def ProjectGadgetGasSZ(snapfile,mass,zmin,zmax,phi=0.0,theta=0.0,psi=0.0): 
    # Projects 3D gas particle data onto a 2D grid
    # Returns both mass data and Xray intensity
    # Euler angles are used to rotate the data if desired
    xpixels=mass.nx
    ypixels=mass.ny
    xy=xpixels*ypixels
    PixelVolume=mass.dx * mass.dy * (zmax-zmin)
    PixelArea = mass.dx * mass.dy
    R=EulerAngles(phi,theta,psi)
    x=zeros([3]) # Position vector
    xp=zeros([3])# Rotated position vector
    # Mass density will go in the input 2d array. 
    # Will create additional 2d arrays for the SZ data
    SZ=Array2d(mass.xmin,mass.xmax,mass.nx,mass.ymin,mass.ymax,mass.ny)

    Npart=snapfile.header.Npart
    Ngas=Npart[0] # Gas particles

    class P(Structure):
        _fields_ = [("x", c_float),("y", c_float),("z", c_float)]

    Pptr = (P*Ngas)()  # Pointer to an array of Pos values
    massptr = (c_float*Ngas)()  # Pointer to an array of input mass values
    hsmlptr = (c_float*Ngas)()  # Pointer to an array of hsml vlues
    quantityptr = (c_float*Ngas)()  # Pointer to an array of input values
    valueptr = (c_float*xy)() # Pointer to an array of output mass values
    valuequantityptr = (c_float*xy)() # Pointer to an array of output values
    desdensngb = 50 # SPH number of neighbors

    for j in range(Ngas):
        massptr[j] = snapfile.data.masses[j]
        hsmlptr[j] = snapfile.data.hsml[j]
	T = snapfile.data.u[j] * BulletConstants.Tconv
	quantityptr[j] =  BulletConstants.SZFactor *  T / PixelArea
        # I think this is the correct way to calculate the SZ Delta T.

        for m in range(3):
            x[m]=snapfile.data.pos[3*j+m]

        xp=dot(R,x) # Rotate the position by the Euler angles
     	Pptr[j].x = xp[0]
     	Pptr[j].y = xp[1]
   	Pptr[j].z = xp[2]

    xmin=c_float(mass.xmin)
    xmax=c_float(mass.xmax)
    ymin=c_float(mass.ymin)
    ymax=c_float(mass.ymax)
    zmin=c_float(zmin)
    zmax=c_float(zmax)

    axis1=0
    axis2=1
    axis3=2 # This is the axis which is projected out
    # I always project out the Z-axis and use Euler angles to rotate

    hmax = c_float(1000.0) # Max Hsml size
    boxsize = c_double(snapfile.header.BoxSize[0]) # Used in Periodic boundary conditions

    start = time.time()

    projectlib.findHsmlAndProject(Ngas, Pptr, hsmlptr, massptr, quantityptr, xmin, xmax, ymin, ymax, zmin, zmax, xpixels, ypixels,                                 desdensngb, axis1, axis2, axis3, hmax, boxsize, valueptr, valuequantityptr)

    for i in range(xpixels):
	for j in range(ypixels):
            mass.data[i,j] = max(1E-12,valueptr[i+j*xpixels])
            SZ.data[i,j] = valuequantityptr[i+j*xpixels] * valueptr[i+j*xpixels]

    elapsed = (time.time()-start)
    print "Elapsed time to run Gadget gas projection = "+str(elapsed)

    return [mass,SZ] #mass contains the mass, SZ the SZ deviation

def ProjectGadgetGasTemp(snapfile,mass,zmin,zmax,phi=0.0,theta=0.0,psi=0.0): 
    # Projects 3D gas particle data onto a 2D grid
    # Returns both mass data andTemp
    # Euler angles are used to rotate the data if desired
    xpixels=mass.nx
    ypixels=mass.ny
    xy=xpixels*ypixels
    PixelVolume=mass.dx * mass.dy * (zmax-zmin)
    PixelArea = mass.dx * mass.dy
    R=EulerAngles(phi,theta,psi)
    x=zeros([3]) # Position vector
    xp=zeros([3])# Rotated position vector
    # Mass density will go in the input 2d array. 
    # Will create additional 2d arrays for the SZ data
    Temp=Array2d(mass.xmin,mass.xmax,mass.nx,mass.ymin,mass.ymax,mass.ny)

    Npart=snapfile.header.Npart
    Ngas=Npart[0] # Gas particles

    class P(Structure):
        _fields_ = [("x", c_float),("y", c_float),("z", c_float)]

    Pptr = (P*Ngas)()  # Pointer to an array of Pos values
    massptr = (c_float*Ngas)()  # Pointer to an array of input mass values
    hsmlptr = (c_float*Ngas)()  # Pointer to an array of hsml vlues
    quantityptr = (c_float*Ngas)()  # Pointer to an array of input values
    valueptr = (c_float*xy)() # Pointer to an array of output mass values
    valuequantityptr = (c_float*xy)() # Pointer to an array of output values
    desdensngb = 50 # SPH number of neighbors

    for j in range(Ngas):
	T = snapfile.data.u[j] * BulletConstants.Tconv 
        massptr[j] = (snapfile.data.rho[j] * snapfile.data.rho[j]) * (exp1(BulletConstants.Emin/T)-exp1(BulletConstants.Emax/T)) / sqrt(T)
        hsmlptr[j] = snapfile.data.hsml[j]
	quantityptr[j] =  T / 1000.0

        for m in range(3):
            x[m]=snapfile.data.pos[3*j+m]

        xp=dot(R,x) # Rotate the position by the Euler angles
     	Pptr[j].x = xp[0]
     	Pptr[j].y = xp[1]
   	Pptr[j].z = xp[2]

    xmin=c_float(mass.xmin)
    xmax=c_float(mass.xmax)
    ymin=c_float(mass.ymin)
    ymax=c_float(mass.ymax)
    zmin=c_float(zmin)
    zmax=c_float(zmax)

    axis1=0
    axis2=1
    axis3=2 # This is the axis which is projected out
    # I always project out the Z-axis and use Euler angles to rotate

    hmax = c_float(1000.0) # Max Hsml size
    boxsize = c_double(snapfile.header.BoxSize[0]) # Used in Periodic boundary conditions

    start = time.time()

    projectlib.findHsmlAndProject(Ngas, Pptr, hsmlptr, massptr, quantityptr, xmin, xmax, ymin, ymax, zmin, zmax, xpixels, ypixels,                                 desdensngb, axis1, axis2, axis3, hmax, boxsize, valueptr, valuequantityptr)

    for i in range(xpixels):
	for j in range(ypixels):
            mass.data[i,j] = max(1E-12,valueptr[i+j*xpixels])
            Temp.data[i,j] = max(1E-12, min(30.0, valuequantityptr[i+j*xpixels])) / (1 + BulletConstants.redshift)

    elapsed = (time.time()-start)
    print "Elapsed time to run Gadget gas projection = "+str(elapsed)

    return Temp


def ProjectGadgetGasSZTemp(snapfile,mass,zmin,zmax,phi=0.0,theta=0.0,psi=0.0): 
    # Projects 3D gas particle data onto a 2D grid
    # Returns both mass data, Xray intensity, SZ data, and Temp data
    # Euler angles are used to rotate the data if desired
    xpixels=mass.nx
    ypixels=mass.ny
    xy=xpixels*ypixels
    PixelVolume=mass.dx * mass.dy * (zmax-zmin)
    PixelArea = mass.dx * mass.dy
    R=EulerAngles(phi,theta,psi)
    x=zeros([3]) # Position vector
    xp=zeros([3])# Rotated position vector
    # Mass density will go in the input 2d array. 
    # Will create additional 2d arrays for the Xray, SZ, and Temp data
    Xray=Array2d(mass.xmin,mass.xmax,mass.nx,mass.ymin,mass.ymax,mass.ny)    
    SZ=Array2d(mass.xmin,mass.xmax,mass.nx,mass.ymin,mass.ymax,mass.ny)    
    Temp=Array2d(mass.xmin,mass.xmax,mass.nx,mass.ymin,mass.ymax,mass.ny)

    Npart=snapfile.header.Npart
    Ngas=Npart[0] # Gas particles

    Nvals = 3 # This routine projects out 3 values, XRay intensity, SZ effect, and Temp

    QNum = Ngas*Nvals
    VQNum = xy*Nvals

    class P(Structure):
        _fields_ = [("x", c_float),("y", c_float),("z", c_float)]

    Pptr = (P*Ngas)()  # Pointer to an array of Pos values
    massptr = (c_float*Ngas)()  # Pointer to an array of input mass values
    hsmlptr = (c_float*Ngas)()  # Pointer to an array of hsml vlues
    quantityptr = (c_float*QNum)()  # Pointer to an array of input values
    valueptr = (c_float*xy)() # Pointer to an array of output mass values
    valuequantityptr = (c_float*VQNum)() # Pointer to an array of output values
    desdensngb = 50 # SPH number of neighbors

    XrayPreFactor = BulletConstants.PreFactor * BulletConstants.TimeFactor * BulletConstants.AreaFactor * BulletConstants.nconv*BulletConstants.nconv
    SZPreFactor = BulletConstants.SZFactor / PixelArea
    for j in range(Ngas):
        massptr[j] = snapfile.data.masses[j]
        hsmlptr[j] = snapfile.data.hsml[j]
	T = snapfile.data.u[j] * BulletConstants.Tconv
	quantityptr[j] =  XrayPreFactor * snapfile.data.rho[j] * sqrt(BulletConstants.me/T) * (exp1(BulletConstants.Emin/T)-exp1(BulletConstants.Emax/T))
	quantityptr[j + Ngas] =  SZPreFactor *  T 
	quantityptr[j + 2 * Ngas] =  T  / 1000.0 # Will return T in keV to compare to Markevitch maps

        for m in range(3):
            x[m]=snapfile.data.pos[3*j+m]

        xp=dot(R,x) # Rotate the position by the Euler angles
     	Pptr[j].x = x[0]
     	Pptr[j].y = x[1]
   	Pptr[j].z = x[2]

    xmin=c_float(mass.xmin)
    xmax=c_float(mass.xmax)
    ymin=c_float(mass.ymin)
    ymax=c_float(mass.ymax)
    zmin=c_float(zmin)
    zmax=c_float(zmax)

    axis1=0
    axis2=1
    axis3=2 # This is the axis which is projected out
    # I always project out the Z-axis and use Euler angles to rotate

    hmax = c_float(1000.0) # Max Hsml size
    boxsize = c_double(snapfile.header.BoxSize[0]) # Used in Periodic boundary conditions

    start = time.time()

    projectnlib.findHsmlAndProjectN(Nvals, Ngas, Pptr, hsmlptr, massptr, quantityptr, xmin, xmax, ymin, ymax, zmin, zmax, xpixels, ypixels, desdensngb, axis1, axis2, axis3, hmax, boxsize, valueptr, valuequantityptr)

    for i in range(xpixels):
	for j in range(ypixels):
            mass.data[i,j] = max(1E-12,valueptr[i+j*xpixels])
            Xray.data[i,j] = max(1E-12,valuequantityptr[i+j*xpixels] * valueptr[i+j*xpixels])
            SZ.data[i,j] = valuequantityptr[i+j*xpixels+xy] * valueptr[i+j*xpixels]
            Temp.data[i,j] = max(1E-12,valuequantityptr[i+j*xpixels+2*xy]) 

    elapsed = (time.time()-start)
    print "Elapsed time to run Gadget gas projection = "+str(elapsed)

    return [mass,Xray,SZ,Temp] #mass contains the mass, SZ the SZ deviation




def ProjectGadgetDM(snapfile,data,parttype,zmin,zmax,phi=0.0,theta=0.0,psi=0.0): 
    # Projects 3D mass data onto a 2D grid
    # Euler angles are used to rotate the data if desired
    xpixels=data.nx
    ypixels=data.ny
    xy=xpixels*ypixels
    data.data=zeros([data.nx,data.ny])# Mass 

    Npart=snapfile.header.Npart
    num=Npart[parttype] # Chosen particle type

    class P(Structure):
        _fields_ = [("x", c_float),("y", c_float),("z", c_float)]

    Pptr = (P*num)()  # Pointer to an array of Pos values
    massptr = (c_float*num)()  # Pointer to an array of input mass values
    valueptr = (c_float*xy)() # Pointer to an array of output mass values

    start = time.time()

    offset=0
    for k in range(parttype):
	offset=offset+Npart[k]

    for n in range(num):
	no=offset+n
     	Pptr[n].x = snapfile.data.pos[3*no]
     	Pptr[n].y = snapfile.data.pos[3*no+1]
   	Pptr[n].z = snapfile.data.pos[3*no+2]
        if snapfile.header.Massarr[parttype]==0:
        	massptr[n]=snapfile.data.masses[no]
        else:
        	massptr[n]=snapfile.header.Massarr[parttype]

    xmin=c_float(data.xmin)
    xmax=c_float(data.xmax)
    ymin=c_float(data.ymin)
    ymax=c_float(data.ymax)
    zmin=c_float(zmin)
    zmax=c_float(zmax)
    phi=c_float(phi)
    theta=c_float(theta)
    psi=c_float(psi)

    projectdmlib.DMProject(num, Pptr, massptr, xmin, xmax, ymin, ymax, zmin, zmax, xpixels, ypixels, theta, phi, psi, valueptr)

    for i in range(xpixels):
	for j in range(ypixels):
            data.data[i,j] = valueptr[i+j*xpixels]

    elapsed = (time.time()-start)
    print "Elapsed time to run Gadget DM projection = "+str(elapsed)

    return data

def ProjectGadgetDMandGalaxies(snapfile,data,zmin,zmax,phi=0.0,theta=0.0,psi=0.0): 
    # Projects 3D mass data onto a 2D grid
    # Euler angles are used to rotate the data if desired
    xpixels=data.nx
    ypixels=data.ny
    xy=xpixels*ypixels
    data.data=zeros([data.nx,data.ny])# Mass 

    Npart=snapfile.header.Npart
    num=Npart[1] + Npart[4]

    class P(Structure):
        _fields_ = [("x", c_float),("y", c_float),("z", c_float)]

    Pptr = (P*num)()  # Pointer to an array of Pos values
    massptr = (c_float*num)()  # Pointer to an array of input mass values
    valueptr = (c_float*xy)() # Pointer to an array of output mass values

    start = time.time()

    offset=Npart[0]

    for n in range(num):
	no=offset+n
     	Pptr[n].x = snapfile.data.pos[3*no]
     	Pptr[n].y = snapfile.data.pos[3*no+1]
   	Pptr[n].z = snapfile.data.pos[3*no+2]
        if snapfile.header.Massarr[1]==0:
        	massptr[n]=snapfile.data.masses[no]
        else:
        	massptr[n]=snapfile.header.Massarr[1]

    xmin=c_float(data.xmin)
    xmax=c_float(data.xmax)
    ymin=c_float(data.ymin)
    ymax=c_float(data.ymax)
    zmin=c_float(zmin)
    zmax=c_float(zmax)
    phi=c_float(phi)
    theta=c_float(theta)
    psi=c_float(psi)

    projectdmlib.DMProject(num, Pptr, massptr, xmin, xmax, ymin, ymax, zmin, zmax, xpixels, ypixels, theta, phi, psi, valueptr)

    for i in range(xpixels):
	for j in range(ypixels):
            data.data[i,j] = valueptr[i+j*xpixels]

    elapsed = (time.time()-start)
    print "Elapsed time to run Gadget DM projection = "+str(elapsed)

    return data


def GetBulletData(filename,data):
	infile=open(filename,'r')
	oneddata=infile.readlines()
	for i in range(data.nx):
		for j in range(data.ny):
			data.data[i,j]=float(oneddata[i+j*data.nx])
	return data

def GetBulletFits(filename,data):
	fits=pyfits.open(filename)
	fitsdata=fits[0].data
	for i in range(data.nx):
		for j in range(data.ny):
			data.data[i,j]=fitsdata[i,j]
	return data

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

def PyramidalKernel3d(deltax,deltay,deltaz):
	if deltax>=1.0 or deltay>=1.0 or deltaz>=1.0:
		return 0.0
	else:
		return (1.0-deltax)*(1.0-deltay)*(1.0-deltaz)

def DataInterpolate3d(data,xprime,yprime,zprime):
	i=int((xprime-data.xmin)/data.dx)
	j=int((yprime-data.ymin)/data.dy)
	k=int((zprime-data.zmin)/data.dz)
	d=0.0
	for m in [i-1,i,i+1]:
	    ml = max(0,min(data.nx-1,m))
            deltax=abs((xprime-data.x[ml])/data.dx)
            for n in [j-1,j,j+1]:
	        nl = max(0,min(data.ny-1,n))
                deltay=abs((yprime-data.y[nl])/data.dy)
		for p in [k-1,k,k+1]:
		        pl = max(0,min(data.nz-1,p))
	                deltaz=abs((zprime-data.z[pl])/data.dz)
	                d=d+PyramidalKernel3d(deltax,deltay,deltaz)*data.data[ml,nl,pl]
	return d

def OneDSlice2(dataA,dataB,sliceA,sliceB,theta):
	imax=int(dataA.data.argmax()/dataA.ny)
	jmax=dataA.data.argmax()-imax*dataA.ny
	for i in range(sliceA.nx):
            xprime=dataA.x[i+1]
            sliceA.x[i]=xprime/cos(theta)
            sliceB.x[i]=sliceA.x[i]
            yprime=dataA.y[jmax]+(i+1-imax)*dataA.dx*sin(theta)
            sliceA.y[i]=yprime
            sliceA.data[i]=DataInterpolate(dataA,xprime,yprime)
            sliceB.data[i]=DataInterpolate(dataB,xprime,yprime)
        return [sliceA,sliceB]

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


def OneDSlice3(dataA,dataB,dataC,sliceA,sliceB,sliceC,theta):
	imax=int(dataA.data.argmax()/dataA.ny)
	jmax=dataA.data.argmax()-imax*dataA.ny
	for i in range(sliceA.nx):
            xprime=dataA.x[i+1]
            sliceA.x[i]=xprime/cos(theta)
            sliceB.x[i]=sliceA.x[i]
            sliceC.x[i]=sliceA.x[i]
            yprime=dataA.y[jmax]+(i+1-imax)*dataA.dx*sin(theta)
            sliceA.data[i]=DataInterpolate(dataA,xprime,yprime)
            sliceB.data[i]=DataInterpolate(dataB,xprime,yprime)
            sliceC.data[i]=DataInterpolate(dataC,xprime,yprime)
        return [sliceA,sliceB,sliceC]

def FindBestShift(data1list,data2list,sigmalist,masklist,align,tol):
    # This routine finds the best alignment between n sets of 2D
    # arrays.  data1 is a list of the larger arrays, and data2 is a list of the smaller,
    # unshifted arrays. sigma is a list of the data standard deviation

    start = time.time()

    numarrays = len(data1list)
    calign=CAlign(align)
    shifteddata1list=list()

    for i in range(numarrays):
	shifteddata1list.append(Array2d(data2list[i].xmin,data2list[i].xmax,data2list[i].nx,data2list[i].ymin,data2list[i].ymax,data2list[i].ny))
    
    arrayset=ArraySet(numarrays,data1list,data2list,shifteddata1list,sigmalist,masklist)

    diff=c_double(0.0)
    gtol=c_double(tol)

    shiftnlib.optimize(byref(arrayset),byref(calign),byref(diff),gtol)
    # c++ subroutine finds best shift

    for i in range(numarrays):
    	ArrayCopyCToPython(arrayset.shifteddata1[i],shifteddata1list[i]) # copy the results back to Python formats

    AlignCopyCToPython(calign,align)
    fom=diff.value

    elapsed = (time.time()-start)
    print "Elapsed time to find optimal alignment = "+str(elapsed)+"\n"

    return [shifteddata1list,align,fom]

def ArrayCopyCToPython(cdata,data):
    # Copies from c to Python
    data.nx=cdata.nx
    data.ny=cdata.ny
    data.xmin=cdata.xmin
    data.xmax=cdata.xmax
    data.ymin=cdata.ymin
    data.ymax=cdata.ymax
    data.dx=cdata.dx
    data.dy=cdata.dy
    
    for j in range(data.ny):
        data.y[j]=cdata.y[j]

    for i in range(data.nx):
        data.x[i]=cdata.x[i]
	for j in range(data.ny):
            data.data[i,j]=cdata.data[i+j*data.nx]

    return data

def AlignCopyCToPython(calign,align):
    # Copies from c to Python
    for i in range(3):
        align.d[i]=calign.d[i]
        align.dmin[i]=calign.dmin[i]
        align.dmax[i]=calign.dmax[i]
    
    return align

def LinearSearch(xold, fold, g, p, stpmax, func):
#This does a linear minimum search. It starts with the full Newton step, 
#then backtracks to find a minimum.  This is adapted from Numerical Recipes.

    #print "Into LinearSearch\n"
    ALF=1.0E-6
    TOLX=1.0E-9
    alam2=0.0
    f2=0.0
    slope=0.0
    summ=0.0
    n=len(xold)
    x=zeros([n])
    xu=zeros([n])
    for i in range(n): summ = summ + p[i]*p[i]
    summ = sqrt(summ)
    if summ>stpmax: 
	for i in range(n): p[i] = p[i] * stpmax/summ
    for i in range(n): slope = slope + g[i]*p[i]
    #print "g=",g," p=",p," slope=",slope,"\n"
    if slope>0.0: print "Roundoff problem in LinearSearch."
    test=0.0;
    for i in range(n):
        temp=abs(p[i])/max(abs(xold[i]),1.0)
        if temp>test: test=temp
    #print "temp=",temp, " test=", test,"\n"
    alamin=TOLX/test
    alam=1.0
    while True:
        #print "alamin=%.6f"%alamin, "alam=%.6f\n"%alam
        for i in range(n): 
            x[i]=xold[i]+alam*p[i]
 	f=func(x)
        if alam<alamin:
            for i in range(n): x[i]=xold[i]
            #print "LinearSearch hit alamin\n"
            return [x,f]
        elif f<=(fold+ALF*alam*slope): 
            #print "Normal return from LinearSearch\n"
            return [x,f]
        else:
            if alam==1.0:tmplam=-slope/(2.0*(f-fold-slope))
            else:
                rhs1=f-fold-alam*slope
                rhs2=f2-fold-alam2*slope
                a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2)
                b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2)
                if a==0.0: tmplam=-slope/(2.0*b)
                else:
                    disc=b*b-3.0*a*slope
                    if disc<0.0: tmplam=0.5*alam
                    elif b<=0.0: tmplam=(-b+sqrt(disc))/(3.0*a)
                    else: tmplam=-slope/(b+sqrt(disc))
                if tmplam>0.5*alam: tmplam=0.5*alam
        alam2=alam
        f2=f
        alam=max(tmplam,0.1*alam)

def FindMinimum(p,  gtol, func, stpmax=.05, eps=.01):
#This finds the minimum of a multidimensional function.
#This is adapted from Numerical Recipes.
  ITMAX=200
  EPS=eps
  TOLX=1.0E-9
  n=len(p)
  hessin=zeros([n,n])
  xi=zeros([n])
  dg=zeros([n])
  hdg=zeros([n])
  fp=func(p)
  g=FunctionGradient(p,fp,EPS,func)
  #print "Finished grad g= ",g,"\n"
  for i in range(n):
      xi[i]=-g[i]
      for j in range(n): hessin[i,j]=0.0
      hessin[i,i]=1.0
  #print "Hessian=" ,hessin, "\n"
  for its in range(ITMAX):
      #print "New iteration xi= ",xi,"\n"
      iter=its
      #print "Into its loop\n"
      [pnew,fret]=LinearSearch(p,fp,g,xi,stpmax,func)
      #print "pnew=", pnew
      fp=fret
      for i in range(n):
          xi[i]=pnew[i]-p[i]
          p[i]=pnew[i]
      test=0.0
      for i in range(n):
          temp=abs(xi[i])/max(abs(p[i]),1.0)
          if temp>test: test=temp
      if test<TOLX: 
          return [p, fret, iter]
      for i in range(n): dg[i]=g[i]
      g=FunctionGradient(p,fret,EPS,func)
      test=0.0
      den=max(abs(fret),1.0)
      for i in range(n):
          temp=abs(g[i])*max(abs(p[i]),1.0)/den
          if temp>test: test=temp
      if test<gtol: return [p, fret, iter]
      for i in range(n):dg[i]=g[i]-dg[i]
      for i in range(n):    
          hdg[i]=0.0
          for j in range(n): hdg[i]=hdg[i]+hessin[i,j]*dg[j]
      fac=fae=sumdg=sumxi=0.0
      for i in range(n):
          fac = fac + dg[i]*xi[i]
          fae = fae + dg[i]*hdg[i]
          sumdg = sumdg + dg[i]*dg[i]
          sumxi = sumxi + xi[i]*xi[i]
      if fac>sqrt(EPS*sumdg*sumxi):
          fac=1.0/fac
          fad=1.0/fae
          for i in range(n): dg[i]=fac*xi[i]-fad*hdg[i]
          for i in range(n):
              for j in range(i,n):
                  hessin[i,j] = hessin[i,j]+fac*xi[i]*xi[j] - fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j]
                  hessin[j,i]=hessin[i,j]
	  #print "New Hessian=" ,hessin, "\n"
	  DetH=linalg.det(hessin)
	  #print "Det H = ",DetH,"\n"
	  if DetH<0.01 :
		  for i in range(n):
		      for j in range(n): hessin[i,j]=0.0
		      hessin[i,i]=1.0
		  #print "DetH too small - resetting Hessian to identity \n"
      for i in range(n):
          xi[i]=0.0
          for j in range(n): 
              xi[i] = xi[i]-hessin[i,j]*g[j]
  print "too many iterations in FindMinimum\n"
  return [p, fret, iter]

def FunctionGradient(x,fold,EPS,func):
  n=len(x)
  df=zeros([n])
  xh=zeros([n])
  for i in range(n): xh[i]=x[i]
  for j in range(n):
      temp=x[j]
      h=EPS*abs(temp)
      if h==0: h=EPS
      xh[j]=temp+h
      fh=func(xh)
      xh[j]=temp
      df[j]=(fh-fold)/h
  return df

					
def FindFom(snapmin,snapmax,Z,phi,theta,psi,RelVel,ConstrainPhi=True,TFudge = 1.0,Mask=(1,0,0,0,0,0,0)):

	# This is a new, faster Fom finding routine that uses only mass and X-ray, and searches for a 
	# minimum in time, psi, and theta.
	# This version first scans through a Theta,Psi matrix to deterine the optimum.

        (dataA,sigmaA,maskA,maskANull,dataB1,sigmaB1,dataB2,sigmaB2,dataB3,sigmaB3,dataC,sigmaC,dataD,sigmaD,maskD,maskDNull) = GetData()

	[dyyA,dxxA] = meshgrid(dataA.y,dataA.x)# Data grid for plots
	[dyyB,dxxB] = meshgrid(dataB1.y,dataB1.x)# Data grid for plots
	[dyyD,dxxD] = meshgrid(dataD.y,dataD.x)# Data grid for plots

        MinVz = BulletConstants.BulletVz - 6.0 * BulletConstants.BulletSigVz # Minimum Bullet radial velocity consistent with observations
        MaxVz = BulletConstants.BulletVz + 6.0 * BulletConstants.BulletSigVz # Maximum Bullet radial velocity consistent with observations
	DeltaPsi = 5 # Using Psi and Theta 100* psi and theta, so I can use integers
	DeltaTheta = 5 # Using DeltaPsi and Delta theta of 0.05 radians or about 3 degrees
	PsiMax = int(100*psi) + 30
  	PsiMin = int(100*psi) - 30
	ThetaMax = int(100*theta) + 30
	ThetaMin = int(100*theta) - 30
	NPsi = int((PsiMax-PsiMin)/DeltaPsi) + 3
	NTheta = int((ThetaMax-ThetaMin)/DeltaTheta) + 3
	FomMatrix = {}
	for snap in range(snapmin-1,snapmax+2): # First fill the whole array with 1E6 - this IDs the boundaries
		for Psi in range(PsiMin-DeltaPsi, PsiMax+2*DeltaPsi, DeltaPsi):	
			for Theta in range(ThetaMin-DeltaTheta, ThetaMax+2*DeltaTheta, DeltaTheta):
				FomMatrix[snap,Psi,Theta] = (1.0E6,0.0)				
	for snap in range(snapmin,snapmax+1): # Next put -1.0 everywhere except at the boundaries - this IDs conditions not yet run
		for Psi in range(PsiMin, PsiMax+DeltaPsi, DeltaPsi):	
			for Theta in range(ThetaMin, ThetaMax+DeltaTheta, DeltaTheta):
				FomMatrix[snap,Psi,Theta] = (-1.0,0.0)		
	snap = snapmin
        lastsimtime = 1000.0 # Nonsense value for first time through
        LastBulletDMPos = [0.0,0.0,0.0]
        LastMainDMPos = [0.0,0.0,0.0]
        if Mask[1] == 0:
            cenfilename='centroids.out'
        if Mask[1] == 1:
            cenfilename='centroidsx1.out'
        Phi = 0.0
	Theta = int(100*theta)
	Psi = int(100*psi)
	counter = 0
	MaxCounter = (snapmax-snapmin+1) * NPsi * NTheta
	TimeStripeDone = False
	ThetaStripeDone = False
	bestfom = 1.0E6
	bestsnap = snapmin
	bestTheta = ThetaMin
	bestPsi = PsiMin
	while counter < MaxCounter:
		psi = Psi / 100.0
		theta = Theta / 100.0
		counter = counter + 1
		print "%d times through - Snap = %d, Psi = %.2f, Theta = %.2f"%(counter, snap, psi, theta)
                sys.stdout.flush()
                R = EulerAngles(-psi,-theta,0.0)
		Vz = dot(R,RelVel)[2] # This is the observed Z_Velocity of the bullet relative to the CM
		if Vz > MaxVz or Vz < MinVz: # If Vz is outside of observed Mean +/- 4 Sigma
			FomMatrix[snap,Psi,Theta] = (1.0E5,0.0) # If it's outside allowed Vz, give it a large FOM
                        fom = 1.0E5
                        Phi = 0.0
                        print 'Outside allowed Vz, snap = %d, Vz = %.2f, MinVz = %.2f, MaxVz = %.2f\n'%(snap,Vz,MinVz,MaxVz)
                        sys.stdout.flush()
		else:
                    try:
			pf = GetPF(snap)
                        simtime = pf.h.parameters['InitialTime'] * BulletConstants.TimeConversion # Put time in Gyears
                        dmsim=Array2d(2.0*dataA.xmin,2.0*dataA.xmax,2*dataA.nx,2.0*dataA.ymin,2.0*dataA.ymax,2*dataA.ny)
                        masssim=Array2d(2.0*dataA.xmin,2.0*dataA.xmax,2*dataA.nx,2.0*dataA.ymin,2.0*dataA.ymax,2*dataA.ny)
                        sumsim=Array2d(2.0*dataA.xmin,2.0*dataA.xmax,2*dataA.nx,2.0*dataA.ymin,2.0*dataA.ymax,2*dataA.ny)
                        [syyA,sxxA] = meshgrid(dmsim.y,dmsim.x) # Sim grid for plots

			ApecData = ReadLookups(Z) # Reads the APEC lookup tables.
			for grid in pf.h.grids:
				grid.set_field_parameter('ApecData',ApecData)
                                grid.set_field_parameter('TFudge',TFudge)
                        DMProject = False
			[dmsim,masssim,xraysim1,xraysim2,xraysim3,szsim]=ProjectEnzoData(pf,masssim,phi,theta,psi,DMProject=DMProject)
			sumsim.data=gaussian_filter(dmsim.data+masssim.data,2.0)
	
                        align = SetAlign(dataB1, phi, theta, psi, ConstrainPhi = ConstrainPhi)

			# The 0.22 is the approximate alignment of the bullet cluster on the sky. 
			tol=1.0E-7 # Alignment tolerance
			data1list=list((sumsim,xraysim1))
			data2list=list((dataA,dataB1))
			sigmalist=list((sigmaA,sigmaB1))
                        if Mask[1] == 1:
                            masklist=list((maskA,maskA))
                        else:
                            masklist=list((maskA,maskANull))

			[shifteddata1list,align,fom]=FindBestShift(data1list, data2list, sigmalist, masklist, align, tol)
			FomMatrix[snap,Psi,Theta] = (fom, align.d[2]) # Store (fom, phi)
                    except:
			FomMatrix[snap,Psi,Theta] = (1.0E5, 0.0) # If there's an error, give it a large FOM
                        fom = 1.0E5                              
                        Phi = 0.0
		if fom < bestfom:
			bestfom = fom
			bestsnap = snap
			bestTheta = Theta
			bestPsi = Psi

		if snap == snapmax and not TimeStripeDone:
			TimeStripeDone = True
			snap = bestsnap
			Theta = ThetaMin
			print "Time Stripe Finished \n"
                        sys.stdout.flush()
		if not TimeStripeDone and not ThetaStripeDone: # First, complete one stripe through time from snapmin to snapmax
                    try:
                            CenOut = FindEnzoCentroids(pf) # CenOut = [NumPart, Masses, Centroids]
                            BulletDMPos = CenOut[2][0]
                            MainDMPos = CenOut[2][1]
                            dt = simtime - lastsimtime
                            lastsimtime = simtime
                            BulletDMVel = (BulletDMPos-LastBulletDMPos)/dt
                            MainDMVel = (MainDMPos-LastMainDMPos)/dt
                            LastBulletDMPos = BulletDMPos 
                            LastMainDMPos = MainDMPos 
                            RelVel = BulletDMVel - MainDMVel
                            if snap > snapmin: # RelVel not valid for first snap
                                    cenfile = open(cenfilename, 'a')
                                    result = "Time = %.3f BulletDMPos = %.2f %.2f %.2f MainDMPos = %.2f %.2f %.2f RelVel = %.2f %.2f %.2f\n"%(simtime, BulletDMPos[0], BulletDMPos[1], BulletDMPos[2], MainDMPos[0], MainDMPos[1], MainDMPos[2], RelVel[0], RelVel[1], RelVel[2])
                                    print result
                                    sys.stdout.flush()
                                    cenfile.write(result)
                                    cenfile.close()
                            snap = snap + 1
                            continue
                    except:
			snap = snap + 1
                        print 'Centroid finding failed! snap = %d\n'%snap
                        sys.stdout.flush()
			continue

		if TimeStripeDone and not ThetaStripeDone and Theta >= ThetaMax: # Then, complete a theta stripe
			ThetaStripeDone = True
			print "Theta Stripe Finished \n"
                        sys.stdout.flush()
		if TimeStripeDone and not ThetaStripeDone: # Then, complete a theta stripe
			Theta = Theta + DeltaTheta
			continue

		# After completing both stripes, find the best fom so far and go to this point:
		for sn in range(snapmin-1,snapmax+1): 
			for Ps in range(PsiMin-DeltaPsi, PsiMax+DeltaPsi, DeltaPsi):	
				for Th in range(ThetaMin-DeltaTheta, ThetaMax+DeltaTheta, DeltaTheta):
					if FomMatrix[sn,Ps,Th][0] < fom and FomMatrix[sn,Ps,Th][0] > 0:
						fom = FomMatrix[sn,Ps,Th][0]
						snap = sn
                                                Phi = FomMatrix[sn,Ps,Th][1]
						Theta = Th
						Psi = Ps

		# If this point is better than all points around it (within tolerance), we're done
                print "Finished both stripes, fom = %f, snap = %d, Phi = %f, Theta = %d, Psi = %d"%(fom,snap,Phi,Theta,Psi)
                sys.stdout.flush()
		FomMatTol = FomMatrix[snap,Psi,Theta][0] * 0.99
		if  FomMatTol  < FomMatrix[snap+1,Psi,Theta][0] and FomMatTol < FomMatrix[snap-1,Psi,Theta][0] \
		and FomMatTol < FomMatrix[snap,Psi+DeltaPsi,Theta][0] and FomMatTol < FomMatrix[snap,Psi-DeltaPsi,Theta][0] \
		and FomMatTol < FomMatrix[snap,Psi,Theta+DeltaTheta][0] and FomMatTol < FomMatrix[snap,Psi,Theta-DeltaTheta][0]:
			simtime = pf.h.parameters['InitialTime'] * BulletConstants.TimeConversion # Put time in Gyears
			return (fom,snap,Phi,Theta/100.0,Psi/100.0,simtime,counter)

		# If not, find a surrounding point not yet done and run it.
		for (sn,Ps,Th) in [(snap-1,Psi,Theta),(snap+1,Psi,Theta),(snap,Psi+DeltaPsi,Theta),(snap,Psi-DeltaPsi,Theta),\
				(snap,Psi,Theta+DeltaTheta),(snap,Psi,Theta-DeltaTheta)]:
			if FomMatrix[sn,Ps,Th][0] > 999.0: # This means it's a boundary
				continue
			elif FomMatrix[sn,Ps,Th][0] < 0.0: # This means it hasn't yet been run, so move there and run it.
				snap = sn
				Psi = Ps
				Theta = Th
				print "New snap = %d, New Psi = %.2f, New Theta = %.2f"%(snap,Psi,Theta)
                                sys.stdout.flush()
				break
		
	return (1.0E5,0,0.0,0.0,0.0,0.0,0) # Returns garbage if it fails to find an optimum


def GetPF(snap):
	# Retrieves a yt - pf file given a snapshot number
	if snap<10:
	    filename = "DD000"+str(snap)+"/output_000"+str(snap) 
	elif snap<100:
	    filename = "DD00"+str(snap)+"/output_00"+str(snap) 
	elif snap<1000:
	    filename = "DD0"+str(snap)+"/output_0"+str(snap) 
	else:
	    filename = "DD"+str(snap)+"/output_"+str(snap) 
	
	pf = load(filename)
	return pf

def GetKernels():
        # This reads in kernels used for image manipulation

	toppath = lageconfig.toppath
	datapathszekernel=toppath+'bullet/data/sze/bullet_transfer_rescaled_small.fits'
	datapathkernel60=toppath+'bullet/data/kernel60.dat'
	datapathkernel15=toppath+'bullet/data/kernel15.dat'
	szekernel=Array2d(0.0,400.0,113,0.0,400.0,113)
	szekernel=GetBulletFits(datapathszekernel,szekernel)

	kernel60=Array2d(0.0,1.0,110,0.0,1.0,110)
	kernel60=GetBulletData(datapathkernel60,kernel60)
	kernel15=Array2d(0.0,1.0,110,0.0,1.0,110)
	kernel15=GetBulletData(datapathkernel15,kernel15)
        return (szekernel, kernel60, kernel15)

def GetData():
        # This reads in the image data
	toppath = lageconfig.toppath
	datapathA=toppath+'bullet/data/kappa_25Apr12.dat'
	datapathB1=toppath+'bullet/data/xray_500_2000_29Jun11.dat'
	datapathB2=toppath+'bullet/data/xray_2000_5000_29Jun11.dat'
	datapathB3=toppath+'bullet/data/xray_5000_8000_29Jun11.dat'
	datapathC=toppath+'bullet/data/sze_data.dat'
	datapathD=toppath+'bullet/data/xray_temp_23Jun11.dat'

	datapathsA=toppath+'bullet/data/mass_sigma_11Feb13.dat'
	datapathsB1=toppath+'bullet/data/xray_sigma_500_2000_14Jul11.dat'
	datapathsB2=toppath+'bullet/data/xray_sigma_2000_5000_14Jul11.dat'
	datapathsB3=toppath+'bullet/data/xray_sigma_5000_8000_14Jul11.dat'
	datapathsC=toppath+'bullet/data/sze_sigma.dat'
	datapathsD=toppath+'bullet/data/xray_temp_sigma_23Jun11.dat'

	datapathmA=toppath+'bullet/data/mask_3May12.dat'
	datapathmANull=toppath+'bullet/data/mask_23Jun11_null.dat'
	
	datapathmD=toppath+'bullet/data/xtemp_mask_23Jun11.dat'
	datapathmDNull=toppath+'bullet/data/xtemp_mask_null_29Jun11.dat' # Ignore temp in this case.
	
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

	dataA=Array2d(dxminA,dxmaxA,xpixelsA,dyminA,dymaxA,ypixelsA)
	dataA=GetBulletData(datapathA,dataA)# Mass density measured data
	sigmaA=Array2d(dxminA,dxmaxA,xpixelsA,dyminA,dymaxA,ypixelsA)
	sigmaA=GetBulletData(datapathsA,sigmaA)#  Mass sigma
	maskA=Array2d(dxminA,dxmaxA,xpixelsA,dyminA,dymaxA,ypixelsA)
	maskA=GetBulletData(datapathmA,maskA)#  Mask
	maskANull=Array2d(dxminA,dxmaxA,xpixelsA,dyminA,dymaxA,ypixelsA)
	maskANull=GetBulletData(datapathmANull,maskANull)#  Null Mask

	xscaleB=yscaleB=BulletConstants.AngularScale*3600*BulletConstants.XRayDataScale 
	# Pixel size in kpc. 4.413 kpc/" is the angular scale at the bullet cluster
	# 3600 "/degree, 9.86E-4 degrees/pixel is the data scale

	xpixelsB=ypixelsB=110
	sxpixelsB=sypixelsB=220
	dxmaxB = xscaleB*xpixelsB/2
	dymaxB = yscaleB*ypixelsB/2
	sxmaxB = dxmaxB*2
	symaxB = dymaxB*2
	dataB1=Array2d(-dxmaxB,dxmaxB,xpixelsB,-dymaxB,dymaxB,ypixelsB)
	dataB1=GetBulletData(datapathB1,dataB1)# XRay measured data
	sigmaB1=Array2d(-dxmaxB,dxmaxB,xpixelsB,-dymaxB,dymaxB,ypixelsB)
	sigmaB1=GetBulletData(datapathsB1,sigmaB1)# XRay sigma
	dataB2=Array2d(-dxmaxB,dxmaxB,xpixelsB,-dymaxB,dymaxB,ypixelsB)
	dataB2=GetBulletData(datapathB2,dataB2)# XRay measured data
	sigmaB2=Array2d(-dxmaxB,dxmaxB,xpixelsB,-dymaxB,dymaxB,ypixelsB)
	sigmaB2=GetBulletData(datapathsB2,sigmaB2)# XRay sigma
	dataB3=Array2d(-dxmaxB,dxmaxB,xpixelsB,-dymaxB,dymaxB,ypixelsB)
	dataB3=GetBulletData(datapathB3,dataB3)# XRay measured data
	sigmaB3=Array2d(-dxmaxB,dxmaxB,xpixelsB,-dymaxB,dymaxB,ypixelsB)
	sigmaB3=GetBulletData(datapathsB3,sigmaB3)# XRay sigma

	#dataB.data=gaussian_filter(dataB.data,0.5) 
	# 0.5-sigma smoothing of X-Ray data
	dataC=Array2d(dxminA,dxmaxA,xpixelsA,dyminA,dymaxA,ypixelsA)
	dataC=GetBulletData(datapathC,dataC)# Measured SZE data
	sigmaC=Array2d(dxminA,dxmaxA,xpixelsA,dyminA,dymaxA,ypixelsA)
	sigmaC=GetBulletData(datapathsC,sigmaC)# SZE Sigma

	dataC.data = abs(dataC.data*1E6) #convert to micro Kelvin
	sigmaC.data = sigmaC.data*1E6 #convert to micro Kelvin

	xscaleD=yscaleD=BulletConstants.AngularScale*3600*BulletConstants.XRayTempDataScale 
	# Pixel size in kpc. 4.413 kpc/" is the angular scale at the bullet cluster
	# 3600 "/degree, 1.09333 degrees/pixel is the data scale

	xpixelsD=ypixelsD=100
	sxpixelsD=sypixelsD=200
	dxmaxD = xscaleD*xpixelsD/2
	dymaxD = yscaleD*ypixelsD/2
	sxmaxD = dxmaxD*2
	symaxD = dymaxD*2
	dataD=Array2d(-dxmaxD,dxmaxD,xpixelsD,-dymaxD,dymaxD,ypixelsD)
	dataD=GetBulletData(datapathD,dataD)# XRay measured data

	sigmaD=Array2d(-dxmaxD,dxmaxD,xpixelsD,-dymaxD,dymaxD,ypixelsD)
	sigmaD=GetBulletData(datapathsD,sigmaD)# Xray Temp Sigma
	maskD=Array2d(-dxmaxD,dxmaxD,xpixelsD,-dymaxD,dymaxD,ypixelsD)
	maskD=GetBulletData(datapathmD,maskD)#  Mask
	maskDNull=Array2d(-dxmaxD,dxmaxD,xpixelsD,-dymaxD,dymaxD,ypixelsD)
	maskDNull=GetBulletData(datapathmDNull,maskDNull)#  Null Mask

        return (dataA,sigmaA,maskA,maskANull,dataB1,sigmaB1,dataB2,sigmaB2,dataB3,sigmaB3,dataC,sigmaC,dataD,sigmaD,maskD,maskDNull)


def SetAlign(dataB1, phi, theta, psi, ConstrainPhi = False):
        # This sets the alignment parameters
	align=Align()
	align.dmax[0]=(dataB1.xmax)
	align.dmax[1]=(dataB1.ymax)
	align.dmin[0]=-align.dmax[0]
	align.dmin[1]=-align.dmax[1]
	align.dmin[2]=0.0
	align.dmax[2]=2.0*pi
	
        align.d[0] = 0.0
        align.d[1] = 0.0

        if ConstrainPhi:
                align.dmin[0] = -150.0
                align.dmax[0] =  150.0
                align.dmin[1] = -150.0
                align.dmax[1] =  150.0

        if cos(psi)>=0.0:
                align.d[2] = 0.22  - arctan( cos(theta)*tan(psi) ) # Seed phi close to final result
                if ConstrainPhi:
                        align.dmin[2] = 0.12  - arctan( cos(theta)*tan(psi) )
                        align.dmax[2] = 0.32  - arctan( cos(theta)*tan(psi) )
        else:
                align.d[2] = 0.22  + pi - arctan( cos(theta)*tan(psi) ) # Seed phi close to final result
                if ConstrainPhi:
                        align.dmin[2] = 0.12  + pi - arctan( cos(theta)*tan(psi) ) # Seed psi close to final result
                        align.dmax[2] = 0.32  + pi - arctan( cos(theta)*tan(psi) ) # Seed psi close to final result

        return align




def BestPlot(snap,Z,phi=0.0,theta=0.0,psi=0.0,PlotSuffix='New_Kappa',FomLocate=False,TFudge=1.0,ConstrainPhi=False, Mask=(1,0,0,0,0,0)):

        (dataA,sigmaA,maskA,maskANull,dataB1,sigmaB1,dataB2,sigmaB2,dataB3,sigmaB3,dataC,sigmaC,dataD,sigmaD,maskD,maskDNull) = GetData()
        (szekernel, kernel60, kernel15) = GetKernels()
	[dyyA,dxxA] = meshgrid(dataA.y,dataA.x)# Data grid for plots
	[dyyB,dxxB] = meshgrid(dataB1.y,dataB1.x)# Data grid for plots
	[dyyD,dxxD] = meshgrid(dataD.y,dataD.x)# Data grid for plots

        align = SetAlign(dataB1, phi, theta, psi, ConstrainPhi = ConstrainPhi)
	pp=PdfPages('Graph_'+PlotSuffix+'.pdf')
	pf = GetPF(snap)
	simtime = pf.h.parameters['InitialTime'] * BulletConstants.TimeConversion # Put time in Gyears

	dmsim=Array2d(2.0*dataA.xmin,2.0*dataA.xmax,2*dataA.nx,2.0*dataA.ymin,2.0*dataA.ymax,2*dataA.ny)
	masssim=Array2d(2.0*dataA.xmin,2.0*dataA.xmax,2*dataA.nx,2.0*dataA.ymin,2.0*dataA.ymax,2*dataA.ny)
	sumsim=Array2d(2.0*dataA.xmin,2.0*dataA.xmax,2*dataA.nx,2.0*dataA.ymin,2.0*dataA.ymax,2*dataA.ny)
	tempmass=Array2d(2.0*dataD.xmin,2.0*dataD.xmax,2*dataD.nx,2.0*dataD.ymin,2.0*dataD.ymax,2*dataD.ny)
	[syyA,sxxA] = meshgrid(dmsim.y,dmsim.x) # Sim grid for plots
	[syyD,sxxD] = meshgrid(tempmass.y,tempmass.x) # Sim grid for plots

	ApecData = ReadLookups(Z) # Reads the APEC lookup tables.
	for grid in pf.h.grids:
		grid.set_field_parameter('ApecData',ApecData)
                grid.set_field_parameter('TFudge',TFudge)

        DMProject = False
	[dmsim,masssim,xraysim1,xraysim2,xraysim3,szsim]=ProjectEnzoData(pf,masssim,phi,theta,psi,DMProject=DMProject)
	sumsim.data=gaussian_filter(dmsim.data+masssim.data,2.0)
	szsim.data=convolve(abs(szsim.data*1E6),szekernel.data,mode='constant',cval=0.0)#Smooth simulation with given transfer function 
	[tempsim,bmagsim]=ProjectEnzoTemp(pf,tempmass,phi,theta,psi)
	# Projected Enzo sim data

	tol=1.0E-7 # Alignment tolerance

	data1list=list((sumsim,xraysim1,xraysim2,xraysim3,szsim,tempsim))
	data2list=list((dataA,dataB1,dataB2,dataB3,dataC,dataD))
	sigmalist=list((sigmaA,sigmaB1,sigmaB2,sigmaB3,sigmaC,sigmaD))
        masklist=list()
        # This code appends the necessary masks. If Mask[i]==0, this data is not included in the FOM
        for i in range(5):
            if Mask[i] == 0:
                masklist.append((maskANull))
                print "i = %d, Mask[i] = %d, appending maskAnull\n"%(i,Mask[i])
            else:
                masklist.append((maskA))
                print "i = %d, Mask[i] = %d, appending maskA\n"%(i,Mask[i])

        if Mask[5] == 0:
            masklist.append((maskDNull))
            print "i = 5, Mask[i] = %d, appending maskDNull\n"%Mask[i]
        else:
            masklist.append((maskD))
            print "i = 5, Mask[i] = %d, appending maskD\n"%Mask[i]

        sys.stdout.flush()
	[shifteddata1list,align,fom]=FindBestShift(data1list, data2list, sigmalist, masklist, align, tol)

	shiftedmsim=shifteddata1list[0]
	shiftedxsim1=shifteddata1list[1]
	shiftedxsim2=shifteddata1list[2]
	shiftedxsim3=shifteddata1list[3]
	shiftedszsim=shifteddata1list[4]
	shiftedtempsim=shifteddata1list[5]
	try:
	    dmlevels=linspace(0,5.0,25) 
	    for i in range(25): dmlevels[i]=int(dmlevels[i]*10)/10.0

	    xlevels=linspace(0.0,10.0,25) 
	    for i in range(25): xlevels[i]=int(xlevels[i]*10)/10.0

	    # First the mass plots
	    figure()
	    suptitle('Bullet Cluster Collision, T = %.2f Gy'%simtime, fontsize=16)
	    subplots_adjust(hspace=0.4, wspace=0.4)

	    subplot("221",aspect=1)
	    title("Mass data")
	    contourf(dxxA,dyyA,dataA.data)
	    for c in contourf(dxxA,dyyA,dataA.data).collections:
		c.set_linewidth(0.1)
		c.set_alpha(1.0)
	    gca().set_xticks([-800.0,-400.0,0.0,800.0,400.0])
	    colorbar()

	    subplot("222",aspect=1)
	    title("Mass simulation")
	    contourf(dxxA,dyyA,shiftedmsim.data)
	    for c in contourf(dxxA,dyyA,shiftedmsim.data).collections:
		c.set_linewidth(0.1)
		c.set_alpha(1.0)
	    colorbar()
	    gca().set_xticks([-800.0,-400.0,0.0,800.0,400.0])
	    
	    [dataslice,simslice] = OneDSlice2Fixed(dataA,shiftedmsim,0.26,80,53)

	    subplot("223",aspect=1)

	    title('Overlay: $\chi^2$ = %.2f '%fom)

	    contourf(dxxA,dyyA,dataA.data)
	    for c in contourf(dxxA,dyyA,dataA.data).collections:
		c.set_linewidth(0.1)
		c.set_alpha(1.0)
	    contour(dxxA,dyyA,shiftedmsim.data, colors='k')
	    plot(dataslice.x,dataslice.y,'w:')# Plot the 1D slice
	    gca().set_xticks([-800.0,-400.0,0.0,800.0,400.0])
	    xlabel('kpc')
	    ylabel('kpc')

	    subplot("224")
	    title('1D Slice through data max')
	    plot(dataslice.x,dataslice.data,label="Data",linewidth=3.0)
	    plot(simslice.x,simslice.data,label="Sim",linewidth=3.0)

	    xlabel('kpc')
	    legend(loc=1,bbox_to_anchor=[1.3, 1.05])
	    ltext = gca().get_legend().get_texts()
	    setp(ltext, fontsize = 9, color = 'k')
	    pp.savefig()
	    
	    # Next the Xray plots 1 

	    figure()
	    suptitle('Bullet Cluster Collision, T = %.2f Gy'%simtime, fontsize=16)
	    subplots_adjust(hspace=0.4, wspace=0.4)

	    subplot("221",aspect=1)
	    title("Log XRay data-500-2000eV")
	    XrayLogMin = -10.0
	    XrayLogMax = -4.0
	    Xraylevels=linspace(XrayLogMin,XrayLogMax,25) 
	    for i in range(25): Xraylevels[i]=int(Xraylevels[i]*10)/10.0
	    XraylevelsCon=linspace(XrayLogMin,XrayLogMax,10) 
	    contourf(dxxB,dyyB,log10(dataB1.data),Xraylevels,cmap=cm.spectral)
	    for c in contourf(dxxB,dyyB,log10(dataB1.data),Xraylevels,cmap=cm.spectral).collections:
		c.set_linewidth(0.1)
		c.set_alpha(1.0)
	    gca().set_xticks([-800.0,-400.0,0.0,800.0,400.0])
	    colorbar()

	    subplot("222",aspect=1)
	    title("Log XRay simulation-500-2000eV")
	    contourf(dxxB,dyyB,log10(shiftedxsim1.data),Xraylevels,cmap=cm.spectral)
	    for c in contourf(dxxB,dyyB,log10(shiftedxsim1.data),Xraylevels,cmap=cm.spectral).collections:
		c.set_linewidth(0.1)
		c.set_alpha(1.0)
	    colorbar()
	    gca().set_xticks([-800.0,-400.0,0.0,800.0,400.0])

	    [dataslice,simslice] = OneDSlice2Fixed(dataB1,shiftedxsim1,0.24,70,55)

	    subplot("223",aspect=1)

	    title('Overlay: $\chi^2$ = %.2f '%fom)

	    contourf(dxxB,dyyB,log10(dataB1.data),Xraylevels,cmap=cm.spectral)
	    for c in contourf(dxxB,dyyB,log10(dataB1.data),Xraylevels,cmap=cm.spectral).collections:
		c.set_linewidth(0.1)
		c.set_alpha(1.0)
	    contour(dxxB,dyyB,log10(shiftedxsim1.data),XraylevelsCon,linestyle='-',colors='k')
	    plot(dataslice.x,dataslice.y,'w:')# Plot the 1D slice
	    gca().set_xticks([-800.0,-400.0,0.0,800.0,400.0])
	    xlabel('kpc')
	    ylabel('kpc')

	    subplot("224")
	    title('1D Slice through data max')
	    plot(dataslice.x,dataslice.data,label="Data",linewidth=3.0)
	    plot(simslice.x,simslice.data,label="Sim",linewidth=3.0)

	    xlabel('kpc')
	    legend(loc=1,bbox_to_anchor=[0.45, 1.05])
	    ltext = gca().get_legend().get_texts()
	    setp(ltext, fontsize = 9, color = 'k')
	    pp.savefig()
	    
	    # Next the Xray plots 2 

	    figure()
	    suptitle('Bullet Cluster Collision, T = %.2f Gy'%simtime, fontsize=16)
	    subplots_adjust(hspace=0.4, wspace=0.4)

	    subplot("221",aspect=1)
	    title("Log XRay data-2000-5000eV")
	    XrayLogMin = -10.0
	    XrayLogMax = -4.0
	    Xraylevels=linspace(XrayLogMin,XrayLogMax,25) 
	    for i in range(25): Xraylevels[i]=int(Xraylevels[i]*10)/10.0
	    XraylevelsCon=linspace(XrayLogMin,XrayLogMax,10) 
	    contourf(dxxB,dyyB,log10(dataB2.data),Xraylevels,cmap=cm.spectral)
	    for c in contourf(dxxB,dyyB,log10(dataB2.data),Xraylevels,cmap=cm.spectral).collections:
		c.set_linewidth(0.1)
		c.set_alpha(1.0)
	    gca().set_xticks([-800.0,-400.0,0.0,800.0,400.0])
	    colorbar()

	    subplot("222",aspect=1)
	    title("Log XRay simulation-2000-5000eV")
	    contourf(dxxB,dyyB,log10(shiftedxsim2.data),Xraylevels,cmap=cm.spectral)
	    for c in contourf(dxxB,dyyB,log10(shiftedxsim2.data),Xraylevels,cmap=cm.spectral).collections:
		c.set_linewidth(0.1)
		c.set_alpha(1.0)
	    colorbar()
	    gca().set_xticks([-800.0,-400.0,0.0,800.0,400.0])

	    [dataslice,simslice] = OneDSlice2Fixed(dataB2,shiftedxsim2,0.24,70,55)

	    subplot("223",aspect=1)

	    title('Overlay: $\chi^2$ = %.2f '%fom)

	    contourf(dxxB,dyyB,log10(dataB2.data),Xraylevels,cmap=cm.spectral)
	    for c in contourf(dxxB,dyyB,log10(dataB2.data),Xraylevels,cmap=cm.spectral).collections:
		c.set_linewidth(0.1)
		c.set_alpha(1.0)
	    contour(dxxB,dyyB,log10(shiftedxsim2.data),XraylevelsCon,linestyle='-',colors='k')
	    plot(dataslice.x,dataslice.y,'w:')# Plot the 1D slice
	    gca().set_xticks([-800.0,-400.0,0.0,800.0,400.0])
	    xlabel('kpc')
	    ylabel('kpc')

	    subplot("224")
	    title('1D Slice through data max')
	    plot(dataslice.x,dataslice.data,label="Data",linewidth=3.0)
	    plot(simslice.x,simslice.data,label="Sim",linewidth=3.0)

	    xlabel('kpc')
	    legend(loc=1,bbox_to_anchor=[0.45, 1.05])
	    ltext = gca().get_legend().get_texts()
	    setp(ltext, fontsize = 9, color = 'k')
	    pp.savefig()

	    # Next the Xray plots 3

	    figure()
	    suptitle('Bullet Cluster Collision, T = %.2f Gy'%simtime, fontsize=16)
	    subplots_adjust(hspace=0.4, wspace=0.4)

	    subplot("221",aspect=1)
	    title("Log XRay data-5000-8000eV")
	    XrayLogMin = -10.0
	    XrayLogMax = -4.0
	    Xraylevels=linspace(XrayLogMin,XrayLogMax,25) 
	    for i in range(25): Xraylevels[i]=int(Xraylevels[i]*10)/10.0
	    XraylevelsCon=linspace(XrayLogMin,XrayLogMax,10) 
	    contourf(dxxB,dyyB,log10(dataB3.data),Xraylevels,cmap=cm.spectral)
	    for c in contourf(dxxB,dyyB,log10(dataB3.data),Xraylevels,cmap=cm.spectral).collections:
		c.set_linewidth(0.1)
		c.set_alpha(1.0)
	    gca().set_xticks([-800.0,-400.0,0.0,800.0,400.0])
	    colorbar()

	    subplot("222",aspect=1)
	    title("Log XRay simulation-5000-8000eV")
	    contourf(dxxB,dyyB,log10(shiftedxsim3.data),Xraylevels,cmap=cm.spectral)
	    for c in contourf(dxxB,dyyB,log10(shiftedxsim3.data),Xraylevels,cmap=cm.spectral).collections:
		c.set_linewidth(0.1)
		c.set_alpha(1.0)
	    colorbar()
	    gca().set_xticks([-800.0,-400.0,0.0,800.0,400.0])

	    [dataslice,simslice] = OneDSlice2Fixed(dataB3,shiftedxsim3,0.24,70,55)

	    subplot("223",aspect=1)

	    title('Overlay: $\chi^2$ = %.2f '%fom)

	    contourf(dxxB,dyyB,log10(dataB3.data),Xraylevels,cmap=cm.spectral)
	    for c in contourf(dxxB,dyyB,log10(dataB3.data),Xraylevels,cmap=cm.spectral).collections:
		c.set_linewidth(0.1)
		c.set_alpha(1.0)
	    contour(dxxB,dyyB,log10(shiftedxsim3.data),XraylevelsCon,linestyle='-',colors='k')
	    plot(dataslice.x,dataslice.y,'w:')# Plot the 1D slice
	    gca().set_xticks([-800.0,-400.0,0.0,800.0,400.0])
	    xlabel('kpc')
	    ylabel('kpc')

	    subplot("224")
	    title('1D Slice through data max')
	    plot(dataslice.x,dataslice.data,label="Data",linewidth=3.0)
	    plot(simslice.x,simslice.data,label="Sim",linewidth=3.0)

	    xlabel('kpc')
	    legend(loc=1,bbox_to_anchor=[0.45, 1.05])
	    ltext = gca().get_legend().get_texts()
	    setp(ltext, fontsize = 9, color = 'k')
	    pp.savefig()

	    # Next the Xray1 Shock plots 

	    figure()
	    suptitle('Bullet Cluster Collision, T = %.2f Gy'%simtime, fontsize=16)
	    subplots_adjust(hspace=0.4, wspace=0.4)


	    theta1 = -0.54
	    theta2 = 0.24
	    theta3 = 1.02

	    [dataslice1,simslice1] = OneDSlice2Fixed(dataB1,shiftedxsim1,theta1,78,56)
	    [dataslice2,simslice2] = OneDSlice2Fixed(dataB1,shiftedxsim1,theta2,78,56)
	    [dataslice3,simslice3] = OneDSlice2Fixed(dataB1,shiftedxsim1,theta3,78,56)

	    subplot("221",aspect=1)

	    title('Overlay: $\chi^2$ = %.2f '%fom)

	    contourf(dxxB,dyyB,log10(dataB1.data),Xraylevels,cmap=cm.spectral)
	    for c in contourf(dxxB,dyyB,log10(dataB1.data),Xraylevels,cmap=cm.spectral).collections:
		c.set_linewidth(0.1)
		c.set_alpha(1.0)
	    contour(dxxB,dyyB,log10(shiftedxsim1.data),XraylevelsCon,linestyle='-',colors='k')
	    plot(dataslice1.x,dataslice1.y,'w:')# Plot the 1D slice
	    plot(dataslice2.x,dataslice2.y,'y:')# Plot the 1D slice
	    plot(dataslice3.x,dataslice3.y,'k:')# Plot the 1D slice
	    gca().set_xticks([-800.0,-400.0,0.0,800.0,400.0])
	    xlabel('kpc')
	    ylabel('kpc')

	    subplot("222")
	    title('1D Slice through data max')
	    plot(dataslice1.x,log10(dataslice1.data),label="Data1",linewidth=3.0)
	    plot(simslice1.x,log10(simslice1.data),label="Sim1",linewidth=3.0)

	    xlabel('kpc')
	    legend(loc=1,bbox_to_anchor=[0.45, 1.05])
	    ltext = gca().get_legend().get_texts()
	    setp(ltext, fontsize = 9, color = 'k')

	    subplot("223")
	    title('1D Slice through data max')
	    plot(dataslice2.x,log10(dataslice2.data),label="Data2",linewidth=3.0)
	    plot(simslice2.x,log10(simslice2.data),label="Sim2",linewidth=3.0)

	    xlabel('kpc')
	    legend(loc=1,bbox_to_anchor=[0.45, 1.05])
	    ltext = gca().get_legend().get_texts()
	    setp(ltext, fontsize = 9, color = 'k')

	    subplot("224")
	    title('1D Slice through data max')
	    plot(dataslice3.x,log10(dataslice3.data),label="Data3",linewidth=3.0)
	    plot(simslice3.x,log10(simslice3.data),label="Sim3",linewidth=3.0)

	    xlabel('kpc')
	    legend(loc=1,bbox_to_anchor=[0.45, 1.05])
	    ltext = gca().get_legend().get_texts()
	    setp(ltext, fontsize = 9, color = 'k')
	    pp.savefig()


	    # Play with edge detection plots 

	    figure()
	    suptitle('Bullet Cluster Collision, T = %.2f Gy'%simtime, fontsize=16)
	    subplots_adjust(hspace=0.4, wspace=0.4)
	    subplot("221",aspect=1)
	    title("Log XRay data-500-2000eV")
	    XrayLogMin = -10.0
	    XrayLogMax = -4.0
	    Xraylevels=linspace(XrayLogMin,XrayLogMax,25) 
	    for i in range(25): Xraylevels[i]=int(Xraylevels[i]*10)/10.0
	    XraylevelsCon=linspace(XrayLogMin,XrayLogMax,10) 
	    contourf(dxxB,dyyB,log10(dataB1.data),Xraylevels,cmap=cm.spectral)
	    for c in contourf(dxxB,dyyB,log10(dataB1.data),Xraylevels,cmap=cm.spectral).collections:
		c.set_linewidth(0.1)
		c.set_alpha(1.0)
	    gca().set_xticks([-800.0,-400.0,0.0,800.0,400.0])
	    colorbar()

	    dataedge=Array2d(dataB1.xmin,dataB1.xmax,dataB1.nx,dataB1.ymin,dataB1.ymax,dataB1.ny)
	    tempsimedge=Array2d(dataD.xmin,dataD.xmax,dataD.nx,dataD.ymin,dataD.ymax,dataD.ny)
	    dataBtrunc=Array2d(dataB1.xmin,dataB1.xmax,dataB1.nx,dataB1.ymin,dataB1.ymax,dataB1.ny)

	    xraymax = pow(10, -6.5)
	    xraymin = pow(10, -8.0)
	    for ii in range(dataB1.nx):
		for jj in range(dataB1.ny):
			dataBtrunc.data[ii,jj] = max(min(xraymax, dataB1.data[ii,jj]),xraymin)

	    dataedge.data = hypot(sobel(dataBtrunc.data, 0),sobel(dataBtrunc.data, 1))

	    #dataedge.data = hypot(sobel(dataB1.data, 0),sobel(dataB1.data, 1))
	    tempsimedge.data = hypot(sobel(shiftedtempsim.data, 0),sobel(shiftedtempsim.data, 1))

	    subplot("222",aspect=1)
	    title("Data Edge Detection - 500-2000eV")
	    contourf(dxxB,dyyB,dataedge.data,cmap=cm.spectral)
	    for c in contourf(dxxB,dyyB,dataedge.data,cmap=cm.spectral).collections:
		c.set_linewidth(0.1)
		c.set_alpha(1.0)
	    colorbar()
	    gca().set_xticks([-800.0,-400.0,0.0,800.0,400.0])

	    subplot("223",aspect=1)
	    title("Sim Edge Detection - Temp")
	    contourf(dxxD,dyyD,tempsimedge.data,cmap=cm.spectral)
	    for c in contourf(dxxD,dyyD,tempsimedge.data,cmap=cm.spectral).collections:
		c.set_linewidth(0.1)
		c.set_alpha(1.0)
	    colorbar()
	    gca().set_xticks([-800.0,-400.0,0.0,800.0,400.0])

	    subplot("224",aspect=1)

	    title('Overlay: $\chi^2$ = %.2f '%fom)

	    contourf(dxxB,dyyB,dataedge.data,cmap=cm.spectral)
	    for c in contourf(dxxB,dyyB,dataedge.data,cmap=cm.spectral).collections:
		c.set_linewidth(0.1)
		c.set_alpha(1.0)
	    contour(dxxD,dyyD,tempsimedge.data,linestyle='-',colors='k')
	    gca().set_xticks([-800.0,-400.0,0.0,800.0,400.0])
	    xlabel('kpc')
	    ylabel('kpc')
	    pp.savefig()


	    # Next the Temp plots
	    
	    figure()
	    suptitle('Bullet Cluster Collision, T = %.2f Gy'%simtime, fontsize=16)
	    subplots_adjust(hspace=0.4, wspace=0.4)
	    		    
	    subplot("221",aspect=1)
	    title("Temp data")
	    TempMin = 0.0
	    TempMax = 20.0
	    Templevels=linspace(TempMin,TempMax,25) 
	    for i in range(25): Templevels[i]=int(Templevels[i]*10)/10.0
	    TemplevelsCon=linspace(TempMin,TempMax,10) 
	    contourf(dxxD,dyyD,dataD.data,Templevels,cmap=cm.spectral)
	    for c in contourf(dxxD,dyyD,dataD.data,Templevels,cmap=cm.spectral).collections:
		c.set_linewidth(0.1)
		c.set_alpha(1.0)
	    gca().set_xticks([-800.0,-400.0,0.0,800.0,400.0])
	    colorbar()
	    
	    subplot("222",aspect=1)
	    title("Temp simulation")
	    contourf(dxxD,dyyD,shiftedtempsim.data,Templevels,cmap=cm.spectral)
	    for c in contourf(dxxD,dyyD,shiftedtempsim.data,Templevels,cmap=cm.spectral).collections:
		c.set_linewidth(0.1)
		c.set_alpha(1.0)
	    colorbar()
	    gca().set_xticks([-800.0,-400.0,0.0,800.0,400.0])

	    [dataslice,simslice] = OneDSlice2Fixed(dataD,shiftedtempsim,0.24,70,52)

	    subplot("223",aspect=1)

	    title('Overlay: $\chi^2$ = %.2f '%fom)

	    contourf(dxxD,dyyD,dataD.data,Templevels,cmap=cm.spectral)
	    for c in contourf(dxxD,dyyD,dataD.data,Templevels,cmap=cm.spectral).collections:
		c.set_linewidth(0.1)
		c.set_alpha(1.0)
	    contour(dxxD,dyyD,shiftedtempsim.data,TemplevelsCon,linestyle='-',colors='k')
	    plot(dataslice.x,dataslice.y,'w:')# Plot the 1D slice
	    gca().set_xticks([-800.0,-400.0,0.0,800.0,400.0])
	    xlabel('kpc')
	    ylabel('kpc')
	    
	    subplot("224")
	    title('1D Slice through data max')
	    plot(dataslice.x,dataslice.data,label="Data",linewidth=3.0)
	    plot(simslice.x,simslice.data,label="Sim",linewidth=3.0)
	    gca().set_xticks([-800.0,-400.0,0.0,800.0,400.0])
	    #ylim(0.0,20.0)
	    xlabel('kpc')
	    legend(loc=1,bbox_to_anchor=[1.3, 1.05])
	    ltext = gca().get_legend().get_texts()
	    setp(ltext, fontsize = 9, color = 'k')	    
	    pp.savefig()
	    
	    # Now the SZE plots
	    figure()
	    suptitle('Bullet Cluster Collision, T = %.2f Gy'%simtime, fontsize=16)
	    subplots_adjust(hspace=0.4, wspace=0.4)

	    SZLevels=linspace(0,400,9)

	    subplot("221",aspect=1)
	    title("SZE data ($\mu K$)")
	    contourf(dxxA,dyyA,dataC.data,SZLevels)
	    for c in contourf(dxxA,dyyA,dataC.data,SZLevels).collections:
		c.set_linewidth(0.1)
		c.set_alpha(1.0)
	    gca().set_xticks([-800.0,-400.0,0.0,800.0,400.0])
	    colorbar()
	    
	    subplot("222",aspect=1)
	    title("SZE simulation ($\mu K$)")
	    contourf(dxxA,dyyA,shiftedszsim.data,SZLevels)
	    for c in contourf(dxxA,dyyA,shiftedszsim.data,SZLevels).collections:
		c.set_linewidth(0.1)
		c.set_alpha(1.0)
	    colorbar()
	    gca().set_xticks([-800.0,-400.0,0.0,800.0,400.0])

	    [dataslice,simslice] = OneDSlice2Fixed(dataC,shiftedszsim,0.24,70,54)

	    subplot("223",aspect=1)

	    title('Overlay: $\chi^2$ = %.2f '%fom)

	    contourf(dxxA,dyyA,dataC.data,SZLevels)
	    for c in contourf(dxxA,dyyA,dataC.data,SZLevels).collections:
		c.set_linewidth(0.1)
		c.set_alpha(1.0)
	    contour(dxxA,dyyA,shiftedszsim.data,SZLevels, colors='k')
	    plot(dataslice.x,dataslice.y,'w:')# Plot the 1D slice
	    gca().set_xticks([-800.0,-400.0,0.0,800.0,400.0])
	    xlabel('kpc')
	    ylabel('kpc')

	    subplot("224")
	    title('1D Slice through data max ($\mu K$)')
	    plot(dataslice.x,dataslice.data,label="Data",linewidth=3.0)
	    plot(simslice.x,simslice.data,label="Sim",linewidth=3.0)

	    xlabel('kpc')
	    legend(loc=1,bbox_to_anchor=[1.3, 1.05])
	    ltext = gca().get_legend().get_texts()
	    setp(ltext, fontsize = 9, color = 'k')

	    pp.savefig()
	    
	    # FOM location plots
	    if FomLocate:

		    numpixels = dataA.nx * dataA.ny
		    masschisquared = pow((dataA.data-shiftedmsim.data)/sigmaA.data,2) * maskA.data
		    xray1chisquared = pow((dataB1.data-shiftedxsim1.data)/sigmaB1.data,2) * maskA.data
		    xray2chisquared = pow((dataB2.data-shiftedxsim2.data)/sigmaB2.data,2) * maskA.data
		    xray3chisquared = pow((dataB3.data-shiftedxsim3.data)/sigmaB3.data,2) * maskA.data
		    szechisquared = pow((dataC.data-shiftedszsim.data)/sigmaC.data,2) * maskA.data
		    tempchisquared = pow((dataD.data-shiftedtempsim.data)/sigmaD.data,2) * maskD.data
		    masscontrib = masschisquared.sum()/(maskA.data.sum())
		    xray1contrib = xray1chisquared.sum()/(maskA.data.sum())
		    xray2contrib = xray2chisquared.sum()/(maskA.data.sum())
		    xray3contrib = xray3chisquared.sum()/(maskA.data.sum())
		    szecontrib = szechisquared.sum()/(maskA.data.sum())
		    tempcontrib = tempchisquared.sum()/(maskD.data.sum())

	   	    FomLevels=linspace(0,80,21)
		    
		    figure()
		    suptitle('Bullet Cluster Collision, T = %.2f Gy'%simtime, fontsize=16)
		    subplots_adjust(hspace=0.7, wspace=0.5)

		    subplot(2,3,1)
		    title("Mass - $\chi^2$ = %.3f"%masscontrib,fontsize=10)
		    contourf(dxxA,dyyA,masschisquared)
		    for c in contourf(dxxA,dyyA,masschisquared).collections:
			c.set_linewidth(0.1)
			c.set_alpha(1.0)
		    gca().set_xticks([-500.0,0.0,500.0])
		    xlabel('kpc')
		    ylabel('kpc')
		    colorbar()
		    
		    subplot(2,3,2)
		    title("XRay-500-2000eV - $\chi^2$ = %.3f"%xray1contrib,fontsize=10)
		    contourf(dxxA,dyyA,xray1chisquared)
		    for c in contourf(dxxA,dyyA,xray1chisquared).collections:
			c.set_linewidth(0.1)
			c.set_alpha(1.0)
		    gca().set_xticks([-500.0,0.0,500.0])
		    xlabel('kpc')
		    ylabel('kpc')
		    colorbar()

		    subplot(2,3,3)
		    title("XRay-2000-5000eV - $\chi^2$ = %.3f"%xray2contrib,fontsize=10)
		    contourf(dxxA,dyyA,xray2chisquared)
		    for c in contourf(dxxA,dyyA,xray2chisquared).collections:
			c.set_linewidth(0.1)
			c.set_alpha(1.0)
		    gca().set_xticks([-500.0,0.0,500.0])
		    xlabel('kpc')
		    ylabel('kpc')
		    colorbar()

		    subplot(2,3,4)
		    title("XRay-5000-8000eV - $\chi^2$ = %.3f"%xray3contrib,fontsize=10)
		    contourf(dxxA,dyyA,xray3chisquared)
		    for c in contourf(dxxA,dyyA,xray3chisquared).collections:
			c.set_linewidth(0.1)
			c.set_alpha(1.0)
		    gca().set_xticks([-500.0,0.0,500.0])
		    xlabel('kpc')
		    ylabel('kpc')
		    colorbar()

		    subplot(2,3,5)
		    title("Temp - $\chi^2$ = %.3f"%tempcontrib,fontsize=10)
		    contourf(dxxD,dyyD,tempchisquared)
		    for c in contourf(dxxD,dyyD,tempchisquared).collections:
			c.set_linewidth(0.1)
			c.set_alpha(1.0)
		    gca().set_xticks([-500.0,0.0,500.0])
		    xlabel('kpc')
		    ylabel('kpc')
		    colorbar()

		    subplot(2,3,6)
		    title("SZE - $\chi^2$ = %.3f"%szecontrib,fontsize=10)
		    contourf(dxxA,dyyA,szechisquared)
		    for c in contourf(dxxA,dyyA,szechisquared).collections:
			c.set_linewidth(0.1)
			c.set_alpha(1.0)
		    gca().set_xticks([-500.0,0.0,500.0])
		    xlabel('kpc')
		    ylabel('kpc')
		    colorbar()
		    pp.savefig()

	    plt.clf()
	    print "Finished plots for snapshot file", snap,"\n"
            sys.stdout.flush()
	except:
	    print "Unexpected error in plot routines.\n"
	    print"sys.exc_info=",sys.exc_info(),"\n"
	    pp.close()
	    return
	pp.close()
	return 

def SimpleFom(snap,Z,phi=0.0,theta=0.0,psi=0.0,ConstrainPhi=False,TFudge=1.0):

    # This is a new, faster Fom finding routine that uses only mass and X-ray, and searches for a 
    # minimum in time, phi, and theta.

    toppath = lageconfig.toppath
    (dataA,sigmaA,maskA,maskANull,dataB1,sigmaB1,dataB2,sigmaB2,dataB3,sigmaB3,dataC,sigmaC,dataD,sigmaD,maskD,maskDNull) = GetData()

    [dyyA,dxxA] = meshgrid(dataA.y,dataA.x)# Data grid for plots
    [dyyB,dxxB] = meshgrid(dataB1.y,dataB1.x)# Data grid for plots
    [dyyD,dxxD] = meshgrid(dataD.y,dataD.x)# Data grid for plots
    align = SetAlign(dataB1, phi, theta, psi, ConstrainPhi = ConstrainPhi)
    try:
            pf = GetPF(snap)

            dmsim=Array2d(2.0*dataA.xmin,2.0*dataA.xmax,2*dataA.nx,2.0*dataA.ymin,2.0*dataA.ymax,2*dataA.ny)
            masssim=Array2d(2.0*dataA.xmin,2.0*dataA.xmax,2*dataA.nx,2.0*dataA.ymin,2.0*dataA.ymax,2*dataA.ny)
            sumsim=Array2d(2.0*dataA.xmin,2.0*dataA.xmax,2*dataA.nx,2.0*dataA.ymin,2.0*dataA.ymax,2*dataA.ny)
            [syyA,sxxA] = meshgrid(dmsim.y,dmsim.x) # Sim grid for plots

            ApecData = ReadLookups(Z) # Reads the APEC lookup tables.
            for grid in pf.h.grids:
                    grid.set_field_parameter('ApecData',ApecData)
                    grid.set_field_parameter('TFudge',TFudge)

            [dmsim,masssim,xraysim1,xraysim2,xraysim3,szsim]=ProjectEnzoData(pf,masssim,phi,theta,psi)
            sumsim.data=gaussian_filter(dmsim.data+masssim.data,2.0)

            data1list=list((sumsim,xraysim1))
            data2list=list((dataA,dataB1))
            sigmalist=list((sigmaA,sigmaB1))
            masklist=list((maskANull,maskA))
            tol = 1.0E-7

            [shifteddata1list,align,fomx1]=FindBestShift(data1list, data2list, sigmalist, masklist, align, tol)

            data1list=list((sumsim,xraysim1))
            data2list=list((dataA,dataB1))
            sigmalist=list((sigmaA,sigmaB1))
            masklist=list((maskA,maskANull))
            tol = 1.0E-7

            [shifteddata1list,align,fommass]=FindBestShift(data1list, data2list, sigmalist, masklist, align, tol)

            return (fomx1,fommass)
    except:
        print "Error in SimpleFom routine",sys.exc_info()[0]
        return (1.0E5,1.0E5) # If there's an error, give it a large FOM

def MoviePlotsMHD1(snapmin,snapmax,xyarray,zmin,zmax,Z,phi=0.0,theta=0.0,psi=0.0,TFudge=1.0):
    # This creates a movie of a 2D slice through the snapshot file, using x, y limits,
    # from the input xyarray, with the slice having z extents from zmin to zmax.
    # This adds dark matter and plots gas temperature, xray intensity, dm density and an overlay.
    # This version uses a 'bright' color scheme.
    # this version uses log temperature and has a more 'zoomed-in' version.
    # This is the same as 5, but with only three stripes.

    nsnaps=snapmax-snapmin

    for snap in range(snapmin,snapmax): # Which snapshot files


	pf = GetPF(snap)
	simtime = pf.h.parameters['InitialTime'] * BulletConstants.TimeConversion # Put time in Gyears

	ApecData = ReadLookups(Z) # Reads the APEC lookup tables.
	for grid in pf.h.grids:
		grid.set_field_parameter('ApecData',ApecData)
                grid.set_field_parameter('TFudge',TFudge)


        start = time.time()
	gasrho=Array2d(xyarray.xmin,xyarray.xmax,xyarray.nx,xyarray.ymin,xyarray.ymax,xyarray.ny)
	DMrho=Array2d(xyarray.xmin,xyarray.xmax,xyarray.nx,xyarray.ymin,xyarray.ymax,xyarray.ny)
 
	simtime = pf.h.parameters['InitialTime'] * BulletConstants.TimeConversion # Put time in Gyears

	[DMrho,gasrho,xraysim1,xraysim2,xraysim3,szsim]=ProjectEnzoData(pf,gasrho,phi,theta,psi)


	[tempsim,bmagsim]=ProjectEnzoTemp(pf,gasrho,phi,theta,psi)
        PixelVolume=xyarray.dx * xyarray.dy * (zmax - zmin)
	DensityFactor = BulletConstants.cm_per_kpc**3 * PixelVolume / (BulletConstants.g_per_Msun * 1E10)
	gasrho.data = gasrho.data / DensityFactor
	#for i in range(100):
		#DMrho.data[:,i] = DMrho2X.data[:,150-i]
	DMrho.data = gaussian_filter(DMrho.data,0.5) / DensityFactor # Smooth and normalize 
	#print "gasrhomin = %.4g, gasrhomax = %.4g, dmrhomin = %.4g, dmrhomax = %.4g, bmin = %.4g, bmax = %.4g\n"%(gasrho.data.min(), gasrho.data.max(), DMrho.data.min(), DMrho.data.max(), bmagsim.data.min(), bmagsim.data.max())
	#sys.exit()
	LogRhoMin = -27.5
	LogRhoMax = -24.5
	LogTMin = 0.0
	LogTMax = 1.5
	LogBMin = -7.5
	LogBMax = -5.5
	
	for i in range(xyarray.nx):
		for j in range(xyarray.ny):
			gasrho.data[i,j]=max(gasrho.data[i,j],1.01*10**LogRhoMin)
			DMrho.data[i,j]=max(DMrho.data[i,j],1.01*10**LogRhoMin)
			tempsim.data[i,j]=max(tempsim.data[i,j],1.01*10**LogTMin)
			bmagsim.data[i,j]=max(bmagsim.data[i,j],1.01*10**LogBMin)
	
        [dyy,dxx] = meshgrid(xyarray.y,xyarray.x)# Data grid for plots
 
        figure(figsize=(8,9))
        suptitle('Bullet Cluster Collision, T = %.2f Gy'%simtime, fontsize=16)
        subplots_adjust(hspace=0.3, wspace=0.2)

	subplot("311",aspect='equal')
	title("Log Density ($g/cm^3$) DM - Contours", fontsize=12)
	Gaslevels=linspace(LogRhoMin,LogRhoMax,200) 
	for i in range(200): Gaslevels[i]=int(Gaslevels[i]*200)/200.0
	contourf(dxx,dyy,log10(gasrho.data),Gaslevels,cmap=cm.spectral)
	for c in contourf(dxx,dyy,log10(gasrho.data),Gaslevels,cmap=cm.spectral).collections:
		c.set_linewidth(0.1)
		c.set_alpha(1.0)
	gca().set_xticks([-1000.0,0.0,1000.0])
	colorbar()
	DMRholevels=linspace(-26.0,-23.0,12)  
        contour(dxx,dyy,log10(DMrho.data),DMRholevels,linestyle='-',colors='w')
	ylabel('kpc')

        subplot("312",aspect='equal')
        title('Log Gas Temperature(keV)', fontsize=12)
	Templevels=linspace(LogTMin,LogTMax,200) 
	for i in range(200): Templevels[i]=int(Templevels[i]*200)/200.0
        contourf(dxx,dyy,log10(tempsim.data),Templevels,cmap=cm.spectral)
	for c in contourf(dxx,dyy,log10(tempsim.data),Templevels,cmap=cm.spectral).collections:
    		c.set_linewidth(0.1)
    		c.set_alpha(1.0)
        colorbar()
	gca().set_xticks([-1000.0,0.0,1000.0])
	ylabel('kpc')

        subplot("313",aspect='equal')
        title('Log Magnetic Field (Gauss)', fontsize=12)
	Blevels=linspace(LogBMin,LogBMax,200) 
	for i in range(200): Blevels[i]=int(Blevels[i]*200)/200.0
        contourf(dxx,dyy,log10(bmagsim.data),Blevels,cmap=cm.spectral)
	for c in contourf(dxx,dyy,log10(bmagsim.data),Blevels,cmap=cm.spectral).collections:
    		c.set_linewidth(0.1)
    		c.set_alpha(1.0)
        colorbar()
	gca().set_xticks([-1000.0,0.0,1000.0])
	xlabel('kpc')
	ylabel('kpc')

        filename = 'movie/collision'+str('%03d'%snap)+'.png'
	plt.savefig(filename)#, size=(2400,1800))
        print 'Wrote file', filename
        plt.clf()

    
    return 


