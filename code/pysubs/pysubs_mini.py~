#!/usr/bin/env python
#Author: Craig Lage, NYU; 
#Date: 14-Sep-12

'''
These subroutines perform various reading, writing, and plotting functions for Gadget 3 and Enzo files for simulating the Bullet cluster.

'''
import lageconfig # system specific path information
from pylab import *
from subprocess import *
import time
import array
from scipy.ndimage import gaussian_filter, convolve
from scipy.ndimage.filters import sobel
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d.axes3d import Axes3D
from scipy.special import exp1,gamma # Exponential integral
import h5py
import sys
import pyfits
import pickle

import BulletConstants # Constants used in Bullet cluster simulations

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
        self.d=zeros([5])
        self.dmin=zeros([5])
        self.dmax=zeros([5])

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
    Synch=Array2d(mass.xmin,mass.xmax,mass.nx,mass.ymin,mass.ymax,mass.ny)

    add_field('TXRay', function = _EnzoTXRay, units=r"\rm{cm}^{-3}\rm{s}^{-1}", projected_units=r"\rm{cm}^{-2}\rm{s}^{-1}", validators=[ValidateParameter('ApecData')], take_log=False)
    add_field('XRay', function = _EnzoXRay, units=r"\rm{cm}^{-3}\rm{s}^{-1}", projected_units=r"\rm{cm}^{-2}\rm{s}^{-1}", validators=[ValidateParameter('ApecData')], take_log=False)
    add_field('BMag', function = _EnzoBMag,take_log=False)
    add_field('Synch', function = _EnzoSynch, units=r"\rm{cm}^{-3}\rm{s}^{-1}", projected_units=r"\rm{cm}^{-2}\rm{s}^{-1}",take_log=False)
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
    BFactor    = 1.0 / (BulletConstants.cm_per_kpc * width[2]) 
    SynchFactor = XFactor * BulletConstants.microJanskys_per_CGS
    projcam=ProjectionCamera(center,normal_vector,width,resolution,"TXRay",north_vector=north_vector,pf=pf, interpolated=True)
    Temp.data = projcam.snapshot()[:,:,0]
    projcam=ProjectionCamera(center,normal_vector,width,resolution,"XRay",north_vector=north_vector,pf=pf, interpolated=True)
    Temp.data = Temp.data / projcam.snapshot()[:,:,0]
    projcam=ProjectionCamera(center,normal_vector,width,resolution,"BMag",north_vector=north_vector,pf=pf, interpolated=True)
    BMag.data = projcam.snapshot()[:,:,0] * BFactor
    projcam=ProjectionCamera(center,normal_vector,width,resolution,"Synch",north_vector=north_vector,pf=pf, interpolated=True)
    Synch.data   = projcam.snapshot()[:,:,0] * SynchFactor
    elapsed = (time.time()-start)
    print "Elapsed time to run Enzo temp projection = "+str(elapsed)
    return [Temp,BMag,Synch]


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

def _EnzoSynch(field,data):
    P = data.get_field_parameter('SpectralIndex')
    Eta = 1.00 # Fudge Factor
    PFactor = (P-2.0)/(P+1.0)*gamma(P/4.0+19.0/12.0)*gamma(P/4.0-1.0/12.0)*gamma((P+5.0)/4.0)/gamma((P+7.0)/4.0)
    B = sqrt(data['Bx']**2 + data['By']**2 + data['Bz']**2)
    omegac = BulletConstants.OmegaCPreFactor * B
    omega =  2.0 * pi * BulletConstants.RadioFrequency * (1 + BulletConstants.redshift)
    return  Eta * BulletConstants.RadioPreFactor * PFactor * pow(B,3.0) * pow(omega/omegac,-(P-1.0)/2.0)

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

def FindEnzoCentroidsObserved(pf, phi=0.0, theta = 0.0, psi = 0.0): 

    enzparam = ReadEnzoParameterFile("AMRTest.enzo")
    EnzoLength     = enzparam.LengthUnits
    EnzoDensity    = enzparam.DensityUnits
    EnzoMass       = EnzoDensity * EnzoLength**3.0
    R250kpc = 250.0 / EnzoLength * BulletConstants.cm_per_kpc
    print "R250kpc = %.4g"%R250kpc 
    start = time.time()
    #print 'Starting Centroid finding\n'
    R=EulerAngles(phi,theta,psi)
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
                X = [pf.h.grids[i]['particle_position_x'][j],pf.h.grids[i]['particle_position_y'][j],pf.h.grids[i]['particle_position_z'][j]]
                XP = dot(R,X)
     		x.append(XP[0])
     		y.append(XP[1])
   		z.append(XP[2])
		NPart = NPart + 1

    Masses=zeros([2])
    NumPart=zeros([2])
    Centroids=zeros([2,3])
    MassWithin250=zeros([2])
    NumWithin250=zeros([2])
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

    for k in range(2): # k denotes the bullet or main particles
        for m in range(3):
            Centroids[k,m] = Centroids[k,m] / NumPart[k]

    for j in range(NPart): # n cycles through the number of particles of this type
        if abs(mass[j] - Masses[0]) < 1E-8:	
            k = 0
        else:
            k = 1
        radius = sqrt((Centroids[k,0]-x[j])**2 + (Centroids[k,1]-y[j])**2)
        if radius < R250kpc:
            MassWithin250[k] = MassWithin250[k] + mass[j]
            NumWithin250[k] = NumWithin250[k] + 1

    for k in range(2):
        MassWithin250[k] = MassWithin250[k] * EnzoMass / BulletConstants.g_per_Msun
        print "For cluster %d, N within 250 = %d"%(k,int(NumWithin250[k]))
    elapsed = (time.time()-start)
    print "Elapsed time to locate centroids = "+str(elapsed)

    return [NumPart, Masses, Centroids, MassWithin250]

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
            if m < 0 or m > data.nx -1:
                continue
            deltax=abs((xprime-data.x[m])/data.dx)
            for n in [j-1,j,j+1]:
                if n < 0 or n > data.ny -1:
                    continue
                deltay=abs((yprime-data.y[n])/data.dy)
                d=d+PyramidalKernel(deltax,deltay)*data.data[m,n]
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


def GetKernels():
        # This reads in kernels used for image manipulation

	toppath = lageconfig.toppath
	datapathszekernel=toppath+'bullet/data/sze/bullet_transfer_rescaled_small.fits'
	datapathkernel60=toppath+'bullet/data/kernel60.dat'
	datapathkernel15=toppath+'bullet/data/kernel15.dat'
	datapathkernel20=toppath+'bullet/data/kernel20.dat'
	szekernel=Array2d(0.0,400.0,113,0.0,400.0,113)
	szekernel=GetBulletFits(datapathszekernel,szekernel)

	kernel60=Array2d(0.0,1.0,110,0.0,1.0,110)
	kernel60=GetBulletData(datapathkernel60,kernel60)
	kernel15=Array2d(0.0,1.0,110,0.0,1.0,110)
	kernel15=GetBulletData(datapathkernel15,kernel15)
	kernel20=Array2d(0.0,1.0,110,0.0,1.0,110)
	kernel20=GetBulletData(datapathkernel20,kernel20)
        return (szekernel, kernel60, kernel15, kernel20)

def GetData():
        # This reads in the image data
	toppath = lageconfig.toppath
	datapathA=toppath+'bullet/data/kappa_25Apr12.dat'
	datapathB1=toppath+'bullet/data/xray_500_2000_29Jun11.dat'
	datapathB2=toppath+'bullet/data/xray_2000_5000_29Jun11.dat'
	datapathB3=toppath+'bullet/data/xray_5000_8000_29Jun11.dat'
	datapathC=toppath+'bullet/data/sze_data_13Mar13.dat'
	datapathD=toppath+'bullet/data/xray_temp_23Jun11.dat'
        datapathE=toppath+'bullet/data/radio_24Apr13.dat'

	datapathsA=toppath+'bullet/data/mass_sigma_11Feb13.dat'
	datapathsB1=toppath+'bullet/data/xray_sigma_500_2000_12Mar13.dat'
	datapathsB2=toppath+'bullet/data/xray_sigma_2000_5000_12Mar13.dat'
	datapathsB3=toppath+'bullet/data/xray_sigma_5000_8000_12Mar13.dat'
	datapathsC=toppath+'bullet/data/sze_sigma_13Mar13.dat'
	datapathsD=toppath+'bullet/data/xray_temp_sigma_23Jun11.dat'
        datapathsE=toppath+'bullet/data/radio_sigma_5Oct11.dat'

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

        dataE=Array2d(-dxmaxA,dxmaxA,xpixelsA,-dymaxA,dymaxA,ypixelsA)
        sigmaE=Array2d(-dxmaxA,dxmaxA,xpixelsA,-dymaxA,dymaxA,ypixelsA)
        dataE=GetBulletData(datapathE,dataE)# Radio measured data                                                                                     
        sigmaE=GetBulletData(datapathsE,sigmaE)# Radio Sigma                                                                       
        return (dataA,sigmaA,maskA,maskANull,dataB1,sigmaB1,dataB2,sigmaB2,dataB3,sigmaB3,dataC,sigmaC,dataD,sigmaD,maskD,maskDNull,dataE,sigmaE)


