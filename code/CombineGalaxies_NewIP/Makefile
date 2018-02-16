
#OPT    +=  -DORBITCORRECTION   # this will correct the orbit if the halos already overlap
                                # the dark halos are assumed to follow a Hernquist profile in this

#OPT    +=  -DVELINCREASE=1.055

EXEC   = CombineGalaxies

OBJS   = main.o  globvars.o \
	 load.o save.o \
	 move.o turn.o \
	 nrsrc/nrutil.o



INCL   = globvars.h nrsrc/nrutil.h Makefile


#SYSTYPE="MPA"
#SYSTYPE="Mako"
#SYSTYPE="Regatta"
#SYSTYPE="RZG_LinuxCluster"
SYSTYPE="RZG_LinuxCluster-gcc"
#SYSTYPE="OPA-Cluster64"

CC       =    gcc        # sets the C-compiler (default)
OPTIMIZE =   -O   # optimization and warning flags (default)


ifeq ($(SYSTYPE),"MAKO")
GSL_INCL = #-I/home/install/app/gsl/include
GSL_LIBS = #-L/home/install/app/gsl/lib
endif


ifeq ($(SYSTYPE),"MPA")
GSL_INCL = -I/usr/common/pdsoft/include
GSL_LIBS = -L/usr/common/pdsoft/lib
endif


ifeq ($(SYSTYPE),"OPA-Cluster64")
CC       =   mpiccg -m64
OPTIMIZE =  -O2 -Wall
GSL_INCL =  -I/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-0.9.3/include
GSL_LIBS =  -L/afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-0.9.3/lib  -Wl,"-R /afs/rzg/home/v/vrs/Libs/opteron64/mvapich2-0.9.3/lib"
endif



ifeq ($(SYSTYPE),"Regatta")
CC       =   cc
OPTIMIZE =   -g -q64 #  -O5 -qstrict -qipa 
GSL_INCL = -I/afs/ipp-garching.mpg.de/u/vrs/gsl_psi64/include  -I/afs/rzg/@sys/lib/include
GSL_LIBS = -L/afs/ipp-garching.mpg.de/u/vrs/gsl_psi64/lib                
endif

ifeq ($(SYSTYPE),"RZG_LinuxCluster")
CC       =   icc
OPTIMIZE =   -O3 
GSL_INCL = -I/home/install/app/gsl/include
GSL_LIBS = -L/home/install/app/gsl/lib                -static
endif

ifeq ($(SYSTYPE),"RZG_LinuxCluster-gcc")
CC       =   gcc
OPTIMIZE =   -O3 
GSL_INCL = -I/home/install/app/gsl/include
GSL_LIBS = -L/home/install/app/gsl/lib               
endif



LIBS   =    -lm  $(GSL_LIBS) -lgsl -lgslcblas 

CFLAGS =    $(OPTIONS) $(OPT)  $(OPTIMIZE)  $(GSL_INCL)


$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS)  -o $(EXEC)  

$(OBJS): $(INCL) 

clean:
	 rm -f $(OBJS) $(EXEC)


