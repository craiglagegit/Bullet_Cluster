EXEC   = Div_Clean.so

OBJS   =  divclean.o kvector.o functions.o fftOperators.o div_clean.o

%.o: %.cc
	$(CC) $(CFLAGS) -c $< -o $@

#GSL_INCL =  -I/export/ursa1/csl336/Software/gsl-1.9/include
#GSL_LIBS =  -L/export/ursa1/csl336/Software/gsl-1.9/lib
#FFTW_INCL=  -I/export/ursa1/csl336/Software/fftw-3.3/include
#FFTW_LIBS=  -L/export/ursa1/csl336/Software/fftw-3.3/lib

GSL_INCL =  -I/home/install/app/gsl-1.9/include
GSL_LIBS =  -L/home/install/app/gsl-1.9/lib
FFTW_INCL=  -I/home/install/app/fftw-3.2.2/include
FFTW_LIBS=  -L/home/install/app/fftw-3.2.2/lib

CFLAGS =  -O3 $(OPTIONS) $(GSL_INCL) $(FFTW_INCL)   -g -fPIC -Wall  -lstdc++
LNKCMD =  g++ -shared 

LFLAGS =  $(GSL_LIBS) $(FFTW_LIBS) 

LIBS =  -lfftw3 -lgsl -lgslcblas -lm -g

CC     =  g++ -g

$(EXEC): $(OBJS)
	 $(LNKCMD) $(OBJS) $(LIBS) $(LFLAGS) -o $(EXEC)

$(OBJS): $(INCL) 

clean:
	 rm -f $(OBJS) $(EXEC)

