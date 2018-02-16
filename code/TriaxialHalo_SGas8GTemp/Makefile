
EXEC   = TriaxialHalo

OBJS   = triaxial.o utils.o fileio.o models.o 


INCL   =  Makefile


CFLAGS = -g 

LIBS   =  -L/usr/lib -lm -lgsl -lgslcblas

CC     =  gcc 


$(EXEC): $(OBJS) 
	$(CC)  $(CFLAGS) $(OBJS) -o $(EPATH)$(EXEC)  $(LIBS)

$(OBJS): $(INCL)

clean:
	rm -f $(OBJS) $(EXEC)




