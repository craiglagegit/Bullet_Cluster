EXEC   = HsmlAndGrid.so

OBJS   = main.o tree.o peano.o allvars.o

#OPTIONS += -DPERIODIC


INCL   = proto.h allvars.h

CFLAGS =  -O3 $(OPTIONS)  -g -fPIC -Wall
LNKCMD =  ld -shared

LIBS   =  -lm -g

CC     =  gcc

$(EXEC): $(OBJS)
	 $(LNKCMD)  $(OBJS) $(LIBS)  -o $(EXEC)


$(OBJS): $(INCL)

clean:
	 rm -f $(OBJS) $(EXEC)




