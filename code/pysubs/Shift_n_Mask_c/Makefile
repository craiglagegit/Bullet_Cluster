EXEC   = Shift_n_Mask.so

OBJS   = shift_n_mask.o

#OPTIONS += -DPERIODIC


INCL   = 

CFLAGS =  -O3 $(OPTIONS)  -g -fPIC -Wall -fno-exceptions -lstdc++ 
LNKCMD =  g++ -shared 

LIBS   =  -lm -g

CC     =  g++

$(EXEC): $(OBJS)
	 $(LNKCMD)  $(OBJS) $(LIBS)  -o $(EXEC)


$(OBJS): $(INCL)

clean:
	 rm -f $(OBJS) $(EXEC)




