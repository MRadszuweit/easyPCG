CFLAGS = -Wall
CC = gcc
OPTIONS = $(OPT)

ifeq ($(strip $(OPT)),debug)

	OPTIONS = -g -pg	

endif

ifeq ($(strip $(OPT)),)
		
	OPTIONS = -march=native	

endif

ifeq ($(strip $(OPT)),O3)
		
	OPTIONS = -ffast-math -O3

endif

LDFLAGS = -lm 
CFILES = testPCG.c easyPCG.c linalg_stuff.c
OFILES = $(CFILES:%.c=%.o)

compile: $(OFILES)

	$(CC) $(CFLAGS) $(OPTIONS) -o testPCG $(OFILES) $(LDFLAGS)

%.o : %.c
	$(CC) $(CFLAGS) $(OPTIONS) -c $< $(LDFLAGS)

clean:
	rm -f *.o testPCG

