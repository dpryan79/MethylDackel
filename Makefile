PREFIX = /usr/local/bin #This can be changed!
CC = gcc
OPTS = -Wall -g -O3

.PHONY: all clean htslib install clean-all

.SUFFIXES:.c .o

all: PileOMeth MBias

.c.o:
	$(CC) -c $(OPTS) -Ihtslib $< -o $@

htslib: 
	$(MAKE) -C htslib

PileOMeth: htslib bed.o
	$(CC) $(OPTS) -Ihtslib -o PileOMeth PileOMeth.c bed.o htslib/libhts.a -lz -lpthread

MBias: htslib svg.o
	$(CC) $(OPTS) -Ihtslib -o MBias MBias.c svg.o htslib/libhts.a -lm -lz -lpthread

clean:
	rm -f *.o PileOMeth MBias

clean-all: clean
	make --directory=htslib clean

install: PileOMeth
	install PileOMeth MBias $(PREFIX)
