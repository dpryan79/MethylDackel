prefix= /usr/local/bin #This can be changed"
CC = gcc
OPTS = -Wall -g -O3

OBJS = pileup.o bed.o

.PHONY: all clean htslib install clean-all

.SUFFIXES:.c .o

OBJS = bed.o svg.o pileup.o

all: PileOMeth PileOMethMBias

.c.o:
	$(CC) -c $(OPTS) -Ihtslib $< -o $@

htslib: 
	$(MAKE) -C htslib

PileOMeth: htslib $(OBJS)
	$(CC) $(OPTS) -Ihtslib -o PileOMeth PileOMeth.c bed.o pileup.o htslib/libhts.a -lz -lpthread

PileOMethMBias: htslib $(OBJS)
	$(CC) $(OPTS) -Ihtslib -o PileOMethMBias MBias.c svg.o bed.o htslib/libhts.a -lm -lz -lpthread

clean:
	rm -f *.o PileOMeth PileOMethMBias

clean-all: clean
	make --directory=htslib clean

install: PileOMeth
	install PileOMeth PileOMethMBias $(prefix)
