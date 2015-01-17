prefix= /usr/local/bin #This can be changed"
CC = gcc
OPTS = -Wall -g -O3

OBJS = pileup.o

.PHONY: all clean htslib install clean-all

.SUFFIXES:.c .o

all: PileOMeth

.c.o:
	$(CC) -c $(OPTS) -Ihtslib $< -o $@

htslib: 
	$(MAKE) -C htslib

PileOMeth: htslib $(OBJS)
	$(CC) $(OPTS) -Ihtslib -o PileOMeth bed.c PileOMeth.c htslib/libhts.a pileup.o -lz -lpthread

clean:
	rm -f *.o PileOMeth

clean-all: clean
	make --directory=htslib clean

install: PileOMeth
	install PileOMeth $(prefix)
