PREFIX = /usr/local/bin #This can be changed"
CC = gcc
OPTS = -Wall -g -O3

OBJS = alignmentHeap.o bloomFilter.o graph.o murmur3.o TargetCreator.o realigner.o SemiGlobal.o threads.o

.PHONY: all clean htslib install clean-all

.SUFFIXES:.c .o

all: PileOMeth

.c.o:
	$(CC) -c $(OPTS) -I$(INCLUDE_DIRS) $< -o $@

htslib: 
	$(MAKE) -C htslib

PileOMeth: htslib
	$(CC) $(OPTS) -Ihtslib -o PileOMeth bed.c PileOMeth.c htslib/libhts.a -lz -lpthread

clean:
	rm -f *.o PileOMeth

clean-all: clean
	make --directory=htslib clean

install: PileOMeth
	install PileOMeth $(PREFIX)
