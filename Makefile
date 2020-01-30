prefix ?= /usr/local/bin #This can be changed
CC ?= gcc
LIBS ?=  # e.g., -L$PREFIX/lib, or where ever htslib is
LIBBIGWIG ?=
CFLAGS ?= -Wall -g -O3 -pthread

.PHONY: all clean install version.h

.SUFFIXES:.c .o

all: MethylDackel

OBJS = common.o bed.o svg.o pileup.o extract.o MBias.o mergeContext.o perRead.o
VERSION = 0.4.0

version.h:
	echo '#define VERSION "$(VERSION)"' > $@

.c.o:
	$(CC) -c $(CFLAGS) $(LIBS) -IlibBigWig $< -o $@

#libbw: 
#	$(MAKE) -C libBigWig all

libMethylDackel.a: version.h $(OBJS)
	-@rm -f $@
	$(AR) -rcs $@ $(OBJS)

lib: libMethylDackel.a

#MethylDackel: libbw version.h libMethylDackel.a $(OBJS)
MethylDackel: libMethylDackel.a version.h $(OBJS)
	$(CC) $(CFLAGS) $(LIBS) -o MethylDackel $(OBJS) main.c libMethylDackel.a $(LIBBIGWIG) -lm -lz -lpthread -lhts -lcurl

test: MethylDackel 
	python tests/test.py

clean:
	rm -f *.o MethylDackel libMethylDackel.a

clean-all: clean
	make --directory=libBigWig clean
	rm -f *.o MethylDackel

install: MethylDackel
	install MethylDackel $(prefix)
