prefix ?= /usr/local/bin #This can be changed
CC ?= gcc
AR ?= ar
RANLIB ?= ranlib
OPTS ?= -Wall -g -O3 -pthread

.PHONY: all clean htslib install clean-all version.h

.SUFFIXES:.c .o

all: lib MethylDackel

OBJS = common.o bed.o svg.o pileup.o extract.o MBias.o mergeContext.o
VERSION = 0.3.0

#If we're building from a git repo, then append the most recent tag
ifneq "$(wildcard .git)" ""
VERSION := $(shell git describe --tags --always --dirty)
endif

version.h:
	echo '#define VERSION "$(VERSION)"' > $@

.c.o:
	$(CC) -c $(OPTS) -Ihtslib -IlibBigWig $< -o $@

htslib: 
	$(MAKE) -C htslib

libbw: 
	$(MAKE) -C libBigWig all

libMethylDackel.a: version.h $(OBJS)
	-@rm -f $@
	$(AR) -rcs $@ $(OBJS)

lib: libMethylDackel.a

MethylDackel: htslib libbw version.h libMethylDackel.a
	$(CC) $(OPTS) -Ihtslib -IlibBigWig -o MethylDackel main.c libMethylDackel.a htslib/libhts.a libBigWig/libBigWig.so -lm -lz -lpthread

test: MethylDackel 
	python tests/test.py

clean:
	rm -f *.o MethylDackel libMethylDackel.a

clean-all: clean
	make --directory=htslib clean
	make --directory=libBigWig clean

install: MethylDackel
	install MethylDackel $(prefix)
