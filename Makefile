prefix ?= /usr/local/bin #This can be changed
CC ?= gcc
AR ?= ar
RANLIB ?= ranlib
OPTS ?= -Wall -g -O3

.PHONY: all clean htslib install clean-all version.h

.SUFFIXES:.c .o

all: lib PileOMeth

OBJS = common.o bed.o svg.o pileup.o extract.o MBias.o mergeContext.o
VERSION = 0.1.10

#If we're building from a git repo, then append the most recent tag
ifneq "$(wildcard .git)" ""
VERSION := $(shell git describe --tags --always --dirty)
endif

version.h:
	echo '#define VERSION "$(VERSION)"' > $@

.c.o:
	$(CC) -c $(OPTS) -Ihtslib $< -o $@

htslib: 
	$(MAKE) -C htslib

libPileOMeth.a: version.h $(OBJS)
	-@rm -f $@
	$(AR) -rcs $@ $(OBJS)

lib: libPileOMeth.a

PileOMeth: htslib version.h libPileOMeth.a
	$(CC) $(OPTS) -Ihtslib -o PileOMeth main.c libPileOMeth.a htslib/libhts.a -lm -lz -lpthread

test: PileOMeth
	python tests/test.py

clean:
	rm -f *.o PileOMeth libPileOMeth.a

clean-all: clean
	make --directory=htslib clean

install: PileOMeth
	install PileOMeth $(prefix)
