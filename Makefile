# makefile for phynder, developed on Richard's Mac.

#CFLAGS= -O3
CFLAGS= -g				# for debugging
HTSDIR=../htslib
CPPFLAGS=-I$(HTSDIR)
HTSLIB=$(HTSDIR)/libhts.a
LDLIBS=$(HTSLIB) -lpthread -lm -lcurl -lz -lbz2 -llzma 

all: phynder

clean:
	$(RM) -r phynder *.o *~ *test *.dSYM

install: phynder
	cp phynder ~/bin

### header dependencies

UTILS_HEADERS=utils.h array.h hash.h dict.h

tree.h: $(UTILS_HEADERS)

newick.h: tree.h

vcf.h: $(UTILS_HEADERS)

### object files

UTILS_OBJS=dict.o array.o hash.o utils.o
$(UTILS_OBJS): utils.h $(UTILS_HEADERS)

tree.o: tree.h

newick.o: newick.h

vcf.o: vcf.h

### programs

phynder: phynder.c tree.o newick.o vcf.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

newick-test: newick.c tree.o $(UTILS_OBJS)
	$(CC) -DTEST $(CFLAGS) $^ -o $@ $(LDLIBS)

### end of file
