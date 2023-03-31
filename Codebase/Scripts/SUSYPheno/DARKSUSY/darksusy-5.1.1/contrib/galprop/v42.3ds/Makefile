#CXX = g++
#FC = g77
CFITSIO = ../../cfitsio
CPPFLAGS = -O3 -Wno-deprecated -I${CFITSIO}/include
CPPFLAGS = -g -O -Wno-deprecated -I${CFITSIO}/include
FFLAGS =
LIBS = -lm -lg2c
FITSLIB = -L$(CFITSIO)/lib -lcfitsio -Wl,-rpath,$(CFITSIO)/lib
LDFLAGS = $(FITSLIB) $(LIBS)                     
FOBJS := $(patsubst %.f,%.o,$(wildcard *.f))                 
CCOBJS := $(patsubst %.cc,%.o,$(wildcard *.cc))
galprop:  ${FOBJS} ${CCOBJS}
#	$(CXX) *.o -o $@ ${LDFLAGS}
	rm -f main.o
	ar r ../../../lib/libdarksusy.a *.o

clean:
	rm -f *.o galprop

.PHONY: galprop
