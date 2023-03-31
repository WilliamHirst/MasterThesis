OBJS = hdecay.o haber.o feynhiggs.o hsqsq.o susylha.o hgaga.o dmb.o elw.o hgg.o

OBJS0 = hdecay0.o haber.o feynhiggs.o hsqsq.o susylha.o hgaga.o dmb.o elw.o hgg.o

#OBJS = hdecay.o haber.o hsqsq.o

#FFLAGS = -std=gnu

FC=gfortran

#FFLAGS = -fno-emulate-complex -fno-automatic -ffixed-line-length-none -ffast-math -march=pentiumpro -Wall -fno-silent

#FC=g77

#FFLAGS = -fno-emulate-complex -fno-automatic -ffixed-line-length-none -ffast-math -march=pentiumpro -malign-double -Wall -fno-silent

#FFLAGS = -Wall -fno-silent

#FC=f77

#FFLAGS= -pc 64 -g77libs

#FC=pgf77

.f.o:
	$(FC) -c $(FFLAGS) $*.f

hdecay: $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o run

hdecay0: $(OBJS0)
	$(FC) $(FFLAGS) $(OBJS0) -o run0

clean:
	rm -f $(OBJS)
