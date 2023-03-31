OBJS = suspect2.o twoloophiggs.o bsg.f sdecay.o

FC=g77

.f.o: 
	$(FC) -c -finit-local-zero $*.f

sdecay: $(OBJS)
	$(FC) $(OBJS) -o run
