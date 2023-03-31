* table.F
* this program prints out a decay table
* last modified 25 Apr 06 th


#include "PDG.h"

	program SLHADemo
	implicit none

#include "SLHA.h"

#define Strip(s) s(1:lnblnk(s))

	integer error
	double complex slhadata(nslhadata)
	double precision width(1024)
	integer id(0:7,1024)
	integer nparticles, nchannels, particle, channel
	character*(PDGLen) name(8)
	character*(8*PDGLen) decay_name

	character*80 argv

	integer iargc, lnblnk
	external iargc

	if( iargc() .ne. 1 ) then
	  call getarg(0, argv)
	  print *, "USAGE:"
	  print *, Strip(argv), " FILE"
	  print *, "reads FILE and displays its decay table"
	  stop
	endif

	call getarg(1, argv)

	call SLHARead(error, slhadata, argv, 1)
	if( error .ne. 0 ) stop "Read error"

	nchannels = SLHADecayTable(slhadata, 0, width, id, 7, 1024)
	do channel = 1, nchannels
	  nparticles = id(0,channel)

	  do particle = 1, nparticles
	    call SLHAPDGName(id(particle,channel), name(particle))
	  enddo

	  if( nparticles .eq. 1 ) then
100	    format("Gamma(", A, ")", T48, F13.5, " GeV")
	    print 100, Strip(name(1)), width(channel)
	  else
101	    format(A, " ->", 7 (:," ",A))
102	    format("  BR(", A, ")", T57, G10.5)
	    write(decay_name, 101) (Strip(name(particle)),
     &        particle = 1, nparticles)
	    print 102, Strip(decay_name), width(channel)
	  endif
	enddo

	end

