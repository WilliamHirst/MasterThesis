* demo.F
* this program demonstrates use of the SLHA library
* last modified 21 Apr 06 th


	program SLHADemo
	implicit none

#include "SLHA.h"

#define Strip(s) s(1:lnblnk(s))

	integer error
	double complex slhadata(nslhadata)

	character*80 argv

	integer iargc, lnblnk
	external iargc

	if( iargc() .ne. 1 ) then
	  call getarg(0, argv)
	  print *, "USAGE:"
	  print *, Strip(argv), " FILE"
	  print *, "reads FILE and writes it out again as FILE.out"
	  stop
	endif

	call getarg(1, argv)

	call SLHARead(error, slhadata, argv, 1)
	if( error .ne. 0 ) stop "Read error"

	call SLHAInfo(BlockSPInfo(1), 1, "SLHADemo")
	call SLHAInfo(BlockSPInfo(1), 2, "1.0")

	call SLHAWrite(error, slhadata, Strip(argv)//".out")
	if( error .ne. 0 ) stop "Write error"
	end

