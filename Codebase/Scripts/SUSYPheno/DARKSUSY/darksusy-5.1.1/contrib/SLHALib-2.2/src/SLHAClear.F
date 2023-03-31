* SLHAClear.F
* clear (i.e. set to some bogus value) all data in the
* SUSY Les Houches Accord data structure (see hep-ph/0311123)
* last modified 21 Apr 06 th

#include "SLHADefs.h"


	subroutine SLHAClear(slhadata)
	implicit none
	double complex slhadata(*)

	integer i

	do i = 1, OffsetDecays
	  slhadata(i) = invalid
	enddo

	do i = OffsetDecays + 1, nslhadata
	  slhadata(i) = DCMPLX(invalid, invalid)
	enddo
	end

