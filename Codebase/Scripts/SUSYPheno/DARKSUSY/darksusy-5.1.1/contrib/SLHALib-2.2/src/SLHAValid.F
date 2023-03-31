* SLHAValid.F
* determine whether a block has all valid entries
* last modified 19 Dec 06 th

#include "SLHADefs.h"


	integer function SLHAValid(slhablock, length)
	implicit none
	double complex slhablock(*)
	integer length

	integer i

	SLHAValid = 0
	do i = 1, length
	  if( DBLE(slhablock(i)) .eq. invalid ) return
	enddo
	SLHAValid = 1
	end

