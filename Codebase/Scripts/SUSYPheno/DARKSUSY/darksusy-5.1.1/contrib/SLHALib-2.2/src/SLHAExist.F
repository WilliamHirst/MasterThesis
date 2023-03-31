* SLHAExist.F
* determine whether a block has any valid entries
* last modified 21 Apr 06 th

#include "SLHADefs.h"


	integer function SLHAExist(slhablock, length)
	implicit none
	double complex slhablock(*)
	integer length

	integer i

	SLHAExist = 0
	do i = 1, length
	  if( DBLE(slhablock(i)) .ne. invalid ) then
	    if( DIMAG(slhablock(i)) .ne. 0 ) then
	      SLHAExist = 2
	      return
	    endif
	    SLHAExist = 1
	  endif
	enddo
	end

