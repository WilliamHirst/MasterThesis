* SLHACopyString.F
* utility routine for SLHAMessage
* last modified 27 Nov 06 th

#include "SLHADefs.h"


	subroutine SLHACopyString(n, s, d)
	implicit none
	integer n
	character*1 s(n), d(SPInfo_Len)

	integer m, i

	m = min(n, SPInfo_Len)
	do i = 1, m
	  d(i) = s(i)
	enddo
	do i = m + 1, SPInfo_Len
	  d(i) = ' '
	enddo
	end

