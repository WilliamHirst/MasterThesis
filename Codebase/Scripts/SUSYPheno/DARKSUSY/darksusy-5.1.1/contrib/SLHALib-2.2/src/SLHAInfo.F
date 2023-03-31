* SLHAInfo.F
* record informational messages for write-out to block SP/DCINFO
* last modified 11 Oct 06 th

#include "SLHADefs.h"

* the following two definitions must match the ones in SLHA2.m
#define MAXMSG 15
#define MSGLEN 5

* codes are, according to the Accord,
*   1 = program name
*   2 = program version
*   3 = warning messages
*   4 = error messages


	subroutine SLHAInfo(slhablock, code, text)
	implicit none
	double complex slhablock(*)
	integer code
	character*(*) text

	integer n

	n = int(slhablock(2))
	if( n .eq. invalid ) then
	  n = 0
	else if( n .ge. MAXMSG ) then
	  print *, "No space left for messages"
	  return
	endif

	slhablock(2) = n + 1
	slhablock(n + 3) = code
	slhablock(1) = max(int(slhablock(1)), code)
	call SLHACopyString(len(text), text,
     &    slhablock(n*MSGLEN + MAXMSG + 3))
	end

