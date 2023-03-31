***********************************************************************
*** This routine selects which set of accelerator constraints to use
*** when dsacbnd is called. For available options, see dsacbnd.f
***********************************************************************
      subroutine dsacset(a)
      implicit none
      include 'dsaccom.h'
      character*(*) a
      aclabel = a
      return
      end

