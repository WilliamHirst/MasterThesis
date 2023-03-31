**********************************************************************
*** function dsepkt calculates the integrated flux of positrons
*** between energy ea and eb from neutralino annihilation in the halo.
*** NOTE. This routine uses the Kamionkowski and Turner expressions.
*** from prd 43(1991)1774. Only included for comparison.
*** units: cm^-2 sec^-1 sr^-1
*** author: joakim edsjo, edsjo@physto.se
*** date: 98-02-10
**********************************************************************

      real*8 function dsepkt(ea,eb,istat)
      implicit none

      include 'dshmcom.h'
      include 'dshrcom.h'
      include 'dshacom.h'

c----------------------------------------------------------------------
      real*8 ea,eb,dsepktdiff,dsf_int
      integer istat
      external dsepktdiff

c----------------------------------------------------------------------

      hristat=0
      dsepkt=dsf_int(dsepktdiff,ea,min(eb,hamwimp),1.d-2)
      istat=hristat

      return
      end
