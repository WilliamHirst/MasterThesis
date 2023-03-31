

**********************************************************************
*** dsntearthmassint gives the mass of the earth in units of kg within
*** a sphere with of radius r meters.
*** in this routine, the actual integration is performed. for speed,
*** use dsntearthmass instead which uses a tabulation of this result.
*** author: joakim edsjo (edsjo@physto.se)
*** date: april 1, 1999
**********************************************************************

      real*8 function dsntearthmassint(r)
      implicit none

      real*8 r,dsntedfunc,dsf_int
      external dsntedfunc

c...integrate earth density

      dsntearthmassint=dsf_int(dsntedfunc,0.0d0,r,1.d-2)

      return
      end
