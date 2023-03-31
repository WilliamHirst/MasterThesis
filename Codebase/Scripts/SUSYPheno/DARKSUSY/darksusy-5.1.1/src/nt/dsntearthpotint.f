**********************************************************************
*** dsntearthpotint gives the gravitational potential inside and outside
*** of the earth as a function of the radius r (in meters).
*** in this routine, the actual integration is performed. for speed,
*** use dsntearthpot instead which uses a tabulation of this result.
*** author: joakim edsjo (edsjo@physto.se)
*** date: april 1, 1999
**********************************************************************

      real*8 function dsntearthpotint(r)
      implicit none

      real*8 r,dsntepfunc,dsf_int,gn,dsntearthmass
      parameter(gn=6.67259d-11) ! m^3 kg^-1 s^-1
      external dsntepfunc

c...integrate earth density

      if (r.lt.6378.140d3) then
        dsntearthpotint=
     &    dsf_int(dsntepfunc,100.0d0,max(r,110.0d0),1.d-2)-
     &    1.123782d8
      else
        dsntearthpotint=-dsntearthmass(6378.140d3)*gn/(max(r,100.0d0))
      endif

      return
      end
