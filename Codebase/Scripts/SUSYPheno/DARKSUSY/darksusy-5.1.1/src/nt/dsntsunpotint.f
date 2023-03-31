**********************************************************************
*** dsntsunpotint gives the gravitational potential inside and outside
*** of the sun as a function of the radius r (in meters).
*** in this routine, the actual integration is performed. for speed,
*** use dsntsunpot instead which uses a tabulation of this result.
*** author: joakim edsjo (edsjo@physto.se)
*** date: april 1, 1999
**********************************************************************

      real*8 function dsntsunpotint(r)
      implicit none
      include 'dssun.h'

      real*8 r,dsntspfunc,dsf_int,gn,dsntsunmass
      parameter(gn=6.67259d-11) ! m^3 kg^-1 s^-1
      external dsntspfunc

c...integrate sun density

      if (r.lt.r_sun) then
        dsntsunpotint=
     &    -gn * m_sun/r_sun  ! surface potential, -1.9069d11 m^2 s^-2
     &    -dsf_int(dsntspfunc,max(r,100.0d0),r_sun,1.d-3)
      else
        dsntsunpotint=-m_sun*gn/(max(r,100.0d0))
      endif

      return
      end
