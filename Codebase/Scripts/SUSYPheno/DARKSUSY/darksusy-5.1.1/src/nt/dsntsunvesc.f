

**********************************************************************
*** dsntsunvesc gives the escape velocity in km/s as a function of
*** the radius r (in meters) from the sun's core.
*** author: joakim edsjo (edsjo@physto.se)
*** input: radius in m
*** output escape velocity in km/s
*** date: 2003-11-26
**********************************************************************

      real*8 function dsntsunvesc(r)
      implicit none

      real*8 r,dsntsunpot,vescsurf,phisurf
      parameter(vescsurf=617.57d0)  ! km/s
      parameter(phisurf=-1.9069d11) ! m^2 s^-2

      dsntsunvesc=vescsurf*sqrt(abs(dsntsunpot(r)/phisurf))

      return
      end
