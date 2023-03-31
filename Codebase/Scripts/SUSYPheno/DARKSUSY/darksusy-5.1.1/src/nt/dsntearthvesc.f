

**********************************************************************
*** dsntearthvesc gives the escape velocity in km/s as a function of
*** the radius r (in meters) from the earth's core.
*** author: joakim edsjo (edsjo@physto.se)
*** input: radius in m
*** output escape velocity in km/s
*** date: april 1, 1999
**********************************************************************

      real*8 function dsntearthvesc(r)
      implicit none

      real*8 r,dsntearthpot,vescsurf,phisurf
      parameter(vescsurf=11.2d0)
      parameter(phisurf=-.627705d8)

c      phisurf=dsntearthpot(6378.140d3)

      dsntearthvesc=vescsurf*sqrt(abs(dsntearthpot(r)/phisurf))

c...The values below are approximate average values in the core and mantle
c...These are the values that are used in the jkg review. Using these instead
c...gives similar capture rates (within about 1%). There is no reason using
c...these approximations here though, they are just included for comparison
c...with the jkg approximations.
c      if (r.le.3483.0d3) then
c        dsntearthvesc=11.2d0*sqrt(1.6d0)
c      else
c        dsntearthvesc=11.2d0*sqrt(1.2d0)
c      endif

      return
      end
