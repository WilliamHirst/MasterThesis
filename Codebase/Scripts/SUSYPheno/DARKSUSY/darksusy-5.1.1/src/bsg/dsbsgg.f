      function dsbsgg(t)

***********************************************************************
* Function G(t) in (E.8) of Gambino and Misiak,                       *
* Nucl. Phys. B611 (2001) 338                                         *
* t must be a positive number                                         *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-03-05                   *
***********************************************************************

      implicit none
      real*8 t
      complex*16 dsbsgg,fct
      real*8 pi
      parameter (pi=3.141592653589793238d0)

      if (t.lt.4.d0) then
       fct=dcmplx(-2.d0*(atan(dsqrt(t/(4.d0-t))))**2,0.d0)
      else
       fct=dcmplx(-pi**2/2.d0
     &                 +2.d0*(log((dsqrt(t)+dsqrt(t-4))/2.d0))**2,
     &                 -2.d0*pi*log((dsqrt(t)+dsqrt(t-4))/2.d0))

      endif
      dsbsgg=fct
      return
      end



