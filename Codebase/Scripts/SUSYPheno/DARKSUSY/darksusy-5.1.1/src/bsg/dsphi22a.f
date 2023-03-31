      function dsphi22a(t)
***********************************************************************
* The first integrand of \phi_22 in (E.2) of Gambino and Misiak,      *
* Nucl. Phys. B611 (2001) 338                                         *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-03-10                   *
***********************************************************************
      implicit none
      real*8 z,t
      real*8 dsphi22a
      real*8 dsabsq
      complex*16 dsbsgg

c     z=m_c^2/m_b^2
      z=0.22d0**2
      
      dsphi22a=(1.d0-z*t)*dsabsq(dsbsgg(t)/dcmplx(t)+dcmplx(0.5d0))

      return
      end



