      function dsphi27b(t)
***********************************************************************
* The 2nd integrand of \phi_27 in (E.3) of Gambino and Misiak,        *
* Nucl. Phys. B611 (2001) 338                                         *
* for the calculation of b --> s gamma                                *
* author:Mia Schelke, schelke@physto.se, 2003-03-10                   *
***********************************************************************
      implicit none
      real*8 z,t
      real*8 dsphi27b
      complex*16 dsbsgg

c     z=m_c^2/m_b^2
      z=0.22d0**2
      

      dsphi27b=(1.d0-z*t)*dreal(dsbsgg(t)+dcmplx(t/2))

      return
      end



