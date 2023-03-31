**********************************************************************
*** function dssigmav returns the annihilation cross section
*** sigma v at p=0 for neutralino-neutralino annihilation.
*** if partch=0, the full sigma v is obtained and if partch>0, the
*** cross section to channel partch is obtained, where is defined
*** as follows:
***
***   partch   process
***   ------   -------
***        0   All processes, i.e. the full annihilation cross section
***        1   H1 H1
***        2   H1 H2
***        3   H2 H2
***        4   H3 H3
***        5   H1 H3
***        6   H2 H3
***        7   H- H+
***        8   H1 Z
***        9   H2 Z
***       10   H3 Z
***       11   W- H+ and W+ H-
***       12   Z0 Z0
***       13   W+ W-
***       14   nu_e nu_e-bar
***       15   e+ e-
***       16   nu_mu nu_mu-bar
***       17   mu+ mu-
***       18   nu_tau nu_tau-bar
***       19   tau+ tau-
***       20   u u-bar
***       21   d d-bar
***       22   c c-bar
***       23   s s-bar
***       24   t t-bar
***       25   b b-bar
***       26   gluon gluon (1-loop)
***       27   q q gluon (not implemented yet, put to zero)
***       28   gamma gamma (1-loop)
***       29   Z gamma (1-loop)
***
*** Units of returned cross section: cm^3 s^-1
**********************************************************************

      function dssigmav(partch)
      implicit none
      include 'dsmssm.h'
      real*8 dssigmav
      integer partch
      include 'dsandwcom.h'
      include 'dshacom.h'
      include 'dsidtag.h'
      include 'dsprep.h'
      real*8 dsandwdcosnn
      integer i
      
      if (newmodelsigmav) then
         wtot=dsandwdcosnn(0.0d0,0.0d0,kn(1),kn(1))
c JE fix: take away below eventually         
c         abr(1) = prtial(22)/wtot ! cc-bar
c         abr(2) = prtial(25)/wtot ! bb-bar
c         abr(3) = prtial(24)/wtot ! tt-bar
c         abr(4) = prtial(19)/wtot ! tau+ tau-
c         abr(5) = prtial(13)/wtot ! w+ w-
c         abr(6) = prtial(12)/wtot ! z0 z0
c         abr(7) = prtial(5)/wtot ! h10 h30
c         abr(8) = prtial(8)/wtot ! z0 h10
c         abr(9) = prtial(9)/wtot ! z0 h20
c         abr(10)= prtial(11)/wtot ! w+ h- / w- h+
c         abr(11)= prtial(6)/wtot ! h20 h30
c         abr(12)= prtial(26)/wtot ! gluon gluon
c         abr(13)= prtial(17)/wtot ! mu+ mu-
c         abr(14)= prtial(29)/wtot ! z gamma
         mx=mass(kn(1))
c... sigma v = w / (4*E_1^2) = wtot / (2*mx**2) since integration
c... over cos theta gives factor of 2.
         sigmav = 0.38937966d-27*3.d10*wtot/(2.d0*mx**2) ! in cm^3/s
         do i=1,29
            sigv(i) = prtial(i)/wtot*sigmav
         enddo
         newmodelsigmav=.false.
      endif
      
      dssigmav=sigmav
      if (partch.gt.0) dssigmav=sigv(partch)
      return

      end
