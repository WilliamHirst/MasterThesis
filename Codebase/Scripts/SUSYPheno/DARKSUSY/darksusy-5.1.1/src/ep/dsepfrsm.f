***********************************************************************
*** real*8 function dsepfrsm solar modulates the positron fraction at a
*** given energy (eep). only gives results for a+ cycle.
*** input:  epfr - interstellar solar modulation fraction
***         eep - positron energy in gev
***         qa - solar modulation cycle. >0 for positrons in a+ cycle and
***             <0 for a- cycle.
*** output: dsepfrsm - solar modulated positron fraction
*** ref: clem et al, apj 464 (1997) 507.
*** author: joakim edsjo, edsjo@physto.se
***********************************************************************

      real*8 function dsepfrsm(epfr,eep,qa)
      implicit none

      real*8 epfr,eep,dseprsm,r
      integer qa

      r=dseprsm(eep)
      if (qa.gt.0) then
        dsepfrsm=(epfr**2*(r+1.0d0)-epfr)/(r*(2.0d0*epfr-1.0d0))
      else
        dsepfrsm=(epfr**2*(r+1.0d0)-epfr*r)/(2.0d0*epfr-1.0d0)
      endif

      return
      end
