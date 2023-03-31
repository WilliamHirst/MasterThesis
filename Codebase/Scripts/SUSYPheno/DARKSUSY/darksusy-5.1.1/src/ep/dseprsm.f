***********************************************************************
*** real function dseprsm gives the ratio of electro+positron flux in an
*** a+ cycle to that in an a- cycle as a function of positron energy
*** in gev.
*** from clem et al, apj 464 (1996) 507.
*** author: joakim edsjo, edsjo@physto.se
***********************************************************************

      real*8 function dseprsm(eep)
      implicit none

      real*8 eep

      dseprsm=max(min(0.45d0+0.17d0*log(eep),1.0d0),0.18d0)

      return
      end
