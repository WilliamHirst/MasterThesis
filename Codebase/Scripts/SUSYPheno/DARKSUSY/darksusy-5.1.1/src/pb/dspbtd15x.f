      real*8 function dspbtd15x(tp)
**********************************************************************
*** antiproton propagation according to various models
*** dspbtd15x is containment time in 10^15 sec
*** inputs:
***     tp - antiproton kinetic energy (gev)
*** from common blocks
***     pbpropmodel - 0 leaky box with energy dependent esc. time
***                   1 chardonnet et al diffusion
***                   2 bergstrom,edsjo,ullio diffusion
***                   3 bergstrom,edsjo,ullio diffusion
***                       but with the DC-like setup as in moskalenko 
***                       et al. ApJ 565 (2002) 280
*** author: paolo gondolo 99-07-13
*** modified: piero ullio 00-07-13
*** modified: piero ullio 04-01-22
**********************************************************************
      implicit none
      include 'dspbcom.h'
      include 'dsmpconst.h'
      include 'dspbprivate.h'
      real*8 tp,pp,difftime,dspbtd15char,dspbtd15beu

                            ! leaky box with energy dep. esc. time
      if (pbpropmodel.eq.0) then
c is this ok ??????????????????????? i took it from the old
c dshrpbardiff. which paper does it come from? pu
        pp=dsqrt(2*m_p*tp+tp**2)
        difftime=6.d0*(1+pp/3.d0)**(-0.6d0) ! 10^15 seconds

                            ! chardonnet et al 1996, bottino et al 1998
      elseif (pbpropmodel.eq.1) then
         difftime = dspbtd15char(tp)

                            ! bergstrom, edsjo & ullio 1999
      elseif (pbpropmodel.eq.2.or.pbpropmodel.eq.3) then
         difftime = dspbtd15beu(tp)


      else
         write(*,*) 'error in dspbtd15x: pbpropmodel out of range: ',
     &        pbpropmodel
         stop
      endif
      dspbtd15x=difftime
      return
      end












