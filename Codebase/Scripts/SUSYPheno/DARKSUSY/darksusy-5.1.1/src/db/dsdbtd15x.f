      real*8 function dsdbtd15x(tp)
**********************************************************************
*** antideuteron propagation according to various models
*** dspbtd15x is containment time in 10^15 sec
*** inputs:
***     tp - kinetic energy per nucleon (gev)
*** from common blocks
***     pbpropmodel - 2 bergstrom,edsjo,ullio diffusion
***                   3 bergstrom,edsjo,ullio diffusion
***                       but with the DC-like setup as in moskalenko 
***                       et al. ApJ 565 (2002) 280
***
*** author: paolo gondolo 99-07-13
*** modified: piero ullio 00-07-13
**********************************************************************
      implicit none
      include 'dspbcom.h'
      include 'dspbprivate.h'
      real*8 tp,dsdbtd15beu,dsdbtd15beum,difftime

      if (pbpropmodel.eq.2.or.pbpropmodel.eq.3) then
         difftime = dsdbtd15beu(tp)
 
      else
         write(*,*) 'error in dsdbtd15x: pbpropmodel out of range: ',
     &        pbpropmodel
         stop
      endif
      dsdbtd15x=difftime
      return
      end












