*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                            dsntcom.h                               ***
***         this piece of code is needed as a separate file          ***
c----------------------------------------------------------------------c
c  author: joakim edsjo (edsjo@physto.se), 2000-08-16

      real*8 tausu,csu,tauea,cea,ntarateea,ntaratesu,ceadk,ntmx,gtot10
      common /ntres/tausu,csu,tauea,cea,ntarateea,ntaratesu,ceadk,
     &  ntmx,gtot10

      real*8 veout
      integer ntcalcmet,ntlambda,nttab
      common /ntpara/veout,ntcalcmet,ntlambda,nttab

      save /ntres/,/ntpara/

*************************** end of dsntcom.h *****************************






