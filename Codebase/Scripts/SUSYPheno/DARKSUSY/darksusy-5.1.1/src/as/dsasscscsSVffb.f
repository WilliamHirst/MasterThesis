**************************************************************
*** SUBROUTINE dsasscscsSVffb                              ***
*** computes ASx and ASy coefficients for                  ***
*** scalar(1) + scalar*(2) -> fermion(3) + fermionbar(4)   ***
*** for a vector boson in the S channel                    ***
***                                                        ***
*** AUTHOR: Piero Ullio, ullio@sissa.it                    ***
*** Date: 01-03-03                                         ***
**************************************************************

      SUBROUTINE dsasscscsSVffb(kp1,kp2,kp3,kp4,kpv)
      implicit none
      include 'dsmssm.h'
      include 'dsascom.h'
      integer kp1,kp2,kp3,kp4,kpv
      complex*16 dsasdepro,den
      den=dsasdepro(Svar,kpv)
***** V in s-channel **+***************************************
      ASxpl(1)=ASxpl(1)
     &  +den*gl(kpv,kp2,kp1)*gl(kpv,kp3,kp4)
      ASxpl(2)=ASxpl(2)
     &  -den*gl(kpv,kp2,kp1)*gl(kpv,kp3,kp4)
      ASxpr(1)=ASxpr(1)
     &  +den*gr(kpv,kp2,kp1)*gr(kpv,kp3,kp4)
      ASxpr(2)=ASxpr(2)
     &  -den*gr(kpv,kp2,kp1)*gr(kpv,kp3,kp4)
      if(kpv.ne.kgamma) then
        ASxpl(1)=ASxpl(1)
     &   -den*gl(kpv,kp2,kp1)*gl(kpv,kp3,kp4)
     &    *(mass1**2-mass2**2)/mass(kpv)**2   
        ASxpl(2)=ASxpl(2)
     &   -den*gl(kpv,kp2,kp1)*gl(kpv,kp3,kp4)
     &    *(mass1**2-mass2**2)/mass(kpv)**2   
        ASxpr(1)=ASxpr(1)
     &   -den*gr(kpv,kp2,kp1)*gr(kpv,kp3,kp4)
     &    *(mass1**2-mass2**2)/mass(kpv)**2   
        ASxpr(2)=ASxpr(2)
     &   -den*gr(kpv,kp2,kp1)*gr(kpv,kp3,kp4)
     &    *(mass1**2-mass2**2)/mass(kpv)**2   
      endif  
      return
      end
