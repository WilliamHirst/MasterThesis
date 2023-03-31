**************************************************************
*** SUBROUTINE dsasscscsSVffbcol                           ***
*** computes ASx and ASy coefficients for                  ***
*** scalar(1) + scalar*(2) -> fermion(3) + fermionbar(4)   ***
*** for a vector boson in the S channel                    ***
***                                                        ***
*** AUTHOR: Piero Ullio, ullio@sissa.it                    ***
*** Date: 01-03-03                                         ***
**************************************************************

      SUBROUTINE dsasscscsSVffbcol(kp1,kp2,kp3,kp4,kpv)
      implicit none
      include 'dsmssm.h'
      include 'dsascom.h'
      integer kp1,kp2,kp3,kp4,kpv
      complex*16 dsasdepro,den
      den=dsasdepro(Svar,kpv)
      if(kpv.eq.kgluon) then
***** gluon in s-channel **************************************
        ASxplc(1,1)=ASxplc(1,1)
     &    +den*gl(kpv,kp2,kp1)*gl(kpv,kp3,kp4)
        ASxplc(1,2)=ASxplc(1,2)
     &    -den*gl(kpv,kp2,kp1)*gl(kpv,kp3,kp4)
        ASxprc(1,1)=ASxprc(1,1)
     &    +den*gr(kpv,kp2,kp1)*gr(kpv,kp3,kp4)
        ASxprc(1,2)=ASxprc(1,2)
     &    -den*gr(kpv,kp2,kp1)*gr(kpv,kp3,kp4)
      else  
***** V in s-channel **+***************************************
        ASxplc(4,1)=ASxplc(4,1)
     &    +den*gl(kpv,kp2,kp1)*gl(kpv,kp3,kp4)
        ASxplc(4,2)=ASxplc(4,2)
     &    -den*gl(kpv,kp2,kp1)*gl(kpv,kp3,kp4)
        ASxprc(4,1)=ASxprc(4,1)
     &    +den*gr(kpv,kp2,kp1)*gr(kpv,kp3,kp4)
        ASxprc(4,2)=ASxprc(4,2)
     &    -den*gr(kpv,kp2,kp1)*gr(kpv,kp3,kp4)
        if(kpv.ne.kgamma) then
          ASxplc(4,1)=ASxplc(4,1)
     &     -den*gl(kpv,kp2,kp1)*gl(kpv,kp3,kp4)
     &      *(mass1**2-mass2**2)/mass(kpv)**2   
          ASxplc(4,2)=ASxplc(4,2)
     &     -den*gl(kpv,kp2,kp1)*gl(kpv,kp3,kp4)
     &      *(mass1**2-mass2**2)/mass(kpv)**2   
          ASxprc(4,1)=ASxprc(4,1)
     &     -den*gr(kpv,kp2,kp1)*gr(kpv,kp3,kp4)
     &      *(mass1**2-mass2**2)/mass(kpv)**2   
          ASxprc(4,2)=ASxprc(4,2)
     &     -den*gr(kpv,kp2,kp1)*gr(kpv,kp3,kp4)
     &      *(mass1**2-mass2**2)/mass(kpv)**2   
        endif  
      endif  
      return
      end
