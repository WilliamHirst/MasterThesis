**************************************************************
*** SUBROUTINE dsasscscTCff                                ***
*** computes ASx and ASy coefficients for                  ***
*** scalar(1) + scalar(2) -> fermion(3) + fermion(4)       ***
*** for a neutralino or chargino in the T channel          ***
***                                                        ***
*** AUTHOR: Piero Ullio, ullio@sissa.it                    ***
*** Date: 01-03-03                                         ***
**************************************************************

      SUBROUTINE dsasscscTCff(kp1,kp2,kp3,kp4,kpchi)
      implicit none
      include 'dsmssm.h'
      include 'dsascom.h'
      integer kp1,kp2,kp3,kp4,kpchi
      complex*16 dsasdepro,den
      den=dsasdepro(Tvar,kpchi)
***** chi in t-channel **+*************************************
      ASxpl(1)=ASxpl(1)
     &  -den*gr(kp1,kp3,kpchi)*gl(kp2,kp4,kpchi)
      ASxpl(3)=ASxpl(3)
     &  +den*gr(kp1,kp3,kpchi)*gl(kp2,kp4,kpchi)
      ASxpr(1)=ASxpr(1)
     &  -den*gl(kp1,kp3,kpchi)*gr(kp2,kp4,kpchi)
      ASxpr(3)=ASxpr(3)
     &  +den*gl(kp1,kp3,kpchi)*gr(kp2,kp4,kpchi)
      ASyl=ASyl
     &  +mass(kpchi)*den*gl(kp1,kp3,kpchi)
     &              *gl(kp2,kp4,kpchi)
      ASyr=ASyr
     &  +mass(kpchi)*den*gr(kp1,kp3,kpchi)
     &              *gr(kp2,kp4,kpchi)
      return
      end
