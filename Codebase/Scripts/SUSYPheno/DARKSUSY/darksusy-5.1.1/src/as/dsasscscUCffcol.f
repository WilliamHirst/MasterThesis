**************************************************************
*** SUBROUTINE dsasscscTCffcol                             ***
*** computes ASx and ASy coefficients for                  ***
*** scalar(1) + scalar(2) -> fermion(3) + fermion(4)       ***
*** for a neutralino or chargino in the U channel          ***
***                                                        ***
*** AUTHOR: Piero Ullio, ullio@sissa.it                    ***
*** Date: 01-03-03                                         ***
**************************************************************

      SUBROUTINE dsasscscUCffcol(kp1,kp2,kp3,kp4,kpchi)
      implicit none
      include 'dsmssm.h'
      include 'dsascom.h'
      integer kp1,kp2,kp3,kp4,kpchi
      complex*16 dsasdepro,den
      den=dsasdepro(Uvar,kpchi)
      if(kpchi.eq.kgluin) then
***** qluino in u-channel *************************************
        ASxplc(3,1)=ASxplc(3,1)
     &    +den*gr(kp2,kp3,kpchi)*gl(kp1,kp4,kpchi)
        ASxplc(3,4)=ASxplc(3,4)
     &    -den*gr(kp2,kp3,kpchi)*gl(kp1,kp4,kpchi)
        ASxprc(3,1)=ASxprc(3,1)
     &    +den*gl(kp2,kp3,kpchi)*gr(kp1,kp4,kpchi)
        ASxprc(3,4)=ASxprc(3,4)
     &    -den*gl(kp2,kp3,kpchi)*gr(kp1,kp4,kpchi)
        ASylc(3)=ASylc(3)
     &    +mass(kpchi)*den*gl(kp2,kp3,kpchi)
     &              *gl(kp1,kp4,kpchi)
        ASyrc(3)=ASyrc(3)
     &    +mass(kpchi)*den*gr(kp2,kp3,kpchi)
     &              *gr(kp1,kp4,kpchi)
      else  
***** chi in t-channel **+*************************************
        ASxplc(6,1)=ASxplc(6,1)
     &    +den*gr(kp2,kp3,kpchi)*gl(kp1,kp4,kpchi)
        ASxplc(6,4)=ASxplc(6,4)
     &    -den*gr(kp2,kp3,kpchi)*gl(kp1,kp4,kpchi)
        ASxprc(6,1)=ASxprc(6,1)
     &    +den*gl(kp2,kp3,kpchi)*gr(kp1,kp4,kpchi)
        ASxprc(6,4)=ASxprc(6,4)
     &    -den*gl(kp2,kp3,kpchi)*gr(kp1,kp4,kpchi)
        ASylc(6)=ASylc(6)
     &    +mass(kpchi)*den*gl(kp2,kp3,kpchi)
     &              *gl(kp1,kp4,kpchi)
        ASyrc(6)=ASyrc(6)
     &    +mass(kpchi)*den*gr(kp2,kp3,kpchi)
     &              *gr(kp1,kp4,kpchi)
      endif
      return
      end
