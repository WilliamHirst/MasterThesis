**************************************************************
*** SUBROUTINE dsasscscsTCffbcol                           ***
*** computes ASx and ASy coefficients for                  ***
*** scalar(1) + scalar*(2) -> fermion(3) + fermionbar(4)   ***
*** for a neutralino or chargino in the T channel          ***
***                                                        ***
*** AUTHOR: Piero Ullio, ullio@sissa.it                    ***
*** Date: 01-03-03                                         ***
**************************************************************

      SUBROUTINE dsasscscsTCffbcol(kp1,kp2,kp3,kp4,kpchi)
      implicit none
      include 'dsmssm.h'
      include 'dsascom.h'
      integer kp1,kp2,kp3,kp4,kpchi
      complex*16 dsasdepro,den
      den=dsasdepro(Tvar,kpchi)
      if(kpchi.eq.kgluin) then
***** qluino in t-channel *************************************
        ASxplc(2,1)=ASxplc(2,1)
     &    +den*gr(kp1,kp3,kpchi)*conjg(gr(kp2,kp4,kpchi))
        ASxplc(2,3)=ASxplc(2,3)
     &    -den*gr(kp1,kp3,kpchi)*conjg(gr(kp2,kp4,kpchi))
        ASxprc(2,1)=ASxprc(2,1)
     &    +den*gl(kp1,kp3,kpchi)*conjg(gl(kp2,kp4,kpchi))
        ASxprc(2,3)=ASxprc(2,3)
     &    -den*gl(kp1,kp3,kpchi)*conjg(gl(kp2,kp4,kpchi))
        ASylc(2)=ASylc(2)
     &    -mass(kpchi)*den*gl(kp1,kp3,kpchi)
     &              *conjg(gr(kp2,kp4,kpchi))
        ASyrc(2)=ASyrc(2)
     &    -mass(kpchi)*den*gr(kp1,kp3,kpchi)
     &              *conjg(gl(kp2,kp4,kpchi))
      else  
***** chi in t-channel ****************************************
        ASxplc(5,1)=ASxplc(5,1)
     &    +den*gr(kp1,kp3,kpchi)*conjg(gr(kp2,kp4,kpchi))
        ASxplc(5,3)=ASxplc(5,3)
     &    -den*gr(kp1,kp3,kpchi)*conjg(gr(kp2,kp4,kpchi))
        ASxprc(5,1)=ASxprc(5,1)
     &    +den*gl(kp1,kp3,kpchi)*conjg(gl(kp2,kp4,kpchi))
        ASxprc(5,3)=ASxprc(5,3)
     &    -den*gl(kp1,kp3,kpchi)*conjg(gl(kp2,kp4,kpchi))
        ASylc(5)=ASylc(5)
     &    -mass(kpchi)*den*gl(kp1,kp3,kpchi)
     &              *conjg(gr(kp2,kp4,kpchi))
        ASyrc(5)=ASyrc(5)
     &    -mass(kpchi)*den*gr(kp1,kp3,kpchi)
     &              *conjg(gl(kp2,kp4,kpchi))
      endif
      return
      end
