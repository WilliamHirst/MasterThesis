**************************************************************
*** SUBROUTINE dsasscscsSHffbcol                           ***
*** computes ASx and ASy coefficients for                  ***
*** scalar(1) + scalar*(2) -> fermion(3) + fermionbar(4)   ***
*** for a Higgs boson in the S channel                     ***
***                                                        ***
*** AUTHOR: Piero Ullio, ullio@sissa.it                    ***
*** Date: 01-03-03                                         ***
**************************************************************

      SUBROUTINE dsasscscsSHffbcol(kp1,kp2,kp3,kp4,kph)
      implicit none
      include 'dsmssm.h'
      include 'dsascom.h'
      integer kp1,kp2,kp3,kp4,kph
      complex*16 dsasdepro,den
      den=dsasdepro(Svar,kph)
***** H in s-channel **+***************************************
      ASylc(4)=ASylc(4)
     &  -den*gl(kph,kp1,kp2)*gl(kph,kp3,kp4)
      ASyrc(4)=ASyrc(4)
     &  -den*gr(kph,kp1,kp2)*gr(kph,kp3,kp4)
      return
      end
