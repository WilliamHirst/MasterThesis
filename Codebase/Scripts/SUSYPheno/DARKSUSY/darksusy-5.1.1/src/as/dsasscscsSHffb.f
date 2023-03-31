**************************************************************
*** SUBROUTINE dsasscscsSHffb                              ***
*** computes ASx and ASy coefficients for                  ***
*** scalar(1) + scalar*(2) -> fermion(3) + fermionbar(4)   ***
*** for a Higgs boson in the S channel                     ***
***                                                        ***
*** AUTHOR: Piero Ullio, ullio@sissa.it                    ***
*** Date: 01-03-03                                         ***
**************************************************************

      SUBROUTINE dsasscscsSHffb(kp1,kp2,kp3,kp4,kph)
      implicit none
      include 'dsmssm.h'
      include 'dsascom.h'
      integer kp1,kp2,kp3,kp4,kph
      complex*16 dsasdepro,den
      den=dsasdepro(Svar,kph)
***** H in s-channel **+***************************************
      ASyl=ASyl
     &  -den*gl(kph,kp1,kp2)*gl(kph,kp3,kp4)
      ASyr=ASyr
     &  -den*gr(kph,kp1,kp2)*gr(kph,kp3,kp4)
      return
      end
