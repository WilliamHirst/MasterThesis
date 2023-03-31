**************************************************************
*** subroutine dsaskinset1                                 ***
*** sets: ep1, ep2, and Svar                               ***
*** input: mass1, mass2, p12                               ***
***                                                        ***
*** AUTHOR: Piero Ullio, ullio@sissa.it                    ***
*** Date: 01-02-28                                         ***
**************************************************************

      subroutine dsaskinset1
      implicit none
      include 'dsascom.h'
****
      ep1=dsqrt(mass1**2+p12**2)
      ep2=dsqrt(mass2**2+p12**2)
      Svar=mass1**2+mass2**2+2.d0*(ep1*ep2+p12**2)
****
      return
      end






