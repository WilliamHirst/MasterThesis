**************************************************************
*** subroutine dsaskinset2                                 ***
*** sets: k34, ek3, ek4, Tvar, Uvar                        ***
*** you must call dsaskinset1 before calling dsaskinset2   *** 
*** input: mass3, mass4, costheta                          ***
***                                                        ***
*** AUTHOR: Piero Ullio, ullio@sissa.it                    ***
*** Date: 01-02-28                                         ***
**************************************************************

      subroutine dsaskinset2
      implicit none
      include 'dsascom.h'
****
      k34=dsqrt(dabs(Svar-(mass3+mass4)**2))
     &    *dsqrt(dabs(Svar-(mass3-mass4)**2))/dsqrt(Svar)/2.d0
      ek3=dsqrt(mass3**2+k34**2)
      ek4=dsqrt(mass4**2+k34**2)
      Tvar=mass1**2+mass3**2-2.d0*(ep1*ek3-p12*k34*costheta)
      Uvar=mass1**2+mass4**2-2.d0*(ep1*ek4+p12*k34*costheta)
****
      return
      end      








