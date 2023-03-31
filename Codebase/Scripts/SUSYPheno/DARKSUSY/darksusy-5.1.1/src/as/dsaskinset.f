**************************************************************
*** subroutine dsaskinset                                  ***
*** set askinder variables: the kinematic variables and    ***
*** the scalar products of the four vectors                ***
*** p(1), p(2), k(3) and k(4) related by                   ***
*** p(1) + p(2) = k(3) + k(4)                              ***
*** in the center of mass frame                            ***
*** input: asparmass, askin variables                      ***
***                                                        ***
*** AUTHOR: Piero Ullio, ullio@sissa.it                    ***
*** Date: 01-02-28                                         ***
**************************************************************

      subroutine dsaskinset
      implicit none
      include 'dsascom.h'
****
      ep1=dsqrt(mass1**2+p12**2)
      ep2=dsqrt(mass2**2+p12**2)
      Svar=mass1**2+mass2**2+2.d0*(ep1*ep2+p12**2)
      k34=dsqrt(dabs(Svar-(mass3+mass4)**2))
     &    *dsqrt(dabs(Svar-(mass3-mass4)**2))/dsqrt(Svar)/2.d0
      ek3=dsqrt(mass3**2+k34**2)
      ek4=dsqrt(mass4**2+k34**2)
****
      scal(1,1)=mass1**2
      scal(1,2)=(Svar-(mass1**2+mass2**2))/2.d0
      scal(1,3)=ep1*ek3-p12*k34*costheta
      scal(1,4)=ep1*ek4+p12*k34*costheta
      scal(2,1)=scal(1,2)
      scal(2,2)=mass2**2
      scal(2,3)=ep2*ek3+p12*k34*costheta
      scal(2,4)=ep2*ek4-p12*k34*costheta
      scal(3,1)=scal(1,3)
      scal(3,2)=scal(2,3)
      scal(3,3)=mass3**2
      scal(3,4)=(Svar-(mass3**2+mass4**2))/2.d0
      scal(4,1)=scal(1,4)
      scal(4,2)=scal(2,4)
      scal(4,3)=scal(3,4)
      scal(4,4)=mass4**2
****
      Tvar=mass1**2+mass3**2-2.d0*scal(1,3)
      Uvar=mass1**2+mass4**2-2.d0*scal(1,4)
****
      return
      end      








