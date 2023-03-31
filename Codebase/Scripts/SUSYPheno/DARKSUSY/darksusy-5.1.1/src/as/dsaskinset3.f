**************************************************************
*** subroutine dsaskinset3                                 ***
*** sets the scalar products of the four vectors           ***
*** p(1), p(2), k(3) and k(4) related by                   ***
*** p(1) + p(2) = k(3) + k(4)                              ***
*** in the center of mass frame                            ***
*** you must call dsaskinset1 and dsaskinset2              *** 
*** before calling dsaskinset3                             *** 
***                                                        ***
*** AUTHOR: Piero Ullio, ullio@sissa.it                    ***
*** Date: 01-02-28                                         ***
**************************************************************

      subroutine dsaskinset3
      implicit none
      include 'dsascom.h'
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
      return
      end      








