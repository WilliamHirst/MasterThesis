C-----------------------------------------------------------------
      SUBROUTINE GAMMASM(T,CSM,F)            
C-----------------------------------------------------------------
C     Eq. 24 of Anlauf
      IMPLICIT NONE
      REAL*8 T,CSM(9),F(9)
      COMMON /G3/ G3         
      REAL*8 G3
      SAVE /G3/         
      
      REAL*8 PI,FAC
      PI=4.*ATAN(1.d0)
      FAC=G3**2/(8.d0*PI**2)
C
      F(1)=FAC*(20./3.d0*CSM(1)-8.d0*CSM(2)+6.d0*CSM(4)
     $         +4.d0*CSM(5)+2.d0*CSM(6))
      F(2)=FAC*(CSM(1)+2./3.d0*CSM(2)+2.d0*CSM(4)
     $         +3./2.d0*CSM(5)+CSM(6)+1./2.d0*CSM(7))
      F(3)=FAC*(-2.d0*CSM(1)+4./3.d0*CSM(2)+16./3.d0*CSM(3)
     $         -CSM(4)+CSM(6)+2.d0*CSM(7)+4./3.d0*CSM(9))
      F(4)=FAC*(2./3.d0*CSM(4)-113./36.d0*CSM(5)-2.*CSM(6)
     $         -113./36.d0*CSM(7))
      F(5)=FAC*(2.d0*CSM(4)+137./18.d0*CSM(5)+2.d0*CSM(6)
     $         +89./18.d0*CSM(7))
      F(6)=FAC*(-2.d0*CSM(4)-113./36.d0*CSM(5)+2./3.d0*CSM(6)
     $         -113./36.d0*CSM(7))
      F(7)=FAC*(-2.d0*CSM(4)-4./3.d0*CSM(5)-2.d0*CSM(6)
     $         +4./3.d0*CSM(7))
      F(8)=FAC*(9./4.d0*CSM(5)+9./4.d0*CSM(7))
      F(9)=0.d0
      
      RETURN
      END
