C-----------------------------------------------------------------
      SUBROUTINE GAMMAWB1(T,CI1,F)            
C-----------------------------------------------------------------
C
C     Right hand side of 1-loop RGEs for evolution of 
C     Wilson coefficients CI1 below M_W scale
C          dCI1_i/dT = F_i(CI1) ~ \gamma_ji * CI1_j
C
      IMPLICIT NONE
      REAL*8 T,CI1(8),F(8)
      COMMON / G3/ G3
      REAL*8 G3
      SAVE /G3/
      INTEGER NF,NU,ND
      REAL*8 PI,FAC
      
      NF=5
      NU=2
      ND=3
C     
      PI=4.*ATAN(1.d0)
      FAC=G3**2/(16.d0*PI**2)
C
      F(1)=FAC*(-2.d0*CI1(1)+6.d0*CI1(2))       
      F(2)=FAC*(6.d0*CI1(1)-2.d0*CI1(2))
      F(3)=FAC*(-2./9.d0*CI1(2)-22./9.d0*CI1(3)
     $         -2./9.d0*NF*CI1(4)-2./9.d0*NF*CI1(6))
      F(4)=FAC*(2./3.d0*CI1(2)+22./3.d0*CI1(3)
     $          +(-2.d0+2./3.d0*NF)*CI1(4)+2./3.d0*NF*CI1(6))    
      F(5)=FAC*(-2./9.d0*CI1(2)-4./9.d0*CI1(3)
     $         -2./9.d0*NF*CI1(4)+2.d0*CI1(5)-2./9.d0*NF*CI1(6))
      F(6)=FAC*(2./3.d0*CI1(2)+4./3.d0*CI1(3)+2./3.d0*NF*CI1(4)
     $          -6.d0*CI1(5)+(-16.d0+2./3.d0*NF)*CI1(6))
      F(7)=FAC*(416./81.d0*CI1(2)-464./81.d0*CI1(3)
     $          +(416./81.d0*NU-232./81.d0*ND)*CI1(4)
     $          +32./9.d0*CI1(5)
     $          +(-448./81.d0*NU+200./81.d0*ND)*CI1(6)
     $          +32./3.d0*CI1(7)-32./9.d0*CI1(8))
      F(8)=FAC*(3.d0*CI1(1)+70./27.d0*CI1(2)
     $          +(140./27.d0+3.d0*NF)*CI1(3)
     $          +(6.d0+70./27.d0*NF)*CI1(4)
     $          +(-14./3.d0-3.d0*NF)*CI1(5)
     $          +(-4.d0-119./27.d0*NF)*CI1(6)
     $          +28./3.d0*CI1(8))        
      RETURN
      END
