C-----------------------------------------------------------------
      SUBROUTINE GAMMAC1(T,CPOA,F)            
C-----------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 T,CPOA(12),F(12)
      COMMON /G3/ G3
      REAL*8 G3
      SAVE /G3/         
      COMMON /CHRGN/ MCHA(2),MSQU(3),MSQD(3),SMIX(6),RMIX(2),ABMIX(2)
      REAL*8 MCHA,MSQU,MSQD,SMIX,RMIX,ABMIX             
      SAVE /CHRGN/
      COMMON /SCALEQ/Q
      REAL*8 Q
      SAVE /SCALEQ/
      COMMON/BSGSM/ MZ,MW,MB,MC,MS,MT,MTAU,XW,S12,S23,S13,ALFAEM,SN2THW
      REAL*8 MZ,MW,MB,MC,MS,MT,MTAU,XW,S12,S23,S13,ALFAEM,SN2THW
      SAVE /BSGSM/
C
      INTEGER NS,NF,NFT
      REAL*8 PI,FAC,BB

      PI=4.*ATAN(1.d0)
      FAC=G3**2/(8.d0*PI**2)
      
      NS=0
      IF(MSQU(1).LT.Q) THEN
       NS=NS+1
      ENDIF
      IF(MSQU(2).LT.Q) THEN
       NS=NS+1
      ENDIF
      IF(MSQU(3).LT.Q) THEN
       NS=NS+10
      ENDIF             
      IF(Q.GT.MT) THEN 
       NF=6
       NFT=3-2*3
      ELSE
       NF=5
       NFT=3-2*2
      ENDIF
      BB=-11./2.d0+1./3.d0*DBLE(NF)+1./12.d0*DBLE(NS)
C
      F(1)=FAC*(20./3.d0*CPOA(1)-8.d0*CPOA(2)+6.d0*CPOA(4)
     $         +4.d0*CPOA(5)+2.d0*CPOA(6))
      F(2)=FAC*(CPOA(1)+2./3.d0*CPOA(2)+2.d0*CPOA(4)
     $         +3./2.d0*CPOA(5)+CPOA(6)+1./2.d0*CPOA(7))
      F(3)=FAC*(-2.d0*CPOA(1)+4./3.d0*CPOA(2)+16./3.d0*CPOA(3)
     $         -CPOA(4)+CPOA(6)+2.d0*CPOA(7)+4./3.d0*CPOA(9)
     $         +3.d0*CPOA(11))
      F(4)=FAC*(2./3.d0*CPOA(4)-113./36.d0*CPOA(5)-2.d0*CPOA(6)
     $         -113./36.d0*CPOA(7))    
      F(5)=FAC*(2.d0*CPOA(4)+137./18.d0*CPOA(5)+2.d0*CPOA(6)
     $         +89./18.*CPOA(7))
      F(6)=FAC*(-2.d0*CPOA(4)-113./36.d0*CPOA(5)+2./3.d0*CPOA(6)
     $         -113./36.d0*CPOA(7))
      F(7)=FAC*(-2.d0*CPOA(4)-4./3.d0*CPOA(5)-2.d0*CPOA(6)
     $          +4./3.d0*CPOA(7))
      F(8)=FAC*(9./4.d0*CPOA(5)+9./4.d0*CPOA(7)
     $         +2.d0*CPOA(12))
      F(9)=0.d0
      F(10)=FAC*(-2.d0*BB*CPOA(10))
      F(11)=FAC*(16./3.d0-2.d0*BB)*CPOA(11)
      F(12)=FAC*(-2.d0*BB*CPOA(12))   
      
      RETURN
      END
