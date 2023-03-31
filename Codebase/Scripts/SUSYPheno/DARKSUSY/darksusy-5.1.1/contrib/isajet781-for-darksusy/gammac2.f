C-----------------------------------------------------------------
      SUBROUTINE GAMMAC2(T,CPOB,F)            
C-----------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 T,CPOB(17),F(17)
      COMMON / G3/ G3
      REAL*8 G3
      SAVE /G3/         
      COMMON /CHRGN/ MCHA(2),MSQU(3),MSQD(3),SMIX(6),RMIX(2),ABMIX(2)
      REAL*8 MCHA,MSQU,MSQD,SMIX,RMIX,ABMIX             
      SAVE /CHRGN/
      COMMON /SCALEQ/Q
      REAL*8 Q
      SAVE /SCALEQ/
      COMMON /BSGSM/ MZ,MW,MB,MC,MS,MT,MTAU,XW,S12,S23,S13,ALFAEM,SN2THW
      REAL*8 MZ,MW,MB,MC,MS,MT,MTAU,XW,S12,S23,S13,ALFAEM,SN2THW
      SAVE /BSGSM/
      INTEGER EN,NS,NF,NFT
      REAL*8 PI,FAC

      PI=4.*ATAN(1.d0)
      FAC=G3**2/(8.d0*PI**2)
      
      EN=3
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
C
      F(1)=FAC*(20./3.d0*CPOB(1)-8.d0*CPOB(2)+6.d0*CPOB(4)
     $         +4.d0*CPOB(5)+2.d0*CPOB(6))
      F(2)=FAC*(CPOB(1)+2./3.d0*CPOB(2)+2.*CPOB(4)
     $         +3./2.d0*CPOB(5)+CPOB(6)+1./2.d0*CPOB(7)
     $         +(EN**2-2)/8.d0/EN*CPOB(10)+1./4.d0*CPOB(11)
     $         +(EN**2-2)/8.d0/EN*CPOB(12)
     $         -((EN**2-2)/16.d0/EN+(EN**2+2)/72.d0/EN)*CPOB(13)
     $         -((EN**2-2)/16.d0/EN-(EN**2+2)/72.d0/EN)*CPOB(14)
     $         +1./4.d0*CPOB(15)-1./8.d0*CPOB(16)-1./8.d0*CPOB(17))
      F(3)=FAC*(-2.d0*CPOB(1)+4./3.d0*CPOB(2)+16./3.d0*CPOB(3)
     $         -CPOB(4)+CPOB(6)+2.*CPOB(7)+4./3.d0*CPOB(9)
     $         -2.*(EN**2-1)/4.d0/EN*CPOB(10)
     $         -2.*(EN**2-1)/4.d0/EN*CPOB(12)
     $         +10./18.*(EN**2-1)/2./EN*CPOB(13)
     $         +8./18.*(EN**2-1)/2./EN*CPOB(14))
      F(4)=FAC*(2./3.d0*CPOB(4)-113./36.d0*CPOB(5)-2.d0*CPOB(6)
     $          -113./36.d0*CPOB(7))    
      F(5)=FAC*(2.d0*CPOB(4)+137./18.d0*CPOB(5)+2.d0*CPOB(6)
     $         +89./18.d0*CPOB(7))
      F(6)=FAC*(-2.d0*CPOB(4)-113./36.d0*CPOB(5)+2./3.d0*CPOB(6)
     $          -113./36.d0*CPOB(7))
      F(7)=FAC*(-2.*CPOB(4)-4./3.d0*CPOB(5)-2.*CPOB(6)
     $          +4./3.d0*CPOB(7))
      F(8)=FAC*(9./4.d0*CPOB(5)+9./4.d0*CPOB(7))
      F(9)=0.
      F(10)=FAC*(1./2.d0*CPOB(10)+9./2.d0*CPOB(12)-1./2.d0*CPOB(13)
     $           -4.d0*CPOB(14)+3./2.d0*CPOB(16)-3./2.d0*CPOB(17))
      F(11)=FAC*(-3./2.d0*CPOB(10)-4.d0*CPOB(11)-3./2.d0*CPOB(12)
     $          +3./2.d0*CPOB(13)-1./2.d0*CPOB(16)+1./2.d0*CPOB(17))
      F(12)=FAC*(-4.d0*CPOB(12))
      F(13)=FAC*(-25./6.d0*CPOB(13)+1./6.d0*CPOB(14)
     $          +1./3.d0*SMIX(2)+2./3.d0*SMIX(3) 
     $          +NF/3.d0*SMIX(4)+NF/3.d0*SMIX(6)
     $          +1./12.d0*ABMIX(1)-1./12.d0*ABMIX(2))   
      F(14)=FAC*(-25./6.d0*CPOB(14)+1./6.d0*CPOB(13)   
     $          -1./3.d0*SMIX(2)-2./3.d0*SMIX(3) 
     $          -NF/3.d0*SMIX(4)-NF/3.d0*SMIX(6)   
     $          -1./12.d0*ABMIX(1)+1./12.d0*ABMIX(2))   
      F(15)=FAC*(-4.d0*CPOB(15))
      F(16)=FAC*(13./18.d0*CPOB(13)-13./18.d0*CPOB(14)
     $          -2.d0*CPOB(16)-2.d0*CPOB(17)
     $          -1./9.d0*SMIX(2)-2./9.d0*SMIX(3) 
     $          -NF/9.d0*SMIX(4)-NF/9.d0*SMIX(6)
     $          -1./36.d0*ABMIX(1)+1./36.d0*ABMIX(2))   
         F(17)=FAC*(-13./18.d0*CPOB(13)+13./18.d0*CPOB(14)
     $          -2.d0*CPOB(16)-2.d0*CPOB(17)
     $          +1./9.d0*SMIX(2)+2./9.d0*SMIX(3) 
     $          +NF/9.d0*SMIX(4)+NF/9.d0*SMIX(6)   
     $          +1./36.d0*ABMIX(1)-1./36.d0*ABMIX(2))   
      
      RETURN
      END
