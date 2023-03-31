C---------------------------------------------------------------------- 
      SUBROUTINE CHARGINO(GK,HK)
C----------------------------------------------------------------------
C     Computes chargino-quark-squark interaction matrices as given by
C     Eq(50) of Anlauf.
C
      IMPLICIT NONE
      REAL*8 GK(2,3),HK(2,3)
      COMMON / G3/ G3
      REAL*8 G3
      SAVE /G3/
      COMMON /CHRGN/ MCHA(2),MSQU(3),MSQD(3),SMIX(6),RMIX(2),ABMIX(2)
      REAL*8 MCHA,MSQU,MSQD,SMIX,RMIX,ABMIX        
      SAVE /CHRGN/        
      COMMON /GGN/ M1,M2,M3,ABOT,ATOP,ATAU
      REAL*8 M1,M2,M3,ABOT,ATOP,ATAU
      SAVE /GGN/      
      COMMON/BSGSM/ MZ,MW,MB,MC,MS,MT,MTAU,XW,S12,S23,S13,ALFAEM,SN2THW
      REAL*8 MZ,MW,MB,MC,MS,MT,MTAU,XW,S12,S23,S13,ALFAEM,SN2THW
      SAVE /BSGSM/
      COMMON /BSGSUG/ TANB,V,VP,MSUSY,MU,MSTP1,MSTP2,MSCHL,
     &                MSCHR,MSUPL,MSEL,MSW1,MGLU,MHPLUS,MHA0,AMZISS(4),
     &                ZMIXSS(4,4),AMW1SS,AMW2SS,GAMMAL,GAMMAR,THETAT,
     &                MTQ,MBQ,MSTLQ,MSTRQ,MGUT,FNGUT,FTMT,XRHNIN(21),
     &                XGMIN(14),GGUTSS,XNUSUG(20),XAMIN(7),EPSNU,
     &                FTRHLD(3),MHUSMG,MHDSMG,
     &                INUHM,IAL3UN,LND5ON
      REAL*8 TANB,V,VP,MSUSY,MU,MSTP1,MSTP2,MSCHL,MSCHR,MSUPL,
     &       MSEL,MSW1,MGLU,MHPLUS,MHA0,AMZISS,ZMIXSS,AMW1SS,AMW2SS,
     &       GAMMAL,GAMMAR,THETAT,MTQ,MBQ,MSTLQ,MSTRQ,MGUT,
     &       FNGUT,FTMT,XRHNIN,XGMIN,GGUTSS,XNUSUG,XAMIN,EPSNU,FTRHLD,
     &       MHUSMG,MHDSMG
      INTEGER IAL3UN,INUHM
      LOGICAL LND5ON
      SAVE /BSGSUG/
C
      REAL*8 BETA,SB,CB,ALT,ALB,UM(2,2),VM(2,2),TM(2,2)
      
      BETA=ATAN(TANB)
      SB=SIN(BETA)
      CB=COS(BETA)       
      ALT=MTQ/(SQRT(2.)*MW*SB)
      ALB=MBQ/(SQRT(2.)*MW*CB)
c
c  Convert chargino mixing matrices to H-K
c
      VM(1,1) =-SIN(GAMMAR)                       *SIGN(1.d0,AMW1SS)
      VM(1,2) = COS(GAMMAR)                       *SIGN(1.d0,AMW1SS)
      VM(2,1) =-COS(GAMMAR)*SIGN(1.d0,TAN(GAMMAR))*SIGN(1.d0,AMW2SS)
      VM(2,2) =-SIN(GAMMAR)*SIGN(1.d0,TAN(GAMMAR))*SIGN(1.d0,AMW2SS) 
      UM(1,1) =-SIN(GAMMAL)                       
      UM(1,2) = COS(GAMMAL)                       
      UM(2,1) =-COS(GAMMAL)*SIGN(1.d0,TAN(GAMMAL))
      UM(2,2) =-SIN(GAMMAL)*SIGN(1.d0,TAN(GAMMAL))

      MCHA(1) = ABS(AMW1SS)
      MCHA(2) = ABS(AMW2SS)
c
c  Convert stop mixing matrix to H-K
c
      TM(1,1) = COS(THETAT) 
      TM(1,2) =-SIN(THETAT)
      TM(2,1) = SIN(THETAT)
      TM(2,2) = COS(THETAT)
      
      MSQU(1)=MSTP1
      MSQU(2)=MSTP2
      MSQU(3)=1./2.d0*(MSCHL+MSCHR) 
C....
      GK(1,1)= VM(1,1)*TM(1,1)-ALT*VM(1,2)*TM(1,2)
      GK(1,2)= VM(1,1)*TM(2,1)-ALT*VM(1,2)*TM(2,2)
      GK(1,3)= (ALT*VM(1,2))**2-(GK(1,1))**2-(GK(1,2))**2
      GK(2,1)= VM(2,1)*TM(1,1)-ALT*VM(2,2)*TM(1,2)
      GK(2,2)= VM(2,1)*TM(2,1)-ALT*VM(2,2)*TM(2,2)
      GK(2,3)= (ALT*VM(2,2))**2-(GK(2,1))**2-(GK(2,2))**2
      
      HK(1,1)= ALB*UM(1,2)*TM(1,1)
      HK(1,2)= ALB*UM(1,2)*TM(2,1)
      HK(1,3)=-GK(1,1)*HK(1,1)-GK(1,2)*HK(1,2)
      HK(2,1)= ALB*UM(2,2)*TM(1,1)
      HK(2,2)= ALB*UM(2,2)*TM(2,1)
      HK(2,3)=-GK(2,1)*HK(2,1)-GK(2,2)*HK(2,2)
C
      RETURN
      END
