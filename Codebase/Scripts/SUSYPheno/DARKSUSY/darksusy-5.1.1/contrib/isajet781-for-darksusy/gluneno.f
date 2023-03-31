C---------------------------------------------------------------------- 
      SUBROUTINE GLUNENO(G,GAM,AGE,AHE)
C----------------------------------------------------------------------
c    Computes 6x6 down squark diagonalization matrix and 
C    neutralino-(d)quark-(d)squark coupling matrices.
C
c    Ref: H.Anlauf  hep-ph/9406286;
c         S.Bertolini et al  NPB353, 591 (1991).
c
c    Created: H.Baer and M.Brhlik
C 
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C   Input:
C     G(i) - vector of RGE evolved parameters
C   Output:
C     GAM       - 6x6 down squark rotation matrix
C     AGE, AHE  - 4x6x6 neutralino-quark-squark matrices 
c                 defined by Eq(C.9) of Bertolini et al.
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv            
      IMPLICIT NONE
      REAL*8 G(111)
      REAL*8 GAM(6,6),AGE(4,6,6),AHE(4,6,6)
      COMMON /GLNN/ MCH0(4),MDSB(6)
      REAL*8 MCH0,MDSB
      SAVE /GLNN/
      COMMON /GGN/ M1,M2,M3,ABOT,ATOP,ATAU
      REAL*8 M1,M2,M3,ABOT,ATOP,ATAU
      SAVE /GGN/      
      COMMON /BSGSM/MZ,MW,MB,MC,MS,MT,MTAU,XW,S12,S23,S13,ALFAEM,SN2THW
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
      REAL*8 DM(6,6),ESQ(3,3),ESD(3,3),EMD(3,3),EAD(3,3),UNI(3,3),
     &       EDY(3,3),XORK(6),ANEU(4,4),BNEU(4,4)
c  ANEU(4,4) - neutralino mixing matrix on (photino,zino,higgs1,higgs2) basis
c              as defined by Haber-Kane
c  BNEU(4,4) - ANEU up to sign change in rows with negative mass eigenvalues
c
      REAL*8 BETA,SB,CB,S2B,C2B,VEV,SW,CW
      REAL*8 YDdiag(3,3),VDL(3,3),VDR(3,3),TD(3,3),M2Q(3,3),M2D(3,3)
      INTEGER I,J,K,L,IERR
      DATA UNI/1.d0,0.d0,0.d0, 0.d0,1.d0,0.d0, 0.d0,0.d0,1.d0/   ! 3x3 identity matrix
      
      BETA=ATAN(TANB)
      SB=SIN(BETA)
      CB=COS(BETA)       
      S2B=SIN(2.*BETA)
      C2B=COS(2.*BETA) 
      VEV=SQRT(V**2+VP**2)          
      SW=SQRT(XW)      
      CW=SQRT(1.d0-XW)
c
c  Convert neutralino mixing matrix to H-K
c
      DO I=1,4
        ANEU(I,1)= ZMIXSS(4,I)
        ANEU(I,2)= ZMIXSS(3,I)
        ANEU(I,3)=-ZMIXSS(2,I)
        ANEU(I,4)=-ZMIXSS(1,I)
	MCH0(I)=ABS(AMZISS(I))
      ENDDO
      DO I=1,4
        DO J=1,4
          BNEU(I,J)=ANEU(I,J)*SIGN(1.d0,-AMZISS(I))
        ENDDO
      ENDDO
c
c  6x6 d-squark mass matrix
c
      CALL YUKDIAG(G,13,YDdiag,VDL,VDR)
      K=0
      DO I=1,3
        DO J=1,3
	  TD(J,I)=G(43+K)
	  M2Q(I,J)=G(65+K)
	  M2D(I,J)=G(92+K)
	  K=K+1
        ENDDO
      ENDDO
c...Rotate matrices to superCKM basis      
      DO I=1,3
        DO J=1,3
          EMD(I,J)=YDdiag(I,J)*VEV*CB
	  EAD(I,J)=0.d0
	  ESQ(I,J)=0.d0
	  ESD(I,J)=0.d0
	  DO K=1,3
	    DO L=1,3
	      EAD(I,J)=EAD(I,J)+VDR(I,K)*TD(L,K)*VDL(J,L)*CB*VEV
              ESQ(I,J)=ESQ(I,J)+VDL(I,K)*M2Q(K,L)*VDL(J,L)
              ESD(I,J)=ESD(I,J)+VDR(I,K)*M2D(K,L)*VDR(J,L)
	    ENDDO
	  ENDDO
	ENDDO
      ENDDO
c...Coupling lambda given by Eq(47) of Anlauf
      EDY(1,1)=EMD(1,1)/(SQRT(2.d0)*CB*MW)
      EDY(1,2)=0.d0
      EDY(1,3)=0.d0
      EDY(2,1)=0.d0
      EDY(2,2)=EMD(2,2)/(SQRT(2.d0)*CB*MW)
      EDY(2,3)=0.d0
      EDY(3,1)=0.d0
      EDY(3,2)=0.d0
      EDY(3,3)=EMD(3,3)/(SQRT(2.d0)*CB*MW)
C...Construct 6x6 down squark mass^2 matrix
      DO I=1,3
      DO J=1,3
        DM(I,J)=0.d0
        DO K=1,3
          DM(I,J)=DM(I,J)+ESQ(I,K)*UNI(K,J)+EMD(K,I)*EMD(K,J)-
     $   	 -MZ**2*(1./2.d0-1./3.d0*XW)*C2B*UNI(I,K)*UNI(K,J)
        ENDDO
c	DM(I,J+3)=EAD(I,J)
c	DO K=1,3
c	  DM(I,J+3)=DM(I,J+3)-EMD(I,K)*MU*TANB*UNI(K,J)
c	ENDDO
c	DM(I+3,J)=EAD(J,I)
c	DO K=1,3
c	  DM(I+3,J)=DM(I+3,J)-MU*TANB*UNI(I,K)*EMD(J,K)
c	ENDDO
c sign in trilinear and mu fixed
	DM(I,J+3)=-EAD(I,J)
	DO K=1,3
	  DM(I,J+3)=DM(I,J+3)+EMD(I,K)*MU*TANB*UNI(K,J)
	ENDDO
	DM(I+3,J)=-EAD(J,I)
	DO K=1,3
	  DM(I+3,J)=DM(I+3,J)+MU*TANB*UNI(I,K)*EMD(J,K)
	ENDDO
        DM(I+3,J+3)=0.d0
        DO K=1,3
          DM(I+3,J+3)=DM(I+3,J+3)+ESD(I,K)*UNI(K,J)+EMD(K,I)*EMD(K,J)
     $  	     +MZ**2*1./3.d0*XW*C2B*UNI(I,K)*UNI(K,J)
        ENDDO
      ENDDO
      ENDDO
c...Diagonalization
      CALL EIGSYS(6,6,DM,MDSB,GAM,IERR,XORK)
      IF (IERR.NE.0) THEN
        PRINT*, 'EISRS1 ERROR IN GLUNENO, IERR=',IERR
        STOP99
      END IF
       
      DO I=1,6
        MDSB(I)=SQRT(abs(MDSB(I)))
      ENDDO
c
c  Neutralino-quark-squark coupling matrices
c
      DO I=1,4
      DO J=1,3
      DO K=1,6
        AGE(I,J,K)=((-1./3.d0)*SW*(ANEU(I,1)*CW+ANEU(I,2)*SW)
     $              +1./CW*(-1./2.d0+1./3.d0*XW)
     $               *(-ANEU(I,1)*SW+ANEU(I,2)*CW))*GAM(J,K)
        AGE(I,J+3,K)=((-1./3.d0)*SW*(BNEU(I,1)*CW+BNEU(I,2)*SW)
     $                +1./CW*(1./3.*XW)
     $                 *(-BNEU(I,1)*SW+BNEU(I,2)*CW))*GAM(J+3,K)
        AHE(I,J,K)=0.d0
        AHE(I,J+3,K)=0.d0
        DO L=1,3
          AHE(I,J,K)=AHE(I,J,K)+BNEU(I,3)*EDY(J,L)*GAM(L,K)
          AHE(I,J+3,K)=AHE(I,J,K)+ANEU(I,3)*EDY(J,L)*GAM(L+3,K)
        ENDDO
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END
