C----------------------------------------------------------------------
      SUBROUTINE RGE157(T,G,F)
C----------------------------------------------------------------------
C
C  Right hand side of the full set of 157 renormalization 
c  group equations
C              dG_i/dT = F_i(G), i=1..157
c  Sparticles decouple at muliple scales according to Castano et al.
C  This is SURG157 subrotine of ISAJET-M with some minor modifications.
c
C  NOTE: This sibroutine follows convention in which fermion masses 
c  	 are proportional to Yukawas, i.e. m ~ Y, while the rest 
c  	 uses convention with fermion masses proportional to transposed
c  	 Yukawas in gauge eigenbasis, i.e.  m ~ Y^T,
C  	 One has to transpose Yukasas and SSB trilinears when passing
C  	 into or from this subroutine.
C
c  Ref.: Martin & Vaughn PRD 50, 2282 (1994);
c        Castano, Piard, Ramond PRD 49, 4882 (1994)
c        with errors in Eqs (B15) and B(16) fixed;
c        S.Antusch et al hep-ph/0501272;
C        Casas & Ibarra hep-ph/0103065;
c        Casas et al PRD 63,097302 (2001).
c
c  Author: Azar Mustafayev
c  Created: 05/29/07
c  Modified: 13/12/07 by Azar Mustafayev 
c             - to speedup matrix multiplications moved to subroutines
C
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C  Local variables:
C     YU  = Y_u - up Yukawa coupling matrix
c     YD  = Y_d - down Yukawa coupling matrix
C     YE  = Y_e - lepton Yukawa coupling matrix 
C     YN  = Y_nu- neutrino Yukawa coupling matrix
C     MRHN= M_RHN - RH neutrino mass matrix 
C     KA  = \kappa- effective dim-5 operator from integration out RH neutrinos 
C     TU  = h_u - up soft-breaking trilinear coupling matrix
C     TD  = h_d - down soft-breaking trilinear coupling matrix
c     TE  = h_e - lepton soft-breaking trilinear coupling matrix 
C     TN  = h_nu- neutrino soft-breaking trilinear coupling matrix 
C     M2Q  = m^2_Q  - squark doublet mass^2 matrix
C     M2L  = m^2_L  - slepton doublet mass^2 matrix
C     M2U  = m^2_u  - up-squark mass^2 matrix
C     M2D  = m^2_d  - down-squark mass^2 matrix
C     M2E  = m^2_e  - slepton singlet mass^2 matrix
C     M2N  = m^2_nu - sneutrino singlet mass^2 matrix
C     YUD = (Y_u)^dagger - Hermitian conjugated Y_u
C     YDD = (Y_d)^dagger 
C     YED = (Y_e)^dagger 
C     YND = (Y_nu)^dagger 
C     TUD = (h_u)^dagger 
C     TDD = (h_d)^dagger 
C     TED = (h_e)^dagger 
C     TND = (h_nu)^dagger 
C
C     G(1) = g_1
C     G(2) = g_2
C     G(3) = g_3
C     G(4) = Y_u(1,1)
C     G(5) = Y_u(1,2)
C     G(6) = Y_u(1,3)
C     ...    ...
C     G(12) = Y_u(3,3)
C     G(13)-G(21) = Y_d
C     G(22)-G(30) = Y_e
C     G(31) = M_1
C     G(32) = M_2
C     G(33) = M_3
C     G(34)-G(42) = h_u
C     G(43)-G(51) = h_d
C     G(52)-G(60) = h_e
C     G(61)= mu
C     G(62)= B*mu
C     G(63)= m^2_Hu
C     G(64)= m^2_Hd
C     G(65)-G(73) = m^2_Q
C     G(74)-G(82) = m^2_L
C     G(83)-G(91) = m^2_u
C     G(92)-G(100)= m^2_d
C     G(101)-G(109)=m^2_e
C     G(110)= v_u
C     G(111)= v_d
C     G(112)-G(120) = Y_nu  
C     G(121)-G(129) = M_RHN 
C     G(130)-G(138) = h_nu
C     G(139)-G(147) = m^2_nu
C     G(148)-G(156) = \kappa
C     G(157) = \lambda   - SM quartic higgs coupling
C
C     X - auxiliary 250x3x3 matrix
C     X(200)-(250) - used for neutrino sector
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      IMPLICIT NONE
      REAL*8 T,G(157),F(157)
      REAL*8 TRACE,TR3X3
c
      COMMON /BSGDEC/MSQDEC(3),MSLDEC(3),MSUDEC(3),MSDDEC(3),
     &               MSEDEC(3),MRNDEC(3),IRHN
      REAL*8 MSQDEC,MSLDEC,MSUDEC,MSDDEC,MSEDEC,MRNDEC
      INTEGER IRHN
      SAVE /BSGDEC/
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
      COMMON /GGN/ M1,M2,M3,ABOT,ATOP,ATAU
c        M1,M2,M3  - gaugino masses at MZ
c        ABOT,ATOP,ATAU  - soft trilinear scalar couplings at MZ 
      REAL*8 M1,M2,M3,ABOT,ATOP,ATAU
      SAVE /GGN/      
c
      REAL*8 PI,FAC,S,X(250,3,3),F2(157),B1,B2,B3,
     &       BETA,SINB,COSB,Q,SP,SIG1,SIG2,SIG3,
     &       T1(3,3),T2(3,3),T3(3,3),T4(3,3),T5(3,3),T6(3,3),I3(3,3),
     &       YU(3,3),YD(3,3),YE(3,3),YN(3,3),TU(3,3),TD(3,3),TE(3,3),
     &       TN(3,3),M2Q(3,3),M2L(3,3),M2U(3,3),M2D(3,3),M2E(3,3),
     &       M2N(3,3),YUD(3,3),YDD(3,3),YED(3,3),YND(3,3),TUD(3,3),
     &       TDD(3,3),TED(3,3),TND(3,3),MRHN(3,3),KA(3,3)
      INTEGER STPSHU,STPHD,STPSHL,STPSHH,
     &        STPHL,STPHH,STPSB,STPSW,STPSG,
     &        STPSQ(3),STPSU(3),STPSD(3),STPSL(3),STPSE(3),ITWOLP
C   Step functions (=1 above and =0 below the scale):
C     STPSHU  = theta_sHu  - up higgsino
C     STPSHD  = theta_sHd  - down higgsino
C     STPSHL  = theta_sh   - light higgsino
C     STPSHH  = theta_sH   - heavy higgsino
C     STPHL   = theta_h    - light Higgs
C     STPHH   = theta_H    - heavy Higgs
C     STPSB   = theta_sB   - bino
C     STPSW   = theta_sW   - wino
C     STPSG   = theta_sg   - gluino
C     STPSQ(i)= theta_sQi  - i-th squark doublet
C     STPSU(i)= theta_sUi  - i-th up squark singlet
C     STPSD(i)= theta_sDi  - i-th down squark singlet
C     STPSL(i)= theta_sLi  - i-th slepton doublet
C     STPSE(i)= theta_sEi  - i-th slepton singlet
C     
C     ITWOLP  - switch for the second loop
c              =0 - 1-loop
c              =1 - only gauge and Yukawas at 2-loop
c              =2 - full 2-loop
      INTEGER NSQ,NSU,NSD,NSL,NSE,NSH,NH,NU,NE,ND,NN
      INTEGER I,J,K,L,M,N
      DATA I3/1.d0,0.d0,0.d0, 0.d0,1.d0,0.d0, 0.d0,0.d0,1.d0/   ! 3x3 identity matrix

      
      PI=4.d0*DATAN(1.d0)
      FAC=16.d0*PI**2

      ITWOLP=2    ! full 2-loop
      
      BETA=ATAN(TANB)
      SINB=SIN(BETA)
      COSB=SQRT(1.d0-SINB**2)
      Q=MGUT*DEXP(T)

c-----assemble matrices------------------------------------------------      
      K=0
      DO I=1,3
        DO J=1,3
	  YU(I,J)=G(4+K)
	  YD(I,J)=G(13+K)
	  YE(I,J)=G(22+K)
	  YN(I,J)=G(112+K)
	  MRHN(I,J)=G(121+K)
	  KA(I,J)=G(148+K)
	  TU(I,J)=G(34+K)
	  TD(I,J)=G(43+K)
	  TE(I,J)=G(52+K)
	  TN(I,J)=G(130+K)
	  M2Q(I,J)=G(65+K)
	  M2L(I,J)=G(74+K)
	  M2U(I,J)=G(83+K)
	  M2D(I,J)=G(92+K)
	  M2E(I,J)=G(101+K)
	  M2N(I,J)=G(139+K)
	  YUD(J,I)=YU(I,J)
	  YDD(J,I)=YD(I,J)
	  YED(J,I)=YE(I,J)
	  YND(J,I)=YN(I,J)
	  TUD(J,I)=TU(I,J)
	  TDD(J,I)=TD(I,J)
	  TED(J,I)=TE(I,J)
	  TND(J,I)=TN(I,J)
	  K=K+1
        ENDDO
      ENDDO
C----Reset auxilliary matrices
      DO K=1,250
        DO I=1,3
          DO J=1,3
            X(K,I,J)=0.d0
	  ENDDO
	ENDDO
      ENDDO
c---- CALCULATE 1-LOOP THRESHOLD EFFECTS ----------
      ND=3
      NE=3
      NN=3
      IF (Q.GE.MT) THEN
        NU=3
      ELSE
        NU=2
      END IF

      DO i=1,3
	IF (Q.GE.MSQDEC(i)) THEN
	  STPSQ(i)=1
        ELSE
	  STPSQ(i)=0
        ENDIF
	IF (Q.GE.MSLDEC(i)) THEN
	  STPSL(i)=1
        ELSE
	  STPSL(i)=0
        ENDIF
	IF (Q.GE.MSUDEC(i)) THEN
	  STPSU(i)=1
        ELSE
	  STPSU(i)=0
        ENDIF
	IF (Q.GE.MSDDEC(i)) THEN
	  STPSD(i)=1
        ELSE
	  STPSD(i)=0
        ENDIF
	IF (Q.GE.MSEDEC(i)) THEN
	  STPSE(i)=1
        ELSE
	  STPSE(i)=0
        ENDIF
      ENDDO
      NSU=STPSU(1)+STPSU(2)+STPSU(3)
      NSD=STPSD(1)+STPSD(2)+STPSD(3)
      NSQ=STPSQ(1)+STPSQ(2)+STPSQ(3)
      NSE=STPSE(1)+STPSE(2)+STPSE(3)
      NSL=STPSL(1)+STPSL(2)+STPSL(3)
      
      IF (Q.GE.DBLE(ABS(MU))) THEN
	STPSHL=1
	STPSHH=1
      ELSE 
	STPSHL=0
	STPSHH=0
      END IF
      NSH=STPSHL+STPSHH

      STPHL=1
      IF (Q.GT.MHPLUS) THEN
        STPHH=1
      ELSE
        STPHH=0
      END IF
      NH=STPHL+STPHH

      IF (Q.GE.MGLU) THEN
        STPSG=1
      ELSE
        STPSG=0
      END IF
	
      IF (Q.GE.M1) THEN
        STPSB=1
      ELSE
        STPSB=0
      END IF
      
      IF (Q.GE.M2) THEN
        STPSW=1
      ELSE
        STPSW=0
      END IF
	    
c-----evaluate auxiliary matrices-------------------------------------      
      CALL MPROD2X(X,1,YUD,YU)
      CALL MPROD2X(X,2,YDD,YD)
      CALL MPROD2X(X,3,YED,YE)
      CALL MPROD2X(X,11,YUD,TU)
      CALL MPROD2X(X,16,YDD,TD)
      CALL MPROD2X(X,17,YED,TE)
      CALL MPROD2X(X,24,TUD,TU)
      CALL MPROD2X(X,29,TDD,TD)
      CALL MPROD2X(X,30,TED,TE)
      CALL MPROD2X(X,37,TU,TUD)
      CALL MPROD2X(X,41,TD,TDD)
      CALL MPROD2X(X,45,TE,TED)
      CALL MPROD2X(X,46,YU,YUD)
      CALL MPROD2X(X,47,YD,YDD)
      CALL MPROD2X(X,48,YE,YED)
      CALL MPROD3X(X,4,YU,YUD,YU)
      CALL MPROD3X(X,5,YU,YDD,YD)
      CALL MPROD3X(X,6,YD,YDD,YD)
      CALL MPROD3X(X,7,YD,YUD,YU)
      CALL MPROD3X(X,8,YE,YED,YE)
      CALL MPROD3X(X,9,TU,YUD,YU)
      CALL MPROD3X(X,10,TU,YDD,YD)
      CALL MPROD3X(X,12,YU,YUD,TU)
      CALL MPROD3X(X,13,YU,YDD,TD)
      CALL MPROD3X(X,14,TD,YDD,YD)
      CALL MPROD3X(X,15,TD,YUD,YU)
      CALL MPROD3X(X,18,YD,YDD,TD)
      CALL MPROD3X(X,19,YD,YUD,TU)
      CALL MPROD3X(X,20,TE,YED,YE)
      CALL MPROD3X(X,21,YE,YED,TE)
      CALL MPROD3X(X,22,M2Q,YUD,YU)
      CALL MPROD3X(X,23,YUD,M2U,YU)
      CALL MPROD3X(X,25,M2Q,YDD,YD)
      CALL MPROD3X(X,26,YDD,M2D,YD)
      CALL MPROD3X(X,27,M2L,YED,YE)
      CALL MPROD3X(X,28,YED,M2E,YE)
      CALL MPROD3X(X,31,YUD,YU,M2Q)
      CALL MPROD3X(X,32,YDD,YD,M2Q)
      CALL MPROD3X(X,33,YED,YE,M2L)
      CALL MPROD3X(X,34,M2U,YU,YUD)
      CALL MPROD3X(X,35,YU,M2Q,YUD)
      CALL MPROD3X(X,36,YU,YUD,M2U)
      CALL MPROD3X(X,38,M2D,YD,YDD)
      CALL MPROD3X(X,39,YD,M2Q,YDD)
      CALL MPROD3X(X,40,YD,YDD,M2D)
      CALL MPROD3X(X,42,M2E,YE,YED)
      CALL MPROD3X(X,43,YE,M2L,YED)
      CALL MPROD3X(X,44,YE,YED,M2E)
C
C  1-loop part
C

c...gauge couplings
      B1=2.d0*(17.*NU/12.d0+5.*ND/12.d0+5.*NE/4.d0+1.*NN/4.d0)/5.d0
     &   +1.*NSQ/30.d0+4.*NSU/15.d0+1.*NSD/15.d0+1.*NSL/10.d0
     &   +1.*NSE/5.d0+1.*NSH/5.d0+1.*NH/10.d0
      B2=-22./3.d0+0.5d0*(NU+ND)+1.d0*(NE+NN)/6.d0
     $   +1.*NSQ/2.d0+1.*NSL/6.d0+1.*NSH/3.d0+1.*NH/6.d0+4.*STPSW/3.d0
      B3=2.*(NU+ND)/3.d0+1.*NSQ/3.d0+1.*NSU/6.d0+1.*NSD/6.d0
     $   +2.d0*STPSG-11.d0

      F(1)=G(1)**3*B1/FAC
      F(2)=G(2)**3*B2/FAC
      F(3)=G(3)**3*B3/FAC
c...Yukawa couplings
      DO 230 I=1,3
        DO 230 J=1,3
	  T1(I,J)=0.d0
	  T2(I,J)=0.d0
	  T3(I,J)=0.d0
	  T4(I,J)=0.d0
	  T5(I,J)=0.d0
	  T6(I,J)=0.d0
	  DO 230 K=1,3
	    T1(I,J)=T1(I,J)+YU(I,K)*STPSQ(K)*YUD(K,J)
	    T2(I,J)=T2(I,J)+YUD(I,K)*STPSU(K)*YU(K,J)
	    T3(I,J)=T3(I,J)+YDD(I,K)*STPSD(K)*YD(K,J)
	    T4(I,J)=T4(I,J)+YD(I,K)*STPSQ(K)*YDD(K,J)
	    T5(I,J)=T5(I,J)+YE(I,K)*STPSL(K)*YED(K,J)
	    T6(I,J)=T6(I,J)+YED(I,K)*STPSE(K)*YE(K,J)
230   CONTINUE	
      CALL MPROD2X(X,180,T1,YU)
      CALL MPROD2X(X,181,YU,T2)
      CALL MPROD2X(X,182,YU,T3)
      CALL MPROD2X(X,183,T4,YD)
      CALL MPROD2X(X,184,YD,T3)
      CALL MPROD2X(X,185,YD,T2)
      CALL MPROD2X(X,186,T5,YE)
      CALL MPROD2X(X,187,YE,T6)
      CALL MPROD3X(X,188,YU,YUD,YU)
      CALL MPROD3X(X,189,YU,YDD,YD)
      CALL MPROD3X(X,190,YD,YUD,YU)

      K=0
      DO 250 I=1,3
        DO 250 J=1,3
      F(4+K)=(YU(I,J)
     &        *(-3./5.d0*G(1)**2
     &  	 *(17./12.d0+3./4.d0*STPSHL
     &  	   -(1./36.d0*STPSQ(J)+4./9.d0*STPSU(I)
     &  	     +1./4.d0*STPSHL)*STPSB)
     &  	-G(2)**2*(9./4.d0+9./4.d0*STPSHL
     &  		  -3./4.d0*(STPSQ(J)+STPSHL)*STPSW)
     &  	-G(3)**2*(8.d0-4./3.d0*(STPSQ(J)+STPSU(I))*STPSG)
     &  	       +((SINB**2*STPHL+COSB**2)*3.d0*TRACE(X,1)
     &  	  +COSB**2*(STPHL-1)*3.d0*TRACE(X,2)
     &  	  +COSB**2*(STPHL-1)*TRACE(X,3)))
     &        +3./2.d0*(SINB**2*STPHL+COSB**2*STPHH)*X(188,I,J)
     &        +1./2.d0*(SINB**2*STPSHL+COSB**2*STPSHH)
     &         *(2.d0*X(180,I,J)+X(181,I,J))
     &        +1./2.d0*(COSB**2*(STPHL-4*(STPHL-STPHH))+SINB**2*STPHH)
     &         *X(189,I,J)
     &        +1./2.d0*(COSB**2*STPSHL+SINB**2*STPSHH)*X(182,I,J)
     &        )/FAC
      F(13+K)=(YD(I,J)
     &         *(-3./5.d0*G(1)**2
     &  	  *(5./12.d0+3./4.d0*STPSHL
     &  	    -(1./36.d0*STPSQ(J)+1./9.d0*STPSD(I)
     &  	      +1./4.d0*STPSHL)*STPSB)
     &  	 -G(2)**2*(9./4.d0+9./4.d0*STPSHL
     &  		   -3./4.d0*(STPSQ(J)+STPSHL)*STPSW)
     &  	 -G(3)**2*(8.d0-4./3.d0*(STPSQ(J)+STPSD(I))*STPSG)
     &  	 +(SINB**2*(STPHL-1)*3.d0*TRACE(X,1)
     &  	  +(COSB**2*STPHL+SINB**2)*3.d0*TRACE(X,2)
     &  	  +(COSB**2*STPHL+SINB**2)*TRACE(X,3)))
     &         +3./2.d0*(COSB**2*STPHL+SINB**2*STPHH)*X(6,I,J)
     &         +1./2.d0*(COSB**2*STPSHL+SINB**2*STPSHH)
     &         *(2.d0*X(183,I,J)+X(184,I,J))
     &         +1./2.d0*(SINB**2*(STPHL-4*(STPHL-STPHH))+COSB**2*STPHH)
     &  	*X(190,I,J)
     &         +1./2.d0*(SINB**2*STPSHL+COSB**2*STPSHH)*X(185,I,J)
     &         )/FAC
      F(22+K)=(YE(I,J)
     &         *(-3./5.d0*G(1)**2
     &  	  *(15./4.d0+3./4.d0*STPSHL
     &  	    -(1./4.d0*STPSL(J)+1.d0*STPSE(I)
     &  	      +1./4.d0*STPSHL)*STPSB)
     &  	 -G(2)**2*(9./4.d0+9./4.d0*STPSHL
     &  		   -3./4.d0*(STPSL(J)+STPSHL)*STPSW)
     &  	 +(SINB**2*(STPHL-1)*3.d0*TRACE(X,1)
     &  	  +(COSB**2*STPHL+SINB**2)*3.d0*TRACE(X,2)
     &  	  +(COSB**2*STPHL+SINB**2)*TRACE(X,3)))
     &  	+3./2.d0*(COSB**2*STPHL+SINB**2*STPHH)*X(8,I,J)
     &  	+1./2.d0*(COSB**2*STPSHL+SINB**2*STPSHH)
     &  	*(2.d0*X(186,I,J)+X(187,I,J))
     &         )/FAC
          K=K+1
250   CONTINUE
c...compute convenient quantity
      S=G(63)-G(64)+TR3X3(M2Q)-TR3X3(M2L)-2.d0*TR3X3(M2U)
     &  +TR3X3(M2D)+TR3X3(M2E)
      
c...gaugino masses      
      F(31)=2.d0*B1*G(1)**2*G(31)/FAC
      F(32)=2.d0*B2*G(2)**2*G(32)/FAC
      F(33)=2.d0*B3*G(3)**2*G(33)/FAC
c...higgs mixing parameters      
      F(61)=G(61)*(3.d0*TRACE(X,1)+3.d0*TRACE(X,2)+TRACE(X,3)
     &             -3.d0*G(2)**2-3./5.d0*G(1)**2)/FAC
      F(62)=(G(62)*(3.d0*TRACE(X,1)+3.d0*TRACE(X,2)+TRACE(X,3)
     &              -3.d0*G(2)**2-3./5.d0*G(1)**2)
     &       +G(61)*(6.d0*TRACE(X,11)+6.d0*TRACE(X,16)+2.d0*TRACE(X,17)
     &               +6.d0*G(2)**2*G(32)+6./5.d0*G(1)**2*G(31)))/FAC
c...higgs mass^2
      F(63)=(6.d0*(G(63)*TRACE(X,1)+TRACE(X,31)+TRACE(X,34)+TRACE(X,24))
     &       -6.d0*G(2)**2*G(32)**2-6./5.d0*G(1)**2*G(31)**2
     &       +3./5.d0*S*G(1)**2)/FAC
      F(64)=(6.d0*(G(64)*TRACE(X,2)+TRACE(X,25))+6.d0*TRACE(X,26)
     &       +2.d0*(G(64)*TRACE(X,3)+TRACE(X,27))+2.d0*TRACE(X,28)
     &       +6.d0*TRACE(X,29)+2.d0*TRACE(X,30)
     &       -6.d0*G(2)**2*G(32)**2-6./5.d0*G(1)**2*G(31)**2
     &       -3./5.d0*S*G(1)**2)/FAC

      K=0
      DO 300 I=1,3
        DO 300 J=1,3
c...soft trilinear scalar couplings
      F(34+K)=(TU(I,J)*3.d0*TRACE(X,1)+5.d0*X(9,I,J)+X(10,I,J)
     &         -TU(I,J)*(16./3.d0*G(3)**2+3.*G(2)**2+13./15.d0*G(1)**2)
     &         +YU(I,J)*6.d0*TRACE(X,11)+4.d0*X(12,I,J)+2.d0*X(13,I,J)
     &         +YU(I,J)*(32./3.d0*G(3)**2*G(33)+6.*G(2)**2*G(32)
     &                   +26./15.d0*G(1)**2*G(31)))/FAC
      F(43+K)=(TD(I,J)*(3.*TRACE(X,2)+TRACE(X,3))
     &         +5.d0*X(14,I,J)+X(15,I,J)
     &         -TD(I,J)*(16./3.d0*G(3)**2+3.*G(2)**2+7./15.d0*G(1)**2)
     &         +YD(I,J)*(6.*TRACE(X,16)+2.*TRACE(X,17))
     &         +4.d0*X(18,I,J)+2.d0*X(19,I,J)
     &         +YD(I,J)*(32./3.d0*G(3)**2*G(33)+6.d0*G(2)**2*G(32)
     &                   +14./15.d0*G(1)**2*G(31)))/FAC
      F(52+K)=(TE(I,J)*(3.d0*TRACE(X,2)+TRACE(X,3))+5.d0*X(20,I,J)
     &         -TE(I,J)*(3.d0*G(2)**2+9./5.d0*G(1)**2)
     &         +YE(I,J)*(6.*TRACE(X,16)+2.*TRACE(X,17))+4.*X(21,I,J)
     &         +YE(I,J)*(6.*G(2)**2*G(32)+18./5.d0*G(1)**2*G(31)))/FAC
c...the rest of SSB mass^2
      F(65+K)=(X(22,I,J)+2.*G(63)*X(1,I,J)+X(25,I,J)+2.*G(64)*X(2,I,J)
     &         +X(31,I,J)+X(32,I,J)+2.d0*X(23,I,J)+2.d0*X(26,I,J)
     &         +2.d0*X(24,I,J)+2.d0*X(29,I,J)
     &         +I3(I,J)*(-32./3.d0*G(3)**2*G(33)**2-6.*G(2)**2*G(32)**2
     &                   -2./15.d0*G(1)**2*G(31)**2
     &                   +G(1)**2*S/5.d0))/FAC
      F(74+K)=(X(27,I,J)+2.*G(64)*X(3,I,J)+2.*X(28,I,J)+X(33,I,J)
     &         +2.d0*X(30,I,J)
     &         +I3(I,J)*(-6.*G(2)**2*G(32)**2-6./5.d0*G(1)**2*G(31)**2
     &                   -3./5.d0*G(1)**2*S))/FAC
      F(83+K)=(2.*X(34,I,J)+4.*G(63)*X(46,I,J)+4.*X(35,I,J)
     &         +2.d0*X(36,I,J)+4.d0*X(37,I,J)
     &         +I3(I,J)*(-32./3.d0*G(3)**2*G(33)**2
     &                   -32./15.d0*G(1)**2*G(31)**2
     &                   -4./5.d0*G(1)**2*S))/FAC
      F(92+K)=(2.*X(38,I,J)+4.*G(64)*X(47,I,J)+4.*X(39,I,J)
     &         +2.*X(40,I,J)+4.*X(41,I,J)
     &         +I3(I,J)*(-32./3.d0*G(3)**2*G(33)**2
     &                   -8./15.d0*G(1)**2*G(31)**2
     &                   +2./5.*G(1)**2*S))/FAC
      F(101+K)=(2.*X(42,I,J)+4.*G(64)*X(48,I,J)+4.*X(43,I,J)
     &          +2.*X(44,I,J)+4.*X(45,I,J)
     &          +I3(I,J)*(-24./5.d0*G(1)**2*G(31)**2
     &                    +6./5.d0*G(1)**2*S))/FAC
	  K=K+1
300   CONTINUE

c...VEVs
      F(110)=-G(110)*(3.*TRACE(X,1)-3./4.d0*(G(2)**2+G(1)**2/5.d0))/FAC
      F(111)=-G(111)*(3.*TRACE(X,2)+TRACE(X,3)
     &  	      -3./4.d0*(G(2)**2+G(1)**2/5.d0))/FAC

c--------- neutrino sector ------------------------------
      IF (IRHN.GT.0) THEN
        IF(Q.LT.MIN(ABS(MRNDEC(1)),ABS(MRNDEC(2)),ABS(MRNDEC(3)))) THEN
	  DO I=112,147
	    F(I)=0.d0     ! all Majorana RHN are decoupled
	  ENDDO
        ELSE
	  CALL MPROD2X(X,200,YND,YN)
	  CALL MPROD2X(X,201,YN,YND)
	  CALL MPROD2X(X,215,YND,TN)
	  CALL MPROD2X(X,225,TN,TND)
	  CALL MPROD2X(X,229,TND,TN)
	  CALL MPROD3X(X,202,YE,YND,YN)
	  CALL MPROD3X(X,203,YN,YND,YN)
	  CALL MPROD3X(X,204,YN,YED,YE)
	  CALL MPROD3X(X,205,YN,YND,MRHN)
	  CALL MPROD3X(X,216,YN,YND,TN)
	  CALL MPROD3X(X,217,TN,YND,YN)
	  CALL MPROD3X(X,218,YN,YED,TE)
	  CALL MPROD3X(X,219,TN,YED,YE)
	  CALL MPROD3X(X,220,YE,YND,TN)
	  CALL MPROD3X(X,221,TE,YND,YN)
	  CALL MPROD3X(X,222,M2N,YN,YND)
	  CALL MPROD3X(X,223,YN,YND,M2N)
	  CALL MPROD3X(X,224,YN,M2L,YND)
	  CALL MPROD3X(X,226,M2L,YND,YN)
	  CALL MPROD3X(X,227,YND,YN,M2L)
	  CALL MPROD3X(X,228,YND,M2N,YN)
	  CALL MPROD3X(X,230,YND,M2L,YN)
	  
	  K=0
          DO 320 I=1,3
            DO 320 J=1,3
c...contributions to up-quark and charged lepton Yukawas
      F(4+K)=F(4+K) + YU(I,J)*TRACE(X,200)/FAC
      F(22+K)=F(22+K) + X(202,I,J)/FAC
c...neutrino Yukawa matrix
      F(112+K)=(3.d0*X(203,I,J)+X(204,I,J)
     &          +YN(I,J)*(TRACE(X,200)+3.d0*TRACE(X,1)
     &                    -3./5.d0*G(1)**2-3.d0*G(2)**2))/FAC
c...Majorana mass matrix
      F(121+K)=2.d0*(X(205,I,J)+X(205,J,I))/FAC
c...neutrino soft trilinear scalar couplings
      F(130+K)=(TN(I,J)*(3.d0*TRACE(X,1)+TRACE(X,200)
     &                   -3./5.d0*G(1)**2-3.d0*G(2)**2)
     &         +2.d0*YN(I,J)*(3.d0*TRACE(X,11)+TRACE(X,215)
     &                       -3./5.d0*G(1)**2*G(31)-3.d0*G(2)**2*G(32))
     &         +4.d0*X(216,I,J)+5.d0*X(217,I,J)+2.d0*X(218,I,J)
     &         +X(219,I,J))/FAC
c...contributions to the other soft trilinear scalar couplings
      F(34+K)=F(34+K)
     &        +(TU(I,J)*TRACE(X,200)+2.d0*YU(I,J)*TRACE(X,215))/FAC
      F(52+K)=F(52+K)+(2.d0*X(220,I,J)+X(221,I,J))/FAC
c...neutrino mass^2
      F(139+K)=(2.d0*(X(222,I,J)+X(223,I,J))
     &          +4.d0*(X(224,I,J)+G(63)*X(201,I,J)+X(225,I,J)))/FAC
c...contributions to slepton doublet mass^2
      F(74+K)= F(74+K)
     &         +(X(226,I,J)+X(227,I,J)
     &           +2.d0*(X(228,I,J)+G(63)*X(200,I,J)+X(229,I,J)))/FAC
	      K=K+1
320       CONTINUE        
C...contribution to up-higgs mass^2
        F(63)=F(63)+(2.d0*(TRACE(X,230)+TRACE(X,228)+G(63)*TRACE(X,200)
     &                     +TRACE(X,229)))/FAC
        ENDIF
c...effective dim-5 neutrino operator
        IF(LND5ON.AND.
     &     Q.LE.MAX(ABS(MRNDEC(1)),ABS(MRNDEC(2)),ABS(MRNDEC(3)))) THEN
          CALL MPROD3X(X,231,KA,YED,YE)
          CALL MPROD3X(X,232,KA,YND,YN)
          CALL MPROD4X(X,233,YED,YE,YED,YE)
          CALL MPROD4X(X,234,YND,YN,YND,YN)
          CALL MPROD4X(X,235,YDD,YD,YDD,YD)
          CALL MPROD4X(X,236,YUD,YU,YUD,YU)

	  IF (Q.LT.MSUSY) THEN     !  use SM RGEs
	    K=0
            DO I=1,3
              DO J=1,3
      F(148+K)=(-3./2.d0*(X(231,J,I)+X(231,I,J))*COSB**2
     &          +1./2.d0*(X(232,J,I)+X(232,I,J))*SINB**2
     &          +KA(I,J)*(2.d0*TRACE(X,3)*COSB**2
     &                    +2.d0*TRACE(X,200)*SINB**2
     &                    +6.d0*TRACE(X,1)*SINB**2
     &                    +6.d0*TRACE(X,3)*COSB**2
     &                    -3.d0*G(2)**2+G(157)))/FAC
	        K=K+1
	      ENDDO
	    ENDDO
c...SM quartic higgs coupling
      F(157)=(6.d0*G(157)**2-3.d0*G(157)*(3.d0*G(2)**2+3./5.d0*G(1)**2)
     &        +3.d0*G(2)**4+3./2.d0*(3./5.d0*G(1)**2+G(2)**2)**2
     &        +4.d0*G(157)*(TRACE(X,3)+TRACE(X,200)+3.d0*TRACE(X,2)
     &                      +3.d0*TRACE(X,1))
     &        -8.d0*(TRACE(X,233)+TRACE(X,234)+3.d0*TRACE(X,235)
     &               +3.d0*TRACE(X,236)))/FAC
	  ELSE                     !  use MSSM RGEs
	    K=0
            DO I=1,3
              DO J=1,3
      F(148+K)=(X(231,J,I)+X(231,I,J)+X(232,J,I)+X(232,I,J)
     &         +KA(I,J)*(2.d0*TRACE(X,200)+6.d0*TRACE(X,1)
     &                   -6./5.d0*G(1)**2-6.d0*G(2)**2))/FAC
	        K=K+1
	      ENDDO
	    ENDDO
	    F(157)=0.d0
	  ENDIF
	ELSE
	  DO I=148,157
	    F(I)=0.d0
	  ENDDO
	ENDIF
      ENDIF
C
C  2-loop part for gauge and Yukawa couplings
C
      IF (ITWOLP.GE.1) THEN
c...evaluate auxiliary matrices
        CALL MPROD4X(X,49,YU,YUD,YU,YUD)
        CALL MPROD4X(X,50,YU,YDD,YD,YUD)
        CALL MPROD4X(X,54,YD,YDD,YD,YDD)
        CALL MPROD4X(X,55,YE,YED,YE,YED)
        CALL MPROD5X(X,51,YU,YUD,YU,YUD,YU)
        CALL MPROD5X(X,52,YU,YDD,YD,YDD,YD)
        CALL MPROD5X(X,53,YU,YDD,YD,YUD,YU)
        CALL MPROD5X(X,56,YD,YDD,YD,YDD,YD)
        CALL MPROD5X(X,57,YD,YUD,YU,YUD,YU)
        CALL MPROD5X(X,58,YD,YUD,YU,YDD,YD)
        CALL MPROD5X(X,59,YE,YED,YE,YED,YE)

        IF (Q.LT.MSUSY) THEN	  ! use SM RGEs below M_SUSY
c...Gauge couplings
      F2(1)=G(1)**3
     &      *(199./50.d0*G(1)**2+27./10.d0*G(2)**2+44./5.d0*G(3)**2
     &        -17./10.d0*SINB**2*TRACE(X,1)-1./2.d0*COSB**2*TRACE(X,2)
     &        -3./2.d0*COSB**2*TRACE(X,3))
      F2(2)=G(2)**3
     &      *(9./10.d0*G(1)**2+35./6.d0*G(2)**2+12.d0*G(3)**2
     &        -3./2.d0*SINB**2*TRACE(X,1)-3./2.d0*COSB**2*TRACE(X,2)
     &        -1./2.d0*COSB**2*TRACE(X,3))
      F2(3)=G(3)**3
     &      *(11./10.d0*G(1)**2+9./2.d0*G(2)**2-26.d0*G(3)**2
     &        -2.d0*SINB**2*TRACE(X,1)-2.d0*COSB**2*TRACE(X,2))
c...Yukawa couplings
          K=0
	  DO I=1,3
            DO J=1,3
              F2(4+K)=0.d0	 ! FOR NOW a la standard ISAJET
              F2(13+K)=0.d0
	      F2(22+K)=0.d0
	      K=K+1
	    ENDDO
	  ENDDO	
        ELSE                       ! use MSSM RGEs above M_SUSY
c...Gauge couplings
      F2(1)=G(1)**3
     &      *(199./25.d0*G(1)**2+27./5.d0*G(2)**2+88./5.d0*G(3)**2
     &        -26./5.d0*TRACE(X,1)-14./5.d0*TRACE(X,2)
     &        -18./5.d0*TRACE(X,3))
      F2(2)=G(2)**3
     &      *(9./5.d0*G(1)**2+25.d0*G(2)**2+24.d0*G(3)**2
     &        -6.d0*TRACE(X,1)-6.d0*TRACE(X,2)-2.d0*TRACE(X,3))
      F2(3)=G(3)**3
     &      *(11./5.d0*G(1)**2+9.d0*G(2)**2+14.d0*G(3)**2
     &        -4.d0*TRACE(X,1)-4.d0*TRACE(X,2))
c...Yukawa couplings
          K=0
	  DO 325 I=1,3
            DO 325 J=1,3
      F2(4+K)=YU(I,J)*(-3.d0)*(3.d0*TRACE(X,49)+TRACE(X,50))
     &        -X(5,I,J)*(3.d0*TRACE(X,2)+TRACE(X,3))
     &        -9.d0*X(4,I,J)*TRACE(X,1)
     &        -4.d0*X(51,I,J)-2.d0*X(52,I,J)-2.d0*X(53,I,J)
     &        +YU(I,J)*(16.d0*G(3)**2+4./5.d0*G(1)**2)*TRACE(X,1)
     &        +(6.d0*G(2)**2+2./5.d0*G(1)**2)*X(4,I,J)
     &        +2./5.d0*G(1)**2*X(5,I,J)
     &        +YU(I,J)*(-16./9.d0*G(3)**4+8.d0*G(3)**2*G(2)**2
     &          	+136./45.d0*G(3)**2*G(1)**2
     &          	+15./2.d0*G(2)**4+G(2)**2*G(1)**2
     &          	+2743./450.d0*G(1)**4)
      F2(13+K)=YD(I,J)*(-3.d0)*(3.d0*TRACE(X,54)+TRACE(X,50)
     &          		+TRACE(X,55))
     &         -3.d0*X(7,I,J)*TRACE(X,1)
     &         -3.d0*X(6,I,J)*(3.d0*TRACE(X,2)+TRACE(X,3))
     &         -4.d0*X(56,I,J)-2.d0*X(57,I,J)-2.d0*X(58,I,J)
     &         +YD(I,J)*((16.d0*G(3)**2-2./5.d0*G(1)**2)*TRACE(X,2)
     &                   +6./5.d0*G(1)**2*TRACE(X,3))
     &         +4./5.d0*G(1)**2*X(7,I,J)
     &         +(6.d0*G(2)**2+4./5.d0*G(1)**2)*X(6,I,J)
     &         +YD(I,J)*(-16./9.d0*G(3)**4+8.d0*G(3)**2*G(2)**2
     &          	 +8./9.d0*G(3)**2*G(1)**2
     &          	 +15./2.d0*G(2)**4+G(2)**2*G(1)**2
     &          	 +287./90.d0*G(1)**4)
      F2(22+K)=YE(I,J)*(-3.d0)*(3.d0*TRACE(X,54)+TRACE(X,50)
     &          		+TRACE(X,55))
     &         -3.d0*X(8,I,J)*(3.d0*TRACE(X,2)+TRACE(X,3))
     &         -4.d0*X(59,I,J)
     &         +YE(I,J)*((16.d0*G(3)**2-2./5.d0*G(1)**2)*TRACE(X,2)
     &                   +6./5.d0*G(1)**2*TRACE(X,3))
     &         +6.d0*G(2)**2*X(8,I,J)
     &         +YE(I,J)*(15./2.d0*G(2)**4+9./5.d0*G(2)**2*G(1)**2
     &  		 +27./2.d0*G(1)**4)
	      K=K+1
325       continue
	ENDIF
c--------- neutrino sector ------------------------------
        IF(IRHN.GT.0.AND.
     &     Q.GE.MIN(ABS(MRNDEC(1)),ABS(MRNDEC(2)),ABS(MRNDEC(3)))) THEN
          CALL MPROD4X(X,209,YN,YED,YE,YND)
          CALL MPROD4X(X,210,YN,YND,YN,YND)
          CALL MPROD5X(X,206,YN,YED,YE,YED,YE)
          CALL MPROD5X(X,207,YN,YED,YE,YND,YN)
          CALL MPROD5X(X,208,YN,YND,YN,YND,YN)
          CALL MPROD5X(X,211,YE,YND,YN,YED,YE)
          CALL MPROD5X(X,212,YE,YND,YN,YND,YN)
          CALL MPROD5X(X,213,YN,YED,YE,YND,MRHN)
          CALL MPROD5X(X,214,YN,YND,YN,YND,MRHN)

c...contributions to gauge couplings
          F2(1)=F2(1)-6./5.d0*G(1)**3*TRACE(X,200)
          F2(2)=F2(2)-2.d0*G(2)**3*TRACE(X,200)
      
	  K=0
          DO 330 I=1,3
            DO 330 J=1,3
c...neutrino Yukawa matrix
      F2(112+K)=-2.d0*X(206,I,J)-2.d0*X(207,I,J)-4.d0*X(208,I,J)
     &          -X(204,I,J)*(3.d0*TRACE(X,2)+TRACE(X,3))
     &          -X(203,I,J)*(3.d0*TRACE(X,201)+9.d0*TRACE(X,1))
     &          +6./5.d0*G(1)**2*(X(204,I,J)+X(203,I,J))
     &          +6.d0*G(2)**2*X(203,I,J)
     &          -YN(I,J)*(TRACE(X,209)+3.d0*TRACE(X,210)
     &                    +3.d0*TRACE(X,50)+9.d0*TRACE(X,49)
     &                    +(4./5.d0*G(1)**2+16.d0*G(3)**2)*TRACE(X,1)
     &                    +207./50.d0*G(1)**4+9./5.d0*G(1)**2*G(2)**2
     &                    +15./2.d0*G(2)**4)
c...contributions to the other Yukawas
      F2(4+K)=F2(4+K)-3.d0*X(4,I,J)*TRACE(X,201)
     &        -YU(I,J)*(TRACE(X,209)+3.d0*TRACE(X,210))
      F2(13+K)=F2(13+K)-YD(I,J)*TRACE(X,209)-X(7,I,J)*TRACE(X,201)
      F2(22+K)=F2(22+K)-2.d0*X(211,I,J)-2.d0*X(212,I,J)
     &         -X(202,I,J)*(TRACE(X,201)+3.d0*TRACE(X,1))
     &         -YE(I,J)*TRACE(X,209)
	      K=K+1
330       continue        
	  DO I=112,120
	    F(I)=F(I)+F2(I)/FAC**2
	  ENDDO
	ENDIF
	DO I=1,30
	  F(I)=F(I)+F2(I)/FAC**2
	ENDDO
      ENDIF
C
C  2-loop part for SSB terms
C
      IF (ITWOLP.GE.2) THEN
c...evaluate auxiliary matrices
        CALL MPROD2X(X,112,TUD,YU)
        CALL MPROD2X(X,113,TDD,YD)
        CALL MPROD2X(X,114,TED,YE)
        CALL MPROD2X(X,145,TU,YUD)
        CALL MPROD2X(X,146,TD,YDD)
        CALL MPROD2X(X,147,TE,YED)
        CALL MPROD2X(X,148,YU,TUD)
        CALL MPROD2X(X,149,YD,TDD)
        CALL MPROD2X(X,150,YE,TED)
        CALL MPROD4X(X,63,TU,YUD,YU,YUD)
        CALL MPROD4X(X,64,TU,YDD,YD,YUD)
        CALL MPROD4X(X,65,TD,YUD,YU,YDD)
        CALL MPROD4X(X,75,TD,YDD,YD,YDD)
        CALL MPROD4X(X,76,TE,YED,YE,YED)
        CALL MPROD4X(X,92,TUD,TU,YUD,YU)
        CALL MPROD4X(X,93,TUD,YU,YUD,TU)
        CALL MPROD4X(X,94,TDD,TD,YUD,YU)
        CALL MPROD4X(X,95,YDD,YD,TUD,TU)
        CALL MPROD4X(X,96,TDD,YD,YUD,TU)
        CALL MPROD4X(X,97,YDD,TD,TUD,YU)
        CALL MPROD4X(X,102,TDD,TD,YDD,YD)
        CALL MPROD4X(X,103,TDD,YD,YDD,TD)
        CALL MPROD4X(X,104,TED,TE,YED,YE)
        CALL MPROD4X(X,105,TED,YE,YED,TE)
        CALL MPROD4X(X,115,YUD,YU,TUD,TU)
        CALL MPROD4X(X,116,YUD,TU,TUD,YU)
        CALL MPROD4X(X,117,YDD,YD,TDD,TD)
        CALL MPROD4X(X,118,YDD,TD,TDD,YD)
        CALL MPROD4X(X,122,YED,YE,TED,TE)
        CALL MPROD4X(X,123,YED,TE,TED,YE)
        CALL MPROD4X(X,125,YUD,YU,YUD,YU)
        CALL MPROD4X(X,126,YDD,YD,YDD,YD)
        CALL MPROD4X(X,127,YED,YE,YED,YE)
        CALL MPROD4X(X,137,TU,TUD,YU,YUD)
        CALL MPROD4X(X,138,YU,YUD,TU,TUD)
        CALL MPROD4X(X,139,TU,YUD,YU,TUD)
        CALL MPROD4X(X,140,YU,TUD,TU,YUD)
        CALL MPROD4X(X,141,TU,TDD,YD,YUD)
        CALL MPROD4X(X,142,YU,YDD,TD,TUD)
        CALL MPROD4X(X,143,TU,YDD,YD,TUD)
        CALL MPROD4X(X,144,YU,TDD,TD,YUD)
        CALL MPROD4X(X,157,YD,YUD,YU,YDD)
        CALL MPROD4X(X,162,TD,TDD,YD,YDD)
        CALL MPROD4X(X,163,YD,YDD,TD,TDD)
        CALL MPROD4X(X,164,TD,YDD,YD,TDD)
        CALL MPROD4X(X,165,YD,TDD,TD,YDD)
        CALL MPROD4X(X,166,TD,TUD,YU,YDD)
        CALL MPROD4X(X,167,YD,YUD,TU,TDD)
        CALL MPROD4X(X,168,TD,YUD,YU,TDD)
        CALL MPROD4X(X,169,YD,TUD,TU,YDD)
        CALL MPROD4X(X,175,TE,TED,YE,YED)
        CALL MPROD4X(X,176,YE,YED,TE,TED)
        CALL MPROD4X(X,177,TE,YED,YE,TED)
        CALL MPROD4X(X,178,YE,TED,TE,YED)
        CALL MPROD5X(X,60,TU,YUD,YU,YUD,YU)
        CALL MPROD5X(X,61,TU,YDD,YD,YDD,YD)
        CALL MPROD5X(X,62,TU,YDD,YD,YUD,YU)
        CALL MPROD5X(X,66,YU,YUD,YU,YUD,TU)
        CALL MPROD5X(X,67,YU,YUD,TU,YUD,YU)
        CALL MPROD5X(X,68,YU,YDD,YD,YDD,TD)
        CALL MPROD5X(X,69,YU,YDD,TD,YDD,YD)
        CALL MPROD5X(X,70,YU,YDD,YD,YUD,TU)
        CALL MPROD5X(X,71,YU,YDD,TD,YUD,YU)
        CALL MPROD5X(X,72,TD,YDD,YD,YDD,YD)
        CALL MPROD5X(X,73,TD,YUD,YU,YUD,YU)
        CALL MPROD5X(X,74,TD,YUD,YU,YDD,YD)
        CALL MPROD5X(X,77,YD,YDD,YD,YDD,TD)
        CALL MPROD5X(X,78,YD,YDD,TD,YDD,YD)
        CALL MPROD5X(X,79,YD,YUD,TU,YUD,YU)
        CALL MPROD5X(X,80,YD,YUD,YU,YUD,TU)
        CALL MPROD5X(X,81,YD,YUD,TU,YDD,YD)
        CALL MPROD5X(X,82,YD,YUD,YU,YDD,TD)
        CALL MPROD5X(X,83,TE,YED,YE,YED,YE)
        CALL MPROD5X(X,84,YE,YED,YE,YED,TE)
        CALL MPROD5X(X,85,YE,YED,TE,YED,YE)
        CALL MPROD5X(X,86,M2Q,YUD,YU,YUD,YU)
        CALL MPROD5X(X,87,YUD,M2U,YU,YUD,YU)
        CALL MPROD5X(X,88,M2Q,YUD,YU,YDD,YD)
        CALL MPROD5X(X,89,YUD,M2U,YU,YDD,YD)
        CALL MPROD5X(X,90,YUD,YU,M2Q,YDD,YD)
        CALL MPROD5X(X,91,YUD,YU,YDD,M2D,YD)
        CALL MPROD5X(X,98,M2Q,YDD,YD,YDD,YD)
        CALL MPROD5X(X,99,YDD,M2D,YD,YDD,YD)
        CALL MPROD5X(X,100,M2L,YED,YE,YED,YE)
        CALL MPROD5X(X,101,YED,M2E,YE,YED,YE)
        CALL MPROD5X(X,106,YUD,YU,M2Q,YUD,YU)
        CALL MPROD5X(X,107,YUD,YU,YUD,M2U,YU)
        CALL MPROD5X(X,108,YUD,YU,YUD,YU,M2Q)
        CALL MPROD5X(X,109,YDD,YD,M2Q,YDD,YD)
        CALL MPROD5X(X,110,YDD,YD,YDD,M2D,YD)
        CALL MPROD5X(X,111,YDD,YD,YDD,YD,M2Q)
        CALL MPROD5X(X,119,YED,YE,M2L,YED,YE)
        CALL MPROD5X(X,120,YED,YE,YED,M2E,YE)
        CALL MPROD5X(X,121,YED,YE,YED,YE,M2L)
        CALL MPROD5X(X,124,M2U,YU,YUD,YU,YUD)
        CALL MPROD5X(X,128,YU,M2Q,YUD,YU,YUD)
        CALL MPROD5X(X,129,YU,YUD,M2U,YU,YUD)
        CALL MPROD5X(X,130,YU,YUD,YU,M2Q,YUD)
        CALL MPROD5X(X,131,YU,YUD,YU,YUD,M2U)
        CALL MPROD5X(X,132,M2U,YU,YDD,YD,YUD)
        CALL MPROD5X(X,133,YU,M2Q,YDD,YD,YUD)
        CALL MPROD5X(X,134,YU,YDD,M2D,YD,YUD)
        CALL MPROD5X(X,135,YU,YDD,YD,M2Q,YUD)
        CALL MPROD5X(X,136,YU,YDD,YD,YUD,M2U)
        CALL MPROD5X(X,151,M2D,YD,YDD,YD,YDD)
        CALL MPROD5X(X,152,YD,M2Q,YDD,YD,YDD)
        CALL MPROD5X(X,153,YD,YDD,M2D,YD,YDD)
        CALL MPROD5X(X,154,YD,YDD,YD,M2Q,YDD)
        CALL MPROD5X(X,155,YD,YDD,YD,YDD,M2D)
        CALL MPROD5X(X,156,M2D,YD,YUD,YU,YDD)
        CALL MPROD5X(X,158,YD,M2Q,YUD,YU,YDD)
        CALL MPROD5X(X,159,YD,YUD,M2U,YU,YDD)
        CALL MPROD5X(X,160,YD,YUD,YU,M2Q,YDD)
        CALL MPROD5X(X,161,YD,YUD,YU,YDD,M2D)
        CALL MPROD5X(X,170,M2E,YE,YED,YE,YED)
        CALL MPROD5X(X,171,YE,M2L,YED,YE,YED)
        CALL MPROD5X(X,172,YE,YED,M2E,YE,YED)
        CALL MPROD5X(X,173,YE,YED,YE,M2L,YED)
        CALL MPROD5X(X,174,YE,YED,YE,YED,M2E)
c...compute some convenient quantities
	SP=-(3.d0*G(63)*TRACE(X,1)+TRACE(X,31))+4.d0*TRACE(X,34)
     &     +3.d0*G(64)*TRACE(X,2)-TRACE(X,25)-2.d0*TRACE(X,26)
     &     +G(64)*TRACE(X,3)+TRACE(X,27)-2.d0*TRACE(X,28)
     &     +(3./2.d0*G(2)**2+3./10.d0*G(1)**2)
     &      *(G(63)-G(64)-TR3X3(M2L))
     &     +(8./3.d0*G(3)**2+3./2.d0*G(2)**2+G(1)**2/30.d0)*TR3X3(M2Q)
     &     -(16./3.d0*G(3)**2+16./15.d0*G(1)**2)*TR3X3(M2U)
     &     +(8./3.d0*G(3)**2+2./15.d0*G(1)**2)*TR3X3(M2D)
     &     +6./5.d0*G(1)**2*TR3X3(M2E)
        SIG1=G(1)**2/5.d0*(3.d0*(G(63)+G(64))+TR3X3(M2Q)
     &			   +3.d0*TR3X3(M2L)+8.d0*TR3X3(M2U)
     &			   +2.d0*TR3X3(M2D)+6.d0*TR3X3(M2E))
	SIG2=G(2)**2*(G(63)+G(64)+3.d0*TR3X3(M2Q)+TR3X3(M2L))
	SIG3=G(3)**2*(2.d0*TR3X3(M2Q)+TR3X3(M2U)+TR3X3(M2D))
c...gaugino masses      
      F2(31)=2.d0*G(1)**2
     &	     *(199./25.d0*G(1)**2*(G(31)+G(31))
     &	       +27./5.d0*G(2)**2*(G(31)+G(32))
     &	       +88./5.d0*G(3)**2*(G(31)+G(33))
     &	       +26./5.d0*(TRACE(X,11)-G(31)*TRACE(X,1))
     &	       +14./5.d0*(TRACE(X,16)-G(31)*TRACE(X,2))
     &	       +18./5.d0*(TRACE(X,17)-G(31)*TRACE(X,3)))
      F2(32)=2.d0*G(2)**2
     &	     *(9./5.d0*G(1)**2*(G(32)+G(31))
     &	       +25.d0*G(2)**2*(G(32)+G(32))
     &	       +24.d0*G(3)**2*(G(32)+G(33))
     &	       +6.d0*(TRACE(X,11)-G(32)*TRACE(X,1))
     &	       +6.d0*(TRACE(X,16)-G(32)*TRACE(X,2))
     &	       +2.d0*(TRACE(X,17)-G(32)*TRACE(X,3)))
      F2(33)=2.d0*G(3)**2
     &	     *(11./5.d0*G(1)**2*(G(33)+G(31))
     &	       +9.d0*G(2)**2*(G(33)+G(32))
     &	       +14.d0*G(3)**2*(G(33)+G(33))
     &	       +4.d0*(TRACE(X,11)-G(33)*TRACE(X,1))
     &	       +4.d0*(TRACE(X,16)-G(33)*TRACE(X,2)))
c...higgs mixing parameters      
      F2(61)=G(61)*(-3.*(3.d0*TRACE(X,49)+3.d0*TRACE(X,54)
     &			 +2.d0*TRACE(X,50)+TRACE(X,55))
     &  	    +(16.d0*G(3)**2+4./5.d0*G(1)**2)*TRACE(X,1)
     &  	    +(16.d0*G(3)**2-2./5.d0*G(1)**2)*TRACE(X,2)
     &  	    +6./5.d0*G(1)**2*TRACE(X,3)
     &  	    +15./2.d0*G(2)**4+9./5.d0*G(1)**2*G(2)**2
     &  	    +207./50.d0*G(1)**4)
      F2(62)=G(62)*(-3.*(3.d0*TRACE(X,49)+3.d0*TRACE(X,54)
     &	        	 +2.d0*TRACE(X,50)+TRACE(X,55))
     &              +(16.d0*G(3)**2+4./5.d0*G(1)**2)*TRACE(X,1)
     &              +(16.d0*G(3)**2-2./5.d0*G(1)**2)*TRACE(X,2)
     &              +6./5.d0*G(1)**2*TRACE(X,3)
     &              +15./2.d0*G(2)**4+9./5.d0*G(1)**2*G(2)**2
     &              +207./50.d0*G(1)**4)
     &      +G(61)*(-12.*(3.d0*TRACE(X,63)+3.d0*TRACE(X,75)
     &      		  +TRACE(X,64)+TRACE(X,65)+TRACE(X,76))
     &      	    +(32.d0*G(3)**2+8./5.d0*G(1)**2)*TRACE(X,11)
     &      	    +(32.d0*G(3)**2-4./5.d0*G(1)**2)*TRACE(X,16)
     &      	    +12./5.d0*G(1)**2*TRACE(X,17)
     &      	    -(32.d0*G(3)**2*G(33)+8./5.d0*G(1)**2*G(31))
     &      	     *TRACE(X,1)
     &      	    -(32.d0*G(3)**2*G(33)-4./5.d0*G(1)**2*G(31))
     &      	     *TRACE(X,2)
     &      	    -12./5.d0*G(1)**2*G(31)*TRACE(X,3)
     &      	    -30.d0*G(2)**4*G(32)
     &      	    -18./5.d0*G(1)**2*G(2)**2*(G(31)+G(32))
     &      	    -414./25.d0*G(1)**4*G(31))
c...higgs mass^2
      F2(63)=-6.*(6.*(G(63)*TRACE(X,49)+TRACE(X,86))+6.*TRACE(X,87)
     &       	  +(G(63)+G(64))*TRACE(X,50)+TRACE(X,88)+TRACE(X,89)
     &       	  +TRACE(X,90)+TRACE(X,91)+6.*TRACE(X,92)
     &       	  +6.*TRACE(X,93)+TRACE(X,94)+TRACE(X,95)+TRACE(X,96)
     &       	  +TRACE(X,97))
     &       +(32.*G(3)**2+8./5.d0*G(1)**2)
     &        *(G(63)*TRACE(X,1)+TRACE(X,31)+TRACE(X,34)+TRACE(X,24))
     &       +32.*G(3)**2*(2.*G(33)**2*TRACE(X,1)-G(33)*TRACE(X,11)
     &       		   -G(33)*TRACE(X,112))
     &       +8./5.d0*G(1)**2*(2.*G(31)**2*TRACE(X,1)-G(31)*TRACE(X,11)
     &       		       -G(31)*TRACE(X,112))
     &       +6./5.d0*G(1)**2*SP
     &       +33.*G(2)**4*G(32)**2
     &       +18./5.d0*G(2)**2*G(1)**2*(G(32)**2+G(31)**2+G(31)*G(32))
     &       +621./25.d0*G(1)**4*G(31)**2
     &       +3.*G(2)**2*SIG2+3./5.d0*G(1)**2*SIG1
      F2(64)=-6.*(6.*(G(64)*TRACE(X,54)+TRACE(X,98))+6.*TRACE(X,99)
     &       	  +(G(63)+G(64))*TRACE(X,50)+TRACE(X,88)+TRACE(X,89)
     &       	  +TRACE(X,90)+TRACE(X,91)+2.*G(64)*TRACE(X,55)
     &       	  +2.*TRACE(X,100)+2.*TRACE(X,101)+6.*TRACE(X,102)
     &       	  +6.*TRACE(X,103)+TRACE(X,95)+TRACE(X,94)+TRACE(X,97)
     &       	  +TRACE(X,96)+2.*TRACE(X,104)+2.*TRACE(X,105))
     &       +(32.*G(3)**2-4./5.d0*G(1)**2)
     &        *(G(64)*TRACE(X,2)+TRACE(X,25)+TRACE(X,26)+TRACE(X,29))
     &       +32.*G(3)**2*(2.*G(33)**2*TRACE(X,2)-G(33)*TRACE(X,16)
     &       		   -G(33)*TRACE(X,113))
     &       -4./5.d0*G(1)**2*(2.*G(31)**2*TRACE(X,2)-G(31)*TRACE(X,16)
     &       		       -G(31)*TRACE(X,113))
     &       +12./5.d0*G(1)**2*(G(64)*TRACE(X,3)+TRACE(X,27)
     &                          +TRACE(X,28)+TRACE(X,45)
     &                          +2.*G(31)**2*TRACE(X,3)
     &                          -G(31)*TRACE(X,114)-G(31)*TRACE(X,17))
     &       -6./5.d0*G(1)**2*SP
     &       +33.*G(2)**4*G(32)**2
     &       +18./5.d0*G(2)**2*G(1)**2*(G(32)**2+G(31)**2+G(31)*G(32))
     &       +621./25.d0*G(1)**4*G(31)**2
     &       +3.*G(2)**2*SIG2+3./5.d0*G(1)**2*SIG1

        K=0
	DO 400 I=1,3
          DO 400 J=1,3
C...soft trilinear scalar couplings
      F2(34+K)=TU(I,J)*(-3.)*(3.*TRACE(X,49)+TRACE(X,50))
     &         -X(10,I,J)*(3.*TRACE(X,2)+TRACE(X,3))
     &         -15.d0*X(9,I,J)*TRACE(X,1)
     &         -6.d0*X(60,I,J)-2.d0*X(61,I,J)-4.d0*X(62,I,J)
     &         +TU(I,J)*(16.*G(3)**2+4./5.d0*G(1)**2)*TRACE(X,1)
     &         +12.d0*G(2)**2*X(9,I,J)
     &         +2./5.d0*G(1)**2*X(10,I,J)
     &         +TU(I,J)*(-16./9.d0*G(3)**4+8.*G(3)**2*G(2)**2
     &                   +136./45.d0*G(3)**2*G(1)**2+15./2.d0*G(2)**4
     &                   +G(2)**2*G(1)**2+2743./450.d0*G(1)**4)
     &         +YU(I,J)*(-6.)*(6.*TRACE(X,63)+TRACE(X,64)+TRACE(X,65))
     &         -18.d0*X(4,I,J)*TRACE(X,11)
     &         -X(5,I,J)*(6.d0*TRACE(X,16)+2.d0*TRACE(X,17))
     &         -12.d0*X(12,I,J)*TRACE(X,1)
     &         -X(13,I,J)*(6.d0*TRACE(X,2)+2.d0*TRACE(X,3))
     &         -6.*X(66,I,J)-8.*X(67,I,J)-4.*X(68,I,J)-4.*X(69,I,J)
     &         -2.*X(70,I,J)-4.*X(71,I,J)
     &         +YU(I,J)*(32.*G(3)**2+8./5.d0*G(1)**2)*TRACE(X,11)
     &         +(6.*G(2)**2+6./5.d0*G(1)**2)*X(12,I,J)
     &         +4./5.d0*G(1)**2*X(13,I,J)
     &         -YU(I,J)*(32.*G(3)**2*G(33)+8./5.d0*G(1)**2*G(31))
     &          *TRACE(X,1)
     &         -(12.*G(2)**2*G(32)+4./5.d0*G(1)**2*G(31))*X(4,I,J)
     &         -4./5.d0*G(1)**2*G(31)*X(5,I,J)
     &         +YU(I,J)*(64./9.d0*G(3)**4*G(33)
     &                   -16.*G(3)**2*G(2)**2*(G(33)+G(32))
     &                   -272./45.d0*G(3)**2*G(1)**2*(G(33)+G(31))
     &                   -30.*G(2)**4*G(32)
     &                   -2.*G(2)**2*G(1)**2*(G(32)+G(31))
     &                   -5486./225.d0*G(1)**4*G(31))
      F2(43+K)=TD(I,J)*(-3.)*(3.*TRACE(X,54)+TRACE(X,50)+TRACE(X,55))
     &         -3.*X(15,I,J)*TRACE(X,1)
     &         -5.*X(14,I,J)*(3.*TRACE(X,2)+TRACE(X,3))
     &         -6.*X(72,I,J)-2.*X(73,I,J)-4.*X(74,I,J)
     &         +TD(I,J)*((16.*G(3)**2-2./5.d0*G(1)**2)*TRACE(X,2)
     &                   +6./5.d0*G(1)**2*TRACE(X,3))
     &         +4./5.d0*G(1)**2*X(15,I,J)
     &         +(12.*G(2)**2+6./5.d0*G(1)**2)*X(14,I,J)
     &         +TD(I,J)*(-16./9.d0*G(3)**4+8.*G(3)**2*G(2)**2
     &                   +8./9.d0*G(3)**2*G(1)**2+15./2.d0*G(2)**4
     &                   +G(2)**2*G(1)**2+287./90.d0*G(1)**4)
     &         +YD(I,J)*(-6.)*(6.*TRACE(X,75)+TRACE(X,64)+TRACE(X,65)
     &                         +2.*TRACE(X,76))
     &         -6.*X(7,I,J)*TRACE(X,11)
     &         -6.*X(6,I,J)*(3.*TRACE(X,16)+TRACE(X,17))
     &         -6.*X(19,I,J)*TRACE(X,1)
     &         -4.*X(18,I,J)*(3.*TRACE(X,2)+TRACE(X,3))
     &         -6.*X(77,I,J)-8.*X(78,I,J)-4.*X(79,I,J)-4.*X(80,I,J)
     &         -4.*X(81,I,J)-2.*X(82,I,J)
     &         +YD(I,J)*((32.*G(3)**2-4./5.d0*G(1)**2)*TRACE(X,16)
     &                   +12./5.d0*G(1)**2*TRACE(X,17))
     &         +8./5.d0*G(1)**2*X(19,I,J)
     &         +(6.*G(2)**2+6./5.d0*G(1)**2)*X(18,I,J)
     &         -YD(I,J)*((32.*G(3)**2*G(33)-4./5.d0*G(1)**2*G(31))
     &                    *TRACE(X,2)
     &                   +12./5.d0*G(1)**2*G(31)*TRACE(X,3))
     &         -(12.*G(2)**2*G(32)+8./5.d0*G(1)**2*G(31))*X(6,I,J)
     &         -8./5.d0*G(1)**2*G(31)*X(7,I,J)
     &         +YD(I,J)*(64./9.d0*G(3)**4*G(33)
     &                   -16.*G(3)**2*G(2)**2*(G(33)+G(32))
     &                   -16./9.d0*G(3)**2*G(1)**2*(G(33)+G(31))
     &                   -30.*G(2)**4*G(32)
     &                   -2.*G(2)**2*G(1)**2*(G(32)+G(31))
     &                   -574./45.d0*G(1)**4*G(31))
      F2(52+K)=TE(I,J)*(-3.)*(3.*TRACE(X,54)+TRACE(X,50)+TRACE(X,55))
     &         -5.*X(20,I,J)*(3.*TRACE(X,2)+TRACE(X,3))
     &         -6.*X(83,I,J)
     &         +TE(I,J)*((16.*G(3)**2-2./5.d0*G(1)**2)*TRACE(X,2)
     &                   +6./5.d0*G(1)**2*TRACE(X,3))
     &         +(12.*G(2)**2-6./5.d0*G(1)**2)*X(20,I,J)
     &         +TE(I,J)*(15./2.d0*G(2)**4+9./5.d0*G(2)**2*G(1)**2
     &                   +27./2.d0*G(1)**4)
     &         +YE(I,J)*(-6.)*(6.*TRACE(X,75)+TRACE(X,64)+TRACE(X,65)
     &                         +2.*TRACE(X,76))
     &         -4.*X(21,I,J)*(3.*TRACE(X,2)+TRACE(X,3))
     &         -6.*X(8,I,J)*(3.*TRACE(X,16)+TRACE(X,17))
     &         -6.*X(84,I,J)-8.*X(85,I,J)
     &         +YE(I,J)*((32.*G(3)**2-4./5.d0*G(1)**2)*TRACE(X,16)
     &                   +12./5.d0*G(1)**2*TRACE(X,17))
     &         +(6.*G(2)**2+6./5.d0*G(1)**2)*X(21,I,J)
     &         -YE(I,J)*((32.*G(3)**2*G(33)-4./5.d0*G(1)**2*G(31))
     &                    *TRACE(X,2)
     &                   +12./5.d0*G(1)**2*G(31)*TRACE(X,3))
     &         -12.*G(2)**2*G(32)*X(8,I,J)
     &         +YE(I,J)*(-30.*G(2)**4*G(32)
     &                   -18./5.d0*G(2)**2*G(1)**2*(G(31)+G(32))
     &                   -54.*G(1)**4*G(31))
c...Squark doublet mass^2
      F2(65+K)=-(2.*X(86,I,J)+8.*G(63)*X(125,I,J))-4.*X(87,I,J)
     &         -4.*X(106,I,J)-4.*X(107,I,J)-2.*X(108,I,J)
     &         -(2.*X(98,I,J)+8.*G(64)*X(126,I,J))-4.*X(99,I,J)
     &         -4.*X(109,I,J)-4.*X(110,I,J)-2.*X(111,I,J)
     &         -(X(22,I,J)+4.*G(63)*X(1,I,J)+2.*X(23,I,J)+X(31,I,J))
     &          *3.*TRACE(X,1)
     &         -(X(25,I,J)+4.*G(64)*X(2,I,J)+2.*X(26,I,J)+X(32,I,J))
     &          *(3.*TRACE(X,2)+TRACE(X,3))
     &         -6.*X(1,I,J)*(TRACE(X,22)+TRACE(X,23))
     &         -X(2,I,J)*(6.*TRACE(X,25)+6.*TRACE(X,26)+2.*TRACE(X,27)
     &          	  +2.*TRACE(X,28))
     &         -4.*(X(115,I,J)+X(92,I,J)+X(116,I,J)+X(93,I,J))
     &         -4.*(X(117,I,J)+X(102,I,J)+X(118,I,J)+X(103,I,J))
     &         -X(24,I,J)*6.*TRACE(X,1)-X(1,I,J)*6.*TRACE(X,24)
     &         -X(112,I,J)*6.*TRACE(X,11)-X(11,I,J)*6.*TRACE(X,112)
     &         -X(29,I,J)*(6.*TRACE(X,2)+2.*TRACE(X,3))
     &         -X(2,I,J)*(6.*TRACE(X,29)+2.*TRACE(X,30))
     &         -X(113,I,J)*(6.*TRACE(X,16)+2.*TRACE(X,17))
     &         -X(16,I,J)*(6.*TRACE(X,113)+2.*TRACE(X,114))
     &         +2./5.d0*G(1)**2*(2.*X(22,I,J)+4.*G(63)*X(1,I,J)
     &          	        +4.*X(23,I,J)+2.*X(31,I,J)+4.*X(24,I,J)
     &          	        -4.*G(31)*X(112,I,J)-4.*G(31)*X(11,I,J)
     &          	        +8.*G(31)**2*X(1,I,J)
     &                          +X(25,I,J)+2.*G(64)*X(2,I,J)
     &                          +2.*X(26,I,J)+X(32,I,J)+2.*X(29,I,J)
     &          	        -2.*G(31)*X(113,I,J)-2.*G(31)*X(16,I,J)
     &          	        +4.*G(31)**2*X(2,I,J))
     &         +I3(I,J)*(2./5.d0*G(1)**2*SP
     &          	 -128./3.*G(3)**4*G(33)**2
     &          	 +32.*G(3)**2*G(2)**2
     &          	  *(G(33)**2+G(32)**2+G(32)*G(33))
     &          	 +32./45.d0*G(3)**2*G(1)**2
     &          	  *(G(33)**2+G(31)**2+G(31)*G(33))
     &          	 +33.*G(2)**4*G(32)**2
     &          	 +2./5.d0*G(2)**2*G(1)**2
     &          	  *(G(32)**2+G(31)**2+G(31)*G(32))
     &          	 +199./75.d0*G(1)**4*G(31)**2
     &          	 +16./3.d0*G(3)**2*SIG3+3.*G(2)**2*SIG2
     &          	 +1./15.d0*G(1)**2*SIG1)
c...Slepton doublet mass^2
      F2(74+K)=-(2.*X(100,I,J)+8.*G(64)*X(127,I,J))-4.*X(101,I,J)
     &         -4.*X(119,I,J)-4.*X(120,I,J)-2.*X(121,I,J)
     &         -(X(27,I,J)+4.*G(64)*X(3,I,J)+2.*X(28,I,J)+X(33,I,J))
     &          *(3.*TRACE(X,2)+TRACE(X,3))
     &         -X(3,I,J)*(6.*X(25,I,J)+6.*X(26,I,J)+2.*X(27,I,J)
     &          	  +2.*X(28,I,J))
     &         -4.*(X(122,I,J)+X(104,I,J)+X(123,I,J)+X(105,I,J))
     &         -X(30,I,J)*(6.*TRACE(X,2)+2.*TRACE(X,3))
     &         -X(3,I,J)*(6.*TRACE(X,29)+2.*TRACE(X,30))
     &         -X(114,I,J)*(6.*TRACE(X,16)+2.*TRACE(X,17))
     &         -X(17,I,J)*(6.*TRACE(X,113)+2.*TRACE(X,114))
     &         +6./5.d0*G(1)**2*(X(27,I,J)+2.*G(64)*X(3,I,J)
     &          	        +2.*X(28,I,J)+X(33,I,J)+2.*X(30,I,J)
     &          	        -2.*G(31)*X(114,I,J)-2.*G(31)*X(17,I,J)
     &          	        +4.*G(31)**2*X(3,I,J))
     &         +I3(I,J)*(-6./5.d0*G(1)**2*SP
     &          	 +33.*G(2)**4*G(32)**2
     &          	 +18./5.d0*G(2)**2*G(1)**2
     &          	  *(G(32)**2+G(31)**2+G(31)*G(32))
     &          	 +621./25.d0*G(1)**4*G(31)**2
     &          	 +3.*G(2)**2*SIG2+3./5.d0*G(1)**2*SIG1)
c...Up squark singlet mass^2
      F2(83+K)=-(2.*X(124,I,J)+8.*G(63)*X(49,I,J))-4.*X(128,I,J)
     &         -4.*X(129,I,J)-4.*X(130,I,J)-2.*X(131,I,J)
     &         -(2.*X(132,I,J)+4.*(G(63)+G(64))*X(50,I,J))
     &         -4.*X(133,I,J)-4.*X(134,I,J)-4.*X(135,I,J)
     &         -2.*X(136,I,J)
     &         -(X(34,I,J)+4.*G(63)*X(46,I,J)+2.*X(35,I,J)+X(36,I,J))
     &          *6.*TRACE(X,1)
     &         -12.*X(46,I,J)*(TRACE(X,22)+TRACE(X,23))
     &         -4.*(X(137,I,J)+X(138,I,J)+X(139,I,J)+X(140,I,J))
     &         -4.*(X(141,I,J)+X(142,I,J)+X(143,I,J)+X(144,I,J))
     &         -12.*(X(37,I,J)*TRACE(X,1)+X(46,I,J)*TRACE(X,24)
     &               +X(145,I,J)*TRACE(X,112)+X(148,I,J)*TRACE(X,11))
     &         +(6.*G(2)**2-2./5.d0*G(1)**2)
     &          *(X(34,I,J)+2.*G(63)*X(46,I,J)+2.*X(35,I,J)+X(36,I,J)
     &            +2.*X(37,I,J))
     &         +12.*G(2)**2*(2.*G(32)**2*X(46,I,J)-G(32)*X(145,I,J)
     &          	     -G(32)*X(148,I,J))
     &         -4./5.d0*G(1)**2*(2.*G(31)**2*X(46,I,J)-G(31)*X(145,I,J)
     &          	         -G(31)*X(148,I,J))
     &         +I3(I,J)*(-8./5.d0*G(1)**2*SP
     &          	 -128./3.d0*G(3)**4*G(33)**2
     &          	 +512./45.d0*G(3)**2*G(1)**2
     &          	  *(G(33)**2+G(31)**2+G(31)*G(33))
     &          	 +3424./75.d0*G(1)**4*G(31)**2
     &          	 +16./3.d0*G(3)**2*SIG3
     &          	 +16./15.d0*G(1)**2*SIG1)
c...Down squark singlet mass^2
      F2(92+K)=-(2.*X(151,I,J)+8.*G(64)*X(54,I,J))-4.*X(152,I,J)
     &         -4.*X(153,I,J)-4.*X(154,I,J)-2.*X(155,I,J)
     &         -(2.*X(156,I,J)+4.*(G(63)+G(64))*X(157,I,J))
     &         -4.*X(158,I,J)-4.*X(159,I,J)-4.*X(160,I,J)
     &         -2.*X(161,I,J)
     &         -(X(38,I,J)+4.*G(64)*X(47,I,J)+2.*X(39,I,J)+X(40,I,J))
     &          *(6.*TRACE(X,2)+2.*TRACE(X,3))
     &         -4.*X(47,I,J)*(3.*TRACE(X,25)+3.*TRACE(X,26)
     &          	      +TRACE(X,27)+TRACE(X,28))
     &         -4.*(X(162,I,J)+X(163,I,J)+X(164,I,J)+X(165,I,J))
     &         -4.*(X(166,I,J)+X(167,I,J)+X(168,I,J)+X(169,I,J))
     &         -4.*X(41,I,J)*(3.*TRACE(X,2)+TRACE(X,3))
     &         -4.*X(47,I,J)*(3.*TRACE(X,29)+TRACE(X,30))
     &         -4.*X(146,I,J)*(3.*TRACE(X,113)+TRACE(X,114))
     &         -4.*X(149,I,J)*(3.*TRACE(X,16)+TRACE(X,17))
     &         +(6.*G(2)**2+2./5.d0*G(1)**2)
     &          *(X(38,I,J)+2.*G(64)*X(47,I,J)+2.*X(39,I,J)
     &            +X(40,I,J)+2.*X(41,I,J))
     &         +12.*G(2)**2*(2.*G(32)**2*X(47,I,J)-G(32)*X(146,I,J)
     &          	     -G(32)*X(149,I,J))
     &         +4./5.d0*G(1)**2*(2.*G(31)**2*X(47,I,J)-G(31)*X(146,I,J)
     &          	         -G(31)*X(149,I,J))
     &         +I3(I,J)*(4./5.d0*G(1)**2*SP
     &          	 -128./3.d0*G(3)**4*G(33)**2
     &          	 +128./45.d0*G(3)**2*G(1)**2
     &          	  *(G(33)**2+G(31)**2+G(31)*G(33))
     &          	 +808./75.d0*G(1)**4*G(31)**2
     &          	 +16./3.d0*G(3)**2*SIG3
     &          	 +4./15.d0*G(1)**2*SIG1)
c...Charged slepton singlet mass^2
      F2(101+K)=-(2.*X(170,I,J)+8.*G(64)*X(55,I,J))-4.*X(171,I,J)
     &          -4.*X(172,I,J)-4.*X(173,I,J)-2.*X(174,I,J)
     &          -(X(42,I,J)+4.*G(64)*X(48,I,J)+2.*X(43,I,J)+X(44,I,J))
     &           *(6.*TRACE(X,2)+2.*TRACE(X,3))
     &          -4.*X(48,I,J)*(3.*TRACE(X,25)+3.*TRACE(X,26)
     &          	       +TRACE(X,27)+TRACE(X,44))
     &          -4.*(X(175,I,J)+X(176,I,J)+X(177,I,J)+X(178,I,J))
     &          -4.*X(45,I,J)*(2.*TRACE(X,2)+TRACE(X,3))
     &          -4.*X(48,I,J)*(3.*TRACE(X,29)+TRACE(X,30))
     &          -4.*X(147,I,J)*(3.*TRACE(X,113)+TRACE(X,114))
     &          -4.*X(150,I,J)*(3.*TRACE(X,16)+TRACE(X,17))
     &          +(6.*G(2)**2-6./5.d0*G(1)**2)
     &           *(X(42,I,J)+2.*G(64)*X(48,I,J)+2.*X(43,I,J)+X(44,I,J)
     &             +2.*X(45,I,J))
     &          +12.*G(2)**2*(2.*G(32)**2*X(48,I,J)-G(32)*X(147,I,J)
     &          	      -G(32)*X(150,I,J))
     &          -12./5.d0*G(1)**2*(2.*G(31)**2*X(48,I,J)
     &          		   -G(31)*X(147,I,J)-G(31)*X(150,I,J))
     &          +I3(I,J)*(12./5.d0*G(1)**2*SP
     &          	  +2808./25.d0*G(1)**4*G(31)**2
     &          	  +12./5.d0*G(1)**2*SIG1)
	    K=K+1
400     CONTINUE
c...VEVs
        F2(110)=0.d0
	F2(111)=0.d0
c--------- neutrino sector ------------------------------
        IF (IRHN.GT.0) THEN
          IF(
     &Q.GE.MIN(ABS(MRNDEC(1)),ABS(MRNDEC(2)),ABS(MRNDEC(3)))) THEN
	    K=0
            DO I=1,3
              DO J=1,3
c...contribution to gaugino masses
            F2(31)=F2(31)+2.d0*G(1)**2 
     &                    *6./5.d0*(TRACE(X,215)-G(31)*TRACE(X,200))
            F2(32)=F2(32)+2.d0*G(2)**2
     &                    *2.d0*(TRACE(X,215)-G(32)*TRACE(X,200))
c...Majorana mass matrix
      F2(121+K)=-2.d0*X(213,J,I)-2.d0*X(214,J,I)
     &          -2.d0*X(213,I,J)-2.d0*X(214,I,J)
     &          -(X(205,J,I)+X(205,I,J))
     &           *(6.d0*TRACE(X,1)+2.d0*TRACE(X,201)
     &             -6./5.d0*G(1)**2-6.d0*G(2)**2)
	        K=K+1
             ENDDO
	    ENDDO
	    DO I=121,129
	      F(I)=F(I)+F2(I)/FAC**2
	    ENDDO
	  ENDIF
c...effective dim-5 neutrino operator
          IF (LND5ON.AND.Q.GE.DBLE(MSUSY).AND.
     &   Q.LE.MAX(ABS(MRNDEC(1)),ABS(MRNDEC(2)),ABS(MRNDEC(3)))) THEN
            CALL MPROD5X(X,237,KA,YED,YE,YED,YE)
            CALL MPROD5X(X,238,KA,YND,YN,YND,YN)
	    K=0
            DO I=1,3
              DO J=1,3
      F2(148+K)=KA(I,J)*(-6.d0*TRACE(X,50)-18.d0*TRACE(X,49)
     &                   -2.d0*TRACE(X,209)-6.d0*TRACE(X,210)
     &                   +(8./5.d0*G(1)**2+32.d0*G(3)**2)*TRACE(X,1)
     &                   +207./25.d0*G(1)**4+18./5.d0*G(1)**2*G(2)**2
     &                   +15.d0*G(2)**4)
     &           -2.d0*X(237,J,I)-2.d0*X(238,J,I)
     &           -2.d0*X(237,I,J)-2.d0*X(238,I,J)
     &           -(X(232,J,I)+X(232,I,J))*(TRACE(X,201)+3.*TRACE(X,1))
     &           -(X(231,J,I)+X(231,I,J))
     &            *(-6./5.d0*G(1)**2+TRACE(X,3)+3.*TRACE(X,2))            
	        K=K+1
              ENDDO
	    ENDDO
	    DO I=148,156
	      F(I)=F(I)+F2(I)/FAC**2
	    ENDDO
	  ENDIF
	ENDIF
	DO I=31,111
	  F(I)=F(I)+F2(I)/FAC**2
	ENDDO
      ENDIF
      
      RETURN
      END
