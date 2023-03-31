C---------------------------------------------------------------------
      SUBROUTINE BSGGUT(YUGUT,YDGUT,YEGUT,YNGUT,G,IMODEL,M0,MHF,A0)            
C----------------------------------------------------------------------
C
C     Sets GUT-scale boundary conditions for G(157).
C
c    Created: 03/13/07 by Azar Mustafayev
c    Modified: 6/12/07 by Azar Mustafayev - compatible with ISAJET 7.75
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IMPLICIT NONE
      REAL*8 G(157),YUGUT(3,3),YDGUT(3,3),YEGUT(3,3),YNGUT(3,3)
      INTEGER IMODEL
      REAL M0,MHF,A0
c
      COMMON /BSGDEC/MSQDEC(3),MSLDEC(3),MSUDEC(3),MSDDEC(3),
     &               MSEDEC(3),MRNDEC(3),IRHN
      REAL*8 MSQDEC,MSLDEC,MSUDEC,MSDDEC,MSEDEC,MRNDEC
      INTEGER IRHN
      SAVE /BSGDEC/
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
c     
      REAL*8 DDILOG
      REAL*8 A(3),XLAMGM,XMESGM,XN5GM,XLM,THRF,THRG,BLHAT,BBHAT,BTHAT,
     &       XC,DY,DG(31),PI
      INTEGER I,J,K

      PI=4.d0*DATAN(1.d0)

C...Compute gauge mediated threshold functions
      IF (IMODEL.EQ.2) THEN
        XLAMGM=M0
        XMESGM=MHF
        XN5GM=A0
        XLM=XLAMGM/XMESGM
        THRF=((1.D0+XLM)*(LOG(1.D0+XLM)-2*DDILOG(XLM/(1.D0+XLM))+
     ,        .5*DDILOG(2*XLM/(1.D0+XLM)))+
     ,       (1.D0-XLM)*(LOG(1.D0-XLM)-2*DDILOG(-XLM/(1.D0-XLM))+
     ,        .5*DDILOG(-2*XLM/(1.D0-XLM))))/XLM**2
        THRG=((1.D0+XLM)*LOG(1.D0+XLM)+(1.D0-XLM)*LOG(1.D0-XLM))/XLM**2
        XC=2*THRF*XLAMGM**2
        DY=SQRT(3./5.)*G(1)*XGMIN(11)
      ENDIF

c...gauge couplings
      G(1)=G(1)
      G(2)=G(2)
      IF (IMODEL.EQ.1.AND.IAL3UN.NE.0) THEN
        G(3)=(G(1)+G(2))/2.d0
      ELSE
        G(3)=G(3)
      ENDIF
c...Yukawa couplings      
      CALL MAT2VEC(G,4,YUGUT,-1)
      CALL MAT2VEC(G,13,YDGUT,-1)
      CALL MAT2VEC(G,22,YEGUT,-1)
c...gaugino masses        
      IF (IMODEL.EQ.1) THEN
        DO J=1,3
          IF (XNUSUG(J).LT.1.d19) THEN
C       Set possible non-universal boundary conditions
            G(J+30)=XNUSUG(J)
          ELSE
	    G(J+30)= DBLE(MHF)
	  END IF
        ENDDO
      ELSEIF (IMODEL.EQ.2) THEN
	DO J=1,3
          G(J+30)=XGMIN(11+J)*XGMIN(8)*THRG*(G(J)/4./PI)**2*XLAMGM
        ENDDO
      ELSEIF (IMODEL.EQ.7.OR.IMODEL.EQ.9) THEN
       G(31)= 33.d0*DBLE(MHF)*G(1)**2/5./16./PI**2
       G(32)= DBLE(MHF)*G(2)**2/16./PI**2
       G(33)=-3.d0*DBLE(MHF)*G(3)**2/16./PI**2
      ENDIF

C...Compute anomaly mediated functions
      IF (IMODEL.EQ.7.OR.IMODEL.EQ.9) THEN
        BLHAT=G(30)*(-9*G(1)**2/5.-3*G(2)**2+3*G(21)**2+4*G(30)**2)
        BBHAT=G(21)*(-7*G(1)**2/15.-3*G(2)**2-16*G(3)**2/3.+
     ,             G(12)**2+6*G(21)**2+G(30)**2)
        BTHAT=G(12)*(-13*G(1)**2/15.-3*G(2)**2-16*G(3)**2/3.+
     ,             6*G(12)**2+G(21)**2)
      ENDIF

c...soft trilinar couplings      
      IF (IMODEL.EQ.1) THEN
        DO J=1,3
          IF (XNUSUG(J+3).LT.1.d19) THEN
C       Set possible non-universal boundary conditions
            A(J)=XNUSUG(J+3)
          ELSE
	    A(J)= DBLE(A0)
          ENDIF
        ENDDO
      ELSE IF (IMODEL.EQ.2) THEN
        DO J=1,3
          A(J)=0.
	ENDDO
      ELSEIF (IMODEL.EQ.7.OR.IMODEL.EQ.9) THEN
        A(3)=-BLHAT*MHF/G(30)/16./PI**2
        A(2)=-BBHAT*MHF/G(21)/16./PI**2
        A(1)=-BTHAT*MHF/G(12)/16./PI**2
      ENDIF
      DO I=1,9
        G(33+i) = G(3+i)*A(1)
        G(42+i) = G(12+i)*A(2)
        G(51+i) = G(21+i)*A(3)
      ENDDO
c...higgs mixing parameters      
      G(61)= G(61)
      G(62)= G(62)
c...higgs  mass^2
      IF (IMODEL.EQ.1) THEN
        DO J=1,2
          IF (XNUSUG(J+6).LT.1.E19) THEN
C       Set possible non-universal boundary conditions
C-FP    ANSI/gfortran fix
            G(J+62)=SIGN(1.D0,XNUSUG(J+6))*(XNUSUG(J+6))**2
          ELSE
            G(J+62)=DBLE(M0**2)
          ENDIF
        ENDDO
        IF (INUHM.EQ.1) THEN
          G(63) = MHUSMG
          G(64) = MHDSMG
        ENDIF
      ELSEIF (IMODEL.EQ.2) THEN
        G(64)=XC*(.75*XGMIN(13)*(G(2)/4./PI)**4+.6*.25*
     ,             XGMIN(12)*(G(1)/4./PI)**4)+XGMIN(9)-DY
        G(63)=XC*(.75*XGMIN(13)*(G(2)/4./PI)**4+.6*.25*
     ,             XGMIN(12)*(G(1)/4./PI)**4)+XGMIN(10)+DY
      ELSEIF (IMODEL.EQ.7.OR.IMODEL.EQ.9) THEN
        G(64)=(-99*G(1)**4/50.-3*G(2)**4/2.+3*G(21)*BBHAT+G(30)*BLHAT)*
     ,        MHF**2/(16*PI**2)**2+XAMIN(6)*M0**2
        G(63)=(-99*G(1)**4/50.-3*G(2)**4/2.+3*G(12)*BTHAT)*
     ,        MHF**2/(16*PI**2)**2+XAMIN(7)*M0**2
      ENDIF
c
c   Set the rest of SSB mass^2      
c
c...Reset mass^2 elements
      DO I=65,147
  	IF(I.LE.109.OR.I.GE.121) G(I)=0.d0
      ENDDO

c...diagonal elements of m^2_Q
      IF (IMODEL.EQ.1) THEN
        IF (XNUSUG(13).LT.1.d19) THEN
C-FP      ANSI/gfortran fix
          G(65)=SIGN(1.D0,XNUSUG(13))*(XNUSUG(13))**2
        ELSE
          G(65)= DBLE(M0**2)
        ENDIF
        G(69)=G(65)   ! degernerate 1st and 2nd generations
        IF (XNUSUG(18).LT.1.d19) THEN
C-FP      ANSI/gfortran fix
          G(73)=SIGN(1.D0,XNUSUG(18))*(XNUSUG(18))**2
        ELSE
          G(73)= DBLE(M0**2)
        ENDIF
      ELSEIF (IMODEL.EQ.2) THEN
        G(65)=XC*(4*XGMIN(14)*(G(3)/4./PI)**4/3.+.75*XGMIN(13)*
     ,        (G(2)/4./PI)**4+.6*XGMIN(12)*(G(1)/4./PI)**4/36.)+DY/3.
        G(69)=G(65)
        G(73)=G(65)
      ELSEIF (IMODEL.EQ.7.OR.IMODEL.EQ.9) THEN
        G(65)=(-11*G(1)**4/50.-3*G(2)**4/2.+8*G(3)**4)*
     ,         MHF**2/(16*PI**2)**2+XAMIN(1)*M0**2
        G(69)=G(65)
        G(73)=(-11*G(1)**4/50.-3*G(2)**4/2.+8*G(3)**4+G(21)*BBHAT+
     ,         G(12)*BTHAT)*MHF**2/(16*PI**2)**2+XAMIN(1)*M0**2
      ENDIF  
c...diagonal elements of m^2_L
      IF (IMODEL.EQ.1) THEN
        IF (XNUSUG(10).LT.1.d19) THEN
C-FP      ANSI/gfortran fix
          G(74)=SIGN(1.D0,XNUSUG(10))*(XNUSUG(10))**2
        ELSE
          G(74)= DBLE(M0**2)
        END IF
        G(78)=G(74)   ! degernerate 1st and 2nd generations
        IF (XNUSUG(15).LT.1.d19) THEN
C-FP      ANSI/gfortran fix
          G(82)=SIGN(1.D0,XNUSUG(15))*(XNUSUG(15))**2
        ELSE
          G(82)= DBLE(M0**2)
        END IF
      ELSEIF (IMODEL.EQ.2) THEN
       G(74)=XC*(.75*XGMIN(13)*(G(2)/4./PI)**4+.6*.25*
     ,        XGMIN(12)*(G(1)/4./PI)**4)-DY
        G(78)=G(74) 
        G(82)=G(74) 
      ELSEIF (IMODEL.EQ.7.OR.IMODEL.EQ.9) THEN
        G(74)=(-99*G(1)**4/50.-3*G(2)**4/2.)*MHF**2/(16*PI**2)**2
     ,+XAMIN(4)*M0**2
        G(78)=G(74) 
        G(82)=(-99*G(1)**4/50.-3*G(2)**4/2.+G(30)*BLHAT)*
     ,        MHF**2/(16*PI**2)**2+XAMIN(4)*M0**2
      ENDIF  
c...diagonal elements of m^2_u
      IF (IMODEL.EQ.1) THEN
        IF (XNUSUG(12).LT.1.d19) THEN
C-FP      ANSI/gfortran fix
          G(83)=SIGN(1.D0,XNUSUG(12))*(XNUSUG(12))**2
        ELSE
          G(83)= DBLE(M0**2)
        END IF
        G(87)=G(83)   ! degernerate 1st and 2nd generations
        IF (XNUSUG(17).LT.1.E19) THEN
C-FP      ANSI/gfortran fix
          G(91)=SIGN(1.D0,XNUSUG(17))*(XNUSUG(17))**2
        ELSE
          G(91)= DBLE(M0**2)
        END IF
      ELSEIF (IMODEL.EQ.2) THEN
        G(83)=XC*(4*XGMIN(14)*(G(3)/4./PI)**4/3.+.6*4*XGMIN(12)*
     ,        (G(1)/4./PI)**4/9.)-4*DY/3.
        G(87)=G(83)
        G(91)=G(83)
      ELSEIF (IMODEL.EQ.7.OR.IMODEL.EQ.9) THEN
        G(83)=(-88*G(1)**4/25.+8*G(3)**4)*MHF**2/(16*PI**2)**2+
     ,XAMIN(3)*M0**2
        G(87)=G(83)
        G(91)=(-88*G(1)**4/25.+8*G(3)**4+2*G(12)*BTHAT)*
     , MHF**2/(16*PI**2)**2+XAMIN(3)*M0**2
      ENDIF  
c...diagonal elements of m^2_d
      IF (IMODEL.EQ.1) THEN
        IF (XNUSUG(11).LT.1.d19) THEN
C-FP      ANSI/gfortran fix
          G(92)=SIGN(1.D0,XNUSUG(11))*(XNUSUG(11))**2
        ELSE
          G(92)= DBLE(M0**2)
        END IF
        G(96)=G(92)   ! degernerate 1st and 2nd generations
        IF (XNUSUG(16).LT.1.d19) THEN
C-FP      ANSI/gfortran fix
          G(100)=SIGN(1.D0,XNUSUG(16))*(XNUSUG(16))**2
        ELSE
          G(100)= DBLE(M0**2)
        END IF
      ELSEIF (IMODEL.EQ.2) THEN
        G(92)=XC*(4*XGMIN(14)*(G(3)/4./PI)**4/3.+.6*XGMIN(12)*
     , (G(1)/4./PI)**4/9.)+2*DY/3.
        G(96)=G(92)  
        G(100)=G(92) 
      ELSEIF (IMODEL.EQ.7.OR.IMODEL.EQ.9) THEN
        G(92)=(-22*G(1)**4/25.+8*G(3)**4)*MHF**2/(16*PI**2)**2+
     ,XAMIN(2)*M0**2
        G(96)=G(92)  
        G(100)=(-22*G(1)**4/25.+8*G(3)**4+2*G(21)*BBHAT)*
     , MHF**2/(16*PI**2)**2+XAMIN(2)*M0**2
      ENDIF  
c...diagonal elements of m^2_e
      IF (IMODEL.EQ.1) THEN
        IF (XNUSUG(9).LT.1.d19) THEN
C-FP      ANSI/gfortran fix
          G(101)=SIGN(1.D0,XNUSUG(9))*(XNUSUG(9))**2
        ELSE
          G(101)= DBLE(M0**2)
        END IF
        G(105)=G(101)   ! degernerate 1st and 2nd generations
        IF (XNUSUG(14).LT.1.d19) THEN
C-FP      ANSI/gfortran fix
          G(109)=SIGN(1.D0,XNUSUG(16))*(XNUSUG(14))**2
        ELSE
          G(109)= DBLE(M0**2)
        END IF
      ELSEIF (IMODEL.EQ.2) THEN
        G(101)=XC*(.6*XGMIN(12)*(G(1)/4./PI)**4)+2*DY
        G(105)=G(101) 
        G(109)=G(101) 
      ELSEIF (IMODEL.EQ.7.OR.IMODEL.EQ.9) THEN
        G(101)=(-198*G(1)**4/25.)*MHF**2/(16*PI**2)**2+XAMIN(5)*M0**2
        G(105)=G(101) 
        G(109)=(-198*G(1)**4/25.+2*G(30)*BLHAT)*MHF**2/(16*PI**2)**2+
     ,XAMIN(5)*M0**2
      ENDIF  

c...additional contributions for mixed moduli-AMSB model
      IF (IMODEL.EQ.9) THEN
        DO I=1,3
          DG(I)=G(I)
        ENDDO
	DG(4)=G(30)
	DG(5)=G(21)
	DG(6)=G(12)
	DG(7)=G(31)
	DG(8)=G(32)
	DG(9)=G(33)
	DG(9)=G(33)
	DG(10)=G(60)/G(30)
	DG(11)=G(51)/G(21)
	DG(12)=G(42)/G(12)
	DG(13)=G(64)
	DG(14)=G(63)
	DG(16)=G(101)
	DG(17)=G(74)
	DG(18)=G(92)
	DG(19)=G(83)
	DG(20)=G(109)
	DG(21)=G(82)
	DG(22)=G(100)
	DG(23)=G(91)
	DG(24)=G(73)
	
        CALL MMAMSB(M0,MHF,DG)
        
	G(31)=DG(7)
	G(32)=DG(8)
	G(33)=DG(9)
	G(33)=DG(9)
	G(60)=DG(10)*G(30)
	G(51)=DG(11)*G(21)
	G(42)=DG(12)*G(12)
	G(64)=DG(13)
	G(63)=DG(14)
	G(101)=DG(16)
	G(74)=DG(17)
	G(92)=DG(18)
	G(83)=DG(19)
	G(109)=DG(20)
	G(82)=DG(21)
	G(100)=DG(22)
	G(91)=DG(23)
	G(73)=DG(24)
      ENDIF

c...higgs VEVs
      G(110) = G(110)
      G(111) = G(111)
C
C  neutrino sector
C
      IF (IRHN.GT.0) THEN
c...Yukawa matrix
        IF (IRHN.EQ.2) THEN     ! impose YN-YU unification with degree EPSNU
          DO I=1,3
	    DO J=1,3
	      YNGUT(I,J)=YUGUT(I,J)*EPSNU
	    ENDDO
	  ENDDO  
        ELSEIF (IRHN.EQ.3) THEN
          K=0
          DO I=1,3
            DO J=1,3
	      YNGUT(J,I)=XRHNIN(10+K)
	      K=K+1
            ENDDO
          ENDDO
        ENDIF
	CALL MAT2VEC(G,112,YNGUT,-1)
c...RHN majorana mass matrix
        DO I=1,9
          G(120+I) = XRHNIN(I) 
        ENDDO
c...soft trilinear coupling
        DO I=0,8
          G(130+I) = G(112+I)*XRHNIN(19)
        ENDDO
c...diagonal elements of SSB mass^2
C-FP    ANSI/gfortran fix
        G(139) = SIGN(1.D0,XRHNIN(11))*(XRHNIN(20))**2
        G(143) = G(139)             ! degernerate 1st and 2nd generations
        G(147) = SIGN(1.D0,XRHNIN(12))*(XRHNIN(21))**2
c...dim-5 related terms
        DO I=148,157
	  G(I)=0.d0
	ENDDO
      ENDIF
      
      RETURN
      END
