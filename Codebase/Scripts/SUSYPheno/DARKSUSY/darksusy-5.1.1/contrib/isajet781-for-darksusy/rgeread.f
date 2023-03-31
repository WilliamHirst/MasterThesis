!
      SUBROUTINE RGEREAD
!
!Purpose: To read in the inputs needed by RGEFLAV and populate the
!         common blocks.
!
!If the user decides to enter complex variables, they must be entered as
!"(x,y)" where z=x+iy.
!
      IMPLICIT NONE
!
      COMMON/SMSAVED/KM,MWEAK,MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      DOUBLE COMPLEX KM(3,3)
      DOUBLE PRECISION MWEAK(6),MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      SAVE/SMSAVED/
!
      COMMON/MYSUGRA/M0,M12,A0,TANB,SIGNMU,MT
      DOUBLE PRECISION M0,M12,A0,TANB,SIGNMU,MT
      SAVE/MYSUGRA/
!
!Contains the unitary matrices that are used in the rotation
!
      COMMON/UNITARY/VLU,VRU,VLD,VRD,SVLQ
      DOUBLE COMPLEX VLU(3,3),VRU(3,3),VLD(3,3),VRD(3,3)
      INTEGER SVLQ
      SAVE/UNITARY/
!
      COMMON/RGEIN/MHIGH,PHASEMU,ACC,COMP,SUG,UNI
      DOUBLE PRECISION MHIGH,PHASEMU
      INTEGER ACC,COMP,SUG,UNI
      SAVE/RGEIN/
      COMMON /RGEFNM/ FNRGE
      CHARACTER*128 FNRGE,FNINPUT,STRADD
!
!High scale SUSY inputs
!
      COMMON/HISUSY/CM1,CM2,CM3,CM1P,CM2P,CM3P
     $             ,CMQ0,CMU0,CMD0,CML0,CME0,CRU,CRD,CRE,CSU,CSD,CSE
     $             ,CTQ,CTU,CTD,CTL,CTE
     $             ,CAU0,CAD0,CAE0,CWU,CWD,CWE,CXU,CXD,CXE
     $             ,CZU,CZD,CZE,CMHU,CMHD,CU,CD,CE
      DOUBLE COMPLEX CM1,CM2,CM3,CM1P,CM2P,CM3P
     $              ,CMQ0,CMU0,CMD0,CML0,CME0,CRU,CRD,CRE,CSU,CSD,CSE
     $              ,CTQ(3,3),CTU(3,3),CTD(3,3),CTL(3,3),CTE(3,3)
     $              ,CAU0,CAD0,CAE0,CWU,CWD,CWE,CXU,CXD,CXE
     $              ,CZU(3,3),CZD(3,3),CZE(3,3),CMHU,CMHD
      INTEGER CU,CD,CE
      SAVE/HISUSY/
!
      COMMON/DEC/NEWTH,ISADEC,BELOW,NSTEPTHRESH,NLTMT,
     $           THSQ,THSU,THSD,THSL,THSE
      DOUBLE PRECISION NEWTH(20)
      INTEGER ISADEC,BELOW(20),NSTEPTHRESH(19),NLTMT
      INTEGER THSQ(3),THSU(3),THSD(3),THSL(3),THSE(3)
      SAVE/DEC/
!
      DOUBLE COMPLEX COMBM1,COMBM2,COMBM3,DAGKM(3,3),CMATMUL
      DOUBLE PRECISION ALP,BET,GAM,DEL,ID(3,3)
      DOUBLE PRECISION DM1,DM2,DM3
     $                ,DMQ0,DMU0,DMD0,DML0,DME0,DRU,DRD,DRE,DSU,DSD,DSE
     $                ,DTQ(3,3),DTU(3,3),DTD(3,3),DTL(3,3),DTE(3,3)
     $                ,DAU0,DAD0,DAE0,DWU,DWD,DWE,DXU,DXD,DXE
     $                ,DZU(3,3),DZD(3,3),DZE(3,3),DMHU,DMHD
      INTEGER I,J,VRID,VLRIN,VLKM,TMP,ERRNO,ERRNOMR,HERMBAD
      CHARACTER COMPTEST
!
      DATA ID(1,1)/1.D0/,ID(1,2)/0.D0/,ID(1,3)/0.D0/
      DATA ID(2,1)/0.D0/,ID(2,2)/1.D0/,ID(2,3)/0.D0/
      DATA ID(3,1)/0.D0/,ID(3,2)/0.D0/,ID(3,3)/1.D0/
!
      ERRNO=0
      ERRNOMR=0
!
      DO I=1,3
        DO J=1,3
          DTQ(I,J)=0.D0
          DTU(I,J)=0.D0
          DTD(I,J)=0.D0
          DTL(I,J)=0.D0
          DTE(I,J)=0.D0
          DZU(I,J)=0.D0
          DZD(I,J)=0.D0
          DZE(I,J)=0.D0
        END DO
      END DO
!
!M_Z inputs. Start with ALPHAEM and SIN2THW from the PDG
!Pole masses from PDG rather than running masses.
!Updated from PDG 2008. For ALPHAEM, see Sec.10: EW Model and
!Constraints on New Physics
!
      ALPHAEM=1.D0/127.925D0
      ALPHASMSB=.1176D0
      XWMSB=.23119D0
      MZ=91.1876D0!Used in later subroutines, could be passed
                  !from an earlier programme
      MW=80.403D0
!
!Now read from the input file
!
      FNINPUT=STRADD(FNRGE,'.rgein')
      OPEN(12,FILE=FNINPUT,STATUS='UNKNOWN')
      READ(12,*)
      READ(12,*)ACC
      IF(ACC.EQ.0)WRITE(*,*)'USING ONE LOOP RGES'
      DO I=1,2
        READ(12,*)
      END DO
      READ(12,*)COMP
      IF(COMP.EQ.0)WRITE(*,*)'USING REAL KM AND RUNNING'
      DO I=1,2
        READ(12,*)
      END DO
      READ(12,*)TMP
      IF(TMP.EQ.0.OR.COMP.EQ.0)THEN
        DO I=1,3
          READ(12,*)
        END DO
      ELSE
        DO I=1,2
          READ(12,*)
        END DO
        READ(12,*)PHASEMU
      END IF
      DO I=1,2
        READ(12,*)
      END DO
      READ(12,*)ISADEC
!
!Now we know COMP, we can call the KM routine
!
      CALL KMIN
      CALL CDAGGER(KM,DAGKM)
!
      DO I=1,2
        READ(12,*)
      END DO
      READ(12,*)SUG
!
      IF(SUG.EQ.1)THEN
        UNI=1
        DO I=1,169
          READ(12,*)
        END DO
      ELSE
!
!Set all the inputs to follow to their mSUGRA defaults.
!It is only if the user chooses non-mSURA conditions that
!these defaults won't be used.
!
        UNI=1
        CM1=DCMPLX(M12,0.D0)
        CM2=DCMPLX(M12,0.D0)
        CM3=DCMPLX(M12,0.D0)
        CM1P=(0.D0,0.D0)
        CM2P=(0.D0,0.D0)
        CM3P=(0.D0,0.D0)
        CMHU=DCMPLX(M0)
        CMHD=DCMPLX(M0)
        CMQ0=DCMPLX(M0)
        CMU0=DCMPLX(M0)
        CMD0=DCMPLX(M0)
        CML0=DCMPLX(M0)
        CME0=DCMPLX(M0)
        CU=1
        CD=1
        CE=1
        CRU=(0.D0,0.D0)
        CRD=(0.D0,0.D0)
        CRE=(0.D0,0.D0)
        CSU=(0.D0,0.D0)
        CSD=(0.D0,0.D0)
        CSE=(0.D0,0.D0)
        CAU0=DCMPLX(A0)
        CAD0=DCMPLX(A0)
        CAE0=DCMPLX(A0)
        CWU=(0.D0,0.D0)
        CWD=(0.D0,0.D0)
        CWE=(0.D0,0.D0)
        CXU=(0.D0,0.D0)
        CXD=(0.D0,0.D0)
        CXE=(0.D0,0.D0)
        DO I=1,3
          DO J=1,3
            CTQ(I,J)=(0.D0,0.D0)
            CTU(I,J)=(0.D0,0.D0)
            CTD(I,J)=(0.D0,0.D0)
            CTL(I,J)=(0.D0,0.D0)
            CTE(I,J)=(0.D0,0.D0)
            CZU(I,J)=(0.D0,0.D0)
            CZD(I,J)=(0.D0,0.D0)
            CZE(I,J)=(0.D0,0.D0)
          END DO
        END DO
!
        DO I=1,4
          READ(12,*)
        END DO
        READ(12,*)TMP
        IF(TMP.EQ.0)THEN
          UNI=0
          DO I=1,2
            READ(12,*)
          END DO
          READ(12,*)MHIGH
        ELSE
          DO I=1,3
            READ(12,*)
          END DO
        END IF
!
!First the complex read statements. Each variable can either be the
!same as the mSUGRA values or a user defined value
!
        IF(COMP.EQ.1)THEN
          DO I=1,2
            READ(12,*)
          END DO
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,2
             READ(12,*)
            END DO
!
!The following three lines check if the user has input the numbers
!in the correct form. If not, there is an error message and the
!programme stops
!
            READ(12,*)COMPTEST
            IF(COMPTEST.NE.'(')THEN
              ERRNO=1
              GOTO 50
            END IF
            BACKSPACE 12 !NOW GO BACK TO READ THE LINE AGAIN
            READ(12,*)COMBM1,COMBM2,COMBM3
            CM1=DCMPLX(DBLE(COMBM1),0.D0)
            CM2=DCMPLX(DBLE(COMBM2),0.D0)
            CM3=DCMPLX(DBLE(COMBM3),0.D0)
            CM1P=DCMPLX(DIMAG(COMBM1),0.D0)
            CM2P=DCMPLX(DIMAG(COMBM2),0.D0)
            CM3P=DCMPLX(DIMAG(COMBM3),0.D0)
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,5
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.NE.'(')THEN
              ERRNO=2
              GOTO 50
            END IF
            BACKSPACE 12
            READ(12,*)CMHU,CMHD
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,5
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,3
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=3
              ERRNOMR=1
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DMQ0
            CMQ0=DCMPLX(DMQ0)
            DO I=1,2
              READ(12,*)
            END DO
            DO I=1,3
              READ(12,*)COMPTEST
              IF(COMPTEST.NE.'(')THEN
                ERRNO=3
                ERRNOMR=4
                GOTO 50
              END IF
              BACKSPACE 12
              READ(12,*)CTQ(I,1),CTQ(I,2),CTQ(I,3)
            END DO
            CALL HERMTEST(CTQ,HERMBAD)
            IF(HERMBAD.EQ.1)WRITE(*,*)'WARNING: T_Q IS NOT HERMITIAN'
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,11
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,4
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=4
              ERRNOMR=1
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DMU0
            CMU0=DCMPLX(DMU0)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)CU
            IF(CU.NE.1.AND.CU.NE.0)THEN
              WRITE(*,*)'ERROR READING CU'
              GOTO 52
            END IF
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=4
              ERRNOMR=2
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DRU
            CRU=DCMPLX(DRU)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=4
              ERRNOMR=3
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DSU
            CSU=DCMPLX(DSU)
            DO I=1,2
              READ(12,*)
            END DO
            DO I=1,3
              READ(12,*)COMPTEST
              IF(COMPTEST.NE.'(')THEN
                ERRNO=4
                ERRNOMR=4
                GOTO 50
              END IF
              BACKSPACE 12
              READ(12,*)CTU(I,1),CTU(I,2),CTU(I,3)
            END DO
            CALL HERMTEST(CTU,HERMBAD)
            IF(HERMBAD.EQ.1)WRITE(*,*)'WARNING: T_U IS NOT HERMITIAN'
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,21
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,4
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=5
              ERRNOMR=1
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DMD0
            CMD0=DCMPLX(DMD0)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)CD
            IF(CD.NE.1.AND.CD.NE.0)THEN
              WRITE(*,*)'ERROR READING CD'
              GOTO 52
            END IF
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=5
              ERRNOMR=2
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DRD
            CRD=DCMPLX(DRD)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=5
              ERRNOMR=3
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DSD
            CSD=DCMPLX(DSD)
            DO I=1,2
              READ(12,*)
            END DO
            DO I=1,3
              READ(12,*)COMPTEST
              IF(COMPTEST.NE.'(')THEN
                ERRNO=5
                ERRNOMR=4
                GOTO 50
              END IF
              BACKSPACE 12
              READ(12,*)CTD(I,1),CTD(I,2),CTD(I,3)
            END DO
            CALL HERMTEST(CTD,HERMBAD)
            IF(HERMBAD.EQ.1)WRITE(*,*)'WARNING: T_D IS NOT HERMITIAN'
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,21
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,3
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=6
              ERRNOMR=1
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DML0
            CML0=DCMPLX(DML0)
            DO I=1,2
              READ(12,*)
            END DO
            DO I=1,3
              READ(12,*)COMPTEST
              IF(COMPTEST.NE.'(')THEN
                ERRNO=6
                ERRNOMR=4
                GOTO 50
              END IF
              BACKSPACE 12
              READ(12,*)CTL(I,1),CTL(I,2),CTL(I,3)
            END DO
            CALL HERMTEST(CTL,HERMBAD)
            IF(HERMBAD.EQ.1)WRITE(*,*)'WARNING: T_L IS NOT HERMITIAN'
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,11
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,4
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=7
              ERRNOMR=1
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DME0
            CME0=DCMPLX(DME0)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)CE
            IF(CE.NE.1.AND.CE.NE.0)THEN
              WRITE(*,*)'ERROR READING CE'
              GOTO 52
            END IF
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=7
              ERRNOMR=2
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DRE
            CRE=DCMPLX(DRE)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=7
              ERRNOMR=3
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DSE
            CSE=DCMPLX(DSE)
            DO I=1,2
              READ(12,*)
            END DO
            DO I=1,3
              READ(12,*)COMPTEST
              IF(COMPTEST.NE.'(')THEN
                ERRNO=7
                ERRNOMR=4
                GOTO 50
              END IF
              BACKSPACE 12
              READ(12,*)CTE(I,1),CTE(I,2),CTE(I,3)
            END DO
            CALL HERMTEST(CTE,HERMBAD)
            IF(HERMBAD.EQ.1)WRITE(*,*)'WARNING: T_E IS NOT HERMITIAN'
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,21
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,4
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.NE.'(')THEN
              ERRNO=8
              GOTO 50
            END IF
            BACKSPACE 12
            READ(12,*)CAU0
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.NE.'(')THEN
              ERRNO=8
              GOTO 50
            END IF
            BACKSPACE 12
            READ(12,*)CWU
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.NE.'(')THEN
              ERRNO=8
              GOTO 50
            END IF
            BACKSPACE 12
            READ(12,*)CXU
            DO I=1,2
              READ(12,*)
            END DO
            DO I=1,3
              READ(12,*)COMPTEST
              IF(COMPTEST.NE.'(')THEN
                ERRNO=8
                GOTO 50
              END IF
              BACKSPACE 12
              READ(12,*)CZU(I,1),CZU(I,2),CZU(I,3)
            END DO
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,18
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,4
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.NE.'(')THEN
              ERRNO=9
              GOTO 50
            END IF
            BACKSPACE 12
            READ(12,*)CAD0
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.NE.'(')THEN
              ERRNO=9
              GOTO 50
            END IF
            BACKSPACE 12
            READ(12,*)CWD
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.NE.'(')THEN
              ERRNO=9
              GOTO 50
            END IF
            BACKSPACE 12
            READ(12,*)CXD
            DO I=1,2
              READ(12,*)
            END DO
            DO I=1,3
              READ(12,*)COMPTEST
              IF(COMPTEST.NE.'(')THEN
                ERRNO=9
                GOTO 50
              END IF
              BACKSPACE 12
              READ(12,*)CZD(I,1),CZD(I,2),CZD(I,3)
            END DO
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,18
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,4
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.NE.'(')THEN
              ERRNO=10
              GOTO 50
            END IF
            BACKSPACE 12
            READ(12,*)CAE0
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.NE.'(')THEN
              ERRNO=10
              GOTO 50
            END IF
            BACKSPACE 12
            READ(12,*)CWE
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.NE.'(')THEN
              ERRNO=10
              GOTO 50
            END IF
            BACKSPACE 12
            READ(12,*)CXE
            DO I=1,2
              READ(12,*)
            END DO
            DO I=1,3
              READ(12,*)COMPTEST
              IF(COMPTEST.NE.'(')THEN
                ERRNO=10
                GOTO 50
              END IF
              BACKSPACE 12
              READ(12,*)CZE(I,1),CZE(I,2),CZE(I,3)
            END DO
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,18
              READ(12,*)
            END DO
          END IF
!
!Now the real read statements. All variables are converted to complex
!at the end
!
        ELSE
          DO I=1,2
            READ(12,*)
          END DO
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,2
             READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=1
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DM1,DM2,DM3
            CM1=DCMPLX(DM1,0.D0)
            CM2=DCMPLX(DM2,0.D0)
            CM3=DCMPLX(DM3,0.D0)
            CM1P=(0.D0,0.D0)
            CM2P=(0.D0,0.D0)
            CM3P=(0.D0,0.D0)
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,5
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=2
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DMHU,DMHD
            CMHU=DCMPLX(DMHU,0.D0)
            CMHD=DCMPLX(DMHD,0.D0)
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,5
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,3
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=3
              ERRNOMR=1
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DMQ0
            CMQ0=DCMPLX(DMQ0)
            DO I=1,2
              READ(12,*)
            END DO
            DO I=1,3
              READ(12,*)COMPTEST
              IF(COMPTEST.EQ.'(')THEN
                ERRNO=3
                ERRNOMR=4
                GOTO 51
              END IF
              BACKSPACE 12
              READ(12,*)DTQ(I,1),DTQ(I,2),DTQ(I,3)
            END DO
            CALL SYMMTEST(DTQ,HERMBAD)
            IF(HERMBAD.EQ.1)WRITE(*,*)'WARNING: T_Q IS NOT SYMMETRIC'
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,11
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,4
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=4
              ERRNOMR=1
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DMU0
            CMU0=DCMPLX(DMU0)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)CU
            IF(CU.NE.1.AND.CU.NE.0)THEN
              WRITE(*,*)'ERROR READING CU'
              GOTO 52
            END IF
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=4
              ERRNOMR=2
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DRU
            CRU=DCMPLX(DRU)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=4
              ERRNOMR=3
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DSU
            CSU=DCMPLX(DSU)
            DO I=1,2
              READ(12,*)
            END DO
            DO I=1,3
              READ(12,*)COMPTEST
              IF(COMPTEST.EQ.'(')THEN
                ERRNO=4
                ERRNOMR=4
                GOTO 51
              END IF
              BACKSPACE 12
              READ(12,*)DTU(I,1),DTU(I,2),DTU(I,3)
            END DO
            CALL SYMMTEST(DTU,HERMBAD)
            IF(HERMBAD.EQ.1)WRITE(*,*)'WARNING: T_U IS NOT SYMMETRIC'
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,21
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,4
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=5
              ERRNOMR=1
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DMD0
            CMD0=DCMPLX(DMD0)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)CD
            IF(CD.NE.1.AND.CD.NE.0)THEN
              WRITE(*,*)'ERROR READING CD'
              GOTO 52
            END IF
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=5
              ERRNOMR=2
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DRD
            CRD=DCMPLX(DRD)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=5
              ERRNOMR=3
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DSD
            CSD=DCMPLX(DSD)
            DO I=1,2
              READ(12,*)
            END DO
            DO I=1,3
              READ(12,*)COMPTEST
              IF(COMPTEST.EQ.'(')THEN
                ERRNO=5
                ERRNOMR=4
                GOTO 51
              END IF
              BACKSPACE 12
              READ(12,*)DTD(I,1),DTD(I,2),DTD(I,3)
            END DO
            CALL SYMMTEST(DTD,HERMBAD)
            IF(HERMBAD.EQ.1)WRITE(*,*)'WARNING: T_D IS NOT SYMMETRIC'
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,21
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,3
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=6
              ERRNOMR=1
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DML0
            CML0=DCMPLX(DML0)
            DO I=1,2
              READ(12,*)
            END DO
            DO I=1,3
              READ(12,*)COMPTEST
              IF(COMPTEST.EQ.'(')THEN
                ERRNO=6
                ERRNOMR=4
                GOTO 51
              END IF
              BACKSPACE 12
              READ(12,*)DTL(I,1),DTL(I,2),DTL(I,3)
            END DO
            CALL SYMMTEST(DTL,HERMBAD)
            IF(HERMBAD.EQ.1)WRITE(*,*)'WARNING: T_L IS NOT SYMMETRIC'
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,11
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,4
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=7
              ERRNOMR=1
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DME0
            CME0=DCMPLX(DME0)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)CE
            IF(CE.NE.1.AND.CE.NE.0)THEN
              WRITE(*,*)'ERROR READING CE'
              GOTO 52
            END IF
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=7
              ERRNOMR=2
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DRE
            CRE=DCMPLX(DRE)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=7
              ERRNOMR=3
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DSE
            CSE=DCMPLX(DSE)
            DO I=1,2
              READ(12,*)
            END DO
            DO I=1,3
              READ(12,*)COMPTEST
              IF(COMPTEST.EQ.'(')THEN
                ERRNO=7
                ERRNOMR=4
                GOTO 51
              END IF
              BACKSPACE 12
              READ(12,*)DTE(I,1),DTE(I,2),DTE(I,3)
            END DO
            CALL SYMMTEST(DTE,HERMBAD)
            IF(HERMBAD.EQ.1)WRITE(*,*)'WARNING: T_E IS NOT SYMMETRIC'
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,21
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,4
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=8
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DAU0
            CAU0=DCMPLX(DAU0)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=8
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DWU
            CWU=DCMPLX(DWU)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=8
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DXU
            CXU=DCMPLX(DXU)
            DO I=1,2
              READ(12,*)
            END DO
            DO I=1,3
              READ(12,*)COMPTEST
              IF(COMPTEST.EQ.'(')THEN
                ERRNO=8
                GOTO 51
              END IF
              BACKSPACE 12
              READ(12,*)DZU(I,1),DZU(I,2),DZU(I,3)
            END DO
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,18
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,4
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=9
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DAD0
            CAD0=DCMPLX(DAD0)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=9
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DWD
            CWD=DCMPLX(DWD)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=9
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DXD
            CXD=DCMPLX(DXD)
            DO I=1,2
              READ(12,*)
            END DO
            DO I=1,3
              READ(12,*)COMPTEST
              IF(COMPTEST.EQ.'(')THEN
                ERRNO=9
                GOTO 51
              END IF
              BACKSPACE 12
              READ(12,*)DZD(I,1),DZD(I,2),DZD(I,3)
            END DO
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,18
              READ(12,*)
            END DO
          END IF
!
          READ(12,*)TMP
          IF(TMP.EQ.0)THEN
            DO I=1,4
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=10
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DAE0
            CAE0=DCMPLX(DAE0)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=10
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DWE
            CWE=DCMPLX(DWE)
            DO I=1,2
              READ(12,*)
            END DO
            READ(12,*)COMPTEST
            IF(COMPTEST.EQ.'(')THEN
              ERRNO=10
              GOTO 51
            END IF
            BACKSPACE 12
            READ(12,*)DXE
            CXE=DCMPLX(DXE)
            DO I=1,2
              READ(12,*)
            END DO
            DO I=1,3
              READ(12,*)COMPTEST
              IF(COMPTEST.EQ.'(')THEN
                ERRNO=10
                GOTO 51
              END IF
              BACKSPACE 12
              READ(12,*)DZE(I,1),DZE(I,2),DZE(I,3)
            END DO
            DO I=1,2
              READ(12,*)
            END DO
          ELSE
            DO I=1,18
              READ(12,*)
            END DO
          END IF
          DO I=1,3
            DO J=1,3
              CTQ(I,J)=DCMPLX(DTQ(I,J))
              CTU(I,J)=DCMPLX(DTU(I,J))
              CTD(I,J)=DCMPLX(DTD(I,J))
              CTL(I,J)=DCMPLX(DTL(I,J))
              CTE(I,J)=DCMPLX(DTE(I,J))
              CZU(I,J)=DCMPLX(DZU(I,J))
              CZD(I,J)=DCMPLX(DZD(I,J))
              CZE(I,J)=DCMPLX(DZE(I,J))
            END DO
          END DO
        END IF
      END IF
!
!Now deal with the rotation matrices
!
      READ(12,*)SVLQ !SVLQ=1 for up quark diagonal basis
      IF(SVLQ.NE.0.AND.SVLQ.NE.1)THEN
        WRITE(*,*)'ERROR: UNKNOWN CHOICE FOR QUARK DOUBLET ROTATION'
        WRITE(*,*)'OUTPUT WILL BE IN UP QUARK DIAGONAL BASIS'
        SVLQ=1
      END IF
      DO I=1,2
        READ(12,*)
      END DO
      READ(12,*)VLRIN !Does user input Vs?
      DO I=1,2
        READ(12,*)
      END DO
      READ(12,*)VLKM !Is V^U_L=KM or KM^dagger?
      DO I=1,2
        READ(12,*)
      END DO
      READ(12,*)VRID !Are V_Rs identity?
!
!Deal with user input Vs
!
      IF(VLRIN.EQ.1)THEN
        DO I=1,5
          READ(12,*)
        END DO
        READ(12,*)ALP,BET,GAM,DEL
        IF(COMP.EQ.0)DEL=0.D0 !This line makes sure that the V''s have
                              !no phase if the user chose real KM
        CALL VGEN(ALP,BET,GAM,DEL,VLU)
        DO I=1,2
          READ(12,*)
        END DO
        READ(12,*)ALP,BET,GAM,DEL
        IF(COMP.EQ.0)DEL=0.D0
        CALL VGEN(ALP,BET,GAM,DEL,VRU)
        DO I=1,2
          READ(12,*)
        END DO
        READ(12,*)ALP,BET,GAM,DEL
        IF(COMP.EQ.0)DEL=0.D0
        CALL VGEN(ALP,BET,GAM,DEL,VRD)
        DO I=1,3
          DO J=1,3
            VLD(I,J)=CMATMUL(0,VLU,KM,I,J)
          END DO
        END DO
!
!Set the Vs if they are not input by the user
!
      ELSE
!
!For the Left-handed matrices, just V^U_L is an input.
!
        IF(VLKM.EQ.0)THEN
          DO I=1,3
            DO J=1,3
              VLU(I,J)=DCMPLX(ID(I,J))
            END DO
          END DO
        ELSE
          DO I=1,3
            DO J=1,3
              VLU(I,J)=KM(I,J)
            END DO
          END DO
        END IF
!
!V^D_L is set so that (V^U_L)^dagger * V^D_L = KM
!
        DO I=1,3
          DO J=1,3
            VLD(I,J)=CMATMUL(0,VLU,KM,I,J)
          END DO
        END DO
!
!The V_Rs are arbitrary. If VRID=0 then the KM is used to set them
!
        IF(VRID.EQ.1)THEN
          DO I=1,3
             DO J=1,3
              VRU(I,J)=DCMPLX(ID(I,J))
              VRD(I,J)=DCMPLX(ID(I,J))
            END DO
          END DO
        ELSE
          DO I=1,3
            DO J=1,3
              VRU(I,J)=DAGKM(I,J)
              VRD(I,J)=KM(I,J)
            END DO
          END DO
        END IF
      END IF
!
      RETURN
!
  50  CONTINUE
      WRITE(*,*)'REAL NUMBERS READ WHERE COMPLEX NUMBERS EXPECTED'
      GOTO 52
  51  CONTINUE
      WRITE(*,*)'COMPLEX NUMBERS READ WHERE REAL NUMBERS EXPECTED'
  52  CONTINUE
      IF(ERRNO.EQ.1)THEN
        WRITE(*,*)'ERROR WAS IN M1,M2,M3'
      ELSE IF(ERRNO.EQ.2)THEN
        WRITE(*,*)'ERROR WAS IN M_{H_U}, M_{H_D}'
      ELSE IF(ERRNO.EQ.3)THEN
        WRITE(*,*)'ERROR WAS IN M_Q'
      ELSE IF(ERRNO.EQ.4)THEN
        WRITE(*,*)'ERROR WAS IN M_U'
      ELSE IF(ERRNO.EQ.5)THEN
        WRITE(*,*)'ERROR WAS IN M_D'
      ELSE IF(ERRNO.EQ.6)THEN
        WRITE(*,*)'ERROR WAS IN M_L'
      ELSE IF(ERRNO.EQ.7)THEN
        WRITE(*,*)'ERROR WAS IN M_E'
      ELSE IF(ERRNO.EQ.8)THEN
        WRITE(*,*)'ERROR WAS IN A_U'
      ELSE IF(ERRNO.EQ.9)THEN
        WRITE(*,*)'ERROR WAS IN A_D'
      ELSE IF(ERRNO.EQ.10)THEN
        WRITE(*,*)'ERROR WAS IN A_E'
      END IF
      IF(ERRNOMR.EQ.1)THEN
        WRITE(*,*)'M_0 WAS COMPLEX'
      ELSE IF(ERRNOMR.EQ.2)THEN
        WRITE(*,*)'R WAS COMPLEX'
      ELSE IF(ERRNOMR.EQ.3)THEN
        WRITE(*,*)'S WAS COMPLEX'
      ELSE IF(ERRNOMR.EQ.4)THEN
        WRITE(*,*)'T WAS OF THE WRONG KIND'
      END IF
      STOP99
!
      END
