CDECK  ID>, ST1CNEU.
      SUBROUTINE ST1CNEU(G,Q)
!
!Purpose: To find the stop1 decay rate to charm+neutralino according
!         to the one-step approximation and according to the more
!         precise calculation from the full solution of the RGEs.
!
      IMPLICIT NONE
!
      COMMON/RGEIN/MHIGH,PHASEMU,ACC,COMP,SUG,UNI
      DOUBLE PRECISION MHIGH,PHASEMU
      INTEGER ACC,COMP,SUG,UNI
      SAVE/RGEIN/
!
      COMMON/MYDECAY/MQQMASS,MUQMASS,MDQMASS,MLQMASS,MEQMASS,
     $             OFFMAXQVAL,OFFMAXUVAL,OFFMAXDVAL,OFFMAXLVAL,
     $             OFFMAXEVAL,OFFMAXQ,OFFMAXU,OFFMAXD,OFFMAXL,OFFMAXE
      DOUBLE COMPLEX MQQMASS(3,3),MUQMASS(3,3),MDQMASS(3,3),
     $               MLQMASS(3,3),MEQMASS(3,3)
      DOUBLE COMPLEX OFFMAXQVAL,OFFMAXUVAL,OFFMAXDVAL,OFFMAXLVAL,
     $               OFFMAXEVAL
      INTEGER OFFMAXQ(2),OFFMAXU(2),OFFMAXD(2),OFFMAXL(2),OFFMAXE(2)
      SAVE/MYDECAY/
!
      COMMON/DECCALC/T1EVE,T1EVA,USQM,COSTHT,SINTHT,GHIK,MST1,MST2,GAMMA
      DOUBLE COMPLEX T1EVE(6,6),T1EVA(6),USQM(6,6),COSTHT,SINTHT
     $              ,GHIK(601)
      DOUBLE PRECISION MST1,MST2,GAMMA
      SAVE/DECCALC/
!
C          SUSY parameters
C          AMGLSS               = gluino mass
C          AMULSS               = up-left squark mass
C          AMELSS               = left-selectron mass
C          AMERSS               = right-slepton mass
C          AMNiSS               = sneutrino mass for generation i
C          TWOM1                = Higgsino mass = - mu
C          RV2V1                = ratio v2/v1 of vev's
C          AMTLSS,AMTRSS        = left,right stop masses
C          AMT1SS,AMT2SS        = light,heavy stop masses
C          AMBLSS,AMBRSS        = left,right sbottom masses
C          AMB1SS,AMB2SS        = light,heavy sbottom masses
C          AMLLSS,AMLRSS        = left,right stau masses
C          AML1SS,AML2SS        = light,heavy stau masses
C          AMZiSS               = signed mass of Zi
C          ZMIXSS               = Zi mixing matrix
C          AMWiSS               = signed Wi mass
C          GAMMAL,GAMMAR        = Wi left, right mixing angles
C          AMHL,AMHH,AMHA       = neutral Higgs h0, H0, A0 masses
C          AMHC                 = charged Higgs H+ mass
C          ALFAH                = Higgs mixing angle
C          AAT                  = stop trilinear term
C          THETAT               = stop mixing angle
C          AAB                  = sbottom trilinear term
C          THETAB               = sbottom mixing angle
C          AAL                  = stau trilinear term
C          THETAL               = stau mixing angle
C          AMGVSS               = gravitino mass
C          MTQ                  = top mass at MSUSY
C          MBQ                  = bottom mass at MSUSY
C          MLQ                  = tau mass at MSUSY
C          FBMA                 = b-Yukawa at mA scale
C          VUQ                  = Hu vev at MSUSY
C          VDQ                  = Hd vev at MSUSY
C          SGNM3                = sign of gaugino mass M3
      COMMON/SSPAR/GORGE,AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,AMSLSS
     $,AMSRSS,AMCLSS,AMCRSS,AMBLSS,AMBRSS,AMB1SS,AMB2SS
     $,AMTLSS,AMTRSS,AMT1SS,AMT2SS,AMELSS,AMERSS,AMMLSS,AMMRSS
     $,AMLLSS,AMLRSS,AML1SS,AML2SS,AMN1SS,AMN2SS,AMN3SS
     $,TWOM1,RV2V1,AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,ZMIXSS(4,4)
     $,AMW1SS,AMW2SS
     $,GAMMAL,GAMMAR,AMHL,AMHH,AMHA,AMHC,ALFAH,AAT,THETAT
     $,AAB,THETAB,AAL,THETAL,AMGVSS,MTQ,MBQ,MLQ,FBMA,
     $VUQ,VDQ,SGNM3
      REAL AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,AMSLSS
     $,AMSRSS,AMCLSS,AMCRSS,AMBLSS,AMBRSS,AMB1SS,AMB2SS
     $,AMTLSS,AMTRSS,AMT1SS,AMT2SS,AMELSS,AMERSS,AMMLSS,AMMRSS
     $,AMLLSS,AMLRSS,AML1SS,AML2SS,AMN1SS,AMN2SS,AMN3SS
     $,TWOM1,RV2V1,AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,ZMIXSS
     $,AMW1SS,AMW2SS
     $,GAMMAL,GAMMAR,AMHL,AMHH,AMHA,AMHC,ALFAH,AAT,THETAT
     $,AAB,THETAB,AAL,THETAL,AMGVSS,MTQ,MBQ,MLQ,FBMA,VUQ,VDQ,SGNM3
      LOGICAL GORGE
      REAL AMZISS(4)
      EQUIVALENCE (AMZISS(1),AMZ1SS)
      SAVE /SSPAR/
!
      COMMON/THRESH/QTHSORT,QTHQL,QTHUR,QTHDR,QTHLL,QTHER,QNSH,QNSG,
     $              QNH,QTHSB,QTHSW,EPS,LOCMH
      DOUBLE PRECISION QTHSORT(20),QTHQL(3),QTHUR(3),QTHDR(3),QTHLL(3),
     $                QTHER(3),QNSH,QNSG,QNH,QTHSB,QTHSW,EPS
                                     !EPSILON IS USED WHEN Q=THRESHOLD
      INTEGER LOCMH
      SAVE/THRESH/
!
      DOUBLE PRECISION VECTZ1(4),MZ,MC,MTMZ,MZ1,ALSQ,BESQ,LMASS
      DOUBLE COMPLEX KMTB,KMCB,DELTAL,DELTAR,EPSILON
      DOUBLE COMPLEX G(601),CROSS,IM,A(3),B(3),AL,BE
      DOUBLE COMPLEX FU(3,3),FD(3,3),AD(3,3),MQ(3,3),MD(3,3),MHD
      DOUBLE COMPLEX GTPQ(3,3),GTPU(3,3),GTQ(3,3),FTUQ(3,3),FTUU(3,3)
      DOUBLE COMPLEX FTUUT(3,3),FTUQS(3,3)
      DOUBLE COMPLEX GORIG(601)
      DOUBLE PRECISION PI,TANB,SB,Q,HKGAMMA,MW,QHIK
      INTEGER I,J,THETA1
!
!Rotate GHIK which for the rest of this must be in the quark mass
!basis. Must first rotate from the basis where the up Yukawas are
!diagonal to the original current basis, and then to the quark mass
!basis.
!
      CALL STROTATE(GHIK,GORIG,0)
      CALL STROTBACK(GORIG,GHIK,1)
      QHIK=Q
      IF(Q.LT.QNH)QHIK=QNH
!
      MHD=GHIK(62)-GHIK(108)**2
!
      IM=(0.D0,1.D0)
      PI=4.D0*DATAN(1.D0)
!
!Accuracy is not important for the two KM elements used in this
!subroutine, so they are entered here as constants rather than
!passed from RGEFLAV
!
      KMTB=(0.9991D0,0.D0)
      KMCB=(0.0413D0,0.D0)
!
      MZ1=ABS(DBLE(AMZ1SS))
      DO I=1,4
        VECTZ1(I)=DBLE(ZMIXSS(I,1))
            !The order of (row,column) has been confirmed with isajet
      END DO
!
      DO I=1,3
        DO J=1,3
!
!This first set is used for the HK calculation.
!We need the quark mass basis terms.
!
          FU(J,I)=GHIK(3+(I-1)*3+J)
          FD(I,J)=GHIK(12+(I-1)*3+J)
          AD(I,J)=GHIK(42+(I-1)*3+J)
          MQ(I,J)=MQQMASS(I,J)
          MD(I,J)=MDQMASS(I,J)
!
!This set is used for the full RGE calculation.
!
          GTPQ(I,J)=G(138+(I-1)*3+J)
          GTPU(I,J)=G(156+(I-1)*3+J)
          GTQ(I,J)=G(185+(I-1)*3+J)
          FTUQ(I,J)=G(232+(I-1)*3+J)
          FTUU(I,J)=G(259+(I-1)*3+J)
        END DO
      END DO
      DO I=1,3
        DO J=1,3
          FTUUT(I,J)=FTUU(J,I)
          FTUQS(I,J)=CONJG(FTUQ(I,J))
        END DO
      END DO
!
      TANB=DBLE(G(110)/G(111))
      SB=SQRT(TANB**2/(1.D0+TANB**2))
      MC=0.677D0 !This is a weak scale value not the value at MS
!
      IF(AMZ1SS.GT.0)THEN
        THETA1=0
      ELSE
        THETA1=1
      END IF
!
!HK (single-step) gamma first
!
      DELTAL=-DBLE(LOG(REAL(MHIGH)/REAL(QHIK)))/16.D0/PI**2*KMTB
     $       *CONJG(KMCB)*((MQ(2,2)+MQ(3,3)+2.D0*MHD+2.D0*MD(3,3))
     $       *CONJG(FD(3,3))*FD(3,3)+2.D0*AD(3,3)*CONJG(AD(3,3)))
      DELTAR=2.D0*DBLE(LOG(REAL(MHIGH)/REAL(QHIK)))/16.D0/PI**2*KMTB
     $       *CONJG(KMCB)*CONJG(FD(3,3))*GHIK(110)*FU(3,3)*AD(3,3)
      EPSILON=(CONJG(DELTAL)*COSTHT-CONJG(DELTAR)*SINTHT
     $          )/(MST1**2-USQM(2,2))
      HKGAMMA=1.D0/16.D0/PI*EPSILON*CONJG(EPSILON)*MST1
     $        *(1.D0-(MZ1**2/MST1**2))**2*1.D0/2.D0
     $        *(DBLE(G(2))*VECTZ1(3)+1.D0/3.D0*DBLE(G(1))
     $        *DSQRT(3.D0/5.D0)*VECTZ1(4))**2
!
!Now my gamma
!First check the kinematics
!
      IF(MST1.LT.(MC+MZ1))THEN
        GAMMA=0.D0
        WRITE(*,*)'DECAY IS KINEMATICALLY FORBIDDEN'
        RETURN
      END IF
!
      DO I=1,3
        A(I)=(-IM)**(THETA1-1)/DSQRT(2.D0)*(GTQ(I,2)*VECTZ1(3)
     $       +GTPQ(I,2)*VECTZ1(4)/3.D0)
        B(I)=4.D0/(3.D0*DSQRT(2.D0))*(CONJG(GTPU(2,I))*(IM)**(THETA1-1)
     $       *VECTZ1(4))
      END DO
!
      LMASS=MST1**4+MZ1**4+MC**4-2.D0*MST1**2*MZ1**2-2.D0*MST1**2
     $        *MC**2-2.D0*MZ1**2*MC**2
!
      AL=(0.D0,0.D0)
      BE=(0.D0,0.D0)
      DO I=1,3
        AL=AL+IM*A(I)*CONJG(T1EVE(I,1))
     $        -FTUUT(I,2)*VECTZ1(1)*(-IM)**THETA1*CONJG(T1EVE(I+3,1))
        BE=BE+IM*B(I)*CONJG(T1EVE(I+3,1))
     $        -FTUQS(I,2)*VECTZ1(1)*(IM)**THETA1*CONJG(T1EVE(I,1))
      END DO
!
      ALSQ=DBLE(AL*CONJG(AL))
      BESQ=DBLE(BE*CONJG(BE))
      CROSS=AL*CONJG(BE)+BE*CONJG(AL)
!
      GAMMA=1.D0/16.D0/PI/MST1**3*((ALSQ+BESQ)*(MST1**2-MC**2-MZ1**2)
     $      -2.D0*MC*MZ1*DBLE(CROSS))*SQRT(LMASS)
!
      WRITE(*,*)'ONE-STEP ESTIMATE OF GAMMA(TP1->Z1SS CH) IS: ',HKGAMMA
      WRITE(*,*)'FULL RGE CALCULATION OF GAMMA(TP1->Z1SS CH) IS: ',GAMMA
!
      RETURN
      END
