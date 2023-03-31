!
      SUBROUTINE DOWNMSCOND
!
!Purpose: Apply various conditions at M_SUSY
!         The EWSB relations are imposed to derive various values.
!         Two versions depending on if m_H is less than m_SUSY
!
      IMPLICIT NONE
!
      COMMON/COUPLINGS/G,DG
      DOUBLE COMPLEX G(601)
      DOUBLE PRECISION DG(601)
      SAVE/COUPLINGS/
!
      COMMON/RGEMS/VEVMH,RGEMS,RGEMU
      DOUBLE COMPLEX VEVMH
      DOUBLE PRECISION RGEMS,RGEMU
      SAVE/RGEMS/
!
      COMMON/SMSAVED/KM,MWEAK,MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      DOUBLE COMPLEX KM(3,3)
      DOUBLE PRECISION MWEAK(6),MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      SAVE/SMSAVED/
!
      COMMON/THRESH/QTHSORT,QTHQL,QTHUR,QTHDR,QTHLL,QTHER,QNSH,QNSG,
     $              QNH,QTHSB,QTHSW,EPS,LOCMH
      DOUBLE PRECISION QTHSORT(20),QTHQL(3),QTHUR(3),QTHDR(3),QTHLL(3),
     $                QTHER(3),QNSH,QNSG,QNH,QTHSB,QTHSW,EPS
                                     !EPSILON IS USED WHEN Q=THRESHOLD
      INTEGER LOCMH
      SAVE/THRESH/
!
      COMMON /BSG/GISA(31),MSQISA(3),MSLISA(3),MSUISA(3),MSDISA(3),
     &            MSEISA(3),MRNISA(3),YNFRZ(3,3),MNFRZ(3,3),TNFRZ(3,3),
     &            RTISA,RBISA,RLISA
c     MSxDEC(i) - decoupling scale of i-th generation of type x sfermion
c     MRNDEC(i) - decoupling scale of i-th RH neutrino
      REAL*8 GISA,MSQISA,MSLISA,MSUISA,MSDISA,MSEISA,MRNISA,
     &       YNFRZ,MNFRZ,TNFRZ
      REAL RTISA,RBISA,RLISA
      SAVE /BSG/
!
      COMMON/SQEIG/MQVE,MQVA,MUPVE,MUPVA,MDVE,MDVA,MLVE,MLVA,MEVE,MEVA
      DOUBLE COMPLEX MQVE(3,3),MUPVE(3,3),MDVE(3,3),MLVE(3,3),MEVE(3,3)
     $               ,MQVA(3),MUPVA(3),MDVA(3),MLVA(3),MEVA(3)
      SAVE/SQEIG/
!
      COMMON/SFMFRZ/MQSAV,MUPSAV,MDSAV,MLSAV,MESAV
      DOUBLE COMPLEX MQSAV(3,4,3),MUPSAV(3,4,3),MDSAV(3,4,3)
      DOUBLE COMPLEX MLSAV(3,4,3),MESAV(3,4,3)
      SAVE/SFMFRZ/
!
      COMMON/MYSUGRA/M0,M12,A0,TANB,SIGNMU,MT
      DOUBLE PRECISION M0,M12,A0,TANB,SIGNMU,MT
      SAVE/MYSUGRA/
!
      COMMON/EWSBSAV/CSHSH,CSHLH,CLHLH,CULSHUR,CDLSHDR,CELSHER,CULLHUR,
     $               CDLLHDR,CELLHER
      DOUBLE COMPLEX CSHSH,CSHLH,CLHLH,CULSHUR(3,3),CDLSHDR(3,3),
     $               CELSHER(3,3),CULLHUR(3,3),CDLLHDR(3,3),CELLHER(3,3)
      SAVE/EWSBSAV/
!
      COMMON/RGEIN/MHIGH,PHASEMU,ACC,COMP,SUG,UNI
      DOUBLE PRECISION MHIGH,PHASEMU
      INTEGER ACC,COMP,SUG,UNI
      SAVE/RGEIN/
!
      REAL G0(601),QIN,MUIN,SIG1,SIG2 !The variables used by SUGEFFFL
      DOUBLE COMPLEX BMHM1SQ,BMHM2SQ,BMHB,BMHAU(3,3),BMHAD(3,3)
      DOUBLE COMPLEX BMHAE(3,3),BMHMTSFU(3,3),BMHMTSFD(3,3)
      DOUBLE COMPLEX BMHMTSFE(3,3),GSAV(601),MDIFF
      DOUBLE COMPLEX BMHM1SQNEW,BMHM2SQNEW
      DOUBLE PRECISION QMSQ(3),UMSQ(3),DMSQ(3),LMSQ(3),EMSQ(3)
      DOUBLE PRECISION SINB,COSB,VU,VD,VUSQ,VDSQ,VSM,MUMSQ
      INTEGER I,J
!
      SINB=DSQRT(TANB**2/(1+TANB**2))
      COSB=SINB/TANB
      VU=DBLE(G(110))
      VD=DBLE(G(111))
      VSM=DBLE(G(428))
      VUSQ=VU**2
      VDSQ=VD**2
!
!Save the original G entries
!
      DO I=1,601
        GSAV(I)=G(I)
      END DO
!
!Solve the simultaneous equations
!
      IF(RGEMS.LT.QNH)THEN
        CSHSH=G(427)
        DO I=1,3
          DO J=1,3
            CULSHUR(I,J)=G(399+(I-1)*3+J)
            CDLSHDR(I,J)=G(408+(I-1)*3+J)
            CELSHER(I,J)=G(417+(I-1)*3+J)
          END DO
        END DO
        BMHM1SQ=COSB**2*CSHSH-2.D0*SINB*COSB*CSHLH+SINB**2*CLHLH
        BMHM2SQ=SINB**2*CSHSH+2.D0*SINB*COSB*CSHLH+COSB**2*CLHLH
        BMHB=SINB*COSB*(-CSHSH+CLHLH)+(SINB**2-COSB**2)*CSHLH
        DO I=1,3
          DO J=1,3
            BMHAU(I,J)=SINB*CULSHUR(I,J)+COSB*CULLHUR(I,J)
            BMHAD(I,J)=COSB*CDLSHDR(I,J)-SINB*CDLLHDR(I,J)
            BMHAE(I,J)=COSB*CELSHER(I,J)-SINB*CELLHER(I,J)
            BMHMTSFU(I,J)=-COSB*CULSHUR(I,J)+SINB*CULLHUR(I,J)
            BMHMTSFD(I,J)=-SINB*CDLSHDR(I,J)-COSB*CDLLHDR(I,J)
            BMHMTSFE(I,J)=-SINB*CELSHER(I,J)-COSB*CELLHER(I,J)
!
            G(33+(I-1)*3+J)=BMHAU(I,J)
            G(42+(I-1)*3+J)=BMHAD(I,J)
            G(51+(I-1)*3+J)=BMHAE(I,J)
            G(429+(I-1)*3+J)=BMHMTSFU(I,J)
            G(438+(I-1)*3+J)=BMHMTSFD(I,J)
            G(447+(I-1)*3+J)=BMHMTSFE(I,J)
          END DO
        END DO
      END IF
!
!Find the squark mass eigenvalues in the quark mass basis
!
      CALL MASSSQM(G)
      CALL DIAGSQM(QMSQ,UMSQ,DMSQ,LMSQ,EMSQ)
!
!Apply the finite shifts to the third generation Yukawas
!
      CALL ROTBACK(1)
      G(12)=G(12)*(1.D0-DBLE(RTISA))
      G(21)=G(21)*(1.D0-DBLE(RBISA))
      G(30)=G(30)*(1.D0-DBLE(RLISA))
      IF(RGEMS.LE.QNH)THEN
        G(120)=G(120)*(1.D0-DBLE(RTISA))
        G(129)=G(129)*(1.D0-DBLE(RBISA))
        G(138)=G(138)*(1.D0-DBLE(RLISA))
      END IF
!
!Set up and call the 1-loop radiative corrections to the scalar potential
!Since Isajet has differing sign conventions, we flip the signs
!of \mu and the a-parameters.
!
      QIN=REAL(RGEMS)
      MUIN=-REAL(REAL(G(108)))
      DO I=1,601
        G0(I)=REAL(REAL(G(I))) !The important quantities have Im()~0
      END DO
      DO I=1,27
        G0(33+I)=-G0(33+I)
        G0(323+I)=-G0(323+I)
        G0(399+I)=-G0(399+I)
        G0(429+I)=-G0(429+I)
      END DO
!
!Fix the values of the sfermion masses
!
      DO I=1,3
        G0(62+(I-1)*3+I)=REAL(QMSQ(I))
        G0(80+(I-1)*3+I)=REAL(UMSQ(I))
        G0(89+(I-1)*3+I)=REAL(DMSQ(I))
        G0(71+(I-1)*3+I)=REAL(LMSQ(I))
        G0(98+(I-1)*3+I)=REAL(EMSQ(I))
      END DO
!
      IF(RGEMS.LT.QNH)THEN
        G0(110)=REAL(SINB*VSM)
        G0(111)=REAL(COSB*VSM)
      END IF
!
!NOTE: Following Isajet's notation, SIG1 is the correction
!      to the down-type Higgs field.
!
      CALL SUGEFFFL(QIN,G0,MUIN,SIG1,SIG2)
!
!Now use EWSB conditions to fix b and (m^2_1+m^2_2).
!Also fix b and mu for the MSSM running couplings.
!
      IF(RGEMS.LT.QNH)THEN
        MDIFF=BMHM1SQ-BMHM2SQ
        BMHM1SQNEW=((-(MDIFF+SIG1-SIG2)*1.D0/2.D0*(1.D0+TANB**2)+1.D0/4.D0
     $          *(G(1)**2+G(2)**2)*(SINB**2-COSB**2)/COSB**2*VSM**2)
     $          /(1.D0-TANB**2))+1.D0/2.D0*(MDIFF+SIG1-SIG2)-SIG1
        BMHM2SQNEW=BMHM1SQNEW-MDIFF
!
!Set the new value of b with the old m1 and m2
!
        BMHB=(BMHM1SQ+BMHM2SQ+SIG1+SIG2)*SINB*COSB
!
!Now replace the values of m1 and m2
!
        BMHM1SQ=BMHM1SQNEW
        BMHM2SQ=BMHM2SQNEW
!
!The MSSM mu parameter is fixed next
!
        MUMSQ=((DBLE(G(351))+SIG2)*TANB**2-DBLE(G(352))-SIG1)
     $       /(1.D0-TANB**2)-1.D0/4.D0*(DBLE(G(291))**2+DBLE(G(292))**2)
     $       *VSM**2
        IF(MUMSQ.LT.0.D0)MUMSQ=MZ**2
!
!Now reset the MSSM mu and b
!
        G(399)=(G(351)+G(352)+SIG1+SIG2+2.D0*G(398)**2)*SINB*COSB
        G(398)=DSQRT(MUMSQ)*EXP((0.D0,1.D0)*PHASEMU)
!
!Reset the frozen parameters
!
        CSHSH=SINB**2*BMHM2SQ+COSB**2*BMHM1SQ-2.D0*SINB*COSB*BMHB
        CSHLH=SINB*COSB*BMHM2SQ-SINB*COSB*BMHM1SQ+(SINB**2-COSB**2)*BMHB
        CLHLH=COSB**2*BMHM2SQ+SINB**2*BMHM1SQ+2.D0*SINB*COSB*BMHB
!
!Reset the running parameter
!
        G(427)=CSHSH
      ELSE
!
!Carry out the greater than m_H conditions in a similar manner
!
        MDIFF=G(62)-G(61)
        BMHM1SQNEW=((-(MDIFF+SIG1-SIG2)*1.D0/2.D0*(1.D0+VUSQ/VDSQ)
     $       +1.D0/4.D0*(G(1)**2+G(2)**2)*(VUSQ-VDSQ)*(VUSQ+VDSQ)/VDSQ)
     $       /(1.D0-VUSQ/VDSQ))+1.D0/2.D0*(MDIFF+SIG1-SIG2)-SIG1
        BMHM2SQNEW=BMHM1SQNEW-MDIFF
        G(109)=(G(62)+G(61)+SIG1+SIG2)*VU*VD/(VUSQ+VDSQ)
        G(62)=BMHM1SQNEW
        G(61)=BMHM2SQNEW
        MUMSQ=((DBLE(G(351))+SIG2)*VUSQ/VDSQ-DBLE(G(352))-SIG1)
     $         /(1.D0-VUSQ/VDSQ)+1.D0/4.D0*(DBLE(G(291))**2
     $         +DBLE(G(292))**2)*(VUSQ+VDSQ)
        IF(MUMSQ.LT.0.D0)MUMSQ=MZ**2
        G(399)=(G(351)+G(352)+SIG1+SIG2+2.D0*G(398)**2)
     $                                              *VU*VD/(VUSQ+VDSQ)
        G(398)=DSQRT(MUMSQ)*EXP((0.D0,1.D0)*PHASEMU)
      END IF
      CALL ROTATE(1)
!
!Replace the original G entries, excluding the fixed couplings
!
      DO I=1,3
        G(I)=GSAV(I)
      END DO
      DO I=31,60
        G(I)=GSAV(I)
      END DO
      DO I=63,108
        G(I)=GSAV(I)
      END DO
      DO I=110,111
        G(I)=GSAV(I)
      END DO
      DO I=139,397
        G(I)=GSAV(I)
      END DO
      DO I=400,426
        G(I)=GSAV(I)
      END DO
      DO I=428,601
        G(I)=GSAV(I)
      END DO
!
      RETURN
      END
