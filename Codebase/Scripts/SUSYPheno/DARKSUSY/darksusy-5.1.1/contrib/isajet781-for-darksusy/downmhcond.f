!
      SUBROUTINE DOWNMHCOND
!
!Purpose: Apply the matching conditions at m_H when running down.
!         If m_H is greater than m_SUSY, the values of decoupled
!         operators are frozen here for use later.
!
      IMPLICIT NONE
!
      COMMON/COUPLINGS/G,DG
      DOUBLE COMPLEX G(601)
      DOUBLE PRECISION DG(601)
      SAVE/COUPLINGS/
!
      COMMON/SMSAVED/KM,MWEAK,MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      DOUBLE COMPLEX KM(3,3)
      DOUBLE PRECISION MWEAK(6),MZ,MW,ALPHAEM,ALPHASMSB,XWMSB
      SAVE/SMSAVED/
!
      COMMON/RGEMS/VEVMH,RGEMS,RGEMU
      DOUBLE COMPLEX VEVMH
      DOUBLE PRECISION RGEMS,RGEMU
      SAVE/RGEMS/
!
      COMMON/MYSUGRA/M0,M12,A0,TANB,SIGNMU,MT
      DOUBLE PRECISION M0,M12,A0,TANB,SIGNMU,MT
      SAVE/MYSUGRA/
!
      COMMON/THRESH/QTHSORT,QTHQL,QTHUR,QTHDR,QTHLL,QTHER,QNSH,QNSG,
     $              QNH,QTHSB,QTHSW,EPS,LOCMH
      DOUBLE PRECISION QTHSORT(20),QTHQL(3),QTHUR(3),QTHDR(3),QTHLL(3),
     $                QTHER(3),QNSH,QNSG,QNH,QTHSB,QTHSW,EPS
                                     !EPSILON IS USED WHEN Q=THRESHOLD
      INTEGER LOCMH
      SAVE/THRESH/
!
      COMMON/EWSBSAV/CSHSH,CSHLH,CLHLH,CULSHUR,CDLSHDR,CELSHER,CULLHUR,
     $               CDLLHDR,CELLHER
      DOUBLE COMPLEX CSHSH,CSHLH,CLHLH,CULSHUR(3,3),CDLSHDR(3,3),
     $               CELSHER(3,3),CULLHUR(3,3),CDLLHDR(3,3),CELLHER(3,3)
      SAVE/EWSBSAV/
!
      DOUBLE PRECISION SINB,COSB
      DOUBLE COMPLEX B,VU,VD
      INTEGER I,J
!
      SINB=DSQRT(TANB**2/(1+TANB**2))
      COSB=SINB/TANB
      B=G(109)
      VU=G(110)
      VD=G(111)
!
      G(429)=(3.D0/5.D0*(G(1))**2+(G(2))**2)/4.D0
     $                                  *((1+TANB**2)/(1-TANB**2))**2
      G(428)=DSQRT(DBLE(VU*CONJG(VU))+DBLE(VD*CONJG(VD)))
      DO I=1,3
        DO J=1,3
          G(111+(I-1)*3+J)=SINB*G(3+(I-1)*3+J)
          G(120+(I-1)*3+J)=COSB*G(12+(I-1)*3+J)
          G(129+(I-1)*3+J)=COSB*G(21+(I-1)*3+J)
          G(399+(I-1)*3+J)=SINB*G(33+(I-1)*3+J)-COSB*G(429+(I-1)*3+J)
          G(408+(I-1)*3+J)=COSB*G(42+(I-1)*3+J)-SINB*G(438+(I-1)*3+J)
          G(417+(I-1)*3+J)=COSB*G(51+(I-1)*3+J)-SINB*G(447+(I-1)*3+J)
        END DO
      END DO
      G(287)=SINB*G(184)
      G(288)=COSB*G(185)
      G(289)=SINB*G(204)
      G(290)=COSB*G(205)
      G(427)=SINB**2*G(61)+COSB**2*G(62)-2.D0*SINB*COSB*B
!
!Save the combinations necessary for the EWSB conditions if m_H > m_SUSY
!
      IF(RGEMS.LT.QNH)THEN
        CSHSH=SINB**2*G(61)+COSB**2*G(62)-2.D0*SINB*COSB*B
        CSHLH=SINB*COSB*G(61)-SINB*COSB*G(62)+(SINB**2-COSB**2)*B
        CLHLH=COSB**2*G(61)+SINB**2*G(62)+2.D0*SINB*COSB*B
        DO I=1,3
          DO J=1,3
            CULSHUR(I,J)=SINB*G(33+(I-1)*3+J)-COSB*G(429+(I-1)*3+J)
            CDLSHDR(I,J)=COSB*G(42+(I-1)*3+J)-SINB*G(438+(I-1)*3+J)
            CELSHER(I,J)=COSB*G(51+(I-1)*3+J)-SINB*G(447+(I-1)*3+J)
            CULLHUR(I,J)=COSB*G(33+(I-1)*3+J)+SINB*G(429+(I-1)*3+J)
            CDLLHDR(I,J)=-SINB*G(42+(I-1)*3+J)-COSB*G(438+(I-1)*3+J)
            CELLHER(I,J)=-SINB*G(51+(I-1)*3+J)-COSB*G(447+(I-1)*3+J)
          END DO
        END DO
      END IF
!
      RETURN
      END
