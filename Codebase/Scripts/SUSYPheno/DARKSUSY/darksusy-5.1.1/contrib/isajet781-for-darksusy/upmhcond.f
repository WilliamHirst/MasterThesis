!
      SUBROUTINE UPMHCOND(GCON)
!
!Purpose: Apply the matching conditions at m_H when running up
!
      IMPLICIT NONE
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
      DOUBLE COMPLEX GCON(215)
      DOUBLE PRECISION SINB,COSB
      INTEGER I
!
      SINB=DSQRT(TANB**2/(1+TANB**2))
      COSB=SINB/TANB
!
!This uses the input value of tanb - i.e. tanb at m_H.
!
      VEVMH=GCON(61)
      GCON(33)=DCMPLX(DSQRT((DBLE(VEVMH))**2/(1.D0+TANB**2)))
      GCON(32)=TANB*GCON(33)
!
      DO I=4,12
        GCON(I)=GCON(I+30)/SINB
      END DO
      DO I=13,30
        GCON(I)=GCON(I+30)/COSB
      END DO
!
      RETURN
      END
