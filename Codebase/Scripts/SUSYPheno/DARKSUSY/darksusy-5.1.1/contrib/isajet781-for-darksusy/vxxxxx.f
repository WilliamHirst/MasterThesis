C
C
C       Subroutine returns the value of evaluated
C       helicity basis boson polarisation wavefunction.
C       Replaces the HELAS routine VXXXXX
C
C       Adam Duff,  1992 September 3
C       <duff@phenom.physics.wisc.edu>
C
      SUBROUTINE VXXXXX(P,VMASS,NHEL,NSV,VC)
C
C          P            IN: BOSON FOUR MOMENTUM
C          VMASS        IN: BOSON MASS
C          NHEL         IN: BOSON HELICITY
C          NSV          IN: INCOMING (-1) OR OUTGOING (+1)
C          VC           OUT: BOSON WAVEFUNCTION
C
C declare input/output variables
C
      IMPLICIT NONE
      COMPLEX*16 VC(6)
      INTEGER*4 NHEL, NSV
      REAL*8 P(0:3), VMASS
C
C declare local variables
C
      REAL*8 RXZERO, RXONE, RXTWO
      PARAMETER( RXZERO=0.0D0, RXONE=1.0D0, RXTWO=2.0D0 )
      REAL*8 PLAT, PABS, RS2, RPLAT, RPABS, RDEN
      COMPLEX*16 CXZERO
      LOGICAL FIRST
      SAVE CXZERO,FIRST
      DATA FIRST/.TRUE./
C
C          Fix compilation with g77
      IF(FIRST) THEN
        FIRST=.FALSE.
        CXZERO=DCMPLX( RXZERO, RXZERO )
      ENDIF
C
C define internal/external momenta
C
      IF ( NSV**2 .NE. 1 ) THEN
         STOP 'VXXXXX:  NSV IS NOT ONE OF -1, +1'
      END IF
C
      RS2 = SQRT( RXONE / RXTWO )
      VC(5) = DCMPLX( P(0), P(3) ) * NSV
      VC(6) = DCMPLX( P(1), P(2) ) * NSV
      PLAT = SQRT( P(1)**2 + P(2)**2 )
      PABS = SQRT( P(1)**2 + P(2)**2 + P(3)**2 )
C
C calculate polarisation four vectors
C
      IF ( NHEL**2 .EQ. 1 ) THEN
         IF ( (PABS .EQ. RXZERO) .OR. (PLAT .EQ. RXZERO) ) THEN
            VC(1) = CXZERO
            VC(2) = DCMPLX( -NHEL * RS2 * DSIGN( RXONE, P(3) ), RXZERO )
            VC(3) = DCMPLX( RXZERO, NSV * RS2 )
            VC(4) = CXZERO
         ELSE
            RPLAT = RXONE / PLAT
            RPABS = RXONE / PABS
            VC(1) = CXZERO
            VC(2) = DCMPLX( -NHEL * RS2 * RPABS * RPLAT * P(1) * P(3),
     &                     -NSV * RS2 * RPLAT * P(2) )
            VC(3) = DCMPLX( -NHEL * RS2 * RPABS * RPLAT * P(2) * P(3),
     &                     NSV * RS2 * RPLAT * P(1) )
            VC(4) = DCMPLX( NHEL * RS2 * RPABS * PLAT,
     &                     RXZERO )
         END IF
      ELSE IF ( NHEL .EQ. 0 ) THEN
         IF ( VMASS .GT. RXZERO ) THEN
            IF ( PABS .EQ. RXZERO ) THEN
               VC(1) = CXZERO
               VC(2) = CXZERO
               VC(3) = CXZERO
               VC(4) = DCMPLX( RXONE, RXZERO )
            ELSE
               RDEN = P(0) / ( VMASS * PABS )
               VC(1) = DCMPLX( PABS / VMASS, RXZERO )
               VC(2) = DCMPLX( RDEN * P(1), RXZERO )
               VC(3) = DCMPLX( RDEN * P(2), RXZERO )
               VC(4) = DCMPLX( RDEN * P(3), RXZERO )
            END IF
         ELSE
            STOP  'VXXXXX: NHEL = 0 IS ONLY FOR MASSIVE BOSONS'
         END IF
      ELSE IF ( NHEL .EQ. 4 ) THEN
         IF ( VMASS .GT. RXZERO ) THEN
            RDEN = RXONE / VMASS
            VC(1) = DCMPLX( RDEN * P(0), RXZERO )
            VC(2) = DCMPLX( RDEN * P(1), RXZERO )
            VC(3) = DCMPLX( RDEN * P(2), RXZERO )
            VC(4) = DCMPLX( RDEN * P(3), RXZERO )
         ELSEIF (VMASS .EQ. RXZERO) THEN
            RDEN = RXONE / P(0)
            VC(1) = DCMPLX( RDEN * P(0), RXZERO )
            VC(2) = DCMPLX( RDEN * P(1), RXZERO )
            VC(3) = DCMPLX( RDEN * P(2), RXZERO )
            VC(4) = DCMPLX( RDEN * P(3), RXZERO )
         ELSE
            STOP 'VXXXXX: NHEL = 4 IS ONLY FOR M>=0'
         END IF
      ELSE
         STOP 'VXXXXX:  NHEL IS NOT ONE OF -1, 0, 1 OR 4'
      END IF
C
C done
C
      RETURN
      END
