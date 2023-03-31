C
C
C
C       Subroutine returns the desired fermion or
C       anti-fermion anti-spinor. ie., <f|
C       A replacement for the HELAS routine OXXXXX
C
C       Adam Duff,  1992 August 31
C       <duff@phenom.physics.wisc.edu>
C
      SUBROUTINE OXXXXX(P,FMASS,NHEL,NSF,FO)
C
C          P            IN: FOUR VECTOR MOMENTUM
C          FMASS        IN: FERMION MASS
C          NHEL         IN: ANTI-SPINOR HELICITY, -1 OR 1
C          NSF          IN: -1=ANTIFERMION, 1=FERMION
C          FO           OUT: FERMION WAVEFUNCTION
C
C declare input/output variables
C
      IMPLICIT NONE
      COMPLEX*16 FO(6)
      INTEGER*4 NHEL, NSF
      REAL*8 P(0:3), FMASS
C
C declare local variables
C
      REAL*8 RXZERO, RXONE, RXTWO
      PARAMETER( RXZERO=0.0D0, RXONE=1.0D0, RXTWO=2.0D0 )
      REAL*8 PLAT, PABS, OMEGAP, OMEGAM, RS2PA, SPAZ
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
C define kinematic parameters
C
      FO(5) = DCMPLX( P(0), P(3) ) * NSF
      FO(6) = DCMPLX( P(1), P(2) ) * NSF
      PLAT = SQRT( P(1)**2 + P(2)**2 )
      PABS = SQRT( P(1)**2 + P(2)**2 + P(3)**2 )
      OMEGAP = SQRT( P(0) + PABS )
C
C do massive fermion case
C
      IF ( FMASS .NE. RXZERO ) THEN
         OMEGAM = FMASS / OMEGAP
         IF ( NSF .EQ. 1 ) THEN
            IF ( NHEL .EQ. 1 ) THEN
               IF ( P(3) .GE. RXZERO ) THEN
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FO(1) = DCMPLX( OMEGAP, RXZERO )
                     FO(2) = CXZERO
                     FO(3) = DCMPLX( OMEGAM, RXZERO )
                     FO(4) = CXZERO
                  ELSE
                     RS2PA = RXONE / SQRT( RXTWO * PABS )
                     SPAZ = SQRT( PABS + P(3) )
                     FO(1) = OMEGAP * RS2PA
     &                     * DCMPLX( SPAZ, RXZERO )
                     FO(2) = OMEGAP * RS2PA / SPAZ
     &                     * DCMPLX( P(1), -P(2) )
                     FO(3) = OMEGAM * RS2PA
     &                     * DCMPLX( SPAZ, RXZERO )
                     FO(4) = OMEGAM * RS2PA / SPAZ
     &                     * DCMPLX( P(1), -P(2) )
                  END IF
               ELSE
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FO(1) = CXZERO
                     FO(2) = DCMPLX( OMEGAP, RXZERO )
                     FO(3) = CXZERO
                     FO(4) = DCMPLX( OMEGAM, RXZERO )
                  ELSE
                     RS2PA = RXONE / SQRT( RXTWO * PABS )
                     SPAZ = SQRT( PABS - P(3) )
                     FO(1) = OMEGAP * RS2PA / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                     FO(2) = OMEGAP * RS2PA * SPAZ / PLAT
     &                     * DCMPLX( P(1), -P(2) )
                     FO(3) = OMEGAM * RS2PA / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                     FO(4) = OMEGAM * RS2PA * SPAZ / PLAT
     &                     * DCMPLX( P(1), -P(2) )
                  END IF
               END IF
            ELSE IF ( NHEL .EQ. -1 ) THEN
               IF ( P(3) .GE. RXZERO ) THEN
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FO(1) = CXZERO
                     FO(2) = DCMPLX( OMEGAM, RXZERO )
                     FO(3) = CXZERO
                     FO(4) = DCMPLX( OMEGAP, RXZERO )
                  ELSE
                     RS2PA = RXONE / SQRT( RXTWO * PABS )
                     SPAZ = SQRT( PABS + P(3) )
                     FO(1) = OMEGAM * RS2PA / SPAZ
     &                     * DCMPLX( -P(1), -P(2) )
                     FO(2) = OMEGAM * RS2PA
     &                     * DCMPLX( SPAZ, RXZERO )
                     FO(3) = OMEGAP * RS2PA / SPAZ
     &                     * DCMPLX( -P(1), -P(2) )
                     FO(4) = OMEGAP * RS2PA
     &                     * DCMPLX( SPAZ, RXZERO )
                  END IF
               ELSE
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FO(1) = DCMPLX( -OMEGAM, RXZERO )
                     FO(2) = CXZERO
                     FO(3) = DCMPLX( -OMEGAP, RXZERO )
                     FO(4) = CXZERO
                  ELSE
                     RS2PA = RXONE / SQRT( RXTWO * PABS )
                     SPAZ = SQRT( PABS - P(3) )
                     FO(1) = OMEGAM * RS2PA * SPAZ / PLAT
     &                     * DCMPLX( -P(1), -P(2) )
                     FO(2) = OMEGAM * RS2PA / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                     FO(3) = OMEGAP * RS2PA * SPAZ / PLAT
     &                     * DCMPLX( -P(1), -P(2) )
                     FO(4) = OMEGAP * RS2PA / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                  END IF
               END IF
            ELSE
               STOP 'OXXXXX:  FERMION HELICITY MUST BE +1,-1'
            END IF
         ELSE IF ( NSF .EQ. -1 ) THEN
            IF ( NHEL .EQ. 1 ) THEN
               IF ( P(3) .GE. RXZERO ) THEN
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FO(1) = CXZERO
                     FO(2) = DCMPLX( OMEGAM, RXZERO )
                     FO(3) = CXZERO
                     FO(4) = DCMPLX( -OMEGAP, RXZERO )
                  ELSE
                     RS2PA = RXONE / SQRT( RXTWO * PABS )
                     SPAZ = SQRT( PABS + P(3) )
                     FO(1) = OMEGAM * RS2PA / SPAZ
     &                     * DCMPLX( -P(1), -P(2) )
                     FO(2) = OMEGAM * RS2PA
     &                     * DCMPLX( SPAZ, RXZERO )
                     FO(3) = -OMEGAP * RS2PA / SPAZ
     &                     * DCMPLX( -P(1), -P(2) )
                     FO(4) = -OMEGAP * RS2PA
     &                     * DCMPLX( SPAZ, RXZERO )
                  END IF
               ELSE
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FO(1) = DCMPLX( -OMEGAM, RXZERO )
                     FO(2) = CXZERO
                     FO(3) = DCMPLX( OMEGAP, RXZERO )
                     FO(4) = CXZERO
                  ELSE
                     RS2PA = RXONE / SQRT( RXTWO * PABS )
                     SPAZ = SQRT( PABS - P(3) )
                     FO(1) = OMEGAM * RS2PA * SPAZ / PLAT
     &                     * DCMPLX( -P(1), -P(2) )
                     FO(2) = OMEGAM * RS2PA / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                     FO(3) = -OMEGAP * RS2PA * SPAZ / PLAT
     &                     * DCMPLX( -P(1), -P(2) )
                     FO(4) = -OMEGAP * RS2PA / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                  END IF
               END IF
            ELSE IF ( NHEL .EQ. -1 ) THEN
               IF ( P(3) .GE. RXZERO ) THEN
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FO(1) = DCMPLX( -OMEGAP, RXZERO )
                     FO(2) = CXZERO
                     FO(3) = DCMPLX( OMEGAM, RXZERO )
                     FO(4) = CXZERO
                  ELSE
                     RS2PA = RXONE / SQRT( RXTWO * PABS )
                     SPAZ = SQRT( PABS + P(3) )
                     FO(1) = -OMEGAP * RS2PA
     &                     * DCMPLX( SPAZ, RXZERO )
                     FO(2) = -OMEGAP * RS2PA / SPAZ
     &                     * DCMPLX( P(1), -P(2) )
                     FO(3) = OMEGAM * RS2PA
     &                     * DCMPLX( SPAZ, RXZERO )
                     FO(4) = OMEGAM * RS2PA / SPAZ
     &                     * DCMPLX( P(1), -P(2) )
                  END IF
               ELSE
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FO(1) = CXZERO
                     FO(2) = DCMPLX( -OMEGAP, RXZERO )
                     FO(3) = CXZERO
                     FO(4) = DCMPLX( OMEGAM, RXZERO )
                  ELSE
                     RS2PA = RXONE / SQRT( RXTWO * PABS )
                     SPAZ = SQRT( PABS - P(3) )
                     FO(1) = -OMEGAP * RS2PA / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                     FO(2) = -OMEGAP * RS2PA * SPAZ / PLAT
     &                     * DCMPLX( P(1), -P(2) )
                     FO(3) = OMEGAM * RS2PA / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                     FO(4) = OMEGAM * RS2PA * SPAZ / PLAT
     &                     * DCMPLX( P(1), -P(2) )
                  END IF
               END IF
            ELSE
               STOP 'OXXXXX:  FERMION HELICITY MUST BE +1,-1'
            END IF
         ELSE
            STOP 'OXXXXX:  FERMION TYPE MUST BE +1,-1'
         END IF
C
C do massless case
C
      ELSE
         IF ( NSF .EQ. 1 ) THEN
            IF ( NHEL .EQ. 1 ) THEN
               IF ( P(3) .GE. RXZERO ) THEN
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FO(1) = DCMPLX( OMEGAP, RXZERO )
                     FO(2) = CXZERO
                     FO(3) = CXZERO
                     FO(4) = CXZERO
                  ELSE
                     SPAZ = SQRT( PABS + P(3) )
                     FO(1) = DCMPLX( SPAZ, RXZERO )
                     FO(2) = RXONE / SPAZ
     &                     * DCMPLX( P(1), -P(2) )
                     FO(3) = CXZERO
                     FO(4) = CXZERO
                  END IF
               ELSE
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FO(1) = CXZERO
                     FO(2) = DCMPLX( OMEGAP, RXZERO )
                     FO(3) = CXZERO
                     FO(4) = CXZERO
                  ELSE
                     SPAZ = SQRT( PABS - P(3) )
                     FO(1) = RXONE / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                     FO(2) = SPAZ / PLAT
     &                     * DCMPLX( P(1), -P(2) )
                     FO(3) = CXZERO
                     FO(4) = CXZERO
                  END IF
               END IF
            ELSE IF ( NHEL .EQ. -1 ) THEN
               IF ( P(3) .GE. RXZERO ) THEN
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FO(1) = CXZERO
                     FO(2) = CXZERO
                     FO(3) = CXZERO
                     FO(4) = DCMPLX( OMEGAP, RXZERO )
                  ELSE
                     SPAZ = SQRT( PABS + P(3) )
                     FO(1) = CXZERO
                     FO(2) = CXZERO
                     FO(3) = RXONE / SPAZ
     &                     * DCMPLX( -P(1), -P(2) )
                     FO(4) = DCMPLX( SPAZ, RXZERO )
                  END IF
               ELSE
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FO(1) = CXZERO
                     FO(2) = CXZERO
                     FO(3) = DCMPLX( -OMEGAP, RXZERO )
                     FO(4) = CXZERO
                  ELSE
                     SPAZ = SQRT( PABS - P(3) )
                     FO(1) = CXZERO
                     FO(2) = CXZERO
                     FO(3) = SPAZ / PLAT
     &                     * DCMPLX( -P(1), -P(2) )
                     FO(4) = RXONE / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                  END IF
               END IF
            ELSE
               STOP 'OXXXXX:  FERMION HELICITY MUST BE +1,-1'
            END IF
         ELSE IF ( NSF .EQ. -1 ) THEN
            IF ( NHEL .EQ. 1 ) THEN
               IF ( P(3) .GE. RXZERO ) THEN
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FO(1) = CXZERO
                     FO(2) = CXZERO
                     FO(3) = CXZERO
                     FO(4) = DCMPLX( -OMEGAP, RXZERO )
                  ELSE
                     SPAZ = SQRT( PABS + P(3) )
                     FO(1) = CXZERO
                     FO(2) = CXZERO
                     FO(3) = -RXONE / SPAZ
     &                     * DCMPLX( -P(1), -P(2) )
                     FO(4) = DCMPLX( -SPAZ, RXZERO )
                  END IF
               ELSE
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FO(1) = CXZERO
                     FO(2) = CXZERO
                     FO(3) = DCMPLX( OMEGAP, RXZERO )
                     FO(4) = CXZERO
                  ELSE
                     SPAZ = SQRT( PABS - P(3) )
                     FO(1) = CXZERO
                     FO(2) = CXZERO
                     FO(3) = -SPAZ / PLAT
     &                     * DCMPLX( -P(1), -P(2) )
                     FO(4) = -RXONE / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                  END IF
               END IF
            ELSE IF ( NHEL .EQ. -1 ) THEN
               IF ( P(3) .GE. RXZERO ) THEN
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FO(1) = DCMPLX( -OMEGAP, RXZERO )
                     FO(2) = CXZERO
                     FO(3) = CXZERO
                     FO(4) = CXZERO
                  ELSE
                     SPAZ = SQRT( PABS + P(3) )
                     FO(1) = DCMPLX( -SPAZ, RXZERO )
                     FO(2) = -RXONE / SPAZ
     &                     * DCMPLX( P(1), -P(2) )
                     FO(3) = CXZERO
                     FO(4) = CXZERO
                  END IF
               ELSE
                  IF ( PLAT .EQ. RXZERO ) THEN
                     FO(1) = CXZERO
                     FO(2) = DCMPLX( -OMEGAP, RXZERO )
                     FO(3) = CXZERO
                     FO(4) = CXZERO
                  ELSE
                     SPAZ = SQRT( PABS - P(3) )
                     FO(1) = -RXONE / SPAZ
     &                     * DCMPLX( PLAT, RXZERO )
                     FO(2) = -SPAZ / PLAT
     &                     * DCMPLX( P(1), -P(2) )
                     FO(3) = CXZERO
                     FO(4) = CXZERO
                  END IF
               END IF
            ELSE
               STOP 'OXXXXX:  FERMION HELICITY MUST BE +1,-1'
            END IF
         ELSE
            STOP 'OXXXXX:  FERMION TYPE MUST BE +1,-1'
         END IF
      END IF
C
C done
C
      RETURN
      END
