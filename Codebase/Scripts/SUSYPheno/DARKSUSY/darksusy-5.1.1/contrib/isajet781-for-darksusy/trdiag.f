C----------------------------------------------------------------------
      SUBROUTINE TRDIAG(NM,N,A,D,E,Z)
C----------------------------------------------------------------------
C     Reduces real symmetric matrix A to trigiagonal form.
C     This is a DOUBLE PRECISION version of the TRED2 subroutine (F220)
C     from CERN program library.
C
C     Created: 02/23/07 by Azar Mustafayev
C
      IMPLICIT NONE
      INTEGER I,J,K,L,N,II,NM,JP1
      REAL*8 A(NM,N),D(N),E(N),Z(NM,N)
      REAL*8 F,G,H,HH,SCALE
      
      DO 100 I = 1, N
         DO 100 J = 1, I
            Z(I,J) = A(I,J)
100   CONTINUE
      IF (N .EQ. 1) GO TO 320
      DO 300 II = 2, N
         I = N + 2 - II
         L = I - 1
         H = 0.D0
         SCALE = 0.D0
         IF (L .LT. 2) GO TO 130
         DO 120 K = 1, L
           SCALE = SCALE + DABS(Z(I,K))
120      CONTINUE
         IF (SCALE .NE. 0.D0) GO TO 140
130      E(I) = Z(I,L)
         GO TO 290
140      DO 150 K = 1, L
            Z(I,K) = Z(I,K) / SCALE
            H = H + Z(I,K) * Z(I,K)
150      CONTINUE
         F = Z(I,L)
         G = -DSIGN(DSQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         Z(I,L) = F - G
         F = 0.D0
         DO 240 J = 1, L
            Z(J,I) = Z(I,J) / (SCALE * H)
            G = 0.D0
            DO 180 K = 1, J
              G = G + Z(J,K) * Z(I,K)
180         CONTINUE
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
            DO 200 K = JP1, L
              G = G + Z(K,J) * Z(I,K)
200         CONTINUE
220         E(J) = G / H
            F = F + E(J) * Z(I,J)
240      CONTINUE
         HH = F / (H + H)
         DO 260 J = 1, L
            F = Z(I,J)
            G = E(J) - HH * F
            E(J) = G
            DO 260 K = 1, J
               Z(J,K) = Z(J,K) - F * E(K) - G * Z(I,K)
260      CONTINUE
         DO 280 K = 1, L
           Z(I,K) = SCALE * Z(I,K)
280      CONTINUE
290        D(I) = H
300   CONTINUE
320   D(1) = 0.D0
      E(1) = 0.D0
      DO 500 I = 1, N
         L = I - 1
         IF (D(I) .EQ. 0.D0) GO TO 380
         DO 360 J = 1, L
            G = 0.D0
            DO 340 K = 1, L
              G = G + Z(I,K) * Z(K,J)
340         CONTINUE
            DO 360 K = 1, L
               Z(K,J) = Z(K,J) - G * Z(K,I)
360      CONTINUE
380      D(I) = Z(I,I)
         Z(I,I) = 1.D0
         IF (L .LT. 1) GO TO 500
         DO 400 J = 1, L
            Z(I,J) = 0.D0
            Z(J,I) = 0.D0
400      CONTINUE
500   CONTINUE
      
      RETURN
      END
