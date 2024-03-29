C----------------------------------------------------------------------
      SUBROUTINE TQLEIG(NM,N,D,E,Z,IERR)
C----------------------------------------------------------------------
C     Computes eigenvalues and eigenvectors of the symmetric tridiagonal 
C     matrix using QL-algorithm.
C     
C     This is a DOUBLE PRECISION version of the TQL2 subroutine (F220)
C     from CERN program library.
C
C     Created: 02/23/07 by Azar Mustafayev
C
      IMPLICIT NONE
      INTEGER I,J,K,L,M,N,II,NM,MML,IERR
      REAL*8 D(N),E(N),Z(NM,N)
      REAL*8 B,C,F,G,H,P,R,S,MACHEP
      
      MACHEP=2.D0**(-23)
      MACHEP=2.D0**(-47)
      IERR = 0
      
      IF (N .EQ. 1) GO TO 1001
      DO 100 I = 2, N
        E(I-1) = E(I)
100   CONTINUE
      F = 0.D0
      B = 0.D0
      E(N) = 0.D0
      DO 240 L = 1, N
         J = 0
         H = MACHEP * (DABS(D(L)) + DABS(E(L)))
         IF (B .LT. H) B = H
         DO 110 M = L, N
            IF (DABS(E(M)) .LE. B) GO TO 120
110      CONTINUE
120      IF (M .EQ. L) GO TO 220
130      IF (J .EQ. 30) GO TO 1000
         J = J + 1
         P = (D(L+1) - D(L)) / (2.D0 * E(L))
         R = DSQRT(P*P+1.D0)
         H = D(L) - E(L) / (P + DSIGN(R,P))
         DO 140 I = L, N
           D(I) = D(I) - H
140      CONTINUE
         F = F + H
         P = D(M)
         C = 1.D0
         S = 0.D0
         MML = M - L
         DO 200 II = 1, MML
            I = M - II
            G = C * E(I)
            H = C * P
            IF (DABS(P) .LT. DABS(E(I))) GO TO 150
            C = E(I) / P
            R = DSQRT(C*C+1.D0)
            E(I+1) = S * P * R
            S = C / R
            C = 1.D0 / R
            GO TO 160
150         C = P / E(I)
            R = DSQRT(C*C+1.D0)
            E(I+1) = S * E(I) * R
            S = 1.D0 / R
            C = C * S
160         P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
180         CONTINUE
200      CONTINUE
         E(L) = S * P
         D(L) = C * P
         IF (DABS(E(L)) .GT. B) GO TO 130
220      D(L) = D(L) + F
240   CONTINUE
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
260      CONTINUE
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
280      CONTINUE
300   CONTINUE
      GO TO 1001
1000  IERR = L

1001  RETURN
      END
