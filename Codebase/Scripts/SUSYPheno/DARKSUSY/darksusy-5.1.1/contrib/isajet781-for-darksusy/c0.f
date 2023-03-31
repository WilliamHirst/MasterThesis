C
C
C**********************************************************************      
C Given by Pierce, (Eq. C19)
C
C----------------------------------------------------------------------
      FUNCTION C0(M1,M2,M3)
C
      REAL C0, M1, M2, M3
      REAL*8 X, Y, Z, X2, Y2, Z2, FRACXY, FRACXZ, FRACYZ
C
      X = M1*M1
      Y = M2*M2
      Z = M3*M3
C
      FRACXY = ABS(X-Y)/(X+Y)
      FRACXZ = ABS(X-Z)/(X+Z)
      FRACYZ = ABS(Y-Z)/(Y+Z)
      IF(FRACXY.LT.1.D-4) Y = X
      IF(FRACXZ.LT.1.D-4) Z = X
      IF(FRACYZ.LT.1.D-4) Y = Z
C
      IF(X.NE.0.AND.Y.NE.0.AND.Z.NE.0)THEN
         IF(X.NE.Y.AND.X.NE.Z.AND.Y.NE.Z)THEN
            C0 = 1.D0/(Y-Z)*(Y/(X-Y)*DLOG(Y/X)-Z/(X-Z)*DLOG(Z/X))
         ELSE
            IF(X.EQ.Y.AND.X.NE.Z)THEN
               C0 = -(1.D0 - Z/(Y-Z)*DLOG(Y/Z))/(Y-Z)
            ELSE
               IF(X.EQ.Y.AND.X.EQ.Z)THEN
                  C0 = -1.D0/2.D0/Z
               END IF
            END IF
         END IF
         IF(X.NE.Y)THEN
            IF(X.EQ.Z)THEN
               C0 = -(1.D0 - Y/(X-Y)*DLOG(X/Y))/(X-Y)
            ELSE
               IF(Y.EQ.Z)THEN
                  C0 = -(1.D0 - X/(Y-X)*DLOG(Y/X))/(Y-X)
               END IF
            END IF
         END IF
      ELSE
         WRITE(*,*) 'C0 COM ARGUMENTO NULO!' 
      END IF
      RETURN
      END
