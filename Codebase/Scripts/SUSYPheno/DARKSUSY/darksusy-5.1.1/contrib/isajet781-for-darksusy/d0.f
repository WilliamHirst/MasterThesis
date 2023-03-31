C
C----------------------------------------------------------------------
      FUNCTION D0(M1,M2,M3,M4)
      REAL M1, M2, M3, M4, X, Y, Z, W, C0, FRACXY, FRACXZ, FRACXW
      INTEGER FLAG
C
      X = M1*M1
      Y = M2*M2
      Z = M3*M3
      W = M4*M4
      FRACXY = (X-Y)/(X+Y)
      FRACXZ = (X-Z)/(X+Z)
      FRACXW = (X-W)/(X+W)
C
      IF(FRACXY.LT.1.E-2) Y = X
      IF(FRACXZ.LT.1.E-2) Z = X
      IF(FRACXW.LT.1.E-2) W = X
C
      FLAG = 0
      IF(X.NE.0.AND.Y.NE.0.AND.Z.NE.0.AND.W.NE.0)THEN
         IF(X.NE.Y)THEN
            D0 = (C0(X,Z,W) - C0(Y,Z,W))/(X-Y)
            FLAG = 1
         END IF
         IF(X.EQ.Y.AND.X.NE.Z.AND.X.NE.W)THEN
            D0 = -(-(Y**2-Z*W)/(Y-Z)**2/(Y-W)**2*LOG(Y/W) +
     &           Z/(Y-Z)**2/(Z-W)*LOG(Z/W) + 1./(Y-Z)/(Y-W))
            FLAG = 1
         END IF
         IF(X.EQ.Y.AND.X.EQ.Z.AND.X.NE.W)THEN
            D0 = -(1./2./Y - 1./(Y-W) + W/(Y-W)**2*LOG(Y/W))/(Y-W)
            FLAG = 1
         END IF
         IF(X.EQ.Y.AND.X.NE.Z.AND.X.EQ.W)THEN
            D0 = -(1./(X-Z)*(-1./2./X - Z/X/(X-Z)) + 
     &           Z/(X-Z)**3*LOG(X/Z))
            FLAG = 1
         END IF
         IF(X.EQ.Y.AND.X.EQ.Z.AND.X.EQ.W)THEN
            D0 = 1./6./W**2
            FLAG = 1
         END IF
      END IF
C
      IF(FLAG.EQ.0) WRITE(*,*) 'D0: ARGUMENTOS ESTRANHOS'
      RETURN
      END
