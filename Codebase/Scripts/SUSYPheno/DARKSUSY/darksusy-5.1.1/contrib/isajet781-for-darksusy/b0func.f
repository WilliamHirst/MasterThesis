C
C**********************************************************************
C FUNCTION B0 FROM PIERCE, PAG. 12
C
      FUNCTION B0FUNC(SCALE,M1,M2)
      REAL B0FUNC, SCALE, M1, M2, MC2, ML2
C
      IF(M1.GT.M2)THEN
         MC2 = M1*M1
         ML2 = M2*M2
      ELSE
         MC2 = M2*M2
         ML2 = M1*M1
      END IF
C
      B0FUNC = -LOG(MC2/SCALE**2) + 1. + ML2/(ML2-MC2)*LOG(MC2/ML2)
C
      RETURN
      END
