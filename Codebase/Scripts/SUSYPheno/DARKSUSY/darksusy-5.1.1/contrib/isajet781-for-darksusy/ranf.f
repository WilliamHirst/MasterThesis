CDECK  ID>, RANFLUX.
      REAL FUNCTION RANF()
C
C          Call RANLUX instead of 48-bit congruental generator RANF
C          Dummy RANFGT/RANFST/RANFMT
C
      DIMENSION X(1)
      CALL RANLUX(X,1)
      RANF=X(1)
      RETURN
      END
