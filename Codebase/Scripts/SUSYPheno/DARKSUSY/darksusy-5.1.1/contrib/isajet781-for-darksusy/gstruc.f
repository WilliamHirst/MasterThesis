CDECK  ID>, GSTRUC.
      REAL FUNCTION GSTRUC(X,QS)
C
C     THIS IS PHOTON PARTON DISTRIBUTION FUNCTION;
C     Improved WEIZACKER-WILLIAMS; NOTE! GSTRUC=0 FOR X<0.0000001
C     See Frixione et al. hep-ph/9310350
      IMPLICIT NONE
      REAL AL,PI,AME,QS,X,DEL
C
      AL=1./128.
      PI=4*ATAN(1.)
      AME=.511E-3
      DEL=AME**2/QS
      IF (X.LT.0.0000001.OR.X.EQ.1.) THEN
        GSTRUC=0.
      ELSE
        GSTRUC = (AL/2./PI)*(LOG((1.-X)/X/X/DEL)*(1.+(1.-X)**2)/X-
     ,           2*(1.-X-DEL*X*X)/X)
      END IF
      RETURN
      END
