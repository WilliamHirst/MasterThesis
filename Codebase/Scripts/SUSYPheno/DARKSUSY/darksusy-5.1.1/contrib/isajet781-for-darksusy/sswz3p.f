CDECK  ID>, SSWZ3P.
        REAL FUNCTION SSWZ3P(Q2)
C-----------------------------------------------------------------------
C          SSWIBF: w1ss -> z1ss pi+ pi0 pi0 or z1ss pi+ pi- pi+
C          Taken from Chen, Drees, Gunion hep-ph/9902309 Eq. A2c
C-----------------------------------------------------------------------
      IMPLICIT NONE
C          Temporary parameters for functions
      COMMON/SSTMP/TMP(10),ITMP(10)
      REAL TMP
      INTEGER ITMP
      SAVE /SSTMP/
C
      REAL Q2
      DOUBLE PRECISION MW1,MZ1,OL,OR,MA2,GA2,QS,
     $RQS,FQSS,MPI,WID,SSDLAM,MR,GQS
      COMPLEX FQS,ZI
      DATA ZI/(0.,1.)/,MA2/1.318D0/,GA2/.107D0/
      DATA MPI/.140D0/,MR/.773D0/
C
      QS=Q2
      RQS=SQRT(QS)
      MW1=TMP(1)
      MZ1=TMP(2)
      OL=TMP(3)
      OR=TMP(4)
      FQS=MA2**2/(MA2**2-QS-ZI*RQS*GA2)
      FQSS=FQS*CONJG(FQS)
      IF (QS.LT.(MR+MPI)) THEN
        GQS=4.1*(QS-9*MPI**2)**3*(1.D0-3.3*(QS-9*MPI**2)+
     ,      5.8*(QS-9*MPI**2)**2)
      ELSE
        GQS=QS*(1.623D0+10.38/QS-9.32/QS/QS+.65/QS**3)
      END IF
      WID=FQSS*SQRT(SSDLAM(MW1**2,MZ1**2,QS))*GQS
     ,*((OL**2+OR**2)*(MW1**2+MZ1**2-2*QS+(MW1**2-MZ1**2)**2/QS)
     ,-12*OL*OR*MW1*MZ1)
      SSWZ3P=1.35*WID
      RETURN
      END
