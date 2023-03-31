CDECK  ID>, SSWZ2P.
        REAL FUNCTION SSWZ2P(Q2)
C-----------------------------------------------------------------------
C          SSWIBF: w1ss -> z1ss pi+ pi0
C          Taken from Chen, Drees, Gunion hep-ph/9902309 Eq. A2b
C-----------------------------------------------------------------------
      IMPLICIT NONE
C          Temporary parameters for functions
      COMMON/SSTMP/TMP(10),ITMP(10)
      REAL TMP
      INTEGER ITMP
      SAVE /SSTMP/
C
      REAL Q2
      DOUBLE PRECISION MW1,MZ1,OL,OR,MR,GR,MRP,GRP,QS,
     $RQS,FQSS,MPI,WID,SSDLAM,B
      COMPLEX FQS,ZI
      DATA ZI/(0.,1.)/,MR/.773D0/,GR/.145D0/,MRP/1.37D0/,GRP/.510D0/
      DATA MPI/.140D0/,B/-.145D0/
C
      QS=Q2
      RQS=SQRT(QS)
      MW1=TMP(1)
      MZ1=TMP(2)
      OL=TMP(3)
      OR=TMP(4)
      FQS=(MR**2/(MR**2-QS-ZI*RQS*GR)+B*
     ,(MRP**2/(MRP**2-QS-ZI*RQS*GRP)))/(1.D0+B)
      FQSS=FQS*CONJG(FQS)
      WID=FQSS*(1.D0-4*MPI**2/QS)**1.5D0*SQRT(SSDLAM(MW1**2,MZ1**2,QS))
     ,*((OL**2+OR**2)*(QS*(MW1**2+MZ1**2-2*QS)+(MW1**2-MZ1**2)**2)
     ,-12*OL*OR*QS*MW1*MZ1)
      SSWZ2P=WID
      RETURN
      END
