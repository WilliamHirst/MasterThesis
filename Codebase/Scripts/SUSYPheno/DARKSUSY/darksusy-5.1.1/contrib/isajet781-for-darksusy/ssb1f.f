CDECK  ID>, SSB1F.
       COMPLEX*16 FUNCTION SSB1F(XS,XMI,XMJ)
C          Implemented by Javier 9/8/05 to remove the Log
C          thresholds already implemented in the RGEs
C          through step by step decoupling
      COMMON/SSINF/XLAM
      DOUBLE PRECISION XLAM
       REAL XS,XMI,XMJ
       DOUBLE PRECISION S,MI,MJ
       COMPLEX*16 SSB1
       S=XS
       MI=XMI
       MJ=XMJ
       SSB1F=SSB1(XS,XMI,XMJ)+LOG(MAX(MI,MJ))-XLAM/2.D0
       RETURN
       END
