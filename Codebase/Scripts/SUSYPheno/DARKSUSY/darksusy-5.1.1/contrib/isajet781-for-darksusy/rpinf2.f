CDECK  ID>, RPINF2
*CMZ :-        -24/09/02  14:59:17  by  Peter Richardson
*-- Author :     Peter Richardson  
C----------------------------------------------------------------------- 
       FUNCTION RPINF2(MI1,MI2,MI3,MI4,WIDTH1,WIDTH2,RESM1,RESM2,
     &                 AIN1,AIN2,BIN1,BIN2,TYPE)
C-----------------------------------------------------------------------
C      FUNCTION TO CALCULATE THE INTERFERENCE TERMS FOR 3-BODY
C      RPARITY VIOLATING DECAY RATES
C-----------------------------------------------------------------------
       IMPLICIT NONE
       DOUBLE PRECISION LOW,UPP
       REAL MI1,MI2,MI3,MI4,WIDTH1,WIDTH2,RESM1,RESM2,AIN1,AIN2
     &      ,BIN1,BIN2,RPINF2
       INTEGER I,TYPE
       DOUBLE PRECISION SSDINT,RPINT2,RPINT3
       EXTERNAL SSDINT,RPINT2,RPINT3
C--common block to pass the masses,etc
       COMMON /INFTRM/ M(4),MSQ(4),GAM(2),MR(2),MRSQ(2),A(2),B(2)
       DOUBLE PRECISION M,GAM,MR,MRSQ,A,B,MSQ
C--Set the masses and couplings in the integration routine
       M(1)   = DBLE(MI1)
       M(2)   = DBLE(MI2)
       M(3)   = DBLE(MI3)
       M(4)   = DBLE(MI4)
       A(1)   = DBLE(AIN1)
       A(2)   = DBLE(AIN2)
       B(1)   = DBLE(BIN1)
       B(2)   = DBLE(BIN2)
       GAM(1) = DBLE(WIDTH1)
       GAM(2) = DBLE(WIDTH2)
       MR(1)  = DBLE(RESM1)
       MR(2)  = DBLE(RESM2)
       DO I=1,4
         MSQ(I)=M(I)**2
       ENDDO
       DO I=1,2
         MRSQ(I)=MR(I)**2
       ENDDO
C--Perform the smoothing
       LOW = ATAN(((M(1)+M(2))**2-MR(1)**2)/(GAM(1)*ABS(MR(1))))
     &       /(GAM(1)*ABS(MR(1)))
       UPP = ATAN(((M(4)-M(3))**2-MR(1)**2)/(GAM(1)*ABS(MR(1))))
     &       /(GAM(1)*ABS(MR(1)))
C--Do the outer integral
       IF(TYPE.EQ.1) RPINF2 = REAL(SSDINT(LOW,RPINT2,UPP))
       IF(TYPE.EQ.2) RPINF2 = REAL(SSDINT(LOW,RPINT3,UPP))
       END
