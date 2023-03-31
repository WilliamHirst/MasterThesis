CDECK  ID>, RPINF1
*CMZ :-        -24/09/02  14:59:17  by  Peter Richardson
*-- Author :     Peter Richardson  
C----------------------------------------------------------------------- 
       FUNCTION RPINF1(MI1,MI2,MI3,MI4,WIDTH,RESM,AIN,BIN)
C-----------------------------------------------------------------------
C      FUNCTION TO RETURN THE INTEGRATED AMPLITUDE SQUARE PIECE OF A
C      THREE BODY DECAY MATRIX ELEMENT
C-----------------------------------------------------------------------
       DOUBLE PRECISION LOW,UPP
       REAL MI1,MI2,MI3,MI4,WIDTH,RESM,AIN,BIN,RPINF1
       DOUBLE PRECISION SSDINT,RPINT1
       EXTERNAL SSDINT,RPINT1
C--common block to pass the masses,etc
       COMMON /SQUTRM/ M(4),GAM,MR,A,B
       DOUBLE PRECISION M,GAM,MR,A,B
C--Set the masses and couplings in the integration routine
       M(1) = DBLE(MI1)
       M(2) = DBLE(MI2)
       M(3) = DBLE(MI3)
       M(4) = DBLE(MI4)
       A    = DBLE(AIN)
       B    = DBLE(BIN)
       GAM  = DBLE(WIDTH)
       MR   = DBLE(ABS(RESM))
C--Perform the smoothing
       LOW = ATAN(((M(1)+M(2))**2-MR**2)/(GAM*MR))/(GAM*MR)
       UPP = ATAN(((M(4)-M(3))**2-MR**2)/(GAM*MR))/(GAM*MR)
C--Do the outer integral
       RPINF1 = 0.5E0*REAL(SSDINT(LOW,RPINT1,UPP))
       END
