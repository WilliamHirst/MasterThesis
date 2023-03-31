C
C ----------------------------------------------------------------------
C
      SUBROUTINE ROTXXX(P,Q , PROT)
C
C this subroutine performs the spacial rotation of a four-momentum.     
C the momentum p is assumed to be given in the frame where the spacial  
C component of q points the positive z-axis.  prot is the momentum p    
C rotated to the frame where q is given.                                
C                                                                       
C input:                                                                
C       real    p(0:3)         : four-momentum p in q(1)=q(2)=0 frame   
C       real    q(0:3)         : four-momentum q in the rotated frame   
C                                                                       
C output:                                                               
C       real    prot(0:3)      : four-momentum p in the rotated frame   
C
      IMPLICIT NONE
      REAL*8    P(0:3),Q(0:3),PROT(0:3),QT2,QT,PSGN,QQ,P1
C
      REAL*8 RXZERO, RXONE
      PARAMETER( RXZERO=0.0D0, RXONE=1.0D0 )
C
      PROT(0) = P(0)
C
      QT2=Q(1)**2+Q(2)**2
C
      IF ( QT2 .EQ. RXZERO ) THEN
          IF ( Q(3) .EQ. RXZERO ) THEN
             PROT(1) = P(1)
             PROT(2) = P(2)
             PROT(3) = P(3)
          ELSE
             PSGN=DSIGN(RXONE,Q(3))
             PROT(1) = P(1)*PSGN
             PROT(2) = P(2)*PSGN
             PROT(3) = P(3)*PSGN
          ENDIF
      ELSE
          QQ=SQRT(QT2+Q(3)**2)
          QT=SQRT(QT2)
          P1=P(1)
          PROT(1) = Q(1)*Q(3)/QQ/QT*P1 -Q(2)/QT*P(2) +Q(1)/QQ*P(3)
          PROT(2) = Q(2)*Q(3)/QQ/QT*P1 +Q(1)/QT*P(2) +Q(2)/QQ*P(3)
          PROT(3) =          -QT/QQ*P1               +Q(3)/QQ*P(3)
      ENDIF
C
      RETURN
      END
