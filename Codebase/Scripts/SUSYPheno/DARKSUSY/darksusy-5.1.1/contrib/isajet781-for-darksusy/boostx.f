CDECK  ID>, DHELAS.
C  *********************************************************************
C  ***                                                               ***
C  ***               coded by H. Murayama & I. Watanabe              ***
C  *** For the formalism and notations, see the following reference: ***
C  ***           H. Murayama, I. Watanabe and K. Hagiwara            ***
C  ***           "HELAS: HELicity Amplitude Subroutines              ***
C  ***               for Feynman diagram evaluation"                 ***
C  ***               KEK Report 91-11, December 1991                 ***
C  ***                                                               ***
C  *********************************************************************
C
C  Converted to double precision by W. Long and T. Seltzer for MadGraph.
C
C  Minor changes for portability by FEP, July 1999. The code is not ANSI
C  standard, but that cannot be helped if MadGraph compatibility is to 
C  be maintained.
C
C ======================================================================
C
      SUBROUTINE BOOSTX(P,Q , PBOOST)
C
C this subroutine performs the lorentz boost of a four-momentum.  the   
C momentum p is assumed to be given in the rest frame of q.  pboost is  
C the momentum p boosted to the frame in which q is given.  q must be a 
C timelike momentum.                                                    
C                                                                       
C input:                                                                
C       real    p(0:3)         : four-momentum p in the q rest  frame   
C       real    q(0:3)         : four-momentum q in the boosted frame   
C                                                                       
C output:                                                               
C       real    pboost(0:3)    : four-momentum p in the boosted frame   
C
      IMPLICIT NONE
      REAL*8    P(0:3),Q(0:3),PBOOST(0:3),PQ,QQ,M,LF
      REAL*8 RXZERO
      PARAMETER( RXZERO=0.0D0 )
C
      QQ=Q(1)**2+Q(2)**2+Q(3)**2
C
      IF ( QQ .NE. RXZERO ) THEN
         PQ=P(1)*Q(1)+P(2)*Q(2)+P(3)*Q(3)
         M=SQRT(Q(0)**2-QQ)
         LF=((Q(0)-M)*PQ/QQ+P(0))/M
         PBOOST(0) = (P(0)*Q(0)+PQ)/M
         PBOOST(1) =  P(1)+Q(1)*LF
         PBOOST(2) =  P(2)+Q(2)*LF
         PBOOST(3) =  P(3)+Q(3)*LF
      ELSE
         PBOOST(0)=P(0)
         PBOOST(1)=P(1)
         PBOOST(2)=P(2)
         PBOOST(3)=P(3)
      ENDIF
C
      RETURN
      END
