C ----------------------------------------------------------------------
C
      SUBROUTINE HSSSXX(S1,S2,S3,G,SMASS,SWIDTH , HSSS)
C
C This subroutine computes an off-shell scalar current from the four-   
C scalar coupling.                                                      
C                                                                       
C INPUT:                                                                
C       complex S1(3)          : first  scalar                        S1
C       complex S2(3)          : second scalar                        S2
C       complex S3(3)          : third  scalar                        S3
C       real    G              : coupling constant                 GHHHH
C       real    SMASS          : mass  of OUTPUT scalar S'              
C       real    SWIDTH         : width of OUTPUT scalar S'              
C                                                                       
C OUTPUT:                                                               
C       complex HSSS(3)        : scalar current           J(S':S1,S2,S3)
C
      IMPLICIT NONE
      COMPLEX*16 S1(3),S2(3),S3(3),HSSS(3),DG
      REAL*8     Q(0:3),G,SMASS,SWIDTH,Q2
C
      HSSS(2) = S1(2)+S2(2)+S3(2)
      HSSS(3) = S1(3)+S2(3)+S3(3)
C
      Q(0)=DBLE( HSSS(2))
      Q(1)=DBLE( HSSS(3))
      Q(2)=DIMAG(HSSS(3))
      Q(3)=DIMAG(HSSS(2))
      Q2=Q(0)**2-(Q(1)**2+Q(2)**2+Q(3)**2)
C
      DG=-G/DCMPLX( Q2-SMASS**2,MAX(SIGN(SMASS*SWIDTH ,Q2),0.D0))
C
      HSSS(1) = DG * S1(1)*S2(1)*S3(1)
C
      RETURN
      END
