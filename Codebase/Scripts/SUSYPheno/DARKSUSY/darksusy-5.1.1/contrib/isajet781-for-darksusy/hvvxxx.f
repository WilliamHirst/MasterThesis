C
C ----------------------------------------------------------------------
C
      SUBROUTINE HVVXXX(V1,V2,G,SMASS,SWIDTH , HVV)
C
C this subroutine computes an off-shell scalar current from the vector- 
C vector-scalar coupling.                                               
C                                                                       
C input:                                                                
C       complex v1(6)          : first  vector                        v1
C       complex v2(6)          : second vector                        v2
C       real    g              : coupling constant                  gvvh
C       real    smass          : mass  of output scalar s               
C       real    swidth         : width of output scalar s               
C                                                                       
C output:                                                               
C       complex hvv(3)         : off-shell scalar current     j(s:v1,v2)
C
      IMPLICIT NONE
      COMPLEX*16 V1(6),V2(6),HVV(3),DG
      REAL*8    Q(0:3),G,SMASS,SWIDTH,Q2
      REAL*8 RXZERO
      PARAMETER( RXZERO=0.0D0 )
C
      HVV(2) = V1(5)+V2(5)
      HVV(3) = V1(6)+V2(6)
C
      Q(0)=DBLE( HVV(2))
      Q(1)=DBLE( HVV(3))
      Q(2)=DIMAG(HVV(3))
      Q(3)=DIMAG(HVV(2))
      Q2=Q(0)**2-(Q(1)**2+Q(2)**2+Q(3)**2)
C
      DG=-G/DCMPLX( Q2-SMASS**2 , MAX(SIGN( SMASS*SWIDTH ,Q2),RXZERO) )
C
      HVV(1) = DG*(V1(1)*V2(1)-V1(2)*V2(2)-V1(3)*V2(3)-V1(4)*V2(4))
C
      RETURN
      END
