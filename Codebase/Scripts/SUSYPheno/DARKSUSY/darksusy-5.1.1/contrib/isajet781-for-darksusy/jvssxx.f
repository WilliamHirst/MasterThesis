C ----------------------------------------------------------------------
C
      SUBROUTINE JVSSXX(VC,S1,S2,G,VMASS,VWIDTH , JVSS)
C
C This subroutine computes an off-shell vector current from the vector- 
C vector-scalar-scalar coupling.  The vector propagator is given in     
C Feynman gauge for a massless vector and in unitary gauge for a massive
C vector.                                                               
C                                                                       
C INPUT:                                                                
C       complex VC(6)          : input  vector                        V 
C       complex S1(3)          : first  scalar                        S1
C       complex S2(3)          : second scalar                        S2
C       real    G              : coupling constant                 GVVHH
C       real    VMASS          : mass  of OUTPUT vector V'              
C       real    VWIDTH         : width of OUTPUT vector V'              
C                                                                       
C OUTPUT:                                                               
C       complex JVSS(6)        : vector current         J^mu(V':V,S1,S2)
C
      IMPLICIT NONE
      COMPLEX*16 VC(6),S1(3),S2(3),JVSS(6),DG
      REAL*8    Q(0:3),G,VMASS,VWIDTH,Q2,VK,VM2
C
      JVSS(5) = VC(5)+S1(2)+S2(2)
      JVSS(6) = VC(6)+S1(3)+S2(3)
C
      Q(0)=DBLE( JVSS(5))
      Q(1)=DBLE( JVSS(6))
      Q(2)=DIMAG(JVSS(6))
      Q(3)=DIMAG(JVSS(5))
      Q2=Q(0)**2-(Q(1)**2+Q(2)**2+Q(3)**2)
      VM2=VMASS**2
C
      IF (VMASS.EQ.0.) GOTO 10
C
      DG=G*S1(1)*S2(1)/DCMPLX( Q2-VM2,MAX(SIGN( VMASS*VWIDTH,Q2),0.D0))
C  For the running width, use below instead of the above DG.
C      DG=G*S1(1)*S2(1)/CMPLX( Q2-VM2 , MAX( VWIDTH*Q2/VMASS ,0.))
C
      VK=(Q(0)*VC(1)-Q(1)*VC(2)-Q(2)*VC(3)-Q(3)*VC(4))/VM2
C
      JVSS(1) = DG*(VC(1)-VK*Q(0))
      JVSS(2) = DG*(VC(2)-VK*Q(1))
      JVSS(3) = DG*(VC(3)-VK*Q(2))
      JVSS(4) = DG*(VC(4)-VK*Q(3))
C
      RETURN
C
  10  DG= G*S1(1)*S2(1)/Q2
C
      JVSS(1) = DG*VC(1)
      JVSS(2) = DG*VC(2)
      JVSS(3) = DG*VC(3)
      JVSS(4) = DG*VC(4)
C
      RETURN
      END
