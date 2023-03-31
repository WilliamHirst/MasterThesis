C
C
C ----------------------------------------------------------------------
C
      SUBROUTINE JVSXXX(VC,SC,G,VMASS,VWIDTH , JVS)
C
C this subroutine computes an off-shell vector current from the vector- 
C vector-scalar coupling.  the vector propagator is given in feynman    
C gauge for a massless vector and in unitary gauge for a massive vector.
C                                                                       
C input:                                                                
C       complex vc(6)          : input vector                          v
C       complex sc(3)          : input scalar                          s
C       real    g              : coupling constant                  gvvh
C       real    vmass          : mass  of output vector v'              
C       real    vwidth         : width of output vector v'              
C                                                                       
C output:                                                               
C       complex jvs(6)         : vector current             j^mu(v':v,s)
C
      IMPLICIT NONE
      COMPLEX*16 VC(6),SC(3),JVS(6),DG,VK
      REAL*8    Q(0:3),VMASS,VWIDTH,Q2,VM2,G
C
      JVS(5) = VC(5)+SC(2)
      JVS(6) = VC(6)+SC(3)
C
      Q(0)=DBLE( JVS(5))
      Q(1)=DBLE( JVS(6))
      Q(2)=DIMAG(JVS(6))
      Q(3)=DIMAG(JVS(5))
      Q2=Q(0)**2-(Q(1)**2+Q(2)**2+Q(3)**2)
      VM2=VMASS**2
C
      IF (VMASS.EQ.0.) GOTO 10
C
      DG=G*SC(1)/DCMPLX( Q2-VM2 , MAX(DSIGN( VMASS*VWIDTH ,Q2),0.D0) )
C  for the running width, use below instead of the above dg.
C      dg=g*sc(1)/dcmplx( q2-vm2 , max( vwidth*q2/vmass ,0.) )
C
      VK=(-Q(0)*VC(1)+Q(1)*VC(2)+Q(2)*VC(3)+Q(3)*VC(4))/VM2
C
      JVS(1) = DG*(Q(0)*VK+VC(1))
      JVS(2) = DG*(Q(1)*VK+VC(2))
      JVS(3) = DG*(Q(2)*VK+VC(3))
      JVS(4) = DG*(Q(3)*VK+VC(4))
C
      RETURN
C
  10  DG=G*SC(1)/Q2
C
      JVS(1) = DG*VC(1)
      JVS(2) = DG*VC(2)
      JVS(3) = DG*VC(3)
      JVS(4) = DG*VC(4)
C
      RETURN
      END
