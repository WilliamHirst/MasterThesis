C ----------------------------------------------------------------------
C
      SUBROUTINE JSSXXX(S1,S2,G,VMASS,VWIDTH , JSS)
C
C This subroutine computes an off-shell vector current from the vector- 
C scalar-scalar coupling.  The coupling is absent in the minimal SM in  
C unitary gauge.  The propagator is given in Feynman gauge for a        
C massless vector and in unitary gauge for a massive vector.            
C                                                                       
C INPUT:                                                                
C       complex S1(3)          : first  scalar                        S1
C       complex S2(3)          : second scalar                        S2
C       real    G              : coupling constant (S1 charge)          
C       real    VMASS          : mass  of OUTPUT vector V               
C       real    VWIDTH         : width of OUTPUT vector V               
C                                                                       
C Examples of the coupling constant G for SUSY particles are as follows:
C   -----------------------------------------------------------         
C   |    S1    | (Q,I3) of S1  ||   V=A   |   V=Z   |   V=W   |         
C   -----------------------------------------------------------         
C   | nu~_L    | (  0  , +1/2) ||   ---   |  GZN(1) |  GWF(1) |         
C   | e~_L     | ( -1  , -1/2) ||  GAL(1) |  GZL(1) |  GWF(1) |         
C   | u~_L     | (+2/3 , +1/2) ||  GAU(1) |  GZU(1) |  GWF(1) |         
C   | d~_L     | (-1/3 , -1/2) ||  GAD(1) |  GZD(1) |  GWF(1) |         
C   -----------------------------------------------------------         
C   | e~_R-bar | ( +1  ,  0  ) || -GAL(2) | -GZL(2) | -GWF(2) |         
C   | u~_R-bar | (-2/3 ,  0  ) || -GAU(2) | -GZU(2) | -GWF(2) |         
C   | d~_R-bar | (+1/3 ,  0  ) || -GAD(2) | -GZD(2) | -GWF(2) |         
C   -----------------------------------------------------------         
C where the S1 charge is defined by the flowing-OUT quantum number.     
C                                                                       
C OUTPUT:                                                               
C       complex JSS(6)         : vector current            J^mu(V:S1,S2)
C
      IMPLICIT NONE
      COMPLEX*16 S1(3),S2(3),JSS(6),DG,ADG
      REAL*8  PP(0:3),PA(0:3),Q(0:3),G,VMASS,VWIDTH,Q2,VM2,MP2,MA2,M2D
C
      JSS(5) = S1(2)+S2(2)
      JSS(6) = S1(3)+S2(3)
C
      Q(0)=DBLE( JSS(5))
      Q(1)=DBLE( JSS(6))
      Q(2)=DIMAG(JSS(6))
      Q(3)=DIMAG(JSS(5))
      Q2=Q(0)**2-(Q(1)**2+Q(2)**2+Q(3)**2)
      VM2=VMASS**2
C
      IF (VMASS.EQ.0.) GOTO 10
C
      DG=G/DCMPLX( Q2-VM2, MAX(SIGN( VMASS*VWIDTH ,Q2),0.D0))
C  For the running width, use below instead of the above DG.
C      DG=G/dCMPLX( Q2-VM2 , MAX( VWIDTH*Q2/VMASS ,0.) )
C
      ADG=DG*S1(1)*S2(1)
C
      PP(0)=DBLE( S1(2))
      PP(1)=DBLE( S1(3))
      PP(2)=DIMAG(S1(3))
      PP(3)=DIMAG(S1(2))
      PA(0)=DBLE( S2(2))
      PA(1)=DBLE( S2(3))
      PA(2)=DIMAG(S2(3))
      PA(3)=DIMAG(S2(2))
      MP2=PP(0)**2-(PP(1)**2+PP(2)**2+PP(3)**2)
      MA2=PA(0)**2-(PA(1)**2+PA(2)**2+PA(3)**2)
      M2D=MP2-MA2
C
      JSS(1) = ADG*( (PP(0)-PA(0)) - Q(0)*M2D/VM2)
      JSS(2) = ADG*( (PP(1)-PA(1)) - Q(1)*M2D/VM2)
      JSS(3) = ADG*( (PP(2)-PA(2)) - Q(2)*M2D/VM2)
      JSS(4) = ADG*( (PP(3)-PA(3)) - Q(3)*M2D/VM2)
C
      RETURN
C
  10  ADG=G*S1(1)*S2(1)/Q2
C
      JSS(1) = ADG*DBLE( S1(2)-S2(2))
      JSS(2) = ADG*DBLE( S1(3)-S2(3))
      JSS(3) = ADG*DIMAG(S1(3)-S2(3))
      JSS(4) = ADG*DIMAG(S1(2)-S2(2))
C
      RETURN
      END
