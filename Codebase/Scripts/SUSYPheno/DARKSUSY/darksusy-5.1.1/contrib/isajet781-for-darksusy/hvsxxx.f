C
C ======================================================================
C ----------------------------------------------------------------------
C
      SUBROUTINE HVSXXX(VC,SC,G,SMASS,SWIDTH , HVS)
C
C this subroutine computes an off-shell scalar current from the vector- 
C scalar-scalar coupling.  the coupling is absent in the minimal sm in  
C unitary gauge.                                                        
C                                                                       
C input:                                                                
C       complex vc(6)          : input vector                          v
C       complex sc(3)          : input scalar                          s
C       complex g              : coupling constant (s charge)           
C       real    smass          : mass  of output scalar s'              
C       real    swidth         : width of output scalar s'              
C                                                                       
C examples of the coupling constant g for susy particles are as follows:
C   -----------------------------------------------------------         
C   |    s1    | (q,i3) of s1  ||   v=a   |   v=z   |   v=w   |         
C   -----------------------------------------------------------         
C   | nu~_l    | (  0  , +1/2) ||   ---   |  gzn(1) |  gwf(1) |         
C   | e~_l     | ( -1  , -1/2) ||  gal(1) |  gzl(1) |  gwf(1) |         
C   | u~_l     | (+2/3 , +1/2) ||  gau(1) |  gzu(1) |  gwf(1) |         
C   | d~_l     | (-1/3 , -1/2) ||  gad(1) |  gzd(1) |  gwf(1) |         
C   -----------------------------------------------------------         
C   | e~_r-bar | ( +1  ,  0  ) || -gal(2) | -gzl(2) | -gwf(2) |         
C   | u~_r-bar | (-2/3 ,  0  ) || -gau(2) | -gzu(2) | -gwf(2) |         
C   | d~_r-bar | (+1/3 ,  0  ) || -gad(2) | -gzd(2) | -gwf(2) |         
C   -----------------------------------------------------------         
C where the sc charge is defined by the flowing-out quantum number.     
C                                                                       
C output:                                                               
C       complex hvs(3)         : scalar current                j(s':v,s)
C
      IMPLICIT NONE
      COMPLEX*16 VC(6),SC(3),HVS(3),DG,QVV,QPV,G
      REAL*8    QV(0:3),QP(0:3),QA(0:3),SMASS,SWIDTH,Q2
C
      HVS(2) = VC(5)+SC(2)
      HVS(3) = VC(6)+SC(3)
C
      QV(0)=DBLE(  VC(5))
      QV(1)=DBLE(  VC(6))
      QV(2)=DIMAG( VC(6))
      QV(3)=DIMAG( VC(5))
      QP(0)=DBLE(  SC(2))
      QP(1)=DBLE(  SC(3))
      QP(2)=DIMAG( SC(3))
      QP(3)=DIMAG( SC(2))
      QA(0)=DBLE( HVS(2))
      QA(1)=DBLE( HVS(3))
      QA(2)=DIMAG(HVS(3))
      QA(3)=DIMAG(HVS(2))
      Q2=QA(0)**2-(QA(1)**2+QA(2)**2+QA(3)**2)
C
      DG=-G/DCMPLX( Q2-SMASS**2 , MAX(DSIGN( SMASS*SWIDTH ,Q2),0D0) )
      QVV=QV(0)*VC(1)-QV(1)*VC(2)-QV(2)*VC(3)-QV(3)*VC(4)
      QPV=QP(0)*VC(1)-QP(1)*VC(2)-QP(2)*VC(3)-QP(3)*VC(4)
C
      HVS(1) = DG*(2D0*QPV+QVV)*SC(1)
C
      RETURN
      END
