C
C ======================================================================
C
      SUBROUTINE VSSXXX(VC,S1,S2,G , VERTEX)
C
C this subroutine computes an amplitude from the vector-scalar-scalar   
C coupling.  the coupling is absent in the minimal sm in unitary gauge. 
C                                                                       
C       complex vc(6)          : input  vector                        v 
C       complex s1(3)          : first  scalar                        s1
C       complex s2(3)          : second scalar                        s2
C       complex g              : coupling constant (s1 charge)          
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
C where the s1 charge is defined by the flowing-out quantum number.     
C                                                                       
C output:                                                               
C       complex vertex         : amplitude                gamma(v,s1,s2)
C
      IMPLICIT NONE
      COMPLEX*16 VC(6),S1(3),S2(3),VERTEX,G
      REAL*8    P(0:3)
C
      P(0)=DBLE( S1(2)-S2(2))
      P(1)=DBLE( S1(3)-S2(3))
      P(2)=DIMAG(S1(3)-S2(3))
      P(3)=DIMAG(S1(2)-S2(2))
C
      VERTEX = G*S1(1)*S2(1)
     &        *(VC(1)*P(0)-VC(2)*P(1)-VC(3)*P(2)-VC(4)*P(3))
C
      RETURN
      END
