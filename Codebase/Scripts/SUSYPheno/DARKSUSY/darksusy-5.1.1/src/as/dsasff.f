**************************************************************
*** SUBROUTINE dsasff                                      ***
*** computes the amplitude squared of the process          ***
*** scalar(1) + scalar(2) -> fermion(3) + fermion(4)       ***
***                                                        *** 
*** input:                                                 ***
*** asparmass, askin, askinder variables                   ***
*** complex vectors ASxpl(i), ASxpr(i), ASyl, ASyr         *** 
*** AUTHOR: Piero Ullio, ullio@sissa.it                    ***
*** Date: 01-02-28                                         ***
**************************************************************

      SUBROUTINE dsasff(ampl2)
      implicit none
      include 'dsascom.h'
      real*8 ampl2,dsabsq
      real*8 term1,term2,term3,term4,term5
      complex*16 xlk2,xlk1,xrk2,xrk1,xlxlsk1k2,xrxrsk1k2,
     & xlxrs
      integer i,j
***** some preliminar quantities
      xlk2=0.d0
      xlk1=0.d0
      xrk2=0.d0 
      xrk1=0.d0
      xlxlsk1k2=0.d0
      xrxrsk1k2=0.d0
      xlxrs=0.d0
      do i=1,4
        xlk2=xlk2+ASxpl(i)*scal(4,i) 
        xlk1=xlk1+ASxpl(i)*scal(3,i) 
        xrk2=xrk2+ASxpr(i)*scal(4,i) 
        xrk1=xrk1+ASxpr(i)*scal(3,i) 
        do j=1,4
          xlxlsk1k2=xlxlsk1k2
     &      +ASxpl(i)*dconjg(ASxpl(j))*scal(i,j)*scal(3,4)
          xrxrsk1k2=xrxrsk1k2
     &      +ASxpr(i)*dconjg(ASxpr(j))*scal(i,j)*scal(3,4)
          xlxrs=xlxrs
     &      +ASxpl(i)*dconjg(ASxpr(j))*scal(i,j)
        enddo
      enddo  
***** 5 terms in the amplitude squared expression
      term1=2.d0*(2.d0*dreal(xlk2*dconjg(xlk1))
     & -xlxlsk1k2+2.d0*dreal(xrk2*dconjg(xrk1))
     & -xrxrsk1k2)
*****
      term2=(-1.d0)**(q3+q4)*mass3*mass4*4.d0*
     & dreal(ASyl*dconjg(ASyr)+xlxrs)
*****
      term3=(-1.d0)**(q4)*mass4*4.d0
     & *dreal(xlk1*dconjg(ASyl)+xrk1*dconjg(ASyr))
*****
      term4=(-1.d0)**(q3)*mass3*4.d0
     & *dreal(xlk2*dconjg(ASyr)+xrk2*dconjg(ASyl))
*****
      term5=2.d0*scal(3,4)*(dsabsq(ASyl)+dsabsq(ASyr))
*****      
***** sum the terms:
      ampl2=term1+term2+term3+term4+term5
*****
      return
      end
