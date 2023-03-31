**************************************************************
*** SUBROUTINE dsasffcol                                   ***
*** computes the amplitude squared of the process          ***
*** scalar(1) + scalar(2) -> fermion(3) + fermion(4)       ***
***                                                        *** 
*** input:                                                 ***
*** asparmass, askin, askinder variables                   ***
*** complex vectors:                                       ***
***   ASxplc(j,i), ASxprc(j,i), ASylc(j), ASyrc(j)         *** 
*** AUTHOR: Piero Ullio, ullio@sissa.it                    ***
*** Date: 01-02-28                                         ***
**************************************************************


      SUBROUTINE dsasffcol(ampl2)
      implicit none
      include 'dsascom.h'
      real*8 ampl2
      real*8 term1,term2,term3,term4,term5
      complex*16 xlk2(6),xlk1(6),xrk2(6),xrk1(6),xlxlsk1k2(6,6)
     & ,xrxrsk1k2(6,6),xlxrs(6,6)
      integer i,j,j1,j2,k
***** some preliminar quantities
      do j=1,6
        xlk2(j)=0.d0
        xlk1(j)=0.d0
        xrk2(j)=0.d0 
        xrk1(j)=0.d0
      enddo  
      do j1=1,6
      do j2=1,6   
        xlxlsk1k2(j1,j2)=0.d0
        xrxrsk1k2(j1,j2)=0.d0
        xlxrs(j1,j2)=0.d0
      enddo  
      enddo
      do j=1,6
        do i=1,4
          xlk2(j)=xlk2(j)+ASxplc(j,i)*scal(4,i) 
          xlk1(j)=xlk1(j)+ASxplc(j,i)*scal(3,i) 
          xrk2(j)=xrk2(j)+ASxprc(j,i)*scal(4,i) 
          xrk1(j)=xrk1(j)+ASxprc(j,i)*scal(3,i)
        enddo
      enddo  
      do j1=1,6
      do j2=1,6   
        do i=1,4
        do k=1,4
          xlxlsk1k2(j1,j2)=xlxlsk1k2(j1,j2)
     &      +ASxplc(j1,i)*dconjg(ASxplc(j2,k))*scal(i,k)*scal(3,4)
          xrxrsk1k2(j1,j2)=xrxrsk1k2(j1,j2)
     &      +ASxprc(j1,i)*dconjg(ASxprc(j2,k))*scal(i,k)*scal(3,4)
          xlxrs(j1,j2)=xlxrs(j1,j2)
     &      +ASxplc(j1,i)*dconjg(ASxprc(j2,k))*scal(i,k)
        enddo
        enddo
      enddo  
      enddo
***** 5 terms in the amplitude squared expression
      term1=0.d0
      do j1=1,6
      do j2=1,6   
        term1=term1+colfactor(j1,j2)*2.d0*
     &    (2.d0*dreal(xlk2(j1)*dconjg(xlk1(j2)))
     &     -xlxlsk1k2(j1,j2)+2.d0*dreal(xrk2(j1)*dconjg(xrk1(j2)))
     &     -xrxrsk1k2(j1,j2))
      enddo
      enddo
*****
      term2=0.d0
      do j1=1,6
      do j2=1,6   
        term2=term2+colfactor(j1,j2)*
     &    dreal(ASylc(j1)*dconjg(ASyrc(j2))+xlxrs(j1,j2))
      enddo
      enddo
      term2=(-1.d0)**(q3+q4)*mass3*mass4*4.d0*term2
*****
      term3=0.d0
      do j1=1,6
      do j2=1,6   
        term3=term3+colfactor(j1,j2)*
     &    dreal(xlk1(j1)*dconjg(ASylc(j2))+xrk1(j1)*dconjg(ASyrc(j2)))
      enddo
      enddo
      term3=(-1.d0)**(q4)*mass4*4.d0*term3
*****
      term4=0.d0
      do j1=1,6
      do j2=1,6   
        term4=term4+colfactor(j1,j2)*
     &    dreal(xlk2(j1)*dconjg(ASyrc(j2))+xrk2(j1)*dconjg(ASylc(j2)))
      enddo
      enddo
      term4=(-1.d0)**(q3)*mass3*4.d0*term4
*****
      term5=0.d0
      do j1=1,6
      do j2=1,6   
        term5=term5+colfactor(j1,j2)*
     &        (dreal(ASylc(j1)*dconjg(ASylc(j2)))
     &         +dreal(ASyrc(j1)*dconjg(ASyrc(j2))))
      enddo
      enddo
      term5=2.d0*scal(3,4)*term5
*****      
***** sum the terms:
      ampl2=term1+term2+term3+term4+term5
*****
      return
      end
