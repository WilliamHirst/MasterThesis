      subroutine dsg4setc1234(kp1,kp2,kp3,kp4,vrtx)
c_______________________________________________________________________
c  auxiliary subroutine to dsvertx for quartic couplings 
c  set the value of the 4-particle vertex to vrtx
c  case of four neutral particles
c  author: paolo gondolo (pxg26@po.cwru.edu) 2001
c=======================================================================
      implicit none
      integer kp1,kp2,kp3,kp4
      complex*16 vrtx
      call dsg4setc(kp1,kp2,kp3,kp4,vrtx)
      if (kp1.ne.kp2) call dsg4setc(kp2,kp1,kp3,kp4,vrtx)
      if (kp2.ne.kp3) call dsg4setc(kp1,kp3,kp2,kp4,vrtx)
      return
      end
