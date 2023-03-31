      subroutine dsg4setc13(kp1,kp2,kp3,kp4,vrtx)
c_______________________________________________________________________
c  auxiliary subroutine to dsvertx for quartic couplings 
c  set the value of the 4-particle vertex to vrtx
c  case of two neutral particles (1 and 3)
c  author: paolo gondolo (pxg26@po.cwru.edu) 2001
c=======================================================================
      implicit none
      integer kp1,kp2,kp3,kp4
      complex*16 vrtx
      call dsg4setc(kp1,kp2,kp3,kp4,vrtx)
      if (kp1.ne.kp2) call dsg4set(kp2,kp1,kp3,kp4,0.d0,0.d0)
      if (kp2.ne.kp3) call dsg4setc(kp1,kp3,kp2,kp4,conjg(vrtx))
      return
      end
