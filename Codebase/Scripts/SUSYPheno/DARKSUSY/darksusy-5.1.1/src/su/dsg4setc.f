      subroutine dsg4setc(kp1,kp2,kp3,kp4,vrtx)
c_______________________________________________________________________
c  auxiliary subroutine to dsvertx for quartic couplings 
c  set the value of the 4-particle vertex to vrtx
c  author: paolo gondolo (pxg26@po.cwru.edu) 2001
c=======================================================================
      implicit none
      integer kp1,kp2,kp3,kp4,q,dsqindx
      complex*16 vrtx
      complex*16 dsg4pa(123256)
      common /nvrtxs/ dsg4pa
      save /nvrtxs/
      q=dsqindx(kp1,kp2,kp3,kp4)
      if (q.gt.0) then
         dsg4pa(q)=vrtx
c         write (*,1000) ' ',q,dreal(dsg4pa(q)),dimag(dsg4pa(q)),
c     &        pname(kp1),pname(kp2),pname(kp3),pname(kp4)
      else if (q.lt.0) then
         dsg4pa(-q)=conjg(vrtx)
c         write (*,1000) '-',-q,dreal(dsg4pa(-q)),dimag(dsg4pa(-q)),
c     &        pname(kp1),pname(kp2),pname(kp3),pname(kp4)
      endif
      return
c 1000 format('dsg4set.... ',a1,'q=',i6.6,' (',f9.5,',',f9.5,
c     &     ')   [',a7,',',a7,',',a7,',',a7,']')
      end
