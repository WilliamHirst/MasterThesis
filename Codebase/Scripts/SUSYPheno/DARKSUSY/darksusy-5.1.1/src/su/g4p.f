      function g4p(kp1,kp2,kp3,kp4)
c_______________________________________________________________________
c  function returning the 4-particle vertex
c  author: paolo gondolo (pxg26@po.cwru.edu) 2001
c=======================================================================
      implicit none
      complex*16 dsg4pa(123256)
      common /nvrtxs/ dsg4pa
      save /nvrtxs/
      complex*16 g4p,ans
      integer kp1,kp2,kp3,kp4,q,dsqindx
      q=dsqindx(kp1,kp2,kp3,kp4)
c      write (*,*) 'get.... q=',q
      if (q.gt.0) then
         ans=dsg4pa(q)
      else if (q.lt.0) then
         ans=dconjg(dsg4pa(-q))
      else
         ans=dcmplx(0.d0,0.d0)
      endif
      g4p=ans
      return
      end
