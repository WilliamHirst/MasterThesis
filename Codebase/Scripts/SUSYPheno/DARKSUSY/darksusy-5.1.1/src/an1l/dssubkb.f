c====================================================================
c
c   auxiliary function used in:  
c   dsi_41.f  dsi_42.f
c
c   author: piero ullio (piero@tapir.caltech.edu)
c
c____________________________________________________________________

      subroutine dssubkb(r1,r2,delta,res1,res2)
      implicit none
      real*8 r1,r2,delta,res1,res2,int,tst1,tst2
      if(r1.ge.0.d0.and.r2.ge.0.d0) then
        int=log(abs((r1+sqrt(delta))/(r2+sqrt(delta))))
        res1=int*sqrt(delta)
        res2=-sqrt(delta)/2.d0
      else if(r1.lt.0.d0.and.r2.lt.0.d0) then
        int=-log(abs((r1-sqrt(delta))/(r2-sqrt(delta))))
        res1=int*sqrt(delta)
        res2=sqrt(delta)/2.d0
      else if(r1.ge.0.d0.and.r2.lt.0.d0) then
        tst1=abs(abs(r1)-sqrt(delta))
        tst2=abs(abs(r2)-sqrt(delta))
        if(tst1.le.tst2) then
          int=log(abs((r1+sqrt(delta))/(r2+sqrt(delta))))
          res1=int*sqrt(delta)
          res2=-sqrt(delta)/2.d0
        else
          int=-log(abs((r1-sqrt(delta))/(r2-sqrt(delta))))
          res1=int*sqrt(delta)
          res2=sqrt(delta)/2.d0
        endif
      else
        tst1=abs(abs(r1)-sqrt(delta))
        tst2=abs(abs(r2)-sqrt(delta))
        if(tst1.le.tst2) then
          int=-log(abs((r1-sqrt(delta))/(r2-sqrt(delta))))
          res1=int*sqrt(delta)
          res2=sqrt(delta)/2.d0
        else
          int=log(abs((r1+sqrt(delta))/(r2+sqrt(delta))))
          res1=int*sqrt(delta)
          res2=-sqrt(delta)/2.d0
        endif
      endif
      return
      end
