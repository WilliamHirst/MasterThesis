c====================================================================
c
c   auxiliary function used in:  
c   dsi_42.f
c
c   author: piero ullio (piero@tapir.caltech.edu)
c
c____________________________________________________________________

      subroutine dssubkc(r1,delta,res1,res2)
      implicit none
      real*8 r1,delta,res1,res2,int
      if(r1.ge.0.d0) then
        int=log(abs((r1+sqrt(delta))/2.d0))
        res1=int*sqrt(delta)
        res2=-sqrt(delta)/2.d0
      else
        int=-log(abs((r1-sqrt(delta))/2.d0))
        res1=int*sqrt(delta)
        res2=sqrt(delta)/2.d0
      endif
      return
      end
