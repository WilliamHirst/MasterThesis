c====================================================================
c
c   auxiliary function used in:  
c   dslp2.f
c
c   author: piero ullio (piero@tapir.caltech.edu)
c
c____________________________________________________________________

      real*8 function dsilp2(x)
      implicit none
      real*8 x,cc1,cc2,arg_1
      common/dslp2c/cc1,cc2
      if(x.eq.0) then
        dsilp2=cc1
      else
        arg_1=cc1*x+cc2*x*x
        if(dabs(arg_1).le.1.d-40) then
          dsilp2=arg_1/x
        else
          dsilp2=1.d0/x*dlog(1.d0+arg_1)
        endif
      endif
      return
      end
