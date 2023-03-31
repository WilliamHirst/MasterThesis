c====================================================================
c
c   auxiliary function used in:  
c   dsi_41.f  dsi_42.f
c
c   author: piero ullio (piero@tapir.caltech.edu)
c
c____________________________________________________________________


      real*8 function dssubka(r,delta)
      implicit none
      real*8 r,delta
      if(delta.ge.0.d0) then
        write(*,*) 'dssubka called with wrong argument'
        stop
      else
        dssubka=sqrt(-delta)*atan(r/sqrt(-delta))
      endif
      return
      end
