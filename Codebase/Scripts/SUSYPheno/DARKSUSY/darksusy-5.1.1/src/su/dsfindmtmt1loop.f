      subroutine dsfindmtmt1loop
      implicit none
      include 'dsmssm.h'
      real*8 ck,dsmqpole1loop
c
      mtmt=mass(kt)
 100  ck=mass(kt)-dsmqpole1loop(mtmt)
      if(ck.lt.0.d0) then
        mtmt=mtmt-0.1d0
        goto 100
      endif
c      write(*,*) 'for alph3 and mass(kt) = ',alph3,mass(kt)
c      write(*,*) 'mtmt = ',mtmt
      return
      end  
