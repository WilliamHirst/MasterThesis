**********************************************************************
*** function called in dspbcharpar1
*** 
*** author: piero ullio (piero@tapir.caltech.edu)
*** date: 00-07-13
**********************************************************************

      real*8 function dspbcharpar2(y)
      implicit none
      real*8 radcoord,vertcoord,dshmaxirho,dshmaxiprob,val,xx,y
      include 'dspbcom.h'
      include 'dspbprivate.h'
      include 'dshmcom.h'
      common /pbchar/xx
      if(hclumpy.eq.1) then                     ! smooth halo
        radcoord=xx*pbrh
        vertcoord=y*pbhh
        val=dshmaxirho(radcoord,vertcoord)
        val=val/rho0
        val=val**2
      elseif(hclumpy.eq.2)then                  ! clumpy halo
        radcoord=xx*pbrh
        vertcoord=y*pbhh
        val=dshmaxiprob(radcoord,vertcoord)
        val=val/rho0
      else
        write(*,*) 'dspbcharpar2 called with wrong hclumpy=',hclumpy
      endif
      dspbcharpar2=val*dsinh(zero*pbhh/pbrh*(1.d0-y))
     &          /dsinh(zero*pbhh/pbrh)
      return
      end
