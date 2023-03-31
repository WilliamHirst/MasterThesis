**********************************************************************
*** function called in dspbbeuparm
*** it is integrated in the cylindrical coordinate z from 0 to
*** pbhg/pbhh (linear change of variables such that z=1 => z=pbhh
*** (half height of the diffusion box)) - part associated with sinh
*** version valid in case of constant galactic wind in the z direction
*** 
*** author: piero ullio (piero@tapir.caltech.edu)
*** date: 00-07-13
*** modified: 04-01-22 (pu)
**********************************************************************

      real*8 function dspbbeupargs(z)
      implicit none
      real*8 z,radcoord,vertcoord,dshmaxirho,dshmaxiprob,val,r,kpc
     & ,arg
      include 'dspbcom.h'
      include 'dspbprivate.h'
      include 'dshmcom.h'
      common /czint/r
      parameter(kpc=3.08567802d0)
      if(hclumpy.eq.1) then                     ! smooth halo
        radcoord=r*pbrh
        vertcoord=z*pbhh
        val=dshmaxirho(radcoord,vertcoord)
        val=val/rho0
        val=val**2
      elseif(hclumpy.eq.2) then                 ! clumpy halo
        radcoord=r*pbrh
        vertcoord=z*pbhh
        val=dshmaxiprob(radcoord,vertcoord)
        val=val/rho0
      else
        write(*,*) 'dspbbeupargs called with wrong hclumpy=',hclumpy
        stop
      endif
      if(z.lt.0.d0.or.z.gt.pbhg/pbhh) then
        write(*,*) 'dspbbeupargs called with wrong z=',z
        stop
      else
        arg=lambdag*pbhh*kpc*(pbhg/pbhh-z)
        if(arg.lt.20.d0) then
          dspbbeupargs=dsinh(lambdag*pbhh*kpc*(pbhg/pbhh-z))
     &              /dsinh(lambdag*pbhg*kpc)
        else
          dspbbeupargs=dcosh(lambdag*pbhh*kpc*z)
     &        -dsinh(lambdag*pbhh*kpc*z)/dtanh(lambdag*pbhg*kpc)
        endif  
        dspbbeupargs=dspbbeupargs*val*pbhh*kpc
      endif
      return
      end
