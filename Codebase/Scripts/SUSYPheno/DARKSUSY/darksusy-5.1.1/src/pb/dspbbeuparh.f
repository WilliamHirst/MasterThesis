**********************************************************************
*** function called in dspbbeuparm
*** it is integrated in the cylindrical coordinate z from 
*** pbhg/pbhh to 1 (linear change of variables such that z=1 => z=pbhh
*** (half height of the diffusion box))
*** version valid in case of constant galactic wind in the z direction
*** 
*** author: piero ullio (piero@tapir.caltech.edu)
*** date: 00-07-13
*** modified: 04-01-22 (pu)
**********************************************************************

      real*8 function dspbbeuparh(z)
      implicit none
      real*8 z,radcoord,vertcoord,dshmaxirho,dshmaxiprob,val,r,kpc,
     & arg
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
        write(*,*) 'dspbbeuparh called with wrong hclumpy=',hclumpy
        stop
      endif
      if(z.lt.pbhg/pbhh.or.z.gt.1.d0) then
        write(*,*) 'dspbbeuparh called with wrong z=',z
        stop
      else
        arg=lambdah*pbhh*kpc*(1.d0-z)
        if(arg.lt.37.d0) then
          dspbbeuparh=dsinh(lambdah*pbhh*kpc*(1.d0-z))
     &             /dsinh(lambdah*(pbhh-pbhg)*kpc)
        else
          dspbbeuparh=dexp(lambdah*kpc*(pbhg-z*pbhh))
        endif
        dspbbeuparh=dspbbeuparh*dexp(gamh*pbhh*kpc*(pbhg/pbhh-z))
     &     *val*pbhh*kpc   
      endif
      return
      end
