**********************************************************************
*** auxiliary function needed in dspbtd15beucl
***
*** diffusion constant in units of 10^27 cm^2 s^-1
*** axec in mb*10^10 cm s^-1
*** lambdag, lambdah in 10^-21 cm^-1
*** addterm in 1/(10^27 cm^2 s^-1 * 10^-21 cm^-1) 
***           = 1/10^6 s/cm
**********************************************************************

      real*8 function dspbaddterm(k,nusk,Jklocal,Jkplus1squared)
      implicit none
      include 'dspbcom.h'
      include 'dspbprivate.h'
      integer k
      real*8 nusk,Jklocal,Jkplus1squared,coeff,dsbessjw,parcoeff
     & ,arg,factors,factorc,rebiggs, rebiggc
      real*8 kpc
      parameter (kpc=3.08567802d0)
      lambdag=dsqrt(axsec*pbng/dg/100.d0+(nusk/pbrh/kpc)**2)
      lambdah=dsqrt(gamh**2+axsec*pbnh/dh/100.d0
     &                +(nusk/pbrh/kpc)**2)
      if(k.eq.0) then
        coeff=1/2.d0
      else
        coeff=(dsin(k*(thetacl+0.5d0*deltathetacl))
     &         -dsin(k*(thetacl-0.5d0*deltathetacl)))/(k*deltathetacl)
      endif
      coeff=coeff/Jkplus1squared*dsbessjw(k,nusk*rcl/pbrh)*Jklocal
      if(zcl.gt.pbhg) then
        coeff=coeff/(dg*lambdag
     &     +(dh*gamh+dh*lambdah/dtanh(lambdah*(pbhh-pbhg)*kpc))
     &      /dtanh(lambdag*pbhg*kpc))
        coeff=coeff/dsinh(lambdag*pbhg*kpc)
        coeff=coeff*0.5d0*dexp(gamh*kpc*(pbhg-zcl))
        if(lambdah*kpc*pbhh.lt.37.d0) then
          coeff=coeff*
     &            dsinh(lambdah*kpc*(pbhh-zcl))
     &             /dsinh(lambdah*kpc*(pbhh-pbhg))
        else  
          coeff=coeff*dexp(lambdah*kpc*(pbhg-zcl))
        endif
      else
        coeff=coeff*0.5d0
        rebiggs=1.d0
     &  /(dg*lambdag+(dh*gamh+dh*lambdah/dtanh(lambdah*(pbhh-pbhg)*kpc))
     &  /dtanh(lambdag*pbhg*kpc)) 
        rebiggs=rebiggs/dg/lambdag*(dh*gamh
     &     +dh*lambdah/dtanh(lambdah*(pbhh-pbhg)*kpc))
        rebiggc=1.d0
     &  /(dg*lambdag+(dh*gamh+dh*lambdah/dtanh(lambdah*(pbhh-pbhg)*kpc))
     &  /dtanh(lambdag*pbhg*kpc))
        arg=lambdag*kpc*(pbhg-zcl)
        if(arg.lt.20.d0) then
          factors=dsinh(lambdag*kpc*(pbhg-zcl))
     &              /dsinh(lambdag*pbhg*kpc)
        else
          factors=dcosh(lambdag*kpc*zcl)
     &        -dsinh(lambdag*kpc*zcl)/dtanh(lambdag*pbhg*kpc)
        endif
        if(arg.lt.20.d0) then
          factorc=dcosh(lambdag*kpc*(pbhg-zcl))
     &              /dsinh(lambdag*pbhg*kpc)
        else  
          factorc=
     &      dcosh(lambdag*kpc*zcl)/dtanh(lambdag*pbhg*kpc)
     &      -dsinh(lambdag*kpc*zcl)
        endif
        coeff=coeff*(factors*rebiggs+factorc*rebiggc)
      endif
      dspbaddterm=coeff
      return
      end
