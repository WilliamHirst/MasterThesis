**********************************************************************
*** function which approximates the function dspbtd15comp by 
*** estimating that diffusion time term supposing to have a point 
*** source located at the galactic center but then weighting it with 
*** the emission over a whole cylinder of radius scale and 
*** height 2*scale, i.e. rho2int (to be given in kpc^3). 
*** The goodness of the approximation should be checked by comparing
***    dspbtd15point(rho2int,tp) with
***    dspbtd15comp(tp) for different value of tp and scale, 
*** and depending on the halo profile chosen and level of precision 
*** required. The comparison has to be performed but setting 
*** rho2int=dshmrho2cylint(scale,scale) and each 
*** pbrcy and pbzcy pair equal to scalebefore calling dspbtd15comp, 
*** possibly resetting the parameter clspset as well, see the header
*** of the function dspbtd15beuclsp
*** 
*** input: rho2int=dshmrho2cylint(scale,scale) in kpc^3, tp in GeV
*** output in 10^15 s
***
*** author: piero ullio (ullio@sissa.it)
*** date: 04-01-22
**********************************************************************

      real*8 function dspbtd15point(rho2int,tp)
      implicit none
      include 'dshmcom.h'
      include 'dspbcom.h'
      real*8 rho2int,tp,dspbtd15beucl,clumpterm,thetaL,L,Lhor
c
      zcl=0.05d0
      thetaL=2.d0*datan(1.d0)
      deltathetacl=datan(1.d0)/100.d0
      L=pbr0
      Lhor=dsqrt(L**2-zcl**2)
      rcl=dsqrt(Lhor**2+pbr0**2-2.d0*pbr0*Lhor*dcos(thetaL))
      thetacl=dasin(Lhor/rcl*dsin(thetaL))
      if(Lhor*dcos(thetaL).gt.pbr0.and.dcos(thetaL).gt.0.d0) then
        thetacl=4.d0*datan(1.d0)-thetacl
      endif
c
      clumpterm=dspbtd15beucl(tp)
      dspbtd15point=rho2int*clumpterm
      return
      end
