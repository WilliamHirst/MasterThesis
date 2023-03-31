      real*8 function dsepspecm(eps,how)
****************************************************************
***                                                          ***
*** function which computes the differential flux of         ***
*** positrons for the energy eps as a result of              ***
*** neutralino annihilation in the halo.                     ***
*** input: eps - positron energy in gev                      ***
***        how = 2 - diffusion model is tabulated on first   ***
***                  call, and then interpolated             ***
***              3 - as 2, but also write the table to disk  ***
***                  at the first call                       ***
***              4 - read table from disk on first call, and ***
***                  use the subsequent calls. If the file   ***
***                  does not exist, it will be created      ***
***                  (as in 3). (default)                    ***
***                                                          ***
*** rescaling is not included                                ***
***                                                          ***
*** Author: Piero Ullio                                      ***
*** Date: 2004-02-03                                         ***
*** Modified: Joakim Edsjo, modifications to file loading    ***
****************************************************************
      implicit none
      include 'dsepcom.h'
      include 'dsprep.h'
      include 'dshmcom.h'
      include 'dshacom.h'
      include 'dsge.h'
      include 'dsdirver.h'
      integer i,how
      real*8 lowloc,uploc,checkres,par
      real*8 epsabs
      integer ndigits
      common/epspecsetcom/epsabs,ndigits

      integer dhow
      common/ephow/dhow
c
      real*8 eps,v,dsepvofeps,factor,vep
      real*8 deltavprime,dsepintgreentab,kpc
      parameter(kpc=3.08567802d0) 
c
      real*8 up,low
      real*8 abserr,alist,blist,elist,epsrel,rlist,result
      integer ier,iord,last,limit,neval
      dimension alist(5000),blist(5000),elist(5000),iord(5000),
     & rlist(5000)
      external dsepspecmint
c
      real*8 vint
      common/epclintcom/vint
c

c...Check if we should be here
      if (dsepdiffmethod.ne.1) then
        write(*,*) 'ERROR in dsepspecm.f:'
        write(*,*) 'You have called this routine with dsepdiffmethod=',
     &    dsepdiffmethod
        write(*,*) 'which is not supported here. Call dsepdiff instead.'
        stop
      endif

      if (eps.gt.hamwimp) then ! outside kinematic region (JE 090109)
         dsepspecm=0.d0
         return
      endif

c...Transfer how flag to common block
      dhow=how

c...Set up integration accuracy
      if (epcl_load) then   ! first call with new diff. model
        ndigits=5
        epsabs=1.d-10
      endif  

      epsrel=10.d0**(-ndigits)
      limit=5000

c...Perform energy integration (done in v)
      v=dsepvofeps(eps)
      vint=v
c continuum contribution
      low=dsepvofeps(hamwimp)
      up=v
      par=0.d0
      do i=10,1,-1
      lowloc=low+(up-low)/10.d0*(i-1)
      uploc=low+(up-low)/10.d0*i
 20   call dqagse(dsepspecmint,lowloc,uploc,epsabs,epsrel,limit,
     &   result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      checkres=par+result
      if(dabs(checkres)/10.d0**ndigits.lt.epsabs) then
        epsabs=dabs(checkres)/10.d0**(ndigits+2) !too low accuracy
        goto 20
      elseif(dabs(checkres)/10.d0**(ndigits+2).gt.epsabs) then
        epsabs=dabs(checkres)/10.d0**(ndigits+2)  !too high accuracy
      endif
      par=par+result
      enddo
      dsepspecm=par*hasv
c sum the line contribution
      deltavprime=up-low
      DeltaVprime=deltavprime*4.0d0*k27*tau16*10.d0/kpc**2 
      dsepspecm=dsepspecm+habr(15)*hasv*dsepintgreentab(DeltaVprime)
c...dsepspecmint(low)
      factor=(rhox/hamwimp)**2*0.5d0  ! rho0-> rhox on 06-08-29
      dsepspecm=dsepspecm*factor*tau16*1.0d16/eps**2  !kpc**-3 cm**-3 GeV**-1
c...convert density to flux
      vep=2.99792458d10 !cm s**-1
c dsepspec_cl in kpc**-3 cm**-2 s**-1 GeV**-1 sr**-1
      dsepspecm=dsepspecm*vep/(4.0d0*pi)
      return
      end
