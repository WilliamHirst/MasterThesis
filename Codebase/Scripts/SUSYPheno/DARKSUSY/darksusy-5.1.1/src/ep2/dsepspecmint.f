      real*8 function dsepspecmint(vprime)
      implicit none
      include 'dshacom.h'
      include 'dsepcom.h'
      real*8 vprime
      real*8 dshaloyield,epsprime,deltavprime,dsepintgreentab,
     &  dsepepsofv,dsepwofu,uprime,dsepuofeps
      integer istat
      real*8 vint
      common/epclintcom/vint
      real*8 kpc
      parameter(kpc=3.08567802d0) 
      deltavprime=vint-vprime
      DeltaVprime=deltavprime*4.0d0*k27*tau16*10.d0/kpc**2 
      dsepspecmint=dsepintgreentab(DeltaVprime)
      epsprime=dsepepsofv(vprime)
      uprime=dsepuofeps(epsprime)
c dsepspecpuint in GeV**-1
      hasmooth=2  ! 2 to smooth well
      dsepspecmint=dsepspecmint*dsepwofu(uprime)
     &   *dshaloyield(epsprime,151,istat)
      hasmooth=0
      return
      end
