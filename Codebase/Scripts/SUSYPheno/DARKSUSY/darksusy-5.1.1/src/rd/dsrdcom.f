c      block data  dsrdcom
      subroutine  dsrdcom
      implicit none
      include 'dsrdcom.h'

      data cosmin,waccd,dpminr,dpthr,wdiffr,wdifft,
     &  hstep,hmin,compeps,xinit,xfinal,umax,cfr
     &  /0.996195d0,0.005d0,1d-4,5d-4,0.05d0,0.02d0,
     &  0.01d0,1.0d-9,0.01d0,2.0d0,200.0d0,10.0d0,0.5d0/
c...0.996195d0 - 5 degrees, 0.999048d0 - 2.5 degrees
      data thavint/1/  ! dgadap
      data rdluerr,rdlulog/6,6/
      data rdprt/0/
      data pdivr,dpres/2.0d0,0.5d0/
      data nlow,nhigh,npres,nthup,cthtest,spltest
     &  /20,10,4,4,0,1/
      end

