***********************************************************************
*** The function dsmhmcut returns the typical mass scale of the smallest
*** gravitationally bound objects [in units of M_sun]
***
***  input: m0  - DM mass [in GeV]
***         tkd - kinetic decoupling temperature [in MeV]
***         how = 1,2,3
***               1 - maximum of 2,3 [default]
***               2 - horizon mass at T_kd
***               3 - cutoff scale associated to free-streaming
***
*** author: torsten bringmann (troms@physto.se), 2010-01-23
***********************************************************************

      real*8 function dsmhmcut(m0,tkd,how)
      implicit none

      include 'dsmhcom.h'

      integer how
      real*8  m0,tkd

      real*8 tmpres1,tmpres2,sqrtg,sqrtgt
      real*8 pi
      parameter (pi=3.141592653589793238d0)

      tmpres1=0d0
      tmpres2=0d0

      call dsmhdof(tkd,sqrtg,sqrtgt)
      tmpres1=3.38D-6*(50d0/tkd/sqrt(sqrtg))**3
      tmpres2=2.90d-6*(1+log(sqrt(sqrtg)*tkd/30.)/18.56)**3/
     -              (m0/100.*sqrtg*tkd/30.)**1.5d0

      if (how.eq.2) dsmhmcut=tmpres1
      if (how.eq.3) dsmhmcut=tmpres2
      if (how.ne.2.and.how.ne.3) then    ! default
        if (tmpres1.ge.tmpres2) dsmhmcut=tmpres1
        if (tmpres1.lt.tmpres2) dsmhmcut=tmpres2
      endif


      return
      end





        

