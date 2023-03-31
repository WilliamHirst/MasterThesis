***********************************************************************
*** auxialiary function for the quick determination of Tkd
*** (integrand for the thermal average). Called by dsmhTkd.
***
*** author: torsten bringmann (troms@physto.se), 2010-01-23
***********************************************************************

      real*8 function dsmhcint2(y)
      implicit none

      include 'dsmssm.h'
      include 'dsmhcom.h'
      include 'dsidtag.h'

      real*8 m0,mf,x,y,tmpres
      real*8 cn,dsmhm2
      integer SMtype,n

      SMtype = 1
      tmpres = 0d0
      dsmhcint2=0d0    

      m0=mass(kn(1))

 10   if (SMtype.le.3) mf=0d0
      if (SMtype.eq.4) mf=mass(ke)
      if (SMtype.eq.5) mf=mass(kmu)
      if (SMtype.eq.6) mf=mass(ktau)
      if (y.lt.mf/Tint/10.) goto 20
      x = sqrt(y**2+(mf/Tint)**2)

      call dsmhm2simp(m0,SMtype,cn,n)
      tmpres = tmpres+cn*y**2*x**n/(1.+exp(x))

 20   SMtype = SMtype+1
      if (SMtype.le.6) goto 10


      dsmhcint2 = tmpres
      
      return
      end

