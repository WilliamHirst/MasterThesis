**********************************************************************
*** returns effective degrees of freedom during kinetic decoupling
***
***   input:   T  -- Temperature in MeV
***   output:  g^{1/2}, {\tilde g}^{1/2}
***
*** author: Torsten Bringmann (troms@physto.se), 2009-01-05
***         (modified version of Paolo Gondolos dsrddof)
***********************************************************************

      subroutine dsmhdof(t,sqrtg,sqrtgt)

      implicit none
      real*8 t,sqrtg,sqrtgt
      integer k,dk

      include 'dsmhcom.h'

      logical initialized
      data initialized / .false. /
      save initialized

      if (.not.initialized) then
        klo = 1
        khi = nf
        initialized = .true.
      endif

      if (t.ge.tmev(nf)) then
         sqrtg = sqrtgtab(nf)
         sqrtgt = sqrtgttab(nf)
      elseif (t.le.tmev(1)) then
         sqrtg = sqrtgtab(1)
         sqrtgt = sqrtgttab(1)
      else
 
c... This check should no longer be necessary with above
c... initialization check unless klo or khigh are modified in
c... another routine
        if (klo.lt.1.or.klo.gt.nf) klo=1
        if (khi.gt.nf.or.khi.lt.1) khi=nf

c... find an interval [tmev(klo),tmev(khi)] around t, 
c... starting with values from previous call
         dk=1
 100     if (tmev(khi).lt.t) then
            khi=khi+dk
            dk=2*dk
            if (khi.lt.nf) goto 100
            khi=nf
         endif
 110     if (tmev(klo).gt.t) then
            klo=klo-dk
            dk=2*dk
            if (klo.gt.1) goto 110
            klo=1
         endif

c... and then narrow it down
 120     if (khi-klo.gt.1) then
            k=(khi+klo)/2
            if (tmev(k).lt.t) then
               klo=k
            else
               khi=k
            endif
            goto 120
         endif
         
c... interpolate between closest tabulated points
         sqrtg = sqrtgtab(klo)+(sqrtgtab(khi)-sqrtgtab(klo))
     &          *(t-tmev(klo))/(tmev(khi)-tmev(klo))
         sqrtgt = sqrtgttab(klo)+(sqrtgttab(khi)-sqrtgttab(klo))
     &          *(t-tmev(klo))/(tmev(khi)-tmev(klo))

      endif

      return
      end
