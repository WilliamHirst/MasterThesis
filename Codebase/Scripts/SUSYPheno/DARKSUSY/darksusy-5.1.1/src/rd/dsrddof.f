      subroutine dsrddof(t,sqrtg,gentr,genergy,gpressure)
c_______________________________________________________________________
c     return effective degrees of freedom at temperature t
c     common:
c       'dsrdcom.h' - included common blocks
c  author: paolo gondolo (paolo@physics.utah.edu) 2005
c  modified: paolo gondolo (paolo@physics.utah.edu) 2008
c=======================================================================
      implicit none
      real*8 t,sqrtg,gentr,genergy,gpressure
      include 'dsrdcom.h'
      integer k,kk

c short version without genergy, gpressure
      genergy=0.d0
      gpressure=0.d0

      if (t.ge.tgev(1)) then

         sqrtg = fg(1)
         gentr = fh(1)
c         genergy = fe(1)
c         gpressure = fp(1)

      else

         kk=1
 100     if (tgev(khi).gt.t) then
            khi=khi+kk
            kk=2*kk
            if (khi.le.nf) goto 100
            khi=nf
         endif
 110     if (tgev(klo).lt.t) then
            klo=klo-kk
            kk=2*kk
            if (klo.ge.1) goto 110
            klo=1
         endif
 120     if (khi-klo.gt.1) then
            k=(khi+klo)/2
            if (tgev(k).lt.t) then
               khi=k
            else
               klo=k
            endif
            goto 120
         endif
         
         sqrtg = (t*(fg(khi)-fg(klo))+
     &        fg(klo)*tgev(khi)-fg(khi)*tgev(klo))/(tgev(khi)-tgev(klo))
         gentr = (t*(fh(khi)-fh(klo))+
     &        fh(klo)*tgev(khi)-fh(khi)*tgev(klo))/(tgev(khi)-tgev(klo))
c         genergy = (t*(fe(khi)-fe(klo))+
c     &        fe(klo)*tgev(khi)-fe(khi)*tgev(klo))/(tgev(khi)-tgev(klo))
c         gpressure = (t*(fp(khi)-fp(klo))+
c     &        fp(klo)*tgev(khi)-fp(khi)*tgev(klo))/(tgev(khi)-tgev(klo))

      endif

      end
