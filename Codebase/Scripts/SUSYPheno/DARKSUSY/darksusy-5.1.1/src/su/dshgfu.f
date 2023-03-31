      subroutine dshgfu(ma,tanb,mq,mur,md,mtop,at,ab,mu,vh,
     &     stop1,stop2,v,mz,alpha1,alpha2,alpha3z)
c
c carena, quiros, wagner gfun -- adapted to darksusy by gondolo
c
      implicit real*8 (a-h,l,m,o-z)
      include 'dsge.h'
      dimension vh(2,2),vh1(2,2),vh2(2,2),vh3t(2,2),vh3b(2,2),al(2,2)
      if(dabs(mu).lt.0.000001) mu = 0.000001
      mq2 = mq**2
      mur2 = mur**2
      md2 = md**2
      tanba = tanb
      sinba = tanba/(tanba**2+1.)**.5
      cosba = sinba/tanba
      sinb = tanb/(tanb**2+1.)**.5
      cosb = sinb/tanb
      g2 = (alpha2*4.*pi)**.5
      g12 = (alpha1*4.*pi)
      g1 = g12**.5
      mw = (g2**2*v**2/2.)**.5
      alpha3 = alpha3z/(1.+23/12./pi*0.12*log(mtop**2/mz**2))
      mb = 3.
      if(mq.gt.mur) then
         mst = mq
      else
         mst = mur
      endif
      msusyt = (mst**2  + mtop**2)**.5
      if(mq.gt.md) then
         msb = mq
      else
         msb = md
      endif
      msusyb = (msb**2 + mb**2)**.5
c     if(mq.gt.md.and.mq.gt.mur) msusyb = msusyt
      tt = log(msusyt**2/mtop**2)
      tb = log(msusyb**2/mtop**2)
c     if(mq.gt.mur) t = log((mq2 + mtop**2)/mtop**2)
c     if(mur.gt.mq.or.mq.eq.mur) t = log((mur2+mtop**2)/mtop**2)
      rmtop = mtop/(1.+4.*alpha3/3./pi)
      ht = rmtop/(v*sinb)
      htst = rmtop/v
      hb = mb/v/cosb
      g32 = alpha3*4.*pi
      bt2 = -(8.*g32 - 9.*ht**2/2. - hb**2/2.)/(4.*pi)**2
      bb2 = -(8.*g32 - 9.*hb**2/2. - ht**2/2.)/(4.*pi)**2
      al2 = 3./8./pi**2*ht**2
      bt2st = -(8.*g32 - 9.*htst**2/2.)/(4.*pi)**2
      alst = 3./8./pi**2*htst**2
      al1 = 3./8./pi**2*hb**2
      al(1,1) = al1
      al(1,2) = (al2+al1)/2.
      al(2,1) = (al2+al1)/2.
      al(2,2) = al2
      mtop4 = rmtop**4.*(1.+2.*bt2*tt- al2*tt)
      mtop2 = mtop4**.5
      mbot4 = mb**4.*(1.+2.*bb2*tb - al1*tb)
      mbot2 = mbot4**.5
      if(ma.gt.mtop) then
         vi = v*(1. + 3./32./pi**2*htst**2*log(mtop**2/ma**2))
         h1i = vi* cosba
         h2i = vi*sinba
         h1t = h1i*(1.+3./8./pi**2*hb**2*log(ma**2/msusyt**2))**.25
         h2t = h2i*(1.+3./8./pi**2*ht**2*log(ma**2/msusyt**2))**.25
         h1b = h1i*(1.+3./8./pi**2*hb**2*log(ma**2/msusyb**2))**.25
         h2b = h2i*(1.+3./8./pi**2*ht**2*log(ma**2/msusyb**2))**.25
      else
         vi = v
         h1i = vi*cosb
         h2i = vi*sinb
         h1t = h1i*(1.+3./8./pi**2*hb**2*log(mtop**2/msusyt**2))**.25
         h2t = h2i*(1.+3./8./pi**2*ht**2*log(mtop**2/msusyt**2))**.25
         h1b = h1i*(1.+3./8./pi**2*hb**2*log(mtop**2/msusyb**2))**.25
         h2b = h2i*(1.+3./8./pi**2*ht**2*log(mtop**2/msusyb**2))**.25
      end if
      tanbst = h2t/h1t
      sinbt = tanbst/(1.+tanbst**2)**.5
      cosbt = sinbt/tanbst
      tanbsb = h2b/h1b
      sinbb = tanbsb/(1.+tanbsb**2)**.5
      cosbb = sinbb/tanbsb
      stop12 = (mq2 + mur2)*.5 + mtop2
     &     +1./8.*(g2**2+g1**2)*(h1t**2-h2t**2)
     &     +(((g2**2-5.*g1**2/3.)/4.*(h1t**2-h2t**2) +
     &     mq2 - mur2)**2*0.25 + mtop2*(at-mu/tanbst)**2)**.5
      stop22 = (mq2 + mur2)*.5 + mtop2
     &     +1./8.*(g2**2+g1**2)*(h1t**2-h2t**2)
     &     - (((g2**2-5.*g1**2/3.)/4.*(h1t**2-h2t**2) +
     &     mq2 - mur2)**2*0.25
     &     + mtop2*(at-mu/tanbst)**2)**.5
      if(stop22.lt.0.) then
         do i =1,2
            do j = 1,2
               vh(i,j) = -1.d+15
            enddo
         enddo
         return
      endif
      sbot12 = (mq2 + md2)*.5
     &     - 1./8.*(g2**2+g1**2)*(h1b**2-h2b**2)
     &     + (((g1**2/3.-g2**2)/4.*(h1b**2-h2b**2) +
     &     mq2 - md2)**2*0.25 + mbot2*(ab-mu*tanbsb)**2)**.5
      sbot22 = (mq2 + md2)*.5
     &     - 1./8.*(g2**2+g1**2)*(h1b**2-h2b**2)
     &     - (((g1**2/3.-g2**2)/4.*(h1b**2-h2b**2) +
     &     mq2 - md2)**2*0.25 + mbot2*(ab-mu*tanbsb)**2)**.5
      if(sbot22.lt.0.) then
         do i =1,2
            do j = 1,2
               vh(i,j) = -1.d+15
            enddo
         enddo
         return
      endif
      stop1 = stop12**.5
      stop2 = stop22**.5
      sbot1 = sbot12**.5
      sbot2 = sbot22**.5
      vh1(1,1) = 1./tanbst
      vh1(2,1) = -1.
      vh1(1,2) = -1.
      vh1(2,2) = tanbst
      vh2(1,1) = tanbst
      vh2(1,2) = -1.
      vh2(2,1) = -1.
      vh2(2,2) = 1./tanbst
cccccccccccccccccccccccccccccccc
ccc   d-terms
cccccccccccccccccccccccccccccccc
      stw=.2320
      f1t=(mq2-mur2)/(stop12-stop22)*(.5-4./3.*stw)*
     &     log(stop1/stop2)
     &     +(.5-2./3.*stw)*log(stop1*stop2/(mq2+mtop2))
     &     + 2./3.*stw*log(stop1*stop2/(mur2+mtop2))
      f1b=(mq2-md2)/(sbot12-sbot22)*(-.5+2./3.*stw)*
     &     log(sbot1/sbot2)
     &     +(-.5+1./3.*stw)*log(sbot1*sbot2/(mq2+mbot2))
     &     - 1./3.*stw*log(sbot1*sbot2/(md2+mbot2))
      f2t=mtop2**.5*(at-mu/tanbst)/(stop12-stop22)*
     &     (-.5*log(stop12/stop22)
     &     +(4./3.*stw-.5)*(mq2-mur2)/(stop12-stop22)*
     &     dshlf2(stop12,stop22))
      f2b=mbot2**.5*(ab-mu*tanbsb)/(sbot12-sbot22)*
     &     (.5*log(sbot12/sbot22)
     &     +(-2./3.*stw+.5)*(mq2-md2)/(sbot12-sbot22)*
     &     dshlf2(sbot12,sbot22))
      vh3b(1,1) = mbot4/(cosbb**2)*(log(sbot1**2*sbot2**2/
     &     (mq2+mbot2)/(md2+mbot2))
     &     + 2.*(ab*(ab-mu*tanbsb)/(sbot1**2-sbot2**2))*
     &     log(sbot1**2/sbot2**2)) +
     &     mbot4/(cosbb**2)*(ab*(ab-mu*tanbsb)/
     &     (sbot1**2-sbot2**2))**2*dshlf2(sbot12,sbot22)
      vh3t(1,1) =
     &     mtop4/(sinbt**2)*(mu*(-at+mu/tanbst)/(stop1**2
     &     -stop2**2))**2*dshlf2(stop12,stop22)
      vh3b(1,1)=vh3b(1,1)+
     &     mz**2*(2*mbot2*f1b-mbot2**.5*ab*f2b)
      vh3t(1,1) = vh3t(1,1) +
     &     mz**2*(mtop2**.5*mu/tanbst*f2t)
      vh3t(2,2) = mtop4/(sinbt**2)*(log(stop1**2*stop2**2/
     &     (mq2+mtop2)/(mur2+mtop2))
     &     + 2.*(at*(at-mu/tanbst)/(stop1**2-stop2**2))*
     &     log(stop1**2/stop2**2)) +
     &     mtop4/(sinbt**2)*(at*(at-mu/tanbst)/
     &     (stop1**2-stop2**2))**2*dshlf2(stop12,stop22)
      vh3b(2,2) =
     &     mbot4/(cosbb**2)*(mu*(-ab+mu*tanbsb)/(sbot1**2
     &     -sbot2**2))**2*dshlf2(sbot12,sbot22)
      vh3t(2,2)=vh3t(2,2)+
     &     mz**2*(-2*mtop2*f1t+mtop2**.5*at*f2t)
      vh3b(2,2) = vh3b(2,2) -mz**2*mbot2**.5*mu*tanbsb*f2b
      vh3t(1,2) = -
     &     mtop4/(sinbt**2)*mu*(at-mu/tanbst)/
     &     (stop1**2-stop2**2)*(log(stop1**2/stop2**2) + at*
     &     (at - mu/tanbst)/(stop1**2-stop2**2)*dshlf2(stop12,stop22))
      vh3b(1,2) =
     &     - mbot4/(cosbb**2)*mu*(at-mu*tanbsb)/
     &     (sbot1**2-sbot2**2)*(log(sbot1**2/sbot2**2) + ab*
     &     (ab - mu*tanbsb)/(sbot1**2-sbot2**2)*dshlf2(sbot12,sbot22))
      vh3t(1,2)=vh3t(1,2) +
     &     mz**2*(mtop2/tanbst*f1t-mtop2**.5*(at/tanbst+mu)/2.*f2t)
      vh3b(1,2)=vh3b(1,2)
     &     +mz**2*(-mbot2*tanbsb*f1b+mbot2**.5*(ab*tanbsb+mu)/2.*f2b)
      vh3t(2,1) = vh3t(1,2)
      vh3b(2,1) = vh3b(1,2)
      tq = log((mq2 + mtop2)/mtop2)
      tu = log((mur2+mtop2)/mtop2)
      tqd = log((mq2 + mb**2)/mb**2)
      td = log((md2+mb**2)/mb**2)
      do i = 1,2
         do j = 1,2
            vh(i,j) =
     &           6./(8.*pi**2*(h1t**2+h2t**2))
     &           *vh3t(i,j)*0.5*(1.-al(i,j)*tt/2.) +
     &           6./(8.*pi**2*(h1b**2+h2b**2))
     &           *vh3b(i,j)*0.5*(1.-al(i,j)*tb/2.)
         enddo
      enddo
      return
      end
