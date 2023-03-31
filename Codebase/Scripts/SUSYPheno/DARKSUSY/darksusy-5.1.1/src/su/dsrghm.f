      subroutine dsrghm(mchi,ma,tanb,mq,mur,md,mtop,au,ad,mu,
     &     mhp,hmp,sa,ca,tanba,v,mz,alpha1,alpha2,alpha3z,
     &     lambda,prtlevel)
c
c carena, quiros, wagner -- adapted to darksusy by gondolo
c
      implicit real*8(a-h,l,m,o-z)
      include 'dsge.h'
      real*8 lambda(7)
      integer prtlevel
      dimension vh(2,2),m2(2,2),m2p(2,2)

      tanba = tanb
      tanbt = tanb

c     mbottom(mtop) = 3. gev
      mb = 3.
      alpha3 = alpha3z/(1. +(11. - 10./3.)/4./pi*alpha3z*
     &     log(mtop**2/mz**2))

c     rmtop= running top quark mass
      rmtop = mtop/(1.+4.*alpha3/3./pi)
      tq = log((mq**2+mtop**2)/mtop**2)
      tu = log((mur**2 + mtop**2)/mtop**2)
      td = log((md**2 + mtop**2)/mtop**2)
      sinb = tanb/((1. + tanb**2)**.5)
      cosb = sinb/tanb
      if(ma.gt.mtop)
     &     tanba = tanb*(1.-3./32./pi**2*
     &     (rmtop**2/v**2/sinb**2-mb**2/v**2/cosb**2)*
     &     log(ma**2/mtop**2))
      if(ma.lt.mtop.or.ma.eq.mtop) tanbt = tanba
      sinb = tanbt/((1. + tanbt**2)**.5)
      cosb = 1./((1. + tanbt**2)**.5)
      cos2b = (tanbt**2 - 1.)/(tanbt**2 + 1.)
      g1 = (alpha1*4.*pi)**.5
      g2 = (alpha2*4.*pi)**.5
      g3 = (alpha3*4.*pi)**.5
      hu = rmtop/v/sinb
      hd =  mb/v/cosb

      call dshgfu(ma,tanba,mq,mur,md,mtop,au,ad,mu,vh,stop1,stop2,
     &     v,mz,alpha1,alpha2,alpha3z)

      if(mq.gt.mur) then
         tp = tq - tu
      else
         tp = tu - tq
      endif
      if(mq.gt.mur) then
         tdp = tu
      else
         tdp = tq
      endif
      if(mq.gt.md) then
         tpd = tq - td
      else
         tpd = td - tq
      endif
      if(mq.gt.md) then
         tdpd = td
      else
         tdpd = tq
      endif

      if(mq.gt.md) then
         dlambda1 = 6./96./pi**2*g1**2*hd**2*tpd
      else
         dlambda1 = 3./32./pi**2*hd**2*(g1**2/3.+g2**2)*tpd
      endif

      if(mq.gt.mur) then
         dlambda2 =12./96./pi**2*g1**2*hu**2*tp
      else
         dlambda2 = 3./32./pi**2*hu**2*(-g1**2/3.+g2**2)*tp
      endif

      dlambda3 = 0.
      dlambda4 = 0.

      if(mq.gt.md) then
         dlambda3 = -1./32./pi**2*g1**2*hd**2*tpd
      else
         dlambda3 = 3./64./pi**2*hd**2*(g2**2-g1**2/3.)*tpd
      endif


      if(mq.gt.mur) then
         dlambda3 = dlambda3 - 1./16./pi**2*g1**2*hu**2*tp
      else
         dlambda3 = dlambda3 + 3./64./pi**2*hu**2*(g2**2+g1**2/3.)*tp
      endif

      if(mq.lt.mur) then
         dlambda4 = -3./32./pi**2*g2**2*hu**2*tp
      else
         dlambda4 = dlambda4 - 3./32./pi**2*g2**2*hd**2*tpd
      endif

      lambda(1) = ((g1**2 + g2**2)/4.)*
     &     (1.-3.*hd**2*(tpd + tdpd)/8./pi**2)
     &     +(3.*hd**4./16./pi**2) *tpd*(1.
     &     + (3.*hd**2/2. + hu**2/2.
     &     - 8.*g3**2) * (tpd + 2.*tdpd)/16./pi**2)
     &     +(3.*hd**4./8./pi**2) *tdpd*(1.  + (3.*hd**2/2. + hu**2/2.
     &     - 8.*g3**2) * tdpd/16./pi**2) + dlambda1
      lambda(2) = ((g1**2 + g2**2)/4.)*(1.-3.*hu**2*
     &     (tp + tdp)/8./pi**2)
     &     +(3.*hu**4./16./pi**2) *tp*(1.
     &     + (3.*hu**2/2. + hd**2/2.
     &     - 8.*g3**2) * (tp + 2.*tdp)/16./pi**2)
     &     +(3.*hu**4./8./pi**2) *tdp*(1. + (3.*hu**2/2. + hd**2/2.
     &     - 8.*g3**2) * tdp/16./pi**2) + dlambda2
      lambda(3) = ((g2**2 - g1**2)/4.)*(1.-3.*
     &     (hu**2)*(tp + tdp)/16./pi**2 -3.*
     &     (hd**2)*(tpd + tdpd)/16./pi**2) +dlambda3
      lambda(4) = (- g2**2/2.)*(1.
     &     -3.*(hu**2)*(tp + tdp)/16./pi**2
     &     -3.*(hd**2)*(tpd + tdpd)/16./pi**2) +dlambda4

      lambda(5) = 0.
      lambda(6) = 0.
      lambda(7) = 0.

      m2(1,1) = 2.*v**2*(lambda(1)*cosb**2+2.*lambda(6)*
     &     cosb*sinb + lambda(5)*sinb**2) + ma**2*sinb**2

      m2(2,2) = 2.*v**2*(lambda(5)*cosb**2+2.*lambda(7)*
     &     cosb*sinb + lambda(2)*sinb**2) + ma**2*cosb**2
      m2(1,2) = 2.*v**2*(lambda(6)*cosb**2+(lambda(3)+lambda(4))*
     &     cosb*sinb + lambda(7)*sinb**2) - ma**2*sinb*cosb

      m2(2,1) = m2(1,2)
ccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   this is the contribution from light charginos/neutralinos
ccccccccccccccccccccccccccccccccccccccccccccccccc

      mssusy=(.5*(mq**2+mur**2)+mtop**2)**.5

      if(mchi.gt.mssusy)goto 3790
      if(mchi.lt.mtop) mchi=mtop

      tchar=log(mssusy**2/mchi**2)

      deltal12=(9./64./pi**2*g2**4+5./192./pi**2*g1**4)*tchar
      deltal3p4=(3./64./pi**2*g2**4+7./192./pi**2*g1**4
     &     +4./32/pi**2*g1**2*g2**2)*tchar

      deltam112=2.*deltal12*v**2*cosb**2
      deltam222=2.*deltal12*v**2*sinb**2
      deltam122=2.*deltal3p4*v**2*sinb*cosb

      m2(1,1)=m2(1,1)+deltam112
      m2(2,2)=m2(2,2)+deltam222
      m2(1,2)=m2(1,2)+deltam122
      m2(2,1)=m2(2,1)+deltam122

 3790 continue

ccccccccccccccccccccccccccccccccccccccccccc
ccc   end of charginos/neutralinos
ccccccccccccccccccccccccccccccccccccccccccc

      do 9800 i = 1,2
         do 9801 j = 1,2
            m2p(i,j) = m2(i,j) + vh(i,j)
 9801    continue
 9800 continue

      if (prtlevel.ge.2) then
         write (*,'('' higgs mass^2 matrix = '',2(f12.2))')
     &        m2p(1,1),m2p(1,2)
         write (*,'(''                       '',2(f12.2))')
     &        m2p(2,1),m2p(2,2)
      endif

      trm2p = m2p(1,1) + m2p(2,2)
      detm2p = m2p(1,1)*m2p(2,2) - m2p(1,2)*m2p(2,1)

      mh2p = (trm2p - (trm2p**2 - 4.* detm2p)**.5)/2.
      hm2p = (trm2p + (trm2p**2 - 4.* detm2p)**.5)/2.
      hmp = hm2p**.5
      if(mh2p.lt.0.) goto 5555
      mhp = mh2p**.5
      sin2alpha = 2.*m2p(1,2)/(trm2p**2-4.*detm2p)**.5
      cos2alpha = (m2p(1,1)-m2p(2,2))/(trm2p**2-4.*detm2p)**.5
      if(cos2alpha.gt.0.) then
         alpha = asin(sin2alpha)/2.
      else if (cos2alpha.lt.0.) then
         alpha = -pi/2.-asin(sin2alpha)/2.
      else
         alpha = 0.d0
      endif
      sa = sin(alpha)
      ca = cos(alpha)
      sqbma = (sinb*ca - cosb*sa)**2
 5555 xin = 1.
 2242 return
      end
