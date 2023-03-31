      subroutine dspole(mchi,ma,tanb,mq,mur,mdr,
     &     mtop,at,ab,mu,mh,mhp,hm,hmp,amp,sa,ca,
     &     v,mz,alpha1,alpha2,alpha3z,sint,lambda,prtlevel,ierr)
c
c carena, quiros, wagner -- adapted to darksusy by gondolo
c
      implicit real*8(a-h,m,o-z)
      include 'dsge.h'
      real*8 lambda(7)
      integer prtlevel
      dimension delta(2,2),coupt(2,2),t(2,2),sstop2(2),
     &     ssbot2(2),b(2,2),coupb(2,2),
     &     hcoupt(2,2),hcoupb(2,2),
     &     acoupt(2,2),acoupb(2,2),pr(3), polar(3)

      ihiggs = 3

      ierr = 0
      delta(1,1) = 1.d0
      delta(2,2) = 1.d0
      delta(1,2) = 0.d0
      delta(2,1) = 0.d0
      alpha3=1.d0/(1.d0/alpha3z+23.d0/6.d0/pi*log(mtop/mz))

      rmtop = mtop/(1.d0+4.d0*alpha3/3.d0/pi)

      ht = rmtop /v
      call dsrghm(mchi,ma,tanb,mq,mur,mdr,mtop,at,ab,
     &     mu,mh,hm,sa,ca,tanba,v,mz,alpha1,alpha2,alpha3z,
     &     lambda,prtlevel)
      sinb = tanb/dsqrt(tanb**2+1.d0)
      cosb = 1.d0/dsqrt(tanb**2+1.d0)
      cos2b = sinb**2 - cosb**2
      sinbpa = sinb*ca + cosb*sa
      cosbpa = cosb*ca - sinb*sa
      rmbot = 3.d0
      mq2 = mq**2
      mur2 = mur**2
      mdr2 = mdr**2
      mst11 = rmtop**2 + mq2  - 0.35d0*mz**2*cos2b
      mst22 = rmtop**2 + mur2 - 0.15d0*mz**2*cos2b
      msb11 = rmbot**2 + mq2  + 0.42*mz**2*cos2b
      msb22 = rmbot**2 + mdr2 + 0.08*mz**2*cos2b
      wmst11 = rmtop**2 + mq2
      wmst22 = rmtop**2 + mur2
      mst12 = rmtop*(at - mu/tanb)
      msb12 = rmbot*(ab - mu*tanb)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     stop eigenvalues calculation
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      stop12 = 0.5*(mst11+mst22) +
     &     0.5*((mst11+mst22)**2 -
     &     4.*(mst11*mst22 - mst12**2))**.5
      stop22 = 0.5*(mst11+mst22) -
     &     0.5*((mst11+mst22)**2 - 4.*(mst11*mst22 - mst12**2))**.5
      if(stop22.lt.0.d0) then
         ierr=1
         call dswrite(1,1,'mstop^2<0 in dspole')
         return
      endif
      sstop2(1) = stop12
      sstop2(2) = stop22
      stop1 = stop12**.5
      stop2 = stop22**.5
      if(mst12.eq.0.) xst11 = 1.
      if(mst12.eq.0.) xst12 = 0.
      if(mst12.eq.0.) xst21 = 0.
      if(mst12.eq.0.) xst22 = 1.
      if(mst12.ne.0.) then
         xst11 = mst12/(mst12**2+(mst11-stop12)**2)**.5
         xst12 = - (mst11-stop12)/(mst12**2+(mst11-stop12)**2)**.5
         xst21 = mst12/(mst12**2+(mst11-stop22)**2)**.5
         xst22 = - (mst11-stop22)/(mst12**2+(mst11-stop22)**2)**.5
      endif
      t(1,1) = xst11
      t(2,2) = xst22
      t(1,2) = xst12
      t(2,1) = xst21
      sbot12 = 0.5*(msb11+msb22) +
     &     0.5*((msb11+msb22)**2 -
     &     4.*(msb11*msb22 - msb12**2))**.5
      sbot22 = 0.5*(msb11+msb22) -
     &     0.5*((msb11+msb22)**2 - 4.*(msb11*msb22 - msb12**2))**.5
      if(sbot22.lt.0.d0) then
         ierr=2
         call dswrite(1,1,'msbottom^2<0 in dspole')
         return
      endif
      sbot1 = sbot12**.5
      sbot2 = sbot22**.5
      ssbot2(1) = sbot12
      ssbot2(2) = sbot22
      if(msb12.eq.0.) xsb11 = 1.
      if(msb12.eq.0.) xsb12 = 0.
      if(msb12.eq.0.) xsb21 = 0.
      if(msb12.eq.0.) xsb22 = 1.
      if(msb12.ne.0.) then
         xsb11 = msb12/(msb12**2+(msb11-sbot12)**2)**.5
         xsb12 = - (msb11-sbot12)/(msb12**2+(msb11-sbot12)**2)**.5
         xsb21 = msb12/(msb12**2+(msb11-sbot22)**2)**.5
         xsb22 = - (msb11-sbot22)/(msb12**2+(msb11-sbot22)**2)**.5
      endif
      b(1,1) = xsb11
      b(2,2) = xsb22
      b(1,2) = xsb12
      b(2,1) = xsb21

      sqr = dsqrt(2.d0)
      vp = v*sqr

cccccccccccccccccccccccccccccccccccccccc

      if(ihiggs.eq.0) return

cccccccccccccccccccccccccccccccccccc
ccc   starting of light higgs
cccccccccccccccccccccccccccccccccccc

      do i = 1,2
         do j = 1,2
            coupt(i,j) =
     &           sint*mz**2*2.*sqr/v/3.*sinbpa*(delta(i,j) +
     &           (3. - 8.*sint)/4./sint*t(1,i)*t(1,j))
     &           -rmtop**2/v**2*vp/sinb*ca*delta(i,j)
     &           -rmtop/vp/sinb*(at*ca + mu*sa)*(t(1,i)*t(2,j) +
     &           t(1,j)*t(2,i))
         enddo
      enddo

      do i = 1,2
         do j = 1,2
            coupb(i,j) =
     &           -sint*mz**2*2.*sqr/v/6.*sinbpa*(delta(i,j) +
     &           (3. - 4.*sint)/2./sint*b(1,i)*b(1,j))
     &           +rmbot**2/v**2*vp/cosb*sa*delta(i,j)
     &           +rmbot/vp/cosb*(ab*sa + mu*ca)*(b(1,i)*b(2,j) +
     &           b(1,j)*b(2,i))
         enddo
      enddo

      prun = mh
      epss = 1.d-4*prun
      iter = 0
 7007 iter = iter + 1
      do i3 = 1,3
         pr(i3)=prun+(i3-2)*epss/2
         p2=pr(i3)**2
         polt = 0.
         do i = 1,2
            do j = 1,2
               polt = polt + coupt(i,j)**2*3.*
     &              dshlf3(p2,sstop2(i),sstop2(j))/16./pi**2
            enddo
         enddo
         polb = 0.
         do i = 1,2
            do j = 1,2
               polb = polb + coupb(i,j)**2*3.*
     &              dshlf3(p2,ssbot2(i),ssbot2(j))/16./pi**2
            enddo
         enddo
         rmtop2 = rmtop**2
         mtop2=mtop**2
c
         poltt =
     &        3.*rmtop**2/8./pi**2/  v  **2*
     &        ca**2/sinb**2 *
     &        (-2.*mtop**2+.5*p2)*
     &        dshlf3(p2,mtop2,mtop2)
c
         pol = polt + polb + poltt
         polar(i3) = p2 - mh**2 - pol
      enddo
      deriv = (polar(3)-polar(1))/epss
      drun = - polar(2)/deriv
      prun = prun + drun
      p2 = prun**2
      if ( abs(drun) .ge. 1.d-4 ) goto 7007

      mhp = p2**.5

cccccccccccccccccccccccccccccccccccccccc
ccc   end of light higgs
cccccccccccccccccccccccccccccccccccccccc

      if(ihiggs.eq.1) return

ccccccccccccccccccccccccccccccccccccccccc
ccc   starting of heavy higgs
cccccccccccccccccccccccccccccccccccccccccc

      do i = 1,2
         do j = 1,2
            hcoupt(i,j) =
     &           -sint*mz**2*2.*sqr/v/3.*cosbpa*(delta(i,j) +
     &           (3. - 8.*sint)/4./sint*t(1,i)*t(1,j))
     &           -rmtop**2/v**2*vp/sinb*sa*delta(i,j)
     &           -rmtop/vp/sinb*(at*sa - mu*ca)*(t(1,i)*t(2,j) +
     &           t(1,j)*t(2,i))
         enddo
      enddo

      do i = 1,2
         do j = 1,2
            hcoupb(i,j) =
     &           sint*mz**2*2.*sqr/v/6.*cosbpa*(delta(i,j) +
     &           (3. - 4.*sint)/2./sint*b(1,i)*b(1,j))
     &           -rmbot**2/v**2*vp/cosb*ca*delta(i,j)
     &           -rmbot/vp/cosb*(ab*ca - mu*sa)*(b(1,i)*b(2,j) +
     &           b(1,j)*b(2,i))
            hcoupb(i,j)=0.
         enddo
      enddo

      prun = hm
      epss = 1.d-4*prun
      iter = 0
 1001 iter = iter + 1
      do i3 = 1,3
         pr(i3)=prun+(i3-2)*epss/2
         hp2=pr(i3)**2
         hpolt = 0.
         do i = 1,2
            do j = 1,2
               hpolt = hpolt + hcoupt(i,j)**2*3.*
     &              dshlf3(hp2,sstop2(i),sstop2(j))/16./pi**2
            enddo
         enddo
         hpolb = 0.
         do i = 1,2
            do j = 1,2
               hpolb = hpolb + hcoupb(i,j)**2*3.*
     &              dshlf3(hp2,ssbot2(i),ssbot2(j))/16./pi**2
            enddo
         enddo
         rmtop2 = rmtop**2
         mtop2  = mtop**2
         hpoltt =
     &        3.*rmtop**2/8./pi**2/  v  **2*
     &        sa**2/sinb**2 *
     &        (-2.*mtop**2+.5*hp2)*
     &        dshlf3(hp2,mtop2,mtop2)
         hpol = hpolt + hpolb + hpoltt
         polar(i3) =hp2-hm**2-hpol
      enddo
      deriv = (polar(3)-polar(1))/epss
      drun = - polar(2)/deriv
      prun = prun + drun
      hp2 = prun**2
      if ( abs(drun) .ge. 1.d-4 ) goto 1001

      hmp = hp2**.5

ccccccccccccccccccccccccccccccccccccccccccc
ccc   end of heavy higgs
cccccccccccccccccccccccccccccccccccccccccccc

      if(ihiggs.eq.2) return

cccccccccccccccccccccccccccccccccccccccccccc
ccc   beginning of pseudoscalar higgs
cccccccccccccccccccccccccccccccccccccccccccc

      do i = 1,2
         do j = 1,2
            acoupt(i,j) =
     &           -rmtop/vp/sinb*(at*cosb + mu*sinb)*
     &           (t(1,i)*t(2,j) -t(1,j)*t(2,i))
         enddo
      enddo
      do i = 1,2
         do j = 1,2
            acoupb(i,j) =
     &           rmbot/vp/cosb*(ab*sinb + mu*cosb)*
     &           (b(1,i)*b(2,j) -b(1,j)*b(2,i))
         enddo
      enddo
      prun = ma
      epss = 1.d-4*prun
      iter = 0
 6006 iter = iter + 1
      do i3 = 1,3
         pr(i3)=prun+(i3-2)*epss/2
         ap2=pr(i3)**2
         apolt = 0.
         do i = 1,2
            do j = 1,2
               apolt = apolt + acoupt(i,j)**2*3.*
     &              dshlf3(ap2,sstop2(i),sstop2(j))/16./pi**2
            enddo
         enddo
         apolb = 0.
         do i = 1,2
            do j = 1,2
               apolb = apolb + acoupb(i,j)**2*3.*
     &              dshlf3(ap2,ssbot2(i),ssbot2(j))/16./pi**2
            enddo
         enddo
         rmtop2 = rmtop**2
         mtop2=mtop**2
         apoltt =
     &        3.*rmtop**2/8./pi**2/  v  **2*
     &        cosb**2/sinb**2 *
     &        (-.5*ap2)*
     &        dshlf3(ap2,mtop2,mtop2)
         apol = apolt + apolb + apoltt
         polar(i3) = ap2 - ma**2 -apol
      enddo
      deriv = (polar(3)-polar(1))/epss
      drun = - polar(2)/deriv
      prun = prun + drun
      ap2 = prun**2
      if ( abs(drun) .ge. 1.d-4 ) goto 6006

      amp = ap2**.5

ccccccccccccccccccccccccccccccccccccccccccc
ccc   end of pseudoscalar higgs
cccccccccccccccccccccccccccccccccccccccccccc

      if(ihiggs.eq.3) return

      end
