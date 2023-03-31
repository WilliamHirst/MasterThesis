      subroutine dshigsct(unphys,hwarning)
c_______________________________________________________________________
c  higgs bosons masses and mixings
c  called by dssusy.
c  needs dssfesct.
c  higloop =  0  tree-level
c             1  brignole-ellis-ridolfi-zwirner eff. pot.
c             2  drees-nojiri eff. pot.
c             3  carena-espinosa-quiros-wagner rg-impr. eff. pot.
c                    (uses subh.f.)
c             4  carena-quiros-wagner impr. eff. pot.
c                    (uses subhpole2.f)
c             5  use FeynHiggs by Heinemeyer, Hollik and Weiglein
c                requires full FeynHiggs to be installed (see below)
c             6  use FeynHiggsFast by Heinemeyer, Hollik and Weiglein
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994,1995
c  modified by: joakim edsjo, edsjo@physto.se, 2000-09-01
c  modified by: paolo gondolo, 1999-2000
c=======================================================================
      implicit none
      include 'dsmssm.h'
      include 'dsio.h'
      real*8 cos2al,sin2al,atop,abot,mt2,mb2,mt4,mb4,ltop,gtop,rtop,
     &     lbot,gbot,rbot,d11,d12,d22,coeff,a,b,c,d,msqh2,msqh1,msqz,
     &     msqw,msqa,msqst1,msqst2,msqsb1,msqsb2,msqhc,cb,sb,cb2,sb2,
     &     cb3,sb3,cb4,sb4,c2b,tw2,mq2,md2,mu2,g,b012,bd12,b12,c012,
     &     cd12,c12,d012,dd12,h1,h2,h3,h4,q2,ftop,fbot,deltasq,ht,hb,
     &     ub(2,2),ut(2,2),x(2,2),y(2,2),cc(2,2),dd(2,2),dsg0loop,
     &     mst2(2),msb2(2),delta0,ht2,hb2,au,ad,mqtl,mqtr,mqbl,mqbr,
     &     msusy2,mchi,mhp,hmp,mdr,amp,alpha1,
     &     alpha2,alpha3,alpha3z,pi2,rmbot,rmtop,ms2,t,vsq,g1sq,g2sq,
     &     g3sq,hu2,hd2,hu4,hd4,musq,xau,xad,aud,ms4,v1,v2,lambda(7),
     &     mhch,halpha
      integer i,j,k,l,ierr,unphys,hwarning,msbarselec
      real*8 tanb,mq,mur,mtop,mh,hm,sa,ca,v,mz ! for dspole
      real*8 pi
      parameter (pi=3.141592653589793238d0)
      character*180 message
      character*8 dsi2s

      hwarning=0          ! warnings from the Higgs calculation
      delrho=0.0d0        ! MSSM correction to rho parameter

      atop = -asoftu(3)   ! je convention change 960525
      abot = -asoftd(3)   ! je convention change 960525
      msqa = ma**2
      msqz = mass(kz)**2
      msqw = mass(kw)**2
      cb = cosbe
      cb2 = cosbe*cb
      cb3 = cosbe*cb2
      cb4 = cosbe*cb3
      sb = sinbe
      sb2 = sinbe*sb
      sb3 = sinbe*sb2
      sb4 = sinbe*sb3
      c2b = cb2-sb2
      tw2 = (sinthw/costhw)**2
      mt2 = mass(kt)**2
      mb2 = mass(kb)**2
      mq2 = mass2q(3)
      md2 = mass2d(3)
      mu2 = mass2u(3)
      v = sqrt(2.0d0)*mass(kw)/g2weak
      v1 = v*cosbe
      v2 = v*sinbe
      vsq = v*v
      lam1 = 0.25d0*(g2weak**2+gyweak**2)
      lam2 = lam1
      lam3 = 0.25d0*(g2weak**2-gyweak**2)
      lam4 = -0.5d0*g2weak**2
      lam5 = 0.d0
      lam6 = 0.d0
      lam7 = 0.d0

c----------------------------------------------- tree level higgs sector
      if (higloop.eq.0) then
        if (ma.eq.0.d0) then
          mass(kh1) = mass(kz)
          mass(kh2) = 0.0d0
          mass(kh3) = 0.0d0
          mass(khc) = mass(kw)
          alpha = -atan(tanbe)
        else
c        r3 = r2 * (1.0d0-r2) / (cos2be*cos2be-r2)
c          r3 = msqa/msqz
c          auxa = r3+1.0d0
c          auxb = auxa*auxa-4.0d0*r3*cos2be**2
c          if (auxb .le. 0.0d0) then
c            unphys = -2
c            return
c          else
c            r2 = 0.5d0*(auxa-sqrt(auxb))
c            r1 = r3/r2 * cos2be**2
c            mass(kh2) = mass(kz) * sqrt(r2)
c            mass(kh1) = mass(kz) * sqrt(r1)
c            mass(khc) = mass(kz) * sqrt(r3+costhw**2)
c            cos2al = cos2be *
c     &        ( (cos2be**2+r2**2-2.0d0*r2) /
c     &        (cos2be**2+r2**2-2.0d0*r2*cos2be**2) )
c            alpha = -0.5d0*acos(cos2al)
c          endif
          a = msqz*cb2 + msqa*sb2
          b = -sb*cb*(msqz+msqa)
          c = msqz*sb2 + msqa*cb2
          if (prtlevel.ge.2) then
             write (*,'('' higgs mass^2 matrix = '',2(f12.2))')
     &            a,b
             write (*,'(''                       '',2(f12.2))')
     &            b,c
          endif
          d = sqrt( (a-c)**2+4.0d0*b**2 )
          msqh2 = 0.5d0*(a+c-d)
          if (msqh2.lt.0.0d0) then
            unphys = -2
            return
          endif
          msqh1 = 0.5d0*(a+c+d)
          sin2al = 2.0d0*b/d
          cos2al = (a-c)/d
          if (sin2al.eq.0.0d0) then
            alpha = 0.0d0
          else
            alpha = sign(0.5d0*acos(cos2al),sin2al)
          endif
          mass(kh1) = sqrt(msqh1)
          mass(kh2) = sqrt(msqh2)
          mass(kh3) = sqrt(msqa)
          mass(khc) = sqrt(msqa+msqw)

        endif

c---------------- one-loop higgs sector (brignole-ellis-ridolfi-zwirner)
      else if (higloop.eq.1) then

        if (mass(kst(1)).eq.0.0)  call dssfesct(unphys)
        if (unphys.ne.0) return

        ! h1, h2 in ellis-ridolfi-zwirner plb257(91)83; plb262(91)477

        mt4 = mass(kt)**4
        mb4 = mass(kb)**4
        msqst1 = mass(kst(1))**2
        msqst2 = mass(kst(2))**2
        msqsb1 = mass(ksb(1))**2
        msqsb2 = mass(ksb(2))**2

        if (msqst1-msqst2.eq.0.d0) then
           rtop = 0.d0
           ltop = 0.d0
           gtop = 0.d0
        else
           rtop = (atop+mu/tanbe) / (msqst1-msqst2)
           ltop = log(msqst1/msqst2)
           gtop = 2.0d0 -
     &          (msqst1+msqst2)/(msqst1-msqst2) * ltop
        endif
        if (msqsb1-msqsb2.eq.0.d0) then
           rbot = 0.d0
           lbot = 0.d0
           gbot = 0.d0
        else
           rbot = (abot+mu*tanbe) / (msqsb1-msqsb2)
           lbot = log(msqsb1/msqsb2)
           gbot = 2.0d0 -
     &          (msqsb1+msqsb2)/(msqsb1-msqsb2) * lbot
        endif

        d11 = mb4/cb2 *
     &    ( log(msqsb1*msqsb2/mb4) +
     &      2.0d0 * abot * rbot * lbot ) +
     &    mb4/cb2 * abot**2 * rbot**2 * gbot +
     &    mt4/sb2 * mu**2 * rtop**2 * gtop
        d22 = mt4/sb2 *
     &    ( log(msqst1*msqst2/mt4) +
     &      2.0d0 * atop * rtop * ltop ) +
     &    mt4/sb2 * atop**2 * rtop**2 * gtop +
     &    mb4/cb2 * mu**2 * rbot**2 * gbot
        d12 = mt4/sb2 * mu * rtop *
     &    (ltop + atop * rtop * gtop) +
     &    mb4/cb2 * mu * rbot *
     &    (lbot + abot * rbot * gbot)

        ! d-corrections in brignole are missing

        coeff = 3.0d0*g2weak**2/(16.0d0*pi*pi*msqw)

        a = msqz*cb2 + msqa*sb2 + coeff*d11
        b = -sb*cb*(msqz+msqa) + coeff*d12
        c = msqz*sb2 + msqa*cb2 + coeff*d22
        if (prtlevel.ge.2) then
           write (*,'('' higgs mass^2 matrix = '',2(f12.2))')
     &          a,b
           write (*,'(''                       '',2(f12.2))')
     &          b,c
        endif

        d = sqrt( (a-c)**2+4.0d0*b**2 )
        msqh2 = 0.5d0*(a+c-d)
        msqh1 = 0.5d0*(a+c+d)
        sin2al = 2.0d0*b/d
        cos2al = (a-c)/d
        if (sin2al.eq.0.0d0) then
          alpha = 0.0d0
        else
          alpha = sign(0.5d0*acos(cos2al),sin2al)
        endif

        ! hc in brignole-ellis-ridolfi-zwirner plb271(91)123
        ! (their mu in eq 4 has opposite sign to mine -- april 3, 1996)
        ! (their mu in eq 5 has same sign as mine -- april 3, 1996)
        ! (their appendix has the same sign -- april 7, 1996)

        coeff = g2weak**2/(2*msqw)
        g = mb2*sb2 + mt2*cb2 - msqw*cb2*sb2
        b012 = coeff/(cb2*sb2) *
     &       ( abot*mb2*mu*sb2 - mb2*mt2*sb*cb +
     &         atop*mt2*mu*cb2 - mt2*mb2*sb*cb )
        bd12 = g2weak**2/sin2be * g
        c012 = coeff/(cb3*sb3) *
     &       ( mu**2*mb4*sb4 -
     &         (mq2+mu2+2*mt2)*mu*abot*mb2*sb3*cb +
     &         (mq2+mu2+2*mt2+mu**2-atop*(atop+abot))*
     &           mt2*mb2*sb2*cb2 +
     &         mu**2*mt4*cb4 -
     &         (mq2+md2+2*mb2)*mu*atop*mt2*cb3*sb +
     &         (mq2+md2+2*mb2+mu**2-abot*(abot+atop))*
     &           mb2*mt2*cb2*sb2 )
        cd12 = -g2weak**2/(2*cb2*sb2) *
     &       ( - mu*c2b/(2*costhw**2)*
     &           (atop*mt2*cb2-abot*mb2*sb2) +
     &         sinbe*cosbe*(mb2*sb2*(mu**2-abot**2)+
     &                      mt2*cb2*(mu**2-atop**2)) +
     &         sinbe*cosbe*(mu2+md2+mb2+mt2+msqw/3*tw2*c2b)*g)
        d012 = -coeff/(cb3*sb3) *
     &       ( mu**2*(mu2+mt2)*mb4*sb4 +
     &         (atop*mt2*(mu**2+2*mb2+atop*abot)-
     &           abot*(mu2+mt2)*(mq2+mt2))*mu*mb2*cb*sb3 +
     &         (md2*mu2/2+mb2*(mq2+mu2+mt2+2*mu**2)-2*mq2*atop*abot-
     &           abot**2*(mu2+mt2)+(mq2+mu**2+atop*abot)**2/2)*
     &           mt2*mb2*sb2*cb2 +
     &         mu**2*(md2+mb2)*mt4*cb4 +
     &         (abot*mb2*(mu**2+2*mt2+abot*atop)-
     &           atop*(md2+mb2)*(mq2+mb2))*mu*mt2*sb*cb3 +
     &         (mu2*md2/2+mt2*(mq2+md2+mb2+2*mu**2)-2*mq2*abot*atop-
     &           atop**2*(md2+mb2)+(mq2+mu**2+abot*atop)**2/2)*
     &           mb2*mt2*cb2*sb2 )
        dd12 = g2weak**2/(72*cb3*sb3) *
     &       ( msqw*tw2**2*sb2*cb2*c2b**2*(7*mt2*mb2-8*msqw*g) +
     &         6*tw2*cb2*sb2*c2b*(-2*msqw*g*(mu2+mt2-2*md2-2*mb2)+
     &           2*mt2*mb2*(mu2+mq2-2*md2)+3*mt2*mb2*(mt2-mb2)) +
     &         9*cb2*sb2*(mt2*mb2*msqw*c2b**2+2*c2b*mt2*mb2*(mt2-mb2)+
     &           4*(mu2+mt2)*(md2+mb2)*g) ) +
     &       g2weak**2/(6*cb3*sb3) *
     &       ( -tw2*cb2*sb2*c2b*(2*abot**2*mb2*(msqw*sb2-mt2)-
     &           atop**2*mt2*(msqw*cb2-mb2)+atop*abot*mt2*mb2) -
     &         3*sb2*cb2*(abot**2*mb2*sb2*(mu2+mt2)+
     &           atop**2*mt2*cb2*(md2+mb2)-abot*atop*mt2*mb2) ) -
     &       g2weak**2*mu/(36*cb2*sb2) *
     &       ( tw2**2*c2b**2*(2*abot*mb2*sb2-atop*mt2*cb2)*msqw -
     &         tw2*c2b*(abot*mb2*sb2*(6*msqw*c2b+9*mt2+12*mq2-3*mu2)+
     &           atop*mt2*cb2*(3*msqw*c2b-9*mb2-6*mq2-3*md2)) -
     &         9*c2b*abot*mb2*sb2*(mu2+mt2) +
     &         9*c2b*atop*mt2*cb2*(md2+mb2) -
     &         36*cb2*sb2*mt2*mb2*(abot+atop) ) -
     &       g2weak**2*mu**2/(12*cb3*sb3) *
     &       ( tw2*c2b*(msqw*cb2*sb2*(2*mt2*cb2-4*mb2*sb2) -
     &           2*mt4*cb4+4*mb4*sb4-2*cb2*sb2*mt2*mb2) -
     &         6*mb2*sb4*cb2*(mu2+2*mt2) -
     &         6*mt2*cb4*sb2*(md2+2*mb2) )
        b12 = b012 + bd12
        c12 = c012 + cd12
        d12 = d012 + dd12
        if ((msqst1-msqst2)*(msqst1-msqsb1)*(msqst1-msqsb2).eq.0.d0)
     &       then
           h1 = 0.d0
        else
           h1 = - (b12*msqst1**2+c12*msqst1+d12) /
     &          ((msqst1-msqst2)*(msqst1-msqsb1)*(msqst1-msqsb2)) +
     &          g2weak**2/(2*sb2)*(mt2/msqw)*(atop*mu/(msqst1-msqst2))
        endif
        if ((msqst2-msqst1)*(msqst2-msqsb1)*(msqst2-msqsb2).eq.0.d0)
     &       then
           h2 = 0.d0
        else
           h2 = - (b12*msqst2**2+c12*msqst2+d12) /
     &          ((msqst2-msqst1)*(msqst2-msqsb1)*(msqst2-msqsb2)) -
     &          g2weak**2/(2*sb2)*(mt2/msqw)*(atop*mu/(msqst1-msqst2))
        endif
        if ((msqsb1-msqst1)*(msqsb1-msqst2)*(msqsb1-msqsb2).eq.0.d0)
     &       then
           h3 = 0.d0
        else
           h3 = - (b12*msqsb1**2+c12*msqsb1+d12) /
     &          ((msqsb1-msqst1)*(msqsb1-msqst2)*(msqsb1-msqsb2)) +
     &          g2weak**2/(2*cb2)*(mb2/msqw)*(abot*mu/(msqsb1-msqsb2))
        endif
        if ((msqsb2-msqst1)*(msqsb2-msqst2)*(msqsb2-msqsb1).eq.0.d0)
     &       then
           h4 = 0.d0
        else
           h4 = - (b12*msqsb2**2+c12*msqsb2+d12) /
     &          ((msqsb2-msqst1)*(msqsb2-msqst2)*(msqsb2-msqsb1)) -
     &          g2weak**2/(2*cb2)*(mb2/msqw)*(abot*mu/(msqsb1-msqsb2))
        endif

        q2 = msqw
        ftop = 2*mt2*(log(mt2/q2)-1.0)
        fbot = 2*mb2*(log(mb2/q2)-1.0)
        deltasq =
     &       -3*g2weak**2/(64*pi**2) *
     &          (mt2/sb2+mb2/cb2-msqw) *
     &          ((ftop-fbot)/(mt2-mb2)) +
     &       3/(16*pi**2*sb*cb*(mt2-mb2)) *
     &          (h1*msqst1*(mt2*log(msqst1/mt2)-mb2*log(msqst1/mb2))+
     &           h2*msqst2*(mt2*log(msqst2/mt2)-mb2*log(msqst2/mb2))+
     &           h3*msqsb1*(mt2*log(msqsb1/mt2)-mb2*log(msqsb1/mb2))+
     &           h4*msqsb2*(mt2*log(msqsb2/mt2)-mb2*log(msqsb2/mb2)))
        msqhc = msqw+msqa+deltasq

        if (msqh2.le.0.0d0) then
          unphys = -3
          call dswrite(1,1,'negative h2 mass squared')
          return
        endif
        if (msqh1.le.0.0d0) then
          unphys = -4
          call dswrite(1,1,'negative h1 mass squared')
          return
        endif
        if (msqhc.le.0.0d0) then
          unphys = -5
          call dswrite(1,1,'negative h+ mass squared')
          return
        endif
        mass(kh1) = sqrt(msqh1)
        mass(kh2) = sqrt(msqh2)
        mass(kh3) = sqrt(msqa)
        mass(khc) = sqrt(msqhc)

c-------- one-loop higgs sector (drees-nojiri, prd45(92)2482 + my errata)
      else if (higloop.eq.2) then

        if (mass(kst(1)).eq.0.0)  call dssfesct(unphys)
        if (unphys.ne.0) return

        ! h1, h2

        mt2 = mass(kt)**2
        mb2 = mass(kb)**2
        mt4 = mass(kt)**4
        mb4 = mass(kb)**4
        msqst1 = mass(kst(1))**2
        msqst2 = mass(kst(2))**2
        msqsb1 = mass(ksb(1))**2
        msqsb2 = mass(ksb(2))**2
        ht = yukawa(kt)
        hb = yukawa(kb)
        ht2 = ht**2
        hb2 = hb**2

        if (msqst1-msqst2.eq.0.d0) then
           rtop = 0.d0
           ltop = 0.d0
           gtop = 0.d0
        else
           rtop = (atop+mu/tanbe) / (msqst1-msqst2)
           ltop = log(msqst1/msqst2)
           gtop = 2.0d0 -
     &          (msqst1+msqst2)/(msqst1-msqst2) * ltop
        endif
        if (msqsb1-msqsb2.eq.0.d0) then
           rbot = 0.d0
           lbot = 0.d0
           gbot = 0.d0
        else
           rbot = (abot+mu*tanbe) / (msqsb1-msqsb2)
           lbot = log(msqsb1/msqsb2)
           gbot = 2.0d0 -
     &          (msqsb1+msqsb2)/(msqsb1-msqsb2) * lbot
        endif

        d11 = 3.d0/(8.d0*pi**2) * (
     &    mb2 * hb2 * ( log(msqsb1*msqsb2/mb4) +
     &                  2.d0 * abot * rbot * lbot +
     &                  abot**2 * rbot**2 * gbot ) +
     &    mt2 * ht2 * mu**2 * rtop**2 * gtop )
        d22 = 3.d0/(8.d0*pi**2) * (
     &    mt2 * ht2 * ( log(msqst1*msqst2/mt4) +
     &                  2.d0 * atop * rtop * ltop +
     &                  atop**2 * rtop**2 * gtop ) +
     &    mb2 * hb2 * mu**2 * rbot**2 * gbot )
        d12 = 3.d0/(8.d0*pi**2) * mu * (
     &    mb2 * hb2 * rbot * ( lbot + abot * rbot * gbot ) +
     &    mt2 * ht2 * rtop * ( ltop + atop * rtop * gtop ) )

        a = msqz*cb2 + msqa*sb2 + d11
        b = -sb*cb*(msqz+msqa) + d12
        c = msqz*sb2 + msqa*cb2 + d22
        if (prtlevel.ge.2) then
           write (*,'('' higgs mass^2 matrix = '',2(f12.2))')
     &          a,b
           write (*,'(''                       '',2(f12.2))')
     &          b,c
        endif

        d = sqrt( (a-c)**2+4.0d0*b**2 )
        msqh2 = 0.5d0*(a+c-d)
        msqh1 = 0.5d0*(a+c+d)
        sin2al = 2.0d0*b/d
        cos2al = (a-c)/d
        if (sin2al.eq.0.0d0) then
          alpha = 0.0d0
        else
          alpha = sign(0.5d0*acos(cos2al),sin2al)
        endif

        ! hc

        ut(1,1) = squlmx(3,3) ! ut(1,l)
        ut(1,2) = squrmx(3,3) ! ut(1,r)
        ut(2,1) = squlmx(6,3) ! ut(2,l)
        ut(2,2) = squrmx(6,3) ! ut(2,r)
        ub(1,1) = sqdlmx(3,3) ! ub(1,l)
        ub(1,2) = sqdrmx(3,3) ! ub(1,r)
        ub(2,1) = sqdlmx(6,3) ! ub(2,l)
        ub(2,2) = sqdrmx(6,3) ! ub(2,r)
        x(1,1) = -hb2*v1+0.5d0*g2weak**2*v1 ! d-term added (april 4, 96)
        x(1,2) = -ht*mu
        x(2,1) = hb*abot ! sign changed my check (april 4, 96)
        x(2,2) = -ht*hb*v2 ! sign changed my check (april 4, 96)
        y(1,1) = -ht2*v2+0.5d0*g2weak**2*v2 ! d-term added (april 4, 96)
        y(1,2) = ht*atop
        y(2,1) = -hb*mu ! sign changed my check (april 4, 96)
        y(2,2) = -ht*hb*v1 ! sign changed my check (april 4, 96)

        do i=1,2
          do j=1,2
            cc(i,j) = 0.d0
            dd(i,j) = 0.d0
            do k=1,2
              do l=1,2
                cc(i,j) = cc(i,j) + ub(i,l) * x(l,k) * ut(j,k)
                dd(i,j) = dd(i,j) + ub(i,l) * y(l,k) * ut(j,k)
              enddo
            enddo
          enddo
        enddo

        q2 = msqw
        ftop = 2.d0*mt2*(log(mt2/q2)-1.0d0)
        fbot = 2.d0*mb2*(log(mb2/q2)-1.0d0)
        mst2(1) = msqst1
        mst2(2) = msqst2
        msb2(1) = msqsb1
        msb2(2) = msqsb2

        delta0 = 1.d0/(32.d0*pi**2) * ( ! 64->32 correction april 7, 1996
     &    3.d0 * ht2 * atop*mu * dsg0loop(q2,msqst1,msqst2) +
     &    3.d0 * hb2 * abot*mu * dsg0loop(q2,msqsb1,msqsb2) )

        deltasq = 2.d0*delta0/sin2be +
     &    3.d0/(16.d0*pi**2*sin2be) * (
     &      - 2.d0 * (ftop-fbot) * ht2*hb2*v1*v2/(mt2-mb2) )
        do i=1,2
          do j=1,2
            deltasq = deltasq + 3.d0/(16.d0*pi**2*sin2be) * (
     &        cc(i,j)*dd(i,j) * dsg0loop(q2,mst2(j),msb2(i)) )
          enddo
        enddo

        msqhc = msqw+msqa+deltasq

        if (msqh2.le.0.0d0) then
          unphys = -3
          call dswrite(1,1,'negative h2 mass squared')
          return
        endif
        if (msqh1.le.0.0d0) then
          unphys = -4
          call dswrite(1,1,'negative h1 mass squared')
          return
        endif
        if (msqhc.le.0.0d0) then
          unphys = -5
          call dswrite(1,1,'negative h+ mass squared')
          return
        endif
        mass(kh1) = sqrt(msqh1)
        mass(kh2) = sqrt(msqh2)
        mass(kh3) = sqrt(msqa)
        mass(khc) = sqrt(msqhc)

c------------------------------ rge improved 1-loop effective potential
c-------------------- carena-espinosa-quiros-wagner, plb 355 (1995) 209
      else if (higloop.eq.3) then

        if (mass(kst(1)).eq.0.0d0)  call dssfesct(unphys)
        if (unphys.ne.0) return
        au=asoftu(3)  ! je should be + here, convention of 960525
        ad=asoftd(3)  ! je should be + here, convention of 960525
        alpha2=alphem/s2thw
        alpha1=alpha2*s2thw/(1.0d0-s2thw)
        alpha3z=alph3mz ! alph3->alph3mz 020912 (JE)
        pi2 = pi**2
        msqa = ma**2
        rmbot = 3.d0            ! mbottom(mtop) = 3 gev
        alpha3 = alpha3z/(1.d0 +(11.d0 - 10.d0/3.d0)/4.d0/pi*alpha3z*
     &       log(mt2/msqz))
        rmtop = mass(kt)/(1.d0+4.d0*alpha3/3.d0/pi)
        ms2 = 0.5d0*(mass(kst(1))**2 + mass(kst(2))**2) + mt2
        ms4 = ms2**2
        t = log(ms2/mt2)

        g1sq = alpha1*4.d0*pi
        g2sq = alpha2*4.d0*pi
        g3sq = alpha3*4.d0*pi
        hu2 = (rmtop/v/sinbe)**2
        hd2 = (rmbot/v/cosbe)**2
        hu4 = hu2**2
        hd4 = hd2**2
        musq = mu**2

        xau = (2.d0*au**2/ms2)*(1.d0 - au**2/12.d0/ms2)
        xad = (2.d0*ad**2/ms2)*(1.d0 - ad**2/12.d0/ms2)
        aud = (-6.d0*musq/ms2 - ( musq- ad*au)**2/ms4
     &       + 3.d0*(au + ad)**2/ms2)/6.d0

        lam1 = ((g1sq+g2sq)/4.d0)*(1.d0-3.d0*hd2*t/8.d0/pi2)
     &       +(3.d0*hd4/8.d0/pi2)*(t+xad/2.d0+(3.d0*hd2/2.d0+hu2/2.d0
     &       -8.d0*g3sq)*(xad*t+t**2)/16.d0/pi2)
     &       -(3.d0*hu4*musq**2/96.d0/pi2/ms4)*(1.d0+(9.d0*hu2-5.d0*hd2
     &       -16.d0*g3sq)*t/16.d0/pi2)
        lam2 = ((g1sq+g2sq)/4.d0)*(1.d0-3.d0*hu2*t/8.d0/pi2)
     &       +(3.d0*hu4/8.d0/pi2)*(t+xau/2.d0+(3.d0*hu2/2.d0+hd2/2.d0
     &       -8.d0*g3sq)*(xau*t+t**2)/16.d0/pi2)
     &       -(3.d0*hd4*musq**2/96.d0/pi2/ms4)*(1.d0+(9.d0*hd2-5.d0*hu2
     &       -16.d0*g3sq)*t/16.d0/pi2)
        lam3 = ((g2sq-g1sq)/4.d0)*(1.d0-3.d0*(hu2+hd2)*t/16.d0/pi2)
     &       +(6.d0*hu2*hd2/16.d0/pi2)*(t+aud/2.d0+(hu2+hd2
     &       -8.d0*g3sq)*(aud*t+t**2)/16.d0/pi2)
     &       +(3.d0*hu4/96.d0/pi2)*(3.d0*musq/ms2-musq*au**2/ms4)*
     &       (1.d0+(6.d0*hu2-2.d0*hd2/2.d0-16.d0*g3sq)*t/16.d0/pi2)
     &       +(3.d0*hd4/96.d0/pi2)*(3.d0*musq/ms2-musq*ad**2/
     &       ms4)*(1.d0+(6.d0*hd2-2.d0*hu2-16.d0*g3sq)*t/16.d0/pi2)
        lam4 = (-g2sq/2.d0)*(1.d0-3.d0*(hu2+hd2)*t/16.d0/pi2)
     &       -(6.d0*hu2*hd2/16.d0/pi2)*(t+aud/2.d0+(hu2+hd2
     &       - 8.d0*g3sq)*(aud*t+t**2)/16.d0/pi2)
     &       +(3.d0*hu4/96.d0/pi2)*(3.d0*musq/ms2-musq*au**2/ms4)*
     &       (1.d0+(6.d0*hu2-2.d0*hd2-16.d0*g3sq)*t/16.d0/pi2)
     &       +(3.d0*hd4/96.d0/pi2)*(3.d0*musq/ms2-musq*ad**2/ms4)*
     &       (1.d0+(6.d0*hd2-2.d0*hu2/2.d0-16.d0*g3sq)*t/16.d0/pi2)
        lam5 = -(3.d0*hu4*musq*au**2/96.d0/pi2/ms4)*
     &       (1.d0-(2.d0*hd2-6.d0*hu2+16.d0*g3sq)*t/16.d0/pi2)
     &       -(3.d0*hd4*musq*ad**2/96.d0/pi2/ms4)*
     &       (1.d0-(2.d0*hu2-6.d0*hd2+16.d0*g3sq)*t/16.d0/pi2)
        lam6 = mu*(3.d0*hu4*musq*au/96.d0/pi2/ms4)*(1.d0-
     &       (7.d0*hd2/2.d0-15.d0*hu2/2.d0+16.d0*g3sq)*t/16.d0/pi2)
     &       +(3.d0*hd4*mu*(ad**3/ms2-6.d0*ad)/ms2/96.d0/pi2)*
     &       (1.d0-(hu2/2.d0-9.d0*hd2/2.+16.d0*g3sq)*t/16.d0/pi2)
        lam7 = mu*(3.d0*hd4*musq*ad/96.d0/pi2/ms4)*(1.d0-
     &       (7.d0*hu2/2.d0-15.d0*hd2/2.d0+16.d0*g3sq)*t/16.d0/pi2)
     &       +(3.d0*hu4*mu*au*(au**2/ms2-6.d0)/96.d0/pi2/ms2)*
     &       (1.d0-(hd2/2.d0-9.d0*hu2/2.d0+16.d0*g3sq)*t/16.d0/pi2)

        msqhc = msqa+(lam5-lam4)*vsq

        a = 2.d0*vsq*(lam1*cb2+2.d0*lam6*sb*cb+lam5*sb2)+
     &       msqa*sb2
        b = 2.d0*vsq*(sb*cb*(lam3+lam4)+lam6*cb2+lam7*sb2)-
     &       msqa*sb*cb
        c = 2.d0*vsq*(lam2*sb2+2.d0*lam7*sb*cb+lam5*cb2)+
     &       msqa*cb2
        if (prtlevel.ge.2) then
           write (*,'('' higgs mass^2 matrix = '',2(f12.2))')
     &          a,b
           write (*,'(''                       '',2(f12.2))')
     &          b,c
        endif

        d = dsqrt( (a-c)**2+4.0d0*b**2 )
        msqh2 = 0.5d0*(a+c-d)
        msqh1 = 0.5d0*(a+c+d)
        sin2al = 2.0d0*b/d
        cos2al = (a-c)/d
        if (sin2al.eq.0.0d0) then
           alpha = 0.0d0
        else
           alpha = sign(0.5d0*acos(cos2al),sin2al)
        endif

        if (msqh2.le.0.0d0) then
           unphys = -3
           call dswrite(1,1,'negative h2 mass squared')
           return
        endif
        if (msqh1.le.0.0d0) then
           unphys = -4
           call dswrite(1,1,'negative h1 mass squared')
           return
        endif
        if (msqhc.le.0.0d0) then
           unphys = -5
           call dswrite(1,1,'negative h+ mass squared')
           return
        endif
        mass(kh1) = sqrt(msqh1)
        mass(kh2) = sqrt(msqh2)
        mass(kh3) = sqrt(msqa)
        mass(khc) = sqrt(msqhc)

        msqst1 = mass(kst(1))**2
        msqst2 = mass(kst(2))**2
        msqsb1 = mass(ksb(1))**2
        msqsb2 = mass(ksb(2))**2
        msusy2 = 0.5d0*(msqst1+msqst2)

        if (abs(msqst2-msqst1).gt.0.5*(msqst2+msqst1)) then
           call dswrite(1,1,
     &          'carena et al expressions are out of range of validity')
           call dswrite(1,1,'|mst2^2-mst1^2|/(mst1^2+mst2^2)>0.5')
           hwarning=ibset(hwarning,0)
        endif
        if (abs(msqsb2-msqsb1).gt.0.5*(msqsb2+msqsb1)) then
           call dswrite(1,1,
     &          'carena et al expressions are out of range of validity')
           call dswrite(1,1,'|msb2^2-msb1^2|/(msb1^2+msb2^2)>0.5')
           hwarning=ibset(hwarning,1)
        endif

c------------------------------ rge improved 1-loop effective potential
c-------------------- carena-quiros-wagner, nucl. phys. b461 (1996) 407
      else if (higloop.eq.4) then
        tanb=tanbe
        mq=sqrt(mass2q(3))
        mur=sqrt(mass2u(3))
        mdr=sqrt(mass2d(3))
        mtop=mass(kt)
        au=asoftu(3)  ! je should be + here, convention of 960525
        ad=asoftd(3)  ! je should be + here, convention of 960525
        mchi=max(mass(kcha(1)),mass(kcha(2)))
        v = dsqrt(vsq)
        mz = mass(kz)
        alpha2=alphem/s2thw
        alpha1=alpha2*s2thw/(1.0d0-s2thw)
        alpha3z=alph3mz   ! alph3->alph3mz 020912 (JE)
 	call dspole(mchi,ma,tanb,mq,mur,mdr,mtop,au,ad,mu,mh,mhp,hm,hmp,
     &       amp,sa,ca,v,mz,alpha1,alpha2,alpha3z,s2thw,lambda,
     &       prtlevel,ierr)
        if (ierr.ne.0) then
           unphys = -6
           call dswrite(1,1,'error in dspole')
           return
        endif
        mass(kh1)=hmp
        mass(kh2)=mhp
        mass(kh3)=amp
        alpha=asin(sa)
        lam1=lambda(1)
        lam2=lambda(2)
        lam3=lambda(3)
        lam4=lambda(4)
        lam5=lambda(5)
        lam6=lambda(6)
        lam7=lambda(7)

        msqhc = msqa+(lam5-lam4)*vsq
        if (msqhc.le.0.0d0) then
           unphys = -5
           call dswrite(1,1,'negative h+ mass squared')
           return
        endif
        mass(khc) = sqrt(msqhc)

c-------------------------------- FeynHiggs: Heinemeyer, Hollik, Weiglein

c...This is for the full FeynHiggs package
      else if (higloop.eq.5) then
 
c...Below is for FeynHiggs 1.2.2       
c        mtop=Mass(kt)
c        Au=Asoftu(3)
c        Ad=Asoftd(3)
c        Mqtl=dsqrt(mass2q(3))
c        Mqtr=dsqrt(mass2u(3))
c        Mqbl=dsqrt(mass2q(3))
c        Mqbr=dsqrt(mass2d(3))
c        msbarselec = 1 ! 1 if Mqtl, ... are on-shell, 2 if they are MSbar
c
c        call dsfeynhiggs(hwarning,HM,mh,mhch,halpha,delrho,
c     &    Asoftu(3),Asoftd(3),mass(kt),mass(kb),
c     &    mu,M2,Mqtl,Mqtr,Mqbl,Mqbr,mass(kgluin),ma,tanbe,
c     &    msbarselec)
c
c        Mass(kh1)=HM
c        Mass(kh2)=mh
c        Mass(kh3)=ma
c        Mass(khc)=mhch
c        Alpha=halpha

         call dsfeynhiggs(hwarning) ! FeynHiggs

c----------------------------FeynHiggsFast: Heinemeyer, Hollik, Weiglein
c... This implementation should be redone to send our masses to
c... FeynHiggsFast instead of letting FeynHiggsFast calculate them.

      else if (higloop.eq.6) then
        
        mtop=Mass(kt)
        Au=Asoftu(3)
        Ad=Asoftd(3)
        Mqtl=dsqrt(mass2q(3))
        Mqtr=dsqrt(mass2u(3))
        Mqbl=dsqrt(mass2q(3))
        Mqbr=dsqrt(mass2d(3))

        call dsfeynhiggsfast(hwarning,HM,mh,mhch,halpha,delrho,
     &    Asoftu(3),Asoftd(3),mass(kt),mass(kb),
     &    mu,M2,Mqtl,Mqtr,Mqbl,Mqbr,mass(kgluin),ma,tanbe)

        Mass(kh1)=HM
        Mass(kh2)=mh
        Mass(kh3)=ma
        Mass(khc)=mhch
        Alpha=halpha

c---------------------------------------------- unrecognized higloop
      else
        call dswrite(0,1,'unrecognized option higloop='//dsi2s(higloop))
        stop
      endif

      if (prtlevel.ge.2) then
        write (message,'('' lambda(1..7) ='',7(1x,f10.7))')
     &        lam1,lam2,lam3,lam4,lam5,lam6,lam7
        call dswrite(2,0,message)
      endif

      end











