      subroutine dsacbnd2(excl)
c_______________________________________________________________________
c  check if accelerator data exclude the present point.
c  output:
c    excl - code of the reason for exclusion (integer); 0 if allowed
c    if not allowed, the reasons are coded as follows
c      bit set   dec.   oct.   reason
c      -------   ----   ----   ------
c            0      1      1   chargino mass
c            1      2      2   gluino mass
c            2      4      4   squark mass
c            3      8     10   slepton mass
c            4     16     20   invisible z width
c            5     32     40   higgs mass
c            6     64    100   neutralino mass
c            7    128    200   b -> s gamma
c            8    256    400   rho parameter
c  author: paolo gondolo 1994-1999
c  history:
c     940407 first version paolo gondolo
c     950323 update paolo gondolo
c     971200 partial update joakim edsjo
c     980428 update joakim edsjo
c     990719 update paolo gondolo
c     000310 update piero ullio
c     000424 added delrho joakim edsjo
c=======================================================================
      implicit none
      include 'dsmssm.h'
      include 'dsio.h'
      include 'dsaccom.h'
      integer i,j,excl
      logical excluded
      real*8 temp,mi,mj,ei,ej,p2,mz,mz2,dsabsq,
     &  beta,y,mh2lim,pi
      parameter (pi=3.141592653589793238d0)
      complex*16 gzij
      real*8 diffneucha,absmacha,absmachalim,arrmh2(25),
     & arrsin2ba(25),arrmh2cd(13),arrtanbe(13),lgtgbe,ylim
      integer iiikkk

      data (arrmh2(i),i=1,25)/
     & 0.914500d+02,0.924500d+02,0.936000d+02,0.944500d+02,
     & 0.955000d+02,0.960000d+02,0.979000d+02,0.990500d+02,
     & 0.100750d+03,0.101050d+03,0.102050d+03,0.102800d+03,
     & 0.103300d+03,0.103900d+02,0.104350d+03,0.105100d+03,
     & 0.105800d+03,0.106250d+03,0.106700d+03,0.106850d+03,
     & 0.107000d+03,0.107100d+03,0.107350d+03,0.107450d+03,
     & 0.107700d+03/

      data (arrsin2ba(i),i=1,25)/
     & 0.000d0,0.050d0,0.150d0,0.200d0,0.213d0,0.200d0,0.187d0,
     & 0.187d0,0.200d0,0.213d0,0.250d0,0.300d0,0.350d0,0.400d0,
     & 0.450d0,0.500d0,0.550d0,0.600d0,0.650d0,0.700d0,0.750d0,
     & 0.800d0,0.850d0,0.900d0,1.000d0/

      data (arrtanbe(i),i=1,13)/46.d0,48.d0,49.d0,50.d0,52.d0,53.d0
     & ,55.d0,55.8d0,56.d0,57.d0,58.d0,59.d0,60.d0/

      data (arrmh2cd(i),i=1,13)/90.d0,92.d0,94.d0,95.d0,96.d0,98.d0
     & ,100.d0,102.d0,104.d0,106.d0,108.d0,110.d0,112.d0/

      excl=0

c----------------------------------------------bounds from b -> s gamma
      call dsbsgamma(bsg,bsgqcd)
      if (bsg.lt.1.00d-4.or.bsg.gt.4.0d-4) then
         excl=ibset(excl,7)
      endif

c------------------------------------------lower bound on chargino mass
c...pu update 99-12-22 from l3 talk at lepc meeting in november 1999
c...http://l3www.cern.ch/analysis/latestresults.html
      diffneucha=abs(mass(kn(1))-mass(kcha(2)))
      absmacha=abs(mass(kcha(2)))
      if (absmacha.lt.85.9d0) excl=ibset(excl,0)
      if (absmacha.lt.94.2d0.and.diffneucha.le.0.2d-1)
     &  excl=ibset(excl,0)
      if (absmacha.lt.94.2d0.and.diffneucha.gt.10.d0)
     &  excl=ibset(excl,0)
      if (absmacha.lt.94.2d0.and.
     &  diffneucha.gt.4.2d0.and.diffneucha.le.10.d0) then
        absmachalim=92.3d0+(94.2d0-92.3d0)*(diffneucha-4.2d0)
     &   /(10.d0-4.2d0)
        if(absmacha.lt.absmachalim) excl=ibset(excl,0)
      endif
      if (absmacha.lt.92.3d0.and.
     &  diffneucha.gt.3.5d0.and.diffneucha.le.4.2d0) then
        absmachalim=85.9d0+(92.3d0-85.9d0)*(diffneucha-3.5d0)
     &   /(4.2d0-3.5d0)
        if(absmacha.lt.absmachalim) excl=ibset(excl,0)
      endif
c        pdg 94 p 1801 iv hidaka91 cdf
      if (abs(mass(kcha(1))).lt.99.0d0) excl=ibset(excl,0)

c--------------------------------------------lower bound on gluino mass
c...pu update 00-03-10   pdg 2000 edition update  http://pdg.lbl.gov
      if (mass(kgluin).lt.190.0d0) excl=ibset(excl,1)
      excluded=.false.
      do i=1,6
        excluded = excluded .or. (mass(ksqu(i)).lt.mass(kgluin)) .or.
     &        (mass(ksqd(i)).lt.mass(kgluin))
      enddo
      if (excluded .and. mass(kgluin).lt.260.0d0) excl=ibset(excl,1)

c------------------------------------------lower bound on squark masses
c        pdg 94 p 1802 i abe92 cdf is 90 with mgluino<410
      excluded=.false.
      do i=1,6
        excluded = excluded .or. (mass(ksqu(i)).lt.90.0d0) .or.
     &        (mass(ksqd(i)).lt.90.0d0)
      enddo
      if (excluded .and. mass(kgluin).lt.410.0d0) excl=ibset(excl,2)
c        pdg 94 abachi95 d0 prl75,618 is 176 with mgluino<300
      excluded=.false.
      do i=1,6
        excluded = excluded .or. (mass(ksqu(i)).lt.176.0d0) .or.
     &        (mass(ksqd(i)).lt.176.0d0)
      enddo
      if (excluded .and. mass(kgluin).lt.300.0d0) excl=ibset(excl,2)
c        pdg 94 p 1802 iii abe92 cdf was 180
c        pdg 99 abe96 cdf prl76,2006 is 224
      excluded=.false.
      do i=1,6
        excluded = excluded .or.
     &    (mass(ksqu(i)).lt.224.d0 .and. mass(ksqu(i)).gt.mass(kgluin))
     &    .or.
     &    (mass(ksqu(i)).lt.224.d0 .and. mass(ksqd(i)).gt.mass(kgluin))
      enddo
      if (excluded) excl=ibset(excl,2)

c---------------------------------lower bound on charged slepton masses
c        selectron pdg 94 p 1801 ii decamp92 aleph was 45 with 41
c                  pdg 99 barate98 aleph plb433,176 is 78 with 73
      if (mass(ksl(1)).lt.78.0d0 .and. abs(mass(kn(kln))).lt.73.d0)
     &     excl=ibset(excl,3)
      if (mass(ksl(2)).lt.78.0d0 .and. abs(mass(kn(kln))).lt.73.d0)
     &     excl=ibset(excl,3)
c        smuon pdg 94 p 1802 ii decamp92 aleph was 45 with 41
c              pdg 99 barate98 aleph plb433,176 is 71 with 66
      if (mass(ksl(3)).lt.71.0d0 .and. abs(mass(kn(kln))).lt.66.d0)
     &     excl=ibset(excl,3)
      if (mass(ksl(4)).lt.71.0d0 .and. abs(mass(kn(kln))).lt.66.d0)
     &     excl=ibset(excl,3)
c        stau pdg 94 p 1802 ii decamp92 aleph was 45 with 38
c             pdg 99 barate98 aleph plb433,176 is 65 with 55
      if (mass(ksl(5)).lt.65.0d0 .and. abs(mass(kn(kln))).lt.55.d0)
     &     excl=ibset(excl,3)
      if (mass(ksl(6)).lt.65.0d0 .and. abs(mass(kn(kln))).lt.55.d0)
     &     excl=ibset(excl,3)

c---------------------------------------lower bound on sneutrino masses
c        pdg 94 p 1801 ii adriani93 l3 was 37.1
c        pdg 99 their very own limit 44.4
c...pu update 00-03-10   pdg 2000 edition update  http://pdg.lbl.gov
      if (mass(ksnu(1)).lt.43.d0) excl=ibset(excl,3)
      if (mass(ksnu(2)).lt.43.d0) excl=ibset(excl,3)
      if (mass(ksnu(3)).lt.43.d0) excl=ibset(excl,3)

c-------------------------------------lower bounds on neutralino masses
c     assumes neutralinos are ordered by increasing mass
c...pu update 00-03-10 limit implemented from cern-ep/99-123
c...g. abbiendi et al. (opal collab.) europ. phys. j. c (2000) accepted
c...http://www.cern.ch/opal/pubs/paper/html/pr291.html
c do we want to implement the limits as a function of tan(beta)???
      if (abs(mass(kn(1))).lt.31.6d0) excl=ibset(excl,6)
      if (abs(mass(kn(2))).lt.55.9d0.and.
     &    abs(abs(mass(kn(2)))-abs(mass(kn(1)))).gt.10.d0)
     &    excl=ibset(excl,6)
      if (abs(mass(kn(3))).lt.106.0d0.and.
     &    abs(abs(mass(kn(2)))-abs(mass(kn(1)))).gt.10.d0)
     &    excl=ibset(excl,6)


c        pdg 94 p 1799 i decamp92 phys rep 216, 253 aleph
c        pdg 99 acciarri95 l3 plb350,109
c      if (abs(mass(kn(1))).lt.23.0d0 .and.
c     &     tanbe.gt.3.d0) excl=ibset(excl,6)
c      if (abs(mass(kn(1))).lt.20.0d0 .and.
c     &     tanbe.gt.2.d0) excl=ibset(excl,6)
c        pdg 99 buskulic96 zphys c72,549 aleph
c      if ( (abs(mass(kn(1))).lt.12.8.and.mass(ksnu(1)).gt.200.d0).or.
c     &     (abs(mass(kn(1))).lt.12.8.and.mass(ksnu(2)).gt.200.d0).or.
c     &     (abs(mass(kn(1))).lt.12.8.and.mass(ksnu(3)).gt.200.d0) )
c     &     excl=ibset(excl,6)
c        pdg 99 acciarri98 l3 epj c4,207
c      if (abs(mass(kn(1))).lt.10.9d0) excl=ibset(excl,6)

c        pdg 99 abbiendi99 opal epj c8,255
      if (abs(mass(kn(2))).lt.44.0d0) excl=ibset(excl,6)
c        pdg 99 abbiendi99 opal epj c8,255
      if (abs(mass(kn(3))).lt.102.d0) excl=ibset(excl,6)
c        pdg 99 acciarri95 l3 plb350,109
      if (abs(mass(kn(4))).lt.127.0d0) excl=ibset(excl,6)

c------------------------------------------------bounds on higgs masses
c h3:
c     approx lep  from phys. rep. 216(1992)253, p. 296
c     je correction 97-02-21
      lgtgbe = log10(tanbe)
      if (mass(kh3).lt.45.and.
     &     lgtgbe.gt.(5.5/29.+sqrt(15./1038.))) then
         excl=ibset(excl,5)
      else if (mass(kh3).lt.(30.0+lgtgbe*29.*5.).and.
     &     lgtgbe.lt.(5.5/29.-sqrt(25./1038.))) then
         excl=ibset(excl,5)
      else if (mass(kh3).lt.(60.-1038.*(lgtgbe-5.5/29.)**2).and.
     &     lgtgbe.le.(5.5/29.+sqrt(15./1038.)).and.
     &     lgtgbe.ge.(5.5/29.-sqrt(25./1038.))) then
         excl=ibset(excl,5)
      endif
c        h3 pdg 94 p 1370 iii adriani92 l3
      if (mass(kh3).lt.22.d0 .and. tanbe.gt.1.0.and.tanbe.lt.50.0)
     &   excl=ibset(excl,5)
c...pu update 00-03-10 limit implemented from opal document pn426
c...http://www.cern.ch/opal/pubs/latest_docs.html
      if(mass(kh3).lt.80.1d0) excl=ibset(excl,5)


c h2:
c...pu update 00-03-10 limit implemented from aleph 2000-006
c...http://alephwww.cern.ch/alpub/oldconf/oldconf_00.html
      beta=atan(tanbe)
      y=sin(beta-alpha)**2
      if(mass(kh2).lt.arrmh2(1)) then
        excl=ibset(excl,5)
        goto 440
      elseif(mass(kh2).gt.arrmh2(25)) then
        goto 440
      else
        ylim=0.0d0
        do iiikkk=1,24
          if(mass(kh2).lt.arrmh2(iiikkk+1).and
     &        .mass(kh2).ge.arrmh2(iiikkk)) then
            ylim=arrsin2ba(iiikkk)+(arrsin2ba(iiikkk+1)
     &        -arrsin2ba(iiikkk))*(mass(kh2)-arrmh2(iiikkk))
     &        /(arrmh2(iiikkk+1)-arrmh2(iiikkk))
            goto 444
          endif
        enddo
 444    if(y.lt.ylim) excl=ibset(excl,5)
      endif
 440  continue


c...for large tan(beta) check against the cdf limit as well
c...fermilab-conf-99/263-e presented at the hep conference in
c...in tampere, july 99
c...approximate limit in the plane tan(beta) mh2 in case of maximal
c...mixing, nearly coincident (for tan(beta)<60) with no-mixing case
      mh2lim=0.0d0
      if(tanbe.ge.arrtanbe(1).and.tanbe.lt.arrtanbe(13)) then
        do iiikkk=1,12
        if(tanbe.lt.arrtanbe(iiikkk+1)
     &         .and.tanbe.ge.arrtanbe(iiikkk)) then
          mh2lim=arrmh2cd(iiikkk)+(arrmh2cd(iiikkk+1)
     &            -arrmh2cd(iiikkk))*(tanbe-arrtanbe(iiikkk))
     &     /(arrtanbe(iiikkk+1)-arrtanbe(iiikkk))
          goto 555
        endif
      enddo
 555  if (mass(kh2).lt.mh2lim) excl=ibset(excl,5)
      endif


c hc:
c...pu update 00-03-10. limit implemented from aleph 2000-011
c...http://alephwww.cern.ch/alpub/oldconf/oldconf_00.html
      if (mass(khc).lt.77.7d0) excl=ibset(excl,5)

c-----------------------------------------bound from gamma_z(invisible)
      mz = mass(kz)
      mz2 = mz*mz
      ! 3 neutrinos
      gzinv = (g2weak/2/costhw)**2 * mz/(24.d0*pi)
      ! 3 sneutrinos
      do i=1,3
         temp = 1.d0-4.d0*mass(ksnu(i))**2/mz2
         if (temp.gt.0.d0) then
            temp = temp**(1.5d0) *
     &           (g2weak/2/costhw)**2 * mz/(48.d0*pi)
            gzinv = gzinv + temp
         endif
      enddo
      ! neutralinos
c...pg fixed june 30, 2000
      do i=1,4
         mi = abs(mass(kn(i)))
         do j=i,4
            mj = abs(mass(kn(j)))
            p2 = (mz2-(mi+mj)**2)*(mz2-(mi-mj)**2)/(4*mz2)
            ei = (mz2-mj**2+mi**2)/(2*mz)
            ej = (mz2-mi**2+mj**2)/(2*mz)
            if (p2.gt.0.d0.and.ei.gt.0.0d0.and.ej.gt.0.0d0) then
               gzij = g2weak/(2.d0*costhw)*
     &              (conjg(neunmx(i,3))*neunmx(j,3)-
     &               conjg(neunmx(i,4))*neunmx(j,4))
               temp=sqrt(p2)/(2.d0*pi*mz2) * (
     &              dsabsq(gzij) * (ei*ej+p2/3.0d0)-
     &              dreal(gzij**2)*mi*mj)
               if (i.eq.j) temp=0.5d0*temp
               gzinv = gzinv + temp
            endif
         enddo
      enddo
c        pdg 94 p 1358 (498.2+4.2)
      if (gzinv.gt.0.5024d0) excl=ibset(excl,4)

c   reason for exclusion
 1000 continue
      if (prtlevel.ge.1) call dswexcl(6,excl)
      return

c--------------------------------------------limit on the rho parameter
c...je addition 2000-04-24
c...uses delta rho as calculated by feynhiggs
      if (delrho.gt.1.3d-3) then
        excl=ibset(excl,8)
      endif


      end

