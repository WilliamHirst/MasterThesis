      subroutine dsacbnd1(excl)
c_______________________________________________________________________
c  check if accelerator data exclude the present point.
c  output:
c    excl - code of the reason for exclusion (integer); 0 if allowed
c  author: paolo gondolo 1994-1999
c  history:
c     940407 first version paolo gondolo
c     950323 update paolo gondolo
c     971200 partial update joakim edsjo
c     980428 update joakim edsjo
c     990719 update paolo gondolo
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
      excl=0

c----------------------------------------------bounds from b -> s gamma
      call dsbsgamma(bsg,bsgqcd)
      if (bsg.lt.1.00d-4.or.bsg.gt.4.0d-4) then
         excl=ibset(excl,7)
      endif

c------------------------------------------lower bound on chargino mass
c...je update 98-04-28 from aleph, talk by j. carr, march 31, 1998
c...http://alephwww.cern.ch/alpub/seminar/carrlepc98/index.html
      if (abs(mass(kcha(2))).lt.91.0d0.and.
     &  abs(mass(kn(1))-mass(kcha(2))).gt.4.0d0)
     &  excl=ibset(excl,0)
c        pdg 99 acciarri 96 l3 plb377,289
      if (abs(mass(kcha(2))).lt.64.d0 .and.
     &     abs(mass(kn(kln))).lt.43.0d0 .and.
     &     abs(mass(kcha(2))).lt.abs(mass(kn(2))))
     &     excl=ibset(excl,0)
c        pdg 94 p 1801 iii decamp92 aleph
      if (abs(mass(kcha(2))).lt.47.d0 .and.
     &     abs(mass(kn(kln))).lt.41.0d0)
     &     excl=ibset(excl,0)
c        pdg 94 p 1801 iv hidaka91 cdf
      if (abs(mass(kcha(1))).lt.99.0d0) excl=ibset(excl,0)

c--------------------------------------------lower bound on gluino mass
c        pdg 94 p 1803 i abe92 cdf was 218 (90%cl)
c        pdg 99 abachi95 d0 prl75,618 is 212 (95%cl)
      excluded=.false.
      do i=1,6
        excluded = excluded .or. (mass(ksqu(i)).lt.mass(kgluin)) .or.
     &        (mass(ksqd(i)).lt.mass(kgluin))
      enddo
      if (excluded .and. mass(kgluin).lt.212.0d0) excl=ibset(excl,1)
c        pdg 94 p 1803 ii abe92 cdf was 100
c        pdg 99 abe97 cdf prd56, r1357 is 173-10 (see note)
      if (mass(kgluin).lt.163.0d0) excl=ibset(excl,1)

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
      if (mass(ksnu(1)).lt.44.4)
     &     excl=ibset(excl,3)
      if (mass(ksnu(2)).lt.44.4)
     &     excl=ibset(excl,3)
      if (mass(ksnu(3)).lt.44.4)
     &     excl=ibset(excl,3)

c-------------------------------------lower bounds on neutralino masses
c     assumes neutralinos are ordered by increasing mass
c        pdg 94 p 1799 i decamp92 phys rep 216, 253 aleph
c        pdg 99 acciarri95 l3 plb350,109
      if (abs(mass(kn(1))).lt.23.0d0 .and.
     &     tanbe.gt.3.d0) excl=ibset(excl,6)
      if (abs(mass(kn(1))).lt.20.0d0 .and.
     &     tanbe.gt.2.d0) excl=ibset(excl,6)
c        pdg 99 buskulic96 zphys c72,549 aleph
      if ( (abs(mass(kn(1))).lt.12.8.and.mass(ksnu(1)).gt.200.d0).or.
     &     (abs(mass(kn(1))).lt.12.8.and.mass(ksnu(2)).gt.200.d0).or.
     &     (abs(mass(kn(1))).lt.12.8.and.mass(ksnu(3)).gt.200.d0) )
     &     excl=ibset(excl,6)
c        pdg 99 acciarri98 l3 epj c4,207
      if (abs(mass(kn(1))).lt.10.9d0) excl=ibset(excl,6)
c        pdg 99 abbiendi99 opal epj c8,255
      if (abs(mass(kn(2))).lt.44.0d0) excl=ibset(excl,6)
c        pdg 99 abbiendi99 opal epj c8,255
      if (abs(mass(kn(3))).lt.102.d0) excl=ibset(excl,6)
c        pdg 99 acciarri95 l3 plb350,109
      if (abs(mass(kn(4))).lt.127.0d0) excl=ibset(excl,6)

c------------------------------------------------bounds on higgs masses

c      pdg 99 abbiendi99 opal epj c7,407
      if (mass(khc).lt.59.5d0) excl=ibset(excl,5)

c...gao, gay aleph 99-053 (conf 99-029) = hep99 tampere
c...http://alephwww.cern.ch/alpub
      beta=atan(tanbe)
      y=(sin(beta-alpha))**2
      mh2lim=82.5d0+(93.d0-82.5d0)*y
      if (mass(kh2).lt.mh2lim) excl=ibset(excl,5)

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

      end

