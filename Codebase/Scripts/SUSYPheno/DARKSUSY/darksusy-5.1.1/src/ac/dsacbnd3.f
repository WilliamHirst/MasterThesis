      subroutine dsacbnd3(excl)
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
c     000904 update according to pdg2000 lars bergstrom      
c     010214 mh2 limits corrected, joakim edsjo
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
      real*8 diffneucha,absmacha,arrmh2(25),
     & arrsin2ba(25),arrmh2cd(13),arrtanbe(13),ylim
      integer iiikkk

      data (arrmh2(i),i=1,25)/
     & 0.914500d+02,0.924500d+02,0.936000d+02,0.944500d+02,
     & 0.955000d+02,0.960000d+02,0.979000d+02,0.990500d+02,
     & 0.100750d+03,0.101050d+03,0.102050d+03,0.102800d+03,
     & 0.103300d+03,0.103900d+03,0.104350d+03,0.105100d+03, ! JE corr 010212
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
c  ... lbe update 000904 to pdg2000             
      diffneucha=abs(mass(kn(1))-mass(kcha(2)))
      absmacha=abs(mass(kcha(2)))
c ...  lbe, pdg2000, l3 limit:
      if (absmacha.lt.67.7d0) excl=ibset(excl,0) 
c ... the following is from pdg2000, abreu et al (delphi):      
      if (absmacha.lt.88.4d0.and.diffneucha.ge.3.d0
     & .and.tanbe.ge.1.d0.and.mass(ksnu(1)).ge.mass(kcha(1)))
     &  excl=ibset(excl,0)
c        pdg 94 p 1801 iv hidaka91 cdf 
c      if (abs(mass(kcha(1))).lt.99.0d0) excl=ibset(excl,0)

c--------------------------------------------lower bound on gluino mass
c  ... lbe update 000904 to pdg2000; use only limit with cascade decays:
      if (mass(kgluin).lt.173.d0) excl=ibset(excl,1)
c      excluded=.false.
c      do i=1,6
c        excluded = excluded .or. (mass(ksqu(i)).lt.mass(kgluin)) .or.
c     &        (mass(ksqd(i)).lt.mass(kgluin))
c      enddo
c      if (excluded .and. mass(kgluin).lt.260.0d0) excl=ibset(excl,1)

c------------------------------------------lower bound on squark masses
c  ... lbe update 000904 to pdg2000             
      excluded=.false.
      do i=1,6 ! aleph limit, pdg2000
        excluded = excluded .or. (mass(ksqu(i)).lt.92.d0) .or.
     &        (mass(ksqd(i)).lt.92.d0)
      enddo
      if (excluded .and. mass(kgluin)-mass(kn(kln)).lt.10.0d0)
     &  excl=ibset(excl,2)
c        pdg 94 abachi95 d0 prl75,618 is 176 with mgluino<300
      excluded=.false.
      do i=1,6 ! cdf limit, pdg2000
        excluded = excluded .or.
     &    (mass(ksqu(i)).lt.224.d0 .and. mass(ksqu(i)).gt.mass(kgluin))
     &    .or.
     &    (mass(ksqu(i)).lt.224.d0 .and. mass(ksqd(i)).gt.mass(kgluin))
      enddo
      if (excluded) excl=ibset(excl,2)

c---------------------------------lower bound on charged slepton masses
c  ... lbe update 000904 to pdg2000
      excluded=.false.
      do i=1,6 ! lep, pdg2000
        excluded = excluded .or. (mass(ksl(i)).lt.41.d0)
        enddo
c ... opal limits, pdg2000:        
      if (mass(ksl(1)).lt.87.1d0 .and. abs(mass(kn(kln))).lt.82.1d0)
     &     excl=ibset(excl,3)
      if (mass(ksl(2)).lt.87.1d0 .and. abs(mass(kn(kln))).lt.82.1d0)
     &     excl=ibset(excl,3)
      if (mass(ksl(3)).lt.82.3d0 .and. abs(mass(kn(kln))).lt.79.3d0)
     &     excl=ibset(excl,3)
      if (mass(ksl(4)).lt.82.3d0 .and. abs(mass(kn(kln))).lt.79.3d0)
     &     excl=ibset(excl,3)
      if (mass(ksl(5)).lt.81.0d0 .and. abs(mass(kn(kln))).lt.73.d0)
     &     excl=ibset(excl,3)
      if (mass(ksl(6)).lt.81.0d0 .and. abs(mass(kn(kln))).lt.73.d0)
     &     excl=ibset(excl,3)

c---------------------------------------lower bound on sneutrino masses
c ... lbe update 000904 to pdg2000 "unpublished lep limits":
      if (mass(ksnu(1)).lt.43.7d0) excl=ibset(excl,3)
      if (mass(ksnu(2)).lt.43.7d0) excl=ibset(excl,3)
      if (mass(ksnu(3)).lt.43.7d0) excl=ibset(excl,3)

c-------------------------------------lower bounds on neutralino masses
c     assumes neutralinos are ordered by increasing mass
c ... lbe update 000904, pdg2000
      if (abs(mass(kn(1))).lt.32.5d0) excl=ibset(excl,6)
      if (abs(mass(kn(2))).lt.55.9d0.and.
     &    abs(abs(mass(kn(2)))-abs(mass(kn(1)))).gt.10.d0.and.
     &    tanbe.ge.1.5d0)
     &    excl=ibset(excl,6)
      if (abs(mass(kn(3))).lt.106.6d0.and.
     &    abs(abs(mass(kn(2)))-abs(mass(kn(1)))).gt.10.d0.and.
     &    tanbe.gt.1.5d0)
     &    excl=ibset(excl,6)
c        pdg 99 abbiendi99 opal epj c8,255
      if (abs(mass(kn(2))).lt.44.0d0) excl=ibset(excl,6)
c        pdg 99 abbiendi99 opal epj c8,255
      if (abs(mass(kn(3))).lt.102.d0) excl=ibset(excl,6)
c        pdg 99 acciarri95 l3 plb350,109
      if (abs(mass(kn(4))).lt.127.0d0) excl=ibset(excl,6)

c------------------------------------------------bounds on higgs masses
c h3:
c ... lbe update, pdg2000 delphi limit:
      if(mass(kh3).lt.84.1d0) excl=ibset(excl,5)


c h2:
c...pu update 00-03-10 limit implemented from aleph 2000-006
c...http://alephwww.cern.ch/alpub/oldconf/oldconf_00.html
c...Corrected, JE 2001-02-14
      beta=atan(TanBe)
      y=sin(beta-alpha)**2
      if(mass(kh2).lt.arrmh2(1)) then
        Excl=ibset(Excl,5)
        goto 440
      elseif(mass(kh2).gt.arrmh2(25)) then
        goto 440
      else  
        do iiikkk=1,24
          if(mass(kh2).lt.arrmh2(iiikkk+1).and
     &        .mass(kh2).ge.arrmh2(iiikkk)) then
            ylim=arrsin2ba(iiikkk)+(arrsin2ba(iiikkk+1)
     &        -arrsin2ba(iiikkk))*(mass(kh2)-arrmh2(iiikkk))
     &        /(arrmh2(iiikkk+1)-arrmh2(iiikkk)) 
            goto 444
          endif
        enddo
 444    if(y.gt.ylim) Excl=ibset(Excl,5)   ! JE corr 2001-02-14
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
c      if (mass(khc).lt.77.7d0) excl=ibset(excl,5)
c ... lbe update, pdg2000, d0 limit:
      if (mass(khc).lt.82.8d0) excl=ibset(excl,5)      
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

