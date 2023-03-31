      subroutine dsacbnd7(excl)
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
c     020927 higgs limits update according to pdg2002 mia schelke
c     021001 susy part. mass limits update to pdg2002 je/ms
c     031204 standard model higgs like mh2 limit for msugra models
c     070529 calling new bsg routine+new experim bsg value, ms
c     090105 Fixed incorrect do-loop for mh3 mass bounds
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
      real*8 diffneucha,absmacha,arrmh3min(6),arrh3mintanbe(6),
     & arrmh3(18),arrh3lowtanbe(18),arrmh3cdf(10),arrh3hightanbe(10),
     & arrmh2(14),arrh2lowtanbe(14),
     & arrmh2cdf(10),arrh2hightanbe(10),ylim
      integer iiikkk

      data (arrmh3min(i),i=1,6)/3.8800d+02,3.9300d+02,4.0800d+02,
     & 4.2000d+02,4.6900d+02,5.0000d+02/

      data (arrh3mintanbe(i),i=1,6)/0.400d0,0.420d0,0.440d0,
     & 0.450d0,0.479d0,0.480d0/

      data (arrmh3(i),i=1,18)/0.9190d+02,0.9189d+02,0.9188d+02,
     & 0.9300d+02,
     & 0.9800d+02,1.0300d+02,1.0500d+02,1.1000d+02,1.1500d+02,
     & 1.2500d+02,1.3900d+02,1.5400d+02,1.6900d+02,1.9500d+02,
     & 2.3600d+02,2.8600d+02,3.5900d+02,5.0000d+02/

      data (arrh3lowtanbe(i),i=1,18)/47.500d0,30.000d0,15.240d0,
     & 10.000d0,
     & 6.660d0,7.040d0,8.000d0,8.570d0,8.000d0,   
     & 6.120d0,4.820d0,4.070d0,3.540d0,3.210d0,
     & 2.910d0,2.670d0,2.490d0,2.400d0/

      data (arrmh3cdf(i),i=1,10)/0.9190d+02,0.9800d+02,1.0100d+02,
     & 1.1000d+02,1.2000d+02,1.3000d+02,1.4000d+02,1.5100d+02,
     & 2.0100d+02,2.4700d+02/

      data (arrh3hightanbe(i),i=1,10)/47.500d0,53.000d0,54.500d0,
     & 59.000d0,68.500d0,76.500d0,87.000d0,83.500d0,
     & 96.000d0,100.000d0/



      data (arrmh2(i),i=1,14)/
     & 0.9000d+02,0.8999d+02,0.9270d+02,0.9230d+02,
     & 0.9450d+02,0.9680d+02,1.0000d+02,1.0320d+02,
     & 1.0820d+02,1.1090d+02,1.1180d+02,1.1319d+02,
     & 1.1320d+02,1.1321d+02/

      data (arrh2lowtanbe(i),i=1,14)/
     & 45.500d0,30.001d0,30.000d0,9.000d0,
     &  7.000d0,6.660d0,7.450d0,8.330d0,
     &  8.690d0,8.000d0,7.000d0,5.000d0,
     &  2.400d0,0.500d0/

      data (arrmh2cdf(i),i=1,10)/0.9000d+02,0.9720d+02,1.0060d+02,
     & 1.1500d+02,1.1830d+02,1.2000d+02,1.2220d+02,1.2380d+02,
     & 1.2500d+02,1.2560d+02/

      data (arrh2hightanbe(i),i=1,10)/45.500d0,53.000d0,55.500d0,
     & 61.500d0,65.000d0,70.000d0,83.500d0,90.000d0,
     & 91.000d0,100.000d0/



      excl=0

c----------------------------------------------bounds from b -> s gamma
c...Heavy Flavour Averaging Group 
c...http://www.slac.stanford.edu/xorg/hfag/ or arXiv:0704:3575v1
c...Barberio et al.:
c...The experimental world average
c...BR(b->s gamma) = (3.55 +- 0.24 +0.09-0.1 +- 0.03) e-4
c...               = (3.55 +- 0.26) e-4
c...For the theoretical error:
c...SM calculation, Misiak et al (hep-ph/0609232 and hep-ph/0609241)
c...+ M. Misiak and M. Steinhauser, private communication:
c...SM error: +- 0.23 e-4
c...Let's assume theoretical SUSY error is also +- 0.23 e-4
c...The total error (experiment+theory) 0.42 e-4. 
c...The two sigma confidence level is then (2.71 -> 4.39 )e-4 
c...
      call dsbsg2007full(bsg,1)  ! bsg w/ SUSY contribution
      if (bsg.lt.2.71d-4.or.bsg.gt.4.39d-4) then
         excl=ibset(excl,7)
      endif

c------------------------------------------lower bound on chargino mass
c...pdg2002 at http://pdg.lbl.gov
c...go to particle listings/other searches/susy part.
c...see the appendix in this (p.64ff) 
c...(citation:Phys.Rev.D66(2002)010001)
c...the following limits still hold 
c  ... lbe update 000904 to pdg2000    
      absmacha=min(mass(kcha(1)),mass(kcha(2)))
      diffneucha=absmacha-mass(kn(kln))
c ...  lbe, pdg2000, l3 limit:
      if (absmacha.lt.67.7d0) excl=ibset(excl,0) 
c ... the following is from pdg2000, abreu et al (delphi):      
      if (absmacha.lt.88.4d0.and.diffneucha.ge.3.d0
     & .and.tanbe.ge.1.d0.and.mass(ksnu(1)).ge.absmacha)
     &  excl=ibset(excl,0)
c        pdg 94 p 1801 iv hidaka91 cdf 
c      if (abs(mass(kcha(1))).lt.99.0d0) excl=ibset(excl,0)

c--------------------------------------------lower bound on gluino mass
c...je/ms update 02-10-01 to
c...pdg2002 at http://pdg.lbl.gov
c...go to particle listings/other searches/susy part.
c...see the appendix in this (p.64ff) 
c...(citation:Phys.Rev.D66(2002)010001)
      if (mass(kgluin).lt.195.d0) excl=ibset(excl,1)
c  ... lbe update 000904 to pdg2000; use only limit with cascade decays:
c      if (mass(kgluin).lt.173.d0) excl=ibset(excl,1)
c      excluded=.false.
c      do i=1,6
c        excluded = excluded .or. (mass(ksqu(i)).lt.mass(kgluin)) .or.
c     &        (mass(ksqd(i)).lt.mass(kgluin))
c      enddo
c      if (excluded .and. mass(kgluin).lt.260.0d0) excl=ibset(excl,1)

c------------------------------------------lower bound on squark masses
c...je/ms update 02-10-01 to
c...pdg2002 at http://pdg.lbl.gov
c...go to particle listings/other searches/susy part.
c...see the appendix in this (p.64ff) 
c...(citation:Phys.Rev.D66(2002)010001)

      excluded=.false.
      do i=1,2 ! DO limit pdg2002 
        excluded = excluded .or.
     &    mass(ksqu(i)).lt.250.d0 .or.
     &    mass(ksqd(i)).lt.250.d0        !fixed u->d, ms 021001
      enddo      
      do i=1,2 ! DO limit pdg2002 
        excluded = excluded .or.
     &    mass(ksqu(i+3)).lt.250.d0 .or.
     &    mass(ksqd(i+3)).lt.250.d0      !fixed u->d, ms 021001
      enddo                                                 
      if (excluded) excl=ibset(excl,2)
c...je/ms update 02-10-01 news:stop/sbottom limits included
      if (mass(ksb1).lt.91.0d0
     &  .and.(mass(ksb1)-mass(kn(kln))).gt.8.0d0)
     &  excl=ibset(excl,2)
      if (mass(ksb2).lt.91.0d0
     &  .and.(mass(ksb2)-mass(kn(kln))).gt.8.0d0)
     &  excl=ibset(excl,2)
      if (mass(kst1).lt.86.4d0
     &  .and.(mass(kst1)-mass(kn(kln))).gt.5.0d0)
     &  excl=ibset(excl,2)
      if (mass(kst2).lt.86.4d0
     &  .and.(mass(kst2)-mass(kn(kln))).gt.5.0d0)
     &  excl=ibset(excl,2)

c.....................??????????

c  ... lbe update 000904 to pdg2000             
c      excluded=.false.
c      do i=1,6 ! aleph limit, pdg2000
c        excluded = excluded .or. (mass(ksqu(i)).lt.92.d0) .or.
c     &        (mass(ksqd(i)).lt.92.d0)
c      enddo
c      if (excluded .and. mass(kgluin)-mass(kn(kln)).lt.10.0d0)
c     &  excl=ibset(excl,2)
c        pdg 94 abachi95 d0 prl75,618 is 176 with mgluino<300
c      excluded=.false.
c      do i=1,6 ! cdf limit, pdg2000
c        excluded = excluded .or.
c     &    (mass(ksqu(i)).lt.224.d0 .and. mass(ksqu(i)).gt.mass(kgluin))
c     &    .or.
c     &    (mass(ksqd(i)).lt.224.d0 .and. mass(ksqd(i)).gt.mass(kgluin)) !fixed u->d 
c      enddo                                                             !ms 021001
c      if (excluded) excl=ibset(excl,2)
c---------------------------------lower bound on charged slepton masses
c...je/ms update 02-10-01 to
c...pdg2002 at http://pdg.lbl.gov
c...go to particle listings/other searches/susy part.
c...see the appendix in this (p.64ff) 
c...(citation:Phys.Rev.D66(2002)010001)
c  ... lbe update 000904 to pdg2000 (still o.k. pdg2002)
      excluded=.false.
      do i=1,6 ! lep, pdg2000
        excluded = excluded .or. (mass(ksl(i)).lt.41.d0)
      enddo
c...021001update:ksl numbers corrected,limits+massdiff.updated,ALEP added 
c ... opal limits, pdg2002:        
      if (mass(ksl(1)).lt.87.1d0
     &     .and.(mass(ksl(1))-mass(kn(kln))).gt.5.0d0)
     &     excl=ibset(excl,3)
      if (mass(ksl(4)).lt.87.1d0
     &     .and.(mass(ksl(4))-mass(kn(kln))).gt.5.0d0)
     &     excl=ibset(excl,3)
      if (mass(ksl(2)).lt.82.3d0
     &     .and.(mass(ksl(2))-mass(kn(kln))).gt.3.0d0)
     &     excl=ibset(excl,3)
      if (mass(ksl(5)).lt.82.3d0
     &     .and.(mass(ksl(5))-mass(kn(kln))).gt.3.0d0)
     &     excl=ibset(excl,3)
      if (mass(ksl(3)).lt.73.0d0
     &     .and.(mass(ksl(3))-mass(kn(kln))).gt.10.0d0)
     &     excl=ibset(excl,3)
      if (mass(ksl(6)).lt.73.0d0
     &     .and.(mass(ksl(6))-mass(kn(kln))).gt.10.0d0)
     &     excl=ibset(excl,3)
c ... alep limits, pdg2002:        
      if (mass(ksl(1)).lt.95.0d0
     &   .and.(mass(ksl(1))-mass(kn(kln))).gt.15.0d0)
     &     excl=ibset(excl,3)
      if (mass(ksl(4)).lt.95.0d0
     &     .and.(mass(ksl(4))-mass(kn(kln))).gt.15.0d0)
     &     excl=ibset(excl,3)
      if (mass(ksl(2)).lt.88.0d0
     &     .and.(mass(ksl(2))-mass(kn(kln))).gt.15.0d0)
     &     excl=ibset(excl,3)
      if (mass(ksl(5)).lt.88.0d0
     &     .and.(mass(ksl(5))-mass(kn(kln))).gt.15.0d0)
     &     excl=ibset(excl,3)
      if (mass(ksl(3)).lt.76.0d0
     &     .and.(mass(ksl(3))-mass(kn(kln))).gt.15.0d0)
     &     excl=ibset(excl,3)
      if (mass(ksl(6)).lt.76.0d0
     &     .and.(mass(ksl(6))-mass(kn(kln))).gt.15.0d0)
     &     excl=ibset(excl,3)
c---------------------------------------lower bound on sneutrino masses
c...je/ms update 021001 to pdg2002(see ref above);the following still holds
c ... lbe update 000904 to pdg2000 "unpublished lep limits":
      if (mass(ksnu(1)).lt.43.7d0) excl=ibset(excl,3)
      if (mass(ksnu(2)).lt.43.7d0) excl=ibset(excl,3)
      if (mass(ksnu(3)).lt.43.7d0) excl=ibset(excl,3)

c-------------------------------------lower bounds on neutralino masses
c     assumes neutralinos are ordered by increasing mass
c...je/ms update 02-10-01 to
c...pdg2002 at http://pdg.lbl.gov
c...go to particle listings/other searches/susy part.
c...see the appendix in this (p.64ff) 
c...(citation:Phys.Rev.D66(2002)010001)
      if (abs(mass(kn(1))).lt.37.0d0) excl=ibset(excl,6)
      if (abs(mass(kn(2))).lt.62.4d0.and.tanbe.ge.1.0d0
     &    .and.tanbe.le.40.0d0)
     &    excl=ibset(excl,6)
      if (abs(mass(kn(3))).lt.99.9d0.and.tanbe.ge.1.0d0
     &    .and.tanbe.le.40.0d0)
     &    excl=ibset(excl,6)
      if (abs(mass(kn(4))).lt.116.0d0.and.tanbe.ge.1.0d0
     &    .and.tanbe.le.40.0d0)
     &    excl=ibset(excl,6)

c------------------------------------------------bounds on higgs masses
c h3:
c ... lbe update, pdg2000 delphi limit:
c      if(mass(kh3).lt.84.1d0) excl=ibset(excl,5)

c...ms update 02-09-27    
c...using limits stated for the maximal mixing scenario
c...pdg2002 at http://pdg.lbl.gov
c...look for higgs bosons under particle listings
c...see sec. 4 and fig 7 in this  
c...(citation:Phys.Rev.D66(2002)010001) 
c...for more detailed figures see
c...low tanbe:LEP preliminary limits LHWG Note 2001-04 at
c...http://lephiggs.web.cern.ch/LEPHIGGS/papers/index.html
c...we use fig 4_3
c...high tanbe:CDF limits, Phys.Rev.Lett. 86(2001)4472
c...(hep-ex/0010052) we use fig 3 (p.18)

      if(tanbe.lt.arrh3mintanbe(6)) then
       if(mass(kh3).lt.arrmh3min(1)) then
         Excl=ibset(Excl,5)
         goto 110
       elseif(mass(kh3).gt.arrmh3min(6)) then
         goto 110
       else
         ylim=0.0d0  
         do iiikkk=1,6
           if(mass(kh3).lt.arrmh3min(iiikkk+1).and
     &          .mass(kh3).ge.arrmh3min(iiikkk)) then
             ylim=arrh3mintanbe(iiikkk)+
     &          (arrh3mintanbe(iiikkk+1)-arrh3mintanbe(iiikkk))
     &          *(mass(kh3)-arrmh3min(iiikkk))
     &          /(arrmh3min(iiikkk+1)-arrmh3min(iiikkk)) 
             goto 111
           endif
         enddo
 111     if(tanbe.gt.ylim) Excl=ibset(Excl,5)  
       endif
      endif  
 110  continue

      if(tanbe.ge.arrh3mintanbe(6).and.tanbe.le.arrh3lowtanbe(18)) then      
         excl=ibset(excl,5)
      endif

      if(tanbe.lt.arrh3lowtanbe(1).and.tanbe.gt.arrh3lowtanbe(18)) then
       if(mass(kh3).lt.arrmh3(1)) then
         Excl=ibset(Excl,5)
         goto 220
       elseif(mass(kh3).gt.arrmh3(18)) then
         goto 220
       else
         ylim=0.0d0  
         do iiikkk=1,17
           if(mass(kh3).lt.arrmh3(iiikkk+1).and
     &          .mass(kh3).ge.arrmh3(iiikkk)) then
             ylim=arrh3lowtanbe(iiikkk)+
     &          (arrh3lowtanbe(iiikkk+1)-arrh3lowtanbe(iiikkk))
     &          *(mass(kh3)-arrmh3(iiikkk))
     &          /(arrmh3(iiikkk+1)-arrmh3(iiikkk)) 
             goto 222
           endif
         enddo
 222     if(tanbe.lt.ylim) Excl=ibset(Excl,5)  
       endif
      endif  
 220  continue
   

      if(tanbe.ge.arrh3hightanbe(1)
     &   .and.tanbe.lt.arrh3hightanbe(10)) then
       if(mass(kh3).lt.arrmh3cdf(1)) then
         Excl=ibset(Excl,5)
         goto 330
       elseif(mass(kh3).gt.arrmh3cdf(10)) then
         goto 330
       else
         ylim=0.0d0  
         do iiikkk=1,10
           if(mass(kh3).lt.arrmh3cdf(iiikkk+1).and
     &          .mass(kh3).ge.arrmh3cdf(iiikkk)) then
             ylim=arrh3hightanbe(iiikkk)+
     &          (arrh3hightanbe(iiikkk+1)-arrh3hightanbe(iiikkk))
     &          *(mass(kh3)-arrmh3cdf(iiikkk))
     &          /(arrmh3cdf(iiikkk+1)-arrmh3cdf(iiikkk)) 
             goto 333
           endif
         enddo
 333     if(tanbe.gt.ylim) Excl=ibset(Excl,5)  
       endif
      endif  
 330  continue
   

c h2:
c...pu update 00-03-10 limit implemented from aleph 2000-006
c...http://alephwww.cern.ch/alpub/oldconf/oldconf_00.html
c...Corrected, JE 2001-02-14

c...ms update 02-09-27    
c...using limits stated for the maximal mixing scenario
c...pdg2002 at http://pdg.lbl.gov
c...look for higgs bosons under particle listings 
c...see sec 4 and fig 7 in this (use this fig for tanbe~30)
c...(citation:Phys.Rev.D66(2002)010001) 
c...for more detailed figures see
c...for low tanbe:LEP preliminary limits LHWG Note 2001-04 at
c...http://lephiggs.web.cern.ch/LEPHIGGS/papers/index.html
c...we use fig 4_2
c...for high tanbe:CDF limits, Phys.Rev.Lett. 86(2001)4472
c...(hep-ex/0010052) we use fig 2 (p.17)



      if(tanbe.lt.arrh2lowtanbe(1)) then
       if(mass(kh2).lt.arrmh2(1)) then
         Excl=ibset(Excl,5)
         goto 440
       elseif(mass(kh2).gt.arrmh2(14)) then
         goto 440
       else
         ylim=0.0d0  
         do iiikkk=1,14
           if(mass(kh2).lt.arrmh2(iiikkk+1).and
     &          .mass(kh2).ge.arrmh2(iiikkk)) then
             ylim=arrh2lowtanbe(iiikkk)+
     &          (arrh2lowtanbe(iiikkk+1)-arrh2lowtanbe(iiikkk))
     &          *(mass(kh2)-arrmh2(iiikkk))
     &          /(arrmh2(iiikkk+1)-arrmh2(iiikkk)) 
             goto 444
           endif
         enddo
 444     if(tanbe.lt.ylim) Excl=ibset(Excl,5)  
       endif
      endif  
 440  continue


      if(tanbe.ge.arrh2hightanbe(1)
     &   .and.tanbe.lt.arrh2hightanbe(10)) then
       if(mass(kh2).lt.arrmh2cdf(1)) then
         Excl=ibset(Excl,5)
         goto 550
       elseif(mass(kh2).gt.arrmh2cdf(10)) then
         goto 550
       else
         ylim=0.0d0  
         do iiikkk=1,10
           if(mass(kh2).lt.arrmh2cdf(iiikkk+1).and
     &          .mass(kh2).ge.arrmh2cdf(iiikkk)) then
             ylim=arrh2hightanbe(iiikkk)+
     &          (arrh2hightanbe(iiikkk+1)-arrh2hightanbe(iiikkk))
     &          *(mass(kh2)-arrmh2cdf(iiikkk))
     &          /(arrmh2cdf(iiikkk+1)-arrmh2cdf(iiikkk)) 
             goto 555
           endif
         enddo
 555     if(tanbe.gt.ylim) Excl=ibset(Excl,5)  
       endif
      endif  
 550  continue

c change this when we set a global variable to say whether we are in the
c Mssm or in mSUGRA
c limit on a standard model like h2, standard model limit from 
c CERN-EP-2003-011 at: http://lephiggs.web.cern.ch/LEPHIGGS/papers/index.html
c

      if(dabs(dabs(sgnmuvar)-1.d0).lt.1.d-10) then
        if(mass(kh2).lt.114.4d0) Excl=ibset(Excl,5) 
      endif
 

c hc:
c...ms update 02-09-27    
c...pdg2002 at http://pdg.lbl.gov
c...look for higgs bosons under particle listings 
c...see sec 5 (citation:Phys.Rev.D66(2002)010001) 
c...for more details see
c...LEP preliminary limits LHWG Note 2001-05 at
c...http://lephiggs.web.cern.ch/LEPHIGGS/papers/index.html
      if (mass(khc).lt.78.6d0) excl=ibset(excl,5)   

c...pu update 00-03-10. limit implemented from aleph 2000-011
c...http://alephwww.cern.ch/alpub/oldconf/oldconf_00.html      
c      if (mass(khc).lt.77.7d0) excl=ibset(excl,5)
c ... lbe update, pdg2000, d0 limit:
c      if (mass(khc).lt.82.8d0) excl=ibset(excl,5)      
c-----------------------------------------bound from gamma_z(invisible)
      mz = mass(kz)
      mz2 = mz*mz
      ! 3 neutrinos, corrected 2002-10-01, factor of 3 added
      ! (Thanks to Dan Hooper)
      gzinv = 3.0d0 * (g2wmz/2/cwmz)**2 * mz/(24.d0*pi)
      ! 3 sneutrinos
      do i=1,3
         temp = 1.d0-4.d0*mass(ksnu(i))**2/mz2
         if (temp.gt.0.d0) then
            temp = temp**(1.5d0) *
     &           (g2wmz/2/cwmz)**2 * mz/(48.d0*pi)
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
               gzij = g2wmz/(2.d0*cwmz)*
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
c        pdg 2002 (499 +- 1.5) => 2 sigma limit, 499 +- 3
      if (gzinv.gt.0.502d0) excl=ibset(excl,4)

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

