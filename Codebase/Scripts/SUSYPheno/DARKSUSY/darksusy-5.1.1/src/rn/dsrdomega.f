      real*8 function dsrdomega(omtype,fast,xf,ierr,iwar,nfc)

**********************************************************************
*** function dsrdomega calculates omega h^2 for the mssm neutralino
*** uses the mssm routines and the relic density routines
*** input:
***   omtype = 0 - no coann
***            1 - include all relevant coannihilations (charginos, 
***                neutralinos and sleptons)
***            2 - include only coannihilations betweeen charginos
***                and neutralinos
***            3 - include only coannihilations between sfermions
***                and the lightest neutralino
***   fast =   0 - standard accurate calculation (accuracy better than 1%)
***            1 - faster calculation: (recommended unless extreme accuracy
***                is needed).
***                * requires less accuracy in tabulation of w_eff
***            2 - quick and dirty method, i.e. expand the annihilation
***                cross section in x (not recommended)
*** output:
***   dsrdomega - omega h^2 for the neutralino
***   xf = m_WIMP / T_f where T_f is the freeze-out temperature
***   ierr = error from dsrdens or dsrdqad
***   iwar = warning from dsrdens of dsrdqad
***   nfc - number of function calls to the effective annihilation
***         cross section
*** authors: joakim edsjo, paolo gondolo
*** date: 98-03-03
*** modified: 98-03-03
***           99-07-30 pg
***           02-02-27 joakim edsjo: including sfermion coanns
***           06-02-22 paolo gondolo: streamlined inclusion of coanns
**********************************************************************

      implicit none
      include 'dsmssm.h'
      include 'dsio.h'
      include 'dsidtag.h'
      include 'dsandwcom.h'
      include 'dsrdcom.h'
c      include 'dsrncom.h'

      integer comax
      parameter(comax=50)
      real*8 oh2,mcoann(comax),dof(comax),dsanwx,tm(comax),
     &     rm(comax),rw(comax),xf,tf,h1width
      integer ncoann,ierr,iwar,nthr,omtype,nfc,fast,fastsav,i
      integer kpart(comax)
      external dsanwx
      logical first
      data first/.true./,fastsav/0/
      save first,fastsav

c----------------------------------------------------------------------

c...Check which coannihilation processes to include
c...The coproc switch determines which processes are included
c...Bit Decimal Coannihilation processes included
c...  0       1 neu_1, neu_2, neu_3, neu_4, cha_1, cha_2
c...  1       2 sel_1,sel_2,smu_1,smu_2,stau_1,stau_2
      if (omtype.eq.0) then      ! no coannihilations
        coproc=0
      elseif (omtype.eq.1) then  ! include all coannihilations
        coproc=3
      elseif (omtype.eq.2) then  ! include only charg. and neu coanns
        coproc=1
      elseif (omtype.eq.3) then  ! include only slepton coanns
        coproc=2
      else
        write(*,*) 'ERROR in dsrdomega: invalid omtype=',omtype
        stop
      endif

      oh2=0.0d0
      xf=0.0d0
      nfc=0

c...temporarily change width of h1 higgs boson
      h1width=width(kh1)
      width(kh1)=max(width(kh1),0.1d0) ! je, dec 11, 1998

      if (first) then
         call dsrdinit
         first=.false.
      endif

      if (fast.eq.1) then       ! fast calculation
c         dwopt=.true.
         waccd=0.05d0
         dpminr=5.0d-4
         dpthr=2.5d-3
         wdiffr=0.5d0
         wdifft=0.1d0
         mcofr=1.5d0
c         brmin=1.0d-4           ! 1.0d-3 doesn't make a big speed difference
      else
         mcofr=2.1d0
      endif

      if (fast.eq.0.and.fastsav.eq.1) then
         write (*,*) 'error in dsrdomega:',
     &    ' can''t go back to fast=0 when fast=1 had been used.'
c        write(*,*) 'rewrite dsrdomega if this is desired.'
         stop
      endif

      fastsav=fast
c--------------------------------------initialize relic density routines
      rdprt = prtlevel

      coann=0
      copart=0

c...Add lightest neutralino to annihilation list
      ncoann=1
      mcoann(1)=mass(kn(1))
      dof(1)=kdof(kn(1))
      kpart(1)=kn(1)

c...Add other neutralinos and charginos if chosen to be included
      if (and(coproc,1).ne.0) then
         coann=1
         copart=ibset(copart,0) ! neu_1 added to list of particles
         mcoann(2)=mass(kn(2))
         mcoann(3)=mass(kn(3))
         mcoann(4)=mass(kn(4))
         dof(2)=kdof(kn(2))
         dof(3)=kdof(kn(3))
         dof(4)=kdof(kn(4))
         if (mcoann(2)/mcoann(1).le.mcofr) then
            ncoann=2
            kpart(ncoann)=kn(2)
            copart=ibset(copart,1)
         endif
         if (mcoann(3)/mcoann(1).le.mcofr) then
            ncoann=3
            kpart(ncoann)=kn(3)
            copart=ibset(copart,2)
         endif
         if (mcoann(4)/mcoann(1).le.mcofr) then
            ncoann=4
            kpart(ncoann)=kn(4)
            copart=ibset(copart,3)
         endif
         if (mass(kcha(1))/mcoann(1).le.mcofr) then
            ncoann=ncoann+1
            mcoann(ncoann)=mass(kcha(1))
            dof(ncoann)=kdof(kcha(1))
            kpart(ncoann)=kcha(1)
            copart=ibset(copart,4)
         endif
         if (mass(kcha(2))/mcoann(1).le.mcofr) then
            ncoann=ncoann+1
            mcoann(ncoann)=mass(kcha(2))
            dof(ncoann)=kdof(kcha(2))
            kpart(ncoann)=kcha(2)
            copart=ibset(copart,5)
         endif
      endif

c...Add slepton coannihilations if chosen to be included
      if(and(coproc,2).ne.0) then

c...Sleptons
         do i=1,6
            if (mass(ksl(i))/mcoann(1).le.mcofr) then
c               write(*,*) idtag,': ','Slepton ',i,' added..'
               ncoann=ncoann+1
               mcoann(ncoann)=mass(ksl(i))
               dof(ncoann)=kdof(ksl(i))
               kpart(ncoann)=ksl(i)
               copart=ibset(copart,i+5)
            endif
         enddo

c...Sneutrinos
         do i=1,3
            if (mass(ksnu(i))/mcoann(1).le.mcofr) then
c               write(*,*) idtag,': ','Sneutrino ',i,' added..'
               ncoann=ncoann+1
               mcoann(ncoann)=mass(ksnu(i))
               dof(ncoann)=kdof(ksnu(i))
               kpart(ncoann)=ksnu(i)
               copart=ibset(copart,i+11)
            endif
         enddo
c...Up squarks
         do i=1,6
            if (mass(ksqu(i))/mcoann(1).le.mcofr) then
c               write(*,*) idtag,': ','Up squark ',i,' added..'
               ncoann=ncoann+1
               mcoann(ncoann)=mass(ksqu(i))
               dof(ncoann)=kdof(ksqu(i))
               kpart(ncoann)=ksqu(i)
               copart=ibset(copart,i+14)
            endif
         enddo

c...Down squarks
         do i=1,6
            if (mass(ksqd(i))/mcoann(1).le.mcofr) then
c               write(*,*) idtag,': ','Down squark ',i,' added..'
               ncoann=ncoann+1
               mcoann(ncoann)=mass(ksqd(i))
               dof(ncoann)=kdof(ksqd(i))
               kpart(ncoann)=ksqd(i)
               copart=ibset(copart,i+20)
            endif
         enddo

      endif

c...add resonances for new call
cpg      call dsrdres(ncoann,mcoann,nrs,rm,rw)
      nres=0
c...z resonance
      if (mass(kz).gt.mcoann(1)*2.0d0) then
        nres=nres+1
        rm(nres)=mass(kz)
        rw(nres)=width(kz)
      endif
c...h1 resonance
      if (mass(kh1).gt.mcoann(1)*2.0d0) then
        nres=nres+1
        rm(nres)=mass(kh1)
        rw(nres)=width(kh1)
      endif
c...h2 resonance
      if (mass(kh2).gt.mcoann(1)*2.0d0) then
        nres=nres+1
        rm(nres)=mass(kh2)
        rw(nres)=width(kh2)
      endif
c...h3 resonance
      if (mass(kh3).gt.mcoann(1)*2.0d0) then
        nres=nres+1
        rm(nres)=mass(kh3)
        rw(nres)=width(kh3)
      endif
c...coannihilation-resonances
      if (ncoann.gt.1) then
        if (mass(khc).gt.mcoann(1)*2.0d0) then   ! sufficient condition
          nres=nres+1
          rm(nres)=mass(khc)
          rw(nres)=width(khc)
        endif
        if (mass(kw).gt.mcoann(1)*2.0d0) then ! sufficient condition
          nres=nres+1
          rm(nres)=mass(kw)
          rw(nres)=width(kw)
        endif
      endif

c...add thresholds for new call
cpg      call dsrdthr(ncoann,mcoann,nt,tm)
      nthr=0
c...ww-threshold
      if (mass(kw).gt.mcoann(1)) then
        nthr=nthr+1
        tm(nthr)=2.0d0*mass(kw)
      endif
c...zz-threshold
      if (mass(kz).gt.mcoann(1)) then
        nthr=nthr+1
        tm(nthr)=2.0d0*mass(kz)
      endif
c...tt-bar-threshold
      if (mass(kt).gt.mcoann(1)) then
        nthr=nthr+1
        tm(nthr)=2.0d0*mass(kt)
      endif
c...note that coannihilation thresholds are automatically added in
c...dsrdstart (called from dsrdens).

      rdtag=idtag

      if (fast.eq.0.or.fast.eq.1) then  ! standard accurate method
         call dsrdens(dsanwx,ncoann,kpart,mcoann,dof,nres,rm,rw,
     &     nthr,tm,oh2,tf,ierr,iwar)
         nfc=nr
      elseif (fast.eq.2) then   ! quick and dirty method
         call dsrdqad(dsanwx,mcoann(1),oh2,ierr)
         iwar=0
         nfc=2
      else 
         write(*,*) 'ERROR in dsrdomega: invalid option fast=',fast 
         stop
      endif

c...  the following two lines are added by je for test
      if (prtlevel.eq.-1) then
         call dsrdwres            ! je test
         write (6,*) '============== new point =============='
         call dsrdwrate(6,6,1)
         call dsrdwrate(6,6,2)
         call dsrdwrate(6,6,3)
      endif

c...  the following lines are added by pg for test 050528
      if (prtlevel.eq.-111) call dsrdwintrp(dsanwx,50)

      dsrdomega = oh2
      xf=mco(1)/tf

      width(kh1)=h1width        ! change back h1 width

      end
