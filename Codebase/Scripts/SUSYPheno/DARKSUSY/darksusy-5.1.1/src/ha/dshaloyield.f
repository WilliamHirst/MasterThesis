*****************************************************************************
***   function dshaloyield gives the total yield of positrons, cont. gammas
***   or neutrinos coming from neutralino annihilation in the halo
***   the yields are given as number / annihilation. the energy egev
***   is the threshold for integrated yields and the energy for
***   differential yields. the yields are
***     yieldk =  51: integrated positron yields
***     yieldk =  52: integrated cont. gammas
***     yieldk =  53: integrated muon neutrinos
***     yieldk =  54: integrated antiproton yields
***     yieldk =  71: integrated neutrino yields (same as 53)
***     yieldk =  72: integrated muon yields at creation
***     yieldk =  73: integrated muon yields in ice
***     yieldk = above+100: differential in energy
*** the annihilation branching ratios and scalar parameters
*** are extracted from dshacom.h. Note. You need to call dshasetup prior
*** to calling this routine.
*** istat will be zero in case of no errors. Otherwise its bits are set as
***   bit  decimal   reason
***     0        1   some inaccesible parts the differential muon spectra
***                  has been wanted, and the returned yield should then
***                  be treated as a lower bound.
***     1        2   energetically forbidden annihilation channels have been
***                  wanted.
***     2        4   problems with dsIBf_intdxdy integration for IB yields
***                  (only used if option haib is set to 'susy')
***     3        8   problems with dsIBf_intdy integration for IB yields
***                  (only used if option haib is set to 'susy')

*** author: joakim edsjo  edsjo@physics.berkeley.edu
*** date: 98-01-29
*** modified: 98-04-15
*** modified: 2007-05-01 Torsten Bringmann (IB contribution added)
*** modified: 2008-01-15 Joakim Edsjo, made more modular
*****************************************************************************

      real*8 function dshaloyield(egev,yieldk,istat)
      implicit none
      include 'dshacom.h'
      include 'dshrcom.h'
      include 'dsidtag.h'

c------------------------ functions ------------------------------------

      real*8 dshayield,dshaib

c------------------------ variables ------------------------------------

      real*8 egev,yield,sv
      integer ch,istat,yieldk

c----------------------------------------------- set-up common variables

      if (.not.dshasetupcalled) then
        write(*,*) 'DS: Error in dshaloyield: dshasetup must be called',
     &    ' before any halo yield calculation'
        write(*,*) 'begins for every new model. Program stopping.'
        stop
      endif

c...loop through different channels and calculate the yield above threshold
c...for each channel.

c      write(*,*)
c      write(*,*) 'model: ',idtag,'  eth = ',egev
      haistat=0
      yield=0.0d0
      do 100 ch=1,29
        if (habr(ch).gt.0.0d0) then
          yield=yield+habr(ch)*dshayield(hamwimp,egev,
     &      ch,yieldk,istat)
          haistat=or(haistat,istat)
        endif
  100 continue


c...add IB contribution
      yield=yield+dshaib(egev,yieldk,istat)
      haistat=or(haistat,istat*4)

      dshaloyield=yield
      hristat=or(hristat,haistat)
      istat=haistat

      end



















