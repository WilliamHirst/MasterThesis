*****************************************************************************
***   function dsIByield gives the yield from internal
***   bremsstrahlung (IB) for SUSY models. 
***   Input: egev - energy in GeV
***          yieldk - which yield to calculate
***             (see dshaloyield for an explanation, 
***              currently only photon and positron yield are implemented)
***   Output: yield (annihilation)**1, differential also GeV**-1
***           istat is set as follows in case of errors
***   bit  decimal  reason
***     0        1  dsIBf_intdy (integration for photon yield) failed
***     1        2  dsIBf_intdy2 (integration for positron yield) failed
***                              -- for direct annihilation into positrons (channel eeg)
***     2        4  dsIBf_intdy2 (integration for positron yield) failed
***                              -- for annihilation channel different from eeg
*** Author: Joakim Edsjo, edsjo@physto.se
***         Torsten Bringmann, bringman@sissa.it
*** Date: 2008-01-15
*****************************************************************************

      real*8 function dsIByield(egev,yieldk,istat)
      implicit none
      include 'dshacom.h'
      include 'dshrcom.h'
      include 'dsidtag.h'
      include 'dsprep.h'
      include 'dsibcom.h'

c------------------------ functions ------------------------------------

      real*8 dsIByieldone

c------------------------ variables ------------------------------------

      real*8 egev,yield,sv,dssigmav
      real*8 IBabr(12)
      integer IBch,istat,yieldk,itmp

c----------------------------------------------- set-up common variables


      istat=0

      yield=0.0d0

c... add IB contribution

      sv=dssigmav(0)
      if (yieldk.eq.52.or.yieldk.eq.152.or.
     &    yieldk.eq.51.or.yieldk.eq.151) then 
        IBabr(1)=dssigmav(13)/sv    ! W+W-
        IBabr(2)=dssigmav(11)/sv    ! W+H- and W-H+
        IBabr(3)=1.d0               ! H+H- : for neutralino annihilations, the 
                                    ! lowest order contribution vanishes 
                                    ! identically for v->0, but we still want to 
                                    ! keep this channel.The IB contribution is 
                                    ! therefore instead normalized to the *total* 
                                    ! annihilation cross section.
        IBabr(4)=dssigmav(15)/sv    ! e+ e-
        IBabr(5)=dssigmav(17)/sv    ! mu+ mu-
        IBabr(6)=dssigmav(19)/sv    ! tau+ tau-
        IBabr(7)=dssigmav(20)/sv    ! u u-bar
        IBabr(8)=dssigmav(21)/sv    ! d d-bar
        IBabr(9)=dssigmav(22)/sv    ! c c-bar
        IBabr(10)=dssigmav(23)/sv   ! s s-bar
        IBabr(11)=dssigmav(24)/sv   ! t t-bar
        IBabr(12)=dssigmav(25)/sv   ! b b-bar

        if (ibhow.eq.2) call dsibselect ! dynamically choose important channels
        do 200 IBch=1,12     
          if ((IBflag(IBch).ne.0).and.(IBabr(IBch).gt.0.0d0)) then  
            yield=yield+IBabr(IBch)*
     &                   dsIByieldone(egev,IBch,yieldk,itmp)
            istat=or(istat,itmp)

           endif
  200   continue
      endif

      dsIByield=yield

      end
