*****************************************************************************
*** function dshayieldf calculates the yield above threshold (or differential
*** at that energy) for the annihilation channel ch and the
*** fluxtype given by yieldk, according to the following table.
***
*** particle       integrated yield     differential yield
*** --------       ----------------     ------------------
*** positron                     51                    151
*** cont. gamma                  52                    152
*** nu_mu and nu_mu-bar          53                    153
*** antiproton                   54                    154
*** cont. gamma w/o pi0          55                    155
*** nu_e and nu_e-bar            56                    156
*** nu_tau and nu_tau-bar        57                    157
*** pi0                          58                    158
*** nu_mu and nu_mu-bar          71                    171 (same as 53/153)
*** muons from nu at creation    72                    172
*** muons from nu at detector    73                    173
***
*** only channels chi = 1-8 are supported.
*** chi =  1 - c c-bar
***        2 - b b-bar
***        3 - t t-bar
***        4 - tau+ tau-
***        5 - W+ W-
***        6 - Z0 Z0
***        7 - mu+ mu-
***        8 - gluon gluon
*** units: (annihilation)**-1  integrated
*** units: gev**-1 (annihilation)**1 differential
***
*** Note: if this routine is called directly, without calling dshayield
*** first, one needs to load the correct yield tables manually first
*** with a call to dshainit(yieldk) (only needs to be done once).
*** If this routine is called outside of the kinematical regions
*** where tables exist, the following is done:
***   - for lower energies, than the lowest simulted ones,
***     extrapolations are used
***   - for masses below the lowest simulated, extrapolations
***     are used
***   - for energies above the highest simulated, the results
***     for the highest energy simulated are used
***
*** Modifications: 2010-11-05 (JE): better extrapolation below lowest
***   simulated energy (now dN/dz is taken as constant below)
*****************************************************************************

      real*8 function dshayieldf(mneu,emuthr,ch,yieldk,istat)
      implicit none
      include 'dshacom.h'
      include 'dsidtag.h'
      include 'dsio.h'

c------------------------ variables ------------------------------------

      real*8 mneu,emuthr,phi1,phi2,mp1,mp2,zpl,mn,z,
     &  tmp,lge,dshayieldget,scalefactor
      integer ch,istat,zi,m1i,m2i,fltype,fi,yieldk
      logical wb
      parameter(lge=0.434294481903d0)
      external dshacom  ! set up common block variables

c-----------------------------------------------------------------------

      dshayieldf=0.d0

c      write(*,*) 'dshayieldf called with: mneu = ',mneu
c      write(*,*) '  emuthr = ',emuthr
      call dshadec(yieldk,fltype,fi)

      wb=.true.

      if (emuthr.ge.mneu) then
        dshayieldf=0.0d0
        return
      endif

      mn = mneu
c...take care of the case where mneu is between the mass of the annihilation
c...product and the lower bound of the simulations (there might be a small
c...gap of less than a gev) for tt-bar,ww and zz.
      scalefactor=1 ! use yields directly as they are
      if (mn.lt.lb(ch)) then
         if (mn.lt.0.98d0*lb(ch)) istat=(istat/2)*2+1
         mn=lb(ch)
         scalefactor=mneu/mn ! rescale energies, i.e. assume dN/dz is OK
      endif

c      if (ch.eq.3.or.ch.eq.5.or.ch.eq.6) then
c        if (mn.lt.lb(ch).and.mn.gt.(0.99*lb(ch))) then
c          mn=lb(ch)
c        endif
c      endif

c...check if mneu within correct bounds
c...this part is obsolete, as we extrapolate now
c      if (mn.lt.lb(ch).and.emuthr.lt.lb(ch)) then
c        if (prtlevel.gt.0) then
c          write(6,*)
c          write(6,*) 'warning in dshayieldf for model: ',idtag
c          write(*,5000) 'a neutralino mass of ',mn,
c     +      ' gev wants to be used,'
c          write(6,5010) 'while the lower bound for channel ',ch,
c     +      ' is ',lb(ch),' gev.'
c          write(6,*) 'the yield is put to 0.0 for these too low masses.'
c          write(6,*) 'the results can thus only be trusted as a',
c     +      ' lower bound.'
c        endif
c        wb=.false.
c        dshayieldf=0.0d0
c        istat=(istat/2)*2+1
c     endif

      if (mn.gt.ub(ch)) then
        if (prtlevel.gt.0) then
          write(6,*)
          write(6,*) 'warning in dshayieldf for model: ',idtag
          write(6,5000) 'a neutralino mass of ',mn,
     +      ' gev wants to be used,'
          write(6,5010) 'while the upper bound for channel ',ch,
     +      ' is ',ub(ch),' gev.'
          write(6,5020) 'a neutralino mass of ',ub(ch),' gev is used',
     +      ' instead for these too high masses.'
          write(6,*) 'the results can thus only be trusted as a',
     +      ' lower bound.'
        endif
        istat=(istat/2)*2+1
      endif


c---------------------------------------------------- integrated yields
      if (wb.and.fltype.eq.1) then
c...determine which entries in phiint to use and how
        if (fi.ge.1.and.fi.le.20) then ! log tabulated yields
          z=(log10(emuthr/mneu)+ndec)/ndec
        else
          z=emuthr/mneu
        endif
        if (z.lt.0.0d0) then
          dshayieldf=0.0d0
          return
        endif
        call dshaifind(z,zindex(0,1),zpl,zi,0,zn-1)

        if (zi.eq.-5.or.zi.ge.zn) then
          dshayieldf=0.0d0
          return
        endif

        call dshaifind(mn,mi(1),tmp,m1i,1,17) ! mn here
        mp1=mi(m1i)
        m2i=m1i+1
        mp2=mi(m2i)

        if (mn.ge.mi(18)) then
          m1i=18
          m2i=18
          mp1=mi(18)
          mp2=mp1

          dshayieldf =
     &      (1.0-zpl)*dshayieldget(zi,m1i,ch,fi,fltype,istat)+
     &      zpl*dshayieldget(zi+1,m1i,ch,fi,fltype,istat)
        else
          phi1 =
     &      (1.0-zpl)*dshayieldget(zi,m1i,ch,fi,fltype,istat)+
     &      zpl*dshayieldget(zi+1,m1i,ch,fi,fltype,istat)
          phi2 =
     &      (1.0-zpl)*dshayieldget(zi,m2i,ch,fi,fltype,istat)+
     &      zpl*dshayieldget(zi+1,m2i,ch,fi,fltype,istat)
          dshayieldf = phi1 + (phi2-phi1)*(mn-mp1)*(mn+mp1)/
     &      ((mp2-mp1)*(mp2+mp1))
        endif
c        if (mn.ne.mneu) then
c          dshayieldf=mneu**2/mn**2*dshayieldf
c        endif
      endif

c-------------------------------------------------- differential yields
      if (wb.and.fltype.eq.2) then
c...determine which entries in phidiff to use and how
        if (fi.ge.1.and.fi.le.20) then ! log tabulated yields
          z=(log10(emuthr/mneu)+ndec)/ndec
        else
          z=emuthr/mneu
        endif
        if (z.lt.0.0d0) then
          dshayieldf=0.0d0
          return
        endif
        call dshaifind(z,zindex(-1,2),zpl,zi,-1,zn-1)

        if (zi.eq.-5.or.zi.ge.zn) then
          dshayieldf=0.0d0
          return
        endif

        if ((emuthr.gt.mn).or.(emuthr.le.0.0).or.zi.eq.-5) then
          dshayieldf=0.0d0
          return
        endif

        call dshaifind(mn,mi(1),tmp,m1i,1,17) ! here we use mn
        mp1=mi(m1i)
        m2i=m1i+1
        mp2=mi(m2i)

        if (mn.ge.mi(18)) then
          m1i=18
          m2i=18
          mp1=mi(18)
          mp2=mp1
          dshayieldf =
     &      (1.0-zpl)*dshayieldget(zi,m1i,ch,fi,fltype,istat)+
     &      zpl*dshayieldget(zi+1,m1i,ch,fi,fltype,istat)
        else
          phi1 =
     &      (1.0-zpl)*dshayieldget(zi,m1i,ch,fi,fltype,istat)+
     &      zpl*dshayieldget(zi+1,m1i,ch,fi,fltype,istat)
          phi2 =
     &      (1.0-zpl)*dshayieldget(zi,m2i,ch,fi,fltype,istat)+
     &      zpl*dshayieldget(zi+1,m2i,ch,fi,fltype,istat)
          dshayieldf = phi1 + (phi2-phi1)*(mn-mp1)*(mn+mp1)/
     &      ((mp2-mp1)*(mp2+mp1))
        endif
c        if (mn.ne.mneu) then
c          dshayieldf=mneu**2/mn**2*dshayieldf
c        endif
c... convert from dyield/dz or dyield/dx to dyield/de
        if (fi.ge.1.and.fi.le.20) then ! log tabulated yields
          dshayieldf=dshayieldf*lge/(ndec*emuthr)
        else
          dshayieldf=dshayieldf/mneu
        endif
      endif


 5000 format(' ',a,f8.2,a)
 5010 format(' ',a,i2,a,f8.2,a)
 5020 format(' ',a,f8.2,a,a)

      end

