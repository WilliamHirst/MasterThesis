      subroutine dsntrates(emuth0,thmax0,rtype,ptype,rateea,
     &  ratesu,istat)
c_______________________________________________________________________
c
c            n e u t r a l i n o   b r a n c h i n g   r a t i o s
c            a n d  c a p t u r e  r a t e  i n  t h e   s u n
c            m u o n   f l u x  c a l c u l a t e d
c    november, 1995
c    uses routines by p. gondolo and j. edsjo
c    modified by l. bergstrom and j. edsjo
c    capture rate routines are written by l. bergstrom
c    input:  emuth0  - muon energy threshold in gev
c            thmax0  - muon angel cut in degrees
c            rtype   - 1 = flux of neutrinos (nu_mu and/or nu_mu-bar) km^-2 yr^-1
c                      2 = contained mu- and/or mu+ events km^-3 yr^-1
c                      3 = through-going mu- and/or mu+ events km^-2 yr^-1
c            ptype   - 1 = particles only (nu_mu or mu-)
c                      2 = anti-particles only (nu_mu-bar or mu+)
c                      3 = summed rates (both particles and anti-particles)
c    hidden input: ntcalcmet - 1 use jkg approximations
c                              2 use jkg for sun, full gould for earth
c                              3 use jkg for sun, full gould+dk for earth
c                              4 use full numerical calculations for Sun, Earth
c    output: rateea  - events from earth ann. per km^2(3) per yr
c            ratesu  - events from sun ann. per km^2(3) per yr
c    slightly modified by j. edsjo.
c    modified by j. edsjo 97-05-15 to match new inv. rate convention
c    modified by j. edsjo 97-12-03 to match muflux3.21 routines.
c    modified by p. gondolo 98-03-04 to detach dsntannrate from susy
c    routines.
c    modified by j. edsjo 98-09-07 to fix istat bug.
c    modified by j. edsjo 98-09-23 to use damour-krauss distributions
c      and full earth formulas.
c    modified by j. edsjo 99-03-17 to include better damour-krauss
c      velocity distributions and numerical capture rate integrations
c      for these non-gaussian distributions
c    modified by p. scott 11-04-23 to allow rates to be calculated
c      for only particles or only antiparticles
c
c=======================================================================
      implicit none
      include 'dshmcom.h'
      include 'dsntcom.h'
      include 'dswacom.h'

      real*8 emuth0,thmax0,rateea,ratesu,arateea,aratesu,
     &  yield,tmp1,tmp2
      real*8 dsntmuonyield
      integer istat,itmp,rtype,ptype

c ----------------------------------------- zero common block data

      tausu=0.0d0
      csu=0.0d0
      tauea=0.0d0
      cea=0.0d0
      ceadk=0.0d0
      gtot10=0.0d0
      ntarateea=0.0d0
      ntaratesu=0.0d0

      aratesu=0.0d0
      arateea=0.0d0

      rateea=0.0d0
      ratesu=0.0d0

c --------------------------------------------- start calculations

      if (rhox.eq.0.0d0) then
        rateea=0.0d0
        ratesu=0.0d0
        return
      endif

c ------------------------------------ annihilation rates at rest:

      

c **************************************************************

      if (wamwimp.gt.emuth0) then

        if (ntcalcmet.eq.1.or.ntcalcmet.eq.2
     &    .or.ntcalcmet.eq.4) then ! jkg and/or gould w/ Gauss or full dist.
          call dsntannrate(wamwimp,wasigsip,wasigsdp,wasv,
     &       arateea,aratesu)

        elseif (ntcalcmet.eq.3) then ! jkg for sun, gould+dk for earth
          call dsntdkannrate(wamwimp,wasigsip,wasigsdp,wasv,arateea,
     &      aratesu)
        else
          write(*,*) 'error in dsntrates: invalid option,',
     &      ' ntcalcmet = ',ntcalcmet
          return
        endif
c...arateea and aratesu in units of 10^24 yr^-1

        yield=dsntmuonyield(emuth0,thmax0,'su',1,rtype,ptype,istat)
c...yield is in units of 10^-30 m^-2(3).
        ratesu=yield*aratesu
c...we now have units km^-2 yr^-1 if rtype=3.
c...one more m^-1 -> km^-1 still to go for rtype=2
        if (rtype.eq.2) then
          ratesu=ratesu*1.0d3
        endif
        itmp=istat
        yield=dsntmuonyield(emuth0,thmax0,'ea',1,rtype,ptype,istat)
c...yield is in units of 10^-30 m^-2(3).
        rateea=yield*arateea
c...we now have units km^-2 yr^-1 if rtype=3.
c...one more m^-1 -> km^-1 still to go for rtype=2
        if (rtype.eq.2) then
          rateea=rateea*1.0d3
        endif
        istat=or(istat,itmp)
      else
        rateea=0.d0
        ratesu=0.d0
      endif

      ntarateea=arateea
      ntaratesu=aratesu

      return

      end




