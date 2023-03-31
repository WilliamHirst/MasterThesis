      subroutine dsddset(sisd,cs)
c...set parameters for scattering cross section
c...  c - character string specifying choice to be made
c...author: paolo gondolo 2000-07-07
c...modified by Gintaras Duda 2007-06-27 for new FF options
c...modified by Paolo Gondolo 2008-02-18
      implicit none
      include 'dsddcom.h'
      character*(*) sisd,cs
      logical goodoption

      goodoption=.false.

c...help
      if (sisd.eq.'help') then
         goodoption=.true.
         write (*,*) 'dsddset: use "call dsddset(key,value)" where'
         write (*,*) ' key=''help'' this help'
         write (*,*) ' key=''si'' for spin-independent scattering'
         write (*,*) ' key=''sd'' for spin-dependent scattering'
         write (*,*) ' for help on a specific key, use e.g.',
     &        ' "call dsddset(''si'',''help'')" '
      endif
      if (sisd.eq.'si'.and.cs.eq.'help') then
         goodoption=.true.
         write (*,*) 'dsddset: use "call dsddset(''si'',value)" where'
         write (*,*) ' value=''gls91'' for <mq q qbar> as in ',
     &        'Gaisser, Leutwyler & Sainio 1991 (GLS)'
         write (*,*) ' value=''dn93'' for <mq q qbar> as in ',
     &        'Drees & Nojiri 1993 (=GLS)'
         write (*,*) ' value=''bg96'' for <mq q qbar> as in ',
     &        'Bergstrom & Gondolo 1996 (=GLS)'
         write (*,*) ' value=''pole'' to include squark pole in ',
     &        'scattering amplitude as in Drees & Nojiri 1993'
         write (*,*) ' value=''nopole'' to skip squark pole in ',
     &        'scattering amplitude as in Griest 1988 and Gondolo'
         write (*,*) ' value=''dn1'' to treat heavy quarks and gluons ',
     &        'as in Drees & Nojiri 1993'
         write (*,*) ' value=''best'' to use the best available ',
     &        'form factor'
         write (*,*) ' value=''no-ff'' to use no ',
     &        'form factor (i.e. q=0)'
         write (*,*) ' value=''L-S'' to use the Helm ',
     &        'form factor with best-fit Lewin-Smith parameters'
         write (*,*) ' value=''ds4.1'' to use the Helm ',
     &        'form factor as in DarkSUSY 4.1'
         write (*,*) ' value=''gauss'' to use the gaussian ',
     &        'form factor'
         write (*,*) ' value=''SOG'' to use the sum over gaussians ',
     &        'form factor'
         write (*,*) ' value=''FB'' to use the Fourier-Bessel ',
     &        'form factor'
         write (*,*) ' value=''Fermi'' to use the Helm ',
     &        'form factor with parameters from a Fermi function'
      endif
      if (sisd.eq.'sd'.and.cs.eq.'help') then
         goodoption=.true.
         write (*,*) 'dsddset: use "call dsddset(''sd'',value)" where'
         write (*,*) ' value=''smc'' for Delta q as measured in SMC'
         write (*,*) ' value=''emc'' for Delta q as measured in EMC'
         write (*,*) ' value=''dn93'' for Delta q as in ',
     &        'Drees & Nojiri 1993 (=EMC)'
         write (*,*) ' value=''bg96'' for Delta q as in ',
     &        'Bergstrom & Gondolo 1996 (=SMC)'
         write (*,*) ' value=''pole'' to include squark pole in ',
     &        'scattering amplitude as in Drees & Nojiri 1993'
         write (*,*) ' value=''nopole'' to skip squark pole in ',
     &        'scattering amplitude as in Griest 1988 and Gondolo'
         write (*,*) ' value=''default'' to use the default ',
     &        'spin form factors'
         write (*,*) ' value=''ISM'' to use various ',
     &        'interacting shell model spins'
         write (*,*) ' value=''ISMR'' to use the ',
     &        'interacting shell model spins for Ge by Ressell et al'
         write (*,*) ' value=''SPSM'' to use the single ',
     &        'particle shell model spins'
         write (*,*) ' value=''OddG'' to use the ',
     &        'odd group model spins'
         write (*,*) ' value=''no-ff'' to use the otherwise specified ',
     &        'spins but no form factor (i.e. q=0)'
         write (*,*) ' value=''gauss'' to use the otherwise specified ',
     &        'spins and a gaussian form factor'
      endif

c...spin-independent

c...gaisser,leutwyler,sainio,plb253(91)252
c...drees,nojiri,prd48(93)3483
c...bergstrom,gondolo,astrop.ph.
      if (sisd.eq.'si'.and.(
     &     cs.eq.'gls91'.or.cs.eq.'dn93'.or.cs.eq.'bg96'
     &     .or.cs.eq.'default')) then
         goodoption=.true.
         ftp(7)=0.023
         ftp(8)=0.034
         ftp(9)=0.0595
         ftp(10)=0.14
         ftp(11)=0.0595
         ftp(12)=0.0595
         ftn(7)=0.019
         ftn(8)=0.041
         ftn(9)=0.0592
         ftn(10)=0.14
         ftn(11)=0.0592
         ftn(12)=0.0592
         dddn = 0
         ddpole = 0
      endif
            
c...spin-dependent
c...drees,nojiri,prd48(93)3483
      if (sisd.eq.'sd'.and.(cs.eq.'dn93'.or.cs.eq.'emc')) then
         goodoption=.true.
         delu=0.77
         deld=-0.49
         dels=-0.15
      endif
            
c...spin-dependent
c...bergstrom,gondolo,astrop.ph.
c...d.adams et al cern-ppe/94-57
      if (sisd.eq.'sd'.and.(
     &     cs.eq.'bg96'.or.cs.eq.'smc'.or.cs.eq.'default')) then
         goodoption=.true.
         delu=0.74
         deld=-0.40
         dels=-0.12
      endif
      
c...allow for pole in scattering
      if (cs.eq.'pole') then
         goodoption=.true.
         ddpole=1
      endif
      if (cs.eq.'nopole'.or.cs.eq.'default') then
         goodoption=.true.
         ddpole=0
      endif

c...drees-nojiri treatment of heavy quarks
      if (sisd.eq.'si'.and.cs.eq.'dn1') then
         goodoption=.true.
         dddn=1
      endif
      if (sisd.eq.'si'.and.cs.eq.'default') then
         goodoption=.true.
         dddn=0
      endif

c...form factors
      if (sisd.eq.'si'.and.(cs.eq.'best'.or.cs.eq.'default')) then
         goodoption=.true.
         ddffsi = 'best'
      endif
      if (sisd.eq.'si'.and.(cs.eq.'no-ff')) then
         goodoption=.true.
         ddffsi = 'no-ff'
      endif
      if (sisd.eq.'si'.and.(cs.eq.'L-S'.or.
     &     cs.eq.'Helm')) then
         goodoption=.true.
         ddffsi = 'L-S'
      endif
      if (sisd.eq.'si'.and.(cs.eq.'ds4.1')) then
         goodoption=.true.
         ddffsi = 'ds4.1'
      endif
      if (sisd.eq.'si'.and.cs.eq.'gauss') then
         goodoption=.true.
         ddffsi = 'gauss'
      endif
      if (sisd.eq.'si'.and.cs.eq.'SOG') then
         goodoption = .true.
         ddffsi = 'SOG'
      endif
      if (sisd.eq.'si'.and.cs.eq.'FB') then
         goodoption = .true.
         ddffsi = 'FB'
      endif
      if (sisd.eq.'si'.and.cs.eq.'Fermi') then
         goodoption = .true.
         ddffsi = 'Fermi'
      endif

      if (sisd.eq.'sd'.and.(cs.eq.'default'.or.cs.eq.'best')) then
         goodoption=.true.
         ddffsd = ' best'
      endif
      if (sisd.eq.'sd'.and.cs.eq.'ISMR') then
         goodoption=.true.
         ddffsd = ddffsd(1:1)//'ISMR'
      endif
      if (sisd.eq.'sd'.and.cs.eq.'ISM') then
         goodoption=.true.
         ddffsd = ddffsd(1:1)//'ISM'
      endif
      if (sisd.eq.'sd'.and.cs.eq.'OddG') then
         goodoption=.true.
         ddffsd = ddffsd(1:1)//'OddG'
      endif
      if (sisd.eq.'sd'.and.cs.eq.'SPSM') then
         goodoption=.true.
         ddffsd = ddffsd(1:1)//'SPSM'
      endif
      if (sisd.eq.'sd'.and.cs.eq.'no-ff') then
         goodoption=.true.
         ddffsd(1:1)='1'
      endif
      if (sisd.eq.'sd'.and.cs.eq.'gauss') then
         goodoption=.true.
         ddffsd(1:1) = 'g'
      endif

c...invalid choice
      if (.not.goodoption) then
         write (*,*) 'dsddset: unrecognized options ',sisd,cs
         stop
      endif
      
      return
      end
