      program dsmain
c
c     This program is a sample DarkSUSY main program. It asks for model
c     parameters interactively and calculates the relic density and
c     various rates for each given model. This program can be used
c     as is for quick calculations, or be modified to your liking for
c     more advanced uses.
c
c     
c-----This line is 72 columns long--------------------------------------
c
      implicit none

      real*8 sumofwidths

      real*8 oh2,xf,dsrdomega                            ! relic density
      real*8 sigsip,sigsin,sigsdp,sigsdn                 ! direct detection
      integer i,nnuc,a(10),z(10),stoich(10)              ! direct detection
      real*8 e,si(10),sd(10),siff(10),sdff(10),t,rsi,rsd ! direct detection
      real*8 phiep,dsepdiff                              ! positron flux
      real*8 fluxgacdiff,fluxgac,fluxgaga,fluxgaz        ! gamma-ray flux
      real*8 jpsi,cospsi0,delta,fdelta,dshmjave,dshmj    ! gamma-ray flux
      real*8 nsigvgaga,nsigvgaz,nsigvgacont,nsigvgacdiff ! gamma-rays
      integer istat                                      ! gamma-ray flux, etc
      real*8 egam,dshrgacontdiff,egath,dshrgacont        ! gamma-ray flux
      real*8 tpbess(3),pb_a,pb_b,pb_c,dshrpbardiff       ! pbar flux
      real*8 eth,thmax,rateea,ratesu,energy,theta        ! neutrino telescopes
      integer rtype                                      ! neutrino telescopes
      real*8 gm2amp,dsgm2muon                            ! g-2 amplitude
      integer unphys,excl,hwarning,acceptable,iend,ierr,iwar,nfc
      integer opt
      real*8 dshrmuhalo,dshrmudiff,nmusigmav,phimuhalo,dsabsq,tp,phidb,
     &  dshrdbardiff,dsepspecm
      real*8 amu,am2,ama,atanbe,amsq,atm,abm,
     &  am0,amhf,aa0,asgnmu
      character message*80
      character slhafile*128,slhaout*128
      integer modchoice
      logical first
      data first/.true./

c     ----- ATLAS     
      character line*80
      character ATLASlha*128
      real*8 oh2with
      integer dooh2with
      character accelbit*128
      character acceltxt*128

c
c     Here we include the file dssusy.h which contains common block
c     definitions for particle masses, susy parameters, and results.
c     We also include dsidtag.h which gives the possibility to tag
c     every model with a 12 character id tag (variable idtag). This
c     id tag is printed in case of errors for that model.
c
      include 'dsmssm.h'    ! new
      include 'dsio.h'      ! new
c      include 'dssusy.h'   ! new
      include 'dsidtag.h'

c
c     This call initializes the DarkSUSY package. This call should
c     be the first call in any program using DarkSUSY. This call initializes
c     some global variables and calls various other modules to initialize
c     them with their default values. Check out src/ini/dsinit.f if you
c     want to see what it does.
c
      call dsinit

c
c     The amount of output onto standard output can be controlled by the
c     variable prtlevel. Setting prtlevel=0 suppresses all messages,
c     except fatal errors that terminate the program. Setting prtlevel=1
c     displays informational messages. Setting prtlevel to higher values
c     displays more and more information. In this example, we set
c     prtlevel=1 to have a message from dsreadpar. We reset it to 0
c     later. Had we set prtlevel=2 here, we would also have messages for
c     each variable which is (re)set in dsreadpar. Caution: if prtlevel
c     is reset in dsreadpar, the new value applies to subsequent
c     messages.
c
      prtlevel=1
c
c     This is an example of redefining DarkSUSY parameters reading them
c     from an input file, in this case the file dsdefault.dat which
c     contains the default values (this is done for test purposes only,
c     because the default values are already included in the code). The
c     user can copy the file dsdefault.dat to a new file, change the
c     values of the variables, and read in the new file.  Not all
c     variables need to be assigned in the new file, but only those
c     whose value differs from the default value. Note: dsdefault.dat
c     resets prtlevel to 0, ie no messages.
c     Note: the dsdefault.dat is not upto date, but a call like the
c     following can be used to reset parameters to your own choices.
c
c      open (unit=10,file='dsdefault.dat',status='old')
c      call dsreadpar(10)
c      close (10)

c
c     Before embarking on the scan in supersymemtric parameter space, we
c     calculate some quantities related to the model of the galactic halo
c     which do not depend on the specific supersymmetric model. 
c
c     First of all, we decide on the density profile of the dark
c     halo. We initialize the spherically symmetric density profile
c     provided with DarkSUSY (see dshmset.f for details).
c     The halo model is set with a call to dshmset, like below.
c     A call to dshmset('default') is done in dsinit, so you only have
c     to call dshmset if you want another profile than the default
c     Navarro Frenk and White profile. See dshmet.f for a list of available
c     standard options.
c    
c      call dshmset('isosm')  ! modified isothermal sphere


c
c     For the fluxes of gamma-rays and neutrinos with the chosen halo
c     profile, we calculate, once and for all, the line of sight
c     integration factor j in the direction of observation, which we
c     define as the direction which forms an angle psi0 with respect to
c     the direction of the galactic centre (e.g. cospsi0 = 1 is the
c     galactic center direction, cospsi0 = -1 the antigalactic centre).
c
      cospsi0=1.d0
c
c     There are two options: you can take into account that your 
c     detector has an angular resolution delta (in sr) and compute
c     j averaged over that solid angle. E.g. :
c
      delta=1.d-3
      write(*,*) 'Calculating jpsi, please be patient...'
      jpsi=dshmjave(cospsi0,delta)
      write(*,*) 'jpsi =',jpsi
c
c     or just take j in the direction psi0. The call would then be
c     (remove comments):
c     
c      jpsi=dshmj(cospsi0) 
c

c
c     Now we are ready to scan the supersymmetric parameter space. We say
c     that we define a supersymmetric 'model'. In this example program
c     we ask the user for the model parameters and then go on with the rest
c     of the calculations. Note, that this is only one way of assigning
c     parameters. For more heavy use, you most likely want to read model
c     parameters from a file or generate them randomly or on a grid. See
c     dstest.f for an example of reading the paramters from a file.
c     The user is free to assign parameters in the way
c     he/she prefers. Caution: all of the following supersymmetric
c     parameters must be assigned: the gaugino masses m1, m2, m3; the
c     Higgs pseudoscalar mass ma; the ratio of Higgs vacuum expectation
c     values tanbe; the square of the squark and slepton mass parameters
c     mass2q(i), mass2u(i), mass2d(i), mass2l(i), mass2e(i); the
c     trilinear soft parameters asofte(i), asoftu(i), asoftd(i). Here
c     i=1,2,3 is a generation index. In the current version of DarkSUSY,
c     all these parameters are at the weak scale. See dsgive_model for 
c     an example of how these are assigned for MSSM-7.
c

c     SOME ATLAS INITIALISATIONS
      dooh2with = 1
      oh2with = -1.
      unphys = 0
      hwarning = 0
      ATLASlha = 'ATLASlha'
      write(*,*) "INFO ATLAS: ATLASlha = ",ATLASlha
      !write(*,*) ATLASlha

c ----- ATLAS: read steer file (minimal, minimal, only to automatise)
      open (33,file='ATLASsteer',status='old')
      read(33,*) modchoice  ! easier, but structure below can be expanded
      read(33,*) dooh2with

      write(*,*) "INFO ATLAS: modchoice = ",modchoice
      write(*,*) "INFO ATLAS: dooh2with = ",dooh2with
      if(modchoice.eq.4) goto 1010 ! skip interactive step



c 10   read (unit,'(a)',end=99) line
c      if (line(1:1).eq.'!'.or.line(1:1).eq.'#'.or.
c     &  len(line).eq.0) goto 10
c      if (line(1:1) == '4') then 
c        modchoice = 4
c        goto 1010  !skipping interactive step
c      endif

c ----- 


 1000 write(*,*) ' '
      write(*,*) 'What kind of model do you want to look at?'
      write(*,*) '   0 = quit'
      write(*,*) '   1 = MSSM-7'
      write(*,*) '   2 = mSUGRA'
      write(*,*) '   3 = as read from an SLHA2 file'
      read(5,*) modchoice
      if (modchoice.eq.0) goto 2000

 1010 continue
      
c...MSSM
      if (modchoice.eq.1) then ! MSSM-7
         write(*,*) 'Enter mu (GeV): '
         read(5,*) amu
         write(*,*) 'Enter M2 (GeV): '
         read(5,*) am2
         write(*,*) 'Enter MA (mass of CP-odd Higgs) (GeV): '
         read(5,*) ama
         write(*,*) 'Enter tan(beta): '
         read(5,*) atanbe
         write(*,*) 'Enter m0 (common sfermion mass scale) (GeV): '
         read(5,*) amsq
         write(*,*) 'Enter At/m0 [-3:3]: '
         read(5,*) atm
         write(*,*) 'Enter Ab/m0 [-3:3]: '
         read(5,*) abm

c...Now transfer these variables to the DarkSUSY common blocks, properly
c...setting up all relevant low-energy parameters         
         call dsgive_model(amu,am2,ama,atanbe,amsq,atm,abm)

c...For larger runs, it is also convenient to define a unique 12-character
c...idtag for each model. This one will be printed by DarkSUSY if there is
c...any problem with a specific model. For now, just set it to something.

         idtag='INTMOD000001'

c...mSUGRA
      elseif (modchoice.eq.2) then

         write(*,*) 'Enter m0 (scalar mass parameter @ GUT) (GeV): '
         read(5,*) am0
         write(*,*) 'Enter m_1/2 (gaugino mass param @ GUT) (GeV): '
         read(5,*) amhf
         write(*,*) 'Enter A0 (trilinear coupling @ GUT) (GeV): '
         read(5,*) aa0
         write(*,*) 'Enter sign of mu (+1 or -1): '
         read(5,*) asgnmu
         write(*,*) 'Enter tan(beta): '
         read(5,*) atanbe

c...Now transfer these parameters to the DarkSUSY common blocks, preparing
c...for RGE running
         call dsgive_model_isasugra(am0,amhf,aa0,asgnmu,atanbe)
      

      elseif (modchoice.eq.3) then
         write(*,*) 'Enter the filename for your SLHA2 file: '
         read(5,'(A)') slhafile
         call dsSLHAread(slhafile,0) ! 0=no warnings, 1=print warnings

      elseif (modchoice.eq.4) then
         write(*,*) 'Reading default SLHA2 file: ',ATLASlha
         !read(slhafile,'(A)') ATLASlha
         !slhafile = 'ATLASlha'
         !call dsSLHAread(slhafile,0) ! 0=no warnings, 1=print warnings
         call dsSLHAread(ATLASlha,0) ! 0=no warnings, 1=print warnings

      else
         write(*,*) 'Not a valid choice of model: ',modchoice
         stop
      endif

      write(*,*) ' '
      write(*,*) '***** MODEL: ',idtag,' *****'

c
c     Now we are ready to calculate the supersymmetric particle spectrum
c     and the 3-particle vertices. We use the routine dssusy. In output,
c     the arguments of dssusy are non-zero if one of the flags unphys or
c     hwarning is set. The flag unphys<0 means that the model is
c     theoretically inconsistent, unphys=1 means that the neutralino is
c     not the LSP. Check the meaning of unphys by calling
c     dswunph(u,unphys) where u is an I/O unit number. The flag hwarning
c     is non-zero when some of the approximations in the radiative
c     corrections in the Higgs sector are applied outside their range of
c     validity. Check the meaning of hwarning by calling
c     dswhwarn(u,hwarning) where u is an I/O unit number.
c
c
      if (modchoice.eq.1) then
        call dssusy(unphys,hwarning)  ! MSSM-7 or any EW scale MSSM model

      elseif (modchoice.eq.2) then
        call dssusy_isasugra(unphys,hwarning) ! mSUGRA

      elseif (modchoice.eq.3.or.modchoice.eq.4) then
c...Note: For SLHA files, we should NOT call dssusy as we would then
c...overwrite things read from the SLHA file. Hence, call dsprep directly.
        call dsprep
        open (unit=30,file='dsmain1.tmp')
        write(30,*) ' '
        write(30,*) '***** MODEL: ',idtag,' *****'
        call dswspectrum(30)
        call dswvertx(30)
        close (30)
      endif


c
c     Then we check if the model satisfies current accelerator limits.
c     We call the routine dsacbnd, which returns the variable excl.
c     This variable excl is non-zero when an accelerator limit is
c     violated. Use a call to dswexcl(u,excl), where u is an I/O unit
c     number, to find which limit(s) is(are) violated. Which set of
c     accelerator bounds to use is fixed by a call to dsacset(aclabel)
c     where aclabel is 'pdg2000', 'pdg1999', etc. The default is
c     'pdg2002c' and is set by a call to dsacset('default'). This call is
c     done automatically in dsinit, so you only need to call dsacset if
c     you want another set of accelerator bounds. Check dsacset.f for
c     the available set of options. If you want to implement your own
c     accelerator bounds, you can of course replace this call with one
c     to your own routines instead.
c
      if (unphys.ne.0.or.hwarning.ne.0) then
         acceptable=-1
         excl=0
         !write(*,*)'here\n'
         !if (unphys.lt.0) then
         !  write(30,*) '---------- more on unphys'
         !  call dswunph(30,unphys)
         !endif

         goto 150  !FIX
      endif

      call dsacbnd(excl)
c
c     As a test, we check if the flags are non-zero, in which case we
c     output the reasons for the theoretical inconsistency or the
c     experimental exclusion.
c
      if (excl.eq.0.and.unphys.eq.0.and.hwarning.eq.0) then
         acceptable=0
      else
         acceptable=-1
      endif

 150  write (message,*) 'acceptable =',acceptable, ' (0=OK, -1=not OK)'
      call dswrite(0,0,message)
c
c     If the model is not acceptable, call dswunph and dswexcl to print
c     out why. Also, if warnings were issued for the Higgs calculation,
c     print out what the warnings were.
c
      if (acceptable.ne.0) then
         call dswunph(6,unphys)
         call dswexcl(6,excl)
         !BKGgoto 1000 ! next model
      endif
      if (hwarning.ne.0) then
         call dswhwarn(6,hwarning)
      endif
c
c     As a test, for the first model we write the particle spectrum and
c     the 3-particle vertices into the file dstest1.tmp.
c     Uncomment the following lines if you want this
c
c      if (first) then
c        open (unit=30,file='dstest1.tmp')
c        write(30,*) ' '
c        write(30,*) '***** MODEL: ',idtag,' *****'
c        call dswspectrum(30)
c        call dswvertx(30)
c        close (30)
c        write (message,*) 'Spectrum and vertices written to dstest1.tmp'
c        call dswrite(0,0,message)
cc     We now set the print level to 0 to get severe error messages only.
c        prtlevel=0
c        first=.false.
c      endif

c
c     We now have all the masses and couplings calculated. Let's print
c     some out to see what they are. dsabsq is just the square of the
c     absolute value of it's complex*16 argument.

      write(*,*) '  Neutralino mass = ',mass(kn(1))
      write(*,*) '  Gaugino fraction = ',
     &  dsabsq(neunmx(1,1))+dsabsq(neunmx(1,2))
      write(*,*) '  H1 mass =  ',mass(kh1),width(kh1)
      write(*,*) '  H2 mass =  ',mass(kh2),width(kh2)
      write(*,*) '  H3 mass =  ',mass(kh3),width(kh3)
      write(*,*) '  H+- mass = ',mass(khc),width(khc)

c
c     We are now ready to calculate the relic density and different
c     direct or indirect detection rates. We will start with the relic
c     density. The relic density is calculated by the function dsrdomega.
c     The first argument determines if coannihilations should be included.
c     If the argument is 0, no coannihilations are included, if it is 1,
c     all relevant coannihilations are included (charginos, neutralinos
c     and sfermions if light enough). One can also choose to only
c     include a partial set of coannihilations: argument=2, only chargino
c     and neutralino coannihilaitons, argument=3, only sfermion 
c     coannihiltions.
c     The second argument determines the
c     accuracy, if set to 1 a faster calculation (e.g. only including
c     coannihilations if the mass difference is less then 30%) is
c     performed and if set to 0 a more accurate calculation is performed.
c     In practice, the fast option is more than adequate (better than
c     5% accuracy). The function returns the relic density, omega h^2,
c     the freeze-out temperature, xf (x=m/T), one error flag, ierr and
c     one warning flag, iwar, both of which are 0 if everything is OK.
c     nfc is the number points in momentum where the cross section had
c     to be calculated. Note that the omega calculation can be very
c     time consuming, needing up to several minutes for tricky models
c     (with many resonances, thresholds, coannihilations etc).
c     
c     The next line selects the effective degrees of freedom in the early
c     Universe. The default is Hindmarsh & Philipsen, table B. If you want
c     something else, e.g. Gondolo and Gelmini 1991 (150 MeV) d.o.f.,
c     uncomment the call to dsrdset below.
c     Use call dsrdset('dof','help') to print out the other possibilities.
c     For example, '3' gives Hindmarsch & Philipsen's table B.

c      call dsrdset('dof','1') ! choose GG (150 MeV) dof instead of HP-B default

      write(*,*) 'Calculating omega h^2 without coannihilations,',
     &     ' please be patient...'
      oh2=dsrdomega(0,1,xf,ierr,iwar,nfc)
      write(*,*) '  without coannihilations Oh2 = ',oh2,ierr,iwar
      write(*,*) 'Calculating omega h^2 with coannihilations,',
     &     ' please be patient...'
      if(dooh2with.eq.1) then 
        oh2with=dsrdomega(1,1,xf,ierr,iwar,nfc)
        write(*,*) '  with coannihilations Oh2 = ',oh2with,ierr,iwar
      else
        write(*,*) '  skipwith coannihilations Oh2 = ',-1.,-1,-1
      endif

c
c     Now we continue to calculate different direct and indirect
c     detection rates. One word of caution applies here regarding the
c     local halo density. If the relic density is so small that the
c     neutralinos cannot contribute 100% to the local halo density, we
c     have to rescale the local halo density of neutralinos with the
c     relic density. In most cases, this would be simple to do later by
c     the user, except for the neutrino telescope routines, where non
c     -linear equations enter. To be consistent, we have thus chosen to
c     rescale all the rates if the relic density is too low. However,
c     it is up to the user to make sure this rescaling is done by calling
c     dshmrescale_rho(oh2,oh2min), where the first argument is the
c     relic density and the second argument is the minimum relic
c     density required for the neutralinos to contribute 100% to the
c     local halo density. We often choose this to 0.025. The rescaled
c     local halo density is then given by rhox=rho0*(oh2/oh2min) if
c     oh2<oh2min and rhox=rho0 otherwise. Note however, that only rates
c     are rescaled, if you calculate cross sections, you get the cross
c     sections without any rescaling.

      call dshmrescale_rho(oh2,0.025d0)

c
c     We are now ready to calculate rates.  Let's start with scattering
c     cross sections for direct detection experiments by calling
c     dsddneunuc which returns the spin-independent and the spin-dependent
c     scattering cross sections off protons and neutrons. The cross
c     sections are returned in units of cm^2.
c

      write (*,*) 'Calculating scattering cross sections...'
      call dsddneunuc(sigsip,sigsin,sigsdp,sigsdn)
      write(*,*) '  sigsip (pb) = ',sigsip*1.0d36
      write(*,*) '  sigsin (pb) = ',sigsin*1.0d36
      write(*,*) '  sigsdp (pb) = ',sigsdp*1.0d36
      write(*,*) '  sigsdn (pb) = ',sigsdn*1.0d36
      call dsddsigma(1,1,1,sigsip,sigsdp)
      write(*,*) ' proton: sigsi (pb) = ',sigsip*1.0d36,
     &     ' sigsd (pb) = ',sigsdp*1.0d36


c
c     Now take a more complex example with scattering off a compound
c
      nnuc = 5
      a(1) = 23                 ! Na-23
      z(1) = 11                 ! Na-23
      a(2) = 127                ! I-127
      z(2) = 53                 ! I-127
      a(3) = 73                 ! Ge-73
      z(3) = 32                 ! Ge-73
      a(4) = 1                  ! p
      z(4) = 1                  ! p
      a(5) = 1                  ! n
      z(5) = 0                  ! n
      call dsddsigma(nnuc,a,z,si,sd)
      e = 10.0 ! recoil energy in keV
      call dsddsigmaff(e,nnuc,a,z,siff,sdff)
      do i=1,nnuc
         write(*,*) '  A=',a(i),' Z=',z(i),
     &        ' sigsi (pb) = ',si(i)*1.0d36,
     &        ' sigsd (pb) = ',sd(i)*1.0d36
         write(*,*) ' sigsi*ff(pb)=',
     &        siff(i)*1.0d36,' sigsd*ff(pb)=',sdff(i)*1.0d36
      enddo
      stoich(1)=1
      stoich(2)=1
      t = 651.3d0
      call dsdddrde(t,e,2,a,z,stoich,rsi,rsd,1)
      write(*,*) ' NaI :  dRdE [counts/kg-day-keV] (SI)=',rsi,
     &     ' (SD)=',rsd
      
c
c     Next we compute cosmic-ray rates from neutralino annihilations 
c     in the galactic halo. There are several cosmic ray species.
c
c
c     The gamma-ray flux contains a factor which depends just on susy
c     parameters (these are halo independent quantities which one may
c     want to save as well). We calculate it here and then multiply by
c     the line-of-sight integral j computed before scanning in parameter
c     space.  There are the following options:
c
c     1) gamma-ray flux with continuum energy spectrum at a given energy
c     egam (GeV). E.g. :
c

      write (*,*) 'Calculating gamma ray fluxes...'
      egam=20.d0
      fluxgacdiff=dshrgacontdiff(egam,jpsi*delta,istat) !ph cm^-2 s^-1 GeV^-1
c
c     if jpsi was computed without averaging over the angular 
c     acceptance of the detector, one should have:
c
c      fluxgacdiff=dshrgacontdiff(egam,jpsi,istat) !ph cm^-2 s^-1 sr^-1 GeV^-1

c
c     2) gamma-ray flux with continuum energy spectrum integrated above
c        some given threshold egath (GeV). E.g. :
c

      egath=1.d0
      fluxgac=dshrgacont(egath,jpsi*delta,istat) !ph cm^-2 s^-1 
c
c     if jpsi was computed without averaging over the angular 
c     acceptance of the detector, one should have:
c
c      fluxgac=dshrgacont(egam,jpsi,istat) !ph cm^-2 s^-1 sr^-1

c
c     3) monochromatic gamma-ray flux induced by 1-loop annihilation
c        processes into a 2-body final state containing a photon. There
c        are two such final states: the 2 photon final state (with 
c        energy of the photons equal to the neutralino mass mx) and
c        the final state with a photon and a Z boson (with energy of 
c        the photon equal to mx * (1 - mz**2 / (4 * mx**2)) and mz
c        the mass of the Z boson)
c

      call dshrgaline(jpsi*delta,fluxgaga,fluxgaz) !ph cm^-2 s^-1 

c
c     if jpsi was computed without averaging over the angular 
c     acceptance of the detector, one should have:
c
c      call dshrgaline(jpsi,fluxgaga,fluxgaz) !ph cm^-2 s^-1 sr^-1

c
c     If we want to do the gymnastics ourselves converting cross
c     sections to fluxes, we can also call dsnsigvgacont,
c     dsnsigvgacdiff and dsnsigvgaline to get the number of photons 
c     (2 for gamma gamma, 1 for Z gamma and whatever the number is for
c     continuous gammas) times the cross section. 
c     The number returned is the dimensionless number
c        N_gammas * (sigma v) / (10^-29 cm^3 s^-1)
c     with (sigma v) in unit of cm^3 s^-1. For dsnsigvgacdiff the units
c     are of course GeV^-1 instead.

      egam=1.0d0
      call dsnsigvgacont(egam,nsigvgacont)     ! dimensionless
      call dsnsigvgacdiff(egam,nsigvgacdiff)   ! GeV^-1
      call dsnsigvgaline(nsigvgaga,nsigvgaz)   ! dimensionless
      
c
c     had one picked hclumpy = 2, the rescaling factor fdelta is
c     still missing. One may have some theoretical estimate for it
c     or assume some value and then check that constraints are not 
c     violated. E.g.:
c
c      fdelta=100.d0
c      fluxgacdiff=fluxgacdiff*fdelta
c      fluxgac=fluxgac*fdelta
c      fluxgaga=fluxgaga*fdelta
c      fluxgaz=fluxgaz*fdelta
c

      write(*,*) '  fluxgacdiff = ',fluxgacdiff,' ph/(cm^2 s GeV)'
      write(*,*) '      fluxgac = ',fluxgac,' ph/(cm^2 s)'
      write(*,*) '     fluxgaga = ',fluxgaga,' ph/(cm^2 s)'
      write(*,*) '      fluxgaz = ',fluxgaz,' ph/(cm^2 s)'
      write(*,*) '  nsigvgacont = ',nsigvgacont
      write(*,*) ' nsigvgacdiff = ',nsigvgacdiff,' GeV^-1'
      write(*,*) '    nsigvgaga = ',nsigvgaga
      write(*,*) '     nsigvgaz = ',nsigvgaz

c
c     Antiproton rates, interstellar, solarmodulated
c
c     Now we come to the antiproton routines. We get the differential
c     spectrum in units of GeV^-1 cm^-2 s^-1 sr^-1 by calling
c     dshrpbardiff with the antiproton kinetic energy as the first
c     argument. The second argument determines if the flux is solar
c     modulated (a la Perko) or not (0=no solar modulation, 1=solar
c     modulation). The third argument determines how the diffusion
c     equations are solved. If the argument is
c          1 = they are solved for the requested energy only
c          2 = they are tabulated for a large range of energies for the first
c              call and the table is used for subsequent calls.
c          3 = as 2, but the table is also written to disk in the current
c              directory on the first call.
c          4 = the table is read from disk on first call, and is used for
c              that and subsequent calls.
c     The tabulation can take several minutes (or hours), 
c     so if you are only interested in a few models, don't tabulte,
c     but if you want to calculate many models, do tabulate 
c     (i.e. using any of the options 2-4). Some standard tables are 
c     included in the distribution and are available in files
c     of the form dat/pbtd-*.dat. You can change the default propagation
c     model with a call to dspbset (see that routine for details).



      write (*,*) 'Calculating antiproton fluxes...'
      tpbess(1)=0.35d0
      tpbess(2)=1.76d0
      tpbess(3)=3.00d0
         


c     Here are two example calls, the first one uses the routines directly,
c     whereas the second uses tables. Tables are slower to startup, especially
c     if the table needs to be recreacted (could take hours), but is
c     faster after that initial slowdown. However, if you use tables and
c     change properties of the diffusion model that requires a retabulation
c     it is your responsibility to retabulate (easiest done by deleting
c     the appropriate tabulated file share/DarkSUSY/pbtd-*.dat).

c      pb_a=dshrpbardiff(tpbess(1),1,1) ! no tables, slower
      pb_a=dshrpbardiff(tpbess(1),1,4) ! use tables, faster after slower startup
      write(*,*) '  solar modulated pbar flux at 0.35 GeV = ',pb_a,
     &  ' GeV^-1 cm^-2 s^-1 sr^-1'
      pb_b=dshrpbardiff(tpbess(2),1,4)
      write(*,*) '  solar modulated pbar flux at 1.76 GeV = ',pb_b,
     &  ' GeV^-1 cm^-2 s^-1 sr^-1'
      pb_c=dshrpbardiff(tpbess(3),1,4)
      write(*,*) '  solar modulated pbar flux at 3.00 GeV = ',pb_c,
     &  ' GeV^-1 cm^-2 s^-1 sr^-1'


c
c     had one picked hclumpy = 2, the rescaling factor fdelta is
c     still missing. One may have some theoretical estimate for it
c     or assume some value and then check that constraints are not 
c     violated. E.g.:
c
c      fdelta=100.d0
c      pb_a=pb_a*fdelta
c      pb_b=pb_b*fdelta
c      pb_c=pb_c*fdelta
c

c
c      We can also calculate the flux of antideutrons. The call is very
c      similar to that for pbar. The dbar's use the same diffusion model
c      as the pbar and can thus be changed with a call to dspbset.
c

      write (*,*) 'Calculating antideutron fluxes...'
      tp=1.0d0

      phidb=dshrdbardiff(tp,1,4)   ! use tables, faster after slower startup
c      phidb=dshrdbardiff(tp,1,1)   ! no tables, slower
      write(*,*) '  solar modulated dbar flux at 1.00 GeV = ',phidb,
     &  ' GeV^-1 cm^-2 s^-1 sr^-1'

c
c     We now continue with the rates of positrons from neutralino
c     halo annihilation. With the default propagation model 'esu04'
c     (Baltz and Edsjo propgation with new parameters), you get
c     the differential flux with a call like below. The first argument
c     is the positron energy (GeV) and the second argument deterimes
c     if the diffusion model should be read from disk or created on the
c     spot. Default is to read from disk, and if the file does not exist,
c     it will recreate it. Some standard files are available in 
c     share/DarkSUSY/eptab-*.dat and share/DarkSUSY/epgretab-*.dat.
c

      write (*,*) 'Calculating positron fluxes at 1 GeV...'
c      phiep=dsepspecm(1.0d0,1) ! no tables, slower
c      phiep=dsepspecm(1.0d0,4) ! use tables
c      write (*,*) '  phiep=',phiep,
c     &  ' GeV^-1 cm^-2 s^-1 sr^-1'

c...The routine below is an alternative routine (completely independent
c...implementation, but with different approximations).
      phiep=dsepdiff(1.0d0,1,4) ! use tables
      write (*,*) '  phiep=',phiep,
     &  ' GeV^-1 cm^-2 s^-1 sr^-1'



c     Analogously with the antiprotons, you can change the diffusion model
c     with a call to dsepset. Note, however, if you choose any other model
c     then epdiffmethod=1 (see dsepset), you should call dsepdiff(e,2,4)
c     instead of dsepspecm.



c
c     Neutrino telescope rates
c
c     To calculate the rates in neutrino telescopes we call dsntrates. We
c     can either calculate the flux of neutrinos, the conversion rate
c     (per volume element) of neutrinos to muons or the muon flux. This
c     is determined by the argument rtype. We can also choose the
c     energy threshold and the maximal half-aperture angle from the
c     center of the Earth/Sun we are interested in. The rates for both
c     the earth and the sun are returned in units of km^-2 yr^-1 for
c     rtype=1,3 and in units of km^-3 yr^-1 for rtype=2. If some
c     warnings were issued, the flag istat is non-zero.
c
c     The defualt calculation method is to use the full expressions by
c     Gould and numerically integrate them over the velocity distribution
c     as specified by the halo profile, e.g. a gaussian for an
c     isothermal sphere. This numerical integration is rather slow and
c     there is thus an option to use tables (read from disk, recreated
c     if absent) instead. This is the default. Changing to numerical 
c     integration directly is handled by a call to dsntset (see dsntset.f
c     for details). Also other (approximate) formulae or the Damour
c     Krauss population are available and can be chosen by a call to dsntset.
c     Note: For the earth, which captures from a population of WIMPs bound
c     in the solar system, a new estimate (Lundberg and Edsjo, 
c     astro-ph/0401113) is used as a default.
c
c     

      write(*,*) 'Calculating rates in neutrino telescopes'
      eth=1.0d0      ! energy threshold (of neutrino/muon), GeV
      thmax=30.0d0   ! the maximum half-aperture angle, degrees
      rtype=3        ! 1=neutrino flux
                     ! 2=neutrino to muon conversion rate
                     ! 3=muon flux
      call dsntrates(eth,thmax,rtype,rateea,ratesu,istat)

      write(*,*) '  Flux from the Earth = ',rateea,
     &  ' km^-2 yr^-1'
      write(*,*) '  Flux from the Sun =   ',ratesu,
     &  ' km^-2 yr^-1'

c     If you want differential rates, you instead call dsntdiffrates. The
c     units of the rates are then km^-2 yr^-1 GeV^-1 degree^-1 for
c     rtype=1 or rtype=3 and km^-3 yr^-1 GeV^-1 degree^-1 for
c     rtype=2. Uncomment the following lines if you want differential
c     rates.

c      energy=10.0d0  ! energy (of neutrino/muon), GeV
c      theta=30.0d0   ! angle from center of Earth/Sun, degrees
c      rtype=3        ! 1=neutrino flux
c                     ! 2=neutrino to muon conversion rate
c                     ! 3=muon flux
c      call dsntdiffrates(energy,theta,rtype,rateea,ratesu,istat)

c
c     Muon rates from the halo
c     
c     It is also possible to calculate the neutrino-induced muon rates
c     in neutrino telescopes that would occur from neutralino
c     annihilations in our galactic halo. The neutrino-induced muon
c     flux above a threshold Eth is given by
c     

      write(*,*) 
     &  'Calculating neutrino-induced muon fluxes from the halo...'
      eth=1.0d0
      phimuhalo=dshrmuhalo(jpsi*delta,eth,nmusigmav,istat)
          ! km^-2 yr^-1

c
c     If we instead wanted the differential flux, we would call
c     (uncomment if you want this)
c
c      phimuhalo=dshrmudiff(jpsi*delta,eth,nmusigmav,istat)
c          ! km^-2 yr^-1 GeV^-1

      write(*,*) '  Muon flux from halo = ',phimuhalo,
     &  ' km^-2 yr^-1'

c
c     Now we calculate the MSSM contribution to the g-2 amplitude
c
      gm2amp=dsgm2muon()
      write(*,*) 'g-2 amplitude [a_mu = (g-2)/2 = ]: ',gm2amp


      write(*,*) 
     &  'Do you want to write out an SLHA2 file for your model?'
      write(*,*) '  0 = no'
      write(*,*) '  1 = yes, with full 6x6 sfermion mixing'
      write(*,*) '  2 = yes, with minimal flavour violation'

      if(modchoice.eq.4) then
        opt = 1
        slhaout = 'ATLASout_6x6.slha2'
        call dsSLHAwrite('ATLASout_6x6.slha2', 1)
        call dsSLHAwrite('ATLASout.slha2', 2)
      else 
        read(*,*) opt
        if(opt.eq.0) goto 1000
        write(*,*) 'Give SLHA2 file name:'
        read(5,'(A)') slhaout
        call dsSLHAwrite(slhaout,opt) 
      endif
      
c
c     end of loop
c
      if(modchoice.ne.4) goto 1000  ! BKG

 2000 continue
      close (10)
        close (30)
c     
c     Last lines of the main test program.
c
      write (*,*) 'The DarkSUSY example program is now finished.'
      write (*,*) '----------------------------------------------'


      ! oh2
      ! oh2with
      ! gm2amp 

      ! details in log: higgs, 
      ! details in log: Tevatron, LEP

      !write(*,*) accelbit

      write(*,*) 'Non-format: ',acceptable, unphys, excl

      print 18, acceptable, unphys, excl, oh2, oh2with, gm2amp

 18   format('ATLAS_res0 | tot: ',I2,' | unphys: ',I2,' | accel: ',I3,
     & ' | cdm= ',F6.3,1x,F6.3,
     & ' | (g-2)/2= ',E10.3)

      write(*,*)'sum_excl: ',excl
      if(btest(excl,0)) write(*,*)'ATLAS_excl: 0  chargino ',mass(kcha1)
      if(btest(excl,1)) write(*,*)'ATLAS_excl: 1  gluino ',mass(kgluin)
      if(btest(excl,2)) write(*,*)'ATLAS_excl: 2  squark ',mass(ksqd(1))
    !     &,mass(ksqu1),mass(ksqu1)
      if(btest(excl,3)) write(*,*)'ATLAS_excl: 3  slepton '
      if(btest(excl,4)) write(*,*)'ATLAS_excl: 4  gamma_z(inv) '
      if(btest(excl,5)) write(*,*)'ATLAS_excl: 5  h2 ',mass(kh2)
      if(btest(excl,6)) write(*,*)'ATLAS_excl: 6  neutralino ',mass(kn1)
      if(btest(excl,7)) write(*,*)'ATLAS_excl: 7  b->s+gamma '
      if(btest(excl,8)) write(*,*)'ATLAS_excl: 8  delta_rho '

      ! shows: dL,sL,b1, dR,sR,b2
      write(*,*)'TEST sq(d):'
     &     ,mass(ksqd(1)), mass(ksqd(2)), mass(ksqd(3))
     &     ,mass(ksqd(4)), mass(ksqd(5)), mass(ksqd(6))

      ! shows: uL,cL,t1, uR,cR,t2
      write(*,*)'TEST sq(u):'
     &     ,mass(ksqu(1)), mass(ksqu(2)), mass(ksqu(3))
     &     ,mass(ksqu(4)), mass(ksqu(5)), mass(ksqu(6))



      stop
      end


c
c     Subroutines that are not used by this example program, but that 
c     may prove useful when reading/writing larger sets of models
c


      subroutine read_model(lunit,nmodel,iend,ierr)
c
c     To read in a model from a file
c
c     NOTE: We here show how you can read in a model file and define all
c     DarkSUSY model parameters. For MSSM models, there exists a routine
c     src/su/dsgive_model.f that makes the definitions for you. 
c     We here read the model parameters and call that routine.

      implicit none
c      include 'dssusy.h'  ! new
      include 'dsidtag.h'
      include 'dsmssm.h'   ! new
      integer nmodel,lunit,iend,ierr
      real*8 at,ab,mqtild
      integer i
      character*40 message
 2000 format (1x,a12,7(1x,e14.8))
      ierr=0
c... When nmodel<0, skip -nmodel lines
      if (nmodel.lt.0) then
         do i=1,-nmodel
            read (lunit,*,end=1000,err=3000)
         enddo
         return
      endif
c... If nmodel>0, read n-th model (assumes header line)
      if (nmodel.gt.0) then
         do i=1,nmodel
            read (lunit,*,end=1000,err=3000)
         enddo
      endif
c... If nmodel=0, read next model
      read (lunit,2000,end=1000,err=3000) 
     &     idtag,mu,m2,ma,tanbe,mqtild,at,ab
c... modify/set additional parameters
c      higloop=5  ! 5 = Full FeynHiggs;  6 = FeynHiggsFast
      ! W mass for unitarity of tree-level annihilation amplitudes
      call dsgive_model(mu,m2,ma,tanbe,mqtild,at,ab)
      return

 1000 continue
      iend=1
      ierr=1
      write (message,*) 'End of model file (unit=',lunit,')'
      call dswrite(1,0,message)
 3000 continue
      ierr=1
      stop
      end



      subroutine random_model()
c
c     To generate model parameters in a random way
c
      implicit none
c      include 'dssusy.h'  ! new
      include 'dsmssm.h'   ! new
      include 'dsidtag.h'
      real*8 dsrndlog,dsrndlin,dsrndsgn
      integer first,n,idum,i
      real*8 msq,atm,abm
      real*8 mumin,mumax,m1min,m1max,m2min,m2max,mamin,mamax,
     &     tbmin,tbmax,msqmin,msqmax,atmmin,atmmax,abmmin,abmmax
      real*8 dsgf2s2thw
      data first/0/, n/0/, idum/0/
      save first,n,idum
c... First time set random number seed
      if (first.eq.0) then
         idum = -3782932
         n=0
         first=1
      endif
c...  Limits of the Mu-M2-MA-tan(beta)-Msq-At-Ab region explored
      mumin = 10.d0
      mumax = 10000.d0
      m2min = 10.d0
      m2max = 10000.d0
      m1min = 10.d0
      m1max = 10000.d0
      mamin = 10.d0
      mamax = 1000.d0
      tbmin = 1.001d0
      tbmax = 60.d0
      msqmin = 50.d0
      msqmax = 1000.d0
      atmmin = -3.d0
      atmmax = 3.d0
      abmmin = -3.d0
      abmmax = 3.d0
c... random point
      n=n+1
      write (idtag,'(a,i8.8)') 'rndm',n
      mu = dsrndsgn(idum) * dsrndlog(idum,abs(mumin),abs(mumax))
      m2 = dsrndsgn(idum) * dsrndlog(idum,abs(m2min),abs(m2max))
      s2wmz=dsgf2s2thw(GFermi,alphem,mass(kz),mass(kt),3)
      m1 = m2 * 5.0d0/3.0d0 * s2wmz/(1.d0-s2wmz)
      ma = dsrndlog(idum,mamin,mamax)
      tanbe = dsrndlog(idum,tbmin,tbmax)
      msq = dsrndlog(idum,msqmin,msqmax)
      atm = dsrndlin(idum,atmmin,atmmax)
      abm = dsrndlin(idum,abmmin,abmmax)
c... further assumptions
      ! W mass for unitarity of tree-level annihilation amplitudes
      mass(kw)=mass(kz)*sqrt(1.d0-s2thw)
      ! GUT relation on gaugino mass parameters
      m3 =  m2 * (alph3mz*s2wmz)/alphem    ! alph3->alph3mz 020912 (JE)
      ! Slepton and squark mass matrices at the weak scale
      mass2q(3) = msq**2
      mass2u(3) = msq**2
      mass2d(3) = 2.0d0*mass2q(3)-mass2u(3)
      mass2l(3) = mass2d(3)
      mass2e(3) = mass2d(3)
      asofte(3) = 0.d0
      asoftu(3) = atm*msq
      asoftd(3) = abm*msq
      do i=1,2
         mass2q(i) = mass2d(3)
         mass2u(i) = mass2d(3)
         mass2d(i) = mass2d(3)
         mass2l(i) = mass2d(3)
         mass2e(i) = mass2d(3)
         asofte(i) = 0.d0
         asoftu(i) = 0.d0
         asoftd(i) = 0.d0
      enddo
      return
      end


      subroutine write_model(lunit)
c
c     To write model parameters to unit lunit
c
      implicit none
c      include 'dssusy.h'   ! new
      include 'dsidtag.h'
      include 'dsmssm.h'    ! new
      real*8 at,ab,mtop,mqtild
      integer lunit
      integer first
      save first
      data first/0/
      if (first.eq.0) then
         write (lunit,1001)
         first = 1
      endif
      mtop = mass(kt)
      mqtild=dsqrt(mass2q(3))
      at = asoftu(3)/mqtild
      ab = asoftd(3)/mqtild
      write (lunit,2001) 
     &     idtag,mu,m2,ma,tanbe,mqtild,at,ab
      return
 1001 format ('#','.....id.....',1x,'......mu......',1x,
     &     '......m2......',1x,'......ma......',1x,
     &     '.....tanbe....',1x,'......m.......',1x,
     &     '.....at/m.....',1x,'.....ab/m.....')
 2001 format (1x,a12,7(1x,e14.8))
      end
