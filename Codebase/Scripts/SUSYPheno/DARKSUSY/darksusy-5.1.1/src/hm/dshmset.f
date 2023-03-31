      subroutine dshmset(c)
****************************************************************
*** subroutine dshmset:                                      ***
*** initialize the density profile and or the small clump    ***
*** probability distribution                                 ***
*** type of halo:                                            ***
***   hclumpy=1 smooth, hclumpy=2 clumpy                     ***
***                                                          ***
*** a few sample cases are given; specified values of the    ***
*** local halo density 'rho0' and of the length scale        ***
*** parameter 'a' should be considered just indicative       ***
***                                                          ***
*** author: piero ullio (piero@tapir.caltech.edu)            ***
*** date: 00-07-13                                           ***
*** small modif: paolo gondolo 00-07-19                      ***
*** mod: 03-11-19 je, 04-01-13 pu, 09-05-07 ps, 09-08-08 ps  ***
****************************************************************
      implicit none
      include 'dshmcom.h'
      include 'dsdirver.h'
      character*(*) c
      real*8 dshmaxirho
      real*8 pi
      parameter (pi=3.141592653589793238d0)
      logical firstdd
      data firstdd/.true./
      save firstdd

c...Set default flags
      udfload=.true.  ! load udffile (if chosen by veldf='user') on next
                      ! call to dshmudf.f
      udfearthload=.true. ! load udfearthfile (if chosen by veldfearth='user')
                          ! on next call to dshmudfearth.f

      isodfload=.true.   ! load isodf file on next call to dshmisotrnum.f

c...Three-dimensional velocity dispersion of the WIMPs in the halo
c...Note that in a simple isothermal sphere, the circular speed
c...(approx. solar speed v_sun) is sqrt(2/3)*vd_3d
      vd_3d = 270.0d0     ! WIMP 3D velocity dispersion, km/s
      v_sun = 220.0d0     ! circular velocity at the Sun location, km/s
      v_obs = v_sun       ! Velocity of observer
      vgalesc = 600.d0          ! galactic escape speed in km/s
c...Observer speed w.r.t. the halo, i.e. including Sun + Earth speed,
c...yearly average,
      vobs = 264.d0       ! observer speed w.r.t. the halo in km/s
      v_earth=29.78d0    ! Earth speed in the solar system in km/s

c a few options for the density profile:

c ... modified isothermal profile, smooth profile
      if (c.eq.'isosm') then
        hclumpy=1
        r_0=8.5d0           ! sun galactocentric distance (kpc)
        rho0=0.3d0          ! local halo density (gev cm^-3)
        haloshape='spherical' !here you choose between spherical and axisymm.
        halotype='albega'  !this picks the alpha-beta-gamma profile, set by:
        alphah=2.d0        
        betah=2.d0
        gammah=0.d0
        ah=3.5d0             ! length scale (kpc)
        Rref=r_0
        rhoref=rho0         ! normalizing the profile to the local halo density
        rhcut=1.d-10        ! cut radius (kpc) there is no cut radius for this
                            ! case but this value avoids numerical problems
        haloid='isosm'      ! tag used to identify pbar and dbar files
        veldf='gauss'       ! tag used to identify velocity profile
        veldfearth='sdbest' ! tag to identify earth vel. profile

c...modified isothermal profile, smooth, with dark disk
c...NOTE: Dark disk part only implemented in neutrino telescope routines so far
      elseif (c.eq.'isosmdd') then
        if (firstdd) then
           write(*,*) 
     &       'DS INFO: You have requested halo profile isosmdd,'
           write(*,*) 
     &       'i.e. a modified isothermal sphere with a dark disk.'
           write(*,*)
     &      'Note that currently only the velocity distribution is used'
           write(*,*)
     &      'and this is implemented for neutrino telescope rates'
           write(*,*)
     &      'from the Sun and direct detection only.'
           firstdd=.false.
        endif

        hclumpy=1
        r_0=8.5d0           ! sun galactocentric distance (kpc)
        rho0=0.3d0          ! local halo density (gev cm^-3)
        haloshape='spherical' !here you choose between spherical and axisymm.
        halotype='albega'  !this picks the alpha-beta-gamma profile, set by:
        alphah=2.d0        
        betah=2.d0
        gammah=0.d0
        ah=3.5d0             ! length scale (kpc)
        Rref=r_0
        rhoref=rho0         ! normalizing the profile to the local halo density
        rhcut=1.d-10        ! cut radius (kpc) there is no cut radius for this
                            ! case but this value avoids numerical problems
        haloid='isosm'    ! tag used to identify pbar and dbar files
        veldf='gaussdd'     ! tag used to identify velocity profile
        veldfearth='sdbest' ! tag to identify earth vel. profile

        v_lag_dd=10.d-3 ! velocity of disk rel us (km/s, default is ~0)
        vd_3d_dd=sqrt(3.d0)*50.d0 ! 3D vel. disp. of dark disk
        delta_dd=1.d0 ! density of dark disk with respect to halo local density

c ... modified isothermal profile, clumpy profile
      elseif (c.eq.'isocl') then
        hclumpy=2
        r_0=8.5d0           !sun galactocentric distance (kpc)
        rho0=0.3d0          !local halo density (gev cm^-3)  
        probshape='spherhalo' !here you choose which probability distribution
                              !you have for the small clumps, spherhalo is
                              !the option to have the same distribution as
                              !the spherical density profile, then you have
                              !to set this profile
        halotype='albega'   !this picks the alpha-beta-gamma profile, set by:
        alphah=2.d0        
        betah=2.d0
        gammah=0.d0
        ah=5.d0             ! length scale (kpc)
        Rref=r_0
        rhoref=rho0         ! normalizing the profile to the local halo density
        rhcut=1.d-10        ! cut radius (kpc) there is no cut radius for this
        haloid='isocl'      ! tag used to identify pbar and dbar files
        veldf='gauss'       ! tag used to identify velocity profile
        veldfearth='sdbest' ! tag to identify earth vel. profile

c ... navarro-frenk-white profile, smooth profile
      elseif (c.eq.'nfwsm'.or.c.eq.'default') then
        hclumpy=1
        r_0=8.d0           ! sun galactocentric distance (kpc)
        rho0=0.3d0          ! local halo density (gev cm^-3)
        haloshape='spherical' !here you choose between spherical and axisymm.
        halotype='albega' !this picks the alpha-beta-gamma profile, set by:
        alphah=1.d0        
        betah=3.d0
        gammah=1.d0
        ah=20.d0             ! length scale (kpc)
        Rref=r_0
        rhoref=rho0         ! normalizing the profile to the local halo density
        rhcut=1.d-5         ! cut radius (kpc) 
        haloid='nfwsm'      ! tag used to identify pbar and dbar files
        veldf='gauss'       ! tag used to identify velocity profile
        veldfearth='sdbest' ! tag to identify earth vel. profile


c ... moore et al. profile, smooth profile
      elseif (c.eq.'moosm') then
        hclumpy=1
        r_0=8.d0           ! sun galactocentric distance (kpc)
        rho0=0.3d0          ! local halo density (gev cm^-3)
        haloshape='spherical' !here you choose between spherical and axisymm.
        halotype='albega'  !this picks the alpha-beta-gamma profile, set by:
        alphah=1.5d0        
        betah=3.d0
        gammah=1.5d0
        ah=28.d0             ! length scale (kpc)
        Rref=r_0
        rhoref=rho0         ! normalizing the profile to the local halo density
        rhcut=1.d-5         ! cut radius (kpc) 
        haloid='moosm'      ! tag used to identify pbar and dbar files
        veldf='gauss'       ! tag used to identify velocity profile
        veldfearth='sdbest' ! tag to identify earth vel. profile

c ... burkert et al. profile, smooth profile
      elseif (c.eq.'burksm') then
        hclumpy=1
        r_0=8.d0           ! sun galactocentric distance (kpc)
        haloshape='spherical' ! here you choose between spherical and axisymm.
        halotype='burkert'  ! this picks the burkert profile, set by:
c model with virial mass and virial radius:
c        Mvir=0.130000d+13
c        cvir=0.160000d+02
c this corresponds to:
        ah=11.6839499d0
        rho0=0.339382133d0
        Rref=r_0
        rhoref=rho0          ! normalizing the profile at the local radius
                             ! it can be anything else
        haloid='burksm'      ! tag used to identify pbar and dbar files
        rhcut=0.0d0
        veldf='numc'       ! tag used to identify velocity profile
        udfnumfile=dsinstall
        call dscharadd(udfnumfile,'share/DarkSUSY/vdf-burksm.dat') ! sun vel. prof.
        isodfnumfile=dsinstall
        call dscharadd(isodfnumfile,'share/DarkSUSY/df-burksm.dat') ! sun vel df (3D)
        isodf='num'         ! tag to use the above file
        veldfearth='sdbest' ! tag to identify earth vel. profile
        v_sun = 216.2856d0     ! circular velocity at the Sun location, km/s


c ... adiabatically contracted profile, smooth profile loaded from file
      elseif (c.eq.'adiabsm') then
        hclumpy=1
c model with virial mass and virial radius:
c        Mvir=0.180000d+13
c        cvir=0.120000d+02
        r_0=8.d0           ! sun galactocentric distance (kpc)
        haloshape='spherical' ! here you choose between spherical and axisymm.
        halotype='numerical'  ! this picks the profile from file:
        hmrhofile=dsinstall
        call dscharadd(hmrhofile,'share/DarkSUSY/rho-adiabsm.dat')
                        ! this must be a file with a set of pairs
                        ! dlog(radialcoord[kpc]) - dlog(rho[GeV cm^-3])
                        ! with the format (2(1x,e14.8))
        rho0=dshmaxirho(r_0,0.d0) ! set the local halo density
        haloid='adiabsm'      ! tag used to identify pbar and dbar files
        rhcut=0.0d0
        veldf='numc'   ! tag used to identify velocity profile
        udfnumfile=dsinstall
        call dscharadd(udfnumfile,'share/DarkSUSY/vdf-adiabsm.dat') ! sun vel. prof.
        isodfnumfile=dsinstall
        call dscharadd(isodfnumfile,'share/DarkSUSY/df-adiabsm.dat') ! sun vel df (3D)
        isodf='num'         ! tag to use the above file
        veldfearth='sdbest'  ! tag to identify earth vel. profile
        v_sun = 220.9832d0     ! circular velocity at the Sun location, km/s

      elseif (c.eq.'boer05') then ! de Boer et al, astro-ph/0508617
        hclumpy=1
        r_0=8.3d0           ! sun galactocentric distance (kpc)
        rho0=0.5d0          ! local halo density (gev cm^-3)
        haloshape='boer' !here you choose between spherical and axisymm.
        halotype='boer'  !
        alphah=2.d0        
        betah=2.d0
        gammah=0.d0
        ah=5.0d0             ! length scale (kpc)
        Rref=r_0
        rhoref=rho0         ! normalizing the profile to the local halo density
        rhcut=1.d-10        ! cut radius (kpc) there is no cut radius for this
                            ! case but this value avoids numerical problems
        haloid='deboer'      ! tag used to identify pbar and dbar files
        veldf='gauss'       ! tag used to identify velocity profile
        veldfearth='sdbest' ! tag to identify earth vel. profile
c...The remaining things are for the triaxiality and rings in de Boer's profile

        eps_xy=0.8d0   ! eccentricity in y-direction
        eps_z=0.75d0   ! eccentricity in z-direction
        phi_gc=90.0d0/180.0d0*pi  ! angle between earth-gc and x-axis (rads)
        rho1=4.5d0     ! density in ring 1, GeV / cm^3
        R1=4.15d0      ! radius of ring 1, kpc
        sigr1=4.15d0   ! width of ring 1 r-direction, kpc
        sigz1=0.17d0   ! width of ring 1 z-direction, kpc
        eps_xy1=0.8d0  ! eccentricity of ring 1
        phi1=-70.0d0/180.0d0*pi   ! angle between earth-gc and x-axis (rads)
        rho2=1.85d0    ! density in ring 2, GeV / cm^3
        R2=12.9d0      ! radius of ring 2, kpc
        sigr2=3.3d0    ! width of ring 2 r-direction, kpc
        sigz2=1.7d0    ! width of ring 2 z-direction, kpc
        eps_xy2=0.95d0 ! eccentricity of ring 2
        phi2=-20.0d0/180.0d0*pi   ! angle between earth-gc and x-axis (rads)
        dn=4.0d0       ! cut-off distance for inner side of ring 2
c        rho0=1.92d0             ! JE TMP 

c ... adiabatically contracted profile, smooth profile loaded from file
      elseif (c.eq.'nfwlssm') then ! Strigari NFW profile
        hclumpy=1
c model with virial mass and virial radius:
        r_0=8.5d0           ! sun galactocentric distance (kpc)
        rho0=0.32d0          ! local halo density (gev cm^-3)
        haloshape='spherical' !here you choose between spherical and axisymm.
        halotype='albega' !this picks the alpha-beta-gamma profile, set by:
        alphah=1.d0        
        betah=3.d0
        gammah=1.d0
        ah=21.d0             ! length scale (kpc)
        Rref=r_0
        rhoref=rho0         ! normalizing the profile to the local halo density
        rhcut=1.d-5         ! cut radius (kpc) 
        haloid='nfwlssm'      ! tag used to identify pbar and dbar files

        veldf='numc'   ! tag used to identify velocity profile
        udfnumfile=dsinstall
        call dscharadd(udfnumfile,'share/DarkSUSY/vdf-nfw-strigari.dat') ! sun vel. prof.
        isodfnumfile=dsinstall
        call dscharadd(isodfnumfile,
     &  'share/DarkSUSY/df-nfw-strigari.dat') ! sun vel df (3D)
        isodf='num'         ! tag to use the above file
        veldfearth='sdbest'  ! tag to identify earth vel. profile
        v_sun = 220.d0     ! circular velocity at the Sun location, km/s

c......general einasto profile
      elseif (c.eq.'einasto') then
        hclumpy=1
        r_0=8.d0            ! sun galactocentric distance (kpc)
        rho0=0.3d0          ! local halo density (gev cm^-3)
        haloshape='spherical' !here you choose between spherical and axisymm.
        halotype='einasto'  !this picks the einasto profile, set by:
        rs_e=an03           ! scale radius (kpc)
        rhos_e=rho0/dexp(-2.d0/alphan03*((r_0/an03)**alphan03-1.d0))    ! scale density (GeV cm^-3)          
        n_e=1.d0/alphan03   ! Einasto index
        rhcut=0.d0          ! cut radius (kpc) 
        haloid='einasto'      ! tag used to identify pbar and dbar files
        veldf='gauss'       ! tag used to identify velocity profile
        veldfearth='sdbest' ! tag to identify earth vel. profile

c......einasto profile for segue 1 dwarf galaxy, as per arXiv:0902.4715
      elseif (c.eq.'segue_einasto') then
        hclumpy=1
        r_0=23.d0           ! sun galactocentric distance (kpc)
        rho0=0.3d0          ! local halo density (gev cm^-3)
        haloshape='spherical' !here you choose between spherical and axisymm.
        halotype='einasto'  !this picks the einasto profile, set by:
        rs_e=0.07d0         ! scale radius (kpc)
        rhos_e=3.8d0        ! scale density (GeV cm^-3)          
        n_e=3.3d0           ! Einasto index
        rhcut=0.d0          ! cut radius (kpc) 
        haloid='segue-einasto'      ! tag used to identify pbar and dbar files
        veldf='gauss'       ! tag used to identify velocity profile
        veldfearth='sdbest' ! tag to identify earth vel. profile

c......einasto profile for draco dwarf galaxy, as per arXiv:0902.4715
      elseif (c.eq.'draco_einasto') then
        hclumpy=1
        r_0=80.d0           ! sun galactocentric distance (kpc)
        rho0=0.3d0          ! local halo density (gev cm^-3)
        haloshape='spherical' !here you choose between spherical and axisymm.
        halotype='einasto'  !this picks the einasto profile, set by:
        rs_e=1.0d0          ! scale radius (kpc)
        rhos_e=1.14d0       ! scale density (GeV cm^-3)          
        n_e=3.3d0           ! Einasto index
        rhcut=0.d0          ! cut radius (kpc) 
        haloid='draco-einasto'      ! tag used to identify pbar and dbar files
        veldf='gauss'       ! tag used to identify velocity profile
        veldfearth='sdbest' ! tag to identify earth vel. profile

      else
         write (*,*) 'dshmset: unrecognized option ',c
         stop
      endif


c...unrescaled local density
      rhox = rho0

      return
      end
