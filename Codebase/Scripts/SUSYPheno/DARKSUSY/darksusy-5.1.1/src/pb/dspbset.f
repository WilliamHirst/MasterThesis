      subroutine dspbset(c)
c...set parameters for antiproton routines
c...  c - character string specifying choice to be made
c...author: paolo gondolo 1999-07-14
      implicit none
      include 'dshmcom.h'
      include 'dspbcom.h'
      include 'dspbprivate.h'
      real*8 dshmrho2cylint
      character*(*) c
      integer first
      data first/0/
      save first


      if (first.eq.0) then
c     call routine to initialize zeros of bessel j0
c        call dspbzerobessj0
      endif

c...chardonnet et al 1996
      if (c.eq.'cha96') then
         pbc='cha96'   ! Used if a pbtd-file is used
         pbpropmodel=1
         pbr0=8.0d0
         pbnzero=30
         pbhg=0.1d0             !(half thickness of the disk in kpc)
         pbhh=3.d0              !(height of the diffusion box in kpc)
         pbrh=20.d0             !(radius of the diffusion box in kpc)
         pbrig0=3.d0            !(proton rigidity in diffusion constant gev)
         pbkg=0.d0              !(diffusion coeff in gas in 10^27 cm^2/s)
         pbkh=6.d0              !(diffusion coeff in halo in 10^27 cm^2/s
         pbng=1.d0              !(density of hydrogen in the disk in cm^-3)
         pbdelta=0.6d0          !(exponent of rigidity dep in diff. koeff.)
         pbnh=0.d0              !(density of hydrogen in the halo in cm^-3)
         pbcvel=0.d0            !(galactic wind velocity in km/s)
         smod=0.320d0           !(solar minimum)
c...   numerical parameters
         pbrcy=0.0d0            ! inner region for split integration
         pbzcy=0.0d0            ! inner region for split integration

c...bottino et al 1998
      else if (c.eq.'bot98') then
         pbc='bot98'   ! Used if a pbtd-file is used
         pbpropmodel=1
         pbr0=8.0d0
         pbnzero=30
         pbhg=0.1d0             !(half thickness of the disk in kpc)
         pbhh=3.d0              !(height of the diffusion box in kpc)
         pbrh=20.d0             !(radius of the diffusion box in kpc)
         pbrig0=1.d0            !(proton rigidity in diffusion constantgev)
         pbkg=0.d0              !(diffusion coeff in gas in 10^27 cm^2/s)
         pbkh=6.d0              !(diffusion coeff in halo in 10^27 cm^2/s
         pbdelta=0.6d0          !(exponent of rigidity dep in diff. koeff.)
         pbng=1.d0              !(density of hydrogen in the disk in cm^-3)
         pbnh=0.d0              !(density of hydrogen in the halo in cm^-3)
         pbcvel=0.d0            !(galactic wind velocity in km/s)
         smod=0.320d0           !(solar minimum)
c...   numerical parameters
         pbrcy=0.0d0            ! inner region for split integration
         pbzcy=0.0d0            ! inner region for split integration


c...bergstrom et al 1999
      else if (c.eq.'beu99') then
         pbc='beu99'   ! Used if a pbtd-file is used
         pbpropmodel=2
         pbr0=8.5d0
         pbnzero=30
         pbrh=20.d0             !(radius of the diffusion box in kpc)
         pbhg=0.1d0             !(half thickness of the disk in kpc)
         pbhh=3.d0              !(height of the diffusion box in kpc)
         pbng=1.d0              !(density of hydrogen in the disk in cm^-3)
         pbnh=0.d0              !(density of hydrogen in the halo in cm^-3)
         pbrig0=3.d0            !(proton rigidity in diffusion constant gev)
         pbkg=6.d0              !(diffusion coeff in gas in 10^27 cm^2/s)
         pbkh=6.d0              !(diffusion coeff in halo in 10^27 cm^2/s)
         pbdelta=0.6d0          !(exponent of rigidity dep in diff. koeff.)
         pbcvel=0.d0            !(galactic wind velocity in km/s)
         smod=0.320d0           !(solar minimum)
c...   numerical parameters
         pbrcy=0.0d0            ! inner region for split integration
         pbzcy=0.0d0            ! inner region for split integration

c...bergstrom et al 1999
c...but with DC-like setup in moskalenko et al. ApJ 565 (2002) 280
      else if (c.eq.'beu04'.or.c.eq.'default') then
         pbc='beu04'   ! Used if a pbtd-file is used
         pbpropmodel=3
         pbr0=r_0
         pbnzero=30
         pbrh=30.d0             !(radius of the diffusion box in kpc)
         pbhg=0.1d0             !(half thickness of the disk in kpc)
         pbhh=4.d0              !(height of the diffusion box in kpc)
         pbng=1.d0              !(density of hydrogen in the disk in cm^-3)
         pbnh=0.d0              !(density of hydrogen in the halo in cm^-3)
         pbrig0=4.d0            !(proton rigidity in diffusion constant gev)
         pbkg=25.d0              !(diffusion coeff in gas in 10^27 cm^2/s)
         pbkh=25.d0              !(diffusion coeff in halo in 10^27 cm^2/s)
         pbdelta=0.6d0          !(exponent of rigidity dep in diff. koeff.)
         pbcvel=10.d0            !(galactic wind velocity in km/s)
         smod=0.320d0           !(solar minimum)
c...   numerical parameters
         write(*,*) 'WARNING in dspbset: using halo model ',haloid
         write(*,*) 'for the setup. If you change halo model, ',
     &     'make sure'
         write(*,*) 'to call dspbset again.'
         if (haloid.eq.'adiabsm') then
           pbrcy=1.0d0            ! inner region for split integration
           pbzcy=1.0d0            ! inner region for split integration
           pbrho2int=dshmrho2cylint(pbrcy,pbzcy)
         else
           pbrcy=0.0d0            ! inner region for split integration
           pbzcy=0.0d0            ! inner region for split integration
         endif

      else if (c.eq.'beu05') then
         pbc='beu05'   ! Used if a pbtd-file is used
         pbpropmodel=3
         pbr0=r_0
         pbnzero=30
         pbrh=30.d0             !(radius of the diffusion box in kpc)
         pbhg=0.1d0             !(half thickness of the disk in kpc)
         pbhh=4.d0              !(height of the diffusion box in kpc)
         pbng=1.d0              !(density of hydrogen in the disk in cm^-3)
         pbnh=0.d0              !(density of hydrogen in the halo in cm^-3)
         pbrig0=4.d0            !(proton rigidity in diffusion constant gev)
         pbkg=25.d0              !(diffusion coeff in gas in 10^27 cm^2/s)
         pbkh=25.d0              !(diffusion coeff in halo in 10^27 cm^2/s)
         pbdelta=0.6d0          !(exponent of rigidity dep in diff. koeff.)
         pbcvel=10.d0            !(galactic wind velocity in km/s)
         smod=0.610d0           !(BESS 98 value)
c...   numerical parameters
         write(*,*) 'WARNING in dspbset: using halo model ',haloid
         write(*,*) 'for the setup. If you change halo model, ',
     &     'make sure'
         write(*,*) 'to call dspbset again.'
         if (haloid.eq.'adiabsm') then
           pbrcy=1.0d0            ! inner region for split integration
           pbzcy=1.0d0            ! inner region for split integration
           pbrho2int=dshmrho2cylint(pbrcy,pbzcy)
         else
           pbrcy=0.0d0            ! inner region for split integration
           pbzcy=0.0d0            ! inner region for split integration
         endif

c... Donato,Fornengo,Maurin,Salati,Tailet (astro-ph/0306207)
      else if (c.eq.'salati03med') then !Michael Gustafsson 2006-01-18
         pbc='salati03med'   ! Used if a pbtd-file is used
         pbpropmodel=3
         pbr0=r_0     
         pbnzero=30   !Besselfunction
         pbrh=20.d0   !ok          !(radius of the diffusion box in kpc)
         pbhg=0.1d0   !ok          !(half thickness of the disk in kpc)
         pbhh=4.d0    !=L          !(height of the diffusion box in kpc)
         pbng=1.d0    !ok          !(density of hydrogen in the disk in cm^-3)
         pbnh=0.d0    !ok          !(density of hydrogen in the halo in cm^-3)
         pbrig0=1.d0  !?           !(proton rigidity in diffusion constant gev)
         pbkg=3.379d0 !=K_0        !(diffusion coeff in gas in 10^27 cm^2/s)
         pbkh=3.379d0 !=K_0        !(diffusion coeff in halo in 10^27 cm^2/s)
         pbdelta=0.70d0            !(exponent of rigidity dep in diff. koeff.)
         pbcvel=12.d0 !=V_c in z-dir, V_A notused?    !(galactic wind velocity in km/s)
         smod=0.610d0 !           !(BESS 98 value)
c...   numerical parameters
         write(*,*) 'WARNING in dspbset: using halo model ',haloid
         write(*,*) 'for the setup. If you change halo model, ',
     &     'make sure'
         write(*,*) 'to call dspbset again.'
         if (haloid.eq.'adiabsm') then
           pbrcy=1.0d0            ! inner region for split integration
           pbzcy=1.0d0            ! inner region for split integration
           pbrho2int=dshmrho2cylint(pbrcy,pbzcy)
         else
           pbrcy=0.0d0            ! inner region for split integration
           pbzcy=0.0d0            ! inner region for split integration
         endif

      else if (c.eq.'salati03max') then !Michael Gustafsson 2006-01-18
         pbc='salati03max'   ! Used if a pbtd-file is used
         pbpropmodel=3
         pbr0=r_0     
         pbnzero=30   
         pbrh=20.d0             !(radius of the diffusion box in kpc)
         pbhg=0.1d0             !(half thickness of the disk in kpc)
         pbhh=15.d0   !=L       !(height of the diffusion box in kpc)
         pbng=1.d0              !(density of hydrogen in the disk in cm^-3)
         pbnh=0.d0              !(density of hydrogen in the halo in cm^-3)
         pbrig0=1.d0            !(proton rigidity in diffusion constant gev)
         pbkg=23.08d0 !=K_0     !(diffusion coeff in gas in 10^27 cm^2/s)
         pbkh=23.08d0 !=K_0     !(diffusion coeff in halo in 10^27 cm^2/s)
         pbdelta=0.46d0         !(exponent of rigidity dep in diff. koeff.)
         pbcvel=5.d0  !=V_c     !(galactic wind velocity in km/s)
         smod=0.610d0           !(BESS 98 value)
c...   numerical parameters
         write(*,*) 'WARNING in dspbset: using halo model ',haloid
         write(*,*) 'for the setup. If you change halo model, ',
     &     'make sure'
         write(*,*) 'to call dspbset again.'
         if (haloid.eq.'adiabsm') then
           pbrcy=1.0d0            ! inner region for split integration
           pbzcy=1.0d0            ! inner region for split integration
           pbrho2int=dshmrho2cylint(pbrcy,pbzcy)
         else
           pbrcy=0.0d0            ! inner region for split integration
           pbzcy=0.0d0            ! inner region for split integration
         endif

      else if (c.eq.'salati03min') then !Michael Gustafsson 2006-01-18
         pbc='salati03min'   ! Used if a pbtd-file is used
         pbpropmodel=3
         pbr0=r_0     
         pbnzero=30   
         pbrh=20.d0              !(radius of the diffusion box in kpc)
         pbhg=0.1d0              !(half thickness of the disk in kpc)
         pbhh=1.0d0     !=L       !(height of the diffusion box in kpc)
         pbng=1.d0               !(density of hydrogen in the disk in cm^-3)
         pbnh=0.d0               !(density of hydrogen in the halo in cm^-3)
         pbrig0=1.d0             !(proton rigidity in diffusion constant gev)
         pbkg=0.4827d0 !=K_0      !(diffusion coeff in gas in 10^27 cm^2/s)
         pbkh=0.4827d0 !=K_0      !(diffusion coeff in halo in 10^27 cm^2/s)
         pbdelta=0.85d0          !(exponent of rigidity dep in diff. koeff.)
         pbcvel=13.5d0 !=V_c      !(galactic wind velocity in km/s)
         smod=0.610d0             !(BESS 98 value)
c...   numerical parameters
         write(*,*) 'WARNING in dspbset: using halo model ',haloid
         write(*,*) 'for the setup. If you change halo model, ',
     &     'make sure'
         write(*,*) 'to call dspbset again.'
         if (haloid.eq.'adiabsm') then
           pbrcy=1.0d0            ! inner region for split integration
           pbzcy=1.0d0            ! inner region for split integration
           pbrho2int=dshmrho2cylint(pbrcy,pbzcy)
         else
           pbrcy=0.0d0            ! inner region for split integration
           pbzcy=0.0d0            ! inner region for split integration
         endif

      else if (c.eq.'galprop') then ! use Green's functions from GALPROP
         pbc='galprop'

c...invalid choice
      else
         write (*,*) 'dspbset: unrecognized option ',c
         stop
      endif


c...To find good values of pbrcy and pbzcy, the following test may
c...prove useful.
c      do i=1,4
c        t=tpcheck(i) 
c        res1=dspbtd15point(pbrho2int,t)
c        write(*,*) 'res1 = ',res1
c        res2=dspbtd15comp(t)
c        write(*,*) 'res1, res2, res1/res2 = ',res1,res2,res1/res2 
c      enddo  


      return
      end





