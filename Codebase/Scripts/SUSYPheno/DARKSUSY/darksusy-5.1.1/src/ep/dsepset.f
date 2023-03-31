      subroutine dsepset(c)
c...set parameters for positron routines
c...  c - character string specifying choice to be made
c...author: joakim edsjo, 2000-07-09
c...modified: paolo gondolo, 2000-07-19
c...          joakim edsjo,  2000-08-15
c...          Joakim Edsjo, 2008-04-17 (reduced vtol: 0.01 -> 0.001)
      implicit none
      include 'dshmcom.h'
      include 'dsepcom.h'
      character*(*) c


c...Make sure to (re)load tables upon next call to dsepdiff
      epload=.true.      ! for src/ep routines
      epcl_load=.true.   ! for src/ep2 routines

c... define some internal parameters
      sumtol=0.01d0
      rtol=1.0d0
      vtol=0.001d0
      rwid=3.0d0

c... fetch some parameters from the halo model common block
      r_e=r_0           ! galactocentric distance (kpc)
      n_c = 1.0d0       ! we now have the normalization in the halo
                        ! density dshmrho directly, hence n_c=1.0

c...baltz and edsjo 1999, without low-energy cut-off in k_diff
      if (c.eq.'be99nocut') then
         dsepdiffmethod=0
         l_h=3.0d0              ! height of diffusion zone (kpc)
         tau16=1.0d0            ! energy loss = (e/gev)^2/tau (tau in 10^16 s)
         alphexp=0.6d0          ! diffusion const = k0 (e/gev)^alphexp
         k27=3.15d0             ! k0 (10^27 cm^2/s)
         epc='be99nocut'        ! idtag for model identification

c...baltz and edsjo 1999, with low-energy cut-off in k_diff
      else if (c.eq.'be99') then
         dsepdiffmethod=1
         l_h=3.0d0              ! height of diffusion zone (kpc)
         tau16=1.0d0            ! energy loss = (e/gev)^2/tau (tau in 10^16 s)
         alphexp=0.6d0          ! diffusion const = k0 (e/gev)^alphexp
         k27=6.08952d0          ! k0 (10^27 cm^2/s)
         epc='be99'             ! idtag for model identification

c...baltz and edsjo 1999, with low-energy cut-off in k_diff (default)
      else if (c.eq.'esu04'.or.c.eq.'default') then
         dsepdiffmethod=1
         l_h=4.0d0              ! height of diffusion zone (kpc)
         tau16=1.0d0            ! energy loss = (e/gev)^2/tau (tau in 10^16 s)

c in MS D_xx= beta * D_0 *(Rig/Rig0)**delta with Rig0 = 4 GV, 
c delta=0.6 iff Rig>Rig0, 0 iff Rig<Rig0, D_0 = 2.5d28 cm^2/s
c and beta = momentum/energy 
c the following should match the high energy behaviour 
c (but not the low energy one)
         alphexp=0.6d0          ! diffusion const = k0 (e/gev)^alphexp
c
c NOTE: this was wrong because h(eps) in the code is not as written in
c eqs (7) in the paper but it is instead 1+(eps/3)**0.6
c
         k27=25.d0*(3.d0/4.d0)**alphexp             ! k0 (10^27 cm^2/s)
         epc='esu04'            ! idtag for model identification


c...kamionkowski and turner, energy dependent escape time
      else if (c.eq.'kt91') then
         dsepdiffmethod=2
         epc='kt91'             ! idtag for model identification


c...moskalenko and strong, with reacceleration
      else if (c.eq.'ms99') then
         dsepdiffmethod=3
         data mselo/1.03d0,2.06d0,5.15d0,10.30d0,25.76d0,51.52d0,
     +        103.0d0,206.1d0,412.1d0,824.3d0/
         data msatab/-1.9732d0,-2.4853d0,-2.6365d0,-1.9555d0,-1.1684d0,
     +        -0.8469d0,-0.6979d0,-0.6173d0,-0.5337d0,-0.4255d0/
         data msbtab/2.3448d0,3.2517d0,4.4743d0,4.4101d0,3.7535d0,
     +        3.3985d0,3.3043d0,3.3187d0,3.2338d0,2.9458d0/
         data msctab/-0.1340d0,-0.6564d0,-1.6189d0,-2.2099d0,-2.6853d0,
     +        -3.0180d0,-3.4776d0,-4.0118d0,-4.4403d0,-4.6095d0/
         
         data msehi/1.03d0,2.06d0,5.15d0/
         data mswtab/-4.9292d0,-7.2475d0,-8.9618d0/
         data msxtab/-0.8786d0,1.0942d0,2.1785d0/
         data msytab/-0.1115d0,0.4742d0,3.1988d0/
         epc='ms99'            ! idtag for model identification
         
c...galprop green's functions
      else if (c.eq.'galprop') then
         dsepdiffmethod=4
         epc='galprop'

c...invalid choice
      else
         write (*,*) 'dsepset: unrecognized option ',c
         stop
      endif

      return
      end
