***********************************************************************
*** dsntICbounds returns counts, likelihoods and p-values from IceCube
*** input: liketype	Sets combination of data to use in likelihood
***			 calculations
***			1 => Number of events only
***			2 => Number of events and event arrival angles
***			3 => Number of events and energy estimator (nchan 
***			     = number of lit DOMs)
***			4 => Number of events, event arrival angles and 
***			     energy estimator
***        pvalFromRef  T => calculate the p-value with reference to a user-
***                          specified likelihood, assuming that the log-
***                          likelihood follows a chi^2 distriubution in
***                          the vicinity of the reference point (i.e.
***                          essentially assume that the reference corresponds
***                          to the best-fit log likelihood, and that Wilks' 
***                          theorem is satisfied).
***                     F => calculate the p-value with reference to the
***                          observed background
***        referenceLike reference value of the (natural) log likelihood to
***                      use for the calculation of the p-value if 
***                      pvalFromRef=T; ignored otherwise.
***        dof          degrees of freedom to assume in the calculation of the
***                      p-value in the case that pvalFromRef=T; ignored 
***                      otherwise.  This should normally be set to the
***                      number of parameters over which the likelihood has
***                      not been maximised when the solution is perturbed 
***                      from the minimum (i.e. the number of free parameters
***                      in a scan minus the number profiled over).
***                              
*** output: Nsignal_predicted	Predicted number of signal events within 
***				phi_cut (includes solar coronal BG)
***	    NBG_expected	Expected number of (non-solar) background
***				events insided phi_cut
***         Ntotal_observed	Observed number of events inside phi_cut
***         lnlike		natural log of chosen likelihood
***	    pvalue		derived p-value for chosen model
***        
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: Mar 20, 2011
*** Updated: Jul 21, 2011
***********************************************************************


      subroutine dsntICbounds(Nsignal_predicted, NBG_expected, 
     & Ntotal_observed, lnlike, pvalue, liketype, 
     & pvalFromRef, referenceLike, dof)

      implicit none
      include 'dswacom.h'
      include 'dsntcom.h'
      include 'dsntIC.h'

      integer Ntotal_observed, liketype, i, j, k, l
      integer counted1, counted2, countrate
      real*8 Nsignal_predicted, NBG_predicted, NBG_expected, dsntICpval
      real*8 lnlike, pvalue, referenceLike, dof, DGAMIC, DGAMMA
      real*8 nLikelihood, angularLikelihood, spectralLikelihood, scaling
      real*8 theta_tot, f_S, dsntICanglike, dsntICspeclike, dsntICnlike
      real*8 erval1, erval2, BGpvalChi2, deltalnlike
      real*8 dsntcapsun, dsntcapsunnum, dsntcapsuntab, ca, tt_sun
      logical pvalFromRef, dsntICspeclike_reset, doProfiling
      character (len=*) pref,f1,f2,f3,f4
      !Hidden option for doing speed profiling
      parameter (doProfiling = .false.)
      parameter (pref = 'share/DarkSUSY/IC_data/')
      parameter (f1 = pref//'events_10deg_IC79.dat')
      parameter (f2 = pref//'energy_histograms_IC79.dat')
      parameter (f3 = pref//'BG_distributions_IC79.dat')
      parameter (f4 = pref//'nuEffArea_IC79.dat')
        
      if (doProfiling) call system_clock(counted1,countrate)

      !Check validity of user-provided liketype
      if (liketype .lt. 1 .or. liketype .gt. 4) then
        write(*,*) 'Invalid liketype provided to dsntICbounds.'
        write(*,*) 'Quitting...'
        stop
      endif
     
      !if dsntICinit has not yet been called, call it with the default IC-79 options
      if (.not. dsntICinit_called) then
        call dsntICinit(f1,f2,f3,f4,20.d0,0.05d0,.true.,.true.)
      endif

      !Set internal wimp mass
      log10mwimp = dlog10(wamwimp)

      !Calculate capture and annihlation rates once and for all
      ! csu: capture rate (s^-1)
      if (ntcalcmet.eq.1) then      ! jkg approximation
        csu=dsntcapsun(wamwimp,wasigsip,wasigsdp)
      elseif (ntcalcmet.eq.2) then  ! Full Gould, analytic Gaussian approx.
        csu=dsntcapsun(wamwimp,wasigsip,wasigsdp)
      elseif (ntcalcmet.eq.4) then  ! Full Gould, num. integrating vel.dist.
        if (nttab.eq.0) then        ! use numerical integration everytime
          csu=dsntcapsunnum(wamwimp,wasigsip,wasigsdp)
        else
          csu=dsntcapsuntab(wamwimp,wasigsip,wasigsdp)
        endif
      endif
      ca=wasv/6.6d28*(wamwimp/20.d0)**(3./2.)
      tausu=1.0d0/dsqrt(csu*ca)                      !Equilibration time (s)
      tt_sun=1.5d17*dsqrt(csu*ca)
      annrateIC=csu*0.5d0*tanh(tt_sun)**2            !Annihilation rate (s^-1)

      if (doProfiling) then
        call system_clock(counted2,countrate)
        write(*,*) 'Elapsed time initialising (s): ', 
     &   real(counted2 - counted1)/real(countrate)
      endif  

      !Calculate signal counts and spectrum
      call dsntICsignal

      if (doProfiling) then
        call system_clock(counted1,countrate)
        write(*,*) 'Elapsed time on signal calc (s): ', 
     &   real(counted1 - counted2)/real(countrate)
      endif  

      !Calculate likelihood
      lnlike = 0.d0
      Nsignal_predicted = 0.d0
      NBG_expected = 0.d0
      !Step through the superbins
      do i = 1, nBinsEAError
        angularLikelihood = 0.d0
        spectralLikelihood = 0.d0
        theta_tot = theta_BG(i) + theta_S(i)
        f_S = theta_S(i) / theta_tot
        if (liketype .ne. 1) then
          !Reset the saved spectral likelihoods next time dsntICspeclike is run
          dsntICspeclike_reset = .true.
          !Step through the events in each superbin
          do j = 1, nEvents_inEAErrBins(i)
            !Add in angular likelihood for this event
            if (liketype .eq. 2 .or. liketype .eq. 4) 
     &       angularLikelihood = angularLikelihood + dsntICanglike(
     &        events_cosphi(i,j),events_cosphiErr(i,j),f_S)
            !Add in spectral likelihood for this event
            if (liketype .eq. 3 .or. liketype .eq. 4) 
     &       spectralLikelihood = spectralLikelihood + dsntICspeclike(
     &        events_nchan(i,j),theta_S(i),f_S,dsntICspeclike_reset,
     &        EAlogE_inEAErrBins(1,i),EAlogE_inEAErrBins(2,i))
          enddo
        endif

        if (doProfiling) then
          call system_clock(counted2,countrate)
          write(*,*) 'Elapsed time on ang/spec likelihood calc (s): ', 
     &     real(counted2 - counted1)/real(countrate)
        endif  

        !Calculate the number likelihood for this superbin
        nLikelihood = dsntICnlike(nEvents_inEAErrBins(i),
     &   theta_tot,theta_S(i),EAErr_inEAErrBins(i),theoryErr)

        if (doProfiling) then
          call system_clock(counted1,countrate)
          write(*,*) 'Elapsed time on number likelihood calc (s): ', 
     &     real(counted1 - counted2)/real(countrate)
        endif  

        !Put together the number, angular and spectral likelihoods
        lnlike = lnlike + nLikelihood + angularLikelihood + 
     &   spectralLikelihood
      enddo 

      !Work out total counts
      NBG_expected = sum(theta_BG)
      Nsignal_predicted = theta_S_total
      Ntotal_observed = sum(nEvents_inEAErrBins)
      
      !Calculate pvalue
      if (pValFromRef) then !compute p-value with reference to input log likelihood

        deltalnlike = max(0.d0, referenceLike-lnlike)
        pvalue = DGAMIC(dof*0.5d0,deltalnlike)/DGAMMA(dof*0.5d0)

      else                  !compute p-value with reference to background

        theta_tot = NBG_expected + theta_S_total
        !p-value from Poissonian statistics
        if (.not. pvalBGPoisComputed) call dsntICbglikeprecomp
        pvalue = dsntICpval(Ntotal_observed, theta_tot, theta_S_total)
        pvalue = pvalue / BGpvalPoissonian
        
      endif

      if (doProfiling) then
        call system_clock(counted2,countrate)
        write(*,*) 'Elapsed time on pval calc (s): ',
     &   real(counted2 - counted1)/real(countrate)
      endif  

      end subroutine dsntICbounds

