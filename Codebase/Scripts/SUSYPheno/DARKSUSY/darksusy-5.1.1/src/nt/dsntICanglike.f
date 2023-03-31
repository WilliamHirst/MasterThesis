***********************************************************************
*** dsntICanglike returns the contribution of a single event to the 
*** unbinned likelihood, based on the reconstructed angle between
*** its arrival direction in the detector and the direction of the Sun.
***
*** input:  cosphi     cos(reconstructed angle from solar position)
***         parabsigma paraboloid (or other) sigma corresponding to 
***                    1 sigma angular error when considering a single
***                    dimension of a 2D Gaussian PDF (i.e. 39,3% C.L.,
***                    not 68%). (degrees)
***         f_S        signal fraction; percentage of predicted counts
***                    within analysis window (cut cone and superbin)
***                    expected to be due to signal rather than back-
***                    ground.
*** output:            ln(Likelihood / degrees^-1)
***       
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: May 6, 2011
***********************************************************************

      double precision function dsntICanglike(cosphi,
     & parabsigma,f_S)

      implicit none
      include 'dsntIC.h'
      include 'dsge.h'

      real*8 cosphi, parabsigma, f_S
      real*8 signalpartiallike,bgpartiallike,dsntICpsf,dsntICbgangpdf

      !Return low likelihood automatically if cosphi = 1.0 exactly
      if (cosphi .eq. 1.d0) then
        dsntICanglike = bigBadLike
        return
      endif

      !Calculate the signal part of the likelihood
      signalpartiallike = f_S*dsntICpsf(acos(cosphi)*180.d0/pi, 
     & 0.d0, phi_max_deg, parabsigma)

      !Calculate the background part of the likelihood
      bgpartiallike = (1.d0-f_S) * dsntICbgangpdf(cosphi)

      !Put them together
      dsntICanglike = log(signalpartiallike + bgpartiallike)

      end function dsntICanglike
