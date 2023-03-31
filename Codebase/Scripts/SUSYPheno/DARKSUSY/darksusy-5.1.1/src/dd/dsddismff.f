      subroutine dsddismff(q,a,z,l2jjpp,l2jjnn,l2jjpn)
c_______________________________________________________________________
c  Spin-dependent structure functions for direct detection.
c  Full interacting shell model (ISM) calculations
c  input:
c    q : real*8  : momentum transfer in GeV ( q=sqrt(2ME) )
c    a : integer : mass number
c    z : integer : atomic number
c  output:
c    l2jjnn, l2jjpp, l2jjpn : the factors \lambda^2 J (J+1) made functions
c       of the momentum transfer q, for pp, nn, and pn; they are
c       related to Spp, Snn, and Spn through
c           S = \lambda^2 J (J+1) (2J+1)/\pi
c       In turn, Spp, Snn, and Spn are related to the spin-dependent 
c       structure functions S_{00}(q), S_{01}(q), S_{11}(q) defined in
c       Engel (PLB264,114,1991) via
c           Spp = S00+S11+S01
c           Snn = S00+S11-S01
c           Spn = 2*(S00-S11)
c  author: paolo gondolo (paolo@physics.utah.edu) 2008
c=======================================================================
      implicit none
      include 'dsddcom.h'
      include 'dsmpconst.h'
      include 'dsio.h'
      real*8 q
      integer a,z
      real*8 l2jjpp,l2jjnn,l2jjpn
      real*8 s00,s01,s11,spp,snn,spn
      real*8 y,j
c     j = total nuclear spin
      
      if (a.eq.29.and.z.eq.14) then ! Si-29
         ! ressell et al, PRD 48 (1993) 5519
         y = (1.75d0*q*fermiGeV/2.d0)**2
         j = 0.5d0
         s00 = 0.00818d0*exp(-4.428d0*y)
         s11 = 0.00818d0*1.06d0*exp(-6.264d0*y)
         s01 = 0.00818d0*(-2.06d0)*exp(-5.413d0*y)
      else if (a.eq.73.and.z.eq.32) then ! Ge-73
         j = 4.5d0
         y = (2.04d0*q*fermiGeV/2.d0)**2
         if (ddffsd.eq.'ISMR') then
            ! ressell et al, PRD 48 (1993) 5519
            s00 = 0.20313d0*1.102d0*exp(-7.468d0*y)
            s11 = 0.20313d0*exp(-8.856d0*y)
            s01 = 0.20313d0*(-2.099d0)*exp(-8.191d0*y)
         else
            ! dimitrov, engel, pittel, hep-ph/9408246
            s00 = (0.1606d0+y*(-1.1052d0+y*(3.2320d0+y*
     &           (-4.9245d0+y*(4.1229d0+y*(-1.8016d0+y*0.3211d0))))))
            s11 = (0.1164d0+y*(-0.9228d0+y*(2.9753d0+y*
     &           (-4.8709d0+y*(4.3099d0+y*(-1.9661d0+y*0.3624d0))))))
            s01 = (-0.2736d0+y*(2.0374d0+y*(-6.2803d0+y*
     &           (9.9426d0+y*(-8.5710+y*(3.8310d0-y*0.6948d0))))))
         endif
      else if (a.eq.27.and.z.eq.13) then ! Al-27
         ! engel, ressell, towner, ormand, hep-ph/9504322
         y = (1.73d0*q*fermiGeV/2.d0)**2
         j = 2.5d0
         s00 = (0.0929516d0+y*(-0.472059d0+y*(1.05996d0+y*
     &        (-1.01148d0))))
         s11 = (0.0657232d0+y*(-0.449840d0+y*(1.35041d0+y*
     &        (-1.68508d0))))
         s01 = (0.1563300d0+y*(-0.935958d0+y*(2.45779d0+y*
     &        (-2.72621))))
      else if (a.eq.39.and.z.eq.19) then ! K-39
         ! engel, ressell, towner, ormand, hep-ph/9504322
         y = (1.73d0*q*fermiGeV/2.d0)**2
         j = 1.5d0
         s00 = (0.0094999d0+y*(-0.0619718d0+y*(0.162844d0+y*
     &        (-0.194282d0+y*0.0891054d0))))
         s11 = (0.0298127d0+y*(-0.2176360d0+y*(0.623646d0+y*
     &        (-0.814418d0+y*0.4050270d0))))
         s01 = (0.332044d0+y*(-0.2319430d0+y*(0.638528d0+y*
     &        (-0.798523d0+y*0.3809750d0))))
      else if (a.eq.23.and.z.eq.11) then ! Na-23
         ! ressell, dean, hep-ph/9702290
         y = (1.6864d0*q*fermiGeV/2.d0)**2
         j = 1.5d0
         s00 = (0.0379935d0+y*(-0.174341d0+y*(0.378299d0+y*
     &        (-0.342962d0))))
         s01 = (0.0646525d0+y*(-0.350289d0+y*(0.910031d0+y*
     &        (-0.985833d0))))
         s11 = (0.0275013d0+y*(-0.169641d0+y*(0.507868d0+y*
     &        (-0.617985d0))))
      else if (a.eq.127.and.z.eq.53) then ! I-127
         ! ressell, dean, hep-ph/9702290 Bonn A
         y = (2.282d0*q*fermiGeV/2.d0)**2
         j = 2.5d0
         s00 = exp(-2.d0*y)*(0.0983393d0+y*(-0.489096d0+y*
     &        (1.1402d0+y*(-1.47168d0+y*(1.1717d0+y*
     &        (-0.564574d0+y*(0.158287d0+y*(-0.0238874d0+y*
     &        (0.00154252d0)))))))))
         s11 = exp(-2.d0*y)*(0.0365709d0+y*(-0.194994d0+y*
     &        (0.504876d0+y*(-0.747451d0+y*(0.704334d0+y*
     &        (-0.393018d0+y*(0.121881d0+y*(-0.0191881d0+y*
     &        (0.00121021d0)))))))))
         s01 = exp(-2.d0*y)*(0.11994d0+y*(-0.618424d0+y*
     &        (1.50893d0+y*(-2.07367d0+y*(1.77307d0+y*
     &        (-0.903597d0+y*(0.26002d0+y*(-0.0387025d0+y*
     &        (0.00235675d0)))))))))
      else if (a.eq.129.and.z.eq.54) then ! Xe-129
         ! ressell, dean, hep-ph/9702290 Bonn A
         y = (2.282d0*q*fermiGeV/2.d0)**2
         j = 0.5d0
         s00 = exp(-2.d0*y)*(0.0713238d0+y*(-0.344779d0+y*
     &        (0.755895d0+y*(-0.933448d0+y*(0.690061d0+y*
     &        (-0.302476d0+y*(0.0765282d0+y*(-0.0103169d0+y*
     &        (0.000573919d0)))))))))
         s11 = exp(-2.d0*y)*(-2.05825d0+y*(1.80756d0+y*
     &        (-1.27746d0+y*(0.654589d0+y*(-0.221971d0+y*
     &        (0.0454635d0+y*(-0.00425694d0+y*(-0.000136779d0+y*
     &        (0.00004396d0))))))))+2.11016d0/(1.d0+y))
         s01 = exp(-2.d0*y)*(-0.12166d0+y*(0.644351d0+y*
     &        (-1.52732d0+y*(2.02061d0+y*(-1.57689d0+y*
     &        (0.723976d0+y*(-0.190399d0+y*(0.0263823d0+y*
     &        (-0.00148593d0)))))))))
      else if (a.eq.131.and.z.eq.54) then ! Xe-131
         ! ressell, dean, hep-ph/9702290 Bonn A
         y = (2.282d0*q*fermiGeV/2.d0)**2
         j = 1.5d0
         s00 = exp(-2.d0*y)*(0.0296421d0+y*(-0.133427d0+y*
     &        (0.377987d0+y*(-0.579614d0+y*(0.578896d0+y*
     &        (-0.345562d0+y*(0.115952d0+y*(-0.0201178d0+y*
     &        (0.00141793d0)))))))))
         s11 = exp(-2.d0*y)*(0.0250994d0+y*(-0.137716d0+y*
     &        (0.366609d0+y*(-0.53851d0+y*(0.492545d0+y*
     &        (-0.269903d0+y*(0.0836943d0+y*(-0.0133959d0+y*
     &        (0.000868668d0)))))))))
         s01 = exp(-2.d0*y)*(-0.0545474d0+y*(0.271757d0+y*
     &        (-0.723023d0+y*(1.0545d0+y*(-0.971333d0+y*
     &        (0.538422d0+y*(-0.168988d0+y*(0.027416d0+y*
     &        (-0.00180527d0)))))))))
      else
         if (prtlevel.gt.0) write (*,*)
     &        'dsddismff: ISM form factor unavailable for A=',
     &        a,' Z=',z
         l2jjpp = -1.d0
         l2jjnn = -1.d0
         l2jjpn = -1.d0
         return
      endif

      spp = s00+s11+s01
      snn = s00+s11-s01
      spn = 2.d0*(s00-s11)
      l2jjpp = spp*pi/(2.d0*j+1.d0)
      l2jjnn = snn*pi/(2.d0*j+1.d0)
      l2jjpn = spn*pi/(2.d0*j+1.d0)

      return
      end
