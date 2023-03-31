
!!**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
!! * antiproton.f *                                galprop package * 2001/05/11
!!**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

      real*8 function ANTIPROTON(key,Pap1,Pp1,NZ1,NA1,NZ2,NA2)        ! IMOS20010511
c***********************************************************************
c             ###   I.Moskalenko   ###    version of 22.06.1998     ###
c Antiproton (+antineitron) production spectrum vs. momentum [barn c/GeV] for 
c pp-, pA-, Ap-, and AA-collisions (Pp1, Pap1 fixed) per 1 target nucleus/cm^3. 
c Refs: Moskalenko I.V. et al. 2002, ApJ 565, 280
c (Tan & Ng 1983, J.Phys.G:Nucl.Phys.9,227; ibid.,p.1289;
c Letaw et al.1983,ApJS,51,271; Gaisser & Schaeffer 1992,ApJ,394,174;
c Westfall et al.1979,PRC,19,1309)
c
c Pap1 [GeV/c] - secondary anti-p momentum; Pp1 [GeV/c] -beam momentum/nucleus
c NA1 & NA2 are the atomic numbers of beam and target nuclei, correspondingly
c (NA1=NA2=1 for pp-collisions)
c***********************************************************************
      implicit real*8 (a-h,m,o-z)
      common/gpmass/M0,G,Mp,Mn,Md,Mpi0,Mpi1,Mmu,MK
     2      /en/Egam,Pp,Ppi,Pe,Pap !GeV, GeV/c
      external TAN_NG

      Pi = 3.1415926
      ANTIPROTON = 0.d0
      Pap = Pap1
      Pp = Pp1/NA1
c      Pp = Pp1  ! use if Pp1 is the beam momentum per nucleon
      s2 = 2.d0*Mp*(dsqrt(Pp**2+Mp**2)+Mp)!square of the total CMS energy
      s1 = dsqrt(s2)                      !Total pp-energy in CMS
      if(s1 .le. 4.d0*Mp) return          !anti-p production threshold

      MX = 3.d0*Mp                        !pp->ppp(anti-p) channel
      gam_cms = s1/2.d0/Mp                !CMS Lorentz factor (Lf)
      betgam_cms = dsqrt(gam_cms**2-1.d0) !CMS beta*gamma
      Eap = dsqrt(Pap**2+Mp*Mp)           !anti-p LS energy
c###
      if(NA2 .gt. 1) then                 ! IMOS20010215
         Eap = Eap+0.06                   ! gives more pbars at low energies
         Pap = dsqrt(Eap**2-Mp*Mp)        ! in case of nucleus target
      endif 
c###
      Eap_max = (s2-MX*MX+Mp*Mp)/2.d0/s1  !max anti-p CMS energy
      Ek = dsqrt(Pp*Pp+Mp*Mp)-Mp          !kinetic energy per nucleon
      COSmax =-1.d0
      if(betgam_cms*Pap .gt. 0.d0)
     1   COSmax =(gam_cms*Eap-Eap_max)/betgam_cms/Pap
      if(COSmax .lt. -1.d0) COSmax = -1.d0      
      AI = 0.d0
      if(COSmax .lt. 1.d0)
     #   call SIM1(1.d0,COSmax,1.d-3,1.d-4,1.d-40,TAN_NG,AI)
c dF/dEap - production spectrum vs. energy; factor 2 accounts for antineutrons
      ANTIPROTON = -AI*Pap  *2.d0  *2.d0*Pi*1.d-3  
     #   *Pap/Eap                         !tranformation to dF/dPap
      if(NA1*NA2 .eq. 1) return           !return if pp-collisions

c pA-, Ap-, and AA-collisions: 
      if(Ek .le. 0.3) return              !kinetic energy must be>0.3GeV
      call NUCLEON_CS(key,Ek,1,NZ2,NA2
     #       ,PP_inel,PA_inel,aPP_non,aPA_non,aPP_ann,aPA_ann)  ! IMOS20010511 IMOS20000606
c nuclei: pA total inelastic cross section, mb (Letaw et al.1983,ApJS,51,271)
      CS2 = PP_inel
      if(NA2 .gt. 1) CS2 = PA_inel
      CS1 = PP_inel
      if(NA1 .gt. 1) then 
         call NUCLEON_CS(key,Ek,1,NZ1,NA1
     #       ,PP_inel,PA_inel,aPP_non,aPA_non,aPP_ann,aPA_ann)  ! IMOS20010511 IMOS20000606
         CS1 = PA_inel
      endif
c multiplicity of anti-p in pA-,Ap-,AA-coll.; Gaisser&Schaeffer 1992,ApJ,394,174
c###
      MULT=1.2*  (NA1*CS2+NA2*CS1)/2./PP_inel   ! IMOS20010223 put factor 1.2
c###                                            ! from comparison w. Simon et al. 
      ANTIPROTON = ANTIPROTON*MULT

      return
      end

      real*8 function TAN_NG(COSX)
c***********************************************************************
c The invar. cross section of the inclusive anti-proton production 
c [mbarn c^3/GeV^2] from p+p->(anti-p)+X reaction (for Pp,(Pap,COSX) fixed);
c (Tan & Ng 1983,J.Phys.G:Nucl.Phys.9,227; ibid., p.1289)
c COSX = cos(theta) is the LS polar angle.
c***********************************************************************
      implicit real*8 (a-h,m,o-z)
      common/gpmass/M0,G,Mp,Mn,Md,Mpi0,Mpi1,Mmu,MK
     2      /en/Egam,Pp,Ppi,Pe,Pap !GeV, GeV/c

      TAN_NG = 0.d0
      s2 = 2.d0*Mp*(dsqrt(Pp**2+Mp**2)+Mp)!square of the total CMS energy
      s1 = dsqrt(s2)                      !Total pp-energy in CMS
      MX = 3.d0*Mp                        !pp->ppp(anti-p) channel
      if(s1 .le. 4.d0*Mp) return          !anti-p production threshold

      gam_cms = s1/2.d0/Mp                !CMS Lorentz factor (Lf)
      betgam_cms = dsqrt(gam_cms**2-1.d0) !CMS beta*gamma
      Eap = dsqrt(Pap**2+Mp*Mp)           !anti-p LS energy
      SINX = 0.d0
      if(1.d0-COSX**2 .gt. 0.d0) SINX = dsqrt(1.d0-COSX**2)
      PT = Pap*SINX                       !anti-p transverse momentum
c radial variable: XR=(E*)/(Emax*), where E*,Emax* - anti-p CMS energies
      XR=(gam_cms*Eap-betgam_cms*Pap*COSX)*2.d0*s1/(s2-MX*MX+Mp*Mp)
      if(XR .gt. 1.d0) return

c      goto 999
c >> parametrization used in Tan & Ng 1983,J.Phys.G:Nucl.Phys.9,1289 <<
c inclusive cross section at s^(1/2) >~ 10 GeV; pp.1296-8
888   A = 0.465*dexp(-3.70d-2*XR)+2.31*dexp(1.40d-2*XR)   ! c/GeV
      B = 3.02d-2*dexp(-3.19*(XR+0.399))*(XR+0.399)**8.39 !(c/GeV)**2
      F = (3.15-1.05d-4)*(1.d0-XR)**7.90                  !mbarn c^3/GeV^2
      if(XR .le. 0.5) F = F+1.05d-4*dexp(-10.1*XR)
      TAN_NG = F*dexp(-(A*PT+B*PT*PT))
      if(s1 .ge. 1.d1) return

c inclusive cross section at s^(1/2) <~ 10 GeV; pp.1303-5
      XT =((gam_cms*Eap-betgam_cms*Pap*COSX)-Mp)
     #   /((s2-MX*MX+Mp*Mp)/2.d0/s1-Mp)   !XT=(T*)/(Tmax*)
      Q = s1-4.d0*Mp                      !4Mp = anti-p production threshold
      A = 0.306*dexp(-0.120*XT)
      B = 0.0552*dexp(2.72*XT)
      C = 0.758-0.680*XT+1.54*XT*XT
      D = 0.594*dexp(2.87*XT)
      delta = 1.-dexp(-(1.-dexp(-A*Q**B))*dexp(C*Q-D))
      TAN_NG = TAN_NG/delta
      return

c ----------------------------------------------------------------------
c >> parametrization used in Tan & Ng 1983,J.Phys.G:Nucl.Phys.9,227 <<
c inclusive cross section at s^(1/2) >~ 10 GeV
999   A = 3.95*dexp(-2.76*XR)             ! c/GeV
      B = 40.5*dexp(-3.21*XR)*XR**2.13    !(c/GeV)**2
      F = 2.10*(1.d0-XR)**7.80            !mbarn c^3/GeV^2
      if(XR .le. 0.5) F = F+3.34*dexp(-17.6*XR)
      TAN_NG = F*dexp(-(A*PT+B*PT*PT))
      if(s1 .ge. 1.d1) return

c inclusive cross section at s^(1/2) <~ 10 GeV
      DXR =XR -dsqrt(PT**2+Mp**2) *2.d0*s1/(s2-MX*MX+Mp*Mp) !DXR=XR-XRmin
      Q = s1-4.d0*Mp                      !4Mp = anti-p production threshold
      delta = 6.25d-3*(dexp(-0.592*Q)+493.*dexp(-5.40*Q))
     #   *(dexp(6.08+2.57*DXR+7.95*DXR**2)-1.d0)*dexp(3.00*DXR*(3.09-Q))
      TAN_NG = (delta+1.d0)*TAN_NG
      return
      end

