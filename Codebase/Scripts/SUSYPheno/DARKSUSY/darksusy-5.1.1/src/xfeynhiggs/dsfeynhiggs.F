       subroutine dsfeynhiggs(hwar)

c -----------------------------------------------------------------       
c Implementation of FeynHiggs in DarkSUSY by
c Author: J. Edsjo, edsjo@physto.se
c Date: 2008-06-24.
c Modified: 2008-07-02 by J. Edsjo, small bug fixes.
c Based on implementation for FeynHiggs 2.4.1 by P. Gondolo
c Modified:2009-03-16 by E. Lundstrom, extraction of (partial)
c widhts and production cross section ratios (needed in Higgs
c Bounds) from FeynArts.
c 2010-09-24: Updated to work with FeynHiggs 2.7.2
c 2011-09-04: PG update to FeynHiggs 2.8.5
c 2013-02-09: PG update to FeynHiggs 2.9.4
c
c Output: mh is lighter scalar Higgs mass, HM is heavier Higgs mass 
c -----------------------------------------------------------------

c --------------------------------------------------------------
c
c Warnings implemented by Joakim Edsjo
c
c                   Bit  Value
c  Bits of hwar:    0 -  1:   Error in FHSetFlags
c                   1 -  2:   Error in FHSetSMPara
c                   2 -  4:   Error in FHSetPara
c                   3 -  8:   Error in FHHiggsCorr
c                   4 -  16:  Error in FHConstraints
c                   5 -  32:  Error in FHSelectUZ
c                   6 -  64:  Error in FHCouplings
c                   7 -  128: Error in FHHiggsProd !Currently not used
c                   8 -  256: Other error (NaN returned)
c---------------------------------------------------------------

      implicit none

c...convert human-readable FeynHiggs indexing
#include "FHCouplings.h"


c...darksusy declarations
      include 'dsmssm.h'
      include 'dsidtag.h'
      integer hwar
      integer jehwar,prlev
      character*12 jehtag
      common /jehiggs/ jehwar,prlev,jehtag

      logical dsisnan
      
c...feynhiggs declarations
      integer mssmpart, fieldren, tanbren, higgsmix, p2approx
      integer looplevel, running_mt, tl_bot_resum,tl_cplx_approx
      real*8 AlfasMZ, ME, MUU, MD, MM, MC, MS, ML, MB, MT, MZ, MW
      real*8 TB, MA0, MHp
      real*8 MSusy
      real*8 M3SL, M2SL, M1SL
      real*8 M3SE, M2SE, M1SE
      real*8 M3SQ, M2SQ, M1SQ
      real*8 M3SU, M2SU, M1SU
      real*8 M3SD, M2SD, M1SD
      complex*16 Atau, At, Ab
      complex*16 Amu, Ac, As
      complex*16 Ae, Au, Ad
      real*8 Qtau, Qt, Qb
      complex*16 MUE, M_1, M_2, M_3
      real*8 scalefactor
      integer error
      real*8 MHiggs(4)
      complex*16 SAeff, UHiggs(3,3), ZHiggs(3,3)
      real*8 gm2
      real*8 deltarho, MWMSSM, MWSM, SW2MSSM, SW2SM
      real*8 edmeTh, edmn, edmHg
      real*8 ccb
      real*8 CKMtheta12, CKMtheta23, CKMtheta13, CKMdelta13
      real*8 invAlfa,GF
      real*8 CKMlam,CKMA,CKMrho,CKMeta
      complex*16 aux
      integer uzint, uzext, mfeff
      complex*16 couplings(ncouplings)
      complex*16 couplingsms(ncouplingsms)
      real*8 gammas(ngammas)
      real*8 gammasms(ngammasms)
      integer excl
c      real*8 prodxs(nprodxs) !Currently not used

c...miscellaneous local declarations
      integer i,j,k,l,m,h,p
      real*8 widthsum

c...set darksusy variables
      hwar=0
      jehtag=idtag
      prlev=0                   ! 0=don't print things, 1=do print things

c...set feynhiggs variables
      mssmpart = 4 ! full MSSM
      fieldren = 0 ! DRbar
      tanbren = 0  ! DRbar
      higgsmix = 2 ! rMSSM, needs to be changed for complex MSSM
      p2approx = 0
      looplevel = 2
      running_mt = 1 ! running mt (updated to FH 2.7.2)
      tl_bot_resum = 1
      tl_cplx_approx = 0 ! rMSSM setting, needs to be changed for complex MSSM
      AlfasMZ = alph3mz
      ME = mass(ke)
      MUU = mass(ku)
      MD = mass(kd)
      MM = mass(kmu)
      MC = mcmc ! check if running, etc.
      MS = mass(ks)
      ML = mass(ktau)
      MB = mbmb ! check if running, etc.
      MT = mass(kt)
      MW = -1 ! mass(kw) ! let FeynHiggs calculate this instead of us giving it
      MZ = mass(kz) 
      invAlfa=1.d0/alphem
      GF=GFermi
      CKMtheta12 = dasin(ckms12)
      CKMtheta23 = dasin(ckms23)
      CKMtheta13 = dasin(ckms13)
      CKMdelta13 = ckmdelta

      CKMlam=ckms12
      CKMA=ckms23/CKMlam**2
      aux=ckms13*exp(dcmplx(0.d0,1.d0)*ckmdelta)/ckmA/ckmlam**3
      CKMrho=dreal(aux)
      CKMeta=dimag(aux)
      

      TB = tanbe
      MA0 = ma
      mhp=-1.d0 ! not used for input, but must be given
      MSusy = -1.d0 ! not used
      M3SL = dsqrt(mass2l(3))
      M3SE = dsqrt(mass2e(3))
      M3SQ = dsqrt(mass2q(3))
      M3SU = dsqrt(mass2u(3))
      M3SD = dsqrt(mass2d(3))
      M2SL = dsqrt(mass2l(2))
      M2SE = dsqrt(mass2e(2))
      M2SQ = dsqrt(mass2q(2))
      M2SU = dsqrt(mass2u(2))
      M2SD = dsqrt(mass2d(2))
      M1SL = dsqrt(mass2l(1))
      M1SE = dsqrt(mass2e(1))
      M1SQ = dsqrt(mass2q(1))
      M1SU = dsqrt(mass2u(1))
      M1SD = dsqrt(mass2d(1))
      At = asoftu(3)
      Ab = asoftd(3)
      Atau = asofte(3)

      Ac = asoftu(2)
      As = asoftd(2)
      Amu = asofte(2)

      Au = asoftu(1)
      Ad = asoftd(1)
      Ae = asofte(1)
      Qtau = 0.d0 ! 0 = use on-shell 
      Qt = 0.d0   ! 0 = use on-shell 
      Qb = 0.d0   ! 0 = use on-shell 
c      Qtau = -1.d0 ! -1 = sqrt(mstop1*mstop2)
c      Qt = -1.d0   ! -1 = sqrt(mstop1*mstop2)
c      Qb = -1.d0   ! -1 = sqrt(mstop1*mstop2)
      MUE = mu
      M_1 = m1
      M_2 = m2
      M_3 = m3
      scalefactor = 1

c...set miscellaneous local variables
      widthsum=0.0d0

c...call feynhiggs
      call FHSetFlags(error, mssmpart, fieldren, tanbren,
     &     higgsmix, p2approx, looplevel, running_mt, tl_bot_resum,
     &     tl_Cplx_Approx)
      if( error .ne. 0 ) then
         write (*,*) 'Error in FeynHiggs FHSetFlags. error=',error
         hwar=ibset(hwar,0)
      endif

      call FHSetSMPara(error,invAlfa,AlfasMZ,GF,
     &  ME,MUU,MD,MM,MC,MS,ML,MB,MW,MZ,CKMlam,CKMA,CKMrho,CKMeta) ! PG updated to FH 2.8.5; JE CHECK: Wolfenstein OK?
      if( error .ne. 0 ) then
         write (*,*) 'Error in FeynHiggs FHSetSMPara. error=',error
         hwar=ibset(hwar,1)
      endif

c      write(*,*) 'XXX: ',
c     &     MT,TB,MA0,MHP,
c     &     M3SL, M3SE, M3SQ, M3SU, M3SD,
c     &     M2SL, M2SE, M2SQ, M2SU, M2SD,
c     &     M1SL, M1SE, M1SQ, M1SU, M1SD,
c     &     MUE,
c     &     Atau, At, Ab,
c     &     Amu, Ac, As,
c     &     Ae, Au, Ad,
c     &     M_1, M_2, M_3,
c     &     Qtau, Qt, Qb

      call FHSetPara(error, scalefactor,
     &     MT,TB,MA0,MHP,
     &     M3SL, M3SE, M3SQ, M3SU, M3SD,
     &     M2SL, M2SE, M2SQ, M2SU, M2SD,
     &     M1SL, M1SE, M1SQ, M1SU, M1SD,
     &     MUE,
     &     Atau, At, Ab,
     &     Amu, Ac, As,
     &     Ae, Au, Ad,
     &     M_1, M_2, M_3,
     &     Qtau, Qt, Qb)
      if( error .ne. 0 ) then
         write (*,*) 'Error in FeynHiggs FHSetPara. error=',error
         hwar=ibset(hwar,2)
      endif

      call FHHiggsCorr(error, MHiggs, SAeff, UHiggs, ZHiggs)

      if( error .ne. 0 ) then
         mass(kh1)=0.d0
         mass(kh2)=0.d0
         mass(kh3)=0.d0
         mass(khc)=0.d0
         hwar=ibset(hwar,3)
         write (*,*) 'Error in FeynHiggs FHHiggsCorr. error=',error
      else
         mass(kh2) = MHiggs(1)    ! H_2 mass (light)
         mass(kh1) = MHiggs(2)    ! H_1 mass (heavy)
         mass(kh3) = MHiggs(3)    ! H30 (A0) mass
         mass(khc) = MHiggs(4)    ! H+- mass
      endif

c...Add extra safecheck
      if (dsisnan(mass(kh1)).or.dsisnan(mass(kh2))
     &  .or.dsisnan(mass(kh3)).or.dsisnan(mass(khc))
     &  .or.mass(kh1).eq.0.d0.or.mass(kh2).eq.0.d0
     &  .or.mass(kh3).eq.0.d0.or.mass(khc).eq.0.d0) then
         write(*,*) 'Error in FeynHiggs NaN returned'
         hwar=ibset(hwar,8)
         return
      endif

c...extract widths and production cross section ratios from FeynHiggs
c...(default, but can be switched off by setting higwid=1 in dsinit)
      if( higwid.eq.5 ) then

c..... no FeynHiggs API below here: we dwelve FeynHiggs guts helped by HiggsBounds
c..... but the HB 3.8.0 example programs HBwithFH are obsolete for FH 2.9.4 because
c..... the 2.9.4 FHSelectUZ requires 4 arguments and no longer 3.

c...following example file in HiggsBounds we set uzint=uzext from here
         uzint=2   !FH-default: uzint=1
         uzext=2   !FH-default: uzext=2
         mfeff=0   !FH-default: mfeff=0
         call FHSelectUZ(error, uzint, uzext, mfeff)
         if ( error .ne. 0) then
            write (*,*) 'Error in FeynHiggs FHSelectUZ. error=',error
            hwar=ibset(hwar,5)
         endif
   
         call FHCouplings(error, couplings, couplingsms,
     &                    gammas, gammasms, 0) ! JE updated to FH 2.7.2; PG to 2.8.3

c...using the following preprocessor macros below 
c...(as defined in FHCouplings.h) 
c... GammaTot(h), Gamma(c), H0HH(h,h1,h2), H0HV(h,hv), H0FF(h,t,g1,g2), H0VV(h,vv),
c... H0SfSf(h,s1,s2,t,g), H0NeuNeu(h,n1,n2), H0ChaCha(h,c1,c2), HpFF(p,g1,g2), 
c... HpHV(hv), HpSfSf(s1,s2,p,g1,g2), HpNeuCha(n1,c2)
         if ( error. ne. 0) then
            write (*,*) 'Error in FeynHiggs FHCouplings. error=',error
            hwar=ibset(hwar,6)
         else
            do i=1,4
               do j=1,32
                  hdwidth(j,i)=0.0d0
               enddo
            enddo
c...total widths
            width(kh2)=GammaTot(1) ! H_2 width (light)            
            width(kh1)=GammaTot(2) ! H_1 width (heavy)   
            width(kh3)=GammaTot(3) ! H30 (A0) width           
            width(khc)=GammaTot(4) ! H+- width           
c...partial widths - neutral higgses
            do i=1,3
               if (i.eq.1) h=2 ! Higgs1 in DS is Higgs2 in FH
               if (i.eq.2) h=1 ! Higgs2 in DS is Higgs1 in FH
               if (i.eq.3) h=3 ! Higgs3 in DS is Higgs3 in FH
               hdwidth(1,i)=  Gamma(H0HH(h,2,2)) !h1 h1
               hdwidth(2,i)=  Gamma(H0HH(h,2,1)) !h1 h2
               hdwidth(3,i)=  Gamma(H0HH(h,1,1)) !h2 h2
               hdwidth(4,i)=  Gamma(H0HH(h,3,3)) !h3 h3
               hdwidth(5,i)=  Gamma(H0HH(h,2,3)) !h1 h3
               hdwidth(6,i)=  Gamma(H0HH(h,1,3)) !h2 h3
               hdwidth(7,i)=  Gamma(H0HH(h,4,4)) !h- h+
               hdwidth(8,i)=  Gamma(H0HV(h,2))   !z h1
               hdwidth(9,i)=  Gamma(H0HV(h,1))   !z h2
               hdwidth(10,i)= Gamma(H0HV(h,3))   !z h3
               hdwidth(11,i)= 0.0d0 !w- h+ and w+ h- (Not existing within FeynHiggs)
               hdwidth(12,i)= Gamma(H0VV(h,3))   !z z
               hdwidth(13,i)= Gamma(H0VV(h,4))   !w+ w-
c...H0FF updated to FynHiggs 2.7.2 below
               hdwidth(14,i)= Gamma(H0FF(h,1,1,1)) !nu_e nu_e-bar
               hdwidth(15,i)= Gamma(H0FF(h,2,1,1)) !e+ e-
               hdwidth(16,i)= Gamma(H0FF(h,1,2,2)) !nu_mu nu_mu-bar
               hdwidth(17,i)= Gamma(H0FF(h,2,2,2)) !mu+ mu-
               hdwidth(18,i)= Gamma(H0FF(h,1,3,3)) !nu_tau nu_tau-bar
               hdwidth(19,i)= Gamma(H0FF(h,2,3,3)) !tau+ tau-
               hdwidth(20,i)= Gamma(H0FF(h,3,1,1)) !u u-bar
               hdwidth(21,i)= Gamma(H0FF(h,4,1,1)) !d d-bar
               hdwidth(22,i)= Gamma(H0FF(h,3,2,2)) !c c-bar
               hdwidth(23,i)= Gamma(H0FF(h,4,2,2)) !s s-bar
               hdwidth(24,i)= Gamma(H0FF(h,3,3,3)) !t t-bar
               hdwidth(25,i)= Gamma(H0FF(h,4,3,3)) !b b-bar
               hdwidth(26,i)= Gamma(H0VV(h,5))   !gluon gluon
               hdwidth(27,i)= 0.0d0 !q q gluon (Not existing within FeynHiggs)
               hdwidth(28,i)= Gamma(H0VV(h,1))   !gamma gamma
               hdwidth(29,i)= Gamma(H0VV(h,2))   !z gamma
               widthsum=0.0d0
               do j=1,2
                  do k=1,2
                     do l=1,4
                        do m=1,3         
                           widthsum=widthsum+
     &                        Gamma(H0SfSf(h,j,k,l,m))
                        enddo
                     enddo
                  enddo
               enddo               
               hdwidth(30,i)= widthsum !sfermions
               widthsum=0.0d0
               do j=1,4
                  do k=1,4
                     widthsum=widthsum+
     &                  Gamma(H0NeuNeu(h,j,k))
                  enddo
               enddo               
               hdwidth(31,i)= widthsum !neutralinos
               widthsum=0.0d0
               do j=1,2
                  do k=1,2
                     widthsum=widthsum+
     &                  Gamma(H0ChaCha(h,j,k))
                  enddo
               enddo               
               hdwidth(32,i)= widthsum !charginos
            enddo
c...partial widths - charged higgs ! JE CHECK: Is ordering for cross-family OK?
            hdwidth(1,4)=  Gamma(HpFF(2,1,1)) !u d-bar
            hdwidth(2,4)=  Gamma(HpFF(2,1,2)) !u s-bar
            hdwidth(3,4)=  Gamma(HpFF(2,1,3)) !u b-bar
            hdwidth(4,4)=  Gamma(HpFF(2,2,1)) !c d-bar
            hdwidth(5,4)=  Gamma(HpFF(2,2,2)) !c s-bar
            hdwidth(6,4)=  Gamma(HpFF(2,2,3)) !c b-bar
            hdwidth(7,4)=  Gamma(HpFF(2,3,1)) !t d-bar
            hdwidth(8,4)=  Gamma(HpFF(2,3,2)) !t s-bar
            hdwidth(9,4)=  Gamma(HpFF(2,3,3)) !t b-bar
            hdwidth(10,4)= Gamma(HpFF(1,1,1)) !nu_e e+
            hdwidth(11,4)= Gamma(HpFF(1,2,2)) !nu_mu mu+
            hdwidth(12,4)= Gamma(HpFF(1,3,3)) !nu_tau tau+
            hdwidth(13,4)= Gamma(HpHV(2))   !w+ h1
            hdwidth(14,4)= Gamma(HpHV(1))   !w+ h2
            hdwidth(15,4)= Gamma(HpHV(3))   !w+ h3
            widthsum=0.0d0
            do j=1,2
               do k=1,2
                  do p=1,2
                     do l=1,2
                        do m=1,3
                           widthsum=widthsum+
     &                          Gamma(HpSfSf(j,k,p,l,m))
                        enddo
                     enddo
                  enddo
               enddo
            enddo               
            hdwidth(20,4)= widthsum !sfermions
            widthsum=0.0d0
            do j=1,4
               do k=1,2
                  widthsum=widthsum+
     &               Gamma(HpNeuCha(j,k))
               enddo
            enddo               
            hdwidth(21,4)= widthsum !neutralinos and charginos

c...LEP production cross section ratios
c...(needed in HiggsBounds, see e.g. HBwithFH.F)
c...(PG: no more needed with HiggsBounds 3.4.0)
c            do i=1,3
c               if (i.eq.1) h=2 ! Higgs1 in DS is Higgs2 in FH
c               if (i.eq.2) h=1 ! Higgs2 in DS is Higgs1 in FH
c               if (i.eq.3) h=3 ! Higgs3 in DS is Higgs3 in FH
c               xsec_lep_hjZ_ratio(i)=
c     &              dreal( Coupling(H0VV(h,3))
c     &                   / CouplingSM(H0VV(h,3)) )**2.0D0
c     &            + dimag( Coupling(H0VV(h,3))
c     &                   / CouplingSM(H0VV(h,3)) )**2.0D0
c            enddo
c            do j=1,3
c               do i=1,3
c                  if (j.eq.1) k=2 ! Higgs1 in DS is Higgs2 in FH
c                  if (j.eq.2) k=1 ! Higgs2 in DS is Higgs1 in FH
c                  if (j.eq.3) k=3 ! Higgs3 in DS is Higgs3 in FH
c                  if (i.eq.1) h=2 ! Higgs1 in DS is Higgs2 in FH
c                  if (i.eq.2) h=1 ! Higgs2 in DS is Higgs1 in FH
c                  if (i.eq.3) h=3 ! Higgs3 in DS is Higgs3 in FH
c                  xsec_lep_hjhi_ratio(j,i)=(
c     &               dreal( Coupling(H0HV(k,h)) )**2.0D0
c     &             + dimag( Coupling(H0HV(k,h)) )**2.0D0 )
c     &             /(GF*sqrt(2.0d0)*MZ**2.0d0)
c               enddo
c            enddo
         endif
 
c         call FHHiggsProd(error, prodxs)
c         if ( error. ne. 0 ) then
c            write (*,*) 'Error in FeynHiggs FHHiggsProd. error=',error
c            hwar=ibset(hwar,7)
c         else
c            !Get Higgs production cross sections here
c            !Currently not needed for anything in DS
c         endif

      endif

      call FHConstraints(error, gm2, deltarho, MWMSSM, MWSM, SW2MSSM, 
     &     SW2SM, edmeTh, edmn, edmHg, ccb) ! PG updated to FH 2.9.4
      if( error .ne. 0 ) then
         alpha=0.d0
         delrho=0.d0
         write (*,*) 'Error in FeynHiggs FHConstraints. error=',error
         hwar=ibset(hwar,4)
      else
         alpha = dasin(dreal(SAeff)) ! alpha
         delrho=deltarho
      endif
      
      end

