**********************************************************************
*** subroutine dstoSLHA transfers DarkSUSY data to the SUSY Les
*** Houches Accord format (into the array slhadata)
*** with the help of SLHALIB by Thomas Hahn.
*** Input: opt = 1 - an SLHA2 file with full 6x6 squark and slepton
***                  mixings and 3x3 trilinear couplings is written.
***                  For minimal flavour violation
***                  this is a bit of an overkill, but more consistent
***                  with how things are stored internally in DarkSUSY
***              2 - a minimal flavour violation SLHA2 file is written
***                  with only third generation mixings and trilinear
***                  couplings. This option will in principle be an
***                  approximation to the internal setup in DarkSUSY,
***                  but most likely good enough for most scenarios.
***              3 - an mSUGRA input file (i.e. only with the SM and
***                  mSUGRA input parameters). This option could be
***                  of interest when the RGEs should be solved with
***                  an external program.
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: 2008-06-24
**********************************************************************
      subroutine dstoSLHA(opt)
      implicit none
c...Note: These include files need to be included with the precompiler
c...directive for include and not the fortran one. Hence, this file
c...has to end on .F instead of .f as well.

#include "SLHA.h"

#include "PDG.h"

      include 'dsmssm.h'
      include 'dsidtag.h'
      include 'dsdirver.h'
#include "dsslha.h"

      real*8 dsmsbarmass
      complex*16 aux
      integer opt
      integer error,g,n,c,n1,n2,c1,c2,s1,s2,g1,g2

c...Clear slhadata
      call SLHAClear(slhadata)

c...Block SPINFO
      if (opt.ne.3) then
         call SLHAPutInfo(BlockSPInfo(1),1,"DarkSUSY")
         call SLHAPutInfo(BlockSPInfo(1),2,dsversion)
         call SLHAPutInfo(BlockSPInfo(1),3,"Model: "//idtag)
      endif


C...Block MODSEL
c      ModSel_Model=...
      ModSel_Model=0 ! 0=MSSM
      if (modeltype.eq.3) then
         ModSel_Model=1 ! 1=mSUGRA
      endif
      if (opt.ne.3) then
         ModSel_RPV = 0 ! no R-partity violation
         ModSel_CPV = 1 ! 1=Only CKM-CP violation
      endif
c...For the flavour violation flag, we set the flag to 3, even in the MFV
c...scenario. The reason for this is that we want to give the full
c...6x6 slepton and squark mixings below
      if (opt.eq.2) then
         ModSel_FV = 0  ! 0=no flavour violation, 6x6 SQMIX etc should be absent
      elseif (opt.eq.1) then
         ModSel_FV = 3  ! 3=flavour violation in quarks and leptons
      elseif (opt.eq.3) then
         continue ! do nothing in this case
      else
         write(*,*) 'DS Error in dstoslha: invalid option: ',opt
         ModSel_FV = 3  ! 3=flavour violation in quarks and leptons
         opt=1
      endif


c...Block SMINPUTS
      SMInputs_invAlfaMZ=1.d0/alphem
      SMInputs_GF=GFermi
      SMInputs_AlfasMZ = alph3mz
      SMInputs_MZ = mass(kz)

c...JE FIX: running masses at 2 GeV
      SMInputs_Mf(3,3) = mass(kt)
      SMInputs_Mf(4,3) = mbmb
      SMInputs_Mf(2,3) = mass(ktau)
      if (opt.ne.3) then ! add more stuff
        SMInputs_Mf(2,1) = mass(ke)
        SMInputs_Mf(2,2) = mass(kmu)
        SMInputs_Mf(3,1) = mass(ku)
        SMInputs_Mf(3,2) = mcmc
        SMInputs_Mf(4,1) = mass(kd)
        SMInputs_Mf(4,2) = mass(ks)
      endif


c...Give MINPAR block if we started with mSUGRA
      if (modeltype.eq.3) then
         MinPar_M0=m0var
         MinPar_M12=mhfvar
         MinPar_A=a0var
         MinPar_signMUE=sgnmuvar
         MinPar_TB=tgbetavar
      endif

      if (opt.eq.3) return ! no extra blocks needed for mSUGRA input SLHA file

c...Block EXTPAR
c      ExtPar_Q = -1 ! this gives scale sqrt(m_stop1*m_stop2), OK?
      ExtPar_Q = mass(kz) ! Input scale is mZ.
      ExtPar_M1 = M1
      ExtPar_M2 = M2
      ExtPar_M3 = M3
      ExtPar_Atau=asofte(3)
      ExtPar_At=asoftu(3)
      ExtPar_Ab=asoftd(3)
      ExtPar_MUE = mu
c      ExtPar_MA02 = ...  ! this is tree-level mass, not pole mass. 
      ExtPar_TB = tanbe
      ExtPar_MA0=mass(kh3) ! this is pole mass
      do g=1,3
          ExtPar_MSL(g)=sqrt(mass2l(g))
          ExtPar_MSE(g)=sqrt(mass2e(g))
          ExtPar_MSQ(g)=sqrt(mass2q(g))
          ExtPar_MSU(g)=sqrt(mass2u(g))
          ExtPar_MSD(g)=sqrt(mass2d(g))
      enddo

c...Block MSoft
      MSoft_Q=mass(kz) ! Input scale is mz
      MSoft_M1 = M1
      MSoft_M2 = M2
      MSoft_M3 = M3
      do g=1,3
         MSoft_MSL(g)=sqrt(mass2l(g))
         MSoft_MSE(g)=sqrt(mass2e(g))
         MSoft_MSQ(g)=sqrt(mass2q(g))
         MSoft_MSU(g)=sqrt(mass2u(g))
         MSoft_MSD(g)=sqrt(mass2d(g))
      enddo

c...Block MASS
      do g=1,3
         Mass_Mf(1,g) = mass(knu(g))
         Mass_Mf(2,g) = mass(kl(g))
         Mass_Mf(3,g) = mass(kqu(g))
         Mass_Mf(4,g) = mass(kqd(g))
         Mass_MSf(1,1,g) = mass(ksnu(g))
c         Mass_MSf(2,1,g) = mass(ksnu(g+3)) ! No sneutrino_2's
         Mass_MSf(1,2,g) = mass(ksl(g))   
         Mass_MSf(2,2,g) = mass(ksl(g+3))
         Mass_MSf(1,3,g) = mass(ksqu(g))
         Mass_MSf(2,3,g) = mass(ksqu(g+3))
         Mass_MSf(1,4,g) = mass(ksqd(g))
         Mass_MSf(2,4,g) = mass(ksqd(g+3))
      enddo
      Mass_MZ = mass(kz)
      Mass_MW = mass(kw)
      Mass_Mh0 = mass(kh2)
      Mass_MHH = mass(kh1)
      Mass_MA0 = mass(kh3)
      Mass_MHp = mass(khc)

c...Neutralinos
      do n=1,4
         if (abs(dimag(NeuNmx(n,1))).gt.0.d0) then ! complex
            Mass_MNeu(n)=-mass(kn(n))
         else ! real
            Mass_MNeu(n)=mass(kn(n))
         endif
      enddo

c...Charginos
      do c=1,2
         Mass_MCha(c)=mass(kcha(c))
      enddo

c...Gluino
      Mass_MGl=mass(kgluin)

c...Gravitino
c      Mass_MGrav=...  ! We don't include the gravitino

c...Block NMIX - neutralino mixing
      do n=1,4
         if (abs(dimag(NeuNmx(n,1))).gt.0.d0) then ! complex, pick complex part
            do n2=1,4
               NMix_ZNeu(n,n2)=dcmplx(0.d0,-1.d0)*NeuNmx(n,n2)
            enddo
         else ! pick real part
            do n2=1,4
               NMix_ZNeu(n,n2)=NeuNmx(n,n2)
            enddo
         endif
      enddo

c...Block UMIX - U matrix for chargino mixing
      do c=1,2
         do c2=1,2
            UMix_UCha(c,c2)=chaumx(c,c2)
         enddo
      enddo

c...Block VMIX - V matrix for chargino mixing
      do c=1,2
         do c2=1,2
            VMix_VCha(c,c2)=chavmx(c,c2)
         enddo
      enddo

c...BLOCK alpha
      Alpha_Alpha = alpha

      if (opt.eq.2) then ! simplified MFV output

c...Block STAUMIX
         StauMix_USf(1,1)=sldlmx(3,3) 
         StauMix_USf(1,2)=sldrmx(3,3) 
         StauMix_USf(2,1)=sldlmx(6,3) 
         StauMix_USf(2,2)=sldrmx(6,3)

c...Block STOPMIX
         StopMix_USf(1,1)=squlmx(3,3) 
         StopMix_USf(1,2)=squrmx(3,3) 
         StopMix_USf(2,1)=squlmx(6,3) 
         StopMix_USf(2,2)=squrmx(6,3)

c...Block SBOTMIX
         SbotMix_USf(1,1)=sqdlmx(3,3) 
         SbotMix_USf(1,2)=sqdrmx(3,3) 
         SbotMix_USf(2,1)=sqdlmx(6,3) 
         SbotMix_USf(2,2)=sqdrmx(6,3)

c...NOTE: For the blocks AE, AU and AD, these things are 
c...in principle defined in the block ExtPar. We here give
c...them anyway for consistency. For the full thing
c...we should use the blocks TE, TU and TD below instead.
c...Block AE
c...Block AU
c...Block AD
         Ae_Q = mass(kz)
         Au_Q = mass(kz)
         Ad_Q = mass(kz)

cc      do g1=1,3
cc         do g2=1,3
cc            Ae_Af(g1,g2)=dcmplx(0.0d0,0.0d0)
cc            Au_Af(g1,g2)=dcmplx(0.0d0,0.0d0)
cc            Ad_Af(g1,g2)=dcmplx(0.0d0,0.0d0)
cc         enddo
cc         Ae_Af(g1,g1)=asofte(g1)
cc         Au_Af(g1,g1)=asoftu(g1)
cc         Ad_Af(g1,g1)=asoftd(g1)
cc      enddo
cc...Note: only (3,3) component should be given, hence do not use above      
         Ae_Af(3,3)=asofte(3)
         Au_Af(3,3)=asoftu(3)
         Ad_Af(3,3)=asoftd(3)

      else ! opt = 1 , full 6x6

c...Block TE
c...Block TU
c...Block TD
c...Missing Yukawas corrected 080704
         Te_Q = mass(kz)
         Tu_Q = mass(kz)
         Td_Q = mass(kz)

         do g1=1,3
            do g2=1,3
               Te_Tf(g1,g2)=dcmplx(0.0d0,0.0d0)
               Tu_Tf(g1,g2)=dcmplx(0.0d0,0.0d0)
               Td_Tf(g1,g2)=dcmplx(0.0d0,0.0d0)
            enddo
            Te_Tf(g1,g1)=asofte(g1)*yukawa(kl(g1))
            Tu_Tf(g1,g1)=asoftu(g1)*yukawa(kqu(g1))
            Td_Tf(g1,g1)=asoftd(g1)*yukawa(kqd(g1))
         enddo


c...Block SNMIX - Sneutrino mixing matrices
         do s1=1,3
            do s2=1,3
               SnuMix_UASf(s1,s2) = slulmx(s1,s2)
            enddo
         enddo

c...Block SLMIX - Charged slepton mixing matrices
         do s1=1,6
            do s2=1,3
               SelMix_UASf(s1,s2)   = sldlmx(s1,s2)
               SelMix_UASf(s1,s2+3) = sldrmx(s1,s2)
            enddo
         enddo

c...Block USQMIX - Up squark mixing matrices
         do s1=1,6
            do s2=1,3
               USqMix_UASf(s1,s2)   = squlmx(s1,s2)
               USqMix_UASf(s1,s2+3) = squrmx(s1,s2)
            enddo
         enddo

c...Block DSQMIX - Down squark mixing matrices
         do s1=1,6
            do s2=1,3
               DSqMix_UASf(s1,s2)   = sqdlmx(s1,s2)
               DSqMix_UASf(s1,s2+3) = sqdrmx(s1,s2)
            enddo
         enddo

      endif ! end of opt=1 (full 6x6)

c...Block YE
c...Block YU
c...Block YD
         Ye_Q = mass(kz)
         Yu_Q = mass(kz)
         Yd_Q = mass(kz)
c      do g1=1,3
c         do g2=1,3
c            Ye_Yf(g1,g2)=dcmplx(0.0d0,0.0d0)
c            Yu_Yf(g1,g2)=dcmplx(0.0d0,0.0d0)
c            Yd_Yf(g1,g2)=dcmplx(0.0d0,0.0d0)
c         enddo
c         Ye_Yf(g1,g1)=yukawa(kl(g1))
c         Yu_Yf(g1,g1)=yukawa(kqu(g1))
c         Yd_Yf(g1,g1)=yukawa(kqd(g1))
c      enddo
c...Note: only (3,3) component should be given, hence do not use the above      
      Ye_Yf(3,3)=yukawa(kl(3))
      Yu_Yf(3,3)=yukawa(kqu(3))
      Yd_Yf(3,3)=yukawa(kqd(3))

    

c...Block MSL2
c...Block MSE2
c...Block MSQ2
c...Block MSU2
c...Block MSD2
      MSL2_Q=mass(kz)
      MSE2_Q=mass(kz)
      MSQ2_Q=mass(kz)
      MSU2_Q=mass(kz)
      MSD2_Q=mass(kz)
      do g1=1,3
         do g2=1,3
            MSL2_MSL2(g1,g2)=dcmplx(0.0d0,0.0d0)
            MSE2_MSE2(g1,g2)=dcmplx(0.0d0,0.0d0)
            MSQ2_MSQ2(g1,g2)=dcmplx(0.0d0,0.0d0)
            MSU2_MSU2(g1,g2)=dcmplx(0.0d0,0.0d0)
            MSD2_MSD2(g1,g2)=dcmplx(0.0d0,0.0d0)
         enddo
         MSL2_MSL2(g1,g1)=mass2l(g1)
         MSE2_MSE2(g1,g1)=mass2e(g1)
         MSQ2_MSQ2(g1,g1)=mass2q(g1)
         MSU2_MSU2(g1,g1)=mass2u(g1)
         MSD2_MSD2(g1,g1)=mass2d(g1)
      enddo
      

c...Block VCKMINPUTS - CKM mixing - don't use inputs like this
c      VCKMIn_lambda = ckms12
c      VCKMIn_A=ckms23/VCKMIn_lambda**2
c      aux=ckms13*exp(dcmplx(0.d0,1.d0)*ckmdelta)
c     &   /VCKMIn_A/VCKMIn_lambda**3
c      VCKMIn_rho=dreal(aux)
c      VCKMIn_eta=dimag(aux)
 

c...Block VCKM
      VCKM_Q = mass(kz)
      do g1=1,3
         do g2=1,3
            VCKM_VCKM(g1,g2)=ckm(g1,g2)
         enddo
      enddo

c...Now give some extra auxiliary blocks
c...Block HMIX
      HMix_Q = mass(kz)
      HMix_MUE = mu
      HMix_TB = tanbe
c...Note, last sqrt(2) comes from going from Gunion and Haber convention
c...to SLHA convention.
      HMix_VEV = sqrt(2.0d0)*mass(kw)/g2weak * sqrt(2.0d0) 
c...JE FIX, we here give pole mass, but should give tree-level mass
      HMix_MA02 = mass(kh3)**2

c...Block GAUGE
      Gauge_Q = mass(kz)
      Gauge_g1 = gywmz
      Gauge_g2 = g2wmz
      Gauge_g3 = g3stro

      return
      end

