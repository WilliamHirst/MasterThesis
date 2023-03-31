**********************************************************************
*** subroutine dsfromSLHA transfers SLHA data in the array slhadata
*** to DarkSUSY. The SLHA things are read with with the help of
*** SLHALIB by Thomas Hahn.
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: 2008-07-01
*** TO FIX: Better error handling
***         More options on how to treat redundant information, e.g.
***         should Yukawas from SLHA file be used instead of letting
***         DarkSUSY run them (as is currently done).
**********************************************************************
      subroutine dsfromSLHA
      implicit none
c...Note: These include files need to be included with the precompiler
c...directive for include and not the fortran one. Hence, this file
c...has to end on .F instead of .f as well.


#include "SLHA.h"

#include "PDG.h"

      include 'dsmssm.h'
      include 'dsidtag.h'
#include "dsslha.h"

      integer plev
      parameter(plev=0) ! 0=only errors, 1=errors and warnings
      integer unphys,hwarning,i,opt
      integer error,g,n,c,n1,n2,c1,c2,s1,s2,g1,g2
      real*8 dsmsbarmass,rtmp,dsgf2s2thw,tmp,dsabsq,dszarg,
     &  dsmqpole4loop
      complex*16 aux
      character*80 model
c...mSUGRA parameters
      real*8 m0,mhf,a0,sgnmu,tgbeta

c...Zero things here
      call dsmssmzero

c...If this model comes from DarkSUSY, try extracting the model idtag
c      call SLHACopyString(SPInfo_Len,SPInfo_Text(1,1),model)
c      write(*,*) 'Model: ',model
c      write(*,*) SPInfo_NLines

c...Block MODSEL
      if (ModSel_Model.ne.invalid) then
      if (ModSel_Model.eq.0) then !
         modeltype=0 ! full MSSM
      elseif (ModSel_Model.eq.1) then !
         modeltype=3 ! mSUGRA
      else
         write(*,*) 
     &  'DS Error: Can currently not read an SLHA2 file with',
     &  ' ModSel_Model = ',ModSel_Model
         stop
      endif
      endif

c...Set higwid parameter
      higwid=1

      if (ModSel_RPV.ne.invalid) then
      if (ModSel_RPV.ne.0) then
         write(*,*) 
     &  'DS Error: Can currently not read an SLHA2 file with',
     &  ' ModSel_RPV = ',ModSel_RPV
         stop
      endif
      endif

      if (ModSel_CPV.ne.invalid) then
      if (ModSel_CPV.ne.1) then
         write(*,*) ('-',i=1,70)
         write(*,*) 
     &  'DS Error: From the SLHA2 file I was expecting',
     &  ' ModSel_CPV = 1, but found ModSel_CPV = ',ModSel_CPV
         write(*,*) ('-',i=1,70)
         stop
      endif
      endif
      
c...ModSelFV
      opt=2 ! no flavor violation (2x2 sfermion matices...)
      if (ModSel_FV.eq.3) then
         opt=1 ! flavor violation (6x6 sfermion matrices...)
      endif

c...Block SMINPUTS
      if (SMInputs_invAlfaMZ.ne.invalid) then
         alphem=1.d0/SMInputs_invAlfaMZ
      endif
      
      call dssetfromslha(SMInputs_GF,GFermi,0,'')

      call dssetfromslha(SMInputs_AlfasMZ,alph3mz,0,'')
      call dssetfromslha(SMInputs_MZ,mass(kz),0,'')

      if (SMInputs_Mf(2,1).ne.invalid) then
         mass(ke)=SMInputs_Mf(2,1)
      endif
      if (SMInputs_Mf(2,2).ne.invalid) then
         mass(kmu)=SMInputs_Mf(2,2)
      endif
      if (SMInputs_Mf(2,3).ne.invalid) then
         mass(ktau)=SMInputs_Mf(2,3)
      endif

      if (SMInputs_Mf(3,1).ne.invalid) then
         mu2gev=SMInputs_Mf(3,1)
      endif

      if (SMInputs_Mf(3,2).ne.invalid) then
         mcmc=SMInputs_Mf(3,2)
      endif
      
      if (SMInputs_Mf(3,3).ne.invalid) then
         mass(kt)=SMInputs_Mf(3,3)
      endif

      if (SMInputs_Mf(4,1).ne.invalid) then
         md2gev=SMInputs_Mf(4,1)
      endif

      if (SMInputs_Mf(4,2).ne.invalid) then
         ms2gev=SMInputs_Mf(4,2)
      endif

      if (SMInputs_Mf(4,3).ne.invalid) then
         mbmb=SMInputs_Mf(4,3)
      endif

c...Calculate sin^2(theta_w) at Z scale
      s2thw=dsgf2s2thw(GFermi,alphem,mass(kz),mass(kt),1)

c mass(ku), mass(kd), mass(ks) are running masses at 2 GeV
c mass(kc), mass(kb), mass(kt) are pole masses
      call dsfindmtmt ! must come before dsmsbarmass
      mass(ku)     =  mu2gev
      mass(kd)     =  md2gev
      mass(ks)     =  ms2gev
      mass(kc)     =  dsmqpole4loop(kc,mcmc)
      mass(kb)     =  dsmqpole4loop(kb,mbmb)


c...Block VCKMINPUTS - CKM mixing
c...We will first read in CKM inputs (if present) and set up the CKM
c...matrix. If the block VCKM is present, it will override below
      call dssetfromslha(VCKMIn_lambda,ckms12,0,'')
      if (VCKMIn_A.ne.invalid) then
         ckms23=VCKMIn_A*ckms12**2
      endif
      if (VCKMIn_rhobar.ne.invalid) then
         if (VCKMIn_etabar.ne.invalid) then
            tmp=VCKMIn_etabar
         else
            tmp=ckmdelta ! DS default
         endif
         aux=VCKMIn_rhobar+dcmplx(0.d0,tmp)
         aux=aux*(ckms23/ckms12**2)*ckms12**3
         ckms13=dsabsq(aux)
         ckmdelta=dszarg(aux)
      endif

      call dssuconst_ckm ! set up CKM matrix from these inputs

c...Now replace what we have just calculated, in case the CKM block
c...is present
c...Block VCKM
      do g1=1,3
         do g2=1,3
            call dssetfromslhac(VCKM_VCKM(g1,g2),
     &         ckm(g1,g2),0,'')
         enddo
      enddo

c...Block MINPAR
c...Read mSUGRA parameters if mSUGRA
      if (modeltype.eq.3) then ! mSUGRA
      if (MinPar_M0.ne.invalid) then
         m0=MinPar_M0
         m0var=m0
      else
         write(*,*) 'DS ERROR: No m0 parameter given for',
     &      ' mSUGRA SLHA file'
         stop
      endif

      if (MinPar_M12.ne.invalid) then
         mhf=MinPar_M12
         mhfvar=mhf
      else
         write(*,*) 'DS ERROR: No mhf parameter given for',
     &      ' mSUGRA SLHA file'
         stop
      endif

      if (MinPar_TB.ne.invalid) then
         tgbeta=MinPar_TB
         tgbetavar=tgbeta
      else
         write(*,*) 'DS ERROR: No tanbeta parameter given for',
     &      ' mSUGRA SLHA file'
         stop
      endif

      if (MinPar_A.ne.invalid) then
         a0=MinPar_A
         a0var=a0
      else
         write(*,*) 'DS ERROR: No A parameter given for',
     &      ' mSUGRA SLHA file'
         stop
      endif

      if (MinPar_signMUE.ne.invalid) then
         sgnmu=MinPar_signMUE
         sgnmuvar=sgnmu
      else
         write(*,*) 'DS ERROR: No sign(mu) parameter given for',
     &      ' mSUGRA SLHA file'
         stop
      endif
      endif

c...Now mSUGRA parameters are read, if no low-energy blocks are present
c...we will regard this as an INPUT file and run ISASUGRA.
c...If low-energy blocks are present, we will regard this as an
c...OUTPUT file from a spectrum calculator and use the values as given
      if (modeltype.eq.3.and.Mass_MNeu(1).eq.invalid) then
c... regard as input file
         write(*,*) 'SLHA mSUGRA file regarded as input file,',
     &     ' will run ISASUGRA.'
         call dsgive_model_isasugra(m0,mhf,a0,sgnmu,tgbeta)
         call dsrge_isasugra(unphys,hwarning)
         call dsmodelsetup_isasugra
         return
      else
         write(*,*) 'SLHA mSUGRA file regarded as output file,',
     &     ' will read low-energy values from file.'
      endif


c...Block EXTPAR, MSoft and HMIX
      if ((ExtPar_Q.eq.invalid).or.
     &   (dsabsq(ExtPar_Q-mass(kz)).gt.2.0d0)) then 
         if (prl.ge.1) then
         write(*,*) ('-',i=1,70)
         write(*,*)
     &  'DS Warning: DarkSUSY wants to have the values in the EXTPAR'
         write(*,*)
     &  'block at the electroweak scale (MZ), but they are given'
         write(*,*)
     &  'at the scale ',ExtPar_Q,' instead. DarkSUSY will not run'
         write(*,*)
     &  'the EXTPAR parameters but use them at face value anyway'
         write(*,*) ('-',i=1,70)
         endif
      endif

      if ((MSoft_Q.eq.invalid).or.
     &   (dsabsq(MSoft_Q-mass(kz)).gt.2.0d0)) then 
         if (prl.ge.1) then
         write(*,*) ('-',i=1,70)
         write(*,*)
     &  'DS Warning: DarkSUSY wants to have the values in the MSoft'
         write(*,*)
     &  'block at the electroweak scale (MZ), but they are given'
         write(*,*)
     &  'at the scale ',MSoft_Q,' instead. DarkSUSY will not run'
         write(*,*)
     &  'the MSoft parameters but use them at face value anyway.'
         write(*,*) ('-',i=1,70)
         endif
      endif

      
      call dssetfromslha2(MSoft_M1,ExtPar_M1,M1,2,'M1')
      call dssetfromslha2(MSoft_M2,ExtPar_M2,M2,2,'M2')
      call dssetfromslha2(MSoft_M3,ExtPar_M3,M3,2,'M3')

      call dssetfromslha(ExtPar_Atau,asofte(3),0,'')
      call dssetfromslha(ExtPar_At,asoftu(3),0,'')
      call dssetfromslha(ExtPar_Ab,asoftd(3),0,'')

      call dssetfromslha2(HMix_MUE,ExtPar_MUE,mu,2,'MU')
      call dssetfromslha2(HMix_TB,ExtPar_TB,tanbe,2,'TanBe')
      call dssetfromslha2(ExtPar_MA0,Mass_MA0,ma,2,'MA0')
      mass(kh3)=ma

      do g=1,3
         if (MSoft_MSL(g).ne.invalid) then
            mass2l(g)=MSoft_MSL(g)**2
         else
            if (ExtPar_MSL(g).ne.invalid) then
               mass2l(g)=ExtPar_MSL(g)**2
            else
             write(*,*) 
     &       'DS Error: MSL(',g,') is not given in the SLHA file'
             write(*,*) 
     &       'but it is needed.'
             stop
            endif
         endif

         if (MSoft_MSE(g).ne.invalid) then
            mass2e(g)=MSoft_MSE(g)**2
         else
            if (ExtPar_MSE(g).ne.invalid) then
               mass2e(g)=ExtPar_MSE(g)**2
            else
             write(*,*) 
     &       'DS Error: MSE(',g,') is not given in the SLHA file'
             write(*,*) 
     &       'but it is needed.'
             stop
            endif
         endif

         if (MSoft_MSQ(g).ne.invalid) then
            mass2q(g)=MSoft_MSQ(g)**2
         else
            if (ExtPar_MSQ(g).ne.invalid) then
               mass2q(g)=ExtPar_MSQ(g)**2
            else
             write(*,*) 
     &       'DS Error: MSQ(',g,') is not given in the SLHA file'
             write(*,*) 
     &       'but it is needed.'
             stop
            endif
         endif

         if (MSoft_MSU(g).ne.invalid) then
            mass2u(g)=MSoft_MSU(g)**2
         else
            if (ExtPar_MSU(g).ne.invalid) then
               mass2u(g)=ExtPar_MSU(g)**2
            else
             write(*,*) 
     &       'DS Error: MSU(',g,') is not given in the SLHA file'
             write(*,*) 
     &       'but it is needed.'
             stop
            endif
         endif

         if (MSoft_MSD(g).ne.invalid) then
            mass2d(g)=MSoft_MSD(g)**2
         else
            if (ExtPar_MSD(g).ne.invalid) then
               mass2d(g)=ExtPar_MSD(g)**2
            else
             write(*,*) 
     &       'DS Error: MSD(',g,') is not given in the SLHA file'
             write(*,*) 
     &       'but it is needed.'
             stop
            endif
         endif

      enddo

c...Now set up things (dssusy) to have defaults for the SUSY sector
c...The things calculated here, will presumably be superseeded below

      call dssuconst
      call dsspectrum(unphys,hwarning)
      

c...Block MASS
      do g=1,3
         call dssetfromslha(Mass_Mf(1,g),mass(knu(g)),0,'')
         call dssetfromslha(Mass_Mf(2,g),mass(kl(g)),0,'')
         call dssetfromslha(Mass_Mf(3,g),mass(kqu(g)),0,'')
         call dssetfromslha(Mass_Mf(4,g),mass(kqd(g)),0,'')

c...For the sparticles, we rely on things being calculated from
c...the EW scale parameters above. If they are not, the values in
c...the MASS block will supersede those values. Hence, only replace
c...if the values are given in the SLHA file.

         call dssetfromslha(Mass_MSf(1,1,g),mass(ksnu(g)),0,'')
c         call dsetfromslha(Mass_MSf(2,1,g),mass(ksnu(g+3)),0,'') ! No sneutrino_2's
         call dssetfromslha(Mass_MSf(1,2,g),mass(ksl(g)),0,'')
         call dssetfromslha(Mass_MSf(2,2,g),mass(ksl(g+3)),0,'')
         call dssetfromslha(Mass_MSf(1,3,g),mass(ksqu(g)),0,'')
         call dssetfromslha(Mass_MSf(2,3,g),mass(ksqu(g+3)),0,'')
         call dssetfromslha(Mass_MSf(1,4,g),mass(ksqd(g)),0,'')
         call dssetfromslha(Mass_MSf(2,4,g),mass(ksqd(g+3)),0,'')
      enddo
      call dssetfromslha(Mass_MZ,mass(kz),0,'')
      if (Mass_MW.ne.invalid) then
c...OK, we now have to enforce the tree-level condidtion for unitarity
c...We then have a choice of calculating both sin^2 theta_W and MW
c...from alpha,MZ and GF as we normally do in DarkSUSY. The line below
c...would enforce that.
c      mass(kw)=mass(kz)*sqrt(1.d0-s2thw)
c...however, it is more prudent to take the value of MW from the SLHA file
c...if given, and instead enforce the tree-level condition by redefining
c...sin^2 theta_W. That we do here:
         mass(kw)=Mass_MW
         s2thw=1.d0-mass(kw)**2/mass(kz)**2
      else
c...In case MW is not given in the SLHA file, we calcualte it from
c...the input parameters.
c...We thus enforce the tree-level condition for unitarity
         mass(kw)=mass(kz)*sqrt(1.d0-s2thw)
      endif

      call dssetfromslha(Mass_MW,mass(kw),0,'')
      call dssetfromslha(Mass_Mh0,mass(kh2),0,'')
      call dssetfromslha(Mass_MHH,mass(kh1),0,'')
      call dssetfromslha(Mass_MA0,mass(kh3),0,'')
      call dssetfromslha(Mass_MHp,mass(khc),0,'')



c...Neutralinos
      do n=1,4
         call dssetfromslha(Mass_MNeu(n),mass(kn(n)),0,'')
      enddo

c...Charginos
      do c=1,2
         call dssetfromslha(Mass_MCha(c),mass(kcha(c)),0,'')
      enddo

c...Gluino
      call dssetfromslha(Mass_MGl,mass(kgluin),0,'')

c...Gravitino
c      Mass_MGrav=...  ! We don't include the gravitino

c...Block NMIX - neutralino mixing
      do n=1,4
         do n2=1,4
             call dssetfromslhac(NMix_ZNeu(n,n2),NeuNmx(n,n2),0,'')
         enddo
      enddo

c...Now rewrite NMix in the complex form we want it, i.e. use complex
c...eigenvectors and positive mass eigenvalues instead of real eigenvectors
c...and signed eigenvalues
      do n=1,4
         do n2=1,4
            if (mass(kn(n)).lt.0.d0) then
               NeuNmx(n,n2)=dcmplx(0.d0,1.d0)*NeuNmx(n,n2)
            endif
          enddo
          mass(kn(n))=abs(mass(kn(n)))
       enddo

c...Block UMIX - U matrix for chargino mixing
      do c=1,2
         do c2=1,2
            call dssetfromslhac(UMix_UCha(c,c2),chaumx(c,c2),0,'')
         enddo
      enddo

c...Block VMIX - V matrix for chargino mixing
      do c=1,2
         do c2=1,2
            call dssetfromslhac(VMix_VCha(c,c2),chavmx(c,c2),0,'')
         enddo
      enddo

c...BLOCK alpha
      call dssetfromslha(Alpha_Alpha,alpha,0,'')


c...Block YE
c...Block YU
c...Block YD
c...First set up defaults from given masses, then override if these
c...blocks are present in the SLHA file
      call dssuconst_yukawa_running ! calculate yukawas from masses

      if (plev.gt.0.and.Ye_Q.ne.invalid.and.
     &   dsabsq(Ye_Q-mass(kz)).gt.2.d0) then
         if (prl.ge.1) then
         write(*,*) ('-',i=1,70)
         write(*,*)
     &   'DS Warning: SLHA scale for block YE is not MZ as expected,'
         write(*,*)
     &    'it is instead ',Ye_Q,'. Will use this anyway'
         write(*,*) ('-',i=1,70)
         endif
      endif

      if (plev.gt.0.and.Yu_Q.ne.invalid.and.
     &   dsabsq(Yu_Q-mass(kz)).gt.2.d0) then
         if (prl.ge.1) then
         write(*,*) ('-',i=1,70)
         write(*,*)
     &   'DS Warning: SLHA scale for block YU is not MZ as expected,'
         write(*,*)
     &    'it is instead ',Yu_Q,'. Will use this anyway'
         write(*,*) ('-',i=1,70)
         endif
      endif

      if (plev.gt.0.and.Yd_Q.ne.invalid.and.
     &   dsabsq(Yd_Q-mass(kz)).gt.2.d0) then
         if (prl.ge.1) then
         write(*,*) ('-',i=1,70)
         write(*,*)
     &   'DS Warning: SLHA scale for block YD is not MZ as expected,'
         write(*,*)
     &    'it is instead ',Yd_Q,'. Will use this anyway'
         write(*,*) ('-',i=1,70)
         endif
      endif
      
c...JE FIX: Note: the following lines will mean that we instead of
c...using running Yukawas from DarkSUSY will use whatever Yukawas that
c...are specified in the SLHA file. In principle it is more consistent
c...to use the Yukawas from the file, but they might not be at the 
c...scale we want them.
      do g1=1,3
         call dssetfromslha(Ye_Yf(g1,g1),yukawa(kl(g1)),0,'')
         call dssetfromslha(Yu_Yf(g1,g1),yukawa(kqu(g1)),0,'')
         call dssetfromslha(Yd_Yf(g1,g1),yukawa(kqd(g1)),0,'')
      enddo


c...Block MSL2
c...Block MSE2
c...Block MSQ2
c...Block MSU2
c...Block MSD2
      if (plev.gt.0.and.MSL2_Q.ne.invalid.and.
     &   dsabsq(MSL2_Q-mass(kz)).gt.2.d0) then
         if (prl.ge.1) then
         write(*,*) ('-',i=1,70)
         write(*,*)
     &   'DS Warning: SLHA scale for block MSL2 is not MZ as expected,'
         write(*,*)
     &    'it is instead ',MSL2_Q,'. Will use this anyway'
         write(*,*) ('-',i=1,70)
         endif
      endif

      if (plev.gt.0.and.MSE2_Q.ne.invalid.and.
     &   dsabsq(MSE2_Q-mass(kz)).gt.2.d0) then
         if (prl.ge.1) then
         write(*,*) ('-',i=1,70)
         write(*,*)
     &   'DS Warning: SLHA scale for block MSE2 is not MZ as expected,'
         write(*,*)
     &    'it is instead ',MSE2_Q,'. Will use this anyway'
         write(*,*) ('-',i=1,70)
         endif
      endif

      if (plev.gt.0.and.MSQ2_Q.ne.invalid.and.
     &   dsabsq(MSQ2_Q-mass(kz)).gt.2.d0) then
         if (prl.ge.1) then
         write(*,*) ('-',i=1,70)
         write(*,*)
     &   'DS Warning: SLHA scale for block MSQ2 is not MZ as expected,'
         write(*,*)
     &    'it is instead ',MSQ2_Q,'. Will use this anyway'
         write(*,*) ('-',i=1,70)
         endif
      endif

      if (plev.gt.0.and.MSU2_Q.ne.invalid.and.
     &   dsabsq(MSU2_Q-mass(kz)).gt.2.d0) then
         if (prl.ge.1) then
         write(*,*) ('-',i=1,70)
         write(*,*)
     &   'DS Warning: SLHA scale for block MSU2 is not MZ as expected,'
         write(*,*)
     &    'it is instead ',MSU2_Q,'. Will use this anyway'
         write(*,*) ('-',i=1,70)
         endif
      endif

      if (plev.gt.0.and.MSD2_Q.ne.invalid.and.
     &   abs(MSD2_Q-mass(kz)).gt.2.d0) then
         if (prl.ge.1) then
         write(*,*) ('-',i=1,70)
         write(*,*)
     &   'DS Warning: SLHA scale for block MSD2 is not MZ as expected,'
         write(*,*)
     &    'it is instead ',MSD2_Q,'. Will use this anyway'
         write(*,*) ('-',i=1,70)
         endif
      endif

      do g1=1,3
         call dssetfromslha(MSL2_MSL2(g1,g1),mass2l(g1),0,'')
         call dssetfromslha(MSE2_MSE2(g1,g1),mass2e(g1),0,'')
         call dssetfromslha(MSQ2_MSQ2(g1,g1),mass2q(g1),0,'')
         call dssetfromslha(MSU2_MSU2(g1,g1),mass2u(g1),0,'')
         call dssetfromslha(MSD2_MSD2(g1,g1),mass2d(g1),0,'')
      enddo


c...Now read mixings, first in minimal flavour violation setup (opt=2)
      if (opt.eq.2) then
c...Set first two generations diagonal         
         do g1=1,2
            do g2=1,2
               sldlmx(g1,g2)=0.d0
               sldlmx(g1+3,g2)=0.d0
               sldrmx(g1,g2)=0.d0
               sldrmx(g1+3,g2)=0.d0
               squlmx(g1,g2)=0.d0
               squlmx(g1+3,g2)=0.d0
               squrmx(g1,g2)=0.d0
               squrmx(g1+3,g2)=0.d0
               sqdlmx(g1,g2)=0.d0
               sqdlmx(g1+3,g2)=0.d0
               sqdrmx(g1,g2)=0.d0
               sqdrmx(g1+3,g2)=0.d0
            enddo
         enddo
         sldlmx(1,1)=1.d0
         sldlmx(2,2)=1.d0
         sldrmx(4,1)=1.d0
         sldrmx(5,2)=1.d0

         squlmx(1,1)=1.d0
         squlmx(2,2)=1.d0
         squrmx(4,1)=1.d0
         squrmx(5,2)=1.d0

         sqdlmx(1,1)=1.d0
         sqdlmx(2,2)=1.d0
         sqdrmx(4,1)=1.d0
         sqdrmx(5,2)=1.d0
         
c... Now read in third generation mixings
c...Block STAUMIX
         call dssetfromslhac(StauMix_USf(1,1),sldlmx(3,3),2,'STAUMIX')
         call dssetfromslhac(StauMix_USf(1,2),sldrmx(3,3),2,'STAUMIX') 
         call dssetfromslhac(StauMix_USf(2,1),sldlmx(6,3),2,'STAUMIX')
         call dssetfromslhac(StauMix_USf(2,2),sldrmx(6,3),2,'STAUMIX')

c...Block STOPMIX
         call dssetfromslhac(StopMix_USf(1,1),squlmx(3,3),2,'STOPMIX')
         call dssetfromslhac(StopMix_USf(1,2),squrmx(3,3),2,'STOPMIX')
         call dssetfromslhac(StopMix_USf(2,1),squlmx(6,3),2,'STOPMIX')
         call dssetfromslhac(StopMix_USf(2,2),squrmx(6,3),2,'STOPMIX')

c...Block SBOTMIX
         call dssetfromslhac(SbotMix_USf(1,1),sqdlmx(3,3),2,'SBOTMIX')
         call dssetfromslhac(SbotMix_USf(1,2),sqdrmx(3,3),2,'SBOTMIX')
         call dssetfromslhac(SbotMix_USf(2,1),sqdlmx(6,3),2,'SBOTMIX')
         call dssetfromslhac(SbotMix_USf(2,2),sqdrmx(6,3),2,'SBOTMIX')


c...Block AE
c...Block AU
c...Block AD
      if (plev.gt.0.and.Ae_Q.ne.invalid.and.
     &   dsabsq(Ae_Q-mass(kz)).gt.2.d0) then
         if (prl.ge.1) then
         write(*,*) ('-',i=1,70)
         write(*,*)
     &   'DS Warning: SLHA scale for block AE is not MZ as expected,'
         write(*,*)
     &    'it is instead ',Ae_Q,'. Will use this anyway'
         write(*,*) ('-',i=1,70)
         endif
      endif

      if (plev.gt.0.and.Au_Q.ne.invalid.and.
     &   dsabsq(Au_Q-mass(kz)).gt.2.d0) then
         if (prl.ge.1) then
         write(*,*) ('-',i=1,70)
         write(*,*)
     &   'DS Warning: SLHA scale for block AU is not MZ as expected,'
         write(*,*)
     &    'it is instead ',Au_Q,'. Will use this anyway'
         write(*,*) ('-',i=1,70)
         endif
      endif

      if (plev.gt.0.and.Ad_Q.ne.invalid.and.
     &   dsabsq(Ad_Q-mass(kz)).gt.2.d0) then
         if (prl.ge.1) then
         write(*,*) ('-',i=1,70)
         write(*,*)
     &   'DS Warning: SLHA scale for block AD is not MZ as expected,'
         write(*,*)
     &    'it is instead ',Ad_Q,'. Will use this anyway'
         write(*,*) ('-',i=1,70)
         endif
      endif

c...Trilinear couplings
      do g1=1,3
         asofte(g1)=0.d0
         asoftu(g1)=0.d0
         asoftd(g1)=0.d0
         call dssetfromslha(Ae_Af(g1,g1),asofte(g1),0,'')
         call dssetfromslha(Au_Af(g1,g1),asoftu(g1),0,'')
         call dssetfromslha(Ad_Af(g1,g1),asoftd(g1),0,'')
      enddo

      else ! opt = 1, full 6x6

c...Block TE
c...Block TU
c...Block TD
c...Currently, we only use the diagonal elements
      do g1=1,3
         if (Te_Tf(g1,g1).ne.invalid) then
            asofte(g1)=Te_Tf(g1,g1)/yukawa(kl(g1))
         endif

         if (Tu_Tf(g1,g1).ne.invalid) then
            asoftu(g1)=Tu_Tf(g1,g1)/yukawa(kqu(g1))
         endif

         if (Td_Tf(g1,g1).ne.invalid) then
            asoftd(g1)=Td_Tf(g1,g1)/yukawa(kqd(g1))
         endif
      enddo

c...Block SNMIX - Sneutrino mixing matrices
c...JE CHECK: Is this OK for sneutrinos? s1,s2 = 1,...,6
      do s1=1,3
         do s2=1,3
            call dssetfromslhac(SnuMix_UASf(s1,s2),
     &         slulmx(s1,s2),0,'')
         enddo
      enddo

c...Block SLMIX - Charged slepton mixing matrices
      do s1=1,6
         do s2=1,3
            call dssetfromslhac(SelMix_UASf(s1,s2),
     &        sldlmx(s1,s2),0,'')
            call dssetfromslhac(SelMix_UASf(s1,s2+3),
     &        sldrmx(s1,s2),0,'')
         enddo
      enddo

c...Block USQMIX - Up squark mixing matrices
      do s1=1,6
         do s2=1,3
            call dssetfromslhac(USqMix_UASf(s1,s2),
     &         squlmx(s1,s2),0,'')
            call dssetfromslhac(USqMix_UASf(s1,s2+3),
     &         squrmx(s1,s2),0,'')
         enddo
      enddo

c...Block DSQMIX - Down squark mixing matrices
      do s1=1,6
         do s2=1,3
            call dssetfromslhac(DSqMix_UASf(s1,s2),
     &         sqdlmx(s1,s2),0,'')
            call dssetfromslhac(DSqMix_UASf(s1,s2+3),
     &         sqdrmx(s1,s2),0,'')
         enddo
      enddo

      endif ! end of opt=1


c...Now we need to make sure that things needed for further calculations
c...are added. In particular, we need to set vertices and (s)particle
c...widths. In the future, we should read these from the SLHA file instead
c...and only calculate here if they are not present in the file.

      call dsvertx
      call dshigwid
      call dsspwid


      return
      end

