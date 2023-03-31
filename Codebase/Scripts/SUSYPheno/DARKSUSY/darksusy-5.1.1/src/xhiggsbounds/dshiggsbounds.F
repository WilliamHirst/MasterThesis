c This is a modified version of the HBwithFH.F file supplied with
c HiggsBounds 2.8.5.

c This file assumes that some FeyhHiggs subroutines have been
c called previously in dsfeynhiggs. They are marked below with 
c 'done in dsfeynhiggs'.

* This is modified version of the demo.F file which is supplied with 
* FeynHiggs (vers 2.8.0, downloaded 16.05.11, www.feynhiggs.de)
* This file is part of the HiggsBounds distribution.

* The phrase 'currently' as used in this file means 'correct
* on the date 16.05.11'.
* See www.feynhiggs.de for the latest version of FeynHiggs
* and the corresponding documentation.


	subroutine dshiggsbounds(HBverdict)
	implicit none

	include 'dsmpconst.h' ! after removing Pi below

	integer HBverdict

	integer error

c used by FHHiggsCorr
	double precision MHiggs(4)
	double complex SAeff, UHiggs(3,3), ZHiggs(3,3)

c used by FHSelectUZ:
	integer uzint, uzext

c used by FHCouplings:
#include "FHCouplings.h"
        double complex couplings(ncouplings), couplingsms(ncouplingsms)
        double precision gammas(ngammas), gammasms(ngammasms)
        integer fast

c used by FHHiggsProd:
        double precision sqrts, prodxs(nprodxs)

c used by FHGetPara:
	integer nmfv
	double precision MASf(6,4), MCha(2), MNeu(4)
	double complex UASf(6,6,4)
	double complex UCha(2,2), VCha(2,2), ZNeu(4,4)
	double complex DeltaMB
	double precision MGl
	double precision MHtree(4), SAtree

c used by FHRetrieveSMPara:	
        double precision invAlfa, AlfasMZ, GF
	double precision ME, MU, MD, MM, MC, MS, ML, MB
	double precision MW, MZ
	double precision CKMlambda, CKMA, CKMrhobar, CKMetabar

c used by initialize_HiggsBounds
        integer nHiggsneut,nHiggsplus
        character(LEN=5) whichanalyses

c used by HiggsBounds_neutral_input_part
        double precision Mh(3),GammaTotal_hj(3)
        integer CP_value(3)
        double precision  CS_lep_hjZ_ratio(3),          
     &          CS_lep_bbhj_ratio(3),CS_lep_tautauhj_ratio(3),
     &          CS_lep_hjhi_ratio_nHbynH(3,3),               
     &          CS_gg_hj_ratio(3),CS_bb_hj_ratio(3),  
     &          CS_bg_hjb_ratio(3),                       
     &          CS_ud_hjWp_ratio(3),CS_cs_hjWp_ratio(3),
     &          CS_ud_hjWm_ratio(3),CS_cs_hjWm_ratio(3), 
     &          CS_gg_hjZ_ratio(3),
     &          CS_dd_hjZ_ratio(3),CS_uu_hjZ_ratio(3),
     &          CS_ss_hjZ_ratio(3),CS_cc_hjZ_ratio(3), 
     &          CS_bb_hjZ_ratio(3),                        
     &          CS_tev_vbf_ratio(3),CS_tev_tthj_ratio(3),
     &          CS_lhc7_vbf_ratio(3),CS_lhc7_tthj_ratio(3),
     &          BR_hjss(3),BR_hjcc(3),                         
     &          BR_hjbb(3),BR_hjmumu(3),BR_hjtautau(3),                     
     &          BR_hjWW(3),BR_hjZZ(3),BR_hjZga(3),                     
     &          BR_hjgaga(3),BR_hjgg(3),
     &          BR_hjinvisible(3),BR_hjhihi_nHbynH(3,3)

c used by HiggsBounds_charged_input
        double precision Mhplus(1),GammaTotal_Hpj(1),
     &          CS_lep_HpjHmj_ratio(1),                    
     &          BR_tWpb,BR_tHpjb(1),                     
     &          BR_Hpjcs(1),BR_Hpjcb(1),BR_Hpjtaunu(1) 

c used by run_HiggsBounds
        integer HBresult,chan,ncombined
        double precision obsratio

c misc:
        integer i,j,as,t
        double precision norm,CW2
        double precision
     &          g2hjbb(3),g2hjWW(3),g2hjZZ(3),                
     &          g2hjgg(3),g2hjhiZ_nHbynH(3,3)
        double precision g2hjbb_s(3),g2hjbb_p(3)
        double precision g2hjtautau_s(3),g2hjtautau_p(3)
        integer sneutrino_lspcandidate_number
        logical invisible_lsp
        double precision lspcandidate_mass 

*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *

c Set number of neutral and charged Higgs bosons in the MSSM:
        nHiggsneut=3
        nHiggsplus=1

c The string 'whichanalyses' determines which subset of experimental 
c results are used.
c In this example, we've used the option 'LandH',
c which instructs HiggsBounds to use tables of results
c from LEP, Tevatron and LHC (i.e. the full set of
c results supplied with HiggsBounds).
c        whichanalyses='LandH'
        whichanalyses='onlyL'

c The subroutine initialize_HiggsBounds reads in all necessary
c tables etc.
c It must be called before any of the other HiggsBounds subroutines.

c done in dsinit
c        call initialize_HiggsBounds(nHiggsneut,nHiggsplus,whichanalyses)

c If you would like to perform scans over variables, the subroutine
c initialize_HiggsBounds (and finish_HiggsBounds) should be called
c outside the do-loops in order to save time.

*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *

c calls to FH subroutines:

c the subroutines setFlags,setPara,setSLHA
c are also contained in this file

c done in dsfeynhiggs
c	call setFlags

* either use setPara to set the parameters directly
* or use setSLHA to read them from an SLHA file
c done in dsfeynhiggs
c	call setPara
cc	call setSLHA("tmp.slha")

	call FHGetPara(error, nmfv, MASf, UASf,
     &    MCha, UCha, VCha, MNeu, ZNeu, DeltaMB, MGl,
     &    MHtree, SAtree)
	if( error .ne. 0 ) stop

c done in dsfeynhiggs
c	call FHHiggsCorr(error, MHiggs, SAeff, UHiggs, ZHiggs)
c	if( error .ne. 0 ) stop

c done in dsfeynhiggs
c NOTE: we are setting uzint=uzext
c        uzint=2
c        uzext=2
c        call FHSelectUZ(error, uzint, uzext)c
	if( error .ne. 0 ) stop

c done in dsfeynhiggs
c        fast=1
c	call FHCouplings(error,
c     &    couplings, couplingsms, gammas, gammasms, fast)
c	if( error .ne. 0 ) stop

c We would like FH to calculate LHC cross sections
        sqrts=14.0D0
        call FHHiggsProd(error, sqrts, prodxs)
	if( error .ne. 0 ) stop 

        call FHRetrieveSMPara(error,
     &    invAlfa, AlfasMZ, GF,
     &    ME, MU, MD, MM, MC, MS, ML, MB,
     &    MW, MZ,
     &    CKMlambda, CKMA, CKMrhobar, CKMetabar)

*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *      

c Set variables needed by HiggsBounds (using results from FeynHiggs).
c See HiggsBounds documentation for definition of variables used
c as arguments to HiggsBounds_neutral_input_part and run_HiggsBounds
c and FeynHiggs documentation for all other variables.

c Note: It is slightly more accurate to use the subroutine HiggsBounds_neutral_input_part
c rather than the subroutine HiggsBounds_neutral_input_effC because the SM branching ratios
c used internally in HiggsBounds (from HDecay) are not identical to the SM branching
c ratios used in FeynHiggs


        do i=1,3

         Mh(i)=MHiggs(i) 
         GammaTotal_hj(i) = GammaTot(i)

         BR_hjss(i)       = BR(H0FF(i,4,2,2))
         BR_hjcc(i)       = BR(H0FF(i,3,2,2))
         BR_hjbb(i)       = BR(H0FF(i,4,3,3))
         BR_hjmumu(i)     = BR(H0FF(i,2,2,2))
         BR_hjtautau(i)   = BR(H0FF(i,2,3,3)) 

         BR_hjWW(i)     = BR(H0VV(i,4))
         BR_hjgaga(i)   = BR(H0VV(i,1))
         BR_hjZga(i)    = BR(H0VV(i,2))
         BR_hjZZ(i)     = BR(H0VV(i,3))
         BR_hjgg(i)     = BR(H0VV(i,5))

         if(GammaSM(H0FF(i,4,3,3)).le.0.0D0)then
          g2hjbb(i)=0.0D0
         else
          g2hjbb(i)=Gamma(H0FF(i,4,3,3)) 
     &      /GammaSM(H0FF(i,4,3,3))
         endif

c Note that this is currently equivalent to
c         g2hjbb(i)= bbh(i)/bbhSM(i)
c         g2hjbb(i)= btagbh(i)/btagbhSM(i)
c as long as MH>80 GeV

         CS_bg_hjb_ratio(i) = g2hjbb(i)
         CS_bb_hj_ratio(i)  = g2hjbb(i)
   
         g2hjbb_s(i)=(abs(RCoupling(H0FF(i,4,3,3))
     &                   /RCouplingSM(H0FF(i,4,3,3))+ 
     &                    LCoupling(H0FF(i,4,3,3))
     &                   /LCouplingSM(H0FF(i,4,3,3)))/2.0D0)**2.0D0
         g2hjbb_p(i)=(abs(RCoupling(H0FF(i,4,3,3))
     &                   /RCouplingSM(H0FF(i,4,3,3))- 
     &                    LCoupling(H0FF(i,4,3,3))
     &                   /LCouplingSM(H0FF(i,4,3,3)))/2.0D0)**2.0D0

         g2hjtautau_s(i)=(abs(RCoupling(H0FF(i,2,3,3))
     &                   /RCouplingSM(H0FF(i,2,3,3))+ 
     &                    LCoupling(H0FF(i,2,3,3))
     &                   /LCouplingSM(H0FF(i,2,3,3)))/2.0D0)**2.0D0
         g2hjtautau_p(i)=(abs(RCoupling(H0FF(i,2,3,3))
     &                   /RCouplingSM(H0FF(i,2,3,3))- 
     &                    LCoupling(H0FF(i,2,3,3))
     &                   /LCouplingSM(H0FF(i,2,3,3)))/2.0D0)**2.0D0

         if(       g2hjbb_p(i).lt.1.0D-10)then
            CP_value(i) = 1
         elseif(   g2hjbb_s(i).lt.1.0D-10)then
            CP_value(i) = -1
         else
            CP_value(i) = 0
         endif

         CS_lep_bbhj_ratio(i)     = g2hjbb_s(i)+g2hjbb_p(i)
         CS_lep_tautauhj_ratio(i) = g2hjtautau_s(i)+g2hjtautau_p(i)

         g2hjWW(i)= dble(  Coupling(H0VV(i,4)) 
     &                   / CouplingSM(H0VV(i,4)) )**2.0D0
     &            + dimag( Coupling(H0VV(i,4)) 
     &                   / CouplingSM(H0VV(i,4)) )**2.0D0
c Note that this is currently equivalent to
c         g2hjWW(i)= WhTev(i)/WhTevSM(i
c	  g2hjWW(i)= qqhTev(i)/qqhTevSM(i)
c as long as MH>80 GeV and uzint=uzext

         g2hjZZ(i)= dble(  Coupling(H0VV(i,3)) 
     &                   / CouplingSM(H0VV(i,3)) )**2.0D0
     &            + dimag( Coupling(H0VV(i,3)) 
     &                   / CouplingSM(H0VV(i,3)) )**2.0D0
c Note that this is currently equivalent to
c         g2hjZZ(i)= ZhTev(i)/ZhTevSM(i)
c as long as MH>80 GeV and uzint=uzext
c It is also equivalent to g2hjWW(i)

         CS_lep_hjZ_ratio(i)        = g2hjZZ(i)

         CS_gg_hjZ_ratio(i)     = 0.0D0
         CS_dd_hjZ_ratio(i)     = g2hjZZ(i)
         CS_uu_hjZ_ratio(i)     = g2hjZZ(i)
         CS_ss_hjZ_ratio(i)     = g2hjZZ(i)
         CS_cc_hjZ_ratio(i)     = g2hjZZ(i)
         CS_bb_hjZ_ratio(i)     = g2hjZZ(i)

         CS_ud_hjWp_ratio(i)    = g2hjZZ(i)
         CS_cs_hjWp_ratio(i)    = g2hjZZ(i)
         CS_ud_hjWm_ratio(i)    = g2hjZZ(i)
         CS_cs_hjWm_ratio(i)    = g2hjZZ(i)

         CS_tev_vbf_ratio(i)     = g2hjZZ(i)
         CS_lhc7_vbf_ratio(i)    = g2hjZZ(i)


         if(tthSM(i).gt.0.0D0)then
           CS_tev_tthj_ratio(i)    = tth(i)/tthSM(i)  
         else
           CS_tev_tthj_ratio(i)    = 0.0D0
         endif

         CS_lhc7_tthj_ratio(i) = CS_tev_tthj_ratio(i)

c tevatron gluon fusion XS is not calculated in FH is MH<90 geV
         if(Mh(i).gt.90.0001D0)then      
           if(gghSM(i).gt.0.0D0)then
            CS_gg_hj_ratio(i) = ggh(i)/gghSM(i)
           else
            CS_gg_hj_ratio(i) = 0.0D0
           endif 
         else
           if(GammaSM(H0VV(i,5)).le.0.0D0)then
            CS_gg_hj_ratio(i)=0.0D0
           else
            CS_gg_hj_ratio(i)= Gamma(H0VV(i,5))/GammaSM(H0VV(i,5))
          endif
         endif

        enddo
  
        norm=GF*sqrt(2.0D0)*MZ**2.0D0

        do j=1,3
          do i=1,3   
            g2hjhiZ_nHbynH(j,i)= (
     &          dble(  Coupling(H0HV(j,i)) )**2.0D0
     &        + dimag( Coupling(H0HV(j,i)) )**2.0D0
     &          )
     &       /norm
            
            CS_lep_hjhi_ratio_nHbynH(j,i) = g2hjhiZ_nHbynH(j,i)

            BR_hjhihi_nHbynH(j,i)=BR(H0HH(j,i,i))
          enddo
        enddo

c higgs->neutralino1 neutralino1 contributes the invisible Higgs decay width
c when neutralino1 or sneutrino is the LSP

        do i=1,3  
         sneutrino_lspcandidate_number=0
         invisible_lsp=.True.

c first determine whether lightest sneutrino is lighter than the lightest neutralino
c
c sneutrino_lspcandidate_number=0 indicates that lightest neutralino is 
c lighter than all the sneutrinos
         lspcandidate_mass=MNeu(1)
         do as=1,3
            if( MASf(as,1) .lt. lspcandidate_mass )then
               lspcandidate_mass=MASf(as,1)
               sneutrino_lspcandidate_number=as
            endif
         enddo

         if(     MCha(1) .lt. lspcandidate_mass )then 
           invisible_lsp=.False.
         elseif( MGl     .lt. lspcandidate_mass )then
           invisible_lsp=.False.
         else
          do as=1,6
           do t=2,4
            if( MASf(as,t) .lt. lspcandidate_mass )then
               invisible_lsp=.False.
            endif
           enddo
          enddo      
         endif

         if(invisible_lsp)then
          if(sneutrino_lspcandidate_number.eq.0)then
            BR_hjinvisible(i) = BR(H0NeuNeu(i,1,1))
          else
            BR_hjinvisible(i) = BR(H0SfSf(i,1,1,1,as))
          endif
         else
          BR_hjinvisible(i) = 0.0D0
         endif

        enddo

*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *

        Mhplus(1)              = MHiggs(4) 
        GammaTotal_Hpj(1)      = GammaTot(4)  
        CS_lep_HpjHmj_ratio(1) = 1.0D0
        BR_tWpb                = BR( tBF(1) ) 
        BR_tHpjb(1)            = BR( tBF(2) )
        BR_Hpjcs(1)            = BR( HpFF(2,2,2) ) 
        BR_Hpjcb(1)            = BR( HpFF(2,2,3) )
        BR_Hpjtaunu(1)         = BR( HpFF(1,3,3) )

*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
c calls to HiggsBounds_neutral_input_part,HiggsBounds_charged_input,
c which give input to HiggsBounds

c        write (*,*) 'DEBUG> calling HiggsBounds_neutral_input_part'
  
        call HiggsBounds_neutral_input_part(Mh,GammaTotal_hj,CP_value, 
     &          CS_lep_hjZ_ratio,                            
     &          CS_lep_bbhj_ratio,CS_lep_tautauhj_ratio,     
     &          CS_lep_hjhi_ratio_nHbynH,                    
     &          CS_gg_hj_ratio,CS_bb_hj_ratio,       
     &          CS_bg_hjb_ratio,                         
     &          CS_ud_hjWp_ratio,CS_cs_hjWp_ratio,    
     &          CS_ud_hjWm_ratio,CS_cs_hjWm_ratio,    
     &          CS_gg_hjZ_ratio,     
     &          CS_dd_hjZ_ratio,CS_uu_hjZ_ratio,     
     &          CS_ss_hjZ_ratio,CS_cc_hjZ_ratio,     
     &          CS_bb_hjZ_ratio,                         
     &          CS_tev_vbf_ratio,CS_tev_tthj_ratio,    
     &          CS_lhc7_vbf_ratio,CS_lhc7_tthj_ratio,    
     &          BR_hjss,BR_hjcc,                             
     &          BR_hjbb,BR_hjmumu,BR_hjtautau,               
     &          BR_hjWW,BR_hjZZ,BR_hjZga, BR_hjgaga,BR_hjgg,  
     &          BR_hjinvisible,BR_hjhihi_nHbynH              )
     
c        write (*,*) 'DEBUG> calling HiggsBounds_charged_input'
        call HiggsBounds_charged_input(Mhplus,GammaTotal_Hpj, 
     &          CS_lep_HpjHmj_ratio,                        
     &          BR_tWpb,BR_tHpjb,                           
     &          BR_Hpjcs,BR_Hpjcb,BR_Hpjtaunu)


*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
c call to run_HiggsBounds
        call run_HiggsBounds( HBresult,chan, obsratio, ncombined )

c        write (*,*) 'DEBUG>  '      
c        write (*,*) 'DEBUG> *************    HiggsBounds Results  **************'
c        write (*,*) 'DEBUG>  '    
c        write (*,*) 'DEBUG> Is this parameter point excluded by LEP, Tevatron'
c        write (*,*) 'DEBUG> or LHC data?'         
c        print*, HBresult, ',  where'
c        write (*,*) 'DEBUG>                0 = yes, it is excluded'
c        write (*,*) 'DEBUG>                1 = no, it has not been excluded'
c        write (*,*) 'DEBUG>               -1 = invalid parameter set'    
c        write (*,*) 'DEBUG>  '    
c        write (*,*) 'DEBUG> The process with the highest statistical sensitivity'
c        write (*,*) 'DEBUG> is'
c        print*, chan,'(see Key.dat)'
c        write (*,*) 'DEBUG> This process has a theoretical rate vs. limit of'
c        print*, obsratio
c        write (*,*) 'DEBUG>  '    
c        write (*,*) 'DEBUG> The number of Higgs which have contributed to the'
c        write (*,*) 'DEBUG> theoretical rate of this process was'
c        print*, ncombined
c        write (*,*) 'DEBUG>  '    
c        write (*,*) 'DEBUG> See HiggsBounds documentation for more information.'
c        write (*,*) 'DEBUG> ****************************************************'
c        write (*,*) 'DEBUG>  '    


*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
c deallocates arrays used by HiggsBounds:

c should be called at end of darksusy, but it is not (memory leakage?)
c        call finish_HiggsBounds

	HBverdict = HBresult

	end
