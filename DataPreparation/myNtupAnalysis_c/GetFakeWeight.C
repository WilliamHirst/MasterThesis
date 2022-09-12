#include <iostream>
#include <math.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TLorentzVector.h>
#include "MatrixMethod.h"
#include <map>
#include <TH1.h>
#include <TFile.h>
#include <TString.h>
#include <TH2.h>
#include <math.h>
#include <numeric>
#include <unordered_map>
#include "GetFakeWeight.h"
 

using namespace std;



//double zmass = 91.1876;



void GetFakeWeight::initializeHist(TString LTDef, TString eff_file, TString frac_file, bool doThreeLep, TString var, bool is1516, bool is17, bool is18, std::vector<TString> uncert_keys, std::vector<std::string> bdt_vec){

  TFile *MMfile_real       = new TFile("./MMinputfiles/FNPntuple_FCL_JUNE09_fromDATA_BDTscores_R21.root");
  TFile *MMfile       = new TFile(eff_file);
  TFile *MMfile_frac  = new TFile(frac_file);
  TFile *MMfile_frac_extra  = new TFile("./MMinputfiles/MMinput_frac_onlyFAKE2L21_FRAC_JUNE10.root");
  TString appendtostring = "R21";
  TString isokey;

  Bool_t use2Dreal = false;
  
  if(LTDef.Contains("2L2J"))isokey = LTDef;
  else isokey = "FCL";

  cout<<"LTDef = "<<LTDef.Data()<<endl;
  cout<<"isokey = "<<isokey.Data()<<endl;

  setFracVar(var);
  
  if(doThreeLep){
    for(int i = 4; i<=9; i++){
      if(!(i==4))continue;
      cout<<"Adding 3L: "<<i<<endl;
      h_el_frac_hf[Form("EEE_0%i_%s",i,isokey.Data())] = (TH1F*)MMfile_frac->Get(Form("h_fract_h_lep_"+var+"_TCF_nL_EEE_E_trueFAKE_REAL2L0%i_"+isokey+"_HF%s",i,appendtostring.Data()));
      h_el_frac_hf[Form("EEM_0%i_%s",i,isokey.Data())] = (TH1F*)MMfile_frac->Get(Form("h_fract_h_lep_"+var+"_TCF_nL_EEM_E_trueFAKE_REAL2L0%i_"+isokey+"_HF%s",i,appendtostring.Data()));
      h_el_frac_hf[Form("MME_0%i_%s",i,isokey.Data())] = (TH1F*)MMfile_frac->Get(Form("h_fract_h_lep_"+var+"_TCF_nL_MME_E_trueFAKE_REAL2L0%i_"+isokey+"_HF%s",i,appendtostring.Data()));

      h_el_frac_lf[Form("EEE_0%i_%s",i,isokey.Data())] = (TH1F*)MMfile_frac->Get(Form("h_fract_h_lep_"+var+"_TCF_nL_EEE_E_trueFAKE_REAL2L0%i_"+isokey+"_LF%s",i,appendtostring.Data()));
      h_el_frac_lf[Form("EEM_0%i_%s",i,isokey.Data())] = (TH1F*)MMfile_frac->Get(Form("h_fract_h_lep_"+var+"_TCF_nL_EEM_E_trueFAKE_REAL2L0%i_"+isokey+"_LF%s",i,appendtostring.Data()));
      h_el_frac_lf[Form("MME_0%i_%s",i,isokey.Data())] = (TH1F*)MMfile_frac->Get(Form("h_fract_h_lep_"+var+"_TCF_nL_MME_E_trueFAKE_REAL2L0%i_"+isokey+"_LF%s",i,appendtostring.Data()));

      h_el_frac_co[Form("EEE_0%i_%s",i,isokey.Data())] = (TH1F*)MMfile_frac->Get(Form("h_fract_h_lep_"+var+"_TCF_nL_EEE_E_trueFAKE_REAL2L0%i_"+isokey+"_CO%s",i,appendtostring.Data()));
      h_el_frac_co[Form("EEM_0%i_%s",i,isokey.Data())] = (TH1F*)MMfile_frac->Get(Form("h_fract_h_lep_"+var+"_TCF_nL_EEM_E_trueFAKE_REAL2L0%i_"+isokey+"_CO%s",i,appendtostring.Data()));
      h_el_frac_co[Form("MME_0%i_%s",i,isokey.Data())] = (TH1F*)MMfile_frac->Get(Form("h_fract_h_lep_"+var+"_TCF_nL_MME_E_trueFAKE_REAL2L0%i_"+isokey+"_CO%s",i,appendtostring.Data()));

    }
  }else{
    for(int i = 2; i<=13; i++){
      if(isokey.Contains("2L2J") && i > 12)continue;
      cout<<"Adding 2L: "<<(Form("h_fract_h_lep_"+var+"_TCF_nL_EEOS_E_trueFAKE_REAL2L%02i_"+isokey+"_BHadronDecay%s",i,appendtostring.Data()))<<endl;
      h_el_frac_hf[Form("EMOS_%02i_%s",i,isokey.Data())] = (TH1F*)MMfile_frac->Get(Form("h_fract_h_lep_"+var+"_TCF_nL_EMOS_E_trueFAKE_REAL2L%02i_"+isokey+"_BHadronDecay%s",i,appendtostring.Data()));
      h_el_frac_hf[Form("EMSS_%02i_%s",i,isokey.Data())] = (TH1F*)MMfile_frac->Get(Form("h_fract_h_lep_"+var+"_TCF_nL_EMSS_E_trueFAKE_REAL2L%02i_"+isokey+"_BHadronDecay%s",i,appendtostring.Data()));
      h_el_frac_hf[Form("EEOS_%02i_%s",i,isokey.Data())] = (TH1F*)MMfile_frac->Get(Form("h_fract_h_lep_"+var+"_TCF_nL_EEOS_E_trueFAKE_REAL2L%02i_"+isokey+"_BHadronDecay%s",i,appendtostring.Data()));
      h_el_frac_hf[Form("EESS_%02i_%s",i,isokey.Data())] = (TH1F*)MMfile_frac->Get(Form("h_fract_h_lep_"+var+"_TCF_nL_EESS_E_trueFAKE_REAL2L%02i_"+isokey+"_BHadronDecay%s",i,appendtostring.Data()));

      h_el_frac_cf[Form("EMOS_%02i_%s",i,isokey.Data())] = (TH1F*)MMfile_frac->Get(Form("h_fract_h_lep_"+var+"_TCF_nL_EMOS_E_trueFAKE_REAL2L%02i_"+isokey+"_CHadronDecay%s",i,appendtostring.Data()));
      h_el_frac_cf[Form("EMSS_%02i_%s",i,isokey.Data())] = (TH1F*)MMfile_frac->Get(Form("h_fract_h_lep_"+var+"_TCF_nL_EMSS_E_trueFAKE_REAL2L%02i_"+isokey+"_CHadronDecay%s",i,appendtostring.Data()));
      h_el_frac_cf[Form("EEOS_%02i_%s",i,isokey.Data())] = (TH1F*)MMfile_frac->Get(Form("h_fract_h_lep_"+var+"_TCF_nL_EEOS_E_trueFAKE_REAL2L%02i_"+isokey+"_CHadronDecay%s",i,appendtostring.Data()));
      h_el_frac_cf[Form("EESS_%02i_%s",i,isokey.Data())] = (TH1F*)MMfile_frac->Get(Form("h_fract_h_lep_"+var+"_TCF_nL_EESS_E_trueFAKE_REAL2L%02i_"+isokey+"_CHadronDecay%s",i,appendtostring.Data()));

      h_el_frac_lf[Form("EMOS_%02i_%s",i,isokey.Data())] = (TH1F*)MMfile_frac->Get(Form("h_fract_h_lep_"+var+"_TCF_nL_EMOS_E_trueFAKE_REAL2L%02i_"+isokey+"_LightFlavorDecay%s",i,appendtostring.Data()));
      h_el_frac_lf[Form("EMSS_%02i_%s",i,isokey.Data())] = (TH1F*)MMfile_frac->Get(Form("h_fract_h_lep_"+var+"_TCF_nL_EMSS_E_trueFAKE_REAL2L%02i_"+isokey+"_LightFlavorDecay%s",i,appendtostring.Data()));
      h_el_frac_lf[Form("EEOS_%02i_%s",i,isokey.Data())] = (TH1F*)MMfile_frac->Get(Form("h_fract_h_lep_"+var+"_TCF_nL_EEOS_E_trueFAKE_REAL2L%02i_"+isokey+"_LightFlavorDecay%s",i,appendtostring.Data()));
      h_el_frac_lf[Form("EESS_%02i_%s",i,isokey.Data())] = (TH1F*)MMfile_frac->Get(Form("h_fract_h_lep_"+var+"_TCF_nL_EESS_E_trueFAKE_REAL2L%02i_"+isokey+"_LightFlavorDecay%s",i,appendtostring.Data()));

      h_el_frac_co[Form("EMOS_%02i_%s",i,isokey.Data())] = (TH1F*)MMfile_frac->Get(Form("h_fract_h_lep_"+var+"_TCF_nL_EMOS_E_trueFAKE_REAL2L%02i_"+isokey+"_PromptPhotonConversion%s",i,appendtostring.Data()));
      h_el_frac_co[Form("EMSS_%02i_%s",i,isokey.Data())] = (TH1F*)MMfile_frac->Get(Form("h_fract_h_lep_"+var+"_TCF_nL_EMSS_E_trueFAKE_REAL2L%02i_"+isokey+"_PromptPhotonConversion%s",i,appendtostring.Data()));
      h_el_frac_co[Form("EEOS_%02i_%s",i,isokey.Data())] = (TH1F*)MMfile_frac->Get(Form("h_fract_h_lep_"+var+"_TCF_nL_EEOS_E_trueFAKE_REAL2L%02i_"+isokey+"_PromptPhotonConversion%s",i,appendtostring.Data()));
      h_el_frac_co[Form("EESS_%02i_%s",i,isokey.Data())] = (TH1F*)MMfile_frac->Get(Form("h_fract_h_lep_"+var+"_TCF_nL_EESS_E_trueFAKE_REAL2L%02i_"+isokey+"_PromptPhotonConversion%s",i,appendtostring.Data()));

      h_el_frac_em[Form("EMOS_%02i_%s",i,isokey.Data())] = (TH1F*)MMfile_frac->Get(Form("h_fract_h_lep_"+var+"_TCF_nL_EMOS_E_trueFAKE_REAL2L%02i_"+isokey+"_ElectronFromMuon%s",i,appendtostring.Data()));
      h_el_frac_em[Form("EMSS_%02i_%s",i,isokey.Data())] = (TH1F*)MMfile_frac->Get(Form("h_fract_h_lep_"+var+"_TCF_nL_EMSS_E_trueFAKE_REAL2L%02i_"+isokey+"_ElectronFromMuon%s",i,appendtostring.Data()));
      h_el_frac_em[Form("EEOS_%02i_%s",i,isokey.Data())] = (TH1F*)MMfile_frac->Get(Form("h_fract_h_lep_"+var+"_TCF_nL_EEOS_E_trueFAKE_REAL2L%02i_"+isokey+"_ElectronFromMuon%s",i,appendtostring.Data()));
      h_el_frac_em[Form("EESS_%02i_%s",i,isokey.Data())] = (TH1F*)MMfile_frac->Get(Form("h_fract_h_lep_"+var+"_TCF_nL_EESS_E_trueFAKE_REAL2L%02i_"+isokey+"_ElectronFromMuon%s",i,appendtostring.Data()));

    }

    /**
    h_el_frac_hf[Form("ALL_%02i_%s",14,isokey.Data())] = (TH1F*)MMfile_frac_extra->Get(Form("h_fract_h_lep_"+var+"_TCF_nL_ALL_E_FAKE2L21_trueFAKE_"+isokey+"_BHadronDecay%s",appendtostring.Data()));
    h_el_frac_cf[Form("ALL_%02i_%s",14,isokey.Data())] = (TH1F*)MMfile_frac_extra->Get(Form("h_fract_h_lep_"+var+"_TCF_nL_ALL_E_FAKE2L21_trueFAKE_"+isokey+"_CHadronDecay%s",appendtostring.Data()));
    h_el_frac_lf[Form("ALL_%02i_%s",14,isokey.Data())] = (TH1F*)MMfile_frac_extra->Get(Form("h_fract_h_lep_"+var+"_TCF_nL_ALL_E_FAKE2L21_trueFAKE_"+isokey+"_LightFlavorDecay%s",appendtostring.Data()));
    h_el_frac_co[Form("ALL_%02i_%s",14,isokey.Data())] = (TH1F*)MMfile_frac_extra->Get(Form("h_fract_h_lep_"+var+"_TCF_nL_ALL_E_FAKE2L21_trueFAKE_"+isokey+"_PromptPhotonConversion%s",appendtostring.Data()));
    h_el_frac_em[Form("ALL_%02i_%s",14,isokey.Data())] = (TH1F*)MMfile_frac_extra->Get(Form("h_fract_h_lep_"+var+"_TCF_nL_ALL_E_FAKE2L21_trueFAKE_"+isokey+"_ElectronFromMuon%s",appendtostring.Data()));
    */
    
  }
  for(map< TString, TH1F* >::iterator map_it = h_el_frac_hf.begin(); map_it != h_el_frac_hf.end(); ++map_it) {
    cout<<"Fixing hf "<<map_it->first<<endl;
    h_el_frac_hf[map_it->first]->SetDirectory(0);
    h_el_frac_cf[map_it->first]->SetDirectory(0);
    cout<<"Fixing lf "<<map_it->first<<endl;
    h_el_frac_lf[map_it->first]->SetDirectory(0);
    cout<<"Fixing co "<<map_it->first<<endl;
    h_el_frac_co[map_it->first]->SetDirectory(0);
    cout<<"Fixing em "<<map_it->first<<endl;
    h_el_frac_em[map_it->first]->SetDirectory(0);
    //h_el_frac_cf[map_it->first]->SetDirectory(0);
  }

  
  if(LTDef.Contains("2L2J")){
    //NEW 2L0J Rates
    auto objPtr = MMfile->Get<TH1F>("eff_pt_el_2015_FAKE2L21_data15-16");
    if (objPtr)h_el_fake_hf["2015_2016_"+LTDef] = objPtr;
    objPtr = MMfile->Get<TH1F>("eff_pt_el_2017_FAKE2L21_data17");
    if (objPtr)h_el_fake_hf["2017_"+LTDef]      = objPtr;
    objPtr = MMfile->Get<TH1F>("eff_pt_el_2018_FAKE2L21_data18");
    if (objPtr)h_el_fake_hf["2018_"+LTDef]      = objPtr;


    objPtr = MMfile->Get<TH1F>("eff_pt_el_2015_FAKE2L21_xsecup_data15-16");
    if (objPtr)h_el_fake_hf_xsecup["2015_2016_"+LTDef] = objPtr;
    objPtr = MMfile->Get<TH1F>("eff_pt_el_2017_FAKE2L21_xsecup_data17");
    if (objPtr)h_el_fake_hf_xsecup["2017_"+LTDef]      = objPtr; 
    objPtr = MMfile->Get<TH1F>("eff_pt_el_2018_FAKE2L21_xsecup_data18");
    if (objPtr)h_el_fake_hf_xsecup["2018_"+LTDef]      = objPtr;

    objPtr = MMfile->Get<TH1F>("eff_pt_el_2015_FAKE2L21_xsecdw_data15-16");
    if (objPtr)h_el_fake_hf_xsecdw["2015_2016_"+LTDef] = objPtr;
    objPtr = MMfile->Get<TH1F>("eff_pt_el_2017_FAKE2L21_xsecdw_data17");
    if (objPtr)h_el_fake_hf_xsecdw["2017_"+LTDef]      = objPtr; 
    objPtr = MMfile->Get<TH1F>("eff_pt_el_2018_FAKE2L21_xsecdw_data18");
    if (objPtr)h_el_fake_hf_xsecdw["2018_"+LTDef]      = objPtr;


    objPtr = MMfile->Get<TH1F>("eff_pt_el_2015_FAKE2L23_data15-16");
    if (objPtr)h_el_fake_co["2015_2016_"+LTDef] = objPtr;
    objPtr = MMfile->Get<TH1F>("eff_pt_el_2017_FAKE2L23_data17");
    if (objPtr)h_el_fake_co["2017_"+LTDef]      = objPtr;
    objPtr = MMfile->Get<TH1F>("eff_pt_el_2018_FAKE2L23_data18");
    if (objPtr)h_el_fake_co["2018_"+LTDef]      = objPtr;

    objPtr = MMfile->Get<TH1F>("eff_pt_el_2015_FAKE2L05_mc16a");
    if (objPtr)h_el_fake_lf["2015_2016_"+LTDef] = objPtr;
    objPtr = MMfile->Get<TH1F>("eff_pt_el_2017_FAKE2L05_mc16d");
    if (objPtr)h_el_fake_lf["2017_"+LTDef]      = objPtr;
    objPtr = MMfile->Get<TH1F>("eff_pt_el_2018_FAKE2L05_mc16e");
    if (objPtr)h_el_fake_lf["2018_"+LTDef]      = objPtr;

    //if(!LTDef.Contains("tight_iso"))
    auto objPtr2D = MMfile->Get<TH2F>("eff_pt_eta_el_2015_REAL2L04_mc16a");
    printf("Adding %s\n",("2015_2016_"+LTDef).Data());
    if (objPtr2D)h_el_real["2015_2016_"+LTDef] = objPtr2D;
    else cout<<"Could not find "<<"eff_pt_eta_el_2015_REAL2L04_mc16a"<<endl;
    objPtr2D = MMfile->Get<TH2F>("eff_pt_eta_el_2017_REAL2L04_mc16d");
    if (objPtr2D) h_el_real["2017_"+LTDef] = objPtr2D;
    objPtr2D = MMfile->Get<TH2F>("eff_pt_eta_el_2018_REAL2L04_mc16e");
    if (objPtr2D) h_el_real["2018_"+LTDef] = objPtr2D;
    //h_el_real["2017_"+LTDef]      = (TH2F*)MMfile->Get("eff_el_2017_REAL2L04");
    //h_el_real["2018_"+LTDef]      = (TH2F*)MMfile->Get("eff_el_2018_REAL2L04");

    // add systematic histograms
    for(std::vector< TString >::iterator it = uncert_keys.begin(); it != uncert_keys.end(); ++it){
      if((*it).Contains("MUON"))continue;
      printf("Adding %s\n",("2015_2016_"+LTDef+"_"+(*it)).Data());
      objPtr = MMfile->Get<TH1F>("frac_pt_el_2015_FAKE2L05_"+(*it)+"_mc16a");
      if (objPtr)h_el_fake_lf["2015_2016_"+LTDef+"_"+(*it)] = objPtr;
      objPtr = MMfile->Get<TH1F>("frac_pt_el_2017_FAKE2L05_"+(*it)+"_mc16d");
      if (objPtr)h_el_fake_lf["2017_"+LTDef+"_"+(*it)]      = objPtr;
      objPtr = MMfile->Get<TH1F>("frac_pt_el_2018_FAKE2L05_"+(*it)+"_mc16e");
      if (objPtr)h_el_fake_lf["2018_"+LTDef+"_"+(*it)]      = objPtr;
      
      auto objPtr2D = MMfile->Get<TH2F>("frac_pt_eta_el_2015_REAL2L04_"+(*it)+"_mc16a");
      if (objPtr2D)h_el_real["2015_2016_"+LTDef+"_"+(*it)] = objPtr2D;
      else cout<<"Could not find "<<"2015_2016_"+LTDef+"_"+(*it)<<endl;
      objPtr2D = MMfile->Get<TH2F>("frac_pt_eta_el_2017_REAL2L04_"+(*it)+"_mc16d");
      if (objPtr2D) h_el_real["2017_"+LTDef+"_"+(*it)] = objPtr2D;
      objPtr2D = MMfile->Get<TH2F>("frac_pt_eta_el_2018_REAL2L04_"+(*it)+"_mc16e");
      if (objPtr2D) h_el_real["2018_"+LTDef+"_"+(*it)] = objPtr2D; 
    }
    
    //if(!LTDef.EqualTo("tight_iso")){
    objPtr = MMfile->Get<TH1F>("eff_pt_mu_2015_FAKE2L21_data15-16");
    if (objPtr)h_mu_fake_hf["2015_2016_"+LTDef] = objPtr;
    objPtr = MMfile->Get<TH1F>("eff_pt_mu_2017_FAKE2L21_data17");
    if (objPtr)h_mu_fake_hf["2017_"+LTDef]      = objPtr; 
    objPtr = MMfile->Get<TH1F>("eff_pt_mu_2018_FAKE2L21_data18");
    if (objPtr)h_mu_fake_hf["2018_"+LTDef]      = objPtr;

    objPtr = MMfile->Get<TH1F>("eff_pt_mu_2015_FAKE2L21_xsecup_data15-16");
    if (objPtr)h_mu_fake_hf_xsecup["2015_2016_"+LTDef] = objPtr;
    objPtr = MMfile->Get<TH1F>("eff_pt_mu_2017_FAKE2L21_xsecup_data17");
    if (objPtr)h_mu_fake_hf_xsecup["2017_"+LTDef]      = objPtr; 
    objPtr = MMfile->Get<TH1F>("eff_pt_mu_2018_FAKE2L21_xsecup_data18");
    if (objPtr)h_mu_fake_hf_xsecup["2018_"+LTDef]      = objPtr;

    objPtr = MMfile->Get<TH1F>("eff_pt_mu_2015_FAKE2L21_xsecdw_data15-16");
    if (objPtr)h_mu_fake_hf_xsecdw["2015_2016_"+LTDef] = objPtr;
    objPtr = MMfile->Get<TH1F>("eff_pt_mu_2017_FAKE2L21_xsecdw_data17");
    if (objPtr)h_mu_fake_hf_xsecdw["2017_"+LTDef]      = objPtr; 
    objPtr = MMfile->Get<TH1F>("eff_pt_mu_2018_FAKE2L21_xsecdw_data18");
    if (objPtr)h_mu_fake_hf_xsecdw["2018_"+LTDef]      = objPtr;

    objPtr2D = MMfile->Get<TH2F>("eff_pt_eta_mu_2015_REAL2L04_mc16a");
    if (objPtr2D)h_mu_real["2015_2016_"+LTDef] = objPtr2D;
    objPtr2D = MMfile->Get<TH2F>("eff_pt_eta_mu_2017_REAL2L04_mc16d");
    if (objPtr2D)h_mu_real["2017_"+LTDef]      = objPtr2D;
    objPtr2D = MMfile->Get<TH2F>("eff_pt_eta_mu_2018_REAL2L04_mc16e");
    if (objPtr2D)h_mu_real["2018_"+LTDef]      = objPtr2D;
    //}

    for(std::vector< TString >::iterator it = uncert_keys.begin(); it != uncert_keys.end(); ++it){
      if((*it).Contains("EL"))continue;
      printf("Adding %s\n",("2015_2016_"+LTDef+"_"+(*it)).Data());
      objPtr2D = MMfile->Get<TH2F>("frac_pt_eta_mu_2015_REAL2L04_"+(*it)+"_mc16a");
      if (objPtr2D)h_mu_real["2015_2016_"+LTDef+"_"+(*it)] = objPtr2D;
      else cout<<"Could not find "<<"frac_pt_eta_mu_2015_REAL2L04_"+(*it)+"_mc16a"<<endl;
      objPtr2D = MMfile->Get<TH2F>("frac_pt_eta_mu_2017_REAL2L04_"+(*it)+"_mc16d");
      if (objPtr2D)h_mu_real["2017_"+LTDef+"_"+(*it)]      = objPtr2D;
      objPtr2D = MMfile->Get<TH2F>("frac_pt_eta_mu_2018_REAL2L04_"+(*it)+"_mc16e");
      if (objPtr2D)h_mu_real["2018_"+LTDef+"_"+(*it)]      = objPtr2D;
    }
    
  }else{
    //NEW 2L0J Rates
    auto objPtr = MMfile->Get<TH1F>("eff_pt_e_"+LTDef+"_2015_FAKE2L21_data15-16");
    if (objPtr)h_el_fake_hf["2015_2016_"+LTDef] = objPtr;
    objPtr = MMfile->Get<TH1F>("eff_pt_e_"+LTDef+"_2017_FAKE2L21_data17");
    if (objPtr)h_el_fake_hf["2017_"+LTDef]      = objPtr;
    objPtr = MMfile->Get<TH1F>("eff_pt_e_"+LTDef+"_2018_FAKE2L21_data18");
    if (objPtr)h_el_fake_hf["2018_"+LTDef]      = objPtr;


    objPtr = MMfile->Get<TH1F>("eff_pt_e_"+LTDef+"_2015_FAKE2L21_xsecup_data15-16");
    if (objPtr)h_el_fake_hf_xsecup["2015_2016_"+LTDef] = objPtr;
    objPtr = MMfile->Get<TH1F>("eff_pt_e_"+LTDef+"_2017_FAKE2L21_xsecup_data17");
    if (objPtr)h_el_fake_hf_xsecup["2017_"+LTDef]      = objPtr; 
    objPtr = MMfile->Get<TH1F>("eff_pt_e_"+LTDef+"_2018_FAKE2L21_xsecup_data18");
    if (objPtr)h_el_fake_hf_xsecup["2018_"+LTDef]      = objPtr;

    objPtr = MMfile->Get<TH1F>("eff_pt_e_"+LTDef+"_2015_FAKE2L21_xsecdw_data15-16");
    if (objPtr)h_el_fake_hf_xsecdw["2015_2016_"+LTDef] = objPtr;
    objPtr = MMfile->Get<TH1F>("eff_pt_e_"+LTDef+"_2017_FAKE2L21_xsecdw_data17");
    if (objPtr)h_el_fake_hf_xsecdw["2017_"+LTDef]      = objPtr; 
    objPtr = MMfile->Get<TH1F>("eff_pt_e_"+LTDef+"_2018_FAKE2L21_xsecdw_data18");
    if (objPtr)h_el_fake_hf_xsecdw["2018_"+LTDef]      = objPtr;


    objPtr = MMfile->Get<TH1F>("eff_pt_e_"+LTDef+"_2015_FAKE2L23_data15-16");
    if (objPtr)h_el_fake_co["2015_2016_"+LTDef] = objPtr;
    objPtr = MMfile->Get<TH1F>("eff_pt_e_"+LTDef+"_2017_FAKE2L23_data17");
    if (objPtr)h_el_fake_co["2017_"+LTDef]      = objPtr;
    objPtr = MMfile->Get<TH1F>("eff_pt_e_"+LTDef+"_2018_FAKE2L23_data18");
    if (objPtr)h_el_fake_co["2018_"+LTDef]      = objPtr;

    objPtr = MMfile->Get<TH1F>("eff_pt_e_"+LTDef+"_2015_FAKE2L05_mc16a");
    if (objPtr)h_el_fake_lf["2015_2016_"+LTDef] = objPtr;
    objPtr = MMfile->Get<TH1F>("eff_pt_e_"+LTDef+"_2017_FAKE2L05_mc16d");
    if (objPtr)h_el_fake_lf["2017_"+LTDef]      = objPtr;
    objPtr = MMfile->Get<TH1F>("eff_pt_e_"+LTDef+"_2018_FAKE2L05_mc16e");
    if (objPtr)h_el_fake_lf["2018_"+LTDef]      = objPtr;

    //if(!LTDef.Contains("tight_iso"))
    //if(use2Dreal){
    // From data
    // auto objPtr2D = MMfile_real->Get<TH2F>("eff_pt_eta_e_"+LTDef+"_2015_REAL2L01_data15-16");
    // printf("Adding %s\n",("2015_2016_"+LTDef).Data());
    // if (objPtr2D)h_el_real["2015_2016_"+LTDef] = objPtr2D;
    // else cout<<"Could not find "<<"eff_pt_eta_e_"+LTDef+"_2015_REAL2L01_data15-16"<<endl;
    // objPtr2D = MMfile_real->Get<TH2F>("eff_pt_eta_e_"+LTDef+"_2017_REAL2L01_data17");
    // if (objPtr2D) h_el_real["2017_"+LTDef] = objPtr2D;
    // objPtr2D = MMfile_real->Get<TH2F>("eff_pt_eta_e_"+LTDef+"_2018_REAL2L01_data18");
    // if (objPtr2D) h_el_real["2018_"+LTDef] = objPtr2D;
    
    // Original from MC
    auto objPtr2D = MMfile->Get<TH2F>("eff_pt_eta_e_"+LTDef+"_2015_REAL2L04_mc16a");
    printf("Adding %s\n",("2015_2016_"+LTDef).Data());
    if (objPtr2D)h_el_real["2015_2016_"+LTDef] = objPtr2D;
    else cout<<"Could not find "<<"eff_pt_eta_e_"+LTDef+"_2015_REAL2L04_mc16a"<<endl;
    objPtr2D = MMfile->Get<TH2F>("eff_pt_eta_e_"+LTDef+"_2017_REAL2L04_mc16d");
    if (objPtr2D) h_el_real["2017_"+LTDef] = objPtr2D;
    objPtr2D = MMfile->Get<TH2F>("eff_pt_eta_e_"+LTDef+"_2018_REAL2L04_mc16e");
    if (objPtr2D) h_el_real["2018_"+LTDef] = objPtr2D;
      
      
    //}else{
    for (std::string bdt_score : bdt_vec) {
      TString key = bdt_score;	
      auto objPtr = MMfile_real->Get<TH1F>("eff_BDT_"+key+"_e_"+LTDef+"_2015_REAL2L01_data15-16");
      printf("Adding %s\n",("2015_2016_"+LTDef+"_"+key).Data());
      if (objPtr)h_el_real_1D["2015_2016_"+LTDef+"_"+key] = objPtr;
      else cout<<"Could not find "<<"eff_BDT_"+key+"_e_"+LTDef+"_2015_REAL2L01_data15-16"<<endl;
      objPtr = MMfile_real->Get<TH1F>("eff_BDT_"+key+"_e_"+LTDef+"_2017_REAL2L01_data17");
      if (objPtr) h_el_real_1D["2017_"+LTDef+"_"+key] = objPtr;
      objPtr = MMfile_real->Get<TH1F>("eff_BDT_"+key+"_e_"+LTDef+"_2018_REAL2L01_data18");
      printf("---> \t Adding %s\n",("2018_"+LTDef+"_"+key).Data());
      if (objPtr) h_el_real_1D["2018_"+LTDef+"_"+key] = objPtr;
      else cout<<"Could not find "<<"eff_BDT_"+key+"_e_"+LTDef+"_2018_REAL2L01_data18"<<endl;
    }
    //}

    // add systematic histograms
    for(std::vector< TString >::iterator it = uncert_keys.begin(); it != uncert_keys.end(); ++it){
      if((*it).Contains("MUON"))continue;
      printf("Adding %s\n",("2015_2016_"+LTDef+"_"+(*it)).Data());
      objPtr = MMfile->Get<TH1F>("frac_pt_e_"+LTDef+"_2015_FAKE2L05_"+(*it)+"_mc16a");
      if (objPtr)h_el_fake_lf["2015_2016_"+LTDef+"_"+(*it)] = objPtr;
      objPtr = MMfile->Get<TH1F>("frac_pt_e_"+LTDef+"_2017_FAKE2L05_"+(*it)+"_mc16d");
      if (objPtr)h_el_fake_lf["2017_"+LTDef+"_"+(*it)]      = objPtr;
      objPtr = MMfile->Get<TH1F>("frac_pt_e_"+LTDef+"_2018_FAKE2L05_"+(*it)+"_mc16e");
      if (objPtr)h_el_fake_lf["2018_"+LTDef+"_"+(*it)]      = objPtr;
      
      auto objPtr2D = MMfile->Get<TH2F>("frac_pt_eta_e_"+LTDef+"_2015_REAL2L04_"+(*it)+"_mc16a");
      if (objPtr2D)h_el_real["2015_2016_"+LTDef+"_"+(*it)] = objPtr2D;
      else cout<<"Could not find "<<"2015_2016_"+LTDef+"_"+(*it)<<endl;
      objPtr2D = MMfile->Get<TH2F>("frac_pt_eta_e_"+LTDef+"_2017_REAL2L04_"+(*it)+"_mc16d");
      if (objPtr2D) h_el_real["2017_"+LTDef+"_"+(*it)] = objPtr2D;
      objPtr2D = MMfile->Get<TH2F>("frac_pt_eta_e_"+LTDef+"_2018_REAL2L04_"+(*it)+"_mc16e");
      if (objPtr2D) h_el_real["2018_"+LTDef+"_"+(*it)] = objPtr2D; 
    }
    
    //if(!LTDef.EqualTo("tight_iso")){
    objPtr = MMfile->Get<TH1F>("eff_pt_m_"+LTDef+"_2015_FAKE2L21_data15-16");
    if (objPtr)h_mu_fake_hf["2015_2016_"+LTDef] = objPtr;
    objPtr = MMfile->Get<TH1F>("eff_pt_m_"+LTDef+"_2017_FAKE2L21_data17");
    if (objPtr)h_mu_fake_hf["2017_"+LTDef]      = objPtr; 
    objPtr = MMfile->Get<TH1F>("eff_pt_m_"+LTDef+"_2018_FAKE2L21_data18");
    if (objPtr)h_mu_fake_hf["2018_"+LTDef]      = objPtr;

    objPtr = MMfile->Get<TH1F>("eff_pt_m_"+LTDef+"_2015_FAKE2L21_xsecup_data15-16");
    if (objPtr)h_mu_fake_hf_xsecup["2015_2016_"+LTDef] = objPtr;
    objPtr = MMfile->Get<TH1F>("eff_pt_m_"+LTDef+"_2017_FAKE2L21_xsecup_data17");
    if (objPtr)h_mu_fake_hf_xsecup["2017_"+LTDef]      = objPtr; 
    objPtr = MMfile->Get<TH1F>("eff_pt_m_"+LTDef+"_2018_FAKE2L21_xsecup_data18");
    if (objPtr)h_mu_fake_hf_xsecup["2018_"+LTDef]      = objPtr;

    objPtr = MMfile->Get<TH1F>("eff_pt_m_"+LTDef+"_2015_FAKE2L21_xsecdw_data15-16");
    if (objPtr)h_mu_fake_hf_xsecdw["2015_2016_"+LTDef] = objPtr;
    objPtr = MMfile->Get<TH1F>("eff_pt_m_"+LTDef+"_2017_FAKE2L21_xsecdw_data17");
    if (objPtr)h_mu_fake_hf_xsecdw["2017_"+LTDef]      = objPtr; 
    objPtr = MMfile->Get<TH1F>("eff_pt_m_"+LTDef+"_2018_FAKE2L21_xsecdw_data18");
    if (objPtr)h_mu_fake_hf_xsecdw["2018_"+LTDef]      = objPtr;

    //if(use2Dreal){
    // From data
    /**
       auto objPtr2D = MMfile_real->Get<TH2F>("eff_pt_eta_m_"+LTDef+"_2015_REAL2L01_data15-16");
       if (objPtr2D)h_mu_real["2015_2016_"+LTDef] = objPtr2D;
       objPtr2D = MMfile_real->Get<TH2F>("eff_pt_eta_m_"+LTDef+"_2017_REAL2L01_data17");
       if (objPtr2D)h_mu_real["2017_"+LTDef]      = objPtr2D;
       objPtr2D = MMfile_real->Get<TH2F>("eff_pt_eta_m_"+LTDef+"_2018_REAL2L01_data18");
       if (objPtr2D)h_mu_real["2018_"+LTDef]      = objPtr2D;
    */
    // Original from MC
      
    objPtr2D = MMfile->Get<TH2F>("eff_pt_eta_m_"+LTDef+"_2015_REAL2L04_mc16a");
    if (objPtr2D)h_mu_real["2015_2016_"+LTDef] = objPtr2D;
    objPtr2D = MMfile->Get<TH2F>("eff_pt_eta_m_"+LTDef+"_2017_REAL2L04_mc16d");
    if (objPtr2D)h_mu_real["2017_"+LTDef]      = objPtr2D;
    objPtr2D = MMfile->Get<TH2F>("eff_pt_eta_m_"+LTDef+"_2018_REAL2L04_mc16e");
    if (objPtr2D)h_mu_real["2018_"+LTDef]      = objPtr2D;
      
    //}else{
    for(std::string bdt_score : bdt_vec){
      TString key = bdt_score;	
      auto objPtr = MMfile_real->Get<TH1F>("eff_BDT_"+key+"_m_"+LTDef+"_2015_REAL2L01_data15-16");
      printf("Adding %s\n",("2015_2016_"+LTDef+"_"+key).Data());
      if (objPtr)h_mu_real_1D["2015_2016_"+LTDef+"_"+key] = objPtr;
      else cout<<"Could not find "<<"eff_BDT_"+key+"_m_"+LTDef+"_2015_REAL2L01_data15-16"<<endl;
      objPtr = MMfile_real->Get<TH1F>("eff_BDT_"+key+"_m_"+LTDef+"_2017_REAL2L01_data17");
      if (objPtr) h_mu_real_1D["2017_"+LTDef+"_"+key] = objPtr;
      objPtr = MMfile_real->Get<TH1F>("eff_BDT_"+key+"_m_"+LTDef+"_2018_REAL2L01_data18");
      printf("Adding %s\n",("2018_"+LTDef+"_"+key).Data());
      if (objPtr) h_mu_real_1D["2018_"+LTDef+"_"+key] = objPtr;
      else cout<<"Could not find "<<"eff_BDT_"+key+"_m_"+LTDef+"_2018_REAL2L01_data18"<<endl;
    }
    //}
    //}

    for(std::vector< TString >::iterator it = uncert_keys.begin(); it != uncert_keys.end(); ++it){
      if((*it).Contains("EL"))continue;
      printf("Adding %s\n",("2015_2016_"+LTDef+"_"+(*it)).Data());
      auto objPtr2D = MMfile->Get<TH2F>("frac_pt_eta_m_"+LTDef+"_2015_REAL2L04_"+(*it)+"_mc16a");
      if (objPtr2D)h_mu_real["2015_2016_"+LTDef+"_"+(*it)] = objPtr2D;
      objPtr2D = MMfile->Get<TH2F>("frac_pt_eta_m_"+LTDef+"_2017_REAL2L04_"+(*it)+"_mc16d");
      if (objPtr2D)h_mu_real["2017_"+LTDef+"_"+(*it)]      = objPtr2D;
      objPtr2D = MMfile->Get<TH2F>("frac_pt_eta_m_"+LTDef+"_2018_REAL2L04_"+(*it)+"_mc16e");
      if (objPtr2D)h_mu_real["2018_"+LTDef+"_"+(*it)]      = objPtr2D;
    }
    
  }

  for(map< TString, TH2F* >::iterator map_it = h_el_real.begin(); map_it != h_el_real.end(); ++map_it) {
    cout<<"map_it->first = "<<map_it->first<<endl;

    /**
    if(dynamic_cast<TH1F*>(h_el_fake_lf[map_it->first]) == nullptr && dynamic_cast<TH2F*>(h_el_fake_lf[map_it->first]) == nullptr){
      cout<<"INFO \t Could not find fake rates/efficiencies for key "<<map_it->first<<endl;
      if(is18 && (map_it->first).EqualTo("2018"))cout<<"FATAL \t Running over 2018 and no information"<<endl;
      if(is17 && (map_it->first).EqualTo("2017"))cout<<"FATAL \t Running over 2017 and no information"<<endl;
      if(is1516 && (map_it->first).EqualTo("2015_16"))cout<<"FATAL \t Running over 2015/16 and no information"<<endl;
      continue;
    }
    */
    //printf("upper edge of h_el_fake_hf = %.0f\n",h_el_fake_hf[map_it->first]->GetXaxis()->GetBinUpEdge(h_el_fake_hf[map_it->first]->GetNbinsX()));
    if(!(h_el_fake_hf.find(map_it->first) == h_el_fake_hf.end()))
      upper_val_pt["el_fake_hf_"+map_it->first] = h_el_fake_hf[map_it->first]->GetXaxis()->GetBinUpEdge(h_el_fake_hf[map_it->first]->GetNbinsX());
    if(!(h_el_fake_co.find(map_it->first) == h_el_fake_co.end()))
      upper_val_pt["el_fake_co_"+map_it->first] = h_el_fake_co[map_it->first]->GetXaxis()->GetBinUpEdge(h_el_fake_co[map_it->first]->GetNbinsX());
    if(!(h_el_fake_lf.find(map_it->first) == h_el_fake_lf.end()))
      upper_val_pt["el_fake_lf_"+map_it->first] = h_el_fake_lf[map_it->first]->GetXaxis()->GetBinUpEdge(h_el_fake_lf[map_it->first]->GetNbinsX());

   
    cout<<"Doing key "<<map_it->first<<" for:"<<endl;
    cout<<"h_el_fake_hf"<<endl;
    if(!(h_el_fake_hf.find(map_it->first) == h_el_fake_hf.end())) {
      h_el_fake_hf[map_it->first]->SetDirectory(0);
    }
    if(!(h_el_fake_hf_xsecup.find(map_it->first) == h_el_fake_hf_xsecup.end())) {
      cout<<"h_el_fake_hf_xsecup"<<endl;
      h_el_fake_hf_xsecup[map_it->first]->SetDirectory(0);
    }
    if(!(h_el_fake_hf_xsecdw.find(map_it->first) == h_el_fake_hf_xsecdw.end())) {
      cout<<"h_el_fake_hf_xsecdw"<<endl;
      h_el_fake_hf_xsecdw[map_it->first]->SetDirectory(0);
    }

    if(!(h_el_fake_lf.find(map_it->first) == h_el_fake_lf.end())) {
      cout<<"h_el_fake_lf"<<endl;
      h_el_fake_lf[map_it->first]->SetDirectory(0);
    }
    if(!(h_el_fake_co.find(map_it->first) == h_el_fake_co.end())) {
      cout<<"h_el_fake_co"<<endl;
      h_el_fake_co[map_it->first]->SetDirectory(0);
    }
    
    if(!(h_el_real.find(map_it->first) == h_el_real.end())) {
      cout<<"h_el_real "<<map_it->first<<endl;
      h_el_real[map_it->first]->SetDirectory(0);
      cout<<"b"<<endl;
      upper_val_pt["el_real_"+map_it->first]    = h_el_real[map_it->first]->GetXaxis()->GetBinUpEdge(h_el_real[map_it->first]->GetNbinsX());
      cout<<"c"<<endl;
      upper_val_eta["el_real_"+map_it->first]    = h_el_real[map_it->first]->GetYaxis()->GetBinUpEdge(h_el_real[map_it->first]->GetNbinsY());
      cout<<"d"<<endl;
    }
    
    
  }
  cout<<"a"<<endl;
  
  for(map< TString, TH2F* >::iterator map_it = h_mu_real.begin(); map_it != h_mu_real.end(); ++map_it) {
    if(dynamic_cast<TH1F*>(h_mu_real[map_it->first]) == nullptr && dynamic_cast<TH2F*>(h_mu_real[map_it->first]) == nullptr){
      cout<<"INFO \t Could not find fake rates/efficiencies for key "<<map_it->first<<endl;
      if(is18 && (map_it->first).EqualTo("2018"))cout<<"FATAL \t Running over 2018 and no information"<<endl;
      if(is17 && (map_it->first).EqualTo("2017"))cout<<"FATAL \t Running over 2017 and no information"<<endl;
      if(is1516 && (map_it->first).EqualTo("2015_16"))cout<<"FATAL \t Running over 2015/16 and no information"<<endl;
      continue;
    }

    if (!(h_mu_fake_hf.find(map_it->first) == h_mu_fake_hf.end())) {
      cout<<"h_mu_fake_hf"<<endl;
      h_mu_fake_hf[map_it->first]->SetDirectory(0);
      upper_val_pt["mu_fake_hf_"+map_it->first] = h_mu_fake_hf[map_it->first]->GetXaxis()->GetBinUpEdge(h_mu_fake_hf[map_it->first]->GetNbinsX());
    }
    if (!(h_mu_real.find(map_it->first) == h_mu_real.end())) {
      cout<<"h_mu_real = "<<map_it->first<<endl;
      h_mu_real[map_it->first]->SetDirectory(0);
      upper_val_pt["mu_real_"+map_it->first]    = h_mu_real[map_it->first]->GetXaxis()->GetBinUpEdge(h_mu_real[map_it->first]->GetNbinsX());
      upper_val_eta["mu_real_"+map_it->first]    = h_mu_real[map_it->first]->GetYaxis()->GetBinUpEdge(h_mu_real[map_it->first]->GetNbinsY());
    }
    
    if(!(h_mu_fake_hf_xsecup.find(map_it->first) == h_mu_fake_hf_xsecup.end())) {
      cout<<"h_mu_fake_hf_xsecup"<<endl;
      h_mu_fake_hf_xsecup[map_it->first]->SetDirectory(0);
    }
    if(!(h_mu_fake_hf_xsecdw.find(map_it->first) == h_mu_fake_hf_xsecdw.end())) {
      cout<<"h_mu_fake_hf_xsecdw"<<endl;
      h_mu_fake_hf_xsecdw[map_it->first]->SetDirectory(0);
    }
    
    cout<<"aa"<<endl;

  }

  cout<<"bb"<<endl;
  
  delete MMfile;
  delete MMfile_frac;
}



// For those who doesn't need the TLorentzVectors...
int GetFakeWeight::setFrac(float var_dep, float mll, int nJet, int n_bjet, double zmass, float met_Sign){
  TLorentzVector a;
  TLorentzVector b;
  TLorentzVector c;
  return setFrac(var_dep, mll, nJet, n_bjet, zmass, met_Sign, a, b, c);
}

int GetFakeWeight::setFrac(float var_dep, float mll, int nJet, int n_bjet, double zmass, float met_Sign, TLorentzVector met_tlv, TLorentzVector lep1, TLorentzVector lep2){

  TString this_typ_old = this_typ;
  
  if(this_typ.EqualTo("MMOS") or this_typ.EqualTo("MMSS")){
    
    hf_frac_nom = 1.0;
    lf_frac_nom = 0.0;
    co_frac_nom = 0.0;
    
    hf_frac_max = 1.0;
    lf_frac_max = 0.0;
    co_frac_max = 0.0;
    
    hf_frac_min = 1.0;
    lf_frac_min = 0.0;
    co_frac_min = 0.0;

    return 1;
  }


  if(verbosity)cout<<"Type is "<<this_typ.Data()<<endl;
  
  std::vector<double> hf;
  std::vector<double> lf;
  std::vector<double> co;

  std::vector<double> hf_up;
  std::vector<double> lf_up;
  std::vector<double> co_up;

  std::vector<double> hf_dw;
  std::vector<double> lf_dw;
  std::vector<double> co_dw;

  
  std::vector<int> frac_regions;
  
  int mt2bin = 1;
  if(frac_var.EqualTo("mT2")){
    if(var_dep >= 0 && var_dep < 20)mt2bin = 1;
    else if(var_dep >= 20 && var_dep < 40)mt2bin = 2;
    else if(var_dep >= 40 && var_dep < 60)mt2bin = 3;
    else if(var_dep >= 60 && var_dep < 100)mt2bin = 4;
    else if(var_dep >= 100)mt2bin = 5;
  }else if(frac_var.EqualTo("metsig")){
    if(var_dep >= 0 && var_dep < 2)mt2bin = 2;
    else if(var_dep >= 2 && var_dep < 4)mt2bin = 2;
    else if(var_dep >= 4 && var_dep < 6)mt2bin = 3;
    else if(var_dep >= 6 && var_dep < 8)mt2bin = 4;
    else if(var_dep >= 8 && var_dep < 10)mt2bin = 5;
    else if(var_dep >= 10)mt2bin = 6;
  }else{
    printf("FATAL \t Does not regognize variable %s\n",frac_var.Data());
    return -1;
  }
  if(verbosity)cout<<frac_var.Data()<<" bin "<<mt2bin<< " for value "<< var_dep << endl;

  //for(int i = 4; i<=11; i++){
  if(this_ana.EqualTo("2L2J")){
    frac_regions.push_back(4);
    if(fabs((mll-zmass)) < 20.0 && n_bjet <= 1){frac_regions.push_back(10);}
    if((fabs((mll-zmass)) < 20.0 && n_bjet >= 2)){frac_regions.push_back(11);}
    if((fabs((mll-zmass)) < 10.0 && n_bjet >= 1 && nJet >= 2)){frac_regions.push_back(9);}
    if((fabs((mll-zmass)) < 10.0 && n_bjet == 0 && nJet == 2)){frac_regions.push_back(8);}
    if(((mll > 12 && mll < 71) && n_bjet == 0 && nJet >= 2)){frac_regions.push_back(7);}
    if((mll > 12 && fabs((mll-zmass)) < 10.0 && nJet >= 2)){frac_regions.push_back(6);}
    if((mll > 12 && nJet >= 2)){frac_regions.push_back(5);}
  }else{
    /**
    frac_regions.push_back(14);
    this_typ = "ALL";
    
    if(fabs((mll-zmass)) < 10.0){
      frac_regions.push_back(2);
      if(nJet <= 2)frac_regions.push_back(3);
    }
    */
    if(((this_typ.Contains("MM") || this_typ.Contains("EE")) && fabs((mll-zmass)) > 15.0) ||
       (this_typ.Contains("EM") || (this_typ.Contains("ME")))){
      if(n_bjet >= 1){
	frac_regions.push_back(13);
	if(nJet == 0)frac_regions.push_back(12);
      }else{
	frac_regions.push_back(4);
	if(nJet == 0)frac_regions.push_back(11);
	//if(nJet == 1)frac_regions.push_back(12);
	if(lep1.DeltaPhi(lep2) > 2.2)frac_regions.push_back(5);
	if(TMath::ATan(fabs(lep1.Eta()-lep2.Eta())/2.)<0.2)frac_regions.push_back(6);
	if(lep1.Pt() > 140 && lep2.Pt() > 20)frac_regions.push_back(7);
	if(mll > 60)frac_regions.push_back(8);
	if((met_tlv+lep1+lep2).Pt() < 5)frac_regions.push_back(9);
	if((met_tlv.DeltaPhi(lep1)>2.2))frac_regions.push_back(10);
      }
    }
    
  }

  // In case no evens passed use preselection, event is likely to not be used anyway.
  if(frac_regions.size() == 0)frac_regions.push_back(4);
  
  for(std::vector< int >::iterator it = frac_regions.begin(); it != frac_regions.end(); it++){
    //if((*it) == 4 || (*it) == 12)continue;

    if(verbosity){
      printf("Using region %i\n",(*it));
      cout<<"Getting frac "<<Form("%s_%02i_%s",this_typ.Data(),(*it),this_ana.Data())<<endl;
    }

    if ( h_el_frac_hf.find(Form("%s_%02i_%s",this_typ.Data(),(*it),this_ana.Data())) == h_el_frac_hf.end() ) {
      //      cout<<"Could not find frac for "<<Form("%s_%02i_%s",this_typ.Data(),(*it),this_ana.Data())<<endl;
      continue;
    }
    hf.push_back(h_el_frac_hf[Form("%s_%02i_%s",this_typ.Data(),(*it),this_ana.Data())]->GetBinContent(mt2bin)+h_el_frac_cf[Form("%s_%02i_%s",this_typ.Data(),(*it),this_ana.Data())]->GetBinContent(mt2bin));
    lf.push_back(h_el_frac_lf[Form("%s_%02i_%s",this_typ.Data(),(*it),this_ana.Data())]->GetBinContent(mt2bin));
    co.push_back(h_el_frac_co[Form("%s_%02i_%s",this_typ.Data(),(*it),this_ana.Data())]->GetBinContent(mt2bin)/**+h_el_frac_em[Form("%s_%02i_%s",this_typ.Data(),(*it),this_ana.Data())]->GetBinContent(mt2bin)*/);

    double sum = (hf.back()+lf.back()+co.back());
    if(fabs(sum-1) > 0.05){
      hf.back() *= 1.0/sum;
      lf.back() *= 1.0/sum;
      co.back() *= 1.0/sum;
    }

    hf_up.push_back(hf.back()+(h_el_frac_hf[Form("%s_%02i_%s",this_typ.Data(),(*it),this_ana.Data())]->GetBinError(mt2bin)+
			       h_el_frac_cf[Form("%s_%02i_%s",this_typ.Data(),(*it),this_ana.Data())]->GetBinError(mt2bin)));
    lf_up.push_back(lf.back()+ h_el_frac_lf[Form("%s_%02i_%s",this_typ.Data(),(*it),this_ana.Data())]->GetBinError(mt2bin));
    co_up.push_back(co.back()+(h_el_frac_co[Form("%s_%02i_%s",this_typ.Data(),(*it),this_ana.Data())]->GetBinError(mt2bin)/**+
															     h_el_frac_em[Form("%s_%02i_%s",this_typ.Data(),(*it),this_ana.Data())]->GetBinError(mt2bin)*/));

    hf_dw.push_back(hf.back()-(h_el_frac_hf[Form("%s_%02i_%s",this_typ.Data(),(*it),this_ana.Data())]->GetBinError(mt2bin)+  
		   	       h_el_frac_cf[Form("%s_%02i_%s",this_typ.Data(),(*it),this_ana.Data())]->GetBinError(mt2bin)));
    lf_dw.push_back(lf.back()- h_el_frac_lf[Form("%s_%02i_%s",this_typ.Data(),(*it),this_ana.Data())]->GetBinError(mt2bin));  
    co_dw.push_back(co.back()-(h_el_frac_co[Form("%s_%02i_%s",this_typ.Data(),(*it),this_ana.Data())]->GetBinError(mt2bin)/**+  
															     h_el_frac_em[Form("%s_%02i_%s",this_typ.Data(),(*it),this_ana.Data())]->GetBinError(mt2bin)*/));

    
    hf.back() = hf.back() < 0.0 ? 0.0 : hf.back();
    lf.back() = lf.back() < 0.0 ? 0.0 : lf.back();
    co.back() = co.back() < 0.0 ? 0.0 : co.back();
    
    hf_up.back() = hf_up.back() < 0.0 ? 0.0 : hf_up.back();
    lf_up.back() = lf_up.back() < 0.0 ? 0.0 : lf_up.back();
    co_up.back() = co_up.back() < 0.0 ? 0.0 : co_up.back();

    hf_dw.back() = hf_dw.back() < 0.0 ? 0.0 : hf_dw.back();
    lf_dw.back() = lf_dw.back() < 0.0 ? 0.0 : lf_dw.back();
    co_dw.back() = co_dw.back() < 0.0 ? 0.0 : co_dw.back();    
    
  }
  
  hf_frac_max = *max_element(std::begin(hf), std::end(hf));
  lf_frac_max = *max_element(std::begin(lf), std::end(lf));
  co_frac_max = *max_element(std::begin(co), std::end(co));

  hf_frac_min = *min_element(std::begin(hf), std::end(hf));
  lf_frac_min = *min_element(std::begin(lf), std::end(lf));
  co_frac_min = *min_element(std::begin(co), std::end(co));

  hf_frac_nom = accumulate( hf.begin(), hf.end(), 0.0)/hf.size();
  lf_frac_nom = accumulate( lf.begin(), lf.end(), 0.0)/lf.size();
  co_frac_nom = accumulate( co.begin(), co.end(), 0.0)/co.size();

  hf_frac_statdw = accumulate( hf_dw.begin(), hf_dw.end(), 0.0)/hf_dw.size();
  lf_frac_statdw = accumulate( lf_dw.begin(), lf_dw.end(), 0.0)/lf_dw.size();
  co_frac_statdw = accumulate( co_dw.begin(), co_dw.end(), 0.0)/co_dw.size();

  hf_frac_statup = accumulate( hf_up.begin(), hf_up.end(), 0.0)/hf_up.size();
  lf_frac_statup = accumulate( lf_up.begin(), lf_up.end(), 0.0)/lf_up.size();
  co_frac_statup = accumulate( co_up.begin(), co_up.end(), 0.0)/co_up.size();

  /**
  if((hf_frac_nom + lf_frac_nom + co_frac_nom) != 1){
    hf_frac_nom += (1-(hf_frac_nom + lf_frac_nom + co_frac_nom));
    hf_frac_max += (1-(hf_frac_nom + lf_frac_nom + co_frac_nom));
    hf_frac_statup  += (1-(hf_frac_nom + lf_frac_nom + co_frac_nom));

    hf_frac_min += (1-(hf_frac_nom + lf_frac_nom + co_frac_nom));
    hf_frac_statdw  += (1-(hf_frac_nom + lf_frac_nom + co_frac_nom));
    
  }
  */
  if(verbosity){
    printf("Systematics: \n");
    printf("HF frac = %.2f ; max = %.2f; min = %.2f\n",hf_frac_nom,hf_frac_max,hf_frac_min);
    printf("LF frac = %.2f ; max = %.2f; min = %.2f\n",lf_frac_nom,lf_frac_max,lf_frac_min);
    printf("CO frac = %.2f ; max = %.2f; min = %.2f\n",co_frac_nom,co_frac_max,co_frac_min);
    printf("SUM frac = %.2f (nom), %.2f (max), %.2f (min)\n",hf_frac_nom+lf_frac_nom+co_frac_nom,hf_frac_max+lf_frac_max+co_frac_max,hf_frac_min+lf_frac_min+co_frac_min);
    printf("Statistsics: \n");
    printf("HF frac = %.2f ; max = %.2f; min = %.2f\n",hf_frac_nom,hf_frac_statup,hf_frac_statdw);
    printf("LF frac = %.2f ; max = %.2f; min = %.2f\n",lf_frac_nom,lf_frac_statup,lf_frac_statdw);
    printf("CO frac = %.2f ; max = %.2f; min = %.2f\n",co_frac_nom,co_frac_statup,co_frac_statdw);
    printf("SUM frac = %.2f (nom), %.2f (max), %.2f (min)\n",hf_frac_nom+lf_frac_nom+co_frac_nom,hf_frac_statup+lf_frac_statup+co_frac_statup,hf_frac_statdw+lf_frac_statdw+co_frac_statdw);
  }

  this_typ = this_typ_old;
  
  return 1;
  
}

void GetFakeWeight::setVariables(TString year, TString ana, TString eventStr){

  this_year = year;

  this_ana = ana;

  this_typ = eventStr;

}


double GetFakeWeight::getFake(float pT, float eta, int leptype, TString triggermatch, uncertainty u){

  std::pair <Double_t,Double_t> upcut;
  Double_t pT_hf = pT;
  Double_t pT_lf = pT;
  Double_t pT_co = pT;
  TString key;
  
  if(triggermatch.Length() <= 0)triggermatch = this_ana;

  key = this_year+"_"+triggermatch;
  
  
  
  if(verbosity){
    cout<<"-----------------------------------------------------------------------------------------"<<endl;
    cout<<"FAKE: Getting "<<key<< " for "<< (leptype == 1 ? "el" : "mu") << " with pT "<< pT<< endl;
    //cout << "hf_frac_nom = "<<hf_frac_nom<<endl;
    //cout << "lf_frac_nom = "<<lf_frac_nom<<endl;
    //cout << "co_frac_nom = "<<co_frac_nom<<endl;

    // for (auto const& pair: h_el_fake_hf) {
    //   std::cout << "{" << pair.first << ": " << pair.second << "}\n";
    // }
  }
  
  switch(leptype){
  case 1:

    upcut = getUpperCut("el_fake_hf_"+key);
    pT_hf = pT < upcut.first ? pT : upcut.first;
    upcut = getUpperCut("el_fake_lf_"+key);
    pT_lf = pT < upcut.first ? pT : upcut.first;
    upcut = getUpperCut("el_fake_co_"+key);
    pT_co = pT < upcut.first ? pT : upcut.first;
    
    
    if(verbosity){
      
      printf("Electron: pT = %.2f (hf: %.2f, lf: %.2f, co: %.2f), eta = %.2f, f = %.4f\n",pT,pT_hf,pT_lf,pT_co,eta,
	     getFakeUnc(key,leptype,pT_hf,pT_lf,pT_co,h_el_fake_hf,h_el_fake_lf,h_el_fake_co,h_el_fake_hf_xsecup,h_el_fake_hf_xsecdw,u));
    }

    //    getFakeStatUp(key, pT_hf, pT_lf, pT_co)
    
    return getFakeUnc(key,leptype,pT_hf,pT_lf,pT_co,h_el_fake_hf,h_el_fake_lf,h_el_fake_co,h_el_fake_hf_xsecup,h_el_fake_hf_xsecdw,u);
    break;
  case 2:

    upcut = getUpperCut("mu_fake_hf_"+key);
    pT_hf = pT < upcut.first ? pT : upcut.first;
    
    if(verbosity)printf("Muon:     pT = %.2f (%.2f), eta = %.2f, f = %.4f\n",pT,pT_hf,eta,getFakeUnc(key,leptype,pT_hf,pT_lf,pT_co,h_mu_fake_hf,h_mu_fake_hf,h_mu_fake_hf,h_mu_fake_hf_xsecup,h_mu_fake_hf_xsecdw,u));
    return getFakeUnc(key,leptype,pT_hf,pT_lf,pT_co,h_mu_fake_hf,h_mu_fake_hf,h_mu_fake_hf,h_mu_fake_hf_xsecup,h_mu_fake_hf_xsecdw,u);
    //return h_mu_fake_hf[key]->GetBinContent(h_mu_fake_hf[key]->FindBin(pT_hf));
    break;
  default:
    cout<<"Could not find lepton type "<<leptype<<endl;
    return -1;
  }
}

double GetFakeWeight::getRealUnc(TString key, int leptype, double pT, double eta, th2Map h_real, uncertainty u, th1Map h_real_1D, TString BDTname, double BDTvalue){

  double systw = 1.0;
  TString systematic;// = leptype == 1 ? "EL_"+unckey : "MUON_"+unckey;
  
  switch(u){
  case NOM:
    if(BDTvalue < 0){
      return h_real[key]->GetBinContent(h_real[key]->FindBin(pT,eta));
    }else{
      key += "_"+BDTname;
      //cout<<"Getting key"<< key.Data()<<endl;
      return h_real_1D[key]->GetBinContent(h_real_1D[key]->FindBin(BDTvalue));
    }
    break;
  case STATUP:
    return h_real[key]->GetBinContent(h_real[key]->FindBin(pT,eta))+h_real[key]->GetBinError(h_real[key]->FindBin(pT,eta));
    break;
  case STATDW:
    return h_real[key]->GetBinContent(h_real[key]->FindBin(pT,eta))-h_real[key]->GetBinError(h_real[key]->FindBin(pT,eta));
    break;
  case SYSTUP:
    systematic = (leptype == 1 ? "EL_"+uncertainty_keys.at(SYSTUP) : "MUON_"+uncertainty_keys.at(SYSTUP));
    // If not available return nominal
    if(h_real.find(key+"_"+systematic) == h_real.end()){
      printf("ERROR \t REAL: Could not find histogram for systematics %s_%s\n",key.Data(),systematic.Data());
      systw = 1.0;
    }else{
      systw = h_real[key+"_"+systematic]->GetBinContent(h_real[key+"_"+systematic]->FindBin(pT,eta));
    }
    return systw*h_real[key]->GetBinContent(h_real[key]->FindBin(pT,eta));;
    
    break;
  case SYSTDW:
    systematic = leptype == 1 ? "EL_"+uncertainty_keys.at(SYSTDW) : "MUON_"+uncertainty_keys.at(SYSTDW);
    //systematic += "_1down";
    // If not available return nominal
    if(h_real.find(key+"_"+systematic) == h_real.end()){
      printf("ERROR \t REAL: Could not find histogram for systematics %s_%s and lep %i\n",key.Data(),systematic.Data(),leptype);
      systw = 1.0;
    }else{
      systw = h_real[key+"_"+systematic]->GetBinContent(h_real[key+"_"+systematic]->FindBin(pT,eta));
    }
    return systw*h_real[key]->GetBinContent(h_real[key]->FindBin(pT,eta));
    
    break;
  case TRIGUP:
    systematic = leptype == 1 ? "EL_"+uncertainty_keys.at(TRIGUP) : "MUON_"+uncertainty_keys.at(TRIGUP);
    if(leptype == 1)systematic = systematic.ReplaceAll("TRIG","Trigger");
    //systematic += "_1up";
    // If not available return nominal
    if(h_real.find(key+"_"+systematic) == h_real.end()){
      printf("ERROR \t REAL: Could not find histogram for systematics %s_%s and lep %i\n",key.Data(),systematic.Data(),leptype);
      systw = 1.0;
    }else{
      systw = h_real[key+"_"+systematic]->GetBinContent(h_real[key+"_"+systematic]->FindBin(pT,eta));
    }
    return systw*h_real[key]->GetBinContent(h_real[key]->FindBin(pT,eta));;
    
    break;
  case TRIGDW:
    systematic = leptype == 1 ? "EL_"+uncertainty_keys.at(TRIGDW) : "MUON_"+uncertainty_keys.at(TRIGDW);
    if(leptype == 1)systematic = systematic.ReplaceAll("TRIG","Trigger");
    //systematic += "_1down";
    // If not available return nominal
    if(h_real.find(key+"_"+systematic) == h_real.end()){
      printf("ERROR \t REAL: Could not find histogram for systematics %s_%s\n",key.Data(),systematic.Data());
      systw = 1.0;
    }else{
      systw = h_real[key+"_"+systematic]->GetBinContent(h_real[key+"_"+systematic]->FindBin(pT,eta));
    }
    return systw*h_real[key]->GetBinContent(h_real[key]->FindBin(pT,eta));
    
    break;
  default:
    return h_real[key]->GetBinContent(h_real[key]->FindBin(pT,eta));
  }
}

double GetFakeWeight::getWeightedFake(double hf_nom, double lf_nom, double co_nom, double hf_var, double lf_var, double co_var, th1Map h_hf, th1Map h_lf, th1Map h_co, double pT_hf, double pT_lf, double pT_co, TString key, TString syst) {

  double ch  = 0.0;
  std::vector<double> maxf;
  double newlf;
  double newhf;
  double newco;
  double sf = 1.0;

  

  bool ismax = hf_var > hf_nom;

  if(ismax  && (lf_var < lf_nom || co_var < co_nom))printf("Is max, but lf: nom =  %.2f > var = %.2f and co: nom = %.2f > var = %.2f\n",lf_nom,lf_var,co_nom,co_var);
  if(!ismax && (lf_var > lf_nom || co_var > co_nom))printf("Is min, but lf: nom =  %.2f < var = %.2f and co: nom = %.2f < var = %.2f\n",lf_nom,lf_var,co_nom,co_var);
  
  ch = ismax ? hf_var-hf_nom : hf_nom-hf_var;
  
  sf = (lf_nom+co_nom);
  newlf = ismax ? (lf_nom-(lf_nom/sf*ch)) : (lf_nom+(lf_nom/sf*ch));
  newco = ismax ? (co_nom-(co_nom/sf*ch)) : (co_nom+(co_nom/sf*ch));

  if(fabs((hf_var+newlf+newco)-1) > 0.05){
    printf("ERROR \t hf = %.4f, lf = %.4f, co = %.4f. SUM = %.2f\n",hf_var,newlf,newco,(hf_var+newlf+newco));
    cout<<"Systematics is "<<syst<<endl;
    printf("%20s : hf = %.4f, lf = %.4f, co = %.4f. SUM = %.2f\n","NOMINAL",hf_nom,lf_nom,co_nom,(hf_nom+lf_nom+co_nom));
    printf("%20s : hf = %.4f, lf = %.4f, co = %.4f. SUM = %.2f\n",syst.Data(),hf_var,lf_var,co_var,(hf_var+lf_var+co_var));
  }

  if(!isnan((hf_var*h_hf[key]->GetBinContent(h_hf[key]->FindBin(pT_hf)) +
	     newlf*h_lf[key]->GetBinContent(h_lf[key]->FindBin(pT_lf)) +
	     newco*h_co[key]->GetBinContent(h_co[key]->FindBin(pT_co))))){
    maxf.push_back((hf_var*h_hf[key]->GetBinContent(h_hf[key]->FindBin(pT_hf)) +
		    newlf*h_lf[key]->GetBinContent(h_lf[key]->FindBin(pT_lf)) +
		    newco*h_co[key]->GetBinContent(h_co[key]->FindBin(pT_co))));
  }
  ch = ismax ? lf_var-lf_nom : lf_nom-lf_var;
  sf = (hf_nom+co_nom);
  
  newhf = ismax ? (hf_nom-(hf_nom/sf*ch)) : (hf_nom+(hf_nom/sf*ch));
  newco = ismax ? (co_nom-(co_nom/sf*ch)) : (co_nom+(co_nom/sf*ch));
  
  if(!isnan((newhf*h_hf[key]->GetBinContent(h_hf[key]->FindBin(pT_hf)) +
	     lf_var*h_lf[key]->GetBinContent(h_lf[key]->FindBin(pT_lf)) +
	     newco*h_co[key]->GetBinContent(h_co[key]->FindBin(pT_co))))){
    maxf.push_back((newhf*h_hf[key]->GetBinContent(h_hf[key]->FindBin(pT_hf)) +
		    lf_var*h_lf[key]->GetBinContent(h_lf[key]->FindBin(pT_lf)) +
		    newco*h_co[key]->GetBinContent(h_co[key]->FindBin(pT_co))));
  }
  
  if(fabs((lf_var+newhf+newco)-1) > 0.05){
    printf("ERROR \t hf = %.4f, lf = %.4f, co = %.4f. SUM = %.2f\n",newhf,lf_var,newco,(newhf+lf_var+newco));
    cout<<"Systematics is "<<syst<<endl;
    printf("%20s : hf = %.4f, lf = %.4f, co = %.4f. SUM = %.2f\n","NOMINAL",hf_nom,lf_nom,co_nom,(hf_nom+lf_nom+co_nom));
    printf("%20s : hf = %.4f, lf = %.4f, co = %.4f. SUM = %.2f\n",syst.Data(),hf_var,lf_var,co_var,(hf_var+lf_var+co_var));
  }
    
  ch = ismax ? co_var-co_nom : co_nom-co_var;
  sf = (hf_nom+lf_nom);
  newhf = ismax ? (hf_nom-(hf_nom/sf*ch)) : (hf_nom+(hf_nom/sf*ch));
  newlf = ismax ? (lf_nom-(lf_nom/sf*ch)) : (lf_nom+(lf_nom/sf*ch));

  if(fabs((co_var+newhf+newlf)-1) > 0.05){
    printf("ERROR \t hf = %.4f, lf = %.4f, co = %.4f. SUM = %.2f\n",newhf,newlf,co_var,(newhf+newlf+co_var));
    cout<<"Systematics is "<<syst<<endl;
    printf("%20s : hf = %.4f, lf = %.4f, co = %.4f. SUM = %.2f\n","NOMINAL",hf_nom,lf_nom,co_nom,(hf_nom+lf_nom+co_nom));
    printf("%20s : hf = %.4f, lf = %.4f, co = %.4f. SUM = %.2f\n",syst.Data(),hf_var,lf_var,co_var,(hf_var+lf_var+co_var));
  }

  if(!isnan((newhf*h_hf[key]->GetBinContent(h_hf[key]->FindBin(pT_hf)) +
	     newlf*h_lf[key]->GetBinContent(h_lf[key]->FindBin(pT_lf)) +
	     co_var*h_co[key]->GetBinContent(h_co[key]->FindBin(pT_co))))){
    
    maxf.push_back((newhf*h_hf[key]->GetBinContent(h_hf[key]->FindBin(pT_hf)) +
		    newlf*h_lf[key]->GetBinContent(h_lf[key]->FindBin(pT_lf)) +
		    co_var*h_co[key]->GetBinContent(h_co[key]->FindBin(pT_co))));
  }

  if(verbosity){
    printf("Size of maxf = %i\n",maxf.size());
    for(const auto& value: maxf) {
      std::cout << value << "\n";
    }
    printf("max = %.2f\n",*max_element(std::begin(maxf), std::end(maxf)));
    printf("min = %.2f\n",*min_element(std::begin(maxf), std::end(maxf)));
  }
  
  if(ismax)return *max_element(std::begin(maxf), std::end(maxf));
  else return *min_element(std::begin(maxf), std::end(maxf));


}

double GetFakeWeight::getFakeUnc(TString key, int leptype, double pT_hf, double pT_lf, double pT_co, th1Map h_hf, th1Map h_lf, th1Map h_co, th1Map h_hf_xsecup, th1Map h_hf_xsecdw, uncertainty u){
  double systw = 1.0;

  TString systematic;// = leptype == 1 ? "EL_"+unckey : "MUON_"+unckey;

  double hf_temp_frac_nom = leptype == 1 ? hf_frac_nom : 1.0;
  double lf_temp_frac_nom = leptype == 1 ? lf_frac_nom : 0.0;
  double co_temp_frac_nom = leptype == 1 ? co_frac_nom : 0.0;

  double hf_temp_frac_max = leptype == 1 ? hf_frac_max : 1.0;
  double lf_temp_frac_max = leptype == 1 ? lf_frac_max : 0.0;
  double co_temp_frac_max = leptype == 1 ? co_frac_max : 0.0;

  double hf_temp_frac_min = leptype == 1 ? hf_frac_min : 1.0;
  double lf_temp_frac_min = leptype == 1 ? lf_frac_min : 0.0;
  double co_temp_frac_min = leptype == 1 ? co_frac_min : 0.0;

  double hf_temp_frac_statdw = leptype == 1 ? hf_frac_statdw : 1.0;
  double lf_temp_frac_statdw = leptype == 1 ? lf_frac_statdw : 0.0;
  double co_temp_frac_statdw = leptype == 1 ? co_frac_statdw : 0.0;
                    		  
  double hf_temp_frac_statup = leptype == 1 ? hf_frac_statup : 1.0;
  double lf_temp_frac_statup = leptype == 1 ? lf_frac_statup : 0.0;
  double co_temp_frac_statup = leptype == 1 ? co_frac_statup : 0.0;

  double nominal = (hf_temp_frac_nom*h_hf[key]->GetBinContent(h_hf[key]->FindBin(pT_hf)) +
		    lf_temp_frac_nom*h_lf[key]->GetBinContent(h_lf[key]->FindBin(pT_lf)) +
		    co_temp_frac_nom*h_co[key]->GetBinContent(h_co[key]->FindBin(pT_co)));

  double xsecup_temp;
  double xsecdw_temp;
  
  if(verbosity)printf("Using hf = %.2f, lf = %.2f, co = %.2f\n",hf_temp_frac_nom,lf_temp_frac_nom,co_temp_frac_nom);
  
  switch(u){
  case NOM:

    if(verbosity)printf("NOM \t f_hf = %.4f, f_lf = %.4f, f_co = %.4f\n",h_hf[key]->GetBinContent(h_hf[key]->FindBin(pT_hf)),h_lf[key]->GetBinContent(h_lf[key]->FindBin(pT_lf)),h_co[key]->GetBinContent(h_co[key]->FindBin(pT_co)));
    
    
    return (hf_temp_frac_nom*h_hf[key]->GetBinContent(h_hf[key]->FindBin(pT_hf)) +
	    lf_temp_frac_nom*h_lf[key]->GetBinContent(h_lf[key]->FindBin(pT_lf)) +
	    co_temp_frac_nom*h_co[key]->GetBinContent(h_co[key]->FindBin(pT_co)));
    break;
  case STATUP:
    return (hf_temp_frac_nom*(h_hf[key]->GetBinContent(h_hf[key]->FindBin(pT_hf))+h_hf[key]->GetBinError(h_hf[key]->FindBin(pT_hf))) +
	    lf_temp_frac_nom*(h_lf[key]->GetBinContent(h_lf[key]->FindBin(pT_lf))+h_lf[key]->GetBinError(h_lf[key]->FindBin(pT_lf))) +
	    co_temp_frac_nom*(h_co[key]->GetBinContent(h_co[key]->FindBin(pT_co))+h_co[key]->GetBinError(h_co[key]->FindBin(pT_co))));
    break;
  case STATDW:
    return (hf_temp_frac_nom*(h_hf[key]->GetBinContent(h_hf[key]->FindBin(pT_hf))-h_hf[key]->GetBinError(h_hf[key]->FindBin(pT_hf))) +
	    lf_temp_frac_nom*(h_lf[key]->GetBinContent(h_lf[key]->FindBin(pT_lf))-h_lf[key]->GetBinError(h_lf[key]->FindBin(pT_lf))) +
	    co_temp_frac_nom*(h_co[key]->GetBinContent(h_co[key]->FindBin(pT_co))-h_co[key]->GetBinError(h_co[key]->FindBin(pT_co))));
    break;
  case WEIGHTUP:
    if(leptype != 1)return nominal;

    return getWeightedFake(hf_temp_frac_nom, lf_temp_frac_nom, co_temp_frac_nom, hf_temp_frac_statup,  lf_temp_frac_statup, co_temp_frac_statup, h_hf, h_lf, h_co, pT_hf, pT_lf, pT_co,key,uncertainty_keys.at(WEIGHTUP));
    
    // return (hf_temp_frac_statup*h_hf[key]->GetBinContent(h_hf[key]->FindBin(pT_hf)) +
    // 	    lf_temp_frac_statup*h_lf[key]->GetBinContent(h_lf[key]->FindBin(pT_lf)) +
    // 	    co_temp_frac_statup*h_co[key]->GetBinContent(h_co[key]->FindBin(pT_co)));
    break;
  case WEIGHTDW:
    if(leptype != 1)return nominal;

    return getWeightedFake(hf_temp_frac_nom, lf_temp_frac_nom, co_temp_frac_nom, hf_temp_frac_statdw,  lf_temp_frac_statdw, co_temp_frac_statdw, h_hf, h_lf, h_co, pT_hf, pT_lf, pT_co,key,uncertainty_keys.at(WEIGHTDW));
    
    // return (hf_temp_frac_statdw*h_hf[key]->GetBinContent(h_hf[key]->FindBin(pT_hf)) +
    // 	    lf_temp_frac_statdw*h_lf[key]->GetBinContent(h_lf[key]->FindBin(pT_lf)) +
    // 	    co_temp_frac_statdw*h_co[key]->GetBinContent(h_co[key]->FindBin(pT_co)));
    break;
  case WEIGHTSYSTUP:
    if(leptype != 1)return nominal;
    //if(verbosity)printf("WEIGHTSYSTUP \t hf = %.4f, lf = %.4f, co = %.4f\n",hf_temp_frac_max,lf_temp_frac_max,co_temp_frac_max);
    
    return getWeightedFake(hf_temp_frac_nom, lf_temp_frac_nom, co_temp_frac_nom, hf_temp_frac_max,  lf_temp_frac_max, co_temp_frac_max, h_hf, h_lf, h_co, pT_hf, pT_lf, pT_co,key,uncertainty_keys.at(WEIGHTSYSTUP));
      
      //return *max_element(std::begin(maxf), std::end(maxf));
    //if(verbosity)printf("WEIGHTSYSTUP \t hf = %.4f, lf = %.4f, co = %.4f. SUM = %.2f\n",hf_temp_frac_max,(lf_temp_frac_nom-(lf_temp_frac_nom*ch)),(co_temp_frac_nom-(co_temp_frac_nom*ch)),hf_temp_frac_max+(lf_temp_frac_nom-(lf_temp_frac_nom*ch))+(co_temp_frac_nom-(co_temp_frac_nom*ch)));
    
    //#  (lf_temp_frac_nom-(lf_temp_frac_nom*ch))*h_lf[key]->GetBinContent(h_lf[key]->FindBin(pT_lf)) +
    //  (co_temp_frac_nom-(co_temp_frac_nom*ch))*h_co[key]->GetBinContent(h_co[key]->FindBin(pT_co)
    
    //return (hf_temp_frac_max*h_hf[key]->GetBinContent(h_hf[key]->FindBin(pT_hf)) +
    //	    lf_temp_frac_max*h_lf[key]->GetBinContent(h_lf[key]->FindBin(pT_lf)) +
    //	    co_temp_frac_max*h_co[key]->GetBinContent(h_co[key]->FindBin(pT_co)));
    break;
  case WEIGHTSYSTDW:
    if(leptype != 1)return nominal;
    if(verbosity)printf("WEIGHTSYSTDW \t hf = %.4f, lf = %.4f, co = %.4f\n",hf_temp_frac_min,lf_temp_frac_min,co_temp_frac_min);

    return getWeightedFake(hf_temp_frac_nom, lf_temp_frac_nom, co_temp_frac_nom, hf_temp_frac_min,  lf_temp_frac_min, co_temp_frac_min, h_hf, h_lf, h_co, pT_hf, pT_lf, pT_co, key,uncertainty_keys.at(WEIGHTSYSTDW));
    
    
    //return (hf_temp_frac_min*h_hf[key]->GetBinContent(h_hf[key]->FindBin(pT_hf)) +
    //	    lf_temp_frac_min*h_lf[key]->GetBinContent(h_lf[key]->FindBin(pT_lf)) +
    //	    co_temp_frac_min*h_co[key]->GetBinContent(h_co[key]->FindBin(pT_co)));
    break;
    // Not that XSECDW changes the fake rate such that FNP estimate is expected to increase, that's why the names are changed 
  case XSECDW:

    if(verbosity)printf("XSECUP \t f_hf = %.4f, f_lf = %.4f, f_co = %.4f\n",h_hf_xsecup[key]->GetBinContent(h_hf_xsecup[key]->FindBin(pT_hf)),h_lf[key]->GetBinContent(h_lf[key]->FindBin(pT_lf)),h_co[key]->GetBinContent(h_co[key]->FindBin(pT_co)));

    xsecup_temp = (hf_temp_frac_nom*h_hf_xsecup[key]->GetBinContent(h_hf_xsecup[key]->FindBin(pT_hf)) +
    		  lf_temp_frac_nom*h_lf[key]->GetBinContent(h_lf[key]->FindBin(pT_lf)) +
    		  co_temp_frac_nom*h_co[key]->GetBinContent(h_co[key]->FindBin(pT_co)));
    xsecdw_temp = (hf_temp_frac_nom*h_hf_xsecdw[key]->GetBinContent(h_hf_xsecdw[key]->FindBin(pT_hf)) +
		   lf_temp_frac_nom*h_lf[key]->GetBinContent(h_lf[key]->FindBin(pT_lf)) +
		   co_temp_frac_nom*h_co[key]->GetBinContent(h_co[key]->FindBin(pT_co)));
    
    return min(xsecup_temp,xsecdw_temp);
    
    //return (hf_temp_frac_nom*h_hf_xsecup[key]->GetBinContent(h_hf_xsecup[key]->FindBin(pT_hf)) +
    //    lf_temp_frac_nom*h_lf[key]->GetBinContent(h_lf[key]->FindBin(pT_lf)) +
    //	    co_temp_frac_nom*h_co[key]->GetBinContent(h_co[key]->FindBin(pT_co)));
    break;
    // Not that XSECUP changes the fake rate such that FNP estimate is expected to decrease, that's why the names are changed 
  case XSECUP:

    if(verbosity)printf("XSECDW \t f_hf = %.4f, f_lf = %.4f, f_co = %.4f\n",h_hf_xsecdw[key]->GetBinContent(h_hf_xsecdw[key]->FindBin(pT_hf)),h_lf[key]->GetBinContent(h_lf[key]->FindBin(pT_lf)),h_co[key]->GetBinContent(h_co[key]->FindBin(pT_co)));
    
    xsecup_temp = (hf_temp_frac_nom*h_hf_xsecup[key]->GetBinContent(h_hf_xsecup[key]->FindBin(pT_hf)) +
		   lf_temp_frac_nom*h_lf[key]->GetBinContent(h_lf[key]->FindBin(pT_lf)) +
		   co_temp_frac_nom*h_co[key]->GetBinContent(h_co[key]->FindBin(pT_co)));
    xsecdw_temp = (hf_temp_frac_nom*h_hf_xsecdw[key]->GetBinContent(h_hf_xsecdw[key]->FindBin(pT_hf)) +
		   lf_temp_frac_nom*h_lf[key]->GetBinContent(h_lf[key]->FindBin(pT_lf)) +
		   co_temp_frac_nom*h_co[key]->GetBinContent(h_co[key]->FindBin(pT_co)));
    
    return max(xsecup_temp,xsecdw_temp);
    
    //return (hf_temp_frac_nom*h_hf_xsecdw[key]->GetBinContent(h_hf_xsecdw[key]->FindBin(pT_hf)) +
    //	    lf_temp_frac_nom*h_lf[key]->GetBinContent(h_lf[key]->FindBin(pT_lf)) +
    //	    co_temp_frac_nom*h_co[key]->GetBinContent(h_co[key]->FindBin(pT_co)));
    break;
  case SYSTUP:
    if(leptype != 1)return nominal;
    systematic = leptype == 1 ? "EL_"+uncertainty_keys.at(SYSTUP) : "MUON_"+uncertainty_keys.at(SYSTUP);
    //systematic += "_1up";
    // If not available return nominal
    if(h_lf.find(key+"_"+systematic) == h_lf.end()){
      printf("ERROR \t FAKE: Could not find histogram for systematics %s_%s\n",key.Data(),systematic.Data());
      systw = 1.0;
    }else{
      systw = h_lf[key+"_"+systematic]->GetBinContent(h_lf[key+"_"+systematic]->FindBin(pT_lf));
    }
    return (hf_temp_frac_nom*h_hf[key]->GetBinContent(h_hf[key]->FindBin(pT_hf)) +
	    systw*lf_temp_frac_nom*h_lf[key]->GetBinContent(h_lf[key]->FindBin(pT_lf)) +
	    co_temp_frac_nom*h_co[key]->GetBinContent(h_co[key]->FindBin(pT_co)));
    break;
  case SYSTDW:
    if(leptype != 1)return nominal;
    systematic = leptype == 1 ? "EL_"+uncertainty_keys.at(SYSTDW) : "MUON_"+uncertainty_keys.at(SYSTDW);
    //systematic += "_1down";
    // If not available return nominal
    if(h_lf.find(key+"_"+systematic) == h_lf.end()){
      printf("ERROR \t FAKE: Could not find histogram for systematics %s_%s\n",key.Data(),systematic.Data());
      systw = 1.0;
    }else{
      systw = h_lf[key+"_"+systematic]->GetBinContent(h_lf[key+"_"+systematic]->FindBin(pT_lf));
    }
    return (hf_temp_frac_nom*h_hf[key]->GetBinContent(h_hf[key]->FindBin(pT_hf)) +
	    systw*lf_temp_frac_nom*h_lf[key]->GetBinContent(h_lf[key]->FindBin(pT_lf)) +
	    co_temp_frac_nom*h_co[key]->GetBinContent(h_co[key]->FindBin(pT_co)));
    break;
  case TRIGUP:
    if(leptype != 1)return nominal;
    systematic = leptype == 1 ? "EL_"+uncertainty_keys.at(TRIGUP) : "MUON_"+uncertainty_keys.at(TRIGUP);
    if(leptype == 1)systematic = systematic.ReplaceAll("TRIG","Trigger");
    //systematic += "_1up";
    // If not available return nominal
    if(h_lf.find(key+"_"+systematic) == h_lf.end()){
      printf("ERROR \t FAKE: Could not find histogram for systematics %s_%s\n",key.Data(),systematic.Data());
      systw = 1.0;
    }else{
      systw = h_lf[key+"_"+systematic]->GetBinContent(h_lf[key+"_"+systematic]->FindBin(pT_lf));
    }
    return (hf_temp_frac_nom*h_hf[key]->GetBinContent(h_hf[key]->FindBin(pT_hf)) +
	    systw*lf_temp_frac_nom*h_lf[key]->GetBinContent(h_lf[key]->FindBin(pT_lf)) +
	    co_temp_frac_nom*h_co[key]->GetBinContent(h_co[key]->FindBin(pT_co)));
    break;
  case TRIGDW:
    if(leptype != 1)return nominal;
    systematic = leptype == 1 ? "EL_"+uncertainty_keys.at(TRIGDW) : "MUON_"+uncertainty_keys.at(TRIGDW);
    if(leptype == 1)systematic = systematic.ReplaceAll("TRIG","Trigger");
    //systematic += "_1down";
    // If not available return nominal
    if(h_lf.find(key+"_"+systematic) == h_lf.end()){
      printf("ERROR \t FAKE: Could not find histogram for systematics %s_%s\n",key.Data(),systematic.Data());
      systw = 1.0;
    }else{
      systw = h_lf[key+"_"+systematic]->GetBinContent(h_lf[key+"_"+systematic]->FindBin(pT_lf));
    }
    return (hf_temp_frac_nom*h_hf[key]->GetBinContent(h_hf[key]->FindBin(pT_hf)) +
	    systw*lf_temp_frac_nom*h_lf[key]->GetBinContent(h_lf[key]->FindBin(pT_lf)) +
	    co_temp_frac_nom*h_co[key]->GetBinContent(h_co[key]->FindBin(pT_co)));
    break;
  default:
    return (hf_temp_frac_nom*h_hf[key]->GetBinContent(h_hf[key]->FindBin(pT_hf)) +
	    lf_temp_frac_nom*h_lf[key]->GetBinContent(h_lf[key]->FindBin(pT_lf)) +
	    co_temp_frac_nom*h_co[key]->GetBinContent(h_co[key]->FindBin(pT_co)));
  }
  
}

double GetFakeWeight::getReal(float pT, float eta, int leptype, TString triggermatch, uncertainty u, TString BDTname, double BDTvalue){


  std::pair <Double_t,Double_t> upcut;
  Double_t pT_real = pT;
  Double_t eta_real = eta;
  TString key;
  // If no triggermatch, set the trigger to the ana key
  if(triggermatch.Length() <= 0)triggermatch = this_ana;

  key = this_year+"_"+triggermatch;
  
  
  if(verbosity){
    cout<<"-----------------------------------------------------------------------------------------"<<endl; 
    cout<<"REAL: Getting "<<key<< " for "<< (leptype == 1 ? "el" : "mu") << " with pT "<< pT<< ", eta "<<eta<<endl;
  }

  eta = fabs(eta);
  
  switch(leptype){
  case 1:

    if(BDTvalue < 0){
      upcut = getUpperCut("el_real_"+key);
      pT_real = pT < upcut.first ? pT : upcut.first;
      eta_real = eta < upcut.second ? eta : upcut.second;
    }
    
    if(verbosity)printf("Electron: pT = %.2f (%.2f), eta = %.2f (%.2f), r = %.4f\n",pT,pT_real,eta,eta_real,getRealUnc(key, leptype, pT_real, eta_real, h_el_real, u, h_el_real_1D, BDTname, BDTvalue));
    return getRealUnc(key, leptype, pT_real, eta_real, h_el_real, u, h_el_real_1D, BDTname, BDTvalue);
    //return h_el_real[key]->GetBinContent(h_el_real[key]->FindBin(pT_real,eta_real));
    break;
  case 2:
    
    if(BDTvalue < 0){ 
      upcut = getUpperCut("mu_real_"+key);
      pT_real = pT < upcut.first ? pT : upcut.first;
      eta_real = eta < upcut.second ? eta : upcut.second;
    }
    
    if(verbosity)printf("Muon   : pT = %.2f (%.2f), eta = %.2f (%.2f), r = %.4f\n",pT,pT_real,eta,eta_real,getRealUnc(key, leptype, pT_real, eta_real, h_mu_real, u, h_mu_real_1D, BDTname, BDTvalue));
    return getRealUnc(key, leptype, pT_real, eta_real, h_mu_real, u, h_mu_real_1D, BDTname, BDTvalue);
    //return h_mu_real[key]->GetBinContent(h_mu_real[key]->FindBin(pT_real,eta_real));
  default:
    cout<<"Could not find lepton type "<<leptype<<endl;
    return -1;
  }
}

std::pair <Double_t,Double_t> GetFakeWeight::getUpperCut(TString key){
  std::pair <Double_t,Double_t> ret (999,999);
  if((upper_val_pt.find(key) == upper_val_pt.end()) && (upper_val_eta.find(key) == upper_val_eta.end())){
    printf("ERROR \t Could not find upper pt nor eta cut for %s\n",key.Data());
    return ret;
  }
  if(!(upper_val_pt.find(key) == upper_val_pt.end())){
    ret.first = upper_val_pt[key]-0.5;
  }else{
    printf("ERROR \t Could not find upper pt cut for %s\n",key.Data());
  }
  if(!(upper_val_eta.find(key) == upper_val_eta.end())){
    ret.second = upper_val_eta[key]-0.01;
  }else if(key.Contains("_real_")){
    printf("ERROR \t Could not find upper eta cut for %s\n",key.Data());
  } 
  return ret; 
}



