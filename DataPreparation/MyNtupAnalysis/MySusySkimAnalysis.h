//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue May  8 09:28:22 2018 by ROOT version 6.10/04
// from TChain tree_NoSys/
//////////////////////////////////////////////////////////

#ifndef MySusySkimAnalysis_h
#define MySusySkimAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include "MatrixMethod.h"
#include "GetFakeWeight.h"
#include <map>
#include <TH1.h>
#include <TH2.h>
#include "TLorentzVector.h"
#include "MatrixMethod.cxx"
#include "HistFitterTree.h"
#include "make2L2JTree.h"
#include "makeMYTree.h"
#include "TRandom.h"

// Headers needed by this particular selector
#include <vector>

//enum truthENUM { CF = 0, HF = 1, LF = 2, CO = 3, REAL = 4, UNKNOWN = 5, DATA = 6};
//std::vector<TString> truthVEC = {"CF", "HF", "LF", "CO", "REAL", "UNKNOWN", "DATA"}; 


enum truthENUMIFF {
  Unknown = 0,
  KnownUnknown = 1,
  IsoElectron = 2,
  ChargeFlipIsoElectron = 3,
  PromptMuon = 4,
  PromptPhotonConversion = 5,
  ElectronFromMuon = 6,
  TauDecay = 7,
  BHadronDecay = 8,
  CHadronDecay = 9,
  LightFlavorDecay = 10,
  DATA = 11,
};

std::vector<TString> truthVECIFF = {
  "Unknown",
  "KnownUnknown",
  "IsoElectron",
  "ChargeFlipIsoElectron",
  "PromptMuon",
  "PromptPhotonConversion",
  "ElectronFromMuon",
  "TauDecay",
  "BHadronDecay",
  "CHadronDecay",
  "LightFlavorDecay",
  "DATA",
};

typedef map< TString, TH1F* > th1Map;
typedef map< TString, TH2F* > th2Map;

class MySusySkimAnalysis : public TSelector {
public :
  TTreeReader     fReader;  //!the tree reader
  TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain
  TSelectorList  *fOutput;

  TTree *tree3;
  TTree *exactcopy;
  TTree *newtree;
  TFile *ntupf;

  TChain *chain;

  Double_t extrafac;

  std::map<TString,Double_t> b_FNPweight;
  
  TFile *ttree_infile;
  TFile *ttree_newfile;
  TTree *centraltree; 
  Long64_t central_nentries;
  Double_t MMweight;
  Double_t MMweight_up;
  Double_t MMweight_dw;
  TBranch *bMMweight;    
  TBranch *bMMweight_up; 
  TBranch *bMMweight_dw;

  Double_t FNPest;
  Double_t FNPest_up;
  Double_t FNPest_dw;

  std::map<TString,Float_t> el_trg_th;
  std::map<TString,Float_t> mu_trg_th;
  
  std::map<Int_t,std::vector<Bool_t>> lep1isSignal;
  std::map<Int_t,std::vector<Bool_t>> lep2isSignal;
  std::map<Int_t,std::vector<Float_t>> lepton1pT;
  std::map<Int_t,std::vector<Float_t>> lepton2pT;

  std::map<Int_t,std::vector<Int_t>> glob_njet;
  std::map<Int_t,std::vector<Int_t>> glob_nbjet;
  std::map<Int_t,std::vector<Bool_t>> glob_isSF;  
  std::map<Int_t,std::vector<Bool_t>> glob_isOS;  
  std::map<Int_t,std::vector<Float_t>> glob_METsig;
  
  std::map<Int_t,std::vector<Double_t>> glob_FNP_TOTAL_UP;
  std::map<Int_t,std::vector<Double_t>> glob_FNP_TOTAL_DOWN;
  std::map<Int_t,std::vector<Double_t>> glob_FNP_WEIGHTS;

  std::map<Int_t,std::map<TString,std::vector<double>>> BDTscores;


  int fnpidx;
  std::map<TString,std::vector<double>> BDTscoremap;
  std::map<Int_t,std::vector<Long64_t>> evnums;
  std::vector<std::string> bdt_vec;
  std::map<TString,Double_t> BDTweight;
  std::map<Int_t,std::vector<TString>> evtype;
  vector<int> done_fnpidx;
  vector<int> years;

  int ev_was_found;
 
  Int_t yr;
  
  Long64_t toprocess;
  Long64_t firstentry; 

  Int_t nlep_base_OR;

  std::map<TString, std::vector<TTreeReaderValue<Float_t>> > uncert;
  vector<TString> uncert_keys;
  
  vector<TString> triglist;
  vector< TTreeReaderValue<Bool_t> > trigval;
  vector< TTreeReaderValue<vector<bool>> > trigvallep;
  std::map<int,std::vector< TTreeReaderValue<bool>>> triglist_e;
  std::map<int,std::vector< TTreeReaderValue<bool>>> triglist_m;
  std::map<int,std::vector<TTreeReaderValue<vector<bool>>>> trigmatch_e;
  std::map<int,std::vector<TTreeReaderValue<vector<bool>>>> trigmatch_m;
  std::map<int, std::vector<TString> > trigstr_e;
  std::map<int, std::vector<TString> > trigstr_m;

  vector<ULong64_t> DuplicationCheck;
  vector<Int_t> DSIDcheck;

  int nDup;
  bool doDupCheck;

  TH1F* h_triggermatched_el;
  TH1F* h_triggermatched_mu;


  TString cutflowstr;
  std::vector<TString> cutfloworder;
  map< TString, float > cutflow;
  map< TString, float > cutflow_wgt;
  map< TString,TString > trigcat;
  std::vector<TString> trig_order;
  map< TString,TString > trigcategories;
  std::vector<TString> categories;
  std::vector<TString> triggerDef;

  std::vector<ULong64_t> mirto_ev;

  HistFitterTree *HFT;
  make2L2JTree *twoLtwoJ;
  makeMYTree *MY;

  TString LT_postfix;
  
  bool makeHistograms;
  //bool computeRates;
  bool fullScan;

  ofstream myfile;
  int nmirto;

  TRandom *rnd;

  bool doFakeEst;
  bool makeNtup;
  bool make2L2JNtup;
  bool makeMYNtup;
  bool doRJR;
  bool allEventsPass;
  bool doCoreRegions;
  bool doAnalysis;

  bool is1516;
  bool is17;
  bool is18;

  Double_t glob_r1;
  Double_t glob_r2;
  Double_t glob_f1;
  Double_t glob_f2;

  std::map<TString,Double_t> glob_r_wgt;
  std::map<TString,Double_t> glob_rf_wgt;
  std::map<TString,Double_t> glob_fr_wgt;
  std::map<TString,Double_t> glob_f_wgt ;
  std::map<TString,Double_t> b_fake_wgt ;
  

  TString date_folder;
  int num_of_leptons;
  TString fname_old;
  TString fname;
  char fnameoutputntupfilename[200];
  char fnameoutputhistfilename[200];
  char fnameoutputcopyfilename[200];
  TString rnum;
  TString fullpath;
  TLorentzVector met_tlv;
  TLorentzVector met_tlv2;
  Double_t metrel_Et;
  std::vector<TString> fillregions;
  std::vector<TString> cutsregions;
  std::vector<TString> fillregions_wcuts;
  std::vector<TString> looseTightDef;

  TH1F *h_min_dphi_lep_met;
  
  TH1F *h_fake_weights;
  TH1F *h_runNumber18;

  TH1F *h_passed_trigMatch_2LTrig_fakeVR;
  
  TH1F* h_probe_el_hf_deltaR_ORjet;
  TH1F* h_probe_el_hf_deltaR_ORbjet;
  TH1F* h_probe_mu_hf_deltaR_ORjet;
  TH1F* h_probe_mu_hf_deltaR_ORbjet;

  TH1F *h_ptllboost;

  TH1F *h_bornMass;

  double themt2;
  double mtw;
  int n_clj;
  int n_clj30;
  int n_bjet;
  int n_bjet77;
  int n_bjet85;
  int n_sgj;
  double themll;

  th1Map h_evt_BDT_TT;
  th1Map h_evt_BDT_Tl;
  th1Map h_evt_BDT_lT;
  th1Map h_evt_BDT_ll;

  th1Map h_real1_diff;
  th1Map h_real2_diff;

  th2Map h_syst_MM;

  th1Map h_lep_nbljetpreOR;

  th1Map h_lep_njet30_nT;
  th1Map h_lep_njet30_nL;

  th1Map h_lep_pT0_nT;
  th1Map h_lep_pT0_nL;

  th1Map h_lep_MET_nL;
  th1Map h_lep_MET_nT;

  th1Map h_lep_nbjet_nL;
  th1Map h_lep_nbjet_nT;

  th1Map h_lep_pT_nT;
  th1Map h_lep_pT_nL;
  th1Map h_lep_BDT_nT;
  th1Map h_lep_BDT_nL;
  th1Map h_lep_eta_nT;
  th1Map h_lep_eta_nL;
  th1Map h_lep_TCF_nL;
  th1Map h_lep_TCF_nT;
  th1Map h_lep_trig_nT;
  th1Map h_lep_trig_nL;

  th1Map h_lep_pT0;
  th1Map h_lep_pT1;
  th1Map h_lep_pT2;
  th1Map h_lep_mt2;
  th1Map h_lep_mtw;
  th1Map h_lep_met;
  th1Map h_lep_mll;
  th1Map h_lep_njet30;
  th1Map h_lep_nbj77;
  th1Map h_lep_nbj85;

  th2Map h_el_fake_hf;
  th2Map h_el_fake_hf_xsecup;
  th2Map h_el_fake_hf_xsecdw;
  th1Map h_el_fake_co;
  th2Map h_el_real;
  th2Map h_mu_real;
  th2Map h_el_real_mc;
  th2Map h_mu_real_mc;

  th1Map h_mu_fake_hf;
  th1Map h_mu_fake_hf_xsecup;
  th1Map h_mu_fake_hf_xsecdw;
  th1Map h_el_fake_lf;

  th1Map h_lep_metsig;
  /* th1Map h_el_fake_lf_xsecup; */
  /* th1Map h_el_fake_lf_xsecdw; */

  th2Map h_lep_pT_eta_nT;
  th2Map h_lep_pT_eta_nL;
  th2Map h_lep_mT2_TCF_nT;
  th2Map h_lep_mT2_TCF_nL;
  th2Map h_lep_mTW_TCF_nT;
  th2Map h_lep_mTW_TCF_nL;
  th2Map h_lep_njet_TCF_nT;
  th2Map h_lep_njet_TCF_nL;

  th2Map h_lep_metsig_TCF_nL;
  th2Map h_lep_metsig_TCF_nT;

  th2Map h_lep_pT_TCF_nL;
  th2Map h_lep_pT_TCF_nT;

  th2Map h_lep_type_origin_nT;
  th2Map h_lep_type_origin_nL;
  th2Map h_lep_pT_origin_nL;
  th2Map h_lep_pT_origin_nT;
  th2Map h_lep_pT_type_nL;
  th2Map h_lep_pT_type_nT;
  th2Map h_lep_pT_firstEgMotherPdgId_nL;
  th2Map h_lep_pT_firstEgMotherPdgId_nT;
  th2Map h_lep_pT_firstEgMotherO_nL;
  th2Map h_lep_pT_firstEgMotherO_nT;
  th2Map h_lep_pT_firstEgMotherT_nL;
  th2Map h_lep_pT_firstEgMotherT_nT;

  th2Map h_lep_mu_pT_nL;
  th2Map h_lep_mu_pT_nT;
  
  th1Map h_el_frac_hf;
  th1Map h_el_frac_lf;
  th1Map h_el_frac_co;
  th1Map h_el_frac_cf;

  th1Map h_MMevtype;

  TH1F *h_HFreg_mindr_jet_probe;
  TH1F *h_HFreg_nprobes;
  TH2F *h_HFreg_ntags_nprobes;
  TH1F *h_lep_lepECIDS;
   
  Long64_t nentries;
   
  int nev;
  int all_nev;

  MatrixMethod *mtx;
  GetFakeWeight *gfw;

  double wgt;
  double wgt_loose;
  bool isData;


  TH1F *h_wgts_leptonWeight;
  TH1F *h_wgts_baselineleptonWeight;
  TH1F *h_wgts_globalDiLepTrigSF;
  TH1F *h_wgts_singleLepTrigSF;
  TH1F *h_wgts_globalBaselineDiLepTrigSF;
  TH1F *h_wgts_singleBaselineLepTrigSF;

  TH1F *h_wgts_lepBLRecoSF_lep1;
  TH1F *h_wgts_lepBLRecoSF_lep2;
  TH1F *h_wgts_lepRecoSF_lep1;
  TH1F *h_wgts_lepRecoSF_lep2;

  TH1F *h_wgts_lepTrigSF_lep1;
  TH1F *h_wgts_lepTrigSF_lep2;
  TH1F *h_wgts_lepBLTrigSF_lep1;
  TH1F *h_wgts_lepBLTrigSF_lep2;

  TH1F *h_wgts_wgt;
  TH1F *h_wgts_wgtloose;

  TH1F *h_wgts_infilllandt_wgt;
  TH1F *h_wgts_infilllandt_wgtloose;

  
  // REL22 ntuples
  // Readers to access the data (delete the ones you do not need).
  /*
    TTreeReaderValue<vector<bool>> lepHLT_e24_lhmedium_L1EM20VH = {fReader, "lepHLT_e24_lhmedium_L1EM20VH"};
    TTreeReaderValue<vector<bool>> lepHLT_e60_lhmedium = {fReader, "lepHLT_e60_lhmedium"};
    TTreeReaderValue<vector<bool>> lepHLT_e120_lhloose = {fReader, "lepHLT_e120_lhloose"};
    TTreeReaderValue<vector<bool>> lepHLT_mu20_iloose_L1MU15 = {fReader, "lepHLT_mu20_iloose_L1MU15"};
    TTreeReaderValue<vector<bool>> lepHLT_2e12_lhloose_L12EM10VH = {fReader, "lepHLT_2e12_lhloose_L12EM10VH"};
    TTreeReaderValue<vector<bool>> lepHLT_mu18_mu8noL1 = {fReader, "lepHLT_mu18_mu8noL1"};
    TTreeReaderValue<vector<bool>> lepHLT_e17_lhloose_mu14 = {fReader, "lepHLT_e17_lhloose_mu14"};
    TTreeReaderValue<vector<bool>> lepHLT_e7_lhmedium_mu24 = {fReader, "lepHLT_e7_lhmedium_mu24"};
    TTreeReaderValue<vector<bool>> lepHLT_e24_lhtight_nod0_ivarloose = {fReader, "lepHLT_e24_lhtight_nod0_ivarloose"};
    TTreeReaderValue<vector<bool>> lepHLT_e24_lhmedium_nod0_L1EM20VH = {fReader, "lepHLT_e24_lhmedium_nod0_L1EM20VH"};
    TTreeReaderValue<vector<bool>> lepHLT_e60_medium = {fReader, "lepHLT_e60_medium"};
    TTreeReaderValue<vector<bool>> lepHLT_mu40 = {fReader, "lepHLT_mu40"};
    TTreeReaderValue<vector<bool>> lepHLT_mu24_iloose_L1MU15 = {fReader, "lepHLT_mu24_iloose_L1MU15"};
    TTreeReaderValue<vector<bool>> lepHLT_mu24_ivarloose_L1MU15 = {fReader, "lepHLT_mu24_ivarloose_L1MU15"};
    TTreeReaderValue<vector<bool>> lepHLT_mu24_ivarmedium = {fReader, "lepHLT_mu24_ivarmedium"};
    TTreeReaderValue<vector<bool>> lepHLT_mu24_imedium = {fReader, "lepHLT_mu24_imedium"};
    TTreeReaderValue<vector<bool>> lepHLT_mu26_imedium = {fReader, "lepHLT_mu26_imedium"};
    TTreeReaderValue<vector<bool>> lepHLT_2e15_lhvloose_nod0_L12EM13VH = {fReader, "lepHLT_2e15_lhvloose_nod0_L12EM13VH"};
    TTreeReaderValue<vector<bool>> lepHLT_2e17_lhvloose_nod0 = {fReader, "lepHLT_2e17_lhvloose_nod0"};
    TTreeReaderValue<vector<bool>> lepHLT_2mu10 = {fReader, "lepHLT_2mu10"};
    TTreeReaderValue<vector<bool>> lepHLT_2mu14 = {fReader, "lepHLT_2mu14"};
    TTreeReaderValue<vector<bool>> lepHLT_mu20_mu8noL1 = {fReader, "lepHLT_mu20_mu8noL1"};
    TTreeReaderValue<vector<bool>> lepHLT_mu22_mu8noL1 = {fReader, "lepHLT_mu22_mu8noL1"};
    TTreeReaderValue<vector<bool>> lepHLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1 = {fReader, "lepHLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1"};
    TTreeReaderValue<vector<bool>> lepHLT_e26_lhtight_nod0_ivarloose = {fReader, "lepHLT_e26_lhtight_nod0_ivarloose"};
    TTreeReaderValue<vector<bool>> lepHLT_e26_lhtight_nod0 = {fReader, "lepHLT_e26_lhtight_nod0"};
    TTreeReaderValue<vector<bool>> lepHLT_e60_lhmedium_nod0 = {fReader, "lepHLT_e60_lhmedium_nod0"};
    TTreeReaderValue<vector<bool>> lepHLT_e140_lhloose_nod0 = {fReader, "lepHLT_e140_lhloose_nod0"};
    TTreeReaderValue<vector<bool>> lepHLT_e300_etcut = {fReader, "lepHLT_e300_etcut"};
    TTreeReaderValue<vector<bool>> lepHLT_mu26_ivarmedium = {fReader, "lepHLT_mu26_ivarmedium"};
    TTreeReaderValue<vector<bool>> lepHLT_mu50 = {fReader, "lepHLT_mu50"};
    TTreeReaderValue<vector<bool>> lepHLT_mu60_0eta105_msonly = {fReader, "lepHLT_mu60_0eta105_msonly"};
    TTreeReaderValue<vector<bool>> lepHLT_2e17_lhvloose_nod0_L12EM15VHI = {fReader, "lepHLT_2e17_lhvloose_nod0_L12EM15VHI"};
    TTreeReaderValue<vector<bool>> lepHLT_2e24_lhvloose_nod0 = {fReader, "lepHLT_2e24_lhvloose_nod0"};
    TTreeReaderValue<vector<bool>> lepHLT_e17_lhloose_nod0_mu14 = {fReader, "lepHLT_e17_lhloose_nod0_mu14"};
    TTreeReaderValue<vector<bool>> lepHLT_e26_lhmedium_nod0_mu8noL1 = {fReader, "lepHLT_e26_lhmedium_nod0_mu8noL1"};
    TTreeReaderValue<vector<bool>> lepHLT_e7_lhmedium_nod0_mu24 = {fReader, "lepHLT_e7_lhmedium_nod0_mu24"};
    TTreeReaderValue<Bool_t> trigMatch_1L2LTrig = {fReader, "trigMatch_1L2LTrig"};
    TTreeReaderValue<Bool_t> trigMatch_1LTrig = {fReader, "trigMatch_1LTrig"};
    TTreeReaderValue<Bool_t> trigMatch_1L2LTrigOR = {fReader, "trigMatch_1L2LTrigOR"};
    TTreeReaderValue<Bool_t> trigMatch_1LTrigOR = {fReader, "trigMatch_1LTrigOR"};
    TTreeReaderValue<Bool_t> trigMatch_2LTrig = {fReader, "trigMatch_2LTrig"};
    TTreeReaderValue<Bool_t> trigMatch_2LTrigOR = {fReader, "trigMatch_2LTrigOR"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_e24_lhmedium_L1EM20VH = {fReader, "trigMatch_HLT_e24_lhmedium_L1EM20VH"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_e60_lhmedium = {fReader, "trigMatch_HLT_e60_lhmedium"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_e120_lhloose = {fReader, "trigMatch_HLT_e120_lhloose"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_mu20_iloose_L1MU15 = {fReader, "trigMatch_HLT_mu20_iloose_L1MU15"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_2e12_lhloose_L12EM10VH = {fReader, "trigMatch_HLT_2e12_lhloose_L12EM10VH"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_mu18_mu8noL1 = {fReader, "trigMatch_HLT_mu18_mu8noL1"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_e17_lhloose_mu14 = {fReader, "trigMatch_HLT_e17_lhloose_mu14"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_e7_lhmedium_mu24 = {fReader, "trigMatch_HLT_e7_lhmedium_mu24"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_e24_lhtight_nod0_ivarloose = {fReader, "trigMatch_HLT_e24_lhtight_nod0_ivarloose"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_e24_lhmedium_nod0_L1EM20VH = {fReader, "trigMatch_HLT_e24_lhmedium_nod0_L1EM20VH"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_e60_medium = {fReader, "trigMatch_HLT_e60_medium"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_mu40 = {fReader, "trigMatch_HLT_mu40"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_mu24_iloose_L1MU15 = {fReader, "trigMatch_HLT_mu24_iloose_L1MU15"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_mu24_ivarloose_L1MU15 = {fReader, "trigMatch_HLT_mu24_ivarloose_L1MU15"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_mu24_ivarmedium = {fReader, "trigMatch_HLT_mu24_ivarmedium"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_mu24_imedium = {fReader, "trigMatch_HLT_mu24_imedium"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_mu26_imedium = {fReader, "trigMatch_HLT_mu26_imedium"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_2e15_lhvloose_nod0_L12EM13VH = {fReader, "trigMatch_HLT_2e15_lhvloose_nod0_L12EM13VH"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_2e17_lhvloose_nod0 = {fReader, "trigMatch_HLT_2e17_lhvloose_nod0"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_2mu10 = {fReader, "trigMatch_HLT_2mu10"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_2mu14 = {fReader, "trigMatch_HLT_2mu14"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_mu20_mu8noL1 = {fReader, "trigMatch_HLT_mu20_mu8noL1"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_mu22_mu8noL1 = {fReader, "trigMatch_HLT_mu22_mu8noL1"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1 = {fReader, "trigMatch_HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_e26_lhtight_nod0_ivarloose = {fReader, "trigMatch_HLT_e26_lhtight_nod0_ivarloose"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_e26_lhtight_nod0 = {fReader, "trigMatch_HLT_e26_lhtight_nod0"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_e60_lhmedium_nod0 = {fReader, "trigMatch_HLT_e60_lhmedium_nod0"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_e140_lhloose_nod0 = {fReader, "trigMatch_HLT_e140_lhloose_nod0"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_e300_etcut = {fReader, "trigMatch_HLT_e300_etcut"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_mu26_ivarmedium = {fReader, "trigMatch_HLT_mu26_ivarmedium"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_mu50 = {fReader, "trigMatch_HLT_mu50"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_mu60_0eta105_msonly = {fReader, "trigMatch_HLT_mu60_0eta105_msonly"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_2e17_lhvloose_nod0_L12EM15VHI = {fReader, "trigMatch_HLT_2e17_lhvloose_nod0_L12EM15VHI"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_2e24_lhvloose_nod0 = {fReader, "trigMatch_HLT_2e24_lhvloose_nod0"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_e17_lhloose_nod0_mu14 = {fReader, "trigMatch_HLT_e17_lhloose_nod0_mu14"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_e26_lhmedium_nod0_mu8noL1 = {fReader, "trigMatch_HLT_e26_lhmedium_nod0_mu8noL1"};
    TTreeReaderValue<Bool_t> trigMatch_HLT_e7_lhmedium_nod0_mu24 = {fReader, "trigMatch_HLT_e7_lhmedium_nod0_mu24"};
    TTreeReaderValue<Float_t> mu = {fReader, "mu"};
    TTreeReaderValue<Float_t> avg_mu = {fReader, "avg_mu"};
    TTreeReaderValue<Float_t> actual_mu = {fReader, "actual_mu"};
    TTreeReaderValue<Int_t> nVtx = {fReader, "nVtx"};
    TTreeReaderValue<Int_t> channel = {fReader, "channel"};
    TTreeReaderValue<Int_t> nLep_base = {fReader, "nLep_base"};
    TTreeReaderValue<Int_t> nLep_signal = {fReader, "nLep_signal"};
    TTreeReaderArray<int> lepFlavor = {fReader, "lepFlavor"};
    TTreeReaderArray<int> lepCharge = {fReader, "lepCharge"};
    TTreeReaderArray<int> lepAuthor = {fReader, "lepAuthor"};
    TTreeReaderArray<float> lepPt = {fReader, "lepPt"};
    TTreeReaderArray<float> lepEta = {fReader, "lepEta"};
    TTreeReaderArray<float> lepPhi = {fReader, "lepPhi"};
    TTreeReaderArray<float> lepM = {fReader, "lepM"};
    TTreeReaderArray<float> lepD0 = {fReader, "lepD0"};
    TTreeReaderArray<float> lepD0Sig = {fReader, "lepD0Sig"};
    TTreeReaderArray<float> lepZ0 = {fReader, "lepZ0"};
    TTreeReaderArray<float> lepZ0SinTheta = {fReader, "lepZ0SinTheta"};
    TTreeReaderValue<vector<bool>> lepPassOR = {fReader, "lepPassOR"};
    TTreeReaderArray<int> lepType = {fReader, "lepType"};
    TTreeReaderArray<int> lepOrigin = {fReader, "lepOrigin"};
    TTreeReaderArray<int> lepIFFClass = {fReader, "lepIFFClass"};
    TTreeReaderArray<int> lepEgMotherType = {fReader, "lepEgMotherType"};
    TTreeReaderArray<int> lepEgMotherOrigin = {fReader, "lepEgMotherOrigin"};
    TTreeReaderArray<int> lepEgMotherPdgId = {fReader, "lepEgMotherPdgId"};
    TTreeReaderArray<int> lepECIDS = {fReader, "lepECIDS"};
    TTreeReaderValue<vector<bool>> lepPassBL = {fReader, "lepPassBL"};
    TTreeReaderValue<vector<bool>> lepVeryLoose = {fReader, "lepVeryLoose"};
    TTreeReaderValue<vector<bool>> lepLoose = {fReader, "lepLoose"};
    TTreeReaderValue<vector<bool>> lepMedium = {fReader, "lepMedium"};
    TTreeReaderValue<vector<bool>> lepTight = {fReader, "lepTight"};
    TTreeReaderValue<vector<bool>> lepIsoFCHighPtCaloOnly = {fReader, "lepIsoFCHighPtCaloOnly"};
    TTreeReaderValue<vector<bool>> lepIsoFCLoose = {fReader, "lepIsoFCLoose"};
    TTreeReaderValue<vector<bool>> lepIsoFCTight = {fReader, "lepIsoFCTight"};
    TTreeReaderValue<vector<bool>> lepIsoFCLoose_FixedRad = {fReader, "lepIsoFCLoose_FixedRad"};
    TTreeReaderValue<vector<bool>> lepIsoFCTight_FixedRad = {fReader, "lepIsoFCTight_FixedRad"};
    TTreeReaderValue<vector<bool>> lepIsoHighPtCaloOnly = {fReader, "lepIsoHighPtCaloOnly"};
    TTreeReaderValue<vector<bool>> lepIsoTightTrackOnly_VarRad = {fReader, "lepIsoTightTrackOnly_VarRad"};
    TTreeReaderValue<vector<bool>> lepIsoTightTrackOnly_FixedRad = {fReader, "lepIsoTightTrackOnly_FixedRad"};
    TTreeReaderValue<vector<bool>> lepIsoLoose_VarRad = {fReader, "lepIsoLoose_VarRad"};
    TTreeReaderValue<vector<bool>> lepIsoTight_VarRad = {fReader, "lepIsoTight_VarRad"};
    TTreeReaderValue<vector<bool>> lepIsoPLVLoose = {fReader, "lepIsoPLVLoose"};
    TTreeReaderValue<vector<bool>> lepIsoPLVTight = {fReader, "lepIsoPLVTight"};
    TTreeReaderValue<vector<bool>> lepTruthMatched = {fReader, "lepTruthMatched"};
    TTreeReaderArray<int> lepTruthCharge = {fReader, "lepTruthCharge"};
    TTreeReaderArray<float> lepTruthPt = {fReader, "lepTruthPt"};
    TTreeReaderArray<float> lepTruthEta = {fReader, "lepTruthEta"};
    TTreeReaderArray<float> lepTruthPhi = {fReader, "lepTruthPhi"};
    TTreeReaderArray<float> lepTruthM = {fReader, "lepTruthM"};
    TTreeReaderArray<float> lepTrigSF = {fReader, "lepTrigSF"};
    TTreeReaderArray<float> lepRecoSF = {fReader, "lepRecoSF"};
    TTreeReaderArray<float> lepBLTrigSF = {fReader, "lepBLTrigSF"};
    TTreeReaderArray<float> lepBLRecoSF = {fReader, "lepBLRecoSF"};
    TTreeReaderValue<Int_t> nJet30 = {fReader, "nJet30"};
    TTreeReaderValue<Int_t> nJet20 = {fReader, "nJet20"};
    TTreeReaderValue<Int_t> nBJet20_MV2c10_FixedCutBEff_77 = {fReader, "nBJet20_MV2c10_FixedCutBEff_77"};
    TTreeReaderArray<float> jetPt = {fReader, "jetPt"};
    TTreeReaderArray<float> jetEta = {fReader, "jetEta"};
    TTreeReaderArray<float> jetPhi = {fReader, "jetPhi"};
    TTreeReaderArray<float> jetM = {fReader, "jetM"};
    TTreeReaderArray<float> jetJVT = {fReader, "jetJVT"};
    TTreeReaderValue<vector<bool>> jetPassOR = {fReader, "jetPassOR"};
    TTreeReaderValue<vector<bool>> jetSignal = {fReader, "jetSignal"};
    TTreeReaderValue<Float_t> mjj = {fReader, "mjj"};
    TTreeReaderArray<float> jetTileEnergy = {fReader, "jetTileEnergy"};
    TTreeReaderArray<float> jetMV2c10 = {fReader, "jetMV2c10"};
    TTreeReaderArray<float> jetdl1 = {fReader, "jetdl1"};
    TTreeReaderValue<Float_t> met_Et = {fReader, "met_Et"};
    TTreeReaderValue<Float_t> met_Sign = {fReader, "met_Sign"};
    TTreeReaderValue<Float_t> met_Phi = {fReader, "met_Phi"};
    TTreeReaderValue<Float_t> mll = {fReader, "mll"};
    TTreeReaderValue<Double_t> pileupWeight = {fReader, "pileupWeight"};
    TTreeReaderValue<Double_t> leptonWeight = {fReader, "leptonWeight"};
    TTreeReaderValue<Double_t> baselineleptonWeight = {fReader, "baselineleptonWeight"};
    TTreeReaderValue<Double_t> eventWeight = {fReader, "eventWeight"};
    TTreeReaderValue<Double_t> bTagWeight = {fReader, "bTagWeight"};
    TTreeReaderValue<Double_t> jvtWeight = {fReader, "jvtWeight"};
    TTreeReaderValue<Double_t> globalDiLepTrigSF = {fReader, "globalDiLepTrigSF"};
    TTreeReaderValue<Double_t> globalBaselineDiLepTrigSF = {fReader, "globalBaselineDiLepTrigSF"};
    TTreeReaderValue<Double_t> flavSymWeight = {fReader, "flavSymWeight"};
    TTreeReaderValue<Double_t> genWeight = {fReader, "genWeight"};
    TTreeReaderValue<Double_t> genWeightUp = {fReader, "genWeightUp"};
    TTreeReaderValue<Double_t> genWeightDown = {fReader, "genWeightDown"};
    TTreeReaderValue<ULong64_t> PRWHash = {fReader, "PRWHash"};
    TTreeReaderValue<ULong64_t> EventNumber = {fReader, "EventNumber"};
    TTreeReaderValue<Float_t> xsec = {fReader, "xsec"};
    TTreeReaderValue<Float_t> GenHt = {fReader, "GenHt"};
    TTreeReaderValue<Float_t> GenMET = {fReader, "GenMET"};
    TTreeReaderValue<Int_t> DatasetNumber = {fReader, "DatasetNumber"};
    TTreeReaderValue<Int_t> RunNumber = {fReader, "RunNumber"};
    TTreeReaderValue<Int_t> RandomRunNumber = {fReader, "RandomRunNumber"};
    TTreeReaderValue<Int_t> FS = {fReader, "FS"};
  */
  /**
     TTreeReaderArray<float> LHE3Weights = {fReader, "LHE3Weights"};
     TTreeReaderArray<TString> LHE3WeightNames = {fReader, "LHE3WeightNames"};
  */
  // Readers to access the data (delete the ones you do not need).

  //My REL21 ANALYSIS STARTS HERE!
  TTreeReaderValue<Float_t> bornMass = {fReader,"mu"};// = {fReader, "bornMass"};
  
  TTreeReaderValue<vector<bool>> lepHLT_e24_lhvloose_nod0_2e12_lhvloose_nod0_L1EM20VH_3EM10VH = {fReader, "lepHLT_e24_lhvloose_nod0_2e12_lhvloose_nod0_L1EM20VH_3EM10VH"};
   TTreeReaderValue<vector<bool>> lepHLT_e12_lhloose_nod0_2mu10 = {fReader, "lepHLT_e12_lhloose_nod0_2mu10"};
   TTreeReaderValue<vector<bool>> lepHLT_2e12_lhloose_nod0_mu10 = {fReader, "lepHLT_2e12_lhloose_nod0_mu10"};
   TTreeReaderValue<vector<bool>> lepHLT_mu20_2mu4noL1 = {fReader, "lepHLT_mu20_2mu4noL1"};
   TTreeReaderValue<vector<bool>> lepHLT_3mu6 = {fReader, "lepHLT_3mu6"};
   TTreeReaderValue<vector<bool>> lepHLT_3mu6_msonly = {fReader, "lepHLT_3mu6_msonly"};
   TTreeReaderValue<vector<bool>> lepHLT_e17_lhloose_nod0_2e9_lhloose_nod0 = {fReader, "lepHLT_e17_lhloose_nod0_2e9_lhloose_nod0"};
   TTreeReaderValue<vector<bool>> lepHLT_e17_lhloose_nod0_2e10_lhloose_nod0_L1EM15VH_3EM8VH = {fReader, "lepHLT_e17_lhloose_nod0_2e10_lhloose_nod0_L1EM15VH_3EM8VH"};
   TTreeReaderValue<vector<bool>> lepHLT_e17_lhloose_2e9_lhloose = {fReader, "lepHLT_e17_lhloose_2e9_lhloose"};
   TTreeReaderValue<vector<bool>> lepHLT_mu18_2mu4noL1 = {fReader, "lepHLT_mu18_2mu4noL1"};
   TTreeReaderValue<vector<bool>> lepHLT_2e12_lhloose_mu10 = {fReader, "lepHLT_2e12_lhloose_mu10"};
   TTreeReaderValue<vector<bool>> lepHLT_e12_lhloose_2mu10 = {fReader, "lepHLT_e12_lhloose_2mu10"};
   TTreeReaderValue<vector<bool>> lepHLT_e24_lhmedium_L1EM20VH = {fReader, "lepHLT_e24_lhmedium_L1EM20VH"};
   TTreeReaderValue<vector<bool>> lepHLT_e60_lhmedium = {fReader, "lepHLT_e60_lhmedium"};
   TTreeReaderValue<vector<bool>> lepHLT_e120_lhloose = {fReader, "lepHLT_e120_lhloose"};
   TTreeReaderValue<vector<bool>> lepHLT_mu20_iloose_L1MU15 = {fReader, "lepHLT_mu20_iloose_L1MU15"};
   TTreeReaderValue<vector<bool>> lepHLT_2e12_lhloose_L12EM10VH = {fReader, "lepHLT_2e12_lhloose_L12EM10VH"};
   TTreeReaderValue<vector<bool>> lepHLT_mu18_mu8noL1 = {fReader, "lepHLT_mu18_mu8noL1"};
   TTreeReaderValue<vector<bool>> lepHLT_e17_lhloose_mu14 = {fReader, "lepHLT_e17_lhloose_mu14"};
   TTreeReaderValue<vector<bool>> lepHLT_e7_lhmedium_mu24 = {fReader, "lepHLT_e7_lhmedium_mu24"};
   TTreeReaderValue<vector<bool>> lepHLT_e24_lhtight_nod0_ivarloose = {fReader, "lepHLT_e24_lhtight_nod0_ivarloose"};
   TTreeReaderValue<vector<bool>> lepHLT_e24_lhmedium_nod0_L1EM20VH = {fReader, "lepHLT_e24_lhmedium_nod0_L1EM20VH"};
   TTreeReaderValue<vector<bool>> lepHLT_e60_medium = {fReader, "lepHLT_e60_medium"};
   TTreeReaderValue<vector<bool>> lepHLT_mu40 = {fReader, "lepHLT_mu40"};
   TTreeReaderValue<vector<bool>> lepHLT_mu24_iloose_L1MU15 = {fReader, "lepHLT_mu24_iloose_L1MU15"};
   TTreeReaderValue<vector<bool>> lepHLT_mu24_ivarloose_L1MU15 = {fReader, "lepHLT_mu24_ivarloose_L1MU15"};
   TTreeReaderValue<vector<bool>> lepHLT_mu24_ivarmedium = {fReader, "lepHLT_mu24_ivarmedium"};
   TTreeReaderValue<vector<bool>> lepHLT_mu24_imedium = {fReader, "lepHLT_mu24_imedium"};
   TTreeReaderValue<vector<bool>> lepHLT_mu26_imedium = {fReader, "lepHLT_mu26_imedium"};
   TTreeReaderValue<vector<bool>> lepHLT_2e15_lhvloose_nod0_L12EM13VH = {fReader, "lepHLT_2e15_lhvloose_nod0_L12EM13VH"};
   TTreeReaderValue<vector<bool>> lepHLT_2e17_lhvloose_nod0 = {fReader, "lepHLT_2e17_lhvloose_nod0"};
   TTreeReaderValue<vector<bool>> lepHLT_2mu10 = {fReader, "lepHLT_2mu10"};
   TTreeReaderValue<vector<bool>> lepHLT_2mu14 = {fReader, "lepHLT_2mu14"};
   TTreeReaderValue<vector<bool>> lepHLT_mu20_mu8noL1 = {fReader, "lepHLT_mu20_mu8noL1"};
   TTreeReaderValue<vector<bool>> lepHLT_mu22_mu8noL1 = {fReader, "lepHLT_mu22_mu8noL1"};
   TTreeReaderValue<vector<bool>> lepHLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1 = {fReader, "lepHLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1"};
   TTreeReaderValue<vector<bool>> lepHLT_e26_lhtight_nod0_ivarloose = {fReader, "lepHLT_e26_lhtight_nod0_ivarloose"};
   TTreeReaderValue<vector<bool>> lepHLT_e26_lhtight_nod0 = {fReader, "lepHLT_e26_lhtight_nod0"};
   TTreeReaderValue<vector<bool>> lepHLT_e60_lhmedium_nod0 = {fReader, "lepHLT_e60_lhmedium_nod0"};
   TTreeReaderValue<vector<bool>> lepHLT_e140_lhloose_nod0 = {fReader, "lepHLT_e140_lhloose_nod0"};
   TTreeReaderValue<vector<bool>> lepHLT_e300_etcut = {fReader, "lepHLT_e300_etcut"};
   TTreeReaderValue<vector<bool>> lepHLT_mu26_ivarmedium = {fReader, "lepHLT_mu26_ivarmedium"};
   TTreeReaderValue<vector<bool>> lepHLT_mu50 = {fReader, "lepHLT_mu50"};
   TTreeReaderValue<vector<bool>> lepHLT_mu60_0eta105_msonly = {fReader, "lepHLT_mu60_0eta105_msonly"};
   TTreeReaderValue<vector<bool>> lepHLT_2e17_lhvloose_nod0_L12EM15VHI = {fReader, "lepHLT_2e17_lhvloose_nod0_L12EM15VHI"};
   TTreeReaderValue<vector<bool>> lepHLT_2e24_lhvloose_nod0 = {fReader, "lepHLT_2e24_lhvloose_nod0"};
   TTreeReaderValue<vector<bool>> lepHLT_e17_lhloose_nod0_mu14 = {fReader, "lepHLT_e17_lhloose_nod0_mu14"};
   TTreeReaderValue<vector<bool>> lepHLT_e26_lhmedium_nod0_mu8noL1 = {fReader, "lepHLT_e26_lhmedium_nod0_mu8noL1"};
   TTreeReaderValue<vector<bool>> lepHLT_e7_lhmedium_nod0_mu24 = {fReader, "lepHLT_e7_lhmedium_nod0_mu24"};
   TTreeReaderValue<Bool_t> trigMatch_1L2LTrig = {fReader, "trigMatch_1L2LTrig"};
   TTreeReaderValue<Bool_t> trigMatch_3LTrig = {fReader, "trigMatch_3LTrig"};
   TTreeReaderValue<Bool_t> trigMatch_1LTrig = {fReader, "trigMatch_1LTrig"};
   TTreeReaderValue<Bool_t> trigMatch_3LTrigOR = {fReader, "trigMatch_3LTrigOR"};
   TTreeReaderValue<Bool_t> trigMatch_1L2LTrigOR = {fReader, "trigMatch_1L2LTrigOR"};
   TTreeReaderValue<Bool_t> trigMatch_1LTrigOR = {fReader, "trigMatch_1LTrigOR"};
   TTreeReaderValue<Bool_t> trigMatch_2LTrig = {fReader, "trigMatch_2LTrig"};
   TTreeReaderValue<Bool_t> trigMatch_2LTrigOR = {fReader, "trigMatch_2LTrigOR"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_e24_lhmedium_L1EM20VH = {fReader, "trigMatch_HLT_e24_lhmedium_L1EM20VH"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_e60_lhmedium = {fReader, "trigMatch_HLT_e60_lhmedium"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_e120_lhloose = {fReader, "trigMatch_HLT_e120_lhloose"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_mu20_iloose_L1MU15 = {fReader, "trigMatch_HLT_mu20_iloose_L1MU15"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_2e12_lhloose_L12EM10VH = {fReader, "trigMatch_HLT_2e12_lhloose_L12EM10VH"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_mu18_mu8noL1 = {fReader, "trigMatch_HLT_mu18_mu8noL1"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_e17_lhloose_mu14 = {fReader, "trigMatch_HLT_e17_lhloose_mu14"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_e7_lhmedium_mu24 = {fReader, "trigMatch_HLT_e7_lhmedium_mu24"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_e24_lhtight_nod0_ivarloose = {fReader, "trigMatch_HLT_e24_lhtight_nod0_ivarloose"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_e24_lhmedium_nod0_L1EM20VH = {fReader, "trigMatch_HLT_e24_lhmedium_nod0_L1EM20VH"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_e60_medium = {fReader, "trigMatch_HLT_e60_medium"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_mu40 = {fReader, "trigMatch_HLT_mu40"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_mu24_iloose_L1MU15 = {fReader, "trigMatch_HLT_mu24_iloose_L1MU15"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_mu24_ivarloose_L1MU15 = {fReader, "trigMatch_HLT_mu24_ivarloose_L1MU15"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_mu24_ivarmedium = {fReader, "trigMatch_HLT_mu24_ivarmedium"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_mu24_imedium = {fReader, "trigMatch_HLT_mu24_imedium"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_mu26_imedium = {fReader, "trigMatch_HLT_mu26_imedium"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_2e15_lhvloose_nod0_L12EM13VH = {fReader, "trigMatch_HLT_2e15_lhvloose_nod0_L12EM13VH"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_2e17_lhvloose_nod0 = {fReader, "trigMatch_HLT_2e17_lhvloose_nod0"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_2mu10 = {fReader, "trigMatch_HLT_2mu10"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_2mu14 = {fReader, "trigMatch_HLT_2mu14"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_mu20_mu8noL1 = {fReader, "trigMatch_HLT_mu20_mu8noL1"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_mu22_mu8noL1 = {fReader, "trigMatch_HLT_mu22_mu8noL1"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1 = {fReader, "trigMatch_HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_e26_lhtight_nod0_ivarloose = {fReader, "trigMatch_HLT_e26_lhtight_nod0_ivarloose"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_e26_lhtight_nod0 = {fReader, "trigMatch_HLT_e26_lhtight_nod0"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_e60_lhmedium_nod0 = {fReader, "trigMatch_HLT_e60_lhmedium_nod0"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_e140_lhloose_nod0 = {fReader, "trigMatch_HLT_e140_lhloose_nod0"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_e300_etcut = {fReader, "trigMatch_HLT_e300_etcut"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_mu26_ivarmedium = {fReader, "trigMatch_HLT_mu26_ivarmedium"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_mu50 = {fReader, "trigMatch_HLT_mu50"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_mu60_0eta105_msonly = {fReader, "trigMatch_HLT_mu60_0eta105_msonly"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_2e17_lhvloose_nod0_L12EM15VHI = {fReader, "trigMatch_HLT_2e17_lhvloose_nod0_L12EM15VHI"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_2e24_lhvloose_nod0 = {fReader, "trigMatch_HLT_2e24_lhvloose_nod0"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_e17_lhloose_nod0_mu14 = {fReader, "trigMatch_HLT_e17_lhloose_nod0_mu14"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_e26_lhmedium_nod0_mu8noL1 = {fReader, "trigMatch_HLT_e26_lhmedium_nod0_mu8noL1"};
   TTreeReaderValue<Bool_t> trigMatch_HLT_e7_lhmedium_nod0_mu24 = {fReader, "trigMatch_HLT_e7_lhmedium_nod0_mu24"};
   TTreeReaderValue<Float_t> mu = {fReader, "mu"};
   TTreeReaderValue<Float_t> avg_mu = {fReader, "avg_mu"};
   TTreeReaderValue<Float_t> actual_mu = {fReader, "actual_mu"};
   TTreeReaderValue<Int_t> nVtx = {fReader, "nVtx"};
   TTreeReaderValue<Int_t> channel = {fReader, "channel"};
   TTreeReaderValue<Int_t> nLep_base = {fReader, "nLep_base"};
   TTreeReaderValue<Int_t> nLep_signal = {fReader, "nLep_signal"};
   TTreeReaderArray<int> lepFlavor = {fReader, "lepFlavor"};
   TTreeReaderArray<int> lepCharge = {fReader, "lepCharge"};
   TTreeReaderArray<int> lepAuthor = {fReader, "lepAuthor"};
   TTreeReaderArray<float> lepPt = {fReader, "lepPt"};
   TTreeReaderArray<float> lepEta = {fReader, "lepEta"};
   TTreeReaderArray<float> lepPhi = {fReader, "lepPhi"};
   TTreeReaderArray<float> lepM = {fReader, "lepM"};
   TTreeReaderArray<float> lepD0 = {fReader, "lepD0"};
   TTreeReaderArray<float> lepD0Sig = {fReader, "lepD0Sig"};
   TTreeReaderArray<float> lepZ0 = {fReader, "lepZ0"};
   TTreeReaderArray<float> lepZ0SinTheta = {fReader, "lepZ0SinTheta"};
   TTreeReaderArray<float> lepPtcone20 = {fReader, "lepPtcone20"};
   TTreeReaderArray<float> lepPtcone30 = {fReader, "lepPtcone30"};
   TTreeReaderArray<float> lepPtcone40 = {fReader, "lepPtcone40"};
   TTreeReaderArray<float> lepTopoetcone20 = {fReader, "lepTopoetcone20"};
   TTreeReaderArray<float> lepTopoetcone30 = {fReader, "lepTopoetcone30"};
   TTreeReaderArray<float> lepTopoetcone40 = {fReader, "lepTopoetcone40"};
   TTreeReaderArray<float> lepPtvarcone20 = {fReader, "lepPtvarcone20"};
   TTreeReaderArray<float> lepPtvarcone30 = {fReader, "lepPtvarcone30"};
   TTreeReaderArray<float> lepPtvarcone40 = {fReader, "lepPtvarcone40"};
   TTreeReaderArray<float> lepPtvarcone30_TightTTVA_pt1000 = {fReader, "lepPtvarcone30_TightTTVA_pt1000"};
   TTreeReaderArray<float> lepPtvarcone30_TightTTVA_pt500 = {fReader, "lepPtvarcone30_TightTTVA_pt500"};
   TTreeReaderArray<float> lepPtvarcone20_TightTTVA_pt1000 = {fReader, "lepPtvarcone20_TightTTVA_pt1000"};
   TTreeReaderArray<float> lepPtvarcone20_TightTTVA_pt500 = {fReader, "lepPtvarcone20_TightTTVA_pt500"};
   TTreeReaderArray<float> lepPtcone20_TightTTVALooseCone_pt500 = {fReader, "lepPtcone20_TightTTVALooseCone_pt500"};
   TTreeReaderArray<float> lepPtcone20_TightTTVALooseCone_pt1000 = {fReader, "lepPtcone20_TightTTVALooseCone_pt1000"};
   TTreeReaderArray<float> lepPtcone20_TightTTVA_pt500 = {fReader, "lepPtcone20_TightTTVA_pt500"};
   TTreeReaderArray<float> lepPtcone20_TightTTVA_pt1000 = {fReader, "lepPtcone20_TightTTVA_pt1000"};
   TTreeReaderArray<float> lepNeflowisol20 = {fReader, "lepNeflowisol20"};
   TTreeReaderValue<vector<bool>> lepPassOR = {fReader, "lepPassOR"};
   TTreeReaderArray<int> lepType = {fReader, "lepType"};
   TTreeReaderArray<int> lepOrigin = {fReader, "lepOrigin"};
   TTreeReaderArray<int> lepIFFClass = {fReader, "lepIFFClass"};
   TTreeReaderArray<int> lepEgMotherType = {fReader, "lepEgMotherType"};
   TTreeReaderArray<int> lepEgMotherOrigin = {fReader, "lepEgMotherOrigin"};
   TTreeReaderArray<int> lepEgMotherPdgId = {fReader, "lepEgMotherPdgId"};
   TTreeReaderArray<int> lepECIDS = {fReader, "lepECIDS"};
   TTreeReaderValue<vector<bool>> lepPassBL = {fReader, "lepPassBL"};
   TTreeReaderValue<vector<bool>> lepVeryLoose = {fReader, "lepVeryLoose"};
   TTreeReaderValue<vector<bool>> lepLoose = {fReader, "lepLoose"};
   TTreeReaderValue<vector<bool>> lepMedium = {fReader, "lepMedium"};
   TTreeReaderValue<vector<bool>> lepTight = {fReader, "lepTight"};
   TTreeReaderValue<vector<bool>> lepHighPt = {fReader, "lepHighPt"};
   TTreeReaderValue<vector<bool>> lepIsoFCHighPtCaloOnly = {fReader, "lepIsoFCHighPtCaloOnly"};
   TTreeReaderValue<vector<bool>> lepIsoFCLoose = {fReader, "lepIsoFCLoose"};
   TTreeReaderValue<vector<bool>> lepIsoFCTight = {fReader, "lepIsoFCTight"};
   TTreeReaderValue<vector<bool>> lepIsoFCLoose_FixedRad = {fReader, "lepIsoFCLoose_FixedRad"};
   TTreeReaderValue<vector<bool>> lepIsoFCTight_FixedRad = {fReader, "lepIsoFCTight_FixedRad"};
   TTreeReaderValue<vector<bool>> lepIsoHighPtCaloOnly = {fReader, "lepIsoHighPtCaloOnly"};
   TTreeReaderValue<vector<bool>> lepIsoTightTrackOnly_VarRad = {fReader, "lepIsoTightTrackOnly_VarRad"};
   TTreeReaderValue<vector<bool>> lepIsoTightTrackOnly_FixedRad = {fReader, "lepIsoTightTrackOnly_FixedRad"};
   TTreeReaderValue<vector<bool>> lepIsoLoose_VarRad = {fReader, "lepIsoLoose_VarRad"};
   TTreeReaderValue<vector<bool>> lepIsoTight_VarRad = {fReader, "lepIsoTight_VarRad"};
   TTreeReaderValue<vector<bool>> lepIsoPLVLoose = {fReader, "lepIsoPLVLoose"};
   TTreeReaderValue<vector<bool>> lepIsoPLVTight = {fReader, "lepIsoPLVTight"};
   TTreeReaderValue<vector<bool>> lepTruthMatched = {fReader, "lepTruthMatched"};
   TTreeReaderArray<int> lepTruthCharge = {fReader, "lepTruthCharge"};
   TTreeReaderArray<float> lepTruthPt = {fReader, "lepTruthPt"};
   TTreeReaderArray<float> lepTruthEta = {fReader, "lepTruthEta"};
   TTreeReaderArray<float> lepTruthPhi = {fReader, "lepTruthPhi"};
   TTreeReaderArray<float> lepTruthM = {fReader, "lepTruthM"};
   TTreeReaderArray<float> lepTrigSF = {fReader, "lepTrigSF"};
   TTreeReaderArray<float> lepRecoSF = {fReader, "lepRecoSF"};
   TTreeReaderArray<float> lepBLTrigSF = {fReader, "lepBLTrigSF"};
   TTreeReaderArray<float> lepBLRecoSF = {fReader, "lepBLRecoSF"};
   TTreeReaderValue<Int_t> nJet30 = {fReader, "nJet30"};
   TTreeReaderValue<Int_t> nJet20 = {fReader, "nJet20"};
   TTreeReaderArray<float> jetPt = {fReader, "jetPt"};
   TTreeReaderArray<float> jetEta = {fReader, "jetEta"};
   TTreeReaderArray<float> jetPhi = {fReader, "jetPhi"};
   TTreeReaderArray<float> jetM = {fReader, "jetM"};
   TTreeReaderArray<float> jetJVT = {fReader, "jetJVT"};
   TTreeReaderValue<vector<bool>> jetPassOR = {fReader, "jetPassOR"};
   TTreeReaderValue<vector<bool>> jetSignal = {fReader, "jetSignal"};
   TTreeReaderValue<Float_t> mjj = {fReader, "mjj"};
   TTreeReaderArray<float> jetTileEnergy = {fReader, "jetTileEnergy"};
   TTreeReaderArray<float> jetMV2c10 = {fReader, "jetMV2c10"};
   TTreeReaderArray<float> jetdl1 = {fReader, "jetdl1"};
   TTreeReaderArray<float> jetdl1r = {fReader, "jetdl1r"};
   TTreeReaderValue<Float_t> met_Et = {fReader, "met_Et"};
   TTreeReaderValue<Float_t> met_Sign = {fReader, "met_Sign"};
   TTreeReaderValue<Float_t> met_Phi = {fReader, "met_Phi"};
   TTreeReaderValue<Float_t> mll = {fReader, "mll"};
   TTreeReaderValue<Double_t> pileupWeight = {fReader, "pileupWeight"};
   TTreeReaderValue<Double_t> leptonWeight = {fReader, "leptonWeight"};
   TTreeReaderValue<Double_t> baselineleptonWeight = {fReader, "baselineleptonWeight"};
   TTreeReaderValue<Double_t> eventWeight = {fReader, "eventWeight"};
   TTreeReaderValue<Double_t> bTagWeight = {fReader, "bTagWeight"};
   TTreeReaderValue<Double_t> jvtWeight = {fReader, "jvtWeight"};
   TTreeReaderValue<Double_t> globalDiLepTrigSF = {fReader, "globalDiLepTrigSF"};
   TTreeReaderValue<Double_t> globalBaselineDiLepTrigSF = {fReader, "globalBaselineDiLepTrigSF"};
   TTreeReaderValue<Double_t> singleLepTrigSF = {fReader, "singleLepTrigSF"};
   TTreeReaderValue<Double_t> singleBaselineLepTrigSF = {fReader, "singleBaselineLepTrigSF"};
   TTreeReaderValue<Double_t> flavSymWeight = {fReader, "flavSymWeight"};
   TTreeReaderValue<Double_t> genWeight = {fReader, "genWeight"};
   TTreeReaderValue<Double_t> genWeightUp = {fReader, "genWeightUp"};
   TTreeReaderValue<Double_t> genWeightDown = {fReader, "genWeightDown"};
   TTreeReaderValue<ULong64_t> PRWHash = {fReader, "PRWHash"};
   TTreeReaderValue<ULong64_t> EventNumber = {fReader, "EventNumber"};
   TTreeReaderValue<Float_t> xsec = {fReader, "xsec"};
   TTreeReaderValue<Float_t> GenHt = {fReader, "GenHt"};
   TTreeReaderValue<Float_t> GenMET = {fReader, "GenMET"};
   TTreeReaderValue<Int_t> DatasetNumber = {fReader, "DatasetNumber"};
   TTreeReaderValue<Int_t> RunNumber = {fReader, "RunNumber"};
   TTreeReaderValue<Int_t> RandomRunNumber = {fReader, "RandomRunNumber"};
   TTreeReaderValue<Int_t> FS = {fReader, "FS"};
   TTreeReaderArray<float> LHE3Weights = {fReader, "LHE3Weights"};
   TTreeReaderArray<TString> LHE3WeightNames = {fReader, "LHE3WeightNames"};
  

  //MY REL21 ANALYSIS STOPS HERE
  // >----

  /**
     TTreeReaderValue<Float_t> leptonWeight_EL_CHARGEID_STAT__1down = {fReader, "leptonWeight_EL_CHARGEID_STAT__1down"};
     TTreeReaderValue<Float_t> leptonWeight_EL_CHARGEID_STAT__1up = {fReader, "leptonWeight_EL_CHARGEID_STAT__1up"};
     TTreeReaderValue<Float_t> leptonWeight_EL_CHARGEID_SYStotal__1down = {fReader, "leptonWeight_EL_CHARGEID_SYStotal__1down"};
     TTreeReaderValue<Float_t> leptonWeight_EL_CHARGEID_SYStotal__1up = {fReader, "leptonWeight_EL_CHARGEID_SYStotal__1up"};
     TTreeReaderValue<Float_t> leptonWeight_EL_EFF_ChargeIDSel_TOTAL_1NPCOR_PLUS_UNCOR__1down = {fReader, "leptonWeight_EL_EFF_ChargeIDSel_TOTAL_1NPCOR_PLUS_UNCOR__1down"};
     TTreeReaderValue<Float_t> leptonWeight_EL_EFF_ChargeIDSel_TOTAL_1NPCOR_PLUS_UNCOR__1up = {fReader, "leptonWeight_EL_EFF_ChargeIDSel_TOTAL_1NPCOR_PLUS_UNCOR__1up"};
     TTreeReaderValue<Float_t> leptonWeight_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR__1down = {fReader, "leptonWeight_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR__1down"};
     TTreeReaderValue<Float_t> leptonWeight_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR__1up = {fReader, "leptonWeight_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR__1up"};
     TTreeReaderValue<Float_t> leptonWeight_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR__1down = {fReader, "leptonWeight_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR__1down"};
     TTreeReaderValue<Float_t> leptonWeight_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR__1up = {fReader, "leptonWeight_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR__1up"};
     TTreeReaderValue<Float_t> leptonWeight_EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR__1down = {fReader, "leptonWeight_EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR__1down"};
     TTreeReaderValue<Float_t> leptonWeight_EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR__1up = {fReader, "leptonWeight_EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR__1up"};
     TTreeReaderValue<Float_t> leptonWeight_MUON_EFF_ISO_STAT__1down = {fReader, "leptonWeight_MUON_EFF_ISO_STAT__1down"};
     TTreeReaderValue<Float_t> leptonWeight_MUON_EFF_ISO_STAT__1up = {fReader, "leptonWeight_MUON_EFF_ISO_STAT__1up"};
     TTreeReaderValue<Float_t> leptonWeight_MUON_EFF_ISO_SYS__1down = {fReader, "leptonWeight_MUON_EFF_ISO_SYS__1down"};
     TTreeReaderValue<Float_t> leptonWeight_MUON_EFF_ISO_SYS__1up = {fReader, "leptonWeight_MUON_EFF_ISO_SYS__1up"};
     TTreeReaderValue<Float_t> leptonWeight_MUON_EFF_RECO_STAT__1down = {fReader, "leptonWeight_MUON_EFF_RECO_STAT__1down"};
     TTreeReaderValue<Float_t> leptonWeight_MUON_EFF_RECO_STAT__1up = {fReader, "leptonWeight_MUON_EFF_RECO_STAT__1up"};
     TTreeReaderValue<Float_t> leptonWeight_MUON_EFF_RECO_STAT_LOWPT__1down = {fReader, "leptonWeight_MUON_EFF_RECO_STAT_LOWPT__1down"};
     TTreeReaderValue<Float_t> leptonWeight_MUON_EFF_RECO_STAT_LOWPT__1up = {fReader, "leptonWeight_MUON_EFF_RECO_STAT_LOWPT__1up"};
     TTreeReaderValue<Float_t> leptonWeight_MUON_EFF_RECO_SYS__1down = {fReader, "leptonWeight_MUON_EFF_RECO_SYS__1down"};
     TTreeReaderValue<Float_t> leptonWeight_MUON_EFF_RECO_SYS__1up = {fReader, "leptonWeight_MUON_EFF_RECO_SYS__1up"};
     TTreeReaderValue<Float_t> leptonWeight_MUON_EFF_RECO_SYS_LOWPT__1down = {fReader, "leptonWeight_MUON_EFF_RECO_SYS_LOWPT__1down"};
     TTreeReaderValue<Float_t> leptonWeight_MUON_EFF_RECO_SYS_LOWPT__1up = {fReader, "leptonWeight_MUON_EFF_RECO_SYS_LOWPT__1up"};
     TTreeReaderValue<Float_t> leptonWeight_MUON_EFF_TTVA_STAT__1down = {fReader, "leptonWeight_MUON_EFF_TTVA_STAT__1down"};
     TTreeReaderValue<Float_t> leptonWeight_MUON_EFF_TTVA_STAT__1up = {fReader, "leptonWeight_MUON_EFF_TTVA_STAT__1up"};
     TTreeReaderValue<Float_t> leptonWeight_MUON_EFF_TTVA_SYS__1down = {fReader, "leptonWeight_MUON_EFF_TTVA_SYS__1down"};
     TTreeReaderValue<Float_t> leptonWeight_MUON_EFF_TTVA_SYS__1up = {fReader, "leptonWeight_MUON_EFF_TTVA_SYS__1up"};
     TTreeReaderValue<Float_t> trigWeight_EL_EFF_TriggerEff_TOTAL_1NPCOR_PLUS_UNCOR__1down = {fReader, "trigWeight_EL_EFF_TriggerEff_TOTAL_1NPCOR_PLUS_UNCOR__1down"};
     TTreeReaderValue<Float_t> trigWeight_EL_EFF_TriggerEff_TOTAL_1NPCOR_PLUS_UNCOR__1up = {fReader, "trigWeight_EL_EFF_TriggerEff_TOTAL_1NPCOR_PLUS_UNCOR__1up"};
     TTreeReaderValue<Float_t> trigWeight_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1down = {fReader, "trigWeight_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1down"};
     TTreeReaderValue<Float_t> trigWeight_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1up = {fReader, "trigWeight_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1up"};
     TTreeReaderValue<Float_t> trigWeight_MUON_EFF_TrigStatUncertainty__1down = {fReader, "trigWeight_MUON_EFF_TrigStatUncertainty__1down"};
     TTreeReaderValue<Float_t> trigWeight_MUON_EFF_TrigStatUncertainty__1up = {fReader, "trigWeight_MUON_EFF_TrigStatUncertainty__1up"};
     TTreeReaderValue<Float_t> trigWeight_MUON_EFF_TrigSystUncertainty__1down = {fReader, "trigWeight_MUON_EFF_TrigSystUncertainty__1down"};
     TTreeReaderValue<Float_t> trigWeight_MUON_EFF_TrigSystUncertainty__1up = {fReader, "trigWeight_MUON_EFF_TrigSystUncertainty__1up"};
     TTreeReaderValue<Float_t> bTagWeight_FT_EFF_B_systematics__1down = {fReader, "bTagWeight_FT_EFF_B_systematics__1down"};
     TTreeReaderValue<Float_t> bTagWeight_FT_EFF_B_systematics__1up = {fReader, "bTagWeight_FT_EFF_B_systematics__1up"};
     TTreeReaderValue<Float_t> bTagWeight_FT_EFF_C_systematics__1down = {fReader, "bTagWeight_FT_EFF_C_systematics__1down"};
     TTreeReaderValue<Float_t> bTagWeight_FT_EFF_C_systematics__1up = {fReader, "bTagWeight_FT_EFF_C_systematics__1up"};
     TTreeReaderValue<Float_t> bTagWeight_FT_EFF_Light_systematics__1down = {fReader, "bTagWeight_FT_EFF_Light_systematics__1down"};
     TTreeReaderValue<Float_t> bTagWeight_FT_EFF_Light_systematics__1up = {fReader, "bTagWeight_FT_EFF_Light_systematics__1up"};
     TTreeReaderValue<Float_t> bTagWeight_FT_EFF_extrapolation__1down = {fReader, "bTagWeight_FT_EFF_extrapolation__1down"};
     TTreeReaderValue<Float_t> bTagWeight_FT_EFF_extrapolation__1up = {fReader, "bTagWeight_FT_EFF_extrapolation__1up"};
     TTreeReaderValue<Float_t> bTagWeight_FT_EFF_extrapolation_from_charm__1down = {fReader, "bTagWeight_FT_EFF_extrapolation_from_charm__1down"};
     TTreeReaderValue<Float_t> bTagWeight_FT_EFF_extrapolation_from_charm__1up = {fReader, "bTagWeight_FT_EFF_extrapolation_from_charm__1up"};
     TTreeReaderValue<Float_t> jvtWeight_JET_JvtEfficiency__1down = {fReader, "jvtWeight_JET_JvtEfficiency__1down"};
     TTreeReaderValue<Float_t> jvtWeight_JET_JvtEfficiency__1up = {fReader, "jvtWeight_JET_JvtEfficiency__1up"};
     TTreeReaderValue<Float_t> jvtWeight_JET_fJvtEfficiency__1down = {fReader, "jvtWeight_JET_fJvtEfficiency__1down"};
     TTreeReaderValue<Float_t> jvtWeight_JET_fJvtEfficiency__1up = {fReader, "jvtWeight_JET_fJvtEfficiency__1up"};
     TTreeReaderValue<Float_t> pileupWeightUp = {fReader, "pileupWeightUp"};
     TTreeReaderValue<Float_t> pileupWeightDown = {fReader, "pileupWeightDown"};
  */
  /**
     TTreeReaderArray<float> LHE3Weights = {fReader, "LHE3Weights"};
     TTreeReaderArray<TString> LHE3WeightNames = {fReader, "LHE3WeightNames"};
  */
  
  
  // Readers to access the data (delete the ones you do not need).
  /**
     TTreeReaderValue<vector<bool>> lepHLT_e24_lhmedium_L1EM20VH = {fReader, "lepHLT_e24_lhmedium_L1EM20VH"};
     TTreeReaderValue<vector<bool>> lepHLT_e60_lhmedium = {fReader, "lepHLT_e60_lhmedium"};
     TTreeReaderValue<vector<bool>> lepHLT_e120_lhloose = {fReader, "lepHLT_e120_lhloose"};
     TTreeReaderValue<vector<bool>> lepHLT_mu20_iloose_L1MU15 = {fReader, "lepHLT_mu20_iloose_L1MU15"};
     TTreeReaderValue<vector<bool>> lepHLT_2e12_lhloose_L12EM10VH = {fReader, "lepHLT_2e12_lhloose_L12EM10VH"};
     TTreeReaderValue<vector<bool>> lepHLT_mu18_mu8noL1 = {fReader, "lepHLT_mu18_mu8noL1"};
     TTreeReaderValue<vector<bool>> lepHLT_e17_lhloose_mu14 = {fReader, "lepHLT_e17_lhloose_mu14"};
     TTreeReaderValue<vector<bool>> lepHLT_e7_lhmedium_mu24 = {fReader, "lepHLT_e7_lhmedium_mu24"};
     TTreeReaderValue<vector<bool>> lepHLT_e24_lhtight_nod0_ivarloose = {fReader, "lepHLT_e24_lhtight_nod0_ivarloose"};
     TTreeReaderValue<vector<bool>> lepHLT_e24_lhmedium_nod0_L1EM20VH = {fReader, "lepHLT_e24_lhmedium_nod0_L1EM20VH"};
     TTreeReaderValue<vector<bool>> lepHLT_e60_medium = {fReader, "lepHLT_e60_medium"};
     TTreeReaderValue<vector<bool>> lepHLT_mu40 = {fReader, "lepHLT_mu40"};
     TTreeReaderValue<vector<bool>> lepHLT_mu24_iloose_L1MU15 = {fReader, "lepHLT_mu24_iloose_L1MU15"};
     TTreeReaderValue<vector<bool>> lepHLT_mu24_ivarloose_L1MU15 = {fReader, "lepHLT_mu24_ivarloose_L1MU15"};
     TTreeReaderValue<vector<bool>> lepHLT_mu24_ivarmedium = {fReader, "lepHLT_mu24_ivarmedium"};
     TTreeReaderValue<vector<bool>> lepHLT_mu24_imedium = {fReader, "lepHLT_mu24_imedium"};
     TTreeReaderValue<vector<bool>> lepHLT_mu26_imedium = {fReader, "lepHLT_mu26_imedium"};
     TTreeReaderValue<vector<bool>> lepHLT_2e15_lhvloose_nod0_L12EM13VH = {fReader, "lepHLT_2e15_lhvloose_nod0_L12EM13VH"};
     TTreeReaderValue<vector<bool>> lepHLT_2e17_lhvloose_nod0 = {fReader, "lepHLT_2e17_lhvloose_nod0"};
     TTreeReaderValue<vector<bool>> lepHLT_2mu10 = {fReader, "lepHLT_2mu10"};
     TTreeReaderValue<vector<bool>> lepHLT_2mu14 = {fReader, "lepHLT_2mu14"};
     TTreeReaderValue<vector<bool>> lepHLT_mu20_mu8noL1 = {fReader, "lepHLT_mu20_mu8noL1"};
     TTreeReaderValue<vector<bool>> lepHLT_mu22_mu8noL1 = {fReader, "lepHLT_mu22_mu8noL1"};
     TTreeReaderValue<vector<bool>> lepHLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1 = {fReader, "lepHLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1"};
     TTreeReaderValue<vector<bool>> lepHLT_e26_lhtight_nod0_ivarloose = {fReader, "lepHLT_e26_lhtight_nod0_ivarloose"};
     TTreeReaderValue<vector<bool>> lepHLT_e26_lhtight_nod0 = {fReader, "lepHLT_e26_lhtight_nod0"};
     TTreeReaderValue<vector<bool>> lepHLT_e60_lhmedium_nod0 = {fReader, "lepHLT_e60_lhmedium_nod0"};
     TTreeReaderValue<vector<bool>> lepHLT_e140_lhloose_nod0 = {fReader, "lepHLT_e140_lhloose_nod0"};
     TTreeReaderValue<vector<bool>> lepHLT_e300_etcut = {fReader, "lepHLT_e300_etcut"};
     TTreeReaderValue<vector<bool>> lepHLT_mu26_ivarmedium = {fReader, "lepHLT_mu26_ivarmedium"};
     TTreeReaderValue<vector<bool>> lepHLT_mu50 = {fReader, "lepHLT_mu50"};
     TTreeReaderValue<vector<bool>> lepHLT_mu60_0eta105_msonly = {fReader, "lepHLT_mu60_0eta105_msonly"};
     TTreeReaderValue<vector<bool>> lepHLT_2e17_lhvloose_nod0_L12EM15VHI = {fReader, "lepHLT_2e17_lhvloose_nod0_L12EM15VHI"};
     TTreeReaderValue<vector<bool>> lepHLT_2e24_lhvloose_nod0 = {fReader, "lepHLT_2e24_lhvloose_nod0"};
     TTreeReaderValue<vector<bool>> lepHLT_e17_lhloose_nod0_mu14 = {fReader, "lepHLT_e17_lhloose_nod0_mu14"};
     TTreeReaderValue<vector<bool>> lepHLT_e26_lhmedium_nod0_mu8noL1 = {fReader, "lepHLT_e26_lhmedium_nod0_mu8noL1"};
     TTreeReaderValue<vector<bool>> lepHLT_e7_lhmedium_nod0_mu24 = {fReader, "lepHLT_e7_lhmedium_nod0_mu24"};
     TTreeReaderValue<Bool_t> trigMatch_1L2LTrig = {fReader, "trigMatch_1L2LTrig"};
     TTreeReaderValue<Bool_t> trigMatch_1LTrig = {fReader, "trigMatch_1LTrig"};
     TTreeReaderValue<Bool_t> trigMatch_1L2LTrigOR = {fReader, "trigMatch_1L2LTrigOR"};
     TTreeReaderValue<Bool_t> trigMatch_1LTrigOR = {fReader, "trigMatch_1LTrigOR"};
     TTreeReaderValue<Bool_t> trigMatch_2LTrig = {fReader, "trigMatch_2LTrig"};
     TTreeReaderValue<Bool_t> trigMatch_2LTrigOR = {fReader, "trigMatch_2LTrigOR"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_e24_lhmedium_L1EM20VH = {fReader, "trigMatch_HLT_e24_lhmedium_L1EM20VH"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_e60_lhmedium = {fReader, "trigMatch_HLT_e60_lhmedium"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_e120_lhloose = {fReader, "trigMatch_HLT_e120_lhloose"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_mu20_iloose_L1MU15 = {fReader, "trigMatch_HLT_mu20_iloose_L1MU15"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_2e12_lhloose_L12EM10VH = {fReader, "trigMatch_HLT_2e12_lhloose_L12EM10VH"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_mu18_mu8noL1 = {fReader, "trigMatch_HLT_mu18_mu8noL1"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_e17_lhloose_mu14 = {fReader, "trigMatch_HLT_e17_lhloose_mu14"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_e7_lhmedium_mu24 = {fReader, "trigMatch_HLT_e7_lhmedium_mu24"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_e24_lhtight_nod0_ivarloose = {fReader, "trigMatch_HLT_e24_lhtight_nod0_ivarloose"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_e24_lhmedium_nod0_L1EM20VH = {fReader, "trigMatch_HLT_e24_lhmedium_nod0_L1EM20VH"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_e60_medium = {fReader, "trigMatch_HLT_e60_medium"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_mu40 = {fReader, "trigMatch_HLT_mu40"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_mu24_iloose_L1MU15 = {fReader, "trigMatch_HLT_mu24_iloose_L1MU15"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_mu24_ivarloose_L1MU15 = {fReader, "trigMatch_HLT_mu24_ivarloose_L1MU15"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_mu24_ivarmedium = {fReader, "trigMatch_HLT_mu24_ivarmedium"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_mu24_imedium = {fReader, "trigMatch_HLT_mu24_imedium"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_mu26_imedium = {fReader, "trigMatch_HLT_mu26_imedium"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_2e15_lhvloose_nod0_L12EM13VH = {fReader, "trigMatch_HLT_2e15_lhvloose_nod0_L12EM13VH"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_2e17_lhvloose_nod0 = {fReader, "trigMatch_HLT_2e17_lhvloose_nod0"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_2mu10 = {fReader, "trigMatch_HLT_2mu10"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_2mu14 = {fReader, "trigMatch_HLT_2mu14"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_mu20_mu8noL1 = {fReader, "trigMatch_HLT_mu20_mu8noL1"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_mu22_mu8noL1 = {fReader, "trigMatch_HLT_mu22_mu8noL1"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1 = {fReader, "trigMatch_HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_e26_lhtight_nod0_ivarloose = {fReader, "trigMatch_HLT_e26_lhtight_nod0_ivarloose"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_e26_lhtight_nod0 = {fReader, "trigMatch_HLT_e26_lhtight_nod0"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_e60_lhmedium_nod0 = {fReader, "trigMatch_HLT_e60_lhmedium_nod0"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_e140_lhloose_nod0 = {fReader, "trigMatch_HLT_e140_lhloose_nod0"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_e300_etcut = {fReader, "trigMatch_HLT_e300_etcut"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_mu26_ivarmedium = {fReader, "trigMatch_HLT_mu26_ivarmedium"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_mu50 = {fReader, "trigMatch_HLT_mu50"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_mu60_0eta105_msonly = {fReader, "trigMatch_HLT_mu60_0eta105_msonly"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_2e17_lhvloose_nod0_L12EM15VHI = {fReader, "trigMatch_HLT_2e17_lhvloose_nod0_L12EM15VHI"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_2e24_lhvloose_nod0 = {fReader, "trigMatch_HLT_2e24_lhvloose_nod0"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_e17_lhloose_nod0_mu14 = {fReader, "trigMatch_HLT_e17_lhloose_nod0_mu14"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_e26_lhmedium_nod0_mu8noL1 = {fReader, "trigMatch_HLT_e26_lhmedium_nod0_mu8noL1"};
     TTreeReaderValue<Bool_t> trigMatch_HLT_e7_lhmedium_nod0_mu24 = {fReader, "trigMatch_HLT_e7_lhmedium_nod0_mu24"};

     TTreeReaderValue<Float_t> mu = {fReader, "mu"};
     TTreeReaderValue<Float_t> avg_mu = {fReader, "avg_mu"};
     TTreeReaderValue<Float_t> actual_mu = {fReader, "actual_mu"};
     TTreeReaderValue<Int_t> nVtx = {fReader, "nVtx"};
     TTreeReaderValue<Int_t> channel = {fReader, "channel"};
     TTreeReaderValue<Int_t> nLep_base = {fReader, "nLep_base"};
     TTreeReaderValue<Int_t> nLep_signal = {fReader, "nLep_signal"};
     TTreeReaderArray<int> lepFlavor = {fReader, "lepFlavor"};
     TTreeReaderArray<int> lepCharge = {fReader, "lepCharge"};
     TTreeReaderArray<int> lepAuthor = {fReader, "lepAuthor"};
     TTreeReaderArray<float> lepPt = {fReader, "lepPt"};
     TTreeReaderArray<float> lepEta = {fReader, "lepEta"};
     TTreeReaderArray<float> lepPhi = {fReader, "lepPhi"};
     TTreeReaderArray<float> lepM = {fReader, "lepM"};
     TTreeReaderArray<float> lepD0 = {fReader, "lepD0"};
     TTreeReaderArray<float> lepD0Sig = {fReader, "lepD0Sig"};
     TTreeReaderArray<float> lepZ0 = {fReader, "lepZ0"};
     TTreeReaderArray<float> lepZ0SinTheta = {fReader, "lepZ0SinTheta"};
     TTreeReaderValue<vector<bool>> lepPassOR = {fReader, "lepPassOR"};
     TTreeReaderArray<int> lepType = {fReader, "lepType"};
     TTreeReaderArray<int> lepOrigin = {fReader, "lepOrigin"};
     TTreeReaderArray<int> lepIFFClass = {fReader, "lepIFFClass"};
     TTreeReaderArray<int> lepEgMotherType = {fReader, "lepEgMotherType"};
     TTreeReaderArray<int> lepEgMotherOrigin = {fReader, "lepEgMotherOrigin"};
     TTreeReaderArray<int> lepEgMotherPdgId = {fReader, "lepEgMotherPdgId"};
     TTreeReaderArray<int> lepECIDS = {fReader, "lepECIDS"};
     TTreeReaderValue<vector<bool>> lepPassBL = {fReader, "lepPassBL"};
     TTreeReaderValue<vector<bool>> lepVeryLoose = {fReader, "lepVeryLoose"};
     TTreeReaderValue<vector<bool>> lepLoose = {fReader, "lepLoose"};
     TTreeReaderValue<vector<bool>> lepMedium = {fReader, "lepMedium"};
     TTreeReaderValue<vector<bool>> lepTight = {fReader, "lepTight"};
     TTreeReaderValue<vector<bool>> lepIsoFCHighPtCaloOnly = {fReader, "lepIsoFCHighPtCaloOnly"};
     TTreeReaderValue<vector<bool>> lepIsoFixedCutHighPtTrackOnly = {fReader, "lepIsoFixedCutHighPtTrackOnly"};
     TTreeReaderValue<vector<bool>> lepIsoGradient = {fReader, "lepIsoGradient"};
     TTreeReaderValue<vector<bool>> lepIsoFCLoose = {fReader, "lepIsoFCLoose"};
     TTreeReaderValue<vector<bool>> lepIsoFCTight = {fReader, "lepIsoFCTight"};
     TTreeReaderValue<vector<bool>> lepIsoFCTightTrackOnly = {fReader, "lepIsoFCTightTrackOnly"};
     TTreeReaderValue<vector<bool>> lepTruthMatched = {fReader, "lepTruthMatched"};
     TTreeReaderArray<int> lepTruthCharge = {fReader, "lepTruthCharge"};
     TTreeReaderArray<float> lepTruthPt = {fReader, "lepTruthPt"};
     TTreeReaderArray<float> lepTruthEta = {fReader, "lepTruthEta"};
     TTreeReaderArray<float> lepTruthPhi = {fReader, "lepTruthPhi"};
     TTreeReaderArray<float> lepTruthM = {fReader, "lepTruthM"};
     TTreeReaderValue<Int_t> nJet30 = {fReader, "nJet30"};
     TTreeReaderValue<Int_t> nJet20 = {fReader, "nJet20"};
     TTreeReaderValue<Int_t> nBJet20_MV2c10_FixedCutBEff_77 = {fReader, "nBJet20_MV2c10_FixedCutBEff_77"};
     TTreeReaderArray<float> jetPt = {fReader, "jetPt"};
     TTreeReaderArray<float> jetEta = {fReader, "jetEta"};
     TTreeReaderArray<float> jetPhi = {fReader, "jetPhi"};
     TTreeReaderArray<float> jetM = {fReader, "jetM"};
     TTreeReaderArray<float> jetJVT = {fReader, "jetJVT"};
     TTreeReaderValue<vector<bool>> jetPassOR = {fReader, "jetPassOR"};
     TTreeReaderValue<vector<bool>> jetSignal = {fReader, "jetSignal"};
     TTreeReaderValue<Float_t> mjj = {fReader, "mjj"};
     TTreeReaderArray<float> jetTileEnergy = {fReader, "jetTileEnergy"};
     TTreeReaderArray<float> jetMV2c10 = {fReader, "jetMV2c10"};
     TTreeReaderArray<float> jetdl1 = {fReader, "jetdl1"};
     TTreeReaderValue<Float_t> met_Et = {fReader, "met_Et"};
     TTreeReaderValue<Float_t> met_Sign = {fReader, "met_Sign"};
     TTreeReaderValue<Float_t> met_Phi = {fReader, "met_Phi"};
     TTreeReaderValue<Float_t> mll = {fReader, "mll"};
     TTreeReaderValue<Double_t> pileupWeight = {fReader, "pileupWeight"};
     TTreeReaderValue<Double_t> leptonWeight = {fReader, "leptonWeight"};
     TTreeReaderValue<Double_t> eventWeight = {fReader, "eventWeight"};
     TTreeReaderValue<Double_t> bTagWeight = {fReader, "bTagWeight"};
     TTreeReaderValue<Double_t> jvtWeight = {fReader, "jvtWeight"};
     TTreeReaderValue<Double_t> globalDiLepTrigSF = {fReader, "globalDiLepTrigSF"};
     TTreeReaderValue<Double_t> flavSymWeight = {fReader, "flavSymWeight"};
     TTreeReaderValue<Double_t> genWeight = {fReader, "genWeight"};
     TTreeReaderValue<Double_t> genWeightUp = {fReader, "genWeightUp"};
     TTreeReaderValue<Double_t> genWeightDown = {fReader, "genWeightDown"};
  */
  /**
     TTreeReaderValue<Float_t> leptonWeight_EL_CHARGEID_STAT__1down = {fReader, "leptonWeight_EL_CHARGEID_STAT__1down"};
     TTreeReaderValue<Float_t> leptonWeight_EL_CHARGEID_STAT__1up = {fReader, "leptonWeight_EL_CHARGEID_STAT__1up"};
     TTreeReaderValue<Float_t> leptonWeight_EL_CHARGEID_SYStotal__1down = {fReader, "leptonWeight_EL_CHARGEID_SYStotal__1down"};
     TTreeReaderValue<Float_t> leptonWeight_EL_CHARGEID_SYStotal__1up = {fReader, "leptonWeight_EL_CHARGEID_SYStotal__1up"};
     TTreeReaderValue<Float_t> leptonWeight_EL_EFF_ChargeIDSel_TOTAL_1NPCOR_PLUS_UNCOR__1down = {fReader, "leptonWeight_EL_EFF_ChargeIDSel_TOTAL_1NPCOR_PLUS_UNCOR__1down"};
     TTreeReaderValue<Float_t> leptonWeight_EL_EFF_ChargeIDSel_TOTAL_1NPCOR_PLUS_UNCOR__1up = {fReader, "leptonWeight_EL_EFF_ChargeIDSel_TOTAL_1NPCOR_PLUS_UNCOR__1up"};
     TTreeReaderValue<Float_t> leptonWeight_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR__1down = {fReader, "leptonWeight_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR__1down"};
     TTreeReaderValue<Float_t> leptonWeight_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR__1up = {fReader, "leptonWeight_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR__1up"};
     TTreeReaderValue<Float_t> leptonWeight_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR__1down = {fReader, "leptonWeight_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR__1down"};
     TTreeReaderValue<Float_t> leptonWeight_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR__1up = {fReader, "leptonWeight_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR__1up"};
     TTreeReaderValue<Float_t> leptonWeight_EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR__1down = {fReader, "leptonWeight_EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR__1down"};
     TTreeReaderValue<Float_t> leptonWeight_EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR__1up = {fReader, "leptonWeight_EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR__1up"};
     TTreeReaderValue<Float_t> leptonWeight_MUON_EFF_ISO_STAT__1down = {fReader, "leptonWeight_MUON_EFF_ISO_STAT__1down"};
     TTreeReaderValue<Float_t> leptonWeight_MUON_EFF_ISO_STAT__1up = {fReader, "leptonWeight_MUON_EFF_ISO_STAT__1up"};
     TTreeReaderValue<Float_t> leptonWeight_MUON_EFF_ISO_SYS__1down = {fReader, "leptonWeight_MUON_EFF_ISO_SYS__1down"};
     TTreeReaderValue<Float_t> leptonWeight_MUON_EFF_ISO_SYS__1up = {fReader, "leptonWeight_MUON_EFF_ISO_SYS__1up"};
     TTreeReaderValue<Float_t> leptonWeight_MUON_EFF_RECO_STAT__1down = {fReader, "leptonWeight_MUON_EFF_RECO_STAT__1down"};
     TTreeReaderValue<Float_t> leptonWeight_MUON_EFF_RECO_STAT__1up = {fReader, "leptonWeight_MUON_EFF_RECO_STAT__1up"};
     TTreeReaderValue<Float_t> leptonWeight_MUON_EFF_RECO_STAT_LOWPT__1down = {fReader, "leptonWeight_MUON_EFF_RECO_STAT_LOWPT__1down"};
     TTreeReaderValue<Float_t> leptonWeight_MUON_EFF_RECO_STAT_LOWPT__1up = {fReader, "leptonWeight_MUON_EFF_RECO_STAT_LOWPT__1up"};
     TTreeReaderValue<Float_t> leptonWeight_MUON_EFF_RECO_SYS__1down = {fReader, "leptonWeight_MUON_EFF_RECO_SYS__1down"};
     TTreeReaderValue<Float_t> leptonWeight_MUON_EFF_RECO_SYS__1up = {fReader, "leptonWeight_MUON_EFF_RECO_SYS__1up"};
     TTreeReaderValue<Float_t> leptonWeight_MUON_EFF_RECO_SYS_LOWPT__1down = {fReader, "leptonWeight_MUON_EFF_RECO_SYS_LOWPT__1down"};
     TTreeReaderValue<Float_t> leptonWeight_MUON_EFF_RECO_SYS_LOWPT__1up = {fReader, "leptonWeight_MUON_EFF_RECO_SYS_LOWPT__1up"};
     TTreeReaderValue<Float_t> leptonWeight_MUON_EFF_TTVA_STAT__1down = {fReader, "leptonWeight_MUON_EFF_TTVA_STAT__1down"};
     TTreeReaderValue<Float_t> leptonWeight_MUON_EFF_TTVA_STAT__1up = {fReader, "leptonWeight_MUON_EFF_TTVA_STAT__1up"};
     TTreeReaderValue<Float_t> leptonWeight_MUON_EFF_TTVA_SYS__1down = {fReader, "leptonWeight_MUON_EFF_TTVA_SYS__1down"};
     TTreeReaderValue<Float_t> leptonWeight_MUON_EFF_TTVA_SYS__1up = {fReader, "leptonWeight_MUON_EFF_TTVA_SYS__1up"};
     TTreeReaderValue<Float_t> trigWeight_EL_EFF_TriggerEff_TOTAL_1NPCOR_PLUS_UNCOR__1down = {fReader, "trigWeight_EL_EFF_TriggerEff_TOTAL_1NPCOR_PLUS_UNCOR__1down"};
     TTreeReaderValue<Float_t> trigWeight_EL_EFF_TriggerEff_TOTAL_1NPCOR_PLUS_UNCOR__1up = {fReader, "trigWeight_EL_EFF_TriggerEff_TOTAL_1NPCOR_PLUS_UNCOR__1up"};
     TTreeReaderValue<Float_t> trigWeight_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1down = {fReader, "trigWeight_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1down"};
     TTreeReaderValue<Float_t> trigWeight_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1up = {fReader, "trigWeight_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1up"};
     TTreeReaderValue<Float_t> trigWeight_MUON_EFF_TrigStatUncertainty__1down = {fReader, "trigWeight_MUON_EFF_TrigStatUncertainty__1down"};
     TTreeReaderValue<Float_t> trigWeight_MUON_EFF_TrigStatUncertainty__1up = {fReader, "trigWeight_MUON_EFF_TrigStatUncertainty__1up"};
     TTreeReaderValue<Float_t> trigWeight_MUON_EFF_TrigSystUncertainty__1down = {fReader, "trigWeight_MUON_EFF_TrigSystUncertainty__1down"};
     TTreeReaderValue<Float_t> trigWeight_MUON_EFF_TrigSystUncertainty__1up = {fReader, "trigWeight_MUON_EFF_TrigSystUncertainty__1up"};
     TTreeReaderValue<Float_t> bTagWeight_FT_EFF_B_systematics__1down = {fReader, "bTagWeight_FT_EFF_B_systematics__1down"};
     TTreeReaderValue<Float_t> bTagWeight_FT_EFF_B_systematics__1up = {fReader, "bTagWeight_FT_EFF_B_systematics__1up"};
     TTreeReaderValue<Float_t> bTagWeight_FT_EFF_C_systematics__1down = {fReader, "bTagWeight_FT_EFF_C_systematics__1down"};
     TTreeReaderValue<Float_t> bTagWeight_FT_EFF_C_systematics__1up = {fReader, "bTagWeight_FT_EFF_C_systematics__1up"};
     TTreeReaderValue<Float_t> bTagWeight_FT_EFF_Light_systematics__1down = {fReader, "bTagWeight_FT_EFF_Light_systematics__1down"};
     TTreeReaderValue<Float_t> bTagWeight_FT_EFF_Light_systematics__1up = {fReader, "bTagWeight_FT_EFF_Light_systematics__1up"};
     TTreeReaderValue<Float_t> bTagWeight_FT_EFF_extrapolation__1down = {fReader, "bTagWeight_FT_EFF_extrapolation__1down"};
     TTreeReaderValue<Float_t> bTagWeight_FT_EFF_extrapolation__1up = {fReader, "bTagWeight_FT_EFF_extrapolation__1up"};
     TTreeReaderValue<Float_t> bTagWeight_FT_EFF_extrapolation_from_charm__1down = {fReader, "bTagWeight_FT_EFF_extrapolation_from_charm__1down"};
     TTreeReaderValue<Float_t> bTagWeight_FT_EFF_extrapolation_from_charm__1up = {fReader, "bTagWeight_FT_EFF_extrapolation_from_charm__1up"};
     TTreeReaderValue<Float_t> jvtWeight_JET_JvtEfficiency__1down = {fReader, "jvtWeight_JET_JvtEfficiency__1down"};
     TTreeReaderValue<Float_t> jvtWeight_JET_JvtEfficiency__1up = {fReader, "jvtWeight_JET_JvtEfficiency__1up"};
     TTreeReaderValue<Float_t> jvtWeight_JET_fJvtEfficiency__1down = {fReader, "jvtWeight_JET_fJvtEfficiency__1down"};
     TTreeReaderValue<Float_t> jvtWeight_JET_fJvtEfficiency__1up = {fReader, "jvtWeight_JET_fJvtEfficiency__1up"};
     TTreeReaderValue<Float_t> pileupWeightUp = {fReader, "pileupWeightUp"};
     TTreeReaderValue<Float_t> pileupWeightDown = {fReader, "pileupWeightDown"};
  */
  
  // TTreeReaderValue<ULong64_t> PRWHash = {fReader, "PRWHash"};
  // TTreeReaderValue<ULong64_t> EventNumber = {fReader, "EventNumber"};
  // TTreeReaderValue<Float_t> xsec = {fReader, "xsec"};
  // TTreeReaderValue<Float_t> GenHt = {fReader, "GenHt"};
  // TTreeReaderValue<Float_t> GenMET = {fReader, "GenMET"};
  // TTreeReaderValue<Int_t> DatasetNumber = {fReader, "DatasetNumber"};
  // TTreeReaderValue<Int_t> RunNumber = {fReader, "RunNumber"};
  // TTreeReaderValue<Int_t> RandomRunNumber = {fReader, "RandomRunNumber"};
  // TTreeReaderValue<Int_t> FS = {fReader, "FS"};
  
  
  
  // 2L2j NTUPLES
  // Readers to access the data (delete the ones you do not need).
  /**
     TTreeReaderValue<Double_t> trigWeight_2LTrig = {fReader, "trigWeight_2LTrig"};
     TTreeReaderValue<Bool_t> trigMatch_2LTrig = {fReader, "trigMatch_2LTrig"};
     TTreeReaderValue<Float_t> mu = {fReader, "mu"};
     TTreeReaderValue<Float_t> avg_mu = {fReader, "avg_mu"};
     TTreeReaderValue<Float_t> actual_mu = {fReader, "actual_mu"};
     TTreeReaderValue<Int_t> nVtx = {fReader, "nVtx"};
     TTreeReaderValue<Int_t> channel = {fReader, "channel"};
     TTreeReaderValue<Int_t> nLep_combi = {fReader, "nLep_combi"};
     TTreeReaderValue<Int_t> nLep_base = {fReader, "nLep_base"};
     TTreeReaderValue<Int_t> nLep_signal = {fReader, "nLep_signal"};
     TTreeReaderArray<int> lepFlavor = {fReader, "lepFlavor"};
     TTreeReaderArray<int> lepCharge = {fReader, "lepCharge"};
     TTreeReaderArray<int> lepAuthor = {fReader, "lepAuthor"};
     TTreeReaderArray<float> lepPt = {fReader, "lepPt"};
     TTreeReaderArray<float> lepEta = {fReader, "lepEta"};
     TTreeReaderArray<float> lepPhi = {fReader, "lepPhi"};
     TTreeReaderArray<float> lepM = {fReader, "lepM"};
     TTreeReaderArray<float> lepD0 = {fReader, "lepD0"};
     TTreeReaderArray<float> lepD0Sig = {fReader, "lepD0Sig"};
     TTreeReaderArray<float> lepZ0 = {fReader, "lepZ0"};
     TTreeReaderArray<float> lepZ0SinTheta = {fReader, "lepZ0SinTheta"};
     TTreeReaderValue<vector<bool>> lepPassOR = {fReader, "lepPassOR"};
     TTreeReaderArray<int> lepType = {fReader, "lepType"};
     TTreeReaderArray<int> lepOrigin = {fReader, "lepOrigin"};
     TTreeReaderArray<int> lepIFFClass = {fReader, "lepIFFClass"};
     TTreeReaderArray<int> lepEgMotherType = {fReader, "lepEgMotherType"};
     TTreeReaderArray<int> lepEgMotherOrigin = {fReader, "lepEgMotherOrigin"};
     TTreeReaderArray<int> lepEgMotherPdgId = {fReader, "lepEgMotherPdgId"};
     TTreeReaderArray<int> lepECIDS = {fReader, "lepECIDS"};
     TTreeReaderArray<int> lepNPix = {fReader, "lepNPix"};
     TTreeReaderValue<vector<bool>> lepPassBL = {fReader, "lepPassBL"};
     TTreeReaderValue<vector<bool>> lepSignal = {fReader, "lepSignal"};
     TTreeReaderValue<vector<bool>> lepTruthMatched = {fReader, "lepTruthMatched"};
     TTreeReaderArray<int> lepTruthCharge = {fReader, "lepTruthCharge"};
     TTreeReaderArray<float> lepTruthPt = {fReader, "lepTruthPt"};
     TTreeReaderArray<float> lepTruthEta = {fReader, "lepTruthEta"};
     TTreeReaderArray<float> lepTruthPhi = {fReader, "lepTruthPhi"};
     TTreeReaderArray<float> lepTruthM = {fReader, "lepTruthM"};
     TTreeReaderValue<Int_t> nJet30 = {fReader, "nJet30"};
     TTreeReaderValue<Int_t> nJet20 = {fReader, "nJet20"};
     TTreeReaderValue<Int_t> nBJet20_MV2c10_FixedCutBEff_77 = {fReader, "nBJet20_MV2c10_FixedCutBEff_77"};
     TTreeReaderArray<float> jetPt = {fReader, "jetPt"};
     TTreeReaderArray<float> jetEta = {fReader, "jetEta"};
     TTreeReaderArray<float> jetPhi = {fReader, "jetPhi"};
     TTreeReaderArray<float> jetM = {fReader, "jetM"};
     TTreeReaderValue<Float_t> mjj = {fReader, "mjj"};
     TTreeReaderValue<Float_t> mjj_minDPhiZMET = {fReader, "mjj_minDPhiZMET"};
     TTreeReaderValue<Float_t> Rjj = {fReader, "Rjj"};
     TTreeReaderValue<Float_t> Rjj_minDPhiZMET = {fReader, "Rjj_minDPhiZMET"};
     TTreeReaderArray<float> jetTileEnergy = {fReader, "jetTileEnergy"};
     TTreeReaderArray<float> jetMV2c10 = {fReader, "jetMV2c10"};
     TTreeReaderValue<Float_t> vectorSumJetsPt = {fReader, "vectorSumJetsPt"};
     TTreeReaderValue<Float_t> vectorSumJetsEta = {fReader, "vectorSumJetsEta"};
     TTreeReaderValue<Float_t> vectorSumJetsPhi = {fReader, "vectorSumJetsPhi"};
     TTreeReaderValue<Float_t> vectorSumJetsM = {fReader, "vectorSumJetsM"};
     TTreeReaderValue<Float_t> met_Et = {fReader, "met_Et"};
     TTreeReaderValue<Float_t> met_Sign = {fReader, "met_Sign"};
     TTreeReaderValue<Float_t> met_Phi = {fReader, "met_Phi"};
     TTreeReaderValue<Float_t> TST_Et = {fReader, "TST_Et"};
     TTreeReaderValue<Float_t> TST_Phi = {fReader, "TST_Phi"};
     TTreeReaderValue<Float_t> deltaPhi_MET_TST_Phi = {fReader, "deltaPhi_MET_TST_Phi"};
     TTreeReaderValue<Float_t> Ht30 = {fReader, "Ht30"};
     TTreeReaderValue<Float_t> mbb = {fReader, "mbb"};
     TTreeReaderValue<Float_t> Rbb = {fReader, "Rbb"};
     TTreeReaderValue<Float_t> METOverPtZ = {fReader, "METOverPtZ"};
     TTreeReaderValue<Float_t> METOverPtW = {fReader, "METOverPtW"};
     TTreeReaderValue<Float_t> PtISR = {fReader, "PtISR"};
     TTreeReaderValue<Float_t> METOverPtISR = {fReader, "METOverPtISR"};
     TTreeReaderValue<Float_t> METOverHT = {fReader, "METOverHT"};
     TTreeReaderValue<Float_t> METOverJ1pT = {fReader, "METOverJ1pT"};
     TTreeReaderValue<Float_t> DPhiJ1Met = {fReader, "DPhiJ1Met"};
     TTreeReaderValue<Float_t> DPhiJ2Met = {fReader, "DPhiJ2Met"};
     TTreeReaderValue<Float_t> DPhiJ3Met = {fReader, "DPhiJ3Met"};
     TTreeReaderValue<Float_t> DPhiJ4Met = {fReader, "DPhiJ4Met"};
     TTreeReaderValue<Float_t> minDPhi2JetsMet = {fReader, "minDPhi2JetsMet"};
     TTreeReaderValue<Float_t> minDPhi4JetsMet = {fReader, "minDPhi4JetsMet"};
     TTreeReaderValue<Float_t> minDPhiAllJetsMet = {fReader, "minDPhiAllJetsMet"};
     TTreeReaderValue<Float_t> dPhiPjjMet = {fReader, "dPhiPjjMet"};
     TTreeReaderValue<Float_t> dPhiPjjMet_minDPhiZMET = {fReader, "dPhiPjjMet_minDPhiZMET"};
     TTreeReaderValue<Float_t> dPhiMetISR = {fReader, "dPhiMetISR"};
     TTreeReaderValue<Float_t> dPhiMetJet1 = {fReader, "dPhiMetJet1"};
     TTreeReaderValue<Float_t> METOverHTLep = {fReader, "METOverHTLep"};
     TTreeReaderValue<Float_t> mll = {fReader, "mll"};
     TTreeReaderValue<Float_t> mlll = {fReader, "mlll"};
     TTreeReaderValue<Float_t> mllll = {fReader, "mllll"};
     TTreeReaderValue<Float_t> Rll = {fReader, "Rll"};
     TTreeReaderValue<Float_t> Ptll = {fReader, "Ptll"};
     TTreeReaderValue<Float_t> absEtall = {fReader, "absEtall"};
     TTreeReaderValue<Float_t> dPhiPllMet = {fReader, "dPhiPllMet"};
     TTreeReaderValue<Float_t> dPhill = {fReader, "dPhill"};
     TTreeReaderValue<Float_t> mt2leplsp_0 = {fReader, "mt2leplsp_0"};
     TTreeReaderValue<Double_t> pileupWeight = {fReader, "pileupWeight"};
     TTreeReaderValue<Double_t> leptonWeight = {fReader, "leptonWeight"};
     TTreeReaderValue<Double_t> eventWeight = {fReader, "eventWeight"};
     TTreeReaderValue<Double_t> bTagWeight = {fReader, "bTagWeight"};
     TTreeReaderValue<Double_t> jvtWeight = {fReader, "jvtWeight"};
     TTreeReaderValue<Double_t> globalDiLepTrigSF = {fReader, "globalDiLepTrigSF"};
     TTreeReaderValue<Double_t> genWeight = {fReader, "genWeight"};
     TTreeReaderValue<Double_t> genWeightUp = {fReader, "genWeightUp"};
     TTreeReaderValue<Double_t> genWeightDown = {fReader, "genWeightDown"};
     TTreeReaderValue<Double_t> truthMll = {fReader, "truthMll"};
     TTreeReaderValue<Double_t> winoBinoMllWeight = {fReader, "winoBinoMllWeight"};
     TTreeReaderValue<Double_t> winoBinoXsecWeight = {fReader, "winoBinoXsecWeight"};
     TTreeReaderValue<Double_t> winoBinoBrFracWeight = {fReader, "winoBinoBrFracWeight"};
     TTreeReaderValue<Double_t> winoBinoWeight = {fReader, "winoBinoWeight"};
     TTreeReaderValue<Double_t> ttbarNNLOWeight = {fReader, "ttbarNNLOWeight"};
     TTreeReaderValue<Double_t> ttbarNNLOWeightUp = {fReader, "ttbarNNLOWeightUp"};
     TTreeReaderValue<Double_t> ttbarNNLOWeightDown = {fReader, "ttbarNNLOWeightDown"};
     TTreeReaderValue<Float_t> truthTopPt = {fReader, "truthTopPt"};
     TTreeReaderValue<Float_t> truthAntiTopPt = {fReader, "truthAntiTopPt"};
     TTreeReaderValue<Float_t> truthTtbarPt = {fReader, "truthTtbarPt"};
     TTreeReaderValue<Float_t> truthTtbarM = {fReader, "truthTtbarM"};
     TTreeReaderValue<Float_t> x1 = {fReader, "x1"};
     TTreeReaderValue<Float_t> x2 = {fReader, "x2"};
     TTreeReaderValue<Float_t> pdf1 = {fReader, "pdf1"};
     TTreeReaderValue<Float_t> pdf2 = {fReader, "pdf2"};
     TTreeReaderValue<Float_t> scalePDF = {fReader, "scalePDF"};
     TTreeReaderValue<Int_t> id1 = {fReader, "id1"};
     TTreeReaderValue<Int_t> id2 = {fReader, "id2"};
     TTreeReaderValue<Int_t> nLeps_RJ = {fReader, "nLeps_RJ"};
     TTreeReaderValue<Int_t> nJets_RJ = {fReader, "nJets_RJ"};
     TTreeReaderValue<Int_t> nBtagJets_RJ = {fReader, "nBtagJets_RJ"};
     TTreeReaderArray<float> jetPt_RJ = {fReader, "jetPt_RJ"};
     TTreeReaderArray<float> jetEta_RJ = {fReader, "jetEta_RJ"};
     TTreeReaderArray<float> jetPhi_RJ = {fReader, "jetPhi_RJ"};
     TTreeReaderArray<float> jetM_RJ = {fReader, "jetM_RJ"};
     TTreeReaderArray<float> lepPt_RJ = {fReader, "lepPt_RJ"};
     TTreeReaderArray<float> lepEta_RJ = {fReader, "lepEta_RJ"};
     TTreeReaderArray<float> lepPhi_RJ = {fReader, "lepPhi_RJ"};
     TTreeReaderArray<float> lepE_RJ = {fReader, "lepE_RJ"};
     TTreeReaderArray<float> lepsign_RJ = {fReader, "lepsign_RJ"};
     TTreeReaderValue<Bool_t> is2Lep2Jet = {fReader, "is2Lep2Jet"};
     TTreeReaderValue<Bool_t> is2L2JInt = {fReader, "is2L2JInt"};
     TTreeReaderValue<Bool_t> is3Lep = {fReader, "is3Lep"};
     TTreeReaderValue<Bool_t> is3LInt = {fReader, "is3LInt"};
     TTreeReaderValue<Bool_t> is3Lep2Jet = {fReader, "is3Lep2Jet"};
     TTreeReaderValue<Bool_t> is3Lep3Jet = {fReader, "is3Lep3Jet"};
     TTreeReaderValue<Bool_t> is4Lep2Jet = {fReader, "is4Lep2Jet"};
     TTreeReaderValue<Bool_t> is4Lep3Jet = {fReader, "is4Lep3Jet"};
     TTreeReaderValue<Float_t> mll_RJ = {fReader, "mll_RJ"};
     TTreeReaderValue<Float_t> H2PP = {fReader, "H2PP"};
     TTreeReaderValue<Float_t> H5PP = {fReader, "H5PP"};
     TTreeReaderValue<Float_t> RPT_HT5PP = {fReader, "RPT_HT5PP"};
     TTreeReaderValue<Float_t> R_minH2P_minH3P = {fReader, "R_minH2P_minH3P"};
     TTreeReaderValue<Float_t> R_H2PP_H5PP = {fReader, "R_H2PP_H5PP"};
     TTreeReaderValue<Float_t> dphiVP = {fReader, "dphiVP"};
     TTreeReaderValue<Float_t> minDphi = {fReader, "minDphi"};
     TTreeReaderValue<Float_t> mTW = {fReader, "mTW"};
     TTreeReaderValue<Float_t> H4PP = {fReader, "H4PP"};
     TTreeReaderValue<Float_t> RPT_HT4PP = {fReader, "RPT_HT4PP"};
     TTreeReaderValue<Float_t> R_HT4PP_H4PP = {fReader, "R_HT4PP_H4PP"};
     TTreeReaderValue<Float_t> PTISR = {fReader, "PTISR"};
     TTreeReaderValue<Float_t> RISR = {fReader, "RISR"};
     TTreeReaderValue<Float_t> PTI = {fReader, "PTI"};
     TTreeReaderValue<Float_t> dphiISRI = {fReader, "dphiISRI"};
     TTreeReaderValue<Float_t> PTCM = {fReader, "PTCM"};
     TTreeReaderValue<Int_t> NjS = {fReader, "NjS"};
     TTreeReaderValue<Int_t> NjISR = {fReader, "NjISR"};
     TTreeReaderValue<Float_t> MZ = {fReader, "MZ"};
     TTreeReaderValue<Float_t> MJ = {fReader, "MJ"};
     TTreeReaderValue<Float_t> mTl3 = {fReader, "mTl3"};
     TTreeReaderValue<Float_t> lept1Pt_VR = {fReader, "lept1Pt_VR"};
     TTreeReaderValue<Float_t> lept1sign_VR = {fReader, "lept1sign_VR"};
     TTreeReaderValue<Float_t> lept2Pt_VR = {fReader, "lept2Pt_VR"};
     TTreeReaderValue<Float_t> lept2sign_VR = {fReader, "lept2sign_VR"};
     TTreeReaderValue<Float_t> mll_RJ_VR = {fReader, "mll_RJ_VR"};
     TTreeReaderValue<Float_t> H2PP_VR = {fReader, "H2PP_VR"};
     TTreeReaderValue<Float_t> H5PP_VR = {fReader, "H5PP_VR"};
     TTreeReaderValue<Float_t> RPT_HT5PP_VR = {fReader, "RPT_HT5PP_VR"};
     TTreeReaderValue<Float_t> R_minH2P_minH3P_VR = {fReader, "R_minH2P_minH3P_VR"};
     TTreeReaderValue<Float_t> R_H2PP_H5PP_VR = {fReader, "R_H2PP_H5PP_VR"};
     TTreeReaderValue<Float_t> dphiVP_VR = {fReader, "dphiVP_VR"};
     TTreeReaderValue<Float_t> PTISR_VR = {fReader, "PTISR_VR"};
     TTreeReaderValue<Float_t> RISR_VR = {fReader, "RISR_VR"};
     TTreeReaderValue<Float_t> PTI_VR = {fReader, "PTI_VR"};
     TTreeReaderValue<Float_t> dphiISRI_VR = {fReader, "dphiISRI_VR"};
     TTreeReaderValue<Float_t> PTCM_VR = {fReader, "PTCM_VR"};
     TTreeReaderValue<Int_t> NjS_VR = {fReader, "NjS_VR"};
     TTreeReaderValue<Int_t> NjISR_VR = {fReader, "NjISR_VR"};
     TTreeReaderValue<Float_t> MZ_VR = {fReader, "MZ_VR"};
     TTreeReaderValue<Float_t> MJ_VR = {fReader, "MJ_VR"};
     TTreeReaderValue<ULong64_t> PRWHash = {fReader, "PRWHash"};
     TTreeReaderValue<ULong64_t> EventNumber = {fReader, "EventNumber"};
     TTreeReaderValue<Float_t> xsec = {fReader, "xsec"};
     TTreeReaderValue<Float_t> GenHt = {fReader, "GenHt"};
     TTreeReaderValue<Float_t> GenMET = {fReader, "GenMET"};
     TTreeReaderValue<Int_t> DatasetNumber = {fReader, "DatasetNumber"};
     TTreeReaderValue<Int_t> RunNumber = {fReader, "RunNumber"};
     TTreeReaderValue<Int_t> RandomRunNumber = {fReader, "RandomRunNumber"};
     TTreeReaderValue<Int_t> FS = {fReader, "FS"};
  */  
  
  MySusySkimAnalysis(TTree * /*tree*/ =0) { }
  virtual ~MySusySkimAnalysis() { }
  virtual Int_t   Version() const { return 2; }
  virtual void    Begin(TTree *tree);
  virtual void    SlaveBegin(TTree *tree);
  virtual void    Init(TTree *tree);
  virtual Bool_t  Notify();
  virtual Bool_t  Process(Long64_t entry);
  virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
  virtual void    SetOption(const char *option) { fOption = option; }
  virtual void    SetObject(TObject *obj) { fObject = obj; }
  virtual void    SetInputList(TList *input) { fInput = input; }
  virtual TList  *GetOutputList() const { return fOutput; }
  virtual void    SlaveTerminate();
  virtual void    Terminate();
  std::vector<TString> getTriggerCat(int id, int yr);

  

   
  //  truthENUM TruthClassification(int T, int O, int firstEgMotherT, int firstEgMotherO, int firstEgMotherPdgId, int lepCharge, int lepFlav);
  void WriteToFile();
  std::vector<TString> getFillVector(int idx1, int idx2, int idx0);
  bool isBLjet(int idx);
  bool isBjet(int idx);
  bool isSGjet(int idx);
  bool isL(int idx, bool checkOR = false);
  bool isT(int idx, bool checkOR = false);
  void isConvRegion(int co_mu_1, int co_mu_2, int co_el_1, int n_bjet, TString LT_postfix);
  void isHFRegion(std::vector<int> HF_probe_lep, std::vector<int> HF_tag_mu, int n_bjet, TString LT_postfix, bool isTightTag = false);
  void isLFRegion(std::vector<int> HF_probe_el, std::vector<int> HF_tag_el, TString LT_postfix);
  void fillLandT(int idx, std::vector<TString> regions, TString ananame, bool isFake = true, bool trigger = true, TString tname = "");
  void isREALRegion(int idx0, int idx1, int idx2, TString LT_postfix, bool isNOOR = false);
  std::vector<double> getFakeWeight(int idx1, int idx2, int idx0, TString year, TString LT_postfix, int vb);
  void classifyEvent(int idx1, int idx2, int& isTT, int& isTL, int& isLT, int& isLL);
  std::pair <TString,TString> getFillString(int idx1, int idx2, int idx0);
  bool ProcessCut(Long64_t entry);
  void fillFakeRateHist(TString key, int idx, bool lepIsTight);
  std::vector<double> getFakeWeight1D(int idx1, int idx2, int idx0, TString year, TString LT_postfix, int vb);
  void fill2L2JTree(int idx1, int idx2, std::map<TString, double> fake_weight);
  void fillMYTree(int idx1, int idx2, std::map< TString, std::vector<double> > fake_wgt_vec);
  void DefineCutFlow();
  void printCutflow();
  void checkTrigMatch(TString trigname, int idx, TString key, bool lepIsTight, TString trigtype = "L");
  void fillTruthInfo(TString& key, int idx, bool lepistight, double wgt);
  Float_t getPtThresholdTrigger(TString tname, TString tobj, int vb = 0);
  std::vector<TString> checkTriggerMatch(int id, int yr);
  void printInfo();

  ClassDef(MySusySkimAnalysis,0);

};

#endif

#ifdef MySusySkimAnalysis_cxx
void MySusySkimAnalysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
   //fChain = (TChain*)fReader.GetTree();
}

Bool_t MySusySkimAnalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef MySusySkimAnalysis_cxx
