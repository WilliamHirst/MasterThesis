//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Oct 23 13:02:19 2020 by ROOT version 6.18/00
// from TChain myTree/
//////////////////////////////////////////////////////////

#ifndef myTree_h
#define myTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include "MatrixMethod.h"
#include "MatrixMethod.cxx"
#include "GetFakeWeight.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/RInterface.hxx"

// Headers needed by this particular selector
#include <vector>



class myTree : public TSelector {
public :
  TTreeReader     fReader;  //!the tree reader
  TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

  int nentries;
  int n_entry;
  GetFakeWeight *gfw;
  MatrixMethod *mtx;
  std::vector<TString> looseTightDef;
  vector<TString> uncert_keys;
  TString filename_prev;
  TString filename;

  int fnpidx;
  std::map<TString,Double_t> BDTweight;
  std::vector<std::string> bdt_vec;
  std::map<Int_t,std::vector<Long64_t>> evnums;

  

  //ROOT::RDataFrame rdf;

  
  TH1F *h_costhetastar;
  TH1F *h_triggermatched_mu;
  TH1F *h_triggermatched_el;
  TH2F *h_triggermatched_el1_el2;
  TH2F *h_triggermatched_mu1_mu2;
  TH2F *h_triggermatched_el_mu;
  TH1F *h_fake_weight;
  TH1F *h_cutflow;
  TH2F *h_lep_count;
  TH2F *h_lep_tm;

  TH2F *h_r1_f1_TT;
  TH2F *h_r1_f1_Tl;
  TH2F *h_r1_f1_lT;
  TH2F *h_r1_f1_ll;

  TH2F *h_r2_f2_TT;
  TH2F *h_r2_f2_Tl;
  TH2F *h_r2_f2_lT;
  TH2F *h_r2_f2_ll;

  TH1F *h_eta_smalldiff;
  TH1F *h_tmatch_smalldiff;

  th1Map h_real1_diff;
  th1Map h_real2_diff;

  bool isDF_1bjet;
  bool isCR_Top;
  std::map<TString,std::vector<double>> BDTscoremap;
  std::map<Int_t,std::map<TString,std::vector<double>>> BDTscores;

  TH2F *h_lep_fake_trig_cat;

  map< TString, TH1F* > h_lep_BDT_n;
  map< TString, TH2F* > h_lep_BDT_r;
  TH1F *h_preselection_met_OS;
  TH1F *h_preselection_met_OS_FNP;
  TH1F *h_preselection_lep2pt_EE_OS;
  TH1F *h_preselection_lep2pt_EE_OS_FNP;
  TH1F *h_preselection_lep2pt_EM_OS;
  TH1F *h_preselection_lep2pt_EM_OS_FNP; 
  TH1F *h_preselection_met_SS;     
  TH1F *h_preselection_met_SS_FNP;

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
  
  std::map<Int_t,std::vector<TString>> evtype;
 
  std::map<TString,Double_t> b_fake_wgt;
  int yr;
  // TTree *TF;
  // TFile *ff;
  // TTree *T;
  // TBranch *bpt;
  // TBranch *newBranch; 

  // Float_t b_WeightEvents;

  make2L0JTree *twoLzeroJ;
  
  std::map<int,std::vector< TTreeReaderValue<int> >> triglist_e;
  std::map<int,std::vector< TTreeReaderValue<int> >> triglist_m;
  std::map<int,std::vector<ROOT::Internal::TTreeReaderArrayBase*> > trigmatch_e;
  std::map<int,std::vector<ROOT::Internal::TTreeReaderArrayBase*> > trigmatch_m;
  std::map<TString,TString> trigcat;
  std::map<int, std::vector<TString> > trigstr_e;
  std::map<int, std::vector<TString> > trigstr_m;

  std::vector<TString> trig_order;
  
  
  // Readers to access the data (delete the ones you do not need).
  TTreeReaderValue<ULong64_t> EventNumber = {fReader, "EventNumber"};
  TTreeReaderValue<UInt_t> RunNumber = {fReader, "RunNumber"};
  TTreeReaderValue<Float_t> GenFiltMET = {fReader, "GenFiltMET"};
  TTreeReaderValue<Float_t> GenFiltHT = {fReader, "GenFiltHT"};
  //TTreeReaderValue<UInt_t> LumiBlock = {fReader, "LumiBlockNumber"};
  TTreeReaderValue<Double_t> xsec = {fReader, "xsec"};
  TTreeReaderValue<Float_t> WeightEvents = {fReader, "WeightEvents"};
  TTreeReaderValue<Float_t> WeightEventsPU_Nominal = {fReader, "WeightEventsPU_Nominal"};
  TTreeReaderValue<Float_t> WeightEventselSF_Nominal = {fReader, "WeightEventselSF_Nominal"};
  TTreeReaderValue<Float_t> WeightEventsmuSF_Nominal = {fReader, "WeightEventsmuSF_Nominal"};
  TTreeReaderValue<Float_t> WeightEventsSF_global_Nominal = {fReader, "WeightEventsSF_global_Nominal"};
  TTreeReaderValue<Float_t> WeightEvents_trigger_single_Nominal = {fReader, "WeightEvents_trigger_single_Nominal"};
  TTreeReaderValue<Float_t> WeightEvents_trigger_global_Nominal = {fReader, "WeightEvents_trigger_global_Nominal"};
  TTreeReaderValue<Float_t> WeightEventsbTag_Nominal = {fReader, "WeightEventsbTag_Nominal"};
  TTreeReaderValue<Float_t> WeightEventsJVT_Nominal = {fReader, "WeightEventsJVT_Nominal"};
  TTreeReaderValue<Float_t> AverageInteractionsPerCrossing = {fReader, "AverageInteractionsPerCrossing"};
  TTreeReaderValue<Int_t> NumPrimaryVertices = {fReader, "NumPrimaryVertices"};
  TTreeReaderValue<Int_t> DecayType = {fReader, "DecayType"};
  TTreeReaderValue<Int_t> HFclass = {fReader, "HFclass"};
  TTreeReaderValue<Bool_t> passMETtrig = {fReader, "passMETtrig"};
  TTreeReaderValue<Int_t> TruthPDGID1 = {fReader, "TruthPDGID1"};
  TTreeReaderValue<Int_t> TruthPDGID2 = {fReader, "TruthPDGID2"};
  TTreeReaderValue<Int_t> TruthPDFID1 = {fReader, "TruthPDFID1"};
  TTreeReaderValue<Int_t> TruthPDFID2 = {fReader, "TruthPDFID2"};
  TTreeReaderValue<Float_t> TruthX1 = {fReader, "TruthX1"};
  TTreeReaderValue<Float_t> TruthX2 = {fReader, "TruthX2"};
  TTreeReaderValue<Float_t> TruthXF1 = {fReader, "TruthXF1"};
  TTreeReaderValue<Float_t> TruthXF2 = {fReader, "TruthXF2"};
  TTreeReaderValue<Float_t> TruthQ = {fReader, "TruthQ"};
  TTreeReaderArray<int> labeljets_reco = {fReader, "labeljets_reco"};
  TTreeReaderArray<float> ptjets_Nominal = {fReader, "ptjets_Nominal"};
  TTreeReaderArray<float> etajets_Nominal = {fReader, "etajets_Nominal"};
  TTreeReaderArray<float> phijets_Nominal = {fReader, "phijets_Nominal"};
  TTreeReaderArray<float> massjets_Nominal = {fReader, "massjets_Nominal"};
  TTreeReaderArray<int> isBjets_Nominal = {fReader, "isBjets_Nominal"};
  TTreeReaderArray<float> BtagWeightjets_Nominal = {fReader, "BtagWeightjets_Nominal"};
  TTreeReaderArray<float> JVTjets_Nominal = {fReader, "JVTjets_Nominal"};
  TTreeReaderValue<vector<bool>> passORjet_Nominal = {fReader, "passORjet_Nominal"};
  TTreeReaderValue<vector<bool>> isHiggsjet_Nominal = {fReader, "isHiggsjet_Nominal"};
  TTreeReaderArray<int> isWhichjet_Nominal = {fReader, "isWhichjet_Nominal"};
  TTreeReaderArray<float> ptleptons_Nominal = {fReader, "ptleptons_Nominal"};
  TTreeReaderArray<float> etaleptons_Nominal = {fReader, "etaleptons_Nominal"};
  TTreeReaderArray<float> phileptons_Nominal = {fReader, "phileptons_Nominal"};
  TTreeReaderArray<float> massleptons_Nominal = {fReader, "massleptons_Nominal"};
  TTreeReaderArray<float> ptcone20_Nominal = {fReader, "ptcone20_Nominal"};
  TTreeReaderArray<float> topoetcone20_Nominal = {fReader, "topoetcone20_Nominal"};
  TTreeReaderArray<float> ptvarcone20_Nominal = {fReader, "ptvarcone20_Nominal"};
  TTreeReaderArray<float> ptvarcone30_TightTTVA_pt1000_Nominal = {fReader, "ptvarcone30_TightTTVA_pt1000_Nominal"};
  TTreeReaderArray<float> ptcone20_TightTTVA_pt1000_Nominal = {fReader, "ptcone20_TightTTVA_pt1000_Nominal"};
  TTreeReaderArray<int> flavlep_Nominal = {fReader, "flavlep_Nominal"};
  TTreeReaderValue<vector<bool>> passORlep_Nominal = {fReader, "passORlep_Nominal"};
  TTreeReaderValue<vector<bool>> isSignallep_Nominal = {fReader, "isSignallep_Nominal"};
  TTreeReaderValue<vector<bool>> isHighPtlep_Nominal = {fReader, "isHighPtlep_Nominal"};
  TTreeReaderValue<vector<bool>> isMediumlep_Nominal = {fReader, "isMediumlep_Nominal"};
  TTreeReaderValue<vector<bool>> isTightlep_Nominal = {fReader, "isTightlep_Nominal"};
  TTreeReaderValue<vector<bool>> LepIsoFCHighPtCaloOnly_Nominal = {fReader, "LepIsoFCHighPtCaloOnly_Nominal"};
  TTreeReaderValue<vector<bool>> LepIsoFixedCutHighPtTrackOnly_Nominal = {fReader, "LepIsoFixedCutHighPtTrackOnly_Nominal"};
  TTreeReaderValue<vector<bool>> LepIsoGradient_Nominal = {fReader, "LepIsoGradient_Nominal"};
  TTreeReaderValue<vector<bool>> LepIsoFCTightTrackOnly_Nominal = {fReader, "LepIsoFCTightTrackOnly_Nominal"};
  TTreeReaderValue<vector<bool>> LepIsoFCLoose_Nominal = {fReader, "LepIsoFCLoose_Nominal"};
  TTreeReaderValue<vector<bool>> LepIsoFCTight_Nominal = {fReader, "LepIsoFCTight_Nominal"};
  TTreeReaderArray<float> PromptLeptonIso_Nominal = {fReader, "PromptLeptonIso_Nominal"};
  TTreeReaderArray<float> PromptLeptonVeto_Nominal = {fReader, "PromptLeptonVeto_Nominal"};
  TTreeReaderArray<float> d0sig_Nominal = {fReader, "d0sig_Nominal"};
  TTreeReaderArray<float> z0sinTheta_Nominal = {fReader, "z0sinTheta_Nominal"};
  TTreeReaderValue<Bool_t> cleaningVeto_Nominal = {fReader, "cleaningVeto_Nominal"};
  TTreeReaderValue<Float_t> EtMiss_tstPhi_Nominal = {fReader, "EtMiss_tstPhi_Nominal"};
  TTreeReaderValue<Float_t> EtMiss_tst_Nominal = {fReader, "EtMiss_tst_Nominal"};
  TTreeReaderValue<Float_t> EtMiss_SigObj_Nominal = {fReader, "EtMiss_SigObj_Nominal"};
  //TTreeReaderValue<Float_t> Etmiss_PVSoftTrkPhi = {fReader, "Etmiss_PVSoftTrkPhi"};
  //TTreeReaderValue<Float_t> Etmiss_PVSoftTrk = {fReader, "Etmiss_PVSoftTrk"};
  TTreeReaderValue<Int_t> HLT_2e12_lhloose_L12EM10VH_Nominal = {fReader, "HLT_2e12_lhloose_L12EM10VH_Nominal"};
  TTreeReaderArray<int> HLT_2e12_lhloose_L12EM10VH_match_Nominal = {fReader, "HLT_2e12_lhloose_L12EM10VH_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_2e15_lhvloose_nod0_L12EM13VH_Nominal = {fReader, "HLT_2e15_lhvloose_nod0_L12EM13VH_Nominal"};
  TTreeReaderArray<int> HLT_2e15_lhvloose_nod0_L12EM13VH_match_Nominal = {fReader, "HLT_2e15_lhvloose_nod0_L12EM13VH_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_2e17_lhvloose_nod0_L12EM15VHI_Nominal = {fReader, "HLT_2e17_lhvloose_nod0_L12EM15VHI_Nominal"};
  TTreeReaderArray<int> HLT_2e17_lhvloose_nod0_L12EM15VHI_match_Nominal = {fReader, "HLT_2e17_lhvloose_nod0_L12EM15VHI_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_2e17_lhvloose_nod0_Nominal = {fReader, "HLT_2e17_lhvloose_nod0_Nominal"};
  TTreeReaderArray<int> HLT_2e17_lhvloose_nod0_match_Nominal = {fReader, "HLT_2e17_lhvloose_nod0_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_2e24_lhvloose_nod0_Nominal = {fReader, "HLT_2e24_lhvloose_nod0_Nominal"};
  TTreeReaderArray<int> HLT_2e24_lhvloose_nod0_match_Nominal = {fReader, "HLT_2e24_lhvloose_nod0_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_2g20_tight_icalovloose_L12EM15VHI_Nominal = {fReader, "HLT_2g20_tight_icalovloose_L12EM15VHI_Nominal"};
  TTreeReaderArray<int> HLT_2g20_tight_icalovloose_L12EM15VHI_match_Nominal = {fReader, "HLT_2g20_tight_icalovloose_L12EM15VHI_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_2g22_tight_L12EM15VHI_Nominal = {fReader, "HLT_2g22_tight_L12EM15VHI_Nominal"};
  TTreeReaderArray<int> HLT_2g22_tight_L12EM15VHI_match_Nominal = {fReader, "HLT_2g22_tight_L12EM15VHI_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_2g25_loose_g15_loose_Nominal = {fReader, "HLT_2g25_loose_g15_loose_Nominal"};
  TTreeReaderArray<int> HLT_2g25_loose_g15_loose_match_Nominal = {fReader, "HLT_2g25_loose_g15_loose_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_2g25_tight_L12EM20VH_Nominal = {fReader, "HLT_2g25_tight_L12EM20VH_Nominal"};
  TTreeReaderArray<int> HLT_2g25_tight_L12EM20VH_match_Nominal = {fReader, "HLT_2g25_tight_L12EM20VH_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_2g50_loose_L12EM20VH_Nominal = {fReader, "HLT_2g50_loose_L12EM20VH_Nominal"};
  TTreeReaderArray<int> HLT_2g50_loose_L12EM20VH_match_Nominal = {fReader, "HLT_2g50_loose_L12EM20VH_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_2g20_tight_Nominal = {fReader, "HLT_2g20_tight_Nominal"};
  TTreeReaderArray<int> HLT_2g20_tight_match_Nominal = {fReader, "HLT_2g20_tight_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_2g20_loose_g15_loose_Nominal = {fReader, "HLT_2g20_loose_g15_loose_Nominal"};
  TTreeReaderArray<int> HLT_2g20_loose_g15_loose_match_Nominal = {fReader, "HLT_2g20_loose_g15_loose_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_2g22_tight_Nominal = {fReader, "HLT_2g22_tight_Nominal"};
  TTreeReaderArray<int> HLT_2g22_tight_match_Nominal = {fReader, "HLT_2g22_tight_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_2mu10_Nominal = {fReader, "HLT_2mu10_Nominal"};
  TTreeReaderArray<int> HLT_2mu10_match_Nominal = {fReader, "HLT_2mu10_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_2mu14_Nominal = {fReader, "HLT_2mu14_Nominal"};
  TTreeReaderArray<int> HLT_2mu14_match_Nominal = {fReader, "HLT_2mu14_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_e120_lhloose_Nominal = {fReader, "HLT_e120_lhloose_Nominal"};
  TTreeReaderArray<int> HLT_e120_lhloose_match_Nominal = {fReader, "HLT_e120_lhloose_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_e140_lhloose_nod0_Nominal = {fReader, "HLT_e140_lhloose_nod0_Nominal"};
  TTreeReaderArray<int> HLT_e140_lhloose_nod0_match_Nominal = {fReader, "HLT_e140_lhloose_nod0_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_e17_lhloose_2e9_lhloose_Nominal = {fReader, "HLT_e17_lhloose_2e9_lhloose_Nominal"};
  TTreeReaderArray<int> HLT_e17_lhloose_2e9_lhloose_match_Nominal = {fReader, "HLT_e17_lhloose_2e9_lhloose_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_e17_lhloose_mu14_Nominal = {fReader, "HLT_e17_lhloose_mu14_Nominal"};
  TTreeReaderArray<int> HLT_e17_lhloose_mu14_match_Nominal = {fReader, "HLT_e17_lhloose_mu14_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_e17_lhloose_nod0_mu14_Nominal = {fReader, "HLT_e17_lhloose_nod0_mu14_Nominal"};
  TTreeReaderArray<int> HLT_e17_lhloose_nod0_mu14_match_Nominal = {fReader, "HLT_e17_lhloose_nod0_mu14_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_e24_lhmedium_L1EM20VH_Nominal = {fReader, "HLT_e24_lhmedium_L1EM20VH_Nominal"};
  TTreeReaderArray<int> HLT_e24_lhmedium_L1EM20VH_match_Nominal = {fReader, "HLT_e24_lhmedium_L1EM20VH_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_e24_lhmedium_nod0_L1EM20VH_Nominal = {fReader, "HLT_e24_lhmedium_nod0_L1EM20VH_Nominal"};
  TTreeReaderArray<int> HLT_e24_lhmedium_nod0_L1EM20VH_match_Nominal = {fReader, "HLT_e24_lhmedium_nod0_L1EM20VH_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1_Nominal = {fReader, "HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1_Nominal"};
  TTreeReaderArray<int> HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1_match_Nominal = {fReader, "HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_e24_lhtight_nod0_ivarloose_Nominal = {fReader, "HLT_e24_lhtight_nod0_ivarloose_Nominal"};
  TTreeReaderArray<int> HLT_e24_lhtight_nod0_ivarloose_match_Nominal = {fReader, "HLT_e24_lhtight_nod0_ivarloose_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_e24_lhvloose_nod0_2e12_lhvloose_nod0_L1EM20VH_3EM10VH_Nominal = {fReader, "HLT_e24_lhvloose_nod0_2e12_lhvloose_nod0_L1EM20VH_3EM10VH_Nominal"};
  TTreeReaderArray<int> HLT_e24_lhvloose_nod0_2e12_lhvloose_nod0_L1EM20VH_3EM10VH_match_Nominal = {fReader, "HLT_e24_lhvloose_nod0_2e12_lhvloose_nod0_L1EM20VH_3EM10VH_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_e26_lhmedium_nod0_mu8noL1_Nominal = {fReader, "HLT_e26_lhmedium_nod0_mu8noL1_Nominal"};
  TTreeReaderArray<int> HLT_e26_lhmedium_nod0_mu8noL1_match_Nominal = {fReader, "HLT_e26_lhmedium_nod0_mu8noL1_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_e26_lhtight_nod0_ivarloose_Nominal = {fReader, "HLT_e26_lhtight_nod0_ivarloose_Nominal"};
  TTreeReaderArray<int> HLT_e26_lhtight_nod0_ivarloose_match_Nominal = {fReader, "HLT_e26_lhtight_nod0_ivarloose_match_Nominal"};
  // TTreeReaderValue<Int_t> HLT_e300_etcut_Nominal = {fReader, "HLT_e300_etcut_Nominal"};
  // TTreeReaderArray<int> HLT_e300_etcut_match_Nominal = {fReader, "HLT_e300_etcut_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_e60_lhmedium_nod0_Nominal = {fReader, "HLT_e60_lhmedium_nod0_Nominal"};
  TTreeReaderArray<int> HLT_e60_lhmedium_nod0_match_Nominal = {fReader, "HLT_e60_lhmedium_nod0_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_e60_lhmedium_Nominal = {fReader, "HLT_e60_lhmedium_Nominal"};
  TTreeReaderArray<int> HLT_e60_lhmedium_match_Nominal = {fReader, "HLT_e60_lhmedium_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_e7_lhmedium_mu24_Nominal = {fReader, "HLT_e7_lhmedium_mu24_Nominal"};
  TTreeReaderArray<int> HLT_e7_lhmedium_mu24_match_Nominal = {fReader, "HLT_e7_lhmedium_mu24_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_e7_lhmedium_nod0_mu24_Nominal = {fReader, "HLT_e7_lhmedium_nod0_mu24_Nominal"};
  TTreeReaderArray<int> HLT_e7_lhmedium_nod0_mu24_match_Nominal = {fReader, "HLT_e7_lhmedium_nod0_mu24_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_g140_loose_Nominal = {fReader, "HLT_g140_loose_Nominal"};
  TTreeReaderArray<int> HLT_g140_loose_match_Nominal = {fReader, "HLT_g140_loose_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_g300_etcutHLT_g35_medium_g25_medium_L12EM20VH_Nominal = {fReader, "HLT_g300_etcutHLT_g35_medium_g25_medium_L12EM20VH_Nominal"};
  TTreeReaderArray<int> HLT_g300_etcutHLT_g35_medium_g25_medium_L12EM20VH_match_Nominal = {fReader, "HLT_g300_etcutHLT_g35_medium_g25_medium_L12EM20VH_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_j80_xe80_Nominal = {fReader, "HLT_j80_xe80_Nominal"};
  TTreeReaderArray<int> HLT_j80_xe80_match_Nominal = {fReader, "HLT_j80_xe80_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_mu18_mu8noL1_Nominal = {fReader, "HLT_mu18_mu8noL1_Nominal"};
  TTreeReaderArray<int> HLT_mu18_mu8noL1_match_Nominal = {fReader, "HLT_mu18_mu8noL1_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_mu20_iloose_L1MU15_Nominal = {fReader, "HLT_mu20_iloose_L1MU15_Nominal"};
  TTreeReaderArray<int> HLT_mu20_iloose_L1MU15_match_Nominal = {fReader, "HLT_mu20_iloose_L1MU15_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_mu20_mu8noL1_Nominal = {fReader, "HLT_mu20_mu8noL1_Nominal"};
  TTreeReaderArray<int> HLT_mu20_mu8noL1_match_Nominal = {fReader, "HLT_mu20_mu8noL1_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_mu22_mu8noL1_calotag_0eta010_Nominal = {fReader, "HLT_mu22_mu8noL1_calotag_0eta010_Nominal"};
  TTreeReaderArray<int> HLT_mu22_mu8noL1_calotag_0eta010_match_Nominal = {fReader, "HLT_mu22_mu8noL1_calotag_0eta010_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_mu22_mu8noL1_Nominal = {fReader, "HLT_mu22_mu8noL1_Nominal"};
  TTreeReaderArray<int> HLT_mu22_mu8noL1_match_Nominal = {fReader, "HLT_mu22_mu8noL1_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_mu24_ivarloose_L1MU15_Nominal = {fReader, "HLT_mu24_ivarloose_L1MU15_Nominal"};
  TTreeReaderArray<int> HLT_mu24_ivarloose_L1MU15_match_Nominal = {fReader, "HLT_mu24_ivarloose_L1MU15_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_mu24_ivarloose_Nominal = {fReader, "HLT_mu24_ivarloose_Nominal"};
  TTreeReaderArray<int> HLT_mu24_ivarloose_match_Nominal = {fReader, "HLT_mu24_ivarloose_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_mu24_ivarmedium_Nominal = {fReader, "HLT_mu24_ivarmedium_Nominal"};
  TTreeReaderArray<int> HLT_mu24_ivarmedium_match_Nominal = {fReader, "HLT_mu24_ivarmedium_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_mu26_ivarmedium_Nominal = {fReader, "HLT_mu26_ivarmedium_Nominal"};
  TTreeReaderArray<int> HLT_mu26_ivarmedium_match_Nominal = {fReader, "HLT_mu26_ivarmedium_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_mu40_Nominal = {fReader, "HLT_mu40_Nominal"};
  TTreeReaderArray<int> HLT_mu40_match_Nominal = {fReader, "HLT_mu40_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_mu50_Nominal = {fReader, "HLT_mu50_Nominal"};
  TTreeReaderArray<int> HLT_mu50_match_Nominal = {fReader, "HLT_mu50_match_Nominal"};
  // TTreeReaderValue<Int_t> HLT_mu60_0eta105_msonly_Nominal = {fReader, "HLT_mu60_0eta105_msonly_Nominal"};
  // TTreeReaderArray<int> HLT_mu60_0eta105_msonly_match_Nominal = {fReader, "HLT_mu60_0eta105_msonly_match_Nominal"};
  TTreeReaderValue<Int_t> HLT_xe100_L1XE50_Nominal = {fReader, "HLT_xe100_L1XE50_Nominal"};
  TTreeReaderValue<Int_t> HLT_xe100_mht_L1XE50_Nominal = {fReader, "HLT_xe100_mht_L1XE50_Nominal"};
  TTreeReaderValue<Int_t> HLT_xe100_pufit_L1XE50_Nominal = {fReader, "HLT_xe100_pufit_L1XE50_Nominal"};
  TTreeReaderValue<Int_t> HLT_xe100_pufit_L1XE55_Nominal = {fReader, "HLT_xe100_pufit_L1XE55_Nominal"};
  TTreeReaderValue<Int_t> HLT_xe100_tc_em_L1XE50_Nominal = {fReader, "HLT_xe100_tc_em_L1XE50_Nominal"};
  TTreeReaderValue<Int_t> HLT_xe110_mht_L1XE50_Nominal = {fReader, "HLT_xe110_mht_L1XE50_Nominal"};
  TTreeReaderValue<Int_t> HLT_xe110_pueta_L1XE50_Nominal = {fReader, "HLT_xe110_pueta_L1XE50_Nominal"};
  TTreeReaderValue<Int_t> HLT_xe110_pufit_L1XE50_Nominal = {fReader, "HLT_xe110_pufit_L1XE50_Nominal"};
  TTreeReaderValue<Int_t> HLT_xe110_pufit_L1XE55_Nominal = {fReader, "HLT_xe110_pufit_L1XE55_Nominal"};
  TTreeReaderValue<Int_t> HLT_xe110_pufit_xe65_L1XE50_Nominal = {fReader, "HLT_xe110_pufit_xe65_L1XE50_Nominal"};
  TTreeReaderValue<Int_t> HLT_xe110_pufit_xe70_L1XE50_Nominal = {fReader, "HLT_xe110_pufit_xe70_L1XE50_Nominal"};
  TTreeReaderValue<Int_t> HLT_xe120_pueta_Nominal = {fReader, "HLT_xe120_pueta_Nominal"};
  TTreeReaderValue<Int_t> HLT_xe120_pufit_L1XE50_Nominal = {fReader, "HLT_xe120_pufit_L1XE50_Nominal"};
  TTreeReaderValue<Int_t> HLT_xe120_pufit_Nominal = {fReader, "HLT_xe120_pufit_Nominal"};
  TTreeReaderValue<Int_t> HLT_xe120_tc_lcw_L1XE50_Nominal = {fReader, "HLT_xe120_tc_lcw_L1XE50_Nominal"};
  TTreeReaderValue<Int_t> HLT_xe70_mht_Nominal = {fReader, "HLT_xe70_mht_Nominal"};
  TTreeReaderValue<Int_t> HLT_xe70_tc_lcw_Nominal = {fReader, "HLT_xe70_tc_lcw_Nominal"};
  TTreeReaderValue<Int_t> HLT_xe80_tc_lcw_L1XE50_Nominal = {fReader, "HLT_xe80_tc_lcw_L1XE50_Nominal"};
  TTreeReaderValue<Int_t> HLT_xe90_mht_L1XE50_Nominal = {fReader, "HLT_xe90_mht_L1XE50_Nominal"};
  TTreeReaderValue<Int_t> HLT_xe90_mht_wEFMu_L1XE50_Nominal = {fReader, "HLT_xe90_mht_wEFMu_L1XE50_Nominal"};
  TTreeReaderValue<Int_t> HLT_xe90_pufit_L1XE50_Nominal = {fReader, "HLT_xe90_pufit_L1XE50_Nominal"};
  TTreeReaderValue<Int_t> HLT_xe90_tc_lcw_wEFMu_L1XE50_Nominal = {fReader, "HLT_xe90_tc_lcw_wEFMu_L1XE50_Nominal"};


  myTree(TTree * /*tree*/ =0) { }
  virtual ~myTree() { }
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
  bool lepL(int idx);
  bool lepT(int idx, int vb = 0);
  std::pair <TString,TString> getFillString(int idx1, int idx2);
  void classifyEvent(int idx1, int idx2, int& isTT, int& isTL, int& isLT, int& isLL, int vb = 0);
  bool checkEventTrigger(int id, int yr);
  std::vector<TString> getTriggerCat(int id, int yr);
  std::vector<TString> checkTriggerMatch(int id, int yr);
  void fill2L0JTree(int idx1, int idx2, std::map<TString,double> fake_weight, Int_t n_bjets, Int_t n_ljets, Int_t n_siglep, Float_t mll, bool isEE, bool isEM, bool isMM, bool isOS, bool isSS);
  bool bjet(int idx);
  bool jetT(int idx);
  bool jetL(int idx);
  Int_t getPtThresholdTrigger(TString tname, int vb);

  
  ClassDef(myTree,0);

};

#endif

#ifdef myTree_cxx
void myTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t myTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef myTree_cxx
