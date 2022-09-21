//------------------------------------------
//   HistFitterTrees.h
//   definitions for core analysis
//
//
//   author: Lukas Marti
//-------------------------------------------
#ifndef HistFitterTree_h
#define HistFitterTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TH1F.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "TLorentzVector.h"
#include "TParameter.h"
//#include "HistFitterTree/MultiLepEvent.h"


class HistFitterTree  {

public:

  // Constructor
  // - MCID => will be used to name the tree
  // - syst => will be used to name the file
  // That way all files are hadd'able
  // Optional arguments fileName and treeName can override default naming convention
  HistFitterTree(TString MCID, TString syst, TString fileName="", TString treeName="");
  ~HistFitterTree();

  TFile* file; 
  TTree* tree;
  TString mymcid;

  UInt_t DatasetNumber;
  UInt_t RunNumber;
  ULong64_t EventNumber;
  Double_t lept1Pt;
  Double_t lept2Pt;
  Double_t lept3Pt;
  Double_t lept1Eta;
  Double_t lept2Eta;
  Double_t lept3Eta;
  Double_t lept1Phi;
  Double_t lept2Phi;
  Double_t lept3Phi;
  Double_t lept1E;
  Double_t lept2E;
  Double_t lept3E;
  Double_t lept1Flav;
  Double_t lept2Flav;
  Double_t lept3Flav;
  Double_t metPhi;
  Double_t MET;
  Double_t METsig;
  Int_t nLeps;
  Int_t nSigLeps;
  Float_t WeightEvents;
  Float_t WeightEventsEff;
  Float_t WeightEventsSF;
  Float_t WeightEventsbTag;
  Float_t WeightEventsJVT;
  Float_t WeightEventselSF;
  Float_t WeightEventsmuSF;
  Double_t L2Mll;
  Double_t L2MT2;
  Bool_t L2isEMU;
  Bool_t L2isEE;
  Bool_t L2isMUMU;
  Int_t L2nCentralLightJet;
  Int_t L2nCentralBJet;
  Int_t L2nForwardJet;
  Int_t L2nCentralLightJet30;
  Int_t L2nCentralBJet30;
  Int_t L2nForwardJet30;
  Int_t L2nCentralLightJet40;
  Int_t L2nCentralBJet40;
  Int_t L2nForwardJet40;
  Int_t L2nCentralLightJet50;
  Int_t L2nCentralBJet50;
  Int_t L2nForwardJet50;
  Int_t L2nCentralLightJet60;
  Int_t L2nCentralBJet60;
  Int_t L2nForwardJet60;
  Bool_t L2SiLepTrigger;
  Bool_t L2DiLepTrigger;
  Double_t jet1Pt;
  Double_t jet2Pt;
  Double_t jet3Pt;
  Double_t jet1Eta;
  Double_t jet2Eta;
  Double_t jet3Eta;
  Double_t jet1Phi;
  Double_t jet2Phi;
  Double_t jet3Phi;
  Float_t b_mu;

  std::vector<Double_t> eventWeight;
  std::vector<Double_t> syst_MM_up;
  std::vector<Double_t> syst_MM_down;
  std::vector<TString> MM_key;

  /**
  Double_t eventWeight;
  Double_t syst_MM_up;
  Double_t syst_MM_down;
  */

  Double_t syst_FT_EFF_B_down;
  Double_t syst_FT_EFF_B_up;
  Double_t syst_FT_EFF_C_down;
  Double_t syst_FT_EFF_C_up;
  Double_t syst_FT_EFF_Light_down;
  Double_t syst_FT_EFF_Light_up;
  Double_t syst_FT_EFF_extrapolation_down;
  Double_t syst_FT_EFF_extrapolation_up;
  Double_t syst_FT_EFF_extrapolation_from_charm_down;
  Double_t syst_FT_EFF_extrapolation_from_charm_up;
  Double_t syst_EL_EFF_ID_TOTAL_UncorrUncertainty_down;
  Double_t syst_EL_EFF_ID_TOTAL_UncorrUncertainty_up;
  Double_t syst_EL_EFF_Iso_TOTAL_UncorrUncertainty_down;
  Double_t syst_EL_EFF_Iso_TOTAL_UncorrUncertainty_up;
  Double_t syst_EL_EFF_Reco_TOTAL_UncorrUncertainty_down;
  Double_t syst_EL_EFF_Reco_TOTAL_UncorrUncertainty_up;
  Double_t syst_EL_EFF_Trigger_TOTAL_UncorrUncertainty_down;
  Double_t syst_EL_EFF_Trigger_TOTAL_UncorrUncertainty_up;
  Double_t syst_MUON_EFF_BADMUON_STAT_down;
  Double_t syst_MUON_EFF_BADMUON_STAT_up;
  Double_t syst_MUON_EFF_BADMUON_SYS_down;
  Double_t syst_MUON_EFF_BADMUON_SYS_up;
  Double_t syst_MUON_EFF_ISO_STAT_down;
  Double_t syst_MUON_EFF_ISO_STAT_up;
  Double_t syst_MUON_EFF_ISO_SYS_down;
  Double_t syst_MUON_EFF_ISO_SYS_up;
  Double_t syst_MUON_EFF_RECO_STAT_down;
  Double_t syst_MUON_EFF_RECO_STAT_up;
  Double_t syst_MUON_EFF_RECO_STAT_LOWPT_down;
  Double_t syst_MUON_EFF_RECO_STAT_LOWPT_up;
  Double_t syst_MUON_EFF_RECO_SYS_down;
  Double_t syst_MUON_EFF_RECO_SYS_up;
  Double_t syst_MUON_EFF_RECO_SYS_LOWPT_down;
  Double_t syst_MUON_EFF_RECO_SYS_LOWPT_up;
  Double_t syst_MUON_EFF_TTVA_STAT_down;
  Double_t syst_MUON_EFF_TTVA_STAT_up;
  Double_t syst_MUON_EFF_TTVA_SYS_down;
  Double_t syst_MUON_EFF_TTVA_SYS_up;
  Double_t syst_MUON_EFF_TrigStatUncertainty_down;
  Double_t syst_MUON_EFF_TrigStatUncertainty_up;
  Double_t syst_MUON_EFF_TrigSystUncertainty_down;
  Double_t syst_MUON_EFF_TrigSystUncertainty_up;
  Double_t syst_JET_JvtEfficiency_1down;
  Double_t syst_JET_JvtEfficiency_1up;
  Double_t syst_JET_fJvtEfficiency_1down;
  Double_t syst_JET_fJvtEfficiency_1up;

  //virtual void InitializeOutput(TFile** file, TString filename,TTree** tree, TString treename );
  void ClearOutputBranches();

  void setSumOfMcWeights(double sumOfMcWeights);


  void WriteTree();

};
#endif // #ifdef HistFitterTrees_cxx
