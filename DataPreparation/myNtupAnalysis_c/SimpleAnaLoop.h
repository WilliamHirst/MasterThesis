//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri May 20 14:34:33 2022 by ROOT version 6.24/02
// from TChain Zmmjets_NoSys/
//////////////////////////////////////////////////////////

#ifndef SimpleAnaLoop_h
#define SimpleAnaLoop_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include <vector>
#include <map>

#include "TString.h"



class SimpleAnaLoop : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain
  
  std::map<Int_t,Float_t> dsidvec;
  Int_t nev;
  Int_t nentries;
  // Readers to access the data (delete the ones you do not need).
  /**
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
  */
   TTreeReaderValue<Float_t> xsec = {fReader, "xsec"};
  /**
TTreeReaderValue<Float_t> GenHt = {fReader, "GenHt"};
   TTreeReaderValue<Float_t> GenMET = {fReader, "GenMET"};
  */
   TTreeReaderValue<Int_t> DatasetNumber = {fReader, "DatasetNumber"};
  /**
  TTreeReaderValue<Int_t> RunNumber = {fReader, "RunNumber"};
   TTreeReaderValue<Int_t> RandomRunNumber = {fReader, "RandomRunNumber"};
   TTreeReaderValue<Int_t> FS = {fReader, "FS"};
   TTreeReaderArray<float> LHE3Weights = {fReader, "LHE3Weights"};
   TTreeReaderArray<TString> LHE3WeightNames = {fReader, "LHE3WeightNames"};
  */

   SimpleAnaLoop(TTree * /*tree*/ =0) { }
   virtual ~SimpleAnaLoop() { }
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

   ClassDef(SimpleAnaLoop,0);

};

#endif

#ifdef SimpleAnaLoop_cxx
void SimpleAnaLoop::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t SimpleAnaLoop::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef SimpleAnaLoop_cxx
