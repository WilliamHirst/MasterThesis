//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Sep 22 10:04:17 2021 by ROOT version 6.18/00
// from TChain mini/
//////////////////////////////////////////////////////////

#ifndef OpenDataLooper_h
#define OpenDataLooper_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include <vector>



class OpenDataLooper : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

  Int_t nentry;
  TString newfile;
  TFile *ntupf;
  TTree *newtree;
  TChain *chain;
  Long64_t nentries;
  TString fullpath;
  TString previous_fullpath;
  char fnameoutputcopyfilename[200];
  int nfilled;
  bool isreported;

  TBranch *bwgt; 
  std::map<Int_t,Double_t> lpx_weights;

  Float_t newxsec;
  Float_t newsumofw;
  double lpx_wgt;
   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Int_t> runNumber = {fReader, "runNumber"};
   TTreeReaderValue<Int_t> eventNumber = {fReader, "eventNumber"};
   TTreeReaderValue<Int_t> channelNumber = {fReader, "channelNumber"};
   TTreeReaderValue<Float_t> mcWeight = {fReader, "mcWeight"};
   TTreeReaderValue<Float_t> scaleFactor_PILEUP = {fReader, "scaleFactor_PILEUP"};
   TTreeReaderValue<Float_t> scaleFactor_ELE = {fReader, "scaleFactor_ELE"};
   TTreeReaderValue<Float_t> scaleFactor_MUON = {fReader, "scaleFactor_MUON"};
   TTreeReaderValue<Float_t> scaleFactor_PHOTON = {fReader, "scaleFactor_PHOTON"};
   TTreeReaderValue<Float_t> scaleFactor_TAU = {fReader, "scaleFactor_TAU"};
   TTreeReaderValue<Float_t> scaleFactor_BTAG = {fReader, "scaleFactor_BTAG"};
   TTreeReaderValue<Float_t> scaleFactor_LepTRIGGER = {fReader, "scaleFactor_LepTRIGGER"};
   TTreeReaderValue<Float_t> scaleFactor_PhotonTRIGGER = {fReader, "scaleFactor_PhotonTRIGGER"};
   TTreeReaderValue<Bool_t> trigE = {fReader, "trigE"};
   TTreeReaderValue<Bool_t> trigM = {fReader, "trigM"};
   TTreeReaderValue<Bool_t> trigP = {fReader, "trigP"};
   TTreeReaderValue<UInt_t> lep_n = {fReader, "lep_n"};
   TTreeReaderValue<vector<bool>> lep_truthMatched = {fReader, "lep_truthMatched"};
   TTreeReaderValue<vector<bool>> lep_trigMatched = {fReader, "lep_trigMatched"};
   TTreeReaderArray<float> lep_pt = {fReader, "lep_pt"};
   TTreeReaderArray<float> lep_eta = {fReader, "lep_eta"};
   TTreeReaderArray<float> lep_phi = {fReader, "lep_phi"};
   TTreeReaderArray<float> lep_E = {fReader, "lep_E"};
   TTreeReaderArray<float> lep_z0 = {fReader, "lep_z0"};
   TTreeReaderArray<int> lep_charge = {fReader, "lep_charge"};
   TTreeReaderArray<unsigned int> lep_type = {fReader, "lep_type"};
   TTreeReaderValue<vector<bool>> lep_isTightID = {fReader, "lep_isTightID"};
   TTreeReaderArray<float> lep_ptcone30 = {fReader, "lep_ptcone30"};
   TTreeReaderArray<float> lep_etcone20 = {fReader, "lep_etcone20"};
   TTreeReaderArray<float> lep_trackd0pvunbiased = {fReader, "lep_trackd0pvunbiased"};
   TTreeReaderArray<float> lep_tracksigd0pvunbiased = {fReader, "lep_tracksigd0pvunbiased"};
   TTreeReaderValue<Float_t> met_et = {fReader, "met_et"};
   TTreeReaderValue<Float_t> met_phi = {fReader, "met_phi"};
   TTreeReaderValue<UInt_t> jet_n = {fReader, "jet_n"};
   TTreeReaderArray<float> jet_pt = {fReader, "jet_pt"};
   TTreeReaderArray<float> jet_eta = {fReader, "jet_eta"};
   TTreeReaderArray<float> jet_phi = {fReader, "jet_phi"};
   TTreeReaderArray<float> jet_E = {fReader, "jet_E"};
   TTreeReaderArray<float> jet_jvt = {fReader, "jet_jvt"};
   TTreeReaderArray<int> jet_trueflav = {fReader, "jet_trueflav"};
   TTreeReaderValue<vector<bool>> jet_truthMatched = {fReader, "jet_truthMatched"};
   TTreeReaderArray<float> jet_MV2c10 = {fReader, "jet_MV2c10"};
   TTreeReaderValue<UInt_t> photon_n = {fReader, "photon_n"};
   TTreeReaderValue<vector<bool>> photon_truthMatched = {fReader, "photon_truthMatched"};
   TTreeReaderValue<vector<bool>> photon_trigMatched = {fReader, "photon_trigMatched"};
   TTreeReaderArray<float> photon_pt = {fReader, "photon_pt"};
   TTreeReaderArray<float> photon_eta = {fReader, "photon_eta"};
   TTreeReaderArray<float> photon_phi = {fReader, "photon_phi"};
   TTreeReaderArray<float> photon_E = {fReader, "photon_E"};
   TTreeReaderValue<vector<bool>> photon_isTightID = {fReader, "photon_isTightID"};
   TTreeReaderArray<float> photon_ptcone30 = {fReader, "photon_ptcone30"};
   TTreeReaderArray<float> photon_etcone20 = {fReader, "photon_etcone20"};
   TTreeReaderArray<int> photon_convType = {fReader, "photon_convType"};
   TTreeReaderValue<UInt_t> tau_n = {fReader, "tau_n"};
   TTreeReaderArray<float> tau_pt = {fReader, "tau_pt"};
   TTreeReaderArray<float> tau_eta = {fReader, "tau_eta"};
   TTreeReaderArray<float> tau_phi = {fReader, "tau_phi"};
   TTreeReaderArray<float> tau_E = {fReader, "tau_E"};
   TTreeReaderValue<vector<bool>> tau_isTightID = {fReader, "tau_isTightID"};
   TTreeReaderValue<vector<bool>> tau_truthMatched = {fReader, "tau_truthMatched"};
   TTreeReaderValue<vector<bool>> tau_trigMatched = {fReader, "tau_trigMatched"};
   TTreeReaderArray<int> tau_nTracks = {fReader, "tau_nTracks"};
   TTreeReaderArray<float> tau_BDTid = {fReader, "tau_BDTid"};
   TTreeReaderValue<Float_t> ditau_m = {fReader, "ditau_m"};
   TTreeReaderArray<float> lep_pt_syst = {fReader, "lep_pt_syst"};
   TTreeReaderValue<Float_t> met_et_syst = {fReader, "met_et_syst"};
   TTreeReaderArray<float> jet_pt_syst = {fReader, "jet_pt_syst"};
   TTreeReaderArray<float> photon_pt_syst = {fReader, "photon_pt_syst"};
   TTreeReaderArray<float> tau_pt_syst = {fReader, "tau_pt_syst"};
   TTreeReaderValue<Float_t> XSection = {fReader, "XSection"};
   TTreeReaderValue<Float_t> SumWeights = {fReader, "SumWeights"};
   TTreeReaderValue<UInt_t> largeRjet_n = {fReader, "largeRjet_n"};
   TTreeReaderArray<float> largeRjet_pt = {fReader, "largeRjet_pt"};
   TTreeReaderArray<float> largeRjet_eta = {fReader, "largeRjet_eta"};
   TTreeReaderArray<float> largeRjet_phi = {fReader, "largeRjet_phi"};
   TTreeReaderArray<float> largeRjet_E = {fReader, "largeRjet_E"};
   TTreeReaderArray<float> largeRjet_m = {fReader, "largeRjet_m"};
   TTreeReaderArray<float> largeRjet_truthMatched = {fReader, "largeRjet_truthMatched"};
   TTreeReaderArray<float> largeRjet_D2 = {fReader, "largeRjet_D2"};
   TTreeReaderArray<float> largeRjet_tau32 = {fReader, "largeRjet_tau32"};
   TTreeReaderArray<float> largeRjet_pt_syst = {fReader, "largeRjet_pt_syst"};
   TTreeReaderArray<int> tau_charge = {fReader, "tau_charge"};


   OpenDataLooper(TTree * /*tree*/ =0) { }
   virtual ~OpenDataLooper() { }
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
  void WriteToFile();

   ClassDef(OpenDataLooper,0);

};

#endif

#ifdef OpenDataLooper_cxx
void OpenDataLooper::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t OpenDataLooper::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef OpenDataLooper_cxx
