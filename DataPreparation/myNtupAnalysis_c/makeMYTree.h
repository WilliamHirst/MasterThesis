//------------------------------------------
//   makeMYTrees.h
//   definitions for core analysis
//
//
//   author: Lukas Marti
//-------------------------------------------
#ifndef makeMYTree_h
#define makeMYTree_h

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
//#include "makeMYTree/MultiLepEvent.h"


class makeMYTree  {

 public:

  // Constructor
  // - MCID => will be used to name the tree
  // - syst => will be used to name the file
  // That way all files are hadd'able
  // Optional arguments fileName and treeName can override default naming convention
  makeMYTree(TString MCID, TString syst, TString fileName="", TString treeName="");
  ~makeMYTree();

  TFile* file; 
  TTree* tree;
  TString mymcid;

  std::vector<Double_t> bMY_MM_weight;
  std::vector<Double_t> bMY_syst_MM_up;
  std::vector<Double_t> bMY_syst_MM_down;
  std::vector<TString>  bMY_MM_key;

  Bool_t                                     bMY_trigMatch_1L2LTrig;                              
  Bool_t                                     bMY_trigMatch_1LTrig;                                     
  Bool_t                                     bMY_trigMatch_1L2LTrigOR;                                  
  Bool_t                                     bMY_trigMatch_1LTrigOR;                                    
  Bool_t                                     bMY_trigMatch_2LTrig;                                      
  Bool_t                                     bMY_trigMatch_2LTrigOR;                                    
  Bool_t                                     bMY_trigMatch_HLT_e24_lhmedium_L1EM20VH;                   
  Bool_t                                     bMY_trigMatch_HLT_e60_lhmedium;                            
  Bool_t                                     bMY_trigMatch_HLT_e120_lhloose;                            
  Bool_t                                     bMY_trigMatch_HLT_mu20_iloose_L1MU15;                      
  Bool_t                                     bMY_trigMatch_HLT_2e12_lhloose_L12EM10VH;                  
  Bool_t                                     bMY_trigMatch_HLT_mu18_mu8noL1;                            
  Bool_t                                     bMY_trigMatch_HLT_e17_lhloose_mu14;                        
  Bool_t                                     bMY_trigMatch_HLT_e7_lhmedium_mu24;                        
  Bool_t                                     bMY_trigMatch_HLT_e24_lhtight_nod0_ivarloose;              
  Bool_t                                     bMY_trigMatch_HLT_e24_lhmedium_nod0_L1EM20VH;              
  Bool_t                                     bMY_trigMatch_HLT_e60_medium;                              
  Bool_t                                     bMY_trigMatch_HLT_mu40;                                   
  Bool_t                                     bMY_trigMatch_HLT_mu24_iloose_L1MU15;                      
  Bool_t                                     bMY_trigMatch_HLT_mu24_ivarloose_L1MU15;                  
  Bool_t                                     bMY_trigMatch_HLT_mu24_ivarmedium;                         
  Bool_t                                     bMY_trigMatch_HLT_mu24_imedium;                            
  Bool_t                                     bMY_trigMatch_HLT_mu26_imedium;                            
  Bool_t                                     bMY_trigMatch_HLT_2e15_lhvloose_nod0_L12EM13VH;            
  Bool_t                                     bMY_trigMatch_HLT_2e17_lhvloose_nod0;                      
  Bool_t                                     bMY_trigMatch_HLT_2mu10;                                   
  Bool_t                                     bMY_trigMatch_HLT_2mu14;                                   
  Bool_t                                     bMY_trigMatch_HLT_mu20_mu8noL1;                            
  Bool_t                                     bMY_trigMatch_HLT_mu22_mu8noL1;                            
  Bool_t                                     bMY_trigMatch_HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1;     
  Bool_t                                     bMY_trigMatch_HLT_e26_lhtight_nod0_ivarloose;              
  Bool_t                                     bMY_trigMatch_HLT_e26_lhtight_nod0;                        
  Bool_t                                     bMY_trigMatch_HLT_e60_lhmedium_nod0;                       
  Bool_t                                     bMY_trigMatch_HLT_e140_lhloose_nod0;                       
  Bool_t                                     bMY_trigMatch_HLT_e300_etcut;                              
  Bool_t                                     bMY_trigMatch_HLT_mu26_ivarmedium;                        
  Bool_t                                     bMY_trigMatch_HLT_mu50;                                    
  Bool_t                                     bMY_trigMatch_HLT_mu60_0eta105_msonly;                     
  Bool_t                                     bMY_trigMatch_HLT_2e17_lhvloose_nod0_L12EM15VHI;          
  Bool_t                                     bMY_trigMatch_HLT_2e24_lhvloose_nod0;                      
  Bool_t                                     bMY_trigMatch_HLT_e17_lhloose_nod0_mu14;                   
  Bool_t                                     bMY_trigMatch_HLT_e26_lhmedium_nod0_mu8noL1;               
  Bool_t                                     bMY_trigMatch_HLT_e7_lhmedium_nod0_mu24;                   
  Float_t                                    bMY_mu;                                                    
  Float_t                                    bMY_avg_mu;                                                
  Float_t                                    bMY_actual_mu;                                             
  Int_t                                      bMY_nVtx;                                                  
  Int_t                                      bMY_channel;                                               
  Int_t                                      bMY_nLep_base;                                             
  Int_t                                      bMY_nLep_signal;                                           
  vector<int>                                bMY_lepFlavor;                                             
  vector<int>                                bMY_lepCharge;                                             
  vector<int>                                bMY_lepAuthor;                                             
  vector<float>                              bMY_lepPt;                                                 
  vector<float>                              bMY_lepEta;                                                
  vector<float>                              bMY_lepPhi;                                                
  vector<float>                              bMY_lepM;                                                  
  vector<float>                              bMY_lepD0;                                                 
  vector<float>                              bMY_lepD0Sig;                                              
  vector<float>                              bMY_lepZ0;                                                 
  vector<float>                              bMY_lepZ0SinTheta;                                         
  vector<float>                              bMY_lepTopoetcone20;                                       
  vector<float>                              bMY_lepTopoetcone30;                                       
  vector<float>                              bMY_lepPtvarcone20;                                        
  vector<float>                              bMY_lepPtvarcone30;                                        
  vector<float>                              bMY_lepCorrTopoetcone20;                                   
  vector<float>                              bMY_lepCorrPtvarcone20;                                    
  vector<float>                              bMY_lepCorrPtvarcone30;                                    
  vector<bool>                               bMY_lepPassOR;                                          
  vector<int>                                bMY_lepType;                                               
  vector<int>                                bMY_lepOrigin;                                             
  vector<int>                                bMY_lepEgMotherType;                                       
  vector<int>                                bMY_lepEgMotherOrigin;                                     
  vector<int>                                bMY_lepEgMotherPdgId;                                      
  vector<int>                                bMY_lepECIDS;                                              
  vector<bool>                               bMY_lepIsHF;                                               
  vector<bool>                               bMY_lepIsLF;                                               
  vector<bool>                               bMY_lepIsCO;                                               
  vector<bool>                               bMY_lepIsCF;                                               
  vector<bool>                               bMY_lepIsUK;                                               
  vector<bool>                               bMY_lepIsPR;                                               
  vector<bool>                               bMY_lepPassBL;                                             
  vector<bool>                               bMY_lepVeryLoose;                                          
  vector<bool>                               bMY_lepLoose;                                              
  vector<bool>                               bMY_lepMedium;                                             
  vector<bool>                               bMY_lepTight;                                              
  vector<bool>                               bMY_lepIsoFCHighPtCaloOnly;	              
  vector<bool>                               bMY_lepIsoFCLoose;		              
  vector<bool>                               bMY_lepIsoFCTight;		              
  vector<bool>                               bMY_lepIsoFCLoose_FixedRad;	              
  vector<bool>                               bMY_lepIsoFCTight_FixedRad;	              
  vector<bool>                               bMY_lepIsoHighPtCaloOnly;	              
  vector<bool>                               bMY_lepIsoTightTrackOnly_VarRad;	              
  vector<bool>                               bMY_lepIsoTightTrackOnly_FixedRad;             
  vector<bool>                               bMY_lepIsoLoose_VarRad;		              
  vector<bool>                               bMY_lepIsoTight_VarRad;		              
  vector<bool>                               bMY_lepIsoPLVLoose;		              
  vector<bool>                               bMY_lepIsoPLVTight;                                                
  vector<bool>                               bMY_lepTruthMatched;                                       
  vector<int>                                bMY_lepTruthCharge;                                        
  vector<float>                              bMY_lepTruthPt;                                            
  vector<float>                              bMY_lepTruthEta;                                           
  vector<float>                              bMY_lepTruthPhi;                                          
  vector<float>                              bMY_lepTruthM;                                             
  Int_t                                      bMY_nJet30;                                                
  Int_t                                      bMY_nJet20;                                               
  Int_t                                      bMY_nBJet20_MV2c10_FixedCutBEff_77;                        
  vector<float>                              bMY_jetPt;                                                 
  vector<float>                              bMY_jetEta;                                                
  vector<float>                              bMY_jetPhi;                                                
  vector<float>                              bMY_jetM;                                                
  vector<float>                              bMY_jetJVT;                                                
  vector<bool>                               bMY_jetPassOR;                                             
  vector<bool>                               bMY_jetSignal;                                             
  Float_t                                    bMY_mjj;                                                   
  vector<float>                              bMY_jetTileEnergy;                                         
  vector<float>                              bMY_jetMV2c10;                                             
  Float_t                                    bMY_met_Et;                                               
  Float_t                                    bMY_met_Sign;                                              
  Float_t                                    bMY_met_Phi;                                              
  Float_t                                    bMY_met_Et_loose;                                          
  Float_t                                    bMY_met_Et_tighter;                                       
  Float_t                                    bMY_met_Et_tenacious;                                     
  Float_t                                    bMY_met_Phi_loose;                                        
  Float_t                                    bMY_met_Phi_tighter;                                       
  Float_t                                    bMY_met_Phi_tenacious;                                     
  Float_t                                    bMY_mll;                                                   
  Double_t                                   bMY_pileupWeight;                                          
  Double_t                                   bMY_leptonWeight;                                          
  Double_t                                   bMY_eventWeight;                                           
  Double_t                                   bMY_bTagWeight;                                            
  Double_t                                   bMY_jvtWeight;                                             
  Double_t                                   bMY_globalDiLepTrigSF;                                     
  Double_t                                   bMY_flavSymWeight;                                         
  Double_t                                   bMY_genWeight;                                             
  Double_t                                   bMY_genWeightUp;                                           
  Double_t                                   bMY_genWeightDown;                                         
  ULong64_t                                  bMY_PRWHash;                                               
  ULong64_t                                  bMY_EventNumber;                                           
  Float_t                                    bMY_xsec;                                                  
  Float_t                                    bMY_GenHt;                                                 
  Float_t                                    bMY_GenMET;                                                
  Int_t                                      bMY_DatasetNumber;                                         
  Int_t                                      bMY_RunNumber;                                             
  Int_t                                      bMY_RandomRunNumber;                                       
  Int_t                                      bMY_FS;                                                

                                               

  //virtual void InitializeOutput(TFile** file, TString filename,TTree** tree, TString treename );
  void ClearOutputBranches();

  void setSumOfMcWeights(double sumOfMcWeights);


  void WriteTree();

};
#endif // #ifdef makeMYTrees_cxx
