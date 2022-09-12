//------------------------------------------
//   make2L2JTrees.h
//   definitions for core analysis
//
//
//   author: Lukas Marti
//-------------------------------------------
#ifndef make2L2JTree_h
#define make2L2JTree_h

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
//#include "make2L2JTree/MultiLepEvent.h"


class make2L2JTree  {

public:

  // Constructor
  // - MCID => will be used to name the tree
  // - syst => will be used to name the file
  // That way all files are hadd'able
  // Optional arguments fileName and treeName can override default naming convention
  make2L2JTree(TString MCID, TString syst, TString fileName="", TString treeName="", std::vector<TString> unc_keys = {});
  ~make2L2JTree();

  TFile* file; 
  TTree* tree;
  TString mymcid;

  std::map<TString,Double_t> b2L2J_FNPweight;
  Double_t b2L2J_FNP_TOTAL_UP;
  Double_t b2L2J_FNP_TOTAL_DOWN;
  //Double_t b2L2J_fake_wgt_dw;
  Double_t b2L2J_trigWeight_2LTrig;                                          
  Bool_t b2L2J_trigMatch_2LTrig;                                             
  Float_t b2L2J_mu;                                                          
  Float_t b2L2J_avg_mu;                                                      
  Float_t b2L2J_actual_mu;                                                   
  Int_t b2L2J_nVtx;                                                          
  Int_t b2L2J_channel;
  Int_t b2L2J_nLep_combi;                                                     
  Int_t b2L2J_nLep_base;                                                     
  Int_t b2L2J_nLep_signal;                                                   
  vector<int> b2L2J_lepFlavor;                                               
  vector<int> b2L2J_lepCharge;                                               
  vector<int> b2L2J_lepAuthor;                                               
  vector<float> b2L2J_lepPt;                                                 
  vector<float> b2L2J_lepEta;                                                
  vector<float> b2L2J_lepPhi;                                                
  vector<float> b2L2J_lepM;                                                  
  vector<float> b2L2J_lepD0;                                                 
  vector<float> b2L2J_lepD0Sig;                                              
  vector<float> b2L2J_lepZ0;                                                 
  vector<float> b2L2J_lepZ0SinTheta;                                         
  vector<bool> b2L2J_lepPassOR;                                              
  vector<int> b2L2J_lepType;                                                 
  vector<int> b2L2J_lepOrigin;
  vector<int> b2L2J_lepIFFClass;           
  vector<int> b2L2J_lepEgMotherType;                                         
  vector<int> b2L2J_lepEgMotherOrigin;                                       
  vector<int> b2L2J_lepEgMotherPdgId;                                        
  vector<int> b2L2J_lepECIDS;                                                
  vector<int> b2L2J_lepNPix;                                                                         
  vector<bool> b2L2J_lepPassBL;                                                           
  vector<bool> b2L2J_lepSignal;                                              
  vector<bool> b2L2J_lepTruthMatched;                                        
  vector<int> b2L2J_lepTruthCharge;                                          
  vector<float> b2L2J_lepTruthPt;                                            
  vector<float> b2L2J_lepTruthEta;                                           
  vector<float> b2L2J_lepTruthPhi;                                           
  vector<float> b2L2J_lepTruthM;                                             
  Int_t b2L2J_nJet30;                                                        
  Int_t b2L2J_nJet20;                                                             
  Int_t b2L2J_nBJet20_MV2c10_FixedCutBEff_77;                                          
  vector<float> b2L2J_jetPt;                                                 
  vector<float> b2L2J_jetEta;                                                
  vector<float> b2L2J_jetPhi;                                                
  vector<float> b2L2J_jetM;                                                  
  Float_t b2L2J_mjj;                                                         
  Float_t b2L2J_mjj_minDPhiZMET;                                             
  Float_t b2L2J_Rjj;                                                         
  Float_t b2L2J_Rjj_minDPhiZMET;                                             
  vector<float> b2L2J_jetTileEnergy;                                         
  vector<float> b2L2J_jetMV2c10;                                             
  Float_t b2L2J_vectorSumJetsPt;                                             
  Float_t b2L2J_vectorSumJetsEta;                                            
  Float_t b2L2J_vectorSumJetsPhi;                                            
  Float_t b2L2J_vectorSumJetsM;                                              
  Float_t b2L2J_met_Et;                                                      
  Float_t b2L2J_met_Sign;                                                    
  Float_t b2L2J_met_Phi;                                                     
  Float_t b2L2J_TST_Et;                                                      
  Float_t b2L2J_TST_Phi;                                                                                 
  Float_t b2L2J_deltaPhi_MET_TST_Phi;                                                                  
  Float_t b2L2J_Ht30;
  Float_t b2L2J_mbb;
  Float_t b2L2J_Rbb;
  Float_t b2L2J_METOverPtZ;                                                  
  Float_t b2L2J_METOverPtW;                                                  
  Float_t b2L2J_PtISR;                                                       
  Float_t b2L2J_METOverPtISR;                                                
  Float_t b2L2J_METOverHT;                                                   
  Float_t b2L2J_METOverJ1pT;                                                
  Float_t b2L2J_DPhiJ1Met;                                                   
  Float_t b2L2J_DPhiJ2Met;                                                   
  Float_t b2L2J_DPhiJ3Met;                                                   
  Float_t b2L2J_DPhiJ4Met;                                                   
  Float_t b2L2J_minDPhi2JetsMet;                                             
  Float_t b2L2J_minDPhi4JetsMet;                                             
  Float_t b2L2J_minDPhiAllJetsMet;                                           
  Float_t b2L2J_dPhiPjjMet;                                                  
  Float_t b2L2J_dPhiPjjMet_minDPhiZMET;                                      
  Float_t b2L2J_dPhiMetISR;                                                  
  Float_t b2L2J_dPhiMetJet1;                                                 
  Float_t b2L2J_METOverHTLep;                                                
  Float_t b2L2J_mll;
  Float_t b2L2J_mlll;
  Float_t b2L2J_mllll;
  Float_t b2L2J_Rll;                                                         
  Float_t b2L2J_Ptll;                                                        
  Float_t b2L2J_absEtall;                                                    
  Float_t b2L2J_dPhiPllMet;                                                 
  Float_t b2L2J_dPhill;                                                                                                
  Float_t b2L2J_mt2leplsp_0;                                                 
  Double_t b2L2J_pileupWeight;                                               
  Double_t b2L2J_leptonWeight;                                               
  Double_t b2L2J_eventWeight;                                                
  Double_t b2L2J_bTagWeight;                                                 
  Double_t b2L2J_jvtWeight;                                                  
  Double_t b2L2J_globalDiLepTrigSF;                                          
  Double_t b2L2J_genWeight;                                                  
  Double_t b2L2J_genWeightUp;                                                
  Double_t b2L2J_genWeightDown;                                              
  Double_t b2L2J_truthMll;                                                   
  Double_t b2L2J_winoBinoMllWeight;                                          
  Double_t b2L2J_winoBinoXsecWeight;                                         
  Double_t b2L2J_winoBinoBrFracWeight;                                       
  Double_t b2L2J_winoBinoWeight;                                                           
  Double_t b2L2J_ttbarNNLOWeight;                                            
  Double_t b2L2J_ttbarNNLOWeightUp;                                          
  Double_t b2L2J_ttbarNNLOWeightDown;                                        
  Float_t b2L2J_truthTopPt;                                                  
  Float_t b2L2J_truthAntiTopPt;                                              
  Float_t b2L2J_truthTtbarPt;                                                
  Float_t b2L2J_truthTtbarM;                                                 
  Float_t b2L2J_x1;                                                          
  Float_t b2L2J_x2;                                                          
  Float_t b2L2J_pdf1;                                                        
  Float_t b2L2J_pdf2;                                                        
  Float_t b2L2J_scalePDF;                                                    
  Int_t b2L2J_id1;                                                           
  Int_t b2L2J_id2;                                                           
  Int_t b2L2J_nLeps_RJ;
  Int_t b2L2J_nJets_RJ;
  Int_t b2L2J_nBtagJets_RJ;
  vector<Float_t> b2L2J_jetPt_RJ;
  vector<Float_t> b2L2J_jetEta_RJ;
  vector<Float_t> b2L2J_jetPhi_RJ;
  vector<Float_t> b2L2J_jetM_RJ;
  vector<Float_t> b2L2J_lepPt_RJ;  
  vector<Float_t> b2L2J_lepEta_RJ; 
  vector<Float_t> b2L2J_lepPhi_RJ; 
  vector<Float_t> b2L2J_lepE_RJ;   
  vector<Float_t> b2L2J_lepsign_RJ;
  Bool_t b2L2J_is2Lep2Jet;
  Bool_t b2L2J_is2L2JInt;
  Bool_t b2L2J_is3Lep;
  Bool_t b2L2J_is3LInt;
  Bool_t b2L2J_is3Lep2Jet;
  Bool_t b2L2J_is3Lep3Jet;
  Bool_t b2L2J_is4Lep2Jet;
  Bool_t b2L2J_is4Lep3Jet;
  Float_t b2L2J_mll_RJ;
  Float_t b2L2J_H2PP;
  Float_t b2L2J_H5PP;
  Float_t b2L2J_RPT_HT5PP;                                       
  Float_t b2L2J_R_minH2P_minH3P;                                 
  Float_t b2L2J_R_H2PP_H5PP;                                
  Float_t b2L2J_dphiVP;                               
  Float_t b2L2J_minDphi;                              
  Float_t b2L2J_mTW;                             
  Float_t b2L2J_H4PP;                            
  Float_t b2L2J_RPT_HT4PP;                           
  Float_t b2L2J_R_HT4PP_H4PP;                          
  Float_t b2L2J_PTISR;                         
  Float_t b2L2J_RISR;                        
  Float_t b2L2J_PTI;                       
  Float_t b2L2J_dphiISRI;                      
  Float_t b2L2J_PTCM;                     
  Int_t b2L2J_NjS;                    
  Int_t b2L2J_NjISR;                   
  Float_t b2L2J_MZ;                  
  Float_t b2L2J_MJ;                 
  Float_t b2L2J_mTl3;                
  Float_t b2L2J_lept1Pt_VR;               
  Float_t b2L2J_lept1sign_VR;              
  Float_t b2L2J_lept2Pt_VR;             
  Float_t b2L2J_lept2sign_VR;            
  Float_t b2L2J_mll_RJ_VR;           
  Float_t b2L2J_H2PP_VR;          
  Float_t b2L2J_H5PP_VR;         
  Float_t b2L2J_RPT_HT5PP_VR;        
  Float_t b2L2J_R_minH2P_minH3P_VR;       
  Float_t b2L2J_R_H2PP_H5PP_VR;      
  Float_t b2L2J_dphiVP_VR;     
  Float_t b2L2J_PTISR_VR;    
  Float_t b2L2J_RISR_VR;   
  Float_t b2L2J_PTI_VR;  
  Float_t b2L2J_dphiISRI_VR; 
  Float_t b2L2J_PTCM_VR;
  Int_t   b2L2J_NjS_VR;
  Int_t   b2L2J_NjISR_VR;
  Float_t b2L2J_MZ_VR;
  Float_t b2L2J_MJ_VR;
  ULong64_t b2L2J_PRWHash;
  ULong64_t b2L2J_EventNumber;
  Float_t b2L2J_xsec;
  Float_t b2L2J_GenHt;
  Float_t b2L2J_GenMET;
  Int_t b2L2J_DatasetNumber;
  Int_t b2L2J_RunNumber;
  Int_t b2L2J_RandomRunNumber;
  Int_t b2L2J_FS; 

  Int_t b2L2J_isTT;    
  Int_t b2L2J_isTl;    
  Int_t b2L2J_islT;    
  Int_t b2L2J_isll;    
  Double_t b2L2J_r1;    
  Double_t b2L2J_r2;    
  Double_t b2L2J_f1;    
  Double_t b2L2J_f2;        

  Int_t b2L2J_nLep_base_OR;                                                   
  
  //virtual void InitializeOutput(TFile** file, TString filename,TTree** tree, TString treename );
  void ClearOutputBranches();

  void setSumOfMcWeights(double sumOfMcWeights);


  void WriteTree();

};
#endif // #ifdef make2L2JTrees_cxx
