#define makeMYTree_cxx
#include "makeMYTree.h"
#include "TParameter.h"
#include "TString.h"

using namespace std;



makeMYTree::makeMYTree(TString MCID, TString syst, TString fileName, TString treeName)
{
  if(fileName==""){
    fileName = syst+"_"+MCID+".root";
  }
  if(treeName==""){
    treeName = "id_" + MCID;
  }
  // Setup a TTree in a output file
  file = TFile::Open(fileName, "RECREATE");
  file->cd();
  //tree = new TTree("id_"+MCID, "id_"+MCID);
  tree = new TTree(treeName, treeName);
  tree->SetAutoSave(10000000);
  TTree::SetBranchStyle(1);
  tree->SetDirectory(file);

  mymcid = MCID;


  tree->Branch("trigMatch_1L2LTrig",                                 &bMY_trigMatch_1L2LTrig);                                   
  tree->Branch("trigMatch_1LTrig",                                   &bMY_trigMatch_1LTrig);                                 
  tree->Branch("trigMatch_1L2LTrigOR",                               &bMY_trigMatch_1L2LTrigOR);                             
  tree->Branch("trigMatch_1LTrigOR",                                 &bMY_trigMatch_1LTrigOR);                               
  tree->Branch("trigMatch_2LTrig",                                   &bMY_trigMatch_2LTrig);                                 
  tree->Branch("trigMatch_2LTrigOR",                                 &bMY_trigMatch_2LTrigOR);                               
  tree->Branch("trigMatch_HLT_e24_lhmedium_L1EM20VH",                &bMY_trigMatch_HLT_e24_lhmedium_L1EM20VH);              
  tree->Branch("trigMatch_HLT_e60_lhmedium",                         &bMY_trigMatch_HLT_e60_lhmedium);                       
  tree->Branch("trigMatch_HLT_e120_lhloose",                         &bMY_trigMatch_HLT_e120_lhloose);                       
  tree->Branch("trigMatch_HLT_mu20_iloose_L1MU15",                   &bMY_trigMatch_HLT_mu20_iloose_L1MU15);                 
  tree->Branch("trigMatch_HLT_2e12_lhloose_L12EM10VH",               &bMY_trigMatch_HLT_2e12_lhloose_L12EM10VH);             
  tree->Branch("trigMatch_HLT_mu18_mu8noL1",                         &bMY_trigMatch_HLT_mu18_mu8noL1);                       
  tree->Branch("trigMatch_HLT_e17_lhloose_mu14",                     &bMY_trigMatch_HLT_e17_lhloose_mu14);                   
  tree->Branch("trigMatch_HLT_e7_lhmedium_mu24",                     &bMY_trigMatch_HLT_e7_lhmedium_mu24);                   
  tree->Branch("trigMatch_HLT_e24_lhtight_nod0_ivarloose",           &bMY_trigMatch_HLT_e24_lhtight_nod0_ivarloose);         
  tree->Branch("trigMatch_HLT_e24_lhmedium_nod0_L1EM20VH",           &bMY_trigMatch_HLT_e24_lhmedium_nod0_L1EM20VH);         
  tree->Branch("trigMatch_HLT_e60_medium",                           &bMY_trigMatch_HLT_e60_medium);                         
  tree->Branch("trigMatch_HLT_mu40",                                 &bMY_trigMatch_HLT_mu40);                               
  tree->Branch("trigMatch_HLT_mu24_iloose_L1MU15",                   &bMY_trigMatch_HLT_mu24_iloose_L1MU15);                 
  tree->Branch("trigMatch_HLT_mu24_ivarloose_L1MU15",                &bMY_trigMatch_HLT_mu24_ivarloose_L1MU15);              
  tree->Branch("trigMatch_HLT_mu24_ivarmedium",                      &bMY_trigMatch_HLT_mu24_ivarmedium);                    
  tree->Branch("trigMatch_HLT_mu24_imedium",                         &bMY_trigMatch_HLT_mu24_imedium);                       
  tree->Branch("trigMatch_HLT_mu26_imedium",                         &bMY_trigMatch_HLT_mu26_imedium);                       
  tree->Branch("trigMatch_HLT_2e15_lhvloose_nod0_L12EM13VH",         &bMY_trigMatch_HLT_2e15_lhvloose_nod0_L12EM13VH);       
  tree->Branch("trigMatch_HLT_2e17_lhvloose_nod0",                   &bMY_trigMatch_HLT_2e17_lhvloose_nod0);                 
  tree->Branch("trigMatch_HLT_2mu10",                                &bMY_trigMatch_HLT_2mu10);                              
  tree->Branch("trigMatch_HLT_2mu14",                                &bMY_trigMatch_HLT_2mu14);                              
  tree->Branch("trigMatch_HLT_mu20_mu8noL1",                         &bMY_trigMatch_HLT_mu20_mu8noL1);                       
  tree->Branch("trigMatch_HLT_mu22_mu8noL1",                         &bMY_trigMatch_HLT_mu22_mu8noL1);                       
  tree->Branch("trigMatch_HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1",  &bMY_trigMatch_HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1);
  tree->Branch("trigMatch_HLT_e26_lhtight_nod0_ivarloose",           &bMY_trigMatch_HLT_e26_lhtight_nod0_ivarloose);         
  tree->Branch("trigMatch_HLT_e26_lhtight_nod0",                     &bMY_trigMatch_HLT_e26_lhtight_nod0);                   
  tree->Branch("trigMatch_HLT_e60_lhmedium_nod0",                    &bMY_trigMatch_HLT_e60_lhmedium_nod0);                  
  tree->Branch("trigMatch_HLT_e140_lhloose_nod0",                    &bMY_trigMatch_HLT_e140_lhloose_nod0);                  
  tree->Branch("trigMatch_HLT_e300_etcut",                           &bMY_trigMatch_HLT_e300_etcut);                         
  tree->Branch("trigMatch_HLT_mu26_ivarmedium",                      &bMY_trigMatch_HLT_mu26_ivarmedium);                    
  tree->Branch("trigMatch_HLT_mu50",                                 &bMY_trigMatch_HLT_mu50);                               
  tree->Branch("trigMatch_HLT_mu60_0eta105_msonly",                  &bMY_trigMatch_HLT_mu60_0eta105_msonly);                
  tree->Branch("trigMatch_HLT_2e17_lhvloose_nod0_L12EM15VHI",        &bMY_trigMatch_HLT_2e17_lhvloose_nod0_L12EM15VHI);      
  tree->Branch("trigMatch_HLT_2e24_lhvloose_nod0",                   &bMY_trigMatch_HLT_2e24_lhvloose_nod0);                 
  tree->Branch("trigMatch_HLT_e17_lhloose_nod0_mu14",                &bMY_trigMatch_HLT_e17_lhloose_nod0_mu14);              
  tree->Branch("trigMatch_HLT_e26_lhmedium_nod0_mu8noL1",            &bMY_trigMatch_HLT_e26_lhmedium_nod0_mu8noL1);          
  tree->Branch("trigMatch_HLT_e7_lhmedium_nod0_mu24",                &bMY_trigMatch_HLT_e7_lhmedium_nod0_mu24);              
  tree->Branch("mu",                                                 &bMY_mu);                                               
  tree->Branch("avg_mu",                                             &bMY_avg_mu);                                           
  tree->Branch("actual_mu",                                          &bMY_actual_mu);                                        
  tree->Branch("nVtx",                                               &bMY_nVtx);                                             
  tree->Branch("channel",                                            &bMY_channel);                                          
  tree->Branch("nLep_base",                                          &bMY_nLep_base);                                        
  tree->Branch("nLep_signal",                                        &bMY_nLep_signal);                                      
  tree->Branch("lepFlavor",                                          &bMY_lepFlavor);                                        
  tree->Branch("lepCharge",                                          &bMY_lepCharge);                                        
  tree->Branch("lepAuthor",                                          &bMY_lepAuthor);                                        
  tree->Branch("lepPt",                                              &bMY_lepPt);                                            
  tree->Branch("lepEta",                                             &bMY_lepEta);                                           
  tree->Branch("lepPhi",                                             &bMY_lepPhi);                                           
  tree->Branch("lepM",                                               &bMY_lepM);                                             
  tree->Branch("lepD0",                                              &bMY_lepD0);                                            
  tree->Branch("lepD0Sig",                                           &bMY_lepD0Sig);                                         
  tree->Branch("lepZ0",                                              &bMY_lepZ0);                                            
  tree->Branch("lepZ0SinTheta",                                      &bMY_lepZ0SinTheta);                                    
  tree->Branch("lepTopoetcone20",                                    &bMY_lepTopoetcone20);                                  
  tree->Branch("lepTopoetcone30",                                    &bMY_lepTopoetcone30);                                  
  tree->Branch("lepPtvarcone20",                                     &bMY_lepPtvarcone20);                                   
  tree->Branch("lepPtvarcone30",                                     &bMY_lepPtvarcone30);                                   
  tree->Branch("lepCorrTopoetcone20",                                &bMY_lepCorrTopoetcone20);                              
  tree->Branch("lepCorrPtvarcone20",                                 &bMY_lepCorrPtvarcone20);                               
  tree->Branch("lepCorrPtvarcone30",                                 &bMY_lepCorrPtvarcone30);                               
  tree->Branch("lepPassOR",                                          &bMY_lepPassOR);                                        
  tree->Branch("lepType",                                            &bMY_lepType);                                          
  tree->Branch("lepOrigin",                                          &bMY_lepOrigin);                                        
  tree->Branch("lepEgMotherType",                                    &bMY_lepEgMotherType);                                  
  tree->Branch("lepEgMotherOrigin",                                  &bMY_lepEgMotherOrigin);                                
  tree->Branch("lepEgMotherPdgId",                                   &bMY_lepEgMotherPdgId);                                 
  tree->Branch("lepECIDS",                                           &bMY_lepECIDS);                                         
  tree->Branch("lepIsHF",                                            &bMY_lepIsHF);                                          
  tree->Branch("lepIsLF",                                            &bMY_lepIsLF);                                          
  tree->Branch("lepIsCO",                                            &bMY_lepIsCO);                                          
  tree->Branch("lepIsCF",                                            &bMY_lepIsCF);                                          
  tree->Branch("lepIsUK",                                            &bMY_lepIsUK);                                          
  tree->Branch("lepIsPR",                                            &bMY_lepIsPR);                                          
  tree->Branch("lepPassBL",                                          &bMY_lepPassBL);                                        
  tree->Branch("lepVeryLoose",                                       &bMY_lepVeryLoose);                                     
  tree->Branch("lepLoose",                                           &bMY_lepLoose);                                         
  tree->Branch("lepMedium",                                          &bMY_lepMedium);                                        
  tree->Branch("lepTight",                                           &bMY_lepTight);                                         
  tree->Branch("lepIsoFCHighPtCaloOnly",                             &bMY_lepIsoFCHighPtCaloOnly);
  tree->Branch("lepIsoFCLoose",                                      &bMY_lepIsoFCLoose);
  tree->Branch("lepIsoFCTight",                                      &bMY_lepIsoFCTight);
  tree->Branch("lepIsoFCLoose_FixedRad",                             &bMY_lepIsoFCLoose_FixedRad);
  tree->Branch("lepIsoFCTight_FixedRad",                             &bMY_lepIsoFCTight_FixedRad);
  tree->Branch("lepIsoHighPtCaloOnly",                               &bMY_lepIsoHighPtCaloOnly);
  tree->Branch("lepIsoTightTrackOnly_VarRad",                        &bMY_lepIsoTightTrackOnly_VarRad);
  tree->Branch("lepIsoTightTrackOnly_FixedRad",                      &bMY_lepIsoTightTrackOnly_FixedRad);
  tree->Branch("lepIsoLoose_VarRad",                                 &bMY_lepIsoLoose_VarRad);
  tree->Branch("lepIsoTight_VarRad",                                 &bMY_lepIsoTight_VarRad);
  tree->Branch("lepIsoPLVLoose",                                     &bMY_lepIsoPLVLoose);
  tree->Branch("lepIsoPLVTight",                                     &bMY_lepIsoPLVTight);                
  tree->Branch("lepTruthMatched",                                    &bMY_lepTruthMatched);                                  
  tree->Branch("lepTruthCharge",                                     &bMY_lepTruthCharge);                                   
  tree->Branch("lepTruthPt",                                         &bMY_lepTruthPt);                                       
  tree->Branch("lepTruthEta",                                        &bMY_lepTruthEta);                                      
  tree->Branch("lepTruthPhi",                                        &bMY_lepTruthPhi);                                      
  tree->Branch("lepTruthM",                                          &bMY_lepTruthM);                                        
  tree->Branch("nJet30",                                             &bMY_nJet30);                                           
  tree->Branch("nJet20",                                             &bMY_nJet20);                                           
  tree->Branch("nBJet20_MV2c10_FixedCutBEff_77",                     &bMY_nBJet20_MV2c10_FixedCutBEff_77);                   
  tree->Branch("jetPt",                                              &bMY_jetPt);                                            
  tree->Branch("jetEta",                                             &bMY_jetEta);                                           
  tree->Branch("jetPhi",                                             &bMY_jetPhi);                                           
  tree->Branch("jetM",                                               &bMY_jetM);                                             
  tree->Branch("jetJVT",                                             &bMY_jetJVT);                                           
  tree->Branch("jetPassOR",                                          &bMY_jetPassOR);                                        
  tree->Branch("jetSignal",                                          &bMY_jetSignal);                                        
  tree->Branch("mjj",                                                &bMY_mjj);                                              
  tree->Branch("jetTileEnergy",                                      &bMY_jetTileEnergy);                                    
  tree->Branch("jetMV2c10",                                          &bMY_jetMV2c10);                                        
  tree->Branch("met_Et",                                             &bMY_met_Et);                                           
  tree->Branch("met_Sign",                                           &bMY_met_Sign);                                         
  tree->Branch("met_Phi",                                            &bMY_met_Phi);                                          
  tree->Branch("met_Et_loose",                                       &bMY_met_Et_loose);                                     
  tree->Branch("met_Et_tighter",                                     &bMY_met_Et_tighter);                                   
  tree->Branch("met_Et_tenacious",                                   &bMY_met_Et_tenacious);                                 
  tree->Branch("met_Phi_loose",                                      &bMY_met_Phi_loose);                                    
  tree->Branch("met_Phi_tighter",                                    &bMY_met_Phi_tighter);                                  
  tree->Branch("met_Phi_tenacious",                                  &bMY_met_Phi_tenacious);                                
  tree->Branch("mll",                                                &bMY_mll);                                              
  tree->Branch("pileupWeight",                                       &bMY_pileupWeight);                                     
  tree->Branch("leptonWeight",                                       &bMY_leptonWeight);                                     
  tree->Branch("eventWeight",                                        &bMY_eventWeight);                                      
  tree->Branch("bTagWeight",                                         &bMY_bTagWeight);                                       
  tree->Branch("jvtWeight",                                          &bMY_jvtWeight);                                        
  tree->Branch("globalDiLepTrigSF",                                  &bMY_globalDiLepTrigSF);                                
  tree->Branch("flavSymWeight",                                      &bMY_flavSymWeight);                                    
  tree->Branch("genWeight",                                          &bMY_genWeight);                                        
  tree->Branch("genWeightUp",                                        &bMY_genWeightUp);                                      
  tree->Branch("genWeightDown",                                      &bMY_genWeightDown);                                    
  tree->Branch("PRWHash",                                            &bMY_PRWHash);                                          
  tree->Branch("EventNumber",                                        &bMY_EventNumber);                                      
  tree->Branch("xsec",                                               &bMY_xsec);                                             
  tree->Branch("GenHt",                                              &bMY_GenHt);                                            
  tree->Branch("GenMET",                                             &bMY_GenMET);                                           
  tree->Branch("DatasetNumber",                                      &bMY_DatasetNumber);                                    
  tree->Branch("RunNumber",                                          &bMY_RunNumber);                                        
  tree->Branch("RandomRunNumber",                                    &bMY_RandomRunNumber);                                  
  tree->Branch("FS",                                                 &bMY_FS);                    

  tree->Branch("MM_weight",                                         &bMY_MM_weight); 	             
  tree->Branch("syst_MM_down",                                      &bMY_syst_MM_down); 					               
  tree->Branch("syst_MM_up",                                        &bMY_syst_MM_up); 
  tree->Branch("MM_key",                                            &bMY_MM_key); 

  ClearOutputBranches();
}


void makeMYTree::ClearOutputBranches(void)
{

  bMY_MM_key.clear();
  bMY_MM_weight.clear();
  bMY_syst_MM_down.clear();
  bMY_syst_MM_up.clear();
  bMY_lepFlavor.clear();                                              
  bMY_lepCharge.clear();                                              
  bMY_lepAuthor.clear();                                              
  bMY_lepPt.clear();                                                  
  bMY_lepEta.clear();                                                 
  bMY_lepPhi.clear();                                                 
  bMY_lepM.clear();                                                   
  bMY_lepD0.clear();                                                  
  bMY_lepD0Sig.clear();                                               
  bMY_lepZ0.clear();                                                  
  bMY_lepZ0SinTheta.clear();                                          
  bMY_lepTopoetcone20.clear();                                        
  bMY_lepTopoetcone30.clear();                                        
  bMY_lepPtvarcone20.clear();                                         
  bMY_lepPtvarcone30.clear();                                         
  bMY_lepCorrTopoetcone20.clear();                                    
  bMY_lepCorrPtvarcone20.clear();                                     
  bMY_lepCorrPtvarcone30.clear();                                     
  bMY_lepPassOR.clear();                                           
  bMY_lepType.clear();                                                
  bMY_lepOrigin.clear();                                              
  bMY_lepEgMotherType.clear();                                        
  bMY_lepEgMotherOrigin.clear();                                      
  bMY_lepEgMotherPdgId.clear();                                       
  bMY_lepECIDS.clear();                                               
  bMY_lepIsHF.clear();                                                
  bMY_lepIsLF.clear();                                                
  bMY_lepIsCO.clear();                                                
  bMY_lepIsCF.clear();                                                
  bMY_lepIsUK.clear();                                                
  bMY_lepIsPR.clear();                                                
  bMY_lepPassBL.clear();                                              
  bMY_lepVeryLoose.clear();                                           
  bMY_lepLoose.clear();                                               
  bMY_lepMedium.clear();                                              
  bMY_lepTight.clear();                                               
  bMY_lepIsoFCHighPtCaloOnly.clear();	    
  bMY_lepIsoFCLoose.clear();		            
  bMY_lepIsoFCTight.clear();		            
  bMY_lepIsoFCLoose_FixedRad.clear();	    
  bMY_lepIsoFCTight_FixedRad.clear();	    
  bMY_lepIsoHighPtCaloOnly.clear();	            
  bMY_lepIsoTightTrackOnly_VarRad.clear();	    
  bMY_lepIsoTightTrackOnly_FixedRad.clear();     
  bMY_lepIsoLoose_VarRad.clear();		    
  bMY_lepIsoTight_VarRad.clear();		    
  bMY_lepIsoPLVLoose.clear();		    
  bMY_lepIsoPLVTight.clear();                    
  bMY_lepTruthMatched.clear();                                        
  bMY_lepTruthCharge.clear();                                         
  bMY_lepTruthPt.clear();                                             
  bMY_lepTruthEta.clear();                                            
  bMY_lepTruthPhi.clear();                                           
  bMY_lepTruthM.clear();                                                     
  bMY_jetPt.clear();                                                  
  bMY_jetEta.clear();                                                 
  bMY_jetPhi.clear();                                                 
  bMY_jetM.clear();                                                 
  bMY_jetJVT.clear();                                                 
  bMY_jetPassOR.clear();                                              
  bMY_jetSignal.clear();                                              
  bMY_jetTileEnergy.clear();                                          
  bMY_jetMV2c10.clear();                                              


  
  return;
}





//-----------------------------------------------------------------------------------------------------------
makeMYTree::~makeMYTree()
{
    // Write out the output tree and close the output file
  file->Write();
  file->Close();
  delete file;
}


void makeMYTree::WriteTree()
{
  //file->cd();
    tree->Fill();
    //tree->Write();
    //file->Write();
    //file->Close();
    ClearOutputBranches();
}




void makeMYTree::setSumOfMcWeights(double sumOfMcWeights)
{
    // Define histogram
    TH1D *sumwhist = new TH1D("sumOfMcWeights_"+mymcid,"sumOfMcWeights_"+mymcid,1,0.,1.);

    // Fill histogram
    sumwhist -> Fill( 0. , sumOfMcWeights ) ;

    // Write intLumi to file
    file->cd();
    sumwhist->SetDirectory(file);
    sumwhist->Write();
    sumwhist->SetDirectory(0);

    delete sumwhist;
}



