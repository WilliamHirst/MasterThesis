#define HistFitterTree_cxx
#include "HistFitterTree.h"
#include "TParameter.h"

using namespace std;



HistFitterTree::HistFitterTree(TString MCID, TString syst, TString fileName, TString treeName) 
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

  tree->Branch("DatasetNumber",                                     &DatasetNumber); 	
  tree->Branch("RunNumber",                                         &RunNumber); 						                                              
  tree->Branch("EventNumber",                                       &EventNumber); 						                                              
  tree->Branch("lept1Pt",                                           &lept1Pt);					      
  tree->Branch("lept2Pt",                                           &lept2Pt);   						      
  tree->Branch("lept3Pt",                                           &lept3Pt);   						      
  tree->Branch("lept1Eta",                                          &lept1Eta);  						      
  tree->Branch("lept2Eta",  					  &lept2Eta);  						  
  tree->Branch("lept3Eta",                                          &lept3Eta);  						      
  tree->Branch("lept1Phi",                                          &lept1Phi);  						      
  tree->Branch("lept2Phi",                                          &lept2Phi);  						      
  tree->Branch("lept3Phi",                                          &lept3Phi);  						      
  tree->Branch("lept1E",                                            &lept1E);    						                                                  
  tree->Branch("lept2E",                                            &lept2E);    						      
  tree->Branch("lept3E",                                            &lept3E);    						      
  tree->Branch("lept1Flav",                                         &lept1Flav); 						      
  tree->Branch("lept2Flav",                                         &lept2Flav); 						      
  tree->Branch("lept3Flav",                                         &lept3Flav); 						      
  tree->Branch("metPhi",                                            &metPhi);    						      
  tree->Branch("MET",                                               &MET);       						      
  tree->Branch("METsig",                                            &METsig);    						      
  tree->Branch("nLeps",                                             &nLeps);     						      
  tree->Branch("nSigLeps",                                          &nSigLeps);     						      
  tree->Branch("eventWeight",                                       &eventWeight); 						        
  tree->Branch("WeightEvents", 					  &WeightEvents); 						  
  tree->Branch("WeightEventsEff",                                   &WeightEventsEff); 					            
  tree->Branch("WeightEventsSF",                                    &WeightEventsSF); 						           
  tree->Branch("WeightEventsbTag",                                  &WeightEventsbTag);					             
  tree->Branch("WeightEventsJVT",                                   &WeightEventsJVT); 					            
  tree->Branch("WeightEventselSF",                                  &WeightEventselSF); 					             
  tree->Branch("WeightEventsmuSF",                                  &WeightEventsmuSF); 					             
  tree->Branch("L2Mll",                                             &L2Mll);     						      
  tree->Branch("L2MT2",                                             &L2MT2);     						      
  tree->Branch("L2isEMU",                                           &L2isEMU);   						      
  tree->Branch("L2isEE",                                            &L2isEE);    						      
  tree->Branch("L2isMUMU",                                          &L2isMUMU);  						      
  tree->Branch("L2nCentralLightJet",                                &L2nCentralLightJet); 					               
  tree->Branch("L2nCentralBJet",                                    &L2nCentralBJet); 						           
  tree->Branch("L2nForwardJet",                                     &L2nForwardJet); 						          
  tree->Branch("L2nCentralLightJet30",                              &L2nCentralLightJet30); 					                 
  tree->Branch("L2nCentralBJet30",                                  &L2nCentralBJet30); 					             
  tree->Branch("L2nForwardJet30",                                   &L2nForwardJet30); 					            
  tree->Branch("L2nCentralLightJet40",                              &L2nCentralLightJet40); 					                 
  tree->Branch("L2nCentralBJet40",                                  &L2nCentralBJet40); 					             
  tree->Branch("L2nForwardJet40",                                   &L2nForwardJet40); 					            
  tree->Branch("L2nCentralLightJet50",                              &L2nCentralLightJet50); 					                 
  tree->Branch("L2nCentralBJet50",                                  &L2nCentralBJet50); 					             
  tree->Branch("L2nForwardJet50",                                   &L2nForwardJet50); 					            
  tree->Branch("L2nCentralLightJet60",                              &L2nCentralLightJet60); 					                 
  tree->Branch("L2nCentralBJet60",                                  &L2nCentralBJet60); 					             
  tree->Branch("L2nForwardJet60",                                   &L2nForwardJet60); 					            
  tree->Branch("L2SiLepTrigger",                                    &L2SiLepTrigger); 						           
  tree->Branch("L2DiLepTrigger",                                    &L2DiLepTrigger); 						           
  tree->Branch("jet1Pt",                                            &jet1Pt);    						      
  tree->Branch("jet2Pt",                                            &jet2Pt);    						      
  tree->Branch("jet3Pt",                                            &jet3Pt);    						      
  tree->Branch("jet1Eta",                                           &jet1Eta);   						      
  tree->Branch("jet2Eta",   					  &jet2Eta);   						  
  tree->Branch("jet3Eta",                                           &jet3Eta);   						      
  tree->Branch("jet1Phi",   					  &jet1Phi);   						  
  tree->Branch("jet2Phi",                                           &jet2Phi);   						      
  tree->Branch("jet3Phi",                                           &jet3Phi);   						      
  tree->Branch("mu",                                                &b_mu);        						      
  tree->Branch("syst_FT_EFF_B_down",                                &syst_FT_EFF_B_down); 					               
  tree->Branch("syst_FT_EFF_B_up",                                  &syst_FT_EFF_B_up); 					             
  tree->Branch("syst_FT_EFF_C_down",                                &syst_FT_EFF_C_down); 					               
  tree->Branch("syst_FT_EFF_C_up",                                  &syst_FT_EFF_C_up); 					             
  tree->Branch("syst_FT_EFF_Light_down",                            &syst_FT_EFF_Light_down); 					                   
  tree->Branch("syst_FT_EFF_Light_up",                              &syst_FT_EFF_Light_up); 					                 
  tree->Branch("syst_FT_EFF_extrapolation_down",                    &syst_FT_EFF_extrapolation_down); 				                           
  tree->Branch("syst_FT_EFF_extrapolation_up",                      &syst_FT_EFF_extrapolation_up); 				                         
  tree->Branch("syst_FT_EFF_extrapolation_from_charm_down",         &syst_FT_EFF_extrapolation_from_charm_down); 		                                      
  tree->Branch("syst_FT_EFF_extrapolation_from_charm_up",           &syst_FT_EFF_extrapolation_from_charm_up); 		                                    
  tree->Branch("syst_EL_EFF_ID_TOTAL_UncorrUncertainty_down",       &syst_EL_EFF_ID_TOTAL_UncorrUncertainty_down); 		                                        
  tree->Branch("syst_EL_EFF_ID_TOTAL_UncorrUncertainty_up",         &syst_EL_EFF_ID_TOTAL_UncorrUncertainty_up); 		                                      
  tree->Branch("syst_EL_EFF_Iso_TOTAL_UncorrUncertainty_down",      &syst_EL_EFF_Iso_TOTAL_UncorrUncertainty_down); 		                                         
  tree->Branch("syst_EL_EFF_Iso_TOTAL_UncorrUncertainty_up",        &syst_EL_EFF_Iso_TOTAL_UncorrUncertainty_up); 		                                       
  tree->Branch("syst_EL_EFF_Reco_TOTAL_UncorrUncertainty_down",     &syst_EL_EFF_Reco_TOTAL_UncorrUncertainty_down); 		                                          
  tree->Branch("syst_EL_EFF_Reco_TOTAL_UncorrUncertainty_up",       &syst_EL_EFF_Reco_TOTAL_UncorrUncertainty_up); 		                                        
  tree->Branch("syst_EL_EFF_Trigger_TOTAL_UncorrUncertainty_down",  &syst_EL_EFF_Trigger_TOTAL_UncorrUncertainty_down); 	                                             
  tree->Branch("syst_EL_EFF_Trigger_TOTAL_UncorrUncertainty_up",    &syst_EL_EFF_Trigger_TOTAL_UncorrUncertainty_up); 		                                           
  tree->Branch("syst_MUON_EFF_BADMUON_STAT_down",                   &syst_MUON_EFF_BADMUON_STAT_down); 			                            
  tree->Branch("syst_MUON_EFF_BADMUON_STAT_up",                     &syst_MUON_EFF_BADMUON_STAT_up); 				                          
  tree->Branch("syst_MUON_EFF_BADMUON_SYS_down",                    &syst_MUON_EFF_BADMUON_SYS_down); 				                           
  tree->Branch("syst_MUON_EFF_BADMUON_SYS_up",                      &syst_MUON_EFF_BADMUON_SYS_up); 				                         
  tree->Branch("syst_MUON_EFF_ISO_STAT_down",                       &syst_MUON_EFF_ISO_STAT_down); 				                        
  tree->Branch("syst_MUON_EFF_ISO_STAT_up",                         &syst_MUON_EFF_ISO_STAT_up); 				                      
  tree->Branch("syst_MUON_EFF_ISO_SYS_down",                        &syst_MUON_EFF_ISO_SYS_down); 				                       
  tree->Branch("syst_MUON_EFF_ISO_SYS_up",                          &syst_MUON_EFF_ISO_SYS_up); 				                     
  tree->Branch("syst_MUON_EFF_RECO_STAT_down",                      &syst_MUON_EFF_RECO_STAT_down); 				                         
  tree->Branch("syst_MUON_EFF_RECO_STAT_up",                        &syst_MUON_EFF_RECO_STAT_up); 				                       
  tree->Branch("syst_MUON_EFF_RECO_STAT_LOWPT_down",                &syst_MUON_EFF_RECO_STAT_LOWPT_down); 			                               
  tree->Branch("syst_MUON_EFF_RECO_STAT_LOWPT_up",                  &syst_MUON_EFF_RECO_STAT_LOWPT_up); 			                             
  tree->Branch("syst_MUON_EFF_RECO_SYS_down",                       &syst_MUON_EFF_RECO_SYS_down); 				                        
  tree->Branch("syst_MUON_EFF_RECO_SYS_up",                         &syst_MUON_EFF_RECO_SYS_up); 				                      
  tree->Branch("syst_MUON_EFF_RECO_SYS_LOWPT_down",                 &syst_MUON_EFF_RECO_SYS_LOWPT_down); 			                              
  tree->Branch("syst_MUON_EFF_RECO_SYS_LOWPT_up",                   &syst_MUON_EFF_RECO_SYS_LOWPT_up); 			                            
  tree->Branch("syst_MUON_EFF_TTVA_STAT_down",                      &syst_MUON_EFF_TTVA_STAT_down); 				                         
  tree->Branch("syst_MUON_EFF_TTVA_STAT_up",                        &syst_MUON_EFF_TTVA_STAT_up); 				                       
  tree->Branch("syst_MUON_EFF_TTVA_SYS_down",                       &syst_MUON_EFF_TTVA_SYS_down); 				                        
  tree->Branch("syst_MUON_EFF_TTVA_SYS_up",                         &syst_MUON_EFF_TTVA_SYS_up); 				                      
  tree->Branch("syst_MUON_EFF_TrigStatUncertainty_down",            &syst_MUON_EFF_TrigStatUncertainty_down); 			                                   
  tree->Branch("syst_MUON_EFF_TrigStatUncertainty_up",              &syst_MUON_EFF_TrigStatUncertainty_up); 			                                 
  tree->Branch("syst_MUON_EFF_TrigSystUncertainty_down",            &syst_MUON_EFF_TrigSystUncertainty_down); 			                                   
  tree->Branch("syst_MUON_EFF_TrigSystUncertainty_up",              &syst_MUON_EFF_TrigSystUncertainty_up); 			                                 
  tree->Branch("syst_JET_JvtEfficiency_1down",                      &syst_JET_JvtEfficiency_1down); 				                         
  tree->Branch("syst_JET_JvtEfficiency_1up",                        &syst_JET_JvtEfficiency_1up); 				                       
  tree->Branch("syst_JET_fJvtEfficiency_1down",                     &syst_JET_fJvtEfficiency_1down); 				                          
  tree->Branch("syst_JET_fJvtEfficiency_1up",                       &syst_JET_fJvtEfficiency_1up);
  tree->Branch("syst_MM_down",                                      &syst_MM_down); 					               
  tree->Branch("syst_MM_up",                                        &syst_MM_up); 
  tree->Branch("MM_key",                                            &MM_key); 


  ClearOutputBranches();
}


void HistFitterTree::ClearOutputBranches(void) 
{

  MM_key.clear();
  eventWeight.clear();
  syst_MM_down.clear();
  syst_MM_up.clear();
  DatasetNumber = 0;
  RunNumber = 0;
  EventNumber = 0;
  lept1Pt     = 0;
  lept2Pt     = 0;
  lept3Pt   = 0;
  lept1Eta  = 0;
  lept2Eta  = 0;
  lept3Eta  = 0;
  lept1Phi  = 0;
  lept2Phi  = 0;
  lept3Phi  = 0;
  lept1E    = 0;
  lept2E    = 0;
  lept3E    = 0;
  lept1Flav = 0;
  lept2Flav = 0;
  lept3Flav = 0;
  metPhi    = 0;
  MET       = 0;
  METsig      = 0;
  nLeps     				  = 0;
  nSigLeps     				  = 0;
  //eventWeight 				  = 0;
  WeightEvents 				  = 0;
  WeightEventsEff 			  = 0;
  WeightEventsSF 				  = 0;
  WeightEventsbTag 			  = 0;
  WeightEventsJVT 			  = 0;
  WeightEventselSF 			  = 0;
  WeightEventsmuSF 			  = 0;
  L2Mll     				  = 0;
  L2MT2     				  = 0;
  L2isEMU   				  = 0;
  L2isEE    				  = 0;
  L2isMUMU  				  = 0;
  L2nCentralLightJet                        = 0;
  L2nCentralBJet 				  = 0;
  L2nForwardJet 				  = 0;
  L2nCentralLightJet30 			  = 0;
  L2nCentralBJet30 			  = 0;
  L2nForwardJet30 			  = 0;
  L2nCentralLightJet40 			  = 0;
  L2nCentralBJet40 			  = 0;
  L2nForwardJet40 			  = 0;
  L2nCentralLightJet50 			  = 0;
  L2nCentralBJet50 			  = 0;
  L2nForwardJet50 			  = 0;
  L2nCentralLightJet60 			  = 0;
  L2nCentralBJet60 			  = 0;
  L2nForwardJet60 			  = 0;
  L2SiLepTrigger                            = 0;
  L2DiLepTrigger 				  = 0;
  jet1Pt    				  = 0;
  jet2Pt    				  = 0;
  jet3Pt    				  = 0;
  jet1Eta   				  = 0;
  jet2Eta   				  = 0;
  jet3Eta   				  = 0;
  jet1Phi   				  = 0;
  jet2Phi   				  = 0;
  jet3Phi   				  = 0;
  b_mu        				  = 0;
  //syst_MM_down                            = 0;
  //syst_MM_up                              = 0;
  syst_FT_EFF_B_down 			  = 0;
  syst_FT_EFF_B_up 			  = 0;
  syst_FT_EFF_C_down 			  = 0;
  syst_FT_EFF_C_up    = 0;
  syst_FT_EFF_Light_down 				  = 0;
  syst_FT_EFF_Light_up 				  = 0;
  syst_FT_EFF_extrapolation_down 			  = 0;
  syst_FT_EFF_extrapolation_up 			  = 0;
  syst_FT_EFF_extrapolation_from_charm_down 	  = 0;
  syst_FT_EFF_extrapolation_from_charm_up 	  = 0;
  syst_EL_EFF_ID_TOTAL_UncorrUncertainty_down 	  = 0;
  syst_EL_EFF_ID_TOTAL_UncorrUncertainty_up 	  = 0;
  syst_EL_EFF_Iso_TOTAL_UncorrUncertainty_down 	  = 0;
  syst_EL_EFF_Iso_TOTAL_UncorrUncertainty_up 	  = 0;
  syst_EL_EFF_Reco_TOTAL_UncorrUncertainty_down 	  = 0;
  syst_EL_EFF_Reco_TOTAL_UncorrUncertainty_up 	  = 0;
  syst_EL_EFF_Trigger_TOTAL_UncorrUncertainty_down  = 0;
  syst_EL_EFF_Trigger_TOTAL_UncorrUncertainty_up 	  = 0;
  syst_MUON_EFF_BADMUON_STAT_down                   = 0;
  syst_MUON_EFF_BADMUON_STAT_up 			  = 0;
  syst_MUON_EFF_BADMUON_SYS_down 			  = 0;
  syst_MUON_EFF_BADMUON_SYS_up 			  = 0;
  syst_MUON_EFF_ISO_STAT_down 			  = 0;
  syst_MUON_EFF_ISO_STAT_up 			  = 0;
  syst_MUON_EFF_ISO_SYS_down 			  = 0;
  syst_MUON_EFF_ISO_SYS_up 			  = 0;
  syst_MUON_EFF_RECO_STAT_down 			  = 0;
  syst_MUON_EFF_RECO_STAT_up 			  = 0;
  syst_MUON_EFF_RECO_STAT_LOWPT_down 		  = 0;
  syst_MUON_EFF_RECO_STAT_LOWPT_up 		  = 0;
  syst_MUON_EFF_RECO_SYS_down 			  = 0;
  syst_MUON_EFF_RECO_SYS_up 			  = 0;
  syst_MUON_EFF_RECO_SYS_LOWPT_down 		  = 0;
  syst_MUON_EFF_RECO_SYS_LOWPT_up                 = 0; 
  syst_MUON_EFF_TTVA_STAT_down 			  = 0;
  syst_MUON_EFF_TTVA_STAT_up 			  = 0;
  syst_MUON_EFF_TTVA_SYS_down 			  = 0;
  syst_MUON_EFF_TTVA_SYS_up 			  = 0;
  syst_MUON_EFF_TrigStatUncertainty_down 		  = 0;
  syst_MUON_EFF_TrigStatUncertainty_up 		  = 0;
  syst_MUON_EFF_TrigSystUncertainty_down 		  = 0;
  syst_MUON_EFF_TrigSystUncertainty_up 		  = 0;
  syst_JET_JvtEfficiency_1down 			  = 0;
  syst_JET_JvtEfficiency_1up 			  = 0;
  syst_JET_fJvtEfficiency_1down 			  = 0;
  syst_JET_fJvtEfficiency_1up 			  = 0;

}





//-----------------------------------------------------------------------------------------------------------
HistFitterTree::~HistFitterTree()
{
    // Write out the output tree and close the output file
  file->Write();
  file->Close();
  delete file;
}


void HistFitterTree::WriteTree()
{
  //file->cd();
    tree->Fill();
    //tree->Write();
    //file->Write();
    //file->Close();
    ClearOutputBranches();
}




void HistFitterTree::setSumOfMcWeights(double sumOfMcWeights) 
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



