#define make2L2JTree_cxx
#include "make2L2JTree.h"
#include "TParameter.h"

using namespace std;



make2L2JTree::make2L2JTree(TString MCID, TString syst, TString fileName, TString treeName, std::vector<TString> unc_keys)
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

  //std::vector<double> b2L0J_FNPweight;
  for(std::vector< TString >::iterator itLT = unc_keys.begin(); itLT != unc_keys.end(); itLT++){
    b2L2J_FNPweight[(*itLT)] = 0.0;
    if((*itLT).EqualTo("NOM")){
      tree->Branch("FNP_WEIGHTS",&b2L2J_FNPweight[*(itLT)]);
    }else{
      tree->Branch("FNP_"+(*itLT),&b2L2J_FNPweight[*(itLT)]);
    }
  }
  tree->Branch("FNP_TOTAL_UP",&b2L2J_FNP_TOTAL_UP);
  tree->Branch("FNP_TOTAL_DOWN",&b2L2J_FNP_TOTAL_DOWN);
  
  
  tree->Branch("trigWeight_2LTrig",                                     			  &b2L2J_trigWeight_2LTrig);                                          
  tree->Branch("trigMatch_2LTrig",                                        			  &b2L2J_trigMatch_2LTrig);                                             
  tree->Branch("mu",                                                     			  &b2L2J_mu);                                                          
  tree->Branch("avg_mu",                                                 			  &b2L2J_avg_mu);                                                      
  tree->Branch("actual_mu",                                              			  &b2L2J_actual_mu);                                                   
  tree->Branch("nVtx",                                                     			  &b2L2J_nVtx);                                                          
  tree->Branch("channel",                                                  			  &b2L2J_channel);                                                       
  tree->Branch("nLep_combi",                                                			  &b2L2J_nLep_combi);                                                     
  tree->Branch("nLep_base",                                                			  &b2L2J_nLep_base);                                                     
  tree->Branch("nLep_signal",                                              			  &b2L2J_nLep_signal);                                                   
  tree->Branch("lepFlavor",                                          			  &b2L2J_lepFlavor);                                               
  tree->Branch("lepCharge",                                          			  &b2L2J_lepCharge);                                               
  tree->Branch("lepAuthor",                                          			  &b2L2J_lepAuthor);                                               
  tree->Branch("lepPt",                                            			  &b2L2J_lepPt);                                                 
  tree->Branch("lepEta",                                           			  &b2L2J_lepEta);                                                
  tree->Branch("lepPhi",                                           			  &b2L2J_lepPhi);                                                
  tree->Branch("lepM",                                             			  &b2L2J_lepM);                                                  
  tree->Branch("lepD0",                                            			  &b2L2J_lepD0);                                                 
  tree->Branch("lepD0Sig",                                         			  &b2L2J_lepD0Sig);                                              
  tree->Branch("lepZ0",                                            			  &b2L2J_lepZ0);                                                 
  tree->Branch("lepZ0SinTheta",                                    			  &b2L2J_lepZ0SinTheta);                                                                   
  tree->Branch("lepPassOR",                                         			  &b2L2J_lepPassOR);                                              
  tree->Branch("lepType",                                            			  &b2L2J_lepType);                                                 
  tree->Branch("lepOrigin",                                          			  &b2L2J_lepOrigin);
  tree->Branch("lepIFFClass",                                           		  &b2L2J_lepIFFClass);       
  tree->Branch("lepEgMotherType",                                    			  &b2L2J_lepEgMotherType);                                         
  tree->Branch("lepEgMotherOrigin",                                  			  &b2L2J_lepEgMotherOrigin);                                       
  tree->Branch("lepEgMotherPdgId",                                   			  &b2L2J_lepEgMotherPdgId);                                        
  tree->Branch("lepECIDS",                                           			  &b2L2J_lepECIDS);                                                
  tree->Branch("lepNPix",                                            			  &b2L2J_lepNPix);                                                                             
  tree->Branch("lepPassBL",                                         			  &b2L2J_lepPassBL);                                              
  tree->Branch("lepSignal",                                         			  &b2L2J_lepSignal);                                              
  tree->Branch("lepTruthMatched",                                   			  &b2L2J_lepTruthMatched);                                        
  tree->Branch("lepTruthCharge",                                     			  &b2L2J_lepTruthCharge);                                          
  tree->Branch("lepTruthPt",                                       			  &b2L2J_lepTruthPt);                                            
  tree->Branch("lepTruthEta",                                      			  &b2L2J_lepTruthEta);                                           
  tree->Branch("lepTruthPhi",                                      			  &b2L2J_lepTruthPhi);                                           
  tree->Branch("lepTruthM",                                        			  &b2L2J_lepTruthM);                                             
  tree->Branch("nJet30",                                                   		  &b2L2J_nJet30);                                                        
  tree->Branch("nJet20",                                                   		  &b2L2J_nJet20);                                                                                                     
  tree->Branch("nBJet20_MV2c10_FixedCutBEff_77",                           		  &b2L2J_nBJet20_MV2c10_FixedCutBEff_77);                                  
  tree->Branch("jetPt",                                            			  &b2L2J_jetPt);                                                 
  tree->Branch("jetEta",                                           			  &b2L2J_jetEta);                                                
  tree->Branch("jetPhi",                                           			  &b2L2J_jetPhi);                                                
  tree->Branch("jetM",                                             			  &b2L2J_jetM);  
  tree->Branch("jetTileEnergy",                                    			  &b2L2J_jetTileEnergy);                                         
  tree->Branch("jetMV2c10",                                        			  &b2L2J_jetMV2c10);
  
  tree->Branch("mjj",                                                    			  &b2L2J_mjj);                                                         
  tree->Branch("mjj_minDPhiZMET",                                        			  &b2L2J_mjj_minDPhiZMET);                                             
  tree->Branch("Rjj",                                                    			  &b2L2J_Rjj);                                                         
  tree->Branch("Rjj_minDPhiZMET",                                        			  &b2L2J_Rjj_minDPhiZMET);                                                                                     
  tree->Branch("vectorSumJetsPt",                                        			  &b2L2J_vectorSumJetsPt);                                             
  tree->Branch("vectorSumJetsEta",                                       			  &b2L2J_vectorSumJetsEta);                                            
  tree->Branch("vectorSumJetsPhi",                                       			  &b2L2J_vectorSumJetsPhi);                                            
  tree->Branch("vectorSumJetsM",                                         			  &b2L2J_vectorSumJetsM);                                              
  tree->Branch("met_Et",                                                 			  &b2L2J_met_Et);                                                      
  tree->Branch("met_Sign",                                               			  &b2L2J_met_Sign);                                                    
  tree->Branch("met_Phi",                                                			  &b2L2J_met_Phi);                                                     
  tree->Branch("TST_Et",                                                 			  &b2L2J_TST_Et);                                                      
  tree->Branch("TST_Phi",                                                			  &b2L2J_TST_Phi);                                                                                            
  tree->Branch("deltaPhi_MET_TST_Phi",                                   			  &b2L2J_deltaPhi_MET_TST_Phi);                                                                          
  tree->Branch("Ht30",                                                   			  &b2L2J_Ht30);
  tree->Branch("Rbb",                                                    			  &b2L2J_Rbb);                                                         
  tree->Branch("mbb",                                                    			  &b2L2J_mbb);                                                         
  tree->Branch("METOverPtZ",                                             			  &b2L2J_METOverPtZ);                                                  
  tree->Branch("METOverPtW",                                             			  &b2L2J_METOverPtW);                                                  
  tree->Branch("PtISR",                                                  			  &b2L2J_PtISR);                                                       
  tree->Branch("METOverPtISR",                                           			  &b2L2J_METOverPtISR);                                                
  tree->Branch("METOverHT",                                              			  &b2L2J_METOverHT);                                                   
  tree->Branch("METOverJ1pT",                                            			  &b2L2J_METOverJ1pT);                                                                          
  tree->Branch("DPhiJ1Met",                                              			  &b2L2J_DPhiJ1Met);                                                   
  tree->Branch("DPhiJ2Met",                                              			  &b2L2J_DPhiJ2Met);                                                   
  tree->Branch("DPhiJ3Met",                                              			  &b2L2J_DPhiJ3Met);                                                   
  tree->Branch("DPhiJ4Met",                                              			  &b2L2J_DPhiJ4Met);                                                   
  tree->Branch("minDPhi2JetsMet",                                        			  &b2L2J_minDPhi2JetsMet);                                             
  tree->Branch("minDPhi4JetsMet",                                        			  &b2L2J_minDPhi4JetsMet);                                             
  tree->Branch("minDPhiAllJetsMet",                                      			  &b2L2J_minDPhiAllJetsMet);                                                           
  tree->Branch("dPhiPjjMet",                                             			  &b2L2J_dPhiPjjMet);                                                  
  tree->Branch("dPhiPjjMet_minDPhiZMET",                                 			  &b2L2J_dPhiPjjMet_minDPhiZMET);                                      
  tree->Branch("dPhiMetISR",                                             			  &b2L2J_dPhiMetISR);                                                  
  tree->Branch("dPhiMetJet1",                                            			  &b2L2J_dPhiMetJet1);                                                 
  tree->Branch("METOverHTLep",                                           			  &b2L2J_METOverHTLep);                                                                     
  tree->Branch("mll",                                                    			  &b2L2J_mll);                                                         
  tree->Branch("mlll",                                                    			  &b2L2J_mlll);                                                         
  tree->Branch("mllll",                                                    			  &b2L2J_mllll);                                                         
  tree->Branch("Rll",                                                    			  &b2L2J_Rll);                                                         
  tree->Branch("Ptll",                                                   			  &b2L2J_Ptll);                                                        
  tree->Branch("absEtall",                                               			  &b2L2J_absEtall);                                                    
  tree->Branch("dPhiPllMet",                                             			  &b2L2J_dPhiPllMet);                                                  
  tree->Branch("dPhill",                                                 			  &b2L2J_dPhill);                                                                                               
  tree->Branch("mt2leplsp_0",                                            			  &b2L2J_mt2leplsp_0);                                                 
  tree->Branch("pileupWeight",                                          			  &b2L2J_pileupWeight);                                               
  tree->Branch("leptonWeight",                                          			  &b2L2J_leptonWeight);                                               
  tree->Branch("eventWeight",                                           			  &b2L2J_eventWeight);                                                
  tree->Branch("bTagWeight",                                            			  &b2L2J_bTagWeight);                                                 
  tree->Branch("jvtWeight",                                             			  &b2L2J_jvtWeight);                                                  
  tree->Branch("globalDiLepTrigSF",                                     			  &b2L2J_globalDiLepTrigSF);                                          
  tree->Branch("genWeight",                                             			  &b2L2J_genWeight);                                                  
  tree->Branch("genWeightUp",                                           			  &b2L2J_genWeightUp);                                                
  tree->Branch("genWeightDown",                                         			  &b2L2J_genWeightDown);                                              
  tree->Branch("truthMll",                                              			  &b2L2J_truthMll);                                                   
  tree->Branch("winoBinoMllWeight",                                     			  &b2L2J_winoBinoMllWeight);                                          
  tree->Branch("winoBinoXsecWeight",                                    			  &b2L2J_winoBinoXsecWeight);                                         
  tree->Branch("winoBinoBrFracWeight",                                  			  &b2L2J_winoBinoBrFracWeight);                                       
  tree->Branch("winoBinoWeight",                                        			  &b2L2J_winoBinoWeight);                                                        
  tree->Branch("ttbarNNLOWeight",                                       			  &b2L2J_ttbarNNLOWeight);                                            
  tree->Branch("ttbarNNLOWeightUp",                                     			  &b2L2J_ttbarNNLOWeightUp);                                          
  tree->Branch("ttbarNNLOWeightDown",                                   			  &b2L2J_ttbarNNLOWeightDown);                                        
  tree->Branch("truthTopPt",                                             			  &b2L2J_truthTopPt);                                                  
  tree->Branch("truthAntiTopPt",                                         			  &b2L2J_truthAntiTopPt);                                              
  tree->Branch("truthTtbarPt",                                           			  &b2L2J_truthTtbarPt);                                                
  tree->Branch("truthTtbarM",                                            			  &b2L2J_truthTtbarM);                                                 
  tree->Branch("x1",                                                     			  &b2L2J_x1);                                                          
  tree->Branch("x2",                                                     			  &b2L2J_x2);                                                          
  tree->Branch("pdf1",                                                   			  &b2L2J_pdf1);                                                        
  tree->Branch("pdf2",                                                   			  &b2L2J_pdf2);                                                        
  tree->Branch("scalePDF",                                               			  &b2L2J_scalePDF);                                                    
  tree->Branch("id1",                                                      			  &b2L2J_id1);                                                           
  tree->Branch("id2",                                                      			  &b2L2J_id2);
  tree->Branch("nLeps_RJ",                                &b2L2J_nLeps_RJ);		  
  tree->Branch("nJets_RJ",				  &b2L2J_nJets_RJ);		  
  tree->Branch("nBtagJets_RJ",				  &b2L2J_nBtagJets_RJ);		  
  tree->Branch("jetPt_RJ",				  &b2L2J_jetPt_RJ);		  
  tree->Branch("jetEta_RJ",				  &b2L2J_jetEta_RJ);		  
  tree->Branch("jetPhi_RJ",				  &b2L2J_jetPhi_RJ);		  
  tree->Branch("jetM_RJ",				  &b2L2J_jetM_RJ);		  
  tree->Branch("lepPt_RJ",  				  &b2L2J_lepPt_RJ);  		  
  tree->Branch("lepEta_RJ", 				  &b2L2J_lepEta_RJ); 		  
  tree->Branch("lepPhi_RJ", 				  &b2L2J_lepPhi_RJ); 		  
  tree->Branch("lepE_RJ",   				  &b2L2J_lepE_RJ);   		  
  tree->Branch("lepsign_RJ",				  &b2L2J_lepsign_RJ);		  
  tree->Branch("is2Lep2Jet",				  &b2L2J_is2Lep2Jet);		  
  tree->Branch("is2L2JInt",				  &b2L2J_is2L2JInt);		  
  tree->Branch("is3Lep",				  &b2L2J_is3Lep);			  
  tree->Branch("is3LInt",				  &b2L2J_is3LInt);		  
  tree->Branch("is3Lep2Jet",				  &b2L2J_is3Lep2Jet);		  
  tree->Branch("is3Lep3Jet",				  &b2L2J_is3Lep3Jet);		  
  tree->Branch("is4Lep2Jet",				  &b2L2J_is4Lep2Jet);		  
  tree->Branch("is4Lep3Jet",				  &b2L2J_is4Lep3Jet);		  
  tree->Branch("mll_RJ",				  &b2L2J_mll_RJ);			  
  tree->Branch("H2PP",					  &b2L2J_H2PP);			  
  tree->Branch("H5PP",					  &b2L2J_H5PP);			  
  tree->Branch("RPT_HT5PP",                               &b2L2J_RPT_HT5PP);                        
  tree->Branch("R_minH2P_minH3P",                         &b2L2J_R_minH2P_minH3P);                  
  tree->Branch("R_H2PP_H5PP",                             &b2L2J_R_H2PP_H5PP);                 
  tree->Branch("dphiVP",                               	  &b2L2J_dphiVP);                   
  tree->Branch("minDphi",                              	  &b2L2J_minDphi);                  
  tree->Branch("mTW",                             	  &b2L2J_mTW);                      
  tree->Branch("H4PP",                            	  &b2L2J_H4PP);                     
  tree->Branch("RPT_HT4PP",                           	  &b2L2J_RPT_HT4PP);                
  tree->Branch("R_HT4PP_H4PP",                            &b2L2J_R_HT4PP_H4PP);             
  tree->Branch("PTISR",                         	  &b2L2J_PTISR);                    
  tree->Branch("RISR",                        		  &b2L2J_RISR);                     
  tree->Branch("PTI",                       		  &b2L2J_PTI);                      
  tree->Branch("dphiISRI",                      	  &b2L2J_dphiISRI);                 
  tree->Branch("PTCM",                     		  &b2L2J_PTCM);                     
  tree->Branch("NjS",                    		  &b2L2J_NjS);                      
  tree->Branch("NjISR",                   		  &b2L2J_NjISR);                    
  tree->Branch("MZ",                  			  &b2L2J_MZ);                  	  
  tree->Branch("MJ",                 			  &b2L2J_MJ);                 	  
  tree->Branch("mTl3",                			  &b2L2J_mTl3);                	  
  tree->Branch("lept1Pt_VR",               		  &b2L2J_lept1Pt_VR);               
  tree->Branch("lept1sign_VR",              		  &b2L2J_lept1sign_VR);             
  tree->Branch("lept2Pt_VR",             		  &b2L2J_lept2Pt_VR);               
  tree->Branch("lept2sign_VR",            		  &b2L2J_lept2sign_VR);             
  tree->Branch("mll_RJ_VR",           			  &b2L2J_mll_RJ_VR);           	  
  tree->Branch("H2PP_VR",          			  &b2L2J_H2PP_VR);          	  
  tree->Branch("H5PP_VR",         			  &b2L2J_H5PP_VR);         	  
  tree->Branch("RPT_HT5PP_VR",        			  &b2L2J_RPT_HT5PP_VR);        	  
  tree->Branch("R_minH2P_minH3P_VR",       		  &b2L2J_R_minH2P_minH3P_VR);       
  tree->Branch("R_H2PP_H5PP_VR",      			  &b2L2J_R_H2PP_H5PP_VR);      	  
  tree->Branch("dphiVP_VR",     			  &b2L2J_dphiVP_VR);     		  
  tree->Branch("PTISR_VR",    				  &b2L2J_PTISR_VR);    		  
  tree->Branch("RISR_VR",   				  &b2L2J_RISR_VR);   		  
  tree->Branch("PTI_VR",  				  &b2L2J_PTI_VR);  		  
  tree->Branch("dphiISRI_VR", 				  &b2L2J_dphiISRI_VR); 		  
  tree->Branch("PTCM_VR",				  &b2L2J_PTCM_VR);		  
  tree->Branch("NjS_VR",				  &b2L2J_NjS_VR);			  
  tree->Branch("NjISR_VR",				  &b2L2J_NjISR_VR);		  
  tree->Branch("MZ_VR",					  &b2L2J_MZ_VR);			  
  tree->Branch("MJ_VR",					  &b2L2J_MJ_VR);			  
  tree->Branch("PRWHash",				  &b2L2J_PRWHash);		  
  tree->Branch("EventNumber",				  &b2L2J_EventNumber);		  
  tree->Branch("xsec",					  &b2L2J_xsec);			  
  tree->Branch("GenHt",					  &b2L2J_GenHt);			  
  tree->Branch("GenMET",				  &b2L2J_GenMET);			  
  tree->Branch("DatasetNumber",				  &b2L2J_DatasetNumber);		  
  tree->Branch("RunNumber",				  &b2L2J_RunNumber);		  
  tree->Branch("RandomRunNumber",			  &b2L2J_RandomRunNumber);	  
  tree->Branch("FS",                       		  &b2L2J_FS);                                                                                    



  // tree->Branch("fake_wgt",                                                                            &b2L2J_fake_wgt);    
  // tree->Branch("fake_wgt_up",                                                                            &b2L2J_fake_wgt_up);    
  // tree->Branch("fake_wgt_dw",                                                                            &b2L2J_fake_wgt_dw);    

  tree->Branch("isTT",                                                                            &b2L2J_isTT);    
  tree->Branch("isTl",                                                                            &b2L2J_isTl);    
  tree->Branch("islT",                                                                            &b2L2J_islT);    
  tree->Branch("isll",                                                                            &b2L2J_isll);    
  tree->Branch("r1",                                                                            &b2L2J_r1);    
  tree->Branch("r2",                                                                            &b2L2J_r2);    
  tree->Branch("f1",                                                                            &b2L2J_f1);    
  tree->Branch("f2",                                                                            &b2L2J_f2);    

  tree->Branch("nLep_base_OR", &b2L2J_nLep_base_OR);   
  
  ClearOutputBranches();
}


void make2L2JTree::ClearOutputBranches(void)
{

  b2L2J_lepFlavor.clear();                                  
  b2L2J_lepCharge.clear();                                  
  b2L2J_lepAuthor.clear();                                  
  b2L2J_lepPt.clear();                                      
  b2L2J_lepEta.clear();                                     
  b2L2J_lepPhi.clear();
b2L2J_lepIFFClass.clear();
  b2L2J_lepM.clear();                                       
  b2L2J_lepD0.clear();                                      
  b2L2J_lepD0Sig.clear();                                   
  b2L2J_lepZ0.clear();                                      
  b2L2J_lepZ0SinTheta.clear();                                   
  b2L2J_lepPassOR.clear();                                  
  b2L2J_lepType.clear();                                    
  b2L2J_lepOrigin.clear();                                  
  b2L2J_lepEgMotherType.clear();                            
  b2L2J_lepEgMotherOrigin.clear();                          
  b2L2J_lepEgMotherPdgId.clear();                           
  b2L2J_lepECIDS.clear();                                   
  b2L2J_lepNPix.clear();                                                       
  b2L2J_lepPassBL.clear();                                     
  b2L2J_lepSignal.clear();                                  
  b2L2J_lepTruthMatched.clear();                            
  b2L2J_lepTruthCharge.clear();                             
  b2L2J_lepTruthPt.clear();                                 
  b2L2J_lepTruthEta.clear();                                
  b2L2J_lepTruthPhi.clear();                                
  b2L2J_lepTruthM.clear();                                  
  b2L2J_jetPt.clear();        
  b2L2J_jetEta.clear();       
  b2L2J_jetPhi.clear();       
  b2L2J_jetM.clear();  	
  b2L2J_jetTileEnergy.clear();
  b2L2J_jetMV2c10.clear();

  b2L2J_jetPt_RJ.clear();		  
  b2L2J_jetEta_RJ.clear();		  
  b2L2J_jetPhi_RJ.clear();		  
  b2L2J_jetM_RJ.clear();		  
  b2L2J_lepPt_RJ.clear();  		  
  b2L2J_lepEta_RJ.clear(); 		  
  b2L2J_lepPhi_RJ.clear(); 		  
  b2L2J_lepE_RJ.clear();   		  
  b2L2J_lepsign_RJ.clear();	

  return;
}





//-----------------------------------------------------------------------------------------------------------
make2L2JTree::~make2L2JTree()
{
    // Write out the output tree and close the output file
  file->Write();
  file->Close();
  delete file;
}


void make2L2JTree::WriteTree()
{
  //file->cd();
    tree->Fill();
    //tree->Write();
    //file->Write();
    //file->Close();
    ClearOutputBranches();
}




void make2L2JTree::setSumOfMcWeights(double sumOfMcWeights)
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



