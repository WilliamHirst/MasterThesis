#define myTree_cxx
// The class definition in myTree.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("myTree.C")
// root> T->Process("myTree.C","some options")
// root> T->Process("myTree.C+")
//

#include "make2L0JTree.h"
#include "myTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <iostream>
#include <cmath>
#include <TLorentzVector.h>
#include <stdio.h>
#include <cstdlib>
#include <fstream>
std::string str;
#define pointTo(x,y) do{ str = #y; x=&y;   }while(0)

double sumofweights = 0;
double sumofweights_up = 0;
double sumofweights_dw = 0;
int nevpass = 0;
int nevpassnbjet1 = 0;
int nevpassem = 0;

double zmass = 91.1876;

bool doExtraVar = false;

#define PRINTER(name) printer(#name, (name))

void printer(char *name, int value) {
  printf("name: %s\tvalue: %d\n", name, value);
}



void myTree::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

 
  nentries = 0;
  n_entry = 0;
  looseTightDef.push_back("notrigm");
  looseTightDef.push_back("med");
  looseTightDef.push_back("loose");
  looseTightDef.push_back("tight_iso");
  looseTightDef.push_back("no");

  //uncert_keys.push_back("MUON_STAT_1up");
  //uncert_keys.push_back("MUON_STAT_1down");
  uncert_keys.push_back("MUON_EFF_TOTAL_1up");
  uncert_keys.push_back("MUON_EFF_TOTAL_1down");
  uncert_keys.push_back("MUON_EFF_TRIG_TOTAL_1up");
  uncert_keys.push_back("MUON_EFF_TRIG_TOTAL_1down");
  uncert_keys.push_back("EL_EFF_TRIG_TOTAL_1down");
  //uncert_keys.push_back("EL_STAT_1up");
  //uncert_keys.push_back("EL_STAT_1down");
  uncert_keys.push_back("EL_EFF_TOTAL_1up");
  uncert_keys.push_back("EL_EFF_TOTAL_1down");
  uncert_keys.push_back("EL_EFF_Trigger_TOTAL_1up");
  uncert_keys.push_back("EL_EFF_Trigger_TOTAL_1down");

  TString option = GetOption();

  yr = option.Atoi();


  if(0){
    ROOT::RDataFrame rdf("FNP_WEIGHTS","FNP_270521/FNP_EirikTest_020621.root",{"RunNumber","EventNumber","BDTDeltaM100_90","BDTVVDeltaM100_90","BDTtopDeltaM100_90","BDTothersDeltaM100_90",
									       "lepflav1","lepflav2","isSignalLep1","isSignalLep2","lep1pT","lep2pT",
									       "njet","nbjet","METsig","isSF","isOS","FNP_TOTAL_UP","FNP_TOTAL_DOWN","FNP_WEIGHTS"});
  
    auto filter = rdf.Filter("1");//is1516 ? ("RunNumber >=  276262 && RunNumber <= 311481") : (is17 ? ("RunNumber >=  325713 && RunNumber <= 340453") : ("RunNumber >= 348885")));
  
    auto evnum = filter.Take<Long64_t>("EventNumber");
    auto runnum = filter.Take<Int_t>("RunNumber");

    auto njet = filter.Take<Int_t>("njet");
    auto nbjet = filter.Take<Int_t>("nbjet");
    auto isSF = filter.Take<Bool_t>("isSF");
    auto isOS = filter.Take<Bool_t>("isOS");
    auto METsig = filter.Take<Float_t>("METsig");

    auto FNP_TOTAL_UP = filter.Take<Double_t>("FNP_TOTAL_UP");
    auto FNP_TOTAL_DOWN = filter.Take<Double_t>("FNP_TOTAL_DOWN");
    auto FNP_WEIGHTS = filter.Take<Double_t>("FNP_WEIGHTS");
    
    auto vecto_njet = njet.GetValue();
    auto vecto_nbjet = nbjet .GetValue();
    auto vecto_isSF = isSF.GetValue();
    auto vecto_isOS = isOS.GetValue();
    auto vecto_METsig = METsig.GetValue();

    auto vecto_FNP_TOTAL_UP = FNP_TOTAL_UP.GetValue();
    auto vecto_FNP_TOTAL_DOWN = FNP_TOTAL_DOWN.GetValue();
    auto vecto_FNP_WEIGHTS = FNP_WEIGHTS.GetValue();

    auto lepflav1 = filter.Take<Int_t>("lepflav1");
    auto lepflav2 = filter.Take<Int_t>("lepflav2");

    auto lepsignal1 = filter.Take<Bool_t>("isSignalLep1");
    auto lepsignal2 = filter.Take<Bool_t>("isSignalLep2");

    auto lep1pT = filter.Take<Float_t>("lep1pT");
    auto lep2pT = filter.Take<Float_t>("lep2pT");

    auto vecto_lepflav1 = lepflav1.GetValue();
    auto vecto_lepflav2 = lepflav2.GetValue();

    auto vecto_lepsignal1 = lepsignal1.GetValue();
    auto vecto_lepsignal2 = lepsignal2.GetValue();

    auto vecto_lep1pT = lep1pT.GetValue();
    auto vecto_lep2pT = lep2pT.GetValue();

    auto vecto_evnum = evnum.GetValue();
    auto vecto_runnum = runnum.GetValue();

    bdt_vec = {"BDTDeltaM100_90","BDTVVDeltaM100_90","BDTtopDeltaM100_90","BDTothersDeltaM100_90"};//"BDTDeltaM100_90","BDTVVDeltaM100_90","BDTtopDeltaM100_90","BDTothersDeltaM100_90","BDTDeltaM30","BDTVVDeltaM30","BDTtopDeltaM30","BDTothersDeltaM30"};
    for (std::string bdt_score : bdt_vec) {
      auto vect = filter.Take<double>(bdt_score);
      TString key = bdt_score;
      BDTscoremap[key] = vect.GetValue();
      h_real1_diff[key] = new TH1F(Form("h_real1_diff_%s",key.Data()),Form("h_real1_diff_%s",key.Data()),400,-1,1); fOutput->AddLast(h_real1_diff[key]);
      h_real2_diff[key] = new TH1F(Form("h_real2_diff_%s",key.Data()),Form("h_real2_diff_%s",key.Data()),400,-1,1); fOutput->AddLast(h_real2_diff[key]);
      h_lep_BDT_n[key] = new TH1F(Form("h_lep_BDT_n_%s",key.Data()),Form("h_lep_BDT_n_%s",key.Data()),20,0,1); fOutput->AddLast(h_lep_BDT_n[key]);
      h_lep_BDT_r[key] = new TH2F(Form("h_lep_BDT_r_%s",key.Data()),Form("h_lep_BDT_r_%s",key.Data()),20,0,1,100,0,1); fOutput->AddLast(h_lep_BDT_r[key]);
    }

    std::cout<<"number of events = "<<vecto_evnum.size()<<endl;
    std::cout<<"number of events = "<<vecto_runnum.size()<<endl;
    int i = 0;
    for(const auto& value: vecto_runnum) {
      evnums[value].push_back(vecto_evnum.at(i));
      //std::cout << "Value = "<< value << "\n";
      //if(i > 100)break;
      evtype[value].push_back(Form("%i%i",vecto_lepflav1.at(i),vecto_lepflav2.at(i)));
      lep1isSignal[value].push_back(vecto_lepsignal1.at(i));
      lep2isSignal[value].push_back(vecto_lepsignal2.at(i));
      lepton1pT[value].push_back(vecto_lep1pT.at(i));
      lepton2pT[value].push_back(vecto_lep2pT.at(i));

      glob_njet[value].push_back(vecto_njet.at(i));
      glob_nbjet[value].push_back(vecto_nbjet.at(i));
      glob_isSF[value].push_back(vecto_isSF.at(i));
      glob_isOS[value].push_back(vecto_isOS.at(i));
      glob_METsig[value].push_back(vecto_METsig.at(i));

      glob_FNP_TOTAL_UP[value].push_back(vecto_FNP_TOTAL_UP.at(i));
      glob_FNP_TOTAL_DOWN[value].push_back(vecto_FNP_TOTAL_DOWN.at(i));
      glob_FNP_WEIGHTS[value].push_back(vecto_FNP_WEIGHTS.at(i));
      
      // bdt_vec = {"BDTDeltaM100_90","BDTVVDeltaM100_90","BDTtopDeltaM100_90","BDTothersDeltaM100_90","BDTDeltaM30","BDTVVDeltaM30","BDTtopDeltaM30","BDTothersDeltaM30"};
      for (std::string bdt_score : bdt_vec) {
	TString key = bdt_score;
	BDTscores[value][key].push_back(BDTscoremap[key].at(i));
	
      }
      
      i = i+1;
    }

    
  }
  h_r1_f1_TT = new TH2F("h_r1_f1_TT","h_r1_f1_TT",100,0,1,100,0,1); fOutput->AddLast(h_r1_f1_TT);
  h_r1_f1_Tl = new TH2F("h_r1_f1_Tl","h_r1_f1_Tl",100,0,1,100,0,1); fOutput->AddLast(h_r1_f1_Tl);
  h_r1_f1_lT = new TH2F("h_r1_f1_lT","h_r1_f1_lT",100,0,1,100,0,1); fOutput->AddLast(h_r1_f1_lT);
  h_r1_f1_ll = new TH2F("h_r1_f1_ll","h_r1_f1_ll",100,0,1,100,0,1); fOutput->AddLast(h_r1_f1_ll);
  h_r2_f2_TT = new TH2F("h_r2_f2_TT","h_r2_f2_TT",100,0,1,100,0,1); fOutput->AddLast(h_r2_f2_TT);
  h_r2_f2_Tl = new TH2F("h_r2_f2_Tl","h_r2_f2_Tl",100,0,1,100,0,1); fOutput->AddLast(h_r2_f2_Tl);
  h_r2_f2_lT = new TH2F("h_r2_f2_lT","h_r2_f2_lT",100,0,1,100,0,1); fOutput->AddLast(h_r2_f2_lT);
  h_r2_f2_ll = new TH2F("h_r2_f2_ll","h_r2_f2_ll",100,0,1,100,0,1); fOutput->AddLast(h_r2_f2_ll);

  h_tmatch_smalldiff = new TH1F("h_tmatch_smalldiff","h_tmatch_smalldiff",looseTightDef.size(),0,looseTightDef.size()); fOutput->AddLast(h_tmatch_smalldiff);
  int ibin = 1;
  for(std::vector< TString >::iterator itLT = looseTightDef.begin(); itLT != looseTightDef.end(); itLT++){
    h_tmatch_smalldiff->GetXaxis()->SetBinLabel(ibin,(*itLT).Data());
    ibin += 1;
  }
  h_eta_smalldiff = new TH1F("h_eta_smalldiff","h_eta_smalldiff",320,-2.6,2.6); fOutput->AddLast(h_eta_smalldiff);
  
  // if(0){
  //   ROOT::RDataFrame rdf("FNP_WEIGHTS","FNP_270521/FNP_EirikTest_020621.root",{"RunNumber","EventNumber","BDTDeltaM100_90","BDTVVDeltaM100_90","BDTtopDeltaM100_90","BDTothersDeltaM100_90","lepflav1","lepflav2"});
  
  //   auto filter = rdf.Filter((yr == 2015 || yr == 2016) ? ("RunNumber >=  276262 && RunNumber <= 311481") : ((yr == 2017) ? ("RunNumber >=  325713 && RunNumber <= 340453") : ("RunNumber >= 348885")));
  
  //   auto evnum = filter.Take<Long64_t>("EventNumber");
  //   auto runnum = filter.Take<Int_t>("RunNumber");
  //   /**
  //      auto lepflav1 = filter.Take<Int_t>("lepflav1");
  //      auto lepflav2 = filter.Take<Int_t>("lepflav2");

  //      auto vecto_lepflav1 = lepflav1.GetValue();
  //      auto vecto_lepflav2 = lepflav2.GetValue();
  //   */
  //   auto vecto_evnum = evnum.GetValue();
  //   auto vecto_runnum = runnum.GetValue();

  //   std::cout<<"number of events = "<<vecto_evnum.size()<<endl;
  //   std::cout<<"number of events = "<<vecto_runnum.size()<<endl;
  //   int i = 0;
  //   for(const auto& value: vecto_runnum) {
  //     evnums[value].push_back(vecto_evnum.at(i));
  //     //std::cout << "Value = "<< value << "\n";
  //     //if(i > 100)break;
  //     //evtype.push_back(Form("%i%i",vecto_lepflav1.at(i),vecto_lepflav2.at(i)));
  //     i = i+1;
  //   }

  //   bdt_vec = {"BDTDeltaM100_90","BDTVVDeltaM100_90","BDTtopDeltaM100_90","BDTothersDeltaM100_90","BDTDeltaM30","BDTVVDeltaM30","BDTtopDeltaM30","BDTothersDeltaM30"};
  //   for (std::string bdt_score : bdt_vec) {
  //     auto vect = filter.Take<double>(bdt_score);
  //     TString key = bdt_score;
  //     BDTscores[key] = vect.GetValue();
  //     h_lep_BDT_n[key] = new TH1F(Form("h_lep_BDT_n_%s",key.Data()),Form("h_lep_BDT_n_%s",key.Data()),20,0,1);
  //     h_lep_BDT_r[key] = new TH2F(Form("h_lep_BDT_r_%s",key.Data()),Form("h_lep_BDT_r_%s",key.Data()),20,0,1,100,0,1);
  //   }
  // }

  filename_prev = "";
  filename = "filename";
  gfw = new GetFakeWeight();
  gfw->setVB(0);
  std::vector<TString> unc_keys;
  for (auto u : GetFakeWeight::all_unc) {
    unc_keys.push_back(gfw->getUncKey(u));
  }
  twoLzeroJ = new make2L0JTree(option,"central",Form("FNP_BKG_%s.root",option.Data()),"myTree",unc_keys, doExtraVar);
  for(std::vector< TString >::iterator itLT = looseTightDef.begin(); itLT != looseTightDef.end(); itLT++){ 
    gfw->initializeHist((*itLT), "./MMinputfiles/FNPntuple_FCL_FEB10_final_v2.root","./MMinputfiles/MMinput_frac_NEW_FRAC_JUNE18.root" , 0, "metsig", true, true, true, uncert_keys);   
  }

  //"./MMinputfiles/MMinput_frac_NEW_FRAC_JUNE18.root"

  // final_trigger_setup  = {"e_medloose":{"triggers":["e24_lhmedium_L1EM20VH","e24_lhmedium_nod0_L1EM20VH","e60_lhmedium_nod0","e120_lhloose","e140_lhloose_nod0"]},
  //                         "e_tight_iso":{"triggers":["e24_lhtight_nod0_ivarloose","e26_lhtight_nod0_ivarloose"]},
  //                         "e_no":{"triggers":["e300_etcut"]},
  //                         "m_medloose":{"triggers":["mu26_ivarmedium"]},
  //                         "m_no":{"triggers":["mu60_0eta105_msonly","mu50","mu20_iloose_L1MU15"]}
  //                       }

  // 2015
  triglist_e[2015].push_back(HLT_e24_lhmedium_L1EM20VH_Nominal);
  triglist_e[2015].push_back(HLT_e60_lhmedium_Nominal);
  triglist_e[2015].push_back(HLT_e120_lhloose_Nominal);
  
  triglist_m[2015].push_back(HLT_mu20_iloose_L1MU15_Nominal);
  //triglist_m[2015].push_back(HLT_mu60_0eta105_msonly_Nominal);
  triglist_m[2015].push_back(HLT_mu50_Nominal);

  // 2016
  triglist_e[2016].push_back(HLT_e24_lhmedium_nod0_L1EM20VH_Nominal);
  triglist_e[2016].push_back(HLT_e24_lhtight_nod0_ivarloose_Nominal);
  triglist_e[2016].push_back(HLT_e26_lhtight_nod0_ivarloose_Nominal);
  triglist_e[2016].push_back(HLT_e60_lhmedium_nod0_Nominal);
  triglist_e[2016].push_back(HLT_e140_lhloose_nod0_Nominal);
  // triglist_e[2016].push_back(HLT_e300_etcut_Nominal);
    
  // triglist_m[2016].push_back(HLT_mu60_0eta105_msonly_Nominal);
  triglist_m[2016].push_back(HLT_mu26_ivarmedium_Nominal);
  triglist_m[2016].push_back(HLT_mu50_Nominal);

  // 2017
  triglist_e[2017].push_back(HLT_e26_lhtight_nod0_ivarloose_Nominal);
  triglist_e[2017].push_back(HLT_e60_lhmedium_nod0_Nominal);
  triglist_e[2017].push_back(HLT_e140_lhloose_nod0_Nominal);
  //triglist_e[2017].push_back(HLT_e300_etcut_Nominal);
  
  triglist_m[2017].push_back(HLT_mu26_ivarmedium_Nominal);
  //triglist_m[2017].push_back(HLT_mu60_0eta105_msonly_Nominal);
  triglist_m[2017].push_back(HLT_mu50_Nominal);

  // 2018
  triglist_e[2018].push_back(HLT_e26_lhtight_nod0_ivarloose_Nominal);
  triglist_e[2018].push_back(HLT_e60_lhmedium_nod0_Nominal);
  triglist_e[2018].push_back(HLT_e140_lhloose_nod0_Nominal);
  //triglist_e[2018].push_back(HLT_e300_etcut_Nominal);
  
  triglist_m[2018].push_back(HLT_mu26_ivarmedium_Nominal);
  //triglist_m[2018].push_back(HLT_mu60_0eta105_msonly_Nominal);
  triglist_m[2018].push_back(HLT_mu50_Nominal);

  // TriggerMatch
  // 2015
  trigmatch_e[2015].push_back(&HLT_e24_lhmedium_L1EM20VH_match_Nominal);
  trigmatch_e[2015].push_back(&HLT_e60_lhmedium_match_Nominal);
  trigmatch_e[2015].push_back(&HLT_e120_lhloose_match_Nominal);
  
  trigmatch_m[2015].push_back(&HLT_mu20_iloose_L1MU15_match_Nominal);
  //trigmatch_m[2015].push_back(&HLT_mu60_0eta105_msonly_match_Nominal);
  trigmatch_m[2015].push_back(&HLT_mu50_match_Nominal);

  // 2016
  trigmatch_e[2016].push_back(&HLT_e24_lhmedium_nod0_L1EM20VH_match_Nominal);
  trigmatch_e[2016].push_back(&HLT_e24_lhtight_nod0_ivarloose_match_Nominal);
  trigmatch_e[2016].push_back(&HLT_e26_lhtight_nod0_ivarloose_match_Nominal);
  trigmatch_e[2016].push_back(&HLT_e60_lhmedium_nod0_match_Nominal);
  trigmatch_e[2016].push_back(&HLT_e140_lhloose_nod0_match_Nominal);
  //trigmatch_e[2016].push_back(&HLT_e300_etcut_match_Nominal);
    
  // trigmatch_m[2016].push_back(&HLT_mu60_0eta105_msonly_match_Nominal);
  trigmatch_m[2016].push_back(&HLT_mu26_ivarmedium_match_Nominal);
  trigmatch_m[2016].push_back(&HLT_mu50_match_Nominal);

  // 2017
  trigmatch_e[2017].push_back(&HLT_e26_lhtight_nod0_ivarloose_match_Nominal);
  trigmatch_e[2017].push_back(&HLT_e60_lhmedium_nod0_match_Nominal);
  trigmatch_e[2017].push_back(&HLT_e140_lhloose_nod0_match_Nominal);
  //trigmatch_e[2017].push_back(&HLT_e300_etcut_match_Nominal);
  
  trigmatch_m[2017].push_back(&HLT_mu26_ivarmedium_match_Nominal);
  //trigmatch_m[2017].push_back(&HLT_mu60_0eta105_msonly_match_Nominal);
  trigmatch_m[2017].push_back(&HLT_mu50_match_Nominal);

  // 2018
  trigmatch_e[2018].push_back(&HLT_e26_lhtight_nod0_ivarloose_match_Nominal);
  trigmatch_e[2018].push_back(&HLT_e60_lhmedium_nod0_match_Nominal);
  trigmatch_e[2018].push_back(&HLT_e140_lhloose_nod0_match_Nominal);
  //trigmatch_e[2018].push_back(&HLT_e300_etcut_match_Nominal);
  
  trigmatch_m[2018].push_back(&HLT_mu26_ivarmedium_match_Nominal);
  //  trigmatch_m[2018].push_back(&HLT_mu60_0eta105_msonly_match_Nominal);
  trigmatch_m[2018].push_back(&HLT_mu50_match_Nominal);

  // TriggerName (IMPORTANT: must have same order as trigmatch (no so important with triglist)!!!)
  // 2015
  trigstr_e[2015].push_back("HLT_e24_lhmedium_L1EM20VH");
  trigstr_e[2015].push_back("HLT_e60_lhmedium");
  trigstr_e[2015].push_back("HLT_e120_lhloose");
  
  trigstr_m[2015].push_back("HLT_mu20_iloose_L1MU15");
  //  trigstr_m[2015].push_back("HLT_mu60_0eta105_msonly");
  trigstr_m[2015].push_back("HLT_mu50");

  // 2016
  trigstr_e[2016].push_back("HLT_e24_lhmedium_nod0_L1EM20VH");
  trigstr_e[2016].push_back("HLT_e24_lhtight_nod0_ivarloose");
  trigstr_e[2016].push_back("HLT_e26_lhtight_nod0_ivarloose");
  trigstr_e[2016].push_back("HLT_e60_lhmedium_nod0");
  trigstr_e[2016].push_back("HLT_e140_lhloose_nod0");
  // trigstr_e[2016].push_back("HLT_e300_etcut");
    
  //  trigstr_m[2016].push_back("HLT_mu60_0eta105_msonly");
  trigstr_m[2016].push_back("HLT_mu26_ivarmedium");
  trigstr_m[2016].push_back("HLT_mu50");

  // 2017
  trigstr_e[2017].push_back("HLT_e26_lhtight_nod0_ivarloose");
  trigstr_e[2017].push_back("HLT_e60_lhmedium_nod0");
  trigstr_e[2017].push_back("HLT_e140_lhloose_nod0");
  //trigstr_e[2017].push_back("HLT_e300_etcut");
  
  trigstr_m[2017].push_back("HLT_mu26_ivarmedium");
  //  trigstr_m[2017].push_back("HLT_mu60_0eta105_msonly");
  trigstr_m[2017].push_back("HLT_mu50");

  // 2018
  trigstr_e[2018].push_back("HLT_e26_lhtight_nod0_ivarloose");
  trigstr_e[2018].push_back("HLT_e60_lhmedium_nod0");
  trigstr_e[2018].push_back("HLT_e140_lhloose_nod0");
  // trigstr_e[2018].push_back("HLT_e300_etcut");
  
  trigstr_m[2018].push_back("HLT_mu26_ivarmedium");
  //  trigstr_m[2018].push_back("HLT_mu60_0eta105_msonly");
  trigstr_m[2018].push_back("HLT_mu50");

  // TriggerCategory
  // 2015
  trigcat["HLT_e24_lhmedium_L1EM20VH"] = "med";
  trigcat["HLT_e60_lhmedium"] = "med";
  trigcat["HLT_e120_lhloose"] = "loose";
  
  trigcat["HLT_mu20_iloose_L1MU15"] = "loose";
  //  trigcat["HLT_mu60_0eta105_msonly"] = "no";
  trigcat["HLT_mu50"] = "no";

  // 2016
  trigcat["HLT_e24_lhmedium_nod0_L1EM20VH"] = "med";
  trigcat["HLT_e24_lhtight_nod0_ivarloose"] = "tight_iso";
  trigcat["HLT_e26_lhtight_nod0_ivarloose"] = "tight_iso";
  trigcat["HLT_e60_lhmedium_nod0"] = "med";
  trigcat["HLT_e140_lhloose_nod0"] = "loose";
  //trigcat["HLT_e300_etcut"] = "no";
    
  //  trigcat["HLT_mu60_0eta105_msonly"] = "no";
  trigcat["HLT_mu26_ivarmedium"] = "med";
  trigcat["HLT_mu50"] = "no";

  // 2017
  trigcat["HLT_e26_lhtight_nod0_ivarloose"] = "tight_iso";
  trigcat["HLT_e60_lhmedium_nod0"] = "med";
  trigcat["HLT_e140_lhloose_nod0"] = "loose";
  //trigcat["HLT_e300_etcut"] = "no";
  
  trigcat["HLT_mu26_ivarmedium"] = "med";
  //  trigcat["HLT_mu60_0eta105_msonly"] = "no";
  trigcat["HLT_mu50"] = "no";

  // 2018
  trigcat["HLT_e26_lhtight_nod0_ivarloose"] = "tight_iso";
  trigcat["HLT_e60_lhmedium_nod0"] = "med";
  trigcat["HLT_e140_lhloose_nod0"] = "loose";
  //trigcat["HLT_e300_etcut"] = "no";
  
  trigcat["HLT_mu26_ivarmedium"] = "med";
  //  trigcat["HLT_mu60_0eta105_msonly"] = "no";
  trigcat["HLT_mu50"] = "no";

  trig_order.push_back("tight_iso");
  trig_order.push_back("med");
  trig_order.push_back("loose");
  trig_order.push_back("no");

  h_triggermatched_el1_el2 = new TH2F("h_triggermatched_el1_el2","",2,0,2,2,0,2);fOutput->AddLast(h_triggermatched_el1_el2);
  h_triggermatched_el1_el2->GetXaxis()->SetBinLabel(1,"NO");
  h_triggermatched_el1_el2->GetXaxis()->SetBinLabel(2,"YES");
  h_triggermatched_el1_el2->GetYaxis()->SetBinLabel(1,"NO");
  h_triggermatched_el1_el2->GetYaxis()->SetBinLabel(2,"YES");

  h_triggermatched_mu1_mu2 = new TH2F("h_triggermatched_mu1_mu2","",2,0,2,2,0,2);fOutput->AddLast(h_triggermatched_mu1_mu2);
  h_triggermatched_mu1_mu2->GetXaxis()->SetBinLabel(1,"NO");
  h_triggermatched_mu1_mu2->GetXaxis()->SetBinLabel(2,"YES");
  h_triggermatched_mu1_mu2->GetYaxis()->SetBinLabel(1,"NO");
  h_triggermatched_mu1_mu2->GetYaxis()->SetBinLabel(2,"YES");

  h_triggermatched_el_mu = new TH2F("h_triggermatched_el_mu","",2,0,2,2,0,2); fOutput->AddLast(h_triggermatched_el_mu);
  h_triggermatched_el_mu->GetXaxis()->SetBinLabel(1,"NO");
  h_triggermatched_el_mu->GetXaxis()->SetBinLabel(2,"YES");
  h_triggermatched_el_mu->GetYaxis()->SetBinLabel(1,"NO");
  h_triggermatched_el_mu->GetYaxis()->SetBinLabel(2,"YES");

  h_triggermatched_el = new TH1F("h_triggermatched_el","",trigstr_e[yr].size()+1,0,trigstr_e[yr].size()+1);fOutput->AddLast(h_triggermatched_el);
  for(unsigned int i = 0; i < trigstr_e[yr].size(); i++){h_triggermatched_el->GetXaxis()->SetBinLabel(i+2,trigstr_e[yr].at(i));}
  h_triggermatched_mu = new TH1F("h_triggermatched_mu","",trigstr_m[yr].size()+1,0,trigstr_m[yr].size()+1);fOutput->AddLast(h_triggermatched_mu);
  for(unsigned int i = 0; i < trigstr_m[yr].size(); i++){h_triggermatched_mu->GetXaxis()->SetBinLabel(i+2,trigstr_m[yr].at(i));}
 
  h_fake_weight = new TH1F("h_fake_weight","",4000,-20,20);fOutput->AddLast(h_fake_weight);
 
  h_cutflow = new TH1F("h_cutflow","",20,0,10);fOutput->AddLast(h_cutflow);
  h_cutflow->GetXaxis()->SetBinLabel(1,"No cuts");
  h_cutflow->GetXaxis()->SetBinLabel(2,"cleaningVeto");
  h_cutflow->GetXaxis()->SetBinLabel(3,"baseline>=2");
  h_cutflow->GetXaxis()->SetBinLabel(4,"Event trig");
  h_cutflow->GetXaxis()->SetBinLabel(5,"Presel");
  h_cutflow->GetXaxis()->SetBinLabel(6,"No trig match");
  h_cutflow->GetXaxis()->SetBinLabel(7,"2nd TM, pt < 27");
  h_cutflow->GetXaxis()->SetBinLabel(8,"njets == 1");
  h_cutflow->GetXaxis()->SetBinLabel(9,"nbets == 0");
  h_cutflow->GetXaxis()->SetBinLabel(10,"SF");
  h_cutflow->GetXaxis()->SetBinLabel(11,"OS");
  h_cutflow->GetXaxis()->SetBinLabel(12,"pT1 > 100");
  h_cutflow->GetXaxis()->SetBinLabel(13,"pT2 > 50");
  h_cutflow->GetXaxis()->SetBinLabel(14,"etsig > 7");
  h_cutflow->GetXaxis()->SetBinLabel(15,"mll > 60");
  h_cutflow->GetXaxis()->SetBinLabel(16,"zveto15");
  h_cutflow->GetXaxis()->SetBinLabel(17,"costhetastar<0.1");
  h_cutflow->GetXaxis()->SetBinLabel(18,"delphi>2.8");
  h_cutflow->GetXaxis()->SetBinLabel(19,"fakes");

  h_lep_fake_trig_cat = new TH2F("h_lep_fake_trig_cat","",6,0,6,6,0,6);fOutput->AddLast(h_lep_fake_trig_cat);
  int ii = 1;
  for (TString key : looseTightDef){
    h_lep_fake_trig_cat->GetXaxis()->SetBinLabel(ii,key);
    h_lep_fake_trig_cat->GetYaxis()->SetBinLabel(ii,key);
    ii += 1;
  }
  
 
  h_lep_tm = new TH2F("h_lep_tm","",2,0,2,10,0,10);fOutput->AddLast(h_lep_tm);
  for(int i = 1; i< 10; i++)h_lep_tm->GetYaxis()->SetBinLabel(i,Form("%i",i));
  h_lep_tm->GetXaxis()->SetBinLabel(1,"NOT TM");
  h_lep_tm->GetXaxis()->SetBinLabel(2,"TM");

  h_lep_count = new TH2F("h_lep_count","",4,0,4,10,0,10);fOutput->AddLast(h_lep_count);

  for(int i = 1; i< 10; i++)h_lep_count->GetYaxis()->SetBinLabel(i,Form("%i",i));
 
  h_lep_count->GetXaxis()->SetBinLabel(1,"BL");
  h_lep_count->GetXaxis()->SetBinLabel(2,"BL + TM");
  h_lep_count->GetXaxis()->SetBinLabel(3,"SG");
  h_lep_count->GetXaxis()->SetBinLabel(4,"SG + TM");

  h_costhetastar = new TH1F("h_costhetastar","",200,0,2);fOutput->AddLast(h_costhetastar);

  h_preselection_met_OS     = new TH1F("h_preselection_met_OS",""    ,30,0,30); fOutput->AddLast(h_preselection_met_OS);
  h_preselection_met_OS_FNP = new TH1F("h_preselection_met_OS_FNP","",30,0,30); fOutput->AddLast(h_preselection_met_OS_FNP);
  h_preselection_lep2pt_EE_OS     = new TH1F("h_preselection_lep2pt_EE_OS",""    ,200,0,200); fOutput->AddLast(h_preselection_lep2pt_EE_OS);
  h_preselection_lep2pt_EE_OS_FNP = new TH1F("h_preselection_lep2pt_EE_OS_FNP","",200,0,200); fOutput->AddLast(h_preselection_lep2pt_EE_OS_FNP);
  h_preselection_lep2pt_EM_OS     = new TH1F("h_preselection_lep2pt_EM_OS",""    ,200,0,200); fOutput->AddLast(h_preselection_lep2pt_EM_OS);
  h_preselection_lep2pt_EM_OS_FNP = new TH1F("h_preselection_lep2pt_EM_OS_FNP","",200,0,200); fOutput->AddLast(h_preselection_lep2pt_EM_OS_FNP);
  h_preselection_met_SS     = new TH1F("h_preselection_met_SS",""    ,30,0,30); fOutput->AddLast(h_preselection_met_SS);
  h_preselection_met_SS_FNP = new TH1F("h_preselection_met_SS_FNP","",30,0,30); fOutput->AddLast(h_preselection_met_SS_FNP);
 
}

void myTree::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t myTree::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // When processing keyed objects with PROOF, the object is already loaded
  // and is available via the fObject pointer.
  //
  // This function should contain the \"body\" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.
  int yr_num = -1;
  
  int idx1 = -1;
  int idx2 = -1;
  
  int nTT = 0;
  int nTl = 0;
  int nlT = 0;
  int nll = 0;

  int n_bjets = 0;
  int n_baseline_jets = 0;
  int n_signal_jets = 0;

  double r1 = 0.0;
  double f1 = 0.0;
  double r2 = 0.0;
  double f2 = 0.0;

  bool isEE = false;
  bool isEM = false;
  bool isMM = false;

  bool isOS = false;
  bool isSS = false;

  TString eventStr = "";

  std::map<TString,double> r_wgt  ;
  std::map<TString,double> rf_wgt ;
  std::map<TString,double> fr_wgt ;
  std::map<TString,double> f_wgt  ;
   
  int n_baseline = 0;
  int n_signal   = 0;
  std::vector<int> lep_idx;
  std::vector<int> lep_sg_idx;
  std::vector<int> lep_sg_tm_idx;
  std::vector<int> lep_tm_idx;
  
  TString yr_str;
  std::vector<TString> lep1_tmatch; 
  std::vector<TString> lep2_tmatch;

  TString triggercat1;
  TString triggercat2;
   

  isDF_1bjet = false;
  isCR_Top = false;
  fReader.SetLocalEntry(entry);

  h_cutflow->AddBinContent(1,1);
  fnpidx = -1;
  if(evnums.size() > 0){
    if ( evnums.find(*RunNumber) == evnums.end() ) {
      return kTRUE;
    }else{
      auto it = find(evnums[*RunNumber].begin(), evnums[*RunNumber].end(), *EventNumber);
      if(it == evnums[*RunNumber].end())return kTRUE;
      else fnpidx = it - evnums[*RunNumber].begin();
    }

    for(std::map<TString,std::vector<double>>::iterator iter = BDTscores[*RunNumber].begin(); iter != BDTscores[*RunNumber].end(); ++iter){
      //printf("Getting %s at index %i = %.2f\n",iter->first.Data(),fnpidx,(iter->second).at(fnpidx));
      BDTweight[iter->first] = (iter->second).at(fnpidx);
    }
    
    if(glob_njet[*RunNumber].at(fnpidx) == 1 && glob_nbjet[*RunNumber].at(fnpidx) == 1 && glob_METsig[*RunNumber].at(fnpidx) > 8 && !glob_isSF[*RunNumber].at(fnpidx) && !glob_isOS[*RunNumber].at(fnpidx)){
      isDF_1bjet = true;
    }
    if(glob_njet[*RunNumber].at(fnpidx) == 1 && glob_nbjet[*RunNumber].at(fnpidx) == 1 && glob_METsig[*RunNumber].at(fnpidx) > 8 &&
       ((!glob_isSF[*RunNumber].at(fnpidx) && BDTweight["BDTDeltaM100_90"] > 0.5  && BDTweight["BDTDeltaM100_90"] <= 0.7) ||
	(glob_isSF[*RunNumber].at(fnpidx) && BDTweight["BDTDeltaM100_90"] > 0.7  && BDTweight["BDTDeltaM100_90"] <= 0.75)) && glob_isOS[*RunNumber].at(fnpidx)){
      isCR_Top = true;
    }
    if(!isCR_Top)return kTRUE;
  }
    
  if(!nentries){
    nentries = ((TChain*)fReader.GetTree())->GetEntries();
  }
  n_entry += 1;

  if(n_entry%100000 == 0)std::cout<<"Entry "<<n_entry<<"/"<<nentries<<std::endl;

  // if(n_entry > 100000){
  //   Terminate();
  //   Abort("Enough!");
  // }
  if(!(*RunNumber == 276262 && (*EventNumber == 9853975 or *EventNumber == 15421209 or *EventNumber == 16207721 or *EventNumber == 17546845 or *EventNumber == 22437345 )))return kTRUE;
  
  //if((*RunNumber) != 307732)return kTRUE;
  //if((*EventNumber) != 1325337478)return kTRUE;
 
 if(*cleaningVeto_Nominal)return kTRUE;
 
 h_cutflow->AddBinContent(2,1);
 
 if((*RunNumber) >= 325713 && (*RunNumber) < 348885){yr_str = "2017";yr_num=2017;}
 else if ((*RunNumber) <= 284484){yr_str = "2015_2016";yr_num=2015;} // 2015
 else if ((*RunNumber) >= 297730 && (*RunNumber) <= 311481){yr_str = "2015_2016";yr_num=2016;} // 2016
 else if ((*RunNumber) >= 348885){yr_str = "2018";yr_num=2018;}
 else cout<<"ERROR Year for run "<<(*RunNumber)<<" is unknown"<<endl;

 if(yr != yr_num)return kTRUE;

   
  for(unsigned int i = 0; i < ptleptons_Nominal.GetSize(); i++){
    if(!lepL(i))continue;
    n_baseline += 1;
    lep_idx.push_back(i);
    //cout<<"Lepton "<<i<<" = "<<ptleptons_Nominal[i]<<endl;
    if(!lepT(i))continue;
    n_signal += 1;
    lep_sg_idx.push_back(i);
  }

  std::vector<int> lep_idx_sorted;
  std::vector<int> temp(lep_idx.size());
  std::copy(lep_idx.begin(), lep_idx.end(), temp.begin());
  //cout<<"temp size is "<< temp.size()<<endl;
  while(temp.size()>0){
    double maxpt = -999;
    int idx = -1;
    for(unsigned int i = 0; i<lep_idx.size(); i++){
      //cout<<"i = "<<i<<endl;
      if(std::count(lep_idx_sorted.begin(), lep_idx_sorted.end(), lep_idx.at(i)))continue;
      if(ptleptons_Nominal[lep_idx.at(i)] > maxpt){
	maxpt = ptleptons_Nominal[lep_idx.at(i)];
	idx = lep_idx.at(i);
      }
    }
    //cout<<"Trying to remove "<<idx <<endl;
    temp.erase(std::remove(temp.begin(), temp.end(), idx), temp.end());
    lep_idx_sorted.push_back(idx);
    //cout<<"temp size is "<< temp.size()<<endl;
  }

  if(lep_idx_sorted.size() != lep_idx.size())cout<<"Sorted size "<<lep_idx_sorted.size()<< " original size = "<<lep_idx.size()<<endl;

  //lep_idx.clear();
  std::copy(lep_idx_sorted.begin(), lep_idx_sorted.end(), lep_idx.begin());

  if(lep_idx_sorted.size() != lep_idx.size())cout<<"Sorted size "<<lep_idx_sorted.size()<< " original size = "<<lep_idx.size()<<endl;

  std::map<int, std::vector<TString> > trigmatches;
  for(unsigned int i = 0; i < lep_idx.size(); i++){
    if(i < lep_idx.size()-1)if(ptleptons_Nominal[lep_idx.at(i)]<ptleptons_Nominal[lep_idx.at(i+1)])cout<<"-----> ERROR \t Leptons not sorted!!!"<<endl;

    trigmatches[lep_idx.at(i)] = checkTriggerMatch(lep_idx.at(i),yr_num);
    
    if(trigmatches[lep_idx.at(i)].size() > 0){
      lep_tm_idx.push_back(lep_idx.at(i));
      h_lep_tm->Fill(1.0,(double)i);
    }else{
      h_lep_tm->Fill(0.0,(double)i);
    }
  }

  /**
  for(std::map<int, std::vector<TString> >::iterator iter = trigmatches.begin(); iter != trigmatches.end(); ++iter)
    {
      printf("Lepton %i (pT = %.2f, typ = %s) matched to \n",iter->first,ptleptons_Nominal[iter->first]/1000.,fabs(flavlep_Nominal[iter->first]) == 13 ? "mu" : "el");
      for(std::vector<TString>::iterator iter2 = (iter->second).begin(); iter2 != (iter->second).end(); ++iter2){
	printf("%s\n",(*iter2).Data());
      }
      printf("\n"); 
    }
  printf("\n");
  */
  
  for(unsigned int i = 0; i < lep_sg_idx.size(); i++){
    if((checkTriggerMatch(lep_sg_idx.at(i),yr_num)).size() > 0)lep_sg_tm_idx.push_back(lep_sg_idx.at(i));
  }
  
  h_lep_count->Fill(0.0,(double)lep_idx.size());
  h_lep_count->Fill(1.0,(double)lep_tm_idx.size());
  h_lep_count->Fill(2.0,(double)n_signal);
  h_lep_count->Fill(3.0,(double)lep_sg_tm_idx.size());
  
  if(n_baseline < 2)return kTRUE;

  h_cutflow->AddBinContent(3,1);


  idx1 = lep_idx.at(0);
  idx2 = lep_idx.at(1);

  for(unsigned int i = 0; i < ptjets_Nominal.GetSize(); i++){
    if(bjet(i)){
      n_bjets += 1;
      continue;
    }
    if(!jetL(i))continue;
    n_baseline_jets += 1;
    if(!jetT(i))continue;
    n_signal_jets += 1;
  }
  

  //cout<<"year = "<<yr_str.Data()<<endl;

  if(checkEventTrigger(idx1,yr_num))h_cutflow->AddBinContent(4,1);
  else if(checkEventTrigger(idx2,yr_num))h_cutflow->AddBinContent(4,1);

  lep1_tmatch = getTriggerCat(idx1,yr_num);
  lep2_tmatch = getTriggerCat(idx2,yr_num);
  if(lep1_tmatch.size() > 1){
    // cout<<"WARNING \t Lepton matched in "<< lep1_tmatch.size() <<" categories : \n";
    // for(unsigned int i = 0; i < lep1_tmatch.size(); i++)cout<<lep1_tmatch.at(i)<<endl;
    for(unsigned int i = 0; i < trig_order.size(); i++){
      if(std::count(lep1_tmatch.begin(), lep1_tmatch.end(), trig_order.at(i))){
	triggercat1 = trig_order.at(i);
	break;
      }
    }
    //cout<<"Choosing "<<triggercat1.Data()<<endl;
  }else if(lep1_tmatch.size() == 1){
    triggercat1 = lep1_tmatch.at(0);
  }else{
    lep1_tmatch.push_back("notrigm");
    triggercat1 = lep1_tmatch.at(0);
  }
  
  if(lep2_tmatch.size() > 1){
    //cout<<"WARNING \t Lepton matched in "<< lep2_tmatch.size() <<" categories : \n";
    // for(unsigned int i = 0; i < lep2_tmatch.size(); i++)cout<<lep2_tmatch.at(i)<<endl;
    for(unsigned int i = 0; i < trig_order.size(); i++){
      if(std::count(lep2_tmatch.begin(), lep2_tmatch.end(), trig_order.at(i))){
	triggercat2 = trig_order.at(i);
	break;
      }
    }
    //cout<<"Choosing "<<triggercat2.Data()<<endl;
  }else if(lep2_tmatch.size() == 1){
    triggercat2 = lep2_tmatch.at(0);
  }else{
    lep2_tmatch.push_back("notrigm");
    triggercat2 = lep2_tmatch.at(0);
  }
  // if(lep2_tmatch.size() > 1){
  //   cout<<"WARNING \t Lepton matched in "<< lep2_tmatch.size() <<" categories :";
  //   for(unsigned int i = 0; i < lep2_tmatch.size(); i++)cout<<lep2_tmatch.at(i)<<endl;
  // }
  
  // if(lep1_tmatch.size() == 0){
  // //cout<<"No triggers matched lep1. Adding notrigm "<<endl;
  //   lep1_tmatch.push_back("notrigm");
  // }
  // if(lep2_tmatch.size() == 0){
  //   //cout<<"No triggers matched lep2. Adding notrigm "<<endl;
  //   lep2_tmatch.push_back("notrigm");
  // }

  TLorentzVector lep1;
  TLorentzVector lep2;
  // In GEV!!
  lep1.SetPtEtaPhiM(ptleptons_Nominal[idx1]/1000.,etaleptons_Nominal[idx1],phileptons_Nominal[idx1],massleptons_Nominal[idx1]/1000.);
  lep2.SetPtEtaPhiM(ptleptons_Nominal[idx2]/1000.,etaleptons_Nominal[idx2],phileptons_Nominal[idx2],massleptons_Nominal[idx2]/1000.);
  Double_t mll = (lep1+lep2).M();

  TLorentzVector met_tlv;
  met_tlv.SetPtEtaPhiM((*EtMiss_tst_Nominal)/1000.,0.0,(*EtMiss_tstPhi_Nominal),0);

  std::pair <TString,TString> ret = getFillString(idx1, idx2);
  eventStr = ret.first + ret.second;

  if(eventStr.Contains("EE"))isEE = true;
  else if(eventStr.Contains("MM"))isMM = true;
  else if(eventStr.Contains("EM"))isEM = true;

  if(!isEE && !isMM && !isEM)printf("ERROR \t Event is not configured into EE, EM or MM!!\n");

  if(eventStr.Contains("OS"))isOS = true;
  else if(eventStr.Contains("SS"))isSS = true;

  if(!isOS && !isSS)printf("ERROR \t Event is not configured into OS or SS!!\n");

  if(isEE)h_triggermatched_el1_el2->Fill(checkTriggerMatch(idx1,yr_num).size() > 0 ? 1 : 0,checkTriggerMatch(idx2,yr_num).size() > 0 ? 1 : 0,1);
  else if(isMM)h_triggermatched_mu1_mu2->Fill(checkTriggerMatch(idx1,yr_num).size() > 0 ? 1 : 0,checkTriggerMatch(idx2,yr_num).size() > 0 ? 1 : 0,1);
  else h_triggermatched_el_mu->Fill(checkTriggerMatch(idx1,yr_num).size() > 0 ? 1 : 0,checkTriggerMatch(idx2,yr_num).size() > 0 ? 1 : 0,1);

  
  
  if(!(mll > 11 && lep1.Pt() > 27 && lep2.Pt() > 9 && *EtMiss_SigObj_Nominal > 3))return kTRUE;
  h_cutflow->AddBinContent(5,1);
  if((checkTriggerMatch(idx1,yr_num).size() <= 0 && checkTriggerMatch(idx2,yr_num).size() <= 0))return kTRUE;
  h_cutflow->AddBinContent(6,1);
  if(checkTriggerMatch(idx1,yr_num).size() <= 0 && checkTriggerMatch(idx2,yr_num).size() > 0 && ptleptons_Nominal[idx2] <= 27000)return kTRUE;
  h_cutflow->AddBinContent(7,1);
  
  classifyEvent(idx1, idx2, nTT, nTl, nlT, nll);
  
  gfw->setVariables(yr_str, "FCL", eventStr);
  gfw->setFrac(*EtMiss_SigObj_Nominal, mll, n_signal_jets, n_bjets, zmass, *EtMiss_SigObj_Nominal, lep1, lep2, met_tlv);

  //getTriggerMatch(idx1,idx2,yr_num);
  

  //TString triggerm = "notrigm";
  //for(std::vector< TString >::iterator it = uncert_keys.begin(); it != uncert_keys.end(); ++it){
  //for (E e = E::Begin; e != E::End; ++e) {

  int trig_idx_1 = 5;
  int trig_idx_2 = 5;
  auto trig_it_1 = find(looseTightDef.begin(), looseTightDef.end(), triggercat1);
  auto trig_it_2 = find(looseTightDef.begin(), looseTightDef.end(), triggercat2);
  if(trig_it_1 != looseTightDef.end())trig_idx_1 = trig_it_1 - looseTightDef.begin();
  if(trig_it_2 != looseTightDef.end())trig_idx_2 = trig_it_2 - looseTightDef.begin();
  
  h_lep_fake_trig_cat->Fill(trig_idx_1, trig_idx_2, 1.0);

  if(/**triggercat1.Contains("tight_iso") || triggercat2.Contains("tight_iso") || */(fabs(flavlep_Nominal[idx1]) == 11 && (fabs(etaleptons_Nominal[idx1]) > 1.37 && fabs(etaleptons_Nominal[idx1]) < 1.52)) || (fabs(flavlep_Nominal[idx2]) == 11 && (fabs(etaleptons_Nominal[idx2]) > 1.37 && fabs(etaleptons_Nominal[idx2]) < 1.52)))return kTRUE;
  
  for (auto u : GetFakeWeight::all_unc) {
    
    r1 = gfw->getReal(ptleptons_Nominal[idx1]/1000., etaleptons_Nominal[idx1], fabs(flavlep_Nominal[idx1]) == 11 ? 1 : 2, triggercat1,u);
    r2 = gfw->getReal(ptleptons_Nominal[idx2]/1000., etaleptons_Nominal[idx2], fabs(flavlep_Nominal[idx2]) == 11 ? 1 : 2, triggercat2,u);
    f1 = gfw->getFake(ptleptons_Nominal[idx1]/1000., etaleptons_Nominal[idx1], fabs(flavlep_Nominal[idx1]) == 11 ? 1 : 2, triggercat1,u);
    f2 = gfw->getFake(ptleptons_Nominal[idx2]/1000., etaleptons_Nominal[idx2], fabs(flavlep_Nominal[idx2]) == 11 ? 1 : 2, triggercat2,u);

    f1 = f1 <= 0.0 ? 0.01 : f1;
    f2 = f2 <= 0.0 ? 0.01 : f2;
    r1 = r1 >= 1.0 ? 0.99 : r1;
    r2 = r2 >= 1.0 ? 0.00 : r2;

    if((gfw->getUncKey(u)).EqualTo("NOM")){
      for (TString bdt_name : bdt_vec) {
	h_lep_BDT_n[bdt_name]->Fill(BDTweight[bdt_name],1.0);
	h_lep_BDT_r[bdt_name]->Fill(BDTweight[bdt_name],r1,1.0);
	h_lep_BDT_r[bdt_name]->Fill(BDTweight[bdt_name],r2,1.0);
      }

      if(nTT){h_r1_f1_TT->Fill(r1,f1,1.0);h_r2_f2_TT->Fill(r2,f2,1.0);}
      if(nTl){h_r1_f1_Tl->Fill(r1,f1,1.0);h_r2_f2_Tl->Fill(r2,f2,1.0);}
      if(nlT){h_r1_f1_lT->Fill(r1,f1,1.0);h_r2_f2_lT->Fill(r2,f2,1.0);}
      if(nll){h_r1_f1_ll->Fill(r1,f1,1.0);h_r2_f2_ll->Fill(r2,f2,1.0);}

      if((r1-f1) < 0.2){
	/**
	printf("r1 - f1 = %.2f - %.2f = %.2f\n",r1,f1,(r1-f1));
	printf("Lepton %i (pT = %.2f, eta = %.2f, typ = %s) matched to %s\n",idx1,ptleptons_Nominal[idx1]/1000.,etaleptons_Nominal[idx1],fabs(flavlep_Nominal[idx1]) == 13 ? "mu" : "el",triggercat1.Data());
	printf("nTT = %i, nTl = %i, nlT = %i, nll = %i\n",nTT,nTl,nlT,nll);
	*/
	h_eta_smalldiff->Fill(etaleptons_Nominal[idx1]);
	auto it = find(looseTightDef.begin(), looseTightDef.end(), triggercat1);
	if(it != looseTightDef.end())
	  {
	    int index = it - looseTightDef.begin();
	    h_tmatch_smalldiff->AddBinContent(index+1,1.0);
	  }
      }
      if((r2-f2) < 0.2){
	/**
	printf("r1 - f1 = %.2f - %.2f = %.2f\n",r2,f2,(r2-f2));
	printf("Lepton %i (pT = %.2f, eta = %.2f, typ = %s) matched to %s\n",idx2,ptleptons_Nominal[idx2]/1000.,etaleptons_Nominal[idx2],fabs(flavlep_Nominal[idx2]) == 13 ? "mu" : "el",triggercat2.Data());
	printf("nTT = %i, nTl = %i, nlT = %i, nll = %i\n",nTT,nTl,nlT,nll);
	*/
	h_eta_smalldiff->Fill(etaleptons_Nominal[idx2]);
	auto it = find(looseTightDef.begin(), looseTightDef.end(), triggercat2);
	if(it != looseTightDef.end())
	  {
	    int index = it - looseTightDef.begin();
	    h_tmatch_smalldiff->AddBinContent(index+1,1.0);
	  }
      }
      
      
    }
    if(((f1 > r1) || (f2 > r2)) && !(gfw->getUncKey(u)).EqualTo("NOM")){
      //printf("ERROR \t f1 = %.2f > r1 = %.2f\n",f1,r1);
      //printf("ERROR \t f2 = %.2f > r2 = %.2f\n",f2,r2);
      r_wgt[gfw->getUncKey(u)]  = r_wgt["NOM"];
      rf_wgt[gfw->getUncKey(u)] = rf_wgt["NOM"];
      fr_wgt[gfw->getUncKey(u)] = fr_wgt["NOM"];
      f_wgt[gfw->getUncKey(u)]  = f_wgt["NOM"];
    }else{
      r_wgt[gfw->getUncKey(u)]  = mtx->N4_RR_TT(r1,f1,r2,f2,nTT,nTl,nlT,nll);
      rf_wgt[gfw->getUncKey(u)] = mtx->N4_RF_TT(r1,f1,r2,f2,nTT,nTl,nlT,nll);
      fr_wgt[gfw->getUncKey(u)] = mtx->N4_FR_TT(r1,f1,r2,f2,nTT,nTl,nlT,nll);
      f_wgt[gfw->getUncKey(u)]  = mtx->N4_FF_TT(r1,f1,r2,f2,nTT,nTl,nlT,nll);
    }
    
    // *WeightEvents = b_fake_wgt;
    b_fake_wgt[gfw->getUncKey(u)] = f_wgt[gfw->getUncKey(u)]+fr_wgt[gfw->getUncKey(u)]+rf_wgt[gfw->getUncKey(u)];


    //cout<<"fake_wgt["<<gfw->getUncKey(u)<<"] = "<<b_fake_wgt[gfw->getUncKey(u)]<<endl;

    //if(((f1 > r1) || (f2 > r2)) && (gfw->getUncKey(u)).EqualTo("NOM")){
    if(1){
    //if(*RunNumber == 338377 && *EventNumber == 1634701019){
    //if(gfw->getUncKey(u).EqualTo("NOM")){
    //if((fabs(b_fake_wgt[gfw->getUncKey(u)]-b_fake_wgt["NOM"])/b_fake_wgt["NOM"]) > 1 && (ptleptons_Nominal[idx1] > 27000 && ptleptons_Nominal[idx2] > 9000 && mll > 11 && n_bjets == 0 && *EtMiss_SigObj_Nominal > 3)){
      gfw->setVB(1);
      cout<<"----------------------- !! START: LARGE FAKE WEIGHT "<<gfw->getUncKey(u)<<" !! ---------------------------"<<endl;

      cout<<"RunNumber = "<<*RunNumber<<endl;
      cout<<"EventNumber = "<<*EventNumber<<endl;
      
      cout<<"fake_wgt["<<gfw->getUncKey(u)<<"] = "<<b_fake_wgt[gfw->getUncKey(u)]<<endl;
      cout<<"fake_wgt[NOMINAL] = "<<b_fake_wgt["NOM"]<<endl;
 
      cout<<"Uncertainty = "<<(fabs(b_fake_wgt[gfw->getUncKey(u)]-b_fake_wgt["NOM"])/b_fake_wgt["NOM"])*100.<<"%%"<<endl;
      
      classifyEvent(idx1, idx2, nTT, nTl, nlT, nll, 1);

      printf("triggercat1 = %s\n",triggercat1.Data());
      printf("triggercat2 = %s\n",triggercat2.Data());
    
      r1 = gfw->getReal(ptleptons_Nominal[idx1]/1000., etaleptons_Nominal[idx1], fabs(flavlep_Nominal[idx1]) == 11 ? 1 : 2, triggercat1,u);
      r2 = gfw->getReal(ptleptons_Nominal[idx2]/1000., etaleptons_Nominal[idx2], fabs(flavlep_Nominal[idx2]) == 11 ? 1 : 2, triggercat2,u);
      f1 = gfw->getFake(ptleptons_Nominal[idx1]/1000., etaleptons_Nominal[idx1], fabs(flavlep_Nominal[idx1]) == 11 ? 1 : 2, triggercat1,u);
      f2 = gfw->getFake(ptleptons_Nominal[idx2]/1000., etaleptons_Nominal[idx2], fabs(flavlep_Nominal[idx2]) == 11 ? 1 : 2, triggercat2,u);

      f1 = f1 <= 0.0 ? 0.01 : f1;
      f2 = f2 <= 0.0 ? 0.01 : f2;
      r1 = r1 >= 1.0 ? 0.99 : r1;
      r2 = r2 >= 1.0 ? 0.99 : r2;

      printf("f1 = %.4f, f2 = %.4f, r1 = %.4f, r2 = %.4f\n",f1,f2,r1,r2);
      
      r_wgt[gfw->getUncKey(u)]  = mtx->N4_RR_TT(r1,f1,r2,f2,nTT,nTl,nlT,nll);
      rf_wgt[gfw->getUncKey(u)] = mtx->N4_RF_TT(r1,f1,r2,f2,nTT,nTl,nlT,nll);
      fr_wgt[gfw->getUncKey(u)] = mtx->N4_FR_TT(r1,f1,r2,f2,nTT,nTl,nlT,nll);
      f_wgt[gfw->getUncKey(u)]  = mtx->N4_FF_TT(r1,f1,r2,f2,nTT,nTl,nlT,nll);

      printf("Lepton %i (pT = %.2f, typ = %s) matched to \n",idx1,ptleptons_Nominal[idx1]/1000.,fabs(flavlep_Nominal[idx1]) == 13 ? "mu" : "el");
      for(std::vector<TString>::iterator iter2 = trigmatches[idx1].begin(); iter2 != trigmatches[idx1].end(); ++iter2){
	printf("%s\n",(*iter2).Data());
      }
      printf("Lepton %i (pT = %.2f, typ = %s) matched to \n",idx2,ptleptons_Nominal[idx2]/1000.,fabs(flavlep_Nominal[idx2]) == 13 ? "mu" : "el");
      for(std::vector<TString>::iterator iter2 = trigmatches[idx2].begin(); iter2 != trigmatches[idx2].end(); ++iter2){
	printf("%s\n",(*iter2).Data());
      }

      if(f1 > r1){
	printf("ERROR \t f1 = %.2f > r1 = %.2f\n",f1,r1);
      }
      if(f2 > r2){
	printf("ERROR \t f2 = %.2f > r2 = %.2f\n",f2,r2);
      }
    
      printf("f_hf = %.2f, f_lf = %.2f, f_co = %.2f\n",gfw->getFracHF(),gfw->getFracLF(),gfw->getFracCO());
      printf("r_wgt = %.2f, rf_wgt = %.2f, fr_wgt = %.2f, f_wgt = %.2f\n",r_wgt[gfw->getUncKey(u)],rf_wgt[gfw->getUncKey(u)],fr_wgt[gfw->getUncKey(u)],f_wgt[gfw->getUncKey(u)]);
      printf("nTT = %i, nTl = %i, nlT = %i, nll = %i\n",nTT,nTl,nlT,nll);
  
      cout<<"----------------------- !! STOP: LARGE FAKE WEIGHT "<<gfw->getUncKey(u)<<" !! ---------------------------"<<endl; 
      gfw->setVB(0);
    }
  }
  // b_WeightEvents = b_fake_wgt;
  // //cout<<"Filling "<<endl; 
  // bpt->Fill();
  // newBranch->Fill(); 
 
  // if (TF->GetEntryWithIndex(*RunNumber,*EventNumber) > 0) {
  //TF->Fill();
// }
  sumofweights += b_fake_wgt["NOM"];
  nevpass += 1;
  if(n_bjets == 1){
    nevpassnbjet1 += 1;
  }//else{
   // cout<<"n-bjet = "<<n_bjets<<endl;
  //}
  if(isEM)nevpassem += 1;
  //else cout<<"isEE = "<<(int)isEE<<", isMM = "<<(int)isMM<<endl;
  
  h_fake_weight->Fill(b_fake_wgt["NOM"]);

  if(ptleptons_Nominal[idx1] > 27000 && ptleptons_Nominal[idx2] > 9000 && mll > 11 && *EtMiss_SigObj_Nominal > 3){
    //if(isEM || ((isEE || isMM) && (mll > (91.2+15) || mll < (91.2-15)))){
    fill2L0JTree(idx1, idx2, b_fake_wgt,n_bjets,n_signal_jets,n_signal,mll,isEE,isEM,isMM,isOS,isSS);
  }
  //}

  Float_t costhetastar = TMath::ATan(fabs(lep1.Eta()-lep2.Eta())/2.);

  h_costhetastar->Fill(fabs(costhetastar));

  if(ptleptons_Nominal[idx1] > 27000 && ptleptons_Nominal[idx2] > 9000 && mll > 11 && n_signal_jets == 0 && n_bjets == 0 && *EtMiss_SigObj_Nominal > 3){
    if(((isEE || isMM) && fabs(mll-zmass) > 15) || isEM){
      if(isOS){
	h_preselection_met_OS->Fill(*EtMiss_SigObj_Nominal,1);
	h_preselection_met_OS_FNP->Fill(*EtMiss_SigObj_Nominal,b_fake_wgt["NOM"]);
	if(isEE){
	  if(n_signal == 2)h_preselection_lep2pt_EE_OS->Fill(ptleptons_Nominal[idx2]/1000,1);
	  h_preselection_lep2pt_EE_OS_FNP->Fill(ptleptons_Nominal[idx2]/1000,b_fake_wgt["NOM"]);
	}
	if(isEM){
	  if(n_signal == 2)h_preselection_lep2pt_EM_OS->Fill(ptleptons_Nominal[idx2]/1000,1);
	  h_preselection_lep2pt_EM_OS_FNP->Fill(ptleptons_Nominal[idx2]/1000,b_fake_wgt["NOM"]);
	}
      }
      if(isSS){
	h_preselection_met_SS->Fill(*EtMiss_SigObj_Nominal,1);
	h_preselection_met_SS_FNP->Fill(*EtMiss_SigObj_Nominal,b_fake_wgt["NOM"]);
      }
    }
  }

  if(n_signal_jets == 1){
    h_cutflow->AddBinContent(8,1);
    if(n_bjets == 0){
      h_cutflow->AddBinContent(9,1);
      if((isEE || isMM)){
	h_cutflow->AddBinContent(10,1);
	if(isOS){
	  h_cutflow->AddBinContent(11,1);
	  if(ptleptons_Nominal[idx1] > 100000){
	    h_cutflow->AddBinContent(12,1);
	    if(ptleptons_Nominal[idx2] > 50000){
	      h_cutflow->AddBinContent(13,1);
	      if(*EtMiss_SigObj_Nominal > 7){
		h_cutflow->AddBinContent(14,1);
		if(mll > 60){
		  h_cutflow->AddBinContent(15,1);
		  if(fabs(mll-zmass) > 15){
		    h_cutflow->AddBinContent(16,1);
		    if(costhetastar < 0.1){
		      h_cutflow->AddBinContent(17,1);
		      if(lep1.DeltaPhi(lep2)>2.8){
			h_cutflow->AddBinContent(18,1);
			h_cutflow->AddBinContent(19,b_fake_wgt["NOM"]);
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  // if(n_signal_jets == 1 && n_bjets == 0 && (isEM || isMM) && isOS && ptleptons_Nominal[idx1] > 100000 &&
  //    ptleptons_Nominal[idx2] > 50000 && *EtMiss_SigObj_Nominal > 7 && mll > 60 && fabs(mll-zmass) > 15 &&
  //    costhetastar < 0.1 && lep1.DeltaPhi(lep2)){
  //   h_cutflow->AddBinContent(7,1);
    
  // }
   
  return kTRUE;
}

void myTree::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void myTree::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  // If last file need to save branch here
  cout<<"Terminate: Opening "<<filename<<" for updating ttree"<<endl;
  // TFile *ff = new TFile(filename,"update");
  // cout<<"Wrting branch to ttree"<<endl;
  // T->Write();
  // cout<<"Closing and deleteing file"<<endl;
  // ff->Close();
  // delete ff;
  //cout<<"Flush one basket"<<endl;
  //bpt->ResetAddress();


  cout<<"sumofweights = "<<sumofweights<<endl;
  cout<<"sumofweights_up = "<<sumofweights_up<<endl;
  cout<<"sumofweights_dw = "<<sumofweights_dw<<endl;
  cout<<"nevpass = "<<nevpass<<endl;
  cout<<"nevpassnbjet1 = "<<nevpassnbjet1<<endl;
  cout<<"nevpassem = "<<nevpassem<<endl;

  cout<<"Starting writing histograms to file"<<endl;
  TFile *fff = new TFile(Form("histograms_%i.root",yr),"RECREATE");
  fff->cd();

  TIterator *nextObject = (TIterator*)fOutput->MakeIterator();
  TObject *tobj;
  int i = 0;

  while((tobj= (TObject*)nextObject->Next())) {
      
    TString name = tobj->GetName();
      
    //cout<<"Writing "<<name.Data()<<endl;
      
    tobj->Write();
      
    i += 1;
      
  }
    
  /**
     for(TString bdt_name : bdt_vec) {
     h_lep_BDT_r[bdt_name]->Write();
     h_lep_BDT_n[bdt_name]->Write();
     }
     h_triggermatched_mu->Write();
     h_triggermatched_el->Write();
     h_fake_weight->Write();
     h_triggermatched_el1_el2->Write();
     h_triggermatched_mu1_mu2->Write();
     h_triggermatched_el_mu->Write();
     h_cutflow->Write();
     h_lep_count->Write();
     h_lep_tm->Write();
     h_costhetastar->Write();
     h_preselection_met_OS->Write();
     h_preselection_met_OS_FNP->Write();
     h_preselection_lep2pt_EE_OS->Write();
     h_preselection_lep2pt_EE_OS_FNP->Write();
     h_preselection_lep2pt_EM_OS->Write();
     h_preselection_lep2pt_EM_OS_FNP->Write();
     h_preselection_met_SS->Write();    
     h_preselection_met_SS_FNP->Write();
     h_lep_fake_trig_cat->Write();
  */
  fff->Close();

  delete twoLzeroJ;
}

bool myTree::jetL(int idx){
  if(ptjets_Nominal[idx] <= 20000)return false;
  if(fabs(etajets_Nominal[idx]) > 2.8)return false;
  if(!passORjet_Nominal->at(idx))return false;
  return true;
}

bool myTree::jetT(int idx){
  if(fabs(etajets_Nominal[idx]) > 2.4)return false;
  if(ptjets_Nominal[idx] < 60000 && JVTjets_Nominal[idx] <= 0.91)return false;
  return true;
}

bool myTree::bjet(int idx){
  if(!jetL(idx) || fabs(etajets_Nominal[idx]) > 2.4)return false;
  if(!isBjets_Nominal[idx])return false;
  return true;
}

bool myTree::lepL(int idx){
  if(ptleptons_Nominal[idx] <= 9000)return false;
  if(!passORlep_Nominal->at(idx))return false;
  if(fabs(z0sinTheta_Nominal[idx]) > 0.5)return false;
  if(fabs(flavlep_Nominal[idx]) == 13 && fabs(etaleptons_Nominal[idx]) >= 2.6)return false;
  if(fabs(flavlep_Nominal[idx]) == 11 && fabs(etaleptons_Nominal[idx]) >= 2.47)return false;
  return true;
}

bool myTree::lepT(int idx, int vb){
  if(!lepL(idx)){if(vb){cout<<"failed bl selection"<<endl;}return false;}
  if(!isSignallep_Nominal->at(idx))return false;
  // As of January 2021 these are no longer needed. We use the signal lept flag instead.
  //if(fabs(flavlep_Nominal[idx]) == 13 && (!isMediumlep_Nominal->at(idx) || !LepIsoFCLoose_Nominal->at(idx) || d0sig_Nominal[idx] >= 3)){if(vb){printf("isMediumlep = %i, LepIsoFCLoose = %i, d0sig = %.2f\n",(int)isMediumlep_Nominal->at(idx),(int)LepIsoFCLoose_Nominal->at(idx),d0sig_Nominal[idx]);}return false;}
  //if(fabs(flavlep_Nominal[idx]) == 11 && (!isTightlep_Nominal->at(idx) || !LepIsoFCLoose_Nominal->at(idx) || d0sig_Nominal[idx] >= 5)){if(vb){printf("isTightlep_Nominal = %i, LepIsoFCLoose = %i, d0sig = %.2f\n",(int)isTightlep_Nominal->at(idx),(int)LepIsoFCLoose_Nominal->at(idx),d0sig_Nominal[idx]);}return false;}
  //if(fabs(flavlep_Nominal[idx]) == 13 && fabs(etaleptons_Nominal[idx]) >= 2.4){if(vb){cout<<"failed eta selection"<<endl;}return false;}
  return true;
}

std::pair <TString,TString> myTree::getFillString(int idx1, int idx2)
{
  TString chstr = "";
  TString sgnstr = "";

  
  if(fabs(flavlep_Nominal[idx1]) == 11 && fabs(flavlep_Nominal[idx2]) == 11){
    chstr = "EE";
  }else if(fabs(flavlep_Nominal[idx1]) == 13 && fabs(flavlep_Nominal[idx2]) == 13){
    chstr = "MM";
  }else if(fabs(flavlep_Nominal[idx1]) == 11 && fabs(flavlep_Nominal[idx2]) == 13){
    chstr = "EM";
  }else if(fabs(flavlep_Nominal[idx1]) == 13 && fabs(flavlep_Nominal[idx2]) == 11){
    chstr = "EM";
  }else{
    cout<<"ERROR \t Strange error flav1 = "<<fabs(flavlep_Nominal[idx1])<<", flav2 = "<<fabs(flavlep_Nominal[idx2])<<endl;
  }
  

  if(flavlep_Nominal[idx1]*flavlep_Nominal[idx2] > 0)sgnstr = "SS";
  else if(flavlep_Nominal[idx1]*flavlep_Nominal[idx2] < 0)sgnstr = "OS";
  
  
  std::pair <TString,TString> ret (chstr,sgnstr);
  
  return ret;
}

void myTree::classifyEvent(int idx1, int idx2, int& isTT, int& isTL, int& isLT, int& isLL, int vb)
{

  isTT = 0;
  isTL = 0;
  isLT = 0;
  isLL = 0;
  if(lepT(idx1,vb) && lepT(idx2,vb))isTT = 1;
  else if(lepT(idx1,vb) && !lepT(idx2,vb))isTL = 1;
  else if(!lepT(idx1,vb) && lepT(idx2,vb))isLT = 1;
  else if(!lepT(idx1,vb) && !lepT(idx2,vb))isLL = 1;

}


bool myTree::checkEventTrigger(int id, int yr){
  //if(fabs(flavlep_Nominal[id]) == 11){
  for(unsigned int i=0; i<triglist_e[yr].size();i++){
    //printf("%s : %s\n",trigstr_e[yr].at(i).Data(),*triglist_e[yr].at(i) ? "True" : "False");
    if(*triglist_e[yr].at(i))return true;
  }
  //}else if(fabs(flavlep_Nominal[id]) == 13){
  for(unsigned int i=0; i<triglist_m[yr].size();i++){
    //printf("%s : %s\n",trigstr_m[yr].at(i).Data(),*triglist_m[yr].at(i) ? "True" : "False");
    if(*triglist_m[yr].at(i))return true;
  }
  //}
  return false;
}




std::vector<TString> myTree::checkTriggerMatch(int id, int yr){
  std::vector<TString> matched_triggers; 
  if(fabs(flavlep_Nominal[id]) == 11){
    for(unsigned int i=0; i<trigmatch_e[yr].size();i++){
      //printf("Threshold for trigger %s is %i\n",trigstr_e[yr].at(i).Data(),getPtThresholdTrigger(trigstr_e[yr].at(i),0));
      //printf("El(%.2f): %s : %s\n",ptleptons_Nominal[id]/1000.,trigstr_e[yr].at(i).Data(),((TTreeReaderArray<int>*)trigmatch_e[yr].at(i))->At(id) ? "True" : "False");
      if(((TTreeReaderArray<int>*)trigmatch_e[yr].at(i))->At(id) && ptleptons_Nominal[id]/1000. > getPtThresholdTrigger(trigstr_e[yr].at(i),0))matched_triggers.push_back(trigstr_e[yr].at(i));
    }
  }else if(fabs(flavlep_Nominal[id]) == 13){
    for(unsigned int i=0; i<trigmatch_m[yr].size();i++){
      //printf("Threshold for trigger %s is %i\n",trigstr_m[yr].at(i).Data(),getPtThresholdTrigger(trigstr_m[yr].at(i),0));
      
      //printf("Mu(%.2f): %s : %s\n",ptleptons_Nominal[id]/1000.,trigstr_m[yr].at(i).Data(),((TTreeReaderArray<int>*)trigmatch_m[yr].at(i))->At(id) ? "True" : "False");
      if(((TTreeReaderArray<int>*)trigmatch_m[yr].at(i))->At(id) && ptleptons_Nominal[id]/1000. > getPtThresholdTrigger(trigstr_m[yr].at(i),0))matched_triggers.push_back(trigstr_m[yr].at(i));
    }
  }
  return matched_triggers;
}

std::vector<TString> myTree::getTriggerCat(int id, int yr){
  std::vector<TString> trigc; 
  if(fabs(flavlep_Nominal[id]) == 11){
    for(unsigned int i=0; i<trigmatch_e[yr].size();i++){
      if(((TTreeReaderArray<int>*)trigmatch_e[yr].at(i))->At(id) && ptleptons_Nominal[id]/1000. > getPtThresholdTrigger(trigstr_e[yr].at(i),0)){
	//printf("El(%.2f): %s : %s\n",ptleptons_Nominal[id]/1000.,trigstr_e[yr].at(i).Data(),((TTreeReaderArray<int>*)trigmatch_e[yr].at(i))->At(id) ? "True" : "False");
	h_triggermatched_el->AddBinContent(i+2,1);
      }
      if(((TTreeReaderArray<int>*)trigmatch_e[yr].at(i))->At(id) && ptleptons_Nominal[id]/1000. > getPtThresholdTrigger(trigstr_e[yr].at(i),0) && !(std::find(trigc.begin(), trigc.end(), trigcat[trigstr_e[yr].at(i)]) != trigc.end()))trigc.push_back(trigcat[trigstr_e[yr].at(i)]);
    }
    if(trigc.size() == 0)h_triggermatched_el->AddBinContent(1,1);
  }else if(fabs(flavlep_Nominal[id]) == 13){
    for(unsigned int i=0; i<trigmatch_m[yr].size();i++){
      if(((TTreeReaderArray<int>*)trigmatch_m[yr].at(i))->At(id) && ptleptons_Nominal[id]/1000. > getPtThresholdTrigger(trigstr_m[yr].at(i),0)){
	//printf("Mu(%.2f): %s : %s\n",ptleptons_Nominal[id]/1000.,trigstr_m[yr].at(i).Data(),((TTreeReaderArray<int>*)trigmatch_m[yr].at(i))->At(id) ? "True" : "False");
	h_triggermatched_mu->AddBinContent(i+2,1);
      }
      if(((TTreeReaderArray<int>*)trigmatch_m[yr].at(i))->At(id)  && ptleptons_Nominal[id]/1000. > getPtThresholdTrigger(trigstr_m[yr].at(i),0) && !(std::find(trigc.begin(), trigc.end(), trigcat[trigstr_m[yr].at(i)]) != trigc.end()))trigc.push_back(trigcat[trigstr_m[yr].at(i)]);
    }
    if(trigc.size() == 0)h_triggermatched_mu->AddBinContent(1,1);   
  }
  return trigc;
}

Int_t myTree::getPtThresholdTrigger(TString tname, int vb){
  
  TString trig_ptcut = ((TObjString *)(tname.Tokenize("_"))->At(1))->String();
  if(trig_ptcut.Contains("e"))trig_ptcut.Remove(trig_ptcut.First("e"),1);
  else if(trig_ptcut.Contains("mu"))trig_ptcut.Remove(trig_ptcut.First("mu"),2);
  if(vb)cout<<"pT cut for trigger "<< tname.Data() <<" is "<< trig_ptcut.Data() << endl;
  Float_t trig_ptcut_val = trig_ptcut.Atoi();
  if(vb)cout<<"pT cut for trigger "<< tname.Data() <<" is "<< trig_ptcut_val << endl;

  return trig_ptcut_val;
  
}


void myTree::fill2L0JTree(int idx1, int idx2, std::map<TString,double> fake_weight, Int_t n_bjets, Int_t n_ljets, Int_t n_siglep, Float_t mll, bool isEE, bool isEM, bool isMM, bool isOS, bool isSS){

  //int i = 0;
  double tot_err_up = 0.0;
  double tot_err_dw = 0.0;
  TString dwn_str = "";
  TString up_str = "";
  double nominal_weight = fake_weight["NOM"];
  for(std::map<TString,double>::iterator iter = fake_weight.begin(); iter != fake_weight.end(); ++iter)
    {
      twoLzeroJ->b2L0J_FNPweight[iter->first] = iter->second;

      // if((iter->first).Contains("WEIGHT")){
      // 	if(fabs(iter->second) < fabs(nominal_weight)){
      // 	  tot_err_dw += fabs((iter->second)-nominal_weight)/nominal_weight;
      // 	  printf("Adding %s in lineraly DOWN\n",iter->first.Data());
      // 	  if((iter->first).Contains("up") || (iter->first).Contains("UP"))
      // 	    printf("ERROR \t %s added to DW: unc = %.2f, nom = %.2f\n",iter->first.Data(),iter->second,nominal_weight); 
      // 	}else if(fabs(iter->second) > fabs(nominal_weight)){
      // 	  tot_err_up += fabs((iter->second)-nominal_weight)/nominal_weight;;
      // 	  printf("Adding %s in lineraly UP\n",iter->first.Data());
      // 	  if((iter->first).Contains("down") || (iter->first).Contains("DW"))
      // 	    printf("ERROR \t %s added to UP: unc = %.2f, nom = %.2f\n",iter->first.Data(),iter->second,nominal_weight);
      // 	}
      // }else 
      if(!(iter->first).EqualTo("NOM")){
	if(iter->second < nominal_weight){
	  tot_err_dw += ((fabs((iter->second)-nominal_weight)/nominal_weight)*(fabs((iter->second)-nominal_weight)/nominal_weight));
	  //printf("Adding %s in quadrature DOWN\n",iter->first.Data());
	  //if((iter->first).Contains("up") || (iter->first).Contains("UP"))
	  //printf("ERROR \t %s added to DW: unc = %.2f, nom = %.2f\n",iter->first.Data(),iter->second,nominal_weight);
	  dwn_str += Form("%.4f *",(((fabs((iter->second)-nominal_weight)/nominal_weight)*(fabs((iter->second)-nominal_weight)/nominal_weight))));
	}else if(iter->second > nominal_weight){
	  tot_err_up += ((fabs((iter->second)-nominal_weight)/nominal_weight)*(fabs((iter->second)-nominal_weight)/nominal_weight));
	  //printf("Adding %s in quadrature UP\n",iter->first.Data());
	  //if((iter->first).Contains("down") || (iter->first).Contains("DW"))
	  //printf("ERROR \t %s added to UP: unc = %.2f, nom = %.2f\n",iter->first.Data(),iter->second,nominal_weight);
	  up_str += Form("%.4f *",(((fabs((iter->second)-nominal_weight)/nominal_weight)*(fabs((iter->second)-nominal_weight)/nominal_weight))));
	}
      }
      
      //i += 1;
    }

  /**
  if(sqrt(tot_err_up) > 3){
    cout<<"EventNumber: "<<*EventNumber<<endl;
    cout<<"UP: "<<up_str.Data()<<endl;
    cout<<"DW: "<<dwn_str.Data()<<endl;
    printf("TOT ERR UP = %.2f\n",sqrt(tot_err_up));
    for(std::map<TString,double>::iterator iter = fake_weight.begin(); iter != fake_weight.end(); ++iter)
      {
	printf("DEBUG \t %s =  %.2f\n",iter->first.Data(),iter->second);
      }
  }
  */
  sumofweights_up += nominal_weight + fabs(sqrt(tot_err_up)*nominal_weight);
  sumofweights_dw += nominal_weight - fabs(sqrt(tot_err_dw)*nominal_weight);


  if(doExtraVar){

    for(std::map<TString,std::vector<double>>::iterator iter = BDTscores[*RunNumber].begin(); iter != BDTscores[*RunNumber].end(); ++iter){
      //printf("Getting %s at index %i = %.2f\n",iter->first.Data(),fnpidx,(iter->second).at(fnpidx));
      twoLzeroJ->b2L0J_BDTweight[iter->first] = (iter->second).at(fnpidx);
    }
    

    /**
    std::vector<std::string> v_bkg {"BDTDeltaM100_90","BDTVVDeltaM100_90","BDTtopDeltaM100_90","BDTothersDeltaM100_90","BDTDeltaM30","BDTVVDeltaM30","BDTtopDeltaM30","BDTothersDeltaM30"};

    auto filt = rdf.Filter(Form("EventNumber == %llu && RunNumber == %u",(*EventNumber,*RunNumber)));
    cout<<"val = "<<filt.Count().GetValue()<<endl;
    if(filt.Count().GetValue() == 1){
      for (std::string bdt_score : {"BDTDeltaM100_90","BDTVVDeltaM100_90","BDTtopDeltaM100_90","BDTothersDeltaM100_90","BDTDeltaM30","BDTVVDeltaM30","BDTtopDeltaM30","BDTothersDeltaM30"}) {
	TString key = bdt_score;
	cout<<"key = "<<key<<endl;
	std::cout<<bdt_score <<" = "<<filt.Take<float>(bdt_score)->at(0)<<endl;
	twoLzeroJ->b2L0J_BDTweight[key] = filt.Take<float>(bdt_score)->at(0);
      }
    }else if(filt.Count().GetValue() > 1){
      cout<<"Found "<<filt.Count().GetValue() <<" events matching criteria"<<endl;
    }
    */
    twoLzeroJ->b2L0J_lep1Pt = ptleptons_Nominal[idx1] >= ptleptons_Nominal[idx2] ? ptleptons_Nominal[idx1] : ptleptons_Nominal[idx2];
    twoLzeroJ->b2L0J_lep2Pt = ptleptons_Nominal[idx1] >= ptleptons_Nominal[idx2] ? ptleptons_Nominal[idx2] : ptleptons_Nominal[idx1];
  
    twoLzeroJ->b2L0J_n_bl_lep = 2;
    twoLzeroJ->b2L0J_mll = mll;
    twoLzeroJ->b2L0J_n_sg_lep = n_siglep;	
    twoLzeroJ->b2L0J_n_bjet = n_bjets;	
    twoLzeroJ->b2L0J_n_ljet = n_ljets;	
    twoLzeroJ->b2L0J_isSF = isEM ? 0 : 1; //; (fabs(flavlep_Nominal[idx1]) == fabs(flavlep_Nominal[idx2])) ? 1 : 0;	
    twoLzeroJ->b2L0J_isOS = isOS ? 1 : 0;//(flavlep_Nominal[idx1]*flavlep_Nominal[idx2]) < 0 ? 1 : 0;	
    twoLzeroJ->b2L0J_isEE = isEE ? 1 : 0;//((fabs(flavlep_Nominal[idx1]) == fabs(flavlep_Nominal[idx2])) && fabs(flavlep_Nominal[idx1]) == 11) ? 1 : 0;	
    twoLzeroJ->b2L0J_isEM = isEM ? 1 : 0;//((fabs(flavlep_Nominal[idx1]) == fabs(flavlep_Nominal[idx2]))) ? 0 : 1;	    
    twoLzeroJ->b2L0J_isMM = isMM ? 1 : 0;//((fabs(flavlep_Nominal[idx1]) == fabs(flavlep_Nominal[idx2])) && fabs(flavlep_Nominal[idx1]) == 13) ? 1 : 0;   
  }
  
  twoLzeroJ->b2L0J_FNP_TOTAL_UP   = nominal_weight + fabs(sqrt(tot_err_up)*nominal_weight);
  twoLzeroJ->b2L0J_FNP_TOTAL_DOWN = nominal_weight - fabs(sqrt(tot_err_dw)*nominal_weight);
  
  twoLzeroJ->b2L0J_EventNumber                                                            = (*EventNumber);                                
  twoLzeroJ->b2L0J_RunNumber                                                              = (*RunNumber);                                   
  twoLzeroJ->b2L0J_GenFiltMET                                                             = (*GenFiltMET);                                  
  twoLzeroJ->b2L0J_GenFiltHT                                                              = (*GenFiltHT);                                   
  twoLzeroJ->b2L0J_xsec                                                                   = (*xsec);                                        
  twoLzeroJ->b2L0J_WeightEvents                                                           = (*WeightEvents);                                
  twoLzeroJ->b2L0J_WeightEventsPU_Nominal                                                 = (*WeightEventsPU_Nominal);                      
  twoLzeroJ->b2L0J_WeightEventselSF_Nominal                                               = (*WeightEventselSF_Nominal);                    
  twoLzeroJ->b2L0J_WeightEventsmuSF_Nominal                                               = (*WeightEventsmuSF_Nominal);                    
  twoLzeroJ->b2L0J_WeightEventsSF_global_Nominal                                          = (*WeightEventsSF_global_Nominal);               
  twoLzeroJ->b2L0J_WeightEvents_trigger_single_Nominal                                    = (*WeightEvents_trigger_single_Nominal);         
  twoLzeroJ->b2L0J_WeightEvents_trigger_global_Nominal                                    = (*WeightEvents_trigger_global_Nominal);         
  twoLzeroJ->b2L0J_WeightEventsbTag_Nominal                                               = (*WeightEventsbTag_Nominal);                    
  twoLzeroJ->b2L0J_WeightEventsJVT_Nominal                                                = (*WeightEventsJVT_Nominal);                     
  twoLzeroJ->b2L0J_AverageInteractionsPerCrossing                                         = (*AverageInteractionsPerCrossing);              
  twoLzeroJ->b2L0J_NumPrimaryVertices                                                     = (*NumPrimaryVertices );                         
  twoLzeroJ->b2L0J_DecayType                                                              = (*DecayType             );                      
  twoLzeroJ->b2L0J_HFclass                                                                = (*HFclass                  );                     
  twoLzeroJ->b2L0J_passMETtrig                                                            = (*passMETtrig              );                  
  twoLzeroJ->b2L0J_TruthPDGID1                                                            = (*TruthPDGID1              );                   
  twoLzeroJ->b2L0J_TruthPDGID2                                                            = (*TruthPDGID2              );                   
  twoLzeroJ->b2L0J_TruthPDFID1                                                            = (*TruthPDFID1              );                  
  twoLzeroJ->b2L0J_TruthPDFID2                                                            = (*TruthPDFID2              );                   
  twoLzeroJ->b2L0J_TruthX1                                                                = (*TruthX1                    );                
  twoLzeroJ->b2L0J_TruthX2                                                                = (*TruthX2                    );                
  twoLzeroJ->b2L0J_TruthXF1                                                               = (*TruthXF1                   );                
  twoLzeroJ->b2L0J_TruthXF2                                                               = (*TruthXF2                   );                
  twoLzeroJ->b2L0J_TruthQ                                                                 = (*TruthQ                     );

  for(unsigned int i = 0; i < labeljets_reco.GetSize(); i++){
    twoLzeroJ->b2L0J_labeljets_reco.push_back(labeljets_reco[i]);                
    twoLzeroJ->b2L0J_ptjets_Nominal.push_back(ptjets_Nominal[i]);                 
    twoLzeroJ->b2L0J_etajets_Nominal.push_back(etajets_Nominal[i]);                
    twoLzeroJ->b2L0J_phijets_Nominal.push_back(phijets_Nominal[i]);                
    twoLzeroJ->b2L0J_massjets_Nominal.push_back(massjets_Nominal[i]);               
    twoLzeroJ->b2L0J_isBjets_Nominal.push_back(isBjets_Nominal[i]);                
    twoLzeroJ->b2L0J_BtagWeightjets_Nominal.push_back(BtagWeightjets_Nominal[i]);         
    twoLzeroJ->b2L0J_JVTjets_Nominal.push_back(JVTjets_Nominal[i]);                    
    twoLzeroJ->b2L0J_passORjet_Nominal.push_back(passORjet_Nominal->at(i));                 
    twoLzeroJ->b2L0J_isHiggsjet_Nominal.push_back(isHiggsjet_Nominal->at(i));                
    twoLzeroJ->b2L0J_isWhichjet_Nominal.push_back(isWhichjet_Nominal[i]);              
  }
  for(unsigned int i = 0; i < ptleptons_Nominal.GetSize(); i++){
    twoLzeroJ->b2L0J_ptleptons_Nominal.push_back(ptleptons_Nominal[i]);                      
    twoLzeroJ->b2L0J_etaleptons_Nominal.push_back(etaleptons_Nominal[i]);                     
    twoLzeroJ->b2L0J_phileptons_Nominal.push_back(phileptons_Nominal[i]);                     
    twoLzeroJ->b2L0J_massleptons_Nominal.push_back(massleptons_Nominal[i]);                    
    twoLzeroJ->b2L0J_ptcone20_Nominal.push_back(ptcone20_Nominal[i]);                       
    twoLzeroJ->b2L0J_topoetcone20_Nominal.push_back(topoetcone20_Nominal[i]);                   
    twoLzeroJ->b2L0J_ptvarcone20_Nominal.push_back(ptvarcone20_Nominal[i]);                    
    twoLzeroJ->b2L0J_ptvarcone30_TightTTVA_pt1000_Nominal.push_back(ptvarcone30_TightTTVA_pt1000_Nominal[i]);   
    twoLzeroJ->b2L0J_ptcone20_TightTTVA_pt1000_Nominal.push_back(ptcone20_TightTTVA_pt1000_Nominal[i]);      
    twoLzeroJ->b2L0J_flavlep_Nominal.push_back(flavlep_Nominal[i]);                        
    twoLzeroJ->b2L0J_passORlep_Nominal.push_back(passORlep_Nominal->at(i));                      
    twoLzeroJ->b2L0J_isSignallep_Nominal.push_back(isSignallep_Nominal->at(i));                    
    twoLzeroJ->b2L0J_isHighPtlep_Nominal.push_back(isHighPtlep_Nominal->at(i));                    
    twoLzeroJ->b2L0J_isMediumlep_Nominal.push_back(isMediumlep_Nominal->at(i));                    
    twoLzeroJ->b2L0J_isTightlep_Nominal.push_back(isTightlep_Nominal->at(i));                     
    twoLzeroJ->b2L0J_LepIsoFCHighPtCaloOnly_Nominal.push_back(LepIsoFCHighPtCaloOnly_Nominal->at(i));         
    twoLzeroJ->b2L0J_LepIsoFixedCutHighPtTrackOnly_Nominal.push_back(LepIsoFixedCutHighPtTrackOnly_Nominal->at(i));   
    twoLzeroJ->b2L0J_LepIsoGradient_Nominal.push_back(LepIsoGradient_Nominal->at(i));                 
    twoLzeroJ->b2L0J_LepIsoFCTightTrackOnly_Nominal.push_back(LepIsoFCTightTrackOnly_Nominal->at(i));         
    twoLzeroJ->b2L0J_LepIsoFCLoose_Nominal.push_back(LepIsoFCLoose_Nominal->at(i));                  
    twoLzeroJ->b2L0J_LepIsoFCTight_Nominal.push_back(LepIsoFCTight_Nominal->at(i));                  
    twoLzeroJ->b2L0J_PromptLeptonIso_Nominal.push_back(PromptLeptonIso_Nominal[i]);                 
    twoLzeroJ->b2L0J_PromptLeptonVeto_Nominal.push_back(PromptLeptonVeto_Nominal[i]);                
    twoLzeroJ->b2L0J_d0sig_Nominal.push_back(d0sig_Nominal[i]);                           
    twoLzeroJ->b2L0J_z0sinTheta_Nominal.push_back(z0sinTheta_Nominal[i]);


    twoLzeroJ->b2L0J_HLT_2e12_lhloose_L12EM10VH_match_Nominal.push_back(HLT_2e12_lhloose_L12EM10VH_match_Nominal[i]      );    
    twoLzeroJ->b2L0J_HLT_2e15_lhvloose_nod0_L12EM13VH_match_Nominal.push_back(HLT_2e15_lhvloose_nod0_L12EM13VH_match_Nominal[i]      );    
    twoLzeroJ->b2L0J_HLT_2e17_lhvloose_nod0_L12EM15VHI_match_Nominal.push_back(HLT_2e17_lhvloose_nod0_L12EM15VHI_match_Nominal[i]       );    
    twoLzeroJ->b2L0J_HLT_2e17_lhvloose_nod0_match_Nominal.push_back(HLT_2e17_lhvloose_nod0_match_Nominal[i]         );  
    twoLzeroJ->b2L0J_HLT_2e24_lhvloose_nod0_match_Nominal.push_back(HLT_2e24_lhvloose_nod0_match_Nominal[i]         );
    twoLzeroJ->b2L0J_HLT_2mu10_match_Nominal.push_back(HLT_2mu10_match_Nominal[i]                      );  
    twoLzeroJ->b2L0J_HLT_2mu14_match_Nominal.push_back(HLT_2mu14_match_Nominal[i]                      );  
    twoLzeroJ->b2L0J_HLT_e120_lhloose_match_Nominal.push_back(HLT_e120_lhloose_match_Nominal[i]                 );
    twoLzeroJ->b2L0J_HLT_e140_lhloose_nod0_match_Nominal.push_back(HLT_e140_lhloose_nod0_match_Nominal[i]            );
    twoLzeroJ->b2L0J_HLT_e17_lhloose_2e9_lhloose_match_Nominal.push_back(HLT_e17_lhloose_2e9_lhloose_match_Nominal[i]           );
    twoLzeroJ->b2L0J_HLT_e17_lhloose_mu14_match_Nominal.push_back(HLT_e17_lhloose_mu14_match_Nominal[i]         );    
    twoLzeroJ->b2L0J_HLT_e17_lhloose_nod0_mu14_match_Nominal.push_back(HLT_e17_lhloose_nod0_mu14_match_Nominal[i]       );    
    twoLzeroJ->b2L0J_HLT_e24_lhmedium_L1EM20VH_match_Nominal.push_back(HLT_e24_lhmedium_L1EM20VH_match_Nominal[i]       );    
    twoLzeroJ->b2L0J_HLT_e24_lhmedium_nod0_L1EM20VH_match_Nominal.push_back(HLT_e24_lhmedium_nod0_L1EM20VH_match_Nominal[i]         );  
    twoLzeroJ->b2L0J_HLT_e24_lhtight_nod0_ivarloose_match_Nominal.push_back(HLT_e24_lhtight_nod0_ivarloose_match_Nominal[i]         );  
    twoLzeroJ->b2L0J_HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1_match_Nominal.push_back(HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1_match_Nominal[i]         );  
    twoLzeroJ->b2L0J_HLT_e24_lhvloose_nod0_2e12_lhvloose_nod0_L1EM20VH_3EM10VH_match_Nominal.push_back(HLT_e24_lhvloose_nod0_2e12_lhvloose_nod0_L1EM20VH_3EM10VH_match_Nominal[i]           );
    twoLzeroJ->b2L0J_HLT_e26_lhmedium_nod0_mu8noL1_match_Nominal.push_back(HLT_e26_lhmedium_nod0_mu8noL1_match_Nominal[i]            );
    twoLzeroJ->b2L0J_HLT_e26_lhtight_nod0_ivarloose_match_Nominal.push_back(HLT_e26_lhtight_nod0_ivarloose_match_Nominal[i]            );
    //twoLzeroJ->b2L0J_HLT_e300_etcut_match_Nominal.push_back(HLT_e300_etcut_match_Nominal[i]               );    
    twoLzeroJ->b2L0J_HLT_e60_lhmedium_nod0_match_Nominal.push_back(HLT_e60_lhmedium_nod0_match_Nominal[i]        );    
    twoLzeroJ->b2L0J_HLT_e60_lhmedium_match_Nominal.push_back(HLT_e60_lhmedium_match_Nominal[i]             );    
    twoLzeroJ->b2L0J_HLT_e7_lhmedium_mu24_match_Nominal.push_back(HLT_e7_lhmedium_mu24_match_Nominal[i]           );  
    twoLzeroJ->b2L0J_HLT_e7_lhmedium_nod0_mu24_match_Nominal.push_back(HLT_e7_lhmedium_nod0_mu24_match_Nominal[i]         );
    twoLzeroJ->b2L0J_HLT_mu18_mu8noL1_match_Nominal.push_back(HLT_mu18_mu8noL1_match_Nominal[i]                 );
    twoLzeroJ->b2L0J_HLT_mu20_iloose_L1MU15_match_Nominal.push_back(HLT_mu20_iloose_L1MU15_match_Nominal[i]       );    
    twoLzeroJ->b2L0J_HLT_mu20_mu8noL1_match_Nominal.push_back(HLT_mu20_mu8noL1_match_Nominal[i]             );    
    twoLzeroJ->b2L0J_HLT_mu22_mu8noL1_calotag_0eta010_match_Nominal.push_back(HLT_mu22_mu8noL1_calotag_0eta010_match_Nominal[i]      );    
    twoLzeroJ->b2L0J_HLT_mu22_mu8noL1_match_Nominal.push_back(HLT_mu22_mu8noL1_match_Nominal[i]               );  
    twoLzeroJ->b2L0J_HLT_mu24_ivarloose_L1MU15_match_Nominal.push_back(HLT_mu24_ivarloose_L1MU15_match_Nominal[i]         );  
    twoLzeroJ->b2L0J_HLT_mu24_ivarloose_match_Nominal.push_back(HLT_mu24_ivarloose_match_Nominal[i]             );  
    twoLzeroJ->b2L0J_HLT_mu24_ivarmedium_match_Nominal.push_back(HLT_mu24_ivarmedium_match_Nominal[i]              );
    twoLzeroJ->b2L0J_HLT_mu26_ivarmedium_match_Nominal.push_back(HLT_mu26_ivarmedium_match_Nominal[i]              );
    twoLzeroJ->b2L0J_HLT_mu40_match_Nominal.push_back(HLT_mu40_match_Nominal[i]                         );
    twoLzeroJ->b2L0J_HLT_mu50_match_Nominal.push_back(HLT_mu50_match_Nominal[i]                     );    
    //twoLzeroJ->b2L0J_HLT_mu60_0eta105_msonly_match_Nominal.push_back(HLT_mu60_0eta105_msonly_match_Nominal[i]       );

    
  }

  for(unsigned int i = 0; i < HLT_g140_loose_match_Nominal.GetSize(); i++){
    twoLzeroJ->b2L0J_HLT_g140_loose_match_Nominal.push_back(HLT_g140_loose_match_Nominal[i]);                
    twoLzeroJ->b2L0J_HLT_g300_etcutHLT_g35_medium_g25_medium_L12EM20VH_match_Nominal.push_back(HLT_g300_etcutHLT_g35_medium_g25_medium_L12EM20VH_match_Nominal[i]);          
    twoLzeroJ->b2L0J_HLT_2g20_tight_icalovloose_L12EM15VHI_match_Nominal.push_back(HLT_2g20_tight_icalovloose_L12EM15VHI_match_Nominal[i]);          
    twoLzeroJ->b2L0J_HLT_2g22_tight_L12EM15VHI_match_Nominal.push_back(HLT_2g22_tight_L12EM15VHI_match_Nominal[i]);     
    twoLzeroJ->b2L0J_HLT_2g25_loose_g15_loose_match_Nominal.push_back(HLT_2g25_loose_g15_loose_match_Nominal[i]);     
    twoLzeroJ->b2L0J_HLT_2g25_tight_L12EM20VH_match_Nominal.push_back(HLT_2g25_tight_L12EM20VH_match_Nominal[i]);     
    twoLzeroJ->b2L0J_HLT_2g50_loose_L12EM20VH_match_Nominal.push_back(HLT_2g50_loose_L12EM20VH_match_Nominal[i]);     
    twoLzeroJ->b2L0J_HLT_2g20_tight_match_Nominal.push_back(HLT_2g20_tight_match_Nominal[i]);     
    twoLzeroJ->b2L0J_HLT_2g20_loose_g15_loose_match_Nominal.push_back(HLT_2g20_loose_g15_loose_match_Nominal[i]);     
    twoLzeroJ->b2L0J_HLT_2g22_tight_match_Nominal.push_back(HLT_2g22_tight_match_Nominal[i]);                  

  }

  for(unsigned int i = 0; i < HLT_j80_xe80_match_Nominal.GetSize(); i++){
    twoLzeroJ->b2L0J_HLT_j80_xe80_match_Nominal.push_back(HLT_j80_xe80_match_Nominal[i]);  
  }

  twoLzeroJ->b2L0J_cleaningVeto_Nominal                                                   = (*cleaningVeto_Nominal                         );  
  twoLzeroJ->b2L0J_EtMiss_tstPhi_Nominal                                                  = (*EtMiss_tstPhi_Nominal                          );
  twoLzeroJ->b2L0J_EtMiss_tst_Nominal                                                     = (*EtMiss_tst_Nominal                             );
  twoLzeroJ->b2L0J_EtMiss_SigObj_Nominal                                                  = (*EtMiss_SigObj_Nominal                          );
  // twoLzeroJ->b2L0J_Etmiss_PVSoftTrkPhi                                                    = (*Etmiss_PVSoftTrkPhi                            );
  // twoLzeroJ->b2L0J_Etmiss_PVSoftTrk                                                       = (*Etmiss_PVSoftTrk                               );

  
  twoLzeroJ->b2L0J_HLT_2e12_lhloose_L12EM10VH_Nominal                                     = (*HLT_2e12_lhloose_L12EM10VH_Nominal             );
  twoLzeroJ->b2L0J_HLT_2e15_lhvloose_nod0_L12EM13VH_Nominal                               = (*HLT_2e15_lhvloose_nod0_L12EM13VH_Nominal      );    
  twoLzeroJ->b2L0J_HLT_2e17_lhvloose_nod0_L12EM15VHI_Nominal                              = (*HLT_2e17_lhvloose_nod0_L12EM15VHI_Nominal       );    
  twoLzeroJ->b2L0J_HLT_2e17_lhvloose_nod0_Nominal                                         = (*HLT_2e17_lhvloose_nod0_Nominal             );    
  twoLzeroJ->b2L0J_HLT_2e24_lhvloose_nod0_Nominal                                         = (*HLT_2e24_lhvloose_nod0_Nominal               );  


  twoLzeroJ->b2L0J_HLT_2g20_tight_icalovloose_L12EM15VHI_Nominal                          = (*HLT_2g20_tight_icalovloose_L12EM15VHI_Nominal         );  
  twoLzeroJ->b2L0J_HLT_2g22_tight_L12EM15VHI_Nominal                                      = (*HLT_2g22_tight_L12EM15VHI_Nominal            );  
  twoLzeroJ->b2L0J_HLT_2g25_loose_g15_loose_Nominal                                       = (*HLT_2g25_loose_g15_loose_Nominal               );
  twoLzeroJ->b2L0J_HLT_2g25_tight_L12EM20VH_Nominal                                       = (*HLT_2g25_tight_L12EM20VH_Nominal               );
  twoLzeroJ->b2L0J_HLT_2g50_loose_L12EM20VH_Nominal                                       = (*HLT_2g50_loose_L12EM20VH_Nominal               );
  twoLzeroJ->b2L0J_HLT_2g20_tight_Nominal                                                 = (*HLT_2g20_tight_Nominal                     );    
  twoLzeroJ->b2L0J_HLT_2g20_loose_g15_loose_Nominal                                       = (*HLT_2g20_loose_g15_loose_Nominal           );    
  twoLzeroJ->b2L0J_HLT_2g22_tight_Nominal                                                 = (*HLT_2g22_tight_Nominal                     );    

  twoLzeroJ->b2L0J_HLT_2mu10_Nominal                                                      = (*HLT_2mu10_Nominal                            );  
  twoLzeroJ->b2L0J_HLT_2mu14_Nominal                                                      = (*HLT_2mu14_Nominal                            );  
  twoLzeroJ->b2L0J_HLT_e120_lhloose_Nominal                                               = (*HLT_e120_lhloose_Nominal                     );  
  twoLzeroJ->b2L0J_HLT_e140_lhloose_nod0_Nominal                                          = (*HLT_e140_lhloose_nod0_Nominal                  );
  twoLzeroJ->b2L0J_HLT_e17_lhloose_2e9_lhloose_Nominal                                    = (*HLT_e17_lhloose_2e9_lhloose_Nominal            );
  twoLzeroJ->b2L0J_HLT_e17_lhloose_mu14_Nominal                                           = (*HLT_e17_lhloose_mu14_Nominal                   );
  twoLzeroJ->b2L0J_HLT_e17_lhloose_nod0_mu14_Nominal                                      = (*HLT_e17_lhloose_nod0_mu14_Nominal          );    
  twoLzeroJ->b2L0J_HLT_e24_lhmedium_L1EM20VH_Nominal                                      = (*HLT_e24_lhmedium_L1EM20VH_Nominal          );    
  twoLzeroJ->b2L0J_HLT_e24_lhmedium_nod0_L1EM20VH_Nominal                                 = (*HLT_e24_lhmedium_nod0_L1EM20VH_Nominal       );    
  twoLzeroJ->b2L0J_HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1_Nominal                        = (*HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1_Nominal         );  
  twoLzeroJ->b2L0J_HLT_e24_lhtight_nod0_ivarloose_Nominal                                 = (*HLT_e24_lhtight_nod0_ivarloose_Nominal         );  
  twoLzeroJ->b2L0J_HLT_e24_lhvloose_nod0_2e12_lhvloose_nod0_L1EM20VH_3EM10VH_Nominal      = (*HLT_e24_lhvloose_nod0_2e12_lhvloose_nod0_L1EM20VH_3EM10VH_Nominal         );  
  twoLzeroJ->b2L0J_HLT_e26_lhmedium_nod0_mu8noL1_Nominal                                  = (*HLT_e26_lhmedium_nod0_mu8noL1_Nominal           );
  twoLzeroJ->b2L0J_HLT_e26_lhtight_nod0_ivarloose_Nominal                                 = (*HLT_e26_lhtight_nod0_ivarloose_Nominal           );
  //twoLzeroJ->b2L0J_HLT_e300_etcut_Nominal                                                 = (*HLT_e300_etcut_Nominal                         );
  twoLzeroJ->b2L0J_HLT_e60_lhmedium_nod0_Nominal                                          = (*HLT_e60_lhmedium_nod0_Nominal              );    
  twoLzeroJ->b2L0J_HLT_e60_lhmedium_Nominal                                               = (*HLT_e60_lhmedium_Nominal                   );    
  twoLzeroJ->b2L0J_HLT_e7_lhmedium_mu24_Nominal                                           = (*HLT_e7_lhmedium_mu24_Nominal               );    
  twoLzeroJ->b2L0J_HLT_e7_lhmedium_nod0_mu24_Nominal                                      = (*HLT_e7_lhmedium_nod0_mu24_Nominal            );  

  twoLzeroJ->b2L0J_HLT_g140_loose_Nominal                                                 = (*HLT_g140_loose_Nominal                       );  
  twoLzeroJ->b2L0J_HLT_g300_etcutHLT_g35_medium_g25_medium_L12EM20VH_Nominal              = (*HLT_g300_etcutHLT_g35_medium_g25_medium_L12EM20VH_Nominal         );  

  twoLzeroJ->b2L0J_HLT_j80_xe80_Nominal                                                   = (*HLT_j80_xe80_Nominal                           );

  twoLzeroJ->b2L0J_HLT_mu18_mu8noL1_Nominal                                               = (*HLT_mu18_mu8noL1_Nominal                       );
  twoLzeroJ->b2L0J_HLT_mu20_iloose_L1MU15_Nominal                                         = (*HLT_mu20_iloose_L1MU15_Nominal                 );
  twoLzeroJ->b2L0J_HLT_mu20_mu8noL1_Nominal                                               = (*HLT_mu20_mu8noL1_Nominal                   );    
  twoLzeroJ->b2L0J_HLT_mu22_mu8noL1_calotag_0eta010_Nominal                               = (*HLT_mu22_mu8noL1_calotag_0eta010_Nominal      );    
  twoLzeroJ->b2L0J_HLT_mu22_mu8noL1_Nominal                                               = (*HLT_mu22_mu8noL1_Nominal                   );    
  twoLzeroJ->b2L0J_HLT_mu24_ivarloose_L1MU15_Nominal                                      = (*HLT_mu24_ivarloose_L1MU15_Nominal            );  
  twoLzeroJ->b2L0J_HLT_mu24_ivarloose_Nominal                                             = (*HLT_mu24_ivarloose_Nominal                   );  
  twoLzeroJ->b2L0J_HLT_mu24_ivarmedium_Nominal                                            = (*HLT_mu24_ivarmedium_Nominal                  );  
  twoLzeroJ->b2L0J_HLT_mu26_ivarmedium_Nominal                                            = (*HLT_mu26_ivarmedium_Nominal                    );
  twoLzeroJ->b2L0J_HLT_mu40_Nominal                                                       = (*HLT_mu40_Nominal                               );
  twoLzeroJ->b2L0J_HLT_mu50_Nominal                                                       = (*HLT_mu50_Nominal                               );
  //twoLzeroJ->b2L0J_HLT_mu60_0eta105_msonly_Nominal                                        = (*HLT_mu60_0eta105_msonly_Nominal            );    


  twoLzeroJ->b2L0J_HLT_xe100_L1XE50_Nominal                                               = (*HLT_xe100_L1XE50_Nominal                   );    
  twoLzeroJ->b2L0J_HLT_xe100_mht_L1XE50_Nominal                                           = (*HLT_xe100_mht_L1XE50_Nominal               );    
  twoLzeroJ->b2L0J_HLT_xe100_pufit_L1XE50_Nominal                                         = (*HLT_xe100_pufit_L1XE50_Nominal             );    
  twoLzeroJ->b2L0J_HLT_xe100_pufit_L1XE55_Nominal                                         = (*HLT_xe100_pufit_L1XE55_Nominal               );  
  twoLzeroJ->b2L0J_HLT_xe100_tc_em_L1XE50_Nominal                                         = (*HLT_xe100_tc_em_L1XE50_Nominal               );  
  twoLzeroJ->b2L0J_HLT_xe110_mht_L1XE50_Nominal                                           = (*HLT_xe110_mht_L1XE50_Nominal                 );  
  twoLzeroJ->b2L0J_HLT_xe110_pueta_L1XE50_Nominal                                         = (*HLT_xe110_pueta_L1XE50_Nominal               );  
  twoLzeroJ->b2L0J_HLT_xe110_pufit_L1XE50_Nominal                                         = (*HLT_xe110_pufit_L1XE50_Nominal               );  
  twoLzeroJ->b2L0J_HLT_xe110_pufit_L1XE55_Nominal                                         = (*HLT_xe110_pufit_L1XE55_Nominal               );  
  twoLzeroJ->b2L0J_HLT_xe110_pufit_xe65_L1XE50_Nominal                                    = (*HLT_xe110_pufit_xe65_L1XE50_Nominal            );
  twoLzeroJ->b2L0J_HLT_xe110_pufit_xe70_L1XE50_Nominal                                    = (*HLT_xe110_pufit_xe70_L1XE50_Nominal            );
  twoLzeroJ->b2L0J_HLT_xe120_pueta_Nominal                                                = (*HLT_xe120_pueta_Nominal                        );
  twoLzeroJ->b2L0J_HLT_xe120_pufit_L1XE50_Nominal                                         = (*HLT_xe120_pufit_L1XE50_Nominal                 );
  twoLzeroJ->b2L0J_HLT_xe120_pufit_Nominal                                                  = (*HLT_xe120_pufit_Nominal                        );
  twoLzeroJ->b2L0J_HLT_xe120_tc_lcw_L1XE50_Nominal                                        = (*HLT_xe120_tc_lcw_L1XE50_Nominal                );
  twoLzeroJ->b2L0J_HLT_xe70_mht_Nominal                                                   = (*HLT_xe70_mht_Nominal                       );    
  twoLzeroJ->b2L0J_HLT_xe70_tc_lcw_Nominal                                                = (*HLT_xe70_tc_lcw_Nominal                    );    
  twoLzeroJ->b2L0J_HLT_xe80_tc_lcw_L1XE50_Nominal                                         = (*HLT_xe80_tc_lcw_L1XE50_Nominal             );    
  twoLzeroJ->b2L0J_HLT_xe90_mht_L1XE50_Nominal                                            = (*HLT_xe90_mht_L1XE50_Nominal                );    
  twoLzeroJ->b2L0J_HLT_xe90_mht_wEFMu_L1XE50_Nominal                                      = (*HLT_xe90_mht_wEFMu_L1XE50_Nominal          );    
  twoLzeroJ->b2L0J_HLT_xe90_pufit_L1XE50_Nominal                                          = (*HLT_xe90_pufit_L1XE50_Nominal              );    
  twoLzeroJ->b2L0J_HLT_xe90_tc_lcw_wEFMu_L1XE50_Nominal                                   = (*HLT_xe90_tc_lcw_wEFMu_L1XE50_Nominal         );

  twoLzeroJ->WriteTree();
}
