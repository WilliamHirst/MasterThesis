#define MySusySkimAnalysis_cxx
// The class definition in MySusySkimAnalysis.h has been generated automatically
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
// root> T->Process("MySusySkimAnalysis.C")
// root> T->Process("MySusySkimAnalysis.C","some options")
// root> T->Process("MySusySkimAnalysis.C+")
//
#define _USE_MATH_DEFINES
#include "./CalcGenericMT2/CalcGenericMT2/MT2_ROOT.h"
#include "MySusySkimAnalysis.h"
#include "HistFitterTree.h"
#include "make2L2JTree.h"
#include "makeMYTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TObjString.h>
#include <time.h>
#include <ctime>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include <numeric>
#include <fstream>
#include <iostream>
//#include "ROOT/RDataFrame.hxx"
//#include "ROOT/RVec.hxx"
//#include "ROOT/RDF/RInterface.hxx"
double gev = 1000.;
double zmass = 91.1876;

clock_t c_begin;
clock_t c_end;
int mypid;

bool doTwoLep = true;
bool doThreeLep = false;

void MySusySkimAnalysis::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
  std::vector< TString >::iterator it1;
  std::vector< TString >::iterator it2;
  std::vector< TString >::iterator it3;
  std::vector< TString >::iterator it4;
  
  std::map< TString, std::vector<TTreeReaderValue<Float_t>> >::iterator mapit1;
  
  cout<<"11"<<endl;
  
  mypid = getpid();

  fOutput = new TSelectorList();

  cutflowstr = "";//CRTOP";//"CRWW";
  DefineCutFlow();

  ev_was_found = 0;

  FNPest = 0.0;
  FNPest_up = 0.0;
  FNPest_dw = 0.0;

  /**
     std::ifstream infile;
     infile.open("eventnumbers_mirto.txt");
     ULong64_t evnum;
     int nev_mirto = 0;
     while (infile >> evnum)
     {
     mirto_ev.push_back(evnum);
     nev_mirto += 1;
     }
     nmirto = 0;
     infile.close();
     cout<<"Loaded "<<nev_mirto<<" events"<<endl;
  */

  //myfile.open("skip_evnums.txt");
  //myfile <<"DSID"<<" "<<"Event"<<" "<<"IFFCLASS"<<" "<<"MotherO"<<" "
  // 	 <<"motherPDGID"<<" "<<"TYPE"<<" "<<"ORIGIN"
  // 	 <<" "<<"MotherT"<<" "<<"px"<<" "<<"py"<<" "<<"pz"<<" "<<"pt"<<" "<<"eta"<<" "<<"passOR"<<" "<<"charge"<<" "<<"truecharge"<<"\n"; 
  
  /**
     TFile::Open("test.root","RECREATE");
     TFile *otfile = new TFile("/scratch3/eirikgr/2L2JInputs/v1.6/SUSY2/SUSY2_Data/data15-16_merged_processed.root");
     TTree *newtree_midl=(TTree*)otfile->Get("data15-16");
     newtree = newtree_midl->CloneTree(0);
     newtree->SetDirectory(0);
     fOutput->Add(newtree);
     otfile->Close();
  */

  //fOutput->Add(newtree);
  doDupCheck = false;
  nDup = 0;

  rnd = new TRandom(time(NULL));

  int nBDT = 20;
  int maxBDT = 1;

  int nPt = 200;
  int maxPt = 1000;

  int nMT2 = 200;
  int maxMT2 = 1000;

  int nMET = 300;
  int maxMET = 1500;

  int nMETSIG = 30;
  int maxMETSIG = 30;

  int nEta = 75;
  int maxEta = 3;

  nev = 0;
  all_nev = 0;

  n_bjet = 0;
  n_sgj = 0;
  n_clj = 0;
  n_clj30 = 0;


  TString option = GetOption();


  TObjArray *tx = option.Tokenize(" ");
  
  is1516 = false;
  is17 = false;
  is18    = false;
  cout<<"111->"<<endl;
  if(tx->GetEntries() >= 2){
    if((((TObjString *)tx->At(1))->String()).EqualTo("year1516")){is1516 = true; years.push_back(2015); years.push_back(2016);} 
    else if((((TObjString *)tx->At(1))->String()).EqualTo("year17")){is17 = true; years.push_back(2017);}
    else if((((TObjString *)tx->At(1))->String()).EqualTo("year18")){is18 = true; years.push_back(2018);}
    else{
      cout<<"ERROR No information about data period!"<<endl;
      return;
    }
  }else{
    cout<<"ERROR No year specified in option - assuming mc16e"<<endl;
    is18 = true;
    years.push_back(2018);
  }
  cout<<"112->"<<endl;
  if(tx->GetEntries() >= 1){
    if((((TObjString *)tx->At(0))->String()).EqualTo("data")){
      isData = true;
      cout<<"This is DATA"<<endl;
    }else{
      isData = false;
      cout<<"This is MC"<<endl;
    }
  }else{
    cout<<"ERROR No type specified in option - assuming mc"<<endl;
  }
  cout<<"113->"<<endl;
  /**
     if(0){
     ROOT::RDataFrame rdf("FNP_WEIGHTS","FNP_270521/FNP_EirikTest_020621.root",{"RunNumber","EventNumber","BDTDeltaM100_90","BDTVVDeltaM100_90","BDTtopDeltaM100_90","BDTothersDeltaM100_90",
     "lepflav1","lepflav2","isSignalLep1","isSignalLep2","lep1pT","lep2pT",
     "njet","nbjet","METsig","isSF","isOS","FNP_TOTAL_UP","FNP_TOTAL_DOWN","FNP_WEIGHTS"});
  
     auto filter = rdf.Filter(is1516 ? ("RunNumber >=  276262 && RunNumber <= 311481") : (is17 ? ("RunNumber >=  325713 && RunNumber <= 340453") : ("RunNumber >= 348885")));
  
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
  */
  cout<<"12"<<endl;
  nentries = 0;
  toprocess = 0;
  firstentry = 0;
  
  if(tx->GetEntries() >= 4){
    toprocess  = (((TObjString *)tx->At(2))->String()).Atoll();
    firstentry = (((TObjString *)tx->At(3))->String()).Atoll();
    cout<<"INFO \t Will process "<<toprocess<<" events, starting from "<<firstentry<<endl;
    nentries = toprocess;
  }

  if(tx->GetEntries() >= 1){
    date_folder = (((TObjString *)tx->At(tx->GetLast()))->String());
  }else{
    cout<<"ERROR \t No folder specfied, creatng dummy"<<endl;
    date_folder = "./dummy";
  }
  delete tx;
  
  mtx = new MatrixMethod();

  uncert["NOMINAL"].push_back(xsec);

  /**
     if(!isData){
     // el eff ID
     uncert["EL_EFF_ID_TOTAL_1down"].push_back(leptonWeight_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR__1down);
     uncert["EL_EFF_ID_TOTAL_1up"].push_back(leptonWeight_EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR__1up);
     // el eff ISO
     uncert["EL_EFF_Iso_TOTAL_1down"].push_back(leptonWeight_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR__1down);
     uncert["EL_EFF_Iso_TOTAL_1up"].push_back(leptonWeight_EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR__1up);
     // el eff RECO
     uncert["EL_EFF_Reco_TOTAL_1down"].push_back(leptonWeight_EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR__1down);
     uncert["EL_EFF_Reco_TOTAL_1up"].push_back(leptonWeight_EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR__1up);

     // mu eff ISO
     uncert["MUON_EFF_ISO_STAT_1down"].push_back(leptonWeight_MUON_EFF_ISO_STAT__1down);
     uncert["MUON_EFF_ISO_STAT_1up"].push_back(leptonWeight_MUON_EFF_ISO_STAT__1up);
     uncert["MUON_EFF_ISO_SYS_1down"].push_back(leptonWeight_MUON_EFF_ISO_SYS__1down);
     uncert["MUON_EFF_ISO_SYS_1up"].push_back(leptonWeight_MUON_EFF_ISO_SYS__1up);

     // mu eff RECO
     uncert["MUON_EFF_RECO_STAT_1down"].push_back(leptonWeight_MUON_EFF_RECO_STAT__1down);
     uncert["MUON_EFF_RECO_STAT_1up"].push_back(leptonWeight_MUON_EFF_RECO_STAT__1up); 
     uncert["MUON_EFF_RECO_STAT_LOWPT_1down"].push_back(leptonWeight_MUON_EFF_RECO_STAT_LOWPT__1down);
     uncert["MUON_EFF_RECO_STAT_LOWPT_1up"].push_back(leptonWeight_MUON_EFF_RECO_STAT_LOWPT__1up);
  
     uncert["MUON_EFF_RECO_SYS_1down"].push_back(leptonWeight_MUON_EFF_RECO_SYS__1down); 
     uncert["MUON_EFF_RECO_SYS_1up"].push_back(leptonWeight_MUON_EFF_RECO_SYS__1up);
     uncert["MUON_EFF_RECO_SYS_LOWPT_1down"].push_back(leptonWeight_MUON_EFF_RECO_SYS_LOWPT__1down);
     uncert["MUON_EFF_RECO_SYS_LOWPT_1up"].push_back(leptonWeight_MUON_EFF_RECO_SYS_LOWPT__1up);

     uncert["MUON_EFF_TTVA_STAT_1down"].push_back(leptonWeight_MUON_EFF_TTVA_STAT__1down);
     uncert["MUON_EFF_TTVA_STAT_1down"].push_back(leptonWeight_MUON_EFF_TTVA_STAT__1up);
     uncert["MUON_EFF_TTVA_STAT_1up"].push_back(leptonWeight_MUON_EFF_TTVA_SYS__1down);
     uncert["MUON_EFF_TTVA_STAT_1down"].push_back(leptonWeight_MUON_EFF_TTVA_SYS__1up);

     // e trig
     //uncert["EL_EFF_TriggerEff_TOTAL_1down"].push_back(trigWeight_EL_EFF_TriggerEff_TOTAL_1NPCOR_PLUS_UNCOR__1down);
     //uncert["EL_EFF_TriggerEff_TOTAL_1up"].push_back(trigWeight_EL_EFF_TriggerEff_TOTAL_1NPCOR_PLUS_UNCOR__1up);
     uncert["EL_EFF_Trigger_TOTAL_1down"].push_back(trigWeight_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1down);
     uncert["EL_EFF_Trigger_TOTAL_1up"].push_back(trigWeight_EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1up);

     // mu trig
     uncert["MUON_EFF_TrigStat_1down"].push_back(trigWeight_MUON_EFF_TrigStatUncertainty__1down);
     uncert["MUON_EFF_TrigStat_1up"].push_back(trigWeight_MUON_EFF_TrigStatUncertainty__1up);
     uncert["MUON_EFF_TrigSyst_1down"].push_back(trigWeight_MUON_EFF_TrigSystUncertainty__1down);
     uncert["MUON_EFF_TrigSyst_1up"].push_back(trigWeight_MUON_EFF_TrigSystUncertainty__1up);

    
     }
  */
  

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

   
  trigval.push_back(trigMatch_1LTrig);
  trigval.push_back(trigMatch_2LTrig);
  trigval.push_back(trigMatch_1L2LTrig);

  triglist.push_back("1LTrig");
  triglist.push_back("2LTrig");
  triglist.push_back("1L2LTrig");

  // Dummy, should not be used.
  trigvallep.push_back(lepHLT_e24_lhmedium_L1EM20VH);
  trigvallep.push_back(lepHLT_e24_lhmedium_L1EM20VH);
  trigvallep.push_back(lepHLT_e24_lhmedium_L1EM20VH);


  trigval.push_back(trigMatch_HLT_e7_lhmedium_mu24);
  trigval.push_back(trigMatch_HLT_e17_lhloose_mu14);
  trigval.push_back(trigMatch_HLT_e26_lhmedium_nod0_mu8noL1);
  triglist.push_back("e7_lhmedium_mu24");
  triglist.push_back("e17_lhloose_mu14");
  triglist.push_back("e26_lhmedium_nod0_mu8noL1");
  trigvallep.push_back(lepHLT_e7_lhmedium_mu24);
  trigvallep.push_back(lepHLT_e17_lhloose_mu14);
  trigvallep.push_back(lepHLT_e26_lhmedium_nod0_mu8noL1);
  
  triglist_m[2015].push_back(trigMatch_HLT_e7_lhmedium_mu24);
  triglist_m[2016].push_back(trigMatch_HLT_e7_lhmedium_mu24);
  triglist_m[2017].push_back(trigMatch_HLT_e7_lhmedium_mu24);
  triglist_m[2018].push_back(trigMatch_HLT_e7_lhmedium_mu24);
  
  triglist_e[2015].push_back(trigMatch_HLT_e7_lhmedium_mu24);
  triglist_e[2016].push_back(trigMatch_HLT_e7_lhmedium_mu24);
  triglist_e[2017].push_back(trigMatch_HLT_e7_lhmedium_mu24);
  triglist_e[2018].push_back(trigMatch_HLT_e7_lhmedium_mu24);

  triglist_m[2015].push_back(trigMatch_HLT_e17_lhloose_mu14);
  triglist_m[2016].push_back(trigMatch_HLT_e17_lhloose_mu14);
  triglist_m[2017].push_back(trigMatch_HLT_e17_lhloose_mu14);
  triglist_m[2018].push_back(trigMatch_HLT_e17_lhloose_mu14);
  
  triglist_e[2015].push_back(trigMatch_HLT_e17_lhloose_mu14);
  triglist_e[2016].push_back(trigMatch_HLT_e17_lhloose_mu14);
  triglist_e[2017].push_back(trigMatch_HLT_e17_lhloose_mu14);
  triglist_e[2018].push_back(trigMatch_HLT_e17_lhloose_mu14);

  triglist_m[2015].push_back(trigMatch_HLT_e26_lhmedium_nod0_mu8noL1);
  triglist_m[2016].push_back(trigMatch_HLT_e26_lhmedium_nod0_mu8noL1);
  triglist_m[2017].push_back(trigMatch_HLT_e26_lhmedium_nod0_mu8noL1);
  triglist_m[2018].push_back(trigMatch_HLT_e26_lhmedium_nod0_mu8noL1);
  
  triglist_e[2015].push_back(trigMatch_HLT_e26_lhmedium_nod0_mu8noL1);
  triglist_e[2016].push_back(trigMatch_HLT_e26_lhmedium_nod0_mu8noL1);
  triglist_e[2017].push_back(trigMatch_HLT_e26_lhmedium_nod0_mu8noL1);
  triglist_e[2018].push_back(trigMatch_HLT_e26_lhmedium_nod0_mu8noL1);

  trigmatch_m[2015].push_back(lepHLT_e7_lhmedium_mu24);
  trigmatch_m[2016].push_back(lepHLT_e7_lhmedium_mu24);
  trigmatch_m[2017].push_back(lepHLT_e7_lhmedium_mu24);
  trigmatch_m[2018].push_back(lepHLT_e7_lhmedium_mu24);
  
  trigmatch_e[2015].push_back(lepHLT_e7_lhmedium_mu24);
  trigmatch_e[2016].push_back(lepHLT_e7_lhmedium_mu24);
  trigmatch_e[2017].push_back(lepHLT_e7_lhmedium_mu24);
  trigmatch_e[2018].push_back(lepHLT_e7_lhmedium_mu24);

  trigmatch_m[2015].push_back(lepHLT_e17_lhloose_mu14);
  trigmatch_m[2016].push_back(lepHLT_e17_lhloose_mu14);
  trigmatch_m[2017].push_back(lepHLT_e17_lhloose_mu14);
  trigmatch_m[2018].push_back(lepHLT_e17_lhloose_mu14);
  
  trigmatch_e[2015].push_back(lepHLT_e17_lhloose_mu14);
  trigmatch_e[2016].push_back(lepHLT_e17_lhloose_mu14);
  trigmatch_e[2017].push_back(lepHLT_e17_lhloose_mu14);
  trigmatch_e[2018].push_back(lepHLT_e17_lhloose_mu14);

  trigmatch_m[2015].push_back(lepHLT_e26_lhmedium_nod0_mu8noL1);
  trigmatch_m[2016].push_back(lepHLT_e26_lhmedium_nod0_mu8noL1);
  trigmatch_m[2017].push_back(lepHLT_e26_lhmedium_nod0_mu8noL1);
  trigmatch_m[2018].push_back(lepHLT_e26_lhmedium_nod0_mu8noL1);
  
  trigmatch_e[2015].push_back(lepHLT_e26_lhmedium_nod0_mu8noL1);
  trigmatch_e[2016].push_back(lepHLT_e26_lhmedium_nod0_mu8noL1);
  trigmatch_e[2017].push_back(lepHLT_e26_lhmedium_nod0_mu8noL1);
  trigmatch_e[2018].push_back(lepHLT_e26_lhmedium_nod0_mu8noL1);

  trigstr_m[2015].push_back("e7_lhmedium_mu24");
  trigstr_m[2016].push_back("e7_lhmedium_mu24");
  trigstr_m[2017].push_back("e7_lhmedium_mu24");
  trigstr_m[2018].push_back("e7_lhmedium_mu24");
  
  trigstr_e[2015].push_back("e7_lhmedium_mu24");
  trigstr_e[2016].push_back("e7_lhmedium_mu24");
  trigstr_e[2017].push_back("e7_lhmedium_mu24");
  trigstr_e[2018].push_back("e7_lhmedium_mu24");

  trigstr_m[2015].push_back("e17_lhloose_mu14");
  trigstr_m[2016].push_back("e17_lhloose_mu14");
  trigstr_m[2017].push_back("e17_lhloose_mu14");
  trigstr_m[2018].push_back("e17_lhloose_mu14");
  
  trigstr_e[2015].push_back("e17_lhloose_mu14");
  trigstr_e[2016].push_back("e17_lhloose_mu14");
  trigstr_e[2017].push_back("e17_lhloose_mu14");
  trigstr_e[2018].push_back("e17_lhloose_mu14");

  trigstr_m[2015].push_back("e26_lhmedium_nod0_mu8noL1");
  trigstr_m[2016].push_back("e26_lhmedium_nod0_mu8noL1");
  trigstr_m[2017].push_back("e26_lhmedium_nod0_mu8noL1");
  trigstr_m[2018].push_back("e26_lhmedium_nod0_mu8noL1");
  
  trigstr_e[2015].push_back("e26_lhmedium_nod0_mu8noL1");
  trigstr_e[2016].push_back("e26_lhmedium_nod0_mu8noL1");
  trigstr_e[2017].push_back("e26_lhmedium_nod0_mu8noL1");
  trigstr_e[2018].push_back("e26_lhmedium_nod0_mu8noL1");
  
  // 2015
  trigval.push_back(trigMatch_HLT_mu26_imedium);
  triglist.push_back("mu26_imedium");
  trigvallep.push_back(lepHLT_mu26_imedium);
  
  triglist_m[2015].push_back(trigMatch_HLT_mu26_imedium);
  trigmatch_m[2015].push_back(lepHLT_mu26_imedium);
  trigstr_m[2015].push_back("mu26_imedium");
  

  // 2015-->
  trigval.push_back(trigMatch_HLT_mu50);
  triglist.push_back("mu50");
  trigvallep.push_back(lepHLT_mu50);
  
  triglist_m[2015].push_back(trigMatch_HLT_mu50);
  triglist_m[2016].push_back(trigMatch_HLT_mu50);
  triglist_m[2017].push_back(trigMatch_HLT_mu50);
  triglist_m[2018].push_back(trigMatch_HLT_mu50);

  trigmatch_m[2015].push_back(lepHLT_mu50);
  trigmatch_m[2016].push_back(lepHLT_mu50);
  trigmatch_m[2017].push_back(lepHLT_mu50);
  trigmatch_m[2018].push_back(lepHLT_mu50);

  
  trigstr_m[2015].push_back("mu50");
  trigstr_m[2016].push_back("mu50");
  trigstr_m[2017].push_back("mu50");
  trigstr_m[2018].push_back("mu50");
 

  // 2016-->
  trigval.push_back(trigMatch_HLT_mu26_ivarmedium);
  triglist.push_back("mu26_ivarmedium");
  trigvallep.push_back(lepHLT_mu26_ivarmedium);
   
  triglist_m[2016].push_back(trigMatch_HLT_mu26_ivarmedium);
  triglist_m[2017].push_back(trigMatch_HLT_mu26_ivarmedium);
  triglist_m[2018].push_back(trigMatch_HLT_mu26_ivarmedium);

  trigmatch_m[2016].push_back(lepHLT_mu26_ivarmedium);
  trigmatch_m[2017].push_back(lepHLT_mu26_ivarmedium);
  trigmatch_m[2018].push_back(lepHLT_mu26_ivarmedium);
  
  trigstr_m[2016].push_back("mu26_ivarmedium");
  trigstr_m[2017].push_back("mu26_ivarmedium");
  trigstr_m[2018].push_back("mu26_ivarmedium");

  // 2015
  trigval.push_back(trigMatch_HLT_2e12_lhloose_L12EM10VH);
  triglist.push_back("2e12_lhloose_L12EM10VH");
  trigvallep.push_back(lepHLT_2e12_lhloose_L12EM10VH);

  triglist_e[2015].push_back(trigMatch_HLT_2e12_lhloose_L12EM10VH);
  trigmatch_e[2015].push_back(lepHLT_2e12_lhloose_L12EM10VH);
  trigstr_e[2015].push_back("2e12_lhloose_L12EM10VH");

  // 2016
  trigval.push_back(trigMatch_HLT_2e17_lhvloose_nod0);
  triglist.push_back("2e17_lhvloose_nod0");
  trigvallep.push_back(lepHLT_2e17_lhvloose_nod0);
  
  triglist_e[2016].push_back(trigMatch_HLT_2e17_lhvloose_nod0);
  trigmatch_e[2016].push_back(lepHLT_2e17_lhvloose_nod0);
  trigstr_e[2016].push_back("2e17_lhvloose_nod0");

  //2017 (325713–326833) and (328394–340453) and 2018
  trigval.push_back(trigMatch_HLT_2e24_lhvloose_nod0);
  triglist.push_back("2e24_lhvloose_nod0");
  trigvallep.push_back(lepHLT_2e24_lhvloose_nod0);
  
  triglist_e[2017].push_back(trigMatch_HLT_2e24_lhvloose_nod0);
  triglist_e[2018].push_back(trigMatch_HLT_2e24_lhvloose_nod0);
  trigmatch_e[2017].push_back(lepHLT_2e24_lhvloose_nod0);
  trigmatch_e[2018].push_back(lepHLT_2e24_lhvloose_nod0);
  trigstr_e[2017].push_back("2e24_lhvloose_nod0");
  trigstr_e[2018].push_back("2e24_lhvloose_nod0");
 
  trigval.push_back(trigMatch_HLT_2e17_lhvloose_nod0_L12EM15VHI);
  triglist.push_back("2e17_lhvloose_nod0_L12EM15VHI");
  trigvallep.push_back(lepHLT_2e17_lhvloose_nod0_L12EM15VHI);
  
  triglist_e[2017].push_back(trigMatch_HLT_2e17_lhvloose_nod0_L12EM15VHI);
  triglist_e[2018].push_back(trigMatch_HLT_2e17_lhvloose_nod0_L12EM15VHI);
  trigmatch_e[2017].push_back(lepHLT_2e17_lhvloose_nod0_L12EM15VHI);
  trigmatch_e[2018].push_back(lepHLT_2e17_lhvloose_nod0_L12EM15VHI);
  trigstr_e[2017].push_back("2e17_lhvloose_nod0_L12EM15VHI");
  trigstr_e[2018].push_back("2e17_lhvloose_nod0_L12EM15VHI");
  
  // 2017 (326834–328393)
  // trigval.push_back(trigMatch_HLT_2e24_lhvloose_nod0);
  // triglist.push_back("2e24_lhvloose_nod0");
  // trigvallep.push_back(lepHLT_2e24_lhvloose_nod0);
  // triglist_e[2017].push_back(trigMatch_HLT_2e24_lhvloose_nod0);


  trigcat["mu26_ivarmedium"] = "ivarmed";
  trigcat["mu26_imedium"] = "imed";
  trigcat["mu50"] = "no";

  trigcat["2e12_lhloose_L12EM10VH"] = "lhlooseL1";
  trigcat["2e17_lhvloose_nod0"] = "lhvloose";
  trigcat["2e24_lhvloose_nod0"] = "lhvloose";
  trigcat["2e17_lhvloose_nod0_L12EM15VHI"] = "lhvlooseL1";

  trigcat["e7_lhmedium_mu24"] = "lhmed";	  
  trigcat["e17_lhloose_mu14"] = "lhlooseL1";
  trigcat["e26_lhmedium_nod0_mu8noL1"] = "lhmed";
  
  trigcategories["mu26_ivarmedium"] = "1M";		  
  trigcategories["mu26_imedium"] = "1M";			  
  trigcategories["mu50"] = "1M";				  
  
  trigcategories["2e12_lhloose_L12EM10VH"] = "2E";	  
  trigcategories["2e17_lhvloose_nod0"] = "2E";		  
  trigcategories["2e24_lhvloose_nod0"] = "2E";		  
  trigcategories["2e17_lhvloose_nod0_L12EM15VHI"] = "2E";

  trigcategories["e7_lhmedium_mu24"] = "1E1M";	  
  trigcategories["e17_lhloose_mu14"] = "1E1M";
  trigcategories["e26_lhmedium_nod0_mu8noL1"] = "1E1M";
  

  trig_order.push_back("no");
  trig_order.push_back("lhlooseL1");
  trig_order.push_back("lhvloose");
  trig_order.push_back("lhvlooseL1");
  trig_order.push_back("imed");
  trig_order.push_back("ivarmed");

  categories = {"e_lhlooseL1","e_lhvloose","e_lhvlooseL1","m_imed","e_lhmed","m_ivarmed","m_no"};

  triggerDef.push_back("notrigm");
  triggerDef.push_back("no");
  triggerDef.push_back("ivarmed");
  triggerDef.push_back("imed");
  triggerDef.push_back("lhvlooseL1");
  triggerDef.push_back("lhvloose");
  triggerDef.push_back("lhlooseL1");
   
  
  // 1L analysis
  /**
     trigval.push_back(trigMatch_HLT_e24_lhmedium_L1EM20VH);
     //trigval.push_back(trigMatch_HLT_e24_lhmedium_nod0_L1EM20VH);  
     trigval.push_back(trigMatch_HLT_mu20_iloose_L1MU15);    
     trigval.push_back(trigMatch_HLT_e60_lhmedium);              
     //  trigval.push_back(trigMatch_HLT_mu60_0eta105_msonly);    
     trigval.push_back(trigMatch_HLT_e120_lhloose);              

     trigval.push_back(trigMatch_HLT_e24_lhtight_nod0_ivarloose);   
     trigval.push_back(trigMatch_HLT_e26_lhtight_nod0_ivarloose);   
     trigval.push_back(trigMatch_HLT_e60_lhmedium_nod0);             
     trigval.push_back(trigMatch_HLT_mu26_ivarmedium);         
     trigval.push_back(trigMatch_HLT_e140_lhloose_nod0);             
     trigval.push_back(trigMatch_HLT_mu50);         
     //  trigval.push_back(trigMatch_HLT_e300_etcut);

     //1E2015 {"HLT_e24_lhmedium_L1EM20VH","HLT_e60_lhmedium","HLT_e120_lhloose"} )
     // 1M2015 {"HLT_mu20_iloose_L1MU15","HLT_mu50","HLT_mu60_0eta105_msonly"} )
     // 1E2016 {"HLT_e24_lhtight_nod0_ivarloose","HLT_e24_lhmedium_nod0_L1EM20VH" ,"HLT_e26_lhtight_nod0_ivarloose","HLT_e60_lhmedium_nod0","HLT_e140_lhloose_nod0","HLT_e300_etcut"}
     // 1M2016 {"HLT_mu26_ivarmedium","HLT_mu50","HLT_mu60_0eta105_msonly"}
     // 1E2017 {"HLT_e26_lhtight_nod0_ivarloose","HLT_e60_lhmedium_nod0","HLT_e140_lhloose_nod0","HLT_e300_etcut"}
     // 1M2017 {"HLT_mu26_ivarmedium","HLT_mu50","HLT_mu60_0eta105_msonly"} 

     triglist.push_back("e24_lhmedium_L1EM20VH"); // ok
     triglist.push_back("mu20_iloose_L1MU15");     //ok
     triglist.push_back("e60_lhmedium");              // ok
     triglist.push_back("e120_lhloose");              //ok

     triglist.push_back("e24_lhtight_nod0_ivarloose");  // ok
     triglist.push_back("e26_lhtight_nod0_ivarloose");   
     triglist.push_back("e60_lhmedium_nod0");             
     triglist.push_back("mu26_ivarmedium");         
     triglist.push_back("e140_lhloose_nod0");             
     triglist.push_back("mu50");         

     // 2015
     triglist_e[2015].push_back(trigMatch_HLT_e24_lhmedium_L1EM20VH);
     triglist_e[2015].push_back(trigMatch_HLT_e60_lhmedium);
     triglist_e[2015].push_back(trigMatch_HLT_e120_lhloose);
  
     triglist_m[2015].push_back(trigMatch_HLT_mu20_iloose_L1MU15);
     //triglist_m[2015].push_back(trigMatch_HLT_mu60_0eta105_msonly);
     triglist_m[2015].push_back(trigMatch_HLT_mu50);

     // 2016
     triglist_e[2016].push_back(trigMatch_HLT_e24_lhmedium_nod0_L1EM20VH);
     triglist_e[2016].push_back(trigMatch_HLT_e24_lhtight_nod0_ivarloose);
     triglist_e[2016].push_back(trigMatch_HLT_e26_lhtight_nod0_ivarloose);
     triglist_e[2016].push_back(trigMatch_HLT_e60_lhmedium_nod0);
     triglist_e[2016].push_back(trigMatch_HLT_e140_lhloose_nod0);
     // triglist_e[2016].push_back(trigMatch_HLT_e300_etcut);
    
     // triglist_m[2016].push_back(trigMatch_HLT_mu60_0eta105_msonly);
     triglist_m[2016].push_back(trigMatch_HLT_mu26_ivarmedium);
     triglist_m[2016].push_back(trigMatch_HLT_mu50);

     // 2017
     triglist_e[2017].push_back(trigMatch_HLT_e26_lhtight_nod0_ivarloose);
     triglist_e[2017].push_back(trigMatch_HLT_e60_lhmedium_nod0);
     triglist_e[2017].push_back(trigMatch_HLT_e140_lhloose_nod0);
     //triglist_e[2017].push_back(trigMatch_HLT_e300_etcut);
  
     triglist_m[2017].push_back(trigMatch_HLT_mu26_ivarmedium);
     //triglist_m[2017].push_back(trigMatch_HLT_mu60_0eta105_msonly);
     triglist_m[2017].push_back(trigMatch_HLT_mu50);                    

     // 2018
     triglist_e[2018].push_back(trigMatch_HLT_e26_lhtight_nod0_ivarloose);
     triglist_e[2018].push_back(trigMatch_HLT_e60_lhmedium_nod0);
     triglist_e[2018].push_back(trigMatch_HLT_e140_lhloose_nod0);
     //triglist_e[2018].push_back(trigMatch_HLT_e300_etcut);
  
     triglist_m[2018].push_back(trigMatch_HLT_mu26_ivarmedium);
     //triglist_m[2018].push_back(trigMatch_HLT_mu60_0eta105_msonly);
     triglist_m[2018].push_back(trigMatch_HLT_mu50);

     // TriggerMatch
     // 2015
     trigmatch_e[2015].push_back(lepHLT_e24_lhmedium_L1EM20VH);
     trigmatch_e[2015].push_back(lepHLT_e60_lhmedium);
     trigmatch_e[2015].push_back(lepHLT_e120_lhloose);
  
     trigmatch_m[2015].push_back(lepHLT_mu20_iloose_L1MU15);
     //trigmatch_m[2015].push_back(lepHLT_mu60_0eta105_msonly);
     trigmatch_m[2015].push_back(lepHLT_mu50);

     // 2016
     trigmatch_e[2016].push_back(lepHLT_e24_lhmedium_nod0_L1EM20VH);
     trigmatch_e[2016].push_back(lepHLT_e24_lhtight_nod0_ivarloose);
     trigmatch_e[2016].push_back(lepHLT_e26_lhtight_nod0_ivarloose);
     trigmatch_e[2016].push_back(lepHLT_e60_lhmedium_nod0);
     trigmatch_e[2016].push_back(lepHLT_e140_lhloose_nod0);
     //trigmatch_e[2016].push_back(lepHLT_e300_etcut);
    
     // trigmatch_m[2016].push_back(lepHLT_mu60_0eta105_msonly);
     trigmatch_m[2016].push_back(lepHLT_mu26_ivarmedium);
     trigmatch_m[2016].push_back(lepHLT_mu50);

     // 2017
     trigmatch_e[2017].push_back(lepHLT_e26_lhtight_nod0_ivarloose);
     trigmatch_e[2017].push_back(lepHLT_e60_lhmedium_nod0);
     trigmatch_e[2017].push_back(lepHLT_e140_lhloose_nod0);
     //trigmatch_e[2017].push_back(lepHLT_e300_etcut);
  
     trigmatch_m[2017].push_back(lepHLT_mu26_ivarmedium);
     //trigmatch_m[2017].push_back(lepHLT_mu60_0eta105_msonly);
     trigmatch_m[2017].push_back(lepHLT_mu50);

     // 2018
     trigmatch_e[2018].push_back(lepHLT_e26_lhtight_nod0_ivarloose);
     trigmatch_e[2018].push_back(lepHLT_e60_lhmedium_nod0);
     trigmatch_e[2018].push_back(lepHLT_e140_lhloose_nod0);
     //trigmatch_e[2018].push_back(lepHLT_e300_etcut);
  
     trigmatch_m[2018].push_back(lepHLT_mu26_ivarmedium);
     //  trigmatch_m[2018].push_back(lepHLT_mu60_0eta105_msonly);
     trigmatch_m[2018].push_back(lepHLT_mu50);

     // TriggerName (IMPORTANT: must have same order as trigmatch (no so important with triglist)!!!)
     // 2015
     trigstr_e[2015].push_back("e24_lhmedium_L1EM20VH");
     trigstr_e[2015].push_back("e60_lhmedium");
     trigstr_e[2015].push_back("e120_lhloose");
  
     trigstr_m[2015].push_back("mu20_iloose_L1MU15");
     //  trigstr_m[2015].push_back("HLT_mu60_0eta105_msonly");
     trigstr_m[2015].push_back("mu50");

     // 2016
     trigstr_e[2016].push_back("e24_lhmedium_nod0_L1EM20VH");
     trigstr_e[2016].push_back("e24_lhtight_nod0_ivarloose");
     trigstr_e[2016].push_back("e26_lhtight_nod0_ivarloose");
     trigstr_e[2016].push_back("e60_lhmedium_nod0");
     trigstr_e[2016].push_back("e140_lhloose_nod0");
     // trigstr_e[2016].push_back("e300_etcut");
    
     //  trigstr_m[2016].push_back("mu60_0eta105_msonly");
     trigstr_m[2016].push_back("mu26_ivarmedium");
     trigstr_m[2016].push_back("mu50");

     // 2017
     trigstr_e[2017].push_back("e26_lhtight_nod0_ivarloose");
     trigstr_e[2017].push_back("e60_lhmedium_nod0");
     trigstr_e[2017].push_back("e140_lhloose_nod0");
     //trigstr_e[2017].push_back("e300_etcut");
  
     trigstr_m[2017].push_back("mu26_ivarmedium");
     //  trigstr_m[2017].push_back("mu60_0eta105_msonly");
     trigstr_m[2017].push_back("mu50");

     // 2018
     trigstr_e[2018].push_back("e26_lhtight_nod0_ivarloose");
     trigstr_e[2018].push_back("e60_lhmedium_nod0");
     trigstr_e[2018].push_back("e140_lhloose_nod0");
     // trigstr_e[2018].push_back("e300_etcut");
  
     trigstr_m[2018].push_back("mu26_ivarmedium");
     //  trigstr_m[2018].push_back("mu60_0eta105_msonly");
     trigstr_m[2018].push_back("mu50");


  
     // Dummy, should not be used.
     trigvallep.push_back(lepHLT_e24_lhmedium_L1EM20VH);
     trigvallep.push_back(lepHLT_e24_lhmedium_L1EM20VH);
     trigvallep.push_back(lepHLT_e24_lhmedium_L1EM20VH);

     trigvallep.push_back(lepHLT_e24_lhmedium_L1EM20VH);
     //trigvallep.push_back(lepHLT_e24_lhmedium_nod0_L1EM20VH);  
     trigvallep.push_back(lepHLT_mu20_iloose_L1MU15);    
     trigvallep.push_back(lepHLT_e60_lhmedium);              
     //  trigvallep.push_back(lepHLT_mu60_0eta105_msonly);    
     trigvallep.push_back(lepHLT_e120_lhloose);              

     trigvallep.push_back(lepHLT_e24_lhtight_nod0_ivarloose);   
     trigvallep.push_back(lepHLT_e26_lhtight_nod0_ivarloose);   
     trigvallep.push_back(lepHLT_e60_lhmedium_nod0);             
     trigvallep.push_back(lepHLT_mu26_ivarmedium);         
     trigvallep.push_back(lepHLT_e140_lhloose_nod0);             
     trigvallep.push_back(lepHLT_mu50);         
     //  trigvallep.push_back(lepHLT_e300_etcut);


     // TriggerCategory
     // 2015
     trigcat["e24_lhmedium_L1EM20VH"] = "med";
     trigcat["e60_lhmedium"] = "med";
     trigcat["e120_lhloose"] = "loose";
  
     trigcat["mu20_iloose_L1MU15"] = "loose";
     //  trigcat["mu60_0eta105_msonly"] = "no";
     trigcat["mu50"] = "no";

     // 2016
     trigcat["e24_lhmedium_nod0_L1EM20VH"] = "med";
     trigcat["e24_lhtight_nod0_ivarloose"] = "tight_iso";
     trigcat["e26_lhtight_nod0_ivarloose"] = "tight_iso";
     trigcat["e60_lhmedium_nod0"] = "med";
     trigcat["e140_lhloose_nod0"] = "loose";
     //trigcat["e300_etcut"] = "no";
    
     //  trigcat["mu60_0eta105_msonly"] = "no";
     trigcat["mu26_ivarmedium"] = "med";
     trigcat["mu50"] = "no";

     // 2017
     trigcat["e26_lhtight_nod0_ivarloose"] = "tight_iso";
     trigcat["e60_lhmedium_nod0"] = "med";
     trigcat["e140_lhloose_nod0"] = "loose";
     //trigcat["e300_etcut"] = "no";
  
     trigcat["mu26_ivarmedium"] = "med";
     //  trigcat["mu60_0eta105_msonly"] = "no";
     trigcat["mu50"] = "no";

     // 2018
     trigcat["e26_lhtight_nod0_ivarloose"] = "tight_iso";
     trigcat["e60_lhmedium_nod0"] = "med";
     trigcat["e140_lhloose_nod0"] = "loose";
     //trigcat["e300_etcut"] = "no";
  
     trigcat["mu26_ivarmedium"] = "med";
     //  trigcat["mu60_0eta105_msonly"] = "no";
     trigcat["mu50"] = "no";
    

     trigcategories["e24_lhmedium_L1EM20VH"] = "1E";
     trigcategories["mu20_iloose_L1MU15"] = "1M";    
     trigcategories["e60_lhmedium"] = "1E";              
     trigcategories["e120_lhloose"] = "1E";              

     trigcategories["e24_lhtight_nod0_ivarloose"] = "1E";   
     trigcategories["e26_lhtight_nod0_ivarloose"] = "1E";   
     trigcategories["e60_lhmedium_nod0"] = "1E";             
     trigcategories["mu26_ivarmedium"] = "1M";         
     trigcategories["e140_lhloose_nod0"] = "1E";             
     trigcategories["mu50"] = "1M";         
     //  trigcategories["e300_etcut"] = "1E";

     trig_order.push_back("tight_iso");
     trig_order.push_back("med");
     trig_order.push_back("loose");
     trig_order.push_back("no");
 
  */
  
  /**
     trigcategories["e24_lhmedium_L1EM20VH"] = "e_med";
     //trigcategories["e24_lhmedium_nod0_L1EM20VH"] = "e_med";
     trigcategories["e60_lhmedium_nod0"] = "e_med";
     trigcategories["e60_lhmedium"] = "e_med";
  
     trigcategories["e120_lhloose"] = "e_loose";
     trigcategories["e140_lhloose_nod0"] = "e_loose";
  
     trigcategories["e24_lhtight_nod0_ivarloose"] = "e_tight_iso";
     trigcategories["e26_lhtight_nod0_ivarloose"] = "e_tight_iso";
  
     trigcategories["mu26_ivarmedium"] = "m_med";
  
     trigcategories["mu20_iloose_L1MU15"] = "m_loose";
  
     trigcategories["mu50"] = "m_no";
  */

  /**
     categories = {"e_med","e_loose","e_tight_iso","m_med","m_loose","m_no"};

     triggerDef.push_back("notrigm");
     triggerDef.push_back("med");
     triggerDef.push_back("loose");
     triggerDef.push_back("tight_iso");
     triggerDef.push_back("no");

  */
  std::vector<TString> label_bins_e;
  std::vector<TString> label_bins_m;
  for(std::vector< int >::iterator it_int = years.begin(); it_int != years.end(); it_int++){
    cout<<"yr = "<<(*it_int)<<endl;
    cout<<"trigstr_e[yr].size()+1 = "<<trigstr_e[(*it_int)].size()+1<<endl;
    cout<<"trigstr_m[yr].size()+1 = "<<trigstr_m[(*it_int)].size()+1<<endl;
    for(unsigned int i = 0; i < trigstr_e[(*it_int)].size(); i++){
      label_bins_e.push_back(trigstr_e[(*it_int)].at(i));
    }
    for(unsigned int i = 0; i < trigstr_m[(*it_int)].size(); i++){
      label_bins_m.push_back(trigstr_m[(*it_int)].at(i));
    }
  }

  h_triggermatched_el = new TH1F("h_triggermatched_el","",label_bins_m.size()+1,0,label_bins_m.size()+1);
  h_triggermatched_mu = new TH1F("h_triggermatched_mu","",label_bins_m.size()+1,0,label_bins_m.size()+1);
  fOutput->AddLast(h_triggermatched_el);
  fOutput->AddLast(h_triggermatched_mu);
  for(unsigned int i = 0; i < label_bins_e.size(); i++){h_triggermatched_el->GetXaxis()->SetBinLabel(i+2,label_bins_e.at(i));}
  for(unsigned int i = 0; i < label_bins_m.size(); i++){h_triggermatched_mu->GetXaxis()->SetBinLabel(i+2,label_bins_m.at(i));}

  h_min_dphi_lep_met = new TH1F("h_min_dphi_lep_met","",160,-8,8);
  fOutput->AddLast(h_min_dphi_lep_met);

  vector<TString> MM_regions;
  vector<TString> controlregions;
  vector<TString> realregions;
  vector<TString> fakeregions;
  vector<TString> lepregions;
  vector<TString> ptregions;
  vector<TString> ncljregions;
  vector<TString> nbjregions;
  //vector<TString> metsigregions;
  vector<TString> metregions;
  vector<TString> mt2regions;
  vector<TString> mllregions;
  vector<TString> trigregions;

  cout<<"14"<<endl;
  doRJR          = false;

  if(!doRJR){
    // looseTightDef.push_back("HP_MT_FCL_FCHPT");
    // looseTightDef.push_back("HP_0M_FCL_FCHPT");
    // looseTightDef.push_back("HP_MT_FCT");
    // looseTightDef.push_back("HP_0M_FCT");
    //looseTightDef.push_back("ODA");
    //looseTightDef.push_back("2L2J");
    //looseTightDef.push_back("FCL_HPTC");
    //looseTightDef.push_back("2L2J");
    //looseTightDef.push_back("FCT");
    //looseTightDef.push_back("FCT");
    //looseTightDef.push_back("FCT_MT");
    //looseTightDef.push_back("ODA");
    //looseTightDef.push_back("FCL");
    looseTightDef.push_back("ZMET");
    //looseTightDef.push_back("GRAD");
    // looseTightDef.push_back("GRAD");
    // looseTightDef.push_back("FCL_FCTT");
    // looseTightDef.push_back("FCL_FCHPT");
  }else{
    //looseTightDef.push_back("2L2J");
    looseTightDef.push_back("FCT");
    //looseTightDef.push_back("FCT_MT");
    doThreeLep = true;
  }
  

  if(false){ // doFakeEstimation
    allEventsPass  = false;
    doCoreRegions  = false;
    doAnalysis     = true;
    doFakeEst      = true;
    makeHistograms = true;
    makeNtup       = true;
  }else{ // Get fake rates/real eff.
    allEventsPass = false;
    doCoreRegions  = false;
    doAnalysis     = false;
    doFakeEst      = false;
    makeHistograms = true;
    makeNtup       = true;
  }
  
  if(!isData)doFakeEst      = false;
  if(!makeNtup){
    make2L2JNtup   = false;
    makeMYNtup   = false;
  }else{
    make2L2JNtup   = false;
    makeMYNtup   = false;
  }
  if(doThreeLep){
    //computeRates   = true;
    fullScan       = false;
  }else{
    //computeRates   = true;
    fullScan       = false;
  }

  // Remove triggers from list if 2L-triggers are used!
  if(looseTightDef.size() == 1 and (looseTightDef.at(0).EqualTo("2L2J"))){// or looseTightDef.at(0).EqualTo("ZMET"))){
    triglist.clear();
    triggerDef.clear();
    triggerDef.push_back(looseTightDef.at(0));
  }
  if(looseTightDef.size() > 1 and !doAnalysis)doCoreRegions  = true;
  if(!makeHistograms)return;
  
  if(!doRJR){
    if(doThreeLep)num_of_leptons = 3;
    else num_of_leptons = 2;
  }else{
    num_of_leptons = 23;
    doThreeLep = true;
  }

  // November 2019: done in python script retrieved through options
  /**
     time_t now = time(0);
     char* dt = ctime(&now);
     string dt_str = string(dt);
     vector<string> tokens; // Create vector to hold our words
     string buf;
     stringstream ss(dt_str); // Insert the string into a stream
     while (ss >> buf)
     tokens.push_back(buf);
     //sprintf(date_folder,"/scratch2/eirikgr/%s_%s_%s_%s_%02dL_NTUP_FCT_%s",tokens.at(0).c_str(),tokens.at(1).c_str(),tokens.at(2).c_str(),tokens.at(tokens.size()-1).c_str(),num_of_leptons,
     //is1516 ? "data1516" : (is17 ? "data17" : "data18"));
     sprintf(date_folder,"/scratch/eirikgr/%s_%s_%s_%s_%02dL_NTUP_FCT_%s",tokens.at(0).c_str(),tokens.at(1).c_str(),tokens.at(2).c_str(),tokens.at(tokens.size()-1).c_str(),num_of_leptons,
     is1516 ? "data1516" : (is17 ? "data17" : "data18"));
     //sprintf(date_folder,"./%s_%s_%s_%s_%iL_NTUP",tokens.at(0).c_str(),tokens.at(1).c_str(),tokens.at(2).c_str(),tokens.at(tokens.size()-1).c_str(),num_of_leptons);
     struct stat st;
     //if(stat(date_folder,&st) == 0){
     if(!((stat(date_folder, &st) == 0) && (((st.st_mode) & S_IFMT) == S_IFDIR))) {
     //printf("Folder %s is present\n",date_folder);
     printf("Creating directory %s\n",date_folder);
     mkdir(date_folder,0750);
     }
  */
  
  

  //trigregions.push_back("1LTrig");
 
  trigregions.push_back("2LTrig");
  /**
     trigregions.push_back("1L2LTrig");
     trigregions.push_back("1LTrig");
  */
  trigregions.push_back("");

  controlregions.push_back("ALL");
  if(doThreeLep){
    controlregions.push_back("EEE");
    controlregions.push_back("MMM");
    controlregions.push_back("EEM");
    controlregions.push_back("MME");
  }
  if(!doThreeLep || doRJR){
    // remember to uncomment in getFillString
    controlregions.push_back("EEOS");
    //TESTING -->
   
    controlregions.push_back("EESS");
    controlregions.push_back("EMOS");
    controlregions.push_back("EMSS");
    controlregions.push_back("MMOS");
    //TESTING -->//controlregions.push_back("MMSS");
    // */

    // remember to uncomment line in getFillVector
    //controlregions.push_back("SFOS");
    //controlregions.push_back("SFSS");
    // controlregions.push_back("EE");
    // controlregions.push_back("EM");
    // controlregions.push_back("MM");
  }

  // if(computeRates){
  realregions.push_back("2L01");
  realregions.push_back("2L02");
  realregions.push_back("2L12");
  realregions.push_back("2L13");
  //if(doThreeLep)
  realregions.push_back("2L03");
  realregions.push_back("2L04");
  //realregions.push_back("2L05");
  //realregions.push_back("2L06");
  if(!doAnalysis){

    
    // TESTING -->
    realregions.push_back("2L05");
    realregions.push_back("2L06");
    realregions.push_back("2L07");
    realregions.push_back("2L08");
    realregions.push_back("2L09");
    realregions.push_back("2L10");
    realregions.push_back("2L11");
    // }
    
    
    if(!doCoreRegions){
      fakeregions.push_back("2L01");
      fakeregions.push_back("2L02");
      fakeregions.push_back("2L04");
      /**
      fakeregions.push_back("2L20");
      fakeregions.push_back("2L30");
      */
    }
    
    //fakeregions.push_back("2L05");
    
    
    // if(doThreeLep){
    //   fakeregions.push_back("2L03");
    //   fakeregions.push_back("2L05");
    //   fakeregions.push_back("2L06");
    // }
    fakeregions.push_back("2L21");
    

    /**
     // Removed for quick running
     fakeregions.push_back("2L49");
     fakeregions.push_back("2L40");
     fakeregions.push_back("2L41");
     fakeregions.push_back("2L42");
     fakeregions.push_back("2L43");
     fakeregions.push_back("2L44");
     fakeregions.push_back("2L45");
    */
    // fakeregions.push_back("2L31");
    // fakeregions.push_back("2L32");
    // fakeregions.push_back("2L33");
    //    fakeregions.push_back("2L23");

    /**
     // Removed for quick running
     fakeregions.push_back("2L24");
     fakeregions.push_back("2L25");
     fakeregions.push_back("2L26");
     fakeregions.push_back("2L27");
    */

    
    lepregions.push_back("E");
    lepregions.push_back("M");
  }

  if(fullScan){
    metregions.push_back("met_min_50");
    metregions.push_back("met_50_100");
    metregions.push_back("met_100_150");
    metregions.push_back("met_150_200");
    metregions.push_back("met_200_max");

    mt2regions.push_back("mt2_min_50");
    mt2regions.push_back("mt2_50_100");
    mt2regions.push_back("mt2_100_150");
    mt2regions.push_back("mt2_150_200");
    mt2regions.push_back("mt2_200_max");
    /**
       metsigregions.push_back("metsig_min_10");
       metsigregions.push_back("metsig_10_20");
       metsigregions.push_back("metsig_20_30");
       metsigregions.push_back("metsig_30_50");
       metsigregions.push_back("metsig_50_max");
    */
    mllregions.push_back("mll_min_50");
    mllregions.push_back("mll_50_80");
    mllregions.push_back("mll_80_110");
    mllregions.push_back("mll_110_200");
    mllregions.push_back("mll_200_300");
    mllregions.push_back("mll_300_max");

    ncljregions.push_back("nclj_0_1");
    ncljregions.push_back("nclj_1_2");
    ncljregions.push_back("nclj_2_max");

    nbjregions.push_back("nbj_0_1");
    nbjregions.push_back("nbj_1_max");

    ptregions.push_back("pT_min_30");
    ptregions.push_back("pT_30_40");
    ptregions.push_back("pT_40_50");
    ptregions.push_back("pT_50_60");
    ptregions.push_back("pT_60_80");
    ptregions.push_back("pT_80_max");
  }

  if(!doRJR/** && !doCoreRegions*/){

    for(std::vector< TString >::iterator itLT = looseTightDef.begin(); itLT != looseTightDef.end(); itLT++){
      if((*itLT).EqualTo("2L2J")){
	if ( std::find(cutsregions.begin(), cutsregions.end(), "njet30_geq_2") != cutsregions.end() )break;
	cutsregions.push_back("njet30_geq_2");
	cutsregions.push_back("zveto20");
	cutsregions.push_back("zveto20_bjet_eq_0");
	cutsregions.push_back("zveto20_bjet_geq_1");
	cutsregions.push_back("zveto20_bjet_leq_1"); // new
	cutsregions.push_back("met_g_40");
	cutsregions.push_back("met_g_40_njet30_geq_2");
	cutsregions.push_back("met_g_40_njet30_geq_2_ecids");
	cutsregions.push_back("met_g_100_metsig_g_6_bjet_leq_1_zveto20_njet30_geq_2");
      }else if((*itLT).EqualTo("ZMET")){
	cutsregions.push_back("SR1");
	cutsregions.push_back("SR2");
	cutsregions.push_back("SR3");
	cutsregions.push_back("VR_Z");
	cutsregions.push_back("VR_Top");
	cutsregions.push_back("VR_Dib");
	cutsregions.push_back("CR_Z");
	cutsregions.push_back("CR_Top");
	cutsregions.push_back("CR_Dib");
      }else{
	if ( std::find(cutsregions.begin(), cutsregions.end(), "CR_Dib") != cutsregions.end() )break;
      	// cutsregions.push_back("zveto20_bjet_eq_0");
      	// cutsregions.push_back("zveto20_bjet_geq_1");
      	// cutsregions.push_back("zveto20_njet20_eq_0");
      	// cutsregions.push_back("zveto20_njet20_l_2");
	cutsregions.push_back("CR_Dib");
	cutsregions.push_back("CR_Top");
	cutsregions.push_back("DF_1bjet");
	cutsregions.push_back("DF_0jet");
      }
    }
	
	
    


    // if 2L0J
    //cutsregions.push_back("metsig_g_3");
    //cutsregions.push_back("metsig_g_3_nsgjet_geq_0");	
    //cutsregions.push_back("metsig_g_3_nsgjet_eq_0");
    //cutsregions.push_back("metsig_g_3_nbjet_geq_1");

    /**
       cutsregions.push_back("SR-SF-1J");
       cutsregions.push_back("SR-DF-1J");
       cutsregions.push_back("SR-SF-0J");
       cutsregions.push_back("SR-DF-0J");

       cutsregions.push_back("CRVZ");
       cutsregions.push_back("CRWW");
       cutsregions.push_back("CRTOP");
       cutsregions.push_back("CRTOPL");
    */


    // cutsregion.push_back("VRWW-0J");
    // cutsregion.push_back("VRWW-1J");
    // cutsregion.push_back("VRVZ");
    // cutsregion.push_back("VRTOPL");
    // cutsregion.push_back("VRTOPH");
    // cutsregion.push_back("VRTOPWW");
  }
  if(looseTightDef.size() <= 1){
   
    
    cutsregions.push_back("FAKE2L04");
    //if(!doCoreRegions){
    if(doAnalysis){
      cutsregions.push_back("M_REAL2L02");
      cutsregions.push_back("M_FAKE2L20");
      cutsregions.push_back("M_FAKE2L21");
      cutsregions.push_back("E_FAKE2L21");
      cutsregions.push_back("E_REAL2L02");
      cutsregions.push_back("E_FAKE2L23");
      cutsregions.push_back("E_FAKE2L20");
    }
    //}
   
    //if(!doCoreRegions){
    //    cutsregions.push_back("E_FAKE2L27");
    
    
    //}
  }
  //cutsregions.push_back()
  // if(doThreeLep){
  //   cutsregions.push_back("3LRJR_STD");
  //   cutsregions.push_back("3LRJR_STD_noB");
  //   cutsregions.push_back("3LRJR_CMP");
  // }
  // if(!doThreeLep || doRJR){
  //   cutsregions.push_back("2L2J-VR-com");
  // }
  char reg[80];
  TH2F *h2;

  Float_t bins[] = { 25, 30, 35, 40, 45, 50, 60, 100, 1000 };
  Int_t  binnum = 9;//sizeof(bins)/sizeof(Float_t) - 1; // or just = 9

  for(it1 = controlregions.begin(); it1 != controlregions.end(); it1++){

    if(doFakeEst){
      for(std::vector< TString >::iterator itLT = looseTightDef.begin(); itLT != looseTightDef.end(); itLT++){
	TString LTDef = (*itLT);//looseTightDef.at(0);
	cout << "Creating" <<Form("%s_%s_stat_up",(*it1).Data(),LTDef.Data())<<endl;
	h2 = new TH2F(Form("%s_%s_stat_up",(*it1).Data(),LTDef.Data()),"",binnum, bins,binnum, bins); fOutput->AddLast(h2);    h_syst_MM[(*it1)+ "_" + LTDef+"_stat_up"]   = h2;
	h2 = new TH2F(Form("%s_%s_stat_dw",(*it1).Data(),LTDef.Data()),"",binnum, bins,binnum, bins); fOutput->AddLast(h2);    h_syst_MM[(*it1)+ "_" + LTDef+"_stat_dw"]   = h2;
	h2 = new TH2F(Form("%s_%s_wgt_up",(*it1).Data(),LTDef.Data()),"",binnum, bins,binnum, bins);  fOutput->AddLast(h2);     h_syst_MM[(*it1)+ "_" + LTDef+"_wgt_up"]    = h2;
	h2 = new TH2F(Form("%s_%s_wgt_dw",(*it1).Data(),LTDef.Data()),"",binnum, bins,binnum, bins);  fOutput->AddLast(h2);     h_syst_MM[(*it1)+ "_" + LTDef+"_wgt_dw"]    = h2;
	h2 = new TH2F(Form("%s_%s_mc_up",(*it1).Data(),LTDef.Data()),"",binnum, bins,binnum, bins);   fOutput->AddLast(h2);      h_syst_MM[(*it1)+ "_" + LTDef+"_mc_up"]     = h2;
	h2 = new TH2F(Form("%s_%s_mc_dw",(*it1).Data(),LTDef.Data()),"",binnum, bins,binnum, bins);   fOutput->AddLast(h2);      h_syst_MM[(*it1)+ "_" + LTDef+"_mc_dw"]     = h2;
	h2 = new TH2F(Form("%s_%s_xsec_up",(*it1).Data(),LTDef.Data()),"",binnum, bins,binnum, bins); fOutput->AddLast(h2);    h_syst_MM[(*it1)+ "_" + LTDef+"_xsec_up"]   = h2;
	h2 = new TH2F(Form("%s_%s_xsec_dw",(*it1).Data(),LTDef.Data()),"",binnum, bins,binnum, bins); fOutput->AddLast(h2);    h_syst_MM[(*it1)+ "_" + LTDef+"_xsec_dw"]   = h2;
	h2 = new TH2F(Form("%s_%s_total_up",(*it1).Data(),LTDef.Data()),"",binnum, bins,binnum, bins); fOutput->AddLast(h2);    h_syst_MM[(*it1)+ "_" + LTDef+"_total_up"]   = h2;
	h2 = new TH2F(Form("%s_%s_total_dw",(*it1).Data(),LTDef.Data()),"",binnum, bins,binnum, bins); fOutput->AddLast(h2);    h_syst_MM[(*it1)+ "_" + LTDef+"_total_dw"]   = h2;
      }
    }

    cout<<"controlregions "<<(*it1).Data()<<endl;

    if(!doCoreRegions){
      sprintf(reg,"%s",(*it1).Data());
      MM_regions.push_back(reg);
      if(doFakeEst){
	sprintf(reg,"%s_estfake",(*it1).Data());
	MM_regions.push_back(reg);
	sprintf(reg,"%s_estfake_UP",(*it1).Data());
	MM_regions.push_back(reg);
	sprintf(reg,"%s_estfake_DW",(*it1).Data());
	MM_regions.push_back(reg);
      }
      if(!isData){
	//MM_regions.push_back(reg);
	sprintf(reg,"%s_trueFAKE",(*it1).Data());
	MM_regions.push_back(reg);
	sprintf(reg,"%s_trueREAL",(*it1).Data());
	MM_regions.push_back(reg);
      }

      for(it2 = cutsregions.begin(); it2 != cutsregions.end(); it2++){

	if((*it2).Contains("-SF-") && (*it1).Contains("EM"))continue;
	if((*it2).Contains("-DF-") && !(*it1).Contains("EM"))continue;

	sprintf(reg,"%s_%s",(*it1).Data(),(*it2).Data());
	//if(!(std::count(MM_regions.begin(), MM_regions.end(), reg))) {
	MM_regions.push_back(reg);
	//}
	if(doFakeEst){
	  //cout<<"Adding "<<reg<<endl;
	  sprintf(reg,"%s_%s_estfake",(*it1).Data(),(*it2).Data());
	  MM_regions.push_back(reg);
	  sprintf(reg,"%s_%s_estfake_UP",(*it1).Data(),(*it2).Data());
	  MM_regions.push_back(reg);
	  sprintf(reg,"%s_%s_estfake_DW",(*it1).Data(),(*it2).Data());
	  MM_regions.push_back(reg);
	}
	if(!isData){
	  sprintf(reg,"%s_%s_trueFAKE",(*it1).Data(),(*it2).Data());
	  MM_regions.push_back(reg);
	  sprintf(reg,"%s_%s_trueREAL",(*it1).Data(),(*it2).Data());
	  MM_regions.push_back(reg);
	  sprintf(reg,"%s_%s_trueCF",(*it1).Data(),(*it2).Data());
	  MM_regions.push_back(reg);
	}
      }
    }//else{
     // if((*it1).EqualTo("MMSS") && !(*it2).Contains("2L4"))continue;
    //}
    //if(!((*it1).EqualTo("EE") || ((*it1).EqualTo("MM"))))continue;
    for(it3 = lepregions.begin(); it3 != lepregions.end(); it3++){

      cout<<"lepregions "<<(*it3).Data()<<endl;

      if((*it3).EqualTo("E") && !((*it1).Contains("EE") || (*it1).Contains("EM") || (*it1).Contains("ALL") || (*it1).Contains("EEE") || (*it1).Contains("MME") || (*it1).Contains("EEM")))continue;
      if((*it3).EqualTo("M") && !((*it1).Contains("MM") || (*it1).Contains("EM") || (*it1).Contains("ALL") || (*it1).Contains("MMM") || (*it1).Contains("MME") || (*it1).Contains("EEM")))continue;

      cout<<"Doing "<<(*it3).Data()<<endl;

      for(it2 = fakeregions.begin(); it2 != fakeregions.end(); it2++){
	if(doCoreRegions){
	  if((*it1).EqualTo("MMSS") && !(*it2).Contains("2L4"))continue;
	  if((*it3).EqualTo("E") && !((*it2).Contains("2L40") || (*it2).Contains("2L41") || (*it2).Contains("2L42") || (*it2).Contains("2L43") || (*it2).Contains("2L44") || (*it2).Contains("2L45") || (*it2).Contains("2L20") || (*it2).Contains("2L21") || (*it2).Contains("2L30") || (*it2).Contains("2L49") || (*it2).Contains("2L02") || (*it2).Contains("2L04") || (*it2).Contains("2L05") || (*it2).Contains("2L23") || (*it2).Contains("2L24") || (*it2).Contains("2L25") || (*it2).Contains("2L26") || (*it2).Contains("2L27")))continue;
	  if((*it3).EqualTo("M") && !((*it2).Contains("2L49") || (*it2).Contains("2L43") || (*it2).Contains("2L44") || (*it2).Contains("2L45") || (*it2).Contains("2L20") || (*it2).Contains("2L21") || (*it2).Contains("2L30") || (*it2).Contains("2L02") || (*it2).Contains("2L04") || (*it2).Contains("2L05")))continue;
	  //if(!((*it1).EqualTo("ALL") || (*it1).EqualTo("EESS")))continue;
	  if(((*it2).Contains("2L4") && (*it1).Contains("OS")))continue;
	}

	
	
	sprintf(reg,"%s_%s_FAKE%s",(*it1).Data(),(*it3).Data(),(*it2).Data());
	MM_regions.push_back(reg);
	if(!isData){
	  sprintf(reg,"%s_%s_trueLightFlavorDecay_FAKE%s",(*it1).Data(),(*it3).Data(),(*it2).Data());
	  MM_regions.push_back(reg);
	  sprintf(reg,"%s_%s_trueChargeFlipIsoElectron_FAKE%s",(*it1).Data(),(*it3).Data(),(*it2).Data());
	  /**if(!doCoreRegions)*/MM_regions.push_back(reg);
	  sprintf(reg,"%s_%s_truePromptPhotonConversion_FAKE%s",(*it1).Data(),(*it3).Data(),(*it2).Data());
	  /**if(!doCoreRegions)*/MM_regions.push_back(reg);
	  sprintf(reg,"%s_%s_trueBHadronDecay_FAKE%s",(*it1).Data(),(*it3).Data(),(*it2).Data());
	  /**if(!doCoreRegions)*/MM_regions.push_back(reg);
	  sprintf(reg,"%s_%s_CHadronDecay_FAKE%s",(*it1).Data(),(*it3).Data(),(*it2).Data());
	  /**if(!doCoreRegions)*/MM_regions.push_back(reg);
	  sprintf(reg,"%s_%s_ElectronFromMuon_FAKE%s",(*it1).Data(),(*it3).Data(),(*it2).Data());
	  /**if(!doCoreRegions)*/MM_regions.push_back(reg);
	  sprintf(reg,"%s_%s_trueFAKE_FAKE%s",(*it1).Data(),(*it3).Data(),(*it2).Data());
	  MM_regions.push_back(reg);
	  sprintf(reg,"%s_%s_trueREAL_FAKE%s",(*it1).Data(),(*it3).Data(),(*it2).Data());
	  MM_regions.push_back(reg);
	  sprintf(reg,"%s_%s_trueUnknown_FAKE%s",(*it1).Data(),(*it3).Data(),(*it2).Data());
	  /**if(!doCoreRegions)*/MM_regions.push_back(reg);
	  sprintf(reg,"%s_%s_KnownUnknown_FAKE%s",(*it1).Data(),(*it3).Data(),(*it2).Data());
	  /**if(!doCoreRegions)*/MM_regions.push_back(reg);
	}

	if((*it2).Contains("2L05") && (*it1).EqualTo("ALL")){
	  for(mapit1 = uncert.begin(); mapit1 != uncert.end(); mapit1++){
	    TString unc_name = mapit1->first;
	    if((*it3).EqualTo("E") && unc_name.Contains("MUON_"))continue;
	    if((*it3).EqualTo("M") && unc_name.Contains("EL_"))continue;
	    sprintf(reg,"%s_%s_FAKE%s_%s",(*it1).Data(),(*it3).Data(),(*it2).Data(),unc_name.Data());
	    MM_regions.push_back(reg);
	    if(!isData){
	      sprintf(reg,"%s_%s_trueLightFlavorDecay_FAKE%s_%s",(*it1).Data(),(*it3).Data(),(*it2).Data(),unc_name.Data());
	      MM_regions.push_back(reg);
	      sprintf(reg,"%s_%s_trueChargeFlipIsoElectron_FAKE%s_%s",(*it1).Data(),(*it3).Data(),(*it2).Data(),unc_name.Data());
	      /**if(!doCoreRegions)*/MM_regions.push_back(reg);
	      sprintf(reg,"%s_%s_truePromptPhotonConversion_FAKE%s_%s",(*it1).Data(),(*it3).Data(),(*it2).Data(),unc_name.Data());
	      /**if(!doCoreRegions)*/MM_regions.push_back(reg);
	      sprintf(reg,"%s_%s_trueBHadronDecay_FAKE%s_%s",(*it1).Data(),(*it3).Data(),(*it2).Data(),unc_name.Data());
	      /**if(!doCoreRegions)*/MM_regions.push_back(reg);
	      sprintf(reg,"%s_%s_CHadronDecay_FAKE%s_%s",(*it1).Data(),(*it3).Data(),(*it2).Data(),unc_name.Data());
	      /**if(!doCoreRegions)*/MM_regions.push_back(reg);
	      sprintf(reg,"%s_%s_ElectronFromMuon_FAKE%s_%s",(*it1).Data(),(*it3).Data(),(*it2).Data(),unc_name.Data());
	      /**if(!doCoreRegions)*/MM_regions.push_back(reg);
	      sprintf(reg,"%s_%s_trueFAKE_FAKE%s_%s",(*it1).Data(),(*it3).Data(),(*it2).Data(),unc_name.Data());
	      MM_regions.push_back(reg);
	      sprintf(reg,"%s_%s_trueREAL_FAKE%s_%s",(*it1).Data(),(*it3).Data(),(*it2).Data(),unc_name.Data());
	      MM_regions.push_back(reg);
	      sprintf(reg,"%s_%s_trueUnknown_FAKE%s_%s",(*it1).Data(),(*it3).Data(),(*it2).Data(),unc_name.Data());
	      /**if(!doCoreRegions)*/MM_regions.push_back(reg);
	      sprintf(reg,"%s_%s_KnownUnknown_FAKE%s_%s",(*it1).Data(),(*it3).Data(),(*it2).Data(),unc_name.Data());
	      /**if(!doCoreRegions)*/MM_regions.push_back(reg);
	    }

	  }	
	}
      }
      for(it2 = realregions.begin(); it2 != realregions.end(); it2++){
	if(doCoreRegions){
	  if((*it1).EqualTo("MMOS")  && !((*it2).Contains("2L01") ||(*it2).Contains("2L02") || (*it2).Contains("2L03") || (*it2).Contains("2L12") || (*it2).Contains("2L13") || (*it2).Contains("2L04") || (*it2).Contains("2L05") || (*it2).Contains("2L06") ||
					  (*it2).Contains("2L07") || (*it2).Contains("2L08") || (*it2).Contains("2L09") || (*it2).Contains("2L10") || (*it2).Contains("2L11") || (*it2).Contains("2L12") || (*it2).Contains("2L13") ))continue;
	  if(((*it1).Contains("EE") || (*it1).Contains("EM")) && !((*it2).Contains("2L01") || (*it2).Contains("2L02") || (*it2).Contains("2L03") || (*it2).Contains("2L12") || (*it2).Contains("2L13") ||(*it2).Contains("2L04") || (*it2).Contains("2L05") ||
								   (*it2).Contains("2L06") || (*it2).Contains("2L07") || (*it2).Contains("2L08") || (*it2).Contains("2L09") || (*it2).Contains("2L10") || (*it2).Contains("2L11") ))continue;
	  if(((*it1).Contains("EEE") || (*it1).Contains("EEM") || (*it1).Contains("MME") || (*it1).Contains("MMM")) && !((*it2).Contains("2L01") || (*it2).Contains("2L02") || (*it2).Contains("2L12") || (*it2).Contains("2L13") ||
															 (*it2).Contains("2L04") || (*it2).Contains("2L05") ||
															 (*it2).Contains("2L06") || (*it2).Contains("2L07") ||
															 (*it2).Contains("2L09") || (*it2).Contains("2L10") || 
															 (*it2).Contains("2L11") || (*it2).Contains("2L12") || (*it2).Contains("2L13")))continue;
	  if((*it1).EqualTo("ALL"))continue;
	}
	
	//cout<<"Adding "<<(*it1).Data()<<(*it3).Data()<<(*it2).Data()<<endl;
	sprintf(reg,"%s_%s_REAL%s",(*it1).Data(),(*it3).Data(),(*it2).Data());
	MM_regions.push_back(reg);
	if(!isData){
	  /**
	     sprintf(reg,"%s_%s_trueLF_REAL%s",(*it1).Data(),(*it3).Data(),(*it2).Data());
	     MM_regions.push_back(reg);
	     sprintf(reg,"%s_%s_trueCF_REAL%s",(*it1).Data(),(*it3).Data(),(*it2).Data());
	     if(!doCoreRegions)MM_regions.push_back(reg);
	     sprintf(reg,"%s_%s_trueCO_REAL%s",(*it1).Data(),(*it3).Data(),(*it2).Data());
	     if(!doCoreRegions)MM_regions.push_back(reg);
	     sprintf(reg,"%s_%s_trueHF_REAL%s",(*it1).Data(),(*it3).Data(),(*it2).Data());
	     if(!doCoreRegions)MM_regions.push_back(reg);
	  */
	  sprintf(reg,"%s_%s_trueFAKE_REAL%s",(*it1).Data(),(*it3).Data(),(*it2).Data());
	  MM_regions.push_back(reg);
	  sprintf(reg,"%s_%s_trueREAL_REAL%s",(*it1).Data(),(*it3).Data(),(*it2).Data());
	  MM_regions.push_back(reg);
	  sprintf(reg,"%s_%s_trueCF_REAL%s",(*it1).Data(),(*it3).Data(),(*it2).Data());
	  MM_regions.push_back(reg);
	  //sprintf(reg,"%s_%s_trueUNKNOWN_REAL%s",(*it1).Data(),(*it3).Data(),(*it2).Data());
	  //if(!doCoreRegions)MM_regions.push_back(reg);
	}


	if((*it2).Contains("2L04") && ((*it1).EqualTo("EEOS") || (*it1).EqualTo("MMOS"))){
	  for(mapit1 = uncert.begin(); mapit1 != uncert.end(); mapit1++){
	    TString unc_name = mapit1->first;
	    if((*it1).EqualTo("EEOS") && unc_name.Contains("MUON_"))continue;
	    if((*it1).EqualTo("MMOS") && unc_name.Contains("EL_"))continue;
	    sprintf(reg,"%s_%s_REAL%s_%s",(*it1).Data(),(*it3).Data(),(*it2).Data(),unc_name.Data());
	    MM_regions.push_back(reg);
	    if(!isData){
	      sprintf(reg,"%s_%s_trueFAKE_REAL%s_%s",(*it1).Data(),(*it3).Data(),(*it2).Data(),unc_name.Data());
	      MM_regions.push_back(reg);
	      sprintf(reg,"%s_%s_trueREAL_REAL%s_%s",(*it1).Data(),(*it3).Data(),(*it2).Data(),unc_name.Data());
	      MM_regions.push_back(reg);
	      sprintf(reg,"%s_%s_trueCF_REAL%s_%s",(*it1).Data(),(*it3).Data(),(*it2).Data(),unc_name.Data());
	      MM_regions.push_back(reg);
	      //sprintf(reg,"%s_%s_trueUNKNOWN_REAL%s",(*it1).Data(),(*it3).Data(),(*it2).Data());
	      //if(!doCoreRegions)MM_regions.push_back(reg);
	    }
	  }
	}	
      }
    }
  }
  
  TH1F *h;

  h_lep_lepECIDS = new TH1F("h_lep_lepECIDS","h_lep_lepECIDS",2,0,2); fOutput->AddLast(h_lep_lepECIDS);

  h_runNumber18 = new TH1F("h_runNumber18","h_runNumber18",15407,348885,364292); fOutput->AddLast(h_runNumber18);

  h_fake_weights = new TH1F("h_fake_weights","h_fake_weights",2000,-100,100); fOutput->AddLast(h_fake_weights);
  h_HFreg_mindr_jet_probe = new TH1F("h_HFreg_mindr_jet_probe","h_HFreg_mindr_jet_probe",100,0,1); fOutput->AddLast(h_HFreg_mindr_jet_probe);
  h_HFreg_nprobes = new TH1F("h_HFreg_nprobes","h_HFreg_nprobes",10,0,10); fOutput->AddLast(h_HFreg_nprobes);
  h_HFreg_ntags_nprobes = new TH2F("h_HFreg_ntags_nprobes","h_HFreg_ntags_nprobes",10,0,10,10,0,10); fOutput->AddLast(h_HFreg_ntags_nprobes);

  h_passed_trigMatch_2LTrig_fakeVR = new TH1F("h_passed_trigMatch_2LTrig_fakeVR","h_passed_trigMatch_2LTrig_fakeVR",2,0,2); fOutput->AddLast(h_passed_trigMatch_2LTrig_fakeVR);

  h_probe_el_hf_deltaR_ORjet  = new TH1F("h_probe_el_hf_deltaR_ORjet","h_probe_el_hf_deltaR_ORjet",100,0,1); fOutput->AddLast(h_probe_el_hf_deltaR_ORjet);
  h_probe_el_hf_deltaR_ORbjet  = new TH1F("h_probe_el_hf_deltaR_ORbjet","h_probe_el_hf_deltaR_ORbjet",100,0,1); fOutput->AddLast(h_probe_el_hf_deltaR_ORbjet);

  h_probe_mu_hf_deltaR_ORjet  = new TH1F("h_probe_mu_hf_deltaR_ORjet","h_probe_mu_hf_deltaR_ORjet",100,0,1); fOutput->AddLast(h_probe_mu_hf_deltaR_ORjet);
  h_probe_mu_hf_deltaR_ORbjet  = new TH1F("h_probe_mu_hf_deltaR_ORbjet","h_probe_mu_hf_deltaR_ORbjet",100,0,1); fOutput->AddLast(h_probe_mu_hf_deltaR_ORbjet);

  h_ptllboost = new TH1F("h_ptllboost","h_ptllboost",150,0,150); fOutput->AddLast(h_ptllboost);


  h_wgts_leptonWeight = new TH1F("h_wgts_leptonWeight","h_wgts_leptonWeight",10000,-5,5); fOutput->AddLast(h_wgts_leptonWeight );
  h_wgts_baselineleptonWeight = new TH1F("h_wgts_baselineleptonWeight","h_wgts_baselineleptonWeight",10000,-5,5); fOutput->AddLast(h_wgts_baselineleptonWeight );
  h_wgts_globalDiLepTrigSF = new TH1F("h_wgts_globalDiLepTrigSF","h_wgts_globalDiLepTrigSF",10000,-5,5); fOutput->AddLast(h_wgts_globalDiLepTrigSF );
  h_wgts_singleLepTrigSF = new TH1F("h_wgts_singleLepTrigSF","h_wgts_singleLepTrigSF",10000,-5,5); fOutput->AddLast(h_wgts_singleLepTrigSF );
  h_wgts_globalBaselineDiLepTrigSF = new TH1F("h_wgts_globalBaselineDiLepTrigSF","h_wgts_globalBaselineDiLepTrigSF",10000,-5,5); fOutput->AddLast(h_wgts_globalBaselineDiLepTrigSF);
  h_wgts_singleBaselineLepTrigSF = new TH1F("h_wgts_singleBaselineLepTrigSF","h_wgts_singleBaselineLepTrigSF",10000,-5,5); fOutput->AddLast(h_wgts_singleBaselineLepTrigSF );

  h_wgts_lepBLRecoSF_lep1 = new TH1F("h_wgts_lepBLRecoSF_lep1","h_wgts_lepBLRecoSF_lep1",10000,-5,5); fOutput->AddLast(h_wgts_lepBLRecoSF_lep1 );
  h_wgts_lepBLRecoSF_lep2 = new TH1F("h_wgts_lepBLRecoSF_lep2","h_wgts_lepBLRecoSF_lep2",10000,-5,5); fOutput->AddLast(h_wgts_lepBLRecoSF_lep2 );

  h_wgts_lepRecoSF_lep1 = new TH1F("h_wgts_lepRecoSF_lep1","h_wgts_lepRecoSF_lep1",10000,-5,5); fOutput->AddLast(h_wgts_lepRecoSF_lep1 );
  h_wgts_lepRecoSF_lep2 = new TH1F("h_wgts_lepRecoSF_lep2","h_wgts_lepRecoSF_lep2",10000,-5,5); fOutput->AddLast(h_wgts_lepRecoSF_lep2 );

  h_wgts_lepBLTrigSF_lep1 = new TH1F("h_wgts_lepBLTrigSF_lep1","h_wgts_lepBLTrigSF_lep1",10000,-5,5); fOutput->AddLast(h_wgts_lepBLTrigSF_lep1 );
  h_wgts_lepBLTrigSF_lep2 = new TH1F("h_wgts_lepBLTrigSF_lep2","h_wgts_lepBLTrigSF_lep2",10000,-5,5); fOutput->AddLast(h_wgts_lepBLTrigSF_lep2 );

  h_wgts_lepTrigSF_lep1 = new TH1F("h_wgts_lepTrigSF_lep1","h_wgts_lepTrigSF_lep1",10000,-5,5); fOutput->AddLast(h_wgts_lepTrigSF_lep1 );
  h_wgts_lepTrigSF_lep2 = new TH1F("h_wgts_lepTrigSF_lep2","h_wgts_lepTrigSF_lep2",10000,-5,5); fOutput->AddLast(h_wgts_lepTrigSF_lep2 );
  
  h_wgts_wgt = new TH1F("h_wgts_wgt","h_wgts_wgt",10000,-5,5); fOutput->AddLast(h_wgts_wgt );
  h_wgts_wgtloose = new TH1F("h_wgts_wgtloose","h_wgts_wgtloose",10000,-5,5); fOutput->AddLast(h_wgts_wgtloose );

  h_wgts_infilllandt_wgt = new TH1F("h_wgts_infilllandt_wgt","h_wgts_infilllandt_wgt",10000,-5,5); fOutput->AddLast(h_wgts_infilllandt_wgt );
  h_wgts_infilllandt_wgtloose = new TH1F("h_wgts_infilllandt_wgtloose","h_wgts_infilllandt_wgtloose",10000,-5,5); fOutput->AddLast(h_wgts_infilllandt_wgtloose );
 

  h_bornMass = new TH1F("h_bornMass","h_bornMass",1000,0,1000); fOutput->AddLast(h_bornMass);

  
  int nit = 0;
  for(it1 = MM_regions.begin(); it1 != MM_regions.end(); it1++){
    nit += 1;
    if(nit%100 == 0)cout<<"Number of historgrams stored "<<fOutput->GetSize()<<"  after " << nit << " of " << MM_regions.size() << " regions" << endl;
    for(it4 = looseTightDef.begin(); it4 != looseTightDef.end(); it4++){
      if(!doCoreRegions){
      
	if((isData && doFakeEst) || !isData){ //(!((*it1).Contains("FAKE") || (*it1).Contains("REAL")) || (((*it1).Contains("trueFAKE") || (*it1).Contains("trueREAL")) && !(*it1).Contains("_REAL") && !(*it1).Contains("_FAKE")))){

	


	  //cout<<"Filling region "<<((*it1)+"_"+(*it4)).Data()<<endl; 
	  h = new TH1F(Form("h_lep_pT0_%s",((*it1)+"_"+(*it4)).Data()),"",nPt,0,maxPt);   fOutput->AddLast(h);    h_lep_pT0[(*it1)+"_"+(*it4)]   = h;
	  h = new TH1F(Form("h_lep_pT1_%s",((*it1)+"_"+(*it4)).Data()),"",nPt,0,maxPt);   fOutput->AddLast(h);    h_lep_pT1[(*it1)+"_"+(*it4)]   = h;
	  h = new TH1F(Form("h_lep_pT2_%s",((*it1)+"_"+(*it4)).Data()),"",nPt,0,maxPt);   fOutput->AddLast(h);    h_lep_pT2[(*it1)+"_"+(*it4)]   = h;
	  h = new TH1F(Form("h_lep_mt2_%s",((*it1)+"_"+(*it4)).Data()),"",nMT2,0,maxMT2); fOutput->AddLast(h);    h_lep_mt2[(*it1)+"_"+(*it4)]   = h;
	  if(doThreeLep){h = new TH1F(Form("h_lep_mtw_%s",((*it1)+"_"+(*it4)).Data()),"",nMT2,0,maxMT2); fOutput->AddLast(h);    h_lep_mtw[(*it1)+"_"+(*it4)]   = h;  }
	  h = new TH1F(Form("h_lep_met_%s",((*it1)+"_"+(*it4)).Data()),"",nPt,0,maxPt);   fOutput->AddLast(h);    h_lep_met[(*it1)+"_"+(*it4)]   = h;
	  h = new TH1F(Form("h_lep_metsig_%s",((*it1)+"_"+(*it4)).Data()),"",nMETSIG,0,maxMETSIG);   fOutput->AddLast(h);    h_lep_metsig[(*it1)+"_"+(*it4)]   = h;
	  h = new TH1F(Form("h_lep_mll_%s",((*it1)+"_"+(*it4)).Data()),"",nPt,0,maxPt);   fOutput->AddLast(h);    h_lep_mll[(*it1)+"_"+(*it4)]   = h;
	  h = new TH1F(Form("h_lep_njet30_%s",((*it1)+"_"+(*it4)).Data()),"",10,0,10);   fOutput->AddLast(h);    h_lep_njet30[(*it1)+"_"+(*it4)]   = h;
	  //if(!((*it1)+"_"+(*it4)).Contains("estfake")){

	  h = new TH1F(Form("h_lep_nbj77_%s",((*it1)+"_"+(*it4)).Data()),"",10,0,10);     fOutput->AddLast(h);    h_lep_nbj77[(*it1)+"_"+(*it4)] = h;
	  h = new TH1F(Form("h_lep_nbj85_%s",((*it1)+"_"+(*it4)).Data()),"",10,0,10);     fOutput->AddLast(h);    h_lep_nbj85[(*it1)+"_"+(*it4)] = h;

	  //}
	  if(/**isData && */doFakeEst && ((*it1)+"_"+(*it4)).Contains("estfake") && (!((*it1)+"_"+(*it4)).Contains("estfake_UP") && !((*it1)+"_"+(*it4)).Contains("estfake_DW"))){
	    h = new TH1F(Form("h_MMevtype_%s",((*it1)+"_"+(*it4)).Data()),"",4,0,4);   fOutput->AddLast(h);    h_MMevtype[(*it1)+"_"+(*it4)]   = h;
	    h_MMevtype[(*it1)+"_"+(*it4)]->GetXaxis()->SetBinLabel(1,"TT");
	    h_MMevtype[(*it1)+"_"+(*it4)]->GetXaxis()->SetBinLabel(2,"Tl");
	    h_MMevtype[(*it1)+"_"+(*it4)]->GetXaxis()->SetBinLabel(3,"lT");
	    h_MMevtype[(*it1)+"_"+(*it4)]->GetXaxis()->SetBinLabel(4,"ll");
	  }

	}
      
      }
      if(((*it1)).Contains("estfake"))continue;


      if((((*it1).Contains("FAKE2L20") || (*it1).Contains("FAKE2L21")) && (*it1).Contains("ALL")) || 
	 (((*it1).Contains("FAKE2L27") || (*it1).Contains("FAKE2L25")) && (*it1).Contains("ALL")) || 
	 (((*it1).Contains("FAKE2L05") || (*it1).Contains("FAKE2L04")) && (*it1).Contains("EE")) || 
	 ((*it1).Contains("REAL2L02") && ((*it1).Contains("EEOS") || (*it1).Contains("MMOS")))){
	//cout<<"Fillin h_lep_nbljetpreOR with "<<(*it1).Data()<<endl;
	h = new TH1F(Form("h_lep_nbljetpreOR_%s",((*it1)+"_"+(*it4)).Data()),"",50,0,50);   fOutput->AddLast(h);    h_lep_nbljetpreOR[(*it1)+"_"+(*it4)]   = h;
      }
    } 
    for(it3 = trigregions.begin(); it3 != trigregions.end(); it3++){
      
      
      TString key1 = (*it1);
      TString key  = (*it1);
      //if(doCoreRegions && !(key.Contains("REAL2L02") || key.Contains("FAKE2L20") || key.Contains("FAKE2L21") || key.Contains("FAKE2L04") || key.Contains("FAKE2L23")))continue;
      
      TString add = "";
      if((*it3).Length() > 0){
	add  = "_"+(*it3);
	key1 += "_"+(*it3);

      }
      key = key1;
      for(it4 = looseTightDef.begin(); it4 != looseTightDef.end(); it4++){

	key = key1 + "_"+(*it4);
	//cout<<"add = "<<key.Data()<<endl;
	//	cout<<"key = "<<key.Data()<<endl;



	if(!isData){
	  //if(!(*it1).Contains("true")){
	  //if(!(*it1).Contains("trueFAKE"))continue;
	  h = new TH1F(Form("h_lep_TCF_nL_%s",(key).Data()),"",truthVECIFF.size(),0,truthVECIFF.size()); fOutput->AddLast(h);    h_lep_TCF_nL[key]   = h;
	  h = new TH1F(Form("h_lep_TCF_nT_%s",(key).Data()),"",truthVECIFF.size(),0,truthVECIFF.size()); fOutput->AddLast(h);    h_lep_TCF_nT[key]   = h;

	  if(!key.Contains("_")){
	    for(it2 = ptregions.begin(); it2 != ptregions.end(); it2++){
	      h2 = new TH2F(Form("h_lep_type_origin_nL_%s_%s",(key).Data(),(*it2).Data()),"",41,0,41,46,0,46); fOutput->AddLast(h2);    h_lep_type_origin_nL[key+"_"+*it2]   = h2;
	      h2 = new TH2F(Form("h_lep_type_origin_nT_%s_%s",(key).Data(),(*it2).Data()),"",41,0,41,46,0,46); fOutput->AddLast(h2);    h_lep_type_origin_nT[key+"_"+*it2]   = h2;
	    }
	  }

	  h2 = new TH2F(Form("h_lep_type_origin_nL_%s",(key).Data()),"",41,0,41,46,0,46); fOutput->AddLast(h2);    h_lep_type_origin_nL[key]   = h2;
	  h2 = new TH2F(Form("h_lep_type_origin_nT_%s",(key).Data()),"",41,0,41,46,0,46); fOutput->AddLast(h2);    h_lep_type_origin_nT[key]   = h2;

	  //cout<<" Now filling "<<(key).Data()<<endl;
	  
	  h2 = new TH2F(Form("h_lep_pT_origin_nL_%s",(key).Data()),"",nPt,0,maxPt,46,0,46); fOutput->AddLast(h2);    h_lep_pT_origin_nL[key]   = h2;
	  h2 = new TH2F(Form("h_lep_pT_origin_nT_%s",(key).Data()),"",nPt,0,maxPt,46,0,46); fOutput->AddLast(h2);    h_lep_pT_origin_nT[key]   = h2;

	  h2 = new TH2F(Form("h_lep_pT_firstEgMotherPdgId_nL_%s",(key).Data()),"",nPt,0,maxPt,82,-41,41); fOutput->AddLast(h2);    h_lep_pT_firstEgMotherPdgId_nL[key]   = h2;
	  h2 = new TH2F(Form("h_lep_pT_firstEgMotherPdgId_nT_%s",(key).Data()),"",nPt,0,maxPt,82,-41,41); fOutput->AddLast(h2);    h_lep_pT_firstEgMotherPdgId_nT[key]   = h2;

	  h2 = new TH2F(Form("h_lep_pT_firstEgMotherO_nL_%s",(key).Data()),"",nPt,0,maxPt,41,0,41); fOutput->AddLast(h2);    h_lep_pT_firstEgMotherO_nL[key]   = h2;
	  h2 = new TH2F(Form("h_lep_pT_firstEgMotherO_nT_%s",(key).Data()),"",nPt,0,maxPt,41,0,41); fOutput->AddLast(h2);    h_lep_pT_firstEgMotherO_nT[key]   = h2;

	  h2 = new TH2F(Form("h_lep_pT_type_nL_%s",(key).Data()),"",nPt,0,maxPt,41,0,41); fOutput->AddLast(h2);    h_lep_pT_type_nL[key]   = h2;
	  h2 = new TH2F(Form("h_lep_pT_type_nT_%s",(key).Data()),"",nPt,0,maxPt,41,0,41); fOutput->AddLast(h2);    h_lep_pT_type_nT[key]   = h2;

	  h2 = new TH2F(Form("h_lep_pT_firstEgMotherT_nL_%s",(key).Data()),"",nPt,0,maxPt,41,0,41); fOutput->AddLast(h2);    h_lep_pT_firstEgMotherT_nL[key]   = h2;
	  h2 = new TH2F(Form("h_lep_pT_firstEgMotherT_nT_%s",(key).Data()),"",nPt,0,maxPt,41,0,41); fOutput->AddLast(h2);    h_lep_pT_firstEgMotherT_nT[key]   = h2;


	  h2 = new TH2F(Form("h_lep_mT2_TCF_nL_%s",(key).Data()),"",nMT2,0,maxMT2,truthVECIFF.size(),0,truthVECIFF.size()); fOutput->AddLast(h2);    h_lep_mT2_TCF_nL[key]   = h2;
	  h2 = new TH2F(Form("h_lep_mT2_TCF_nT_%s",(key).Data()),"",nMT2,0,maxMT2,truthVECIFF.size(),0,truthVECIFF.size()); fOutput->AddLast(h2);    h_lep_mT2_TCF_nT[key]   = h2;
	  
	  if(doThreeLep){
	    h2 = new TH2F(Form("h_lep_mTW_TCF_nL_%s",(key).Data()),"",nMT2,0,maxMT2,truthVECIFF.size(),0,truthVECIFF.size()); fOutput->AddLast(h2);    h_lep_mTW_TCF_nL[key]   = h2;
	    h2 = new TH2F(Form("h_lep_mTW_TCF_nT_%s",(key).Data()),"",nMT2,0,maxMT2,truthVECIFF.size(),0,truthVECIFF.size()); fOutput->AddLast(h2);    h_lep_mTW_TCF_nT[key]   = h2;
	  }
	  h2 = new TH2F(Form("h_lep_njet_TCF_nL_%s",(key).Data()),"",10,0,10,truthVECIFF.size(),0,truthVECIFF.size()); fOutput->AddLast(h2);    h_lep_njet_TCF_nL[key]   = h2;
	  h2 = new TH2F(Form("h_lep_njet_TCF_nT_%s",(key).Data()),"",10,0,10,truthVECIFF.size(),0,truthVECIFF.size()); fOutput->AddLast(h2);    h_lep_njet_TCF_nT[key]   = h2;

	  h2 = new TH2F(Form("h_lep_metsig_TCF_nL_%s",(key).Data()),"",nMETSIG,0,maxMETSIG,truthVECIFF.size(),0,truthVECIFF.size()); fOutput->AddLast(h2);    h_lep_metsig_TCF_nL[key]   = h2;
	  h2 = new TH2F(Form("h_lep_metsig_TCF_nT_%s",(key).Data()),"",nMETSIG,0,maxMETSIG,truthVECIFF.size(),0,truthVECIFF.size()); fOutput->AddLast(h2);    h_lep_metsig_TCF_nT[key]   = h2;


	  h2 = new TH2F(Form("h_lep_pT_TCF_nL_%s",(key).Data()),"",nPt,0,maxPt,truthVECIFF.size(),0,truthVECIFF.size()); fOutput->AddLast(h2);    h_lep_pT_TCF_nL[key]   = h2;
	  h2 = new TH2F(Form("h_lep_pT_TCF_nT_%s",(key).Data()),"",nPt,0,maxPt,truthVECIFF.size(),0,truthVECIFF.size()); fOutput->AddLast(h2);    h_lep_pT_TCF_nT[key]   = h2;



	  int ibin = 1;
	  for(it2 = truthVECIFF.begin(); it2 != truthVECIFF.end(); it2++){
	    h_lep_TCF_nL[(key)]->GetXaxis()->SetBinLabel(ibin,(*it2));
	    h_lep_TCF_nT[(key)]->GetXaxis()->SetBinLabel(ibin,(*it2));
	    h_lep_mT2_TCF_nL[(key)]->GetYaxis()->SetBinLabel(ibin,(*it2));
	    h_lep_mT2_TCF_nT[(key)]->GetYaxis()->SetBinLabel(ibin,(*it2));
	    if(doThreeLep){
	      h_lep_mTW_TCF_nL[(key)]->GetYaxis()->SetBinLabel(ibin,(*it2));
	      h_lep_mTW_TCF_nT[(key)]->GetYaxis()->SetBinLabel(ibin,(*it2));
	    }
	    h_lep_njet_TCF_nL[(key)]->GetYaxis()->SetBinLabel(ibin,(*it2));
	    h_lep_njet_TCF_nT[(key)]->GetYaxis()->SetBinLabel(ibin,(*it2));
	    h_lep_metsig_TCF_nL[(key)]->GetYaxis()->SetBinLabel(ibin,(*it2));
	    h_lep_metsig_TCF_nT[(key)]->GetYaxis()->SetBinLabel(ibin,(*it2));
	    h_lep_pT_TCF_nL[(key)]->GetYaxis()->SetBinLabel(ibin,(*it2));
	    h_lep_pT_TCF_nT[(key)]->GetYaxis()->SetBinLabel(ibin,(*it2));
	    ibin += 1;
	  }
	}

	// don't need to do SS/OS for the REAL/FAKE regions
	//if((/**(key).Contains("REAL") ||*/ (key).Contains("FAKE2L")) && ((key).Contains("OS") || (key).Contains("SS")))continue;

	if(doThreeLep && doRJR){
	  //cout<<"key = "<<(*it1).Data()<<endl;
	  h = new TH1F(Form("h_lep_pT0_nT_%s",(key).Data()),"",nPt,0,maxPt); fOutput->AddLast(h);    h_lep_pT0_nT[(key)]   = h;
	  h = new TH1F(Form("h_lep_pT0_nL_%s",(key).Data()),"",nPt,0,maxPt); fOutput->AddLast(h);    h_lep_pT0_nL[(key)]   = h;
	}

	h = new TH1F(Form("h_lep_MET_nL_%s",(key).Data()),"",200,0,1000); fOutput->AddLast(h);    h_lep_MET_nL[(key)]   = h;
	h = new TH1F(Form("h_lep_MET_nT_%s",(key).Data()),"",200,0,1000); fOutput->AddLast(h);    h_lep_MET_nT[(key)]   = h;
	
	h = new TH1F(Form("h_lep_nbjet_nL_%s",(key).Data()),"",10,0,10); fOutput->AddLast(h);    h_lep_nbjet_nL[(key)]   = h;
	h = new TH1F(Form("h_lep_nbjet_nT_%s",(key).Data()),"",10,0,10); fOutput->AddLast(h);    h_lep_nbjet_nT[(key)]   = h;

	h = new TH1F(Form("h_lep_pT_nT_%s",(key).Data()),"",nPt,0,maxPt); fOutput->AddLast(h);    h_lep_pT_nT[key]   = h;
	h = new TH1F(Form("h_lep_pT_nL_%s",(key).Data()),"",nPt,0,maxPt); fOutput->AddLast(h);    h_lep_pT_nL[key]   = h;
      
	for (TString bdt_name : bdt_vec) {
	  h = new TH1F(Form("h_evt_BDT_TT_%s_%s",bdt_name.Data(),(key).Data()),"",nBDT,0,maxBDT); fOutput->AddLast(h);    h_evt_BDT_TT[bdt_name+"_"+key]   = h;
	  h = new TH1F(Form("h_evt_BDT_Tl_%s_%s",bdt_name.Data(),(key).Data()),"",nBDT,0,maxBDT); fOutput->AddLast(h);    h_evt_BDT_Tl[bdt_name+"_"+key]   = h;
	  h = new TH1F(Form("h_evt_BDT_lT_%s_%s",bdt_name.Data(),(key).Data()),"",nBDT,0,maxBDT); fOutput->AddLast(h);    h_evt_BDT_lT[bdt_name+"_"+key]   = h;
	  h = new TH1F(Form("h_evt_BDT_ll_%s_%s",bdt_name.Data(),(key).Data()),"",nBDT,0,maxBDT); fOutput->AddLast(h);    h_evt_BDT_ll[bdt_name+"_"+key]   = h;
	  
	  h = new TH1F(Form("h_lep_BDT_nT_%s_%s",bdt_name.Data(),(key).Data()),"",nBDT,0,maxBDT); fOutput->AddLast(h);    h_lep_BDT_nT[bdt_name+"_"+key]   = h;
	  h = new TH1F(Form("h_lep_BDT_nL_%s_%s",bdt_name.Data(),(key).Data()),"",nBDT,0,maxBDT); fOutput->AddLast(h);    h_lep_BDT_nL[bdt_name+"_"+key]   = h;
	}
	
	h2 = new TH2F(Form("h_lep_pT_eta_nT_%s",(key).Data()),"",nPt,0,maxPt,nEta,-1*maxEta,maxEta); fOutput->AddLast(h2);    h_lep_pT_eta_nT[key]   = h2;
	h2 = new TH2F(Form("h_lep_pT_eta_nL_%s",(key).Data()),"",nPt,0,maxPt,nEta,-1*maxEta,maxEta); fOutput->AddLast(h2);    h_lep_pT_eta_nL[key]   = h2;

	//	cout<<"Adding key "<<(key).Data()<<endl;
	h = new TH1F(Form("h_lep_njet30_nT_%s",(key).Data()),"",10,0,10); fOutput->AddLast(h);    h_lep_njet30_nT[key]   = h;
	h = new TH1F(Form("h_lep_njet30_nL_%s",(key).Data()),"",10,0,10); fOutput->AddLast(h);    h_lep_njet30_nL[key]   = h;

	h2 = new TH2F(Form("h_lep_mu_pT_nL_%s",(key).Data()),"",80,0,80,100,0,1000); fOutput->AddLast(h2);    h_lep_mu_pT_nL[key]   = h2;
	h2 = new TH2F(Form("h_lep_mu_pT_nT_%s",(key).Data()),"",80,0,80,100,0,1000); fOutput->AddLast(h2);    h_lep_mu_pT_nT[key]   = h2;

	h = new TH1F(Form("h_lep_eta_nT_%s",(key).Data()),"",nEta,-1*maxEta,maxEta); fOutput->AddLast(h);    h_lep_eta_nT[key]   = h;
	h = new TH1F(Form("h_lep_eta_nL_%s",(key).Data()),"",nEta,-1*maxEta,maxEta); fOutput->AddLast(h);    h_lep_eta_nL[key]   = h;

	
	if(triglist.size() > 0){
	  h = new TH1F(Form("h_lep_trig_nT_lepnotmatched_%s",(key).Data()),"",triglist.size(),0,triglist.size()); fOutput->AddLast(h);    h_lep_trig_nT["lepnotmatched_"+key]   = h;
	  h = new TH1F(Form("h_lep_trig_nL_lepnotmatched_%s",(key).Data()),"",triglist.size(),0,triglist.size()); fOutput->AddLast(h);    h_lep_trig_nL["lepnotmatched_"+key]   = h;
	  h = new TH1F(Form("h_lep_trig_nT_lepmatched_%s",(key).Data()),"",triglist.size(),0,triglist.size()); fOutput->AddLast(h);    h_lep_trig_nT["lepmatched_"+key]   = h;
	  h = new TH1F(Form("h_lep_trig_nL_lepmatched_%s",(key).Data()),"",triglist.size(),0,triglist.size()); fOutput->AddLast(h);    h_lep_trig_nL["lepmatched_"+key]   = h;

	  
	  
	  for(int i=1; i<h_lep_trig_nT["lepnotmatched_"+key]->GetNbinsX()+1; i++){


	    //cout<<"Now filling key = "<<Form("h_lep_eta_nT_%s_lepmatched_%s",triglist.at(i-1).Data(),(key).Data())<<endl;
	    
	    
	    h_lep_trig_nT["lepnotmatched_"+key]->GetXaxis()->SetBinLabel(i,triglist.at(i-1));
	    //cout<<"e"<<endl;
	    h_lep_trig_nL["lepnotmatched_"+key]->GetXaxis()->SetBinLabel(i,triglist.at(i-1));
	    h_lep_trig_nT["lepmatched_"+key]->GetXaxis()->SetBinLabel(i,triglist.at(i-1));
	    h_lep_trig_nL["lepmatched_"+key]->GetXaxis()->SetBinLabel(i,triglist.at(i-1));
	    //cout<<"d"<<endl;
	    h = new TH1F(Form("h_lep_eta_nT_%s_lepmatched_%s",triglist.at(i-1).Data(),(key).Data()),"",nEta,-1*maxEta,maxEta); fOutput->AddLast(h);    h_lep_eta_nT[triglist.at(i-1)+"_lepmatched_"+key]   = h;
	    h = new TH1F(Form("h_lep_eta_nL_%s_lepmatched_%s",triglist.at(i-1).Data(),(key).Data()),"",nEta,-1*maxEta,maxEta); fOutput->AddLast(h);    h_lep_eta_nL[triglist.at(i-1)+"_lepmatched_"+key]   = h;
	    h = new TH1F(Form("h_lep_eta_nT_%s_lepnotmatched_%s",triglist.at(i-1).Data(),(key).Data()),"",nEta,-1*maxEta,maxEta); fOutput->AddLast(h);    h_lep_eta_nT[triglist.at(i-1)+"_lepnotmatched_"+key]   = h;
	    h = new TH1F(Form("h_lep_eta_nL_%s_lepnotmatched_%s",triglist.at(i-1).Data(),(key).Data()),"",nEta,-1*maxEta,maxEta); fOutput->AddLast(h);    h_lep_eta_nL[triglist.at(i-1)+"_lepnotmatched_"+key]   = h;

	    for (TString bdt_name : bdt_vec) {
	      h = new TH1F(Form("h_lep_BDT_nT_%s_%s_lepmatched_%s",bdt_name.Data(),triglist.at(i-1).Data(),(key).Data()),"",nBDT,0,maxBDT); fOutput->AddLast(h);    h_lep_BDT_nT[bdt_name+"_"+triglist.at(i-1)+"_lepmatched_"+key]   = h;
	      h = new TH1F(Form("h_lep_BDT_nL_%s_%s_lepmatched_%s",bdt_name.Data(),triglist.at(i-1).Data(),(key).Data()),"",nBDT,0,maxBDT); fOutput->AddLast(h);    h_lep_BDT_nL[bdt_name+"_"+triglist.at(i-1)+"_lepmatched_"+key]   = h;
	      h = new TH1F(Form("h_lep_BDT_nT_%s_%s_lepnotmatched_%s",bdt_name.Data(),triglist.at(i-1).Data(),(key).Data()),"",nBDT,0,maxBDT); fOutput->AddLast(h);    h_lep_BDT_nT[bdt_name+"_"+triglist.at(i-1)+"_lepnotmatched_"+key]   = h;
	      h = new TH1F(Form("h_lep_BDT_nL_%s_%s_lepnotmatched_%s",bdt_name.Data(),triglist.at(i-1).Data(),(key).Data()),"",nBDT,0,maxBDT); fOutput->AddLast(h);    h_lep_BDT_nL[bdt_name+"_"+triglist.at(i-1)+"_lepnotmatched_"+key]   = h;
	    }
	    
	    h = new TH1F(Form("h_lep_pT_nT_%s_lepmatched_%s",triglist.at(i-1).Data(),(key).Data()),"",nPt,0,maxPt); fOutput->AddLast(h);    h_lep_pT_nT[triglist.at(i-1)+"_lepmatched_"+key]   = h;
	    h = new TH1F(Form("h_lep_pT_nL_%s_lepmatched_%s",triglist.at(i-1).Data(),(key).Data()),"",nPt,0,maxPt); fOutput->AddLast(h);    h_lep_pT_nL[triglist.at(i-1)+"_lepmatched_"+key]   = h;
	    h = new TH1F(Form("h_lep_pT_nT_%s_lepnotmatched_%s",triglist.at(i-1).Data(),(key).Data()),"",nPt,0,maxPt); fOutput->AddLast(h);    h_lep_pT_nT[triglist.at(i-1)+"_lepnotmatched_"+key]   = h;
	    h = new TH1F(Form("h_lep_pT_nL_%s_lepnotmatched_%s",triglist.at(i-1).Data(),(key).Data()),"",nPt,0,maxPt); fOutput->AddLast(h);    h_lep_pT_nL[triglist.at(i-1)+"_lepnotmatched_"+key]   = h;

	    h2 = new TH2F(Form("h_lep_pT_eta_nT_%s_lepmatched_%s",triglist.at(i-1).Data(),(key).Data()),"",nPt,0,maxPt,nEta,-1*maxEta,maxEta); fOutput->AddLast(h2);    h_lep_pT_eta_nT[triglist.at(i-1)+"_lepmatched_"+key]   = h2;
	    h2 = new TH2F(Form("h_lep_pT_eta_nL_%s_lepmatched_%s",triglist.at(i-1).Data(),(key).Data()),"",nPt,0,maxPt,nEta,-1*maxEta,maxEta); fOutput->AddLast(h2);    h_lep_pT_eta_nL[triglist.at(i-1)+"_lepmatched_"+key]   = h2;
	    h2 = new TH2F(Form("h_lep_pT_eta_nT_%s_lepnotmatched_%s",triglist.at(i-1).Data(),(key).Data()),"",nPt,0,maxPt,nEta,-1*maxEta,maxEta); fOutput->AddLast(h2);    h_lep_pT_eta_nT[triglist.at(i-1)+"_lepnotmatched_"+key]   = h2;
	    h2 = new TH2F(Form("h_lep_pT_eta_nL_%s_lepnotmatched_%s",triglist.at(i-1).Data(),(key).Data()),"",nPt,0,maxPt,nEta,-1*maxEta,maxEta); fOutput->AddLast(h2);    h_lep_pT_eta_nL[triglist.at(i-1)+"_lepnotmatched_"+key]   = h2;


	    if(triglist.at(i-1).Contains("Trig")){
	      h = new TH1F(Form("h_lep_pT_nT_%s_evnotrig_%s",triglist.at(i-1).Data(),(key).Data()),"",nPt,0,maxPt); fOutput->AddLast(h);    h_lep_pT_nT[triglist.at(i-1)+"_evnotrig_"+key]   = h;
	      h = new TH1F(Form("h_lep_pT_nL_%s_evnotrig_%s",triglist.at(i-1).Data(),(key).Data()),"",nPt,0,maxPt); fOutput->AddLast(h);    h_lep_pT_nL[triglist.at(i-1)+"_evnotrig_"+key]   = h;
	      h = new TH1F(Form("h_lep_eta_nT_%s_evnotrig_%s",triglist.at(i-1).Data(),(key).Data()),"",nEta,-1*maxEta,maxEta); fOutput->AddLast(h);    h_lep_eta_nT[triglist.at(i-1)+"_evnotrig_"+key]   = h;
	      h = new TH1F(Form("h_lep_eta_nL_%s_evnotrig_%s",triglist.at(i-1).Data(),(key).Data()),"",nEta,-1*maxEta,maxEta); fOutput->AddLast(h);    h_lep_eta_nL[triglist.at(i-1)+"_evnotrig_"+key]   = h;
	      h2 = new TH2F(Form("h_lep_pT_eta_nT_%s_evnotrig_%s",triglist.at(i-1).Data(),(key).Data()),"",nPt,0,maxPt,nEta,-1*maxEta,maxEta); fOutput->AddLast(h2);    h_lep_pT_eta_nT[triglist.at(i-1)+"_evnotrig_"+key]   = h2;
	      h2 = new TH2F(Form("h_lep_pT_eta_nL_%s_evnotrig_%s",triglist.at(i-1).Data(),(key).Data()),"",nPt,0,maxPt,nEta,-1*maxEta,maxEta); fOutput->AddLast(h2);    h_lep_pT_eta_nL[triglist.at(i-1)+"_evnotrig_"+key]   = h2;
	    }
	    if(trigcategories[triglist.at(i-1)].Contains("E") && !(key).Contains("_E_"))continue;
	    if(trigcategories[triglist.at(i-1)].Contains("M") && !(key).Contains("_M_"))continue;

	    //cout<<"a"<<endl;

	    /**
	       for (TString bdt_name : bdt_vec) {
	       h = new TH1F(Form("h_lep_BDT_nT_%s_%s_%s",bdt_name.Data(),triglist.at(i-1).Data(),(key).Data()),"",nBDT,0,maxBDT); fOutput->AddLast(h);    h_lep_BDT_nT[bdt_name+"_"+triglist.at(i-1)+"_"+key]   = h;
	       h = new TH1F(Form("h_lep_BDT_nL_%s_%s_%s",bdt_name.Data(),triglist.at(i-1).Data(),(key).Data()),"",nBDT,0,maxBDT); fOutput->AddLast(h);    h_lep_BDT_nL[bdt_name+"_"+triglist.at(i-1)+"_"+key]   = h;	      
	       }
	       h = new TH1F(Form("h_lep_pT_nT_%s_%s",triglist.at(i-1).Data(),(key).Data()),"",nPt,0,maxPt); fOutput->AddLast(h);    h_lep_pT_nT[triglist.at(i-1)+"_"+key]   = h;
	       h = new TH1F(Form("h_lep_pT_nL_%s_%s",triglist.at(i-1).Data(),(key).Data()),"",nPt,0,maxPt); fOutput->AddLast(h);    h_lep_pT_nL[triglist.at(i-1)+"_"+key]   = h;
	       h = new TH1F(Form("h_lep_eta_nT_%s_%s",triglist.at(i-1).Data(),(key).Data()),"",nEta,-1*maxEta,maxEta); fOutput->AddLast(h);    h_lep_eta_nT[triglist.at(i-1)+"_"+key]   = h;
	       h = new TH1F(Form("h_lep_eta_nL_%s_%s",triglist.at(i-1).Data(),(key).Data()),"",nEta,-1*maxEta,maxEta); fOutput->AddLast(h);    h_lep_eta_nL[triglist.at(i-1)+"_"+key]   = h;
	       h2 = new TH2F(Form("h_lep_pT_eta_nT_%s_%s",triglist.at(i-1).Data(),(key).Data()),"",nPt,0,maxPt,nEta,-1*maxEta,maxEta); fOutput->AddLast(h2);    h_lep_pT_eta_nT[triglist.at(i-1)+"_"+key]   = h2;
	       h2 = new TH2F(Form("h_lep_pT_eta_nL_%s_%s",triglist.at(i-1).Data(),(key).Data()),"",nPt,0,maxPt,nEta,-1*maxEta,maxEta); fOutput->AddLast(h2);    h_lep_pT_eta_nL[triglist.at(i-1)+"_"+key]   = h2;
	    */
	  }
	}
	
	/**
	   h2 = new TH2F(Form("h_lep_pT_eta_nT_%s",(key).Data()),"",nPt,0,maxPt,nEta,-1*maxEta,maxEta); fOutput->AddLast(h2);    h_lep_pT_eta_nT[key]   = h2;
	   h2 = new TH2F(Form("h_lep_pT_eta_nL_%s",(key).Data()),"",nPt,0,maxPt,nEta,-1*maxEta,maxEta); fOutput->AddLast(h2);    h_lep_pT_eta_nL[key]   = h2;
	*/
	/**
	   for(it2 = categories.begin(); it2 != categories.end(); it2++){

	   //if((*it2).Contains("e_") && (key.Contains("_M_") || key.Contains("MM")))continue;
	   //if((*it2).Contains("m_") && (key.Contains("_E_") || key.Contains("EE")))continue;
	   //cout<<"b"<<endl;
	   h = new TH1F(Form("h_lep_eta_nT_%s_lepmatched_%s",(*it2).Data(),(key).Data()),"",nEta,-1*maxEta,maxEta); fOutput->AddLast(h);    h_lep_eta_nT[(*it2)+"_lepmatched_"+key]   = h;
	   h = new TH1F(Form("h_lep_eta_nL_%s_lepmatched_%s",(*it2).Data(),(key).Data()),"",nEta,-1*maxEta,maxEta); fOutput->AddLast(h);    h_lep_eta_nL[(*it2)+"_lepmatched_"+key]   = h;
	   h = new TH1F(Form("h_lep_eta_nT_%s_lepnotmatched_%s",(*it2).Data(),(key).Data()),"",nEta,-1*maxEta,maxEta); fOutput->AddLast(h);    h_lep_eta_nT[(*it2)+"_lepnotmatched_"+key]   = h;
	   h = new TH1F(Form("h_lep_eta_nL_%s_lepnotmatched_%s",(*it2).Data(),(key).Data()),"",nEta,-1*maxEta,maxEta); fOutput->AddLast(h);    h_lep_eta_nL[(*it2)+"_lepnotmatched_"+key]   = h;

	   for (TString bdt_name : bdt_vec) {
	   h = new TH1F(Form("h_lep_BDT_nT_%s_%s_lepmatched_%s",bdt_name.Data(),(*it2).Data(),(key).Data()),"",nBDT,0,maxBDT); fOutput->AddLast(h);    h_lep_BDT_nT[bdt_name+"_"+(*it2)+"_lepmatched_"+key]   = h;
	   h = new TH1F(Form("h_lep_BDT_nL_%s_%s_lepmatched_%s",bdt_name.Data(),(*it2).Data(),(key).Data()),"",nBDT,0,maxBDT); fOutput->AddLast(h);    h_lep_BDT_nL[bdt_name+"_"+(*it2)+"_lepmatched_"+key]   = h;
	   h = new TH1F(Form("h_lep_BDT_nT_%s_%s_lepnotmatched_%s",bdt_name.Data(),(*it2).Data(),(key).Data()),"",nBDT,0,maxBDT); fOutput->AddLast(h);    h_lep_BDT_nT[bdt_name+"_"+(*it2)+"_lepnotmatched_"+key]   = h;
	   h = new TH1F(Form("h_lep_BDT_nL_%s_%s_lepnotmatched_%s",bdt_name.Data(),(*it2).Data(),(key).Data()),"",nBDT,0,maxBDT); fOutput->AddLast(h);    h_lep_BDT_nL[bdt_name+"_"+(*it2)+"_lepnotmatched_"+key]   = h;
	   h = new TH1F(Form("h_lep_BDT_nT_%s_evnotrig_%s_%s",bdt_name.Data(),(*it2).Data(),(key).Data()),"",nBDT,0,maxBDT); fOutput->AddLast(h);    h_lep_BDT_nT[bdt_name+"_"+(*it2)+"_evnotrig_"+key]   = h;
	   h = new TH1F(Form("h_lep_BDT_nL_%s_evnotrig_%s_%s",bdt_name.Data(),(*it2).Data(),(key).Data()),"",nBDT,0,maxBDT); fOutput->AddLast(h);    h_lep_BDT_nL[bdt_name+"_"+(*it2)+"_evnotrig_"+key]   = h;
	   }
	   h = new TH1F(Form("h_lep_pT_nT_%s_lepmatched_%s",(*it2).Data(),(key).Data()),"",nPt,0,maxPt); fOutput->AddLast(h);    h_lep_pT_nT[(*it2)+"_lepmatched_"+key]   = h;
	   h = new TH1F(Form("h_lep_pT_nL_%s_lepmatched_%s",(*it2).Data(),(key).Data()),"",nPt,0,maxPt); fOutput->AddLast(h);    h_lep_pT_nL[(*it2)+"_lepmatched_"+key]   = h;
	   h = new TH1F(Form("h_lep_pT_nT_%s_lepnotmatched_%s",(*it2).Data(),(key).Data()),"",nPt,0,maxPt); fOutput->AddLast(h);    h_lep_pT_nT[(*it2)+"_lepnotmatched_"+key]   = h;
	   h = new TH1F(Form("h_lep_pT_nL_%s_lepnotmatched_%s",(*it2).Data(),(key).Data()),"",nPt,0,maxPt); fOutput->AddLast(h);    h_lep_pT_nL[(*it2)+"_lepnotmatched_"+key]   = h;

	   h2 = new TH2F(Form("h_lep_pT_eta_nT_%s_lepmatched_%s",(*it2).Data(),(key).Data()),"",nPt,0,maxPt,nEta,-1*maxEta,maxEta); fOutput->AddLast(h2);    h_lep_pT_eta_nT[(*it2)+"_lepmatched_"+key]   = h2;
	   h2 = new TH2F(Form("h_lep_pT_eta_nL_%s_lepmatched_%s",(*it2).Data(),(key).Data()),"",nPt,0,maxPt,nEta,-1*maxEta,maxEta); fOutput->AddLast(h2);    h_lep_pT_eta_nL[(*it2)+"_lepmatched_"+key]   = h2;
	   h2 = new TH2F(Form("h_lep_pT_eta_nT_%s_lepnotmatched_%s",(*it2).Data(),(key).Data()),"",nPt,0,maxPt,nEta,-1*maxEta,maxEta); fOutput->AddLast(h2);    h_lep_pT_eta_nT[(*it2)+"_lepnotmatched_"+key]   = h2;
	   h2 = new TH2F(Form("h_lep_pT_eta_nL_%s_lepnotmatched_%s",(*it2).Data(),(key).Data()),"",nPt,0,maxPt,nEta,-1*maxEta,maxEta); fOutput->AddLast(h2);    h_lep_pT_eta_nL[(*it2)+"_lepnotmatched_"+key]   = h2;
	 
	   h = new TH1F(Form("h_lep_pT_nT_%s_evnotrig_%s",(*it2).Data(),(key).Data()),"",nPt,0,maxPt); fOutput->AddLast(h);    h_lep_pT_nT[(*it2)+"_evnotrig_"+key]   = h;
	   h = new TH1F(Form("h_lep_pT_nL_%s_evnotrig_%s",(*it2).Data(),(key).Data()),"",nPt,0,maxPt); fOutput->AddLast(h);    h_lep_pT_nL[(*it2)+"_evnotrig_"+key]   = h;
	   h = new TH1F(Form("h_lep_eta_nT_%s_evnotrig_%s",(*it2).Data(),(key).Data()),"",nEta,-1*maxEta,maxEta); fOutput->AddLast(h);    h_lep_eta_nT[(*it2)+"_evnotrig_"+key]   = h;
	   h = new TH1F(Form("h_lep_eta_nL_%s_evnotrig_%s",(*it2).Data(),(key).Data()),"",nEta,-1*maxEta,maxEta); fOutput->AddLast(h);    h_lep_eta_nL[(*it2)+"_evnotrig_"+key]   = h;
	   h2 = new TH2F(Form("h_lep_pT_eta_nT_%s_evnotrig_%s",(*it2).Data(),(key).Data()),"",nPt,0,maxPt,nEta,-1*maxEta,maxEta); fOutput->AddLast(h2);    h_lep_pT_eta_nT[(*it2)+"_evnotrig_"+key]   = h2;
	   h2 = new TH2F(Form("h_lep_pT_eta_nL_%s_evnotrig_%s",(*it2).Data(),(key).Data()),"",nPt,0,maxPt,nEta,-1*maxEta,maxEta); fOutput->AddLast(h2);    h_lep_pT_eta_nL[(*it2)+"_evnotrig_"+key]   = h2;
	   }
	*/
	
      }
	
      }
    //cout<<"c"<<endl;
    for(it2 = metregions.begin(); it2 != metregions.end(); it2++){
      h = new TH1F(Form("h_lep_pT_nT_%s_%s",(*it1).Data(),(*it2).Data()),"",nPt,0,maxPt); fOutput->AddLast(h);    h_lep_pT_nT[*it1+"_"+*it2]   = h;
      h = new TH1F(Form("h_lep_pT_nL_%s_%s",(*it1).Data(),(*it2).Data()),"",nPt,0,maxPt); fOutput->AddLast(h);    h_lep_pT_nL[*it1+"_"+*it2]   = h;
    }
    for(it2 = mt2regions.begin(); it2 != mt2regions.end(); it2++){
      h = new TH1F(Form("h_lep_pT_nT_%s_%s",(*it1).Data(),(*it2).Data()),"",nPt,0,maxPt); fOutput->AddLast(h);    h_lep_pT_nT[*it1+"_"+*it2]   = h;
      h = new TH1F(Form("h_lep_pT_nL_%s_%s",(*it1).Data(),(*it2).Data()),"",nPt,0,maxPt); fOutput->AddLast(h);    h_lep_pT_nL[*it1+"_"+*it2]   = h;
    }
    /**
       for(it2 = metsigregions.begin(); it2 != metsigregions.end(); it2++){
       h = new TH1F(Form("h_lep_pT_nT_%s_%s",(*it1).Data(),(*it2).Data()),"",nPt,0,maxPt); fOutput->AddLast(h);    h_lep_pT_nT[*it1+"_"+*it2]   = h;
       h = new TH1F(Form("h_lep_pT_nL_%s_%s",(*it1).Data(),(*it2).Data()),"",nPt,0,maxPt); fOutput->AddLast(h);    h_lep_pT_nL[*it1+"_"+*it2]   = h;
       }
    */
    for(it2 = mllregions.begin(); it2 != mllregions.end(); it2++){
      h = new TH1F(Form("h_lep_pT_nT_%s_%s",(*it1).Data(),(*it2).Data()),"",nPt,0,maxPt); fOutput->AddLast(h);    h_lep_pT_nT[*it1+"_"+*it2]   = h;
      h = new TH1F(Form("h_lep_pT_nL_%s_%s",(*it1).Data(),(*it2).Data()),"",nPt,0,maxPt); fOutput->AddLast(h);    h_lep_pT_nL[*it1+"_"+*it2]   = h;
    }
    for(it2 = ncljregions.begin(); it2 != ncljregions.end(); it2++){
      h = new TH1F(Form("h_lep_pT_nT_%s_%s",(*it1).Data(),(*it2).Data()),"",nPt,0,maxPt); fOutput->AddLast(h);    h_lep_pT_nT[*it1+"_"+*it2]   = h;
      h = new TH1F(Form("h_lep_pT_nL_%s_%s",(*it1).Data(),(*it2).Data()),"",nPt,0,maxPt); fOutput->AddLast(h);    h_lep_pT_nL[*it1+"_"+*it2]   = h;
    }
    for(it2 = nbjregions.begin(); it2 != nbjregions.end(); it2++){
      h = new TH1F(Form("h_lep_pT_nT_%s_%s",(*it1).Data(),(*it2).Data()),"",nPt,0,maxPt); fOutput->AddLast(h);    h_lep_pT_nT[*it1+"_"+*it2]   = h;
      h = new TH1F(Form("h_lep_pT_nL_%s_%s",(*it1).Data(),(*it2).Data()),"",nPt,0,maxPt); fOutput->AddLast(h);    h_lep_pT_nL[*it1+"_"+*it2]   = h;
    }

    //}
  }

  cout<<"Number of historgrams stored "<<fOutput->GetSize()<<endl;
  
  if(doFakeEst){
    gfw = new GetFakeWeight();
    gfw->setVB(0);
    for(std::vector< TString >::iterator itLT = triggerDef.begin(); itLT != triggerDef.end(); itLT++){ 
      gfw->initializeHist((*itLT), "./MMinputfiles/FNPntuple_FCL_FEB10_final_v2.root", "./MMinputfiles/MMinput_frac_NEW_FRAC_MAY26.root", 0, "metsig", true, true, true, uncert_keys, bdt_vec);
      //gfw->initializeHist((*itLT), "./MMinputfiles/FNPntuple_2L2J_FEB10_final_v2.root", "./MMinputfiles/MMinput_frac_2L2J_FEB10.root", 0, "metsig", true, true, true, uncert_keys, bdt_vec);  
    }
  }
}

void MySusySkimAnalysis::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

bool MySusySkimAnalysis::ProcessCut(Long64_t entry)
{

  //if(*nJet20<=1)return false;
  //  if(*nBJet20_MV2c10_FixedCutBEff_77 != 0)return false;
  //if((*met_Et) <= 110)return false; 
  //if((*met_Sign) <= 10)return false; 
  //if(themt2 < 60)continue;
  
  
  return true;
  
}


Bool_t MySusySkimAnalysis::Process(Long64_t entry)
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


  // std::vector<TString> convregions;
  // std::vector<TString> HFregions;
  // std::vector<TString> LFregions;

  fillregions.clear();
  fillregions.push_back("ALL");
  cutsregions.clear();
  fillregions_wcuts.clear();

  std::vector<int> baseline_lep;
  std::vector<int> baseline_mu;
  std::vector<int> baseline_el;
  std::vector<int> baseline_mu_OR;
  std::vector<int> baseline_el_OR;
  std::vector<int> signal_mu;
  std::vector<int> signal_el;
  std::vector<int> noOR_baseline_lep;
  std::vector<int> noOR_baseline_mu;
  std::vector<int> noOR_baseline_el;

  std::vector<TString> lep1_tmatch; 
  std::vector<TString> lep2_tmatch;

  TString triggercat1;
  TString triggercat2;

  std::vector<TString> vec;

  std::vector< TString >::iterator it;
  std::vector< TString >::iterator it2;
  std::vector< TString >::iterator itend;
  std::vector< int >::iterator int_it;

  Bool_t isCR_Dib = false;
  Bool_t isCR_Top = false;
  Bool_t isDF_1bjet = false;
  Bool_t isDF_0jet = false;

  int nTT = 0;
  int nTl = 0;
  int nlT = 0;
  int nll = 0;
  //uncert.erase ( uncert.begin(), uncert.end() );    // erasing by range
  
  

  //fReader.SetLocalEntry(entry);
 
  if(nev == 0){

    //    ntupf = new TFile("ntup.root","RECREATE");
    //tree3 = ((TChain*)fReader.GetTree())->CloneTree(-1);
    //((TChain*)fReader.GetTree())->GetListOfClones()->Remove(tree3);
    //tree3->SetBranchStatus("*",0);
    //tree3->ResetBranchAddresses();
    //tree3->SetBranchAddress("jet_pt",&jetPt);
    
    
    fullpath = ((TChain*)fReader.GetTree())->GetFile()->GetName();//TString(((TChain*)(MySusySkimAnalysis::fChain))->GetFile()->GetName());

    /**
       ttree_infile  = new TFile("./data15-16_merged_processed_copy.root");
       centraltree = (TTree*)ttree_infile->Get("data15-16");//(is1516 ? "data15-16" : (is17 ? "data17" : "data18")));
       central_nentries = 0;
       ttree_newfile = new TFile("./FAKES.root","recreate");
       centraltree->SetBranchStatus("*",1);
       centraltree->SetBranchStatus("eventWeight",0);
       //centraltree->SetBranchStatus("eventWeight",0);
       newtree = centraltree->CloneTree(-1);
       //newtree->SetDirectory(0);
       //ttree_infile->Close();
       */

    //ttree_newfile = new TFile("/scratch3/eirikgr/2L2JInputs/v1.6/SUSY2/data15-16_merged_processed.root_FAKES.root");//Form("%s_FAKES.root",fullpath.ReplaceAll(".root","").Data()),"recreate");
    // bMMweight    = centraltree->Branch("eventWeight",&MMweight,"eventWeight/D");
    // bMMweight_up = centraltree->Branch("MMweight_up",&MMweight_up,"MMweight_up/D");
    // bMMweight_dw = centraltree->Branch("MMweight_dw",&MMweight_dw,"MMweight_dw/D");



    TString fname_temp;
    int start = fullpath.Last('/');
    fname_temp = fullpath(start,fullpath.Length());
    int start2 = fname_temp.Last('/');
    fname = fname_temp(start2+1,fname_temp.Length()-6);
    TObjArray *tx = fname.Tokenize("_");
    rnum = ((TObjString *)tx->At(0))->String();
    delete tx;

    //int fi = 0;
    //while(true){
    sprintf(fnameoutputhistfilename,"%s/%s_HISTOGRAMS_%lld_%lld.root",date_folder.Data(),fname.Data(),toprocess,firstentry);
    struct stat buffer;   
    if(stat (fnameoutputhistfilename, &buffer) == 0){
      cout<<"WARNING \t Histograms already exists for " << fname.Data()<< " with index "<<toprocess<<"_"<<firstentry<<" in "<< date_folder.Data()<<endl;
      Abort("Histograms exists in folder");
    }
    // }else{
    // 	break;
    // }
    //fi += 1;
    // }

    //cout<<"INFO \t Running over "<<rnum.Data()<<endl;

    if(rnum.Contains("data"))isData = true;
    else isData = false;

    if(fname.Contains("Zeejets") || fname.Contains("Zmmjets") || fname.Contains("Zttjets")){bornMass = {fReader, "bornMass"};}

    //cout<<"This is data = "<<isData<<endl;
    if(!nentries){
      nentries = ((TChain*)fReader.GetTree())->GetEntries();
    }

    
    
    if(makeNtup){
      std::vector<TString> unc_keys;
      if(doFakeEst){
	for (auto u : GetFakeWeight::all_unc) {
	  unc_keys.push_back(gfw->getUncKey(u));
	}
      }
      char ofile[200];
      sprintf(ofile,"%s/%s_MVA_%lld_%lld.root",date_folder.Data(),fname.Data(),toprocess,firstentry);
      TString pidst = string(ofile);
      HFT = new HistFitterTree(rnum,"central",pidst,"MM_CENTRAL");
      if(make2L2JNtup){
	sprintf(ofile,"%s/%s_2L2J_%lld_%lld.root",date_folder.Data(),fname.Data(),toprocess,firstentry);
	//sprintf(fnameoutputhistfilename,"%s/%s_HISTOGRAMS_%lld_%lld.root",date_folder.Data(),fname.Data(),toprocess,firstentry);
	TString pidst = string(ofile);
	twoLtwoJ = new make2L2JTree(rnum,"central",pidst,"MM_CENTRAL",unc_keys);
      }
      if(makeMYNtup){
	sprintf(ofile,"%s/%s_MY_%lld_%lld.root",date_folder.Data(),fname.Data(),toprocess,firstentry);
	//sprintf(fnameoutputhistfilename,"%s/%s_HISTOGRAMS_%lld_%lld.root",date_folder.Data(),fname.Data(),toprocess,firstentry);
	TString pidst = string(ofile);
	MY = new makeMYTree(rnum,"central",pidst,"MM_CENTRAL");
      }

      /**
      cout<<"Creating copy of ntuples..."<<endl;
      sprintf(fnameoutputntupfilename,"%s/%s_FNP_%lld_%lld.root",date_folder.Data(),fname.Data(),toprocess,firstentry);
      ntupf = new TFile(fnameoutputntupfilename,"RECREATE");


      //tree3 = ((TChain*)fReader.GetTree())->AddFriend("FNPfriend",fnameoutputntupfilename);

      tree3 = new TTree("FNPfriend","FNPfriend");
      */

      printf("Adding tree %s from file %s\n",(((TChain*)fReader.GetTree())->GetName()),(((TChain*)fReader.GetTree())->GetFile()->GetName()));
      chain = new TChain(((TChain*)fReader.GetTree())->GetName());
      chain->Add(((TChain*)fReader.GetTree())->GetFile()->GetName());

      
      sprintf(fnameoutputcopyfilename,"%s/%s_COPY_%lld_%lld.root",date_folder.Data(),fname.Data(),toprocess,firstentry);
      //TFile gh(fnameoutputcopyfilename,"RECREATE");
      ntupf = new TFile(fnameoutputcopyfilename,"RECREATE");
      newtree=chain->CloneTree(0);
      //fOutput->Add(newtree);
      //gh.Close();

     
      //exactcopy = ((TChain*)fReader.GetTree())->GetFile()->Get(((TChain*)fReader.GetTree())->GetName());

      
      /**
      cout<<"1"<<endl;
      tree3 = ((TChain*)fReader.GetTree())->CloneTree(0,"fast");
      tree3->ResetBranchAddresses();
      cout<<"2"<<endl;
      ((TChain*)fReader.GetTree())->GetListOfClones()->Remove(tree3);
      */
      /**
	 TObjArray *branchList  = ((TChain*)fReader.GetTree())->GetListOfBranches();
	 Int_t nBranch = ((TChain*)fReader.GetTree())->GetNbranches();
	 TString varnames[nBranch];
	 TBranch *var[nBranch];
	 for(int i=0;i<nBranch;i++){
	 cout<<"branch name = "<<branchList->At(i)->GetName()<<endl;
	 varnames[i] = branchList->At(i)->GetName();
	 cout<<"var"<<endl;
	 var[i] = ((TChain*)fReader.GetTree())->GetBranch(varnames[i]);//->GetLeaf(varnames[i]);
	 //var[i] = fChain->GetBranch(varnames[i])->GetLeaf(varnames[i]);
	 //value_i = var[i]->GetValue();
	 cout<<"setting adress"<<endl;
	 //h1_master[i]->Fill(value_i,weight);
	 tree3->SetBranchAddress(varnames[i],var[i]->GetAddress());
	 }
      */
     
      cout<<"Adding FNP weights + systematics"<<endl;
      for(std::vector< TString >::iterator itLT = unc_keys.begin(); itLT != unc_keys.end(); itLT++){
	b_fake_wgt[(*itLT)] = 0.0;
	cout<<"Adding barnch "<<(*(itLT)).Data()<<endl;
	if((*itLT).EqualTo("NOM")){
	  newtree->Branch("FNP_WEIGHTS",&b_fake_wgt[*(itLT)]);
	}else{
	  newtree->Branch("FNP_"+(*itLT),&b_fake_wgt[*(itLT)]);
	}
      }
      cout<<"Done"<<endl;
     
      //tree3->SetBranchStatus("*",0);
      //tree3->SetBranchStatus("GEN*",0);
      newtree->Branch("isTT",&nTT);
      newtree->Branch("isTl",&nTl);
      newtree->Branch("islT",&nlT);
      newtree->Branch("isll",&nll);



      
    }
     
  }
     
 fReader.SetLocalEntry(entry);
  
  if(nev == 0){
    cout<<"Starting with event "<<*EventNumber<<endl;

  }
  //      cout<<"bornMass = "<<(*bornMass)<<endl;
  
  //cout<<"aa"<<endl;
  //cout<<"Test = "<<( fReader.GetCurrentEntry() )<<endl;
  //if(makeNtup)
  // for(unsigned int i = 0; i<lepPt.GetSize(); i++){
  //   Bool_t med = lepMedium->at(i);
  //   Bool_t tig = lepTight->at(i);
    
  // }
  //chain->Fill();
  //cout<<"index = "<<((TChain*)fReader.GetTree())->GetTreeIndex()<<endl;

  //tree3->Branch("lepPt" , &lepPt ); 
  //tree3->Fill();

  nev += 1;
  all_nev += 1;

  if(is18 && isData && makeHistograms)h_runNumber18->Fill(*RunNumber);


  //skip non-nominal ttbar samples
  if((!isData) && (((*DatasetNumber) >= 410633 && (*DatasetNumber) <= 410637) || (!is18 && (*DatasetNumber) == 410472)))return true;

  //skip non-2.2.11 Z+jets samples
  if((!isData) && ((*DatasetNumber) >= 308092 && (*DatasetNumber) <= 308093))return true;

  
  
  
  int printevery = 1000000;
  if(all_nev%printevery == 0){
    double elapsed_secs = 0;
    double ev_per_sec = 0;
    double time_left = 0.0;
    double minutesRemainder = 0.0;
    double secondsRemainder = 0.0;
    int hours = 0;
    int minutes = 0;
    int seconds = 0;
    if(all_nev>0){
      c_end = clock();
      elapsed_secs = double(c_end - c_begin) / CLOCKS_PER_SEC;
      //cout<<"elapsed_secs = "<<elapsed_secs<<endl;
      if(elapsed_secs > 0){
	ev_per_sec = ((double)printevery)/elapsed_secs;
	if(ev_per_sec > 0){
	  time_left = (nentries-all_nev)/ev_per_sec;
	  time_left = time_left/(60.*60);
	  hours = time_left;
	  minutesRemainder = (time_left - hours) * 60;
	  minutes = minutesRemainder;
	  secondsRemainder = (minutesRemainder - minutes) * 60;
	  seconds = secondsRemainder;
	}
      }
    }
    cout<<"PID("<<mypid<<"): "<<"Processing file for runnum "<<rnum.Data()<< " with events "<< all_nev << "/"
	<< nentries << " ("<<((float)all_nev/(float)nentries)*100.<<"%)" << " (ev/sec = "<<ev_per_sec<<"). "<<"Estimated time left = "<<hours<<"h"<<minutes<<"m"<<seconds<<"s"<<endl;
    c_begin = clock();

    //printf("Reading %lld bytes in %d transactions\n",MMEST_file->GetBytesRead(), MMEST_file->GetReadCalls());
    //MMEST_tree->PrintCacheStats();
    //cout<<"Cache size = "<<MMEST_tree->GetCacheSize()<<endl;
  }

  // if(!ProcessCut(entry)){
  //   all_nev += 1;
  //   return true;
  // }
  //if(all_nev < 100000)return kTRUE;

  

  //if((*EventNumber)!=64235104)return true;;

  // if(!(std::find(mirto_ev.begin(), mirto_ev.end(), (*EventNumber)) != mirto_ev.end()))return true;

  //if(*DatasetNumber == 410442)return kTRUE;

  // returning if not connected to the period we're interested in
  if(!allEventsPass){
    if(is1516 &&  (*RandomRunNumber) > 320000)return kTRUE;
    if(is17   && ((*RandomRunNumber) < 320000 || (*RandomRunNumber) > 348000))return kTRUE;
    if(is18   &&  (*RandomRunNumber) < 348000)return kTRUE;
  }

  // if(*RunNumber == 284420 && (*EventNumber == 674058450 || *EventNumber == 670579280 || *EventNumber == 679813952 || *EventNumber == 492727254 || *EventNumber == 494103246 )){

  //   cout<<"event "<<*EventNumber<<" found before"<<endl;

  // }else{
  //   return kTRUE;
  // }
  if((*DatasetNumber >= 700320 && *DatasetNumber < 700500)){
    
    if((*DatasetNumber >= 700320 && *DatasetNumber <= 700328) && *bornMass > 120e3){
      return kTRUE;
    }
    if(makeHistograms)h_bornMass->Fill(*bornMass/1000.);
  }

  
  fnpidx = -1;
  if(evnums.size() > 0){
    if ( evnums.find(*RunNumber) == evnums.end() ) {
      return kTRUE;
    }else{
      auto it = find(evnums[*RunNumber].begin(), evnums[*RunNumber].end(), *EventNumber);
      if(it == evnums[*RunNumber].end())return kTRUE;
      else fnpidx = it - evnums[*RunNumber].begin();
    }
    
    /**
       printf("Me:    Event %llu and Run %llu passed\n",*EventNumber,*RunNumber);
       printf("Group: Event %llu and Run %llu passed\n",evnums[*RunNumber].at(fnpidx),*RunNumber);

    
       auto itt = find(done_fnpidx.begin(), done_fnpidx.end(), fnpidx);
       if(itt == done_fnpidx.end())done_fnpidx.push_back(fnpidx);
       else cout<<"FNP idx "<<fnpidx<<" was already used"<<endl;
  
  

       if(evtype.at(fnpidx).Contains("-1111") || evtype.at(fnpidx).Contains("11-11") ||
       evtype.at(fnpidx).Contains("-1313") || evtype.at(fnpidx).Contains("13-13")){
       printf("Event type is %s\n",evtype.at(fnpidx).Data());
       printf("Event %llu and Run %llu passed\n",*EventNumber,*RunNumber);
       }
    */
  
    ev_was_found += 1;
  
    for(std::map<TString,std::vector<double>>::iterator iter = BDTscores[*RunNumber].begin(); iter != BDTscores[*RunNumber].end(); ++iter){
      //printf("Getting %s at index %i = %.2f\n",iter->first.Data(),fnpidx,(iter->second).at(fnpidx));
      BDTweight[iter->first] = (iter->second).at(fnpidx);
    }
  
    /**
       printf("glob_njet[*RunNumber].at(fnpidx) = %i \n",glob_njet[*RunNumber].at(fnpidx));
       printf("glob_nbjet[*RunNumber].at(fnpidx) = %i \n",glob_nbjet[*RunNumber].at(fnpidx));
       printf("glob_METsig[*RunNumber].at(fnpidx) = %.2f \n",glob_METsig[*RunNumber].at(fnpidx));
       printf("glob_isSF[*RunNumber].at(fnpidx) = %i \n",(int)glob_isSF[*RunNumber].at(fnpidx));
    */
    if(glob_njet[*RunNumber].at(fnpidx) == 0 && glob_nbjet[*RunNumber].at(fnpidx) == 0 && glob_METsig[*RunNumber].at(fnpidx) > 8 && 
       BDTweight["BDTDeltaM100_90"]  > 0.2 && BDTweight["BDTDeltaM100_90"]  <= 0.65 && BDTweight["BDTVVDeltaM100_90"] > 0.2 &&
       /**BDTweight["BDTtopDeltaM100_90"] < 0.1 &&*/ (!glob_isSF[*RunNumber].at(fnpidx) || (glob_isSF[*RunNumber].at(fnpidx) && BDTweight["BDTothersDeltaM100_90"] < 0.01))){
      isCR_Dib = true;
    }else if(glob_njet[*RunNumber].at(fnpidx) == 1 && glob_nbjet[*RunNumber].at(fnpidx) == 1 && glob_METsig[*RunNumber].at(fnpidx) > 8/** &&
	     ((!glob_isSF[*RunNumber].at(fnpidx) && BDTweight["BDTDeltaM100_90"] > 0.5  && BDTweight["BDTDeltaM100_90"] <= 0.7) ||
		(glob_isSF[*RunNumber].at(fnpidx) && BDTweight["BDTDeltaM100_90"] > 0.7  && BDTweight["BDTDeltaM100_90"] <= 0.75))*/){
      isCR_Top = true;
      //
      if(fabs(glob_FNP_TOTAL_UP[*RunNumber].at(fnpidx)/glob_FNP_WEIGHTS[*RunNumber].at(fnpidx)) > 5 || fabs(glob_FNP_TOTAL_DOWN[*RunNumber].at(fnpidx)/glob_FNP_WEIGHTS[*RunNumber].at(fnpidx) > 5))
	printf("FNP = %.2f + %.2f - %.2f\n",glob_FNP_WEIGHTS[*RunNumber].at(fnpidx),glob_FNP_TOTAL_UP[*RunNumber].at(fnpidx),glob_FNP_TOTAL_DOWN[*RunNumber].at(fnpidx));
      FNPest += glob_FNP_WEIGHTS[*RunNumber].at(fnpidx);
      FNPest_up += glob_FNP_TOTAL_UP[*RunNumber].at(fnpidx);
      FNPest_dw += glob_FNP_TOTAL_DOWN[*RunNumber].at(fnpidx);
    }else if(glob_njet[*RunNumber].at(fnpidx) == 1 && glob_nbjet[*RunNumber].at(fnpidx) == 1 && glob_METsig[*RunNumber].at(fnpidx) > 8 && !glob_isSF[*RunNumber].at(fnpidx) && !glob_isOS[*RunNumber].at(fnpidx)){
      isDF_1bjet = true;
    }else if(glob_njet[*RunNumber].at(fnpidx) == 0 && glob_nbjet[*RunNumber].at(fnpidx) == 0 && glob_METsig[*RunNumber].at(fnpidx) > 8 && !glob_isSF[*RunNumber].at(fnpidx) && !glob_isOS[*RunNumber].at(fnpidx)){
      isDF_0jet = true;  
    }
    if(!isCR_Top && !isCR_Dib && !isDF_1bjet && !isDF_0jet)return kTRUE;
  }
  // if(*RunNumber == 284420 && (*EventNumber == 674058450 || *EventNumber == 670579280 || *EventNumber == 679813952 || *EventNumber == 492727254 || *EventNumber == 494103246 )){
    
  //cout<<"event "<<*EventNumber<<" found atd index "<<fnpidx<<endl;
  
  // }
  
  //if((*DatasetNumber != 364100) /**&& (*EventNumber) != 168843958*/)return kTRUE;
  /**
   if(nev>500000){
     if(makeNtup)delete HFT;
     if(makeNtup && make2L2JNtup)delete twoLtwoJ;
     if(makeNtup && makeMYNtup)delete MY;
     WriteToFile();
     printCutflow();
     for (unsigned int i=0; i<DSIDcheck.size(); i++){
       cout<<"DSID : "<<DSIDcheck[i]<<endl;
     }
     Abort("Enough!");
   }
  */
  // lumi15 = 3219.56
  // lumi16 = 32988.1
  // lumi17 = 44307.4
  // lumi18 = 58450.1
  yr = -1;
  float scalelumi = 1.0;
  wgt = 1.0;
  wgt_loose = 1.0;
  if(!isData){
    //if(is1516 || is17){
    scalelumi = (*RandomRunNumber) < 320000 ? 36207.65 : (((*RandomRunNumber) > 320000 && (*RandomRunNumber) < 348000) ? 44307.4 : 58450.1);
    //}else if(is18){
    //scalelumi = 58450.1;
    //}
    extrafac = 1.0;
    if(*DatasetNumber >= 700300 && *DatasetNumber <= 700500)extrafac = 1.0e3;
    
    if((*RandomRunNumber) <= 284500)yr = 2015;
    else if((*RandomRunNumber) > 284500 && (*RandomRunNumber) < 320000)yr = 2016;
    else if((*RandomRunNumber) > 320000 && (*RandomRunNumber) < 348000)yr = 2017;
    else if((*RandomRunNumber) > 348000)yr = 2018;
    else cout<<"ERROR \t Could not find year for RandomRunNumber "<<(*RandomRunNumber)<<endl; 

    // <------------------------------------- COMMENT BACK IN. NOT AVAILABLE IN 2L2J OFFICIIAL NTUPLES
    wgt = ((*genWeight) * (*eventWeight) * /**(*leptonWeight) * */ (*jvtWeight) * (*bTagWeight) * (*pileupWeight) * scalelumi * extrafac);// * ((*itLT).Contains("2L2J") ? (*globalDiLepTrigSF) : (*singleLepTrigSF)));
    wgt_loose = ((*genWeight) * (*eventWeight) * /**(*baselineleptonWeight) * */ (*jvtWeight) * (*bTagWeight) * (*pileupWeight) * scalelumi * extrafac);// * ((*itLT).Contains("2L2J") ? (*globalBaselineDiLepTrigSF) : (*singleBaselineLepTrigSF)));	   
    // <------------------------------------- COMMENT BACK IN. NOT AVAILABLE IN 2L2J OFFICIIAL NTUPLES

    //wgt = ((*genWeight) * (*eventWeight) * (*leptonWeight) * (*jvtWeight) * (*bTagWeight) * (*pileupWeight) * scalelumi); // trigger gets multiplied further down.
    //wgt_loose = wgt;
    
    //wgt       = ((*genWeight) * (*eventWeight) * (*jvtWeight) * (*bTagWeight) * (*pileupWeight) * scalelumi);
    //wgt_loose = ((*genWeight) * (*eventWeight) * (*jvtWeight) * (*bTagWeight) * (*pileupWeight) * scalelumi);
    if(makeHistograms){
      h_wgts_leptonWeight->Fill(*leptonWeight);
      h_wgts_baselineleptonWeight->Fill(*baselineleptonWeight);
      h_wgts_globalDiLepTrigSF->Fill(*globalDiLepTrigSF);
      h_wgts_singleLepTrigSF->Fill(*singleLepTrigSF);
      h_wgts_globalBaselineDiLepTrigSF->Fill(*globalBaselineDiLepTrigSF);
      h_wgts_singleBaselineLepTrigSF->Fill(*singleBaselineLepTrigSF);
    }
    
  }else{
    wgt = 1.0;
    wgt_loose = 1.0; 
    if((*RunNumber) >= 325713 && (*RunNumber) <= 340453)yr = 2017;
    else if ((*RunNumber) >= 297730 && (*RunNumber) <= 311481)yr = 2016;
    else if ((*RunNumber) >= 276262 && (*RunNumber) <= 284484)yr = 2015;
    else if ((*RunNumber) >= 348885)yr = 2018;
    else cout<<"ERROR Year for run "<<(*RunNumber)<<" is unknown"<<endl;
    
  }

  /**
  for (unsigned int i=0; i<DSIDcheck.size(); i++){
    if (DSIDcheck[i] == *DatasetNumber){
      return kTRUE;
    }
  }
  DSIDcheck.push_back(*DatasetNumber);   
  */
  if(doDupCheck && !isData){
    for (unsigned int i=0; i<DuplicationCheck.size(); i++){
      if (DuplicationCheck[i] == *EventNumber){
	cout << "DUPLICATION: " << *EventNumber << endl;
	nDup += 1;
	return kTRUE;
      }
    }
    DuplicationCheck.push_back(*EventNumber); 
  }

  
  //if(*RunNumber != 357713 || *EventNumber != 1413589294)return kTRUE;
  /**
     if(evnums.find(*RunNumber) == evnums.end() ) {
     return kTRUE;
     }else{
     if((std::count(evnums[*RunNumber].begin(), evnums[*RunNumber].end(), *EventNumber))){
     //cout<<"pass"<<endl;
     //evnums.remove(evnums.begin(), evnums.end(), *EventNumber);
     }else{
     return kTRUE;
     }
    
     }
  */
  
  /**
     if((std::count(evnums.begin(), evnums.end(), *EventNumber))){
     //cout<<"pass"<<endl;
     //evnums.remove(evnums.begin(), evnums.end(), *EventNumber);
     }else{
     return kTRUE;
     }
  */
  
  //if(nev>500000){
  //  WriteToFile();
  //  Abort("Enough!");
  //}
  std::map<TString,vector<double>> fake_wgt_vec;
  for(std::vector< TString >::iterator itLT = looseTightDef.begin(); itLT != looseTightDef.end(); itLT++){

    // Use di-lepton triggers for 2L2J analysis
    if(!allEventsPass && ((*itLT).Contains("2L2J")/** || (*itLT).Contains("ZMET")*/)  && !(*trigMatch_2LTrigOR))continue;

    // These are the jets used for the real analysis
    if((*itLT).Contains("2L2J"))n_clj = *nJet30;

    if(!isData){
      wgt *= (*itLT).Contains("2L2J") ? (*globalDiLepTrigSF) : (*singleLepTrigSF);  //<-------------- COMMENT BACK IN. NOT AVAILABLE IN 2L2J NTUPLES!
      wgt_loose *= (*itLT).Contains("2L2J") ? (*globalBaselineDiLepTrigSF) : (*singleBaselineLepTrigSF);//<-------------- COMMENT BACK IN. NOT AVAILABLE IN 2L2J NTUPLES!
    }

    
    noOR_baseline_lep.clear();
    baseline_lep.clear();
    baseline_mu.clear();
    baseline_el.clear();
    baseline_mu_OR.clear();
    baseline_el_OR.clear();
    signal_mu.clear();
    signal_el.clear();
    noOR_baseline_mu.clear();
    noOR_baseline_el.clear();

    vec.clear();

    fillregions.clear();
    fillregions.push_back("ALL");
    cutsregions.clear();
    fillregions_wcuts.clear();

    LT_postfix = (*itLT);

    //cout<<"LT_postfix = "<<LT_postfix.Data()<<", nev = "<<nev<<endl;
  
    for(unsigned int i = 0; i<lepPt.GetSize(); i++){

      //printf("lep%i) me: %.2f group 1: %.2f group 2: %.2f\n",i,lepPt[i],lepton1pT[*RunNumber].at(fnpidx),lepton2pT[*RunNumber].at(fnpidx));
    
      
      //cout<<"i = "<<i<<endl;
      //cout<<"lepPt["<<i<<"] = "<<lepPt[i]<<endl;
      if(lepPassOR->at(i) && isL(i))baseline_lep.push_back(i);
      if(isL(i))noOR_baseline_lep.push_back(i);

      if(lepPassOR->at(i) && lepFlavor[i] == 2 && isL(i))baseline_mu.push_back(i);
      else if(lepPassOR->at(i) && lepFlavor[i] == 1 && isL(i))baseline_el.push_back(i);

      if(lepPassOR->at(i) && lepFlavor[i] == 2 && isL(i))baseline_mu_OR.push_back(i);
      else if(lepPassOR->at(i) && lepFlavor[i] == 1 && isL(i))baseline_el_OR.push_back(i);

      if(lepPassOR->at(i) && lepFlavor[i] == 2 && isT(i))signal_mu.push_back(i);
      else if(lepPassOR->at(i) && lepFlavor[i] == 1 && isT(i))signal_el.push_back(i);

      if(lepFlavor[i] == 2 && isL(i))noOR_baseline_mu.push_back(i);
      else if(lepFlavor[i] == 1 && isL(i))noOR_baseline_el.push_back(i);

    }

    //if(baseline_lep.size() != 3)return true;
    //if(lepCharge[baseline_lep.at(0)]*lepCharge[baseline_lep.at(1)]<0)return true;
    int n_bljpreOR = 0;

    n_bjet = 0;
    n_sgj = 0;
    int n_cbj = 0;
    int n_cbj30 = 0;
    int n_cbj40 = 0;
    int n_cbj50 = 0;
    int n_cbj60 = 0;

    n_bjet85 = 0;
    n_bjet77 = 0;

    int n_blj = 0;
    int n_blj30 = 0;
    int n_blj40 = 0;
    int n_blj50 = 0;
    int n_blj60 = 0;

    n_clj = 0;
    n_clj30 = 0;
    int n_clj40 = 0;
    int n_clj50 = 0;
    int n_clj60 = 0;

    int n_fwj = 0;
    int n_fwj30 = 0;
    int n_fwj40 = 0;
    int n_fwj50 = 0;
    int n_fwj60 = 0;

    //cout<<"jetPt.GetSize() = "<<jetPt.GetSize()<<endl;

    for(unsigned int i=0; i<jetPt.GetSize(); i++){

      // cout<<"jetpt["<<i<<"] = "<<jetPt[i]<<endl;

      if(isSGjet(i)){
	n_sgj += 1;
	//}

	//if(fabs(jetEta[i])<2.4){
      }
      if(!isBLjet(i))continue;
      n_bljpreOR += 1;
      n_blj += 1;
      if(jetPt[i] > 30)n_blj30 += 1;
      if(jetPt[i] > 40)n_blj40 += 1;
      if(jetPt[i] > 50)n_blj50 += 1;
      if(jetPt[i] > 60)n_blj60 += 1;
      if(isSGjet(i)){
	// increment counters
	if(jetMV2c10[i]>0.11)n_bjet85 += 1;
	if(jetMV2c10[i]>0.64)n_bjet77 += 1;
      }
      if(isSGjet(i) && !isBjet(i)){
	n_clj += 1;
	if(jetPt[i] > 30)n_clj30 += 1;
	if(jetPt[i] > 40)n_clj40 += 1;
	if(jetPt[i] > 50)n_clj50 += 1;
	if(jetPt[i] > 60)n_clj60 += 1;
      }else if(isSGjet(i) && fabs(jetEta[i])<2.4 && isBjet(i)){
	n_cbj += 1;
	if(jetPt[i] > 30)n_cbj30 += 1;
	if(jetPt[i] > 40)n_cbj40 += 1;
	if(jetPt[i] > 50)n_cbj50 += 1;
	if(jetPt[i] > 60)n_cbj60 += 1;
      }else if(fabs(jetEta[i])>=2.4){
	n_fwj += 1;
	if(jetPt[i] > 30)n_fwj30 += 1;
	if(jetPt[i] > 40)n_fwj40 += 1;
	if(jetPt[i] > 50)n_fwj50 += 1;
	if(jetPt[i] > 60)n_fwj60 += 1;
      }
      // }

      if(isBjet(i)){
	n_bjet += 1;
      }
    }
    

    
    
    //cout<<"a"<<endl;
    //  Double_t Met_x = (*met_Et)*cos(*met_Phi);
    // Double_t Met_y = (*met_Et)*sin(*met_Phi);
    met_tlv.SetPtEtaPhiM((*met_Et),0.0,(*met_Phi),0);
    //met_tlv2.SetPtEtaPhiM((*met_Et),0.0,(*met_Phi),0);

    std::vector<int> LF_tag_el;
    std::vector<int> LF_probe_el;
    // Find the probe (i.e. lepton from W decay)
    for(int_it = baseline_el.begin(); int_it != baseline_el.end(); int_it++){
      TLorentzVector el;
      if(!lepTight->at(*int_it))continue; // <------------------------------------- COMMENT BACK IN. NOT AVAILABLE IN 2L2J OFFICIIAL NTUPLES
      if(!isT((*int_it)))continue;

      el.SetPtEtaPhiM(lepPt[*int_it],lepEta[*int_it],lepPhi[*int_it],lepM[*int_it]);
      float mindr_jet = 9999;
      for(unsigned int i=0; i<jetPt.GetSize(); i++){
	if(!isBLjet(i))continue;
	TLorentzVector jet;
	jet.SetPtEtaPhiM(jetPt[i],jetEta[i],jetPhi[i],jetM[i]);
	// check if this muon can be used as probe
	if(mindr_jet>el.DeltaR(jet)){
	  mindr_jet = el.DeltaR(jet);
	}
	if(mindr_jet < 0.4)break;
      }
      if(mindr_jet < 0.4)continue;
      LF_tag_el.push_back((*int_it));
    }
    // Check if there are any other leptons (must be fakes!)
    for(int_it = baseline_el.begin(); int_it != baseline_el.end(); int_it++){
      if(!(std::find(LF_tag_el.begin(), LF_tag_el.end(), (*int_it)) != LF_tag_el.end())) {
	LF_probe_el.push_back((*int_it));
      }
    }
    if(LF_tag_el.size() == 1 && LF_probe_el.size() == 1)isLFRegion(LF_probe_el,LF_tag_el,(*itLT));


    int n_bjet_HF = 0;
    std::vector<int> HF_tag_mu;
    std::vector<int> HF_tag_mu_tight;
    std::vector<int> HF_probe_lep;
    std::vector<int> HF_probe_lep_tight;
    HF_probe_lep.clear();
    HF_tag_mu.clear();
    int closest_to_bjet = -1;
    float mindr_bjet = 9999;
    for(int_it = noOR_baseline_mu.begin(); int_it != noOR_baseline_mu.end(); int_it++){
      //bool isTag = false;
      float drbjet = 9999;
      TLorentzVector mu;
      float mindr = 9999;
      mu.SetPtEtaPhiM(lepPt[*int_it],lepEta[*int_it],lepPhi[*int_it],lepM[*int_it]);
      // loop over all jets
      n_bjet_HF = 0;
      for(unsigned int i=0; i<jetPt.GetSize(); i++){
	if(!isBLjet(i))continue;
	TLorentzVector jet;
	jet.SetPtEtaPhiM(jetPt[i],jetEta[i],jetPhi[i],jetM[i]);
	// check if this muon can be used as probe
	if(mindr_bjet>mu.DeltaR(jet) && isBjet(i)){
	  mindr_bjet = mu.DeltaR(jet);
	  closest_to_bjet = (*int_it);
	  //isTag = true;
	}
	if(isBjet(i))n_bjet_HF += 1;
	
	// require jets to pass OR if considered in the search for probe
	if(!jetPassOR->at(i))continue; // <------------------------------------- COMMENT BACK IN. NOT AVAILABLE IN 2L2J OFFICIIAL NTUPLES
	
	if(mindr>mu.DeltaR(jet)){
	  mindr = mu.DeltaR(jet);
	}
	// check if this muons can be used as tag
	
	//printf("DeltaR(mu,jet) = %.6f\n",mu.DeltaR(jet));
	/*if(mu.DeltaR(jet)<0.4)isTag = true;*/
	
      }
      if(mindr > 0.4)HF_probe_lep.push_back((*int_it));
    }
    if(closest_to_bjet >= 0 && mindr_bjet <= 0.3){
      HF_tag_mu.push_back(closest_to_bjet);
      HF_probe_lep.erase(remove(HF_probe_lep.begin(), HF_probe_lep.end(), closest_to_bjet), HF_probe_lep.end());
      if(makeHistograms)h_HFreg_mindr_jet_probe->Fill(mindr_bjet,wgt);
      if(isT(closest_to_bjet)){
	HF_tag_mu_tight.push_back(closest_to_bjet);
	HF_probe_lep_tight.erase(remove(HF_probe_lep_tight.begin(), HF_probe_lep_tight.end(), closest_to_bjet), HF_probe_lep_tight.end());
      }
    }

    //baseline_el
    // trying with electrons removed by OR
    //for(int_it = baseline_el.begin(); int_it != baseline_el.end(); int_it++){
    for(int_it = noOR_baseline_el.begin(); int_it != noOR_baseline_el.end(); int_it++){
      float mindr = 9999;
      
      TLorentzVector el;
      el.SetPtEtaPhiM(lepPt[*int_it],lepEta[*int_it],lepPhi[*int_it],lepM[*int_it]);
      // loop over all jets (after OR)
      for(unsigned int i=0; i<jetPt.GetSize(); i++){
	if(!jetPassOR->at(i))continue; // <------------------------------------- COMMENT BACK IN. NOT AVAILABLE IN 2L2J OFFICIIAL NTUPLES
	if(!isBLjet(i))continue;
	TLorentzVector jet;
	jet.SetPtEtaPhiM(jetPt[i],jetEta[i],jetPhi[i],jetM[i]);
	
	// check if this electron can be used as probe
	if(mindr>el.DeltaR(jet))mindr = el.DeltaR(jet);
      }
      //cout<<"Try: Electron is probe with flavor "<<lepFlavor[*int_it]<<" and pT = "<<lepPt[*int_it]<<" and mindr = "<<mindr<<endl; 
      if(mindr > 0.4){
	HF_probe_lep.push_back((*int_it));
	//cout<<"Succ: Electron is probe with flavor "<<lepFlavor[*int_it]<<" and pT = "<<lepPt[*int_it]<<" and mindr = "<<mindr<<endl; 
      }
    }


    
    //if(computeRates){
    //if(noOR_baseline_mu.size() == 2 && noOR_baseline_el.size() == 1)isConvRegion(noOR_baseline_mu.at(0),noOR_baseline_mu.at(1),noOR_baseline_el.at(0),n_bjet);
    if(HF_tag_mu.size() == 1 && HF_probe_lep.size() == 1){
      for(int_it = HF_probe_lep.begin(); int_it != HF_probe_lep.end(); int_it++){// <------------------------------------- COMMENT BACK IN. NOT AVAILABLE IN 2L2J OFFICIIAL NTUPLES
      	if(HF_tag_mu.at(0) == (*int_it))continue;
      	TLorentzVector el;
      	el.SetPtEtaPhiM(lepPt[*int_it],lepEta[*int_it],lepPhi[*int_it],lepM[*int_it]);
      	//loop over all jets (after OR)
      	for(unsigned int i=0; i<jetPt.GetSize(); i++){
      	  TLorentzVector jet;
      	  jet.SetPtEtaPhiM(jetPt[i],jetEta[i],jetPhi[i],jetM[i]);
      	  if(!jetPassOR->at(i)){ 
      	    if(lepFlavor[*int_it] == 1){
      	      if(makeHistograms){
		h_probe_el_hf_deltaR_ORjet->Fill(jet.DeltaR(el));
		if(jetMV2c10[i]>0.64)h_probe_el_hf_deltaR_ORbjet->Fill(jet.DeltaR(el));
	      }
      	    }else{
      	      if(makeHistograms){
		h_probe_mu_hf_deltaR_ORjet->Fill(jet.DeltaR(el));
		if(jetMV2c10[i]>0.64)h_probe_mu_hf_deltaR_ORbjet->Fill(jet.DeltaR(el));
	      }
      	    }
      	  } 
      	}
      }// <------------------------------------- COMMENT BACK IN. NOT AVAILABLE IN 2L2J OFFICIIAL NTUPLES
      isHFRegion(HF_probe_lep,HF_tag_mu,n_bjet_HF,(*itLT));
    }
    
    if(HF_tag_mu_tight.size() == 1 && HF_probe_lep_tight.size() == 1)isHFRegion(HF_probe_lep_tight,HF_tag_mu_tight,n_bjet_HF,(*itLT),true);

    /**
    auto it  = find(cutsregions.begin(), cutsregions.end(), "M_FAKE2L21");
    auto itt = find(cutsregions.begin(), cutsregions.end(), "E_FAKE2L21");
    if(it == cutsregions.end() && itt == cutsregions.end())return kTRUE;
    //}
    */
    if(makeHistograms){
      h_HFreg_nprobes->Fill(HF_probe_lep.size(),wgt);
      h_HFreg_ntags_nprobes->Fill(HF_tag_mu.size(),HF_probe_lep.size(),wgt);
    }
    
    //    if(noOR_baseline_mu.size() == 2 && noOR_baseline_el.size() == 1)isConvRegion(noOR_baseline_mu.at(0),noOR_baseline_mu.at(1),noOR_baseline_el.at(0),n_bjet,(*itLT));

    //continue;
    
    //if(!doAnalysis && !doFakeEst)continue;

    int idx0;
    int idx1;
    int idx2;
    //return true;
    /**
       if((noOR_baseline_lep.size() >= 3 ? isT(noOR_baseline_lep.at(0)) : true) && ((!doThreeLep && (noOR_baseline_lep.size() == 2)) || (doThreeLep && (noOR_baseline_lep.size() == 3)))){
       idx0 = noOR_baseline_lep.size() >= 3 ? noOR_baseline_lep.at(0) : -1;
       idx1 = noOR_baseline_lep.size() >= 3 ? noOR_baseline_lep.at(1) : noOR_baseline_lep.at(0);
       idx2 = noOR_baseline_lep.size() >= 3 ? noOR_baseline_lep.at(2) : noOR_baseline_lep.at(1);
       fillregions = getFillVector(idx1,idx2,idx0);
       isREALRegion(idx0,idx1,idx2,LT_postfix,true);
       }
    */
      
    
    // cout<<"n BL = "<<(baseline_mu.size()+baseline_el.size())<<endl;
    // for(unsigned int i = 0; i<lepPt.GetSize(); i++){
    //   cout<<"lepPt["<<i<<"] = "<<lepPt[i]<<endl;
    //   cout<<"lepEta["<<i<<"] = "<<lepEta[i]<<endl;
    //   cout<<"lepPassOR["<<i<<"] = "<<lepPassOR->at(i)<<endl;
    //   cout<<"lepFlavor["<<i<<"] = "<<lepFlavor[i]<<endl;
    //   cout<<"fabs(lepZ0SinTheta["<<i<<"]) = " << fabs(lepZ0SinTheta[i])<<endl;
    // }
    
    
    if((baseline_mu.size()+baseline_el.size())<2){
      //printf("n(mu) = %i, n(el) = %i\n",baseline_mu.size(),baseline_el.size());
      //myfile << *RunNumber << " " << *EventNumber << "\n";
      /**
	 cout<<" baseline_lep.size() = "<< baseline_lep.size()<<endl;
	 cout<<" baseline_mu.size() = "<< baseline_mu.size()<<endl;
	 cout<<" baseline_el.size() = "<< baseline_el.size()<<endl;
	 cout<<" signal_mu.size() = "<< signal_mu.size()<<endl;
	 cout<<" signal_el.size() = "<< signal_el.size()<<endl;
	 printInfo();
      */
      continue;
    }


    //cout<<" (baseline_mu.size()+baseline_el.size()) = "<<(baseline_mu.size()+baseline_el.size())<<endl;
    

    // if RJR analysis - let all events pass!
    if(!doRJR && !allEventsPass){
      if(doTwoLep && ((baseline_mu.size()+baseline_el.size()) != 2)/** && (signal_mu.size()+signal_el.size()) == 2)*/)continue;
      else if(doThreeLep && ((baseline_mu.size()+baseline_el.size()) != 3)/** && (signal_mu.size()+signal_el.size()) == 2)*/)continue;
    }else if(doRJR){
      doThreeLep = ((baseline_mu.size()+baseline_el.size()) == 3) ? true : false;
    }

    
    if(cutflowstr.EqualTo("CRTOP") || cutflowstr.EqualTo("CRWW")){cutflow["== 2 baseline"] += 1;cutflow_wgt["== 2 baseline"] += wgt;}

    
    //if(lepPt.GetSize() < 2)return true;
    // lep1Type = TruthClassification(lepType[idx1], lepOrigin[idx1], lepEgMotherType[idx1], lepEgMotherOrigin[idx1], lepEgMotherPdgId[idx1], lepCharge[idx1]);
    // lep2Type = TruthClassification(lepType[idx2], lepOrigin[idx2], lepEgMotherType[idx2], lepEgMotherOrigin[idx2], lepEgMotherPdgId[idx2], lepCharge[idx2]);
    //  cout<<"-----"<<endl;
    idx0 = doThreeLep ? baseline_lep.at(0) : -1;
    idx1 = doThreeLep ? baseline_lep.at(1) : baseline_lep.at(0);
    idx2 = doThreeLep ? baseline_lep.at(2) : baseline_lep.at(1);

    
    // if(!allEventsPass && !(*itLT).Contains("2L2J")){
    //   int lepistriggermatched = 0;
    //   for(unsigned int i=0; i<trigvallep.size();i++){
    // 	if(triglist.at(i).Contains("Trig"))continue;
    // 	//if(*trigval.at(i)){
    // 	if(trigvallep.at(i)->at(idx1))lepistriggermatched += 1;
    // 	if(trigvallep.at(i)->at(idx2))lepistriggermatched += 1;
    // 	if(lepistriggermatched)break;
    //   }
    //   if(!lepistriggermatched)continue;
    // }
    if(!allEventsPass && (*itLT).Contains("2L2J")){
      if(doThreeLep && !(lepPt[idx0] > 25 && lepPt[idx1] > 25 && lepPt[idx2] > 20))continue;
      else if(!(lepPt[idx1] > 25 && lepPt[idx2] > 25))continue;
    }else if(!allEventsPass && (*itLT).Contains("ZMET")){
      //if(lepFlavor[idx1] != lepFlavor[idx2])continue;
      int n_match_1 = checkTriggerMatch(idx1,yr).size();
      int n_match_2 = checkTriggerMatch(idx2,yr).size();
      // cout<<"checkTriggerMatch("<<idx1<<","<<yr<<").size() = "<<n_match_1<<endl;
      // cout<<"checkTriggerMatch("<<idx2<<","<<yr<<").size() = "<<n_match_2<<endl;
      if((n_match_1 <= 0 && n_match_2 <= 0))continue;
      if((lepFlavor[idx1] == 2 && lepFlavor[idx2] == 2) && !(lepPt[idx1] > 27 && lepPt[idx2] > 20))continue;
      else if((lepFlavor[idx1] == 1 && lepFlavor[idx2] == 1) && !(lepPt[idx1] > 25 && lepPt[idx2] > 25))continue;
      // cout<<"lepPt["<<idx1<<"]  = "<<lepPt[idx1] <<endl;
      // cout<<"lepPt["<<idx2<<"]  = "<<lepPt[idx2] <<endl;
    }else if(!allEventsPass /**&& !(*itLT).Contains("2L2J")*/){
      if(lepPt[idx1] < 27 || lepPt[idx2] < 9)continue;
      if((checkTriggerMatch(idx1,yr).size() <= 0 && checkTriggerMatch(idx2,yr).size() <= 0))continue;
      if(checkTriggerMatch(idx1,yr).size() <= 0 && checkTriggerMatch(idx2,yr).size() > 0 && lepPt[idx2] <= 27000)continue;
    }

    

    triggercat1 = "";
    triggercat2 = "";
    if(!(*itLT).Contains("2L2J")){
      // Find triger match for the leptons
      lep1_tmatch = getTriggerCat(idx1,yr);
      lep2_tmatch = getTriggerCat(idx2,yr);
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
    }
    
    //printf("cat 1 = %s, cat2 = %s\n",triggercat1.Data(),triggercat2.Data());
      
    vector<int> myIdx;
    myIdx.clear();
    myIdx.push_back(idx0);
    myIdx.push_back(idx1);
    myIdx.push_back(idx2);
    bool lep0isTight = doThreeLep ? isT(idx0) : true;
    bool lep1isTight = isT(idx1);
    bool lep2isTight = isT(idx2);


    if(!isData){
      //<------------------------------------- COMMENT BACK IN. NOT AVAILABLE IN 2L2J OFFICIIAL NTUPLES
      if(doThreeLep){
	wgt *= lep0isTight ? lepRecoSF[idx0] : lepBLRecoSF[idx0];
	wgt_loose *= lepBLRecoSF[idx0];
      }
    
      wgt *= lep1isTight ? lepRecoSF[idx1] : 1.0;//lepBLRecoSF[idx1];
      wgt *= lep2isTight ? lepRecoSF[idx2] : 1.0;//lepBLRecoSF[idx2];
    
      wgt_loose *= lepBLRecoSF[idx1];
      wgt_loose *= lepBLRecoSF[idx2];
    
      if(!(*itLT).Contains("2L2J")){
      
	if(doThreeLep){
	  wgt *= lep0isTight ? lepTrigSF[idx0] : 1.0;
	  wgt_loose *= lepBLTrigSF[idx0];
	}     
	wgt *= lep1isTight ? lepTrigSF[idx1] : 1.0;
	wgt *= lep2isTight ? lepTrigSF[idx2] : 1.0;
      
	wgt_loose *= lepBLTrigSF[idx1];
	wgt_loose *= lepBLTrigSF[idx2];
      
      }
      if(makeHistograms){
	h_wgts_lepBLRecoSF_lep1->Fill(lepBLRecoSF[idx1]);
	h_wgts_lepBLRecoSF_lep2->Fill(lepBLRecoSF[idx2]);
	h_wgts_lepBLTrigSF_lep1->Fill(lepBLTrigSF[idx1]);
	h_wgts_lepBLTrigSF_lep2->Fill(lepBLTrigSF[idx2]);

	h_wgts_lepRecoSF_lep1->Fill(lepRecoSF[idx1]);
	h_wgts_lepRecoSF_lep2->Fill(lepRecoSF[idx2]);
	h_wgts_lepTrigSF_lep1->Fill(lepTrigSF[idx1]);
	h_wgts_lepTrigSF_lep2->Fill(lepTrigSF[idx2]);

	h_wgts_wgt->Fill(wgt);
	h_wgts_wgtloose->Fill(wgt_loose);
      }
    }
    
    //<------------------------------------- COMMENT BACK IN. NOT AVAILABLE IN 2L2J OFFICIIAL NTUPLES
    
    
    //wgt_loose = wgt;

    // printf("Wgt    = %.2f\n",wgt);
    // printf("Wgt BL = %.2f\n",wgt_loose);
      
    TLorentzVector lep0;
    TLorentzVector lep1;
    TLorentzVector lep2;
    vector<TLorentzVector> myLeptons;
    myLeptons.clear();
    if(doThreeLep)lep0.SetPtEtaPhiM(lepPt[idx0],lepEta[idx0],lepPhi[idx0],lepM[idx0]);
    lep1.SetPtEtaPhiM(lepPt[idx1],lepEta[idx1],lepPhi[idx1],lepM[idx1]);
    lep2.SetPtEtaPhiM(lepPt[idx2],lepEta[idx2],lepPhi[idx2],lepM[idx2]);
    myLeptons.push_back(lep0);
    myLeptons.push_back(lep1);
    myLeptons.push_back(lep2);

    Double_t min_dphi_lep_met = 9999;
    for(auto l : myLeptons){
      Double_t dphi = fabs(l.DeltaPhi(met_tlv));
      if(dphi < min_dphi_lep_met){
	min_dphi_lep_met = dphi;
      }
    }
    if(makeHistograms)h_min_dphi_lep_met->Fill(min_dphi_lep_met);
    
    metrel_Et = (min_dphi_lep_met < M_PI/2.0) ? (*met_Et)*sin(min_dphi_lep_met) : (*met_Et);

    
    std::pair <TString,TString> ret = getFillString(idx1, idx2, idx0);

    std::vector< TString > vec;    
    if(lepCharge[idx1]*lepCharge[idx2] > 0 && lepECIDS[idx1] == 1 && lepECIDS[idx2] == 1){
      if(n_bjet == 0)vec.push_back("2L40");
      if(n_bjet == 1)vec.push_back("2L41");
      if(n_bjet >= 1)vec.push_back("2L42");
    }
    if(lepCharge[idx1]*lepCharge[idx2] > 0){
      if(n_bjet == 0)vec.push_back("2L43");
      if(n_bjet == 1)vec.push_back("2L44");
      if(n_bjet >= 1)vec.push_back("2L45");
    }	
    if(vec.size() > 0){
      fillLandT(idx1,vec,(*itLT),true);
      fillLandT(idx2,vec,(*itLT),true);

      /**
	 fillLandT(idx1,vec,true,*trigMatch_2LTrigOR,"2LTrig");
	 fillLandT(idx2,vec,true,*trigMatch_2LTrigOR,"2LTrig");

	 fillLandT(idx1,vec,true,*trigMatch_1L2LTrigOR,"1L2LTrig");
	 fillLandT(idx2,vec,true,*trigMatch_1L2LTrigOR,"1L2LTrig");

	 fillLandT(idx1,vec,true,*trigMatch_1LTrigOR,"1LTrig");
	 fillLandT(idx2,vec,true,*trigMatch_1LTrigOR,"1LTrig");
      */
    }
    vec.clear();
  
    
    // int nTT = 0;
    // int nTl = 0;
    // int nlT = 0;
    // int nll = 0;
    // //int nBL_afterOR = 0; 
    
    // classifyEvent(idx1, idx2, nTT, nTl, nlT, nll);
    //myfile <<(*RunNumber)<<" "<<(*EventNumber)<<" "<<(baseline_mu.size()+baseline_el.size())<<" "<< 
    //  (signal_mu.size()+signal_el.size()) <<" "<< lepPassOR->at(idx1)<<" "<< lepPassOR->at(idx2)<<" "<<ret.first<<" "<<ret.second<< " " << 
    //  nTT << " " << nTl << " " << nlT << " " << nll << " " << (baseline_mu_OR.size() + baseline_el_OR.size()) << endl;

    nlep_base_OR = (baseline_mu_OR.size() + baseline_el_OR.size());
    
    
    
    if((cutflowstr.EqualTo("CRTOP") || cutflowstr.EqualTo("CRWW")) && (ret.first).EqualTo("EM")){cutflow["DF"] += 1;cutflow_wgt["DF"] += wgt;}
    else if((cutflowstr.EqualTo("CRTOP") || cutflowstr.EqualTo("CRWW"))) continue;
    if((cutflowstr.EqualTo("CRTOP") || cutflowstr.EqualTo("CRWW")) && (ret.second).EqualTo("OS")){cutflow["OS"] += 1;cutflow_wgt["OS"] += wgt;}
    else if((cutflowstr.EqualTo("CRTOP") || cutflowstr.EqualTo("CRWW"))) continue;
    if((cutflowstr.EqualTo("CRTOP") || cutflowstr.EqualTo("CRWW")) && (signal_mu.size()+signal_el.size()) == 2){cutflow["== 2 signal"] += 1;cutflow_wgt["== 2 signal"] += wgt;}
    else if((cutflowstr.EqualTo("CRTOP") || cutflowstr.EqualTo("CRWW"))) continue;
    if((cutflowstr.EqualTo("CRTOP") || cutflowstr.EqualTo("CRWW")) && (lep1.Pt() > 25 && lep2.Pt() > 9)){cutflow["pT > 25"] += 1;cutflow_wgt["pT > 25"] += wgt;}
    else if((cutflowstr.EqualTo("CRTOP") || cutflowstr.EqualTo("CRWW"))) continue;
    if((cutflowstr.EqualTo("CRTOP") || cutflowstr.EqualTo("CRWW")) && n_clj == 0){cutflow["n_clj == 0"] += 1;cutflow_wgt["n_clj == 0"] += wgt;}
    else if((cutflowstr.EqualTo("CRTOP") || cutflowstr.EqualTo("CRWW"))) continue;

    if(cutflowstr.EqualTo("CRWW") && n_bjet == 0){cutflow["n_bjet == 0"] += 1;cutflow_wgt["n_bjet == 0"] += wgt;}
    else if(cutflowstr.EqualTo("CRWW")) continue;

    if(cutflowstr.EqualTo("CRTOP") && n_bjet == 1){cutflow["n_bjet == 1"] += 1;cutflow_wgt["n_bjet == 1"] += wgt;}
    else if(cutflowstr.EqualTo("CRTOP")) continue;
    


    //cout<<"4"<<endl;
    // find the lepton pair closest to Z-mass
    if(doThreeLep){
      double diff = 10000000000.0;
      int Zlep1 = -99;
      int Zlep2 = -99;
      int Wlep1 = -999;
      double Zmass = -999.0;
      bool foundSFOS = false;
      mtw = -1;
      themll = -1;
      for(unsigned int i=0; i<myLeptons.size(); i++)
	{
	  for(unsigned int j=i+1; j<myLeptons.size(); j++)
	    {
	      //Opposite-Sign
	      if(lepCharge[myIdx.at(i)]*lepCharge[myIdx.at(j)]<0)
		{
		  //Same-Flavor
		  if(abs(lepCharge[myIdx.at(i)])==abs(lepCharge[myIdx.at(j)]))
		    {
		      double mass = (myLeptons.at(i)+myLeptons.at(j)).M();
		      double massdiff = fabs(mass-91187.6);
		      if(massdiff<diff)
			{
			  diff=massdiff;
			  Zmass=mass;
			  Zlep1 = myIdx.at(i);
			  Zlep2 = myIdx.at(j);
			  foundSFOS = true;
			}
		    }
		}
	    }

	}
      if(!doRJR && !foundSFOS && !allEventsPass)continue;
      else if(doRJR && foundSFOS){//return true;
	if((Zlep1==myIdx.at(0) && Zlep2==myIdx.at(1)) || (Zlep1==myIdx.at(1) && Zlep2==myIdx.at(0)) ) Wlep1=myIdx.at(2);
	else if((Zlep1==myIdx.at(0) && Zlep2==myIdx.at(2)) || (Zlep1==myIdx.at(2) && Zlep2==myIdx.at(0)) ) Wlep1=myIdx.at(1);
	else if((Zlep1==myIdx.at(1) && Zlep2==myIdx.at(2)) || (Zlep1==myIdx.at(2) && Zlep2==myIdx.at(1)) ) Wlep1=myIdx.at(0);
	TLorentzVector lepW;
	lepW.SetPtEtaPhiM(lepPt[Wlep1],lepEta[Wlep1],lepPhi[Wlep1],lepM[Wlep1]);
	double wlepMetphi = lepW.DeltaPhi(met_tlv);
	mtw = sqrt(2*lepW.Pt()*met_tlv.Pt()*(1-cos(wlepMetphi)));
	themll = Zmass;
      }
    }else{
      themll = (lep1+lep2).M();
      mtw = -1;
    }
    ComputeMT2 mycalc = ComputeMT2(lep1,lep2,met_tlv,0.,0.);
    themt2 = mycalc.Compute();
    //ComputeMT2 mycalc2 = ComputeMT2(lep1,lep2,met_tlv2,0.,0.);
    //double themt2_2 = mycalc2.Compute();
    truthENUMIFF lep1Type = DATA; // default for data
    truthENUMIFF lep2Type = DATA; // default for data
    truthENUMIFF lep0Type = DATA; // default for data

    
    if((cutflowstr.EqualTo("CRTOP") || cutflowstr.EqualTo("CRWW")) && themll > 11){cutflow["mll > 11"] += 1;cutflow_wgt["mll > 11"] += wgt;}
    else if((cutflowstr.EqualTo("CRTOP") || cutflowstr.EqualTo("CRWW"))) continue;


    // if(cutflowstr.EqualTo("CRTOP") && (*met_Et) > 110){cutflow["met_Et > 110"] += 1;cutflow_wgt["met_Et > 110"] += wgt;}
    // else if(cutflowstr.EqualTo("CRTOP")) continue;
    if(cutflowstr.EqualTo("CRTOP") && (*met_Sign) > 8){cutflow["met_Sign > 8"] += 1;cutflow_wgt["met_Sign > 8"] += wgt;}
    else if(cutflowstr.EqualTo("CRTOP")) continue;
    // if(cutflowstr.EqualTo("CRTOP") && themt2 > 80){cutflow["themt2 > 80"] += 1;cutflow_wgt["themt2 > 80"] += wgt;}
    // else if(cutflowstr.EqualTo("CRTOP")) continue;

    /**
       if(cutflowstr.EqualTo("CRWW") && (*met_Et) <= 60 && (*met_Et) >= 100){cutflow["60 <= met_Et <= 100"] += 1;cutflow_wgt["60 <= met_Et <= 100"] += wgt;}
       else if(cutflowstr.EqualTo("CRWW")) continue;
       if(cutflowstr.EqualTo("CRWW") && (*met_Sign) >= 5 && (*met_Sign) <= 10){cutflow["5 <= met_Sign <= 10"] += 1;cutflow_wgt["5 <= met_Sign <= 10"] += wgt;}
       else if(cutflowstr.EqualTo("CRWW")) continue;
       if(cutflowstr.EqualTo("CRWW") && themt2 >= 60 && themt2 <= 65){cutflow["60 <= themt2 <= 65"] += 1;cutflow_wgt["60 <= themt2 <= 65"] += wgt;}
       else if(cutflowstr.EqualTo("CRWW")) continue;
    */
  
    // if(std::find(mirto_ev.begin(), mirto_ev.end(), (*EventNumber)) != mirto_ev.end()) {
    //   cout<<"Found event"<<(*EventNumber)<<endl;
    //   printf("pT1 = %.2f\n",lep1.Pt());
    //   printf("pT2 = %.2f\n",lep2.Pt());
    //   printf("mll = %.2f\n",themll);
    //   printf("MT2 = %.2f\n",themt2);
    //   printf("fl1 = %i\n",(int)lepFlavor[idx1]);
    //   printf("fl2 = %i\n",(int)lepFlavor[idx2]);
    //   printf("nbj = %i\n",n_cbj);
    //   printf("ncj = %i\n",n_clj);
    //   printf("sig = %.2f\n",(*met_Sign));
    //   printf("met = %.2f\n",(*met_Et));
    //   printf("met (tlv) = %.2f\n",met_tlv.Pt());
    // }else{
    //   return true;
    // }


    //  if(!(lep1.Pt()>25. && lep2.Pt()>25. && themll>25. && fabs(themll-zmass) < 20 && ((lepFlavor[idx1] == 1 && lepFlavor[idx2] == 1) || (lepFlavor[idx1] == 2 && lepFlavor[idx2] == 2))))return true;
    //if(!(lep1.Pt()>25. && lep2.Pt()>25. && themll>25. && themt2>100 && (signal_mu.size()+signal_el.size())==2 && ((lepFlavor[idx1] == 1 && lepFlavor[idx2] == 2) || (lepFlavor[idx1] == 2 && lepFlavor[idx2] == 1)) && n_cbj==1 && n_clj==1 && (*met_Sign)>5 && (*met_Et)>110))return true;

    // if(themt2<=60)return true;
    // if(lepCharge[idx1]*lepCharge[idx2] < 0)return true;


    //cout<<"Event passed!"<<(*EventNumber)<<endl;
    // if(!(std::find(mirto_ev.begin(), mirto_ev.end(), (*EventNumber)) != mirto_ev.end())){
    //     myfile <<(*EventNumber)<<"\n";
    //     nmirto += 1;
    //     cout<<"Found now "<<nmirto<<" events"<<endl;
    //     cout<<"Found event "<<(*EventNumber)<<endl;
    //     cout<<"Run "<<(*RunNumber)<<endl;
    //     cout<<"trigMatch_1L2LTrigOR = "<<(*trigMatch_1L2LTrigOR)<<endl;
    //     printf("pT1 = %.2f\n",lep1.Pt());
    //     printf("pT2 = %.2f\n",lep2.Pt());
    //     printf("mll = %.2f\n",themll);
    //     printf("MT2 = %.2f\n",themt2);
    //     printf("fl1 = %i\n",(int)lepFlavor[idx1]);
    //     printf("fl2 = %i\n",(int)lepFlavor[idx2]);
    //     printf("nbj = %i\n",n_cbj);
    //     printf("ncj = %i\n",n_clj);
    //     printf("sig = %.2f\n",(*met_Sign));
    //     printf("met = %.2f\n",(*met_Et));
    //     printf("met (tlv) = %.2f\n",met_tlv.Pt());
    //     printf("metphi1 = %.2f\n",met_tlv.Phi());

    // }
    // if(nmirto>20){
    //   WriteToFile();
    //   Abort("Enough!");
    // }

    if(!isData){
      if(doThreeLep)lep0Type = (truthENUMIFF)lepIFFClass[idx0];//TruthClassification(lepType[idx0], lepOrigin[idx0], lepEgMotherType[idx0], lepEgMotherOrigin[idx0], lepEgMotherPdgId[idx0], lepCharge[idx0], lepFlavor[idx0]);
      lep1Type = (truthENUMIFF)lepIFFClass[idx1];//TruthClassification(lepType[idx1], lepOrigin[idx1], lepEgMotherType[idx1], lepEgMotherOrigin[idx1], lepEgMotherPdgId[idx1], lepCharge[idx1], lepFlavor[idx1]);
      lep2Type = (truthENUMIFF)lepIFFClass[idx2];//TruthClassification(lepType[idx2], lepOrigin[idx2], lepEgMotherType[idx2], lepEgMotherOrigin[idx2], lepEgMotherPdgId[idx2], lepCharge[idx2], lepFlavor[idx2]);
    }

    fillregions.clear();
    fillregions = getFillVector(idx1,idx2,idx0);
    

    if(((*itLT).Contains("2L2J") && (*trigMatch_2LTrigOR))){
      
      if(makeHistograms)h_passed_trigMatch_2LTrig_fakeVR->Fill((*trigMatch_2LTrigOR) ? 1.0 : 0.0 );
      
      if(!doRJR && !doCoreRegions){

	

	if(lepPt[idx1] > 25 && lepPt[idx2] > 25){

	  if(*met_Et > 100 && (*met_Sign) > 6 && fabs(zmass-themll) > 20 && n_bjet <= 1 && *nJet30 >= 2)cutsregions.push_back("met_g_100_metsig_g_6_bjet_leq_1_zveto20_njet30_geq_2");
	  // if(themll > 100 && n_clj == 0 && themt2 > 80 && (*met_Et) > 110 && (*met_Sign) > 10 && n_bjet == 1)cutsregions.push_back("CRTOP");
	  // if(themll > 100 && n_clj == 0 && themt2 > 80 && (*met_Et) > 110 && (*met_Sign) > 10 && n_bjet == 1 && (*trigMatch_1L2LTrigOR))cutsregions.push_back("CRTOPL");
	  // if(themll > 100 && n_clj == 0 && themt2 >= 60 && themt2 <= 65 && (*met_Et) >= 60 && (*met_Et) <= 100 && (*met_Sign) >= 5 && (*met_Sign) <= 10 && n_bjet == 0)cutsregions.push_back("CRWW");
	  // if(themll >= 61.2 && themll <= 121.2 && n_clj == 0 && themt2 > 120 && (*met_Et) > 110 && (*met_Sign) > 10 && n_bjet == 0)cutsregions.push_back("CRVZ");
	  // if((*met_Et) > 110 && (*met_Sign) > 10 && n_bjet == 0){
	  //   if(n_clj == 1 && themll > 121.2)cutsregions.push_back("SR-SF-1J");
	  //   if(n_clj == 1 && themll > 100)cutsregions.push_back("SR-DF-1J");
	  //   if(n_clj == 0 && themll > 121.2)cutsregions.push_back("SR-SF-0J");
	  //   if(n_clj == 0 && themll > 100)cutsregions.push_back("SR-DF-0J");
	  
	  //}
	  if(*nJet30 >= 2)cutsregions.push_back("njet30_geq_2");
	  if(fabs(zmass-themll) > 20){
	    cutsregions.push_back("zveto20");
	    if(n_bjet == 0)cutsregions.push_back("zveto20_bjet_eq_0");
	    if(n_bjet >= 1)cutsregions.push_back("zveto20_bjet_geq_1");
	    if(n_bjet <= 1)cutsregions.push_back("zveto20_bjet_leq_1");
	    
	  }
	  if(*met_Et > 40){
	    cutsregions.push_back("met_g_40");
	    if(*nJet30 >= 2){
	      cutsregions.push_back("met_g_40_njet30_geq_2");
	      if(!lepECIDS[idx1] && !lepECIDS[idx2]){
		cutsregions.push_back("met_g_40_njet30_geq_2_ecids");
	      }
	      
	    }
	  }
	  
	}
      }
    }else if(((*itLT).Contains("ZMET"))/** && (*trigMatch_2LTrigOR))*/){
      if(themll >= 110 && themll <= 180 && n_bjet == 0 && met_tlv.DeltaPhi((lep1+lep2)) > 2.5 && metrel_Et > 20 && metrel_Et < 50)cutsregions.push_back("CR_Z");
      if(themll >= 70 && themll <= 110 && n_bjet == 0 && met_tlv.DeltaPhi((lep1+lep2)) > 2.5 && (*met_Sign) > 10)cutsregions.push_back("CR_Dib");
      if(themll >= 110 && themll <= 180 && met_tlv.DeltaPhi((lep1+lep2)) > 2.5 && *met_Et > 50)cutsregions.push_back("CR_Top");

      if(themll >= 110){
	if(themll < 180 && n_bjet == 0 && met_tlv.DeltaPhi((lep1+lep2)) > 2.5 && metrel_Et >= 50 && metrel_Et <= 80)cutsregions.push_back("VR_Z");
	if(themll < 180 && n_bjet == 0 && (*met_Sign) > 10)cutsregions.push_back("VR_Dib");
	if(themll > 180 && metrel_Et >= 50 && met_tlv.DeltaPhi((lep1+lep2)) > 2.5)cutsregions.push_back("VR_Top");
      }else if(themll >= 180 && n_bjet == 0 && met_tlv.DeltaPhi((lep1+lep2)) > 2.5){
	if(metrel_Et >= 50 && metrel_Et <= 100)cutsregions.push_back("SR1");
	if(metrel_Et >= 100 && metrel_Et <= 150)cutsregions.push_back("SR2");
	if(metrel_Et >= 150)cutsregions.push_back("SR3");
      } 
    }else if(!((*itLT).Contains("2L2J")) && !((*itLT).Contains("ZMET"))){
      if(isCR_Top)cutsregions.push_back("CR_Top");
      if(isCR_Dib)cutsregions.push_back("CR_Dib");
      if(isDF_1bjet)cutsregions.push_back("DF_1bjet");
      if(isDF_0jet)cutsregions.push_back("DF_0jet");
      
      /**
	 if(lepPt[idx1] > 30 && lepPt[idx2] > 9 && themll > 11 && (*met_Sign) > 3){
	 if(lepFlavor[idx1] == lepFlavor[idx2]){
	 if(fabs(zmass-themll) > 20){
	 if(n_bjet == 0)cutsregions.push_back("zveto20_bjet_eq_0");
	 if(n_bjet >= 1)cutsregions.push_back("zveto20_bjet_geq_1");
	 if((*nJet20) == 0)cutsregions.push_back("zveto20_njet20_eq_0");
	 if((*nJet20) < 2)cutsregions.push_back("zveto20_njet20_l_2");
	 }
	 }else{
	 if(n_bjet == 0)cutsregions.push_back("zveto20_bjet_eq_0");
	 if(n_bjet >= 1)cutsregions.push_back("zveto20_bjet_geq_1");
	 if((*nJet20) == 0)cutsregions.push_back("zveto20_njet20_eq_0");
	 if((*nJet20) < 2)cutsregions.push_back("zveto20_njet20_l_2");
	 }	  
	 }
      */
    }
    
    // NEW RECO:
    /**
       if(!((*itLT).Contains("2L2J")) || ((*itLT).Contains("2L2J") && (*trigMatch_2LTrigOR))){
      
       h_passed_trigMatch_2LTrig_fakeVR->Fill((*trigMatch_2LTrigOR) ? 1.0 : 0.0 );
      
       if(!doRJR && !doCoreRegions){
       if((!((*itLT).Contains("2L2J")) && lepPt[idx] > 30 && lepPt[idx2] > 9) || ((*itLT).Contains("2L2J") && (lepPt[idx1] > 25 && lepPt[idx2] > 25))){
       // if(themll > 100 && n_clj == 0 && themt2 > 80 && (*met_Et) > 110 && (*met_Sign) > 10 && n_bjet == 1)cutsregions.push_back("CRTOP");
       // if(themll > 100 && n_clj == 0 && themt2 > 80 && (*met_Et) > 110 && (*met_Sign) > 10 && n_bjet == 1 && (*trigMatch_1L2LTrigOR))cutsregions.push_back("CRTOPL");
       // if(themll > 100 && n_clj == 0 && themt2 >= 60 && themt2 <= 65 && (*met_Et) >= 60 && (*met_Et) <= 100 && (*met_Sign) >= 5 && (*met_Sign) <= 10 && n_bjet == 0)cutsregions.push_back("CRWW");
       // if(themll >= 61.2 && themll <= 121.2 && n_clj == 0 && themt2 > 120 && (*met_Et) > 110 && (*met_Sign) > 10 && n_bjet == 0)cutsregions.push_back("CRVZ");
       // if((*met_Et) > 110 && (*met_Sign) > 10 && n_bjet == 0){
       //   if(n_clj == 1 && themll > 121.2)cutsregions.push_back("SR-SF-1J");
       //   if(n_clj == 1 && themll > 100)cutsregions.push_back("SR-DF-1J");
       //   if(n_clj == 0 && themll > 121.2)cutsregions.push_back("SR-SF-0J");
       //   if(n_clj == 0 && themll > 100)cutsregions.push_back("SR-DF-0J");
	  
       //}
	  
       if(fabs(zmass-themll) > 20){
       cutsregions.push_back("zveto20");
       if(n_bjet == 0)cutsregions.push_back("zveto20_bjet_eq_0");
       if(n_bjet >= 1)cutsregions.push_back("zveto20_bjet_geq_1");
       }
       if((*itLT).Contains("2L2J")){
       if(*met_Et > 40){
       cutsregions.push_back("met_g_40");
       if(*nJet30 >= 2){
       cutsregions.push_back("met_g_40_njet30_geq_2");
       if(!lepECIDS[idx1] && !lepECIDS[idx2]){
       cutsregions.push_back("met_g_40_njet30_geq_2_ecids");
       }
       if(*met_Sign > 5)cutsregions.push_back("met_g_40_njet30_geq_2_metSig_g_5");
       }
       }else if(!((*itLT).Contains("2L2J"))){
       if(*met_Sign > 3){
       cutsregions.push_back("metsig_g_3");
       if(n_sgj>0){
       cutsregions.push_back("metsig_g_3_nsgjet_geq_0");
       }else if(n_sgj == 0){
       cutsregions.push_back("metsig_g_3_nsgjet_eq_0");
       }
       if(n_bjet >= 1)cutsregions.push_back("metsig_g_3_nbjet_geq_1");
       }
       }
       }
       }
       }
       }
    */
    //if(computeRates){
    // printf("noOR_baseline_mu = %i \n",noOR_baseline_mu.size());
    // printf("noOR_baseline_el = %i \n",noOR_baseline_el.size());
    // printf("baseline_lep = %i \n",baseline_lep.size());
    // printf("n_bjet_HF = %i \n",n_bjet_HF);
    // printf("HF_tag_mu = %i \n",HF_tag_mu.size());
    // printf("HF_probe_lep = %i \n",HF_probe_lep.size());
    if(!doCoreRegions){
      if(lep0isTight && (*met_Et)<40 && n_bjet > 0 )vec.push_back("2L01");
      if(lep0isTight &&  n_bjet == 0)vec.push_back("2L02");
      if(lep0isTight && (*met_Et)<40 && *nJet30 >= 2){ vec.push_back("2L04"); cutsregions.push_back("FAKE2L04");}
    }
    //if(lep0isTight/** && (*met_Et)<40 && *nJet30 >= 2 && !lepECIDS[idx1] && !lepECIDS[idx2]*/){vec.push_back("2L05");cutsregions.push_back("FAKE2L05");}
    //if(doThreeLep && lep0isTight){
    //if(baseline_lep.size() == 3 && (themll < 75 || themll > 105) && lepPt[idx0] > 25 && lepPt[idx1] > 25 && lepPt[idx2] > 20 && n_bjet == 0 && *H4PP > 150)vec.push_back("2L03");
    //if(baseline_lep.size() == 3 && (themll < 75 || themll > 105) && lepPt[idx0] > 25 && lepPt[idx1] > 25 && lepPt[idx2] > 20 && n_bjet == 1 && *H4PP > 150)vec.push_back("2L05");
    //if(baseline_lep.size() == 3 && (themll < 75 || themll > 105) && lepPt[idx0] > 25 && lepPt[idx1] > 25 && lepPt[idx2] > 20 && n_bjet == 0 && *PTISR > 50)vec.push_back("2L06");
    //}
    //if((*met_Et)<40 && (lepCharge[idx1]*lepCharge[idx2] > 0) && n_bjet > 0 && (lepFlavor[idx1]+lepFlavor[idx2]) == 3)vec.push_back("2L03");
    fillLandT(idx1,vec,(*itLT),true);

    fillLandT(idx2,vec,(*itLT),true);
      
    fillLandT(idx1,vec,(*itLT),true,*trigMatch_2LTrigOR,"2LTrig");
    fillLandT(idx2,vec,(*itLT),true,*trigMatch_2LTrigOR,"2LTrig");


    /**
    fillLandT(idx1,cutsregions,(*itLT),true);

    fillLandT(idx2,cutsregions,(*itLT),true);
      
    fillLandT(idx1,cutsregions,(*itLT),true,*trigMatch_2LTrigOR,"2LTrig");
    fillLandT(idx2,cutsregions,(*itLT),true,*trigMatch_2LTrigOR,"2LTrig");
    */
    /**
       fillLandT(idx1,vec,true,*trigMatch_1L2LTrigOR,"1L2LTrig");
       fillLandT(idx2,vec,true,*trigMatch_1L2LTrigOR,"1L2LTrig");

       fillLandT(idx1,vec,true,*trigMatch_1LTrigOR,"1LTrig");
       fillLandT(idx2,vec,true,*trigMatch_1LTrigOR,"1LTrig");
    */

    // fillLandT(idx1,vec,true,*trigMatch_1LTrigOR,"1LTrig");
    // fillLandT(idx2,vec,true,*trigMatch_1LTrigOR,"1LTrig");
    //cout<<"baseline_lep.size() = "<<baseline_lep.size()<<endl;
    
      if(/**(baseline_mu.size() == 2 || baseline_el.size() == 2) &&*/ lep0isTight && ((!doThreeLep && baseline_lep.size() == 2) || (doThreeLep && baseline_lep.size() == 3))){
	isREALRegion(idx0,idx1,idx2,LT_postfix);
	//cout<<"Filling REAL regions!"<<endl;
      }
      //}

      

    if(!makeHistograms/** || doCoreRegions*/)continue;

    //cout<<"a"<<endl;
    /**
       if((lepCharge[idx1]*lepCharge[idx2]>0 && lepCharge[idx1]==lepCharge[idx2])  ||
       (lepCharge[idx1]*lepCharge[idx0]>0 && lepCharge[idx1]==lepCharge[idx0])  ||
       (lepCharge[idx2]*lepCharge[idx0]>0 && lepCharge[idx2]==lepCharge[idx0]))cout<<"3L charges"<<endl;
       if(lepPt[idx1] > 25 && lepPt[idx2] > 25 && lepPt[idx0] > 20)cout<<"3L leppts"<<endl;
       if(n_bjet == 0)cout<<"nbjet"<<endl;
       if((themll < 75 || themll > 105))cout<<"mll"<<endl;
    */
    // if(doThreeLep && (((lepCharge[idx1]*lepCharge[idx2]>0 && lepCharge[idx1]==lepCharge[idx2])  ||
    // 		     (lepCharge[idx1]*lepCharge[idx0]>0 && lepCharge[idx1]==lepCharge[idx0])  ||
    // 		     (lepCharge[idx2]*lepCharge[idx0]>0 && lepCharge[idx2]==lepCharge[idx0])) &&
    // 		    (lepPt[idx1] > 25 && lepPt[idx2] > 25 && lepPt[idx0] > 20) &&
    // 		    (n_bjet == 0 && (themll < 75 || themll > 105)))){
    //   cutsregions.push_back("RJRstd");
    // }

    // if(doThreeLep){
    //   if(baseline_lep.size() == 3 && (themll < 75 || themll > 105) && lepPt[idx0] > 25 && lepPt[idx1] > 25 && lepPt[idx2] > 20 && n_bjet == 0 && *H4PP > 150) cutsregions.push_back("3LRJR_STD");
    //   if(baseline_lep.size() == 3 && (themll < 75 || themll > 105) && lepPt[idx0] > 25 && lepPt[idx1] > 25 && lepPt[idx2] > 20 && *H4PP > 150)cutsregions.push_back("3LRJR_STD_noB");
    //   if(baseline_lep.size() == 3 && (themll < 75 || themll > 105) && lepPt[idx0] > 25 && lepPt[idx1] > 25 && lepPt[idx2] > 20 && n_bjet == 0 && *PTISR > 50)cutsregions.push_back("3LRJR_CMP");
    // }else{
    //   if(lepPt[idx1] > 25 && lepPt[idx2] > 25 /*&& *trigMatch_1L2LTrigOR && (*nJet30) >= 2*/ && (*met_Et) > 100 && (*met_Et) < 200 /*&& ((*mjj) < 60 || (*mjj)>100)*/) cutsregions.push_back("2L2J-VR-com");
    //   //if(lepPt[idx1] > 25 && lepPt[idx2] > 25 && *trigMatch_1L2LTrigOR && (*nJet30) >= 2 && (*met_Et) > 100 && (*met_Et) < 200 && ((*mjj) < 60 || (*mjj)>100)) cutsregions.push_back("2L2J-VR-com");
    // }

    it    = fillregions.begin();
    itend = fillregions.end();
    for(; it != itend; it++){
      fillregions_wcuts.push_back((*it));
      for(it2 = cutsregions.begin(); it2 != cutsregions.end(); it2++){
	//cout << "Region "<<((*it)+"_"+(*it2))<<endl;
	fillregions_wcuts.push_back((*it)+"_"+(*it2)+"_"+LT_postfix);
      }
    }
  
    // for(it = fillregions.begin(); it != fillregions.end(); it++){
    //   printf("Does contain region %s\n",(*it).Data());
    // }
    Double_t fake_wgt = 1.0;
    Double_t fake_wgt_up = 1.0;
    Double_t fake_wgt_dw = 1.0;
    int vb = 0;
    //cout<<"a"<<endl;
    if(doFakeEst){
      TString yr_str;
      if(isData){
	if((*RunNumber) >= 325713 && (*RunNumber) < 348885)yr_str = "2017";
	else if ((*RunNumber) < 325713)yr_str = "2015_2016";
	else if ((*RunNumber) >= 348885)yr_str = "2018";
	else cout<<"ERROR Year for run "<<(*RunNumber)<<" is unknown"<<endl;
      }else{
	if(is1516)yr_str = "2015_2016";
	else if(is17)yr_str = "2017";
	else if(is18)yr_str = "2018";
	else cout<<"ERROR Year for run "<<(*RunNumber)<<" is unknown"<<endl;
      }

      
      
            
      TString eventStr;
      
      classifyEvent(idx1, idx2, nTT, nTl, nlT, nll);
      
      if(vb)printf("nTT = %i, nTl = %i, nlT = %i, nll = %i\n",nTT, nTl, nlT, nll);
      
      std::pair <TString,TString> ret = getFillString(idx1, idx2, idx0);
      eventStr = ret.first + ret.second;


      //       RunNumber | EventNumber | lepflav1 | lepflav2 | 
      // 284420    | 713531919   | -13      | 13       | 
      // 284420    | 305529041   | 13       | -13      | 
      // 283429    | 1018382351  | 13       | -13      | 
      // 283429    | 3392236231  | -11      | 11       | 
      // 283429    | 3394612049  | 11       | -11      | 

      
      if(*RunNumber == 284420 && *EventNumber == 713531919)printf("eventStr = %s\n",eventStr.Data());

      if(doFakeEst){
	gfw->setVariables(yr_str, LT_postfix, eventStr);
	//
	gfw->setFrac(*met_Sign, *mll, *nJet30, n_bjet, zmass, *met_Sign, lep1, lep2, met_tlv);
      
      
	//cout<<"year "<<yr_str.Data()<<", "<<(*itLT).Data()<<endl;
	if(((*itLT).Contains("2L2J") && lepPt[idx1] > 25 && lepPt[idx2] > 25) || !(*itLT).Contains("2L2J")){
	  for (auto u : GetFakeWeight::all_unc) {

	    double acc_r1 = 0.0;
	    double acc_r2 = 0.0;
	    double BDT_r1 = 0.0;
	    double BDT_r2 = 0.0;
	    
	    glob_r1 = gfw->getReal(lepPt[idx1], lepEta[idx1], fabs(lepFlavor[idx1]), triggercat1, u);
	    glob_r2 = gfw->getReal(lepPt[idx2], lepEta[idx2], fabs(lepFlavor[idx2]), triggercat2, u);
	    glob_f1 = gfw->getFake(lepPt[idx1], lepEta[idx1], fabs(lepFlavor[idx1]), triggercat1, u);
	    glob_f2 = gfw->getFake(lepPt[idx2], lepEta[idx2], fabs(lepFlavor[idx2]), triggercat2, u);
	    
	    if(gfw->getUncKey(u).EqualTo("NOM")){
	      int nBDT = 0;
	      for (std::string bdt_score : bdt_vec) {
		TString key = bdt_score;
		if(!key.EqualTo("BDTVVDeltaM100_90"))continue;
		BDT_r1 = gfw->getReal(lepPt[idx1], lepEta[idx1], fabs(lepFlavor[idx1]), triggercat1, u, key, BDTweight[key]);
		BDT_r2 = gfw->getReal(lepPt[idx2], lepEta[idx2], fabs(lepFlavor[idx2]), triggercat2, u, key, BDTweight[key]);

		h_real1_diff[key]->Fill(BDT_r1-glob_r1,wgt);
		h_real2_diff[key]->Fill(BDT_r2-glob_r2,wgt);
		  
		acc_r1 += BDT_r1;
		acc_r2 += BDT_r1;
		
		nBDT += 1;
		//printf("%s = %.2f, r1 = %.2f, r2 = %.2f\n",key.Data(),BDTweight[key],glob_r1,glob_r2);
	      }

	      // If one wants to use the real effs from the BDTscores
	      // glob_r1 = acc_r1/nBDT;
	      // glob_r2 = acc_r2/nBDT;
	      //printf("Average: r1 = %.2f, r2 = %.2f\n",acc_r1/bdt_vec.size(),acc_r2/bdt_vec.size());
	    }

	    glob_f1 = glob_f1 <= 0.0 ? 0.01 : glob_f1;
	    glob_f2 = glob_f2 <= 0.0 ? 0.01 : glob_f2;
	    glob_r1 = glob_r1 >= 1.0 ? 0.99 : glob_r1;
	    glob_r2 = glob_r2 >= 1.0 ? 0.00 : glob_r2;
	  
	    glob_r_wgt[gfw->getUncKey(u)]  = mtx->N4_RR_TT(glob_r1,glob_f1,glob_r2,glob_f2,nTT,nTl,nlT,nll);
	    glob_rf_wgt[gfw->getUncKey(u)] = mtx->N4_RF_TT(glob_r1,glob_f1,glob_r2,glob_f2,nTT,nTl,nlT,nll);
	    glob_fr_wgt[gfw->getUncKey(u)] = mtx->N4_FR_TT(glob_r1,glob_f1,glob_r2,glob_f2,nTT,nTl,nlT,nll);
	    glob_f_wgt[gfw->getUncKey(u)]  = mtx->N4_FF_TT(glob_r1,glob_f1,glob_r2,glob_f2,nTT,nTl,nlT,nll);
	  
  
  
	    // *WeightEvents = b_fake_wgt;
	    b_fake_wgt[gfw->getUncKey(u)] = glob_f_wgt[gfw->getUncKey(u)]+glob_fr_wgt[gfw->getUncKey(u)]+glob_rf_wgt[gfw->getUncKey(u)];
	  
	    //cout<<"fake_wgt["<<gfw->getUncKey(u)<<"] = "<<b_fake_wgt[gfw->getUncKey(u)]<<endl;

	    //if(*RunNumber == 350067 && *EventNumber == 1047493169){
	    //if(gfw->getUncKey(u).EqualTo("NOM")){
	    if(fabs(b_fake_wgt[gfw->getUncKey(u)]) > 15 && (lepPt[idx1] > 25 && lepPt[idx2] > 9 && themll > 11 && *nJet30 >= 2 && n_bjet == 0 && *met_Sign > 3)){
	      gfw->setVB(1);
	      cout<<"----------------------- !! START: LARGE FAKE WEIGHT "<<gfw->getUncKey(u)<<" !! ---------------------------"<<endl;
	    
	      cout<<"RunNumber = "<<*RunNumber<<endl;
	      cout<<"EventNumber = "<<*EventNumber<<endl;
	    
	      cout<<"fake_wgt["<<gfw->getUncKey(u)<<"] = "<<b_fake_wgt[gfw->getUncKey(u)]<<endl;
	    
	      classifyEvent(idx1, idx2, nTT, nTl, nlT, nll);

	      printf("triggercat1 = %s\n",triggercat1.Data());
	      printf("triggercat2 = %s\n",triggercat2.Data());
	    
	      glob_r1 = gfw->getReal(lepPt[idx1], lepEta[idx1], fabs(lepFlavor[idx1]), triggercat1, u);
	      glob_r2 = gfw->getReal(lepPt[idx2], lepEta[idx2], fabs(lepFlavor[idx2]), triggercat2, u);
	      glob_f1 = gfw->getFake(lepPt[idx1], lepEta[idx1], fabs(lepFlavor[idx1]), triggercat1, u);
	      glob_f2 = gfw->getFake(lepPt[idx2], lepEta[idx2], fabs(lepFlavor[idx2]), triggercat2, u);
	    
	      glob_f1 = glob_f1 <= 0.0 ? 0.01 : glob_f1;
	      glob_f2 = glob_f2 <= 0.0 ? 0.01 : glob_f2;
	      glob_r1 = glob_r1 >= 1.0 ? 0.99 : glob_r1;
	      glob_r2 = glob_r2 >= 1.0 ? 0.99 : glob_r2;
	    
	      printf("f1 = %.4f, f2 = %.4f, r1 = %.4f, r2 = %.4f\n",glob_f1,glob_f2,glob_r1,glob_r2);
	    
	      glob_r_wgt[gfw->getUncKey(u)]  = mtx->N4_RR_TT(glob_r1,glob_f1,glob_r2,glob_f2,nTT,nTl,nlT,nll);
	      glob_rf_wgt[gfw->getUncKey(u)] = mtx->N4_RF_TT(glob_r1,glob_f1,glob_r2,glob_f2,nTT,nTl,nlT,nll);
	      glob_fr_wgt[gfw->getUncKey(u)] = mtx->N4_FR_TT(glob_r1,glob_f1,glob_r2,glob_f2,nTT,nTl,nlT,nll);
	      glob_f_wgt[gfw->getUncKey(u)]  = mtx->N4_FF_TT(glob_r1,glob_f1,glob_r2,glob_f2,nTT,nTl,nlT,nll);
	    
	      printf("Lepton %i (pT = %.2f, typ = %s) matched to \n",idx1,lepPt[idx1]/1000.,fabs(lepFlavor[idx1]) == 13 ? "mu" : "el");
	      // for(std::vector<TString>::iterator iter2 = trigmatches[idx1].begin(); iter2 != trigmatches[idx1].end(); ++iter2){
	      //   printf("%s\n",(*iter2).Data());
	      // }
	      printf("Lepton %i (pT = %.2f, typ = %s) matched to \n",idx2,lepPt[idx2]/1000.,fabs(lepFlavor[idx2]) == 13 ? "mu" : "el");
	      // for(std::vector<TString>::iterator iter2 = trigmatches[idx2].begin(); iter2 != trigmatches[idx2].end(); ++iter2){
	      //   printf("%s\n",(*iter2).Data());
	      // }
	    
	      if(glob_f1 > glob_r1)printf("ERROR \t f1 = %.2f > r1 = %.2f\n",glob_f1,glob_r1);
	      if(glob_f2 > glob_r2)printf("ERROR \t f2 = %.2f > r2 = %.2f\n",glob_f2,glob_r2);
	    
	      printf("f_hf = %.2f, f_lf = %.2f, f_co = %.2f\n",gfw->getFracHF(),gfw->getFracLF(),gfw->getFracCO());
	      printf("r_wgt = %.2f, rf_wgt = %.2f, fr_wgt = %.2f, f_wgt = %.2f\n",glob_r_wgt[gfw->getUncKey(u)],glob_rf_wgt[gfw->getUncKey(u)],glob_fr_wgt[gfw->getUncKey(u)],glob_f_wgt[gfw->getUncKey(u)]);
	      printf("nTT = %i, nTl = %i, nlT = %i, nll = %i\n",nTT,nTl,nlT,nll);
	    
	      cout<<"----------------------- !! STOP: LARGE FAKE WEIGHT "<<gfw->getUncKey(u)<<" !! ---------------------------"<<endl; 
	      gfw->setVB(0);
	    }
	  }

	  // fake_wgt_vec[LT_postfix] = getFakeWeight1D(idx1,idx2,idx0,yr_str,LT_postfix,1);
	  // fake_wgt                 = fake_wgt_vec[LT_postfix].at(0);
	  // fake_wgt_up  = fake_wgt_vec[LT_postfix].size() > 1 ? fake_wgt_vec[LT_postfix].at(1) : 1.0;
	  // fake_wgt_dw  = fake_wgt_vec[LT_postfix].size() > 2 ? fake_wgt_vec[LT_postfix].at(2) : 1.0;
	  // h_fake_weights->Fill(fake_wgt);
	
	  // if(fabs(fake_wgt-1.0)<0.001){
	  //   printf("wgt = %.2f : %s : %.2f : %.2f : nlep = %lu\n", fake_wgt, yr_str.Data(),lepPt[idx1],lepPt[idx2],(baseline_mu.size()+baseline_el.size()));
	  //   printf("r1 = %.4f, r2 = %.4f, f1 = %.4f, f1 = %.4f\n",glob_r1,glob_r2,glob_f1,glob_f2);
	  //   printf("f_hf = %.2f, f_lf = %.2f, f_co = %.2f",gfw->getFracHF(),gfw->getFracLF(),gfw->getFracCO());
	  //   printf("r_wgt = %.2f, rf_wgt = %.2f, fr_wgt = %.2f, f_wgt = %.2f\n",glob_r_wgt,glob_rf_wgt,glob_fr_wgt,glob_f_wgt);
	  //   gfw->setVB(1);
	  //   fake_wgt_vec[LT_postfix] = getFakeWeight1D(idx1,idx2,idx0,yr_str,LT_postfix,1);
	  //   //gfw->setVB(0);
	  //   for(unsigned int i = 0; i < fillregions.size(); i++){
	  //     printf("\n%s\n",fillregions.at(i).Data());
	  //   }
	  // }
	}
	fake_wgt = b_fake_wgt["NOM"];
      }
      
      /**
	 MMweight = fake_wgt;
	 MMweight_up = fake_wgt_up;
	 MMweight_dw = fake_wgt_dw;

	 centraltree->GetEntry(central_nentries);
	 centraltree->Fill();
	 central_nentries += 1;
      */

      //cout<<"a"<<endl;
      // check
      // if(fake > 0 && fake_wgt_up < 1)cout<<"Fake is "<<fake_wgt<<" and Wgt up is smaller than one "<<fake_wgt_up<<endl;
      // else if(fake < 0 && fake_wgt_dw > 1)cout<<"Fake is "<<fake_wgt<<" and wgtt dw is larger  than one "<<fake_wgt_dw<<endl;
      // else if(fake > 0 && fake_wgt_dw < 1)cout<<"Fake is "<<fake_wgt<<" and wgtt dw is smaller  than one "<<fake_wgt_dw<<endl;
      // 	      else if(fake < 0 && fake_wgt_up < 1)cout<<"Fake is "<<fake_wgt<<" and wgtt dw is smaller  than one "<<fake_wgt_dw<<endl;


      /** }else{
	  if((*RandomRunNumber) >= 320000)yr_str = "2017";
	  else if ((*RandomRunNumber) < 320000)yr_str = "2015_2016";
	  else yr_str = "2018";
	  //cout<<"year "<<yr_str.Data()<<endl;
	  fake_wgt = getFakeWeight(idx1,idx2,idx0,yr_str,0);
	  }
      */

    }

       //cout<<"b"<<endl;
    
    if(makeNtup && (*itLT).EqualTo(looseTightDef.back())){

      //cout<<"Key is " << (*itLT).Data() << " and last entry is " << looseTightDef.back().Data()<<endl;

      Double_t lepflav[3] = {-1,-1,-1};

      if(doThreeLep){
	lepflav[0] = lepFlavor[idx0] == 1 ? 11.0*lepCharge[idx0] : 13.0*lepCharge[idx0];
	lepflav[1] = lepFlavor[idx1] == 1 ? 11.0*lepCharge[idx1] : 13.0*lepCharge[idx1];
	lepflav[2] = lepFlavor[idx2] == 1 ? 11.0*lepCharge[idx2] : 13.0*lepCharge[idx2];
      }else{
	lepflav[0] = lepFlavor[idx1] == 1 ? 11.0*lepCharge[idx1] : 13.0*lepCharge[idx1];
	lepflav[1] = lepFlavor[idx2] == 1 ? 11.0*lepCharge[idx2] : 13.0*lepCharge[idx2];
	lepflav[2] = -1;
      }

      HFT->RunNumber = (*RandomRunNumber);
      HFT->DatasetNumber = (*DatasetNumber);
      HFT->EventNumber = (*EventNumber);
      HFT->lept1Pt = doThreeLep ? lepPt[idx0]*gev : lepPt[idx1]*gev;
      HFT->lept2Pt = doThreeLep ? lepPt[idx1]*gev : lepPt[idx2]*gev;
      HFT->lept3Pt = doThreeLep ? lepPt[idx2]*gev : -1;
      HFT->lept1Eta = doThreeLep ? lepEta[idx0] : lepEta[idx1];
      HFT->lept2Eta = doThreeLep ? lepEta[idx1] : lepEta[idx2];
      HFT->lept3Eta = doThreeLep ? lepEta[idx2] : -1;
      HFT->lept1Phi  = doThreeLep ? lepPhi[idx0] : lepPhi[idx1];
      HFT->lept2Phi  = doThreeLep ? lepPhi[idx1] : lepPhi[idx2];
      HFT->lept3Phi  = doThreeLep ? lepPhi[idx2] : -1;
      HFT->lept1E = doThreeLep ? lep0.E()*gev : lep1.E()*gev;
      HFT->lept2E = doThreeLep ? lep1.E()*gev : lep2.E()*gev;
      HFT->lept3E = doThreeLep ? lep2.E()*gev : -1;
      HFT->lept1Flav  = lepflav[0];
      HFT->lept2Flav  = lepflav[1];
      HFT->lept3Flav  = lepflav[2];
      HFT->metPhi = met_tlv.Phi();
      HFT->MET = (*met_Et)*gev;
      HFT->METsig = (*met_Sign);
      HFT->nLeps = baseline_lep.size();
      HFT->nSigLeps = (signal_mu.size()+signal_el.size());
      for(std::vector< TString >::iterator itLT2 = looseTightDef.begin(); itLT2 != looseTightDef.end(); itLT2++){
	//if(fake_wgt_vec[*itLT2].size() > 0)cout<<"Filling weight for "<<(*itLT2).Data()<<" with weight "<<fake_wgt_vec[*itLT2].at(0)<<endl;
	//else cout<< "No weights for "<< (*itLT2) << endl;
	HFT->eventWeight.push_back(fake_wgt_vec[*itLT2].size() > 0 ? fake_wgt_vec[*itLT2].at(0) : 1.0);
	HFT->syst_MM_up.push_back(fake_wgt_vec[*itLT2].size() > 1 ? fake_wgt_vec[*itLT2].at(1): 1.0);
	HFT->syst_MM_down.push_back(fake_wgt_vec[*itLT2].size() > 2 ? fake_wgt_vec[*itLT2].at(2) : 1.0);
	HFT->MM_key.push_back(*itLT2);
      }
      HFT->WeightEvents = 1.0;
      HFT->WeightEventsEff = 1.0;
      HFT->WeightEventsSF = 1.0;
      HFT->WeightEventsbTag = 1.0;
      HFT->WeightEventsJVT = 1.0;
      HFT->WeightEventselSF = 1.0;
      HFT->WeightEventsmuSF = 1.0;
      HFT->L2Mll = themll*gev;
      HFT->L2MT2 = themt2*gev;
      HFT->L2isEMU = ((lepFlavor[idx1] == 1 && lepFlavor[idx2] == 2) || (lepFlavor[idx1] == 2 && lepFlavor[idx2] == 1));
      HFT->L2isEE = (lepFlavor[idx1] == 1 && lepFlavor[idx2] == 1);
      HFT->L2isMUMU = (lepFlavor[idx1] == 2 && lepFlavor[idx2] == 2);
      HFT->L2nCentralLightJet = n_clj;
      HFT->L2nCentralBJet = n_cbj;
      HFT->L2nForwardJet = n_fwj;
      HFT->L2nCentralLightJet30 = n_clj30;
      HFT->L2nCentralBJet30 = n_cbj30;
      HFT->L2nForwardJet30 = n_fwj30;
      HFT->L2nCentralLightJet40 = n_clj40;
      HFT->L2nCentralBJet40 = n_cbj40;
      HFT->L2nForwardJet40 = n_fwj40;
      HFT->L2nCentralLightJet50 = n_clj50;
      HFT->L2nCentralBJet50 = n_cbj50;
      HFT->L2nForwardJet50 = n_fwj50;
      HFT->L2nCentralLightJet60 = n_clj60;
      HFT->L2nCentralBJet60 = n_cbj60;
      HFT->L2nForwardJet60 = n_fwj60;
      HFT->L2SiLepTrigger = true;//(*trigMatch_1LTrigOR);
      HFT->L2DiLepTrigger = (*trigMatch_2LTrigOR);
      HFT->jet1Pt = jetPt.GetSize() > 0 ? jetPt[0]*gev : -1;
      HFT->jet2Pt = jetPt.GetSize() > 1 ? jetPt[1]*gev : -1;
      HFT->jet3Pt = jetPt.GetSize() > 2 ? jetPt[2]*gev : -1;
      HFT->jet1Eta  = jetEta.GetSize() > 0 ? jetEta[0] : -999;
      HFT->jet2Eta  = jetEta.GetSize() > 1 ? jetEta[1] : -999;
      HFT->jet3Eta  = jetEta.GetSize() > 2 ? jetEta[2] : -999;
      HFT->jet1Phi = jetPhi.GetSize() > 0 ? jetPhi[0] : -999;
      HFT->jet2Phi = jetPhi.GetSize() > 1 ? jetPhi[1] : -999;
      HFT->jet3Phi = jetPhi.GetSize() > 2 ? jetPhi[2] : -999;
      HFT->b_mu = (*mu);

      // HFT->eventweight = wgt;
      // HFT->runNumber = (*RandomRunNumber);
      // HFT->eventNumber = (*EventNumber);
      // HFT->mcChannel = rnum.Atoi();
      // HFT->NumPrimaryVertices = (*nVtx);
      // HFT->eventweightNoPRW = ((*genWeight) * (*eventWeight) * (*leptonWeight) * (*jvtWeight) * (*bTagWeight) * /**(*FFWeight) **/ ((*RandomRunNumber) > 320000 ? 43800 : 36200));
      // HFT->xsecUp = (*xsec);
      // HFT->xsecDown = (*xsec);
      // HFT->njets = n_blj;
      // HFT->nbjets = n_bjet;
      // HFT->njets30 = n_blj30;
      // HFT->njets40 = n_blj40;
      // HFT->njets50 = n_blj50;
      // HFT->njets60 = n_blj60;
      // HFT->L2nCentralLightJets = n_clj;
      // HFT->L2nCentralLightJets30 = n_clj30;
      // HFT->L2nCentralLightJets40 = n_clj40;
      // HFT->L2nCentralLightJets50 = n_clj50;
      // HFT->L2nCentralLightJets60 = n_clj60;
      // HFT->L2nCentralBJets = n_cbj;
      // HFT->L2nForwardJets = n_fwj;
      // HFT->L2nForwardJets30 = n_fwj30;
      // HFT->L2nForwardJets40 = n_fwj40;
      // HFT->L2nForwardJets50 = n_fwj50;
      // HFT->L2nForwardJets60 = n_fwj60;
      //     HFT->lept1Type = lepType[idx1];
      //     HFT->lept2Type = lepType[idx2];
      //     HFT->lept1Origin = lepOrigin[idx1];
      //     HFT->lept2Origin = lepOrigin[idx2];
      //     HFT->lept1Pt = lepPt[idx1];
      //     HFT->lept2Pt = lepPt[idx2];
      //     HFT->lept1Eta = lepEta[idx1];
      //     HFT->lept2Eta = lepEta[idx2];
      //     HFT->lept1Phi = lepPhi[idx1];
      //     HFT->lept2Phi = lepPhi[idx2];
      //     HFT->lept1Flav = lepFlavor[idx1];
      //     HFT->lept2Flav = lepFlavor[idx2];
      //     HFT->lept1q = lepCharge[idx1];
      //     HFT->lept2q = lepCharge[idx2];
      //     HFT->lept1OR = lepPassOR->at(idx1);
      //     HFT->lept2OR = lepPassOR->at(idx2);
      //     HFT->lept1baseline = (lepPassOR->at(idx1) && isL(idx1));
      //     HFT->lept2baseline = (lepPassOR->at(idx2) && isL(idx2));
      //     HFT->lept1signal = (lepPassOR->at(idx1) && isT(idx1));
      //     HFT->lept2signal = (lepPassOR->at(idx2) && isT(idx2));
      //     HFT->lept1MCtruthClass = (int)lep1Type;
      //     HFT->lept2MCtruthClass = (int)lep2Type;
      //     HFT->metPhi = met_tlv.Phi();
      //     HFT->MET = (*met_Et);
      //     HFT->METsig = 0.0;
      //     HFT->L2dPhib = met_tlv.DeltaPhi(lep1+lep2+met_tlv);
      //     HFT->L2dPhi_MET = met_tlv.DeltaPhi(lep1+lep2);
      //     HFT->L2LtMET = lep1.Pt()+lep2.Pt()+met_tlv.Pt();
      //     HFT->nSigLep = signal_mu.size() + signal_el.size();
      //     HFT->L2Mll = themll;
      //     HFT->L2MT2 = themt2;
      //     HFT->pbllmod = (lep1+lep2+met_tlv).Pt();
      //     HFT->L2isEMU  = ((lepFlavor[idx1] == 1 && lepFlavor[idx2] == 2) || (lepFlavor[idx1] == 2 && lepFlavor[idx2] == 1));
      //     HFT->L2isMUMU = (lepFlavor[idx1] == 2 && lepFlavor[idx2] == 2);
      //     HFT->L2isEE   = (lepFlavor[idx1] == 1 && lepFlavor[idx2] == 1);
      //     HFT->L2dileptonpt = (lep1+lep2).Pt();
      //     HFT->nBaseLep   = baseline_mu.size() + baseline_el.size();
      //     HFT->isHF = (lep1Type == 1 && lep2Type == 1) ? true : false;
      //     HFT->isLF = (lep1Type == 2 && lep2Type == 2) ? true : false;
      //     HFT->isCF = (lep1Type == 0 && lep2Type == 0) ? true : false;
      //     HFT->isCO = (lep1Type == 3 && lep2Type == 3) ? true : false;
      //     HFT->isRL = (lep1Type == 4 && lep2Type == 4) ? true : false;
      //     HFT->isUN = (lep1Type == 5 && lep2Type == 5) ? true : false;

      //     if(lep1Type == 4 && lep2Type == 4)HFT->MMcat = 0;
      //     else if(lep1Type == 4 && lep2Type != 4)HFT->MMcat = 1;
      //     else if(lep1Type != 4 && lep2Type == 4)HFT->MMcat = 2;
      //     else if(lep1Type != 4 && lep2Type != 4)HFT->MMcat = 3;
      //     else HFT->MMcat = 4;

      //     if((lep1Type == 1 || lep2Type == 1))HFT->eventCat = 1;
      //     else if((lep1Type == 2 || lep2Type == 2))HFT->eventCat = 2;
      //     else if((lep1Type == 3 || lep2Type == 3))HFT->eventCat = 3;
      //     else if((lep1Type == 4 || lep2Type == 4))HFT->eventCat = 4;
      //     else if((lep1Type == 5 || lep2Type == 5))HFT->eventCat = 5;
      //     else if((lep1Type == 0 || lep2Type == 0))HFT->eventCat = 0;
      //     else HFT->eventCat = 6;

      //     //HFT->isTraining = rnd->Uniform() <= 0.5 ? true : false;
      //     if(!isData){
      // HFT->isTraining = rnd->Uniform() <= 0.5 ? true : false;
      //     }else{
      //       HFT->isTraining = false;
      //     }
      //     // if((lep1Type >= 0 && lep1Type <= 5) && (lep2Type >= 0 && lep2Type <= 5) && (lep1Type == lep2Type)){
      //     //   HFT->isTraining = rnd->Uniform() <= 0.5 ? true : false;
      //     // }else{
      //     //   HFT->isTraining = false;
      //     // }
      //cout<<"filling"<<endl;
      //tree3->Fill();
      chain->GetEntry( fReader.GetCurrentEntry() );
      newtree->Fill();
      HFT->WriteTree();
      if(makeMYNtup)fillMYTree(idx1, idx2, fake_wgt_vec);
      if(make2L2JNtup && looseTightDef.size()==1)fill2L2JTree(idx1, idx2, b_fake_wgt);
     
    }
    /**
       cout<<"Will look at following regions:"<<endl;
       for(it = fillregions_wcuts.begin(); it != fillregions_wcuts.end(); it++){
       cout<<(*it).Data()<<endl;
       }
       cout<<"------------------"<<endl;
    */
    
    // cout<<"lep1isTight = "<<(lep1isTight ? "True" : "False");
    // cout<<endl;
    // cout<<"lep2isTight = "<<(lep2isTight ? "True" : "False");
    // cout<<endl;

    for(it = fillregions_wcuts.begin(); it != fillregions_wcuts.end(); it++){

      //cout<<"Checking "<<(*it).Data()<<endl;
      if(!(h_lep_pT_nL.find((*it)) == h_lep_pT_nL.end())){
	//cout<<"Filling "<<(*it).Data()<<endl;
	fillFakeRateHist((*it), idx1, lep1isTight);
	fillFakeRateHist((*it), idx2, lep2isTight);
	// if(doThreeLep){
	//   h_lep_pT0_nL[(*it)]->Fill(lepPt[idx0],wgt);
	//   if(lep0isTight)h_lep_pT0_nT[(*it)]->Fill(lepPt[idx0],wgt);
	// }

      

	//cout<<"idx 0 = "<<idx0<<endl;
	if(doThreeLep)fillTruthInfo((*it),idx0,lep0isTight,wgt);
	//cout<<"b"<<endl;
	if(!lep0isTight)continue;
       
	for (TString bdt_name : bdt_vec) {
	  h_lep_BDT_nL[bdt_name+"_"+(*it)]->Fill(BDTweight[bdt_name],wgt);
	  if(nTT)h_evt_BDT_TT[bdt_name+"_"+(*it)]->Fill(BDTweight[bdt_name],wgt);
	  if(nTl)h_evt_BDT_Tl[bdt_name+"_"+(*it)]->Fill(BDTweight[bdt_name],wgt);
	  if(nlT)h_evt_BDT_lT[bdt_name+"_"+(*it)]->Fill(BDTweight[bdt_name],wgt);
	  if(nll)h_evt_BDT_ll[bdt_name+"_"+(*it)]->Fill(BDTweight[bdt_name],wgt);
	}
	//cout<<"Filling h_lep_pT_nL for "<< (*it) << " in main body" << endl;
	/**
	h_lep_pT_nL[(*it)]->Fill(lepPt[idx1],wgt);
	h_lep_njet30_nL[(*it)]->Fill(*nJet30,wgt);
	h_lep_eta_nL[(*it)]->Fill(lepEta[idx1],wgt);
	h_lep_pT_eta_nL[(*it)]->Fill(lepPt[idx1],lepEta[idx1],wgt);
	h_lep_mu_pT_nL[(*it)]->Fill(*mu,lepPt[idx1],wgt);
	*/
	if(lep1isTight){
	  for (TString bdt_name : bdt_vec) {
	    h_lep_BDT_nT[bdt_name+"_"+(*it)]->Fill(BDTweight[bdt_name],wgt);
	  }
	  /**
	  h_lep_pT_nT[(*it)]->Fill(lepPt[idx1],wgt);
	  h_lep_njet30_nT[(*it)]->Fill(*nJet30,wgt);
	  h_lep_eta_nT[(*it)]->Fill(lepEta[idx1],wgt);
	  h_lep_pT_eta_nT[(*it)]->Fill(lepPt[idx1],lepEta[idx1],wgt);
	  h_lep_mu_pT_nT[(*it)]->Fill(*mu,lepPt[idx1],wgt);
	  */
	}
	for (TString bdt_name : bdt_vec) {
	  h_lep_BDT_nL[bdt_name+"_"+(*it)]->Fill(BDTweight[bdt_name],wgt);
	}
	/**
	h_lep_pT_nL[(*it)]->Fill(lepPt[idx2],wgt);
	h_lep_njet30_nL[(*it)]->Fill(*nJet30,wgt);
	h_lep_eta_nL[(*it)]->Fill(lepEta[idx2],wgt);
	h_lep_pT_eta_nL[(*it)]->Fill(lepPt[idx2],lepEta[idx2],wgt);
	h_lep_mu_pT_nL[(*it)]->Fill(*mu,lepPt[idx2],wgt);
	*/
	if(lep2isTight){
	  for (TString bdt_name : bdt_vec) {
	    h_lep_BDT_nT[bdt_name+"_"+(*it)]->Fill(BDTweight[bdt_name],wgt);
	  }
	  /**
	  h_lep_pT_nT[(*it)]->Fill(lepPt[idx2],wgt);
	  h_lep_njet30_nT[(*it)]->Fill(*nJet30,wgt);
	  h_lep_eta_nT[(*it)]->Fill(lepEta[idx2],wgt);
	  h_lep_pT_eta_nT[(*it)]->Fill(lepPt[idx2],lepEta[idx2],wgt);
	  h_lep_mu_pT_nT[(*it)]->Fill(*mu,lepPt[idx2],wgt);
	  */
	}
      }//else{
      // 	printf("WARNING \t Did not find %s for h_lep_pT_nL\n",((*it)).Data());

      // }
      //
  
      if(!(h_lep_nbljetpreOR.find((*it)) == h_lep_nbljetpreOR.end())){
	//cout<<"Region "<<(*it)<<endl;
	h_lep_nbljetpreOR[(*it)]->Fill(n_bljpreOR,wgt);
      }
      //cout<<"Region after "<<(*it)<<endl;
      
      if(h_lep_pT1.find((*it)) == h_lep_pT1.end()){
	//	printf("WARNING \t Did not find %s in h_lep_pT1\n",(*it).Data());
	continue;
      }

      //cout<<"Region after "<<(*it)<<" and isData = "<<isData<<endl;
      
      if(lep1isTight && lep2isTight){
	//cout<<"1) "<<(*it).Data()<<endl;
	if(doThreeLep)h_lep_pT0[(*it)]->Fill(lepPt[idx0],wgt);
	h_lep_pT1[(*it)]->Fill(lepPt[idx1],wgt);
	h_lep_pT2[(*it)]->Fill(lepPt[idx2],wgt);
	h_lep_mt2[(*it)]->Fill(themt2,wgt);
	h_lep_metsig[(*it)]->Fill(*met_Sign,wgt);
	if(doThreeLep)h_lep_mtw[(*it)]->Fill(mtw,wgt);
	h_lep_met[(*it)]->Fill(*met_Et,wgt);
	h_lep_mll[(*it)]->Fill(themll,wgt);
	h_lep_njet30[(*it)]->Fill(*nJet30,wgt);
	h_lep_nbj77[(*it)]->Fill(n_bjet77,wgt);
	h_lep_nbj85[(*it)]->Fill(n_bjet85,wgt);
      }
      if(doFakeEst){
	//cout<<"2) "<<(*it).Data()<<"_estfake with weight "<<wgt*fake_wgt<<endl;
	if(doThreeLep)h_lep_pT0[(*it)+"_estfake"]->Fill(lepPt[idx0],wgt*fake_wgt);
	h_lep_pT1[(*it)+"_estfake"]->Fill(lepPt[idx1],wgt*fake_wgt);
	h_lep_pT2[(*it)+"_estfake"]->Fill(lepPt[idx2],wgt*fake_wgt);
	h_lep_mt2[(*it)+"_estfake"]->Fill(themt2,wgt*fake_wgt);
	h_lep_metsig[(*it)+"_estfake"]->Fill(*met_Sign,wgt*fake_wgt);
	if(doThreeLep)h_lep_mtw[(*it)+"_estfake"]->Fill(mtw,wgt*fake_wgt);
	h_lep_met[(*it)+"_estfake"]->Fill(*met_Et,wgt*fake_wgt);
	h_lep_mll[(*it)+"_estfake"]->Fill(themll,wgt*fake_wgt);
	h_lep_njet30[(*it)+"_estfake"]->Fill(*nJet30,wgt*fake_wgt);
	h_lep_nbj77[(*it)+"_estfake"]->Fill(n_bjet77,wgt*fake_wgt);
	h_lep_nbj85[(*it)+"_estfake"]->Fill(n_bjet85,wgt*fake_wgt);

	if(doThreeLep)h_lep_pT0[(*it)+"_estfake_UP"]->Fill(lepPt[idx0],wgt*fake_wgt_up*fake_wgt);
	h_lep_pT1[(*it)+"_estfake_UP"]->Fill(lepPt[idx1],wgt*fake_wgt_up*fake_wgt);
	h_lep_pT2[(*it)+"_estfake_UP"]->Fill(lepPt[idx2],wgt*fake_wgt_up*fake_wgt);
	h_lep_mt2[(*it)+"_estfake_UP"]->Fill(themt2,wgt*fake_wgt_up*fake_wgt);
	h_lep_metsig[(*it)+"_estfake_UP"]->Fill(*met_Sign,wgt*fake_wgt_up*fake_wgt);
	if(doThreeLep)h_lep_mtw[(*it)+"_estfake_UP"]->Fill(mtw,wgt*fake_wgt_up*fake_wgt);
	h_lep_met[(*it)+"_estfake_UP"]->Fill(*met_Et,wgt*fake_wgt_up*fake_wgt);
	h_lep_mll[(*it)+"_estfake_UP"]->Fill(themll,wgt*fake_wgt_up*fake_wgt);
	h_lep_njet30[(*it)+"_estfake_UP"]->Fill(*nJet30,wgt*fake_wgt_up*fake_wgt);
	h_lep_nbj77[(*it)+"_estfake_UP"]->Fill(n_bjet77,wgt*fake_wgt_up*fake_wgt);
	h_lep_nbj85[(*it)+"_estfake_UP"]->Fill(n_bjet85,wgt*fake_wgt_up*fake_wgt);

	if(doThreeLep)h_lep_pT0[(*it)+"_estfake_DW"]->Fill(lepPt[idx0],wgt*fake_wgt_dw*fake_wgt);
	h_lep_pT1[(*it)+"_estfake_DW"]->Fill(lepPt[idx1],wgt*fake_wgt_dw*fake_wgt);
	h_lep_pT2[(*it)+"_estfake_DW"]->Fill(lepPt[idx2],wgt*fake_wgt_dw*fake_wgt);
	h_lep_mt2[(*it)+"_estfake_DW"]->Fill(themt2,wgt*fake_wgt_dw*fake_wgt);
	h_lep_metsig[(*it)+"_estfake_DW"]->Fill(*met_Sign,wgt*fake_wgt_dw*fake_wgt);
	if(doThreeLep)h_lep_mtw[(*it)+"_estfake_DW"]->Fill(mtw,wgt*fake_wgt_dw*fake_wgt);
	h_lep_met[(*it)+"_estfake_DW"]->Fill(*met_Et,wgt*fake_wgt_dw*fake_wgt);
	h_lep_mll[(*it)+"_estfake_DW"]->Fill(themll,wgt*fake_wgt_dw*fake_wgt);
	h_lep_njet30[(*it)+"_estfake_DW"]->Fill(*nJet30,wgt*fake_wgt_dw*fake_wgt);
	h_lep_nbj77[(*it)+"_estfake_DW"]->Fill(n_bjet77,wgt*fake_wgt_dw*fake_wgt);
	h_lep_nbj85[(*it)+"_estfake_DW"]->Fill(n_bjet85,wgt*fake_wgt_dw*fake_wgt);

	int nTT = 0;
	int nTl = 0;
	int nlT = 0;
	int nll = 0;
	classifyEvent(idx1, idx2, nTT, nTl, nlT, nll);
	h_MMevtype[(*it)+"_estfake"]->AddBinContent(1,nTT);
	h_MMevtype[(*it)+"_estfake"]->AddBinContent(2,nTl);
	h_MMevtype[(*it)+"_estfake"]->AddBinContent(3,nlT);
	h_MMevtype[(*it)+"_estfake"]->AddBinContent(4,nll);
      }
      if(isData)continue;
    
      
      TString truthkey0 = doThreeLep ? truthVECIFF.at(lep1Type) : "trueREAL";
      TString truthkey1 = truthVECIFF.at(lep1Type);
      TString truthkey2 = truthVECIFF.at(lep2Type);
      // June 11th configuration
      TString addkey1 = (*it)+"_";
      TString addkey2 = "";
      if(h_lep_pT_origin_nL.find(((*it)+"_trueREAL")) == h_lep_pT_origin_nL.end()){
	//printf("WARNING \t Did not find %s in h_lep_pT1\n",(*it).Data());
	
	//cout<<"Switching from "<<((*it)+"_trueREAL").Data();
	
	//ALL_E_REAL2L02_trueREAL_ZMET
	TObjArray *tx = ((*it)).Tokenize("_");

	
	if(tx->GetEntries() >= 2){

	  addkey1 = "";
	  
	  //ALL_E_trueREAL_REAL2L02_ZMET
	  int i_idx = 0;
	  for(int i = 0; i<tx->GetEntries(); i++){
	    if((((TObjString *)tx->At(i))->String()).EqualTo(LT_postfix) || (((TObjString *)tx->At(i))->String()).Contains("REAL")
	       || (((((TObjString *)tx->At(i))->String()).Contains("_M_") || (((TObjString *)tx->At(i))->String()).Contains("_E_")) && (((TObjString *)tx->At(i))->String()).Contains("FAKE")))break;
	    addkey1 += (((TObjString *)tx->At(i))->String())+"_";
	    i_idx += 1;
	  }
	  addkey2 = "";
	  
	  for(int i = (i_idx); i<tx->GetEntries(); i++){
	    addkey2 += "_"+(((TObjString *)tx->At(i))->String());
	  }

	  
	}
	delete tx;
	
	//cout<<" to "<<addkey1.Data()<<"trueREAL"<<addkey2.Data()<<endl;;
    
    
	//continue;
      }

      if(doThreeLep && (!truthkey0.EqualTo("ChargeFlipIsoElectron") && !truthkey0.EqualTo("IsoElectron") && !truthkey0.EqualTo("PromptMuon") && !truthkey0.Contains("Unknown")))truthkey0 = (addkey1+"trueFAKE"+addkey2);
      else if(doThreeLep && !truthkey0.EqualTo("ChargeFlipIsoElectron")) truthkey0 = (addkey1+"trueREAL"+addkey2);
      else if(doThreeLep) truthkey0 = (addkey1+"trueCF"+addkey2);

      if(!truthkey1.EqualTo("ChargeFlipIsoElectron") && !truthkey1.EqualTo("IsoElectron") && !truthkey1.EqualTo("PromptMuon") && !truthkey1.Contains("Unknown"))truthkey1 = (addkey1+"trueFAKE"+addkey2);
      else if(!truthkey1.EqualTo("ChargeFlipIsoElectron")) truthkey1 = (addkey1+"trueREAL"+addkey2);
      else truthkey1 = (addkey1+"trueCF"+addkey2);

      if(!truthkey2.EqualTo("ChargeFlipIsoElectron") && !truthkey2.EqualTo("IsoElectron") && !truthkey2.EqualTo("PromptMuon") && !truthkey2.Contains("Unknown"))truthkey2 = (addkey1+"trueFAKE"+addkey2);
      else if(!truthkey2.EqualTo("ChargeFlipIsoElectron")) truthkey2 = (addkey1+"trueREAL"+addkey2);
      else truthkey2 = (addkey1+"trueCF"+addkey2);

 
      // if(!truthkey0.Contains("REAL") && !truthkey0.Contains("CF"))truthkey0 = ((*it)+"_trueFAKE");
      // else if(truthkey0.Contains("REAL") || truthkey0.Contains("CF"))truthkey0 = ((*it)+"_trueREAL");
      // if(!truthkey1.Contains("REAL") && !truthkey1.Contains("CF"))truthkey1 = ((*it)+"_trueFAKE");
      // else if(truthkey1.Contains("REAL") || truthkey1.Contains("CF"))truthkey1 = ((*it)+"_trueREAL");
      // if(!truthkey2.Contains("REAL") && !truthkey2.Contains("CF"))truthkey2 = ((*it)+"_trueFAKE");
      // else if(truthkey2.Contains("REAL") || truthkey2.Contains("CF"))truthkey2 = ((*it)+"_trueREAL");
      // June 8th configuration
      //if(!truthkey1.Contains("REAL") && !truthkey1.Contains("UNKNOWN"))truthkey1 = ((*it)+"_trueFAKE");
      //else if(truthkey1.Contains("REAL"))truthkey1 = ((*it)+"_trueREAL");

      //if(!truthkey2.Contains("REAL") && !truthkey2.Contains("UNKNOWN"))truthkey2 = ((*it)+"_trueFAKE");
      //else if(truthkey2.Contains("REAL"))truthkey2 = ((*it)+"_trueREAL");
      // if(lep1isTight && lep2isTight){
      //   //cout<<"3) "<<truthkey1.Data()<<endl;
      //   if(doThreeLep)h_lep_pT0[(*it)]->Fill(lepPt[idx0],wgt);
      //   h_lep_pT1[(*it)]->Fill(lepPt[idx1],wgt);
      //   h_lep_pT2[(*it)]->Fill(lepPt[idx2],wgt);
      //   h_lep_mt2[(*it)]->Fill(themt2,wgt);
      //   h_lep_mtw[(*it)]->Fill(mtw,wgt);
      //   h_lep_met[(*it)]->Fill(*met_Et,wgt);
      //   h_lep_mll[(*it)]->Fill(themll,wgt);
      // }
      if(truthkey0.Contains("trueREAL") && doThreeLep)fillTruthInfo(truthkey0,idx0,lep0isTight,wgt);
      if(truthkey1.Contains("trueREAL")              )fillTruthInfo(truthkey1,idx1,lep1isTight,wgt);
      if(truthkey2.Contains("trueREAL")              )fillTruthInfo(truthkey2,idx2,lep2isTight,wgt);

      // cout<<"truthkey1 = "<<truthkey1.Data()<<endl;
      // cout<<"truthkey2 = "<<truthkey2.Data()<<endl;
      /**
      cout<<"4) "<<truthkey1.Data()<<endl;
      printf("lep0isTight = %s\n",lep0isTight ? "True" : "False");
      printf("lep1isTight = %s\n",lep1isTight ? "True" : "False");
      printf("lep2isTight = %s\n",lep2isTight ? "True" : "False");
      */
      if(lep1isTight && lep2isTight && truthkey0.Contains("trueREAL") && truthkey1.Contains("trueREAL") && truthkey2.Contains("trueREAL")){
	if(doThreeLep)h_lep_pT0[truthkey0]->Fill(lepPt[idx0],wgt);
	h_lep_pT1[truthkey1]->Fill(lepPt[idx1],wgt);
	
	h_lep_pT2[truthkey1]->Fill(lepPt[idx2],wgt);
	
	h_lep_mt2[truthkey1]->Fill(themt2,wgt);
	
	h_lep_metsig[truthkey1]->Fill(*met_Sign,wgt);
	
	if(doThreeLep)h_lep_mtw[truthkey1]->Fill(mtw,wgt);
	h_lep_met[truthkey1]->Fill(*met_Et,wgt);
	
	h_lep_mll[truthkey1]->Fill(themll,wgt);
	
	h_lep_njet30[truthkey1]->Fill(*nJet30,wgt);
	
	h_lep_nbj77[truthkey1]->Fill(n_bjet77,wgt);
	
	h_lep_nbj85[truthkey1]->Fill(n_bjet85,wgt);
	
	
      }
      
      h_lep_TCF_nL[(*it)]->AddBinContent(lep1Type+1,wgt);
      
      h_lep_mT2_TCF_nL[(*it)]->Fill(themt2,lep1Type,wgt);
      if(doThreeLep)h_lep_mTW_TCF_nL[(*it)]->Fill(mtw,lep1Type,wgt);
      h_lep_njet_TCF_nL[(*it)]->Fill(n_clj,lep1Type,wgt);
      h_lep_metsig_TCF_nL[(*it)]->Fill(*met_Sign,lep1Type,wgt);
      h_lep_pT_TCF_nL[(*it)]->Fill(lepPt[idx1],lep1Type,wgt);

      
      fillTruthInfo((*it),idx1,lep1isTight,wgt);
      
      if(truthkey1.Contains("trueFAKE")){
	fillTruthInfo(truthkey1,idx1,lep1isTight,wgt);
	h_lep_TCF_nL[truthkey1]->AddBinContent(lep1Type+1,wgt);
	h_lep_mT2_TCF_nL[truthkey1]->Fill(themt2,lep1Type,wgt);
	if(doThreeLep)h_lep_mTW_TCF_nL[truthkey1]->Fill(mtw,lep1Type,wgt);
	h_lep_njet_TCF_nL[truthkey1]->Fill(n_clj,lep1Type,wgt);
	h_lep_metsig_TCF_nL[truthkey1]->Fill(*met_Sign,lep1Type,wgt);
	h_lep_pT_TCF_nL[truthkey1]->Fill(lepPt[idx1],lep1Type,wgt);
	
	      
      }
      if(lep1isTight){
	h_lep_TCF_nT[(*it)]->AddBinContent(lep1Type+1,wgt);
	h_lep_mT2_TCF_nT[(*it)]->Fill(themt2,lep1Type,wgt);
	if(doThreeLep)h_lep_mTW_TCF_nT[(*it)]->Fill(mtw,lep1Type,wgt);
	h_lep_njet_TCF_nT[(*it)]->Fill(n_clj,lep1Type,wgt);
	h_lep_metsig_TCF_nT[(*it)]->Fill(*met_Sign,lep1Type,wgt);
	h_lep_pT_TCF_nT[(*it)]->Fill(lepPt[idx1],lep1Type,wgt);
	if(truthkey1.Contains("trueFAKE")){
	  h_lep_TCF_nT[truthkey1]->AddBinContent(lep1Type+1,wgt);
	  h_lep_mT2_TCF_nT[truthkey1]->Fill(themt2,lep1Type,wgt);
	  if(doThreeLep)h_lep_mTW_TCF_nT[truthkey1]->Fill(mtw,lep1Type,wgt);
	  h_lep_njet_TCF_nT[truthkey1]->Fill(n_clj,lep1Type,wgt);
	  h_lep_metsig_TCF_nT[truthkey1]->Fill(*met_Sign,lep1Type,wgt);
	  h_lep_pT_TCF_nT[truthkey1]->Fill(lepPt[idx1],lep1Type,wgt);
	}
      }
      h_lep_TCF_nL[(*it)]->AddBinContent(lep2Type+1,wgt);
      h_lep_mT2_TCF_nL[(*it)]->Fill(themt2,lep2Type,wgt);
      if(doThreeLep)h_lep_mTW_TCF_nL[(*it)]->Fill(mtw,lep2Type,wgt);
      h_lep_njet_TCF_nL[(*it)]->Fill(n_clj,lep2Type,wgt);
      h_lep_metsig_TCF_nL[(*it)]->Fill(*met_Sign,lep2Type,wgt);
      h_lep_pT_TCF_nL[(*it)]->Fill(lepPt[idx2],lep2Type,wgt);
      fillTruthInfo((*it),idx2,lep2isTight,wgt);
            
      if(truthkey2.Contains("trueFAKE")){
	fillTruthInfo(truthkey2,idx2,lep2isTight,wgt);
	h_lep_TCF_nL[truthkey2]->AddBinContent(lep2Type+1,wgt);
	h_lep_mT2_TCF_nL[truthkey2]->Fill(themt2,lep2Type,wgt);
	if(doThreeLep)h_lep_mTW_TCF_nL[truthkey2]->Fill(mtw,lep2Type,wgt);
	h_lep_njet_TCF_nL[truthkey2]->Fill(n_clj,lep2Type,wgt);
	h_lep_metsig_TCF_nL[truthkey2]->Fill(*met_Sign,lep2Type,wgt);
	h_lep_pT_TCF_nL[truthkey2]->Fill(lepPt[idx2],lep2Type,wgt);
	
	      
      }
      if(lep2isTight){
	h_lep_TCF_nT[(*it)]->AddBinContent(lep2Type+1,wgt);
	h_lep_mT2_TCF_nT[(*it)]->Fill(themt2,lep2Type,wgt);
	if(doThreeLep)h_lep_mTW_TCF_nT[(*it)]->Fill(mtw,lep2Type,wgt);
	h_lep_njet_TCF_nT[(*it)]->Fill(n_clj,lep2Type,wgt);
	h_lep_metsig_TCF_nT[(*it)]->Fill(*met_Sign,lep2Type,wgt);
	h_lep_pT_TCF_nT[(*it)]->Fill(lepPt[idx2],lep2Type,wgt);
	if(truthkey2.Contains("trueFAKE")){
	  h_lep_TCF_nT[truthkey2]->AddBinContent(lep2Type+1,wgt);
	  h_lep_mT2_TCF_nT[truthkey2]->Fill(themt2,lep2Type,wgt);
	  if(doThreeLep)h_lep_mTW_TCF_nT[truthkey2]->Fill(mtw,lep2Type,wgt);
	  h_lep_njet_TCF_nT[truthkey2]->Fill(n_clj,lep2Type,wgt);
	  h_lep_metsig_TCF_nT[truthkey2]->Fill(*met_Sign,lep2Type,wgt);
	  h_lep_pT_TCF_nT[truthkey2]->Fill(lepPt[idx2],lep2Type,wgt);
	}
      }
    }

  }// end loop over LT definitions

  //cout<<"Done with this event"<<endl;

  return kTRUE;
}

void MySusySkimAnalysis::fillTruthInfo(TString& key, int idx, bool lepistight, double wgt)
{

  if(isData)return;

  if(h_lep_pT_origin_nL.find(key) == h_lep_pT_origin_nL.end()){
    //printf("WARNING \t Did not find %s in h_lep_pT1\n",(*it).Data());

    //cout<<"Switching from "<<key.Data();

    //ALL_E_REAL2L02_trueREAL_ZMET
    TObjArray *tx = key.Tokenize("_");


    //ALL_E_trueREAL_REAL2L02_ZMET
  
    key = (((TObjString *)tx->At(0))->String())+"_"+(((TObjString *)tx->At(1))->String())+"_"+(((TObjString *)tx->At(3))->String())+"_"+(((TObjString *)tx->At(2))->String());

    for(int i = 4; i<tx->GetEntries(); i++){
      key += "_"+(((TObjString *)tx->At(i))->String());
    }

    delete tx;
    
    //    cout<<" to "<<key.Data()<<endl;;

    
    
    //continue;
  }

    
  //cout<<"key = "<<key<<endl;
  h_lep_pT_origin_nL[key]->Fill(lepPt[idx],lepOrigin[idx],wgt);
  if(lepistight)h_lep_pT_origin_nT[key]->Fill(lepPt[idx],lepOrigin[idx],wgt);
  //cout<<"a"<<endl;

  h_lep_pT_type_nL[key]->Fill(lepPt[idx],lepType[idx],wgt);
  if(lepistight)h_lep_pT_type_nT[key]->Fill(lepPt[idx],lepType[idx],wgt);
  //cout<<"b"<<endl;

  h_lep_pT_firstEgMotherO_nL[key]->Fill(lepPt[idx],lepEgMotherOrigin[idx],wgt);
  if(lepistight)h_lep_pT_firstEgMotherO_nT[key]->Fill(lepPt[idx],lepEgMotherOrigin[idx],wgt);
  //cout<<"c"<<endl;

  h_lep_pT_firstEgMotherT_nL[key]->Fill(lepPt[idx],lepEgMotherType[idx],wgt);
  if(lepistight)h_lep_pT_firstEgMotherT_nT[key]->Fill(lepPt[idx],lepEgMotherType[idx],wgt);
  //cout<<"d"<<endl;

  h_lep_pT_firstEgMotherPdgId_nL[key]->Fill(lepPt[idx],lepEgMotherOrigin[idx],wgt);
  if(lepistight)h_lep_pT_firstEgMotherPdgId_nT[key]->Fill(lepPt[idx],lepEgMotherOrigin[idx],wgt);
  //cout<<"e"<<endl;

  h_lep_type_origin_nL[key]->Fill(lepType[idx],lepOrigin[idx],wgt);
  if(lepistight)h_lep_type_origin_nT[key]->Fill(lepType[idx],lepOrigin[idx],wgt);
  //cout<<"f"<<endl;

  
}



void MySusySkimAnalysis::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void MySusySkimAnalysis::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

  
  cout<<"Processed "<<ev_was_found<<" events, now finishing"<<endl;
  cout<<"FNPest = "<<FNPest<<endl;
  cout<<"FNPest_up = "<<FNPest_up<<endl;
  cout<<"FNPest_dw = "<<FNPest_dw<<endl;
  
  if(makeNtup)delete HFT;
  cout<<"1"<<endl;
  if(makeNtup && make2L2JNtup)delete twoLtwoJ;
  cout<<"2"<<endl; 
  if(makeNtup && makeMYNtup)delete MY;
  cout<<"3"<<endl;
  if(!makeHistograms)return;
  cout<<"4"<<endl;
  
  WriteToFile();
  cout<<"5"<<endl;
  
  printCutflow();
  cout<<"6"<<endl;
  
  /**
  for (unsigned int i=0; i<DSIDcheck.size(); i++){
    cout<<"DSID : "<<DSIDcheck[i]<<endl;
  }
  */

  //return;
  
}

void MySusySkimAnalysis::fillLandT(int idx, std::vector<TString> regions, TString ananame, bool isFake, bool trigger, TString tname)
{
  //if(doAnalysis)return;
  if(!trigger)return;
  std::vector< TString >::iterator it;
  std::vector< TString >::iterator it2;
  std::vector< TString >::iterator it3;
  std::vector< TString >::iterator it4;
  std::map< TString, std::vector<TTreeReaderValue<Float_t>> >::iterator  mapit;
  bool lepistriggermatched = false;
  TString lep = "";
  TString fakeORreal = isFake ? "FAKE" : "REAL";
  bool lepIsTight = isT(idx);
  TString unckey = "";
  if(lepFlavor[idx] == 1)lep = "E";
  else if(lepFlavor[idx] == 2)lep = "M";
  else printf("ERROR /t Flavour %i of lepton is unknown!!\n",lepFlavor[idx]);
  double old_wgt = wgt;
  double old_wgt_loose = wgt_loose;
  float scalelumi = 1.0;
  if(!isData)scalelumi = (*RandomRunNumber) < 320000 ? 36207.65 : (((*RandomRunNumber) > 320000 && (*RandomRunNumber) < 348000) ? 44307.4 : 58450.1);
  //cout<<"-------------------"<<*EventNumber<<"-----------------------------"<<endl;
  for(mapit = uncert.begin(); mapit != uncert.end(); ++mapit){
    wgt = old_wgt;
    wgt_loose = old_wgt_loose;
    unckey = "";
    if(!(mapit->first).EqualTo("NOMINAL")){
      unckey = "_"+mapit->first;
      if(unckey.Contains("Trig")){
	wgt = ((*genWeight) * (*eventWeight) * (*leptonWeight) * (*jvtWeight) * (*bTagWeight) * (*pileupWeight) * scalelumi * extrafac * (*((mapit->second).at(0)))/** (*FFWeight)*/);
	//printf("INFO \t Uncertainty %s - changing weight from %.6f to %.6f \n",unckey.Data(),old_wgt,wgt);
      }else if(unckey.Contains("MUON") || unckey.Contains("EL")){
	wgt = ((*genWeight) * (*eventWeight) * (*((mapit->second).at(0))) * (*jvtWeight) * (*bTagWeight) * (*pileupWeight) * scalelumi * extrafac * (*globalDiLepTrigSF)/** (*FFWeight)*/);
	//printf("INFO \t Uncertainty %s - changing weight from %.6f to %.6f \n",unckey.Data(),old_wgt,wgt);
      }else{
	printf("ERROR \t Could not identify class of uncertainty %s \n",unckey.Data());
	continue;
      }
    }else if(!isData){
      wgt       = ((lepRecoSF[idx])   * (*genWeight) * (*eventWeight) * (*jvtWeight) * (*bTagWeight) * (*pileupWeight) * scalelumi * extrafac * (ananame.Contains("2L2J") ? (*globalDiLepTrigSF) : lepTrigSF[idx]));
      wgt_loose = ((lepBLRecoSF[idx]) * (*genWeight) * (*eventWeight) * (*jvtWeight) * (*bTagWeight) * (*pileupWeight) * scalelumi * extrafac * (ananame.Contains("2L2J") ? (*globalBaselineDiLepTrigSF) : lepBLTrigSF[idx]));
      if(makeHistograms){
	h_wgts_infilllandt_wgt->Fill(wgt);
	h_wgts_infilllandt_wgtloose->Fill(wgt_loose);
      }
      // printf("::fillLandT Wgt    = %.8f\n",wgt);
      // printf("::fillLandT Wgt BL = %.8f\n",wgt_loose);
    }
     
    //cout<<"Uncertainty key is "<<unckey<<endl; 
    for(it = fillregions.begin(); it != fillregions.end(); it++){
      //cout<<"Looking at "<<(*it).Data()<<endl;
      if(isFake && ((*it).Contains("OS") /**|| (*it).Contains("SS")*/))continue; // skip doing the sign dep. since this is already defined in the fake/real regions
      for(it2 = regions.begin(); it2 != regions.end(); it2++){
	TString key = (*it2).BeginsWith("2L") ? (*it)+"_"+lep+"_"+fakeORreal+(*it2) : (*it)+"_"+lep+"_"+(*it2);
	TString add = (*it2);
	if(tname.Length() > 0){
	  key += "_"+tname;
	  add += "_"+tname;
	}
	key += unckey+"_"+LT_postfix;
	add += unckey+"_"+LT_postfix;
	//printf("INFO \t Key %s - changing weight from %.6f to %.6f \n",key.Data(),old_wgt,wgt);
	if(h_lep_pT_nL.find(key) == h_lep_pT_nL.end()){
	  // cout<<"Could not find key "<<key<<endl;
	  continue;
	}
	//cout<<"key is "<<key.Data()<<" and wgt_loose = " <<wgt_loose << " and wgt = "<< wgt<<endl;
	if(makeHistograms){
	  h_lep_njet30_nL[key]->Fill(*nJet30,wgt_loose);
	  h_lep_MET_nL[key]->Fill(*met_Et,wgt_loose);
	  h_lep_nbjet_nL[key]->Fill(n_bjet,wgt_loose);
	}
	for (TString bdt_name : bdt_vec) {
	  if(makeHistograms)h_lep_BDT_nL[bdt_name+"_"+key]->Fill(BDTweight[bdt_name],wgt_loose);
	}
	//cout<<"fillLandT :: Filling h_lep_pT_nL for "<< key << endl;
	if(makeHistograms){
	  h_lep_pT_nL[key]->Fill(lepPt[idx],wgt_loose);
	  h_lep_eta_nL[key]->Fill(lepEta[idx],wgt_loose);
	  h_lep_pT_eta_nL[key]->Fill(lepPt[idx],lepEta[idx],wgt_loose);
	  h_lep_mu_pT_nL[key]->Fill(*mu,lepPt[idx],wgt_loose);
	}
	if(triglist.size()>0){
	  //if(*trigMatch_1L2LTrig)
	  checkTrigMatch("1L2LTrig",idx,key,lepIsTight,(lepFlavor[idx] == 1 ? "E" : "M"));
	  //if(*trigMatch_1LTrig)
	  checkTrigMatch("1LTrig",idx,key,lepIsTight,(lepFlavor[idx] == 1 ? "1E" : "1M")); //<------------------------------------- COMMENT BACK IN. NOT AVAILABLE IN 2L2J OFFICIIAL NTUPLES
	  //if(*trigMatch_2LTrigOR)
	  checkTrigMatch("2LTrig",idx,key,lepIsTight,(lepFlavor[idx] == 1 ? "2E" : "2M"));

	  
	  std::map< TString, TString >::iterator  mapit;
	  for(mapit = trigcat.begin(); mapit != trigcat.end(); ++mapit){
	    checkTrigMatch(mapit->first,idx,key,lepIsTight,mapit->second);
	  }
	  
	  // for(unsigned int i=0; i<trigvallep.size();i++){
	  //   if(triglist.at(i).Contains("Trig"))continue;
	  //   //if(*trigval.at(i)){
	  //   Float_t pt_th = getPtThresholdTrigger(triglist.at(i));
	    
	  //   if(trigvallep.at(i)->at(idx) && *trigval.at(i) && lepPt[idx] >= pt_th &&
	  //      ((lepFlavor[idx] == 1 && (std::count(trigstr_e[yr].begin(), trigstr_e[yr].end(), triglist.at(i)))) ||
	  // 	(lepFlavor[idx] == 2 && (std::count(trigstr_m[yr].begin(), trigstr_m[yr].end(), triglist.at(i)))))){
	      
	  //     //if((lepFlavor[idx] == 1 && !(std::count(trigstr_e[yr].begin(), trigstr_e[yr].end(), triglist.at(i)))) ||
	  //     //	 (lepFlavor[idx] == 2 && !(std::count(trigstr_m[yr].begin(), trigstr_m[yr].end(), triglist.at(i)))))continue;
	  //     //Float_t pt_th = getPtThresholdTrigger(triglist.at(i));
	  //     //if(lepPt[idx] < pt_th)continue;
	  //     lepistriggermatched = true;
	  //     for (TString bdt_name : bdt_vec) {
	  // 	//if(key.Contains("EEOS") && key.Contains("REAL2L01"))
	  // 	//  cout<<"Filling "<<bdt_name+"_"+triglist.at(i)+"_"+key<<endl;
	  // 	h_lep_BDT_nL[bdt_name+"_"+triglist.at(i)+"_lepmatched_"+key]->Fill(BDTweight[bdt_name],wgt);
	  //     }
	  //     h_lep_pT_nL[triglist.at(i)+"_lepmatched_"+key]->Fill(lepPt[idx],wgt);
	  //     h_lep_eta_nL[triglist.at(i)+"_lepmatched_"+key]->Fill(lepEta[idx],wgt);
	  //     h_lep_pT_eta_nL[triglist.at(i)+"_lepmatched_"+key]->Fill(lepPt[idx],lepEta[idx],wgt);
	  //     h_lep_trig_nL[key]->AddBinContent(i+1,wgt);
	  //     if(lepIsTight){
	  // 	for (TString bdt_name : bdt_vec) {
	  // 	  h_lep_BDT_nT[bdt_name+"_"+triglist.at(i)+"_lepmatched_"+key]->Fill(BDTweight[bdt_name],wgt);
	  // 	}
	  // 	h_lep_trig_nT[key]->AddBinContent(i+1,wgt);
	  // 	h_lep_pT_nT[triglist.at(i)+"_lepmatched_"+key]->Fill(lepPt[idx],wgt);
	  // 	h_lep_eta_nT[triglist.at(i)+"_lepmatched_"+key]->Fill(lepEta[idx],wgt);
	  // 	h_lep_pT_eta_nT[triglist.at(i)+"_lepmatched_"+key]->Fill(lepPt[idx],lepEta[idx],wgt);
	  //     }
	  //   }
	  // }
	  
	  for(unsigned int i=0; i<trigval.size();i++){
	    if(triglist.at(i).Contains("Trig"))continue;

	    Float_t pt_th = getPtThresholdTrigger(triglist.at(i),lepFlavor[idx] == 1 ? "e" : "mu");
	    
	    if(trigvallep.at(i)->at(idx) && *trigval.at(i) && lepPt[idx] >= pt_th &&
	       ((lepFlavor[idx] == 1 && (std::count(trigstr_e[yr].begin(), trigstr_e[yr].end(), triglist.at(i)))) ||
		(lepFlavor[idx] == 2 && (std::count(trigstr_m[yr].begin(), trigstr_m[yr].end(), triglist.at(i)))))){

	      //printf("Key %s : Event triggered by %s and %s %i (pT = %.2f) matched to it\n",(triglist.at(i)+"_lepmatched_"+key).Data(),triglist.at(i).Data(),lepFlavor[idx] == 1 ? "electron" : "muon",idx,lepPt[idx]);
	      
	      for (TString bdt_name : bdt_vec) {
		if(makeHistograms)h_lep_BDT_nL[bdt_name+"_"+triglist.at(i)+"_lepmatched_"+key]->Fill(BDTweight[bdt_name],wgt_loose);
	      }
	      if(makeHistograms){
		h_lep_pT_nL[triglist.at(i)+"_lepmatched_"+key]->Fill(lepPt[idx],wgt_loose);
		h_lep_eta_nL[triglist.at(i)+"_lepmatched_"+key]->Fill(lepEta[idx],wgt_loose);
		h_lep_pT_eta_nL[triglist.at(i)+"_lepmatched_"+key]->Fill(lepPt[idx],lepEta[idx],wgt_loose);
		h_lep_trig_nL["lepmatched_"+key]->AddBinContent(i+1,wgt_loose);
	      }
	      if(lepIsTight){
		for (TString bdt_name : bdt_vec) {
		  if(makeHistograms)h_lep_BDT_nT[bdt_name+"_"+triglist.at(i)+"_lepmatched_"+key]->Fill(BDTweight[bdt_name],wgt);
		}
		if(makeHistograms){
		  h_lep_trig_nT["lepmatched_"+key]->AddBinContent(i+1,wgt);
		  h_lep_pT_nT[triglist.at(i)+"_lepmatched_"+key]->Fill(lepPt[idx],wgt);
		  h_lep_eta_nT[triglist.at(i)+"_lepmatched_"+key]->Fill(lepEta[idx],wgt);
		  h_lep_pT_eta_nT[triglist.at(i)+"_lepmatched_"+key]->Fill(lepPt[idx],lepEta[idx],wgt);
		}
	      }
	    }else if(*trigval.at(i)){
	      //printf("Key %s : Event triggered by %s, but %s %i (pT = %.2f) NOT matched to it\n",(triglist.at(i)+"_lepnotmatched_"+key).Data(),triglist.at(i).Data(),lepFlavor[idx] == 1 ? "electron" : "muon",idx,lepPt[idx]);
	      for (TString bdt_name : bdt_vec) {
		if(makeHistograms)h_lep_BDT_nL[bdt_name+"_"+triglist.at(i)+"_lepnotmatched_"+key]->Fill(BDTweight[bdt_name],wgt_loose);
	      }
	      if(makeHistograms){
	      h_lep_pT_nL[triglist.at(i)+"_lepnotmatched_"+key]->Fill(lepPt[idx],wgt_loose);
	      h_lep_eta_nL[triglist.at(i)+"_lepnotmatched_"+key]->Fill(lepEta[idx],wgt_loose);
	      h_lep_pT_eta_nL[triglist.at(i)+"_lepnotmatched_"+key]->Fill(lepPt[idx],lepEta[idx],wgt_loose);
	      h_lep_trig_nL["lepnotmatched_"+key]->AddBinContent(i+1,wgt_loose);
	      }
	      if(lepIsTight){
		for (TString bdt_name : bdt_vec) {
		  if(makeHistograms)h_lep_BDT_nT[bdt_name+"_"+triglist.at(i)+"_lepnotmatched_"+key]->Fill(BDTweight[bdt_name],wgt);
		}
		if(makeHistograms){
		h_lep_pT_nT[triglist.at(i)+"_lepnotmatched_"+key]->Fill(lepPt[idx],wgt);
		h_lep_eta_nT[triglist.at(i)+"_lepnotmatched_"+key]->Fill(lepEta[idx],wgt);
		h_lep_pT_eta_nT[triglist.at(i)+"_lepnotmatched_"+key]->Fill(lepPt[idx],lepEta[idx],wgt);
		h_lep_trig_nT["lepnotmatched_"+key]->AddBinContent(i+1,wgt);
		}
	      }
	    }
	  }
	}
	
	//       if(!lepistriggermatched){
	// 	for (TString bdt_name : bdt_vec) {
	// 	  h_lep_BDT_nL[bdt_name+"_"+triglist.at(i)+"_lepnotmatched_"+key]->Fill(BDTweight[bdt_name],wgt);
	// 	}
	// 	h_lep_pT_nL[triglist.at(i)+"_lepnotmatched_"+key]->Fill(lepPt[idx],wgt);
	// 	h_lep_eta_nL[triglist.at(i)+"_lepnotmatched_"+key]->Fill(lepEta[idx],wgt);
	// 	h_lep_pT_eta_nL[triglist.at(i)+"_lepnotmatched_"+key]->Fill(lepPt[idx],lepEta[idx],wgt);
	//       }
	//       if(lepIsTight){
	// 	for (TString bdt_name : bdt_vec) {
	// 	  h_lep_BDT_nT[bdt_name+"_"+triglist.at(i)+"_eventTrig_"+key]->Fill(BDTweight[bdt_name],wgt);
	// 	}
	// 	h_lep_trig_nT["eventTrig_"+key]->AddBinContent(i+1,wgt);
	// 	h_lep_pT_nT[triglist.at(i)+"_eventTrig_"+key]->Fill(lepPt[idx],wgt);
	// 	h_lep_eta_nT[triglist.at(i)+"_eventTrig_"+key]->Fill(lepEta[idx],wgt);
	// 	h_lep_pT_eta_nT[triglist.at(i)+"_eventTrig_"+key]->Fill(lepPt[idx],lepEta[idx],wgt);
	// 	if(!lepistriggermatched){
	// 	  for (TString bdt_name : bdt_vec) {
	// 	    h_lep_BDT_nT[bdt_name+"_"+triglist.at(i)+"_lepnotmatched_"+key]->Fill(BDTweight[bdt_name],wgt);
	// 	  }
	// 	  h_lep_pT_nT[triglist.at(i)+"_lepnotmatched_"+key]->Fill(lepPt[idx],wgt);
	// 	  h_lep_eta_nT[triglist.at(i)+"_lepnotmatched_"+key]->Fill(lepEta[idx],wgt);
	// 	  h_lep_pT_eta_nT[triglist.at(i)+"_lepnotmatched_"+key]->Fill(lepPt[idx],lepEta[idx],wgt);
	// 	}
	      
	//       }
	//     }
	//   }	
	// }
	fillFakeRateHist(key, idx, lepIsTight);
	if(lepIsTight){
	  if(makeHistograms)h_lep_njet30_nT[key]->Fill(*nJet30,wgt);
	  for (TString bdt_name : bdt_vec) {
	    if(makeHistograms)h_lep_BDT_nT[bdt_name+"_"+key]->Fill(BDTweight[bdt_name],wgt);
	  }
	  if(makeHistograms)h_lep_pT_nT[key]->Fill(lepPt[idx],wgt);
	  if(makeHistograms)h_lep_eta_nT[key]->Fill(lepEta[idx],wgt);
	  if(makeHistograms)h_lep_pT_eta_nT[key]->Fill(lepPt[idx],lepEta[idx],wgt);
	  if(makeHistograms)h_lep_mu_pT_nT[key]->Fill(*mu,lepPt[idx],wgt);
	  if(makeHistograms)h_lep_MET_nT[key]->Fill(*met_Et,wgt);
	  if(makeHistograms)h_lep_nbjet_nT[key]->Fill(n_bjet,wgt);
	}
	if(isData)continue;
      

	std::vector<TString> truthkey;
	truthkey.clear();
	truthENUMIFF tcf = (truthENUMIFF)lepIFFClass[idx];//TruthClassification(lepType[idx], lepOrigin[idx], lepEgMotherType[idx], lepEgMotherOrigin[idx], lepEgMotherPdgId[idx], lepCharge[idx], lepFlavor[idx]);
	if(makeHistograms)h_lep_TCF_nL[key]->AddBinContent(tcf+1,wgt_loose);
	if(makeHistograms)h_lep_mT2_TCF_nL[key]->Fill(themt2,tcf,wgt_loose);
	if(doThreeLep)if(makeHistograms)h_lep_mTW_TCF_nL[key]->Fill(mtw,tcf,wgt_loose);
	if(makeHistograms)h_lep_njet_TCF_nL[key]->Fill(n_clj,tcf,wgt_loose);
	if(makeHistograms)h_lep_metsig_TCF_nL[key]->Fill(*met_Sign,tcf,wgt_loose);
	if(makeHistograms)h_lep_pT_TCF_nL[key]->Fill(lepPt[idx],tcf,wgt_loose);
	if(isT(idx)){
	  if(makeHistograms)h_lep_TCF_nT[key]->AddBinContent(tcf+1,wgt);
	  if(makeHistograms)h_lep_mT2_TCF_nT[key]->Fill(themt2,tcf,wgt);
	  if(doThreeLep)if(makeHistograms)h_lep_mTW_TCF_nT[key]->Fill(mtw,tcf,wgt);
	  if(makeHistograms)h_lep_njet_TCF_nT[key]->Fill(n_clj,tcf,wgt);
	  if(makeHistograms)h_lep_metsig_TCF_nT[key]->Fill(*met_Sign,tcf,wgt);
	  if(makeHistograms)h_lep_pT_TCF_nT[key]->Fill(lepPt[idx],tcf,wgt);
	}
	//cout<<"g"<<endl;
	TString lep2Type = truthVECIFF.at(tcf);
	//truthkey.push_back((*it)+"_"+lep+"_true"+lep2Type+"_"+fakeORreal+(*it2));
	//if(!doCoreRegions || lep2Type.Contains("REAL") || lep2Type.Contains("LF")){
	truthkey.push_back((*it)+"_"+lep+"_true"+lep2Type+"_"+fakeORreal+add);
	//}
	//if(!lep2Type.Contains("REAL") && !lep2Type.Contains("CF"))truthkey.push_back((*it)+"_"+lep+"_trueFAKE_"+fakeORreal+(*it2));
      //      if(!lep2Type.Contains("REAL") && !lep2Type.Contains("CF"))truthkey.push_back((*it)+"_"+lep+"_trueFAKE_"+fakeORreal+add);

      if(!lep2Type.EqualTo("ChargeFlipIsoElectron") && !lep2Type.EqualTo("IsoElectron") && !lep2Type.EqualTo("PromptMuon") && !lep2Type.Contains("Unknown"))truthkey.push_back((*it)+"_"+lep+"_trueFAKE_"+fakeORreal+add);
      else if(!lep2Type.EqualTo("ChargeFlipIsoElectron")) truthkey.push_back((*it)+"_"+lep+"_trueREAL_"+fakeORreal+add);
      else truthkey.push_back((*it)+"_"+lep+"_trueCF_"+fakeORreal+add);

	
      
      for(it3 = truthkey.begin(); it3 != truthkey.end(); it3++){

	if(h_lep_pT_nL.find((*it3)) == h_lep_pT_nL.end()){
	  //printf("WARNING \t Could not find key h_lep_pT_nL_%s\n",(*it3).Data());
	  continue;
	}
	// Fill trigger also for truth MC
	if(triglist.size()>0){
	  //if(*trigMatch_1L2LTrig)
	  checkTrigMatch("1L2LTrig",idx,(*it3),lepIsTight,(lepFlavor[idx] == 1 ? "E" : "M"));
	  //if(*trigMatch_1LTrig)
	  checkTrigMatch("1LTrig",idx,(*it3),lepIsTight,(lepFlavor[idx] == 1 ? "1E" : "1M")); //<------------------------------------- COMMENT BACK IN. NOT AVAILABLE IN 2L2J OFFICIIAL NTUPLES
	  //if(*trigMatch_2LTrigOR)
	  checkTrigMatch("2LTrig",idx,(*it3),lepIsTight,(lepFlavor[idx] == 1 ? "2E" : "2M"));


	  std::map< TString, TString >::iterator  mapit;
	  
	  for(mapit = trigcat.begin(); mapit != trigcat.end(); ++mapit){
	    checkTrigMatch(mapit->first,idx,(*it3),lepIsTight,mapit->second);
	  }

	  /**
	  for(it4 = categories.begin(); it4 != categories.end(); it4++){
	    checkTrigMatch((*it4),idx,(*it3),lepIsTight,(*it4));
	  }
	  */

	  lepistriggermatched = false;
	  for(unsigned int i=0; i<trigvallep.size();i++){
	    if(triglist.at(i).Contains("Trig"))continue;

	    Float_t pt_th = getPtThresholdTrigger(triglist.at(i),lepFlavor[idx] == 1 ? "e" : "mu");

	    // Trigger is matched and trigger has fired
	    if(trigvallep.at(i)->at(idx) && *trigval.at(i) && lepPt[idx] >= pt_th &&
	       ((lepFlavor[idx] == 1 && (std::count(trigstr_e[yr].begin(), trigstr_e[yr].end(), triglist.at(i)))) ||
		(lepFlavor[idx] == 2 && (std::count(trigstr_m[yr].begin(), trigstr_m[yr].end(), triglist.at(i)))))){
    
	      //f((lepFlavor[idx] == 1 && !(std::count(trigstr_e[yr].begin(), trigstr_e[yr].end(), triglist.at(i)))) ||
	      //	 (lepFlavor[idx] == 2 && !(std::count(trigstr_m[yr].begin(), trigstr_m[yr].end(), triglist.at(i)))))continue;
	      
	      //if(lepPt[idx] < pt_th)continue;
	      lepistriggermatched = true;
	      //printf("Key %s : Event triggered by %s and %s %i (pT = %.2f) matched to it\n",(triglist.at(i)+"_lepmatched_"+(*it3)).Data(),triglist.at(i).Data(),lepFlavor[idx] == 1 ? "electron" : "muon",idx,lepPt[idx]);
	      for (TString bdt_name : bdt_vec) {
		h_lep_BDT_nL[bdt_name+"_"+triglist.at(i)+"_lepmatched_"+(*it3)]->Fill(BDTweight[bdt_name],wgt_loose);
	      }
	      h_lep_pT_nL[triglist.at(i)+"_lepmatched_"+(*it3)]->Fill(lepPt[idx],wgt_loose);
	      h_lep_eta_nL[triglist.at(i)+"_lepmatched_"+(*it3)]->Fill(lepEta[idx],wgt_loose);
	      h_lep_pT_eta_nL[triglist.at(i)+"_lepmatched_"+(*it3)]->Fill(lepPt[idx],lepEta[idx],wgt_loose);
	      h_lep_trig_nL["lepmatched_"+(*it3)]->AddBinContent(i+1,wgt_loose);
	      if(lepIsTight){
		for (TString bdt_name : bdt_vec) {
		  h_lep_BDT_nT[bdt_name+"_lepmatched_"+triglist.at(i)+"_"+(*it3)]->Fill(BDTweight[bdt_name],wgt);
		}
		h_lep_trig_nT["lepmatched_"+(*it3)]->AddBinContent(i+1,wgt);
		h_lep_pT_nT[triglist.at(i)+"_lepmatched_"+(*it3)]->Fill(lepPt[idx],wgt);
		h_lep_eta_nT[triglist.at(i)+"_lepmatched_"+(*it3)]->Fill(lepEta[idx],wgt);
		h_lep_pT_eta_nT[triglist.at(i)+"_lepmatched_"+(*it3)]->Fill(lepPt[idx],lepEta[idx],wgt);
	      }
	      // Trigger has fired, but not matched to this lepton
	    }else if(*trigval.at(i)){
	      //printf("Key %s : Event triggered by %s, but %s %i (pT = %.2f) NOT matched to it\n",(triglist.at(i)+"_lepnotmatched_"+(*it3)).Data(),triglist.at(i).Data(),lepFlavor[idx] == 1 ? "electron" : "muon",idx,lepPt[idx]);
	      for (TString bdt_name : bdt_vec) {		 
		h_lep_BDT_nL[bdt_name+"_"+triglist.at(i)+"_lepnotmatched_"+(*it3)]->Fill(BDTweight[bdt_name],wgt_loose);
	      }
	      h_lep_trig_nL["lepnotmatched_"+(*it3)]->AddBinContent(i+1,wgt_loose);
	      h_lep_pT_nL[triglist.at(i)+"_lepnotmatched_"+(*it3)]->Fill(lepPt[idx],wgt_loose);
	      h_lep_eta_nL[triglist.at(i)+"_lepnotmatched_"+(*it3)]->Fill(lepEta[idx],wgt_loose);
	      h_lep_pT_eta_nL[triglist.at(i)+"_lepnotmatched_"+(*it3)]->Fill(lepPt[idx],lepEta[idx],wgt_loose);
	      if(lepIsTight){
		for (TString bdt_name : bdt_vec) {		 
		  h_lep_BDT_nT[bdt_name+"_"+triglist.at(i)+"_lepnotmatched_"+(*it3)]->Fill(BDTweight[bdt_name],wgt);
		}
		h_lep_trig_nT["lepnotmatched_"+(*it3)]->AddBinContent(i+1,wgt);
		h_lep_pT_nT[triglist.at(i)+"_lepnotmatched_"+(*it3)]->Fill(lepPt[idx],wgt);
		h_lep_eta_nT[triglist.at(i)+"_lepnotmatched_"+(*it3)]->Fill(lepEta[idx],wgt);
		h_lep_pT_eta_nT[triglist.at(i)+"_lepnotmatched_"+(*it3)]->Fill(lepPt[idx],lepEta[idx],wgt);
		
	      }
	    }
	  }
	}
      
	  /**
	  for(unsigned int i=0; i<trigval.size();i++){
	    if(triglist.at(i).Contains("Trig"))continue;
	    if((!(std::count(trigstr_e[yr].begin(), trigstr_e[yr].end(), triglist.at(i)))) && !(std::count(trigstr_m[yr].begin(), trigstr_m[yr].end(), triglist.at(i))))continue;
	    //if(!trigvallep.at(i)->at(idx)){
	    if(*trigval.at(i)){
	      for (TString bdt_name : bdt_vec) {		 
		h_lep_BDT_nL[bdt_name+"_"+triglist.at(i)+"_eventTrig_"+(*it3)]->Fill(BDTweight[bdt_name],wgt);
	      }
	      h_lep_pT_nL[triglist.at(i)+"_eventTrig_"+(*it3)]->Fill(lepPt[idx],wgt);
	      h_lep_eta_nL[triglist.at(i)+"_eventTrig_"+(*it3)]->Fill(lepEta[idx],wgt);
	      h_lep_pT_eta_nL[triglist.at(i)+"_eventTrig_"+(*it3)]->Fill(lepPt[idx],lepEta[idx],wgt);
	      h_lep_trig_nL["eventTrig_"+(*it3)]->AddBinContent(i+1,wgt);
	      if(!lepistriggermatched){
		for (TString bdt_name : bdt_vec) {		 
		  h_lep_BDT_nL[bdt_name+"_"+triglist.at(i)+"_lepnotmatched_"+(*it3)]->Fill(BDTweight[bdt_name],wgt);
		}
		h_lep_pT_nL[triglist.at(i)+"_lepnotmatched_"+(*it3)]->Fill(lepPt[idx],wgt);
		h_lep_eta_nL[triglist.at(i)+"_lepnotmatched_"+(*it3)]->Fill(lepEta[idx],wgt);
		h_lep_pT_eta_nL[triglist.at(i)+"_lepnotmatched_"+(*it3)]->Fill(lepPt[idx],lepEta[idx],wgt);
	      }
	      if(lepIsTight){
		for (TString bdt_name : bdt_vec) {		 
		  h_lep_BDT_nT[bdt_name+"_"+triglist.at(i)+"_eventTrig_"+(*it3)]->Fill(BDTweight[bdt_name],wgt);
		}
		h_lep_trig_nT["eventTrig_"+(*it3)]->AddBinContent(i+1,wgt);
		h_lep_pT_nT[triglist.at(i)+"_eventTrig_"+(*it3)]->Fill(lepPt[idx],wgt);
		h_lep_eta_nT[triglist.at(i)+"_eventTrig_"+(*it3)]->Fill(lepEta[idx],wgt);
		h_lep_pT_eta_nT[triglist.at(i)+"_eventTrig_"+(*it3)]->Fill(lepPt[idx],lepEta[idx],wgt);
		if(!lepistriggermatched){
		  for (TString bdt_name : bdt_vec) {		 
		    h_lep_BDT_nT[bdt_name+"_"+triglist.at(i)+"_lepnotmatched_"+(*it3)]->Fill(BDTweight[bdt_name],wgt);
		  }
		  h_lep_pT_nT[triglist.at(i)+"_lepnotmatched_"+(*it3)]->Fill(lepPt[idx],wgt);
		  h_lep_eta_nT[triglist.at(i)+"_lepnotmatched_"+(*it3)]->Fill(lepEta[idx],wgt);
		  h_lep_pT_eta_nT[triglist.at(i)+"_lepnotmatched_"+(*it3)]->Fill(lepPt[idx],lepEta[idx],wgt);
		}
	      }
	    }
	  }
	}
	  */
	  
	//if((*it3).Contains("trueLF"))cout<<"key 1 "<<(*it3).Data()<<endl;
	
	fillFakeRateHist((*it3), idx, lepIsTight);
	for (TString bdt_name : bdt_vec) {
	  h_lep_BDT_nL[bdt_name+"_"+(*it3)]->Fill(BDTweight[bdt_name],wgt_loose);
	}
	h_lep_pT_nL[(*it3)]->Fill(lepPt[idx],wgt_loose);
	h_lep_eta_nL[(*it3)]->Fill(lepEta[idx],wgt_loose);
	h_lep_pT_eta_nL[(*it3)]->Fill(lepPt[idx],lepEta[idx],wgt_loose);
	h_lep_TCF_nL[(*it3)]->AddBinContent(tcf+1,wgt_loose);
	h_lep_mT2_TCF_nL[(*it3)]->Fill(themt2,tcf,wgt_loose);
	if(doThreeLep)h_lep_mTW_TCF_nL[(*it3)]->Fill(mtw,tcf,wgt_loose);
	h_lep_njet_TCF_nL[(*it3)]->Fill(n_clj,tcf,wgt_loose);
	h_lep_metsig_TCF_nL[(*it3)]->Fill(*met_Sign,tcf,wgt_loose);
	h_lep_pT_TCF_nL[(*it3)]->Fill(lepPt[idx],tcf,wgt_loose);
	h_lep_pT_origin_nL[(*it3)]->Fill(lepPt[idx],lepOrigin[idx],wgt_loose);
	h_lep_pT_type_nL[(*it3)]->Fill(lepPt[idx],lepType[idx],wgt_loose);
	h_lep_pT_firstEgMotherO_nL[(*it3)]->Fill(lepPt[idx],lepEgMotherOrigin[idx],wgt_loose);
	h_lep_pT_firstEgMotherT_nL[(*it3)]->Fill(lepPt[idx],lepEgMotherType[idx],wgt_loose);
	h_lep_pT_firstEgMotherPdgId_nL[(*it3)]->Fill(lepPt[idx],lepEgMotherOrigin[idx],wgt_loose);
	h_lep_type_origin_nL[(*it3)]->Fill(lepType[idx],lepOrigin[idx],wgt_loose);
	h_lep_mu_pT_nL[key]->Fill(*mu,lepPt[idx],wgt_loose);
	if(fullScan){
	  if(lepPt[idx] < 30)     h_lep_type_origin_nL[(*it3)+"_"+"pT_min_30"]->Fill(lepType[idx],lepOrigin[idx],wgt_loose);
	  else if(lepPt[idx] < 40)h_lep_type_origin_nL[(*it3)+"_"+"pT_30_40" ]->Fill(lepType[idx],lepOrigin[idx],wgt_loose);
	  else if(lepPt[idx] < 50)h_lep_type_origin_nL[(*it3)+"_"+"pT_40_50" ]->Fill(lepType[idx],lepOrigin[idx],wgt_loose);
	  else if(lepPt[idx] < 60)h_lep_type_origin_nL[(*it3)+"_"+"pT_50_60" ]->Fill(lepType[idx],lepOrigin[idx],wgt_loose);
	  else if(lepPt[idx] < 80)h_lep_type_origin_nL[(*it3)+"_"+"pT_60_80" ]->Fill(lepType[idx],lepOrigin[idx],wgt_loose);
	  else                    h_lep_type_origin_nL[(*it3)+"_"+"pT_80_max"]->Fill(lepType[idx],lepOrigin[idx],wgt_loose);
	}
	if(lepIsTight){
	  for (TString bdt_name : bdt_vec) {
	    h_lep_BDT_nT[bdt_name+"_"+(*it3)]->Fill(BDTweight[bdt_name],wgt);
	  }
	  h_lep_pT_nT[(*it3)]->Fill(lepPt[idx],wgt);
	  h_lep_eta_nT[(*it3)]->Fill(lepEta[idx],wgt);
	  h_lep_pT_eta_nT[(*it3)]->Fill(lepPt[idx],lepEta[idx],wgt);
	  h_lep_TCF_nT[(*it3)]->AddBinContent(tcf+1,wgt);
	  h_lep_mT2_TCF_nT[(*it3)]->Fill(themt2,tcf,wgt);
	  if(doThreeLep)h_lep_mTW_TCF_nT[(*it3)]->Fill(mtw,tcf,wgt);
	  h_lep_njet_TCF_nT[(*it3)]->Fill(n_clj,tcf,wgt);
	  h_lep_metsig_TCF_nT[(*it3)]->Fill(*met_Sign,tcf,wgt);
	  h_lep_pT_TCF_nT[(*it3)]->Fill(lepPt[idx],tcf,wgt);
	  h_lep_pT_origin_nT[(*it3)]->Fill(lepPt[idx],lepOrigin[idx],wgt);
	  h_lep_pT_type_nT[(*it3)]->Fill(lepPt[idx],lepType[idx],wgt);
	  h_lep_pT_firstEgMotherO_nT[(*it3)]->Fill(lepPt[idx],lepEgMotherOrigin[idx],wgt);
	  h_lep_pT_firstEgMotherT_nT[(*it3)]->Fill(lepPt[idx],lepEgMotherType[idx],wgt);
	  h_lep_pT_firstEgMotherPdgId_nT[(*it3)]->Fill(lepPt[idx],lepEgMotherOrigin[idx],wgt);
	  h_lep_type_origin_nT[(*it3)]->Fill(lepType[idx],lepOrigin[idx],wgt);
	  h_lep_mu_pT_nT[key]->Fill(*mu,lepPt[idx],wgt);
	  if(fullScan){
	    if(lepPt[idx] < 30)     h_lep_type_origin_nT[(*it3)+"_"+"pT_min_30"]->Fill(lepType[idx],lepOrigin[idx],wgt);
	    else if(lepPt[idx] < 40)h_lep_type_origin_nT[(*it3)+"_"+"pT_30_40" ]->Fill(lepType[idx],lepOrigin[idx],wgt);
	    else if(lepPt[idx] < 50)h_lep_type_origin_nT[(*it3)+"_"+"pT_40_50" ]->Fill(lepType[idx],lepOrigin[idx],wgt);
	    else if(lepPt[idx] < 60)h_lep_type_origin_nT[(*it3)+"_"+"pT_50_60" ]->Fill(lepType[idx],lepOrigin[idx],wgt);
	    else if(lepPt[idx] < 80)h_lep_type_origin_nT[(*it3)+"_"+"pT_60_80" ]->Fill(lepType[idx],lepOrigin[idx],wgt);
	    else                    h_lep_type_origin_nT[(*it3)+"_"+"pT_80_max"]->Fill(lepType[idx],lepOrigin[idx],wgt);
	  }
	}
      }
      }// end regions
    }// end fillregions
  }
  // re-set wgt before returning
  wgt = old_wgt;
  wgt_loose = old_wgt_loose;
}

void MySusySkimAnalysis::isREALRegion(int idx0, int idx1, int idx2, TString LT_postfix, bool isNOOR)
{
  std::vector< TString > vec;
  TLorentzVector lep1;
  TLorentzVector lep2;
  lep1.SetPtEtaPhiM(lepPt[idx1],lepEta[idx1],lepPhi[idx1],lepM[idx1]);
  lep2.SetPtEtaPhiM(lepPt[idx2],lepEta[idx2],lepPhi[idx2],lepM[idx2]);
  Double_t mll = (lep1+lep2).M();



  if(doThreeLep){
    if(fabs((mll-zmass)) < 15.0){
      if(isNOOR)vec.push_back("2L12");
      else vec.push_back("2L02");
    }
    if(!isNOOR){
      //if(fabs((mll-zmass)) < 15.0 && n_clj < 4 && n_bjet == 0)vec.push_back("2L03");
      //if(lepPt[idx1] > 25 && lepPt[idx2] > 20)vec.push_back("2L01");
      if(!doRJR /**&& lepPt[idx0] > 25 && lepPt[idx1] > 25 && lepPt[idx2] > 20*/ && n_bjet == 0)vec.push_back("2L04");
      else if(doRJR && lepPt[idx0] > 25 && lepPt[idx1] > 25 && lepPt[idx2] > 20 && n_bjet == 0){
	vec.push_back("2L04");
	if(n_clj == 0)vec.push_back("2L05");
	if(n_clj == 1)vec.push_back("2L06");
	if(fabs((mll-zmass)) > 30)vec.push_back("2L07");
	//if((*met_Sign) > 10)vec.push_back("2L08");
	if((*met_Et) > 80)vec.push_back("2L09");
      }
    }
  }else{
    vec.push_back("2L01");
    if(fabs((mll-zmass)) < 10.0){
      if(isNOOR){
	vec.push_back("2L13");
      }else{
	if((*nJet30 >= 2 && LT_postfix.EqualTo("2L2J")) || n_clj < 2)vec.push_back("2L03");	
	vec.push_back("2L02");	
	if(lepFlavor[idx1] == 2 && lepFlavor[idx2] == 2)cutsregions.push_back("M_REAL2L02");
	if(lepFlavor[idx1] == 1 && lepFlavor[idx2] == 1)cutsregions.push_back("E_REAL2L02");
      }
    }
    if(!isNOOR){
      //if(fabs((mll-zmass)) < 10.0 && lepCharge[idx1]*lepCharge[idx2]<0 && n_clj < 2)vec.push_back("2L03");
      //if(mll > 25.0 && themt2 > 50 && lepPt[idx1] > 25 && lepPt[idx2] > 20)vec.push_back("2L01");
      if(!doRJR && LT_postfix.EqualTo("2L2J") /**&& mll > 25.0*/ && lepPt[idx1] > 25 && lepPt[idx2] > 25 && *nJet30 >= 1){
	vec.push_back("2L04");
	
	// high mass
	if(fabs((mll-zmass)) < 20.0 && n_bjet <= 1)vec.push_back("2L10");
	if(fabs((mll-zmass)) < 20.0 && n_bjet >= 2)vec.push_back("2L11");
	
	// low/intermediate mass
	if(fabs((mll-zmass)) < 10.0 && n_bjet >= 1 && *nJet30 >= 2)vec.push_back("2L09");

	// low mass
	if(fabs((mll-zmass)) < 10.0 && n_bjet == 0 && *nJet30 == 2)vec.push_back("2L08");

	// off-shell
	if((mll > 12 && mll < 71) && n_bjet == 0 && *nJet30 >= 2)vec.push_back("2L07");


	//strong
	if(mll > 12 && fabs((mll-zmass)) < 10.0 && *nJet30 >= 2)vec.push_back("2L06");
	if(mll > 12 && *nJet30 >= 2)vec.push_back("2L05");
      }else if(!doRJR && LT_postfix.EqualTo("ZMET") && mll > 70 && ((lepFlavor[idx1] == 2 && lepPt[idx1] > 27 && lepPt[idx2] > 20) || (lepFlavor[idx1] == 1 && lepPt[idx1] > 25))){
	  vec.push_back("2L04");
	  if(n_bjet == 0)vec.push_back("2L05");
	  if(mll > 180)vec.push_back("2L06");
	  if(metrel_Et > 50)vec.push_back("2L07");
	  if(met_tlv.DeltaPhi((lep1+lep2)) > 2.5)vec.push_back("2L08");
	  if(n_bjet == 0 && mll > 180 && metrel_Et > 50 && met_tlv.DeltaPhi((lep1+lep2)) > 2.5)vec.push_back("2L09");
	  if(metrel_Et > 50 && metrel_Et <= 100)vec.push_back("2L10");
	  if(metrel_Et > 100 && metrel_Et <= 150)vec.push_back("2L11");
	  if(metrel_Et > 150)vec.push_back("2L12");
      }else if(!doRJR && !LT_postfix.EqualTo("2L2J") && !LT_postfix.EqualTo("ZMET") && mll > 11 && *met_Sign > 3 && lepPt[idx1] > 27 && lepPt[idx2] > 9 && n_clj <  2){// && mll > 25.0 && lepPt[idx1] > 25 && lepPt[idx2] > 25){
	h_ptllboost->Fill((met_tlv+lep1+lep2).Pt(),wgt);
	if(((lepFlavor[idx1] == lepFlavor[idx2]) && fabs(mll-zmass) > 15.0) || (lepFlavor[idx1] != lepFlavor[idx2])){
	  vec.push_back("2L04");
	  if(n_bjet >=  1){
	    vec.push_back("2L13");
	    if(n_clj == 0)vec.push_back("2L12");
	    //if(n_clj == 1)vec.push_back("2L12");
	  }else if(n_bjet == 0){
	    if(n_clj == 0)vec.push_back("2L11");
	    //if(n_clj == 1)vec.push_back("2L12");
	    if(lep1.DeltaPhi(lep2) > 2.2)vec.push_back("2L05");
	    if(TMath::ATan(fabs(lep1.Eta()-lep2.Eta())/2.)<0.2)vec.push_back("2L06");
	    if(lep1.Pt() > 140 && lep2.Pt() > 20)vec.push_back("2L07");
	    if(mll > 60)vec.push_back("2L08"); 
	    if((met_tlv+lep1+lep2).Pt() < 5)vec.push_back("2L09");
	    if((met_tlv.DeltaPhi(lep1)>2.2))vec.push_back("2L10");
	  }
	}
      }else if(doRJR /**&& lepPt[idx1] > 25 && lepPt[idx2] > 25*/ && *nJet30 >= 2 && n_bjet == 0){
	vec.push_back("2L04");
	if(n_clj == 0)vec.push_back("2L05");
	if(n_clj == 1)vec.push_back("2L06");
	if(fabs((mll-zmass)) > 30)vec.push_back("2L07");
	//if((*met_Sign) > 10)vec.push_back("2L08");
	if((*met_Et) > 80)vec.push_back("2L09");
      }
    }
  }
  
  if(vec.size()>0){
    fillLandT(idx1,vec,LT_postfix,false);
    fillLandT(idx2,vec,LT_postfix,false);
    fillLandT(idx1,vec,LT_postfix,false,*trigMatch_2LTrigOR,"2LTrig");
    fillLandT(idx2,vec,LT_postfix,false,*trigMatch_2LTrigOR,"2LTrig");
/**
    fillLandT(idx1,vec,false,*trigMatch_1L2LTrigOR,"1L2LTrig");
    fillLandT(idx2,vec,false,*trigMatch_1L2LTrigOR,"1L2LTrig");
    
    fillLandT(idx1,vec,false,*trigMatch_1LTrigOR,"1LTrig");
    fillLandT(idx2,vec,false,*trigMatch_1LTrigOR,"1LTrig");
*/


    // fillLandT(idx1,vec,false,*trigMatch_1LTrigOR,"1LTrig");
    // fillLandT(idx2,vec,false,*trigMatch_1LTrigOR,"1LTrig");
  }
}

void MySusySkimAnalysis::isLFRegion(std::vector<int> LF_probe_el, std::vector<int> LF_tag_el, TString LT_postfix)
{
  std::vector< TString > vec;
  TLorentzVector lep_probe;
  TLorentzVector lep_tag;
  int lf_el_probe = LF_probe_el.at(0);
  int lf_el_tag = LF_tag_el.at(0);

  lep_probe.SetPtEtaPhiM(lepPt[lf_el_probe],lepEta[lf_el_probe],lepPhi[lf_el_probe],lepM[lf_el_probe]);

  lep_tag.SetPtEtaPhiM(lepPt[lf_el_tag],lepEta[lf_el_tag],lepPhi[lf_el_tag],lepM[lf_el_tag]);

  TLorentzVector tot_2lep = (TLorentzVector)lep_probe+lep_probe;
  Float_t invMass_2L_LF = tot_2lep.M();
  Float_t mt = (lep_tag.Mt() + met_tlv.Mt())*(lep_tag.Mt() + met_tlv.Mt()) - (lep_tag+met_tlv).Perp2();
  mt = (mt >= 0.) ? sqrt(mt) : sqrt(-mt);

  if(/***nBJet20_MV2c10_FixedCutBEff_77 == 0 && */mt > 40 && *met_Et > 25 && fabs(invMass_2L_LF-zmass)>10)vec.push_back("2L49");
  if(vec.size()>0)fillLandT(lf_el_probe,vec,LT_postfix,true);
}

void MySusySkimAnalysis::isHFRegion(std::vector<int> HF_probe_lep, std::vector<int> HF_tag_mu, int n_bjet, TString LT_postfix, bool isTightTag)
{
  std::vector< TString > vec;
  std::vector< int >::iterator int_it;
  int idx_hf_probe_lep = -1;
  int hf_mu_1;
  int hf_lep_2;
  TLorentzVector lep1;
  TLorentzVector lep2;
  for(int_it = HF_probe_lep.begin(); int_it != HF_probe_lep.end(); int_it++){
    idx_hf_probe_lep += 1;
    if(HF_tag_mu.at(0) == (*int_it))continue;
    hf_mu_1 = HF_tag_mu.at(0);

    if(!lepPassOR->at(*int_it))continue;

    lep1.SetPtEtaPhiM(lepPt[hf_mu_1],lepEta[hf_mu_1],lepPhi[hf_mu_1],lepM[hf_mu_1]);
    lep2.SetPtEtaPhiM(lepPt[(*int_it)],lepEta[(*int_it)],lepPhi[(*int_it)],lepM[(*int_it)]);
    TLorentzVector tot_2lep = (TLorentzVector)lep1+lep2;
    Float_t invMass_2L_HF = tot_2lep.M();
    Float_t mt = (lep2.Mt() + met_tlv.Mt())*(lep2.Mt() + met_tlv.Mt()) - (lep2+met_tlv).Perp2();
    mt = (mt >= 0.) ? sqrt(mt) : sqrt(-mt);

    if(n_bjet == 1){
      if(mt < 50 && (*met_Et) < 50){
	if(fabs(lepFlavor[(*int_it)]) == 2 && fabs(invMass_2L_HF-zmass)>10){
	  if(!isTightTag){
	    //vec.push_back("2L20");
	    cutsregions.push_back("M_FAKE2L20");
	    if(*nJet30 >= 2){
	      vec.push_back("2L21");
	      cutsregions.push_back("M_FAKE2L21");
	    }
	  }// else{
	  //   vec.push_back("2L30");
	  // }
	}else if(fabs(lepFlavor[(*int_it)]) == 1){
	  if(!isTightTag){
	    //vec.push_back("2L20");
	    cutsregions.push_back("E_FAKE2L20");
	    if(*nJet30 >= 2){
	      vec.push_back("2L21");
	      cutsregions.push_back("E_FAKE2L21");
	    }
	  }// else{
	  //   vec.push_back("2L30");
	  // }
	}
	
      }
    }
    if(vec.size()>0){
      fillLandT((*int_it),vec,LT_postfix,true);
      fillLandT((*int_it),vec,LT_postfix,true,*trigMatch_2LTrigOR,"2LTrig");

      /**
      fillLandT((*int_it),vec,true,*trigMatch_1L2LTrigOR,"1L2LTrig");
      fillLandT((*int_it),vec,true,*trigMatch_1LTrigOR,"1LTrig");
     
      */
      //fillLandT((*int_it),vec,true,*trigMatch_1LTrigOR,"1LTrig");
    }
    vec.clear();
  }
}

void MySusySkimAnalysis::isConvRegion(int co_mu_1, int co_mu_2, int co_el_1, int n_bjet, TString LT_postfix)
{
  if(!makeHistograms)return;
  std::vector<TString> vec;
  if(!lepPassOR->at(co_el_1))return;
  
  if(lepCharge[co_mu_1]*lepCharge[co_mu_2] < 0/** && lepPassOR->at(co_el_1)*/){
    TLorentzVector mu1;
    TLorentzVector mu2;
    TLorentzVector el1;
    mu1.SetPtEtaPhiM(lepPt[co_mu_1],lepEta[co_mu_1],lepPhi[co_mu_1],lepM[co_mu_1]);
    mu2.SetPtEtaPhiM(lepPt[co_mu_2],lepEta[co_mu_2],lepPhi[co_mu_2],lepM[co_mu_2]);
    el1.SetPtEtaPhiM(lepPt[co_el_1],lepEta[co_el_1],lepPhi[co_el_1],lepM[co_el_1]);
    TLorentzVector tot_3L = (TLorentzVector)mu1+mu2+el1;
    TLorentzVector tot_2mu = (TLorentzVector)mu1+mu2;
    Double_t invMass_3L = tot_3L.M();
    Double_t invMass_2mu = tot_2mu.M();


    TLorentzVector truthel1;
    truthel1.SetPtEtaPhiM(lepTruthPt[co_el_1],lepTruthEta[co_el_1],lepTruthPhi[co_el_1],lepTruthM[co_el_1]);

    
    
    // printf("fabs(invMass_3L-zmass) = %.2f\n",fabs(invMass_3L-zmass));
    // printf("met = %.2f\n",*met_Et);
    if(fabs(invMass_3L-zmass)<10 && fabs(lepEta[co_el_1]) < 1.4){
      vec.push_back("2L26");
    }

    if(n_bjet == 0 && fabs(invMass_3L-zmass)<10 && *met_Et<50.) {
      // if(truthel1.Px() != 0){
      // 	/**if(lepIFFClass[co_el_1] == 6)*/
      // 	myfile <<(*DatasetNumber)<<" "<<(*EventNumber)<<" "<<lepIFFClass[co_el_1]<<" "<<lepEgMotherOrigin[co_el_1]<<" "
      // 	       <<lepEgMotherPdgId[co_el_1]<<" "<<lepType[co_el_1]<<" "<<lepOrigin[co_el_1]
      // 	       <<" "<<lepEgMotherType[co_el_1]<<" "<<truthel1.Px()<<" "<<truthel1.Py()<<" "
      // 	       <<truthel1.Pz()<<" "<<lepTruthPt[co_el_1]<<" "<<lepTruthEta[co_el_1]<<" "
      // 	       <<lepPassOR->at(co_el_1)<<" "<<lepCharge[co_el_1]<<" "<<lepTruthCharge[co_el_1]<<"\n";
      // }else{
      // 	/**if(lepIFFClass[co_el_1] == 6)*/
      // 	myfile <<(*DatasetNumber)<<" "<<(*EventNumber)<<" "<<lepIFFClass[co_el_1]<<" "<<lepEgMotherOrigin[co_el_1]<<" "
      // 	       <<lepEgMotherPdgId[co_el_1]<<" "<<lepType[co_el_1]<<" "<<lepOrigin[co_el_1]
      // 	       <<" "<<lepEgMotherType[co_el_1]<<" "<<truthel1.Px()<<" "<<truthel1.Py()<<" "
      // 	       <<truthel1.Pz()<<" "<<lepTruthPt[co_el_1]<<" "<<lepTruthEta[co_el_1]<<" "
      // 	       <<lepPassOR->at(co_el_1)<<" "<<lepCharge[co_el_1]<<" "<<lepTruthCharge[co_el_1]<<"\n";
      // }
      vec.push_back("2L23");
      h_lep_lepECIDS->Fill(lepECIDS[co_el_1]);
      if(!lepECIDS[co_el_1]){
	vec.push_back("2L24");
	if(*nJet30 >= 2){
	  vec.push_back("2L27");
	  cutsregions.push_back("E_FAKE2L27");
	}
      }else{
	vec.push_back("2L25");
	cutsregions.push_back("E_FAKE2L25");
      }
    }
  }

  if(vec.size()>0){
    fillLandT(co_el_1,vec,LT_postfix);
    fillLandT(co_el_1,vec,LT_postfix,true,*trigMatch_2LTrigOR,"2LTrig");
    /**
    fillLandT(co_el_1,vec,true,*trigMatch_1L2LTrigOR,"1L2LTrig");
    fillLandT(co_el_1,vec,true,*trigMatch_1LTrigOR,"1LTrig");
    */
    
    //fillLandT(co_el_1,vec,true,*trigMatch_1LTrigOR,"1LTrig");
  }
  vec.clear();
}

std::pair <TString,TString> MySusySkimAnalysis::getFillString(int idx1, int idx2, int idx0)
{
  TString chstr = "";
  TString sgnstr = "";

  if(doThreeLep){
    if(lepFlavor[idx1] == 1 && lepFlavor[idx2] == 1 && lepFlavor[idx0] == 1){
      chstr = "EEE";
    }else if(lepFlavor[idx1] == 2 && lepFlavor[idx2] == 2 && lepFlavor[idx0] == 2){
      chstr = "MMM";
    }else if((lepFlavor[idx1] == 2 && lepFlavor[idx2] == 2 && lepFlavor[idx0] == 1) ||
	     (lepFlavor[idx1] == 1 && lepFlavor[idx2] == 2 && lepFlavor[idx0] == 2) ||
	     (lepFlavor[idx1] == 2 && lepFlavor[idx2] == 1 && lepFlavor[idx0] == 2)){
      chstr = "MME";
    }else if((lepFlavor[idx1] == 1 && lepFlavor[idx2] == 1 && lepFlavor[idx0] == 2) ||
	     (lepFlavor[idx1] == 2 && lepFlavor[idx2] == 1 && lepFlavor[idx0] == 1) ||
	     (lepFlavor[idx1] == 1 && lepFlavor[idx2] == 2 && lepFlavor[idx0] == 1)){
      chstr = "EEM";
    }else{
      cout<<"ERROR \t Strange error flav1 = "<<lepFlavor[idx1]<<", flav2 = "<<lepFlavor[idx2]<<", flav3 = "<<lepFlavor[idx0]<<endl;
    }
    sgnstr = "";
  }else{
    if(lepFlavor[idx1] == 1 && lepFlavor[idx2] == 1){
      chstr = "EE";
    }else if(lepFlavor[idx1] == 2 && lepFlavor[idx2] == 2){
      chstr = "MM";
    }else if(lepFlavor[idx1] == 1 && lepFlavor[idx2] == 2){
      chstr = "EM";
    }else if(lepFlavor[idx1] == 2 && lepFlavor[idx2] == 1){
      chstr = "EM";
    }else{
      cout<<"ERROR \t Strange error flav1 = "<<lepFlavor[idx1]<<", flav2 = "<<lepFlavor[idx2]<<endl;
    }


    if(lepCharge[idx1]*lepCharge[idx2] > 0)sgnstr = "SS";
    else if(lepCharge[idx1]*lepCharge[idx2] < 0)sgnstr = "OS";

  }

  std::pair <TString,TString> ret (chstr,sgnstr);

  return ret;
}

std::vector<TString> MySusySkimAnalysis::getFillVector(int idx1, int idx2, int idx0)
{

  TString chstr = "";
  TString sgnstr = "";
  std::vector<TString> regions;
  regions.push_back("ALL");


  std::pair <TString,TString> ret = getFillString(idx1, idx2, idx0);

  chstr = ret.first;
  sgnstr = ret.second;

  if(doThreeLep){
    regions.push_back(chstr);
  }
  //if(!doThreeLep && (chstr.Contains("EE") || chstr.Contains("MM")))regions.push_back("SF"+sgnstr);
  regions.push_back(chstr+sgnstr);
  return regions;

}

void MySusySkimAnalysis::classifyEvent(int idx1, int idx2, int& isTT, int& isTL, int& isLT, int& isLL )
{

  isTT = 0;
  isTL = 0;
  isLT = 0;
  isLL = 0;
  if(isT(idx1) && isT(idx2))isTT = 1;
  else if(isT(idx1) && !isT(idx2))isTL = 1;
  else if(!isT(idx1) && isT(idx2))isLT = 1;
  else if(!isT(idx1) && !isT(idx2))isLL = 1;

}

std::vector<double> MySusySkimAnalysis::getFakeWeight1D(int idx1, int idx2, int idx0, TString year, TString LT_postfix, int vb){

std::vector<double> vec;
return vec;

}
//   int nTT = 0;
//   int nTl = 0;
//   int nlT = 0;
//   int nll = 0;

//   double r_wgt  = 0.0;
//   double rf_wgt = 0.0;
//   double fr_wgt = 0.0;
//   double f_wgt  = 0.0;
 
//   TString eventStr;

//   classifyEvent(idx1, idx2, nTT, nTl, nlT, nll);

//   if(vb)printf("nTT = %i, nTl = %i, nlT = %i, nll = %i\n",nTT, nTl, nlT, nll);
  
//   std::pair <TString,TString> ret = getFillString(idx1, idx2, idx0);
//   eventStr = ret.first + ret.second;

//   if(vb)printf("eventStr = %s\n",eventStr.Data());

//   gfw->setVariables(year, LT_postfix, eventStr);
//   //
//   gfw->setFrac(*met_Sign, *mll, *nJet30, n_bjet, zmass, *met_Sign);
  
//   if(vb)printf("metsig = %.2f, mll = %.2f, njet30 = %i, nbjet = %i, zmass = %.2f\n",*met_Sign,*mll,*nJet30,n_bjet,zmass);

//   // a few hacks
//   double pT1_fake    = lepPt[idx1] < 200 ? lepPt[idx1] : 199.0;
//   double pT1_fake_co = lepPt[idx1] < 100 ? lepPt[idx1] : 99.0;
//   double pT1_real    = pT1_fake;
//   double pT2_fake    = lepPt[idx2] < 200 ? lepPt[idx2] : 199.0;
//   double pT2_fake_co = lepPt[idx2] < 100 ? lepPt[idx2] : 99.0;
//   double pT2_real    = pT2_fake;

//   if(vb){
//     printf("pt1fake = %.2f, pt2fake = %.2f\n",pT1_fake,pT2_fake);
//     printf("pt1real = %.2f, pt2real = %.2f\n",pT1_real,pT2_real);
//   }
  
//   r_wgt  = mtx->N4_RR_TT(gfw->getReal(pT1_real, lepEta[idx1], lepFlavor[idx1]),
// 			 gfw->getFake(pT1_fake, lepEta[idx1], lepFlavor[idx1]),
// 			 gfw->getReal(pT2_real, lepEta[idx2], lepFlavor[idx2]),
// 			 gfw->getFake(pT2_fake, lepEta[idx2], lepFlavor[idx2]),
// 			 nTT,nTl,nlT,nll);
//   if(vb)printf("--> r_wgt = %.6f\n",r_wgt);
//   rf_wgt = mtx->N4_RF_TT(gfw->getReal(pT1_real, lepEta[idx1], lepFlavor[idx1]),
// 			 gfw->getFake(pT1_fake, lepEta[idx1], lepFlavor[idx1]),
// 			 gfw->getReal(pT2_real, lepEta[idx2], lepFlavor[idx2]),
// 			 gfw->getFake(pT2_fake, lepEta[idx2], lepFlavor[idx2]),
// 			 nTT,nTl,nlT,nll);
//   if(vb)printf("--> rf_wgt = %.6f\n",rf_wgt);
//   fr_wgt = mtx->N4_FR_TT(gfw->getReal(pT1_real, lepEta[idx1], lepFlavor[idx1]),
// 			 gfw->getFake(pT1_fake, lepEta[idx1], lepFlavor[idx1]),
// 			 gfw->getReal(pT2_real, lepEta[idx2], lepFlavor[idx2]),
// 			 gfw->getFake(pT2_fake, lepEta[idx2], lepFlavor[idx2]),
// 			 nTT,nTl,nlT,nll);
//   if(vb)printf("--> fr_wgt = %.6f\n",fr_wgt);
//   f_wgt  = mtx->N4_FF_TT(gfw->getReal(pT1_real, lepEta[idx1], lepFlavor[idx1]),
// 			 gfw->getFake(pT1_fake, lepEta[idx1], lepFlavor[idx1]),
// 			 gfw->getReal(pT2_real, lepEta[idx2], lepFlavor[idx2]),
// 			 gfw->getFake(pT2_fake, lepEta[idx2], lepFlavor[idx2]),
// 			 nTT,nTl,nlT,nll);
//   if(vb)printf("--> f_wgt = %.6f\n",f_wgt);

//   double f_wgt_syst  = mtx->N4_FF_TT(gfw->getReal(pT1_real, lepEta[idx1], lepFlavor[idx1],"",GetFakeWeight::XSECUP),
// 				     gfw->getFake(pT1_fake, lepEta[idx1], lepFlavor[idx1],"",GetFakeWeight::XSECUP),
// 				     gfw->getReal(pT2_real, lepEta[idx2], lepFlavor[idx2],"",GetFakeWeight::XSECUP),
// 				     gfw->getFake(pT2_fake, lepEta[idx2], lepFlavor[idx2],"",GetFakeWeight::XSECUP),
// 				     nTT,nTl,nlT,nll);
//   if(vb)printf("--> f_wgt_syst = %.6f\n",f_wgt_syst);


//   std::vector<double> retvec;
//   double fake = rf_wgt+fr_wgt+f_wgt;
//   retvec.push_back(fake);
//   retvec.push_back(fake);// > 0 ? (1+frac_total_up) : (1-frac_total_up));
//   retvec.push_back(fake);// > 0 ? (1-frac_total_dw) : (1+frac_total_dw));

//   // setting some global variables (for ntuple)
//   glob_r1 = gfw->getReal(lepPt[idx1], lepEta[idx1], lepFlavor[idx1]);
//   glob_r2 = gfw->getReal(lepPt[idx2], lepEta[idx2], lepFlavor[idx2]);
//   glob_f1 = gfw->getFake(lepPt[idx1], lepEta[idx1], lepFlavor[idx1]);
//   glob_f2 = gfw->getFake(lepPt[idx2], lepEta[idx2], lepFlavor[idx2]);

//   glob_r_wgt = r_wgt;
//   glob_rf_wgt = rf_wgt;
//   glob_fr_wgt = fr_wgt;
//   glob_f_wgt = f_wgt;
  

//   return retvec;

// }



  

// std::vector<double> MySusySkimAnalysis::getFakeWeight1D(int idx1, int idx2, int idx0, TString year, TString LT_postfix, int vb)
// {

//   int nTT = 0;
//   int nTl = 0;
//   int nlT = 0;
//   int nll = 0;

//   TLorentzVector lep1;
//   TLorentzVector lep2;
//   lep1.SetPtEtaPhiM(lepPt[idx1],lepEta[idx1],lepPhi[idx1],lepM[idx1]);
//   lep2.SetPtEtaPhiM(lepPt[idx2],lepEta[idx2],lepPhi[idx2],lepM[idx2]);
//   Double_t mll = (lep1+lep2).M();

//   TString year_const = year;

//   double r1 = 0.0;
//   double f1 = 0.0;
//   double r2 = 0.0;
//   double f2 = 0.0;

//   TString eventStr = "";
//   TString fulleventStr = "";

//   double r_wgt  = 0.0;
//   double rf_wgt = 0.0;
//   double fr_wgt = 0.0;
//   double f_wgt  = 0.0;
//   double fake_wgt = 0.0;

//   std::vector<double> r_wgt_var;
//   std::vector<double> rf_wgt_var;
//   std::vector<double> fr_wgt_var;
//   std::vector<double> f_wgt_var;

//   std::vector<double> r_wgt_stat;
//   std::vector<double> rf_wgt_stat;
//   std::vector<double> fr_wgt_stat;
//   std::vector<double> f_wgt_stat;
//   std::vector<double> fake_wgt_stat;

//   std::vector<double> r_wgt_syst;
//   std::vector<double> rf_wgt_syst;
//   std::vector<double> fr_wgt_syst;
//   std::vector<double> f_wgt_syst;
//   std::vector<double> fake_wgt_syst;

//   std::vector<double> r_wgt_xsec;
//   std::vector<double> rf_wgt_xsec;
//   std::vector<double> fr_wgt_xsec;
//   std::vector<double> f_wgt_xsec;
//   std::vector<double> fake_wgt_xsec;

//   double r_wgt_mc;
//   double rf_wgt_mc;
//   double fr_wgt_mc;
//   double f_wgt_mc;
//   double fake_wgt_mc;



//   double hf_frac_nom = 0.;
//   double lf_frac_nom = 0.;
//   double co_frac_nom = 0.;

//   double hf_frac_max = 0.;
//   double lf_frac_max = 0.;
//   double co_frac_max = 0.;

//   double hf_frac_min = 0.;
//   double lf_frac_min = 0.;
//   double co_frac_min = 0.;

//   std::vector<double> hf;
//   std::vector<double> lf;
//   std::vector<double> co;

//   std::vector<double> f1_wgt_var;
//   std::vector<double> f2_wgt_var;

//   double f1_xsecup_var = 0.0;
//   double f2_xsecup_var = 0.0;
//   double f1_xsecdw_var = 0.0;
//   double f2_xsecdw_var = 0.0;

//   std::vector<double> r1_stat_var;
//   std::vector<double> r2_stat_var;
//   std::vector<double> f1_stat_var;
//   std::vector<double> f2_stat_var;

//   double r1_mc = 1.0;
//   double r2_mc = 1.0;

//   vector<double>::const_iterator it;

//   classifyEvent(idx1, idx2, nTT, nTl, nlT, nll);
//   std::pair <TString,TString> ret = getFillString(idx1, idx2, idx0);
//   eventStr = ret.first + ret.second;
//   fulleventStr = ret.first + ret.second + "_" + LT_postfix;

//   // a few hacks
//   double pT1_fake    = lepPt[idx1] < 200 ? lepPt[idx1] : 199.0;
//   double pT1_fake_co = lepPt[idx1] < 100 ? lepPt[idx1] : 99.0;
//   double pT1_real    = pT1_fake;
//   double pT2_fake    = lepPt[idx2] < 200 ? lepPt[idx2] : 199.0;
//   double pT2_fake_co = lepPt[idx2] < 100 ? lepPt[idx2] : 99.0;
//   double pT2_real    = pT2_fake;      


//   double eta1_real    = fabs(lepEta[idx1]) < 2.6 ? fabs(lepEta[idx1]) : 2.5;
//   double eta2_real    = fabs(lepEta[idx2]) < 2.6 ? fabs(lepEta[idx2]) : 2.5;

 
//   if(LT_postfix.Contains("HP_")){
//     eta1_real    = fabs(lepEta[idx1]) < 2.5 ? fabs(lepEta[idx1]) : 2.45;
//     eta2_real    = fabs(lepEta[idx2]) < 2.5 ? fabs(lepEta[idx2]) : 2.45;
//   }


//   pT1_fake    = lepPt[idx1] >= 25 ? pT1_fake    : 25.0;	
//   pT1_fake_co = lepPt[idx1] >= 25 ? pT1_fake_co : 25.0;	
//   pT1_real    = pT1_fake;				
//   pT2_fake    = lepPt[idx2] >= 25 ? pT2_fake    : 25.0;	
//   pT2_fake_co = lepPt[idx2] >= 25 ? pT2_fake_co : 25.0;	
//   pT2_real    = pT2_fake;       

                         


//   if(vb > 0){
//     printf("Year is %s, Events is %s, Number is %llu\n",year.Data(),eventStr.Data(),*EventNumber);
//     printf("Lep1 is %i: pT = %.2f (%.2f), |eta| = %.2f\n",lepFlavor[idx1], lepPt[idx1],pT1_fake,eta1_real);
//     printf("Lep2 is %i: pT = %.2f (%.2f), |eta| = %.2f\n",lepFlavor[idx2], lepPt[idx2],pT2_fake,eta2_real);
//     printf("TT = %i, TL = %i, LT = %i, LL = %i, \n",nTT, nTl, nlT, nll);
//   }

//   int mt2bin = 1;
//   float var_dep;
//   if(doThreeLep)var_dep = mtw;
//   else var_dep = themt2;
//   if(var_dep >= 0 && var_dep < 20)mt2bin = 1;
//   else if(var_dep >= 20 && var_dep < 40)mt2bin = 2;
//   else if(var_dep >= 40 && var_dep < 60)mt2bin = 3;
//   else if(var_dep >= 60 && var_dep < 100)mt2bin = 4;
//   else if(var_dep >= 100)mt2bin = 5;
//   if(!(eventStr.EqualTo("MMOS") || eventStr.EqualTo("MMSS") || eventStr.EqualTo("MMM"))){
//     //if(eventStr.Contains("OS") || eventStr.Contains("SS")){
//     for(int i = 4; i<=11; i++){

//       if(!doRJR && LT_postfix.EqualTo("2L2J")){
// 	  if((fabs((mll-zmass)) < 20.0 && n_bjet <= 1) && (i != 4 && i != 10))continue;
// 	  if((fabs((mll-zmass)) < 20.0 && n_bjet >= 2) && (i != 4 && i != 11))continue;
// 	  if((fabs((mll-zmass)) < 10.0 && n_bjet >= 1 && *nJet30 >= 2) && (i != 4 && i != 9))continue;
// 	  if((fabs((mll-zmass)) < 10.0 && n_bjet == 0 && *nJet30 == 2) && (i != 4 && i != 8))continue;
// 	  if(((mll > 12 && mll < 71) && n_bjet == 0 && *nJet30 >= 2) && (i != 4 && i != 7))continue;
// 	  if((mll > 12 && fabs((mll-zmass)) < 10.0 && *nJet30 >= 2) && (i != 4 && i != 6))continue;
// 	  if((mll > 12 && *nJet30 >= 2) && (i != 4 && i != 5))continue;
// 	}
//       //cout<<"Using region "<<i<<endl;
	
// 	//if(h_el_frac_hf.find(Form("%s_%02i_%s",eventStr.Data(),i,LT_postfix.Data())) == h_el_frac_hf.end())continue;
// 	hf.push_back(h_el_frac_hf[Form("%s_%02i_%s",eventStr.Data(),i,LT_postfix.Data())]->GetBinContent(mt2bin)+h_el_frac_cf[Form("%s_%02i_%s",eventStr.Data(),i,LT_postfix.Data())]->GetBinContent(mt2bin));
// 	lf.push_back(h_el_frac_lf[Form("%s_%02i_%s",eventStr.Data(),i,LT_postfix.Data())]->GetBinContent(mt2bin));
// 	co.push_back(h_el_frac_co[Form("%s_%02i_%s",eventStr.Data(),i,LT_postfix.Data())]->GetBinContent(mt2bin));
//     }
//     //}// else{
//     //   hf.push_back(h_el_frac_hf[Form("%s",eventStr.Data())]->GetBinContent(mt2bin));
//     //   lf.push_back(h_el_frac_lf[Form("%s",eventStr.Data())]->GetBinContent(mt2bin));
//     //   co.push_back(h_el_frac_co[Form("%s",eventStr.Data())]->GetBinContent(mt2bin));
//     // }

//     hf_frac_max = *max_element(std::begin(hf), std::end(hf));
//     lf_frac_max = *max_element(std::begin(lf), std::end(lf));
//     co_frac_max = *max_element(std::begin(co), std::end(co));

//     hf_frac_min = *min_element(std::begin(hf), std::end(hf));
//     lf_frac_min = *min_element(std::begin(lf), std::end(lf));
//     co_frac_min = *min_element(std::begin(co), std::end(co));

//     hf_frac_nom = accumulate( hf.begin(), hf.end(), 0.0)/hf.size();
//     lf_frac_nom = accumulate( lf.begin(), lf.end(), 0.0)/lf.size();
//     co_frac_nom = accumulate( co.begin(), co.end(), 0.0)/co.size();

//     if(vb > 0){
//       printf("mt2 = %.2f, mt2bin = %i, wgt(hf) = %.2f, wgt(lf) = %.2f, wgt(co) = %.2f, wgt(sum) = %.2f\n",var_dep,mt2bin,
// 	     hf_frac_nom,
// 	     lf_frac_nom,
// 	     co_frac_nom,
// 	     hf_frac_nom+lf_frac_nom+co_frac_nom);
//       printf("mt2 = %.2f, mt2bin = %i, wgt_up(hf) = %.2f, wgt_up(lf) = %.2f, wgt_up(co) = %.2f, wgt_up(sum) = %.2f\n",var_dep,mt2bin,
// 	     hf_frac_max,
// 	     lf_frac_max,
// 	     co_frac_max,
// 	     hf_frac_max+lf_frac_max+co_frac_max);
//       printf("mt2 = %.2f, mt2bin = %i, wgt_dw(hf) = %.2f, wgt_dw(lf) = %.2f, wgt_dw(co) = %.2f, wgt_dw(sum) = %.2f\n",var_dep,mt2bin,
// 	     hf_frac_min,
// 	     lf_frac_min,
// 	     co_frac_min,
// 	     hf_frac_min+lf_frac_min+co_frac_min);
//     }
//   }
//   if(lepFlavor[idx1] == 1){

//     //r1 = h_el_real[year+"_"+LT_postfix]->GetBinContent(h_el_real[year+"_"+LT_postfix]->FindBin(pT1_real));
//     r1 = h_el_real[year+"_"+LT_postfix]->GetBinContent(h_el_real[year+"_"+LT_postfix]->FindBin(pT1_real,eta1_real));

//     //r1_mc = h_el_real_mc[year+"_"+LT_postfix]->GetBinContent(h_el_real_mc[year+"_"+LT_postfix]->FindBin(pT1_real));
//     r1_mc = h_el_real_mc[year+"_"+LT_postfix]->GetBinContent(h_el_real_mc[year+"_"+LT_postfix]->FindBin(pT1_real,eta1_real));

//     r1_stat_var.push_back(h_el_real[year+"_"+LT_postfix]->GetBinContent(h_el_real[year+"_"+LT_postfix]->FindBin(pT1_real,eta1_real)) + h_el_real[year+"_"+LT_postfix]->GetBinError(h_el_real[year+"_"+LT_postfix]->FindBin(pT1_real,eta1_real)));
//     r1_stat_var.push_back(h_el_real[year+"_"+LT_postfix]->GetBinContent(h_el_real[year+"_"+LT_postfix]->FindBin(pT1_real,eta1_real)) - h_el_real[year+"_"+LT_postfix]->GetBinError(h_el_real[year+"_"+LT_postfix]->FindBin(pT1_real,eta1_real)));

//     f1 = (hf_frac_nom*h_el_fake_hf[year+"_"+LT_postfix]->GetBinContent(h_el_fake_hf[year+"_"+LT_postfix]->FindBin(pT1_fake)) +
//     	  lf_frac_nom*h_el_fake_lf[year+"_"+LT_postfix]->GetBinContent(h_el_fake_lf[year+"_"+LT_postfix]->FindBin(pT1_fake)) +
//     	  co_frac_nom*h_el_fake_co[year+"_"+LT_postfix]->GetBinContent(h_el_fake_co[year+"_"+LT_postfix]->FindBin(pT1_fake_co))
//     	  );

//     f1_stat_var.push_back(hf_frac_nom*(h_el_fake_hf[year+"_"+LT_postfix]->GetBinContent(h_el_fake_hf[year+"_"+LT_postfix]->FindBin(pT1_fake)) + h_el_fake_hf[year+"_"+LT_postfix]->GetBinError(h_el_fake_hf[year+"_"+LT_postfix]->FindBin(pT1_fake))) +
// 			  lf_frac_nom*(h_el_fake_lf[year+"_"+LT_postfix]->GetBinContent(h_el_fake_lf[year+"_"+LT_postfix]->FindBin(pT1_fake)) + h_el_fake_lf[year+"_"+LT_postfix]->GetBinError(h_el_fake_lf[year+"_"+LT_postfix]->FindBin(pT1_fake))) +
// 			  co_frac_nom*(h_el_fake_co[year+"_"+LT_postfix]->GetBinContent(h_el_fake_co[year+"_"+LT_postfix]->FindBin(pT1_fake_co)) + h_el_fake_co[year+"_"+LT_postfix]->GetBinError(h_el_fake_co[year+"_"+LT_postfix]->FindBin(pT1_fake_co)))
// 			  );
//     f1_stat_var.push_back(hf_frac_nom*(h_el_fake_hf[year+"_"+LT_postfix]->GetBinContent(h_el_fake_hf[year+"_"+LT_postfix]->FindBin(pT1_fake)) - h_el_fake_hf[year+"_"+LT_postfix]->GetBinError(h_el_fake_hf[year+"_"+LT_postfix]->FindBin(pT1_fake))) +
// 			  lf_frac_nom*(h_el_fake_lf[year+"_"+LT_postfix]->GetBinContent(h_el_fake_lf[year+"_"+LT_postfix]->FindBin(pT1_fake)) - h_el_fake_lf[year+"_"+LT_postfix]->GetBinError(h_el_fake_lf[year+"_"+LT_postfix]->FindBin(pT1_fake))) +
// 			  co_frac_nom*(h_el_fake_co[year+"_"+LT_postfix]->GetBinContent(h_el_fake_co[year+"_"+LT_postfix]->FindBin(pT1_fake_co)) - h_el_fake_co[year+"_"+LT_postfix]->GetBinError(h_el_fake_co[year+"_"+LT_postfix]->FindBin(pT1_fake_co)))
// 			  );

//     for(unsigned int i = 0; i<hf.size(); i++){
//       f1_wgt_var.push_back((hf.at(i)*h_el_fake_hf[year+"_"+LT_postfix]->GetBinContent(h_el_fake_hf[year+"_"+LT_postfix]->FindBin(pT1_fake)) +
// 			    lf.at(i)*h_el_fake_lf[year+"_"+LT_postfix]->GetBinContent(h_el_fake_lf[year+"_"+LT_postfix]->FindBin(pT1_fake)) +
// 			    co.at(i)*h_el_fake_co[year+"_"+LT_postfix]->GetBinContent(h_el_fake_co[year+"_"+LT_postfix]->FindBin(pT1_fake_co))
// 			    ));
//     }

//     f1_xsecup_var = ((hf_frac_nom*h_el_fake_hf_xsecup[year+"_"+LT_postfix]->GetBinContent(h_el_fake_hf_xsecup[year+"_"+LT_postfix]->FindBin(pT1_fake)) +
// 		      lf_frac_nom*h_el_fake_lf[year+"_"+LT_postfix]->GetBinContent(h_el_fake_lf[year+"_"+LT_postfix]->FindBin(pT1_fake)) +
// 		      co_frac_nom*h_el_fake_co[year+"_"+LT_postfix]->GetBinContent(h_el_fake_co[year+"_"+LT_postfix]->FindBin(pT1_fake_co))
// 		      ));
//     f1_xsecdw_var = ((hf_frac_nom*h_el_fake_hf_xsecdw[year+"_"+LT_postfix]->GetBinContent(h_el_fake_hf_xsecdw[year+"_"+LT_postfix]->FindBin(pT1_fake)) +
// 		      lf_frac_nom*h_el_fake_lf[year+"_"+LT_postfix]->GetBinContent(h_el_fake_lf[year+"_"+LT_postfix]->FindBin(pT1_fake)) +
// 		      co_frac_nom*h_el_fake_co[year+"_"+LT_postfix]->GetBinContent(h_el_fake_co[year+"_"+LT_postfix]->FindBin(pT1_fake_co))
// 		      ));


//   }else if(lepFlavor[idx1] == 2){

//     // if(doRJR){
//     //   if(year.Contains("2017") && pT1_fake >= 60)year = "2018";
//     // }

//     if(LT_postfix.Contains("HP_")){
//       pT1_fake = lepPt[idx1] < 50 ? lepPt[idx1] : 49.0;
//     }else{
//       pT1_fake = lepPt[idx1] < 80 ? lepPt[idx1] : 79.0;
//     }
//     // if(doRJR){
//     //   if(year.Contains("2017"))year = "2018";
//     // }


//     r1 = h_mu_real[year+"_"+LT_postfix]->GetBinContent(h_mu_real[year+"_"+LT_postfix]->FindBin(pT1_real,eta1_real));

//     r1_mc = h_mu_real_mc[year+"_"+LT_postfix]->GetBinContent(h_mu_real_mc[year+"_"+LT_postfix]->FindBin(pT1_real,eta1_real));

//     r1_stat_var.push_back(h_mu_real[year+"_"+LT_postfix]->GetBinContent(h_mu_real[year+"_"+LT_postfix]->FindBin(pT1_real,eta1_real)) + h_mu_real[year+"_"+LT_postfix]->GetBinError(h_mu_real[year+"_"+LT_postfix]->FindBin(pT1_real,eta1_real)));
//     r1_stat_var.push_back(h_mu_real[year+"_"+LT_postfix]->GetBinContent(h_mu_real[year+"_"+LT_postfix]->FindBin(pT1_real,eta1_real)) - h_mu_real[year+"_"+LT_postfix]->GetBinError(h_mu_real[year+"_"+LT_postfix]->FindBin(pT1_real,eta1_real)));



//     f1 = h_mu_fake_hf[year+"_"+LT_postfix]->GetBinContent(h_mu_fake_hf[year+"_"+LT_postfix]->FindBin(pT1_fake));

//     f1_stat_var.push_back(h_mu_fake_hf[year+"_"+LT_postfix]->GetBinContent(h_mu_fake_hf[year+"_"+LT_postfix]->FindBin(pT1_fake)) + h_mu_fake_hf[year+"_"+LT_postfix]->GetBinError(h_mu_fake_hf[year+"_"+LT_postfix]->FindBin(pT1_fake)));
//     f1_stat_var.push_back(h_mu_fake_hf[year+"_"+LT_postfix]->GetBinContent(h_mu_fake_hf[year+"_"+LT_postfix]->FindBin(pT1_fake)) - h_mu_fake_hf[year+"_"+LT_postfix]->GetBinError(h_mu_fake_hf[year+"_"+LT_postfix]->FindBin(pT1_fake)));

//     f1_wgt_var.push_back(f1);
//     //cout<<"b"<<endl;
//     f1_xsecup_var = (h_mu_fake_hf_xsecup[year+"_"+LT_postfix]->GetBinContent(h_mu_fake_hf_xsecup[year+"_"+LT_postfix]->FindBin(pT1_fake)));
//     f1_xsecdw_var = (h_mu_fake_hf_xsecdw[year+"_"+LT_postfix]->GetBinContent(h_mu_fake_hf_xsecdw[year+"_"+LT_postfix]->FindBin(pT1_fake)));
//     //cout<<"c"<<endl;

//     year = year_const;

//   }

//   if(lepFlavor[idx2] == 1){

//     r2 = h_el_real[year+"_"+LT_postfix]->GetBinContent(h_el_real[year+"_"+LT_postfix]->FindBin(pT2_real,eta2_real));

//     r2_mc = h_el_real_mc[year+"_"+LT_postfix]->GetBinContent(h_el_real_mc[year+"_"+LT_postfix]->FindBin(pT2_real,eta2_real));

//     r2_stat_var.push_back(h_el_real[year+"_"+LT_postfix]->GetBinContent(h_el_real[year+"_"+LT_postfix]->FindBin(pT2_real,eta2_real)) + h_el_real[year+"_"+LT_postfix]->GetBinError(h_el_real[year+"_"+LT_postfix]->FindBin(pT2_real,eta2_real)));
//     r2_stat_var.push_back(h_el_real[year+"_"+LT_postfix]->GetBinContent(h_el_real[year+"_"+LT_postfix]->FindBin(pT2_real,eta2_real)) - h_el_real[year+"_"+LT_postfix]->GetBinError(h_el_real[year+"_"+LT_postfix]->FindBin(pT2_real,eta2_real)));

//     f2 = (hf_frac_nom*h_el_fake_hf[year+"_"+LT_postfix]->GetBinContent(h_el_fake_hf[year+"_"+LT_postfix]->FindBin(pT2_fake)) +
//     	  lf_frac_nom*h_el_fake_lf[year+"_"+LT_postfix]->GetBinContent(h_el_fake_lf[year+"_"+LT_postfix]->FindBin(pT2_fake)) +
// 	  co_frac_nom*h_el_fake_co[year+"_"+LT_postfix]->GetBinContent(h_el_fake_co[year+"_"+LT_postfix]->FindBin(pT2_fake_co))
//     	  );

//     f2_stat_var.push_back(hf_frac_nom*(h_el_fake_hf[year+"_"+LT_postfix]->GetBinContent(h_el_fake_hf[year+"_"+LT_postfix]->FindBin(pT2_fake)) + h_el_fake_hf[year+"_"+LT_postfix]->GetBinError(h_el_fake_hf[year+"_"+LT_postfix]->FindBin(pT2_fake))) +
// 			  lf_frac_nom*(h_el_fake_lf[year+"_"+LT_postfix]->GetBinContent(h_el_fake_lf[year+"_"+LT_postfix]->FindBin(pT2_fake)) + h_el_fake_lf[year+"_"+LT_postfix]->GetBinError(h_el_fake_lf[year+"_"+LT_postfix]->FindBin(pT2_fake))) +
// 			  co_frac_nom*(h_el_fake_co[year+"_"+LT_postfix]->GetBinContent(h_el_fake_co[year+"_"+LT_postfix]->FindBin(pT2_fake_co)) + h_el_fake_co[year+"_"+LT_postfix]->GetBinError(h_el_fake_co[year+"_"+LT_postfix]->FindBin(pT2_fake_co)))
// 			  );
//     f2_stat_var.push_back(hf_frac_nom*(h_el_fake_hf[year+"_"+LT_postfix]->GetBinContent(h_el_fake_hf[year+"_"+LT_postfix]->FindBin(pT2_fake)) - h_el_fake_hf[year+"_"+LT_postfix]->GetBinError(h_el_fake_hf[year+"_"+LT_postfix]->FindBin(pT2_fake))) +
// 			  lf_frac_nom*(h_el_fake_lf[year+"_"+LT_postfix]->GetBinContent(h_el_fake_lf[year+"_"+LT_postfix]->FindBin(pT2_fake)) - h_el_fake_lf[year+"_"+LT_postfix]->GetBinError(h_el_fake_lf[year+"_"+LT_postfix]->FindBin(pT2_fake))) +
// 			  co_frac_nom*(h_el_fake_co[year+"_"+LT_postfix]->GetBinContent(h_el_fake_co[year+"_"+LT_postfix]->FindBin(pT2_fake_co)) - h_el_fake_co[year+"_"+LT_postfix]->GetBinError(h_el_fake_co[year+"_"+LT_postfix]->FindBin(pT2_fake_co)))
// 			  );

//     for(unsigned int i = 0; i<hf.size(); i++){
//       f2_wgt_var.push_back((hf.at(i)*h_el_fake_hf[year+"_"+LT_postfix]->GetBinContent(h_el_fake_hf[year+"_"+LT_postfix]->FindBin(pT2_fake)) +
// 			    lf.at(i)*h_el_fake_lf[year+"_"+LT_postfix]->GetBinContent(h_el_fake_lf[year+"_"+LT_postfix]->FindBin(pT2_fake)) +
// 			    co.at(i)*h_el_fake_co[year+"_"+LT_postfix]->GetBinContent(h_el_fake_co[year+"_"+LT_postfix]->FindBin(pT2_fake_co))
// 			    ));
//     }

//     f2_xsecup_var = ((hf_frac_nom*h_el_fake_hf_xsecup[year+"_"+LT_postfix]->GetBinContent(h_el_fake_hf_xsecup[year+"_"+LT_postfix]->FindBin(pT2_fake)) +
// 		      lf_frac_nom*h_el_fake_lf[year+"_"+LT_postfix]->GetBinContent(h_el_fake_lf[year+"_"+LT_postfix]->FindBin(pT2_fake)) +
// 		      co_frac_nom*h_el_fake_co[year+"_"+LT_postfix]->GetBinContent(h_el_fake_co[year+"_"+LT_postfix]->FindBin(pT2_fake_co))
// 		      ));
//     f2_xsecdw_var = ((hf_frac_nom*h_el_fake_hf_xsecdw[year+"_"+LT_postfix]->GetBinContent(h_el_fake_hf_xsecdw[year+"_"+LT_postfix]->FindBin(pT2_fake)) +
// 		      lf_frac_nom*h_el_fake_lf[year+"_"+LT_postfix]->GetBinContent(h_el_fake_lf[year+"_"+LT_postfix]->FindBin(pT2_fake)) +
// 		      co_frac_nom*h_el_fake_co[year+"_"+LT_postfix]->GetBinContent(h_el_fake_co[year+"_"+LT_postfix]->FindBin(pT2_fake_co))
// 		      ));


//   }else if(lepFlavor[idx2] == 2){

//     // if(doRJR){
//     //   if(year.Contains("2017") && pT2_fake >= 60)year = "2018";
//     // }
//     if(LT_postfix.Contains("HP_")){
//       pT2_fake = lepPt[idx2] < 50 ? lepPt[idx2] : 49.0;
//     }else{
//       pT2_fake = lepPt[idx2] < 80 ? lepPt[idx2] : 79.0;
//     }


//     // if(doRJR){
//     //   if(year.Contains("2017"))year = "2018";
//     // }


//     r2 = h_mu_real[year+"_"+LT_postfix]->GetBinContent(h_mu_real[year+"_"+LT_postfix]->FindBin(pT2_real,eta2_real));

//     r2_mc = h_mu_real_mc[year+"_"+LT_postfix]->GetBinContent(h_mu_real_mc[year+"_"+LT_postfix]->FindBin(pT2_real,eta2_real));

//     r2_stat_var.push_back(h_mu_real[year+"_"+LT_postfix]->GetBinContent(h_mu_real[year+"_"+LT_postfix]->FindBin(pT2_real,eta2_real)) + h_mu_real[year+"_"+LT_postfix]->GetBinError(h_mu_real[year+"_"+LT_postfix]->FindBin(pT2_real,eta2_real)));
//     r2_stat_var.push_back(h_mu_real[year+"_"+LT_postfix]->GetBinContent(h_mu_real[year+"_"+LT_postfix]->FindBin(pT2_real,eta2_real)) - h_mu_real[year+"_"+LT_postfix]->GetBinError(h_mu_real[year+"_"+LT_postfix]->FindBin(pT2_real,eta2_real)));


//     f2 = h_mu_fake_hf[year+"_"+LT_postfix]->GetBinContent(h_mu_fake_hf[year+"_"+LT_postfix]->FindBin(pT2_fake));

//     f2_stat_var.push_back(h_mu_fake_hf[year+"_"+LT_postfix]->GetBinContent(h_mu_fake_hf[year+"_"+LT_postfix]->FindBin(pT2_fake)) + h_mu_fake_hf[year+"_"+LT_postfix]->GetBinError(h_mu_fake_hf[year+"_"+LT_postfix]->FindBin(pT2_fake)));
//     f2_stat_var.push_back(h_mu_fake_hf[year+"_"+LT_postfix]->GetBinContent(h_mu_fake_hf[year+"_"+LT_postfix]->FindBin(pT2_fake)) - h_mu_fake_hf[year+"_"+LT_postfix]->GetBinError(h_mu_fake_hf[year+"_"+LT_postfix]->FindBin(pT2_fake)));


//     f2_wgt_var.push_back(f2);

//     f2_xsecup_var = (h_mu_fake_hf_xsecup[year+"_"+LT_postfix]->GetBinContent(h_mu_fake_hf_xsecup[year+"_"+LT_postfix]->FindBin(pT2_fake)));
//     f2_xsecdw_var = (h_mu_fake_hf_xsecdw[year+"_"+LT_postfix]->GetBinContent(h_mu_fake_hf_xsecdw[year+"_"+LT_postfix]->FindBin(pT2_fake)));

//     year = year_const;

//   }
//   if(vb > 0){
//     printf("r1 = %.2f, f1 = %.2f\n",r1,f1);
//     printf("r2 = %.2f, f2 = %.2f\n",r2,f2);
//   }

//   r_wgt  = mtx->N4_RR_TT(r1,f1,r2,f2,nTT,nTl,nlT,nll);
//   rf_wgt = mtx->N4_RF_TT(r1,f1,r2,f2,nTT,nTl,nlT,nll);
//   fr_wgt = mtx->N4_FR_TT(r1,f1,r2,f2,nTT,nTl,nlT,nll);
//   f_wgt  = mtx->N4_FF_TT(r1,f1,r2,f2,nTT,nTl,nlT,nll);
//   for(unsigned int i = 0; i<r1_stat_var.size(); i++){

    

//     for(unsigned int j = 0; j<f1_stat_var.size(); j++){
//       if(vb > 0){
// 	printf("Stat variations: r1 = %.2f, f1 = %.2f, r2 = %.2f, f2 = %.2f\n",
// 	       r1_stat_var.at(i) < 1.0 ? r1_stat_var.at(i) : 0.99, 
// 	       f1_stat_var.at(j) > 0.0 ? f1_stat_var.at(j) : 0.01, 
// 	       r2_stat_var.at(i) < 1.0 ? r2_stat_var.at(i) : 0.99, 
// 	       f2_stat_var.at(j) > 0.0 ? f2_stat_var.at(j) : 0.01);
	  
//       }
//       r_wgt_stat.push_back( mtx->N4_RR_TT(r1_stat_var.at(i) < 1.0 ? r1_stat_var.at(i) : 0.99, 
// 					  f1_stat_var.at(j) > 0.0 ? f1_stat_var.at(j) : 0.01, 
// 					  r2_stat_var.at(i) < 1.0 ? r2_stat_var.at(i) : 0.99, 
// 					  f2_stat_var.at(j) > 0.0 ? f2_stat_var.at(j) : 0.01,
// 					  nTT,nTl,nlT,nll));
//       rf_wgt_stat.push_back(mtx->N4_RF_TT(r1_stat_var.at(i) < 1.0 ? r1_stat_var.at(i) : 0.99, 
// 					  f1_stat_var.at(j) > 0.0 ? f1_stat_var.at(j) : 0.01, 
// 					  r2_stat_var.at(i) < 1.0 ? r2_stat_var.at(i) : 0.99, 
// 					  f2_stat_var.at(j) > 0.0 ? f2_stat_var.at(j) : 0.01,
// 					  nTT,nTl,nlT,nll));
//       fr_wgt_stat.push_back(mtx->N4_FR_TT(r1_stat_var.at(i) < 1.0 ? r1_stat_var.at(i) : 0.99, 
// 					  f1_stat_var.at(j) > 0.0 ? f1_stat_var.at(j) : 0.01, 
// 					  r2_stat_var.at(i) < 1.0 ? r2_stat_var.at(i) : 0.99, 
// 					  f2_stat_var.at(j) > 0.0 ? f2_stat_var.at(j) : 0.01,
// 					  nTT,nTl,nlT,nll));
//       f_wgt_stat.push_back( mtx->N4_FF_TT(r1_stat_var.at(i) < 1.0 ? r1_stat_var.at(i) : 0.99, 
// 					  f1_stat_var.at(j) > 0.0 ? f1_stat_var.at(j) : 0.01, 
// 					  r2_stat_var.at(i) < 1.0 ? r2_stat_var.at(i) : 0.99, 
// 					  f2_stat_var.at(j) > 0.0 ? f2_stat_var.at(j) : 0.01,
// 					  nTT,nTl,nlT,nll));
//       fake_wgt_stat.push_back(rf_wgt_stat.back()+fr_wgt_stat.back()+f_wgt_stat.back());
//     }
//   }
//   for(unsigned int i = 0; i<f1_wgt_var.size(); i++){
//     for(unsigned int j = 0; j<f2_wgt_var.size(); j++){
//       r_wgt_syst.push_back( mtx->N4_RR_TT(r1,f1_wgt_var.at(i) > 0.0 ? f1_wgt_var.at(i) : 0.01, r2, f2_wgt_var.at(j) > 0.0 ? f2_wgt_var.at(j) : 0.01, nTT,nTl,nlT,nll));
//       rf_wgt_syst.push_back(mtx->N4_RF_TT(r1,f1_wgt_var.at(i) > 0.0 ? f1_wgt_var.at(i) : 0.01, r2, f2_wgt_var.at(j) > 0.0 ? f2_wgt_var.at(j) : 0.01, nTT,nTl,nlT,nll));
//       fr_wgt_syst.push_back(mtx->N4_FR_TT(r1,f1_wgt_var.at(i) > 0.0 ? f1_wgt_var.at(i) : 0.01, r2, f2_wgt_var.at(j) > 0.0 ? f2_wgt_var.at(j) : 0.01, nTT,nTl,nlT,nll));
//       f_wgt_syst.push_back( mtx->N4_FF_TT(r1,f1_wgt_var.at(i) > 0.0 ? f1_wgt_var.at(i) : 0.01, r2, f2_wgt_var.at(j) > 0.0 ? f2_wgt_var.at(j) : 0.01, nTT,nTl,nlT,nll));
//       fake_wgt_syst.push_back(rf_wgt_syst.back()+fr_wgt_syst.back()+f_wgt_syst.back());
//     }
//   }

//   r_wgt_xsec.push_back( mtx->N4_RR_TT(r1, f1_xsecup_var > 0.0 ? f1_xsecup_var: 0.01, r2, f2_xsecup_var > 0.0 ? f2_xsecup_var : 0.01, nTT,nTl,nlT,nll));
//   rf_wgt_xsec.push_back(mtx->N4_RF_TT(r1, f1_xsecup_var > 0.0 ? f1_xsecup_var: 0.01, r2, f2_xsecup_var > 0.0 ? f2_xsecup_var : 0.01, nTT,nTl,nlT,nll));
//   fr_wgt_xsec.push_back(mtx->N4_FR_TT(r1, f1_xsecup_var > 0.0 ? f1_xsecup_var: 0.01, r2, f2_xsecup_var > 0.0 ? f2_xsecup_var : 0.01, nTT,nTl,nlT,nll));
//   f_wgt_xsec.push_back( mtx->N4_FF_TT(r1, f1_xsecup_var > 0.0 ? f1_xsecup_var: 0.01, r2, f2_xsecup_var > 0.0 ? f2_xsecup_var : 0.01, nTT,nTl,nlT,nll));
//   fake_wgt_xsec.push_back(rf_wgt_xsec.back()+fr_wgt_xsec.back()+f_wgt_xsec.back());

//   r_wgt_xsec.push_back( mtx->N4_RR_TT(r1,f1_xsecdw_var,r2,f2_xsecdw_var,nTT,nTl,nlT,nll));
//   rf_wgt_xsec.push_back(mtx->N4_RF_TT(r1,f1_xsecdw_var,r2,f2_xsecdw_var,nTT,nTl,nlT,nll));
//   fr_wgt_xsec.push_back(mtx->N4_FR_TT(r1,f1_xsecdw_var,r2,f2_xsecdw_var,nTT,nTl,nlT,nll));
//   f_wgt_xsec.push_back( mtx->N4_FF_TT(r1,f1_xsecdw_var,r2,f2_xsecdw_var,nTT,nTl,nlT,nll));
//   fake_wgt_xsec.push_back(rf_wgt_xsec.back()+fr_wgt_xsec.back()+f_wgt_xsec.back());

//   r_wgt_mc  = mtx->N4_RR_TT(r1_mc < 1.0 ? r1_mc : 0.99, f1 > 0.0 ? f1 : 0.01, r2_mc < 1.0 ? r2_mc : 0.99, f2 > 0.0 ? f2 : 0.01, nTT,nTl,nlT,nll);
//   rf_wgt_mc = mtx->N4_RF_TT(r1_mc < 1.0 ? r1_mc : 0.99, f1 > 0.0 ? f1 : 0.01, r2_mc < 1.0 ? r2_mc : 0.99, f2 > 0.0 ? f2 : 0.01, nTT,nTl,nlT,nll);
//   fr_wgt_mc = mtx->N4_FR_TT(r1_mc < 1.0 ? r1_mc : 0.99, f1 > 0.0 ? f1 : 0.01, r2_mc < 1.0 ? r2_mc : 0.99, f2 > 0.0 ? f2 : 0.01, nTT,nTl,nlT,nll);
//   f_wgt_mc  = mtx->N4_FF_TT(r1_mc < 1.0 ? r1_mc : 0.99, f1 > 0.0 ? f1 : 0.01, r2_mc < 1.0 ? r2_mc : 0.99, f2 > 0.0 ? f2 : 0.01, nTT,nTl,nlT,nll);
//   fake_wgt_mc  = rf_wgt_mc+fr_wgt_mc+f_wgt_mc;

//   fake_wgt = rf_wgt+fr_wgt+f_wgt;

//   if(vb > 0){
//     printf("wgt(r)    =  %.4f\n",r_wgt);
//     printf("wgt(rf)   =  %.4f\n",rf_wgt);
//     printf("wgt(fr)   =  %.4f\n",fr_wgt);
//     printf("wgt(ff)   =  %.4f\n",f_wgt);
//     printf("wgt(fake) =  %f\n",rf_wgt+fr_wgt+f_wgt);


//   }
//   double frac_fake_mc;
//   double max_fake_syst;
//   double min_fake_syst;
//   double max_fake_stat;
//   double min_fake_stat;
//   double max_fake_wgt_syst;
//   double min_fake_wgt_syst;
//   double max_fake_wgt_stat = max(*max_element(std::begin(fake_wgt_stat), std::end(fake_wgt_stat)),*min_element(std::begin(fake_wgt_stat), std::end(fake_wgt_stat)));
//   double min_fake_wgt_stat = min(*max_element(std::begin(fake_wgt_stat), std::end(fake_wgt_stat)),*min_element(std::begin(fake_wgt_stat), std::end(fake_wgt_stat)));
//   double max_fake_wgt_xsec = max(*max_element(std::begin(fake_wgt_xsec), std::end(fake_wgt_xsec)),*min_element(std::begin(fake_wgt_xsec), std::end(fake_wgt_xsec)));
//   double min_fake_wgt_xsec = min(*max_element(std::begin(fake_wgt_xsec), std::end(fake_wgt_xsec)),*min_element(std::begin(fake_wgt_xsec), std::end(fake_wgt_xsec)));
//   if(!(eventStr.EqualTo("MMOS") || eventStr.EqualTo("MMSS"))){
//     max_fake_wgt_syst = max(*max_element(std::begin(fake_wgt_syst), std::end(fake_wgt_syst)),*min_element(std::begin(fake_wgt_syst), std::end(fake_wgt_syst)));
//     min_fake_wgt_syst = min(*max_element(std::begin(fake_wgt_syst), std::end(fake_wgt_syst)),*min_element(std::begin(fake_wgt_syst), std::end(fake_wgt_syst)));
//   }else{
//     max_fake_wgt_syst = rf_wgt+fr_wgt+f_wgt;
//     min_fake_wgt_syst = rf_wgt+fr_wgt+f_wgt;
//   }


//   if(vb > 0){

//     if(max_fake_wgt_syst < fake_wgt)cout<<"ERROR \t Max fake wgt syst "<< max_fake_wgt_syst << " is smaller than nominal "<<fake_wgt<< " for event type "<<fulleventStr.Data()<<endl;
//     if(min_fake_wgt_syst > fake_wgt)cout<<"ERROR \t Min fake wgt syst "<< min_fake_wgt_syst << " is larger  than nominal "<<fake_wgt<< " for event type "<<fulleventStr.Data()<<endl;
//     if(max_fake_wgt_stat < fake_wgt)cout<<"ERROR \t Max fake wgt stat "<< max_fake_wgt_stat << " is smaller than nominal "<<fake_wgt<< " for event type "<<fulleventStr.Data()<<endl;
//     if(min_fake_wgt_stat > fake_wgt)cout<<"ERROR \t Min fake wgt stat "<< min_fake_wgt_stat << " is larger  than nominal "<<fake_wgt<< " for event type "<<fulleventStr.Data()<<endl;


//     printf("wgt_syst_up(fake) =  %f\n",max_fake_wgt_syst);
//     printf("wgt_syst_dw(fake) =  %f\n",min_fake_wgt_syst);

//     printf("wgt_stat_up(fake) =  %f\n",max_fake_wgt_stat);
//     printf("wgt_syst_dw(fake) =  %f\n",min_fake_wgt_stat);

//     printf("wgt_MC(fake)   =  %f\n",fake_wgt_mc);

//  }
//   if(fake_wgt_mc >= fake_wgt){
//     max_fake_syst = sqrt(pow(fabs((max_fake_wgt_syst - fake_wgt)/fake_wgt),2) + pow(fabs((fake_wgt_mc - fake_wgt)/fake_wgt),2) + pow(fabs((max_fake_wgt_xsec - fake_wgt)/fake_wgt),2));
//     min_fake_syst = sqrt(pow(fabs((fake_wgt - min_fake_wgt_syst)/fake_wgt),2) + pow(fabs((max_fake_wgt_xsec - fake_wgt)/fake_wgt),2));
//     frac_fake_mc = fabs((fake_wgt_mc - fake_wgt)/fake_wgt);
//   }else{
//     min_fake_syst = sqrt(pow(fabs((fake_wgt - min_fake_wgt_syst)/fake_wgt),2) + pow(fabs((fake_wgt - fake_wgt_mc)/fake_wgt),2) + pow(fabs((fake_wgt - min_fake_wgt_xsec)/fake_wgt),2));
//     max_fake_syst = sqrt(pow(fabs((max_fake_wgt_syst - fake_wgt)/fake_wgt),2) + pow(fabs((max_fake_wgt_xsec - fake_wgt)/fake_wgt),2));
//     frac_fake_mc = fabs((fake_wgt - fake_wgt_mc)/fake_wgt);
//   }

//   max_fake_stat = fabs((fake_wgt - min_fake_wgt_stat)/fake_wgt);
//   min_fake_stat = fabs((max_fake_wgt_stat - fake_wgt)/fake_wgt);

//   h_syst_MM[fulleventStr+"_stat_up"]->SetBinContent(h_syst_MM[fulleventStr+"_stat_up"]->FindBin(lepPt[idx1],lepPt[idx2]),max_fake_stat*100.);
//   h_syst_MM[fulleventStr+"_stat_dw"]->SetBinContent(h_syst_MM[fulleventStr+"_stat_dw"]->FindBin(lepPt[idx1],lepPt[idx2]),min_fake_stat*100.);

//   h_syst_MM[fulleventStr+"_wgt_up"]->SetBinContent(h_syst_MM[fulleventStr+"_wgt_up"]->FindBin(lepPt[idx1],lepPt[idx2]),fabs((max_fake_wgt_syst - fake_wgt)/fake_wgt)*100.);
//   h_syst_MM[fulleventStr+"_wgt_dw"]->SetBinContent(h_syst_MM[fulleventStr+"_wgt_dw"]->FindBin(lepPt[idx1],lepPt[idx2]),fabs((fake_wgt - min_fake_wgt_syst)/fake_wgt)*100.);

//   if(fake_wgt_mc >= fake_wgt)
//     h_syst_MM[fulleventStr+"_mc_up"]->SetBinContent(h_syst_MM[fulleventStr+"_mc_up"]->FindBin(lepPt[idx1],lepPt[idx2]),frac_fake_mc*100.);
//   else
//     h_syst_MM[fulleventStr+"_mc_dw"]->SetBinContent(h_syst_MM[fulleventStr+"_mc_dw"]->FindBin(lepPt[idx1],lepPt[idx2]),frac_fake_mc*100.);

//   h_syst_MM[fulleventStr+"_xsec_up"]->SetBinContent(h_syst_MM[fulleventStr+"_xsec_up"]->FindBin(lepPt[idx1],lepPt[idx2]),fabs((max_fake_wgt_xsec - fake_wgt)/fake_wgt)*100.);
//   h_syst_MM[fulleventStr+"_xsec_dw"]->SetBinContent(h_syst_MM[fulleventStr+"_xsec_dw"]->FindBin(lepPt[idx1],lepPt[idx2]),fabs((fake_wgt - min_fake_wgt_xsec)/fake_wgt)*100.);


//   double frac_total_up = 1.0;
//   double frac_total_dw = 1.0;


//   if(lepFlavor[idx1] == 2 && lepFlavor[idx2] == 2){
//     frac_total_up = sqrt(pow(max_fake_syst,2) + pow(max_fake_stat,2));
//     frac_total_dw = sqrt(pow(min_fake_syst,2) + pow(min_fake_stat,2));
//   }else{
//     frac_total_up = sqrt(pow(max_fake_syst,2) + pow(max_fake_stat,2));
//     frac_total_dw = sqrt(pow(min_fake_syst,2) + pow(min_fake_stat,2));
//   }

//   h_syst_MM[fulleventStr+"_total_up"]->SetBinContent(h_syst_MM[fulleventStr+"_total_up"]->FindBin(lepPt[idx1],lepPt[idx2]),frac_total_up*100);
//   h_syst_MM[fulleventStr+"_total_dw"]->SetBinContent(h_syst_MM[fulleventStr+"_total_dw"]->FindBin(lepPt[idx1],lepPt[idx2]),frac_total_dw*100);


//   if(vb > 0){

//     printf("frac_syst_up(fake) =  %f\n",fabs((max_fake_wgt_syst - fake_wgt)/fake_wgt));
//     printf("frac_syst_dw(fake) =  %f\n",fabs((fake_wgt - min_fake_wgt_syst)/fake_wgt));

//     printf("frac_stat_up(fake) =  %f\n",max_fake_stat);
//     printf("frac_syst_dw(fake) =  %f\n",min_fake_stat);

//     printf("frac_MC(fake)   =  %f\n",frac_fake_mc);

//     printf("frac_syst_up_total(fake)   =  %f\n",max_fake_syst);
//     printf("frac_syst_dw_total(fake)   =  %f\n",min_fake_syst);

//     printf("frac_up_total(fake)   =  %f\n",frac_total_up);
//     printf("frac_dw_total(fake)   =  %f\n",frac_total_dw);

//  }


//   if(fabs(rf_wgt+fr_wgt+f_wgt)>2.0){
//     printf("a) Year is %s, Events is %s\n",year.Data(),fulleventStr.Data());
//     printf("Lep1 is %i: pT = %.2f (f, %.2f; r, %.2f), |eta| = %.2f\n",lepFlavor[idx1],lepPt[idx1],pT1_fake,pT1_real,eta1_real);
//     printf("Lep2 is %i: pT = %.2f (f, %.2f; r, %.2f), |eta| = %.2f\n",lepFlavor[idx2],lepPt[idx2],pT2_fake,pT2_real,eta2_real);
//     printf("TT = %i, TL = %i, LT = %i, LL = %i, \n",nTT, nTl, nlT, nll);
//     if(vb > 0 && !(eventStr.EqualTo("MMOS") || eventStr.EqualTo("MMSS"))){
//       printf("mt2 = %.2f, mt2bin = %i, wgt(hf) = %.2f + %.2f - %.2f, wgt(lf) = %.2f + %.2f - %.2f, wgt(co) = %.2f + %.2f - %.2f, sum = %.2f + %.2f - %.2f\n",var_dep,mt2bin,
// 	     hf_frac_nom,hf_frac_max,hf_frac_min,
// 	     lf_frac_nom,lf_frac_max,lf_frac_min,
// 	     co_frac_nom,co_frac_max,co_frac_min,
// 	     hf_frac_nom+lf_frac_nom+co_frac_nom, hf_frac_max+lf_frac_max+co_frac_max, hf_frac_min+lf_frac_min+co_frac_min);
//     }
//     printf("r1 = %.2f, f1 = %.2f\n",r1,f1);
//     printf("r2 = %.2f, f2 = %.2f\n",r2,f2);
//     printf("wgt(r)    =  %.4f\n",r_wgt);
//     printf("wgt(rf)   =  %.4f\n",rf_wgt);
//     printf("wgt(fr)   =  %.4f\n",fr_wgt);
//     printf("wgt(ff)   =  %.4f\n",f_wgt);
//     printf("wgt(fake) =  %.4f\n",rf_wgt+fr_wgt+f_wgt);
//   }

//   std::vector<double> retvec;
//   double fake = rf_wgt+fr_wgt+f_wgt;
//   retvec.push_back(fake);
//   retvec.push_back(fake > 0 ? (1+frac_total_up) : (1-frac_total_up));
//   retvec.push_back(fake > 0 ? (1-frac_total_dw) : (1+frac_total_dw));

//   if(vb > 0){
//     for(unsigned int i = 0; i < retvec.size(); i++){
//       cout<<"Wgt("<<i<<") = "<<retvec.at(i)<<endl;
//     }
//   }

//   // setting some global variables (for ntuple)
//   glob_r1 = r1;
//   glob_r2 = r2;
//   glob_f1 = f1;
//   glob_f2 = f2;
  

//   return retvec;

// }

std::vector<double> MySusySkimAnalysis::getFakeWeight(int idx1, int idx2, int idx0, TString year, TString LT_postfix, int vb)
{

  if(doThreeLep){return getFakeWeight1D(idx1,idx2,idx0,year,LT_postfix,vb);}

  int nTT = 0;
  int nTl = 0;
  int nlT = 0;
  int nll = 0;
  double r1 = 0.0;
  double f1 = 0.0;
  double r2 = 0.0;
  double f2 = 0.0;
  TString eventStr = "";
  double r_wgt  = 0.0;
  double rf_wgt = 0.0;
  double fr_wgt = 0.0;
  double f_wgt  = 0.0;

  classifyEvent(idx1, idx2, nTT, nTl, nlT, nll);
  std::pair <TString,TString> ret = getFillString(idx1, idx2, idx0);
  eventStr = ret.first + ret.second + "_" + LT_postfix;

  // a few hacks
  double pT1_fake = lepPt[idx1] < 200 ? lepPt[idx1] : 199.0;
  double pT1_real = pT1_fake;
  double pT2_fake = lepPt[idx2] < 200 ? lepPt[idx2] : 199.0;
  double pT2_real = pT2_fake;
  
  if(vb > 0){
    printf("b) Year is %s, Events is %s\n",year.Data(),eventStr.Data());
    printf("Lep1 is %i: pT = %.2f (%.2f), |eta| = %.2f\n",lepFlavor[idx1], lepPt[idx1],pT1_fake,fabs(lepEta[idx1]));
    printf("Lep2 is %i: pT = %.2f (%.2f), |eta| = %.2f\n",lepFlavor[idx2], lepPt[idx2],pT2_fake,fabs(lepEta[idx2]));
    printf("TT = %i, TL = %i, LT = %i, LL = %i, \n",nTT, nTl, nlT, nll);
  }

  int mt2bin = 1;
  if(themt2 >= 0 && themt2 < 20)mt2bin = 1;
  else if(themt2 >= 20 && themt2 < 40)mt2bin = 2;
  else if(themt2 >= 40 && themt2 < 60)mt2bin = 3;
  else if(themt2 >= 60 && themt2 < 100)mt2bin = 4;
  else if(themt2 >= 100)mt2bin = 5;

  if(vb > 0 && (eventStr.Contains("EEOS") || eventStr.Contains("EMOS") || eventStr.Contains("EESS") || eventStr.Contains("EMSS"))){
    printf("mt2 = %.2f, mt2bin = %i, wgt(hf) = %.2f, wgt(cf) = %.2f, wgt(lf) = %.2f, wgt(cf) = %.2f, wgt(co) = %.2f, sum = %.2f\n",themt2,mt2bin,
	   h_el_frac_hf[eventStr]->GetBinContent(mt2bin),
	   h_el_frac_cf[eventStr]->GetBinContent(mt2bin),
	   h_el_frac_lf[eventStr]->GetBinContent(mt2bin),
	   h_el_frac_co[eventStr]->GetBinContent(mt2bin),
	   h_el_frac_cf[eventStr]->GetBinContent(mt2bin),
	   h_el_frac_hf[eventStr]->GetBinContent(mt2bin)+h_el_frac_cf[eventStr]->GetBinContent(mt2bin)+h_el_frac_lf[eventStr]->GetBinContent(mt2bin)+h_el_frac_co[eventStr]->GetBinContent(mt2bin)+h_el_frac_cf[eventStr]->GetBinContent(mt2bin));
  }


  if(lepFlavor[idx1] == 1){

    // if(year.Contains("2015_2016") && pT1_fake >= 30 && pT1_fake < 35 && fabs(lepEta[idx1]) >= 2.4) pT1_fake = 29.0;
    // if(year.Contains("2017")      && pT1_fake >= 35 && pT1_fake < 40 && fabs(lepEta[idx1]) >= 0.4 && fabs(lepEta[idx1]) < 0.6) pT1_fake = 34.0;

    r1 = h_el_real[year+"_"+LT_postfix]->GetBinContent(h_el_real[year+"_"+LT_postfix]->FindBin(pT1_real,fabs(lepEta[idx1])));
    f1 = ((h_el_frac_hf[eventStr]->GetBinContent(mt2bin)+h_el_frac_cf[eventStr]->GetBinContent(mt2bin))*h_el_fake_hf[year+"_"+LT_postfix]->GetBinContent(h_el_fake_hf[year+"_"+LT_postfix]->FindBin(pT1_fake,fabs(lepEta[idx1]))) +
	  h_el_frac_lf[eventStr]->GetBinContent(mt2bin)*h_el_fake_lf[year+"_"+LT_postfix]->GetBinContent(h_el_fake_lf[year+"_"+LT_postfix]->FindBin(pT1_fake,fabs(lepEta[idx1]))) +
	  (h_el_frac_co[eventStr]->GetBinContent(mt2bin)+h_el_frac_cf[eventStr]->GetBinContent(mt2bin))*h_el_fake_co[year+"_"+LT_postfix]->GetBinContent(h_el_fake_co[year+"_"+LT_postfix]->FindBin(pT1_fake,0))
	  );
  }else if(lepFlavor[idx1] == 2){
    r1 = h_mu_real[year+"_"+LT_postfix]->GetBinContent(h_mu_real[year+"_"+LT_postfix]->FindBin(pT1_real,fabs(lepEta[idx1])));
    f1 = h_mu_fake_hf[year+"_"+LT_postfix]->GetBinContent(h_mu_fake_hf[year+"_"+LT_postfix]->FindBin(pT1_fake));// < 80 ? pT1_fake : 79,0));
  }
  if(lepFlavor[idx2] == 1){


    // if(year.Contains("2015_2016") && pT2_fake >= 30 && pT2_fake < 35 && fabs(lepEta[idx2]) >= 2.4) pT2_fake = 29.0;
    // if(year.Contains("2017")      && pT2_fake >= 35 && pT2_fake < 40 && fabs(lepEta[idx2]) >= 0.4 && fabs(lepEta[idx2]) < 0.6) pT2_fake = 34.0;

    r2 = h_el_real[year+"_"+LT_postfix]->GetBinContent(h_el_real[year+"_"+LT_postfix]->FindBin(pT2_real,fabs(lepEta[idx2])));
    f2 = ((h_el_frac_hf[eventStr]->GetBinContent(mt2bin)+h_el_frac_cf[eventStr]->GetBinContent(mt2bin))*h_el_fake_hf[year+"_"+LT_postfix]->GetBinContent(h_el_fake_hf[year+"_"+LT_postfix]->FindBin(pT2_fake,fabs(lepEta[idx2]))) +
	  h_el_frac_lf[eventStr]->GetBinContent(mt2bin)*h_el_fake_lf[year+"_"+LT_postfix]->GetBinContent(h_el_fake_lf[year+"_"+LT_postfix]->FindBin(pT2_fake,fabs(lepEta[idx2]))) +
	  (h_el_frac_co[eventStr]->GetBinContent(mt2bin)+h_el_frac_cf[eventStr]->GetBinContent(mt2bin))*h_el_fake_co[year+"_"+LT_postfix]->GetBinContent(h_el_fake_co[year+"_"+LT_postfix]->FindBin(pT2_fake,0))
	  );
  }else if(lepFlavor[idx2] == 2){
    r2 = h_mu_real[year+"_"+LT_postfix]->GetBinContent(h_mu_real[year+"_"+LT_postfix]->FindBin(pT2_real,fabs(lepEta[idx1])));
    f2 = h_mu_fake_hf[year+"_"+LT_postfix]->GetBinContent(h_mu_fake_hf[year+"_"+LT_postfix]->FindBin(pT2_fake));// < 80 ? pT2_fake : 79,0));
  }

  if(vb > 0){
    printf("r1 = %.2f, f1 = %.2f\n",r1,f1);
    printf("r2 = %.2f, f2 = %.2f\n",r2,f2);
  }

  r_wgt  = mtx->N4_RR_TT(r1,f1,r2,f2,nTT,nTl,nlT,nll);
  rf_wgt = mtx->N4_RF_TT(r1,f1,r2,f2,nTT,nTl,nlT,nll);
  fr_wgt = mtx->N4_FR_TT(r1,f1,r2,f2,nTT,nTl,nlT,nll);
  f_wgt  = mtx->N4_FF_TT(r1,f1,r2,f2,nTT,nTl,nlT,nll);

  if(vb > 0){
    printf("wgt(r)    =  %.4f\n",r_wgt);
    printf("wgt(rf)   =  %.4f\n",rf_wgt);
    printf("wgt(fr)   =  %.4f\n",fr_wgt);
    printf("wgt(ff)   =  %.4f\n",f_wgt);
    printf("wgt(fake) =  %.4f\n",rf_wgt+fr_wgt+f_wgt);
  }

  if(fabs(rf_wgt+fr_wgt+f_wgt)>10){
    printf("Year is %s, Events is %s\n",year.Data(),eventStr.Data());
    printf("Lep1 is %i: pT = %.2f (f, %.2f; r, %.2f), |eta| = %.2f\n",lepFlavor[idx1],lepPt[idx1],pT1_fake,pT1_real,fabs(lepEta[idx1]));
    printf("Lep2 is %i: pT = %.2f (f, %.2f; r, %.2f), |eta| = %.2f\n",lepFlavor[idx2],lepPt[idx2],pT2_fake,pT2_real,fabs(lepEta[idx2]));
    printf("TT = %i, TL = %i, LT = %i, LL = %i, \n",nTT, nTl, nlT, nll);
    if(vb > 0 && (eventStr.Contains("EEOS") || eventStr.Contains("EMOS") || eventStr.Contains("EESS") || eventStr.Contains("EMSS"))){
      printf("mt2 = %.2f, mt2bin = %i, wgt(hf) = %.2f, wgt(cf) = %.2f, wgt(lf) = %.2f, wgt(cf) = %.2f, wgt(co) = %.2f, sum = %.2f\n",themt2,mt2bin,
	     h_el_frac_hf[eventStr]->GetBinContent(mt2bin),
	     h_el_frac_cf[eventStr]->GetBinContent(mt2bin),
	     h_el_frac_lf[eventStr]->GetBinContent(mt2bin),
	     h_el_frac_co[eventStr]->GetBinContent(mt2bin),
	     h_el_frac_cf[eventStr]->GetBinContent(mt2bin),
	     h_el_frac_hf[eventStr]->GetBinContent(mt2bin)+h_el_frac_cf[eventStr]->GetBinContent(mt2bin)+h_el_frac_lf[eventStr]->GetBinContent(mt2bin)+h_el_frac_co[eventStr]->GetBinContent(mt2bin)+h_el_frac_cf[eventStr]->GetBinContent(mt2bin));
    }
    printf("r1 = %.2f, f1 = %.2f\n",r1,f1);
    printf("r2 = %.2f, f2 = %.2f\n",r2,f2);
    printf("wgt(r)    =  %.4f\n",r_wgt);
    printf("wgt(rf)   =  %.4f\n",rf_wgt);
    printf("wgt(fr)   =  %.4f\n",fr_wgt);
    printf("wgt(ff)   =  %.4f\n",f_wgt);
    printf("wgt(fake) =  %.4f\n",rf_wgt+fr_wgt+f_wgt);
  }

  std::vector<double> retvec;
  double fake = rf_wgt+fr_wgt+f_wgt;
  retvec.push_back(fake);

  return retvec;

}

bool MySusySkimAnalysis::isBLjet(int idx)
{

  bool is2L2J       = false;
  bool isZMET       = false;
  bool isFCT        = false;
  bool isODA        = false;
  bool isFCL        = false;
  bool isFCL_HPTC   = false;
  bool isGRAD      = false;
  bool isFCL_FCTT   = false;
  bool isFCL_FCHPT  = false;

  // for Helen
  bool isHP_MT_FCL_FCHPT = false;
  bool isHP_0M_FCL_FCHPT = false;
  bool isHP_MT_FCT = false;
  bool isHP_0M_FCT = false;

  if(LT_postfix.EqualTo("2L2J"))is2L2J = true;
  else if(LT_postfix.EqualTo("ZMET"))isZMET = true;
  else if(LT_postfix.EqualTo("HP_MT_FCL_FCHPT"))isHP_MT_FCL_FCHPT = true;
  else if(LT_postfix.EqualTo("HP_0M_FCL_FCHPT"))isHP_0M_FCL_FCHPT = true;
  else if(LT_postfix.EqualTo("HP_MT_FCT"))isHP_MT_FCT = true;
  else if(LT_postfix.EqualTo("HP_0M_FCT"))isHP_0M_FCT = true;
  else if(LT_postfix.EqualTo("FCT"))isFCT = true;
  else if(LT_postfix.EqualTo("FCT_MT"))isFCT = true;
  else if(LT_postfix.EqualTo("ODA"))isODA = true;
  else if(LT_postfix.EqualTo("FCL") || LT_postfix.EqualTo("notrigm") || LT_postfix.EqualTo("med") || LT_postfix.EqualTo("loose") || LT_postfix.EqualTo("tight_iso") || LT_postfix.EqualTo("no"))isFCL = true;
  else if(LT_postfix.EqualTo("FCL_HPTC"))isFCL_HPTC = true;
  else if(LT_postfix.EqualTo("GRAD"))isGRAD = true;
  else if(LT_postfix.EqualTo("FCL_FCTT"))isFCL_FCTT = true;
  else if(LT_postfix.EqualTo("FCL_FCHPT"))isFCL_FCHPT = true;
  else {
    cout<<"ERROR \t Unknown tight-loose postfix "<<LT_postfix.Data()<<endl;
    return false;
  }

  if(jetPt[idx]<=20)return false;

  if(is2L2J){
    if(fabs(jetEta[idx])>=2.8)return false;
  }else if(isODA && fabs(jetEta[idx])>=2.5)return false;
  else if(isZMET && fabs(jetEta[idx])>=4.5)return false;
  else if(!is2L2J && !isODA && fabs(jetEta[idx])>=2.8)return false;
  return true;

}



bool MySusySkimAnalysis::isSGjet(int idx)
{

  bool is2L2J       = false;
  bool isZMET       = false;
  bool isFCT        = false;
  bool isODA        = false;
  bool isFCL        = false;
  bool isFCL_HPTC   = false;
  bool isGRAD      = false;
  bool isFCL_FCTT   = false;
  bool isFCL_FCHPT  = false;
  // for Helen
  bool isHP_MT_FCL_FCHPT = false;
  bool isHP_0M_FCL_FCHPT = false;
  bool isHP_MT_FCT = false;
  bool isHP_0M_FCT = false;

  if(LT_postfix.EqualTo("2L2J"))is2L2J = true;
  else if(LT_postfix.EqualTo("ZMET"))isZMET = true;
  else if(LT_postfix.EqualTo("HP_MT_FCL_FCHPT"))isHP_MT_FCL_FCHPT = true;
  else if(LT_postfix.EqualTo("HP_0M_FCL_FCHPT"))isHP_0M_FCL_FCHPT = true;
  else if(LT_postfix.EqualTo("HP_MT_FCT"))isHP_MT_FCT = true;
  else if(LT_postfix.EqualTo("HP_0M_FCT"))isHP_0M_FCT = true;
  else if(LT_postfix.EqualTo("FCT"))isFCT = true;
  else if(LT_postfix.EqualTo("FCT_MT"))isFCT = true;
  else if(LT_postfix.EqualTo("ODA"))isODA = true;
  else if(LT_postfix.EqualTo("FCL") || LT_postfix.EqualTo("notrigm") || LT_postfix.EqualTo("med") || LT_postfix.EqualTo("loose") || LT_postfix.EqualTo("tight_iso") || LT_postfix.EqualTo("no"))isFCL = true;
  else if(LT_postfix.EqualTo("FCL_HPTC"))isFCL_HPTC = true;
  else if(LT_postfix.EqualTo("GRAD"))isGRAD = true;
  else if(LT_postfix.EqualTo("FCL_FCTT"))isFCL_FCTT = true;
  else if(LT_postfix.EqualTo("FCL_FCHPT"))isFCL_FCHPT = true;
  else {
    cout<<"ERROR \t Unknown tight-loose postfix "<<LT_postfix.Data()<<endl;
    return false;
  }

  if(!isBLjet(idx))return false;
  if(!is2L2J && !jetSignal->at(idx))return false; //<------------------------------------- COMMENT BACK IN. NOT AVAILABLE IN 2L2J OFFICIIAL NTUPLES
  if(!jetPassOR->at(idx))return false; //<------------------------------------- COMMENT BACK IN. NOT AVAILABLE IN 2L2J OFFICIIAL NTUPLES
  if((!is2L2J && !isZMET) && fabs(jetEta[idx])>=2.4)return false;
  //else if(isFCT && fabs(jetEta[idx])>=2.5)return false;
  
  //if(fabs(jetEta[idx])>=2.4)return false;
  if(is2L2J && jetPt[idx]<50 && fabs(jetEta[idx]) < 2.5 && fabs(jetJVT[idx]) <= 0.59)return false; //<------------------------------------- COMMENT BACK IN. NOT AVAILABLE IN 2L2J OFFICIIAL NTUPLES
  if(is2L2J){
    //return true;
    if(fabs(jetEta[idx])<2.4 && jetPt[idx] < 120 && jetJVT[idx] <= 0.59)return false;  //<------------------------------------- COMMENT BACK IN. NOT AVAILABLE IN 2L2J OFFICIIAL NTUPLES
    if(fabs(jetEta[idx])>=2.4 && jetPt[idx] < 120 && jetJVT[idx] <= 0.11)return false; //<------------------------------------- COMMENT BACK IN. NOT AVAILABLE IN 2L2J OFFICIIAL NTUPLES
  }
  if((!is2L2J && !isZMET) && jetPt[idx]<60 && fabs(jetJVT[idx]) <= 0.91)return false; //<------------------------------------- COMMENT BACK IN. NOT AVAILABLE IN 2L2J OFFICIIAL NTUPLES

  if(isZMET && jetPt[idx]<60 && fabs(jetEta[idx])<2.4 && fabs(jetJVT[idx]) <= 0.5)return false;
  //  if(isZMET && jetPt[idx]<120 && (fabs(jetEta[idx])>2.5 && fabs(jetEta[idx]) < 4.5) && fabs(jetJVT[idx]) <= 0.5)return false;
  return true;
  
}

bool MySusySkimAnalysis::isBjet(int idx)
{

  bool is2L2J       = false;
  bool isZMET       = false;
  bool isFCT        = false;
  bool isODA        = false;
  bool isFCL        = false;
  bool isFCL_HPTC   = false;
  bool isGRAD      = false;
  bool isFCL_FCTT   = false;
  bool isFCL_FCHPT  = false;
// for Helen
  bool isHP_MT_FCL_FCHPT = false;
  bool isHP_0M_FCL_FCHPT = false;
  bool isHP_MT_FCT = false;
  bool isHP_0M_FCT = false;
  if(LT_postfix.EqualTo("2L2J"))is2L2J = true;
  else if(LT_postfix.EqualTo("ZMET"))isZMET = true;
  else if(LT_postfix.EqualTo("HP_MT_FCL_FCHPT"))isHP_MT_FCL_FCHPT = true;
  else if(LT_postfix.EqualTo("HP_0M_FCL_FCHPT"))isHP_0M_FCL_FCHPT = true;
  else if(LT_postfix.EqualTo("HP_MT_FCT"))isHP_MT_FCT = true;
  else if(LT_postfix.EqualTo("HP_0M_FCT"))isHP_0M_FCT = true;
  else if(LT_postfix.EqualTo("FCT"))isFCT = true;
  else if(LT_postfix.EqualTo("FCT_MT"))isFCT = true;
  else if(LT_postfix.EqualTo("ODA"))isODA = true;
  else if(LT_postfix.EqualTo("FCL") || LT_postfix.EqualTo("notrigm") || LT_postfix.EqualTo("med") || LT_postfix.EqualTo("loose") || LT_postfix.EqualTo("tight_iso") || LT_postfix.EqualTo("no"))isFCL = true;
  else if(LT_postfix.EqualTo("FCL_HPTC"))isFCL_HPTC = true;
  else if(LT_postfix.EqualTo("GRAD"))isGRAD = true;
  else if(LT_postfix.EqualTo("FCL_FCTT"))isFCL_FCTT = true;
  else if(LT_postfix.EqualTo("FCL_FCHPT"))isFCL_FCHPT = true;
  else {
    cout<<"ERROR \t Unknown tight-loose postfix "<<LT_postfix.Data()<<endl;
    return false;
  }

  if(!isBLjet(idx))return false;
  if((!is2L2J && !isZMET) && fabs(jetEta[idx])>=2.4)return false;
  //<------------------------------------- COMMENT BACK IN. NOT AVAILABLE IN 2L2J OFFICIIAL NTUPLES
  if(!doRJR && !is2L2J && jetdl1r[idx]>=0.665)return true; // 85%
  else if(jetdl1r[idx]>=2.195)return true; // 77%
  //<------------------------------------- COMMENT BACK IN. NOT AVAILABLE IN 2L2J OFFICIIAL NTUPLES
  //if(!doRJR && !is2L2J&& jetMV2c10[idx]>0.11)return true; // 85%
  //else if(jetMV2c10[idx]>0.64)return true; // 77%
  return false;

}


bool MySusySkimAnalysis::isL(int idx, bool checkOR)
{


  bool isEl = false;
  bool isMu = false;

  if(lepFlavor[idx] == 1)isEl = true;
  else isMu = true;
  
  bool isZMET       = false;
  bool is2L2J       = false;
  bool isFCT        = false;
  bool isFCT_MT     = false;
  bool isODA        = false;
  bool isFCL        = false;
  bool isFCL_HPTC   = false;
  bool isGRAD      = false;
  bool isFCL_FCTT   = false;
  bool isFCL_FCHPT  = false;

  // for Helen
  bool isHP_MT_FCL_FCHPT = false;
  bool isHP_0M_FCL_FCHPT = false;
  bool isHP_MT_FCT = false;
  bool isHP_0M_FCT = false;
  if(LT_postfix.EqualTo("2L2J"))is2L2J = true;
  else if(LT_postfix.EqualTo("ZMET"))isZMET = true;
  else if(LT_postfix.EqualTo("HP_MT_FCL_FCHPT"))isHP_MT_FCL_FCHPT = true;
  else if(LT_postfix.EqualTo("HP_0M_FCL_FCHPT"))isHP_0M_FCL_FCHPT = true;
  else if(LT_postfix.EqualTo("HP_MT_FCT"))isHP_MT_FCT = true;
  else if(LT_postfix.EqualTo("HP_0M_FCT"))isHP_0M_FCT = true;
  else if(LT_postfix.EqualTo("FCT_MT"))isFCT_MT = true;
  else if(LT_postfix.EqualTo("FCT"))isFCT = true;
  else if(LT_postfix.EqualTo("ODA"))isODA = true;
  else if(LT_postfix.EqualTo("FCL") || LT_postfix.EqualTo("notrigm") || LT_postfix.EqualTo("med") || LT_postfix.EqualTo("loose") || LT_postfix.EqualTo("tight_iso") || LT_postfix.EqualTo("no"))isFCL = true;
  else if(LT_postfix.EqualTo("FCL_HPTC"))isFCL_HPTC = true;
  else if(LT_postfix.EqualTo("GRAD"))isGRAD = true;
  else if(LT_postfix.EqualTo("FCL_FCTT"))isFCL_FCTT = true;
  else if(LT_postfix.EqualTo("FCL_FCHPT"))isFCL_FCHPT = true;
  else {
    cout<<"ERROR \t Unknown tight-loose postfix "<<LT_postfix.Data()<<endl;
    return false;
  }

  if(isODA)return true;

  if(checkOR && !lepPassOR->at(idx))return false;

  // cout<<"lepFlavor["<<idx<<"] =" <<lepFlavor[idx]<<endl;
  // cout<<"lepPt["<<idx<<"] = "<<lepPt[idx]<<endl;
  // cout<<"lepEta["<<idx<<"] = "<<lepEta[idx]<<endl;
  // cout<<"lepPassBL->at("<<idx<<") = "<<lepPassBL->at(idx)<<endl;
  // cout<<"lepMedium->at("<<idx<<") = "<<lepMedium->at(idx)<<endl;

  //if(doRJR && isFCT && lepPt[idx] <= 20){/**cout<<"a"<<endl;*/return false;}
  //if(!doRJR && lepPt[idx] <= 25){/**cout<<"a"<<endl;*/return false;}
  if((fabs(lepZ0SinTheta[idx])>=0.5)){return false;}
  
  if(is2L2J){
    if(!lepPassBL->at(idx))return false;
    //if(lepPt[idx] < 20)return false;
    //if((isMu && fabs(lepD0Sig[idx]) >= 3))return false;
    if((isEl && fabs(lepEta[idx]) >= 2.47) || (isMu && fabs(lepEta[idx]) >= 2.7))return false;  //<------------------------------------- COMMENT BACK IN. NOT AVAILABLE IN 2L2J OFFICIIAL NTUPLES
    if((isEl && !(lepPassBL->at(idx) && lepVeryLoose->at(idx))) || (isMu && !lepMedium->at(idx)))return false;  //<------------------------------------- COMMENT BACK IN. NOT AVAILABLE IN 2L2J OFFICIIAL NTUPLES
  }else if(isZMET){
    if(isEl && (fabs(lepEta[idx]) >= 2.47 || (fabs(lepEta[idx]) >= 1.37 && fabs(lepEta[idx]) <= 1.52)))return false;
    if(isMu && fabs(lepEta[idx]) >= 2.5)return false;
    if((isMu && fabs(lepD0Sig[idx]) >= 3) || (isEl && fabs(lepD0Sig[idx]) >= 5))return false;
    if(isEl && !(lepPassBL->at(idx) && lepLoose->at(idx)))return false;
    if(isMu && !lepHighPt->at(idx))return false;
  }else if( doRJR && isFCT      && ((isEl && fabs(lepEta[idx]) >= 2.47) || (isMu && fabs(lepEta[idx]) >= 2.4))){/**cout<<"b"<<endl;*/return false;
  }else if(isFCL && ((isEl && fabs(lepEta[idx]) >= 2.47) || (isEl && !lepLoose->at(idx)) || (isEl && !lepPassBL->at(idx)) ||
                     (isMu && fabs(lepEta[idx]) >= 2.6 ) || (isMu && !lepMedium->at(idx)))){return false;
  }else if(isFCL_HPTC && ((isEl && fabs(lepEta[idx]) >= 2.47) || (isMu && fabs(lepEta[idx]) >= 2.6))){return false;
  }else if(!doRJR && ((isEl && fabs(lepEta[idx]) >= 2.47) || (isMu && fabs(lepEta[idx]) >= 2.7))){/**cout<<"b"<<endl;*/return false;}
  //if( isODA               && ((isEl && fabs(lepEta[idx]) >= 2.47) || (isMu && fabs(lepEta[idx]) >= 2.5))){/**cout<<"b"<<endl;*/return false;}
  
  // Not neded, all leptons passes these cuts
  //if(isFCL && (isEl && !(lepLoose->at(idx) && lepPassBL->at(idx)))){return false;}

  // if((isFCT_MT || isHP_MT_FCL_FCHPT || isHP_MT_FCT) && !lepMedium->at(idx)){
  //   return false;
  // }else if(isMu && !(isHP_0M_FCL_FCHPT || isHP_0M_FCT) && !lepMedium->at(idx)){return false;}
  
  return true;

}

bool MySusySkimAnalysis::isT(int idx, bool checkOR)
{

  //if(!lepSignal->at(idx))return false;


  bool isEl = false;
  bool isMu = false;

  if(lepFlavor[idx] == 1)isEl = true;
  else isMu = true;

  if(!isL(idx,checkOR))return false;
  if((isMu && fabs(lepD0Sig[idx]) >= 3) || (isEl && fabs(lepD0Sig[idx]) >= 5))return false;

  // Moved this to the definition of tight
  //if(!lepPassOR->at(idx))return false;
  bool isZMET       = false;
  bool is2L2J       = false;
  bool isFCT        = false;
  bool isFCT_MT     = false;
  bool isODA        = false;
  bool isFCL        = false;
  bool isFCL_HPTC   = false;
  bool isGRAD      = false;
  bool isFCL_FCTT   = false;
  bool isFCL_FCHPT  = false;
  // for Helen
  bool isHP_MT_FCL_FCHPT = false;
  bool isHP_0M_FCL_FCHPT = false;
  bool isHP_MT_FCT = false;
  bool isHP_0M_FCT = false;
  if(LT_postfix.EqualTo("2L2J"))is2L2J = true;
  else if(LT_postfix.EqualTo("ZMET"))isZMET = true;
  else if(LT_postfix.EqualTo("HP_MT_FCL_FCHPT"))isHP_MT_FCL_FCHPT = true;
  else if(LT_postfix.EqualTo("HP_0M_FCL_FCHPT"))isHP_0M_FCL_FCHPT = true;
  else if(LT_postfix.EqualTo("HP_MT_FCT"))isHP_MT_FCT = true;
  else if(LT_postfix.EqualTo("HP_0M_FCT"))isHP_0M_FCT = true;
  else if(LT_postfix.EqualTo("FCT_MT"))isFCT_MT = true;
  else if(LT_postfix.EqualTo("FCT"))isFCT = true;
  else if(LT_postfix.EqualTo("ODA"))isODA = true;
  else if(LT_postfix.EqualTo("FCL") || LT_postfix.EqualTo("notrigm") || LT_postfix.EqualTo("med") || LT_postfix.EqualTo("loose") || LT_postfix.EqualTo("tight_iso") || LT_postfix.EqualTo("no"))isFCL = true;
  else if(LT_postfix.EqualTo("FCL_HPTC"))isFCL_HPTC = true;
  else if(LT_postfix.EqualTo("GRAD"))isGRAD = true;
  else if(LT_postfix.EqualTo("FCL_FCTT"))isFCL_FCTT = true;
  else if(LT_postfix.EqualTo("FCL_FCHPT"))isFCL_FCHPT = true;
  else {
    cout<<"ERROR \t Unknown tight-loose postfix "<<LT_postfix.Data()<<endl;
    return false;
  }



  if(is2L2J){
    if((isMu && !lepIsoTightTrackOnly_VarRad->at(idx)) || (isEl && !lepIsoFCTight->at(idx)))return false;//<------------------------------------- COMMENT BACK IN. NOT AVAILABLE IN 2L0J OFFICIIAL NTUPLES
    //if((isMu && !lepIsoFCTight->at(idx)) || (isEl && !lepIsoFCTight->at(idx)))return false;
    if(!lepMedium->at(idx))return false;//<------------------------------------- COMMENT BACK IN. NOT AVAILABLE IN 2L0J OFFICIIAL NTUPLES
    //return false;
    //if(!lepSignal->at(idx))return false; //<------------------------------------- COMMENT BACK IN. NOT AVAILABLE IN 2L0J OFFICIIAL NTUPLES
  }else if(isZMET){
    if(isEl && !lepMedium->at(idx))return false;
    if(isEl && !lepIsoFCTight->at(idx))return false;
    if(isMu && !lepIsoTightTrackOnly_VarRad->at(idx))return false;
  }else if(isFCL){
    //return false;
    if((isMu && !lepMedium->at(idx)) || (isEl && !lepMedium->at(idx)))return false; //<------------------------------------- COMMENT BACK IN. NOT AVAILABLE IN 2L2J OFFICIIAL NTUPLES
    if((isMu && !lepIsoTightTrackOnly_VarRad->at(idx)) || (isEl && !lepIsoFCLoose->at(idx)))return false; //<------------------------------------- COMMENT BACK IN. NOT AVAILABLE IN 2L2J OFFICIIAL NTUPLES
  }
  
  return true;

}

void MySusySkimAnalysis::fillFakeRateHist(TString key, int idx, bool lepIsTight){

  if(!fullScan)return;
  // if(key.Contains("FAKE2L03")){ 
  //   cout<<"Filling histograms with "<<key.Data()<<endl;
  //   cout<<"n_clj = "<<n_clj<<endl;
  // }

  if(*met_Et < 50)h_lep_pT_nL[key+"_met_min_50"]->Fill(lepPt[idx],wgt);
  else if(*met_Et < 100)h_lep_pT_nL[key+"_met_50_100"]->Fill(lepPt[idx],wgt);
  else if(*met_Et < 150)h_lep_pT_nL[key+"_met_100_150"]->Fill(lepPt[idx],wgt);
  else if(*met_Et < 200)h_lep_pT_nL[key+"_met_150_200"]->Fill(lepPt[idx],wgt);
  else h_lep_pT_nL[key+"_met_200_max"]->Fill(lepPt[idx],wgt);
  /**
     if(metsig < 10)h_lep_pT_nL[key+"_metsig_min_10"]->Fill(lepPt[idx],wgt);
     else if(metsig < 20)h_lep_pT_nL[key+"_metsig_10_20"]->Fill(lepPt[idx],wgt);
     else if(metsig < 30)h_lep_pT_nL[key+"_metsig_20_30"]->Fill(lepPt[idx],wgt);
     else if(metsig < 50)h_lep_pT_nL[key+"_metsig_30_50"]->Fill(lepPt[idx],wgt);
     else h_lep_pT_nL[key+"_metsig_50_max"]->Fill(lepPt[idx],wgt);
  */
  if(themll < 50)h_lep_pT_nL[key+"_mll_min_50"]->Fill(lepPt[idx],wgt);
  else if(themll < 80)h_lep_pT_nL[key+"_mll_50_80"]->Fill(lepPt[idx],wgt);
  else if(themll < 110)h_lep_pT_nL[key+"_mll_80_110"]->Fill(lepPt[idx],wgt);
  else if(themll < 200)h_lep_pT_nL[key+"_mll_110_200"]->Fill(lepPt[idx],wgt);
  else if(themll < 300)h_lep_pT_nL[key+"_mll_200_300"]->Fill(lepPt[idx],wgt);
  else h_lep_pT_nL[key+"_mll_300_max"]->Fill(lepPt[idx],wgt);

  if(themt2 < 50)h_lep_pT_nL[key+"_mt2_min_50"]->Fill(lepPt[idx],wgt);
  else if(themt2 < 100)h_lep_pT_nL[key+"_mt2_50_100"]->Fill(lepPt[idx],wgt);
  else if(themt2 < 150)h_lep_pT_nL[key+"_mt2_100_150"]->Fill(lepPt[idx],wgt);
  else if(themt2 < 200)h_lep_pT_nL[key+"_mt2_150_200"]->Fill(lepPt[idx],wgt);
  else h_lep_pT_nL[key+"_mt2_200_max"]->Fill(lepPt[idx],wgt);

  if(n_clj < 1)h_lep_pT_nL[key+"_nclj_0_1"]->Fill(lepPt[idx],wgt);
  else if(n_clj < 2)h_lep_pT_nL[key+"_nclj_1_2"]->Fill(lepPt[idx],wgt);
  else h_lep_pT_nL[key+"_nclj_2_max"]->Fill(lepPt[idx],wgt);

  if(n_bjet < 1)h_lep_pT_nL[key+"_nbj_0_1"]->Fill(lepPt[idx],wgt);
  else h_lep_pT_nL[key+"_nbj_1_max"]->Fill(lepPt[idx],wgt);

  if(lepIsTight){
    if(*met_Et < 50)h_lep_pT_nT[key+"_met_min_50"]->Fill(lepPt[idx],wgt);
    else if(*met_Et < 100)h_lep_pT_nT[key+"_met_50_100"]->Fill(lepPt[idx],wgt);
    else if(*met_Et < 150)h_lep_pT_nT[key+"_met_100_150"]->Fill(lepPt[idx],wgt);
    else if(*met_Et < 200)h_lep_pT_nT[key+"_met_150_200"]->Fill(lepPt[idx],wgt);
    else h_lep_pT_nT[key+"_met_200_max"]->Fill(lepPt[idx],wgt);
    /**
    if(metsig < 10)h_lep_pT_nT[key+"_metsig_min_10"]->Fill(lepPt[idx],wgt);
    else if(metsig < 20)h_lep_pT_nT[key+"_metsig_10_20"]->Fill(lepPt[idx],wgt);
    else if(metsig < 30)h_lep_pT_nT[key+"_metsig_20_30"]->Fill(lepPt[idx],wgt);
    else if(metsig < 50)h_lep_pT_nT[key+"_metsig_30_50"]->Fill(lepPt[idx],wgt);
    else h_lep_pT_nT[key+"_metsig_50_max"]->Fill(lepPt[idx],wgt);
    */
    if(themll < 50)h_lep_pT_nT[key+"_mll_min_50"]->Fill(lepPt[idx],wgt);
    else if(themll < 80)h_lep_pT_nT[key+"_mll_50_80"]->Fill(lepPt[idx],wgt);
    else if(themll < 110)h_lep_pT_nT[key+"_mll_80_110"]->Fill(lepPt[idx],wgt);
    else if(themll < 200)h_lep_pT_nT[key+"_mll_110_200"]->Fill(lepPt[idx],wgt);
    else if(themll < 300)h_lep_pT_nT[key+"_mll_200_300"]->Fill(lepPt[idx],wgt);
    else h_lep_pT_nT[key+"_mll_300_max"]->Fill(lepPt[idx],wgt);

    if(themt2 < 50)h_lep_pT_nT[key+"_mt2_min_50"]->Fill(lepPt[idx],wgt);
    else if(themt2 < 100)h_lep_pT_nT[key+"_mt2_50_100"]->Fill(lepPt[idx],wgt);
    else if(themt2 < 150)h_lep_pT_nT[key+"_mt2_100_150"]->Fill(lepPt[idx],wgt);
    else if(themt2 < 200)h_lep_pT_nT[key+"_mt2_150_200"]->Fill(lepPt[idx],wgt);
    else h_lep_pT_nT[key+"_mt2_200_max"]->Fill(lepPt[idx],wgt);

    if(n_clj < 1)h_lep_pT_nT[key+"_nclj_0_1"]->Fill(lepPt[idx],wgt);
    else if(n_clj < 2)h_lep_pT_nT[key+"_nclj_1_2"]->Fill(lepPt[idx],wgt);
    else h_lep_pT_nT[key+"_nclj_2_max"]->Fill(lepPt[idx],wgt);

    if(n_bjet < 1)h_lep_pT_nT[key+"_nbj_0_1"]->Fill(lepPt[idx],wgt);
    else h_lep_pT_nT[key+"_nbj_1_max"]->Fill(lepPt[idx],wgt);
  }
}

// truthENUM MySusySkimAnalysis::TruthClassification(int T, int O, int firstEgMotherT, int firstEgMotherO, int firstEgMotherPdgId, int lepCharge, int lepFlav)
// {


//   bool C1 = (T==2 || ( T==4  && O==5 && fabs(firstEgMotherPdgId) == 11));

//   if(T==6  && (O==2  || O==10  || O==12  || O==13  || O==14  || O==15  || O==22  || O==43))return REAL;
//   if((C1 && (firstEgMotherPdgId*lepCharge)<0) || (T==4 && O==5 && firstEgMotherO==40) || (T==15 && O==40))return REAL;
//   if(C1 && (firstEgMotherPdgId*lepCharge)>0)return CF;
//   if((T==14  && O==37) ||
//      (T==4 && O==5 && firstEgMotherT==14 && firstEgMotherO==37) ||
//      (T==4 && O==5 && firstEgMotherT==16  && firstEgMotherO==38) ||
//      (T==4 && O==5 && firstEgMotherT==15  && firstEgMotherO==39) ||
//      (T==4 && O==5 && firstEgMotherT==13 && fabs(firstEgMotherPdgId) == 22) ||
//      (T==15 && O == 39))return CO;
//   if((T==4  && (O==6  || O==7  || O==23  || O==24  )) ||
//      (T==4  && O==5 && firstEgMotherT==16  && ( firstEgMotherO==42  || firstEgMotherO==23  || firstEgMotherO==24  )) ||
//      (T==16 && (O==42  || O==23  )) ||
//      (T==8  && ( O==34  || O==35  || O==23  || O==24  )) ||
//      (T==17) )return LF;
//   if((( T==3 || T==15 ) && O==9)  ||
//      ( T==4 && O==5 && firstEgMotherT==15 && firstEgMotherO==9 ) ||
//      (T==7  && O==9)  ||
//      (((T==3  || T==15 ) && O==8 )) ||
//      ( T==4 && O==5 && firstEgMotherT==15 && firstEgMotherO==8))return HF;
//   if((T==3  && (O==26  || O==29  || O==33  )) ||
//      (T==16  && O==26)  ||
//      ( T==4 && O==5 && firstEgMotherT==16 && firstEgMotherO==26 ) ||
//      (( T==6  || T==7 ) && (O==26  || O==29  || O==33 )))return HF;
//   if((T==3  && ( O==25  || O==32  || O==27  )) ||
//      (T==4  && O==27)  ||
//      (T==16  && ( O==25  || O==27  )) ||
//      ( T==4 && O==5 && firstEgMotherT==16 && ( firstEgMotherO==25 || firstEgMotherO==27 )) ||
//      (T==7  && ( O==25  || O==32  || O==27  )) ||
//      (( T==6  || T==8  ) && ( O==27 /** || originbkgLep==27*/ )))return HF;

//   //cout<<"UNKNOWN "<<lepFlav<<endl;
//   //cout<<"T = "<<T<<", O = "<<O<<", firstEgMotherT = "<<firstEgMotherT<<", firstEgMotherO = "<<firstEgMotherO<<", firstEgMotherPdgId = "<<firstEgMotherPdgId <<", lepCharge = "<<lepCharge<<endl;

//   return UNKNOWN;
// }

void MySusySkimAnalysis::checkTrigMatch(TString trigname, int idx, TString key, bool lepIsTight, TString trigtype){
  int verb = 0;
  int trigidx = -1;
  bool lepIsTriggered = false;
  bool evIsTriggered  = false;

  TObjArray *tx = trigtype.Tokenize("_");
  // If trigtype contains e_<cat> or m_<cat>: need to extract category only
  if(tx->GetSize() > 1){
    trigtype = (((TObjString *)tx->Last())->String());
  }
  delete tx;
  
  std::vector<TString>::iterator it  = std::find(triglist.begin(), triglist.end(), trigname); 
  if(it != triglist.end()){
    trigidx = it - triglist.begin();
    if(verb)std::cout << "DEBUG \t 1) Element "<<trigname<<" found in index "<<trigidx<<endl;
  }else{
    std::vector<TString>::iterator it  = std::find(categories.begin(), categories.end(), trigname);
    if(it != categories.end()){
      trigidx = it - categories.begin();
      if(verb)std::cout << "DEBUG \t 2) Element "<<trigname<<" found in index "<<trigidx<<endl;
    }else{
      std::cout << "ERROR \t Element "<<trigname<<" not found."; 
      return;
    }
  }
  for(unsigned int i=0; i<trigvallep.size();i++){
    if(triglist.at(i).Contains("Trig"))continue;
    if(!trigcategories[triglist.at(i)].Contains(trigtype) && !trigcat[triglist.at(i)].Contains(trigtype))continue;
    if((lepFlavor[idx] == 1 && !(std::count(trigstr_e[yr].begin(), trigstr_e[yr].end(), triglist.at(i)))) ||
       (lepFlavor[idx] == 2 && !(std::count(trigstr_m[yr].begin(), trigstr_m[yr].end(), triglist.at(i)))))continue;
    if(verb)std::cout << "DEBUG \t Trigger "<<trigname<<" matching against "<<triglist.at(i)<<endl; 
    if(trigvallep.at(i)->at(idx) && lepPt[idx] >= getPtThresholdTrigger(triglist.at(i),lepFlavor[idx] == 1 ? "e" : "mu"))lepIsTriggered = true;
    if(*trigval.at(i))evIsTriggered = true;  
  }
  //if(!lepIsTriggered){

  

  if(h_lep_pT_nL.find(trigname+"_lepmatched_"+key) == h_lep_pT_nL.end()){
    cout<<"WARNING \t Could not find key "<<trigname+"_lepmatched_"+key<<endl;
    return;
  }
  
  if(!evIsTriggered){
    if(trigname.Contains("Trig")){
      if(verb)std::cout<<"Event not triggered : key = "<<trigname+"_"+key<<std::endl;
      for (TString bdt_name : bdt_vec) {
	h_lep_BDT_nL[bdt_name+"_"+trigname+"_evnotrig_"+key]->Fill(BDTweight[bdt_name],wgt);
      }
      h_lep_pT_nL[trigname+"_evnotrig_"+key]->Fill(lepPt[idx],wgt);
      h_lep_pT_eta_nL[trigname+"_evnotrig_"+key]->Fill(lepPt[idx],lepEta[idx],wgt);
      //h_lep_trig_nL["lepnotrig_"+key]->AddBinContent(trigidx+1,wgt);
      if(lepIsTight){
	//h_lep_trig_nT["lepnotrig_"+key]->AddBinContent(trigidx+1,wgt);
       
	for (TString bdt_name : bdt_vec) {
	  h_lep_BDT_nT[bdt_name+"_"+trigname+"_evnotrig_"+key]->Fill(BDTweight[bdt_name],wgt);
	}
	h_lep_pT_nT[trigname+"_evnotrig_"+key]->Fill(lepPt[idx],wgt);
	h_lep_pT_eta_nT[trigname+"_evnotrig_"+key]->Fill(lepPt[idx],lepEta[idx],wgt);
      }
    }
  }else{
    if(!lepIsTriggered){
       if(verb)std::cout<<"Lep not triggered : key = "<<trigname+"_"+key<<std::endl;
      for (TString bdt_name : bdt_vec) {
	h_lep_BDT_nL[bdt_name+"_"+trigname+"_lepnotmatched_"+key]->Fill(BDTweight[bdt_name],wgt);
      }
      h_lep_pT_nL[trigname+"_lepnotmatched_"+key]->Fill(lepPt[idx],wgt);
      h_lep_pT_eta_nL[trigname+"_lepnotmatched_"+key]->Fill(lepPt[idx],lepEta[idx],wgt);
      if(trigidx>=0)h_lep_trig_nL["lepnotmatched_"+key]->AddBinContent(trigidx+1,wgt);
      if(lepIsTight){
	 
	for (TString bdt_name : bdt_vec) {
	  h_lep_BDT_nT[bdt_name+"_"+trigname+"_lepnotmatched_"+key]->Fill(BDTweight[bdt_name],wgt);
	}
	if(trigidx>=0)h_lep_trig_nT["lepnotmatched_"+key]->AddBinContent(trigidx+1,wgt);
	h_lep_pT_nT[trigname+"_lepnotmatched_"+key]->Fill(lepPt[idx],wgt);
	h_lep_pT_eta_nT[trigname+"_lepnotmatched_"+key]->Fill(lepPt[idx],lepEta[idx],wgt);
      }
    }else{
      if(verb)std::cout<<"Lep triggered : key = "<<trigname+"_"+key<<std::endl;
      for (TString bdt_name : bdt_vec) {
	cout<<bdt_name+"_"+trigname+"_"+key<<endl;
	h_lep_BDT_nL[bdt_name+"_"+trigname+"_lepmatched_"+key]->Fill(BDTweight[bdt_name],wgt);
      }
      h_lep_pT_nL[trigname+"_lepmatched_"+key]->Fill(lepPt[idx],wgt);
      h_lep_pT_eta_nL[trigname+"_lepmatched_"+key]->Fill(lepPt[idx],lepEta[idx],wgt);
      if(trigidx>=0)h_lep_trig_nL["lepmatched_"+key]->AddBinContent(trigidx+1,wgt);
      if(lepIsTight){
	for (TString bdt_name : bdt_vec) {
	  h_lep_BDT_nT[bdt_name+"_"+trigname+"_lepmatched_"+key]->Fill(BDTweight[bdt_name],wgt);
	}
	if(trigidx>=0)h_lep_trig_nT["lepmatched_"+key]->AddBinContent(trigidx+1,wgt);
	h_lep_pT_nT[trigname+"_lepmatched_"+key]->Fill(lepPt[idx],wgt);
	h_lep_pT_eta_nT[trigname+"_lepmatched_"+key]->Fill(lepPt[idx],lepEta[idx],wgt);
      }
    }
  }
  /**else if(trigname.EqualTo("2LTrig") && !lepIsTriggered && key.Contains("EEOS")){
     cout<<"WARNING \t Event is triggered by 2LTrig, but lepton is not matched to any trigger: pT = " <<lepPt[idx]<<endl;
     for(unsigned int i = 0; i<lepPt.GetSize(); i++){
     printf("Lepton %i: \n",i);
     printf("pT %.2f: \n",lepPt[i]);
     printf("2e24 %i: \n",(int)lepHLT_2e24_lhvloose_nod0->at(i));
     printf("2e17 %i: \n",(int)lepHLT_2e17_lhvloose_nod0_L12EM15VHI->at(i));
     printf("e26,mu8 %i: \n",(int)lepHLT_e26_lhmedium_nod0_mu8noL1->at(i));
     printf("e7,mu24 %i: \n",(int)lepHLT_e7_lhmedium_nod0_mu24->at(i));
     printf("e17,mu14 %i: \n",(int)lepHLT_e17_lhloose_nod0_mu14->at(i));
     }
     }*/
}
void MySusySkimAnalysis::DefineCutFlow(){

  if(cutflowstr.EqualTo("CRTOP")){
    cutflow["== 2 baseline"] = 0;
    cutflow_wgt["== 2 baseline"] = 0;
    cutfloworder.push_back("== 2 baseline");
    cutflow["DF"] = 0;
    cutflow_wgt["DF"] = 0;
    cutfloworder.push_back("DF");
    cutflow["OS"] = 0;
    cutflow_wgt["OS"] = 0;
    cutfloworder.push_back("OS");
    cutflow["== 2 signal"] = 0;
    cutflow_wgt["== 2 signal"] = 0;
    cutfloworder.push_back("== 2 signal");
    cutflow["pT > 25"] = 0;
    cutflow_wgt["pT > 25"] = 0;
    cutfloworder.push_back("pT > 25");
    cutflow["n_clj == 0"] = 0;
    cutflow_wgt["n_clj == 0"] = 0;
    cutfloworder.push_back("n_clj == 0");
    cutflow["n_bjet == 1"] = 0;
    cutflow_wgt["n_bjet == 1"] = 0;
    cutfloworder.push_back("n_bjet == 1");
    cutflow["mll > 11"] = 0;
    cutflow_wgt["mll > 11"] = 0;
    cutfloworder.push_back("mll > 11");
    //cutflow["met_Et > 110"] = 0;
    //cutflow_wgt["met_Et > 110"] = 0;
    //cutfloworder.push_back("met_Et > 110");
    cutflow["met_Sign > 8"] = 0;
    cutflow_wgt["met_Sign > 8"] = 0;
    cutfloworder.push_back("met_Sign > 8");
    //cutflow["themt2 > 80"] = 0;
    //cutflow_wgt["themt2 > 80"] = 0;
    //cutfloworder.push_back("themt2 > 80");
  }else if(cutflowstr.EqualTo("CRWW")){
    cutflow["== 2 baseline"] = 0;
    cutflow_wgt["== 2 baseline"] = 0;
    cutfloworder.push_back("== 2 baseline");
    cutflow["DF"] = 0;
    cutflow_wgt["DF"] = 0;
    cutfloworder.push_back("DF");
    cutflow["OS"] = 0;
    cutflow_wgt["OS"] = 0;
    cutfloworder.push_back("OS");
    cutflow["== 2 signal"] = 0;
    cutflow_wgt["== 2 signal"] = 0;
    cutfloworder.push_back("== 2 signal");
    cutflow["pT > 25"] = 0;
    cutflow_wgt["pT > 25"] = 0;
    cutfloworder.push_back("pT > 25");
    cutflow["n_clj == 0"] = 0;
    cutflow_wgt["n_clj == 0"] = 0;
    cutfloworder.push_back("n_clj == 0");
    cutflow["n_bjet == 0"] = 0;
    cutflow_wgt["n_bjet == 0"] = 0;
    cutfloworder.push_back("n_bjet == 0");
    cutflow["mll > 100"] = 0;
    cutflow_wgt["mll > 100"] = 0;
    cutfloworder.push_back("mll > 100");
    cutflow["60 <= met_Et <= 100"] = 0;
    cutflow_wgt["60 <= met_Et <= 100"] = 0;
    cutfloworder.push_back("60 <= met_Et <= 100");
    cutflow["5 <= met_Sign <= 10"] = 0;
    cutflow_wgt["5 <= met_Sign <= 10"] = 0;
    cutfloworder.push_back("5 <= met_Sign <= 10");
    cutflow["60 <= themt2 <= 65"] = 0;
    cutflow_wgt["60 <= themt2 <= 65"] = 0;
    cutfloworder.push_back("60 <= themt2 <= 65");
  }
    
}

void MySusySkimAnalysis::printCutflow(){
  printf("CUTFLOW FOR %s:\n",cutflowstr.Data());
  printf("%-20s %10s %10s\n","CUT","COUNT","WGT");
 
  for(unsigned int i= 0; i<cutfloworder.size(); i++){
    printf("%-20s %10.0f %10.2f\n",cutfloworder.at(i).Data(),cutflow[cutfloworder.at(i)],cutflow_wgt[cutfloworder.at(i)]);
  }

  // for(std::map<TString,float>::iterator iter = cutflow.begin(); iter != cutflow.end(); ++iter)
  //   {
  //     printf("%-20s %10.2f\n",iter->first.Data(),iter->second);
  //   }  
}

void MySusySkimAnalysis::printInfo(){
  printf("RunNumber = %i, EventNumber = %llu\n",*RunNumber,*EventNumber);
  for(unsigned int i = 0; i<lepPt.GetSize(); i++){
    printf("-----------------------------------------------------------------");
    printf("\n");
    printf("Lep is %s \n",(lepFlavor[i] == 1 ? "el": "mu"));
    printf("pT = %.2f, eta = %.2f\n",lepPt[i],lepEta[i]);
    printf("passOR = %s\n",lepPassOR->at(i) ? "yes" : "no");
    printf("z0 = %.2f, d0 = %.2f\n",fabs(lepZ0SinTheta[i]),fabs(lepD0Sig[i]));
    printf("BL = %s, loose = %s, medium = %s, tight = %s\n",lepPassBL->at(i) ? "yes" : "no",lepVeryLoose->at(i) ? "yes" : "no",lepMedium->at(i) ? "yes" : "no",lepTight->at(i) ? "yes" : "no"); //<---- COMMENT BACK IN. NOT AVAILABLE IN 2L2J OFFICIAL NTUPLES
    printf("iso = %s\n",lepIsoFCLoose->at(i) ? "yes" : "no"); //<---- COMMENT BACK IN. NOT AVAILABLE IN 2L2J OFFICIAL NTUPLES
  }
}
void MySusySkimAnalysis::WriteToFile(){

  cout<<nDup<<" duplicated events"<<endl;

  //myfile.close();
  //char ofile[200];
  //sprintf(ofile,"%s/%s_HISTOGRAMS.root",date_folder,fname.Data());
  //printf("Writing to %s/%s_HISTOGRAMS.root\n",date_folder,fname.Data());

  if(makeNtup){
    ntupf->Write();
    ntupf->Close();
    // //delete tree3;
    delete ntupf;
  }
  
  // In terminate;
  /**cout<<"a"<<endl;
  newtree->GetCurrentFile()->Write();
  cout<<"n"<<endl;
  delete newtree->GetCurrentFile(); //
  cout<<"c"<<endl;
  */
  TFile *f = new TFile(fnameoutputhistfilename,"RECREATE");
  TIterator *nextObject = (TIterator*)fOutput->MakeIterator();
  TObject *tobj;
  int i = 0;

  while((tobj= (TObject*)nextObject->Next())) {

    TString name = tobj->GetName();

    //cout<<"Writing "<<name.Data()<<endl;

    tobj->Write();

    i += 1;

  }
  printf("%i histograms written to file %s\n",i,fnameoutputhistfilename);
  f->Close();
  delete f;
  delete tobj;
  delete nextObject;
}


Float_t MySusySkimAnalysis::getPtThresholdTrigger(TString tname, TString tobj, int vb){

  //return 10.0;
  // Check if thresholds already found, use them directly from dictionary
  
  if(tobj.EqualTo("e")){
    if(!(el_trg_th.find(tname) == el_trg_th.end())){
      //cout<<"el_trg_th["<<tname<<"] = "<<el_trg_th[tname]<<endl;
      return el_trg_th[tname];
    }
  }else if(tobj.EqualTo("mu")){
    if(!(mu_trg_th.find(tname) == mu_trg_th.end())){
      //cout<<"mu_trg_th["<<tname<<"] = "<<mu_trg_th[tname]<<endl;
      return mu_trg_th[tname];
    }
  }
  
  if(vb)cout<<"Checking "<<tname.Data()<<endl;
  TObjArray *tx = tname.Tokenize("_");
  Float_t trig_ptcut_val = -1;
  
  for(int i = 0; i<tx->GetEntries(); i++){
    TString trig_ptcut = ((TObjString *)tx->At(i))->String();
    if(vb)cout<<"Checking value "<<i<<" : "<<trig_ptcut.Data()<<endl;
    if(vb)cout<<tobj.Data()<<" trigger object is at idx "<<trig_ptcut.First(tobj.Data())<<endl;
    //if(vb)cout<<"Muon trigger object is at idx "<<trig_ptcut.First("mu")<<endl;
    if(trig_ptcut.Contains(tobj))trig_ptcut.Remove(trig_ptcut.First(tobj.Data())-trig_ptcut.First(tobj.Data()),tobj.Length()+trig_ptcut.First(tobj.Data()));
    else continue;
    
    if(trig_ptcut.Length() > 0 && !trig_ptcut.IsDigit()){
      if(vb)cout<<"Trig branch "<<trig_ptcut.Data()<<" does not contain only numbers"<<endl;
      vector<Float_t> valid_pt_th;
      TString trig_ptcut_copy = trig_ptcut;
      for(int k = 0; k<trig_ptcut.Length(); k++){
	for(int j = 1; j<trig_ptcut.Length(); j++){
	  //TString cutted_str = trig_ptcut_copy.Strip(k,j)
	  //cout<<"Trying "<<trig_ptcut_copy.Remove(k,j).Data()<<endl;
	  if(trig_ptcut_copy.Remove(k,j).IsDigit()){
	    // Need an extra check to avoid concatinating 8 and 1 in strings like mu8noL1
	    if(trig_ptcut.Contains(trig_ptcut_copy)){
	      valid_pt_th.push_back(trig_ptcut_copy.Atoi());
	    }
	  }
	  trig_ptcut_copy = trig_ptcut;	  
	}
      }
      if(!valid_pt_th.size())continue;
      trig_ptcut_val = *max_element(valid_pt_th.begin(), valid_pt_th.end()); // C++11
    }else{
      if(trig_ptcut.Atoi() > trig_ptcut_val){
	trig_ptcut_val = trig_ptcut.Atoi();	
      }
    }
    if(vb)cout<<"pT cut for trigger "<< tname.Data() <<" is "<< trig_ptcut_val << endl;
    
    //else cout<<"ERROR \t Could not find any object matching "<<tobj.Data()<<" in trigger " << tname << endl;
    //else if(trig_ptcut.Contains("mu"))trig_ptcut.Remove(trig_ptcut.First("mu")-trig_ptcut.First("mu"),2+trig_ptcut.First("mu"));
    
  }
  // Update dictinoaries
  if(vb)cout<<"Updating maps for "<<tname.Data()<<" with value "<<trig_ptcut_val<<endl;
  if(tobj.EqualTo("e")){
    el_trg_th[tname] = trig_ptcut_val;
  }else if(tobj.EqualTo("mu")){
    mu_trg_th[tname] = trig_ptcut_val; 
  }
  if(vb)cout<<"Deleteeing tx"<<endl;
  delete tx;
  if(vb)cout<<"Returning value "<<trig_ptcut_val<<endl;
  return trig_ptcut_val;
  
}

void MySusySkimAnalysis::fillMYTree(int idx1, int idx2, std::map< TString, std::vector<double> > fake_wgt_vec){
  /*
  for(std::vector< TString >::iterator itLT2 = looseTightDef.begin(); itLT2 != looseTightDef.end(); itLT2++){
    //if(fake_wgt_vec[*itLT2].size() > 0)cout<<"Filling weight for "<<(*itLT2).Data()<<" with weight "<<fake_wgt_vec[*itLT2].at(0)<<endl;
    //else cout<< "No weights for "<< (*itLT2) << endl;
    MY->bMY_MM_weight.push_back(fake_wgt_vec[*itLT2].size() > 0 ? fake_wgt_vec[*itLT2].at(0) : 1.0);
    MY->bMY_syst_MM_up.push_back(fake_wgt_vec[*itLT2].size() > 1 ? fake_wgt_vec[*itLT2].at(1): 1.0);
    MY->bMY_syst_MM_down.push_back(fake_wgt_vec[*itLT2].size() > 2 ? fake_wgt_vec[*itLT2].at(2) : 1.0);
    MY->bMY_MM_key.push_back(*itLT2);
  }

  
  MY->bMY_trigMatch_1L2LTrig =                                     (*trigMatch_1L2LTrig);                                
  MY->bMY_trigMatch_1LTrig =                                   	   (*trigMatch_1LTrig);                                  
  MY->bMY_trigMatch_1L2LTrigOR =                               	   (*trigMatch_1L2LTrigOR);                              
  MY->bMY_trigMatch_1LTrigOR =                                 	   (*trigMatch_1LTrigOR);                                
  MY->bMY_trigMatch_2LTrig =                                   	   (*trigMatch_2LTrigOR);                                  
  MY->bMY_trigMatch_2LTrigOR =                                 	   (*trigMatch_2LTrigOR);                                
  MY->bMY_trigMatch_HLT_e24_lhmedium_L1EM20VH =                	   (*trigMatch_HLT_e24_lhmedium_L1EM20VH);               
  MY->bMY_trigMatch_HLT_e60_lhmedium =                         	   (*trigMatch_HLT_e60_lhmedium);                        
  MY->bMY_trigMatch_HLT_e120_lhloose =                         	   (*trigMatch_HLT_e120_lhloose);                        
  MY->bMY_trigMatch_HLT_mu20_iloose_L1MU15 =                   	   (*trigMatch_HLT_mu20_iloose_L1MU15);                  
  MY->bMY_trigMatch_HLT_2e12_lhloose_L12EM10VH =               	   (*trigMatch_HLT_2e12_lhloose_L12EM10VH);              
  MY->bMY_trigMatch_HLT_mu18_mu8noL1 =                         	   (*trigMatch_HLT_mu18_mu8noL1);                        
  MY->bMY_trigMatch_HLT_e17_lhloose_mu14 =                     	   (*trigMatch_HLT_e17_lhloose_mu14);                    
  MY->bMY_trigMatch_HLT_e7_lhmedium_mu24 =                     	   (*trigMatch_HLT_e7_lhmedium_mu24);                    
  MY->bMY_trigMatch_HLT_e24_lhtight_nod0_ivarloose =           	   (*trigMatch_HLT_e24_lhtight_nod0_ivarloose);          
  MY->bMY_trigMatch_HLT_e24_lhmedium_nod0_L1EM20VH =           	   (*trigMatch_HLT_e24_lhmedium_nod0_L1EM20VH);          
  MY->bMY_trigMatch_HLT_e60_medium =                           	   (*trigMatch_HLT_e60_medium);                          
  MY->bMY_trigMatch_HLT_mu40 =                                 	   (*trigMatch_HLT_mu40);                                
  MY->bMY_trigMatch_HLT_mu24_iloose_L1MU15 =                   	   (*trigMatch_HLT_mu24_iloose_L1MU15);                  
  MY->bMY_trigMatch_HLT_mu24_ivarloose_L1MU15 =                	   (*trigMatch_HLT_mu24_ivarloose_L1MU15);               
  MY->bMY_trigMatch_HLT_mu24_ivarmedium =                      	   (*trigMatch_HLT_mu24_ivarmedium);                     
  MY->bMY_trigMatch_HLT_mu24_imedium =                         	   (*trigMatch_HLT_mu24_imedium);                        
  MY->bMY_trigMatch_HLT_mu26_imedium =                         	   (*trigMatch_HLT_mu26_imedium);                        
  MY->bMY_trigMatch_HLT_2e15_lhvloose_nod0_L12EM13VH =         	   (*trigMatch_HLT_2e15_lhvloose_nod0_L12EM13VH);        
  MY->bMY_trigMatch_HLT_2e17_lhvloose_nod0 =                   	   (*trigMatch_HLT_2e17_lhvloose_nod0);                  
  MY->bMY_trigMatch_HLT_2mu10 =                                	   (*trigMatch_HLT_2mu10);                               
  MY->bMY_trigMatch_HLT_2mu14 =                                	   (*trigMatch_HLT_2mu14);                               
  MY->bMY_trigMatch_HLT_mu20_mu8noL1 =                         	   (*trigMatch_HLT_mu20_mu8noL1);                        
  MY->bMY_trigMatch_HLT_mu22_mu8noL1 =                         	   (*trigMatch_HLT_mu22_mu8noL1);                        
  MY->bMY_trigMatch_HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1 =  	   (*trigMatch_HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1); 
  MY->bMY_trigMatch_HLT_e26_lhtight_nod0_ivarloose =           	   (*trigMatch_HLT_e26_lhtight_nod0_ivarloose);          
  MY->bMY_trigMatch_HLT_e26_lhtight_nod0 =                     	   (*trigMatch_HLT_e26_lhtight_nod0);                    
  MY->bMY_trigMatch_HLT_e60_lhmedium_nod0 =                    	   (*trigMatch_HLT_e60_lhmedium_nod0);                   
  MY->bMY_trigMatch_HLT_e140_lhloose_nod0 =                    	   (*trigMatch_HLT_e140_lhloose_nod0);                   
  MY->bMY_trigMatch_HLT_e300_etcut =                           	   (*trigMatch_HLT_e300_etcut);                          
  MY->bMY_trigMatch_HLT_mu26_ivarmedium =                      	   (*trigMatch_HLT_mu26_ivarmedium);                     
  MY->bMY_trigMatch_HLT_mu50 =                                 	   (*trigMatch_HLT_mu50);                                
  MY->bMY_trigMatch_HLT_mu60_0eta105_msonly =                  	   (*trigMatch_HLT_mu60_0eta105_msonly);                 
  MY->bMY_trigMatch_HLT_2e17_lhvloose_nod0_L12EM15VHI =        	   (*trigMatch_HLT_2e17_lhvloose_nod0_L12EM15VHI);       
  MY->bMY_trigMatch_HLT_2e24_lhvloose_nod0 =                   	   (*trigMatch_HLT_2e24_lhvloose_nod0);                  
  MY->bMY_trigMatch_HLT_e17_lhloose_nod0_mu14 =                	   (*trigMatch_HLT_e17_lhloose_nod0_mu14);               
  MY->bMY_trigMatch_HLT_e26_lhmedium_nod0_mu8noL1 =            	   (*trigMatch_HLT_e26_lhmedium_nod0_mu8noL1);           
  MY->bMY_trigMatch_HLT_e7_lhmedium_nod0_mu24 =                	   (*trigMatch_HLT_e7_lhmedium_nod0_mu24);               
  MY->bMY_mu =                                                 	   (*mu);                                                
  MY->bMY_avg_mu =                                             	   (*avg_mu);                                            
  MY->bMY_actual_mu =                                          	   (*actual_mu);                                         
  MY->bMY_nVtx =                                               	   (*nVtx);                                              
  MY->bMY_channel =                                            	   (*channel);                                           
  MY->bMY_nLep_base =                                          	   (*nLep_base);                                         
  MY->bMY_nLep_signal =                                        	   (*nLep_signal);   
  for(unsigned int i = 0; i<lepFlavor.GetSize(); i++){                          
    MY->bMY_lepFlavor.push_back                            	   (lepFlavor[i]);                                         
    MY->bMY_lepCharge.push_back                            	   (lepCharge[i]);                                         
    MY->bMY_lepAuthor.push_back                            	   (lepAuthor[i]);                                         
    MY->bMY_lepPt.push_back                                	   (lepPt[i]);                                             
    MY->bMY_lepEta.push_back                               	   (lepEta[i]);                                            
    MY->bMY_lepPhi.push_back                               	   (lepPhi[i]);                                            
    MY->bMY_lepM.push_back                                 	   (lepM[i]);                                              
    MY->bMY_lepD0.push_back                                	   (lepD0[i]);                                             
    MY->bMY_lepD0Sig.push_back                             	   (lepD0Sig[i]);                                          
    MY->bMY_lepZ0.push_back                                	   (lepZ0[i]);                                             
    MY->bMY_lepZ0SinTheta.push_back                        	   (lepZ0SinTheta[i]);                                     
    MY->bMY_lepPassOR.push_back                            	   (lepPassOR->at(i));                                         
    MY->bMY_lepType.push_back                              	   (lepType[i]);                                           
    MY->bMY_lepOrigin.push_back                            	   (lepOrigin[i]);                                         
    MY->bMY_lepEgMotherType.push_back                      	   (lepEgMotherType[i]);                                   
    MY->bMY_lepEgMotherOrigin.push_back                    	   (lepEgMotherOrigin[i]);                                 
    MY->bMY_lepEgMotherPdgId.push_back                     	   (lepEgMotherPdgId[i]);                                  
    MY->bMY_lepECIDS.push_back                             	   (lepECIDS[i]);                                          
    MY->bMY_lepPassBL.push_back                            	   (lepPassBL->at(i));                                         
    MY->bMY_lepVeryLoose.push_back                         	   (lepVeryLoose->at(i));                                      
    MY->bMY_lepLoose.push_back                             	   (lepLoose->at(i));                                          
    MY->bMY_lepMedium.push_back                            	   (lepMedium->at(i));                                         
    MY->bMY_lepTight.push_back                             	   (lepTight->at(i));                                          
    MY->bMY_lepIsoFCHighPtCaloOnly.push_back               	   (lepIsoFCHighPtCaloOnly->at(i));                            
    MY->bMY_lepIsoFixedCutHighPtTrackOnly.push_back        	   (lepIsoFixedCutHighPtTrackOnly->at(i));                     
    MY->bMY_lepIsoGradient.push_back                       	   (lepIsoGradient->at(i));                                    
    MY->bMY_lepIsoFCLoose.push_back                        	   (lepIsoFCLoose->at(i));                                     
    MY->bMY_lepIsoFCTight.push_back                        	   (lepIsoFCTight->at(i));                                     
    MY->bMY_lepIsoFCTightTrackOnly.push_back               	   (lepIsoFCTightTrackOnly->at(i));                            
    MY->bMY_lepTruthMatched.push_back                      	   (lepTruthMatched->at(i));                                   
    MY->bMY_lepTruthCharge.push_back                       	   (lepTruthCharge[i]);                                    
    MY->bMY_lepTruthPt.push_back                           	   (lepTruthPt[i]);                                        
    MY->bMY_lepTruthEta.push_back                          	   (lepTruthEta[i]);                                       
    MY->bMY_lepTruthPhi.push_back                          	   (lepTruthPhi[i]);                                       
    MY->bMY_lepTruthM.push_back                            	   (lepTruthM[i]);   
  }                                      
  MY->bMY_nJet30 =                                             	   (*nJet30);                                            
  MY->bMY_nJet20 =                                             	   (*nJet20);                                            
  MY->bMY_nBJet20_MV2c10_FixedCutBEff_77 =                     	   (*nBJet20_MV2c10_FixedCutBEff_77);  
        
  for(unsigned int i = 0; i<jetPt.GetSize(); i++){     
    MY->bMY_jetPt.push_back                                	   (jetPt[i]);                                             
    MY->bMY_jetEta.push_back                               	   (jetEta[i]);                                            
    MY->bMY_jetPhi.push_back                               	   (jetPhi[i]);                                            
    MY->bMY_jetM.push_back                                 	   (jetM[i]);                                              
    MY->bMY_jetJVT.push_back                               	   (jetJVT[i]);                                            
    MY->bMY_jetPassOR.push_back                            	   (jetPassOR->at(i));                                         
    MY->bMY_jetSignal.push_back                            	   (jetSignal->at(i));   
    MY->bMY_jetTileEnergy.push_back                        	   (jetTileEnergy[i]);                                     
    MY->bMY_jetMV2c10.push_back                            	   (jetMV2c10[i]); 
  }  

  MY->bMY_mjj =                                                	   (*mjj);                                                                                     
  MY->bMY_met_Et =                                             	   (*met_Et);                                            
  MY->bMY_met_Sign =                                           	   (*met_Sign);                                          
  MY->bMY_met_Phi =                                            	   (*met_Phi);                                           
  MY->bMY_mll =                                                	   (*mll);                                               
  MY->bMY_pileupWeight =                                       	   (*pileupWeight);                                      
  MY->bMY_leptonWeight =                                       	   (*leptonWeight);                                      
  MY->bMY_eventWeight =                                        	   (*eventWeight);                                       
  MY->bMY_bTagWeight =                                         	   (*bTagWeight);                                        
  MY->bMY_jvtWeight =                                          	   (*jvtWeight);                                         
  MY->bMY_globalDiLepTrigSF =                                  	   (*globalDiLepTrigSF);                                 
  MY->bMY_flavSymWeight =                                      	   (*flavSymWeight);                                     
  MY->bMY_genWeight =                                          	   (*genWeight);                                         
  MY->bMY_genWeightUp =                                        	   (*genWeightUp);                                       
  MY->bMY_genWeightDown =                                      	   (*genWeightDown);                                     
  MY->bMY_PRWHash =                                            	   (*PRWHash);                                           
  MY->bMY_EventNumber =                                        	   (*EventNumber);                                       
  MY->bMY_xsec =                                               	   (*xsec);                                              
  MY->bMY_GenHt =                                              	   (*GenHt);                                             
  MY->bMY_GenMET =                                             	   (*GenMET);                                            
  MY->bMY_DatasetNumber =                                      	   (*DatasetNumber);                                     
  MY->bMY_RunNumber =                                          	   (*RunNumber);                                         
  MY->bMY_RandomRunNumber =                                    	   (*RandomRunNumber);                                   
  MY->bMY_FS =                                                 	   (*FS);                                                

  MY->WriteTree();
  **/
}


void MySusySkimAnalysis::fill2L2JTree(int idx1, int idx2, std::map<TString, double> fake_weight){

 
  int nTT = 0;
  int nTl = 0;
  int nlT = 0;
  int nll = 0;
  classifyEvent(idx1, idx2, nTT, nTl, nlT, nll);

  twoLtwoJ->b2L2J_isTT = nTT;
  twoLtwoJ->b2L2J_isTl = nTl;
  twoLtwoJ->b2L2J_islT = nlT;
  twoLtwoJ->b2L2J_isll = nll;
  twoLtwoJ->b2L2J_r1 = glob_r1;
  twoLtwoJ->b2L2J_r2 = glob_r2;
  twoLtwoJ->b2L2J_f1 = glob_f1;
  twoLtwoJ->b2L2J_f2 = glob_f2;
  /**     
  twoLtwoJ->b2L2J_trigWeight_2LTrig =                                                           (*trigWeight_2LTrig);
  twoLtwoJ->b2L2J_trigMatch_2LTrig =                                                           (*trigMatch_2LTrigOR);
  twoLtwoJ->b2L2J_mu =                                                           (*mu);
  twoLtwoJ->b2L2J_avg_mu =                                                           (*avg_mu);
  twoLtwoJ->b2L2J_actual_mu =                                                           (*actual_mu);
  twoLtwoJ->b2L2J_nVtx =                                                           (*nVtx);
  twoLtwoJ->b2L2J_channel =                                                           (*channel);
  twoLtwoJ->b2L2J_nLep_combi =                                                           (*nLep_combi);
  twoLtwoJ->b2L2J_nLep_base =                                                           (*nLep_base);
  twoLtwoJ->b2L2J_nLep_base_OR = nlep_base_OR;
  twoLtwoJ->b2L2J_nLep_signal =                                                           (*nLep_signal);
  for(unsigned int i = 0; i<lepFlavor.GetSize(); i++){
    twoLtwoJ->b2L2J_lepFlavor.push_back                                                           (lepFlavor[i]);
    twoLtwoJ->b2L2J_lepCharge.push_back                                                           (lepCharge[i]);
    twoLtwoJ->b2L2J_lepAuthor.push_back                                                           (lepAuthor[i]);
    twoLtwoJ->b2L2J_lepPt.push_back                                                           (lepPt[i]);
    twoLtwoJ->b2L2J_lepEta.push_back                                                           (lepEta[i]);
    twoLtwoJ->b2L2J_lepPhi.push_back                                                           (lepPhi[i]);
    twoLtwoJ->b2L2J_lepIFFClass.push_back                                                           (lepIFFClass[i]);
    twoLtwoJ->b2L2J_lepM.push_back                                                           (lepM[i]);
    twoLtwoJ->b2L2J_lepD0.push_back                                                           (lepD0[i]);
    twoLtwoJ->b2L2J_lepD0Sig.push_back                                                           (lepD0Sig[i]);
    twoLtwoJ->b2L2J_lepZ0.push_back                                                           (lepZ0[i]);
    twoLtwoJ->b2L2J_lepZ0SinTheta.push_back                                                           (lepZ0SinTheta[i]);
    twoLtwoJ->b2L2J_lepPassOR.push_back                                                           (lepPassOR->at(i));
    twoLtwoJ->b2L2J_lepType.push_back                                                           (lepType[i]);
    twoLtwoJ->b2L2J_lepOrigin.push_back                                                           (lepOrigin[i]);
    twoLtwoJ->b2L2J_lepEgMotherType.push_back                                                           (lepEgMotherType[i]);
    twoLtwoJ->b2L2J_lepEgMotherOrigin.push_back                                                           (lepEgMotherOrigin[i]);
    twoLtwoJ->b2L2J_lepEgMotherPdgId.push_back                                                           (lepEgMotherPdgId[i]);
    twoLtwoJ->b2L2J_lepECIDS.push_back                                                           (lepECIDS[i]);
    twoLtwoJ->b2L2J_lepNPix.push_back                                                           (lepNPix[i]);
    twoLtwoJ->b2L2J_lepPassBL.push_back                                                           (lepPassBL->at(i));
    twoLtwoJ->b2L2J_lepSignal.push_back                                                           (lepSignal->at(i));
    twoLtwoJ->b2L2J_lepTruthMatched.push_back                                                           (lepTruthMatched->at(i));
    twoLtwoJ->b2L2J_lepTruthCharge.push_back                                                           (lepTruthCharge[i]);
    twoLtwoJ->b2L2J_lepTruthPt.push_back                                                           (lepTruthPt[i]);
    twoLtwoJ->b2L2J_lepTruthEta.push_back                                                           (lepTruthEta[i]);
    twoLtwoJ->b2L2J_lepTruthPhi.push_back                                                           (lepTruthPhi[i]);
    twoLtwoJ->b2L2J_lepTruthM.push_back                                                           (lepTruthM[i]);
  }
  twoLtwoJ->b2L2J_nJet30 =                                                           (*nJet30);
  twoLtwoJ->b2L2J_nJet20 =                                                           (*nJet20);
  twoLtwoJ->b2L2J_nBJet20_MV2c10_FixedCutBEff_77 =                                                           (*nBJet20_MV2c10_FixedCutBEff_77);
  for(unsigned int i = 0; i<jetPt.GetSize(); i++){
    twoLtwoJ->b2L2J_jetPt.push_back                                                           (jetPt[i]);
    twoLtwoJ->b2L2J_jetEta.push_back                                                           (jetEta[i]);
    twoLtwoJ->b2L2J_jetPhi.push_back                                                           (jetPhi[i]);
    twoLtwoJ->b2L2J_jetM.push_back                                                           (jetM[i]);
    twoLtwoJ->b2L2J_jetTileEnergy.push_back                                                              (jetTileEnergy[i]);
    twoLtwoJ->b2L2J_jetMV2c10.push_back                                                               (jetMV2c10[i]);
  }
  twoLtwoJ->b2L2J_mjj =                                                       			  (*mjj);                                                  
  twoLtwoJ->b2L2J_mjj_minDPhiZMET =                                           			  (*mjj_minDPhiZMET);                                      
  twoLtwoJ->b2L2J_Rjj =                                                       			  (*Rjj);                                                  
  twoLtwoJ->b2L2J_Rjj_minDPhiZMET =                                           			  (*Rjj_minDPhiZMET);                                      
  twoLtwoJ->b2L2J_vectorSumJetsPt =                                           			  (*vectorSumJetsPt);                                      
  twoLtwoJ->b2L2J_vectorSumJetsEta =                                          			  (*vectorSumJetsEta);                                     
  twoLtwoJ->b2L2J_vectorSumJetsPhi =                                          			  (*vectorSumJetsPhi);                                     
  twoLtwoJ->b2L2J_vectorSumJetsM =                                            			  (*vectorSumJetsM);                                       
  twoLtwoJ->b2L2J_met_Et =                                                    			  (*met_Et);                                               
  twoLtwoJ->b2L2J_met_Sign =                                                  			  (*met_Sign);                                             
  twoLtwoJ->b2L2J_met_Phi =                                                   			  (*met_Phi);                                              
  twoLtwoJ->b2L2J_TST_Et =                                                    			  (*TST_Et);                                               
  twoLtwoJ->b2L2J_TST_Phi =                                                   			  (*TST_Phi);                                              
  twoLtwoJ->b2L2J_deltaPhi_MET_TST_Phi =                                      			  (*deltaPhi_MET_TST_Phi);                                 
  twoLtwoJ->b2L2J_Ht30 =                                                      			  (*Ht30);                                                 
  twoLtwoJ->b2L2J_Rbb =                                                       			  (*Rbb);                                                  
  twoLtwoJ->b2L2J_mbb =                                                       			  (*mbb);                                                  
  twoLtwoJ->b2L2J_METOverPtZ =                                                			  (*METOverPtZ);                                           
  twoLtwoJ->b2L2J_METOverPtW =                                                			  (*METOverPtW);                                           
  twoLtwoJ->b2L2J_PtISR =                                                     			  (*PtISR);                                                
  twoLtwoJ->b2L2J_METOverPtISR =                                              			  (*METOverPtISR);                                         
  twoLtwoJ->b2L2J_METOverHT =                                                 			  (*METOverHT);                                            
  twoLtwoJ->b2L2J_METOverJ1pT =                                               			  (*METOverJ1pT);                                          
  twoLtwoJ->b2L2J_DPhiJ1Met =                                                 			  (*DPhiJ1Met);                                            
  twoLtwoJ->b2L2J_DPhiJ2Met =                                                 			  (*DPhiJ2Met);                                            
  twoLtwoJ->b2L2J_DPhiJ3Met =                                                 			  (*DPhiJ3Met);                                            
  twoLtwoJ->b2L2J_DPhiJ4Met =                                                 			  (*DPhiJ4Met);                                            
  twoLtwoJ->b2L2J_minDPhi2JetsMet =                                           			  (*minDPhi2JetsMet);                                      
  twoLtwoJ->b2L2J_minDPhi4JetsMet =                                           			  (*minDPhi4JetsMet);                                      
  twoLtwoJ->b2L2J_minDPhiAllJetsMet =                                         			  (*minDPhiAllJetsMet);                                    
  twoLtwoJ->b2L2J_dPhiPjjMet =                                                			  (*dPhiPjjMet);                                           
  twoLtwoJ->b2L2J_dPhiPjjMet_minDPhiZMET =                                    			  (*dPhiPjjMet_minDPhiZMET);                               
  twoLtwoJ->b2L2J_dPhiMetISR =                                                			  (*dPhiMetISR);                                           
  twoLtwoJ->b2L2J_dPhiMetJet1 =                                               			  (*dPhiMetJet1);                                          
  twoLtwoJ->b2L2J_METOverHTLep =                                              			  (*METOverHTLep);                                         
  twoLtwoJ->b2L2J_mll =                                                       			  (*mll);                                                  
  twoLtwoJ->b2L2J_mlll =                                                       			  (*mlll);                                                 
  twoLtwoJ->b2L2J_mllll =                                                       			  (*mllll);                                                
  twoLtwoJ->b2L2J_Rll =                                                       			  (*Rll);                                                  
  twoLtwoJ->b2L2J_Ptll =                                                      			  (*Ptll);                                                 
  twoLtwoJ->b2L2J_absEtall =                                                  			  (*absEtall);                                             
  twoLtwoJ->b2L2J_dPhiPllMet =                                                			  (*dPhiPllMet);                                           
  twoLtwoJ->b2L2J_dPhill =                                                    			  (*dPhill);                                               
  twoLtwoJ->b2L2J_mt2leplsp_0 =                                               			  (*mt2leplsp_0);                                          
  twoLtwoJ->b2L2J_pileupWeight =                                             			  (*pileupWeight);                                         
  twoLtwoJ->b2L2J_leptonWeight =                                             			  (*leptonWeight);                                         
  twoLtwoJ->b2L2J_eventWeight =                                              			  (*eventWeight);                                          
  twoLtwoJ->b2L2J_bTagWeight =                                               			  (*bTagWeight);                                           
  twoLtwoJ->b2L2J_jvtWeight =                                                			  (*jvtWeight);                                            
  twoLtwoJ->b2L2J_globalDiLepTrigSF =                                        			  (*globalDiLepTrigSF);                                    
  twoLtwoJ->b2L2J_genWeight =                                                			  (*genWeight);                                            
  twoLtwoJ->b2L2J_genWeightUp =                                              			  (*genWeightUp);                                          
  twoLtwoJ->b2L2J_genWeightDown =                                            			  (*genWeightDown);                                        
  twoLtwoJ->b2L2J_truthMll =                                                 			  (*truthMll);                                             
  twoLtwoJ->b2L2J_winoBinoMllWeight =                                        			  (*winoBinoMllWeight);                                    
  twoLtwoJ->b2L2J_winoBinoXsecWeight =                                       			  (*winoBinoXsecWeight);                                   
  twoLtwoJ->b2L2J_winoBinoBrFracWeight =                                     			  (*winoBinoBrFracWeight);                                 
  twoLtwoJ->b2L2J_winoBinoWeight =                                           			  (*winoBinoWeight);                                       
  twoLtwoJ->b2L2J_ttbarNNLOWeight =                                          			  (*ttbarNNLOWeight);                                      
  twoLtwoJ->b2L2J_ttbarNNLOWeightUp =                                        			  (*ttbarNNLOWeightUp);                                    
  twoLtwoJ->b2L2J_ttbarNNLOWeightDown =                                      			  (*ttbarNNLOWeightDown);                                  
  twoLtwoJ->b2L2J_truthTopPt =                                                			  (*truthTopPt);                                           
  twoLtwoJ->b2L2J_truthAntiTopPt =                                            			  (*truthAntiTopPt);                                       
  twoLtwoJ->b2L2J_truthTtbarPt =                                              			  (*truthTtbarPt);                                         
  twoLtwoJ->b2L2J_truthTtbarM =                                               			  (*truthTtbarM);                                          
  twoLtwoJ->b2L2J_x1 =                                                        			  (*x1);                                                   
  twoLtwoJ->b2L2J_x2 =                                                        			  (*x2);                                                   
  twoLtwoJ->b2L2J_pdf1 =                                                      			  (*pdf1);                                                 
  twoLtwoJ->b2L2J_pdf2 =                                                      			  (*pdf2);                                                 
  twoLtwoJ->b2L2J_scalePDF =                                                  			  (*scalePDF);                                             
  twoLtwoJ->b2L2J_id1 =                                                         			  (*id1);                                                  
  twoLtwoJ->b2L2J_id2 = 				                        			  (*id2);						 		 
  twoLtwoJ->b2L2J_nLeps_RJ = 		  							       (*nLeps_RJ);		  								 	 
  twoLtwoJ->b2L2J_nJets_RJ = 		  						    	  (*nJets_RJ);		  								 		 
  twoLtwoJ->b2L2J_nBtagJets_RJ = 		  					    	  (*nBtagJets_RJ);

  for(unsigned int i = 0; i<jetPt_RJ.GetSize(); i++){
  twoLtwoJ->b2L2J_jetPt_RJ.push_back(		  						    	 jetPt_RJ[i]);		  								 		 
  twoLtwoJ->b2L2J_jetEta_RJ.push_back(		  						    	 jetEta_RJ[i]);		  								 		 
  twoLtwoJ->b2L2J_jetPhi_RJ.push_back(		  						    	 jetPhi_RJ[i]);		  								 		 
  twoLtwoJ->b2L2J_jetM_RJ.push_back(		  						    	 jetM_RJ[i]);
  }
  for(unsigned int i = 0; i<lepPt_RJ.GetSize(); i++){
  twoLtwoJ->b2L2J_lepPt_RJ.push_back(  		  						    	  lepPt_RJ[i]);  		  								 		 
  twoLtwoJ->b2L2J_lepEta_RJ.push_back( 		  						    	  lepEta_RJ[i]); 		  								 		 
  twoLtwoJ->b2L2J_lepPhi_RJ.push_back( 		  						    	  lepPhi_RJ[i]); 		  								 		 
  twoLtwoJ->b2L2J_lepE_RJ.push_back(   		  						    	  lepE_RJ[i]);   		  								 		 
  twoLtwoJ->b2L2J_lepsign_RJ.push_back(		  						    	  lepsign_RJ[i]);
  }
  twoLtwoJ->b2L2J_is2Lep2Jet = 		  						    	  (*is2Lep2Jet);		  								 		 
  twoLtwoJ->b2L2J_is2L2JInt = 		  						    	  (*is2L2JInt);		  								 		 
  twoLtwoJ->b2L2J_is3Lep = 			  					    	  (*is3Lep);			  							 		 
  twoLtwoJ->b2L2J_is3LInt = 		  						    	  (*is3LInt);		  								 		 
  twoLtwoJ->b2L2J_is3Lep2Jet = 		  						    	  (*is3Lep2Jet);		  								 		 
  twoLtwoJ->b2L2J_is3Lep3Jet = 		  						    	  (*is3Lep3Jet);		  								 		 
  twoLtwoJ->b2L2J_is4Lep2Jet = 		  						    	  (*is4Lep2Jet);		  								 		 
  twoLtwoJ->b2L2J_is4Lep3Jet = 		  						    	  (*is4Lep3Jet);		  								 		 
  twoLtwoJ->b2L2J_mll_RJ = 			  					    	  (*mll_RJ);			  							 		 
  twoLtwoJ->b2L2J_H2PP = 			  					    	  (*H2PP);			  								 			 
  twoLtwoJ->b2L2J_H5PP = 			  					    	  (*H5PP);			  								 			 
  twoLtwoJ->b2L2J_RPT_HT5PP =                         						       (*RPT_HT5PP);                        							 	 
  twoLtwoJ->b2L2J_R_minH2P_minH3P =                   						       (*R_minH2P_minH3P);                  							 	 
  twoLtwoJ->b2L2J_R_H2PP_H5PP =                  						       (*R_H2PP_H5PP);                 								 		 
  twoLtwoJ->b2L2J_dphiVP =                    						    	  (*dphiVP);                   								 		 
  twoLtwoJ->b2L2J_minDphi =                   						    	  (*minDphi);                  								 		 
  twoLtwoJ->b2L2J_mTW =                       						    	  (*mTW);                      								 		 
  twoLtwoJ->b2L2J_H4PP =                      						    	  (*H4PP);                     								 		 
  twoLtwoJ->b2L2J_RPT_HT4PP =                 						    	  (*RPT_HT4PP);                								 		 
  twoLtwoJ->b2L2J_R_HT4PP_H4PP =              							       (*R_HT4PP_H4PP);             								 	 
  twoLtwoJ->b2L2J_PTISR =                     						    	  (*PTISR);                    								 		 
  twoLtwoJ->b2L2J_RISR =                      						    	  (*RISR);                     								 		 
  twoLtwoJ->b2L2J_PTI =                       						    	  (*PTI);                      								 		 
  twoLtwoJ->b2L2J_dphiISRI =                  						    	  (*dphiISRI);                 								 		 
  twoLtwoJ->b2L2J_PTCM =                      						    	  (*PTCM);                     								 		 
  twoLtwoJ->b2L2J_NjS =                       						    	  (*NjS);                      								 		 
  twoLtwoJ->b2L2J_NjISR =                     						    	  (*NjISR);                    								 		 
  twoLtwoJ->b2L2J_MZ =                   	  					    	  (*MZ);                  	  								 			 
  twoLtwoJ->b2L2J_MJ =                  	  					    	  (*MJ);                 	  								 			 
  twoLtwoJ->b2L2J_mTl3 =                 	  					    	  (*mTl3);                	  								 			 
  twoLtwoJ->b2L2J_lept1Pt_VR =                						    	  (*lept1Pt_VR);               								 		 
  twoLtwoJ->b2L2J_lept1sign_VR =              						    	  (*lept1sign_VR);             								 		 
  twoLtwoJ->b2L2J_lept2Pt_VR =                						    	  (*lept2Pt_VR);               								 		 
  twoLtwoJ->b2L2J_lept2sign_VR =              						    	  (*lept2sign_VR);             								 		 
  twoLtwoJ->b2L2J_mll_RJ_VR =            	  					    	  (*mll_RJ_VR);           	  								 			 
  twoLtwoJ->b2L2J_H2PP_VR =           	  						    	  (*H2PP_VR);          	  								 		 
  twoLtwoJ->b2L2J_H5PP_VR =          	  						    	  (*H5PP_VR);         	  								 		 
  twoLtwoJ->b2L2J_RPT_HT5PP_VR =         	  					    	  (*RPT_HT5PP_VR);        	  								 			 
  twoLtwoJ->b2L2J_R_minH2P_minH3P_VR =        						    	  (*R_minH2P_minH3P_VR);       								 		 
  twoLtwoJ->b2L2J_R_H2PP_H5PP_VR =       	  					    	  (*R_H2PP_H5PP_VR);      	  								 			 
  twoLtwoJ->b2L2J_dphiVP_VR =      		  					    	  (*dphiVP_VR);     		  							 		 
  twoLtwoJ->b2L2J_PTISR_VR =     		  					    	  (*PTISR_VR);    		  								 			 
  twoLtwoJ->b2L2J_RISR_VR =    		  						    	  (*RISR_VR);   		  								 		 
  twoLtwoJ->b2L2J_PTI_VR =   		  						    	  (*PTI_VR);  		  								 		 
  twoLtwoJ->b2L2J_dphiISRI_VR =  		  					    	  (*dphiISRI_VR); 		  								 			 
  twoLtwoJ->b2L2J_PTCM_VR = 		  						    	  (*PTCM_VR);		  								 		 
  twoLtwoJ->b2L2J_NjS_VR = 			  					    	  (*NjS_VR);			  							 		 
  twoLtwoJ->b2L2J_NjISR_VR = 		  						    	  (*NjISR_VR);		  								 		 
  twoLtwoJ->b2L2J_MZ_VR = 			  					    	  (*MZ_VR);			  							 		 
  twoLtwoJ->b2L2J_MJ_VR = 			  					    	  (*MJ_VR);			  							 		 
  twoLtwoJ->b2L2J_PRWHash = 		  						    	  (*PRWHash);		  								 		 
  twoLtwoJ->b2L2J_EventNumber = 		  					    	  (*EventNumber);		  								 			 
  twoLtwoJ->b2L2J_xsec = 			  					    	  (*xsec);			  								 			 
  twoLtwoJ->b2L2J_GenHt = 			  					    	  (*GenHt);			  							 		 
  twoLtwoJ->b2L2J_GenMET = 			  					    	  (*GenMET);			  							 		 
  twoLtwoJ->b2L2J_DatasetNumber = 		  					    	  (*DatasetNumber);		  							 		 
  twoLtwoJ->b2L2J_RunNumber = 		  						    	  (*RunNumber);		  								 		 
  twoLtwoJ->b2L2J_RandomRunNumber = 	  						    	  (*RandomRunNumber);	  								 		 
  twoLtwoJ->b2L2J_FS = 									    	  (*FS);                                                                                    


  double tot_err_up = 0.0;
  double tot_err_dw = 0.0;
  TString dwn_str = "";
  TString up_str = "";
  double nominal_weight = b_fake_wgt["NOM"];
  for(std::map<TString,double>::iterator iter = b_fake_wgt.begin(); iter != b_fake_wgt.end(); ++iter)
    {
      twoLtwoJ->b2L2J_FNPweight[iter->first] = iter->second;
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
  
  if(sqrt(tot_err_up) > 3){
    cout<<"UP: "<<up_str.Data()<<endl;
    cout<<"DW: "<<dwn_str.Data()<<endl;
    printf("TOT ERR UP = %.2f\n",sqrt(tot_err_up));
    for(std::map<TString,double>::iterator iter = b_fake_wgt.begin(); iter != b_fake_wgt.end(); ++iter)
      {
	printf("DEBUG \t %s =  %.2f\n",iter->first.Data(),iter->second);
      }
  }
  
  twoLtwoJ->b2L2J_FNP_TOTAL_UP   = nominal_weight + fabs(sqrt(tot_err_up)*nominal_weight);
  twoLtwoJ->b2L2J_FNP_TOTAL_DOWN = nominal_weight - fabs(sqrt(tot_err_dw)*nominal_weight);

  // twoLtwoJ->b2L2J_fake_wgt     = b_fake_wgt.size() > 0 ? b_fake_wgt.at(0) : 1.0;
  // twoLtwoJ->b2L2J_fake_wgt_up  = b_fake_wgt.size() > 1 ? b_fake_wgt.at(1) : 1.0;
  // twoLtwoJ->b2L2J_fake_wgt_dw  = b_fake_wgt.size() > 2 ? b_fake_wgt.at(2) : 1.0;
  */   
  twoLtwoJ->WriteTree();
  
 
}

std::vector<TString> MySusySkimAnalysis::getTriggerCat(int id, int yr){
  std::vector<TString> trigc;
  // cout<<"Numb of e triggers"<<trigmatch_e[yr].size()<<endl;
  // cout<<"Numb of m triggers"<<trigmatch_m[yr].size()<<endl;
  if(lepFlavor[id] == 1){
    if(trigmatch_e.find(yr) == trigmatch_e.end()){
      cout<<"ERROR \t Could not find "<<yr<<" in trigmatch_e"<<endl;
    }
    for(unsigned int i=0; i<trigmatch_e[yr].size();i++){
      if(i >= trigstr_e[yr].size()){
	cout<<"ERROR \t Entry "<< i << " is not in trigstr_e["<<yr<<"] which has size "<<trigstr_e[yr].size()<<endl;
      }
      if((trigmatch_e[yr].at(i))->at(id) && lepPt[id] > getPtThresholdTrigger(trigstr_e[yr].at(i),"e",0)){
	//printf("El(%.2f): %s : %s\n",lepPt[id],trigstr_e[yr].at(i).Data(),(trigmatch_e[yr].at(i))->at(id) ? "True" : "False");
	if((i+2) > (unsigned int)h_triggermatched_el->GetNbinsX()){
	  cout<<"EL: Adding 1 to bin "<<i+2<<" for histogram with "<<h_triggermatched_el->GetNbinsX()<<" bins"<<endl;
	}
	h_triggermatched_el->AddBinContent(i+2,1);
      }
      if((trigmatch_e[yr].at(i))->at(id) && lepPt[id] > getPtThresholdTrigger(trigstr_e[yr].at(i),"e",0) && !(std::find(trigc.begin(), trigc.end(), trigcat[trigstr_e[yr].at(i)]) != trigc.end())){
	trigc.push_back(trigcat[trigstr_e[yr].at(i)]);
      }
    }
    if(trigc.size() == 0)h_triggermatched_el->AddBinContent(1,1);
  }else if(lepFlavor[id] == 2){
    if(trigmatch_m.find(yr) == trigmatch_m.end()){
      cout<<"ERROR \t Could not find "<<yr<<" in trigmatch_m"<<endl;
    }
    for(unsigned int i=0; i<trigmatch_m[yr].size();i++){
      if(i >= trigstr_m[yr].size()){
	cout<<"ERROR \t Entry "<< i << " is not in trigstr_m["<<yr<<"] which has size "<<trigstr_m[yr].size()<<endl;
      }
      if((trigmatch_m[yr].at(i))->at(id) && lepPt[id] > getPtThresholdTrigger(trigstr_m[yr].at(i),"mu",0)){
	//printf("Mu(%.2f): %s : %s\n",lepPt[id],trigstr_m[yr].at(i).Data(),(trigmatch_m[yr].at(i))->at(id) ? "True" : "False");
	if((i+2) > (unsigned int)h_triggermatched_mu->GetNbinsX()){
	  cout<<"MU: Adding 1 to bin "<<i+2<<" for histogram with "<<h_triggermatched_mu->GetNbinsX()<<" bins"<<endl;
	}
	h_triggermatched_mu->AddBinContent(i+2,1);
      }
      if((trigmatch_m[yr].at(i))->at(id)  && lepPt[id] > getPtThresholdTrigger(trigstr_m[yr].at(i),"mu",0) && !(std::find(trigc.begin(), trigc.end(), trigcat[trigstr_m[yr].at(i)]) != trigc.end()))trigc.push_back(trigcat[trigstr_m[yr].at(i)]);
    }
    if(trigc.size() == 0)h_triggermatched_mu->AddBinContent(1,1);   
  }
  return trigc;
}

std::vector<TString> MySusySkimAnalysis::checkTriggerMatch(int id, int yr){
  std::vector<TString> matched_triggers; 
  if(lepFlavor[id] == 1){
    for(unsigned int i=0; i<trigmatch_e[yr].size();i++){
      Float_t pt_th = getPtThresholdTrigger(trigstr_e[yr].at(i),"e",0);
      //printf("Threshold for trigger %s is %.0f\n",trigstr_e[yr].at(i).Data(),pt_th);
      //printf("El(%.2f): %s : %s\n",lepPt[id],trigstr_e[yr].at(i).Data(),(trigmatch_e[yr].at(i))->at(id) ? "True" : "False");
      if((trigmatch_e[yr].at(i))->at(id) && lepPt[id] > pt_th)matched_triggers.push_back(trigstr_e[yr].at(i));
    }
  }else if(lepFlavor[id] == 2){
    for(unsigned int i=0; i<trigmatch_m[yr].size();i++){
      Float_t pt_th = getPtThresholdTrigger(trigstr_m[yr].at(i),"mu",0);
      //printf("Threshold for trigger %s is %.0f\n",trigstr_m[yr].at(i).Data(),pt_th);
      //printf("Mu(%.2f): %s : %s\n",lepPt[id],trigstr_m[yr].at(i).Data(),(trigmatch_m[yr].at(i))->at(id) ? "True" : "False");
      if((trigmatch_m[yr].at(i))->at(id) && lepPt[id] > pt_th)matched_triggers.push_back(trigstr_m[yr].at(i));
    }
  }

  return matched_triggers;
}
