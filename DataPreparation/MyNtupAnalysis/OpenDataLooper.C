#define OpenDataLooper_cxx
// The class definition in OpenDataLooper.h has been generated automatically
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
// root> T->Process("OpenDataLooper.C")
// root> T->Process("OpenDataLooper.C","some options")
// root> T->Process("OpenDataLooper.C+")
//


#include "OpenDataLooper.h"
#include <TH2.h>
#include <TStyle.h>
#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/RInterface.hxx"


void OpenDataLooper::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

  ROOT::RDataFrame rdf("analysis","/storage/eirikgr/user.egramsta.LPXReweightingTool_grid_test3_ANALYSIS.root/*.root");
  auto evnum = rdf.Take<ULong64_t>("EventNumber");
  auto lpx_weight = rdf.Take<Double_t>("weight");
  //auto truemass = filter.Take<Int_t>("truthmass");
  
  auto vecto_evnum = evnum.GetValue();
  auto vecto_lpx_weight = lpx_weight.GetValue();

  int i = 0;
  for(const auto& value: vecto_evnum) {
    lpx_weights[value] = vecto_lpx_weight.at(i);
    i += 1;
  }
  
  nentry = 0;
  nentries = 0;
  fullpath = "";
  previous_fullpath = "";
  nfilled = 0;
  isreported = false;
  TString option = GetOption();
}

void OpenDataLooper::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t OpenDataLooper::Process(Long64_t entry)
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

  fReader.SetLocalEntry(entry);

  
  // If updating xsection
  //if(*XSection != 1.0 && *SumWeights != 1.0){
  //  return kTRUE;
  //}

  fullpath = ((TChain*)fReader.GetTree())->GetFile()->GetName();

  if(nentry == 0 || !fullpath.EqualTo(previous_fullpath)){

    std::cout<<"Previous file = "<< previous_fullpath.Data()<<std::endl;
    std::cout<<"New file = "<< fullpath.Data()<<std::endl;
    std::cout<<"nentry = "<<nentry<<endl;
    
    if(nentry > 0){
      WriteToFile();
    }
    
    previous_fullpath = fullpath;
    TObjArray *tx = fullpath.Tokenize("/");
    TString fname_temp;
    TString inputfile = ((TObjString *)tx->Last())->String();
    TString path = "/";
  
    for(int i = 0; i<=(tx->GetEntries()-2); i++){
      path += ((TObjString *)tx->At(i))->String();
      path += "/";
    }

    delete tx;
   
    newfile = inputfile.ReplaceAll(".root","_new.root");

    chain = new TChain(((TChain*)fReader.GetTree())->GetName());
    chain->Add(((TChain*)fReader.GetTree())->GetFile()->GetName());
    
    sprintf(fnameoutputcopyfilename,"%s/%s",path.Data(),newfile.Data());
    std::cout<<"Opening file  = "<<fnameoutputcopyfilename<<std::endl;
    ntupf = new TFile(fnameoutputcopyfilename,"RECREATE");
    newtree=chain->CloneTree(0);

    newtree->SetBranchAddress("XSection",&newxsec);
    newtree->SetBranchAddress("SumWeights",&newsumofw);

    bwgt = newtree->Branch("lpxweight",&lpx_wgt);
    //newtree->SetBranchAddress("lpxweight",&lpx_wgt);

    isreported = false;

    //nentries = 0;
    //nentry = 0;
  }

  if(!nentries){
      nentries = ((TChain*)fReader.GetTree())->GetEntries();
    }

  nentry += 1;

  if(nentry%100000 == 0)printf("Doing entry %i/%i\n",nentry,nentries);
  
  // std::cout<<"XSection = "<<*XSection<<std::endl; 
  // std::cout<<"SumWeights = "<<*SumWeights<<std::endl; 

  //if(*XSection == 1.0 || *SumWeights == 1){
    if(newfile.Contains("410155")){
      newxsec = 5.4822e-4*1.096*1000.;
      newsumofw = 4075279.75386;
    }else if(newfile.Contains("410218")){
      newxsec = 3.6864e-5*1.12*1000.;
      newsumofw = 51968.9384584;
    }else if(newfile.Contains("410219")){
      newxsec = 3.6868e-5*1.12*1000.;
      newsumofw = 52007.5311319;
    }else if(!isreported){
      std::cout<<"Could not find information for "<<newfile.Data()<<std::endl;
      isreported = true;
    }
    //}

  chain->GetEntry( fReader.GetCurrentEntry() );

  lpx_wgt = -999;
  if ( lpx_weights.find(*eventNumber) != lpx_weights.end() ) {
    lpx_wgt = lpx_weights[*eventNumber];
  }

  bwgt->Fill();
  newtree->Fill();

  nfilled += 1;

  //if(nentry > 10)// {
  //   Abort("hei");
  //   ntupf->Write();
  //   ntupf->Close();
  //   // //delete tree3;
  //   delete ntupf;
    

  // }



  

  return kTRUE;
}

void OpenDataLooper::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

  
}

void OpenDataLooper::WriteToFile()
{

  if(nfilled>0){
    std::cout<<"Writing "<< nfilled << " events to  = "<<fnameoutputcopyfilename<<std::endl;
    nfilled = 0;
    ntupf->Write();
    ntupf->Close();
    // //delete tree3;
    //delete newtree;
    //delete chain;
    delete ntupf;
  }else{
    std::cout<<"Did not find any ntuples with bad branches :-)"<<std::endl;
  }
  
}

void OpenDataLooper::Terminate()
{
  // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

  WriteToFile();

}
