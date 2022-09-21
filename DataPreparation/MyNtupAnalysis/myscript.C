#include <TROOT.h>
#include <TChain.h>
#include <TObjString.h>
#include <iostream>
#include <string>
#include <filesystem>
#include <unistd.h>
//#include "MySusySkimAnalysis.h"


void myscript(TString indir, TString key, TString year, TString date_folder, int sumentries = 0, int procnum = -1, int procstart = -1){

  gROOT->ProcessLine(".L HistFitterTree.cxx+");
  gROOT->ProcessLine(".L make2L2JTree.cxx+");
  gROOT->ProcessLine(".L makeMYTree.cxx+");
  gROOT->ProcessLine(".L GetFakeWeight.C+");
  char tmp[256];
  getcwd(tmp, 256);
  std::cout<<"Setting macro path to "<<tmp<<std::endl;
  gROOT->SetMacroPath(tmp);
  TChain *c;
  TString physproc = key;
  //if(!key.Contains("data")){
  TObjArray *tx = key.Tokenize("_");
  physproc = (((TObjString *)tx->At(0))->String());
  //}
  
  cout<<"hei"<<endl;

  if(!key.Contains("data")){
    c = new TChain(Form("%s_NoSys",physproc.Data()));
  }else{
    c = new TChain(physproc.Data());
  }
  // if(key.Contains("data")){
  //   c->Add(Form("/storage/eirikgr/ANAntuples/PHYS_Data/%s_*_merged_processed.root",key));
  // }else{
  //   if(year.EqualTo("year1516")){
  std::cout<<Form("%s/%s_*_merged_processed.root",indir.Data(),physproc.Data())<<std::endl;
  c->Add(Form("%s/%s_*merged_processed.root",indir.Data(),physproc.Data()));
  // }else if(year.EqualTo("year17")){
  //   c->Add(Form("/storage/eirikgr/ANAntuples/EXOT0_Bkgs_mc16d/%s_*_merged_processed.root",key));
  // }else if(year.EqualTo("year18")){
  //   c->Add(Form("/storage/eirikgr/ANAntuples/EXOT0_Bkgs_mc16e/%s_*_merged_processed.root",key));
  // }else{
  //   std::cout<<"ERROR \t Unknown year "<<year<<std::endl;
  //   return;
  // }
  //}
  cout<<"hade"<<endl;

  if(c->GetEntries() != sumentries){
    std::cout<<"ERROR \t Inconsistency in number of events in TTree ("<<c->GetEntries() <<") and number of events used in splitting ("<<sumentries<<")"<<std::endl;
    return;
  }

  if(sumentries and procnum and procstart >= 0){
    std::cout<<Form("data %s %i %i %s",year.Data(),procnum,procstart,date_folder.Data())<<std::endl;
    if(key.Contains("data")){
      c->Process("MySusySkimAnalysis.C++", Form("data %s %i %i %s",year.Data(),procnum,procstart,date_folder.Data()),procnum,procstart);
    }else{
      c->Process("MySusySkimAnalysis.C++", Form("mc %s %i %i %s",year.Data(),procnum,procstart,date_folder.Data()),procnum,procstart);
    }
  }else{
    if(key.Contains("data")){
      c->Process("MySusySkimAnalysis.C++", Form("data %s %s",year.Data(),date_folder.Data()),procnum,procstart);
    }else{
      c->Process("MySusySkimAnalysis.C++", Form("mc %s %s",year.Data(),date_folder.Data()),procnum,procstart);
    }
  }
  
}
