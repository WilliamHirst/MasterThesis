#include <TROOT.h>
#include <TChain.h>
#include <TObjString.h>
#include <iostream>
#include <string>
#include <filesystem>
#include <unistd.h>
//#include "MySusySkimAnalysis.h"


void myscript_simple(){
  /**
  gROOT->ProcessLine(".L HistFitterTree.cxx+");
  gROOT->ProcessLine(".L make2L2JTree.cxx+");
  gROOT->ProcessLine(".L makeMYTree.cxx+");
  gROOT->ProcessLine(".L GetFakeWeight.C+");
  char tmp[256];
  getcwd(tmp, 256);
  std::cout<<"Setting macro path to "<<tmp<<std::endl;
  gROOT->SetMacroPath(tmp);
  */  
  TChain *c = new TChain("Wjets_NoSys");
  c->Add("/storage/eirikgr/ANAntuples/EXOT0_Bkgs_mc16e/Wjets_merged_processed.root");
  c->Process("MySusySkimAnalysis.C++");


  
}
