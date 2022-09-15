#ifndef helperFunctions_h
#define helperFunctions_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include <ctime>
#include <dirent.h>
#include "TLorentzVector.h"
#include "TParameter.h"
#include <ROOT/RVec.hxx>
using VecF_t = const ROOT::RVec<float>&;
using VecI_t = const ROOT::RVec<int>&;
using VecB_t = const ROOT::VecOps::RVec<bool>;
using Vec2_t = const ROOT::VecOps::RVec<float>;
using Vec_t = const ROOT::VecOps::RVec<int>;

//TCanvas c("c","x hist");
std::string progressBar;
float nEvents;
int everyN;

clock_t c_begin = 0;
clock_t c_end = 0;
// void updateHistogram(TH1D &h_){		     
//   c.cd();
//   h_.Draw();
//   c.Update();    
// };

void setRunParameters(float n, int en){
  nEvents = n;
  everyN = en;
}

void printProgressBar(TH1D &h_){
  double elapsed_secs = 0;
  double ev_per_sec = 0;
  double time_left = 0.0;
  double minutesRemainder = 0.0;
  double secondsRemainder = 0.0;
  int hours = 0;
  int minutes = 0;
  int seconds = 0;
  double all_nev = progressBar.size()*everyN;
  c_end = clock();
  elapsed_secs = double(c_end - c_begin) / CLOCKS_PER_SEC;
  if(elapsed_secs > 0){
    ev_per_sec = ((double)everyN)/elapsed_secs;
    if(ev_per_sec > 0){
      time_left = (nEvents-all_nev)/ev_per_sec;
      time_left = time_left/(60.*60);
      hours = time_left;
      minutesRemainder = (time_left - hours) * 60;
      minutes = minutesRemainder;
      secondsRemainder = (minutesRemainder - minutes) * 60;
      seconds = secondsRemainder;
    }
    std::cout<<"Events/sec = "<<ev_per_sec<<"). "<<"Estimated time left = "<<hours<<"h"<<minutes<<"m"<<seconds<<"s"<<std::endl;
  }
  std::mutex barMutex; // Only one thread at a time can lock a mutex. Let's use this to avoid concurrent printing.
  // Magic numbers that yield good progress bars for nSlots = 1,2,4,8
  const auto barWidth = nEvents / everyN;
  std::lock_guard<std::mutex> l(barMutex); // lock_guard locks the mutex at construction, releases it at destruction
  progressBar.push_back('#');
  // re-print the line with the progress bar
  std::cout << "\r[" << std::left << std::setw(barWidth) << progressBar << ']' << std::flush;
  c_begin = clock();
}

//auto drawHisto = updateHistogram;

int Zlep1 = -1;
int Zlep2 = -1;
int Wlep1 = -1;
std::pair <double,double> getLeptonsFromZ(VecI_t chlep, VecI_t& fllep, VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e, Float_t met_et, Float_t met_phi);
bool  myfilter(float x);
bool  isOS(const ROOT::VecOps::RVec<int>& chlep);
bool  isSF(VecI_t& fllep);
bool isEE(VecI_t& fllep);
bool isMM(VecI_t& fllep);
double getSF(VecF_t& sf);
int flavourComp3L(VecI_t& fllep);
bool deltaRlepjet(float lpt, float leta, float lphi, float le, VecF_t& jpt, VecF_t& jeta, VecF_t& jphi, VecF_t& je);
float ComputeInvariantMass(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& m);
float OSSFInvariantMass(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& m, VecF_t& ch, VecF_t& tp);
float calcMT2(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e, Float_t met_et, Float_t met_phi, int idx1 = 0, int idx2 = 1);
float ptllboost(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e, Float_t met_et, Float_t met_phi);
float costhetastar(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e);
float deltaPhi_ll(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e);
float deltaPhi_metl(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e, Float_t met_et, Float_t met_phi);
bool  checkPt(VecF_t& pt, float cut1, float cut2, float cut3);
float deltaPhi_metll(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e, Float_t met_et, Float_t met_phi);
std::pair <int,int> num_bl_sg_lep(VecF_t& pt, VecF_t& eta, VecI_t& fllep, VecB_t passOR, VecB_t passLOOSE, VecB_t passMEDIUM, VecB_t passBL, VecF_t z0sinth, VecB_t ISO, VecF_t d0sig, VecF_t passTIGHT);
Double_t getMetRel(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e, Float_t met_et, Float_t met_phi);
int getZlep1(){return Zlep1;}
int getZlep2(){return Zlep2;}
int getWlep1(){return Wlep1;}
Float_t getLumiSF(Int_t randrnum);
float getVar(Vec2_t& var, int idx);
int getTypeTimesCharge(Vec_t& chlep, Vec_t& fllep, int idx);
float getE(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& m, int idx);
float getM(VecF_t& m, int idx);
float getMt(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& m, int idx);
float deltaR(VecF_t& eta, VecF_t& phi, int idx1, int idx2);
float SSHt(VecF_t& pt, VecF_t& ch);
float S_Et_miss(Float_t et_miss, VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& m);
bool isNGoodLeptons(VecI_t& isGood, int nrLeps);






#endif
