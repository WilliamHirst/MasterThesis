#ifndef Cfunctions_h
#define Cfunctions_h

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
#include <ROOT/RVec.hxx>
using Vec2_t = const ROOT::VecOps::RVec<float>;
using Vec_t = const ROOT::VecOps::RVec<int>;

bool  myfilter(float x);
bool  isOS(const ROOT::VecOps::RVec<int>& chlep);
bool  isSF(Vec_t& fllep);
float ComputeInvariantMass(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e);
float calcMT2(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e, Float_t met_et, Float_t met_phi);
float ptllboost(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e, Float_t met_et, Float_t met_phi);
float costhetastar(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e);
float deltaPhi_ll(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e);
float deltaPhi_metl(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e, Float_t met_et, Float_t met_phi);
bool  checkPt(Vec2_t& pt, float cut1, float cut2);

#endif
