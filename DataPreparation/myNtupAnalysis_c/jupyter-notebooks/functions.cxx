#include <ROOT/RVec.hxx>
#include "../../NTUPana_SusysSkim/CalcGenericMT2/CalcGenericMT2/MT2_ROOT.h"

//R.gInterpreter.AddIncludePath("../NTUPana_SusysSkim/CalcGenericMT2/CalcGenericMT2/");
using Vec2_t = const ROOT::VecOps::RVec<float>;
using Vec_t = const ROOT::VecOps::RVec<int>;
  
bool myfilter(float x) {
   return x > 5;
}

bool isOS(const ROOT::VecOps::RVec<int>& chlep) {
  if(chlep[0]*chlep[1] < 0)return kTRUE;
  return kFALSE;
}

bool isSF(Vec_t& fllep) {
    if(fllep[0] == fllep[1])return kTRUE;
    return kFALSE;
}

float ComputeInvariantMass(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e) {
  TLorentzVector p1;
  TLorentzVector p2;
  p1.SetPtEtaPhiM(pt[0], eta[0], phi[0], e[0]);
  p2.SetPtEtaPhiM(pt[1], eta[1], phi[1], e[1]);
  return (p1 + p2).M();
}


float calcMT2(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e, Float_t met_et, Float_t met_phi) {

    TLorentzVector p1;
    TLorentzVector p2;
    TLorentzVector met;
    p1.SetPtEtaPhiM(pt[0], eta[0], phi[0], e[0]);
    p2.SetPtEtaPhiM(pt[1], eta[1], phi[1], e[1]);
    met.SetPtEtaPhiM(met_et, 0.0, met_phi, 0.0);
    return ComputeMT2(p1,p2,met,0.,0.).Compute();
}

float ptllboost(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e, Float_t met_et, Float_t met_phi) {

    TLorentzVector p1;
    TLorentzVector p2;
    TLorentzVector met;
    p1.SetPtEtaPhiM(pt[0], eta[0], phi[0], e[0]);
    p2.SetPtEtaPhiM(pt[1], eta[1], phi[1], e[1]);
    met.SetPtEtaPhiM(met_et, 0.0, met_phi, 0.0);
    return (met+p1+p2).Pt();
}

float costhetastar(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e) {

    TLorentzVector p1;
    TLorentzVector p2;
    p1.SetPtEtaPhiM(pt[0], eta[0], phi[0], e[0]);
    p2.SetPtEtaPhiM(pt[1], eta[1], phi[1], e[1]);
    return TMath::ATan(fabs(p1.Eta()-p2.Eta())/2.);
}

float deltaPhi_ll(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e) {

    TLorentzVector p1;
    TLorentzVector p2;
    p1.SetPtEtaPhiM(pt[0], eta[0], phi[0], e[0]);
    p2.SetPtEtaPhiM(pt[1], eta[1], phi[1], e[1]);
    return p1.DeltaPhi(p2);
}

float deltaPhi_metl(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e, Float_t met_et, Float_t met_phi) {

    TLorentzVector p1;
    TLorentzVector p2;
    TLorentzVector met;
    p1.SetPtEtaPhiM(pt[0], eta[0], phi[0], e[0]);
    p2.SetPtEtaPhiM(pt[1], eta[1], phi[1], e[1]);
    met.SetPtEtaPhiM(met_et, 0.0, met_phi, 0.0);
    
    if(p1.Pt() > p2.Pt()){
        return p1.DeltaPhi(met);
    }else{
        return p2.DeltaPhi(met);
    }
}

bool checkPt(Vec2_t& pt, float cut1, float cut2){
    if((pt[0] > cut1 && pt[1] > cut2) || (pt[1] > cut1 && pt[0] > cut2))return kTRUE;
    return kFALSE;
}
