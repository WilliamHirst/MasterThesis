#define helperFunctions_cxx


#include <ROOT/RVec.hxx>
#include "../CalcGenericMT2/CalcGenericMT2/MT2_ROOT.h"

using VecF_t = const ROOT::RVec<float>&;
using VecI_t = const ROOT::RVec<int>&;
using VecB_t = const ROOT::VecOps::RVec<bool>;

#include "helperFunctions.h"


bool myfilter(float x) {
   return x > 5;
}

auto sum = [](int a, int b) {
        return a + b;
    };


// int num_baseline_lep(VecF_t& pt, VecF_t& eta, VecI_t& fllep, VecB_t passOR, VecB_t passLOOSE, VecB_t passMEDIUM, VecB_t passBL, VecF_t z0sinth){
//   int nbl = 0;
//   for(unsigned int i=0; i<fllep.size(); i++)
//     {
//       if(pt[i] < 9)continue;
//       if((fllep[i] == 1 && fabs(eta[i])>2.47) || ((fllep[i] == 2 && fabs(eta[i])>2.6)))continue;
//       if(!passOR[i])continue
//       if((fllep[i] == 1 && (!passLOOSE[i] || !passBL[i])) || (fllep[i] == 22 && !passMEDIUM[i]))continue;
//       if(fabs(z0sinth)<0.5)continue;
//       nbl += 1;
//     }
//   return nbl;
// }

Double_t getMetRel(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e, Float_t met_et, Float_t met_phi){
  TLorentzVector l;
  Double_t min_dphi_lep_met = 9999;
  TLorentzVector met;
  met.SetPtEtaPhiM(met_et, 0.0, met_phi, 0.0);
  for(unsigned int i=0; i<pt.size(); i++)
    {
      l.SetPtEtaPhiM(pt[i], eta[i], phi[i], e[i]);
      Double_t dphi = fabs(l.DeltaPhi(met));
      if(dphi < min_dphi_lep_met){
	min_dphi_lep_met = dphi;
      }
    }
  return (min_dphi_lep_met < M_PI/2.0) ? (met_et)*sin(min_dphi_lep_met) : (met_et);
 
}

std::pair <int,int> num_bl_sg_lep(VecF_t& pt, VecF_t& eta, VecI_t& fllep, VecB_t passOR, VecB_t passLOOSE, VecB_t passMEDIUM, VecB_t passBL, VecF_t z0sinth, VecB_t ISO, VecF_t d0sig, VecF_t passTIGHT){
  int nbl = 0;
  int nsg = 0;
  std::pair <int,int> nlep;
  for(unsigned int i=0; i<fllep.size(); i++)
    {
      if(pt[i] < 9)continue;
      if((fllep[i] == 1 && fabs(eta[i])>2.47) || ((fllep[i] == 2 && fabs(eta[i])>2.6)))continue;
      if(!passOR[i])continue;
      if((fllep[i] == 1 && (!passLOOSE[i] || !passBL[i])) || (fllep[i] == 2 && !passMEDIUM[i]))continue;
      if(fabs(z0sinth[i])>0.5)continue;
      
      nbl += 1;
      
      if((fllep[i] == 1 && !passTIGHT[i]))continue;
      if(!ISO[i])continue;
      if((fllep[i] == 1 && fabs(d0sig[i])>5) || ((fllep[i] == 2 && fabs(d0sig[i])>3)))continue;
      
      nsg += 1;
    }
  nlep = std::make_pair(nbl,nsg);
  return nlep;
}


std::pair <double,double> getLeptonsFromZ(VecI_t chlep, VecI_t& fllep, VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e, Float_t met_et, Float_t met_phi){
  double diff = 10000000000.0;
  /**
  int Zlep1 = -99;
  int Zlep2 = -99;
  int Wlep1 = -999;
  */
  double Zmass = -1.0;
  double Wmass = -1.0;
  bool foundSFOS = false;
  std::pair <double,double> masses;
  for(unsigned int i=0; i<chlep.size(); i++)
    {
      for(unsigned int j=i+1; j<chlep.size(); j++)
	{
	  //Opposite-Sign
	  if(chlep[i]*chlep[j]<0)
	    {
	      //Same-Flavor
	      if(abs(fllep[i])==abs(fllep[j]))
		{
		  TLorentzVector p1;
		  p1.SetPtEtaPhiM(pt[i], eta[i], phi[i], e[i]);
		  TLorentzVector p2;
		  p2.SetPtEtaPhiM(pt[j], eta[j], phi[j], e[j]);
		  double mass = (p1+p2).M();
		  double massdiff = fabs(mass-91187.6);
		  if(massdiff<diff)
		    {
		      diff=massdiff;
		      Zmass=mass;
		      Zlep1 = i;
		      Zlep2 = j;
		      foundSFOS = true;
		    }
		}
	    }
	}

    }
  
  if(foundSFOS){
    TLorentzVector met;
    met.SetPtEtaPhiM(met_et, 0.0, met_phi, 0.0);
    
    if((Zlep1==0 && Zlep2==1) || (Zlep1==1 && Zlep2==0) ) Wlep1=2;
    else if((Zlep1==0 && Zlep2==2) || (Zlep1==2 && Zlep2==0) ) Wlep1=1;
    else if((Zlep1==1 && Zlep2==2) || (Zlep1==2 && Zlep2==1) ) Wlep1=0;
    
    TLorentzVector lepW;
    lepW.SetPtEtaPhiM(pt[Wlep1], eta[Wlep1], phi[Wlep1], e[Wlep1]);
    double wlepMetphi = lepW.DeltaPhi(met);
    Wmass = sqrt(2*lepW.Pt()*met.Pt()*(1-cos(wlepMetphi)));
  }
  masses = std::make_pair(Zmass,Wmass);
    
  return masses;
}

Float_t getLumiSF(Int_t randrnum){
    Float_t lumi15 = 3219.56;
    Float_t lumi16 = 32988.1;
    Float_t lumi17 = 44307.4;
    Float_t lumi18 = 58450.1;
    if(randrnum < 320000)return lumi15+lumi16;
    else if(randrnum > 320000 && randrnum < 348000)return lumi17;
    else if(randrnum > 348000)return lumi18;
    else{std::cout<<"ERROR \t RandomRunnumber "<<randrnum<<" has no period attached"<<std::endl;}
    return 1.0;
}

double getSF(VecF_t& sf){
  const auto size = sf.size();
  double scalef = 1.0;
  for (size_t i=0; i < size; ++i) {
    scalef *= sf[i];
  }
  return scalef;
}

bool isOS(const ROOT::VecOps::RVec<int>& chlep) {
  if(chlep[0]*chlep[1] < 0)return kTRUE;
  return kFALSE;
}

// bool isTriggerMatched(const ROOT::VecOps::RVec<int>& isTM) {
//   std::vector<int> tm_vec; 
//   const auto size = isTM.size();
//    for (size_t i=0; i < size; ++i) {
//      if(isTM[i])tm_vec.push_back(0);
//      else tm_vec.push_back(1);
//    }
   
//   return kFALSE;
// }

int flavourComp3L(VecI_t& fllep) {
  const auto size = fllep.size();
  if(size>=3){
    //std::cout<<"ERROR \t Vector must be at least 3 long!"<<std::endl;
    if(fllep[0] == 1 && fllep[1] == 1 && fllep[2] == 1)return 0;
    if(fllep[0] == 1 && fllep[1] == 1 && fllep[2] == 2)return 1;
    if(fllep[0] == 1 && fllep[1] == 2 && fllep[2] == 2)return 2;
    if(fllep[0] == 2 && fllep[1] == 2 && fllep[2] == 2)return 3;
    if(fllep[0] == 2 && fllep[1] == 2 && fllep[2] == 1)return 4;
    if(fllep[0] == 2 && fllep[1] == 1 && fllep[2] == 1)return 5;
  }else if(size==2){
    if(fllep[0] == 1 && fllep[1] == 1)return 6;                                                                                                                                                                                                                    
    if(fllep[0] == 2 && fllep[1] == 2)return 7;                                                                                                                                                                                                                    
    if((fllep[0] == 1 && fllep[1] == 2) || (fllep[0] == 2 && fllep[1] == 1))return 8;                                                                                                                                                                                               
  }
  return -1;
}

bool deltaRlepjet(float lpt, float leta, float lphi, float le, VecF_t& jpt, VecF_t& jeta, VecF_t& jphi, VecF_t& je){
  const auto njet = int(jpt.size());
  TLorentzVector p2;
  TLorentzVector p1;
  double deltaR;
  double mindr = 9999;
  //for (size_t i=0; i < nlep; ++i) {
  p1.SetPtEtaPhiM(lpt, leta, lphi, le);
  for (int j=0; j < njet; ++j) {
    p2.SetPtEtaPhiM(jpt[j], jeta[j], jphi[j], je[j]);
    deltaR = p1.DeltaR(p2);
    if(deltaR < mindr){
      mindr = deltaR;
    }
    //}
  }
  //  std::cout<<"mindr = "<<mindr<<std::endl;
  return mindr;  
}


bool isSF(VecI_t& fllep) {
    if(fllep[0] == fllep[1])return kTRUE;
    return kFALSE;
}

bool isEE(VecI_t& fllep) {
    if(fllep[0] == fllep[1] && fllep[0] == 1)return kTRUE;
    return kFALSE;
}

bool isMM(VecI_t& fllep) {
    if(fllep[0] == fllep[1] && fllep[0] == 2)return kTRUE;
    return kFALSE;
}

float ComputeInvariantMass(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& m) {
  TLorentzVector pf;
  const auto size = int(pt.size());
  for(int i = 0; i++; i<size){
      TLorentzVector pi;
      pi.SetPtEtaPhiM(pt[i], eta[i], phi[i], m[i]);
     pf += pi;
  }
  return pf.M();
}

float OSSFInvariantMass(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& m, VecF_t& ch, VecF_t& tp){
    int idx1;
    int idx2;
    if (tp[0] == tp[1] &&  ch[0] != ch[1]){
        idx1 = 0;
        idx2 = 1;
    }
    else if (tp[0] == tp[2] &&  ch[0] != ch[2]){
        idx1 = 0;
        idx2 = 2;
    }
    else if (tp[1] == tp[2] &&  ch[1] != ch[2]){
        idx1 = 1;
        idx2 = 2;
    }
    else{
        return 0;
    }
    TLorentzVector p1;
    TLorentzVector p2;
    p1.SetPtEtaPhiM(pt[idx1], eta[idx1], phi[idx1], m[idx1]);
    p2.SetPtEtaPhiM(pt[idx2], eta[idx2], phi[idx2], m[idx2]);
    return (p1 + p2 ).M();
    
    
}


float calcMT2(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e, Float_t met_et, Float_t met_phi, int idx1, int idx2) {

  const auto size = int(pt.size());
  if(idx1 > size || idx2 > size){
    printf("calcMT2::ERROR \t Indices %i and %i are higher than size of vector %i\n",idx1,idx2,size);
    return -1;
  }
  TLorentzVector p1;
  TLorentzVector p2;
  TLorentzVector met;
  p1.SetPtEtaPhiM(pt[idx1], eta[idx1], phi[idx1], e[idx1]);
  p2.SetPtEtaPhiM(pt[idx2], eta[idx2], phi[idx2], e[idx2]);
  met.SetPtEtaPhiM(met_et, 0.0, met_phi, 0.0);
  return ComputeMT2(p1,p2,met,0.,0.).Compute();
}

float ptllboost(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e, Float_t met_et, Float_t met_phi) {

    TLorentzVector p1;
    TLorentzVector p2;
    TLorentzVector met;
    p1.SetPtEtaPhiM(pt[0], eta[0], phi[0], e[0]);
    p2.SetPtEtaPhiM(pt[1], eta[1], phi[1], e[1]);
    met.SetPtEtaPhiM(met_et, 0.0, met_phi, 0.0);
    return (met+p1+p2).Pt();
}

float costhetastar(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e) {

    TLorentzVector p1;
    TLorentzVector p2;
    p1.SetPtEtaPhiM(pt[0], eta[0], phi[0], e[0]);
    p2.SetPtEtaPhiM(pt[1], eta[1], phi[1], e[1]);
    return TMath::ATan(fabs(p1.Eta()-p2.Eta())/2.);
}

float deltaPhi_ll(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e) {

    TLorentzVector p1;
    TLorentzVector p2;
    p1.SetPtEtaPhiM(pt[0], eta[0], phi[0], e[0]);
    p2.SetPtEtaPhiM(pt[1], eta[1], phi[1], e[1]);
    return p1.DeltaPhi(p2);
}

float deltaPhi_metl(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e, Float_t met_et, Float_t met_phi) {

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

float deltaPhi_metll(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e, Float_t met_et, Float_t met_phi) {

  if(pt.size() < 2){
    return -999;
  }
  
  TLorentzVector p1;
  TLorentzVector p2;
  TLorentzVector dil;
  TLorentzVector met;
  p1.SetPtEtaPhiM(pt[0], eta[0], phi[0], e[0]);
  p2.SetPtEtaPhiM(pt[1], eta[1], phi[1], e[1]);

  dil = (TLorentzVector)(p1+p2);
  met.SetPtEtaPhiM(met_et, 0.0, met_phi, 0.0);
    
  return dil.DeltaPhi(met);

}

bool checkPt(VecF_t& pt, float cut1, float cut2){
    if((pt[0] > cut1 && pt[1] > cut2) || (pt[1] > cut1 && pt[0] > cut2))return kTRUE;
    return kFALSE;
}

float getVar(Vec2_t& var, int idx){
    const auto size = var.size();
    if(idx > size){
    printf("Can not ask for idx %i when there are only %i object(s)", idx, size);
    return -1;       
  }
    return var[idx];  
}

int getTypeTimesCharge(Vec_t& chlep, Vec_t& fllep, int idx){
    const auto sizech = chlep.size();
    const auto sizefl = fllep.size();
    
    if(idx > sizech or idx > sizefl){
    printf("Can not ask for idx %i when there are only %i object(s)", idx, sizech);

    return -1;       
  }
    if(sizech != sizefl){
    printf("Different size of fl %i wand charge %i.", sizefl, sizech);

    return -1;       
  }
    return chlep[idx]*fllep[idx];  
}

float getE(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& m, int idx){
    TLorentzVector lep;
    const auto size = pt.size();
    if(idx > size){
        printf("Can not ask for idx %i when there are only %i object(s)", idx, size);
        return -1;
    }
    lep.SetPtEtaPhiM(pt[idx], eta[idx], phi[idx], m[idx]);
    return lep.E();
}

float getM(VecF_t& m, int idx){
    const auto size = m.size();
    if(idx > size){
        printf("Can not ask for idx %i when there are only %i object(s)", idx, size);
        return -1;
    }
    float m_scaled = m[idx];
    return m_scaled;
}

float getMt(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& m, int idx){
    TLorentzVector lep;
    const auto size = pt.size();
    if(idx > size){
        printf("Can not ask for idx %i when there are only %i object(s)", idx, size);
        return -1;
    }
    lep.SetPtEtaPhiM(pt[idx], eta[idx], phi[idx], m[idx]);
    return lep.Mt();
}


float deltaR(VecF_t& eta, VecF_t& phi, int idx1, int idx2){
    float dEta = eta[idx1] - eta[idx2];
    float dPhi = phi[idx1] - phi[idx2];
    float dR = sqrt(dEta*dEta + dPhi*dPhi);
    return dR;
             
}

float SSHt(VecF_t& pt, VecF_t& ch){
    int idx1;
    int idx2;
    if (ch[0] == ch[1]){
        idx1 = 0;
        idx2 = 1;
    }
    else if (ch[0] == ch[2]){
        idx1 = 0;
        idx2 = 2;
    }
    else if (ch[1] == ch[2]){
        idx1 = 1;
        idx2 = 2;
    }
    else{
         return 0;
    }
            
    float H = pt[idx1] + pt[idx2];
    return H;
}

float S_Et_miss(Float_t et_miss, VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& m){
    TLorentzVector et1;
    TLorentzVector et2;
    TLorentzVector et3;
    et1.SetPtEtaPhiM(pt[0], eta[0], phi[0], m[0]);
    et2.SetPtEtaPhiM(pt[1], eta[1], phi[1], m[1]);
    et3.SetPtEtaPhiM(pt[1], eta[1], phi[1], m[1]);
    float sumEt = et1.Et() + et2.Et() + et3.Et();

    float sig = et_miss/sumEt;
    return sig;
}



