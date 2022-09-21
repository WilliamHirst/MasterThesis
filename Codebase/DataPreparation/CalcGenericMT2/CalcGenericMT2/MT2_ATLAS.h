#include "MT2.h"
#include "TLorentzVector.h"
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <iostream>

using namespace std;

#include "xAODBase/IParticle.h"
#include "xAODMissingET/MissingET.h"
#ifndef CALCGENERIC_MT2_ATLAS_H
#define CALCGENERIC_MT2_ATLAS_H

/*
===========================================================================
This is a class for computing (asymmetric) MT2 using the bisection method. 
Created: October 17, 2014
Author: Christopher Lester
Maintainers: Christopher Lester (lester@hep.phy.cam.ac.uk) and Benjamin Nachman (bnachman@cern.ch)
===========================================================================

===============
Instructions:
===============

For TLorentzVectors visa, visb, metp (the visible momenta and the vector sum of the missed 
particles) an for trial masses ma and mb, 

ComputeaMT2 mt2class = ComputeaMT2(visa,visb,metp,0.,80.);
double mt2 = mt2class.Compute();

you can also do something like

ComputeaMT2 mt2class = ComputeaMT2(visa,visb,metp,"Wboson","Neutrino");
double mt2 = mt2class.Compute();

to use the W boson and neutrino masses as trial masses.  If you use a name that
is not in the database, it will list them for you and quit the program.

If you want to know the actual values of the missing particle momenta that minimize MT2, you can do
std::pair <double,double> sols = mt2class.get_solutions();
double solution1 = sols.first;
double solution2 = sols.first;

*/

class ComputeaMT2 {
 public:
  double particletype(TString part);
  void listparticles();
  ComputeaMT2(TLorentzVector visa, TLorentzVector visb, TLorentzVector MET, double ma=0., double mb=0.);
  ComputeaMT2(TLorentzVector visa, TLorentzVector visb, TLorentzVector MET, TString name1, TString name2);
  ComputeaMT2(vector<TLorentzVector> visa, vector<TLorentzVector> visb, vector<TLorentzVector> MET, double ma=0., double mb=0.);
  ComputeaMT2(const xAOD::IParticle* visa, const xAOD::IParticle* visb, const xAOD::MissingET* met, double ma=0., double mb=0.);
  ComputeaMT2(const xAOD::IParticle* visa, const xAOD::IParticle* visb, TLorentzVector MET, double ma=0., double mb=0.);

  ~ComputeaMT2();
  double Compute() ;
  std::pair <double,double> get_solutions(){
    if (myMT2 < 0){
      std::cout << "MT2 Error: You have not yet run MT2!" << std::endl;
      exit(1);   
    } 
    else{
      std::pair <double,double> out =  ben_findsols(myMT2, V1.Px(), V1.Py(), V1.M(), Ma, V2.Px(), V2.Py(), ET.Px(), ET.Py(), V2.M(),  Mb);
      return out;
    }
  }
 private:
  void Init(TLorentzVector visa, TLorentzVector visb, TLorentzVector MET, double ma, double mb);
  TLorentzVector V1;
  TLorentzVector V2;
  TLorentzVector ET;
  double Ma;
  double Mb;
  double myMT2;
  double sol1;
  double sol2; 
};

ComputeaMT2::ComputeaMT2(TLorentzVector visa, TLorentzVector visb, TLorentzVector MET, double ma, double mb){  
  Init(visa, visb, MET, ma, mb);
}
ComputeaMT2::ComputeaMT2(TLorentzVector visa, TLorentzVector visb, TLorentzVector MET, TString name1, TString name2){
  if (particletype(name1) >= 0 && particletype(name2) >= 0){
    Init(visa, visb, MET, particletype(name1), particletype(name2));
  }
  else{
    std::cout << "MT2 Error: the particle names you give are not on the list (or have the wrong syntax)!" << std::endl;
    std::cout << "Choices are: " << std::endl;
    listparticles();
    exit(1);
  }
}

ComputeaMT2::ComputeaMT2(const xAOD::IParticle* visa, const xAOD::IParticle* visb, TLorentzVector MET, double ma, double mb){  
  Init(visa->p4(), visb->p4(), MET, ma, mb);
}
ComputeaMT2::ComputeaMT2(vector<TLorentzVector> visa, vector<TLorentzVector> visb, vector<TLorentzVector> MET, double ma, double mb){
  TLorentzVector visa_sum;
  TLorentzVector visb_sum;
  TLorentzVector met_sum;
  for (unsigned int i=0; i<visa.size(); i++){
    visa_sum+=visa[i];
  }
  for (unsigned int i=0; i<visb.size();i++){
    visb_sum+=visb[i];
  }
  for (unsigned int i=0; i<MET.size();i++){
    met_sum+=MET[i];
  }
  Init(visa_sum, visb_sum, met_sum, ma, mb);
}
ComputeaMT2::ComputeaMT2(const xAOD::IParticle* visa, const xAOD::IParticle* visb, const xAOD::MissingET* met, double ma, double mb){
  TLorentzVector MET;
  MET.SetXYZM(met->mpx()/1000.,met->mpy()/1000.,0.,0.);
  Init(visa->p4(), visb->p4(), MET, ma, mb);
}

void ComputeaMT2::Init(TLorentzVector visa, TLorentzVector visb, TLorentzVector MET, double ma, double mb) {
  
  if (V1.M() < 0 || V2.M() < 0){
    std::cout << "MT2 Error: Your input particles have negative mass" << std::endl;
    exit(1);
  }

  V1 = visa;
  V2 = visb;
  ET = MET;
  Ma = ma;
  Mb = mb;
  myMT2 = -1;
  sol1 = -1;
  sol2 = -1;

}

ComputeaMT2::~ComputeaMT2() {
}

double ComputeaMT2::Compute()  {

  asymm_mt2_lester_bisect::disableCopyrightMessage();
  myMT2 = asymm_mt2_lester_bisect::get_mT2(V1.M(),V1.Px(),V1.Py(),V2.M(),V2.Px(),V2.Py(),ET.Px(),ET.Py(),Ma,Mb);
  sol1 = 0;
  sol2 = 0;
  return myMT2;

}

double ComputeaMT2::particletype(TString part){

  if (part == "Wboson_PDG2014_GeV") return 80.385;
  if (part == "Neutrino") return 0;
  if (part == "Top_PDG2014_GeV") return 173.21;
  if (part == "Lepton") return 0.;
  return -1;

}

void ComputeaMT2::listparticles(){
  std::cout << "Wboson_PDG2014_GeV, mass = 80.385 GeV" << std::endl;
  std::cout << "Neutrino, mass = 0 GeV" << std::endl;
  std::cout << "Top_PDG2014_GeV, mass = 173.21 GeV" << std::endl;
  std::cout << "Lepton, mass = 0 GeV" << std::endl;
}
#endif


