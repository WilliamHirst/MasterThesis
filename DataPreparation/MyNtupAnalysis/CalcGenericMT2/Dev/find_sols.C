#include "../CalcGenericMT2/MT2_ROOT.h"
#include <iostream>
#include <cmath>
#include <cstdlib> // for exit()
#include <cassert>

#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TH1F.h"
#include "TF2.h"
#include "THStack.h"
#include "TLegend.h"
#include "TCanvas.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include "TLine.h"
#include "TStyle.h"

double pi = 4.*atan(1.);

double fAT(TLorentzVector, TLorentzVector);
vector<double> CA(TLorentzVector visa, TLorentzVector visb, TLorentzVector MET, double ma, double mb);
TLorentzVector perpify(TLorentzVector, TLorentzVector);
vector<TLorentzVector>  genBranch(double mT, double mW, double mB, double mN, double mL);

double MT(TLorentzVector pa, TLorentzVector pb){
  TLorentzVector pap;
  TLorentzVector pbp;
  pap.SetPtEtaPhiM(pa.Pt(),0.,pa.Phi(),pa.M());
  pbp.SetPtEtaPhiM(pb.Pt(),0.,pb.Phi(),pb.M());
  return (pap+pbp).M();
}

void MT2(){

  TH1F *diff_analytic = new TH1F("MT2","MT2",300,-1e-6,1e-6);

  double mT = 172;
  double mW = 80;
  double mN = 0;
  double mL = 0.0001;
  double mB = 10.0;
  int number = 10;

  for (int i = 1; i<= number; i++){

    if (i%100==0) std::cout << i << std::endl;

    vector<TLorentzVector> side1 = genBranch(mT, mW, mB, mB, mB);
    vector<TLorentzVector> side2 =  genBranch(mT, mW, mB, mN, mL);

    //side1
    TLorentzVector p4B1 = side1[0];
    TLorentzVector p4W1 = side1[1];
    TLorentzVector p4L1 = side1[2];
    TLorentzVector p4N1 = side1[3];

    //side2
    TLorentzVector p4B2 = side2[0];
    TLorentzVector p4W2 = side2[1];
    TLorentzVector p4L2 = side2[2];
    TLorentzVector p4N2 = side2[3];

    TLorentzVector MET = p4N2+p4W1;

    TLorentzVector up = -MET + p4B1 + p4B2 + p4L2;
    TLorentzVector ep = perpify(p4L2,up);

    TLorentzVector j1p = perpify(p4B1,up);
    TLorentzVector j2p = perpify(p4B2,up);
    TLorentzVector metp = perpify(MET,up);

    TLorentzVector visa = j2p+ep;
    TLorentzVector visb = j1p;

    double Ma = 0;
    double Mb = 80.;

    vector<double> mt2a_ana = CA(visa,visb,metp,Ma,Mb);
    ComputeMT2 stuff = ComputeMT2(visa,visb,metp,Ma,Mb);
    double lester = stuff.Compute();
    if (mt2a_ana[0]==1) diff_analytic->Fill(lester-mt2a_ana[1]);
    std::cout << lester << " " << mt2a_ana[1] << " " << mt2a_ana[0] << std::endl;

    double mymin = 10000000;
    double xx = 0.;
    double yy = 0.;
    for (double metpx = -100; metpx<=100; metpx+=0.5){
      for (double metpy = -100; metpy<=100; metpy+=0.5){
        TLorentzVector hold = TLorentzVector(metpx,metpy,0.,sqrt(Ma*Ma+metpx*metpx+metpy*metpy));
        TLorentzVector hold2 = TLorentzVector(metp.Px()-metpx,metp.Py()-metpy,0.,sqrt(Mb*Mb+pow(metp.Px()-metpx,2)+pow(metp.Py()-metpy,2)));
        if (max(MT(visa,hold),MT(visb,hold2)) < mymin){
          mymin = max(MT(visa,hold),MT(visb,hold2));
          xx = metpx;
          yy = metpy;
          //std::cout << "here " << mymin << " " << lester << " " << metpx << " " << metpy << " " << MT(visa,hold) << " " << MT(visb,hold)  << std::endl;
        }
        //std::cout << metpx << " " << metpy << " " << MT(visa,hold) << " " << MT(visb,hold2) << " " << max(MT(visa,hold),MT(visb,hold2)) << " " << lester << std::endl;
      }
    }


    double E = sqrt(visa.Pt()*visa.Pt()+visa.M()*visa.M());
    double px = visa.Px();
    double py = visa.Py();
    double M = lester;
    double E2 = E*E;
    double M2 = M*M;
    double M4 = M2*M2;
    double Ma2 = Ma*Ma;
    double Ma4 = Ma2*Ma2;
    double px2 = px*px;
    double py2 = py*py;
    double px4 = px2*px2;
    double py4 = py2*py2;
    double px3 = px2*px;
    double py3 = py2*py;
    double E4 = E2*E2;
    double TermA = E2*px-M2*px+Ma2*px-px2*px-px*py2;
    double TermB = -2.*px*py;
    double TermSqy0 = E4*E2-2.*E4*M2-2.*E4*Ma2-2.*E4*px2-2.*E4*py2+E2*M4-2.*E2*M2*Ma2+2.*E2*M2*px2+2.*E2*M2*py2+E2*Ma4+2.*E2*Ma2*px2-2.*E2*Ma2*py2+E2*px4+2.*E2*px2*py2+E2*py4;
    double TermSqy1 = -4.*E4*py+4.*E2*M2*py-4.*E2*Ma2*py+4.*E2*px2*py+4.*E2*py3;
    double TermSqy2 = -4.*E4+4.*E2*px2+4.*E2*py2;

    //First, determine the range.

    double myx = 0.;
    double myy = 0.;

    if (TermSqy1*TermSqy1-4.*TermSqy0*TermSqy2 < 0){
      //std::cout << "no solutions! are you sure you found MT2???" << std::endl; //unbalanced!
    }
    else{

      double sol1 = (-TermSqy1 - sqrt(TermSqy1*TermSqy1-4.*TermSqy0*TermSqy2))/(2.*TermSqy2);
      double sol2 = (-TermSqy1 + sqrt(TermSqy1*TermSqy1-4.*TermSqy0*TermSqy2))/(2.*TermSqy2);
      double low = sol1;
      double high = sol2;
      if (low > high){
	low = sol2;
	high = sol1;
      }

      double myclose = 99999999.;
      for (double metpy = low; metpy<=high; metpy+=(high-low)/10000.){
	double metpx = -(TermB*metpy+TermA-sqrt(TermSqy0+TermSqy1*metpy+TermSqy2*metpy*metpy))*0.5/(E2-px2);
	double metpx2 = -(TermB*metpy+TermA+sqrt(TermSqy0+TermSqy1*metpy+TermSqy2*metpy*metpy))*0.5/(E2-px2);
	TLorentzVector hold = TLorentzVector(metpx,metpy,0.,sqrt(metpx*metpx+metpy*metpy+Ma*Ma));
	TLorentzVector hold2 = TLorentzVector(metpx2,metpy,0.,sqrt(metpx2*metpx2+metpy*metpy+Ma*Ma));
	double mt1a = MT(visa,hold);
	double mt1b = MT(visa,hold2);
	double metpxb = metp.Px()-metpx;
	double metpx2b = metp.Px()-metpx2;
	TLorentzVector holdb = TLorentzVector(metpxb,metp.Py()-metpy,0.,sqrt(metpxb*metpxb+pow(metp.Py()-metpy,2)+Mb*Mb));
	TLorentzVector hold2b = TLorentzVector(metpx2b,metp.Py()-metpy,0.,sqrt(metpx2b*metpx2b+pow(metp.Py()-metpy,2)+Mb*Mb));
	double mt2a = MT(visb,holdb);
	double mt2b = MT(visb,hold2b);
	if (fabs(mt1a-mt2a) < myclose){
	  myclose = fabs(mt1a-mt2a);
	  myy = metpy;
	  myx = metpx;
	}
	if (fabs(mt1b-mt2b) < myclose){
	  myclose = fabs(mt1b-mt2b);
	  myy = metpy;
	  myx = metpx2;
	}
      }
    }

    std::pair <double,double> sols = stuff.get_solutions();
    double misspx = sols.first;
    double misspy = sols.second;

    std::cout << myx << " " << myy << " " << xx << " " << yy << " " << misspx << " " << misspy << std::endl;
  }

  TCanvas * c1 = new TCanvas("","",500,500);
  diff_analytic->GetXaxis()->SetTitle("Lester - Matchev [GeV]");
  diff_analytic->Draw();
  c1->Print("test.pdf");

}

TLorentzVector perpify(TLorentzVector V, TLorentzVector P){

  TVector3 V3 = V.Vect();
  TVector3 P3 = P.Vect();
  TVector3 hold1 = V3.Cross(P3);
  hold1 = P3.Cross(hold1);
  hold1 = hold1*(1./(P3.Mag()*P3.Mag()));
  TLorentzVector out;
  out.SetXYZM(hold1.Px(),hold1.Py(),V.Pz(),V.M());

  return out;
}

vector<TLorentzVector>  genBranch(double mT, double mW, double mB, double mN, double mL){
    //generates b-l-nu system from top quark at rest

    TRandom3 rand(0); //was 10

    //quantities in the ttbar frame
    double pW = sqrt( (pow(mT,2)-pow((mW-mB),2))*(pow(mT,2)-pow((mW+mB),2)))/(2.*mT);

    double eW = (pow(mT,2)+pow(mW,2)-pow(mB,2))/(2.*mT);
    //double eB = (pow(mT,2)+pow(mB,2)-pow(mW,2))/(2.*mT); //Correct, but not necesary
    double betaWMag = pW/eW;

    //quantities in the W rest frame
    double pL = sqrt( (pow(mW,2)-pow((mL-mN),2))*(pow(mW,2)-pow((mL+mN),2)))/(2.*mW);
    double eL = (pow(mW,2)+pow(mL,2)-pow(mN,2))/(2.*mW);
    //double eN = (pow(mW,2)-pow(mL,2)+pow(mN,2))/(2.*mW); //Correct, but not necesary
    double betatauMag = pL/eL;

    //generate W direction
    double cosThetaW = rand.Uniform(-1.,1.);
    double sinThetaW = sqrt(1.-pow(cosThetaW,2));
    double phiW = rand.Uniform(0.,2.*pi);
    TVector3 betaW = TVector3(betaWMag*sinThetaW*cos(phiW), betaWMag*sinThetaW*sin(phiW), betaWMag*cosThetaW);
    TLorentzVector p4W;
    p4W.SetXYZM(pW*sinThetaW*cos(phiW), pW*sinThetaW*sin(phiW), pW*cosThetaW, mW);

    //b quark in ttbar frame
    TLorentzVector p4B;
    p4B.SetXYZM(-pW*sinThetaW*cos(phiW), -pW*sinThetaW*sin(phiW), -pW*cosThetaW, mB);

    //generate a direction for the tau in the W rest frame
    double cosThetaL = rand.Uniform(-1.,1.);
    double sinThetaL = sqrt(1.-pow(cosThetaL,2));
    double phiL = rand.Uniform(0.,2.*pi);
    TLorentzVector p4L;
    p4L.SetXYZM(pL*sinThetaL*cos(phiL), pL*sinThetaL*sin(phiL), pL*cosThetaL, mL);
    TVector3 betatau = TVector3(betatauMag*sinThetaL*cos(phiL), betatauMag*sinThetaL*sin(phiL), betatauMag*cosThetaL);

    //neutrino in the W rest frame
    TLorentzVector p4N;
    p4N.SetXYZM(-pL*sinThetaL*cos(phiL), -pL*sinThetaL*sin(phiL), -pL*cosThetaL, mN);

    //boost lepton2 and neutrinos to the ttbar frame
    p4L.Boost(betaW);
    p4N.Boost(betaW);

    vector<TLorentzVector> out;
    out.push_back(p4B);
    out.push_back(p4W);
    out.push_back(p4L);
    out.push_back(p4N);

    return out;
}

vector<double> CA(TLorentzVector visa, TLorentzVector visb, TLorentzVector MET, double ma, double mb){

  vector<double> out;

  TLorentzVector V1 = visa;
  TLorentzVector V2 = visb;
  TLorentzVector ET = MET;
  double Ma = ma;
  double Mb = mb;

  TLorentzVector a = V1;
  TLorentzVector b = V2;

  double M2p = 0.5*(Mb*Mb+Ma*Ma);
  double M2m = 0.5*(Mb*Mb-Ma*Ma);

  TVector3 Q = -(a.Vect()+b.Vect());
  double Q2 = Q.Px()*Q.Px()+Q.Py()*Q.Py();
  double parta = Q2+(Ma+Mb)*(Ma+Mb);
  parta = sqrt(std::max(0.,parta));

  double e12 = a.M()*a.M()+a.Pt()*a.Pt();
  double e1 = std::max(e12,0.);
  e1 = sqrt(e1);

  double e22 = b.M()*b.M()+b.Pt()*b.Pt();
  double e2 = std::max(e22,0.);
  e2 = sqrt(e2);

  //First, we check balanced versus unbalanced.

  TVector3 qt0a = (Ma/a.M())*a.Vect();
  TVector3 qt0b = (Mb/b.M())*b.Vect();

  TVector3 qtb = -qt0a+Q;
  TVector3 qta = -qt0b+Q;

  double eatild2 = Ma*Ma+qta.X()*qta.X()+qta.Y()*qta.Y();
  double eatild = std::max(eatild2,0.);
  eatild = sqrt(eatild);

  double ebtild2 = Mb*Mb+qtb.X()*qtb.X()+qtb.Y()*qtb.Y();
  double ebtild = std::max(ebtild2,0.);
  ebtild = sqrt(std::max(0.,ebtild));

  double MTa2 = a.M()*a.M()+Ma*Ma+2.*(e1*eatild-a.Px()*qta.Px()-a.Py()*qta.Py());
  double MTb2 = b.M()*b.M()+Mb*Mb+2.*(e2*ebtild-b.Px()*qtb.Px()-b.Py()*qtb.Py());

  double MTa = sqrt(std::max(0.,MTa2));
  double MTb = sqrt(std::max(0.,MTb2));

  if ((MTa > b.M()+Mb) && (MTb > a.M()+Ma)){ //balanced

    double AT = fAT(a,b);
    double ma2 = a.M2();
    double mb2 = b.M2();
    double commondenom = 2.*AT-ma2-mb2;

    double sqrttermA = 4.*M2p/commondenom;
    double sqrttermB = 2.*M2m/commondenom;
    sqrttermB = sqrttermB*sqrttermB;
    double sqrtterm = 1.+sqrttermA+sqrttermB;
    sqrtterm = sqrt(std::max(0.,sqrtterm));

    double sqrtfact = sqrt(std::max(0.,AT*AT-ma2*mb2));

    //Let's figure out the ambiguity in sign.

    double sqrtsminhat = e1+e2+parta; //sqrt(\hat{s}_T)_{min}

    double sterma = 2.*(e2 - e1)*M2m/commondenom;
    double stermb = ((e2+e1)*AT-(e2*ma2+e1*mb2))/sqrtfact;

    double shatM = e1+e2+sterma-stermb*sqrtterm;

    double sign = -1.;

    if (shatM < sqrtsminhat){
      sign = 1.;
    }

    //Now, back to MT2

    double MT2term1 = M2m*(mb2-ma2)/commondenom;
    double MT2sq = M2p+AT+MT2term1+sign*sqrtterm*sqrtfact;
    double MT2 = sqrt(std::max(0.,MT2sq));

    out.push_back(1);
    out.push_back(MT2);
    return out;

  }

  if ((MTa > b.M()+Mb)){ //unbalanced

    out.push_back(-1);
    out.push_back(a.M()+Ma);
    return out;
  }

  if ((MTb > a.M()+Ma)){ //unbalanced

    out.push_back(-1);
    out.push_back(b.M()+Mb);
    return out;

  }

  out.push_back(-1);
  return out;

}

double fAT(TLorentzVector vis1, TLorentzVector vis2){

  double e12 = vis1.M2()+vis1.Pt()*vis1.Pt();
  double e1 = std::max(e12,0.);
  e1 = sqrt(e1);

  double e22 = vis2.M2()+vis2.Pt()*vis2.Pt();
  double e2 = std::max(e22,0.);
  e2 = sqrt(e2);

  double dotprod = vis1.Px()*vis2.Px()+vis1.Py()*vis2.Py();

  return e1*e2+dotprod;

}



