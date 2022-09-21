#ifndef GetFakeWeight_h
#define GetFakeWeight_h
#include <TH1.h>
#include <TH2.h>


typedef map< TString, TH1F* > th1Map;
typedef map< TString, TH2F* > th2Map;


class GetFakeWeight

{
  
  
 public:

  GetFakeWeight() {}
  ~GetFakeWeight() {}

  enum uncertainty
    {
     NOM = 0,
     STATUP = 1,
     STATDW = 2,
     WEIGHTUP = 3,
     WEIGHTDW = 4,
     XSECUP = 5,
     XSECDW = 6,
     SYSTUP = 7,
     SYSTDW = 8,
     TRIGUP = 9,
     TRIGDW = 10,
     WEIGHTSYSTUP = 11,
     WEIGHTSYSTDW = 12,
    };

  static constexpr std::initializer_list<uncertainty> all_unc = {NOM,
							  STATUP,
							  STATDW,
							  WEIGHTUP,
							  WEIGHTDW,
							  XSECUP,
							  XSECDW,
							  SYSTUP,
							  SYSTDW,
							  TRIGUP,
							  TRIGDW,
							  WEIGHTSYSTUP,
							  WEIGHTSYSTDW};
  std::vector<TString> uncertainty_keys = {"NOM","STATUP","STATDW","WEIGHTUP","WEIGHTDW","XSECUP","XSECDW","EFF_TOTAL_1up","EFF_TOTAL_1down","EFF_TRIG_TOTAL_1up","EFF_TRIG_TOTAL_1down","WEIGHTSYSTUP","WEIGHTSYSTDW"};
  
  
  TString this_year;
  TString this_ana;
  TString this_typ;

  TString frac_var;
  
  map<TString, Double_t> upper_val_pt;
  map<TString, Double_t> upper_val_eta;
  
  double hf_frac_nom = 0.;
  double lf_frac_nom = 0.;
  double co_frac_nom = 0.;

  double hf_frac_max = 0.;
  double lf_frac_max = 0.;
  double co_frac_max = 0.;

  double hf_frac_min = 0.;
  double lf_frac_min = 0.;
  double co_frac_min = 0.;

  double hf_frac_statdw = 0;
  double lf_frac_statdw = 0;
  double co_frac_statdw = 0;
                    		  
  double hf_frac_statup = 0;
  double lf_frac_statup = 0;
  double co_frac_statup = 0;


  th1Map h_el_frac_hf;
  th1Map h_el_frac_lf;
  th1Map h_el_frac_co;
  th1Map h_el_frac_cf;
  th1Map h_el_frac_em;

  th1Map h_el_fake_hf;
  th1Map h_el_fake_hf_xsecup;
  th1Map h_el_fake_hf_xsecdw;
  th1Map h_el_fake_co;

  th1Map h_el_real_1D;
  th1Map h_mu_real_1D;
  th2Map h_el_real;
  th2Map h_mu_real;
  th2Map h_el_real_mc;
  th2Map h_mu_real_mc;

  th1Map h_mu_fake_hf;
  th1Map h_mu_fake_hf_xsecup;
  th1Map h_mu_fake_hf_xsecdw;
  th1Map h_el_fake_lf;

  int verbosity;
  
  void initializeHist(TString LTDef, TString eff_file, TString frac_file, bool doThreeLep, TString variable, bool is1516 = false, bool is17 = false, bool is18 = false, std::vector<TString> uncert_keys = {}, std::vector<std::string> bdt_vec = {});
  int setFrac(float var_dep, float mll, int nJet30, int n_bjet, double zmass, float met_Sign, TLorentzVector met_tlv, TLorentzVector lep1, TLorentzVector lep2);
  int setFrac(float var_dep, float mll, int nJet30, int n_bjet, double zmass, float met_Sign);
  void setVariables(TString year, TString ana, TString eventStr);
  double getFake(float pT, float eta, int leptype, TString triggermatch = "", uncertainty u = NOM);
  double getReal(float pT, float eta, int leptype, TString triggermatch = "", uncertainty u = NOM, TString BDTname = "", double BDTvalue = -1.0);
  void setVB(int vb = 0){verbosity = vb;}
  float getFracHF(){return hf_frac_nom;}
  float getFracLF(){return lf_frac_nom;}
  float getFracCO(){return co_frac_nom;}
  void setFracVar(TString var){frac_var = var;}
  std::pair <Double_t,Double_t> getUpperCut(TString key);
  double getFakeUnc(TString key, int leptype, double pT_hf, double pT_lf, double pT_co, th1Map h_hf, th1Map h_lf, th1Map h_co, th1Map h_hf_xsecup, th1Map h_hf_xsecdw, uncertainty u);
  double getRealUnc(TString key, int leptype, double pT, double eta, th2Map h_real, uncertainty u, th1Map h_real_1D = {}, TString BDTname = "", double BDTvalue = -1.0);
  TString getUncKey(uncertainty u){return uncertainty_keys.at(u);}
  double getWeightedFake(double hf_nom, double lf_nom, double co_nom, double hf_var, double lf_var, double co_var, th1Map h_hf, th1Map h_lf, th1Map h_co, double pT_hf, double pT_lf, double pT_co, TString key, TString syst);
};


#endif
