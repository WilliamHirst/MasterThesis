import ROOT as R
import pandas as pd
from os import listdir
from os.path import isfile, join, isdir
import time
import sys
from samples import configure_samples
import plottingTool as pt
from pyHelperFunctions import *



d_samp,d_type,d_reg = configure_samples()

R.EnableImplicitMT(200)

R.gROOT.ProcessLine(".L helperFunctions.cxx+");
R.gSystem.AddDynamicPath("-I/home/wohirst/MasterThesis/Codebase/DataPreparation/RDataFrameAna")
R.gInterpreter.Declare('#include "helperFunctions.h"') # Header with the definition of the myFilter function
R.gSystem.Load("helperFunctions_cxx.so") # Library with the myFilter function

storage = "/storage/William_Sakarias/William_Data"


hname = "MET_2L_mm"


fldic = {"eee":0,
         "eem":1,
         "emm":2,
         "mem":3,
         "mmm":4,
         "mme":5,
         "mee":6,
         "eme":7,
}


bkgdic = {"Wjets":{"color":R.kMagenta},
          "Zjets2":{"color":R.kBlue-7},
          "diboson2L":{"color":R.kRed-7},
          "diboson3L":{"color":R.kBlue-7},
          "diboson4L":{"color":R.kGreen-2},
          "higgs":{"color":R.kYellow+2},
          "singletop":{"color":R.kOrange-2},
          "topOther":{"color":R.kSpring-9},
          "triboson":{"color":R.kViolet-7},
          "ttbar":{"color":R.kRed+7},
          "data18":{"color":R.kBlack}
}
featdic = {"lep1_Pt"  : {"xlabel":"P_{t}(l_{1}) [GeV]",
                        "nr_bins": 40, "min" : 25, "max" : 300},
           "lep2_Pt"  : {"xlabel":"P_{t}(l_{2}) [GeV]",
                        "nr_bins": 30, "min" : 20, "max" : 250},
           "lep3_Pt"  : {"xlabel":"P_{t}(l_{3}) [GeV]",
                        "nr_bins": 30, "min" : 7, "max" : 120},
           "lep1_E"  : {"xlabel":"E(l_{1}) [GeV]",
                        "nr_bins": 40, "min" : 20, "max" : 500},
           "lep2_E"  : {"xlabel":"E(l_{2}) [GeV]",
                        "nr_bins": 40, "min" : 20, "max" : 500},
           "lep3_E"  : {"xlabel":"E(l_{3}) [GeV]",
                        "nr_bins": 40, "min" : 20, "max" : 500},
           "lep1_Mt"  : {"xlabel":"M_{T}(l_{1}) [GeV]",
                        "nr_bins": 40, "min" : 20, "max" : 250},
           "lep2_Mt"  : {"xlabel":"M_{T}(l_{2}) [GeV]",
                        "nr_bins": 30, "min" : 15, "max" : 200},
           "lep3_Mt"  : {"xlabel":"M_{T}(l_{3}) [GeV]",
                        "nr_bins": 30, "min" : 10, "max" : 120},
           "lep1_Eta" : {"xlabel":"\eta(l_{1})",
                        "nr_bins": 20, "min" : -3, "max" : 3},
           "lep2_Eta" : {"xlabel":"\eta(l_{2})",
                        "nr_bins": 20, "min" : -3, "max" : 3},
           "lep3_Eta" : {"xlabel":"\eta(l_{3})",
                        "nr_bins": 20, "min" : -3, "max" : 3},
           "lep1_Phi" : {"xlabel":"\phi(l_{1})",
                        "nr_bins": 20, "min" : -3.5, "max" : 3.5},
           "lep2_Phi" : {"xlabel":"\phi(l_{2})",
                        "nr_bins": 20, "min" : -3.5, "max" : 3.5},
           "lep3_Phi" : {"xlabel":"\phi(l_{3})",
                        "nr_bins": 20, "min" : -3.5, "max" : 3.5},
           "lep1_Charge" : {"xlabel": "Charge(l_{1})",
                        "nr_bins": 2, "min" : -2, "max" : 2},
           "lep2_Charge" : {"xlabel": "Charge(l_{2})",
                        "nr_bins": 2, "min" : -2, "max" : 2},
           "lep3_Charge" : {"xlabel": "Charge(l_{3})",
                        "nr_bins": 2, "min" : -2, "max" : 2},
           "lep1_Flavor" : {"xlabel": "Flavor(l_{1})",
                        "nr_bins": 2, "min" : 0.5, "max" : 2.5},
           "lep2_Flavor" : {"xlabel": "Flavor(l_{2})",
                        "nr_bins": 2, "min" : 0.5, "max" : 2.5},
           "lep3_Flavor" : {"xlabel": "Flavor(l_{3})",
                        "nr_bins": 2, "min" : 0.5, "max" : 2.5},
           "met_Et"   : {"xlabel":"E_{T}^{miss}[GeV]"},
           "met_Phi"  : {"xlabel":"\phi (miss)"},
           "deltaR"  : {"xlabel":"\Delta R"},
           "mlll"  : {"xlabel":"M_{lll} [GeV]"},
           "mll_OSSF"  : {"xlabel":"M_{ll} (OSSF) [GeV]"},
           "Ht_lll"  : {"xlabel":"H_{t}(lll)[GeV]"},
           "Ht_SS"  : {"xlabel":"H_{t}(SS)[GeV]"},
           "Ht_met_Et"  : {"xlabel":"H_{t}(lll) + E_{T}^{miss}[GeV]"},
           "M_jj"  : {"xlabel":"M_{jj}[GeV]"},
           "met_Sign" : {"xlabel": "S(E_{t}^{miss}) [GeV])"}

}
Nlep = 3
lepv = ["lepPt","lepEta","lepPhi", "lepCharge", "lepFlavor"]



def getTriggerThreshold(tname):
    thr = []
    reg = re.findall(r'_\d*([e]*[mu]*\d{1,})_{0,}',tname)
    for r in reg:
        thr.append(int(re.sub('\D', '', r)))
    return max(thr)



trgdic = {"2015":{"1L":["HLT_e24_lhmedium_L1EM20VH",
                        "HLT_e60_lhmedium",
                        "HLT_e120_lhloose",
                        "HLT_mu20_iloose_L1MU15",
                        "HLT_mu50"],
                  "2L":["HLT_2e12_lhloose_L12EM10VH",
                        "HLT_2mu10",
                        "HLT_mu18_mu8noL1",
                        "HLT_e17_lhloose_mu14",
                        "HLT_e7_lhmedium_mu24"
                  ],
                  "3L":["HLT_e17_lhloose_2e9_lhloose",
                        "HLT_mu18_2mu4noL1",
                        "HLT_2e12_lhloose_mu10",
                        "HLT_e12_lhloose_2mu10"
                  ]},   
          "2016":{"1L":["HLT_e24_lhmedium_nod0_L1EM20VH",
                        "HLT_e24_lhtight_nod0_ivarloose",
                        "HLT_e26_lhtight_nod0_ivarloose",
                        "HLT_e60_lhmedium_nod0",
                        "HLT_e140_lhloose_nod0",
                        "HLT_mu26_ivarmedium",
                        "HLT_mu50"],
                  "2L":["HLT_2e15_lhvloose_nod0_L12EM13VH",
                        "HLT_2e17_lhvloose_nod0",
                        "HLT_2mu10",
                        "HLT_2mu14",
                        "HLT_mu20_mu8noL1",
                        "HLT_mu22_mu8noL1",
                        "HLT_e17_lhloose_nod0_mu14",
                        "HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1",
                        "HLT_e7_lhmedium_nod0_mu24"
                  ],
                  "3L":["HLT_e24_lhvloose_nod0_2e12_lhvloose_nod0_L1EM20VH_3EM10VH",
                        "HLT_e12_lhloose_nod0_2mu10",
                        "HLT_2e12_lhloose_nod0_mu10",
                        "HLT_mu20_2mu4noL1",
                        "HLT_3mu6",
                        "HLT_3mu6_msonly",
                        "HLT_e17_lhloose_nod0_2e10_lhloose_nod0_L1EM15VH_3EM8VH"
                  ]},
          "2017":{"1L":["HLT_e26_lhtight_nod0_ivarloose",
                        "HLT_e60_lhmedium_nod0",  
                        "HLT_e140_lhloose_nod0",    
                        "HLT_e300_etcut",                                
                        "HLT_mu26_ivarmedium",	     
                        "HLT_mu50"]
                  ,"2L":["HLT_2e17_lhvloose_nod0_L12EM15VHI",
                         "HLT_2e24_lhvloose_nod0",
                         "HLT_2mu14",
                         "HLT_mu22_mu8noL1",
                         "HLT_e17_lhloose_nod0_mu14",
                         "HLT_e26_lhmedium_nod0_mu8noL1",
                         "HLT_e7_lhmedium_nod0_mu24"
                  ],
                  "3L":["HLT_e24_lhvloose_nod0_2e12_lhvloose_nod0_L1EM20VH_3EM10VH",
                        "HLT_e12_lhloose_nod0_2mu10",
                        "HLT_2e12_lhloose_nod0_mu10",
                        "HLT_mu20_2mu4noL1",
                        "HLT_3mu6",
                        "HLT_3mu6_msonly"
                  ]},
          "2018":{"1L":["HLT_e26_lhtight_nod0_ivarloose",
                        "HLT_e60_lhmedium_nod0",  
                        "HLT_e140_lhloose_nod0",    
                        "HLT_e300_etcut",                                
                        "HLT_mu26_ivarmedium",	     
                        "HLT_mu50"],
                  "2L":["HLT_2e17_lhvloose_nod0_L12EM15VHI",
                        "HLT_2e24_lhvloose_nod0",
                        "HLT_2mu14",
                        "HLT_mu22_mu8noL1",
                        "HLT_e17_lhloose_nod0_mu14",
                        "HLT_e26_lhmedium_nod0_mu8noL1",
                        "HLT_e7_lhmedium_nod0_mu24"],
                  "3L":["HLT_e24_lhvloose_nod0_2e12_lhvloose_nod0_L1EM20VH_3EM10VH",
                        "HLT_e12_lhloose_nod0_2mu10",
                        "HLT_2e12_lhloose_nod0_mu10",
                        "HLT_mu20_2mu4noL1",
                        "HLT_3mu6"
                  ]},
}


import re
trigstr = {}
evtrigstr = {}
for yr in trgdic.keys():
    for x in trgdic[yr].keys():
        if not len(trgdic[yr][x]): continue
        if not x in trigstr.keys():
            trigstr[x] = {}
            evtrigstr[x] = {}
        if not yr in trigstr[x].keys():
            trigstr[x][yr] = "("
            evtrigstr[x][yr] = "("
        for trigger in trgdic[yr][x]:
            if trigger == "1":
                trigstr[x][yr] += "(1) || "
                evtrigstr[x][yr] += "1 || "
            else:
                trigstr[x][yr] += "(lep%s && lepPt > %i) || "%(trigger,getTriggerThreshold(trigger))
                evtrigstr[x][yr] += "trigMatch_%s || "%(trigger)
        trigstr[x][yr] = trigstr[x][yr][:-4]+")"
        evtrigstr[x][yr] = evtrigstr[x][yr][:-4]+")"




good_runs = []
    
histo = {}
allhisto = []
nEvents = 618282964
nSlots = R.GetThreadPoolSize();
print("Number of slots = %i"%nSlots);
everyN = int(100 * nSlots)

def runANA(mypath_mc, mypath_data, everyN, fldic, histo, allhisto, nEvents = 0):
    nh = 100
    if not isfile("histograms.root"):
        histo = getHistograms("histograms.root")
        return
    else:
        if isdir(mypath_mc):
            df_mc = getDataFrames(mypath_mc)
            print("Loading %s into dataframe with keys %s" %(mypath_mc,",".join(df_mc.keys())))
        else:
            df_mc = {}

        if isdir(mypath_data):
            df_data = getDataFrames(mypath_data)
            print("Loading %s into dataframe with keys %s" %(mypath_data,",".join(df_data.keys())))
        else:
            df_data = {}

        df = {**df_mc,**df_data}
        for k in df.keys():
            isData = "data" in k

            if not isData:
                df[k] = df[k].Define("scaletolumi","(RandomRunNumber) < 320000 ? 36207.65 : (((RandomRunNumber) > 320000 && (RandomRunNumber) < 348000) ? 44307.4 : 58450.1)")

            df[k] = df[k].Define("new_xsec","(DatasetNumber == 308981) ? (0.30649*69.594)/80000. : 0.0")

            # Baseline leptons
            ele_filter = "lepFlavor==1 && lepPassOR > 0 && (lepEta <= 2.47 && lepEta >= -2.47) && ((lepZ0SinTheta)<=0.5 && (lepZ0SinTheta)>=-0.5)"
            muo_filter = "lepFlavor==2 && lepPassOR > 0 && (lepEta <= 2.7  && lepEta >= -2.7) && lepLoose > 0 && ((lepZ0SinTheta)<=0.5 && (lepZ0SinTheta)>=-0.5)"

            df[k] = df[k].Define("ele_BL", ele_filter) 
            df[k] = df[k].Define("muo_BL", muo_filter) 
            df[k] = df[k].Define("nlep_BL","ROOT::VecOps::Sum(ele_BL)+ROOT::VecOps::Sum(muo_BL)")

            # Signal leptons
            df[k] = df[k].Define("ele_SG","ele_BL && lepIsoLoose_VarRad && lepTight && (lepD0Sig <= 5 && lepD0Sig >= -5)") 
            df[k] = df[k].Define("muo_SG","muo_BL && lepIsoLoose_VarRad && (lepD0Sig <= 3 && lepD0Sig >= -3)")
            df[k] = df[k].Define("nlep_SG","ROOT::VecOps::Sum(ele_SG)+ROOT::VecOps::Sum(muo_SG)")
            df[k] = df[k].Filter("nlep_BL == 3")
            df[k] = df[k].Filter("nlep_SG == 3")

            isGoodLepton = f"ele_SG || muo_SG"
            df[k] = df[k].Define("isGoodLepton", isGoodLepton)
            
            if not isData:

                df[k] = df[k].Define("is2015","RandomRunNumber <= 284500")
                df[k] = df[k].Define("is2016","(RandomRunNumber > 284500 && RandomRunNumber < 320000)")
                df[k] = df[k].Define("is2017","(RandomRunNumber > 320000 && RandomRunNumber < 348000)")
                df[k] = df[k].Define("is2018","RandomRunNumber > 348000")

                df[k] = df[k].Define("lepwgt_SG","getSF(lepRecoSF[ele_SG || muo_SG])")

                df[k] = df[k].Define("trgwgt_SG","getSF(lepTrigSF[ele_SG || muo_SG])")

                df[k] = df[k].Define("wgt_SG","(new_xsec ? (new_xsec) : (genWeight))*eventWeight*jvtWeight*bTagWeight*pileupWeight*scaletolumi*lepwgt_SG*trgwgt_SG")

                df[k] = df[k].Define("wgt_EV_SG","(eventWeight*jvtWeight*bTagWeight*pileupWeight*scaletolumi*lepwgt_SG*trgwgt_SG)")

            else:
                df[k] = df[k].Define("is2015","(RunNumber >= 276262 && RunNumber <= 284484)")
                df[k] = df[k].Define("is2016","(RunNumber >= 297730 && RunNumber <= 311481)")
                df[k] = df[k].Define("is2017","(RunNumber >= 325713 && RunNumber <= 340453)")
                df[k] = df[k].Define("is2018","RunNumber >= 348885")

                df[k] = df[k].Define("wgt_SG","1.0")
                df[k] = df[k].Define("wgt_EV","1.0")
            
            if not nEvents:
                this_nEvents = int(df[k].Count().GetValue())
                nEvents += this_nEvents
                print("Loading %s with %.0f events. Now %.0f events"%(k,this_nEvents,nEvents))
            else:
                print("Loading %s"%(k))


            # Triggers            
            trigmatch_2015_2L = "(lepHLT_2e12_lhloose_L12EM10VH[isGoodLepton] && lepPt[isGoodLepton] > 12) || (lepHLT_e17_lhloose_mu14[isGoodLepton] && lepPt[isGoodLepton] > 17) || (lepHLT_mu18_mu8noL1[isGoodLepton] && lepPt[isGoodLepton] > 18)"
            triggered_2015_2L = "(trigMatch_HLT_2e12_lhloose_L12EM10VH || trigMatch_HLT_e17_lhloose_mu14 || trigMatch_HLT_mu18_mu8noL1)"
            trigmatch_2016_2L = "(lepHLT_2e17_lhvloose_nod0[isGoodLepton] && lepPt[isGoodLepton] > 17) || (lepHLT_e17_lhloose_nod0_mu14[isGoodLepton] && lepPt[isGoodLepton] > 17) || (lepHLT_mu22_mu8noL1[isGoodLepton] && lepPt[isGoodLepton] > 22)"
            triggered_2016_2L = "(trigMatch_HLT_2e17_lhvloose_nod0 || trigMatch_HLT_e17_lhloose_nod0_mu14 || trigMatch_HLT_mu22_mu8noL1)"
            trigmatch_2017_2L = "(lepHLT_2e17_lhvloose_nod0_L12EM15VHI[isGoodLepton] && lepPt[isGoodLepton] > 17) || (lepHLT_2e24_lhvloose_nod0[isGoodLepton] && lepPt[isGoodLepton] > 24) || (lepHLT_mu22_mu8noL1[isGoodLepton] && lepPt[isGoodLepton] > 22) || (lepHLT_e17_lhloose_nod0_mu14[isGoodLepton] && lepPt[isGoodLepton] > 17)"
            triggered_2017_2L = "(trigMatch_HLT_2e17_lhvloose_nod0_L12EM15VHI || trigMatch_HLT_2e24_lhvloose_nod0 || trigMatch_HLT_mu22_mu8noL1 || trigMatch_HLT_e17_lhloose_nod0_mu14)"
            trigmatch_2018_2L = "(lepHLT_2e17_lhvloose_nod0_L12EM15VHI[isGoodLepton] && lepPt[isGoodLepton] > 17) || (lepHLT_2e24_lhvloose_nod0[isGoodLepton] && lepPt[isGoodLepton] > 24) || (lepHLT_mu22_mu8noL1[isGoodLepton] && lepPt[isGoodLepton] > 22) || (lepHLT_e17_lhloose_nod0_mu14[isGoodLepton] && lepPt[isGoodLepton] > 17)"
            triggered_2018_2L = "(trigMatch_HLT_2e17_lhvloose_nod0_L12EM15VHI || trigMatch_HLT_2e24_lhvloose_nod0 || trigMatch_HLT_mu22_mu8noL1 || trigMatch_HLT_e17_lhloose_nod0_mu14)"
            
            Trig_year = {"2015" : {"match" : trigmatch_2015_2L,
                                   "triggered" : triggered_2015_2L},
                         "2016" : {"match" : trigmatch_2016_2L,
                                   "triggered" : triggered_2016_2L},
                         "2017" : {"match" : trigmatch_2017_2L,
                                   "triggered" : triggered_2017_2L},
                         "2018" : {"match" : trigmatch_2018_2L,
                                   "triggered" : triggered_2018_2L}
                        }
            
            years = list(Trig_year.keys())
            for year in years:
                y_t = Trig_year[year]
                other_years = [years[i] for i in range(4) if years[i] != year ]
                isTriggered = y_t["triggered"]
                isMatch = y_t["match"]
                notThisYear = f"is{other_years[0]} || is{other_years[1]} || is{other_years[2]}"
                df[k] = df[k].Filter(f"{isTriggered} || {notThisYear}")
                df[k] = df[k].Define(f"lep_trig_{year}", f"{isMatch} ")
                df[k] = df[k].Filter(f"ROOT::VecOps::Sum(lep_trig_{year}) >= 2 || {notThisYear}")

            for i in range(Nlep):
                df[k] = df[k].Define("lep%i_flav"%(i+1),"getTypeTimesCharge(lepCharge[isGoodLepton],lepType[isGoodLepton],%i)"%(i))
                for v in lepv:
                    if "lep" in v:
                        var = v.replace("lep","")
                    else:
                        var = v
                    bins_dic = featdic["lep%i_%s"%(i+1,var)]
                    df[k] = df[k].Define("lep%i_%s"%(i+1,var),"getVar(lep%s[isGoodLepton],%i)"%(var,i))
                    histo["lep%i_%s_%s"%(i+1,var,k)] = df[k].Histo1D(("lep%i_%s_%s"%(i+1,var,k),
                          "lep%i_%s_%s;Feature;Entries"%(i+1,var,k),bins_dic["nr_bins"],bins_dic["min"],bins_dic["max"]),
                          "lep%i_%s"%(i+1,var),"wgt_SG")
                    
            # Stransverse mass       
            df[k] = df[k].Define("MT2_12","calcMT2(lepPt[isGoodLepton], lepEta[isGoodLepton], lepPhi[isGoodLepton], lepM[isGoodLepton], met_Et, met_Phi, 0, 1)")
            df[k] = df[k].Define("MT2_13","calcMT2(lepPt[isGoodLepton], lepEta[isGoodLepton], lepPhi[isGoodLepton], lepM[isGoodLepton], met_Et, met_Phi, 0, 2)")
            df[k] = df[k].Define("MT2_23","calcMT2(lepPt[isGoodLepton], lepEta[isGoodLepton], lepPhi[isGoodLepton], lepM[isGoodLepton], met_Et, met_Phi, 1, 2)")
            
            # Energy
            df[k] = df[k].Define("lep1_E", "getE(lepPt[isGoodLepton], lepEta[isGoodLepton], lepPhi[isGoodLepton], lepM[isGoodLepton], 0)")
            df[k] = df[k].Define("lep2_E", "getE(lepPt[isGoodLepton], lepEta[isGoodLepton], lepPhi[isGoodLepton], lepM[isGoodLepton], 1)")
            df[k] = df[k].Define("lep3_E", "getE(lepPt[isGoodLepton], lepEta[isGoodLepton], lepPhi[isGoodLepton], lepM[isGoodLepton], 2)")
            
            # Transverse Mass
            df[k] = df[k].Define("lep1_Mt", "getMt(lepPt[isGoodLepton], lepEta[isGoodLepton], lepPhi[isGoodLepton], lepM[isGoodLepton], 0)")
            df[k] = df[k].Define("lep2_Mt", "getMt(lepPt[isGoodLepton], lepEta[isGoodLepton], lepPhi[isGoodLepton], lepM[isGoodLepton], 1)")
            df[k] = df[k].Define("lep3_Mt", "getMt(lepPt[isGoodLepton], lepEta[isGoodLepton], lepPhi[isGoodLepton], lepM[isGoodLepton], 2)")
                
            # Jets
            isGoodJet = "jet_BL && (jetPt > 60 || (jetPt <=60 && jetJVT <= 0.91 && jetJVT >= -0.91))"
            df[k] = df[k].Define("jet_BL","jetPt >= 20 && (jetEta <= 2.8 && jetEta >= -2.8)")
            df[k] = df[k].Define("jet_SG","jet_BL && (jetPt > 60 || (jetPt <=60 && jetJVT <= 0.91 && jetJVT >= -0.91))")
            df[k] = df[k].Define("isGoodJet",isGoodJet)
            
            df[k] = df[k].Define("bjet85","isGoodJet && jetdl1r>=0.665")
            df[k] = df[k].Define("bjet77","isGoodJet && jetdl1r>=2.195")
            
            df[k] = df[k].Define("jet_SG_pT","jetPt[isGoodJet]")
            df[k] = df[k].Define("jet_SG_eta","jetEta[isGoodJet]")
            df[k] = df[k].Define("njet_SG","ROOT::VecOps::Sum(jet_SG)")
            
            # Delta R
            df[k] = df[k].Define("deltaR","deltaR(lepEta[isGoodLepton], lepPhi[isGoodLepton], 0, 1)")
            
            # Invariant Mass (lll)
            df[k] = df[k].Define("mlll","ComputeInvariantMass(lepPt[isGoodLepton], lepEta[isGoodLepton], lepPhi[isGoodLepton], lepM[isGoodLepton])")
            
            # Invariant Mass (ll-OSSF)
            df[k] = df[k].Define("mll_OSSF","OSSFInvariantMass(lepPt[isGoodLepton], lepEta[isGoodLepton], lepPhi[isGoodLepton], lepM[isGoodLepton], lepCharge[isGoodLepton], lepFlavor[isGoodLepton])")
            
            # Ht(lll)
            df[k] = df[k].Define("Ht_lll","ROOT::VecOps::Sum(lepPt[isGoodLepton])")
            
            # Ht (SS)
            df[k] = df[k].Define("Ht_SS","SSHt(lepPt[isGoodLepton], lepCharge[isGoodLepton])")
            
            # Ht(lll) + Energy of missing transverse momentum
            df[k] = df[k].Define("Ht_met_Et","Ht_lll + met_Et")
            
            # Invariant mass of leading jet-pair.
            df[k] = df[k].Define("M_jj","getMjj(jetPt[isGoodJet], jetEta[isGoodJet], jetPhi[isGoodJet], jetM[isGoodJet], njet_SG)")

            # nBjets
            df[k] = df[k].Define("nbjet85","ROOT::VecOps::Sum(bjet85)")
            df[k] = df[k].Define("nbjet77","ROOT::VecOps::Sum(bjet77)")

         
            # Flavour combo
            df[k] = df[k].Define("flcomp","flavourComp3L(lepFlavor[ele_BL || muo_BL])")

            # Negative Weights
            df_Neg = df[k].Filter("wgt_SG < 0")

            df_Neg = df_Neg.Define("lep1_Pt_Neg",  "lep1_Pt")
            df_Neg = df_Neg.Define("lep2_Pt_Neg",  "lep2_Pt")
            df_Neg = df_Neg.Define("lep1_Phi_Neg", "lep1_Phi")
            df_Neg = df_Neg.Define("lep2_Phi_Neg", "lep2_Phi")
            df_Neg = df_Neg.Define("wgt_SG_Abs", "abs(wgt_SG)")

            
            # HISTOGRAMS
            histo["flcomp_%s"%(k)] = df[k].Histo1D(("h_%s_%s"%("flcomp",k),"h_%s_%s"%("flcomp",k),len(fldic.keys()),0,len(fldic.keys())),"flcomp","wgt_SG")
            
            #histo["MT2_12_%s"%(k)] = df[k].Histo1D(("h_%s_%s"%("MT2_12",k),"h_%s_%s;m_{T}^{2}(12) [GeV];Entries"%("MT2_12",k),40,0,400),"MT2_12","wgt_SG")
            #histo["MT2_13_%s"%(k)] = df[k].Histo1D(("h_%s_%s"%("MT2_13",k),"h_%s_%s;m_{T}^{2}(13) [GeV];Entries"%("MT2_13",k),40,0,300),"MT2_13","wgt_SG")
            #histo["MT2_23_%s"%(k)] = df[k].Histo1D(("h_%s_%s"%("MT2_23",k),"h_%s_%s;m_{T}^{2}(23) [GeV];Entries"%("MT2_23",k),40,0,250),"MT2_23","wgt_SG")
            
            #bins_dic = featdic["lep1_E"]
            #histo["lep1_E_%s"%(k)] = df[k].Histo1D(("h_%s_%s"%("lep1_E",k),"h_%s_%s;E(l1) [GeV];Entries"%("lep1_E",k),bins_dic["nr_bins"],bins_dic["min"],bins_dic["max"]),"lep1_E","wgt_SG")
            #histo["lep2_E_%s"%(k)] = df[k].Histo1D(("h_%s_%s"%("lep2_E",k),"h_%s_%s;E(l2) [GeV];Entries"%("lep2_E",k),bins_dic["nr_bins"],bins_dic["min"],bins_dic["max"]),"lep2_E","wgt_SG")
            #histo["lep3_E_%s"%(k)] = df[k].Histo1D(("h_%s_%s"%("lep3_E",k),"h_%s_%s;E(l3) [GeV];Entries"%("lep3_E",k),bins_dic["nr_bins"],bins_dic["min"],bins_dic["max"]),"lep3_E","wgt_SG")
            
            bins_dic = featdic["lep1_Mt"]
            histo["lep1_Mt_%s"%(k)] = df[k].Histo1D(("h_%s_%s"%("lep1_Mt",k),"h_%s_%s;Mt(l1) [GeV];Entries"%("lep1_Mt",k),bins_dic["nr_bins"],bins_dic["min"],bins_dic["max"]),"lep1_Mt","wgt_SG")
            histo["lep2_Mt_%s"%(k)] = df[k].Histo1D(("h_%s_%s"%("lep2_Mt",k),"h_%s_%s;Mt(l2) [GeV];Entries"%("lep2_Mt",k),bins_dic["nr_bins"],bins_dic["min"],bins_dic["max"]),"lep2_Mt","wgt_SG")
            histo["lep3_Mt_%s"%(k)] = df[k].Histo1D(("h_%s_%s"%("lep3_Mt",k),"h_%s_%s;Mt(l3) [GeV];Entries"%("lep3_Mt",k),bins_dic["nr_bins"],bins_dic["min"],bins_dic["max"]),"lep3_Mt","wgt_SG")
            
            histo["nbjet85_%s"%(k)] = df[k].Histo1D(("h_%s_%s"%("nbjet85",k),"h_%s_%s"%("nbjet85",k),5,0,4),"nbjet85","wgt_SG")        
            histo["nbjet77_%s"%(k)] = df[k].Histo1D(("h_%s_%s"%("nbjet77",k),"h_%s_%s"%("nbjet77",k),4,0,3),"nbjet77","wgt_SG")
            
            histo["met_Phi_%s"%k] = df[k].Histo1D(("h_%s_%s"%("met_Phi",k),"h_%s_%s; Phi of missing transvere momentum;Entries"%("met_Phi",k),20,-3.5,3.5),"met_Phi","wgt_SG")

            histo["met_Et_%s"%k] = df[k].Histo1D(("h_%s_%s"%("met_Et",k),"h_%s_%s; Energy of missing transverse momentum [GeV];Entries"%("met_Et",k),40,25,350),"met_Et","wgt_SG")
                        
            histo["njet_SG_%s"%k] = df[k].Histo1D(("njet_SG_%s"%k,"njet_SG_%s"%k,10,0,10),"njet_SG","wgt_SG")

            histo["deltaR_%s"%k] = df[k].Histo1D(("deltaR_%s"%k,"deltaR_%s"%k,20,0,6),"deltaR","wgt_SG")

            histo["mlll_%s"%k] = df[k].Histo1D(("mlll_%s"%k,"mlll_%s"%k,40,50,500),"mlll","wgt_SG")

            histo["mll_OSSF_%s"%k] = df[k].Histo1D(("mll_OSSF_%s"%k,"mll_OSSF_%s"%k,40,0,400),"mll_OSSF","wgt_SG")

            histo["Ht_lll_%s"%k] = df[k].Histo1D(("Ht_lll_%s"%k,"Ht_lll_%s"%k,25,40,500),"Ht_lll","wgt_SG")

            histo["Ht_SS_%s"%k] = df[k].Histo1D(("Ht_SS_%s"%k,"Ht_SS_%s"%k,40,20,400),"Ht_SS","wgt_SG")

            histo["Ht_met_Et_%s"%k] = df[k].Histo1D(("Ht_met_Et_%s"%k,"Ht_met_Et_%s"%k,40,40,600),"Ht_met_Et","wgt_SG")
            
            histo["M_jj_%s"%k] = df[k].Histo1D(("M_jj_%s"%k,"M_jj_%s"%k,25,0,800),"M_jj","wgt_SG")
            
            histo["met_Sign_%s"%k] = df[k].Histo1D(("met_Sign_%s"%k,"met_Sign_%s"%k,20,0,20),"met_Sign","wgt_SG")


            bins_dic = featdic["lep1_Pt"]
            histo["lep1_Pt_Neg_%s"%k] = df_Neg.Histo1D(("lep1_Pt_Neg_%s"%k,"lep1_Pt_Neg_%s"%k,bins_dic["nr_bins"],bins_dic["min"],bins_dic["max"]),"lep1_Pt_Neg","wgt_SG_Abs")
            bins_dic = featdic["lep2_Pt"]
            histo["lep2_Pt_Neg_%s"%k] = df_Neg.Histo1D(("lep2_Pt_Neg_%s"%k,"lep2_Pt_Neg_%s"%k,bins_dic["nr_bins"],bins_dic["min"],bins_dic["max"]),"lep2_Pt_Neg","wgt_SG_Abs")
            bins_dic = featdic["lep1_Phi"]
            histo["lep1_Phi_Neg_%s"%k] = df_Neg.Histo1D(("lep1_Phi_Neg_%s"%k,"lep1_Phi_Neg_%s"%k,bins_dic["nr_bins"],bins_dic["min"],bins_dic["max"]),"lep1_Phi_Neg","wgt_SG_Abs")
            bins_dic = featdic["lep2_Phi"]
            histo["lep2_Phi_Neg_%s"%k] = df_Neg.Histo1D(("lep2_Phi_Neg_%s"%k,"lep2_Phi_Neg_%s"%k,bins_dic["nr_bins"],bins_dic["min"],bins_dic["max"]),"lep2_Phi_Neg","wgt_SG_Abs")

        for k in histo.keys():
            allhisto.append(histo[k])


        print("Calculating %i histograms"%len(allhisto))
        start = time.time()
        R.RDF.RunGraphs(allhisto)
        end = time.time()
        print("%10i | %.2f"%(len(allhisto),(end - start)))

        hfile = R.TFile("histograms.root","RECREATE")
        hfile.cd()

        writeHistsToFile(histo, True)

        return histo, df

histo, df = runANA("/storage/shared/data/master_students/William_Sakarias/data/PHYS_3LBkgs_mc16e",
               "/storage/shared/data/master_students/William_Sakarias/data/data18",
               everyN,fldic,histo,allhisto)


if 0:
    newd = dict.fromkeys(histo)
    for key in newd.keys():
        if "_BL_" in key:
            key_SG = key.replace("_BL_","_SG_")
            key_BL = key
            key_EF = key.replace("_BL_","_EF_")
        else:
            continue
        if key_EF in histo.keys():
            continue
        print(key_SG,key_BL,key_EF)
        if key_SG in histo.keys() and key_BL in histo.keys() and not key_EF in histo.keys():
            try: 
                histo[key_EF] = getRatio1D(histo[key_SG],histo[key_BL],1)
            except:
                histo[key_EF] = getRatio1D(histo[key_SG].GetValue(),histo[key_BL].GetValue(),1)
writeHistsToFile(histo, False)

toplot = []
for bkg in bkgdic.keys():
    toplot.append(bkg)

allFeats = histo.keys()
featuresPlot = []
for feat in allFeats:
    if feat[-5:] == "Wjets":
        featuresPlot.append(feat[:-6])

for feature in featuresPlot:
    try:
        print(feature)
        if feature in featdic:
            xlabel = featdic[feature]["xlabel"]
        else:
            xlabel = feature
        p = pt.Plot(histo,feature,toplot,xtext = xlabel)
        p.can.SaveAs(f"../../../thesis/Figures/FeaturesHistograms/{feature}.pdf")
        p.can.Draw()
    except:
        print(f"Was not able to plot histogram for {feature}.")


featuresPlot.append("wgt_SG")
features = []
for feat in featuresPlot:
    if "Neg" not in feat:
        features.append(feat)

exit()
df_s = {}
for k in df.keys():
    print("Transforming " + k + "-ROOT to Numpy.")
    numpy = df[k].AsNumpy(features)
    print("Transforming " + k + "-ROOT to Pandas.")
    pandas = pd.DataFrame(data=numpy)
    del numpy
    pandas.info()
    print("Done")
    
    print("Saving to hdf5...")
    if "data" in k:
        pandas.to_hdf(f"{storage}/3L_data.hdf5","mini")
    else:
        pandas.to_hdf(f"{storage}/" + k +"mc.hdf5","mini")
    print("Done")
    del pandas

