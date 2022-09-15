import ROOT as R
import pandas as pd
from os import listdir
from os.path import isfile, join, isdir
import array
import time
import sys
from samples import configure_samples
import plottingTool as pt



d_samp,d_type,d_reg = configure_samples()#False,False,True,False,False)

R.EnableImplicitMT(200)

R.gROOT.ProcessLine(".L helperFunctions.cxx+");
R.gSystem.AddDynamicPath("-I/home/wohirst/MasterThesis/myNtupAnalysis/RDataFrameAna")
R.gInterpreter.Declare('#include "helperFunctions.h"') # Header with the definition of the myFilter function
R.gSystem.Load("helperFunctions_cxx.so") # Library with the myFilter function

from IPython.display import display, HTML
display(HTML("<style>.container { width:85% !important; }</style>"))



hname = "MET_2L_mm"


fldic = {"eee":0,
         "eem":1,
         "emm":2,
         "mmm":3,
         "mme":4,
         "mee":5,
         "all":9
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
                        "nr_bins": 60, "min" : 0, "max" : 300},
           "lep2_Pt"  : {"xlabel":"P_{t}(l_{1}) [GeV]",
                        "nr_bins": 60, "min" : 0, "max" : 300},
           "lep3_Pt"  : {"xlabel":"P_{t}(l_{1}) [GeV]",
                        "nr_bins": 60, "min" : 0, "max" : 300},
           "lep1_E"  : {"xlabel":"E(l_{1}) [GeV]",
                        "nr_bins": 100, "min" : 0, "max" : 500},
           "lep2_E"  : {"xlabel":"E(l_{2}) [GeV]",
                        "nr_bins": 100, "min" : 0, "max" : 500},
           "lep3_E"  : {"xlabel":"E(l_{3}) [GeV]",
                        "nr_bins": 100, "min" : 0, "max" : 500},
           "lep1_Mt"  : {"xlabel":"M_{T}(l_{1}) [GeV]",
                        "nr_bins": 60, "min" : 0, "max" : 300},
           "lep2_Mt"  : {"xlabel":"M_{T}(l_{2}) [GeV]",
                        "nr_bins": 60, "min" : 0, "max" : 300},
           "lep3_Mt"  : {"xlabel":"M_{T}(l_{3}) [GeV]",
                        "nr_bins": 60, "min" : 0, "max" : 200},
           "lep1_Eta" : {"xlabel":"\eta(l_{1})",
                        "nr_bins": 50, "min" : -3, "max" : 3},
           "lep2_Eta" : {"xlabel":"\eta(l_{2})",
                        "nr_bins": 50, "min" : -3, "max" : 3},
           "lep3_Eta" : {"xlabel":"\eta(l_{3})",
                        "nr_bins": 50, "min" : -3, "max" : 3},
           "lep1_Phi" : {"xlabel":"\phi(l_{1})",
                        "nr_bins": 50, "min" : -3.5, "max" : 3.5},
           "lep2_Phi" : {"xlabel":"\phi(l_{2})",
                        "nr_bins": 50, "min" : -3.5, "max" : 3.5},
           "lep3_Phi" : {"xlabel":"\phi(l_{3})",
                        "nr_bins": 50, "min" : -3.5, "max" : 3.5},
           "lep1_M"   : {"xlabel": "M(l_{1})[MeV]",
                        "nr_bins": 24, "min" : 0, "max" : 120},
           "lep2_M"   : {"xlabel": "M(l_{2})[MeV]",
                        "nr_bins": 24, "min" : 0, "max" : 120},
           "lep3_M"   : {"xlabel": "M(l_{3})[MeV]",
                        "nr_bins": 24, "min" : 0, "max" : 120},
           "lep1_Ptcone30" : {"xlabel": "P_{t}cone30(l_{1})",
                        "nr_bins": 5, "min" : 0, "max" : 50},
           "lep2_Ptcone30" : {"xlabel": "P_{t}cone30(l_{1})",
                        "nr_bins": 5, "min" : 0, "max" : 50},
           "lep3_Ptcone30" : {"xlabel": "P_{t}cone30(l_{1})",
                        "nr_bins": 5, "min" : 0, "max" : 50},
           "lep1_Z0" : {"xlabel": "Z_{0}(l_{1})",
                        "nr_bins": 100, "min" : 0, "max" : 300},
           "lep2_Z0" : {"xlabel": "Z_{0}(l_{2})",
                        "nr_bins": 100, "min" : 0, "max" : 300},
           "lep3_Z0" : {"xlabel": "Z_{0}(l_{3})",
                        "nr_bins": 100, "min" : 0, "max" : 300},
           "met_et"   : {"xlabel":"E_{T}^{miss}[GeV]"},
           "met_phi"  : {"xlabel":"\phi (miss)"},
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
lepv = ["lepPt","lepEta","lepPhi", "lepM", "lepPtcone30", "lepZ0"]



def getTriggerThreshold(tname):
    thr = []
    #print(tname)
    reg = re.findall(r'_\d*([e]*[mu]*\d{1,})_{0,}',tname)
    for r in reg:
        #print(int(re.sub('\D', '', r)))
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


def convertRDFCutflowToTex(cutflow1,cutflow2):
    i = 0
    tabstr = ""
    for c in cutflow1:
        cname = c.GetName()
        c2 = cutflow2.At(cname)
        if i == 0:
            nevc1 = c.GetAll()
            nevc2 = c2.GetAll()
        cname = cname.replace(">","$>$")
        cname = cname.replace("<","$<$")
        tabstr += "%-30s & $%.0f$ & $%.0f$ & $%.2f$ & $%.2f$ & $%.0f$ & $%.0f$ & $%.2f$ & $%.2f$ \\\ \n"%(cname,c.GetPass(),c.GetAll(),c.GetEff(),(c.GetPass()/nevc1)*100.,c2.GetPass(),c2.GetAll(),c2.GetEff(),(c2.GetPass()/nevc2)*100.)
        i += 1
    print(tabstr)


def writeHistsToFile(histo, writetofile = True):
    for k in histo.keys():
        col = -1
        sp = k.split("_")
        typ = ""
        for i in range(len(sp)):
            s = "_".join(sp[i:])
            if s in d_samp.keys():
                typ = s
        if not typ:
            print("Did to find match for key %s"%k)
            continue
        #for plk in d_samp.keys():
        #    if plk == typ:
        #print(typ)
        evtyp = list(fldic.keys())
        if "flcomp" in k:
            for i in range(1,histo[k].GetNbinsX()+1):
                histo[k].GetXaxis().SetBinLabel(i,evtyp[i-1])
        if d_samp[typ]["type"] == "bkg":
            histo[k].SetFillColor(d_samp[typ]["f_color"])
            histo[k].SetLineColor(d_samp[typ]["f_color"])
            histo[k].SetMarkerStyle(0)
            histo[k].SetMarkerSize(0)
        elif d_samp[typ]["type"] == "data":
            histo[k].SetFillColor(d_samp[typ]["f_color"])
            histo[k].SetLineColor(d_samp[typ]["l_color"])
            histo[k].SetMarkerStyle(20)
        elif d_samp[typ]["type"] == "sig":
            histo[k].SetFillColor(0)
            histo[k].SetLineColor(d_samp[typ]["l_color"])
            histo[k].SetMarkerStyle(0)
            histo[k].SetMarkerSize(0)
            histo[k].SetLineStyle(9)
            histo[k].SetLineWidth(2)
        if writetofile:
            histo[k].Write()

def getHistograms(fname):
    histo = {}
    f1 = R.TFile(fname)
    dirlist = f1.GetListOfKeys()
    it = dirlist.MakeIterator()
    key = it.Next()
    while key:
        cl = R.gROOT.GetClass(key.GetClassName());
        if cl.InheritsFrom("TH1D") or cl.InheritsFrom("TH2D"):
            obj = key.ReadObj()
            histo[obj.GetName().replace("h_","")] = obj.Clone()
            histo[obj.GetName().replace("h_","")].SetDirectory(0)
            key = it.Next()
        else:
            key = it.Next()
            continue
    f1.Close()
    return histo

def getTreeName(fname):
    f1 = R.TFile(fname)
    dirlist = f1.GetListOfKeys()
    it = dirlist.MakeIterator()
    key = it.Next()
    while key:
        cl = R.gROOT.GetClass(key.GetClassName());
        if cl.InheritsFrom("TTree"):
            obj = key.ReadObj()
            if obj.GetName() in ["CutBookkeepers","MetaTree"]: 
                key = it.Next()
                continue
            return obj.GetName()
        else:
            key = it.Next()
            continue
    f1.Close()
    return "noname"


def getDataFrames(mypath, nev = 0): 
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    df = {}
    files = {}
    for of in onlyfiles:
        if not "merged" in of or not of.endswith(".root"): continue
        sp = of.split("_")
        typ = ""
        for s in sp:
            if "merged" in s or s.isnumeric(): break
            typ += s
        if not typ in files.keys():
            files[typ] = {"files":[], "treename":""}
        treename = getTreeName(mypath+"/"+of)
        if treename == "noname":
            print("ERROR \t Could not find any TTree in %s"%(mypath+"/"+of))
            continue
        files[typ]["treename"] = treename
        files[typ]["files"].append(mypath+"/"+of)
        
        #print(typ)
        #if not typ == "singleTop": continue
        #df[typ] = R.Experimental.MakeNTupleDataFrame("mini",mypath+"/"+of)#("%s_NoSys"%typ,mypath+"/"+of)
    for typ in files.keys():
        print("Adding %i files for %s"%(len(files[typ]["files"]),typ))
        df[typ] = R.RDataFrame(files[typ]["treename"],files[typ]["files"])
        if nev:
            df[typ] = df[typ].Range(nev)
    return df

def getRatio1D(hT,hL,vb=0):
    asym = R.TGraphAsymmErrors();
    hR = hT.Clone(hT.GetName().replace("hT","hE"))
    hR.Divide(hT,hL,1.,1.,'b')
    if vb: print(":::->Dividing T = %.2f on L = %.2f" %(hT.Integral(),hL.Integral()))
    asym.Divide(hT,hL,"cl=0.683 b(1,1) mode")
    for i in range(0,hR.GetNbinsX()+1):
        hR.SetBinError(i+1,asym.GetErrorY(i))
    return hR


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
            print(k)
            isData = "data" in k
            if isData:
                print(df[k].GetColumnNames())

            if not isData:
                df[k] = df[k].Define("scaletolumi","(RandomRunNumber) < 320000 ? 36207.65 : (((RandomRunNumber) > 320000 && (RandomRunNumber) < 348000) ? 44307.4 : 58450.1)")

            df[k] = df[k].Define("new_xsec","(DatasetNumber == 308981) ? (0.30649*69.594)/80000. : 0.0")

            # Baseline leptons
            ele_filter = "lepFlavor==1 && lepPassOR > 0 && (lepEta <= 2.47 && lepEta >= -2.47) && ((lepZ0SinTheta)<=0.5 && (lepZ0SinTheta)>=-0.5)"
            muo_filter = "lepFlavor==2 && lepPassOR > 0 && (lepEta <= 2.7  && lepEta >= -2.7) && lepLoose > 0 && ((lepZ0SinTheta)<=0.5 && (lepZ0SinTheta)>=-0.5)"
            
            isGoodLepton = f"{ele_filter} || {muo_filter}"
            df[k] = df[k].Define("isGoodLepton", isGoodLepton)

            df[k] = df[k].Define("ele_BL", ele_filter) 
            df[k] = df[k].Define("muo_BL", muo_filter) 
            df[k] = df[k].Define("nlep_BL","ROOT::VecOps::Sum(ele_BL)+ROOT::VecOps::Sum(muo_BL)")

            # Signal leptons
            df[k] = df[k].Define("ele_SG","ele_BL && lepIsoLoose_VarRad && lepTight && (lepD0Sig <= 5 && lepD0Sig >= -5)") 
            df[k] = df[k].Define("muo_SG","muo_BL && lepIsoLoose_VarRad && (lepD0Sig <= 3 && lepD0Sig >= -3)")
            df[k] = df[k].Define("nlep_SG","ROOT::VecOps::Sum(ele_SG)+ROOT::VecOps::Sum(muo_SG)")
            df[k] = df[k].Filter("nlep_BL == 3")
            df[k] = df[k].Filter("nlep_SG == 3")
            
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
                    
            # Momentum cuts 
            df[k] = df[k].Filter("lep1_Pt > 40")
            df[k] = df[k].Filter("lep2_Pt > 40")
            df[k] = df[k].Filter("lep3_Pt > 15")
            
            # Significance of missing transverse energy cut
            df[k] = df[k].Filter("met_Sign>= 5")
            
            # Triggers
            Years = ["2015", "2016", "2017", "2018"]
            
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
                print(year)
                y_t = Trig_year[year]
                other_years = [years[i] for i in range(4) if years[i] != year ]
              
                isTriggered = y_t["triggered"]
                isMatch = y_t["match"]
                df[k] = df[k].Filter(f"{isTriggered} || is{other_years[0]} || is{other_years[1]} || is{other_years[2]}")
                df[k] = df[k].Define(f"lep_trig_{year}", f"{isMatch} +  || is{other_years[0]} || is{other_years[1]} || is{other_years[2]}")
                df[k] = df.Filter(f"ROOT::VecOps::Sum(lep_trig_{year}) == 2")
                    
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
            
            # Energy of missing transverse momentum
            df[k] = df[k].Define("met_et","met_Et")
            
            # Phi of missing transverse momentum
            df[k] = df[k].Define("met_phi","met_Phi")
            
            # Delta R
            df[k] = df[k].Define("deltaR","deltaR(lepEta[isGoodLepton], lepPhi[isGoodLepton], 0, 1)")
            
            # Invariant Mass (lll)
            df[k] = df[k].Define("mlll","ComputeInvariantMass(lepPt[isGoodLepton], lepEta[isGoodLepton], lepPhi[isGoodLepton], lepM[isGoodLepton])")
            
            # Invariant Mass (ll-OSSF)
            df[k] = df[k].Define("mll_OSSF","OSSFInvariantMass(lepPt[isGoodLepton], lepEta[isGoodLepton], lepPhi[isGoodLepton], lepM[isGoodLepton], lepCharge[isGoodLepton], lepType[isGoodLepton])")
            
            # Ht(lll)
            df[k] = df[k].Define("Ht_lll","ROOT::VecOps::Sum(lepPt[isGoodLepton])")
            
            # Ht (SS)
            df[k] = df[k].Define("Ht_SS","SSHt(lepPt[isGoodLepton], lepCharge[isGoodLepton])")
            
            # Ht(lll) + Energy of missing transverse momentum
            df[k] = df[k].Define("Ht_met_Et","Ht_lll + met_Et")
            
            # Invariant mass of leading jet-pair.
            df[k] = df[k].Define("M_jj","mjj")
            

                        

            # nBjets
            df[k] = df[k].Define("nbjet85","ROOT::VecOps::Sum(bjet85)")
            df[k] = df[k].Define("nbjet77","ROOT::VecOps::Sum(bjet77)")
            
            # Flavour combo
            df[k] = df[k].Define("flcomp","flavourComp3L(lepFlavor[ele_BL || muo_BL])")
            
            
            # HISTOGRAMS
            histo["flcomp_%s"%(k)] = df[k].Histo1D(("h_%s_%s"%("flcomp",k),"h_%s_%s"%("flcomp",k),len(fldic.keys()),0,len(fldic.keys())),"flcomp","wgt_SG")
            
            histo["MT2_12_%s"%(k)] = df[k].Histo1D(("h_%s_%s"%("MT2_12",k),"h_%s_%s;m_{T}^{2}(12) [GeV];Entries"%("MT2_12",k),40,0,400),"MT2_12","wgt_SG")
            histo["MT2_13_%s"%(k)] = df[k].Histo1D(("h_%s_%s"%("MT2_13",k),"h_%s_%s;m_{T}^{2}(13) [GeV];Entries"%("MT2_13",k),40,0,400),"MT2_13","wgt_SG")
            histo["MT2_23_%s"%(k)] = df[k].Histo1D(("h_%s_%s"%("MT2_23",k),"h_%s_%s;m_{T}^{2}(23) [GeV];Entries"%("MT2_23",k),40,0,400),"MT2_23","wgt_SG")
            
            bins_dic = featdic["lep1_E"]
            histo["lep1_E_%s"%(k)] = df[k].Histo1D(("h_%s_%s"%("lep1_E",k),"h_%s_%s;E(l1) [GeV];Entries"%("lep1_E",k),bins_dic["nr_bins"],bins_dic["min"],bins_dic["max"]),"lep1_E","wgt_SG")
            histo["lep2_E_%s"%(k)] = df[k].Histo1D(("h_%s_%s"%("lep2_E",k),"h_%s_%s;E(l2) [GeV];Entries"%("lep2_E",k),bins_dic["nr_bins"],bins_dic["min"],bins_dic["max"]),"lep2_E","wgt_SG")
            histo["lep3_E_%s"%(k)] = df[k].Histo1D(("h_%s_%s"%("lep3_E",k),"h_%s_%s;E(l3) [GeV];Entries"%("lep3_E",k),bins_dic["nr_bins"],bins_dic["min"],bins_dic["max"]),"lep3_E","wgt_SG")
            
            bins_dic = featdic["lep1_Mt"]
            histo["lep1_Mt_%s"%(k)] = df[k].Histo1D(("h_%s_%s"%("lep1_Mt",k),"h_%s_%s;Mt(l1) [GeV];Entries"%("lep1_Mt",k),bins_dic["nr_bins"],bins_dic["min"],bins_dic["max"]),"lep1_Mt","wgt_SG")
            histo["lep2_Mt_%s"%(k)] = df[k].Histo1D(("h_%s_%s"%("lep2_Mt",k),"h_%s_%s;Mt(l2) [GeV];Entries"%("lep2_Mt",k),bins_dic["nr_bins"],bins_dic["min"],bins_dic["max"]),"lep2_Mt","wgt_SG")
            histo["lep3_Mt_%s"%(k)] = df[k].Histo1D(("h_%s_%s"%("lep3_Mt",k),"h_%s_%s;Mt(l3) [GeV];Entries"%("lep3_Mt",k),bins_dic["nr_bins"],bins_dic["min"],bins_dic["max"]),"lep3_Mt","wgt_SG")
            
            histo["nbjet85_%s"%(k)] = df[k].Histo1D(("h_%s_%s"%("nbjet85",k),"h_%s_%s"%("nbjet85",k),9,0,8),"nbjet85","wgt_SG")        
            histo["nbjet77_%s"%(k)] = df[k].Histo1D(("h_%s_%s"%("nbjet77",k),"h_%s_%s"%("nbjet77",k),7,0,6),"nbjet77","wgt_SG")
            
            histo["met_phi_%s"%k] = df[k].Filter("nlep_SG == 3").Histo1D(("h_%s_%s"%("met_phi",k),"h_%s_%s; Phi of missing transvere momentum;Entries"%("met_phi",k),100,-3.5,3.5),"met_phi","wgt_SG")

            histo["met_et_%s"%k] = df[k].Filter("nlep_SG == 3").Histo1D(("h_%s_%s"%("met_et",k),"h_%s_%s; Energy of missing transverse momentum [GeV];Entries"%("met_et",k),50,0,500),"met_et","wgt_SG")
            
            histo["nlep_SG_%s"%k] = df[k].Histo1D(("nlep_SG_%s"%k,"nlep_SG_%s"%k,6,0,5),"nlep_SG","wgt_SG")
            
            histo["njet_SG_%s"%k] = df[k].Histo1D(("njet_SG_%s"%k,"njet_SG_%s"%k,10,0,10),"njet_SG","wgt_SG")
            
            histo["deltaR_%s"%k] = df[k].Histo1D(("deltaR_%s"%k,"deltaR_%s"%k,20,0,8),"deltaR","wgt_SG")
            
            histo["mlll_%s"%k] = df[k].Histo1D(("mlll_%s"%k,"mlll_%s"%k,40,0,400),"mlll","wgt_SG")
            
            histo["mll_OSSF_%s"%k] = df[k].Histo1D(("mll_OSSF_%s"%k,"mll_OSSF_%s"%k,40,0,400),"mll_OSSF","wgt_SG")
            
            histo["Ht_lll_%s"%k] = df[k].Histo1D(("Ht_lll_%s"%k,"Ht_lll_%s"%k,50,0,500),"Ht_lll","wgt_SG")

            histo["Ht_SS_%s"%k] = df[k].Histo1D(("Ht_SS_%s"%k,"Ht_SS_%s"%k,40,0,400),"Ht_SS","wgt_SG")

            histo["Ht_met_Et_%s"%k] = df[k].Histo1D(("Ht_met_Et_%s"%k,"Ht_met_Et_%s"%k,80,0,800),"Ht_met_Et","wgt_SG")
            
            histo["M_jj_%s"%k] = df[k].Histo1D(("M_jj_%s"%k,"M_jj_%s"%k,80,0,800),"M_jj","wgt_SG")
            
            histo["met_Sign_%s"%k] = df[k].Histo1D(("met_Sign_%s"%k,"met_Sign_%s"%k,20,0,20),"met_Sign","wgt_SG")


            


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


        

        return histo

histo = runANA("/storage/shared/data/master_students/William_Sakarias/data/PHYS_3LBkgs_mc16e","/storage/shared/data/master_students/William_Sakarias/data/data18",everyN,fldic,histo,allhisto)


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
%jsroot on
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
        p.can.SaveAs(f"../../Figures/FeaturesHistograms/{feature}.pdf");
        p.can.Draw()
    except:
        print(f"Was not able to plot histogram for {feature}.")

#df.to_hdf("/storage/shared/data/master_students/William_Sakarias/" + k +"mc.hdf5","mini")
#df.to_hdf("/storage/shared/data/master_students/William_Sakarias/3L_data.hdf5","mini")