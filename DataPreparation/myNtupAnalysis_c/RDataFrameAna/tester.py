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
R.gSystem.AddDynamicPath("-I/home/wohirst/myNtupAnalysis/RDataFrameAna")
R.gInterpreter.Declare('#include "helperFunctions.h"') # Header with the definition of the myFilter function
R.gSystem.Load("helperFunctions_cxx.so") # Library with the myFilter function

from IPython.display import display, HTML
display(HTML("<style>.container { width:85% !important; }</style>"))

# +
#nh = int(sys.argv[1])
# -

hname = "MET_2L_mm"
#if len(sys.argv) > 2:
    #hname = sys.argv[2]

fldic = {"eee":0,
         "eem":1,
         "emm":2,
         "mmm":3,
         "mme":4,
         "mee":5,
         "ee":6,
         "mm":7,
         "em":8,
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
                        "nr_bins": 100, "min" : 0, "max" : 300},
           "lep2_Pt"  : {"xlabel":"P_{t}(l_{1}) [GeV]",
                        "nr_bins": 100, "min" : 0, "max" : 300},
           "lep3_Pt"  : {"xlabel":"P_{t}(l_{1}) [GeV]",
                        "nr_bins": 100, "min" : 0, "max" : 300},
           "lep1_Eta" : {"xlabel":"\eta(l_{1})",
                        "nr_bins": 100, "min" : -3, "max" : 3},
           "lep2_Eta" : {"xlabel":"\eta(l_{2})",
                        "nr_bins": 100, "min" : -3, "max" : 3},
           "lep3_Eta" : {"xlabel":"\eta(l_{3})",
                        "nr_bins": 100, "min" : -3, "max" : 3},
           "lep1_Phi" : {"xlabel":"\phi(l_{1})",
                        "nr_bins": 100, "min" : -3.5, "max" : 3.5},
           "lep2_Phi" : {"xlabel":"\phi(l_{2})",
                        "nr_bins": 100, "min" : -3.5, "max" : 3.5},
           "lep3_Phi" : {"xlabel":"\phi(l_{3})",
                        "nr_bins": 100, "min" : -3.5, "max" : 3.5},
           "lep1_M"   : {"xlabel": "M(l_{1})[GeV]",
                        "nr_bins": 100, "min" : 0, "max" : 300},
           "lep2_M"   : {"xlabel": "M(l_{2})[GeV]",
                        "nr_bins": 100, "min" : 0, "max" : 300},
           "lep3_M"   : {"xlabel": "M(l_{3})[GeV]",
                        "nr_bins": 100, "min" : 0, "max" : 300},
           "lep1_Ptcone30" : {"xlabel": "P_{t}cone30(l_{1})",
                        "nr_bins": 100, "min" : 0, "max" : 300},
           "lep2_Ptcone30" : {"xlabel": "P_{t}cone30(l_{1})",
                        "nr_bins": 100, "min" : 0, "max" : 300},
           "lep3_Ptcone30" : {"xlabel": "P_{t}cone30(l_{1})",
                        "nr_bins": 100, "min" : 0, "max" : 300},
           "lep1_Z0" : {"xlabel": "Z_{0}(l_{1})",
                        "nr_bins": 100, "min" : 0, "max" : 300},
           "lep2_Z0" : {"xlabel": "Z_{0}(l_{2})",
                        "nr_bins": 100, "min" : 0, "max" : 300},
           "lep3_Z0" : {"xlabel": "Z_{0}(l_{3})",
                        "nr_bins": 100, "min" : 0, "max" : 300},
           "met_et"   : {"xlabel":"E_{T}^{miss}[GeV]"},
           "met_phi"  : {"xlabel":"\phi (miss)"}
           
}
Nlep = 3
lepv = ["lepPt","lepEta","lepPhi", "lepM",
        "lepPtcone30",
        "lepZ0"]



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

# R.gInterpreter.Declare("""
#     TCanvas c("c","x hist");
#     auto drawHisto = [](TH1D &h_){		     
# 		   c.cd();
# 		   h_.Draw();
# 		   c.Update();    
# 		 };


#     """)

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

        #mypath = "/storage/eirikgr/ANAntuples/PHYS_Data/"
        if isdir(mypath_data):
            df_data = getDataFrames(mypath_data)
            print("Loading %s into dataframe with keys %s" %(mypath_data,",".join(df_data.keys())))
        else:
            df_data = {}

        df = {**df_mc,**df_data}
        #print(df["Wjets"].GetColumnNames())
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

            df[k] = df[k].Define("ele_BL", ele_filter) #((lepZ0SinTheta)<=0.5 && (lepZ0SinTheta)>=-0.5) &&
            df[k] = df[k].Define("muo_BL", muo_filter) #((lepZ0SinTheta)<=0.5 && (lepZ0SinTheta)>=-0.5) &&

            df[k] = df[k].Define("nlep_BL","ROOT::VecOps::Sum(ele_BL)+ROOT::VecOps::Sum(muo_BL)")

            # Signal leptons
            df[k] = df[k].Define("ele_SG","ele_BL && lepIsoLoose_VarRad && lepTight && (lepD0Sig <= 5 && lepD0Sig >= -5)") #&& lepTight && (lepD0Sig <= 5 && lepD0Sig >= -5)
            df[k] = df[k].Define("muo_SG","muo_BL && lepIsoLoose_VarRad && (lepD0Sig <= 3 && lepD0Sig >= -3)") #&& (lepD0Sig <= 3 && lepD0Sig >= -3)

            df[k] = df[k].Define("nlep_SG","ROOT::VecOps::Sum(ele_SG)+ROOT::VecOps::Sum(muo_SG)")
            df[k] = df[k].Filter("nlep_SG == 3")
                        
            if not isData:

                df[k] = df[k].Define("is2015","RandomRunNumber <= 284500")
                df[k] = df[k].Define("is2016","(RandomRunNumber > 284500 && RandomRunNumber < 320000)")
                df[k] = df[k].Define("is2017","(RandomRunNumber > 320000 && RandomRunNumber < 348000)")
                df[k] = df[k].Define("is2018","RandomRunNumber > 348000")

                #df[k] = df[k].Define("lepwgt_BL","getSF(lepBLRecoSF[ele_BL || muo_BL])")
                df[k] = df[k].Define("lepwgt_SG","getSF(lepRecoSF[ele_SG || muo_SG])")

                #df[k] = df[k].Define("trgwgt_BL","getSF(lepBLTrigSF[ele_BL || muo_BL])")
                df[k] = df[k].Define("trgwgt_SG","getSF(lepTrigSF[ele_SG || muo_SG])")

                #df[k] = df[k].Define("wgt_BL","(new_xsec ? (new_xsec) : (genWeight))*eventWeight*jvtWeight*bTagWeight*pileupWeight*scaletolumi*lepwgt_BL*trgwgt_BL")
                df[k] = df[k].Define("wgt_SG","(new_xsec ? (new_xsec) : (genWeight))*eventWeight*jvtWeight*bTagWeight*pileupWeight*scaletolumi*lepwgt_SG*trgwgt_SG")

                #df[k] = df[k].Define("wgt_EV_BL","(eventWeight*jvtWeight*bTagWeight*pileupWeight*scaletolumi*lepwgt_BL*trgwgt_BL)")
                df[k] = df[k].Define("wgt_EV_SG","(eventWeight*jvtWeight*bTagWeight*pileupWeight*scaletolumi*lepwgt_SG*trgwgt_SG)")

            else:
                df[k] = df[k].Define("is2015","(RunNumber >= 276262 && RunNumber <= 284484)")
                df[k] = df[k].Define("is2016","(RunNumber >= 297730 && RunNumber <= 311481)")
                df[k] = df[k].Define("is2017","(RunNumber >= 325713 && RunNumber <= 340453)")
                df[k] = df[k].Define("is2018","RunNumber >= 348885")

                #df[k] = df[k].Define("wgt_BL","1.0")
                df[k] = df[k].Define("wgt_SG","1.0")
                df[k] = df[k].Define("wgt_EV","1.0")


            #df[k] = df[k].Filter("eventIsTriggered_1L","1L trigger")
            #df[k] = df[k].Filter("ROOT::VecOps::Sum(lepIsTrigMatched_1L[ele_BL || muo_BL]) > 0","Trigger Matched")
            
            if not nEvents:
                this_nEvents = int(df[k].Count().GetValue())
                nEvents += this_nEvents
                print("Loading %s with %.0f events. Now %.0f events"%(k,this_nEvents,nEvents))
            else:
                print("Loading %s"%(k))
                
            histo["nlep_SG_%s"%k] = df[k].Histo1D(("nlep_SG_%s"%k,"nlep_SG_%s"%k,10,0,10),"nlep_SG","wgt_SG")
            

            #df[k] = df[k].Define("Zcand_mass","getLeptonsFromZ(lepCharge[ele_SG > 0 || muo_SG > 0], lepFlavor[ele_SG > 0 || muo_SG > 0], lepPt[ele_SG > 0 || muo_SG > 0], lepEta[ele_SG > 0 || muo_SG > 0], lepPhi[ele_SG > 0 || muo_SG > 0], lepM[ele_SG > 0 || muo_SG > 0], met_Et, met_Phi).first")
            #df[k] = df[k].Define("Wcand_mass","getLeptonsFromZ(lepCharge[ele_SG > 0 || muo_SG > 0], lepFlavor[ele_SG > 0 || muo_SG > 0], lepPt[ele_SG > 0 || muo_SG > 0], lepEta[ele_SG > 0 || muo_SG > 0], lepPhi[ele_SG > 0 || muo_SG > 0], lepM[ele_SG > 0 || muo_SG > 0], met_Et, met_Phi).second")

            df[k] = df[k].Define("isZlep1","getZlep1()")
            df[k] = df[k].Define("isZlep2","getZlep2()")
            df[k] = df[k].Define("isWlep1","getWlep1()")

            df[k] = df[k].Filter("nlep_SG == 3").Define("MT2_12","calcMT2(lepPt[ele_SG > 0 || muo_SG > 0], lepEta[ele_SG > 0 || muo_SG > 0], lepPhi[ele_SG > 0 || muo_SG > 0], lepM[ele_SG > 0 || muo_SG > 0], met_Et, met_Phi, 0, 1)")
            df[k] = df[k].Filter("nlep_SG == 3").Define("MT2_13","calcMT2(lepPt[ele_SG > 0 || muo_SG > 0], lepEta[ele_SG > 0 || muo_SG > 0], lepPhi[ele_SG > 0 || muo_SG > 0], lepM[ele_SG > 0 || muo_SG > 0], met_Et, met_Phi, 0, 2)")
            df[k] = df[k].Filter("nlep_SG == 3").Define("MT2_23","calcMT2(lepPt[ele_SG > 0 || muo_SG > 0], lepEta[ele_SG > 0 || muo_SG > 0], lepPhi[ele_SG > 0 || muo_SG > 0], lepM[ele_SG > 0 || muo_SG > 0], met_Et, met_Phi, 1, 2)")
            
            for i in range(Nlep):
                df[k] = df[k].Define("lep%i_flav"%(i+1),"getTypeTimesCharge(lepCharge[isGoodLepton],lepType[isGoodLepton],%i)"%(i))
                for v in lepv:
                    if "lep" in v:
                        var = v.replace("lep","")
                    else:
                        var = v
                    print(i, v)
                    bins_dic = featdic["lep%i_%s"%(i+1,var)]
                    df[k] = df[k].Define("lep%i_%s"%(i+1,var),"getVar(lep%s[isGoodLepton],%i)"%(var,i))
                    histo["lep%i_%s_%s"%(i+1,var,k)] = df[k].Histo1D(("lep%i_%s_%s"%(i+1,var,k),
                          "lep%i_%s_%s;Feature;Entries"%(i+1,var,k),bins_dic["nr_bins"],bins_dic["min"],bins_dic["max"]),
                          "lep%i_%s"%(i+1,var),"wgt_SG")
            
            
            # Jets
            df[k] = df[k].Define("jet_BL","jetPt >= 20 && (jetEta <= 2.8 && jetEta >= -2.8)")
            df[k] = df[k].Define("jet_SG","jet_BL && (jetPt > 60 || (jetPt <=60 && jetJVT <= 0.91 && jetJVT >= -0.91))")
            df[k] = df[k].Define("bjet85","jet_SG && jetdl1r>=0.665")
            df[k] = df[k].Define("bjet77","jet_SG && jetdl1r>=2.195")
            df[k] = df[k].Define("jet_SG_pT","jetPt[jet_SG > 0]")
            df[k] = df[k].Define("jet_SG_eta","jetEta[jet_SG > 0]")

            #df[k] = df[k].Define("minDR_jetlep1","deltaRlepjet(lepPt[0],lepEta[0],lepPhi[0],lepM[0],jetPt[jet_BL > 0],jetEta[jet_BL > 0],jetPhi[jet_BL > 0],jetM[jet_BL > 0])")
            
            #histo["minDR_jetlep1_%s"%k] = df[k].Histo1D(("h_%s_%s"%("minDR_jetlep1",k),"h_%s_%s;min #DeltaR(lep1,jet);Entries"%("minDR_jetlep1",k),100,0,2),"minDR_jetlep1","wgt_SG") 

            df[k] = df[k].Define("lep_BL_pT","lepPt[ele_BL > 0 || muo_BL > 0]")
            df[k] = df[k].Define("lep_SG_pT","lepPt[ele_SG > 0 || muo_SG > 0]")
            
            #met_et
            df[k] = df[k].Define("met_et","met_Et")
            histo["met_et_%s"%k] = df[k].Filter("nlep_SG == 3").Histo1D(("h_%s_%s"%("met_et",k),"h_%s_%s; Energy of missing transverse momentum [GeV];Entries"%("met_et",k),100,0,500),"met_et","wgt_SG")
            
            #met_phi
            df[k] = df[k].Define("met_phi","met_Phi")
            histo["met_phi_%s"%k] = df[k].Filter("nlep_SG == 3").Histo1D(("h_%s_%s"%("met_phi",k),"h_%s_%s; Phi of missing transvere momentum;Entries"%("met_phi",k),100,-3.5,3.5),"met_phi","wgt_SG")


            df[k] = df[k].Define("ele_SG_pT","lepPt[ele_SG > 0]")
            df[k] = df[k].Define("ele_SG_eta","lepEta[ele_SG > 0]")
            df[k] = df[k].Define("muo_SG_pT","lepPt[muo_SG > 0]")
            df[k] = df[k].Define("muo_SG_eta","lepEta[muo_SG > 0]")
            df[k] = df[k].Define("njet_SG","ROOT::VecOps::Sum(jet_SG)")
            
            histo["njet_SG_%s"%k] = df[k].Histo1D(("njet_SG_%s"%k,"njet_SG_%s"%k,10,0,10),"njet_SG","wgt_SG")


            df[k] = df[k].Define("nbjet85","ROOT::VecOps::Sum(bjet85)")
            df[k] = df[k].Define("nbjet77","ROOT::VecOps::Sum(bjet77)")

            df[k] = df[k].Define("flcomp","flavourComp3L(lepFlavor[ele_BL || muo_BL])")
            
            histo["flcomp_%s"%(k)] = df[k].Histo1D(("h_%s_%s"%("flcomp",k),"h_%s_%s"%("flcomp",k),len(fldic.keys()),0,len(fldic.keys())),"flcomp","wgt_SG")
            

            etabins = array.array('f',[-2.7,-2.0,-1.8,-1.52,-1.37,-1.2,-1.0,-0.8,-0.6,-0.4,0.0,0.4,0.6,0.8,1.0,1.2,1.37,1.52,1.8,2.0,2.7])

            xbins = array.array('f',[0,10,20,25,30,40,50,60,80,100,200])
            ybins = array.array('f',[0,10,20,25,30,40,50,60,80,100,200])

            for nlep in ["2L"]:#,"2L"]:#,"2L"]:
                trigs = "(1"#(eventIsTriggered_%s"%nlep# && ROOT::VecOps::Sum(lepIsTrigMatched_%s[ele_BL || muo_BL]) > 0"%(nlep,nlep)
                for flk in ["all"]:#fldic.keys():
                    comp = fldic[flk]
                    if comp <= 8:
                        filterstr = "%s && flcomp == %i)"%(trigs,comp)
                    else:
                        filterstr = "%s)"%(trigs)

                    #histo["minDR_jetlep1_%s_%s_%s"%(nlep,flk,k)] = df[k].Filter(filterstr).Histo1D(("h_%s_%s_%s_%s"%("minDR_jetlep1",nlep,flk,k),"h_%s_%s_%s_%s;min #DeltaR(lep1,jet);Entries"%("minDR_jetlep1",nlep,flk,k),11,-1,10),"minDR_jetlep1","wgt_SG") 

                    histo["Zlep1_%s_%s_%s"%(nlep,flk,k)] = df[k].Filter(filterstr).Histo1D(("h_%s_%s_%s_%s"%("Zlep1",nlep,flk,k),"h_%s_%s_%s_%s;Lepton number;Entries"%("Zlep1",nlep,flk,k),11,-1,10),"isZlep1","wgt_SG")
                    histo["Zlep2_%s_%s_%s"%(nlep,flk,k)] = df[k].Filter(filterstr).Histo1D(("h_%s_%s_%s_%s"%("Zlep2",nlep,flk,k),"h_%s_%s_%s_%s;Lepton number;Entries"%("Zlep2",nlep,flk,k),11,-1,10),"isZlep2","wgt_SG")
                    histo["Wlep1_%s_%s_%s"%(nlep,flk,k)] = df[k].Filter(filterstr).Histo1D(("h_%s_%s_%s_%s"%("Wlep1",nlep,flk,k),"h_%s_%s_%s_%s;Lepton number;Entries"%("Wlep1",nlep,flk,k),11,-1,10),"isWlep1","wgt_SG")

                    #histo["Zcand_mass_%s_%s_%s"%(nlep,flk,k)] = df[k].Filter(filterstr).Histo1D(("h_%s_%s_%s_%s"%("Zcand_mass",nlep,flk,k),"h_%s_%s_%s_%s;m_{ll}^{Z-cand} [GeV];Entries"%("Zcand_mass",nlep,flk,k),101,-10,1000),"Zcand_mass","wgt_SG")
                    #histo["Wcand_mass_%s_%s_%s"%(nlep,flk,k)] = df[k].Filter(filterstr).Histo1D(("h_%s_%s_%s_%s"%("Wcand_mass",nlep,flk,k),"h_%s_%s_%s_%s;m_{T}^{W-cand} [GeV];Entries"%("Wcand_mass",nlep,flk,k),101,-10,1000),"Wcand_mass","wgt_SG")

                    #histo["MET_%s_%s_%s"%(nlep,flk,k)] = df[k].Filter(filterstr).Histo1D(("h_%s_%s_%s_%s"%("MET",nlep,flk,k),"h_%s_%s_%s_%s;Missing Transverse Energy [GeV]; Entries"%("MET",nlep,flk,k),500,0,1000),"met_Et","wgt_SG")

                    histo["MT2_12_%s_%s_%s"%(nlep,flk,k)] = df[k].Filter(filterstr).Histo1D(("h_%s_%s_%s_%s"%("MT2_12",nlep,flk,k),"h_%s_%s_%s_%s;m_{T}^{2}(12) [GeV];Entries"%("MT2_12",nlep,flk,k),500,0,400),"MT2_12","wgt_SG")
                    histo["MT2_13_%s_%s_%s"%(nlep,flk,k)] = df[k].Filter(filterstr).Histo1D(("h_%s_%s_%s_%s"%("MT2_13",nlep,flk,k),"h_%s_%s_%s_%s;m_{T}^{2}(13) [GeV];Entries"%("MT2_13",nlep,flk,k),500,0,400),"MT2_13","wgt_SG")
                    histo["MT2_23_%s_%s_%s"%(nlep,flk,k)] = df[k].Filter(filterstr).Histo1D(("h_%s_%s_%s_%s"%("MT2_23",nlep,flk,k),"h_%s_%s_%s_%s;m_{T}^{2}(23) [GeV];Entries"%("MT2_23",nlep,flk,k),500,0,400),"MT2_23","wgt_SG")

                    histo["lepZ0SinTheta_%s_%s_%s"%(nlep,flk,k)] = df[k].Filter(filterstr).Histo1D(("h_%s_%s_%s_%s"%("lepZ0SinTheta",nlep,flk,k),"h_%s_%s_%s_%s;Z_{0}#sin#theta;Entries"%("lepZ0SinTheta",nlep,flk,k),1500,-5,10),"lepZ0SinTheta","wgt_SG")
                    histo["lepD0Sig_%s_%s_%s"%(nlep,flk,k)] = df[k].Filter(filterstr).Histo1D(("h_%s_%s_%s_%s"%("lepD0Sig",nlep,flk,k),"h_%s_%s_%s_%s;d_{0}/#sigma(d_{0});Entries"%("lepD0Sig",nlep,flk,k),400,-20,20),"lepD0Sig","wgt_SG")

                    histo["lepPt_ele_SG_%s_%s_%s"%(nlep,flk,k)] = df[k].Filter(filterstr).Histo1D(("h_%s_%s_%s_%s"%("lepPt_ele_SG",nlep,flk,k),"h_%s_%s_%s_%s;p_{T}^{SG leptons} [GeV];Entries"%("lepPt_ele_SG",nlep,flk,k),len(xbins)-1,xbins),"ele_SG_pT","wgt_SG")


                    if nh == 1: continue

                    

                    histo["lepEta_ele_SG_%s_%s_%s"%(nlep,flk,k)] = df[k].Filter(filterstr).Histo1D(("h_%s_%s_%s_%s"%("lepEta_ele_SG",nlep,flk,k),"h_%s_%s_%s_%s"%("lepEta_ele_SG",nlep,flk,k),len(etabins)-1,etabins),"ele_SG_eta","wgt_SG")
                    histo["lepEta_muo_SG_%s_%s_%s"%(nlep,flk,k)] = df[k].Filter(filterstr).Histo1D(("h_%s_%s_%s_%s"%("lepEta_muo_SG",nlep,flk,k),"h_%s_%s_%s_%s"%("lepEta_muo_SG",nlep,flk,k),len(etabins)-1,etabins),"muo_SG_eta","wgt_SG")

                    
                    if nh == 2: continue

                    histo["lepPt_muo_SG_%s_%s_%s"%(nlep,flk,k)] = df[k].Filter(filterstr).Histo1D(("h_%s_%s_%s_%s"%("lepPt_muo_SG",nlep,flk,k),"h_%s_%s_%s_%s"%("lepPt_muo_SG",nlep,flk,k),len(xbins)-1,xbins),"muo_SG_pT","wgt_SG")


                    if nh == 3: continue

                    
                    histo["nbjet85_%s_%s_%s"%(nlep,flk,k)] = df[k].Filter(filterstr).Histo1D(("h_%s_%s_%s_%s"%("nbjet85",nlep,flk,k),"h_%s_%s_%s_%s"%("nbjet85",nlep,flk,k),20,0,20),"nbjet85","wgt_SG")
                    
                    histo["nbjet77_%s_%s_%s"%(nlep,flk,k)] = df[k].Filter(filterstr).Histo1D(("h_%s_%s_%s_%s"%("nbjet77",nlep,flk,k),"h_%s_%s_%s_%s"%("nbjet77",nlep,flk,k),20,0,20),"nbjet77","wgt_SG")
                    
                    histo["njet_SG_%s_%s_%s"%(nlep,flk,k)] = df[k].Filter(filterstr).Histo1D(("h_%s_%s_%s_%s"%("njet_SG",nlep,flk,k),"h_%s_%s_%s_%s"%("njet_SG",nlep,flk,k),20,0,20),"njet_SG","wgt_SG")

                    if nh == 4: continue
                    
                    histo["jetPt_SG_%s_%s_%s"%(nlep,flk,k)] = df[k].Filter(filterstr).Histo1D(("h_%s_%s_%s_%s"%("jetPt_SG",nlep,flk,k),"h_%s_%s_%s_%s"%("jetPt_SG",nlep,flk,k),len(xbins)-1,xbins),"jet_SG_pT","wgt_SG")

                    if nh == 5: continue
                    

                    if nh == 6: continue
                    
                    
                    histo["lep1Pt_jet1Pt_%s_%s_%s"%(nlep,flk,k)] = df[k].Filter("njet_SG > 0").Define("lep1_pT","ROOT::VecOps::Max(lep_SG_pT)").Define("jet1_pT","ROOT::VecOps::Max(jet_SG_pT)").Filter(filterstr).Histo2D(("h_%s_%s_%s_%s"%("lep1Pt_jet1Pt",nlep,flk,k),"h_%s_%s_%s_%s"%("lep1Pt_jet1Pt",nlep,flk,k),100,0,1000,100,0,1000),"lep1_pT","jet1_pT") #len(xbins)-1,xbins,len(xbins)-1,xbins)

                    if nh == 7: continue
                    
                    histo["r_j1_l1_%s_%s_%s"%(nlep,flk,k)] = df[k].Filter(filterstr).Filter("njet_SG > 0").Define("lep1_pT","ROOT::VecOps::Max(lep_SG_pT)").Define("jet1_pT","ROOT::VecOps::Max(jet_SG_pT)").Define("r_j1_l1","lep1_pT/jet1_pT").Histo1D(("h_%s_%s_%s_%s"%("r_j1_l1",nlep,flk,k),"h_%s_%s_%s_%s"%("r_j1_l1",nlep,flk,k),100,0,20),"r_j1_l1","wgt_SG")

                    print(len(histo.keys()))
                    
        for k in histo.keys():
            allhisto.append(histo[k])

        print("Calculating %i histograms"%len(allhisto))
        #sys.exit()
        start = time.time()
        R.RDF.RunGraphs(allhisto)
        end = time.time()
        print("%10i | %.2f"%(len(allhisto),(end - start)))

        hfile = R.TFile("histograms.root","RECREATE")
        hfile.cd()

        writeHistsToFile(histo, True)


        

        return histo

histo = runANA("/storage/shared/data/master_students/William_Sakarias/data/PHYS_3LBkgs_mc16e","/storage/shared/data/master_students/William_Sakarias/data/data18",everyN,fldic,histo,allhisto)