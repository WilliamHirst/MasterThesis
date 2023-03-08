from genericpath import exists
from os import listdir,walk
import ROOT as R
from os.path import isfile, join, isdir
import array
import sys
from samples import configure_samples
from pyHelperFunctions import *
import re


d_samp,d_type,d_reg = configure_samples()#False,False,True,False,False)

fldic = {"eee":0,
         "eem":1,
         "emm":2,
         "mem":3,
         "mmm":4,
         "mme":5,
         "mee":6,
         "eme":7,
}


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
            histo[k].SetFillColorAlpha(d_samp[typ]["f_color"],0.95)
            histo[k].SetLineColor(d_samp[typ]["f_color"])
            histo[k].SetMarkerStyle(0)
            histo[k].SetMarkerSize(0)
        elif d_samp[typ]["type"] == "data":
            histo[k].SetFillColor(d_samp[typ]["f_color"])
            histo[k].SetLineColor(d_samp[typ]["l_color"])
            histo[k].SetMarkerStyle(20)
        elif d_samp[typ]["type"] == "sig":
            #histo[k].SetFillColor(0)
            # histo[k].SetFillColorAlpha(d_samp[typ]["l_color"],0.3)
            histo[k].SetLineColor(d_samp[typ]["l_color"])
            
            histo[k].SetLineStyle(9)
            histo[k].SetLineWidth(3)
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


def getDataFrames(mypath, nev = 0, storage = ""):
    onlyfiles = []
    for path,dirs,files in walk(mypath):
        #print(path,dirs,files)
        for f in files:
            if isfile(join(path, f)) and f.endswith("_merged_processed.root"):
                #print(join(path,f))
                onlyfiles.append(join(path, f))

 
    df = {}
    files = {}
    for of in onlyfiles:
        if not "merged" in of or not of.endswith(".root"): continue
        sp = of.split("/")[-1].split("_")
        typ = ""
        treename = getTreeName(of)
        for s in sp:
            if "merged" in s: break
            typ += s

        if not typ in files.keys():
            files[typ] = {"files":[], "treename":""}
            
        if removeP0(typ) in files.keys() and "MGPy8EGA" in typ and "p0p0" in typ:
            files[typ] = files[removeP0(typ)]
            treenames = [files[removeP0(typ)]["treename"],treename]
            files[typ]["treename"] = treenames
            files.pop(removeP0(typ))

        elif removeAllP0(typ) in files.keys() and "MGPy8EGA" in typ and "p0" in typ:
            files[typ] = files[removeAllP0(typ)]
            treenames = [files[removeAllP0(typ)]["treename"],treename]
            files[typ]["treename"] = treenames
            files.pop(removeAllP0(typ))
        
        if treename == "noname":
            print("ERROR \t Could not find any TTree in %s"%(of))
            continue
        if files[typ]["treename"] == "":
            files[typ]["treename"] = treename
        files[typ]["files"].append(of)
    

    for typ in files.keys():
        print("Adding %i files for %s"%(len(files[typ]["files"]),typ))
        if isinstance(files[typ]["treename"], str) :
            df[typ] = R.RDataFrame(files[typ]["treename"],files[typ]["files"])
        else:
            files[typ]["files"][0] = f"{storage}/signal/" + files[typ]["files"][0].split("/")[-1]
            df[typ] = R.RDataFrame(files[typ]["treename"][-1],files[typ]["files"])
            # tree1 = files[typ]["files"][0] + "/" +files[typ]["treename"][0]
            # tree2 = files[typ]["files"][1] + "/" +files[typ]["treename"][1]
            # tree3 = files[typ]["files"][2] + "/" +files[typ]["treename"][1]
            # C = R.TChain();
            # C.AddFile(tree1);
            # C.AddFile(tree2);
            # C.AddFile(tree3);
            # df[typ] = R.RDataFrame(C)
        if nev:
            df[typ] = df[typ].Range(nev)
    return df

def addP0(string):
    elem = string.split("WZ")
    mass1 = elem[1][0:3]
    string = string.replace(mass1, f"{mass1}p0")
    mass2 = re.search(f'MGPy8EGA14N23LOC1N2WZ{mass1}p0(.*)3L2L7', string)[1]
    string = string.replace(f"{mass2}3L2", f"{mass2}p03L2")
    return string
def addExtraP0(string):
    string = string.replace("p0", "p0p0")
    return string
def removeP0(string):
    string = string.replace("p0p0", "p0")
    return string
def removeAllP0(string):
    string = string.replace("p0", "")
    return string


def getRatio1D(hT,hL,vb=0):
    asym = R.TGraphAsymmErrors();
    hR = hT.Clone(hT.GetName().replace("hT","hE"))
    hR.Divide(hT,hL,1.,1.,'b')
    if vb: print(":::->Dividing T = %.2f on L = %.2f" %(hT.Integral(),hL.Integral()))
    asym.Divide(hT,hL,"cl=0.683 b(1,1) mode")
    for i in range(0,hR.GetNbinsX()+1):
        hR.SetBinError(i+1,asym.GetErrorY(i))
    return hR

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


def changeTreeName(mypath, storage):
    onlyfiles = []
    for path,dirs,files in walk(mypath):
        for f in files:
            if isfile(join(path, f)) and f.endswith("_merged_processed.root"):
                onlyfiles.append(join(path, f))
 
    df = {}
    files = {}
    for of in onlyfiles:
        if not "merged" in of or not of.endswith(".root"): continue
        sp = of.split("/")[-1].split("_")
        typ = ""
        treename = getTreeName(of)
        for s in sp:
            if "merged" in s: break
            typ += s

        if not typ in files.keys():
            files[typ] = {"files":[], "treename":""}
            
        if removeP0(typ) in files.keys() and "MGPy8EGA" in typ and "p0p0" in typ:
            files[typ] = files[removeP0(typ)]
            treenames = [files[removeP0(typ)]["treename"],treename]
            files[typ]["treename"] = treenames
            files.pop(removeP0(typ))

        elif removeAllP0(typ) in files.keys() and "MGPy8EGA" in typ and "p0" in typ:
            files[typ] = files[removeAllP0(typ)]
            treenames = [files[removeAllP0(typ)]["treename"],treename]
            files[typ]["treename"] = treenames
            files.pop(removeAllP0(typ))
        
        if treename == "noname":
            print("ERROR \t Could not find any TTree in %s"%(of))
            continue
        if files[typ]["treename"] == "":
            files[typ]["treename"] = treename
        files[typ]["files"].append(of)
    

    for typ in files.keys():
        print("Adding %i files for %s"%(len(files[typ]["files"]),typ))
        if not isinstance(files[typ]["treename"], str) :
            name = files[typ]["files"][0].split("/")[-1]  
            oldfile = R.TFile(files[typ]["files"][0], "READ");
            oldtree =  R.TTree();
            oldtree = oldfile.Get(files[typ]["treename"][0].split("/")[-1]);
            oldtree.SetName(files[typ]["treename"][1]);
            newfile = R.TFile.Open(f"{storage}/{name}", "RECREATE");
            newtree = oldtree.CloneTree();
            newfile.Write();
    return 

if __name__ == "__main__":
    data_loc = "/storage/shared/data/master_students/William_Sakarias/data_vOCT2022"
    storage = "/storage/William_Sakarias/William_Data/signal"
    changeTreeName(f"{data_loc}/",storage)

