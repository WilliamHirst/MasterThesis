#!/usr/bin/env python
# coding: utf-8

# In[1]:


import ROOT as R
import glob
import samples
import os
R.EnableImplicitMT(150)



# In[2]:


d_samp, d_type, d_reg = samples.configure_samples()


# In[3]:


astyle = "/home/eirikgr/atlasstyle-00-04-02/"#                                                                                                                                                                                                                                            
astyle ="/mn/felt/u1/eirikgr/bin/atlasstyle-00-04-02/"                                                                                                                                                                                                                                    
astyle ="/home/eirikgr/atlasrootstyle/"                                                                                                                                                                                                                                                   
R.gROOT.SetMacroPath(astyle)                                                                                                                                                                                                                                                                
R.gROOT.LoadMacro("AtlasUtils.C")                                                                                                                                                                                                                                                           
R.gROOT.LoadMacro("AtlasStyle.C")                                                                                                                                                                                                                                                           
R.gROOT.LoadMacro("AtlasLabels.C")                                                                                                                                                                                                                                                          
R.SetAtlasStyle()
R.gStyle.SetErrorX(0.5) 


# In[4]:


os.system(' g++ -shared -fPIC -o helperFunctions.so ./helperFunctions.cxx `root-config --cflags --glibs`')


# In[5]:


R.gSystem.AddDynamicPath("./")
R.gROOT.ProcessLine(".include ./");
R.gInterpreter.AddIncludePath("./");
R.gInterpreter.Declare('#include "./helperFunctions.h"') # Header with the definition of the myFilter function
R.gSystem.Load("./helperFunctions.so") # Library with the myFilter function


# In[6]:


lumi15 = 3219.56
lumi16 = 32988.1
lumi17 = 44307.4
lumi18 = 58450.1
lumi = lumi15+lumi16+lumi17+lumi18
#lumi = lumi18

infiles = []
indir = "/storage/eirikgr/ANAntuples/EXOT0_Bkgs_mc16a/*merged_processed.root"
infiles += glob.glob(indir)
indir = "/storage/eirikgr/ANAntuples/EXOT0_Bkgs_mc16d/*merged_processed.root"
infiles += glob.glob(indir)
indir = "/storage/eirikgr/ANAntuples/EXOT0_Bkgs_mc16e/*merged_processed.root"
infiles += glob.glob(indir)
indir = "/storage/eirikgr/ANAntuples/EXOT0_Data/data15_merged_processed.root"
infiles += glob.glob(indir)
indir = "/storage/eirikgr/ANAntuples/EXOT0_Data/data16_merged_processed.root"
infiles += glob.glob(indir)
indir = "/storage/eirikgr/ANAntuples/EXOT0_Data/data17_merged_processed.root"
infiles += glob.glob(indir)
indir = "/storage/eirikgr/ANAntuples/EXOT0_Data/data18_merged_processed.root"
infiles += glob.glob(indir)


# In[7]:


lumi


# In[8]:


alldic = {}
for infi in infiles:
    sp = infi.split("/")[-1]
    ID = sp.split("_")[0]
    if not ID in alldic:
        alldic[ID] = {"infiles":[]}
    alldic[ID]["infiles"].append(infi)


# In[9]:


dfull = {}
for key in alldic.keys():
    print(key)
    if not 'data' in key:
        dfull[key] = R.RDataFrame("%s_NoSys"%key,alldic[key]["infiles"])
    else:
        dfull[key] = R.RDataFrame(key,alldic[key]["infiles"])


# In[10]:


for key in alldic.keys():
    print("#"*50)
    print("Key %s has %i files:"%(key,len(alldic[key]["infiles"])))
    i = 1
    for f in alldic[key]["infiles"]:
        print("%i) %s"%(i,f)) 
        i+=1


# In[11]:


weight = {}
for key in dfull.keys():
    if 'data' in key:
        weight[key] = '(1.0)'
    elif 'Z' in key and 'jets' in key:
        weight[key] = '(globalDiLepTrigSF*genWeight*eventWeight*leptonWeight*jvtWeight*bTagWeight*pileupWeight*scaletolumi*1000.)'
    else:
        weight[key] = '(globalDiLepTrigSF*genWeight*eventWeight*leptonWeight*jvtWeight*bTagWeight*pileupWeight*scaletolumi)'


# In[12]:


for key in dfull.keys():
    if not "data" in key:
        dfull[key] = dfull[key].Define("scaletolumi","getLumiSF(RandomRunNumber)")
    
    dfull[key] = dfull[key].Define("weight",weight[key])
    dfull[key] = dfull[key].Define("baseline_el","lepPt > 20 && abs(lepEta) < 2.47 && !(abs(lepEta) >= 1.37 && abs(lepEta) <= 1.52) && abs(lepD0Sig) < 5 && abs(lepZ0SinTheta) < 0.5 && lepPassBL && lepLoose && lepPassOR")
    dfull[key] = dfull[key].Define("baseline_mu","lepPt > 20 && abs(lepEta) < 2.50 && abs(lepD0Sig) < 3 && abs(lepZ0SinTheta) < 0.5 && lepHighPt")
    dfull[key] = dfull[key].Define("signal_el","baseline_el && lepMedium && lepIsoFCTight")
    dfull[key] = dfull[key].Define("signal_mu","baseline_mu && lepIsoTightTrackOnly_VarRad")
    dfull[key] = dfull[key].Define("signal_lep","signal_el || signal_mu")
    
    dfull[key] = dfull[key].Define("signal_jet","jetPt > 20 && abs(jetEta) < 4.5 && ((jetJVT > 0.5 && jetPt < 60 && abs(jetEta) < 2.4) || jetPt >= 60) && jetPassOR && jetSignal")
    dfull[key] = dfull[key].Define("signal_bjet","signal_jet && jetdl1r>=0.665")
    dfull[key] = dfull[key].Define("n_bjet","Sum(signal_bjet)")
    dfull[key] = dfull[key].Define("signal_fwdjet","signal_jet && abs(jetEta) > 2.5")
    dfull[key] = dfull[key].Define("JVT_fwdjet","jetJVT[signal_fwdjet]")
    
    
    
    dfull[key] = dfull[key].Define("n_bl_lep","Sum(baseline_mu)+Sum(baseline_el)")
    dfull[key] = dfull[key].Define("n_sg_lep","Sum(signal_lep)")
    dfull[key] = dfull[key].Define("sglep_eta","lepEta[signal_lep]")
    dfull[key] = dfull[key].Define("sgel_eta","lepEta[signal_el]")
    dfull[key] = dfull[key].Define("sgmu_eta","lepEta[signal_mu]")
    dfull[key] = dfull[key].Define("sglep_flavor","lepFlavor[signal_lep]")
    dfull[key] = dfull[key].Define("isOS","isOS(lepCharge[signal_lep])")
    dfull[key] = dfull[key].Define("isSF","isSF(lepFlavor[signal_lep])")
    dfull[key] = dfull[key].Define("isEE","isEE(lepFlavor[signal_lep])")
    dfull[key] = dfull[key].Define("isMM","isMM(lepFlavor[signal_lep])")
    dfull[key] = dfull[key].Define("mll_new","ComputeInvariantMass(lepPt[signal_lep],lepEta[signal_lep],lepPhi[signal_lep],lepM[signal_lep])")


# In[13]:


dpresel = {}
for key in dfull.keys():
    dfull[key] = dfull[key].Filter("(n_bl_lep == 2)","2 baseline leptons")
    dfull[key] = dfull[key].Filter("(n_sg_lep == 2)","2 signal leptons")
    dfull[key] = dfull[key].Filter("checkPt(lepPt[signal_lep],25,20)","pT > 25/20 GeV")
    dfull[key] = dfull[key].Filter("trigMatch_2LTrigOR","matched to trigger")
    dfull[key] = dfull[key].Filter("isOS","Opposite sign")
    dfull[key] = dfull[key].Filter("mll_new > 70","mll > 70 GeV")
    dpresel[key] = dfull[key]


# In[14]:


for key in dfull.keys():
    dpresel[key] = dpresel[key].Define("METrel","getMetRel(lepPt[signal_lep],lepEta[signal_lep],lepPhi[signal_lep],lepM[signal_lep],met_Et,met_Phi)")
    dpresel[key] = dpresel[key].Define("dphimetll","deltaPhi_metll(lepPt[signal_lep],lepEta[signal_lep],lepPhi[signal_lep],lepM[signal_lep],met_Et,met_Phi)")


# In[15]:


dCRZ = {}
dCRDib = {}
dCRTop = {}

for key in dfull.keys():
    dCRZ[key] = dpresel[key].Filter("mll_new > 110 && mll_new < 180","mll in 110-180")
    dCRZ[key] = dpresel[key].Filter("n_bjet == 0","bjet veto")
    dCRZ[key] = dpresel[key].Filter("METrel > 20 && METrel < 50","metrel in 20 - 50")
    
    dCRDib[key] = dpresel[key].Filter("mll_new > 70 && mll_new < 110","mll in 70-110")
    dCRDib[key] = dpresel[key].Filter("n_bjet == 0","bjet veto")
    dCRDib[key] = dpresel[key].Filter("met_Sign > 10","met sign. > 10")
    
    dCRTop[key] = dpresel[key].Filter("!isSF","DF")
    dCRTop[key] = dpresel[key].Filter("mll_new > 110 && mll_new < 180","mll in 110-180")
    dCRTop[key] = dpresel[key].Filter("METrel > 50","metrel > 50")


# In[16]:


frm_dic = {'dphimetll':{"xmin":0,"xmax":3.125,"nbins":25,"tit":"#Delta#Phi(MET,ll) [GeV]"},
            'mll_new':{"xmin":0,"xmax":1500,"nbins":300,"tit":"m_{ll}) [GeV]"},
           #'sglep_eta':{"xmin":-3,"xmax":3,"nbins":60,"tit":"#eta"},
           #'sgel_eta':{"xmin":-3,"xmax":3,"nbins":60,"tit":"#eta"},
           #'sgmu_eta':{"xmin":-3,"xmax":3,"nbins":60,"tit":"#eta"},
           #'sglep_flavor':{"xmin":0,"xmax":3,"nbins":3,"tit":"lep flavor"},
           'METrel':{"xmin":0,"xmax":1500,"nbins":300,"tit":"MET_{rel} [GeV]"},
           'JVT_fwdjet':{"xmin":0,"xmax":1,"nbins":10,"tit":"Fwd_jet JVT"},
          }


# In[17]:


hist = {}
allhistos = []
for key in dpresel.keys():
    if not key in hist.keys():
        hist[key] = {}
    for ch in ["isEE","isMM","!isSF"]:
        for reg in ["CRZ","CRTop","CRDib"]:
            newkey = ch+"_"+reg
            if not newkey in hist[key].keys():
                hist[key][newkey] = {}
            for var in frm_dic.keys():
                if reg == "CRZ" and ch is not "!isSF":
                    hist[key][newkey][var] = dCRZ[key].Filter(ch).Histo1D(("h_%s_%s_%s_%s"%(ch,key,var,reg),"",frm_dic[var]['nbins'],frm_dic[var]['xmin'],frm_dic[var]['xmax']),var,"weight")
                    allhistos.append(hist[key][newkey][var])
                elif reg == "CRDib" and ch is not "!isSF":
                    hist[key][newkey][var] = dCRDib[key].Filter(ch).Histo1D(("h_%s_%s_%s_%s"%(ch,key,var,reg),"",frm_dic[var]['nbins'],frm_dic[var]['xmin'],frm_dic[var]['xmax']),var,"weight")
                    allhistos.append(hist[key][newkey][var])
                elif reg == "CRTop" and ch == "!isSF":
                    hist[key][newkey][var] = dCRTop[key].Filter(ch).Histo1D(("h_%s_%s_%s_%s"%(ch,key,var,reg),"",frm_dic[var]['nbins'],frm_dic[var]['xmin'],frm_dic[var]['xmax']),var,"weight")
                    allhistos.append(hist[key][newkey][var])


# In[18]:


print("Number of histograms = %i"%len(allhistos))


# In[ ]:

R.RDF.RunGraphs(allhistos)
#get_ipython().run_cell_magic('time', '', 'R.RDF.RunGraphs(allhistos)')


# In[ ]:


#get_ipython().run_line_magic('tb', '')


# In[ ]:


ths = {}
datastack = {}
for key in dfull.keys():
    if not key in d_samp.keys():
        print("ERROR \t No info on sample %s"%key)
    for ch in ["isEE","isMM","!isSF"]:
        for reg in ["CRZ","CRTop","CRDib"]:
            newkey = ch+"_"+reg
            if not newkey in ths.keys():
                ths[newkey] = {}
                datastack[newkey] = {}
            for var in frm_dic.keys():
                if not var in hist[key][newkey]: continue
                if not var in ths[newkey].keys():
                    ths[newkey][var] = R.THStack()
                    datastack[newkey][var] = R.THStack()
                # To get the overflow into last bin
                hist[key][newkey][var].SetBinContent(hist[key][newkey][var].GetNbinsX(), hist[key][newkey][var].GetBinContent(hist[key][newkey][var].GetNbinsX()) + hist[key][newkey][var].GetBinContent(hist[key][newkey][var].GetNbinsX() + 1))
                if "data" in key:
                    hist[key][newkey][var].SetMarkerStyle(20)
                    datastack[newkey][var].Add(hist[key][newkey][var].GetValue())
                else:
                    hist[key][newkey][var].SetLineColor(d_samp[key]['f_color'])
                    hist[key][newkey][var].SetFillColor(d_samp[key]['f_color'])
                    ths[newkey][var].Add(hist[key][newkey][var].GetValue())


# In[ ]:


ch = "isMM_CRZ"
var = "mll_new"


can  = R.TCanvas('can','',900,700)
can.Draw()
pad1 = R.TPad('pad1', '', 0.0, 0.40, 1.0, 1.0)
pad2 = R.TPad('pad2', '', 0.0, 0.00, 1.0, 0.4)
pad1.SetTopMargin(0.1)
pad1.SetBottomMargin(0.005)
pad2.SetBottomMargin(0.3)
pad1.SetLeftMargin(0.1)
pad2.SetLeftMargin(0.1)
pad1.Draw()
pad2.Draw()
pad1.cd()
pad1.SetLogy()
ths[ch][var].Draw("hist")
if "CRZ" in ch:
    ths[ch][var].GetXaxis().SetRangeUser(110,180)
ths[ch][var].SetMinimum(1.0)
hsumMC = (ths[ch][var].GetStack().Last()).Clone("hsumMC")
hsumDATA = (datastack[ch][var].GetStack().Last()).Clone("hsumDATA")
hsumDATA.Draw("same e")
hratio = hsumMC.Clone("hratio")
hratio.Reset()
hratio.Divide(hsumDATA,hsumMC)
pad2.cd()
pad2.SetGridx()
pad2.SetGridy()
hratio.SetLineColor(R.kBlack)
hratio.SetMarkerColor(R.kBlack)
hratio.SetMarkerStyle(20)
hratio.SetMinimum(0.6)
hratio.SetMaximum(1.4)
if "CRZ" in ch:
    hratio.GetXaxis().SetRangeUser(110,180)
hratio.Draw("e")
hratio.GetXaxis().SetTitle("pT")
hratio.SetMaximum(1.8)
hratio.SetMinimum(0.4)


# In[ ]:


cols = ["DatasetNumber","weights"]


# In[ ]:


for c in dfull["Diboson"].GetColumnNames():
    print(c)


# In[ ]:


import pandas as pd
pandas = {}
for key in dpresel.keys():
    print(key)
    pandas[key] = pd.DataFrame(data=dpresel[key].Filter("dphimetll > 2.5","dphimetll > 2.5").AsNumpy(cols))


# In[ ]:


for key in pandas.keys():
    print(key)
    print("_"*50)
    print(pandas[key]["DatasetNumber"].value_counts())


# In[ ]:


4./1.25


# In[ ]:




