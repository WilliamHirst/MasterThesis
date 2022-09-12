#!/usr/bin/env python
# coding: utf-8

# In[1]:


import ROOT as R
from os import listdir
from os.path import isfile, join
import sys
sys.path.insert(1,"../plotting/")
from samples import configure_samples
R.EnableImplicitMT()
#get_ipython().run_line_magic('jsroot', 'on')


# In[2]:


d_samp, d_type, d_reg = configure_samples()


# In[3]:


mypath = "/storage/eirikgr/ANAntuples/PHYS_Bkgs_mc16e/"
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
df = {}
for of in onlyfiles:
    if not of.endswith("_merged_processed.root"): continue
    typ = of.replace("_merged_processed.root","")
    print(typ)
    df[typ] = R.RDataFrame("%s_NoSys"%typ,mypath+"/"+of)


# In[ ]:





# In[4]:


mypath = "/storage/eirikgr/ANAntuples/PHYS_Data/"
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
for of in onlyfiles:
    if not of.endswith("_merged_processed.root"): continue
    typ = of.replace("_merged_processed.root","")
    print(typ)
    df[typ] = R.RDataFrame("%s"%typ,mypath+"/"+of)


# In[5]:


mypath = "/scratch/eirikgr/ANAoutput/Tue_Aug_31_2021_02L_NTUP_data18_OLD2L2J/"
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
for of in onlyfiles:
    if not of.endswith("_merged_processed_COPY_0_0.root"): continue
    typ = "FNP"
    treename = of.replace("_merged_processed_COPY_0_0.root","")
    print(typ)
    df[typ] = R.RDataFrame(treename,mypath+"/"+of)


# In[6]:


#R.gROOT.ProcessLine(".L functions.cxx+");
R.gSystem.AddDynamicPath("-I/home/eirikgr/myNtupAnalysis/jupyter-notebooks/")
R.gInterpreter.Declare('#include "Cfunctions.h"') # Header with the definition of the myFilter function
R.gSystem.Load("Cfunctions.so") # Library with the myFilter function


# In[7]:


cols = df["PHWW"].GetColumnNames()
variables = []
blacklist = ["lepIFFClass", "lepEgMotherPdgId"]
plotdic = {}
for c in cols:
    typ = df["PHWW"].GetColumnType(c)
    if "bool" in typ or "Bool_t" in typ: continue
    #print(typ)
    cstr = str(c)
    if cstr in blacklist: continue
    if "lep" in cstr and not ("HLT" in cstr or "Weight" in cstr or "Iso" in 
                              cstr or "Truth" in cstr or "Origin" in cstr or "Type" in cstr
                             or "SF" in cstr):
        variables.append(c)
        print(c)
        plotdic[c] = {"nbin":0,"max":0,"min":0}


# In[8]:


for key in plotdic.keys():
    print(key)


# In[9]:


plotdic["lepFlavor"]["nbin"] = 2
plotdic["lepFlavor"]["nmax"] = 2
plotdic["lepFlavor"]["nmin"] = 0

plotdic["lepCharge"]["nbin"] = 3
plotdic["lepCharge"]["nmax"] = 1
plotdic["lepCharge"]["nmin"] = -1

plotdic["lepAuthor"]["nbin"] = 10
plotdic["lepAuthor"]["nmax"] = 10
plotdic["lepAuthor"]["nmin"] = 0

plotdic["lepPt"]["nbin"] = 100
plotdic["lepPt"]["nmax"] = 1000
plotdic["lepPt"]["nmin"] = 0

plotdic["lepEta"]["nbin"] = 80
plotdic["lepEta"]["nmax"] = 4
plotdic["lepEta"]["nmin"] = -4

plotdic["lepPhi"]["nbin"] = 80
plotdic["lepPhi"]["nmax"] = 4
plotdic["lepPhi"]["nmin"] = -4

plotdic["lepM"]["nbin"] = 100
plotdic["lepM"]["nmax"] = 1000
plotdic["lepM"]["nmin"] = 0

plotdic["lepD0"]["nbin"] = 200
plotdic["lepD0"]["nmax"] = 10
plotdic["lepD0"]["nmin"] = -10

plotdic["lepZ0"]["nbin"] = 200
plotdic["lepZ0"]["nmax"] = 10
plotdic["lepZ0"]["nmin"] = -10

plotdic["lepD0Sig"]["nbin"] = 100
plotdic["lepD0Sig"]["nmax"] = 10
plotdic["lepD0Sig"]["nmin"] = 0

plotdic["lepZ0SinTheta"]["nbin"] = 100
plotdic["lepZ0SinTheta"]["nmax"] = 10
plotdic["lepZ0SinTheta"]["nmin"] = 0

plotdic["lepECIDS"]["nbin"] = 2
plotdic["lepECIDS"]["nmax"] = 2
plotdic["lepECIDS"]["nmin"] = 0

plotdic["mymll"] = {}
plotdic["mymll"]["nbin"] = 200
plotdic["mymll"]["nmax"] = 2000
plotdic["mymll"]["nmin"] = 0

plotdic["mymt2"] = {}
plotdic["mymt2"]["nbin"] = 200
plotdic["mymt2"]["nmax"] = 2000
plotdic["mymt2"]["nmin"] = 0


# In[10]:


for key in df.keys():
    if "data" in key:
        df[key] = df[key].Define("wgt","1.0")
    elif "FNP" in key:
        df[key] = df[key].Define("wgt","FNP_WEIGHTS")
    else:
        df[key] = df[key].Define("scalelumi","(RandomRunNumber) < 320000 ? 36207.65 : (((RandomRunNumber) > 320000 && (RandomRunNumber) < 348000) ? 44307.4 : 58450.1);")
        df[key] = df[key].Define("wgt","((genWeight) * (eventWeight) * (leptonWeight) * (jvtWeight) * (bTagWeight) * (pileupWeight) * scalelumi)")
    df[key] = df[key].Define("baseline_el","lepPt > 9.0 && lepFlavor == 1 && abs(lepEta) < 2.47 && abs(lepEta) < 2.47 && lepLoose && lepPassBL && abs(lepZ0SinTheta) < 0.5")
    df[key] = df[key].Define("baseline_mu","lepPt > 9.0 && lepFlavor == 2 && abs(lepEta) < 2.6 && lepMedium && abs(lepZ0SinTheta) < 0.5")
    df[key] = df[key].Define("baseline_lep","baseline_el || baseline_mu")
    df[key] = df[key].Define("isOS","isOS(lepCharge[baseline_lep])")
    df[key] = df[key].Filter("Sum(baseline_lep) == 2","2 baseline leptons")
    df[key] = df[key].Filter("isOS","opposite sign")
    df[key] = df[key].Define("mymll","ComputeInvariantMass(lepPt[baseline_lep],lepEta[baseline_lep],lepPhi[baseline_lep],lepM[baseline_lep])")    
    df[key] = df[key].Define("mymt2","calcMT2(lepPt[baseline_lep],lepEta[baseline_lep],lepPhi[baseline_lep],lepM[baseline_lep], met_Et, met_Phi)")    


# In[11]:


histo = {}
allhisto = []
variables = ["mymll"]#,"mymt2"]
for v in variables:
    print(v)
    if not v in plotdic.keys():
        print("No information on binning for %s"%v)
        continue
    histo[v] = {}
    for key in df.keys():
        histo[v][key] = df[key].Histo1D(("h_%s_%s"%(v,key),"h_%s_%s"%(v,key),plotdic[v]["nbin"],
                                         plotdic[v]["nmin"],plotdic[v]["nmax"]),v,"wgt")
        allhisto.append(histo[v][key])


# In[12]:


print("Number of histograms = %i"%len(allhisto))
print("Number of runs = %i"%df["FNP"].GetNRuns())
print("Number of slots = %i"%df["FNP"].GetNSlots())
#df["FNP"].Describe()
#%%time
#%%time


# In[13]:


#graph = df["FNP"].SaveGraph(1)


# In[ ]:


#get_ipython().run_cell_magic('time', '', 'R.RDF.RunGraphs(allhisto)')


# In[ ]:


stack = {}
data = {}
legend = {}
for v in histo.keys():
    # Add legend
    legend[v] = R.TLegend(0.60, 0.60, 0.8, 0.85)
    legend[v].SetTextFont(42)
    legend[v].SetFillStyle(0)
    legend[v].SetBorderSize(0)
    legend[v].SetTextSize(0.04)
    stack[v] = R.THStack()
    data[v] = {}
    for key in histo[v].keys():
        if not key in d_samp.keys():
            print("skipping %s"%key)
            continue
        if not "data" in key:
            h = histo[v][key].GetValue()
            h.SetFillColor(d_samp[key]['f_color'])
            h.SetLineColor(d_samp[key]['f_color'])
            stack[v].Add(h)
            legend[v].AddEntry(h,"%-s"%key,"f")
        else:
            data[v][key] = histo[v][key].GetValue()
            data[v][key].SetMarkerStyle(d_samp[key]['m_type'])
            data[v][key].SetLineColor(d_samp[key]['f_color'])
            data[v][key].SetMarkerColor(d_samp[key]['f_color'])
            legend[v].AddEntry(data[v][key],"%-s"%key,"lp")


# In[ ]:


#get_ipython().run_cell_magic('time', '', 'datak = "data18"\nlumi = 58450.1\nsumMC = {}\nratio = {}\nv = "lepD0Sig"\n# Create canvas with pad\nc = R.TCanvas("c", "", 900, 700)\nc.Draw()\nR.gStyle.SetOptStat(0)\npad = R.TPad("upper_pad", "", 0, 0.2, 1, 1.0)\npad2 = R.TPad("lower_pad", "", 0, 0, 1, 0.2)\npad.SetTickx(False)\npad.SetTicky(False)\npad.SetBottomMargin(0.005)\npad.SetLogy()\npad.Draw()\npad2.Draw()\npad.cd()\nstack[v].Draw("HIST")\n\nsumMC[v] = stack[v].GetStack().Last().Clone("sumMC")\nratio[v] = sumMC[v].Clone("ratio")\nratio[v].Clear()\nratio[v].Divide(sumMC[v],data[v][datak])\n\n# Draw stack with MC contributions\nstack[v].GetXaxis().SetLabelSize(0.04)\nstack[v].GetXaxis().SetTitleSize(0.045)\nstack[v].GetXaxis().SetTitleOffset(1.3)\nstack[v].GetXaxis().SetTitle("m_{T}^{W#rightarrow l#nu} [GeV]")\nstack[v].GetYaxis().SetTitle("Events")\nstack[v].GetYaxis().SetLabelSize(0.04)\nstack[v].GetYaxis().SetTitleSize(0.045)\nstack[v].SetMaximum(1e7 * lumi/1000.)\nstack[v].SetMinimum(1)\n\n# Draw data\ndata[v][datak].SetMarkerStyle(20)\ndata[v][datak].SetMarkerSize(1.2)\ndata[v][datak].SetLineWidth(2)\ndata[v][datak].SetLineColor(R.kBlack)\ndata[v][datak].Draw("E SAME")\n\n# Draw legend\nlegend[v].Draw("SAME")\n\n# Add ATLAS label\ntext = R.TLatex()\ntext.SetNDC()\ntext.SetTextFont(72)\ntext.SetTextSize(0.045)\ntext.DrawLatex(0.21, 0.86, "ATLAS")\ntext.SetTextFont(42)\ntext.DrawLatex(0.21 + 0.09, 0.86, "Preliminary")\ntext.SetTextSize(0.04)\ntext.DrawLatex(0.21, 0.80, "#sqrt{{s}} = 13 TeV, {:.1f} fb^{{-1}}".format(lumi / 1000.0))\n\npad2.cd()\npad2.SetGridy()\npad2.SetTopMargin(0.01)\npad2.SetTickx(False)\npad2.SetTicky(False)\nratio[v].SetTitle("")\nratio[v].GetXaxis().SetLabelSize(0.15)\nratio[v].GetYaxis().SetLabelSize(0.15)\nratio[v].SetMaximum(2)\nratio[v].SetMinimum(2)\nratio[v].Draw("ep")')


# In[ ]:





# In[ ]:




