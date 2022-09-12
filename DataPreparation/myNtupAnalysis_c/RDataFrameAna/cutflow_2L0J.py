#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import ROOT as R
get_ipython().run_line_magic('jsroot', 'off')


# In[ ]:


frm_dic = {#'ptllboost':{"xmin":0,"xmax":1000,"nbins":100,"tit":"p_{T,boost}^{ll} [GeV]"},
           #'costhetastar':{"xmin":0,"xmax":1,"nbins":150,"tit":"cos#theta^{*}_{l,l}"},
           #'deltaphill':{"xmin":0,"xmax":4,"nbins":40,"tit":"#Delta#phi_{l,l}"},
           #'deltaphimetl':{"xmin":0,"xmax":4,"nbins":40,"tit":"#Delta#phi_{MET,l}"},
           #'mt2_new':{"xmin":0,"xmax":200,"nbins":20,"tit":"m_{T}^{2}"},
           'n_sg_jet':{"xmin":0,"xmax":10,"nbins":10,"tit":"N_{jets}"},
           #'n_nojvt_jet':{"xmin":0,"xmax":10,"nbins":10,"tit":"N_{jets (no JVT)}"},
           'n_b_jet':{"xmin":0,"xmax":10,"nbins":10,"tit":"N_{b-jets}"},
           #'jetJVT':{"xmin":0,"xmax":1,"nbins":50,"tit":"jet_{JVT}"},
           'jetPt':{"xmin":0,"xmax":1000,"nbins":100,"tit":"jet_{p_{T}}"},
           'jetEta':{"xmin":-4,"xmax":4,"nbins":80,"tit":"jet_{#eta}"},
           'met_Sign':{"xmin":0,"xmax":20,"nbins":40,"tit":"MET sign."},
           'met_Et':{"xmin":0,"xmax":1000,"nbins":100,"tit":"MET"}
          }


# In[ ]:


d = {}
rootdir = "/scratch/eirikgr/ntuples/"
d["r21"] = R.RDataFrame("Zjets_NoSys",rootdir+"/MC/SUSY2_Bkgs_mc16e/DEC20syst/merged/Zjets_merged_processed.root")
d["r22"] = R.RDataFrame("Zjets_NoSys",rootdir+"/PHYS_Bkgs_mc16e/JUN16/merged/Zjets_merged_processed.root")


# In[ ]:


weight = 'genWeight*eventWeight*leptonWeight*jvtWeight*bTagWeight*pileupWeight*globalDiLepTrigSF*58450.1'


# In[ ]:


dfull = {}
for key in d.keys():
    dfull[key] = d[key].Define("weight",weight)
    dfull[key] = dfull[key].Define("baseline_el","lepPt > 9.0 && lepFlavor == 1 && abs(lepEta) < 2.47 && abs(lepEta) < 2.47 && lepLoose && lepPassBL && abs(lepZ0SinTheta) < 0.5")
    dfull[key] = dfull[key].Define("baseline_mu","lepPt > 9.0 && lepFlavor == 2 && abs(lepEta) < 2.6 && lepMedium && abs(lepZ0SinTheta) < 0.5")
    dfull[key] = dfull[key].Define("signal_el","baseline_el && lepPt > 9 && lepTight && lepIsoFCLoose && abs(lepD0Sig) < 5")
    dfull[key] = dfull[key].Define("signal_mu","baseline_mu && lepPt > 9 && lepIsoFCLoose && abs(lepD0Sig) < 3")
    dfull[key] = dfull[key].Define("signal_jet","jetPt > 60 && abs(jetEta) < 2.4 && ((jetPt < 60 && jetJVT > 0.5) || jetPt >= 60)")
    dfull[key] = dfull[key].Define("signal_jet_noJVT","jetPt > 20 && abs(jetEta) < 2.4")
    dfull[key] = dfull[key].Define("b_jet","jetPt > 20 && abs(jetEta) < 2.4 && jetMV2c10 > 0.11")
    dfull[key] = dfull[key].Define("signal_lep","signal_el || signal_mu")


# In[ ]:


col = d["r22"].GetColumnNames()
for c in col:
    print(c)


# In[ ]:


for key in d.keys():
    dfull[key] = dfull[key].Define("n_bl_lep","Sum(baseline_mu)+Sum(baseline_el)")
    dfull[key] = dfull[key].Define("n_sg_lep","Sum(signal_lep)")
    dfull[key] = dfull[key].Define("n_sg_jet","Sum(signal_jet)")
    dfull[key] = dfull[key].Define("n_b_jet","Sum(b_jet)")
    dfull[key] = dfull[key].Define("n_nojvt_jet","Sum(signal_jet_noJVT)")
    


# In[ ]:


#cols = R.vector('string')()
#cols.push_back('beamSpotWeight')
#p = dfull['r22'].Display(cols)
#p.Print()


# In[ ]:


var = "n_sg_jet"
hist = {}
for key in dfull.keys():
    hist[key] = dfull[key].Histo1D(("h_%s"%var,"",frm_dic[var]['nbins'],frm_dic[var]['xmin'],frm_dic[var]['xmax']),var,"weight")


# In[ ]:


c = R.TCanvas()
c.Draw()
hist["r22"].Draw("hist")


# In[ ]:


R.gInterpreter.Declare(
"""
using Vec_t = const ROOT::VecOps::RVec<int>;
bool isOS(Vec_t& chlep) {
    if(chlep[0]*chlep[1] < 0)return kTRUE;
    return kFALSE;
}
""")


# In[ ]:


R.gInterpreter.Declare(
"""
using Vec_t = const ROOT::VecOps::RVec<int>;
bool isSF(Vec_t& fllep) {
    if(fllep[0] == fllep[1])return kTRUE;
    return kFALSE;
}
""")


# In[ ]:


R.gInterpreter.Declare(
"""
using Vec2_t = const ROOT::VecOps::RVec<float>;
float ComputeInvariantMass(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e) {
    ROOT::Math::PtEtaPhiMVector p1(pt[0], eta[0], phi[0], e[0]);
    ROOT::Math::PtEtaPhiMVector p2(pt[1], eta[1], phi[1], e[1]);
    return (p1 + p2).mass();
}
""")


# In[ ]:


R.gInterpreter.AddIncludePath("../NTUPana_SusysSkim/CalcGenericMT2/CalcGenericMT2/");
R.gInterpreter.Declare(
"""
#include "MT2_ROOT.h"
using Vec2_t = const ROOT::VecOps::RVec<float>;
float calcMT2(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e, Float_t met_et, Float_t met_phi) {

    TLorentzVector p1;
    TLorentzVector p2;
    TLorentzVector met;
    p1.SetPtEtaPhiM(pt[0], eta[0], phi[0], e[0]);
    p2.SetPtEtaPhiM(pt[1], eta[1], phi[1], e[1]);
    met.SetPtEtaPhiM(met_et, 0.0, met_phi, 0.0);
    return ComputeMT2(p1,p2,met,0.,0.).Compute();
}
""")


# In[ ]:


R.gInterpreter.Declare(
"""
using Vec2_t = const ROOT::VecOps::RVec<float>;
float ptllboost(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e, Float_t met_et, Float_t met_phi) {

    TLorentzVector p1;
    TLorentzVector p2;
    TLorentzVector met;
    p1.SetPtEtaPhiM(pt[0], eta[0], phi[0], e[0]);
    p2.SetPtEtaPhiM(pt[1], eta[1], phi[1], e[1]);
    met.SetPtEtaPhiM(met_et, 0.0, met_phi, 0.0);
    return (met+p1+p2).Pt();
}
""")


# In[ ]:


R.gInterpreter.Declare(
"""
using Vec2_t = const ROOT::VecOps::RVec<float>;
float costhetastar(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e) {

    TLorentzVector p1;
    TLorentzVector p2;
    p1.SetPtEtaPhiM(pt[0], eta[0], phi[0], e[0]);
    p2.SetPtEtaPhiM(pt[1], eta[1], phi[1], e[1]);
    return TMath::ATan(fabs(p1.Eta()-p2.Eta())/2.);
}
""")


# In[ ]:


R.gInterpreter.Declare(
"""
using Vec2_t = const ROOT::VecOps::RVec<float>;
float deltaPhi_ll(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e) {

    TLorentzVector p1;
    TLorentzVector p2;
    p1.SetPtEtaPhiM(pt[0], eta[0], phi[0], e[0]);
    p2.SetPtEtaPhiM(pt[1], eta[1], phi[1], e[1]);
    return p1.DeltaPhi(p2);
}
""")


# In[ ]:


R.gInterpreter.Declare(
"""
using Vec2_t = const ROOT::VecOps::RVec<float>;
float deltaPhi_metl(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e, Float_t met_et, Float_t met_phi) {

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
""")


# In[ ]:


R.gInterpreter.Declare(
"""
using Vec2_t = const ROOT::VecOps::RVec<float>;
bool checkPt(Vec2_t& pt, float cut1, float cut2){
    if((pt[0] > cut1 && pt[1] > cut2) || (pt[1] > cut1 && pt[0] > cut2))return kTRUE;
    return kFALSE;
}
""")


# In[ ]:


for key in d.keys():
    dfull[key] = dfull[key].Filter("(n_bl_lep == 2)","2 baseline leptons")
    dfull[key] = dfull[key].Filter("(n_sg_lep == 2)","2 signal leptons")


# In[ ]:


for key in d.keys():
    dfull[key] = dfull[key].Filter("checkPt(lepPt[signal_lep],27,9)","pT > 27/9 GeV")


# In[ ]:


for key in d.keys():
    dfull[key] = dfull[key].Define("isOS","isOS(lepCharge[signal_lep])")
    dfull[key] = dfull[key].Define("isSF","isSF(lepFlavor[signal_lep])")


# In[ ]:


for key in d.keys():
    dfull[key] = dfull[key].Filter("trigMatch_1LTrigOR","matched to trigger")
    #dfull[key] = dfull[key].Filter("isSF","Same flavour")
    dfull[key] = dfull[key].Filter("isOS","Opposite sign")


# In[ ]:


for key in d.keys():
    dfull[key] = dfull[key].Define("mll_new","ComputeInvariantMass(lepPt[signal_lep],lepEta[signal_lep],lepPhi[signal_lep],lepM[signal_lep])")


# In[ ]:


for key in d.keys():
    dfull[key] = dfull[key].Filter("mll_new > 11","mll > 11 GeV")


# In[ ]:


for key in d.keys():
    dfull[key] = dfull[key].Filter("n_sg_jet < 2","N(jet) < 2")
    dfull[key] = dfull[key].Filter("n_b_jet == 0","N(bjet) == 0")
    dfull[key] = dfull[key].Filter("met_Sign > 3","MET sign. > 3")
    dfull[key] = dfull[key].Filter("((isSF && abs(mll_new-91.1876) > 15) || !isSF)","Z-veto (only SF)")


# In[ ]:


for key in d.keys():
    dfull[key] = dfull[key].Define("ptllboost","ptllboost(lepPt[signal_lep],lepEta[signal_lep],lepPhi[signal_lep],lepM[signal_lep],met_Et,met_Phi)")
    dfull[key] = dfull[key].Define("costhetastar","costhetastar(lepPt[signal_lep],lepEta[signal_lep],lepPhi[signal_lep],lepM[signal_lep])")
    dfull[key] = dfull[key].Define("deltaphill","deltaPhi_ll(lepPt[signal_lep],lepEta[signal_lep],lepPhi[signal_lep],lepM[signal_lep])")
    dfull[key] = dfull[key].Define("deltaphimetl","deltaPhi_metl(lepPt[signal_lep],lepEta[signal_lep],lepPhi[signal_lep],lepM[signal_lep],met_Et,met_Phi)")
    dfull[key] = dfull[key].Define("mt2_new","calcMT2(lepPt[signal_lep],lepEta[signal_lep],lepPhi[signal_lep],lepM[signal_lep],met_Et,met_Phi)")


# In[ ]:


presel = {}
for key in d.keys():
    presel[key] = dfull[key]


# In[ ]:


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
    


# # Direct Slepton Signal regions

# In[ ]:


for key in d.keys():
    dfull[key] = dfull[key].Filter("n_sg_jet == 0","jet-veto")
    dfull[key] = dfull[key].Filter("checkPt(lepPt[signal_lep],140,20)","pT > 140/20 GeV")
    dfull[key] = dfull[key].Filter("met_Sign > 7","MET sign. > 7")
    dfull[key] = dfull[key].Filter("ptllboost < 5","p_{T,boost}^{ll} < 5 GeV")
    dfull[key] = dfull[key].Filter("costhetastar < 0.2","cos#theta^{*}_{ll} < 0.2")
    dfull[key] = dfull[key].Filter("deltaphill > 2.2","#Delta#phi_{l,l} > 2.2")
    dfull[key] = dfull[key].Filter("deltaphimetl > 2.2","#Delta#phi_{MET,l} > 2.2")


# In[ ]:


allcuts = {}
for key in d.keys():
    allcuts[key] = dfull[key].Report()


# In[ ]:


print("RELEASE 22")
convertRDFCutflowToTex(allcuts['r22'],allcuts['r21'])
allcuts['r22'].Print()


# In[ ]:


print("RELEASE 21")
allcuts['r21'].Print()


# In[ ]:





# In[ ]:


h = {"r21":{},"r22":{}}
for key in d.keys():
    for var in frm_dic.keys():
        #if not "jet" in var: continue# and not "met" in var: continue
        h[key][var] = presel[key].Histo1D(("h_%s"%var,"",frm_dic[var]['nbins'],frm_dic[var]['xmin'],frm_dic[var]['xmax']),var,"weight")
        h[key][var].SetLineColor(R.kPink+5 if key == 'r22' else R.kTeal+5)
        h[key][var].SetFillColor(R.kPink+5 if key == 'r22' else R.kTeal+5)
        h[key][var].SetFillStyle(3004 if key == 'r22' else 3005)
        h[key][var].SetLineWidth(2)


# In[ ]:


cols = R.vector('string')()
cols.push_back('ptllboost')
cols.push_back('costhetastar')
cols.push_back('deltaphill')
cols.push_back('deltaphimetl')
p = dfull['r22'].Display(cols)


# In[ ]:


p.Print()


# In[ ]:


def customise_gPad(top=0.03, bot=0.15, left=0.17, right=0.08):

    R.gPad.Update()

    R.gStyle.SetTitleFontSize(0.0)

    # gPad margins
    R.gPad.SetTopMargin(top)
    R.gPad.SetBottomMargin(bot)
    R.gPad.SetLeftMargin(left)
    R.gPad.SetRightMargin(right)
    
    R.gStyle.SetOptStat(0) # Hide usual stats box
    
    R.gPad.Update()


# In[ ]:


def customise_axes(hist, xtitle, ytitle, scaleFactor=1.1, IsLogY=False, enlargeYaxis=False, scaling=False):

    # Set a universal text size
    text_size = 45

    R.TGaxis.SetMaxDigits(4)

    ##################################
    # X axis
    xax = hist.GetXaxis()

    # Precision 3 Helvetica (specify label size in pixels)
    xax.SetLabelFont(43)
    xax.SetTitleFont(43)
    # xax.SetTitleFont(13) # times
    
    xax.SetTitle(xtitle)
    xax.SetTitleSize(text_size)

    # Top panel
    if 'Events' in ytitle:
        xax.SetLabelSize(0)
        xax.SetLabelOffset(0.02)
        xax.SetTitleOffset(2.0)
        xax.SetTickSize(0.04)
    # Bottom panel
    else:
        xax.SetLabelSize(text_size - 7)
        xax.SetLabelOffset(0.03)
        xax.SetTitleOffset(3.5)
        xax.SetTickSize(0.08)

    # xax.SetRangeUser(0,2000)
    # xax.SetNdivisions(-505)

    R.gPad.SetTickx()

    ##################################
    # Y axis
    yax = hist.GetYaxis()

    # Precision 3 Helvetica (specify label size in pixels)
    yax.SetLabelFont(43)
    yax.SetTitleFont(43)

    yax.SetTitle(ytitle)
    yax.SetTitleSize(text_size)
    yax.SetTitleOffset(1.8)
    
    yax.SetLabelOffset(0.015)
    yax.SetLabelSize(text_size - 7)
    
    ymax = hist.GetMaximum()
    ymin = hist.GetMinimum()

    # if ymin == 0.0:
    #     print 'ymin = 0.0'

    # Top events panel
    if 'Events' in ytitle:
        yax.SetNdivisions(505)
        if IsLogY:
            if enlargeYaxis:
                ymax = 2 * 10 ** 10
                ymin = 0.00001
            else:
                # ymax = 3 * 10 ** 4
                # ymin = 0.5
                ymax = 3 * 10 ** 3
                ymin = 0.00001
            #if scaling:
            #    hist.SetMaximum(1.0)
            #else:
            hist.SetMaximum(ymax)
            hist.SetMinimum(ymin)
        else:
            #if scaling:
            #    hist.SetMaximum(1.0)
            #else:
            hist.SetMaximum(ymax*scaleFactor)
            hist.SetMinimum(0.0)
    # Bottom panel
    elif 'Ratio' in ytitle:
        yax.SetNdivisions(505)
        # Dynamic 
        if ymax*scaleFactor > 5:
            hist.SetMaximum(5)
        else: hist.SetMaximum(ymax*scaleFactor)
        if ymin*0.9 < -1:
            hist.SetMinimum(-2)#ymin*0.9)
        else: hist.SetMinimum(ymin*0.9)
        # Fixed
        #hist.SetMinimum(-0.5) 
        #hist.SetMaximum(2.5)  

    R.gPad.SetTicky()

    R.gPad.Update()


# In[ ]:





# In[ ]:


h["r21"]['ptllboost'].GetValue()


# In[ ]:


print(frm_dic.keys())


# In[ ]:


var = "jetEta"#"ptllboost"


# In[ ]:


h3 = {}
for var in h['r21'].keys():
    
    xtitle = frm_dic[var]['tit']
    ytitle = "Events"
    IsLogY = False
    enlargeYaxis = False
    scaling = True
    uselogY = True
    
    filepre = "presel_ptjet60"
    
    if scaling:
        integral1 = h["r21"][var].Integral()
        scale1 = 1. / integral1
        h["r21"][var].Scale(scale1) 

        integral2 = h["r22"][var].Integral()
        scale2 = 1. / integral2
        h["r22"][var].Scale(scale2) 

    # Third histogram (ratio)
    h3[var] = h["r22"][var].GetValue().Clone("h3_%s"%var)
    h3[var].Sumw2() #
    h3[var].Divide(h["r21"][var].GetValue())
    h3[var].SetLineColor(R.kGray+2)
    h3[var].SetLineWidth(2)
    h3[var].SetMarkerStyle(21)

    can  = R.TCanvas('','',1000,1000)
    customise_gPad()
    pad1 = R.TPad('pad1', '', 0.0, 0.40, 1.0, 1.0)
    pad2 = R.TPad('pad2', '', 0.0, 0.00, 1.0, 0.4)
    # Margings used for the pads
    gpLeft = 0.17
    gpRight = 0.05

    leg = R.TLegend(0.80,0.90,0.85,0.95)
    leg.SetBorderSize(0)
    #leg.SetTextSize(0.033)
    #leg.SetTextFont(42) # Helvetica
    leg.SetTextSize(15)
    leg.SetTextFont(43) # Helvetica
    leg.SetNColumns(1)
    leg.AddEntry(h["r21"][var].GetValue(),"REL21",'lf')
    leg.AddEntry(h["r22"][var].GetValue(),"REL22",'lf')

    #-------
    # PAD1
    #-------
    if uselogY:
        pad1.SetLogy()
    pad1.Draw()
    pad1.cd()
    customise_gPad(top=0.08, bot=0.04, left=gpLeft, right=gpRight)
    can.Draw()
    if h["r21"][var].GetMaximum() > h["r22"][var].GetMaximum():
        customise_axes(h["r21"][var], xtitle, ytitle, 1.1, IsLogY, enlargeYaxis, scaling) # Need to feed in the one that is drawn first
        if uselogY:
            h["r21"][var].GetYaxis().SetRangeUser(0.001,h["r21"][var].GetMaximum())
            h["r22"][var].GetYaxis().SetRangeUser(0.001,h["r21"][var].GetMaximum())
        h["r21"][var].Draw("hist")
        h["r22"][var].Draw("hist same")
    elif h["r21"][var].GetMaximum() < h["r22"][var].GetMaximum():
        customise_axes(h["r22"][var], xtitle, ytitle, 1.1, IsLogY, enlargeYaxis, scaling) 
        h["r21"][var].GetYaxis().SetRangeUser(0.001,h["r22"][var].GetMaximum())
        h["r22"][var].GetYaxis().SetRangeUser(0.001,h["r22"][var].GetMaximum())
        h["r22"][var].Draw("hist")
        h["r21"][var].Draw("hist same")
    

    can.cd()

    leg.Draw()

    pad2.Draw()
    pad2.cd()
    customise_gPad(top=0.05, bot=0.39, left=gpLeft, right=gpRight)
    #customise_gPad(top=0, bot=0.39, left=gpLeft, right=gpRight) # joins upper and lower plot
    pad2.SetGridy()
    h3[var].Draw("ep")
    customise_axes(h3[var], xtitle, "Ratio", 1.1, IsLogY, enlargeYaxis, scaling == 'True')
    can.SaveAs("/mn/felt/u1/eirikgr/talks/SUSYBkgf/%s_%s_%s.png"%(var,"scaled" if scaling else "noscale",filepre))


# In[ ]:





# In[ ]:




