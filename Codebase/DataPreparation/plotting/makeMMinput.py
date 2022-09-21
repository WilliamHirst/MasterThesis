from ROOT import *
import sys, os
from os import listdir
from os.path import isfile, join
import array
from math import sqrt
import socket
from samples import configure_samples
import copy


TH1.SetDefaultSumw2()
TH2.SetDefaultSumw2()

gROOT.SetBatch(False)

d_samp,d_type, d_reg = configure_samples()

if socket.gethostname() == "quark":
  astyle = "/home/eirikgr/atlasrootstyle/"#atlasstyle-00-04-02/"#
elif socket.gethostname() == "hepp01.hpc.uio.no":
  astyle = "/home/eirikgr/atlasrootstyle/"
else:
  astyle = "/mn/felt/u1/eirikgr/bin/atlasstyle-00-04-02/"                                                                                                                                                                                                                                                                                                                                    
gROOT.SetMacroPath(astyle)                                                                                                                                                                                     
gROOT.LoadMacro("AtlasUtils.C")                                                                                                                                                                                
gROOT.LoadMacro("AtlasStyle.C")                                                                                                                                                                                
gROOT.LoadMacro("AtlasLabels.C")                                                                                                                                                                               
SetAtlasStyle()  
gStyle.SetErrorX()
gStyle.SetEndErrorSize() 

gStyle.SetPaintTextFormat(".2f");
#gROOT.ProcessLine("gErrorIgnoreLevel = kFatal")
def getxsecerr(xsec_up,xsec_dw, xsec_nom):

    error_up_x = array.array('f')
    error_low_x = array.array('f')
    error_up_y = array.array('f')
    error_low_y = array.array('f')
    x = array.array('f')
    y = array.array('f')
    for ibin in range(1,xsec_up.GetNbinsX()+1):

        x.append(xsec_nom.GetXaxis().GetBinCenter(ibin))
        y.append(xsec_nom.GetBinContent(ibin))

        error_up_x.append(xsec_up.GetXaxis().GetBinUpEdge(ibin)-xsec_up.GetXaxis().GetBinCenter(ibin))
        error_low_x.append(xsec_dw.GetXaxis().GetBinCenter(ibin)-xsec_dw.GetXaxis().GetBinLowEdge(ibin))

            
        xsecdw = xsec_dw.GetBinContent(ibin)
        if xsec_dw.GetBinContent(ibin) < xsec_nom.GetBinContent(ibin):
            xsecdw = 0

        xsecup = xsec_up.GetBinContent(ibin)
        if xsec_up.GetBinContent(ibin) > xsec_nom.GetBinContent(ibin):
            xsecup = 0

        error_up_y.append(xsecdw-xsec_nom.GetBinContent(ibin))
        error_low_y.append(xsec_nom.GetBinContent(ibin)-xsecup)
    return x,y,error_up_x,error_up_y,error_low_x,error_low_y

def getEff(hTD,hLD,vb = 0):
    is2D = False
    hT = hTD.Clone(hTD.GetName()+"_getEFF")
    hL = hLD.Clone(hLD.GetName()+"_getEFF")
    if type(hT) is TH2F: 
        if vb: print("DEBUG \t Histogram is 2D!")
        is2D = True
    elif type(hT) is int:
        print("WARNING \t Histogram is not a histogram")
        return hT
    if vb: print("Dividing %s on %s"%(hTD.GetName(),hLD.GetName()))
    if not is2D and vb:
        for i in range(0,hT.GetNbinsX()+2):
          print("BEFORE: Bin [%03d-%03d] : %.2f +/- %.2f %.2f +/- %.2f %.2f"%(hT.GetXaxis().GetBinLowEdge(i+1),hT.GetXaxis().GetBinUpEdge(i+1),
                                                                              hT.GetBinContent(i+1),hT.GetBinError(i+1),hL.GetBinContent(i+1),hL.GetBinError(i+1),
                                                                              hT.GetBinContent(i+1)/hL.GetBinContent(i+1) if hL.GetBinContent(i+1) else 0.0))
    h_eff = hTD.Clone(hTD.GetName().replace("nT","eff"))
    h_eff.Divide(hT,hL)
    h_eff.SetDirectory(0)
    if is2D:
        asym = TEfficiency(hT,hL)
        asym.SetStatisticOption(TEfficiency.kBBayesian)
        asym.SetConfidenceLevel(0.68);
        asym.SetBetaAlpha(1.)
        asym.SetBetaBeta(1.) 
        ni = 0
        for i in range(1,h_eff.GetNbinsX()+1):
            for j in range(1,h_eff.GetNbinsY()+1):
                ibin = h_eff.GetBin(i,j)
                if h_eff.GetBinContent(ibin)>1.0 and ("mc" in hT.GetName() or "data" in hT.GetName()):
                    print("ERROR \t Be careful, r = %.4f for %s. Errors may be wrong!"%(h_eff.GetBinContent(ibin),hT.GetName()))
                h_eff.SetBinError(ibin,(asym.GetEfficiencyErrorLow(ibin)+asym.GetEfficiencyErrorUp(ibin))/2.0)
                if h_eff.GetBinError(ibin)>0.3 and ("mc" in hT.GetName() or "data" in hT.GetName()):
                  print("ERROR \t Be careful, err = %.4f for %s. Errors may be wrong!"%(h_eff.GetBinError(ibin),hT.GetName()))
                
    else:
        asym = TGraphAsymmErrors();
        asym.Divide(hT,hL,"cl=0.683 b(1,1) mode")
        ni = 0
        for i in range(0,h_eff.GetNbinsX()):
            if (h_eff.GetBinContent(i)>1.0 or h_eff.GetBinContent(i)<0.0) and ("mc" in hT.GetName() or "data" in hT.GetName()):
                print("ERROR \t Be careful, f = %.4f for %s. Errors may be wrong!"%(h_eff.GetBinContent(i),hT.GetName()))
            if h_eff.GetBinContent(i+1):
                h_eff.SetBinError(i+1,asym.GetErrorY(ni))
                ni += 1
                if vb: print("Bin [%03d-%03d] : %.2f +/- %.2f"%(h_eff.GetXaxis().GetBinLowEdge(i+1),h_eff.GetXaxis().GetBinUpEdge(i+1),h_eff.GetBinContent(i+1),asym.GetErrorY(ni)))
            else:
                h_eff.SetBinError(i+1,0.0)
            if (h_eff.GetBinError(i)>0.3) and ("mc" in hT.GetName() or "data" in hT.GetName()):
              print("ERROR \t Be careful, err = %.4f for %s. Errors may be wrong!"%(h_eff.GetBinError(i),hT.GetName()))
    return h_eff

def doSubt(hT,hL,hT_mc,hL_mc,syst = 1.0,vb = 0):
    # Need to clone, if not changes are visible outside function!
    hTmc = hT_mc.Clone("hTmc")
    hLmc = hL_mc.Clone("hLmc")
    if syst != 1.0:
      if vb: print("INFO \t Scaling MC with %.2f"%syst)
      hTmc.Scale(syst)
      hLmc.Scale(syst)
      hTd = hT.Clone(hT.GetName()+"_%.2f_subt"%syst)
      hLd = hL.Clone(hL.GetName()+"_%.2f_subt"%syst)
    else:
      hTd = hT.Clone(hT.GetName()+"_subt")
      hLd = hL.Clone(hL.GetName()+"_subt")
    
    if vb: print("Subtracting %s from %s" %(hT.GetName(),hT_mc.GetName()))
    #print("Before : %10s %10s"%("Tight","Loose"))
    before = {"T":[],"L":[]}
    for ibin in range(1,hTd.GetNbinsX()+1):
        #print("Bin %02d : %10.0f %10.0f"%(ibin,hTd.GetBinContent(ibin),hLd.GetBinContent(ibin)))
        before["T"].append(hTd.GetBinContent(ibin))
        before["L"].append(hLd.GetBinContent(ibin))
    hTd.Add(hTmc,-1.0)
    hLd.Add(hLmc,-1.0)

    # Finxing errors after subtraction!
    for ibin in range(1,hTd.GetNbinsX()+1):
      err = (hTd.GetBinError(ibin)*hTd.GetBinError(ibin))-(hTmc.GetBinError(ibin)*hTmc.GetBinError(ibin))
      if vb:
        print("Data = %.2f +/- %.2f"%(hTd.GetBinContent(ibin),hTd.GetBinError(ibin)))
        print("MC   = %.2f +/- %.2f"%(hTmc.GetBinContent(ibin),hTmc.GetBinError(ibin)))
        print("T: Error is %s"%sqrt(err))
      if err > 0:
        hTd.SetBinError(ibin,sqrt(err))
      else:
        hTd.SetBinError(ibin,0.0)
      err = (hLd.GetBinError(ibin)*hLd.GetBinError(ibin))-(hLmc.GetBinError(ibin)*hLmc.GetBinError(ibin))
      if vb: print("L: Error is %s"%err)
      if err > 0:
        hLd.SetBinError(ibin,sqrt(err))
      else:
        hLd.SetBinError(ibin,0.0)
    
    if vb: print("Bin           : %21s %21s %4s %4s"%("Before","After", "%%", "%%"))
    for ibin in range(1,hTd.GetNbinsX()+1):
        if not (before["L"][ibin-1] and before["T"][ibin-1]): continue
        if vb: print("Bin [%03d-%03d] : %10.0f %10.0f %10.0f %10.0f %4.0f %4.0f"%(hTd.GetXaxis().GetBinLowEdge(ibin),
                                                                                  hTd.GetXaxis().GetBinUpEdge(ibin),
                                                                                  before["T"][ibin-1],before["L"][ibin-1],hTd.GetBinContent(ibin),
                                                                                  hLd.GetBinContent(ibin),((before["T"][ibin-1]-hTd.GetBinContent(ibin))/before["T"][ibin-1])*100.,
                                                                                  ((before["L"][ibin-1]-hLd.GetBinContent(ibin))/before["L"][ibin-1])*100.))

    if vb: print("Returning %s and %s" %(hTd.GetName(),hLd.GetName()))
    hTd.SetDirectory(0)
    hLd.SetDirectory(0)
    return hTd,hLd

def doRebin(h,xbins=array.array('f',[0,10,20,25,30,40,50,60,80,100,1000]),ybins=array.array('f',[0,10,20,25,30,40,50,60,80,100,1000]),vb = 0):
    if type(h) is TH2F:
        hreb = TH2F(h.GetName()+"_REB",h.GetTitle()+"_REB",len(xbins)-1,xbins,len(ybins)-1,ybins)
        xaxis = h.GetXaxis()
        yaxis = h.GetYaxis()
        for j in range(1,yaxis.GetNbins()+1):
            xfill = -1
            for i in range(1,xaxis.GetNbins()+1):
                binc = h.GetBinContent(i,j)
                if binc < -1:
                  binc = 0.0
                yfill = -1
                #hreb.Fill(xaxis.GetBinCenter(i),yaxis.GetBinCenter(j),binc)
                if(xaxis.GetBinCenter(i)>=xbins[-1]): 
                    xfill = xbins[-1]-1
                    #print("Outside x : Bin [%.0f-%.2f]"%(xaxis.GetBinCenter(i),yaxis.GetBinCenter(j)))
                if(abs(yaxis.GetBinCenter(j))>=ybins[-1]): 
                    yfill = ybins[-1]-0.1
                    #print("Outside y : Bin [%.0f-%.2f]"%(xaxis.GetBinCenter(i),yaxis.GetBinCenter(j)))
                if xfill >= 0 and yfill >= 0:
                    hreb.Fill(xfill,yfill,binc)
                elif yfill >= 0:
                    hreb.Fill(xaxis.GetBinCenter(i),yfill,binc)
                elif xfill >= 0:
                    hreb.Fill(xfill,abs(yaxis.GetBinCenter(j)),binc)
                else:
                    hreb.Fill(xaxis.GetBinCenter(i),abs(yaxis.GetBinCenter(j)),binc)
        for i in range(hreb.GetNbinsX()+1):
            for j in range(hreb.GetNbinsY()+1):
                if hreb.GetBinContent(i,j) > 0:
                    hreb.SetBinError(i,j,sqrt(hreb.GetBinContent(i,j)))
                else:
                    hreb.SetBinError(i,j,0)
    else:
        hreb = TH1F(h.GetName()+"_REB",h.GetTitle()+"_REB",len(xbins)-1,xbins)
        for nbin in range(1,h.GetNbinsX()+1):
            if vb:
              print("nbin = %i, nbins = %i"%(nbin,h.GetNbinsX()))
              print("Low edge = %i, Adding %.2f into bin %i"%(h.GetBinLowEdge(nbin),h.GetBinContent(nbin),hreb.GetBinLowEdge(hreb.FindBin(abs(h.GetBinLowEdge(nbin))))))
            hreb.AddBinContent(hreb.FindBin(abs(h.GetBinLowEdge(nbin))),h.GetBinContent(nbin))
        for nbin in range(1,hreb.GetNbinsX()+1):
            hreb.SetBinError(nbin,sqrt(hreb.GetBinContent(nbin)) if hreb.GetBinContent(nbin) >= 0 else 0.0) 

        # Fix overflow
        nbin = h.GetNbinsX()+1
        if vb: print("nbin = %i, nbins = %i"%(nbin,h.GetNbinsX()))
        hreb.AddBinContent(hreb.GetNbinsX(),h.GetBinContent(nbin))
        hreb.SetBinError(nbin,sqrt(hreb.GetBinContent(nbin)) if hreb.GetBinContent(nbin) >= 0 else 0.0)

    return hreb


vb = 0

trigger_setup = {"e_med":{"triggers":["e24_lhmedium_L1EM20VH","e24_lhmedium_nod0_L1EM20VH","e60_lhmedium_nod0"],"col":kOrange-4,"marker":22},
                 "e_loose":{"triggers":["e120_lhloose","e140_lhloose_nod0"],"col":kYellow-7,"marker":23},
                 "e_tight_iso":{"triggers":["e24_lhtight_nod0_ivarloose","e26_lhtight_nod0_ivarloose"],"col":kRed-3,"marker":21},
                 "m_med":{"triggers":["mu26_ivarmedium"],"col":kOrange-4,"marker":22},
                 "m_loose":{"triggers":["mu20_iloose_L1MU15"],"col":kYellow-7,"marker":23},
                 "m_no":{"triggers":["mu50"],"col":kGreen-7,"marker":21}}





trigger_categories = {"2e12_lhloose_L12EM10VH":"e_lhoose",
                      "2e17_lhvloose_nod0":"e_lhvloose",
                      "2e24_lhvloose_nod0":"e_lhvloose",
                      "2e17_lhvloose_nod0_L12EM15VHI":"e_lhvloose",
                      "mu26_imedium":"m_imed",
                      "mu26_ivarmedium":"m_ivarmed",
                      "mu50":"m_no",
                      "e17_lhloose_mu14":"em_lhloose",
                      "e26_lhmedium_nod0_mu8noL1":"em_lhmed"
                      
                      
}
plot_categories = True

'''
trigger_setup = {#"e_med":{"triggers":["e24_lhmedium_L1EM20VH","e24_lhmedium_nod0_L1EM20VH","e60_lhmedium_nod0"],"col":kOrange-4,"marker":22},
  "e_lhlooseL1":{"triggers":["2e12_lhloose_L12EM10VH"],"col":kYellow-7,"marker":23},
  "e_lhvloose":{"triggers":["2e17_lhvloose_nod0","2e24_lhvloose_nod0"],"col":kYellow-7,"marker":23},
  "e_lhvlooseL1":{"triggers":["2e17_lhvloose_nod0_L12EM15VHI"],"col":kYellow-7,"marker":23},
  #"e_tight_iso":{"triggers":["e24_lhtight_nod0_ivarloose","e26_lhtight_nod0_ivarloose"],"col":kRed-3,"marker":21},
  "m_imed":{"triggers":["mu26_imedium"],"col":kOrange-4,"marker":22},
  "m_ivarmed":{"triggers":["mu26_ivarmedium"],"col":kOrange-4,"marker":22},
  #"m_loose":{"triggers":["mu20_iloose_L1MU15"],"col":kYellow-7,"marker":23},
  "m_no":{"triggers":["mu50"],"col":kGreen-7,"marker":21}}
'''



bkg_suffix = ""
data_suffix = ""
# d_samp = {'ttbar'    :{'type':'bkg', 'leg':'t#bar{t}',                  'f_color':kRed-4,   'path':'ttbar'+bkg_suffix},
#           'singleTop':{'type':'bkg', 'leg':'Single top',                'f_color':kYellow-4,  'path':'singleTop'+bkg_suffix},
#           'Zjets'    :{'type':'bkg', 'leg':'Z+jets',                    'f_color':kOrange+1,  'path':'Zjets'+bkg_suffix},
#           #'Zee'    :{'type':'bkg', 'leg':'Z+jets (ee)',                    'f_color':kOrange+1,  'path':'Zee'+bkg_suffix},
#           #'Zmumu'    :{'type':'bkg', 'leg':'Z+jets (mm)',                    'f_color':kOrange+1,  'path':'Zmumu'+bkg_suffix},
#           #'Ztautau'    :{'type':'bkg', 'leg':'Z+jets (tt)',                    'f_color':kOrange+1,  'path':'Ztautau'+bkg_suffix},
#           'Wjets'    :{'type':'bkg', 'leg':'W+jets',                    'f_color':kTeal+5,  'path':'Wjets'+bkg_suffix},
#           #'PHZZ'     :{'type':'bkg', 'leg':'ZZ (powheg)',        'f_color':kAzure+7, 'path':'PHZZ'+bkg_suffix},
#           #'PHWZ'     :{'type':'bkg', 'leg':'WZ (powheg)',        'f_color':kAzure-7, 'path':'PHWZ'+bkg_suffix},
#           #'PHWW'     :{'type':'bkg', 'leg':'WW (powheg)',        'f_color':kAzure, 'path':'PHWW'+bkg_suffix},
#           'diboson'  :{'type':'bkg', 'leg':'Diboson',                   'f_color':kAzure+7, 'path':'diboson'+bkg_suffix},
#           'triboson' :{'type':'bkg', 'leg':'Triboson',                  'f_color':kGreen-9,'path':'triboson'+bkg_suffix},
#           'higgs'    :{'type':'bkg', 'leg':'Higgs',                     'f_color':kAzure+6,'path':'higgs'+bkg_suffix},
#           'lowMassDY':{'type':'bkg', 'leg':'Low mass DY',               'f_color':kMagenta-7,'path':'lowMassDY'+bkg_suffix},
#           #'DY':{'type':'bkg', 'leg':'Low mass DY',               'f_color':kMagenta-7,'path':'DY'+bkg_suffix},
#           'topOther' :{'type':'bkg', 'leg':'Top other',                 'f_color':kGray,'path':'topOther'+bkg_suffix},            
#           'data15-16'  :{'type':'data','leg':'Data (2015,2016)',        'f_color':kBlack,'l_color':kBlack,  'path': 'data'+data_suffix,'m_type':24},
#           'data17'     :{'type':'data','leg':'Data (2017)',             'f_color':kBlack,'l_color':kBlack,  'path': 'data'+data_suffix,'m_type':25},
#           'data18'     :{'type':'data','leg':'Data (2018)',             'f_color':kBlack,'l_color':kBlack,  'path': 'data'+data_suffix,'m_type':26},
#           'mc16a'            :{'type':'mc','leg':'MC',     'f_color':kBlack,'l_color':kBlack,  'path': 'data'+bkg_suffix,'m_type':23},
#           'mc16d'            :{'type':'mc','leg':'MC',     'f_color':kBlack,'l_color':kBlack,  'path': 'data'+bkg_suffix,'m_type':23},
#           'mc16e'            :{'type':'mc','leg':'MC',     'f_color':kBlack,'l_color':kBlack,  'path': 'data'+bkg_suffix,'m_type':23}}


triglist = {}
triglist["e60_lhmedium_nod0"] = {"color":kOrange+10, "year":[2016,2017,2018]}        
triglist["e60_lhmedium"] = {"color":kOrange+10, "year":[2015]}            
triglist["e24_lhmedium_L1EM20VH"] = {"color":kSpring,       "year":[2015]}  
#triglist["e24_lhmedium_nod0_L1EM20VH"] = {"color":kGray,   "year":[2016]}  
triglist["mu20_iloose_L1MU15"] = {"color":kYellow,          "year":[2015]}    
triglist["e120_lhloose"] = {"color":kBlue,                  "year":[2015]}              
triglist["e24_lhtight_nod0_ivarloose"] = {"color":kCyan,    "year":[2016]}   
triglist["e26_lhtight_nod0_ivarloose"] = {"color":kPink+10, "year":[2016,2017,2018]}   
triglist["mu26_ivarmedium"] = {"color":kGray,        "year":[2016,2017,2018]}         
triglist["e140_lhloose_nod0"] = {"color":kGreen-10,  "year":[2016,2017,2018]}             
triglist["mu50"] = {"color":kViolet+4,               "year":[2015,2016,2017,2018]}         

## ZMET
triglist = {}
triglist["2e12_lhloose_L12EM10VH"] = {"color":kOrange+10, "year":[2015]}
triglist["2e17_lhvloose_nod0"] = {"color":kOrange+10, "year":[2016]}
triglist["2e24_lhvloose_nod0"] = {"color":kOrange+10, "year":[2017,2018]}
triglist["2e17_lhvloose_nod0_L12EM15VHI"] = {"color":kOrange+10, "year":[2017,2018]}
triglist["mu26_imedium"] = {"color":kOrange+10, "year":[2015]}
triglist["mu26_ivarmedium"] = {"color":kOrange+10, "year":[2016,2017,2018]}
triglist["mu50"] = {"color":kOrange+10, "year":[2015,2016,2017,2018]}
#triglist["e17_lhloose_mu14"] = {"color":kOrange-10, "year":[2015,2016,2017,2018],"marker":37}
#triglist["e26_lhmedium_nod0_mu8noL1"] = {"color":kOrange+3, "year":[2015,2016,2017,2018],"marker":39}

## If want to plot trigger-by-trigger: uncomment this + set usetriggerbytrigger = true
usetriggerbytrigger = True
'''
trigger_setup = {}
trigger_setup["e60_lhmedium_nod0"] = {"col":kOrange+10, "year":[2016,2017,2018],"marker":22}        
trigger_setup["e60_lhmedium"] = {"col":kOrange+10, "year":[2015],"marker":23}            
trigger_setup["e24_lhmedium_L1EM20VH"] = {"col":kSpring,       "year":[2015],"marker":32}  
#trigger_setup["e24_lhmedium_nod0_L1EM20VH"] = {"col":kGray,   "year":[2016],"marker":22}  
trigger_setup["mu20_iloose_L1MU15"] = {"col":kYellow,          "year":[2015],"marker":25}    
trigger_setup["e120_lhloose"] = {"col":kBlue,                  "year":[2015],"marker":26}              
trigger_setup["e24_lhtight_nod0_ivarloose"] = {"col":kCyan,    "year":[2016],"marker":27}   
trigger_setup["e26_lhtight_nod0_ivarloose"] = {"col":kPink+10, "year":[2016,2017,2018],"marker":28}   
trigger_setup["mu26_ivarmedium"] = {"col":kGray,        "year":[2016,2017,2018],"marker":29}         
trigger_setup["e140_lhloose_nod0"] = {"col":kGreen-10,  "year":[2016,2017,2018],"marker":30}             
trigger_setup["mu50"] = {"col":kViolet+4,               "year":[2015,2016,2017,2018],"marker":31}         
'''

trigger_setup = {}
trigger_setup["2e12_lhloose_L12EM10VH"] = {"col":kOrange+10, "year":[2015],"marker":22} 
trigger_setup["2e17_lhvloose_nod0"] = {"col":kOrange+10, "year":[2016],"marker":32}  
trigger_setup["2e24_lhvloose_nod0"] = {"col":kOrange+10, "year":[2017,2018],"marker":22}  
trigger_setup["2e17_lhvloose_nod0_L12EM15VHI"] = {"col":kOrange+10, "year":[2017,2018],"marker":25}   
trigger_setup["mu26_imedium"] = {"col":kOrange+10, "year":[2015],"marker":26}  
trigger_setup["mu26_ivarmedium"] = {"col":kOrange+10, "year":[2016,2017,2018],"marker":27}   
trigger_setup["mu50"] = {"col":kOrange+10, "year":[2015,2016,2017,2018],"marker":28}
trigger_setup["e17_lhloose_mu14"] = {"col":kOrange-10, "year":[2015,2016,2017,2018],"marker":37}
trigger_setup["e26_lhmedium_nod0_mu8noL1"] = {"col":kOrange+3, "year":[2015,2016,2017,2018],"marker":39}



#



indir = sys.argv[1]
typ   = sys.argv[2]
lep   = sys.argv[3]
uncert = sys.argv[4]
analy = sys.argv[5]
if len(sys.argv) > 6:
  BDTscore = sys.argv[6]
else:
  BDTscore = ""
is18 = False
is17 = False
is1516 = False
datakey = ""
tryear = []
if "data18"     in indir or "year18" in indir: 
    datakey = "data18"
    is18   = True
    mckey = "mc16e"
    fnkey = "2018"
    tryear = [2018]
elif "data17"   in indir or "year17" in indir: 
    datakey = "data17"
    is17   = True
    mckey = "mc16d"
    fnkey = "2017"
    tryear = [2017]
elif "data1516" in indir  or "year1516" in indir: 
    datakey = "data15-16"
    is1516 = True
    mckey = "mc16a"
    fnkey = "2015"
    tryear = [2015,2016]

if analy == "FCL":
  regions_to_use = {"light":"FAKE2L05",
                    "heavy":"FAKE2L21",
                    "conv" :"FAKE2L23",
                    "real" :"REAL2L01",
                    "real1d" :"REAL2L01",
                    "real1deta" :"REAL2L01"}
elif analy == "ZMET":
  regions_to_use = {"light":"FAKE2L05",
                    "heavy":"FAKE2L21",
                    "conv" :"FAKE2L23",
                    "real" :"REAL2L04",
                    "real1d" :"REAL2L04",
                    "real1deta" :"REAL2L04"}
elif analy == "2L2J":
  triglist = {}#"":{"color":kViolet+4,"year":[2015,2016,2017,2018]}}
  trigger_setup = {"":{"triggers":[""],"col":kOrange-4,"marker":22}}
  regions_to_use = {"light":"FAKE2L05",
                    "heavy":"FAKE2L21",
                    "conv" :"FAKE2L23",
                    "real" :"REAL2L04",
                    "real1d" :"REAL2L04",
                    "real1deta" :"REAL2L04"}
  

infiles = [f for f in listdir(indir) if isfile(join(indir, f))]

notrig = []
subt_notrig = []
trig = {}
subt_trig = {}
rebx = array.array('f',[0,20,25,30,40,50,60,80,100,1000])
reby = array.array('f',[0,20,25,30,40,50,60,80,100,1000])
vari ="pT"
for key in trigger_setup.keys():
    if lep == "el" and key.startswith("m_"): continue
    if lep == "mu" and key.startswith("e_"): continue
    trig[key] = ""
    subt_trig[key] = ""

if "light" in typ:
    rebx = array.array('f',[0,20,25,30,40,60,80,100,1000])
    doPlot = [mckey,"Diboson","ttbar","Zjets","Wjets"]
    ytit = "Light Flavoured Fake Rate (electrons)"
    xtit = "p_{T}^{el} [GeV]"
    for key in triglist.keys():
        isOK = False
        for yr in tryear:
            if yr in triglist[key]["year"]:
                isOK = True
                break
        if not isOK: continue
        notrig.append("h_lep_pT_nT_%slepnotrig_ALL_E_trueLightFlavorDecay_FAKE2L05%s_%s"%("%s_"%key if key else "","_%s"%uncert if uncert else "",analy))
    for key in trigger_setup.keys():
        isOK = False
        for yr in tryear:
            if yr in trigger_setup[key]["year"]:
                isOK = True
                break
        if not isOK: continue
        if key.startswith("m_"): continue
        trig[key] = "h_lep_pT_nT_%sALL_E_trueLightFlavorDecay_FAKE2L05%s_%s"%("%s_"%key if key else "","_%s"%uncert if uncert else "",analy)
elif "conv" in typ:
    doPlot = ["Zjets",datakey]
    ytit = "Conversion Fake Rate (electrons)"
    xtit = "p_{T}^{el} [GeV]"
    for key in triglist.keys():
        isOK = False
        for yr in tryear:
            if yr in triglist[key]["year"]:
                isOK = True
                break
        if not isOK: continue
        notrig.append("h_lep_pT_nT_%slepnotrig_ALL_E_FAKE2L23_%s"%("%s_"%key if key else "",analy))
    for key in trigger_setup.keys():
        isOK = False
        for yr in tryear:
            if yr in trigger_setup[key]["year"]:
                isOK = True
                break
        if not isOK: continue
        trig[key] = "h_lep_pT_nT_%sALL_E_FAKE2L23_%s"%("%s_"%key if key else "",analy)
elif "real1deta" in typ:
    vari = "eta"
    doPlot = [mckey]
    rebx = array.array('f',[0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.7])
    for key in triglist.keys():
        isOK = False
        for yr in tryear:
            if yr in triglist[key]["year"]:
                isOK = True
                break
        if not isOK: continue
        #print("Trigger %s is ok for year %s"%("%s_"%key if key else "",fnkey))
        if lep == "mu":
            notrig.append("h_lep_eta_nT_%slepnotmatched_MMOS_M_%s%s_%s"%("%s_"%key if key else "",regions_to_use[typ],"_%s"%uncert if uncert else "",analy))
        elif lep == "el":
            notrig.append("h_lep_eta_nT_%slepnotmatched_EEOS_E_%s%s_%s"%("%s_"%key if key else "",regions_to_use[typ],"_%s"%uncert if uncert else "",analy))
    for key in trigger_setup.keys():
        isOK = False
        for yr in tryear:
            if yr in trigger_setup[key]["year"]:
                isOK = True
                break
        if not isOK: continue
        if lep == "mu":
            #rebx = array.array('f',[0,15,20,30,40,60,100,1000])
            rebx = array.array('f',[0.0,0.4,0.8,1.2,1.6,2.0,2.4,2.7])
            ytit = "Real Efficiency (muons)"
            xtit = "|#eta^{#mu}|"
            #notrig.append("h_lep_eta_eta_nT_%slepnotrig_MMOS_M_REAL2L04%s_%s"%("%s_"%key if key else "","_%s"%uncert if uncert else ""))
            if (not usetriggerbytrigger and key.startswith("e_")) or (usetriggerbytrigger and key.startswith("e")): continue
            if not usetriggerbytrigger:
              trig[key] = "h_lep_eta_nT_%slepmatched_MMOS_M_%s%s_%s"%("%s_"%key if key else "",regions_to_use[typ],"_%s"%uncert if uncert else "",analy)
            else:
              trig[key] = "h_lep_eta_nT_%slepmatched_MMOS_M_%s%s_%s"%("%s_"%key if key else "",regions_to_use[typ],"_%s"%uncert if uncert else "",analy)
        elif lep == "el":
            if vb: print("key = ",key)
            rebx = array.array('f',[0.0,0.4,0.8,1.2,1.37,1.52,1.6,2.0,2.4,2.7])
            ytit = "Real Efficiency (electrons)"
            xtit = "|#eta^{el}|"
            #notrig.append("h_lep_eta_eta_nT_%slepnotrig_EEOS_E_REAL2L04%s_%s"%("%s_"%key if key else "","_%s"%uncert if uncert else ""))
            if (not usetriggerbytrigger and key.startswith("m_")) or (usetriggerbytrigger and key.startswith("m")): continue
            if not usetriggerbytrigger:
              trig[key] = "h_lep_eta_nT_%slepmatched_EEOS_E_%s%s_%s"%("%s_"%key if key else "",regions_to_use[typ],"_%s"%uncert if uncert else "",analy)
            else:
              trig[key] = "h_lep_eta_nT_%slepmatched_EEOS_E_%s%s_%s"%("%s_"%key if key else "",regions_to_use[typ],"_%s"%uncert if uncert else "",analy)
elif "real1d" in typ:
    vari = "pT"
    doPlot = [mckey]
    if vari == "pT":
      rebx = array.array('f',[0,15,20,30,35,40,60,80,1000])
      #BDTscore = ""
    elif vari == "BDT":
      doPlot = [datakey]
      rebx = array.array('f',[])#0.05,0.1,0.15,0.20,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0])
      #BDTscore = "BDTothersDeltaM30"
      #rebx = array.array('f',[0,15,20,30,35,40,60,80,1000])
    for key in triglist.keys():
        isOK = False
        for yr in tryear:
            if yr in triglist[key]["year"]:
                isOK = True
                break
        if not isOK: continue
        #print("Trigger %s is ok for year %s"%("%s_"%key if key else "",fnkey))
        if lep == "mu":
            notrig.append("h_lep_%s_nT%s_%slepnotmatched_MMOS_M_%s%s_%s"%(vari,"_"+BDTscore if BDTscore else "","%s_"%key if key else "",regions_to_use[typ],"_%s"%uncert if uncert else "",analy))
        elif lep == "el":
            notrig.append("h_lep_%s_nT%s_%slepnotmatched_EEOS_E_%s%s_%s"%(vari,"_"+BDTscore if BDTscore else "","%s_"%key if key else "",regions_to_use[typ],"_%s"%uncert if uncert else "",analy))
    for key in trigger_setup.keys():
        isOK = False
        for yr in tryear:
            if yr in trigger_setup[key]["year"]:
                isOK = True
                break
        if not isOK: continue
        if lep == "mu":
            if vari == "pT": rebx = array.array('f',[0,15,20,30,40,60,1000])
            ytit = "Real Efficiency (muons)"
            xtit = "p_{T}^{#mu} [GeV]"
            #notrig.append("h_lep_%s_eta_nT_%slepnotrig_MMOS_M_%s%s_%s"%("%s_"%key if key else "","_%s"%uncert if uncert else ""))
            if (not usetriggerbytrigger and key.startswith("e_")) or (usetriggerbytrigger and key.startswith("e")): continue
            if not usetriggerbytrigger:
              trig[key] = "h_lep_%s_nT%s_%slepmatched_MMOS_M_%s%s_%s"%(vari,"_"+BDTscore if BDTscore else "","%s_"%key if key else "",regions_to_use[typ],"_%s"%uncert if uncert else "",analy)
            else:
              trig[key] = "h_lep_%s_nT%s_%slepmatched_MMOS_M_%s%s_%s"%(vari,"_"+BDTscore if BDTscore else "","%s_"%key if key else "",regions_to_use[typ],"_%s"%uncert if uncert else "",analy)
        elif lep == "el":
            ytit = "Real Efficiency (electrons)"
            xtit = "p_{T}^{el} [GeV]"
            #notrig.append("h_lep_%s_eta_nT_%slepnotrig_EEOS_E_%s%s_%s"%("%s_"%key if key else "","_%s"%uncert if uncert else ""))
            if (not usetriggerbytrigger and key.startswith("m_")) or (usetriggerbytrigger and key.startswith("m")): continue
            if not usetriggerbytrigger:
              trig[key] = "h_lep_%s_nT%s_%slepmatched_EEOS_E_%s%s_%s"%(vari,"_"+BDTscore if BDTscore else "","%s_"%key if key else "",regions_to_use[typ],"_%s"%uncert if uncert else "",analy)
            else:
              trig[key] = "h_lep_%s_nT%s_%slepmatched_EEOS_E_%s%s_%s"%(vari,"_"+BDTscore if BDTscore else "","%s_"%key if key else "",regions_to_use[typ],"_%s"%uncert if uncert else "",analy)
elif "real" in typ:
    doPlot = [datakey]
    rebx = array.array('f',[0,15,20,30,35,40,60,80,1000])
    reby = array.array('f',[0.0,0.4,1.0,1.37,1.52,2.0,2.7])
    for key in triglist.keys():
        isOK = False
        for yr in tryear:
            if yr in triglist[key]["year"]:
                isOK = True
                break
        if not isOK: continue
        #print("Trigger %s is ok for year %s"%("%s_"%key if key else "",fnkey))
        if lep == "mu":
            notrig.append("h_lep_pT_eta_nT_%slepnotmatched_MMOS_M_%s%s_%s"%("%s_"%key if key else "",regions_to_use[typ],"_%s"%uncert if uncert else "",analy))
        elif lep == "el":
            notrig.append("h_lep_pT_eta_nT_%slepnotmatched_EEOS_E_%s%s_%s"%("%s_"%key if key else "",regions_to_use[typ],"_%s"%uncert if uncert else "",analy))
    for key in trigger_setup.keys():
        isOK = False
        for yr in tryear:
            if yr in trigger_setup[key]["year"]:
                isOK = True
                break
        if not isOK: continue
        if lep == "mu":
            rebx = array.array('f',[0,15,20,30,40,60,100,1000])
            reby = array.array('f',[0.0,0.6,1.2,1.8,2.7])
            ytit = "|#eta^{#mu}|"
            xtit = "p_{T}^{#mu} [GeV]"
            #notrig.append("h_lep_pT_eta_nT_%slepnotrig_MMOS_M_%s%s_%s"%("%s_"%key if key else "",regions_to_use[typ],"_%s"%uncert if uncert else ""))
            if key.startswith("e_"): continue
            trig[key] = "h_lep_pT_eta_nT_%slepmatched_MMOS_M_%s%s_%s"%("%s_"%key if key else "",regions_to_use[typ],"_%s"%uncert if uncert else "",analy)
        elif lep == "el":
            ytit = "|#eta^{el}|"
            xtit = "p_{T}^{el} [GeV]"
            #notrig.append("h_lep_pT_eta_nT_%slepnotrig_EEOS_E_%s%s_%s"%("%s_"%key if key else "","_%s"%uncert if uncert else ""))
            if key.startswith("m_"): continue
            trig[key] = "h_lep_pT_eta_nT_%slepmatched_EEOS_E_%s%s_%s"%("%s_"%key if key else "",regions_to_use[typ],"_%s"%uncert if uncert else "",analy)
elif "heavy" in typ:
    doPlot = [datakey,"ttbar","Wjets"]
    for key in triglist.keys():
        isOK = False
        for yr in tryear:
            if yr in triglist[key]["year"]:
                isOK = True
                break
        if not isOK: continue
        if lep == "mu":
            notrig.append("h_lep_pT_nT_%slepnotmatched_ALL_M_%s_%s"%("%s_"%key if key else "",regions_to_use[typ],analy))
            subt_notrig.append("h_lep_pT_nT_%slepnotmatched_ALL_M_trueREAL_%s_%s"%("%s_"%key if key else "",regions_to_use[typ],analy))
        elif lep == "el":
            notrig.append("h_lep_pT_nT_%slepnotmatched_ALL_E_%s_%s"%("%s_"%key if key else "",regions_to_use[typ],analy))
            subt_notrig.append("h_lep_pT_nT_%slepnotmatched_ALL_E_trueREAL_%s_%s"%("%s_"%key if key else "",regions_to_use[typ],analy))
    for key in trigger_setup.keys():
        isOK = False
        for yr in tryear:
            if yr in trigger_setup[key]["year"]:
                isOK = True
                break
        if not isOK: continue
        if lep == "mu":
            rebx = array.array('f',[0,10,15,20,25,30,35,40,60,80,100,1000])
            ytit = "Heavy Flavoured Fake Rate (muons)"
            xtit = "p_{T}^{#mu} [GeV]"
            #rebx = array.array('f',[0,20,25,30,40,60,80,1000])
            #notrig.append("h_lep_pT_nT_%slepnotmatched_ALL_M_FAKE2L21_%s"%("%s_"%key if key else "",analy))
            #subt_notrig.append("h_lep_pT_nT_%slepnotmatched_ALL_M_trueREAL_FAKE2L21_%s"%("%s_"%key if key else "",analy))
            if (not usetriggerbytrigger and key.startswith("e_")) or (usetriggerbytrigger and key.startswith("e")): continue

            if not usetriggerbytrigger:
              trig[key] = "h_lep_pT_nT_%slepmatched_ALL_M_%s_%s"%("%s_"%key if key else "",regions_to_use[typ],analy)
              subt_trig[key] = "h_lep_pT_nT_%slepmatched_ALL_M_trueREAL_%s_%s"%("%s_"%key if key else "",regions_to_use[typ],analy)
            else:
              trig[key] = "h_lep_pT_nT_%slepmatched_ALL_M_%s_%s"%("%s_"%key if key else "",regions_to_use[typ],analy)
              subt_trig[key] = "h_lep_pT_nT_%slepmatched_ALL_M_trueREAL_%s_%s"%("%s_"%key if key else "",regions_to_use[typ],analy)

        elif lep == "el":
            if vb: print("key = ",key)
            rebx = array.array('f',[0,15,20,25,30,40,50,60,80,100,1000])
            ytit = "Heavy Flavoured Fake Rate (electrons)"
            xtit = "p_{T}^{el} [GeV]"
            #notrig.append("h_lep_pT_nT_%slepnotmatched_ALL_E_FAKE2L21_%s"%("%s_"%key if key else "",analy))
            #subt_notrig.append("h_lep_pT_nT_%slepnotmatched_ALL_E_trueREAL_FAKE2L21_%s"%("%s_"%key if key else "",analy))
            #if (not usetriggerbytrigger and key.startswith("m_")) or (usetriggerbytrigger and key.startswith("m")): continue
          
            if not usetriggerbytrigger:
              trig[key] = "h_lep_pT_nT_%slepmatched_ALL_E_%s_%s"%("%s_"%key if key else "",regions_to_use[typ],analy)
              subt_trig[key] = "h_lep_pT_nT_%slepmatched_ALL_E_trueREAL_%s_%s"%("%s_"%key if key else "",regions_to_use[typ],analy)
            else:
              trig[key] = "h_lep_pT_nT_%slepmatched_ALL_E_%s_%s"%("%s_"%key if key else "",regions_to_use[typ],analy)
              subt_trig[key] = "h_lep_pT_nT_%slepmatched_ALL_E_trueREAL_%s_%s"%("%s_"%key if key else "",regions_to_use[typ],analy)

ht_notrig = {}
hl_notrig = {}

ht_trig = {}
hl_trig = {}

ht_subt_notrig = {}
hl_subt_notrig = {}

ht_subt_trig = {}
hl_subt_trig = {}
didplot = []

h_trig_err_xsec = {}
h_notrig_err_xsec = {}
for f in infiles:
    proc = f.split("_")[0]
    if not "_HISTOGRAMS_" in f: continue
    # if proc in ht_notrig.keys():
    #     print("Process %s has already been read. Check!"%proc)
    #     continue
    if not f.split("_")[1] == "merged" and not f.split("_")[3] == "merged":
        print("File %s has bad formatting. Check!"%f)
        continue
    if not proc in d_samp.keys():
        print("Not plotting %s"%proc)
        continue
    didplot.append(proc)
    if typ == "light" and "data" in proc: continue
    if uncert and "data" in proc: 
      print("Not plotting %s"%proc)
      continue
    if not proc in ht_notrig.keys():
      ht_notrig[proc] = 0
      hl_notrig[proc] = 0
      if typ == "heavy":
        ht_subt_notrig[proc] = 0
        hl_subt_notrig[proc] = 0
    
    mmfile = TFile(indir+"/"+f)

    # -------
    # Histograms for non-triggered
    # -------
    # Get the non-triggered rates
    for nt in notrig:

        if not proc == datakey:
          print("Switching name from %s"%nt)
          nt = nt.replace(regions_to_use[typ],"trueFAKE_%s"%regions_to_use[typ])
          print("to %s for %s"%(nt,proc))
        # Get first the tight histogram
        objPtr = mmfile.Get(nt);
        #for i in range(1,objPtr.GetNbinsX()+1):
        #  print("Bin %i : [%.2f,%.2f]"%(i,objPtr.GetXaxis().GetBinLowEdge(i),objPtr.GetXaxis().GetBinUpEdge(i)))
        if objPtr and type(ht_notrig[proc]) is int: 
            ht_notrig[proc] = objPtr.Clone(objPtr.GetName()+"_%s"%proc)
            ht_notrig[proc].SetDirectory(0)
            if vb: print("Adding for first %s from %s, now %.2f"%(nt,proc,ht_notrig[proc].Integral()))
        elif objPtr:
            ht_notrig[proc].Add(objPtr)
            if vb: print("Adding on top %s from %s, now %.2f"%(nt,proc,ht_notrig[proc].Integral()))
        else:
            print("1) Could not find %s in %s"%(nt,indir+"/"+f))

        # Then get the loose histogram
        objPtr = mmfile.Get(nt.replace("nT","nL"));
        if objPtr and type(hl_notrig[proc]) is int: 
            hl_notrig[proc] = objPtr.Clone(objPtr.GetName()+"_%s"%proc)
            hl_notrig[proc].SetDirectory(0)
            #print("Adding %s from %s, now %.2f"%(nt,proc,hl_notrig[proc].Integral()))
        elif objPtr:
            hl_notrig[proc].Add(objPtr)
            #print("Adding %s from %s, now %.2f"%(nt,proc,hl_notrig[proc].Integral()))
        else:
            print("2) Could not find %s in %s"%(nt.replace("nT","nL"),indir+"/"+f))

    # -------
    # Histograms for subtraction
    # -------
    # Get the non-triggered rates (for subtraction)
    for nt in subt_notrig:
        if "data" in proc: continue
        # Get first the tight histogram
        objPtr = mmfile.Get(nt);
        if objPtr and type(ht_subt_notrig[proc]) is int: 
            ht_subt_notrig[proc] = objPtr.Clone(objPtr.GetName()+"_%s"%proc)
            ht_subt_notrig[proc].SetDirectory(0)
        elif objPtr:
            ht_subt_notrig[proc].Add(objPtr)
            #print("Adding %s from %s, now %.2f"%(nt,proc,ht_subt_notrig[proc].Integral()))
        else:
            print("3) Could not find %s in %s"%(nt,indir+"/"+f))

        # Then get the loose histogram
        objPtr = mmfile.Get(nt.replace("nT","nL"));
        if objPtr and type(hl_subt_notrig[proc]) is int: 
            hl_subt_notrig[proc] = objPtr.Clone(objPtr.GetName()+"_%s"%proc)
            hl_subt_notrig[proc].SetDirectory(0)
        elif objPtr:
            hl_subt_notrig[proc].Add(objPtr)
            #print("Adding %s from %s, now %.2f"%(nt.replace("nT","nL"),proc,hl_subt_notrig[proc].Integral()))
        else:
            print("4) Could not find %s in %s"%(nt.replace("nT","nL"),indir+"/"+f))

    # -------
    # Histograms for triggers
    # -------
    # Get the non-triggered rates
    for key in trig.keys():
        if not trig[key]: continue

        if not proc == datakey and not "trueFAKE" in trig[key]:
          print("Switching name from %s"%nt)
          trig[key] = trig[key].replace(regions_to_use[typ],"trueFAKE_%s"%regions_to_use[typ])
          print("to %s for "%nt)
        if proc == datakey and "trueFAKE" in trig[key]:
          trig[key] = trig[key].replace("trueFAKE_%s"%regions_to_use[typ],regions_to_use[typ])

        extrakey = ""
        if key in trigger_categories.keys():
          extrakey = trigger_categories[key]
          if not extrakey in ht_trig.keys():
            ht_trig[extrakey] = {proc:0}
          elif not proc in ht_trig.keys():
            ht_trig[extrakey][proc] = 0

        if not extrakey:
          print("Did not find category for trigger %s"%key)

        if not key in ht_trig.keys():
            ht_trig[key] = {proc:0}
        elif not proc in ht_trig[key].keys():
            ht_trig[key][proc] = 0

            
        # Get first the tight histogram
        objPtr = mmfile.Get(trig[key]);
        if objPtr and type(ht_trig[key][proc]) is int:
            ht_trig[key][proc] = objPtr.Clone(objPtr.GetName()+"_%s"%proc)
            ht_trig[key][proc].SetDirectory(0)
            if vb: print("Trig :: Adding for %s first %s from %s, now %.2f"%(key,trig[key],proc,ht_trig[key][proc].Integral()))
        elif objPtr:
            ht_trig[key][proc].Add(objPtr)
            if vb: print("Trig :: Adding %s on top %s from %s, now %.2f"%(key,trig[key],proc,ht_trig[key][proc].Integral()))
        else:
            print("5) Could not find %s in %s"%(trig[key],indir+"/"+f))

        if extrakey and objPtr and type(ht_trig[extrakey][proc]) is int:
          ht_trig[extrakey][proc] = objPtr.Clone(objPtr.GetName().replace(trig[key],extrakey)+"_%s"%proc)
          ht_trig[extrakey][proc].SetDirectory(0)
          if vb: print("Trig :: Adding for %s first %s from %s, now %.2f"%(extrakey,trig[key],proc,ht_trig[extrakey][proc].Integral()))
        elif extrakey and objPtr:
          ht_trig[extrakey][proc].Add(objPtr)
          if vb: print("Trig :: Adding %s on top %s from %s, now %.2f"%(extrakey,trig[key],proc,ht_trig[extrakey][proc].Integral()))
        elif extrakey:
          print("6) Could not find %s in %s"%(trig[key],indir+"/"+f))


        extrakey = ""
        if key in trigger_categories.keys():
          extrakey = trigger_categories[key]
          if not extrakey in hl_trig.keys():
            hl_trig[extrakey] = {proc:0}
          elif not proc in hl_trig.keys():
            hl_trig[extrakey][proc] = 0

        if not key in hl_trig.keys():
            hl_trig[key] = {proc:0}
        elif not proc in hl_trig[key].keys():
            hl_trig[key][proc] = 0
            
        # Then get the loose histogram
        objPtr = mmfile.Get(trig[key].replace("nT","nL"));
        if objPtr and type(hl_trig[key][proc]) is int:
            hl_trig[key][proc] = objPtr.Clone(objPtr.GetName()+"_%s"%proc)
            hl_trig[key][proc].SetDirectory(0)
        elif objPtr:
            hl_trig[key][proc].Add(objPtr)
        else:
            print("7) Could not find %s in %s"%(trig[key].replace("nT","nL"),indir+"/"+f))

        if extrakey and objPtr and type(hl_trig[extrakey][proc]) is int:
          hl_trig[extrakey][proc] = objPtr.Clone(objPtr.GetName().replace(trig[key],extrakey)+"_%s"%proc)
          hl_trig[extrakey][proc].SetDirectory(0)
        elif extrakey and objPtr:
          hl_trig[extrakey][proc].Add(objPtr)
        elif extrakey:
          print("8) Could not find %s in %s"%(trig[extrakey].replace("nT","nL"),indir+"/"+f))


        if "data" in proc: continue
        # -------
        # Histograms for subtraction
        # -------
        if typ == "heavy":
          
            if not key in ht_subt_trig.keys():
                ht_subt_trig[key] = {proc:0}
            elif not proc in ht_subt_trig[key].keys():
                ht_subt_trig[key][proc] = 0
                
            # Get first the tight histogram
            objPtr = mmfile.Get(subt_trig[key]);
            if objPtr and type(ht_subt_trig[key][proc]) is int:
                ht_subt_trig[key][proc] = objPtr.Clone(objPtr.GetName()+"_%s"%proc)
                ht_subt_trig[key][proc].SetDirectory(0)
                if vb: print("Trig (subt) :: Adding for %s first %s from %s, now %.2f"%(key,subt_trig[key],proc,ht_subt_trig[key][proc].Integral()))
            elif objPtr:
                ht_subt_trig[key][proc].Add(objPtr)
                if vb: print("Trig (subt) :: Adding %s on top %s from %s, now %.2f"%(extrakey,subt_trig[key],proc,ht_subt_trig[key][proc].Integral()))
            else:
                print("9) Could not find %s in %s"%(subt_trig[key],indir+"/"+f))

            # And for the trigger categories
            if extrakey and not extrakey in ht_subt_trig.keys():
                ht_subt_trig[extrakey] = {proc:0}
            elif extrakey and not proc in ht_subt_trig[extrakey].keys():
                ht_subt_trig[extrakey][proc] = 0
            
            if extrakey and objPtr and type(ht_subt_trig[extrakey][proc]) is int:
                ht_subt_trig[extrakey][proc] = objPtr.Clone(objPtr.GetName().replace(trig[key],extrakey)+"_%s"%proc)
                ht_subt_trig[extrakey][proc].SetDirectory(0)
                if vb: print("Trig (subt) :: Adding for %s first %s from %s, now %.2f"%(extrakey,subt_trig[key],proc,ht_subt_trig[extrakey][proc].Integral()))
            elif extrakey and objPtr:
                ht_subt_trig[extrakey][proc].Add(objPtr)
                if vb: print("Trig (subt) :: Adding %s on top %s from %s, now %.2f"%(extrakey,subt_trig[key],proc,ht_subt_trig[extrakey][proc].Integral()))
            else:
                print("10) Could not find %s in %s for extrakey %s"%(subt_trig[key],indir+"/"+f,extrakey))

            if not key in hl_subt_trig.keys():
                hl_subt_trig[key] = {proc:0}
            elif not proc in hl_subt_trig[key].keys():
                hl_subt_trig[key][proc] = 0
                
            # Then get the loose histogram
            objPtr = mmfile.Get(subt_trig[key].replace("nT","nL"));
            if objPtr and type(hl_subt_trig[key][proc]) is int:
                hl_subt_trig[key][proc] = objPtr.Clone(objPtr.GetName()+"_%s"%proc)
                hl_subt_trig[key][proc].SetDirectory(0)
            elif objPtr:
                hl_subt_trig[key][proc].Add(objPtr)
            else:
                print("11) Could not find %s in %s"%(subt_trig[key],indir+"/"+f))

            # And for the trigger categories
            if not extrakey in hl_subt_trig.keys() and extrakey:
                hl_subt_trig[extrakey] = {proc:0}
            elif extrakey and not proc in hl_subt_trig[extrakey].keys():
                hl_subt_trig[extrakey][proc] = 0
            
            if extrakey and objPtr and type(hl_subt_trig[extrakey][proc]) is int:
                hl_subt_trig[extrakey][proc] = objPtr.Clone(objPtr.GetName().replace(trig[key],extrakey)+"_%s"%proc)
                hl_subt_trig[extrakey][proc].SetDirectory(0)
            elif extrakey and objPtr:
                hl_subt_trig[extrakey][proc].Add(objPtr)
            else:
                print("12) Could not find %s in %s for extrakey %s"%(subt_trig[key],indir+"/"+f,extrakey))

    mmfile.Close() 


# -------
# Non-trig histograms
# -------
# Add sum for tight histogram
ht_notrig_copy = copy.deepcopy(ht_notrig)
for proc in ht_notrig_copy.keys():
  if type(ht_notrig[proc]) is int: continue
  if not "data" in proc: 
    if not mckey in ht_notrig.keys():
      ht_notrig[mckey] = ht_notrig[proc].Clone(ht_notrig[proc].GetName().replace(proc,mckey))
    else:
      ht_notrig[mckey].Add(ht_notrig[proc])
  elif is1516:
    if not datakey in  ht_notrig.keys():
      ht_notrig[datakey] = ht_notrig[proc].Clone(ht_notrig[proc].GetName().replace(proc,datakey))
    elif is1516:
      ht_notrig[datakey].Add(ht_notrig[proc])
# Subtraction
ht_subt_notrig_copy = copy.deepcopy(ht_subt_notrig)
for proc in ht_subt_notrig_copy.keys():
  if type(ht_subt_notrig[proc]) is int: continue
  if not "data" in proc: 
    if not mckey in ht_subt_notrig.keys():
        ht_subt_notrig[mckey] = ht_subt_notrig[proc].Clone(ht_subt_notrig[proc].GetName().replace(proc,mckey))
    else:
        ht_subt_notrig[mckey].Add(ht_subt_notrig[proc])
  elif is1516:
    if not datakey in ht_subt_notrig.keys():
      ht_subt_notrig[datakey] = ht_subt_notrig[proc].Clone(ht_subt_notrig[proc].GetName().replace(proc,datakey))
    else:
      ht_subt_notrig[datakey].Add(ht_subt_notrig[proc])
# Add sum for loose histogram
hl_notrig_copy = copy.deepcopy(hl_notrig)
for proc in hl_notrig_copy.keys():
  if type(hl_notrig[proc]) is int: continue
  if not "data" in proc:
    if not mckey in hl_notrig.keys():
        hl_notrig[mckey] = hl_notrig[proc].Clone(hl_notrig[proc].GetName().replace(proc,mckey))
    else:
        hl_notrig[mckey].Add(hl_notrig[proc])
  elif is1516:
    if not datakey in hl_notrig.keys():
      if vb: print("Data Adding %s"%proc)
      hl_notrig[datakey] = hl_notrig[proc].Clone(hl_notrig[proc].GetName().replace(proc,datakey))
    else:
      if vb: print("Data Adding %s"%proc)
      hl_notrig[datakey].Add(hl_notrig[proc])
# Subtraction
hl_subt_notrig_copy = copy.deepcopy(hl_subt_notrig)
for proc in hl_subt_notrig_copy.keys():
  if type(hl_subt_notrig[proc]) is int: continue
  if not "data" in proc:  
    if not mckey in hl_subt_notrig.keys():
        hl_subt_notrig[mckey] = hl_subt_notrig[proc].Clone(hl_subt_notrig[proc].GetName().replace(proc,mckey))
    else:
        hl_subt_notrig[mckey].Add(hl_subt_notrig[proc])
  elif is1516:
    if not datakey in hl_subt_notrig.keys():
      hl_subt_notrig[datakey] = hl_subt_notrig[proc].Clone(hl_subt_notrig[proc].GetName().replace(proc,datakey))
    else:
      hl_subt_notrig[datakey].Add(hl_subt_notrig[proc])

# -------
# Trig histograms
# -------
# Add sum for tight histogram
ht_trig_copy = copy.deepcopy(ht_trig)
for key in ht_trig_copy.keys():
    for proc in ht_trig_copy[key].keys():
      if type(ht_trig[key][proc]) is int: continue
      if not "data" in proc: 
        if not mckey in ht_trig[key].keys():
            ht_trig[key][mckey] = ht_trig[key][proc].Clone(ht_trig[key][proc].GetName().replace(proc,mckey))
        else:
            ht_trig[key][mckey].Add(ht_trig[key][proc])
      elif is1516:
        if not datakey in ht_trig[key].keys():
          ht_trig[key][datakey] = ht_trig[key][proc].Clone(ht_trig[key][proc].GetName().replace(proc,datakey))
          if vb: print("ht_trig %s starting adding %s, int = %.2f"%(key,proc,ht_trig[key][datakey].Integral()))
        else:
          ht_trig[key][datakey].Add(ht_trig[key][proc])
          if vb: print("ht_trig %s adding %s int = %.2f"%(key,proc,ht_trig[key][datakey].Integral()))
# Subtraction
ht_subt_trig_copy = copy.deepcopy(ht_subt_trig)
for key in ht_subt_trig_copy.keys():
    for proc in ht_subt_trig_copy[key].keys():
      if type(ht_subt_trig[key][proc]) is int: continue
      if not "data" in proc:
        if not mckey in ht_subt_trig[key].keys():
            ht_subt_trig[key][mckey] = ht_subt_trig[key][proc].Clone(ht_subt_trig[key][proc].GetName().replace(proc,mckey))
        else:
            ht_subt_trig[key][mckey].Add(ht_subt_trig[key][proc])
      elif is1516:
        if not datakey in ht_subt_trig[key].keys():
            ht_subt_trig[key][datakey] = ht_subt_trig[key][proc].Clone(ht_subt_trig[key][proc].GetName().replace(proc,datakey))
        else:
            ht_subt_trig[key][datakey].Add(ht_subt_trig[key][proc])
# Add sum for loose histogram
hl_trig_copy = copy.deepcopy(hl_trig)
for key in hl_trig_copy.keys():
    for proc in hl_trig_copy[key].keys():
      if type(hl_trig[key][proc]) is int: continue
      if not "data" in proc: 
        if not mckey in hl_trig[key].keys():
          hl_trig[key][mckey] = hl_trig[key][proc].Clone(hl_trig[key][proc].GetName().replace(proc,mckey))
        else:
          hl_trig[key][mckey].Add(hl_trig[key][proc])
      elif is1516:
        if not datakey in hl_trig[key].keys():
          hl_trig[key][datakey] = hl_trig[key][proc].Clone(hl_trig[key][proc].GetName().replace(proc,datakey))
        else:
          hl_trig[key][datakey].Add(hl_trig[key][proc])
# Subtraction
hl_subt_trig_copy = copy.deepcopy(hl_subt_trig)
for key in hl_subt_trig_copy.keys():
    for proc in hl_subt_trig_copy[key].keys():
      if type(hl_subt_trig[key][proc]) is int: continue
      if not "data" in proc:
        if not mckey in hl_subt_trig[key].keys():
            hl_subt_trig[key][mckey] = hl_subt_trig[key][proc].Clone(hl_subt_trig[key][proc].GetName().replace(proc,mckey))
        else:
            hl_subt_trig[key][mckey].Add(hl_subt_trig[key][proc])
      elif is1516:
        if not datakey in hl_subt_trig[key].keys():
            hl_subt_trig[key][datakey] = hl_subt_trig[key][proc].Clone(hl_subt_trig[key][proc].GetName().replace(proc,datakey))
        else:
            hl_subt_trig[key][datakey].Add(hl_subt_trig[key][proc])


trigger_setup.update({"e_lhoose":{"col":kYellow-7,"marker":23},
                      "e_lhvloose":{"col":kOrange-4,"marker":22},
                      "e_lhvlooseL1":{"col":kGreen-7,"marker":21},
                      "m_imed":{"col":kGreen-7,"marker":24},
                      "m_ivarmed":{"col":kGreen-7,"marker":25},
                      "m_no":{"col":kGreen-7,"marker":26},
                      "em_lhloose":{"col":kGreen-7,"marker":27},
                      "em_lhmed":{"col":kGreen-7,"marker":27}})

## Legend used for plotting nT and nL distributions (to check that everything's looking good!)
diag_leg = TLegend(0.75,0.75,0.99,0.95)
diag_leg.SetBorderSize(0)
diag_leg.SetTextFont(42)
diag_leg.SetTextSize(0.03)
diag_leg.SetTextColor(1)
diag_leg.SetFillStyle(0)
diag_leg.SetLineColor(0)
diag_leg_filled = []

heff_trig = {}
heff_trig_xsecup = {}
heff_trig_xsecdw = {}
ht_trig_xsecup = {}
ht_trig_xsecdw = {}
hl_trig_xsecup = {}
hl_trig_xsecdw = {}
h_lT_trig = {}
h_lL_trig = {}
h_lT_trig_nosub = {}
h_lL_trig_nosub = {}
for key in hl_trig.keys():
    if not key in heff_trig.keys():
        heff_trig[key] = {}
        heff_trig_xsecup[key] = {}
        heff_trig_xsecdw[key] = {}
        ht_trig_xsecup[key] = {}
        ht_trig_xsecdw[key] = {}
        hl_trig_xsecup[key] = {}
        hl_trig_xsecdw[key] = {}
        h_lT_trig[key] = THStack()
        h_lL_trig[key] = THStack()
        h_lT_trig_nosub[key] = {}
        h_lL_trig_nosub[key] = {}
    if not key in ht_trig.keys():
        print("ERROR \t Could not find %s in ht_trig, but exists in hl_trig"%key)
        continue
    for proc in hl_trig[key].keys():
        if(key == "e_tight_iso" and typ == "real"):
          newrebx = array.array('f',[0,15,20,30,35,40,70,100,1000])
          newreby = array.array('f',[0.0,0.4,1.0,1.37,1.52,2.0,2.7])
          ht_trig[key][proc] = doRebin(ht_trig[key][proc],newrebx,newreby)
          hl_trig[key][proc] = doRebin(hl_trig[key][proc],newrebx,newreby)
        elif(key == "m_loose" and typ == "real"):
          newrebx = array.array('f',[0,15,20,30,60,100,1000])
          newreby = array.array('f',[0.0,0.6,1.2,1.8,2.7])
          ht_trig[key][proc] = doRebin(ht_trig[key][proc],newrebx,newreby)
          hl_trig[key][proc] = doRebin(hl_trig[key][proc],newrebx,newreby)
        elif len(rebx) > 0:
          ht_trig[key][proc] = doRebin(ht_trig[key][proc],rebx,reby)
          hl_trig[key][proc] = doRebin(hl_trig[key][proc],rebx,reby)
    # Need to do this after all the rebinning for subtraction to work!
    for proc in hl_trig[key].keys():
        if proc == datakey:
          h_lT_trig_nosub[key][proc] = ht_trig[key][proc]
          h_lL_trig_nosub[key][proc] = hl_trig[key][proc]
        if typ == "heavy" and proc == datakey: 
            ht_trig_xsecup[key][proc], hl_trig_xsecup[key][proc] = doSubt(ht_trig[key][proc],hl_trig[key][proc],ht_subt_trig[key][mckey],hl_subt_trig[key][mckey],1.1)
            ht_trig_xsecdw[key][proc], hl_trig_xsecdw[key][proc] = doSubt(ht_trig[key][proc],hl_trig[key][proc],ht_subt_trig[key][mckey],hl_subt_trig[key][mckey],0.9)
            ht_trig[key][proc], hl_trig[key][proc] = doSubt(ht_trig[key][proc],hl_trig[key][proc],ht_subt_trig[key][mckey],hl_subt_trig[key][mckey])
            heff_trig_xsecup[key][proc] = getEff(ht_trig_xsecup[key][proc],hl_trig_xsecup[key][proc])
            heff_trig_xsecup[key][proc].SetDirectory(0)
            heff_trig_xsecdw[key][proc] = getEff(ht_trig_xsecdw[key][proc],hl_trig_xsecdw[key][proc])
            heff_trig_xsecdw[key][proc].SetDirectory(0)
        heff_trig[key][proc] = getEff(ht_trig[key][proc],hl_trig[key][proc])
        heff_trig[key][proc].SetDirectory(0)
        if not "data" in proc and not mckey in proc:
          if vb: print("Adding %s"%proc)
          ht_trig[key][proc].SetLineColor(d_samp[proc]["f_color"])
          ht_trig[key][proc].SetFillColor(d_samp[proc]["f_color"])
          hl_trig[key][proc].SetLineColor(d_samp[proc]["f_color"])
          hl_trig[key][proc].SetFillColor(d_samp[proc]["f_color"])
          h_lT_trig[key].Add(ht_trig[key][proc])
          h_lL_trig[key].Add(hl_trig[key][proc])
          if not proc in diag_leg_filled:
            diag_leg.AddEntry(ht_trig[key][proc],"%s"%d_samp[proc]['leg'],"f")
            diag_leg_filled.append(proc)
        

heff_notrig = {}
heff_notrig_xsecup = {}
heff_notrig_xsecdw = {}
ht_notrig_xsecup = {}
ht_notrig_xsecdw = {}
hl_notrig_xsecup = {}
hl_notrig_xsecdw = {}
for proc in hl_notrig.keys():
    if type(hl_notrig[proc]) is int:
      continue
    # if not proc in heff_notrig.keys():
    #     heff_notrig[proc] = {}
    #     heff_notrig_xsecup[proc] = {}
    #     heff_notrig_xsecdw[proc] = {}
    #     ht_notrig_xsecup[proc] = {}
    #     ht_notrig_xsecdw[proc] = {}
    #     hl_notrig_xsecup[proc] = {}
    #     hl_notrig_xsecdw[proc] = {}
    if not proc in ht_notrig.keys():
        print("ERROR \t Could not find %s in ht_notrig, but exists in hl_notrig"%proc)
        continue
    if len(rebx) > 0:
      ht_notrig[proc] = doRebin(ht_notrig[proc],rebx,reby)
      hl_notrig[proc] = doRebin(hl_notrig[proc],rebx,reby)

# Need to do this after all the rebinning for subtraction to work!
for proc in hl_notrig.keys():
    if type(hl_notrig[proc]) is int: continue
    if typ == "heavy" and proc == datakey:
        ht_notrig_xsecup[proc], hl_notrig_xsecup[proc] = doSubt(ht_notrig[proc],hl_notrig[proc],ht_subt_notrig[mckey],hl_subt_notrig[mckey],1.1)
        ht_notrig_xsecdw[proc], hl_notrig_xsecdw[proc] = doSubt(ht_notrig[proc],hl_notrig[proc],ht_subt_notrig[mckey],hl_subt_notrig[mckey],0.9)
        ht_notrig[proc], hl_notrig[proc] = doSubt(ht_notrig[proc],hl_notrig[proc],ht_subt_notrig[mckey],hl_subt_notrig[mckey])
        #for i in range(1,ht_notrig[proc].GetNbinsX()+1):
        #  print("%.2f +/- %.2f"%(ht_notrig[proc].GetBinContent(i),ht_notrig[proc].GetBinError(i)))
        heff_notrig_xsecup[proc] = getEff(ht_notrig_xsecup[proc],hl_notrig_xsecup[proc])
        heff_notrig_xsecdw[proc] = getEff(ht_notrig_xsecdw[proc],hl_notrig_xsecdw[proc])
    heff_notrig[proc] = getEff(ht_notrig[proc],hl_notrig[proc])
    if type(heff_notrig[proc]) is int:
        print("WARNING \t Histogram heff_notrig for %s is empty"%proc)
    else:
        heff_notrig[proc].SetDirectory(0)


## Assign file names for saving to MM ntuple
prefix = "eff_"
if typ == "real":
  prefix += "pt_eta"
elif typ=="real1d":
  prefix += "%s%s"%(vari,"" if not BDTscore else "_"+BDTscore)
elif typ == "real1deta":
  prefix += "eta"


prefix += "_"+analy+"_"+fnkey+"_"+regions_to_use[typ]
if uncert:
    prefix += "_"+uncert

##
## Plotting
##
canv = []
c = TCanvas("c_%s_%s_%s%s"%(datakey,lep,typ,"_"+uncert if uncert else "_nominal"),"c_%s_%s_%s_%s"%(datakey,lep,typ,"_"+uncert if uncert else "%s_nominal"%lep),1000,600)
c.cd()
if not "eta" in typ and not vari == "BDT": c.SetLogx()
c.SetGridx()
c.SetGridy()
if typ == "real1deta":
  leg = TLegend(0.20,0.47,0.43,0.70)
elif typ == "light":
  leg = TLegend(0.60,0.80,0.80,0.92)
  leg2 = TLegend(0.75,0.80,0.95,0.92)
  leg2.SetBorderSize(0)
  leg2.SetTextFont(42)
  leg2.SetTextSize(0.03)
  leg2.SetTextColor(1)
  leg2.SetFillStyle(0)
  leg2.SetLineColor(0)
elif not typ == "real1d":
  leg = TLegend(0.18,0.72,0.41,0.94)
elif typ == "real1d":
  leg = TLegend(0.55,0.50,0.79,0.70)
leg.SetBorderSize(0)
leg.SetTextFont(42)
leg.SetTextSize(0.03)
leg.SetTextColor(1)
leg.SetFillStyle(0)
leg.SetLineColor(0)
firstDraw = True
utfil = TFile("FNPntuple_%s_JULY2022_R21.root"%(analy),"UPDATE")
added = False
trigcatinleg = []
mcinleg = []
for key in heff_trig.keys():
    if plot_categories and key in trigger_categories.keys(): continue
    for proc in heff_trig[key].keys():
        if not proc in doPlot: continue
        if not heff_trig[key][proc].Integral():
          #print("Skipping %s because empty!"%key)
          continue
        if typ == "real":
            canv.append(c.DrawClone())
            canv[-1].SetName(c.GetName()+"_%s_%s"%(key,proc))
            canv[-1].SetTitle(c.GetName()+"_%s_%s"%(key,proc))
            canv[-1].cd()
            if not "eta" in typ: canv[-1].SetLogx()
            canv[-1].SetLeftMargin(0.1)
            canv[-1].SetRightMargin(0.15)
            firstDraw = True
        heff_trig[key][proc].SetLineColor(d_samp[proc]["f_color"])
        heff_trig[key][proc].SetMarkerColor(d_samp[proc]["f_color"])
        heff_trig[key][proc].SetMarkerStyle(trigger_setup[key]["marker"])
        if typ == "light" and not key in trigcatinleg:
          leg.AddEntry(heff_trig[key][proc],"%s"%key,"lp")
          trigcatinleg.append(key)
        if typ == "light" and not d_samp[proc]["leg"] in mcinleg:
          leg2.AddEntry(heff_trig[key][proc],"%s"%d_samp[proc]["leg"],"fl")
          mcinleg.append(d_samp[proc]["leg"])
        if not typ == "light":
          leg.AddEntry(heff_trig[key][proc],"%s: "%key+d_samp[proc]["leg"],"flp")
        if firstDraw:
            heff_trig[key][proc].GetXaxis().SetMoreLogLabels()
            heff_trig[key][proc].GetXaxis().SetTitle(xtit)
            heff_trig[key][proc].GetYaxis().SetTitle(ytit)
            if typ == "real":
                heff_trig[key][proc].GetYaxis().SetTitleOffset(0.8)
                heff_trig[key][proc].Draw("colz text e")
                if analy == "FCL": myText(0.10, 0.90, kBlack, "Trig category: %s" %key)
            elif typ == "real1d" or typ == "real1deta":
                heff_trig[key][proc].GetYaxis().SetRangeUser(0.2 if typ == "real1d" else 0.3,1.0)
                heff_trig[key][proc].Draw("e1")
            else:
                if typ == "heavy":
                  heff_trig[key][proc].GetYaxis().SetRangeUser(0,1)
                if typ == "light":
                  heff_trig[key][proc].GetYaxis().SetRangeUser(0,0.8)
                heff_trig[key][proc].Draw("e1")
            firstDraw = False
        else:
            heff_trig[key][proc].Draw("same e1")

        if proc in heff_trig_xsecup[key].keys() and not uncert:
            heff_trig_xsecup[key][proc].Write(prefix.replace(analy,key)+"_xsecup_"+proc if (analy == "FCL" or analy == "ZMET") else prefix.replace(analy,lep)+"_xsecup_"+proc)
            heff_trig_xsecdw[key][proc].Write(prefix.replace(analy,key)+"_xsecdw_"+proc if (analy == "FCL" or analy == "ZMET") else prefix.replace(analy,lep)+"_xsecdw_"+proc)
        heff_trig[key][proc].Write(prefix.replace(analy,key)+"_"+proc if (analy == "FCL" or analy == "ZMET") else prefix.replace(analy,lep)+"_"+proc)
        print("Writing %s"%prefix.replace(analy,key)+"_"+proc if (analy == "FCL" or analy == "ZMET") else prefix.replace(analy,lep)+"_"+proc)
        

        if typ == "heavy" and proc == datakey:
            if not key in h_trig_err_xsec.keys():
                h_trig_err_xsec[key] = {}                
            x,y,error_up_x,error_up_y,error_low_x,error_low_y = getxsecerr(heff_trig_xsecup[key][proc],heff_trig_xsecdw[key][proc],heff_trig[key][proc])
            h_trig_err_xsec[key][proc] = TGraphAsymmErrors(len(error_low_x),x,y,error_low_x,error_up_x,error_low_y,error_up_y)
            h_trig_err_xsec[key][proc].SetLineStyle(kDashed);
            h_trig_err_xsec[key][proc].SetLineColor(kGray);
            h_trig_err_xsec[key][proc].SetFillColor(kGray);
            h_trig_err_xsec[key][proc].SetFillStyle(3018)
            h_trig_err_xsec[key][proc].SetMarkerSize(0)
            h_trig_err_xsec[key][proc].Draw("e2same");
            if not added:
                leg.AddEntry(h_trig_err_xsec[key][proc],"MC subtr. unc.","f")
                added = True

for proc in heff_notrig.keys():
    if not proc in doPlot: continue
    if not heff_notrig[proc].Integral(): continue
    if typ == "real":
        canv.append(c.DrawClone())
        canv[-1].SetName(c.GetName()+"_notrig_%s"%(proc))
        canv[-1].SetTitle(c.GetName()+"_notrig_%s"%(proc))
        canv[-1].cd()
        if not "eta" in typ and not vari == "BDT": canv[-1].SetLogx()
        canv[-1].SetLeftMargin(0.1)
        canv[-1].SetRightMargin(0.15)
        firstDraw = True
    heff_notrig[proc].SetLineColor(d_samp[proc]["f_color"])
    heff_notrig[proc].SetMarkerColor(d_samp[proc]["f_color"])
    heff_notrig[proc].SetMarkerStyle(24)
    if typ == "light" and not "notrig" in trigcatinleg:
      leg.AddEntry(heff_notrig[proc],"notrig","flp")
      trigcatinleg.append("notrig")
    if typ == "light" and not d_samp[proc]["leg"] in mcinleg:
      leg2.AddEntry(heff_notrig[proc],"%s"%d_samp[proc]["leg"],"lp")
      mcinleg.append(d_samp[proc]["leg"])
    if not typ == "light":
      leg.AddEntry(heff_notrig[proc],"m_notrig: "+d_samp[proc]["leg"],"flp")

    if firstDraw:
        heff_notrig[proc].GetXaxis().SetMoreLogLabels()
        heff_notrig[proc].GetXaxis().SetTitle(xtit)
        heff_notrig[proc].GetYaxis().SetTitle(ytit)
        if typ == "real":
            heff_notrig[proc].GetYaxis().SetTitleOffset(0.8)
            heff_notrig[proc].Draw("colz text e")
            if analy == "FCL": myText(0.10, 0.90, kBlack, "Trig category: No trig match")
        elif typ == "real1d" or typ == "real1deta":
            heff_notrig[proc].GetYaxis().SetRangeUser(0.2 if typ == "real1d" else 0.3,1.0)
            heff_notrig[proc].Draw("e1")
        else:
          if typ == "heavy":
            heff_notrig[proc].GetYaxis().SetRangeUser(0,1)
            heff_notrig[proc].Draw("e1")
        firstDraw = False
    else:
        heff_notrig[proc].Draw("same e1")

    if typ == "heavy" and proc == datakey:
        x,y,error_up_x,error_up_y,error_low_x,error_low_y = getxsecerr(heff_notrig_xsecup[proc],heff_notrig_xsecdw[proc],heff_notrig[proc])
        h_notrig_err_xsec[proc] = TGraphAsymmErrors(len(error_low_x),x,y,error_low_x,error_up_x,error_low_y,error_up_y)
        h_notrig_err_xsec[proc].SetLineStyle(kDashed);
        h_notrig_err_xsec[proc].SetLineColor(kGray);
        h_notrig_err_xsec[proc].SetFillColor(kGray);
        h_notrig_err_xsec[proc].SetFillStyle(3018)
        h_notrig_err_xsec[proc].SetMarkerSize(0)
        h_notrig_err_xsec[proc].Draw("e2same");
        if not added:
            leg.AddEntry(h_notrig_err_xsec[proc],"MC subtr. unc.","f")
            added = True
    #print("Writing %s to file %s"%((prefix.replace(analy,("%s_notrigm"%("e" if lep == "el" else "m")))+"_"+proc),utfil.GetName()))
    if proc in heff_notrig_xsecup.keys() and not uncert:
        heff_notrig_xsecup[proc].Write(prefix.replace(analy,("%s_notrigm"%("e" if lep == "el" else "m")))+"_xsecup_"+proc)
        heff_notrig_xsecdw[proc].Write(prefix.replace(analy,("%s_notrigm"%("e" if lep == "el" else "m")))+"_xsecdw_"+proc)
    heff_notrig[proc].Write(prefix.replace(analy,("%s_notrigm"%("e" if lep == "el" else "m")))+"_"+proc)


if typ == "light":
  leg.Draw()
  leg2.Draw()
elif not typ == "real":
    leg.Draw()
    c.Update()
utfil.Close()

if typ == "real":
  for c in canv:
    c.SaveAs("forZMET/%s.pdf"%c.GetTitle())
    c.SaveAs("forZMET/%s.png"%c.GetTitle())
else:
  c.SaveAs("forZMET/%s.pdf"%c.GetTitle())
  c.SaveAs("forZMET/%s.png"%c.GetTitle())

allbkg = []
for key in d_samp.keys():
    if key in ['mc16a','mc16e','mc16d','data15-16','data17','data18']: continue
    allbkg.append(key)

for key in d_samp.keys():
    if key in didplot and key in allbkg:
        allbkg.remove(key)
        

if len(allbkg) > 0:
  print("#"*100)
  print("Did not find following background(s): %s"%",".join(allbkg))
  print("#"*100)




if not typ == "real":
  cc = TCanvas("cc_check_%s_%s_%s%s"%(datakey,lep,typ,"_"+uncert if uncert else "_nominal"),"cc_check_%s_%s_%s_%s"%(datakey,lep,typ,"_"+uncert if uncert else "%s_nominal"%lep),1000,600)
  cc.cd()
  if not "eta" in typ and not vari == "BDT": cc.SetLogx()
  cc.SetGridx()
  cc.SetGridy()
  cc.SetLogy()
  for key in h_lT_trig.keys():
    if h_lT_trig_nosub[key][datakey].Integral() < 10: continue 
    print("Plotting %s"%key)
    h_lT_trig[key].Draw("hist")

    d_max = -1
    if datakey in h_lT_trig_nosub[key].keys():
      d_max = h_lT_trig_nosub[key][datakey].GetMaximum()
      h_lT_trig_nosub[key][datakey].Draw("same ep")

    if d_max > 0:
      h_lT_trig[key].SetMaximum(d_max*1.2)
      h_lT_trig[key].SetMinimum(10)
  diag_leg.Draw()
