from ROOT import *
import sys, os
from os import listdir
from os.path import isfile, join
import array
from math import sqrt
import socket

TH1.SetDefaultSumw2()
TH2.SetDefaultSumw2()

gROOT.SetBatch(True)

if socket.gethostname() == "quark":
  astyle = "/home/eirikgr/atlasrootstyle/"#atlasstyle-00-04-02/"#
else:
  astyle = "/mn/felt/u1/eirikgr/bin/atlasstyle-00-04-02/"                                                                                                                                                                           
gROOT.SetMacroPath(astyle)                                                                                                                                                                                     
gROOT.LoadMacro("AtlasUtils.C")                                                                                                                                                                                
gROOT.LoadMacro("AtlasStyle.C")                                                                                                                                                                                
gROOT.LoadMacro("AtlasLabels.C")                                                                                                                                                                               
SetAtlasStyle()  
gStyle.SetErrorX()
gStyle.SetEndErrorSize() 

gStyle.SetPaintTextFormat(".4f");


infile = sys.argv[1]
year = sys.argv[2]
analy = sys.argv[3]
tfile = TFile(infile,"READ")

if year == "2018": 
    mckey = "mc16e"
    datakey = "data18"
elif year == "2016": 
    year = "2015"
    mckey = "mc16a"
    datakey = "data15-16"
elif year == "2015": 
    mckey = "mc16a"
    datakey = "data15-16"
elif year == "2017": 
    mckey = "mc16d"
    datakey = "data17"

uncert_col = {}

uncert_col["EL_EFF_ID_TOTAL_1down"] =  kRed           
uncert_col["EL_EFF_ID_TOTAL_1up"] =    kRed           
uncert_col["EL_EFF_Iso_TOTAL_1down"] = kOrange           
uncert_col["EL_EFF_Iso_TOTAL_1up"] =          kOrange    
uncert_col["EL_EFF_Reco_TOTAL_1down"] =           kGreen +3
uncert_col["EL_EFF_Reco_TOTAL_1up"] =             kGreen +3

uncert_col["MUON_EFF_ISO_STAT_1down"] =           kCyan
uncert_col["MUON_EFF_ISO_STAT_1up"] =             kCyan
uncert_col["MUON_EFF_ISO_SYS_1down"] =            kOrange
uncert_col["MUON_EFF_ISO_SYS_1up"] =              kOrange
uncert_col["MUON_EFF_RECO_STAT_1down"] =          kBlue
uncert_col["MUON_EFF_RECO_STAT_1up"] =            kBlue
uncert_col["MUON_EFF_RECO_STAT_LOWPT_1down"] =    kOrange+7
uncert_col["MUON_EFF_RECO_STAT_LOWPT_1up"] =      kOrange+7
uncert_col["MUON_EFF_RECO_SYS_1down"] =           kGreen +3
uncert_col["MUON_EFF_RECO_SYS_1up"] =             kGreen +3
uncert_col["MUON_EFF_RECO_SYS_LOWPT_1down"] =     kGray
uncert_col["MUON_EFF_RECO_SYS_LOWPT_1up"] =       kGray


uncert_col["MUON_EFF_ISO_TOTAL_1down"] =            kOrange
uncert_col["MUON_EFF_ISO_TOTAL_1up"] =              kOrange
uncert_col["MUON_EFF_RECO_TOTAL_LOWPT_1down"] =    kOrange+7
uncert_col["MUON_EFF_RECO_TOTAL_LOWPT_1up"] =      kOrange+7
uncert_col["MUON_EFF_RECO_TOTAL_1down"] =           kGreen +3
uncert_col["MUON_EFF_RECO_TOTAL_1up"] =             kGreen +3

uncert_col["EL_EFF_Trigger_TOTAL_1down"] =        kPink +10
uncert_col["EL_EFF_Trigger_TOTAL_1up"] =          kPink +10

uncert_col["MUON_EFF_TrigStat_1down"] =           kPink +10
uncert_col["MUON_EFF_TrigStat_1up"] =             kPink +10
uncert_col["MUON_EFF_TrigSyst_1down"] =           kAzure +10
uncert_col["MUON_EFF_TrigSyst_1up"] =             kAzure +10

uncert_col["MUON_EFF_TRIG_TOTAL_1down"] =           kPink +10
uncert_col["MUON_EFF_TRIG_TOTAL_1up"] =           kPink +10

uncert_col["MUON_EFF_TOTAL_1down"] =           kAzure +7
uncert_col["MUON_EFF_TOTAL_1up"] =             kAzure +7

uncert_col["EL_EFF_TOTAL_1down"] =           kAzure +7
uncert_col["EL_EFF_TOTAL_1up"] =             kAzure +7

uncert_col["MUON_STAT_1up"] = kGreen
uncert_col["MUON_STAT_1down"] = kGreen

uncert_col["EL_STAT_1up"] = kGreen
uncert_col["EL_STAT_1down"] = kGreen

uncert = []

# uncert.append("EL_EFF_ID_TOTAL_1down")            
# uncert.append("EL_EFF_ID_TOTAL_1up")              
# uncert.append("EL_EFF_Iso_TOTAL_1down")           
# uncert.append("EL_EFF_Iso_TOTAL_1up")             
# uncert.append("EL_EFF_Reco_TOTAL_1down")          
# uncert.append("EL_EFF_Reco_TOTAL_1up")            
# uncert.append("MUON_EFF_ISO_STAT_1down")          
# uncert.append("MUON_EFF_ISO_STAT_1up")            
# uncert.append("MUON_EFF_ISO_SYS_1down")           
# uncert.append("MUON_EFF_ISO_SYS_1up")             
# uncert.append("MUON_EFF_RECO_STAT_1down")         
# uncert.append("MUON_EFF_RECO_STAT_1up")           
# uncert.append("MUON_EFF_RECO_STAT_LOWPT_1down")   
# uncert.append("MUON_EFF_RECO_STAT_LOWPT_1up")     
# uncert.append("MUON_EFF_RECO_SYS_1down")          
# uncert.append("MUON_EFF_RECO_SYS_1up")            
# uncert.append("MUON_EFF_RECO_SYS_LOWPT_1down")    
# uncert.append("MUON_EFF_RECO_SYS_LOWPT_1up")      
# uncert.append("EL_EFF_Trigger_TOTAL_1down")       
# uncert.append("EL_EFF_Trigger_TOTAL_1up")         
# uncert.append("MUON_EFF_TrigStat_1down")          
# uncert.append("MUON_EFF_TrigStat_1up")            
# uncert.append("MUON_EFF_TrigSyst_1down")          
# uncert.append("MUON_EFF_TrigSyst_1up")   
trigger_setup = {}
if analys == "FCL":
  trigger_setup = {"e_med":{"triggers":["e24_lhmedium_L1EM20VH","e24_lhmedium_nod0_L1EM20VH","e60_lhmedium_nod0"],"col":kOrange-4,"marker":22},
                   "e_loose":{"triggers":["e120_lhloose","e140_lhloose_nod0"],"col":kYellow-7,"marker":23},
                   "e_tight_iso":{"triggers":["e24_lhtight_nod0_ivarloose","e26_lhtight_nod0_ivarloose"],"col":kRed-3,"marker":21},
                   "m_med":{"triggers":["mu26_ivarmedium"],"col":kOrange-4,"marker":22},
                   "m_loose":{"triggers":["mu20_iloose_L1MU15"],"col":kYellow-7,"marker":23},
                   "m_no":{"triggers":["mu50"],"col":kGreen-7,"marker":21},
                   "m_notrigm":{"triggers":["lepnotrig"],"col":kGreen-7,"marker":21},
                   "e_notrigm":{"triggers":["lepnotrig"],"col":kGreen-7,"marker":21}}
elif analy == "ZMET":
  trigger_setup = {"e_loose":{"triggers":["2e12_lhloose_L12EM10VH","2e17_lhvloose_nod0","2e24_lhvloose_nod0","2e17_lhvloose_nod0_L12EM15VHI"],"col":kYellow-7,"marker":23},
                   "m_med":{"triggers":["mu26_imedium","mu26_ivarmedium"],"col":kOrange-4,"marker":22},
                   "m_no":{"triggers":["mu50"],"col":kGreen-7,"marker":21}}
if analy == "2L2J":
  trigger_setup = {"el":{"triggers":[""],"col":kBlack,"marker":22},
                   "mu":{"triggers":[""],"col":kBlack,"marker":22}}
  

## REAL
regions = ["FAKE2L21","REAL2L04"]#,"FAKE2L05","FAKE2L21","FAKE2L23"]
#region = "FAKE2L05"

if len(regions) > 1: gROOT.SetBatch(True)

finalfile = infile.replace(".root","_final_v2.root")
tfinalfile = TFile(finalfile,"UPDATE")

c = {}
leg = {}
for r in regions:
  c[r] = {}
  leg[r] = {}
for region in regions:

    h_xsecup  = {}
    h_xsecdw  = {}
    h_nominal  = {}
    h_fraction_nominal  = {}
    h_syst     = {}
    h_fraction = {}

    for key in trigger_setup.keys():
        if region in ["FAKE2L05","FAKE2L23"] and key.startswith("m"):
            print("Skipping %s for key %s"%(region,key))
            continue

        if region in ["FAKE2L05","FAKE2L23","FAKE2L21"]: 
            var = "pt"
        else: 
            var = "pt_eta"
        #var = "pt"

        hname = "eff_%s_%s_%s_%s_%s"%(var,key,year,region,mckey if region in ["REAL2L04","FAKE2L05"] else datakey)        
        objPtr = tfile.Get(hname)
        if objPtr: 
          h_nominal[region] = {}
          h_fraction_nominal[region] = {}
          h_nominal[region][key] = objPtr.Clone(objPtr.GetName())
          h_nominal[region][key].SetDirectory(0)
          h_fraction_nominal[region][key] = h_nominal[region][key].Clone(h_nominal[region][key].GetName().replace("eff","frac"))
          h_fraction_nominal[region][key].Divide(h_nominal[region][key],h_nominal[region][key])
          for ibin in range(1,h_fraction_nominal[region][key].GetNbinsX()+1):
            h_fraction_nominal[region][key].SetBinError(ibin,0.0)
          h_fraction_nominal[region][key].SetDirectory(0)
        else:
            print("Could not find %s in %s"%(hname,infile))

        # Thes regions do not have systematics
        if region in ["FAKE2L23","FAKE2L21"] and key in h_nominal.keys():
            h_nominal[region][key].Write("eff_%s_%s_%s_%s_%s"%(var,key,year,region,mckey if region in ["REAL2L04","FAKE2L05"] else datakey))
            if region == "FAKE2L21":
                hname_xsec = "eff_%s_%s_%s_%s_xsecup_%s"%(var,key,year,region,datakey)        
                objPtr = tfile.Get(hname_xsec)
                if objPtr:
                    h_xsecup[key] = objPtr.Clone(objPtr.GetName())
                    h_xsecup[key].SetDirectory(0)
                    h_xsecup[key].Write("eff_%s_%s_%s_%s_xsecup_%s"%(var,key,year,region,datakey))
                hname_xsec = "eff_%s_%s_%s_%s_xsecdw_%s"%(var,key,year,region,datakey)        
                objPtr = tfile.Get(hname_xsec)
                if objPtr:
                    h_xsecdw[key] = objPtr.Clone(objPtr.GetName())
                    h_xsecdw[key].SetDirectory(0)
                    h_xsecdw[key].Write("eff_%s_%s_%s_%s_xsecdw_%s"%(var,key,year,region,datakey))
            continue

        for u in uncert:
            if key.startswith("e") and "MUON" in u: continue
            if key.startswith("m") and "EL" in u: continue
            hname = "eff_%s_%s_%s_%s_%s_%s"%(var,key,year,region,u,mckey)        
            objPtr = tfile.Get(hname)
            if objPtr: 
                if not key in h_syst.keys():
                    h_syst[key] = {}
                    h_fraction[key] = {}
                h_syst[key][u] = objPtr.Clone(objPtr.GetName())
                h_syst[key][u].SetDirectory(0)
            else:
                print("Could not find %s in %s"%(hname,infile))
                continue
            h_fraction[key][u] = h_syst[key][u].Clone(h_syst[key][u].GetName().replace("eff","frac"))
            h_fraction[key][u].Divide(h_fraction[key][u],h_nominal[region][key])

        if not key in h_fraction.keys(): continue

        for u in ["MUON_STAT_1down","MUON_STAT_1up","EL_STAT_1down","EL_STAT_1up"]:
            if key.startswith("m") and "EL" in u: continue
            if key.startswith("e") and "MUON" in u: continue
            h_syst[key][u] = h_nominal[region][key].Clone(h_nominal[region][key].GetName().replace(region,region+"_"+u))
            if type(h_syst[key][u]) is TH2F:
                for xbin in range(1,h_syst[key][u].GetNbinsX()+1):
                    for ybin in range(1,h_syst[key][u].GetNbinsY()+1):
                        h_syst[key][u].AddBinContent(h_syst[key][u].GetBin(xbin,ybin),((-1) if "1down" in u else (1))*h_nominal[region][key].GetBinError(xbin,ybin))
            else:
                for ibin in range(1,h_syst[key][u].GetNbinsX()+1):
                    h_syst[key][u].AddBinContent(ibin,((-1) if "1down" in u else (1))*h_nominal[region][key].GetBinError(ibin))

            h_fraction[key][u] = h_syst[key][u].Clone(h_syst[key][u].GetName().replace("eff","frac"))
            h_fraction[key][u].Divide(h_fraction[key][u],h_nominal[region][key])

        ## Add SYST and STAT in quadrature
        for unc in ["MUON_EFF_TRIG_TOTAL_1down","MUON_EFF_ISO_TOTAL_1down",
                    "MUON_EFF_RECO_TOTAL_1down","MUON_EFF_RECO_TOTAL_LOWPT_1down",
                    "MUON_EFF_ISO_TOTAL_1up","MUON_EFF_RECO_TOTAL_1up",
                    "MUON_EFF_RECO_TOTAL_LOWPT_1up","MUON_EFF_TRIG_TOTAL_1up"]:

            if key.startswith("e"): continue

            statunc = (unc.replace("TOTAL","STAT") if not "TRIG" in unc else unc.replace("TRIG_TOTAL","TrigStat"))
            systunc = (unc.replace("TOTAL","SYS") if not "TRIG" in unc else unc.replace("TRIG_TOTAL","TrigSyst"))

            
            h_fraction[key][unc] = h_fraction[key][statunc].Clone(h_fraction[key][statunc].GetName().replace("STAT","TOTAL"))

            print("statkey = %s, newkey = %s"%(statunc,unc))

            if type(h_fraction[key][unc]) is TH2F:
                for xbin in range(1,h_fraction[key][unc].GetNbinsX()+1):
                    for ybin in range(1,h_fraction[key][unc].GetNbinsY()+1):
                        if "1down" in unc :
                            syst_val = 1-h_fraction[key][systunc].GetBinContent(xbin,ybin)
                            stat_val = 1-h_fraction[key][statunc].GetBinContent(xbin,ybin)
                            h_fraction[key][unc].SetBinContent(h_fraction[key][unc].GetBin(xbin,ybin),1-sqrt((syst_val*syst_val)+(stat_val*stat_val)))
                        if "1up" in unc :
                            syst_val = h_fraction[key][systunc].GetBinContent(xbin,ybin)-1
                            stat_val = h_fraction[key][statunc].GetBinContent(xbin,ybin)-1
                            h_fraction[key][unc].SetBinContent(h_fraction[key][unc].GetBin(xbin,ybin),1+sqrt((syst_val*syst_val)+(stat_val*stat_val)))
            else:
                for xbin in range(1,h_fraction[key][unc].GetNbinsX()+1): 
                    if "1down" in unc :
                        syst_val = 1-h_fraction[key][systunc].GetBinContent(xbin)
                        stat_val = 1-h_fraction[key][statunc].GetBinContent(xbin)
                        h_fraction[key][unc].SetBinContent(xbin,1-sqrt((syst_val*syst_val)+(stat_val*stat_val)))
                    if "1up" in unc :
                        syst_val = h_fraction[key][systunc].GetBinContent(xbin)-1
                        stat_val = h_fraction[key][statunc].GetBinContent(xbin)-1
                        h_fraction[key][unc].SetBinContent(xbin,1+sqrt((syst_val*syst_val)+(stat_val*stat_val)))

        ## Add RECO and ISO lineraly
        for unc in ["MUON_EFF_ISO_TOTAL","MUON_EFF_RECO_TOTAL",
                    "MUON_EFF_RECO_TOTAL_LOWPT","EL_EFF_ID_TOTAL",
                    "EL_EFF_Iso_TOTAL","EL_EFF_Reco_TOTAL"]:

            if key.startswith("e") and "MUON" in unc: continue
            if key.startswith("m") and "EL" in unc: continue

            totkey = unc+"_1down"
            newkey = unc.split("_")[0]+"_"+unc.split("_")[1]+"_"+"_".join(unc.split("_")[3:])+"_1down"

            print("totkey = %s, newkey = %s"%(totkey,newkey))

            if not newkey in h_fraction[key].keys():
              h_fraction[key][newkey] = h_fraction[key][totkey].Clone(h_fraction[key][totkey].GetName().replace("STAT","TOTAL")+"_1down")
              h_fraction[key][newkey].Reset("ICE")
            newkeyup = newkey.replace("_1down","_1up")
            if not newkeyup in h_fraction[key].keys():
              h_fraction[key][newkeyup] = h_fraction[key][totkey].Clone(h_fraction[key][totkey].GetName().replace("STAT","TOTAL")+"_1up")
              h_fraction[key][newkeyup].Reset("ICE")

            if type(h_fraction[key][newkey]) is TH2F:
                for xbin in range(1,h_fraction[key][newkey].GetNbinsX()+1):
                    for ybin in range(1,h_fraction[key][newkey].GetNbinsY()+1):
                        if h_fraction[key][totkey].GetBinContent(xbin,ybin) < 1.0:
                          h_fraction[key][newkey].AddBinContent(h_fraction[key][newkey].GetBin(xbin,ybin),(1-abs(h_fraction[key][totkey].GetBinContent(xbin,ybin))))
                        else:
                          h_fraction[key][newkeyup].AddBinContent(h_fraction[key][newkey].GetBin(xbin,ybin),(1-abs(h_fraction[key][totkey].GetBinContent(xbin,ybin))))
            else:
                for xbin in range(1,h_fraction[key][newkey].GetNbinsX()+1):
                    if h_fraction[key][totkey]:
                      print("%s %s fraction down = %.4f"%(totkey,key,(h_fraction[key][totkey].GetBinContent(xbin))))
                      if (h_fraction[key][totkey].GetBinContent(xbin)) < 1.0:
                        h_fraction[key][newkey].AddBinContent(xbin,(1-abs(h_fraction[key][totkey].GetBinContent(xbin))))
                      else:
                        h_fraction[key][newkeyup].AddBinContent(xbin,(abs(h_fraction[key][totkey].GetBinContent(xbin))-1))

            totkey = unc+"_1up"
            #newkey = unc.split("_")[0]+"_"+unc.split("_")[1]+"_"+"_".join(unc.split("_")[3:])+"_1up"

            if type(h_fraction[key][newkey]) is TH2F:
                for xbin in range(1,h_fraction[key][newkey].GetNbinsX()+1):
                    for ybin in range(1,h_fraction[key][newkey].GetNbinsX()+1):
                        if h_fraction[key][totkey].GetBinContent(xbin,ybin) > 1.0:
                          h_fraction[key][newkeyup].AddBinContent(h_fraction[key][newkey].GetBin(xbin,ybin),(abs(h_fraction[key][totkey].GetBinContent(xbin,ybin))-1))
                        else:
                          h_fraction[key][newkey].AddBinContent(h_fraction[key][newkey].GetBin(xbin,ybin),(1-abs(h_fraction[key][totkey].GetBinContent(xbin,ybin))))
            else:
                for xbin in range(1,h_fraction[key][newkey].GetNbinsX()+1):
                    if h_fraction[key][totkey]:
                      print("%s %s fraction up = %.4f"%(totkey,key,h_fraction[key][totkey].GetBinContent(xbin)))
                      if h_fraction[key][totkey].GetBinContent(xbin) > 1.0:
                        h_fraction[key][newkeyup].AddBinContent(xbin,(abs(h_fraction[key][totkey].GetBinContent(xbin))-1))
                      else:
                        h_fraction[key][newkey].AddBinContent(xbin,(1-abs(h_fraction[key][totkey].GetBinContent(xbin))))

        # Finally getting the totals with respect to 1!
        for unc in ['EL_EFF_TOTAL_1up', 'EL_EFF_TOTAL_1down','MUON_EFF_TOTAL_1up', 'MUON_EFF_TOTAL_1down']:
          if key.startswith("e") and "MUON" in unc: continue
          if key.startswith("m") and "EL" in unc: continue

          if type(h_fraction[key][unc]) is TH2F:
            for xbin in range(1,h_fraction[key][unc].GetNbinsX()+1):
              for ybin in range(1,h_fraction[key][unc].GetNbinsY()+1):
                ibin = h_fraction[key][unc].GetBin(xbin,ybin)
                if "up" in unc:
                  h_fraction[key][unc].SetBinContent(ibin,1+h_fraction[key][unc].GetBinContent(xbin,ybin))
                elif "down" in unc:
                  h_fraction[key][unc].SetBinContent(ibin,1-abs(h_fraction[key][unc].GetBinContent(xbin,ybin)))
          else:
            for xbin in range(1,h_fraction[key][unc].GetNbinsX()+1):
              if "up" in unc:
                print("UP: Setting bin %i content to %.2f"%(xbin,1+h_fraction[key][unc].GetBinContent(xbin)))
                h_fraction[key][unc].SetBinContent(xbin,1+h_fraction[key][unc].GetBinContent(xbin))
              elif "down" in unc:
                print("Unc = %s and cont = %.4f"%(unc,h_fraction[key][unc].GetBinContent(xbin)))
                print("DWON: Setting bin %i content to %.2f"%(xbin,1-abs(h_fraction[key][unc].GetBinContent(xbin))))
                h_fraction[key][unc].SetBinContent(xbin,1-abs(h_fraction[key][unc].GetBinContent(xbin)))


        # Adding the new systeamtics
        new_uncert = uncert[:]
        for u in ['EL_STAT_1up', 'EL_STAT_1down', 'EL_EFF_TOTAL_1up', 'EL_EFF_TOTAL_1down',
                  'MUON_STAT_1up', 'MUON_STAT_1down', 'MUON_EFF_TOTAL_1up', 'MUON_EFF_TOTAL_1down',
                  'MUON_EFF_TRIG_TOTAL_1up', 'MUON_EFF_TRIG_TOTAL_1down']:
          new_uncert.append(u)


        print("key = %s"%key)
        for unc in h_fraction[key].keys():
            print("unc = %s"%(unc))
        if type(h_fraction_nominal[region][key]) is not TH2F:
            h_fraction_nominal[region][key].SetMarkerStyle(20)
            h_fraction_nominal[region][key].SetMarkerColor(kBlack)
            h_fraction_nominal[region][key].GetYaxis().SetRangeUser(0.8,1.2)
            h_fraction_nominal[region][key].GetYaxis().SetTitle("Fake Rate")
            h_fraction_nominal[region][key].GetXaxis().SetTitle("p_{T} [GeV]")
            h_fraction_nominal[region][key].GetXaxis().SetMoreLogLabels()
        h_nominal[region][key].Write("eff_%s_%s_%s_%s_%s"%(var,key,year,region,mckey))
        if not key in c[region].keys():
          c[region][key] = []
          leg[region][key] = []
        nu = -1
        legadded = []
        for unc in new_uncert:
          
            nu += 1;
            print("key = %s, region = %s and unc = %s"%(key,region,unc))
            if key.startswith("e") and "MUON" in unc: continue
            if key.startswith("m") and "EL" in unc: continue
            if type(h_fraction_nominal[region][key]) is TH2F:
                c[region][key].append(TCanvas("c_%s_%s_%s_%s"%(region,key,unc,year),"c_%s_%s_%s_%s"%(region,key,unc,year),1000,600))
                c[region][key][-1].SetLeftMargin(0.1)
                c[region][key][-1].SetRightMargin(0.15)
            elif len(c[region][key]) == 0:
                c[region][key].append(TCanvas("c_%s_%s_%s_%s"%(region,key,year,"all"),"c_%s_%s_%s_%s"%(region,key,year,"all"),1000,600))
                leg[region][key].append(TLegend(0.18,0.72,0.41,0.94))
                leg[region][key][-1].SetBorderSize(0)
                leg[region][key][-1].SetTextFont(42)
                leg[region][key][-1].SetTextSize(0.03)
                leg[region][key][-1].SetTextColor(1)
                leg[region][key][-1].SetFillStyle(0)
                leg[region][key][-1].SetLineColor(0)
            c[region][key][-1].cd()
            c[region][key][-1].SetLogx()
            c[region][key][-1].SetGridx()
            c[region][key][-1].SetGridy()
            h_fraction_nominal[region][key].SetMinimum(0.8)
            h_fraction_nominal[region][key].SetMaximum(1.2)
            if type(h_fraction_nominal[region][key]) is not TH2F and nu == 0: 
              h_fraction_nominal[region][key].Draw("P")
            if unc in uncert_col.keys():
                if type(h_fraction[key][unc]) is TH2F: 
                    h_fraction[key][unc].GetXaxis().SetTitle("p_{T} [GeV]")
                    h_fraction[key][unc].GetYaxis().SetTitle("|#eta|")
                    for xbin in range(1,h_fraction[key][unc].GetNbinsX()+1):
                      for ybin in range(1,h_fraction[key][unc].GetNbinsY()+1):
                        if h_fraction[key][unc].GetBinContent(xbin,ybin) > 2 or h_fraction[key][unc].GetBinContent(xbin,ybin) < 0:
                          h_fraction[key][unc].SetBinContent(xbin,ybin,0.0)
                    if "1up" in unc:
                      h_fraction[key][unc].SetMaximum(1.1)
                      h_fraction[key][unc].SetMinimum(1.0)
                    elif "1down" in unc:
                      h_fraction[key][unc].SetMaximum(1.0)
                      h_fraction[key][unc].SetMinimum(0.9)
                    #h_fraction[key][unc].GetZaxis().SetLimits(0.5,1.5)
                    h_fraction[key][unc].Draw("colz text")
                else:
                    legstr = "_".join(unc.split("_")[:-1])
                    if not legstr in legadded:
                      leg[region][key][-1].AddEntry(h_fraction[key][unc],legstr,"l")
                      legadded.append(legstr)
                    h_fraction[key][unc].SetLineColor(uncert_col[unc])
                    h_fraction[key][unc].Draw("same hist")
            else:
                print("Could not find color for %s"%unc)
            if unc in ['EL_STAT_1up', 'EL_STAT_1down', 'EL_EFF_TOTAL_1up', 'EL_EFF_TOTAL_1down',
                       'MUON_STAT_1up', 'MUON_STAT_1down', 'MUON_EFF_TOTAL_1up', 'MUON_EFF_TOTAL_1down',
                       'MUON_EFF_TRIG_TOTAL_1up', 'MUON_EFF_TRIG_TOTAL_1down','EL_EFF_Trigger_TOTAL_1up','EL_EFF_Trigger_TOTAL_1down']:
                print("Writing %s to %s"%("frac_%s_%s_%s_%s_%s_%s"%(var,key,year,region,unc,mckey),finalfile))
                h_fraction[key][unc].Write("frac_%s_%s_%s_%s_%s_%s"%(var,key,year,region,unc,mckey))
            #if len(c[region][key]) > 5: break
        if type(h_fraction_nominal[region][key]) is not TH2F: h_fraction_nominal[region][key].Draw("same P")

        ni = 0
        for canv in c[region][key]:
            try:
              leg[region][key][ni].Draw()
            except:
              print("No legend available for 2D hist")
            canv.Update()
            canv.SaveAs("for2L0Jpaper/"+canv.GetTitle()+".pdf")
            canv.SaveAs("for2L0Jpaper/"+canv.GetTitle()+".png")
            ni += 1
tfinalfile.Close()

        
