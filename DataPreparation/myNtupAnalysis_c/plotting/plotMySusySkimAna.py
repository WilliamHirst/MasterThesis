from ROOT import *#TFile, TH1, TH2, THStack, TH1F,TLegend, TCanvas, gROOT, gStyle, TH2F, TGraph, TPad, TColor, kWhite, kBlack, kDashed, kYellow, kMagenta, TGraphAsymmErrors, kGray, TLatex, Double
from optparse import OptionParser
from os import listdir
from os.path import isfile, join
from samples import configure_samples
import array
from math import *
import sys,os
import glob
import ctypes
# python -i plotMySusySkimAna.py --doMC 1 --inputdir Fri_Jun_1_2018_2L_NTUP/ --histname h_lep_mT2_TCF_nT_EMOS --doSub 0
doLOG = False

astyle = "/home/eirikgr/atlasstyle-00-04-02/"#
astyle ="/mn/felt/u1/eirikgr/bin/atlasstyle-00-04-02/"
astyle ="/home/eirikgr/atlasrootstyle/"
gROOT.SetMacroPath(astyle)
gROOT.LoadMacro("AtlasUtils.C")
gROOT.LoadMacro("AtlasStyle.C")
gROOT.LoadMacro("AtlasLabels.C")
#SetAtlasStyle()  

#gStyle.SetErrorX(1) 

d_samp,d_type, d_reg = configure_samples()

gStyle.SetPaintTextFormat(".4f");

metadata = {}

TH1.SetDefaultSumw2()
TH2.SetDefaultSumw2()

canv = []

lumi = {}
lumi["data15-16"] = 36.2
lumi["data1516"] = 36.2
lumi["data17"] = 44.3
lumi["data18"] = 59.9
lumi["alldata"] = 36.2 + 44.3 + 59.9


parser = OptionParser()
parser.add_option("--doMC", help="do mc", default=0)
parser.add_option("--doSub", help="subtract MC", default=0)
parser.add_option("--inputdir", help="Location of files", default="")
parser.add_option("--inputdirFAKE", help = "Location of fake files", default = "")
parser.add_option("--doNEV", help="do NEV", default="")
parser.add_option("--datakey", help="data15-16, data17 or data18", default="")
parser.add_option("--histname", help="name of histogram", default="")
parser.add_option("--makeMMfile", help="make MM file", default=0)
parser.add_option("--doEff", help="make fake/real eff plots", default=0)
parser.add_option("--setBATCH", help="turn of batch mode", default=0)
parser.add_option("--estfake", help="estimate fakes from data", default=0)
parser.add_option("--subBkgs", help="which backgrounds to subtract", default="ALL")
parser.add_option("--allData", help="if combine 15-18", default=0)
(options, args) = parser.parse_args()

if int(options.setBATCH):
    gROOT.SetBatch(True)

alldata = int(options.allData)
doMC = int(options.doMC)
inputdir   = str(options.inputdir)
datakey = str(options.datakey)
histname      = str(options.histname).split(",")
doSub = int(options.doSub)
estfake = int(options.estfake)
inputdir_FAKE = str(options.inputdirFAKE)
if str(options.subBkgs) != "ALL":
    subBkgs = str(options.subBkgs).split(",")
else:
    subBkgs = ["ALL"] 

fake_hist_name = "="
data_hist_name = "="

comb1516 = True

varmap = {"pT":"pT",
          "mt2":"MT2",
          "pT1":"pT1",
          "pT2":"pT2",
          "mll":"MLL",
          "met":"MET"}
if estfake and inputdir_FAKE:
    FAKEfiles   = [join(inputdir_FAKE,f) for f in listdir(inputdir_FAKE) if (isfile(join(inputdir_FAKE, f)) and f.endswith("_HISTOGRAMS_0.root"))]
hsumMC = {}
h_eff_sub_pEff = {}
h_eff_sub = {}
h_eff_sub_xsecup = {}
h_eff_sub_xsecdw = {}
h_nT_sub = {}
h_nT_sub_xsecup = {}
h_nT_sub_xsecdw = {}
h_nL_sub = {}
h_nL_sub_xsecup = {}
h_nL_sub_xsecdw = {}
h_eff_err_xsec = {}
h_allMC   = {}
h_allMC_xsecup   = {}
h_allMC_xsecdw   = {}

#MCfiles   = [join(inputdir_MC,f) for f in listdir(inputdir_MC) if (isfile(join(inputdir_MC, f)) and f.endswith("processed_HISTOGRAMS.root"))]
#DATAfiles = [join(inputdir_MC,f) for f in listdir(inputdir_DATA) if (isfile(join(inputdir_DATA, f)) and f.endswith("processed_HISTOGRAMS.root"))]    

uncert = []
uncert.append("EL_EFF_ID_TOTAL_1down")            
uncert.append("EL_EFF_ID_TOTAL_1up")              
uncert.append("EL_EFF_Iso_TOTAL_1down")           
uncert.append("EL_EFF_Iso_TOTAL_1up")             
uncert.append("EL_EFF_Reco_TOTAL_1down")          
uncert.append("EL_EFF_Reco_TOTAL_1up")            
uncert.append("MUON_EFF_ISO_STAT_1down")          
uncert.append("MUON_EFF_ISO_STAT_1up")            
uncert.append("MUON_EFF_ISO_SYS_1down")           
uncert.append("MUON_EFF_ISO_SYS_1up")             
uncert.append("MUON_EFF_RECO_STAT_1down")         
uncert.append("MUON_EFF_RECO_STAT_1up")           
uncert.append("MUON_EFF_RECO_STAT_LOWPT_1down")   
uncert.append("MUON_EFF_RECO_STAT_LOWPT_1up")     
uncert.append("MUON_EFF_RECO_SYS_1down")          
uncert.append("MUON_EFF_RECO_SYS_1up")            
uncert.append("MUON_EFF_RECO_SYS_LOWPT_1down")    
uncert.append("MUON_EFF_RECO_SYS_LOWPT_1up")      
uncert.append("EL_EFF_TriggerEff_TOTAL_1down")    
uncert.append("EL_EFF_TriggerEff_TOTAL_1up")      
uncert.append("EL_EFF_Trigger_TOTAL_1down")       
uncert.append("EL_EFF_Trigger_TOTAL_1up")         
uncert.append("MUON_EFF_TrigStat_1down")          
uncert.append("MUON_EFF_TrigStat_1up")            
uncert.append("MUON_EFF_TrigSyst_1down")          
uncert.append("MUON_EFF_TrigSyst_1up")            
 
def getHistWithoutSyst(hname):
    for u in uncert:
        if u in hname:
            return hname.replace("_"+u,"")
    return hname

text_size = 0.45
def myText(x, y, text, tsize=0.05, color=kBlack, angle=0) :
  
  l = TLatex()
  l.SetTextSize(tsize)
  l.SetNDC()
  l.SetTextColor(color)
  l.SetTextAngle(angle)
  l.DrawLatex(x,y,'#bf{' + text + '}')
  l.SetTextFont(4)


def MakeResidPlot(hmc,hdata,resid,resid_up,resid_dwn,rangeX,asymErr = 0,h_fake = 0):

    hN = hdata.GetName()

    resid = hdata.Clone("resid")
    #resid.Reset()
    #resid_tax=resid.GetXaxis()
    #resid_tax.SetTitle(hmc.GetXaxis().GetTitle())
    resid.GetYaxis().SetTitle("Data/MC");

    # ### construct and draw the hmax hmin (errorband)
    h_max = TH1F()
    h_max = hmc.Clone()
    h_max.Reset()
    h_min = TH1F()
    h_min = hmc.Clone()
    h_min.Reset()



    """ get the max value of the error (symmetric at the moment, so just get the max value)
    in order to set the range of the plot to cover the uncertainty at least. Cannot include
    all the points as the range would be tooo large """

    err_max = 0.

    Nbins = hmc.GetNbinsX()
    for i in range(0,Nbins+1):
        nom = hmc.GetBinContent(i+1)
        
        if asymErr != 0:
            # print "nom",nom
            # print "asymErr.GetErrorYhigh(%i) = %.2f" %(i,asymErr.GetErrorYhigh(i))
            # print "asymErr.GetErrorYlow(%i) = %.2f" %(i,asymErr.GetErrorYlow(i))
            # print "hmc.GetBinContent(%i) = %.2f" %(i+1,hmc.GetBinContent(i+1))
            h_max.SetBinContent(i+1,nom + asymErr.GetErrorYhigh(i))
            if asymErr.GetErrorYlow(i) > 0:
                h_min.SetBinContent(i+1,nom - asymErr.GetErrorYlow(i))
            else:
                h_min.SetBinContent(i+1,nom - h_fake.GetBinContent(i+1))
        else:
            err = hmc.GetBinError(i+1)
            h_max.SetBinContent(i+1,nom+err)
            h_min.SetBinContent(i+1,nom-err)

        temp_err=0.
        try:
            temp_err = err/hmc.GetBinContent(i+1)
        except:
            temp_err = 0.
        if temp_err > err_max:
            err_max = temp_err

    # to avoid too high valules
    if err_max > 3:
        err_max = 3


    # print "Resid:",resid.GetBinWidth(1)
    # print "hMC:",hmc.GetBinWidth(1)
    # print "hDATA:",hdata.GetBinWidth(1)
    # print "h_max:",h_max.GetBinWidth(1)
    # print "h_min:",h_min.GetBinWidth(1)

    #resid.Add(hdata)
    #resid.Divide(resid,hmc,1.,1.,'b')

    #asym = TGraphAsymmErrors();
    #print "Dividing resid = %.2f on hmc = %.2f" %(resid.Integral(),hmc.Integral())
    #resid.Divide(resid,hmc,"cl=0.683 b(1,1) mode")
    resid.Divide(resid,hmc)
    #for i in range(0,resid.GetNbinsX()+1):
    #    resid.SetBinError(i+1,asym.GetErrorY(i))

    resid.GetXaxis().SetTitleSize(0.13)
    resid.GetXaxis().SetTitleOffset(1.13)
    #resid.GetXaxis().SetLabelSize(0.1)
    resid.GetXaxis().SetLabelSize(0.105)
    resid.GetYaxis().SetTitleOffset(0.5);
    resid.GetYaxis().SetTitleSize(0.4);
    resid.GetYaxis().SetLabelSize(0.09)
    resid.GetYaxis().SetTitle("Data/MC");
    resid.SetMarkerStyle(20)
    resid.SetMarkerSize(0.75)


    maxi = resid.GetBinContent(resid.GetMaximumBin())

    if maxi > 3.0:
        resid.SetMaximum(maxi*1.5);
    else: resid.SetMaximum(3.);
    resid.SetMinimum(0.01);
    
    if rangeX > 0:
        resid.GetXaxis().SetRangeUser(0,rangeX)

        
    #h_max.Divide(h_max,hmc,1.,1.,'b')
    #h_min.Divide(h_min,hmc,1.,1.,'b')

    #asym = TGraphAsymmErrors();
    #print "Dividing h_max = %.2f on hmc = %.2f" %(h_max.Integral(),hmc.Integral())
    #h_max.Divide(h_max,hmc,"cl=0.683 b(1,1) mode")
    h_max.Divide(h_max,hmc)
    #for i in range(0,h_max.GetNbinsX()+1):
    #    h_max.SetBinError(i+1,asym.GetErrorY(i))

    #asym = TGraphAsymmErrors();
    #print "Dividing hmin = %.2f on hmc = %.2f" %(h_min.Integral(),hmc.Integral())
    #h_min.Divide(h_min,hmc,"cl=0.683 b(1,1) mode")
    h_min.Divide(h_min,hmc)
    #for i in range(0,h_min.GetNbinsX()+1):
    #    h_min.SetBinError(i+1,asym.GetErrorY(i))

    resid_up = h_max.Clone('resid_up')
    resid_dwn = h_min.Clone('resid_dwn')

    resid_up.SetLineStyle(kDashed);
    resid_up.SetLineColor(kBlack);
    resid_up.SetFillColor(kGray);
    resid_up.SetFillStyle(3001)
    resid_dwn.SetLineStyle(kDashed);
    resid_dwn.SetLineColor(kBlack);
    resid_dwn.SetFillColor(10);
    resid_dwn.SetFillStyle(1001)

    return resid,resid_up,resid_dwn


def set_palette(name='palette', ncontours=999):
    """Set a color palette from a given RGB list
    stops, red, green and blue should all be lists of the same length
    see set_decent_colors for an example"""

    if name == "gray" or name == "grayscale":
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [1.00, 0.84, 0.61, 0.34, 0.00]
        green = [1.00, 0.84, 0.61, 0.34, 0.00]
        blue  = [1.00, 0.84, 0.61, 0.34, 0.00]
    # elif name == "whatever":
        # (define more palettes)
    else:
        # default palette, looks cool
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [0.00, 0.00, 0.87, 1.00, 0.51]
        green = [0.00, 0.81, 1.00, 0.20, 0.00]
        blue  = [0.51, 1.00, 0.12, 0.00, 0.00]

    s = array.array('d', stops)
    r = array.array('d', red)
    g = array.array('d', green)
    b = array.array('d', blue)

    npoints = len(s)
    TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
    gStyle.SetNumberContours(ncontours)

def doRebin2D(old,xbins=array.array('f',[0,10,20,25,30,40,50,60,80,100,200]),ybins=array.array('f',[0,10,20,25,30,40,50,60,80,100,200])):
    # create a new TH2 with your bin arrays spec
    h = TH2F(old.GetName()+"_oldrebin",old.GetTitle(),len(xbins)-1,xbins,len(ybins)-1,ybins)
    xaxis = old.GetXaxis()
    yaxis = old.GetYaxis()
    for j in range(1,yaxis.GetNbins()):
        for i in range(1,xaxis.GetNbins()+2):
            h.Fill(xaxis.GetBinCenter(i),yaxis.GetBinCenter(j),old.GetBinContent(i,j));
        if "TCF" in old.GetName():
            print("a")
            h.GetYaxis().SetBinLabel(j,yaxis.GetBinLabel(j))
    for i in range(h.GetNbinsX()+1):
        for j in range(h.GetNbinsY()+1):
            if h.GetBinContent(i,j) > 0:
                h.SetBinError(i,j,sqrt(h.GetBinContent(i,j)))
            else:
                h.SetBinError(i,j,0)
    return h

def fillMetadata(h):
    metadata[h]["xmin"] = 0
    metadata[h]["xmax"] = 1000

    if "eta" in h: metadata[h]["xtit"] = "|#eta|"
    elif ("pT" in h or "pt" in h) and not "pT1" in h and not "pT2" in h:
        metadata[h]["xtit"] = "p_{T} [GeV]"
        if "ALL_M" in h and "FAKE2L21" in h:
            metadata[h]["xmin"] = 0
            metadata[h]["xmax"] = 999
        else:
            metadata[h]["xmin"] = 0
            metadata[h]["xmax"] = 999
    elif "lep_met_" in h:
        metadata[h]["xmin"] = 0
        metadata[h]["xmax"] = 400
        metadata[h]["xtit"] = "E_{T}^{miss} [GeV]"
    elif "lep_metsig_" in h:
        metadata[h]["xmin"] = 0
        metadata[h]["xmax"] = 30
        metadata[h]["xtit"] = "E_{T}^{miss,sign}"
    elif "lep_njet30_" in h:
        metadata[h]["xmin"] = 0
        metadata[h]["xmax"] = 10
        metadata[h]["xtit"] = "N_{jets}^{p_{T}>30}"
    elif "_mt2_" in h:
        metadata[h]["xmin"] = 0
        metadata[h]["xmax"] = 300
        metadata[h]["xtit"] = "m_{T2} [GeV]"
    elif "_mtw_" in h:
        metadata[h]["xmin"] = 0
        metadata[h]["xmax"] = 300
        metadata[h]["xtit"] = "m_{T}^{W} [GeV]"
    elif "_mll_" in h:
        metadata[h]["xmin"] = 25
        metadata[h]["xmax"] = 1000
        metadata[h]["xtit"] = "m_{ll} [GeV]"
    elif "_pT1_" in h:
        metadata[h]["xmin"] = 25
        metadata[h]["xmax"] = 999
        metadata[h]["xtit"] = "p_{T}^{1} [GeV]"
    elif "_pT2_" in h:
        metadata[h]["xmin"] = 0
        metadata[h]["xmax"] = 350
        metadata[h]["xtit"] = "p_{T}^{2} [GeV]"
    # elif "_nbj77_" in h:
    #     metadata[h]["xmin"] = 0
    #     metadata[h]["xmax"] = 4
    #     metadata[h]["xtit"] = "N(b-jets (85%))"
    elif "_lep_TCF_n" in h:
        metadata[h]["xmin"] = 0
        metadata[h]["xmax"] = 12
        metadata[h]["xtit"] = "Fake Category"

    if "eff" in h: 
        if "FAKE" in h: metadata[h]["ytit"] = "Fake Rate"
        if "REAL2L" in h: metadata[h]["ytit"] = "Real Efficiency"
    else: metadata[h]["ytit"] = "Entries"

    if "EE" in h: metadata[h]["label"] = "Electrons"
    elif "MM" in h: metadata[h]["label"] = "Muons"


def doSubtraction2D(h_nT_MC,h_nL_MC,h_nT_DATA,h_nL_DATA,datap,is2D):

    h_nT_DATA_subt = h_nT_DATA.Clone(h_nT_DATA.GetName()+"_subt")
    h_nL_DATA_subt = h_nL_DATA.Clone(h_nL_DATA.GetName()+"_subt")
    h_eff_DATA_subt = h_nL_DATA.Clone((h_nL_DATA.GetName()).replace("nL","eff")+"_subt")

    for i in range(1,h_nT_MC.GetNbinsX()+2):
        for j in range(1,h_nT_MC.GetNbinsY()+2):
            print("bf scaling h_nT_MC (%i,%i) = %f" %(i,j,h_nT_MC.GetBinContent(i,j)))
            print("bf scaling h_nL_MC (%i,%i) = %f" %(i,j,h_nL_MC.GetBinContent(i,j)))

    # if datap == "data17":
    #     print "Scaling"
    #     h_nT_MC.Scale(0.5)
    #     h_nL_MC.Scale(0.5)


    for i in range(1,h_nT_MC.GetNbinsX()+2):
        for j in range(1,h_nT_MC.GetNbinsY()+2):
            print("h_nT_MC (%i,%i) = %f" %(i,j,h_nT_MC.GetBinContent(i,j)))
            print("h_nL_MC (%i,%i) = %f" %(i,j,h_nL_MC.GetBinContent(i,j)))

            print("h_nT_DATA (%i,%i) = %f" %(i,j,h_nT_DATA.GetBinContent(i,j)))
            print("h_nL_DATA (%i,%i) = %f" %(i,j,h_nL_DATA.GetBinContent(i,j)))

            h_nT_DATA_subt.SetBinContent(i,j,h_nT_DATA.GetBinContent(i,j)-h_nT_MC.GetBinContent(i,j))
            h_nL_DATA_subt.SetBinContent(i,j,h_nL_DATA.GetBinContent(i,j)-h_nL_MC.GetBinContent(i,j))

            if ((h_nT_DATA.GetBinError(i,j)*h_nT_DATA.GetBinError(i,j))-(h_nT_MC.GetBinError(i,j)*h_nT_MC.GetBinError(i,j))) > 0:
                h_nT_DATA_subt.SetBinError(i,j,sqrt((h_nT_DATA.GetBinError(i,j)*h_nT_DATA.GetBinError(i,j))-(h_nT_MC.GetBinError(i,j)*h_nT_MC.GetBinError(i,j))))
            else:
                h_nT_DATA_subt.SetBinError(i,j,0)
            if ((h_nL_DATA.GetBinError(i,j)*h_nL_DATA.GetBinError(i,j))-(h_nL_MC.GetBinError(i,j)*h_nL_MC.GetBinError(i,j))) > 0:
                h_nL_DATA_subt.SetBinError(i,j,sqrt((h_nL_DATA.GetBinError(i,j)*h_nL_DATA.GetBinError(i,j))-(h_nL_MC.GetBinError(i,j)*h_nL_MC.GetBinError(i,j))))
            else:
                h_nL_DATA_subt.SetBinError(i,j,0)


    #h_nT_DATA_subt.Add(h_nT_MC,-1.0)
    #h_nL_DATA_subt.Add(h_nL_MC,-1.0)
    for i in range(1,h_nT_MC.GetNbinsX()+2):
        for j in range(1,h_nT_MC.GetNbinsY()+2):
            print("h_nT_DATA_subt (%i,%i) = %f" %(i,j, h_nT_DATA_subt.GetBinContent(i,j)))
            print("h_nL_DATA_subt (%i,%i) = %f" %(i,j, h_nL_DATA_subt.GetBinContent(i,j)))

    h_eff_DATA_subt.Divide(h_nT_DATA_subt,h_nL_DATA_subt,1.,1.,'b')

    if is2D:
        asym = TEfficiency(h_nT_DATA_subt,h_nL_DATA_subt)
        asym.SetStatisticOption(TEfficiency.kBBayesian)
        asym.SetConfidenceLevel(0.68);
        asym.SetBetaAlpha(1.)
        asym.SetBetaBeta(1.) 
        for i in range(0,h_eff_DATA_subt.GetNbinsX()+1):
            for j in range(0,h_eff_DATA_subt.GetNbinsY()+1):
                ibin = h_eff_DATA_subt.GetBin(i,j)
                h_eff_DATA_subt.SetBinError(ibin,(asym.GetEfficiencyErrorLow(ibin)+asym.GetEfficiencyErrorUp(ibin))/2.0)
    else:
        asym = TGraphAsymmErrors();
        print(":::->Dividing h_nT_DATA_subt = %.2f on h_nL_DATA_subt = %.2f" %(h_nT_DATA_subt.Integral(),h_nL_DATA_subt.Integral()))
        asym.Divide(h_nT_DATA_subt,h_nL_DATA_subt,"cl=0.683 b(1,1) mode")
        for i in range(0,h_eff_DATA_subt.GetNbinsX()+1):
            h_eff_DATA_subt.SetBinError(i+1,asym.GetErrorY(i))
            

    print(type(h_eff_DATA_subt))

    h_eff_DATA_subt.SetDirectory(0)

    return h_eff_DATA_subt
    
def doSubtraction(h_nT_MC,h_nL_MC,h_nT_DATA,h_nL_DATA,datap):

    h_nT_DATA_subt = h_nT_DATA.Clone(h_nT_DATA.GetName()+"_subt")
    h_nL_DATA_subt = h_nL_DATA.Clone(h_nL_DATA.GetName()+"_subt")
    print("Changing name to %s" %((h_nL_DATA.GetName()).replace("nL","eff")+"_subt"))
    h_eff_DATA_subt = h_nL_DATA.Clone((h_nL_DATA.GetName()).replace("nL","eff")+"_subt")

    for i in range(1,h_nT_MC.GetNbinsX()+2):
        print("bf scaling h_nT_MC (%i) = %f" %(i,h_nT_MC.GetBinContent(i)))
        print("bf scaling h_nL_MC (%i) = %f" %(i,h_nL_MC.GetBinContent(i)))

    # if datap == "data17":
    #     print "Scaling"
    #     h_nT_MC.Scale(0.5)
    #     h_nL_MC.Scale(0.5)

    prct_T = -1.0
    prev_prct_T = -1.0
    prct_L = -1.0
    prev_prct_L = -1.0
    for i in range(1,h_nT_MC.GetNbinsX()+1):
        print("h_nT_MC (%i) = %.2f +/- %.2f" %(i,h_nT_MC.GetBinContent(i),h_nT_MC.GetBinError(i)))
        print("h_nL_MC (%i) = %.2f +/- %.2f" %(i,h_nL_MC.GetBinContent(i),h_nL_MC.GetBinError(i)))

        print("h_nT_DATA (%i) = %f" %(i,h_nT_DATA.GetBinContent(i)))
        print("h_nL_DATA (%i) = %f" %(i,h_nL_DATA.GetBinContent(i)))

        if h_nT_DATA.GetBinContent(i)-h_nT_MC.GetBinContent(i) > h_nL_DATA.GetBinContent(i)-h_nL_MC.GetBinContent(i):
            print("Correct from nL = ",h_nL_MC.GetBinContent(i))
            print("Correct from nT = ",h_nT_MC.GetBinContent(i))
            h_nL_MC.SetBinContent(i,h_nL_MC.GetBinContent(i)-h_nL_MC.GetBinError(i))
            h_nT_MC.SetBinContent(i,h_nT_MC.GetBinContent(i)+h_nT_MC.GetBinError(i))
            print("to nL = ",h_nL_MC.GetBinContent(i))
            print("to nT = ",h_nT_MC.GetBinContent(i))

        if h_nT_DATA.GetBinContent(i): print("nT : subtract %.2f %% from data" %(100.0*(h_nT_MC.GetBinContent(i)/h_nT_DATA.GetBinContent(i))))
        if h_nL_DATA.GetBinContent(i): print("nL : subtract %.2f %% from data" %(100.0*(h_nL_MC.GetBinContent(i)/h_nL_DATA.GetBinContent(i))))

        if prct_T < 1:
            prev_prct_T = prct_T
            if h_nT_DATA.GetBinContent(i):
                prct_T      = (h_nT_MC.GetBinContent(i)/h_nT_DATA.GetBinContent(i))
        if prct_L < 1:
            prev_prct_L = prct_L
            if h_nL_DATA.GetBinContent(i):
                prct_L      = (h_nL_MC.GetBinContent(i)/h_nL_DATA.GetBinContent(i))

        if not h_nT_MC.GetBinContent(i) >= h_nT_DATA.GetBinContent(i):
            h_nT_DATA_subt.SetBinContent(i,h_nT_DATA.GetBinContent(i)-h_nT_MC.GetBinContent(i))
            if ((h_nT_DATA.GetBinError(i)*h_nT_DATA.GetBinError(i))-(h_nT_MC.GetBinError(i)*h_nT_MC.GetBinError(i))) >= 0:
                h_nT_DATA_subt.SetBinError(i,sqrt((h_nT_DATA.GetBinError(i)*h_nT_DATA.GetBinError(i))-(h_nT_MC.GetBinError(i)*h_nT_MC.GetBinError(i))))
            else:
                h_nT_DATA_subt.SetBinError(i,0)
        else:
            print("Tight: Scaling with %.2f" %(1.0-prev_prct_T))
            h_nT_DATA_subt.SetBinContent(i,h_nT_DATA_subt.GetBinContent(i)*(1.0-prev_prct_T))
        if not h_nL_MC.GetBinContent(i) >= h_nL_DATA.GetBinContent(i):
            h_nL_DATA_subt.SetBinContent(i,h_nL_DATA.GetBinContent(i)-h_nL_MC.GetBinContent(i))
            if ((h_nL_DATA.GetBinError(i)*h_nL_DATA.GetBinError(i))-(h_nL_MC.GetBinError(i)*h_nL_MC.GetBinError(i))) >= 0:
                h_nL_DATA_subt.SetBinError(i,sqrt((h_nL_DATA.GetBinError(i)*h_nL_DATA.GetBinError(i))-(h_nL_MC.GetBinError(i)*h_nL_MC.GetBinError(i))))
            else:
                h_nL_DATA_subt.SetBinError(i,0)
        else:
            print("Lose: Scaling with %.2f" %(1.0-prev_prct_L))
            h_nT_DATA_subt.SetBinContent(i,h_nT_DATA_subt.GetBinContent(i)*(1.0-prev_prct_L))

    h_eff_DATA_subt.Divide(h_nT_DATA_subt,h_nL_DATA_subt,1.,1.,'b')
    #h_nT_DATA_subt.Add(h_nT_MC,-1.0)
    #h_nL_DATA_subt.Add(h_nL_MC,-1.0)
    for i in range(1,h_nT_MC.GetNbinsX()+1):
        print("h_nT_DATA_subt (%i) = %f" %(i, h_nT_DATA_subt.GetBinContent(i)))
        print("h_nL_DATA_subt (%i) = %f" %(i, h_nL_DATA_subt.GetBinContent(i)))

    asym = TGraphAsymmErrors();
    print("<-:::Dividing h_nT_DATA_subt = %.2f on h_nL_DATA_subt = %.2f" %(h_nT_DATA_subt.Integral(),h_nL_DATA_subt.Integral()))
    for ibin in range(1,h_nT_DATA_subt.GetNbinsX()+1):
        print("Dividing : ibin %i : nT = %.2f, nL = %.2f" %(ibin,h_nT_DATA_subt.GetBinContent(ibin),h_nL_DATA_subt.GetBinContent(ibin)))
    asym.Divide(h_nT_DATA_subt,h_nL_DATA_subt,"cl=0.683 b(1,1) mode")
    ni = 0
    for i in range(0,h_eff_DATA_subt.GetNbinsX()+1):
        print("Bin %i" %i)
        print("Hist high/low edge: %.2f/%.2f" %(h_eff_DATA_subt.GetXaxis().GetBinLowEdge(i+1),h_eff_DATA_subt.GetXaxis().GetBinUpEdge(i+1)))
        print("Asym high/low edge: %.2f/%.2f" %(asym.GetXaxis().GetBinLowEdge(ni),asym.GetXaxis().GetBinUpEdge(ni)))
        if h_eff_DATA_subt.GetBinContent(i+1):
            print("asym err = ", asym.GetErrorY(ni))
            h_eff_DATA_subt.SetBinError(i+1,asym.GetErrorY(ni))
            ni += 1
        else:
            h_eff_DATA_subt.SetBinError(i+1,0.0)

    return h_eff_DATA_subt, h_nT_DATA_subt, h_nL_DATA_subt
    

def doRebin(hT,bins=array.array('f',[0,10,20,25,30,40,50,60,80,100,200])):
    print(("Histogram is %s"%hT.GetName())) 
    VB = 0
    name = ""
    name = hT.GetName()+"_REBIN"
    hT_reb = TH1F(name,name,len(bins)-1,bins)
    for nbin in range(1,hT.GetNbinsX()+2):
        if VB > 1: print(("Adding %.2f from value %.2f into bin %i" %(hT.GetBinContent(nbin),abs(hT.GetBinLowEdge(nbin)),hT_reb.FindBin(abs(hT.GetBinLowEdge(nbin))))))
        hT_reb.AddBinContent(hT_reb.FindBin(abs(hT.GetBinLowEdge(nbin))),hT.GetBinContent(nbin))
    for nbin in range(1,hT_reb.GetNbinsX()+2):
        hT_reb.SetBinError(nbin,sqrt(hT_reb.GetBinContent(nbin)) if hT_reb.GetBinContent(nbin) >= 0 else 0.0)
    return hT_reb

def getEff(hT,hL):
    print(("Histogram %s"%hT.GetName()))
    hEff = hT.Clone(hT.GetName().replace("nT","eff"))
    hEff.Clear();
    print(("nT = ",hT.Integral()))
    print(("nL = ",hL.Integral()))
    hT.ClearUnderflowAndOverflow()
    hL.ClearUnderflowAndOverflow()
    asym = TGraphAsymmErrors();
    print("Dividing hT = %.2f on hL = %.2f" %(hT.Integral(),hL.Integral()))
    for ibin in range(1,hT.GetNbinsX()+2):
        print("hT = ",hT.GetBinContent(ibin))
        print("hL = ",hL.GetBinContent(ibin))
        #if hT.GetBinContent(ibin) > hL.GetBinContent(ibin):
        #    print "----> OBS!! <-----"
        #    hT.SetBinContent(ibin,hT.GetBinContent(ibin-1))
        #    hL.SetBinContent(ibin,hL.GetBinContent(ibin-1))
        #    print "hT (fixed) = ",hT.GetBinContent(ibin)
        #    print "hL (fixed) = ",hL.GetBinContent(ibin)
    hEff.Divide(hT,hL,1.,1.,'b')
    asym.Divide(hT,hL,"cl=0.683 b(1,1) mode")
    print("after")
    ni = 0
    for i in range(0,hEff.GetNbinsX()+1):
        if hT.GetBinContent(i+1) or hL.GetBinContent(i+1):
            print("Setting bin %i: %.2f +/- %.2f" %(i+1,hEff.GetBinContent(i+1),asym.GetErrorY(ni)))
            hEff.SetBinError(i+1,asym.GetErrorY(ni))
            ni += 1
        else:
            hEff.SetBinError(i+1,0)
    hEff.SetDirectory(0)
    hT.SetDirectory(0)
    hL.SetDirectory(0)
    print(type(hEff))
    return hEff

def getSortedList(nev):
    bkg_sorted_nev = []
    while len(bkg_sorted_nev) != len(list(nev.keys())):
        maxi = -999
        for key in list(nev.keys()):
            if 'data' in key: continue
            if key in bkg_sorted_nev: continue
            if nev[key] > maxi:
                maxi = nev[key]
                maxi_key = key
        bkg_sorted_nev.append(maxi_key)
    return bkg_sorted_nev

def getAbs(totals,errors,binlabel,hname):
    h_fract = {}
    for i in range(0,len(totals)):
        for lab in list(totals[i].keys()):
            if lab == "ALL" or lab == "DATA": continue
            if not lab in list(h_fract.keys()):
                h_fract[lab] = TH1F("h_fract_%s_%s"%(hname,lab),"h_fract_%s_%s"%(lab,hname),len(totals),0,len(totals))
                if 'm_type' in list(d_samp[lab].keys()): h_fract[lab].SetMarkerStyle(d_samp[lab]['m_type'])
                if 'l_color' in list(d_samp[lab].keys()):h_fract[lab].SetLineColor(d_samp[lab]['l_color'])
                if 'm_color' in list(d_samp[lab].keys()): h_fract[lab].SetMarkerColor(d_samp[lab]['m_color'])
                if 'f_color' in list(d_samp[lab].keys()): 
                    print("Type is ", lab)
                    print("Setting color ",d_samp[lab]['f_color'])
                    h_fract[lab].SetFillColor(d_samp[lab]['f_color'])
                    h_fract[lab].SetLineColor(d_samp[lab]['f_color'])
            if i < len(binlabel):
                print("i", i+1)
                print("binl", binlabel[i])
                h_fract[lab].GetXaxis().SetBinLabel(i+1,binlabel[i])
            if totals[i]["ALL"] and totals[i][lab]:
                if totals[i][lab] < 0:
                    totals[i]["ALL"] -= totals[i][lab]
                    continue
                h_fract[lab].SetBinContent(i+1,totals[i][lab])
                err = sqrt(totals[i][lab])#/totals[i]["ALL"] * sqrt((errors[i][lab]/(totals[i][lab]*totals[i][lab])) + (errors[i]["ALL"]/(totals[i]["ALL"]*totals[i]["ALL"])))
                h_fract[lab].SetBinError(i+1,err)

    return h_fract


def getFract(totals,errors,binlabel,hname):
    h_fract = {}
    for i in range(0,len(totals)):
        for lab in list(totals[i].keys()):
            if lab == "ALL" or lab == "DATA": continue
            if not lab in list(h_fract.keys()):
                h_fract[lab] = TH1F("h_fract_%s_%s"%(hname,lab),"h_fract_%s_%s"%(lab,hname),len(totals),0,len(totals))
                if 'm_type' in list(d_samp[lab].keys()): h_fract[lab].SetMarkerStyle(d_samp[lab]['m_type'])
                if 'l_color' in list(d_samp[lab].keys()):h_fract[lab].SetLineColor(d_samp[lab]['l_color'])
                if 'm_color' in list(d_samp[lab].keys()): h_fract[lab].SetMarkerColor(d_samp[lab]['m_color'])
                if 'f_color' in list(d_samp[lab].keys()): 
                    print("Type is ", lab)
                    print("Setting color ",d_samp[lab]['f_color'])
                    h_fract[lab].SetFillColor(d_samp[lab]['f_color'])
                    h_fract[lab].SetLineColor(d_samp[lab]['f_color'])
            if i < len(binlabel):
                print("i", i+1)
                print("binl", binlabel[i])
                h_fract[lab].GetXaxis().SetBinLabel(i+1,binlabel[i])
            if totals[i]["ALL"] and totals[i][lab]:
                if totals[i][lab] < 0:
                    totals[i]["ALL"] -= totals[i][lab]
                    continue
                h_fract[lab].SetBinContent(i+1,totals[i][lab]/totals[i]["ALL"])
                err = totals[i][lab]/totals[i]["ALL"] * sqrt((errors[i][lab]/(totals[i][lab]*totals[i][lab])) + (errors[i]["ALL"]/(totals[i]["ALL"]*totals[i]["ALL"])))
                h_fract[lab].SetBinError(i+1,err)

    return h_fract



hist = {}
hist_xsecup = {}
hist_xsecdw = {}
hEff = {}
hEff_xsecup = {}
hEff_xsecdw = {}
evcount = {}

rebarray = array.array('f',[])

set_palette()
is2D = False

addname = ""
if datakey == "data18": addname = "mc16e"
if datakey == "data17": addname = "mc16d"
if datakey == "data15-16": addname = "mc16a"
if datakey == "alldata": addname = "mc16ade"

midl_histname = histname[:]
for h in midl_histname:
    if "trueREAL" in h:# and not estfake:
        histname.append(h.replace("trueREAL","trueFAKE"))
        
        
for h in histname:
    print(h)
    doRB = True
    doRBSimple = False
    if "TCF" in h: 
        doRB = False
    elif "mu_pT" in h:
        print("-------------------------------------------------------->mu pT")
        doRB = False
        reby = array.array('f',[20,30,40,50,1000])
        if datakey == "data15-16":
            rebx = array.array('f',[0,20,25,30,35,55])
        else:
            rebx = array.array('f',[0,20,25,30,35,40,45,50,55,75,100])
            # for i in range(0,100):
        #     rebx.append(i)
    elif "pT_eta" in h:
        print("-------------------------------------------------------->pt_eta")
        rebx = array.array('f',[0,15,25,30,35,40,50,60,70,80,90,100,120,140,200])
        #reby = array.array('f',[0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.7])
        #reby = array.array('f',[0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.37,1.52,1.6,1.8,2.0,2.2,2.3,2.4,2.5,2.7])
        reby = array.array('f',[0.0,0.16,0.32,0.48,0.64,0.8,0.96,1.12,1.28,1.44,1.5,1.66,1.82,1.98,2.14,2.3,2.78])
        if "REAL2L02" in h or "REAL2L08" in h:
            if "_E_" in h: 
                rebx = array.array('f',[0,15,25,30,35,40,50,60,70,80,90,100,120,140,200,1000])
                reby = array.array('f',[0.0,0.4,0.6,0.8,1.0,1.2,1.37,1.52,1.8,2.0,2.7])#[0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.7])
            elif "_M_" in h:
                rebx = array.array('f',[0,15,25,35,40,45,50,60,70,80,90,100,120,140,200,1000])
                reby = array.array('f',[0.0,0.4,0.8,1.4,1.7,2.0,2.7])#[0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.7])
        elif "FAKE2L02" in h or "FAKE2L01" in h and not "trueLight" in h:
            rebx = array.array('f',[25,30,35,40,60,80,200])
            reby = array.array('f',[0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.7])
        #elif "2L04" in h and "trueLight" in h:
        #    rebx = array.array('f',[25,35,50,70,1000])
        #    reby = array.array('f',[0.0,0.8,1.6,2.6])
        elif "FAKE2L20" in h:
            rebx = array.array('f',[25,26,27,28,29,30,32,35,40,45,50,60,80,200])
            reby = array.array('f',[0.0,0.5,1.5,2.0,2.7])
        elif "FAKE2L23" in h or "FAKE2L24" in h or "FAKE2L25" in h: 
            rebx = array.array('f',[25,40,1000])
            reby = array.array('f',[0.0,1.2,2.6])
    elif ("pT_origin" in h or "pT_type" in h) and not "type_origin" in h:
        print("-------------------------------------------------------->pt origin")
        rebx = array.array('f',[25,30,40,50,60,80,200])
        if "pT_origin" in h: reby = array.array('f',[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45])
        if "pT_type" in h: reby = array.array('f',[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40])
    elif "type_origin" in h:
        print("-------------------------------------------------------->type vs origin")
        doRB = False
        doRBSimple = False
    elif ("pT" in h or "pt" in h) and ("nT" in h or "nL" in h) and not estfake: 
        print("-------------------------------------------------------->pt no estfake")
        if "2L23" in h or "2L24" in h or "2L25" in h or "2L26" in h or "2L27" in h or "truePromptPhotonConversion_FAKE2L05" in h:
            if "lepnotrig" in h:
                #if datakey == "data15-16":
                rebarray = array.array('f',[5,1000])#[0,10,20,30,40,50,100])
                #else:
                #    rebarray = array.array('f',[5,10,15,20,25,30,40,50,60,80,140,200,1000])
            else:
                rebarray = array.array('f',[5,10,15,20,25,30,40,50,60,80,140,200,1000])
            doLOG = True
        #elif ("2L04" in h or "2L02" in h) and "trueLight" in h:
        #    rebarray = array.array('f',[25,35,50,80,1000])
        elif "_M_" in h and ("FAKE2L20" in h or "FAKE2L01" in h or "FAKE2L21" in h or estfake):
            if "lepnotrig" in h:
                rebarray = array.array('f',[5,10,15,20,25,30,40,50,60,70,100,1000])#[0,10,20,30,40,50,100])
            elif "mu60_0eta105_msonly" in h or "mu50" in h:
                rebarray = array.array('f',[5,10,15,20,25,30,35,40,45,50,55,60,65,70,80,100,1000])#[10,25,30,40,50,60,80,1000])
                doLOG = True
            elif "mu20_iloose" in h:
                rebarray = array.array('f',[5,20,40,60,100,1000])#[10,25,30,40,50,60,80,1000])
            else:
                #rebarray = array.array('f',[10,25,30,35,40,50,1000])#[25,30,50,100,1000])
                rebarray = array.array('f',[5,10,15,20,25,30,35,40,45,50,60,1000])#[10,25,30,40,50,60,80,1000])
                doLOG = True
        elif ("EE" in h or "ALL" in h or "EM" in h) and ("FAKE2L04" in h or "FAKE2L01" in h or "FAKE2L03" in h or "FAKE2L05" in h or "FAKE2L02" in h):
            if "lepnotrig" in h:
                if datakey == "data15-16":
                    rebarray = array.array('f',[5,10,15,20,25,35,50,80,120,1000])#[0,10,20,30,40,50,100])
                else:
                    rebarray = array.array('f',[5,10,15,20,25,35,50,60,80,110,140,1000])
            else:
                rebarray = array.array('f',[5,10,15,20,25,35,50,60,70,80,90,110,140,1000])
            doLOG = True
        elif "_E_" in h and ("FAKE2L02" in h or "FAKE2L01" in h or "FAKE2L21" in h or "FAKE2L20" in h or "FAKE2L4" in h):
            if "lepnotrig" in h:
                #if datakey == "data15-16":
                rebarray = array.array('f',[5,10,15,20,30,1000])#[0,10,20,30,40,50,100])
                #else:
                #    rebarray = array.array('f',[5,10,15,20,25,30,35,40,45,50,60,70,80,100,200,1000])#[0,10,20,30,40,50,100])
            elif "e60_lhmedium_nod0" in h:
                rebarray = array.array('f',[5,10,15,20,25,30,35,40,45,50,60,70,80,100,200,1000])#[0,10,20,30,40,50,100])
            elif "e300_etcut" in h:
                rebarray = array.array('f',[5,10,15,20,25,30,35,40,45,50,60,70,80,100,200,300,1000])#[0,10,20,30,40,50,100])
            elif "e140_lhloose_nod0" in h or "e120_lhloose" in h:
                rebarray = array.array('f',[5,10,15,20,25,30,35,40,45,50,60,70,80,100,120,140,200,300,1000])#[0,10,20,30,40,50,100])
            elif "e24_lhmedium" in h:
                rebarray = array.array('f',[5,20,25,40,60,100,1000])#[10,25,30,40,50,60,80,1000])
            else:
                rebarray = array.array('f',[5,10,15,20,25,30,35,40,45,50,60,70,80,100,200,1000])#[0,10,20,30,40,50,100])
                doLOG = True
        else:
            print("Here I am!")
            if "lepnotrig" in h or "eventTrig" in h:
                if "REAL2L03" in h or "REAL2L03" in h:
                    rebarray = array.array('f',[5,10,15,30,1000])
                else:
                    rebarray = array.array('f',[5,10,15,30,1000])
            else:
                rebarray = array.array('f',[5,10,15,20,25,35,40,50,60,70,80,120,1000])

    elif "lep_pT1" in h:
        rebarray = array.array('f',[25,30,35,40,45,50,60,70,80,90,100,200,1000])
        #doRB = False
        #doRBSimple = True
        #rebfac = 4
        print("-------------------------------------------------------->pt1!")
        #elif "pT" in h and estfake: 
        #    rebarray = array.array('f',[25,30,40,50,60,100,1000])
    elif "lep_eta" in h:
        print("-------------------------------------------------------->eta")
        #rebarray = array.array('f',[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0])
        rebarray = array.array('f',[0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.37,1.52,1.6,1.8,2.0,2.2,2.3,2.4,2.5,2.7])
        #rebarray = array.array('f',[0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.7])
        if "2L23" in h or "2L24" in h or "2L25" in h or "2L26" in h or "2L27" in h:
            rebarray = array.array('f',[0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,2.0,2.2,2.7])
    elif "lep_eta" in h:
        print("-------------------------------------------------------->eta")
        doRB = False
    elif "lep_mll" in h:
        print("-------------------------------------------------------->mll")
        doRB = False
        doRBSimple = True
        rebfac = 10
    elif "lep_metsig" in h:
        print("-------------------------------------------------------->mll")
        doRB = False
        doRBSimple = False
        rebfac = 2
    elif "lep_pT_" in h:
        doRB = False
        doRBSimple = True
        rebfac = 30
        print("-------------------------------------------------------->here3")
    elif "HFreg" in h :
        print("-------------------------------------------------------->HFreg")
        doRB = False
        doRBSimple = False
        rebfac = 25
    elif "lep_nbj" in h or "lep_njet" in h:
        print("-------------------------------------------------------->nbj")
        doRB = False
        doRBSimple = False
        rebfac = 1
    elif "lep_pT1" in h and estfake:
        #rebarray = array.array('f',[25,30,35,40,45,50,60,100,1000])
        doRB = False
        doRBSimple = True
        rebfac = 4
        print("-------------------------------------------------------->pt1!")
    elif "lep_pT2" in h:# and estfake:
        #rebarray = array.array('f',[25,30,35,40,45,50,60,100,1000])
        doRB = False
        doRBSimple = True
        rebfac = 4
        print("-------------------------------------------------------->here!")
    elif ("lep_met_" in h and estfake) or ("lep_MET" in h):
        #rebarray = array.array('f',[25,30,35,40,45,50,60,100,1000])
        doRB = False
        doRBSimple = True
        rebfac = 4
        print("-------------------------------------------------------->here!")
    elif "lep_mt2" in h:# and estfake:
        rebarray = array.array('f',[0,60,80,100,140,400])
        doRB = False
        doRBSimple = True
        rebfac = 2
        print("-------------------------------------------------------->here!")
    elif "lep_njet30" in h and estfake:
        #rebarray = array.array('f',[0,60,80,100,140,400])
        doRB = False
        doRBSimple = False
        rebfac = 1
        print("-------------------------------------------------------->here!")
    elif "lep_mll" in h and estfake:
        #rebarray = array.array('f',[25,30,35,40,45,50,60,100,1000])
        doRB = False
        doRBSimple = True
        rebfac = 10
        print("-------------------------------------------------------->here!")
    elif "lep_mt2" in h and "CRTOP" in h:
        rebarray = array.array('f',[100,120,140,160,180,200])
        doRB = True
        doRBSimple = False
        rebfac = 20
        print("-------------------------------------------------------->here!")
    else:
        doRB = False
        doRBSimple = False
        rebfac = 1
        print("-------------------------------------------------------->here2!")
    print(h)
    if not h in hist:
        hist[h] = {}
        hist_xsecup[h] = {}
        hist_xsecdw[h] = {}
        evcount[h] = {}
    for typ in list(d_samp.keys()):
        if typ in ["CF","fake","LF","HF","CO","mc16a","mc16c","mc16d","mc16e","Unknown","IsoElectron","ChargeFlipIsoElectron","PromptMuon","PromptPhotonConversion",
                   "ElectronFromMuon","TauDecay","BHadronDecay","CHadronDecay","LightFlavorDecay","all"]: continue
        #if typ not in ['PHZZ','PHWZ','PHWW']: continue
        #if typ not in ['diboson']: continue
        if "TCF" in h and "data" in typ: continue
        if "type_origin" in h and "data" in typ: continue
        if not doMC and not "data" in typ: continue
        if "data" in typ and not datakey and not (alldata or comb1516): continue
        #if "data18" in typ: continue
        #if not "boson" in typ: continue
        #if estfake and "EE" in h and typ in ["Zjets","ttbar"]: continue
        #if estfake and not typ in ["Wjets","singleTop","Vgamma","higgs","ttbar",datakey]: continue
        if alldata: dataper = ["data18","data17","data1516"]
        elif comb1516: dataper = ["data1516"]
        elif datakey in ["data18","data17"]: dataper = [datakey]
        elif datakey == "data15-16": dataper = ["data1516"]
        else: dataper = []
        oldtyp = typ
        print("1) oldtyp = %s, typ = %s" %(oldtyp,typ))
        for fn in dataper:

            typ = oldtyp

            print("2) oldtyp = %s, typ = %s" %(oldtyp,typ))


            if fn in ["data18","data17"]:
                lumikey = fn
            elif fn in ["data1516"]:
                lumikey = "data15-16"

            if typ in ["data18","data17"] and not fn == typ: continue
            if typ in ["data15-16","data15","data16"] and not fn == "data1516": continue

            #print "-->",inputdir.split("_")[-1]
            postfix = ""
            if "test" in inputdir or "oldORDef" in inputdir or "EST" in inputdir or "RATES" in inputdir or "R21" in inputdir or "WGT" in inputdir or "OLD2L2J" in inputdir or "ZMET" in inputdir:
                postfix = inputdir.split("_")[-1].replace("/","")
                inputdir = inputdir.replace(inputdir.split("_")[-2],fn.replace("data","year"))
            else:
                inputdir = inputdir.replace(inputdir.split("_")[-1],fn.replace("data","year")+"/")

            print(typ)
            #fnames = glob.glob("%s/%s_merged_processed_HISTOGRAMS_0_0.root" %(inputdir,typ))
            fnames = glob.glob("%s/%s_*merged_processed_HISTOGRAMS_*.root" %(inputdir,typ))
            print(fnames)
            if len(fnames) <= 0:
                print(("WARNING \t Could not find file for %s in %s" %(typ,inputdir)))
                continue
            else:
                print("DEBUG FOund file for %s"%typ)
            #fname = fnames[0]
            #print "fname is %s" %fname
            for fname in fnames:
                if not os.path.isfile(fname):
                    print(("WARNING \t Could not find file for %s in %s" %(typ,fn)))
                    if (datakey == "data18" or fn == "data18") and typ in ["lowMassDY","Vgamma"]:
                        fname = "%s/%s_merged_processed_HISTOGRAMS.root" %(inputdir.replace("data18","data17"),typ)
                        print(("INFO \t Using file from mc16d for %s in stead %s" %(typ,fname)))
                    else:
                        continue
                    if not os.path.isfile(fname):
                        print(("WARNING \t Could not find file for %s in %s" %(typ,fn)))
                        continue
                if "data" in typ:
                    oldtyp = typ
                    if alldata: typ = "alldata"
                    elif comb1516: typ = "data15-16"
                    print("3) oldtyp = %s, typ = %s" %(oldtyp,typ))
                    #lumi["alldata"] += lumi[lumikey]
                    if "true" in h or "ElectronFromMuon" in h:
                        print("-----> before ", h)
                        #htoget = getHistWithoutSyst(h)
                        htoget = ""
                        for t in h.split("_"):
                            if not "true" in t and not "ElectronFromMuon" in t: htoget += "%s_" %t
                        htoget = htoget[:-1]
                        print(htoget)
                        htoget = getHistWithoutSyst(htoget)
                        print("-----> after ", htoget)
                    else: 
                        #print "-----> before ", htoget
                        htoget = getHistWithoutSyst(h)
                        print("-----> after ", htoget)
                else:
                    if estfake and not "trueREAL" in h and not "trueFAKE" in h:
                        htoget = h+"_trueREAL"
                    else:
                        htoget = h
                print(htoget)            
                #print(typ)
                file = TFile(fname)
                rebisdone = False
                if not typ in list(hist[h].keys()):
                    print("INFO \t Creating %s %s from %s" %(h,typ,fn))
                    if type(file.Get(htoget)) is TH2F:
                        is2D = True
                        hist[h][typ] = TH2F()
                        hist_xsecup[h][typ] = TH2F()
                        hist_xsecdw[h][typ] = TH2F()
                    else:
                        hist[h][typ] = TH1F()
                        hist_xsecup[h][typ] = TH1F()
                        hist_xsecdw[h][typ] = TH1F()
                    try:
                        print("Getting ",htoget)
                        hist[h][typ] = (file.Get(htoget)).Clone(h+"_"+typ)
                    except:
                        print("Could not find histogram %s for type %s" %(htoget,typ))
                        sys.exit()

                    hist_xsecup[h][typ] = (file.Get(htoget)).Clone(h+"_"+typ+"_xsecup")
                    hist_xsecdw[h][typ] = (file.Get(htoget)).Clone(h+"_"+typ+"_xsecdw")
                    print("INFO \t Starting with %.2f events from histogram %s" %(hist[h][typ].Integral(),htoget))
                else:
                    h_temp = ((file.Get(htoget)).Clone(h+"_"+typ))
                    h_temp.SetDirectory(0)
                    print("INFO \t Adding %s %s from %s with %.2f events" %(h,typ,fn,h_temp.Integral()))
                    if doRB: 
                        if is2D: 
                            hist[h][typ] = doRebin2D(h_temp,rebx,reby)
                            hist_xsecup[h][typ] = doRebin2D(h_temp,rebx,reby)
                            hist_xsecdw[h][typ] = doRebin2D(h_temp,rebx,reby)
                        else:
                            #print "Adding %.2f events for histogram %s from %s" %(h_temp.Integral(0,h_temp.GetNbinsX()+1),htoget,typ)
                            hist[h][typ].Add(doRebin(h_temp,rebarray))
                            hist_xsecup[h][typ].Add(doRebin(h_temp,rebarray))
                            hist_xsecdw[h][typ].Add(doRebin(h_temp,rebarray))
                        rebisdone = True
                    elif doRBSimple:
                        h_temp = h_temp.Rebin(rebfac)
                        hist[h][typ].Add(h_temp)
                        hist_xsecup[h][typ].Add(h_temp)
                        hist_xsecdw[h][typ].Add(h_temp)
                        rebisdone = True
                    else:
                        hist[h][typ].Add(h_temp)
                        hist_xsecup[h][typ].Add(h_temp)
                        hist_xsecdw[h][typ].Add(h_temp)
                    if "data" in typ:
                        for i in range(0,hist[h][typ].GetNbinsX()+1):
                            print("BinContent(%i) = %.2f" %(i,hist[h][typ].GetBinContent(i)))
                    if not is2D: print("INFO \t Now with %.2f events" %hist[h][typ].Integral(0,hist[h][typ].GetNbinsX()+1))

                print("a")
                print("rebisdone = ",rebisdone)
                # old_htoget = htoget
                # if "trueREAL" in h and not "data" in typ:
                #     print "Before Adding is ", hist[h][typ].Integral()
                #     for trus in ["trueHF","trueLight"]:
                #         htoget = h.replace("trueREAL",trus)
                #         print "Adding", htoget
                #         hist[h][typ].Add((file.Get(htoget)).Clone(h+"_"+typ))
                #         print "Now is ", hist[h][typ].Integral()
                # htoget = old_htoget
                if type(hist[h][typ]) is TH2F: 
                    is2D = True
                    if not "TCF" in h and not "type_origin" in h:
                        doRB = True
                if not rebisdone:
                    if doRB: 
                        if is2D: 
                            hist[h][typ] = doRebin2D(hist[h][typ],rebx,reby)
                            hist_xsecup[h][typ] = doRebin2D(hist_xsecup[h][typ],rebx,reby)
                            hist_xsecdw[h][typ] = doRebin2D(hist_xsecdw[h][typ],rebx,reby)
                        else: 
                            print("Rebinning!")
                            hist[h][typ] = doRebin(hist[h][typ],rebarray)
                            hist_xsecup[h][typ] = doRebin(hist_xsecup[h][typ],rebarray)
                            hist_xsecdw[h][typ] = doRebin(hist_xsecdw[h][typ],rebarray)
                    elif doRBSimple:
                        print("rebfac = ", rebfac)
                        print("bin width = ", hist[h][typ].GetBinWidth(1))
                        hist[h][typ] = hist[h][typ].Rebin(rebfac)
                        hist_xsecup[h][typ] = hist_xsecup[h][typ].Rebin(rebfac)
                        hist_xsecdw[h][typ] = hist_xsecdw[h][typ].Rebin(rebfac)
                print("b")
                hist[h][typ].SetDirectory(0)
                hist_xsecup[h][typ].SetDirectory(0)
                hist_xsecdw[h][typ].SetDirectory(0)

                if not is2D: evcount[h][typ] = (hist[h][typ].Integral(0,hist[h][typ].GetNbinsX()+1))
                if "nT" in h:# or (estfake and "data" in typ):
                    rebisdone = False
                    if "nT" in h:
                        hnameL = h.replace("nT","nL")
                        htoget = htoget.replace("nT","nL")
                        if not hnameL in list(hist.keys()):
                            hist[hnameL] = {}
                            hist_xsecup[hnameL] = {}
                            hist_xsecdw[hnameL] = {}
                            evcount[hnameL] = {}
                    if not typ in list(hist[hnameL].keys()):
                        print("INFO \t Creating %s %s from %s" %(hnameL,typ,fn))
                        if type(file.Get(htoget)) is TH2F:
                            is2D = True
                            hist[hnameL][typ] = TH2F()
                            hist_xsecup[hnameL][typ] = TH2F()
                            hist_xsecdw[hnameL][typ] = TH2F()
                        else:
                            hist[hnameL][typ] = TH1F()
                            hist_xsecup[hnameL][typ] = TH1F()
                            hist_xsecdw[hnameL][typ] = TH1F()
                        try:
                            hist[hnameL][typ] = (file.Get(htoget)).Clone(hnameL+"_"+typ)
                        except:
                            print("Could not find histogram %s for type %s" %(htoget,typ))
                            sys.exit()

                        hist_xsecup[hnameL][typ] = (file.Get(htoget)).Clone(hnameL+"_"+typ+"_xsecup")
                        hist_xsecdw[hnameL][typ] = (file.Get(htoget)).Clone(hnameL+"_"+typ+"_xsecdw")
                        print("INFO \t Starting with %.2f events" %hist[hnameL][typ].Integral())
                    else:
                        h_temp = ((file.Get(htoget)).Clone(h+"_"+typ))
                        h_temp.SetDirectory(0)
                        #print "INFO \t Adding %s %s from %s with %.2f events" %(h,typ,fn,h_temp.Integral())
                        if doRB: 
                            if is2D: 
                                hist[hnameL][typ] = doRebin2D(h_temp,rebx,reby)
                                hist_xsecup[hnameL][typ] = doRebin2D(h_temp,rebx,reby)
                                hist_xsecdw[hnameL][typ] = doRebin2D(h_temp,rebx,reby)
                            else:
                                print("Adding %.2f events" %h_temp.Integral(0,h_temp.GetNbinsX()+1))
                                hist[hnameL][typ].Add(doRebin(h_temp,rebarray))
                                hist_xsecup[hnameL][typ].Add(doRebin(h_temp,rebarray))
                                hist_xsecdw[hnameL][typ].Add(doRebin(h_temp,rebarray))
                            rebisdone = True
                        elif doRBSimple:
                            h_temp = h_temp.Rebin(rebfac)
                            hist[hnameL][typ].Add(h_temp)
                            hist_xsecup[hnameL][typ].Add(h_temp)
                            hist_xsecdw[hnameL][typ].Add(h_temp)
                            rebisdone = True
                        else:
                            hist[hnameL][typ].Add(h_temp)
                            hist_xsecup[hnameL][typ].Add(h_temp)
                            hist_xsecdw[hnameL][typ].Add(h_temp)
                        if "data" in typ:
                            for i in range(0,hist[hnameL][typ].GetNbinsX()+1):
                                print("BinContent(%i) = %.2f" %(i,hist[hnameL][typ].GetBinContent(i)))
                        if not is2D: 
                            print("INFO \t Now with %.2f events" %hist[hnameL][typ].Integral(0,hist[hnameL][typ].GetNbinsX()+1))

                    print("a")
                    print("rebisdone = ",rebisdone)
                    # old_htoget = htoget
                    # if "trueREAL" in h and not "data" in typ:
                    #     print "Before Adding is ", hist[hnameL][typ].Integral()
                    #     for trus in ["trueHF","trueLight"]:
                    #         htoget = h.replace("trueREAL",trus)
                    #         print "Adding", htoget
                    #         hist[hnameL][typ].Add((file.Get(htoget)).Clone(h+"_"+typ))
                    #         print "Now is ", hist[hnameL][typ].Integral()
                    # htoget = old_htoget
                    if type(hist[hnameL][typ]) is TH2F: 
                        is2D = True
                        if not "TCF" in h and not "type_origin" in h:
                            doRB = True
                    if not rebisdone:
                        if doRB: 
                            if is2D: 
                                hist[hnameL][typ] = doRebin2D(hist[hnameL][typ],rebx,reby)
                                hist_xsecup[hnameL][typ] = doRebin2D(hist_xsecup[hnameL][typ],rebx,reby)
                                hist_xsecdw[hnameL][typ] = doRebin2D(hist_xsecdw[hnameL][typ],rebx,reby)
                            else: 
                                print("Rebinning!")
                                hist[hnameL][typ] = doRebin(hist[hnameL][typ],rebarray)
                                hist_xsecup[hnameL][typ] = doRebin(hist_xsecup[hnameL][typ],rebarray)
                                hist_xsecdw[hnameL][typ] = doRebin(hist_xsecdw[hnameL][typ],rebarray)
                        elif doRBSimple:
                            print("rebfac = ", rebfac)
                            print("bin width = ", hist[hnameL][typ].GetBinWidth(1))
                            hist[hnameL][typ] = hist[hnameL][typ].Rebin(rebfac)
                            hist_xsecup[hnameL][typ] = hist_xsecup[hnameL][typ].Rebin(rebfac)
                            hist_xsecdw[hnameL][typ] = hist_xsecdw[hnameL][typ].Rebin(rebfac)
                    print("b")
                    hist[hnameL][typ].SetDirectory(0)
                    hist_xsecup[hnameL][typ].SetDirectory(0)
                    hist_xsecdw[hnameL][typ].SetDirectory(0)

                    if not is2D: evcount[hnameL][typ] = (hist[hnameL][typ].Integral(0,hist[hnameL][typ].GetNbinsX()+1))

                    # if "nT" in h:
                    #     hnameL = h.replace("nT","nL")
                    #     htoget = htoget.replace("nT","nL")
                    # elif estfake:
                    #     hnameL = h+"_estfake"
                    #     htoget = htoget+"_estfake"
                    #     print "htoget",htoget
                    # if not hnameL in hist.keys():
                    #     hist[hnameL] = {}
                    #     hist_xsecup[hnameL] = {}
                    #     hist_xsecdw[hnameL] = {}
                    #     evcount[hnameL] = {}

                    # if not typ in hist[hnameL].keys():
                    #     print "INFO \t Creating %s %s from %s" %(h,typ,fn)

                    #     if is2D:
                    #         hist[hnameL][typ] = TH2F()
                    #     else:
                    #         hist[hnameL][typ] = TH1F()

                    #     hist[hnameL][typ] = (file.Get(htoget)).Clone(hnameL+"_"+typ)
                    #     if is2D:
                    #         hist_xsecup[hnameL][typ] = TH2F()
                    #     else:
                    #         hist_xsecup[hnameL][typ] = TH1F()
                    #     hist_xsecup[hnameL][typ] = (file.Get(htoget)).Clone(hnameL+"_"+typ+"_xsecup")
                    #     if is2D:
                    #         hist_xsecdw[hnameL][typ] = TH1F()
                    #     else:
                    #         hist_xsecdw[hnameL][typ] = TH2F()
                    #     hist_xsecdw[hnameL][typ] = (file.Get(htoget)).Clone(hnameL+"_"+typ+"_xsecdw")
                    # if doRB: 
                    #      if is2D: 
                    #          hist[hnameL][typ] = doRebin2D(hist[hnameL][typ],rebx,reby)
                    #          hist_xsecup[hnameL][typ] = doRebin2D(hist_xsecup[hnameL][typ],rebx,reby)
                    #          hist_xsecdw[hnameL][typ] = doRebin2D(hist_xsecdw[hnameL][typ],rebx,reby)
                    #      else: 
                    #          hist[hnameL][typ] = doRebin(hist[hnameL][typ],rebarray)
                    #          hist_xsecup[hnameL][typ] = doRebin(hist_xsecup[hnameL][typ],rebarray)
                    #          hist_xsecdw[hnameL][typ] = doRebin(hist_xsecdw[hnameL][typ],rebarray)
                    # elif doRBSimple:
                    #     print "rebfac = ", rebfac
                    #     print "bin width = ", hist[hnameL][typ].GetBinWidth(1);
                    #     hist[hnameL][typ] = hist[hnameL][typ].Rebin(rebfac)
                    #     hist_xsecup[hnameL][typ] = hist_xsecup[hnameL][typ].Rebin(rebfac)
                    #     hist_xsecdw[hnameL][typ] = hist_xsecdw[hnameL][typ].Rebin(rebfac)
                    # hist[hnameL][typ].SetDirectory(0)
                    # hist_xsecup[hnameL][typ].SetDirectory(0)
                    # hist_xsecdw[hnameL][typ].SetDirectory(0)
                    # if not is2D: evcount[hnameL][typ] = (hist[hnameL][typ].Integral(0,hist[hnameL][typ].GetNbinsX()+1))
                    # else: evcount[hnameL][typ] = (hist[hnameL][typ].Integral(0,hist[hnameL][typ].GetNbinsX()+1,0,hist[hnameL][typ].GetNbinsY()+1))
                    # #if "REAL" in h or "FAKE" in h:
                    if "nT" in h and not (alldata or comb1516) and (("_REAL" in h and not "trueFAKE" in h) or ("_FAKE" in h and not "trueREAL" in h) or not "true" in h):
                        #if "_FAKE" in h and "trueREAL" in h: continue
                        #if "trueREAL" in h and not "trueREAL" in hist[h][typ].GetName(): continue

                        print("Stored values", list(hEff.keys()))

                        print("hEff = ",hist[h][typ].GetName().replace("nT","eff"))

                        print("typ", typ)
                        print("Passed", hist[h][typ].GetBinContent(1))
                        print("Total", hist[h.replace("nT","nL")][typ].GetBinContent(1))

                        if is2D:
                            print("--> Dividing %s with %s" %(hist[h][typ].GetName(),hist[h.replace("nT","nL")][typ].GetName()))
                            # Avoid having more T than L (needed in some cases, mostly for Z+jets)
                            for ibin in range(1,hist[h][typ].GetNbinsX()+2):
                                for jbin in range(1,hist[h][typ].GetNbinsX()+2):
                                    print("hT = ",hist[h][typ].GetBinContent(ibin,jbin))
                                    print("hL = ",hist[h.replace("nT","nL")][typ].GetBinContent(ibin,jbin))
                                    #if hist[h][typ].GetBinContent(ibin,jbin) > hist[h.replace("nT","nL")][typ].GetBinContent(ibin,jbin):
                                    #    print "----> OBS!! <-----"
                                    #    hist[h][typ].SetBinContent(ibin,jbin,hist[h][typ].GetBinContent(ibin-1,jbin-1))
                                    #    hist[h.replace("nT","nL")][typ].SetBinContent(ibin,jbin,hist[h.replace("nT","nL")][typ].GetBinContent(ibin-1,jbin-1))
                                    #    print "hT (fixed) = ",hist[h][typ].GetBinContent(ibin,jbin)
                                    #    print "hL (fixed) = ",hist[h.replace("nT","nL")][typ].GetBinContent(ibin,jbin)


                            hEff[typ] = TEfficiency(hist[h][typ],hist[h.replace("nT","nL")][typ])
                            hEff[typ].SetName(hist[h][typ].GetName().replace("nT","eff")+"_"+addname)
                            hEff[typ].SetStatisticOption(TEfficiency.kBBayesian)
                            hEff[typ].SetConfidenceLevel(0.68);
                            hEff[typ].SetBetaAlpha(1.)
                            hEff[typ].SetBetaBeta(1.) 
                            hEff[typ].SetDirectory(0)

                            #canv.append(TCanvas("nL_%s"%typ,"nL_%s"%typ,1))
                            #canv[-1].cd()
                            #hist[h.replace("nT","nL")][typ].Draw("colz text")
                            #canv.append(TCanvas("nT_%s"%typ,"nT_%s"%typ,1))
                            #canv[-1].cd()
                            #hist[h][typ].Draw("colz text")
                            #canv.append(TCanvas("eff_%s"%typ,"eff_%s"%typ,1))
                            #canv[-1].cd()
                            #hEff[typ].Draw("colz text")

                            hEff_xsecup[typ] = TEfficiency(hist_xsecup[h][typ],hist_xsecup[hnameL][typ])
                            hEff_xsecup[typ].SetName(hist_xsecup[h][typ].GetName().replace("nT","eff")+"_xsecup")
                            hEff_xsecup[typ].SetStatisticOption(TEfficiency.kBBayesian)
                            hEff_xsecup[typ].SetConfidenceLevel(0.68);
                            hEff_xsecup[typ].SetBetaAlpha(1.)
                            hEff_xsecup[typ].SetBetaBeta(1.) 
                            hEff_xsecup[typ].SetDirectory(0)

                            hEff_xsecdw[typ] = TEfficiency(hist_xsecdw[h][typ],hist_xsecdw[hnameL][typ])
                            hEff_xsecdw[typ].SetName(hist_xsecdw[h][typ].GetName().replace("nT","eff")+"_xsecdw")
                            hEff_xsecdw[typ].SetStatisticOption(TEfficiency.kBBayesian)
                            hEff_xsecdw[typ].SetConfidenceLevel(0.68);
                            hEff_xsecdw[typ].SetBetaAlpha(1.)
                            hEff_xsecdw[typ].SetBetaBeta(1.) 
                            hEff_xsecdw[typ].SetDirectory(0)
                        else:
                            #hEff[typ] = hist[h][typ].Clone(hist[h][typ].GetName().replace("nT","eff")+"_"+addname)
                            #hEff_xsecup[typ] = hist_xsecup[h][typ].Clone(hist_xsecup[h][typ].GetName().replace("nT","eff")+"_xsecup")
                            #hEff_xsecdw[typ] = hist_xsecdw[h][typ].Clone(hist_xsecdw[h][typ].GetName().replace("nT","eff")+"_xsecdw")

                            htemp = getEff(hist[h][typ],hist[hnameL][typ])
                            htemp.SetDirectory(0)
                            hEff[typ] = htemp.Clone(hist[h][typ].GetName().replace("nT","eff")+"_"+addname)

                            htemp2 = getEff(hist_xsecup[h][typ],hist_xsecup[hnameL][typ])
                            htemp2.SetDirectory(0)
                            hEff_xsecup[typ] = htemp2.Clone(hist_xsecup[h][typ].GetName().replace("nT","eff")+"_xsecup")

                            htemp3 = getEff(hist_xsecdw[h][typ],hist_xsecdw[hnameL][typ])
                            htemp3.SetDirectory(0)
                            hEff_xsecdw[typ] = htemp3.Clone(hist_xsecdw[h][typ].GetName().replace("nT","eff")+"_xsecdw")
                            #hEff_xsecup[typ] = getEff(hist_xsecup[h][typ],hist_xsecup[hnameL][typ])#.Clone(hist_xsecup[h][typ].GetName().replace("nT","eff")+"_xsecup")#.Divide(hist_xsecup[h][typ],hist_xsecup[hnameL][typ],1.,1.,'b')
                            #hEff_xsecdw[typ] = getEff(hist_xsecdw[h][typ],hist_xsecdw[hnameL][typ])#.Clone(hist_xsecdw[h][typ].GetName().replace("nT","eff")+"_xsecdw")#.Divide(hist_xsecdw[h][typ],hist_xsecdw[hnameL][typ],1.,1.,'b')
                            hEff[typ].SetDirectory(0)
                            hEff_xsecup[typ].SetDirectory(0)
                            hEff_xsecdw[typ].SetDirectory(0)
                            hEff[typ].SetLineColor(d_samp[typ]['f_color'])
                            hEff[typ].SetMarkerColor(d_samp[typ]['f_color'])
                            hEff[typ].SetFillColor(d_samp[typ]['f_color'])
                            if "data" in typ and estfake:
                                hEff[typ].SetLineColor(kBlack)
                                hEff[typ].SetMarkerColor(kBlack)
                                hEff[typ].SetFillColor(kBlack)
                            #if typ == "Wjets":
                            #    for i in range(0,hEff[typ].GetNbinsX()+1):
                            #        print "bin ",i
                            #        print "nT = ",hist[h][typ].GetBinContent(i)
                            #        print "nL = ",hist[hnameL][typ].GetBinContent(i)
                            #        print "f  = ",hEff[typ].GetBinContent(i)
                            #asym = TGraphAsymmErrors();
                            #print "-->Dividing hist[%s][%s] = %.2f with %i bins on hist[%s][%s] = %.2f with %i bins" %(h,typ,hist[h][typ].Integral(),hist[h][typ].GetNbinsX(),hnameL,typ,hist[hnameL][typ].Integral(),hist[hnameL][typ].GetNbinsX())
                            #if typ == "Wjets":
                            #for ibin in range(1,hist[h][typ].GetNbinsX()+1):
                            #    print "-->Low edge %.0f, high edge %.0f" %(hist[h][typ].GetXaxis().GetBinLowEdge(ibin),hist[h][typ].GetXaxis().GetBinLowEdge(ibin+1))
                            #    print "-->Dividing %s (bin %i): %.2f / %.2f" %(typ,ibin,hist[h][typ].GetBinContent(ibin),hist[hnameL][typ].GetBinContent(ibin))
                            #hist[h][typ].ClearUnderflowAndOverflow()
                            #hist[hnameL][typ].ClearUnderflowAndOverflow()
                            #asym.Divide(hist[h][typ],hist[hnameL][typ],"cl=0.683 b(1,1) mode")
                            #ni = 0
                            #for i in range(0,hEff[typ].GetNbinsX()+1):
                            #    print "Typ: %s, Bin %i : %.2f" %(typ,i,hEff[typ].GetBinContent(i))
                            #    print "nL = %.2f, nT = %.2f" %(hist[h][typ].GetBinContent(i),hist[hnameL][typ].GetBinContent(i))
                            #    #print "asym", asym.GetErrorY(i)
                            #    #print "asym high", asym.GetErrorYhigh(i)
                            #    #print "asym low", asym.GetErrorYlow(i)
                            #    #print "heff", hEff[typ].GetBinError(i+1)
                            #    if hEff[typ].GetBinContent(i+1):
                            #        print "Setting error for bin %i to %.2f" %(i+1,asym.GetErrorY(ni))
                            #        hEff[typ].SetBinError(i+1,asym.GetErrorY(ni))
                            #        ni += 1
                            #    else:
                            #        hEff[typ].SetBinError(i+1,0.0)

            print("datakey = ",datakey)
            if estfake and "data" in typ:
                olddatakey = datakey
                if alldata: datakey = "alldata"
                elif comb1516: datakey = "data15-16" 
                print("Getting fakes!")
                for addi in ["","UP","DW"]:
                    if addi:
                        fake_hist_name = h +"_estfake_"+addi
                    else:
                        fake_hist_name = h +"_estfake"
                    if "trueREAL" in fake_hist_name: fake_hist_toget = fake_hist_name.replace("_trueREAL","")
                    elif "trueFAKE" in fake_hist_name: fake_hist_toget = fake_hist_name.replace("_trueFAKE","")
                    else: fake_hist_toget = fake_hist_name
                    hist_midl = (file.Get(fake_hist_toget)).Clone("hist_midl")
                    print("Adding fake histogram with name ", fake_hist_name)
                    print("Getting fake histogram with name %s and events %.2f" %(fake_hist_toget,hist_midl.Integral()))
                    if doRB: 
                        hist_midl = doRebin(hist_midl,rebarray)
                    elif doRBSimple:
                        hist_midl = hist_midl.Rebin(rebfac)
                        
                    if not fake_hist_name in list(hist.keys()):
                        hist[fake_hist_name] = {}
                        if not datakey in list(hist[fake_hist_name].keys()):
                            hist[fake_hist_name][datakey] = hist_midl.Clone(fake_hist_name+"_"+typ)
                            hist[fake_hist_name][datakey].SetDirectory(0)
                            print("INFO \t Fakes : Creating %s %s from %s" %(h,typ,fn))
                            print("INFO \t Fakes : Starting with %.2f events" %(hist[fake_hist_name][datakey].Integral()))
                            
                    else:
                        print("INFO \t Fakes : Adding %s %s from %s" %(h,typ,fn))
                        hist[fake_hist_name][datakey].Add(hist_midl.Clone(fake_hist_name+"_"+typ))
                        print("INFO \t Fakes : Now with %.2f events" %(hist[fake_hist_name][datakey].Integral()))
                fake_hist_name = fake_hist_name = h +"_estfake"
                datakey = olddatakey
            print("Closing file %s" %file.GetName())
            file.Close()
    

    nfakefiles = 0
    if estfake and inputdir_FAKE:
        
        for nf in FAKEfiles:
            if datakey == "data15-16" and not ("data16" in nf or "data15" in nf): continue
            if datakey == "data17" and not "data17" in nf: continue
            if datakey == "data18" and not "data18" in nf: continue

            tfile = TFile(nf)

            ch = h.split("_")[3]
            var = h.split("_")[2]
            
            if "zveto20_met40" in h: addstr = "ZV_MET"
            elif "zveto20" in h: addstr = "ZV"
            else: addstr = ""
            if addstr:
                hname = "h_lep_%s_%s_VR2LFAKES_%s_GradientLoose_TRG_NO_fake_pT_eta" %(varmap[var],ch,addstr)
            else:
                hname = "h_lep_%s_%s_VR2LFAKES_GradientLoose_TRG_NO_fake_pT_eta" %(varmap[var],ch)

            fake_hist_name = hname

            for vari in ["","wgtup","wgtdw","statup","statdw"]: #"_fake_pT_eta",
                if vari == "_fake_pT_eta":
                    hname_new = hname.replace(vari,"")
                    data_hist_name = hname_new
                elif vari: hname_new = hname+"_"+vari
                else: hname_new = hname

                #print "getting ",hname_new
                if not hname_new in list(hist.keys()):
                    hist[hname_new] = {}
                if not datakey in list(hist[hname_new].keys()):
                    hist[hname_new][datakey] = TH1F()
                
                hist_midl = tfile.Get(hname_new).Clone("hist_midl")
                if doRB: 
                    if is2D: hist_midl = doRebin2D(hist_midl,rebx,reby)
                    else: hist_midl = doRebin(hist_midl,rebarray)
                elif doRBSimple:
                    hist_midl = hist_midl.Rebin(rebfac)
                
                hist_midl.SetDirectory(0)

                if nfakefiles == 0:
                    hist[hname_new][datakey] = hist_midl.Clone(hname_new)
                else:
                    hist[hname_new][datakey].Add(hist_midl)
                    
                hist[hname_new][datakey].SetDirectory(0)

            tfile.Close()
            
            nfakefiles += 1

        for hname_new in list(hist.keys()):
            if not "VR2LFAKES" in hname_new: continue
            if "EE" in hname_new:
                print("Scalingg", hname_new)
                print("before", hist[hname_new][datakey].Integral())
                #hist[hname_new][datakey].Scale(0.8)
                print("after", hist[hname_new][datakey].Integral())        

    print("Added %i files for histogram %s" %(nfakefiles,h))
   
if alldata or comb1516:
    if alldata: datakey = "alldata"
    elif comb1516: datakey = "data15-16"

    print("Doing eff afterwards!")
    for h in list(hist.keys()):
        if "_REAL" in h and "trueFAKE" in h: continue
        if "_FAKE" in h and "trueREAL" in h: continue
        for typ in list(hist[h].keys()):
            if not "nT" in h: continue
            hnameL = h.replace("nT","nL")
            print("hEff = ",hist[h][typ].GetName().replace("nT","eff"))
    
            print("typ", typ)
            print("Passed", hist[h][typ].GetBinContent(1))
            print("Total", hist[h.replace("nT","nL")][typ].GetBinContent(1))
    
            if is2D:
                hEff[typ] = TEfficiency(hist[h][typ],hist[h.replace("nT","nL")][typ])
                hEff[typ].SetName(hist[h][typ].GetName().replace("nT","eff")+"_"+addname)
                hEff[typ].SetStatisticOption(TEfficiency.kBBayesian)
                hEff[typ].SetConfidenceLevel(0.68);
                hEff[typ].SetBetaAlpha(1.)
                hEff[typ].SetBetaBeta(1.) 

                hEff_xsecup[typ] = TEfficiency(hist_xsecup[h][typ],hist_xsecup[hnameL][typ])
                hEff_xsecup[typ].SetName(hist_xsecup[h][typ].GetName().replace("nT","eff")+"_xsecup")
                hEff_xsecup[typ].SetStatisticOption(TEfficiency.kBBayesian)
                hEff_xsecup[typ].SetConfidenceLevel(0.68);
                hEff_xsecup[typ].SetBetaAlpha(1.)
                hEff_xsecup[typ].SetBetaBeta(1.) 

                hEff_xsecdw[typ] = TEfficiency(hist_xsecdw[h][typ],hist_xsecdw[hnameL][typ])
                hEff_xsecdw[typ].SetName(hist_xsecdw[h][typ].GetName().replace("nT","eff")+"_xsecdw")
                hEff_xsecdw[typ].SetStatisticOption(TEfficiency.kBBayesian)
                hEff_xsecdw[typ].SetConfidenceLevel(0.68);
                hEff_xsecdw[typ].SetBetaAlpha(1.)
                hEff_xsecdw[typ].SetBetaBeta(1.) 
            else:
                hEff[typ] = hist[h][typ].Clone(hist[h][typ].GetName().replace("nT","eff")+"_"+addname)
                hEff_xsecup[typ] = hist_xsecup[h][typ].Clone(hist_xsecup[h][typ].GetName().replace("nT","eff")+"_xsecup")
                hEff_xsecdw[typ] = hist_xsecdw[h][typ].Clone(hist_xsecdw[h][typ].GetName().replace("nT","eff")+"_xsecdw")  
                hEff[typ].Divide(hist[h][typ],hist[hnameL][typ],1.,1.,'b')
                hEff_xsecup[typ].Divide(hist_xsecup[h][typ],hist_xsecup[hnameL][typ],1.,1.,'b')
                hEff_xsecdw[typ].Divide(hist_xsecdw[h][typ],hist_xsecdw[hnameL][typ],1.,1.,'b')
                hEff[typ].SetDirectory(0)
                hEff_xsecup[typ].SetDirectory(0)
                hEff_xsecdw[typ].SetDirectory(0)
                hEff[typ].SetLineColor(d_samp[typ]['f_color'])
                hEff[typ].SetMarkerColor(d_samp[typ]['f_color'])
                hEff[typ].SetFillColor(d_samp[typ]['f_color'])
                if "data" in typ and estfake:
                    hEff[typ].SetLineColor(kBlack)
                    hEff[typ].SetMarkerColor(kBlack)
                    hEff[typ].SetFillColor(kBlack)
                if typ == "Wjets":
                    for i in range(0,hEff[typ].GetNbinsX()+1):
                        print("bin ",i)
                        print("nT = ",hist[h][typ].GetBinContent(i))
                        print("nL = ",hist[hnameL][typ].GetBinContent(i))
                        print("f  = ",hEff[typ].GetBinContent(i))
                asym = TGraphAsymmErrors();
                print("-->Dividing hist[%s][%s] = %.2f with %i bins on hist[%s][%s] = %.2f with %i bins" %(h,typ,hist[h][typ].Integral(),hist[h][typ].GetNbinsX(),hnameL,typ,hist[hnameL][typ].Integral(),hist[hnameL][typ].GetNbinsX()))
                #if typ == "Wjets":
                for ibin in range(1,hist[h][typ].GetNbinsX()+1):
                    print("-->Low edge %.0f, high edge %.0f" %(hist[h][typ].GetXaxis().GetBinLowEdge(ibin),hist[h][typ].GetXaxis().GetBinLowEdge(ibin+1)))
                    print("-->Dividing %s (bin %i): %.2f / %.2f" %(typ,ibin,hist[h][typ].GetBinContent(ibin),hist[hnameL][typ].GetBinContent(ibin)))
                hist[h][typ].ClearUnderflowAndOverflow()
                hist[hnameL][typ].ClearUnderflowAndOverflow()
                asym.Divide(hist[h][typ],hist[hnameL][typ],"cl=0.683 b(1,1) mode")
                ni = 0
                for i in range(0,hEff[typ].GetNbinsX()+1):
                    print("Typ: %s, Bin %i : %.2f" %(typ,i,hEff[typ].GetBinContent(i)))
                    print("nL = %.2f, nT = %.2f" %(hist[h][typ].GetBinContent(i),hist[hnameL][typ].GetBinContent(i)))
                    #print "asym", asym.GetErrorY(i)
                    #print "asym high", asym.GetErrorYhigh(i)
                    #print "asym low", asym.GetErrorYlow(i)
                    #print "heff", hEff[typ].GetBinError(i+1)
                    if hEff[typ].GetBinContent(i+1):
                        print("Setting error for bin %i to %.2f" %(i+1,asym.GetErrorY(ni)))
                        hEff[typ].SetBinError(i+1,asym.GetErrorY(ni))
                        ni += 1
                    else:
                        hEff[typ].SetBinError(i+1,0.0)

    

#for h in hist.keys():
#    for typ in hist[h].keys():
#        if not typ in ['Zeejets', 'Zmmjets', 'Zttjets']: continue
#        print("Scaling %s "%typ)
#        #hist[h][typ].Scale(46.5)
#        print("before", hist[h][typ].Integral())
#        hist[h][typ].Scale(1.0e3)
#        #hist[h][typ].Scale(0.971055475361775)
        #if typ in ["lowMassDY","Vgamma"]:
        #   hist[h][typ].Scale(59937.2/44307.4)
#        print("after", hist[h][typ].Integral())
        # if "data" in typ: continue
        # if datakey == "data15-16": hist[h][typ].Scale(0.8265)
        # elif datakey == "data17": hist[h][typ].Scale(1/0.8265)
        # elif datakey == "data18": hist[h][typ].Scale(0.3)

if 0:#("EE" in h or "MM" in h or "EM" in h or "ALL" in h):# and not estfake:
    for h in list(hist.keys()):
        for typ in list(hist[h].keys()):
            if "data" in typ: continue
            #if "data18" in datakey:  hist[h][typ].Scale(0.30)
            continue
            if "EE" in h or ("ALL" in h and ("2L23" in h or "2L24" in h or "2L25" in h or "2L27" in h)):
                print("scaling ", typ)
                if "data15-16" in datakey:
                    #print "before", hist[h][typ].Integral()
                    if not estfake: hist[h][typ].Scale(0.8417)
                    else: hist[h][typ].Scale(0.8)
                if "data18" in datakey:
                    hist[h][typ].Scale(0.3344)
                    #print "after", hist[h][typ].Integral()
                elif "data17" in datakey: hist[h][typ].Scale(0.2288)
            # elif "MM" in h and not estfake:
            #     if "data15-16" in datakey: hist[h][typ].Scale(0.4975)
            #     elif "data17" in datakey: hist[h][typ].Scale(0.2288)
            elif "EM" in h and not estfake:
                print("scaling")
                if "data15-16" in datakey: hist[h][typ].Scale(0.4700)
                elif "data18" in datakey: hist[h][typ].Scale(1.86)
                elif "data17" in datakey: hist[h][typ].Scale(0.2288)
                    


#sys.exit()
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
leg = {}

pad1 = []
pad2 = []
pad3 = []
stack = {}

totals_mc = []
errors_mc = []

totals = []
errors = []

h_ratio = []
h_scaled = []

h_fract_mc = []
h_fract = []
stck_fk = []
stck = []
leg1 = []
leg2 = []
            
arr = {}
if is2D:
    for h in list(hist.keys()):
        totals_mc = []
        errors_mc = []
        
        totals = []
        errors = []
        
        h_ratio = []
        h_scaled = []
        if "TCF" in h: 
            ntyp = -1
            for typ in list(hist[h].keys()):
                if "data" in typ: continue
                ntyp += 1
                binlabel = []
                if not "_njet_" in h and not "metsig" in h: hist[h][typ] = doRebin2D(hist[h][typ],array.array('f',[0,20,40,60,100,1000]),array.array('f',[0,1,2,3,4,5,6,7,8,9,10,11]))
                elif "metsig" in h: hist[h][typ] = doRebin2D(hist[h][typ],array.array('f',[0,2,4,6,8,10,100]),array.array('f',[0,1,2,3,4,5,6,7,8,9,10,11]))
                for i in range(1,hist[h][typ].GetNbinsX()+1):
                    if ntyp == 0:
                        #totals.append({"CF":0,"CO":0,"HF":0,"LF":0,"ALL":0})
                        #errors.append({"CF":0,"CO":0,"HF":0,"LF":0,"ALL":0})
                        totals.append({
                            "Unknown":0,
                            "KnownUnknown":0,
                            "IsoElectron":0,
                            "ChargeFlipIsoElectron":0,
                            "PromptMuon":0,
                            "PromptPhotonConversion":0,
                            "ElectronFromMuon":0,
                            "TauDecay":0,
                            "BHadronDecay":0,
                            "CHadronDecay":0,
                            "LightFlavorDecay":0,
                            "DATA":0,
                            "ALL":0
                        })
                        errors.append({
                            "Unknown":0,
                            "KnownUnknown":0,
                            "IsoElectron":0,
                            "ChargeFlipIsoElectron":0,
                            "PromptMuon":0,
                            "PromptPhotonConversion":0,
                            "ElectronFromMuon":0,
                            "TauDecay":0,
                            "BHadronDecay":0,
                            "CHadronDecay":0,
                            "LightFlavorDecay":0,
                            "DATA":0,
                            "ALL":0
                        })
                        if len(totals_mc) < i:
                            totals_mc.append({})
                            errors_mc.append({})
                        print("addin ALL")
                        totals_mc[i-1]["ALL"] = 0
                        errors_mc[i-1]["ALL"] = 0
                    totals_mc[i-1][typ] = 0
                    errors_mc[i-1][typ] = 0
                    binlabel.append("[%i-%i>"%(hist[h][typ].GetXaxis().GetBinLowEdge(i),hist[h][typ].GetXaxis().GetBinLowEdge(i+1)))
                    
                    for j in range(1,hist[h][typ].GetNbinsY()+2):
                        print(j)
                        lab = hist[h][typ].GetYaxis().GetBinLabel(j)
                        if not lab in list(totals[i-1].keys()): continue
                        if i == 5 and "nT" in h: print("%s bin (%i,%s) = %.2f" %(typ,i,lab,hist[h][typ].GetBinContent(i,j)))
                        totals[i-1]["ALL"] +=  hist[h][typ].GetBinContent(i,j)
                        errors[i-1]["ALL"] += (hist[h][typ].GetBinError(i,j)*hist[h][typ].GetBinError(i,j))
                        totals[i-1][lab]   +=  hist[h][typ].GetBinContent(i,j)
                        errors[i-1][lab]   += (hist[h][typ].GetBinError(i,j)*hist[h][typ].GetBinError(i,j))    

                        totals_mc[i-1]["ALL"] +=  hist[h][typ].GetBinContent(i,j)
                        errors_mc[i-1]["ALL"] += (hist[h][typ].GetBinError(i,j)*hist[h][typ].GetBinError(i,j))
                        totals_mc[i-1][typ] +=  hist[h][typ].GetBinContent(i,j)
                        errors_mc[i-1][typ] += (hist[h][typ].GetBinError(i,j)*hist[h][typ].GetBinError(i,j))

            max_totals = -999
            for i in range(0,len(totals)):
                for lab in list(totals[i].keys()):
                    if lab == "ALL" or lab == "DATA": continue
                    if totals[i][lab] > max_totals:
                        max_totals = totals[i][lab]

            if "trueFAKE" in h:
                h_fract.append(getFract(totals,errors,binlabel,h))
            else:
                h_fract.append(getAbs(totals,errors,binlabel,h))
            ndraw = 0
            canv.append(TCanvas("c_ratio_%s"%h,"c_ratio_%s"%h,1))
            pad3.append(TPad("pad3_mc","pad3_mc",0.8,0.0,1.0,1.0,21))
            pad3[-1].Draw();
            pad3[-1].SetFillColor(0);
            pad1.append(TPad("pad1_fk","pad1_fk",0.0,0.5,0.8,1.0,21))
            pad1[-1].Draw();
            pad1[-1].SetRightMargin(0.01);
            pad2.append(TPad("pad2_mc","pad2_mc",0.0,0.0,0.8,0.5,22))
            pad2[-1].SetRightMargin(0.01);
            pad2[-1].Draw();
            pad2[-1].SetFillColor(0);
            pad2[-1].SetTopMargin(0.05);
            pad2[-1].SetBottomMargin(0.2);
            pad1[-1].SetFillColor(0);
            pad1[-1].SetBottomMargin(0.1);
            pad1[-1].SetTopMargin(0.09);
            pad1[-1].cd() 
            leg1.append(TLegend(0.02149,0.536025,0.953,0.951477))
            leg1[-1].SetTextSize(0.1)
            leg1[-1].SetTextColor(1)
            leg2.append(TLegend(0.0186246,0.101266,0.95702,0.476793))
            leg2[-1].SetTextSize(0.1)
            if not "trueFAKE" in h:
                pad2[-1].SetLogy(1);
            pad2[-1].SetGridx()
            pad2[-1].SetGridy()
            if not "trueFAKE" in h:
                pad1[-1].SetLogy(1);
            pad1[-1].SetGridx()
            pad1[-1].SetGridy()
            stck_fk.append(THStack("stck_fk_"+h,"stck_fk_"+h))
            for key in sorted(h_fract[-1].keys()):
                pad1[-1].cd()
                #if ndraw == 0: 
                #    h_fract[-1][key].GetYaxis().SetRangeUser(1e-2,1.0)
                #    h_fract[-1][key].GetYaxis().SetTitle("Fraction")
                #    h_fract[-1][key].GetXaxis().SetTitle("mT_{2} bin [GeV]")
                #    h_fract[-1][key].GetXaxis().SetRangeUser(0,5)
                #    h_fract[-1][key].Draw("e1")
                #else: h_fract[-1][key].Draw("e1 same")
                MMfile = TFile("MMinput_frac.root","UPDATE")
                MMfile.cd()
                
                h_fract[-1][key].Write(h_fract[-1][key].GetName()+postfix)
                MMfile.Close()

                stck_fk[-1].Add(h_fract[-1][key])
                leg1[-1].AddEntry(h_fract[-1][key],d_samp[key]['leg'],"f")
                ndraw += 1
            stck_fk[-1].Draw("nostackb hist")
            stck_fk[-1].GetYaxis().SetTitle("Fraction")
            stck_fk[-1].SetMinimum(0)
            if "trueFAKE" in h: 
                stck_fk[-1].SetMaximum(1.0)
            else:
                stck_fk[-1].SetMinimum(1.0)
                stck_fk[-1].SetMaximum(max_totals)
            stck_fk[-1].GetYaxis().SetLabelSize(0.06)
            stck_fk[-1].GetXaxis().SetLabelSize(0.09)
            stck_fk[-1].GetXaxis().SetTitleSize(0.06)
            stck_fk[-1].GetXaxis().SetTitleOffset(1.55)
            stck_fk[-1].GetYaxis().SetTitleOffset(0.6)
            stck_fk[-1].GetYaxis().SetTitleSize(0.06)
            if not "_njet_" in h and not "metsig" in h: stck_fk[-1].GetXaxis().SetRangeUser(0,5)
            #canv[-1].SetLogy()
            xr = pad1[-1].GetX2()-pad1[-1].GetX1()
            for i in range(0,len(totals_mc)):
                xpos = (stck_fk[-1].GetXaxis().GetBinCenter(i+1)-pad1[-1].GetX1())/ xr
                myText(xpos, 0.6, "%.0f"%totals_mc[i]["ALL"], 0.08, kBlack)
            region = ""
            for r in list(d_reg.keys()):
                if r in h:
                    region = r
                    break
            if region:
                myText(0.25, 0.8, "Reg.:"+region.replace("REAL","")+", Ch: "+h.split("_")[5], text_size*0.15, kBlack)
                #myText(0.25, 0.8, "Cuts: "+d_reg[region]["descr"], text_size*0.15, kBlack)
            pad1[-1].Update()
            
            max_totals_mc = -999
            for i in range(0,len(totals_mc)):
                for lab in list(totals_mc[i].keys()):
                    if lab == "ALL" or lab == "DATA": continue
                    if totals_mc[i][lab] > max_totals_mc:
                        max_totals_mc = totals_mc[i][lab]
            
            if "trueFAKE" in h:
                h_fract_mc.append(getFract(totals_mc,errors_mc,binlabel,h))
            else:
                h_fract_mc.append(getAbs(totals_mc,errors_mc,binlabel,h))
            ndraw = 0
            stck.append(THStack("stck_"+h,"stck_"+h))
            for key in list(h_fract_mc[-1].keys()):
                pad2[-1].cd()
                stck[-1].Add(h_fract_mc[-1][key])
                ndraw += 1
                leg2[-1].AddEntry(h_fract_mc[-1][key],d_samp[key]['leg'],"f")
            stck[-1].Draw("nostackb hist")
            #stck[-1].GetYaxis().SetRangeUser(1e-2,1.0)
            stck[-1].SetMinimum(0)
            if "trueFAKE" in h: 
                stck[-1].SetMaximum(1.0)
            else:
                stck[-1].SetMinimum(1.0)
                stck[-1].SetMaximum(max_totals_mc)
            #canv[-1].SetBottomMargin(0.2)
            stck[-1].GetYaxis().SetLabelSize(0.06)
            stck[-1].GetXaxis().SetLabelSize(0.08)
            stck[-1].GetXaxis().SetTitleSize(0.06)
            stck[-1].GetXaxis().SetTitleOffset(1.55)
            stck[-1].GetYaxis().SetTitleOffset(0.6)
            stck[-1].GetYaxis().SetTitleSize(0.06)
            stck[-1].GetYaxis().SetTitle("Fraction")
            if "mT2" in h: stck[-1].GetXaxis().SetTitle("mT_{2} bin [GeV]")
            elif "mTW" in h: stck[-1].GetXaxis().SetTitle("mT_{W} bin [GeV]")
            elif "metsig" in h: stck[-1].GetXaxis().SetTitle("E_{T}^{miss}-sign.")
            if not "_njet_" in h and not "metsig" in h:
                stck[-1].GetXaxis().SetRangeUser(0,5)
            elif "_njet_" in h:
                stck[-1].GetXaxis().SetTitle("N_{light jets}^{p_{T}>20, |#eta|<2.4}")
            #canv[-1].SetLogy()

            pad3[-1].cd()
            leg1[-1].Draw()
            leg2[-1].Draw()

            canv[-1].Update()

        elif "origin" in h or "type" in h:
            ntyp = -1
            calcFR = True
            for typ in list(hist[h].keys()):
                if "data" in typ: continue
                if not ("Vgamma" in typ or "Wjets" in typ or "ttbar" in typ): continue
                ntyp += 1
                canv.append(TCanvas("c_%s_%s"%(typ,h),"c_%s_%s"%(typ,h),1))
                canv[-1].cd()
                
                h_scaled.append(hist[h][typ].Clone(hist[h][typ].GetName()+"_scaled"))
                if not "type_origin" in h:
                    for nx in range(1,h_scaled[-1].GetNbinsX()+2):
                        intgrl = h_scaled[-1].Integral(nx,nx,0,h_scaled[-1].GetNbinsY()+1)
                        if intgrl <= 0: continue
                        for ny in range(1,h_scaled[-1].GetNbinsY()+2):
                            if h_scaled[-1].GetBinContent(h_scaled[-1].GetBin(nx,ny)) < 0: continue
                            h_scaled[-1].SetBinContent(h_scaled[-1].GetBin(nx,ny),h_scaled[-1].GetBinContent(h_scaled[-1].GetBin(nx,ny))/intgrl)
                            h_scaled[-1].SetBinError(h_scaled[-1].GetBin(nx,ny),h_scaled[-1].GetBinError(h_scaled[-1].GetBin(nx,ny))/intgrl)
                elif h_scaled[-1].Integral(0,h_scaled[-1].GetNbinsX()+1,0,h_scaled[-1].GetNbinsY()+1) > 0:
                    h_scaled[-1].Scale(1.0/h_scaled[-1].Integral(0,h_scaled[-1].GetNbinsX()+1,0,h_scaled[-1].GetNbinsY()+1))
                    
                h_scaled[-1].SetMarkerColor(kBlack)
                if not "type_origin" in h:
                    h_scaled[-1].GetXaxis().SetTitle("p_{T} [GeV]")
                    h_scaled[-1].GetZaxis().SetRangeUser(0,0.35)
                    h_scaled[-1].GetYaxis().SetRangeUser(0,25)
                    canv[-1].SetLogx()
                else:
                    h_scaled[-1].GetXaxis().SetTitle("type")
                    h_scaled[-1].SetMarkerSize(3)
                    h_scaled[-1].GetZaxis().SetRangeUser(0,0.8)
                    h_scaled[-1].GetYaxis().SetRangeUser(0,10)
                    h_scaled[-1].GetXaxis().SetRangeUser(0,6)
                    
                h_scaled[-1].GetYaxis().SetTitle("Origin" if 'origin' in h else "Type")
                h_scaled[-1].Draw("colz text")
                h_scaled[-1].SetDirectory(0)
                
                if not "nL" in h: 
                    canv.append(TCanvas("c_%s_%s"%(typ,h.replace("nT","eff")),"c_%s_%s"%(typ,h.replace("nT","eff")),1))
                    canv[-1].cd()
                    h_ratio.append(hist[h][typ].Clone(hist[h][typ].GetName().replace("nT","eff")))
                    h_ratio[-1].Divide(hist[h][typ],hist[h.replace("nT","nL")][typ],1.0,1.0,"b")
                    asym = TGraphAsymmErrors();
                    print("<--Dividing hist[%s][%s] = %.2f on hist[%s][%s] = %.2f" %(h,typ,hist[h][typ].Integral(),h.replace("nT","nL"),typ,hist[h.replace("nT","nL")][typ].Integral()))
                    asym.Divide(hist[h][typ],hist[h.replace("nT","nL")][typ],"cl=0.683 b(1,1) mode")
                    #ni += 1
                    ni = 0
                    for i in range(0,h_ratio[-1].GetNbinsX()+1):
                        if h_ratio[-1].GetBinContent(i+1):
                            h_ratio[-1].SetBinError(i+1,asym.GetErrorY(ni))
                            ni += 1
                        else:
                            h_ratio[-1].SetBinError(i+1,0.0)
                        
                    h_ratio[-1].GetZaxis().SetRangeUser(0,0.7)
                    if not "type_origin" in h:
                        canv[-1].SetLogx()
                        h_ratio[-1].GetYaxis().SetRangeUser(0,25)
                        h_ratio[-1].GetXaxis().SetTitle("p_{T} [GeV]")
                    else:
                        h_ratio[-1].GetXaxis().SetTitle("Type")
                        h_ratio[-1].SetMarkerSize(3)
                        h_ratio[-1].GetYaxis().SetRangeUser(0,10)
                        h_ratio[-1].GetXaxis().SetRangeUser(0,6)
                    h_ratio[-1].SetMarkerColor(kBlack)
                    h_ratio[-1].GetYaxis().SetTitle("Origin" if 'origin' in h else "Type")
                    h_ratio[-1].Draw("colz text")
                    h_ratio[-1].SetDirectory(0)
        elif "pT_eta" in h:
            print("h = ",h) 
            if "nT" in h or "nL" in h: 
                typ = datakey
                print("Plotting ", hist[h][typ].GetName())
                canv.append(TCanvas("c_2d_%s"%h,"c_2d_%s"%h,1))
                canv[-1].SetLogx()
                canv[-1].SetRightMargin(0.14);
                hist[h][typ].Draw("colz text e")
                print("binc = ", hist[h][typ].GetBinContent(15))
                canv[-1].Update()
                print(typ)
                MMfile = TFile("MMinput.root","UPDATE")
                MMfile.cd()
                hist[h][typ].Write(hist[h][typ].GetName()+postfix)
                MMfile.Close()
                #h_ratio.append(TEfficiency(hist[h][typ],hist[h.replace("nT","nL")][typ]))
                #h_ratio[-1].SetStatisticOption(TEfficiency.kBBayesian)
                #h_ratio[-1].SetConfidenceLevel(0.68);
                #h_ratio[-1].SetBetaAlpha(1.)
                #h_ratio[-1].SetBetaBeta(1.) 
                #h_ratio.append(hist[h][typ].Clone(hist[h][typ].GetName().replace("nT","eff")))
                #print "nT = ", hist[h][typ].Integral()
                #print "nL = ", hist[h.replace("nT","nL")][typ].Integral()
                #h_ratio[-1].Divide(hist[h][typ],hist[h.replace("nT","nL")][typ],1.0,1.0,"b")
                #h_ratio[-1].Draw("colz text")
                #h_ratio[-1].SetDirectory(0)
    
    if typ in list(hEff.keys()) and ("pT_eta" in hEff[typ].GetName() or "lep_mu_pT" in hEff[typ].GetName()):
        for typ in list(hEff.keys()):
            canv.append(TCanvas("c_2deff_%s"%typ,"c_2deff_%s"%typ,1))
            if not "lep_mu_pT" in hEff[typ].GetName():
                canv[-1].SetLogx()
            canv[-1].SetRightMargin(0.14);
            #if not datakey in typ: continue
            hEff[typ].Draw("EY")
            canv[-1].Update()
            print(typ)
            hEff[typ+"_hist"] = hEff[typ].GetPaintedHistogram()
            try:
                hEff[typ+"_hist"].SetDirectory(0)
            except:
                print("Could not get eff histogram for %s" %typ)
                continue
            hEff[typ+"_hist"].SetName(hEff[typ].GetName()+"_hist")
            print(type(hEff[typ+"_hist"]))
            for i in range(0,hEff[typ+"_hist"].GetNbinsX()+1):
                for j in range(0,hEff[typ+"_hist"].GetNbinsY()+1):
                    ibin = hEff[typ+"_hist"].GetBin(i,j)
                    hEff[typ+"_hist"].SetBinError(ibin,(hEff[typ].GetEfficiencyErrorLow(ibin)+hEff[typ].GetEfficiencyErrorUp(ibin))/2.0)
            hEff[typ+"_hist"].SetMarkerColor(kBlack)
            if "REAL" in hEff[typ+"_hist"].GetName(): 
                hEff[typ+"_hist"].SetMarkerSize(0.8)
                hEff[typ+"_hist"].Draw("colz text e")
            else: 
                hEff[typ+"_hist"].Draw("colz text e")
            if "lep_mu_pT" in hEff[typ+"_hist"].GetName():
                hEff[typ+"_hist"].GetXaxis().SetTitle("mu")
                hEff[typ+"_hist"].GetYaxis().SetTitle("p_{T} [GeV]")
                hEff[typ+"_hist"].GetYaxis().SetRangeUser(0,50)
            else:
                hEff[typ+"_hist"].GetXaxis().SetTitle("p_{T} [GeV]")
                hEff[typ+"_hist"].GetYaxis().SetTitle("|#eta|")
            if "REAL" in hEff[typ+"_hist"].GetName(): hEff[typ+"_hist"].GetZaxis().SetTitle("Efficiency")
            else: hEff[typ+"_hist"].GetZaxis().SetTitle("Fake Rate")
            if "REAL" in hEff[typ+"_hist"].GetName(): hEff[typ+"_hist"].GetZaxis().SetRangeUser(0.65,1.0)
            else: hEff[typ+"_hist"].GetZaxis().SetRangeUser(0.0,1.0)
            canv[-1].Update()
            MMfile = TFile("MMinput.root","UPDATE")
            MMfile.cd()
            hEff[typ+"_hist"].Write(hEff[typ+"_hist"].GetName()+postfix)
            print("Writing " + hEff[typ+"_hist"].GetName() + " to file MMinput.root")
            MMfile.Close()
        
        h_allMC_2D = {}
        h_allMC_2D_pEff = {}
        nL_name = {}
        nT_name = {}
        for h in list(hist.keys()):
            hname = "nT" if "nT" in h else "nL"
            hname = hname+"_trueREAL" if "trueREAL" in h else (hname+"_trueFAKE" if "trueFAKE" in h else hname)
            h_allMC_2D[hname] = TH2F()
            h_allMC_2D[hname].SetName((h.replace("nT","eff") if "nT" in h else h.replace("nL","eff"))+"_totalMC")
            if "nL" in h: nL_name["trueFAKE" if "trueFAKE" in h else "trueREAL"] = h
            elif "nT" in h: nT_name["trueFAKE" if "trueFAKE" in h else "trueREAL"] = h
            print(h, hname)
            ntyp = 0
            for typ in list(hist[h].keys()):
                if "data" in typ: continue
                if not type(hist[h][typ]) is TH2F: 
                    print("Type %s is 1D" %typ)
                    continue
                if ntyp == 0 and (typ in subBkgs or "ALL" in subBkgs):
                    h_allMC_2D[hname] = hist[h][typ].Clone(hist[h][typ].GetName().replace(typ,addname))
                    ntyp += 1
                    continue
                elif ntyp == 0:
                    continue
                print(type(hist[h][typ]))
                print("Adding %s to total 2D" %typ)
                if typ in subBkgs or "ALL" in subBkgs: h_allMC_2D[hname].Add(hist[h][typ],1.0)
                ntyp += 1
        MMfile = TFile("MMinput.root","UPDATE")
        MMfile.cd()
        for key in list(h_allMC_2D.keys()):
            h_allMC_2D[key].Write(h_allMC_2D[key].GetName()+postfix)
        print("Writing " + h_allMC_2D[key].GetName() + " to file MMinput.root")
        MMfile.Close()

        if "nT_trueFAKE" in list(h_allMC_2D.keys()) and "nL_trueFAKE" in list(h_allMC_2D.keys()):
            canv.append(TCanvas("c_eff_trueFAKE_%s_totalMC"%h_allMC_2D["nT_trueFAKE"].GetName().replace("nT","eff"),"c_eff_trueFAKE_%s_totalMC"%h_allMC_2D["nT_trueFAKE"].GetName().replace("nT","eff"),1))
            canv[-1].cd()
            h_allMC_2D["eff_trueFAKE"] = TEfficiency(h_allMC_2D["nT_trueFAKE"],h_allMC_2D["nL_trueFAKE"])# = h_allMC_2D["nL_trueFAKE"].Clone(h_allMC_2D["nL_trueFAKE"].GetName().replace("nT","eff"))
            h_allMC_2D["eff_trueFAKE"].SetName((h_allMC_2D["nT_trueFAKE"].GetName().replace("nT","eff")))
            #h_allMC_2D["eff_trueFAKE"].SetDirectory(0)
            h_allMC_2D["eff_trueFAKE"].Draw("colz text e")
        if "nL_trueREAL" in list(h_allMC_2D.keys()) and "nL_trueREAL" in list(h_allMC_2D.keys()):
            canv.append(TCanvas("c_eff_trueREAL_%s_totalMC"%h_allMC_2D["nT_trueREAL"].GetName().replace("nT","eff"),"c_eff_trueREAL_%s_totalMC"%h_allMC_2D["nT_trueREAL"].GetName().replace("nT","eff"),1))
            canv[-1].cd()
            for ibin in range(1,h_allMC_2D["nT_trueREAL"].GetNbinsX()+2):
                for jbin in range(1,h_allMC_2D["nT_trueREAL"].GetNbinsX()+2):
                    print("hT = ",h_allMC_2D["nT_trueREAL"].GetBinContent(ibin,jbin))
                    print("hL = ",h_allMC_2D["nL_trueREAL"].GetBinContent(ibin,jbin))
                    #if h_allMC_2D["nT_trueREAL"].GetBinContent(ibin,jbin) > h_allMC_2D["nL_trueREAL"].GetBinContent(ibin,jbin):
                    #    print "----> OBS!! <-----"
                    #    h_allMC_2D["nT_trueREAL"].SetBinContent(ibin,jbin,h_allMC_2D["nT_trueREAL"].GetBinContent(ibin-1,jbin-1))
                    #    h_allMC_2D["nL_trueREAL"].SetBinContent(ibin,jbin,h_allMC_2D["nL_trueREAL"].GetBinContent(ibin-1,jbin-1))
                    #    print "hT (fixed) = ",h_allMC_2D["nT_trueREAL"].GetBinContent(ibin,jbin)
                    #    print "hL (fixed) = ",h_allMC_2D["nL_trueREAL"].GetBinContent(ibin,jbin)
            h_allMC_2D["eff_trueREAL"] = TEfficiency(h_allMC_2D["nT_trueREAL"],h_allMC_2D["nL_trueREAL"])# = h_allMC_2D["nL_trueREAL"].Clone(h_allMC_2D["nL_trueREAL"].GetName().replace("nT","eff"))
            h_allMC_2D["eff_trueREAL"].SetName((h_allMC_2D["nT_trueREAL"].GetName().replace("nT","eff")))
            #h_allMC_2D["eff_trueFAKE"].SetDirectory(0)
            h_allMC_2D["eff_trueFAKE"].Draw("colz text e")

        MMfile = TFile("MMinput.root","UPDATE")
        MMfile.cd()
        for key in list(h_allMC_2D.keys()):
            print("Name is %s" %h_allMC_2D[key].GetName())
            if "lep_mu_pT" in h_allMC_2D[key].GetName():
                h_allMC_2D[key].GetXaxis().SetTitle("mu")
                h_allMC_2D[key].GetYaxis().SetTitle("p_{T} [GeV]")
                h_allMC_2D[key].GetYaxis().SetRangeUser(0,50.)
            else:
                if "REAL" in h_allMC_2D[key].GetName():
                    h_allMC_2D[key].SetTitle(" ; p_{T} [GeV] ; |#eta| ; Efficiency")
                else:
                    h_allMC_2D[key].SetTitle(" ; p_{T} [GeV] ; |#eta| ; Fake Rate")
                #h_allMC_2D[key].GetYaxis().SetTitle("|#eta|")
                #canv[-1].SetLogx()
            #if "REAL" in h_allMC_2D[key].GetName(): h_allMC_2D[key].GetZaxis().SetTitle("Efficiency")
            #else: h_allMC_2D[key].GetZaxis().SetTitle("Fake Rate")
            #if "REAL" in h_allMC_2D[key].GetName(): h_allMC_2D[key].GetZaxis().SetRangeUser(0.65,1.0)
            #else: h_allMC_2D[key].GetZaxis().SetRangeUser(0.0,1.0)      
            
            h_allMC_2D[key].Write(h_allMC_2D[key].GetName()+postfix)
            print("Writing %s to MMinput.root" %(h_allMC_2D[key].GetName()+postfix))

            # if not "trueREAL" in key and not "trueFAKE" in key:
            #     h_eff_sub[key+"_subtREAL"] = doSubtraction2D(h_allMC_2D["nT_trueREAL"],h_allMC_2D["nL_trueREAL"],hist[nT_name][datakey],hist[nL_name][datakey],datakey)
            #     h_eff_sub[key+"_subtREAL"].SetName(h_eff_sub[key+"_subtREAL"].GetName()+"_subtREAL")
            #     h_eff_sub[key+"_subtFAKE"] = doSubtraction2D(h_allMC_2D["nT_trueFAKE"],h_allMC_2D["nL_trueFAKE"],hist[nT_name][datakey],hist[nL_name][datakey],datakey)
            #     h_eff_sub[key+"_subtFAKE"].SetName(h_eff_sub[key+"_subtFAKE"].GetName()+"_subtFAKE")
            #     canv.append(TCanvas("c_eff_%s_subtREAL"%h_allMC_2D["nT_trueREAL"].GetName().replace("nT","eff"),"c_eff_%s_subtREAL"%h_allMC_2D["nT_trueREAL"].GetName().replace("nT","eff"),1))
            #     h_eff_sub[key+"_subtREAL"].Draw("colz text e")
            #     canv.append(TCanvas("c_eff_%s_subtFAKE"%h_allMC_2D["nT_trueFAKE"].GetName().replace("nT","eff"),"c_eff_%s_subtFAKE"%h_allMC_2D["nT_trueFAKE"].GetName().replace("nT","eff"),1))
            #     h_eff_sub[key+"_subtFAKE"].Draw("colz text e")
            #     h_eff_sub[key+"_subtFAKE"].Write(h_eff_sub[key][typ].GetName()+postfix)
            #     h_eff_sub[key+"_subtREAL"].Write(h_eff_sub[key][typ].GetName()+postfix)
            #     print "Writing " + h_eff_sub[key][typ].GetName() + " to file MMinput.root"

        MMfile.Close()
            ## do subtraction

        # h_eff_sub = {}
        # h_eff_sub_pEff = {}
        # #for nhist in range(len(h_allMC_2D["nT"])):
        # canv.append(TCanvas("c_2d_%s"%h_allMC_2D["nT"].GetName(),"c_2d_%s"%h_allMC_2D["nT"].GetName(),1))
        # print("hopla")
        # # Avoid having more T than L (needed in some cases, mostly for Z+jets)
        # for ibin in range(1,h_allMC_2D["nT"].GetNbinsX()+2):
        #     for jbin in range(1,h_allMC_2D["nT"].GetNbinsY()+2):
        #         print("hT = ",h_allMC_2D["nT"].GetBinContent(ibin))
        #         print("hL = ",h_allMC_2D["nL"].GetBinContent(ibin))
        #         if h_allMC_2D["nT"].GetBinContent(ibin,jbin) > h_allMC_2D["nL"].GetBinContent(ibin,jbin):
        #             print("----> OBS!! <-----")
        #             h_allMC_2D["nT"].SetBinContent(ibin,jbin,h_allMC_2D["nT"].GetBinContent(ibin-1,jbin-1))
        #             h_allMC_2D["nL"].SetBinContent(ibin,jbin,h_allMC_2D["nL"].GetBinContent(ibin-1,jbin-1))
        #             print("hT (fixed) = ",h_allMC_2D["nT"].GetBinContent(ibin,jbin))
        #             print("hL (fixed) = ",h_allMC_2D["nL"].GetBinContent(ibin,jbin))
        # #h_allMC_2D["eff"] = h_allMC_2D["nT"].Clone(h_allMC_2D["nT"].GetName().replace("nT","eff"))
        # #h_allMC_2D["eff"].Divide(h_allMC_2D["nT"],h_allMC_2D["nL"],1.,1.,'b')
        # h_allMC_2D_pEff["eff"] = TEfficiency(h_allMC_2D["nT"],h_allMC_2D["nL"])
        # h_allMC_2D_pEff["eff"].SetStatisticOption(TEfficiency.kBBayesian)
        # h_allMC_2D_pEff["eff"].SetConfidenceLevel(0.68);
        # h_allMC_2D_pEff["eff"].SetBetaAlpha(1.)
        # h_allMC_2D_pEff["eff"].SetBetaBeta(1.) 
        # for typ in hist[h].keys():
        #     if not "nT" in h: continue
        #     h_allMC_2D[typ] = hist[h][typ].Clone(hist[h][typ].GetName().replace("nT","eff"))
        #     h_allMC_2D[typ].Divide(hist[h][typ],hist[h.replace("nT","nL")][typ],1.,1.,'b')
        # h_allMC_2D["eff"].Draw("colz text e")
        # canv[-1].SetRightMargin(0.14);
        # canv.append(TCanvas("c_2deff_trueMC_teff","c_2deff_trueMC_teff",1))
        # canv[-1].cd()
        # canv[-1].SetRightMargin(0.14);
        # for i in range(0,h_allMC_2D["eff"].GetNbinsX()+1):
        #     for j in range(0,h_allMC_2D["eff"].GetNbinsY()+1):
        #         ibin = h_allMC_2D["eff"].GetBin(i,j)
        #         h_allMC_2D["eff"].SetBinError(ibin,(h_allMC_2D_pEff["eff"].GetEfficiencyErrorLow(ibin)+h_allMC_2D_pEff["eff"].GetEfficiencyErrorUp(ibin))/2.0)
        # h_allMC_2D_pEff["eff"].Draw("colz text e")
        # if "lep_mu_pT" in hEff[typ].GetName():
        #     h_allMC_2D["eff"].GetXaxis().SetTitle("mu")
        #     h_allMC_2D["eff"].GetYaxis().SetTitle("p_{T} [GeV]")
        #     h_allMC_2D["eff"].GetYaxis().SetRangeUser(0,50.)
        # else:
        #     h_allMC_2D["eff"].GetXaxis().SetTitle("p_{T} [GeV]")
        #     h_allMC_2D["eff"].GetYaxis().SetTitle("|#eta|")
        #     canv[-1].SetLogx()
        # if "REAL" in h_allMC_2D["eff"].GetName(): h_allMC_2D["eff"].GetZaxis().SetTitle("Efficiency")
        # else: h_allMC_2D["eff"].GetZaxis().SetTitle("Fake Rate")
        # if "REAL" in h_allMC_2D["eff"].GetName(): h_allMC_2D["eff"].GetZaxis().SetRangeUser(0.65,1.0)
        # else: h_allMC_2D["eff"].GetZaxis().SetRangeUser(0.0,1.0)

        # MMfile = TFile("MMinput.root","UPDATE")
        # MMfile.cd()
        # h_allMC_2D["eff"].Write(h_allMC_2D["eff"].GetName()+postfix)
        # h_allMC_2D_pEff["eff"].Write(h_allMC_2D["eff"].GetName()+postfix)
        # print("ohlala--- > Writing " + h_allMC_2D["eff"].GetName() + " to file MMinput.root")
        # print("ohlala--- > Writing " + h_allMC_2D_pEff["eff"].GetName() + " to file MMinput.root")
        # MMfile.Close()

        # canv[-1].SetRightMargin(0.14);
            

            # key = h_allMC_2D["eff"][nhist].GetName()
            # if not key in h_eff_sub.keys():
            #     h_eff_sub[key] = {}
            #     h_eff_sub_pEff[key] = {}
            # for typ in hEff:
            #     if "_hist" in typ: continue
            #     if not "data" in typ: continue
            #     canv.append(TCanvas("c_%s_2D_%s_%s"%(key,hEff[typ].GetName(),typ),"c_%s_2D_%s_%s"%(key,hEff[typ].GetName(),typ),1))
            #     canv[-1].cd()


            #     canv[-1].SetRightMargin(0.14);
            #     if doSub:
            #         h_eff_sub[key][typ] = doSubtraction2D(h_allMC_2D["nT"][nhist],h_allMC_2D["nL"][nhist],hist[nT_name][typ],hist[nL_name][typ],typ)
            #     else:
            #         h_eff_sub[key][typ] = hist[nL_name][typ].Clone((hist[nL_name][typ].GetName()).replace("nL","eff")+postfix)
            #         h_eff_sub[key][typ].Divide(hist[nT_name][typ],hist[nL_name][typ],1.,1.,'b')
            #         h_eff_sub_pEff[key]["eff"] = TEfficiency(hist[nT_name][typ],hist[nL_name][typ])
            #         h_eff_sub_pEff[key]["eff"].SetStatisticOption(TEfficiency.kBBayesian)
            #         h_eff_sub_pEff[key]["eff"].SetConfidenceLevel(0.68);
            #         h_eff_sub_pEff[key]["eff"].SetBetaAlpha(1.)
            #         h_eff_sub_pEff[key]["eff"].SetBetaBeta(1.) 
            #         for i in range(0,h_eff_sub[key][typ].GetNbinsX()+1):
            #             for j in range(0,h_eff_sub[key][typ].GetNbinsY()+1):
            #                 ibin = h_eff_sub[key][typ].GetBin(i,j)
            #                 h_eff_sub[key][typ].SetBinError(ibin,(h_eff_sub_pEff[key]["eff"].GetEfficiencyErrorLow(ibin)+h_eff_sub_pEff[key]["eff"].GetEfficiencyErrorUp(ibin))/2.0)
            #     print "after", type(h_eff_sub[key][typ])

            #     if "lep_mu_pT" in hEff[typ].GetName():
            #         h_eff_sub[key][typ].GetXaxis().SetTitle("mu")
            #         h_eff_sub[key][typ].GetYaxis().SetTitle("p_{T} [GeV]")
            #         h_eff_sub[key][typ].GetYaxis().SetRangeUser(20,50.)
            #         h_eff_sub[key][typ].GetXaxis().SetRangeUser(0,75)
            #         h_eff_sub[key][typ].SetMarkerSize(1.6)
            #     else:
            #         h_eff_sub[key][typ].GetXaxis().SetTitle("p_{T}")
            #         h_eff_sub[key][typ].GetYaxis().SetTitle("|#eta|")
            #         canv[-1].SetLogx()
            #     h_eff_sub[key][typ].GetZaxis().SetTitle("Fake rate")
            #     if "REAL2L" in nL_name:
            #         h_eff_sub[key][typ].GetZaxis().SetRangeUser(0.65,1.0)
            #         h_eff_sub[key][typ].SetMarkerSize(0.75)
            #     else:
            #         h_eff_sub[key][typ].GetZaxis().SetRangeUser(0.0,1.0)
            #         h_eff_sub[key][typ].SetMarkerSize(1.0)
            #     h_eff_sub[key][typ].Draw("colz text e")
            #     MMfile = TFile("MMinput.root","UPDATE")
            #     MMfile.cd()
            #     h_eff_sub[key][typ].Write(h_eff_sub[key][typ].GetName()+postfix)
            #     print "Writing " + h_eff_sub[key][typ].GetName() + " to file MMinput.root"
            #     MMfile.Close()
    for c in canv:
       c.cd()
       c.Update()
       c.SaveAs("plots/%s_zoomed.png" %c.GetName())
else:
    stack_resid = {}
    stack_resid_UP = {}
    stack_resid_DOWN = {}
    hist_nofake = {}
    h_err_tot = {}
    for h in list(hist.keys()):
        f_num_of_events = open("num_of_events_%s.txt"%(h),"w")
        print("fake_hist_name",fake_hist_name)
        print("data_hist_name",data_hist_name)
        if "estfake" in h or fake_hist_name in h or data_hist_name in h: continue
        print("a")
        print(h)
        if not h in list(metadata.keys()):
            metadata[h] = {"xtit":"","ytit":"","label":""}
            fillMetadata(h)
        canv.append(TCanvas("c_%s_%s"%(h,datakey),"c_%s_%s"%(h,datakey),1))
        canv[-1].cd()
        doRatio = True
        if "TCF" in h: doRatio = False
        if doRatio:# 'NOCUTS' in hn and 
           
            pad1.append(TPad("pad1_%s"%h,"pad1_%s"%h,0.0,0.3,0.8,1.0,21))
            pad3.append(TPad("pad3_%s"%h,"pad3_%s"%h,0.8,0.0,1.0,1.0,kWhite,1,1))
            pad3[-1].SetLineColor(kWhite)
            pad1[-1].Draw();

            if not "CRTOP" in h:
                pad1[-1].SetLogy(1);
            if not estfake and not "eta" in h and doLOG: pad1[-1].SetLogx(1); #not estfake and 
            
            
            pad2.append(TPad("pad2_%s"%h,"pad2_%s"%h,0.0,0.0,0.8,0.3,22))
            pad2[-1].SetRightMargin(0.03);
            pad2[-1].Draw();
            pad2[-1].SetFillColor(kWhite);
            pad2[-1].SetTopMargin(0.05);
            pad2[-1].SetBottomMargin(0.3);
            if not estfake and not "eta" in h and doLOG: pad2[-1].SetLogx(1); #not estfake and 

            pad1[-1].SetFillColor(kWhite);
            pad1[-1].SetBottomMargin(0.01);
            pad1[-1].SetRightMargin(0.03);
            pad3[-1].Draw()
            pad3[-1].SetFillColor(kWhite)

            pad3[-1].cd()
            l = TLatex() 
            l.SetNDC();
            l.SetTextColor(kBlack);
            l.SetTextSize(0.1)
            l.DrawLatex(0.0658166,0.44916,"#int L = %.1f fb^{-1}" %(lumi[datakey]));

            # l = TLatex() 
            # l.SetNDC();
            # l.SetTextColor(kBlack);
            # l.SetTextSize(0.09)
            # l.DrawLatex(0.00358166,0.64916,"%s (%.1f fb^{-1})" %(datakey,lumi[datakey]));
            
            pad1[-1].cd() 
        #if not "eta" in h and not "TCF" in h and not "mt2" in h and not "met" in h and not "mll" in h:
        #  pad1.SetLogx()
        #    pad2.SetLogx()
        stack[h] = THStack();
        if doRatio: 
            leg[h] = TLegend(0.00358166,0.54916,0.984957,0.915966)
            leg[h].SetTextSize(0.09)
        else: 
            if "_M" in h and "FAKE2L21" in h:
                leg[h] = TLegend(0.123209,0.617647,0.272206,0.983193)
            else:
                leg[h] = TLegend(0.809456,0.613924,0.958453,0.981013)
            leg[h].SetTextSize(0.03)
            canv[-1].SetLogy()
            canv[-1].SetLogx()
        leg[h].SetBorderSize(0)
        leg[h].SetTextFont(42)
        leg[h].SetTextColor(1)
        leg[h].SetFillColor(0)
        leg[h].SetLineColor(0)
        bkg_sorted = []
        err = ctypes.c_double(0)#Double(0)
        err_up = ctypes.c_double(0) #Double(0)
        err_dw = ctypes.c_double(0) #Double(0)
        if doMC: bkg_sorted = getSortedList(evcount[h])
        doneType = []
        sum_intgrl = 0
        for typ in reversed(bkg_sorted):
            if typ in doneType: continue
            print("Adding %s to stack" %typ)
            hist[h][typ].SetFillColor(d_samp[typ]['f_color'])
            hist[h][typ].SetLineColor(d_samp[typ]['f_color'])
            leg[h].AddEntry(hist[h][typ],d_samp[typ]['leg'],"f")
            print("Adding %.2f events to stack from %s" %(hist[h][typ].Integral(),typ))
            stack[h].Add(hist[h][typ])
            if not h in list(hist_nofake.keys()):
                hist_nofake[h] = hist[h][typ].Clone("h_MC_nofake")
            else:
                hist_nofake[h].Add(hist[h][typ])
            doneType.append(typ)
            
            print(type(hist[h][typ]), h, typ)
            if type(hist[h][typ]) is TH2F:
                intgrl = hist[h][typ].IntegralAndError(2,2,2,2,err)#(0,hist[h][typ].GetNbinsX()+1,err)
            else:
                intgrl = hist[h][typ].IntegralAndError(2,2,err)#(0,hist[h][typ].GetNbinsX()+1,err)
            f_num_of_events.write("%10s %4.2f +/- %4.2f\n" %(typ,intgrl,float(err.value)))
            sum_intgrl += intgrl
        if not sum_intgrl: canv[-1].SetLogy(0)
        if estfake:
            hnm = h+"_estfake"
            tnm = datakey
            hist[hnm][datakey].SetFillColor(d_samp['fake']['f_color'])
            hist[hnm][datakey].SetLineColor(d_samp['fake']['f_color'])
            leg[h].AddEntry(hist[hnm][datakey],d_samp['fake']['leg'],"f")
            stack[h].Add(hist[hnm][datakey])
            intgrl = hist[hnm][datakey].IntegralAndError(2,2,err)#0,hist[hnm][datakey].GetNbinsX()+1,err)
            ingrl_up = hist[hnm+"_UP"][datakey].IntegralAndError(2,2,err_up)
            ingrl_dw = hist[hnm+"_DW"][datakey].IntegralAndError(2,2,err_dw)
            f_num_of_events.write("%10s %.2f +/- %.2d + %.2d - %.2d\n" %("Fakes",intgrl,err.value,err_up.value,err_dw.value))
            print("mT2 > 100: %.2f" %hist[hnm][datakey].Integral(hist[hnm][datakey].GetXaxis().FindBin(100),hist[hnm][datakey].GetNbinsX()+1))
            for ibin in range(1,hist[hnm][datakey].GetNbinsX()+1):
                print("[%i-%i] %.2f" %(hist[hnm][datakey].GetBinLowEdge(ibin),hist[hnm][datakey].GetBinLowEdge(ibin+1),hist[hnm][datakey].GetBinContent(ibin)))
        if doMC:
            #if not is2D: 
            try:
                stack[h].Draw("hist")
                stack[h].GetXaxis().SetMoreLogLabels(True)
            except:
                print("Failed")
            if not doRatio: 
                stack[h].GetXaxis().SetTitle(metadata[h]['xtit'])
                stack[h].GetXaxis().SetLabelSize(0.035)
                stack[h].GetXaxis().SetTitleSize(0.03)
                stack[h].GetXaxis().SetTitleOffset(1.8)
            stack[h].GetYaxis().SetTitle(metadata[h]['ytit'])
            stack[h].GetYaxis().SetTitleSize(0.06)
            stack[h].GetYaxis().SetTitleOffset(0.7)
            stack[h].GetXaxis().SetRangeUser(metadata[h]['xmin'],metadata[h]['xmax'])  
            hnm_data = ""
            if datakey in list(hist[h].keys()):
                if not data_hist_name == "=": hnm_data = data_hist_name
                else: hnm_data = h
            if "MMSS" in h and "CRTOP" in h:
                stack[h].SetMaximum(12)
                stack[h].SetMinimum(0)
            elif not "TCF" in h and hnm_data:
                stack[h].SetMaximum(stack[h].GetMaximum()*1.5 if stack[h].GetMaximum() > hist[hnm_data][datakey].GetMaximum() else hist[hnm_data][datakey].GetMaximum()*1.5)
                stack[h].SetMinimum(1)
        
        if doMC:
            hsumMC[h] = (stack[h].GetStack().Last()).Clone(h+"_all_"+addname)
            if type(hsumMC[h]) is TH2F:
                intgrl = hsumMC[h].IntegralAndError(2,2,2,2,err)#(0,hsumMC[h].GetNbinsX()+1,err)
            else:
                intgrl = hsumMC[h].IntegralAndError(2,2,err)#(0,hsumMC[h].GetNbinsX()+1,err)
        else:
            intgrl = 0.0
        f_num_of_events.write("-"*30)
        f_num_of_events.write("\n")
        f_num_of_events.write("%10s %.2f +/- %.2f\n" %("Total",intgrl,err.value))
        f_num_of_events.write("-"*30)
        f_num_of_events.write("\n")
        hnm_data = ""
        if datakey in list(hist[h].keys()):
            if not data_hist_name == "=": hnm_data = data_hist_name
            else: hnm_data = h
            print("Data histname is %s, original histname is %s"%(hnm_data,h))
            if not estfake:
                if "data15-16" in datakey or "alldata" in datakey: hist[hnm_data][datakey].SetMarkerStyle(20)
                if "data17" in datakey: hist[hnm_data][datakey].SetMarkerStyle(4)
                if "data18" in datakey: hist[hnm_data][datakey].SetMarkerStyle(30)
                hist[hnm_data][datakey].SetMarkerColor(d_samp[datakey]['f_color'])
                hist[hnm_data][datakey].SetLineColor(d_samp[datakey]['f_color'])
            else:
                hist[hnm_data][datakey].SetMarkerStyle(20)
                hist[hnm_data][datakey].SetMarkerColor(kBlack)
                hist[hnm_data][datakey].SetLineColor(kBlack)
                
            
            hist[hnm_data][datakey].Draw("same e")
            leg[h].AddEntry(hist[hnm_data][datakey],"%s"%(d_samp[datakey]['leg']),"lp")
            if type(hist[hnm_data][datakey]) is TH2F:
                intgrl = hist[hnm_data][datakey].IntegralAndError(2,2,2,2,err)#(0,hist[hnm_data][datakey].GetNbinsX()+1,err)
            else:
                intgrl = hist[hnm_data][datakey].IntegralAndError(2,2,err)#(0,hist[hnm_data][datakey].GetNbinsX()+1,err)
            f_num_of_events.write("%10s %.2f +/- %.2f" %(datakey,intgrl,err.value))

        if doMC and hnm_data:
            hsumMC[h] = (stack[h].GetStack().Last()).Clone(h+"_all_"+addname)
            tot_scalef = 0
            nbins = 0
            for i in range(1,hsumMC[h].GetNbinsX()+1):
                mc_binc   = hsumMC[h].GetBinContent(i)
                data_binc = hist[hnm_data][datakey].GetBinContent(i)
                print("mc_binc =",mc_binc)
                print("data_binc =",data_binc)
                if mc_binc:
                    print("scalef %i = %.4f" %(i,data_binc/mc_binc))
                    tot_scalef += data_binc/mc_binc
                    nbins += 1


            if nbins > 0: print("Total scalef is %.4f" %(tot_scalef/nbins))

                    

        if estfake and fake_hist_name+"_UP" in list(hist.keys()):

            print("Getting systematics")
            error_up_x = array.array('f')
            error_low_x = array.array('f')
            error_up_y = array.array('f')
            error_low_y = array.array('f')
            x = array.array('f')
            y = array.array('f')
            syst_prct = 0.20
            stat_prct = 0.05
            for ibin in range(1,hist_nofake[h].GetNbinsX()+1):

                nom = hist[hnm][datakey].GetBinContent(ibin)
                nom_fake = hist[fake_hist_name][datakey].GetBinContent(ibin)

                if not nom_fake: 
                    error_up_y.append(0.0)
                    error_low_y.append(0.0)
                    x.append(hsumMC[h].GetXaxis().GetBinCenter(ibin))
                    y.append(hsumMC[h].GetBinContent(ibin))
                    
                    error_up_x.append(hsumMC[h].GetXaxis().GetBinUpEdge(ibin)-hsumMC[h].GetXaxis().GetBinCenter(ibin))
                    error_low_x.append(hsumMC[h].GetXaxis().GetBinCenter(ibin)-hsumMC[h].GetXaxis().GetBinLowEdge(ibin))
                    continue

                #stat_err_y_up  = hist_nofake[h].GetBinError(ibin);
                #stat_err_y_low = hist_nofake[h].GetBinError(ibin);

                stat_err_y_up  = hsumMC[h].GetBinError(ibin);
                stat_err_y_low = hsumMC[h].GetBinError(ibin);

                if fake_hist_name+"_UP" in list(hist.keys()):
                    up_wgt = hist[fake_hist_name+"_UP"][datakey].GetBinContent(ibin)
                    dw_wgt = hist[fake_hist_name+"_DW"][datakey].GetBinContent(ibin)


                # MC syst uncertinty
                nom_mc = hsumMC[h].GetBinContent(ibin);
                mc_unc = 0.1

                # print "stat_err_y_up = ", stat_err_y_up
                # print "stat_err_y_dw = ", stat_err_y_low
                # print "((stat_err_y_up)/nom_mc) = ",((stat_err_y_up)/nom_mc)
                # print "((up_wgt-nom_fake)/nom_fake) = ",((up_wgt-nom_fake)/nom_fake)
                # print "up_wgt = ", up_wgt
                # print "dw_wgt = ", dw_wgt
                # print "nom_fake = ", nom_fake

                # print "Bin %i" %ibin
                # print "Nom MC = %.2f" %nom_mc
                
                #i# f nom_fake > 0:
                #     print "Fake stat = %.2f (sqrt(%.2f) = %.2f)" %(hist[fake_hist_name][datakey].GetBinError(ibin)/nom_fake,nom_fake,sqrt(nom_fake))
                # if hist_nofake[h].GetBinContent(ibin) > 0: 
                #     print "MC stat = %.2f (sqrt(%.2f) = %.2f)" %(hist_nofake[h].GetBinError(ibin)/hist_nofake[h].GetBinContent(ibin),hist_nofake[h].GetBinContent(ibin),sqrt(hist_nofake[h].GetBinContent(ibin)))
                # if hsumMC[h].GetBinContent(ibin) > 0:
                #     if hsumMC[h].GetBinContent(ibin): print "MC + Fake stat = %.2f (sqrt(%.2f) = %.2f)" %(hsumMC[h].GetBinError(ibin)/hsumMC[h].GetBinContent(ibin),hsumMC[h].GetBinContent(ibin),sqrt(hsumMC[h].GetBinContent(ibin)))
                # print "mc_unc = %.2f" %mc_unc

                print("-"*100)
                print("Bin %i" %ibin)
                print("\n")
                print("Stat syst UP = %.2f" %((stat_err_y_up)/nom_mc))
                print("Stat syst DW = %.2f" %((stat_err_y_low)/nom_mc))
                print("\n")
                print("Fake syst UP = %.2f" %((up_wgt-nom_fake)/nom_fake))
                print("Fake syst DW = %.2f" %((nom_fake-dw_wgt)/nom_fake))
                print("\n")
                print("Full UP unc: ", (sqrt(((up_wgt-nom_fake)/nom_fake)**2 + (stat_err_y_up/nom_mc)**2)))
                print("Full DW unc: ", (sqrt(((nom_fake-dw_wgt)/nom_fake)**2 + (stat_err_y_low/nom_mc)**2)))
                print("-"*100)

                #error_up_y.append(up_wgt)#sqrt( ((up_wgt-nom_fake)/nom_fake)**2 + stat_err_y_up**2))
                #error_up_y.append(nom_fake * sqrt(((up_wgt-nom_fake)/nom_fake)**2 + ((stat_err_y_up-nom_fake)/nom_fake)**2))
                #print "err up: ",(sqrt(((up_wgt-nom_fake)/nom_fake)**2 + (stat_err_y_up/nom_mc)**2))# + (mc_unc)**2))

                if 0:
                    error_up_y.append(nom_mc * (sqrt(((up_wgt-nom_fake)/nom_fake)**2 + (stat_err_y_up/nom_mc)**2)))# + (mc_unc)**2)))
                    dw_wgt = dw_wgt if dw_wgt > 0 else 0.0
                    error_low_y.append(nom_mc * (sqrt(((nom_fake-dw_wgt)/nom_fake)**2 + (stat_err_y_low/nom_mc)**2)))## + (mc_unc)**2)))
                
                
                    # if "EMSS" in  h:
                    #     error_low_y[-1] = error_low_y[-1]*0.1
                    #     error_up_y[-1] = error_up_y[-1]*0.1
                    # if "MMSS" in  h:
                    #     print "Reducing uncertainties!"
                    #     error_low_y[-1] = error_low_y[-1]*0.4
                    #     error_up_y[-1] = error_up_y[-1]*0.4

                    print("Total up = %.2f" %error_up_y[-1])
                    print("Total dw = %.2f" %error_low_y[-1])

                else:
                    error_up_y.append(nom_mc*sqrt(0.25*0.25 + (stat_err_y_up/nom_mc)**2))# + (nom_mc*0.25))# + (mc_unc)**2)))
                    dw_wgt = dw_wgt if dw_wgt > 0 else 0.0
                    error_low_y.append(nom_mc*sqrt(0.25*0.25 + (stat_err_y_up/nom_mc)**2))# - (nom_mc*0.25))## + (mc_unc)**2)))


                x.append(hsumMC[h].GetXaxis().GetBinCenter(ibin))
                y.append(hsumMC[h].GetBinContent(ibin))

                error_up_x.append(hsumMC[h].GetXaxis().GetBinUpEdge(ibin)-hsumMC[h].GetXaxis().GetBinCenter(ibin))
                error_low_x.append(hsumMC[h].GetXaxis().GetBinCenter(ibin)-hsumMC[h].GetXaxis().GetBinLowEdge(ibin))

                

            h_err_tot[h] = TGraphAsymmErrors(len(error_low_x),x,y,error_low_x,error_up_x,error_low_y,error_up_y)
            h_err_tot[h].SetLineStyle(kDashed);
            h_err_tot[h].SetLineColor(kGray+3);
            h_err_tot[h].SetFillColor(kGray+3);
            h_err_tot[h].SetFillStyle(3018)
            h_err_tot[h].SetMarkerSize(0)
            h_err_tot[h].Draw("e2same");
            #print "filled systematics"
        
                

            
        if doMC and hnm_data:   
            stack_resid[h] = TH1F()
            stack_resid_UP[h] = TH1F()
            stack_resid_DOWN[h] = TH1F()
            if doRatio:
                pad2[-1].cd().SetGridy(1);
                pad2[-1].cd().Update();
                pad2[-1].cd()
            rangeX = hist[hnm_data][datakey].GetXaxis().GetBinUpEdge(hist[hnm_data][datakey].GetNbinsX())
            if estfake and h in list(h_err_tot.keys()): stack_resid[h],stack_resid_UP[h],stack_resid_DOWN[h] = MakeResidPlot(hsumMC[h],hist[hnm_data][datakey],stack_resid[h],stack_resid_UP[h],stack_resid_DOWN[h],rangeX,h_err_tot[h],hist[fake_hist_name][datakey])
            else: stack_resid[h],stack_resid_UP[h],stack_resid_DOWN[h] = MakeResidPlot(hsumMC[h],hist[hnm_data][datakey],stack_resid[h],stack_resid_UP[h],stack_resid_DOWN[h],rangeX,0)
            if not "FAKE2L21" in h:
                stack_resid[h].GetYaxis().SetRangeUser(0.0,2.0)
            stack_resid[h].GetYaxis().SetNdivisions(7,4,0)
            stack_resid[h].GetXaxis().SetMoreLogLabels(True)
            if doRatio:
                stack_resid[h].Draw("pe0");
            stack_resid[h].GetXaxis().SetTitleOffset(1.04)
            stack_resid[h].GetXaxis().SetTitle(metadata[h]['xtit'])
            stack_resid[h].GetYaxis().SetTitle("Data/MC")
            stack_resid[h].GetYaxis().SetTitleOffset(0.35)
            stack_resid[h].GetYaxis().SetTitleSize(0.12)
            stack_resid[h].GetXaxis().SetRangeUser(metadata[h]['xmin'],metadata[h]['xmax'])   
            if doRatio:
                stack_resid_UP[h].Draw("hist][ same");
                stack_resid_DOWN[h].Draw("hist][ same");
                stack_resid[h].Draw("pe0 same");
                stack_resid[h].Draw("axis same");
                stack_resid[h].Draw("axiG same");

            if doRatio: pad3[-1].cd()
            else: gPad.cd()
            leg[h].Draw()
            if "trueREAL" in h or "trueLight" in h or "truePromptPhotonConversion" in h or "REAL2L0" in h or "REAL2L1" in h or "ElectronFromMuon" in h:
                if not "ALL" in subBkgs:
                    if "nT" in h: 
                        h_allMC["nT"] = TH1F()
                        h_allMC_xsecup["nT"] = TH1F()
                        h_allMC_xsecdw["nT"] = TH1F()
                        hname = "nT"
                        nT_name = h
                    if "nL" in h: 
                        h_allMC["nL"] = TH1F()
                        h_allMC_xsecup["nL"] = TH1F()
                        h_allMC_xsecdw["nL"] = TH1F()
                        hname = "nL"
                        nL_name = h
                    ntyp = 0
                    for typ in list(hist[h].keys()):
                        if "data" in typ: continue
                        if ntyp == 0 and (typ in subBkgs or "ALL" in subBkgs):
                            print("ntyp = %i and subtracting %s" %(ntyp,typ))
                            h_allMC[hname] = hist[h][typ].Clone(hist[h][typ].GetName().replace(typ,addname))
                            h_allMC_xsecup[hname] = hist[h][typ].Clone(hist[h][typ].GetName().replace(typ,addname)+"_xsecup")
                            h_allMC_xsecdw[hname] = hist[h][typ].Clone(hist[h][typ].GetName().replace(typ,addname)+"_xsecdw")
                            ntyp += 1
                            continue
                        if typ in subBkgs or "ALL" in subBkgs: 
                            print("ntyp = %i and subtracting %s" %(ntyp,typ))
                            h_allMC[hname].Add(hist[h][typ],1)
                            h_allMC_xsecup[hname].Add(hist[h][typ],1)
                            h_allMC_xsecdw[hname].Add(hist[h][typ],1)
                            ntyp += 1
                else:
                    print("Subtracting ALL")
                    if "nT" in h:
                        h_allMC["nT"] = (stack[h].GetStack().Last()).Clone(h+"_all_"+addname) 
                        h_allMC_xsecup["nT"] = (stack[h].GetStack().Last()).Clone(h+"_all_"+addname+"_xsecup") 
                        h_allMC_xsecup["nT"].Scale(1.1)
                        h_allMC_xsecdw["nT"] = (stack[h].GetStack().Last()).Clone(h+"_all_"+addname+"_xsecdw") 
                        h_allMC_xsecdw["nT"].Scale(0.9)
                        nT_name = h
                    elif "nL" in h:
                        h_allMC["nL"] = (stack[h].GetStack().Last()).Clone(h+"_all_"+addname) 
                        h_allMC_xsecup["nL"] = (stack[h].GetStack().Last()).Clone(h+"_all_"+addname+"_xsecup") 
                        h_allMC_xsecup["nL"].Scale(1.1)
                        h_allMC_xsecdw["nL"] = (stack[h].GetStack().Last()).Clone(h+"_all_"+addname+"_xsecdw") 
                        h_allMC_xsecdw["nL"].Scale(0.9)
                        nL_name = h

        # MMfile = TFile("MMinput.root","UPDATE")
        # MMfile.cd()
        # if "nT" in h_allMC_xsecup.keys():
        #     h_allMC_xsecup["nT"].Write(h_allMC_xsecup["nT"].GetName()+postfix)
        #     print("Writing %s to file"%h_allMC_xsecup["nT"].GetName()+postfix)
        # if "nL" in h_allMC_xsecup.keys():
        #     h_allMC_xsecup["nL"].Write(h_allMC_xsecup["nL"].GetName()+postfix)
        #     print("Writing %s to file"%h_allMC_xsecup["nL"].GetName()+postfix)
        # if "nT" in h_allMC_xsecdw.keys():
        #     h_allMC_xsecdw["nT"].Write(h_allMC_xsecdw["nT"].GetName()+postfix)
        #     print("Writing %s to file"%h_allMC_xsecdw["nT"].GetName()+postfix)
        # if "nL" in h_allMC_xsecdw.keys():
        #     h_allMC_xsecdw["nL"].Write(h_allMC_xsecdw["nL"].GetName()+postfix)
        #     print("Writing %s to file"%h_allMC_xsecdw["nL"].GetName()+postfix)
       
        # MMfile.Close()
        f_num_of_events.close()
        #canv[-1].Update()
        #canv[-1].SaveAs("plots/%s_zoomed.png" %canv[-1].GetName())
    for h in list(h_err_tot.keys()):
        newfile = open("cutflow_"+h+".txt","w")
        err = ctypes.c_double(0)#Double(0)
        print("Cutflow for histogram: %s" %h)
        print("-"*35)
        newfile.write("Cutflow for histogram: %s\n" %h)
        newfile.write("-"*35)
        newfile.write("\n")
        nbkg = 0
        y_up = 0
        y_dw = 0
        liboffstr = ""
        for i in range(h_err_tot[h].GetN()):
            nominal = hist[h+"_estfake"][datakey].GetBinContent(i+1)
            #print "up", h_err_tot[h.replace("trueREAL","trueFAKE")].GetErrorYhigh(i)
            #print "nm", nominal
            #print "dw", h_err_tot[h.replace("trueREAL","trueFAKE")].GetErrorYlow(i)
            y_up += (h_err_tot[h.replace("trueREAL","trueFAKE")].GetErrorYhigh(i))
            y_dw += (h_err_tot[h.replace("trueREAL","trueFAKE")].GetErrorYlow(i))
        for bkg in bkg_sorted:
            if nbkg == len(bkg_sorted)-1: break
            integr = hist[h][bkg].IntegralAndError(0,hist[h][bkg].GetNbinsX()+1,err)
            liboffstr += "%s|%.2f|+/-|%.2f\n" %(bkg,integr,err.value)
            yieldstr = "%.2f | +/- | %.2f" %(integr,err.value)
            print("%10s | %20s" %(bkg,yieldstr))
            newfile.write("%10s | %20s\n" %(bkg,yieldstr))
            nbkg += 1
        if estfake:
            integr = hist[h+"_estfake"][datakey].IntegralAndError(0,hist[h+"_estfake"][datakey].GetNbinsX()+1,err)
            yieldstr = "%.2f | +/- | %.2f" %(integr,(y_up+y_dw)/2.)
            liboffstr += "%s|%.2f|+/-|%.2f\n" %("FNP",integr,(y_up+y_dw)/2.)
            print("%10s | %20s" %("FNP",yieldstr))
            newfile.write("%10s | %20s\n" %("FNP",yieldstr))
        print("-"*35)
        newfile.write("-"*35)
        newfile.write("\n")
        if doMC:
            integr = hsumMC[h].IntegralAndError(0,hsumMC[h].GetNbinsX()+1,err)
            yieldstr = "%.2f | +/- | %.2f" %(integr,sqrt(err.value*err.value+((y_up+y_dw)/2.)*((y_up+y_dw)/2.)))
            liboffstr += "%s|%.2f|+/-|%.2f\n" %("TOTAL",integr,sqrt(err.value*err.value+((y_up+y_dw)/2.)*((y_up+y_dw)/2.)))
            print("%-10s | %20s" %("Total",yieldstr))
            print("-"*35)
            newfile.write("%-10s | %20s\n" %("Total",yieldstr))
            newfile.write("-"*35)
            newfile.write("\n")
        #print hnm_data
        if hnm_data in list(hist.keys()):
            integr = hist[hnm_data][datakey].IntegralAndError(0,hist[hnm_data][datakey].GetNbinsX()+1,err)
            yieldstr = "%.2f | +/- | %.2f" %(integr,err.value)
            print("%-10s | %20s" %("Observed",yieldstr))
            newfile.write("%-10s | %20s\n" %("Observed",yieldstr))
            liboffstr += "%s|%.2f|\n" %("Observed",integr)
        print(liboffstr)
    
    if "nT" in list(h_allMC.keys()): 
        for ibin in range(1,h_allMC["nT"].GetNbinsX()+1):
            print("%i : %.2f" %(ibin,h_allMC["nT"].GetBinContent(ibin)))
    if not doRatio:
        leg[h].Draw()

    
    if datakey and len(list(hEff.keys()))>0:# and datakey in hEff:
        if not datakey in hEff and len(list(hEff.keys()))>0:
            canv.append(TCanvas("c_%s"%(hEff[list(hEff.keys())[0]].GetName()),"c_%s"%(hEff[list(hEff.keys())[0]].GetName()),1))
            dkey = list(hEff.keys())[0]
        elif datakey in hEff:
            canv.append(TCanvas("c_%s"%(hEff[datakey].GetName()),"c_%s"%(hEff[datakey].GetName()),1))
            dkey = datakey
        canv[-1].cd()
        if not "eta" in hEff[dkey].GetName() and not "trig" in hEff[dkey].GetName():
            if doLOG: canv[-1].SetLogx()
            leg["eff"] = TLegend(0.530086,0.254202,0.826648,0.491597)
        if "FAKE2L23" in hEff[dkey].GetName() or "FAKE2L27" in hEff[dkey].GetName() or "FAKE2L24" in hEff[dkey].GetName() or "FAKE2L25" in hEff[dkey].GetName() or "_truePromptPhotonConversion" in hEff[dkey].GetName() or "ElectronFromMuon" in hEff[dkey].GetName(): 
            if "_pT_" in hEff[dkey].GetName():
                leg["eff"] = TLegend(0.574499,0.181435,0.869628,0.320675)
            elif "_eta_" in hEff[dkey].GetName():
                leg["eff"] = TLegend(0.133238,0.74346,0.428367,0.880591)
        elif "FAKE2L20" in hEff[dkey].GetName(): leg["eff"] = TLegend(0.180516,0.706751,0.475645,0.843882)
        elif "FAKE2L21" in hEff[dkey].GetName() or "_trueLight" in hEff[dkey].GetName(): leg["eff"] = TLegend(0.180516,0.704726,0.475645,0.869198)
        elif "FAKE2L03" in hEff[dkey].GetName(): leg["eff"] = TLegend(0.19341,0.582278,0.489971,0.721519)
        else: leg["eff"] = TLegend(0.530086,0.254202,0.826648,0.491597)
        leg["eff"].SetBorderSize(0)
        leg["eff"].SetTextFont(42)
        leg["eff"].SetTextSize(0.04)
        leg["eff"].SetTextColor(1)
        leg["eff"].SetFillColor(0)
        leg["eff"].SetLineColor(0)
        neff = 0
        for typ in list(hEff.keys()):
            if ("_trueLight" in hEff[typ].GetName() or "_truePromptPhotonConversion" in hEff[typ].GetName() or "ElectronFromMuon" in hEff[typ].GetName()) and "data" in typ: continue
            if doSub and "data" in typ: continue
            
            error_up_x = array.array('f')
            error_low_x = array.array('f')
            error_up_y = array.array('f')
            error_low_y = array.array('f')
            x = array.array('f')
            y = array.array('f')
            
            if not hEff[typ].GetName() in list(metadata.keys()):
                metadata[hEff[typ].GetName()] = {"xtit":"","ytit":"","label":""}
                fillMetadata(hEff[typ].GetName())
            leg["eff"].SetHeader(metadata[hEff[typ].GetName()]["label"])
            if "MM" in hEff[typ].GetName() and "FAKE2L20" in hEff[typ].GetName() and not "data" in typ: continue
            if not "data" in typ:
                #if not ("_trueLight" in hEff[typ].GetName() or "_truePromptPhotonConversion" in hEff[typ].GetName()) and not "REAL2L" in hEff[typ].GetName(): continue
                print("------------------>",typ)
                #if "lowMassDY" in typ: continue
                #if "Vgamma" in typ: continue
                print(hEff[typ].GetName())
                if "_trig_" in hEff[typ].GetName():
                    continue
                #elif "REAL2L0" in hEff[typ].GetName() and not typ in bkg_sorted[:4]: 
                #    print "Not Plotting ", typ
                #    continue
                #if "REAL2L" in hEff[typ].GetName() and (not (typ in bkg_sorted[:3] or typ in ["ttbar"]) or typ in ["Wjets"] or typ in ["Vgamma"]): continue #not typ == "Wjets" and 
                elif ("_trueLight" in hEff[typ].GetName() or "_truePromptPhotonConversion" in hEff[typ].GetName() or "ElectronFromMuon" in hEff[typ].GetName()) and not typ in bkg_sorted: continue #not typ == "Wjets" and 
                elif ("2L21" in hEff[typ].GetName() or "2L4" in hEff[typ].GetName()) and not typ in bkg_sorted[:1]: continue
                elif ("2L23" in hEff[typ].GetName() or "2L27" in hEff[typ].GetName() or "2L24" in hEff[typ].GetName() or "2L25" in hEff[typ].GetName()) and not typ in bkg_sorted[:1]: continue
                elif "2L20" in hEff[typ].GetName() and not typ in ['ttbar']: continue
                elif "2L02" in hEff[typ].GetName() and not "REAL2L" in hEff[typ].GetName() and not "_trueLight" in hEff[typ].GetName() and not typ in ['ttbar','diboson','Vgamma','singleTop']: continue
                #elif "trueLight" in hEff[typ].GetName() and typ not in ['ttbar_nominal']: continue
                print("hei")
                hEff[typ].SetMarkerStyle(21+neff)
            else:
                hEff[typ].SetMarkerStyle(20)
            print("Settting color for ",typ)
            hEff[typ].SetMarkerColor(d_samp[typ]['f_color'])
            hEff[typ].SetLineColor(d_samp[typ]['f_color'])
            MMfile = TFile("MMinput.root","UPDATE")
            MMfile.cd()
            hEff[typ].Write(hEff[typ].GetName()+postfix)
            hEff[dkey].Write(hEff[dkey].GetName()+postfix)
            print("Writing " + hEff[typ].GetName() + " to file MMinput.root")
            print("Writing " + hEff[dkey].GetName() + " to file MMinput.root")
            MMfile.Close()
            if neff == 0:
                hEff[typ].Draw("e")
                try:
                    hEff[typ].GetXaxis().SetMoreLogLabels(True)
                except:
                    print("no")
                hEff[typ].GetXaxis().SetTitle(metadata[hEff[typ].GetName()]['xtit'])
                hEff[typ].GetYaxis().SetTitle(metadata[hEff[typ].GetName()]['ytit'])
                if "REAL2L" in hEff[typ].GetName():
                    if "EE" in hEff[typ].GetName(): hEff[typ].GetYaxis().SetRangeUser(0.4,1.01)
                    if "MM" in hEff[typ].GetName(): hEff[typ].GetYaxis().SetRangeUser(0.4,1.01)
                elif "FAKE2L23" in hEff[typ].GetName() or "FAKE2L24" in hEff[typ].GetName() or "FAKE2L25" in hEff[typ].GetName(): hEff[typ].GetYaxis().SetRangeUser(0.0,1.0)
                elif "FAKE2L2" in hEff[typ].GetName(): hEff[typ].GetYaxis().SetRangeUser(0.0,1.0)
                elif "FAKE2L03" in hEff[typ].GetName(): hEff[typ].GetYaxis().SetRangeUser(0.0,0.8)
                elif "FAKE2L01" in hEff[typ].GetName(): hEff[typ].GetYaxis().SetRangeUser(0.0,0.8)
                elif "_true" in hEff[typ].GetName(): hEff[typ].GetYaxis().SetRangeUser(0.0,0.8)
            else:
                hEff[typ].Draw("same e")
            for ibin in range(1,hEff_xsecup[typ].GetNbinsX()+1):
                x.append(hEff[typ].GetXaxis().GetBinCenter(ibin))
                y.append(hEff[typ].GetBinContent(ibin))
            
                error_up_x.append(hEff_xsecup[typ].GetXaxis().GetBinUpEdge(ibin)-hEff_xsecup[typ].GetXaxis().GetBinCenter(ibin))
                error_low_x.append(hEff_xsecdw[typ].GetXaxis().GetBinCenter(ibin)-hEff_xsecdw[typ].GetXaxis().GetBinLowEdge(ibin))

                error_low_y.append(hEff_xsecup[typ].GetBinContent(ibin))
                error_up_y.append(hEff_xsecdw[typ].GetBinContent(ibin))

            # h_eff_err_xsec[typ] = TGraphAsymmErrors(len(error_low_x),x,y,error_low_x,error_up_x,error_low_y,error_up_y)
            # h_eff_err_xsec[typ].SetLineStyle(kDashed);
            # h_eff_err_xsec[typ].SetLineColor(d_samp[typ]['f_color']);
            # #h_eff_err_xsec[typ].SetFillColor((d_samp[typ]['f_color']));
            # #h_eff_err_xsec[typ].SetFillStyle(3018)
            # h_eff_err_xsec[typ].SetMarkerSize(0)
            # h_eff_err_xsec[typ].Draw("e1 same");
            
            leg["eff"].AddEntry(hEff[typ],d_samp[typ]['leg'],"pl")
            neff += 1
        print("Trying ",h)
        
        for h in list(stack.keys()):
            if not ("LightFla" in hEff[dkey].GetName() or "CO" in hEff[dkey].GetName() or "REAL2L0" in hEff[dkey].GetName() or "REAL2L1" in hEff[dkey].GetName()) or not doMC: continue# and (("_REAL" in h and not "trueFAKE" in h) or ("_FAKE" in h and not "trueREAL" in h)):
            if "nT" in h: continue
            print("Looking at ",h.replace("nL","nT")+"_all_"+addname)
            
            #if not "nT" in h_allMC.keys():
            h_allMC["nT"] = (stack[h.replace("nL","nT")].GetStack().Last()).Clone(h.replace("nL","nT")+"_all_"+addname) 
            #if not "nL" in h_allMC.keys():
            h_allMC["nL"] = (stack[h].GetStack().Last()).Clone(h+"_all_"+addname) 
                
            if not "MC" in list(hEff.keys()):
                hEff["MC"] = []
            hEff["MC"].append(getEff(h_allMC["nT"],h_allMC["nL"]))
            hEff["MC"][-1].SetName(h_allMC["nT"].GetName().replace("nT","eff"))
            canv.append(TCanvas("c_%s"%(hEff["MC"][-1].GetName()),"c_%s"%(hEff["MC"][-1].GetName()),1))
            if doLOG: canv[-1].SetLogx()
            if "MC" in list(hEff_xsecup.keys()):
                hEff_xsecup["MC"].append(getEff(h_allMC_xsecup["nT"],h_allMC_xsecup["nL"]))
                hEff_xsecup["MC"][-1].SetName(h_allMC_xsecup["nT"].GetName().replace("nT","eff"))
                hEff_xsecdw["MC"].append(getEff(h_allMC_xsecdw["nT"],h_allMC_xsecdw["nL"]))
                hEff_xsecdw["MC"][-1].SetName(h_allMC_xsecdw["nT"].GetName().replace("nT","eff"))
            hEff["MC"][-1].SetMarkerStyle(3)
            hEff["MC"][-1].SetMarkerColor(kMagenta)
            hEff["MC"][-1].SetLineColor(kMagenta)
            #if not "trueLight" in h:
            leg["eff"].AddEntry(hEff["MC"][-1],"Total MC","pl")
            hEff["MC"][-1].Draw("same e")
            MMfile = TFile("MMinput.root","UPDATE")
            MMfile.cd()
            hEff["MC"][-1].Write(hEff["MC"][-1].GetName()+postfix)
            #hEff[dkey].Write(hEff[dkey].GetName()+postfix)
            hEff[list(hEff.keys())[0]].Write(hEff[list(hEff.keys())[0]].GetName()+postfix)
            if "MC" in list(hEff_xsecup.keys()):
                hEff_xsecup["MC"][-1].Write(hEff_xsecup["MC"][-1].GetName()+postfix)
                hEff_xsecup[dkey].Write(hEff_xsecup[dkey].GetName()+postfix)
                hEff_xsecdw["MC"][-1].Write(hEff_xsecdw["MC"][-1].GetName()+postfix)
                hEff_xsecdw[dkey].Write(hEff_xsecdw[dkey].GetName()+postfix)
            h_allMC["nT"].Write(h_allMC["nT"].GetName()+postfix)
            h_allMC["nL"].Write(h_allMC["nL"].GetName()+postfix)
            if "nT" in list(h_allMC_xsecup.keys()):
                h_allMC_xsecup["nT"].Write(h_allMC_xsecup["nT"].GetName()+postfix)
                h_allMC_xsecup["nL"].Write(h_allMC_xsecup["nL"].GetName()+postfix)
                h_allMC_xsecdw["nT"].Write(h_allMC_xsecdw["nT"].GetName()+postfix)
                h_allMC_xsecdw["nL"].Write(h_allMC_xsecdw["nL"].GetName()+postfix)
            print("Writing " +  hEff["MC"][-1].GetName() + " to file MMinput.root")
            #print "Writing " +   hEff[dkey].GetName() + " to file MMinput.root"
            MMfile.Close()

        if not doSub:
            leg["eff"].Draw()
    #if not doSub and not is2D:
    #    for c in canv:
    #        c.Update()
    #        c.SaveAs("plots/%s.png" %c.GetName())
    # if "nT" in h_allMC.keys():
    #     print "hei"
    #     nbins = 0
    #     tot_scalef = 0
    #     for i in range(1,h_allMC["nT"].GetNbinsX()+1):

    #         mc_binc   = h_allMC["nT"].GetBinContent(i)
    #         data_binc = hist[nT_name][dkey].GetBinContent(i)

    #         if mc_binc:
    #             print "scalef %i = %.4f" %(i,data_binc/mc_binc)
    #             tot_scalef += data_binc/mc_binc
    #             nbins += 1

                
    #     if nbins: print "Total scalef is %.4f" %(tot_scalef/nbins)
        
        for typ in hEff:
            error_up_x = array.array('f')
            error_low_x = array.array('f')
            error_up_y = array.array('f')
            error_low_y = array.array('f')
            x = array.array('f')
            y = array.array('f')
            
            #if not "trueREAL" in hEff: continue
            print("typ ",typ)
            if not datakey in typ: continue
            if not hEff[typ].GetName() in list(metadata.keys()):
                metadata[hEff[typ].GetName()] = {"xtit":"","ytit":"","label":""}
                fillMetadata(hEff[typ].GetName())
            
            if doSub:
                print("Now, subtracting ", h_allMC["nT"].GetName())
                print("nL_name = ",nL_name)
                h_eff_sub[typ], h_nT_sub[typ], h_nL_sub[typ] = doSubtraction(h_allMC["nT"],h_allMC["nL"],hist[nT_name][typ],hist[nL_name][typ],typ)
                h_eff_sub_xsecup[typ], h_nT_sub_xsecup[typ], h_nL_sub_xsecup[typ] = doSubtraction(h_allMC_xsecup["nT"],h_allMC_xsecup["nL"],hist[nT_name][typ],hist[nL_name][typ],typ)
                h_eff_sub_xsecup[typ].SetName(h_eff_sub_xsecup[typ].GetName()+"_xsecup")
                h_eff_sub_xsecdw[typ], h_nT_sub_xsecdw[typ], h_nL_sub_xsecdw[typ] = doSubtraction(h_allMC_xsecdw["nT"],h_allMC_xsecdw["nL"],hist[nT_name][typ],hist[nL_name][typ],typ)
                h_eff_sub_xsecdw[typ].SetName(h_eff_sub_xsecdw[typ].GetName()+"_xsecdw")
                h_eff_sub[typ].SetMarkerStyle(20)
                h_eff_sub[typ].SetMarkerColor(d_samp[typ]['f_color'])
                h_eff_sub[typ].SetLineColor(d_samp[typ]['f_color'])

                h_nT_sub_xsecup[typ].SetName(h_nT_sub_xsecup[typ].GetName()+"_xsecup")
                h_nL_sub_xsecup[typ].SetName(h_nL_sub_xsecup[typ].GetName()+"_xsecup")
                h_nT_sub_xsecdw[typ].SetName(h_nT_sub_xsecdw[typ].GetName()+"_xsecdw")
                h_nL_sub_xsecdw[typ].SetName(h_nL_sub_xsecdw[typ].GetName()+"_xsecdw")


                for ibin in range(1,h_eff_sub_xsecup[typ].GetNbinsX()+1):

                    print("bin ",ibin)

                    print("h_eff_sub_xsecup[%s] = %.2f" %(typ,h_eff_sub_xsecup[typ].GetBinContent(ibin)))
                    print("h_eff_sub_xsecdw[%s] = %.2f" %(typ,h_eff_sub_xsecdw[typ].GetBinContent(ibin)))
                    print("h_eff_sub[%s] = %.2f" %(typ,h_eff_sub[typ].GetBinContent(ibin)))

                    x.append(h_eff_sub[typ].GetXaxis().GetBinCenter(ibin))
                    y.append(h_eff_sub[typ].GetBinContent(ibin))

                    error_up_x.append(h_eff_sub_xsecup[typ].GetXaxis().GetBinUpEdge(ibin)-h_eff_sub_xsecup[typ].GetXaxis().GetBinCenter(ibin))
                    error_low_x.append(h_eff_sub_xsecdw[typ].GetXaxis().GetBinCenter(ibin)-h_eff_sub_xsecdw[typ].GetXaxis().GetBinLowEdge(ibin))

                    print("x = ",x[-1])
                    print("error_up_x = ",error_up_x[-1])
                    print("error_low_x = ",error_low_x[-1])

                    xsecdw = h_eff_sub_xsecdw[typ].GetBinContent(ibin)
                    if h_eff_sub_xsecdw[typ].GetBinContent(ibin) < h_eff_sub[typ].GetBinContent(ibin):
                        xsecdw = 0

                    xsecup = h_eff_sub_xsecup[typ].GetBinContent(ibin)
                    if h_eff_sub_xsecup[typ].GetBinContent(ibin) > h_eff_sub[typ].GetBinContent(ibin):
                        xsecup = 0

                    error_up_y.append(xsecdw-h_eff_sub[typ].GetBinContent(ibin))
                    error_low_y.append(h_eff_sub[typ].GetBinContent(ibin)-xsecup)

                    print("y = ",y[-1])
                    print("error_up_y = ",error_up_y[-1])
                    print("error_low_y = ",error_low_y[-1])

                h_eff_err_xsec[h] = TGraphAsymmErrors(len(error_low_x),x,y,error_low_x,error_up_x,error_low_y,error_up_y)


                if neff == 0:
                    h_eff_sub[typ].Draw("e")
                    h_eff_sub[typ].GetYaxis().SetRangeUser(0.0,1.0)
                    h_eff_sub[typ].GetXaxis().SetTitle(metadata[hEff[typ].GetName()]['xtit'])
                    h_eff_sub[typ].GetYaxis().SetTitle(metadata[hEff[typ].GetName()]['ytit'])
                else:
                    h_eff_sub[typ].Draw("same e")

                h_eff_err_xsec[h].SetLineStyle(kDashed);
                h_eff_err_xsec[h].SetLineColor(kGray);
                h_eff_err_xsec[h].SetFillColor(kGray);
                h_eff_err_xsec[h].SetFillStyle(3018)
                h_eff_err_xsec[h].SetMarkerSize(0)
                h_eff_err_xsec[h].Draw("e2same");
                leg["eff"].AddEntry(h_eff_sub[typ],d_samp[typ]['leg'],"pl")
                leg["eff"].AddEntry(h_eff_err_xsec[h],"MC subtr. unc.","f")
                neff += 1
                MMfile = TFile("MMinput.root","UPDATE")
                MMfile.cd()

                h_eff_sub[typ].Write(h_eff_sub[typ].GetName()+postfix)
                h_eff_sub_xsecup[typ].Write(h_eff_sub_xsecup[typ].GetName()+postfix)
                h_eff_sub_xsecdw[typ].Write(h_eff_sub_xsecdw[typ].GetName()+postfix)

                h_nT_sub[typ].Write(h_nT_sub[typ].GetName()+postfix)
                h_nT_sub_xsecup[typ].Write(h_nT_sub_xsecup[typ].GetName()+postfix)
                h_nT_sub_xsecdw[typ].Write(h_nT_sub_xsecdw[typ].GetName()+postfix)
                
                h_nL_sub[typ].Write(h_nL_sub[typ].GetName()+postfix)
                h_nL_sub_xsecup[typ].Write(h_nL_sub_xsecup[typ].GetName()+postfix)
                h_nL_sub_xsecdw[typ].Write(h_nL_sub_xsecdw[typ].GetName()+postfix)

                print("Writing " +   h_eff_sub[typ].GetName() + " to file MMinput.root")
                print("Writing " +   h_eff_sub_xsecup[typ].GetName() + " to file MMinput.root")
                print("Writing " +   h_eff_sub_xsecdw[typ].GetName() + " to file MMinput.root")
                MMfile.Close()

        MMfile = TFile("MMinput.root","UPDATE")
        MMfile.cd()
        for key in list(hist.keys()):
            for typ in list(hist[key].keys()):
                hist[key][typ].Write(hist[key][typ].GetName()+postfix)
                print("---> Writing " + hist[key][typ].GetName()+postfix + " to file MMinput.root")
        MMfile.Close()
    


        leg["eff"].Draw()
    #canv[-1].Update()
    #canv[-1].SaveAs("plots/%s_zoomed.png" %(canv[-1].GetName()))
    
for c in canv:
    try:
        c.Update()
        print(c.GetName())
        c.SaveAs("plots/%s_zoomed.png" %c.GetName())
    except:
        print(("Problems saving %s"%c.GetName()))
        continue
#    c.SaveAs("plots/%s.pdf" %c.GetName())
    



if 0:
    for key in list(hist.keys()):
        if '_estfake_DW' in key: dwkey = key
        elif '_estfake_UP' in key: upkey = key
        elif '_estfake' in key: nomkey = key
        elif '_trueFAKE' in key and not '_estfake' in key: dtkey = key
    #print "Bin ## : %10s %10s %10s %10s %10s %10s %10s %10s %10s" %("DATA","MM","MM DW","MM UP","MC","MM max","Pred(MM+MC)","Discr.","[%]")
    print("Bin ## : %10s %10s %10s %10s %10s %10s %10s" %("DATA","MM","MM DW","MM UP","MC","MM max","Pred(MM+MC)"))

    hsumMC_unc = (stack[dtkey.replace('trueFAKE','trueREAL')].GetStack().Last()).Clone("h_lep_met_MMSS_zveto20_trueREAL_all_unccalc")

    for ibin in range(1,hsumMC_unc.GetNbinsX()+1):



        MM_max = hist[nomkey][datakey].GetBinContent(ibin)# + hist[upkey][datakey].GetBinContent(ibin)
        pred =  hsumMC_unc.GetBinContent(ibin)# + hist[upkey][datakey].GetBinContent(ibin)
        discrp = hist[dtkey][datakey].GetBinContent(ibin) - pred
        if MM_max != 0:
            ratio = (discrp/MM_max)*100.
        else:
            ratio = 0.0

        #print "Bin %02d : %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f " %(ibin,hist[dtkey][datakey].GetBinContent(ibin),hist[nomkey][datakey].GetBinContent(ibin),hist[dwkey][datakey].GetBinContent(ibin),hist[upkey][datakey].GetBinContent(ibin),hsumMC_unc.GetBinContent(ibin)-hist[nomkey][datakey].GetBinContent(ibin),MM_max,pred,discrp,ratio)

        print("Bin %02d : %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f " %(ibin,hist[dtkey][datakey].GetBinContent(ibin),hist[nomkey][datakey].GetBinContent(ibin),hist[dwkey][datakey].GetBinContent(ibin),hist[upkey][datakey].GetBinContent(ibin),hsumMC_unc.GetBinContent(ibin)-hist[nomkey][datakey].GetBinContent(ibin),MM_max,pred))
