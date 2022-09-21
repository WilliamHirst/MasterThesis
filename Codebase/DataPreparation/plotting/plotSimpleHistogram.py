from ROOT import TRandom, TFile, TH1, TH2, THStack, TH1F,TLegend, TCanvas, gROOT, gStyle, TH2F, TGraph, TPad, TColor, kWhite, kBlack, kDashed, kYellow, kMagenta, TGraphAsymmErrors, kGray
from optparse import OptionParser
from os import listdir
from os.path import isfile, join
from samples import configure_samples
import array
from math import *
import sys,os

d_samp,d_type,d_reg = configure_samples()
gStyle.SetPaintTextFormat(".2f");
TH1.SetDefaultSumw2()
TH2.SetDefaultSumw2()

parser = OptionParser()
parser.add_option("--doMC", help="do mc", default=0)
parser.add_option("--inputdir", help="Location of files", default="")
parser.add_option("--datakey", help="data15-16 or data17", default="")
parser.add_option("--histname", help="name of histogram", default="")
parser.add_option("--setBATCH", help="turn of batch mode", default=0)
parser.add_option("--estfake", help="estimate fakes from data", default=0)
(options, args) = parser.parse_args()

if int(options.setBATCH):
    gROOT.SetBatch(True)

doMC      = int(options.doMC)
inputdir  = str(options.inputdir)
datakey   = str(options.datakey)
histname  = str(options.histname).split(",") # h_lep_mt2_ALL_RJRstd_estfake
estfake   = int(options.estfake)

def getSortedList(nev):
    bkg_sorted_nev = []
    while len(bkg_sorted_nev) != len(nev.keys()):
        maxi = -999
        for key in nev.keys():
            if 'data' in key: continue
            if key in bkg_sorted_nev: continue
            if nev[key] > maxi:
                maxi = nev[key]
                maxi_key = key
        bkg_sorted_nev.append(maxi_key)
    return bkg_sorted_nev

evcount = {}
hists = {}
for h in histname:
    if not h in hists.keys():
        hists[h] = {}
        evcount[h] = {}
    for typ in d_samp.keys():
        if "data" in typ and not typ == datakey: continue
        fname = "%s/%s_merged_processed_HISTOGRAMS_0_0.root" %(inputdir,typ)
        if not os.path.isfile(fname):
            print("WARNING \t Could not find file for %s" %typ)
            continue
        tfile = TFile(fname)
        if estfake:
            if "data" in typ: hname = h
            else: hname = h + "_trueREAL"
        else: hname = h
        hists[h][typ] = tfile.Get(hname).Clone(h+"_"+typ)
        hists[h][typ] = hists[h][typ].Rebin(5)
        hists[h][typ].SetDirectory(0)
        hists[h][typ].SetLineColor(kBlack)
        #hists[h][typ].SetMarkerColor(d_samp[typ]['f_color'])
        if not "data" in typ:
            hists[h][typ].SetFillColor(d_samp[typ]['f_color'])
            evcount[h][typ] = hists[h][typ].Integral()
            #hists[h][typ].Scale(0.46)
        else:
            hists[h][typ].SetMarkerStyle(20)
            hists[h][typ].SetMarkerColor(kBlack)

        if typ == datakey and estfake:
            hname =  h + "_estfake"
            print hname
            hists[h]["fake"] = tfile.Get(hname).Clone(h+"_"+"fake")
            hists[h]["fake"].SetDirectory(0)
            hists[h]["fake"].SetFillColor(d_samp["fake"]['f_color'])
            hists[h]["fake"].SetLineColor(kBlack)
            evcount[h]["fake"] = hists[h]["fake"].Integral()
            hists[h]["fake"] = hists[h]["fake"].Rebin(5)
            

canv = []
canv.append(TCanvas("c_%s"%histname,"c_%s"%histname,1))
canv[-1].cd()
canv[-1].SetLogy()
hstack = {}
leg = TLegend(0.729226,0.56962,0.855301,0.888186)
leg.SetBorderSize(0)
leg.SetTextFont(42)
leg.SetTextSize(0.025)
leg.SetTextColor(1)
leg.SetFillColor(0)
leg.SetLineColor(0)
nh = 0

r = TRandom();
for histn in histname:

    if len(histname) > 1:
        color = int((113-51)*r.Rndm()+51)
    else:
        color = -1

    #histn = h

    hstack[histn] = THStack()
    
    bkg_sorted = getSortedList(evcount[histn])


    for typ in reversed(bkg_sorted):
        print("Adding",typ)
        hists[histn][typ].SetFillStyle(3305+nh*10)
        if color >= 0:
            hists[histn][typ].SetFillColor(color)
            hists[histn][typ].SetLineColor(color)
        hstack[histn].Add(hists[histn][typ])
        if len(histname) > 1:
            leg.AddEntry(hists[histn][typ],d_samp[typ]['leg']+" "+histn,"f")
        else:
            leg.AddEntry(hists[histn][typ],d_samp[typ]['leg'],"f")
        if nh == 0:
            hstack[histn].Draw("hist")
            hstack[histn].GetXaxis().SetRangeUser(0,250)
            
        else:
            hstack[histn].Draw("same hist")
            
    if datakey in hists[histn].keys():
        hists[histn][datakey].Draw("e p same")
        leg.AddEntry(hists[histn][datakey],d_samp[datakey]['leg'],"lp")
    nh += 1
leg.Draw()


# hsumMC = (hstack.GetStack().Last()).Clone(histname+"_all")
# tot_scalef = 0.0
# nbins = 0.0
# for i in range(1,hsumMC.GetNbinsX()+1):
#     mc_binc   = hsumMC.GetBinContent(i)
#     data_binc = hists[histname][datakey].GetBinContent(i)
#     if mc_binc:
#         print "scalef %i = %.4f" %(i,data_binc/mc_binc)
#         tot_scalef += data_binc/mc_binc
#         nbins += 1


# print "Total scalef is %.4f" %(tot_scalef/nbins)

