import sys
import ROOT as R
from ROOT import kGray
from os import listdir, system
from os.path import isfile, isdir, join
import array
from math import sqrt
astyle = "/home/eirikgr/atlasrootstyle/"#
#astyle ="/mn/felt/u1/eirikgr/bin/atlasstyle-00-04-02/"
R.gROOT.SetMacroPath(astyle)
R.gROOT.LoadMacro("AtlasUtils.C")
R.gROOT.LoadMacro("AtlasStyle.C")
R.gROOT.LoadMacro("AtlasLabels.C")
R.SetAtlasStyle() 

R.gStyle.SetErrorX(0.5) 

R.gROOT.SetBatch(True)

def doRebin(hT,bins=array.array('f',[50,60,70,80,90,100,120,150,200,300,1000])):
    #print("Histogram is %s"%hT.GetName()) 
    VB = 1
    name = ""
    name = hT.GetName()+"_REBIN"
    hT_reb = R.TH1F(name,name,len(bins)-1,bins)
    tot_err = {}
    for nbin in range(1,hT.GetNbinsX()+2):
        newbin = hT_reb.FindBin(abs(hT.GetBinLowEdge(nbin)))
        if VB > 1: print("Adding %.2f from value %.2f into bin %i" %(hT.GetBinContent(nbin),abs(hT.GetBinLowEdge(nbin)),newbin))
        hT_reb.AddBinContent(newbin,hT.GetBinContent(nbin))
        if not newbin in tot_err.keys():
            tot_err[newbin] = 0.0
        tot_err[newbin] += (hT.GetBinError(nbin)*hT.GetBinError(nbin))
    for newbin in tot_err.keys():
        if VB > 1: print("Setting error for bin %i = %1f"%(newbin,sqrt(tot_err[newbin])))
        hT_reb.SetBinError(newbin,sqrt(tot_err[newbin]) if tot_err[newbin] else 0.0)
    return hT_reb


def getTreeName(f1,docheck = False):
    print("Checking %s"%f1.GetName())
    dirlist = f1.GetListOfKeys()
    it = dirlist.MakeIterator()
    key = it.Next()
    while key:
        cl = R.gROOT.GetClass(key.GetClassName());
        if cl.InheritsFrom("TTree"): 
            obj = key.ReadObj()
            nev = obj.Draw("preselection>>hist","preselection == 1","goff")
            #hist = R.gDirectory.Get("hist")
            #nev = hist.GetEntries()
            if not nev:
                print("Tree is empty in %s, skipping"%f1.GetName())
                return "noname"
            return obj.GetName()
        else:
            try:
                key = iter.Next()
            except:
                print("Problems with file %s"%f1.GetName())
                return "noname"
            continue
    return "noname"


R.EnableImplicitMT(60)#EnableImplicitMT()

inputdir = sys.argv[1]

nom_samples = [363356,
               363358,
               363359,
               363360,
               363489,
               364250,
               364253,
               364254,
               364255,
               700320,
               700321,
               700322,
               700323,
               700324,
               700325,
               700326,
               700327,
               700328,
               700452,
               700453,
               700454,
               700455,
               700456,
               700457,
               700458,
               700459,
               700460,
               700467,
               700468,
               700469,
               700470,
               700471,
               700472,
               700473,
               700474,
               700475]

regions = ["SR1","preselection","SR2","SR3","CR_Dib","CR_Top","CR_Z","VR1","VR2","VR3"]


regions = ["VR1","VR2","VR3"]

cutstr = {"SR1":"(metrel > 50 && metrel < 100) && dphi_ll_met > 2.5 && nbjets == 0 && (mll > 180 && mll < 1500)",
          "SR2":"(metrel > 100 && metrel < 150) && dphi_ll_met > 2.5 && nbjets == 0 && (mll > 180 && mll < 1500)",
          "SR3":"(metrel > 150) && dphi_ll_met > 2.5 && nbjets == 0 && (mll > 180 && mll < 1500)",
          "VR1":"(metrel > 50 && metrel < 100)  && dphi_ll_met > 2.5 && nbjets == 0 && mll > 110",# && (mll > 110 && mll < 180)",
          "VR2":"(metrel > 100 && metrel < 150) && dphi_ll_met > 2.5 && nbjets == 0 && mll > 110",# && (mll > 110 && mll < 180)",
          "VR3":"(metrel > 150) && dphi_ll_met > 2.5 && nbjets == 0 && mll > 110",# && (mll > 110 && mll < 180)",
          "CR_Dib":"(mll > 70 && mll < 110) &&  nbjets == 0 && @jet_e.size() == 0 &&  dphi_ll_met > 2.5 && metsign > 10",
          "CR_Top":"(mll > 180) &&  nbjets >= 1 &&  dphi_ll_met > 2.5 && metrel > 50",
          "CR_Z":"(mll > 70 && mll < 110) &&  nbjets == 0 && dphi_ll_met > 2.5 && metrel > 50"}

histo = {}
df = {}
zjt_samples = {}
dib_samples = {}

files = [join(inputdir, f) for f in listdir(inputdir) if (isfile(join(inputdir, f)) and f.endswith(".root"))]


infiles = {}

evcount = {}

for inputfile in files:


    rootfile = inputfile.split("/")[-1]

    print(rootfile)

    
    #print(inputfile.replace(".root",".txt"))
    
        
    #def plotWithRatio(canvas,this_histo,nom_histo):

        #h3.SetDirectory(0)
        #can.Update()
        #return ca
        
    #f1 = R.TFile.Open(inputfile)
    #tname = getTreeName(f1,False)
    #f1.Close()
    tname = "ntuple"
    if not tname == "ntuple":
        print("Skipping %s since tree not found"%(rootfile))
        print("mv %s /scratch/eirikgr/theorySystNtuples/untared/empty/%s."%(inputfile,rootfile))
        system("mv %s /scratch/eirikgr/theorySystNtuples/untared/empty/%s."%(inputfile,rootfile))
        continue
        

    dsid = int(rootfile.split("_")[0])
    descr = rootfile.split(".")[0].replace(str(dsid),"").replace("ZMETdilep","")[1:-1]
    #print(dsid,descr)
    systvar = "nom"
    if not dsid in nom_samples:
        systvar = descr.split("_")[-1]
    non_syst_descr = ""
    if "Sh_2211" in descr:
        strt_idx = 2
        if not systvar in zjt_samples.keys():
            zjt_samples[systvar] = []
        zjt_samples[systvar].append(inputfile)
    elif "Sherpa_22" in descr:
        strt_idx = 3
        if not systvar in dib_samples.keys():
            dib_samples[systvar] = []
        dib_samples[systvar].append(inputfile)
    for d in descr.split("_")[strt_idx:]:
        if d == systvar: break
        if d in ['CFilterBVeto','CVetoBVeto','BFilter']:
            break
        non_syst_descr += "%s_"%d
    non_syst_descr = non_syst_descr[:-1]

    if not non_syst_descr in df.keys():
        df[non_syst_descr] = {}
        histo[non_syst_descr] = {}
        infiles[non_syst_descr] = {}
        evcount[non_syst_descr] = {}
    if not systvar in infiles[non_syst_descr]:
        infiles[non_syst_descr][systvar] = []
        evcount[non_syst_descr][systvar] = 0
    infiles[non_syst_descr][systvar].append(inputfile)

    txtname = inputfile.replace(".root",".txt")
    if isfile(txtname):
        print("Found %s"%txtname)
        lines = [line.rstrip() for line in open(txtname)]
        il = 0
        for l in lines:
            if il == 0:
                il += 1
                continue
            if l.split(",")[0] == "All":
                nev = float(l.split(",")[-1])
                evcount[non_syst_descr][systvar] = nev
                break
            il += 1
            
      
for nsd in infiles.keys():
    for systvar in infiles[nsd].keys():
        #else:
        print("INFO Adding %i files for samples %s and systematics %s"%(len(infiles[nsd][systvar]),nsd,systvar))
        if len(infiles[nsd][systvar]) > 1:
            for n in range(len(infiles[nsd][systvar])):
                print("%i : %s"%(n,infiles[nsd][systvar][n]))
        if not systvar in df[nsd]:
            df[nsd][systvar] = R.RDataFrame("ntuple",infiles[nsd][systvar])
            #df[nsd][systvar].Define("nev",nev)
            histo[nsd][systvar] = {}


allhistograms = []
for proc in df.keys():
    for syst in df[proc].keys():
        for r in regions:
            if r in cutstr.keys():
                cut = cutstr[r]
            else:
                cut = r
            try:
                #print(r,cut)
                histo[proc][syst]["%s"%r] = df[proc][syst].Filter(cut).Histo1D(("metrel_%s_%s_%s"%(proc,syst,r),"metrel_%s_%s_%s"%(proc,syst,r),500,0,1000),"metrel","eventWeight")
                #histo[proc][syst]["%s"%r] = doRebin(histo[proc][syst]["%s"%r])#.Rebin(10,"h_rebin_metrel_%s_%s_%s"%(proc,syst,r))
                allhistograms.append(histo[proc][syst]["%s"%r])
                #histo[proc][syst]["%s"%r].SetDirectory()
            except:
                print("WARNING \t Could not create histogram for %s %s %s"%(proc,syst,r))
            


dib_merged_df = {}
dib_histo = {}
zjt_merged_df = {}
zjt_histo = {}

for syst in dib_samples.keys():
    dib_histo[syst] = {}
    dib_merged_df[syst] = R.RDataFrame("ntuple",dib_samples[syst])
    #print("dib")
    #print("\n".join(dib_samples[syst]))
    for r in regions:
        if r in cutstr.keys():
            cut = cutstr[r]
        else:
            cut = r
        try:
            dib_histo[syst]["%s"%r] = dib_merged_df[syst].Filter(cut).Histo1D(("metrel_dib_%s_%s"%(syst,r),"metrel_dib_%s_%s"%(syst,r),500,0,1000),"metrel","eventWeight")
            allhistograms.append(dib_histo[syst]["%s"%r])
        except:
            print("WARNING \t Problems with %s, %s, %s"%(syst,r,"\n".join(dib_samples[syst])))
for syst in zjt_samples.keys():
    zjt_histo[syst] = {}
    zjt_merged_df[syst] = R.RDataFrame("ntuple",zjt_samples[syst])
    #print("zjets")
    #print("\n".join(zjt_samples[syst]))
    for r in regions:
        if r in cutstr.keys():
            cut = cutstr[r]
        else:
            cut = r
        try:
            zjt_histo[syst]["%s"%r] = zjt_merged_df[syst].Filter(cut).Histo1D(("metrel_zjt_%s_%s"%(syst,r),"metrel_zjt_%s_%s"%(syst,r),500,0,1000),"metrel","eventWeight")
            allhistograms.append(zjt_histo[syst]["%s"%r])
        except:
            print("WARNING \t Problems with %s, %s, %s"%(syst,r,"\n".join(zjt_samples[syst])))
        
#sys.exit()

print("INFO \t Executing %i histograms"%len(allhistograms))
R.RDF.RunGraphs(allhistograms);
descr = "llvv"




#this_histo = histo[descr]["SR1"]
#nom_histo = histo[descr]["SR1"].Clone("to_divide")

plot_dic = {'CKKW15':R.kRed+1, 'nom':R.kBlack, 'QSF4':R.kAzure+7, 'QSF025':R.kSpring-1, 'CKKW30':R.kYellow+1,"QSFDN":R.kMagenta+1,"QSFUP":R.kCyan+1,"CSSKIN":R.kOrange+7}

doregion = ""

plothisto = {}
#plothisto["all"] = dib_histo
plothisto = histo

nsyst = {}
for desc in plothisto.keys():
    print("Background %s: %s"%(desc,",".join(plothisto[desc].keys())))
    nsyst[desc] = len(plothisto[desc].keys())

h_syststack = {}
h_errstack = {}
yaxis = {}
pad1 = {}
pad2 = {}
can = {}
for descr in plothisto.keys():
    if not descr in h_syststack.keys():
        h_syststack[descr] = {}
        h_errstack[descr] = {}
        yaxis[descr] = {}
        can[descr] = {}
        pad1[descr] = {}
        pad2[descr] = {}
    for doregion in regions:
        if not doregion in h_syststack[descr].keys():
            h_syststack[descr][doregion] = R.THStack("hstack_%s_%s"%(descr,doregion),"hstack_%s_%s"%(descr,doregion))
            h_errstack[descr][doregion] = R.THStack("hstack_err_%s_%s"%(descr,doregion),"hstack_err_%s_%s"%(descr,doregion))
            can[descr][doregion]  = R.TCanvas('can_%s_%s'%(descr,doregion),'',900,700)
            pad1[descr][doregion] = R.TPad('pad1_%s_%s'%(descr,doregion), '', 0.0, 0.40, 1.0, 1.0)
            pad2[descr][doregion] = R.TPad('pad2_%s_%s'%(descr,doregion), '', 0.0, 0.00, 1.0, 0.4)
        pad1[descr][doregion].SetBottomMargin(0.005)
        pad2[descr][doregion].SetBottomMargin(0.3)
        pad1[descr][doregion].SetLeftMargin(0.1)
        pad2[descr][doregion].SetLeftMargin(0.1)
        pad1[descr][doregion].Draw()

        pad2[descr][doregion].Draw()



        legend1 = R.TLegend(0.67, 0.60, 0.99, 0.6+(nsyst[descr]*0.05))
        legend1.SetTextFont(42)
        legend1.SetFillStyle(0)
        legend1.SetBorderSize(0)
        legend1.SetTextSize(0.05)

        nh = 0
        h_nom = 0
        pad1[descr][doregion].cd()
        #pad1[descr][doregion].SetLogy()
        #pad1[descr][doregion].SetLogx()
        absmax = -999
        absmin = 999
        
        
        for syst in plothisto[descr].keys():
            if doregion and not doregion in plothisto[descr][syst].keys():
                continue
            if syst in plot_dic.keys():
                color = plot_dic[syst]
            #print("Bin width = %i"%plothisto[descr][syst][doregion].GetBinWidth(1))
            #for ibin in range(1,plothisto[descr][syst][doregion].GetNbinsX()+1):
            #    print("BF RB Bin %i ) %.2f +/- %.2f"%(ibin,plothisto[descr][syst][doregion].GetBinContent(ibin),plothisto[descr][syst][doregion].GetBinError(ibin)))
            plothisto[descr][syst][doregion].Sumw2()
            minx = 0
            maxx = 0
            maxy = 0
            miny = 0
            if doregion in ["VR1","SR1"]:
                bins=array.array('f',[50,100])#[50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100])#[50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100]
                minx = 50
                maxx = 100
                maxy = 3e-4
                miny = 0
            elif doregion in ["VR2","SR2"]:
                bins=array.array('f',[100,150])#[100,110,120,130,140,150]) #[100,110,120,130,140,150]
                minx = 100
                maxx = 150
            elif doregion in ["VR3","SR3"]:
                bins=array.array('f',[150,5000])#[150,200,500])#[150,200,500]
                minx = 150
                maxx = 500
                #miny = 0.75
                #maxy = 1.2
                pad1[descr][doregion].SetLogy(False)
            plothisto[descr][syst][doregion] = doRebin(plothisto[descr][syst][doregion],bins)#.Rebin(10,"h_rebin_metrel_%s_%s_%s"%(proc,syst,r))
            #for ibin in range(1,plothisto[descr][syst][doregion].GetNbinsX()+1):
            #    print("%s,%s,%s : AF RB Bin %i ) %.2f +/- %.2f"%(descr,syst,doregion,ibin,plothisto[descr][syst][doregion].GetBinContent(ibin),plothisto[descr][syst][doregion].GetBinError(ibin)))
            #plothisto[descr][syst][doregion] = plothisto[descr][syst][doregion].Rebin(10,"h_rebin_metrel_%s_%s_%s"%(descr,syst,r))
            plothisto[descr][syst][doregion].SetDirectory(0)
            plothisto[descr][syst][doregion].SetLineColor(color)
            plothisto[descr][syst][doregion].SetMarkerColor(color)
            plothisto[descr][syst][doregion].SetMarkerStyle(1)
            isscaled = False
            #if plothisto[descr][syst][doregion].Integral(0,plothisto[descr][syst][doregion].GetNbinsX()+1) != 0:
            #    plothisto[descr][syst][doregion].Scale(1.0/plothisto[descr][syst][doregion].Integral(0,plothisto[descr][syst][doregion].GetNbinsX()+1))
            plothisto[descr][syst][doregion].Scale(1.0/evcount[descr][syst])
            isscaled = False
            #isscaled = False
            #else:
            #    continue
            #if nh == 0:
            #    if minx and maxx:
            #        plothisto[descr][syst][doregion].GetXaxis().SetRangeUser(minx,maxx)
            #    if miny and maxy:
            #        plothisto[descr][syst][doregion].GetYaxis().SetRangeUser(miny,maxy)
            #    elif isscaled:
            #        plothisto[descr][syst][doregion].GetYaxis().SetRangeUser(0.001,1.0)
            #    #hh = plothisto[descr][syst][doregion].DrawCopy("hist")
            h_syststack[descr][doregion].Add(plothisto[descr][syst][doregion].Clone(plothisto[descr][syst][doregion].GetName()+"_forstack"))
            #    plothisto[descr][syst][doregion].SetMinimum(0.12)
            #    pl_descr = descr
            #    pl_syst = syst
            #    pl_doregion = doregion
            #else:
            #    plothisto[descr][syst][doregion].DrawCopy("hist same")
            plothisto[descr][syst][doregion].SetFillStyle(3018)
            plothisto[descr][syst][doregion].SetFillColor(color)
            #plothisto[descr][syst][doregion].Draw("e2same")
            h_errstack[descr][doregion].Add(plothisto[descr][syst][doregion].Clone(plothisto[descr][syst][doregion].GetName()+"_unc_forstack"))

            newmax = plothisto[descr][syst][doregion].GetMaximum()
            newmin = plothisto[descr][syst][doregion].GetMinimum()
            if newmax > absmax:
                absmax = newmax
            if newmin < absmin:
                absmin = newmin
            #hh.SetFillStyle(3018)
            #hh.SetFillColor(color)
            #hh.Draw("e2same")
            #pad1[descr][doregion].Update()
            #for ibin in range(1,plothisto[descr][syst][doregion].GetNbinsX()+1):
            #    print("AF SC Bin %i ) %.4f +/- %.10f"%(ibin,plothisto[descr][syst][doregion].GetBinContent(ibin),plothisto[descr][syst][doregion].GetBinError(ibin)))
            legend1.AddEntry(plothisto[descr][syst][doregion],syst,"l")
            if "nom" in syst:
                h_nom = plothisto[descr][syst][doregion].Clone("h_nom")
                h_nom.SetDirectory(0)
            nh += 1
            
            #can[descr][doregion].cd()
        #hh.SetMinimum(0.01)
        #hh.SetMaximum(absmax*1.05)
        #print(pl_descr,pl_syst,pl_doregion,absmax)
        #plothisto[pl_descr][pl_syst][pl_doregion].GetYaxis().SetRangeUser(0,absmax)
        #plothisto[pl_descr][pl_syst][pl_doregion].SetMaximum(absmax)
        #pad1[descr][doregion].Update()
        # pad1[descr][doregion].cd()
        # h_syststack[descr][doregion].Draw("hist unstacked")
        # #h_syststack[descr][doregion].Draw("e2same unstacked")
        # yaxis[descr][doregion] = h_syststack[descr][doregion].GetYaxis()
        
        # maxi = (h_syststack[descr][doregion].GetStack().Last()).GetMaximum()

        # yaxis[descr][doregion].SetRangeUser(0,maxi)
        # print(descr,doregion,maxi)
        # #h_syststack[descr][doregion].GetHistogram().SetMaximum(maxi)
        # #h_syststack[descr][doregion].GetYaxis().SetRangeUser(0,maxi)
        # h_syststack[descr][doregion].Draw("hist unstacked")
        # #h_syststack[descr][doregion].Draw("hist unstacked")
        # #h_syststack[descr][doregion].Draw("e2same unstacked")
        # #pad1[descr][doregion].Update()
        # pad1[descr][doregion].Modified()
        # pad1[descr][doregion].Update()
        # legend1.SetHeader("Region : %s"%doregion)
        # legend1.Draw()            

        newplothisto = {}

        legend = R.TLegend(0.10, 0.10, 0.15, 0.15)
        legend.SetTextFont(42)
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        legend.SetTextSize(0.08)
        pad2[descr][doregion].cd()
        pad2[descr][doregion].SetGridy()
        #pad2[descr][doregion].SetLogx()
        nh = 0
        absmax = -999
        absmin = 999
        plotsyst = ""
        total_syst_up = []
        total_syst_dw = []
        summary_str = ""
        for syst in plothisto[descr].keys():
            #print(syst)
            if "nom" in syst: continue
            if not doregion in plothisto[descr][syst].keys(): continue
            try:
                newplothisto[syst].Reset()
                newplothisto["total_up"].Reset()
                newplothisto["total_dw"].Reset()
            except:
                print("INFO \t Hisogram does not exist. Good.")
            newplothisto[syst] = plothisto[descr][syst][doregion].Clone("h_%s_%s_%s"%(descr,syst,doregion))
            #for ibin in range(1,newplothisto[syst].GetNbinsX()+1):
            #    print("%s,%s,%s : BF SUB SYST Bin %i ) %.2f +/- %.2f"%(descr,syst,doregion,ibin,newplothisto[syst].GetBinContent(ibin),newplothisto[syst].GetBinError(ibin)))
            #for ibin in range(1,h_nom.GetNbinsX()+1):
            #    print("%s,%s,%s : AF SUB NOM  Bin %i ) %.2f +/- %.2f"%(descr,syst,doregion,ibin,h_nom.GetBinContent(ibin),h_nom.GetBinError(ibin)))
            
            newplothisto[syst].Add(h_nom,-1.0)
            #for ibin in range(1,newplothisto[syst].GetNbinsX()+1):
            #    print("%s,%s,%s : AF SUB Bin %i ) %.2f +/- %.2f"%(descr,syst,doregion,ibin,newplothisto[syst].GetBinContent(ibin),newplothisto[syst].GetBinError(ibin)))
            newplothisto[syst].Divide(newplothisto[syst],h_nom,1,1,"B")

            
            # asym = R.TGraphAsymmErrors();
            # asym.Divide(newplothisto[syst],h_nom,"cl=0.683 b(1,1) mode")
            # ni = 0
            # for i in range(0,newplothisto[syst].GetNbinsX()+1):
            #     if newplothisto[syst].GetBinContent(i+1):
            #         newplothisto[syst].SetBinError(i+1,asym.GetErrorY(ni))
            #         ni += 1
            #     else:
            #         newplothisto[syst].SetBinError(i+1,0.0)
                    

            
            if syst in plot_dic.keys():
                color = plot_dic[syst]
            newplothisto[syst].SetFillStyle(0)
            newplothisto[syst].SetLineColor(color)
            newplothisto[syst].SetLineStyle(9)
            newmax = newplothisto[syst].GetMaximum()
            newmin = newplothisto[syst].GetMinimum()
            if newmax > absmax:
                absmax = newmax
            if newmin < absmin:
                absmin = newmin
            if nh == 0:
                newplothisto[syst].GetYaxis().SetLabelSize(0.08)
                newplothisto[syst].GetYaxis().SetTitle("#frac{N_{syst}-N_{nom}}{N_{nom}}")
                newplothisto[syst].GetYaxis().SetTitleOffset(0.55)
                newplothisto[syst].GetYaxis().SetTitleSize(0.07)
                newplothisto[syst].GetXaxis().SetTitle("E_{T}^{miss,rel}")
                newplothisto[syst].GetXaxis().SetTitleSize(0.07)
                newplothisto[syst].GetXaxis().SetLabelSize(0.08)
                hplotsyst = newplothisto[syst].DrawCopy("hist")
                plotsyst = syst
            else:
                newplothisto[syst].DrawCopy("hist same")

            newplothisto[syst].SetFillStyle(3018)
            newplothisto[syst].SetFillColor(color)
            newplothisto[syst].Draw("e2same")
            #legend.AddEntry(newplothisto[syst],"%s"%syst,"l")
            summary_str += "-"*60
            summary_str += "\nSystematics %s\n"%syst
            summary_str += "-"*60
            summary_str += "\n"
            for ibin in range(1,newplothisto[syst].GetNbinsX()+1):
                summary_str += "[%.0f-%.0f] %.5f +/- %.2f\n"%(newplothisto[syst].GetBinLowEdge(ibin),newplothisto[syst].GetBinLowEdge(ibin+1),newplothisto[syst].GetBinContent(ibin)*100.0,newplothisto[syst].GetBinError(ibin)*100.0)
        

            if len(total_syst_up) == 0:
                total_syst_up = [0] * (newplothisto[syst].GetNbinsX()+1)
            if len(total_syst_dw) == 0:
                total_syst_dw = [0] * (newplothisto[syst].GetNbinsX()+1)
            for ibin in range(1,newplothisto[syst].GetNbinsX()+1):
                if newplothisto[syst].GetBinContent(ibin) >= 0:
                    total_syst_up[ibin-1] += (newplothisto[syst].GetBinContent(ibin)*newplothisto[syst].GetBinContent(ibin))
                    total_syst_dw[ibin-1] += 0
                else:
                    total_syst_dw[ibin-1] += (newplothisto[syst].GetBinContent(ibin)*newplothisto[syst].GetBinContent(ibin))
                    total_syst_up[ibin-1] += 0
            
            nh += 1
        
        if len(newplothisto.keys()) == 0:
            continue
        key0 = list(newplothisto.keys())[0]
        newplothisto["total_up"] = newplothisto[key0].Clone(newplothisto[key0].GetName()+"_total_up")
        newplothisto["total_dw"] = newplothisto[key0].Clone(newplothisto[key0].GetName()+"_total_dw")
        ibin = 1
        for ts in total_syst_up:
            newplothisto["total_up"].SetBinContent(ibin,sqrt(ts))
            ibin += 1
        ibin = 1
        for ts in total_syst_dw:
            newplothisto["total_dw"].SetBinContent(ibin,-1*sqrt(ts) if ts > 0 else 0.0)
            ibin += 1
        newplothisto["total_up"].SetLineWidth(3)
        newplothisto["total_up"].SetLineColor(R.kBlack)
        #newplothisto["total_up"].Draw("hist same")
        newplothisto["total_dw"].SetLineWidth(3)
        newplothisto["total_dw"].SetLineColor(R.kBlack)
        #newplothisto["total_dw"].Draw("hist same")
        #legend.AddEntry(newplothisto["total_up"],"Total","l")
        #newmax = newplothisto["total_up"].GetMaximum()
        #newmin = newplothisto["total_dw"].GetMinimum()
        #if newmax > absmax:
        #    absmax = newmax
        #if newmin < absmin:
        #    absmin = newmin
        
        if plotsyst:
            if not "VR3" in doregion:
                absmax = 0.4
            else:
                absmax = 0.5
            absmin = -0.4
            #absmax = 0.4
            #absmin = -0.4
            #if absmax < 0.02:
            #    absmax = 0.05
            #if absmin > -0.01:
            #    absmin = -0.05
            hplotsyst.GetYaxis().SetRangeUser(absmin,absmax)
            #newplothisto[plotsyst].SetMinimum(absmin*1.1)
            #print("descr = %s syst = %s doregion = %s, max = %.1f min = %.1f"%(descr,syst,doregion,absmax,absmin))
        line = R.TLine(minx,0,maxx,0)
        line.SetLineWidth(2)
        line.Draw()
        pad2[descr][doregion].Update()
        legend.Draw()


        pad1[descr][doregion].cd()
        h_syststack[descr][doregion].Draw("hist NOSTACK")
        
        yaxis[descr][doregion] = h_syststack[descr][doregion].GetYaxis()
        
        maxi = (h_syststack[descr][doregion].GetStack().Last()).GetMaximum()

        yaxis[descr][doregion].SetRangeUser(0,maxi)
        print(descr,doregion,maxi)
        #h_syststack[descr][doregion].GetHistogram().SetMaximum(maxi)
        #h_syststack[descr][doregion].GetYaxis().SetRangeUser(0,maxi)
        h_syststack[descr][doregion].Draw("hist NOSTACK")
        h_errstack[descr][doregion].Draw("e2same NOSTACK")
        #h_syststack[descr][doregion].Draw("hist unstacked")
        #h_syststack[descr][doregion].Draw("e2same unstacked")
        #pad1[descr][doregion].Update()
        pad1[descr][doregion].Modified()
        pad1[descr][doregion].Update()
        can[descr][doregion].Modified()
        can[descr][doregion].Update()
        
        legend1.SetHeader("Region : %s"%doregion)
        legend1.Draw() 

        
        can[descr][doregion].SaveAs("syst_%s_%s.png"%(descr,doregion))
        can[descr][doregion].SaveAs("syst_%s_%s.pdf"%(descr,doregion))
        #break
        #break
            #h3 = this_plothisto.Clone("h3")
            #h3.Sumw2()
            #h3.Divide(nom_plothisto) # PHYS/SUSY2
            #h3.SetLineColor(kGray+2)
            #h3.SetLineWidth(2)
            ##
            #h3.Draw("e0p")
            #plotWithRatio(can[descr][doregion] ,plothisto[doregion],plothisto[doregion].Clone("to_divide"))
            #can[descr][doregion].Update()
            # for r in regions:
            #     hist[r] = R.TH1F("h_metrel_%s"%r,"h_metrel_%s"%r,500,0,1000)
            #     c.Draw("metrel>>h_metrel_%s"%r,"%s"%r)
        utfil = open("syst_summary_%s_%s.txt"%(descr,doregion),"a")
        utfil.write(summary_str)
        utfil.close()
