import ROOT as R
#from ROOT import *
import os, sys

#from infofile import infos

R.gROOT.SetBatch(1)
R.gStyle.SetOptStat(0);
R.gStyle.SetPadLeftMargin(0.13)
R.gStyle.SetLegendBorderSize(0)
R.gStyle.SetPalette(1)
R.gStyle.SetGridStyle(2)
R.gStyle.SetPadLeftMargin(0.13)
R.TH1.AddDirectory(False)

channel = sys.argv[1] 


# Variables 

variables = ['pt1', 'pt2', 'eta1', 'eta2', 'phi1', 'phi2', 'mll', 'met']

xtitles = {'pt1':'Leading lepton p_{T} (GeV)', 'pt2':'Subleading lepton p_{T} (GeV)', 'eta1':'Leading lepton #eta', 'eta2':'Subleading lepton #eta', 'phi1':'Leading lepton #phi', 'phi2':'Subleading lepton #phi', 'mll':'m_{ll} (GeV)', 'met':'E_{T}^{miss} (GeV)'}


# Backgrounds

backgrounds = ['Zjets', 'Top', 'Diboson', 'Wjets'] 

Zjets = [364100, 364101, 364102, 364103, 364104, 364105, 364106, 364107, 364108, 364109, 364110, 364111, 364112, 364113, 364114, 364115, 364116, 364117, 364118, 364119, 364120, 364121, 364122, 364123, 364124, 
         364125, 364126, 364127, 364128, 364129, 364130, 364131, 364132, 364133, 364134, 364135, 364136, 364137, 364138, 364139, 364140, 364141]

Wjets = [364156, 364157, 364158, 364159, 364160, 364161, 364162, 364163, 364164, 364165, 364166, 364167, 364168, 364169, 364170, 364171, 364172, 364173, 364174, 364175, 364176, 364177, 364178, 364179, 364180, 
         364181, 364182, 364183, 364184, 364185, 364186, 364187, 364188, 364189, 364190, 364191, 364192, 364193, 364194, 364195, 364196, 364197]

Diboson = [363356, 363358, 363359, 363360, 363489, 363490, 363491, 363492, 363493]

Top = [410000, 410011, 410012, 4100013, 410014, 410025, 410026]

# Signals

signals = ['Zprime2000']

Zprime2000 = [301215, 301220];


fileIDs = {'Diboson':Diboson, 'Zjets':Zjets, 'Wjets':Wjets, 'Top':Top, 'Zprime2000':Zprime2000}

hist_bkg = {}
for var in variables:
    hist_bkg[var] = {}
    for bkg in backgrounds:
        hist_bkg[var][bkg] = R.TH1F()

hist_sig = {}
for var in variables:
    hist_sig[var] = {}
    for sig in signals:
        hist_sig[var][sig] = R.TH1F()


colours = dict(Diboson=R.kAzure+1, Top=R.kRed+1, Zjets=R.kOrange-2, Wjets=R.kGray, Zprime2000=R.kBlue) 


# Extract info about cross section and sum of weights from infofile 

# info = {} 
# for key in list(infos.keys()): 
#     ID = infos[key]['DSID']
#     info[ID] = {} 
#     info[ID]['xsec'] = infos[key]['xsec'] 
#     info[ID]['sumw'] = infos[key]['sumw'] 
#     info[ID]['events'] = infos[key]['events']


# Function for making histograms
L = 10.06 # integrated luminosity = (10.06 +/- 0.37) fb^-1 

def fill_hist(h, h_name, key, ID):

    h_midl = infile.Get(h_name).Clone("h_midl")

    #xsec = 1000*info[ID]['xsec']
    #nev = info[ID]['sumw'] 

    #N_mc = xsec*L

    #sf = 1.0;N_mc/nev  

    if not h.GetName(): 
        h=infile.Get(h_name)  
        #h.Scale(sf)
        n = h.GetNbinsX()
        for i in range(n):
            bc = h.GetBinContent(i)
            if bc<0:
                h.SetBinContent(i,0)
        h.SetFillColor(colours[key])
        h.SetLineColor(colours[key])	 
    else:
        #h_midl.Scale(sf)
        n = h_midl.GetNbinsX()
        for i in range(n):
            bc = h_midl.GetBinContent(i)
            if bc < 0: 
                h_midl.SetBinContent(i,0) 
        h.Add(h_midl)

    return h  

## Loop over files in MC directory
print('Looping over MC files')
mc_files=os.listdir('./Histograms/MC')
print(mc_files)
for filename in os.listdir('Histograms/MC/'):
    if '.root' in filename: 
        filepath = 'Histograms/MC/'+filename 
        infile = R.TFile(filepath)
        file_id = int(filename.split('.')[2])
        print("Got MC file: %s" % filepath)
        for var in variables:
            for bkg in backgrounds:
                if file_id in fileIDs[bkg]: 
                    print("Adding %i as %s"%(file_id,bkg))
                    hist_bkg[var][bkg] = fill_hist(hist_bkg[var][bkg], 'h_'+channel+'_'+var, bkg, file_id)
            for sig in signals:
                if file_id in fileIDs[sig]: 
                    hist_sig[var][sig] = fill_hist(hist_sig[var][sig], 'h_'+channel+'_'+var, sig, file_id)


# Get data 

data = R.TFile('Histograms/Data/hist.Data.2016.root')
hist_d ={}

for var in variables:
    hist_d[var] = data.Get('h_'+channel+'_'+var) 
    hist_d[var].SetMarkerStyle(20)
    hist_d[var].SetMarkerSize(0.7)
    hist_d[var].SetLineColor(R.kBlack)
    hist_d[var].GetYaxis().SetTitle("Events")
    hist_d[var].GetXaxis().SetTitle(xtitles[var]) 
    hist_d[var].GetXaxis().SetTitleFont(43)
    hist_d[var].GetXaxis().SetTitleSize(16)
    hist_d[var].GetYaxis().SetTitleFont(43)
    hist_d[var].GetYaxis().SetTitleSize(16)
    hist_d[var].GetXaxis().SetLabelFont(43)
    hist_d[var].GetXaxis().SetLabelSize(16)
    hist_d[var].GetYaxis().SetLabelFont(43)
    hist_d[var].GetYaxis().SetLabelSize(16)
    hist_d[var].GetXaxis().SetTitleOffset(4)
    hist_d[var].GetYaxis().SetTitleOffset(1.5)


# Style histograms, make stack and histograms with full background

stack = {} 
hist_r = {}
hist_mc = {}

for var in variables:
    stack[var] = R.THStack(var, "")  
    hist_mc[var] = R.TH1F()
    hist_r[var] = R.TH1F()
    for bkg in reversed(backgrounds): 
        hist_bkg[var][bkg].GetYaxis().SetTitle("Events")
        hist_bkg[var][bkg].GetXaxis().SetTitle(xtitles[var]) 
        hist_bkg[var][bkg].GetXaxis().SetTitleFont(43)
        hist_bkg[var][bkg].GetXaxis().SetTitleSize(16)
        hist_bkg[var][bkg].GetYaxis().SetTitleFont(43)
        hist_bkg[var][bkg].GetYaxis().SetTitleSize(16)
        hist_bkg[var][bkg].GetXaxis().SetLabelFont(43)
        hist_bkg[var][bkg].GetXaxis().SetLabelSize(16)
        hist_bkg[var][bkg].GetYaxis().SetLabelFont(43)
        hist_bkg[var][bkg].GetYaxis().SetLabelSize(16)
        hist_bkg[var][bkg].GetXaxis().SetTitleOffset(4)
        hist_bkg[var][bkg].GetYaxis().SetTitleOffset(1.5)
        stack[var].Add(hist_bkg[var][bkg]) 
        if not hist_mc[var].GetName(): 
            hist_mc[var] = hist_bkg[var][bkg].Clone()
        else: 
            hist_mc[var].Add(hist_bkg[var][bkg])
        hist_r[var] = hist_d[var].Clone()
        hist_r[var].Divide(hist_mc[var])
        hist_r[var].SetTitle("")
        hist_r[var].GetXaxis().SetTitle(xtitles[var])
        hist_r[var].GetYaxis().SetTitle("Data/#SigmaMC")
        hist_r[var].GetYaxis().SetNdivisions(506)
        hist_r[var].SetMarkerStyle(20)
        hist_r[var].SetMarkerSize(0.7)


# Make plot legend 

leg = R.TLegend(0.70,0.50,0.88,0.88)
leg.SetFillStyle(4000)  
leg.SetFillColor(0)
leg.SetTextFont(42)
leg.SetBorderSize(0)

bkg_labels = {'Zjets':'Z+jets', 'Top':'Top', 'Diboson':'Diboson', 'Wjets':'W+jets'}

sig_labels = {'Zprime2000':"Z' (2 TeV)"}

for bkg in backgrounds: 
    leg.AddEntry(hist_bkg['pt1'][bkg], bkg_labels[bkg], "f")

for sig in signals: 
    leg.AddEntry(hist_sig['pt1'][sig], sig_labels[sig], "f")

leg.AddEntry(hist_d['pt1'],"Data","ple")

selection = ""
if channel == "ee": 
    selection = "ee" 
if channel == "uu": 
    selection = "#mu#mu"

#Create directory for saving plots
if not os.path.exists('./Histograms'):
    os.makedirs('./Histograms')
if not os.path.exists('./Histograms/Plots/'):
    os.makedirs('./Histograms/Plots')


# Make plots
for var in variables: 

    cnv = R.TCanvas("cnv_"+var,"", 500, 500)
    cnv.SetTicks(1,1) 
    cnv.SetLeftMargin(0.13) 
    #cnv.SetLogy()

    p1 = R.TPad("p1", "", 0, 0.35, 1, 1) 
    p2 = R.TPad("p2", "", 0, 0.0, 1, 0.35) 

    p1.SetLogy()
    p1.SetBottomMargin(0.0) 
    p1.Draw() 
    p1.cd()

    stack[var].Draw("hist")
    stack[var].SetMinimum(10E-2)
    stack[var].GetYaxis().SetTitle("Events")
    stack[var].GetYaxis().SetTitleFont(43)
    stack[var].GetYaxis().SetTitleSize(16)
    stack[var].GetYaxis().SetLabelFont(43)
    stack[var].GetYaxis().SetLabelSize(16)
    stack[var].GetYaxis().SetTitleOffset(1.5)
    if var in ['eta1', 'eta2', 'phi1', 'phi2']: 
        maximum = stack[var].GetMaximum() 
        stack[var].SetMaximum(maximum*10E4)

    hist_d[var].Draw("same e0")
    leg.Draw("same")

    for sig in signals:
        hist_sig[var][sig].SetFillColor(0);
        hist_sig[var][sig].Draw("same hist");

    s = R.TLatex()
    s.SetNDC(1);
    s.SetTextAlign(13);
    s.SetTextColor(R.kBlack);
    s.SetTextSize(0.044);
    s.DrawLatex(0.4,0.86,"#font[72]{ATLAS} Open Data");
    s.DrawLatex(0.4,0.81,"#bf{#sqrt{s} = 13 TeV,^{}%.1f^{}fb^{-1}}" % (L));
    s.DrawLatex(0.4,0.76,"#bf{"+selection+" selection}");


    p1.Update() 
    p1.RedrawAxis() 

    cnv.cd() 

    p2.Draw() 
    p2.cd() 

    p2.SetGridy()

    hist_r[var].SetMaximum(1.99) 
    hist_r[var].SetMinimum(0.01) 	
    hist_r[var].Draw("0PZ") 

    p2.SetTopMargin(0) 
    p2.SetBottomMargin(0.35) 
    p2.Update()        
    p2.RedrawAxis() 

    cnv.cd() 
    cnv.Update()
    cnv.Print('./Histograms/Plots/'+channel+'_'+var+'.png') 
    cnv.Close() 
