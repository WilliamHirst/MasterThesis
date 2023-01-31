import ROOT as R
import sys
sys.path.insert(1, "../../")
from Utilities import *
sys.path.insert(1, "../../DataPreparation/RDataFrameAna")
import plottingTool as pt
from pyHelperFunctions import *


def PlotRootHisto(MC, MC_wgt, Channels, title, xlabel, bins, Data = None, CutOff = None):

    R.EnableImplicitMT(200)

    R.gSystem.AddDynamicPath("-I/home/wohirst/MasterThesis/CodebasPreparation/RDataFrameAna")

    df = mergeToRoot(MC, MC_wgt, Channels, Data=Data, CutOff = CutOff)
    if CutOff is None:
        CutOff = 0

    bkgdic = {"Wjets":{"color":R.kYellow+2},
          "Zjets2":{"color":R.kBlue-7},
          "Zeejets":{"color":R.kBlue-7},
          "Zmmjets":{"color":R.kBlue-7},
          "Zttjets":{"color":R.kBlue-7},
          "diboson2L":{"color":R.kRed-7},
          "diboson3L":{"color":R.kBlue-7},
          "diboson4L":{"color":R.kGreen-2},
          "higgs":{"color":R.kYellow+2},
          "singletop":{"color":R.kOrange-2},
          "topOther":{"color":R.kSpring-9},
          "triboson":{"color":R.kYellow+2},
          "ttbar":{"color":R.kRed+7},
          #"data18":{"color":R.kBlack},
          #"data17":{"color":R.kBlack},
          #"data16":{"color":R.kBlack},
          #"data15":{"color":R.kBlack}
    }   
    sigdic = {"ttbarHNLfullLepMLp15":{"color":R.kBlack},
          "ttbarHNLfullLepMLp75":{"color":R.kBlack},
          "ttbarHNLfullLepMLm15":{"color":R.kBlack},
          "ttbarHNLfullLepMLm75":{"color":R.kBlack},
    }
    featdic = {"ML_Val" : {"xlabel": xlabel}}

        
    histo = {}
    allhisto = []


    def runANA(df, histo, allhisto):

        for k in df.keys():
            # HISTOGRAMS
            histo["ML_Val_%s"%k] = df[k].Histo1D(("ML_Val_%s"%k,"ML_Val_%s"%k,bins,CutOff,1),"ML_Val","wgt")
        
        for k in histo.keys():
            allhisto.append(histo[k])
        R.RDF.RunGraphs(allhisto)
        hfile = R.TFile("histograms.root","RECREATE")
        hfile.cd()

        writeHistsToFile(histo, True)

        return histo, df

    histo, df = runANA(df,histo, allhisto)



    toplot = []
    for bkg in bkgdic.keys():
        toplot.append(bkg)

    allFeats = histo.keys()
    featuresPlot = []
    for feat in allFeats:
        if feat[-5:] == "Wjets":
            featuresPlot.append(feat[:-6])
            
    for sig in sigdic.keys():
        toplot.append(sig)

    for feature in featuresPlot:
        if feature in featdic:
            xlabel = featdic[feature]["xlabel"]
        else:
            xlabel = feature
        p = pt.Plot(histo,feature,toplot,xtext = xlabel, mergeSig = False)
        p.can.SaveAs(f"../../../thesis/Figures/MLResults/{title}.pdf")
        p.can.Draw()
    

