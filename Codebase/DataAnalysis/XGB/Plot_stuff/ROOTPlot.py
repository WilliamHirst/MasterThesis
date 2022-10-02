import ROOT as R
from os.path import isfile, isdir
import sys
sys.path.insert(1, "../../")
from Utilities import *
sys.path.insert(1, "../../DataPreparation/RDataFrameAna")
from samples import configure_samples
import plottingTool as pt
from pyHelperFunctions import *


def PlotRootHisto(MC, MC_wgt, Data, Channels, title, xlabel, xmin, xmax, bins):

    R.EnableImplicitMT(200)

    #R.gROOT.ProcessLine(".L helperFunctions.cxx+");
    R.gSystem.AddDynamicPath("-I/home/wohirst/MasterThesis/Codebase/DataPreparation/RDataFrameAna")
    #R.gInterpreter.Declare('#include "helperFunctions.h"') # Header with the definition of the myFilter function
    #R.gSystem.Load("helperFunctions_cxx.so") # Library with the myFilter function

    df = mergeToRoot(MC, MC_wgt, Data, Channels)

    bkgdic = {"Wjets":{"color":R.kYellow+2},
            "Zjets2":{"color":R.kBlue-7},
            "diboson2L":{"color":R.kRed-7},
            "diboson3L":{"color":R.kBlue-7},
            "diboson4L":{"color":R.kGreen-2},
            "higgs":{"color":R.kYellow+2},
            "singletop":{"color":R.kOrange-2},
            "topOther":{"color":R.kSpring-9},
            "triboson":{"color":R.kYellow+2},
            "ttbar":{"color":R.kRed+7},
            "data18":{"color":R.kBlack}
    }
    featdic = {"ML_Val" : {"xlabel": xlabel}}

        
    histo = {}
    allhisto = []


    def runANA(df, histo, allhisto):
        
        for k in df.keys():
            # HISTOGRAMS
            histo["ML_Val_%s"%k] = df[k].Histo1D(("ML_Val_%s"%k,"ML_Val_%s"%k,bins,xmin,xmax),"ML_Val","wgt")
        
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

    for feature in featuresPlot:
        if feature in featdic:
            xlabel = featdic[feature]["xlabel"]
        else:
            xlabel = feature
        p = pt.Plot(histo,feature,toplot,xtext = xlabel)
        p.can.SaveAs(f"../../../thesis/Figures/ML_Results/{title}.pdf")
        p.can.Draw()
    
