import sys
sys.path.insert(1, "../")
from Plot_stuff.ROOTPlot import *

sys.path.insert(1, "../../")
from Utilities import saveLoad



name = "xgb"
signal = "SUSY"

channels = saveLoad(f"../results/XGB/channels_{name}.npy")
predict_sorted = saveLoad(f"../results/XGB/predict_sorted_{name}.npy")
weights_sorted = saveLoad(f"../results/XGB/weights_sorted_{name}.npy")

sigdic = {"MGPy8EGA14N23LOC1N2WZ400p0p0250p0p03L2L7": 0,
          "MGPy8EGA14N23LOC1N2WZ650p0p0400p0p03L2L7": 0,
          "MGPy8EGA14N23LOC1N2WZ700p0p050p0p03L2L7": 0,
          "MGPy8EGA14N23LOC1N2WZ800p0p0400p0p03L2L7": 0}

PlotRootHisto(predict_sorted, 
              weights_sorted, 
              channels, 
              title = f"XGB/SUSY/MLDist/{name}Dist", 
              xlabel = f"XGB-Output", 
              sigdic=sigdic,
              bins = 30,
              noData = True)
PlotRootHisto(predict_sorted, 
              weights_sorted, 
              channels, 
              title = f"XGB/SUSY/MLDist/{name}Dist_C7", 
              xlabel = f"XGB-Output", 
              bins = 30,
              sigdic=sigdic,
              CutOff = 0.9,
              noData = True)
