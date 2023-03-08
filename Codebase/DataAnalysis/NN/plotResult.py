import sys
sys.path.insert(1, "../")
from Plot_stuff.ROOTPlot import *

sys.path.insert(1, "../../")
from Utilities import saveLoad



name = "PNNV2"
signal = "SUSY"

channels = saveLoad(f"results/PNNDistTest/channels_{name}.npy")
predict_sorted = saveLoad(f"results/PNNDistTest/predict_sorted_{name}.npy")
weights_sorted = saveLoad(f"results/PNNDistTest/weights_sorted_{name}.npy")

sigdic = {"MGPy8EGA14N23LOC1N2WZ250p050p03L2L7": 0,
        "MGPy8EGA14N23LOC1N2WZ200p0p0100p0p03L2L7": 0,
        "MGPy8EGA14N23LOC1N2WZ300p0p0200p0p03L2L7": 0,
        "MGPy8EGA14N23LOC1N2WZ500p0p0100p0p03L2L7": 0}

PlotRootHisto(predict_sorted, 
              weights_sorted, 
              channels, 
              title = f"NN/SUSY/MLDist/PNNDistTest/{name}Dist", 
              xlabel = f"{name}(200,300)-Output", 
              sigdic=sigdic,
              bins = 30,
              noData = True)
PlotRootHisto(predict_sorted, 
              weights_sorted, 
              channels, 
              title = f"NN/SUSY/MLDist/PNNDistTest/{name}Dist_C7", 
              xlabel = f"{name}(200,300)-Output", 
              bins = 30,
              sigdic=sigdic,
              CutOff = 0.7,
              noData = True)
