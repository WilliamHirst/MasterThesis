import sys
sys.path.insert(1, "../")
from Plot_stuff.ROOTPlot import *

sys.path.insert(1, "../../")
from Utilities import saveLoad



name = "PNN"
signal = "SUSY"

channels = saveLoad(f"results/channels_{name}.npy")
# predict_data = saveLoad(f"predict_data_test.npy").ravel()
predict_sorted = saveLoad(f"results/predict_sorted_{name}.npy")
weights_sorted = saveLoad(f"results/weights_sorted_{name}.npy")


sigdic = {"MGPy8EGA14N23LOC1N2WZ250p050p03L2L7": 0,
        "MGPy8EGA14N23LOC1N2WZ200p0p0100p0p03L2L7": 0,
        "MGPy8EGA14N23LOC1N2WZ300p0p0200p0p03L2L7": 0,
        "MGPy8EGA14N23LOC1N2WZ500p0p0100p0p03L2L7": 0}

PlotRootHisto(predict_sorted, 
              weights_sorted, 
              channels, 
              title = f"NN/SUSY/MLDist/{name}Dist", 
              xlabel = f"{name}-Output", 
              sigdic=sigdic,
              bins = 30,
              noData = True)
PlotRootHisto(predict_sorted, 
              weights_sorted, 
              channels, 
              title = f"NN/SUSY/MLDist/{name}Dist_C7", 
              xlabel = f"{name}-Output", 
              bins = 30,
              sigdic=sigdic,
              CutOff = 0.99,
              noData = True)
