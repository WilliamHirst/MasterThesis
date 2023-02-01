import sys
sys.path.insert(1, "../")
from Plot_stuff.ROOTPlot import *

sys.path.insert(1, "../../")
from Utilities import saveLoad



name = "MaxOut"
signal = "SUSY"

channels = saveLoad(f"results/{signal}{name}channels_test.npy")
# predict_data = saveLoad(f"{signal}{name}predict_data_test.npy").ravel()
predict_sorted = saveLoad(f"results/{signal}{name}predict_sorted_test.npy")
weights_sorted = saveLoad(f"results/{signal}{name}weights_sorted_test.npy")

sigdic = {"MGPy8EGA14N23LOC1N2WZ750p0p00p0p03L2L7": 0,
        "MGPy8EGA14N23LOC1N2WZ750p0p050p0p03L2L7": 0,
        "MGPy8EGA14N23LOC1N2WZ800p0p00p0p03L2L7": 0,
        "MGPy8EGA14N23LOC1N2WZ800p0p050p0p03L2L7": 0}

PlotRootHisto(predict_sorted, 
              weights_sorted, 
              channels, 
              title = f"NN/SUSY/MLDist/{name}Dist", 
              xlabel = "XGB-Output", 
              sigdic=sigdic,
              bins = 30)
PlotRootHisto(predict_sorted, 
              weights_sorted, 
              channels, 
              title = f"NN/SUSY/MLDist/{name}Dist_C7", 
              xlabel = "XGB-Output", 
              bins = 30,
              sigdic=sigdic,
              CutOff = 0.99)
