import sys
sys.path.insert(1, "../")
from Plot_stuff.ROOTPlot import *

sys.path.insert(1, "../../")
from Utilities import saveLoad



name = "MaxOut"
signal = "SUSY"

channels = saveLoad(f"{signal}{name}channels_test.npy")
# predict_data = saveLoad(f"{signal}{name}predict_data_test.npy").ravel()
predict_sorted = saveLoad(f"{signal}{name}predict_sorted_test.npy")
weights_sorted = saveLoad(f"{signal}{name}weights_sorted_test.npy")


PlotRootHisto(predict_sorted, 
              weights_sorted, 
              channels, 
              title = f"NN/SUSY/MLDist/{name}Dist", 
              xlabel = "XGB-Output", 
              bins = 30)
PlotRootHisto(predict_sorted, 
              weights_sorted, 
              channels, 
              title = f"NN/SUSY/MLDist/{name}Dist_C7", 
              xlabel = "XGB-Output", 
              bins = 30,
              CutOff = 0.7)
