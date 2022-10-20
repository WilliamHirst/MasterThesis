import sys
sys.path.insert(1, "../")
from Plot_stuff.ROOTPlot import *

sys.path.insert(1, "../../")
from Utilities import saveLoad

signal = "ttbar" 
channels = saveLoad("channels_test.npy")
predict_data = saveLoad("predict_data_test.npy").ravel()
predict_sorted = saveLoad("predict_sorted_test.npy")
weights_sorted = saveLoad("weights_sorted_test.npy")


PlotRootHisto(predict_sorted, 
              weights_sorted, 
              predict_data, 
              channels, 
              title = f"NN/{signal}SearchDist", 
              xlabel = "XGB-Output", 
              bins = 30)
PlotRootHisto(predict_sorted, 
              weights_sorted, 
              predict_data, 
              channels, 
              title = f"NN/{signal}SearchDist_C7", 
              xlabel = "XGB-Output", 
              bins = 30,
              CutOff = 0.7)
