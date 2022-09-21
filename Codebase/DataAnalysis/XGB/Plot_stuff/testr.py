import numpy as np
from Codebase.DataAnalysis.XGB.Plot_stuff.RHM import ROOT_Histo_Maker

mu, sigma = 0, 0.1 # mean and standard deviation

s = 0
mc_data = []
weight = []
for i in range(4):
    
    mc_data.append(np.random.poisson(5, 10000*int((i/2 + 1))).ravel())
    weight.append(np.ones(10000*int((i/2 + 1))))
    s += int(10000*(i/2 + 1))

data = np.random.poisson(5, s)


ROOT_Histo_Maker(mc_data, weight, ["Diboson", "ttbar" , "Higgs", "Z+jets"], data, bin_max=20, bin_min = 0,nr_bins=14, y_max=100000, y_min = 0.5, variable_name = r"$E_t[GeV]$")