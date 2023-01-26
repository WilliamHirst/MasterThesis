import matplotlib.pyplot as plt
from Plot_stuff.ROCM import *
import matplotlib.colors as colors
import sys
sys.path.insert(1, "../../")
from Utilities import *

import numpy as np


def HM(model, X, Y, W, columns,name, metric = "Auc", data = None, save = False, mlType="NN"):
    columns_s = columns[Y.to_numpy() == 1] 
    unique_c = columns_s.unique()
    threshold = 0.998

    if mlType == "NN":
        predict_prob = lambda X: model.predict(X, batch_size=8192)
    elif mlType == "XGB":
        predict_prob = lambda X: model.predict_proba(X)[:,1]


    method = name.split('/')[-1]
    print(method)
    if save:
        EmptyJson(metric, method)


    bkg = (Y==0).to_numpy()
    Z, map1, map2, M1, M2 = getGrid(unique_c)

    min_val = 10000

    fig, _ = plt.subplots()
    for c in unique_c:
        txt = c.split("p0")
        m1 = txt[0][21:24]
        m2 = txt[2]

        print(m1,m2)
        index_i = (columns == c).to_numpy() + bkg
        X_i = X[index_i].copy()
        Y_i = Y[index_i].copy()
        W_i = W[index_i].copy()



        if metric == "Auc":
            W_i.loc[(Y_i==0).to_numpy()] *= np.sum(W_i[(Y_i==1).to_numpy()])/np.sum(W_i[(Y_i==0).to_numpy()])
            score = (plotRoc( Y_i, 
                            predict_prob(X_i), 
                            W_i,
                            "",
                            plot = False,
                            return_score = True))
            colorBar = lambda Z: 10**np.array(Z)
            if 100*score < min_val:
                min_val =  score

        elif metric == "Sig":
            sf = np.sum(W_i[(Y_i==1).to_numpy()])/np.sum(W_i[(Y_i==0).to_numpy()])
            W_i.loc[(Y_i==0).to_numpy()] *= sf
            score = Calc_Sig(predict_prob(X_i), Y_i, W_i, sf =sf, best_threshold=threshold)
            plt.text(m1,m2, f"{score:.3f}", color = "white") 
            colorBar = lambda Z: Z
            if score < min_val:
                min_val =  score
        if save:
            saveToJson(score, m1, m2, metric, method)

        if metric == "AUC":
            score = score*100 - 90

        Z[map2[f"{m2}"], map1[f"{m1}"]] = score



    Z = np.where(Z == 0, np.nan, Z)
    cmap = plt.pcolormesh(M1, M2, colorBar(Z), cmap = 'magma')#, levels = np.logspace(np.log10(min_val), np.log10(np.nanmax(Z)),100)), norm = colors.LogNorm(),

    cbar = fig.colorbar(cmap)
    cbar.ax.tick_params(size=0)
    cbar.set_ticks([])
    plt.savefig(f"../../../thesis/Figures/MLResults/{mlType}/{name}{metric}.pdf", bbox_inches="tight")
    plt.show()

def getGrid(col):
    m1 = []
    m2 = []
    for c in col:
        txt = c.split("p0")
        m1.append(int(txt[0][21:24]))
        m2.append(int(txt[2]))
    m1 = np.sort(np.unique(m1))
    m2 = np.sort(np.unique(m2))
    m1 = [f"{m}" for m in m1]
    m2 = [f"{m}" for m in m2]
    Z = np.zeros([len(m2), len(m1)])
    map1 = {}
    map2 = {}
    for i in range(len(m1)):
        map1[f"{m1[i]}"] = i
    for i in range(len(m2)):
        map2[f"{m2[i]}"] = i
    return Z, map1, map2, m1, m2
        

        
        
