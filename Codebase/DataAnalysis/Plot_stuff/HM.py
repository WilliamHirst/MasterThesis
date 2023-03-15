import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patheffects as pe
import re



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
    print(Z)
    print(M1, M2)

    for c in unique_c:
        m1, m2 = getMass(c)

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
     

        elif metric == "Sig":
            sf = np.sum(W_i[(Y_i==1).to_numpy()])/np.sum(W_i[(Y_i==0).to_numpy()])
            W_i.loc[(Y_i==0).to_numpy()] *= sf
            score = Calc_Sig(predict_prob(X_i), Y_i, W_i, sf =sf, best_threshold=threshold,max_sig= np.max(predict_prob(X_i)))

       
        if save:
            saveToJson(score, m1, m2, metric, method)

        if metric == "Auc":
            score = score*100 - 90

        Z[map2[f"{m2}"], map1[f"{m1}"]] = score

    gridPlotter(mlType, name, metric)

def gridPlotter(mlType, name, metric, file_name = "SIG", cut_off = 10):
    import matplotlib
    method = name.split('/')[-1]
    with open(f'../{mlType}/results/{file_name}.json', 'r') as openfile:
        # Reading from json file
        json_object = json.load(openfile)
    scores = json_object[method]
    M1 = []
    M2 = []
    for score in list(scores.keys()):
        M1.append(int(score.split("_")[0]))
        M2.append(int(score.split("_")[1]))
    Z,  M1, M2, map1, map2 = createMap(M1, M2)
    fig, ax = plt.subplots()
    ax.grid(visible = False)
    for score_key in list(scores.keys()):
        elem = scores[score_key]
        score = elem["score"]
        m1_i =elem["m1"]
        m2_i =elem["m2"]
        if score >cut_off:
            scoreString = f"{score:.0f}"
            scale = np.NAN
        else:
            scoreString = f"{score:.2f}"
            scale = 1
        plt.text(map1[f"{m1_i}"]+0.5,map2[f"{m2_i}"]+0.5, scoreString, ha='center', va='center', color = "white", fontsize = 'medium', path_effects=[pe.withStroke(linewidth=1, foreground="black")]) 
        Z[map2[f"{m2_i}"], map1[f"{m1_i}"]] = score*scale
    colorBar = lambda Z: Z
    cmap = matplotlib.cm.magma.copy()
    cmap.set_bad('white',1.)
    cmap.set_under(color= "#eeeeee")
    cmap = plt.pcolormesh(np.arange(len(M1)+1), np.arange(len(M2)+1), colorBar(Z), cmap = cmap,vmin=0.0000001, edgecolors = "whitesmoke", linewidth =0.1)
    cbar = fig.colorbar(cmap)
    cbar.ax.tick_params(size=0)
    cbar.set_ticks([])

    # Set ticks in center of cells
    ax.set_xticks(np.arange(len(M1)) + 0.5, minor=False)
    ax.set_yticks(np.arange(len(M2)) + 0.5, minor=False)
    
    ax.set_xlabel(r"$\tilde{\chi}_2$ [Gev]", x = 1, fontsize =20)
    ax.set_ylabel(r"$\tilde{\chi}_1$ [Gev]", y = 1, fontsize =20,rotation=0, labelpad = -20)

    ax.set_xticklabels(M1,rotation=90)
    ax.set_yticklabels(M2)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
    plt.savefig(f"../../../thesis/Figures/MLResults/{mlType}/SUSY/Grid/{name}{metric}.pdf", bbox_inches="tight")
    plt.show()

def getGrid(col):
    m1 = []
    m2 = []
    for c in col:
        m1_i, m2_i = getMass(c)
        m1.append(int(m1_i))
        m2.append(int(m2_i))
    Z,  m1, m2, map1, map2 = createMap(m1, m2)
    return Z, map1, map2, m1, m2

def createMap(m1, m2):
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
    return Z,  m1, m2, map1, map2

        

def getMass(string):
    elem = string.split("WZ")
    m1 = elem[1][0:3]
    if "p0p0" in string:
        elem = string.split("p0p0")
    else:
        elem = string.split("p0")
    m2 = elem[1]  
    return m1, m2
        
if __name__ == "__main__":
    from ROCM import plotRoc
    from plot_set import *
    # cmap.set_under(color='black')
    
    #gridPlotter(mlType = "NN", name ="Events", metric = "NrEvents", file_name="NrEvents", cut_off=1000)
    # names = ["ChannelOutGrid", "HybridPCALeakyGrid", "HybridPCAMaxOutGrid", "MaxOutGrid", "MaxOutPCAGrid", "NNGrid", "NNPCAGrid", "NNshallowGrid", "PNNGrid", "PNNPCAGrid", "StochChannelOutGrid"]
    names = ["MaxOutPCA_FS_MLMGrid", "MaxOutPCA_FSGrid", "NN_FS_MLMGrid", "NN_FSGrid", "PNNPCA_FS_MLMGrid", "PNNPCA_FSGrid"]
    #names = ["MaxOut_InterpolationGrid", "NN_InterpolationGrid", "NN_OneMass_InterpolationGrid", "NN_OneMass_Overfitting_InterpolationGrid", "NN_OneMass_Overfitting8_InterpolationGrid", "NN_OneMass_Overfitting10_InterpolationGrid", "NN_OneMass_Overfitting15_InterpolationGrid"]
    for name in names:
        gridPlotter(mlType = "NN", name =f"FS/{name}", metric = "Sig")
    #HM(0, 0, 0, 0, 0,0)