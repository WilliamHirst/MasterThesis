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
    print(Z)
    print(M1, M2)


    min_val = 10000

    fig, ax = plt.subplots()
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
            score = Calc_Sig(predict_prob(X_i), Y_i, W_i, sf =sf, best_threshold=threshold,max_sig= np.max(predict_prob(X_i)))

            plt.text(map1[f"{m1}"]+0.5,map2[f"{m2}"]+0.5, f"{score:.3f}", color = "white") 
            colorBar = lambda Z: Z
            if score < min_val:
                min_val =  score
        if save:
            saveToJson(score, m1, m2, metric, method)

        if metric == "Auc":
            score = score*100 - 90

        Z[map2[f"{m2}"], map1[f"{m1}"]] = score


    # colorBar = lambda Z: Z
    # Z = [[0.,         0.,        1.1199161 ,  0.,         0.67562842, 0.5045912 ],
    #      [0.,         0.,        0.,          0.8597455,  0.68695831, 0.50660583],
    #      [0.,         0.,        1.07967766,  0.,         0.64853513, 0.4907103 ],
    #      [0.,         0.,        0.,          0.79560227, 0.43838745, 0.49555962],
    #      [0.,         2.01730255,0.96195477,  0.,         0.61578891, 0.47611924],
    #      [1.38402886, 0.,        0.,          0.69550025, 0.59030954, 0.45754442],
    #      [0.,         0.90806263,0.76242112,  0.,         0.53022775, 0.43298332],
    #      [0.,         0.,        0.,          0.5452668,  0.48655325, 0.40809008],
    #      [0.,         0.,        0.47394461,  0.,         0.40635828, 0.35986222]]

    # M1 = ['400', '450', '650', '700', '750', '800']
    # M2 = ['0', '50', '100', '150', '200', '250', '300', '350', '400']



    #Z = np.where(Z == 0, np.nan, Z)
    cmap = plt.pcolormesh(np.arange(len(M1)+1), np.arange(len(M2)+1), colorBar(Z), cmap = 'magma')

    cbar = fig.colorbar(cmap)
    cbar.ax.tick_params(size=0)
    cbar.set_ticks([])
    # Set ticks in center of cells
    ax.set_xticks(np.arange(len(M1)) + 0.5, minor=False)
    ax.set_yticks(np.arange(len(M2)) + 0.5, minor=False)

    ax.set_xticklabels(M1,rotation=90)
    ax.set_yticklabels(M2)
    print(M2)
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
        

        
        
if __name__ == "__main__":
    HM(0, 0, 0, 0, 0,0)