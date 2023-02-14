import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sys

sys.path.insert(1, "../../")
sys.path.insert(1, "../")
from Utilities import *
from plot_set import *

import numpy as np


def plotComp(metric, plotters =  None, name =  ""):
    if metric == "Auc":
        file = 'AUC'
    else:
        file = 'SIG'

    with open('../NN/results/SIG.json', 'r') as openfile:
        # Reading from json file
        json_object = json.load(openfile)

    fig, ax = plt.subplots(figsize=(10,8))

    methods = list(json_object.keys())

    m1 = []
    m2 = []

    method_0 = json_object[list(methods)[0]]
    map = {}
    for i in method_0:
        m1.append(int(method_0[i]["m1"]))
        m2.append(int(method_0[i]["m2"]))


    m1 = np.sort(np.unique(m1))
    m2 = np.sort(np.unique(m2))
    for i in range(len(m1)):
        for j in range(len(m2)):
            map[f"{m1[i]}_{m2[j]}"] = [i,j]

    for i in method_0:
        scores = []
        m1_i = method_0[i]["m1"]
        m2_i = method_0[i]["m2"]
        
        for method in methods:
            if plotters is not None and method not in plotters:
                continue
            score = json_object[method][f"{m1_i}_{m2_i}"]["score"]
            scores.append(score)

        scores = 100*np.asarray(scores)/np.max(scores)
        scores = np.around(scores)

        draw_pie(scores,  map[f"{m1_i}_{m2_i}"][0], map[f"{m1_i}_{m2_i}"][1], 1750, ax=ax)
    
    if plotters is not None:
        methods = [met for met in methods if met in plotters]

    colorList = [color_cycle(i+2) for i in range(len(methods))]
    method_label = [met[:-4] for met in methods]

    ax.set_xlim([-0.5, len(m1)-0.5])
    ax.set_ylim([-0.5, len(m2)-0.5])

    ax.set_xticks(np.arange(len(m1)) , minor=False)
    ax.set_yticks(np.arange(len(m2)) , minor=False)
    ax.set_xticklabels(m1,rotation=90, fontsize = 18)
    ax.set_yticklabels(m2, fontsize = 18)

    ax.legend(borderpad = 1.25, framealpha = 1,fontsize = 'xx-large',loc =  'lower left', labels = method_label,labelcolor = colorList,)
    ax.set_xlabel(r"$\tilde{\chi}_2$ [Gev]",fontsize =24, loc = "right")
    ax.set_ylabel(r"$\tilde{\chi}_1$ [Gev]",fontsize =24, loc = "top",rotation=0, labelpad = -40)
    fig.savefig(f"../../../thesis/Figures/MLResults/NN/SUSY/Comparison/{name}NetworkComp.pdf")





def draw_pie(dist, 
             xpos, 
             ypos, 
             size, 
             ax=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(10,8))

    # for incremental pie slices
    cumsum = np.cumsum(dist)
    cumsum = cumsum/ cumsum[-1]
    pie = [0] + cumsum.tolist()
    i = 2
    for r1, r2 in zip(pie[:-1], pie[1:]):
        angles = np.linspace(2 * np.pi * r1, 2 * np.pi * r2)
        x = [0] + np.cos(angles).tolist()
        y = [0] + np.sin(angles).tolist()

        xy = np.column_stack([x, y])
        ax.scatter([xpos], [ypos], marker=xy, s=size, color = color_cycle(i))
        i += 1
        
    ax.scatter([xpos], [ypos], s=size, facecolor = 'none', edgecolors=color_cycle(2+np.where(dist ==100)[0][0]), linewidth = 8)
    return ax

if __name__ == "__main__":

    metric = "Sig"
    plotters = ["StochChannelOutGrid", "ChannelOutGrid", "MaxOutGrid"]
    plotters = ["MaxOutGrid", "PNNGrid", "NNGrid"]
    plotters = ["MaxOutPCAGrid", "PNNPCAGrid", "NNPCAGrid"]
    name = "PCA"
    plotComp(metric, plotters = plotters, name = name)
    