import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sys

sys.path.insert(1, "../../")
sys.path.insert(1, "../")
from Utilities import *

import numpy as np


def plotComp(metric):
    if metric == "Auc":
        file = 'AUC'
    else:
        file = 'SIG'

    with open('../NN/results/SIG.json', 'r') as openfile:
        # Reading from json file
        json_object = json.load(openfile)

    fig, ax = plt.subplots(figsize=(10,8))

    methods = json_object.keys()

    m1 = []
    m2 = []

    method_0 = json_object[list(methods)[0]]

    for i in method_0:
        m1.append(int(method_0[i]["m1"]))
        m2.append(int(method_0[i]["m2"]))

    # m1 = np.sort(m1)
    # m2 = np.sort(m2)

    Scores = []
    M1 = []
    M2 = []
    for i in method_0:
        scores = []
        m1_i = method_0[i]["m1"]
        m2_i = method_0[i]["m2"]
        
        for method in methods:
            score = json_object[method][f"{m1_i}_{m2_i}"]["score"]
            scores.append(score)
            # print(score)
        scores = 100*np.asarray(scores)/np.max(scores)
        scores = np.around(scores)
        print(scores)

        # Scores.append(scores_i)
        # M1.append(int(m1_i))
        # M2.append(int(m2_i))

        # ax.scatter(x=int(m1_i), y=int(m2_i), s=0)
        # ax.pie(np.array(Scores), center = (M1, M2), radius = 1)
        print(m1, m2)
        print(scores)
        draw_pie(scores, m1_i, m2_i, 1500, ax=ax)

    fig.savefig("test.pdf")





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

    for r1, r2 in zip(pie[:-1], pie[1:]):
        angles = np.linspace(2 * np.pi * r1, 2 * np.pi * r2)
        x = [0] + np.cos(angles).tolist()
        y = [0] + np.sin(angles).tolist()

        xy = np.column_stack([x, y])
        ax.scatter([xpos], [ypos], marker=xy, s=size)
    return ax

if __name__ == "__main__":

    metric = "Sig"
    plotComp(metric)
    