import matplotlib.pyplot as plt
from Plot_stuff.ROCM import *
import matplotlib.colors as colors

import numpy as np


def HM(model, X, Y, W, columns):
    columns_s = columns[Y.to_numpy() == 1] 
    unique_c = columns_s.unique()


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

        W_i.loc[(Y_i==0).to_numpy()] *= np.sum(W_i[(Y_i==1).to_numpy()])/np.sum(W_i[(Y_i==0).to_numpy()])

        auc = (plotRoc( Y_i, 
                        model.predict(X_i, batch_size=8192), 
                        W_i,
                        "",
                        plot = False,
                        return_score = True))*100 - 90

        print(m1,m2, f"{(((auc+90)/100)):.3f}")
        plt.text(m1,m2, f"{((auc+90)/100):.3f}", color = "white")

        Z[map2[f"{m2}"], map1[f"{m1}"]] = auc

        if 100*auc < min_val:
            min_val =  auc
    Z = np.where(Z == 0, np.nan, Z)
    
    #print(Z)
    
    #print(min_val, np.max(Z))
    # min_val = 7.443651606392535

    M1 = ['400', '450', '650', '700', '750', '800']
    M2 = ['0', '100', '150', '200', '250', '300', '350', '400', '500']
    
    # Z = [[       np.nan,        np.nan ,9.8069044 ,        np.nan ,9.85493751 ,9.83655702],
    #     [       np.nan,        np.nan ,9.76849106,        np.nan ,9.84785144 ,9.8427822 ],
    #     [       np.nan,        np.nan,        np.nan ,9.81122076 ,9.83380099 ,9.8633117 ],
    #     [       np.nan ,8.96158436 ,9.72649121,        np.nan ,9.80999454 ,9.79416771],
    #     [7.32182614,        np.nan,        np.nan ,9.73709742 ,9.75160429 ,9.85670611],
    #     [       np.nan ,7.44365161 ,9.59068287,        np.nan ,9.72308134 ,9.77407459],
    #     [       np.nan,        np.nan,        np.nan ,9.56079261 ,9.70886932 ,9.72227695],
    #     [       np.nan,        np.nan ,9.11776657,        np.nan ,9.64152682 ,9.75732222],
    #     [       np.nan,        np.nan,        np.nan ,9.81110924 ,9.86912759 ,9.85820018]]


    # print(np.logspace(np.log10(min_val), np.log10(np.nanmax(Z)),100))

    cmap = plt.pcolormesh(M1, M2, 10**np.array(Z), cmap = 'magma')#, levels = np.logspace(np.log10(min_val), np.log10(np.nanmax(Z)),100)), norm = colors.LogNorm(),

    cbar = fig.colorbar(cmap)
    cbar.ax.tick_params(size=0)
    cbar.set_ticks([])
    plt.savefig(f"HM_test.pdf", bbox_inches="tight")
    plt.show()

def getGrid(col):
    m1 = []
    m2 = []
    for c in col:
        txt = c.split("p0")
        m1.append(txt[0][21:24])
        m2.append(txt[2])
    m1 = np.sort(np.unique(m1))
    m2 = np.sort(np.unique(m2))
    Z = np.zeros([len(m2), len(m1)])
    map1 = {}
    map2 = {}
    for i in range(len(m1)):
        map1[f"{m1[i]}"] = i
    for i in range(len(m2)):
        map2[f"{m2[i]}"] = i
    return Z, map1, map2, m1, m2
        

        
        
