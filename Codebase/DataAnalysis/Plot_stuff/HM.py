import matplotlib.pyplot as plt
from Plot_stuff.ROCM import *
import numpy as np


def HM(model, X, Y, W, columns):
    columns_s = columns[Y.to_numpy() == 1] 
    unique_c = columns_s.unique()


    print(unique_c)
    bkg = (Y==0).to_numpy()
    Z, map1, map2, M1, M2 = getGrid(unique_c)

    for c in unique_c:
        txt = c.split("p0")
        m1 = txt[0][21:24]
        m2 = txt[2]

        print(m1,m2)
        index_i = (columns == c).to_numpy() + bkg
        X_i = X[index_i]
        Y_i = Y[index_i]
        W_i = W[index_i]
        Z[map2[f"{m2}"], map1[f"{m1}"]] = plotRoc(  Y_i, 
                                                    model.predict(X_i, batch_size=8192), 
                                                    W_i,
                                                    "",
                                                    plot = False,
                                                    return_score = True)
    
    fig, ax = plt.subplots()
    cmap = plt.contourf(M1, M2, Z)

    fig.colorbar(cmap)
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
    m2 = np.sort(np.unique(m2))[::-1]
    Z = np.zeros([len(m2), len(m1)])
    map1 = {}
    map2 = {}
    for i in range(len(m1)):
        map1[f"{m1[i]}"] = i
    for i in range(len(m2)):
        map2[f"{m2[i]}"] = i
    return Z, map1, map2, m1, m2
        

        
        
