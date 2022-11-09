from sklearn.metrics import auc
from sklearn.metrics import det_curve

# Plotting settings
import matplotlib.pyplot as plt

import numpy as np



def plotDET(Y, Y_pred, weights, title, return_score = False, name = None, plot = None ):
    fpr, fnr, thresholds = det_curve(Y,Y_pred, sample_weight = weights, pos_label=1)
    sorted_index = np.argsort(fpr)
    fpr =  np.array(fpr)[sorted_index]
    fnr = np.array(fnr)[sorted_index]
    
    fpr = fpr[np.where(fpr<=1)]
    fnr = fnr[np.where(fpr<=1)]
    
    lw = 2
    if plot:
        plt.figure()
        plt.plot(fpr, fnr, color='darkorange',
                lw=lw, label='ROC curve')
        plt.xlim([-0.005, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(title)
        plt.legend(loc="lower right")
        plt.savefig(f"{name}", bbox_inches="tight")
        plt.show()
 