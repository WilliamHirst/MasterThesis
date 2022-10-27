from sklearn.metrics import auc
from sklearn.metrics import roc_curve

# Plotting settings
import matplotlib.pyplot as plt

import numpy as np



def plotRoc(Y, Y_pred, weights, title, return_score = False, name = None, plot = None ):
    fpr, tpr, thresholds = roc_curve(Y,Y_pred, sample_weight = weights, pos_label=1)
    sorted_index = np.argsort(fpr)
    fpr =  np.array(fpr)[sorted_index]
    tpr = np.array(tpr)[sorted_index]
    
    fpr = fpr[np.where(fpr<=1)]
    tpr = tpr[np.where(fpr<=1)]
    
    roc_auc = auc(fpr,tpr)
    lw = 2
    if plot:
        plt.figure()
        plt.plot(fpr, tpr, color='darkorange',
                lw=lw, label='ROC curve (area = %0.2e)' % (roc_auc))
        plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
        plt.xlim([-0.005, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(title)
        plt.legend(loc="lower right")
        plt.savefig(f"{name}", bbox_inches="tight")
        plt.show()
    if return_score:
        return roc_auc
