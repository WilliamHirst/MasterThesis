import matplotlib.pyplot as plt
import matplotlib.colors as colors

import sys
sys.path.insert(1, "../")
from Plot_stuff.plot_set import *

def THP(history, model, signal):
    """
    Training History Plotter
    """
    training =  history.history['auc']
    validation =  history.history['val_auc']
    model_name = model.split("_")[0]

    fig, ax = plt.subplots(figsize=(8,6))
    plt.plot(training, label = "Training", c = color_cycle[-1])
    plt.plot(validation, label = "Validation", c = color_cycle[6])
    ax.set_ylabel("AUC", fontsize=18)
    ax.set_xlabel("Epochs", fontsize=18)
    ax.grid(color = "whitesmoke", linestyle ='-', linewidth =1.5)
    ax.set_axisbelow(True)
    ax.set_title(f"{model_name}", fontsize = 18, loc = "right")
    plt.legend(fontsize = 18)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
    plt.savefig(f"../../../thesis/Figures/MLResults/NN/{signal}/History/{model}History.pdf", bbox_inches="tight")


if __name__ == "__main__":
    import numpy as np
    class H():
        def __init__(self):
            self.history = dict(auc = np.sin(np.linspace(-2*np.pi, 2*np.pi, 20)), val_auc = np.cos(np.linspace(-2*np.pi, 2*np.pi, 20)))
    history = H()
    THP(history, 'Test', 'SUSY')
