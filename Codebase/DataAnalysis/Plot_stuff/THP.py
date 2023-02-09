import matplotlib.pyplot as plt
import matplotlib.colors as colors


def THP(history, model, signal):
    """
    Training History Plotter
    """
    training =  history.history['auc']
    validation =  history.history['val_auc']

    plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')

    plt.xlabel(r"$Epochs$", fontsize=16)
    plt.ylabel(r"$AUC$", fontsize=16)
    plt.plot(training, label = "Training")
    plt.plot(validation, label = "Validation")
    plt.legend(fontsize = 16)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
    plt.savefig(f"../../../thesis/Figures/MLResults/NN/{signal}/History/{model}History.pdf", bbox_inches="tight")
