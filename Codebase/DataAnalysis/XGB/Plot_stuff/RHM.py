# --------------------------------------------------------------#
# This is a python file used to create histograms in the style #
# of ATLAS ROOT.                                               #
# --------------------------------------------------------------#

import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import seaborn as sns


class ROOT_Histo_Maker:
    """
    *----------------------------------*The ROOT Histogram Maker (RHM)*---------------------------------*
    |   Args:                                                                                           |
    |        MC_Data         --              (2D - Array): An array where each element contains the     |
    |                                                      values of the variable for each event in     |
    |                                                      a given channel.                             |
    |        MC_Weights      --             (2D - Array): An array where each element contains the      |
    |                                                     distribution of weights for a given           |
    |                                                     channel.                                      |
    |        channel_labels  --             (1D - Array): An array containing the labels of each        |
    |                                                     channel.                                      |
    |        Data            --             (1D - Array): Array containing the value for each event     |
    |                                                     in the data.                                  |
    |        nr_bins         --             (int)       : Number of bins in histogram.                  |
    |        bin_max         --             (int)       : Largest value for edge in bins.               |
    |        bin_min         --             (int)       : Smallest value for edge in bins.              |
    |        variable_name   --             (string)    : Name of variable to be plotted.               |
    |        y_scale         --             (string)    : String to set scale of y-axis.                |
    |        y_max           --             (float)     : Float to decide y-axis max value.             |
    |        y_min           --             (float)     : Float to decide y-axis min value.             |
    |        default_pallet  --             (bool)      : Decide to use dedualt plotting settings       |
    |                                                     or not.                                       |
    |        saveAs          --             (string)    : Name to save figure as. If argument not given |
    |                                                     figure will not be saved.                     |
    |        show            --             (bool)      : Boolean to show or not show figure.           |
    |                                                     If alteration to the figure is wanted,        |
    |                                                     show should be set to "False".                |
    *---------------------------------------------------------------------------------------------------*
    """

    def __init__(
        self,
        MC_Data,
        MC_Weights,
        channel_labels,
        Data = None,
        nr_bins=25,
        bin_max=1000,
        bin_min=0,
        variable_name="",
        y_scale="log",
        y_max=None,
        y_min=None,
        default_pallet=True,
        saveAs=None,
        show=True,
    ):

        self.MC_Data = MC_Data
        self.MC_Weights = MC_Weights
        self.channel_labels = channel_labels
        self.Data = Data
        self.nr_bins = nr_bins
        self.bin_max = bin_max
        self.bin_min = bin_min
        self.variable_name = variable_name
        self.y_scale = y_scale
        self.y_max = y_max
        self.y_min = y_min
        self.default_pallet = default_pallet
        self.saveAs = saveAs
        self.show = show

        self.sortData()

        self.fig, (self.ax1, self.ax2) = plt.subplots(
            2,
            1,
            gridspec_kw={"height_ratios": [3, 1]},
            num=0,
            dpi=80,
            facecolor="w",
            edgecolor="k",
            figsize=(7.5, 5.8),
        )

        self.createHisto()

        if self.saveAs != None:
            plt.savefig(self.saveAs, bbox_inches='tight')
        if self.show:
            plt.show()

    def createHisto(self):
        if self.Data is not None:
            N, bins = np.histogram(
                self.Data, bins=self.nr_bins, range=(self.bin_min, self.bin_max)
            )
            x = (np.array(bins[0:-1]) + np.array(bins[1:])) / 2

        self.ax1.set_yscale(self.y_scale)
        self.ax1.set_ylabel("#Events", fontsize=16)
        self.ax1.tick_params(
            axis="x", which="both", bottom=False, top=False, labelbottom=False
        )
        n, bins, patches = self.ax1.hist(
            self.MC_Data,
            weights=self.MC_Weights,
            bins=self.nr_bins,
            range=(self.bin_min, self.bin_max),
            histtype="barstacked",
            stacked=True,
            label=self.channel_labels,
        )
        if self.Data is not None:
            self.ax1.scatter(x, N, c="black", label="Data", zorder=100)

        box1 = self.ax1.get_position()
        self.ax1.set_position([box1.x0, box1.y0, box1.width * 0.7, box1.height])
        self.ax1.legend(fontsize=12, loc='center left', bbox_to_anchor=(1, 0.62), fancybox=True, shadow=True)

        if self.y_max != None:
            self.ax1.set_ylim(top=self.y_max)
        else:
            self.ax1.set_ylim(top = np.max(n)*1e2)

        if self.y_min != None:
            self.ax1.set_ylim(bottom=self.y_min)
        """else:
            self.ax1.set_ylim(top = np.min(n[n>0])*1e-1)"""

        self.ax1.set_xlim([bins[0], bins[-1]])
        n = self.calcN(bins, self.MC_Weights, self.MC_Data)
        if self.Data is not None:
            self.ax2.scatter(x, N / n, c="k", alpha=1, s=20)
        self.ax2.axhline(1, linestyle="--", c="k", alpha=0.7, linewidth=1)
        self.ax2.set_xlabel(self.variable_name, fontsize=16)#, loc="right")
        self.ax2.set_xlim([bins[0], bins[-1]])
        self.ax2.set_ylabel("Data/MC", fontsize=16)
        self.ax2.set_ylim([0.0, 2])
        plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)

    def calcN(self, bins, weight, values):
        s = 0
        n = []
        for j in range(len(bins) - 1):
            s = 0
            for i in range(len(values)):
                s += np.sum(
                    weight[i][(values[i] < bins[j + 1]) * (values[i] > bins[j])]
                )
            n.append(s)
        return n

    def sortData(self):
        nrEvents = [np.sum(w) for w in self.MC_Weights]
        sort_indx = sorted(range(len(nrEvents)), key=lambda k: nrEvents[k])
        self.MC_Data = [self.MC_Data[i] for i in sort_indx]
        self.MC_Weights = [self.MC_Weights[i] for i in sort_indx]
        self.channel_labels = [self.channel_labels[i] for i in sort_indx]
        

