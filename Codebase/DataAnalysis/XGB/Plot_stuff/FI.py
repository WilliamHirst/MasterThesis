import matplotlib.pyplot as plt
import numpy as np

def plotFI(model, features, signal):
    FI = model.feature_importances_
    sort_indx = FI.argsort()
    FI = FI[sort_indx]
    features = features[sort_indx]
    nrFeats = 10

    plt.clf()

    fig, (ax) = plt.subplots(
            1,
            1,
            num=0,
            dpi=80,
            facecolor="w",
            edgecolor="k",
            figsize=(7.5, 5.8),
        )
    ax.barh(features[-nrFeats:], FI[-nrFeats:], align='center')
    plt.yticks(np.arange(nrFeats),
           [rf"${featdic[feature]}$" for feature in features[-nrFeats:]])
    ax.set_xlabel(r"$Feature Importance$", fontsize=16)
    others = [rf"${featdic[features[:nrFeats][i]]}$" + f":{FI[i]:.3f}" for i in reversed(range(len(features[:nrFeats])))]
    others.insert(0,"Feature importance:")
    textstr = '\n'.join(others)
    props = dict(boxstyle='round', 
                 facecolor = u'#eeeeee',
                 edgecolor =  u'#bcbcbc',
                 alpha=0.9)
    ax.text(1.05, 0.9 , textstr, transform=ax.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)


    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)

    plt.savefig(f"../../../thesis/Figures/MLResults/XGB/{signal}FeatureImportance.pdf", bbox_inches='tight')
    plt.show()
    



featdic = {"lep1_Pt"  : r"P_{t}(l_{1}) ",
           "lep2_Pt"  : r"P_{t}(l_{2}) ",
           "lep3_Pt"  : r"P_{t}(l_{3}) ",
           "lep1_E"  : r"E(l_{1}) ",
           "lep2_E"  : r"E(l_{2}) ",
           "lep3_E"  : r"E(l_{3}) ",
           "lep1_Mt"  : r"M_{T}(l_{1}) ",
           "lep2_Mt"  : r"M_{T}(l_{2}) ",
           "lep3_Mt"  : r"M_{T}(l_{3}) ",
           "lep1_Eta" : r"\eta(l_{1})",
           "lep2_Eta" : r"\eta(l_{2})",
           "lep3_Eta" : r"\eta(l_{3})",
           "lep1_Phi" : r"\phi(l_{1})",
           "lep2_Phi" : r"\phi(l_{2})",
           "lep3_Phi" : r"\phi(l_{3})",
           "lep1_Charge" :  "Charge(l_{1})",
           "lep2_Charge" :  "Charge(l_{2})",
           "lep3_Charge" :  "Charge(l_{3})",
           "lep1_Flavor" :  "Flavor(l_{1})",
           "lep2_Flavor" :  "Flavor(l_{2})",
           "lep3_Flavor" :  "Flavor(l_{3})",
           "met_Et"   : r"E_{T}^{miss}",
           "met_Phi"  : r"\phi (miss)",
           "deltaR"  : r"\Delta R",
           "mlll"  : r"M_{lll} ",
           "njet_SG" : r"Nr \ signal \ Jets",
           "met_Sign"  : r"S(E_{t}^{miss})",
           "flcomp"  : r"Flavour Comb",
           "nbjet85" : r"NrB-jets (85)",
           "nbjet77" : r"NrB-jets (77)",
           "mll_OSSF"  : r"M_{ll} (OSSF) ",
           "Ht_lll"  : r"H_{t}(lll)",
           "Ht_SS"  : r"H_{t}(SS)",
           "Ht_met_Et"  : r"H_{t}(lll) + E_{T}^{miss}",
           "M_jj"  : r"M_{jj}",
           "lep1_Pt_Neg" :  r"Negative Weights abs(P_{t}(l_{1})) ",
           "lep2_Pt_Neg" :  r"Negative Weights abs(P_{t}(l_{2})) ",
           "lep1_Phi_Neg" :  r"Negative Weights abs(\phi(l_{1}))",
           "lep2_Phi_Neg" :  r"Negative Weights abs(\phi(l_{2}))"}
