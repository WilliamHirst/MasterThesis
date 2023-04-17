import sys
sys.path.insert(1, "../../")
from Utilities import loadDf,scaleData, PCAData
from plot_set import *
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerTuple
from scipy.spatial import ConvexHull
from scipy import interpolate

import numpy as np




def plotPCA():
    myPath = "/storage/William_Sakarias/William_Data"
    df, _, _, channels = loadDf(myPath,notInc=["LRS", "filtch", "HNL", "MGPy8EGA"])
    channels = df["channel"]
    channelsU = channels.unique() 
    df = df.drop(columns=["wgt_SG", "channel"])
    df = scaleData(df)
    dfPCA = PCAData(df)
    alphas = np.linspace(1.0, .2, len(channelsU))
    skip = False
    skip2 = False
    
    fig, axs = plt.subplots(
            2,
            2,
            num=0,
            dpi=80,
            facecolor="w",
            edgecolor="k",
            figsize=(7.4, 7.4),
        )
    ax1 = axs[0,0]
    ax2 = axs[0,1]
    ax3 = axs[1,0]
    ax4 = axs[1,1]
    plotters = []
    legends = []
    j = -1
    for i in range(len(channelsU)):
        channel = channelsU[i]
        if "Z" in channel and "jets" in channel: 
            if skip:
                continue
            index = np.asarray(channels == "Zmmjets") + np.asarray(channels == "Zeejets") + np.asarray(channels == "Zttjets")
            skip = True
            channel = "Zjets"
        elif "diboson" in channel:
            if skip2:
                continue
            index = np.asarray(channels == "diboson2L") + np.asarray(channels == "diboson3L") + np.asarray(channels == "diboson4L")
            skip2 = True
            channel = "Diboson"
        else:
            index = channels==channel
        j += 1
        dfPCA_0i = dfPCA[index,-2]
        dfPCA_1i = dfPCA[index,-1]
        color = color_cycle[int(j*len(color_cycle)/8)]
        m_1, m_2 = np.mean(dfPCA_0i), np.mean(dfPCA_1i)
        s_1, s_2 = np.std(dfPCA_0i), np.std(dfPCA_1i)
        p_1 = dfPCA_0i[np.abs(dfPCA_0i - m_1) < 0.5*s_1 ]
        p_2 = dfPCA_0i[np.abs(dfPCA_1i - m_2) < 0.5*s_2 ]
        points = np.c_[p_1[:10],p_2[:10]]
        s1= ax1.scatter(points[:,0], points[:,1], color = color)

        interp_x, interp_y = getCurve(points)
        f2, = ax2.fill(interp_x, interp_y, '--', c=color, alpha=alphas[j], label = channel)
        points = np.c_[p_1[:100],p_2[:100]]
        interp_x, interp_y = getCurve(points)
        ax3.fill(interp_x, interp_y, '--', c=color, alpha=alphas[j])
        points = np.c_[p_1[:300],p_2[:300]]
        interp_x, interp_y = getCurve(points)
        ax4.fill(interp_x, interp_y, '--', c=color, alpha=alphas[j])
        
        plotters.append((s1,f2))

    ax2.set_yticks([])
    ax4.set_yticks([])
    ax1.set_xticks([])
    ax2.set_xticks([])
    ax1.grid(False)
    ax2.grid(False)
    ax3.grid(False)
    ax4.grid(False)
    ax4.set_xlabel("PCA Feature Second Last", loc =  'right',fontsize =20, labelpad = 15)
    ax1.set_ylabel("PCA Feature Last", loc =  'top',fontsize =20)
    plt.tight_layout()
    plt.savefig("PCAPlotLast.pdf")
    plt.show()




def getCurve(points):
    hull = ConvexHull(points)
    x_hull = np.append(points[hull.vertices,0],
                    points[hull.vertices,0][0])
    y_hull = np.append(points[hull.vertices,1],
                    points[hull.vertices,1][0])

    # interpolate
    dist = np.sqrt((x_hull[:-1] - x_hull[1:])**2 + (y_hull[:-1] - y_hull[1:])**2)
    dist_along = np.concatenate(([0], dist.cumsum()))
    spline, u = interpolate.splprep([x_hull, y_hull], 
                                    u=dist_along, s=0, per=1)
    interp_d = np.linspace(dist_along[0], dist_along[-1], 50)
    interp_x, interp_y = interpolate.splev(interp_d, spline)
    return interp_x, interp_y






if __name__ == "__main__":
    plotPCA()



