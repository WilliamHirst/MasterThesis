import sys
sys.path.insert(1, "../../")
from Utilities import loadDf,scaleData, PCAData
from plot_set import *
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import numpy as np




def plotPCA():
    myPath = "/storage/William_Sakarias/William_Data"
    df, y, df_data, channels = loadDf(myPath,notInc=["LRS", "filtch", "HNL"])
    channels = df["channel"]
    channelsU = channels.unique() 
    df = df.drop(columns=["wgt_SG", "channel"])
    df = scaleData(df)
    dfPCA = PCAData(df)
    colors = ["khaki", "paleturquoise", "indigo", "pink", "chocolate"]
    alphas = np.linspace(1.0, .2, len(channelsU))
    skip = False
    
    fig, (ax1, ax2) = plt.subplots(
            1,
            2,
            num=0,
            dpi=80,
            facecolor="w",
            edgecolor="k",
            figsize=(7.4, 4.8),
        )
    plotters = []
    legends = []
    for i in range(len(channelsU)):
        channel = channelsU[i]
        if "Z" in channel and "jets" in channel: 
            if skip:
                continue
            index = np.asarray(channels == "Zmmjets") + np.asarray(channels == "Zeejets") + np.asarray(channels == "Zttjets")
            skip = True
            channel = "Zjets"
        else:
            index = channels==channel

        dfPCA_0i = dfPCA[index,0]
        dfPCA_1i = dfPCA[index,1]
        print(channel)
        try:
            color = plt.rcParams["axes.prop_cycle"].by_key()["color"][i]
        except:
            color = colors[i-10]
        m_1, m_2 = np.mean(dfPCA_0i), np.mean(dfPCA_1i)
        s_1, s_2 = np.std(dfPCA_0i), np.std(dfPCA_1i)
        p_1 = dfPCA_0i[np.abs(dfPCA_0i - m_1) < s_1 ][:10]
        p_2 = dfPCA_0i[np.abs(dfPCA_1i - m_2) < s_2 ][:10]
        points = np.c_[p_1,p_2]
        
        s1= ax1.scatter(points[:,0], points[:,1], color = color)

        hull = ConvexHull(points)
        x_hull = np.append(points[hull.vertices,0],
                        points[hull.vertices,0][0])
        y_hull = np.append(points[hull.vertices,1],
                        points[hull.vertices,1][0])
        
        # interpolate
        from scipy import interpolate
        dist = np.sqrt((x_hull[:-1] - x_hull[1:])**2 + (y_hull[:-1] - y_hull[1:])**2)
        dist_along = np.concatenate(([0], dist.cumsum()))
        spline, u = interpolate.splprep([x_hull, y_hull], 
                                        u=dist_along, s=0, per=1)
        interp_d = np.linspace(dist_along[0], dist_along[-1], 50)
        interp_x, interp_y = interpolate.splev(interp_d, spline)
        # plot shape

        f2 = ax2.fill(interp_x, interp_y, '--', c=color, alpha=alphas[i], label = channel)
        ax2.set_yticks([])
        ax1.grid(False)
        ax2.grid(False)
        plotters.append((s1,f2))
        plotters.append(channel)

    plt.legend()
    plt.tight_layout()
    plt.savefig("PCAPlot.pdf")
    plt.show()




    






if __name__ == "__main__":
    plotPCA()



