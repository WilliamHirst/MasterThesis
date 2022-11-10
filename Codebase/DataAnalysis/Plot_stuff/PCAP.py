import sys
sys.path.insert(1, "../../")
from Utilities import loadDf,scaleData, PCAData
from plot_set import *
import matplotlib.pyplot as plt



def plotPCA():
    myPath = "/storage/William_Sakarias/William_Data"

    df, y, df_data, channels = loadDf(myPath,notInc=["LRS", "filtch", "HNL"])
    channels = df["channel"]
    channelsU = channels.unique() 
    df = df.drop(columns=["wgt_SG", "channel"])
    df = scaleData(df)
    dfPCA = PCAData(df)

    
    plt.figure()
    for i in range(len(channelsU)):
        channel = channelsU[i]
        index = channels==channel
        plt.scatter(dfPCA[index,0], dfPCA[index,1])

    plt.savefig("PCAPlot.pdf")
    plt.show()




    






if __name__ == "__main__":
    plotPCA()



