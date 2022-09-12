import ROOT as R
import numpy as np

#df = R.RDataFrame("CollectionTree","/storage/PHYSLITE/D2AOD_PHYSLITE.fromPHYS.pool.root")
df = R.RDataFrame("CollectionTree","/storage/PHYSLITE/DAOD_PHYSLITE.fromAOD.pool.root")

def makeAlias(df):
    goodcols = R.vector('string')()
    cols = df.GetColumnNames()
    for c in cols:
        colname = str(c)
        if "xAOD::" in str(df.GetColumnType(c)) or "BASE" in str(df.GetColumnType(c)) or "string" in str(df.GetColumnType(c)):# or "vector" in str(df.GetColumnType(c)):
            print("INFO \t %s not a valid column with type %s"%(colname,df.GetColumnType(c)))
            continue
        if "char" in str(df.GetColumnType(c)):
            print("INFO \t %s not a valid column with type %s"%(colname,df.GetColumnType(c)))
            continue
        if not "AuxDyn" in colname:
            #print(c)
            #blacklist.push_back(c)
            continue
        if "Link" in colname or "detDescrTags" in colname or "originalTruthParticle" in colname:
            continue
        variable = colname
        colname = colname.replace("AuxDyn","")
        sp = colname.split(".")
        alias = sp[0]+"_"+sp[1]
        if "Analysis" in alias:
            alias = alias.replace("Analysis","")
        #print("Type for %s is %s"%(colname,type(c)))
        try:
            df = df.Alias(alias,variable)
            goodcols.push_back(alias)
        except:
            print("Could not make alias %s for %s"%(alias,variable))
            #blacklist.push_back(c)
            continue

    print("Good columns")
    for gc in goodcols:
        #print(gc)
        #if "xAOD::" in str(df.GetColumnType(gc)):
        try:
            print("Name = %-100s, Type = %-50s"%(str(gc),df.GetColumnType(gc)))
        except:
            print("Problems with %s"%str(gc))
    #print(type(goodcols))
    print("Makes a snapshot")
    return df.Snapshot("CollectionTree","out.root",goodcols)


#R.gSystem.AddDynamicPath("/home/eirikgr/myNtupAnalysis/jupyter-notebooks/")
#R.gInterpreter.Declare('#include "Cfunctions.h"') # Header with the definition of the myFilter function
#R.gSystem.Load("Cfunctions.so") # Library with the myFilter function

#PHYSdf = makeAlias(df)




df = R.RDataFrame("CollectionTree","out.root")

pt = df.AsNumpy(columns=["Jets_phi"])
#AODdf  = makeAlias(df2)

#cols = AODdf.GetColumnNames()
#for c in cols:
#    print(c)

#myHist1 = PHYSdf.Histo1D({"PHYSdf", "PHYSdf", 1000, 0., 1000}, "myColumn");

    
#cols = R.vector('string')()
#cols.push_back('TrigMatch_HLT_e160_lhvloose_nod0')
#p = df.Display(cols)
#p.Print()
