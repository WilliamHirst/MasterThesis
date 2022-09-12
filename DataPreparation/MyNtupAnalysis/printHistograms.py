from ROOT import *
import os, sys
fil = TFile(sys.argv[1])
dirlist = fil.GetListOfKeys()
iter = dirlist.MakeIterator()
key = iter.Next()

utfil = open("histogramsinfiles.txt","w")
utfil_empty = open("histogramsinfiles_empty.txt","w")


while key:
    cl = gROOT.GetClass(key.GetClassName());
    if not cl.InheritsFrom("TH1"): 
        key = iter.Next()
        continue

    obj = key.ReadObj()
    #print key.GetNbytes()
    if obj.GetSumOfWeights() == 0:
        utfil_empty.write(key.GetName()+"\n")
    else:
        utfil.write(key.GetName()+"\n")
    obj.Delete()
    key = iter.Next()
utfil.close()
utfil_empty.close()
