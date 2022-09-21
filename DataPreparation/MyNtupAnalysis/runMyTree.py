from ROOT import TChain, gROOT
import sys
from os import walk

year = str(sys.argv[1])

c = TChain("myTree")
nfiles = 0
for (dirpath, dirnames, filenames) in walk("/storage/eirikgr/melana_output221020/ntuples_v15_08102021_syst/"):
    for d in dirnames:
        #print(d)
        if not d == year: continue
        for (dirpath_2, dirnames_2, filenames_2) in walk(dirpath+d):
            for dd in dirnames_2:
                for (dirpath_3, dirnames_3, filenames_3) in walk(dirpath_2+"/"+dd):
                    print("Adding %i files from %s"%(len(filenames_3),dirpath_3))
                    #print(dirnames_3)
                    for f in filenames_3:
                        if ".root" in f:
                            c.Add(dirpath_3+"/"+f)
                            nfiles += 1
print("Added total of %i files"%(nfiles))
#sys.exit()
#c = TChain("myTree")
#c.Add("/storage/eirikgr/melana_output221020/data/2017/user.jsabater.data17_13TeV.periodF_sleptons_v8syst_output.root/user.jsabater.21174959._00000*.root*")
#c.Add("/storage/eirikgr/user.jsabater.21174959._000003.output.root")
#c.Add("/storage/eirikgr/melana_output221020/data/*/user.jsabater.*/user.jsabater.*.root*")
gROOT.ProcessLine(".L makeMYTree.cxx+")
gROOT.ProcessLine(".L make2L0JTree.cxx+")
gROOT.ProcessLine(".L GetFakeWeight.C+")
c.Process("myTree.C++",year)
