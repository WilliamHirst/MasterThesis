#!/usr/bin/env python
from ROOT import TChain, TFile, gROOT
import os, sys
import shutil
import multiprocessing
from datetime import datetime
c = {}#TChain("tree_NoSys")

gROOT.SetBatch(True)

folders = sys.argv[1]
physproc = sys.argv[2]

numsplit = 0
if "_" in physproc:
    numsplit = int(physproc.split("_")[-1])

physproc = physproc.split("_")[0]

scaletodata = sys.argv[3]

sumentries = 0
procnum = 0
procstart = -1
if len(sys.argv) > 6:
    sumentries = int(sys.argv[4])
    procnum    = int(sys.argv[5])
    procstart  = int(sys.argv[6])
    print("sumentries", sumentries)
    print("procnum", procnum)
    print("procstart", procstart)

if scaletodata not in ["year18","year1516", "year17"]: 
    print("Please specify period to scale to (year18 or year1516 or year17)" )

is1516 = False
is17   = False
is18   = False
if scaletodata == "year1516": is1516 = True
elif scaletodata == "year17": is17 = True
elif scaletodata == "year18": is18 = True

#Create a Proof-lite instance:
#ROOT.TProof.Open("");
#tell the chain that we want to use PROOF
#chain->SetProof();
#And this will now use all your cores!
#chain->Process("MySelector.C+");

def worker(n):                                                                                                                                                                                                     
    if n == 0:                                                                                                                                                                                                     
        c[c.keys()[n]].Process("MySusySkimAnalysis.C++","")                                                                                                                                                             
    else:                                                                                                                                                                                                          
        c[c.keys()[n]].Process("MySusySkimAnalysis.C","")     

for folder in folders.split(","):                                                                                                                                                                    
    for fname in os.listdir(folder):
        #if "FAKES" in fname: continue
        
        if not fname.endswith(".root"): continue
        
        #if not fname.endswith("merged_processed.root"): continue
        if not (fname.startswith(physproc+"_merged")) and not (fname.startswith(physproc) and fname.endswith('_merged_processed.root')) : continue
        
        #if not "Zjets" in fname or "diboson" in fname: continue
        sample = physproc
        if "data" in fname: isData = True
        else: isData = False
        sample = ""
        for sp in fname.split("_"):
            if "merged" in sp or sp.isnumeric(): 
                break
            sample += "%s_"%sp
        sample = sample[:-1]
        #sample = fname.replace("_merged_processed.root","")
        if len(sample.split("_")): physproc = sample
        if not physproc == sample: continue
        #if len(sample.split("_")) > 1: 
        #    print "Skipping ", sample
        #    continue
        #sample = fname.split("_")[0]
        #if "410472_" in fname:
        #    sample = sample.replace("410472_","")
        #    sample += "_nominal"
        #print sample
        if not sample in c.keys():
            #c[sample] = TChain(("tree_NoSys") if not isData else "tree_NoSys") #"data")#
            c[sample] = TChain(("%s_NoSys"%sample) if not isData else sample)#"data")#sample)#"data") #"data")#
        c[sample].Add(folder+"/"+fname)
        print("Added %s/%s to chain for %s" %(folder,fname,sample))
        #break
print("Added %i samples" %len(c.keys()))
gROOT.ProcessLine(".L HistFitterTree.cxx+")
gROOT.ProcessLine(".L make2L2JTree.cxx+")
gROOT.ProcessLine(".L makeMYTree.cxx+")
gROOT.ProcessLine(".L GetFakeWeight.C+")
print("a")

num_of_leptons = 2
today = datetime.today()
datestring = today.strftime('%a_%b_%d_%Y')
date_folder = "/scratch/eirikgr/ANAoutput/%s_%02dL_NTUP_%s_%s" %(datestring,num_of_leptons,"data1516" if is1516 else ("data17" if is17 else "data18"), "REL22" if "REL22" in folders else "TEST")

if not os.path.exists(date_folder):
    print("INFO \t Creating output folder %s" %date_folder)
    os.makedirs(date_folder)
else:
    print("INFO \t Output folder already exists %s" %date_folder)

#sys.exit()

for sample in c.keys():
    print(sample)
    if sumentries and sumentries != c[sample].GetEntries():
        print("ERROR \t Number of events used for job splitting (%i) is not equal events in TTree (%i). Aborting." %(sumentries,c[sample].GetEntries()))
        sys.exit()
    if "data" in sample:
        if sumentries and procnum and procstart >= 0:
            print("%i  :  %i" %(procnum,procstart))
            c[sample].Process("MySusySkimAnalysis.C++","data %s %i %i %s"%(scaletodata,procnum,procstart,date_folder),procnum,procstart)
        else:
            c[sample].Process("MySusySkimAnalysis.C++","data %s %s"%(scaletodata,date_folder))
    else:
        if sumentries and procnum and procstart >= 0:
            c[sample].Process("MySusySkimAnalysis.C++","mc %s %i %i %s"%(scaletodata,procnum,procstart,date_folder),procnum,procstart)   
        else:
            c[sample].Process("MySusySkimAnalysis.C++","mc %s %s"%(scaletodata,date_folder))
# for dirname, dirnames, filenames in os.walk("/scratch3/eirikgr/SusySkimAnaNtuples/MM_2018_05_07/MC/SUSY2_Bkgs_mc16a/MM0705/merged/"):
#     for subdirname in dirnames:
#         dir = os.path.join(dirname, subdirname)
#         for fname in os.listdir(dir):
#             print dir+"/"+fname
#             c.Add(dir+"/"+fname)
#         break
# c.LoadTree(-1)
# jobs = []                                                                                                                                                                                                      
# for i in range(len(c.keys())):                                                                                                                                                                                 
#     if c[c.keys()[i]].GetEntries() == 0: continue                                                                                                                                                                        
#     p = multiprocessing.Process(target=worker, args=(i,))                                                                                                                                                      
#     jobs.append(p)                                                                                                                                                                                             
#     p.start()


#infof = open("infofile.txt",'w')
#infof.write(date_folder)
#infof.close()
