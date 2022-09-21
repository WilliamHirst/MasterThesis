from ROOT import TChain
import sys
from os import listdir, walk
from os.path import isfile, join

c = TChain("mini")

indir = sys.argv[1];

files = [f for f in listdir(indir) if (isfile(join(indir, f)) and f.endswith(".root"))]

nf = 0
for f in files:
    if "_new" in f: continue
    if not "mc_361106.Zee" in f: continue
    c.Add(indir+"/"+f)
    nf += 1

print("Added %i files"%nf)
c.Process("OpenDataLooper.C++")
