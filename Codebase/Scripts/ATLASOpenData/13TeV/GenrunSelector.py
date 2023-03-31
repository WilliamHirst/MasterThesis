from ROOT import TChain, TFile, gROOT
import sys, os, time

t0 = time.time()

gROOT.ProcessLine("gErrorIgnoreLevel = kError")

h_file = sys.argv[1]
c_file = sys.argv[2]
leptype = sys.argv[3]  
datatype = sys.argv[4]

#if type == 'Data':
#        input_dir = '/storage/data/fys5555/ATLAS_opendata/2lep/Data/'
#elif type == 'MC':
#        input_dir = '/storage/data/fys5555/ATLAS_opendata/2lep/MC/'
input_dir = '/storage/data/fys5555/ATLAS_opendata/'+leptype + '/' + datatype +'/'
myChain = TChain('mini')

for filename in os.listdir(input_dir):
        if not '.root' in filename: continue
        print(filename)
        myChain.Add(input_dir+filename)
#if type == 'MC':
#        for filename in os.listdir(input_dir):
#                if not '.root' in filename: continue
#                print(filename)
#                myChain.Add(input_dir+filename)

if not os.path.exists('./Histograms'):
    os.makedirs('./Histograms')
if not os.path.exists('./Histograms/MC/'):
    os.makedirs('./Histograms/MC')
if not os.path.exists('./Histograms/Data/'):
    os.makedirs('./Histograms/Data')

entries = myChain.GetEntries()

print("-------------------------------------------")
#if type == 'Data':
#        print("Running on real data!")
#else:
#        print("Running on Monte Carlo!")
print("Number of events to process: %d" %entries)
print("-------------------------------------------")

os.system('cp ' + h_file + ' ./MySelector.h')
os.system('cp ' + c_file + ' ./MySelector.C')

#if type == 'Data':
myChain.Process('MySelector.C+',datatype)
#else:
#        myChain.Process('MySelector.C+', "MC")

t = int( time.time()-t0 )/60

print("Time spent: %d min" %t)

