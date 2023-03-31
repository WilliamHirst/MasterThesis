from ROOT import TChain, TFile
import sys, os, time

t0 = time.time() 

arg1 = sys.argv[1]  

if arg1 == 'Data': 
        input_dir = '/storage/data/fys5555/ATLAS_opendata/2lep/Data/'
elif arg1 == 'MC': 
        input_dir = '/storage/data/fys5555/ATLAS_opendata/2lep/MC/'


myChain = TChain('mini') 


for filename in os.listdir(input_dir):
        if not '.root' in filename: continue 
        print(filename)  
        myChain.Add(input_dir+filename) 
if arg1 == 'MC': 
        for filename in os.listdir(input_dir):
                if not '.root' in filename: continue 
                print(filename)  
                myChain.Add(input_dir+filename)

if not os.path.exists('./Histograms'):
    os.makedirs('./Histograms')
if not os.path.exists('./Histograms/MC/'):
    os.makedirs('./Histograms/MC')
if not os.path.exists('./Histograms/Data/'):
    os.makedirs('./Histograms/Data')

entries = myChain.GetEntries() 

print("-------------------------------------------")
if arg1 == 'Data': 
        print("Running on real data!")
else: 
        print("Running on Monte Carlo!") 
print("Number of events to process: %d" %entries)
print("-------------------------------------------")

if arg1 == 'Data': 
        myChain.Process("MySelector.C+", "Data")
else: 
        myChain.Process("MySelector.C+", "MC") 

t = int( time.time()-t0 )/60  

print("Time spent: %d min" %t) 
