from ROOT import TChain, TFile, gROOT
import sys, os, time
from optparse import OptionParser
import multiprocessing as mp

gROOT.ProcessLine("gErrorIgnoreLevel = kError")

def wait_completion(pids, sleep):
    #Wait until completion of one of the launched jobs
    while True:
        for pid in pids:
            if not pid.is_alive():
                try:
                    print( ">>> Completed: "+str(pid.name)+" ("+str(pid.pid)+")" )
                except AttributeError:
                    print( ">>> Completed: "+str(pid.pid) )
                pids.remove(pid)
                return
        #print( "Waiting for completion of jobs..." )
        time.sleep(sleep) # wait before retrying
    return


def parallel_exec(target, args_list, nproc, name=None, sleep = 0.01):

        print( "Maximum number of CPUs to be used: %d" %nproc )
        pids = []
        for a in args_list:
                if len(pids) >= nproc:
                        wait_completion(pids, sleep)
                p = mp.Process(target=target, args=(a,))        
                n = str(a) 
                print( ">>> Started: "+str(n) )
                p.name = n
                pids.append(p)
                p.start()
        wait_all(pids,sleep)
        return


def wait_all(pids, sleep):

    while len(pids) > 0:
        wait_completion(pids, sleep)
    print( "All jobs finished !" )
    return


def runFile(path):
        myChain = TChain('mini')
        myChain.Add(path)

        entries = myChain.GetEntries()
        print( "%d events to process for %s" %(entries, path) )

        option = ''
        if '/Data/' in path:
                print("INFO \t Running on real data with %i events"%(entries))
                period = path.split('/')[-1].split('data')[1].split('.')[0];
                option = 'Data'+period
        else:
                print("INFO \t Running on MC with %i events"%(entries))
                option = 'MC'
        arg = option+' '+runcardfile
        myChain.Process("MySelector.C+", arg)
        return


parser = OptionParser()
parser.add_option("--header", help="header file to use", default="")
parser.add_option("--source", help="source file to use", default="")
parser.add_option("--runcard", help="use runcard to specify cuts", default="")
parser.add_option("--doData", help = "run on data", default=0)
parser.add_option("--doMC", help="run on mc", default=0)
parser.add_option("--dataset", help="which data set to run on [1largeRjet1lep,1lep,1lep1tau,2lep,3lep,4lep,exactly2lep,GamGam]", default="3lep")
parser.add_option("--numCores", help="number of cores to run on", default=25)
(options, args) = parser.parse_args()

h_file = str(options.header)
c_file = str(options.source)
runcard = str(options.runcard)
doData = int(options.doData)
doMC = int(options.doMC)
dataset = str(options.dataset)
numCores = int(options.numCores)

input_dir = []
runcardfile = ""
rootpath = "/storage/shared/data/fys5555/ATLAS_opendata_split/"

if h_file:
        print("INFO \t Using uploaded header file")
        os.system('cp ' + h_file + ' ./MySelector.h')
if c_file:
        print("INFO \t Using uploaded source file")
        os.system('cp ' + c_file + ' ./MySelector.C')
if doData:
        print("INFO \t Will run on data")
        input_dir.append('%s/%s/Data/'%(rootpath,dataset))
if doMC:
        print("INFO \t Will run on MC")
        input_dir.append('%s/%s/MC/'%(rootpath,dataset))
if runcard:
        print("INFO \t Will use cuts specified in the uploaded runcard")
        runcardfile = "runcard.txt"
        os.system('cp ' + runcard + ' ./'+runcardfile)

# Check that input exists and prepare list of files to run over
input_list = [];
for idir in input_dir:
        if not os.path.isdir(idir):
                raise Exception("ERROR \t Could not find any data set for %s on %s"%(dataset,"MC" if "/MC/" in idir else "Data"))
        for filename in os.listdir(idir):
                if '.root' in filename:
                        input_list.append(idir+filename);

#Sort inputs from largest to smallest to get an early start on the largest files
input_list = sorted(input_list, key = lambda x: -os.stat(x).st_size)

if not os.path.exists('./Histograms'):
    os.makedirs('./Histograms')
if not os.path.exists('./Histograms/MC/'):
    os.makedirs('./Histograms/MC')
if not os.path.exists('./Histograms/Data/'):
    os.makedirs('./Histograms/Data')

#Compile the code if necessary
gROOT.ProcessLine(".L MySelector.C+");

parallel_exec(runFile,input_list,numCores)

#Merge data outputs
if doData:
        print( "Merging data outputs..." )
        sources = [];
        for f in os.listdir('Histograms/Data/'):
            if 'Data_' in f and '.root' in f:
                sources.append(f);
        cmd = 'hadd -f Histograms/Data/hist.Data.2016.root ';
        for s in sources:
            cmd += 'Histograms/Data/' + s + ' ';
        print( cmd )
        os.system(cmd);
