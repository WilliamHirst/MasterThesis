import pty
import shlex
import os
import time
import subprocess
import sys
import shutil
from ROOT import TChain
from datetime import datetime
import platform

indir = sys.argv[1]
doRun = sys.argv[2].split(",") #  ttbar,Wjets,PHWW,topOther,triboson,lowMassDY,Vgamma,higgs,PHWZ,singleTop,PHZZ,Zjets,diboson

# all:
#ttbar,Wjets,PHWW,topOther,triboson,lowMassDY,Vgamma,higgs,PHWZ,singleTop,PHZZ,Zjets,diboson,
# main:
#ttbar,Wjets,topOther,triboson,lowMassDY,Vgamma,higgs,singleTop,Zjets,diboson

#rel21: higgs,Zmumu,Zee,Ztautau,ttbar,singleTop,QCD,DY
#rel22: diboson,higgs,Zmumu,Zee,Ztautau,ttbar,singleTop,QCD,DY
#data: data17,data18,data16,data15

#diboson,triboson,higgs,lowMassDY,ttbar,topOther,Vgamma,PHWW,Zjets,PHZZ,Wjets,PHWZ,singleTop

#ZMET
#Diboson,ttbar,Zeejets,Zmmjets,Zttjets,Wjets,singletop






#doRun  = ["triboson","lowMassDY","higgs","topOther"] #,"singleTop","data18","data15-16","data17"
#doRun += ["diboson","Zjets","ttbar","Wjets","singleTop"]

#doRun = ["Zjets"]

#doRun  = ["data15-16","data17","data18"]#,"data15-16","data17"]

#doRun  = ["singleTop"]#,"lowMassDY"]

rundic = {}
for r, d, f in os.walk(indir):
    for file in f:
        if "merged_processed" in file and file.endswith(".root"):
            tag = file.split("_")[0]
            if tag not in doRun: 
                print("Skipping %s" %tag)
                continue
            if tag not in rundic.keys():
                rundic[tag] = []
            print("Adding ",os.path.join(r, file))
            rundic[tag].append(os.path.join(r, file))

year = "year18"           
if "mc16e"  in indir: year = "year18"
elif "mc16cd" in indir: year = "year17"
elif "mc16d" in indir: year = "year17"
elif "mc16a"  in indir: year = "year1516"
#else: year = ""

orig_dir = os.getcwd()

# Taring all files needed for execution
cmd = "tar -zcf runpack.tar.gz ./MMinputfiles/ myscript.C runMySusySkimAnalysis.py MySusySkimAnalysis.C MySusySkimAnalysis.h GetFakeWeight.C GetFakeWeight.h MatrixMethod.cxx MatrixMethod.h HistFitterTree.cxx HistFitterTree.h make2L2JTree.cxx make2L2JTree.h makeMYTree.cxx makeMYTree.h CalcGenericMT2/"
os.system(cmd)

if "REL22" in indir: rel = "rel22"
else: rel = "rel21"

rel = "BLWGT"

#sys.exit()
nkey = 0
screenSessions = []
startentry = 0
# If only one file, let's split the file into distinct events
splitFile = False
numsplit = 5
#tag = list(rundic.keys())[0]
doSplit = 1
new_rundic = {}
sumentries = {}
if doSplit:#len(rundic.keys()) == 1:
    for tag in sorted(rundic.keys()):
        #print(tag,rundic[tag])
        splitFile = True
        if not "data" in tag:
            c = TChain("%s_NoSys"%tag)
        elif "data" in tag and ("2L2JNtuples" in indir or "ntuples2L2J" in indir):
            c = TChain("data")
        else:
            c = TChain("%s"%tag)
        for j in range(len(rundic[tag])):
            #print("Adding ",rundic[tag][j])
            #print("Events before ",c.GetEntries())
            c.Add(rundic[tag][j])
            if not tag in sumentries.keys():
                sumentries[tag] = 0
            sumentries[tag] = c.GetEntries()
            #print("After ",c.GetEntries())
        print("INFO \t Splitting %s into %i jobs. Total number of events : %i" %(tag,numsplit,sumentries[tag]))


        for i in range(1,numsplit+1):
            #if not tag+"_"+str(i) in new_rundic.keys():
            #    new_rundic[tag+"_"+str(i)] = []
            new_rundic[tag+"_"+str(i)] = rundic[tag]

num_of_leptons = 2
today = datetime.today()
datestring = today.strftime('%a_%b_%d_%Y')            
#            
#sys.exit()
oldtag = ""
for key in sorted(new_rundic.keys()):
    #print("Running %s" %key)

    i = 0
    if "_" in key:
        i = int(key.split("_")[-1])
        tag = key.split("_")[0]
    else:
        tag = key
    
    if "data18" in key: 
        year = "year18"
        #sumentries = 99249653
    elif "data17" in key: year = "year17"
    elif "data15" in key: year = "year1516"
    elif "data16" in key: year = "year1516"
    elif "data15-16" in key: year = "year1516"

    date_folder = "/scratch/eirikgr/ANAoutput/%s_%02dL_NTUP_%s_%s" %(datestring,num_of_leptons,year,"ZMET")
    if not os.path.exists(date_folder):
        print("INFO \t Creating output folder %s" %date_folder)
        os.makedirs(date_folder)
    else:
        print("INFO \t Output folder already exists %s" %date_folder)
        
    if not year:
        print("ERROR \t Could not find year for %s in file %s" %(key,new_rundic[key]))
        continue

    #for i in range(0,numsplit):
    perproc = 0

    # reset when new process
    if not tag == oldtag:
        oldtag = tag
        startentry = 0
    
    if splitFile and numsplit > 1:
        perproc = int(sumentries[tag]/numsplit)
        #if splitFile: print "Sample %s has %.0f entries. Each of the %i procs will take %i events. A total of %.0f" %(key,sumentries,numsplit,perproc,numsplit*perproc)
        if i == (numsplit):
            procnum = sumentries[tag]-startentry
            procstart = startentry
            print("Proc %i: # events %i, process [%i-%i]" %(i,sumentries[tag]-startentry,startentry,(sumentries[tag]-startentry)+startentry))
        else:
            procnum = perproc
            procstart = startentry
            print("Proc %i: # events %i, process [%i-%i]" %(i,perproc,startentry,startentry+perproc))
        startentry += perproc+1

    ########
    #### THIS IS WHERE JOBS ARE SCHEDULED INTO DIFFERENT PROCESSES
    #######

    # Path to temporary directory where job is executed
    path = "runmidl_%s_%s_%s_%i" %(rel,key,year,i)

    # Create it if it doesn't exist
    if os.path.exists(path) and os.path.isdir(path):
        shutil.rmtree(path)

    try:
        os.mkdir(path)
    except OSError:
        print ("ERROR \t Creation of the directory %s failed" % path)
        continue

    # Copy tarball to temp. directory
    cmd = "cp runpack.tar.gz ./%s/." %path
    os.system(cmd)

    # Go into temp directory
    new_dir  = os.chdir(path) # change the current working directory
    cwd = os.getcwd()

    # Untar tarball and make steering script
    runStr  = "tar -zxf runpack.tar.gz\n"
    if not 'hpc' in platform.node():
        runStr += "export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase\n"
        runStr += "source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh\n"
        runStr += "lsetup root\n"
    else:
        runStr += "module purge\n"
        runStr += "module load DataAnalysis/1.0.3-foss-2019b-Python-3.7.4\n"
    #runStr += "chmod +x ./runMySusySkimAnalysis.py\n"
    runStr += "chmod +x ./myscript.C\n"
    if perproc == 0:
        #runStr += "./runMySusySkimAnalysis.py %s %s %s" %(indir,key,year)
        runStr += "root -l -b -q \'./myscript.C+(\"%s\",\"%s\",\",%s\",\"%s\")\' &> log.out"%(indir,key,year,date_folder)
    else:
        #runStr += "./runMySusySkimAnalysis.py %s %s %s %i %i %i" %(indir,key,year,sumentries[tag],procnum,procstart)
        runStr += "root -l -b -q \'./myscript.C+(\"%s\",\"%s\",\"%s\",\"%s\",%i,%i,%i)\' &> log.out"%(indir,key,year,date_folder,sumentries[tag],procnum,procstart)
    print("Run command : %s"%runStr)
    #runStr += "source bashtest.sh %s" %key

    # Write setup commands to steer file
    runFile = open("runSH.sh","w")
    runFile.write(runStr)
    runFile.close()

    # Issue command
    cmd = "chmod +x ./runSH.sh"
    os.system(cmd)

    # prepare cmd that will start job in a screen session
    # Note: it's actually two nested screen sessions
    cmdString = "screen -dmS dummy%s%s%s%ixc screen -S %s%s%s%ixc ./runSH.sh" %(rel,key,year,i,rel,key,year,i)

    screenSessions.append("%s%s%s%ixc"%(rel,key,year,i))

    # allocate new pseudo-terminal and spawn screen sessions
    (master, slave) = pty.openpty()
    cmdArgs = shlex.split(cmdString)
    p = subprocess.Popen(cmdArgs, close_fds = False, shell=False,
            stdin=slave, stdout=slave, stderr=slave)


    os.chdir(orig_dir) # change the current working directory
    cwd = os.getcwd()

    nkey += 1

print("Followig %i sessions were started:" %len(screenSessions))
for ss in screenSessions:
    print(ss)





