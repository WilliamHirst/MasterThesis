import ROOT
from datetime import datetime
import os
from os import listdir
from os.path import isfile, join

num_of_leptons = 2
today = datetime.today()
datestring = today.strftime('%a_%b_%d_%Y')#"Thu_Oct_08_2020"#
#print(datestring)
#sys.exit()

tomerge = ["HISTOGRAMS"]#,"2L2J","MVA"]

tohadd = {}
for y in ["data1516","data17","data18"]:
    date_folder = "/scratch/eirikgr/ANAoutput/%s_%02dL_NTUP_%s_TEST" %(datestring,num_of_leptons,y)
    if not os.path.isdir(date_folder):
        print("Folder %s does not exist"%date_folder)
        continue

    onlyfiles = [f for f in listdir(date_folder) if isfile(join(date_folder, f))]
    
    for of in onlyfiles:
        if of.endswith("_0_0.root"):
            print("File %s not to be added"%of)
            continue
        key = of.split("_")[0]
        if not key in tohadd.keys():
            tohadd[key] = {}
        ftype = of.split("_")[-3]
        if not ftype in tomerge:
            print("Skipping merge of %s since not scheduled"%ftype)
            continue
        if not ftype in tohadd[key].keys():
            tohadd[key][ftype] = []
        tohadd[key][ftype].append(join(date_folder, of))

    for key in tohadd.keys():
        print("#"*100)
        print("\n%s\n"%key)
        print("#"*100)
        for ftype in tohadd[key].keys():
            print("-"*100)
            print("\n%s\n"%ftype)
            print("-"*100)
            haddstr = "hadd -f "
            newname = date_folder+"/"+key+"_"+"_".join(tohadd[key][ftype][0].split("/")[-1].split("_")[-5:-2])+"_0_0.root"
            haddstr += "%s "%(newname)
            rmstr = "mv "
            nhadd = 0
            for ds in tohadd[key][ftype]:
                if not y in ds: continue
                haddstr += "%s "%(ds)
                rmstr += "%s "%ds
                nhadd += 1
            if not nhadd: continue
            print("\n\n"+haddstr)
            print("nhadd=",nhadd)
            #os.system(haddstr)

            rmstr += " %s/unhadded/." %date_folder
            #os.system(rmstr)
            #print("\n\n"+rmstr)



            

        
        
        #os.system(haddstr)
    

