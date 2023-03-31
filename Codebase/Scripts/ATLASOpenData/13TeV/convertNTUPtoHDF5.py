import uproot3 as uproot
import pandas as pd
import numpy as np
from os import listdir, walk
from os.path import isfile, join
import os
import h5py
import awkward0 as awkward
import time
from itertools import combinations
import uproot3_methods.classes.TLorentzVector
import sys


def initialize(indir):
    filedic = {}
    files = [f for f in listdir(indir) if (isfile(join(indir, f)) and f.endswith(".root"))]
    for f in files:
        ident = f.split(".")[1]
        cat = ""
        typ = ""
        for key in bkg_dsid_toplot.keys():
            if key == ident:
                cat = bkg_dsid_toplot[key]["cat"]
                typ = "bkg"
                break
        for key in sig_dsid_toplot.keys():
            if key == ident:
                cat = sig_dsid_toplot[key]["cat"]
                typ = "signal"
                break
        if "data_" in f:
            cat = "data"
            typ = "data"
        if cat:
            if not cat in filedic:
                filedic[cat] = {"files":[],"type":typ,"dsid":[]}
            filedic[cat]["files"].append(indir+"/"+f)
            filedic[cat]["dsid"].append(f.split(".")[0].split("_")[-1])
        else:
            print("WARNING \t Could not find category for %s"%ident)
    return filedic

def getSignalCategories():
    newdic = {}
    for key in sig_dsid_toplot.keys():
        if not sig_dsid_toplot[key]["cat"] in newdic.keys():
            newdic[sig_dsid_toplot[key]["cat"]] = 0
        newdic[sig_dsid_toplot[key]["cat"]] += 1
    print("#"*31)
    print("#### Background categories ####")
    print("#"*31)
    print("%-20s %10s"%("Category","N(samples)"))
    print("-"*31)
    for key in sorted(newdic.keys()):
        print("%-20s %10s"%(key,newdic[key]))
    return sorted(newdic.keys())

def getBkgCategories():
    newdic = {}
    for key in bkg_dsid_toplot.keys():
        if not bkg_dsid_toplot[key]["cat"] in newdic.keys():
            newdic[bkg_dsid_toplot[key]["cat"]] = 0
        newdic[bkg_dsid_toplot[key]["cat"]] += 1
    print("#"*31)
    print("#### Background categories ####")
    print("#"*31)
    print("%-20s %10s"%("Category","N(samples)"))
    print("-"*31)
    for key in sorted(newdic.keys()):
        print("%-20s %10s"%(key,newdic[key]))
    return sorted(newdic.keys())

def getSamplesInCategory(cat):
    newdic = {}
    for key in bkg_dsid_toplot.keys():
        if bkg_dsid_toplot[key]["cat"] == cat:
            newdic[bkg_dsid_toplot[key]["DSID"]] = key
    for key in sig_dsid_toplot.keys():
        if sig_dsid_toplot[key]["cat"] == cat:
            newdic[sig_dsid_toplot[key]["DSID"]] = key
    print("#"*31)
    print("####### Category %s #######" %cat)
    print("#"*31)
    print("%-20s %10s"%("DSID","Description"))
    print("-"*31)
    for key in sorted(newdic.keys()):
        print("%-20s %10s"%(key,newdic[key]))
    return newdic

def getMCCategory(filedic):
    MCcat = {}
    for cat in filedic:
        for dsid in filedic[cat]["dsid"]:
            MCcat[dsid] = cat
    return MCcat

def calc_sf(xsec,lumi,nev,mcWeight,scaleFactor_PILEUP,scaleFactor_ELE,scaleFactor_MUON,scaleFactor_BTAG,scaleFactor_LepTRIGGER):
    if lumi <= 0: 
        print("Lumi {:d} is not valid".format(lumi)) 
        return 0
    wgt = (mcWeight)*(scaleFactor_PILEUP)*(scaleFactor_ELE)*(scaleFactor_MUON)*(scaleFactor_BTAG)*(scaleFactor_LepTRIGGER)
    return wgt * ((xsec*1000*lumi)/nev)

def skimming(df,events,prev_n,next_n,nlep,lep_ptcut):
    #pt, met = events.arrays(['lep_pt','met_et'],outputtype=tuple,entrystart=prev_n,entrystop=next_n)
    # Vector to hold all the 
    
    pt  = awkward.fromiter(df['lep_pt'])
    met = awkward.fromiter(df['met_et'])
    
    skim = []
    
    # MET > 50 GeV
    skim.append(pd.Series(met > 50000,name='bools').astype(int))
    
    # Makes sure that each sub-vector has the same number of entries (filling -999 if not)
    pt = pt.pad(nlep).fillna(-999)
    
    # Require leptons to have enough pt (according to values set in lep_ptcut)
    for i in range(len(lep_ptcut)):
        se = pt > lep_ptcut[i]
        skim.append(pd.Series(pt[pt > lep_ptcut[i]].counts >= 1).astype(int))
        mask = np.logical_and(pt != pt.max(), pt.max != -999)
        pt = pt[mask]
        
     
    # Make sure we only have exactly nlep (after the pt cuts)
    skim.append(pd.Series(pt[pt > 0].counts == 0).astype(int))
    
        
    # < here one can add additional cuts. Remeber to append the result to the skim vector)
    
    # Make sure that all our entries in the skim vector has value 1 
    # If not this means that one of the cuts above did not pass (i.e. we don't want to keep the event)
    # Adding values from all skim vectors together should give a total equal to the length of the skim vector
    sk_final = skim[0]
    for i in range(1,len(skim)):
        sk_final = sk_final.add(skim[i])
    final_skim = pd.Series(sk_final == len(skim))
    
    # Keep only rows where we have right number of leptons with pT above thresholds
    df = df[final_skim.values]
    
    return df

def jetaugmentation(df,events,prev_n,next_n):
    #pt, eta, phi, e = events.arrays(['jet_pt','jet_eta','jet_phi','jet_E'],outputtype=tuple,entrystart=prev_n,entrystop=next_n)
    pt  = awkward.fromiter(df['jet_pt'])
    eta = awkward.fromiter(df['jet_eta'])
    phi = awkward.fromiter(df['jet_phi'])
    e   = awkward.fromiter(df['jet_E'])
    df['jet_n60'] = pt[pt > 60000].counts
    #df['jet_n60'] = np.count_nonzero(pt > 60000)
    return df

def get_invmass(value):
    #print("Type of value is ",type(value))
    try:
        return(value.mass)
    except:
        #print("Could not return mass for ", value)
        return np.nan #-1
#-----

def get_deltaPhi(x,y):
    #print("Type of value is ",type(x))
    try:
        return(x.delta_phi(y))
    except:
        print("Could not return deltaPhi for ", (x,y))
        return -1
#-----


def get_deltaR(x,y):
    #print(" of value is ",type(x))
    try:
        return(x.delta_r(y))
    except:
        print("Could not return deltaR for ", (x,y))
        return -1
#-----

def lepaugmentation(df,events,prev_n,next_n,nlep):
    #pt, eta, phi, e = events.arrays(['lep_pt','lep_eta','lep_phi','lep_E'],outputtype=tuple,entrystart=prev_n,entrystop=next_n)
    
    #print(df.keys())
    #subset = df[['lep_pt']]#,'lep_eta','lep_phi','lep_E']]
    #pt = tuple(subset.values)
    #pt = fromparents(df['lep_pt'])
    pt  = awkward.fromiter(df['lep_pt'])
    eta = awkward.fromiter(df['lep_eta'])
    phi = awkward.fromiter(df['lep_phi'])
    e   = awkward.fromiter(df['lep_E'])
    print("c")
    pt_org  = awkward.fromiter(df['lep_pt'])
    
    tlv = uproot3_methods.classes.TLorentzVector.TLorentzVectorArray.from_ptetaphim(pt, eta, phi, e)
    
    df["tlv"] = tlv[:]
    
    pt = pt.pad(nlep).fillna(-999)
    eta = eta.pad(nlep).fillna(-999)
    phi = phi.pad(nlep).fillna(-999)
    e = e.pad(nlep).fillna(-999)
    print("b")
    # Make the lepton variables
    for i in range(1,nlep+1):
        df['lep%i_pt'%i]  = pt[pt.argmax(),:].flatten()
        df['lep%i_eta'%i] = eta[pt.argmax(),:].flatten()
        df['lep%i_phi'%i] = phi[pt.argmax(),:].flatten()
        df['lep%i_E'%i]   = e[pt.argmax(),:].flatten()
        df['lep%i_tlv'%i]   = tlv[pt.argmax()].flatten()
        # Remove the hardest and continue to find the next to hardest etc.
        mask = np.logical_and(pt != pt.max(), pt.max != -999)
        pt  = pt[mask]
        eta = eta[mask]
        phi = phi[mask]
        e   = e[mask]
        tlv   =   tlv[mask]
    
    # Compute variables for all combinations of 2 leptons
    pairs = pt_org.argchoose(2)
    left  = pairs.i0
    right = pairs.i1

    #print(df.head())
    #print(df["lep1_tlv"].head())
    print("a")
    for ilep in range(len(left[0])):
        i = left[0][ilep]
        j = right[0][ilep]
        idx1 = left[0][i]
        idx2 = right[0][i]
        
        
        df['mll_%i%i'%(i+1,j+1)]   = (df['lep%i_tlv'%(i+1)]+df['lep%i_tlv'%(j+1)]).apply(get_invmass)
        df['dphi_%i%i'%(i+1,j+1)] = df.apply(lambda x : get_deltaPhi(x['lep%i_tlv'%(i+1)],x['lep%i_tlv'%(j+1)]), axis=1)
        df['dR_%i%i'%(i+1,j+1)]   = df.apply(lambda x : get_deltaR(x['lep%i_tlv'%(i+1)],x['lep%i_tlv'%(j+1)]), axis=1)
    
    dodrop = ['tlv']
    for i in range(1,nlep+1):
        dodrop.append('lep%i_tlv'%i)
    print("d")
    df = df.drop(dodrop,axis=1)
    df = df.select_dtypes(exclude=['int32'])
        
    return df

def createDataFrames(indir,nlep,chunksize,printevry,datatype,skimtag,lep_ptcut,branches,lumi,categories = []):
    # Count how many events and files we write/read
    n_allfiles   = 0
    nfile   = 0
    totev   = 0
    totev_skim = 0
    out_filenum = 1
    #try:
    #    del [result]
    #except:
    #    print("WARNING\t Result does not exists. Good :-)")
           
    root_files = initialize(indir)
 
    # Checking if file exitst (not needed)
    #path = indir+"/%s_%s.h5" %(datatype,skimtag)
    #if os.path.exists(path):
    #    os.remove(path)
    #    print("WARNING\t Removing file {:s}".format(path))
    for cat in root_files:
        
        if len(categories) and not cat in categories:
            print("INFO \t Skipping category %s"%cat)
            continue
        
        isSignal = False
        isBkg = False
        isData = False
        if root_files[cat]['type'] == 'signal':
            isSignal == True
        elif root_files[cat]['type'] == 'bkg':
            isBkg == True
        elif root_files[cat]['type'] == 'data':
            isData == True
        else:
            print("ERROR \t Could not find type of data set for category {:s}. Skipping".format(cat))
        
        for f in root_files[cat]["files"]:
            #if not "410012" in f: continue
            n_allfiles += 1
            # Getting the file and extracting the information from info dictionary
            
            events = uproot.open(f)["mini"]
            
            nentries = events.numentries
            
            print("INFO  \t Opening file {:d}/{:d}: {:s} with {:d} events".format(n_allfiles,len(root_files),f,nentries))
            
            #path = indir+"/%s_%s.h5" %(f.split(".root")[0],skimtag)
            
            path = indir+"/%s_%s_num_%i.h5" %(datatype,skimtag,out_filenum)
            
            
            n = 1
            prev_n = 0

            while True:

                # Measure time to read 
                if n == 1: start = time.time()   
                else:   
                    end = time.time()
                    dur = (end - start)
                    dur_sec = chunksize/dur
                    m, s = divmod((nentries-((n-1)*chunksize))/dur_sec, 60)
                    h, m = divmod(m, 60)
                    print("INFO  \t Event/sec  = {:.0f}. ETC = {:d}h{:02d}m{:02d}s".format(dur_sec,int(h),int(m),int(s)))
                    start = time.time()

                # Get the range of events to read
                next_n = n*chunksize if n*chunksize < nentries else nentries
                print("INFO  \t Reading entries {:d} - {:d} of {:d}. Total so far: {:d}".format(prev_n,next_n,nentries,totev_skim))
                #df = events.arrays(branches)#flatten=False,entrystart=prev_n,entrystop=next_n)
                df = events.pandas.df(branches,flatten=False,entrystart=prev_n,entrystop=next_n)

                df = df.astype({'channelNumber': 'float32'})     
                df = df.astype({'eventNumber': 'float32'})

                totev += len(df.index)

                if df.shape[0] > 0:
                    # If MC: get the scale factor corresponding to the total data luminosity (is set above). 
                    # If data: set scale factor to 1!
                    if not "data" in f:
                        #df['wgt'] = 
                        np.vectorize(calc_sf)(df.XSection,lumi,df.SumWeights,df.mcWeight,
                                                         df.scaleFactor_PILEUP,df.scaleFactor_ELE,
                                                         df.scaleFactor_MUON,df.scaleFactor_BTAG,
                                                         df.scaleFactor_LepTRIGGER)
                    else:
                        df['wgt'] = 1.0

                    # Define if it is signal or background (important for ML classification)
                    if isBkg or isData:
                        df['isSignal'] = 0
                    elif isSignal:
                        df['isSignal'] = 1
                    else:
                        df['isSignal'] = 0

                    df['MCType'] = cat

                    #df = df.astype({'MCType': 'string'})

                    ###########   
                    # Skimming
                    ###########  
                    df = skimming(df,events,prev_n,next_n,nlep,lep_ptcut)

                if df.shape[0] > 0: 

                    ##############
                    # Augumenting
                    ##############
                    df = jetaugmentation(df,events,prev_n,next_n)
                    df = lepaugmentation(df,events,prev_n,next_n,nlep)

                    ###########   
                    # Slimming
                    ###########
                    df = df.drop(['mcWeight','scaleFactor_PILEUP','scaleFactor_ELE','scaleFactor_MUON','scaleFactor_BTAG','scaleFactor_LepTRIGGER'],axis=1)
                    df = df.drop(['jet_n', 'jet_pt', 'jet_eta', 'jet_phi', 'jet_E', 'jet_jvt'],axis=1)
                    df = df.drop(['jet_trueflav', 'jet_truthMatched', 'jet_MV2c10', 'jet_pt_syst'],axis=1)
                    df = df.drop(['lep_n', 'lep_truthMatched', 'lep_trigMatched', 'lep_pt', 'lep_eta'],axis=1)
                    df = df.drop(['lep_phi', 'lep_E', 'lep_z0', 'lep_charge', 'lep_type', 'lep_isTightID'],axis=1)
                    df = df.drop(['lep_ptcone30', 'lep_etcone20', 'lep_trackd0pvunbiased'],axis=1)
                    df = df.drop(['lep_tracksigd0pvunbiased', 'lep_pt_syst'],axis=1)
                    # If first time, create result data frame. If not; concatenate
                    try: 
                        result = pd.concat([result,df],axis=0, ignore_index=True)
                    except:
                        result = df
                        print("WARNING\t Starting a new result panda")

                    print(result.shape)
                    totev_skim += len(df.index)
                 
                    # Delete the temporary data frame
                    del [df]

                if totev_skim > printevry:
                    #result = sortAndIndex(result,nlep)
                    path = indir+"/%s_%s_num_%i.h5" %(datatype,skimtag,out_filenum) 
                    print("INFO  \t Read {:d} events in {:d} files, for which {:d} ({:.2f}%) were written to {:s}"
                        .format(totev,nfile,totev_skim,(float(totev_skim)/float(totev))*100.,path))
                    result.to_hdf(path,key='mini',mode='w',format='table')
                    out_filenum += 1
                    totev = 0
                    totev_skim = 0
                    nfile = 0
                    del [result]

                # If read everything
                if n*chunksize > nentries: break

                # Update counters before continuing
                prev_n = n*chunksize + 1
                n += 1
                #break

            nfile += 1
            del [events]  
            break
    # Make sure we write the last events
    #result = sortAndIndex(result,nlep)
    path = indir+"/%s_%s_num_%i.h5" %(datatype,skimtag,out_filenum)
    result.to_hdf(path,'mini',mode='w',format='table')
    print("INFO  \t Read {:d} events in {:d} files, for which {:d} ({:.2f}%) were written to {:s}"
        .format(totev,nfile,totev_skim,(float(totev_skim)/float(totev))*100.,path))
    #del [result]
    return result

if os.path.isfile("../../../Input/Background_samples_13TeV.txt"):
    infofile = open("../../../Input/Background_samples_13TeV.txt", "r")
else:
    infofile = open(sys.path[1]+"/Input/Background_samples_13TeV.txt", "r")
bkg_dsid_toplot = {}
sig_dsid_toplot = {}
for line in infofile:
    if line.startswith("#"): continue
    try:
        fname,dsid,cat = line.split()
    except:
        continue
    bkg_dsid_toplot[fname] = {"cat":cat,"DSID":int(dsid)}
infofile.close()

if os.path.isfile("../../../Input/Signal_samples_13TeV.txt"):
    infofile = open("../../../Input/Signal_samples_13TeV.txt", "r")
else:
    infofile = open(sys.path[1]+"/Input/Signal_samples_13TeV.txt", "r")
for line in infofile:
    if line.startswith("#"): continue
    try:
        fname,dsid,cat = line.split()
    except:
        continue
    sig_dsid_toplot[fname] = {"cat":cat,"DSID":int(dsid)}
infofile.close()

res = createDataFrames("/scratch/eirikgr/openData_13TeV/2lep/MC/",2, 50000, 1000, "MC", "2L_pt25_25_met50", [25000,25000],
                       ['eventNumber','channelNumber',
                        'mcWeight','scaleFactor_PILEUP','scaleFactor_ELE','scaleFactor_MUON',
                        'scaleFactor_BTAG','scaleFactor_LepTRIGGER',
                        'met_*','jet_*','lep_*','XSection','SumWeights'], 10.6)

