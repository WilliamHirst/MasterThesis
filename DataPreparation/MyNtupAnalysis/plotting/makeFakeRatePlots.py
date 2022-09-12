import os, sys


inputdir = sys.argv[1]
analysis = sys.argv[2]
datakey  = sys.argv[3].split(",")
variable = sys.argv[4].split(",")
tag      = sys.argv[5]
triggers = sys.argv[6].split(",")
alldata  = int(sys.argv[7])
nosubmit = int(sys.argv[8])

if alldata and len(datakey) > 1:
    print "ERROR \t Specify only one datakey when running alldata!"
    sys.exit()

dataperiod = (inputdir.split("/")[-1]).split("_")[6]

hfreg    = ["2L21"]#,"2L40","2L41","2L42","2L43","2L44","2L45"]#"2L20",,"2L30"
convreg  = ["2L23"]#,"2L24","2L25","2L26","2L27"]
realreg  = ["2L02","2L03","2L04","2L05","2L06","2L07","2L08","2L09","2L10","2L11","2L12"]#["2L04"]#["2L02","2L03","2L04","2L05","2L06","2L07","2L08","2L09","2L10","2L11","2L12"]#,"2L09"]#["2L02","2L03","2L04","2L05","2L06","2L07","2L08","2L09"] #"2L02",
lightreg = ["F2L05"]#["F2L05"]#,"F2L02","F2L04","F2L05"]#,"2L49"] #"F2L02","F2L01","F2L02","F2L04",
#convreg  = ["F2L01","F2L02","F2L04","F2L05"]#,"2L49"] #"F2L02",

reg = hfreg + convreg + realreg + lightreg

uncert = []
uncert.append("")
'''
uncert.append("EL_EFF_ID_TOTAL_1down")            
uncert.append("EL_EFF_ID_TOTAL_1up")              
uncert.append("EL_EFF_Iso_TOTAL_1down")           
uncert.append("EL_EFF_Iso_TOTAL_1up")             
uncert.append("EL_EFF_Reco_TOTAL_1down")          
uncert.append("EL_EFF_Reco_TOTAL_1up")            
uncert.append("MUON_EFF_ISO_STAT_1down")          
uncert.append("MUON_EFF_ISO_STAT_1up")            
uncert.append("MUON_EFF_ISO_SYS_1down")           
uncert.append("MUON_EFF_ISO_SYS_1up")             
uncert.append("MUON_EFF_RECO_STAT_1down")         
uncert.append("MUON_EFF_RECO_STAT_1up")           
uncert.append("MUON_EFF_RECO_STAT_LOWPT_1down")   
uncert.append("MUON_EFF_RECO_STAT_LOWPT_1up")     
uncert.append("MUON_EFF_RECO_SYS_1down")          
uncert.append("MUON_EFF_RECO_SYS_1up")            
uncert.append("MUON_EFF_RECO_SYS_LOWPT_1down")    
uncert.append("MUON_EFF_RECO_SYS_LOWPT_1up")      
#uncert.append("EL_EFF_TriggerEff_TOTAL_1down")    
#uncert.append("EL_EFF_TriggerEff_TOTAL_1up")      
uncert.append("EL_EFF_Trigger_TOTAL_1down")       
uncert.append("EL_EFF_Trigger_TOTAL_1up")         
uncert.append("MUON_EFF_TrigStat_1down")          
uncert.append("MUON_EFF_TrigStat_1up")            
uncert.append("MUON_EFF_TrigSyst_1down")          
uncert.append("MUON_EFF_TrigSyst_1up")  
'''

if not nosubmit:
    try:
        os.remove("MMinput.root")
    except:
        print("Error while deleting file MMinput.root")

nr = 0
true = []
for r in reg:
    oldr = r
    nr += 1
    print "Doing region %i/%i" %(nr,len(reg))
    fakeorreal = "FAKE"
    for l in ["E","M"]:
        if r in hfreg:
            if l == "E" and "2L4" in r:
                ch = "EESS"
            elif l == "M" and "2L4" in r:
                ch = "MMSS"
            else:
                ch = "ALL"
            sub = True
            true = ["_trueREAL_"]
        if r in convreg:
            if l == "M": 
                continue
            else: 
                ch = "ALL"
                sub = False
                true = ["_"]
        if r in realreg:
            sub = False
            true = ["_trueREAL_","_"]
            fakeorreal = "REAL"
            if l == "M": 
                ch = "MMOS"
            elif l == "E": 
                ch = "EEOS"
        elif r in lightreg:
            r = r.replace("F","")
            if l == "M": 
                continue
            fakeorreal = "FAKE"
            sub = False
            ch = "ALL"
            true = ["_trueLightFlavorDecay_"]#,"_truePromptPhotonConversion_"] 
        #else:
        #    print "ERROR \t Could not find channel for %s" %reg
        #    continue
            
        for d in datakey:
            new_inputdir = inputdir
            if d not in inputdir:
                new_inputdir = inputdir.replace(dataperiod,d)
            for v in variable:
                if v == "pT_eta" and not fakeorreal == "REAL": continue
                if v == "TCF": sub = False
                for tr in true:
                    for t in triggers:
                        for notrig in ["","_eventTrig"]:#"_lepnotrig",
                            if not t and notrig == "_lepnotrig": continue
                            for unc in uncert:
                                if unc and l == "M" and not "MUON_" in unc: continue
                                if unc and l == "E" and not "EL_" in unc: continue
                                if not r == "2L04" and unc: continue
                                #notrig = ""
                                #if "mu" in t and "E" in l: 
                                #    notrig = "_lepnotrig" 
                                #if "el" in t and "M" in l: 
                                #    notrig = "_lepnotrig" 
                                if t:
                                    histname = "h_lep_%s_%s_%s%s_%s_%s%s%s%s%s_%s" %(v,"nL" if v == "TCF" else "nT",t,notrig,ch,l,tr,fakeorreal,r,"_"+unc if unc else "",analysis)
                                else:
                                    histname = "h_lep_%s_%s_%s_%s%s%s%s%s_%s" %(v,"nL" if v == "TCF" else "nT",ch,l,tr,fakeorreal,r,"_"+unc if unc else "",analysis)
                                basestr = "python plotMySusySkimAna.py --setBATCH 1 --doMC 1 --inputdir %s --histname %s --doSub %i  --datakey %s" %(new_inputdir,histname,1 if sub else 0, "data15-16" if d == "data1516" else d)
                                if alldata:
                                    basestr += " --allData 1"
                                if nosubmit: 
                                    print basestr
                                    #if "mu" in t and "M" in l: 
                                    #    print basestr.replace(t,"%s_lepnotrig"%t)
                                else: 
                                    os.system(basestr)
                                    #if "mu" in t and "M" in l: 
                                    #    os.system(basestr.replace(t,"%s_lepnotrig"%t))
        r = oldr

if not nosubmit:
    try:
        if alldata: 
            os.rename(r'MMinput.root',r'MMinput_%s_alldata.root' %tag)
        else:
            os.rename(r'MMinput.root',r'MMinput_%s.root' %tag)
    except:
        print "Could not find any file with name MMinput.root"

#if alldata: sys.exit()

#if not nosubmit:
#    sys.exit()

#sys.exit()

try:
    os.remove("MMinput_frac.root")
except:
    print("Error while deleting file MMinput_frac.root")

alldata = False
for d in datakey:
    new_inputdir = inputdir
    if d not in inputdir:
        new_inputdir = inputdir.replace(dataperiod,d)
    for i in range(2,14):
        for ch in ["EESS","EMSS","EEOS","EMOS","MMOS"]:
            for v in ["metsig"]:#,"njet"]: #"njet","mT2",
                if not alldata:
                    basestr = "python plotMySusySkimAna.py --setBATCH 1 --doMC 1 --inputdir %s --histname h_lep_%s_TCF_nL_%s_%s_trueFAKE_REAL2L%02d_%s --doSub 0 --datakey %s" %(new_inputdir,v,ch,"M" if "MM" in ch else "E" ,i,analysis, "data15-16" if d == "data1516" else d)
                else:
                    basestr = "python plotMySusySkimAna.py --setBATCH 1 --doMC 1 --inputdir %s --histname h_lep_%s_TCF_nL_%s_%s_trueFAKE_REAL2L%02d_%s --doSub 0 --allData 1" %(new_inputdir,v,ch,"M" if "MM" in ch else "E" ,i,analysis)
                if nosubmit: print basestr
                else: os.system(basestr)
try:
    os.rename(r'MMinput_frac.root',r'MMinput_frac_%s.root' %tag)
except:
    print "Could not find any file with name MMinput_frac.root"
