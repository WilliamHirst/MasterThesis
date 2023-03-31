#!/usr/bin/env python

'''
#########################################
# b.k.gjelsten@fys.uio.no
# Aug 16, 2010: first version
#########################################
# PROGRAM DESCRIPTION


# KNOWN BUGS / "BUGS"


# FEATURE REQUESTS


#########################################
'''

import os,sys
from libExpValues import exp
from kilelib import WriteToFile, GetLinuxDist
from lineup import lineup

from mytestarg_lib import * 


import os
######################################### METHODS

######################################### METHODS
class Empty: pass

def doHELP():
    print("  Usage:  darksusy.py --plha < <list of lha-files piped in>")
    print("  Usage:  cat list_of_lha_files | darksusy.py --plha")
    print("  Ex:  echo susy_DG_emt_TB10_m1m2mu_140_140_140.lha| darksusy.py --plha -vb 1")
    print("  Ex:  ls *.lha| darksusy.py --plha -vb 1")
    print("  Options:")
    print("     -cdmcoann 0/1      (0:without coann, 1:with coann [can be slow])") 
    print("     --cdmwithcoann     (only without)") 
    print("     --cdmwithoutcoann  (does both)") 
    print("     --keeprundir       (keeps ./run_ds/ for inspection)")
    print()
    print("  Ex how to run standalone on list of slha and still get M1,M2,mu into the file:")
    print("    for fnslha in `ls susy_*.slha` ; do vars=`echo ${fnslha} | sed s/'.slha'//| sed s/susy_DGnoSL_TB10_M1M2MU_0// | sed s/_/' '/g` ; ls ${fnslha} | darksusy.py  --plha  -vb 0  -cdmcoann 1  -repl .slha,''  -prevar 'M1 M2 mu'  -preval \"\"\"${vars}\"\"\" ; done   # cd ~/grids_lsp/DGnoSL_TB10_2013-09-16/lhawig_M1eq050   # 2014-01-30")


accel_excl_txt = ['C1','gl','sq','sl','widZ','h','N1','b->sy','rho']
def accel_excl(accel, mode=2):
    if mode==0:
        txt = ""
        for i in range(9):
            if accel & (2**i): txt += " 1"
            else: txt += " 0"
        txt = txt[1:]
        return txt

    if mode == 1: 
        txt = []
        for i in range(9):
            if accel & (2**i): txt.append(accel_excl_txt[i])
        return txt

    if mode == 2: 
        txt = ""
        for i in range(9):
            if accel & (2**i): txt += "  "+accel_excl_txt[i]
        txt = txt[2:]
        return txt


def do_darksusy():

    
    iS = g.iS
    g.res.append({})
    res = g.res[iS]
    
    # set Steer file (or use default, but ensure '4')
    # <not implemented>
    if userATLASsteer:
        f = open('./ATLASsteer','w')
        f.write('4  #modchoice (4: auto)\n')
        for L in userATLASsteer: f.write('%s\n' %(L))   
        f.close()
        if VB>=4: print('user set ATLASsteer: %s' %(userATLASsteer)) 
    else: 
        cmd = 'cp -p %s ./ATLASsteer' %(DARKSUSY_STEER)
        if VB>=4: print(cmd)
        os.system(cmd)

    # set lha file
    if VB>=4: print('fn_lha :  ', fn_lha[iS])
    if VB>=4: print('pwd    :  ', os.getcwd())
    cmd = 'cp -p %s/%s ./ATLASlha' %(DIR_START, fn_lha[iS])
    if VB>=4: print(cmd)
    os.system(cmd)
    #print 'hei'
    

    # execute
    cmd = "%s > %s" %(DARKSUSY_EXE, fn_dsmain_out)
    if VB>=4: print(cmd)
    os.system(cmd)


    # ------------------------------------
    # Some minimal initialisation
    res['ATLAS_res0'] = ""
    res['accel'] = 0
    res['ATLAS_excl_txt'] = ""
    # ------------------------------------
    
    # read log & fill table
    f = open(fn_dsmain_out); lines = f.readlines(); f.close()
    for iL in range(len(lines)):
        L = lines[iL].strip()
        w = L.split()

        # A - Pick results from a summary line
        if L.startswith('ATLAS_res0'):

            #print L
            
            res['ATLAS_res0'] = L
            res['tot']       = int(w[3])
            res['unphys']    = int(w[6])
            res['accel']     = int(w[9])

            try: res['cdm_noco']  = float(w[12])  # Getting ****** for some (extreme) scenarios, 
            except: res['cdm_noco'] = -2.         #   ex: (M1,M2,MU) = (0,210,999)
            try: res['cdm_co']    = float(w[13])  # Eventually look into
            except: res['cdm_co'] = -2.

            res['cdm_best'] = res['cdm_co']
            if res['cdm_best'] < 0 and res['cdm_noco'] >= 0: res['cdm_best'] = res['cdm_noco']

            # VS 1sigma errors
            for zzz in ['cdm_noco', 'cdm_co', 'cdm_best']:
                if res[zzz] >= 0: res[zzz+'_VS1sig'] = (res[zzz] - exp['cdm']) / exp['cdm_error']
                else: res[zzz+'_VS1sig'] = -3.

            # VS EXP value
            for zzz in ['cdm_noco', 'cdm_co', 'cdm_best']:
                if res[zzz] >= 0: res[zzz+'_VSexp'] = res[zzz] / exp['cdm']
                else: res[zzz+'_VSexp'] = -3.

            
            res['(g-2)/2']   = float(w[16])  
            res['(g-2)/2_txt']   = w[16]
            res['amuSUSYadd_txt']   = w[16]
            res['amuSUSYadd']   = float(w[16])
            res['amuVS1sig'] = res['amuSUSYadd'] / exp['amuErr']
            # Note: what is given by darksusy is the addition to amu=(gmu-2)/2 from SUSY
            # Hence it is to be compared to the exp error of amu measurement

            # res[''] = w[1]
            # res[''] = w[1]
            # res[''] = w[1]

        # B - Pick some more results
        if L.startswith('sum_excl:'):
            res['excl_sum'] = int(w[1])
            res['excl_str'] = ""
            res['excl_str2'] = ""  # same but with one space in between each digit
            for i in range(9):
                if res['excl_sum'] & (2**i): res['excl_str'] += '1'
                else: res['excl_str'] += '0'
                if res['excl_sum'] & (2**i): res['excl_str2'] += '1 '
                else: res['excl_str2'] += '0 '
            res['excl_str2'] = res['excl_str2'].strip()

        # C - and some more
        if L.startswith('ATLAS_excl:'):
            zz = ''
            for i in range(len(w)):
                if i>1: zz += " "+w[i]
            if res['ATLAS_excl_txt'] == "": res['ATLAS_excl_txt'] = zz.strip()
            else: res['ATLAS_excl_txt'] += " , "+zz.strip()
            
        

    if VB>=4: print('res.keys(): ',list(res.keys()))
    res['accel_txt'] = accel_excl(res['accel'], mode=2)
    
    # Old 1-liner
    if VB>=2: print("CDM = %6.3f  %6.3f" %(res['cdm_noco'], res['cdm_co']))
    #if VB>=1: print "%s | %s | %s" %(fn_tag[iS], res['ATLAS_res0'], res['accel_txt'])
    outds0 = "[DS] %s | %s " %(fn_tag[iS], 'CGQLZhNbrho')
    outds  = "[DS] %s | %s    | %s | %s" %(fn_tag[iS], res['excl_str'], res['ATLAS_res0'].replace('ATLAS_res0 |','').strip() , res['ATLAS_excl_txt'])
    #print "%s | %s    | %3i | %s" %(fn_tag[iS], res['excl_str'], res['excl_sum'], res['ATLAS_excl_txt'])
    if VB>=1:
        print(outds0)
        print(outds) 

    # New 2-liner
    zhead = preVar + '  DMwoCo  DMwCo  DMbest  DMwoCoVSexp  DMwCoVSexp  DMbestVSexp  DMwoCoVS1sig  DMwCoVS1sig  DMbestVS1sig  amuSUSYadd  amuVS1sig  IsBad  UnPhys  AccBound  CGQLZhNbr C G Q L Z h N bsg rho'
    #zout  = preVal + '  %7.4f  %7.4f  %7.4f  %8.3f  %8.3f  %8.3f  %8.2f  %8.2f  %8.2f  %s  %8.3f  %i  %i  %i  %s  %s' %(res['cdm_noco'], res['cdm_co'], res['cdm_best'], res['cdm_noco_VSexp'], res['cdm_co_VSexp'], res['cdm_best_VSexp'], res['cdm_noco_VS1sig'], res['cdm_co_VS1sig'], res['cdm_best_VS1sig'], res['amuSUSYadd_txt'], res['amuVS1sig'], res['tot'], res['unphys'], res['accel'], res['excl_str'], res['excl_str2'])
    zout  = preVal + '  %7.4f  %7.4f  %7.4f  %9.3e  %11.5e  %9.3e  %8.2f  %8.2f  %8.2f  %s  %8.3f  %i  %i  %i  %s  %s' %(res['cdm_noco'], res['cdm_co'], res['cdm_best'], res['cdm_noco_VSexp'], res['cdm_co_VSexp'], res['cdm_best_VSexp'], res['cdm_noco_VS1sig'], res['cdm_co_VS1sig'], res['cdm_best_VS1sig'], res['amuSUSYadd_txt'], res['amuVS1sig'], res['tot'], res['unphys'], res['accel'], res['excl_str'], res['excl_str2'])
    outs = lineup([zhead,zout])
    WriteToFile(fn=DIR_START+'/'+fn_tag[iS]+'_ds%s.tlt' %(fnadd), outs=outs, VB=VB-1)


    if SAVEds0: f = open(DIR_START+'/'+fn_tag[iS]+fnadd+'.ds0','w'); f.write('%s\n' %(outds0)); f.write('%s\n' %(outds)); f.close()
    if SAVEds: f = open(DIR_START+'/'+fn_tag[iS]+fnadd+'.ds','w');  f.write('%s\n' %(outds));  f.close()
    cmd = 'cp -p %s %s' %(fn_dsmain_out, DIR_START+'/'+fn_tag[iS]+'.dslong')
    if SAVEdslong: os.system(cmd)

    
    #g.outs.append(res['ATLAS_res0'])

    # store output somewhere
    # py-res, log?, lha2-files?, dsmain1.tmp?,



def Main():
    #print "Executing Main()"


    if os.path.exists(DIR_RUN) and DIR_RUN != DIR_RUN0FULL:
        print("WARNING::darksusy.py  DIR_RUN expect not to exist, but does: %s" %(DIR_RUN))
        

    if not os.path.exists(DIR_RUN):
        if VB > 1: print("darksusy.py: INFO Creating dir %s" %(DIR_RUN))
        os.mkdir(DIR_RUN)
        

    # Loop
    for g.iS in range(g.nlha): 
        os.chdir(DIR_RUN)
        do_darksusy()
        os.chdir(DIR_START)
    

    # Summarise?
    if(VB>=3): print("darksusy.py: INFO: %i scenarios done" %(g.nlha))

    # Delete directory
    if '*' not in DIR_RUN and not KEEP_DIR_RUN: # just to avoid disaster if so should try to occur
        os.system('rm -r %s' %(DIR_RUN))
        if VB > 1: print('darksusy.py: INFO Removing dir %s' %(DIR_RUN))


######################################### GLOBAL VARS & INITIALISATION
g = Empty()
g.iS = 0
g.res = []

VB = 1

fn_lha0 = []
fn_lha = []
fn_wig = []
fn_out = []
fn_tag = []
fn_tagReplace = []
fn_tagWith = []

DIR_TMP = "/tmp"

dist = GetLinuxDist()  # 2014-08-25  (added for RHEL6)

if not 'SUSYPHENO_PATH' in list(os.environ.keys()):
    print("$SUSYPHENO_PATH not set. Exiting.")
    sys.exit()
susyphenopath = os.environ['SUSYPHENO_PATH']
DIR_DARKSUSY = susyphenopath+"/DARKSUSY/darksusy-5.1.1/test"  # RHEL5
#if   dist in ['RHEL5']: DIR_DARKSUSY = "/net/abel-evs/cargo/fysepf/epfshare/prog/darksusy/darksusy-5.0.5/test"  # RHEL5
#if   dist in ['RHEL5']: DIR_DARKSUSY = "/net/abel-evs/cargo/fysepf/epfshare/prog/darksusy-5.0.5/test"  # RHEL5
#elif dist in ['RHEL6']: DIR_DARKSUSY = "/net/abel-evs/cargo/fysepf/epfshare/prog/darksusy-5.1.1_RHEL6/test"     # RHEL6
#else:
#    DIR_DARKSUSY = "/net/abel-evs/cargo/fysepf/epfshare/prog/darksusy/darksusy-5.0.5/test"
#    print 'Warning  darksusy  non-recognised distribution: %s  (Proceeding as if RHEL5)' %(dist)

# NOTE: 5.0.5 compiles on RHEL5 but not on RHEL6; 5.1.1 does the opposite   # 2014-08-25
#             Hence: RHEL5 uses 5.0.5 while RHEL6 uses 5.1.1 (could be some very minor changes)

    
DARKSUSY0 = 'dsmain'

DARKSUSY_EXE = "%s/%s" %(DIR_DARKSUSY, DARKSUSY0)  
DARKSUSY_STEER = "%s/%s" %(DIR_DARKSUSY,'ATLASsteer')
#DARKSUSY_LHA = "%s/%s" %(DIR_DARKSUSY,'ATLASlha')

DIRTAG = ''

DIR_START = os.getcwd()
DIR_RUN0 = "run_ds"
DIR_RUN0FULL = '%s/%s' %(DIR_START, DIR_RUN0)  # use to check if dir name is input
fn_dsmain_out = "out_dsmain"
SAVEds     = 1
SAVEds0    = 1
SAVEdslong = 0

KEEP_DIR_RUN = 0
userATLASsteer = []

preVal = ''
preVar = ''

fnadd = ''

######################################### ARGUMENT READING
Argv = list(sys.argv)
dumA = []


if testarg(Argv,'-h') == 1:
    doHELP()
    sys.exit()

if testarg(Argv,'--automan') == 1:
    autoMan()
    sys.exit()

if testarg(Argv,'--vb') == 1: VB = 1
if testargget(Argv,'-vb',dumA) == 1: VB = int(dumA[0])

# --------------------------------------
# FOR 1 FILE ONLY (non-piped)
if testargget(Argv,'-lha',dumA) == 1:
    fn_lha = [dumA[0]]

if testargget(Argv,'-base',dumA) == 1:
    base = dumA[0]
    if os.path.exists(base+".wig"): fn_wig = [base+".wig"]
    if os.path.exists(base+"_wig"): fn_wig = [base+"_wig"]
    if os.path.exists(base+".out"): fn_out = [base+".out"]
    if os.path.exists(base+"_out"): fn_out = [base+"_out"]

    if not (fn_wig and fn_out): sys.exit("darksusy.py: FATAL  base: %s" %(base))

if testargget(Argv,'-wigout',dumA) == 1:
    z = dumA[0].split(',')
    fn_wig = [z[0]]
    fn_out = [z[1]]
    #if not os.path.exists(fn_wig): sys.exit("darksusy: FATAL  non-existent fn_wig: %s" %(fn_wig[0]))
    #if not os.path.exists(fn_out): sys.exit("darksusy: FATAL  non-existent fn_out: %s" %(fn_out[0]))
# --------------------------------------

# FOR LIST OF FILES (piped)

if testarg(Argv,'--plha') == 1:
    lines = sys.stdin.readlines()
    for L in lines:
        fn_lha.append(L.strip())

if testarg(Argv,'--pbase') == 1:
    lines = sys.stdin.readlines()
    for L in lines:
        base = L.strip()
        if os.path.exists(base+".wig"): fn_wig.append(base+".wig")
        if os.path.exists(base+"_wig"): fn_wig.append(base+"_wig")
        if os.path.exists(base+".out"): fn_out.append(base+".out")
        if os.path.exists(base+"_out"): fn_out.append(base+"_out")
        
if testarg(Argv,'--pwig') == 1:
    lines = sys.stdin.readlines()
    for iL in range(len(lines)):
        wig = lines[iL].strip()
        fn_wig.append(wig)
        base = wig.replace("_wig","").replace(".wig","")
        if os.path.exists(base+".out"): fn_out.append(base+".out")
        elif os.path.exists(base+"_out"): fn_out.append(base+"_out")
        else: sys.exit("darksusy: FATAL  No out file matching (%i) wig: %s" %(iL, wig))
                      

if testarg(Argv,'--pwigout') == 1:
    lines = sys.stdin.readlines()
    for L in lines:
        w = L.strip.split()
        fn_wig.append(w)
        fn_out.append(w)
        
if testargget(Argv,'-rundirtag',dumA) == 1:
    DIRTAG = dumA[0]
    if not DIRTAG.startswith('_'): DIRTAG = '_' + DIRTAG

if testargget(Argv,'-rundir',dumA) == 1:
    DIR_RUN0 = dumA[0]


# --------------------------------------
if testargget(Argv,'-repl',dumA) == 1:
    w = dumA[0].split(',')
    nw = len(w)/2
    if len(w) != 2*nw: sys.exit("darksusy.py: non-allowed: repl: %s" %(dumA[0]))
    for i in range(int(nw)):
        fn_tagReplace.append(w[2*i])
        fn_tagWith.append(w[2*i+1])

if testarg(Argv,'--savedslong') == 1:
    SAVEdslong = 1
if testargget(Argv,'-savedslong',dumA) == 1:
    SAVEdslong = int(dumA[0])
if testargget(Argv,'-saveds',dumA) == 1:
    SAVEds = int(dumA[0])
if testargget(Argv,'-saveds0',dumA) == 1:
    SAVEds0 = int(dumA[0])
if testarg(Argv,'--nods') == 1:
    SAVEds = 0
if testarg(Argv,'--nods0') == 1:
    SAVEds0 = 0


if testarg(Argv,'--keeprundir') == 1: 
    KEEP_DIR_RUN = 1


if testarg(Argv,'--cdmwithoutcoann') == 1:
    userATLASsteer.append('0  #dooh2with (without coannihilations')

if testarg(Argv,'--cdmwithcoann') == 1:
    userATLASsteer.append('1  #dooh2with (with coannihilations')

if testargget(Argv,'-cdmcoann',dumA) == 1:
    userATLASsteer.append('%s  #dooh2with (with coannihilations' %(dumA[0]))

if testargget(Argv,'-prevar',dumA) == 1:   # 2014-01-26
    preVar = dumA[0]
            
if testargget(Argv,'-preval',dumA) == 1: 
    preVal = dumA[0]
            
if testargget(Argv,'-fnadd',dumA) == 1: 
    fnadd = dumA[0]
    if not fnadd.startswith('_'): fnadd = '_' + fnadd
            

if len(Argv)>1:
    doHELP()
    print("Unresolved arguments: ", Argv)
    sys.exit()



######################################### CHECKS & PREPS

DIR_RUN = "%s/%s%s" %(DIR_START, DIR_RUN0, DIRTAG)


# Make lha from wig & out (put in /tmp/)
if fn_wig and fn_out and not fn_lha:
    if VB>=3: print("darksusy.py: Making lha file")
    for iS in range(len(fn_wig)): 
        fn_lha0.append(fn_wig[iS].split("/").pop().replace(".wig","").replace("_wig","") + ".lha")
        fn_lha.append("%s/%s" %(DIR_TMP, fn_lha0[iS]))
        cmd = "isa2lha.py %s -out %s > %s" %(fn_wig[iS], fn_out[iS], fn_lha[iS])
        if VB>=4: print("%3i |  %s" %(iS,cmd))
        os.system(cmd)
        if VB>=4: os.system("ls -l %s" %(fn_lha[iS]))

g.nlha = len(fn_lha)

if g.nlha == 0: sys.exit("darksusy.py  FATAL: No lha (or wig/out) files given")

# Make lha0 (lha0 is only the file name, no path)
if fn_lha and not fn_lha0:
    for iS in range(len(fn_lha)):
        fn_lha0.append(fn_lha[iS].split("/").pop())

# Make tag 
if fn_lha:
    if VB>=3: print("darksusy.py: making tag from lha (n=%i)" %(g.nlha))
    for iS in range(len(fn_lha)):
        z = fn_lha0[iS].replace(".lha","")
        for i in range(len(fn_tagReplace)):
            z = z.replace(fn_tagReplace[i],fn_tagWith[i])
        if VB>=3: print("%3i  tag: %s" %(iS, z))
        fn_tag.append(z)


######################################### EXECUTION

#lines = sys.stdin.readlines()

Main()

