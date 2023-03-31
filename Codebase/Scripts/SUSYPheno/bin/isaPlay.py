#!/usr/bin/env python

#
# This is maybe the last version before renaming it as 
# isaPlay.py
#
##############################################################################
# DESCRIPTION
#   This script produces easily-readable 'derived' files of masses and BRs
#   based on isawig output
# 
# USAGE
#       isaPlay.py  <file_isawig>  <BR_cut>  [info_text]
#   Ex: isaPlay.py  su1.new.data  0.001       > su1.new.derived
#   Ex: isaPlay.py  su1.new.data  0.001  SU1  > su1.new.derived
#
# DEVELOPMENT
#      2007 Version 1: Borge Kile Gjelsten <gjelsten@cern.ch>
#      2008 Version 2
#      2010 Version 3: - Name change from StructureISAWIGfiles.py to isaPlay.py
#                      - Some global action put into function
#                      - Add 1leg and 2leg methods (no editing/usage between May3 and Sep14)
#  Sep 2010 Version 4: - Added study3(), a clean structure for quickly assessing the event types for a given scenario
#                      - Fixed some bugs
#                      - Changed sparticle naming from old to new (libISAWIG too)
#
#############################################################################
# Interesting features to add: 
#
#
# Refreshing:
#  - May: used as test work directory: /scratch2/borgeg/cargo/alt/prog/isajet/isa778/scan/cascades_test/dir_casc/
#
#
#############################################################################
# COMMENTS TO STUDY3L (SEP 2010)
# ------------------------------
# - only tested on SU4, expect things to happen in other scenarios: 
#    - wino-LSP will give trouble in that new SM particles enter (pions): will need to extend structure
# - May still be some undesirable "features" due to incorrect copying of list/dict
# - Results not really tested: should find a way to test on a given scenario
# -
#
#############################################################################



# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>  PREPARATIONS
import sys,os,string,re,math


path1 = "/mn/kvant/u1/borgeg/bin/"
file1 = "mytestarg_lib.py"
if os.path.exists(file1): exec(compile(open(file1, "rb").read(), file1, 'exec'))
elif os.path.exists(path1+file1): exec(compile(open(path1+file1, "rb").read(), path1+file1, 'exec'))
else:
    print("FATAL: file not available: %s" %file1)
    sys.exit()



from libISAWIG import *
# libISAWIG import gives output like 'registering stepper ClassicalRK4', ... Apr25
# apparently finds libISAWIG in ~/bin/ 


# Sep16
# MyTimeTools.py needs to lie in a location contained in the PYTHONPATH variable,
# e.g. export PYTHONPATH:$PYTHONPATH:.  and put MyTimeTools.py in the work directory.
# 
from MyTimeTools import *


path1 = "/mn/kvant/u1/borgeg/"
file1 = ".pythonrc2"
if os.path.exists(file1): exec(compile(open(file1, "rb").read(), file1, 'exec'))
elif os.path.exists(path1+file1): exec(compile(open(path1+file1, "rb").read(), path1+file1, 'exec'))
else:
    "WARNING: file not available: %s" %file1
    "May have problem enjoying history and tab-completion functionality in iPython"



# =============================================================================



if(len(sys.argv)<2 or (sys.argv[1]=="-h")):
    print("  NB: Not all help messages are up2date, but STUDY LEPTON examples are")
    print("    TO PRINT DECAYS:")
    print("          isaPlay.py  <wigfilename>  [BRmin]  [COMMENT_ON_TOP]")
    print("      Ex: isaPlay.py  dc1new.txt")
    print("      Ex: isaPlay.py  dc1new.txt  0.001")
    print("      Ex: isaPlay.py  dc1new.txt  0.001  [SU3]")
    print("      Ex: isaPlay.py  dc1new.txt  0.01  short   #show only 'most relevant' sparticles")
    print("      Ex: isaPlay.py  dc1new.txt  0.01  masses  #show only masses")
    print("      Ex: isaPlay.py  dc1new.txt  0.01  lsep    #show 1./2. gen leptons separately")
    print("      Ex: isaPlay.py  dc1new.txt  0.01  qsep    #show 1./2. gen quarks separately")
    print("      Ex: isaPlay.py  dc1new.txt  0.01  -maxmass 1000  #shows only sparticles of mass<1000 GeV")
    print() 
    print("    TO ANALYSE 2LEGS:")
    print("      Ex: isaPlay.py  ../dir_isawig/sps1a_isawig  -2leg  -brmin2leg 1e-4  -prox ../dir_ox/sps1a_7TeV.ox")
    print("      Ex: isaPlay.py  ../dir_isawig/sps1a_isawig  --qsep -2leg  -brmin2leg 1e-4  -prox ../dir_ox/sps1a_7TeV.ox")
    print("      Ex: python -i ~/bin/isaPlay.py ../dir_isawig/sps1a_isawig  -2leg  -prox ../dir_ox/sps1a_7TeV_ox  -brmin2leg 1e-4  -brmin1leg 1e-4")
    print()
    print("    TO STUDY LEPTONS")
    print("      Ex: python -i ~/bin/isaPlay.py  ../dir_isawigsu/msugra_su4_wig  -prox ../dir_oxsu/su4_7TeV.ox  -study3l -brmin2leg 1e-3")
    print("      Ex: python -i ~/bin/isaPlay.py  ../dir_isawigsu/msugra_su4_wig  -prox ../dir_oxsu/su4_hacked_7TeV.ox  -study3l -brmin2leg 1e-5 -SMout l,q ")
    print("      Ex: python -i ~/bin/isaPlay.py  ../dir_isawigsu/msugra_su4_wig  -prox ../dir_oxsu/su4_hacked_7TeV.ox  -study3l -brmin2leg 1e-5 -SMout l -showSMp -VBscreen 0 ")
    print("      Ex: isaPlay.py  ../dir_isawigsu/msugra_su4_wig  -prox ../dir_oxsu/su4_hacked_7TeV.ox  -study3l -brmin2leg 1e-5 -SMout l -showSMp -VBscreen 0 ")
    print("      Ex: ./isaPlay.py  msugra_su4_wig  -prox su4_hacked_7TeV.ox  -study3l -brmin2leg 1e-5 -SMout l -showSMp")
    print("      Ex: ./isaPlay.py  msugra_su4_wig  -prox su4_hacked_7TeV.ox  -study3l -brmin2leg 1e-5 -SMout l -showSMp -VBscreen 0")
    print()
    print("      Ex:    cd ~/grids_lsp/DGnoSL_TB10_2013-09-16/   #  2013-11-08")
    print("             for fnbase in `ls lhawig/ | grep wig | sed s/.wig// ` ; do isaPlay.py lhawig/${fnbase}.wig  -prox PROX/${fnbase}.ox  -scen ${fnbase}  -study3l ; done ")
    print()
    print("      Ex:    isaPlay.py  <wigfile>  -maxmass 900  -brmin 0.0001   # print nice&old der file    2013-12-21")
    print("      Ex:    isaPlay.py  <wigfile>  -maxmass 900  --m             # print nice&old mass file   2013-12-21")
    print() 
    sys.exit()





# ============================================================================= GLOBAL VARS (#1)
# ============================================================================= 
# ============================================================================= STEERING

# Global variables needed (extremely inelegant)
N1b = -1.
N1w = -1.
N1h1 = -1.
N1h2 = -1.
N1h = -1.
N2b = -1.
N2w = -1.
N2h1 = -1.
N2h2 = -1.
N2h = -1.
N3b = -1.
N3w = -1.
N3h1 = -1.
N3h2 = -1.
N3h = -1.
N4b = -1.
N4w = -1.
N4h1 = -1.
N4h2 = -1.
N4h = -1.


singammaR = -1.
cosgammaR = -1.
singammaL = -1.
cosgammaL = -1.
C1wL = -1.
C1hL = -1.
C1wR = -1.
C1hR = -1.
C2wL = -1.
C2hL = -1.
C2wR = -1.
C2hR = -1.

thetat = -1.
thetab = -1.
thetatau = -1.
t1L = -1.
t1R = -1.
t2R = -1.
t2L = -1.
b1L = -1.
b1R = -1.
b2R = -1.
b2L = -1.
tau1L = -1.
tau1R = -1.
tau2R = -1.
tau2L = -1.

At = -1.
Ab = -1.
Atau = -1.

col = -1.
mu = -1.



# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# <><><><><><><><><><><><><>  STEERING  <><><><><><><><><><><><><><><><><><><>

PRINTDECAY  = 1   #D=1  main steer
PRINTMASSES = 0  #D=0  main steer
MODEMASSES = 0

PRINTCOMPOSITION = 2  #1:join hU and hD,  2:separate hU and hD
PRINTMASSESOFDAUGHTERS=3
PRINTFULLPARENTEACHLINE=1
PRINTPARMASS=2 #1:old; 2:new(mass left)
PRINTSUM=1
BRmin=0.001
BRPOS=1  #1:left; 2:right
NOMASSESABOVE10TeV=1  #has no effect on other stuff
PRINTHEAD=0
PRINTUSEINFO=0  #D=1
PRINTLINE=1
NSPACE=1
NDECMAX=99
PRINTARROW=1
massEq="m="
PRINTWIDTH = 1  #0,1
#PRINTCOMP = 2   

MAXMASSPRINT = 99999.

PRINTTOFILE = 0
PRINTTOFILEONLY = 0

PSKIP = []


SHORT=0   #D=0  (isaPlayShort has this 1, is the only difference)

LeptonMode=0
QuarkMode=0

#donotprint=["t","h","H","A","H+"]
donotprint=["t"]
onlyprint=[]

shownotpart = []

COMMENT_ON_TOP=""

#CASCADES = 0  #replaced by oneleganalysis and twoleganalysis (so far)
ONELEGANALYSIS = 0
TWOLEGANALYSIS = 0
STUDY3L = 0

# used with 1legs
BRmin1leg = 1e-5

# used with 2legs
BRmin2leg = 1e-3

scenID = ""
SHOWCRUCIALS = True

PrOX_filename = ""  # PrOX = ProspinoOsloXsec

TWOLEG_HIST = 0

SMout = ['l','q','b','T','y']  #default as of Sep27

VBscreen = 1
VBfile = 0

saveSMp = 0
showSMp = 0

outdir = 'RES_ISAPLAY'
#fnSMp = 'outSMp.txt'
fnSMp0 = 'outSMp'

scenID = ''  # outSMp files will get scenID prepended

outlabel = ''

SMout_reductionPath = ['y','T','b','q']  #this is what currently triggers output files, actually  # leaves only l left at the end
#SMout_reductionPath = ['y','b','q','l']  #this is what currently triggers output files, actually # leaves only T left at the end
# NB: may also need to do changes to libISAWIG:  pOlds = 
SMout_reductionPath0 = list(SMout_reductionPath)  # constant, just to keep

checkmll = 0

LSPid = "N1"


# ============================================================================= 
# ============================================================================= ARGUMENT TESTING
# ============================================================================= 

mode_casc_starttypes = 'all'

Argv=sys.argv
dumA=[]
VERBARG=1


if testarg(Argv,"--automan") == 1:
    autoMan()
    sys.exit()


# this is a very old setup ... these few lines below should probably be put closer to the bottom
filen=Argv[1]
filen2=filen.split("/"); filen2=filen2[len(filen2)-1]
del(Argv[1])

# =============================================================
#print Argv, testarg(Argv,"--automan"), Argv

if testargget(Argv,"-reductionpath",dumA) == 1:  # 2012-06-23
   SMout_reductionPath = dumA[0].split(',')
   print("Old reduction path: %s   -->  New reduction path: %s" %(str(SMout_reductionPath0), str(SMout_reductionPath)))
   
if(testargget(Argv,"-brmin",dumA)==1):
    BRmin=float(dumA[0])
if(testargget(Argv,"-head",dumA)==1):
    COMMENT_ON_TOP=dumA[0]
if(testargget(Argv,"-nspace",dumA)==1):
    NSPACE=int(dumA[0])
if(testargget(Argv,"-ndecmax",dumA)==1):
    NDECMAX=int(dumA[0])
if(testarg(Argv,"--noline")==1):
    PRINTLINE=0
if(testarg(Argv,"--fullpar")==1):
    PRINTFULLPARENTEACHLINE=1
if(testarg(Argv,"--nofullpar")==1):
    PRINTFULLPARENTEACHLINE=0
if(testarg(Argv,"--massdau")==1):
    PRINTMASSESOFDAUGHTERS=1
if(testarg(Argv,"-massdau","1")==1):
    PRINTMASSESOFDAUGHTERS=1
if(testarg(Argv,"-massdau","2")==1):
    PRINTMASSESOFDAUGHTERS=2
if(testarg(Argv,"-massdau","3")==1):
    PRINTMASSESOFDAUGHTERS=3
if(testarg(Argv,"--nomeq")==1):
    massEq=""
if(testarg(Argv,"--arrow")==1 or testarg(Argv,"--to")==1):
    PRINTARROW=1
if(testarg(Argv,"--noarrow")==1 or testarg(Argv,"--noto")==1):
    PRINTARROW=0
if(testargget(Argv,"-brpos",dumA)==1):
    BRPOS=int(dumA[0])

if(testargget(Argv,"-f",dumA)==1):
    PRINTTOFILE = 1
    filenprint = filen.split("/").pop().replace("isawig","derived")+dumA[0]
if(testargget(Argv,"-f0",dumA)==1):
    PRINTTOFILE = 1
    PRINTTOFILEONLY = 1
    filenprint = filen.split("/").pop().replace("isawig","derived")+dumA[0]


if testargget(Argv,"-pskip",dumA)==1:    #skip particles in DECAY-printout
    # ex: '-pskip H+,H,A,h,N4,N3,C2,grav'
    zdum = dumA[0].split(",")
    for zdum2 in zdum: PSKIP.append(zdum2)  

    
if testarg(Argv,"-m1line")==1:
    PRINTDECAY=0
    MODEMASSES = 2
    PRINTMASSES=1

if testarg(Argv,"-m1linem")==1:
    PRINTDECAY=0
    MODEMASSES = 1
    PRINTMASSES=1
    
if(testarg(Argv,"--masses")==1 or testarg(Argv,"--m")==1):
    #SHORT=0
    PRINTDECAY=0
    PRINTMASSES=1
    MODEMASSES = 0 
    #LeptonMode=1  #D=0 (comment out)
    #QuarkMode=1   #D=0
if(testarg(Argv,"--noshort")==1):
    SHORT=0
    #PRINTDECAY=1
    #PRINTMASSES=0
    #PRINTUSEINFO=0
    #PRINTSUM=0
if(testarg(Argv,"--short")==1):
    SHORT=1
    #PRINTDECAY=1
    #PRINTMASSES=0
    #PRINTUSEINFO=0
    #PRINTSUM=0
if(testarg(Argv,"--lsep")==1 or testarg(Argv,"--sleptons")==1): 
    #SHORT=0
    #PRINTDECAY=1
    #PRINTMASSES=0
    #PRINTSUM=1
    LeptonMode=1
if(testarg(Argv,"--qsep")==1 or testarg(Argv,"--squarks")==1): 
    #SHORT=0
    #PRINTDECAY=1
    #PRINTMASSES=0
    #PRINTSUM=1
    QuarkMode=1
if(testargget(Argv,"-part",dumA)==1):
    #SHORT=0
    #PRINTUSEINFO=0
    #PRINTDECAY=1
    #PRINTMASSES=0
    #PRINTSUM=0
    QuarkMode=1
    LeptonMode=1
    #print dumA[0]
    if(dumA[0]=="1"):
        onlyprint += ["N1","C1","N2","eR","muR","ta1"];
    elif(dumA[0]=="2"):
        onlyprint += ["N1","C1","N2","eR","muR","ta1","gl","uL","uR"];
    elif(dumA[0]=="3"):
        onlyprint += ["h","N1","C1","N2","N3","C2","N4","eR","eL","ve","ta1","gl","uL","uR","dL","dR","t1","t2","b1","b2"];
    else:
        #onlyprint += dumA[0].split()  # May3: replaced this 
        onlyprint += dumA[0].split(',')

if(testargget(Argv,"-notpart",dumA)==1):
    shownotpart = dumA[0].split(',')

if(testarg(Argv,"-width")==1): PRINTWIDTH = 1
if(testarg(Argv,"-nowidth")==1): PRINTWIDTH = 0



if(testarg(Argv,"--nohead")==1):
    PRINTHEAD=0
if(testarg(Argv,"--head")==1):
    PRINTHEAD=1
if(testarg(Argv,"--noinfo")==1):
    PRINTUSEINFO=0
if(testarg(Argv,"--info")==1):
    PRINTUSEINFO=1
if(testargget(Argv,"-info",dumA)==1):
    if(dumA[0]=="0"):
        PRINTUSEINFO=0
    elif(dumA[0]=="1"):
        PRINTUSEINFO=1
    else:
        print("Stopping:  '-info "+dum[0]+"'   is not understood.")
        sys.exit()
if(testarg(Argv,"--sum")==1):
    PRINTSUM=1
if(testarg(Argv,"--nosum")==1):
    PRINTSUM=0

if(testargget(Argv,"-maxmass",dumA)==1):
    MAXMASSPRINT = float(dumA[0])

if(testarg(Argv,"-casc")==1):
    #CASCADES = 1
    ONELEGANALYSIS = 1
    PRINTDECAY = 0

if(testarg(Argv,"-1leg")==1):
    #CASCADES = 1
    ONELEGANALYSIS = 1
    PRINTDECAY = 0

if(testarg(Argv,"-2leg")==1):
    #CASCADES = 1
    TWOLEGANALYSIS = 1
    PRINTDECAY = 0

if(testargget(Argv,'-brmin1leg',dumA)==1):   #for use with -casc
    BRmin1leg = float(dumA[0])

if(testargget(Argv,'-brmin2leg',dumA)==1):   #for use with -2leg
    BRmin2leg = float(dumA[0])

if(testargget(Argv,'-scen',dumA)==1):
    scenID = dumA[0]

if(testargget(Argv,'-prox',dumA)==1):
    PrOX_filename = dumA[0]
    if not os.path.exists(PrOX_filename):
        print("FATAL: Non-existent ox file: "+PrOX_filename)
        sys.exit()

if(testarg(Argv,'-2leghist')==1):
    TWOLEG_HIST = 1   #1:use hist, 0:not use [if 1, need to have pyroot]

if(testarg(Argv,'-2leghistoff')==1):
    TWOLEG_HIST = 0   #1:use hist, 0:not use [if 1, need to have pyroot]

if(testargget(Argv,'-smout',dumA)==1):
    SMout = dumA[0].split(',')
    # test
    for sm in SMout:
        if sm not in ['l','q','b','T','y']:
            print("FATAL: Non-allowed SMout particle: %s" %(sm))
            sys.exit()

if(testargget(Argv,'-vbscreen',dumA)==1):
    VBscreen = int(dumA[0])
    
if(testargget(Argv,'-vbfile',dumA)==1):
    VBfile = int(dumA[0])

if(testarg(Argv,'-showsmp')==1):
    showSMp = 1
if(testarg(Argv,'-noshowsmp')==1):
    showSMp = 0

if(testarg(Argv,'-savesmp')==1):
    saveSMp = 1
if(testarg(Argv,'-nosavesmp')==1):
    saveSMp = 0

if(testargget(Argv,'-fnsmp',dumA)==1):
    print("NB: might be unwise to change fnSMp from default %s to %s" %(fnSMp0,dumA[0]))
    fnSMp0 = dumA[0]

if(testargget(Argv,'-scen',dumA)==1):
    scenID = dumA[0]

    
if(testarg(Argv,"-study3l")==1):
    #CASCADES = 1
    STUDY3L = 1
    PRINTDECAY = 0
    LeptonMode = 1
    QuarkMode = 1
    if VBscreen: print("NB: study3L sets LeptonMode = 1 and QuarkMode = 1, equivalent to --qsep --lsep")

if testargget(Argv,"-outdir",dumA)==1:
    print("NB: changing from default outdir %s to %s" %(outdir, dumA[0]))
    outdir = dumA[0]

if testargget(Argv,"-starttypes",dumA)==1:
    mode_casc_starttypes = dumA[0]

if testarg(Argv,"-nodecay") == 1:
    PRINTDECAY = 0

# ============================================================= 2010-12-12
if(testargget(Argv,"-quarkmode",dumA)==1): 
    QuarkMode = int(dumA[0])

if(testargget(Argv,"-leptonmode",dumA)==1): 
    LeptonMode = int(dumA[0])
# =============================================================

if testargget(Argv,"-lsp",dumA)==1:
    LSPid = dumA[0]

if(testarg(Argv,"-checkmll")==1): 
    checkmll = 1
# =============================================================
if(len(Argv)>1):
    print("Argv: ", end=' ')
    print(Argv)
    print("Unresolved arguments in command line. Exiting.")
    sys.exit()

if(not(os.path.exists(filen))):
    print("File "+filen+" does not exist..")
    sys.exit()


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>



ifile=open(filen,"r"); line=ifile.readlines(); ifile.close()

if ONELEGANALYSIS or TWOLEGANALYSIS or STUDY3L:  # Mar4
    if not os.path.exists(outdir): os.mkdir(outdir)

fnSMp = outdir + "/" + fnSMp0
if scenID: fnSMp += "_" + scenID

# ============================================================================= 
# ============================================================================= GLOBAL VARS
# ============================================================================= 

anDec=[]
aDecDau=[] #will get additional [anDec][3]
aDecBR=[]  #will get additional [anDec]
#aDecWidth=[]
aDecnJoin=[]
aDecnDau=[] #how many daughters (2 or 3?) in this decay chain

aIDlog=[]
for i in range(500):
    aIDlog.append(-1)
aIDgen=[-1]
m=["-1"]
w=[-1]
t=[-1]
ct=[-1]
nam=["nam"]
rm=[]  #real mass ... i.e. not 'text mass'

nIDlog = -9

ziterMO = []
iterMO = []
namMO = []


# ============================================================================= BEGIN
# ============================================================================= FUNCTIONS
# =============================================================================


# =============================================================================
# analysis2legs moved to top of file (here) not to confuse with study3l in editing process
def analysis2legs(mode=5):

    if mode < 1: return
    if g.VB >=2: print("analysis2legs(1): make1legs   " + 30*"-")
    make1legs()

    if mode < 2: return
    #if g.VB >=2: print "analysis2legs(2): getPrOX " + 30*"-"
    g.ox = getPrOXsec(g.PrOX_filename)
    g.OX = OX(g.ox[0], g.ox[1], g.ox[2])  # OX-object (useful)

    if mode < 3: return
    if g.VB >=2: print("analysis2legs(3): make2legs " + 30*"-")
    make2legs()
    #g.twolegs.Show(5)

    # ----------------------------------------------------------- BELOW: refining
    if mode < 4: return
    #if g.VB >=2: print "analysis2legs(4): refine2legs "+30*"-"
    
    #refine2legs()
    DB = True

    g.twolegs.Show(0,txt='<--- Just after making all 1leg x 1leg combinations')
    if DB: g.twolegs.Show(fn='an2legs0_brute2leg.txt')


    # NBNBNB check is very relevant input to Duplicates('join')

    #check = ['SU']  
    check = ['SU','SM'] 

    if check == ['SU']: NOTIFY.append("check = %s : only the two SUSY chains are checked. SM information is rendered incorrect when 2legs are joined" %check)
    if check == ['SU','SM']: NOTIFY.append("check = %s : SUSY chains and SM chains are checked. SM information remains correct when 2legs are joined" %check)
    
    # ['LRboth', 'SU', 'SM'/'SMa']  any combination of these three (not SM and SMa together)
    # LRcheckboth = False #True
    # True:  2legs will keep their original LR-order (not BR-order, of course)
    # False: 2legs will get LR-order according to 'value' of SU1 vs SU2


    # interchange the two legs per 2leg if needed to get SU1 < SU2
    if 'LRboth' not in check: 
        g.twolegs.ShortSmallestLeft()

    if DB: g.twolegs.Show(fn='an2legs0a_brute2leg.txt')  #still leptons (_ll) here 

    # join duplicates
    g.twolegs.Duplicates('join',check) #LRcheckboth)
    if DB: g.twolegs.Show(fn='an2legs0b_brute2leg.txt')  #leptons gone...
    g.twolegs.orderWithBR()
    if DB: g.twolegs.Show(fn='an2legs1_ordered.txt')
    if DB: g.twolegs.Show(0,txt='<--- After joining identical chains')

    # CONFORM INITIAL PARTICLES (were left out earlier)
    #g.twolegs.JoinQs(pJoin="qL", start=0, end=1) #now included in 'SQ'
    #g.twolegs.JoinQs(pJoin="qR", start=0, end=1) #now included in 'SQ'
    g.twolegs.JoinQs(pJoin="SQ", start=0, end=1)

    
    g.twolegs.Duplicates('join',check)
    g.twolegs.orderWithBR()
    if DB: g.twolegs.Show(fn='an2legs2_conforminitial.txt')
    if DB: g.twolegs.Show(0,txt='<--- After uL,uR,dL,... -> SQ and join identicals')
    
    # Ni -> N_ 
    if 1: #new
        g.twolegs.JoinNs()
        g.twolegs.Duplicates('join',check)
        g.twolegs.orderWithBR()
        if DB: g.twolegs.Show(fn='an2legs_joinNs.txt')
        if DB: g.twolegs.Show(0,txt='Ni->N')

    # Ci -> C_
    if 1: #new 
        g.twolegs.JoinCs()
        g.twolegs.Duplicates('join',check)
        g.twolegs.orderWithBR()
        if DB: g.twolegs.Show(fn='an2legs_joinCs.txt')
        g.twolegs.Show(0,txt='Ci->C')

    g.twolegs.ShortSmallestLeft() ; g.twolegs.Duplicates('join',check)  #ShortSmallestLeft (has effect here)
    g.twolegs.Show(fn=g.it.next('pres','_')+g.scenID+'_an2legs_gl_SQ_N_C_SL_orderBR.txt')  #<--

    # N_,C_ -> X_
    g.twolegs.JoinNCs()
    g.twolegs.Duplicates('join',check)
    g.twolegs.orderWithBR()
    if DB: g.twolegs.Show(fn='an2legs3_joinNCs.txt')
    if DB: g.twolegs.Show(0,txt='<--- After N,C -> X')
    
    g.twolegs.ShortSmallestLeft() ; g.twolegs.Duplicates('join',check)  #ShortSmallestLeft (has effect here)
    g.twolegs.orderWithBR()
    g.twolegs.Show(fn=g.it.next('pres','_')+g.scenID+'_an2legs_gl_SQ_X_SL_orderBR.txt')  #<--

    # ShortSmallestLeft
    g.twolegs.ShortSmallestLeft() ; g.twolegs.Duplicates('join',check)  #ShortSmallestLeft
    g.twolegs.orderWithBR()
    if DB: g.twolegs.Show(fn='an2legs4_ShortSmallestLeft.txt')
    #if DB: g.twolegs.Show(0,txt='ShortSmallestLeft')
    
    # ----------------------------------------------------------- ABOVE: refining


    g.twolegs.Show(0,txt='<--- 2legs: final')
    #g.twolegs.Show(fn=g.scenID+'_an2legs_final_orderBR.txt')
    g.twolegs.orderWithLegvalue()
    g.twolegs.Show(fn=g.it.next('pres','_')+g.scenID+'_an2legs_final_orderLegvalue.txt')

    # Here refine SMp content .. (not yet implmented)


    # Then act on SMp
    if 1: 
        g.twolegs.FillSMp()
        g.twolegs.Duplicates('join',['SU','SMp']) # NB: now losing SM1 and SM2 relevance (info incorrect after this)
        g.twolegs.Show(5,txt="SU SMp")
        g.twolegs.orderWithLegvalue()
        g.twolegs.Show(fn=g.it.next('pres','_')+g.scenID+'_an2legs_final_orderLegvalueSMp.txt',SM=['+SMp'])  #nohead?

        #g.twolegs.FillSMp()
        g.twolegs.Duplicates('join',['SMp'])  # NB: now losing SU1 and SU2 (info incorrect after this)
        g.twolegs.orderWithSMp()
        #g.twolegs.Show(fn=g.scenID+'_an2legs_final_orderSMp.txt',SM=['+SMp'])


        # Simplify SMq
        g.twolegs.Show(0,txt='Before Expand')

        # Neat construction. ExpandSMp returns True if was expanded. If True, then enters block ('pass') and reruns test. 
        while g.twolegs.ExpandSMp(mode=['tt','t','h','H','A','H+','W','Z','b','T','QQ','vv','TT','Tv','lv','ll']):  # continues until False is returned (i.e. all expansion done)
            #print 'Expanded'
            pass
        g.twolegs.Show(0,txt='Right after Expand')
        g.twolegs.Duplicates('join',check=['SMp'])  # NB: lose whatever SU info is left
        g.twolegs.orderWithSMp()
        g.twolegs.Show(0,txt='Right after Expand')
        #g.twolegs.Show(fn=g.scenID+'_an2legs_final_orderSMpExp.txt',SM=['+SMp'])
        g.twolegs.Show(fn=g.it.next('pres','_')+g.scenID+'_an2legs_final_orderSMpExp2.txt',SM=['SMp'])  #head..

        # Make plot of SMp
        if TWOLEG_HIST: #here assuming ROOT to be available (e.g. athena sourced)
            g.twolegs.FillH()
            fn0 = g.it.next('pres','_')+g.scenID+'_an2legs_final_orderSMpExp2'
            cc,tleg = g.twolegs.ShowH(which=[-1,0,1,2,3,4], fn=[fn0+'.pdf', fn0+'.eps'], silent=True)
            

    if TWOLEG_HIST: return cc,tleg  #to preserve

    #if mode < 5: return
    #if g.VB >=2: print "\nanalysis2legs(5): show2legs " + 30*"-"
    #show2legs()
# =============================================================================
# =============================================================================
# =============================================================================


# =============================================================================
def ReadMassesBRsMixingparsEtc(): 

    iIDlog=1 #leave 0th slot for gluino
    iL=1 #start at line 1

    ##############################
    # A: READ MASSES
    ##############################
    
    nMlines=int(line[0])
    while (iL<nMlines+1):
        col=str.split(line[iL].strip())
        IDspec=int(col[0])
        IDgen=fIDgen(IDspec,QuarkMode,LeptonMode)
        iL+=1
        dum= str(iL)+"  "+str(IDspec)+"  "+fname(IDspec,QuarkMode,LeptonMode)  #com
        if(IDgen!=IDspec): continue
        
        # special for gluino (with this trick gluino occupies 0th place)
        if(IDgen==449):
            aIDlog[IDgen]=0
            aIDgen[0]=IDgen
            nam[0]=fname(IDgen,QuarkMode,LeptonMode)
        
            m[0]=col[1]
            t[0]=float(col[2])   # lifetime
            w[0]=6.582e-25/t[0]  # width
            ct[0]=3.0e8*t[0]     # distance traversed in detector. Should be <1mm or so to not be displaced or stable   
        else:
            aIDlog[IDgen]=iIDlog
            aIDgen.append(IDgen)
            nam.append(fname(IDgen,QuarkMode,LeptonMode))

            # print "%2i  %-4s  IDgen:%2i" %(len(nam)-1, nam[len(nam)-1], IDgen) #helper
        
            m.append(str(abs(float(col[1]))))
            t.append(float(col[2]))
            w.append(6.582e-25/t[len(t)-1])
            ct.append(3.0e8*t[len(t)-1])
            iIDlog+=1

        # print dum

    #print 'her daa',nMlines,iIDlog
    global nIDlog
    nIDlog=iIDlog

    for ii in range(len(m)): rm.append(float(m[ii]))

    if(0):
        for i in range(nIDlog):
            print(str(i).rjust(3)+"  "+nam[i].ljust(5)+"  "+str(aIDgen[i]).rjust(3)+"  "+m[i].ljust(8)+"  "+w[i])


    ################################
    # B: POPULATE BRANCHING RATIOS 
    ################################

    for i in range(nIDlog):
        anDec.append(0)
        aDecBR.append([])
        # aDecWidth.append([])
        aDecDau.append([])
        aDecnJoin.append([]) #used for testing
        aDecnDau.append([])

    dumN=0
    # sys.exit()
    while(1):
        col=str.split(line[iL].strip())
        thisline=line[iL]
        iL+=1

        dumN+=1
        # if(dumN>30): break

        if(len(col)==1): continue

        #if(int(col[0])<100): break  #jumps out of loop when done with sparticles
        # if(int(col[0])==12): break #jumps out of loop when done top+ (then remains top-)  # FRAGILE preMay23
        if len(col) == 2: break  #May23 (replacing the line above)


        IDspec=int(col[0])
        IDgen=fIDgen(IDspec,QuarkMode,LeptonMode)

        if(IDgen!=IDspec): continue  #skip if is a duplicate particle
        
        IDlog=aIDlog[IDgen]
        
        zDau=[fIDgen(int(col[3]),QuarkMode,LeptonMode),fIDgen(int(col[4]),QuarkMode,LeptonMode),fIDgen(int(col[5]),QuarkMode,LeptonMode),fIDgen(int(col[6]),QuarkMode,LeptonMode)]
        zDau.sort(); zDau.reverse()
        zBR=float(col[1])
        # check if decay chain is already there
        ziDec=-1
        for iDec in range(anDec[IDlog]):
            # print "Testing: "+str(zDau[0]).rjust(4)+str(zDau[0]).rjust(4)+str(zDau[0]).rjust(4)+"  VS "+str(aDecDau[IDlog][iDec][0]).rjust(4)+str(aDecDau[IDlog][iDec][1]).rjust(4)+str(aDecDau[IDlog][iDec][2]).rjust(4)
            if(zDau==aDecDau[IDlog][iDec]):
                ziDec=iDec
                break

        # If need to add new chain 
        if(ziDec==-1):
            aDecBR[IDlog].append(0)
            # aDecWidth[IDlog].append(0)
            aDecDau[IDlog].append(zDau)
            # aDecnDau[IDlog].append((zDau[0]>0) + (zDau[1]>0) + (zDau[2]>0))
            aDecnDau[IDlog].append((zDau[0]>0) + (zDau[1]>0) + (zDau[2]>0) + (zDau[3]>0))
            anDec[IDlog]+=1
            ziDec=len(aDecBR[IDlog])-1
            aDecnJoin[IDlog].append(0)

        # print "ziDec="+str(ziDec)
        # Then add current BR (also if new chain)
        aDecBR[IDlog][ziDec]+=zBR
        # aDecWidth[IDlog][ziDec]+=zWidth
        aDecnJoin[IDlog][ziDec]+=1
    
        # print thisline.strip()



    ################################
    # C: READ MIXING PARAMETERS ETC. 
    ################################

    if 0: #preMay23
        iL+=4 #first skip tbar lines # ... if there are new particles below (e.g. light H+), this trick will fail
        col=str.split(line[iL].strip()); iL+=1


    tanb=float(col[0]); alphaH=float(col[1])
    #print 'tanb, alphaH: ', tanb, alphaH

    # NB: notation below confirmed Sep28,2009; stems exactly from comment isajet.f line 42620
    #     But what is higgs1 and higgs2 ?
    #     See from hep-ph/0311123 (reference in isajet) that h1=Hd and h2=Hu
    #     same ref: alphaH is zero, then (h0,H0) = (Hu0,Hd0)   [so lightest higgs is usually Hu0]

    #print line[iL].strip()
    col=str.split(line[iL].strip()); iL+=1
    #print col


    global N1b,N1w,N1h1,N1h2,N1h,N2b,N2w,N2h1,N2h2,N2h,N3b,N3w,N3h1,N3h2,N3h,N4b,N4w,N4h1,N4h2,N4h
    global singammaR,cosgammaR,singammaL,cosgammaL,C1wL,C1hL,C1wR,C1hR,C2wL,C2hL,C2wR,C2hR
    global thetat,thetab,thetatau,t1L,t1R,t2R,t2L,b1L,b1R,b2R,b2L,tau1L,tau1R,tau2R,tau2L,At,Ab,Atau,mu


    N1b=float(col[0])**2; N1w=float(col[1])**2; N1h1=float(col[2])**2; N1h2=float(col[3])**2; N1h=N1h1+N1h2
    col=str.split(line[iL].strip()); iL+=1
    N2b=float(col[0])**2; N2w=float(col[1])**2; N2h1=float(col[2])**2; N2h2=float(col[3])**2; N2h=N2h1+N2h2
    col=str.split(line[iL].strip()); iL+=1
    N3b=float(col[0])**2; N3w=float(col[1])**2; N3h1=float(col[2])**2; N3h2=float(col[3])**2; N3h=N3h1+N3h2
    col=str.split(line[iL].strip()); iL+=1
    N4b=float(col[0])**2; N4w=float(col[1])**2; N4h1=float(col[2])**2; N4h2=float(col[3])**2; N4h=N4h1+N4h2

    col=str.split(line[iL].strip()); iL+=1; singammaR=float(col[0]); cosgammaR=float(col[1])
    col=str.split(line[iL].strip()); iL+=1; singammaL=float(col[0]); cosgammaL=float(col[1])
    C1wL=singammaL**2; C1hL=cosgammaL**2; C1wR=singammaR**2; C1hR=cosgammaR**2; 
    C2wL=cosgammaL**2; C2hL=singammaL**2; C2wR=cosgammaR**2; C2hR=singammaR**2; 

    col=str.split(line[iL].strip()); iL+=1
    thetat=float(col[0]); thetab=float(col[1]); thetatau=float(col[2])
    t1L=math.cos(thetat)**2; t1R=math.sin(thetat)**2
    t2R=math.cos(thetat)**2; t2L=math.sin(thetat)**2
    b1L=math.cos(thetab)**2; b1R=math.sin(thetab)**2
    b2R=math.cos(thetab)**2; b2L=math.sin(thetab)**2
    tau1L=math.cos(thetatau)**2; tau1R=math.sin(thetatau)**2
    tau2R=math.cos(thetatau)**2; tau2L=math.sin(thetatau)**2

    col=str.split(line[iL].strip()); iL+=1
    At=float(col[0]); Ab=float(col[1]); Atau=float(col[2])

    col=str.split(line[iL].strip()); iL+=1
    mu=float(col[0])
    # -------------------------------------------------
    # -------------------------------------------------



# =============================================================================
def orderDecays(): 
    for IDlog in range(nIDlog):
        #print nam[IDlog]+"  "+m[IDlog]+"  "+w[IDlog]
        for iDec in range(anDec[IDlog]-1):
            for jDec in range(anDec[IDlog]):
                if(jDec<=iDec): continue

                #then test & interchange
                if(aDecBR[IDlog][iDec]>aDecBR[IDlog][jDec]): continue  #no change

                #interchange the two:
                dum=aDecBR[IDlog][iDec]; aDecBR[IDlog][iDec]=aDecBR[IDlog][jDec]; aDecBR[IDlog][jDec]=dum;
                #dum=aDecWidth[IDlog][iDec]; aDecWidth[IDlog][iDec]=aDecWidth[IDlog][jDec]; aDecWidth[IDlog][jDec]=dum;
                dum=aDecnJoin[IDlog][iDec]; aDecnJoin[IDlog][iDec]=aDecnJoin[IDlog][jDec]; aDecnJoin[IDlog][jDec]=dum;
                dum=aDecnDau[IDlog][iDec]; aDecnDau[IDlog][iDec]=aDecnDau[IDlog][jDec]; aDecnDau[IDlog][jDec]=dum;
                dum=aDecDau[IDlog][iDec]; aDecDau[IDlog][iDec]=aDecDau[IDlog][jDec]; aDecDau[IDlog][jDec]=dum;
                



# =============================================================================
def printUseInfo(): 
    print("-----------------------------------------------------------")
    print("-----------------------------------------------------------")
    print("-----------------------------------------------------------")
    print("")
    
def printHead(): 
    print(" Based on: "+filen+"  "+COMMENT_ON_TOP)

def printNBs(): 
    print("")
    print("  [NB: using BRmin="+str(BRmin)+"]")
    print("  [NB: decay to e.g. C1+ and C1- are summed]")
    print("  [NB: two first generations are summed, name of first is used, e.g. eL is eL+muL]")
    #print "  [NB: BRs are ordered in size]\n"
    print("  [ BRs are ordered according to size ]\n")
#   ---


# =============================================================================
def printDecay():
    global donotprint #added
    page = []
    if(SHORT):
        donotprint+=["dL","b1","b2","t1","t2","dR","snue","snut","H","A","H+"]
    linebefore=["gl","eL","N1","h"]
    donotprint += shownotpart  #Oct7

    #print 'A', nIDlog
    for IDlog in range(nIDlog):
        if(nam[IDlog] in donotprint): continue;
        if(nam[IDlog] in PSKIP): continue
        
        if float(m[IDlog]) > MAXMASSPRINT: continue
        if(len(onlyprint) != 0 and nam[IDlog] not in onlyprint): continue
        #for ii in range(NSPACE): print ""
        for ii in range(NSPACE): page.append("")
        if(PRINTLINE and nam[IDlog] in linebefore):
            #print "-----------------------------------------------------------"
            page.append("-----------------------------------------------------------")
        if(NOMASSESABOVE10TeV):
            tparent="  "+nam[IDlog].ljust(4)+" ("+massEq+"%7.2f)   "  % float(m[IDlog])
        else:
            tparent="  "+nam[IDlog].ljust(4)+" ("+massEq+"%8.2f)  "  % float(m[IDlog])

        tpar=nam[IDlog].ljust(4)
        tparmass=massEq+"%7.2f"  % float(m[IDlog])
        if(PRINTPARMASS==1):
            tparent="  "+tpar+" ("+tparmass+") "
        elif(PRINTPARMASS==2):
            tparent=" ("+tparmass+")   "+tpar
        tparent0=tparent

        #print "QuarkMode = %i" %(QuarkMode)

        #----- begin(new code)
        zname=nam[IDlog]
        tcomp=""
        if(zname=="ta1"): tcomp=" [L=%4.2f  R=%4.2f]" %(tau1L,tau1R)
        if(zname=="ta2"): tcomp=" [L=%4.2f  R=%4.2f]" %(tau2L,tau2R)
        if(zname=="t1"): tcomp=" [L=%4.2f  R=%4.2f]" %(t1L,t1R)
        if(zname=="t2"): tcomp=" [L=%4.2f  R=%4.2f]" %(t2L,t2R)
        if(zname=="b1"): tcomp=" [L=%4.2f  R=%4.2f]" %(b1L,b1R)
        if(zname=="b2"): tcomp=" [L=%4.2f  R=%4.2f]" %(b2L,b2R)
        if(zname=="C1"): tcomp=" [w=%4.2f  h=%4.2f](L)  [w=%4.2f  h=%4.2f](R)" %(C1wL,C1hL,C1wR,C1hR)
        if(zname=="C2"): tcomp=" [w=%4.2f  h=%4.2f](L)  [w=%4.2f  h=%4.2f](R)" %(C2wL,C2hL,C2wR,C2hR)

        if PRINTCOMPOSITION == 1: 
            if(zname=="N1"): tcomp=" [b=%4.2f  w=%4.2f  h=%4.2f]" %(N1b,N1w,N1h)
            if(zname=="N2"): tcomp=" [b=%4.2f  w=%4.2f  h=%4.2f]" %(N2b,N2w,N2h)
            if(zname=="N3"): tcomp=" [b=%4.2f  w=%4.2f  h=%4.2f]" %(N3b,N3w,N3h)
            if(zname=="N4"): tcomp=" [b=%4.2f  w=%4.2f  h=%4.2f]" %(N4b,N4w,N4h)

        if PRINTCOMPOSITION == 2: 
            if(zname=="N1"): tcomp=" [b=%4.2f  w=%4.2f  hU=%4.2f  hD=%4.2f]" %(N1b,N1w,N1h2,N1h1)
            if(zname=="N2"): tcomp=" [b=%4.2f  w=%4.2f  hU=%4.2f  hD=%4.2f]" %(N2b,N2w,N2h2,N2h1)
            if(zname=="N3"): tcomp=" [b=%4.2f  w=%4.2f  hU=%4.2f  hD=%4.2f]" %(N3b,N3w,N3h2,N3h1)
            if(zname=="N4"): tcomp=" [b=%4.2f  w=%4.2f  hU=%4.2f  hD=%4.2f]" %(N4b,N4w,N4h2,N4h1)


        tcomp0=tcomp
        #----- end(new code)

        #tline=tparent
        #print tparent,

        #if(nam[IDlog]=="N1"):  #changed 2009, jan26 to line below (works with GMSB as well)
        if(len(list(range(anDec[IDlog])))==0): 
            if(PRINTCOMPOSITION):
                tstable=" <stable>                 "+tcomp
            else:
                tstable=" <stable>"
            #tline+=tstable
            #print tstable
            #if(BRPOS==1): print "           "+tparent+tstable
            #if(BRPOS==2): print tparent+tstable+"           "
            if(BRPOS==1): page.append("           "+tparent+tstable)
            if(BRPOS==2): page.append(tparent+tstable+"           ")
            continue
        
        totBR=0
        for iDec in range(anDec[IDlog]):
            if(aDecBR[IDlog][iDec]<BRmin): continue
            if(iDec>=NDECMAX): continue
            totBR+=aDecBR[IDlog][iDec]
            #totWidth+=aDecWidth[IDlog][iDec]

            if(not PRINTFULLPARENTEACHLINE and iDec!=0):
                #tparent="                        "
                #tcomp=""
                tparent="                   " #april20
                tcomp=""
                
            tline=tparent  
            
            #tBR="BR= %6.4f "  % aDecBR[IDlog][iDec]
            tBR="BR= %7.5f  "  % aDecBR[IDlog][iDec]
            tline+=tBR

            partWidth = aDecBR[IDlog][iDec] * w[IDlog]
            if partWidth < 1e-4:
                tWidth = "w=%9.2e  " %partWidth
            else: 
                tWidth = "w=%9.6f  " %partWidth

            #tWidth = "w=%3.1f" %partWidth  #Hack2010,apr20 (large widths)
            #tWidth = " w=%3.1f" %partWidth  #Hack2010,apr20
            
            tarrow=""
            if(PRINTARROW): tarrow=" --> "
            tline+=tarrow
                
            
            #rewrite in terms of dau1 and daurest
            tdau1=""
            tdaurest=""
            tdaumass1=""

            
            
            for iDau in range(aDecnDau[IDlog][iDec]):
                dum=""

                if(iDau==0):
                    tdau1=fname(aDecDau[IDlog][iDec][iDau],QuarkMode,LeptonMode).ljust(4)
                    #tdaumass1="         ("+massEq+"%4.0f)" % float(m[aIDlog[aDecDau[IDlog][iDec][0]]])
                    tdaumass1="         ("+massEq+"%6.1f)" % float(m[aIDlog[aDecDau[IDlog][iDec][0]]])  #Jan8

                else:
                    tdaurest+=" + "+fname(aDecDau[IDlog][iDec][iDau],QuarkMode,LeptonMode).ljust(4)
                    

            # NOW CONSTRUCT THE LINE TO PRINT --------
            tline=tparent+tarrow

            tBRandW = tBR
            if PRINTWIDTH: tBRandW = tBR+tWidth

            if(BRPOS==1):
                tline=tBRandW+tline
            else:
                tline=tline+tBRandW
                
            if(PRINTMASSESOFDAUGHTERS==0):
                tline+=tdau1+tdaurest
            elif(PRINTMASSESOFDAUGHTERS==1):
                tline+=tdaumass1+tdau1+tdaurest
            elif(PRINTMASSESOFDAUGHTERS==2):
                tline+=tdau1+tdaumass1+tdaurest
            elif(PRINTMASSESOFDAUGHTERS==3):
                tline+=tdau1+tdaurest+tdaumass1        
                
            tline+=tcomp

            tct=""
            if(ct[IDlog]>1e-4): tct="  <---- ct=%7.1emm  SEMISTABLE?" %(ct[IDlog]*1000.) #show in mm
            tline+=tct

            #print tline
            page.append(tline)
            
        dum="                    [= %6.4f]" % totBR
        if(PRINTSUM and BRPOS==1): page.append(" [= %6.4f]" % totBR)
        if(PRINTSUM and BRPOS==2): page.append(hline(" ",len(tparent+tarrow))+" [= %6.4f]" % totBR)
            
    page.append("")
    
    if PRINTTOFILE:
        f = open(filenprint,"w")
        for line in page: f.write("%s\n" %line)
        f.close()
    if not PRINTTOFILEONLY:
        for line in page: print(line)
        pass

    #print 'here', PRINTTOFILEONLY

# =============================================================================
def printMasses(mode=MODEMASSES):
    # mode=0: print masses (detailed)
    # mode=1: print names&masses on one line (ordered)
    # mode=2: print names on one line (ordered)
    donotprint=["t"]
    if(SHORT):
        donotprint+=["dL","b1","b2","t1","t2","dR","snue","snut","H","A","H+"]

    #printarray=[]
    printarray2=[]  # 2012-06-02
    for IDlog in range(nIDlog):
        zname=nam[IDlog]
        if float(m[IDlog]) > MAXMASSPRINT: continue  # 2013-12-21
        #if(zname==""): print "hmmmmm"+str(IDlog)
        if(zname in donotprint): continue;
        zline=" %7.2f  " % float(m[IDlog])
        zType=pType(zname)
        if(zType==0): zline+="     %-4s    " % zname    # others (gravitino, ..)
        if(zType==5): zline+="%-2s           " % zname    # higgs
        if(zType==1): zline+="   %-2s        " % zname    # gluino
        if(zType==2): zline+="   %-2s        " % zname    # squark
        if(zType==4): zline+="      %-2s     " % zname    # neutralino/chargino
        if(zType==3): zline+="         %-4s"   % zname    # sleptons
       
        zline2=""
        if(zname=="N1"): zline2=" [b=%4.2f  w=%4.2f  h=%4.2f]" %(N1b,N1w,N1h)
        if(zname=="N2"): zline2=" [b=%4.2f  w=%4.2f  h=%4.2f]" %(N2b,N2w,N2h)
        if(zname=="N3"): zline2=" [b=%4.2f  w=%4.2f  h=%4.2f]" %(N3b,N3w,N3h)
        if(zname=="N4"): zline2=" [b=%4.2f  w=%4.2f  h=%4.2f]" %(N4b,N4w,N4h)
        if(zname=="C1"): zline2=" [w=%4.2f  h=%4.2f](L)  [w=%4.2f  h=%4.2f](R)" %(C1wL,C1hL,C1wR,C1hR)
        if(zname=="C2"): zline2=" [w=%4.2f  h=%4.2f](L)  [w=%4.2f  h=%4.2f](R)" %(C2wL,C2hL,C2wR,C2hR)
        if(zname=="ta1"): zline2=" [L=%4.2f  R=%4.2f]" %(tau1L,tau1R)
        if(zname=="ta2"): zline2=" [L=%4.2f  R=%4.2f]" %(tau2L,tau2R)
        if(zname=="t1"): zline2=" [L=%4.2f  R=%4.2f]" %(t1L,t1R)
        if(zname=="t2"): zline2=" [L=%4.2f  R=%4.2f]" %(t2L,t2R)
        if(zname=="b1"): zline2=" [L=%4.2f  R=%4.2f]" %(b1L,b1R)
        if(zname=="b2"): zline2=" [L=%4.2f  R=%4.2f]" %(b2L,b2R)
        zline+=zline2

        #print zline
        #printarray.append(zline)
        printarray2.append([float(m[IDlog]), zline])  # 2012-06-02 # to avoid C1 being accidentally put lighter than N1 when equal to second decimal

        # onelines:
        #oneline.append(zname

    # Then sort according to masses
    #printarrayS=printarray
    #printarrayS.sort(None,None,True)  #reverse

    printarray2.sort()
    printarray2.reverse()
    printarrayS = []
    for zzz in printarray2: printarrayS.append(zzz[1])

    if mode==0:
        for iL in range(len(printarrayS)):
            print(printarrayS[iL])

    # adhoc for oneline Apr26
    if mode>0:
        print('%s  ' %(scenID), end=' ') 
        for iL in range(len(printarrayS)):
            zm = printarrayS[iL].split()[0]
            zp = printarrayS[iL].split()[1]

            if zp in ['h','H','A','H+','grav']: continue  # skip higgses and gravitino
            

            if mode == 1:
                print('%s(%s) ' %(zp,zm), end=' ')
            if mode == 2:
                print('%s ' %(zp), end=' ')


# ============================================================================= END
# ============================================================================= FUNCTIONS
# =============================================================================


# ===============================================================================================
# - Here ends the traditional 'StructureISAWIGfiles.py' which produces derived files
# - The 5-6 lines above could be put inside a common 'if DERIVED' to signal clearly that
#   they are only run in the derived mode
# - (The old structure is now sufficiently 'objectified'. There are lots of global variables,
#    but that shouldn't pose a problem.) 
# - Below are newer stuff which uses the decay info already read in to
#   - inspect 1legs, in particular di-leptonic legs (with Per & David)
#   - inspect 2legs w/prospino: to address the so-called eigenvalue method
# 
# ===============================================================================================
# ===============================================================================================


    
# ===============================================================================================
# ===============================================================================================
# ===============================================================================================
# ===============================================================================================
# ===============================================================================================
# ====================  BUILD NEW STRUCTURE =====================================================
# ===============================================================================================
# ===============================================================================================
# ===============================================================================================
# ===============================================================================================
# ===============================================================================================

def calcmll(a,b,c):
    #print a, b, c
    return math.sqrt((a*a-b*b)*(b*b-c*c))/b

def checkmllEtc():
    pN2 = part['N2']
    #pN2.showDecays()
    mN2 = part['N2'].mass
    mN1 = part['N1'].mass
    for SL in ['lR','lL','eR','eL']: 
        if pN2.BR(SL) < 0.2: continue
        msl = part[SL].mass
        mll = calcmll(mN2,msl,mN1)
        print("checkmllEtc:  mll(%s) = %6.2f    (N2,%s,N1) = (%.1f, %.1f, %.1f)   %s" %(SL,mll, SL, mN2,msl,mN1, filen))
        

# =============================================================================
# CONSTRUCTING dict part={..}
#              and mass-ordered list partMO=[..]
#              as well as includes Z/W-splitting in decay list
# Based solely on:  aIDgen[iP], nam[iP], m[iP], t[iP], anDec[iP], aDecBR, aDecnDau, aDecDau
#              and: SPLIT_ZWetc
# =============================================================================
def importFullParticleInfo_all(): 
    for iP in range(len(aIDgen)):
        importFullParticleInfo(iP)
    
# =============================================================================
def importFullParticleInfo(iP): 
    zpart = Particle(aIDgen[iP], nam[iP], m[iP], t[iP])
    part[nam[iP]] = zpart

    if zpart.name not in ['h','H','H+','A','t']: partMO.append(zpart)



    # Loop over decays (for this particle)
    #part[nam[iP]].NdecayList = anDec[iP]  #re-set below
    for iD in range(anDec[iP]):
        arg1=aDecBR[iP][iD]; arg2=aDecnDau[iP][iD]; arg3=aDecDau[iP][iD]
        zdecay = Decay(arg1,arg2,arg3)
        # print zdecay.dauName, zdecay.n, zdecay.SM[0]

        dau0 = zdecay.SM[0]
        
        if dau0 not in SPLIT_ZWetc: 
            part[nam[iP]].decayList.append(zdecay)
        else:
            BR0 = zdecay.BR
            if dau0 == 'Z':
                zdecay = Decay(arg1,arg2,arg3); zdecay.SM[0] = 'Z->ll'; zdecay.BR = BR0 * g.PDG.Ztoll ; part[nam[iP]].decayList.append(zdecay)
                zdecay = Decay(arg1,arg2,arg3); zdecay.SM[0] = 'Z->TT'; zdecay.BR = BR0 * g.PDG.ZtoTT ; part[nam[iP]].decayList.append(zdecay)
                zdecay = Decay(arg1,arg2,arg3); zdecay.SM[0] = 'Z->vv'; zdecay.BR = BR0 * g.PDG.Ztovv ; part[nam[iP]].decayList.append(zdecay)
                zdecay = Decay(arg1,arg2,arg3); zdecay.SM[0] = 'Z->qq'; zdecay.BR = BR0 * g.PDG.Ztoqq ; part[nam[iP]].decayList.append(zdecay)
                zdecay = Decay(arg1,arg2,arg3); zdecay.SM[0] = 'Z->bb'; zdecay.BR = BR0 * g.PDG.Ztobb ; part[nam[iP]].decayList.append(zdecay)

            elif dau0 == 'W':
                zdecay = Decay(arg1,arg2,arg3); zdecay.SM[0] = 'W->lv'; zdecay.BR = BR0 * g.PDG.Wtolv ; part[nam[iP]].decayList.append(zdecay)
                zdecay = Decay(arg1,arg2,arg3); zdecay.SM[0] = 'W->Tv'; zdecay.BR = BR0 * g.PDG.WtoTv ; part[nam[iP]].decayList.append(zdecay)
                #zdecay = Decay(arg1,arg2,arg3); zdecay.SM[0] = 'W->QQ'; zdecay.BR = BR0 * g.PDG.WtoQQ ; part[nam[iP]].decayList.append(zdecay)
                zdecay = Decay(arg1,arg2,arg3); zdecay.SM[0] = 'W->qq'; zdecay.BR = BR0 * g.PDG.Wtoqq ; part[nam[iP]].decayList.append(zdecay)

            # Sep16 hmm
            elif dau0 == 't':
                #zdecay = Decay(arg1,arg2,arg3); zdecay.SM[0] = 't->bW'; zdecay.BR = BR0 * g.PDG.Wtolv ; part[nam[iP]].decayList.append(zdecay)
                zdecay = Decay(arg1,arg2,arg3); zdecay.SM[0] = 't->blv'; zdecay.BR = BR0 * g.PDG.Wtolv ; part[nam[iP]].decayList.append(zdecay)
                zdecay = Decay(arg1,arg2,arg3); zdecay.SM[0] = 't->bTv'; zdecay.BR = BR0 * g.PDG.WtoTv ; part[nam[iP]].decayList.append(zdecay)
                zdecay = Decay(arg1,arg2,arg3); zdecay.SM[0] = 't->bqq'; zdecay.BR = BR0 * g.PDG.Wtoqq ; part[nam[iP]].decayList.append(zdecay)
                
            else:
                print("Not yet implemented for ",dau0)

    part[nam[iP]].NdecayList = len(part[nam[iP]].decayList)  # need to be (re)set manually due to 'decay-expansion' above  (~8hour debugging...)


    # testing Apr17 (result is ok & complete)
    if 0 and zpart.name in ['b2']:
        if 0: 
            zpart.decayList.pop(0)
            zpart.decayList.pop(0)
            zpart.decayList.pop(0)
            zpart.decayList.pop(0)
            zpart.NdecayList -= 4
        if 0: 
            aaa = zpart.decayList[1]
            zpart.decayList[1] = zpart.decayList[6]
            zpart.decayList[6] = aaa
        for dec in zpart.decayList:
            print(dec.SUSY, "   ",dec.SM, "  ",dec.BR)

# =============================================================================
def testLSP(): 
    #if not (partMO[0].name=='grav' and partMO[1].name=='N1' or partMO[0].name=='N1'):
    if not (partMO[0].name=='grav' and partMO[1].name==LSPid or partMO[0].name==LSPid):
        txt = partMO[0].name
        for iP in range(1,len(partMO)): txt += " < "+partMO[iP].name
        #print 'N1 not LSP: ',filen
        print('%s not LSP (but was assumed to be): ' %(LSPid, filen))
        print(txt)
        exit()

# =============================================================================
#def isLSP(name): return (name=='N1')
def isLSP(name): return (name==LSPid)

# =============================================================================
def MakeTableToOrderAccordingToMasses():
    for ii in range(len(aIDgen)):
        ziterMO.append([float(m[ii]), ii])
    ziterMO.sort()
    ziterMO.reverse()
    for ii in range(len(ziterMO)):
        iterMO.append(ziterMO[ii][1])
        namMO.append(nam[ziterMO[ii][1]])

# =============================================================================
def NextDaughter(ziter,casc):

    ziter += 1
    tpre = 2*ziter*" "+str(ziter)
    tpre = '%-10s' %tpre
    # daughters

    # 1 First find the last SUSY particle (Particle()) in cascade (leg)
    if casc.Last() == '--': print(60*"c")  # TEST (Never happens) DEBUG
    pLast = part[casc.Last()] 
    
    # 2 Store & Return when reach N1 (decay is then fully developed)
    if isLSP(casc.Last()):
        cascStore.Add(casc)
        return
    
    # 3 Loop over possible decays for the (so far) last particle: cann NextDaughter (recursive)
    for iD in range(pLast.NdecayList):
        decay = pLast.decayList[iD]
        # - Find name of next sparticle 
        # pNewName = decay.SUSY[0]
        
        # - Make a copy cascNew of casc. CascNew is built on below, casc is the one to start of when recursed back
        cascNew = casc.myCopy()

        # - cascNew is built on with the next sparticle in line (decay.SUSY[0]) 
        # cascNew.Add(decay.SUSY[0], iD, decay.SM, decay.BR)
        cascNew.Add(decay.SUSY[0], decay.SM, decay.BR)

        # drop if BR too low
        if cascNew.BRtot<g.BRmin1leg:
            continue 

        # - If the last particle is LSP, then store the resulting chain and continue/recurse back  #
        if isLSP(cascNew.Last()):
            cascStore.Add(cascNew)                #<--- NB(ugly): stores in cascStore
            # if cascNew.isDilep:
            #    cascStoreDilep.append(cascNew)
            continue

        NextDaughter(ziter,cascNew)
        
# =============================================================================
# =============================================================================
def analysis1leg(): 
    # print 10*"CASCADES  "

    # STEP 1: BUILD CASCSTORE
    # For all particles in g.casc_starttypes store all chains with BR>g.BRmin1leg in cascStore
    for itp in range(len(namMO)):
        pnam = namMO[itp]
        p = part[pnam]

        #print "%2i  %4.f  %s" %(itp,p.mass,p.name) #verbosity

        #accept a certain selection of cascade starts types (for now)
        if pnam not in g.casc_starttypes: continue


        # For given start particle, go through decays and build a chain

        # casc = Cascade(pnam, itp)
        
        casc = Cascade(pnam) #starts a cascade (leg) with pnam as initiator  libISAWIG::Cascade  # class has 5-600 lines ... init is simple
        NextDaughter(0,casc)  #<-- NB stores (hidden) in cascStore  # 50 lines
        # NextDaughter() recursively finds all cascades with total BR > g.BRmin1leg,
        # then when LSP has been for a decay chain stores the chains in CascStore

    # at this stage have already isDilep set for Z->ll, not ll, not sleLL, not snuLL


    # STEP 2: REFINE CASCSTORE
    if 1:

        DEBUG = 1  #D=0   #To save list for each operation
        JOIN  = 1  #D=0

        if DEBUG: cascStore.Save("debug_casc1_all_0")

        
        # OPERATION: nue/num/nut->v, ta->T, e->l, m->l, pho->y ;
        # OPERATION: combine multiple SM entries in a decay into 1 entry, i.e. to give: pipipi, lv, ud, ..
        cascStore.SimplifyNotation();
        if JOIN: cascStore.Duplicates('join');
        if DEBUG: cascStore.Save("debug_casc2_all_SimplifyNotation")

        # OPERATION: uL,dL,... -> qL  (not b and t)
        cascStore.JoinQs(pJoin="qL", start=1);  # start=1: do not change the first sparticle (due to xsec consideratoins)
        if JOIN: cascStore.Duplicates('drop'); # drop if head, else join
        if DEBUG: cascStore.Save("debug_casc3_all_JoinQs_qL")

        # OPERATION: uR,dR,... -> qR  (not b and t)
        cascStore.JoinQs(pJoin="qR", start=1);  # start=1: do not change the first sparticle (due to xsec consideratoins)
        if JOIN: cascStore.Duplicates('drop'); # drop if head, else join
        if DEBUG: cascStore.Save("debug_casc4_all_JoinQs_qR")

        # OPERATION: qL,qR,b1,b2 -> SQ (only t1 and t2 left)
        cascStore.JoinQs(pJoin="SQ", start=1); # start=1: do not change the first sparticle (due to xsec consideratoins)
        if JOIN: cascStore.Duplicates('join');
        if DEBUG: cascStore.Save("debug_casc5_all_JoinQs_SQ")


        # OPERATION: ud->QQ, 'Z->qq'->'Z->QQ'
        cascStore.JoinQQ("QQ");
        if JOIN: cascStore.Duplicates('join');
        if DEBUG: cascStore.Save("debug_casc6_all_JoinQQ_QQ")

        # Final
        # OPERATION
        cascStore.Duplicates('join');
        if DEBUG: cascStore.Save("debug_casc7_all_Duplicates_join")
        cascStore.Classify(categAll)

        cascStore.Save("casc_all_"+filen2)

        cascStoreDilep = cascStore.Dileptons(g)
        cascStoreDilep.Save("casc_dilep_"+filen2)


        
    # BELOW 'OLD'
    if 0: 
        cascStore.Show()
        cascStoreDilep = cascStore.Dileptons(g)
        cascStoreDilep.Show()

    if 0:
        cascStore.SimplifyNotation(); #cascStore.Show()
        cascStore.Classify(categAll)
        print(10*"--", "cascClassified")
        cascStore.Show()
        cascStoreDilep2 = cascStore.Dileptons(g)
        print(10*"--", "Dilep2")
        cascStoreDilep2.Show()

    if 0:
        cascStoreDilep.JoinQs("qL"); #cascStoreDilep.Show()
        cascStoreDilep.JoinQs("qR"); #cascStoreDilep.Show()
        cascStoreDilep.JoinQs("SQ"); #cascStoreDilep.Show()
        cascStoreDilep.JoinQQ("QQ"); #cascStoreDilep.Show()

        cascStoreDilep.Duplicates('join'); #cascStoreDilep.Show()
        cascStoreDilep.SimplifyNotation(); #cascStoreDilep.Show()
        cascStoreDilep.Save("casc_dilepton_"+filen2)

    if 0:
        cascStore.JoinQs("qL"); #cascStore.Show()
        cascStore.JoinQs("qR"); #cascStore.Show()
        cascStore.JoinQs("SQ"); #cascStore.Show()
        cascStore.JoinQQ("QQ"); #cascStore.Show()

        cascStore.Duplicates('join'); #cascStore.Show()
        cascStore.SimplifyNotation(); #cascStore.Show()
        cascStore.Save("casc_all_"+filen2)
# =============================================================================
# =============================================================================
def printPrOXsec(xT,x,xInput=[],form='%12.7f'):
    for ikey in range(len(xT)):
        key = xT[ikey]
        print(' %8s   ' %(key) + form %(x[key]))

# =============================================================================
def getPrOXsec(fn):
    #print "<getPrOXsec>"
    xT = []
    x  = {}
    xInput = {}
    if not os.path.exists(fn):
        print('WARNING: getPrOXsec: nonexistent file: %s' %(fn))
        return (xT, x, xInput)
    f=open(fn); l=f.readlines(); f.close()
    for iL in range(len(l)):
        if l[iL].startswith("#"):
            # Look for relevant extra info
            # (some are probably too prospino-specific)
            line = l[iL][1:].strip()
            ww = line.split()
            nw = len(ww)
            if line.startswith('Collider') and nw>=3: xInput['collider'] = line.split()[2]
            if line.startswith('COM Energy') and nw>=4: xInput['comenergy'] = line.split()[3]
            if line.startswith('Calculation') and nw>=3: xInput['calculation'] = line.split()[2]
            if line.startswith('Variable') and nw>=3: xInput['variable'] = line.split()[2]
            if line.startswith('squark mass') and nw>=3: xInput['squarkmass'] = line[line.index(':')+1:].strip()
            if line.startswith('scale') and nw>=3: xInput['squarkmass'] = line[line.index(':')+1:].strip()
            if line.startswith('Units') and nw>=3: xInput['unit'] = line.split()[2]
        else:
            w=l[iL].split()
            if len(w)<3:
                print('WARNING: getPrOXsec: skipping illegal line (below):\n%s' %(l[iL]))
                continue
            key = w[0]+' '+w[1]
            xT.append(key)
            x[key] = float(w[2])
            
    return (xT, x, xInput)
# =============================================================================
def make1legs():
    #print "<make1legs>"
    # USES: 
    #   - cascStore


    # 1 BUILD CASCSTORE
    for itp in range(len(namMO)):
        pnam = namMO[itp]
        p = part[pnam]
        if pnam not in g.casc_starttypes: continue
        cascStore.initial.append(pnam)
        casc = Cascade(pnam) #starts a cascade (leg) with pnam as initiator
        NextDaughter(0,casc)  #<-- NB stores (hidden) in cascStore
    #   ---

    # 2 REFINE CASCSTORE
    DEBUG = 1  #D=0   #To save list for each operation  DEBUG option
    JOIN  = 0  #D=0
    
    if DEBUG: cascStore.Save("debug_1leg1_all_0")
    
    # OPERATION: nue/num/nut->v, ta->T, e->l, m->l, pho->y ;
    # OPERATION: combine multiple SM entries in a decay into 1 entry, i.e. to give: pipipi, lv, ud, ..

    cascStore.Status('<--- 1leg: Just after making all 1legs')
    cascStore.SimplifyNotation();
    if JOIN: cascStore.Duplicates('join');
    if DEBUG: cascStore.Save("debug_1leg2_all_SimplifyNotation")
    cascStore.Save(fn=g.it.next('pres','_')+g.scenID+'_an1legs_noReduction.txt')
    #cascStore.Status('<--- 1leg: Simplify notation (no reduction)')
   
    # OPERATION: uL,dL,... -> qL  (not b and t)
    cascStore.JoinQs(pJoin="qL", start=1);  # start=1: do not change the first sparticle (due to xsec consideratoins)
    if JOIN: cascStore.Duplicates('drop'); # drop if head, else join
    if DEBUG: cascStore.Save("debug_1leg3_all_JoinQs_qL")
    cascStore.Status('<--- 1leg: uL,dL,... -> qL')
    
    # OPERATION: uR,dR,... -> qR  (not b and t)
    cascStore.JoinQs(pJoin="qR", start=1);  # start=1: do not change the first sparticle (due to xsec consideratoins)
    if JOIN: cascStore.Duplicates('drop');  # drop if head, else join
    if DEBUG: cascStore.Save("debug_1leg4_all_JoinQs_qR")
    cascStore.Status('<--- 1leg: uR,dR,... -> qR')
    
    # OPERATION: qL,qR,b1,b2 -> SQ (only t1 and t2 left)
    cascStore.JoinQs(pJoin="SQ", start=1); # start=1: do not change the first sparticle (due to xsec consideratoins)
    if JOIN: cascStore.Duplicates('join');
    if DEBUG: cascStore.Save("debug_1leg5_all_JoinQs_SQ")
    cascStore.Status('<--- 1leg: qL,qR -> SQ')

    # OPERATION: eL,eR,mL,mR -> sl
    cascStore.Joinsl()
    if DEBUG: cascStore.Save("debug_1leg5b_all_Joinsl")
    cascStore.Status('<--- 1leg: eL,eR,mL,mR -> sl')
    
    # OPERATION: eL,mL,mL,mR,T1,T2 -> SL
    cascStore.JoinSL()
    if DEBUG: cascStore.Save("debug_1leg5b_all_JoinSL")
    cascStore.Status('<--- 1leg: sl,T1,T2 -> SL')
    
    # OPERATION: ud->QQ, 'Z->qq'->'Z->QQ'
    cascStore.JoinQQ("QQ");
    if JOIN: cascStore.Duplicates('join');
    if DEBUG: cascStore.Save("debug_1leg6_all_JoinQQ_QQ")
    cascStore.Status('<--- 1leg: diquarks(ud,uu,..) -> QQ')

    # OPERATION: Set any SM particle to 'x' (remove SM differences)
    if 0: 
        cascStore.eraseSM()
        if JOIN: cascStore.Duplicates('join');
        if DEBUG: cascStore.Save("debug_1leg6b_all_eraseSM")
        cascStore.Status("<--- 1leg: set all SM particle combinations (q,b,t,W,Z,l,T,v,ll,TT,vv,lv,Tv) to 'x' (i.e. only SUSY part matter)")
        

    
    # Final
    # OPERATION
    cascStore.Duplicates('join');
    if DEBUG: cascStore.Save("debug_1leg7_all_Duplicates_join")
    cascStore.Classify(categAll)
    
    #cascStore.Save("casc_all_1leg_"+filen2)
    cascStore.Save("casc_all_1leg/casc_all_1leg_"+filen2)
    cascStore.Status("<--- 1leg: final")

# =============================================================================
def make2legs():
    #print "<make2legs>"

    g.twolegs = TwoLegs(g.scenID)
    g.tl = g.twolegs

    # From Apr22: now uses the same as g.casc_starttypes
    #initialParticles = ['gl','uL','dL','sL','cL','uR','dR','sR','cR']
    #initialParticles = ['gl','uL','dL','uR','dR']
    initialParticles = g.casc_starttypes
    if VBscreen: print("initialParticles: ",str(initialParticles))

    g.twolegs.BRcut = g.BRmin2leg
    BRtotall = 0.
    BRtotreq = 0.
    BRdebug0 = 0.
    BRdebug1 = 0.
    BRdebug2 = 0.
    BRdebug0b = 0.
    BRdebug1b = 0.
    BRdebug2b = 0.

    # LOOP OVER ORDERED XSECS (drop [0] which is tot): p1 and p2
    for iOX in range(0,g.OX.N()):

        key = g.OX.xOT[iOX]

        # NB: the total cross section is kept in the same dict (dangerous)
        if key == 'tot':
            #print "hit key = tot,  iOX = %i" %(iOX)
            continue
        
        p1,p2 = key.split()

        # xrel0 is the relative xsec of the given subprocess (key here is e.g. 'sL cL', the two outgoing sparticles)
        xrel0 = g.OX.xO[key]/max(g.OX.tot,1e-9)

        BRtotall += xrel0
        # SKIP IF NOT REQUESTED INITIALPARTICLES
        if not (p1 in initialParticles and p2 in initialParticles): continue

        BRtotreq += xrel0
        BRdebug0 += xrel0

        # SKIP IF TOO LOW (REL) XSEC
        if xrel0 < g.BRmin2leg:
            # print "too low xsec: %f" %xrel0
            continue        

        BRdebug0b += xrel0


        # OUTER LOOP OVER LIST (STOP ONLY AT 1LEGS STARTING WITH P1)
        for i1 in range(cascStore.N()):
            casc1 = cascStore.list[i1]
            if casc1.SUSY[0] != p1: continue

            xrel1 = xrel0 * casc1.BRtot
            BRdebug1 += xrel1  
            #BRp1debug += casc1.BRtot

            # SKIP IF TOO LOW REL XSEC
            if xrel1 < g.BRmin2leg:
                #print "too low xsec: %f" %xrel1
                continue

            BRdebug1b += xrel1  
            

            # INNER LOOP OVER LIST (STOP ONLY AT 1LEGS STARTING WITH P2)
            #BRp2debug = 0.
            for i2 in range(cascStore.N()):
                casc2 = cascStore.list[i2]
                if casc2.SUSY[0] != p2: continue

                xrel2 = xrel1 * casc2.BRtot
                BRdebug2 += xrel2

                # SKIP IF TOO LOW REL XSEC
                if xrel2 < g.BRmin2leg: continue
                #BRp2debug += casc2.BRtot
                
                # IF ARRIVE HERE, THEN HAVE TWO LEGS WITH CORRECT INITIAL PARTICLES
                #                 AND SUFFICIENTLY GOOD XSEC x BR1 x BR2
                BRdebug2b += xrel2

                g.twolegs.Add(casc1.SUSY, casc2.SUSY, casc1.SM, casc2.SM, i1, i2, xrel2, g.OX.tot * xrel2, add2dict(casc1.SMp, casc2.SMp))
                


    # g.twolegs.sumBR() # done per Add
    g.twolegs.BRtotreq = BRtotreq
    g.twolegs.BRtotall = BRtotall
    #g.twolegs.Show(0,txt='blasted!!')
    if VBscreen: 
        print('For all below: full coverage is 1. They are xsec(1leg or 2leg) / xsec(total).')
        print('Coverage: BRdebug0:  %8.6f   Sum of BR of all 1legs for the given initial particles, see casc_starttypes()' %(BRdebug0))
        print('Coverage: BRdebug0b: %8.6f   As above, but enforcing cut (g.BRmin2leg > %8.6f)' %(BRdebug0b, g.BRmin2leg))
        print('Coverage: BRdebug1:  %8.6f   BR(p1,p2,cut) * sum[Dec(p1)]'  %(BRdebug1))
        print('Coverage: BRdebug1b: %8.6f   BR(p1,p2,cut) * sum[Dec(p1,cut)]'  %(BRdebug1b))
        print('Coverage: BRdebug2:  %8.6f   BR(p1,p2,cut) * sum[Dec(p1,cut)*sum[Dec(p2)]]'  %(BRdebug2))
        print('Coverage: BRdebug2b: %8.6f   BR(p1,p2,cut) * sum[Dec(p1,cut)*sum[Dec(p2,cut)]]'  %(BRdebug2b))


# =============================================================================
#def refine2legs():
#    # function moved inside analysis2legs()
#    pass

# =============================================================================
def show2legs():
    #print "<show2legs>"
    g.twolegs.Show(5)
    pass





# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
def study3L_make1legs():

    DEBUG = VBfile  #D=0   #To save lists to disk for each operation  DEBUG option


    # ===== 1  BUILD CASCSTORE
    for itp in range(len(namMO)):
        pnam = namMO[itp]
        #print " ... %s    %s" %(itp,pnam)
        p = part[pnam]
        if pnam not in g.casc_starttypes: continue
        cascStore.initial.append(pnam)
        casc = Cascade(pnam) #starts a cascade (leg) with pnam as initiator
        NextDaughter(0,casc)  #<-- NB stores (hidden) in cascStore
        
    if VBscreen: print("cascStore: %i  %s"  %(cascStore.NuniquePart(),str(cascStore.UniquePart())))
    if DEBUG: cascStore.Save("debug_1leg1")
    if VBscreen: cascStore.Status('  [1leg: Just made all 1legs]')


    # ===== 2  A FIRST SIMPLIFICAION OF SM OUTPUT:  e,m->l, and pi+pi+pi etc. (wino-LSP)
    #                                               NB: e,m -> l   [code is not yet set up to handle e and m separately]
    #                                               NB: occurence of pi's is likely to break the code (may need to generalise a few places)
    cascStore.SimplifyNotationSM();   
    if DEBUG: cascStore.Save("debug_1leg2_SimplifyNotationSM")
        

    # ===== 3  COLLAPSES INTERMEDIATE SPARTICLES, KEEP ONLY INITIAL AND FINAL  (the intermediate sparticles are not of interest for this study)
    cascStore.KeepOnlyFirstAndLastSparticle()
    cascStore.Duplicates('join');    
    if VBscreen: cascStore.Status("  [1leg: KeepFirstAndLast_Join]")
    if DEBUG: cascStore.Save("debug_1leg3_KeepFirstAndLast_Join")
    if VBscreen: mytimer.Add("make1legs: AFTER  KeepFirstAndLast")


    # ===== 4  FROM SM LIST MAKES SMp LIST (WHICH IS WHAT IS USED FROM NOW ONE)
    cascStore.FillSMp()
    if VBscreen: print("AFTER  1leg FillSMp(): N = %i" %cascStore.N())
    #if DEBUG: cascStore.Save("debug_1leg4_FillSMp")
    if DEBUG: cascStore.Save("debug_1leg4_FillSMp_mode2",mode=2)


    # ===== 5  EXPAND NEARLY ALL SM PARTICLES/COMBINATIONS INTO A FEW ANALYSIS-RELEVANT ONES, TYPICALLY SOME OF  l, q, b, T, y
    # PARTICLES/COMBINATIONS TO KEEP (typically some of l,q,b,T,y) SHOULD NOT BE IN THE LIST BELOW:
    # SMout is by default ['l']. Can be set on command line by ex: -SMout l,q,T
    #all_particlecombinations = ['t','h','W','Z','QQ','vv','TT','Tv','lv','ll', 'tb','u','d','g', 'uu','dd','ud','bb',  'b','T','y', 'q','l']
    all_particlecombinations = ['tt','t','h','H','A','H+','W','Z','QQ','vv','TT','Tv','lv','ll', 'tb','u','d','g', 'uu','dd','ud','bb',  'b','T','y', 'q','l']  #2011-01-03
    particles_to_transform = list(all_particlecombinations)
    for sm in SMout:
        if sm in particles_to_transform: particles_to_transform.remove(sm)

    #print particles_to_transform
        
    #particles_to_transform = ['t','h','W','Z','QQ','vv','TT','Tv','lv','ll', 'tb','u','d','g', 'uu','dd','ud','bb',  'T','y']
    while cascStore.ExpandSMp(mode=particles_to_transform):  pass   
    if VBscreen: print("AFTER  1leg ExpandSMp: N = %i" %cascStore.N())
    #if DEBUG: cascStore.Save("debug_1leg5_ExpandSMp")
    if DEBUG: cascStore.Save("debug_1leg5_ExpandSMp_mode2",mode=2)
    

    # ===== 6  TIDY SMp LIST; DROP REMAINING EMPTY, KEEP ONLY A FEW ANALYSIS-RELEVANT ONES
    #cascStore.TidySMp(keep=['l','T','q','b','y'], drop=['v'])
    cascStore.TidySMp(keep=['l','T','q','b','y'], drop=['v'], warntxt=scenID)
    if VBscreen: mytimer.Add("make1legs: After 1leg Tidy: N = %i" %(cascStore.N()),-1)
    #if DEBUG: cascStore.Save("debug_1leg6_TidySMp")
    if DEBUG: cascStore.Save("debug_1leg6_TidySMp_mode2",mode=2)


    # ===== 7  FINAL JOIN OF DUPLICATES FOR 1LEG
    cascStore.Duplicates('join','SMp');  # <-- need to join SU+SMp (not SM) [Is this wrong? Jan10]
    if VBscreen: mytimer.Add("make1legs: After 1leg Duplicates(join): N = %i" %(cascStore.N()),-1)
    #if DEBUG: cascStore.Save("debug_1leg7_Final")
    if DEBUG: cascStore.Save("debug_1leg7_Final_mode2",mode=2)

    

# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# Sep16: study3L() started out as identical copy of analysis2leg() [which is kept unchanged for other studies]

def study3L(mode=5):

    DEBUG = VBfile #True

    ########## A - get xsecs
    if VBscreen: print("\n\n A - get xsecs")
    #if g.VB >=2: print "INFO: getPrOXsec " + 30*"-"
    g.ox = getPrOXsec(g.PrOX_filename)
    g.OX = OX(g.ox[0], g.ox[1], g.ox[2])  # OX-object (useful)
    mytimer.Add("getPrOXsec")
    if VBscreen: print("Number of subprocesses: %i" %(g.OX.N()))
    if VBscreen: print("Total xsec: %.3f pb" %(g.OX.tot))


    ########## B - make1legs
    if VBscreen: print("\n\n B - make1legs")
    study3L_make1legs()
    #print "make1legs: %i onelegs" %(cascStore.N())
    if VBscreen: mytimer.Add("make1legs: %i onelegs" %(cascStore.N()),-1)
    #cascStore.AllPart_Coverage()
    

    ########## C - make2legs
    if VBscreen: print("\n\n C - make2legs")
    make2legs()
    if VBscreen: mytimer.Add("make2legs: %i twolegs" %(g.twolegs.N()),-1)
    if DEBUG: g.twolegs.Show(fn='debug_2leg1_untouched.txt')


    ########## D - refining (very little needed, most is done in 1leg)
    if VBscreen: print("\n\n D - refining")
    #g.twolegs.Show(0,txt='<--- Just after making all 1leg x 1leg combinations')
    g.twolegs.Duplicates('join',['SMp']) 
    if DEBUG: g.twolegs.Show(fn='debug_2leg2_Duplicates_Join.txt')
    if VBscreen: mytimer.Add("make2legs: %i twolegs after Duplicates_Join" %(g.twolegs.N()),-1)
    g.twolegs.orderWithSMp()
    if VBscreen: print("\n(In iPython you can look at results with g.twolegs.ShowSMp())")


    ########## E - output
    if saveSMp:
        #g.twolegs.Show(fn=fnSMp, SM=['SMp'], SMpHead2=1)
        g.twolegs.ShowSMp(fn=fnSMp, keys=SMout)
        if VBscreen: print("Saved SMp output to %s" %(fnSMp))

    if showSMp:
        #g.twolegs.Show(SM=['SMp'], SMpHead2=1)
        g.twolegs.ShowSMp(keys=SMout)
        

    ########## F - reduction path
    
    if SMout_reductionPath:
        
        if VBscreen: print("\n\n F - reductionPath")
        SMout0 = list(SMout)  # store the original

        for iRed in range(-1,len(SMout_reductionPath)):
            thisSMred = SMout_reductionPath[iRed]


            if iRed > -1: 
                while g.twolegs.ExpandSMp(mode=[thisSMred]): pass
                g.twolegs.Duplicates('join',check=['SMp'])
                if thisSMred in SMout: 
                    SMout.remove(thisSMred)
                else:
                    print('ERROR::isaPlay  trying to remove thisSMred: %s from SMout: %s' %(thisSMred, SMout))

            pid = "_"
            for p in SMout: pid += p                

            g.twolegs.orderWithSMp()
            #g.twolegs.ShowSMp(fn=fnSMp+pid, keys=SMout0)
            g.twolegs.ShowSMp(fn=fnSMp + mode_casc_starttypes + pid, keys=SMout0)

            print("Saving: %-30s   (%4i twolegs)" %(fnSMp+pid, g.twolegs.N()))

    



    ########## F - MAKE PLOTS OF SMp (PROBABLY NOT WORKING)
    if TWOLEG_HIST:  #here assuming ROOT to be available (e.g. ATLAS release sourced)
        g.twolegs.FillH()
        fn0 = '2legs_final'
        cc,tleg = g.twolegs.ShowH(which=[-1,0,1,2,3,4], fn=[fn0+'.pdf', fn0+'.eps'], silent=True)
        return cc,tleg  #to preserve for interactive python
    


# =============================================================================




# =============================================================================
# ============================================================================= End of functions
# ============================================================================= for various
# ============================================================================= casc implementations
# =============================================================================


# =============================================================================
# ============================================================================= BEGIN EXECUTION
# =============================================================================


mytimer = MyTimer()



ReadMassesBRsMixingparsEtc()
orderDecays()
MakeTableToOrderAccordingToMasses()    #(not used for derived)

# ============================================================================= FOR PRODUCING 
# ============================================================================= DERIVED FILES
# ============================================================================= 
#print 'hei', PRINTDECAY
if(PRINTDECAY and PRINTUSEINFO): printUseInfo()
if(PRINTHEAD): printHead()
if(PRINTDECAY and PRINTUSEINFO): printNBs()

if(PRINTDECAY): printDecay()
if(PRINTMASSES):  printMasses()





# ============================================================================= FOR ANALYSING 
# ============================================================================= 1LEGS OR 2LEGS
# ============================================================================= 

if not (ONELEGANALYSIS or TWOLEGANALYSIS or STUDY3L or checkmll):
    sys.exit()



# ============================================================================= ADDITIONAL GLOBAL VARIABLES

# --- The main global variable henceforth
g = FreeObject()
g.VB = 2
g.scenID = scenID
#g.PRES = 'pres_'  #will be prepended to the most important file-savings in twoleg-mode (for preserve)

g.it = Incrementer(0)


# ### PARTICLE DATA GROUP
g.PDG = PDG  # from libISAWIG.py
# ########################## 

# list with important runtime messages
NOTIFY = []

# list with important codetime messages
CRUCIALS = []

part = {}
partMO = []

# ============================================================================= 




if (TWOLEGANALYSIS or STUDY3L) and BRmin2leg != BRmin1leg:
    NOTIFY.append("Changing BRmin1leg to be consistent with BRmin2leg:  %7.1e ->  %7.1e" %(BRmin1leg, BRmin2leg))
    BRmin1leg = BRmin2leg * 1.
    
g.BRmin1leg = BRmin1leg
g.BRmin2leg = BRmin2leg


g.PrOX_filename = PrOX_filename

SPLIT_ZWetc = []
#SPLIT_ZWetc = ['Z','W']  #Contains particle to split according to PDG BRs (Z,W,t,H+,...)
#SPLIT_ZWetc = ['Z']  #Contains particle to split according to PDG BRs (Z,W,t,H+,...)
#print "NB: SPLIT_ZWetc = ",str(SPLIT_ZWetc)



# Import all the mass/decay info from the old (bad) structure into the new (somewhat better) structure
importFullParticleInfo_all()


partMO.sort(key=lambda obj: obj.mass)

# ##### SKIP SCENARIOS WHERE N1 IS NOT LSP #####
testLSP()  # may exit


if(checkmll): checkmllEtc()

# #####################################################

g.casc_starttypesD = {}
g.casc_starttypesD['all'] = ['gl','dL','uL','dR','uR','sL','cL','sR','cR','b1','b2','t1','t2', 'N1','N2','N3','N4','C1','C2',  'eR','eL','mL','mR','T1','T2','ve','vm','vT'] # all possible sparticles
g.casc_starttypesD['kromo'] = ['gl','dL','uL','dR','uR','sL','cL','sR','cR','b1','b2','t1','t2']
g.casc_starttypesD['nonkromo'] = ['N1','N2','N3','N4','C1','C2',  'eR','eL','mL','mR','T1','T2','ve','vm','vT'] 
g.casc_starttypesD['test'] = ['t2']

g.casc_starttypes = []

if mode_casc_starttypes not in list(g.casc_starttypesD.keys()):
    #print "NB: Illegal mode_casc_starttypes: '%s'   Set to 'all'" %(mode_casc_starttypes)
    print(mode_casc_starttypes)
    zz1 = mode_casc_starttypes.split(',')
    for zz2 in zz1: g.casc_starttypes.append(zz2)
    mode_casc_starttypes = mode_casc_starttypes.replace(',','_')

    print("NB: manually set starttypes to: ", g.casc_starttypes)

else:     
    g.casc_starttypes = g.casc_starttypesD[mode_casc_starttypes]


#g.casc_starttypes = ['N4','N3','N2','N1', 'C2','C1','uL']
#g.casc_starttypes = ['gl','dL','uL','sL','cL','dR','uR','sR','cR' ,'b1','b2']
#g.casc_starttypes = ['gl','dL','uL','dR','uR','b1','b2']  # Apr17
#g.casc_starttypes = ['gl','dL','uL','dR','uR','b1','b2','t1','t2']  # Apr17

#g.casc_starttypes = ['gl','dL','uL','dR','uR'] # Set1
#g.casc_starttypes = ['gl','dL','uL','dR','uR','sL','cL','sR','cR'] # Set2
#g.casc_starttypes = ['gl','dL','uL','dR','uR','sL','cL','sR','cR','b1','b2','t1','t2'] # Set3
#g.casc_starttypes = ['gl','dL','uL','dR','uR','sL','cL','sR','cR','b1','b2','t1','t2','N2','C1'] # Set4
#g.casc_starttypes = ['gl','dL','uL','dR','uR','sL','cL','sR','cR','b1','b2','t1','t2','N1','N2','N3','N4','C1','C2'] # Set5

#g.casc_starttypes = ['gl','dL','uL','dR','uR','sL','cL','sR','cR','b1','b2','t1','t2','N1','N2','N3','N4','C1','C2','lL','lR','sl','SL','T1']
#g.casc_starttypes = ['gl','dL','uL','dR','uR','sL','cL','sR','cR','b1','b2','t1','t2', 'N1','N2','N3','N4','C1','C2', 'eL','eR','ta1','ta2','snue','snut','lL','lR','vL','T1','T2']
#g.casc_starttypes = ['gl','dL','uL','dR','uR','sL','cL','sR','cR','b1','b2','t1','t2', 'N1','N2','N3','N4','C1','C2',  'eR','eL','mL','mR','T1','T2','ve','vm','vT'] # all possible sparticles


# 

#
# COMMENTS / RECOMMENDATIONS
# <Table of coverage>
# set1: SU4: e-6:0.602 , e-4:0.589 , e-3:0.557 , e-2:0.315
# set2: SU4: e-6:0.665 , e-4:0.647 , e-3:0.591 , e-2:0.315 
# set3: SU4: e-6:0.858 , e-4:0.838 , e-3:0.775 , e-2:0.499
# set4: SU4: e-6:0.986 , e-4:0.965 , e-3:0.896 , e-2:0.596
# set5: SU4: e-6:0.989 , e-4:0.965 , e-3:0.896 , e-2:0.596
#
#  The numbers commented to the right in the lines above (and put in table just above) show the 'coverage' for the given BRminXleg
#  See that having it at 1e-4 nearly gives full coverage.
#  Moving it to 1e-3 starts to affect the coverage, but is still ok. The number of 2legs is typically halved or more (good for overview).
#  See that 1e-2 is too harsh
#  Recommendation: use 1e-3 and set5 (to get full overview, except sleptons)
#                  use 1e-3 and set3 (to get all hard-jet events)

CRUCIALS.append(('SPLIT_ZWetc', SPLIT_ZWetc))


# ==========================================================================
# CHECKING MODES FOR SQUARK AND SLEPTON MERGING
#   QuarkMode == 1: c,u,s and d are counted separately
#   QuarkMode == 0: u=(u,c), d=(d,s)
#   (QuarkMode == 2: q=(u,d,c,s)   ... used for derived, but not here)
#   (QuarkMode == 3: q=(u,d,c,s,b) ... used for derived, but not here)
#   So if c and s are in casc_starttypes, QuarkMode must be 1
#   Similarly for mu,e
#   --qsep from command line ensures QuarkMode=1
#   --lsep from command line ensures LeptonMode=1


if STUDY3L or TWOLEGANALYSIS or ONELEGANALYSIS:
    for pid in ['sL','cL','sR','cR']:
        if pid in g.casc_starttypes and QuarkMode != 1:
            print("%s in g.casc_starttypes (%s) is only valid if the '--qsep' option is set." %(pid,str(g.casc_starttypes)))
            print('Exiting')
            sys.exit()

    for pid in ['mL','mR','vm']: #correct?   No, is incorrect. Anyway, need to fix Leptons first
        if pid in g.casc_starttypes and LeptonMode != 1:
            print("%s in g.casc_starttypes (%s) is only valid if the '--lsep' option is set." %(pid,str(g.casc_starttypes)))
            print('Exiting')
            sys.exit()
        

# ==========================================================================

cascStore = CascadeStore(g.BRmin1leg)  #1leg?
#cascStoreDilep = []

# ========================================================================== TWOLEGANALYSIS

if STUDY3L:
    study3L()


# ========================================================================== TWOLEGANALYSIS

if TWOLEGANALYSIS:
    analysis2legs()

    # some variables for easy access when interactive: 
    xT = g.ox[0]
    x  = g.ox[1]
    xInput = g.ox[2]


# ========================================================================== ONELEGANALYSIS


if ONELEGANALYSIS: analysis1leg()




# ========================================================================== SOME INFO/WARNINGS
# ========================================================================== AT THE END
# ========================================================================== 

if 0 and SHOWCRUCIALS:
    print(50*"=")
    print('CRUCIAL PARAMETERS:')
    for i in range(len(CRUCIALS)):
        print('%20s  :  %s' %(CRUCIALS[i][0], CRUCIALS[i][1]))
    print(50*"=")


NBs = []
NBs.append('SPLIT_ZWetc')
NBs.append('2leg: make1legs: setting JOIN  = 1 causes clear coverage loss in e.g. su1 and su3, but not su4 ... NOT UNDERSTOOD (may be other effects as well)')
NBs.append('Source of errors: copying of list and dict with incorrect methods (i.e. without .extend/.copy)')
NBs.append('Remember: prospino xsecs for smuons are probably not included')
NBs.append('Remember: make a prospino-xsec check of 4 mixing-matrices LHA conventions available')
NBs.append("Bug somewhere in SMp. In 400_000_150_000_150 pres4_* has 4 2legs (1e-4) while pres6 has none. FIXED: was due to missing key 'QQ' in FillSMp etc.")

if 0 and len(NBs) > 0:
    print('\nBEAWARE (variables known to cause problems under certain circumstances) [...features/bugs]')
    for nb in NBs: print('   ***** ',nb)

if 0 and len(NOTIFY) > 0:
    print('\n========== NOTIFICATIONS ==========')
    for noti in NOTIFY: print('   ***** ',noti)


# ################################################################################
# ################################################################################
#
# IMPROVEMENTS TO BE MADE:
#  [x]  u,d,(b) -> q
#  [x]  Split Z into qq and ll and ..., W ...  

#  [x]  Joining of signatures
#        [/] dd and uu join/sum (and rewrite as qq)
#        [/] also rewrite ud as 'qq'
#        [/] Option for rewriting uL,dL as qL and uR,dR as qR
#        [/] Option for adding qL and qR into sq
#  [x]  Also option for including b into sq
#  [ ]  Dropping of signatures with extra leptons (?) [at least put in cascStoreTrilep]

#  [x]  Store the Dilep list of one scenario to disk

#  SECONDARY PROGRAM (COMBINING)
#  [x]  Join Dilep lists, for identical chains pick highest BR
#  [ ]  Classify the cascades in kinematically Unique Cascade Kinematics classes (UCK)
#       [ ] each cascade belong to a UCK and can be labelled a UCK realisation (UCKR)
#       [ ] List the UCKs and their respective UCKR (with BRs): that's a "result"

#  Then analyse:
#  [ ]  Do a scan over M1, M2, mu, (vary tanB and a few others as well)
#       while keep common m(q) and m(gl) slightly above m(N4/C2)
#     [ ] Look at systematics from varying m(q) and m(gl) somewhat 
#     [ ] Look at systematics from varying tanB, ... 
#     - Remember: we are doing 'existence proof' (BR>1e-4), so only care for max BR
#       for a given UCK(R)


#  Later projects/tasks
#
#  [ ]  Classify Dilepton according to neutrinos, top, extra leptons, (b's), ...
#  [ ]  Classify non-Dilepton cascades
#  [ ]  Find and classify all possibilities for endpoint usage (1-lepton, ...)
#  [ ]  Investigate effects of non-degenerate slepton masses (OSOF-subtraction, ...)
# 


############################### Sep20 - REMAINING ITEMS
# DEBUGGING
# [/] Keep filling text of the debug: to understand where coverage is reduced
# [ ] NextDaughter, interchange, .. (and several other places): copy of list/dict: done correctly? 
# [x] Why is N1 not in initial particles?
# [x] Is the xsec correctly covered? (including N1?)
# [/] Drop/test unnecessary SUSY-simplifications in study3L_make1legs ?
# [/] Test with full l,T,q,b,y
#
# FEATURES
# [/] Clean: stop file creation (some 15 files...)
# [ ] control output by some g.VB3L levels (
# [ ] sensible handles to toggle between l, ... , and l,T,q,b,y
#  
# [ ] Make small user guide
# [ ] Make tarball




