#!/usr/bin/env python
'''
Program : signalgrid_loop
Version : 1.0, Dec 2012
Author  : b.k.gjelsten@fys.uio.no

Program description:
 - Loops over M1, M2, mu, tanbeta
 - Tunes higgs mass by changing At, typically done in 4-5 iterations (can probably learn to hit on first)
 - Time: 2s per go (can we do faster?) 
     - with 4-5 Higgsiterations we have ~10s per scenario

---------------
SUSYHIT description [SUSY-HIT = SU(spect)-S(deca)Y-H(decay)-I(n)T(erface)
 - A. Input 1: the *suspect2_lha.in file which specifies EWSB MSSM, Higgs options, MSSM24(27)-input (a-la what isasusy takes)
 - B. This is then sent to SuSpect v.X which does RGE etc. and gets the masses: output: suspect2.out + slhaspectrum_in
 - C. The slhaspectru_in file is then fed to SDecay and HDecay which does decays and produces the final susyhit_slha.out

 - BRIEF: Edit MSSM24-inputpars=suspect2_lha.in (in std lha format)   <-- this is done automatically by signalgrid_loop.py
          Suspect(suspect2_lha.in)             => suspect2.out and slhaspectrum_in
          SDecay and HDecay (slhaspectrum_in)  => susyhit_slha.out
---------------
I am using SUSY-HIT version 1.3, from July 14, 2012 (still the newest as of June 3, 2013) [see it from tar file date]
This includes:
  suspect2.f  2.41 (Aug/Dec 2008)   [As of May 23, 2013 there is a version 2.43] 
  HDecay      3.4  (Nov/Dec 2008)   [As of May 23, 2013 there is a version 5.10] 
  SDecay      1.3b (Dec 2008)       [is the most recent, is stamped July 14, 2012 on homepage]

Homepages: 
SUSPECT  http://www.coulomb.univ-montp2.fr/perso/jean-loic.kneur/Suspect/index.html
HDECAY   http://people.web.psi.ch/spira/proglist.html
SDECAY   http://www.itp.kit.edu/~maggie/SDECAY/
---
SUSYHIT  http://www.itp.kit.edu/~maggie/SUSY-HIT/

------------------------------------------------


OLD feature-wishes [maybe obsolete, or no longer considered very relevant]: 
 o validate
 / Higgs mass
   o add option for table-look-up
 x Isasusy part
   x DarkSUSY 
 o Other constraints
 o DarkSUSY on suspect-lha as well
 o internal parallelisation ... [tja ... do 
 o Option to let At be remembered from previous point .. should then often hit at first go (faster)
 o Verbose:
   o at start output the plans (how many loops + higgs intent + constraints ++)



NEW/STILL: 
 o time estimate:
   o initial dummy loop to check how many exist already
   o adaptive time estimate as we proceed
   o check which of the things take time
   
 o darksusy w/coann: option for breaking off?  <--- or do as bg job? [brute] 

 o memory: increases, is at 2 GB when a loop is at 1732/3332 .. so untenable ... is maybe fixed now



FRAGILE++

x s.require (tuning of h, T1, ...) is done manually in code



TODO (2014-01-21)

o small program to pick out:
  x warnings/errors ("caution") from *_slhaspectrum.in [part of BLOCK SPINFO] [3:warnings(spectrum may still be ok), 4:errors(results should not be used)]
  x fine-tuning parameters  _slhaspectrum.in:BLOCK SU_FINETUNE
  x constrained by exp data  _slhaspectrum.in:BLOCK SU_LOWPAR
o LARGE: Implement Softsusy++ as well


WEAKNESS
o new parameters from slha, e.g. the newly implemented 't' do not work out of the box.
  Need to change: susyhit_tool:susyhit_make_lha_in and pMSSM (need to specifically pick up, do not make general search [should I?])


2014-09: 

ISSUES
o add a 'try' (or other tests) somewhere on micromegas reading-back ; sometimes it fails (found out why, issue warnings)
  ex: fails with -mssmpars heavy=800  -mssmloop M3=2700:Q1=2500:Q2=2500:Q3=800:dR=2500:uR=2500:sR=2500:cR=2500:bR=800:tR=800

o introduced hack s.opt = opt for the prevar/preval stuff ; maybe standardise better the last one

o add some more standard tables (for easier easy-study)


'''

import sys, os, random
from math import sqrt  #,copysign
import traceback
import bkgjelstenArgReader 

import susyhit_tool
import softsusy_tool
import isasusy_tool


from kilelib import IterateNdim, SaveHistory, WriteToFile, tlt2dict
from libSUSY import txt2pdgId_plain
from lineup import lineup

#import pyslha  # hmm



# ##################################################### GLOBAL METHODS
# ##################################################### GLOBAL METHODS
# ##################################################### GLOBAL METHODS


##########
slhacode2txt = {
    'EXTPAR': {
    1:'M1', 
    2:'M2',
    3:'M3',
    11:'At',
    12:'Ab',
    13:'AT',
    23:'mu',
    25:'TB',  # Note, there is also TBatZ in MINPAR below ; that one will be overridden by the current, if present
    26:'mA',
    31:'L1',
    32:'L2',
    33:'L3',
    34:'eR',
    35:'mR',
    36:'TR',
    41:'Q1',
    42:'Q2',
    43:'Q3',
    44:'uR',
    45:'cR',
    46:'tR',
    47:'dR',
    48:'sR',
    49:'bR'
    },
    
    'MINPAR': {
    3:'TBatZ'  # Note TB above (TB will beat TBatZ, as in the SLHA convention)
    },

    'SMINPUTS': {
    6:'t'
    }

    }

##########
def GetInputParsFromSlha(fn_slha):
    par = {}
    
    # Open slha
    slhaObj = pyslha.readSLHAFile(fn_slha)
    

    # Read
    for block in list(slhacode2txt.keys()):
        for slhacode in slhacode2txt[block]:
            var = slhacode2txt[block][slhacode]
            if slhacode not in slhaObj[0][block].entries: continue  # if code not present, e.g. 25 (TB at SUSYEW scale)
            val = slhaObj[0][block].entries[slhacode]
            par[var] = val

    #for var in par: print 'GetInputParsFromSlha   %-4s  %10.3f' %(var, par[var])  # verbosity

    return par

    
##########
class pMSSM:
    # This class can be extended to to various things
    # For now it is only a placekeeper for the parameter dictionary (.par)

    def __init__(s, opt={}):
        #print 'opt: ', opt
        s.par = {}
        #s.heavy = opt.get('heavy',3000.)  # When called by signalgrid opt now always contains 'heavy'  # 2013-12-21
        s.par['heavy'] = opt.get('heavy',3000.)
        
        #print 'DEBUG init: heavy = ', s.par['heavy']
        for z in ['EWSB_scale']:
            if z in opt: s.par[z] = opt[z]  
        s.allowed_nonpMSSMkeys = ['A','mslep','mq','sleprelmpar','Afunnel']  # should probably get this list from the caller
        s.SetPrimaries(opt=opt)
        s.SetSecondaries(opt=opt)
        s.CheckPars(pars=opt)
        

    # ##### 
    def CheckPars(s, pars, VB=1):
        outs_status = []
        for key in list(pars.keys()):
            if key not in s.par and key not in s.allowed_nonpMSSMkeys: 
                outs_status.append('WARNING  pMSSM::CheckPars  Found unknown parameter in par dict:  %-10s   %s' %(str(key), str(pars[key])))
                if VB>0: print(outs_status[-1])

        return outs_status
    
    
    # ##### 
    def UpdateAndGetPars(s, opt={}):
        outs_updatestatus = s.UpdatePars(opt=opt)
        return s.par, outs_updatestatus
    

    # ##### 
    def UpdatePars(s, opt={}):
        outs_updatestatus = s.CheckPars(pars=opt)  # 2014-01-30
        #s.SetPrimaries(opt=opt)
        for o in opt: s.par[o] = opt[o]
        s.SetSecondaries(opt=opt)
        return outs_updatestatus
    
    
    # #####
    def SetPrimaries(s, opt={}): 
        par = s.par
        heavy = s.par['heavy']

        # M3 (~gluino mass): 1 par
        par['M3'] = opt.get('M3',heavy)

        # squark masses: 9 pars
        msq = opt.get('msq', heavy)
        par['Q1'] = opt.get('Q1', msq)
        par['Q2'] = opt.get('Q2', msq)
        par['Q3'] = opt.get('Q3', msq)
        par['dR'] = opt.get('dR', msq)
        par['uR'] = opt.get('uR', msq)
        par['sR'] = opt.get('sR', msq)
        par['cR'] = opt.get('cR', msq)
        par['bR'] = opt.get('bR', msq)
        par['tR'] = opt.get('tR', msq)

        # slepton masses: 6 pars
        mslep = opt.get('mslep', heavy)
        par['L1'] = opt.get('L1', mslep)
        par['L2'] = opt.get('L2', mslep)
        par['L3'] = opt.get('L3', mslep)
        par['eR'] = opt.get('eR', mslep)
        par['mR'] = opt.get('mR', mslep)
        par['TR'] = opt.get('TR', mslep)

        
        # neutralino/chargino sector: 3 pars
        #par['M1'] = opt.get('M1', 100.)
        #par['M2'] = opt.get('M2', 200.)
        #par['mu'] = opt.get('mu', 500.)
        par['M1'] = opt.get('M1', heavy)
        par['M2'] = opt.get('M2', heavy)
        par['mu'] = opt.get('mu', heavy)

        # Higgs sector: 2 pars
        #par['TB'] = opt.get('TB', 10)  # tan(beta)
        if 'TB' in opt: par['TB'] = opt['TB']   # 2014-01-21: replaced one line above with these 3 lines: to allow TB at Z scale too
        elif 'TBatZ' in opt: par['TBatZ'] = opt['TBatZ']
        else: par['TB'] = 10
        #par['mA'] = opt.get('mA', 2000.)  # For TB20_140_300_300 need this at 1400 to avoid (intermediate) tachyon A (odd-higgs) 
        par['mA'] = opt.get('mA', heavy)  # For TB20_140_300_300 need this at 1400 to avoid (intermediate) tachyon A (odd-higgs) 


        # trilinear couplings: 3 pars
        #A = opt.get('A', 3400.)
        A = opt.get('A', heavy)
        par['At'] = opt.get('At', A)
        par['Ab'] = opt.get('Ab', A)
        par['AT'] = opt.get('AT', A)


        # other
        if 't' in opt: par['t'] = opt['t']
        
        # this is s.par (pMSSM.par) being set


    # #####
    def SetSecondaries(s, opt={}):  # make optional? 
        par = s.par
        
        # (trivial) additionals needed e.g. by suspect
        if 'Au' not in par: par['Au'] = par['At']
        if 'Ad' not in par: par['Ad'] = par['Ab']
        if 'Ae' not in par: par['Ae'] = par['AT']

        # Secondaries (mixing, slepton mass, ..)
        if opt.get('mixstop_mode','') == 'max':
            par['At'] = sqrt(6. * par['Q3'] * par['tR']) # + par['mu'] / par['TB']
        if 'mixstop' in opt:
            par['At'] = sqrt(opt['mixstop'] * par['Q3'] * par['tR'])# + par['mu'] / par['TB']
        if 'mixstopreltomax' in opt:
            par['At'] = opt['mixstopreltomax'] * sqrt(6. * par['Q3'] * par['tR'])   # + par['mu'] / par['TB']

        # Note it is not At, but Xt = At-mu/TB which is the crucial variable
        # So At = Xt + mu/TB ;  Xt = At - mu/TB
        # max mixing is for Xt=sqrt(6*t1*t2)
        # So, scaling is maybe best by: Xt = (At-mu/TB) = XtMax * sc = sqrt(6*t1*t2) * sc
        #  => At = mu/TB + sc*sqrt(6*t1*t2)



##########
class signalgrid: 

    def __init__(s, cmd=[]):
        s.myname = sys.argv[0].split('/').pop()
        s.VB = 1  # D=2
        #maiken changed
        #s.HOME = os.getenv('HOME')
        s.HOME = os.getenv('PWD')
        #s.gridname = 'agrid'  # not yet in use

        #s.fnstart = 'susy_DGnoL'
        #s.LoopVarsMode = 'DGnoL_v0'

        s.dict = {}

        s.warn = []
        s.fn_warn = 'warnings_signalgrid_loop.txt'

        s.fnstart = 'susy'  # DEFAULT
        s.fnadd = 'MYGRID' # <--- TO BE SETto be set on command line
        s.LoopVarsMode = 'free'  # DEFAULT

        
        s.LoopVars = {}
        #s.pMSSM_mode = 'DGnoL_v0'
        #s.pMSSM_opt = {}

        s.calculators = ['susyhit','isasusy','softsusy','suspect']
        #s.calculators_dotoo = ['susyhit','isasusy','softsusy']  # well, switch off by default? 
        #s.calculators_dotoo = ['susyhit','isasusy', 'softsusy']
        s.calculators_dotoo = []
        s.calculator = 'susyhit'


        #s.pMSSM_object = pMSSM()  # so far no options, then all is set to heavy(3TeV) except M1,M2,mu,TB,mA
        #s.heavy = 3000.  # GeV  # will be put in s.pars_input_manual

        s.pmssm_vars = ['TB','mu','M1','M2','mA', 'M3', 'L1','L2','L3','eR','mR','TR', 'Q1','Q2','Q3','dR','uR','sR','cR','bR','tR', 'At','Ab','AT']

        s.pars_input_manual = {}
        s.pars_input_manual['heavy'] = 3000.
        s.showpmssmvals = 0
        s.pars_pmssmdefault = {'M1':100., 'M2':200., 'mu':500., 'mA':2000., 'A':3400}  # Note: TB is not in ... (complications with TB, TBatZ)
        #s.pars_pmssmdefault = {}
        s.pars_loop = {}
        s.pars_looporder = []
        s.pars_loop_multival = []   # this contains only the vars actually involved in a grid, i.e. the ones which have multiple input values
        
        # should allow this to be input from command line ; or to select predefined
        s.require = {}
        s.setequalto = {}

        s.runloop = 0  # main parameter to decide if loop (Main) is run or not
        
        s.runcard = ""  # if mssm loop parameters set through a runcard
        

        ##s.require['h'] = {'nom':125.50, 'min':125.45, 'max':125.55, 'tunevar':'At', 'step':300}
        #s.require['h'] = {'nom':125.50, 'delminmax':0.02, 'tunevar':'At', 'step':300}  # DGnoSL
        ##s.require['T1'] = {'nom':{'abs_N1':1-0.5, 'abs_N2':0.5}, 'delmin':0.01, 'delmax':0.01, 'tunevar':'TR', 'step':5}
        s.require_mode = []

        s.configsetups = ['DGnoSL','DGstauR50','DGemtR50','DGemtR50noh','Higgs0.3','Higgs','hfunnelM1','AfunnelUpDS','AfunnelDownDS']
        s.configsetup = []  # Newest addition. Can control several things. For now only s.require

        if not 'SUSYPHENO_PATH' in list(os.environ.keys()):
            print("$SUSYPHENO_PATH not set. Exiting.")
            sys.exit()
        susyphenopath = os.environ['SUSYPHENO_PATH']

        s.darksusy_exe = '%s/bin/darksusy.py' %(susyphenopath)
        s.fn_ds_script = 'script_darksusy_do.sh'
        s.ds_script = []
        s.ds_opt_cdmcoann = 1
        s.ds_donow = 1

        s.force = 0

        # Number of digits in filename
        s.ndig = {}
        s.ndig['default'] = 3
        s.ndig['M1'] = 3
        s.ndig['M2'] = 3
        s.ndig['mu'] = 3
        s.ndig['mA'] = 4
        s.ndig['TB'] = 1
        s.ndig['eR'] = 3
        s.ndig['mR'] = 3
        s.ndig['TR'] = 3
        s.ndig['L1'] = 3
        s.ndig['L2'] = 3
        s.ndig['L3'] = 3
        s.ndig['sleprelmpar'] = 4

        s.format = {}
        s.format['sleprelmpar'] = '%4.2f'

        #s.input = {}
        #s.input['M1'] = []
        #s.input['M2'] = []
        #s.input['mu'] = []
        #s.input['TB'] = []

        s.nBreak = 9999999

        s.M1welltempered = 0      # 0:off, 1:mode1, ..[not implemented more, but can be]
        s.M1welltempered_rel = -1 # in use (and set anew) when above is on 

        s.mode = {}
        #s.mode['fnID'] = 1   # mode of setting filename: 1,2  (See Get_fnID())
        s.mode['fnID'] = -1   # 2014-01-23:  -1:fully general with _ as delimiter, -2:ditto with = as delimiter (tab no good)

        s.useEifabove = 1e9  # 1e4

        s.susyhit_deltmp = 1
        
        s.dict['EWSB_scale'] = -1  # 2014-03-09

        s.dict['slha2wig_py'] = "~/bin/slha2wig.py"
        s.dict['isaplay_py'] = "~/bin/isaPlay.py"


        s.dict['makewig'] = 0
        s.dict['makeder'] = 0
        s.dict['makedermass'] = 0


        s.dict['makeoneliner_suspect'] = 0  # 2014-03-11: default is now off, is mostly OBSOLETE
        s.dict['maketwoliner_scanvars'] = 1
        s.dict['errlowtune_addscanvars'] = 0  # 2014-01-30  ... turn this off by default due to risk of making tlts which cannot be pasted

        s.dict['slha2tlt'] = 2      # 0:none, 1:main calculator, 2:all
        s.dict['makedarksusy'] = 1  # 0:none, 1:main calculator, 2:all
        s.dict['micromegas'] = 1    # 0:none, 1:main calculator, 2:all

        s.dict['slha2tlt_tables'] = 'masses,,extpar,,3genmix,C1mix,Nmix1:::mixing,,elweak,,rest'
        s.dict['slha2tlt_opt'] =  '-dict I,forceremakeflat,1'


        
        s.dict['isaplay_maxmass'] = -1
        s.dict['isaplay_maxmass_relofheavy'] = 0.8
        s.dict['isaplay_brmin'] = 0.0001

        s.dict['slepmass_maxofN2'] = 1.  # corr. to DGemtR50
        s.dict['slepmass_rel'] = 0.5     # corr. to DGemtR50

        s.dict['sleprelmpar_which'] = ['eR','mR','eL','mR','TR','TL']
        s.dict['sleprelmpar_maxrelto'] = 'N2'  # use this mass as max
        s.dict['sleprelmpar_minrelto'] = 'N1'  # use this mass as max

        s.dict['Afunnel_Amin']  = 0.
        s.dict['Afunnel_scale'] = 1.0
        s.dict['Afunnel_add']   = 0.

        s.dict['AfunnelUp_nom']   = 1.
        s.dict['AfunnelUp_delminmax']   = 0.02

        s.dict['AfunnelDown_nom']   = 1.
        s.dict['AfunnelDown_delminmax']   = 0.02

        s.dict['h_delminmax'] = 0.02   # this is default as used for e.g. DGnoSL

        s.dict['derivedAfterGo1_keys'] = ['sleprelmpar', 'Afunnel']
        
        s.dict['hfunnel_add'] = 0.
        s.dict['hfunnel_delminmax'] = 0.1

        s.dict['require_ntriesMax'] = 20   # 10

        #s.dict['stdvars'] = ['TB','M1','M2','mu','mA']   # these will be put in all tables ; facilitates plotting
        s.dict['stdvars'] = ['TB','M1','M2','mu','mA','eR','mR','TR']   # these will be put in all tables ; facilitates plotting
        s.dict['stdvars_format'] = '%.0f'
        #s.stdvals = []

        s.par_steer = {}   # steers the calculators (number of loops etc.) 

        #s.dict['sleprelmpar_max'] =   # use this mass as max

        # mass is given as 
        

        s.derivedAfterGo1 = {}
        s.fn_slhainput = ''
        s.pars_input_fromslha = {}

        s.fn_history = 'history_signalgrid_loop.txt'
        s.fn_history_global = '%s/.history_signalgrid_loop.txt' %(s.HOME)

        s.scanvars_def = 'all'
        s.scanvars_defs = ['all','multi']  # will add 'manual' or 'user'

        s.dict['stopifunknownparameter'] = 1  # 2014-01-30

        
        s.fn_history = 'history_signalgrid_loop.txt'
        s.fn_history_global = '%s/.history_signalgrid_loop.txt' %(s.HOME)
        SaveHistory(fn=s.fn_history, argv=sys.argv, opt=['stripcmd','date'])
        SaveHistory(fn=s.fn_history_global, argv=sys.argv, opt=['stripcmd','date','dirsameline'])

        # -------------
        # -------------
        if 'ReadArg' in cmd: s.ReadArg()


        
        # -------------
        # -------------
        # AfterMath
        s.PostInit()  # Sets up loop, 
        

        if 'Main' in cmd:
            if s.runloop: s.MainLoop()
            else: print('WARNING  MainLoop() not run. Need to specify with -mssmloop and/or -mssmpars.')
        
    
    # ##########
    def showHelp(s):
        print()
        print('DESCRIPTION')
        print(' This python script takes as input pMSSM parameters (and scan type parameters):')
        print('   then runs SUSY calculator(s) [susyhit, softsusy, isasusy] to produce slha file(s)')
        print('     options: autotune Higgs mass, autotune slepton masses, ...')
        print('   then runs secondary programs on the slhas: DarkSUSY, MicrOmegas')
        print('   output: tlt (two-line-table) which has the quantity/value on the first/second line')
        print()
        print(' The tlt-output is a practical base for making tables with one line per scenario (and the quantity on the first line.)')
        print(' Such tables are then appropriate format for plotting with the ROOT-plotter script munch.py')
        print()

        print('\nOUTPUT FILES (tlt)')
        print("   <pointdesc>_scanvars_all.tlt        : just the variables specified on input line with '-mssloop <list of vars>'  (i.e. the ones not taking default values)")
        print('   <pointdesc>_scanvars_multi.tlt      : similar, but only vars which on input line has more than one value')
        print('   <pointdesc>_susyhit_slha_mass.tlt   : input vars + all SUSY masses')
        print('   <pointdesc>_susyhit_slha_extpar.tlt : input vars + ~all Lagrangian parameters')
        print('   <pointdesc>_susyhit_slha_elweak.tlt : input vars + parameters relevant for elweak sector')
        print('   <pointdesc>_susyhit_slha_rest.tlt   : input vars + some additional vars (the rest)')
        print('   <pointdesc>_susyhit_slha_mixing.tlt : input vars + mixing parameters in squark/slepton/neutralino/chargino sector')
        print('   <pointdesc>_susyhit_errlowtune.tlt  : input vars + errors + finetuning parameters + experimentally constrained vars ')
        print('   <pointdesc>_susyhit_MO.tlt          : input vars + MicrOmegas output (DarkMatter++, experimentally constrained vars)')
        print('   <pointdesc>_susyhit_ds.tlt          : input vars + DarkSUSY output (DarkMatter, experimental bounds)')
        
        print() 

        print('\nUSAGE: %s <options>' %(s.myname))
        print()
        print('OPTIONS: ')
        print("   --showpmssmpars  : shows all 24 pmssm parameters [Note: one letter per sparticle; T=tau, m=muon, .. ; TB=tan(beta)]")
        print("   --showpmssmvals  : shows default values of all 24 pmssm parameters (and maybe a few scan variables)")
        print("   -mssmloop <Colon-separated list of pmssm vars for loop>")
        print("      Ex: -mssmloop M1=100,200:M2=300:mu=400,500   : gives 2D scan in M1 x mu  (100,200)x(400,500)")
        print("   -mssmpars <Colon-separated list of pmssm/scan vars>  ; here only one value per variable (for scan variables use -mssmloop)")
        print("      can be specific parameters like mA,eR,M1")
        print("      can also be pseudopars like heavy, msq, mslep (the value of 'heavy', of squarks, of sleptons")
        print("      Ex: -mssmpars heavy=3500,msq=3000,mslep=2000,bR=1500,eR=500 : resets heavy, msq, mslep, but mb1 and meR get their specific values")
        print("   Squark/slepton values: if specific value is set, take it; elif msq/mslep is set, take it; else take heavy")
        print("   For A-coupling, if specific value is set (At,Ab,AT), take it, else take A")
        print() 
        print("   --maxstopmix                : sets pmssm par At to value which gives maximum stop mixing (and ~maximum mass of lightest Higgs)")
        print("   -mixstopreltomax <mixvalue> : amount of mixing: 1 gives max mixing, 0 gives no mixing, intermidiate gives intermediate  (affects At)")
        print("   --cdmcoann     : include coannihilation procs in DarkMatter calculation (is more accurate, but can in some cases take noticeably longer)")
        print("   -configsetup <list of setups> : specify predefined configsetup.")
        print("                Ex: '-configsetup Higgs' involves a tuning of Higgsmass to 125.5 GeV within h_delminmax=%.3f GeV" %(s.dict['h_delminmax']))
        print("                Existing configsetups: %s" %(s.configsetups))
        print("     -dict F,h_delminmax,0.3    # to change h_delminmax to 0.3  (Note: -dict can be used to change other vars as well. Will then need to be separated by ':'")

        print() 
        print("   -calc <spectrum calculator>       : specify spectrum calculater. D=%s ; Allowed ones:%s" %(s.calculator, s.calculators))
        print("   -calc_dotoo <spectrum calculator> : add more calculators to use alongside the main")
        print("   --allspectrumcalc                 : add all to dotoo")
        print("   --force : will recalculate points already calculated (default is to not recalculate) ")
        print() 
        print("   The output filenames are based on fnstart and fnadd and fnIDmode ")
        print("     -fnIDmode <mode:-3/-2/-1/0/1> : D=-1 (uses '_' to separate var and val; with -2/-3 the separator is '='/'eq')   (see code for more details)")
        print("     -fnstart <txt> : common beginning of output filename base (D='susy')")
        print("     -fnadd   <txt> : common end of output filename base (D='MYGRID')")
        print("     -ndig <var1:ndig1,var2:ndig2,..>  : filenames: reset number of digits used for select vars.  (Ex: '-ndig M1:4,mA:5') ")
        print("     -useEifabove <value>              : filenames: use scientific (E) format if values above specified value")
        print() 
        
        
        #print '    Ex: %s  -loopvarsmode free  -fnadd DGnoSL  -TB 10  -M1M2mu 100,120,200  -M1add 45,65,85' %(s.myname)
        #print '    Ex: %s  -loopvarsmode free  -fnadd DGnoSL  -TB 10  -M1welltempered 1,0.5  -M2mu 100,120,150,200,250,300,350,400,450,500,999' %(s.myname)
        #print '    Ex: signalgrid_loop.py  -useEifabove 1e3  -fnIDmode 1  -fnstart susy  -fnadd DGnoSL  -mssmpars mA=1000   -TB 10  -M1 50   -M2mu 100,120,140,150,160,180,200,250,300,350,400,450,500,3e3      # 2013-09-16 DGnoSL official try1'
        #print '    Ex: signalgrid_loop.py  -useEifabove 1e3  -fnIDmode 1  -fnstart susy  -fnadd DGnoSL  -mssmpars mA=1000   -TB 10  -M1 3e3  -M2mu 100,120,140,150,160,180,200,250,300,350,400,450,500,3e3      # 2013-09-16'
        #print '    Ex: signalgrid_loop.py  -useEifabove 1e3  -fnIDmode 1  -fnstart susy  -fnadd DGnoSL  -mssmpars mA=1000   -TB 10  -M1 55  -M2mu 100,120,140,150,160,180,200,250,300,350,400,450,500,3e3'
        #print '    Ex: signalgrid_loop.py  -useEifabove 1e3  -fnIDmode 1  -fnstart susy  -fnadd DGstauR50  -mssmpars mA=2000   -TB 50  -M1 75  -M2mu 200  -configsetup DGstauR50  # 2013-11-22  Does not reproduce exactly. Probably due to slight change in requirements on T1 and h (and/or their combination)'
        print()
        #print '   -mssmparstxt mixstop_mode:max   # to force max stop mixing, and thus roughly max higgs mass'
        #print '   -mssmpars mixstop=6   # to force max stop mixing, and thus roughly max higgs mass'
        #print '   -mssmpars mixstopreltomax=1.   # 1.:to force max stop mixing, 0.5 to have At=AtMax/2 etc. '
        #print '   --maxstopmix   # max stop mixing (yet another way)'
        
        #print '    (NB: seems we cannot have two parallel loops in the same directory (should check))'

        print("   -vb <0/1/...>  : verbosity level; 0 is least verbosity (D=1)")
        print("   -h                : this help message")
        print("   --h               : brute info on implemented options")
        print("   --dict            : dump various default values (which can be changed)")
        print("  (See code for more options/details.)")
        print() 
        print("  Note: For reference previously issued commands are kept in:")
        print("     %s  (this directory)" %(s.fn_history))
        print("     %s  (any directory)" %(s.fn_history_global))
        print()
        print("EXAMPLES: ")
        print('  signalgrid_loop.py  -fnadd _testgrid0  -mssmloop M1=100:M2=200:mu=400:TB=10:mA=2000   # one pmssm point, the non-specified input variables take default values (try --showpmssmvals')
        print('  signalgrid_loop.py  -fnadd _testgrid1  -mssmloop M1=100:M2=100,150,200,250,300,350,400:mu=100,150,200,250,300,350,400:TB=10:mA=2000  # scan: 2D grid in M2xmu, 7x7 points')

        
    # ##########
    def PostInit(s):

        if s.calculator in s.calculators_dotoo: s.calculators_dotoo.remove(s.calculator)

        #if 'heavy' in s.pars_input_manual: sys.exit('heavy already')
        #s.pars_input_manual['heavy'] = s.heavy  # 2013-12-21

        if s.fn_slhainput: 
            s.pars_input_fromslha = GetInputParsFromSlha(s.fn_slhainput)

        if 'TB' in s.pars_input_fromslha and 'TBatZ' in s.pars_input_manual: pars_input_fromslha.pop('TB')  # needed since TB overrides TBatZ
            
        s.pars_input = dict(s.pars_input_fromslha) # could be empty (if no file given)
        s.pars_input.update(s.pars_input_manual)   # manual will override any common values 
        if s.dict['EWSB_scale'] == -1 and 'EWSB_scale' not in s.pars_input:
            s.pars_input['EWSB_scale'] = -1   # hack required for softsusy [susyhit already defaults to -1]

        if s.VB>1:
            print('s.pars_loop: ', s.pars_loop)
            print('s.pars_input_manual: ', s.pars_input_manual)
            print('s.pars_input: ', s.pars_input)
            print('s.pars_loop:  ', s.pars_loop)
            print('s.pars_pmssmdefault: ', s.pars_pmssmdefault)
        
        # Moved default values out of pMSSM
        for var in s.pars_pmssmdefault:
            if var not in s.pars_input:
                s.pars_input[var] = s.pars_pmssmdefault[var]
                if s.VB>1: print('INFO  Setting default for %-5s: %10.3f' %(var, s.pars_pmssmdefault[var]))
            
        

        if s.VB>2:
            for var in s.pars_input:
                zloop = ''
                if var in s.pars_loop: zloop = '   Will be overriden by loop values: %s' %(s.pars_loop[var])

                try:    print('PostInit :  s.pars_input  %-7s  %10.3f %s' %(var, s.pars_input[var], zloop))
                except: print('PostInit :  s.pars_input  %-7s  %s %s' %(var, s.pars_input[var], zloop))
        
        s.pMSSM_object = pMSSM(opt=s.pars_input) 
        #if s.pars_input: s.pMSSM_object.UpdatePars(opt=s.pars_input)

        if s.VB>2 or s.showpmssmvals:
            for var in sorted(s.pMSSM_object.par):
                zloop = ''
                if var in s.pars_loop: zloop = '   Will be overriden by loop values: %s' %(s.pars_loop[var])
                print('PostInit:  s.pMSSM_object.par  %-7s  %10.3f %s' %(var, s.pMSSM_object.par[var], zloop))
        
        # The above moved hereto on 2014-01-21

        if 'DGnoSL' in s.configsetup:
            s.require['h'] = {'nom':125.50, 'delminmax':s.dict['h_delminmax'], 'tunevar':'At', 'step':300}  # DGnoSL
            # s.require['h'] = {'nom':125.50, 'min':125.45, 'max':125.55, 'tunevar':'At', 'step':300}

            
        if 'DGstauR50' in s.configsetup:  # This does not reproduce exactly, but rounded
            s.require['T1'] = {'nom':{'abs_N1':1-0.5, 'abs_N2':0.5}, 'delmin':0.01, 'delmax':0.01, 'tunevar':'TR', 'step':5}
            s.require['h'] = {'nom':125.50, 'delminmax':s.dict['h_delminmax'], 'tunevar':'At', 'step':300}  # DGnoSL

        if 'DGemtR50' in s.configsetup:
            #s.require['T1'] = {'nom':{'abs_N1':1-0.5, 'abs_N2':0.5}, 'delmin':0.01, 'delmax':0.01, 'tunevar':'TR', 'step':5}
            s.require['eR'] = {'nom':{'abs_N1':1-0.5, 'abs_N2':0.5}, 'delmin':0.01, 'delmax':0.01, 'tunevar':'eR', 'step':5}
            s.require['h'] = {'nom':125.50, 'delminmax':s.dict['h_delminmax'], 'tunevar':'At', 'step':300}  # DGnoSL
            s.setequalto['mR'] = 'eR'
            s.setequalto['TR'] = 'eR'  # might not work due to yukawa difference [think it does, but need to recheck]
        
        if 'DGemtR50noh' in s.configsetup:
            #s.require['T1'] = {'nom':{'abs_N1':1-0.5, 'abs_N2':0.5}, 'delmin':0.01, 'delmax':0.01, 'tunevar':'TR', 'step':5}
            s.require['eR'] = {'nom':{'abs_N1':1-0.5, 'abs_N2':0.5}, 'delmin':0.01, 'delmax':0.01, 'tunevar':'eR', 'step':5}
            #s.require['h'] = {'nom':125.50, 'delminmax':s.dict['h_delminmax'], 'tunevar':'At', 'step':300}  # DGnoSL
            s.setequalto['mR'] = 'eR'
            s.setequalto['TR'] = 'eR'  # might not work due to yukawa difference [think it does, but need to recheck]
        

        if 'Higgs' in s.configsetup:
            s.require['h'] = {'nom':125.50, 'delminmax':s.dict['h_delminmax'], 'tunevar':'At', 'step':300}

        if 'Higgs0.3' in s.configsetup:
            s.require['h'] = {'nom':125.50, 'delminmax':0.3, 'tunevar':'At', 'step':300}
            
        #if 'A' in s.configsetup:
        #    s.require['A'] = 

            #s.require['eR'] = {'nom':{'abs_N1':s.dict['slepmass_rel'], 'abs_N2':s.dict['slepmass_maxofN2']*s.dict['slepmass_rel']}, 'delmin':0.01, 'delmax':0.01, 'tunevar':'eR', 'step':5}
            # ex.: slepmass_maxofN2 = 2  means sleptons will have mass between N1 and 2*N2, decided by slepmass_rel 
            #s.setequalto['mR'] = 'eR'
            #s.require['eL'] = {'nom':{'abs_N1':s.dict['slepmass_rel'], 'abs_N2':s.dict['slepmass_maxofN2']*s.dict['slepmass_rel']}, 'delmin':0.01, 'delmax':0.01, 'tunevar':'eL', 'step':5}
            #s.setequalto['mL'] = 'eL'


        if 'hfunnelM1' in s.configsetup:
            s.require['N1'] = {'nom': 125.5/2 + s.dict['hfunnel_add'], 'delminmax':s.dict['hfunnel_delminmax'], 'tunevar':'M1', 'step':2.}


        if 'AfunnelUpDS' in s.configsetup:
            s.require['DMwCoVSexp'] = {'nom': s.dict['AfunnelUp_nom'], 'delminmax':s.dict['AfunnelUp_delminmax'], 'tunevar':'mA', 'step':10., 'maxstepAbs':15}

        if 'AfunnelDownDS' in s.configsetup:
            s.require['DMwCoVSexp'] = {'nom': s.dict['AfunnelDown_nom'], 'delminmax':s.dict['AfunnelDown_delminmax'], 'tunevar':'mA', 'step':-10., 'maxstepAbs':15}


        if s.VB: s.ShowInfo(mode=0)



    # ##########
    def ShowInfo(s,mode=0):
        vars = ['M1','M2','mu','TB']
        vars = s.pars_looporder
        nPoints = 1
        outs = []
        for var in vars: 
            #outs.append('ShowInfo::  %-5s  %3i  %s' %(var, len(s.LoopVars[var]), s.LoopVars[var]))
            #nPoints *= len(s.LoopVars[var])
            outs.append('ShowInfo::  %-5s  %3i  %s' %(var, len(s.pars_loop[var]), s.pars_loop[var]))
            nPoints *= len(s.pars_loop[var])
        secsPerPoint = 2.
        totsecs = int(secsPerPoint*nPoints) 
        hours = totsecs/3600 ; mins = (totsecs-hours*3600)/60 ; secs = (totsecs-hours*3600-mins*60)
        thetime = '%ih%2sm%2ss' %(hours, str(mins).zfill(2), str(secs).zfill(2))
        outs.append('ShowInfo::  Number of points: %i    Time: %s  (%is) assuming %.2fs per point' %(nPoints, thetime, totsecs, secsPerPoint))

        for out in outs: print(out)
        



    # ##########
    def GetValueFromReq(s, reqval, slhamass):
        if type(reqval) in [int,float]:
            val = reqval
        else:
            # assume a dictionary, e.g. {'N1':0.25, 'N2':0.75} or (better) {'abs_N1':0.25, 'abs_N2':0.75}
            # when using 'abs_', abs values are taken, this can be needed e.g. for neutralinos which can have negative mass parameter
            val = 0.
            for zvar0 in reqval:
                zvar = zvar0.replace('abs_','')
                pdgid = txt2pdgId_plain[zvar]
                addval = slhamass[pdgid] * reqval[zvar0]
                if zvar0.startswith('abs_'): addval = abs(addval)
                val += addval

        return val
        
        
    # ##########
    def DerivedAfterGo1(s, slhamass, origSUSYpar):  # 2014-01-20  ... all this is rather adhoc
        mssmChanges = {}

        #print 'slhamass keys: ', slhamass.keys()

        # #####
        # NOTE: NEED TO DOCUMENT THIS A BIT MORE, GIVE AN EXAMPLE
        if 'sleprelmpar' in s.pars_loop: 
            for var in s.dict['sleprelmpar_which']:
                zrel = origSUSYpar['sleprelmpar'] #  s.pars_loop['sleprelmpar']
                zreltomin = txt2pdgId_plain[s.dict['sleprelmpar_minrelto']]
                zreltomax = txt2pdgId_plain[s.dict['sleprelmpar_maxrelto']]
                if s.VB>3: print('DerivedAfterGo1: ', zrel, zreltomin, zreltomax)
                if zrel < 1: 
                    val = abs(slhamass[zreltomin]) * (1-zrel) + abs(slhamass[zreltomax]) * zrel
                else:
                    val = abs(slhamass[zreltomax]) * zrel
                mssmChanges[var] = val

        # #####
        if 'Afunnel' in s.pars_loop:  # the par Afunnel itself is just on or off 
            mssmChanges['mA'] = max(s.dict['Afunnel_Amin'], 2 * abs(slhamass[txt2pdgId_plain['N1']]) * s.dict['Afunnel_scale'] + s.dict['Afunnel_add'])
            # default: Afunnel_Amin=0, Afunnel_scale=1.; Afunnel_add=0 " which gives exactly A = 2*N1, which should be optimal

            
        
        if s.VB>3: print('DerivedAfterGo1: %s' %(mssmChanges))
        return mssmChanges

    
    # ##########

    # ##########
    def TuneVariable(s, var, SUSYparCopy, slhaObj, require, firstIteration, verbose, opt, istuned):

        #var = 'h'
        #if var in s.require:

        req=require[var]

        slhamass = slhaObj[0]['MASS'].entries

        
        # ### INIT

        valNOM = s.GetValueFromReq(req['nom'], slhamass)  # (should work fine with AfunnelUpDS)
        #if type(req['nom']) in [int,float]: valNOM = req['nom']
        #else:
        #    valNOM = 0.
        #    for zvar in req['nom']:
        #        pdgid = txt2pdgId_plain[zvar]
        #        valNOM += slhamass[pdgid] * req['nom'][zvar]


        if 'min' in req: valMIN = req['min']
        elif 'delmin'    in req: valMIN = valNOM - req['delmin']
        elif 'delminmax' in req: valMIN = valNOM - req['delminmax']
        elif 'relmin'    in req: valMIN = valNOM * (1 - req['relmin'])
        elif 'relminmax' in req: valMIN = valNOM * (1 - req['relminmax'])

        if 'max' in req: valMAX = req['max']
        elif 'delmax'    in req: valMAX = valNOM + req['delmax']
        elif 'delminmax' in req: valMAX = valNOM + req['delminmax']
        elif 'relmax'    in req: valMAX = valNOM * (1 + req['relmax'])
        elif 'relminmax' in req: valMAX = valNOM * (1 + req['relminmax'])




        # ### START
        
        tunevar = req['tunevar']   # e.g. At, eR, M1, mA


        #if var == 'mh': thisval = slhaObj[0]['MASS'].entries[25]  # yepp...

        if var in ['DMwCoVSexp']:
            thisval = s.res_ds[var]
        elif var in ['mgen3']: 
            thisval = (slhamass[1000005]+slhamass[2000005]+slhamass[1000006]+slhamass[2000006]) / 4.
        else: 
            # usually var is a particle h, T1, eR, .., N1 [then ok] ; but also DMwCoVSexp[hmm]
            thisval = slhamass[txt2pdgId_plain[var]]

            
        # Hm, does the cutoff below work? (Maybe. Will start working again if istuned[var] changes [due to other parameters]
        if istuned.get(var,False):
            istuned[var] = ( valMIN <= thisval <= valMAX )
            verbose[var] = '%s was in tune last time, trying to stay so' %(var);
            if s.VB>2: print(verbose[var], '      -  Latest istuned[var]: %s  (val=%.3f)' %(istuned[var], thisval))
            return SUSYparCopy[tunevar]
        # ---

        istuned[var] = ( valMIN <= thisval <= valMAX )


        # -------

        y3 = valNOM

        # borrow the previous values
        x1 = req.get('x2',0)
        y1 = req.get('y2',0)


        # The values of the last calculation
        y2 = thisval
        x2 = SUSYparCopy[tunevar]
        req['x2'] = x2
        req['y2'] = y2


        # CONVERGENCE?
        # if (valMIN <= thisval <= valMAX:  # NEW
            # verbose[var] = '%s   Converged %s = %.3f with %s = %.3f' %(opt['verboseLoop1'], var, thisval, tunevar, x2)
            # return 1  # constraint satisfied
        #  sdf
        


        # Not yet ok, so need to iterate

        
        # if status['ntries'] == 1:
        if firstIteration: 
            # The first iteration is in principle a guess
            # Below we improve the guess by getting the sign correct as well as the step (for higgs)
            # Note: the step may be quite dependent
            # NB: becomes a bit messy. Also is not very much validatede

            # need to set x3

            if var == 'h':  # rough and (too?) quick to find the sign (increase or decrease At in first step) # tunevar is mixing
                thesign = 1
                z_Atmax = sqrt(6. * slhaObj[0]['MASS'].entries[1000006] * slhaObj[0]['MASS'].entries[2000006])
                #print 'max: %.2f' %(z_Atmax)
                if x2 < z_Atmax and y2 < valMIN: thesign = 1
                elif x2 > z_Atmax and y2 > valMIN: thesign = 1
                else: thesign = -1
                # Then try to decide the step (probably quite value dependent)
                # del(mh)/del(At) is ~0.002 : so del(At) = del(mh) / 0.002 (roughly). Maybe also very value dependent
                x3 = x2 + thesign * abs(valNOM-y2) / 0.002  
                print('%s startvalue %s : %.3f' %(var, tunevar, x3))
                # x3 = x2 + thesign * req['step']  # default 

                
            elif var in ['T1','eL','eR']:  # then tunevar is ...
                x3 = valNOM = s.GetValueFromReq(req['nom'], slhamass)  # this sets e.g. L1 to N1*0.5+N2*0.5, which is good starting point
                print('%s startvalue %s : %.3f     N1:%.1f  N2:%.1f' %(var, tunevar, x3, slhamass[1000022],slhamass[1000023]))

            elif var == 'N1':  # then tunevar is 'M1'
                x3 = 125.5/2  #+0.5678  # first idea.. (if I put x3=125.5/2, and also have it 
                print('%s startvalue %s : %.3f' %(var, tunevar, x3))

            elif var in ['DMwCoVSexp']:  # AfunnelUpDS, AfunnelDownDS, .. ; tunevar is mA
                x3 = 2 * slhamass[1000022] - 0.1 + req['step']
                print('%s startvalue %s : %.3f' %(var, tunevar, x3))
        
            elif var in ['mgen3']: 
                x3 = (slhamass[1000005]+slhamass[2000005]+slhamass[1000006]+slhamass[2000006]) / 4.

            else: 
                x3 = slhamass[txt2pdgId_plain[var]] + require[var]['step']   # hack
                if s.VB>3: print('DEBUG::  TuneVariable:firstIteration: %s  x3=%.5f' %(var,x3))
                

        else: 

            # Newton's method (very fast) ; the steps beyond the first are automatic
            if s.VB>3 and var in ['h','N1','mA']: print('Hm   x3:  %12s  =  %.5f + ( %.5f - %.5f ) * ( %.5f - %.5f ) / ( %.5f - %.5f ) '   %('',x2,y3,y2,x2,x1,y2,y1))


            #if istuned[var]:  # this didn't work .. maybe because it was slightly off-tuned by other parameters
            #                  # but never got back into the game 
            #    if s.VB>2: print 'TuneVariable: var %s is already tuned, setting delx3 to 0, input value will stay as it was'
            #    delx3 = 0
            
            if y2 - y1 == 0:
                if s.VB>2: print('TuneVariable: var %s : y2 - y1 = 0  (Looks like it is already very tuned.)' %(var))
                delx3 = 0   # 2014-02-19
                #return 'error' # hack: is picked up
                if not istuned[var]: delx3 = req['step']*0.5  # hack to avoid stopping at (bad) point which changes very slowly
            else: 
                delx3 = (y3 - y2) * (x2 - x1) / (y2 - y1)

            

            # hack to try to prevent very quickly changing quantities from diverging
            if 'maxstepAbs' in req and abs(delx3) > req['maxstepAbs']:
                if delx3 > 0: delx3 = req['maxstepAbs']
                else: delx3 = -req['maxstepAbs']

            # If was True last time (and we here again because we are waiting for some other requirement to be satisfied)
            # then should just keep the value as it was. (If, in the end it changes, that will be taken care of in next round, I think.) 
              
            
            x3 = x2 + delx3
            
            if s.VB>3 and var == 'h': print('Hmmm x3:  %12.5f  =  %.5f + ( %.5f - %.5f ) * ( %.5f - %.5f ) / ( %.5f - %.5f ) ' %(x3,x2,y3,y2,x2,x1,y2,y1))
            if s.VB>3: print('TuneVariable end:   [%s=] x1: %.5f  x2: %.5f  x3: %.5f    |  [%s=] y1: %.5f  y2: %.5f  y3: %.5f' %(tunevar, x1,x2,x3,  var,y1,y2,y3))

            # ------- 


        verbose[var] = ' istuned:%-5s | %s in ( %.3f , %.3f )  x1: %.1f  y1: %.3f | x2: %.1f  y2: %.3f | x3: %.1f  y3: %.3f  ' %(istuned[var], var, valMIN, valMAX,  x1,y1, x2,y2, x3,y3)   

        #x3 = 0.1*x3 + 0.9*x2  # Hack to go slowly
        
        return x3   # returns the suggested input value. Will be used if not all are tuned
                    
        
        
    # ##########
    def ExecuteCalculator(s, fn_base_out, SUSYpar, par_steer, fnmode, deltmp):
        if s.VB>5: print('SUSYpar (before Execute): ', SUSYpar)
        # could just use s.par_steer instead of inputting it

        #prevar, preval = s.GetPrevarval(SUSYpar)
        prevar, preval = s.opt['scanvars']['all'], s.opt['scanvals']['all']   # hack


        if s.calculator in ['susyhit']: 
            susyhit_tool.susyhit_execute_in_subdir(fn_base_out=fn_base_out+'_susyhit', SUSYpar=SUSYpar, par_steer=par_steer, fnmode=0, deltmp=s.susyhit_deltmp)
            
        if s.calculator in ['softsusy']:
            #softsusy_tool.softsusy_execute(fn_base_out=fn_base_out+'_softsusy', SUSYpar=SUSYpar, par_steer=par_steer)
            softsusy_tool.softsusy_execute(optD={'fn_base_out':fn_base_out+'_softsusy', 'SUSYpar':SUSYpar, 'par_steer':par_steer, 'prevar':prevar, 'preval':preval})

        if s.calculator in ['isasusy']:
            #fnID_isasusy = fnID + '_isasusy'
            opt = {'isasusy':1, 'isa2lha':2, 'isaplay':1, 'darksusy':0, 'status':1, 'prevar':prevar, 'preval':preval}
            isasusy_tool.isasusy_execute(fnID=fn_base_out+'_isasusy', par=SUSYpar, opt=opt)



    # ##########
    def GetPrevarval(s, SUSYpar):
        stdvals = []
        for zvar in s.dict['stdvars']: stdvals.append(s.dict['stdvars_format'] %(SUSYpar.get(zvar,-1)))
        prevar = ''
        preval = ''
        for zvar in s.dict['stdvars']: prevar += '  %s' %(zvar)
        for zval in stdvals: preval += '  %s' %(zval)
        prevar = prevar.strip() ;
        preval = preval.strip()
        
        return prevar, preval


    # ##########
    def DoScenario(s, SUSYpar, opt={}):
        # if s.VB: print 'DoScenario  %s' %(opt.get('verboseLoop1',''))

        SUSYparCopy = dict(SUSYpar)
        for var in s.setequalto: # 2013-12-21  (have it also below, in the iterations)
            SUSYparCopy[var] = SUSYparCopy[s.setequalto[var]]
        # Calculate spectrum 
        # a) Write slha_in file (for suspect)
        # b) Execute external code
        # c) Change filenames at wish

        #if 'hfunnelM1' in s.configsetup: # commented out: start value is already set in TuneVariable, so not needed to have it here
        #    SUSYparCopy['M1'] = 125.9/2  # HACK: start value. NB: tuning fails if I give same value here as startvalue in TuneVariable
        
        fnID = opt.get('fnID', 'agrid')  # this remains the clean fnID (without any _softsusy or similar)


        require = dict(s.require)  # make a copy because we might change it below

        

        # A0) Make input-parameter 2-liners (good to have before loop, will then still have if things fail, and can use..) 
        if s.dict['maketwoliner_scanvars']: 
            if s.VB>2: print('DEBUG: Now maketwoliner_scanvars ... (first)')
            for scanvars in s.scanvars_defs:
                outs = lineup([opt['scanvars'][scanvars] , opt['scanvals'][scanvars]])
                fn = '%s_scanvars_%s.tlt' %(fnID, scanvars)
                WriteToFile(outs=outs, fn=fn, VB=s.VB-1)
                

        # A - MAIN CALCULATION

        if 1:  # re-indent
            
            if s.VB>4: print('Doing main calculation: %s' %(s.calculator))

            #if s.calculator == 'susyhit': fn_lha = '%s_susyhit_slha.out' %(fnID)   # 
            if s.calculator == 'susyhit'  : fn_lha = '%s_susyhit_susyhit_slha.out' %(fnID)   # 
            if s.calculator == 'softsusy' : fn_lha = '%s_softsusy.slha' %(fnID)
            if s.calculator == 'isasusy'  : fn_lha = '%s_isasusy.slha' %(fnID)
            

            if s.calculator in ['susyhit'] and susyhit_tool.gen12nonequal(SUSYpar, calculator=s.calculator):
                s.warn.append(susyhit_tool.gen12nonequal(SUSYpar, calculator=s.calculator))
                if s.VB>0: print('Warning::signalgrid_loop  %s' %(s.warn[-1]))
                

            # IF REQUIRE==0 AND DERIVEDAFTERGO1==0  ***********
            if (not require) and (not s.derivedAfterGo1):
                # Do first and only calculation
                #susyhit_tool.susyhit_execute_in_subdir(fnID=fnID, SUSYpar=SUSYparCopy, par_steer=par_steer, fnmode=0, deltmp=s.susyhit_deltmp)
                s.ExecuteCalculator(fn_base_out=fnID, SUSYpar=SUSYparCopy, par_steer=s.par_steer, fnmode=0, deltmp=s.susyhit_deltmp)


            # IF REQUIRE==0 AND DERIVEDAFTERGO1==1  ***********
            if (not require) and (s.derivedAfterGo1):
                # Do first calculation
                #susyhit_tool.susyhit_execute_in_subdir(fnID=fnID, SUSYpar=SUSYparCopy, par_steer=par_steer, fnmode=0, deltmp=s.susyhit_deltmp)
                s.ExecuteCalculator(fn_base_out=fnID, SUSYpar=SUSYparCopy, par_steer=s.par_steer, fnmode=0, deltmp=s.susyhit_deltmp)
                #fn_lha = '%s_susyhit_slha.out' %(fnID)   # 

                # Read slha
                slhaObj = pyslha.readSLHAFile(fn_lha)
                slhamass = slhaObj[0]['MASS'].entries
                # Do the changes
                zmssmChanges = s.DerivedAfterGo1(slhamass, origSUSYpar=SUSYpar)
                for var in zmssmChanges: SUSYparCopy[var] = zmssmChanges[var]
                # Do second and final calculation
                #susyhit_tool.susyhit_execute_in_subdir(fnID=fnID, SUSYpar=SUSYparCopy, par_steer=par_steer, fnmode=0, deltmp=s.susyhit_deltmp)
                s.ExecuteCalculator(fn_base_out=fnID, SUSYpar=SUSYparCopy, par_steer=s.par_steer, fnmode=0, deltmp=s.susyhit_deltmp)
                del slhaObj



            # IF REQUIRE==1                         ***********
            if (require): 
                # HERE BEGINS PREPS FOR REQUIRE: PUT INSIDE IF-BLOCK
                status = {}
                # status['isok'] = 0
                ntries = 0
                #ntriesMAX = 10
                # lha = {}  # not needed
                
                isokAll = 0
                istuned = {}
                for var in list(require.keys()): istuned[var] = 0  # init
                hasbeentunedonce = {}
                for var in list(require.keys()): hasbeentunedonce[var] = 0  # init
                if s.VB>3: print('DoScenario: require.keys() = %s' %(list(require.keys())))
                while not isokAll and ntries < s.dict['require_ntriesMax']: 
                    if s.VB: print('DoScenario  ntries = %i' %(ntries))
                    ntries += 1
                    
                    # =================================
                    if 'DMwCoVSexp' in s.require and ntries == 1: SUSYparCopy['mA'] = 4 * min(SUSYparCopy['M1'],SUSYparCopy['M2'],SUSYparCopy['mu'])  # Hack to first do with a non-funnel value of mA (to be analysed ~40 lines below)
                    
                    #susyhit_tool.susyhit_execute_in_subdir(fnID=fnID, SUSYpar=SUSYparCopy, par_steer=par_steer, fnmode=0, deltmp=s.susyhit_deltmp)
                    s.ExecuteCalculator(fn_base_out=fnID, SUSYpar=SUSYparCopy, par_steer=s.par_steer, fnmode=0, deltmp=s.susyhit_deltmp)
                    #fn_lha = '%s_susyhit_slha.out' %(fnID)
                        
                    # -- When failing [pMSSM issue, e.g. M1=0] the line above gives output close to:  
                    #  invalid number: incomprehensible list input
                    #  apparent state: unit 22 named slhaspectrum.in
                    #  last format: list io
                    #  lately reading sequential formatted external IO
                    #  sh: line 1: 13270 Aborted                 /net/abel-evs/cargo/fysepf/borgeg/alt/prog/susyhit/run > /dev/null
                    # -- 
                    
                    
                    # Read slha
                    try: slhaObj = pyslha.readSLHAFile(fn_lha)
                    except:
                        print('problem reading back %s' %(fn_lha))
                        os.system('ls -l %s' %(fn_lha))
                        continue
                    slhamass = slhaObj[0]['MASS'].entries


                    # IF DERIVED  # 2014-01-20
                    if s.derivedAfterGo1:
                        zmssmChanges = s.DerivedAfterGo1(slhamass, origSUSYpar=SUSYpar)
                        for var in zmssmChanges: SUSYparCopy[var] = zmssmChanges[var]
                        

                    # ================================= 
                    # run DS if needed (used in tuning)
                    if 'DMwCoVSexp' in s.require:

                        # 
                        zending = '_susyhit_slha.out'
                        fn = '%s%s' %(fnID, zending) 
                        opts = "  -prevar '%s'  -preval '%s'" %(opt['scanvars']['all'], opt['scanvals']['all'])


                        # Hack: first try, need to place ourself in the funnel region
                        if ntries == 1:

                            # First get and check DM value of the (surely) non-funnel mA already calculate above
                            fnadd = '_tmp%i' %(random.randrange(1e4,1e5))
                            cmd = "ls %s | %s  --plha  -vb 0  -cdmcoann %i  -repl %s,''  -rundirtag %s  %s  -fnadd %s  --nods  --nods0" %(fn, s.darksusy_exe, s.ds_opt_cdmcoann, zending, fnID, opts, fnadd)
                            os.system(cmd)
                            fn_dstlt = '%s_ds%s.tlt' %(fn.replace(zending,''), fnadd)
                            s.res_ds = tlt2dict(fn_dstlt)  # under-the-hood for now
                            if s.VB>3: os.system('cat %s' %(fn_dstlt))
                            if s.VB>3: print('   (with mA = %.3f .. this was the (zeroth) off-funnel calc)' %(SUSYparCopy['mA']))
                            os.remove(fn_dstlt)

                            # PRETEST #1: if DM is already below exp, we will never reach exp
                            if s.res_ds['DMwCoVSexp'] < s.require['DMwCoVSexp']['nom'] + s.require['DMwCoVSexp']['delminmax'] :
                                if s.VB>1: print('INFO::DoScenario %s has DM satisfied even without Afunnel. Skipping.' %(fnID))
                                cmd = 'touch  DMokwithoutAfunnel__%s' %(fnID)
                                os.system(cmd)
                                return -5  # this breaks off 
                            # --- 
                            
                            # Now do the first calculation (the zeroth was done earlier), this is in the funnel
                            SUSYparCopy['mA'] = 2*slhamass[txt2pdgId_plain['N1']] - 0.1  # NB security hack fragile?
                            #susyhit_tool.susyhit_execute_in_subdir(fnID=fnID, SUSYpar=SUSYparCopy, par_steer=par_steer, fnmode=0, deltmp=s.susyhit_deltmp)
                            s.ExecuteCalculator(fn_base_out=fnID, SUSYpar=SUSYparCopy, par_steer=s.par_steer, fnmode=0, deltmp=s.susyhit_deltmp)
                            #fn_lha = .. given earlier
                            try: slhaObj = pyslha.readSLHAFile(fn_lha)
                            except: continue
                            slhamass = slhaObj[0]['MASS'].entries
                        # -----
                        
                        zending = '_susyhit_slha.out'
                        fn = '%s%s' %(fnID, zending) 
                        opts = "  -prevar '%s'  -preval '%s'" %(opt['scanvars']['all'], opt['scanvals']['all'])

                        fnadd = '_tmp%i' %(random.randrange(1e4,1e5))
                        cmd = "ls %s | %s  --plha  -vb 0  -cdmcoann %i  -repl %s,''  -rundirtag %s  %s  -fnadd %s  --nods  --nods0" %(fn, s.darksusy_exe, s.ds_opt_cdmcoann, zending, fnID, opts, fnadd)  # 2013-12-21 added rundirtag: can now run in parallel (earlier not)
                        os.system(cmd)

                        fn_dstlt = '%s_ds%s.tlt' %(fn.replace(zending,''), fnadd)
                        s.res_ds = tlt2dict(fn_dstlt)  # under-the-hood for now
                        if s.VB>3: os.system('cat %s' %(fn_dstlt))
                        if s.VB>3: print('   (with mA = %.3f)' %(SUSYparCopy['mA']))
                        os.remove(fn_dstlt)


                        # PRETEST #2: if DM is not satisfied even with Afunnel, we will never reach exp
                        #   ... This only happens for mu>~3TeV or so, I think, so will not fire often
                        #       Ex: M1=600,M2=1E3,mu=3E3: DM/exp=2 , while for M1=300 DM=0.7 or so
                        #       Should have a similar hook for hfunnelUp/Down
                        if ntries == 1: 
                            if s.res_ds['DMwCoVSexp'] > s.require['DMwCoVSexp']['nom'] + s.require['DMwCoVSexp']['delminmax'] :
                                if s.VB>1: print('INFO::DoScenario %s fails DM requirement in peak-Afunnel. Skipping.' %(fnID))
                                cmd = 'touch  DMnotokwithAfunnel__%s' %(fnID)
                                os.system(cmd)
                                return -6  # this breaks off 
                        # --- 

                    # ================================= 
                    


                    if s.VB > 3: print('\n')
                    
                    isokAll = 1  # start value
                    newval = {}
                    verbose = {}
                    # FIRST go through each requirement.
                    #   Note: changes/improvements in EACH constraint variables is continued until ALL constraints are satisfied.
                    #         Otherwise, parameter changes from the undone constraints can push the tuned one out of tune
                    #         (And then there are problems to restart Newton's method.) 
                    for var in list(require.keys()):
                        if 'okiftunedonce' in s.require_mode  and  hasbeentunedonce[var]: continue

                        if s.VB>2: print('entering TuneVariable for %s' %(var))
                        zrestune = s.TuneVariable(var, SUSYparCopy, slhaObj, require, (ntries==1), verbose, opt, istuned)
                        if s.VB>2: print('ending TuneVariable for %s : zrestune = %s' %(var, zrestune))
                        # Note: TuneVariable does not issue susyhit
                        if str(zrestune) == 'error':
                            print('TUNING FAILED (1) ... returning ... without writing scenario: touching a file bad__%s' %(fnID))
                            cmd = 'touch  bad__%s' %(fnID)
                            os.system(cmd)
                            return -1 # error
                    
                        newval[var] = zrestune
                        isokAll *= istuned[var]
                        if istuned[var]: hasbeentunedonce[var] = 1
                        if s.VB>1: print('INTERMEDIATE  %2i  %-4s  : %s' %(ntries, var, verbose[var]) + '  |||  h: %.3f    At: %.5f   N1: %.3f' %(slhamass[25],SUSYparCopy['At'],slhamass[1000022]))
                        # TuneVariable can alter: require, verbose  (no longer SUSYparCopy, that is done outside, here)


                    # THEN update SUSYparCopy if not all constraints are satisfied
                    if not isokAll:  
                        for var in list(newval.keys()):
                            ztunevar = s.require[var]['tunevar']
                            if s.VB>3: print('inserting for %s new %s : %.5f' %(var, ztunevar, newval[var]))
                            SUSYparCopy[ztunevar] = newval[var]
                            if var in ['mgen3']: SUSYparCopy['bR'] = SUSYparCopy['tR'] = newval[var]   # fairly ugly hack
                        for var in s.setequalto:
                            SUSYparCopy[var] = SUSYparCopy[s.setequalto[var]]  # 2013-12-21, DGemtR50
                   
                    else:
                        if s.VB>1: print('All constraints satisfied.\n')

                        # If the requirements (higgs, .. ) are all ok on the first try,
                        #   then any derivedAfterGo1 are not taken into account.
                        # Then need to do them here and recalculate
                        if s.derivedAfterGo1 and ntries == 1:
                            if s.VB>0: print('INFO:: requirements satisfied on first go, but need s.derivedAfterGo1 (recalculating with susyhit_tool)')
                            #susyhit_tool.susyhit_execute_in_subdir(fnID=fnID, SUSYpar=SUSYparCopy, par_steer=par_steer, fnmode=0, deltmp=s.susyhit_deltmp)
                            s.ExecuteCalculator(fn_base_out=fnID, SUSYpar=SUSYparCopy, par_steer=s.par_steer, fnmode=0, deltmp=s.susyhit_deltmp)
                            # Read the new ? Don't think we need it
                            #try:
                            #    slhaObj = pyslha.readSLHAFile(fn_lha)
                            #except:
                            #    isokAll
                            #slhamass = slhaObj[0]['MASS'].entries
                            
                        
                    # ----- 
                    del slhaObj
                    
                    
                    # VERBOSITY 
                    out = '%s | %i %3i ' %(opt['verboseLoop1'], isokAll, ntries)
                    for var in list(verbose.keys()): out += verbose[var]
                    # if s.VB: print out
                    
                    
                    
                    
                # ---- end(while) - the loop to fulfill requirements is now ended. Typically with success (isok)



                if not isokAll:
                    print('TUNING FAILED (2) ... returning ... without writing scenario: touching a file bad__%s' %(fnID))
                    cmd = 'touch  bad__%s' %(fnID)
                    os.system(cmd)
                    return -2 # error

                # ### HERE ENDS REQUIRE==1 BLOCK



            # ################################### END OF SUSYHIT CALCULATIONS: NOW MAKE DERIVED FILES (MIGHT ALSO MAKE OUTSIDE)

            #fn_slha = '%s.slha' %(fnID)
            fn_slha = '%s_%s.slha' %(fnID, s.calculator)

            
            # d1) set link
            if s.calculator in ['susyhit','suspect']:  # probably different for suspect alone
                cmd = "ln -sf  %s_susyhit_susyhit_slha.out  %s" %(fnID, fn_slha)  # might consider making the .slha one the real and the other the link
                if s.VB>4: print(cmd)
                os.system(cmd)


            if s.calculator in ['softsusy']: 
                #cmd = "ln -sf  %s_softsusy.slha  %s" %(fnID, fn_slha)  # onto itself??
                #if s.VB>4: print cmd
                #print 'now doing it ...'
                #os.system(cmd)
                pass


            if s.calculator in ['isasusy']:  # not yet implemented
                #cmd = "ln -sf  %s_isasusy.slha  %s" %(fnID, fn_slha)  # tja, would be the same .. 
                #if s.VB>4: print cmd
                #os.system(cmd)
                pass

            # ###############
            

            ### Prepare stdvars/stdvals ; lists and text
            #stdvals = []
            #for zvar in s.dict['stdvars']: stdvals.append(s.dict['stdvars_format'] %(SUSYparCopy.get(zvar,-1)))
            #stdvarsT = '' ; stdvalsT = ''
            #for zvar in s.dict['stdvars']: stdvarsT += '  %s' %(zvar)
            #for zval in stdvals: stdvalsT += '  %s' %(zval)
            #stdvarsT = stdvarsT.strip() ; stdvalsT = stdvalsT.strip()
            #s.prevar = stdvarsT
            #s.preval = stdvalsT
            #s.prevar, s.preval = s.GetPrevarval(SUSYparCopy)
            s.prevar, s.preval = opt['scanvars']['all'], opt['scanvals']['all']  # post 2014-10-08
            ###
            


        # ================= BEGIN dotoo

        # B - ISASUSY (ONLY FOR COMPARISON / CHECKS) <-- allow to use different versions
        # ----------
        if 'isasusy' in s.calculators_dotoo:
            if s.VB>2: print('Doing Isasusy too')
            fnID_isasusy = fnID + '_isasusy'
            # this are default, can switch off things (isaplay, isa2lha)
            opt = {'isasusy':1, 'isa2lha':2, 'isaplay':1, 'darksusy':0, 'status':1, 'prevar':s.prevar, 'preval':s.preval}  
            #opt = {}
            isasusy_tool.isasusy_execute(fnID=fnID_isasusy, par=SUSYparCopy, opt=opt)  



        # C - SOFTSUSY 
        # ------------
        if 'softsusy' in s.calculators_dotoo:
            if s.VB>2: print('Doing Softsusy too')
            fnID_softsusy = fnID + '_softsusy'
            #cmd = "softpoint_3.4.0.x leshouches  <  %s.slha  >  %s.slha" %(fnID, fnID_softsusy)
            #os.system(cmd)
            #softsusy_tool.softsusy_execute(fn_base_out=fnID_softsusy, SUSYpar=SUSYparCopy, par_steer=s.par_steer)
            softsusy_tool.softsusy_execute(optD={'fn_base_out':fnID_softsusy, 'SUSYpar':SUSYparCopy, 'par_steer':s.par_steer, 'fn_add':'_softsusy', 'prevar':s.prevar, 'preval':s.preval} )
            # There will be warnings ... maybe need a wrapper
            # Often no spectrum will be produced
            
            

        # D - SUSPECT standalone/newest 
        # ------------
        if 'suspect' in s.calculators_dotoo:
            if s.VB>2: print('Doing Suspect too')
            fnID_suspect = fnID + '_suspect'
            cmd = 'suspect2.43.sh  %s  %s_suspect' %(fn_slha, fnID)   # produces file <fnID>_suspect_out and <fnID>_suspect_slha
            if s.VB>3: print(cmd)
            os.system(cmd)




        # E - SUSYHIT standalone/newest 
        # ------------
        if 'susyhit' in s.calculators_dotoo:
            if s.VB>2: print('Doing Susyhit too')
            fnID_susyhit = fnID + '_susyhit'
            susyhit_tool.susyhit_execute_in_subdir(fn_base_out=fnID_susyhit, SUSYpar=SUSYpar, par_steer=s.par_steer, fnmode=0, deltmp=s.susyhit_deltmp)   # 
            cmd = 'mv %s_susyhit_slha.out %s.slha' %(fnID_susyhit, fnID_susyhit)
            if s.VB>3: print(cmd)
            os.system(cmd)



        # ================= END dotoo


        # BEGIN Derivatives of main calculation  (maybe make separate method)


        # e) Make suspect2 one-liner  (more elegant to import than execute, but not necessary)  # susyhit/suspect
        #if s.calculator in ['susyhit','suspect']: 
        if 'susyhit' in s.calculators_dotoo + [s.calculator]: 
            if s.VB>2: print('DEBUG: Now suspect2errlowtune ...')
            #opts = ''
            #if s.dict['errlowtune_addscanvars']:
            #    opts += "  -prevar '%s'  -preval '%s'" %(opt['scanvars'][s.scanvars_def], opt['scanvals'][s.scanvars_def])
            opts = "  -prevar '%s'  -preval '%s'" %(s.prevar, s.preval)
            if s.VB>3: print(opts)
            cmd = "suspect2errlowtune.py  -f %s_susyhit_suspect2.out  %s" %(fnID, opts)   # works without abs path since is in PATH
            os.system(cmd)


            

        # j) Make tlts
        for calc in [s.calculator] + s.calculators_dotoo :

            thisfn_slha = '%s_%s.slha' %(fnID, calc) 
            
            if s.dict['slha2tlt'] >= 2 or (s.dict['slha2tlt'] >= 1 and calc == s.calculator): 
                if s.VB>2: print('DEBUG: Now slha2tlt %s ...' %(calc))
                #cmd = 'ls %s | slha2tlt.py  -tables %s  %s  -vb 0' %(thisfn_slha, s.dict['slha2tlt_tables'], s.dict['slha2tlt_opt'])
                zprevarvals = "-prevar '%s'  -preval '%s' " %(s.prevar, s.preval)  # pre 2014-10-09
                cmd = "ls %s | slha2tlt.py  -tables %s  %s  %s -vb 0" %(thisfn_slha, s.dict['slha2tlt_tables'], zprevarvals, s.dict['slha2tlt_opt'])
                if s.VB>4: print(cmd)
                os.system(cmd)


            # k) Make micromegas
            if s.dict['micromegas'] >= 2 or (s.dict['micromegas'] >= 1 and calc == s.calculator):
                # for now produces only MO_mini.tlt ... will soon allow also the various contributions
                if s.VB>2: print('DEBUG: Now micromegas %s ...' %(calc))
                MO_opts = "  -prevar '%s'  -preval '%s'" %(opt['scanvars']['all'], opt['scanvals']['all'])
                cmd = 'ls %s | micromegasTool.py -vb 0  %s' %(thisfn_slha, MO_opts)
                if s.VB>4: print(cmd)
                os.system(cmd)



            # f) Calculate additional stuff (DarkSUSY) [or prepare for it]
            if s.dict['makedarksusy'] >= 2 or (s.dict['makedarksusy'] >= 1 and calc == s.calculator): 
                if s.VB>2: print('DEBUG: Now makedarksusy %s ...' %(calc))
                # fn = '%s_susyhit_slha.out' %(fnID)  # 2014-01-26, use instead .slha
                # fn = '%s.slha' %(fnID) 
                opts = ''
                opts += "  -prevar '%s'  -preval '%s'" %(opt['scanvars']['all'], opt['scanvals']['all'])  # use 'multi' instead of s.scanvars_def, which could also mean 'all'
                # print 'OPTS:', opts
                cmd = "ls %s | %s  --plha  -vb 0  -cdmcoann %i  -repl .slha,''  -rundirtag %s  %s" %(thisfn_slha, s.darksusy_exe, s.ds_opt_cdmcoann, fnID, opts)  # 2013-12-21 added rundirtag: can now run in parallel (earlier not)
                s.ds_script.append(cmd)
                if s.ds_donow:
                    if s.VB>2: print(cmd)
                    os.system(cmd)

        # =====

        # d) MakeOneLiner  # susyhit ... somewhat obsolete in times of tlt
        if s.dict['makeoneliner_suspect'] and s.calculator in ['susyhit','suspect']:
            if s.VB>2: print('DEBUG: Now MakeOneLiner_suspect ...')
            s.MakeOneLiner_suspect(fnID)


        # g) Optionally: make wig and out
        if s.dict['makewig'] and s.calculator in ['susyhit']: # (for now) only do for susyhit ; the hacked slha2wig has problems reading slha from softsusy
            if s.VB>2: print('DEBUG: Now makewig ...')
            #fn = '%s.slha' %(fnID)
            cmd = "ls %s | %s -vb 0" %(fn_slha, s.dict['slha2wig_py'])
            if s.VB>4: print(cmd)
            os.system(cmd)
            
            
            # some preps for maxder and maxdermass    # 2013-12-21
            fn_wig = '%s.wig' %(fnID)
            if s.dict['isaplay_maxmass'] > 0: zmaxmass = s.dict['isaplay_maxmass']
            else: zmaxmass = s.dict['isaplay_maxmass_relofheavy'] * s.pars_input_manual['heavy']   #s.heavy
            
            print() 
            
            if s.dict['makeder']:  # 2013-12-21
                fn_der = '%s.der' %(fnID)
                cmd = "%s  %s  -maxmass %.0f  -brmin %.7f  >  %s" %(s.dict['isaplay_py'], fn_wig, zmaxmass, s.dict['isaplay_brmin'], fn_der)
                os.system(cmd)
                
            if s.dict['makedermass']:  # 2013-12-21
                fn_mass = '%s.mass' %(fnID)
                cmd = "%s  %s  -maxmass %.0f  --masses  >  %s" %(s.dict['isaplay_py'], fn_wig, zmaxmass, fn_mass)
                os.system(cmd)


        # h) Prepare for Prospino LO/NLO, ... [can also do afterwards]
        # ...
            


            
        # END Derivatives of main calculation






        # Clean up (and fix the memory problems?)
        del require, SUSYparCopy

        return 0  # 0=ok


    # ##########
    def MakeOneLiner_suspect(s, fnID):  # somewhat obsolete

        #fn_ds = '%s_susyhit_slha.ds' %(fnID)  # 2014-01-27
        fn_ds = '%s.ds' %(fnID)
        txtDS = os.popen('cat %s' %(fn_ds)).readlines()[0].strip()
        txtDS2 = txtDS[txtDS.find('_slha ')+5:]
        
        
        fn_out = '%s_suspect2.out' %(fnID)
        res_out = susyhit_tool.susyhit_read_outfile(fn_out)
        cat = 'SU_LOWPAR'
        txt_suspect = ' | %9.3e  %9.3e  %9.3e ' %(res_out[cat]['delrho'],res_out[cat]['gmu-2'],res_out[cat]['BR(b->sgam)'])
        cat = 'SU_FINETUNE'
        txt_suspect += ' | %7.1f  %7.1f  %7.1f  %7.1f ' %(res_out[cat]['dmZ^2/mZ^2(mu^2)'],res_out[cat]['dmZ^2/mZ^2(B.mu)'],res_out[cat]['dmt/mt(mu)'],res_out[cat]['dmt/mt(B.mu)'])
        
        cat = 'FLAGS'
        txt_suspect += ' | %s ' %(res_out[cat]['flagstxt2'])
        
        
        txt = '%s # %s # %s' %(fnID, txt_suspect, txtDS2)
        
        fn_oneline = '%s_oneline.txt' %(fnID)
        f = open(fn_oneline,'w')
        f.write('%s\n' %(txt))
        f.close()

    

    # ##########
    def DoNonStdParSetting(s):
        # For now just a sequential test of various adhoc settings (to be able to use the same loop structure for everything

        pMSSMpar = s.pMSSM_object.par

        if s.M1welltempered == 1:  # is also possible to rename this to M1derived ; than non-'welltempered' adhocs be included too
            M2 = pMSSMpar['M2']
            mu = pMSSMpar['mu']
            M1 = s.M1welltempered * M2 * mu / sqrt(M2*M2 + mu*mu)
            SUSYpar, outs_updatestatus = s.pMSSM_object.UpdateAndGetPars({'M1':M1})
            if outs_updatestatus: s.warn += outs_updatestatus
            if outs_updatestatus and s.dict['stopifunknownparameter']:
                for out in outs_updatestatus: print(out)
                sys.exit("Fatal::signalgrid_loop  Unknown parameter")
            
            s.mode['fnID'] = 2

        

    # ##########
    def Get_fnID(s):

        pMSSMpar = s.pMSSM_object.par

        # First prepare width of vars
        # Also make textstring scanvars and scanvals
        scanvars = {}
        scanvals = {}
        scanvars['all'] = ''
        scanvals['all'] = ''
        scanvars['multi'] = ''
        scanvals['multi'] = ''
        ztxt = {}
        for var in s.pars_looporder:
            zval = s.format.get(var,'%.0f') %(pMSSMpar[var])  # defaults to int if not specified in s.format
            ztxt[var] = str(zval).zfill(s.ndig.get(var,s.ndig['default']))
            #ztxt[var] = str(int(pMSSMpar[var])).zfill(s.ndig.get(var,1))
            if pMSSMpar[var] >= s.useEifabove:
                ztxt[var] = ('%1.0E' %(pMSSMpar[var])).replace('E+0','E')  # FRAGILE for digits != 3 
            if s.VB>1: print('txt: ', var, ztxt[var], zval, pMSSMpar[var])

            # Then the scanvars
            scanvars['all'] += ('  ' + var)
            scanvals['all'] += ('  ' + ztxt[var])
            if var in s.pars_loop_multival: 
                scanvars['multi'] += ('  ' + var)
                scanvals['multi'] += ('  ' + ztxt[var])

        for z in scanvars, scanvals:
            for key in list(z.keys()):
                z[key] = z[key].strip()   # will also add 'manual' here
                


        if s.mode['fnID'] == -3:  # fully general
            fnID = '%s%s' %(s.fnstart, s.fnadd)
            for var in s.pars_looporder:
                fnID += '_%seq%s' %(var, ztxt[var])

        if s.mode['fnID'] == -2:  # fully general
            fnID = '%s%s' %(s.fnstart, s.fnadd)
            for var in s.pars_looporder:
                fnID += '_%s=%s' %(var, ztxt[var])

        if s.mode['fnID'] == -1:  # fully general
            fnID = '%s%s' %(s.fnstart, s.fnadd)
            for var in s.pars_looporder:
                fnID += '_%s_%s' %(var, ztxt[var])

        if s.mode['fnID'] == 0:
            fnID = '%s%s' %(s.fnstart, s.fnadd)
            for var in s.pars_looporder: fnID += var
            for var in s.pars_looporder: fnID += '_%s' %(ztxt[var])    

        if s.mode['fnID'] == 1:  # legacy (very limited)
            fnID = '%s%s' %(s.fnstart, s.fnadd)
            #fnID += '_TB%i_M1M2MU_%s_%s_%s' %(pMSSMpar['TB'], ztxt['M1'], ztxt['M2'], ztxt['mu'])
            fnID += '_TB%s_M1M2MU_%s_%s_%s' %(ztxt['TB'], ztxt['M1'], ztxt['M2'], ztxt['mu'])

        return fnID,scanvars,scanvals



    
    # ########## 
    def MainLoop(s):
        if s.VB>1: print("INFO::%s  Main" %(s.myname))

        nBad = 0
        nDone = 0
        nSkippedExisting = 0

        # Span the grid beforehand 
        resIteration = IterateNdim(varDict=s.pars_loop, vars=s.pars_looporder) 
        combs = resIteration['combinations']
        if s.VB:
            print('ncombs: %i' %(len(combs)))
            for var in s.pars_looporder: print('  mssmloop: %-5s  %2i  %s' %(var, len(s.pars_loop[var]), s.pars_loop[var]))

        # Do the loop
        for icomb in range(len(combs)):
            mssmChanges = combs[icomb]  # this is a dict of the interating vars, e.g. {'M1':120, 'M2':250, 'mu':200}

            if icomb > s.nBreak: continue

            SUSYpar, outs_updatestatus = s.pMSSM_object.UpdateAndGetPars(mssmChanges)
            if s.VB>4: print('start SUSYpar: ', SUSYpar)


            if outs_updatestatus: s.warn += outs_updatestatus
            if outs_updatestatus and s.dict['stopifunknownparameter']:
                for out in outs_updatestatus: print(out)
                sys.exit("Fatal::signalgrid_loop  Unknown parameter")


                s.DoNonStdParSetting()  # for now only contains M1welltempered


            opt = {}
            opt['fnID'],opt['scanvars'],opt['scanvals'] = s.Get_fnID()
            
            opt['verboseLoop1'] = ' %3i / %i | %s | %s | %s ' %(icomb, len(combs), s.pars_looporder, resIteration['indices'][icomb], opt['fnID'])  # tja
            s.opt = opt   # hack to allow inside ExecuteCalculator

            #if not ( SUSYpar['M1'] <= SUSYpar['M2'] <= SUSYpar['mu'] ): print 'HACK skipping %s' %(opt['fnID']); continue


            # Allow to not do if already done: 
            #fn_test = '%s_susyhit_slha.out' %(opt['fnID'])
            #if not s.force and os.path.exists(fn_test):
            fn_slha = '%s_%s.slha' %(opt['fnID'], s.calculator)
            if not s.force and os.path.exists(fn_slha):   # 2014-02-19 : this must be more correct than the above?
                if s.VB: print('%s   SKIPPING SINCE EXISTS (use --force to force)' %(opt['verboseLoop1']))
                nSkippedExisting += 1
                continue

            # --- 
            if 1: # use with Afunnel
                if not s.force and os.path.exists('DMnotokwithAfunnel__%s' %(opt['fnID'])):
                    if s.VB: print('%s   SKIPPING EXISTING DMnotokwithAfunnel (use --force to force)' %(opt['verboseLoop1']))
                    nSkippedExisting += 1
                    continue

                if not s.force and os.path.exists('DMokwithoutAfunnel__%s' %(opt['fnID'])):
                    if s.VB: print('%s   SKIPPING EXISTING DMokwithoutAfunnel (use --force to force)' %(opt['verboseLoop1']))
                    nSkippedExisting += 1
                    continue
            # ---

            if s.VB: print(opt['verboseLoop1'])

            # Action!
            zstatus = s.DoScenario(SUSYpar=SUSYpar, opt=opt)
            if zstatus <  0: nBad += 1
            if zstatus == 0: nDone += 1


        # Summary
        print()
        print('nBad:             %i' %(nBad))
        print('nDone:            %i' %(nDone))
        print('nSkippedExisting: %i' %(nSkippedExisting))
        print()
        print('Warnings:         %i' %(len(s.warn)))
        print()
        
        WriteToFile(fn=s.fn_warn, outs=s.warn, VB=s.VB-1)
        


    # ##########
    def ReadArg(s): 

        # ################################### ARGUMENT READING
        Arg = bkgjelstenArgReader.ArgReader(sys.argv, VB=0)

        '''
        if Arg.hasget('-alist'):  print 'a string list: ',Arg.list()
        if Arg.hasget('-alisti'): print 'an integer list: ',Arg.listI()
        if Arg.hasget('-alistf'): print 'a float list: ',Arg.listF()
        if Arg.hasget('-x'):  print 'a string: ',Arg.val()
        if Arg.hasget('-xI'): print 'an integer: ',Arg.valI()
        if Arg.hasget('-xF'): print 'a float: ',Arg.valF()
        '''

        if Arg.has(['--dict','--showdict']):
            print('DUMPING DEFAULT VALUES IN THE VARIABLE DICTIONARY, s.dict')
            print('   These can be changed by e.g. -dict I,var1,val1:var2,val2:F,var3,val3:var4,val4')
            print('   where I/F denote integer and float. In the example var2 inherits I from var1,')
            print('   while var4 is neither given as I nor F and therefore is taken as string.')
            print('   You may need to inspect the code, though, to understand how they are used.')
            print() 
            for key in s.dict:
                print('%-20s  %s' %(key, s.dict[key]))
            sys.exit()

        if Arg.has(['--h','---help','--automan']):
            os.system("cat %s | grep 'if Arg\.has'" %(sys.argv[0]))
            print("[ The above shows all implemented options. For full details inspect the code, %s  (Try also '-h') ]" %(sys.argv[0]))
            sys.exit()
            
        if Arg.has(['-h','--help','--h','-help']):
            s.showHelp()
            sys.exit()

        if Arg.hasget('-vb'):
            s.VB = Arg.valI()
            print('Verbosity level: %i' %(s.VB))

        if Arg.hasget('-mssmpars'):  # -mssmpars M1=100,M2=  # Sets non-loop parameters (common for all scenarios in loop)
            zz = Arg.list()
            for z in zz:
                var,val = z.split('=')
                if s.VB: print('mssmpars: %s=%s' %(var,val))
                s.pars_input_manual[var] = float(val)
            s.runloop = 1 

        if Arg.hasget('-mssmloop'):  # -mssmloop M1=100,200,300:M2=100,200:TB=10,30   # can take over from mssmpars as well
            zz = Arg.list(':')
            for z in zz:
                var,vals = z.split('=')
                if s.VB: print('mssmloop %s: %s' %(var, vals))
                s.pars_loop[var] = [float(number) for number in vals.split(',')]
                s.pars_looporder.append(var)
                if len(s.pars_loop[var]) > 1: s.pars_loop_multival.append(var)
                if var in s.dict['derivedAfterGo1_keys']: s.derivedAfterGo1 = 1
            s.runloop = 1

        if Arg.hasget('-runcard'):
            s.runcard =  Arg.val()
            print("using runcard %s"%s.runcard)
            if not os.path.isfile(s.runcard):
                print("Could not find run card %s in current directory"%s.runcard)
            else:
                lines = [line.rstrip() for line in open(s.runcard)]
                for l in lines:
                    if l.startswith("#"): continue
                    var,vals = l.split("=")
                    print(var,vals)
                    if var in s.pars_loop.keys():
                        print("Variable %s already specified with %s, overriding"%(var,",".join(s.pars_loop[var])))
                    s.pars_loop[var] = [float(number) for number in vals.split(',')]
                    s.pars_looporder.append(var)
                    if len(s.pars_loop[var]) > 1: s.pars_loop_multival.append(var)
                    if var in s.dict['derivedAfterGo1_keys']: s.derivedAfterGo1 = 1
            s.runloop = 1

        if Arg.has('--showpmssmpars'):
            print('Allowed pmssm vars: %s' %(s.pmssm_vars))

        if Arg.has('--showpmssmvals'):
            s.showpmssmvals = 1

        # can be dropped (use -mssmpars instead)  # Why?
        if Arg.hasget('-heavy'):   
            #s.heavy = Arg.valF()
            sys.exit('Warning::signalgrid_loop  Option -heavy <heavy> will become obsolete. Instead of -heavy: use -mssmpars heavy=<heavy>')
            #s.pars_input_manual['heavy'] = Arg.valF()
            
        
        if Arg.hasget('-mssmparstxt'):
            zz = Arg.list()
            for z in zz:
                a,b = z.split(':')
                s.pars_input_manual[a] = b

        if Arg.has('--maxstopmix'):
            s.pars_input_manual['mixstopreltomax'] = 1.

        if Arg.hasget('-mixstopreltomax'):
            s.pars_input_manual['mixstopreltomax'] = Arg.valF()
            

        if Arg.hasget('-cdmcoann'):
            s.ds_opt_cdmcoann = Arg.valI()

        if Arg.hasget('-loopvarsmode'): 
            s.LoopVarsMode = Arg.val()

        #if Arg.hasget('-gridname'): 
        #    s.gridname = Arg.val()

        if Arg.hasget('-fnstart'): 
            s.fnstart = Arg.val()

        if Arg.hasget('-fnadd'): 
            s.fnadd = Arg.val()
            if not s.fnadd.startswith('_'): s.fnadd = '_'+s.fnadd

        if Arg.has('--force'):
            s.force = 1

        if Arg.hasget('-force'):
            s.force = Arg.valI()

        if Arg.hasget('-nbreak'):
             s.nBreak = Arg.valI()

        # ------
        if Arg.hasget('-TB'):
            s.pars_loop['TB'] = Arg.listF()
            s.pars_looporder += ['TB']
            
        if Arg.hasget('-M1'):
            s.pars_loop['M1'] = Arg.listF()
            s.pars_looporder += ['M1']
            
        if Arg.hasget('-M2'):
            s.pars_loop['M2'] = Arg.listF()
            s.pars_looporder += ['M2']
            
        #if Arg.hasget(['-mu','MU']):  # 2014-01-30  removed 'MU' altogether and implemented pMSSM::CheckPars()
        if Arg.hasget(['-mu']): 
            s.pars_loop['mu'] = Arg.listF()
            s.pars_looporder += ['mu']

            
            
        if Arg.hasget('-M2mu'):
            s.pars_loop['M2'] = Arg.listF()
            s.pars_loop['mu'] = list(s.pars_loop['M2'])
            s.pars_looporder += ['M2','mu']
            
        if Arg.hasget('-M1M2mu'):
            s.pars_loop['M2'] = Arg.listF()
            s.pars_loop['mu'] = list(s.pars_loop['M2'])
            s.pars_loop['M1'] = list(s.pars_loop['M2'])
            s.pars_looporder += ['M1','M2','mu']

        # the following three allows to add to one (or more) of them
        if Arg.hasget('-M1add'):
            s.pars_loop['M1'] += Arg.listF()
            s.pars_loop['M1'].sort()  # tjok

        if Arg.hasget('-M2add'):
            s.pars_loop['M2'] += Arg.listF()
            s.pars_loop['M2'].sort()  # tjok

        if Arg.hasget('-muadd'):
            s.pars_loop['mu'] += Arg.listF()
            s.pars_loop['mu'].sort()  # tjok

        if Arg.hasget(['-fn_slhainput','-slhainput']):
            s.fn_slhainput = Arg.val()

        if Arg.hasget(['-format']):  # ex: -format TB:%.1f,M1:%5.1f
            zs = Arg.list()
            for z in zs: 
                var,val = z.split(':')
                s.format[var] = val
            
        if Arg.hasget(['-scanvars']): 
            s.scanvars_def = Arg.val()
            if s.scanvars_def not in s.scanvars_defs: sys.exit('Fatal::signalgrid_loop  -scanvars %s not among allowed ones: %s' %(s.scanvars_def, s.scanvars_defs))
            
        #if Arg.hasget('-derivedGo1_which'):
        #    s.which


        # ------



        if Arg.hasget(['-calc']):
            z = Arg.val()
            if z not in s.calculators: sys.exit('FATAL:: specified spectrum calculator %s not among known ones, %s' %(z, s.calculators))
            s.calculator = z
            

        if Arg.hasget(['-calc_dotoo','-calc_too']): 
            w = Arg.list()
            s.calculators_dotoo = []
            for z in w: #s.calculators:
                if z not in s.calculators:
                    sys.exit('FATAL:: known calculators: %s    UNKNOWN: %s' %(s.calculators, z))
                s.calculators_dotoo.append(z)


        if Arg.has(['--allspectrumcalc','--allcalcs']): 
            s.calculators_dotoo = ['susyhit','isasusy','softsusy','suspect']
            #if s.calculator in s.calculators_dotoo: s.calculators_dotoo.remove(s.calculator)  # done in PostInit
            
        
        if Arg.hasget('-ndig'):
            z = Arg.list()  # 2014-08-25: added the z part (allows several to be specified)
            for w in z:
                w = w.split(':')
                s.ndig[w[0]] = int(w[1])

        if Arg.hasget('-ndigTB'):
            s.ndig['TB'] = Arg.valI()

        if Arg.hasget('-ndigM1'):
            s.ndig['M1'] = Arg.valI()

        if Arg.hasget('-ndigM2'):
            s.ndig['M2'] = Arg.valI()
            
        if Arg.hasget('-ndigmu'):
            s.ndig['mu'] = Arg.valI()

        if Arg.hasget('-ndigpars'):
            s.ndig['M1'], s.ndig['M2'], s.ndig['mu'] = Arg.listI()

        if Arg.hasget('-require_mode'): 
            s.require_mode = Arg.list()

        if Arg.hasget('-useEifabove'): 
            s.useEifabove = Arg.valF()

        if Arg.hasget('-require'): 
            # input format: VAR,TUNEVAR,nom,step,delmin[,delmax]  # if delmax not given, the delmin is taken as both delmin and delmax
            zreqs = Arg.list(',,')
            for zreq in zreqs: 
                wreq = zreq.split(',')
                if len(wreq) not in [4,5]: 
                    sys.exit('FATAL:: irregular require-input: %s' %(zreq))
                znom = wreq.pop(0)
                s.require[znom] = {}
                s.require[znom]['tunevar'] = wreq.pop(0)
                s.require[znom]['nom'] = float(wreq.pop(0))
                s.require[znom]['step'] = float(wreq.pop(0))
                s.require[znom]['delmin'] = float(wreq.pop(0))
                if len(wreq) == 0: 
                    s.require[znom]['delmax'] = s.require[znom]['delmin']
                else: 
                    s.require[znom]['delmax'] = float(wreq.pop(0))

                if s.VB>1: print('Require-input: %s : %s' %(znom, s.require[znom]))

        
        if Arg.hasget('-setequalto'):  # Ex: -setequalto mR=eR,TR=eR
            w = Arg.list()
            for z in w:
                var1,var2 = z.split('=')
                s.setequalto[var1] = var2


        #if Arg.hasget('-SUSYpars'):
        #    zz = Arg.list()
        #    for z in zz:
        #        zvar,zval = z.split('=')
                

        if Arg.hasget('-mssmlooporder'):
            z = Arg.list()
            if len(z) != len(s.pars_loop):
                sys.exit('Fatal: mismatch between pars_looporder= %s  and pars_loop= %s' %(s.pars_looporder, s.pars_loop)) 
            s.pars_looporder = z

        if Arg.hasget(['-fnIDmode']):
            s.mode['fnID'] = Arg.valI()
        

        if Arg.hasget('-M1welltempered'):
            z = Arg.list()
            s.M1welltempered = int(z[0])
            s.M1welltempered_rel = float(z[1])
            s.pars_loop['M1'] = [0]  # <-- need to give a value, otherwise the loop will short-cut


        if Arg.hasget('-configsetup'):
            s.configsetup = Arg.list()
            for z in s.configsetup:
                if z not in s.configsetups: sys.exit('Fatal::signalgrid_loop  Unrecognised configsetup %s not in %s' %(z, s.configsetups))

        if Arg.hasget('-susyhit_deltmp'):
            s.susyhit_deltmp = Arg.valI()


        if Arg.hasget('-par_steer'):
            w = Arg.list()
            for z in w:
                var,val = z.split('=')
                s.par_steer[var] = val  # note, is string
            

        if Arg.hasget('-stdvars_add'):  # Ex: -stdvars_add eR,mR,TR
            s.dict['stdvars'] += Arg.list()

        if Arg.has('--nodm'):   
            s.dict['micromegas']   = 0
            s.dict['makedarksusy'] = 0

        
        # ----- The new general procedure for var input (should this be put into the ArgReader?)
        if Arg.hasget('-dict'):
            zs = Arg.list(':')
            # print zs
            for z in zs:
                zw = z.split(',')
                
                # First determine var type (default is string)
                ztype = 'string'
                if zw[0] in ['I']: ztype = zw.pop(0)
                elif zw[0] in ['F']: ztype = zw.pop(0)

                # Then get the key / var name and check
                key = zw.pop(0)
                if key not in s.dict:
                    # this restriction might be dropped
                    print(s.dict)
                    sys.exit('FATAL  non-existing var set with -var: %s  (%s)' %(key, zs))

                if len(zw) == 0: sys.exit('FATAL  non-allowed arg for -var: %s' %(zs))
                # The fill the dict/var 
                s.dict[key] = []  # First make a list. If only one entry, turn list into a plain value (bottom)
                for zw1 in zw:
                    zval = zw1
                    if ztype == 'I': zval = int(zw1)
                    elif ztype == 'F': zval = float(zw1)
                    s.dict[key].append(zval)
                if len(zw) == 1: s.dict[key] = s.dict[key][0]   # if just one entry, don't use list
        # ----- 



        if not Arg.AllOk():
            print('Problems...')
            s.showHelp()
            sys.exit('Fatal')
    
        # ################################### POST-INIT

        
    
############################## EXECUTE IF RUN AS SCRIPT (NOT JUST IMPORTED)
if __name__ == '__main__':
    t = signalgrid(cmd=['ReadArg','Main'])

    print("Cleaning up log-files")
    os.system("cat *_ds.tlt | lineup_removeduplicates.sh > table_ds.txt")
    os.system("cat *_elweak.tlt | lineup_removeduplicates.sh > table_elweak.txt")
    os.system("cat *_masses.tlt | lineup_removeduplicates.sh > table_masses.txt")
    os.system("cat *_errlowtune.tlt | lineup_removeduplicates.sh > table_errlowtune.txt")
    os.system("cat *_MO_mini.tlt | lineup_removeduplicates.sh > table_MO.txt")
    os.system("paste table_elweak.txt table_ds.txt table_masses.txt table_errlowtune.txt table_MO.txt > table_all_info_combined.txt")
############################## 

