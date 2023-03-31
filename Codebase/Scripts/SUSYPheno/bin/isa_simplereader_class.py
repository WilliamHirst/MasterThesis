#!/usr/bin/env python
#########################################
# b.k.gjelsten@fys.uio.no
# Aug 16, 2010: first version
# 2013-09-11: started objectifying the old "global" isa_simplereader & simplifying

#########################################

''' 

OBJECTIVE: 
  o SHOULD DEVELOP THIS INTO A FLEXIBLE, EASY-TO-USE TOOL FOR INVESTIGATING SCENARIOS based on: SLHA, xsecs(various formats)
  o Maybe this will become part of the general "SusySpaceFlyer" (Matrix) Framework

NEED TO:
  o Objectivy and secure all the dependencies, in particular those that read wig/out files (libISABR, ...)
  o Should also make a real slha-BR reader ; a-la what we have (or under the hood go via wig/out)

  o Need to simplify code (e.g. remove list-option? remove vars which relate to plotting and specific 2D-slices)
  o Should add simple/flexible dumper
  o Should add simple/flexible plotter


NOTE:
  x Have replaced class Isapro s.wig and s.isaout with s.dict['wig'] and s.dict['isaout'] for much more flexible handling & inspection : should do also for the other dicts (s.mass, ...) 
  
  

'''


import os,sys, pickle, math
from kilelib import IsInRangesT,PickleToFile,WriteToFile,SaveHistory   # used in isa_ex_DGemt.py
import bkgjelstenArgReader 

import lineup

# hack for pyroot on Nebu
if os.getcwd().startswith('/home/borgeg/'): sys.path.append('/usr/lib/root/')
# https://bugs.launchpad.net/ubuntu/+source/root-system/+bug/512600

HOME = os.getenv("HOME")

#execfile("/mn/kvant/u1/borgeg/.pythonrc2")
#execfile("/mn/kvant/u1/borgeg/bin/mytestarg_lib.py")
#import mytestarg_lib

from libISAWIG import *
import libISABR as isabr
isabr.LeptonMode = 0
isabr.QuarkMode  = 0

#from ROOT import *

from libIsaplot import *
from kilelib_ROOT import *  # for ColourRatio2D

from libDarksusy import ReadDarksusy2  # NB this is not the same routine as isa2root uses, based in it, but altered

from kilelib_ROOT import stdcol as Col
from kilelib_ROOT import NextColour
from kilelib_ROOT import LEPlimit_C1_gaugino_func


######################################### METHODS

# =================================
def inside(a,i0,i1):
    if a >= i0 and a<= i1: return 1
    else: return 0

# =================================

# =================================
def GetSliceInfo(IDs, freedict, repl=[], delim='M1M2MU_',VB=1):   # 2013-07-14

    v1=freedict['v1']
    v2=freedict['v2']
    fixvar=freedict['fixvar']
    fixval=freedict['fixval']
    v1skip = freedict.get('v1skip',[])
    v2skip = freedict.get('v2skip',[])
    v1skipI = []
    v2skipI = []
    for z in v1skip: v1skipI.append(int(z))
    for z in v2skip: v2skipI.append(int(z))
    
    fixvalI = int(fixval)
    
    res = {}
    res['v1'] = v1
    res['v2'] = v2
    res['fixvar'] = fixvar
    res['fixval'] = fixval
    res['IDs'] = IDs
    res['shortIDs'] = []
    
    
    arr = {}
    
    arr['v1'] = []
    arr['v2'] = []
    arr['fixvar'] = [] 
    arr['v1 x v2'] = []
    for iID in range(len(IDs)):
        ID1 = IDs[iID].strip()
        ID2 = ID1.split(delim).pop()
        for r in repl: ID2 = ID2.replace(r,'')

        # should now have an array
        ID3 = ID2.split('_')
        if len(ID3) != 3:
            print('Error::GetSliceInfo  SKIPPING Non-standard allID name: %s  ->  %s' %(ID1, ID3))
            continue

        res['shortIDs'].append(ID2)

        z = {}
        z['M_1'] = int(ID3[0])
        z['M_2'] = int(ID3[1])
        z['MU']  = int(ID3[2])

        # ---
        if z[fixvar] != fixvalI:  # should not happen
            if VB: print('INFO::GetSliceInfo  fixvar %s = %i not on fixval: %i : Skipping this one : %s' %(fixvar, z[fixvar], fixvalI, ID3))
            continue
        # ---

        if z[v1] in v1skipI: continue
        if z[v2] in v2skipI: continue

        if z[v1] not in arr['v1']: arr['v1'].append(z[v1])
        if z[v2] not in arr['v2']: arr['v2'].append(z[v2])
        if z[fixvar] not in arr['fixvar']: arr['fixvar'].append(z[fixvar])  # well, probably not needed
        arr['v1 x v2'].append((z[v1],z[v2]))
        

    # Post treatment
    arr['v1'].sort()
    arr['v2'].sort()
    arr['fixvar'].sort()
    arr['v1 x v2'].sort()

    res['arr'] = arr
    return res

    

# =================================
filestatus_ok_keys = ['ok_fn_out', 'ok_fn_wig', 'ok_fn_ds', 'ok_fn_prox', 'ok_fn_play']    # not so elegant, but this is a global var used by classes Isapro and IsaproScan

class Isapro:
    def __init__(s,ID, fn_isaout='', fn_wig='', fn_prox='', fn_isaplay='', fn_ds='', playmode='l', optD={}):

        s.dict = {}

        s.ID = ID
        s.fn_isout = fn_isaout
        s.fn_wig = fn_wig
        s.fn_prox = fn_prox
        s.fn_isaplay = fn_isaplay
        s.fn_ds = fn_ds
        s.playmode = playmode
        s.optD = optD

        #s.mass = []
        s.dict['mass'] = []
        s.dict['isaout'] = []
        #s.wig  = []
        s.dict['wig'] = []
        s.time  = []
        s.distance  = []
        s.width  = []
        s.xTot = -1.

        


        s.oxok = 0 # hack 2012-06-04
        s.ox = 'None' # Hack 2013-10-01 to make things work without prox

        s.VB = optD.get('VB',0)
        
        # 2012-06-09  # status codes (to keep overview) to be set to 1 in reading procedure if no problem
        s.status = {}
        for ok_key in filestatus_ok_keys: s.status[ok_key] = 0

        if fn_isaout: s.read_isaout(fn_isaout)
        else: s.dict['isaout'] = {}  # empty list when have not isaout (hack 11jan4)
        if fn_wig:
            s.read_isawig(fn_wig)
            
        if fn_wig:
            s.read_isabr(fn_wig)

        # print fn_prox
        if fn_prox and s.optD['USEPROX']:
            s.read_prox(fn_prox)
        
        if fn_isaplay:
            s.read_isaplay(fn_isaplay, playmode)
            
        if fn_ds:
            s.read_ds(fn_ds)

        

        s.tag = []

    def LSP(s):
        mMin = 9999999.
        pMin = ""
        for p in list(s.dict['mass'].keys()):
            if abs(s.dict['mass'][p]) < mMin and p not in ['grav','h','H','A','H+','t','']:  #hmm, what is the '' ? 
                mMin = abs(s.dict['mass'][p])
                pMin = p
        return pMin
        
    def read_isaout(s,fn):
        s.fn_isaout = fn
        s.dict['isaout'] = ReadIsaout(fn)
        if s.dict['isaout']: s.status['ok_fn_out'] = 1

    def read_isawig(s,fn):
        s.fn_wig = fn
        #s.dict['mass'], s.wig, s.time, s.distance, s.width = ReadIsawig(fn, mode=5)
        s.dict['mass'], s.dict['wig'], s.time, s.distance, s.width = ReadIsawig(fn, mode=5)
        if s.dict['mass']: s.status['ok_fn_wig'] = 1

    def read_isabr(s,fn):
        # print 'fn: ', fn
        ifile = open(fn); isabr.line = ifile.readlines(); ifile.close()   # ???? WHAT IS isabr.line ???? 
        #isabr.ExampleBR()
        #isabr.ExampleNewStructure()
        #isabr.GetParticles()
        #s.part = dict(isabr.part)
        if s.VB>1: print('read_isabr: fn = %s' %(fn))
        s.br = isabr.ParticleDecay()  # look in libISABR.py to see what is included
       

    def read_prox(s,fn):
        s.fn_prox = fn
        if os.path.exists(fn): s.oxok = 1
        else:
            s.oxok = 0
            #print 'WARNING::read_prox  non-existing PROX: %s' %(fn)  # 2012-10-10 maybe comment out later
            
        # the things below apparently "work", i.e. don't crash, even if fn doesn't exist
        x,xT,xInput = getOX(fn)
        if x: s.status['ok_fn_prox'] = 1
        s.ox = OX(x,xT,xInput, ID=fn)
        # s.OX = OX(getOX(fn))  #possible?
        s.xTot = s.ox.getTot()
        if s.VB>2:
            print('read_prox: %s' %(fn))
            print('  x: %s' %(x))
            
    def read_isaplay(s,fn,playmode='l'):  #NOT IN USE   #2012-06-01 not?
        s.fn_isaplay = fn
        s.play = IsaPlayRes(fn,playmode)
        s.status['ok_fn_play'] = s.play.ok
        # print 'DEBUG PLAY STATUS: ', s.play.ok
        

    def read_ds(s,fn):  
        s.fn_ds = fn
        #print fn
        s.ds = ReadDarksusy2(fn)  # NB this is not the same routine as isa2root uses, based in it, but altered
        if s.ds: s.status['ok_fn_ds'] = 1   # s.ds will be empty ({}) if file not there
        else:
            print('Warning::isa_simplereader  darksusy file not found: %s' %(fn))
            

    def Tag(s,mode=1):
         if mode == 1:
            mlR = s.dict['mass']['eR']
            if   inside(mlR, abs(s.dict['mass']['N1']), abs(s.dict['mass']['N2'])): s.tag.append('N1<lR<N2')
            elif inside(mlR, abs(s.dict['mass']['N2']), abs(s.dict['mass']['N3'])): s.tag.append('N2<lR<N3')
            elif inside(mlR, abs(s.dict['mass']['N3']), abs(s.dict['mass']['N4'])): s.tag.append('N3<lR<N4')
            elif mlR > abs(s.dict['mass']['N4']): s.tag.append('lR>N4')
            else: s.tag.append('lR<N1')  # ever happen?

    def showNNmix(s): 
        #w = s.wig
        w = s.dict['wig']
        if not w: 
            print('showNNmix: wig empty')
            return 

        print("| N1 |      | %6.3f  %6.3f  %6.3f  %6.3f |     | B  | "  %(w['n11'],w['n12'],w['n13'],w['n14']))
        print("| N2 |      | %6.3f  %6.3f  %6.3f  %6.3f |  .  | W  | "  %(w['n21'],w['n22'],w['n23'],w['n24']))
        print("| N3 |  ==  | %6.3f  %6.3f  %6.3f  %6.3f |     | H1 | "  %(w['n31'],w['n32'],w['n33'],w['n34']))
        print("| N4 |      | %6.3f  %6.3f  %6.3f  %6.3f |     | H2 | "  %(w['n41'],w['n42'],w['n43'],w['n44']))

    def getNcomp(s,iNi,ret=2,show=0): 
        comp1 = []
        comp2 = []
        comp1T = "N%i  " %(iNi)
        comp2T = "N%i  " %(iNi)
        for iNj in range(1,5): 
            mixname = 'n%i%i' %(iNi,iNj)
            #comp1.append(s.wig[mixname])
            comp1.append(s.dict['wig'][mixname])
            comp2.append(comp1[-1]**2)
            comp1T += '%6.3f  ' %(comp1[-1])
            comp2T += '%6.3f  ' %(comp2[-1])
        if show in [1,12]: print('Linear:  %s' %(comp1T))
        if show in [2,12]: print('Square:  %s' %(comp2T))
        if ret == 1: return comp1
        elif ret == 2: return comp2
        elif ret == -1: return comp1T
        elif ret == -2: return comp2T


# =================================
class IsaproScan:
    def __init__(s,dir_wig='dir_wig_dummy/',dir_prox='dir_prox_dummy/',dir_isaout='dir_isaout_dummy/', dir_isaplay='dir_isaplay_dummy/', dir_ds='dir_ds_dummy/', optD={}):
        
        s.L = []
        s.D = {}
        s.dir_wig = dir_wig
        s.dir_prox = dir_prox
        s.dir_isaout = dir_isaout
        s.dir_isaplay = dir_isaplay
        s.dir_ds = dir_ds
        s.optD = optD
        s.VB = optD.get('VB',0)
        
        s.status = {}
        for ok_key in filestatus_ok_keys: s.status[ok_key] = 0

        # print 'dir_prox: %s' %(s.dir_prox)  # this is actually set afterwards


    def Add(s, ID, fn_isaout='', fn_wig='', fn_prox='', fn_isaplay='', fn_ds='', playmode='l'):
        if not fn_isaout and s.optD['USEOUT']: fn_isaout = "%s%s%s" %(s.dir_isaout, ID, s.optD['ending_out'])
        if not fn_wig: fn_wig = "%s%s%s" %(s.dir_wig, ID, s.optD['ending_wig'])
        if not fn_prox and s.optD['withLHAending']: fn_prox = "%s%s.ox" %(s.dir_prox, ID+".lha")  #adhoc; some grids (glsq, M1M2MU) were generated with scenid ending with .lha
        if not fn_prox: fn_prox = "%s%s.ox" %(s.dir_prox, ID)
        if not fn_isaplay and s.optD['USEPLAY']:
            # fn_isaplay = "%soutSMp_%s%s_l" %(s.dir_isaplay, ID, isaplay_combcode)
            # print 'Add: playmode=%s' %(playmode)
            fn_isaplay = "%soutSMp_%s%s_%s" %(s.dir_isaplay, ID, isaplay_combcode, playmode)
            # print 'huhuh fn_isaplay: %s' %(fn_isaplay)
        #print 'hm1: ', fn_ds
        if not fn_ds and s.optD['USEDS']: fn_ds = "%s%s%s" %(s.dir_wig, ID, '.ds')
        #print 'hm2: ', fn_ds
        #print 'hm3: |%s|%s|%s|' %(s.dir_wig,ID,'.ds')

        #print 'Add l313: ', ID, fn_isaout, fn_wig, fn_prox, fn_isaplay, fn_ds, playmode
        #print 'Here fn_prox = %s' %(fn_prox)
        isapro = Isapro(ID=ID, fn_isaout=fn_isaout, fn_wig=fn_wig, fn_prox=fn_prox, fn_isaplay=fn_isaplay, fn_ds=fn_ds, playmode=playmode, optD=s.optD)
        
        #isapro = Isapro(ID, fn_isaout, fn_wig, fn_prox)
        s.L.append(isapro)
        s.D[ID] = isapro

        for ok_key in isapro.status: s.status[ok_key] += isapro.status[ok_key] 
        
        
    def ReadFromIDlist(s,IDs, playmode='l'): 
        # for ID in IDs: iScan.Add(ID)
        for iID in range(len(IDs)):
            ID = IDs[iID]
            if s.VB>0 and (math.fmod(iID,100) == 0): print("INFO::isa_simplereader::IsaproScan::ReadFromIDlist  read %i/%i" %(iID,len(IDs)))
            s.Add(ID, playmode=playmode)
            # print 'ReadFromIDlist: playmode',playmode
            # this is the place to see how many lhawig, prox, play etc. are available (can also have a dedicated one)
        if s.VB: 
            print("Read %i scenarios" %(len(IDs)))
            for ok_key in filestatus_ok_keys:
                print('   %-12s %5i' %(ok_key, s.status[ok_key]))
        

    def Show3Laccu(s,I0=0,I1=4):
        for iL in range(I0,I1+1):
            if iL > len(s.L): continue
            isapro = s.L[iL]
            print("%3i  %6.4f  %s" %(iL, isapro.play.xrelAccu[3], isapro.ID))

    def N(s): return len(s.L)
    
    def Tag(s,mode=1):
        # for isapro in iScan.L: isapro.Tag(mode)
        for isapro in s.L: isapro.Tag(mode)
            
    def ShowScenLeft(s,IDstart):
        for iS in range(s.N()):
            S = s.L[iS]
            if S.ID.startswith(IDstart):
                print("%4i  %s" %(iS,S.ID))
        
    def FindID(s, anID, match='exact', show=1, ret=''): 
        matches = ['exact','grep']
        returns = ['','object','index']
        if match not in matches:
            print('WARNING: Isaproscan::FindID: match mode must be in %s  (using exact)' %(matches))
            match = 'exact'
        iLs = []
        for iL in range(s.N()): 
            if match == 'exact' and s.L[iL].ID == anID: 
                iLs.append(iL)
            if match in ['grep'] and anID in s.L[iL].ID: 
                iLs.append(iL)

        # Both print and return result
        if show: 
            for iL in iLs: 
                # print " %4i  %s" %(iL, iScan.L[iL].ID) 
                print(" %4i  %s" %(iL, s.L[iL].ID)) 
        if 'index' in ret: return iLs
        if 'object' in ret:
            if len(iLs) >= 1: 
                if len(iLs) > 1: print('FindID: several found, returning the first')
                # return iScan.L[iLs[0]]
                return s.L[iLs[0]]
            else: 
                print('FindID: none found')
                
        
# =================================
            

# ##################################################### GLOBAL METHODS
# ##################################################### GLOBAL METHODS
# ##################################################### GLOBAL METHODS

class isa_simplereader:

    
    # ##########
    def SetExcludedScenarios(s):
        s.excludescenDG_C1C1 = [(100,100,100),(140,100,100),(250,100,100)]
        s.excludescenDG_N1C1 = [(100,100,100),(140,100,100),(250,100,100),(100,100,110),(100,110,100),(140,110,100)]
        s.excludescenDG_N1N1 = [(100,100,100),(100,100,110),(100,100,120),(100,110,100),(100,110,110),(100,120,100),(100,140,100),(140,100,100),(140,100,110),(140,110,100),(140,120,100),(250,100,100),(250,100,110),(250,110,100)]  # 2012-01-25 (also defined in isa_ex_DGemt.py)
        s.exclude14scen = excludescenDG_N1N1
        s.excludescenDG_N1 = excludescenDG_N1N1
        s.excludescenDG_C1 = excludescenDG_N1C1
        s.excludescenDG_14 = exclude14scen
        print("INFO  isa_simplereader:  Defined exclude14scen scenarios with the 14 scenarios of DGemt and DGnoL which lack values for certain subprocess cross-sections due to low lying masses")
      

    # ##########
    def ShowExampleOuterLoops(s): 
        print("\nNow you might want to do e.g.:")
        print("   execfile('isa_simpleplotter.py')")
        print("   execfile('isa_ex_plotdict.py')")
        print("   execfile('colourratio_ex1.py')")
        print("   execfile('isa_ex_DGemt.py')")
        print("   execfile('isa_definecontours.py')   # may want to execute this one first")
        # print "   execfile('')"
        print() 


    # #########
    def DGLegacy(s):  # not to be used in general cases
        if 'fixvar' in s.freedict: s.fixvar = s.freedict['fixvar']      # needed?  THESE ARE VERY SPECIFIC
        if 'fixval' in s.freedict: s.fixval = int(s.freedict['fixval'])
        s.fixvarT = s.fixvar.replace('_','')  # M_1->M1, M_2->M2, MU->MU  # 2013-07-13

        if fixvar == '' or fixvarT == '':
            sys.exit('Fatal::isa_simplereader  (fixvar,fixvarT,fixval) = (%s,%s,%i)' %(fixvar,fixvarT,fixval))

        if 'v1' not in freedict: sys.exit('Fatal::isa_simplereader  v1 not in freedict (has %s)' %(list(freedict.keys())))  # 2013-07-14
        if 'v2' not in freedict: sys.exit('Fatal::isa_simplereader  v2 not in freedict (has %s)' %(list(freedict.keys())))

        if s.contours:  # 2012-06-02
            ggCont = {'scanID':s.aninputID, 'fixvar':s.freedict['fixvar'], 'fixval':s.freedict['fixval']}
            if VB: print('Now: GetContourHistos()')
            s.hContours = GetContourHistos(gg=ggCont, contours=s.contours, Scan=s.Scan, flip=s.flip) 
            del ggCont

        s.OptDict_reader['LEPlimit_C1_gaugino'] = LEPlimit_C1_gaugino_func()   # Returns a TGraph ; can us Eval(delM)

        s.sliceInfos = []   # 2013-07-14
        for iID in range(len(s.IDs)):
            s.sliceInfos.append(GetSliceInfo(IDs=s.IDs[iID], freedict=s.freedict))

        s.asliceInfo = sliceInfos[0]
        s.freedict['sliceInfos'] = s.sliceInfos
        s.freedict['asliceInfo'] = s.asliceInfo

        #if outer_py: 
        #    print "isa_simplereader: ADHOC autoexecute outer: %s" %(outer_py)
        #    execfile(outer_py)

        if not os.path.exists(figdir):
            print('INFO::isa_simplereader  Creating figdir: %s' %(figdir))
            os.mkdir(figdir)


    # ##########
    def ExampleUsage(s):
        if s.VB>0: print('ExampleUsage: s.nScans = %i' %(s.nScans) , s.inputID)
        for iScan in range(s.nScans): 
            # Scan[iScan].ReadFromIDlist(IDs[iScan])
            print()
            print('Example output for scan %i : %s' %(iScan, s.inputID[iScan]))

            print('s.Scan: ', s.Scan)

            print("scenario number 0:")
            s.Scan[iScan].L[0].br.part['N2'].showDecays()
            
            print("scenario number %i" %(len(s.IDs[iScan])-1))
            s.Scan[iScan].L[-1].br.part['N2'].showDecays()

            print("Scan[%i].L[0].br.part['N2'].showDecays()" %(iScan))
            
            print("ex.: ipro = nScan.L[22] ")
            print("     pN2 = ipro.br.part['N2']") 

            # Now comes plotting etc. ... interactively?
            # or writing to ROOT
            # iScan.Show3Laccu()


    # ##########
    # ##########
    # ##########
    def __init__(s, cmd=[], optD={}):
        # ====================== PRE INIT
        if 'argv' in optD: s.argv = optD['argv']
        else: s.argv = sys.argv
        s.cmd = cmd
        s.myname = sys.argv[0].split('/').pop()
        s.VB = 0  # 1
        s.HOME = os.getenv('HOME')

        s.fn_globalhistory = '%s/.globalhistory__isa_simplereader_class.txt' %(s.HOME)
        SaveHistory(fn=s.fn_globalhistory, argv=s.argv, opt=['stripcmd','dirsameline','date'])

        s.cwd  = os.getcwd()  # current work directory  
        #s.dir0 = '%s/XXX' %(s.HOME)
        s.dir0 = ''

        s.dict = {}

        # s.dict['test'] = 'testval'
        # s.dict['test2'] = 'testval2'
        # s.dict['testI'] = 3
        # s.dict['testF'] = 4.34

        s.warn = []
        s.fn_warn = 'warnings.txt'
        s.fn_report = 'report'
        s.report = []

        # ##########
        # begin(old)
        # ######################################## INITIALISATION
        # dir_isaout  = 'ISAOUT/'
        # dir_wig     = 'LHAWIG/'
        # dir_prox    = 'PROX/'

        s.vars = {}

        s.dir_isaout  = [] #'isasimple/'
        s.dir_wig     = [] #'isasimple/'
        s.dir_prox    = [] #'isasimple/'
        
        s.fn_IDs = []
        s.scanID = []
        
        s.dir_isaplay = []  #'RES_ISAPLAY/'
        s.dir_ds = []
        
        # isaplay_combcode = "" # run1, run2
        s.isaplay_combcode = "_all"
        s.isaplay_combcode = "all" 
        # isaplay_combcode = "_kromo"

        s.vars['VB'] = s.VB
        s.vars['USEOUT'] = 1
        s.vars['USEPROX'] = 0
        s.vars['USEPLAY'] = 0
        s.vars['USEDS'] = 0
        s.vars['ending_out'] = ".out"
        s.vars['ending_wig'] = ".wig"
        
        s.IDs = []
        s.IDs2 = [] #dummy
        
        s.vars['withLHAending'] = 0
        
        s.inputID = []  #""
        s.fn_add0 = ""
        
        s.readervar = []
        
        s.outer_py = ''
        
        s.fntag = ''   # this can be given on command-line (i.e. here), but is also used hardcoded in various scripts (sometimes replacing the command-line given here, hence use fntag2 to ensure)
        s.fntag2 = ''  # this is meant as command-line exclusive [ though not always sent to get2DhistM1M2MU() ]
        s.fnpre = ''
        s.fnbase = ''
        
        # freedict = {'fixvar':'M_1', 'fixvarT':'M1'}
        s.freedict = {}  # 2013-07-13
        
        s.SetBatch = 0
        
        s.Scan = []
        s.scanTag = []
        s.nIDs = []
        
        s.plotinfo = ''
        s.txtinfo = []
        
        s.fixvar = '' ; fixvarT = '' ; fixval = 99
        
        s.steer = {}
        s.steerList = {}
        s.inputOptDict = {}   # not (yet?) in use ; there is maybe also another usage (intended) in isa_eg_DGemt.py 
        
        s.OptDict_reader = {}  # in use
        
        s.figdir = 'fig1'
        
        s.contours = {}
        s.hContours = {}
        
        s.flip = 0
        
        s.format = {}   # can with this send format messages to the plotting routines, e.g. how many decimals etc. 
        
        s.fn_savehistos = ''
        
        s.gDict = {}
        
        s.playmodes = ['l','q','b','T','y']
        s.playmode = 'l'   # 'l': means looks for outSMp_......all_l  and reads the right column (5)
        # 'T': means looks for outSMp_......all_T  and reads the right column (8)
        # etc. q,b,y


        # Additional (previously gloabl) variables defined after Arg reading ... many are hacks for outer loopers
        s.PLOTTER_INITDONE = 0  # dirty helper variable in leptonscan_plot.py  # probably obsolete
        s.AUTO = 0  # used to auto over isa_simpleplotter.py

        s.COL = NextColour()   # colour iterator
        #s.gTex = TLatex()

        s.hCont = []  # can be filled in isa_definecontours.py

        # end(old)

        # new:
        s.cmds = ['ReadArg','PostInit','Main','DumpScenario','ExampleUsage']

        # the more general structure
        s.dump = {}
        s.dump_headdict = {}
        s.dumpformat = {}
        s.dumpformat['mass'] = '%6.1f'
        s.dumpformat['wig'] = '%4.1f'   # only the decimal counts, the rest is taken care of by lineup
        s.dumpformat['isaout'] = '%6.1f'

        s.dict['delim_between_blocks'] = 1
        s.dict['show_scenarioname'] = 1  # D=1, but for beauty, will often switch off
        s.dict['clean'] = 0

        s.dict['save'] = 1
        s.dict['show'] = 1


        #s.dump_wig = []
        #s.dump_isaout = []

        #s.dump_masses = []
        #s.dump_masses = ['N1','N2']

        s.dump_ds = []

        s.dump_xsecs = []
        #s.dump_xsecs = ['tot', 'C1 C1', 'N2 C1']

        s.dump_brs = []
        #s.dump_brs = [['N2','eR'],['C1','eR']]

        s.dump_brs_all = []
        #s.dump_brs_all = ['N2','C1']  # give particles and all BRs will be dumped (on several lines) (showDecay())  # debug/investigation mode
        
        s.fn_dumpscenario0 = 'DumpScenario'
        
        #s.ox = []  # Hack to make things run without xsecs

        # ##########

        # ====================== READ ARG
        if 'ReadArg' in s.cmd: s.ReadArg()


        # ====================== POST INIT
        #print 'Before PostInit:  s.dir_wig: ', s.dir_wig
        if 'PostInit' in s.cmd: s.PostInit()
        #print 'After PostInit:  s.dir_wig: ', s.dir_wig

        # ====================== EXECUTE 
        if 'Main' in s.cmd: s.Main()



    



    # ##########
    def SetIDs(s):
        for iScan in range(s.nScans):
            s.IDs.append([])
            for ID in s.IDs2: 
                # Here can add a check that the file is ok?
                if s.VB>2: print('  ID: %s' %(ID.strip()))
                s.IDs[iScan].append(ID.strip())
        
            # nIDs = len(IDs)
            if s.VB > 1:
                print('s.IDs2: %s' %(s.IDs2))
                print('s.nScans: %i' %(s.nScans))
                print('iScan : %s' %(iScan))
                print('s.scanTag: %s' %(s.scanTag))
                print('s.fn_IDs: %s' %(s.fn_IDs))
                print('s.fn_IDs[%i]: %s' %(iScan, s.fn_IDs[iScan]))
            s.scanTag.append(s.fn_IDs[iScan].replace("_allID","").replace("/allIDs","").replace("allID",""))
            s.nIDs.append(len(s.IDs))
        #if s.VB>0: print 'SetIDs: s.IDs: %s' %(s.IDs)
        del s.IDs2


    # ##########
    def DefineIsaproScan(s):
        for iScan in range(s.nScans):
            if len(s.inputID) < iScan+1: s.inputID.append('Scan_%i' %(iScan))
            if s.VB>0: print("Creating Scan object %i : %s" %(iScan, s.inputID[iScan]))
    
            # zScan = IsaproScan(dir_wig=dir_wig[iScan], dir_prox=dir_prox[iScan], dir_isaout=dir_isaout[iScan])
            zScan = IsaproScan(optD=s.vars)
            # by default directories are set to dummies. Now set the relevant ones to correct values
            if s.VB>0: print('iScan: %i' %(iScan))
            if s.VB>1: print('lens: ', len(s.dir_wig), len(s.dir_prox), len(s.dir_isaout), len(s.dir_isaplay), len(s.dir_ds))
            if iScan < len(s.dir_wig): zScan.dir_wig = s.dir_wig[iScan]
            if iScan < len(s.dir_prox):
                zScan.dir_prox = s.dir_prox[iScan]
            if iScan < len(s.dir_isaout): zScan.dir_isaout = s.dir_isaout[iScan]
            if iScan < len(s.dir_isaplay): zScan.dir_isaplay = s.dir_isaplay[iScan]
            if iScan < len(s.dir_ds): zScan.dir_ds = s.dir_ds[iScan]
            

            if s.VB>1: print('Define: s.IDs: ' , s.IDs)
            zScan.ReadFromIDlist(s.IDs[iScan], s.playmode)
            s.Scan.append(zScan)

        # Hack. Sets the first scan to the default one (code is not general enough). In 99% of cases we only consider one grid
        s.aScan = s.Scan[0]
        s.aninputID = s.inputID[0]



    # ##########
    def PostInit(s): 
        s.vars['VB'] = s.VB


        if s.dict['clean']:
            s.dict['show_scenarioname'] = 0
            s.dict['delim_between_blocks'] = 0
            

        if s.dir_isaout == []: s.dir_isaout = ['./']    # 2013-11-28: 
        if s.dir_wig    == []: s.dir_wig = ['./']
        if s.dir_prox   == []: s.dir_prox = ['./']


        s.fn_dumpscenario_pickle = '%s%s.pickle' %(s.fn_dumpscenario0, s.fntag)
        s.fn_dumpscenario_txt = '%s%s.txt' %(s.fn_dumpscenario0, s.fntag)


        if s.dir0: s.fn_warn = '%s/%s' %(s.dir0, s.fn_warn)
        if s.dir0: s.fn_report = '%s/%s' %(s.dir0, s.fn_report)

        if not s.IDs2:
            if s.VB>0: print("Reading ID list from stdin ...")
            z = sys.stdin.readlines()
            s.IDs2 = []
            zdirpickup = ''
            for zL0 in z:
                zL = zL0.strip().split('/').pop().replace('.slha','')
                if zdirpickup == '': # 2013-12-20
                    zw = zL0.split('/')
                    zw.pop()
                    for zwz in zw: zdirpickup += zwz + '/'
                    #zdirpickup = zwz.rstrip('/')  # actually requires '/' at the end (old/bad convention)
                    if zdirpickup == '': zdirpickup = './'

                s.IDs2.append(zL)
            s.fn_IDs.append('from_stdin')

            if s.dir_isaout == ['./']: s.dir_isaout = [zdirpickup]  # 2013-12-20
            if s.dir_wig == ['./']: s.dir_wig = [zdirpickup]  # 2013-12-20



        #s.nScans = len(s.IDs2)
        s.nScans = len(s.fn_IDs)
        for i in range(s.nScans): 
            if len(s.inputID) <= i: s.inputID.append('inputID%i' %(i))  # then don't need to specifically set -inputid ${M2MURANGE} (can be dangerous since this variable is used in the histogram setup (old times))

        s.SetIDs()

        s.DefineIsaproScan()


        if s.cmds == []: sys.exit("Fatal::isa_simplereader_class  No commands selected (ex: -cmd DumpScenario)")

        

    # ##########
    def Plotter(s):
        #from ROOT import * 
        import ROOT
        s.ROOT = ROOT
        if 'atlas' in s.plotstyle: 
            s.gROOT.LoadMacro("/mn/kvant/u1/borgeg/scripts/sky_AtlasStyle.C")
            SetAtlasStyle()
        #from kilelib_ROOT import *  # will fail
        if s.SetBatch: s.gROOT.SetBatch()



    # ##################################################### CLASS METHODS
    # ##################################################### CLASS METHODS
    # ##################################################### CLASS METHODS

    # ##########
    def showHelp(s):
        s.oldHELP()
        print(30*'-')
        print(' Usage: %s [options]' %(s.myname))
        print('        %s  -dict test,txt1:I,testI,3:test2,txt2a,txt2b:F,testF,4.14   # for using autodict (NB: vars need to be defined in __init__)' %(s.myname))
        print('   Ex:  cd ~/grids_lsp/res_DGemt_TB6_to700')
        print('        less allIDs/allIDs_M1eq100 | tail -3 | isa_simplereader_class.py -cmd DumpScenario -lhawig lhawig  -prox PROX_7TeV_NLO_n867/  -VB 0  -dump_xsecs N2_C1,C1_C1  -dump_brs N2:eR,C1:eR,N2:N1:Z')
        print('        less allIDs/allIDs_M1eq100 | tail -3 | isa_simplereader_class.py -cmd DumpScenario -lhawig lhawig  -prox PROX_7TeV_NLO_n867/  -VB 0  -dump_xsecs N2_C1,C1_C1  -dump_brs N2:eR,C1:eR,N2:N1:Z  -dump_brs_all N2,C1')
        print('        less allIDs/allIDs_M1eq100 | head -20 | isa_simplereader_class.py -cmd DumpScenario -lhawig lhawig  -prox PROX_7TeV_NLO_n867/  -VB 0  -dump_masses N1,N2,N3,N4,C1,C2,eR,T1,h,t  -dump_xsecs tot,N1_N1,N1_C1,N1_N2,C1_C1,N2_C1,N3_C1,N1_N3,N2_N3  -dump_brs N2:eR,N2:T1,C1:eR,C1:T1,C1:N1:W,C1:N1:nue,C1:N1:nut,N3:eR,N3:T1  -fntag DGemt_TB6_M1eq100')
        print("   Ex:  less allIDs/allIDs_n100 | isa_simplereader_class.py -cmd DumpScenario -lhawig lhawig_n100/  -dump_masses T1,h,N1,N2,N3,N4,C1,C2  -dump_brs N2:T1,N2:N1:e,C1:T1,N3:T1 -fntag test  -vb 0| lineupauto  # 2013-10-01 from ~/grids_lsp/DGstau2013/")
        print("   Ex:  less allIDs/allIDs_n100 | isa_simplereader_class.py -cmd DumpScenario -lhawig lhawig_n100/  -dump_masses T1,h,N1,N2,N3,N4,C1,C2  -dump_brs N2:T1,N2:N1:e,C1:T1,N3:T1 -fntag test  -vb 0  | sed s/susy_DGstauR50_TB50_M1M2MU_// | sed s/Scenario/'M1 M2 MU'/ | sed s/_/' '/g | lineupauto   # 2013-10-01 from ~/grids_lsp_DGstau2013/  (Ready to plug directly into munch.py)")
        print("   Ex:  cat allIDs/allIDs_M1eq050_n196 | head -10 | isa_simplereader_class.py -lhawig lhawig_M1eq050  -cmd DumpScenario -dump_masses N1,N2,N3,N4,C1,C2,h  -fntag test  -vb 1  -dump isaout:M_1/M1,M_2/M2,MU::wig:alfah  # 2013-11-21  (cd ~/grids_lsp/DGnoSL_TB10_2013-09-16/")
        print("   Ex:  cat allIDs/allIDs_M1eq050_n196 | grep -v 3E3 | isa_simplereader_class.py -lhawig lhawig_M1eq050  -cmd DumpScenario -dump_masses N1,N2,N3,N4,C1,C2,h  -fntag to500  -vb 1  -dump isaout:M_2/M2,MU  -dict I,show_scenarioname,0:I,delim_between_blocks,0  -fn_dump DGnoSL  # 2013-11-22  (cd ~/grids_lsp/DGnoSL_TB10_2013-09-16/)")
        print("   Ex:  cat allID | grep -v 3E3 | isa_simplereader_class.py -lhawig lhawig_M1eq050 -cmd DumpScenario  -fn_dump DGnoSL_vardump -fntag to500  -dump isaout:M_2/M2,MU::mass:N1,N2,N3,N4,C1,C2,h  -dict I,show_scenarioname,0:I,delim_between_blocks,0  -dumpformat mass:%6.2f  # 2013-11-23")
        print("   Ex:  ls lhawig_n165/susy_DGstauR50_TB50_M1M2MU_*.slha | isa_simplereader_class.py  -dump isaout:MU,M_2/M2::mass:N1  -dump_brs h:b:b,h:N1:N1   # cd ~/grids_lsp/DGstau2013/  # 2013-12-20")
        

    # ##########
    def oldHELP(s):
        print(" OLDHELP:")
        print("  Usage:  First read info: python -i ~/bin/isa_simplereader.py  -id <file_with_list_of_scenIDs>  -dir <DIR_to_wig/out/lha>  [-prox <DIR_to_ox>")
        print("          Then plot      : execfile('isa_simpleplotter.py')")
        print("          By default eps/pdf figs are put in fig1/")
        print()
        print("   Ex: python -i ~/bin/isa_simplereader.py  -id res_xsec_glsq/allID  -dir BASKET/  -prox res_xsec_glsq/PROX/  --withlhaending")
        print("       python -i ~/bin/isa_simplereader.py  -id res_xsec_M1M2MU/allID  -dir BASKET/  -prox res_xsec_M1M2MU/PROX --withlhaending")
        print("       python -i ~/bin/isa_simplereader.py  -id res_g96/allID_1  -dir res_g96/  -ending .")
        print("       python -i ~/bin/isa_simplereader.py -ending '.'  -id allID_3_wide  -dir LHAWIG/  --play  -combcode 'kromo'")



    # ##########
    def DumpWarnings(s):
        f = open(s.fn_warn,'w')
        for out in s.warn: f.write('%s\n' %(out))
        f.close()

        
    # ##########
    def Main(s):
        if s.VB: print("INFO::%s  Main" %(s.myname))
        # for key in s.dict.keys(): print 'dict:  %-10s  %s' %(key, s.dict[key])

        if 'ExampleUsage' in s.cmd:
            print(s.cmd)
            s.ExampleUsage()

        if 'DumpScenario' in s.cmd: s.DumpScenario()


    # ##########
    def DumpScenario(s):  # could actually be method of IsaproScan
        if s.VB: print(50*'#', 'DumpScenario')

        outs = []
        head = ''
        res = {}
        for iScan in range(s.nScans):
            if iScan > 0: continue

            scan = s.Scan[iScan]
            nScen = len(scan.L)
            if s.VB>1: print('dir_prox: %s' %(scan.dir_prox))

            for iScen in range(nScen):

                # ###################
                # ### Init
                scen = scan.L[iScen]  # this is the Isapro object (the scenario)

                ID = scen.ID
                res[ID] = {'mass':{}, 'xsec':{}, 'br':{}, 'wig':{}, 'isaout':{}}  # could also make dynamic
                
                if s.VB>1: 
                    print() 
                    print('\nscen: %i  %s' %(iScen, scen.ID))
                    print(' ', type(scen), scen)  # Isapro instance
                    
                    print('  fn_prox: %s' %(scen.fn_prox))
                    print('  xTot = %.6f' %(scen.xTot))
                    # print '  x:  %s' %(scen.x)
                    print('  ox: %s' %(scen.ox))

                    #print '  wig: %s' %(scen.wig)
                    print('  wig: %s' %(scen.dict['wig']))
                    print('  isaout: %s' %(scen.dict['isaout']))

                #out = ''   # deprecated (I think); instead outt is used
                outt = ''


                


                # ###################
                # ### wig content ('n12',....,'staumix11','mu','a_tau',...)  # a bit cryptic
                #mkeys = ['isaout','wig']
                #mkeys = ['wig']
                # This general structure works so far for 'mass','wig','isaout'
                for mkey in list(s.dump.keys()): 
                    if s.VB>2: print('mkey = %s' %(mkey))
                    if s.VB>2: print('s.dump = %s' %(s.dump))
                    for ivar in range(len(s.dump[mkey])):
                        var = s.dump[mkey][ivar]

                        if var not in scen.dict[mkey]: sys.exit("Fatal::isa_simplereader_class  var %s not among allowed ones for dict %s: %s" %(var, mkey, list(scen.dict[mkey].keys())))

                        if ivar == 0 and s.dict['delim_between_blocks']:
                            outt += ' | '
                            if iScan == 0 and iScen == 0: head += ' | '
                        if iScan == 0 and iScen == 0: head += ' %s ' %(s.dump_headdict[mkey][var])
                        val = scen.dict[mkey][var]
                        if mkey in ['mass']: val = abs(val)
                        #print mkey, s.dumpformat[mkey]
                        outt += ' %s ' %(s.dumpformat[mkey]%(val))
                        res[ID][mkey][var] = val
                    



                # ###################
                # ### Cross-section                
                ox = scen.ox   # libISAWIG::OX-object; has leafs ordered xsec titles: xOT[],  the xsecs: xO{}, (x,xT,xInput)
                for ixsec in range(len(s.dump_xsecs)):
                    #if ixsec == 0: out += '| XSEC[pb]: '
                    #if ixsec > 0 : out += ','
                    if ixsec == 0 and s.dict['delim_between_blocks']:
                        if iScan == 0 and iScen == 0: head += ' | '
                        outt += ' | '
                    
                    xsec = s.dump_xsecs[ixsec]   # ex: 'N2 C1'
                    xsecT = xsec.replace(' ','') # ex: 'N2C1'
                    val = ox.xO.get(xsec,-1)
                    #out += ' %s= %10.6f ' %(xsecT, val)

                    if iScan == 0 and iScen == 0: head += ' %s ' %(xsecT)
                    outt += ' %10.6f ' %(val)

                    res[ID]['xsec'][xsecT] = val

                    
                # ###################
                # ### BRs
                br = scen.br   # libISABR::ParticleDecay. This object has a leaf part which is a dict of libISAWIG::Particle objects, ex. br.part['N2']
                #print 's.dump_brs: %s' %(s.dump_brs)
                for ibr in range(len(s.dump_brs)):
                    #if ibr == 0: out += '| BR: '
                    #if ibr > 0: out += ','
                    if ibr == 0 and s.dict['delim_between_blocks']:
                        if iScan == 0 and iScen == 0: head += ' | '
                        outt += ' | '
                        
                    
                    wbr = s.dump_brs[ibr]
                    #print 'wbr: %s' %(wbr)
                    
                    SUSYparent=wbr[0]
                    

                    SUSYchild = wbr[1]
                    if len(wbr)==3: SMchildren = wbr[2]
                    else: SMchildren = ''
                    val = br.part[SUSYparent].BR(SUSYchild, SMchildren)
                    BRtxt = '%s->%s' %(SUSYparent, SUSYchild)
                    if SMchildren: BRtxt += ',%s' %(SMchildren)
                    
                    #out += ' %s->%s' %(SUSYparent, SUSYchild)
                    #if SMchildren: out += ',%s' %(SMchildren)
                    #out += ' %s' %(BRtxt)
                    #out += '= %6.4f ' %(val)

                    #head += ' %s->%s' %(SUSYparent, SUSYchild)
                    #if SMchildren: head += ',%s' %(SMchildren)
                    if iScan == 0 and iScen == 0: head += ' %s' %(BRtxt)
                    outt += ' %6.4f ' %(val)

                    res[ID]['br'][BRtxt] = val



                # ###################
                # ### DarkSusy
                #ds = scen.ds
                #print ds
                for ids in range(len(s.dump_ds)):
                    var = s.dump_ds[ids]
                    #val = ds[var]
                    val = scen.ds[var]
                    if ids == 0 and s.dict['delim_between_blocks']: 
                        if iScan == 0 and iScen == 0: head += ' | '
                        outt += ' | '
                        outt += ' %6.4f ' %(val)




                # ###################
                # ### Final fixing
                if outt == '': outt = 'No info dumped. Specify with e.g. -dump_masses N1,N2,N3,N4,C1,C2,eR,T1,h,t  -dump_xsecs N2_C1,C1_C1  -dump_brs N2:eR,C1:eR,N2:N1:Z  -dump_brs_all N2,C1'

                #out = '%s  %s |' %(scen.ID, out)

                if s.dict['show_scenarioname']: 
                    outt = '%s  %s |' %(scen.ID, outt)
                    if iScan == 0 and iScen == 0: head = '%s  %s |' %('Scenario', head)
                
                outs.append(outt)
                
                # output
                #print out

                if 0:
                    if iScan == 0 and iScen == 0: print(head)
                    print(outt)


                # VERBOSE
                for ibr in range(len(s.dump_brs_all)):
                    if ibr == 0: print('===== Showing all decays for particles: %s' %(s.dump_brs_all))
                    SUSYparent = s.dump_brs_all[ibr]
                    br.part[SUSYparent].showDecays()
                    print(30*'-')

            #
        #
        outs.insert(0,head)
        outs = lineup.lineup(outs)

        #fn_dumpscenario_pickle = 'DumpScenario.pickle'
        #fn_dumpscenario_txt = 'DumpScenario.txt'
        

        if s.dict['save']: 
            #fntmp='fn_dumpscenario.tmp'
            PickleToFile(fn=s.fn_dumpscenario_pickle, thepickle=res, VB=s.VB-1)
            WriteToFile(fn=s.fn_dumpscenario_txt, outs=outs, VB=s.VB-1)
            #WriteToFile(fn=fntmp, outs=outs, VB=s.VB-1)
            #os.system('cat %s | kilelinetool.py -lineup_header 0 > %s; rm %s' %(fntmp, s.fn_dumpscenario_txt, fntmp))

        if s.dict['show']:
            for out in outs: print(out)
            
        
                
    # #######################################
    # ####################################### ARGUMENT READING
    # #######################################
    # ##########
    def ReadArg(s): 

        Arg = bkgjelstenArgReader.ArgReader(s.argv, VB=0)

        '''
        if Arg.hasget('-alist'):  print 'a string list: ',Arg.list()
        if Arg.hasget('-alisti'): print 'an integer list: ',Arg.listI()
        if Arg.hasget('-alistf'): print 'a float list: ',Arg.listF()
        if Arg.hasget('-x'):  print 'a string: ',Arg.val()
        if Arg.hasget('-xI'): print 'an integer: ',Arg.valI()
        if Arg.hasget('-xF'): print 'a float: ',Arg.valF()
        '''

        if Arg.has(['help','-h','--help','--h','-help']):
            s.showHelp()
            sys.exit()

        if Arg.hasget(['-vb','-VB']):
            s.VB = Arg.valI()
            if s.VB: print('Verbosity level: %i' %(s.VB))

        # new

        if Arg.hasget('-cmd'):
            zz = Arg.list()
            for z in zz:
                if z in s.cmds: s.cmd.append(z)
                else: sys.exit("Fatal: Command '%s' not among known ones: %s" %(z, s.cmds))  

        if Arg.hasget('-dump'):  # Ex: -dump isawig:M_1/M1,M_2/M2,MU::wig:alphah     #mass:N1,N2,N3
            zdumps = Arg.list('::')
            for zdump in zdumps:
                zdict,vars = zdump.split(':')
                s.dump[zdict] = []
                s.dump_headdict[zdict] = {}
                for varandvarT in vars.split(','):
                    w = varandvarT.split('/')
                    var = w[0]
                    if len(w) == 2: varT = w[1]
                    else: varT = w[0]

                    s.dump[zdict].append(var)
                    s.dump_headdict[zdict][var] = varT  # this is just to replace name in the header
            s.cmd.append("DumpScenario")
            #print s.dump
            #print s.dump_headdict
            
                    


        if Arg.hasget('-fn_dump'):
            s.fn_dumpscenario0 = Arg.val()

        #if Arg.hasget('-dump_masses'):
        #    #s.dump_masses = Arg.list()
        #    s.dump['mass'] = Arg.list()
        #    s.dump_headdict['mass'] = 
            
        if Arg.hasget('-dump_ds'):
            s.dump_ds = Arg.list()
            s.cmd.append("DumpScenario")

        if Arg.hasget('-dump_xsecs'):
            s.dump_xsecs = []
            zz = Arg.list()
            for z in zz: 
                s.dump_xsecs.append(z.replace('_',' ').replace(':',' '))
            s.cmd.append("DumpScenario")

        if Arg.hasget('-dump_brs'):
            s.dump_brs = []
            zz = Arg.list()
            for z in zz:
                s.dump_brs.append(z.split(':'))
            s.cmd.append("DumpScenario")
                
                
        if Arg.hasget('-dump_brs_all'):
            s.dump_brs_all = Arg.list()
            s.cmd.append("DumpScenario")

            
        if Arg.has('--clean'): s.dict['clean'] = 1


        if Arg.has('--nosave'): s.dict['save'] = 0
        if Arg.has('--noshow'): s.dict['show'] = 0


        # begin(old)
        # ######################################## ARGUMENT READING

        #if testarg(Argv,'--automan') == 1:
        #    autoMan()
        #    sys.exit()

        if Arg.has('-h'):
            doHELP()
            sys.exit()

        if Arg.has('--man'):
            doHELP()
            sys.exit()

        if Arg.hasget('-testread'): 
            print("testread: arg = ", Arg.val())



        if Arg.hasget('-id'): 
            zz = Arg.list()
            for iz in range(len(zz)):
                z = zz[iz]
                s.fn_IDs.append(z)
                print(iz, z)
                os.system('pwd')
                # f = open(z); s.IDs2.append(f.readlines()); f.close()  # pre 2013-11-23 (later: allows allIDs with directory and ending with .slha ; both will be removed [untested]
                s.IDs2.append([])
                f = open(z); zLs = f.readlines(); f.close()
                for zL in zLs:
                    zL = zL.strip().split('/').pop().replace('.slha','')
                    s.IDs2[-1].append(zL)


        if Arg.hasget('-ending'):
            z = Arg.val()
            if z == '.':
                s.vars['ending_out'] = '.out'
                s.vars['ending_wig'] = '.wig'
            elif z == '_':
                s.vars['ending_out'] = '_out'
                s.vars['ending_wig'] = '_wig'
            else:
                print('FATAL: Non-allowed ending: %s' %(z))
                sys.exit()

        if Arg.hasget('-endingoutwig'): 
            z = Arg.list()
            s.vars['ending_out'] = z[0]+'out'
            s.vars['ending_wig'] = z[1]+'wig'

        if Arg.hasget(['-dir','-lhawig']):  
            zz = Arg.list()
            
            for iz in range(len(zz)):
                z = zz[iz]
                if not z.endswith("/"): z += "/"
                s.dir_isaout.append(z)
                s.dir_wig.append(z)
                s.dir_prox.append(z)

        if Arg.hasget('-prox'): 
            zz = Arg.list()
            s.dir_prox = []
            for iz in range(len(zz)):
                z = zz[iz]
                if not z.endswith("/"): z += "/"
                s.dir_prox.append(z)
            s.vars['USEPROX'] = 1

        if Arg.hasget(['-dir_play','-play']): 
            zz = Arg.list()
            for iz in range(len(zz)):
                z = zz[iz]
                if not z.endswith("/"): z += "/"
                s.dir_isaplay.append(z)
            s.vars['USEPLAY'] = 1

        if Arg.hasget(['-dir_ds','-ds']): 
            zz = Arg.list()
            for iz in range(len(zz)):
                z = zz[iz]
                if not z.endswith("/"): z += "/"
                s.dir_ds.append(z)
            s.vars['USEDS'] = 1

        if Arg.has('--play'):
            s.vars['USEPLAY'] = 1

        if Arg.hasget('-playmode'):
            s.playmode = Arg.val()
            if s.playmode not in s.playmodes: sys.exit('FATAL  isa_simplereader  playmode %s not in allowed ones: %s' %(s.playmode, s.playmodes))
    
        if Arg.has(['--ds','--darksusy']):
            s.vars['USEDS'] = 1

        #if Arg.has('-ds'): 
        #    s.vars['USEDS'] = Arg.valI()

        if Arg.has('--noout'):
            s.vars['USEOUT'] = 0

        if Arg.hasget('-combcode'): 
            s.isaplay_combcode = Arg.val()

        if Arg.has('--withlhaending'):
            s.vars['withLHAending'] = 1

        if Arg.hasget('-inputid'): 
            zz = Arg.list()
            for iz in range(len(zz)):
                z = zz[iz]
                s.inputID.append(z)

        if Arg.hasget('-fn_add'): 
            s.fn_add0 = Arg.val()

        if Arg.hasget('-cutflow'):  # Mar23: to allow (atlfast) cutflow results
            zz = Arg.list()
            for z in zz:
                print("gridpar/cutflow: execfile('%s')" %(z))
                exec(compile(open(z, "rb").read(), z, 'exec'))  # 2013-09: WILL WORK?
        
        if Arg.hasget('-execfile'):  # Mar23: to allow (atlfast) cutflow results
            zz = Arg.list()
            for z in zz:
                print("execfile: execfile('%s')" %(z))
                exec(compile(open(z, "rb").read(), z, 'exec'))  # 2013-09: WILL WORK?
        
        if Arg.hasget('-execgdict'):  # Mar23: to allow (atlfast) cutflow results
            zz = Arg.list()
            for z in zz:
                z2 = 'dir_dicts/' + z
                if os.path.exists(z2):   # NB: dir_dicts/ is searched before ./ 
                    print("execfile: execfile('%s',gDict)" %(z2))
                    exec(compile(open(z2, "rb").read(), z2, 'exec'), gDict)
                elif os.path.exists(z):
                    print("execfile: execfile('%s',gDict)" %(z))
                    exec(compile(open(z, "rb").read(), z, 'exec'), gDict)
                else:
                    sys.exit('FATAL  isa_simplereader: non-existent file %s in ./ and dir_dicts/  (%s)' %(z,os.getcwd())) 
        
        if Arg.hasget('-readervar'):  # Mar23: to allow (atlfast) cutflow results
            s.readervar = Arg.list()

        if Arg.hasget('-freedict'):    # is apparently not in use by howto_std_plots*.sh
            z = Arg.list()
            for z2 in z: 
                z3 = z2.split(':')
                s.freedict[z3[0]] = z3[1]
                print('freedict: ', s.freedict)

        if Arg.hasget('-outer'):
            z = Arg.val()
            if z != 'none':
                s.outer_py = z
                if not os.path.exists(s.outer_py):
                    s.outer_py = "%s/bin/%s" %(s.HOME,s.outer_py)
                    if not os.path.exists(s.outer_py):
                        sys.exit("FATAL: -outer not existing: %s  (also not in ~/bin/)" %(z))
    

        if Arg.hasget('-fntag'):  
            s.fntag = Arg.val()
            if not s.fntag.startswith('_'): s.fntag = '_' + s.fntag

        if Arg.hasget('-fntag2'):  
            s.fntag2 = Arg.val()

        if Arg.hasget('-fnpre'):  
            s.fnpre = Arg.val()

        if Arg.hasget('-fnbase'):  
            s.fnbase = Arg.val()

        if Arg.has('--setbatch'): 
            s.SetBatch = 1

        if Arg.hasget('-plotinfo'):    # is by default put bottom right (e.g. grid type DGnoL, DGemtR, ..) [colourratio, ..?]
            s.plotinfo = Arg.val()

        if Arg.hasget('-txtinfo'):     # a bit unclear where is put
            s.txtinfo = Arg.list()
    
        if Arg.hasget('-fixval'): 
            s.freedict['fixval'] = Arg.valI()
    
        if Arg.hasget('-fixvar'):   # 2013-07-13
            s.freedict['fixvar'] = Arg.val()


        if Arg.hasget(['-v1','-xvar']):   # 2013-07-13
            s.freedict['v1'] = Arg.val()

        if Arg.hasget(['-v2','-yvar']):   # 2013-07-13
            s.freedict['v2'] = Arg.val()

        if Arg.hasget('-v1skip'):   # 2013-07-13
            s.freedict['v1skip'] = Arg.list()
        if Arg.hasget('-v2skip'):   # 2013-07-13
            s.freedict['v2skip'] = Arg.list()

        if Arg.hasget('-lastval'):   # 2013-07-13
            s.freedict['lastval'] = Arg.valI()
        if Arg.hasget('-lastvalv1'):   # 2013-07-13
            s.freedict['lastvalv1'] = Arg.valI()
        if Arg.hasget('-lastvalv2'):   # 2013-07-13
            s.freedict['lastvalv2'] = Arg.valI()


        if Arg.has('--test'):
            s.fnpre = 'test_'

        if Arg.hasget('-steer'): 
            z = Arg.list()
            for zz in z:
                if zz == '': continue # allow empty (comma-trick)
                w = zz.split(':')
                if len(w) != 2: sys.exit("FATAL::isa_simplereader  illegal steer item: %s   (entire is: '%s')" %(zz, z))
                s.steer[w[0]] = w[1]


        if Arg.hasget('-steerlist'): 
            s.steerList = Arg.list()


    
        if Arg.hasget('-steertxtlist'):   # input a text list (or several) into steer{}
            zz = Arg.list(';')
            for z in zz:
                z = z.split(':')
                key = z[0]  # e.g. 'RANGE'
                zcodes = z[1].split(',')
                s.steer[key] = zcodes


        if Arg.hasget('-steerintlist'):   # input an integer list into steer{}
            zz = Arg.list(';')
            for z in zz:
                z = z.split(':')
                key = z[0]  # e.g. 'RANGE'
                zcodes = z[1].split(',')
                s.steer[key] = []
                for zcode in zcodes:
                    if zcode == '': continue  # skipping empty (can happen)
            
                    if zcode.count('-') == 1: # then have e.g. 101-105    # so e.g. -steerintlist RANGE:101,102,104-109,121  is allowed
                        z0,z1 = zcode.split('-')
                        s.steer[key] += list(range(int(z0), int(z1)+1))
                    else:  # enter in usual way, i.e. code per code
                        s.steer[key].append(int(zcode))

        if Arg.hasget('-format'):   # input an integer list
            zz = Arg.list(';')
            for z in zz:
                za = z.split(',')
                s.format[za[0]] = za[1]

    
        if Arg.hasget('-dumpformat'):
            zz = Arg.list()
            for z in zz:
                w = z.split(':')
                mkey = w[0]
                zform = w[1]
                s.dumpformat[mkey] = zform
            s.cmd.append("DumpScenario")


        if Arg.hasget('-figdir'): 
            s.figdir = Arg.val()

        if Arg.hasget('-contours'): 
            s.contours = DefineContours(Arg.val())
    
        if Arg.hasget('-flip'): 
            s.flip = Arg.valI()
            if s.flip not in [0,1]: sys.exit('Fatal::isa_simplereader  non-allowed flip (%i) ... must be 0 or 1' %(s.flip))

        if Arg.hasget('-fn_savehistos'):
            s.fn_savehistos = Arg.val()
            if not s.fn_savehistos.endswith('.root'): fn_savehistos += '.root'


        # end(old)
            

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
            sys.exit("FATAL  Ending due to problems of arguments")
    
        # ################################### POST-INIT


############################## EXECUTE IF RUN AS SCRIPT (NOT JUST IMPORTED)
if __name__ == '__main__':
    t = isa_simplereader(cmd=['ReadArg','PostInit','Main'])
############################## 






                 






    












