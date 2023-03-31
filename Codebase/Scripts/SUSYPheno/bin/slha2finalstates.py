#!/usr/bin/env python

"""
Program : slha2finalstates.py
Version : 1.0
Author  : b.k.gjelsten@fys.uio.no
Description : 




==== TODO

o Rename xsecRel to BRtot (or vice versa)

o Signed 1legs and 2legs ... 
  o Current 2legs are unsigned: need to think this through. It is probably not correct info
  o Option to do fully signed version? : signed 1legs, and by consequence signed 2legs.
    Then would need new methods to optionally join e.g. N2C1+ and N2C1- (since N2 is fully symmetric in + and -)
    -> NEED TO THINK THROUGH SUBPROC BY SUBPROC
  

o Add "tags":
  o Keep track of signs, most notably SS leptons



==== VALIDATIONS
x BRsummation seems to work fine (checked N2->N1cc + N2->N1ss becomes N2->N1->N1qq with correct BR. 
x DecayPaths1Leg: cascStore is filling up as expected
x Small asymmetries between e and m are checked and understood ; due to 2.8% and 1.7% asymm in T and W decays


==== WARNING, ISSUES
o slha files (from SUSYHIT) are bad for b1 and b2: SOME NEGATIVE BRs (susy_DGnoSL_TB10_M1M2MU_050_200_300.slha)
  - should investigate further, then report to authors of SUSYHIT++

o Cannot tackle SUSY daughters in SM decay (h->C1C1)
  o Is N1s ok?
  o Solution?: include SM decays in the ordinary NextDaughter procedure

o Flavour reduction is currently not entirely under control
  o Default reduction is done in GetParticleDict_ala_IsaPlay() via ParticleName(pdgID, optD=s.pname_optD) ... No.. is empty, and there are no off-shell decays : need to check another scenario
  
  o SM reduction in SM decays:
    - is done purely in allSMcombExpand: this is now steered entirely by e.g. "-particles2keep b,Th,l,q,y"  (GOOD!)
  o SM reduction in SUSY decays: 
    - 

===== PERFORMANCE, SPEED
- Full 30particles in xsec: takes 1m20s for BRmin1leg,1e-3:F,BRmin2leg,1e-3
- The bulk usage is in 1leg::ExpandSMp and Duplicates, and for the scenario, 050_100_300, gl, t2 and t1 took nearly all the time. Don't know why




===== SHARP (Very tuned approach)
  x sort masses: let smallest go first,
  x include also SM
  o Separate conversion tables (1:1 conversions, e.g. e->l .. and there is Th..) : maybe do these even in particle array
  o 1-leg by:
    / Datastructure: one leg: a list of daughters (including SUSY), 1leg['~N2'] = [[N1,v,v],[N1,l,l]] and a dict of desired (dynamical) classifiers
    / if (first) SUSY daughter is known [it will be]: insert the list directly to the larger list, maybe multiply in smart way with an already developed SM-part ; that way a later SM expansion will not be needed




===== SHARP
x performance: early test shows <900 GeV 1leg at 0.5s vs 1.7s for legacy#2
o BUGS: insurmountable problems with Hc->N4N4 and earlier on very slow with Hc->N2N3 : Need to understand... why??


o Make easy-to-use outputs


Eventually:
o Mass to b and the other quarks (decay table)
o dump: intermediate step where particle decay is itself expanded: a) particle conversion only (not h,Z,W,..) ; then also with those ; this could be a precursor to GetOnelegs()


Structure, naming, object types: need to be standardised

    varnames: xsecRel vs BRtot .. etc.: need standards
    
    

2014-01-14: To do:

o Remake Onelegs as list of BBCascades : hmm, that's really tough ; it is really specialised/complex
  Maybe just do at the end of GetOnelegs? 

o Output: 
  x lists 
  o easily readable dicts (&pickle)



FUTURE FEATURES/IMPROVEMENTS
x FS: ability to read existing pickle (well, reading is trivial, that is just the pickle)
x FS: ability to treat further, e.g. joining channels (skipping final particles)

o Develop superclass which contains a list of Mother (typically a scenario) ? 
  - What would be relevant features?
  o Feature: object.finalstate({'l':1, 'q':2}) should give output like rel={'C1+:C1-':0.713, 'C1+:N2':0.22}
  o Feature: generalis & join subprocs: Ex. generic: let C1+,C1- -> C1 or let N2,N3 -> N23, ... 
  

o Develop constructive (non-debugging) verbosity on program ... e.g. is Expand used only once as one should expect
o Develop time estimate: where is time spent
o Unspecialise code: Check if works if: no xsecs, ...

o Opt: Introduce a second particle dict where all 1:1 exchanges [as given as default/in input; ex: u,d,s,c->q, ve,vm,vT->v, ..] are already done

o Opt: use libExpValues for SM BRs 

o FlatDict


------------------------------------------------------------------------------------------
STILL FRAGILE & CONFUSING  <-- CHECK THE FOLLOWING IF THERE IS A NON-UNDERSTANDABLE ERROR

o The use/nonuse of the separate width attribute of Mother. See e.g. how using dropwidth=1 fixed things for thisMother_Onelegs.AddMother


---------------------
o 2-liners
  o masses (1-2file)
  o other quantities (bino part, mixing angles) (1 file)
  o widths & BRs <-- several files,  (several files)






HOW TO MAKE XSECS_TEMPLATE FOR INTERPOLATION USING STANDARD UNCERTAINTY TABLES (OR EQUIVALENT):  2014-02-11
  (more rigorous & dynamic solution to be developed ... in relation to SUSYBANK / XSECBANK)

cd /home/scratch/borgeg/lsp/alt/SignalUncertaintiesUtils/SignalUncertaintiesUtils_2013-09-06/res_DGnoSL_TB10_M1eq50

for sub in ${subs} ; do less Herwigpp_UEEE3_CTEQ6L1_DGnoSL_TB10_X_Y_Z_2L.txt | sed s/' MU '/' mu '/ | grep ' XS \| '${sub}'  |' | TableToTGraph.py  -coords mu,M2 -resvars XS:XS_${sub}  -save dir_template_XS/template_mu_M2__XS_${sub}  -dict I,fnaddaxes,0; done 

python 
import ROOT, os
from kilelib import WritePickle, LoadPickle
fn_pickle="DGnoSL_NLO_3E3__mu_M2.pickle"


# LOAD FIRST TIME (afterwards can just use pickle)
dir_fxsec = "dir_template_XS"
fn_fxsec_start = "template_mu_M2__"
fxsec = {}
fns = os.popen("ls %s/%s*.root" %(dir_fxsec, fn_fxsec_start)).readlines()
for fn in fns: fn = fn.strip(); ff = ROOT.TFile(fn); var = fn.split('/').pop().replace(fn_fxsec_start,"").replace(".root",""); var2 = int(var.replace('XS_','')); z = ff.Get(var); z.SetDirectory(0); fxsec[var2] = z ; ff.Close()




# 2014-08-27 : COMMENTS

- Crucial parameters:

x SM particles (h,Z,W,..) are not "JoinIdenticalled", it seems, while this is fine for SUSY particles ... FIXED


CLASS TYPES: 
 A  class FS:  [could be called FS,fs,finalstate,decay] 
       br         Ex: br=0.34
       children   Ex: children = {'~N1':1, 'l':2, 'b':2, 'Th':1}
     Is a small class: Is just one decay channel, has no info on the parent (except br)
      
                          
 B  class Mother: [could be called FSlist or PartFSs, FSs]
       name = particle name
       mass = mass
       width = total width 
       DecaList: list of FS objects
     Has many functions:
       ShowDecay()
       DevelopFinalstates()  <-- fat
       ExpandDecay()         <-- fat




OBJECT TYPES: 
  particles : dict of FSlist. But it does not comply with particles2keep, is given in e,m,T,...

  with initiators = ['~N2','~eR'] (and a 'x' autoput in front), we build Z,h,~N1 while building ~N2
  

BUG 1: FIXED (FIX1a-c in SlhaTools)
 - 1leg RELacc sums to above 1 for h, N2, ... while C1 and W are ok (for this scenario)
 - if initiators = ~N1,~C1 we have the problem
 - if initiators = ~N1 they sum to 1 .. in this case ~C1 (obviously) but also W are absent from 1leg file
     Why is that? Shouldn't W be generated e.g. in relation with decay of h? How is that taken care of?
 x FIXED: all related to the usage (&change) of the particles dict in ExpandDecay()     



NEXTs: 
 x Document 
 / Set up nice examples & defaults
 / Allow for silent running 
 o Allow for importing and running from other scripts
 o how to input for plotting? Can probably already use the flatdict
 x Allow some printing of tlts? 
 x Add more info to flatdict (in order to plot as function of e.g. M1,M2,...): import slha0 with slha2tlt imports ? 


"""

import sys,os
import bkgjelstenArgReader 
from kilelib import WriteToFile, WritePickle, LoadPickle, ExtractVarsViaFormat, SaveHistory, IterateNdim, GeneaLog, dict2tlt


import pyslha_edited as pyslha
from SlhaTools import ParticleName, pname_generic2signed, PDG, namGeneric

from SlhaTools import Mother, FS, AddSMtoMothers, SumMothers
from SlhaTools import slha2flatdict  # to read slha[0]

from ProspinoLib import func_prospino_txtpair2code, prospino_code2txtpair

from MyTimeTools import MyTimer



# ##################################################### GLOBAL METHODS
# ##################################################### GLOBAL METHODS
# ##################################################### GLOBAL METHODS

def isLSP(name):
    return (name == '~N1' or name == 'N1')
def isSUSY(name): return (name.startswith('~'))


def add2dict(a,b):
    # http://stackoverflow.com/questions/1031199/adding-dictionaries-in-python
    return dict( (n, a.get(n, 0)+b.get(n, 0)) for n in set(a)|set(b) )





# ##################################################### GLOBAL METHODS


class slha2finalstates: 

    def __init__(s, cmd=[], optD={}):

        # ====================== PRE INIT
        if 'argv' in optD: s.argv = optD['argv']
        else: s.argv = sys.argv

        # The genea structure is structure output (to help see what goes on in the program)
        # output to screen is then preferably done like:
        #  s.geneaLog.Add(genea=genea, S='fatal', txt=zzz)  ; with S in ['info','warning','error','fatal']
        #  And at the beginning of each routine, we do genea = s.geneaLog.Genea(genea) - which updates the genea string 
        genea = optD.get('genea','')
        s.geneainit = optD.get('geneainit','')
        s.geneaadd = optD.get('geneaadd','')
        s.geneaLog = GeneaLog(optD={'maxdepth':5, 'maxwidth1':12})
        

        s.cmd = cmd
        s.myname = sys.argv[0].split('/').pop()
        s.VB = 1
        s.HOME = os.getenv('HOME')
        s.cwd  = os.getcwd()  # current work directory  
        #s.dir0 = '%s/XXX' %(s.HOME)
        s.dir0 = '.'
        s.dir1 = '%s' %(s.dir0)

        s.dict = {}

        s.warn = []
        s.fn_warn = 'warnings_slha2finalstates.txt'
        s.fn_report = 'report'
        s.report = []


        # Dicts and Lists for scenario variables
        s.fns = []
        s.slha = {}
        #s.particles = {}  # particles are not global, they belong to a scenario, hence no global variable

        #s.pname_optD = {}
        s.generic = []
        #s.generic ['q','v']
        #s.dict['dumpBR_which'] = []  # empty dumps all (if dump is on)

        s.dict['dumpBR'] = 0
        s.dict['1leg_minmass'] = 0.
        s.dict['1leg_maxmass'] = 900.  # will then skip "infinite"-mass sparticles

        s.dict['dump_minmass'] = 0.
        s.dict['dump_maxmass'] = 900.  # will then skip "infinite"-mass sparticles

        s.dict['dump_nmaxchildren'] = 4
        s.dict['dump_mode'] = 'flat'

        s.DEISM = {}  # DressedEventwithInitiatorsandSM .. this is the final pheno object for the standard full cascade expansion and closing

        s.dict['BRmin1leg'] = 1e-5
        s.dict['BRmin2leg'] = s.dict['BRmin1leg']
        s.dict['BRmin2legFinal'] = 1e-5

        s.dict['BRmaxdev1leg'] = 1e-4  # will give warning if exceeds this
        s.dict['BRmaxdev2leg'] = 1e-3  # will give warning if exceeds this

        s.dict['runmode'] = 1  #  [0:only particle dump?]  1:1legOnly; 2:2legToo  [D=-1]

        s.dict['1legstofile'] = 1

        s.SUSYinitiatorsBrute = []
        s.SUSYinitiators1leg = []
        #s.SUSYinitiators2leg = []  #<-- this needs just be a local variable
        s.subprocs2legBrute = []
        s.subprocs2leg = []
        

        s.particles2keep = ['l','Th','b','q','y','v','~N1']

        s.fn_xsecs = ''
        s.xsecs_dictname = 'X'

        s.fn_xsecs_template = ''
        s.xsecs_template_x = -1
        s.xsecs_template_y = -1


        s.slhaformat = '*'

        s.timer = MyTimer()

        s.dict['BRminDecay'] = 1e-3  # cut on the BR straight from particles dict (slha)
        
        s.fntag = ''

        s.flatdict = {}
        s.flatdict_slha0 = {}
        s.dict['flatdict'] = ['onelegs','twolegs']

        s.dict['flatdict_gteqlt'] = ['>=','=','<=']
        s.dict['flatdict_npartMaxCommon'] = 2
        s.dict['flatdict_npartMax'] = {'L':3,'l':3,'Th':3}   # (pretty arbitrary default)

        s.dict['flatdict_slha0'] = 1  # 1:Put slha[0] part in flatdict too ; 0:no

        s.tlts = []
        s.dict['tlt_delim_name'] = ':::'
        s.dict['tlt_format0'] = '%.4f'  #  (small because of BRs...)
        s.dict['tlt_delim_var'] = ','
        s.dict['tlt_delim_rename'] = '::'

        # THIS IS COPY FROM slha2tlt [inelegant...]
        s.formats = {}
        for z in ['MASS','EXTPAR','mu','mu(Q)','~q4mean','q4mean','GMEANstop','GMEANsbot','Q']: s.formats[z] = '%.1f'
        for z in ['CALC','CALCver','FULLCALC']: s.formats[z] = '%s'
        for z in ['ICALC']: s.formats[z] = '%i'
        s.formats['TB(MX)'] = '%.3f'  # hack
        s.dict['tlt_formats'] = s.formats


        s.dict['fn0'] = 'fs_'

        s.fsinfilename = ''   # this is to add the finalstates to the filenames, e.g. add "_L_h" or "_l_Th_b"


        s.fn_history = 'history_slha2finalstates.txt'
        s.fn_history_global = '%s/.history_slha2finalstates.txt' %(s.HOME)
        SaveHistory(fn=s.fn_history, argv=sys.argv, opt=['stripcmd','date'])
        SaveHistory(fn=s.fn_history_global, argv=sys.argv, opt=['stripcmd','date','dirsameline'])
        

        # ====================== READ ARG
        if 'ReadArg' in s.cmd: s.ReadArg(genea=genea)


        # ====================== POST INIT
        if 'PostInit' in s.cmd: s.PostInit(genea=genea)


        # ====================== EXECUTE 
        if 'Main' in s.cmd: s.Main(genea=genea)

        s.geneaLog.SaveToFile()  # 2014-08-28




    # ##########
    def PostInit(s, genea): 
        genea = s.geneaLog.Genea(genea)
        #s.fn_warn = '%s/%s' %(s.dir1, 'warnings.txt')
        #s.fn_report = '%s/%s' %(s.dir1, 'report')

        s.pname_optD = {'generic':s.generic}

        if s.fsinfilename:
            s.fsinfilename = ''
            for p2keep in s.particles2keep: s.fsinfilename += '_%s' %(p2keep)

        # Pipe in slha files? 
        if not s.fns: 
            if s.VB > 0: print("INFO::slha2tlt  getting filelist from stdin... (otherwise need to specify '-f <file_with_slhalist>'")
            s.fns = sys.stdin.readlines()

        


    # ##################################################### CLASS METHODS
    # ##################################################### CLASS METHODS
    # ##################################################### CLASS METHODS
    # ##########
    def showHelp(s):
        print('DESCRIPTION')
        print("  slha2finalstates.py takes as input a list of slha files and the list of desired finalstate particles (e.g. l,b ; leptons and b-quarks)")
        print("  then for each slha does (some of) the following for all SUSY particles (or the ones selected)")
        print("   0  Reads the branching ratio (BR) table from slha (+ adds some SM decays): ")
        print("   1  Produces so-called 1legs")
        print("        Here the decays of (s)particles are rewritten in terms of the selected finalstate particles (e.g. l,b)")
        print("        Ex: For h there would be 'empty' decay in ~24% of cases, h->2b in ~63% of cases")
        print("        Rewriting the BR tables in terms of finalstate particles can be practical when scrutinizing search channels")
        print("        (Effects like kinematic cuts, experimental identification, resolution etc. are not included, this is the pure theoretical BR-table)")
        print("   2  Produces so-called 2legs (optionally)")
        print("        Here all (selected) 1leg+1leg combinations are made to produce 2leg BR tables.")
        print("        Since SUSY will typically produce 2 legs in the detector, the 2leg BR tables are more relevant than the 1leg tables")
        print("        However, to get fully useful numbers out, we also need to know and use the subprocess cross-sections")
        print("        If these are input, the 2leg BR tables will also contain the cross-section ")
        print("        and a relative importance (cross-section compared to the full SUSY cross-section,") 
        print("        as well as a (most useful) summary 2leg tables which contains the sum of the 2legs.") 
        print("        This summary 2leg table contains the most useful information.") 
        print()    
        print("   3  Produces also a flat dictionary which contains all the results (e.g. to use for plotting with munch.py)")
        print("      Optionally: produces tlts (two-line-tables) of a (small) selection of the variables (see below)")
        print() 
        print()
        print("  (S)PARTICLE NAMING CONVENTIONS:")
        print("    Try 'slha2finalstates --showparticles' to see the hierarchy of (s)particle names")
        print("    (Note: T=tau, v=neutrino ; the rest should be self-explanatory)")
        print("    (NB: sneutrinos are missing.. CHECK)")
        print("    As seen there are basic names (e,d,~eR,..) and more general ones (l, L)")
        print("    In addition comes Th which is a hadronic tau. If Th is specified in particles2keep, all T expanded into Th (and possibly leptons and neutrinos if selected)")
        print("    And there is also W,Z,h,t.")
        print("    (Note: The general names like l,q can only be used for SM particles (usually in particles2keep), not for SUSY particles.)")
        print() 
        print()
        print("  COMMENTS")
        print("    Currently no bugs are known, but could do with some more validation.")
        print("    Is pretty slow when there are many initiators.")
        print("    The algorithm is a bit on the complicated side (to gain speed, which probably didn't happen).")
        print() 
        print() 
        print("OPTIONS  (see '--h' for full list (though without explanations)")
        print("   -f <file_with_list>  : the slha files. Or they can be piped in, a-la 'ls *slha | slha2finalstates.py [options]")
        print("   -particles2keep <commalist> : these are the finalstate particles you want to study the scenario(s) in terms of")
        print("                                  D=l,Th,b,q,y,v,~N1 , ex: l,q (in terms of leptons and quarks), or l,q,h")
        print("                                  See above for info on the (s)particle naming conventions")
        print("   -initiators <commalist>     : [D=empty] The SUSY initiators you want for the 1legs [and 2legs]")
        print("                                 Ex:  ~N2,~C1,~N3  or  ~N2,~C1+  or  ~N2,~eR,~C1-  or ...")
        print(" ")
        print("   -subprocs <commalist>       : [D=empty] List of subprocesses to use for 2legs")
        print("                                 Ex:  ~N2:~C1,~N2:~N1,~N2:~N2,~N2:~N3  or ~N2:~C1+,~C1:~C1  or ...")
        print("                                 If -subprocs are specified, they override the -initiators (which are then superfluous)")
        print("                                 If only -initiators are specified, all their combinations will make up the subprocs")
        print("   -runmode <0/1/2>            : [D=%i]  0:BRdump only ; 1:BRdump,1legs ; 2:BRdump,1legs,2legs" %(s.dict['runmode']))
        print("   -dumpBR  <0/1/2/3>          : [D=%i]  0:none ; 1:to screen ; 2:to file ; 3:to screen and file" %(s.dict['dumpBR']))
        print("   -fsinfilename <0/1>         : [D=%i]  Add the finalstates specification (e.g. 'l_b') to the filenames")
        print("                                 (output filenames are based on the slha-filenames)")
        print("   -xsecs <xsecsfile>             : Attach a xsecfile. [DETAILS ON FORMAT]")
        print("   -xses_template <templatefile>  : Attach a xsec template file:  interpolate to get xsecs [DETAILS ON FORMAT]")
        print() 
        print("   -flatdict_gteqlt <ComparisonConditions>   : Details the quantities to be stored in the output flatdictionary")
        print("        D='>=,=,<=' which means that quantities for e.g. h ->  >=1b + =1l is stored. If only '=' was included, only quantities like h -> =1b + =1l  or h -> =0b + 1l  etc. are stored")
        print("        Ex:  -gteqlt '=' , -gteqlt '>=' , -gteqlt '>=,=' , etc.")
        print("        Note that the arguments will typically need to be protected with ''")
        print("   -flatdict_npartMaxCommon <int>   : [D=%i]  Is the max number of particles of a given type for which to show results (common value)" %(s.dict['flatdict_npartMaxCommon']))
        print("   -flatdict_npartMax <dictlist>    : [D=%s]  Is the max number of particles of a given type for which to show results (specific for given particle type)" %(s.dict['flatdict_npartMax']))
        print() 
        print("   -tlt <commalist_of_vars:::name>  :")
        print("         ex: -tlt M1,M2,mu,BR:h:=2l_=0b,BR:h:=2l_=1b:::tableBR1                   # (BR:h:=2l_=0b is here a variable)")
        print("         ex: -tlt M1,M2,mu,BR:h:=2l_=0b::h_2l,BR:h:=2l_=1b::h_2l_1b:::tableBR1    # Can rename a var with '::'")
        print("         ex: -tlt 'M1,M2,mu,BR:h:=2l_=0b::varx,BR:h:=1l_>=0b:::table1'            # If a var contain e.g. '>', you'll need to protect the tlt definition with ''")
        print("       (See dict vars 'tlt_*' for variables which affect the tlts)  (try '--dict')")

        print() 
        print("  A number of (dict) parameters are accessible via a generic input system: ")
        print("    Ex: -dict I,var1,3:F,var2,3.4:var3,sometext:I,var4,44")
        print("        Variables are separated by colon; a given variable is specified by type,varname,value")
        print("        Above var1 and var4 are integers, var2 is float and var3 is a string.")
        print() 
        print("    dict-parameters:")
        print("      1leg_minmass    [D=%.0f] : particles below this are not used to initiate 1legs" %(s.dict['1leg_minmass']))
        print("      1leg_maxmass    [D=%.0f] : particles above this are not used to initiate 1legs" %(s.dict['1leg_maxmass']))
        print("      dump_minmass    [D=%.0f] : particles below this are not shown in dump" %(s.dict['dump_minmass']))
        print("      dump_maxmass    [D=%.0f] : particles above this are not shown in dump" %(s.dict['dump_maxmass']))
        print("      BRmin1leg       [D=%.6f] : 1legs with BR below this limit are dropped from 1leg lists" %(s.dict['BRmin1leg']))
        print("      BRmin2leg       [D=%.6f] : 2legs with BR below this limit are dropped from 2leg lists" %(s.dict['BRmin1leg']))
        print("      BRmin2legFinal  [D=%.6f] : 2legs with BR below this limit are dropped from 2leg summary list" %(s.dict['BRmin2legFinal']))
        print("      BRmaxdev1leg    [D=%.6f] : if 1-sum(BRs) for a given 1leg is less than BRmaxdev1leg, a warning is given" %(s.dict['BRmaxdev1leg']))
        print("      BRmaxdev2leg    [D=%.6f] : if 1-sum(BRs) for a given 2leg is less than BRmaxdev1leg, a warning is given" %(s.dict['BRmaxdev2leg']))
        print("      1legstofile     [D=%i]  0:no ; 1:yes" %(s.dict['1legstofile']))
        print("      flatdict        [D=['onelegs','twolegs']]  Specify which of onelegs/twolegs are to be included in the flatdict")

        
        print() 
        print("   --silent          : No output to screen")
        print("   -geneaVB <'dump'/'debug','info','warn'>  : D='info'  Verbosity level")
        print()
        print("   -vb <0/1/..>      : change verbosity (D=1)")
        print("   -h                : this help message")
        print("   --h               : brute info on implemented options")
        print("   --dict            : dump various default values (which can be changed)")
        print() 
        print()
        print('EXAMPLES')
        #print '        %s  -dict test,txt1:I,testI,3:test2,txt2a,txt2b:F,testF,4.14   # for using autodict (NB: vars need to be defined in __init__)' %(s.myname)
        #print "    Ex: slha2finalstates.py -f susy_DGnoSL_TB10_M1M2MU_050_200_300.slha  -dict I,dumpBR,1  -vb 2   # 2013-12-01 developing, debugging"
        #print "        slha2finalstates.py -f susy_DGnoSL_TB10M1M2MU_050_200_300.slha  -dict I,dumpBR,1  -vb 2  -particles2keep b,T,l,q,y,v  -vb 2"
        #print  "    Ex: slha2finalstates.py -f susy_DGnoSL_TB10_M1M2MU_050_200_300.slha  -dict I,dumpBR,1  -initiators ~N1,~N2,~N3,~N4,~C1,~C2  -vb 2  -particles2keep b,Th,l,q,y  -vb 2"
        #print  "        slha2finalstates.py -f susy_DGnoSL_TB10_M1M2MU_050_200_300.slha  -dict I,dumpBR,1  -initiators ~N1,~N2,~N3,~N4,~C1,~C2  -vb 2  -particles2keep b,Th,l,q,y  -slhaformat 'susy_DGnoSL_TB10_M1M2MU_*_*_*.slha'  -vb 2"
        #print  "        slha2finalstates.py -f susy_DGnoSL_TB10_M1M2MU_050_200_300.slha  -dict I,dumpBR,1  -subprocs ~N2:~C1+  -vb 2  -particles2keep b,Th,l,q,y  -slhaformat 'susy_DGnoSL_TB10_M1M2MU_*_*_*.slha'  -xsecs xsecsnew.py  -vb 4"
        #print  "        slha2finalstates.py -f susy_DGnoSL_TB10_M1M2MU_050_100_300.slha  -slhaformat 'susy_DGnoSL_TB10_M1M2MU_*_*_*.slha'  -xsecs xsecsnew.py  -particles2keep y,l,Th,b,q  -dict I,dumpBR,3:F,BRmin1leg,1e-3:F,BRmin2leg,1e-3  -subprocs ~N2:~C1+"

        #print  "        slha2finalstates.py -f susy_DGnoSL_TB10_M1M2MU_050_100_300.slha  -slhaformat 'susy_DGnoSL_TB10_M1M2MU_*_*_*.slha'  -xsecs xsecsnew.py  -particles2keep y,l,Th,b,q  -dict I,runmode,1:I,dumpBR,0:F,BRmin1leg,1e-3:F,BRmin2leg,1e-3  -subprocs ~N2:~C1+    # 2013-12-13  testing, developing [cd ~/grids_lsp/DGnoSL_TB10_Fine_2013-11-22/lhawig_M1eq050"

        #print  "        slha2finalstates.py -f susy_DGnoSL_TB10_M1M2MU_050_100_300.slha  -slhaformat 'susy_DGnoSL_TB10_M1M2MU_*_*_*.slha'  -xsecs xsecsnew.py  -particles2keep ~N1,y,l,Th,b,q  -dict I,runmode,1:I,dumpBR,2:F,BRmin1leg,1e-4:F,BRmin2leg,1e-4  -vb 2    # 2014-01-20"
        #print  "        slha2finalstates.py  -f susy_DGnoSL_TB10_M1M2MU_050_177_400.slha  -slhaformat 'susy_DGnoSL_TB10_M1M2MU_*_*_*.slha'  -xsecs_template DGnoSL_NLO_3E3__mu_M2.pickle -dict I,runmode,1:I,dumpBR,2:F,BRminDecay,1e-5:F,BRmin1leg,1e-5:F,BRmin2leg,1e-5:F,BRmaxdev1leg,1e-5:F,BRmaxdev2leg,1e-3:fn0,  -vb  -1  -subprocs all  -particles2keep L,b  --fsinfile  # 2014-02-11"
        #print  "          (Note: can interpolate xsecs with -xsecs_template using dict of TGraph2Ds"
        #print
        #print  " 2014-08-28: (the above ones are old and need to be checked)"
        print(" ls *.slha | slha2finalstates.py  -particles2keep l,b  -initiators ~N2  -runmode 1  -dict I,dumpBR,2 -geneaVB debug                       # 1leg")
        print(" ls *.slha | slha2finalstates.py  -particles2keep l,Th,Q  -initiators ~N2  -runmode 2  -dict I,dumpBR,2 -geneaVB debug  -subprocs ~N2:~C1+   # 2leg")
        print(" ls *.slha | slha2finalstates.py  -particles2keep L,Q  -initiators ~N2  -runmode 2  -dict I,dumpBR,2 -geneaVB debug  -subprocs ~N2:~C1+  -xsecs xsecsnew.py # 2leg")
        print(" ls *.slha | slha2finalstates.py  -particles2keep L,q,b  -initiators ~N2  -runmode 2  -dict I,dumpBR,2 -geneaVB debug  -subprocs ~N2:~C1+   -xsecs_template DGnoSL_NLO_3E3__mu_M2.pickle ")
 
        print() 


    # ##########
    def DumpWarnings(s, genea):
        genea = s.geneaLog.Genea(genea)
        f = open(s.fn_warn,'w')
        for out in s.warn: f.write('%s\n' %(out))
        f.close()
        

    # ########## SHARP
    def GetMotherDict(s, genea, slha1, optD={}):
        genea = s.geneaLog.Genea(genea)
        particles = {}

        s.geneaLog.Add(genea=genea, S='info', txt='now starting')
        s.geneaLog.Add(genea=genea, S='dump', txt=' len(slha1)=%i   keys=%s' %(len(slha1), list(slha1.keys())))
        # Loop over particles in branching ratio dict
        for pdgID in slha1:

            #if pdgID in [6]: continue   # hack, ... there is the problem that top does not come with mass. Need to go to other block

            # 1 First get particle info from slha object
            pname = ParticleName(pdgID, optD=s.pname_optD)

            s.geneaLog.Add(genea=genea, S='dump', txt='DEBUG::GetParticleDict  %s' %(pname))

            # 2 Drop similar ones    FRAGILE: must be used with care... hm..
            if pname in particles:
                s.geneaLog.Add(genea=genea, S='debug', txt='INFO::GetMotherDict  skipping %s since already in particles' %(pname))
                continue 


            # 3a Create the particle object
            zmass = slha1[pdgID].mass
            if zmass is None:
                zmass = optD.get('masses',{}).get(pname,None)   # Need this for the top mass (not so important, though)
                s.geneaLog.Add(genea=genea, S='debug', txt='INFO::GetMotherDict  for %s the mass needed to be taken from slhaInfo (%.3f)' %(pname, zmass))
                
            s.geneaLog.Add(genea=genea, S='debug', txt='DEBUG::GetMotherDict  %-4s  mass: %.3f' %(pname, zmass))
            p = particles[pname] = Mother(pname=pname, mass=zmass, width=slha1[pdgID].totalwidth, particles2keep=s.particles2keep, geneaLog=s.geneaLog, genea=genea)



            # 3b Then fill decays into the particle object [incl. summing of decay products]
            for slhadecay in slha1[pdgID].decays:
                # First get sorted children list
                childrenlist = []
                for zid in slhadecay.ids: childrenlist.append(ParticleName(abs(zid), optD=s.pname_optD))
                # Then add to particle (AdDecay adds to existing fs or adds new)
                p.AddFS_brute(childrenlist, slhadecay.br)  # this 

            # 3c Sort
            p.SortDecaysWithBR()


        return particles



    # ##########
    def GetInfo(s, genea, slha):
        """
        Small routine (for now) to get some vitals. Might expand later. 
        """
        genea = s.geneaLog.Genea(genea)
        res = {}
        res['masses'] = {}
        res['masses']['W'] = slha[0]['MASS'].entries[24]
        res['masses']['b'] = slha[0]['MASS'].entries[5]
        res['masses']['Z'] = slha[0]['SMINPUTS'].entries[4]
        res['masses']['t'] = slha[0]['SMINPUTS'].entries[6] # top
        res['masses']['T'] = slha[0]['SMINPUTS'].entries[7] # tau
        #res['masses']['W'] = slha[0]['SMINPUTS'].entries[5]   

        return res
    
    

    # ##########
    def GetInitiators(s, genea):
        """
        Output: 
        zsubprocs2leg <- used to develop 2legs
        zSUSYinitiators1leg <- used to develop 1legs

        """            

        genea = s.geneaLog.Genea(genea)
        zSUSYinitiators1leg = []
        zsubprocs2leg = []

        # Either specified by the subprocs (signed or not) ; this is in principle economic, non-cartesian
        if s.subprocs2legBrute:  # Ex -subprocs N1:C1+,N2:C1- or also -subprocs N1:C1 (will expand C1)
            for subproc in s.subprocs2legBrute:
                p1s = pname_generic2signed.get(subproc[0], [subproc[0]])
                p2s = pname_generic2signed.get(subproc[1], [subproc[1]])
                for p1 in p1s:
                    for p2 in p2s:
                        orderedtuple = tuple(sorted([p1,p2]))
                        if orderedtuple not in zsubprocs2leg: zsubprocs2leg.append(orderedtuple)  # fill the pairs

                # Unsigned ones
                for p in subproc:
                    pgeneric = p.rstrip('+-_')
                    if pgeneric not in zSUSYinitiators1leg: zSUSYinitiators1leg.append(pgeneric)


        # Or specify by initiators (signed or not) ; this is cartesian
        elif s.SUSYinitiatorsBrute:  #ex -SUSYinitiators N1,N2,C1 [will expand C1]   or N1,N2,C1+  [C1- will not be touched]
            zSUSYinitiators2leg = []
            # 1) Initiators from command line
            for z in s.SUSYinitiatorsBrute:
                # 1leg (unsigned)
                #ps = pname_generic2signed.get(z,[z])
                zgeneric = z.rstrip('+-_') 
                if zgeneric not in zSUSYinitiators1leg: zSUSYinitiators1leg += [zgeneric]
                # 2leg (signed) [tmp]
                zSUSYinitiators2leg += pname_generic2signed.get(z,[z])

            # then fill the pairs
            for ip1 in range(len(zSUSYinitiators2leg)):
                for ip2 in range(ip1,len(zSUSYinitiators2leg)):
                    
                    orderedtuple = tuple(sorted([zSUSYinitiators2leg[ip1], zSUSYinitiators2leg[ip2]]))
                    if orderedtuple not in zsubprocs2leg: zsubprocs2leg.append(orderedtuple)  #

        # Note: inside the actual loop, there will be an xsec test (and maybe other tests) on top of zsubprocs2leg
        
        # If subprocs2leg remains empty after this, all the subprocesses of the xsecs will be done. 
        # If the xsecs are empty, nothing will be developed
        #  Hm: Will allow a few more methods: allow sparticles within a given mass range, .. 
        
        return zSUSYinitiators1leg, zsubprocs2leg


    # ########## SHARP
    def MotherDump(s, genea, fn, particles, which=[], optD={}, fnadd=''):
        genea = s.geneaLog.Genea(genea)
        s.geneaLog.Add(genea=genea, S='info', txt='now starting')
        #if s.VB>1: print 'INFO::MotherDump starting'
        outs = []

        if which == []: whichpart = sorted(list(particles.keys()), key=lambda x: abs(particles[x].mass))
        else: whichpart = list(which)  # probably need to do this, otherwise will change the input list

        for iPart in range(len(whichpart)):
            pname = whichpart[iPart]
            p = particles[pname]
            if not ( s.dict['dump_minmass'] <= abs(p.mass) <= s.dict['dump_maxmass'] ) :
                continue
            outs.append('')
            out = '%-3s  %7.1f  %.2e' %(p.name, abs(p.mass), p.width)
            outs.append(out)

            sumBR = 0.

            for decay in p.FSlist:
                sumBR += decay.br
                out = '  BR   %.5f  %8.5f  %9.2e  %-3s  -> ' %(sumBR, decay.br, decay.br*p.width, pname)
                summassdaughters = 0.

                # would actually make sense to implement these as methods of the BlackBodyDecay class

                if s.dict['dump_mode'] in ['flat','old','legacy']:
                    ChildrenList = decay.ChildrenList()
                    for zname in ChildrenList: 
                        out += ' %-4s' %(zname)
                        if zname in particles: summassdaughters += abs(particles[zname].mass)
                        else: summassdaughters += optD.get('masses',{}).get(zname,0.)
                    for i in range(s.dict['dump_nmaxchildren']-len(ChildrenList)): out += ' %-4s' %('-') 
                    
                else:  
                    for zname in decay.children:
                        if decay.children[zname] == 1: out += ' %4s' %(zname)
                        elif decay.children[zname] > 1: out += ' %i%-3s' %(decay.children[zname], zname)
                        if zname in particles: summassdaughters += abs(particles[zname].mass) * decay.children[zname]
                        else: summassdaughters += optD.get('masses',{}).get(zname,0.) * decay.children[zname]
                    for i in range(s.dict['dump_nmaxchildren']-len(decay.children)): out += ' %4s' %('-') 

                
                out += '  %7.1f' %(abs(p.mass) - summassdaughters)
                #print out
                outs.append(out)

        outs.append('')

        # Final
        if s.dict['dumpBR'] & 1:
            s.geneaLog.Add(genea=genea, S='info', txt=80*'-')
            s.geneaLog.Add(genea=genea, S='info', txt=out)
            s.geneaLog.Add(genea=genea, S='info', txt=80*'-')
            
        if s.dict['dumpBR'] & 2: WriteToFile(fn='%s%s%s%s_BRdump.txt' %(s.dict['fn0'],fn.replace('.txt','').replace('.slha',''), s.fntag, fnadd), outs=outs, VB=s.VB-1)
        
        return outs        





    # ##########################################################################
    # ##########################################################################
    # ##########################################################################
    # ##########################################################################

    # ########## SHARP
    def GetOnelegs(s, genea, fn, ID, particles, SUSYinitiators=[]):
        """
        Use s.SUSYinitiators1leg or (if empty) take all from xsecs : either way, allow part_minmass to skip sparticles
        Set smallest masses first (also SM)
        From particle (LegacyParticle) which has decayList equal to): do recursive 
        """
        genea = s.geneaLog.Genea(genea)
        #print 50*'a' ; particles['~C1'].ShowDecays()
        #for d in particles['~C1'].FSlist: print d.children
        #if s.VB>0: s.timer.Add('GetOnelegs:  start',-1)
        if SUSYinitiators == []: SUSYinitiators = list(particles.keys())
        SUSYinitiators1legSorted = sorted(SUSYinitiators, key=lambda x: abs(particles[x].mass))   # lightest is now first
        s.geneaLog.Add(genea=genea, S='info', txt='now starting')
        s.geneaLog.Add(genea=genea, S='debug', txt='SUSYinitiators1legSorted = %s' %(SUSYinitiators1legSorted))
        s.geneaLog.Add(genea=genea, S='debug', txt='particles.keys() = %s' %(list(particles.keys())))


        Mothers_Onelegs = {}
        # Ex Old:  1legsDict['~C1'] = [ [0.17,'~N1','q','q'], [0.23,'~N1','l','v'], ...]   # will probably need to make more advanced
        # Ex New:  1legsDict['~C1'] = [ [0.17,{'~N1':1, 'q':2}], [0.23,{'~N1':1, 'l':1, 'v':1}], ...]



        if 0:  # This is actually needed. NO, IT IS NO LONGER NEEDED. ALSO NOT DANGEROUS, WILL JUST SKIP IF EXISTS
            zaddthese = []
            zaddthese += ['x']
            #zaddthese += ['~N1','y','v','q','L',  'u','d','c','s','l','Th','b','g','ve','vm','vT', 'e','m'] 
            #zaddthese += ['T','W','Z','h','t']
            #zaddthese = ['W']   # the problem
            SUSYinitiators1legSorted = zaddthese + SUSYinitiators1legSorted
            #for zadd in reversed(zaddthese):
            #    if zadd not in SUSYinitiators1legSorted: 
            #        SUSYinitiators1legSorted = [zadd] + SUSYinitiators1legSorted
            #        print 'added: %s' %(zadd)
            print('\nNB: adding %s to SUSYinitiators1legSorted which  then becomes %s' %(zaddthese, SUSYinitiators1legSorted))
            print('  THIS ACTUALLY SAVES THE DAY (PARTLY) ; NOW THERE IS NO PROBLEM WITH C1 AND N2 ANYMORE')
            print('  HOWEVER THERE ARE PROBLEMS WITH HEAVIER SPARTICLES. IT SEEMS THIS IS DUE TO THEIR DECAY INTO LIGHTER SPARTICLES. (IF THESE LIGHTER SPARTICLES ARE PUT IN PARTICLES2KEEP, THEY SUM TO ONE.)')
            print(' SO SOMEHOW THE PROBLEM IS THE USE OF OTHER PARTICLES WITH COMPLEX DECAY')
            print('\n')

        
        SUSYinitiators1legSorted.insert(0,'x')  # need this to avoid eternal loop in ExpandDecay

        # ### A Loop over initiators
        for iPart in range(len(SUSYinitiators1legSorted)):
            pname_initiator = SUSYinitiators1legSorted[iPart]
            p = particles[pname_initiator]
            if not ( s.dict.get('1leg_minmass',0.) <= abs(p.mass) <= s.dict.get('1leg_maxmass',99999.) ):
                s.geneaLog.Add(genea=genea, S='info', txt='skipping %s of mass %.3f' %(pname_initiator, p.mass))
                continue

            s.geneaLog.Add(genea=genea, S='info', txt='Loop over initiators: iPart=%i  %s' %(iPart, pname_initiator))

            if pname_initiator in Mothers_Onelegs:
                s.geneaLog.Add(genea=genea, S='info', txt=' Note: skipping %s since already on Mothers_Onelegs' %(pname_initiator))
                continue
            
            thisMother_Onelegs = Mother(pname=pname_initiator, particles2keep=s.particles2keep, geneaLog=s.geneaLog, genea=genea)


            if s.VB>1:
                if 'e' in Mothers_Onelegs: s.geneaLog.Add(genea=genea, S='info', txt="(first: Mothers_Onelegs['e']: FSlist: %s" %(Mothers_Onelegs['e'].FSlist))
                s.geneaLog.Add(genea=genea, S='debug', txt='GetOnelegs  ipart: %i  %s   ndecays: %i' %(iPart, pname_initiator, len(p.FSlist)))
                s.geneaLog.Add(genea=genea, S='debug', txt='   Decay before start developing:')
                s.geneaLog.Add(genea=genea, S='debug', txt=p.ShowDecays(opt=['ret','marklast']))   # will contain NBEXTRA since e.g. e,m,T are not yet developed (to l,Th)
                s.geneaLog.Add(genea=genea, S='info', txt='   ---------------')
                

            # I am here sending in the Mothers_Onelegs_existing
            thisMother_Onelegs.DevelopFinalstates(FSlist=p.FSlist, particles=particles, Mothers_Onelegs_existing=Mothers_Onelegs, optD=s.dict, VB=s.VB, genea=genea)

            s.geneaLog.Add(genea=genea, S='debug', txt='Built: %4s : finalstates= %2i    SumBR=%9.7f  (width=%9.7f)' %(thisMother_Onelegs.name, thisMother_Onelegs.N(), thisMother_Onelegs.ShowSumBR(), thisMother_Onelegs.width))

            if 1-thisMother_Onelegs.ShowSumBR() > s.dict['BRmaxdev1leg']:
                s.warn.append('Warning::GetOnelegs  Too low SumBR for %-4s: %.6f' %(pname_initiator, thisMother_Onelegs.ShowSumBR()))
                s.geneaLog.Add(genea=genea, S='warn', txt=s.warn[-1])


            #thisMother_Onelegs.JoinIdenticalFSs()  # 2014-08-27  attempt ...
            Mothers_Onelegs[pname_initiator] = thisMother_Onelegs

            s.geneaLog.Add(genea=genea, S='debug', txt='just built %-4s :  len(Mothers_Onelegs) = %2i  : %s' %(thisMother_Onelegs.name, len(Mothers_Onelegs), list(Mothers_Onelegs.keys())))

            #s.geneaLog.Add(genea=genea, S='dump', txt='\nAA\nAA\nAA')

        # ===== end of loop A over initiators

        return Mothers_Onelegs
    
    



    # ########## SHARP
    def GetTwolegs(s, genea, Mothers_Onelegs, subprocs=[], xsecs={}, fn='fn', ID='ID'):
        """
        The encapsulating structure is mostly taken from Make2legs (legacy) ;
        In the processes of rewriting (restarted 2013-12-22, last 2013-12-16)
        """
        genea = s.geneaLog.Genea(genea)

        s.geneaLog.Add(genea=genea, S='info', txt='now starting')


        # Get subprocess list: either input (from -subprocs <> , or from xsecs, or from xsec-template [in that order])
        subprocs = list(subprocs)  # think I need this in case of empty subprocs and subsequent filling from xsecs.keys()
        # If empty, take all from xsecs 
        if subprocs == []:
            for prospcode in xsecs:
                if prospcode == 0 or prospcode >= 249: continue
                txttuple = tuple(sorted(list(prospino_code2txtpair.get(prospcode,[]))))
                if txttuple == (): s.geneaLog.Add(genea=genea, S='warn', txt='Warning::GetTwolegs  Unrecognised code: %i' %(prospcode)) ; continue
                subprocs.append(txttuple)


        if s.fn_xsecs_template:
            xsecTot = xsecs[0].Interpolate(s.xsecs_template_x, s.xsecs_template_y)
        else: 
            xsecTot = xsecs.get(0,{}).get('xs',-1.)  # use to allow cutting on relative xsec [might also add possibility to calc on the fly]

        #print 'xsecTot: ', xsecTot  # 

        usexsec = 1
        if not s.fn_xsecs_template and not xsecs: usexsec = 0 ; xsecTot = 1.
        

        # Then build the 2leg of each subprocess
        Twolegs = {}
        for isubproc in range(len(subprocs)):
            # Init
            subproctuple = subprocs[isubproc]

            if s.subprocs2leg != [] and subproctuple not in s.subprocs2leg: continue
            
            p1signed,p2signed = subproctuple
            prospcode = func_prospino_txtpair2code(subproctuple)
            p1generic = p1signed.rstrip('+-_')
            p2generic = p2signed.rstrip('+-_')
            if s.fn_xsecs_template:
                xsecSub0 = xsecs[prospcode].Interpolate(s.xsecs_template_x, s.xsecs_template_y)
            elif xsecs: 
                xsecSub0 = xsecs.get(prospcode,{}).get('xs',-1e-12)  # IF MODE XSEC_UNCERT
            else:
                xsecSub0 = 1.
            
            if xsecTot < 0: xsecRel0 = -1
            else: xsecRel0 = xsecSub0/xsecTot


        
            s.geneaLog.Add(genea=genea, S='debug', txt='%-10s  xsecTot: %9.6f   xsecSub0: %9.6f   xsecRel0: %9.6f' %(subproctuple, xsecTot, xsecSub0, xsecRel0))

            # Skip-test #0
            if xsecs:
                if xsecTot > 0 and xsecRel0 < s.dict.get('xsecRelMin',0.) :
                    subproc_skipped[subproctuple] = {'xsecRel': xsecRel0}
                    continue
                if xsecSub0 < s.dict.get('xsecMin',0.) :
                    subproc_skipped[subproctuple] = {'xsec': xsecSub0}
                    continue

            #if not usexsec:
            #    xsecTot = xsecSub0 = xsecRel0 = 1.  # 2014-08-28: hack (should clean up all of these xsec things)
        

            if p1generic not in Mothers_Onelegs:
                s.geneaLog.Add(genea=genea, S='warn', txt='Warning::GetTwolegs  Skipping subprocess %s since %s not in Mothers_Onelegs.keys(): %s' %(subproctuple, p1generic, list(Mothers_Onelegs.keys())))
                continue
            if p2generic not in Mothers_Onelegs:
                s.geneaLog.Add(genea=genea, S='warn', txt='Warning::GetTwolegs  Skipping subprocess %s since %s not in Mothers_Onelegs.keys(): %s' %(subproctuple, p2generic, list(Mothers_Onelegs.keys())))
                #print 'Warning::GetTwolegs  %s not in Mothers_Onelegs.keys(): %s' %(p2generic, Mothers_Onelegs.keys())
                continue

            # For entire subproc: xsecSub0, xsecRel0 (directly from xsec dict)


            # Create the twoleg object .. might prefer a real object
            pair = (p1signed, p2signed)
            pairname = str(pair).replace('(','').replace(')','').replace(',','_').replace(' ','').replace("'","")
            twoleg = Mother(pairname, particles2keep=s.particles2keep, geneaLog=s.geneaLog, genea=genea) # BlackBoxTwoleg object?
            Twolegs[pair] = twoleg
            # Additional attributes (beyond Mother)
            twoleg.pair = tuple(pair)
            twoleg.xsecSub = xsecSub0
            twoleg.xsecRel = xsecRel0
            twoleg.width = twoleg.xsecSub  # 2014-01-14
            #twoleg.totalwidth = twoleg.xsecSub
            s.geneaLog.Add(genea=genea, S='dump', txt='xxxsec: %.9f' %(twoleg.width))


            ### NEED TO CLEAN UP THE USE OF .xsecSub and .xsecRel .. these are not attributes of Mother ..
            
            # Combine all 1legs for these two initiators            
            


            # Loop over decays for p1
            ###for i1 in range(len(Onelegs[p1generic])):  # Onelegs[p1name] is a list of Mother
            for i1 in range(Mothers_Onelegs[p1generic].N()):  # Onelegs[p1name] is a list of Mother  #BBFL#
                #casc1 = Mothers_Onelegs[p1generic][i1]  # is a FS
                casc1 = Mothers_Onelegs[p1generic].FSlist[i1]  # is a FS
                xsecSub = xsecSub0 * casc1.br
                xsecRel = xsecRel0 * casc1.br
                s.geneaLog.Add(genea=genea, S='dump', txt='   i1: %i    %4s  thisBR: %10.8f    sub: %10.8f   rel: %10.8f' %(i1, p1generic, casc1.br, xsecSub, xsecRel))

                # Skip-test #1
                if xsecSub < s.dict.get('xsecMin',0.): continue
                if xsecRel < s.dict.get('xsecRelMin',0.): continue
                if casc1.br < s.dict.get('BRmin2leg',0.): continue

                

                # Loop over decays for p2
                ###for i2 in range(len(Onelegs[p2generic])): 
                for i2 in range(Mothers_Onelegs[p2generic].N()):   #BBFL#
                    # if VB>3: print '   i2: %i' %(i2)
                    ###casc2 = Onelegs[p2generic][i2]
                    casc2 = Mothers_Onelegs[p2generic].FSlist[i2]  #BBFL#
                    xsecSub = xsecSub0 * casc1.br * casc2.br
                    xsecRel = xsecRel0 * casc1.br * casc2.br
                    
                    s.geneaLog.Add(genea=genea, S='dump', txt='     i2: %i   %4s  thisBR: %10.8f    sub: %10.8f   rel: %10.8f' %(i2, p2generic, casc2.br, xsecSub, xsecRel))

                    # Skip-test #2
                    if usexsec: 
                        if xsecSub < s.dict.get('xsecMin',0.): continue
                        if xsecRel < s.dict.get('xsecRelMin',0.): continue
                        if casc1.br*casc2.br < s.dict.get('BRmin2leg',0.): continue
                    
                    
                    # Combine, using Mother as container
                    #twoleg.brs = [casc1.br, casc2.br]
                    #twoleg.pnames = [p1signed,p2signed]
                    # [MORE TO IMPLEMENT]

                    # Create this object and insert into twoleg
                    thisdecay = FS(casc1.children, casc1.br)
                    thisdecay.MultWithFS(casc2)
                    twoleg.AddFS(thisdecay)



            # Need here simplify / join the decays 
            twoleg.JoinIdenticalFSs()
            twoleg.SortDecaysWithBR()
            s.geneaLog.Add(genea=genea, S='dump', txt='twoleg: %.7f' %(twoleg.width))  # width complete

            if 1-twoleg.ShowSumBR() > s.dict['BRmaxdev2leg']:
                s.warn.append('Warning::GetTwolegs  Too low SumBR for %-10s: %.6f' %(pairname, twoleg.ShowSumBR()))
                s.geneaLog.Add(genea=genea, S='warn', txt=s.warn[-1])
            
            

        s.geneaLog.Add(genea=genea, S='info', txt="GetTwolegs  number of Twolegs: %i:  %s" %(len(Twolegs), list(Twolegs.keys())))
        if s.VB>1:
            for pair in Twolegs:
                twoleg = Twolegs[pair]
                s.geneaLog.Add(genea=genea, S='debug', txt="GetTwolegs   %-10s  %8.6f   %10.6f pb | %3i" %(twoleg.name, twoleg.xsecRel, twoleg.xsecSub, twoleg.N()))

        if s.VB>2:
            for pair in Twolegs:
                twoleg = Twolegs[pair]
                s.geneaLog.Add(genea=genea, S='dump', txt="GetTwolegs   %-10s  %8.6f   %10.6f pb | %3i" %(twoleg.name, twoleg.xsecRel, twoleg.xsecSub, twoleg.N()))
                twoleg.ShowDecays(opt=['show'], which=['~N1','b','T', 'l','q','y','v'])
                

        res = {'Twolegs':Twolegs}

        return res

                

    # ########## SHARP
    def FindAllFinalstates(s, genea, fn, xsecs={}, ID=''):  # corresponds to DevelopCascade
        """
        This is for a given scenario (fn, xsecs)
        """
        
        genea = s.geneaLog.Genea(genea)
        #if s.VB>=0: s.timer.Add("FindAllFinalstates: start",-1)
        s.geneaLog.Add(genea=genea, S='info', txt='now starting...')
        fn_clean = fn.replace('.slha','').replace('.txt','')  # maybe replaced by ID
        if ID == '': ID = fn_clean
        
        res = {}

        # #################################################  
        # #################################################   PART 0: particle dict
        # #################################################
        s.geneaLog.Add(genea=genea, S='info',txt=40*'#'+' PART 0: CONSTRUCT DICT particles') 

        
        # 1 Read slha        
        slha = pyslha.readSLHAFile(fn)
        #if s.VB>=0: s.timer.Add("DevelopCascade: after readSLHAFile",-1)


        # 1.5 Get flatdict for SLHA[0] object (as in slh2tlt)
        if s.dict['flatdict_slha0']: 
            s.flatdict_slha0 = slha2flatdict(SLHA=slha)
            s.flatdict.update(s.flatdict_slha0['flat'])
        #print 'len(s.flatdict) = %i' %(len(s.flatdict))
        #print 'flatdict.keys() = %s' %(s.flatdict.keys())
        #print 'flatdict = %s' %(s.flatdict)
        

        # 2 Get some additional info from slha  (mass of W,Z,t,T,b)
        slhaInfo = s.GetInfo(genea=genea, slha=slha)  

        # 3 SlhaBR->isaPlayBR : Make a dict of particles with decaylist
        particles = s.GetMotherDict(genea=genea, slha1=slha[1], optD={'masses':slhaInfo['masses']} )
        # typically does not contain W,Z, .. but has h and t

        # Note: ~N1 exists in particles with non-zero mass 
        
        #if s.VB>=0: s.timer.Add("DevelopCascade: after GetParticleDict_ala_IsaPlay",-1)
        if s.dict['dumpBR']: 
            s.geneaLog.Add(genea=genea, S='debug',txt="Before AddSMtoMothers: ")
            s.MotherDump(genea=genea, fn=fn, particles=particles)  # DEBUG


        # 3b
        AddSMtoMothers(geneaLog=s.geneaLog, genea=genea, particles=particles, PDG=PDG, slhaInfo=slhaInfo, VB=s.VB)
        if s.dict['dumpBR']: 
            s.geneaLog.Add(genea=genea, S='debug',txt="After AddSMtoMothers: ")
            s.MotherDump(genea=genea, fn=fn, particles=particles, fnadd='_addSM')  # DEBUG 
        # W,Z,.. have e.g. 'e' in their final states rather than e.g. 'l' that is in particles2keep
        # h also has this


        # 4 Verbose
        #if s.dict['dumpBR']:
        #    s.ParticleDump(fn=fn, particles=particles, which=s.SUSYinitiators1leg, optD={'masses':slhaInfo['masses']} )
        #    s.timer.Add("DevelopCascade: after dumpBR",-1)

        #if s.VB>3:
        #    for pname in particles: print 'DEBUG::FindAllFinalstates (again)  %-4s  %.3f' %(pname, particles[pname].mass)




        # #################################################  
        # #################################################   PART 1: 1 legs
        # #################################################
        if s.dict['runmode'] >= 1: 

            #particles['W'].DevelopFinalstates()  # fails
            s.geneaLog.Add(genea=genea, S='debug',txt=['BEFORE 1 LEG',particles['W'].ShowDecays(opt=['ret','marklast'])])


            s.geneaLog.Add(genea=genea, S='info',txt=40*'#'+' PART 1: CONSTRUCT 1 LEGS', space=[2,0]) 
            # 5 1legs
            Mothers_Onelegs = s.GetOnelegs(genea=genea, fn=fn, ID=ID, particles=particles, SUSYinitiators=s.SUSYinitiators1leg)
            # s.Onelegs_Show(Mothers_Onelegs=Mothers_Onelegs, fnslha=fn)  # jada, hm

            s.geneaLog.Add(genea=genea, S='debug',txt=['AFTER 1 LEG',particles['W'].ShowDecays(opt=['ret','marklast'])])



            # verbosity (& to file)
            outs1legs = []
            for onelegID in Mothers_Onelegs:
                outs1legs += Mothers_Onelegs[onelegID].ShowDecays(opt=['ret','marklast'], what=['RELacc','REL'])
                outs1legs += ['']
            if s.dict['1legstofile']:
                fn_onelegs = '%s%s%s_onelegs%s.txt' %(s.dict['fn0'], fn_clean, s.fntag, s.fsinfilename)
                WriteToFile(fn=fn_onelegs, outs=outs1legs, VB=s.VB-1)
            if s.VB>0:  # nicify (genea)
                s.geneaLog.Add(genea=genea, S='info', txt=20*'*' + '1legs')
                for out in outs1legs: s.geneaLog.Add(genea=genea, S='info', txt=out)


            # 6 Put 1legs in flat dict [opt]
            if 'onelegs' in s.dict['flatdict']:
                s.flatdict.update(s.PutInFlatDict(genea=genea, legs=Mothers_Onelegs, mode='onelegs', gteqlt=s.dict['flatdict_gteqlt'], npartMax=s.dict['flatdict_npartMax'], npartMaxCommon=s.dict['flatdict_npartMaxCommon']))




        # #################################################  
        # #################################################   PART 2: 2 legs
        # #################################################  
        # 7 2legs
        if s.dict['runmode'] >= 2: 
            s.geneaLog.Add(genea=genea, S='info',txt=40*'#'+' PART 2: CONSTRUCT 2 LEGS') 
            res_twolegs = s.GetTwolegs(genea=genea, Mothers_Onelegs=Mothers_Onelegs, xsecs=xsecs, subprocs=s.subprocs2leg, fn=fn, ID=ID)  # 2014-08-28: added s.subprocs2leg
            twolegs = res_twolegs['Twolegs']
            IDs_subprocs = list(twolegs.keys())



            # 8 sum the 2legs (only if have xsecs) [to be implemented]
            sumID = 'Sum'  # make dynamic? 
            if xsecs: 
                sumMothers = SumMothers(list(twolegs.values()), name=sumID, particles2keep=s.particles2keep)
                sumMothers.RemoveFinalStatesWithLowBR(s.dict['BRmin2leg'])
                twolegs[sumID] = sumMothers


            # 9 put 2legs in flat dict
            if 'twolegs' in s.dict['flatdict']:  # is default
                s.flatdict.update(s.PutInFlatDict(genea=genea, legs=twolegs, mode='twolegs', gteqlt=s.dict['flatdict_gteqlt'], npartMax=s.dict['flatdict_npartMax'], npartMaxCommon=s.dict['flatdict_npartMaxCommon']))


            
            # 10 to txt file
            outs2legs = []
            for twolegID in IDs_subprocs + [sumID]:
                if not twolegID in twolegs: continue  # hack to not crash on sumID if no xsec
                outs2legs += twolegs[twolegID].ShowDecays(opt=['ret','marklast'], formats={'parent':'%-9s'})
                outs2legs += ['']
                s.geneaLog.Add(genea=genea, S='debug', txt='twolegs[ID].SumBR() : %.7f' %(twolegs[twolegID].ShowSumBR()))
                s.geneaLog.Add(genea=genea, S='debug', txt='twolegs[ID].width   : %.7f' %(twolegs[twolegID].width))
                # outs2legs += ['']
                # outs2legs += sumMothers.ShowDecays(opt=['ret'])
                fn_twolegs = '%s%s%s_twolegs%s.txt' %(s.dict['fn0'], fn_clean, s.fntag, s.fsinfilename)
                WriteToFile(fn=fn_twolegs, outs=outs2legs, VB=s.VB-1)

                if s.VB>0: # nicify (genea)
                    s.geneaLog.Add(genea=genea, S='info', txt=80*'*'+'2legs')
                    for out in outs2legs: s.geneaLog.Add(genea=genea, S='info', txt=out)

        
                # 11 to pickle (might put in a superdict with additional info)
                fn_2legs_pickle = '%s%s%s_twolegs%s.pickle' %(s.dict['fn0'], fn_clean, s.fntag, s.fsinfilename)
                WritePickle(fn=fn_2legs_pickle, thepickle=twolegs, VB=s.VB-1)
        

        # #################################################  
        # ################################################# flatdict: 1leg, 2leg
        # #################################################  
        if s.dict['runmode'] >= 1: 

            # 12 store flat dict [1leg, 2leg]
            fn_flatdict = '%s%s%s_flatdict%s.pickle' %(s.dict['fn0'], fn_clean, s.fntag, s.fsinfilename)
            WritePickle(fn=fn_flatdict, thepickle=s.flatdict, VB=s.VB-1)

            # 12b store a brute list of the flatdict (for easy looking up and grepping) [1leg, 2leg]
            fn_flattxt = fn_flatdict.replace('.pickle','.txt')
            outs_flattxt = []
            for key in sorted(s.flatdict.keys()): 

                val = s.flatdict[key]
                if type(val) == 'float': 
                    outs_flattxt.append('%-30s  %14.8f' %(key, val))
                elif type(val) == 'int': 
                    outs_flattxt.append('%-30s  %i' %(key, val))
                else: 
                    outs_flattxt.append('%-30s  %s' %(key, val))

            WriteToFile(fn=fn_flattxt, outs=outs_flattxt, VB=s.VB-1)


        # tlt? 
        for tltdefname in s.tlts: 
            tltdef,tltname = tltdefname.split(s.dict['tlt_delim_name'])
            outstlt = dict2tlt(tltdef=tltdef, thedict=s.flatdict, D=s.dict, formats=s.dict['tlt_formats'], var2block=s.flatdict_slha0.get('var2block',{}), format0=s.dict['tlt_format0'], delim_var=s.dict['tlt_delim_var'], delim_rename=s.dict['tlt_delim_rename'])
            fn_tlt =  '%s%s%s_%s.tlt' %(s.dict['fn0'], fn_clean, s.fntag, tltname)
            WriteToFile(fn=fn_tlt, outs=outstlt, VB=s.VB-1)
            
            



    # ########## (It might be that the effective parts of this should be moved to SlhaTools)
    def PutInFlatDict(s, genea, legs, mode, whichlegs=[], npartMax={'L':3,'l':3,'Th':3}, npartMaxCommon=2, npartMin={}, npartMinCommon=0, gteqlt=['>=','=','<='], opt=[], drop=[], delimDecay='_', delim0=':'):

        genea = s.geneaLog.Genea(genea)
        #if s.VB>0: print 'INFO::PutInFlatDict  Starting...'
        s.geneaLog.Add(genea=genea, S='info', txt='now starting')
        gteqlts = ['>=','=','<=']  # allowed values
        
        res = {}


        
        # INIT: Need to do this in a recursive-like way

        # 1) Build an iterator system with n particles, each going between a min and a max, e.g. [0,1,2]
        iterRangeVars = {}
        for zname in s.particles2keep:
            iterRangeVars[zname] = list(range(npartMin.get(zname, npartMinCommon), npartMax.get(zname, npartMaxCommon) + 1))

        resIterboxVars = IterateNdim(vars=s.particles2keep, varDict=iterRangeVars)
        combsVars = resIterboxVars['combinations']
        s.geneaLog.Add(genea=genea, S='debug', txt='DEBUG::PutInFlatDict  n(combsVars) = %i' %(len(combsVars)))
        if s.VB>2:
            for iz in range(len(combsVars)): s.geneaLog.Add(genea=genea, S='dump', txt='   combsVars  %3i   %s' %(iz, combsVars[iz]))

        # 2) Build an iterator system for >=<
        #gteqltArr = []
        for z in gteqlt:
            if z not in gteqlts:
                sys.exit("Error::PutInFlatDict  unrecognised gteqlt: '%s' (len=%i)  not among allowed ones %s" %(gteqlt, len(gteqlt), gteqlts))
            
        iterRangeGteqlt = {}
        for zname in s.particles2keep:
             # It is here possible to optionally skip certain combinations (e.g. keep only >= for L)
            iterRangeGteqlt[zname] = list(gteqlt) 
        resIterboxGteqlt = IterateNdim(vars=s.particles2keep, varDict=iterRangeGteqlt)
        combsGteqlt = resIterboxGteqlt['combinations']
        


        # A1 Loop over all 1legs (e.g. N2->) / 2legs (e.g. C1+C1- ->)

        if whichlegs: legnames = list(whichlegs)
        else: legnames = sorted(legs.keys())
        for ilegname in range(len(legnames)): 
            legname = legnames[ilegname]
            leg = legs[legname]
            
            
            if leg.name in ['d','u','c','s','b','q','~N1','y','v','ve','vm','vT','Th','L','l','T','e','m','g']: continue
            
            # A2 Loop over all possible var combinations
            for icombVars in range(len(combsVars)):
                combVars = combsVars[icombVars]
                s.geneaLog.Add(genea=genea, S='dump', txt='A2 Loop: combVars: %2i  %s' %(icombVars, combVars))
                
                # initialise variables
                # (for now only one variable ... but this may explode, e.g. to give combinations like
                # BR:N22Th,1+h, BR:N2>2Th,1h, BR:N2>2+Th,1h, BR:N2>2+Th,1+h


                # <could also have the "loop over all decays in this 1leg" here>
                

                # A3a Loop over XS/BR type
                #for BRorXS in ['BR','XS']: # should not have an outer/intermediate loop like this, should just do the loop at the very end
                #    if BRorXS == 'XS' and mode == 'onelegs': continue  
                if 1: 

                    # A3b Loop over Gteqlt combination
                    for icombGteqlt in range(len(combsGteqlt)):
                        combGteqlt = combsGteqlt[icombGteqlt]  # ex: combGteqlt=['>','=','>'] if there are three vars, e.g. Th,l,b
                        s.geneaLog.Add(genea=genea, S='dump', txt='combGteqlt = %s' %(combGteqlt))

                        #var = BRorXS + ':' + leg.name   # ex: BR_~C1 or BR_~
                        var = ':' + leg.name   # ex: BR_~C1 or BR_~
                        for ip in range(len(s.particles2keep)):
                            p2keep = s.particles2keep[ip]
                            if ip == 0: var += delim0
                            else: var += delimDecay
                            var += '%s%i%s' %(combGteqlt[p2keep], combVars[p2keep], p2keep)  # adds e.g. '>2Th' for a given 
                        if '~' in drop: var = var.replace('~','')   # BR_~C1  ->  BR_C1
                        if '=' in drop: var = var.replace('=','')  # allows to drop '=' from var name
                        
                        # Now have the final var name to be put in flatdict
                        # Will then loop over fss and compare the given var-construction
                        
                        s.geneaLog.Add(genea=genea, S='dump', txt='  %2i  %-10s   %s   %s   %s' %(ilegname, leg.name, combVars, combGteqlt, var))
                            

                        # init (could also do below)
                        res['BR'+var] = 0.
                        if mode == 'twolegs': res['XS'+var] = 0.

                        # --- 
                            
                        
                        # A3 Loop over all decays in this 1leg
                        for ifs in range(leg.N()):
                            fs = leg.FSlist[ifs]


                            # Now compare given decay chain with comb and fill
                            # dicts with given setting for each p2keep: combGteqlt ({'Th':'>','l':'='}), combVars ({'Th':2, 'l':1}, fs.children
                            ch = fs.children
                            isOK = 1
                            for p2keep in s.particles2keep:
                                if combGteqlt[p2keep] == '>=' and not (ch.get(p2keep,0) >= combVars[p2keep]): isOK = 0; break
                                if combGteqlt[p2keep] == '='  and not (ch.get(p2keep,0) == combVars[p2keep]): isOK = 0; break
                                if combGteqlt[p2keep] == '<=' and not (ch.get(p2keep,0) <= combVars[p2keep]): isOK = 0; break
                                # Note, for children (ch): need to use get, because key is not in dict if value is zero [could change elsewhere]

                            
                            # Now add to variable
                            if isOK:
                                res['BR'+var] += fs.br
                                if mode == 'twolegs': res['XS'+var] += fs.br * leg.width


                

        return res

        
    # ##########
    def Main(s, genea):
        genea = s.geneaLog.Genea(genea)

        #if s.VB>0: print "INFO::%s  Main" %(s.myname)
        #if s.VB>=0: s.timer.Add("Main: start")
        s.geneaLog.Add(genea=genea, S='info', txt="now starting")
        
        # Get xsec, typically one common file (can have option to get per file as well)
        s.xsecs = {}

        # Standard 1: point dictionary
        if s.fn_xsecs:
            tmpDict = {}
            exec(compile(open(s.fn_xsecs, "rb").read(), s.fn_xsecs, 'exec'), tmpDict)
            s.xsecs = tmpDict[s.xsecs_dictname]
            del tmpDict
            if s.VB>0: 's.xsecs.keys(): %s' %(list(s.xsecs.keys()))

        # Standard 2: interpolation templates
        elif s.fn_xsecs_template:
            s.xsecs = LoadPickle(fn=s.fn_xsecs_template) # Needs to be a plain/flat template with prospino codes as keys


        #if s.VB>=0: s.timer.Add("Main: after xsec",-1)
        #s.geneaLog.Add(genea=genea, S='info', txt="after xsec")   # add timer info ...
        
        # Prepare initiators
        s.SUSYinitiators1leg, s.subprocs2leg = s.GetInitiators(genea=genea)
        s.geneaLog.Add(genea=genea, S='info', txt='s.SUSYiniators1leg: %s' %(s.SUSYinitiators1leg))
        s.geneaLog.Add(genea=genea, S='info', txt='s.subprocs2leg: %s' %(s.subprocs2leg))


        # Loop over scenarios: develop cascades (1leg and 2leg) inside
        if s.VB>5: print('s.fns: %s' %(s.fns))
        for iF in range(len(s.fns)):
            fn = s.fns[iF].strip()

            # Cross-sections [should be inside an if-statement?]
            xseckey = tuple(ExtractVarsViaFormat(fn, s.slhaformat))
            ID = ''
            for z in xseckey: ID = (ID + '_' + z).lstrip('_')  # somewhat limited ... 
            #print 'DEBUG xseckey:', xseckey
            s.geneaLog.Add(genea=genea, S='info', txt='ID=%s  fn=%s' %(ID, fn))

            # VERY UNDER THE HOOD: 
            if s.fn_xsecs_template: 
                s.xsecs_template_x = float(xseckey[1])   # FRAGILE & UNDER THE HOOD
                s.xsecs_template_y = float(xseckey[2])   # FRAGILE

            
            if s.fn_xsecs_template:
                xsecs = s.xsecs   # send along the entire thing (Standard 2)            
            else:
                xsecs = s.xsecs.get(xseckey, {})
                
            if not xsecs: s.warn.append("No xsec dict for %s   - Twolegs will be empty." %(fn))


            
            #if s.dict['runmode'] > 0:  # only option [CLEAN UP]
            s.DEISM[fn] = s.FindAllFinalstates(genea=genea, fn=fn, xsecs=xsecs, ID=ID)
            #else:
            #    print 'Warning::slha2finalstates  NOTHING IS RUN. Need to have runmode > 0   (for now)'
            
            #if s.VB>=0: s.timer.Add("Main: After FindAllFinalstates for %i %s" %(iF, fn.split('/').pop()), -1)


        
            s.geneaLog.Add(genea=genea, S='info', txt='\n%s'%(100*'='))
            s.geneaLog.Add(genea=genea, S='info', txt='NOTE: if not all subprocs in scenario are used, the 2leg_9_final BRaccu will not go to its integer value, but to something lower. (BRsub similary affected since its sum equals BRaccu at the various stages.')
            s.geneaLog.Add(genea=genea, S='info', txt='      Same effect in 2leg_0_untouched, will get same number as in 2leg_9_final.')



        WriteToFile(fn=s.fn_warn, outs=s.warn, VB=s.VB-1)


    # ##########
    def ReadArg(s, genea): 

        genea = s.geneaLog.Genea(genea)
        # ################################### ARGUMENT READING
        Arg = bkgjelstenArgReader.ArgReader(s.argv, VB=0)

        #if Arg.hasget('-alist'):  print 'a string list: ',Arg.list()
        #if Arg.hasget('-alisti'): print 'an integer list: ',Arg.listI()
        #if Arg.hasget('-alistf'): print 'a float list: ',Arg.listF()
        #if Arg.hasget('-x'):  print 'a string: ',Arg.val()
        #if Arg.hasget('-xI'): print 'an integer: ',Arg.valI()
        #if Arg.hasget('-xF'): print 'a float: ',Arg.valF()

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
            os.system("cat %s | grep 'if Arg\.has'" %(s.argv[0]))
            print("[ The above shows all implemented options. For full details inspect the code, %s  (Try also '-h') ]" %(s.argv[0]))
            sys.exit()
            
        if Arg.has(['-h','--help','--h','-help']):
            s.showHelp()
            sys.exit()

        if Arg.hasget('-vb'):
            s.VB = Arg.valI()
            if s.VB>0: print('Verbosity level: %i' %(s.VB))



        if Arg.hasget(['-geneaVB']):   # main steering for structured output ('debug','info','warn','error','fatal')
            z = Arg.val()
            s.geneaLog.strengthVB = s.geneaLog.StrengthI(z)
            print(s.geneaLog.strengthVB)

        if Arg.has('--silent'): 
            s.geneaLog.strengthVB = s.geneaLog.StrengthI('fatal')
            s.VB = 0

        if Arg.has(['--showparticles']):
            zparts = list(namGeneric.keys())
            zparts.sort()
            for zpart in zparts: 
                print('  %-5s  %s' %(zpart, namGeneric[zpart]))
            sys.exit()



        if Arg.hasget(['-genea','-geneaadd']):
            s.geneaadd = Arg.val()
            #print 'DDD: ', s.geneaAdd
        if Arg.hasget(['-geneainit']): 
            s.geneainit = Arg.val()

        if Arg.hasget(['-geneadepth']): 
            s.geneaLog.maxdepth = Arg.valI()

        if Arg.hasget(['-geneawidth1']): 
            s.geneaLog.maxwidth1 = Arg.valI()


        #if Arg.hasget('-f'): s.fns = Arg.list()  # 2014-08-27: can this be correct??
        if Arg.hasget('-f'): s.fns = os.popen('cat %s' %(Arg.val())).readlines()

        if Arg.hasget('-generic'): s.generic = Arg.list()

        if Arg.hasget(['-particles2keep']):
          s.particles2keep = Arg.list()
          # TODO: add test that these all exist

        #if Arg.hasget(['-SUSYinitiators1leg','-initiators1leg','-in1leg']):
        #    s.SUSYinitiators1leg = Arg.list()


        if Arg.hasget(['-SUSYinitiators','-initiators','-in']):  # this one sets also  s.SUSYinitiators1leg if empty
            s.SUSYinitiatorsBrute = Arg.list()

        if Arg.hasget(['-subprocs']):  # sets s.subprocs2legBrute, also sets s.SUSYinitiators1leg if empty
            zz = Arg.list()
            if zz != ['all']: # if user says -subprocs all, then set none (which means all)
                for z in zz:
                    w = z.split(':')
                    w.sort()  # alphatically ordered, that's the default standard (then in principle don't need pair)
                    s.subprocs2legBrute.append(w)
                
            
        if Arg.hasget(['-fn_xsecs','-xsecs']):
            s.fn_xsecs = Arg.val()

        if Arg.hasget(['-fn_xsecs_template','-xsecs_template']):
            s.fn_xsecs_template = Arg.val()

        if Arg.hasget('-slhaformat'):
            s.slhaformat = Arg.val()

        if Arg.hasget('-fntag'):
            s.fntag = Arg.val()

        if Arg.has('--fsinfile'):
            s.fsinfilename = 1
        if Arg.hasget('-fsinfile'):
            s.fsinfilename = Arg.valI()

        if Arg.hasget('-runmode'):
            s.dict['runmode'] = Arg.valI()

        if Arg.hasget('-flatdict_gteqlt'):  # Ex:  -gteqlt '>=,='
            s.dict['flatdict_gteqlt'] = Arg.list()

        if Arg.hasget('-flatdict_npartMaxCommon'):  # Ex:  -flatdict_npartMaxCommon 2
            s.dict['flatdict_npartMaxCommon'] = Arg.valI()

        if Arg.hasget('-flatdict_npartMax'):  # Ex:  -flatdict_npartMax L=3,b=4,q=1
            zs = Arg.list()
            for z in zs: 
                key,val = z.split('=')
                s.dict['flatdict_npartMax'][key] = int(val)

        if Arg.hasget(['-tlt','tlts']): 
            s.tlts = Arg.list(',,')


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
    t = slha2finalstates(cmd=['ReadArg','PostInit','Main'])
############################## 

