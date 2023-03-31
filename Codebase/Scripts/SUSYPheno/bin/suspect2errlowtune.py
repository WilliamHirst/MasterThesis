#!/usr/bin/env python
'''
Program : suspect2errlowtune.py
Version : 1.0
Author  : b.k.gjelsten@fys.uio.no
Description : 
This script just goes through the suspect2.out file and
produces simple output of
  warnings/errors
  finetuning variables
  lowenergy values
  

TODO:
x add how many stddev the given gmuminus2 and bsgam is away from exp value [or other relevant goal]
 
'''

import sys,os,math
from bkgjelstenArgReader import ArgReader
import pyslha_edited as pyslha
from lineup import lineup
from kilelib import WriteToFile, WritePickle
from libExpValues import exp

# ##################################################### GLOBAL METHODS
# ##################################################### GLOBAL METHODS
# ##################################################### GLOBAL METHODS

# ##########


# ##########
suspect_errorcode = {
    # Note: start at 1, like in suspec2 output
    1:   "tachyon 3rd gen. sfermion from RGE"    # ERR
    ,2:  "tachyon 1,2 gen. sfermion from RGE"
    ,3:  "tachyon A    (maybe temporary: see final mass) "  # WARN: "3   Warning:  MA^2(Q) <0 at a scale MZ<Q<EWSB !"
    ,4:  "tachyon 3rd gen. sfermion from mixing"
    ,5:  "mu unstable after many iter"   # WARN
    ,6:  "non-convergent mu from EWSB "
    ,7:  "EWSB maybe inconsistent        (!but RG-improved only check)"
    ,8:  "V_eff maybe UFB or CCB         (!but RG-improved only check)"
    ,9:  "Higgs boson masses are NaN "
    ,10: "RGE problems (non-pert and/or Landau poles)"  # ERR: "4   STOP: non-pert. R.C., or Landau pole in RGE"
    }


# ##########
class suspect2errlowtune: 

    def __init__(s, cmd=[], optD={}):

        # ====================== PRE INIT
        if 'argv' in optD: s.argv = optD['argv']
        else: s.argv = sys.argv
        s.cmd = cmd
        s.myname = sys.argv[0].split('/').pop()
        s.VB = 1
        s.HOME = os.getenv('HOME')
        s.cwd  = os.getcwd()  # current work directory  
        #s.dir0 = '%s/XXX' %(s.HOME)
        s.dir0 = '.'
        s.dir1 = '%s' %(s.dir0)

        s.dict = {}

        # s.dict['test'] = 'testval'
        # s.dict['test2'] = 'testval2'
        # s.dict['testI'] = 3
        # s.dict['testF'] = 4.34

        s.warn = []
        s.fn_warn = 'warnings.txt'
        s.fn_report = 'report'
        s.report = []

        s.fn = ''
        s.ID = ''
        s.pickle = 0

        s.prevar = ''
        s.preval = ''

        s.dryrun = 0
        s.nspaces = 3

        # ====================== READ ARG
        if 'ReadArg' in s.cmd: s.ReadArg()


        # ====================== POST INIT
        if 'PostInit' in s.cmd: s.PostInit()


        # ====================== EXECUTE 
        if 'Main' in s.cmd: s.Main()



    # ##########
    def PostInit(s): 
        s.fn_warn = '%s/%s' %(s.dir1, 'warnings.txt')
        s.fn_report = '%s/%s' %(s.dir1, 'report')

        
    

    # ##################################################### CLASS METHODS
    # ##################################################### CLASS METHODS
    # ##################################################### CLASS METHODS
    # ##########
    def showHelp(s):
        print(' Usage: %s [options]' %(s.myname))
        #print '        %s  -dict test,txt1:I,testI,3:test2,txt2a,txt2b:F,testF,4.14   # for using autodict (NB: vars need to be defined in __init__)' %(s.myname)


    # ##########
    def DumpWarnings(s):
        f = open(s.fn_warn,'w')
        for out in s.warn: f.write('%s\n' %(out))
        f.close()

        
    # ##########
    def ReadSuspectSPINFO(s, fn):
        # Need to also take a look at the slha file to see which are warnings and which are errors
        #slhaObj = pyslha.readSLHAFile(fn_slha)
        f = open(fn) ; lines = f.readlines() ; f.close()

        res = {3:[], 4:[], 'Err':0}  # 0:all ok;  1:warn;  2:err

        BLOCK = ''
        for iL in range(len(lines)): 
            L = lines[iL].strip()
            if L.startswith('#'): continue
            w = L.split()
            
            if L.startswith('BLOCK'):
                if BLOCK == 'SPINFO': return res  # Hack: if leaving SPINFO, just return whatever results we have
                BLOCK = w[1]
                continue

            # if in SPINFO block, collect warnings and errors
            if BLOCK == 'SPINFO':
                code = int(w[0])
                if code in [3,4]:   # 3:warning, 4:error
                    themessage = L[L.find(w[1]):]
                    res[code].append(themessage)
                    if code == 3: res['Err'] = max(res['Err'], 1)  # warning
                    if code == 4: res['Err'] = max(res['Err'], 2)  # error
            

        # shouldn't arrive here
        print('Warning::ReadSuspectSPINFO  should not arrive to the bottom')
        return res


    # ##########
    def ReadSuspectOut(s, fn):
        # for key in s.dict.keys(): print 'dict:  %-10s  %s' %(key, s.dict[key])

        # Read file
        if fn:
            f = open(fn) ; lines = f.readlines() ; f.close()
        else: 
            if s.VB>0: print('suspect2errlowtune::PostInit  Now reading input from stdin since no filename was given (-f)')
            lines = sys.stdin.readlines()



        res = {'caution':0}

        # Loop through and interpret
        
        for iL in range(len(lines)):
            L = lines[iL].rstrip()

            if L.startswith('CAUTION'):
                res['caution'] = 1



            if L.startswith('Low-energy/LEP precision'): 
                iL += 2
                L = lines[iL].rstrip()
                L = L.replace('E-','Eminus').replace('-',' -').replace('Eminus','E-')
                w = L.split()
                if len(w) != 3: s.warn.append('warning::suspect2errlowtune  Irregular line %s' %(L)) ; print(s.warn[-1])
                #res['LowEnergy'] = {'DelRho':float(w[0]), 'gmuminus2':float(w[1]), 'bsgam':float(w[2])}
                #res['LowEnergy'] = {'DelRho':w[0], 'gmuminus2':w[1], 'bsgam':w[2]}
                res['LowEnergy'] = {'DelRho':w[0], 'amuSUSYadd':w[1], 'bsgam':w[2]}
                # NB: what is given is not (amu=)gmuminus2/2, but the addition to the SM of (amu=)gmuminus2/2
                # exp(amu) = 0.00116592089 with error 6.3e-10
                # The amuSUSYadd is to be compared with the error (6.3e-10)

            if L.startswith('Fine-tuning values'):
                iL += 2
                L = lines[iL].rstrip()
                w = L.split()
                #res['Finetuning'] = {'dmZ(mu^2)':float(w[0]), 'dmZ(B.mu)':float(w[1]), 'dmt(mu)':float(w[2]), 'dmt(B.mu)':float(w[3])}
                #res['Finetuning'] = {'dmZ(mu^2)':float(w[0]), 'dmZ(B.mu)':float(w[1]), 'dmt(mu)':float(w[2]), 'dmt(B.mu)':float(w[3])}
                res['Finetuning'] = {'dmZ(mu^2)':w[0], 'dmZ(B.mu)':w[1], 'dmt(mu)':w[2], 'dmt(B.mu)':w[3]}


            if L.startswith('Warning/Error'):
                iL += 2
                L = lines[iL].rstrip()
                L = L.replace('-',' -')
                w = L.split()
                testsString = ''
                ErrTest = ''
                #testsInt = 0
                for iw in range(len(w)):
                    testsString += '%i ' %(abs(int(float(w[iw]))))
                    ErrTest += '%i' %(abs(int(float(w[iw]))))
                    #testsInt += 10**(9-iw) * abs(int(w[iw]))
                    res['E%i' %(iw+1)] = abs(int(float(w[iw])))
                testsString = testsString.rstrip()
                testsArr = []
                for iw in range(len(w)):
                    if w[iw] != 0: testsArr.append(suspect_errorcode[iw+1])
                res['Warning'] = {'ErrTest':ErrTest, 'testsString':testsString, 'testsArr':testsArr}
                


        # Return
        return res


        
    # ##########
    def Main(s):
        #if s.VB: print "INFO::%s  Main" %(s.myname)

        if s.fn: s.ID = s.fn.replace('_suspect2.out','')
        
        res = s.ReadSuspectOut(s.fn)

        if s.fn:
            fn_slha = s.ID + '_slhaspectrum.in'
            res['SPINFO'] = s.ReadSuspectSPINFO(fn_slha)


        # Make resflat
        resflat = {}
        for key in res:
            if type(res[key]) is dict:
                for key2 in res[key]: resflat[key2] = res[key][key2]
            else:
                resflat[key] = res[key]


        # Now make output 
        vars_finetune = ['dmZ(mu^2)','dmZ(B.mu)','dmt(mu)','dmt(B.mu)']
        #vars_lowenergy = ['DelRho','gmuminus2','bsgam']
        vars_lowenergy = ['DelRho','amuSUSYadd','bsgam']
        #vars_err = ['Err','ErrTest']
        #vars_err = ['Err','ErrTest','E1','E2','E3','E4','E5','E6','E7','E8','E9','E10']
        vars_err = ['Err','E1','E2','E3','E4','E5','E6','E7','E8','E9','E10']  # 2014-03-12

        head = s.prevar 
        for var in vars_err + vars_finetune + vars_lowenergy: head += '  ' + var
        out = s.preval
        for var in vars_err + vars_finetune + vars_lowenergy: out += '  ' + str(resflat[var])

        # Add some relative values
        head += '  rhoVS1sig'
        if float(resflat['DelRho']) > 0: out += '  %7.3f' %( float(resflat['DelRho']) / exp['rhoErrPlus'])
        elif float(resflat['DelRho']) <= 0: out += '  %7.3f' %( float(resflat['DelRho']) / exp['rhoErrMinus'])
        
        head += '  amuVS1sig'
        #out += '  %7.3f' %(float(resflat['gmuminus2'])/exp['amuErr'])
        out += '  %7.3f' %(float(resflat['amuSUSYadd'])/exp['amuErr'])
        
        head += '  bsgamVS1sig'
        out += '  %7.3f' %( (float(resflat['bsgam']) - exp['bsgam']) / exp['bsgamErr'] )
        
        outs = [head.strip(),out.strip()]
        outs = lineup(outs, n=s.nspaces)

        fn_out = s.ID + '_errlowtune.tlt'
        if not s.dryrun: WriteToFile(fn=fn_out, outs=outs, VB=s.VB-1)
        else:
            for out in outs: print(out)
        
        if s.pickle and not s.dryrun: 
            fn_pickle=s.ID + '_errlowtune.pickle'
            WritePickle(fn=fn_pickle, thepickle=res, VB=s.VB-1)

        

    # ##########
    def ReadArg(s): 

        # ################################### ARGUMENT READING
        Arg = ArgReader(s.argv, VB=0)

        '''
        if Arg.hasget('-alist'):  print 'a string list: ',Arg.list()
        if Arg.hasget('-alisti'): print 'an integer list: ',Arg.listI()
        if Arg.hasget('-alistf'): print 'a float list: ',Arg.listF()
        if Arg.hasget('-x'):  print 'a string: ',Arg.val()
        if Arg.hasget('-xI'): print 'an integer: ',Arg.valI()
        if Arg.hasget('-xF'): print 'a float: ',Arg.valF()
        '''

        if Arg.has(['-h','--help','--h','-help']):
            s.showHelp()
            sys.exit()

        if Arg.hasget('-vb'):
            s.VB = Arg.valI()
            if s.VB: print('Verbosity level: %i' %(s.VB))

        if Arg.hasget(['-f']):
            s.fn = Arg.val()

        if Arg.hasget(['-ID']):
            s.fn = Arg.val() + '_suspect2.out'  # use to make s.fn; later ID is then (re)derived ..

        if Arg.has(['--pickle']):
            s.pickle = 1

        if Arg.hasget(['-prevar']):  # -prevar 'M1 M2 M3' 
            s.prevar = Arg.val()
            
        if Arg.hasget(['-preval']):  # -preval '100 150 150'   # preval and prevar should match, a space will be inserted before the regular 2-liner
            s.preval = Arg.val()
            
        if Arg.hasget(['-n','-nspaces']):  # passed on to lineup
            s.nspaces = Arg.valI()

        if Arg.has(['--dryrun']):
            s.dryrun = 1

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
    t = suspect2errlowtune(cmd=['ReadArg','PostInit','Main'])
############################## 

