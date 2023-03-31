#!/usr/bin/env python
'''
Program : isasusy_status.py
Version : 1.0
Author  : b.k.gjelsten@fys.uio.no
Description : 

'''

import sys,os
from bkgjelstenArgReader import ArgReader
from lineup import lineup
from kilelib import WriteToFile


# ##################################################### GLOBAL METHODS
# ##################################################### GLOBAL METHODS
# ##################################################### GLOBAL METHODS


# ##########
def isaout_investigate(fn_out, prevar='', preval='', fn_tlt='', VB=0): 
    # pars isaout, investigate, put findings into resdict, including outs_tlt

    res = {}
    vars = ['warn','badconv','sssave','ssxintUNK','warnUNK']
    for var in vars: res[var] = 0

    #res[''] = 0
    #res[''] = 0
    #res[''] = 0
    
    f = open(fn_out); lines = f.readlines(); f.close()

    for iL in range(len(lines)):
        L = lines[iL].rstrip()

        if L.startswith(' WARNING'):
            res['warn'] += 1

            # ###
            if L.startswith(' WARNING IN SSXINT'):

                if L.startswith(' WARNING IN SSXINT: BAD CONVERGENCE'):
                    res['badconv'] += 1
                #elif: ...
                else:
                    res['ssxintUNK'] += 1



            elif L.startswith(' WARNING: SSSAVE'):
                res['sssave'] += 1
                res['warn'] +- 1    # non-important warning



            else: 
                res['warnUNK']


    res['head'] = prevar
    for var in vars: res['head'] += '  ' + var
    res['head'] = res['head'].strip()
    res['line'] = preval
    for var in vars: res['line'] += '  %s' %(res[var])
    res['line'] = res['line'].strip()
    res['tlt_outs'] = lineup([res['head'], res['line']])

    if fn_tlt:
        WriteToFile(fn=fn_tlt, outs=res['tlt_outs'], VB=VB-1)

    return res


# ##########
class isasusy_status: 

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


        s.prevar = optD.get('prevar','')
        s.preval = optD.get('preval','')

        s.fn_out  = ''
        s.fn_tlt = ''


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
        print('        %s  -dict test,txt1:I,testI,3:test2,txt2a,txt2b:F,testF,4.14   # for using autodict (NB: vars need to be defined in __init__)' %(s.myname))


    # ##########
    def DumpWarnings(s):
        f = open(s.fn_warn,'w')
        for out in s.warn: f.write('%s\n' %(out))
        f.close()
        
    # ##########
    def Main(s):
        if s.VB>1: print("INFO::%s  Main" %(s.myname))
        # for key in s.dict.keys(): print 'dict:  %-10s  %s' %(key, s.dict[key])

        if s.fn_out  == '': sys.exit('Fatal::isasusy_status  Need to give isasusy-output filename (-fn_out <fname>)')
        #if s.fn_tlt == '': sys.exit('Fatal::isasusy_status  Need to give isasusy-output filename (-fn_tlt <fname>)')

        if s.fn_tlt == '':
            s.fn_tlt = s.fn_out
            for zend in ['.out','_out']:
                if s.fn_tlt.endswith(zend): s.fn_tlt = s.fn_tlt[:-len(zend)]
            s.fn_tlt += '_status.tlt'

        # parse isasusy.out file: pick up any irregularities, produce tlt
        s.res = isaout_investigate(fn_out=s.fn_out, prevar=s.prevar, preval=s.preval, fn_tlt=s.fn_tlt)
        
        

        


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


        if Arg.hasget(['-f','-fn_out']):  # -prevar 'M1 M2 M3' 
            s.fn_out = Arg.val()

        if Arg.hasget(['-save','-fn_tlt']):  # -prevar 'M1 M2 M3' 
            s.fn_tlt = Arg.val()

        if Arg.hasget(['-prevar']):  # -prevar 'M1 M2 M3' 
            s.prevar = Arg.val()
            
        if Arg.hasget(['-preval']):  # -preval '100 150 150' 
            s.preval = Arg.val()


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
    t = isasusy_status(cmd=['ReadArg','PostInit','Main'])
############################## 

