#!/usr/bin/env python
'''
Program : softsusy_tool.py
Version : 1.0
Author  : b.k.gjelsten@fys.uio.no
Description :
Just a simple wrapper around standalone softsusy to at least create a tlt to show overall status (w.r.t. failure or not)

'''

import sys,os,random
from bkgjelstenArgReader import ArgReader
from kilelib import WriteToFile, readfile_strip
from lineup import lineup

import slhafrompardict

# ##################################################### GLOBAL METHODS
# ##################################################### GLOBAL METHODS
# ##################################################### GLOBAL METHODS

# ##########
#def softsusy_execute(fn_base_out='', fn_slha_in='', SUSYpar={}, par_steer={}, fn_add='_softsusy'):
def softsusy_execute(optD): 
    #obj_softsusy = softsusy_tool(cmd=['PostInit','Main'], optD={'fn_base_out':fn_base_out, 'fn_slha_in':fn_slha_in, 'SUSYpar':SUSYpar, 'par_steer':par_steer, 'fn_add':fn_add})
    obj_softsusy = softsusy_tool(cmd=['PostInit','Main'], optD=optD)
    del obj_softsusy
    
    

# ##########
class softsusy_tool: 

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

        s.SUSYpar = optD.get('SUSYpar',{})   # this includes softsusy steering vars
        s.par_steer = optD.get('par_steer',{})

        s.dict = {}

        # s.dict['test'] = 'testval'
        # s.dict['test2'] = 'testval2'
        # s.dict['testI'] = 3
        # s.dict['testF'] = 4.34

        s.warn = []
        s.fn_warn = 'warnings.txt'
        s.fn_report = 'report'
        s.report = []

        s.dict['exe'] = 'softpoint_3.4.0.x'
        s.fn_slha_in  = optD.get('fn_slha_in','')
        s.fn_base_out = optD.get('fn_base_out','')
        #s.fn_slha_out = ''

        s.dict['fn_add'] = optD.get('fn_add','_softsusy')
        s.dict['slha_ending'] = '.slha'
        s.dict['keep_in2'] = 1

        s.fn_warn_global = '%s/.warning_softsusy_tool' %(s.HOME)

        s.prevar = optD.get('prevar','')
        s.preval = optD.get('preval','')
        
        
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


        # 1a Prepare input file for softsusy execution, prepare output filenames
        if s.fn_slha_in == '':
            if s.SUSYpar == {}:     sys.exit('Fatal::softsusy_tool with empty fn_slha_in you need to provide SUSYpar')
            if s.fn_base_out == '': sys.exit('Fatal::softsusy_tool with empty fn_slha_in you need to provide fn_base_out')
            outs_slha_in2 = slhafrompardict.slhafrompardict(par1=s.SUSYpar, par2=s.par_steer, calculator=['softsusy'])

        else:
            outs_slha_in2 = readfile_strip(s.fn_slha_in, opt=['rstrip'])
            if s.par_steer:
                outs_slha_in2 += slhafrompardict.softsusyblockfrompardict(par=s.par_steer)
            if s.fn_base_out == '':
                s.fn_base_out = s.fn_slha_in.replace('.slha','').replace('_slha','').replace('.lha','') + s.dict['fn_add']

            
        # 1b filenames
        fn_slha_in2 = '%s_softin2.slha' %(s.fn_base_out)  # if to be kept
        fn_slha_intmp = '%s_softin2.slha_tmp%s' %(s.fn_base_out, str(random.randrange(0,1e5)).zfill(5))
        #fn_tlt = s.fn_base_out + '.tlt'
        fn_tlt = s.fn_base_out + '_status.tlt'  # 2014-03-12
        fn_slha_out = s.fn_base_out + s.dict['slha_ending']


        # 2 Write temporary input slha file
        WriteToFile(fn=fn_slha_intmp, outs=outs_slha_in2, VB=s.VB-1)

        
        # 3 Issue the command
        #  --- should here add option to add softsusy options, see. p.12-13 in manual
        cmd = "%s  leshouches  <  %s" %(s.dict['exe'], fn_slha_intmp)
        lines_slha = os.popen(cmd).readlines()

        

        # 4 Check the output: clean output, produce slha, produce tlt

        res = {}
        res['spect'] = 1
        
        lines_drop = []
        lines_warn = []
        warn_global = []
        for iL in range(len(lines_slha)):
            L = lines_slha[iL].rstrip()
            
            if L.startswith("# WARNING: don't recognise block"):
                lines_drop.append(iL)
            
            elif L.startswith("WARNING: did not understand parameter"):
                lines_drop.append(iL)
                lines_warn.append('### %s' %(L))  # add at the end? 

            elif L.startswith("# Declining to write spectrum because of serious problem with point"):
                res['spect'] = 0
            
            elif L.startswith("# SOFTSUSY problem with point"): 
                if "A0 tachyon" in L: res['tachA0'] = 1
                elif "stau tachyon" in L: res['tachT1'] = 1
                elif "selectron tachyon" in L: res['tachsel'] = 1
                elif "smuon tachyon" in L: res['tachsmu'] = 1
                # elif ...
                elif 'tachyon' in L:
                    res['tach'] = 1 
                    warn_global.append(L) # unrecognised problem ... should register somewhere
                else: 
                    res['probl'] = 1  
                    warn_global.append(L) # unrecognised problem ... should register somewhere

                lines_warn.append(L)

            elif L.startswith("# MIXING"):
                wz = L.replace('=',' ').split()
                res['MIX'] = int(wz[2])
                res['accuMax'] = float(wz[5])
                res['accu'] = float(wz[8])


        # 4b Remove lines
        for iL in reversed(lines_drop): lines_slha.pop(iL)



        # 5 Write resulting slha file
        WriteToFile(fn=fn_slha_out, outs=lines_slha, opt=['rstrip'], VB=s.VB-1)


        # 6 Make tlt
        head = s.prevar + " spect  tachA0  tachT1  tach  probl  MIX  accu  accuMax"
        out  = s.preval + " %i  %i  %i  %i  %i  %i  %8.2e  %8.2e" %(res['spect'], res.get('tachA0',0), res.get('tachT1',0), res.get('tach',0), res.get('probl',0), res.get('MIX',-1), res.get('accu',-1), res.get('accuMax',-1))
        outs_tlt = lineup(lines=[head,out])
        WriteToFile(fn=fn_tlt, outs=outs_tlt, VB=s.VB-1)


        # 7 Remove or rename intmp
        if s.dict['keep_in2']:
            cmd = 'mv %s %s' %(fn_slha_intmp, fn_slha_in2)
            os.system(cmd)

        
        # 8 Global warning (where I can check and improve this parsing)
        if warn_global: WriteToFile(fn=s.fn_warn_global, outs=warn_global, wORa='a', VB=s.VB)

        

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

        if Arg.hasget(['-fn_slha', '-fn_slha_in', '-fn_in', '-fn', '-f']):
            s.fn_slha_in = Arg.val()

        if Arg.hasget(['-fn_base_out']): 
            s.fn_base_out = Arg.val()

        if Arg.hasget(['-par','-SUSYpar']):
            w = Arg.list()
            for z in w:
                var,val = z.split('=')
                s.SUSYpar[var] = val  # note, is string

        if Arg.hasget(['-par_steer']):
            w = Arg.list()
            for z in w:
                var,val = z.split('=')
                s.par_steer[var] = val  # note, is string

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
    t = softsusy_tool(cmd=['ReadArg','PostInit','Main'])
############################## 

