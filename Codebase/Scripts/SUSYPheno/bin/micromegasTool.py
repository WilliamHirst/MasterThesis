#!/usr/bin/env python
'''
Program : micromegasTool.py
Version : 1.0
Author  : b.k.gjelsten@fys.uio.no
Description : 

'''

import sys,os
from bkgjelstenArgReader import ArgReader
from libExpValues import exp
from kilelib import WriteToFile, WritePickle, GetLinuxDist
from lineup import lineup
import platform

# ##################################################### GLOBAL METHODS
# ##################################################### GLOBAL METHODS
# ##################################################### GLOBAL METHODS
def ReadMicromegasOutput(lines, fn='', VB=1, opt=[]):

    blocks = {
        "mSUGRA scenario":"SUGRA"
        ,"mSUGRA non-universal Higgs scenario":"SUGRANUH"
        ,"AMSB scenario":"AMBS"
        ,"EWSB scale input":"EWSB"
        ,"SLHA file input":"SLHA"
        ,"MASSES OF HIGGS AND SUSY PARTICLES:":"MASS"
        ,"Physical Constraints:":"CONSTRAINTS"
        ,"Calculation of relic density":"RELICDENSITY"
        ,"Indirect detection":"INDIRECT"
        ,"Calculation of CDM-nucleons amplitudes":"CDM-nucleon amplitudes"
        ,"CDM-nucleon cross sections[pb]":"CDM-nuclean XSECS"
        ,"Direct Detection":"DIRECT"
        ,"Neutrino Telescope":"NEUTRINOTELESCOPE"
        ,"Decays":"DECAYS"
        ,"Calculation of cross section":"XSECS"
        }


    if fn: txtfn = '  [fn=%s] ' %(fn) 

    res = {}
    inBlock = ''
    res['warn'] = 0
    res['warnmsg'] = []
    nL = len(lines)
    

    # Loop through lines
    iL = -1
    #for iL in range(len(lines)):
    while iL < nL-1:
        iL += 1 ; L = lines[iL].strip() ; w = L.split()

        if L == '': continue

        # New block
        #print("L = ",L)
        if L.startswith("====") and ( L.endswith("====") or 'Neutrino Telescope' in L):   # inelegant, but works
            inBlock = ''
            for block in blocks: # loops through keys
                #print("block = ",block)
                if block in L:
                    inBlock = blocks[block]
                    res[inBlock] = {}
                    break

            if inBlock == '': # block not recognised
                res['warnmsg'].append("UNKNOWN BLOCK %s (continuing but things may fail): %s" %(txtfn, L))
                res['warn'] = 1
                if VB: print(res['warnmsg'][-1])
                inBlock = 'UNKNOWN'
                if inBlock not in res: res[inBlock] = {}

        bres = res[inBlock]  # 

        if VB>1: print(' %-20s  %3i  %s' %(inBlock, iL, L))

        # #####
        if inBlock == 'SLHA': 
            if L == "Warnings from spectrum calculator:":
                iL += 1 ; L = lines[iL].strip() ; w = L.split()
                while L != '':  # loop to pick up all warnings, but skip '.....none'
                    if L != '.....none': 
                        if 'warnmsg' not in bres: bres['warnmsg'] = []
                        bres['warnmsg'].append(L)
                        res['warn'] = 1
                    iL += 1 ; L = lines[iL].strip() ; w = L.split()
                # at loop end, L is ''

            if L.startswith('Dark matter candidate'):
                bres['LSP_MO'] = w[4].replace("'","")
                iL += 1 ; L = lines[iL].strip() ; w = L.split()
                iL += 1 ; L = lines[iL].strip() ; w = L.split()
                bres['LSP_bino'] = float(w[2].replace('*bino',''))
                bres['LSP_wino'] = float(w[3].replace('*wino',''))
                bres['LSP_higgsino1'] = float(w[4].replace('*higgsino1',''))
                bres['LSP_higgsino2'] = float(w[5].replace('*higgsino2',''))
                bres['LSP_b'] = bres['LSP_bino']**2
                bres['LSP_w'] = bres['LSP_wino']**2
                bres['LSP_h1'] = bres['LSP_higgsino1']**2
                bres['LSP_h2'] = bres['LSP_higgsino2']**2
                bres['LSP_h'] = bres['LSP_h1'] + bres['LSP_h2']
                bres['LSP_g'] = bres['LSP_b'] + bres['LSP_w']
        

        # #####
        elif inBlock == 'MASS': 
            continue

        
        # #####
        elif inBlock == 'CONSTRAINTS':
            if L.startswith('deltartho='):
                bres['DelRhoMO'] = w[0].replace('deltartho=','')
                if   float(bres['DelRhoMO'])  > 0: bres['rhoVS1sigMO'] = float(bres['DelRhoMO']) / exp['rhoErrPlus']
                elif float(bres['DelRhoMO']) <= 0: bres['rhoVS1sigMO'] = float(bres['DelRhoMO']) / exp['rhoErrMinus']
                
            if L.startswith('gmuon='):
                bres['amuSUSYaddMO'] = w[0].replace('gmuon=','')  # think
                bres['amuVS1sigMO'] = float(bres['amuSUSYaddMO']) / exp['amuErr']
                
            if L.startswith('bsmumu='): bres['bsmumu'] = w[0].replace('bsmumu=','')
            if L.startswith('btaunu'): bres['btaunu'] = w[0].replace('btaunu=','')
            if L.startswith('Rl23='): bres['Rl23'] = w[0].replace('Rl23=','')
            
            if L.startswith('bsgnlo='):
                bres['bsgamMO'] = w[0].replace('bsgnlo=','')
                bres['bsgamSMMO'] = w[3]
                bres['bsgamVS1sigMO'] = ( float(bres['bsgamMO']) - float(bres['bsgamSMMO']) ) / exp['bsgamErr']  # additional

            if L.startswith('dtaunu='): bres['dtaunu'] = w[0].replace('dtaunu=','') ; bres['dmunu'] = w[1].replace('dmunu=','')

            if L.startswith('MassLimits OK'): bres['MassLimitsOK'] = 1
            if L.startswith('WARNING:'):
                if 'warnmsg' not in bres: bres['warnmsg'] = []
                bres['warnmsg'].append(L)
                res['warn'] = 1


        # #####
        elif inBlock == 'RELICDENSITY':
            if L.startswith('Xf='):
                bres['Xf'] = w[0].replace('Xf=','')
                bres['Omega'] = w[1].replace('Omega=','')
                bres['OmegaVS1sig'] = ( float(bres['Omega']) - exp['cdm'] ) / exp['cdm_error']
                bres['OmegaVSexp'] = float(bres['Omega']) / exp['cdm']
                
            if L.startswith('# Relative contributions in % are displayed'): 
                iL += 1 ; L = lines[iL].strip() ; w = L.split()
                # Loop over relative contributions
                #   lines have form: 13% ~o1 ~o1 ->t T
                
                while (iL < nL-1) and (len(w) == 5) and ('%' in L) and ('->' in L):   # FRAGILE? 
                    w2 = L.replace('%','').replace('->','').strip().split()
                    val, in1_OM, in2_OM, out1_OM, out2_OM = w2
                    key = 'DM:%s_%s__%s_%s' %(in1_OM, in2_OM, out1_OM, out2_OM)
                    bres[key] = val
                    iL += 1 ; L = lines[iL].strip() ; w = L.split()
                

        # #####
        else:
            # not picking up these blocks for now
            pass
            

        #print 'res:  ', res
        #print 'bres: ', bres
        #print
        
    # #####
    # POST: remove empty blocks
    keys = list(res.keys())
    for key in keys: 
        if type(res[key]) is dict and res[key] == {}: res.pop(key)


    return res



# #########
class micromegasTool: 

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
        s.d = s.dict

        # s.dict['test'] = 'testval'
        # s.dict['test2'] = 'testval2'
        # s.dict['testI'] = 3
        # s.dict['testF'] = 4.34

        s.warn = []
        s.fn_warn = 'warnings.txt'
        s.fn_report = 'report'
        s.report = []

        
        s.d['dist'] = GetLinuxDist()   # 2014-08-25

        # fixed
        s.d['exes'] = ['micromegas_mini_slha','micromegas_orig_slha','micromegas_mini_ewsb','micromegas_orig_ewsb']
        s.d['exe'] = 'micromegas_mini_slha'  # this is default executable: relic density + SM constraints

        
        s.d['fns'] = []   # input files (piped), can be slha or ewsb

        s.preVar = ''
        s.preVal = ''
        

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

        if not 'SUSYPHENO_PATH' in list(os.environ.keys()):
            print("$SUSYPHENO_PATH not set. Exiting.")
            sys.exit()
        susyphenopath = os.environ['SUSYPHENO_PATH']

        s.d['dir_exe'] = susyphenopath+'/MICROMEGA/micromegas_3.5.5/MSSM'
        #if   s.d['dist'] in ['RHEL5']: s.d['dir_exe'] = '/net/abel-evs/cargo/fysepf/epfshare/prog/micromegas_3.5.5/MSSM'
        #elif s.d['dist'] in ['RHEL6']: s.d['dir_exe'] = '/net/abel-evs/cargo/fysepf/epfshare/prog/micromegas_3.5.5_RHEL6//micromegas_3.5.5/MSSM'   # inelegant structure
        #else:
        #    s.d['dir_exe'] = '/net/abel-evs/cargo/fysepf/epfshare/prog/micromegas/micromegas_3.5.5/MSSM'  # default
        #    print 'Warning  micromegasTool  non-recognised distribution: %s  (Proceeding as if RHEL5)' %(s.d['dist'])


    # ##################################################### CLASS METHODS
    # ##################################################### CLASS METHODS
    # ##################################################### CLASS METHODS
    # ##########
    def showHelp(s):
        print(' Usage: ls *.slha | micromegasTool.py')
        print('    (Note: some options are available in s.dict, more to be added.)')


    # ##########
    def DumpWarnings(s):
        f = open(s.fn_warn,'w')
        for out in s.warn: f.write('%s\n' %(out))
        f.close()
        
    # ##########
    def Main(s):
        #if s.VB: print "INFO::%s  Main" %(s.myname)
        # for key in s.dict.keys(): print 'dict:  %-10s  %s' %(key, s.dict[key])


        # Get list of input files
        if s.VB>0: print("INFO::micromegasTool   Now reading list of input files from stdin")
        s.d['fns'] = sys.stdin.readlines()


        # Loop over slha files
        for iF in range(len(s.d['fns'])):

            fn = s.d['fns'][iF].strip()
            fn_base = fn.replace('.slha','').replace('.ewsb','')

            if s.VB>1: print('DEBUG::micromegasTool  %i/%i  %s' %(iF, len(s.d['fns']), fn))
            
            # Calculate
            cmd = "%s/%s  %s" %(s.d['dir_exe'], s.d['exe'], fn)
            lines = os.popen(cmd).readlines()  # <--- FRAGILE ... should make safer (against errors)


            res = ReadMicromegasOutput(lines)
            #print '\nres\n', res

            resflat = s.FlatDict(res)
            #print '\nresflat\n', resflat

            s.SaveOutput(resflat, fn_base)



    # ##########
    def FlatDict(s, res):
        # Make flat dict
        resflat = {}
        for key in res:
            if not type(res[key]) is dict:
                if type(res[key]) is list: resflat[key] = list(res[key])
                else: resflat[key] = res[key]
            else:
                for key2 in res[key]:
                    if type(res[key][key2]) is list:
                        if key2 not in resflat: resflat[key2] = []
                        resflat[key2] += list(res[key][key2])   # lists will be added, e.g. overall warn and warn inside a block
                    else: resflat[key2] = res[key][key2]
        return resflat
    
            
    # ##########
    def SaveOutput(s, resflat, fn_base, opt=[]):

        formats = {'warn':'%i', 'OmegaVSexp':'%9.3e', 'OmegaVS1sig':'%7.3f', 'rhoVS1sigMO':'%7.3f', 'amuVS1sigMO':'%7.3f','bsgamVS1sigMO':'%7.3f'}

        # mini-tlt
        if 1: 
            outs = []
            vars = ['warn','Omega','OmegaVSexp','OmegaVS1sig', 'bsmumu','btaunu','Rl23','dtaunu', 'dmunu', 'DelRhoMO','amuSUSYaddMO','bsgamMO', 'rhoVS1sigMO', 'amuVS1sigMO','bsgamVS1sigMO']
            head = s.preVar
            out =  s.preVal
            for var in vars:
                head += '  %s' %(var)
                if var not in resflat:
                    out += '  -'
                    print("Warning::micromegasTool  SaveOutput  var %s not in resflat. Set to '-'" %(var))
                else:
                    if type(resflat[var]) is str:
                        out += '  %s' %(resflat[var])
                    else:
                        format = formats.get(var,'%8.2e')
                        out += '  ' + format %(resflat[var]) 

            outs = lineup([head.strip(), out.strip()])

            fn_mini_tlt = '%s_MO_mini.tlt' %(fn_base)
            WriteToFile(fn=fn_mini_tlt, outs=outs, VB=s.VB-2)


        # dict


        # flatdict




            


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

        if Arg.hasget('-dist'):
            s.d['dist'] = Arg.val()  # allowed 'RHEL5', 'RHEL6'

        if Arg.hasget('-prevar'): 
            s.preVar = Arg.val()
            
        if Arg.hasget('-preval'): 
            s.preVal = Arg.val()
            


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
    t = micromegasTool(cmd=['ReadArg','PostInit','Main'])
############################## 

