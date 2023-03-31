#!/usr/bin/env python
'''
Program : slhadumper.py
Version : 1.0
Author  : b.k.gjelsten@fys.uio.no
Description : 


TODO:
o simplifications for some key vars (M1, ...) : should not need to know their integer key
o allow dumping BR (slha[1]) as well





'''

import sys,os
import bkgjelstenArgReader 

import pyslha_edited as pyslha

from lineup import lineup

from SlhaTools import pdg2nam, nam2pdg



# ##################################################### GLOBAL METHODS
# ##################################################### GLOBAL METHODS
# ##################################################### GLOBAL METHODS


class slhadumper: 

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

        s.fns_slha = ''
        s.fn_slhas = []

        s.fn_dump = ''

        s.format = {'mass':'%7.2f'}
        s.block_key_name = []

        s.delim1 = ','
        s.delim2 = ':'
        
        s.dict['screen'] = 1

        s.dict['positivemass'] = 1
        s.dict['dumppyslha'] = 0

        s.showfn = 0

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

        if s.fn_slhas: lines = os.popen('cat %s').readlines()
        else:
            if s.VB: print("INFO::slhadumper  Expecting table in by pipe     [if this is unexpected, you may want to break off and type 'slhadumper.py -help']")
            lines = sys.stdin.readlines()
        for L in lines: s.fn_slhas.append(L.strip())
        

        

    # ##################################################### CLASS METHODS
    # ##################################################### CLASS METHODS
    # ##################################################### CLASS METHODS
    # ##########
    def showHelp(s):
        print('DESCRIPTION')
        print('  Takes as input a list of slha files (piped in)')
        print('  and the variables to be shown, ')
        print('  then produces a two-line-table.')
        print('  The variables are specified with their BLOCK and slha key (and optionally their name)')
        print('    This is not so nice, though, as one needs to know the inner structure of slha.)')
        print('    For masses the key is the pdgID, but the sparticle names can be given instead')
        print('    To see the sparticle name convention:  slhadumper.py --showparticlenames')
        print("  To brute dump all the blocks execept branching ratios, add   '-dict I,dumppyslha,1' ")
        print("  To brute dump also the branching ratios, add                 '-dict I,dumppyslha,2' ")
        print(' ')
        print(' Note1: Making tables (tlt) of branching ratios is not yet implemented.')
        print(' Note2: There may be more advanced scripts to investigate slha files, e.g. slha2tlt.py')
        print() 
        print()
        print('EXAMPLES')
        print('   ls *.slha | slhadumper.py  -vars MASS:1000022:neutralino1,1000024:chargino1')
        print('   ls *.slha | slhadumper.py  -vars MASS:N1,C1  -format MASS:%.3f   # works with sparticle names instead of pdgID code, can set format of vars in MASS block')
        #print "   ls susy_DGnoSL_TB10_M1M2MU_050_500_*.slha | slhadumper.py  -vars EXTPAR:1:M1,2:M2,23:MU,MASS:N1:,N2,N3,N4,C1,C2  -format MASS:%.4f  -vb 0  # 2013-12-10  # slha from SUSYHIT, e.g. cd ~/grids_lsp/DGnoSL_TB10_Fine_2013-11-22/lhawig_M1eq050/"
        print("   ls *.slha | slhadumper.py  -vars EXTPAR:1:M1,2:M2,23:MU,MASS:N1:,N2,N3,N4,C1,C2  -format MASS:%.4f  -vb 0  # 2013-12-10  # slha from SUSYHIT")
        #print "   ls *slha |tail -20 | slhadumper.py  -vars MSOFT:1:M1,2:M2,HMIX:1:MU,MASS:N1:,N2,N3,N4,C1,C2  -format MASS:%.4f  -vb 0   # 2013-12-10  # self-produced slha (with MSOFT, but not EXTPAR), e.g. cd ~/grids_lsp/res_DGemtFine/dir_slha_samename"
        print("   ls *.slha |tail -20 | slhadumper.py  -vars MSOFT:1:M1,2:M2,HMIX:1:MU,MASS:N1:,N2,N3,N4,C1,C2  -format MASS:%.4f  -vb 0    # 2013-12-10  # self-produced slha (with MSOFT, but not EXTPAR)")
        print()
        print("OPTIONS  (see --h for full list (though without explanations)")
        print("   -vb <0/1/..>      : change verbosity (D=1)")
        print("   -h                : this help message")
        print("   --h               : brute info on implemented options")

    # ##########
    def DumpWarnings(s):
        f = open(s.fn_warn,'w')
        for out in s.warn: f.write('%s\n' %(out))
        f.close()

        
    # ##########
    def Main(s):
        #if s.VB: print "INFO::%s  Main" %(s.myname)
        # for key in s.dict.keys(): print 'dict:  %-10s  %s' %(key, s.dict[key])

        outs = []
        for islha in range(len(s.fn_slhas)):
            fn_slha = s.fn_slhas[islha]
            if s.VB>1: print('DEBUG::Main  %i %s' %(islha, fn_slha))

            slha = pyslha.readSLHAFile(fn_slha)
            # slha[0]: general part
            # slha[1]: BR part

            out = ''
            head = ''

            if s.showfn: out = fn_slha ; head = 'fn_slha'


            # ------ BRUTE DUMPING
            if s.dict['dumppyslha'] >= 1:
                print(120*'-')
                print('----- Dumping pyslha[0]: %s' %(fn_slha))
                for block in slha[0]:
                    #print 'BLOCK: %s' %(block)
                    print() ; print(slha[0][block])
                print(120*'-')
                print() 
            if s.dict['dumppyslha'] == 2:
                print('----- Dumping pyslha[1]:')
                #print slha[1]
                for p in slha[1]:
                    print()
                    print('p:    ', p)
                    print('obj:  ', slha[1][p])
                    print('mass: ', slha[1][p].mass)
                    try: print('widt: ', slha[1][p].total_width)
                    except: print('widt: ', 'no total_width')
                    #print 'keys: ', slha[1][p].keys()
                    #print 'entr: ', slha[1][p].entries
                print(120*'-')
                print() 

            if s.dict['dumppyslha'] == 3:  # slightly more structured way of dumping
                # slha[1] is a dict (in pdgID) of pyslha.Particle objects
                # these objects have attributes pid, mass, totalwidth which can be accessed as e.g. slha[1][1000023].mass
                #   and decays, which is a plain list of the possible decays 
                #   Each decay is of class pyslha.Decay with attributes br (a float) and ids (a list of pdgIDs)
                # So ... to access the BRs, we'll need something like
                print('----- Traversing and dumping pyslha[1]: ')
                for pdgID_parent in slha[1]:
                    print() 
                    parentT = str(pdg2nam.get(pdgID_parent,pdgID_parent))
                    parent = slha[1][pdgID_parent]
                    decays = parent.decays
                    for idec in range(len(decays)):
                        dec = decays[idec]
                        br = dec.br
                        childrenT = ''
                        for pdgID_ch in dec.ids:
                            pdgID_ch = abs(pdgID_ch)  # Note: selecting 
                            childrenT += '_' + str(pdg2nam.get(pdgID_ch,pdgID_ch))
                        childrenT = childrenT.lstrip('_')
                        decayT = parentT + '_' + childrenT 
                        print("%6.3f  %s" %(br, decayT))
                print('--- Note: particles/antiparticles are not distinguished in the above dump')
                print('---       This also means there will be two entries for some decays')
                
                

            # ------ MAKE TLT
            for block,key,name in s.block_key_name:
                #if s.VB>2: print 'DEBUG::Main  block,key,name:  %s,%s,%s' %(block,key,name)
                if block not in slha[0]:
                    s.warn.append('Warning   Skipping unknown block %s for %i %s  with blocks %s' %(block, islha, fn_slha, list(slha[0].keys())))
                    print(s.warn[-1])
                    continue
                
                obj = slha[0][block].entries
                val = obj.get(key,'-9')

                if s.VB > 2: print(block, key, name) #, obj, val
                
                form = s.format.get(block,'%.1f')
                if form.endswith('f'): val = float(val)
                if block in ['MASS'] and s.dict['positivemass']: val = abs(val)  # in case of negative mass parameters..
                out += ' ' + form %(val)

                if islha == 0: head += ' ' + name

            # 
            if islha == 0: outs.append(head)
            outs.append(out)

        outs = lineup(outs)

        if s.fn_dump: WriteToFile(fn=s.fn_dump, outs=outs, VB=s.VB-1)
        if s.dict['screen']:
            for out in outs: print(out)
            



    # ##########
    def ReadArg(s): 

        # ################################### ARGUMENT READING
        Arg = bkgjelstenArgReader.ArgReader(s.argv, VB=0)


        #if Arg.hasget('-alist'):  print 'a string list: ',Arg.list()
        #if Arg.hasget('-alisti'): print 'an integer list: ',Arg.listI()
        #if Arg.hasget('-alistf'): print 'a float list: ',Arg.listF()
        #if Arg.hasget('-x'):  print 'a string: ',Arg.val()
        #if Arg.hasget('-xI'): print 'an integer: ',Arg.valI()
        #if Arg.hasget('-xF'): print 'a float: ',Arg.valF()


        if Arg.has(['--h','---help','--automan']):
            os.system("cat %s | grep 'if Arg\.has'" %(s.argv[0]))
            print("[ The above shows all implemented options. For full details inspect the code, %s  (Try also '-h') ]" %(s.argv[0]))
            sys.exit()
            

        if Arg.has(['-h','--help','--h','-help']):
            s.showHelp()
            sys.exit()

        if Arg.hasget('-vb'):
            s.VB = Arg.valI()
            if s.VB: print('Verbosity level: %i' %(s.VB))

        if Arg.has('--showparticlenames'):
            print('Particle naming convention:')
            for nam in nam2pdg: print('  %-6s  %i' %(nam, nam2pdg[nam]))
            sys.exit()

        if Arg.hasget('-format'):
            zz = Arg.list()
            for z in zz:
                w = z.split(':')
                s.format[w[0]] = w[1]

        if Arg.hasget('-delim1'):
            s.delim1 = Arg.val()

        if Arg.hasget('-delim2'):
            s.delim2 = Arg.val()

        if Arg.has('--showfn'): s.showfn = 1
        if Arg.hasget('-showfn'): s.showfn = Arg.valI()
            
        if Arg.hasget('-vars'):
            zz = Arg.list(s.delim1)
            block = ''
            for z in zz:
                name = ''
                w = z.split(s.delim2)
                if len(w) == 1 and block == '': sys.exit('Fatal in arguments: %s' %(zz))
                elif len(w) == 1: key = w[0]
                elif len(w) == 2:
                    if block == '': block,key = w
                    else: key,name = w
                    #if name == '': name = key
                elif len(w) == 3: block,key,name = w
                else:
                    sys.exit('Fatal in arguments: %s' %(zz))

                if name == '': name = str(key)
                # Below: try twice, with and without adding '~'. This allows both N1 and ~N1 to be input
                key = nam2pdg.get('~'+key,key)  
                key = nam2pdg.get(key,key)      
                
                key = int(key)
                s.block_key_name.append([block,key,name])
            
                    

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
    t = slhadumper(cmd=['ReadArg','PostInit','Main'])
############################## 

