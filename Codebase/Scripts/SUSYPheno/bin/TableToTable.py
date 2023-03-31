#!/usr/bin/env python
'''
Program : TableToTable.py
Version : 1.0
Author  : b.k.gjelsten@fys.uio.no
Description : 
 Input: text table with header
 Output: text table with a selection of the quantities in the input table and arithemtic combinations thereof

Usage: from terminal, or imported

'''

import sys,os
import bkgjelstenArgReader 
from kilelib import DictFromTable, GetPlainArray
#from TableToTGraph import GetPlainArray    # will move?
from lineup import lineup


# ##################################################### GLOBAL METHODS
# ##################################################### GLOBAL METHODS
# ##################################################### GLOBAL METHODS




class TableToTable: 

    def __init__(s, cmd=['ReadArg','PostInit','Main'], optD={}):

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

        s.vars  = []
        s.varsT = {}
        s.arraytype = ""  # empty:return list; otherwise return an array.array of given type (typical 'f')

        s.show = 1

        s.format = {'default':"%.2f"}
        s.vars_addoriginal = 1  # add original variables or not

        s.delim_vars = ','
        s.delim_varsT = ':'
        s.fn_in = ''

        s.res = {}


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
        print(' Usage:  cat table.txt | TableToTable.py  [options]')
        #print '        %s  -dict test,txt1:I,testI,3:test2,txt2a,txt2b:F,testF,4.14   # for using autodict (NB: vars need to be defined in __init__)' %(s.myname)
        print('    Ex:  cat table.txt | TableToTable.py   # does nothing (except changes output format)')
        print("    Ex:  cat table.txt | TableToTable.py  -vars 'M2,N2->N1,h:N2_N1_h,(N2-N1)/125.5:relmassdiff' ")
        print('    Ex:  cat table.txt | TableToTable.py')
        print("    Ex:  TableToTable.py  -f DumpScenario.txt -delim_vars ,,   -vars 'M2,,N2->N1,h:N2_N1_h,,N2->N1,Z:N2_N1_Z,,N2->N1,h+N2->N1,Z:N2_N1_horZ'  -vb 0  -format M2:%.1f,,N2_N1_h:%.3f")
        print() 
        print('  The table.txt needs to have the quantity names on the first line')
        print() 
        print('  Options:')
        print('    -vars : which vars or combinations to show.')
        print('            All of +,-,*,/,() can be used in the variable combination, e.g. (N2+N1)/2')
        print('            Quantities can be renamed by adding the new label after a column')
        print("            If the name of a quantity contains the key-characters ',' or ':', these will need to be changed by specifying ex. one/both of -delim_vars ,,  -delim_varsT ::")
        print("            Note: if you need to put vars inside '', e.g. because they contain '>', use an overall '' for the entire list, ex.:  -vars 'M2,N2->N1,h'   (do NOT write -vars M2,'N2->N1,h') ")
        print('    -addvars: same as vars, but the original variables are also kept')
        print('    -show 0/1 : print the output table to screen. (If routine accessed within API, it may make sense not to output to screen')
        print("    -format <key>:<format1>,<key2>:<format2>,...  : specifies format for output (using the new var labels). If unspecified, the default is '%.2f' Changes in delimiter will apply here as well.")
        print("            Ex: -format M2:%.1f,,N2_N1_h:%.3f")
        print('    --silent : no info except the resulting table (if show is on)')
        print('    -vb <0,1,..> : set verbosity level')

        print()
        print('  TableToTable can also be imported into other Python scripts:')
        print('    import TableToTable')
        print("    t = TableToTable.TableToTable(optD={'argv':\"-delim_vars ,,   -vars 'M2,,N2->N1,h:N2_N1_h,,N2->N1,Z:N2_N1_Z,,N2->N1,h+N2->N1,Z:N2_N1_horZ'  -vb 0  -format M2:%.1f,,N2_N1_h:%.3f  -f DumpScenario.txt  -show 0\".split()})")
        print('    Then access results from the dict t.res')
        print()
        

    # ##########
    def DumpWarnings(s):
        f = open(s.fn_warn,'w')
        for out in s.warn: f.write('%s\n' %(out))
        f.close()
        
    # ##########
    def Main(s):
        # if s.VB: print "INFO::%s  Main" %(s.myname)

        # 1 Read input table
        if s.fn_in: 
            s.res['tableText'] = os.popen('cat %s' %(s.fn_in)).readlines()
        else:
            if s.VB>0: print("Getting table from pipe ...")
            s.res['tableText'] = sys.stdin.readlines()
        s.res['tableDict'], s.res['tableHead'] = DictFromTable(s.res['tableText'], returnhead=1)

        # 2 Opt: add original vars
        if s.vars_addoriginal:
            s.vars = s.res['tableHead'] + s.vars
        if s.VB>1: print('DEBUG::vars: %s' %(s.vars))
        

        # Make the dict of resulting values
        s.res['resval'] = {}
        nLines = 0
        for ivar in range(len(s.vars)):
            var = s.vars[ivar]
            varT = s.varsT.get(var, var)
            if s.VB>1: print('DEBUG::loop ', ivar, var, varT)
            s.res['resval'][varT] = GetPlainArray(table=s.res['tableDict'], var=var, arraytype=s.arraytype, VB=s.VB)
            if s.VB>1: print('DEBUG::loop nvals: ', len(s.res['resval'][varT]))
            if s.VB>2: print('DEBUG::loop vals: ', s.res['resval'][varT])
            nLines = len(s.res['resval'][varT])
            #print 'nLines: %i' %(nLines)

        # Make text table
        s.res['outs'] = []
        for iL in range(nLines):
            # The header line
            if iL == 0:
                out = ''
                for var in s.vars: out += '  ' + s.varsT.get(var,var)
                s.res['outs'].append(out.strip())
            # The standard line
            out = ''
            for var in s.vars:
                varT = s.varsT.get(var,var)
                
                out += '  ' + s.format.get(varT,s.format.get('default')) %(s.res['resval'][varT][iL])
            s.res['outs'].append(out)
        s.res['outs'] = lineup(s.res['outs'])


        # Dump to screen? 
        if s.show:
            for out in s.res['outs']: print(out)



    # ##########
    def ReadArg(s): 

        # ################################### ARGUMENT READING
        Arg = bkgjelstenArgReader.ArgReader(s.argv, VB=0)

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

        if Arg.hasget(['-show']): s.show = Arg.valI()

        if Arg.has(['--silent']): s.VB = 0

        if Arg.hasget(['-delim_vars']): s.delim_vars = Arg.val()
        if Arg.hasget(['-delim_varsT']): s.delim_vars = Arg.val()
        if Arg.hasget(['-delim']): s.delim_vars, s.delim_varsT = Arg.list()  # (limited)
        

        if Arg.hasget(['-vars']):  # Ex: -vars N1,N2:N20,N3,N2-N1:N2mN1,N2/N1
            #zz = Arg.list(s.delim_vars)
            zz = Arg.val()  #.split(s.delim_vars)
            #print 'AAA', zz
            #print 'BBB', str(zz)
            zz = zz.strip("'").split(s.delim_vars)  # NB: hack to avoid getting a string like "'kdkfds....fsad'"
            #print 'DEBUG vars: ', zz
            
            for z in zz:
                w = z.split(s.delim_varsT)
                zvar = w[0]
                if len(w)>1: zvarT = w[1]
                else: zvarT = w[0]
                s.vars.append(zvar)
                s.varsT[zvar] = zvarT
            s.vars_addoriginal = 0


        if Arg.hasget(['-addvars']):  # Ex: -vars N1,N2:N20,N3,N2-N1:N2mN1,N2/N1
            #zz = Arg.list(s.delim_vars)
            zz = Arg.val()
            #print 'AAA', zz
            #print 'BBB', str(zz)
            zz = zz.strip("'").split(s.delim_vars)
            for z in zz:
                w = z.split(s.delim_varsT)
                zvar = w[0]
                if len(w)>1: zvarT = w[1]
                else: zvarT = w[0]
                s.vars.append(zvar)
                s.varsT[zvar] = zvarT
            s.vars_addoriginal = 1


        if Arg.hasget(['-format']):
            zz = Arg.list(s.delim_vars)
            for z in zz:
                w = z.split(s.delim_varsT)
                s.format[w[0]] = w[1]


        if Arg.hasget(['-f','-fn','-fn_in']):
            s.fn_in = Arg.val()
            

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
    t = TableToTable(cmd=['ReadArg','PostInit','Main'])
############################## 

