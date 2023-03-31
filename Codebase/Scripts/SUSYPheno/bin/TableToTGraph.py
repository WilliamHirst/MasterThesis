#!/usr/bin/env python
'''
Program : TableToTGraph.py
Author  : b.k.gjelsten@fys.uio.no
Version : 1.0   22.11.2013

Description : 


STATUS 2014-01-30
 - Has given some features to munch.py (which makes plots (.pdf)), which is then uptodate with the latest technologies
 - TableToTGraph.py is still the state-of-the-art TGraph-maker

 - Eventually need to see if it should become part of munch.py


'''

import sys,os, array
import ROOT
from bkgjelstenArgReader import ArgReader
from kilelib import DictFromTable
from kilelib import GetPlainArray

# ##################################################### GLOBAL METHODS
# ##################################################### GLOBAL METHODS
# ##################################################### GLOBAL METHODS



# ###################################################################
def MakeTGraph(table, coords, resvar, coordsT=[], resvarT='', VB=1):
    # table is a dict: table={x:[...], y:[...], N1:[...], ...}
    # Both coord and var (both arrays) can be complex entities, like x-y,2x
    # resvar is just one variable ; there is only one 'z'-var per TGraph (but the root file could contain a list of TGraphs) 

    # 1 Init
    if resvarT == '': resvarT = resvar
    for coord in coords: 
        if coord not in coordsT: coordsT[coord] = coord

    #print 'coordsT: ', coordsT

    # 2 Build plain arrays
    coordval = {}
    resval = []
    for coord in coords:
        coordval[coord] = GetPlainArray(table=table, var=coord, arraytype='f')
        if VB>1: print(('\n', coord, len(coordval[coord]), coordval[coord]))
    resval = GetPlainArray(table=table, var=resvar, arraytype='f')
    if VB>1: print(('\n', resvar, len(resval), resval))


    # 3 Make TGraphs
    if len(coords) == 1:
        tGraph = ROOT.TGraph  (len(resval), coordval[coords[0]], resval)
        print('here')
    elif len(coords) == 2:
        tGraph = ROOT.TGraph2D(len(resval), coordval[coords[0]], coordval[coords[1]], resval)
        if VB>2:
            for i in range(len(resval)): print(('%4i  %9.3f  %9.3f    %15.9f' %(i, coordval[coords[0]][i], coordval[coords[1]][i], resval[i])))
                                                                                                                                         
    else:
        tGraph = 'Error::MakeTGraph  (Returning this error message.)    Irregular coords: %s' %(coords)
        print(tGraph)


    # 4 Post treatment
    tGraph.SetName(resvarT)
    #print "SetName: %s" %(resvarT)
    tGraph.SetTitle(resvarT)
    tax = tGraph.GetXaxis()
    tax.SetTitle(coordsT[coords[0]])
    #print 'tax: ', coordsT[0]
    #print tGraph.GetXaxis().GetTitle()
    if len(coords) >= 2: 
        tay = tGraph.GetYaxis()
        tay.SetTitle(coordsT[coords[1]])
        
    if VB>=0: print(('INFO::MakeTGraph  Created TGraph (%iD)  %s with %i entries' %(len(coords), resvarT, len(resval))))


    # 9 Return
    # print 'ret: ', tGraph.GetTitle()
    return tGraph



# ###################################################################
# ###################################################################
# ###################################################################
class TableToTGraph: 

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


        s.coords = []
        s.coordsT = {}
        
        s.resvars = []
        s.resvarsT = {}

        s.fn_root = 'TableToTGraph.root'
        s.fn_table = ''

        s.table = {}
        
        s.dict['fnaddaxes'] = 1


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

        if s.fn_table: lines = os.popen('cat %s').readlines()
        else:
            if s.VB: print("INFO::TableToTGraph  Expecting table in by pipe     [if this is unexpected, you may want to break off and type 'TableToTGraph.py -help']")
            lines = sys.stdin.readlines()

        if len(lines) <= 2: print(('Warning::TableToTGraph  Table has only %i lines' %(len(lines))))

        s.table = DictFromTable(lines)

        # Hack to add coords to filename
        if s.dict['fnaddaxes'] and len(s.coords) == 2:
            z = '_x=%s_y=%s' %(s.coordsT[s.coords[0]], s.coordsT[s.coords[1]])
            s.fn_root = s.fn_root.replace('.root','%s.root' %(z))
    

    # ##################################################### CLASS METHODS
    # ##################################################### CLASS METHODS
    # ##################################################### CLASS METHODS
    # ##########
    def showHelp(s):
        print((' Usage: %s [options]' %(s.myname)))
        #print '        %s  -dict test,txt1:I,testI,3:test2,txt2a,txt2b:F,testF,4.14   # for using autodict (NB: vars need to be defined in __init__)' %(s.myname)
        print("    Ex: cat DGnoSL_masses.txt | TableToTGraph.py -coords MU,M2  -resvars 'N1,N2,N3,N4,C1,C2,h,N2-N1:N2mN1,C1-N1:C1mN1,N2-N1:N2mN1,(N2+N1)/2:N2plusN1on2,(h-2*N1):hfunnel'  -save DGnoSL_masses")
        print()
        print('    -dict I,fnaddaxes,0     # to NOT add the coordinates to the output filename')
        print('    -resvars: contains the variable(combination)s which get their TArray ')
        print('              Adding a column to a variable simply renames the TGraph ')
        print('     Note: there are problems in setting the titles on the axis, therefore the x and y variable are given specifically in the output root file.')
        print() 
        print('           where DGnoSL_masses.txt looks something like:')
        print('              M2      MU     N1      N2      N3      N4      C1      C2       h')
        print('           100.0   100.0   31.2    70.2   119.6   176.9    56.3   177.3   125.5')
        print('           100.0   120.0   36.0    74.5   138.3   186.5    65.3   186.9   125.5')
        print('           100.0   140.0   39.3    78.8   157.3   197.8    72.8   198.4   125.5')
        print('           100.0   150.0   40.6    80.9   166.8   204.1    75.9   204.7   125.5')
        print('              ...')
        print('           500.0   350.0   47.6   342.1   362.0   549.4   340.9   549.4   125.5')
        print('           500.0   400.0   47.9   386.1   412.1   554.9   385.3   554.9   125.5')
        print('           500.0   450.0   48.1   426.5   462.2   565.3   426.0   565.3   125.5')
        print('           500.0   500.0   48.3   460.1   512.2   581.9   459.9   581.9   125.5')
        print()


    # ##########
    def DumpWarnings(s):
        f = open(s.fn_warn,'w')
        for out in s.warn: f.write('%s\n' %(out))
        f.close()


        
    # ##########
    def Main(s):
        if s.VB>1: print(("INFO::%s  Main" %(s.myname)))
        # for key in s.dict.keys(): print 'dict:  %-10s  %s' %(key, s.dict[key])

        # Init
        if os.path.exists(s.fn_root): os.remove(s.fn_root)
        froot = ROOT.TFile(s.fn_root, 'recreate')
        

        for iresvar in range(len(s.resvars)):
            resvar = s.resvars[iresvar]

            # Make TGraph
            tGraph = MakeTGraph(table=s.table, coords=s.coords, resvar=resvar, coordsT=s.coordsT, resvarT=s.resvarsT[resvar], VB=s.VB)
            #tGraph.Draw("lego")
            tGraph.Write()  # for TGraph2D this happens implicit, but for TGraph it is required to get it in the file
            #print 'out: ', tGraph.GetXaxis().GetTitle()


        # Save ROOT file
        froot.Write()
        froot.Close()
        if s.VB>0: print(('INFO::TableToTGraph  Created resulting file %s' %(s.fn_root)))
       


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


        if Arg.hasget(['-coords']):  # Ex: -coords MU,M_2:M2
            zz = Arg.list()
            for z in zz:
                w = z.split(':')
                zvar = w[0]
                if len(w)>1: zvarT = w[1]
                else: zvarT = w[0]
                s.coords.append(zvar)
                s.coordsT[zvar] = zvarT

                
        if Arg.hasget(['-resvars']):  # Ex: -vars N1,N2:N20,N3,N2-N1:N2mN1,N2/N1
            zz = Arg.list()
            for z in zz:
                w = z.split(':')
                zvar = w[0]
                if len(w)>1: zvarT = w[1]
                else: zvarT = w[0]
                s.resvars.append(zvar)
                s.resvarsT[zvar] = zvarT
                #print 'dict: ', zvar, zvarT, s.resvarsT[zvar]


        if Arg.hasget(['-fnroot','-save']):
            s.fn_root = Arg.val()
            if not s.fn_root.endswith('.root'): s.fn_root += '.root'
            
            
        if Arg.hasget(['-fntable']):
            s.fn_table = Argv.val()


        if Arg.hasget('-vb'):
            s.VB = Arg.valI()
            if s.VB: print(('Verbosity level: %i' %(s.VB)))


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
                    print((s.dict))
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
    t = TableToTGraph(cmd=['ReadArg','PostInit','Main'])
############################## 

