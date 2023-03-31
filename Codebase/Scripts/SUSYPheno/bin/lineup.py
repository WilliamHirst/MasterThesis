#!/usr/bin/env python
'''
Program : lineup.py
Version : 1.0
Author  : b.k.gjelsten@fys.uio.no
Description : 

'''

import sys,os
import bkgjelstenArgReader

# ##################################################### GLOBAL METHODS
# ##################################################### GLOBAL METHODS
# ##################################################### GLOBAL METHODS


def lineup(lines=[], optD={}, cmd=['ReadArg','PostInit','Lineup'], n=-1):
    zcmd = list(cmd)  # 2013-11-20: needed to do this, otherwise the default cmd value is changed once ReadArg is removed below. Strange.
    if ('ReadArg' in zcmd) and (not 'argv' in optD): zcmd.remove('ReadArg')
    if n>-1: optD['n'] = n  # shortcut
    oLineup = Lineup(lines=lines, cmd=zcmd, optD=optD)
    linedup = oLineup.linedup
    just = oLineup.res['just']  # 2013-11-03
    del oLineup
    returnAlignArray = optD.get('align',0)  # allows to return the alignment array as well (useful for e.g. making html table)
    if returnAlignArray:
        return linedup, just
    else: 
        return linedup


class Lineup: 

    def __init__(s, lines=[], cmd=['ReadArg','PostInit','Lineup'], optD={}):

        # ====================== PRE INIT
        if 'argv' in optD: s.argv = optD['argv']
        else: s.argv = sys.argv
        s.cmd = cmd
        s.myname = sys.argv[0].split('/').pop()
        s.VB = 0
        s.HOME = os.getenv('HOME')
        s.cwd  = os.getcwd()  # current work directory  
        #s.dir0 = '%s/XXX' %(s.HOME)
        s.dir0 = ''

        s.dict = {}

        s.warn = []
        s.fn_warn = 'warnings.txt'
        s.fn_report = 'report'
        s.report = []
        s.lines = list(lines) 


        s.lineup_nspaces = 3
        if 'n' in optD: s.lineup_nspaces = optD['n']
        if 'nspaces' in optD: s.lineup_nspaces = optD['nspaces']
        s.lineup_align = ['l','r']
        s.lineup_autoalign = 1  # skips lineup_align
        s.lineup_autoexcept = 1  # skips lines with wrong n(col)
        s.lineup_exceptlines = []
        s.lineup_exceptpattern = ['---------------','===============']
        s.lineup_adjusthline = 1  # lines in lines_except which are hline in '-','=' or '#' get their width adjusted
        s.lineup_defline = -1   
        s.lineup_header = [0]   # header lines ... will not affect the auto-alignment based on whether a column is numbers or not
        s.lineup_headerpattern = []   # header lines ... will not affect the auto-alignment based on whether a column is numbers or not
        s.lineup_wordisnumber = ['None','NaN','NONE','NAN','none','nan','-']  # these words are treated as numbers (i.e. right-aligned): e.g. 'None', 'NaN' 
        s.lineup_maxcol = 999
        s.lineup_colsonly = []
        s.lineup_colsskip = []
        s.lineup_delim = ['']


        # ====================== READ ARG
        if 'ReadArg' in s.cmd: s.ReadArg()

        # ====================== READ ARG
        if 'stdin.readlines' in s.cmd:
            zz = sys.stdin.readlines()
            for z in zz: s.lines.append(z.rstrip())
            

        # ====================== POST INIT
        if 'PostInit' in s.cmd: s.PostInit()


        # ====================== EXECUTE 
        if 'Lineup' in s.cmd:
            s.res = s.Lineup()
            s.linedup = s.res['outs']  #Lineup()

            
            # ====================== PRINT?
            if 'Print' in s.cmd:
                for out in s.linedup: print(out)


            # ====================== RETURN?
            #if 'Return' in s.cmd:
            #    return s.linedup    # apparently not allowed



        # SUMMARY
        if len(s.warn):
            print('WARNINGS (%i):' %(len(s.warn)))
            for warn in s.warn: print('  '+warn)



    # ##########
    def PostInit(s): 
        if s.dir0: s.fn_warn = '%s/%s' %(s.dir0, s.fn_warn)
        if s.dir0: s.fn_report = '%s/%s' %(s.dir0, s.fn_report)

    

    # ##################################################### CLASS METHODS
    # ##################################################### CLASS METHODS
    # ##################################################### CLASS METHODS
    # ##########
    def showHelp(s):
        print(' Usage: cat <txtfile> | %s [options]' %(s.myname))
        #print '        %s  -dict test,txt1:I,testI,3:test2,txt2a,txt2b:F,testF,4.14   # for using autodict (NB: vars need to be defined in __init__)' %(s.myname)
        print("    Ex: cat <txtfile> | %s --clean   # to clean any auto-settings" %(s.myname))
        print() 
        print("   Note: Default is auto with first line as potential header line")
        
        

    # ##########
    def DumpWarnings(s):
        f = open(s.fn_warn,'w')
        for out in s.warn: f.write('%s\n' %(out))
        f.close()

        
    # ##########
    def Lineup(s):

        lines = s.lines
        Nlines = len(lines)
    
        maxwidth = []
        lines_except = []
        autoalign = []

        outs = []


        # ### colsonly, colsskip: replace text in header with column
        for iL in s.lineup_header:  # can make more general by giving real header in separate variable
            if len(lines) <= iL: continue  # fragile? 
            whead = lines[iL].strip().split()
            for lineup_cols in [s.lineup_colsonly, s.lineup_colsskip]:
                notreplaced = []
                for iC in range(len(lineup_cols)):
                    col = lineup_cols[iC]

                    if type(col) is not int:
                        # assume it is string
                        wasreplaced = False
                        for iw in range(len(whead)):
                            headtxt = whead[iw]

                            if headtxt == col: 
                                lineup_cols[iC] = iw  # replace the text with the integer from the header
                                wasreplaced = True
                                break
                        if not wasreplaced:
                            notreplaced.append(col)

                for col in notreplaced:
                    lineup_cols.remove(col)
                    s.warn.append("Warning  Lineup: column '%s' ignored as was not found in header" %(col))
                    

    
        # ### autoexcept: find the number of columns in the table lines, except all others from aligning
        ncols = {}
        if s.lineup_autoexcept:
            for iL in range(Nlines):
                nw = len(lines[iL].strip().split())
                if nw in [0,1]: continue
                if nw not in ncols: ncols[nw] = 0
                ncols[nw] += 1
            # ---
            themax = [-1,-1]
            for nw in ncols:
                if s.VB: print(' %2i: %2i' %(nw,ncols[nw]))
                if ncols[nw] > themax[1]: themax = [nw,ncols[nw]]
            ncoltable = themax[0]
            if ncoltable == -1: autoexcept = 0
            #ncoltable = max(ncols, key=ncols.get)  # method to get key with largest value
            #print 'table: ', ncoltable
            # ---
    
    
        # ### Find maxwidth per column
        for iL in range(Nlines):
            line = lines[iL].strip()
            word = line.split()   #string.split(line)
            if s.VB>1: print('len(word):%i  ncoltable:%i' %(len(word), ncoltable))
            # --- skip some lines (will be printed as were)
            if iL in s.lineup_exceptlines or iL-Nlines in s.lineup_exceptlines or (s.lineup_autoexcept and len(word) != ncoltable):
                lines_except.append(iL)
                if s.VB>1: print("except line %2i: %s" %(iL, line))
                continue
            skip = 0
            for patt in s.lineup_exceptpattern:
                if patt in line:
                    skip = 1
                    break
            if skip:
                lines_except.append(iL)
                continue
            # --- end skipping
    
    
            # allow changes to skip columns:
            if s.lineup_colsonly:
                line2 = ''
                
                #for iw in range(len(word)):
                #    if iw in s.lineup_colsonly: line2 += '  '+word[iw]

                for iw in s.lineup_colsonly:  # rearranging too!
                    if iw < len(word): line2 += '  '+word[iw]
                    
                lines[iL] = line = line2.strip()
                word = line.split()

                
            if s.lineup_colsskip:
                line2 = ''
                for iw in range(len(word)):
                    if iw not in s.lineup_colsskip: line2 += '  '+word[iw]
                lines[iL] = line = line2.strip()
                word = line.split()
            # ---
    
            # Go through the columns
            for iw in range(len(word)):
                # fills out maxwidth upon need
                if(iw==len(maxwidth)):
                   maxwidth.append(0)  # init first time
                   autoalign.append('r')
    
                if(len(word[iw])>maxwidth[iw]): maxwidth[iw]=len(word[iw])
    
                # --- autoalign
                for patt in s.lineup_headerpattern:
                    if patt in line and iL not in s.lineup_header:
                        s.lineup_header.append(iL)
                        # print s.lineup_header
                        break
                if iL not in s.lineup_header:  # and iL not in lineup_exceptlines:
                    try:
                        if word[iw].endswith('%'): float(word[iw][:-1])  # allows percentages
                        else: float(word[iw])
                        # if is float, don't have to do anything, because was already init'ed as float, i.e. 'r'
                    except:
                        if word[iw] not in s.lineup_wordisnumber:
                            # print iL, iw, word[iw]
                            autoalign[iw] = 'l'  # if one time not a number, then left-align
                # ---
    
    
        # ### Combine autoalign and lineup_align (simple for now:one or the other, no combination)
        if s.lineup_autoalign: just = autoalign
        else: just = s.lineup_align
    
        # print just
        njust = len(just)
        maxcol = s.lineup_maxcol
        #rightmargin = sum(maxwidth) + s.lineup_nspaces * (len(maxwidth)-1)   # need to generalise when delim is set
        # if inserting text the total width must be calculated like this:

        # Make s.lineup_delim complete in case it is not
        #if s.lineup_delim: 
        for i in range(len(s.lineup_delim), len(maxwidth)):
            s.lineup_delim.append(s.lineup_nspaces*" ")
        s.lineup_delim.append("")
                
        #if len(s.lineup_delim) > 0:
        rightmargin = len(s.lineup_delim[len(maxwidth)])
        for i in range(len(maxwidth)):
            rightmargin += len(s.lineup_delim[i]) 
            rightmargin += maxwidth[i]
        #if len(s.lineup_delim) > len(maxwidth): rightmargin += len(s.lineup_delim[len(maxwidth)-1])

        
        # print rightmargin
        for iL in range(Nlines):
            line = lines[iL].strip()
            word = line.split()
            outline=""
    
            if iL in lines_except:
                if line != '' and s.lineup_adjusthline:
                    if line == len(line)*'-' or line == len(line)*'=' or line == len(line)*'#':
                        line = rightmargin*line[0]
                #print line
                outs.append(line)
                continue


            for iw in range(len(word)):

                if iw > maxcol-1:  # not sure what the intention of this was. (never/rarely used)
                    outline+=" "+word[iw]
                    continue
    
                #if iw == 0 and s.lineup_delim: outline += s.lineup_delim[iw]
                #if 0 < iw < len(s.lineup_delim): outline += s.lineup_delim[iw] + s.lineup_nspaces*" "
                #if iw>0: outline += s.lineup_nspaces*" "
                outline += s.lineup_delim[iw]
                
                if iw > njust-1:
                    thisjust = just[njust-1]
                else:
                    thisjust = just[iw]
    
                if thisjust in ('l','L'):
                    outline+=word[iw].ljust(maxwidth[iw])
                elif thisjust in ('r','R'):
                    outline+=word[iw].rjust(maxwidth[iw])
                elif thisjust in ('c','C'):
                    # outline+=word[iw].center(maxwidth[iw]+myspace)
                    outline+=word[iw].center(maxwidth[iw])
                else:
                    print("non-allowed code: %s" %thisjust)
    
            outline = outline.rstrip()
            #print outline
            outs.append(outline)


        res = {}
        res['outs'] = outs
        res['just'] = just
        return res
    
 

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


        if Arg.has(['--clean','--c']):
            # this option is to clean any default auto-settings 
            s.lineup_header = []
            s.lineup_wordisnumber = []

        if Arg.hasget(['-delim']):
            s.lineup_delim = Arg.list()
            #print s.lineup_delim
            s.lineup_nspaces = 0  # <-- so removes those ... need to set manually with -nspaces


        if Arg.hasget(['-nspaces','-n']):
            s.lineup_nspaces = Arg.valI()
        
        if Arg.hasget('-header'):
            for z in Arg.list():
                try: s.lineup_header.append(int(z))
                except: s.lineup_headerpattern.append(z)
        
        
        if Arg.hasget('-colsonly'):
            s.lineup_colsonly = Arg.listIif()
        
        if Arg.hasget('-colsskip'):
            s.lineup_colsskip = Arg.listIif()
        
        if Arg.hasget(['-likenumberclean','-isnumberclean','-wordisnumberclean','-treatlikenumberclean']):
            s.lineup_wordisnumber = Arg.list()

        if Arg.hasget(['-likenumber','-isnumber','-wordisnumber','-treatlikenumber']):  # Note, this adds to the existing wordisnumbers
            s.lineup_wordisnumber += Arg.list()
        
        if Arg.hasget('-exceptlines'):
            s.lineup_exceptlines = Arg.listI()
        
        if Arg.hasget('-align'):
            s.lineup_align = Arg.list()
        
        if Arg.hasget('-autoalign'):
            s.lineup_autoalign = Arg.valI()
        
        if Arg.hasget('-autoexcept'):
            s.lineup_autoexcept = Arg.valI()
        
 

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
    t = Lineup(cmd=['ReadArg','stdin.readlines', 'PostInit','Lineup','Print'])
############################## 

