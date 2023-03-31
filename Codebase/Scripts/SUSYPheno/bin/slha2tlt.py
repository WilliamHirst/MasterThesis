#!/usr/bin/env python

"""
Program : slha2tlt.py
Version : 1.0
Author  : b.k.gjelsten@fys.uio.no
Description : Make various tables from slha files  (see slha2tlt.py -h)


TO DO:

o allow printing/saving of a "flat dict table" (a tct, two-column-table)
o maybe also printing/saving a "structured dict table"
o inspect problems, update status (also in -h)
o PROBLEM with reading the pickles, ... there is key '#' somewhere ...  [2014-08: still a problem?]
o to use var2block for format issues, need to '--remake' if exists [should fix this]
o Option to add a common list of vars to be added to each table 

"""

import sys,os 
from bkgjelstenArgReader import ArgReader

from kilelib import WriteToFile, LoadPickle, WritePickle, SaveHistory, dict2tlt
from lineup import lineup
from SlhaTools import slha2flatdict, snames, ShowFlatSlhaVars, SlhaDecay

# ##################################################### GLOBAL METHODS
# ##################################################### GLOBAL METHODS
# ##################################################### GLOBAL METHODS


class slha2tlt:  

    def __init__(s, cmd=[], optD={}):

        # ====================== PRE INIT
        if 'argv' in optD: s.argv = optD['argv']
        else: s.argv = sys.argv
        s.cmd = cmd
        s.myname = sys.argv[0].split('/').pop()
        s.VB = 1

        #maiken change
        #s.HOME = os.getenv('HOME')
        s.HOME = os.getcwd()
        s.cwd  = os.getcwd()  # current work directory  
        #s.dir0 = '%s/XXX' %(s.HOME)
        s.dir0 = '.'
        s.dir1 = '%s' %(s.dir0)

        s.dict = {}

        s.warn = []
        s.fn_warn = 'slha2tlt_warnings.txt'
        s.fn_report = 'report'
        s.report = []

        s.dict['fn_pickleend']   = '_slha2flat.pickle'
        s.dict['fn_flatoutsend'] = '_slha2flat.txt'
        s.dict['fn_slhaend'] = '.slha'
        s.dict['fn_tltend'] = '.tlt'
        #s.dict['val_missing'] = -9
        s.dict['tlt_nspaces'] = 2  # lineup

        s.dict['forceremakeflat'] = 0   # make pickle
        s.dict['saveflat'] = 1
        s.dict['saveflat_outs'] = 1
        s.dict['forcename'] = 1

        s.dict['tlt_varprepend'] = ''  # can prepend or append some string 
        s.dict['tlt_varappend'] = ''   # useful e.g. if comparing two tables intended to be identical

        s.dict['tlt_varnotfound_include'] = 1
        s.dict['tlt_varnotfound_val'] = '0'    # -99997'  0 will stand out as '0' compared to e.g. '0.000', so it is safe enough
        s.dict['tlt_valisnone_val'] = -99995   # 
        s.dict['tlt_format0'] = "%.3f"
        s.dict['tlt_delim_var'] = ','
        s.dict['tlt_delim_rename'] = '::'
        s.dict['tlt_delim_name'] = ':::'
        
        s.dict['BR_groups_use'] = ['v','l','q']
        
        s.optBR = {}  # used to configure decay, e.g. 'groups_use' (l,q,v)

        s.ffn = ''   # file with filenames
        s.fns = []
        s.tables = []
        #s.tables_all = ['MASS','EXTPAR','MIX']
        
        s.vars = {}
        s.vars['masses'] = snames['all']  # SlhaTools
        s.vars['mweakinos'] = snames['weakinos']
        s.vars['msleptons'] = snames['sleptons']
        s.vars['mhiggses']  = snames['higgses'] 
        s.vars['mcolour']   = snames['colour']

        s.vars['extpar_part1'] = ['Minput','Q','M1','M2','M3', 'At','Ab','Au','Ad','Ae', 'mu', 'TB(MX)', 'mA']
        s.vars['extpar_sleptons'] = ['eR','mR','TR','L1','L2','L3']
        s.vars['extpar_squarks'] = ['bR','tR','Q3', 'dR','sR','uR','cR','Q1','Q2']
        s.vars['extpar_sfermions'] = s.vars['extpar_sleptons'] + s.vars['extpar_squarks']
        s.vars['extpar'] = s.vars['extpar_part1'] + s.vars['extpar_sfermions']

        s.vars['extpar_colour'] = ['M3', 'bR','tR','Q3', 'dR','sR','uR','cR','Q1','Q2', 'Ab','At']        
        s.vars['colour'] = s.vars['extpar_colour'] + s.vars['mcolour']

        s.vars['msoft_part1'] = ['M1(Q)','M2(Q)','M3(Q)', 'At(Q)','Ab(Q)','Au(Q)','Ad(Q)','Ae(Q)', 'mu(Q)', 'TB(Q)', 'mA']
        s.vars['msoft_sfermions'] = ['eR(Q)','mR(Q)','TR(Q)','L1(Q)','L2(Q)','L3(Q)','bR(Q)','tR(Q)','Q3(Q)', 'dR(Q)','sR(Q)','uR(Q)','cR(Q)','Q1(Q)','Q2(Q)']
        s.vars['msoft'] = s.vars['msoft_part1'] + s.vars['msoft_sfermions']

        s.vars['Nmix1'] = ['N1b','N1w','N1h', 'N2b','N2w','N2h', 'N3b','N3w','N3h', 'N4b','N4w','N4h']
        s.vars['Nmix2'] = ['N1g','N1h1','N1h2', 'N2g','N2h1','N2h2', 'N3g','N3h1','N3h2', 'N4g','N4h1','N4h2']
        s.vars['Nmix3'] = ['N1b','N1w','N1h','N1h1','N1h2', 'N2b','N2w','N2h','N2h1','N2h2', 'N3b','N3w','N3h','N3h1','N3h2', 'N4b','N4w','N4h','N4h1','N4h2']
        s.vars['C1mix'] = ['C1w_U','C1h_U','C1w_V','C1h_V']
        s.vars['C2mix'] = ['C2w_U','C2h_U','C2w_V','C2h_V']
        s.vars['Cmix'] = s.vars['C1mix'] + s.vars['C2mix']
        s.vars['3genmix'] = ['t1R','b1R','T1R']

        s.vars['mixing1'] = s.vars['3genmix'] + s.vars['C1mix'] + s.vars['Nmix1']
        s.vars['mixing2'] = s.vars['3genmix'] + s.vars['C1mix'] + s.vars['Nmix2']
        s.vars['mixing3'] = s.vars['3genmix'] + s.vars['C1mix'] + s.vars['Nmix3']

        s.vars['scales'] = ['FULLCALC','ICALC','Minput','Q','GMEANstop','~t1','~t2','GMEANsbot','~q4mean','q4mean','TB(Z)','TB(Q)','TB(MX)','mA','mu(Q)','mu(MX)','VEV(Q)']
        s.vars['elweak'] = ['M1','M2','mu','TB(MX)','TB(Q)','mA', 'h','At'] + s.vars['mweakinos'] 
        s.vars['sleptonR1'] = ['eR','~eR']
        s.vars['sleptonR']  = ['eR','mR', '~eR','~mR']
        s.vars['sleptonL1'] = ['L1','~eL','~Ve']
        s.vars['sleptonL']  = ['L1','L2', '~eL','~mL', '~Ve', '~Vm']
        s.vars['stau1']     = ['TR','L3','AT', '~T1'] 
        s.vars['stau']      = ['TR','L3','AT', '~T1','~T2'] 

        # Hm ... alternative would be to let q and l have variable meaning .. (would be less clutter, but less safe for the user)
        s.vars['BR_N2_lq'] = ['N2_N1_h','N2_N1_Z','N2_N1_y','N2_N1_q_q','N2_N1_b_b','N2_N1_l_l','N2_N1_T_T','N2_N1_v_v','N2_C1_W','N2_C1_q_q','N2_C1_l_v','N2_C1_T_v']
        s.vars['BR_N2_lQ'] = ['N2_N1_h','N2_N1_Z','N2_N1_y','N2_N1_Q_Q','N2_N1_l_l','N2_N1_T_T','N2_N1_v_v','N2_C1_W','N2_C1_q_q','N2_C1_l_v','N2_C1_T_v']
        s.vars['BR_N2_Lq'] = ['N2_N1_h','N2_N1_Z','N2_N1_y','N2_N1_q_q','N2_N1_b_b','N2_N1_L_L','N2_N1_v_v','N2_C1_W','N2_C1_q_q','N2_C1_L_v']
        s.vars['BR_N2_LQ'] = ['N2_N1_h','N2_N1_Z','N2_N1_y','N2_N1_Q_Q','N2_N1_L_L','N2_N1_v_v','N2_C1_W','N2_C1_Q_Q','N2_C1_L_v']

        s.vars['BR_C1_lq'] = ['C1_N1_W','C1_N1_q_q','C1_N1_l_v','C1_N1_T_v']
        s.vars['BR_C1_lQ'] = ['C1_N1_W','C1_N1_Q_Q','C1_N1_l_v','C1_N1_T_v']
        s.vars['BR_C1_Lq'] = ['C1_N1_W','C1_N1_q_q','C1_N1_L_v']
        s.vars['BR_C1_LQ'] = ['C1_N1_W','C1_N1_Q_Q','C1_N1_L_v']


        s.formats = {}
        for z in ['MASS','EXTPAR','mu','mu(Q)','~q4mean','q4mean','GMEANstop','GMEANsbot','Q']: s.formats[z] = '%.1f'
        for z in ['CALC','CALCver','FULLCALC']: s.formats[z] = '%s'
        for z in ['ICALC']: s.formats[z] = '%i'
        s.formats['TB(MX)'] = '%.3f'  # hack

        s.vars['rest'] = ['TB','TB(MX)','TB(Q)','t','b','VEV(Q)', 'Yt','Yb','YT']

        s.dict['commonvars'] = []

        s.prevar = ''
        s.preval = ''

        s.fn_finaltablebase = 'table'   # 2014-08-26: '' no longer default
        s.dict['lineup1'] = 'lineup_removeheaders.sh'
        s.dict['lineup2'] = 'lineup_removeheaders_addfilename.sh'
        s.dict['lineup']  = s.dict['lineup1']
        
        s.fn_history = 'history_slha2tlt.txt'
        s.fn_history_global = '%s/.history_slha2tlt.txt' %(s.HOME)
        s.orig_argv = list(sys.argv)
        s.dict['savehistory'] = 1

        
        # ====================== READ ARG
        if 'ReadArg' in s.cmd: s.ReadArg()


        # ====================== POST INIT
        if 'PostInit' in s.cmd: s.PostInit()


        # ====================== EXECUTE 
        if 'Main' in s.cmd: s.Main()



    # ##########
    def PostInit(s): 
        s.fn_warn = '%s/%s' %(s.dir1, 'warnings_slha2tlt.txt')
        s.fn_report = '%s/%s' %(s.dir1, 'report')

        if s.dict['savehistory']: SaveHistory(fn=s.fn_history, argv=s.orig_argv, opt=['stripcmd','date'])
        if s.dict['savehistory']: SaveHistory(fn=s.fn_history_global, argv=s.orig_argv, opt=['stripcmd','date','dirsameline'])

        s.optBR['groups_use'] = s.dict['BR_groups_use']

    

    # ##################################################### CLASS METHODS
    # ##################################################### CLASS METHODS
    # ##################################################### CLASS METHODS
    # ##########
    def showPredefined(s):
        print('Standard table vars:')
        for cat in sorted(s.vars.keys()):
            print('   %-18s :  %s' %(cat, s.vars[cat]))
        
    
    # ##########
    def showHelp(s):
        print('DESCRIPTION')
        print('   slha2tlt.py takes a list of slha-files as input and produces for each: ')
        print('     - a flat dictionary with the variables which is pickled to disk')
        print("       (Note: if the pickle already exists, it will be read instead of the slha file. To remake pickle, add '--remake')")
        print('     - two-line-table (tlt) file(s) per scenario with e.g. masses, elweak parameters, ...')
        print('       (there are predefined suggestions, or you can define your own)')
        print('       (tlt: first line contains the variable names (cannot contain space), second line contains the value)')
        print('     - one summary table per table definition where all the scenarios are combined')
        print('   The intention is to make for easier access of the info in slha files.')
        print('    (The tlt-output is a practical base for making tables with one line per scenario (and the quantity on the first line.)')
        print('    Such tables are in turn appropriate format for plotting with the ROOT-plotter script munch.py)')
        print() 

        print()
        print('USAGE')
        print('  ls *.slha | slha2tlt.py -tables <tabledescriptions>  [options]') 
        print('     <tabledescriptions> is a list of tabledefs separated by ,,')
        print('     The tabledefs can consist of ')
        print('          predefined var lists, see :   slha2tlt.py --showpredefined')
        print('          self-selected block vars  :   slha2tlt.py --showvars   (possibly in combination with predefined var lists)')
        print("    Filename:")
        print("     The filename of a given tlt takes the slha-file as base, replaces '.slha' with '.tlt' and inserts a table description.")
        print("     If the table is predefined the table description takes that name by default.")
        print("     Otherwise the table description needs to be added after the table definition, e.g. '-tables ~T1,h,mA:::mytabledesc'")
        print("     (The ':::mytabledesc' can also be used for predefined tables to give them other names.)")
        print() 
        print()
        print("OPTIONS  (see --h for full list (though without explanations)")
        print("   -dict varprepend,<PREPEND>:varappend,<APPEND>  :  prepend or append (common) string to vars")
        print("     (any of  '-dict varprepend,A'  or  '-dict varappend,B'  or  '-dict varprepend,A:varappend,B'")
        print("   --nosummarytable  : drop the summary tables")
        print("   -ffn <filename>   : can use a file with the list of slha-files instead of piping the list")
        print("   --withfilename    : this adds the tlt-filenames to the summary tables (sometimes useful if the filename contains info not in the table)")
        print("   --nohist          : do not add to history  (default is to update %s and %s)" %(s.fn_history, s.fn_history_global))
        print("   -commonvars <varlist>  : vars to insert in each table selected (e.g. scan variables)")
        print("   -vb <0/1/..>      : change verbosity (D=1)")
        print("   -h                : this help message")
        print("   --h               : brute info on implemented options")
        print("   --dict            : dump various default values (which can be changed)")
        print()
        print('EXAMPLES:')
        print('  ls *.slha | slha2tlt.py  -tables mixing3     # (one predefined table)')
        print("  ls *.slha | slha2tlt.py  -tables M1,M2,mu,'TB(MX)','TB(Q)',mixing3:::mixing     # the predefined table with some relevant additional parameters +  name (:::) to 'mixing'")
        print('  ls *.slha | slha2tlt.py  -tables masses,,elweak,,extpar,,rest,,mixing3     # (using predefined only)')
        print("  ls *.slha | slha2tlt.py  -tables masses,rest:::mytable1,,Nmix1,,elweak     # (predefined, rename the 'rest' table to 'mytable1')")
        print("  ls *.slha | slha2tlt.py  -tables ~T1,TR,AT,h,mA:::tablex,,masses           # (self-selected vars in table 'tablex')   # 2014-01-29")
        print("  ls *.slha | slha2tlt.py  -tables elweak,~T1,TR,AT,h,mA:::tablex,,masses    # (predefined and self-selected vars in table 'tablex')   # 2014-01-29")
        print("  ls *.slha | slha2tlt.py  -tables masses,,extpar,,3genmix,C1mix,Nmix1:::mixing,,elweak,,rest ")
        #print "  ls *.slha | slha2tlt.py  -tables masses,,extpar,,3genmix,C1mix,Nmix1:mixing,,elweak,,rest  --remakeflat   # 2014-02-06  (there is currently a problem with reading back existing flat pickle, fix is to just forceremakeflat. Hm, this is now fixed, so not needed)"
        print("  ls *.slha | slha2tlt.py  -tables FULLCALC,ICALC,Minput,Q,'TB(MX)','TB(Q)',mA,mu,'mu(Q)','VEV(Q)':::scales")   # Some variables need to be protected with ''"
        print("  ls *.slha | slha2tlt.py  -tables BR_N2_lq  -BR_groups_use v,l,q  -commonvars M1,M2,mu  # Predefined BRs of N2 with l=e,m and q=u,d,s,c)")
        print("  ls *.slha | slha2tlt.py  -tables BR_N2_LQ  -BR_groups_use v,L,Q  -commonvars M1,M2,mu  # Predefined BRs of N2 with L=e,m,T and q=u,d,s,c,b)  (can also use Lq and lQ)")
        print("  ls *.slha | slha2tlt.py  -tables BR_C1_lq  -BR_groups_use v,l,q                        # Predefined BRs of C1")
        print("  ls *.slha | slha2tlt.py  -tables BR_C1_lq,N2_N1_l_l,N2_N1_Z:::test  -BR_groups_use v,l,q      # Predefined BRs of C1 some of N2")
        print("  ls *.slha | slha2tlt.py  -tables ~N2_~N1_*,~N2_~C1_*,N2_C1_W:::test  -BR_groups_use v,l,q     # Complex BRs (with '*') on the fly")

        print("     Note: if vars contain e.g. parenthesis, ex. TB(Q), you will need to protect the arguments with '' like shown above.")
        
    
    # ##########
    def Main(s):
        if s.VB>1: print("INFO::%s  Main" %(s.myname))
        # for key in s.dict.keys(): print 'dict:  %-10s  %s' %(key, s.dict[key])

        # Get filenames
        if s.ffn: f = open(s.ffn); s.fns = f.readlines() ; f.close()
        else:
            if s.VB > 0: print('INFO::slha2tlt  getting filelist from stdin...')
            s.fns = sys.stdin.readlines()


        fn_tlts = {}  
        for ifn in range(len(s.fns)):

            # 0) Init
            fn = s.fns[ifn].strip()
            if len(fn.split()) != 1: print('Skipping nonstandard filename:  %s' %(fn)); continue
            if fn.endswith(s.dict['fn_pickleend']): fn_base = fn.replace(s.dict['fn_pickleend'],'')
            elif fn.endswith(s.dict['fn_slhaend']): fn_base = fn.replace(s.dict['fn_slhaend'],'')
            else:
                fn_base = fn
                print('Skipping nonstandard filename:  %s' %(fn))
            fn_slha = fn_base + s.dict['fn_slhaend']
            fn_flat = fn_base + s.dict['fn_pickleend']
            fn_flatout = fn_base + s.dict['fn_flatoutsend']
            
            # 1) Get the flat dict
            if os.path.exists(fn_flat) and not s.dict['forceremakeflat']:
                if s.VB>1: print("INFO::slha2tlt  Loading existing pickle: %s" %(fn_flat))
                resflat = LoadPickle(fn_flat)
                res = {'flat':resflat}  # hack to
            else:
                if s.VB>1 and os.path.exists(fn_flat): print("INFO::slha2tlt  Recreating flat dict (& pickle): %s" %(fn_flat))
                elif s.VB>1: print("INFO::slha2tlt  Creating flat dict (& pickle): %s" %(fn_flat))
                if os.path.exists(fn_flat): os.remove(fn_flat)  # remove if not use 
                res = slha2flatdict(fn=fn, VB=s.VB, opt=s.optBR)
                if res == {}:
                    print('Error::slha2tlt  Skipping problem slha file: %s' %(fn))
                    continue
                
                resflat = res.get('flat',{})
                if s.dict['saveflat']: WritePickle(fn=fn_flat, thepickle=resflat, VB=s.VB-1)  # 2014-02-05 changed from res to res['flat']
                if s.dict['saveflat_outs']: WriteToFile(fn=fn_flatout, outs=res.get('flat_outs'), VB=s.VB-1)
            #F = resflat

            if s.VB>3: print("resflat = %s" %(resflat))

            var2block = res.get('var2block',{})   # if dict is reused, this is empty. To have it need to '--remake'

            # 2) Make std tables   (should try to use kilelib::dict2tlt)
            for table in s.tables:

                # A) First construct the variable list
                # table can composite, can structure like this: table=name1,name2:fnadd
                # this will take vars from cat1 and cat2 and the file will get filename from fnadd
                # if fnadd not given, the first part (full part) will also give filename (a bit odd)
                #w = table.split(':')
                w = table.split(s.dict['tlt_delim_name'])
                tltdef = w.pop(0)  # this is the variable list


                
                if w: fnadd = w.pop()  # is this tablename construction in use??? Yes
                else: fnadd = tltdef   # allows to take name of predefined tables (Hmm)

                if s.dict['tlt_delim_var'] in fnadd:   # this is to avoid that free tables are given very bad names
                    s.warn.append("Skipping table %s since no name is given (no '%s<name>' at the end). To switch off this requirement add '-dict I,forcename,0'" %(table, s.dict['tlt_delim_name']))
                    print(s.warn[-1])
                    continue

                for zvar in reversed(s.dict['commonvars']): tltdef = '%s,%s' %(zvar, tltdef)  # neatly placing commonvars (if any) in front

                # Check if there are complex BRs (format 'BR_...._*')
                complexBRs = []
                for zvar in tltdef.split(','): 
                    #if zvar.startswith('BR_') and 
                    if zvar.endswith('_*'): complexBRs.append(zvar)
                if complexBRs: 
                    #print 'complexBRs: ', complexBRs
                    scen = SlhaDecay(fn=fn_slha, optD=s.optBR)
                    flatcomplexBRs, keys = scen.string2dict(complexBRs)
                    #print 'flatcomplexBRs: ', flatcomplexBRs
                    resflat.update(flatcomplexBRs)
                    del scen
                # ---

                outstlt = dict2tlt(tltdef=tltdef, thedict=resflat, prevar=s.prevar, preval=s.preval, D=s.dict, formats=s.formats, var2block=var2block, format0=s.dict['tlt_format0'], delim_var=s.dict['tlt_delim_var'], delim_rename=s.dict['tlt_delim_rename'], varcollection=s.vars)

                fn_table = '%s_slha_%s%s' %(fn_base, fnadd, s.dict['fn_tltend'])
                WriteToFile(fn=fn_table, outs=outstlt, VB=s.VB-1)
                # put table names in dictlist
                if fnadd not in fn_tlts: fn_tlts[fnadd] = []
                fn_tlts[fnadd].append(fn_table)


                """   # this is now put in kilelib::dict2tlt (but is a bit on the specialised side?)
                cats = tltdef.split(',')
                vars = []
                #print cats
                
                # Check if name given
                if s.dict['forcename'] and len(cats) > 1 and not isnamed:
                    s.warn.append("Skipping table %s since no name is given (no ':<name>' at the end). To switch off this requirement add '-dict I,forcename,0'" %(table))
                    print s.warn[-1]
                    continue
                
                for cat in cats: 
                    if cat in s.vars: vars += s.vars[cat]  # could ensure that there is no duplicates
                    else: vars += [cat]  # assume it is a free variable
                    #else: s.warn.append('%s   (sub)table %s nonexistent. Skipping.' %(fn_base, cat))

                if s.VB>1: print 'table: %s   vars: %s' %(table, vars)
                
                head = s.prevar  # ''
                out  = s.preval  # ''
                for var in vars:
                    if var not in resflat:
                        #s.warn.append('%s  Table %s: var %s not in flatdict. Value set to %s' %(fn_base, table, var, s.dict['val_missing']))
                        if s.dict['tlt_varnotfound_includ']: 
                            head += '  ' + s.dict['tlt_varprepend'] + var + s.dict['tlt_varappend']
                            out  += '  ' + s.dict['tlt_varnotfound_val']
                            s.warn.append("%s  Table %s: var %s not in flatdict. Setting value to %s.  (To skip, use ' -dict I,varnotfound_include,0') " %(fn_base, table, var, s.dict['tlt_varnotfound_val']))
                        else:
                            s.warn.append("%s  Table %s: var %s not in flatdict. Skipping value (fragile)  (To keep, use ' -dict I,varnotfound_include,1')" %(fn_base, table, var))
                        continue
                    if var in s.formats: 
                        format = s.formats[var]  # var-specific (rare)
                    elif var2block.get(var,'asdfasdf') in s.formats:
                        format = s.formats[var2block[var]] 
                    else: format = '%.3f'
                    # Hmm ... the usage of var2block ruins the reusage of pickles ... [STILL THE CASE?]

                    val = resflat[var]
                    if val == None: val = s.dict['tlt_valisnone_val']
                    head += '  ' + s.dict['tlt_varprepend'] + var + s.dict['tlt_varappend']
                    out  += '  ' + format %(val)


                # Lineup and write to file 
                outs = lineup([head,out], n=s.dict['tlt_nspaces'])
                WriteToFile(fn=fn_table, outs=outs, VB=s.VB-1)
                """

            # end of treating this file/scenario

        # 
        WriteToFile(fn=s.fn_warn, outs=s.warn, VB=s.VB-1)
        if len(s.warn) > 0 and s.VB>0: print('Note: there were warnings (%i lines in %s)' %(len(s.warn), s.fn_warn))
            

        # Allow to also make a combined table on the fly
        if s.fn_finaltablebase:
            for fnadd in fn_tlts:
                fn_finaltable = '%s_%s.txt' %(s.fn_finaltablebase, fnadd)
                cmd = "("
                for fn in fn_tlts[fnadd]: cmd += ' ls %s ;' %(fn.strip())
                cmd += ") | %s  -nspaces %s  >  %s" %(s.dict['lineup'], s.dict['tlt_nspaces'], fn_finaltable)
                if s.VB>1: print(cmd)
                os.system(cmd)



    # ##########
    def ReadArg(s): 

        # ################################### ARGUMENT READING
        Arg = ArgReader(s.argv, VB=0)

        if Arg.has(['--dict','--showdict']):
            print('DUMPING DEFAULT VALUES IN THE VARIABLE DICTIONARY, s.dict')
            print('   These can be changed by e.g. -dict I,var1,val1:var2,val2:F,var3,val3:var4,val4')
            print('   where I/F denote integer and float. In the example var2 inherits I from var1,')
            print('   while var4 is neither given as I or F and therefore is taken as string.')
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
            if s.VB: print('Verbosity level: %i' %(s.VB))


        if Arg.hasget('-tables'):
            s.tables = Arg.list(',,')  # Note

        if Arg.hasget(['-ffn','-ff']):
            s.ffn = Arg.val()

        if Arg.hasget(['-nspaces']):
            s.dict['nspaces'] = Arg.valI()


        if Arg.hasget(['-prevar']):  # -prevar 'M1 M2 M3' 
            s.prevar = Arg.val()
            
        if Arg.hasget(['-preval']):  # -preval '100 150 150'   # preval and prevar should match, a space will be inserted before the regular 2-liner
            s.preval = Arg.val()
            
        if Arg.hasget(['-commonvars']): 
            s.dict['commonvars'] = Arg.list()

        if Arg.hasget(['-BR_groups_use']): 
            s.dict['BR_groups_use'] = Arg.list()

        if Arg.has('--showpredefined'):
            s.showPredefined()
            sys.exit()  # not superelegant ... 

        if Arg.has('--showvars'):
            ShowFlatSlhaVars(show=1)
            sys.exit()  # not superelegant ... 


        if Arg.has(['--remakepickle','--remakedict','--remakeflat','--remake']):
            s.dict['forceremakeflat'] = 1

        if Arg.hasget(['-fn_finaltable','-fn_finaltablebase','-finaltable']):
            s.fn_finaltablebase = Arg.val()
            for z in ['.txt']:
                if s.fn_finaltablebase.endswith(z): s.fn_finaltablebase = s.fn_finaltablebase[:-len(z)]

        if Arg.has(['--nosummary','--nosummarytable']):
            s.fn_finaltablebase = ''

        if Arg.has(['--withfilename','--withfilenames','--withfn','--withfns']):
            s.dict['lineup'] = s.dict['lineup2']
            
        if Arg.has(['--nohist']):
            s.dict['savehistory'] = 0

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
    t = slha2tlt(cmd=['ReadArg','PostInit','Main'])
############################## 

