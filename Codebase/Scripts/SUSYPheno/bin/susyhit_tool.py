import os, time, sys
from slhafrompardict import slhafrompardict
from kilelib import GetLinuxDist

# #####
def gen12nonequal(SUSYpar, calculator='susyhit'):
    # to give warning if first two generations are not identical, since susyhit (suspect) treats them as such (only reads parameters for the first generation)
    statusTxt = ''
    for gen12 in [['eR','mR'], ['dR','sR'], ['uR','cR'], ['L1','L2'], ['Q1','Q2']]:
        if SUSYpar[gen12[0]] != SUSYpar[gen12[1]]: statusTxt += '%s!=%s ' %(gen12[0], gen12[1])
    return statusTxt


# #####
def susyhit_execute_in_subdir(SUSYpar={}, par_steer={}, lha_in=[], fn_base_out='', fnmode=0, deltmp=1, silent=1):

    fnID = fn_base_out
    if fnID and not fnID.endswith('_'): fnID += '_'
    
    # this is to allow running of several susyhit instances in parallel
    if not (SUSYpar or lha_in): sys.exit('FATAL::susyhit_execute_in_subdir  Need either SUSYpar or lha_in: none given')
    if (SUSYpar and lha_in): sys.exit('FATAL::susyhit_execute_in_subdir  Need *either* SUSYpar or lha_in: both given')

    #if SUSYpar: lha_in = susyhit_make_lha_in(par=SUSYpar, mode=['ret'])
    if SUSYpar: lha_in = slhafrompardict(par1=SUSYpar, par2=par_steer, calculator=['susyhit'])
                    
    # Then make and enter subdir
    dir0 = 'tmp.susyhit_tool_%.6f' %(time.time())
    os.mkdir(dir0)
    os.chdir(dir0)

    # Write the input file
    fn = 'suspect2_lha.in'
    f = open(fn,'w')
    for out in lha_in: f.write('%s\n' %(out))
    f.close()

    dist = GetLinuxDist()

    if not 'SUSYPHENO_PATH' in list(os.environ.keys()):
        print("$SUSYPHENO_PATH not set. Exiting.")
        sys.exit()
    susyphenopath = os.environ['SUSYPHENO_PATH']

    susyhit_dir = susyphenopath+'/SUSYHIT/'
    # Make necessary links (NOTE: THESE INPUT FILES ARE JUST CONSTANTS, SO DO NOT AFFECT PARALLEL RUNNING)
    #if dist in ['RHEL5']: susyhit_dir = '/net/abel-evs/cargo/fysepf/epfshare/prog/susyhit'
    # now works on both RHEL5 and RHEL6 (may have been just confusion)
    #susyhit_dir = '/net/abel-evs/cargo/fysepf/borgeg/alt/prog/RHEL6/susyhit/susyhit_newHDECAY'  # DDD
    #elif dist in ['RHEL6']: susyhit_dir = '/net/abel-evs/cargo/fysepf/epfshare/prog/susyhit_RHEL6'  # DDD
    #else:
    #    susyhit_dir = '/net/abel-evs/cargo/fysepf/epfshare/prog/susyhit'
    #    print 'Warning  susyhit_tool  non-recognised distribution: %s  (Proceeding as if RHEL5)' %(dist)


        
    cmd = 'ln -s %s/susyhit.in .' %(susyhit_dir)
    os.system(cmd)
    os.system('pwd')  # DDD

    
    # Do the calculation  # NB
    if silent == 2: pipe = '  &>  /dev/null'
    elif silent == 1: pipe = '  >  /dev/null' 
    else: pipe = ''
    cmd = '%s/run %s' %(susyhit_dir, pipe)
    if not silent: print(cmd)
    #print cmd  # debug
    os.system(cmd)


    # Change and move filenames
    if fnmode == 0: # just move the output files out of the tmp dir. If fnID is not set, cannot run in parallel)
        for fn in ['suspect2_lha.in', 'suspect2.out', 'slhaspectrum.in', 'susyhit_slha.out']:
            cmd = 'mv %s ../%s%s' %(fn, fnID, fn)
            os.system(cmd)
    if fnmode == 1: # make nicer output (not implemented)
        pass


    # Go up and remove dir
    os.chdir('../')

    if deltmp: os.system('rm -r %s' %(dir0))

    
    
# #####
def suspect_read_outfile(fnout):
    return susyhit_read_outfile(fnout)
    

# #####
def susyhit_read_outfile(fnout):  # ... it is really the outfile of suspect (suspect2.out), see above
    if not os.path.exists(fnout):
        print('susyhit_tool::susyhit_read_out_file: Warning non-existing file: %s' %(fnout))
        return {}

    lines = os.popen('cat %s' %(fnout)).readlines()

    res = {}
    
    # Quick method to only read Warning/Error flags++ (maybe later expand to read entire out file)
    for iL in range(len(lines)):
        L = lines[iL].strip()
        #print L

        # =====
        # Note most blocks below are not filled (because we anyway can get (most of) them with pyslha's reading of slha

        if L == 'Input values:':
            cat = ''

        
        if L == 'Mass matrices and mixing angles:':
            cat = ''


        if L == 'Input non-universal soft terms at M_EWSB':
            cat = 'INPUTatEWSB'
            res[cat] = {}
            w_var = lines[iL+2].strip().split()
            w_val = lines[iL+3].strip().split()
            for i in range(len(w_var)):
                res[cat][w_var[i]] = float(w_val[i])
            iL += 3
            #<can then continue reading next lines>


        if L == 'Fermion masses and gauge couplings: Q=EWSB': 
            cat = ''


        if L == 'Mass matrices and mixing angles:': 
            cat = ''


        if L == 'Final Higgs and SUSY particle masses:':
            cat = ''



        # NOTE: reading of most of the blocks above are not implemented
        # =====

        
            
        if L == 'Low-energy/LEP precision parameter values:':
            cat = 'SU_LOWPAR'
            res[cat] = {}
            iL += 2
            w = lines[iL].strip().replace('E-','Eminus').replace('-',' -').replace('Eminus','E-').split() # 2013-12-21
            # Above: hack to fix numbers grown together (happens if negative: ex: 0.5758E-03-0.2947E-09 0.3447E-03 )
            #        Note: needed to first take aside 'E-'
            if len(w) < 3: print('Warn: SU_LOWPAR line: %s' %(w))  # 
            res[cat]['delrho'] = float(w[0])
            res[cat]['gmu-2'] = float(w[1])
            res[cat]['BR(b->sgam)'] = float(w[2])
        
        if L == 'Fine-tuning values for info: fine-tuned if >>1':
            cat = 'SU_FINETUNE'
            res[cat] = {}
            iL += 2
            w = lines[iL].strip().split()
            res[cat]['dmZ^2/mZ^2(mu^2)'] = float(w[0])
            res[cat]['dmZ^2/mZ^2(B.mu)'] = float(w[1])
            res[cat]['dmt/mt(mu)'] = float(w[2])
            res[cat]['dmt/mt(B.mu)'] = float(w[3])

        
        if L == 'Warning/Error Flags: errmess(1)-(10):':
            cat = 'FLAGS'
            res[cat] = {}
            iL += 2

            # 'flagstxt': first the txt one-liner for all flags 
            flagstxt = ''
            w = lines[iL].strip().replace('.','. ').split()
            for ww in w: flagstxt += '%3i' %(int(float(ww)))
            res[cat]['flagstxt'] = flagstxt
            # shorter (0:ok, 1:warning)
            res[cat]['flagstxt2'] = ''
            for ww in w: res[cat]['flagstxt2'] += '%i' %(abs(int(float(ww))))

            # 'flags': then each flag. Note that [0] is excluded, starting from [1], a-la fortran
            res[cat]['flags'] = ['note, starting at index 1']
            for ww in w: res[cat]['flags'].append(int(float(ww)))

            # 'flaginfo': then general info on each flat. Note [0] entry is included as general comment
            res[cat]['flaginfo'] = []
            iL += 1
            for i in range(11): 
                iL += 1
                res[cat]['flaginfo'].append(lines[iL].strip())

            # 'flagconclusions': then the remaining comments / recommendations, typically the conclusion
            res[cat]['flagconclusion'] = []
            while iL < len(lines)-1:
                iL += 1
                res[cat]['flagconclusion'].append(lines[iL].strip())
                # print 'DEBUG  Ending %i / %i  %s' %(iL, len(lines), lines[iL].strip())

        # print 'DEBUG  loop iL: %i / %i' %(iL, len(lines))
        # end of loop

        
    return res

# #####
