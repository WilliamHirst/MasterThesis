import os
# ----------------------------------------------------------
def ReadDarksusy(fn):
    f = open(fn); lines = f.readlines(); f.close()
    L = lines[0]
    w = L.split()

    par = {}

    #print L
    w.pop(0); w.pop(0); w.pop(0)
    par['tests'] = w.pop(0)
    w.pop(0); w.pop(0)
    par['ds_oktot'] = int(w.pop(0) == 0)
    w.pop(0); w.pop(0)
    par['ds_unphys'] = int(w.pop(0))  # always ok, i.e. 0 ? 
    w.pop(0); w.pop(0)
    par['ds_accel'] = int(w.pop(0))   # just contains tests in binary form (9 bits), so not so useful
    w.pop(0); w.pop(0)
    
    try: par['ds_cdm'] = float(w.pop(0))
    except: par['ds_cdm'] = -2.
    try: par['ds_cdm2'] = float(w.pop(0))
    except: par['ds_cdm2'] = -2
    
    w.pop(0); w.pop(0)
    par['ds_gminus2'] = float(w.pop(0))
    w.pop(0)

    testsH = ['C','G','Q','L','Z','h','N','bsgam','rho']
    for i in range(9):
        par['ds_ok'+testsH[i]] = int(par['tests'][i] == '0')   # 0 means test is ok (not excluded)

    return par


# ----------------------------------------------------------
def ReadDarksusy2(fn):
    if not os.path.exists(fn): return {}
    f = open(fn); lines = f.readlines(); f.close()
    L = lines[0]
    w = L.split()

    par = {}

    #print L
    w.pop(0); w.pop(0); w.pop(0)
    par['tests'] = w.pop(0)  # this is the entire 9-character test string (is used towards the bottom)
    w.pop(0); w.pop(0)
    par['ds_excl_combined'] = abs(int(w.pop(0)))   # this trick ensures 0:all tests ok;  1:at least one test fails
    w.pop(0); w.pop(0)
    par['ds_unphys'] = int(w.pop(0))  # is (nearly?) always ok, i.e. 0  
    w.pop(0); w.pop(0)
    par['ds_accel'] = int(w.pop(0))   # just contains tests in binary form (9 bits), so not so useful
    w.pop(0); w.pop(0)
    
    try: par['ds_cdm'] = float(w.pop(0))
    except: par['ds_cdm'] = -2.
    try: par['ds_cdm2'] = float(w.pop(0))
    except: par['ds_cdm2'] = -2
    
    w.pop(0); w.pop(0)
    par['ds_gminus2'] = float(w.pop(0))
    w.pop(0)

    testsH = ['C','G','Q','L','Z','h','N','bsgam','rho']
    par['ds_excl_combined_except_h'] = 0
    par['ds_excl_combined_except_Ch'] = 0
    for i in range(9):
        # par['ds_excl_'+testsH[i]] = int(par['tests'][i] == '1')   # '1' means test gives exclusion  # 2012-06-11: swapped 0 and 1 for single-exclusions
        # par[testsH[i]] = int(par['tests'][i] == '0')   # 0 means test is ok (not excluded)  # same as above, just shorter
        
        if testsH[i] not in ['h']:
            if par['tests'][i] != '0': par['ds_excl_combined_except_h'] = 1
        if testsH[i] not in ['C','h']:
            if par['tests'][i] != '0': par['ds_excl_combined_except_Ch'] = 1

        par['ds_excl_'+testsH[i]] = int(par['tests'][i])   # just take what is in the file: '1' means test gives exclusion  # 2012-06-11: swapped 0 and 1 for single-exclusions
        par[testsH[i]] = int(par['tests'][i])    # 1 means excluded by the given test  # same as above, just shorter  # 2012-06-11: swapped
        
        if testsH[i] not in ['h']:
            if par['tests'][i] == '1': par['ds_excl_combined_except_h'] = 1
        if testsH[i] not in ['C','h']:
            if par['tests'][i] == '1': par['ds_excl_combined_except_Ch'] = 1    # 1 is exclusion

    return par
 

# ----------------------------------------------------------
def Darksusy_these_bits_notexcluded(par, bits):
    '''
    This procedure tests the given bits (e.g. bsgam, N, h) in darksusy output
    bits is a comma-separated list of bits, e.g. ['N','C','h']
    '''
    isok = 0
    for bit in bits:
        if bit not in list(par.keys()):
            print("Warning::Darksusy_testbits  bit '%s' not keys: %s" %(bit, str(list(par.keys()))))
            continue
        if par[bit] != 0: return 0

    # have arrived here, all bits were ok, i.e. ==0, i.e. not excluded
    return 1


# ----------------------------------------------------------
