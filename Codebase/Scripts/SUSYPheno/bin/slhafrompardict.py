# ##########
def susyhitblockfrompardict(par):
    outs = []
    outs.append("Block SU_ALGO  # !Optional SUSPECT v>=2.3* block: algorithm control parameters")
    outs.append("# !IF block absent (or if any parameter undefined), defaut values are taken")
    outs.append("     2    21  # 2-loop RGE (defaut, 1-loop RGE is: 11 instead)")
    outs.append("     3    1   # 1: g_1(gut) = g_2(gut) consistently calculated from input")
    outs.append("#   (other possibility is 0: High scale input =HIGH in block EXTPAR below)")
    outs.append("     4    2   # RGE accuracy: 1: moderate, 2: accurate (but slower)")
    # outs.append("     6    1   #  1: M_Hu, M_Hd input (default in constrained models)")
    outs.append("     6    0   #  1: M_Hu, M_Hd input (default in constrained models)")
    outs.append("#        (other possibility 0: MA_pole, MU(EWSB) input instead)")
    outs.append("     7    2   #  choice for sparticles masses rad. corr. (=/= h):")
    outs.append("#               2 ->all (recommended, defaut); 1->no R.C. in squarks & gauginos.")
    outs.append("     8    %i   # 1 (defaut): EWSB scale=(mt_L*mt_R)^(1/2)" %(par.get('EWSB_scale_mode',1)))
    outs.append("#         (Or = 0: arbitrary EWSB scale: give EWSB in Block EXTPAR below)")
    outs.append("     9    2   # Final spectrum accuracy: 1 -> 1% acc.; 2 -> 0.01 % acc.(defaut)")
    outs.append("     10   %i   # Higgs boson masses rad. corr. calculation options:" %(par.get('Hmass_mode',2)))
    outs.append("#             A simple (but very good) approximation (advantage=fast)  : 0")
    outs.append("#             Full one-loop calculation                                : 1")
    outs.append("#             One-loop  + dominant DSVZ 2-loop (defaut,recommended)    : 2")
    outs.append("     11   0   # Higher order Higgs 'scheme' choice in rad. corr. at mZ:")
    outs.append("#          RUNNING DRbar Higgs masses at loop-level at mZ (defaut)    : 0")
    outs.append("#          POLE          Higgs masses at loop-level at mZ             : 1")
    outs.append("#")

    return outs

# ##########
def softsusyblockfrompardict(par, keepempty=1):
    outs = []
    outs.append    ("Block SOFTSUSY       # SOFTSUSY specific inputs")
    if 'TOLERANCE' in par:
        outs.append("  1   %8.3e    # tolerance [D=1.0e-4] : desired fractionaly accuracy in output" %(float(par.get('TOLERANCE',1e-4))))
    if 'MIXING' in par:
        outs.append("  2   %3.1f         # mixing [D=1.0]: 0.0:zero quark mixing, but includes diagonal terms of first two families ; 1.0/2.0: all mixing at Mz to be in up or down quark sector" %(float(par.get('MIXING',1.0))))
    if 'PRINTOUT' in par:
        outs.append("  3   %i              # printout [D=0]: >0:warning flags, mZ, TB(Ms), for each iteration: mu(Ms), m3(Ms), Mz;  >1:finetuning vars, flags for negative mass squared scalar;  >2:subiterations,tachyon details" %(int(par.get('PRINTOUT',0))))
    if 'QEWSB' in par:
        outs.append("  4   %8.3e    # QEWSB: change electroweak symmetry breaking scale" %(float(par.get('QEWSB'))))
    if 'INCLUDE_2_LOOP_SCALAR_CORRECTIONS' in par:
        outs.append("  5   %i              # INCLUDE_2_LOOP_SCALAR_CORRECTIONS [D=?]" %(int(par.get('INCLUDE_2_LOOP_SCALAR_CORRECTIONS'))))
    if 'PRECISION' in par:
        outs.append("  6   %i              # PRECISION [D=?] : number of significant figures in SLHA output" %(int(par.get('PRECISION'))))
    if 'numHiggsLoops' in par:
        outs.append("  7   %i              # numHiggsLoops [D=?]: number of loops in REWSB/mh calculation" %(int(par.get("numHiggsLoops"))))
    if 'forceSlha1' in par:
        outs.append(" 10   %i              # forceSlha1 [D=?]: if =1, tries to force output into SLHA 1 format" %(int(par.get('forceSlha1'))))
    outs.append("#")
    
    if not keepempty and len(outs)==2: return []
    return outs

        
# ##########
def slhafrompardict(par1={}, par2={}, calculator=[], mode=[], fn=''):
    # TODO: implement more general interface (to allow mSUGRA etc.) 

    par = dict(par1)
    par.update(par2)  
    
    #cat suspect2_lha.in | myprepend.py "    outs.append(\"" | myappend.py "\")"
    outs = []
    
    outs.append("# SUSY Les Houches Accord 2.0 - example input file for SUSPECT ver >= 2.4")
    outs.append("Block MODSEL  # Select model (with the second parameter):")
    outs.append("#            General MSSM (arbitrary soft terms) at low scale input:    0")
    outs.append("#            SUGRA (!includes non-univ. soft terms, def. in block EXTPAR):  1")
    outs.append("#            GMSB                                     : 2")
    outs.append("#            AMSB                                    :  3")
    outs.append("#            Bottom-up RGE for general MSSM input at EWSB scale: -1")
    outs.append("#            (a specific SuSpect option)")
    #outs.append("    1    1   #  mSUGRA")
    outs.append("    1    0   #  pMSSM")
    outs.append("#")

    # ##########
    if 'susyhit' in calculator or 'suspect' in calculator: outs += susyhitblockfrompardict(par=par)

    
    # ##########
    if 'softsusy' in calculator: outs += softsusyblockfrompardict(par=par)

    outs.append("Block SMINPUTS   # Standard Model inputs (if any undefined, defaut values taken)")
    outs.append("     1     127.934  # alpha_em^(-1)(MZ) SM MSbar")
    outs.append("#    2     1.16639d-5  # G_F")
    outs.append("     3     0.1172  # alpha_s(mZ) SM, MSbar")
    outs.append("#    4     91.187   # mZ pole mass")
    outs.append("     5     4.25    # Mb(mb) SM MSbar")
    outs.append("     6     %5.1f      # Mtop(pole)" %(par.get('t',172.5)))  # 2014-01-21
    outs.append("     7     1.7771     # Mtau(pole)")
    outs.append("#")
    outs.append("Block MINPAR  # specific model input parameters")
    outs.append("#   input for SUGRA models (! comment (#) all other (GMSB,AMSB) lines):")
    outs.append("#     1    100.    # m0")
    outs.append("#     2    250.    # m1%2")
    outs.append("#     5   -100.      # A0")
    if 'TBatZ' in par: outs.append("      3    %5.2f    # tanbeta(MZ)" %(par['TBatZ']))
    outs.append("#     4      1.0   # sign(mu)")
    outs.append("#   input for GMSB models (! comment (#) all other (mSUGRA,AMSB) lines):")
    outs.append("#     1    100.d3  # Lambda_susy")
    outs.append("#     2    200.d3  # Lambda_mess")
    outs.append("#     3    10     # tanbeta(MZ)")
    outs.append("#     4    1.      # sign(MU)")
    outs.append("#     5    1       # Nl_mes")
    outs.append("#     6    1       # Nq_mes")
    outs.append("#   input for AMSB models (! comment (#) all other (mSUGRA,GMSB) lines):")
    outs.append("#     1    450.    # m0")
    outs.append("#     2    60.d3   # M_3%2 gravitino mass")
    outs.append("#     3    10.     # tanbeta(MZ)")
    outs.append("#     4    1.      # sign(MU)")
    outs.append("#     5    1.      # cQ  : weight of m0 for Q_L (3rd gen.) doublet")
    outs.append("#     6    1.      # cuR : weight of m0 for u_R")
    outs.append("#     7    1.      # cdR : weight of m0 for d_R")
    outs.append("#     8    1.      # cL  : weight of m0 for L (1st, 2d gen.) doublet")
    outs.append("#     9    1.      # ceR : weight of m0 for e_R (1st, 2d gen.)")
    outs.append("#     10   1.      # cHu : weight of m0 for Hu")
    outs.append("#     11   1.      # cHd : weight of m0 for Hd")
    outs.append("#")
    outs.append("Block EXTPAR  # general MSSM input (! IF uncommented, values replace MINPAR ones)")
    """
    outs.append("#         0      4.65294922E+02   # EWSB_scale")
    outs.append("#         1     2.5E+02   # M_1")
    outs.append("#         2     2.5E+02   # M_2")
    outs.append("#         3     2.5E+02   # M_3")
    outs.append("#        11    -1E+02   # A_t")
    outs.append("#        12    -1E+02   # A_b")
    outs.append("#        13    -2.51737263E+02   # A_tau")
    outs.append("#        14    -6.77437438E+02   # A_u")
    outs.append("#        15    -8.59633345E+02   # A_d")
    outs.append("#        16    -2.53493493E+02   # A_e")
    outs.append("#        23     6.44045685E+02   # mu(EWSB)")
    outs.append("#        26     1.570838901E+02   # MA_pole")
    outs.append("#        25     1.00000000E+01   # tanbeta(MZ)")
    outs.append("#        31     1.E+02   # M_eL")
    outs.append("#        32     1.E+02   # M_muL")
    outs.append("#        33     1.94701043E+02   # M_tauL")
    outs.append("#        34     1.35972069E+02   # M_eR")
    outs.append("#        35     1.35972069E+02   # M_muR")
    outs.append("#        36     1.33500446E+02   # M_tauR")
    outs.append("#        41     5.45887520E+02   # M_q1L")
    outs.append("#        42     5.45887520E+02   # M_q2L")
    outs.append("#        43     4.97055448E+02   # M_q3L")
    outs.append("#        44     5.27854642E+02   # M_uR")
    outs.append("#        45     5.27854642E+02   # M_cR")
    outs.append("#        46     4.21596092E+02   # M_tR")
    outs.append("#        47     5.25761034E+02   # M_dR")
    outs.append("#        48     5.25761034E+02   # M_sR")
    outs.append("#        49     5.22462473E+02   # M_bR")
    """
    
    if 'EWSB_scale' in par:
        if par['EWSB_scale'] == -1: outs.append("         0      %i   # EWSB_scale" %(par['EWSB_scale']))
        else: outs.append("         0      %10.3f   # EWSB_scale" %(par['EWSB_scale']))
    else: 
        outs.append("#         0      4.65294922E+02   # EWSB_scale")
    outs.append("         1     %10.3f   # M_1" %(par['M1']))
    outs.append("         2     %10.3f   # M_2" %(par['M2']))
    outs.append("         3     %10.3f   # M_3" %(par['M3']))
    outs.append("        11     %10.3f   # A_t" %(par['At']))
    outs.append("        12     %10.3f   # A_b" %(par['Ab']))
    outs.append("        13     %10.3f   # A_tau" %(par['AT']))
    outs.append("        14     %10.3f   # A_u" %(par['Au']))
    outs.append("        15     %10.3f   # A_d" %(par['Ad']))
    outs.append("        16     %10.3f   # A_e" %(par['Ae']))
    outs.append("        23     %10.3f   # mu(EWSB)" %(par['mu']))
    outs.append("        26     %10.3f   # MA_pole" %(par['mA']))
    if 'TB' in par: outs.append("        25     %10.3f   # tanbeta(EW)" %(par['TB']))  # 2014-01-21 the if and MZ->EW
    outs.append("        31     %10.3f   # M_eL" %(par['L1']))
    outs.append("        32     %10.3f   # M_muL" %(par['L2']))
    outs.append("        33     %10.3f   # M_tauL" %(par['L3']))
    outs.append("        34     %10.3f   # M_eR" %(par['eR']))
    outs.append("        35     %10.3f   # M_muR" %(par['mR']))
    outs.append("        36     %10.3f   # M_tauR" %(par['TR']))
    outs.append("        41     %10.3f   # M_q1L" %(par['Q1']))
    outs.append("        42     %10.3f   # M_q2L" %(par['Q2']))
    outs.append("        43     %10.3f   # M_q3L" %(par['Q3']))
    outs.append("        44     %10.3f   # M_uR" %(par['uR']))
    outs.append("        45     %10.3f   # M_cR" %(par['cR']))
    outs.append("        46     %10.3f   # M_tR" %(par['tR']))
    outs.append("        47     %10.3f   # M_dR" %(par['dR']))
    outs.append("        48     %10.3f   # M_sR" %(par['sR']))
    outs.append("        49     %10.3f   # M_bR" %(par['bR']))


    if fn:
        f = open(fn,'w')
        for out in outs: f.write('%s\n' %(out))
        f.close()
        
    if 'ret' in mode or fn=='': return outs

# ##########
