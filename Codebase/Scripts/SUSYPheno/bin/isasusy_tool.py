import os, time, random
from isasusy_status import isaout_investigate

def isasusy_make_in(par={}, mode=[], fnID='', fnIDtmp='', save=1):
    if fnIDtmp == '': fnIDtmp = fnID
    
    outs = []
    outs.append('%s.out' %(fnIDtmp))  # 2014-03-10: hack to allow long output names: use tmp here, rename to the original ones outside
    outs.append('%s.wig' %(fnIDtmp))
    outs.append('%7.3f' %(par.get('mtop',172.5)))

    outs.append('%9.3f  %9.3f  %9.3f  %9.3f    %s' %(par.get('gl',par.get('M3')), par['mu'], par['mA'], par['TB'], '  # gl  mu  mA  tanbeta'))
    
    outs.append('%9.3f  %9.3f  %9.3f  %9.3f  %9.3f    %s' %(par['Q1'], par['dR'], par['uR'], par['L1'], par['eR'], '  # Q1  dR  uR  L1  eR'))
    
    outs.append('%9.3f  %9.3f  %9.3f  %9.3f  %9.3f  %9.3f  %9.3f  %9.3f    %s' %(par['Q3'], par['bR'], par['tR'], par['L3'], par['TR'], par['At'], par['Ab'], par['AT'], '  # Q3  bR  tR  L3  TR   A_t  A_b  A_tau'))
    
    outs.append('%9.3f  %9.3f  %9.3f  %9.3f  %9.3f    %s' %(par['Q2'], par['sR'], par['cR'], par['L2'], par['mR'], '  # Q2  sR  cR  L2  mR'))
    
    outs.append('%9.3f  %9.3f    %s' %(par['M1'], par['M2'], '  # M1  M2'))
    outs.append('/')

    
    if save: # default is to save
        f = open('%s_in' %(fnID), 'w')
        for out in outs: f.write('%s\n' %(out))
        f.close()

    if 'ret' in mode: return outs
    


def isasusy_execute(fnID, opt={}, par={}):

    DO_ISASUSY = opt.get('isasusy',1)
    DO_ISA2LHA = opt.get('isa2lha',2)
    DO_ISAPLAY = opt.get('isaplay',1)
    DO_DARKSUSY = opt.get('darksusy',0)
    
    DO_isasusy_status = opt.get('status',1)  # 2014-03-12

    prevar = opt.get('prevar','')
    preval = opt.get('preval','')

    isasusy_x = "isasusy.x"
    isa2lha_py = "isa2lha.py"
    isaplay_py = 'isaPlay.py'
    darksusy_py = 'darksusy.py'

    fn_out = fnID + '.out'
    fn_wig = fnID + '.wig'
    fn_lha = fnID + '.lha'

    if par:
        fnIDtmp = 'tmp.isasusy_%s' %(str(random.randrange(0,1e6)).zfill(6))
        isasusy_make_in(par=par, fnID=fnID, fnIDtmp=fnIDtmp)
    

    if DO_ISASUSY:
        for fn in [fn_out, fn_wig]:
            if os.path.exists(fn): os.remove(fn)
        cmd = '%s < %s_in | grep -v ENTER' %(isasusy_x, fnID)
        os.system(cmd)
        if par:  # 2014-03-10: fix to allow long filenames (how come this was not implemented earlier?)
            cmd = 'mv %s.out %s ; mv %s.wig %s' %(fnIDtmp, fn_out, fnIDtmp, fn_wig)
            os.system(cmd)

    if DO_isasusy_status:
        fn_tlt = fn_out.replace('.out','') + '_status.tlt'
        isaout_investigate(fn_out=fn_out, prevar=prevar, preval=preval, fn_tlt=fn_tlt)

            
    
    if DO_ISA2LHA: 
        cmd = "%s  %s  -out %s  -lha %s  --mssm  --official  -vb 0" %(isa2lha_py, fn_wig, fn_wig.replace("wig","out"),fn_wig.replace("wig","lha")) 
        os.system(cmd)

        if DO_ISA2LHA >= 2:  #then make an 'slha' file, which also contains the decays (BAD NOTATION: .slha/.lha  = with/without decays)
            cmd = "%s  %s  -out %s  -lha %s  --mssm  --official  --decay  -vb 0" %(isa2lha_py, fn_wig, fn_wig.replace("wig","out"),fn_wig.replace("wig","slha")) 
            os.system(cmd)


        cmd = "isaPlay.py  %s  -leptonmode 0 -quarkmode 0 -maxmass 800 -brmin 0.001 -part 3 --qsep  >  %s" %(fn_wig, fn_wig.replace("wig","der"))
        cmd = "isaPlay.py  %s  -leptonmode 0 -quarkmode 0 -maxmass 800 -brmin 0.00001 -part 3 --qsep  >  %s" %(fn_wig, fn_wig.replace("wig","der2"))

        #cmd = "isaPlay.py  %s  -leptonmode 0 -quarkmode 0 -maxmass 800 -brmin 0.001 -part 3 --qsep -notpart h > %s" %(fn_wig, fn_wig.replace("wig","der")) 

        os.system(cmd)


    # ISAPLAY: make derived and mass file (for easy reading)
    if DO_ISAPLAY:
        cmd = "isaPlay.py  %s  >  %s" %(fn_wig, fn_wig.replace(".wig",".der"))
        os.system(cmd)
        cmd = "isaPlay.py  %s  --masses  >  %s" %(fn_wig, fn_wig.replace(".wig",".mas"))
        os.system(cmd)
        cmd = "isaPlay.py  %s  -maxmass 1200 >  %s" %(fn_wig, fn_wig.replace(".wig",".der2"))
        os.system(cmd)


    # DARKSUSY: to check against constraints
    if DO_DARKSUSY:
        cmd = "echo %s | %s  --plha  -rundirtag %s" %(fn_lha, darksusy_py, TIMETAG)  # With unique TIMETAG can run several darksusy in same subdir 
        os.system(cmd)


    if 1: # make the .txt copy here
        cmd = "cp -p %s %s" %(fn_lha.replace("lha","wig"), fn_lha.replace("lha","txt"))
        os.system(cmd)

