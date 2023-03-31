#!/usr/bin/env python
'''
Program : lha2isa.py
Version : 1.0
Author  : b.k.gjelsten@fys.uio.no
Description : 

TO FIX

x Need to deal with the decays of antiparticles

o Note also the hack for the top mass ... maybe look up the value

o There are some negative widths ... (from slha file?)

o Finally do the EWSB parameters at the end (can make optional)

o decays are not ordered after particle code here while they are in the original: simple to fix


-- NOTE: tried to use pyslha to write the wig file, but did not work out of the box

'''




import sys,os
import bkgjelstenArgReader 
from kilelib import ReadPipeOrFile, FilenameReplace, WriteToFile
import pyslha_edited as pyslha
import math
#from isa2lha import isa2pdg  # problem with importing 

# ##################################################### GLOBAL METHODS
# ##################################################### GLOBAL METHODS
# ##################################################### GLOBAL METHODS


isaIDs = list(range(401,457+1)) + list(range(203,207+1)) + [458] + [6,12]
badisaIDs = [438,440,442,444,446,448, 458]

# =========================================================================== #  isa2pdg

def PdgID_of_AntiPart(pdg):
    apdg = abs(pdg)

    # standard cases
    if apdg in list(range(1000001,1000006+1)) + list(range(2000001,2000006+1)) + list(range(1000011,1000016+1)) + list(range(2000011,2000016+1)) + [1000024,1000037] + [37] + list(range(1,6+1)) + list(range(11,16+1)) + [24] + [211]:
        return -pdg

    # its own antiparticle
    if pdg in [1000021, 1000022,1000023,1000025,1000035, 1000039, 25, 35, 36, 21,22,23, 111]:
        return pdg

    # should not land here
    out = 'PdgID_of_AntiPart::FATAL: unknown pdg: %i ' %(pdg)
    sys.exit(out)

    
# #####
isa2pdg = {
    401: 1000001,  #dL
    402: 1000002, 
    403: 1000003, 
    404: 1000004, 
    405: 1000005, 
    406: 1000006,
    407: -1000001,  #dLbar
    408: -1000002, 
    409: -1000003, 
    410: -1000004, 
    411: -1000005, 
    412: -1000006,

    413: 2000001,  #dR
    414: 2000002, 
    415: 2000003, 
    416: 2000004, 
    417: 2000005, 
    418: 2000006, 
    419: -2000001,  #dRbar
    420: -2000002, 
    421: -2000003, 
    422: -2000004, 
    423: -2000005, 
    424: -2000006,
    
    425: 1000011,  #eL -
    426: 1000012,  #ve
    427: 1000013,  #mL -
    428: 1000014,  #vm
    429: 1000015,  #TL -
    430: 1000016,  #vT
    431: -1000011,  #eL bar +
    432: -1000012,  #ve 
    433: -1000013,  #mL +
    434: -1000014,  #vm
    435: -1000015,  #TL +
    436: -1000016,  #vT

    
    437: 2000011,  #eR -
#   438: 2000012,  #veR  :)
    439: 2000013,  #mR -
#   440: 2000014,  #vmR
    441: 2000015,  #TR -
#   442: 2000016,  #vTR
    443: -2000011,  #eR bar +
#   444: -2000012,  #veR
    445: -2000013,  #mR +
#   446: -2000014,  #vmR
    447: -2000015,  #TR +
#   448: -2000016,  #vTR

    449: 1000021,  #gl

    450: 1000022,  #N1
    451: 1000023,  #N2
    452: 1000025,  #N3
    453: 1000035,  #N4
    454: 1000024,  #C1 +
    455: 1000037,  #C2 +
    456: -1000024,  #C1 -
    457: -1000037,  #C2 -

    458: 1000039,  #gravitino (G32)  # 2013-06-26: was 10000039 (on too many zeros)

    203: 25,  #h0 (NB SMhiggs has code 25)  
    #203: 26,  #h0 (NB SMhiggs has code 25)  # 2012-06-11: 26 didn't work with herwig++ 
    204: 35,  #H0
    205: 36,  #A0   # was always correct
    206: 37,  #H+
    207: -37,  #H-


    # Standard Model: 

    1: 1,  #d
    2: 2,  #u
    3: 3,  #s
    4: 4,  #c
    5: 5,  #b
    6: 6,  #t
    7: -1,  #dbar
    8: -2,  #ubar
    9: -3,  #sbar
    10: -4,  #cbar
    11: -5,  #bbar
    12: -6,  #tbar

    121: 11,  #e -
    122: 12,  #ve  
    123: 13,  #m -
    124: 14,  #vm  
    125: 15,  #T -
    126: 16,  #vT  
    127: -11,  #ebar +    # 2012-06-28: sign was wrong 
    128: -12,  #vebar     #  -"-
    129: -13,  #mbar +    #  -"-
    130: -14,  #vmbar     #  -"-
    131: -15,  #Tbar +    #  -"-
    132: -16,  #vTbar     #  -"-


    13: 21,  #gluon
    59: 22,  #gamma
    200: 23,  #Z
    198: 24,  #W+
    199: -24, #W-

    21: 111,  #PI0
    38: 211,  #PI+
    30: -211, #PI-  

    }


def pdg2isa(pdg):
    keys = list(isa2pdg.keys())
    for key in keys:
        if isa2pdg[key] == pdg: return key

    print('Warning::pdg2isa   pdg not in isa2pdg: %i' %(pdg))
    return 0
    

class lha2isa: 

    def __init__(s, cmd=[]):
        s.cmd = cmd
        s.myname = sys.argv[0].split('/').pop()
        s.dict = {}
        s.VB = 1
        s.HOME = os.getenv('HOME')
        #s.dir0 = '%s/XXX' %(s.HOME)
        s.dir0 = ''

        s.dict = {}
        
        # s.dict['test'] = 'testval'
        # s.dict['test2'] = 'testval2'
        # s.dict['testI'] = 3
        # s.dict['testF'] = 4.34

        s.pipe = 1
        s.fnfn_lha = ''   # file with filenames (either this or pipe)
        s.fns_lha = []  # need to be set
        s.fns_wig = []  # derived (or set?)

        s.makewig = 1

        s.dict['fn_repl_lha2wig']  = ['_susyhit_slha.out','.wig']
        s.dict['fn_repl_lha2wig'] += ['.slha','.wig']
        s.dict['fn_repl_lha2wig'] += ['.lha','.wig']
        s.dict['fn_repl_lha2wig'] += ['slha','wig']
        s.dict['fn_repl_lha2wig'] += ['lha','wig']

        s.makeout = 1

        s.dict['fn_repl_lha2out']  = ['_susyhit_slha.out','.out']
        s.dict['fn_repl_lha2out'] += ['.slha','.out']
        s.dict['fn_repl_lha2out'] += ['.lha','.out']

        s.fn_wig_in = ''  # can specify on command line, otherwise will use a replacement formular (specified above) 
        s.fn_out_in = ''

        s.forcedir = ''

        s.warn = []
        s.fn_warn = 'warnings.txt'
        if s.dir0: s.fn_warn = '%s/%s' %(s.dir0, s.fn_warn)
        
        if 'ReadArg' in s.cmd: s.ReadArg()
        if 'Main' in s.cmd: s.Main()
        
    

    # ##################################################### CLASS METHODS
    # ##################################################### CLASS METHODS
    # ##################################################### CLASS METHODS
    # ##########
    def showHelp(s):
        print(' Usage: %s [options]' %(s.myname))
        print('        echo myfile.slha | lha2isa.py')
        print('        lha2isa.py -f fn_with_list_of_slhas  # NOT IMPLEMENTED')
        print('        echo myfile.slha | lha2isa.py  -fnwig a.wig  -fnout a.out')
        print('        ls *.slha | lha2isa.py')


    # ##########
    def DumpWarnings(s):
        f = open(s.fn_warn,'w')
        for out in s.warn: f.write('%s\n' %(out))
        f.close()

        
    # ##########
    def translate_lha2out(s, lha, fnout='', ret=0):
        # only superminimal version is made
        # the reading algo, ReadIsaout, is flexible enough to understand and not fail
        outs = []

        if 'EXTPAR' not in list(lha[0].keys()):
            print('FATAL::lha2isa::translate_lha2out   No EXTPAR in lha friend of %s  ... could rewrite translate_lha2out to instead take (most) from MSOFT' %(fnout))
    
        msoft = lha[0]['MSOFT'].entries
        extpar = lha[0]['EXTPAR'].entries
        hmix = lha[0]['HMIX'].entries
        minpar = lha[0]['MINPAR'].entries

        TANB_Z     = minpar.get(3,-1)   # at MZ
        TANB_INPUT = extpar.get(25,-1)  # at EW scale (Q)
        TANB_EW    = hmix.get(2,-1)     # at EW scale (Q)  .. I think TANB_INPUT and EW are identical, while TANB_Z is different

        if TANB_Z > 0: TANB_use = TANB_Z
        elif TANB_INPUT > 0: TANB_use = TANB_INPUT
        elif TANB_EW > 0: TANB_use = TANB_EW
        else:
            TANBuse = -1
            print('warning::slha2wig  DID not find TANBETA value, set to unphysical -1')

            
        outs.append(' TAN(BETA) = %9.3f' %(TANB_use))
                
        MU    = extpar[23]   # or hmix[1]
        MGLSS = extpar[3]    # this is M_3 ...
        mA    = extpar[26]   
        outs.append(' M(GLSS)   = %9.3f   MU        = %9.3f   M(HA)  = %9.3f' %(MGLSS, MU, mA))
        
        M1 = extpar[1]  # or msoft[1]
        M2 = extpar[2]  # or msoft[2]
        outs.append(' M_1       = %9.3f   M_2       = %9.3f' %(M1,M2) )
        

        # Finalise:
        outs.append('')
        outs.append(' PARENT (this is keyword for libISAWIG::ReadIsaout)')  # <-- this is what triggers exit in ReadIsaout

        # TO FILE 
        if fnout: WriteToFile(fn=fnout, outs=outs, VB=s.VB)


        # RETURN
        if ret: return outs
        

    # ##########
    def translate_lha2wig(s, lha, fnwig='', ret=0): 
        outs = []

        if s.VB>2: print(len(isaIDs))

        # 0 INIT
        massobj = lha[0].get('MASS',{})
        massdict = massobj.entries
        nomass = -9e9

        particledict = lha[1]


        # 1 PARTICLE SECTION
        outs.append('%4i' %(len(isaIDs)))
        for iisaID in range(len(isaIDs)):

            isaID = isaIDs[iisaID]
            # if iisaID > 3: continue
            # if s.VB>1: print isaID

            if isaID in badisaIDs:
                outs.append("%5i %11.4f    %11.5e" %(isaID, 0, 1e30))
                continue

            pdgID = isa2pdg.get(isaID,0)
            apdgID = abs(pdgID)
            if pdgID == 0:
                print('WARNING: isaID: %i' %(isaID))
                

            #if mass == nomass:
            #    s.warn.append('WARNING::translate_lha2wig  No mass for isaID %i   pdgID %i' %(isaID, pdgID))
            #    print s.warn[-1]


            if apdgID in [1000022]: # fragile
                mass = massdict.get(apdgID,nomass)
                lifetime = 1e30  # tja
                
            else:
                particle = particledict[apdgID]  # NB: for antiparticle will need to 'invert'
                mass = particle.mass
                width = particle.totalwidth
                if width == 0: lifetime = 1e30  # 2014-01-20
                else: lifetime = 6.582e-25 / width
                
            if pdgID in [-6,6]:
                mass = 172.5   # ever? # fragile


            # print isaID
            # print mass, type(mass)
            # print width
            out = "%5i %11.4f    %11.5e" %(isaID, mass, lifetime)
            outs.append(out)
    



        # 2 DECAY SECTION
        for iisaID in range(len(isaIDs)):
            isaID = isaIDs[iisaID]

            if isaID in badisaIDs:
                outs.append('   0')
                continue

            
            pdgID = isa2pdg.get(isaID,0)
            apdgID = abs(pdgID)

            isAntiPart = False
            if pdgID < 0: isAntiPart = True

                

            try:
                particle = particledict[apdgID]
            except:
                # Necessary for LSP
                outs.append('%4i' %(0))
                continue
            
            decays = particle.decays
            ndecays = len(decays)
            outs.append('%4i' %(ndecays))

            #print 'DEBUG  %3i  %4i  %8i' %(iisaID, isaID, pdgID)
            for idecay in range(len(decays)):
                decay = decays[idecay]

                whatisthis = 0
                if isaID in [6,12]: whatisthis = 100
                
                out = ' %5i %15.8e   %3i' %(isaID, decay.br, whatisthis)
                for ida in range(decay.nda):
                    thispdgID = decay.ids[ida]

                    if isAntiPart: thispdgID = PdgID_of_AntiPart(thispdgID)
                
                    thisisaID = pdg2isa(thispdgID)
                    out += ' %5i' %(thisisaID)
                for idum in range(5-decay.nda): 
                    out += ' %5i' %(0)

                outs.append(out)


        # 3 INPUT SECTION
        tanbeta = lha[0]['MINPAR'].entries.get(3,-1)  # tanbeta_Z: this one is not filled in slha from susyhit
        if tanbeta < 0: tanbeta = lha[0]['HMIX'].entries.get(2,-1)  # tanbeta_EW: ... can deviate from input value
        if tanbeta < 0: tanbeta = lha[0]['EXTPAR'].entries.get(25,-1)  # tanbeta_EW: ... can deviate from input value
        if tanbeta < 0:
            print('Warning::slha2wig  tanbeta not found. set to -1')

        alphah = lha[0]['ALPHA'].entries
        outs.append('    %12.8f    %12.8f' %(tanbeta, alphah))

        nmix = lha[0]['NMIX']
        for ii in [1,2,3,4]:
            outs.append('    %12.8f    %12.8f    %12.8f    %12.8f' %(nmix.entries[ii][1],nmix.entries[ii][2],nmix.entries[ii][3],nmix.entries[ii][4]))

        vmixs = lha[0]['VMIX'].entries
        outs.append('    %12.8f    %12.8f    %12.8f    %12.8f' %(vmixs[1][1],vmixs[1][2],vmixs[2][1],vmixs[2][2]))
        
        umixs = lha[0]['UMIX'].entries
        outs.append('    %12.8f    %12.8f    %12.8f    %12.8f' %(umixs[1][1],umixs[1][2],umixs[2][1],umixs[2][2]))
        
        thetat   = math.acos(lha[0]['STOPMIX'].entries[1][1])
        thetab   = math.acos(lha[0]['SBOTMIX'].entries[1][1])
        thetatau = math.acos(lha[0]['STAUMIX'].entries[1][1])
        outs.append('    %12.8f    %12.8f    %12.8f' %(thetat, thetab, thetatau))

        a_t   = lha[0]['AU'].entries[3][3]
        a_b   = lha[0]['AD'].entries[3][3]
        a_tau = lha[0]['AE'].entries[3][3]
        outs.append(' %15.8f %15.8f %15.8f' %(a_t, a_b, a_tau))

        mu = lha[0]['HMIX'].entries[1]
        outs.append(' %15.8f' %(mu))

        outs.append('    T')



        # TO FILE 
        if fnwig: WriteToFile(fn=fnwig, outs=outs, VB=s.VB)
            

        # RETURN
        if ret: return outs

        

    # ##########
    def Main(s):
        # if s.VB: print "INFO::%s  Main" %(s.myname)
        # for key in s.dict.keys(): print 'dict:  %-10s  %s' %(key, s.dict[key])

        s.fns_lha = ReadPipeOrFile(pipe=s.pipe, f=s.fnfn_lha)
        # for out in s.fns_lha: print out

        if s.VB>0: print("INFO  Number of slha files: %i" %(len(s.fns_lha)))
        if len(s.fns_lha) == 0: print('lha2isa: Zero slha files read. No action. (pipe=%i ; f=%s)' %(s.pipe, s.fnfn_lha))

        for ilha in range(len(s.fns_lha)):
            fn_lha = s.fns_lha[ilha]

            #print s.fns_lha
            lha = pyslha.readSLHAFile(fn_lha)

            #pyslha.writeISAWIGFile(fn_wig, blocks=lha[0], decays=lha[1])  # not out of the box

            # This small structure ensures flexibility in getting the output dir right
            outdir = ''
            if '/' in fn_lha:
                w = fn_lha.split('/')
                w.pop()
                for iw in range(len(w)):
                    #if iw > 0: outdir += '/'
                    ww = w[iw]
                    outdir += '%s/' %(ww)
            if s.forcedir: outdir = s.forcedir + '/'
                    
            if s.makewig:
                if s.fn_wig_in: fn_wig = outdir + s.fn_wig_in
                else: fn_wig = FilenameReplace(fn_lha, repl=s.dict['fn_repl_lha2wig'], safeapp='.wig')
                s.translate_lha2wig(lha, fnwig=fn_wig)


            if s.makeout:
                if s.fn_out_in: fn_out = outdir + s.fn_out_in
                else: fn_out = FilenameReplace(fn_lha, repl=s.dict['fn_repl_lha2out'], safeapp='.out')
                s.translate_lha2out(lha, fnout=fn_out)
        



    # ##########
    def ReadArg(s): 

        # ################################### ARGUMENT READING
        Arg = bkgjelstenArgReader.ArgReader(sys.argv, VB=0)

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

        if Arg.has('--pipe'):
            s.pipe = 1

        if Arg.hasget(['-fnlha','-f']):
            s.fnfn_lha = Arg.val()
            s.pipe = 0


        if Arg.hasget(['-fnwig','-fn_wig']):
            s.fn_wig_in = Arg.val()

        if Arg.hasget(['-fnout','-fn_out']):
            s.fn_out_in = Arg.val()


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
    t = lha2isa(cmd=['ReadArg','Main'])
############################## 

