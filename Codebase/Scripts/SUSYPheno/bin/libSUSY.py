from math import sqrt


pdgId2txt_plain = {1:'d',2:'u',3:'s',4:'c',5:'b',6:'t',
              11:'e',12:'ve',13:'m',14:'vm',15:'T',16:'vT',
              21:'g',22:'y',23:'Z',24:'W',25:'h',35:'H',36:'A',37:'Hc',
              1000001:'dL',1000002:'uL',1000003:'sL',1000004:'cL',1000005:'b1',1000006:'t1',
              2000001:'dR',2000002:'uR',2000003:'sR',2000004:'cR',2000005:'b2',2000006:'t2',
              1000011:'eL',1000013:'mL',1000015:'T1',1000012:'veL',1000014:'vmL',1000016:'vTL',
              2000011:'eR',2000013:'mR',2000015:'T2',
              1000021:'gl',
              1000022:'N1',1000023:'N2',1000024:'C1',1000025:'N3',1000035:'N4',1000037:'C2',
              1000039:'G',  # gravitino
              }


txt2pdgId_plain = dict((v,k) for k, v in pdgId2txt_plain.items())

slha_extpar_key2txtslha = {
    0: 'EWSB',
    1: 'M_1',
    2: 'M_2',
    3: 'M_3',
    11: 'A_t',
    12: 'A_b',
    13: 'A_tau',
    14: 'A_u',
    15: 'A_d',
    16: 'A_e',
    23: 'mu(EWSB)',
    25: 'tanbeta(in)',
    26: 'MA_pole',
    31: 'M_eL',
    32: 'M_muL',
    33: 'M_tauL',
    34: 'M_eR',
    35: 'M_muR',
    36: 'M_tauR',
    41: 'M_q1L',
    42: 'M_q2L',
    43: 'M_q3L',
    44: 'M_uR',
    45: 'M_cR',
    46: 'M_tR',
    47: 'M_dR',
    48: 'M_sR',
    49: 'M_bR'
    }

slha_extpar_txtslha2key = dict((v,k) for k, v in slha_extpar_key2txtslha.items())
     
slha_extpar_key2twoch = {
    0: 'SB',
    1: 'M1',
    2: 'M2',
    3: 'M3',
    11: 'At',
    12: 'Ab',
    13: 'AT',
    14: 'Au',
    15: 'Ad',
    16: 'Ae',
    23: 'mu',
    25: 'TB',
    26: 'mA',
    31: 'L1',
    32: 'L2',
    33: 'L3',
    34: 'eR',
    35: 'mR',
    36: 'TR',
    41: 'Q1',
    42: 'Q2',
    43: 'Q3',
    44: 'uR',
    45: 'cR',
    46: 'tR',
    47: 'dR',
    48: 'sR',
    49: 'bR'
    }

slha_extpar_twoch2key = dict((v,k) for k, v in slha_extpar_key2twoch.items())


prosp_twoch2fivekey = {
    
    #'twoch': (fn,ip1,ip2,isq1,isq2)

    # incomplete, need to add many subprocesses

    'gg': ('gg',9,9,0,0)
    
    ,'dg': ('sg',9,9,1,1)
    ,'ug': ('sg',9,9,2,2)
    ,'sg': ('sg',9,9,3,3)
    ,'cg': ('sg',9,9,4,4)
    ,'dbg': ('sg',9,9,-1,-1)
    ,'ubg': ('sg',9,9,-2,-2)
    ,'sbg': ('sg',9,9,-3,-3)
    ,'cbg': ('sg',9,9,-4,-4)

    #,'': ('',9,9,0,0)
    
    ,'gN1':  ('ng',1,9,0,0)
    ,'gN2':  ('ng',2,9,0,0)
    ,'gN3':  ('ng',3,9,0,0)
    ,'gN4':  ('ng',4,9,0,0)
    ,'gC1+': ('ng',5,9,0,0)
    ,'gC2+': ('ng',6,9,0,0)
    ,'gC1-': ('ng',7,9,0,0)
    ,'gC2-': ('ng',8,9,0,0)
    
    ,'N1N1': ('nn',1,1,0,0)
    ,'N1N2': ('nn',1,2,0,0)
    ,'N1N3': ('nn',1,3,0,0)
    ,'N1N4': ('nn',1,4,0,0)
    ,'N2N2': ('nn',2,2,0,0)
    ,'N2N3': ('nn',2,3,0,0)
    ,'N2N4': ('nn',2,4,0,0)
    ,'N3N3': ('nn',3,3,0,0)
    ,'N3N4': ('nn',3,4,0,0)
    ,'N4N4': ('nn',4,4,0,0)
    
    ,'N1C1+': ('nn',1,5,0,0)
    ,'N1C2+': ('nn',1,6,0,0)
    ,'N1C1-': ('nn',1,7,0,0)
    ,'N1C2-': ('nn',1,8,0,0)
    ,'N2C1+': ('nn',2,5,0,0)
    ,'N2C2+': ('nn',2,6,0,0)
    ,'N2C1-': ('nn',2,7,0,0)
    ,'N2C2-': ('nn',2,8,0,0)
    ,'N3C1+': ('nn',3,5,0,0)
    ,'N3C2+': ('nn',3,6,0,0)
    ,'N3C1-': ('nn',3,7,0,0)
    ,'N3C2-': ('nn',3,8,0,0)
    ,'N4C1+': ('nn',4,5,0,0)
    ,'N4C2+': ('nn',4,6,0,0)
    ,'N4C1-': ('nn',4,7,0,0)
    ,'N4C2-': ('nn',4,8,0,0)

    ,'C1+C1-': ('nn',5,7,0,0)
    ,'C1+C2-': ('nn',5,8,0,0)
    ,'C2+C1-': ('nn',6,7,0,0)
    ,'C2+C2-': ('nn',6,8,0,0)


    ,'ee'  :  ('ll',0,9,0,0)
    ,'eLeL':  ('ll',1,9,0,0)
    ,'eReR':  ('ll',2,9,0,0)
    ,'vEvE':  ('ll',3,9,0,0)
    ,'eL+vE': ('ll',4,9,0,0)
    ,'eL-vE': ('ll',5,9,0,0)
    ,'T1T1':  ('ll',6,9,0,0)
    ,'T2T2':  ('ll',7,9,0,0)
    ,'T1T2':  ('ll',8,9,0,0)  # So this is T1+T2-  +  T1-T2+, I guess ??  <--- Need to check 
    ,'vTvT':  ('ll',9,9,0,0)
    ,'T1+vT': ('ll',10,9,0,0)
    ,'T1-vT': ('ll',11,9,0,0)
    ,'T2+vT': ('ll',12,9,0,0)
    ,'T2-vT': ('ll',13,9,0,0)

    }


prosp_fivekey2twoch  = dict((v,k) for k, v in prosp_twoch2fivekey.items())


def pdgId2txt(pdgId,mode=[]):
    apdgId = abs(pdgId)
    if apdgId in pdgId2txt_plain: return pdgId2txt_plain[apdgId]
    else: return '?'


def isSUSY(pdg):
    if  1000001 <=  pdg <= 2000016: return True
    if  1000001 <= -pdg <= 2000016: return True
    return False



def GetMassFromSLHA(fn, pdgId):
    # This is brute mass-getter
    # (more advanced is to create a SLHA object (not yet implemented) and SLHA.getmass(pdgid)

    # Method: Look for 'Block mass', then look for pdgid (until next block)
    
    f = open(fn); lines = f.readlines(); f.close()

    m = -1.

    inBlockMASS = False
    for iL in range(len(lines)):
        L = lines[iL].rstrip()
        w = L.split()

        if not inBlockMASS:
            if len(w) >= 2 and w[0].lower() == 'block' and w[1].lower() == 'mass':
                inBlockMASS = True
            continue

        # when enter here, we are inside Block MASS
        if L.startswith('#'): continue  # skip comment line
        if w[0].lower() == 'block':  # the block mass is ending, break the loop (no mass gotten)
            inBlockMASS = False
            break

        if len(w) < 2: # should not have less than two entries now
            print('WARNING: GetMassFromSLHA: Irregular line in block mass: ""%s""' %(L))
            continue

        thispdgId = int(w[0])
        
        if thispdgId == abs(pdgId):  #Yepp, this is the one I seek
            thism = float(w[1])
            return thism


    # If reach here, the mass was not found
    print('WARNING: GetMassFromSLHA: Did not find mass of %i in file %s   [RETURNING -1]' %(pdgId, fn))
    return -1


prospino_pid_dict = {'nn':{1:1000022, 2:1000023, 3:1000025, 4:1000035, 5:1000024, 6:1000037, 7:-1000024, 8:-1000037} }


def GetPdgIdsFromProspinoSubprocess(sub, mode='abs'):

    id1 = 0
    id2 = 0
    
    if sub.startswith('nn'):
        i1 = int(sub[2:3])
        i2 = int(sub[3:4])
        id1 = prospino_pid_dict['nn'][i1]
        id2 = prospino_pid_dict['nn'][i2]


    if sub.startswith('ng'):
        i1 = int(sub[2:3])
        id1 = prospino_pid_dict['nn'][i1]
        id2 = 1000021

    if sub.startswith('ns'):  # well, ... using here dL for general squark mass
        i1 = int(sub[2:3])
        id1 = prospino_pid_dict['nn'][i1]
        id2 = 1000001    

    if sub.startswith('ll'):
        i = int(sub[2:])
        if i == 1:  id1,id2 = -1000011, 1000011  #eLeL
        if i == 2:  id1,id2 = -2000011, 2000011  #eReR
        if i == 3:  id1,id2 =  1000012, 1000012  #veve
        if i == 4:  id1,id2 = -1000011, 1000012  #eL+veL
        if i == 5:  id1,id2 =  1000011, 1000012  #eL-veL
         
        if i == 6:  id1,id2 = -1000015, 1000015  #T1T1
        if i == 7:  id1,id2 = -2000015, 2000015  #T2T2
        if i == 8:  id1,id2 = -1000015, 2000015  #T1T2
        
        if i == 9:  id1,id2 =  1000016, 1000016  #vTLvTL
        if i == 10: id1,id2 = -1000015, 1000016  #T1+vTL
        if i == 11: id1,id2 =  1000015, 1000016  #T1-vTL
        if i == 12: id1,id2 = -2000015, 1000016  #T2+vTL
        if i == 13: id1,id2 =  2000015, 1000016  #T2-vTL
        
    if id1 == 0 or id2 == 0:
        print('WARNING: GetSparticleFromSubprocess: irregular subprocess %s' %(sub))

    
    if mode == 'abs': return abs(id1),abs(id2)
    
    return id1,id2


def GetPnamesFromProspinoSubprocess(sub, mode='abs'):
    id1, id2 = GetPdgIdsFromProspinoSubprocess(sub,mode)
    nam1 = pdgId2txt_plain[abs(id1)]
    nam2 = pdgId2txt_plain[abs(id2)]
    return nam1,nam2


prosp_twoch2twoword = {
    
    'gg': 'g g'
    
    ,'dg': 'd g'  # standard?
    ,'ug': 'u g'  
    ,'sg': 's g'
    ,'cg': 'c g'
    ,'dbg': 'db g'
    ,'ubg': 'ub g'
    ,'sbg': 'sb g'
    ,'cbg': 'cb g'

    #,'': ('',9,9,0,0)
    
    ,'gN1':  'g N1'
    ,'gN2':  'g N2'
    ,'gN3':  'g N3'
    ,'gN4':  'g N4'
    ,'gC1+': 'g C1+'
    ,'gC2+': 'g C2+'
    ,'gC1-': 'g C1-'
    ,'gC2-': 'g C2-'
    
    ,'N1N1': 'N1 N1'
    ,'N1N2': 'N1 N2' 
    ,'N1N3': 'N1 N3' 
    ,'N1N4': 'N1 N4' 
    ,'N2N2': 'N2 N2' 
    ,'N2N3': 'N2 N3' 
    ,'N2N4': 'N2 N4' 
    ,'N3N3': 'N3 N3' 
    ,'N3N4': 'N3 N4' 
    ,'N4N4': 'N4 N4' 
    
    ,'N1C1+': 'N1 C1+'
    ,'N1C2+': 'N1 C2+'
    ,'N1C1-': 'N1 C1-'
    ,'N1C2-': 'N1 C2-'
    ,'N2C1+': 'N2 C1+'
    ,'N2C2+': 'N2 C2+'
    ,'N2C1-': 'N2 C1-'
    ,'N2C2-': 'N2 C2-'
    ,'N3C1+': 'N3 C1+'
    ,'N3C2+': 'N3 C2+'
    ,'N3C1-': 'N3 C1-'
    ,'N3C2-': 'N3 C2-'
    ,'N4C1+': 'N4 C1+'
    ,'N4C2+': 'N4 C2+'
    ,'N4C1-': 'N4 C1-'
    ,'N4C2-': 'N4 C2-'

    ,'C1+C1-': 'C1+ C1-'
    ,'C1+C2-': 'C1+ C2-'
    ,'C2+C1-': 'C2+ C1-'
    ,'C2+C2-': 'C2+ C2-'


    ,'ee'  :  'e e'
    ,'eLeL':  'eL eL'
    ,'eReR':  'eR eR'
    ,'vEvE':  'vE vE'
    ,'eL+vE': 'eL+ vE'
    ,'eL-vE': 'eL- vE'
    ,'T1T1':  'T1 T1'
    ,'T2T2':  'T2 T2'
    ,'T1T2':  'T1 T2'
    ,'vTvT':  'vT vT'
    ,'T1+vT': 'T1+ vT'
    ,'T1-vT': 'T1- vT'
    ,'T2+vT': 'T2+ vT'
    ,'T2-vT': 'T2- vT'

    }


prosp_twoword2twoch  = dict((v,k) for k, v in prosp_twoch2twoword.items())


# ########### ########### ########## (2014-02-17)
# ########### ########### ########## Analytic formulas for neutralino masses and the neutralino mixing matrix
# ########### ########### ########## Phys Rev D Vol 45, Number 11, 1 June 1992
# ########### ########### ##########   ... DOESN'T WORK .... I get sqrt of negative .. :( 

def neutralino_mass(M1, M2, mu, TB, mode=1):  #
    M1=float(M1)
    M2=float(M2)
    mu=float(mu)
    TB=float(TB)

    VB=1
    # init
    mZ = 91.2
    sinOw2 = 0.23120  # err=0.00015  (at Z-peak)
    cosOw2 = 1-sinOw2

    # sin2beta = 2*tanbeta / (1+tanbeta^2)
    sin2beta = 2.*TB / (1.+TB**2)
    if VB: print('sin2beta: ', sin2beta)
    
    res = {}

    Mp = 5./3 * M1
    M  = M2

    c2 = (Mp*M - mZ**2 - mu**2) - (3./8)*(Mp+M)**2

    c3 = -(1./8)*(Mp+M)**3 + (1./2)*(Mp+M) * (Mp*M - mZ**2 - mu**2) + (Mp+M)*mu**2 + (Mp*cosOw2 + M*sinOw2)*mZ**2 - mu*mZ**2*sin2beta

    c4 = (Mp*cosOw2 + M*sinOw2) * mZ**2*mu*sin2beta - Mp*M*mu**2 + (1./4)*(Mp+M) * ( (Mp+M)*mu**2 + (Mp*cosOw2 + M*sinOw2)*mZ**2 - mu*mZ**2*sin2beta ) + (1./16)*(Mp*M - mZ**2 - mu**2) * (Mp+M)**2 - (3./256)*(Mp+M)**4

    U = -(1./3)*c2**3 - 4*c4
    S = -c3**2 - (2./27)*c2**3 + (8./3)*c2*c4

    D = -4*U**3 - 27*S**2   # for (100,200,300,10) D is negative, -5e44

    if VB: print('c2: ', c2)
    if VB: print('c3: ', c3)
    if VB: print('c4: ', c4)

    if VB: print('U: ',U)
    if VB: print('S: ',S)
    if VB: print('D: ',D)

    z1 = c3**2 + (2./27)*c2**3 - (8./3)*c2*c4
    #if VB: print 'z_real: ',z_real
    #z2 = sqrt(D/27.)
    z = z1 + complex(0,1) * (complex(D/27.))**(1./2)
    if VB: print('z: ',z)
    
    #z = complex( c3**2 + (2./27)*c2**3 - (8./3)*c2*c4,  sqrt(D/27.) )


    if mode==1: # interpretation 1
        if z.real < 0: z3 = -(abs(z.real)**(1./3))   # hack due to problem if negative [since cannot do analytical]
        else: z3 = z.real**(1./3)
        a = 1./2**(1./3) * z3     # -13067507.9508
    else: # interpretation 2
        a = 1./2**(1./3) * (z**(1./3)).real  # 6533753.97539

    if VB: print('a: ',a)


    # Some intermediate values 
    A = sqrt( (1./2)*a - (1./6)*c2 )  # TYPICALLY THIS ONE BECOMES IMAGINARY
    if VB: print('A: ', A)
    B = -(1./2)*a -(1./3)*c2
    if VB: print('B: ', B)
    C = c3 / sqrt(8*a - (8./3)*c2)
    if VB: print('C: ', C)


    # The the neutralino masses
    N1 = - A + sqrt( B + C ) + (1./4)*(Mp+M)
    N2 = + A - sqrt( B - C ) + (1./4)*(Mp+M)
    N3 = - A - sqrt( B + C ) + (1./4)*(Mp+M)
    N4 = + A + sqrt( B - C ) + (1./4)*(Mp+M)

    aN1=abs(N1); aN2=abs(N2); aN3=abs(N3); aN4=abs(N4)
    
    res['N1'] = N1
    res['N2'] = N2
    res['N3'] = N3
    res['N4'] = N4
    
    if VB: print('(M1,M2,mu,TB) = (%.2f,%.2f,%.2f,%.2f) ===>>> (N1,N2,N3,N4)=(%.2f,%.2f,%.2f,%.2f)' %(M1,M2,mu,TB, aN1,aN2,aN3,aN4))
    
    return res

# ##########

