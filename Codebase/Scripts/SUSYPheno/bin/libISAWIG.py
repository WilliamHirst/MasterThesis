import sys, os
from math import cos, sin

################################################################################
################################################################################
# DANGEROUS HACKS:
#  - Expansion of h,H,A,H+ is set 100% to bb regardless of actual decay pattern (no warnings are given)
################################################################################
################################################################################


# (ta1,ta2,snue,snut,eL,eR) for what in the new notation is (T1,T2,vl,vT,lL,lR) [joined sleptons]
# The non-joined leptons e and m are: (eL,eR,mL,mR,ve,vm)

#from ROOT import *  # makes irritating output: registering stepper CashKarpRKF45, etc.
#from ROOT import TH1F, TCanvas, TLegend, gStyle, gROOT  # sufficient?
#from ROOT import cos, sin  # seem to need only these


HOME = "/mn/kvant/u1/borgeg"
StructureISAWIGfiles = HOME+"/bin/StructureISAWIGfiles.py"


from math import fmod

# =============================================================================
# =============================================================================
# =============================================================================
def isGluino(p): return (p in ['gl'])
def isqL(p): return (p in ['~dL','~uL','~sL','~cL'])
def isqR(p): return (p in ['~dR','~uR','~sR','~cR'])
def issq(p): return (p in ['~dL','~uL','~sL','~cL','~dR','~uR','~sR','~cR'])
def isSQ(p): return (p in ['~dL','~uL','~sL','~cL','~dR','~uR','~sR','~cR','~b1','~b2'])
def isN(p): return (p in ['~N1','~N2','~N3','~N4', '~N', '~N_'])
def isC(p): return (p in ['~C1','~C2', '~C', '~C_'])
def isNC(p): return (p in ['~N1','~N2','~N3','~N4','~C1','~C2', '~N','~C','~X','~N_','~C_'])
def isSlepton(p): return (p in ['~lR','~lL','~eR','~eL','~mR','~mL','~sl'])
def isSnu(p): return (p in ['~ve','~vm','~vl','~vT'])
def isStau(p): return (p in ['~T1','~T2'])

def isq(p): return (p in ['d','u','s','c'])
def isQ(p): return (p in ['d','u','s','c','b'])
def isLepton(p): return (p in ['l','e','m'])
def isTau(p): return (p in ['T'])


setGluino = ['gl']
setqL = ['dL','uL','sL','cL']
setqR = ['dR','uR','sR','cR']
setsq = ['dL','uL','sL','cL','dR','uR','sR','cR']
setSQ = ['dL','uL','sL','cL','dR','uR','sR','cR','b1','b2']
setN = ['N1','N2','N3','N4']
setC = ['C1','C2']
setNC = ['N1','N2','N3','N4','C1','C2']

setq = ['d','u','s','c']
setQ = ['d','u','s','c','b']
setj = ['d','u','s','c','b','T']

# setMassiveDaughter = ['Z->qq','Z->bb','Z->QQ', 'W->qq','W->QQ', 'h','H','A','Hp']
setMassiveDaughter = ['Z->qq','Z->bb','Z->QQ', 'W->qq','W->QQ', 'h','H','A','Hp', 'Z->TT']


categAll=['shead','top','pi','exlep','nu','dilep','ditau','extau','y']  #used in Classify

# ################################################### OLD (STILL IN USE)

def addDict_slow(Dicts):
    # adds dictionaries ... must be a more elegant/faster/automatic way??
    # There is: look at addDict
    if len(Dicts) == 0: return {}
    newDict = Dicts[0]
    for iDict in range(1, len(Dicts)):
        zDict = Dicts[iDict]
        for key in list(zDict.keys()):
            if not key in list(newDict.keys()): newDict[key] = zDict[key]
            else: newDict[key] += zDict[key]

    return newDict


def add2dict(a,b):
    # http://stackoverflow.com/questions/1031199/adding-dictionaries-in-python
    return dict( (n, a.get(n, 0)+b.get(n, 0)) for n in set(a)|set(b) )


# =============================================================================
pdg2nam = {

    1: "d",
    2: "u",
    3: "s",
    4: "c",
    5: "b",
    6: "t",

    11: "el",
    13: "mu",
    15: "ta",
    
    12: "vel",
    14: "vmu",
    16: "vta",

    21: "g", 
    22: "y", 
    23: "Z",
    24: "W",
    25: "h",
    26: "h", # 2012-05-06
    35: "H",
    36: "A",   # 2014-01-29 fixed; 36<->37
    37: "H+",  # 2014-01-29 fixed; 36<->37
    1000001: "dL", 
    1000002: "uL", 
    1000003: "sL", 
    1000004: "cL", 
    1000005: "b1", 
    1000006: "t1", 
    2000001: "dR", 
    2000002: "uR", 
    2000003: "sR", 
    2000004: "cR", 
    2000005: "b2", 
    2000006: "t2", 

    1000011: "eL", 
    1000013: "mL", 
    1000015: "T1", 

    2000011: "eR",
    2000013: "mR", 
    2000015: "T2",
    # 2000012: "mR",  # 2012-05-06  Was wrong!! (in use somewhere?)
    # 2000014: "T2",  # 2012-05-06  Was wrong!! (in use somewhere?)
    
    1000012: "ve", 
    1000014: "vm", 
    1000016: "vT", 

    1000021: "gl",
    
    1000022: "N1", 
    1000023: "N2", 
    1000024: "C1", 
    1000025: "N3", 
    1000035: "N4", 
    1000037: "C2", 
}

# =============================================================================
nam2pdg = {
    "Z":  23, 
    "W":  24, 
    "h":  25, 
    "H":  35, 
    "A":  36,  # 2014-01-29 fixed; 36<->37
    "H+": 37,  # 2014-01-29 fixed; 36<->37
    "dL": 1000001, 
    "uL": 1000002, 
    "sL": 1000003, 
    "cL": 1000004, 
    "b1": 1000005, 
    "t1": 1000006, 
    "dR": 2000001, 
    "uR": 2000002, 
    "sR": 2000003, 
    "cR": 2000004, 
    "b2": 2000005, 
    "t2": 2000006, 

    "eL": 1000011, 
    "mL": 1000013, 
    "T1": 1000015, 
    "eR": 2000011, 
    "mR": 2000013, 
    "T2": 2000015, 

    "ve": 1000012, 
    "vm": 1000014, 
    "vT": 1000016, 

    #"muL": 1000013, 
    #"ta1": 1000015, 
    #"muR": 2000013, 
    #"ta2": 2000015, 
    #"snue": 1000012, 
    #"snum": 1000014, 
    #"snut": 1000016, 

    "gl": 1000021, 
    "N1": 1000022, 
    "N2": 1000023, 
    "C1": 1000024, 
    "N3": 1000025, 
    "N4": 1000035, 
    "C2": 1000037, 
}




# =============================================================================
def fname(IDgen,QuarkMode=1,LeptonMode=1):

    IDgen=fIDgen(IDgen,QuarkMode,LeptonMode)
    
    if(IDgen==401): return "dL"
    if(IDgen==402): return "uL"
    if(IDgen==403): return "sL"
    if(IDgen==404): return "cL"
    if(IDgen==405): return "b1"
    if(IDgen==406): return "t1"
    if(IDgen==413): return "dR"
    if(IDgen==414): return "uR"
    if(IDgen==415): return "sR"
    if(IDgen==416): return "cR"
    if(IDgen==417): return "b2"
    if(IDgen==418): return "t2"
    
    if(IDgen==425): return "eL"
    if(IDgen==426): return "ve"
    if(IDgen==427): return "mL"
    if(IDgen==428): return "vm"
    if(IDgen==429): return "T1"
    if(IDgen==430): return "vT"
    if(IDgen==437): return "eR"
    if(IDgen==439): return "mR"
    if(IDgen==441): return "T2"
    
    if(IDgen==449): return "gl"
    
    if(IDgen==450): return "N1"
    if(IDgen==451): return "N2"
    if(IDgen==452): return "N3"
    if(IDgen==453): return "N4"
    
    if(IDgen==454): return "C1"
    if(IDgen==455): return "C2"
    
    if(IDgen==203): return "h"
    if(IDgen==204): return "H"
    if(IDgen==205): return "A"
    if(IDgen==206): return "H+"

    if(IDgen==1): return "d"
    if(IDgen==2): return "u"
    if(IDgen==3): return "s"
    if(IDgen==4): return "c"
    if(IDgen==5): return "b"
    if(IDgen==6): return "t"

    if(IDgen==13): return "g"
    if(IDgen==59): return "y" #"pho"

    if(IDgen==121): return "e"
    if(IDgen==122): return "nue"
    if(IDgen==123): return "mu"
    if(IDgen==124): return "numu"
    if(IDgen==125): return "ta"
    if(IDgen==126): return "nut"

    if(IDgen==198): return "W"
    if(IDgen==200): return "Z"

    if(IDgen==-1): return ""

    if(IDgen==38):  return "PI+-" #38=PI+; 30=PI-
    if(IDgen==21):  return "PI0"

    if(IDgen==458): return "grav"

    # Some extras
    if(IDgen==501): return "qL"
    if(IDgen==502): return "qR"

    if(IDgen==601): return "q"


    # if arrives here, have unknown identity so stops
    res="NBNB: UNKNOWN:"+str(IDgen)
    print(res + " ... ENDING ***********************");
    print(res + " ... ENDING ***********************");
    print(res + " ... ENDING ***********************");
    sys.exit()
    #return res



# =============================================================================
def fIDgen(IDspec, QuarkMode=0, LeptonMode=0):
    
    if(IDspec==407): return 401  # dL,uL,sL,cL
    if(IDspec==408): return 402
    if(QuarkMode==0):  # c contained in u, s contained in d  [more correct: s contained in d]
        if(IDspec==403 or IDspec==409): return 401
        if(IDspec==404 or IDspec==410): return 402
    if(QuarkMode==1):  # c, u, s and d are counted separately
        if(IDspec==403 or IDspec==409): return 403
        if(IDspec==404 or IDspec==410): return 404
    if(QuarkMode==2 or QuarkMode==3):  # q=u,d,c,s
        if(IDspec==401 or IDspec==407): return 501
        if(IDspec==402 or IDspec==408): return 501
        if(IDspec==403 or IDspec==409): return 501
        if(IDspec==404 or IDspec==410): return 501
    if(QuarkMode==3):  # q=u,d,c,s,b
        if(IDspec==405 or IDspec==411): return 501


    if(IDspec==411): return 405
    if(IDspec==412): return 406


    if(IDspec==419): return 413  # dR,uR,sR,cR
    if(IDspec==420): return 414
    if(QuarkMode==0):
        if(IDspec==415 or IDspec==421): return 413
        if(IDspec==416 or IDspec==422): return 414
    if(QuarkMode==1):
        if(IDspec==415 or IDspec==421): return 415
        if(IDspec==416 or IDspec==422): return 416
    if(QuarkMode==2 or QuarkMode==3):
        if(IDspec==413 or IDspec==419): return 502
        if(IDspec==414 or IDspec==420): return 502
        if(IDspec==415 or IDspec==421): return 502
        if(IDspec==416 or IDspec==422): return 502
    if(QuarkMode==3):
        if(IDspec==417 or IDspec==423): return 502

    if(IDspec==423): return 417
    if(IDspec==424): return 418
    
    if(IDspec==431): return 425  # eL,nu_e,muL,nu_mu
    if(IDspec==432): return 426
    if(LeptonMode==0):
        if(IDspec==427 or IDspec==433): return 425
        if(IDspec==428 or IDspec==434): return 426
    if(LeptonMode==1):
        if(IDspec==427 or IDspec==433): return 427
        if(IDspec==428 or IDspec==434): return 428
        
    if(IDspec==435): return 429
    if(IDspec==436): return 430

    if(IDspec==443): return 437  # eR,muR
    if(LeptonMode==0):
        if(IDspec==439 or IDspec==445): return 437
    if(LeptonMode==1):
        if(IDspec==439 or IDspec==445): return 439
    if(IDspec==447): return 441

    if(IDspec==456): return 454
    if(IDspec==457): return 455

    if(IDspec==207): return 206


    if(IDspec==7): return 1  # d,u,s,c
    if(IDspec==8): return 2
    if(QuarkMode==0):
        if(IDspec==3 or IDspec==9):  return 1
        if(IDspec==4 or IDspec==10): return 2
    if(QuarkMode==1):
        if(IDspec==3 or IDspec==9):  return 3
        if(IDspec==4 or IDspec==10): return 4
    if(QuarkMode==2 or QuarkMode==3):
        if(IDspec==1 or IDspec==7):  return 601
        if(IDspec==2 or IDspec==8):  return 601
        if(IDspec==3 or IDspec==9):  return 601
        if(IDspec==4 or IDspec==10): return 601
    if(QuarkMode==3):
        if(IDspec==5 or IDspec==11): return 601
    
    if(IDspec==11): return 5
    if(IDspec==12): return 6

    if(IDspec==127): return 121 # e,nue,mu,numu
    if(IDspec==128): return 122
    if(LeptonMode==0): # count mu=e and nu_mu=nu_e (mSUGRA-like)
        if(IDspec==123 or IDspec==129): return 121
        if(IDspec==124 or IDspec==130): return 122
    if(LeptonMode==1): # count mu and e separately as well as nu_mu and nu_e separately
        if(IDspec==123 or IDspec==129): return 123
        if(IDspec==124 or IDspec==130): return 124

    if(IDspec==131): return 125
    if(IDspec==132): return 126

    if(IDspec==199): return 198

    if(IDspec==30): return 38

    #nonparticles=[438,440,442,444,446,448, 458]
    nonparticles=[438,440,442,444,446,448]
    if(IDspec in nonparticles): return -1

    #if not found and not in nonparticles, then return input value
    return IDspec


# =============================================================================
def pType(name):
    c=name
    if(c=="gl"): return 1
    if(c=="dR" or c=="uR" or c=="dL" or c=="uL" or
       c=="sR" or c=="cR" or c=="sL" or c=="cL" or
       c=="b1" or c=="b2" or c=="t1" or c=="t2"): return 2
    #if(c=="eL" or c=="eR" or c=="muL" or c=="muR" or c=="snue" or c=="snum" or c=="ta1" or c=="ta2" or c=="snut"): return 3
    if(c=="eL" or c=="eR" or c=="mL" or c=="mR" or c=="ve" or c=="vm" or c=="T1" or c=="T2" or c=="vT"): return 3
    if(c=="N1" or c=="N2" or c=="N3" or c=="N4" or c=="C1" or c=="C2"): return 4
    if(c=="h" or c=="H" or c=="A" or c=="H+"): return 5
    if(c=="grav"): return 0

    print("Type not defined for: "+name)
    return 0









# ################################################## NEW


class Incrementer:
    def __init__(self,it = 0):
        self.it = it

    def next(self,txt1='',txt2=''):
        self.it += 1
        return txt1+str(self.it)+txt2
    
# =============================================================================
class FreeObject:
    pass


# =============================================================================
class IsaPlayRes:
    # Sep28: so far only implemented for lepton-only results (mode == 'l')
    # 2012-06-23: implemented for any 'singly-surviving SM particle'; l, q, b, T or y (done to get results for T)
    
    def __init__(self,fn='',playmode='l'):
        self.x = {}
        self.xrel = {}
        self.xTot = 0.
        self.xrelTot = 0.
        self.xAccu = {}
        self.xrelAccu = {}
        self.ok = 0
        
        if fn:
            self.init(fn,playmode=playmode)
            self.calcMore()


    def init(self,fn, playmode='l'):
        if not os.path.exists(fn): 
            print('WARNING  libISAWIG::IsaPlayRes  Non-existing res_isaplay file %s  (continuing)' %(fn)) 
            # self.ok = 0
            return
        #print 'BUT NOW ALL IS FINE'
        self.ok = 1
        f = open(fn); lines = f.readlines(); f.close()

        dummy_xsecAccu = {}  

        #print 'playmode: %s' %(playmode)
        colID = 5  # is default, is 'l'
        if playmode == 'q': colID = 6
        if playmode == 'b': colID = 7
        if playmode == 'T': colID = 8
        if playmode == 'y': colID = 9

        
        for iL in range(len(lines)):
            L = lines[iL]
            if iL < 2: continue
            w = L.strip().split()
            if len(w) < 10:
                print("WARNING IsaPlayRes: too short line:  %s" %(L.strip()))

            Nl = int(w[colID])  # number of leptons (thisnL)
            self.x[Nl] = float(w[1])     # XSEC         for nlep == thisnL
            self.xrel[Nl] = float(w[4])  # relXSEC (BR) for nLep == thisnL

            dummy_xsecAccu[Nl] = float(w[0])


        # Hack since XSEC is given with only 3 decimals while BR is given with 6 decimals in the play files,
        # is useful to redo all XSEC anew and get them from XSECtotal * BR(nL=thisnL) ; since XSECtotal is the largest and hence has the lowest uncertainty
        for Nl in list(self.x.keys()):
            self.x[Nl] = self.xrel[Nl] * dummy_xsecAccu[0]

        del dummy_xsecAccu

            
    
    def calcMore(self):
        # First fill upper lepton multiplicities (with zero)
        for nLep in range(0,9):
            if nLep not in self.x:
                self.x[nLep] = 0.
                self.xrel[nLep] = 0.
                
        # Then ...
        for key in list(self.x.keys()):
            # Fill Aot
            self.xTot += self.x[key]
            self.xrelTot += self.xrel[key]

            # Fill Accu
            self.xAccu[key] = 0.
            self.xrelAccu[key] = 0.
            # For each key, loop over all and add the ones with nLep >= this_nLep
            for key2 in list(self.x.keys()):
                if key2 >= key:
                    self.xAccu[key] += self.x[key2]
                    self.xrelAccu[key] += self.xrel[key2]


# =============================================================================
class OXcategories:
    def __init__(self):
        self.part = {}
        p = self.part
        
        p['g'] = ['gl']
        p['u'] = ['uL','uR']
        p['d'] = ['dL','dR']
        p['s'] = ['sL','sR']
        p['c'] = ['cL','cR']
        p['b'] = ['b1','b2']
        p['t'] = ['t1','t2']
        p['q4'] = p['u']+p['d']+p['s']+p['c']
        p['q5'] = p['q4']+p['b']
        p['q6'] = p['q5']+p['t']

        p['qR'] = ['dR','uR','sR','cR']
        p['qL'] = ['dL','uL','sL','cL']

        p['gq4'] = p['q4']+p['g']
        p['gq5'] = p['q5']+p['g']
        p['gq6'] = p['q6']+p['g']


        # Human readables
        p['gq'] = p['gq6']
        p['k'] = p['gq6']  # k = kromo: strongly interacting sparticles
        
        p['N'] = ['N1','N2','N3','N4']
        p['C'] = ['C1','C2']
        p['X'] = p['N']+p['C']
        
        p['e'] = ['eL','eR']
        p['m'] = ['mL','mR']
        p['T'] = ['T1','T2']
        p['l'] = p['e']+p['m']

        p['v'] = ['ve','vm','vT']   # could be relevant to split out vT
        p['vl'] = ['ve','vm']
        p['vT'] = ['vT']
        
        p['l3'] = p['l']+p['T']  #NB
        p['L'] = p['l3']+p['v']

        p['LX'] = p['X']+p['L']


        # ---------------
        self.comb = {}
        c = self.comb
        c['gl gl'] = ['gl gl']
        c['gg'] = c['gl gl']
        
        c['gl q4'] = self.combine1with2s('gl','q4') 
        c['gl qL'] = self.combine1with2s('gl','qL')
        c['gl qR'] = self.combine1with2s('gl','qR')
        c['gq'] = c['gl q4']
        c['gqL'] = c['gl qL']
        c['gqR'] = c['gl qR']


        c['q4 q4'] = []
        for key1 in self.part['q4']:
            for key2 in self.part['q4']:
                if self.qOrderOk(key1,key2): c['q4 q4'].append(key1+" "+key2)

        c['qL qL'] = []
        for key1 in self.part['qL']:
            for key2 in self.part['qL']:
                if self.qOrderOk(key1,key2): c['qL qL'].append(key1+" "+key2)  # not tested

        c['qR qR'] = []
        for key1 in self.part['qR']:
            for key2 in self.part['qR']:
                if self.qOrderOk(key1,key2): c['qR qR'].append(key1+" "+key2)  # not tested

        c['qL qR'] = []
        for key1 in self.part['qL']:
            for key2 in self.part['qR']:
                if self.qOrderOk(key1,key2): c['qL qR'].append(key1+" "+key2)  # not tested

        c['qLqL'] = c['qL qL']
        c['qRqR'] = c['qR qR']
        c['qLqR'] = c['qL qR']

        c['q5 q5'] = c['q4 q4']+['b1 b1','b2 b2']  #,'b1 b2'  #Hm
        c['q6 q6'] = c['q5 q5']+['t1 t1','t2 t2']  #,'t1 t2'  #Hm
        c['q4q4'] = c['q4 q4']
        c['q5q5'] = c['q5 q5']
        c['q6q6'] = c['q6 q6']
        c['qq'] = c['q6 q6']  # NB, using q as q6

        c['kk'] = c['gg']+c['gq']+c['qq']
        

        c['gl N'] = ['gl N1','gl N2','gl N3','gl N4']
        c['gl C'] = ['gl C1','gl C2']
        c['gl X'] = c['gl N']+c['gl C'] 
        c['gN'] = c['gl N']
        c['gC'] = c['gl C']
        c['gX'] = c['gl X']

        c['q4 N'] = self.combine1swith2s('q4','N')
        c['q4 C'] = self.combine1swith2s('q4','C')
        c['q4 X'] = c['q4 N']+c['q4 C']
        c['qN'] = c['q4 N']
        c['qC'] = c['q4 C']
        c['qX'] = c['q4 X']

        c['kX'] = c['gX']+c['qX']

        c['N N'] = ['N1 N1','N1 N2','N1 N3','N1 N4', 'N2 N2','N2 N3','N2 N4', 'N3 N3','N3 N4', 'N4 N4']
        c['C C'] = ['C1 C1','C1 C2','C2 C2']
        c['N C'] = ['N1 C1','N1 C2','N2 C1','N2 C2','N3 C1','N3 C2','N4 C1','N4 C2']
        c['X X'] = c['N N']+c['C C']+c['N C']
        c['N1 N'] = ['N1 N1','N1 N2','N1 N3','N1 N4']
        c['N1 C'] = ['N1 C1','N1 C2']
        c['N1 X'] = c['N1 N'] + c['N1 C']
        c['NN'] = c['N N']
        c['CC'] = c['C C']
        c['NC'] = c['N C']
        c['XX'] = c['X X']

        c['e e'] = ['eL eL','eR eR']
        c['m m'] = ['mL mL','mR mR']
        c['lL lL'] = ['eL eL','mL mL']
        c['lR lR'] = ['eR eR','mR mR']
        c['l l'] = c['e e']+c['m m']
        c['T T'] = ['T1 T1','T1 T2','T2 T2']
        c['l3 l3'] = c['l l']+c['T T']

        c['l v'] = ['eL ve','mL vm']
        c['T v'] = ['T1 vT','T2 vT']
        c['l3 v'] = c['l v']+c['T v']

        c['vl vl'] = ['ve ve','vm vm']
        c['vT vT'] = ['vT vT']
        c['v v'] = c['vl vl']+c['vT vT']
        
        c['L L'] = c['l3 l3']+c['l3 v']+c['v v']

        # the ones below may not be needed...
        c['ee'] = c['e e']
        c['mm'] = c['m m']
        c['lLlL'] = c['lL lL']
        c['lRlR'] = c['lR lR']
        c['ll'] = c['l l']
        c['TT'] = c['T T']
        c['l3l3'] = c['l3 l3']
        c['lv'] = c['l v']
        c['Tv'] = c['T v']
        c['l3v'] = c['l3 v']
        c['vlvl'] = c['vl vl']
        c['vTvT'] = c['vT vT']
        c['vv'] = c['v v']
        c['LL'] = c['L L']

        c['XL'] = c['XX']+c['LL']  # rather non-standard notion, but ... 
        

    def qOrderOk(self,p1,p2):
        LROrder = {'L':1, 'R':2}
        fOrder = {'u':1, 'd':2, 's':3, 'c': 4}
        f1 = p1[0]; f2 = p2[0]
        LR1 = p1[1]; LR2 = p2[1]
        if (LROrder[LR1] < LROrder[LR2]) or (LROrder[LR1] == LROrder[LR2] and fOrder[f1] <= fOrder[f2]): return 1
        return 0

        
    def combine1swith2s(self,p1,p2):
        res = []
        for key1 in self.part[p1]:
            #for key2 in self.part[p2]: res[p1+"_"+p2].append(key1+' '+key2)
            for key2 in self.part[p2]: res.append(key1+' '+key2)
        return res

    def combine1with2s(self,p1,p2):
        res = []
        for key2 in self.part[p2]:
            #print '%s  %s  %s' %(p1,p2,str(key2))
            #res[p1+"_"+p2].append(p1+" "+key2)
            res.append(p1+" "+key2)
        return res
        
        
        

# =============================================================================
class OX:
    def __init__(self,xT,x,xInput, ID=''):
        self.x = x
        self.xT = xT
        self.xInput = xInput
        self.ID = ID

        self.addTot()
        self.tot = self.x['tot']
        self.xOT = []  #xsec-ordered
        self.xO = {}   #xsec-ordered
        self.makeOrdered()
        self.CombineXsecs()  #produces self.xComb
        

    def N(self):
        return len(self.xO)

    def Sum(self,List,vb=0):
        #if len(self.x) <= 1: return 0.  # NB: return 0 when non-existent   # Hacky? Fragile? (inserted and removed 2012-06-09)
        Sum = 0.
        for proc in List:
            zz = 0.
            if proc in self.xComb: Sum += self.xComb[proc]
            elif proc in self.x: Sum += self.x[proc]
            else:
                if vb: print("OX:Sum: (%s)  Non-existent proc: %s  %s" %(self.ID, proc, list(self.x.keys())))   # this is no problem
        return Sum
            

    def getx(self,key,verb=0): 
        w = key.split()
        if len(w)!=2:
            print("OX:getx: nonallowed %s" %key)
            return 0.
        p1 = w[0]
        p2 = w[1]
        p1p2 = w[0]+" "+w[1]
        p2p1 = w[1]+" "+w[0]
        if p1p2 in list(self.x.keys()): return self.x[p1p2]
        if p2p1 in list(self.x.keys()): return self.x[p2p1]
        if verb: print("Illegal xsec key: %s" %key)
        return 0.
        
    def makeOrdered(self):
        zxO = []
        for key in self.xT:
            zxO.append('%12.7f  %s' %(self.x[key], key))
        zxO.sort()
        for thiszxO in zxO:
            w = thiszxO.split()
            if len(w)==3:
                key = w[1]+" "+w[2]   # all standard ones
            else:
                key = w[1] # tot
                
            self.xOT.append(key)
            self.xO[key] = float(w[0])
        

    def getTot(self):
        xtot = 0.
        for k in self.xT:
            if k!='tot': xtot += self.x[k]
        return xtot
    

    def addTot(self):
        xtot = 0.
        for k in self.xT:
            if k!='tot': xtot += self.x[k]
        if 'tot' not in self.xT:
            self.xT.append('tot')
            self.x['tot'] = xtot
        if 'tot' in self.xT and abs(self.x['tot']-xtot) > 0.001*xtot :
            print("Variation in 'tot': original=%12.7f  recalculated=%12.7f  (%s)" %(self.x['tot'], xtot, self.ID))


    def Show(self):
        for k in self.xT:
            print('%6s : %12.7f' %(k, self.x[k]))

    
    def CombineXsecs(self,vb=0):
        oxcat = OXcategories()
        xsecComb = {}
        for comb in list(oxcat.comb.keys()):
            n = 0
            xsecComb[comb] = 0.
            for subp in oxcat.comb[comb]:
                #print comb
                #print subp
                #print self.x[subp] #crashes on vT vT
                if subp in list(self.x.keys()):
                    xsecComb[comb] += self.x[subp]
                else:
                    xsecComb[comb] += 0.  # this allows for using truncated ox-files (default)
                    
                n += 1
            if vb: print("%6s  %2i  %8.3f" %(comb,n,xsecComb[comb]))
        
        del oxcat
        self.xComb = xsecComb
        

# =============================================================================
def getOX(fn, VB=0):
    #print "<getOsloXsec>"
    xT = []
    x  = {}
    xInput = {}
    if not os.path.exists(fn):
        if VB: print('WARNING  libISAWIG::getOX  getOsloXsec: nonexistent file: %s' %(fn))
        return (xT, x, xInput)
    xT.append('tot')  #tot
    x['tot'] = 0.  #tot
    f=open(fn); l=f.readlines(); f.close()
    for iL in range(len(l)):
        if l[iL].startswith("#"):
            # Look for relevant extra info
            # (some are probably too prospino-specific)
            line = l[iL][1:].strip()
            ww = line.split()
            nw = len(ww)
            if line.startswith('Collider') and nw>=3: xInput['collider'] = line.split()[2]
            if line.startswith('COM Energy') and nw>=4: xInput['comenergy'] = line.split()[3]
            if line.startswith('Calculation') and nw>=3: xInput['calculation'] = line.split()[2]
            if line.startswith('Variable') and nw>=3: xInput['variable'] = line.split()[2]
            if line.startswith('squark mass') and nw>=3: xInput['squarkmass'] = line[line.index(':')+1:].strip()
            if line.startswith('scale') and nw>=3: xInput['squarkmass'] = line[line.index(':')+1:].strip()
            if line.startswith('Units') and nw>=3: xInput['unit'] = line.split()[2]
        else:
            w=l[iL].split()
            if len(w)<3:
                print('WARNING: getOsloXsec: skipping illegal line (below):\n%s' %(l[iL]))
                continue
            key = w[0]+' '+w[1]
            xT.append(key)
            x[key] = float(w[2])
            x['tot'] += x[key]   #tot

    # xT     : list of xsec names (only required to keep a natural order of the xsecs)
    # x      : dict of the xsecs with keys as in xT
    # xInput : some extra (calculation) info availble in the ox file, like energy,
    
    return (xT, x, xInput)

# =============================================================================
        

#SMp_keys = ['h','Z','W','t','b','q','l','T','v']
SMp_keys =  ['l','T','q','b','h','Z','W','t']
# SMp_keys2 = ['QQ','TT','vT','lv','ll']  #Sep20: fixed vT->Tv, added 'vv'
SMp_keys2 = ['QQ','vv','TT','Tv','lv','ll']
# =============================================================================
class TwoLeg:
    def __init__(self,SU1,SU2,SM1,SM2,id1,id2,BR,xsec, SMp={},):
        self.SU1 = SU1
        self.SU2 = SU2
        self.SM1 = SM1
        self.SM2 = SM2
        self.id1 = id1
        self.id2 = id2
        self.BR = BR
        self.xsec = xsec
        self.SMp = {}
        if SMp: self.SMp = SMp.copy()   # is this dict treated correctly (is a new instance created?)
        else: self.FillSMp()
                
    def Copy(self):
        # (this is correct copying, although SU1 = list(g.SU1) is easier/more intuitive
        SU1 = []; SU1.extend(self.SU1)
        SU2 = []; SU2.extend(self.SU2)
        SM1 = []; SM1.extend(self.SM1)
        SM2 = []; SM2.extend(self.SM2)
        return TwoLeg(SU1, SU2, SM1, SM2, self.id1, self.id2, self.BR, self.xsec,  self.SMp.copy())

    def FillSMp(self):
        #print self.SM1+self.SM2
        self.SMp = {}
        for ip in self.SM1+self.SM2:
            if len(ip)>1: print('FillSMp: long ip: %s' %(str(ip)))  #could make more general. Expect to never happen. 
            ipp = ip[0] 
            if ipp == '-': continue
            #print ipp
            #print self.SMp.keys()
            if ipp not in list(self.SMp.keys()): self.SMp[ipp] = 0
            self.SMp[ipp] += 1
                
    def txtSMp(self,mode=[],keys=[],showzero=0):
        txt = ""
        headtxt = ""
        space = " " #"_"
        zero  = "  " #"0"

        SMp_keys12 = SMp_keys + SMp_keys2  # uses globals
        if keys: SMp_keys12 = keys  # allow to pipe in
        for ikey in range(len(SMp_keys12)):
            key = SMp_keys12[ikey]
            if ikey!=0:
                txt += space
                headtxt += space

            dtxt = zero
            if key in list(self.SMp.keys()):
                if self.SMp[key]!=0 or showzero:  #(also not print 0)
                    dtxt = '%2i' %(self.SMp[key])
            txt += dtxt
            headtxt += '%2s' %(key)

        if 'head' in mode: return headtxt
        else: return txt


    def Show(self, ret=False, SM=['+SM'], head=False):
        # SM: [], ['SM'], ['SMp'], ['+SM'], ['+SMp']

        txt = '%6.4f  %-30s  %-30s' %(self.BR, self.strS(self.SU1), self.strS(self.SU2))

        if 'SM' in SM:
            txt = '%6.4f  %-30s  %-30s' %(self.BR, self.strS(self.SM1), self.strS(self.SM2))


        if '+SM' in SM:
            txt = '%6.4f  %-30s  %-30s  %-30s  %-30s' %(self.BR, self.strS(self.SU1), self.strS(self.SU2), self.strS(self.SM1), self.strS(self.SM2))


        if 'SMp' in SM:
            if head:
                txt =   '%6s  %-30s' %('', self.txtSMp(mode=['head']))
            else: 
                txt = '%6.4f  %-30s' %(self.BR, self.txtSMp())
        

        if '+SMp' in SM:
            if head:
                txt =   '%6s  %-30s  %-30s  %-30s' %('', '', '', self.txtSMp(mode=['head']))
            else: 
                txt = '%6.4f  %-30s  %-30s  %-30s' %(self.BR, self.strS(self.SU1), self.strS(self.SU2), self.txtSMp())



        #txt = "%6.4f  %-40s  %-40s" %(self.BR, str(self.SU1), str(self.SU2))
        #if SM:
        #    txt = "%6.4f  %-40s  %-40s  %-20s  %-20s" %(self.BR, str(self.SU1), str(self.SU2), str(self.SM1), str(self.SM2))
        if ret:
            return txt
        else:
            print(txt)

    def strS(self, SS):
        txt = ""
        for iS in range(len(SS)):
            if iS != 0: txt += "__"
            if type(SS[iS]) is str: txt += SS[iS]   #SU
            if type(SS[iS]) is list: txt += SS[iS][0] #SM
        return txt

    def Legvalue(self):
        #return (len(self.SU1), len(self.SU2), str(self.SU1), str(self.SU2))
        return (len(self.SU1), str(self.SU1), len(self.SU2), str(self.SU2))

    def SMpvalue(self):
        return self.txtSMp(self.SMp)  # .. is it correct to insert self.SMp ? 

    def ShortSmallestLeft(self):
        if (len(self.SU1), str(self.SU1)) > (len(self.SU2), str(self.SU2)):         
            x = self.SU1; self.SU1 = self.SU2; self.SU2 = x
            x = self.SM1; self.SM1 = self.SM2; self.SM2 = x
            x = self.id1; self.id1 = self.id2; self.id2 = x
            


    def identicalTo(self, tl2, check=[]): 
        sameSU = 1*(self.SU1==tl2.SU1 and self.SU2==tl2.SU2) | 2*(self.SU1==tl2.SU2 and self.SU2==tl2.SU1)
        sameSM = 1*(self.SM1==tl2.SM1 and self.SM2==tl2.SM2) | 2*(self.SM1==tl2.SM2 and self.SM2==tl2.SM1)
        #sameSMp= (self.SMp==tl2.SMp)
        sameSMp= (self.SMpvalue()==tl2.SMpvalue()) #only checks relevant particles (those in SMp_keys)
        doboth = ('LRboth' in check)
        doone  = ('LRboth' not in check)

        #print check
        fullOK = 1
        
        if 'SU' in check and 'SM' in check:
            if doone  and not ( sameSU&1 and sameSM&1 ): fullOK *= 0
            if doboth and not (sameSU&sameSM ): fullOK *= 0

        elif 'SU' in check and 'SMp' in check:
            if doone and not  ( sameSU&1 and sameSMp ): fullOK *= 0
            if doboth and not ( sameSU and sameSMp ): fullOK *= 0

        elif 'SU' in check:
            #print 'here'
            if doone  and not ( sameSU&1 ): fullOK *= 0
            if doboth and not ( sameSU ): fullOK *= 0

        elif 'SM' in check:
            if doone  and not ( sameSM&1 ): fullOK *= 0
            if doboth and not ( sameSM ): fullOK *= 0

        elif 'SMp' in check:
            if not ( sameSMp ): fullOK *= 0

        else:
            print('identicalTo: not supposed to enter here: check = %s' %(str(check)))
            print('Exiting')
            sys.exit()

        return fullOK


    def JoinNs(self):
        SUs = [self.SU1, self.SU2]
        for SU in SUs:
            for iS in range(len(SU)):
                if isN(SU[iS]): SU[iS] = 'N_'

    def JoinCs(self):
        SUs = [self.SU1, self.SU2]
        for SU in SUs:
            for iS in range(len(SU)):
                if isC(SU[iS]): SU[iS] = 'C_'

    def JoinNCs(self):
        SUs = [self.SU1, self.SU2]
        for SU in SUs:
            for iS in range(len(SU)):
                if isNC(SU[iS]): SU[iS] = 'X_'

    def JoinlL(self):
        SUs = [self.SU1, self.SU2]
        for SU in SUs:
            for iS in range(len(SU)):
                if SU[iS] in ['eL','mL']: SU[iS] = 'lL'

    def JoinlR(self):
        SUs = [self.SU1, self.SU2]
        for SU in SUs:
            for iS in range(len(SU)):
                if SU[iS] in ['eR','mR']: SU[iS] = 'lR'

    def Joinsl(self):
        SUs = [self.SU1, self.SU2]
        for SU in SUs:
            for iS in range(len(SU)):
                if SU[iS] in ['eL','mL','eR','mR','lR','lL']: SU[iS] = 'sl'

    def Joinsl(self):
        SUs = [self.SU1, self.SU2]
        for SU in SUs:
            for iS in range(len(SU)):
                if SU[iS] in ['eL','mL','eR','mR','lR','lL']: SU[iS] = 'sl'

    def JoinQs(self,pJoin='SQ', start=0, end=1, keep=[]):  # based on Cascade.JoinQs
        # setq = ['d','u','s','c']
        # setQ = ['d','u','s','c','b']
        # setj = ['d','u','s','c','b','T']

        if start<0: start = 0
        if end<start: end = self.N()
        
        # 1) Rewrite,  2) Join
        aJoin = {}
        aJoin["qL"] = ["dL","uL","sL","cL"]
        aJoin["qR"] = ["dR","uR","sR","cR"]
        # aJoin["sq"] = ["qL","qR"]
        aJoin["sq"] = ["qL","qR","dL","uL","sL","cL","dR","uR","sR","cR"]
        aJoin["sb"] = ["b1","b2"]
        #aJoin["SQ"] = ["qL","qR","b1","b2","sq"]
        #aJoin["SQ"] = ["sq", "qL","qR","b1","b2", "dL","uL","sL","cL","dR","uR","sR","cR"] #preApr23
        aJoin["SQ"] = ["sq", "qL","qR","b1","b2", "dL","uL","sL","cL","dR","uR","sR","cR","t1","t2"]
        # NB enlarged compared to original (can then do 'SQ' in one step,
        #    otherwise needed to do e.g. 'qL','qR','SQ'
        #    Should be no difference

        SUs = [self.SU1, self.SU2]
        SMs = [self.SM1, self.SM2]

        for iS in range(len(SUs)):
            SU = SUs[iS]
            SM = SMs[iS]
        
            nSU = len(SU)
            for iS in range(min(start,nSU), min(end,nSU)):
            #for iS in range(min(start,nSU), max(end,nSU)):
                if SU[iS] in aJoin[pJoin]:
                    SU[iS] = pJoin

                    # adhoc for sb->top  ... took out Apr20 for 2leg
                    if ('sbtop' in keep) and pJoin in ['SQ'] and (self.SM[iS][0]=='t' or self.SM[iS+1][0]=='t'):
                        SU[iS] = 'sb'


                    # THEN DO SM:

                    if SM[iS][0] in setq: SM[iS]=['q']
                    if pJoin in ['sb','SQ'] and SM[iS][0] in setQ: SM[iS] = ['q']


                    # The lines below (either/both) do some things ...
                    if SM[iS+1][0] in setq:
                        #print '||| %s   %s  %i' %(str(SM), SM[iS+1][0], iS)
                        SM[iS+1]=['q']
                        #print '    %s   %s  %i' %(str(SM), SM[iS+1][0], iS)
                       
                    if pJoin in ['sb','SQ'] and SM[iS+1][0] in setQ: SM[iS+1] = ['q']


                    '''

                    # adhoc for sb->top  ... took out Apr20 for 2leg
                    if ('sbtop' in keep) and pJoin in ['SQ'] and (SM[iS][0]=='t' or SM[iS+1][0]=='t'):
                        SU[iS] = 'tb'
                    if ('sttop' in keep) and pJoin in ['SQ'] and (SM[iS][0]=='t' or SM[iS+1][0]=='t'):
                        SU[iS] = 'tb'
                    '''


# =============================================================================
class TwoLegs:
    def __init__(self,scenID): 
        self.list = []
        self.BRcut = 0.
        self.BRtot = 0.
        self.BRtotall = 0.
        self.BRtotreq = 0.
        self.xsec = 0.
        self.scenID = scenID
        self.h = {}


    def FillH(self,hmax = 0.5, hmin = 1e-3):
        #col = {-1:kBlack, 0:kBlue, 1:kRed, 2:kGreen, 3:kMagenta, 4:kPink, 5:kYellow, 6:kBlack, 7:kBlack, 8:kBlack}
        col = {-1:1, 0:600, 1:632, 2:416, 3:616, 4:900, 5:400, 6:1, 7:1, 8:1}  #Sep14
        wid = {-1:3, 0:2, 1:2, 2:2, 3:2, 4:2, 5:2, 6:2, 7:2, 8:2}
        sty = {-1:2, 0:1, 1:1, 2:1, 3:1, 4:1, 5:1, 6:1, 7:1, 8:1}
        # Fills histograms according to content of SMp
        # h[0]: 0-lepton
        from ROOT import TH1F  # hm, might make more general
        L = -1

        # First create all histos (until/including l=8)
        for l in range(-1,9):
            zh = TH1F('l'+str(l), 'l'+str(l), 11,-0.5, 10.5)
            zh.SetLineColor(col[l])
            zh.SetLineWidth(wid[l])
            zh.SetLineStyle(sty[l])
            self.h[l] = zh
            if hmax!=0: self.h[l].SetMaximum(hmax)
            if hmin!=0: self.h[l].SetMinimum(hmin)
            
        # Then loop over all twolegs and fill
        for tl in self.list:
            l = 0
            if 'l' in list(tl.SMp.keys()): l = tl.SMp['l']
            q = 0
            if 'q' in list(tl.SMp.keys()): q = tl.SMp['q']

            if l>6:
                print('FillH: too many leptons!! n(l) = ',str(l))
                continue

            #if l not in self.h.keys():
            #    zh = TH1F('l'+str(l), 'l'+str(l), 11,-0.5, 10.5)
            #    zh.SetLineColor(col[l])
            #    zh.SetLineWidth(2)
            #    self.h[l] = zh
            #if L not in self.h.keys():
            #    zh = TH1F('l'+str(L), 'l'+str(L), 11,-0.5, 10.5)
            #    zh.SetLineColor(col[L])
            #    zh.SetLineWidth(3)
            #    zh.SetLineStyle(2)
            #    self.h[L] = zh
            
            self.h[l].Fill(q, tl.BR)
            self.h[L].Fill(q, tl.BR)


    #def ShowH(self,cc,tleg,lA=0,lB=-1):
    def ShowH(self, which=[-1,0,1,2], fn=[], silent=False):
        # TO DO FROM COMMAND LINE: 
        from ROOT import TCanvas, TLegend, gStyle, gROOT
        
        if silent: gROOT.SetBatch(True)  # doesn't show TCanvas 
        
        #g.twolegs.FillH()
        #g.twolegs.ShowH(cc,tleg,0,2)
        #
        #tleg.Clear()
        gStyle.SetPadColor(0)
        gStyle.SetStatColor(0)
        gStyle.SetCanvasColor(0)
        gStyle.SetFillColor(0)
        gStyle.SetOptStat(0)
        gStyle.SetPadBorderSize(1)
        gStyle.SetStatBorderSize(1)
        gStyle.SetCanvasBorderSize(1)
        gStyle.SetLegendBorderSize(1)

        cc = TCanvas('cc','',0,0,600,400)
        tleg = TLegend(0.7,0.7,0.9,0.9)
        cc.SetLogy()
        cc.cd()

        #tleg = TLegend(0.7,0.7,0.9,0.9)
        L = -1
        for l in which: 
            px = 'hist'
            hh = self.h[l]
            if l!=-1:
                px += 'same'
                tleg.AddEntry(hh,'%i lepton' %l)
            else:
                hh.SetTitle('')
                tax = hh.GetXaxis()
                tax.SetTitle('Number of jets')
                tax.CenterTitle()
                tay = hh.GetYaxis()
                tay.SetTitle('Branching ratio')
                tay.CenterTitle()
                
                tleg.AddEntry(hh,'Sum leptons')
                
            hh.Draw(px)
        tleg.Draw()

        # Save canvas to file
        for fname in fn:
            cc.SaveAs(fname)
            
        return cc, tleg
        
    def ExpandSMp(self,mode=['h']): 
        # mode h: h->bb (1) [assumption]
        # mode W: h->W->lv,Tv,QQ (3)
        # mode Z: h->Z->ll,TT,vv,QQ
        # mode t: bW (so might want to call before mode W)
        # mode b: b->q
        # mode T: T->l,Q,..
        ch = {}
        #ch['t'] = [(PDG.Wtolv,'b','l','v'), (PDG.WtoTv,'b','T','v'), (PDG.Wtoqq,'b','q','q')]

        ch['h'] =  [[1.,'b','b']]  #NB: approx-hack
        ch['H'] =  [[1.,'b','b']]  #NB: approx-hack
        ch['A'] =  [[1.,'b','b']]  #NB: approx-hack
        ch['H+'] = [[1.,'t','b']]  #NB: approx-hack

        ch['t'] = [[1.,'b','W']]
        ch['W'] = [[PDG.Wtolv,'l','v'], [PDG.WtoTv,'T','v'], [PDG.Wtoqq,'q','q']]

        ch['l'] = [[1.]]  # 2012-06-23
        
        # ###### NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB
        # ###### This is TwoLegs::ExpandSMp() 
        # ###### Most particles are probably already expanded in the 1-leg step, i.e. in CascadeStore::ExpandSMp()
        # ###### You need to implement changes there (in this file)
        # ###### NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB
        
        ch['Z'] = [[PDG.Ztoll,'l','l'], [PDG.ZtoTT,'T','T'], [PDG.Ztovv,'v','v'], [PDG.ZtoQQ,'q','q']]
        ch['QQ'] = [[1.,'q','q']] #added Apr26
        #ch['Q'] = [[1.,'q']]      #added Apr26  In use?
        ch['TT'] = [[1.,'T','T']] #added Apr26
        ch['ll'] = [[1.,'l','l']] #added Apr26
        ch['vv'] = [[1.,'v','v']] #added Sep20

        #Sep19: keep v
        ch['lv'] = [[1.,'l','v']] #added Apr26
        ch['Tv'] = [[1.,'T','v']] #added Apr26
        # Sep16
        ch['blv'] = [[1.,'b','l','v']] #added Apr26
        ch['bTv'] = [[1.,'b','T','v']] #added Apr26
        ch['bqq'] = [[1.,'b','q','q']] #added Apr26

        
        ch['b'] = [[1.,'q']]
        ch['T'] = [[PDG.Ttolvv,'l','v','v'], [PDG.Ttovh,'q','v']]

        #ch['q'] = [[1.,'v']]  # Hack
        ch['q'] = [[1.]]


        expanded = False
        KeepOn = True
        while KeepOn: # this outermost loop is now probably not needed

            KeepOn = False
            iTL = 0
            # Iterate through list. Will continuously expand, so use dynamical self.N()
            while iTL < self.N():
                ##print 'while %i < %i' %(iTL, self.N())

                #pOlds = SMp_keys2 + ['t','h','W','Z','b','T']  #Apr26  "all" possibilities to transform
                #pOlds = SMp_keys2 + ['t','h','W','Z','b','T','q']  #Sep19  "all" possibilities to transform
                pOlds = SMp_keys2 + ['t','h','W','Z','b','T','q','l']  # 2012-06-23

                for pOld in pOlds: 
                    if pOld in mode:
                        
                        if pOld in list(self.list[iTL].SMp.keys()): 
                            if self.list[iTL].SMp[pOld] == 0: continue #adhoc

                            expanded = True
                            TL = self.list[iTL]

                            # Remove list item (but still have it in TL)
                            self.list.pop(iTL)
                            # Change TL
                            ##print 'TL: %s  %i -> ' %(pOld,TL.SMp[pOld]),
                            ###print "OLD:   ", TL.SMp, '   ', pOld
                            TL.SMp[pOld] -= 1
                            ##print '%s' %(TL.SMp[pOld])
                        
                            # Loop over transformations: create new objects, adjust, insert
                            for dec in ch[pOld]:
                                # new object
                                ##print 'TL   : %s %i' %(pOld,TL.SMp[pOld])
                                TLnew = TL.Copy()
                                #print "%3s" %(pOld)
                                ##print 'TL   : %s %i' %(pOld,TL.SMp[pOld])
                                ##print 'TLnew: %s %i' %(pOld,TLnew.SMp[pOld])
                                # adjust BR
                                TLnew.BR *= dec[0]
                                TLnew.xsec *= dec[0]  #Sep27
                                # add the new particles (of decay/transformation)
                                for ipnew in range(1,len(dec)):
                                    pnew = dec[ipnew]
                                    if pnew in list(TLnew.SMp.keys()): TLnew.SMp[pnew] += 1
                                    else: TLnew.SMp[pnew] = 1
                                # insert in list (a possibility is to append ... might be faster)
                                ###print "NEW:   ", TLnew.SMp
                                self.list.insert(iTL,TLnew)
                            # continue without increasing iTL, i.e. redo the same entry (for possible additional changes)
                            ###print
                            del TL
                            continue

                        # (now do next pOld happening to be in mode)

                    # (now try next pOld)

                #
                iTL += 1

        return expanded  # so that I know if an expansion took place


    def Add(self,SU1,SU2,SM1,SM2,i1,i2,BR,xsec, SMp={}):  #Sep20 added SMp
        #anew = TwoLeg(SU1,SU2,SM1,SM2,i1,i2,BR, SMp)
        anew = TwoLeg(list(SU1),list(SU2),list(SM1),list(SM2),i1,i2,BR,xsec, dict(SMp)) #Sep20: fix
        # with the preSep20 solution, the lists in cascStore were changed when twolegs were changed/refined
        self.list.append(anew)
        self.BRtot += BR
        self.xsec += xsec

    def N(self):
        return len(self.list)

    def Status(self, ret=False, showID=True):
        line = 'BRtot = %6.4f [of %6.4f] ,  N(2legs) = %4i  ,  BRcut = %7.1e'  %(self.BRtot, self.BRtotreq, self.N(), self.BRcut)
        if showID: line += '    %s' %self.scenID
        if ret: return line
        else: print(line)

    def ShowSMp(self, i2=-1, i1=0, SMpHead2=1, fn="", showzero=1, keys = ['l','q','T','b','y']):
        # nearly same as Show, but redid to keep compatibility with any old usage. Now show xsec too.
        if i2<0 or i2>self.N(): i2 = self.N()
        if i1<0 or i1>i2: i1 = 0

        BRacc = 0.
        XSECacc = 0.
        l = []
        

        # make list
        for iC in range(i1,i2):
            tl  = self.list[iC]
            if (SMpHead2==1 and iC==i1 or SMpHead2==2) and not fmod(i1-iC, 20): 
                line = '%8s %8s | %8s  %8s  %s' %('XSECaccu','XSECsub','BRaccu','BRsub',tl.txtSMp(mode=['head'],keys=keys))  #tl.Show(ret=True, SM=['SMp'], head=True))
                l.append(line)
                l.append('-'*len(line))
            BRacc += tl.BR
            XSECacc += tl.xsec
            line = '%8.3f %8.3f | %8.6f  %8.6f  %s' %(XSECacc, tl.xsec, BRacc, tl.BR, tl.txtSMp(keys=keys, showzero=showzero))

            l.append(line)

        # print/save list
        if not fn.strip():
            print()
            for line in l: print(line)
        else:
            f=open(fn.strip(),'w')
            for line in l: f.write(line+'\n')
            f.close()
        


    def Show(self, i2=-1, i1=0, fn="", status=True, txt="", opt=['showBRsum'], SM=['+SM'], SMpHead2=1):
        if i2<0 or i2>self.N(): i2 = self.N()
        if i1<0 or i1>i2: i1 = 0
        l=[]
        BRsofar = 0.

        # ADD STATUS?
        if status:
            l.append(self.Status(ret=True)+"   "+txt)

        # THE 1LEGS
        for iC in range(i1,i2):

            # header if SMp-option
            if (SMpHead2==1 and iC==i1 or SMpHead2==2) and ('SMp' in SM or '+SMp' in SM) and not fmod(i1-iC, 20): 
                line = '%13s%s' %('',self.list[iC].Show(ret=True, SM=SM, head=True))
                l.append(line)
                
            # list    
            line = self.list[iC].Show(ret=True,SM=SM)
            if 'showBRsum' in opt:
                BRsofar += self.list[iC].BR
                line = 'sum=%6.4f   %s' %(BRsofar, line)
            l.append(line)

        # Output to file (if fn is set) or screen:
        if not fn.strip():
            for line in l: print(line)
        else:
            f=open(fn.strip(),'w')
            for line in l: f.write(line+'\n')
            f.close()
        
            
    def sumBR(self,i1=0,i2=-1): 
        if i2<0 or i2>self.N(): i2 = self.N()
        if i1<0 or i1>i2: i1 = 0
        self.BRtot = 0.
        for iC in range(i1,i2):
            self.BRtot += self.list[iC].BR
        
    def ShortSmallestLeft(self):
        for twoleg in self.list: twoleg.ShortSmallestLeft()

    def Duplicates(self,mode='join',check=['SU']): #, LRcheckboth=False):
        # check: 'SU':     check also SU
        # check: 'SMp':    check also SMp (intermediate)
        # check: 'SM':     check also SM1 and SM2 (strictest)
        # check: 'LRboth': check 1==1 and 2==2 as well as 1==2 and 2==1 (not needed if ShortSmallestLeft() has been executed
        iC1 = 0
        while iC1 < len(self.list)-1: 
        #for iC1 in range(len(self.list)-1):
            duplicates = []
            for iC2 in range(iC1+1,len(self.list)):
                if self.list[iC1].identicalTo(self.list[iC2],check): #LRcheckboth):
                    duplicates.append(iC2)

            # Reverse list to keep list indices correct when popping has started
            duplicates.reverse()
            # Handle duplicates (join) to iC1
            for iC2 in duplicates:
                if mode=='join':
                    self.list[iC1].BR += self.list[iC2].BR
                    self.list[iC1].xsec += self.list[iC2].xsec
                self.list.pop(iC2)

            iC1 += 1

    def interchange(self,iA,iB):
        L = self.list
        SU1 = L[iA].SU1 ; L[iA].SU1 = L[iB].SU1 ; L[iB].SU1 = SU1
        SU2 = L[iA].SU2 ; L[iA].SU2 = L[iB].SU2 ; L[iB].SU2 = SU2
        SM1 = L[iA].SM1 ; L[iA].SM1 = L[iB].SM1 ; L[iB].SM1 = SM1
        SM2 = L[iA].SM2 ; L[iA].SM2 = L[iB].SM2 ; L[iB].SM2 = SM2
        SMp = L[iA].SMp ; L[iA].SMp = L[iB].SMp ; L[iB].SMp = SMp
        id1 = L[iA].id1 ; L[iA].id1 = L[iB].id1 ; L[iB].id1 = id1
        id2 = L[iA].id2 ; L[iA].id2 = L[iB].id2 ; L[iB].id2 = id2
        BR  = L[iA].BR  ; L[iA].BR  = L[iB].BR  ; L[iB].BR  = BR          
        xsec= L[iA].xsec; L[iA].xsec= L[iB].xsec; L[iB].xsec= xsec          

    def JoinlR(self):
        for iC in range(len(self.list)): self.list[iC].JoinlR()

    def JoinlL(self):
        for iC in range(len(self.list)): self.list[iC].JoinlR()

    def Joinsl(self):
        for iC in range(len(self.list)): self.list[iC].Joinsl()

    def JoinQs(self,pJoin='SQ', start=0, end=1, keep=[]):
        for iC in range(len(self.list)):
            self.list[iC].JoinQs(pJoin,start)

    def JoinNs(self):
        for iT in range(len(self.list)): self.list[iT].JoinNs()

    def JoinCs(self):
        for iT in range(len(self.list)): self.list[iT].JoinCs()

    def JoinNCs(self):
        for iT in range(len(self.list)): self.list[iT].JoinNCs()

    def FillSMp(self):
        for iT in range(len(self.list)): self.list[iT].FillSMp()

    def orderWithBR(self):
        for iC1 in range(len(self.list)-1):
            for iC2 in range(iC1+1,len(self.list)):
                if self.list[iC1].BR < self.list[iC2].BR:
                    self.interchange(iC1,iC2)

    def orderWithLegvalue(self):
        for iC1 in range(len(self.list)-1):
            for iC2 in range(iC1+1,len(self.list)):
                if self.list[iC1].Legvalue() > self.list[iC2].Legvalue():
                   self.interchange(iC1,iC2)
                
    def orderWithSMp(self):
        for iC1 in range(len(self.list)-1):
            for iC2 in range(iC1+1,len(self.list)):
                if self.list[iC1].SMpvalue() < self.list[iC2].SMpvalue():  # with leptons to the left we select on leptons
                   self.interchange(iC1,iC2)

        
# =============================================================================
class Cascade:
    # def __init__(self,thisSUSY,thisCODE):
    def __init__(self,thisSUSY): #,thisCODE):
        #self.BRtot = 1.
        #self.N = 0
        self.rate = 0 #not yet in use
        self.BR = []
        self.BRcum = []
        self.BRtot = 1.
        self.BRtotBefDup = 1.
        self.SUSY = [] #initial particle
        self.SM = []
        self.SMp = {}
        # self.CODE = []
        
        # self.Add(thisSUSY,thisCODE)
        self.Add(thisSUSY)

        self.dilepX = []
        self.ditauX = []
        self.casctype = {}
        self.casctypeI = {} #not yet in use
        

    # shead = ""
    # top = 0
    # pi = 0
    # nu = 0
    # exlep = 0
    # dilep = 0
    fullSUSY = ""
    fullSM = ""

    isDilep = 0
    isSleptonLL = 0
    isSnuLL = 0
    isZll = 0 
    is3bodyLL = 0 

    isDitau = 0
    isStauTT = 0
    isSnuTT = 0
    isZTT = 0 
    is3bodyTT = 0 

    def N(self): return len(self.SUSY)

    # Sep16
    def nSM(self): return len(self.SM)
    def nSUSY(self): return len(self.SUSY)
    def nSMp(self): return len(self.SMp)


    def CopySep19(self):
        # newCasc = Cascade(self.SUSY[0], self.CODE[0])
        newCasc = Cascade(self.SUSY[0])
        for iP in range(1,self.N()):
            newCasc.SUSY.append(self.SUSY[iP])
        newCasc.SM = list(self.SM)
        newCasc.SMp = dict(self.SMp)
        newCasc.BRtot = self.BRtot
        # "minimal copy", do I need more? 
        
        return newCasc


    def FillSMp(self):
        for ip in self.SM:
            if len(ip)>1: print('FillSMp: long ip: %s' %(str(ip)))  #could make more general. Expect to never happen. 
            ipp = ip[0] 
            if ipp == '-': continue
            #print ipp
            #print self.SMp.keys()
            if ipp not in list(self.SMp.keys()): self.SMp[ipp] = 0
            self.SMp[ipp] += 1
        

    def Add(self,thisSUSY,thisSM=['-'],thisBR=1.):
        self.SUSY.append(thisSUSY)
        # self.CODE.append(thisCODE)
        self.SM.append(thisSM)
        self.BR.append(thisBR)
        #self.N += 1
        N = self.N()
        if N==1:
            self.BRcum.append(thisBR)
        else:
            #print "aaa",self.N-2
            self.BRcum.append(thisBR * self.BRcum[N-2])
            #print "bbb"
        self.BRtot = self.BRcum[N-1]

        # self.StatDilep()
        

    def RedoBR(self): #NB this will fail if executed after JoinDuplicates
        self.BRcum[0] = self.BR[0] #needed?
        for iP in range(1,self.N()):
            self.BRcum[iP] = self.BRcum[iP-1] * self.BR[iP]
        self.BRtot = self.BRcum[self.N()-1]

    def strSUSY(self):
        txt = ""
        for iS in range(self.N()):
            if iS!=0: txt += "__"
            txt += self.SUSY[iS]
        return txt

    def strSM(self):
        txt = ""
        for iS in range(self.nSM()):
            if iS!=0: txt += "__"
            txt += self.SM[iS][0]  # NB: assuming SM-particles are joined (but vector kept)
            if len(self.SM[iS]) >= 2: txt += '+'+self.SM[iS][1]  #Apr17
        return txt

    def strSMp(self,mode=0):
        txt = ''
        keys = list(self.SMp.keys())
        keys.sort()
        for key in keys:
            if mode==1: txt += '%s: ' %(key)
            txt += '%i   ' %(self.SMp[key])
        return txt
            

    def eraseSM(self):
        for iS in range(self.N()):
            self.SM[iS] = 'x'  # maybe need ['x'] NBNB

    def RedoStrSUSYSM(self):
        self.fullSUSY = self.strSUSY()
        self.fullSM = self.strSM()

    # ### Below: classifications

    def sameCasctype(self,testtype={}):
        for var in list(testtype.keys()):
            if testtype[var] != self.casctype[var]: return False
        return True

    def isOfType(self,casctypeA={}): #NB casctypeA is an array, its entities can be arrays
        for key in casctypeA:
            # Put testrequirements into list form if not already
            if type(casctypeA[key]) is list:
                varA = casctypeA[key]
            else:
                varA = [casctypeA[key]]

            # Put self.casctype into list form for vars not already so (for uniformity)
            #print type(self.casctype[key])
            if type(self.casctype[key]) is list:
                #print "yes"
                thistypeA = self.casctype[key]
            else:
                thistypeA = [self.casctype[key]]

            #print self.casctype[key], '  -->  ', thistypeA


            # Test all (e.g. if self has dilepX=['ZLL','SleLL'] then will require both)
            for thistype in thistypeA: 
                #print key, thistypeA, " ----- VS ----- ", varA
                if thistype not in varA: return 0
                
        return 1
            

    def Classify(self,categ=[]):
        SUSY = self.SUSY
        SM = self.SM[0]
        #fullSM = ""
        #for zSM in SM: fullSM+=zSM+"__"
        fullSM = self.strSM()
        # ISshead = ['gl','SQ','sq','sb','st']
        sheadLab = {'gl':'G', 'SQ':'Q','sq':'Q','sb':'Q','st':'Q','SHEAD':'SHEAD'}
        
        if 'shead' in categ:  #assuming for now StrongHEAD only
            shead = ""
            iS = 0
            while SUSY[iS] in list(sheadLab.keys()):
                shead += sheadLab[SUSY[iS]]
                iS += 1
            self.casctype['shead']=shead

        if 'top' in categ:
            self.casctype['top'] = fullSM.count("t")

        if 'pi' in categ:
            self.casctype['pi'] = fullSM.count("PI")

        if 'nu' in categ:
            self.casctype['nu'] = fullSM.count("v")

        if 'dilep' in categ:
            self.casctype['dilep'] = self.getnDilep()
            self.casctype['dilepX'] = self.strDilep

        if 'exlep' in categ:
            self.casctype['exlep'] = fullSM.count("l") - 2*self.isDilep

        if 'ditau' in categ:
            self.casctype['ditau'] = self.getnDitau()
            self.casctype['ditauX'] = self.strDitau

        if 'extau' in categ:
            self.casctype['extau'] = fullSM.count("T") - 2*self.isDitau

        if 'y' in categ:
            self.casctype['y'] = fullSM.count("y")


    # ##################
    
    def getnSleptonLL(self):
        n = 0
        for iP in range(self.N()-1):
            if isSlepton(self.SUSY[iP]) and isLepton(self.SM[iP][0]) and isLepton(self.SM[iP+1][0]):
                n += 1
        self.isSleptonLL = n
        return n

    def getnSnuLL(self):
        n = 0
        for iP in range(self.N()-1):
            if isSnu(self.SUSY[iP]) and isLepton(self.SM[iP][0]) and isLepton(self.SM[iP+1][0]):
                n += 1
        self.isSnuLL = n 
        return n

    def getnZll(self):
        fullSM = self.strSM()
        self.isZll =  fullSM.count("_Z->ll")
        return self.isZll  #fullSM.count("_Z->ll")

    def getn3bodyLL(self):
        fullSM = self.strSM()
        self.is3bodyLL = fullSM.count("_ll")
        return self.is3bodyLL

    def getnDilep(self):
        strDilep = []
        for iLL in range(self.getnSleptonLL()): strDilep.append("sleLL")
        for iLL in range(self.getnSnuLL()): strDilep.append("snuLL")
        for iLL in range(self.getnZll()): strDilep.append("ZLL")
        for iLL in range(self.getn3bodyLL()): strDilep.append("3bLL")
        self.strDilep = strDilep
        self.isDilep = len(strDilep) 
        return len(strDilep)
    # ##################
    
    def getnStauTT(self):
        n = 0
        for iP in range(self.N()-1):
            if isStau(self.SUSY[iP]) and isTau(self.SM[iP][0]) and isTau(self.SM[iP+1][0]):
                n += 1
        self.isStauTT = n
        return n

    def getnSnuTT(self):
        n = 0
        for iP in range(self.N()-1):
            if isSnu(self.SUSY[iP]) and isTau(self.SM[iP][0]) and isTau(self.SM[iP+1][0]):
                n += 1
        self.isSnuTT = n 
        return n

    def getnZTT(self):
        fullSM = self.strSM()
        self.isZTT =  fullSM.count("_Z->TT")
        return self.isZTT 

    def getn3bodyTT(self):
        fullSM = self.strSM()
        self.is3bodyTT = fullSM.count("_TT")
        return self.is3bodyTT

    def getnDitau(self):
        strDitau = []
        for iTT in range(self.getnStauTT()): strDitau.append("staTT")
        for iTT in range(self.getnSnuTT()): strDitau.append("snuTT")
        for iTT in range(self.getnZTT()): strDitau.append("ZTT")
        for iTT in range(self.getn3bodyTT()): strDitau.append("3bTT")
        self.strDitau = strDitau
        self.isDitau = len(strDitau) 
        return len(strDitau)
    
    # ##################
    

    # ### Above: classifications
    def CollapseShead(self,mode=0):
        # Count shead
        nshead = 0
        while self.SUSY[nshead] in ['SHEAD','G','gl','SQ','Q','sq','sb','st']:
            nshead += 1

        if nshead == 0:
            print("NB: had no SHEAD", str(self.SUSY))
            return
        
        # Replace with one common
        if mode==0: 
            for iS in range(nshead-1):
                del self.SUSY[0]
                del self.SM[0]
            self.SUSY[0] = 'SHEAD'
            self.SM[0] = ['-']
            self.SM[1] = ['-']  # drop all info from SHEAD (except that there is one)
        # Or remove entirely
        if mode==1: 
            for iS in range(nshead):
                del self.SUSY[0]
                del self.SM[0]
            self.SM[0] = ['-']


    def Last(self):
        if self.N() == 0:
            return '--'
        else:
            return self.SUSY[self.N()-1]


    def myCopy(self):
        # newCasc = Cascade(self.SUSY[0], self.CODE[0])
        newCasc = Cascade(self.SUSY[0])
        for iP in range(1,self.N()):
            # newCasc.Add(self.SUSY[iP], self.CODE[iP], self.SM[iP], self.BR[iP])
            newCasc.Add(self.SUSY[iP], self.SM[iP], self.BR[iP])
        newCasc.getnDilep()  #not entire classification, but ...
        return newCasc


    def Show(self,mode=0,pure=0):
        if mode==0:
            if pure==0: txt = "%6.4f  %-40s  %-40s" %(self.BRtot, str(self.SUSY), str(self.SM))
            if pure==1: txt = " %-40s  %-40s" %(str(self.SUSY), str(self.SM))
        if mode==1:
            zSM=[]
            for iD in range(len(self.SM)): zSM.append(self.SM[iD][0])
            if pure==0: txt = "%7.5f  %-40s  %-40s" %(self.BRtot, str(self.SUSY), str(zSM))
            if pure==1: txt = " %-40s  %-40s" %(str(self.SUSY), str(zSM))
        if mode==2:
            if pure==0: txt = "%7.5f  %-40s  %-40s" %(self.BRtot, str(self.SUSY), self.strSMp(mode=1))
            if pure==1: txt = " %-40s  %-40s" %(str(self.SUSY), self.strSMp(mode=1))
            
            
        return txt

    def JoinlL(self):
        for iS in range(self.N()):
            if self.SUSY[iS] in ['eL','mL']: self.SUSY[iS] = 'lL'

    def JoinlR(self):
        for iS in range(self.N()):
            if self.SUSY[iS] in ['eR','mR']: self.SUSY[iS] = 'lR'

    def Joinsl(self):
        for iS in range(self.N()):
            if self.SUSY[iS] in ['eL','mL','eR','mR','lR','lL']: self.SUSY[iS] = 'sl'

    def JoinSL(self):
        for iS in range(self.N()):
            if self.SUSY[iS] in ['eL','mL','eR','mR','lR','lL','T1','T2','sl']: self.SUSY[iS] = 'SL'

    def JoinQs(self,pJoin='SQ', start=0, keep=[]):  # for 2-body decay;
        # start parameter allows for e.g. not changing the initial sparticle, of relevance to folding in of xsec
        # 1) Rewrite,  2) Join
        aJoin = {}
        aJoin["qL"] = ["dL","uL","sL","cL"]
        aJoin["qR"] = ["dR","uR","sR","cR"]
        # aJoin["sq"] = ["qL","qR"]
        aJoin["sq"] = ["qL","qR","dL","uL","sL","cL","dR","uR","sR","cR"]
        aJoin["sb"] = ["b1","b2"]
        #aJoin["SQ"] = ["qL","qR","b1","b2","sq"]
        aJoin["SQ"] = ["qL","qR","b1","b2","sq","t1","t2"]


        for iS in range(start,self.N()):
            if self.SUSY[iS] in aJoin[pJoin]:
                self.SUSY[iS] = pJoin

                # THEN DO SM: 

                # ############################# (DEBUG)
                if 0 and iS>0:
                    parent = self.SUSY[iS-1]
                    if isGluino(parent):
                        self.SM[iS] = ['q']
                    elif isNC(parent):
                        self.SM[iS] = ['q']
                    elif isSquark(parent):
                        pass
                    else: 
                        print("Interesting parent: ", self.SM[iS-1])
                # #############################

                # NB:b-removal
                
                if self.SM[iS][0] in setq: self.SM[iS]=['q']
                if pJoin in ['sb','SQ'] and self.SM[iS][0] in setQ: self.SM[iS] = ['q']

                if self.SM[iS+1][0] in setq: self.SM[iS+1]=['q']
                if pJoin in ['sb','SQ'] and self.SM[iS+1][0] in setQ: self.SM[iS+1] = ['q']


                # adhoc for sb->top  ... took out Apr20 for 2leg
                if ('sbtop' in keep) and pJoin in ['SQ'] and (self.SM[iS][0]=='t' or self.SM[iS+1][0]=='t'):
                    self.SUSY[iS] = 'tb'
                if ('sttop' in keep) and pJoin in ['SQ'] and (self.SM[iS][0]=='t' or self.SM[iS+1][0]=='t'):
                    self.SUSY[iS] = 'tb'
                    
        

    def JoinQQ(self,pJoin='QQ',mode=['3body','Z']):  # for 3-body decay 
        # setq = ['d','u','s','c']
        # setQ = ['d','u','s','c','b']
        # setj = ['d','u','s','c','b','T']
        for iS in range(self.N()):
            if '3body' in mode and len(self.SM[iS][0])==2: 
                if self.SM[iS][0] in ['H+']: continue #adhoc: non-3body with 2 characters
                #print self.SM[iS]
                p1 = self.SM[iS][0][0] 
                p2 = self.SM[iS][0][1] #this works because is string
                if pJoin=='qq' and p1 in setq and p2 in setq:
                    self.SM[iS][0] = 'qq'
                if pJoin=='QQ' and p1 in setQ and p2 in setQ:
                    self.SM[iS][0] = 'QQ'
                if pJoin=='jj' and p1 in setj and p2 in setj:
                    self.SM[iS][0] = 'jj'

            # if 'Z' in mode and len(self.SM[iS])==1 and 'Z->' in self.SM[iS][0]:
            if 'Z' in mode and 'Z->' in self.SM[iS][0]:
                # if pJoin=='qq': no action
                if pJoin=='QQ' and self.SM[iS][0]=='Z->qq': self.SM[iS][0]='Z->QQ'
                if pJoin=='QQ' and self.SM[iS][0]=='Z->bb': self.SM[iS][0]='Z->QQ'
                if pJoin=='jj' and self.SM[iS][0]=='Z->QQ': self.SM[iS][0]='Z->jj'
                if pJoin=='jj' and self.SM[iS][0]=='Z->qq': self.SM[iS][0]='Z->jj'
                if pJoin=='jj' and self.SM[iS][0]=='Z->bb': self.SM[iS][0]='Z->jj'
                if pJoin=='jj' and self.SM[iS][0]=='Z->TT': self.SM[iS][0]='Z->jj'

            # Not needed for W since there is nothing to join (starts as QQ and Tv is different)
            # if 'W' in mode and 'W->' in self.SM[iS][0]: 

            

    def JoinQQ_obsolete(self,pJoin='QQ',mode=['3body','Z']):  # for 3-body decay 
        for iS in range(self.N()):
            if '3body' in mode and len(self.SM[iS])==2:
                p1 = self.SM[iS][0]
                p2 = self.SM[iS][1]
                if pJoin=='qq' and p1 in setq and p2 in setq:
                    self.SM[iS][0] = 'q'
                    self.SM[iS][1] = 'q'
                if pJoin=='QQ' and p1 in setQ and p2 in setQ:
                    self.SM[iS][0] = 'Q'
                    self.SM[iS][1] = 'Q'
                if pJoin=='jj' and p1 in setj and p2 in setj:
                    self.SM[iS][0] = 'j'
                    self.SM[iS][1] = 'j'

            if 'Z' in mode and len(self.SM[iS])==1 and 'Z->' in self.SM[iS][0]:
                # if pJoin=='qq': no action
                if pJoin=='QQ' and self.SM[iS][0]=='Z->qq': self.SM[iS][0]='Z->QQ'
                if pJoin=='QQ' and self.SM[iS][0]=='Z->bb': self.SM[iS][0]='Z->QQ'
                if pJoin=='jj' and self.SM[iS][0]=='Z->QQ': self.SM[iS][0]='Z->jj'
                if pJoin=='jj' and self.SM[iS][0]=='Z->qq': self.SM[iS][0]='Z->jj'
                if pJoin=='jj' and self.SM[iS][0]=='Z->bb': self.SM[iS][0]='Z->jj'
                if pJoin=='jj' and self.SM[iS][0]=='Z->TT': self.SM[iS][0]='Z->jj'


    def JoinTtoj(self):
        for iS in range(self.N()):
            for iD in range(len(self.SM[iS])): 
                if self.SM[iS][iD] == 'T': self.SM[iS][iD] = 'j'

    def JoinNs(self):
        for iS in range(self.N()):
            if isN(self.SUSY[iS]): self.SUSY[iS] = "N"

    def JoinCs(self):
        for iS in range(self.N()):
            if isC(self.SUSY[iS]): self.SUSY[iS] = "C"

    def JoinNCs(self):
        for iS in range(self.N()):
            if isNC(self.SUSY[iS]): self.SUSY[iS] = "X"

    def JoinMassive(self): # assume SimplifyNotation() has been done
        for iS in range(self.N()):
            if self.SM[iS][0] in setMassiveDaughter: self.SM[iS][0]='M'

    def SimplifyNotation(self):
        self.SimplifyNotationSM()
        self.SimplifyNotationSUSY()

    def SimplifyNotationSM(self): #Sep21: split SM and SUSY 
        #print self.SUSY, "         ", self.SM   #debugging

        # SM-particles
        for iS in range(self.N()):
            for ii in range(len(self.SM[iS])):
                if self.SM[iS][ii] in ['nue','num','nut']: self.SM[iS][ii] = 'v'
                if self.SM[iS][ii] in ['ta']: self.SM[iS][ii] = 'T'
                if self.SM[iS][ii] in ['e']: self.SM[iS][ii] = 'l'  #assume e/m join
                if self.SM[iS][ii] in ['m']: self.SM[iS][ii] = 'l'
                if self.SM[iS][ii] in ['pho']: self.SM[iS][ii] = 'y'  #poor-man's 'gamma'


            if len(self.SM[iS])==3:  #typical for PI
                self.SM[iS][0] = self.SM[iS][0]+self.SM[iS][1]+self.SM[iS][2]
                self.SM[iS].pop(2)
                self.SM[iS].pop(1)
                
            if len(self.SM[iS])==2:
                if self.SM[iS][0] == 'v':
                    self.SM[iS][0] = self.SM[iS][1]+self.SM[iS][0]  # to get 'ev' not 've' etc
                else: 
                    self.SM[iS][0] = self.SM[iS][0]+self.SM[iS][1]
                self.SM[iS].pop(1)
        
    def SimplifyNotationSUSY(self):
        # SUSY-particles # Apr17 (not fully tested for possible dependence on snue etc.)
        for iS in range(self.N()):
            if self.SUSY[iS] in ['eL','mL']: self.SUSY[iS] = 'lL'
            if self.SUSY[iS] in ['ta1']: self.SUSY[iS] = 'T1'
            if self.SUSY[iS] in ['ta2']: self.SUSY[iS] = 'T2'
            if self.SUSY[iS] in ['snue','snum']: self.SUSY[iS] = 'vl'
            if self.SUSY[iS] in ['snut']: self.SUSY[iS] = 'vT'


    def RenameSleptons(self,mode=0):  # might need refinement and renaming to RenameSleptons  In use as of Apr17 ?
        for iS in range(self.N()):
            if mode >= 0:  # different naming above ... 
                #if self.SUSY[iS] == 'snue': self.SUSY[iS] = 'SN'
                if self.SUSY[iS] in ['ve','vm','vl','vT']: self.SUSY[iS] = 'SN'
                if self.SUSY[iS] in ['eR','eL','mL','mR','lL','lR']: self.SUSY[iS] = 'SL'
                #if self.SUSY[iS] == 'eL': self.SUSY[iS] = 'SL'

            if mode >= 1:  # join slepton/sneutrino
                if self.SUSY[iS] == 'SL': self.SUSY[iS] = 'LN'  
                if self.SUSY[iS] == 'SN': self.SUSY[iS] = 'LN'
        

    # Sep16
    def KeepOnlyFirstAndLastSparticle(self):
        for iS in range(self.N()-1-1,-1+1,-1): self.SUSY.pop(iS)  #backwards, skipping first and last, popping all in between

            
    def TidySMp(self,keep,drop=[],warntxt=""):
        # First add
        for key in keep:
            if key not in list(self.SMp.keys()): self.SMp[key] = 0

        #drop = ['v','y']
        #drop = ['v']
        #drop = []

        # Then find and remove zeros
        remove = []
        for key in list(self.SMp.keys()):
            if self.SMp[key] == 0 and key not in keep:
                remove.append(key)
            if key in drop and key not in keep:  #allow for removing 'v' (and 'y'?)
                remove.append(key)

        # the actual removing
        for key in remove: self.SMp.pop(key)

        # Check:
        if 1:
            #print "\nkeep: %s" %str(set(keep))
            #print "SMp:  %s" %str(set(self.SMp))
            if set(keep) != set(self.SMp): 
                print("WARNING [%s]: 1leg TidySMp: SMp = %s" %(warntxt,set(self.SMp)))
            
            

class Particle:
    def __init__(self,ID, name, massT, lifetime):
        self.ID = ID
        self.name = name
        self.massT = massT
        self.mass = float(massT)
        self.lifetime = lifetime
        if lifetime > 0: # hack: if lifetime is -1 (or <=0), then the particle is stable  # 2013-11-30
            self.width = 6.582e-25/self.lifetime
            self.distance = 3.0e8*lifetime  #t[0]
        else:
            self.width = -1.
            self.distance = -1.
 
        self.decayList = []
        self.NdecayList = 0 #set manually later

    def showDecays(self):
        BRsum = 0.
        for iD in range(len(self.decayList)):
            BRsum += self.decayList[iD].BR
            print(' %2i  %5.3f [%5.3f]  %4s  ->  %4s + %s' %(iD, self.decayList[iD].BR, BRsum, self.name, self.decayList[iD].SUSY[0], str(self.decayList[iD].SM)))
            
    def findDecay(self,children):
        hits = []
        ch = children.split(",")
        #print ch
        su = ch.pop(0)
        #print su
        sus = su.split("|")
        if ch: sms = ch[0].split("|")
        else: sms = []
        #print sms
        
        for iD in range(len(self.decayList)):
            # match susy part
            
            #if self.decayList[iD].SUSY[0] == su:
            if sus == [] or self.decayList[iD].SUSY[0] in sus:

                # NB: checking only first SMchild (more not needed, I think)
                if self.decayList[iD].SM[0] in sms or sms == []:

                    hits.append(iD)
                    

                # Old way: loop over children (new is to use only first child, can then OR)
#                 SM = self.decayList[iD].SM
#                 # match SM part
#                 matchall = 1
                
#                 # loop over requested children
#                 for ich in range(len(ch)):
#                     matchthis = 0
#                     # test against children in given decay
#                     for iSM in range(len(SM)): 
#                         if ch[ich] == SM[iSM]:
#                             matchthis = 1
#                             break
#                     matchall *= matchthis

#                 if matchall: hits.append(iD)
        ###
        return hits  #can be empty, can have several
    
        
    def BR(self,children,SMseparate=""):
        # can take arg as (SU,SM) or both in one
        # several SU can be given e.g. as 'N1|N2' (then summed)
        # For SM
        if SMseparate: children += "," + SMseparate
        hits = self.findDecay(children)
            
        BRtot = 0.
        for hit in hits:
            BRtot += self.decayList[hit].BR
        #if 1 and len(hits) > 1: print "Particle:BR: found %i matching decays" %(len(hits))
        return BRtot
    

# Uh, ugly hack to have QuarkMode and LeptonMode available   ... can put inside Decay() ?? 
QuarkMode = 0
LeptonMode = 0

class Decay: 
    def __init__(self,BR, n, dauGen):
        self.BR = BR
        self.n = n
        self.dauGen = dauGen  #is an ordered array of daughter codes
        self.dauName = []
        self.SUSY = []
        self.SM = []
        self.stable = (n==0)
        for iD in range(n):
            self.dauName.append(fname(dauGen[iD],QuarkMode,LeptonMode))
            if iD==0: self.SUSY.append(self.dauName[iD])
            if iD> 0: self.SM.append(self.dauName[iD])

        
class CascadeStore:
    def __init__(self,BRcut=-1.):
        self.list = []
        self.fn_dilep = "casc_dilep.txt"
        self.BRtot = 0.   # sum of all 1leg-BRs. It should equal the number of initial sparticles allowed minus loss due to BRcutoff. 
        self.BRcut = BRcut
        self.initial = []

    def AllPart_Coverage(self,parts=[]):
        useparts = list(self.UniquePart())
        if parts: useparts = list(set(useparts) & set(parts))
        # easytest:
        for part in parts:
            if part not in useparts:
                print("NB: Not in cascStore: particle %s" %part)

        BRs = {}
        N = {}
        # CALC
        for casc in self.list:
            p = casc.SUSY[0]
            if p in useparts:
                if p not in list(BRs.keys()):
                    BRs[p] = 0.
                    N[p] = 0
                    
                BRs[p] += casc.BRtot
                N[p] += 1

        # SHOW
        for p in useparts:
            print("%2s  %2i  %8.6f" %(p,N[p],BRs[p]))
        

    def NuniquePart(self):
        part = []
        for casc in self.list:
            if casc.SUSY[0] not in part: part.append(casc.SUSY[0])
        return len(part)

    def UniquePart(self):
        part = []
        for casc in self.list:
            if casc.SUSY[0] not in part: part.append(casc.SUSY[0])
        return part


    def GivenPart_N(self,p):
        n = 0
        for casc in self.list:
            if casc.SUSY[0] == p: n += 1
        return n

 
    def GivenPart_BR(self,p):
        BR = 0.
        for casc in self.list:
            if casc.SUSY[0] == p: BR += casc.BRtot
        return BR

    def GivenPart_Show(self,p):
        n = 0
        for casc in self.list:
            if casc.SUSY[0] == p: print(casc.Show())
        return n
               
        

    def BRaverage(self):
        if self.Ninitial() == 0:
            print("Zero initial particles")
            return 0.
        return self.BRtot / self.Ninitial()

    def Ninitial(self):
        return len(self.initial)

    def N(self):
        return len(self.list)


    def Status(self,txt="",ret=False):
        line = 'BRtot = %6.4f ,  N(1legs) = %4i  ,  BRcut = %7.1e   %s'  %(self.BRtot, self.N(), self.BRcut, txt)
        if ret: return line
        else: print(line)
        

    pureShow = 0
    def Show(self,i0=0,i1=-1,mode=0):
        if i1<0: i1=len(self.list)
        i1 = min(i1,len(self.list))
        print("     BR      SUSY CHAIN                                OUTGOING SM PARTICLES")
        print(95*"-")
        for iC in range(i0,i1):
            c = self.list[iC]
            if self.pureShow==0: print("%3i  %s" %(iC,c.Show(mode=mode,pure=self.pureShow)))
            if self.pureShow==1: print("%s" %(c.Show(mode=mode,pure=self.pureShow)))



        

    def Save(self,fn="",mode=0):
        if fn=="": fn=self.fn_dilep
        f = open(fn,"w")
        for iC in range(len(self.list)):
            c=self.list[iC]

            if mode==0: txt = "%4i   %7.5f   %-40s  %-40s" %(iC,c.BRtot,c.strSUSY(),c.strSM())
            if mode==2: txt = "%4i   %7.5f   %-40s  %-40s" %(iC,c.BRtot,c.strSUSY(),c.strSMp(mode=1))
            f.write("%s\n" %(txt))
            # print tex
        f.close()
    
    def Read(self,fn="",verb=0):  # Reads *one* casc_*_<wig> file (need to use outer loop)
        if fn=="": fn=self.fn_dilep
        f = open(fn); lines = f.readlines(); f.close()
        if verb==1: print("%4i cascades in %s" %(len(lines),fn))
        added = 0
        for iL in range(len(lines)):
            word = lines[iL].strip().split()
            if len(word)!=4: print("Warning in Read(): lines[%i] = %s" %(iL,lines[iL].strip()))
            num = int(word[0])
            BRtot = float(word[1])
            zSUSY = word[2]
            zSM = word[3]
            SUSY = zSUSY.split("__")
            SM = zSM.split("__")
            
            cascNew = Cascade(SUSY[0])
            for iS in range(1,len(SUSY)): cascNew.Add(SUSY[iS],[SM[iS]])  #NB legacy for SM
            cascNew.BRtot = BRtot
            cascNew.getnDilep()  #not entire classification, but ...

            added += self.AddIfNewOrBetter(cascNew)
            
        if verb==1: print("%4i added" %(added))


    def Add(self,casc):
        self.list.append(casc)   #not yet in use
        self.BRtot += casc.BRtot


    def AddIfNewOrBetter(self,cascNew):
        N = len(self.list)
        for iC in range(N):
            if self.list[iC].SUSY == cascNew.SUSY and self.list[iC].SM == cascNew.SM:
                if self.list[iC].BRtot < cascNew.BRtot:
                    self.list[iC].BRtot = cascNew.BRtot #replace BRtot, keep rest
                    #print "better"
                    return 1
                return 0
        # if arrive here, cascNew is new
        self.Add(cascNew)
        #print "new"
        return 1
        
    

    def RedoStrSUSYSM(self):
        for iC in range(self.N()): self.list[iC].RedoStrSUSYSM()

    def OrderByBR(self): 
        self.list.sort(key=lambda obj: -obj.BRtot)

    def OrderBySUSY(self):  #doesn't work ...
        self.list.sort(key=lambda obj: obj.fullSUSY)
        
    def OrderBySM(self):  #fails
        self.list.sort(key=lambda obj: -obj.fullSM)

        
    def ShowOrdered(self,i0=0,i1=-1):
        if i1<0: i1=len(self.list)
        i1 = min(i1,len(self.list))
        listOrdered = []
        for L in self.list: listOrdered.append(L)
        listOrdered.sort(key=lambda obj: -obj.BRtot)
        for iC in range(i0,i1):
            c = listOrdered[iC]
            print("%3i  %s" %(iC,c.Show()))
        # could now choose to keep listOrdered
        

    def SubSet(self,i0=0,i1=20):
        i0 = max(0,i0)
        i1 = min(i1,self.N())
        storeSubSet = CascadeStore()
        for iC in range(i0,i1):
            storeSubSet.Add(self.list[iC])
        return storeSubSet

    def Dileptons(self,g):  # work-in-progress
        storeDilep = CascadeStore()
        for iC in range(len(self.list)):
            zcasc = self.list[iC]
            if zcasc.isDilep > 0: 

                if zcasc.BRtot < g.BRtotMin:
                    continue

                storeDilep.Add(self.list[iC])
            
        return storeDilep

    def JoinlL(self):
        for iC in range(len(self.list)): self.list[iC].JoinlL()

    def JoinlR(self):
        for iC in range(len(self.list)): self.list[iC].JoinlR()

    def Joinsl(self):
        for iC in range(len(self.list)): self.list[iC].Joinsl()

    def JoinSL(self):
        for iC in range(len(self.list)): self.list[iC].JoinSL()

    def JoinQs(self,pJoin='SQ', start=0):  # start: allows for not changing first sparticle (of relevance to xsec)
        for iC in range(len(self.list)):
            self.list[iC].JoinQs(pJoin,start)            

    def JoinQQ(self,pJoin='QQ'):
        for iC in range(len(self.list)):
            self.list[iC].JoinQQ(pJoin)

    def JoinTtoj(self):
        for iC in range(len(self.list)):
            self.list[iC].JoinTtoj()

    def JoinNs(self):
        for iC in range(len(self.list)): self.list[iC].JoinNs()
        
    def JoinCs(self):
        for iC in range(len(self.list)): self.list[iC].JoinCs()
        
    def JoinNCs(self):
        for iC in range(len(self.list)): self.list[iC].JoinNCs()
        
    def JoinMassive(self):
        for iC in range(len(self.list)): self.list[iC].JoinMassive()
        
    def eraseSM(self):
        for iC in range(len(self.list)): self.list[iC].eraseSM()
            

    def Duplicates(self,mode='join',modeSM='SM'):  #Sep20: modeSM
        iC1 = 0
        while iC1 < len(self.list)-1:
            # print iC1
            # for iC1 in range(0,len(self.list)):
            self.list[iC1].BRtotBefDupl = self.list[iC1].BRtot  # not used
            duplicates = []
            for iC2 in range(iC1+1,len(self.list)):

                if modeSM == 'SM' and self.list[iC1].SUSY == self.list[iC2].SUSY and self.list[iC1].SM == self.list[iC2].SM: 
                     duplicates.append(iC2)
                if modeSM == 'SMp' and self.list[iC1].SUSY == self.list[iC2].SUSY and self.list[iC1].SMp == self.list[iC2].SMp: 
                     duplicates.append(iC2)

            duplicates.reverse()

            for iC2 in duplicates: 
                # print self.list[iC2].Show()+'     DUPLICATE'
                if mode=='join':  self.list[iC1].BRtot += self.list[iC2].BRtot
                if mode=='drop':  self.list[iC1].BRtot = max(self.list[iC1].BRtot,self.list[iC2].BRtot) 
                self.list.pop(iC2)

            iC1 += 1
            

    def SimplifyNotation(self,n=-1):
        for iC in range(len(self.list)):
            if n>=0 and iC>n: continue       # for debugging
            self.list[iC].SimplifyNotation()

    def SimplifyNotationSM(self,n=-1):
        for iC in range(len(self.list)):
            if n>=0 and iC>n: continue       # for debugging
            self.list[iC].SimplifyNotationSM()

    def SimplifyNotationSUSY(self,n=-1):
        for iC in range(len(self.list)):
            if n>=0 and iC>n: continue       # for debugging
            self.list[iC].SimplifyNotationSUSY()

    def RenameSleptons(self,mode=0):
        for iC in range(self.N()):
            self.list[iC].RenameSleptons(mode)
            

    def CollapseShead(self):
        for iC in range(self.N()):
            self.list[iC].CollapseShead()


    def Classify(self,categ=[]):
        # do base classification of all cascades in cascStore
        for iC in range(self.N()):
            casc = self.list[iC]
            casc.Classify(categ)
            

    def countCasctypes(self):
        # find all casctypes
        ctypes = {}
        alltypes = []
        alltypesN = []
        for iC in range(self.N()):
            casc = self.list[iC]

            # Fill ctypes
            for key in list(casc.casctype.keys()):
                if key in ctypes:
                    if casc.casctype[key] not in ctypes[key]:
                        ctypes[key].append(casc.casctype[key])
                else:
                    ctypes[key] = [casc.casctype[key]]

            # Fill alltypes (and occurence)
            if casc.casctype not in alltypes:
                alltypes.append(casc.casctype)
                alltypesN.append(1)
            else:
                for iT in range(len(alltypes)):
                    if casc.casctype == alltypes[iT]:
                        alltypesN[iT] += 1
                        break

        # Not very elegant...
        self.alltypes = alltypes
        self.alltypesN = alltypesN
                
        # sort
        for key in list(ctypes.keys()):
            ctypes[key].sort()


            
    def cascOfType(self,casctypeA=[]):
        zcascStore = CascadeStore()
        for iC in range(self.N()):
            casc = self.list[iC]
            if casc.isOfType(casctypeA):
                zcascStore.Add(casc)
                #print "yes"
        return zcascStore
            

    def KeepOnlyFirstAndLastSparticle(self):
        for iC in range(len(self.list)):
            casc = self.list[iC]
            casc.KeepOnlyFirstAndLastSparticle()
            

    def ExpandWZt(self,pOlds=['t','W','Z']):

        expanded = False
        KeepOn = True
        while KeepOn: 
            KeepOn = False
            iC = 0
            while iC < self.N():

                for pOld in pOlds:

                    if pOld in self.list[iC].SM:
                        # Have found t, W or Z in this 1leg
                        expanded = True
                        
                        casc = self.list[iC]
                        # Remove this 1leg from list
                        self.list.pop(iC)

                        # Loop over the possible decays: 
                        for dec in ch[pOld]:
                            cascNew = casc.Copy()
                            cascNew.BR *= dec[0]

                            # Replace the first instance of pOld with the daughter(s)

                            for ipnew in range(1,len(dec)):
                                pass #unfinished


    def FillSMp(self):
        for iC in range(len(self.list)):
            casc = self.list[iC]
            casc.FillSMp()
            
    def TidySMp(self,keep,drop=[],warntxt=""):
        for iC in range(len(self.list)):
            casc = self.list[iC]
            casc.TidySMp(keep,drop,warntxt)
        
        
    def ExpandSMp(self,mode=['h']):  #Sep19: this is CascadeStore (there is another one, first made, at TwoLegs
        # mode h: h->bb (1) [assumption]
        # mode W: h->W->lv,Tv,QQ (3)
        # mode Z: h->Z->ll,TT,vv,QQ
        # mode t: bW (so might want to call before mode W)
        # mode b: b->q
        # mode T: T->l,Q,..

        # ###### NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB 
        # ######
        # ###### WE ARE NOW IN CascadeStore::ExpandSMp()  [making of 1-legs]
        # ###### THIS IS THE CORRECT PLACE TO MAKE (TEST) CHANGES IN WHAT/HOW TO EXPAND
        # ######
        # ###### THERE IS ALSO A TwoLegs::ExpandSMp()
        # ###### CHANGES THERE ARE PERFORMED AFTER CHANGES IN THE 1LEGS, 
        # ###### AND SO WILL IN GENERAL NOT HAVE EFFECT
        # ######
        # ###### NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB NB 
        
        ch = {}
        #ch['t'] = [(PDG.Wtolv,'b','l','v'), (PDG.WtoTv,'b','T','v'), (PDG.Wtoqq,'b','q','q')]

        ch['h'] =  [[1.,'b','b']]  #NB: approx-hack
        ch['H'] =  [[1.,'b','b']]  #NB: approx-hack
        ch['A'] =  [[1.,'b','b']]  #NB: approx-hack
        ch['H+'] = [[1.,'t','b']]  #NB: approx-hack

        ch['t'] = [[1.,'b','W']]
        ch['W'] = [[PDG.Wtolv,'l','v'], [PDG.WtoTv,'T','v'], [PDG.Wtoqq,'q','q']]
        
        #ch['W'] = [[1.,'q','q']]  # Hack to assess importance of W decay for number of leptons [2011-01-10]
        #                         # In SU4 this is sufficient to isolate the effect of stop production
        #                         # In other scenarios you may need to approach it differently since there might be other W sources
        ch['Z'] = [[PDG.Ztoll,'l','l'], [PDG.ZtoTT,'T','T'], [PDG.Ztovv,'v','v'], [PDG.ZtoQQ,'q','q']]
        ch['QQ'] = [[1.,'q','q']] #added Apr26
        #ch['Q'] = [[1.,'q']]      #added Apr26  In use?
        ch['TT'] = [[1.,'T','T']] #added Apr26
        ch['ll'] = [[1.,'l','l']] #added Apr26
        ch['vv'] = [[1.,'v','v']] #added Sep20
        #if 0: 
        #    ch['lv'] = [[1.,'l']] #added Apr26
        #    ch['Tv'] = [[1.,'T']] #added Apr26
        #    # Sep16
        #    ch['blv'] = [[1.,'b','l']] #added Apr26
        #    ch['bTv'] = [[1.,'b','T']] #added Apr26
        #    ch['bqq'] = [[1.,'b','q','q']] #added Apr26

        #Sep19: keep v
        ch['lv'] = [[1.,'l','v']] #added Apr26
        ch['Tv'] = [[1.,'T','v']] #added Apr26
        ch['blv'] = [[1.,'b','l','v']] #added Apr26
        ch['bTv'] = [[1.,'b','T','v']] #added Apr26
        ch['bqq'] = [[1.,'b','q','q']] #added Apr26

        # Sep20
        ch['tb'] = [[1.,'t','b']] 
        ch['u'] = [[1.,'q']] 
        ch['d'] = [[1.,'q']] 
        ch['g'] = [[1.,'q']] 

        
        ch['b'] = [[1.,'q']]
        ch['T'] = [[PDG.Ttolvv,'l','v','v'], [PDG.Ttovh,'q','v']]

        #ch['q'] = [[1.,'v']]  # Hack
        ch['q'] = [[1.]]

        #Sep21: further simplifications (to skip some of the old simplification routines)
        ch['uu'] = [[1.,'q','q']]
        ch['dd'] = [[1.,'q','q']]
        ch['ud'] = [[1.,'q','q']]
        ch['bb'] = [[1.,'b','b']]

        ch['tt'] = [[1.,'t','t']]  #2011-01-03

        ch['l'] = [[1.]]  # 2012-06-23

        expanded = False
        KeepOn = True
        while KeepOn: # this outermost loop is now probably not needed

            KeepOn = False
            iTL = 0
            # Iterate through list. Will continuously expand, so use dynamical self.N()
            while iTL < self.N():
                ##print 'while %i < %i' %(iTL, self.N())

                #pOlds = SMp_keys2 + ['t','h','W','Z','b','T']  #Apr26  "all" possibilities to transform
                #pOlds = SMp_keys2 + ['t','h','W','Z','b','T','q']  #Sep19  "all" possibilities to transform
                pOlds = SMp_keys2 + ['tt','t','h','A','H','H+','W','Z','b','T','q']  #Sep19  "all" possibilities to transform  #2011-01-03
                pOlds += ['tb','u','d','g']  #Sep20
                pOlds += ['uu','dd','ud','bb'] #Sep21

                for pOld in pOlds: 
                    if pOld in mode:
                        
                        if pOld in list(self.list[iTL].SMp.keys()): 
                            if self.list[iTL].SMp[pOld] == 0: continue #adhoc

                            expanded = True
                            TL = self.list[iTL]

                            # Remove list item (but still have it in TL)
                            self.list.pop(iTL)
                            # Change TL
                            ##print 'TL: %s  %i -> ' %(pOld,TL.SMp[pOld]),
                            ###print "OLD:   ", TL.SMp, '   ', pOld
                            TL.SMp[pOld] -= 1
                            ##print '%s' %(TL.SMp[pOld])
                        
                            # Loop over transformations: create new objects, adjust, insert
                            for dec in ch[pOld]:
                                # new object
                                ##print 'TL   : %s %i' %(pOld,TL.SMp[pOld])
                                #TLnew = TL.Copy()
                                TLnew = TL.CopySep19()
                                #print "%3s" %(pOld)
                                ##print 'TL   : %s %i' %(pOld,TL.SMp[pOld])
                                ##print 'TLnew: %s %i' %(pOld,TLnew.SMp[pOld])
                                # adjust BR
                                #TLnew.BR *= dec[0]
                                TLnew.BRtot *= dec[0]
                                # add the new particles (of decay/transformation)
                                for ipnew in range(1,len(dec)):
                                    pnew = dec[ipnew]
                                    if pnew in list(TLnew.SMp.keys()): TLnew.SMp[pnew] += 1
                                    else: TLnew.SMp[pnew] = 1
                                # insert in list (a possibility is to append ... might be faster)
                                ###print "NEW:   ", TLnew.SMp
                                self.list.insert(iTL,TLnew)
                                #self.list.append(TLnew)
                            # continue without increasing iTL, i.e. redo the same entry (for possible additional changes)
                            ###print
                            del TL
                            continue

                        # (now do next pOld happening to be in mode)

                    # (now try next pOld)

                #
                iTL += 1

        return expanded  # so that I know if an expansion took place

    #def MakeSMlist(self,mode):
    #    while KeepOn: # this most-outside loop is now probably not needed
    #        KeepOn = False

        
        

# #################################################

PDG = FreeObject()

# ### PARTICLE DATA GROUP
PDG.Ztoee = 3.363e-2
PDG.Ztomm = 3.366e-2
PDG.ZtoTT = 3.370e-2
PDG.Ztovv = 20.00e-2
#PDG.ZtoQQ = 69.91e-2
PDG.ZtoQQ = 1 - (PDG.Ztoee + PDG.Ztomm + PDG.ZtoTT + PDG.Ztovv)   #=69.901e-2  #Sep20
PDG.Ztobb = 15.12e-2

PDG.Ztoll = PDG.Ztoee + PDG.Ztomm
PDG.Ztoqq = PDG.ZtoQQ - PDG.Ztobb

PDG.Wtoev = 10.75e-2
PDG.Wtomv = 10.57e-2
PDG.WtoTv = 11.25e-2
#PDG.WtoQQ = 67.60e-2
#PDG.Wtoqq = 67.60e-2  # the b-contribution should be negligible
PDG.Wtoqq = 67.43e-2  # the b-contribution should be negligible
# Sep20:had to manually take Wtoqq down to acheive Sum(W) = 1.0


PDG.Wtolv = PDG.Wtoev + PDG.Wtomv

PDG.Ttoevv = 17.85e-2 #+ 1.75e-2  #2011-01-10: small bug discovered. The small numbers are included in the large (but I added them)
PDG.Ttomvv = 17.36e-2 #+ 0.36e-2
PDG.Ttolvv = PDG.Ttoevv + PDG.Ttomvv
PDG.Ttovh  = 1. - PDG.Ttolvv

PDG.ttobW = 1.

# ########################## 



# =========================================================================== May3: copy
# =========================================================================== method from
# =========================================================================== isa2lha.py

def ReadIsawig(fn_isawig, mode=2):
    #from ROOT import cos, sin  #oh well
    #from math import cos, sin

    wig = {}
    mass = {}
    time = {}
    width = {}
    distance = {}

    if not os.path.exists(fn_isawig): # inelegant but ok
        print("WARNING  libISAWIG::ReadIsawig   non-existent file %s" %(fn_isawig))
        if mode==2: return mass, wig
        if mode==3: return mass, wig, time 
        if mode==5: return mass, wig, time, distance, width
        if mode==9: return {'mass':mass, 'wig':wig, 'time':time, 'distance':distance, 'width':width}
        
    
    f_isawig=open(fn_isawig); lines=f_isawig.readlines(); f_isawig.close()

    QuarkMode = 1
    LeptonMode = 1

    # these are in wig-listing, but have zero mass and inf width : 
    nonparticles = [438,440,442,444,446,448]  # is a hack ... also used earlier

    NP = 0

    # 1. READ MASSES
    np = int(lines[NP].strip())  #wig file states how many particle codes there are
    NP += 1
    for iL in range(NP,NP+np): 
        word = lines[iL].strip().split()

        if int(word[0]) in nonparticles: continue  # 2012-05-06

        pname = fname(int(word[0]), QuarkMode, LeptonMode)

        if pname not in list(mass.keys()):
            mass[pname] = float(word[1])
            time[pname] = float(word[2])
            width[pname] = 6.582e-25/time[pname]

            distance[pname] = 3.0e8*time[pname]


    
    # 2. READ DECAYS (first entry per particle is a single-integer-line)
    #print 'DEBUG: ', fn_isawig
    # Method 1
    NP = np
    while 1:
        NP += 1
        word = lines[NP].split()
        if len(word) == 2: break
        

    # Method 2
    while 0:
        print('DEBUG  NP:%i  np:%i' %(NP,np))
        NP += np
        word = lines[NP].split()
        #print "linestart: ",lines[NP],NP

        if len(word) == 2: break   # All decays read. Jump out

        np = int(word[0])
        NP += 1
        #for iL in range(NP+1,NP+np+1):
        for iL in range(NP,NP+np):

            # Read/interpret decay chain
            # (not yet implemented)
            # print 'READING ... ',lines[iL],
            continue



    # 3. READ PARAMETERS: MIXING ETC. (land here when done with decays)
    wig['tan(beta)'] = float(word[0])
    wig['alfah'] = float(word[1])

    if 1: # 2013-06-26
        # print wig['tan(beta)']
        word = lines[NP+1].split()
        wig['n11'] = float(word[0])
        wig['n12'] = float(word[1])
        wig['n13'] = float(word[2])
        wig['n14'] = float(word[3])
        word = lines[NP+2].split()
        wig['n21'] = float(word[0])
        wig['n22'] = float(word[1])
        wig['n23'] = float(word[2])
        wig['n24'] = float(word[3])
        word = lines[NP+3].split()
        wig['n31'] = float(word[0])
        wig['n32'] = float(word[1])
        wig['n33'] = float(word[2])
        wig['n34'] = float(word[3])
        word = lines[NP+4].split()
        wig['n41'] = float(word[0])
        wig['n42'] = float(word[1])
        wig['n43'] = float(word[2])
        wig['n44'] = float(word[3])
        word = lines[NP+5].split()
        wig['vmix11'] = float(word[0])
        wig['vmix12'] = float(word[1])
        wig['vmix21'] = float(word[2])
        wig['vmix22'] = float(word[3])
        word = lines[NP+6].split()
        wig['umix11'] = float(word[0])
        wig['umix12'] = float(word[1])
        wig['umix21'] = float(word[2])
        wig['umix22'] = float(word[3])
        word = lines[NP+7].split()
        wig['thetat']   = float(word[0])
        wig['thetab']   = float(word[1])
        wig['thetatau'] = float(word[2])
        wig['stopmix11']  =  cos(wig['thetat'])
        wig['stopmix12']  =  sin(wig['thetat'])
        wig['stopmix21']  = -sin(wig['thetat'])
        wig['stopmix22']  =  cos(wig['thetat'])
        wig['sbotmix11']  =  cos(wig['thetab'])
        wig['sbotmix12']  =  sin(wig['thetab'])
        wig['sbotmix21']  = -sin(wig['thetab'])
        wig['sbotmix22']  =  cos(wig['thetab'])
        wig['staumix11']  =  cos(wig['thetatau'])
        wig['staumix12']  =  sin(wig['thetatau'])
        wig['staumix21']  = -sin(wig['thetatau'])
        wig['staumix22']  =  cos(wig['thetatau'])
    
        word = lines[NP+8].split()
        wig['a_t']   = float(word[0])
        wig['a_b']   = float(word[1])
        wig['a_tau'] = float(word[2])
        word = lines[NP+9].split()
        wig['mu'] = float(word[0])

    
    if mode==2: return mass, wig
    if mode==3: return mass, wig, time 
    if mode==5: return mass, wig, time, distance, width

    if mode==9: return {'mass':mass, 'wig':wig, 'time':time, 'distance':distance, 'width':width}
    
# =========================================================================== May3: copy
# =========================================================================== method from
# =========================================================================== isa2lha.py

def ReadIsaout(fn_isaout, VB=0):
    if not os.path.exists(fn_isaout):
        if VB: print("WARNING  libISAWIG::ReadIsaout: non-existent file %s" %(fn_isaout))
        return {}
    
    param = {}
    f_isaout=open(fn_isaout); lines=f_isaout.readlines(); f_isaout.close()
    iL=0
    while lines[iL].find("PARENT") < 0:

        # print "%i %s"  %(iL,lines[iL]),
        word = lines[iL].strip().split()
        if len(word)==0:
            iL += 1
            continue

        if lines[iL].startswith(" *  ISAJET"): 
            param['PROGRAM'] = word[1]
            param['PROGRAMVERSION'] = lines[iL].replace("ISAJET","").replace("*","").strip()

        if lines[iL].startswith(" INPUTS FOR"):
            param['PROGRAMMODE'] = word[2].replace(":","")

        if word[0] == 'M(TP)': 
            param['M(TP)']  = float(word[2])


        if word[0] == 'M(GLSS)':
            param['M(GLSS)'] = float(word[2])
            param['MU'] = float(word[5])
            param['M(HA)'] = float(word[8])

        if word[0] == 'TAN(BETA)':
            param['TAN(BETA)'] = float(word[2])
 
        if word[0] == 'M(Q1':
            param['M(Q1)'] = float(word[3])
            param['M(DR)'] = float(word[6])
            param['M(UR)'] = float(word[9])
 
        if word[0] == 'M(L1)':
            param['M(L1)'] = float(word[2])
            param['M(ER)'] = float(word[5])
            param['M(Q3)'] = float(word[8])
 
        if word[0] == 'M(BR)':
            param['M(BR)'] = float(word[2])
            param['M(TR)'] = float(word[5])
            param['M(L3)'] = float(word[8])
 
        if word[0] == 'M(LR)':
            param['M(LR)'] = float(word[2])   #this is m_tauR (Mar10, 2010)
            param['A_T'] = float(word[5])
            param['A_B'] = float(word[8])
 
        if word[0] == 'A_TAU':
            param['A_TAU'] = float(word[2])
 
        if word[0] == 'ALPHAEM':
            param['ALPHAEM'] = float(word[2])
            param['SIN2(THW)'] = float(word[5])
            param['ALPHA3'] = float(word[8])
 
        if word[0] == 'M(Q2)':
            param['M(Q2)'] = float(word[2])
            param['M(SR)'] = float(word[5])
            param['M(CR)'] = float(word[8])
 
        if word[0] == 'M(L2)':
            param['M(L2)'] = float(word[2])
            param['M(MR)'] = float(word[5])
 
        if word[0] == 'M_1':
            param['M_1'] = float(word[2])
            param['M_2'] = float(word[5])
 
        if word[0] == 'NEUTRALINO':
            param['M(N1)s'] = float(word[4])
            param['M(N2)s'] = float(word[5])
            param['M(N3)s'] = float(word[6])
            param['M(N4)s'] = float(word[7])
            param['M(N1)'] = abs(float(word[4]))
            param['M(N2)'] = abs(float(word[5]))
            param['M(N3)'] = abs(float(word[6]))
            param['M(N4)'] = abs(float(word[7]))

        if word[0] == 'EIGENVECTOR' and word[1] == "1":
            param['NCOMP1'] = [float(word[3]),float(word[4]),float(word[5]),float(word[6])]

        if word[0] == 'EIGENVECTOR' and word[1] == "2":
            param['NCOMP2'] = [float(word[3]),float(word[4]),float(word[5]),float(word[6])]

        if word[0] == 'EIGENVECTOR' and word[1] == "3":
            param['NCOMP3'] = [float(word[3]),float(word[4]),float(word[5]),float(word[6])]

        if word[0] == 'EIGENVECTOR' and word[1] == "4":
            param['NCOMP4'] = [float(word[3]),float(word[4]),float(word[5]),float(word[6])]
  
        if word[0] == 'CHARGINO':
            param['M(C1)s'] = float(word[4])
            param['M(C2)s'] = float(word[5])
            param['M(C1)'] = abs(float(word[4]))
            param['M(C2)'] = abs(float(word[5]))
 
        if word[0] == 'GAMMAL,':
            param['GAMMAL'] = float(word[3])
            param['GAMMAR'] = float(word[4])
 
        if word[0] == 'M(HL)':
            param['M(HL)'] = float(word[2])
            param['M(HH)'] = float(word[5])
            param['M(H+)'] = float(word[8])
            param['ALFAH'] = float(word[9].split("=")[1])
 
        if word[0] == 'M(T1)':
            param['M(T1)'] = float(word[2])
            param['M(T2)'] = float(word[5])
            param['THETAT'] = float(word[7])
 
        if word[0] == 'M(B1)':
            param['M(B1)'] = float(word[2])
            param['M(B2)'] = float(word[5])
            param['THETAB'] = float(word[7])
 
        if word[0] == 'M(TAU1)':
            param['M(TAU1)'] = float(word[2])
            param['M(TAU2)'] = float(word[5])
            param['THETATAU'] = float(word[7])


        iL += 1

    return param
# =========================================================================== 
# =========================================================================== 
# =========================================================================== 

