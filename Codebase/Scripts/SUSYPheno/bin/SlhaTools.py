#import pyslha_edited as pyslha
import pyslha_edited as pyslha
import math
from kilelib import GeneaLog, IterateNdim

""" pyslha creates the following structure
slha = pyslha.ReadSLHAFile(fn_slha)
slha is a tuple with 2 entries: 

slha[0] is dict with all the blocks except decays
  slha[0].keys() = ['YE', 'YD', 'MODSEL', 'YU', 'SPINFO', 'STAUMIX', 'NMIX', 'DCINFO', 'MSOFT', 'MASS', 'GAUGE', 'STOPMIX', 'UMIX', 'AE', 'AD', 'SBOTMIX', 'AU', 'SMINPUTS', 'MINPAR', 'EXTPAR', 'VMIX', 'ALPHA', 'HMIX']

slha[1] is a dict with the decays: 
  slha[1].keys() = [2000001, 2000002, 2000003, 2000004, 2000005, 2000006, 2000011, 2000013, 2000015, 25, 35, 36, 6, 1000001, 1000002, 1000003, 1000004, 1000005, 1000006, 1000011, 1000012, 1000013, 1000014, 1000015, 1000016, 1000021, 1000022, 1000023, 1000024, 1000025, 37, 1000035, 1000037]

The objects in slha[0] are pyslha::Block objects: 
    Object representation of any BLOCK elements read from the SLHA file.  Blocks
    have a name, may have an associated Q value, and then a collection of data
    entries, stored as a recursive dictionary. Types in the dictionary are
    numeric (int or float) when a cast from the string in the file has been
    possible.
The block has a name: slha[0][block].name ; ex: slha[0]['MASS'].name (='MASS') ... not very useful
The block has a dict of values: slha[0][block].entries, ex: slha[0]['MASS'].entries[2000001] is mass of dR

The entries dict in MASS block uses the pdgID if the particles as key
The entries dict in all other blocks uses slha-defined keys, ex: slha[0]['MSOFT'].entries[1] is M1 value, [2] is M2




slha[1] (is a dict with pdgID as keys)
type(slha[1][2000001]) is pyslha.Particle : has attributes pid, totalwidth, mass and decays[]
type(slha[1][2000001].decays) is a list of pyslha.Decay objects which has attributes: parentid (unfilled), br, nda, ids[]
  ids is the list of decay products  (and nda is just the number of decay products)


So the Particle and Decay class has a bit of the functionality I need
  The BRs are already read in
  Only need to have routines to treat them
    FullDecays() <-- starting to do some isaPlay stuff
    FindDecay()  FindFullDecay
    BR()

Usable:
libISAWIG:CAscade?

"""


# #############################################################
class FreeObject: pass

PDG = FreeObject()

# ### PARTICLE DATA GROUP
PDG.Ztoee = 3.363e-2
PDG.Ztomm = 3.366e-2   # ee is smaller by 0.09% = 1e-4 ; probably within errors : So Z asymmetry is absent
PDG.ZtoTT = 3.370e-2
PDG.Ztovv = 20.00e-2
#PDG.ZtoQQ = 69.91e-2
PDG.ZtoQQ = 1 - (PDG.Ztoee + PDG.Ztomm + PDG.ZtoTT + PDG.Ztovv)   #=69.901e-2  #Sep20
PDG.Ztobb = 15.12e-2

PDG.Ztoll = PDG.Ztoee + PDG.Ztomm
PDG.Ztoqq = PDG.ZtoQQ - PDG.Ztobb

PDG.Wtoev = 10.75e-2
PDG.Wtomv = 10.57e-2  # ev is 1.7% larger than mv   So W asymmetry is noticeable .. 
PDG.WtoTv = 11.25e-2
#PDG.WtoQQ = 67.60e-2
#PDG.Wtoqq = 67.60e-2  # the b-contribution should be negligible
PDG.Wtoqq = 67.43e-2  # the b-contribution should be negligible
# Sep20:had to manually take Wtoqq down to acheive Sum(W) = 1.0


PDG.Wtolv = PDG.Wtoev + PDG.Wtomv

PDG.Ttoevv = 17.85e-2 #+ 1.75e-2  #2011-01-10: small bug discovered. The small numbers are included in the large (but I added them)
PDG.Ttomvv = 17.36e-2 #+ 0.36e-2     # evv is 2.8% larger than mvv ; noticeable asymmetry
PDG.Ttolvv = PDG.Ttoevv + PDG.Ttomvv
PDG.Ttovh  = 1. - PDG.Ttolvv

PDG.ttobW = 1.




# ##########
pdg2namSigned = { 

    #-1: "d+",
    #-2: "u-",
    #-3: "s+",
    #-4: "c-",
    #-5: "b+",
    #-6: "t-",
    
    -1: "d_",
    -2: "u_",
    -3: "s_",
    -4: "c_",
    -5: "b_",
    -6: "t_",

    -11: "e+",
    -13: "m+",
    -15: "T+",

    #1: "d-",
    #2: "u+",
    #3: "s-",
    #4: "c+",
    #5: "b-",
    #6: "t+",
    
    1: "d",
    2: "u",
    3: "s",
    4: "c",
    5: "b",
    6: "t",

    11: "e-",
    13: "m-",
    15: "T-",
    
    12: "ve",
    14: "vm",
    16: "vT",

    -12: "ve_",
    -14: "vm_",
    -16: "vT_",

    21: "g", 
    22: "y", 
    23: "Z",
    24: "W+",
    -24: "W-",
    
    25: "h",
    26: "h", # 2012-05-06
    35: "H",
    36: "A",     # 2014-01-29 fixed (36 and 37 were swapped)
    37: "Hc+",   # 2014-01-29 fixed
    -37: "Hc-",  # 2014-01-29 fixed

    
    #1000001: "~dL-", 
    #1000002: "~uL+", 
    #1000003: "~sL-", 
    #1000004: "~cL+", 
    #1000005: "~b1-", 
    #1000006: "~t1+", 
    #2000001: "~dR-", 
    #2000002: "~uR+", 
    #2000003: "~sR-", 
    #2000004: "~cR+", 
    #2000005: "~b2-", 
    #2000006: "~t2+",
    
    #-1000001: "~dL+", 
    #-1000002: "~uL-", 
    #-1000003: "~sL+", 
    #-1000004: "~cL-", 
    #-1000005: "~b1+", 
    #-1000006: "~t1-", 
    #-2000001: "~dR+", 
    #-2000002: "~uR-", 
    #-2000003: "~sR+", 
    #-2000004: "~cR-", 
    #-2000005: "~b2+", 
    #-2000006: "~t2-", 


    1000001: "~dL", 
    1000002: "~uL", 
    1000003: "~sL", 
    1000004: "~cL", 
    1000005: "~b1", 
    1000006: "~t1", 
    2000001: "~dR", 
    2000002: "~uR", 
    2000003: "~sR", 
    2000004: "~cR", 
    2000005: "~b2", 
    2000006: "~t2",
    
    -1000001: "~dL_", 
    -1000002: "~uL_", 
    -1000003: "~sL_", 
    -1000004: "~cL_", 
    -1000005: "~b1_", 
    -1000006: "~t1_", 
    -2000001: "~dR_", 
    -2000002: "~uR_", 
    -2000003: "~sR_", 
    -2000004: "~cR_", 
    -2000005: "~b2_", 
    -2000006: "~t2_", 

    1000011: "~eL-", 
    1000013: "~mL-", 
    1000015: "~T1-", 

    2000011: "~eR-",
    2000013: "~mR-", 
    2000015: "~T2-",

    -1000011: "~eL+", 
    -1000013: "~mL+", 
    -1000015: "~T1+", 

    -2000011: "~eR+",
    -2000013: "~mR+", 
    -2000015: "~T2+",

    
    1000012: "~Ve", 
    1000014: "~Vm", 
    1000016: "~VT", 
    -1000012: "~Ve_", 
    -1000014: "~Vm_", 
    -1000016: "~VT_", 

    1000021: "~gl",
    
    1000022: "~N1", 
    1000023: "~N2", 
    1000024: "~C1+", 
    -1000024: "~C1-", 
    1000025: "~N3", 
    1000035: "~N4", 
    1000037: "~C2+", 
    -1000037: "~C2-", 


}

namSigned2pdg = dict( (v,k) for k,v in pdg2namSigned.items() )

pname_generic2signed = {

    'W': ['W+','W-'],
    'Hc': ['Hc+','Hc-'],
    'e':['e+','e-'],
    'm':['m+','m-'],
    'T':['T+','T-'],
    
    #'~dL':['~dL','~dL_'],  # hmm
    
    '~eL':['~eL+','~eL-'],
    '~mL':['~mL+','~mL-'],
    '~TL':['~TL+','~TL-'],
    '~eR':['~eR+','~eR-'],
    '~mR':['~mR+','~mR-'],
    '~TR':['~TR+','~TR-'],

    '~C1':['~C1+','~C1-'],
    '~C2':['~C2+','~C2-'],
}
    
    


# ##########
pdg2nam = {   # include negative too? 

    1: "d",
    2: "u",
    3: "s",
    4: "c",
    5: "b",
    6: "t",

    11: "e",
    13: "m",
    15: "T",  # new: T for all taus
    
    12: "ve",
    14: "vm",
    16: "vT",

    21: "g", 
    22: "y", 
    23: "Z",
    24: "W",
    
    25: "h",
    26: "h", # 2012-05-06
    35: "H",
    36: "A",   # 2014-01-29 fixed (36 and 37 were swapped)
    37: "Hc",  # 2014-01-29 fixed (36 and 37 were swapped)
    
    1000001: "~dL", 
    1000002: "~uL", 
    1000003: "~sL", 
    1000004: "~cL", 
    1000005: "~b1", 
    1000006: "~t1", 
    2000001: "~dR", 
    2000002: "~uR", 
    2000003: "~sR", 
    2000004: "~cR", 
    2000005: "~b2", 
    2000006: "~t2", 

    1000011: "~eL", 
    1000013: "~mL", 
    1000015: "~T1", 

    2000011: "~eR",
    2000013: "~mR", 
    2000015: "~T2",
    
    #1000012: "~Ve",  # could use ~ve ... but keep for now
    #1000014: "~Vm", 
    #1000016: "~VT", 
    1000012: "~ve",   # 2014-08-29
    1000014: "~vm", 
    1000016: "~vT", 

    1000021: "~gl",
    
    1000022: "~N1", 
    1000023: "~N2", 
    1000024: "~C1", 
    1000025: "~N3", 
    1000035: "~N4", 
    1000037: "~C2", 
}

nam2pdg = dict( (v,k) for k,v in pdg2nam.items() )

snames = {}
snames['weakinos']  = ['~N1','~C1','~N2','~N3','~C2','~N4']
#snames['sleptons']  = ['~eR','~mR','~T1','~eL','~mL','~T2', '~Ve','~Vm','~VT']
snames['sleptons']  = ['~eR','~mR','~T1','~eL','~mL','~T2', '~ve','~vm','~vT']   # 2014-08-29
snames['squarks3']  = ['~b1','~b2','~t1','~t2']
snames['squarks12'] = ['~dR','~sR','~uR','~cR','~dL','~sL', '~uL','~cL']
snames['colour'] = ['~gl'] + snames['squarks3'] + snames['squarks12']
snames['higgses']   = ['h','H','A','Hc']
#snames['all'] = snames['weakinos'] + snames['sleptons'] + ['~gl'] + snames['colour']  + snames['higgses']  # 2014-09-02: removed ~gl since is already in 'colour' (gave two entries for ~gl with slha2tlt)
snames['all'] = snames['weakinos'] + snames['sleptons'] + snames['colour']  + snames['higgses']

# ##########
namGeneric = {
    
    'l': ['e','m'], 
    'L': ['e','m','T'],  # l or L
    'q': ['d','u','s','c'], 
    'Q': ['d','u','s','c','b'],  # q or Q
    'v': ['ve','vm','vT'], 
    #'v': ['ve','vm','vt'],   # 2014-02-04  (was this ever in use?)

    '~e' : ['~eL','~eR'],  # some of these can be combined, most not
    '~m' : ['~mL','~mR'],
    '~T' : ['~T1','~T2'],
    '~lL': ['~eL','~mL'],
    '~lR': ['~eR','~mR'],
    '~l' : ['~eL','~mL','~eR','~mR'],
    '~L' : ['~eL','~mL','~eR','~mR', '~T1', '~T2'], 
    '~v' : ['~ve','~vm','~vT'],   # 2014-08-29

    '~dL': ['~dL','~sL'], 
    '~dR': ['~dR','~sR'], 
    '~uL': ['~uL','~cL'], 
    '~uR': ['~uR','~cR'], 

    '~qL': ['~dL','~sL','~uL','~cL'],   # unreliable
    '~qR': ['~dR','~sR','~uR','~cR'],   # unreliable
    '~q' : ['~dL','~sL','~uL','~cL','~dR','~sR','~uR','~cR'],   # unreliable
    '~b' : ['~b1','~b2'],
    '~Q' : ['~dL','~sL','~uL','~cL','~dR','~sR','~uR','~cR','~b1','~b2'],   # unreliable
    
    }


# ---------
"""
SMrewrite = {}

SMrewrite['u'] = [[1.,'q']] 
SMrewrite['d'] = [[1.,'q']] 
SMrewrite['g'] = [[1.,'q']]  # NB: controversial? Rather use j than q? 

SMrewrite['b'] = [[1.,'q']]  # Note: only the ones requested are changed (so by default b is not turned to q)

SMrewrite['c'] = [[1.,'q']]  # 2013-12-06 required after dynamical Higgs expansion
SMrewrite['s'] = [[1.,'q']]  # 2013-12-06 required after dynamical Higgs expansion
SMrewrite['e'] = [[1.,'l']]  # 2013-12-06 required after dynamical Higgs expansion
SMrewrite['m'] = [[1.,'l']]  # 2013-12-06 required after dynamical Higgs expansion
SMrewrite['ve'] = [[1.,'v']] 
SMrewrite['vm'] = [[1.,'v']] 
SMrewrite['vT'] = [[1.,'v']] 


SMrewrite['q'] = [[1.]]
SMrewrite['l'] = [[1.]]  # 2012-06-23
SMrewrite['y'] = [[1.]] 
SMrewrite['v'] = [[1.]] 
"""
# ---------



# ##########
def ParticleName(x, optD={}):  # Ex: p = ParticleName(1000004, optD={'generic'=['l','q']})

    # First find the basic pname
    if type(x) is int: 
        if optD.get('dict','') == 'signed': pname = pdg2namSigned.get(x,'')
        else: pname = pdg2nam.get(x,'')
        if pname == '':
            pname = '?%i?' %(x)
            print('Warning::ParticleName  Unknown pdg %i. Returning %s' %(x, pname))
    else:
        pname = x


    # Then generalise according to options
    for gen in optD.get('generic',[]):
        if pname in namGeneric[gen]:
            pname = gen
            break

    return pname



# ##########
def AddSMtoMothers(geneaLog, genea, particles, PDG, slhaInfo={}, SMrewrite=1, VB=1): 
    # Done undert-the-hood with particles object
    genea = geneaLog.Genea(genea)
    geneaLog.Add(genea=genea, S='info', txt='now starting')
    
    pname = 'W'
    if pname not in particles:
        mass = slhaInfo.get('masses',{}).get(pname, 80.483)
        width = slhaInfo.get('width',{}).get(pname, 2.085)  # this is not in slhaInfo (nor in slha)
        p = Mother(pname, mass, width)
        p.AddFS_brute(['e','v'], PDG.Wtoev)
        p.AddFS_brute(['m','v'], PDG.Wtomv)
        p.AddFS_brute(['T','v'], PDG.WtoTv)
        p.AddFS_brute(['q','q'], PDG.Wtoqq)  # these are light q, nothing to b of course
        particles[pname] = p
        geneaLog.Add(genea=genea, S='info', txt='AddSMtoMothers for %s' %(pname))
        geneaLog.Add(genea=genea, S='info', txt=p.ShowDecays(ret=1))   # text list sent to geneaLog

    pname = 'Z'
    if pname not in particles:
        mass = slhaInfo.get('masses',{}).get(pname, 91.187)
        width = slhaInfo.get('width',{}).get(pname, 2.4952)  # this is not in slhaInfo (nor in slha)
        p = Mother(pname, mass, width)
        p.AddFS_brute(['e','e'], PDG.Ztoee)
        p.AddFS_brute(['m','m'], PDG.Ztomm)
        p.AddFS_brute(['T','T'], PDG.ZtoTT)
        p.AddFS_brute(['v','v'], PDG.Ztovv)
        p.AddFS_brute(['b','b'], PDG.Ztobb) 
        p.AddFS_brute(['q','q'], PDG.Ztoqq) 
        particles[pname] = p
        #if VB>0: print '\nAddSMtoMothers for %s' %(pname) ; p.ShowDecays() 
        geneaLog.Add(genea=genea, S='info', txt='AddSMtoMothers for %s' %(pname))
        geneaLog.Add(genea=genea, S='info', txt=p.ShowDecays(ret=1))   # text list sent to geneaLog

    #  ## ifTfundamental
    pname = 'T' 
    if pname not in particles:
        mass = slhaInfo.get('masses',{}).get(pname, 1.776)
        width = slhaInfo.get('width',{}).get(pname, 0.)  # lifetime 2.906e-13 s # calculate # s.lifetime = 6.582e-25/s.width # check units
        p = Mother(pname, mass, width)
        p.AddFS_brute(['e','v','v'], PDG.Ttoevv)
        p.AddFS_brute(['m','v','v'], PDG.Ttomvv)
        p.AddFS_brute(['Th','v'], PDG.Ttovh)
        particles[pname] = p
        geneaLog.Add(genea=genea, S='info', txt='AddSMtoMothers for %s' %(pname))
        geneaLog.Add(genea=genea, S='info', txt=p.ShowDecays(ret=1))   # text list sent to geneaLog
    #

    # ### Then the 1:1 rewriting
    if SMrewrite:
        rewrite =  [ ['u','q'], ['d','q'], ['s','q'], ['c','q'], ['g','q'], ['b','q'], ['e','l'], ['m','l'], ['ve','v'], ['vm','v'], ['vT','v'] ]


        ##rewrite += [ ['l','L'], ['T','L'] ]  ## ifTfundamental
        rewrite += [ ['l','L'], ['Th','L'] ]  ## ifTfundamental
        #rewrite += [ ['l',''], ['Th',''] ]  # 2014-02-03   # let l,Th -> L (if desired)

        rewrite += [ ['L','x'] ]
        rewrite += [ ['v','x'], ['y','x'], ['q','x'] ]
        rewrite += [ ['~N1','x'] ]  # 2014-01-20  Fixed: including this lead to inclusion of oneleg h->N1N1 also when N1 not in p2keep [and fixed the onelegs which included h]

        # 'x' is just 'anything' ; somehow there are problems if L,v,y,v,~N1 decay to '', but if they decay to 'x' which in turn goes to '', it seems to work (to be further tested)

        rewrite += [ ['x',''] ]
        
        for pname,pnamegeneric in rewrite:
            if pname in particles:
                p = particles[pname]
                geneaLog.Add(genea=genea, S='warn', txt='AddSMtoMothers  %s already exists, so using this (preserve e.g. mass:%.3f, width:%.3f), but add new decay (was N()=%i)' %(pname, p.mass, p.width, p.N()))
                p.ShowDecays()
            else:
                p = Mother(pname)
                particles[pname] = p
                
            if pnamegeneric == '':
                p.AddFS_brute(['x'], 1.)  # hack
            else:
                p.AddFS_brute([pnamegeneric], 1.)

            geneaLog.Add(genea=genea, S='debug', txt=p.ShowDecays(ret=1, what=['RELacc','REL']))

    # does not return anything; the change is done in the input dict: particles
    # end



# ##########
def SumMothers(ListOfMothers, name='sum', particles2keep=[]):
    sumMothers = Mother(pname=name, width=1., particles2keep=particles2keep)  # a bit non-intuitive to use width=1.
    SUMwidth = 0
    for mother in ListOfMothers: 
        # add the cascades of this subproc with appropriate xsec scaling:
        # (internal BRs should add up to the subprocces cross-section)
        sumMothers.AddMother(mother)  # Note: the width is passed into the FSlist
        SUMwidth += mother.width  # This preserves the total width in the sum
        
    sumMothers.NormaliseDecays(totBRsum=SUMwidth)  # this normalises the BRs in FSlist and puts the sum in width
    sumMothers.JoinIdenticalFSs()
    sumMothers.SortDecaysWithBR()
    
    return sumMothers


# ##########
def MotherfromMultiplyingTwoFSlists(fsList1, fsList2, name='2leg', sc=1., particles2keep=[]):
    mother = Mother(pname=name, particles2keep=particles2keep)
    # print "DEBUG::MotherfromMultiplyingTwoFSlists  len(list1)=%i  len(list2)=%i" %(len(fsList1),len(fsList2))
    mother.MultiplyTwoFSlists(fsList1, fsList2, sc=sc, opt=['join'])
    return mother


# ##########
class Mother: 
    # #####
    def __init__(s, pname='noname', mass=0., width=0., FSlist=[], particles2keep=[], VB=1, geneaLog=0, genea=''):
        if geneaLog: s.geneaLog = geneaLog
        else: s.geneaLog = GeneaLog(optD={'maxdepth':5, 'maxwidth1':12})
        genea = s.geneaLog.Genea(genea)
        
        s.name = pname
        s.mass = mass
        s.width = width  # maybe rename to sc or scale ... that could be used for both BR and xsec

        #s.totalwidth = 0.    # <-- This should guard the theoretical total ; so taking the difference with width gives estimate of loss (from finite BRcut)

        s.FSlist = list(FSlist)     # list of FS objects, each of which is a full(y developed) decay channel

        if particles2keep == []: s.particles2keep = ['l','Th','b','q','y','v','~N1']
        else: s.particles2keep = list(particles2keep)

        s.VB = VB
        

    # #####
    def N(s):
        return len(s.FSlist)
    
    # #####
    def ShowDecays1line(s):
        txt = ''
        for fs in s.FSlist: txt += ' (%s),' %(fs.ShowBrute())
        return txt.rstrip(',').strip()

    # #####
    def SortDecaysWithBR(s):
        s.FSlist = sorted(s.FSlist, key=lambda x: x.br, reverse=True)

    # #####
    def SortDecaysWithChildren(s):
        s.FSlist = sorted(s.FSlist, key=lambda x: x.children)

        
    # #####
    def AddFS_brute(s, childrenlist, br):  # add *one* fs [same as the below]
        fs = FS(childrenlist, br)
        s.AddFS(fs)


    # #####
    def AddFS(s, fs, sc=1.): # add *one* fs 
        added = False
        for existingdecay in s.FSlist:
            if existingdecay.children == fs.children:   # dict comparison
                existingdecay.br += fs.br * sc
                added = True
                break
        if not added: 
            s.FSlist.append(fs.Copy())   # over-safe? 
            s.FSlist[-1].br *= sc  # usually no scaling
            

    # #####
    #def AddFSList(s, bbFinalstateList, sc=1.):
    #    for fs in bbFinalstateList.FSlist: 
    #        s.AddFS(fs, sc=fs.br * bbFinalstateList.width)
            
    # #####
    def AddMother(s, mother, dropwidth=0):  #, sc=1.):   # Hmm Note: 
        for fs in mother.FSlist:
            if dropwidth: zsc = 1.
            else: zsc = mother.width
            s.AddFS(fs, sc=zsc) # * sc)   # Hmm: if BBFL.width = 0 (default), these will be added with zero BR
            

    # #####
    def Copy(s, name='copyname', sc=1.):
        #print 'HERE1, length = %i' %(len(s.FSlist))
        copy = Mother(name, mass=s.mass, width=s.width)
        #copy.AddFSList(s) #.FSlist, sc=sc)
        #print 'HERE2, length = %i' %(len(s.FSlist))
        copy.FSlist = s.CopyFSlist()
        #print 'HERE3, length = %i,  types are  %s and %s' %(len(s.FSlist), type(copy), type(copy.FSlist))
        return copy


    # #####
    def CopyFSlist(s):
        fsList = []
        for fs in s.FSlist:
            fsList.append(fs.Copy())
        
        return fsList


    # #####
    def MultiplyTwoFSlists(s, MotherLeg1, MotherLeg2, sc=1., opt=['join']):  # Typical when 1leg x 1leg -> 2leg ; sc would be xsec
        # So this should give sum of brs equal to 1 (or below depending on BRcuts) and have the xsec in s.width
        s.width = sc 
        for iD1 in range(len(MotherLeg1)): 
            fsLeg1 = MotherLeg1[iD1]

            for iD2 in range(len(MotherLeg2)):
                fsLeg2 = MotherLeg2[iD2].Copy()

                fsLeg2.MultWithFS(fsLeg1)

                s.AddFS(fsLeg2)

        if 'join' in opt: s.JoinIdenticalFSs()
    


    # #####
    def JoinIdenticalFSs(s):
        s.SortDecaysWithChildren()  # Note, a sort is done
        nStart = s.N()
        for iD in range(nStart-1, 0, -1):
            if s.FSlist[iD].children == s.FSlist[iD-1].children:
                s.FSlist[iD-1].br += s.FSlist[iD].br
                s.FSlist.pop(iD)

    # #####
    def RemoveFinalStatesWithLowBR(s, BRcut, opt=[]): 
        if 'sort' in opt: s.SortDecaysWithBR()
        
        iD_firstwithtoosmallBR = -1
        for iD in range(len(s.FSlist)): 
            if s.FSlist[iD].br < BRcut: 
                iD_firstwithtoosmallBR = iD
                break
        if iD_firstwithtoosmallBR > -1:
            s.FSlist = s.FSlist[:iD_firstwithtoosmallBR]
        
        

    # #####
    def NormaliseDecays(s, totBRsum=0.):   # Typically used after 1leg1leg->2leg
        # Loop #1 to find of BRsum
        BRsum = 0
        for fs in s.FSlist:
            BRsum += fs.br  # use if not input (then will get acc=1, which is not correct)
            #if s.VB>2: print 'DEBUG::NormaliseDecays  %10.6f  %10.6f' %(BRsum, fs.br)
        # Loop #2 to normalise to BRsum=1
        if totBRsum > 0: BRsum = totBRsum  # NB: input (new)
        
        for fs in s.FSlist:
            fs.br /= BRsum
        # Adjust overall scale (typically = cross-section)
        s.width *= BRsum  # old
        if totBRsum > 0: s.width = totBRsum  # new



    # #####
    def ShowSumBR(s):
        sumBR = 0.
        for fs in s.FSlist: sumBR += fs.br
        return sumBR
        

    # #####
    def ShowDecays(s, what=['ABSacc','RELacc','ABS','REL'], opt=['show', 'marklast'], which=[], pretext='', head=1, formats={}, txt_TO='  ->  ', ret=0):  # inspired by PhenoTool::Onelegs_Show
        
        whats = ['ABSacc','ABS','RELacc','REL']

        if which == []: which = s.particles2keep
        if ret:  # can also specify ret in this way
            if 'ret' not in opt: opt.append('ret')
            if 'show' in opt: opt.remove('show')
        
        # INIT
        if 'ABS' not in formats: formats['ABS'] = "%11.6f"  # %9.4f
        if 'ABSacc' not in formats: formats['ABSacc'] = formats['ABS']
        if 'REL' not in formats: formats['REL'] = "%9.7f"   # %7.5f
        if 'RELacc' not in formats: formats['RELacc'] = formats['REL']
        if 'parent' not in formats: formats['parent'] = "%-4s"

        N = s.N()
        outs = []

        RELacc = 0.
        ABSacc = 0.
        val = {}
        for zwhat in whats: val[zwhat] = 0.
        for idecay in range(N):
            fs = s.FSlist[idecay]
            #print 'DEBBUG: ', idecay, N, fs

            val['REL'] = fs.br
            val['ABS'] = fs.br * s.width
            val['RELacc'] += val['REL']
            val['ABSacc'] += val['ABS']

            # Build line
            out   = pretext
            zhead = pretext

            for zwhat in what:
                out   += formats[zwhat] %(val[zwhat]) + '  '
                #zformat = formats[zwhat].split('.')[0]+'s'
                zhead += (formats[zwhat].split('.')[0]+'s').replace('%','%') %(zwhat) + '  '


            if not 'noparent' in opt:  # (can opt not to show parent)
                out   += '  ' + formats['parent'] %(s.name) + txt_TO
                zhead += ('  %%-%is' %(max(len(s.name),int(formats['parent'].replace('%','').replace('-','').replace('s','')))) )  %('PAR') + txt_TO  # (to make a simple thing diffcult)
            

            zhead += '%s' %(fs.Show(mode='dict', which=which, head=1)) 
            out   += '%s' %(fs.Show(mode='dict', which=which))

            if idecay == N-1 and 'marklast' in opt: out += '   (last)'

            if idecay == 0 and head > 0 : outs.append(zhead)
            outs.append(out)


        if 'show' in opt:
            for out in outs: print(out)

        if 'ret' in opt: return outs

            
    # #####
    # ##### This is the problem routine, complex like hell
    def DevelopFinalstates(s, genea='', FSlist=[], particles=[], Mothers_Onelegs_existing={}, optD={}, particles2keep=[], VB=1):
        genea = s.geneaLog.Genea(genea)
        if particles2keep != []: s.particles2keep = list(particles2keep)

        if FSlist == []: FSlist = list(s.FSlist)  # use/make separate list
        s.FSlist = []

        BRmin1leg  = optD.get('BRmin1leg',1e-3)
        BRminDecay = optD.get('BRminDecay',1e-3)

        thisMother_Onelegs = s   # hack for now
        
        if 1: 

            # B Loop over initial decays [E.g.: N2 has decay list [N1+Z,N1+h,C1+W,C1+lv,C1+qq,..]; loop over these
            for iDecay in range(len(FSlist)):   # p is BBfinalstates (of a particle)
                fs = FSlist[iDecay]
                
                if fs.br < BRminDecay: continue
                

                # I think a decay box can be represented by a Mother (a BBFL)...
                # Should become BBFLs:
                #   sons_decays{}          -> Mother_ofson[]
                #   decay_combinations{}   -> Mother_accumulated_ofson[]
                #   Mothers_Onelegs[]


                # The start object is one FS: fs
                # Will now create a decay box for each particle (son) in this decay: Mothers_ofson[0], Mothers_ofson[1], ... 
                #   [a decay box = a fully expanded decay list for each particle]
                # The N decay boxes are then combined/multipled to the the full list of fully expanded decays
                # Ex: for a decay 
                #   Mother_accumulated_ofson[0] = Mothers_ofson[0]
                #   Mother_accumulated_ofson[1] = Mothers_ofson[1] x Mothers_ofson[0]
                #   Mother_accumulated_ofson[2] = Mothers_ofson[2] x Mothers_ofson[1]
                #   ...
                #   Mother_accumulated_ofson[-1]  will hold the final result
                # So the final list is Mother_accumulated_ofson[N-1] : this list will then be added to Mothers_Onelegs[pname_initiator]


                Mothers_accumulated_ofson = {}  # Note: keys are 0,1,..
                Mothers_ofson = {}              # Note: keys are 0,1,.. ; values are lists of Fss (-> #BBFL#)
                pname_sons = fs.ChildrenList()  # ex: ['~N1','l','l'] or ['~N1','c','s'],  maybe later:['~N1','q','q']
                    
                #if VB>2: print '\n\n', 80*'-', '\nDevelopFinalstates::B::loopoverinitialdecays    %-3s    iDecay:%i  %s  BR=%7.5f' %(s.name, iDecay, pname_sons, fs.br)
                s.geneaLog.Add(genea=genea, S='dump', txt='loopoverinitialdecays    %-3s    iDecay:%i  %s  BR=%7.5f' %(s.name, iDecay, pname_sons, fs.br), space=[1,0])

                # C Loop over sons and build a fully expanded decay list for each son
                for ison in range(len(pname_sons)):

                    # 1) First create the BBFL for this son: Mothers_ofson[ison]  # Note: using index as key
                    pname_son = pname_sons[ison]

                    if pname_son in s.particles2keep:  # this is required
                        ###Mothers_ofson[ison] = [FS({pname_son:1}, br=1.)] # works? #BBFL#
                        s.geneaLog.Add(genea=genea, S='dump', txt='C::loopoversons:  %-3s is in s.particles2keep. Creating Mother(%s)' %(pname_son, pname_son))
                        Mothers_ofson[ison] = Mother(pname=pname_son, particles2keep=s.particles2keep, geneaLog=s.geneaLog, genea=genea);
                        Mothers_ofson[ison].AddFS_brute({pname_son:1}, br=1.)
                        
                    elif pname_son in Mothers_Onelegs_existing:
                        Mothers_ofson[ison] = Mothers_Onelegs_existing[pname_son]  # list of Fss #BBFL#
                        s.geneaLog.Add(genea=genea, S='dump', txt='loopoversons:  %-3s now being picked from Mothers_Onelegs:' %(pname_son))
                        s.geneaLog.Add(genea=genea, S='dump', txt='   mass:%3.f  width:%.3f ' %(Mothers_ofson[ison].mass, Mothers_ofson[ison].width))
                            
                    else:
                        s.geneaLog.Add(genea=genea, S='debug', txt='C::loopoversons:  %-3s (since in decay of %s) now being created in ExpandDecay:' %(pname_son, s.name))
                        s.geneaLog.Add(genea=genea, S='dump', txt=['For comparison show the entry in particles dict (before):', particles[pname_son].ShowDecays(ret=1)])  # at this stage particles[W] is still fine
                        # So the conclusion must be that ExpandDecay alters/corrupts particles
                        Mothers_ofson[ison] = s.ExpandDecay(pname=pname_son, particles=particles, Mothers_Onelegs_existing=Mothers_Onelegs_existing, genea=genea)  # list of Fss ; Might want to simplify/join at this stage..
                        Mothers_Onelegs_existing[pname_son] = Mothers_ofson[ison]  # 2014-02-03  Is this the bugfix?  Yes. Though I don't understand why it would fail earlier. I would think it should just redo and redo, but be correct each time. 

                        s.geneaLog.Add(genea=genea, S='dump', txt='.. and there it (%s) was created: ' %(pname_son))
                        s.geneaLog.Add(genea=genea, S='dump', txt=Mothers_Onelegs_existing[pname_son].ShowDecays(ret=1))
                      
                        s.geneaLog.Add(genea=genea, S='dump', txt=['For comparison show the entry in particles dict:', particles[pname_son].ShowDecays(ret=1)])  # at this stage particles[W] is ruptured 
                        # this seems to work at least for sons like 'd','e','~N1' which go to 'q','l',[] ..
                        ### BBFL: check


                    if VB>3: 
                        print('  loopoversons (done)  ison=%i :  (len=%i)  %s' %(ison, Mothers_ofson[ison].N(), Mothers_ofson[ison]))
                        for z in Mothers_ofson[ison].FSlist:
                            print('                 children %s' %(z.children))
 
                        
                    # 2) Then add this son's "decay box (list of decays to end particles)" to all entries of the Mothers_accumulated_ofson
                    if ison == 0:
                        # if first
                        ###Mothers_accumulated_ofson[ison] = list(Mothers_ofson[ison])  # pointer to list of Fss
                        Mothers_accumulated_ofson[ison] = Mothers_ofson[ison].Copy()  #BBFL#  (pointer/object?)
                        s.geneaLog.Add(genea=genea, S='dump', txt='isonA: %i  len(Mothers_accumulated_ofson[%i])=%i' %(ison, ison, Mothers_accumulated_ofson[ison].N()))
                        s.geneaLog.Add(genea=genea, S='dump', txt='decaybox0  ison=0')
                        s.geneaLog.Add(genea=genea, S='dump', txt=Mothers_accumulated_ofson[ison].ShowDecays(ret=1))

                    else:    
                        
                        ###Mothers_accumulated_ofson[ison] = []
                        Mothers_accumulated_ofson[ison] = Mother(particles2keep=s.particles2keep, geneaLog=s.geneaLog, genea=genea)  #BBFL# 
                        s.geneaLog.Add(genea=genea, S='dump', txt='    ison=%i  %s' %(ison, Mothers_accumulated_ofson))
                        s.geneaLog.Add(genea=genea, S='dump', txt='    isonB: %i  len(Mothers_accumulated_ofson[%i])=%i' %(ison, ison-1, Mothers_accumulated_ofson[ison-1].N()))


                        if Mothers_accumulated_ofson[ison-1].N() == 0 and Mothers_ofson[ison].N() == 0:  # IF BOTH ACCU AND CURRENT PARTICLE HAVE EMPTY FL
                            s.geneaLog.Add(genea=genea, S='dump', txt='    decayboxA  both zero: pass')
                            pass # ... keep too? Maybe ... ex: ['~N1','v','l'] with only l to be kept ... TEST
                        
                        elif Mothers_accumulated_ofson[ison-1].N() == 0: # IF ACCU EMPTY
                            ###for fs in Mothers_ofson[ison]:
                            ###    Mothers_accumulated_ofson[ison].append(FS(fs.children, fs.br))  #*fs.br))
                            Mothers_accumulated_ofson[ison] = Mothers_ofson[ison].Copy()  # #BBFL#  might also do without the copy
                            if VB>3: print('Mothers_ofson[%i]:' %(ison)); Mothers_ofson[ison].ShowDecays()
                            if VB>3: print('Mothers_accumulated_ofson[%i]:' %(ison)); Mothers_accumulated_ofson[ison].ShowDecays()
                            s.geneaLog.Add(genea=genea, S='dump', txt='-- the above should not be zero')
                            
                        elif Mothers_ofson[ison].N() == 0: # IF CURRENT PARTICLE HAS EMPTY FL
                            ###for existing_fs in Mothers_accumulated_ofson[ison-1]:
                            ###    Mothers_accumulated_ofson[ison].append(existing_fs.children, existing_fs.br) #*fs.br)
                            Mothers_accumulated_ofson[ison] = Mothers_accumulated_ofson[ison-1].Copy()  # #BBFL# might also do without the copy 
                            s.geneaLog.Add(genea=genea, S='dump', txt='    decayboxC (zero current)')

                        else: 
                            s.geneaLog.Add(genea=genea, S='debug', txt='     BeforeMult: ison=%i   (ndecay=%i)' %(ison, Mothers_accumulated_ofson[ison].N()))
                            Mothers_accumulated_ofson[ison] = MotherfromMultiplyingTwoFSlists(Mothers_accumulated_ofson[ison-1].FSlist, Mothers_ofson[ison].FSlist, particles2keep=s.particles2keep, name=s.name)  #BBFL#  #correct with s.name?
                            s.geneaLog.Add(genea=genea, S='debug', txt='     AfterMult:  ison=%i   (ndecay=%i)' %(ison, Mothers_accumulated_ofson[ison].N()))

                        #


                    # done with Mothers_accumulated_ofson[ison] = Mothers_accumulated_of_son[ison-1] x Mothers_ofson[ison] to BBFLs
                    
                
                # Should now have complete list of fully-expanded Fss for the given starter decay for the given particle
                # add them to the decay list for this given particle
                #print 70*'=',"DDD len(Mothers_accumulated_ofson) = %i , type: %s" %(len(Mothers_accumulated_ofson), type(Mothers_accumulated_ofson))

                # Note: for the end-states, ~N1,l,Th,y,v,q, the loop above is not done at all since pname_sons = []
                
                nson = len(Mothers_accumulated_ofson)
                s.geneaLog.Add(genea=genea, S='dump', txt='nson=%i, name=%s' %(nson, s.name))
                #if VB>2: print 10*'=',"DDD the given decay  (%i  %s -> %s) was expanded into len(Mothers_accumulated_ofson) = %i decays, type: %s" %(iDecay, s.name, pname_sons, Mothers_accumulated_ofson[nson-1].N(), type(Mothers_accumulated_ofson[nson-1]))


                if VB>2:
                    for ison in range(nson):
                        #print 25*'+'
                        s.geneaLog.Add(genea=genea, S='dump', txt=['Mothers_ofson[%i]: %s' %(ison, Mothers_ofson[ison].name), Mothers_ofson[ison].ShowDecays(ret=1), 'Mothers_accumulated_ofson[%i]' %(ison), Mothers_accumulated_ofson[ison].ShowDecays(ret=1)])
                        #print 25*'~'
                
                ###for zfs in Mothers_accumulated_ofson[len(Mothers_accumulated_ofson)-1]: zfs.br *= fs.br  # to include the initial br
                s.geneaLog.Add(genea=genea, S='dump', txt='DDEBUG: pname_sons=%s  len(Mothers_accumulated_ofson)=%i  keys=%s' %(pname_sons, len(Mothers_accumulated_ofson), list(Mothers_accumulated_ofson.keys())))

                # Correction of BRs
                s.geneaLog.Add(genea=genea, S='dump', txt='Now to correction of BRs  : %s' %(s.name))
                s.geneaLog.Add(genea=genea, S='dump', txt=s.ShowDecays(ret=1))  # this shows nothing (decay of initial sparticle, e.g. ~C1)
                if len(pname_sons) > 0: # 2014-01-20 put inside if-statement
                    for zfs in Mothers_accumulated_ofson[len(Mothers_accumulated_ofson)-1].FSlist:
                        s.geneaLog.Add(genea=genea, S='debug', txt='  correcting ...  %7.5f x %7.5f[overallforthischain] = %7.5f    (%s)' %(zfs.br, fs.br, zfs.br*fs.br, zfs.children))
                        zfs.br *= fs.br  #BBFL# to include the initial br
                    # ##Mothers_Onelegs[pname_initiator] += Mothers_accumulated_ofson[len(Mothers_accumulated_ofson)-1]  ### #BBFL#
                    #if VB>2: print 'there : '
                    s.geneaLog.Add(genea=genea, S='debug', txt=s.ShowDecays(ret=1)) # this shows nothing, again nothing (s is e.g. ~C1)
                    thisMother_Onelegs.AddMother(Mothers_accumulated_ofson[len(Mothers_accumulated_ofson)-1], dropwidth=1)  #NB width
                    s.geneaLog.Add(genea=genea, S='debug', txt='Show last accumulated:  ')
                    #geneaTxt = s.geneaLog.Add(genea=genea, S='debug', ret=1, store=0, txt='')
                    s.geneaLog.Add(genea=genea, S='debug', txt=Mothers_accumulated_ofson[len(Mothers_accumulated_ofson)-1].ShowDecays(ret=1))  # HERE BR=1
                    s.geneaLog.Add(genea=genea, S='debug', txt='Show thisMother_Onelegs (%s) in the building:' %(thisMother_Onelegs.name))
                    s.geneaLog.Add(genea=genea, S='debug', txt=thisMother_Onelegs.ShowDecays(ret=1))  # HERE BR=0 (absolute BR is zero)

                s.geneaLog.Add(genea=genea, S='dump', txt=['before delete : ', s.ShowDecays(ret=1)])
                del Mothers_accumulated_ofson, Mothers_ofson  # deeper delete preferable/needed? 
                s.geneaLog.Add(genea=genea, S='dump', txt=['after delete : ', s.ShowDecays(ret=1)])
                #if VB>2: print 'after delete :  '; s.ShowDecays()
                # ===== end of treatment for this given decay, the result is a BBFL (possibly several decays) which is added to Mothers_Onelegs[pname_initiator]


            # ===== end of loop B over unexpanded decays

            s.geneaLog.Add(genea=genea, S='dump', txt=s.ShowDecays(ret=1))
            thisMother_Onelegs.JoinIdenticalFSs()
            thisMother_Onelegs.SortDecaysWithBR()
            thisMother_Onelegs.RemoveFinalStatesWithLowBR(BRcut=BRmin1leg)

            # The treatment is now complete
            
        #del FSlist  # would this remove the particles[pname_initiator].FSlist ? 
        # end



    # #####       
    def ExpandDecay(s, pname, particles, Mothers_Onelegs_existing, VB=1, genea=''):  # this is just a helper routine .. maybe should be elsewhere
        # THIS IS ONLY USED WHEN PARTICLES FOR WHICH THERE IS NO 1LEG ARE FOUND IN SOME DECAY 
        # if pname in s.particles2keep: return pname
        # returns a 1legsDict object (which currently is just a decay list)
        # BBFL: decays -> Mother_expanded
        # 2014-08-29: FIX1a-c : Now the nonconservation of BR is fixed (since particles dict not used&changed)

        genea = s.geneaLog.Genea(genea)
        
        decays = []

        if 1 and pname == 'x':  # 2014-08-29  FIX1c [needed with FIX1b to avoid infinite loop for 'x']
            zFS = FS({},1.)
            resBBFSList = Mother(pname, width=1., FSlist=[zFS], particles2keep=s.particles2keep, geneaLog=s.geneaLog, genea=genea)
            return resBBFSList


        if pname in particles:
            s.geneaLog.Add(genea=genea, S='dump', txt='  ExpandDecay begin %s  len(decays@start) = %i' %(pname, len(decays)))
            #decays += particles[pname].FSlist   # list of FS objects (.br and .children{}) 
            for bbfs in particles[pname].FSlist: decays.append(bbfs.Copy())   # FIX1b now do not particles dict (since deep copy)
            s.geneaLog.Add(genea=genea, S='dump', txt='  ExpandDecay begin %s  from particles : len(decays@start) = %i' %(pname, len(decays)))
        else:
            decays = [ FS(children=[pname], br=1.) ]
            s.geneaLog.Add(genea=genea, S='dump', txt='  ExpandDecay begin %s  from scratch   : len(decays@start) = %i' %(pname, len(decays)))

        # Then expand
        KeepOn = True
        s.geneaLog.Add(genea=genea, S='dump', txt='  ExpandDecay begin %s  len(decays@start) = %i' %(pname, len(decays)))
        while KeepOn: 
            KeepOn = False
            iTL = 0

            while iTL < len(decays):  # [will loop forever if 'x' is not done first in 1leg]
                # 
                changeinthisdecay = False

                childrenlist = decays[iTL].ChildrenList()

                s.geneaLog.Add(genea=genea, S='dump', txt='  iTL: %i  childrenlist: %s' %(iTL, childrenlist))
                
                # Traverse one specific decay channel
                for pSM in childrenlist:   # SMcombExpand: 
                    if pSM in s.particles2keep:
                        s.geneaLog.Add(genea=genea, S='dump', txt='  pSM %s not to expand' %(pSM))
                        continue

                    #if pSM in decays[iTL].children:  # satisfied by construction (could rewrite since this is now the case)
                   

                    # If here, the given child is to be expanded
                    KeepOn = True
                    changeinthisdecay = True
                    
                    thisfs = decays[iTL]   # pick out the given decay
                    zBRsub0 = thisfs.br  # this is e.g. 

                    s.geneaLog.Add(genea=genea, S='dump', txt='par:%s  changeinthisdecay: len=%i  iTL=%i  %s' %(pname, len(decays), iTL, decays[iTL].Show()))

                    s.geneaLog.Add(genea=genea, S='dump', txt='par:%s  len=%i  iTL=%i  expanding %s  1->%i  : %s -> %s ' %(pname, len(decays), iTL, pSM,  particles[pSM].N(), pSM, particles[pSM].ShowDecays1line()))
                    

                    ### thisfs = Mother_expanded.FSlist[iTL]  #BBFL#
                    # The first replacement decay is inserted in the same position, then need not pop the initial decay, just change it
                    thisfs.children[pSM] -= 1    # change it: pop away (one instance of the) expanding particle
                    if thisfs.children[pSM] == 0: thisfs.children.pop(pSM)
                    brOrig = thisfs.br
                    childrenOrig = dict(thisfs.children)  # original children, except the expanded one
                    
                    if particles[pSM].N() == 0:
                        s.geneaLog.Add(genea=genea, S='warn', txt='Warning::ExpandDecay  Hack  Empty decay for %s (If LSP, I think it is ok..)' %(pSM))
                        del childrenOrig
                        continue
                    
                    if 1 and pSM in Mothers_Onelegs_existing:   
                        replaceFSlist = Mothers_Onelegs_existing[pSM].FSlist  # 2014-08-29: FIX1a: this fixes the BR:1->1.07 for h (and more)
                    else: 
                        replaceFSlist = particles[pSM].FSlist


                    thisfs.MultWithFS(replaceFSlist[0]) # NOT CORRECT #1  HMM Why not? 

                    zBRsubs = thisfs.br
                    s.geneaLog.Add(genea=genea, S='dump', txt='par:%s  BRtest  old: %9.6f   newsum: %9.6f      thisnew %i/%i %9.6f ' %(pname, zBRsub0, zBRsubs, 0, len(particles[pSM].FSlist), thisfs.br))



                    # Further replacement decays are inserted at the end, using the Orig values
                    ifs = 0
                    for zfs in replaceFSlist[1:]:  # if only 1 entry, then empty, also empty if 0 entry
                        ifs += 1 
                        newdecay = FS(childrenOrig, brOrig)
                        newdecay.MultWithFS(zfs) # NOT CORRECT #2  (Why not?)
                        decays.append(newdecay)

                        zBRsubs += newdecay.br
                        s.geneaLog.Add(genea=genea, S='dump', txt='par:%s  BRtest  old: %9.6f   newsum: %9.6f      thisnew %i/%i %9.6f ' %(pname, zBRsub0, zBRsubs, ifs, len(replaceFSlist), newdecay.br))


                    # Test if BR is conserved (used to not be due to unintended changes and use of particles dict ; should be fine now)
                    if (zBRsub0+zBRsubs) != 0 and abs(zBRsub0-zBRsubs) / (zBRsub0+zBRsubs) > 1e-5: 
                        s.geneaLog.Add(genea=genea, S='warn', txt='par:%s  BR not conserved in expansion %3i/%i   %8.6f -> %8.6f   : %s -> ' %(pname, iTL, len(decays), zBRsub0, zBRsubs, pSM)) 
                    else: 
                        s.geneaLog.Add(genea=genea, S='dump', txt='match: %8.6f' %((zBRsub0-zBRsubs) / (zBRsub0+zBRsubs)))

                    del childrenOrig
                    

                #
                if not changeinthisdecay: iTL += 1
            #


        s.geneaLog.Add(genea=genea, S='dump', txt='  ExpandDecay  %-3s end:  len(decays@end) = %i' %(pname, len(decays)))
        for ibb in range(len(decays)): s.geneaLog.Add(genea=genea, S='dump', txt='  ExpandDecay  %-3s end:  %i  %s' %(pname, ibb, decays[ibb].ChildrenList()))

        
        ###return decays
        resBBFSList = Mother(pname, width=1., FSlist=decays, particles2keep=s.particles2keep, geneaLog=s.geneaLog, genea=genea)  # Minimal #BBFL#
        resBBFSList.ShowDecays()
        s.geneaLog.Add(genea=genea, S='dump', txt='n(decays before joinidentical) = %i     SumBR=%9.6f' %(resBBFSList.N(), resBBFSList.ShowSumBR()))
        resBBFSList.JoinIdenticalFSs()  # 2014-08-28: previously was not done for h,W,Z (moved inside
        s.geneaLog.Add(genea=genea, S='dump', txt='n(decays after joinidentical)  = %i     SumBR=%9.6f' %(resBBFSList.N(), resBBFSList.ShowSumBR()))

        return resBBFSList


    # #####


# ##########
class FS: # Is just one decay channel
    # does not keep track of where things come from
    def __init__(s,children={}, br=1.):  # input children can be dict or list; both treated in AddToChildren
        s.br = br
        s.children = {}
        # ex: s.children = {'~N1':1, 'l':2, 'b':2, 'Th':1}
        s.AddChildren(children)


    # #####
    def Copy(s):
        thechildren = dict(s.children)
        acopy = FS(thechildren, s.br)
        return acopy 
        

    # #####
    def AddChildren(s,children, br=1.):  # initialises and adds / sets the children
        s.br *= br
        if type(children) is dict: s.children = dict(children)
        if type(children) is list:
            s.children = {}
            for p in children:
                if p not in s.children: s.children[p] = 0
                s.children[p] += 1

    # ##### 
    def MultWithFS(s, fs, sc=1.):
        s.AddToChildren(fs.ChildrenList(), fs.br * sc)

        
    # #####
    def AddToChildren(s,childrenlist, brreduction=1.):  # add to existing
        if type(childrenlist) is dict: raise Exception("Fatal::FS::AddToChildren  childrenlist must be a list, not a dict")
        s.br *= brreduction
        for pname in childrenlist:  
            if pname not in s.children: s.children[pname] = 0
            s.children[pname] += 1

    # #####
    def ChildrenList(s):
        childrenlist = []
        for zname in s.children:
            for i in range(s.children[zname]): childrenlist.append(zname.replace('~','00000~')) # hack: SUSY children first
        childrenlist.sort()
        for i in range(len(childrenlist)):
            childrenlist[i] = childrenlist[i].replace('00000~','~') # hack: SUSY children first
        return childrenlist

    # #####
    def ShowBrute(s):
        zlist = []
        for zpart in s.children:
            for ziter in range(s.children[zpart]): zlist.append(zpart)
        zlist.sort()
        zout = ''
        for zpart in zlist: zout += ' ' + zpart
        return zout.strip()
    
        
    # #####
    def Show(s, mode='dict', which=['~N1','b','Th', 'l','q','y','v'], nmax=4, nwid=5, zero='.', head=0):
        # construct correctly ordered plist
        chlist = list(which)
        chlist2 = []
        for ch in s.children:
            if ch not in which: chlist2.append(ch)


        out = ''
        headT = ''
        # Mode nameddict
        if mode in ['nameddict', 'dict']:
            format = ' %%-%is' %(nwid)  # e.g. "%4s"
            for ch in chlist:
                zout = '%s' %(str(s.children.get(ch,zero)))
                if mode in ['nameddict']: zout += '%s' %(ch)
                out += format %(zout)
                headT += format %(ch)
            #for i in range(nmax-len(chlist)): outs += format %()

            out2 = ''
            for ch in chlist2:
                out2 += '%s:%i ' %(ch, s.children[ch])
            if out2: out += '  NBEXTRA: %s' %(out2.rstrip(' '))   # this is to clearly show that there are finalstate particles not anticipated in the input setup
            
            if head: return headT
            else: return out

        if mode in ['flat']:
            pass
        

# #############################################################
# #############################################################
# #############################################################






# #############################################################
# #############################################################
# #############################################################

slhablockcode2name = {
    #'MODSEL': {
    #3:'',
    #4:'',
    #5:'',
    #6:''
    #} , 

    'SPINFO':{
    1:'CALC',
    2:'CALCver'
    } ,
    

    'SMINPUTS':{
    1:'al_emInv(Z)',
    2:'G_F',      #  [GeV^-2]
    3:'al_s(Z)',  # alpha_S(M_Z)^MSbar
    4:'Z',        # pole
    5:'b',        # mb(mb)^MSbar
    6:'t',        # mt pole mass
    7:'T',        # mtau pole mass
    } ,

    'MINPAR': { 
    1:'m0',      # this is mSUGRA ; for GMSB/AMSB the parameters have different meaning
    2:'m1/2',
    #3:['TB','TB(Z)'],      # (common for all models)
    3:['TB(Z)'],      # (common for all models)
    4:'sign(mu)',
    5:'A'

    } , 

    #'AD': {1:'Ad', 2:'As', 3:'Ab'}, 
    #'AU': {1:'Au', 2:'Ac', 3:'At'}, 
    #'AE': {1:'Ae', 2:'Am', 3:'AT'}, 
    #'YD': {1:'Yd', 2:'Ys', 3:'Yb'}, 
    #'YU': {1:'Yu', 2:'Yc', 3:'Yt'}, 
    #'YE': {1:'Ye', 2:'Ym', 3:'YT'}, 

    # The pyslha structure is crazy: 
    'AD': {1:{1:'Ad'}, 2:{2:'As'}, 3:{3:'Ab'}}, 
    'AU': {1:{1:'Au'}, 2:{2:'Ac'}, 3:{3:'At'}}, 
    'AE': {1:{1:'Ae'}, 2:{2:'Am'}, 3:{3:'AT'}}, 
    'YD': {1:{1:'Yd'}, 2:{2:'Ys'}, 3:{3:'Yb'}}, 
    'YU': {1:{1:'Yu'}, 2:{2:'Yc'}, 3:{3:'Yt'}}, 
    'YE': {1:{1:'Ye'}, 2:{2:'Ym'}, 3:{3:'YT'}}, 

    'ALPHA':'h_alpha',

    #'MASS': use SlhaTools.pdg2nam(pdgid) directly


    'HMIX': {
    #1:['mu(Q)','mu','MU'], 
    1:'mu(Q)', 
    2:'TB(Q)',
    3:'VEV(Q)',
    4:'mA^2(Q)'
    } , 


    'MASS': pdg2nam,
    

    'EXTPAR': {
    0:['Minput'],  # was (wrongly) EWSB [only true if Minput=EWSB, which, true, is often specified]
    1:['M1','M1(MX)'],
    2:['M2','M2(MX)'],
    3:['M3','M3(MX)'],

    11:['At','At(MX)'],  # duplicates of part of AD etc above; but many-to-one is no problem
    12:['Ab','Ab(MX)'],
    13:['AT','AT(MX)'],
    14:['Au','Au(MX)'],
    15:['Ad','Ad(MX)'],
    16:['Ae','Ae(MX)'],
    
    21:['mH1^2(MX)'],  
    22:['mH2^2(MX)'],  
    23:['mu','mu(MX)'], 
    24:['mA^2(MX)'],  
    25:['TB(MX)'],
    26:['mA','mA(pole)'],   # (pole) [typically given instead of mA^2(MX)]
    27:['mH+','mH+(pole)'],  # (pole) can be given instead of mA
    
    31:['L1','L1(MX)'],
    32:['L2','L2(MX)'],
    33:['L3','L3(MX)'],
    34:['eR','eR(MX)'],
    35:['mR','mR(MX)'],
    36:['TR','TR(MX)'],
    
    41:['Q1','Q1(MX)'],
    42:['Q2','Q2(MX)'],
    43:['Q3','Q3(MX)'],
    44:['uR','uR(MX)'],
    45:['cR','cR(MX)'],
    46:['tR','tR(MX)'],
    47:['dR','dR(MX)'],
    48:['sR','sR(MX)'],
    49:['bR','bR(MX)'],

    51:['Nmsg1','Nmsg1(MX)'],  # GMSBonly:  U(1) messenger index
    52:['Nmsg2','Nmsg2(MX)'],  # GMSBonly: SU(2) messenger index
    53:['Nmsg3','Nmsg3(MX)'],  # GMSBonly: SU(3) messenger index


    # further parameters are for NMSSM++
    } ,


    
    'MSOFT': {

    1:['M1(Q)'],
    2:['M2(Q)'],
    3:['M3(Q)'],

    #11:['At(Q)'],  # do the A-parameters exist in MSOFT? Not according to LH manual
    #12:['Ab(Q)'],  
    #13:['AT(Q)'],
    #14:['Au(Q)'],
    #15:['Ad(Q)'],
    #16:['Ae(Q)'],
    
    21:['mH1^2(Q)'], 
    22:['mH2^2(Q)'], 
    #23:['mu(Q)','mu'],   # these entries do not exist according to LH manual
    #24:['mA^2(Q)'],   # 
    #25:['TB(Q)'],
    #26:['mA','mA(pole)'],   # (pole)
    #27:['mH+','mH+(pole)'],  # (pole) can be given instead of mA
    
    31:['L1(Q)'],
    32:['L2(Q)'],
    33:['L3(Q)'],
    34:['eR(Q)'],
    35:['mR(Q)'],
    36:['TR(Q)'],
    
    41:['Q1(Q)'],
    42:['Q2(Q)'],
    43:['Q3(Q)'],
    44:['uR(Q)'],
    45:['cR(Q)'],
    46:['tR(Q)'],
    47:['dR(Q)'],
    48:['sR(Q)'],
    49:['bR(Q)'],

    #51:'Nmsg1(Q)',  # GMSBonly:  U(1) messenger index
    #52:'Nmsg2(Q)',  # GMSBonly: SU(2) messenger index
    #53:'Nmsg3(Q)',  # GMSBonly: SU(3) messenger index


    # further parameters are for NMSSM++
    } ,



    'GAUGE':{1:"g'(Q)", 2:'g(Q)', 3:'g3(Q)'} ,  

    'NMIX': {
    1:{1:'N11', 2:'N12', 3:'N13', 4:'N14'} , 
    2:{1:'N21', 2:'N22', 3:'N23', 4:'N24'} , 
    3:{1:'N31', 2:'N32', 3:'N33', 4:'N34'} , 
    4:{1:'N41', 2:'N42', 3:'N43', 4:'N44'}       
    } , 

    'UMIX': {1:{1:'U11',2:'U12'}, 2:{1:'U21',2:'U22'}} ,   # chargino

    'VMIX': {1:{1:'V11',2:'V12'}, 2:{1:'V21',2:'V22'}} ,   # chargino

    'SBOTMIX':{1:{1:'b11',2:'b12'}, 2:{1:'b21',2:'b22'}} , 

    'STOPMIX':{1:{1:'t11',2:'t12'}, 2:{1:'t21',2:'t22'}} , 
    
    'STAUMIX':{1:{1:'T11',2:'T12'}, 2:{1:'T21',2:'T22'}} , 

    
    }



# ##########
slhablocks = [  # order as given by SUSYHIT (for 'intuitive' printing)
    'DCINFO'
    ,'SPINFO'
    ,'MODSEL'
    ,'SMINPUTS'
    ,'MINPAR'
    ,'EXTPAR'
    ,'MASS'
    ,'NMIX'
    ,'UMIX'
    ,'VMIX'
    ,'STOPMIX'
    ,'SBOTMIX'
    ,'STAUMIX'
    ,'ALPHA'
    ,'HMIX'
    ,'GAUGE'
    ,'AU'
    ,'AD'
    ,'AE'
    ,'Yu'
    ,'Yd'
    ,'Ye'
    ,'MSOFT'
]


# ##########
def ShowFlatSlhaVars(show=1,ret=0): 
    D = slhablockcode2name
    blocks = sorted(D.keys())

    outs = ['  %-8s   %s' %('BLOCK','VAR')]
        
    for block in blocks:
        d1 = D[block]

        if not type(d1) is dict:
            outs.append('  %-8s   %s' %(block, d1))

        else: # dict
            keys1 = sorted(d1.keys())
            for key1 in keys1:
                d2 = d1[key1]
                
                if not type(d2) is dict:
                    outs.append('  %-8s   %s' %(block, d2))
                else: # dict[dict] 

                    keys2 = sorted(d2.keys())
                    for key2 in keys2:
                        d3 = d2[key2]

                        outs.append('  %-8s   %s' %(block, d3))

                    
    
    if show:
        for out in outs: print(out)

    if ret: return outs



        
    

# ##########
def slha2flatdict(SLHA=0, fn='', opt={}, VB=0):  # NOTE: name misleading: the flatdict is not returned, but is a in key 'flat'
    # returns in flat dict the content of slha (some quantities are left out)
    # returns also a (dict)list of the order of the variables
    # returns a list of infos/warnings

    # #####  an internal function for cleaner code
    def AddToResDict(res, block, var, val):
        if block == 'MASS': val = abs(val)  # hack to ensure positive neutralino-mass 
        if type(var) is list: vars = var # to allow several vars for same (e.g. mu,MU)
        else: vars = [var]
        for zvar in vars:
            res['flat'][zvar] = val
            res['vars'][block].append(zvar)
            res['var2block'][zvar] = block
    # #####


    # 1) LOAD (if needed)
    if fn:
        try:
            SLHA = pyslha.readSLHAFile(fn)
        except:
            print('Error::slha2flatdict  Problem reading slha file %s' %(fn))
            return {}
        
    if not SLHA:
        print('Warning::slha2flatdict need to give SLHA object (SLHA=) or filename (fn=). Returning empty dictionary')
        return {}



    # 2) INIT
    res = {'flat':{}, 'warn':[], 'vars':{}, 'var2block':{}}




    # 3) FILL FLAT DICT FROM SLHA OBJECT
    for block in SLHA[0]: 
        if VB>2: print('block: ', block)
        if block not in slhablockcode2name:
            res['warn'].append('Not in slhablockcode2name: block=%s' %(block))
            if VB>1: print('Warning::SlhaTools  ', res['warn'][-1])
            continue

        slha = SLHA[0][block].entries
        c2n = slhablockcode2name[block]
        
        res['vars'][block] = []

        # --- adhoc structure to pick up Q (scale) of each block ; some have Q=None, others have the same value (optimally)
        try:
            
            #print 'block: ', block, SLHA[0][block]
            zQ = SLHA[0][block].q
            #if zQ == None: print 'is none for block %s' %(block)
            res['flat']['Q_'+block] = zQ

            if 'Q' in res['flat'] and zQ != None:
                if res['flat']['Q'] == None:
                    res['flat']['Q'] = zQ
                elif res['flat']['Q'] != zQ:
                    res['warn'].append('Different Q replacing old:   %s has %s != %s already stored' %(block, zQ, res['flat']['Q']))
                    res['flat']['Q'] = zQ
                else:
                    pass  # new and old are identical, no action

            else:
                res['flat']['Q'] = zQ  # do not exist, so fill in (even if None)
                
            #if VB>=0: print 'DEBUG Q found for block %-10s : %s' %(block, zQ)
            #print 'flat: ', res['flat']
        except:
            pass
        # ---
        
        if type(slha) is dict:  # dict

            if VB>3: print('slha is dict (block=%s)' %(block))
            
            for key in list(slha.keys()): 

                if key not in c2n:
                    res['warn'].append('Not in slhablockcode2name: block=%s   key=%s' %(block, key))
                    continue 


                # dict in dict (ex: mixing matrices)
                if type(slha[key]) is dict:
                    if VB>3: print('dict in dict: ', block, key)
                    for key2 in list(slha[key].keys()):

                        if VB>3: print('key2: ', block, key, key2, type(key2), c2n, c2n[key], type(c2n[key]))
                        if key2 not in c2n:
                            res['warn'].append('Not in slhablockcode2name: block=%s   key=%s  key2=%s' %(block, key, key2))
                            continue

                        AddToResDict(res=res, block=block, var=c2n[key][key2], val=slha[key][key2])
                        

                else:
                    if VB>3: print('just dict, key=%s' %(key))

                    AddToResDict(res=res, block=block, var=c2n[key], val=slha[key])


        else:
            if VB>3: print('slha is not dict (block=%s)' %(block))



    # 4) DERIVED VARIABLES #1
    if 'CALCver' in res['flat'] and type(res['flat']['CALCver']) is not str: res['flat']['CALCver'] = str(res['flat']['CALCver'])  # Hack, version number could be returned as float, want it as string
    
    if 'CALC' in res['flat']:

        # Combine CALC and CALCver
        if 'CALCver' in res['flat']: res['flat']['FULLCALC'] = res['flat']['CALC'] + '_' + res['flat']['CALCver']
        else:  res['flat']['FULLCALC'] = res['flat']['CALC']

        # Turn CALC into a ICALC
        if   'isasusy'  in res['flat']['CALC'].lower(): res['flat']['ICALC'] = 1
        elif 'suspect'  in res['flat']['CALC'].lower(): res['flat']['ICALC'] = 2
        elif 'softsusy' in res['flat']['CALC'].lower(): res['flat']['ICALC'] = 3
        elif 'spheno'   in res['flat']['CALC'].lower(): res['flat']['ICALC'] = 4
        else: res['flat']['ICALC'] = -1  # unknown



    # 4.5) REDUCED OBJECT RETURN: Relevant for slha files which only contain input (EXTPAR etc.) 
    if not 'MASS' in SLHA[0]:
        res['warn'].append('Returning reduced object: no MASS block in slha file.')
        if VB>0: print('WARNING::SlhaTools.slha2flatdict  ' + res['warn'][-1])
        return res


    # 5) DERIVED VARIABLES #2
    # neutralino composition, mixing angles
    F = res['flat']
    
    # a) neutralino composition: N1b,..,N4b, N1w, ... N1g,... N1h1, ... N1h2, ... N1h, ... N1g
    for i in ['1','2','3','4']: 
        F['N'+i+'b']  = F['N'+i+'1']**2
        F['N'+i+'w']  = F['N'+i+'2']**2
        F['N'+i+'g']  = norm2([F['N'+i+'1'],F['N'+i+'2']])
        F['N'+i+'h1'] = F['N'+i+'3']**2
        F['N'+i+'h2'] = F['N'+i+'4']**2
        F['N'+i+'h']  = norm2([F['N'+i+'3'],F['N'+i+'4']])

        for z in ['b','w','g','h1','h2','h']:
            res['var2block']['N'+i+z] = 'NMIX'
            res['vars']['NMIX'].append('N'+i+z)
        
    # b) L/R content of stop, sbottom, stau
    F['t1L'] = F['t11']**2
    F['t1R'] = F['t12']**2
    F['t2L'] = F['t21']**2
    F['t2R'] = F['t22']**2
    for z in ['t1L','t1R','t2L','t2R']:
        res['var2block'][z] = 'STOPMIX'
        res['vars']['STOPMIX'].append(z)

    F['b1L'] = F['b11']**2
    F['b1R'] = F['b12']**2
    F['b2L'] = F['b21']**2
    F['b2R'] = F['b22']**2    
    for z in ['b1L','b1R','b2L','b2R']:
        res['var2block'][z] = 'SBOTMIX'
        res['vars']['SBOTMIX'].append(z)

    F['T1L'] = F['T11']**2
    F['T1R'] = F['T12']**2
    F['T2L'] = F['T21']**2
    F['T2R'] = F['T22']**2
    for z in ['T1L','T1R','T2L','T2R']:
        res['var2block'][z] = 'STAUMIX'
        res['vars']['STAUMIX'].append(z)


    # implement mixing angles as well


    # c) wino/higgsino content of chargino ... 
    F['C1w_U'] = F['U11']**2
    F['C1h_U'] = F['U12']**2
    F['C2w_U'] = F['U21']**2
    F['C2h_U'] = F['U22']**2

    F['C1w_V'] = F['V11']**2
    F['C1h_V'] = F['V12']**2
    F['C2w_V'] = F['V21']**2
    F['C2h_V'] = F['V22']**2
    for i in ['1','2']:
        for wh in ['w','h']:
            for UV in ['U','V']: res['var2block']['C%s%s_%s' %(i,wh,UV)] = 'UMIX'  # and 'VMIX'
        

    # d) some additional variables
    F['GMEANstop'] = math.sqrt(F['~t1']*F['~t2'])
    F['GMEANsbot'] = math.sqrt(F['~b1']*F['~b2'])

    F['~q4Rmean'] = (F['~dR'] + F['~uR'] + F['~sR'] + F['~cR']) / 4.   # masses
    F['~q4Lmean'] = (F['~dL'] + F['~uL'] + F['~sL'] + F['~cL']) / 4.
    F['~q4mean'] = (F['~q4Rmean'] + F['~q4Lmean']) / 2. 

    F['q4Rmean'] = (F['dR'] + F['uR'] + F['sR'] + F['cR']) / 4.     # input mass parameters
    F['q4Lmean'] = (F['Q1'] + F.get('Q2',F['Q1'])) / 2.
    F['q4mean'] = (F['q4Rmean'] + F['q4Lmean']) / 2. 



    # 6 BRANCHING RATIOS
    decay = SlhaDecay(optD=opt)  #{'groups_use':opt['BR_groups_use']})
    decay.ReadFromSlha(slha1=SLHA[1])
    BRflatDict,BRflatOuts = decay.FlatDict()
    F.update(BRflatDict)


    # 7 FLAT OUTPUT (NON-STRUCTURED, ONLY ALPHABETICALLY)
    outs = []
    for var in sorted(F): 
        out = ' %-20s  %s' %(var, F[var])
        if type(F[var]) is float: out = ' %-20s  %15.9f' %(var, F[var])
        elif type(F[var]) is int: out = ' %-20s  %5i' %(var, F[var])
        else: out = ' %-20s  %s' %(var, F[var])
        outs.append(out)
    res['flat_outs'] = outs

    return res


# ##########
def norm(a):
    if not type(a) is list: a = [a]
    res = 0.
    for z in a: res += z*z
    return math.sqrt(res)

# ##########
def norm2(a):
    if not type(a) is list: a = [a]
    res = 0.
    for z in a: res += z*z
    return res


# ##########
import pickle
class SlhaDecay:
    '''
    NEXTs: 
    o headers... huh
    / move to separate script (& allow usage from prompt?)
    o Some negative BRs??
    o Validate sum
    o Add also possibility for e.g. N2_N1_l+q_* to slha2tlt
    '''

    
    # #####
    def __init__(s, fn='', optD={}, optL=[]):
        s.groups = {}
        s.groups['Q'] = ['d','u','s','c','b']
        s.groups['q'] = ['d','u','s','c']
        s.groups['q5'] = ['d','u','s','c','b']  # disfavoured?
        s.groups['q4'] = ['d','u','s','c']      # disfavoured?
        s.groups['l'] = ['e','m']
        s.groups['L'] = ['e','m','T']
        s.groups['v'] = ['ve','vm','vT']

        s.fn_slha = fn
        s.fn_base = s.fn_slha
        for ending in ['.slha','.lha','.SLHA','.LHA','.txt','.TXT']: 
            if s.fn_base.endswith(ending): s.fn_base = s.fn_base[:-len(ending)]

        s.optL = list(optL)
        s.optD = dict(optD)
        s.formatBR = optD.get('formatBR','%6.4f')
        s.minBR = optD.get('minBR',1e-3)

        s.BRkeystart = optD.get('BRkeystart','')
        if s.BRkeystart and not s.BRkeystart.endswith('_'): s.BRkeystart += '_'

        for group in optD.get('groups',[]):
            s.groups[group] = optD['groups'][group]

        s.groups_use = optD.get('groups_use',['q','v','l'])
        
        s.dec = {}

        if s.fn_slha: 
            s.ReadFromSlha(s.fn_slha)
        if 'flat' in optL: 
            s.FlatDict(fn_base=s.fn_base)
        if 'stdout' in optL:
            s.ShowDecay(fn_base=s.fn_base)



    # #####
    def StdKey(s,key): 
        # sort all, but put ~sparticles before SM (via replace-trick)
        keyL = []
        for zchild in key: keyL.append(zchild.replace('~','000~'))
        keyL.sort()
        for ichild,zchild in enumerate(keyL): keyL[ichild] = zchild.replace('000~','~')
        return tuple(keyL)


    # #####
    def UniqueKey(s, nam, key): 
        stdkey = s.StdKey(key)
        for k in list(s.dec[nam]['decay'].keys()):
            if sorted(stdkey) == sorted(k): return k,1
        return stdkey,0  # returns found-status (1 if found, 0 if not found)
       

    # #####
    def StarKey(s, nam, key): 
        match = []
        listkey = list(key)
        if '*' in listkey: listkey.remove('*')
        #print 'listkey: ', listkey
        for k in list(s.dec[nam]['decay'].keys()):
            listk = list(k)
            ismatch = 1
            #print 'listk: ',listk
            for element in listkey: 
                if element not in listk: ismatch = 0; break  #print 'break'; break
                else: listk.remove(element)
            if ismatch: match.append([nam, k])
        #print 'StarKey: match: ', match
        return match


    # #####
    def FlatDict(s, fn_base='', fn_pickle='', fn_outs=''): 
        
        flatdict = {}
        flatouts = []
        parent = sorted(s.dec, key=lambda name: abs(s.dec[name]['mass']), reverse=True)
        for ipar,zpar in enumerate(parent):
            decay = s.dec[zpar]['decay']
            for zchildren in sorted(decay, key=decay.__getitem__, reverse=True): 
                
                # Build full decay key
                BRkey = '%s%s' %(s.BRkeystart, zpar)
                for zchild in zchildren: BRkey += '_%s' %(zchild)
                if s.optD.get('droptilde',1): BRkey = BRkey.replace('~','')
                
                flatdict[BRkey] = decay[zchildren]
                flatouts.append(' %-20s %8.5f' %(BRkey, decay[zchildren]))

        if fn_base: fn_pickle = fn_base + '_BRflat.pickle' ; fn_outs = fn_base + '_BRflat.txt'
        if fn_pickle: 
            f = open(fn_pickle,'w')
            pickle.dump(flatdict,f)
            f.close()
        if fn_outs: 
            f = open(fn_outs,'w')
            for out in flatouts: f.write('%s\n' %(out))
            f.close()

        return flatdict, flatouts


    # #####
    def string2dict(s, decays): 
        res = {}
        keys = []
        if type(decays) is str: decays = decays.split(',,')
        for idec,zdec in enumerate(decays): 
            zparts = zdec.split('_')
            zparent = zparts.pop(0).split(',')
            zchildren = []
            for zpart in zparts: 
                zchildren.append(zpart.split(','))
            zchildren = tuple(zchildren)

            decres = s.FindDecays(zparent, zchildren)
            
            key = '%s' %(str(zparent).replace('[','').replace(']','').replace(' ','').replace("'",''))
            for zchild in zchildren: key += '_%s' %(str(zchild).replace('[','').replace(']','').replace(' ','').replace("'",''))
            if s.optD.get('tilde',0): key = key.replace('~','')

            res[key] = decres['totbr']
            keys.append(key)

        return res,keys  # keep the original order in keys


    # #####
    def ShowDecay(s, parent=[], children=[], matches=[], fn_base='', fn='', ret=0): 

        if fn_base: fn = '%s_BRstdout.txt' %(fn_base)

        if parent == []: 
            parent = sorted(s.dec, key=lambda name: abs(s.dec[name]['mass']), reverse=True)  # sorts according to mass

        if children: matches = s.FindDecays(parent, children)['match']

        zbrsum = 0.
        import operator
        outs = []
        if matches or children: # input, or by specifying children
            for match in matches:
                #print 'match: ', match
                zpar,zchildren,zbr = match
                zbrsum += zbr
                if zbr < s.minBR: continue
                out = '  ' + s.formatBR %(zbr) 
                txt = '%-3s -> ' %(zpar)
                for zchild in zchildren: txt += '%-3s ' %(zchild)
                if s.optD.get('tltbrname',0): txt = txt.strip().replace(' ->','').replace('   ',' ').replace('  ',' ').replace(' ','_')
                if s.optD.get('tltbrname',0) and s.optD.get('tilde',1): txt = txt.replace('~','')
                #for i in range(3-len(zchildren)): out += '%-3s ' %('-')    # optional
                outs.append('%s   %s' %(out,txt))

            # Hack to add summary line for complex requests
            if len(outs) > 1: 
                txt = 'SUM'
                if children: 
                    txt = '%-3s -> ' %(str(parent).replace('[','').replace(']','').replace(' ','').replace("'",''))
                    for zchild in children: txt += '%-3s ' %(str(zchild).replace('[','').replace(']','').replace(' ','').replace("'",''))
                if s.optD.get('tltbrname',0): txt = txt.strip().replace(' ->','').replace('   ',' ').replace('  ',' ').replace(' ','_')
                if s.optD.get('tltbrname',0) and s.optD.get('tilde',1): txt = txt.replace('~','')
                outs.append('  ' + s.formatBR %(zbrsum) + '   %s  [SUM]' %(txt))
            # --- 

        else: 
            # Show all decays for this given parents (or all)
            #sortedparents = sorted(s.dec, key=lambda name: s.dec[name]['mass'], reverse=True)
            for ipar,zpar in enumerate(parent):
                if zpar in s.optD.get('not2show',['t','W','Z','T','b']): continue 
                if not s.optD.get('minmass',0) <= abs(s.dec[zpar]['mass']) <= s.optD.get('maxmass',99999): continue
                if ipar > 0: outs.append('')
                if s.optD.get('showparthead',1): 
                    outs.append('PARENT: %-3s    MASS: %3.1f    WIDTH: ' %(zpar, s.dec.get(zpar,{}).get('mass',0)))
                    zwidth = s.dec.get(zpar,{}).get('width',0)      # simplistic hack to present width meaningfully
                    if zwidth > 0.2: outs[-1] += '%3.1f' %(zwidth)
                    else: outs[-1] += '%.1e' %(zwidth)
                sumbr = 0.
                decay = s.dec[zpar]['decay']
                for zchildren in sorted(decay, key=decay.__getitem__, reverse=True): 
                    dm = abs(s.dec.get(zpar,{}).get('mass',0))
                    if decay[zchildren] < s.minBR: continue
                    sumbr += decay[zchildren]
                    out = ''
                    if s.optD.get('showsum',1): out += '  ' + s.formatBR %(sumbr)
                    out += '  ' + s.formatBR %(decay[zchildren]) + '   %-3s -> ' %(zpar)
                    #print 'DEBUG:ShowDecay ', zpar, zchildren
                    mout = '%4.0f '%(abs(s.dec.get(zpar,{}).get('mass',0)))
                    for zchild in zchildren: 
                        out += '%-3s ' %(zchild)
                        #print zchild, s.dec['~N2']['mass']
                        zm = s.dec.get(zchild,{}).get('mass',0)
                        dm -= abs(zm)
                        mout += '%4.0f '%(abs(zm))
                    for i in range(3-len(zchildren)): out += '%-3s ' %('')    # optional, no, not anymore
                    if s.optD.get('showmassdiff',1): out += ' %6s' %('D=%.0f' %(dm))     # optional
                    if s.optD.get('showmass',1): out += '   %s' %(mout)       # optional
                    outs.append(out)
                    
        if fn: 
            f = open(fn,'w')
            for out in outs: f.write('%s\n' %(out))
            f.close()
        elif ret: 
            return outs
        else: 
            for out in outs: print(out)


    # #####
    def FindDecays(s, parent, children): 
        # this procedure is terribly complex ... and probably slow ; should also have a simplistic-fast one? 
        if type(parent) is list: parents = list(parent)
        else: parents = [parent]
        if type(children) is list: childrenlist = list(children)
        else: childrenlist = [children]
        # Ex: parents = ['~N2','~N3']  # usually has little meaning to use list of parents
        # Ex: childrenlist = [ (['~N1','~N2'],'q'), ('~C1','L','*') ]

        res = {'totbr':0, 'match':[]}

        for zparent in parents: 
            for zchildren in childrenlist: 
                # zchildren might also expand into many decays... 
                # here need to multiply out all combinations
                # could also implement simplifications, ex: '~NX' -> ['~N1','~N2','~N3','~N4'] ... 
                childD = {}
                for ichild,zchild in enumerate(zchildren): 
                    if not type(zchild) is list: zchild = [zchild]
                    childD[ichild] = zchild
                #print 'AAA: ', childD  
                zres = IterateNdim(varDict=childD)

                #print 'combinations: ', zres['combinations']

                # Loop over combinations
                for comb in zres['combinations']: 
                    zchnams = list(comb.values())
                    tchnams = tuple(zchnams)
                    #print tchnams

                    match = []
                    if '*' not in tchnams: 
                        key,keyexists = s.UniqueKey(nam=zparent, key=tchnams)
                        if keyexists: match = [[zparent,key]]
                    else: 
                        match = s.StarKey(zparent, tchnams)  # is a list (can be several keys)
                        
                    res['match'] += match
                    #zbr = s.dec[zparent]['decay'][key]  #[tchnams]
                    #res['totbr'] += zbr
                    #res['match'].append([zbr, zparent, tchnams])
                        

        # Now need to loop over matches and find 'totbr' (or rather)
        for i in range(len(res['match'])-1,-1,-1):
            if res['match'][i] in res['match'][0:i]: 
                #print 'popping ', i, res['match'][i]
                res['match'].pop(i)

        # Finally sum the results
        for imatch in range(len(res['match'])):
            zparent,key = res['match'][imatch]
            zbr = s.dec[zparent]['decay'][key]
            res['match'][imatch].append(zbr)  
            res['totbr'] += zbr
            #print imatch, key
            
        return res
                                          

    # #####
    def GetBR(s, parent, children): 
        res = FindDecays(parent, children)
        return res['totbr']


    # #####
    def ReadFromSlha(s, fn='', slha1=''):
        # Get slha1 object
        if fn:
            slha = pyslha.readSLHAFile(fn)
            slha1 = slha[1]
            s.fn_slha = fn
        if not slha1: print('No slha given') ; return {}

        # Traverse slha1 object
        #   for each mother, fill the br array
        pdgIDs = sorted(slha1)
        for pdgID in pdgIDs:

            nam = pdg2nam.get(pdgID,'')
            if nam == '':
                print('Skipping unknown pdgID: %i' %(pdgID))
                continue

            if nam in s.dec:
                print('Skipping %s (%i) since already exists (antiparticle?)' %(nam, pdgID))
                continue

            #if nam != 'A': continue  # debug

            p = s.dec[nam] = {}

            # inelegant hacks .. 
            try: p['mass'] = slha1[pdgID].mass
            except: p['mass'] = -999
            try: p['width'] = slha1[pdgID].totalwidth
            except: p['width'] = -9
            if p['mass'] is None: p['mass'] = -999   # top
            if p['width'] is None: p['width'] = -9   # top
            if p['mass'] < 0:  # hack
                if nam == 't': p['mass'] = 172.5   # APPROXIMATE. NOTE THIS IS PROBABLY IN SLHA AS WELL
                if nam == 'W': p['mass'] =  80.382 ; p['width'] = 2.085   # these are not included in the decay array, so will not happen
                if nam == 'Z': p['mass'] =  91.187 ; p['width'] = 2.4952  # better do something like in slha2finalstates
                if nam == 'T': p['mass'] =   1.776 ; p['width'] = 0.

            p['pdgID'] = pdgID


            # Then read in the decays (and join e.g. u,d,s,c,b into q if specified in s.groups_use)
            p['decay'] = {}
            for zdecay in slha1[pdgID].decays:
                zbr = zdecay.br
                zchids = zdecay.ids
                zchnams = []
                
                for zchid in zchids: 
                    zchnam = pdg2nam.get(abs(zchid),str(abs(zchid))) 

                    # Now gather (sign is already gone)
                    for group in s.groups_use:
                        if zchnam in s.groups[group]: zchnam = group
                    
                    zchnams.append(zchnam)
                    
                tchnams = tuple(zchnams)

                key,status = s.UniqueKey(nam=nam, key=tchnams)
                if status == 1: 
                    p['decay'][key] += zbr
                else:
                    p['decay'][key]  = zbr

                #print nam, key, p['decay'][key]
                                           
            # Then 

        # Then add SMs if not there: (not decays, only 
        if 't' in s.dec and s.dec['t']['mass'] < 0: 
            s.dec['t']['mass'] = 172.5   # Not necessarily correct
        if 'W' not in p: 
            s.dec['W'] = {'pdgID':24, 'mass':80.483, 'width':2.085, 'decay':[]}
        if 'Z' not in p: 
            s.dec['Z'] = {'pdgID':23, 'mass':91.187, 'width':2.4952, 'decay':[]}
        if 'T' not in p: 
            s.dec['T'] = {'pdgID':15, 'mass': 1.776, 'width':0., 'decay':[]}


    # #####

# ##########
