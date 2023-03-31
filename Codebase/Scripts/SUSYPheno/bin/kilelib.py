###############################
# library of generic functions
# b.k.gjelsten@fys.uio.no
###############################

import os, sys, time, stat, string, pickle, datetime
import resource
import inspect
import array


# ###################################################################
veryglobalDict = {}
def ReadDict(fn, name): # is there a sligthly more elegant way? 
    exec(compile(open(fn, "rb").read(), fn, 'exec'), veryglobalDict)
    return veryglobalDict[dictname]
# ###################################################################
class GeneaLog:
    # todo: timemode
    # todo: memmode
    
    def __init__(s, optD={}):
        s.strengthVB = optD.get('strengthVB',2)

        # Genealogy
        s.maxwidth1 = optD.get('maxwidth1', 10)  # this can usefully be set in the init stage 
        s.maxdepth = optD.get('maxdepth',5)   # This can usefully be set in the init stage. Depends on code. 
        s.overflowtxt = optD.get('overflowtxt','..')
        s.genea_delim = optD.get('genea_delim',':')
        s.width_strength = optD.get('width_strength',10)

        s.line_globalleft = optD.get('line_globalleft','')
        s.line_delim = optD.get('line_delim',' : ')
        s.line_beforetxt = optD.get('line_beforetxt','   ')

        s.dir = optD.get('dir','.')
        s.fn_base = optD.get('fn_base','geneaLog')

        s.strengthTxt = {0:'dump', 1:'debug', 2:'info', 3:'warning', 4:'error', 5:'fatal'}
        s.tosave = optD.get('tosave',['info','debug','dump'])
        if 'all' in s.tosave: s.tosave = list(s.strengthTxt.values())

        s.upperlower = optD.get('upperlower','upper')  # default is capitalise all 
        
        #s.realmaxwidth1 = s.maxwidth1 - len(s.overflowtxt)

        s.timemode = optD.get('timemode',0)
        s.memmode = optD.get('memmode',0)
        s.format_t = optD.get('format_t','%6.1f')
        

        s.mode = optD.get('mode',1)

        s.flush = optD.get('flush',0)  # if 1, will print when genea is changed in new algo
        s.indent = optD.get('indent','  ')

        s.t0_s = time.time()
        s.t0_dt = datetime.datetime.now()

        
        s.log = {}
        s.time = {}
        for strengthI in range(6): 
            s.log[strengthI] = []
            s.log[s.strengthTxt[strengthI]] = s.log[strengthI]  # link to the same
            #print strengthI, s.strengthTxt[strengthI]
            s.time[strengthI] = []
            s.time[s.strengthTxt[strengthI]] = s.time[strengthI]  # link to the same
            


    # #####
    def StrengthI(s, S='info', VB=0):
        strength = str(S).lower()
        res = 2
        if   strength in ['0','dump']: res = 0 #'dump'
        elif strength in ['1','d','debug']: res = 1 #'debug'
        elif strength in ['2','i','info']: res = 2 #'info'
        elif strength in ['3','w','warn','warning']: res = 3 #'warning'
        elif strength in ['4','e','err','error']: res = 4 #'error'
        elif strength in ['5','f','fatal']: res = 5 #'fatal'
        else:
            res = 2 #'info'
            print("GeneaLog: non-recognised strength level: %s.  Allowed ones: 'dump'=0, 'debug'=1, 'info'=2, 'warn'=3, 'error'=4, 'fatal'=5.  Setting to 'info'=2. " %(S))
        return res




    # #####
    def Add(s, txt, S='info', genea='', G='', store=1, ret=0, space=[0,0]):
        # store can be set to 0 (then not stored in log array)
        # ret can be set to 1 (then returned, not printed)
        # G and genea is the same (allow G for compact code)

        strengthI = s.StrengthI(S=S)
        if space[0] and not ret and strengthI >= s.strengthVB: print((space[0]-1)*' \n')

        tnow_s = time.time() - s.t0_s  # seconds since start
            
        # Allow sending in a list of texts
        if type(txt) is list:
            for ztxt in txt: s.Add(txt=ztxt, S=S, genea=genea, G=G)
            if space[1] and not ret and strengthI >= s.strengthVB: print((space[1]-1)*' \n')
            return 

        if genea == '': genea = G  # allow use of G for compact
        
        strengthTxt = s.strengthTxt[strengthI]
        if s.upperlower in ['upper']: strengthTxt = strengthTxt.upper()
        if s.upperlower in ['lower']: strengthTxt = strengthTxt.lower()
        if s.upperlower in ['capitalisefirst']: strengthTxt = strengthTxt.capitalize()  # first letter only


        # Make the line
        wG = genea.split(s.genea_delim)
        out = s.line_globalleft
        
        if s.mode == 1: # no other modes yet (probably never)

            # Genealogy
            # print 'maxdepth: ', s.maxdepth
            # print 'len(wG):  ', len(wG)
            if len(wG) > s.maxdepth:
                s.Add(genea='', S='warn', txt='WARNING WARNING  GeneaLog s.maxdepth need be increased from %i to %i' %(s.maxdepth, len(wG)))

            for iD in range(max(s.maxdepth, len(wG))):
                if iD < len(wG):
                    if iD > 0: out += s.line_delim
                    out += '%%-%is' %(s.maxwidth1) %(wG[iD])
                else:
                    out += (len(s.line_delim) + s.maxwidth1)*' '

            # Can add time & mem info here.. (to be implemented)
            if s.timemode == 1:
                tTxt = '  ' + s.format_t %(tnow_s) + 's'
                out += tTxt
                # can change s.format_t dynamically
                # can also implement more advanced time-prints

            if s.memmode > 0:
                pass

            # Strength
            
            out += '   %%-%is' %(s.width_strength) %(strengthTxt)

            # The message
            out += s.line_beforetxt
            out += len(wG) * s.indent
            out += txt
            



        # Store in arrays
        if store:  # D=1
            for istr in range(6): 
                if strengthI >= istr:                    
                    s.log[istr].append(out)
                    s.time[istr].append(tnow_s)

        # print? (or return)
        if ret: # D=0
            return out  # returns regardless of strength
        else:
            if strengthI >= s.strengthVB: print(out)   # 2014-08-28 (added test)
            if space[1] and strengthI >= s.strengthVB: print((space[1]-1)*' \n')


    # #####
    def Genea(s, genea=''):
        # useful to have inside class because use the same delim and maxwidth (Note: wouldn't need same maxwidth..)
        z = inspect.stack()[1][3]
        if len(z) > s.maxwidth1: z = z[0:s.maxwidth1-len(s.overflowtxt)] + s.overflowtxt
        genea2 = genea + s.genea_delim + z
        if genea2.startswith(s.genea_delim): genea2 = genea2[1:]

        if s.flush or 0: print(genea2)
        
        return genea2
        

    # #####
    def SaveToFile(s):
        for tosave in s.tosave: 
            fn = '%s/%s_%s.txt' %(s.dir, s.fn_base, tosave)
            tosaveI = s.StrengthI(tosave)
            WriteToFile(fn=fn, outs=s.log[tosaveI], VB=0)
        
    

# ###################################################################
def Genealogy(genea='', delim=':', maxwidth1=10, overflowtxt='..'):  # this is a standalone version, there is also one inside class kilelib::GeneaLog
    # useless because command needs to be issued where it is, cannot go into this subdir (hm, could have used [1][3]
    z = inspect.stack()[1][3]
    if len(z) > maxwidth1: z = z[0:maxwidth1-len(overflowtxt)] + overflowtxt
    genea2 = genea + delim + z
    if genea2.startswith(delim): genea2 = genea2[1:]
    return genea2


# ###################################################################
class Log:
    #import MemUsage
    # 2012-11-28
    def __init__(self, fn='', mode=['stdout'], timed=0, optD={}): 
        self.ctime0 = time.ctime()
        self.timed = optD.get('timed',timed)  # eventually may want to remove variable timed
        self.mem = optD.get('mem',0)
        self.memformat = optD.get('memformat','Mem: %4i MB')
        
        self.fn = fn
        self.mode = mode  # allowed entries: 'stdout','file','silent'
        self.rstrip = 1
        if self.fn: self.mode.append('file')


    def nostdout(self):
        while 'stdout' in self.mode: self.mode.remove('stdout')
        

    def setfile(self,fn):
        self.fn = fn
        if 'file' not in self.mode: self.mode.append('file')

    def log(self,txt):
        if type(txt) is str: txt = [txt]

        TIME=''
        if self.timed == 1: TIME = '%s | ' %(time.ctime())  # very simple

        MEM=''
        if self.mem == 1: MEM = self.memformat %(MemUsage('MB')) + ' | ' 

        if 'file' in self.mode: f = open(self.fn, 'a')

        for iT in range(len(txt)):
            T = txt[iT]
            if self.rstrip: T = T.rstrip()
            T = TIME + MEM + T
        
            if 'file' in self.mode: f.write('%s\n' %(T))
            if 'stdout' in self.mode: print(T)

        if 'file' in self.mode: f.close()

# ###################################################################
def CleanList(List, strip='rstrip'):
    inList = list(List)
    outList = []
    for L in inList:
        if strip == 'rstrip': outList.append(L.rstrip())  # default
        if strip == 'strip': outList.append(L.strip())
        if strip == 'lstrip': outList.append(L.lstrip())
    return outList
    
# ###################################################################
def PutOnOneLine(lines, delim='  |  '):
    oneline = ''
    for iL in range(len(lines)): 
        if iL > 0: oneline += '  |  '
        oneline += '%s' %(lines[iL])
    return oneline
     
# ###################################################################
def AlterPath(path, rel='../'):
    if not rel.endswith('/'): rel += '/'  # always? 
    if path.startswith('/'): path = path
    else: path = '../' + path
    return path
# ###################################################################
def FilenameReplace(fn, repl=[], safeapp='_safeapp', safepre='', app='', pre=''):
    fnout = fn
    nr = len(repl)
    if 2*(nr/2) != nr: return fn + '_incorrectReplaceString'

    for ir in range(nr/2):
        replthis = repl[2*ir]
        withthis = repl[2*ir+1]
        fnout = fnout.replace(replthis, withthis)

    fnout = pre + fnout + app
    
    if fnout == fn: fnout = safepre + fnout + safeapp
    if fnout == '': fnout = safepre + fnout + safeapp
    if fnout == '':
        print('Error::FilenameReplace  resulting filename is empty ... fn: %s  , repl: %s' %(fn, repl))
        fnout = 'fnwasempty'
    
 
    return fnout
# ###################################################################
def ReadPipeOrFile(pipe=0, f='', VB=0, excludestart=[], excludeempty=1, lstrip=0):
    out = []
    if pipe:
        if VB: print('INFO::ReadPipeOrFile:: Now reading from stdin')
        z = sys.stdin.readlines()
        #print z
    else:
        if os.path.exists(f):
            ff = open(f); z = ff.readlines(); ff.close()
            #for zz in z: out.append(zz.rstrip())
        else:
            print('ERROR::ReadPipeOrFile  Returning empty list; No such file %s' %(f))
            return []

    
    for zz in z:
        cont = 0
        zz = zz.rstrip()
        if lstrip: zz = zz.lstrip()
        if excludeempty and zz == '': cont = 1
        for beg in excludestart:
            if zz.startswith(beg): cont = 1
        if cont: continue
        out.append(zz)
        
    return out
    
# ###################################################################
def readfile_strip(fn, VB=0, txt='', opt=['strip']):
    if txt: print(txt)
    elif VB: print('readfile_strip:  %s' %(fn))
    f = open(fn,"r"); lines = f.readlines(); f.close()
    outs = []
    for l in lines:
        if 'strip' in opt: l = l.strip()
        if 'rstrip' in opt: l = l.rstrip()
        if 'lstrip' in opt: l = l.lstrip()
        outs.append(l)
    del f, lines
    return outs


# ###################################################################
def readlines_popen(cmd):
    a = os.popen(cmd).readlines()
    out = []
    for L in a: out.append(L.strip())
    del a
    return out


# ###################################################################
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


# ###################################################################
def numberofdecimals(s,VB=1):
    w = s.split(',')
    if len(w) == 1: return 0
    if len(w) == 2: return len(w[1])
    if VB>0: 'print not a number: %s' %(s)
    return -1


# ###################################################################
def getninitialspaces(a, VB=1):
    if not type(a) is str:
        if VB: print('getninitialspace: not a string: %s' %(a))
        return -1
    
    nspace = 0
    for i in range(len(a)):
        if a[i:i+1] != ' ': return nspace
        nspace += 1

    # can arrive here if string is only blanks
    return nspace


# ###################################################################
def getsplitspaces(L,VB=0):
    '''
    So for a string L you can regain the string by
    Lw = L.split()
    Ls = getsplitspaces(L)
    t = ''
    for i in range(len(Lw)):
      t += "%s%s" %(Ls[i],Lw[i])
    t += Ls[-1]
    '''
    
    w = L.split()
    L2 = L
    s = []

    for i in range(len(w)):
        s.append( getninitialspaces(L2) * ' ' )
        #s.append( getninitialspaces(L2) )
        L2 = L2.lstrip()
        L2 = L2[len(w[i]):]

        if VB>1:
            print(s)

    if VB>1: print('L2 is now: |%s|' %(L2))
    s.append( getninitialspaces(L2) * ' ' ) # include the last
    
    return s


# ###################################################################
def safeGet(thedict, thekey, ifnotindict=-1): 
    if thekey in thedict: return thedict[thekey]
    else: return ifnotindict



# ###################################################################
def IsInRangesT(RangesT, I, Imin0=-999999, Imax0=999999):
    # ex: RangesT = '101-105,121,201-205'
    if RangesT == '': return 0
    for Item in RangesT.split(','):
        
        # A) Check range
        if '-' in Item:
            inrange = 1
            Range = Item.strip().split('-')
            if Range[0] and I < float(Range[0]): inrange = 0
            if Range[1] and I > float(Range[1]): inrange = 0
            if inrange: return 1  # success

        # B) Check specific number
        else:
            if I == float(Item): return 1  # success

    # if arrive here, I is not in RangesT
    return 0


# ###################################################################
def ReadFileOrStdin(fn, vb=1, dorstrip=1):
    lines = []
    if fn:
        if os.path.exists(fn):
            f = open(fn)
            lines = f.readlines()
            f.close()
        else:
            print('FATAL::ReadFileOrStdin   non-existent file %s' %(fn))

    else: 
        if vb: print('INFO::ReadFileOrStdin  Getting input from stdin ...')
        lines = sys.stdin.readlines()

    linesout = lines
    if dorstrip:
        linesout = []
        for L in lines: linesout.append(L.rstrip())

    return linesout


# ###################################################################
def fileage(pathname,unit='s'):
    if not os.path.exists(pathname): return -1  # if no file, return negative
    res = time.time() - os.stat(pathname)[stat.ST_MTIME]
    if unit == 'm': res /= 60.
    if unit == 'h': res /= 3600.
    if unit == 'd': res /= 86400.
    return res  # default is unit='s' (seconds)


# ###################################################################
def LoadPickle(fn, VB=0, txt=''):
    if txt: print(txt)
    elif VB: print('LoadPickle:  %s' %(fn))
    f = open(fn)
    res = pickle.load(f)
    f.close()
    return res
# ###################################################################
def WritePickle(fn, thepickle, VB=1):
    PickleToFile(fn, thepickle, VB=VB)
def PickleToFile(fn, thepickle, VB=1):
    f = open(fn,'wb')
    pickle.dump(thepickle, f)
    f.close()
    if VB>0: print('INFO::PickleToFile  created file %s with %i entries' %(fn, len(thepickle)))
# ###################################################################
def WriteToFile(fn, outs, wORa='w', opt=[], VB=1):
    if not type(outs) is list: outs2 = [outs]
    else: outs2 = list(outs)
    if wORa not in ['w','a']:
        print("Warning  WriteToFile  Illegal wORa '%s', allowed ones are 'w','a'. Using 'w'." %(wORa))
        wORa = 'w'
    f = open(fn,wORa)
    for out in outs2:
        if 'rstrip' in opt: out = out.rstrip()
        f.write('%s\n' %(out))
    f.close()
    #if VB>0 and len(outs2) == 0: print 'INFO::WriteToFile  intended to create/append to file %s, but ' %(fn)
    if VB>0 and wORa == 'w': print('INFO::WriteToFile  created file %s with %i lines' %(fn, len(outs2)))
    if VB>0 and wORa == 'a': print('INFO::WriteToFile  appended %i lines to file %s' %(len(outs2), fn))


# ###################################################################
def AlignColumns(lines, align='r', n=2, maxcol=999, rstrip=0, lspace=0, insertL=[], insertTopical={}):
    # should allow to also shuffle columns
    # insertTopical={'keycol':1}  # minimum
    # insertTopcial={'keycol':0, 'keysubstr':[0,3], 'insert':['apace','head','hline','hline'], 'hlinechar':'-' }  # MAXIMUM

    # insertL=[['-',0,2,999]]  # puts line above and below title as well as at the end
    #   Note that the numbers given are the rows of the lines in the final table, i.e. including the lines themselves
    # lines = lines[0:4]
    outs = []

    myspace = n
    align = align.replace('  ','').replace(' ','')  # remove any white space inserted for overview
    just = align.split(',')
    #if just == []: just = ['r']
    njust = len(just)

    Nlines = len(lines)
    
    maxwidth = []

    # print 'DEBUG njust: ', njust

    # find maxwidth
    for iL in range(Nlines):
        line = lines[iL].strip()
        word = string.split(line)
        for iw in range(len(word)):
            # fills out maxwidth upon need
            if(iw==len(maxwidth)):
                maxwidth.append(0)
            if(len(word[iw])>maxwidth[iw]):
                maxwidth[iw]=len(word[iw])

    # print 'DEBUG maxwidth: ', maxwidth

    for iL in range(Nlines):
        line = lines[iL].strip()
        word = string.split(line)
        outline = ""
        for iw in range(len(word)):

            if iw > maxcol-1:
                outline+=" "+word[iw]
                continue
        
            if iw > 0: outline += myspace*" "
        
            if iw > njust-1: 
                thisjust = just[njust-1]
            else:
                thisjust = just[iw]

            # print 'DEBUG thisjust: ', thisjust
                
            if thisjust in ('l','L'):
                outline+=word[iw].ljust(maxwidth[iw])
            elif thisjust in ('r','R'):
                outline+=word[iw].rjust(maxwidth[iw])
            elif thisjust in ('c','C'):
                outline+=word[iw].center(maxwidth[iw])
            else:
                print("WARNING::AlignColumns: line/col %i/%i : non-allowed code: %s  [line=%s]" %(iL, iw, thisjust, line))
                # ... could make this more flexible..?
        
        if rstrip: outline = outline.rstrip()

        outs.append(lspace*' ' + outline)


    # ---- Insert hlines advanced
    if insertTopical and len(outs) > 0:
        Top = insertTopical
        oldkey = ''
        wid = len(outs[0])
        insertBase = {}
        insertBase['hline'] = lspace*' ' + (wid-lspace)*Top.get('hlinechar','-')
        insertBase['head'] = outs[Top.get('headcol',0)]
        insertBase['space'] = ''
        insertCode = Top.get('insert',['hline'])
        #for z in Top.get('insert',['hline']):
        #    toinsert.append(insertBase[z])
            
        for iL in range(len(outs)-1, 0-1, -1):
            L = outs[iL]
            w = L.split()
            if w[Top['keycol']] < len(w)-1: continue
            if iL < len(outs)-1: oldkey = key
            key = w[Top['keycol']]
            if 'keysubstr' in Top: # allows to choose substring
                key = key[ Top['keysubstr'][0] : min(len(key),Top['keysubstr'][1]) ]

            if iL == len(outs)-1: continue  # do nothing on last line
            if iL == 0 and not 'alsofirstline' in Top.get('opt',[]): continue # do not insert between line 0 and 1 unless 'alsofirstline' is in Top['opt'] - this protects from duplication of the header

            if key != oldkey:
                # here do the insertion
                for iz in range(len(insertCode)):
                    z = insertBase[insertCode[iz]]
                    outs.insert(iL+iz+1, z)



    # ---- Insert hlines (if any)
    if insertL and len(outs) > 0: 
        wid = len(outs[0])
        for zz in insertL:
            hline = lspace*' ' + (wid-lspace) * zz.pop(0)   # faulty? 
            for z in zz:
                outs.insert(int(z), hline)
    # ----

    return outs


# ###################################################################
def ReconstructCommandLine(sysargv, optL=['popcommand','stamp'], stampin=''):   # From prod_requester ... even more specific there
    if 'stampin' == '': 
        t = list(time.localtime())
        for iT in range(len(t)): t[iT] = str(t[iT])
        stamp = t[0].zfill(4)+'-'+t[1].zfill(2)+'-'+t[2].zfill(2)+'_'+t[3].zfill(2)+'-'+t[4].zfill(2)+'-'+t[5].zfill(2)
    else:
        stamp = stampin

    if 'stamp' in optL: 
        commandline = 'stamp=%s ; ' %(stamp)
    else:
        commandline = ''
        
    for iz in range(len(sysargv)):
        z = sysargv[iz]
        if 'popcommand' in optL and iz == 0: z = z.split('/').pop()
        if z.startswith('-'): commandline += ' '
        commandline += ' %s' %(z)
        
    #if 'argv' in optD: commandline += '    # as object'  # a bit too specific, maybe
    return commandline


# ###################################################################
def SmartInterval(i1,i2):
    # Will turn e.g. (172004,16) to (172004,172016)
    ti1 = str(i1)
    ti2 = str(i2)
    zdecdiff = len(ti1) - len(ti2)
    ndectot = len(ti1)
    if zdecdiff > 0: 
        z1 = ti1[:zdecdiff]
        i2add = int(z1) * 10**(ndectot-zdecdiff)
        i2new = i2add + i2

    i2 = i2new
    if i2 < i1:
        print('Warning::SmartInterval  i1 > i2 (%i>%i)  Probably wrong input' %(i1,i2))
        
    return [i1,i2]


# ###################################################################
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
# ###################################################################
def ReplaceText(txtin, repl):
    # txtin can be string or list of strings
    if type(txtin) is list: txtin2 = txtin
    else: txtin2 = [txtin]
    
    txtout = []
    for L in txtin2:
        for key in repl:
            if key in L: L = L.replace(key, repl[key])
            
        txtout.append(L)

    del txtin2
        
    if type(txtin) is list: return txtout
    else: return txtout[0]
        

# ###################################################################
def CommandExists(cmd, VB=0):
    # Note: doesn't find aliases (ex. os.popen('command -v ll') gives no output inside python)
    z = os.popen('command -v %s' %(cmd)).readlines()
    if VB>0: print(z)  # only for debugging
    if len(z) == 0: return False
    else: return True


# ###################################################################
def GetFileInfo(fn):
    import os
    if not os.path.exists(fn): return {}
    
    res = {}
    
    L = os.popen('ls -l %s' %(fn)).readlines().pop(0).strip()
    w = L.split()
    res['user'] = w[2]
    res['chmode'] = w[0]
    
    res['stat'] = os.stat(fn)
    res['st_uid'] = res['stat'].st_uid
    res['size'] = res['stat'].st_size
    try:
        import pwd
        pw = pwd.getpwuid(res['stat'].st_uid)
        res['user'] = pw.pw_name  # overwrites
        res['gecos'] = pw.pw_gecos
        res['name'] = res['gecos'].split(',')[0]
        res['name_nospace'] = res['name'].replace(' ','')
    except:
        pass
    
    return res


# ###################################################################
def HumanRound(x, delim='', factor=1., decimal=0, hack=0):
    format = "%%.%if" %(decimal)
    unit = ''
    if   x > factor * 1e24: res = format %(x/1e24) ; unit = 'Y'
    elif   x > factor * 1e21: res = format %(x/1e21) ; unit = 'Z'
    elif   x > factor * 1e18: res = format %(x/1e18) ; unit = 'E'
    elif   x > factor * 1e15: res = format %(x/1e15) ; unit = 'P'
    elif   x > factor * 1e12: res = format %(x/1e12) ; unit = 'T'
    elif   x > factor * 1e9 : res = format %(x/1e9 ) ; unit = 'G'
    elif   x > factor * 1e6 : res = format %(x/1e6 ) ; unit = 'M'
    elif   x > factor * 1e3 : res = format %(x/1e3 ) ; unit = 'k'
    else: res = format %(x) ; unit = ''

    if hack == 1:
        if decimal > 0 and float(res.split('.')[0]) >= 10: res = res.split('.')[0]
    if hack == 2: 
        if decimal > 0 and not res.startswith('0.'): res = res.split('.')[0]

    return res + delim + unit
    
# ###################################################################
def MemUsage(unit='kB'):
    res = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if unit == 'B': res *= 1024
    elif unit == 'kB': res *= 1.
    elif unit == 'MB': res /= 1024.
    elif unit == 'GB': res /= (1024.*1024)
    elif unit == 'TB': res /= (1024.*1024*1024)
    else:
        print('Warning: unsupported unit %s, will return in kB' %(unit))
    res = int(res)
    return res

# ###################################################################
def AreNumbers(z, delim=' ', nanisnumber=1):  # Note nan,NaN,.. is treated as a number
    if type(z) is list: w = z
    elif type(z) is str: w = z.strip().split(delim)
    else:
        print('Warning IsNumbers unknown input type: %s  %s' %(type(z), z))
        return False

    # testing
    for ww in w:
        try:
            float(ww)
            if not nanisnumber and ww.lower() in ['nan']: return False
        except: return False
    return True
    
# ###################################################################
def RecursiveTraverseAll(varDict, varKeys=[], indices=[], combinations=[]):
    '''

    '''
    # --- Init
    if varKeys == []:
        varKeys = list(varDict.keys())
        varKeys.sort()
    if indices == []:
        indices = len(varKeys)*[0]
    
    # --- 


# ###################################################################
#def IterateNext(index, indexMax):   # Not needed, done inside the Iterate
#    index[-1] += 1
#    for ipos in reversed(range(len(index))):
#        if index[ipos] == indexMax[ipos]:
#            index[ipos] = 0
#            index[ipos-1] += 1
#    return index
#    #if index == indexMaxMinus1:
#       

# ###################################################################
def IterateNdim(varDict, vars=[]):
    # Usage: res = IterateNdim(varDict) where e.g. varDict = {'M1':[100,200,300], 'M2':[110,310], 'M3':[120,320]}
    # Option: can insert list of vars, e.g. vars=['M1','M3','M2'] to get them iterated in the desired order
    # Output is a result dict: with all combinations, both in index form and dict form
    # Created to be used by signalgrid_loop.py [2013-09-16]
    # Should rewrite such that it also can take a plain ndim list as input 

    # --- Init, max the maxlist, e.g. [5,3,3,5,12] and the start index [0,0,0,0,0]
    if vars == []:
        vars = list(varDict.keys())
        vars.sort()

    indices = []
    index = len(vars)*[0]
    indexMax = []
    indexMaxMinus1 = []
    for ivar in range(len(vars)):
        indexMax.append(len(varDict[vars[ivar]]))
        indexMaxMinus1.append(len(varDict[vars[ivar]])-1)

    # --- Loop
    indices.append(list(index))
    while index != indexMaxMinus1:
        index[-1] += 1
        for ipos in reversed(list(range(len(index)))):
            if index[ipos] == indexMax[ipos]:
                index[ipos] = 0
                index[ipos-1] += 1
        indices.append(list(index))


    # --- Results
    combinations = []
    for index in indices:
        D = {}
        for ivar in range(len(index)):
            var = vars[ivar]
            ind = index[ivar]
            D[var] = varDict[var][ind]
        combinations.append(D)

    res = {'indices':indices, 'combinations':combinations, 'vars':vars}
    return res

# ###################################################################
def ReadFileOrStdin(fn='', stdin=0, opt=['rstrip'], skip=[], VB=0):
    # read files
    # if no files given or if stdin==1: read from stdin
    lines = []

    # Read file(s)
    for zfn in fn.split(','): 
        if zfn != '':
            if not os.path.exists(zfn):
                print("WARNING  ReadFileOrStdin  Non-existent file %s" (zfn))
            else: 
                f = open(zfn)
                lines += f.readlines()
                f.close()


    # Read stdin
    if fn == '' or stdin:
        if VB: print("INFO::ReadFileOrStdin  Reading from stdin ...")
        lines += sys.stdin.readlines()


    # Treat input
    Lines = []
    for L in lines:
        # Fix line
        if 'rstrip' in opt: L = L.rstrip()
        if 'lstrip' in opt: L = L.lstrip()
        if 'strip' in opt: L = L.strip()

        # Skip lines
        if 'empty' in skip and L.strip() == '' : continue
        if 'beginswith#' in skip and len(L)>0 and L[0] == '#' : continue

        Lines.append(L)

    return Lines

# ###################################################################
def GetSubsec(ndigit=6):
    T = time.time()
    ms = str(T-int(T)).replace('0.','')  # now we have a string
    if len(ms) >= ndigit: return ms[:ndigit]
    else: return ms + (ndigit-len(ms))*'0'
    
    
# ###################################################################
def TimeStamp(opt=['-'], drop=[], subsecdigits=6):
    # 2014-09-11:  MANY/MOST OF THESE CAN BE DONE WITH DATETIME, E.G. timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    # Ex: mystamp_s = timestamp()['sec']


    subsec = GetSubsec(ndigit=subsecdigits)  # note is string, but is ok here
    T0 = time.localtime()
    T = []
    for iT in range(6): T.append(T0[iT])
    T.append(subsec)
    
    res = {}

    # Init
    times = ['year','month','day','hour','min','sec','subsec']
    unit = ['-','-','_','h','m','s','']
    #Zfill = {'year':4, 'subsec':subsecdigits}
    
    # Fill res according to opt
    fullstamp = ''
    for iT in range(len(times)):
        t = times[iT]
        if t not in drop:
            #fullstamp += str(T[iT]).zfill(Zfill.get(t,2)) + unit[iT]
            fullstamp += str(T[iT]).zfill(2) + unit[iT]
            res[t] = fullstamp

    # Fix style according to opt
    for k in res:
        for zopt in opt: 
            if zopt in ['shortest','short','-',':']: res[k] = res[k].replace('h',':').replace('m',':').replace('s','_')
            if zopt in ['shortest','short','-']: res[k] = res[k].replace(':','-')
            if zopt in ['shortest','short']: res[k] = res[k].replace('-','')
            if zopt in ['shortest']: res[k] = res[k].replace('_','')

        if res[k][-1] in ['_','-',':']: res[k]=res[k][:-1]

    return res


# ###################################################################
def ExtractVarsViaFormat(txt, format):
    w = format.split('*')
    vals = []
    for iw in range(len(w)):
        #print 'before: txt: %s' %(txt)
        zi = txt.find(w[iw])
        #print 'txt: %s   zi: %i   vals: %s' %(txt, zi, vals)
        if zi>0: vals.append(txt[:zi])
        txt = txt[zi+len(w[iw]):]
        iw += 1

    if txt: vals.append(txt)

    return vals
# ###################################################################
def fromTimestamp(secs, now='', optD={}):
    return fromSecsSince1970(secs=secs, now=now, optD=optD)
# ###################################################################
def fromSecsSince1970(secs, now='', optD={}):
    import datetime
    if not now: now = datetime.datetime.now()
    res = {}
    res['now'] = now

    format_ndays = optD.get('format_ndays','%.1f')
    
    dt = datetime.datetime.utcfromtimestamp(int(secs))
    res['tdel'] = now - dt
    res['ndays'] = res['tdel'].days
    res['nhours'] = res['tdel'].seconds/3600.  # in addition to the days
    #res['NDAYS'] = res['ndays'] + res['nhours']/24.
    res['ndaysTxt'] = format_ndays %(res['ndays'] + res['nhours']/24.)  # can make 
    res['timestring_m'] = '%4s-%2s-%2s_%2s-%2s' %(str(dt.year).zfill(4), str(dt.month).zfill(2), str(dt.day).zfill(2), str(dt.hour).zfill(2), str(dt.minute).zfill(2))

    del dt

    if 'ret' in optD: 
        rets = ()
        for key in optD['ret']: rets += (res[key],)
        del res
        return rets

    return res

# ###################################################################
def UpdateIndexHtml(Dir='.', lsadd='', pre=[], post=[], fn_index='index.html', VB=1):
    outsD = {}
    zfns = os.popen("ls -l --time-style=long-iso  %s/ %s" %(Dir, lsadd)).readlines()
    for z in zfns:
        w = z.strip().split()
        if len(w) != 8:
            print('WARNING  UpdateIndexHtml  Remaking dir %s   Skipping irregular ls -l line: %s' %(Dir, z.strip()))
            continue
        zsize, zdate, ztime, zfn = w[4:8]
        outsD[zfn] = '%s  %s  %s  %s' %(zsize,zdate,ztime,zfn)
        outsD[zfn] = {'size':zsize, 'date':zdate, 'time':ztime, 'fn':zfn}
    outs = pre 

    zfns2 = list(outsD.keys())
    zfns2.sort()
    zfns2.reverse()
    for zfn in zfns2:
        #outs.append("<a href='%s'> %s </a> <br>" %(zfn, outsD[zfn]))
        D = outsD[zfn]
        outs.append("%s &nbsp; %s &nbsp; <a href='%s'> %s </a> &nbsp;&nbsp; (size=%s)<br>" %(D['date'], D['time'], D['fn'], D['fn'], D['size']))

    outs += post

    fn_index = '%s/%s' %(Dir, fn_index)
    WriteToFile(fn=fn_index, outs=outs, VB=VB)


# ###################################################################
def AgeOfLast(lscmd, unit):  # Fragile, ex: AgeOfLast("ls -tr dir_hist/* | grep pickle", unit="d")
    z = os.popen(lscmd).readlines()
    if z == []: return -2
    return fileage(pathname=z.pop().strip(), unit=unit)
    
# ###################################################################
def IndicesOfDuplicates(w):
    duplicate = []
    for iw in range(len(w)-1, -1, -1):
        if w[iw] in w[:iw]: duplicate.append(iw)
    return duplicate

# ###################################################################
def DictlistsFromListOfDicts(dictORfnpickleList, first=1, vars2pick=[], warn=[], consistencycheck=1):
    # vars2pick: if nonempty, this list specifies which variables to allow into the Dictlist

    VB=0 # debug
    dictlist = {}
    
    for idictORfnpickle in range(len(dictORfnpickleList)):
        dictORfnpickle = dictORfnpickleList[idictORfnpickle].strip()
        loaded = 0
        if type(dictORfnpickle) is str:
            #print dictORfnpickle
            d = LoadPickle(dictORfnpickle)
            loaded = 1
            if VB>1: print('DDD: loaded: d = ', d)
            
        # Then go through dict and add to Dictlists
        for key in list(d.keys()):
            if vars2pick and key not in vars2pick:
                if VB>1: print('DEBUG: since not in vars2pick, dropping key %s' %(key))
                continue   # do not keep this var
            
            if idictORfnpickle == 0 and first == 1:
                if VB>1: print('DEBUG: start the list for key %s' %(key))
                dictlist[key] = []   # start the list

            # check
            if key not in dictlist: 
                warn.append("Warning::GetDictlistFromListOfDicts  new variable found in (first=%i)  dict %2i :  %s   (SKIPPING)" %(first, idictORfnpickle, key))
                print(warn[-1])
                continue
                
            dictlist[key].append(d[key])

        #if loaded: del d
        

    if VB>1: print('DEBUG: dictlist: %s' %(dictlist))        


    # Consistency check at the end
    if consistencycheck:
        keys = sorted(dictlist.keys())
        nvals = {}
        for key in keys: nvals[key] = len(dictlist[key])
        if max(nvals.values()) != min(nvals.values()):
            warn.append('Error::GetDictlistFromListOfDicts: Inconsistent number per array min(len)=%i VS max(len)=%i' %(min(nvals.values()), max(nvals.values())))
            print(warn[-1])
            for key in keys: print('  %-10s  len: %4i' %(key, nvals[key]))
            print(warn[-1])  # print once more

            
    return dictlist


# ###################################################################
def DictFromTable(table, vartype='float', returnhead=0):
    # returns e.g. res = { 'M1':[100,100,100,100], 'M2':[200,300,200,300], 'mu':[200,200,300,300], 'N2':[220,270,273,330] }
    #         and optionally the header as given (to preserve the order) [usually not relevant]
    # the lists are ordered according to the table list; simple, straightforward
    # robust: if there should happen to be several columns with the same header, this is now taken care of
    # limitation: vartype is the same for all vars
    
    res = {}
    head = []
    for iL in range(len(table)): 
        L = table[iL].strip()
        w = L.split()
        if iL == 0:
            head = w
            #if len(set(head)) != len(head): print 'Warning::DictFromTable   Table has non-unique headers: %s' %(head)
            duplicates = IndicesOfDuplicates(head)
        for iw in range(len(w)):
            if iw in duplicates: continue   # takes care of tables where the same header category appears several times 
            if iL == 0: res[head[iw]] = []
            else:
                val = w[iw]
                if type in ['float','f','F']: val = float(val)
                if type in ['int','i','I']: val = int(val)
                res[head[iw]].append(val)

    if returnhead: return res, head
    return res
    
# ###################################################################
def GetPlainArray(table, var, arraytype='', ind=-1, VB=1, scale=1., protection=0):
    # The var is a complex variable
    vals = []
    tablekeys = list(table.keys())

    splits = {} 
    splits[0] = {'+':' + ', '-':' - ', '*':' * ', '/':' / ', '(':' ( ', ')':' ) '}  # default replacement
    splits[1] = {'[+]':' + ', '[-]':' - ', '[*]':' * ', '[/]':' / ', '[(]':' ( ', '[)]':' ) '}  # protected1, needed for XS:C1+N2
    splits[2] = {'[[+]]':' + ', '[[-]]':' - ', '[[*]]':' * ', '[[/]]':' / ', '[[(]]':' ( ', '[[)]]':' ) '}  # protected 2 (ever needed?
   

    if VB>1: print("GetPlainArray: var = %s" %(var))
    
    secure = {'->': '___>'}

    # Secure some standard combinations: part 1    
    for z in secure: var = var.replace(z,secure[z]) 

    # Insert whites (depending on protection mode)
    #for z in ['+','-','*','/','(',')']: var = var.replace(z,' '+z+' ')
    for z in splits[protection]: var = var.replace(z, splits[protection][z])
    if VB>1: print('var: %s' %(var))

    # Secure some standard combinations: part 2
    for z in secure: var = var.replace(secure[z],z) 


    # Find relevant dicts, find (&check) length of arrays
    wvar = var.split()
    if VB>1: print('var: %s    wvar: %s' %(var, wvar))
    tablelength = {}
    for z in wvar: 
        if z in tablekeys:
            #var2 += table[z][i]
            tablelength[z] = len(table[z])
    if len(tablelength) == 0:
        print('\nERROR::GetPlainArray : In table.keys(): %s \n   var %s is not found\n[WILL NOW CRASH]\n' %(tablekeys, var))
    thelength = min(tablelength.values())
    if min(tablelength.values()) != max(tablelength.values()):
        print("Warning::GetPlainArray  table arrays have different lengths: %s     (var=%s)" %(tablelength, var))

    # Loop to create array
    for ival in range(thelength):
        if ind > -1 and ind != ival: continue  # additional feature
        # Create expression that can be evaluated
        valT = ''
        for iz in range(len(wvar)):
            z = wvar[iz]
            if z in tablekeys:
                #valT += str(table[z][ival])
                valT += str(float(table[z][ival]))  # 2014-02-11: to allow e.g. '0158' as input (fails in eval)
                #print 'here: ', str(table[z][ival]), float(str(table[z][ival]))
            else: valT += ' ' + z + ' '

            #print ival, iz, z, valT
        val = eval(valT)

        vals.append(val * scale)  # scale allows to scale the value by some factor
        #print ival, var, valT, val

    # Finished
    if arraytype != '': vals = array.array(arraytype,vals)   # transform list->array if desired (default: no)
    
    if ind > -1 and len(vals) == 1: return vals[0]  # additional feature
    return vals

# ###################################################################
def SaveHistory(fn, argv, opt=['stripcmd'], VB=0):
    #opt=['dir','space','stripcmd','date','dirsameline']
    outs = []

    # dir
    if 'dir' in opt: 
        outs.append('cd  %s/' %(os.getcwd()))
        
    # command line [w/date]
    out = ''
    for iz in range(len(argv)):
        z = argv[iz]
        if iz == 0 and 'stripcmd' in opt: z = z.split('/').pop()
        if z.startswith('-'): out += "  %s" %(z)
        else: out += " %s" %(z)
    out = out.strip()
    # date
    if 'dirsameline' in opt: out += '        # cd  %s/' %(os.getcwd())
    if 'date' in opt: out += '       # %s' %(TimeStamp()['sec'])
    outs.append(out)

    # space
    if 'space' in opt:
        outs.append('')
    

    WriteToFile(fn=fn, outs=outs, wORa='a', VB=VB)


# ###################################################################
def MakeDirIfNeeded(DIRs, VB=1, pretext=''):
    if type(DIRs) is str: DIRs = [DIRs]

    for DIR in DIRs:
        if not os.path.exists(DIR):
            if VB>0: print('%sMakeDirIfNeeded  Creating dir %s' %(pretext, DIR))
            os.mkdir(DIR)
    

# ###################################################################
def tlt2dict(fn, valtype='f'):
    res = {}
    f = open(fn) ; lines = f.readlines(); f.close()
    if len(lines) != 2: print('Error::tlt2dict  tlt %s does not have 2 lines. Returning {}' %(fn)); return res
    wvar = lines[0].strip().split()
    wval = lines[1].strip().split()
    if len(wvar) != len(wval): print('Error::tlt2dict  tlt %s has unequal header and value lines. Returning {}' %(fn)); return res
    for i in range(len(wvar)):
        var = wvar[i]
        if valtype in ['f','F']: val = float(wval[i])
        elif valtype in ['i','I']: val = int(float(wval[i]))
        else: val = wval[i]  # string
        res[var] = val
    return res
        
# ###################################################################
def dict2tlt(tltdef, thedict, D={}, formats={}, var2block={}, varcollection={}, ret=1, fn_table='', prevar='', preval='', VB=1, delim_var=',', delim_rename='::', format0='%.3f'):
    from lineup import lineup
    # (copyinspired by slha2tlt)
    # this code is too specialised to be here (the direct use of var2block)
    warn = []
    # A) First construct the variable list
    # table can composite, can structure like this: table=name1,name2:fnadd
    # this will take vars from cat1 and cat2 and the file will get filename from fnadd
    # if fnadd not given, the first part (full part) will also give filename (a bit odd)

    vars = []
    cats = tltdef.split(delim_var)


    rename = {}
    for cat in cats:
        z = cat.split(delim_rename)
        if len(z) == 2: 
            cat = z[0]
            rename[cat] = z[1]
        if cat in varcollection: 
            vars += varcollection[cat]  # could ensure that there is no duplicates
        else: 
            vars += [cat]  # assume it is a free variable

    #if VB>1: print 'tltdef: %s   vars: %s' %(tltdef, vars)
                
    head = prevar 
    out  = preval 
    for var in vars:
        if var not in thedict:
            #s.warn.append('%s  Table %s: var %s not in flatdict. Value set to %s' %(fn_base, table, var, D['val_missing']))
            if D.get('tlt_varnotfound_include',1): 
                head += '  ' + D.get('tlt_varprepend','') + rename.get(var,var) + D.get('tlt_varappend','')
                out  += '  ' + D.get('tlt_varnotfound_val','-99997')
                warn.append("%s  Table %s: var %s not in flatdict. Setting value to %s.  (To skip, use ' -dict I,varnotfound_include,0') " %('fn_base', tltdef, var, D.get('tlt_varnotfound_val','-99997')))
            else:
                warn.append("%s  Table %s: var %s not in thedict. Skipping value (fragile)  (To keep, use ' -dict I,varnotfound_include,1')" %(fn_base, tltdef, var))
            continue
        #print var
        #print var2block.keys()
        if var in formats: 
            format = formats[var]  # var-specific (rare)
        elif var2block.get(var,'asdfasdf') in formats:
            #print 'here'
            format = formats[var2block[var]] 
        else: format = format0
        # Hmm ... the usage of var2block ruins the reusage of pickles ... [STILL THE CASE?]
        
        val = thedict[var]
        if val == None: val = D.get('tlt_valisnone_val',-99995)
        head += '  ' + D.get('tlt_varprepend','') + rename.get(var,var) + D.get('tlt_varappend','')
        out  += '  ' + format %(val)
        

    # Lineup and return/write to file 
    outs = lineup([head,out], n=D.get('tlt_nspaces',2))
    
    
    if ret: return outs
    if fn_table: WriteToFile(fn=fn_table, outs=outs, VB=VB-1)




# ###################################################################
def TableToHtml(outs, remove_cols=[], insert_cols=[], insert_lists=[]):
    # (this algo is taken from prod_requestMonitor and very slightly generalised)
    from lineup import lineup

    # Test:
    for List in insert_lists:
        if len(List) != len(outs): print('Warning::TableToHtml  len(outs) = %i  !=  %i = len(an insert_list)' %(len(outs), len(List)))

    
    outsHtml = ["<table border='1'>"]

    outsHtml0, alignHtml = lineup(outs, optD={'n':1, 'align':1})

    for iL in range(len(outsHtml0)):
        out0 = outsHtml0[iL]
        w = out0.split()

        # 0) Remove '|' entries, need not in html table
        for iw in range(len(w)-1, -1, -1):
            if w[iw] == '|':
                w.pop(iw)
                if iL == 0: alignHtml.pop(iw)
        

        # 1) Remove columns (if any)
        for col in sorted(remove_cols, reverse=True): 
            w.pop(col)
            if iL == 0: alignHtml.pop(col)


        # 2) Insert columns (if any)
        for icol in range(len(insert_cols)):
            col = insert_cols[icol]
            txt = insert_lists[icol][iL]

            w.insert(col,txt)
            if iL == 0: alignHtml.insert(col,'l')  # hardcoded


        # 3) Turn array back to text
        out = "<tr> "
        for iw in range(len(w)):
            zout = w[iw]
            if iL == 0: thtd = 'th'
            else: thtd = 'td'
            if alignHtml[iw] in ['l']: zalign = 'left'
            elif alignHtml[iw] in ['r']: zalign = 'right'
            elif alignHtml[iw] in ['c']: zalign = 'center'
            else:
                zalign = 'center'
                    
            out += "<%s align='%s'> %s <%s>" %(thtd, zalign, zout, thtd)

        out += " </tr>"
        outsHtml.append(out)


    outsHtml += ["</table>"]
    return outsHtml


# ###################################################################
def GetLinuxDist():
    ## distro is used from python3.8 >
    #import distro
    import platform as distro
    dist = 'RHEL5'  # default


    try: distT = distro.linux_distribution()
    except: distT = distro.dist()

    try: ver = float(distT[1])
    except:
        print('ERROR  GetLinuxDist   Version not a number: %s  [Returning dist=%s]' %(distT, dist))
        return dist
        

    if distT[0] in ['redhat','Red Hat Enterprise Linux Workstation']   or 'Red Hat' in distT[0]:
        if   5 <= ver < 6: dist = 'RHEL5' 
        elif 6 <= ver < 7: dist = 'RHEL6'  # 2014-08 
        elif 7 <= ver < 8: dist = 'RHEL7'
        elif 8 <= ver < 9: dist = 'RHEL8'
        else:
            print('ERROR  GetLinuxDist  Version not recognised: %s  [Returning dist=%s]' %(distT, dist))
            #return dist

    else:
        print('INFO  GetLinuxDist  Distribution not recognised: %s     (Need to add to GetLinuxDist())   [Returning dist=%s] ' %(distT, dist))
        #return dist
        
    return dist
    

# ###################################################################
# ###################################################################
