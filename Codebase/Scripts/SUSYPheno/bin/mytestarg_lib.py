import random
VERBARG=0

def testargget(zArgv,ztest,zget):
    #  1: test string found, zArgv and zget changed
    #  0: test string not found, zArgv and zget untouched
    # -2: test string found but not enough arguments in argv, zArgv, zget untouched [warn]

    #(Reusing code from testget(), so code may be a bit bulky

    for ii in range(len(zget)): zget.pop()  #Fix jan20: this erases whatever was in the list ...
    znarg=1 #len(args)
    znArgv=len(zArgv)

    stat=0
    for iHit in range(znArgv):
        if(zArgv[iHit].lower()==ztest):
            stat=1

            if(znArgv-iHit-1<znarg):
                print("myWARN: testarg(): only "+str(znArgv-iHit+1)+" arguments following "+ztest+", required is "+str(znarg-2))
                stat=-2
                return stat

            # Have enough args: get variable 
            zget.append(zArgv[iHit+1])  #need to use array to get back values
            #zget=zArgv[iHit+1]
                    
            # Strip global variable (if arrive here, stat==1, so don't have to test)
            if(VERBARG>=5):
                print("testarg: argv before: ", end=' ')
                print(zArgv)
                
            zdel=zArgv[iHit]
            del(zArgv[iHit])
            zdel+=" "+zArgv[iHit]
            del(zArgv[iHit]) 

            if(VERBARG>=4):
                print("testarg: argv after:  ", end=' ')
                print(zArgv)
            if(VERBARG>=2): print("testarg: argv remove:  '"+zdel+"'")
 
            return stat

    # If arrive here, argument was not found
    return stat #(stat==0)

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

 
# testargget2 is more elegant than testargget, not requireing the dumA-structure
# however it is not yet in use (as of March 12, 2010)
#
def testargget2(zArgv,ztest):  
    #  1: test string found, zArgv and zget changed
    #  0: test string not found, zArgv and zget untouched
    # -2: test string found but not enough arguments in argv, zArgv, zget untouched [warn]

    #(Reusing code from testget(), so code may be a bit bulky

    ###for ii in range(len(zget)): zget.pop()  #Fix jan20: this erases whatever was in the list ...
    znarg=1 #len(args)
    znArgv=len(zArgv)
    zget = -999

    stat=0
    for iHit in range(znArgv):
        if(zArgv[iHit].lower()==ztest):
            stat=1

            if(znArgv-iHit-1<znarg):
                print("myWARN: testarg(): only "+str(znArgv-iHit+1)+" arguments following "+ztest+", required is "+str(znarg-2))
                stat=-2
                return stat, zget

            # Have enough args: get variable 
            ###zget.append(zArgv[iHit+1])  #need to use array to get back values
            #zget=zArgv[iHit+1]
            zget = zArgv[iHit+1]
                    
            # Strip global variable (if arrive here, stat==1, so don't have to test)
            if(VERBARG>=5):
                print("testarg: argv before: ", end=' ')
                print(zArgv)
                
            zdel=zArgv[iHit]
            del(zArgv[iHit])
            zdel+=" "+zArgv[iHit]
            del(zArgv[iHit]) 

            if(VERBARG>=4):
                print("testarg: argv after:  ", end=' ')
                print(zArgv)
            if(VERBARG>=2): print("testarg: argv remove:  '"+zdel+"'")
 
            return stat, zget

    # If arrive here, argument was not found
    return stat, zget #(stat==0)

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 


def testarg(zArgv,ztest,*args):
    # RETURN VALUES:
    #  1: test string found, zArgv stripped
    #  0: test string not found, zArgv untouched
    # -1: test string found but additional arguments not found, zArgv untouched 
    # -2: test string found but not enough arguments in argv, zArgv untouched [warn]
    
    znarg=len(args)
    znArgv=len(zArgv)

    stat=0
    for iHit in range(znArgv):
        if(zArgv[iHit].lower()==ztest):
            stat=1

            if(znArgv-iHit-1<znarg):
                print("myWARN: testarg(): only "+str(znArgv-iHit+1)+" arguments following "+ztest+", required is "+str(znarg-2))
                stat=-2
                return stat

            # Have enough args: continue test 
            for ia in range(znarg):
                if(zArgv[iHit+1+ia]!=args[ia]):
                    stat=-1
                    return stat
                    
            # Strip global variable (if arrive here, stat==1, so don't have to test)
            if(VERBARG>=5):
                print("testarg: argv before: ", end=' ')
                print(zArgv)
                
            zdel=zArgv[iHit]
            del(zArgv[iHit])
            #del(zArgv[iHit])
            for ia in range(znarg):
                zdel+=" "+zArgv[iHit]
                del(zArgv[iHit]) 
            if(VERBARG>=4):
                print("testarg: argv after:  ", end=' ')
                print(zArgv)
            if(VERBARG>=2): print("testarg: argv remove:  '"+zdel+"'")
 
            return stat

    # If arrive here, argument was not found
    return stat #(stat==0)

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

def gettimestring(dec=0, rndm=0):
    # returns string yyyymmdd_hhmmss
    from time import localtime as mylocaltime
    t = mylocaltime()
    timestring =  str(t[0])+str(t[1]).zfill(2)+str(t[2]).zfill(2)
    timestring += "_"+str(t[3]).zfill(2)+str(t[4]).zfill(2)+str(t[5]).zfill(2)
    del mylocaltime

    # adhoc to get more detailed than seconds (if needed)
    if dec > 0:
        from time import time as mytimez
        t2 = mytimez()
        s1 = "%.7f" %(t2)
        s2 = s1.split(".").pop()
        s3 = s2[0:min(7,dec)]
        timestring += "_%s" %(s3)
        del mytimez

    if rndm: timestring += getranfilen('_',rndm)  # add a random string (sometimes needed)
    
    return timestring

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

def getranfilen(base,n):
    return base+str(int(random.uniform(0,10**n))).zfill(n)

def randomsleep(t):
    from time import sleep as mysleep
    t = mysleep(t * random.random())
                   
def insertZero(s1,wid1,wid2):
    #(only first entry is replaced, more rigorous implementation would test on occurences)
    # insertZero("aa_p4",1,2) -> "aa_p04"
    # insertZero("aa_p4",1,3) -> "aa_p004"
    # insertZero("aa_44003_3bc",1,2) -> "aa_44003_03bc"
    
    n=len(s1)
    i=0
    while(i<n-wid1+1): 
        # Find a digit
        if(s1[i].isdigit()):
            # Check how many digits follow this one
            j=i+1
            while(j<n and s1[j].isdigit()): j+=1
            
            #print s1[i:j]
            # Test if acceptable region (will allow for more general usage)
            if(j-i>=wid1 and j-i<=wid2):
                #then replace
                s2=s1[0:i]+s1[i:j].zfill(wid2)+s1[j:n]
                return s2
            else:
                i=j
                continue
        i+=1
        
    #if arrive here, then no digit within [wid1,wid2] is found
    return s1


def minNotZero(a,b):
    if(max(a,b)<=0): return 0;
    if(a==0): return b
    if(b==0): return a
    return min(a,b)

def hline(marker,n):
    nIter=n/len(marker)+1
    zstr=""
    for ii in range(nIter):
       zstr+=marker
    return zstr[0:n]

def namecut(name,wid,mode,add):
    n=len(name)
    if(len(name)<=wid): return name
    if(mode=="cutR"): return name[0:wid-len(add)]+add
    if(mode=="cutL"): return add+name[n-wid+len(add):n]
    print("non-understood mode: "+mode)
    return name

def autoMan():
    print(50*"-")
    import os,sys
    os.system("cat %s |grep 'if testarg\|if(testarg\|if (testarg'" %(sys.argv[0]))
    print(50*"-")
