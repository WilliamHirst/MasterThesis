#!/usr/bin/env python

#########################################
# ArgReader
# b.k.gjelsten@fys.uio.no
# Sep 25, 2011: second version 
#########################################
# DESCRIPTION
#
#
# TO BE IMPLEMENTED / FIXED
#  - Option for case sensitivity
#  
#########################################


class ArgReader:
    
    # Nomenclature:
    #   args: options and values (all the text sent along)
    #   opts: options, the keywords, the codes tested against
    #   val(ue): the value following the keyword and to be returned

    def __init__(self, args, VB=0):
        self.args = list(args)
        # self.argsOrig = list(args)  # untouched (currently not in use)
        self.value = None
        self.exception = []
        self.lastOpt = None

        self.VB = VB


    # ########### Test if is empty, noexception, allok
    ### Test if is empty
    def IsEmpty(self,vb=1):  # Legacy
        status = 1
        if len(self.args) > 1:
            status = 0
            if vb: print('ERROR::ArgReader: Non-empty args: %s' %(self.args[1:]))
        return status

    ### Test if exceptions (fixed 2013-10-08 ; never worked earlier)
    def NoException(self,vb=1):  # Legacy
        status = 1
        if len(self.exception) > 0:
            status = 0
            if vb:
                for L in self.exception: print(L)
        return status


    ### Test if all ok
    def AllOk(self,vb=1):
        return self.NoException(vb) * self.IsEmpty(vb)


    # ############ New test. This allows to get error messages after 
    ### 
    def NotEmpty(self):
        status = []
        if len(self.args) > 1: status = ['ERROR::ArgReader: Non-empty args: %s' %(self.args[1:])]
        return status

    ###
    def HaveException(self):
        status = []
        if len(self.exception) > 0:
            status = 0
            if vb:
                for L in self.exception: status.append(L)
        return status
        
    ### Test if all ok
    def ErrorMessages(self):
        status  = self.NotEmpty()
        status += self.HaveException()
        return status




    # ########### Methods for retrieving values

    def val(self): return self.value

    def valI(self):
        try:
            return int(self.value)
        except:
            self.exception.append('Exception: Last option: %s   Value valI not a number: %s' %(self.lastOpt, self.value))
            print(self.exception[-1])
            return None
        
    

    def valF(self):
        try: return float(self.value)
        except: 
            self.exception.append('Exception: Last option: %s   Value valF not a number: %s' %(self.lastOpt, self.value))
            print(self.exception[-1])
            return None


    def list(self,delim=','):
        if delim == '': L = self.value.split()
        else: L = self.value.split(delim)
        return L


    def listI(self,delim=','):
        L = self.list(delim)
        L2 = []
        for i in range(len(L)):
            try: L2.append(int(L[i]))
            except:
                L2.append(None)
                self.exception.append('Exception: Last option: %s   listI: %ith list element not a number: %s' %(self.lastOpt, i, L[i]))
                print(self.exception[-1])
        return L2


    def listIif(self,delim=','):
        L = self.list(delim)
        L2 = []
        for i in range(len(L)):
            try: L2.append(int(L[i]))
            except: L2.append(L[i])  # returns the string if not int
        return L2


    def listF(self,delim=','):
        L = self.list(delim)
        L2 = []
        for i in range(len(L)):
            try: L2.append(float(L[i]))
            except:
                L2.append(None)
                self.exception.append('Exception: Last option: %s   listF: %ith list element not a number: %s' %(self.lastOpt, i, L[i]))
                print(self.exception[-1])
        return L2


    def listFif(self,delim=','):
        L = self.list(delim)
        L2 = []
        for i in range(len(L)):
            try: L2.append(float(L[i]))
            except: L2.append(L[i])  # returns the string if not float
        return L2


    # ##########
    def AddArg(self,toadd):
        # used to hack / tamper with arguments while checking
        if type(toadd) is list: self.args += toadd
        else: self.args.append(toadd)


    # ########## No value retrieval
    def has(self,opt):
        if self.VB: print('DEBUG: .has:  args=',self.args, '   opt=',opt)
        
        if type(opt) == list: opts = opt
        else: opts = [opt]

        toremove = []
        status = 0

        # check if opt in arg and toremove and return status
        for iarg in range(len(self.args)):
            arg = self.args[iarg]
            
            for iopt in range(len(opts)):
                opt = opts[iopt]
                
                if opt == arg: 
                    # mark arg for removal
                    
                    toremove.append(iarg)
                    # self.args.toremove(opt)
                    # update status, but do not break off code if opt found,
                    #                but rather continue to toremove others if there
                    status += 1
                    self.lastOpt = arg

        # toremove.sort()  # is already sorted
        toremove.reverse()
        for i in toremove:
            self.args.pop(i)

        if self.VB: print('DEBUG: .has:  args=',self.args, '   (after test)')

        return status



    # ########## With value retrieval
    def hasget(self,opt):
        if self.VB: print('DEBUG: .has:  args=',self.args, '   opt=',opt)
        
        if type(opt) == list: opts = opt
        else: opts = [opt]

        toremove = []
        status = 0
        self.value = None
        
        for iarg in range(len(self.args)):
            arg = self.args[iarg]
            if iarg == len(self.args) - 1: continue # cannot be last opt since need a return value

            for iopt in range(len(opts)):
                opt = opts[iopt]
                
                if opt == arg: 
            
                    toremove.append(iarg)
                    toremove.append(iarg+1)
                    self.value = self.args[iarg+1]
                    status += 1
                    self.lastOpt = arg
                    break  # here need to break (don't want to collect several values)  

                
        toremove.reverse()
        for i in toremove:
            self.args.pop(i)

        if self.VB: print('DEBUG: .has:  args=',self.args, '   (after test)')
        return status

        
# ##########
    

