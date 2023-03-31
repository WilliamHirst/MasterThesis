from array import array

#from ROOT import TColor,gStyle
#from ROOT import TDatime

#from ROOT import *    # 2011-11-11: import all  # This causes problem with using '-h'
import ROOT    # 2013-11-29: changed from import * [due to problem above]
#from ROOT import TColor
from math import fmod,ceil
import random

# =============================================================================
def set_palette(name='coolhot', ncontours=999):
    # http://ultrahigh.org/2007/08/20/making-pretty-root-color-palettes/
    """Set a color palette from a given RGB list
    stops, red, green and blue should all be lists of the same length
    see set_decent_colors for an example"""

    # NB: leftmost column defines start value (red[0],green[0],blue[0] is first colour)


    if name == "gray" or name == "grayscale":
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [1.00, 0.84, 0.61, 0.34, 0.00]
        green = [1.00, 0.84, 0.61, 0.34, 0.00]
        blue  = [1.00, 0.84, 0.61, 0.34, 0.00]

    elif name == 'myred':
        stops = [0.00, 1.00]
        red   = [1.00, 1.00]
        green = [0.00, 0.00]
        blue  = [1.00, 0.00]

    elif name == 'coolhot1':
        # default palette, looks cool
        stops = [0.00, 0.20, 0.50, 0.95, 1.00]  #sets the 'length' of each color 

        red   = [0.00, 0.00, 0.87, 1.00, 1.00]  #here five colours defined, the first has (r,g,b)=(0,0,1.), etc.
        green = [0.00, 0.81, 1.00, 0.20, 0.00]
        blue  = [1.00, 1.00, 0.12, 0.00, 0.00]

    elif name == 'coolhot1b': # not satisfied but replaces coolhot1 as default for now (to see text in zero)
        # default palette, looks cool
        stops = [0.00, 0.38, 0.93, 1.00]

        red   = [0.00, 0.87, 1.00, 1.00]
        green = [0.30, 1.00, 0.20, 0.00]
        blue  = [1.00, 0.12, 0.00, 0.00]

    elif name == 'coolhot1c': # soft blue at bottom: easiy to read the loest
        stops = [0.00, 0.38, 0.93, 1.00]

        red   = [0.00, 0.87, 1.00, 1.00]
        green = [0.70, 1.00, 0.20, 0.00]
        blue  = [1.00, 0.12, 0.00, 0.00]


    elif name == 'coolhot2':
        stops = [0.00, 0.20, 0.95, 1.00]

        red   = [0.00, 0.00, 1.00, 1.00]
        green = [0.00, 0.81, 0.20, 0.00]
        blue  = [1.00, 1.00, 0.00, 0.00]

    elif name == 'coolhot3':
        stops = [0.00, 0.50, 1.00]

        red   = [0.00, 1.00, 1.00]
        green = [0.00, 1.00, 0.00]
        blue  = [1.00, 0.22, 0.00]

    elif name == 'coolhot4':
        stops = [0.00, 1.00]  

        red   = [0.00, 1.00]  
        green = [0.00, 0.00]
        blue  = [1.00, 0.00]

    elif name == 'coolhot5':
        stops = [0.00, 1.00]  

        red   = [1.00, 1.00]  
        green = [1.00, 0.00]
        blue  = [0.20, 0.00]

    else:
        # default palette, looks cool
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [0.00, 0.00, 0.87, 1.00, 0.51]
        green = [0.00, 0.81, 1.00, 0.20, 0.00]
        blue  = [0.51, 1.00, 0.12, 0.00, 0.00]


    s = array('d', stops)
    r = array('d', red)
    g = array('d', green)
    b = array('d', blue)

    npoints = len(s)
    ROOT.TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)   # <--- This one cannot be issued many times: Slows down more and more. Why? 
    #return 0
    ROOT.gStyle.SetNumberContours(ncontours)


# =============================================================================
def TDatimeFromStamp(stamp):
  # 1. Dechiffer the stamp
  # assume here stamp has format of type 2011-12-12_12:14:13 (mode0) (and h, m and s as well)
  W = stamp.replace('-',' ').replace('_',' ').replace(':',' ').replace('h',' ').replace('m',' ').replace('s','').split()
  
  z = []
  nw = len(W)
  for iw in range(6):
    if iw<nw: z.append(int(W[iw]))
    else: z.append(0)

  # 3. the result
  dt0 = TDatime(z[0],z[1],z[2],z[3],z[4],z[5])
  return dt0


  
# =============================================================================
def ColourRatio2D(hs, hsInfo=[], opt={}, VB=1, meta=0, fillmode=2): 
    # fillmodes: 1: BottomLeft, horisontal
    #            2: TopLeft, horisontal (new default)

    DEBUG = 0  # can set for some more printout

    zmeta = {}
    zopt = dict(opt)  # make a copy
    zhsInfo = list(hsInfo)
    
    # ---------------------- CHECKS
    nhs = len(hs)
    if nhs == 0:
        print("No histograms given")
        return 

    # ---------------------- INIT

    # --- hsInfo



    # --- opt
    
    #max = ['sum']
    if 'max' not in zopt:
        zopt['max'] = 'sum'
        if VB: print("INFO::ColorRatio2D: setting 'max' to 'sum'")
        
    plotmodes = ['random1','chrono']
    if 'plotmode' not in zopt or zopt['plotmode'] not in plotmodes:
        zopt['plotmode'] = 'random1'
        if VB: print("INFO::ColorRatio2D: setting 'plotmode' to '%s'" %(zopt['plotmode']))

    tag = str(int(random.uniform(10000,99999)))

    if 'hnam' not in zopt: #  or 1:
        zopt['hnam'] = 'hCR_%s' %(tag)
        if VB>1: print("INFO::ColorRatio2D: setting 'hnam' to %s" %(zopt['hnam']))

    if 'htit' not in zopt:
        zopt['htit'] = zopt['hnam']
        if VB: print("INFO::ColorRatio2D: setting 'htit' to %s" %(zopt['htit']))

    if 'xx' not in zopt:
        zopt['xx'] = 3
        if VB: print("INFO::ColorRatio2D: setting 'xx' to %i" %(zopt['xx']))

    if 'yy' not in opt:
        zopt['yy'] = 3
        if VB: print("INFO::ColorRatio2D: setting 'yy' to %i" %(zopt['yy']))
        
        
    xx = zopt['xx']
    yy = zopt['yy']
    hnam = zopt['hnam']
    htit = zopt['htit']
    
    
    h0 = hs[0]
    tax = h0.GetXaxis()
    tay = h0.GetYaxis()
    x = h0.GetNbinsX()
    y = h0.GetNbinsY()


    # Make (flexible) binning:
    binsX = array('f')
    for ix in range(x):
        for ixx in range(xx):
            binsX.append( tax.GetBinLowEdge(ix+1) + ixx*1.*tax.GetBinWidth(ix+1)/xx )
    binsX.append(tax.GetXmax())
    binsY = array('f')
    for iy in range(y):
        for iyy in range(yy):
            binsY.append( tay.GetBinLowEdge(iy+1) + iyy*1.*tay.GetBinWidth(iy+1)/yy )
    binsY.append(tay.GetXmax())

    #print 'binsX: ',binsX
    #print 'binsY: ',binsY
    
    #xmin = h0.GetXaxis().GetXmin()
    #xmax = h0.GetXaxis().GetXmax()
    #ymin = h0.GetYaxis().GetXmin()
    #ymax = h0.GetYaxis().GetXmax()
    
    X = x*xx
    Y = y*yy

    #print 'Now init %s  (%s)' %(hnam, tag)
    #hCR = ROOT.TH2D(hnam, htit, X, xmin, xmax, Y, xmin, xmax)
    hCR = ROOT.TH2D(hnam, htit, X, binsX, Y, binsY)

    tax = hCR.GetXaxis()
    tay = hCR.GetYaxis()

    #print 'DEBUG kilelib_ROOT: ', zopt.get('xtit','xtit'), '   ', zopt.get('ytit','ytit')  axis title set here !!!
    if 'xtit' in zopt: tax.SetTitle(zopt['xtit'])
    if 'ytit' in zopt: tay.SetTitle(zopt['ytit'])
    tax.CenterTitle()
    tay.CenterTitle()
    tax.SetTitleOffset(zopt.get('xoff',1.15))  # 2012-10-10
    #tay.SetTitleOffset(zopt.get('xoff',1.4))   # 2012-10-10   
    tay.SetTitleOffset(zopt.get('yoff',1.4))   # 2012-10-10    # 2013-07-16 (untested)
    


    # make box
    box = xx*yy*[0]
    #for iyy in range(yy): box.append(xx*[0])

    # make ninbox[nhs] and rest[nhs]
    ninbox = nhs*[0]
    remain = nhs*[0]

    ninboxAllbins = nhs*[0]
    ninboxAllbinsNone = 0

    nzerosum = 0
    nzerosum_intermedbin = 0
    nzerosum_real = 0

    # ---------------------- LOOP
    for ix in range(x):
        binx = ix+1
        # if ix not in [0]: continue
        for iy in range(y):
            biny = iy+1
            # if iy not in [14]: continue

            isintermediatebin = ( (ix/2)*2 == ix or (iy/2)*2 == iy )  # <-- maybe fragile (if we enter without intermediate bins, or if they somehow follow a different structure)

            # init the box
            for ib in range(xx*yy): box[ib] = 0
            # for ib in range(xx*yy): remain[ib] = 0


            # a) Find max for this (ix,iy)
            thissum = 0.
            for ih in range(nhs):
                thissum += hs[ih].GetBinContent(binx,biny)

            if DEBUG:
                out = 'DEBUG: ix:%i  iy:%i  sum = %6.4f | ' %(ix,iy, thissum)
                for ih in range(nhs):
                    out += "%6.4f  " %(hs[ih].GetBinContent(binx,biny))
                print(out)
                    
            if thissum == 0:
                nzerosum += 1
                if isintermediatebin:
                    nzerosum_intermedbin += 1
                    continue
                else: 
                    nzerosum_real += 1
                    out = "kilelib_ROOT::ColourRatio2D  WARNING: zero sum for ix=%2i  iy=%2i, (no: skipping this bin) [THIS IS NOT AN INTERMEDIATE BIN]" %(ix,iy)
                    if VB and nzerosum_real < 3: print(out)
                    if VB and nzerosum_real == 3: print(out + '   [SUPPRESSING FURTHER WARNINGS]')
                    #if VB>1: print "WARNING: zero sum for ix=%2i  iy=%2i, (no: skipping this bin) [could well be an 'intermediate-bin' (between points)]" %(ix,iy)
                    # thishmax = 1.
                    

            if zopt['max'] == 'sum':
                thishmax = thissum
                if DEBUG: print('DEBUG using sum as max')
            else:
                thishmax = float(zopt['max'])  # use a constant
                if DEBUG: print('DEBUG taking max from opt: %.3f' %(float(zopt['max'])))

            # b) Fill#1 according to hist & max
            ninboxNone = 0
            ninboxTot = 0
            for ih in range(nhs):
                if thishmax == 0.:
                    if DEBUG: print('DEBUG zero thishmax') 
                    z0 = 0.
                else: 
                    z0 = hs[ih].GetBinContent(binx,biny) / thishmax

                z = z0 * xx*yy  # float
                ninbox[ih] = int(z)
                remain[ih] = z-int(z)
                ninboxTot += ninbox[ih]
                #if 1 and (binx==11 and biny==14) or (binx==9 and biny==9): #should be MU=240 and M2=300  or MU=M2=200
                if 0 and (binx==1 and biny==3): 
                    print('DBUG %-10s  %6.4f  ninbox[%i]=%3i  %6.4f' %(zhsInfo[ih]['labelT'], z/(xx*yy), ih, ninbox[ih], remain[ih]))

            if DEBUG: print('DEBUG  kilelib_ROOT: ninboxTot=%i' %(ninboxTot))

            # c) Fill#2: Now deal with the remains (utjamningsmandat) in random fashion (these can (and will) also land on empty if sum < 1)
            iIter = 0
            while ninboxTot < xx*yy:
                iIter += 1
                # i) First sum the remain
                remainTot = 0.
                for ih in range(nhs): remainTot += remain[ih]

                # ii) Get a random from a uniform distribution in the appropriate range
                rndm = random.uniform(0,remainTot)
                # print 'DEBUG: ninboxTot= %3i  random(0,%.3f): %.3f' %(ninboxTot, remainTot, rndm) 
                
                # iii) Then check the remain again and see who won
                remainSumSoFar = 0.
                hit = 0
                for ih in range(nhs):
                    remainSumSoFar += remain[ih]
                    # print 'DEBUG  rndm=%.3f  vs  remainSumSoFar:%.3f' %(rndm, remainSumSoFar)
                    if rndm < remainSumSoFar:
                        # first hit is the correct one
                        ninbox[ih] += 1
                        ninboxTot += 1
                        remain[ih] = 0
                        #print 'DEBUG hit ih = %i' %(ih)
                        hit = 1
                        break

                # if arrive here, then rndm has landed in the slot of missing
                if not hit: 
                    ninboxNone += 1
                    ninboxTot += 1
                
            
                # DEBUG (should never fire)
                if iIter > 10000:
                    print('BUG in Fill#2 (the remains)')
                    return 

            if DEBUG: print('DEBUG  kilelib_ROOT::ColourRatio2D ninboxNone + sum vs ninboxTot vs xx*yy: %5i + %5i vs %5i vs %5i' %(ninboxNone, sum(ninbox), ninboxTot, xx*yy))

            # d) uptdate ninboxAllbins (possibly used in labeling)
            for ih in range(nhs):
                #print 'DEBUG  end ninbox[%i] = %3i' %(ih, ninbox[ih])
                ninboxAllbins[ih] += ninbox[ih]
            ninboxAllbinsNone += ninboxNone


            # e) Fill box according to plotmode
            if zopt['plotmode'] == 'chrono':
                nbin = 0
                for ih in range(nhs):
                    if DEBUG: print('DOBUG: i=%i  ninbox[%i]=%i' %(ih,ih,ninbox[ih]))
                    for i in range(ninbox[ih]):
                        if nbin+i >= len(box): print("FAILS  kilelib_ROOT::ColourRatio2D  (ix,iy)=(%i,%i)  nhs:%i, nbin:%i, len(box):%i   i:%i" %(ix,iy,nhs, nbin, len(box), i))
                        if nbin+i >= len(box):
                            print("FATAL::kilelib_ROOT::ColourRatio2D  nbin+i (%i) >= len(box) %i. This typically happens if optmax is not set to 'sum' in colourratio_ex1.py and the sum of histos (for this bin) exceeds 1 (i.e. if histos are not already ratios) " %(nbin+i, len(box)))
                        box[nbin+i] = ih+1        # <-- colour
                    nbin += ninbox[ih]
                    
            elif zopt['plotmode'] == 'random1':
                for ih in range(nhs):
                    # Go randomly through the box until an empty slot is found
                    iIter = 0
                    ihputinbox = 0
                    while ihputinbox < ninbox[ih]: 
                        rndm = random.randint(0,xx*yy-1)
                        if box[rndm] == 0:
                            box[rndm] = ih+1      # <-- colour
                            ihputinbox += 1
                            if VB>1: print('debug box=%s  %i / %i (ih=%i) (%i,%i)' %(box, ihputinbox, ninbox[ih], ih, ix, iy))
                        iIter += 1
                        if VB>1: print('debug box=%s  %i / %i (ih=%i) (%i,%i)  <-' %(box, ihputinbox, ninbox[ih], ih, ix, iy))
                        if iIter > 10000:
                            print('BUG: plotmode random: never finds empty slot (ix=%2i iy=%2i  ih=%i, rndm=%2i, ihputinbox=%i, ninbox[ih]=%i, box=%s)' %(ix, iy, ih, rndm, ihputinbox,ninbox[ih],box))
                            return


            # f) Put box in histogram
            #    (ix and iy are the ordinary/big bins
            if fillmode == 1: # from bottom left, horisontal
                for ib in range(xx*yy):
                    iyy = ib/xx
                    ixx = ib - iyy*xx
                    
                    iX = ix*xx + ixx
                    iY = iy*yy + iyy

                    #print 'ib,  iyy, ixx = %5i  %3i %3i |  iX: %5i  iY: %5i' %(ib,iyy,ixx, iX, iY)
                    hCR.SetBinContent(iX+1, iY+1, box[ib])

            elif fillmode == 2:  # from top left, horisontal  # new default since it fits with the legend order
                for ib in range(xx*yy):
                    iyy = ib/xx
                    ixx = ib - iyy*xx

                    iX = ix*xx + ixx
                    iY = iy*yy + yy-1-iyy

                    #print 'ib,  iyy, ixx = %5i  %3i %3i |  iX: %5i  iY: %5i' %(ib,iyy,ixx, iX, iY)
                    hCR.SetBinContent(iX+1, iY+1, box[ib])

                    if ib >= sum(ninbox):
                        if VB>1: print('here ib: %3i   box[ib]: %3i' %(ib, box[ib]))
                        hCR.SetBinContent(iX+1, iY+1, 1e-5)
                

    if VB:
        nxReal = (x+1)/2; nyReal = (y+1)/2
        if VB>0 or x*y-nxReal*nyReal != nzerosum:
            print('INFO::Colourratio  (x,xx,X):(%i,%i,%i)  (y,yy,Y)=(%i,%i,%i)  totbin: %i   totrealbin: %i    diff: %i   nzerosum: %i  ' %(x,xx,X, y,yy,Y, x*y, nxReal*nyReal, x*y-nxReal*nyReal, nzerosum), end=' ') 
            if x*y-nxReal*nyReal == nzerosum: print('[EQUAL]')
            else: print('[NOT EQUAL!!]')


    # Construct meta (which might be sent back)
    zmeta['ninboxAllbins'] = ninboxAllbins  # can be used to show/not show Legend
    zmeta['ninboxAllbinsNone'] = ninboxAllbinsNone


    # --- loop done
    if not meta: return hCR

    return hCR, meta


# =============================================================================
def DrawLine(x1,y1,x2,y2,opt=[]):  # doesn't work that splendidly -- have to return line for it to draw
    if 'col' in opt: col = opt['col']
    else: col = ROOT.kBlack
    if 'sty' in opt: sty = opt['sty']
    else: sty = 1
    if 'wid' in opt: wid = opt['wid']
    else: wid = 1
    
    line = ROOT.TLine(x1,y1,x2,y2)
    line.SetLineColor(col)
    line.SetLineStyle(sty)
    line.SetLineWidth(wid)
    line.Draw()
    return line

# =============================================================================
def VerticalLines(L, y1, y2, xs):
    for x in xs: L.DrawLine(x,y1,x,y2)
# =============================================================================
def HorisontalLines(L, x1, x2, ys):
    for y in ys: L.DrawLine(x1,y,x2,y)
# =============================================================================
def Getbins(axis, opt=[]):
    bins = array('f')
    for i in range(axis.GetNbins()):
        bins.append(axis.GetBinLowEdge(i+1))
    bins.append(axis.GetXmax())
    if 'drop0' in opt and bins: bins.pop(0)
    if 'dropN' in opt and bins: bins.pop()
    return bins
    
# =============================================================================
def DrawHistGrid(L, h, col=1, sty=2, wid=1):
    # Also nice: DrawHistGrid(gLine, h, col=0, sty=1)  # abit more discrete
    L.SetLineColor(col)
    L.SetLineStyle(sty)
    L.SetLineWidth(wid)
    
    xbins = Getbins(h.GetXaxis(),['drop0','dropN'])
    ybins = Getbins(h.GetYaxis(),['drop0','dropN'])

    xmin = h.GetXaxis().GetXmin()
    xmax = h.GetXaxis().GetXmax()
    ymin = h.GetYaxis().GetXmin()
    ymax = h.GetYaxis().GetXmax()
    
    VerticalLines(L,ymin,ymax, xbins)
    HorisontalLines(L,xmin,xmax, ybins)

# =============================================================================

def FlipAxes(h): 
    # create new histo
    hn = h.GetName() + '_flip'
    htit = h.GetTitle()
    # (x and y below refers to h-histo)
    xbins = Getbins(h.GetXaxis())
    ybins = Getbins(h.GetYaxis())
    nx = len(xbins)-1
    ny = len(ybins)-1

    h2 = ROOT.TH2F(hn,htit,ny,ybins,nx,xbins)  # flips x and y

    h2.GetXaxis().ImportAttributes(h.GetYaxis())  # flips x and y
    h2.GetYaxis().ImportAttributes(h.GetXaxis())  # flips x and y 

    
    for ix in range(nx+2):   # includes under- and overflow bin 
        for iy in range(ny+2):    # includes under- and overflow bin 
            h2.SetBinContent(iy,ix, h.GetBinContent(ix,iy))  # flips x and y
            h2.SetBinError(iy,ix, h.GetBinError(ix,iy))  # flips x and y

    return h2

# =============================================================================

def GetNoncontBins(val, width, VB=0): 
    binsL = {}
    for var in list(val.keys()): 
        binsL[var] = []
        for v in val[var]: 
            previous = -99999
            if binsL[var]: previous = binsL[var][-1]
            thislower = v - 0.5*width
            if previous != thislower: binsL[var].append(v - 0.5*width) # conditional because can overlap prevous
            if previous > thislower: 
                print('FATAL: GetNoncontBins: the bins are overlapping (maybe related to binwidth..), returning []')
                return []
            binsL[var].append(v + 0.5*width)

    if VB: 
        for key in list(binsL.keys()): 
            print('%s: %s' %(key, binsL[key]))

    return binsL


# =============================================================================
# ============================================================================= Below: methods for Inter/Extrapolation
# =============================================================================


def FirstNonZeroBinI(h,ibinx,ibiny,step,xory,nret=1,nstepmax=-1,VB=0):
    # xory: 'x' or 'y': horisontal or vertical search
    # nret: how many values to return: for interpolation, one is usual
    #       for extrapolation, two is needed
    # nstepmax: max number of bins away from the bin in question (usually 1 will do). -1:off

    if xory == 'x': 
        ibMax = h.GetNbinsX()
        ib = ibinx
    else: 
        ibMax = h.GetNbinsY()
        ib = ibiny
    ibMin = 1
    okbins = []
    nstep = 0
    while ib >= ibMin and ib <= ibMax: 
        if xory == 'x': z = h.GetBinContent(ib,ibiny)
        else: z = h.GetBinContent(ibinx,ib)

        if z != 0: okbins.append(ib)

        if len(okbins) >= nret: break
        nstep += 1
        if nstep > nstepmax > -1: break  # tests also for nstepmax > -1 (on/off)

        ib += step
            
    if len(okbins) < nret: return -1
    elif nret > 1: return okbins[0:nret]  # returns an array
    else: return okbins[0]  # returns a value


def PolateBin(h,ibinx,ibiny,xory,ib1,ib2): 
    #print 'ib1=%i ib2=%i' %(ib1,ib2)
    if xory == 'x': 
        tax = h.GetXaxis()
        x = tax.GetBinCenter(ibinx)
        zX1 = h.GetBinContent(ib1, ibiny)
        zX2 = h.GetBinContent(ib2, ibiny)
        xX1 = tax.GetBinCenter(ib1)
        xX2 = tax.GetBinCenter(ib2)
        rateX = (zX2-zX1)/(1.*xX2-xX1)
        delX = x-xX1
        z = zX1 + rateX * (x-xX1)

    if xory == 'y': 
        tay = h.GetYaxis()
        y = tay.GetBinCenter(ibiny)
        zY1 = h.GetBinContent(ibinx, ib1)
        zY2 = h.GetBinContent(ibinx, ib2)
        yY1 = tay.GetBinCenter(ib1)
        yY2 = tay.GetBinCenter(ib2)
        rateY = (zY2-zY1)/(1.*yY2-yY1)
        delY = y-yY1
        z = zY1 + rateY * (y-yY1)
    
    return z


def biniContainsbinxs(ax,ibin,binxs): 
    Low = ax.GetBinLowEdge(ibin)
    High = ax.GetBinUpEdge(ibin)
    hit = 0
    for binx in binxs: 
        if Low < binx < High: 
            return 1
    return 0
    
    
def GetBinContentXorY(h, ibinx, ibiny, ibinreplace, xory): 
    if xory == 'x': return h.GetBinContent(ibinreplace, ibiny)
    if xory == 'y': return h.GetBinContent(ibinx, ibinreplace)


def LinearPolation(h,mode=['inter','extra'],nstepmax=1,bins2polateX=[],bins2polateY=[], hmode='new', col=1, VB=1, zpolmax=None, zpolmin=None, zpolmaxmode='', zpolminmode='', xorys=['x','y'], steps=[-1,1]): 
    # hmode: 'old':fill the old; 'new':new(clone) and fill; 'newdiff':new and fill only extra
    # zpolmin(max)mode: if outside zpolmin(max): '': use the closest value  'absmin(max)': use the zpolmin(max) value

    if hmode in ['new','newdiff']: 
        hn = h.GetName() + '_polate'
        h2 = h.Clone(hn)
        if hmode == 'newdiff': h2.Reset()
        h2.SetMarkerColor(col)
        h2.SetLineColor(col)
    elif hmode in ['old']: 
        h2 = h
    else: 
        print("Warning: LinearPolation: Non-allowed hmode: %s   (USING 'new')" %(hmode))
        hmode = 'new'

    tax = h.GetXaxis()
    tay = h.GetYaxis()
    nbinx = h.GetNbinsX()
    nbiny = h.GetNbinsY()

    skippedX = []
    skippedY = []

    # Interpolate
    if 'inter' in mode: 
        for ibinx in range(1,nbinx+1): 

            binx = tax.GetBinCenter(ibinx)
            # ---------------------------- pinpoint exactly which bins to polate
            if bins2polateX and not biniContainsbinxs(tax, ibinx, bins2polateX): 
                skippedX.append((ibinx,binx))
                continue

            for ibiny in range(1,nbiny+1): 

                biny = tay.GetBinCenter(ibiny)
                # ---------------------------- pinpoint exactly which bins to polate
                if bins2polateY and not biniContainsbinxs(tay, ibiny, bins2polateY): 
                    skippedY.append((ibiny,biny))
                    continue

                bintag = "(%2i,%2i) = (%4.0f,%4.0f)" %(ibinx,ibiny, binx,biny)


                z = h.GetBinContent(ibinx,ibiny)
                if z != 0: continue
                nAv = 0
                zInt = 0.


                # horisontal
                if ibinx not in (1,nbinx): 
                    ibinx1 = FirstNonZeroBinI(h=h,ibinx=ibinx,ibiny=ibiny,step=-1,xory='x',nret=1,nstepmax=-1)
                    ibinx2 = FirstNonZeroBinI(h=h,ibinx=ibinx,ibiny=ibiny,step= 1,xory='x',nret=1,nstepmax=-1)
                    if ibinx1 > 0 and ibinx2 > 0: 
                        zX = PolateBin(h,ibinx,ibiny,xory='x',ib1=ibinx1,ib2=ibinx2)
                        zInt += zX
                        nAv += 1

                # vertical 
                if ibiny not in (1,nbiny): 
                    ibiny1 = FirstNonZeroBinI(h=h,ibinx=ibinx,ibiny=ibiny,step=-1,xory='y',nret=1,nstepmax=-1)
                    ibiny2 = FirstNonZeroBinI(h=h,ibinx=ibinx,ibiny=ibiny,step= 1,xory='y',nret=1,nstepmax=-1)
                    if ibiny1 > 0 and ibiny2 > 0: 
                        zY = PolateBin(h,ibinx,ibiny,xory='y',ib1=ibiny1,ib2=ibiny2)
                        zInt += zY
                        nAv += 1

                 
                ### Fill 
                if zInt: 
                    zInt /= (1.*nAv)
                    h2.SetBinContent(ibinx,ibiny, zInt)
                    if VB: print('Interpolating (ibinx,ibiny) = (%i,%i): %.3f' %(ibinx,ibiny,zInt))
                else: 
                    if VB>1: print('No interpolation for (%2i,%2i) = (%4.0f,%4.0f)' %(ibinx,ibiny, binx,biny))



    # Extrapolate
    if 'extra' in mode: 
        for ibinx in range(1,nbinx+1): 
            binx = tax.GetBinCenter(ibinx)
            # ---------------------------- pinpoint exactly which bins to polate
            if bins2polateX and not biniContainsbinxs(tax, ibinx, bins2polateX): 
                skippedX.append((ibinx,binx))
                continue

            for ibiny in range(1,nbiny+1): 
                biny = tay.GetBinCenter(ibiny)
                # ---------------------------- pinpoint exactly which bins to polate
                if bins2polateY and not biniContainsbinxs(tay, ibiny, bins2polateY): 
                    skippedY.append((ibiny,biny))
                    continue
                
                bintag = "(%2i,%2i) = (%4.0f,%4.0f)" %(ibinx,ibiny, binx,biny)
                
                z = h.GetBinContent(ibinx,ibiny)
                if z != 0: continue
                nAv = 0
                zExt = 0.

                ### Extrapolate:   ... should maybe move this double loop outside of the double bin loop
                for xory in xorys: #['x','y']: 
                    for step in steps: #[-1,1]: 
                        
                        if nAv > 0: break  # with this line the first is picked, so no averaging
                        ibins = FirstNonZeroBinI(h=h,ibinx=ibinx,ibiny=ibiny,step=step,xory=xory,nret=2,nstepmax=-1)
                        if ibins != -1: 
                            if step == -1: zzExt = PolateBin(h,ibinx,ibiny,xory=xory,ib1=ibins[1],ib2=ibins[0])
                            if step ==  1: zzExt = PolateBin(h,ibinx,ibiny,xory=xory,ib1=ibins[0],ib2=ibins[1])
                            if VB: print('Extrapolating %s in %s with step %i : %.3f' %(bintag, xory, step, zzExt))

                            # =====
                            # allowing a check on min/max (e.g. to allow dropping negative extrapolations)
                            if zpolmax and zzExt > zpolmax: 
                                zpolmaxval = GetBinContentXorY(h=h,ibinx=ibinx,ibiny=ibiny,ibinreplace=ibins[0], xory=xory)
                                if zpolmaxmode == 'absmax': 
                                    zpolmaxval = zpolmax
                                print('     NB: zzExt of %.3f  set to zpolmaxval = %.3f' %(zzExt,zpolmaxval))
                                zzExt = zpolmaxval
                            if zpolmin and zzExt < zpolmin: 
                                zpolminval = GetBinContentXorY(h=h,ibinx=ibinx,ibiny=ibiny,ibinreplace=ibins[0], xory=xory)
                                if zpolminmode == 'absmin': zpolminval = zpolmin
                                print('     NB: zzExt of %.3f  set to zpolminval = %.3f' %(zzExt,zpolminval))
                                zzExt = zpolminval
                            # =====

                            zExt += zzExt
                            nAv += 1
                        

                ### Fill (general code to allow for averaging over several extrapolations)
                if zExt: 
                    zExt /= (1.*nAv)

                    '''
                    # allowing a check on min/max (e.g. to allow dropping negative extrapolations)
                    if zpolmax and zExt > zpolmax: 
                        print 'INFO: zExt of %.3f  set to zpolmax = %.3f' %(zExt,zpolmax)
                        zExt = zpolmax
                    if zpolmin and zExt < zpolmin: 
                        print 'INFO: zExt of %.3f  set to zpolmin = %.3f' %(zExt,zpolmin)
                        zExt = zpolmin
                    '''

                    h2.SetBinContent(ibinx,ibiny, zExt)
                    #if VB: print 'Extrapolating %s: %.3f' %(bintag,zExt)
                else: 
                    if VB: print('No extrapolation for %s' %(bintag))

   
    # Some verbosity
    if skippedX and VB: print('INFO: LinearPolate: skippedX = %s' %(skippedX))
    if skippedY and VB: print('INFO: LinearPolate: skippedY = %s' %(skippedY))

    # Finalise
    if hmode in ['new','newdiff']: 
        return h2

###############################################

# =============================================================================
#def GetBinFromValue(ax,x): #use  TAxis::FindBin(value) instead (returns the appropriate bin number)
#    nbin = ax.GetNbins()
#    for ibin in range(1,nbin+1): 
#        ax.GetBinLowEdge(ibin) < x < ax.GetBinUpEdge(ibin): return ibin
#    return -1



# =============================================================================
def HistoEdit(h,hn='',opt={}): 
    if not hn: hn = h.GetName() + '_edit'
    h2 = h.Clone(hn)
    nx = h.GetNbinsX()
    ny = h.GetNbinsY()

    

    for ix in range(1,nx+1):
        for iy in range(1,ny+1): 

            hval = h.GetBinContent(ix,iy)
            # max
            if 'hmax' in opt and hval > opt['hmax']:
                hval = opt['hmax']
                h2.SetBinContent(ix,iy, hval)

            # rounding
            if 'Round' in opt: 
                for k in list(opt['Round'].keys()):
                    zmin=k[0]; zmax=k[1]
                    unit = opt['Round'][k]
                    # print 'zmin: %.0f   zmax: %.0f   unit: %.0f' %(zmin, zmax, unit)
                    if zmin < hval < zmax: 
                        rest = fmod(hval, unit)
                        hval2 = hval - rest + (rest > unit*0.1) * unit
                        h2.SetBinContent(ix,iy, hval2)


    return h2
# =============================================================================

class ColTest:
    def __init__(self,opt={}):
        # This method creats a marker, gives appropriate (default/input) attributes, plots and returns

        # Create and set marker
        self.mk = ROOT.TMarker()
        self.mk.SetNDC()

        # Default settings
        self.mk.SetX(0.02)
        self.mk.SetY(0.98)
        self.mk.SetMarkerColor(kRed)
        self.mk.SetMarkerSize(3.)

        # Additional settings
        self.SetAtt(opt=opt)


    def SetAtt(self,opt={}):
        # allow simplified for col change (only a number, not a dict)
        if type(opt) == int: opt = {'col': opt}
        if 'x' in opt: self.mk.SetX(opt['x'])
        if 'y' in opt: self.mk.SetY(opt['y'])
        if 'col' in opt: self.mk.SetMarkerColor(opt['col'])
        if 'msize' in opt: self.mk.SetMarkerSize(opt['msize'])


    def Draw(self,opt={}):
        self.SetAtt(opt)
        self.mk.Draw()

# =============================================================================
class NextColourStyle:
    def __init__(self, styles=[], Ncol=-1, cset=0, colFirst=1, argtext=''):
        if argtext:  # can use this to set/override Ncol, colFirst, styles
            z = argtext.split(':')
            Ncol = int(z.pop(0))
            styles = [int(zz) for zz in z.pop(0).split(',')]
            if len(z)>0: zcolFirst = int(z[0])
            else: zcolFirst = 1
        # -----

        self.nextColour = NextColour(cset=cset, Nmax=Ncol)
        self.Ncol = self.nextColour.Nmax  # is the max one admitted if not set (or set to -1)

        self.styles = styles
        self.colFirst = colFirst
        # colFirst=1 means go through all colours with a given style
        # then increment the style (and restart on the colours)
        # colFirst=0 first goes through the styles, then increment colour (and restart on the styles)
        self.iCol = 0
        self.iStyle = 0

        # A bit fragile on the starting procedure
        if colFirst: self.iCol = -1
        else: self.iStyle = -1


    def currentColour(self):
        return self.nextColour.col(self.iCol)

    def currentStyle(self):
        return self.styles[self.iStyle]

    def __next__(self):
        if self.colFirst: 
            self.iCol += 1
            if self.iCol == self.Ncol:
                self.iCol = 0
                self.iStyle += 1

                # check if have come all around, if so: restart
                if self.iStyle == len(self.styles): self.iStyle = 0  

        else:
            self.iStyle += 1
            if self.iStyle == len(self(styles)):
                self.iStyle = 0
                self.iCol += 1
                
                # check if have come all around, if so: restart
                if self.iCol == self.Ncol: self.iCol = 0 

        return [self.currentColour(), self.currentStyle()]



# =============================================================================
class NextColour:
    def __init__(self, cset=0, Nmax=-1):

        self.cols = {}
        self.cols[0] = [2,64,209,223,220, 622,227,416, 791,794, 15,18, 217,606]
        
        if cset < len(self.cols): self.cset = cset
        else:
            print('Warning::NextColour  no colour set %i. Using the default (0) out of %i sets' %(cset, len(self.cols)))
            self.cset = 0

        if Nmax > -1: self.Nmax = Nmax
        else: self.Nmax = len(self.cols[self.cset])


        self.iter = -1

        
    def reset(self,cset=-1): 
        self.iter = -1
        if cset >= 0: self.cset = cset

    def prev(self):
        if self.iter < 0: self.iter = 0   # in case there is no previous, and hence iterator is at -1, increment to 0
        return self.cols[self.cset][self.iter]

    def col(self,indx): 
        if indx < len(self.cols[self.cset]): return self.cols[self.cset][indx]
        else:
            print('Warning::NextColour  Requested colour index %i is outside array length of %i, using colour[0]' %(int(mode), len(self.cols[self.cset])))
            return self.cols[self.cset][0]

    def next(self,mode=''): 
        self.iter += 1            
        if self.iter >= self.Nmax:
            if self.iter > len(self.cols[self.cset]):
                print('Warning::NextColour  In set %i of %i colours we now go back to 0' %(self.cset, len(self.cols[self.cset])))
            self.iter = 0
        return self.cols[self.cset][self.iter]

# =============================================================================
stdcol = {

    '0L' : 922 # grey, dark
    ,'1L' : 920 # grey, light
    
    #,'2L' : 414 # green
    #,'3L' : 57  # ~2, blue (dark)
    ,'2L' : 207 #  red  ### 865  # blue
    ,'3L' : 66  #  blue ### 100  # red
    ,'4L' : 415  # green, dark
    #,'4L' : 614  # purple

    #,'ge2L': 412  # green, light
    #,'ge3L': 66   # blue, intermediate
    ,'ge2L': 623  # red, light
    ,'ge3L': 851   # blue, light  ## 851  # blue, light
    #,'ge4L': 93   # red, light (orange-like)
    #,'ge4L': 413  # green ## #414 #609  # purple, light (pinkish)
    ,'ge4L': 596  # green ## #414 #609  # purple, light (pinkish)

    ,'0T' : 922 # grey, dark
    ,'1T' : 920 # grey, light
    
    ,'2T' : 414 # green
    ,'3T' : 66  # blue
    ,'4T' : 100 # red
    ,'5T' : 614 # purple

    ,'ge2T': 412  # green, light
    ,'ge3T': 70   # blue, light
    ,'ge4T': 93   # red, light (orange-like)
    ,'ge5T': 609  # purple, light (pinkish)



    ,'h'  :876  # lilla
    ,'Z'  :413  # groen
    ,'Z*' :409  # groen, lysmatt
    
    ,'W'  :220  # gul (matt)
    ,'W*' :91  # orange (90->96)[yellow->orange->red]


    # ~2L
    ,'sl+l' :2   # raud  [Ni->sl+l] (sl=lL/lR)
    ,'lL+l' :2   # raud  [Ni->lL+l]
    ,'lL+l b' :623   # raud  [Ni->lL+l]
    ,'lR+l' :623 # raud lys  [Ni->lR+l]
    ,'lR+l b' :2 # raud lys  [Ni->lR+l]
    ,'ll'   :622 # raud lysare  [Ni->Nj+ll,  Ci->Cj+ll]
    ,'X+ll' :622 # raud lysare  [Ni->Nj+ll,  Ci->Cj+ll]
    ,'ll+X' :622 # raud lysare  [Ni->Nj+ll,  Ci->Cj+ll] (pink,rosa)


    # ~1L
    ,'sl+v' :634 # raud moerk [C->sl+v]
    ,'lL+v' :634 # raud moerk [C->lL+v]
    ,'lR+v' :634 # raud moerk [C->lR+v]
    ,'lv'   :634 # raud moerk (->brun) [633]  [Ni->Cj+lv,  Cj->Ni+lv]  off-shell, 1l
    ,'X+lv' :634 # raud moerk (->brun) [633]  [Ni->Cj+lv,  Cj->Ni+lv]  off-shell, 1l
    ,'lv+X' :634 # raud moerk (->brun) [633]  [Ni->Cj+lv,  Cj->Ni+lv]  off-shell, 1l

    ,'vL+L' :634 # raud moerk  [Ci->vL+L]  (L=l,T) # needed?
    ,'vl+l' :634 # raud moerk [Ci->vl+l]
    

    # ~2T 
    ,'Ti+T' :227  # cyan  [Ni->Ti+T]
    ,'Ti+T b' :225  # cyan  [Ni->Ti+T]
    ,'Ti+T c' :226  # cyan  [Ni->Ti+T]  # hmm
    ,'TT'   :422  # cyan light  [Ni->Nj+TT,  Ci->Cj+TT]
    ,'X+TT' :422  # cyan light  [Ni->Nj+TT,  Ci->Cj+TT]
    ,'TT+X' :422  # cyan light  [Ni->Nj+TT,  Ci->Cj+TT]


    # ~1T
    ,'Ti+v' :227  # cyan  [Ci->Ti+v]
    ,'Tv'   :226  # cyan semi-dark [Ni->Cj+Tv,  Cj->Ni+Tv]
    ,'X+Tv' :226  # cyan semi-dark [Ni->Cj+Tv,  Cj->Ni+Tv]
    ,'Tv+X' :226  # cyan semi-dark [Ni->Cj+Tv,  Cj->Ni+Tv]

    ,'vT+T' :226  # 624  [ci->vT+T]


    # ~0T - "useless" (so far) decays into sneutrino 
    ,'vL+v' :43 # kaki-brun  [Ni->vL+v]  # needed?
    ,'vl+v' :43 # kaki-brun  [Ni->vl+v]
    ,'vl+v b' :41 # kaki-brun  [Ni->vl+v]
    ,'vT+v' :42  # 624  [Ni->vT+v]
    ,'vT+v b' :47  # 624  [Ni->vT+v]
    ,'vv'   :46  # 624  [Ni->Nj+vv, Ci->Cj+vv]
    ,'X+vv' :46  # 624  [Ni->Nj+vv, Ci->Cj+vv]
    ,'vv+X' :46  # 624  [Ni->Nj+vv, Ci->Cj+vv]
    ,'vv+X__v2' :622 # raud lysare  [Ni->Nj+ll,  Ci->Cj+ll] (pink,rosa)


    # quarks
    ,'qq'  :922 # semidark grey (920,921,922)
    ,'C1qq':921 # hmm
    ,'C1qq__v2':920 # hmm

    ,'PIs' : 218 # groen, militaer (217-8)
    ,'PIs__v2' : 800    # Using ROOT.kOrange actually gives the -h error ... (why?) Could somehow be related to this being a global dict
    
    # Higgses
    ,'HA' : 901 # lilla, lys
    ,'H'  : 901 # lilla, lys
    ,'A'  : 907 # lilla->rosa

    # Other
    ,'y'  : 798 # kOrange-2 (tja)
    ,'y__v2' :93  # orange (90->96)[yellow->orange->red]


    # Suprocs
    ,'N1N1': 922  # grey,dark
    ,'N1C1': 920  # grey,light
    
    ,'N1N2': 92   # 2L ; orange (91-96: yellow->orange->red) to be similar to C1C1)
    ,'C1C1': 220  # 2L ; yellow (matt, same as W)
    
    ,'N2C1': 413  # green

    ,'N2N2': 803  # brownish;  hardly ever occurs (adhoc)
    ,'N2N3': 64   # not decided
    ,'N2N4': 3    # not set
    ,'N1N4': 2    # not set

    ,'NCheavy': 2    # red
    ,'NNheavy': 877  #223, purple
    ,'CCheavy': 609  # pink
    ,'heavy':   877  # purple
    
    ,'slsl': 64      # blue
    ,'TT'  : 227     # cyan,
    
    }


# =============================================================================
def FromDict(mode=['getTGraph2D'], Dict={}, xFromKeyIndex='', yFromKeyIndex='', xyFromKeyIndices=[], FixKeyIndices={}, zKey='', xtit='x', ytit='y', tit='A quantity', zsc=1.0, xoff=0.8, yoff=1.0, xcent=1, ycent=1, VB=0, FORM='%.0f', txtsize=0.02, extrapolateflat={}):

    ############################ STATUS, usage
    # 
    # o the 'drawtext' mode is a hack, need to return a tex; works fine, but is inelegant
    # o txtsize: should have some automatic way of determining an appropriate size (e.g. when extrapolated, the text typically gets a bit too large) 
    #
    # xoff, yoff: default is 1.0   For 2D plotting that is nice, for 3D plotting (surf,lego) 1.7 is better
    #
    # ###########################
    # Ex:
    # cd ~/grids_lsp/res_DGemt_TB6_to700/dir_dicts
    # from ROOT import *
    # from kilelib_ROOT import FromDict
    # execfile('geneff_total.py')
    # gr2D = FromDict(mode=['getTGraph2D'], Dict=subproceff, xyFromKeyIndices=[2,1], FixKeyIndices={0:100}, zKey='eff', xtit='MU [GeV]', ytit='M2 [GeV]', tit='Generator Filter Efficiency [%]', zsc=100.)
    # gr2D.Draw('zcol')
    # tex = FromDict(mode=['drawtext'], Dict=subproceff, xyFromKeyIndices=[2,1], FixKeyIndices={0:100}, zKey='eff', xtit='MU [GeV]', ytit='M2 [GeV]', tit='Generator Filter Efficiency [%]', zsc=100.,txtsize=0.03)
    # c1.Update()
    # # gr2D.Draw('lego2')
    # # gr2D.Draw('surf1')
    #############################

    
    if type(mode) != list:
        print('kilelib_ROOT::FromDict  mode must be a list')
        return 'ERROR'

    modes = ['getTGraph2D', 'drawtext']

    # 1 --- PREPARATIONS 

    x = array('d',[])
    y = array('d',[])
    z = array('d',[])

    if xyFromKeyIndices:
        xFromKeyIndex = xyFromKeyIndices[0]
        yFromKeyIndex = xyFromKeyIndices[1]

    if not (xFromKeyIndex and yFromKeyIndex):
        print('WARNING  kilelib_ROOT::GetTGraph2D_fromDict   indices not given, x: %s    y: %s' %(xFromKeyIndex, yFromKeyIndex))
        return 0


    if 'drawtext' in mode: 
        tex = ROOT.TLatex()
        tex.SetTextAlign(22)
        tex.SetTextSize(txtsize)
    


    # 2 --- LOOP 
    dkeys = list(Dict.keys())
    dkeys.sort()
    for dkey in dkeys:
        dkeyT = str(dkey).replace(' ','')

        # Check if key to be skipped or not (if is not in FixKeyIndices)
        usethis = 1
        for fkey in list(FixKeyIndices.keys()):
            if dkey[fkey] != FixKeyIndices[fkey]: usethis = 0
        if not usethis: continue
        
        xval = ''
        yval = ''
        zval = ''
        
        if xFromKeyIndex: xval = float(dkey[xFromKeyIndex])   # 'if' not needed, is checked beforehand, 
        if yFromKeyIndex: yval = float(dkey[yFromKeyIndex])   #  ... but keep for the possibility to add more complex methods, 
        #                                                 #  ... e.g. x and/or y inside the dict, not in the key

        if zKey: zval = Dict[dkey][zKey]
        else: zval = Dict[dkey]   # if no z-key is given, assume it is a simple dict (i.e. one whose values are numbers, not dicts)

        zval *= zsc  # allow to scale up/down


        # is here possible to just plot the number found, e.g. to put text onto another
        if 'drawtext' in mode: 
            ztxt = FORM %(zval)
            #tex.SetX(xval)
            #tex.SetY(yval)
            #tex.SetTitle()
            #tex.DrawClone()
            tex.DrawText(xval,yval,ztxt)
        

        x.append(xval)
        y.append(yval)
        z.append(zval)

        #if VB: print '%4i  %s  (MU,M2) = (%3i,%3i) : eff = %7.4f' %(len(z), dkeyT, xval,yval,zval)
        if VB: print('%4i  %s  (x,y) = (%3i,%3i) : eff = %7.4f' %(len(z), dkeyT, xval,yval,zval))



    # 3 --- POST TREATMENT

    if 'getTGraph2D' in mode:

        # -------
        if extrapolateflat:
            ExtrapolateFlat(x,y,z, extrapolateflat)  # this changes x,y,z (adds points according to extrapolateflat

                

        # -------
        n = len(z)


        gr2D = ROOT.TGraph2D(n,x,y,z)
        
        gr2D.SetTitle(tit)
        
        tax = gr2D.GetXaxis()
        tay = gr2D.GetYaxis()
        
        tax.SetTitle(xtit)
        tay.SetTitle(ytit)
        
        tax.SetTitleOffset(xoff)
        tay.SetTitleOffset(yoff)
        
        tax.CenterTitle(xcent)
        tay.CenterTitle(ycent)



    # The actual return procedure should be improved

    if 'getTGraph2D' in mode: return gr2D

    if mode == 'plot': return tex


# =============================================================================
def ExtrapolateFlat(x,y,z, extrapolateflat):
    pol = extrapolateflat

    # First make a dummy TGraph2D, will help
    gr = ROOT.TGraph2D(len(z),x,y,z)
    tax = gr.GetXaxis()
    tay = gr.GetYaxis()

    delta = 1e-5  # apparently the very last value is outside, so need to add a delta
    
    # --- FIRST extend the blocks in the requrested directions

    if 'xmax' in pol:
        for ibiny in range(0,tay.GetNbins()+1):  # NB takes both edges
            yval = tay.GetBinUpEdge(ibiny)
            x.append(pol['xmax']+delta)
            y.append(yval)
            z.append(gr.Interpolate(gr.GetXmax(), yval))
            
    if 'xmin' in pol:
        for ibiny in range(0,tay.GetNbins()+1):  # NB takes both edges
            yval = tay.GetBinUpEdge(ibiny)
            x.append(pol['xmin']-delta)
            y.append(yval)
            z.append(gr.Interpolate(gr.GetXmin(), yval))

    if 'ymax' in pol:
        for ibinx in range(0,tax.GetNbins()+1):  # NB takes both edges
            xval = tax.GetBinUpEdge(ibinx)
            x.append(xval)
            y.append(pol['ymax']+delta)
            z.append(gr.Interpolate(xval, gr.GetYmax()))

    if 'ymin' in pol:
        for ibinx in range(0,tax.GetNbins()+1):  # NB takes both edges
            xval = tax.GetBinUpEdge(ibinx)
            x.append(xval)
            y.append(pol['ymin']-delta)
            z.append(gr.Interpolate(xval, gr.GetYmin()))


    # --- THEN, if needed add the corners (e.g. if xmax and ymax, then the new (xmax,ymax) corner is currently empty

    if 'xmax' in pol and 'ymax' in pol:  # top-right
        x.append(pol['xmax'])
        y.append(pol['ymax'])
        z.append(gr.Interpolate(gr.GetXmax(), gr.GetYmax()))

    if 'xmax' in pol and 'ymin' in pol:  # bottom-right
        x.append(pol['xmax'])
        y.append(pol['ymin'])
        z.append(gr.Interpolate(gr.GetXmax(), gr.GetYmin()))

    if 'xmin' in pol and 'ymax' in pol:  # top-right
        x.append(pol['xmin'])
        y.append(pol['ymax'])
        z.append(gr.Interpolate(gr.GetXmin(), gr.GetYmax()))

    if 'xmin' in pol and 'ymin' in pol:  # bottom-right
        x.append(pol['xmin'])
        y.append(pol['ymin'])
        z.append(gr.Interpolate(gr.GetXmin(), gr.GetYmin()))


    del tax,tay,gr

    # don't have to return anything, the changes to x,y,z are already done 

# =============================================================================

                 

LEPlimit_C1_gaugino_list = [
    # [<mC1-mN1>, <C1 limit>]
    [0.001,103.0],  # asymptotic value (hack)
    [0.10, 103.0],
    [0.14, 102.8],
    [0.145, 93.0],
    [0.15,  92.0],
    [0.20,  93.4],
    [0.21,  93.5],
    [0.22,  95.7],
    [0.30,  96.2],
    [0.35,  95.4],
    [0.40,  95.7],
    [0.50,  96.0],
    [0.60,  96.3],
    [0.70,  96.5],
    [0.80,  96.7],
    [0.90,  96.9],
    [1.00,  97.0],
    [1.50,  97.0],
    [2.00,  96.1],
    [3.00, 101.0],
    [3.75, 102.1],
    [4.00, 102.6],
    [5.00, 102.8],
    [6.00, 102.4],
    [7.00, 102.7],
    [8.00, 102.8],
    [9.00, 102.8],
    [10.0, 102.9],
    [110 , 102.9]  # asymptotic value (hack)
    ]


def LEPlimit_C1_gaugino_func():
    # Array above was taken from LEP plot
    # delM = delM0 x exp(y/y0) with delM0 = 0.10 GeV and y0 = 0.20 cm and y is the distance down to y = 10^-1 ('zero')
    # The formula above gives then the correct delM value in between
    

    # First load list
    x = array('d')
    y = array('d')
    zlist = list(LEPlimit_C1_gaugino_list)
    zlist.sort()
    for pair in zlist:
        x.append(pair[0])  # mass difference  C1-N1
        y.append(pair[1])  # LEP limit on C1


    # Then create TGraph
    C1_limit = ROOT.TGraph(len(x),x,y)

    # Can then get C1-limit as a function of delM=m(C1)-m(N1) by C1_limit.Eval(delM)
    return C1_limit


# =============================================================================

