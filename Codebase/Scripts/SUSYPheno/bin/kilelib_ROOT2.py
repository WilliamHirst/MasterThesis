'''
The main ROOT library is kilelib_ROOT.py
A potential problem with that file is that it does from ROOT import * at start.
Hence, you need to have ROOT sourced first ...
Sometimes it is relevant to use such methods without having ROOT "explicitly" (globally?) sourced
Those cases should be put in this file, for now at least

'''
from math import fmod, ceil


# #####
def BLegends(toTLegend, objarr=[], titarr=[], optD={}):

    # The basic box
    x1 = optD.get('x1', 0.12)
    if 'xwid' in optD and not 'x2' in optD:
        xwid = optD.get('xwid')
        x2 = x1+xwid
    elif 'x2' in optD and not 'xwid' in optD:
        x2 = optD['x2']
        xwid = x2-x1
    elif 'x2' in optD and 'xwid' in optD: 
        x2 = optD['x2']
        xwid = optD['xwid']
        x1 = xwid - x2
    else:
        xwid = 0.13   # default (could make dynamic somehow)
        x2 = x1+xwid

        
    y2 = optD.get('y2',0.895)
    DY = optD.get('DY',0.035)

    # Attributes
    sty = optD.get('sty','')   # 'L' seems like a sensible default, but what is? 
    fillcol = optD.get('fillcol',0)
    bordersize = optD.get('bordersize',1)
    maxchar = optD.get('maxchar',9999)
    maxcharoverflow = optD.get('maxcharoverflow','')


    nGraphs = len(objarr)
    if len(objarr) != len(titarr):
        print('Warning: Legends len(objarr) = %i  !=  %i = len(titarr)' %(len(objarr), len(titarr)))


    # The number of boxes : can specify both (cartesian forced) or one (flexible)
    nbox = optD.get('nbox',-1)
    nperbox = optD.get('nperbox',-1)
    if nperbox == -1 and nbox == -1:
        nbox = 1
        nperbox = nGraphs
    elif nperbox == -1:
        nperbox = ceil(1.0 * nGraphs / nbox)  # ceiling  (math.ceil)   [also exists: floor, math.floor]
    elif nbox == -1:
        nbox = ceil(1.0 * nGraphs / nperbox)


    # for rescaling hmax if desired
    margtop = optD.get('margtop',0.10)
    margbot = optD.get('margbot',0.10)
    

    y1min = 1.

    legs = []  # array of legends
    # Loop over array and produce Legends (to be returned)  [Fill first box first, then next etc.] [can implement horisontal too]
    hmax = -99999
    for iG in range(nGraphs):

        # Need to make new box? 
        if fmod(iG, nperbox) == 0:
            thisx1 = x1 + len(legs)*xwid
            thisx2 = x2 + len(legs)*xwid
            thisy2 = y2
            ninthisbox = min( nGraphs-len(legs)*nperbox, nperbox )   # the number of graphs in this box
            thisy1 = y2 - DY*ninthisbox
            y1min = min(y1min, thisy1)

            # print thisx1,thisy1,thisx2,thisy2, xwid

            leg = toTLegend(thisx1,thisy1,thisx2,thisy2)
            leg.SetFillColor(fillcol)
            leg.SetBorderSize(bordersize)
            legs.append(leg)
        # -----
        hmax = max(hmax, objarr[iG].GetMaximum())

        #thisname = titarr[iG][0:min(len(titarr[iG]), maxchar)]
        if len(titarr[iG]) > maxchar:
            thisname = titarr[iG][0:maxchar-len(maxcharoverflow)] + maxcharoverflow
        else:
            thisname = titarr[iG]

        #print 'DEBUG: ', titarr[iG], len(titarr[iG]), maxchar, thisname, min(len(titarr[iG]), maxchar), min(2,3), type(maxchar)
        
        if 'sty' in optD: leg.AddEntry(objarr[iG], thisname, sty)
        else: leg.AddEntry(objarr[iG], thisname)  # do like this because don't know what the default sty-value is (it is not '')  # ... I actually think it is 'lpf' ... if you have the 'f' (fill?), you will get a border


    oldactiveY = 1. - margtop - margbot
    newactiveY = y1min - margbot
    hmaxnew = hmax * oldactiveY / newactiveY


    # Default is legend(s) starting from upper left
    # Here is a hack to optionally move them to the right and/or to the bottom
    for index,leg in enumerate(reversed(legs)): 
        if 'right' in optD:  # ex: optD['right'] = 0.88
            if index == 0: delX = optD['right'] - leg.GetX2()  # amount to be moved to the right
            leg.SetX2(leg.GetX2()+delX)
            leg.SetX1(leg.GetX1()+delX)
        if 'bottom' in optD: # ex: optD['bottom'] = 0.12
            if index == 0: delY = leg.GetY1() - optD['bottom']
            leg.SetY1(leg.GetY1()-delY)
            leg.SetY2(leg.GetY2()-delY)


    res = {'legends':legs, 'y1min':y1min, 'hmaxrescale':oldactiveY/newactiveY, 'hmaxscaleincrease':oldactiveY/newactiveY-1., 'hmax':hmax, 'hmaxnew':hmaxnew}  # knowing y1min allows the user to increase ymax of the plot to avoid the legend cover the curves

    return res

# #####

