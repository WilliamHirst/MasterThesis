#from ROOT import *   # 2014-01-24 replaced by the below (to protect isa_simplereader_class and munch
import ROOT  
from array import *
import random
# gROOT.LoadMacro("/mn/kvant/u1/borgeg/scripts/sky_AtlasStyle.C")
# SetAtlasStyle()
from kilelib_ROOT import *
from kilelib import safeGet
import kilelib_physics as mypdg
from kilelib_physics import Z_LLR

import os

from libDarksusy import Darksusy_these_bits_notexcluded
from kilelib import IsInRangesT

from kilelib_ROOT import LEPlimit_C1_gaugino_func


# =================================
def DefineContours(txtline):
    # Example of txtline (one only)
    #   mC1/105,[103.5,110],2,7,2    # id/code,[contourValues],colour,style,width
    #   105,[103.5,110],2,7,2        # id dropped (becomes str(code))
    #   std:[mC1,mN2]                # std: looks through the 'list': these are then predefined/standard 

    # Any kind of the above contour-descriptions can then be combined with a ';'
    
    contours = {}
    z = txtline.split(';')

    # Loop over all the contour descriptions (separated by ';')
    for iz in range(len(z)):
        
        # First see if there are any standard (predefined) contours
        if z[iz].startswith('std:'):
            zz = z[iz][4:].split(',')

            for zzz in zz:

                # C1:RED
                if zzz == 'mC1':  
                    contours['mC1'] = {'yis':105, 'val':[103.5], 'col':2, 'sty':7, 'wid':2}
                elif zzz == 'mC1a':  # identical to mC1
                    contours['mC1a'] = {'yis':105, 'val':[103.5], 'col':2, 'sty':7, 'wid':2}  
                elif zzz == 'mC1b':
                    contours['mC1b'] = {'yis':105, 'val':[96], 'col':623, 'sty':7, 'wid':2}
                elif zzz == 'mC1c':
                    contours['mC1c'] = {'yis':105, 'val':[92], 'col':820, 'sty':7, 'wid':2}
                elif zzz == 'mC1ab':
                    contours['mC1ab'] = {'yis':105, 'val':[96,103.5], 'col':2, 'sty':7, 'wid':2}


                # C1-N1:BLUE
                elif zzz == 'C1mN1a': 
                    contours['C1mN1a'] = {'yis':401, 'val':[3.0], 'col':4, 'sty':7, 'wid':2}
                elif zzz == 'C1mN1b': 
                    contours['C1mN1b'] = {'yis':401, 'val':[0.15], 'col':7, 'sty':7, 'wid':2}
                elif zzz == 'C1mN1ab': 
                    contours['C1mN1ab'] = {'yis':401, 'val':[0.15,3.0], 'col':4, 'sty':7, 'wid':2}

                elif zzz == 'meR': 
                    contours['meR'] = {'yis':121, 'val':[100.], 'col':1, 'sty':2, 'wid':2}
                elif zzz == 'meR95100': 
                    contours['meR95100'] = {'yis':121, 'val':[95.,100.], 'col':1, 'sty':2, 'wid':2}
                elif zzz == 'meL': 
                    contours['meL'] = {'yis':122, 'val':[100.], 'col':1, 'sty':2, 'wid':2}
                elif zzz == 'mT1': 
                    contours['mT1'] = {'yis':124, 'val':[90.,95.], 'col':1, 'sty':2, 'wid':2}
                
                elif zzz == 'mC1range100': 
                    contours[zzz] = {'yis':105, 'val':list(range(100,1001,100)), 'col':1, 'sty':2, 'wid':2}
                elif zzz == 'mN2range100': 
                    contours[zzz] = {'yis':102, 'val':list(range(100,1001,100)), 'col':1, 'sty':2, 'wid':2}
                elif zzz == 'mN1range50': 
                    contours[zzz] = {'yis':101, 'val':list(range(50,1001,50)), 'col':1, 'sty':2, 'wid':2}

                # should here add some on darkSUSY, mass differences, 
                    
                else:
                    print('WARNING::DefineContours  Skipping unknown predefined id: %s' %(zzz))

                    
            continue

        # Then do the improvised/on-the-fly contours
        zz = z[iz].split(',')
        
        zIDandCode = zz.pop(0).split('/')
        zid = zIDandCode.pop(0)
        if zIDandCode: zcode = int(zIDandCode.pop(0))
        else: zcode = int(zid)
        #if zid not in contourIDs.keys(): sys.exit('FATAL::isa_simplereader  unknown contour code: %s' %(zid))
        
        contours[zid] = {'yis':zcode, 'val':[], 'col':2, 'sty':7, 'wid':2}
        if zz:
            zvals = zz.pop(0).replace('[','').replace(']','').split(',')
            for zval in zvals: contours[zid]['val'].append(float(zval))
        if zz: contours[zid]['col'] = int(zz.pop(0))
        if zz: contours[zid]['sty'] = int(zz.pop(0))
        if zz: contours[zid]['wid'] = int(zz.pop(0))
        if zz: sys.exit('FATAL: incorrect contour format: %s' %(txtline))
    return contours


# =================================
def GetContourHistos(gg, contours, Scan, flip=0):
    hContours = {}

    # Loop over the contours
    for conid in list(contours.keys()): 
        con = contours[conid]
        # First get the histogram of the quantity in question
        htmp = get2DhistM1M2MU(gg=gg, yis=con['yis'], Scan=Scan, ret='h', flip=flip)
        # Then make contour and set attributes
        conts = array('d',con['val'])
        htmp.SetContour(len(conts),conts)
        htmp.SetLineColor(con['col'])
        htmp.SetLineStyle(con['sty'])
        htmp.SetLineWidth(con['wid'])
        # Finally put in dictionary
        hContours[conid] = htmp
        del htmp

    # return the dictionary with contours, these can then elsewhere be plotted with e.g. Draw('cont3same')
    return hContours




# =================================
def GetAxisRangesAuto(val, loweredge='', upperedge='', firstval='', lastval=''):
    binsL = []

    zval = list(val)
    if type(firstval) in [int,float]: zval[0] = firstval
    if type(lastval)  in [int,float]: zval[-1] = lastval

    # print '\n zval: ', zval

    # left edge
    if type(loweredge) in [int,float]: binsL.append(loweredge)
    else: binsL.append(zval[0] - (zval[1]-zval[0])/2.)

    # all intermediate
    for i in range(len(zval)-1): binsL.append((zval[i] + zval[i+1])/2.)

    
    # right edge
    if type(upperedge) in [int,float]: binsL.append(upperedge)
    else: binsL.append(zval[-1] + (zval[-1] - zval[-2])/2.)

    return binsL

    
# =================================
def GetVarText(varin, outtype):
    outtypes = ['var','varTxt','varLatex']
    if outtype not in outtypes: 
        print('Error::GetVarText  Illegal type %s not in %s' %(outtype, outtypes))
        return 'UnknownVarFor'+varin
    vars = {}
    vars['var'] = {'M_1':'M_1','M1':'M_1','M_{1}':'M_1', 'M_2':'M_2','M2':'M_2','M_{2}':'M_2', 'MU':'MU','#mu':'MU'}
    vars['varTxt'] = {'M_1':'M1','M1':'M1','M_{1}':'M1', 'M_2':'M2','M2':'M2','M_{2}':'M2', 'MU':'MU','#mu':'MU'}
    vars['varLatex'] = {'M_1':'M_{1}','M1':'M_{1}','M_{1}':'M_{1}', 'M_2':'M_{2}','M2':'M_{2}','M_{2}':'M2', 'MU':'#mu','#mu':'#mu'}

    if varin not in vars[outtype]: 
        print('Error::GetVarText  Var %s not in known vars ' %(varin, list(vars[outtype].keys())))
        return 'UnknownVarFor'+varin

    return vars[outtype][varin]
# =================================
def VarTranslateOld(v1, v2):
    vars = {}
    vars['v1'] = v1
    vars['v2'] = v2
    vars['V1'] = GetVarText(v1,'varTxt')
    vars['V2'] = GetVarText(v2,'varTxt')
    vars['v1T'] = GetVarText(v1,'varLatex') + ' [GeV]'
    vars['v2T'] = GetVarText(v2,'varLatex') + ' [GeV]'

    return vars


# =================================
def GetAxisRangesEtcNew(gg, figname_class='', id2=''):

    sliceID = gg['scanID']
    sliceInfo = gg['asliceInfo']
    
    setup = {}
    val = {}
    binsL = {}

    lastval = gg.get('lastval','')
    lastvalv1 = gg.get('lastvalv1',lastval)
    lastvalv2 = gg.get('lastvalv2',lastval)

    # intermediate (for simplicity of reference)
    v1 = sliceInfo['v1']
    v2 = sliceInfo['v2']
    fixvar = sliceInfo['fixvar']
    fixval = sliceInfo['fixval']
    

    zvars = VarTranslateOld(v1, v2) 
    for key in list(zvars.keys()): setup[key] = zvars[key]
    #print zvars
    V1 = zvars['V1']
    V2 = zvars['V2']

    val['v1'] = sliceInfo['arr']['v1']
    val['v2'] = sliceInfo['arr']['v2']
    if lastvalv1 != '': val['v1'][-1] = lastvalv1  # Direct and actual change 2013-07-16
    if lastvalv2 != '': val['v2'][-1] = lastvalv2
    binsL['v1'] = GetAxisRangesAuto(val['v1'],lastval=lastvalv1)  # With the direct dchange in val above, the use of lastval is not needed here
    binsL['v2'] = GetAxisRangesAuto(val['v2'],lastval=lastvalv2)
    #print '\nbinsL: ', binsL['v1']
    if 1:
        val[V1] = val['v1']
        val[V2] = val['v2']
        binsL[V1] = binsL['v1']
        binsL[V2] = binsL['v2']

    #print 'KEYS: ', binsL.keys()


    setup['val'] = val
    setup['binsL'] = binsL
    setup['v1'] = v1
    setup['v2'] = v2

    fnn0 = "%s_NOTUSED_%seq%s" %(sliceID, fixvar, str(int(fixval)).zfill(3))
    if figname_class: fnn0 += '_' + figname_class
    setup['fnn0'] = fnn0

    return setup


# =================================
def GetAxisRangesEtc(sliceID, fixvar, fixval, figname_class='', id2=''):

    setup = {}
    val   = {}
    binsL = {}

    if 1: # 2012-06-01: inserted as default

        if fixvar == 'M_1':  # this is actually the only valid one tof fox for plane-plotting
            v1 = 'MU';  V1 = 'mu'; v1T = 'MU [GeV]'
            v2 = 'M_2'; V2 = 'M2'; v2T = 'M2 [GeV]'

            fnn0 = "%s_MUM2_M1eq%s" %(sliceID, str(int(fixval)).zfill(3))
        else: 
            print('FATAL  libIsaplot::get2DhistM1M2MU  ... %s not made for fixing %s' %(sliceID, fixvar))
            return 



    # -----------------------------------------------------------------
    # Just set these as defaults   2012-06-29
    v1 = 'MU';  V1 = 'mu'; v1T = 'MU [GeV]'
    v2 = 'M_2'; V2 = 'M2'; v2T = 'M2 [GeV]'
    fnn0 = "%s_MUM2_M1eq%s" %(sliceID, str(int(fixval)).zfill(3))
    # -----------------------------------------------------------------


    if sliceID in ['xsecXX4all']: 
        val['M1'] = [80,90,100,110,120,140,160,180,200, 220, 240,260,280,300,999]
        val['M2'] = [80,90,100,110,120,140,160,180,200, 220, 240,260,280,300,999]
        val['mu'] = [80,90,100,110,120,140,160,180,200, 220, 240,260,280,300,999]
        binsL['M1'] = [75,85,95,105,115,130,150,170,190,210,230,250,270,290,310,350]
        binsL['M2'] = list(binsL['M1'])
        binsL['mu'] = list(binsL['M1'])


    if sliceID in ['xsecXX4all']:   # might be more general than the above one
        if fixvar == 'MU': 
            v1 = 'M_1'; V1 = 'M1'; v1T = 'M1 [GeV]'
            v2 = 'M_2'; V2 = 'M2'; v2T = 'M2 [GeV]'
            fnn0 = "XX4all_M1M2_MUeq%.0f" %fixval
        elif fixvar == 'M_1':
            v1 = 'MU';  V1 = 'mu'; v1T = 'MU [GeV]'
            v2 = 'M_2'; V2 = 'M2'; v2T = 'M2 [GeV]'
            fnn0 = "XX4all_MUM2_M1eq%.0f" %fixval
        elif fixvar == 'M_2':
            v1 = 'M_1'; V1 = 'M1'; v1T = 'M1 [GeV]'
            v2 = 'MU';  V2 = 'mu'; v2T = 'MU [GeV]'
            fnn0 = "XX4all_M1MU_M2eq%.0f" %fixval
        else:
            print("FATAL  libisaplot::get2DhistM1M2MU  incorrect fixvar: %s" %(fixvar))



    if sliceID in ['DGemt','DGnoL','DG_n147']: #Mar15
        val['M1'] = [100,140,250]
        val['M2'] = [100,110,120,140,160,180,250]
        val['mu'] = list(val['M2'])
        binsL['M1'] = [80, 120, 200,300]  # not in use (M1 is always fixed)
        binsL['M2'] = [95,105,115,130,150,170,200,260]
        binsL['mu'] = list(binsL['M2'])



    if sliceID in ['DGemt','DGnoL','DGemt_n300','DGemt_n768','DGemtFine','DGemtLR_n300','DG_n300','DG_n147','DGemt_to700','DGemtFine_to700','DGemtFine_to500']: 

        if fixvar == 'M_1':  # this is actually the only valid one tof fox for plane-plotting
            v1 = 'MU';  V1 = 'mu'; v1T = 'MU [GeV]'
            v2 = 'M_2'; V2 = 'M2'; v2T = 'M2 [GeV]'

            fnn0 = "%s_MUM2_M1eq%s" %(sliceID, str(int(fixval)).zfill(3))
            if 0: 
                fnn0 = "BKGjelsten_%s_MUM2_M1eq%s" %(sliceID.split('_')[0], str(int(fixval)).zfill(3))  # HACK FOR 3L 2012-02-04
                print("INFO  libIsaplot::get2DhistM1M2MU  Filename hack for 3L note: fnn0 = %s" %(fnn0))

        else: 
            print('FATAL  libIsaplot::get2DhistM1M2MU  ... %s not made for fixing %s' %(sliceID, fixvar))
            return 


    if sliceID in ['M2xMU=100-700','M2xMU=100-700_18x18']:
        val['M2'] =    [100,110,120,140,160, 180.0, 210.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0] #, 1500.0]   # overflow for 1500
        binsL['M2'] = [95,105,115,130,150, 170, 195.0, 230.0, 275.0, 325.0, 375.0, 425.0, 475.0, 525.0, 575.0, 625.0, 675.0, 725.0] #, 775.0]
        val['mu'] = list(val['M2'])
        binsL['mu'] = list(binsL['M2'])


    if sliceID in ['M2xMU=100-500','M2xMU=100-500_13x13']:
        val['M2'] =    [100,110,120,140,160, 180.0, 210.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0] # overflow for 1500
        binsL['M2'] = [95,105,115,130,150, 170, 195.0, 230.0, 275.0, 325.0, 375.0, 425.0, 475.0, 525.0]
        val['mu'] = list(val['M2'])
        binsL['mu'] = list(binsL['M2'])


    if sliceID in ['M2xMU=100-350','M2xMU=100-350_10x10']:
        val['M2'] =    [100,110,120,140,160, 180.0, 210.0, 250.0, 300.0, 350.0]
        binsL['M2'] = [95,105,115,130,150, 170, 195.0, 230.0, 275.0, 325.0, 375.0]
        val['mu'] = list(val['M2'])
        binsL['mu'] = list(binsL['M2'])


    if sliceID in ['M2xMU=100-250','M2xMU=100-250_8x8']:
        val['M2'] =    [100,110,120,140,160, 180.0, 210.0, 250.0]
        binsL['M2'] = [95,105,115,130,150, 170, 195.0, 230.0, 275.0]
        val['mu'] = list(val['M2'])
        binsL['mu'] = list(binsL['M2'])


    if sliceID in ['DGemt_to700']:  #May12
        val['M1'] = [100,140,250]
        val['M2'] =    [100,110,120,140,160, 180.0, 210.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0]
        binsL['M2'] = [95,105,115,130,150, 170, 195.0, 230.0, 275.0, 325.0, 375.0, 425.0, 475.0, 525.0, 575.0, 625.0, 675.0, 725.0]
        val['mu'] = list(val['M2'])
        binsL['M1'] = [80, 120, 200,300]  # not in use (M1 is always fixed)
        binsL['mu'] = list(binsL['M2'])



    # ---------- Begin: without 110 
    if sliceID in ['M2xMU=100-700b']: 
        val['M2'] =    [100,120,140,160, 180.0, 210.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0] #, 1500.0]   # overflow for 1500
        binsL['M2'] = [90,110,130,150, 170, 195.0, 230.0, 275.0, 325.0, 375.0, 425.0, 475.0, 525.0, 575.0, 625.0, 675.0, 725.0] #, 775.0]
        val['mu'] = list(val['M2'])
        binsL['mu'] = list(binsL['M2'])


    if sliceID in ['M2xMU=100-500b']: 
        val['M2'] =    [100,120,140,160, 180.0, 210.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0] # overflow for 1500
        binsL['M2'] = [90,110,130,150, 170, 195.0, 230.0, 275.0, 325.0, 375.0, 425.0, 475.0, 525.0]
        val['mu'] = list(val['M2'])
        binsL['mu'] = list(binsL['M2'])


    if sliceID in ['M2xMU=100-350b']: 
        val['M2'] =    [100,120,140,160, 180.0, 210.0, 250.0, 300.0, 350.0]
        binsL['M2'] = [90,110,130,150, 170, 195.0, 230.0, 275.0, 325.0, 375.0]
        val['mu'] = list(val['M2'])
        binsL['mu'] = list(binsL['M2'])


    if sliceID in ['M2xMU=100-250b']: 
        val['M2'] =    [100,120,140,160, 180.0, 210.0, 250.0]
        binsL['M2'] = [90,110,130,150, 170, 195.0, 230.0, 275.0]
        val['mu'] = list(val['M2'])
        binsL['mu'] = list(binsL['M2'])


    if sliceID in ['DGemt_to700b']:  
        val['M1'] = [100,140,250]
        val['M2'] =    [100,120,140,160, 180.0, 210.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0]
        binsL['M2'] = [90,110,130,150, 170, 195.0, 230.0, 275.0, 325.0, 375.0, 425.0, 475.0, 525.0, 575.0, 625.0, 675.0, 725.0]
        val['mu'] = list(val['M2'])
        binsL['M1'] = [80, 120, 200,300]  
        binsL['mu'] = list(binsL['M2'])

    # ---------- End: without 110 
        

        
    if sliceID in ['DGemt_n300','DGemtLR_n300','DG_n300']: # Nov24
        val['M1'] = [100,140,250]
        val['M2'] =    [100,110,120,140,160, 180.0, 210.0, 250.0, 300.0, 350.0]
        binsL['M2'] = [95,105,115,130,150, 170, 195.0, 230.0, 275.0, 325.0, 375.0]
        val['mu'] = list(val['M2'])
        binsL['M1'] = [80, 120, 200,300]  # not in use (M1 is always fixed)
        binsL['mu'] = list(binsL['M2'])


    if sliceID in ['DGemt_n768']: # Nov24
        val['M1'] = [100,140,250]
        val['M2'] =    [100,110,120,140,160, 180.0, 195.0, 210.0, 230.0, 250.0, 275.0, 300.0, 325.0, 350.0, 375.0, 400.0]
        binsL['M2'] = [95,105,115,130,150,170,  187.5, 202.5, 220.0, 240.0, 262.5, 287.5, 312.5, 337.5, 362.5, 387.5, 412.5]
        val['mu'] = list(val['M2'])
        binsL['M1'] = [80, 120, 200,300]  # not in use (M1 is always fixed)
        binsL['mu'] = list(binsL['M2'])


    if sliceID in ['DGemtFine']: # Nov24  # NB: this is meant to do together with allIDs/allIDs_DGemtFine_M1_*bn1369  NOT _n1521 
        val['M1'] = [100,140,250]
        val['M2'] =    [100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200, 210,220,230,240,250,260,270,280,290,300,310,320,330,340,350]
        val['mu'] = list(val['M2'])
        # quicky for binsL (but is anyway replaced if 'binwidth' is set
        binsL['M2'] = [95]
        for i in range(len(val['M2'])-1): binsL['M2'].append( (val['M2'][i]+val['M2'][i+1])/2 )
        binsL['M2'].append(355)
        binsL['mu'] = list(binsL['M2'])
        #print 'val   : ' , val
        #print 'binsL : ', binsL

    if sliceID in ['DGemtFine_to700']: # Nov24  # NB: this is meant to do together with allIDs/allIDs_DGemtFine_M1_*bn1369  NOT _n1521 
        val['M1'] = [100,140,250]
        val['M2'] = [100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200, 210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,400,450,500,550,600,650,700]
        val['mu'] = list(val['M2'])
        # quicky for binsL (but is anyway replaced if 'binwidth' is set
        binsL['M2'] = [95]
        for i in range(len(val['M2'])-1): binsL['M2'].append( (val['M2'][i]+val['M2'][i+1])/2 )
        binsL['M2'].append(725)
        binsL['mu'] = list(binsL['M2'])

    if sliceID in ['DGemtFine_to500']: # 
        val['M1'] = [100,140,250]
        val['M2'] = [100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200, 210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,400,450,500]
        val['mu'] = list(val['M2'])
        # quicky for binsL (but is anyway replaced if 'binwidth' is set
        binsL['M2'] = [95]
        for i in range(len(val['M2'])-1): binsL['M2'].append( (val['M2'][i]+val['M2'][i+1])/2 )
        binsL['M2'].append(525)
        binsL['mu'] = list(binsL['M2'])


    if sliceID in ['M2xMU=100-500newDGnoSL']: # 2013-06-26
        val['M1'] = [45]  # anything
        #val['M2'] = [100,120,150,175,200,225,250,275,300,350,400,450,500,999]
        val['M2'] = [100,120,150,175,200,225,250,275,300,350,400,450,500]
        val['M2'] = [100,120,150,175,200,225,250,275,300,350,400,450,500, 550]
        val['mu'] = list(val['M2'])
        binsL['M2'] = [40]
        for i in range(len(val['M2'])-1): binsL['M2'].append( (val['M2'][i]+val['M2'][i+1])/2 )
        binsL['M2'].append(val['M2'][-1] + 25)
        binsL['mu'] = list(binsL['M2'])

    if figname_class: fnn0 += '_' + figname_class


    if 1: # 2012-10-05  # 2013-07-14: available in isa_simplereader::GetVarText
        if v1T == 'M1 [GeV]': v1T = 'M_{1} [GeV]'
        if v1T == 'M2 [GeV]': v1T = 'M_{2} [GeV]'
        if v1T == 'MU [GeV]': v1T = '#mu [GeV]'
        if v2T == 'M1 [GeV]': v2T = 'M_{1} [GeV]'
        if v2T == 'M2 [GeV]': v2T = 'M_{2} [GeV]'
        if v2T == 'MU [GeV]': v2T = '#mu [GeV]'


    setup['val'] = val
    setup['binsL'] = binsL
    setup['v1'] = v1
    setup['v2'] = v2
    setup['V1'] = V1
    setup['V2'] = V2
    setup['v1T'] = v1T
    setup['v2T'] = v2T
    setup['fnn0'] = fnn0

    # print 'SETUP: ',setup
    return setup



# =================================
def get_ipro_values(ipro, yis, fnn0='', keyID='none', OptDict={}, yis2=-1, yis3=-1, dict2={}, dict3={}): 

    # print 'start ipro: fnn0 = %s' %(fnn0)

    _ZERO = 0.  # values used to make safer


    if type(yis) != int: print('WARNING  libIsaplot::get_ipro_values  type(yis) is ', type(yis))

    res = {'status':'ok'}
    y = 0.
    htit = ''
    fnn = ''
    hMax = -998.
    hMin = -999.
        

    if 1: # (to stay with the original indentation)

        ph  = ipro.br.part['h']
        pN2 = ipro.br.part['N2']
        pC1 = ipro.br.part['C1']
        pN3 = ipro.br.part['N3']
        
        pT1 = ipro.br.part['T1']

        pgl = ipro.br.part['gl']
        pdL = ipro.br.part['dL']
        puL = ipro.br.part['uL']
        pdR = ipro.br.part['dR']
        puR = ipro.br.part['uR']
        


        # ===================================================== y-values start
        # ===================================================== y-values start
        # ===================================================== y-values start

        
        # if type(yis) is int:  # 2012-06-04 
        # try: 
            
        # OVERVIEW
        # ->100    XSECs
        #   100->  Masses
        #   200->  BRs


        # -------------------------------------------------- 1 - 99
        if 1 <= yis <= 99: 

            if 50 <= yis <= 73 and not ipro.oxok: return {'status':'continue','warn':"50 <= yis <= 73 and not ipro.oxok"}
            if IsInRangesT('10-15,20-49,72',yis):
                try: ipro.play.ok
                except:
                    #print 'No play object'
                    return {'status':'continue','warn':"IsInRangesT('10-15,20-49,72',yis) and not ipro.play.ok"}  # 2011-11-14


            # N-LEP-XSEC (RELATIVE == nL)
            if yis == 10: 
                htit = "Event fraction with == 0 leptons"; fnn = "%s_rel_eq0L" %(fnn0)
                y = ipro.play.xrel[0] / ipro.play.xrelAccu[0]
            
            if yis == 11:
                htit = "Event fraction with == 1 leptons"; fnn = "%s_rel_eq1L" %(fnn0)
                y =  ipro.play.xrel[1] / ipro.play.xrelAccu[0]
            
            if yis == 12: 
                htit = "Event fraction with == 2 leptons"; fnn = "%s_rel_eq2L" %(fnn0)
                y = ipro.play.xrel[2] / ipro.play.xrelAccu[0]

            if yis == 13: 
                htit = "Event fraction with == 3 leptons"; fnn = "%s_rel_eq3L" %(fnn0)
                y = ipro.play.xrel[3] / ipro.play.xrelAccu[0]

            if yis == 14: 
                htit = "Event fraction with == 4 leptons"; fnn = "%s_rel_eq4L" %(fnn0)
                y = ipro.play.xrel[4] / ipro.play.xrelAccu[0]

            if yis == 15: 
                htit = "Event fraction with == 5 leptons"; fnn = "%s_rel_eq5L" %(fnn0)
                y = ipro.play.xrel[5] / ipro.play.xrelAccu[0]


            # N-LEP-XSEC (RELATIVE >= nL)

            if yis == 20: 
                htit = "Event fraction with #geq 0 leptons"; fnn = "%s_rel_ge0L" %(fnn0)
                y = ipro.play.xrelAccu[0] / ipro.play.xrelAccu[0]

            if yis == 21: 
                htit = "Event fraction with #geq 1 leptons"; fnn = "%s_rel_ge1L" %(fnn0)
                y = ipro.play.xrelAccu[1] / ipro.play.xrelAccu[0]

            if yis == 22: 
                htit = "Event fraction with #geq 2 leptons"; fnn = "%s_rel_ge2L" %(fnn0)
                y = ipro.play.xrelAccu[2] / ipro.play.xrelAccu[0]

            if yis == 23: 
                htit = "Event fraction with #geq 3 leptons"; fnn = "%s_rel_ge3L" %(fnn0)
                y = ipro.play.xrelAccu[3] / ipro.play.xrelAccu[0]

            if yis == 24: 
                htit = "Event fraction with #geq 4 leptons"; fnn = "%s_rel_ge4L" %(fnn0)
                y = ipro.play.xrelAccu[4] / ipro.play.xrelAccu[0]

            #if yis == 25: 
            #    htit = "Event fraction with #geq 5 leptons"; fnn = "%s_rel_ge5L" %(fnn0)
            #    y = ipro.play.xrelAccu[5] / ipro.play.xrelAccu[0]



            if yis in [25,26,27,28,29]:
                lumiT = '%.1f' %(yis2)
                if yis2 == int(yis2): lumiT = '%.0f' %(yis2)
                
            if yis == 25:
                htit = "nev at %s fb^{-1} with #geq 0 leptons" %(lumiT); fnn = "%s_nev%sifb_ge0L" %(fnn0,lumiT); y = ipro.play.xAccu[0] * 1000 * yis2

            if yis == 26: 
                htit = "nev at %s fb^{-1} with #geq 1 leptons" %(lumiT); fnn = "%s_nev%sifb_ge1L" %(fnn0,lumiT); y = ipro.play.xAccu[1] * 1000 * yis2

            if yis == 27: 
                htit = "nev at %s fb^{-1} with #geq 2 leptons" %(lumiT); fnn = "%s_nev%sifb_ge2L" %(fnn0,lumiT); y = ipro.play.xAccu[2] * 1000 * yis2

            if yis == 28: 
                htit = "nev at %s fb^{-1} with #geq 3 leptons" %(lumiT); fnn = "%s_nev%sifb_ge3L" %(fnn0,lumiT); y = ipro.play.xAccu[3] * 1000 * yis2

            if yis == 29: 
                htit = "nev at %s fb^{-1} with #geq 4 leptons" %(lumiT); fnn = "%s_nev%sifb_ge4L" %(fnn0,lumiT); y = ipro.play.xAccu[4] * 1000 * yis2


            #if yis == 35: 
            #    htit = "Cross-section [pb] with == 5 leptons"; fnn = "%s_xsec_eq5L" %(fnn0); y = ipro.play.x[5] 
            
            if yis == 30: 
                htit = "Cross-section [pb] with #geq 0 leptons"; fnn = "%s_xsec_ge0L" %(fnn0); y = ipro.play.xAccu[0] 
            if yis == 31: 
                htit = "Cross-section [pb] with #geq 1 leptons"; fnn = "%s_xsec_ge1L" %(fnn0); y = ipro.play.xAccu[1] 
            if yis == 32: 
                htit = "Cross-section [pb] with #geq 2 leptons"; fnn = "%s_xsec_ge2L" %(fnn0); y = ipro.play.xAccu[2] 
            if yis == 33: 
                htit = "Cross-section [pb] with #geq 3 leptons"; fnn = "%s_xsec_ge3L" %(fnn0); y = ipro.play.xAccu[3] 
            if yis == 34: 
                htit = "Cross-section [pb] with #geq 4 leptons"; fnn = "%s_xsec_ge4L" %(fnn0); y = ipro.play.xAccu[4] 

            if yis == 39: 
                htit = "Cross-section [pb]"; fnn = "%s_xsec_tot" %(fnn0); y = ipro.play.xAccu[0] 

            if yis == 40: 
                htit = "Cross-section [pb] with #geq 0 leptons"; fnn = "%s_xsec_ge0L" %(fnn0); y = ipro.play.xAccu[0] 
            if yis == 41: 
                htit = "Cross-section [pb] with #geq 1 leptons"; fnn = "%s_xsec_ge1L" %(fnn0); y = ipro.play.xAccu[1] 
            if yis == 42: 
                htit = "Cross-section [pb] with #geq 2 leptons"; fnn = "%s_xsec_ge2L" %(fnn0); y = ipro.play.xAccu[2] 
            if yis == 43: 
                htit = "Cross-section [pb] with #geq 3 leptons"; fnn = "%s_xsec_ge3L" %(fnn0); y = ipro.play.xAccu[3] 
            if yis == 44: 
                htit = "Cross-section [pb] with #geq 4 leptons"; fnn = "%s_xsec_ge4L" %(fnn0); y = ipro.play.xAccu[4] 

            #if yis == 45: 
            #    htit = "Cross-section [pb] with #geq 5 leptons"; fnn = "%s_xsec_ge5L" %(fnn0); y = ipro.play.xAccu[5] 
            
            # 2012-10-06: reordered these to allow for 0L and 1L in fb
            if yis == 45: 
                htit = "Cross-section [fb] with #geq 0 leptons"; fnn = "%s_xsec_fb_ge0L" %(fnn0); y = ipro.play.xAccu[0] * 1000. 
            if yis == 46: 
                htit = "Cross-section [fb] with #geq 1 leptons"; fnn = "%s_xsec_fb_ge1L" %(fnn0); y = ipro.play.xAccu[1] * 1000. 
            if yis == 47: 
                htit = "Cross-section [fb] with #geq 2 leptons"; fnn = "%s_xsec_fb_ge2L" %(fnn0); y = ipro.play.xAccu[2] * 1000. 
            if yis == 48: 
                htit = "Cross-section [fb] with #geq 3 leptons"; fnn = "%s_xsec_fb_ge3L" %(fnn0); y = ipro.play.xAccu[3] * 1000. 
            if yis == 49: 
                htit = "Cross-section [fb] with #geq 4 leptons"; fnn = "%s_xsec_fb_ge4L" %(fnn0); y = ipro.play.xAccu[4] * 1000. 


            ### N-LEP-XSEC (RELATIVE == nL)
            #print ipro.ox.x.keys()

            #if 50 <= yis <= 72:
            #    try: ipro.ox
            #    except:
            #        print "Seems like ipro.ox problem ",  ipro.ID, iSS
                    #continue
            
            if yis == 50: 
                htit = "#sigma(#tilde{N}_{1}#tilde{N}_{1}) / #sigma(#tilde{#chi}#tilde{#chi})"
                y = ipro.ox.x.get('N1 N1',_ZERO) / ipro.ox.xComb['X X'];  fnn = "%s_xsec_N1N1onXX" %(fnn0);
                '''
                try:
                    y = ipro.ox.x['N1 N1'] / ipro.ox.xComb['X X'];  fnn = "%s_xsec_N1N1onXX" %(fnn0);
                    print 'yis=50: success', ipro.ID, iSS
                except:
                    print "Seems like ipro.ox problem; may crash ; yis=50  ", ipro.ID, iSS
                    y = -1.
                '''
                    
            if yis == 51: 
                htit = "#sigma(#tilde{N}_{1}#tilde{C}_{1}) / #sigma(#tilde{#chi}#tilde{#chi})"
                # print iSS, ipro.ID, ipro.ox.x.keys()
                #y = ipro.ox.x['N1 C1'];  fnn = "%s_xsec_N1C1onXX" %(fnn0);
                y = ipro.ox.x.get('N1 C1',_ZERO) / ipro.ox.xComb['X X'];  fnn = "%s_xsec_N1C1onXX" %(fnn0);

            if yis == 52: 
                htit = "#sigma(#tilde{N}_{1}#tilde{N}_{2}) / #sigma(#tilde{#chi}#tilde{#chi})"
                y = ipro.ox.x['N1 N2'] / ipro.ox.xComb['X X'];  fnn = "%s_xsec_N1N2onXX" %(fnn0);

            if yis == 53: 
                htit = "#sigma(#tilde{C}_{1}#tilde{C}_{1}) / #sigma(#tilde{#chi}#tilde{#chi})"
                y = ipro.ox.x.get('C1 C1',_ZERO) / ipro.ox.xComb['X X'];  fnn = "%s_xsec_C1C1onXX" %(fnn0);

            if yis == 54: 
                htit = "#sigma(#tilde{N}_{2}#tilde{C}_{1}) / #sigma(#tilde{#chi}#tilde{#chi})"
                y = ipro.ox.x['N2 C1'] / ipro.ox.xComb['X X'];  fnn = "%s_xsec_N2C1onXX" %(fnn0);

            if yis == 55: 
                htit = "#sigma(#tilde{N}_{3}#tilde{C}_{1}) / #sigma(#tilde{#chi}#tilde{#chi})"
                y = ipro.ox.x['N3 C1'] / ipro.ox.xComb['X X'];  fnn = "%s_xsec_N3C1onXX" %(fnn0);
                
            if yis == 56: 
                htit = "#sigma(#tilde{N}_{2}#tilde{N}_{2}) / #sigma(#tilde{#chi}#tilde{#chi})"
                y = ipro.ox.x['N2 N2'] / ipro.ox.xComb['X X'];  fnn = "%s_xsec_N2N2onXX" %(fnn0);

            if yis == 57: 
                htit = "#sigma(#tilde{N}_{4}#tilde{C}_{1}) / #sigma(#tilde{#chi}#tilde{#chi})"
                y = ipro.ox.x['N4 C1'] / ipro.ox.xComb['X X'];  fnn = "%s_xsec_N4C1onXX" %(fnn0);

            if yis == 58:
                # print ipro.ox.x.keys()
                htit = "#sigma(#tilde{N}_{i}#tilde{C}_{2}) / #sigma(#tilde{#chi}#tilde{#chi})"
                y = (ipro.ox.x['N1 C2'] + ipro.ox.x['N2 C2'] + ipro.ox.x['N3 C2'] + ipro.ox.x['N4 C2']) / ipro.ox.xComb['X X'];  fnn = "%s_xsec_NiC2onXX" %(fnn0);

            if yis == 61: 
                htit = "#sigma(#tilde{C}_{1,2}#tilde{C}_{2}) / #sigma(#tilde{#chi}#tilde{#chi})"
                #y = 0.01;  fnn = "%s_xsec_N4C1onXX" %(fnn0);
                y = (ipro.ox.x['C1 C2'] + ipro.ox.x['C2 C2']) / ipro.ox.xComb['X X'];  fnn = "%s_xsec_N4C1onXX" %(fnn0);

            if yis == 62: 
                htit = "#sigma(#tilde{N}_{34}#tilde{N}_{34}) / #sigma(#tilde{#chi}#tilde{#chi})"
                #print ipro.ox.x.keys()
                y = (ipro.ox.x['N3 N3'] + ipro.ox.x['N3 N4'] + ipro.ox.x['N4 N4']) / ipro.ox.xComb['X X'];  fnn = "%s_xsec_N34N34onXX" %(fnn0);

            if yis == 63: 
                htit = "#sigma(#tilde{N}_{12}#tilde{N}_{34}) / #sigma(#tilde{#chi}#tilde{#chi})"
                y = (ipro.ox.x['N1 N3'] + ipro.ox.x['N1 N4'] + ipro.ox.x['N2 N3'] + ipro.ox.x['N2 N4']) / ipro.ox.xComb['X X'];  fnn = "%s_xsec_N12N34onXX" %(fnn0);

            if yis == 64: 
                htit = "#sigma(all with heavy) / #sigma(#tilde{#chi}#tilde{#chi})"
                y = (ipro.ox.x['N1 N3'] + ipro.ox.x['N1 N4'] + ipro.ox.x['N2 N3'] + ipro.ox.x['N2 N4'] + ipro.ox.x['N3 N3'] + ipro.ox.x['N3 N4'] + ipro.ox.x['N4 N4'] + ipro.ox.x['C1 C2'] + ipro.ox.x['C2 C2'] + ipro.ox.x['N1 C2'] + ipro.ox.x['N2 C2'] + ipro.ox.x['N3 C2'] + ipro.ox.x['N4 C2'] + ipro.ox.x['N3 C1'] + ipro.ox.x['N4 C1']) / ipro.ox.xComb['X X'];  fnn = "%s_xsec_NCwithheavyonXX" %(fnn0);

            if yis == 65: 
                htit = "#sigma(#tilde{N}_{3,4}#tilde{C}_{1}) / #sigma(#tilde{#chi}#tilde{#chi})"
                #y = 0.01;  fnn = "%s_xsec_N4C1onXX" %(fnn0);
                y = (ipro.ox.x['N3 C1'] + ipro.ox.x['N4 C1']) / ipro.ox.xComb['X X'];  fnn = "%s_xsec_N34C1onXX" %(fnn0);


            if yis == 59: # special dirty one: the amount of 3L from N2C1, i.e.: relN2C1 * BR(N2->ll) * BR(C1->l)
                htit = "(Total) Event fraction ==3L from #tilde{N}_{2}#tilde{C}_{1} alone"
                #y = ( ipro.ox.x['N2 C1'] / ipro.ox.xComb['X X'] ) * ( pN2.BR('eR') + pN2.BR('eL') + pN2.BR('N1','e') + pN2.BR('T1') + pN2.BR('T2') + pN2.BR('N1','ta') ) * ( pC1.BR('ve') + pC1.BR('eR') + pC1.BR('eL') + pC1.BR('N1','nue') + ( pC1.BR('vT') + pC1.BR('T1') + pC1.BR('T2') + pC1.BR('N1','nut') ) * mypdg.TtoL + pC1.BR('N1','W') * mypdg.WtoLtot )
                zbrN2 = ( pN2.BR('eR') + pN2.BR('eL') + pN2.BR('N1','e') + ( pN2.BR('T1') + pN2.BR('T2') + pN2.BR('N1','ta') ) * mypdg.TtoL * mypdg.TtoL + pN2.BR('N1','Z') * mypdg.ZtoLLtot )
                zbrC1 = ( pC1.BR('ve') + pC1.BR('eR') + pC1.BR('eL') + pC1.BR('N1','nue') + ( pC1.BR('vT') + pC1.BR('T1') + pC1.BR('T2') + pC1.BR('N1','nut') ) * mypdg.TtoL + pC1.BR('N1','W') * mypdg.WtoLtot )
                y = zbrN2 * zbrC1 
                if 0: 
                    print('y59: ', y, brN2, brC1, pN2.BR('N1','e'), pN2.BR('N1','ta'), brN2*brC1)
                    

                fnn = "%s_xsec_N3C1onXX_to3L" %(fnn0);

            #x_LL

            if yis == 70:
                htit = "#sigma(ll) [pb]"
                y = ( ipro.ox.getx('eR eR') + ipro.ox.getx('mR mR') + ipro.ox.getx('T1 T1') + ipro.ox.getx('eL eL') + ipro.ox.getx('mL mL') + ipro.ox.getx('T2 T2')  + ipro.ox.getx('T1 T2') )  
                fnn = "%s_xsec_ll" %(fnn0)
                
            if yis == 71:
                htit = "#sigma(ll) / #sigma(#tilde{#chi}#tilde{#chi})"
                y = ( ipro.ox.getx('eR eR') + ipro.ox.getx('mR mR') + ipro.ox.getx('T1 T1') + ipro.ox.getx('eL eL') + ipro.ox.getx('mL mL') + ipro.ox.getx('T2 T2')  + ipro.ox.getx('T1 T2') ) / (ipro.ox.xComb['X X'])
                fnn = "%s_xsec_llonXX" %(fnn0)


            if yis == 72:  # NB: requires play 
                htit = "#sigma(ll) / #sigma_{2L}(#tilde{#chi}#tilde{#chi})"
                den = ipro.play.xAccu[2]
                y = 0.
                if den > 0: y = ( ipro.ox.getx('eR eR') + ipro.ox.getx('mR mR') + ipro.ox.getx('T1 T1') + ipro.ox.getx('eL eL') + ipro.ox.getx('mL mL') + ipro.ox.getx('T2 T2')  + ipro.ox.getx('T1 T2') ) / den
                fnn = "%s_xsec_llonXXge2L" %(fnn0)

                
            if yis == 73: #naming ... this is LL rather than ll (filename, also above)
                htit = "#sigma(#tilde{L}#tilde{L}) / ( #sigma(#tilde{#chi}#tilde{#chi}) + #sigma(#tilde{L}#tilde{L}) )"
                num = ( ipro.ox.getx('eR eR') + ipro.ox.getx('mR mR') + ipro.ox.getx('T1 T1') + ipro.ox.getx('eL eL') + ipro.ox.getx('mL mL') + ipro.ox.getx('T2 T2')  + ipro.ox.getx('T1 T2') )
                y = num / (num + ipro.ox.xComb['X X'])               
                fnn = "%s_xsec_llonXL" %(fnn0)


            # szMar = 1.0; gStyle.SetPaintTextFormat('.1f'); 
            if 0:
                gStyle.SetPaintTextFormat('.0f'); 
                htit += " (in %)"
                hMax *= 100
                hMax = 100
                DRAWTEXT = 1 #2
                y *= 100



        # -------------------------------------------------- 101 - 199
        if 701 <= yis <= 799: 
            if not ipro.oxok: return {'status':'continue'}  # 2011-11-14
            
            # DIRECT NEUTRALINOS/CHARGINOS
            # print ipro.ox.x.keys()
            
            if yis == 701:
                y = ipro.ox.x.get('N1 N1',_ZERO); htit = '#sigma(#tilde{N}_{1}#tilde{N}_{1})  [pb]'; fnn = "%s_xsecN1N1" %(fnn0)
            if yis == 702:
                y = ipro.ox.x['N1 N2']; htit = '#sigma(#tilde{N}_{1}#tilde{N}_{2})  [pb]'; fnn = "%s_xsecN1N2" %(fnn0)
            if yis == 703:
                y = ipro.ox.x['N1 N3']; htit = '#sigma(#tilde{N}_{1}#tilde{N}_{3})  [pb]'; fnn = "%s_xsecN1N3" %(fnn0)
            if yis == 704:
                y = ipro.ox.x['N1 N4']; htit = '#sigma(#tilde{N}_{1}#tilde{N}_{4})  [pb]'; fnn = "%s_xsecN1N4" %(fnn0)
            if yis == 705:
                y = ipro.ox.x.get('N1 C1',_ZERO); htit = '#sigma(#tilde{N}_{1}#tilde{C}_{1})  [pb]'; fnn = "%s_xsecN1C1" %(fnn0)
            if yis == 706:
                y = ipro.ox.x['N1 C2']; htit = '#sigma(#tilde{N}_{1}#tilde{C}_{2})  [pb]'; fnn = "%s_xsecN1C2" %(fnn0)

            if yis == 707:
                y = ipro.ox.x['N2 N2']; htit = '#sigma(#tilde{N}_{2}#tilde{N}_{2})  [pb]'; fnn = "%s_xsecN2N2" %(fnn0)
            if yis == 708:
                y = ipro.ox.x['N2 N3']; htit = '#sigma(#tilde{N}_{2}#tilde{N}_{3})  [pb]'; fnn = "%s_xsecN2N3" %(fnn0)
            if yis == 709:
                y = ipro.ox.x['N2 N4']; htit = '#sigma(#tilde{N}_{2}#tilde{N}_{4})  [pb]'; fnn = "%s_xsecN2N4" %(fnn0)
            if yis == 710:
                y = ipro.ox.x['N2 C1']; htit = '#sigma(#tilde{N}_{2}#tilde{C}_{1})  [pb]'; fnn = "%s_xsecN2C1" %(fnn0)
            if yis == 711:
                y = ipro.ox.x['N2 C2']; htit = '#sigma(#tilde{N}_{2}#tilde{C}_{2})  [pb]'; fnn = "%s_xsecN2C2" %(fnn0)

            if yis == 712:
                y = ipro.ox.x['N3 N3']; htit = '#sigma(#tilde{N}_{3}#tilde{N}_{3})  [pb]'; fnn = "%s_xsecN3N3" %(fnn0)
            if yis == 713:
                y = ipro.ox.x['N3 N4']; htit = '#sigma(#tilde{N}_{3}#tilde{N}_{4})  [pb]'; fnn = "%s_xsecN3N4" %(fnn0)
            if yis == 714:
                y = ipro.ox.x['N3 C1']; htit = '#sigma(#tilde{N}_{3}#tilde{C}_{1})  [pb]'; fnn = "%s_xsecN3C1" %(fnn0)
            if yis == 715:
                y = ipro.ox.x['N3 C2']; htit = '#sigma(#tilde{N}_{3}#tilde{C}_{2})  [pb]'; fnn = "%s_xsecN3C2" %(fnn0)

            if yis == 716:
                y = ipro.ox.x['N4 N4']; htit = '#sigma(#tilde{N}_{4}#tilde{N}_{4})  [pb]'; fnn = "%s_xsecN4N4" %(fnn0)
            if yis == 717:
                y = ipro.ox.x['N4 C1']; htit = '#sigma(#tilde{N}_{4}#tilde{C}_{1})  [pb]'; fnn = "%s_xsecN4C1" %(fnn0)
            if yis == 718:
                y = ipro.ox.x['N4 C2']; htit = '#sigma(#tilde{N}_{4}#tilde{C}_{2})  [pb]'; fnn = "%s_xsecN4C2" %(fnn0)

            if yis == 719:
                y = ipro.ox.x.get('C1 C1',_ZERO); htit = '#sigma(#tilde{C}_{1}#tilde{C}_{1})  [pb]'; fnn = "%s_xsecC1C1" %(fnn0)
            if yis == 720:
                y = ipro.ox.x['C1 C2']; htit = '#sigma(#tilde{C}_{1}#tilde{C}_{2})  [pb]'; fnn = "%s_xsecC1C2" %(fnn0)

            if yis == 721:
                #print ipro.ox.x.keys()
                y = ipro.ox.x['C2 C2']; htit = '#sigma(#tilde{C}_{2}#tilde{C}_{2})  [pb]'; fnn = "%s_xsecC2C2" %(fnn0)

                
            if yis == 722:
                y = ipro.ox.Sum(['N2 C1','N2 C2','N3 C1','N3 C2','N4 C1','N4 C2']); htit = '#sigma(#tilde{N}_{234}#tilde{C}_{12})  [pb]'; fnn = "%s_xsecN234C12" %(fnn0)

            ### XX total sums
            if yis == 723:
                y = ipro.ox.xComb['X X']; htit = '#sigma(#tilde{#chi}#tilde{#chi})  [pb]'; fnn = "%s_xsecXX" %(fnn0) 
                
            if yis == 724:
                y = ipro.ox.xComb['X X'] - ipro.ox.x.get('N1 N1',_ZERO); htit = '#sigma(#tilde{#chi}#tilde{#chi} - #tilde{N}_{1}#tilde{N}_{1})  [pb]'; fnn = "%s_xsecXXa" %(fnn0)
                
            if yis == 725:
                y = ipro.ox.xComb['X X'] - ipro.ox.x.get('N1 N1',_ZERO) - ipro.ox.x.get('N1 C1',_ZERO); htit = '#sigma(#tilde{#chi}#tilde{#chi} - #tilde{N}_{1}#tilde{N}_{1} - #tilde{N}_{1}#tilde{C}_{1})  [pb]'; fnn = "%s_xsecXXb" %(fnn0) 
                

            ### XX partial sums
            if yis == 726: 
                #htit = "#sigma(all with heavy) / #sigma(#tilde{#chi}#tilde{#chi})"  # 2012-07-12 wrong title fixed
                htit = "#sigma(#tilde{#chi}#tilde{#chi} some heavy)  [pb]"
                y = ipro.ox.x['N1 N3'] + ipro.ox.x['N1 N4'] + ipro.ox.x['N2 N3'] + ipro.ox.x['N2 N4'] + ipro.ox.x['N3 N3'] + ipro.ox.x['N3 N4'] + ipro.ox.x['N4 N4'] + ipro.ox.x['C1 C2'] + ipro.ox.x['C2 C2'] + ipro.ox.x['N1 C2'] + ipro.ox.x['N2 C2'] + ipro.ox.x['N3 C2'] + ipro.ox.x['N4 C2'] + ipro.ox.x['N3 C1'] + ipro.ox.x['N4 C1'];  fnn = "%s_xsecXXwithheavy" %(fnn0);

            if yis == 727: 
                htit = "#sigma(#tilde{N}#tilde{N} some heavy)  [pb]"
                y = ipro.ox.x['N1 N3'] + ipro.ox.x['N1 N4'] + ipro.ox.x['N2 N3'] + ipro.ox.x['N2 N4'] + ipro.ox.x['N3 N3'] + ipro.ox.x['N3 N4'] + ipro.ox.x['N4 N4'] ; fnn = "%s_xsecNNsomeheavy" %(fnn0);

            if yis == 728: 
                htit = "#sigma(#tilde{C}#tilde{C} some heavy)  [pb]"
                y = ipro.ox.x['C1 C2'] + ipro.ox.x['C2 C2'] ; fnn = "%s_xsecCCsomeheavy" %(fnn0);

            if yis == 729: 
                htit = "#sigma(#tilde{N}#tilde{C} some heavy)  [pb]"
                y = ipro.ox.x['N1 C2'] + ipro.ox.x['N2 C2'] + ipro.ox.x['N3 C2'] + ipro.ox.x['N4 C2'] + ipro.ox.x['N3 C1'] + ipro.ox.x['N4 C1'];  fnn = "%s_xsecNCsomeheavy" %(fnn0);


            # might not be needed
            if yis == 730: 
                htit = "#sigma(#tilde{N}_{34}#tilde{N}_{34})  [pb]"
                y = (ipro.ox.x['N3 N3'] + ipro.ox.x['N3 N4'] + ipro.ox.x['N4 N4']);  fnn = "%s_xsecN34N34" %(fnn0);

            if yis == 731: 
                htit = "#sigma(#tilde{N}_{12}#tilde{N}_{34})  [pb]"
                y = (ipro.ox.x['N1 N3'] + ipro.ox.x['N1 N4'] + ipro.ox.x['N2 N3'] + ipro.ox.x['N2 N4']);  fnn = "%s_xsecN12N34" %(fnn0);

                

            # DIRECT SLEPTONS
            if yis == 741:
                y = 2 * ipro.ox.x.get('eL eL',0); htit = '#sigma(#tilde{l}_{L}#tilde{l}_{L})  [pb]'; fnn = "%s_xseclLlL" %(fnn0)
            if yis == 742:
                y = 2 * ipro.ox.x.get('eR eR',0); htit = '#sigma(#tilde{l}_{R}#tilde{l}_{R})  [pb]'; fnn = "%s_xseclRlR" %(fnn0)
            if yis == 743:
                y = 2 * ipro.ox.x.get('ve ve',0); htit = '#sigma(#tilde{v}_{l}#tilde{v}_{l})  [pb]'; fnn = "%s_xsecvlvl" %(fnn0)
            if yis == 744:
                y = 2 * ipro.ox.x.get('eL ve',0); htit = '#sigma(#tilde{l}_{L}#tilde{v}_{l})  [pb]'; fnn = "%s_xseclLvl" %(fnn0)

            if yis == 745:
                y = ipro.ox.x.get('T1 T1',0); htit = '#sigma(#tilde{#tau}_{1}#tilde{#tau}_{1})  [pb]'; fnn = "%s_xsecT1T1" %(fnn0)
            if yis == 746:
                y = ipro.ox.x.get('T1 T2',0); htit = '#sigma(#tilde{#tau}_{2}#tilde{#tau}_{2})  [pb]'; fnn = "%s_xsecT1T2" %(fnn0)
            if yis == 747:
                y = ipro.ox.x.get('T2 T2',0); htit = '#sigma(#tilde{#tau}_{2}#tilde{#tau}_{2})  [pb]'; fnn = "%s_xsecT2T2" %(fnn0)
            if yis == 748:
                y = ipro.ox.x.get('vT vT',0); htit = '#sigma(#tilde{#nu}_{#tau}#tilde{#nu}_{#tau})  [pb]'; fnn = "%s_xsecvTvT" %(fnn0)
            if yis == 749:
                y = ipro.ox.x.get('T1 vT',0); htit = '#sigma(#tilde{#tau}_{1}#tilde{#nu}_{#tau})  [pb]'; fnn = "%s_xsecT1vT" %(fnn0)
            if yis == 750:
                y = ipro.ox.x.get('T2 vT',0); htit = '#sigma(#tilde{#tau}_{2}#tilde{#nu}_{#tau})  [pb]'; fnn = "%s_xsecT2vT" %(fnn0)


            # Various slep combinations
            if yis == 751:
                y = 2 * (ipro.ox.x.get('eL eL',0) + ipro.ox.x.get('eR eR',0)); htit = '#sigma(#tilde{l}#tilde{l})  [pb]'; fnn = "%s_xsecslsl" %(fnn0)
            if yis == 752:
                y = ipro.ox.x.get('T1 T1',0) + ipro.ox.x.get('T1 T2',0) + ipro.ox.x.get('T2 T2',0); htit = '#sigma(#tilde{#tau}#tilde{#tau})  [pb]'; fnn = "%s_xsecTiTj" %(fnn0)

            if yis == 753:
                y = ipro.ox.x.get('T1 vT',0) + ipro.ox.x.get('T2 vT',0) ; htit = '#sigma(#tilde{#tau}#tilde{#nu}_{#tau})  [pb]'; fnn = "%s_xsecTivT" %(fnn0)

            if yis == 754:
                y = 2 * ipro.ox.x.get('eL eL',0) + 2 * ipro.ox.x.get('eR eR',0) + ipro.ox.x.get('T1 T1',0) + ipro.ox.x.get('T1 T2',0) + ipro.ox.x.get('T2 T2',0); htit = '#sigma(#tilde{L}#tilde{L})  [pb]'; fnn = "%s_xsecSLSL" %(fnn0)

            if yis == 755:
                y = 2 * ipro.ox.x.get('eL ve',0) + ipro.ox.x.get('T1 vT',0) + ipro.ox.x.get('T2 vT',0); htit = '#sigma(#tilde{L} #tilde{#nu})  [pb]'; fnn = "%s_xsecSLvL" %(fnn0)
  
            if yis == 756:
                y = 2 * ipro.ox.x.get('ve ve',0) + ipro.ox.x.get('vT vT',0); htit = '#sigma(#tilde{#nu} #tilde{#nu})  [pb]'; fnn = "%s_xsecvLvL" %(fnn0)

            if yis == 757: # LL
                # y = 2 * ( ipro.ox.x['eL eL'] + ipro.ox.x['eR eR'] + ipro.ox.x['ve ve'] + ipro.ox.x['eL ve'] ) + ipro.ox.x['T1 T1'] + ipro.ox.x['T1 T2'] + ipro.ox.x['T2 T2'] + ipro.ox.x['T1 vT'] + ipro.ox.x['T2 vT'] + ipro.ox.x['vT vT']; htit = '#sigma(direct sleptons)  [pb]'; fnn = "%s_xsecSlep" %(fnn0)
                y = 2 * (ipro.ox.Sum(['eL eL','eR eR','ve ve','eL ve']) ) + ipro.ox.Sum(['T1 T1','T1 T2','T2 T2','T1 vT','T2 vT','vT vT']); htit = '#sigma(DirSlep)  [pb]'; fnn = "%s_xsecdirslep" %(fnn0)


            ### XX and SLEP
            if yis == 758: # XX + LL
                # y = 2 * ( ipro.ox.x['eL eL'] + ipro.ox.x['eR eR'] + ipro.ox.x['ve ve'] + ipro.ox.x['eL ve'] ) + ipro.ox.x['T1 T1'] + ipro.ox.x['T1 T2'] + ipro.ox.x['T2 T2'] + ipro.ox.x['T1 vT'] + ipro.ox.x['T2 vT'] + ipro.ox.x['vT vT']; htit = '#sigma(direct sleptons)  [pb]'; fnn = "%s_xsecSlep" %(fnn0)
                y = ipro.ox.xComb['X X'] + 2 * (ipro.ox.Sum(['eL eL','eR eR','ve ve','eL ve']) ) + ipro.ox.Sum(['T1 T1','T1 T2','T2 T2','T1 vT','T2 vT','vT vT']); htit = '#sigma(#tilde{#chi}#tilde{#chi}+DirSlep)  [pb]'; fnn = "%s_xsecXL" %(fnn0)  #2012-07-12: fixed by adding XX




        # =============== MASSES
        msle_min = min(abs(ipro.mass['eR']),abs(ipro.mass['eL']),abs(ipro.mass['T1']),abs(ipro.mass['T2']))
        msnu_min = min(abs(ipro.mass['ve']),abs(ipro.mass['vT']))
        # mX_min = min(abs(ipro.mass['N1']),abs(ipro.mass['N2']),abs(ipro.mass['N3']),abs(ipro.mass['N4']),abs(ipro.mass['C1']),abs(ipro.mass['C2']))   # Not in use
        # -------------------------------------------------- 101 - 199
        if 101 <= yis <= 199: 

            if yis == 101: 
                htit = 'm(#tilde{N}_{1})  [GeV]'; fnn = "%s_mN1"%(fnn0); y = abs(ipro.mass['N1']); hMax = 100
            if yis == 102:
                htit = 'm(#tilde{N}_{2})  [GeV]'; fnn = "%s_mN2"%(fnn0); y = abs(ipro.mass['N2']); hMax = 100
                #print 'N2 fnn: %s' %(fnn)
            if yis == 103: 
                htit = 'm(#tilde{N}_{3})  [GeV]'; fnn = "%s_mN3"%(fnn0); y = abs(ipro.mass['N3']); hMax = 100
            if yis == 104: 
                htit = 'm(#tilde{N}_{4})  [GeV]'; fnn = "%s_mN4"%(fnn0); y = abs(ipro.mass['N4']); hMax = 100
            if yis == 105: 
                htit = 'm(#tilde{C}_{1})  [GeV]'; fnn = "%s_mC1"%(fnn0); y = abs(ipro.mass['C1']); hMax = 100
            if yis == 106: 
                htit = 'm(#tilde{C}_{2})  [GeV]'; fnn = "%s_mC2"%(fnn0); y = abs(ipro.mass['C2']); hMax = 100

            if yis == 121: 
                htit = 'm(#tilde{l}_{R})  [GeV]'; fnn = "%s_meR"%(fnn0); y = abs(ipro.mass['eR']); hMax = 100
            if yis == 122: 
                htit = 'm(#tilde{l}_{L})  [GeV]'; fnn = "%s_meL"%(fnn0); y = abs(ipro.mass['eL']); hMax = 100
            if yis == 123: 
                htit = 'm(#tilde{#nu}_{l})  [GeV]'; fnn = "%s_mvl"%(fnn0); y = abs(ipro.mass['ve']); hMax = 100
            if yis == 124: 
                htit = 'm(#tilde{#tau}_{1})  [GeV]'; fnn = "%s_mT1"%(fnn0); y = abs(ipro.mass['T1']); hMax = 100
            if yis == 125: 
                htit = 'm(#tilde{#tau}_{2})  [GeV]'; fnn = "%s_mT2"%(fnn0); y = abs(ipro.mass['T2']); hMax = 100
            if yis == 126: 
                htit = 'm(#tilde{#nu}_{#tau})  [GeV]'; fnn = "%s_mvT"%(fnn0); y = abs(ipro.mass['vT']); hMax = 100
                
            
            


        # -------------------------------------------------- 401 - 499
        if 401 <= yis <= 499:
            
            # === MASS DIFFERENCES
            if yis == 401:
                htit = 'm(#tilde{C}_{1}) - m(#tilde{N}_{1})  [GeV]'; fnn = "%s_C1mN1"%(fnn0); y = abs(ipro.mass['C1'])-abs(ipro.mass['N1']); hmax = 100
            if yis == 402:
                htit = 'm(#tilde{N}_{2}) - m(#tilde{N}_{1})  [GeV]'; fnn = "%s_N2mN1"%(fnn0); y = abs(ipro.mass['N2'])-abs(ipro.mass['N1']); hmax = 100
            if yis == 403:
                htit = 'm(#tilde{N}_{3}) - m(#tilde{N}_{1})  [GeV]'; fnn = "%s_N3mN1"%(fnn0); y = abs(ipro.mass['N3'])-abs(ipro.mass['N1']); hmax = 100
            if yis == 404:
                htit = 'm(#tilde{N}_{4}) - m(#tilde{N}_{1})  [GeV]'; fnn = "%s_N4mN1"%(fnn0); y = abs(ipro.mass['N4'])-abs(ipro.mass['N1']); hmax = 100
            if yis == 405:
                htit = 'm(#tilde{C}_{2}) - m(#tilde{N}_{1})  [GeV]'; fnn = "%s_C2mN1"%(fnn0); y = abs(ipro.mass['C2'])-abs(ipro.mass['N1']); hmax = 100

            if yis == 406:
                htit = 'm(#tilde{N}_{2}) - m(#tilde{C}_{1})  [GeV]'; fnn = "%s_N2mC1"%(fnn0); y = abs(ipro.mass['N2'])-abs(ipro.mass['C1']); hmax = 100
            if yis == 407:
                htit = 'm(#tilde{N}_{3}) - m(#tilde{C}_{1})  [GeV]'; fnn = "%s_N3mC1"%(fnn0); y = abs(ipro.mass['N3'])-abs(ipro.mass['C1']); hmax = 100
            if yis == 408:
                htit = 'm(#tilde{N}_{4}) - m(#tilde{C}_{1})  [GeV]'; fnn = "%s_N4mC1"%(fnn0); y = abs(ipro.mass['N4'])-abs(ipro.mass['C1']); hmax = 100
            if yis == 409:
                htit = 'm(#tilde{C}_{2}) - m(#tilde{C}_{1})  [GeV]'; fnn = "%s_C2mC1"%(fnn0); y = abs(ipro.mass['C2'])-abs(ipro.mass['C1']); hmax = 100

            if yis == 410:
                htit = 'm(#tilde{N}_{3}) - m(#tilde{N}_{2})  [GeV]'; fnn = "%s_N3mN2"%(fnn0); y = abs(ipro.mass['N3'])-abs(ipro.mass['N2']); hmax = 100
            if yis == 411:
                htit = 'm(#tilde{N}_{4}) - m(#tilde{N}_{2})  [GeV]'; fnn = "%s_N4mN2"%(fnn0); y = abs(ipro.mass['N4'])-abs(ipro.mass['N2']); hmax = 100
            if yis == 412:
                htit = 'm(#tilde{C}_{2}) - m(#tilde{N}_{2})  [GeV]'; fnn = "%s_C2mN2"%(fnn0); y = abs(ipro.mass['C2'])-abs(ipro.mass['N2']); hmax = 100


            if yis == 413:
                htit = 'm(#tilde{C}_{1}) - m(#tilde{l}_{R})  [GeV]'; fnn = "%s_C1mlR"%(fnn0); y = abs(ipro.mass['C1'])-abs(ipro.mass['eR']); hmax = 100
            if yis == 414:
                htit = 'm(#tilde{N}_{2}) - m(#tilde{l}_{R})  [GeV]'; fnn = "%s_N2mlR"%(fnn0); y = abs(ipro.mass['N2'])-abs(ipro.mass['eR']); hmax = 100
            if yis == 415:
                htit = 'm(#tilde{N}_{3}) - m(#tilde{l}_{R})  [GeV]'; fnn = "%s_N3mlR"%(fnn0); y = abs(ipro.mass['N3'])-abs(ipro.mass['eR']); hmax = 100
            if yis == 416:
                htit = 'm(#tilde{N}_{4}) - m(#tilde{l}_{R})  [GeV]'; fnn = "%s_N4mlR"%(fnn0); y = abs(ipro.mass['N4'])-abs(ipro.mass['eR']); hmax = 100
            if yis == 417:
                htit = 'm(#tilde{C}_{2}) - m(#tilde{l}_{R})  [GeV]'; fnn = "%s_C2mlR"%(fnn0); y = abs(ipro.mass['C2'])-abs(ipro.mass['eR']); hmax = 100

            if yis == 418:
                htit = 'm(#tilde{C}_{1}) - m(#tilde{l}_{L})  [GeV]'; fnn = "%s_C1mlL"%(fnn0); y = abs(ipro.mass['C1'])-abs(ipro.mass['eL']); hmax = 100
            if yis == 419:
                htit = 'm(#tilde{N}_{2}) - m(#tilde{l}_{L})  [GeV]'; fnn = "%s_N2mlL"%(fnn0); y = abs(ipro.mass['N2'])-abs(ipro.mass['eL']); hmax = 100
            if yis == 420:
                htit = 'm(#tilde{N}_{3}) - m(#tilde{l}_{L})  [GeV]'; fnn = "%s_N3mlL"%(fnn0); y = abs(ipro.mass['N3'])-abs(ipro.mass['eL']); hmax = 100
            if yis == 421:
                htit = 'm(#tilde{N}_{4}) - m(#tilde{l}_{L})  [GeV]'; fnn = "%s_N4mlL"%(fnn0); y = abs(ipro.mass['N4'])-abs(ipro.mass['eL']); hmax = 100
            if yis == 422:
                htit = 'm(#tilde{C}_{2}) - m(#tilde{l}_{L})  [GeV]'; fnn = "%s_C2mlL"%(fnn0); y = abs(ipro.mass['C2'])-abs(ipro.mass['eL']); hmax = 100

            if yis == 423:
                htit = 'm(#tilde{C}_{1}) - m(#tilde{#tau}_{1})  [GeV]'; fnn = "%s_C1mT1"%(fnn0); y = abs(ipro.mass['C1'])-abs(ipro.mass['T1']); hmax = 100
            if yis == 424:
                htit = 'm(#tilde{N}_{2}) - m(#tilde{#tau}_{1})  [GeV]'; fnn = "%s_N2mT1"%(fnn0); y = abs(ipro.mass['N2'])-abs(ipro.mass['T1']); hmax = 100
            if yis == 425:
                htit = 'm(#tilde{N}_{3}) - m(#tilde{#tau}_{1})  [GeV]'; fnn = "%s_N3mT1"%(fnn0); y = abs(ipro.mass['N3'])-abs(ipro.mass['T1']); hmax = 100
            if yis == 426:
                htit = 'm(#tilde{N}_{4}) - m(#tilde{#tau}_{1})  [GeV]'; fnn = "%s_N4mT1"%(fnn0); y = abs(ipro.mass['N4'])-abs(ipro.mass['T1']); hmax = 100
            if yis == 427:
                htit = 'm(#tilde{C}_{2}) - m(#tilde{#tau}_{1})  [GeV]'; fnn = "%s_C2mT1"%(fnn0); y = abs(ipro.mass['C2'])-abs(ipro.mass['T1']); hmax = 100

            if yis == 428:
                htit = 'm(#tilde{C}_{1}) - m(#tilde{#tau}_{2})  [GeV]'; fnn = "%s_C1mT2"%(fnn0); y = abs(ipro.mass['C1'])-abs(ipro.mass['T2']); hmax = 100
            if yis == 429:
                htit = 'm(#tilde{N}_{2}) - m(#tilde{#tau}_{2})  [GeV]'; fnn = "%s_N2mT2"%(fnn0); y = abs(ipro.mass['N2'])-abs(ipro.mass['T2']); hmax = 100
            if yis == 430:
                htit = 'm(#tilde{N}_{3}) - m(#tilde{#tau}_{2})  [GeV]'; fnn = "%s_N3mT2"%(fnn0); y = abs(ipro.mass['N3'])-abs(ipro.mass['T2']); hmax = 100
            if yis == 431:
                htit = 'm(#tilde{N}_{4}) - m(#tilde{#tau}_{2})  [GeV]'; fnn = "%s_N4mT2"%(fnn0); y = abs(ipro.mass['N4'])-abs(ipro.mass['T2']); hmax = 100
            if yis == 432:
                htit = 'm(#tilde{C}_{2}) - m(#tilde{#tau}_{2})  [GeV]'; fnn = "%s_C2mT2"%(fnn0); y = abs(ipro.mass['C2'])-abs(ipro.mass['T2']); hmax = 100

            if yis == 433:
                htit = 'm(#tilde{C}_{1}) - m(#tilde{#nu}_{l})  [GeV]'; fnn = "%s_C1mvl"%(fnn0); y = abs(ipro.mass['C1'])-abs(ipro.mass['ve']); hmax = 100
            if yis == 434:
                htit = 'm(#tilde{N}_{2}) - m(#tilde{#nu}_{l})  [GeV]'; fnn = "%s_N2mvl"%(fnn0); y = abs(ipro.mass['N2'])-abs(ipro.mass['ve']); hmax = 100
            if yis == 435:
                htit = 'm(#tilde{N}_{3}) - m(#tilde{#nu}_{l})  [GeV]'; fnn = "%s_N3mvl"%(fnn0); y = abs(ipro.mass['N3'])-abs(ipro.mass['ve']); hmax = 100
            if yis == 436:
                htit = 'm(#tilde{N}_{4}) - m(#tilde{#nu}_{l})  [GeV]'; fnn = "%s_N4mvl"%(fnn0); y = abs(ipro.mass['N4'])-abs(ipro.mass['ve']); hmax = 100
            if yis == 437:
                htit = 'm(#tilde{C}_{2}) - m(#tilde{#nu}_{l})  [GeV]'; fnn = "%s_C2mvl"%(fnn0); y = abs(ipro.mass['C2'])-abs(ipro.mass['ve']); hmax = 100

            if yis == 438:
                htit = 'm(#tilde{C}_{1}) - m(#tilde{#nu}_{l})  [GeV]'; fnn = "%s_C1mvl"%(fnn0); y = abs(ipro.mass['C1'])-abs(ipro.mass['ve']); hmax = 100
            if yis == 439:
                htit = 'm(#tilde{N}_{2}) - m(#tilde{#nu}_{l})  [GeV]'; fnn = "%s_N2mvl"%(fnn0); y = abs(ipro.mass['N2'])-abs(ipro.mass['ve']); hmax = 100
            if yis == 440:
                htit = 'm(#tilde{N}_{3}) - m(#tilde{#nu}_{l})  [GeV]'; fnn = "%s_N3mvl"%(fnn0); y = abs(ipro.mass['N3'])-abs(ipro.mass['ve']); hmax = 100
            if yis == 441:
                htit = 'm(#tilde{N}_{4}) - m(#tilde{#nu}_{l})  [GeV]'; fnn = "%s_N4mvl"%(fnn0); y = abs(ipro.mass['N4'])-abs(ipro.mass['ve']); hmax = 100
            if yis == 442:
                htit = 'm(#tilde{C}_{2}) - m(#tilde{#nu}_{l})  [GeV]'; fnn = "%s_C2mvl"%(fnn0); y = abs(ipro.mass['C2'])-abs(ipro.mass['ve']); hmax = 100

            if yis == 443:
                htit = 'm(#tilde{C}_{1}) - m(#tilde{#nu}_{#tau})  [GeV]'; fnn = "%s_C1mvT"%(fnn0); y = abs(ipro.mass['C1'])-abs(ipro.mass['vT']); hmax = 100
            if yis == 444:
                htit = 'm(#tilde{N}_{2}) - m(#tilde{#nu}_{#tau})  [GeV]'; fnn = "%s_N2mvT"%(fnn0); y = abs(ipro.mass['N2'])-abs(ipro.mass['vT']); hmax = 100
            if yis == 445:
                htit = 'm(#tilde{N}_{3}) - m(#tilde{#nu}_{#tau})  [GeV]'; fnn = "%s_N3mvT"%(fnn0); y = abs(ipro.mass['N3'])-abs(ipro.mass['vT']); hmax = 100
            if yis == 446:
                htit = 'm(#tilde{N}_{4}) - m(#tilde{#nu}_{#tau})  [GeV]'; fnn = "%s_N4mvT"%(fnn0); y = abs(ipro.mass['N4'])-abs(ipro.mass['vT']); hmax = 100
            if yis == 447:
                htit = 'm(#tilde{C}_{2}) - m(#tilde{#nu}_{#tau})  [GeV]'; fnn = "%s_C2mvT"%(fnn0); y = abs(ipro.mass['C2'])-abs(ipro.mass['vT']); hmax = 100


            if yis == 448:
                htit = 'm(#tilde{C}_{1}) - m(#tilde{l}_{min})  [GeV]'; fnn = "%s_C1ml1"%(fnn0); y = abs(ipro.mass['C1'])-msle_min; hmax = 100  # lightest slep = l1
            if yis == 449:
                htit = 'm(#tilde{N}_{2}) - m(#tilde{l}_{min})  [GeV]'; fnn = "%s_N2ml1"%(fnn0); y = abs(ipro.mass['N2'])-msle_min; hmax = 100 
            if yis == 450:
                htit = 'm(#tilde{N}_{3}) - m(#tilde{l}_{min})  [GeV]'; fnn = "%s_N3ml1"%(fnn0); y = abs(ipro.mass['N3'])-msle_min; hmax = 100 
            if yis == 451:
                htit = 'm(#tilde{N}_{4}) - m(#tilde{l}_{min})  [GeV]'; fnn = "%s_N4ml1"%(fnn0); y = abs(ipro.mass['N4'])-msle_min; hmax = 100 
            if yis == 452:
                htit = 'm(#tilde{C}_{2}) - m(#tilde{l}_{min})  [GeV]'; fnn = "%s_C2ml1"%(fnn0); y = abs(ipro.mass['C2'])-msle_min; hmax = 100 


            if yis == 453:
                htit = 'm(#tilde{C}_{1}) - m(#tilde{#nu}_{min})  [GeV]'; fnn = "%s_C1mv1"%(fnn0); y = abs(ipro.mass['C1'])-msnu_min; hmax = 100  # lightest snu  = v1
            if yis == 454:
                htit = 'm(#tilde{N}_{2}) - m(#tilde{#nu}_{min})  [GeV]'; fnn = "%s_N2mv1"%(fnn0); y = abs(ipro.mass['N2'])-msnu_min; hmax = 100 
            if yis == 455:
                htit = 'm(#tilde{N}_{3}) - m(#tilde{#nu}_{min})  [GeV]'; fnn = "%s_N3mv1"%(fnn0); y = abs(ipro.mass['N3'])-msnu_min; hmax = 100 
            if yis == 456:
                htit = 'm(#tilde{N}_{4}) - m(#tilde{#nu}_{min})  [GeV]'; fnn = "%s_N4mv1"%(fnn0); y = abs(ipro.mass['N4'])-msnu_min; hmax = 100 
            if yis == 457:
                htit = 'm(#tilde{C}_{2}) - m(#tilde{#nu}_{min})  [GeV]'; fnn = "%s_C2mv1"%(fnn0); y = abs(ipro.mass['C2'])-msnu_min; hmax = 100 


            if yis == 458:
                htit = 'm(#tilde{l}_{R}) - m(#tilde{N}_{1})  [GeV]'; fnn = "%s_lRmN1"%(fnn0); y = abs(ipro.mass['eR'])-abs(ipro.mass['N1']); hmax = 100 
            if yis == 459:
                htit = 'm(#tilde{l}_{L}) - m(#tilde{N}_{1})  [GeV]'; fnn = "%s_lLmN1"%(fnn0); y = abs(ipro.mass['eL'])-abs(ipro.mass['N1']); hmax = 100 
            if yis == 460:
                htit = 'm(#tilde{#tau}_{1}) - m(#tilde{N}_{1})  [GeV]'; fnn = "%s_T1mN1"%(fnn0); y = abs(ipro.mass['T1'])-abs(ipro.mass['N1']); hmax = 100 
            if yis == 461:
                htit = 'm(#tilde{#nu}_{l}) - m(#tilde{N}_{1})  [GeV]'; fnn = "%s_vlmN1"%(fnn0); y = abs(ipro.mass['ve'])-abs(ipro.mass['N1']); hmax = 100 
            if yis == 462:
                htit = 'm(#tilde{#nu}_{T}) - m(#tilde{N}_{1})  [GeV]'; fnn = "%s_vTmN1"%(fnn0); y = abs(ipro.mass['vT'])-abs(ipro.mass['N1']); hmax = 100 

            if yis == 463:
                htit = 'm(#tilde{l}_{1}) - m(#tilde{N}_{1})  [GeV]'; fnn = "%s_l1mN1"%(fnn0); y = msle_min-abs(ipro.mass['N1']); hmax = 100 
            if yis == 464:
                htit = 'm(#tilde{v}_{1}) - m(#tilde{N}_{1})  [GeV]'; fnn = "%s_v1mN1"%(fnn0); y = msnu_min-abs(ipro.mass['N1']); hmax = 100 
            if yis == 465:
                htit = 'min( m(#tilde{l}_{1}, m(#tilde{#nu}_{1})) ) - m(#tilde{N}_{1})  [GeV]'; fnn = "%s_L1mN1"%(fnn0); y = msnu_min-abs(ipro.mass['N1']); hmax = 100  # L1 = min(v1,l1)



            if yis == 466:
                htit = 'm(#tilde{l}) - m(#tilde{#nu}_{l})  [GeV]'; fnn = "%s_SLmSN"%(fnn0); y = abs(ipro.mass['eL'])-abs(ipro.mass['ve']); hmax = 100
            if yis == 467:
                htit = 'm(#tilde{#tau}_{1}) - m(#tilde{#nu}_{#tau})  [GeV]'; fnn = "%s_T1mvT"%(fnn0); y = abs(ipro.mass['T1'])-abs(ipro.mass['vT']); hmax = 100
            if yis == 468:
                htit = 'm(#tilde{#tau}_{2}) - m(#tilde{#nu}_{#tau})  [GeV]'; fnn = "%s_T2mvT"%(fnn0); y = abs(ipro.mass['T2'])-abs(ipro.mass['vT']); hmax = 100
            if yis == 469:
                htit = 'm(#tilde{#tau}_{2}) - m(#tilde{#tau}_{1})  [GeV]'; fnn = "%s_T2mT1"%(fnn0); y = abs(ipro.mass['T2'])-abs(ipro.mass['T1']); hmax = 100
            

            # =============== BRs


            
        # -------------------------------------------------- 201 - 299 (or 10201 - 10299)  BRs
        if 201 <= yis <= 299 or 10201 <= yis <= 10299: 

            # DECAY: LIGHTEST HIGGS
            if yis == 201: 
                htit = "BR: h -> #gamma#gamma  [#times 10^{-4}]"; fnn = "%s_BR_h_yy" %(fnn0);  y = ph.BR('y','y')*1e4; hMax = 100
                
            if yis == 202: 
                htit = "BR: h -> bb  [%]"; fnn = "%s_BR_h_bb" %(fnn0);  y = ph.BR('b','b')*1e2; hMax = 100

            if yis == 203: 
                htit = "BR: h -> gg  [%]"; fnn = "%s_BR_h_gg" %(fnn0);  y = ph.BR('g','g')*1e2; hMax = 100

            if yis == 204: 
                htit = "BR: h -> qq  [%]"; fnn = "%s_BR_h_qq" %(fnn0);  y = (ph.BR('u','u') + ph.BR('d','d'))*1e2; hMax = 100

            if yis == 205: 
                htit = "BR: h -> #tau#tau  [%]"; fnn = "%s_BR_h_tautau" %(fnn0);  y = ph.BR('ta','ta')*1e2; hMax = 100

            if yis == 206: 
                htit = "BR: h -> ZZ*  [%]"; fnn = "%s_BR_h_Zz" %(fnn0);  y = ph.BR('Z')*1e2; hMax = 100

            if yis == 207: 
                htit = "BR: h -> WW*  [%]"; fnn = "%s_BR_h_Ww" %(fnn0);  y = ph.BR('W')*1e2; hMax = 100

            if yis == 208: 
                htit = "BR: h -> WW  [%]"; fnn = "%s_BR_h_WW" %(fnn0);  y = ph.BR('W','W')*1e2; hMax = 100
                
            #if yis == 209: 
            #    htit = "BR: h -> ZZ  [%]"; fnn = "%s_BR_h_ZZ" %(fnn0);  y = ph.BR('Z','Z')*1e2; hMax = 100
                

            if yis == 210:
                htit = "BR: h -> N1N1  [%]"; fnn = "%s_BR_h_N1N1" %(fnn0);  y = ph.BR('N1','N1')*1e2; hMax = 100

            if yis == 211: 
                htit = "BR: h -> N1N2  [%]"; fnn = "%s_BR_h_N1N2" %(fnn0);  y = ph.BR('N1','N2')*1e2; hMax = 100

            if yis == 212: 
                htit = "BR: h -> C1C1  [%]"; fnn = "%s_BR_h_C1C1" %(fnn0);  y = ph.BR('C1','C1')*1e2; hMax = 100

            if yis == 213: 
                htit = "BR: h -> C1C2  [%]"; fnn = "%s_BR_h_C1C2" %(fnn0);  y = ph.BR('C1','C2')*1e2; hMax = 100

            BRhSM =  ph.BR('b','b') +  ph.BR('ta','ta') + ph.BR('u','u') + ph.BR('d','d') +  ph.BR('y','y') + ph.BR('g','g') + ph.BR('Z') + ph.BR('W')
            BRhSUSY = ph.BR('N1','N1') + ph.BR('N1','N2') + ph.BR('C1','C1')

            if yis == 217: 
                htit = "BR: h -> SUSY  [%]"; fnn = "%s_BR_h_SUSY" %(fnn0);  y = BRhSUSY*1e2; hMax = 100

            if yis == 218: 
                htit = "BR: h -> SM  [%]"; fnn = "%s_BR_h_SM" %(fnn0);  y = BRhSM*1e2; hMax = 100

            if yis == 219: 
                htit = "BR: h -> SMSUSY  [%]"; fnn = "%s_BR_h_SMSUSY" %(fnn0);  y = (BRhSUSY+BRhSM)*1e2; hMax = 100


            
            # DECAY: N2
            if yis == 220:
                htit = "BR: N2 -> h"; fnn = "%s_BR_N2_h" %(fnn0);  y = pN2.BR('N1','h')
            if yis == 221:
                htit = "BR: N2 -> Z"; fnn = "%s_BR_N2_Z" %(fnn0);  y = pN2.BR('N1','Z')
            if yis == 222: 
                htit = "BR: N2 -> qq"; fnn = "%s_BR_N2_qq" %(fnn0);  y = pN2.BR('N1','u') + pN2.BR('N1','d') + pN2.BR('N1','b')

            if yis == 223: 
                htit = "BR: N2 -> C1W"; fnn = "%s_BR_N2_C1W" %(fnn0);  y = pN2.BR('C1','W')
            if yis == 224: 
                htit = "BR: N2 -> C1lv"; fnn = "%s_BR_N2_C1lv" %(fnn0);  y = pN2.BR('C1','nue') + pN2.BR('C1','nut')
            if yis == 225: 
                htit = "BR: N2 -> C1qq"; fnn = "%s_BR_N2_C1qq" %(fnn0);  y = pN2.BR('C1','u') + pN2.BR('C1','d') + pN2.BR('C1','b')
            if yis == 226: 
                htit = "BR: N2 -> C1Ws"; fnn = "%s_BR_N2_C1Ws" %(fnn0);  y = pN2.BR('C1','nue') + pN2.BR('C1','nut') +  pN2.BR('C1','u') + pN2.BR('C1','d') + pN2.BR('C1','b')
            if yis == 227: 
                htit = "BR: N2 -> N1Zs"; fnn = "%s_BR_N2_N1Zs" %(fnn0);  y = pN2.BR('N1','e') + pN2.BR('N1','ta') + pN2.BR('N1','nue') + pN2.BR('N1','nut') +  pN2.BR('N1','u') + pN2.BR('N1','d') + pN2.BR('N1','b')


            if yis == 230:
                htit = "BR: N2 -> ll"; fnn = "%s_BR_N2_ll" %(fnn0);  y = pN2.BR('eR') + pN2.BR('eL') + pN2.BR('N1','e') + pN2.BR('T1') + pN2.BR('T2') + pN2.BR('N1','ta'); 
            if yis == 231:
                htit = "BR: N2 -> ll (sl)"; fnn = "%s_BR_N2_ll_sl" %(fnn0);  y = pN2.BR('eR') + pN2.BR('eL') + pN2.BR('T1') + pN2.BR('T2')
            if yis == 232:
                htit = "BR: N2 -> ll (*)"; fnn = "%s_BR_N2_ll_s" %(fnn0);  y = pN2.BR('N1','e') + pN2.BR('N1','ta')


            if yis == 233:
                htit = "BR: N2 -> ee"; fnn = "%s_BR_N2_ee" %(fnn0);  y = pN2.BR('eR') + pN2.BR('eL') + pN2.BR('N1','e'); 
            if yis == 234:
                htit = "BR: N2 -> ee (sl)"; fnn = "%s_BR_N2_ee_sl" %(fnn0);  y = pN2.BR('eR') + pN2.BR('eL') 
            if yis == 235:
                htit = "BR: N2 -> ee (R)"; fnn = "%s_BR_N2_ee_R" %(fnn0);  y = pN2.BR('eR') 
            if yis == 236:
                htit = "BR: N2 -> ee (L)"; fnn = "%s_BR_N2_ee_L" %(fnn0);  y = pN2.BR('eL') 
            if yis == 237:
                htit = "BR: N2 -> ee (*)"; fnn = "%s_BR_N2_ee_s" %(fnn0);  y = pN2.BR('N1','e') 

            if yis == 238:
                htit = "BR: N2 -> TT"; fnn = "%s_BR_N2_TT" %(fnn0);  y = pN2.BR('T1') + pN2.BR('T2') + pN2.BR('N1','ta'); 
            if yis == 239:
                htit = "BR: N2 -> TT (sl)"; fnn = "%s_BR_N2_TT_sl" %(fnn0);  y = pN2.BR('T1') + pN2.BR('T2')
            if yis == 240:
                htit = "BR: N2 -> TT (1)"; fnn = "%s_BR_N2_TT_1" %(fnn0);  y = pN2.BR('T1')
            if yis == 241:
                htit = "BR: N2 -> TT (2)"; fnn = "%s_BR_N2_TT_2" %(fnn0);  y = pN2.BR('T2')
            if yis == 242:
                htit = "BR: N2 -> TT (*)"; fnn = "%s_BR_N2_TT_s" %(fnn0);  y = pN2.BR('N1','ta')

            if yis == 243:
                htit = "BR: N2 -> H/A"; fnn = "%s_BR_N2_HA" %(fnn0);  y = pN2.BR('N1','H') + pN2.BR('N1','A')
            if yis == 244:
                htit = "BR: N2 -> y"; fnn = "%s_BR_N2_y" %(fnn0);  y = pN2.BR('N1','y')
            if yis == 245:
                htit = "BR: N2 -> vv"; fnn = "%s_BR_N2_vv" %(fnn0);  y = pN2.BR('ve') + pN2.BR('vT') + pN2.BR('N1','nue') + pN2.BR('N1','nut')
            if yis == 246:
                htit = "BR: N2 -> vv (*)"; fnn = "%s_BR_N2_vv" %(fnn0);  y = pN2.BR('N1','nue') + pN2.BR('N1','nut')

            if yis == 247:
                htit = "BR: N2 -> vv (sn)"; fnn = "%s_BR_N2_vv_sn" %(fnn0);  y = pN2.BR('ve') + pN2.BR('vT')

            if yis == 248:
                htit = "BR: N2 -> vev (sn)"; fnn = "%s_BR_N2_vev_sn" %(fnn0);  y = pN2.BR('ve')

            if yis == 249:
                htit = "BR: N2 -> vev / eLe"; fnn = "%s_BR_N2_vevONeLe" %(fnn0);
                if pN2.BR('eL') != 0: y = pN2.BR('ve') / pN2.BR('eL')
                else: y = -1



            if yis == 10220:
                htit = "BR: N2 -> 2L"; fnn = "%s_BR_N2_2L" %(fnn0);  y = pN2.BR('eR|eL|T1|T2') + pN2.BR('N1','e|ta') + pN2.BR('N1','Z')*(mypdg.ZtoLL+mypdg.ZtoTT)

            if yis == 10221:
                htit = "BR: N2 -> 1L"; fnn = "%s_BR_N2_1L" %(fnn0);  y = pN2.BR('C1','nue|nut') + pN2.BR('C1','W')*(mypdg.WtoL+mypdg.WtoT)
 
            if yis == 10222:
                htit = "BR: N2 -> 0L"; fnn = "%s_BR_N2_0L" %(fnn0);  y = pN2.BR('C1','W')*mypdg.Wtoq + pN2.BR('C1','u|d|b') + pN2.BR('N1','nue|nut|u|d|b') + pN2.BR('N1','Z')*( 1-(mypdg.ZtoLL+mypdg.ZtoTT) ) + pN2.BR('N1','h')

            if yis == 10224: 
                htit = "BR: N2 -> C1ev"; fnn = "%s_BR_N2_C1ev" %(fnn0);  y = pN2.BR('C1','nue') 
            if yis == 10225: 
                htit = "BR: N2 -> C1Tv"; fnn = "%s_BR_N2_C1Tv" %(fnn0);  y = pN2.BR('C1','nut') 


            # DECAY: C1

            if yis == 251: 
                 htit = "BR: C1 -> W"; fnn = "%s_BR_C1_W" %(fnn0);  y = pC1.BR('N1','W'); fnn = "%s_BR_C1_W" %(fnn0) 

            if yis == 10252:
                htit = "BR: C1 -> N1Ws"; fnn = "%s_BR_C1_N1Ws" %(fnn0);  y = pC1.BR('N1','nue') + pC1.BR('N1','nut') +  pC1.BR('N1','u') + pC1.BR('N1','d') + pC1.BR('N1','b')


            if yis == 252: 
                htit = "BR: C1 -> qq"; fnn = "%s_BR_C1_qq" %(fnn0);  y = pC1.BR('N1','u') + pC1.BR('N1','d') + pC1.BR('N1','b')


            if yis == 253: 
                htit = "BR: C1 -> lv";         fnn = "%s_BR_C1_lv" %(fnn0);  y = pC1.BR('vT') + pC1.BR('ve') + pC1.BR('eR') + pC1.BR('eL') + pC1.BR('T1') + pC1.BR('T2') + pC1.BR('N1','nue') + pC1.BR('N1','nut')

            if yis == 254: 
                htit = "BR: C1 -> lv (sl)"; fnn = "%s_BR_C1_lv_sl" %(fnn0);  y = pC1.BR('vT') + pC1.BR('ve') + pC1.BR('eR') + pC1.BR('eL') + pC1.BR('T1') + pC1.BR('T2')

            if yis == 255: 
                htit = "BR: C1 -> lv (*)"; fnn = "%s_BR_C1_lv_s" %(fnn0);  y = pC1.BR('N1','nue') + pC1.BR('N1','nut')



            if yis == 256: 
                htit = "BR: C1 -> ev"; fnn = "%s_BR_C1_ev" %(fnn0);  y = pC1.BR('ve') + pC1.BR('eR') + pC1.BR('eL') + pC1.BR('N1','nue')
            if yis == 257: 
                htit = "BR: C1 -> Tv"; fnn = "%s_BR_C1_Tv" %(fnn0);  y = pC1.BR('vT') + pC1.BR('T1') + pC1.BR('T2') + pC1.BR('N1','nut')

            if yis == 258: 
                htit = "BR: C1 -> ev (sl)"; fnn = "%s_BR_C1_ev_sl" %(fnn0);  y = pC1.BR('ve') + pC1.BR('eR') + pC1.BR('eL')
            if yis == 259: 
                htit = "BR: C1 -> Tv (sl)"; fnn = "%s_BR_C1_Tv_sl" %(fnn0);  y = pC1.BR('vT') + pC1.BR('T1') + pC1.BR('T2')

            if yis == 260: 
                htit = "BR: C1 -> ev (*)"; fnn = "%s_BR_C1_ev_s" %(fnn0);  y = pC1.BR('N1','nue')
            if yis == 261: 
                htit = "BR: C1 -> Tv (*)"; fnn = "%s_BR_C1_Tv_s" %(fnn0);  y = pC1.BR('N1','nut')

            if yis == 268: 
                htit = "BR: C1 -> N1+PIs"; fnn = "%s_BR_C1_N1PIs" %(fnn0);  y = pC1.BR('N1','PI+-') + pC1.BR('N1','PI0')  # seems that PI+- starts

                
            if yis == 269: 
                htit = "BR: C1 -> e (incl decay of #tau, W)"; fnn = "%s_BR_C1_eTot" %(fnn0);  y = pC1.BR('ve') + pC1.BR('eR') + pC1.BR('eL') + pC1.BR('N1','nue') + ( pC1.BR('vT') + pC1.BR('T1') + pC1.BR('T2') + pC1.BR('N1','nut') ) * mypdg.TtoL + pC1.BR('N1','W') * mypdg.WtoLtot

                


            # DECAY: N3
            if yis == 270: 
                htit = "BR: N3 -> sl"; fnn = "%s_BR_N3_sl" %(fnn0);  y = pN3.BR('eR') + pN3.BR('eL') + pN3.BR('N1','e'); fnn = "%s_BR_N3_sl" %(fnn0)

            
            if yis == 284:
                htit = "BR: N3 -> ee (sl)"; fnn = "%s_BR_N3_ee_sl" %(fnn0);  y = pN3.BR('eR') + pN3.BR('eL') 
            if yis == 289:
                htit = "BR: N3 -> TT (sl)"; fnn = "%s_BR_N3_TT_sl" %(fnn0);  y = pN3.BR('T1') + pN3.BR('T2')
            if yis == 273: 
                htit = "BR: N3 -> C1W"; fnn = "%s_BR_N3_C1W" %(fnn0);  y = pN3.BR('C1','W')
            if yis == 276: 
                htit = "BR: N3 -> C1Ws"; fnn = "%s_BR_N3_C1Ws" %(fnn0);  y = pN3.BR('C1','nue') + pN3.BR('C1','nut') +  pN3.BR('C1','u') + pN3.BR('C1','d') + pN3.BR('C1','b')

            if yis == 271:
                htit = "BR: N3 -> Z"; fnn = "%s_BR_N3_Z" %(fnn0);  y = pN3.BR('N1','Z')
            if yis == 272: 
                htit = "BR: N3 -> qq"; fnn = "%s_BR_N3_qq" %(fnn0);  y = pN3.BR('N1','u') + pN3.BR('N1','d') + pN3.BR('N1','b')
                print('DDDEBUG  udb: %9.6f   q: %9.6f' %(pN3.BR('N1','u') + pN3.BR('N1','d') + pN3.BR('N1','b'), pN3.BR('N1','q')))   # So 'q' doesn't work
 
            # DECAY: stau
            if yis == 291:
                htit = "BR: T1 -> N1"; fnn = "%s_BR_T1_N1" %(fnn0);  y = pT1.BR('N1')
            if yis == 292:
                htit = "BR: T1 -> C1"; fnn = "%s_BR_T1_C1" %(fnn0);  y = pT1.BR('C1')
            

            if yis == 299:
                htit = 'freeBR'
                fnn = '%s_freeBR' %(fnn0)
                # print 'yis2: ',yis2
                if 'htit' in list(yis2.keys()): htit = yis2['htit']
                if 'fnnadd' in list(yis2.keys()): fnn = '%s_%s' %(fnn0, yis2['fnnadd'])
                elif 'fnn' in list(yis2.keys()): fnn = yis2['fnn']

                brs = yis2['BRs']  # is a list of lists, e.g brs = [ [-1.,'N4','N1','Z'], ['N4','N1','h'] ]
                #                     when a number is given first for a list, the BR will be scaled by this (may be useful feature)

                # print yis2['BRs'], '  AND  ', brs

                y = 0. 

                for ibrs in range(len(brs)):
                    br = list(brs[ibrs])
                    br_orig = brs[ibrs]
                    if len(br) < 1:
                        print('libIsaplot  yis=299  Warning len(br)=%i  skipping..' %(len(br)))
                        continue

                    # Any scale factor?
                    zsc = 1. # can scale BRs (e.g. subtract)
                    if ( type(br[0]) is float or type(br[0]) is int ): zsc = float(br.pop(0))

                    if len(br) == 0: # hack used to insert constants, e.g. 1 to be used in 1-BR(xx)-BR(yy).. as in 'the rest'
                        y += zsc
                        continue

                    # Then identify the decaying particle
                    pParentT = br.pop(0)
                    if pParentT in list(ipro.br.part.keys()): pParent = ipro.br.part[pParentT]
                    else:
                        print('libIsaplot  Warning: parent not found: %s  (skipping)' %(pParentT))
                        continue

                    # Then get the decay
                    #print 'DEBUG br: ',br
                    if len(br) == 0:
                        print('libIsaplot  Warning  br list not appropriate: ', br_orig)
                        continue
                        
                    zd1 = br.pop(0)
                    zd2 = ''
                    if br: zd2 = br.pop(0)
                    
                    zy = pParent.BR(zd1,zd2)

                    #print 'DEBUG: ',pParentT, zd1, zd2, zy

                    if br: print('libIsaplot:  Warning  there are additional SM? particles not considered', br_orig)
                    y += zy*zsc
                    
                del brs
                # print 'y = ', y
                

        # -------------------------------------------------- 301 - 399 : COMP (++)
        if 301 <= yis <= 399: 

            # ======================================= COMP

            if yis == 301:
                htit = "N1compB";  fnn = "%s_N1compB" %(fnn0);  y = ipro.wig['n11']**2
            if yis == 302:
                htit = "N1compW";  fnn = "%s_N1compW" %(fnn0);  y = ipro.wig['n12']**2
            if yis == 303:
                htit = "N1compH1"; fnn = "%s_N1compH1" %(fnn0); y = ipro.wig['n13']**2
            if yis == 304:
                htit = "N1compH2"; fnn = "%s_N1compH2" %(fnn0); y = ipro.wig['n14']**2

            if yis == 306:
                htit = "N2compB";  fnn = "%s_N2compB" %(fnn0);  y = ipro.wig['n21']**2
            if yis == 307:
                htit = "N2compW";  fnn = "%s_N2compW" %(fnn0);  y = ipro.wig['n22']**2
            if yis == 308:
                htit = "N2compH1"; fnn = "%s_N2compH1" %(fnn0); y = ipro.wig['n23']**2
            if yis == 309:
                htit = "N2compH2"; fnn = "%s_N2compH2" %(fnn0); y = ipro.wig['n24']**2

            if yis == 311:
                htit = "N3compB";  fnn = "%s_N3compB" %(fnn0);  y = ipro.wig['n31']**2
            if yis == 312:
                htit = "N3compW";  fnn = "%s_N3compW" %(fnn0);  y = ipro.wig['n32']**2
            if yis == 313:
                htit = "N3compH1"; fnn = "%s_N3compH1" %(fnn0); y = ipro.wig['n33']**2
            if yis == 314:
                htit = "N3compH2"; fnn = "%s_N3compH2" %(fnn0); y = ipro.wig['n34']**2

            if yis == 316:
                htit = "N4compB";  fnn = "%s_N4compB" %(fnn0);  y = ipro.wig['n41']**2
            if yis == 317:
                htit = "N4compW";  fnn = "%s_N4compW" %(fnn0);  y = ipro.wig['n42']**2
            if yis == 318:
                htit = "N4compH1"; fnn = "%s_N4compH1" %(fnn0); y = ipro.wig['n43']**2
            if yis == 319:
                htit = "N4compH2"; fnn = "%s_N4compH2" %(fnn0); y = ipro.wig['n44']**2


            #if yis == 321:
            #    htit = "C1compW"; fnn = "%s_N4compH1" %(fnn0); y = ipro.wig['n43']**2

            # #### overall xsecs, filtereffs etc. (based on dict)



        # -------------------------------------------------- 501 - 599
        if 501 <= yis <= 599: 
            if yis in range(500,530) and not dict2: 
                print('WARNING  libIsaplot::get2DhistM1M2MU  For yis = %i need to provide dict2' %(yis))
                return {'status':'continue'}

            if yis == 501: 
                htit = "xsecLO"; y = safeGet(dict2,(M1,M2,mu)); fnn = '%s' %(fnn0)
                

            if yis == 502: 

                if (M1,M2,mu,yis2) in list(dict2.keys()): 
                    y = dict2[M1,M2,mu,yis2][yis3]; 
                    if y == -1: return {'status':'continue'}  # is -1 when there are no generated events
                    if y == 0: y = 1e-7  # 0 is not plotted with 'text' but 1e-7 is

                    if yis3 == 'eff': 
                        Nbeg = dict2[M1,M2,mu,yis2]['Nbeg']
                        Nend = dict2[M1,M2,mu,yis2]['Nend']
                        if Nbeg>0: yerr = sqrt(Nend)/Nbeg
                        else: yerr = -1
                        

                    fnn = '%s_%s_%s' %(fnn0, gg['subT'], yis3)
                    htit = 'Subprocess %s: %s' %(gg['subT'], yis3)
                else: 
                    return {'status':'continue'}


 
            if yis == 503: 
                # htit = "xsecLO"; 
                if (M1,M2,mu,yis2) in list(dict2.keys()): 
                    denom = dict2[M1,M2,mu, 999 ][yis3]
                    if denom == 0: y = 0.
                    else: y = 1. * dict2[M1,M2,mu,yis2][yis3] / denom
                else: return {'status':'continue'}
 


            if yis == 514: # DGemtKat_2ifb: nSR1mc
                zkey = (int(M1),int(M2),int(mu),yis2)
                if zkey in list(dict2.keys()):
                    keyvar = 'n_mc_SR1'
                    if keyvar in dict2[zkey]: 
                        y = dict2[zkey][keyvar]
                    else: 
                        print('WARNING  libIsaplot::get2DhistM1M2MU  get2D: yis: %i  zkey: %s  has no keyvar %s' %(yis, zkey, keyvar))

                else:
                    print('WARNING  libIsaplot::get2DhistM1M2MU  get2D: yis = %i   Not found: %s' %(yis, str(zkey)))
                    y = 0.
                fnn = '%s_nSRmc' %(fnn0)
                htit = 'nSR1mc (RAW) (kat)'


            if yis == 515: # DGemtKat_2ifb: nFiltmc
                zkey = (int(M1),int(M2),int(mu),yis2)
                if zkey in list(dict2.keys()):
                    keyvar = 'n_mc_Filt'
                    if keyvar in dict2[zkey]: 
                        y = dict2[zkey][keyvar]
                    else: 
                        print('WARNING  libIsaplot::get2DhistM1M2MU  get2D: yis: %i  zkey: %s  has no keyvar %s' %(yis, zkey, keyvar))

                else:
                    print('WARNING  libIsaplot::get2DhistM1M2MU  Not found: ',zkey)
                    y = 0.
                fnn = '%s_nFiltmc' %(fnn0)
                htit = 'nFiltmc (RAW) (kat)'




            if yis == 516: # DGemtKat_2ifb: nFiltmc_nonzeroSR1
                zkey = (int(M1),int(M2),int(mu),yis2)
                if zkey in list(dict2.keys()):
                    keyvar = 'n_mc_Filt_nonzeroSR1'
                    if keyvar in dict2[zkey]: 
                        y = dict2[zkey][keyvar]
                    else: 
                        print('WARNING  libIsaplot::get2DhistM1M2MU  get2D: yis: %i  zkey: %s  has no keyvar %s' %(yis, zkey, keyvar))

                else:
                    print('Not found: ',zkey)
                    y = 0.
                fnn = '%s_nFiltmc_nonzeroSR1' %(fnn0)
                htit = 'nFiltmc_nonzeroSR1 (RAW) (kat)'



            if yis == 517: # DGemtKat_2ifb: accFilt2SR1_mc
                y = 0.
                zkey = (int(M1),int(M2),int(mu),yis2)
                if zkey in list(dict2.keys()):
                    keynumer = 'n_mc_SR1'
                    keydenom = 'n_mc_Filt'
                    if keydenom in dict2[zkey]: 
                        denom = dict2[zkey][keydenom]  # should never happen
                        if denom == 0: print('get2D: yis: %i  zkey: %s  has zero denom (nFilt)' %(yis, zkey))
                        if denom > 0 and keynumer in dict2[zkey]: 
                            y = 1. * dict2[zkey][keynumer] / denom
                    else: 
                        print('get2D: yis: %i  zkey: %s  has no keyvar %s' %(yis, zkey, keydenom))

                else:
                    print('Not found: ',zkey)
                    y = 0.
                fnn = '%s_nFiltmc' %(fnn0)
                htit = 'acc (nSR1_mc/nFilt_mc) (RAW) (kat)'
                


            if yis == 510: # DGemtKat_2ifb: SR1(lum)
                zkey = (int(M1),int(M2),int(mu),yis2)
                if zkey in list(dict2.keys()):
                    keyvar = 'n_lum_SR1'
                    if keyvar in dict2[zkey]: y = dict2[zkey][keyvar]

                    keyerr = 'ERRn_lum_SR1'
                    if keyerr in dict2[zkey]:
                        yerr = dict2[zkey][keyerr]
                        #print 'setting yerr: %.2f' %(yerr)
                else:
                    print('Not found: ',zkey)
                    y = 0.
                fnn = '%s_2ifb' %(fnn0)
                htit = 'SR1 (DGemt_Katarina)'


            if yis == 512: # DGemtKat_2ifb
                zkey = (int(M1),int(M2),int(mu),yis2)
                if not zkey in list(dict2.keys()):
                    print('WINFO: not in dict2: ', zkey)
                    return {'status':'continue'}
                if not zkey in list(dict3.keys()):
                    print('WINFO: not in dict3: ', zkey)
                    return {'status':'continue'}
                if dict3[zkey] == 0:
                    print('WINFO: denom (dict3) is zero for key', zkey)
                    return {'status':'continue'}

                # Now all should be fine
                keyvar = 'n_lum_SR1'
                keyerr = 'ERRn_lum_SR1'
                yfull = dict2[zkey][keyvar]
                yfast = dict3[zkey][keyvar]
                y = 1. * yfull / yfast
                errfull = dict2[zkey][keyerr]
                errfast = dict3[zkey][keyerr]
                syst = 0.

                fnn = '%s_2ifb' %(fnn0)                
                htit = 'SR1 (DGemt_Katarina)'
                
                if 'relsysterr' in OptDict:
                    syst = OptDict['relsysterr'] * y
                    if OptDict['relsysterr'] > 0: 
                        fnn += '_relsyst%.2f' %(OptDict['relsysterr'])
                        htit += ' RelSystErr: %.2f' %(OptDict['relsysterr'])
                yerr = sqrt( 1. * errfast**2 * yfull**2 / yfast**4 + 1. * errfull**2 / yfast**2 + syst**2) 
                


            if yis == 511: # total acceptance = DGemtKat_2ifb /  theory nev >=3L at 2ifb  [to check if ratio has some reliability]
                zkey = (int(M1),int(M2),int(mu),yis2)
                if zkey in list(dict2.keys()):
                    keyvar = 'n_lum_SR1'
                    if keyvar in dict2[zkey]: y = dict2[zkey][keyvar]
                else:
                    print('libIsaplot yis= %i   Not found: %s' %(yis,zkey))
                    y = 0.

                den = (ipro.play.xAccu[3] * 2000.) 
                if den != 0: y /= den  # <-- Here divide with the theoretical 2000 ipb min value (for 3+ leptons): giving a "total acceptance"
                else:
                    print(' %s  hmm, denom=0, setting to -1: yis = %i' %(M1M2mu, yis))
                    y = -1
                fnn = '%s_totalacceptance3L' %(fnn0)
                htit = 'Total Acceptance SR1 (DGemt_Katarina)'



            if yis == 513: # DGemtKat_2ifb_polated =  total_acceptance * theory nev >=3L at 2ifb  [to interpolate/extrapolate to new points]
                y = 0.
                
                yTruth = (ipro.play.xAccu[3] * 2000.) 
                if 'hmap' in dict2: 
                    yTotAcc = dict2['hmap'].GetBinContent(binx,biny)
                    y = yTotAcc * yTruth
                    

                fnn = '%s_SR_from_totalacc_and_truth' %(fnn0)
                htit = 'SR1 modelled from Tot Acc and Truth'


            if yis == 518: # polated:  nFilt = desired_nSR1 / accFilt2SR1_mc_polate
                y = 0.
                nSR1 = 0.
                if 'hmap' in dict2 and 'nSR1' in dict2: 
                    denom = dict2['hmap'].GetBinContent(binx,biny)
                    nSR1 = dict2['nSR1']
                    y = nSR1 / denom
                else: 
                    print("Need 'nSR1' and 'hmap' in dict2")
                    
                fnn = '%s_Expected_nFilt_with_nSR1_eq_%.0f' %(fnn0, nSR1)
                htit = 'Expected nFilt to get nSR1=%.0f' %(nSR1)


            if yis == 519: # polated:  nSR1 = nFilt(input) * accFilt2SR1_mc_polate
                y = 0.
                nFilt = 0.
                if 'hmap' in dict2 and 'nFilt' in dict2: 
                    numer = dict2['hmap'].GetBinContent(binx,biny)
                    nFilt = dict2['nFilt']
                    y = nFilt * numer
                    
                fnn = '%s_Expected_nSR1_with_nFilt_eq_%.0f' %(fnn0, nFilt)
                htit = 'Expected nSR1 to get nFilt=%.0f' %(nFilt)


            if yis == 520: # polated:  nSR1 = nFilt(input) * accFilt2SR1_mc_polate
                y = 0.
                if 'extpri' in dict2: 
                    if keyM1M2mu in list(dict2['extpri'].keys()): 
                        thispri = dict2['extpri'][keyM1M2mu]
                        if 'prieq' in OptDict: y = thispri == OptDict['prieq']
                        else: y = thispri
                    else: 
                        #print keyM1M2mu; print dict2['extpri'].keys()
                        pass

                else: 
                    #print 'opt:', opt
                    #print 'dict2.keys():', dict2.keys()
                    pass
                if 'prieq' in OptDict: prieqT = str(OptDict['prieq'])
                else: prieqT = ''
                fnn = '%s_pri%s' %(fnn0, prieqT)
                htit = 'Extended points Pri %s' %(prieqT)
                

            if yis == 521: # 
                if not IkeyM1M2mu in dict2['phys2id']: 
                    if VB: print('yis=%i  phys2id misses Ikey %s' %(yis,IkeyM1M2mu))
                    return {'status':'continue'}
                if not IkeyM1M2mu in dict2['extpri']: 
                    if VB: print('yis=%i  extpri misses Ikey %s' %(yis,IkeyM1M2mu))
                    return {'status':'continue'}
                out = '%s  %s  %2i  %5i' %(M1M2mu, dict2['phys2id'][IkeyM1M2mu], dict2['extpri'][IkeyM1M2mu], dict2['hnev'].GetBinContent(binx,biny))
                print(out)
                outs.append(out)



            if yis == 531:  # effective cross-section (need to remove some by hand)
                #print 'asdfasdf'
                htit = "Filter Cross-section [pb]"
                fnn = "%s_effectivexsec" %(fnn0)
                zkey = (intM1,intM2,intmu,yis2)
                #print zkey
                if zkey not in list(dict2.keys()):
                    if VB: print('yis=%i  effective xsec missing for %s' %(yis, zkey))
                    #print dict2.keys()
                    return {'status':'continue'}
                y = dict2[zkey]


            #print yis, yis2
            if yis == 551:  # tina
                
                htit = "DGemtCut: %s" %(yis2)
                fnn = "%s_%s" %(fnn0, yis2)
                y = 0.
                subprocs = [999]
                if 'subprocs' in list(dict3.keys()): subprocs = list(dict3['subprocs'])
                for sub in subprocs: 
                    zkey = (intM1,intM2,intmu, sub)
                    if zkey in list(dict2.keys()): 
                        y += dict2[zkey][yis2]

                
            if yis == 552:  # tina: ratio not from C1N2
                
                #htit = "DGemtCut: %s" %(yis2)
                #htit = "N2C1/XX(tot)  at %s" %(yis2)
                htit = "#tilde{#chi}^{0}_{2}#tilde{#chi}^{#pm}_{1} / #tilde{#chi}#tilde{#chi}  at  %s" %(yis2)
                fnn = "%s_%s_N2C1onTOT" %(fnn0, yis2)
                y = 0.
                if 'subprocs' in list(dict3.keys()): subprocs = list(dict3['subprocs'])
                zkey = (intM1,intM2,intmu, 999)
                if not zkey in list(dict2.keys()):
                    pass
                elif dict2[zkey][yis2] == 0:
                    y = -1.
                else:
                    den = dict2[zkey][yis2]
                    num = 0.
                    for sub in subprocs: 
                        zkey = (intM1,intM2,intmu, sub)
                        if zkey in list(dict2.keys()): 
                            num += dict2[zkey][yis2]
                    y = num / den


            if yis == 561:  # darkSUSY  # yis2 e.g. 'tot', if 1:excluded, if 0:non-excluded
                htit = "Exclusion (%s)" %(yis2)
                fnn = "%s_excl%s" %(fnn0, yis2)
                if keyM1M2mu in list(dict2.keys()):
                    y = dict2[keyM1M2mu][yis2]
                    y = fabs(y) 


            if yis == 562:  # darkSUSY  ; looking only at key excl, but with yis2 eq. to 'C1', 'SL', 'N1' or 'Zwidth', which are the relevant exclusions for DGemt
                zdum = {'C1':'chargino mass', 'SL':'slepton mass', 'N1':'LSP mass', 'Zwidth':'Z width'}
                htit = "Exclusion (%s)" %(zdum[yis2])
                fnn = "%s_excl%s" %(fnn0, yis2)
                if keyM1M2mu in list(dict2.keys()):
                    if yis2 in dict2[keyM1M2mu]['excl']: y = 1  # excluded by yis2 (e.g. C1 mass)
                    else: y = 0



            if yis == 571:  #eirik's 2L-cutflows
                htit = "DGemt2L: %s" %(yis3)
                fnn = "%s_%s" %(fnn0, yis2[0])  #FRAGILE
                y = 0. 
                if keyM1M2mu in list(dict2.keys()): 
                    y = 0.
                    for yis2b in yis2: #use yis2 as an array
                        if yis2b in list(dict2[keyM1M2mu].keys()): 
                            y += dict2[keyM1M2mu][yis2b]
                            y += 1e-7  #(to plot zeros)
                        else: 
                            print('Warning: for set %s  key %s is not found' %(str(keyM1M2mu), yis2b)) 


            if yis in [572,573,574,575,576]:  #eirik's 2L-cutflows
                htit = "DGemt2L:  Z_{LLR} of %s" %(yis2)
                fnn = "%s_ZLLR_%s" %(fnn0, yis2)
                y = 0. 
                if keyM1M2mu in list(dict2.keys()): 
                    y = 0.
                    yis2b = 's' + yis2
                    if yis2b in list(dict2[keyM1M2mu].keys()): 
                        s = dict2[keyM1M2mu][yis2b]
                        if s == 0: 
                            s = 1e-5
                        b = dict3[yis3][yis2]
                        delb = dict3['err'][yis2]
                        if yis == 572: y = Z_LLR(s,b)
                        if yis == 573: y = Z_LLR(s,b+delb)
                        if yis == 574: y = Z_LLR(s,b+delb*2)

                        if yis == 575: y = s / sqrt(b)
                        if yis == 576: y = s / sqrt(b+delb*delb)

                    else: 
                        print('Warning: for set %s  key %s is not found' %(str(keyM1M2mu), yis2)) 



        # -------------------------------------------------- 601 - 699  DARKSUSY
        if 601 <= yis <= 699: 

            

            if IsInRangesT('601-631',yis) and not ipro.ds:
                return {'status':'continue'}
            
            #    hMin = 0.
            #    hMax = 1.

            if yis == 601:
                # print yis2
                # print ipro.ds
                codeT = ''
                for z in yis2: codeT += z
                htit = "DS excluded: %s" %(codeT)
                fnn = "%s_ds_excl_%s" %(fnn0, codeT)
                bits = yis2  # format bits = ['Z','C',..] for the things to test
                y = Darksusy_these_bits_notexcluded(ipro.ds, bits)

            if yis == 602:
                # print 'keys: ', ipro.ds.keys()
                # print 'ipro.ds: ', ipro.ds
                # print fnn0
                
                dskey = 'ds_excl_combined'
                htit = "DS excluded: combined" 
                fnn = "%s_%s" %(fnn0, dskey)
                y = ipro.ds[dskey]  # NB: 0 is ok, 1 is not ok

            if yis == 603:
                dskey = 'ds_excl_combined_except_h'
                htit = "DS excluded except h"
                fnn = "%s_%s" %(fnn0, dskey)
                y = ipro.ds[dskey]  # NB: 0 is ok, 1 is not ok
                
            if yis == 604:  # <--- kind of a standard in combination with C1 curve
                dskey = 'ds_excl_combined_except_Ch'
                htit = "DS excluded except h,C1"
                fnn = "%s_%s" %(fnn0, dskey)
                y = ipro.ds[dskey]  
                
            if yis == 611:
                dskey = 'ds_excl_h'; htit = "DS excluded: h"; fnn = "%s_%s" %(fnn0, dskey); y = ipro.ds[dskey] 
            if yis == 612:
                dskey = 'ds_excl_C'; htit = "DS excluded: C"; fnn = "%s_%s" %(fnn0, dskey); y = ipro.ds[dskey] 
            if yis == 613:
                dskey = 'ds_excl_Z'; htit = "DS excluded: Z"; fnn = "%s_%s" %(fnn0, dskey); y = ipro.ds[dskey] 
            if yis == 614:
                dskey = 'ds_excl_L'; htit = "DS excluded: sl"; fnn = "%s_%s" %(fnn0, dskey); y = ipro.ds[dskey] 
            if yis == 615:
                dskey = 'ds_excl_N'; htit = "DS excluded: N"; fnn = "%s_%s" %(fnn0, dskey); y = ipro.ds[dskey] 
            if yis == 616:
                dskey = 'ds_excl_bsgam'; htit = "DS excluded: bsgam"; fnn = "%s_%s" %(fnn0, dskey); y = ipro.ds[dskey] 
            if yis == 617:
                dskey = 'ds_excl_rho'; htit = "DS excluded: rho"; fnn = "%s_%s" %(fnn0, dskey); y = ipro.ds[dskey] 
            if yis == 618:
                dskey = 'ds_excl_G'; htit = "DS excluded: gl"; fnn = "%s_%s" %(fnn0, dskey); y = ipro.ds[dskey] 
            if yis == 619:
                dskey = 'ds_excl_Q'; htit = "DS excluded: sq"; fnn = "%s_%s" %(fnn0, dskey); y = ipro.ds[dskey] 

            if yis == 621:
                dskey = 'ds_gminus2'; htit = "DS gminus2"; fnn = "%s_%s" %(fnn0, dskey); y = ipro.ds[dskey]
                
            if yis == 622:
                dskey = 'ds_cdm'; htit = "DS cdm"; fnn = "%s_%s" %(fnn0, dskey); y = ipro.ds[dskey]
                
            if yis == 623:
                dskey = 'ds_cdm2'; htit = "DS cdm (detailed)"; fnn = "%s_%s" %(fnn0, dskey); y = ipro.ds[dskey]
                
            if yis == 624:
                htit = "DS cdm (best)"; fnn = "%s_ds_cdmbest" %(fnn0);
                if ipro.ds['ds_cdm2'] >= 0: y = ipro.ds['ds_cdm2']
                else: y = ipro.ds['ds_cdm']
                

            # should maybe add a combination of 622 and 623 for cases when 623 is not calculated
            # should also add limit cuts (one and two-sided at 1, 2 and 3 sigma)
            # cdm = 0.1123 with error 0.0035 (http://en.wikipedia.org/wiki/Lambda-CDM_model)
            # so varying a few sigmas will have practically no effect
            # (the theoretical error, on the DS-calculation is probably higher ; can I get to it?)

            if 625 <= yis <= 631:  # <-- FRAGILE (remember to updated 631 if more added)
                cdm_susy = ipro.ds['ds_cdm2']
                if cdm_susy < 0: y = -9  # use this to avoid plotting if non-calculated

                if cdm_susy < 0: cdm_susy = ipro.ds['ds_cdm']

                if yis == 625:
                    y = int(cdm_susy > mypdg.cdm + 3*mypdg.cdm_error)
                    htit = "DS cdm(susy) > cdm(wmap) [3#sigma]"
                    fnn = "%s_ds_cdm_gt_wmap_3sig" %(fnn0);
                    
                if yis == 626:
                    y = int(cdm_susy > mypdg.cdm + 2*mypdg.cdm_error)
                    htit = "DS cdm(susy) > cdm(wmap) [2#sigma]"
                    fnn = "%s_ds_cdm_gt_wmap_2sig" %(fnn0);
                    
                if yis == 627:
                    y = int(cdm_susy > mypdg.cdm + 1*mypdg.cdm_error)
                    htit = "DS cdm(susy) > cdm(wmap) [1#sigma]"
                    fnn = "%s_ds_cdm_gt_wmap_1sig" %(fnn0);
                    
                if yis == 628: # within 3 sigma (both ways)
                    y = int( abs(cdm_susy - mypdg.cdm) < 3*mypdg.cdm_error )
                    htit = "DS | cdm(susy) - cdm(wmap) | < 3#sigma"
                    fnn = "%s_ds_cdm_within_3sig" %(fnn0);

                if yis == 629: # within factor 2 (both ways)
                    y = int( mypdg.cdm/2 < cdm_susy < 2*mypdg.cdm ) 
                    htit = "DS cdm:  wmap/2 < susy < wmap*2" 
                    fnn = "%s_ds_cdm_within_factor2" %(fnn0);

                if yis == 630: # within factor 2 (one ways)
                    y = int( 0 < cdm_susy < 2*mypdg.cdm ) 
                    htit = "DS cdm(susy) < 2*cdm(wmap)" 
                    fnn = "%s_ds_cdm_below_factor2" %(fnn0);

                if yis == 631: # take ratio - can have contours at e.g. 0.5, 1., 1.5, .. whatever
                    y = cdm_susy / mypdg.cdm  
                    htit = "DS cdm(susy) / cdm(wmap)" 
                    fnn = "%s_ds_cdmratio_susyVSwmap" %(fnn0);

            if yis == 632:
                SM_gminus2 = 0.00116592080  # http://en.wikipedia.org/wiki/Muon
                dskey = 'ds_gminus2'; htit = "DS gminus2 / SM value"; fnn = "%s_%s_over_SM" %(fnn0, dskey); y = ipro.ds[dskey] / SM_gminus2

            
            if ( (601 <= yis <= 619) or (625 <= yis <= 630) ) and y == 0: y = 1e-7   # Hack to plot the zeros ... 


        # -------------------------------------------------- 901 - dictionaries
        if 901 <= yis <= 999:

            if yis == 901:
                htit = "Generator Filter Efficiency (GFE) [%]"
                fnn = "%s_geneff" %(fnn0)
                try:
                    y = OptDict['gDict']['subproceff'][keyID]['eff'] * 100. 
                except:
                    y = 0
 
            if yis == 902:
                htit = "Generated 7 TeV: nev [kev]"
                fnn = "%s_nevgen7TeV" %(fnn0)
                try:
                    y = OptDict['gDict']['nevgentot'][keyID] / 1000.   # number in kev
                except:
                    y = 0

            if yis == 904:
                htit = "Request 8 TeV: nev [kev]"
                fnn = "%s_nevgen8TeV" %(fnn0)
                try:
                    y = OptDict['gDict']['nevgen8TeV'][keyID] / 1000.   # number in kev
                except:
                    y = 0

            if yis == 905:
                htit = "Request 8 TeV: point priority"
                fnn = "%s_pointpri8TeV" %(fnn0)
                try:
                    y = OptDict['gDict']['pri8TeV'][keyID]
                except:
                    y = 0


            if 910 <= yis <= 919:  # tja
                try: zdict = OptDict['gDict']['request']
                except: return {'status':'continue'}
                
            if yis == 910:
                ztt='nFmin'; htit = "Request: %s [kev]" %(ztt); fnn = "%s_%s" %(fnn0,ztt)
                try: y = zdict[keyID][ztt]/1000
                except: return {'status':'continue'}
            if yis == 911:
                ztt='nFmax'; htit = "Request: %s [kev]" %(ztt); fnn = "%s_%s" %(fnn0,ztt)
                try: y = zdict[keyID][ztt]/1000
                except: return {'status':'continue'}
            if yis == 912:
                ztt='nFmaxCut'; htit = "Request: %s [kev]" %(ztt); fnn = "%s_%s" %(fnn0,ztt)
                try: y = zdict[keyID][ztt]/1000
                except: return {'status':'continue'}

            if yis == 920:
                zSR='3LSR1a'
                ztt='sig'; htit = "%s: %s" %(zSR,ztt); fnn = "%s_%s_%s" %(fnn0,zSR,ztt)
                #print OptDict['gDict']['limit'].keys()
                try: y = OptDict['gDict']['limit'][zSR][keyID][ztt]
                except: return {'status':'continue'}
            if yis == 921:
                zSR='3LSR1a'
                ztt='sigerr_stat_vs_sig'; htit = "%s: %s [%%]" %(zSR,ztt); fnn = "%s_%s_%s" %(fnn0,zSR,ztt)
                #print OptDict['gDict']['limit'].keys()
                try: y = OptDict['gDict']['limit'][zSR][keyID][ztt] * 100. 
                except: return {'status':'continue'}
            if yis == 922:
                zSR='3LSR1a'
                ztt='sigerr_stat_vs_toterr'; htit = "%s: %s [%%]" %(zSR,ztt); fnn = "%s_%s_%s" %(fnn0,zSR,ztt)
                #print OptDict['gDict']['limit'].keys()
                try: y = OptDict['gDict']['limit'][zSR][keyID][ztt] * 100.
                except: return {'status':'continue'}

            if yis == 923:
                zSR='2LSR4'
                ztt='sig'; htit = "%s: %s" %(zSR,ztt); fnn = "%s_%s_%s" %(fnn0,zSR,ztt)
                #print OptDict['gDict']['limit'].keys()
                try: y = OptDict['gDict']['limit'][zSR][keyID][ztt]
                except: return {'status':'continue'}
            if yis == 924:
                zSR='2LSR4'
                ztt='sigerr_stat_vs_sig'; htit = "%s: %s [%%]" %(zSR,ztt); fnn = "%s_%s_%s" %(fnn0,zSR,ztt)
                #print OptDict['gDict']['limit'].keys()
                try: y = OptDict['gDict']['limit'][zSR][keyID][ztt] * 100.
                except: return {'status':'continue'}
            if yis == 925:
                zSR='2LSR4'
                ztt='sigerr_stat_vs_toterr'; htit = "%s: %s [%%]" %(zSR,ztt); fnn = "%s_%s_%s" %(fnn0,zSR,ztt)
                #print OptDict['gDict']['limit'].keys()
                try: y = OptDict['gDict']['limit'][zSR][keyID][ztt] * 100.
                except: return {'status':'continue'}


            # raw & scaled numbers from Tina, June17 2012, for 2011 4.7 fb-1 ; in preps of 2011-extension
            if yis == 930:
                zSR='3LSR1a'
                zt1='raw'; zt2='nev'; htit = "%s: %s %s" %(zSR,zt1,zt2); fnn = "%s_%s_%s%s" %(fnn0,zSR,zt1,zt2)
                #print OptDict['gDict']['limit'].keys()
                try: y = OptDict['gDict']['nev3LSR'][zt1][zSR][keyID][999] # sum
                except: return {'status':'continue'}
            if yis == 931:
                zSR='3LSR1b'
                zt1='raw'; zt2='nev'; htit = "%s: %s %s" %(zSR,zt1,zt2); fnn = "%s_%s_%s%s" %(fnn0,zSR,zt1,zt2)
                #print OptDict['gDict']['limit'].keys()
                try: y = OptDict['gDict']['nev3LSR'][zt1][zSR][keyID][999] # sum
                except: return {'status':'continue'}
            if yis == 932:
                zSR='3LSR2' 
                zt1='raw'; zt2='nev'; htit = "%s: %s %s" %(zSR,zt1,zt2); fnn = "%s_%s_%s%s" %(fnn0,zSR,zt1,zt2)
                #print OptDict['gDict']['limit'].keys()
                try: y = OptDict['gDict']['nev3LSR'][zt1][zSR][keyID][999] # sum
                except: return {'status':'continue'}
            if yis == 933:
                zSR='3LSR1a'
                zt1='scaled'; zt2='nev'; htit = "%s: %s %s" %(zSR,zt1,zt2); fnn = "%s_%s_%s%s" %(fnn0,zSR,zt1,zt2)
                #print OptDict['gDict']['limit'].keys()
                try: y = OptDict['gDict']['nev3LSR'][zt1][zSR][keyID][999] # sum
                except: return {'status':'continue'}
            if yis == 934:
                zSR='3LSR1b'
                zt1='scaled'; zt2='nev'; htit = "%s: %s %s" %(zSR,zt1,zt2); fnn = "%s_%s_%s%s" %(fnn0,zSR,zt1,zt2)
                #print OptDict['gDict']['limit'].keys()
                try: y = OptDict['gDict']['nev3LSR'][zt1][zSR][keyID][999] # sum
                except: return {'status':'continue'}
            if yis == 935:
                zSR='3LSR2'
                zt1='scaled'; zt2='nev'; htit = "%s: %s %s" %(zSR,zt1,zt2); fnn = "%s_%s_%s%s" %(fnn0,zSR,zt1,zt2)
                #print OptDict['gDict']['limit'].keys()
                try: y = OptDict['gDict']['nev3LSR'][zt1][zSR][keyID][999] # sum
                except: return {'status':'continue'}




            if yis == 940:
                htit = "LEPlimit on C1 (gaugino)"
                fnn = "%s_LEPlimit_C1_gaugino" %(fnn0)
                delM = abs(ipro.mass['C1']) - abs(ipro.mass['N1'])
                if 'LEPlimit_C1_gaugino' in OptDict:
                    LEPlimit_C1 = OptDict['LEPlimit_C1_gaugino'].Eval(delM)
                    if abs(ipro.mass['C1']) < LEPlimit_C1: y = 1  # excluded (returns 1)
                    else: y = 0  # non-excluded
                else:
                    sys.exit("FATAL  libIsaplot  Need to have OptDict['LEPlimit_C1_gaugino'] defined to use yis=940")



        # -------------------------------------------------- 9001 - adhoc / tests
        if 9001 <= yis <= 9099:
            
            if yis == 9001:
                y = ( 100.-abs(ipro.mass['N1']) ) / (abs(ipro.mass['N2'])-abs(ipro.mass['N1']) ) * 100
                htit = "Lep ins. (% betw #tilde{N}_{1} and #tilde{N}_{2}) to get m(#tilde{l})=100 GeV" 
                fnn = "%s_LepAtXtogetmslep100" %(fnn0);

            if yis == 9002:
                y = ( 90.-abs(ipro.mass['N1']) ) / (abs(ipro.mass['N2'])-abs(ipro.mass['N1']) ) * 100
                htit = "Lep ins. (% betw #tilde{N}_{1} and #tilde{N}_{2}) to get m(#tilde{l})=90 GeV" 
                fnn = "%s_LepAtXtogetmslep90" %(fnn0);

            if yis == 9005:
                y = ( abs(ipro.mass['C1'])-abs(ipro.mass['N1']) ) / (abs(ipro.mass['N2'])-abs(ipro.mass['N1']) ) * 100
                htit = "m(#tilde{C}_{1}  (in % betw #tilde{N}_{1} and #tilde{N}_{2})" 
                fnn = "%s_C1posbetweenN1N2" %(fnn0);


            #if yis == 9010:
            #    y = 0
            #    if 'coverage' in OptDict:
            #        if 'truth' in OptDict['coverage']:
            #            if keyM1M2mu in OptDict['coverage']['truth']: y = 1
                
                

        # --------------------------------------------------
                    
        # -------------
        else: # can here implement more general structures (started by yis not being int)
            pass
        
        #===============
    
    res['y'] = y
    res['htit'] = htit
    res['fnn'] = fnn
    res['hMin'] = hMin
    res['hMax'] = hMax
    
    return res


# =================================
def get2DhistM1M2MU(gg, Scan, yis, VB=0, OPTs=[], ret='h', texs=[], sc=1.,  dict2={}, dict3={}, yis2=0, yis3='', flip=0, palette='coolhot1c', drawtext={}, grid={}, dopt='colz', opt={}, hCont=[], exclude=[], filetype='pdf', stdgrid={}, moveoverflow=0, debug=0, iScan=0, dict_free2d={}): 

    agridID = opt.get('fnbase','')

    free2d_whatx = ''
    if 'free2d_whatx' in opt:
        free2d_whatx = int(opt['free2d_whatx'])
        #print 'yesss'
    #print 'noooo'
    #print opt.keys()
    #print free2d_whatx, type(free2d_whatx)
        
    # whatx is a hack to allow plotting ordinary histos
    xfree2d = array('d')
    yfree2d = array('d')

    STDFORM = 1  # NB, this one is meant to be more implemented

    # standardise (will later rename OPTs and opt in argument list)
    OptList = OPTs
    OptDict = opt


    if 'drawtext_text' in OptDict: drawtext['text'] = OptDict['drawtext_text']  # 2012-06-17
    if 'drawtext_size' in OptDict: drawtext['size'] = OptDict['drawtext_size']
    if 'drawtext_form' in OptDict: drawtext['form'] = OptDict['drawtext_form']


    txtThis = 'libIsaplot::get2DhistM1M2MU  '
    txtINFO = 'INFO  '+txtThis
    txtWARNING = 'WARNING  '+txtThis
    txtHACK = 'HACK  '+txtThis

    
    if 'figdir' in OptDict:
        #print txtINFO + "from OptDict setting figdir to '%s'" %OptDict['figdir']
        figdir = OptDict['figdir']
    elif 'figdir' in gg:
        print(txtINFO + "from gg setting figdir to '%s'" %gg['figdir'])
        figdir = gg['figdir']
    else:
        figdir = 'fig1'


    if 0: 
        OptDict['fntag'] = '_TB10'
        print("HACK  libIsaplot::get2DhistM1M2MU ---------------> NB  HARDCODING   OptDict['fntag'] = %s  inside" %(OptDict['fntag']))


    aScan = Scan[iScan] 

    RATIO = 0
    if 'RATIO' in gg and gg['RATIO'] == 1:
        if len(Scan) <= 1:
            print("WARNING  libIsaplot::get2DhistM1M2MU:  'RATIO' is on, but len(Scan) == %i   ... PLOTTING NORMAL" %(len(Scan)))
        else: 
            RATIO = 1
    
    
    MOVEOVERFLOW = 0
    if 'useasymptotic' in OptList: MOVEOVERFLOW = 1
    if 'asymptotic' in OptDict: MOVEOVERFLOW = int(OptDict['asymptotic'])

    if 1:
        MOVEOVERFLOW = 1
        print('NB: MANUAL ASYMPTOTIC HACK: setting MOVEOVERFLOW = %i' %(MOVEOVERFLOW))

    ticklength = 0.003
    if 'ticklength' in gg: ticklength = gg['ticklength']


    if palette:
        if not ('palettestatus' in list(gg.keys()) and gg['palettestatus'] == 'set'):  # 2012-06-07: found that set_palette slows down. Hack. 
            set_palette(palette) # MEM this one causes memory(?) problems (gets slower and slower) 
            gg['palettestatus'] = 'set'

    #if 'palette:coolhot1' in OptList: set_palette('coolhot1')

    DEBUG = debug

    zmeta = {}

    POSx1x2 = 1
    #if 'POSx1x2' in OPTs.keys(): POSx1x2 = OPTs['POSx1x2']

    #linlog = 'float'
    #if 'linlog' in gg: linlog = gg['linlog']

    # ##################################################################
    # INIT (from input)

    fixvar = gg['fixvar']
    fixval = gg['fixval']
    # print 'DEBUG: OptDict', OptDict
    figname_class = ''
    if 'figname_class' in OptDict:
        figname_class = OptDict['figname_class']
        if figname_class.startswith('_'): figname_class = figname_class[1:]  # (overkill?)
        if figname_class.endswith('_'): figname_class = figname_class[:-1]
        # print 'hit ',figname_class
    

    #setup = GetAxisRangesEtc(sliceID=gg['scanID'], fixvar=fixvar, fixval=fixval, figname_class=figname_class)   # 2012-06-09
    setup = GetAxisRangesEtcNew(gg=gg, figname_class=figname_class)   # 2012-06-09
    
    val = setup['val']
    binsL = setup['binsL']
    #for key in binsL.keys(): print 'setup %-5s  %s' %(key, binsL[key])
    v1 = setup['v1']
    v2 = setup['v2']
    V1 = setup['V1']
    V2 = setup['V2']
    v1T = setup['v1T']
    v2T = setup['v2T']
    fnn0 = setup['fnn0']
    #print 'HERE fnn0: %s' %(fnn0)



    if 'binwidth' in list(gg.keys()):    # 2013-07-14 : IT ALSO WORKS FINE (PRODUCES PLOTS) WITH THE STANDARD PROCEDURE
        #print 'BINWIDTH is ', gg['binwidth']
        if 1:  # HACKHACK
            zhacked = 0
            if len(val['v1']) < len(val['v2']): val['v1'] = list(val['v2']) ; zhacked = 1
            if len(val['v2']) < len(val['v1']): val['v2'] = list(val['v1']) ; zhacked = 1
            if len(val[V1]) < len(val[V2]): val[V1] = list(val[V2]) ; zhacked = 1
            if len(val[V2]) < len(val[V1]): val[V2] = list(val[V1]) ; zhacked = 1
            
            #if len(val['v1']) < len(val['v2']): val['v1'] = list(val['v2'])
            #if len(val['v2']) < len(val['v1']): val['v2'] = list(val['v1'])
            if zhacked:
                print('NBNBNB: MANUAL HACK TO SYMMETRISE THE SLICE (FOR USE WHEN FIXVAR IS M2 OR MU) : val = ', val)
                
        binsL = GetNoncontBins(val, gg['binwidth'])   # in kilelib_ROOT.py 
        if VB>0: print("INFO  libIsaplot::get2DhistM1M2MU making noncont bins (this sets binsL anew ...)") 


    if 'fnbase' in OptDict and OptDict['fnbase']:
        #fnn0 = "%sMUM2_M1eq%s" %(OptDict['fnbase'], str(int(fixval)).zfill(3))  #HACK 2012-04-17
        fnn0 = "%s%seq%s" %(OptDict['fnbase'], GetVarText(fixvar,'varTxt'), str(int(fixval)).zfill(3))  #HACK 2013-07-12  # removed MUM2_ to be on same footing as colratio
        # NB: filename of simple is set here
        if figname_class: fnn0 += '_' + figname_class

    # ------------------------------------
    #if flip == 1:   # default as of 2010-04-08 is that flip=0 corresponds to (x,y)=(MU,M2)
    #    v1old = v1; v1=v2; v2=v1old
    #    V1old = V1; V1=V2; V2=V1old
    #    v1Told = v1T; v1T=v2T; v2T=v1Told
    
    
    # print 'DEBUG: ', v1, v2, v1T, v2T

    bins = {}
    nbins = {}
    for key in binsL:
        bins[key] = array('f')
        bins[key].fromlist(binsL[key])
        nbins[key] = len(binsL[key]) - 1


    
    hn = "%s_vs_%s_%s" %(V1,V2,str(random.randint(1e4,1e5-1)))
    ht = V1+'vs'+V2

    if 0:
        print('DEBUG')
        print(V1, V2)
        print(list(nbins.keys()))
        print(list(bins.keys()))

    if not (V1 in nbins and V2 in nbins):
        print("ERROR::libIsaplot  Will crash. Have not defined nbins[] / bins[] ; probably due to non-recognised gg['scanID']: %s, or could be related to gg['binwidth'] which fires GetNoncontBins (fails e.g. for DGemtFine)" %(gg['scanID']))
        if 'binwidth' in gg: print('Here binwidth is %.3f' %(gg['binwidth']))
        else: print("There is no gg['binwidth']")
        print("bins  : ", bins)
        print("binsL : ", binsL)
        print("nbins : ", nbins)
        print("bins  : ", bins)

    print('owl: ', V1,V2, list(bins.keys()), list(nbins.keys()))

    zh = ROOT.TH2F(hn,hn,nbins[V1],bins[V1], nbins[V2],bins[V2])
    #return zh  MEM already slow

    tax = zh.GetXaxis()
    tay = zh.GetYaxis()
    taz = zh.GetZaxis()

    # print 'DEBUG1: v1T:%s  v2T:%s' %(v1T,v2T)   # axis titles for simple (not ColourRatio2D) are set here
    tax.SetTitle(v1T)
    tay.SetTitle(v2T)
    tax.CenterTitle()
    tay.CenterTitle()
    tax.SetTitleOffset(OptDict.get('xoff',1.15))
    tay.SetTitleOffset(OptDict.get('yoff',1.4))
    # print 'LabelSize: ', tax.GetLabelSize()
    zlabelsize = 0.045
    tax.SetLabelSize(zlabelsize)  # 2013-07-12
    tay.SetLabelSize(zlabelsize)
    if 'xoff' in OptDict: print("DEBBUG  A;LKJ;FLKAJ;SLKDJFA;LSF  xoff set to %.3f" %(OptDict['xoff']))
    #print 'asdfasdfasf;alkj;ljs;dflasdf'

    tax.SetTickLength(ticklength)
    tay.SetTickLength(ticklength)
    

    ymaxmax = 0.
    
    fnn = 'dummy'
    htit = 'htit_dummy'

    yerr = -999.

    outs = []
    
    denom_titleremove = [' [pb]', ' [GeV]']   # to remove from title when 'denom'-option is set
        
    isRelativeQuantity = 0
    if IsInRangesT('10-15,20-24,50-58,59,61-64,73,201-299,10201-10299,301-319',yis): isRelativeQuantity = 1

    # ##################################################################
    # LOOP & FILLING
    iS = 0
    for iSS in range(aScan.N()): 
        iS += 1
        ipro = aScan.L[iSS]
        #if iSS > 10: continue  #DEBUG


        #print ipro.mass, v1, v2
        ### Choose x1 and x2
        if v1 in ipro.isaout: x1 = ipro.isaout[v1]  #<-- VARY
        elif v1 in ipro.mass: x1 = ipro.mass[v1]

        if v2 in ipro.isaout: x2 = ipro.isaout[v2]  #<-- VARY
        elif v2 in ipro.mass: x2 = ipro.mass[v2]

        # 
        if DEBUG: print('DEBUG  libIsaplot::get2DhistM1M2MU  %4.0f  %4.0f ' %(x1,x2))

        if POSx1x2: 
            x1 = abs(x1)
            x2 = abs(x2)

        binx = tax.FindBin(x1)
        biny = tay.FindBin(x2)



        # ---------------------------------------- SPECIFIC FOR (M1,M2,MU)-SCENARIO 
        # Choose M1, M2, mu from file
        M1 = abs(ipro.isaout['M_1'])
        M2 = abs(ipro.isaout['M_2'])
        mu = abs(ipro.isaout['MU'])
        
        intM1 = int(M1)
        intM2 = int(M2)
        intmu = int(mu)

        #print intM1,intM2,intmu 

        #M1M2mu = "%3s_%3s_%3s" %(str(M1).zfill(3), str(M2).zfill(3), str(mu).zfill(3))
        M1M2mu = "%3s_%3s_%3s" %(str(intM1).zfill(3), str(intM2).zfill(3), str(intmu).zfill(3))
        #keyM1M2mu = (M1,M2,mu)
        IkeyM1M2mu = (int(M1),int(M2),int(mu))
        keyM1M2mu = IkeyM1M2mu
        # Dec: suddenly(?) had to insert int (above) 

        if IkeyM1M2mu in exclude:
            print("INFO: libIsaplot::get2DhistM1M2MU  Deliberately excluding key %s (in exclude list)" %(str(IkeyM1M2mu)))
            continue
        # ----------------------------------------



        # ----- Testing N1 LSP
        lighterthanN1 = []
        for z9 in ['eR','eL','T1','T2','ve','vT','N2','N3','N4','C1','C2']:
            if abs(ipro.mass[z9]) < abs(ipro.mass['N1']): lighterthanN1.append([abs(ipro.mass[z9]), z9])
        if lighterthanN1:
            lighterthanN1.sort()
            out9 = 'WARNING::libIsaplot iSS:%i  %s   sparticles lighter than N1 (%.3f) : ' %(iSS, M1M2mu, abs(ipro.mass['N1']))
            for z9 in lighterthanN1: 
                out9 += '    %s (%.3f)' %(z9[1], z9[0])
            print(out9)
            # might here also skip the point
        # -----




        res = get_ipro_values(ipro, yis, fnn0, keyM1M2mu, OptDict, yis2, yis3, dict2, dict3)  # 2012-06-09
        # print 'HERE res: ', res

        #if 'warn' in res: print "%3s  %2i %2i  (%3i,%3i)  | %s" %(yis, iSS, iS, x1, x2, res['warn'])

        if res['status'] == 'continue': continue # semi-hack 

        # ----------
        # hack to plot ordinary 2d
        if free2d_whatx:
            resx = get_ipro_values(ipro, free2d_whatx, fnn0, keyM1M2mu, OptDict, yis2, yis3, dict2, dict3)  # 2012-07-29
            if resx['status'] == 'continue': continue

            xfree2d.append(resx['y'])
            yfree2d.append(res['y'])

            continue
        # ----------

        y = res['y']
        htit = res['htit']
        fnn = res['fnn']
        hMax = res['hMax']
        hMin = res['hMin']

        # possibility to prepend to htit
        htit = OptDict.get('htit_prepend', '') + htit + OptDict.get('htit_append','')
        if 'htit_remove' in OptDict:
            for rrr in OptDict['htit_remove']: htit = htit.replace(rrr,'')

        
        # ---------------------------------- Now allow for picking a denominator  # 2012-06-09
        if 'denom_yis' in OptDict:
            yisDenom = int(OptDict['denom_yis'])
            # ... here add options when needed (if need to change yis2, dict2, ...)
            # ... 
            res = get_ipro_values(ipro, yisDenom, fnn0, yis2, yis3, dict2, dict3)
            denom = res['y']

            # Re-set y-value
            if denom != 0:
                y /= denom
            else:
                print('WARNING: libIsaplot::denom  zero denominator, orig y = %.3f' %(y))
                y = -1.

            # well, htit and fnn is now set for every bin, but ok, no harm done
            # Re-set title: 
            htit = '%s / %s' %(htit, res['htit'])  # res['htit'] container the denominator
            
            for r in denom_titleremove: htit = htit.replace(r,'')
            if 'denom_htit' in OptDict: htit = OptDict['denom_htit']  # allows to simply set title by hand [NOT IN USE]
            #print 'htit: ', htit, ' ASDF ', res['htit']
            # Re-set filename
            #print 'DEBUG1: fnn: ',fnn,fnn0

            #fnn = fnn.replace(fnn0,'')  # 2012-07-12: removed
            fnn = "%s__RELATIVE_TO_%s" %(fnn, res['fnn'].replace(fnn0,''))
            #print 'DEBUG2: fnn: ',fnn
            if 'denom_fnn' in OptDict: fnn = OptDict['denom_fnn']  # allows to simply set title by hand
            #print 'DEBUG3: fnn: ',fnn

            isRelativeQuantity = 1

        # ----------------------------------


        # ---------------------------------- Now allow for scaling to lumi  # 2012-06-09
        if 'nev@ifb' in OptDict:
            lumi = float(OptDict['nev@ifb'])
            y *= lumi * 1000.   # xsec is in pb, lumi is in ifb so need to scale by additional 1000

            lumiT = '%.1f' %(lumi)
            if int(lumi) == lumi: lumiT = '%.0f' %(lumi)

            htit = htit.replace('Cross-section [pb]', 'NEV at %s fb^{-1}' %(lumiT))
            htit = htit.replace('#sigma', 'NEV at %s fb^{-1} ' %(lumiT))
            htit = htit.replace('[pb]','')
            htit = htit.replace('  ',' ')  # yes?

            fnn = fnn.replace('xsec','nevat%sifb_'%(lumiT))
            fnn = fnn.replace('__','_')  # yes? Fragile (Could change other parts of the filename)
            

        # ===================================================== y-values end
        # ===================================================== y-values end
        # ===================================================== y-values end




        # ===================================================== PLOTTING
        # ===================================================== PLOTTING
        # ===================================================== PLOTTING

        if ret in ['outs','none']: continue  # no plotting
        
        
        ##### CHECK RANGE ######### 
        if 1: 
            if binx > tax.GetNbins(): 
                if MOVEOVERFLOW:
                    binx = tax.GetNbins()
                    if DEBUG: print('WARNING: binx: %i > %i  (set to n (MOVEOVERFLOW==1))' %(binx, tax.GetNbins()))
                else:
                    if DEBUG: print('WARNING: binx: %i > %i  (SKIPPING, MOVEOVERFLOW==0)' %(binx, tax.GetNbins()))
                    continue
            if biny > tay.GetNbins():
                if MOVEOVERFLOW: 
                    if DEBUG: print('WARNING: biny: %i > %i  (set to n (MOVEOVERFLOW==1))' %(biny, tay.GetNbins()))
                    biny = tay.GetNbins()
                else: 
                    if DEBUG: print('WARNING: biny: %i > %i  (SKIPPING, MOVEOVERFLOW==0)' %(biny, tay.GetNbins()))
                    continue

        if DEBUG: print("%4i  %4i   %5.0f  %5.0f" %(binx, biny, x1, x2))


        ##### CHECKs #########
        # can here insert checks to drop point; 
        #   N1 not LSP
        #   bounds on C1, slepton, ...
        #   ...


        y *= sc   # here scale with input variable sc

        if 'percentage' in OptList and isRelativeQuantity: y *= 100.



        # CHECK DOUBLE-FILLING
        if 1:
            if zh.GetBinContent(binx,biny) != 0: print("WARNING: Trying to refill bin (%i,%i) [%i,%i]  old: %7.3f   new: %7.3f" %(binx,biny, x1,x2, zh.GetBinContent(binx,biny), y))

        zh.SetBinContent(binx,biny,y)

        if yerr != -999: zh.SetBinError(binx,biny,yerr)


        ymaxmax = max(ymaxmax, y)
        
        #print 'DEBUG v1:%-3s binx:%3i   v2:%-3s biny:%3i  y=%5.3f' %(v1,binx,v2,biny,y) 


    # ################################################################## Here continues/ends the looping over grid points
    # POST (nicify histo etc.) 

    # ---------- hack 2012-07-29
    if free2d_whatx:
        # define histo / tgraph
        # xfree2d_max = max(xfree2d)
        # xfree2d_min = min(xfree2d)
        # yfree2d_max = max(yfree2d)
        # yfree2d_min = min(yfree2d)
        # if not 'nx' in dict_free2d: dict_free2d['nx'] = 100
        # if not 'ny' in dict_free2d: dict_free2d['ny'] = 100

        gr = ROOT.TGraph(len(xfree2d), xfree2d, yfree2d) 

        return gr

        # return 


    # 


    if gg['scanID'] in ['DGemt','DGnoL','DGemt_n300','DGemt_n768', 'DG_n300','DG_n147'] and flip:   #  A WARNING!!
        for i in range(3): print(80*'+')
        print('   NBNBNB:  FLIP IS SET. THIS MEANS M2 IS PUT ON THE X-AXIS ... (OLD STANDARD)')
        for i in range(3): print(80*'+')


    if ret == 'outs': return outs  # later?
    if ret == 'none': return 


    if 'hMax' in list(gg.keys()): 
        hMax = gg['hMax']
        zh.SetMaximum(hMax)

    if 'hMin' in list(gg.keys()): 
        hMin = gg['hMin']
        zh.SetMinimum(hMin)

    if 'percentage' in OptList and isRelativeQuantity:
        zh.SetMinimum(0.)
        zh.SetMaximum(100.)
        htit += ' [%]'
        drawtext['form'] = '.0f'

    #if 'hMax' in OPTs and hMax != -998: 
    #    zh.SetMaximum(hMax)

    #if 'hMin' in OPTs and hMin != -999: 
    #    zh.SetMinimum(hMin)


    # Hack 2012-06-11
    if IsInRangesT('601-619,625-630',yis):
        hMin = 0.
        hMax = 1.
        zh.SetMinimum(0.)
        zh.SetMaximum(1.)



    # Hack 2012-10-10
    #print agridID
    if 1 and 'DGstauRfix95' in agridID:
        #print 'ASDFASDFASDF'
        htit = htit.replace('lepton','tau')
        htit = htit.replace('2L','2#tau')


    if VB>1: 
        print("Ex: Draw('colz')")

    if flip == 2: zzh = FlipAxes(zh)  # if ==1, do by hack above (easier and probably safer)
    else: zzh = zh

    # can here fill metainfo
    zmeta['fnn'] = fnn
    if 'fntag' in OptDict: zmeta['fnn'] += OptDict['fntag']
    if 'fntag2' in OptDict: zmeta['fnn'] += OptDict['fntag2']
    if 'fnpre' in OptDict: zmeta['fnn'] = OptDict['fnpre'] + zmeta['fnn']

    if '[pb]' in htit and sc == 1000: htit = htit.replace('[pb]','[fb]')
    elif sc > 10.: htit += '  [x %.0f]' %(sc)
    elif sc > 1.: htit += '  [x %.1f]' %(sc)
    zmeta['htit'] = htit

    themin = zh.GetMinimum()
    themax = zh.GetMaximum()

    #if linlog == 'float' and themin != 0 and themax/themin < 10.: 
    #    print 'INFO:het2DhistM1M2MU: themax/themin = %.3f  so using linear z-scale'
    #    c.SetLogy(0)

    if RATIO:
        ggTwo = dict(gg)
        ggTwo.pop('RATIO')
        hTwo = get2DhistM1M2MU(gg=ggTwo, Scan=Scan, yis=yis, VB=VB, OPTs=OptList, ret='h', texs=texs, sc=sc,  dict2=dict2, dict3=dict3, yis2=yis2, yis3=yis3, flip=flip, palette=palette, drawtext=drawtext, grid=grid, dopt=dopt, opt=OptDict, hCont=hCont, exclude=exclude, filetype=filetype, stdgrid=stdgrid, moveoverflow=moveoverflow, debug=debug, iScan=1)

        # hOne = h.Clone('hOne')
        # hOne.Divide(hTwo)

        if not gg['RATIO_invert']: 
            zzh.Divide(hTwo)

        else: 
            hTwo.Divide(zzh)
            #del zzh
            zzh = hTwo.Clone()

        if 'RATIO_text' in gg: drawtext['text'] = gg['RATIO_text']  #new
        if 'RATIO_form' in gg: drawtext['form'] = gg['RATIO_form']
        else: drawtext['form'] = '.2f'
        if 'RATIO_size' in gg: drawtext['size'] = gg['RATIO_size']
        else: drawtext['size'] = 1.0
        if 'RATIO_hMax' in gg: zzh.SetMaximum(gg['RATIO_hMax'])
        else: zzh.SetMaximum(2.)
        if 'RATIO_hMin' in gg: zzh.SetMinimum(gg['RATIO_hMin'])
        else: zzh.SetMinimum(0.)

        if 'RATIO_fnnTag' in gg: zmeta['fnn'] += gg['RATIO_fnnTag']
        else: zmeta['fnn'] += '_RATIO'

    #print '\nand now meta:    ', zmeta
    #print '\nand now optDict: ', OptDict

    if ret == 'h': return zzh
    if ret == 'meta': return zzh, zmeta
    if ret == 'c': return zzh, GetCanvas(zzh, meta=zmeta, texs=texs, drawtext=drawtext, grid=grid, dopt=dopt, hCont=hCont, filetype=filetype, stdgrid=stdgrid, figdir=figdir, yis=yis, OptDict=OptDict)
    if ret == 'hmc': return zzh, zmeta, GetCanvas(zzh, meta=zmeta, texs=texs, drawtext=drawtext, grid=grid, dopt=dopt, hCont=hCont, filetype=filetype, stdgrid=stdgrid, figdir=figdir, yis=yis, OptDict=OptDict)
    print('get2DhistM1M2MU: invalid ret: %s' %(ret))
    return zzh

#############
def GetCanvas(h, meta={}, dopt='colz', htit='', grid={'on':1}, drawtext={}, save=1, texs=[], hCont=[], filetype='pdf', stdgrid={}, epstopdf=True, figdir='fig1', yis=-1, OptDict={}): 

    STDFORM = 1  # NB, this one is meant to be more implemented
    setlog = OptDict.get('setlog',[])

    stamp = random.randint(1000,9999)
    cnam = 'c%i' %(stamp)
    c = ROOT.TCanvas(cnam,'',10,10,700,700)
    c.SetTopMargin(0.10)
    c.SetBottomMargin(0.15)
    c.SetLeftMargin(0.15)
    c.SetRightMargin(0.15)

    h.Draw(dopt)

    zgLine = ROOT.TLine()
    if grid: 
        gridcol = 1
        gridsty = 2
        gridwid = 2
        if 'col' in grid: gridcol = grid['col']
        if 'sty' in grid: gridsty = grid['sty']
        if 'wid' in grid: gridwid = grid['wid']
        DrawHistGrid(zgLine, h, col=gridcol, sty=gridsty, wid=gridwid)

    if hCont:
        for con in list(hCont.keys()): hCont[con].Draw('cont3same')

    if stdgrid:
        if 'x' in stdgrid: c.SetGrid(stdgrid['x'])
        if 'y' in stdgrid: c.SetGrid(stdgrid['y'])

    if drawtext: 
        if 'text' in drawtext: text = drawtext['text']
        else: text = 'text'
        text += ' same'

        if 'form' in drawtext: form = drawtext['form']
        else: form = '.2f'

        #if STDFORM and (yis in [401]): #2012-10-08 removed
        #    form = '.1f'

        gStyle.SetPaintTextFormat(form)


        if 'size' in drawtext: sz = drawtext['size']
        else: sz = 1.0
        h.SetMarkerSize(sz)

        if 'col' in drawtext: sz = drawtext['col']
        else: col = 1
        h.SetMarkerColor(col)

        h.Draw(text)

    if not htit and 'htit' in meta: htit = meta['htit']
    if htit: 
        txt = ROOT.TLatex()
        txt.SetNDC()
        txt.SetTextAlign(21)
        txt.SetTextSize(0.050)
        txt.DrawLatex(0.50,0.93,htit)

    for t in texs: t.Draw()


    if 'coverageList' in OptDict:
        DrawCoverage(List=OptDict['coverageList'],  opt={'col':2, 'sty':24, 'size':[2.3,2.4,2.5]}, M1=OptDict['fixval'])

    
    if 'x' in setlog: c.SetLogx(1)  # Hack
    if 'y' in setlog: c.SetLogy(1)  # Hack
    if 'z' in setlog: c.SetLogz(1)  # Hack

    if save:
       fnn = 'test'
       if 'fnn' in list(meta.keys()): fnn = meta['fnn']
       fn = '%s/%s.%s' %(figdir, fnn,filetype) # not in use now
       c.SaveAs('%s/%s.%s' %(figdir, fnn,'eps'))  # 2012-10-05: hack to get both pdf and eps
       os.system('epstopdf %s/%s.%s' %(figdir, fnn,'eps'))  # 2012-10-05: hack to get both pdf and eps
       #c.SaveAs(fn)
       #if filetype == 'eps' and epstopdf: 
       #    os.system('cd %s ; epstopdf %s.%s' %(figdir, fnn,filetype))

    
    return c

#############
def get2DhistSubXsec(gg, aScan, yis, VB=0, OPTs=[], ret='h', texs=[], sc=1.,  dict2={}, dict3={}, yis2=0, yis3='', flip=0, palette='coolhot1c', drawtext={}, grid={}, dopt='colz', opt={}, nbinx1=100, nbinx2=100, binx1m=0, binx1M=400, binx2m=0,binx2M=400): 
    '''
    gg = {}; gg['scanID']='xsecXX4all'; gg['fixvar']='M_1'; gg['fixval']=140
    # gg = {}; gg['scanID']='DGemt'; gg['fixvar']='M_1'; gg['fixval']=140

    zh = get2DhistM1M2MU(gg=gg,aScan=aScan,yis=105,VB=1)

    HMs:
      - fnn not in use

    '''
    OptList = OPTs
    OptDict = opt


    ticklength = 0.003
    if 'ticklength' in gg: ticklength = gg['ticklength']
    if ticklength != 0.03: print('INFO:get2DhistM1M2MU: ticklength: %.5f ' %(ticklength))

    # this one (fntag) is not really in use: update? 
    fntag = ''
    if 'fntag' in OptDict: fntag = OptDict['fntag']

    if palette:
        if not ('palettestatus' in list(gg.keys()) and gg['palettestatus'] == 'set'):  # 2012-06-07: found that set_palette slows down. Hack. 
            set_palette(palette)
            gg['palettestatus'] = 'set'
        
    #if 'palette:coolhot1' in OptList: set_palette('coolhot1')

    DEBUG = 0

    zmeta = {}

    POSx1x2 = 1
    #if 'POSx1x2' in OPTs.keys(): POSx1x2 = OPTs['POSx1x2']

    #linlog = 'float'
    #if 'linlog' in gg: linlog = gg['linlog']

    # ##################################################################
    # INIT (from input)

    fixvar = gg['fixvar']
    fixval = gg['fixval']
    val = {}
    binsL = {}


    #binx1m = 0
    #binx1M = 400
    #binx2m = 0
    #binx2M = 400

    #nbinx1 = 1000
    #nbinx2 = 1000

    binx1d = (binx1M-binx1m)/(1.*nbinx1)
    binx2d = (binx2M-binx2m)/(1.*nbinx2)

    binsx1 = array('f')
    binsx2 = array('f')

    # Uniform (can be more flexible)
    for i in range(0,nbinx1+1): binsx1.append(binx1m + i*binx1d)
    for i in range(0,nbinx2+1): binsx2.append(binx2m + i*binx2d)

    v1T = OPTs['v1T']  #<--- 2012-06-09: this will of course fail (is list): means that this routine is obsolete
    v2T = OPTs['v2T']
    strv1T = str(v1T)
    strv2T = str(v2T)
    hn = "xsec_%s_%s" %(strv1T,strv2T)
    zh = ROOT.TH2F(hn,hn, nbinx1,binsx1, nbinx2,binsx2)

    tax = zh.GetXaxis()
    tay = zh.GetYaxis()
    taz = zh.GetZaxis()

    # print 'DEBUG2: strv1T:%s  strv2T:%s' %(strv1T,strv2T)
    tax.SetTitle(strv1T)
    tay.SetTitle(strv2T)
    tax.CenterTitle()
    tay.CenterTitle()


    # ##################################################################
    # LOOP & FILLING
    nfilledbins = 0
    nduplicatebins = 0

    iS = 0
    for iSS in range(aScan.N()): 
        iS += 1
        ipro = aScan.L[iSS]

        # array (allows for several values per scenario to be filled, e.g. N masses)
        valsx1 = []
        valsx2 = []
        ys = []

        # VB: 
        VB1 = 'VB: iSS: %4i  iS: %4i ' %(iSS,iS)
        VBdet = "  %s  " %(ipro.ID)

        txtdet = VB1 + VBdet  # will be renamed

        # make array
        varsx1 = list(v1T)
        varsx2 = list(v2T)


        # Fill var array (xsecs)
        for varx1 in varsx1: 
            valx1 = abs(ipro.mass[varx1])
            if VB>=2: print(VB1+'varx1: %s  %6.2f' %(varx1,valx1))
            if valx1 < binx1m or valx1 > binx1M: continue
            iNi = int(varx1.replace('N',''))

            for varx2 in varsx2: 

                valx2 = abs(ipro.mass[varx2])
                if VB>=2: print(VB1+'                  varx2: %s  %6.2f' %(varx2,valx2))
                if valx2 < binx2m or valx2 > binx2M: continue
                iNj = int(varx2.replace('N',''))

                txtdet = '%3s  %3s  %6.2f  %6.2f' %(varx1,varx2,valx1,valx2) 

                xsecT = '%s %s' %(varx1, varx2)
                if xsecT in list(ipro.ox.x.keys()): 
                    y = ipro.ox.x[xsecT] 
                    valsx1.append(valx1)
                    valsx2.append(valx2)
                    if yis in [1,2]: y = 1.
                    ys.append(y)
                    #txtdet += '  %10.6f  %s  %s' %(y, VBdet, VB1)
                    txtdet += '  xsec = %10.6f  || %s - %s || %s' %(y, ipro.getNcomp(iNi,ret=-2), ipro.getNcomp(iNj,ret=-2), VBdet)
                    if VB == -3 or VB >= 2: print(txtdet)
                else: 
                    # y=-1
                    if VB >= 2: print('NB no xsec for: iSS = %3i  iS = %3i   xsecT = %s' %(iSS, iS, xsecT))
                

        # Fill histogram (inside loop)
        if 1: # if histogram
            for iy in range(len(ys)): 
                binx = tax.FindBin(valsx1[iy])
                biny = tay.FindBin(valsx2[iy])
                
                yold = zh.GetBinContent(binx,biny)
                if yold != 0: 
                    nduplicatebins += 1
                    print('Warning: bin (%i,%i)  [%6.2f,%6.2f] already filled. Old: %8.5f   New: %8.5f' %(binx,biny, valsx1[iy],valsx2[iy], yold,ys[iy]))
                    # now refilling with the last entry

                if yis == 2: zh.SetBinContent(binx,biny,ys[iy]+yold)
                else: zh.SetBinContent(binx,biny,ys[iy])
                #print 'FILLING:', txtdet
                nfilledbins += 1

    # Summary
    print('----------------------')
    print('nfilledbins:    %5i' %(nfilledbins))
    print('nduplicatebins: %5i  (last one plotted)' %(nduplicatebins))
    print('widx = %.1f   widy = %.1f' %(binx1d, binx2d))
    # Return
    return zh


#############
def DrawCoverage(List=[], M1=-1, opt={}, txt=''):  # List sent in is e.g. mc11_7TeV['full']
    VB = 1
    
    mrk_col = 2
    mrk_size = 3.5
    mrk_sty = 24
    # OPEN ONES:   24:circle, 25:rectangle, 26:triangle, 27:diamond, 28:cross, 30:star
    # FILLED ONES: 20:circle, 21:rectangle, 22:triangle, 23:triangleupsidedown, 29:star
    # STREK:        2:plus, 3:"star", 5:cross
    mrk_sizes = []

    if 'sty' in opt: mrk_sty = opt['sty']
    if 'col' in opt: mrk_col = opt['col']
    if 'size' in opt:
        if type(opt['size']) in [int,float]: mrk_size = opt['size']
        if type(opt['size']) is list: mrk_sizes = opt['size']


    mk = ROOT.TMarker()
    mk.SetMarkerSize(mrk_size)
    mk.SetMarkerStyle(mrk_sty)
    mk.SetMarkerColor(mrk_col)
    #mk.SetMarkerAlign(22)

    points = List
    ninslice = 0

    pointsdone = []
    
    for p in points:
        zM1,zM2,zMU = p[0],p[1],p[2]

        if zM1 != M1: continue
        if p in pointsdone: continue  # this allows to input as list Dict['full']+Dict['fast']

        if mrk_sizes: # Hack to allow drawing fatter lines :/
            for sz in mrk_sizes:
                mk.SetMarkerSize(sz)
                mk.DrawMarker(zMU,zM2)
        else: 
            mk.DrawMarker(zMU,zM2)


        ninslice += 1
        pointsdone.append(p)


    if VB:
        out = 'DrawCoverage  M1=%i:  npoints = %i  ' %(M1, ninslice)
        if txt: out += '  (%s)' %(txt)
        print(out)


def DrawCoverageAll(Dict={}, M1=-1, mode=1):
    if mode == 0: 
        DrawCoverage(List=Dict['full'],  opt={'col':1, 'sty':20, 'size':1.5}, M1=M1)
        DrawCoverage(List=Dict['fast'],  opt={'col':1, 'sty':20, 'size':1.5}, M1=M1)
        DrawCoverage(List=Dict['full'],  opt={'col':2, 'sty':24, 'size':2.0}, M1=M1, txt='full')
        DrawCoverage(List=Dict['fast'],  opt={'col':3, 'sty':24, 'size':2.5}, M1=M1, txt='fast')
        DrawCoverage(List=Dict['truth'], opt={'col':1, 'sty':25, 'size':2.8}, M1=M1, txt='truth')

    if mode == 1: 
        DrawCoverage(List=Dict['full']+Dict['fast'],  opt={'col':1, 'sty':24, 'size':1.5}, M1=M1)


#############
def getTGraph2D(aScan, what, minmax=(100,350), delta=1e-7): 
    # Quicky which might be more standardised and implemented in the general


    x = array('d')
    y = array('d')
    z = array('d')
    for iSS in range(aScan.N()):
        ipro = aScan.L[iSS]
        # M1 = abs(ipro.isaout['M_1'])
        M2 = abs(ipro.isaout['M_2'])
        MU = abs(ipro.isaout['MU'])

        if 1:
            if not minmax[0] <= MU <= minmax[1]: continue
            if not minmax[0] <= M2 <= minmax[1]: continue

        if 1: 
            # NB: hack to include border values in Interpolate(x,y)
            if M2 == minmax[0]: M2 -= delta
            if MU == minmax[0]: MU -= delta
            if M2 == minmax[1]: M2 += delta
            if MU == minmax[1]: MU += delta
            
        x.append(float(MU))
        y.append(float(M2))

        if 0:
            pass
        elif what == 'N1': z.append(abs(ipro.mass['N1']))
        elif what == 'N2': z.append(abs(ipro.mass['N2']))
        elif what == 'N3': z.append(abs(ipro.mass['N3']))
        elif what == 'N4': z.append(abs(ipro.mass['N4']))
        elif what == 'C1': z.append(abs(ipro.mass['C1']))
        elif what == 'C2': z.append(abs(ipro.mass['C2']))
        elif what in ['C1mN1','C1-N1']: z.append(abs(ipro.mass['C1'])-abs(ipro.mass['N1']))
        elif what in ['N2mN1','N2-N1']: z.append(abs(ipro.mass['N2'])-abs(ipro.mass['N1']))
        else:
            print('WARNING: libIsaplot: getTGraph2D  illegal what: %s' %(what))
            z.append(-1)

    gr = ROOT.TGraph2D(len(x),x,y,z)
    gr.SetName(what); gr.SetTitle(what)
    
    return gr
        

#############
def getLEPlimitsC1(aScan,minmax=(100,350), storeall=0, dMU=1., dM2=1., do2D=0, do1D=1, M1=999): 


    grC1 = getTGraph2D(aScan,what='C1',minmax=minmax)
    grC1mN1 = getTGraph2D(aScan,what='C1-N1',minmax=minmax)
    LEPlimitC1 = LEPlimit_C1_gaugino_func()

    # dMU = 4. # GeV
    # dM2 = 4. # GeV
    

    MU_0 = minmax[0]
    M2_0 = minmax[0]
    MU_1 = minmax[1]
    M2_1 = minmax[1]

    delta = 1e-7

    MU = MU_0 - dMU
    #M2 = M2_0 - dM2  # since inc before

    x = array('d')
    y = array('d')
    z = array('d')


    xx = array('d')
    yy = array('d')

    while MU < MU_1: 
        MU += dMU
        #print 'MU = %.0f' %(MU)
        M2 = M2_0 - dM2
        
        # 2D action    
        while do2D and (M2 < M2_1): 
            M2 += dM2

            delM = grC1mN1.Interpolate(MU,M2)
            LEPlimC1 = LEPlimitC1.Eval(delM)
            C1 = grC1.Interpolate(MU,M2)

            isexcl = 0
            if C1 < LEPlimC1: isexcl = 1

            if storeall or (C1-10 < LEPlimC1 < C1+10): 
                x.append(MU)
                y.append(M2)  
                z.append(isexcl)



        # 1D action: find the y-value (M2-value) where we go from exclusion to non-exclusion)
        cont = True
        ddM2 = 0.1
        M2 = M2_0 - ddM2
        while cont and M2 < M2_1 :
            M2 += ddM2
            delM = grC1mN1.Interpolate(MU,M2)
            LEPlimC1 = LEPlimitC1.Eval(delM)
            C1 = grC1.Interpolate(MU,M2)
            isexcl = 0
            if C1 < LEPlimC1: isexcl = 1

            if not isexcl:
                xx.append(MU)
                yy.append(M2 - ddM2/2)
                cont = False


    # hack: close upper-left
    if 1:
        xx.insert(0, xx[0] - dMU*0.2)
        yy.insert(0, M2_1)

    res = {}
    res['xx'] = xx
    res['yy'] = yy

        
    
    if do2D: res['gr2D'] = ROOT.TGraph2D(len(z), x,y,z)
    
    res['gr1D'] = ROOT.TGraph(len(xx), xx,yy)
    res['gr1D'].SetLineWidth(2)
    res['gr1D'].SetName('LEPlimit_C1_M1eq%.0f' %(M1))
    res['gr1D'].SetTitle('LEPlimit_C1_M1eq%.0f' %(M1))
    res['gr1D'].SetMarkerStyle(2)
    res['gr1D'].SetMarkerSize(0.4)
    
    tax = res['gr1D'].GetXaxis()
    tay = res['gr1D'].GetYaxis()
    tax.SetTitle("MU [GeV]")
    tay.SetTitle("M2 [GeV]")
    #tax.SetTitle("#mu [GeV]")  # 2012-10-05
    #tay.SetTitle("M_2 [GeV]")
    
    return res


############# implementation not finalised
def getFree2D(aScan, nx=100,ny=100, whatx=1, whaty=101): 
    # both whatx and whaty takes the same values 


    
    # loop over and gather
    x = array('d')
    y = array('d')
    
    iS = 0
    for iSS in range(aScan.N()): 
        iS += 1
        ipro = aScan.L[iSS]

        

    
    # prepare histo / TGraph

    # loop and fill


    # return 


#############
#############
