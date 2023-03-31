#!/usr/bin/env python
'''
Program : munch.py
Version : 1.0  2013-09-15
Author  : b.k.gjelsten@fys.uio.no
Description : 

Simple 1D-plotter / 2D-plotter

  input:
    - stdin or filelist(s) (various formats) [or possibly x[,y],z if imported]
      -> builds x[,y],z arrays for plotting

    - config: canvas/plot 
    - config: per curve


  output:
    - pdf


Potential features
  - include also 


====
x,y,z: n-dim lists of floats, where n is the number of quantities/curves



NEXT-TO-DO

  o option: axes on right side as well



2D:

 o generalise the 2D-plotting (currently much is hardcoded)
 o axes are set on hDef; relies on 'same' to be used for plotting [fragile]
 o optionally to specify contour values / palette
 o set contour on top of the default zcont4
 o lowest values of zcont4 go white (top right)




1D:
 x config marker (hide, set size, type, colour, ... default is one in common, same colour as plot)

 

general:
 x Allow combining columns, similar to what is done in TableToTArray. [COULD MAYBE REUSE MOST OF THAT CODE]
   - this is now done in TableToTable.py before munch.py 
   
 o Remove the possibility of using several input files (does that just cause confusion?) ? 

 o Allow coordinates to be other than the first column(s)
   - can already fix externally with using TableToTable




2014-01-29: WILL REMAKE (join with TableToTGraph and extend)

NOW: Histo from TGraph2D
x several plots in series
x contours

o input: ltl or dict  [lists of ltl or dicts to be handled prior to munch]
  o might allow usage of fn_base list, then add endings for as many dicts/tlts as desired [separate algorithm]
o option for markers (typically instead of text)


PROBLEMS: 
o Why can I not get to the outermost edge with colour? Check.


EVENTUALLY
o auto-label contours (a: TLegend, b:manually along some axis
o Move text-algorithm to lib

o (1D, 2D) x (THf,TGraph)



# How to run interactively: 
See lsp:/home/scratch/borgeg/master/susyscans/scan_test1/testinteractive.py



STATUS 2014-01-30:
 - munch.py is now as-good&new-as TableToTGraph (has copied from the input section)
 - remains to test input from several pasted tables and to test pickles (and several pickles)
 - then need to add more features as above + understand why the outer rim is not coloured
 - Then revisit 1D and do the upgrades there as well
 - Then implement colourratio plot to be used by munch

 - Currently munch.py and TableToTGraph.py are both uptospeed: will need to see if TableToTGraph can become a part of munch


FRAGILE:
o reading from table: should chech that all lines have the same number of columns


ISSUES:
o grid no plotted on top of surf


o Allow xyrange to go in and actually change x,y and z arrays (otherwise contours are apparently screwed up)


======= 2014-09-01
1D:
[x] s.VAR -> s.dict['VAR']
[ ]  mcol,msty,msize : could do as for lcol etc. : use as a list, one value per graph
  


'''

import sys,os
import array
import bkgjelstenArgReader
from kilelib import AreNumbers, SaveHistory, DictFromTable, DictlistsFromListOfDicts, GetPlainArray, SaveHistory, WriteToFile
from kilelib_ROOT import NextColourStyle   # this one is causing problems with -h
from kilelib_ROOT2 import BLegends

import ROOT

# ##################################################### GLOBAL METHODS
# ##################################################### GLOBAL METHODS
# ##################################################### GLOBAL METHODS

# ##########
def SetupROOT():
    ROOT.gStyle.SetCanvasColor(0)
    ROOT.gStyle.SetTitleFillColor(0)
    ROOT.gStyle.SetOptTitle(1)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetTitleBorderSize(0)
    ROOT.gStyle.SetPalette(60)
    # palette = ROOT.MakeNullPointer(ROOT.Int_t)
    # palette = 5*(ROOT.Int_t)
    # palette[0] = 15
    # palette[1] = 20
    # palette[2] = 23
    # palette[3] = 30
    # palette[4] = 32
    
    # ROOT.gStyle.SetPalette(10000)

# ##########
def TextsFromRaw(plottexts_raw, ROOT, VB=0):
    # Input format is:  <text>,x,y[,size[D=],color[D=1],font[D=]]:<text>,...
    plottexts = []
    if VB>2: print('plottexts_raw: %s' %(plottexts_raw))
    for zA in plottexts_raw:
        zB = zA.split(',')
        if len(zB) < 3: sys.exit('FATAL: nonallowed -text: %s' %(zA))
        thistext = ROOT.TLatex()
        thistext.SetNDC()
        ztit = zB[0]
        if ztit.startswith("'") and ztit.endswith("'"): ztit = ztit[1:-1]  # hack to remove enclosing 's if any
        thistext.SetTitle(ztit)
        thistext.SetX(float(zB[1]))
        thistext.SetY(float(zB[2]))
        # optional parameters
        # print len(zB), zB
        if len(zB) >= 4: thistext.SetTextSize(float(zB[3]))  # TextSize
        if len(zB) >= 5: thistext.SetTextColor(int(zB[4]))   # TextColor
        if len(zB) >= 6: thistext.SetTextFont(int(zB[5]))    # TextFont
        # then put in array
        plottexts.append(thistext)
    return plottexts

# ##########
class munch: 

    def __init__(s, cmd=[], optD={}):

        # ====================== PRE INIT
        if 'argv' in optD: s.argv = optD['argv']
        else: s.argv = sys.argv


        s.HOME = os.getenv('HOME')
        s.fn_globalhistory = '%s/.globalhistory__munch.txt' %(s.HOME)
        SaveHistory(fn=s.fn_globalhistory, argv=s.argv, opt=['dir','space','stripcmd','date'])
        
        s.cmd = cmd
        s.myname = sys.argv[0].split('/').pop()
        s.VB = 1
        s.cwd  = os.getcwd()  # current work directory  
        #s.dir0 = '%s/XXX' %(s.HOME)
        s.dir0 = ''

        s.undefined = 9999999

        s.dict = {}

        # s.dict['test'] = 'testval'
        # s.dict['test2'] = 'testval2'
        # s.dict['testI'] = 3
        # s.dict['testF'] = 4.34

        s.warn = []
        s.fn_warn = 'warnings.txt'
        s.fn_report = 'report'
        s.report = []

        

        s.rootbatch = 1

        s.fn_table = []
        s.fn_pickles = []
        s.ffn_pickle = ''
        s.readfrompickle = 0

        s.vars2pick = []
        s.moredicts = {}
        

        s.Linput = []  # obsolete

        s.arraytype = {'x':'f', 'y':'f', 'z':'f'}

        

        # ===================== Plotting defaults: 
        
        s.labelx = 'labelx'
        s.labely = 'labely'
        s.labelz = []  # filled by first line in input files, or separately in -labelz
        s.labelz_byhand = []

        s.xtitle = '' # titlex'
        s.ytitle = '' # titley'
        s.ztitle = '' # titlez'

        s.xrangemin = s.undefined
        s.xrangemax = s.undefined
        s.yrangemin = s.undefined
        s.yrangemax = s.undefined
        s.zrangemin = s.undefined
        s.zrangemax = s.undefined

        s.xmin = []
        s.xmax = []
        s.ymin = []
        s.ymax = []
        s.zmin = []
        s.zmax = []

        s.dict['zmin_air'] = 0.10  # "margin" from real max to top of canvas (relative)
        s.dict['zmax_air'] = 0.10  # "margin" from real min to bottom of canvas (relative)

        s.dict['zmin_air_log'] = 2.  # "margin" from max to top of canvas    (division)
        s.dict['zmax_air_log'] = 2.  # "margin" from min to bottom of canvas (multiplication)

        s.ndim = 1
        s.dict['line1is'] = 'auto'


        # === canvas
        s.dict['canv_xpos'] = 10
        s.dict['canv_ypos'] = 10
        s.dict['canv_xwidth'] = 700
        s.dict['canv_ywidth'] = 700

        s.dict['xtitleoffset'] = 1.2
        s.dict['ytitleoffset'] = 1.4

        s.dict['tickyright'] = 1  # 0:no, 1:ticks, 2:ticks&label 
        s.dict['ticktop'] = 0

        s.dict['xLog'] = 0
        #s.dict['xLogMin'] = 0.1
        s.dict['yLog'] = 0
        s.dict['yLogMin'] = 0.1
        s.dict['zLog'] = 0
        s.dict['zLogMin'] = 0.1

        s.dict['gridx'] = 0
        s.dict['gridy'] = 0

        s.dict['label_x1'] = 0.12
        s.dict['label_xwid'] = 0.13
        s.dict['label_y2'] = 0.895
        s.dict['label_DY'] = 0.035
        s.dict['label_maxchar'] = 99
        s.dict['label_maxcharoverflow'] = '..'

        s.dict['legend2right']  = -9. # no effect if negative
        s.dict['legend2bottom'] = -9. 
        s.dict['1D_reorder'] = 1   # D=1;  1:largest first,  -1:smallest first,  0:keep order 

        s.title = ''
        s.plottexts_raw = []
        s.plottexts = []

        s.dict['textprecision'] = 1  # 1 or 2: the number of digits for the text on the canvas
        s.dict['marker'] = 0  # 0:not use ; 1:use
        s.dict['marker_size']   = 1.5
        s.dict['marker_style']  = 24
        s.dict['marker_colour']  = 1

        s.numbermarker = ''  # if set, it will replace the number with given ascii ; the size and colour remain the same
        
        # === curves
        s.gr  = []
        s.hgr = []
        s.gr_cont = []
        s.hgr_cont = {}
        
        s.dict['drawstyle_hdef'] = ''
        s.dict['drawstyle_gr1D'] = 'lp'
        #s.drawstyle_gr2D = 'zcont4same'  # cannot plot contours on top of this  (2014-01-30 changed to surf2)
        s.drawstyle_gr2D = 'zsurf2same'  # is the same as zcont4same if you also do c1.SetTheta(90); c1.SetPhi(0)
        s.dict['CanvasTheta'] = 90
        s.dict['CanvasPhi'] = 0
        s.dict['drawstyle_cont'] = 'cont2same'
        #s.dict['drawstyle_cont'] = 'cont5same'  # test

        s.dict['practicalZERO'] = 1e-15  # useful for log plots

        """
        s.dict['legDY'] = 0.035 # notyetinuse
        s.dict['legX1'] = 0.12
        s.dict['legX2'] = 0.25
        s.dict['legY1'] = 0.60  # Not used. Is rather generated from mlegDY and number of graphs
        s.dict['legY2'] = 0.895
        """

        s.dict['autoColorStyleMode'] = '8:1,3,7:1'
        s.dict['autoColorStyleUse'] = 1
        s.dict['autoColorStyle'] = None

        s.dict['width_gr1D'] = 3

        s.dict['numbersonplot'] = 0

        s.dict['topmargin'] = 0.10
        s.dict['bottommargin'] = 0.10
        s.dict['leftmargin'] = 0.10
        s.dict['rightmargin'] = 0.10

        s.dict['ticklengthx'] = 0.03  #D=0.03
        s.dict['ticklengthy'] = 0.03
        # ===


        s.dict['fn_save'] = 'munch.pdf'
        #s.cols = []

        s.dict['scale'] =  1.
        s.dict['mcol']  = -1
        s.dict['msize'] = -1.  # default size is 1
        s.dict['msty']  = -1

        s.dict['lcol'] = []   # line colour (manual; if set, it overrides autoColorStyle)
        s.dict['lsty'] = []   # line style  ( -- " -- )
        s.dict['lwid'] = []   # line width  ( -- " -- )

        s.x = []  #array(s.arraytype['x'])
        s.y = []  #array(s.arraytype['y'])
        s.z = []  #array(s.arraytype['z'])

        # (some structures taken from TableToTGraph)
        s.xcoord  = ''
        s.xcoordT = ''
        s.ycoord  = ''
        s.ycoordT = ''
        
        s.conts  = []
        s.contsT = {}
        s.contsDef = {}  # contains dict of type {'vals': zvals, 'col':zcol, 'sty':zsty, 'wid':zwid}


        s.resvars  = []
        s.resvarsT = {}

        s.dict['numbers_angle'] = 0
        s.dict['numbers_size'] = 0.030
        s.dict['numbers_colour'] = 1
        s.contval = {}

        s.dict['delim_resvar_comma'] = ','
        s.dict['delim_resvar_colon'] = ':'

        s.dict['saveCtoo'] = 0
        s.dict['operationprotection'] = 0  # 
        s.dict['cutrange'] = 1  # allow cutting range with xyrange
        s.dict['tableverbose_resvar'] = 0  # 1:outputs to screen the effective table plotted (for inspection/post-treatment/...)
        #s.dict['fn_tableverbose_resvar'] = ''

        s.dict['fn_history'] = 'history_munch.txt'
        s.dict['fn_history_global'] = '%s/.history_munch.txt' %(s.HOME)
        SaveHistory(fn=s.dict['fn_history'], argv=sys.argv, opt=['stripcmd','date'])
        SaveHistory(fn=s.dict['fn_history_global'], argv=sys.argv, opt=['stripcmd','date','dirsameline'])

        s.dict['plot_formats'] = ['pdf','eps','ps','jpg','gif','png']
        s.dict['commonpdf'] = 1


        # ====================== READ ARG
        if 'ReadArg' in s.cmd: s.ReadArg()


        # ====================== POST INIT
        if 'PostInit' in s.cmd: s.PostInit()


        # ====================== EXECUTE 
        if 'Main' in s.cmd: s.Main()



    # ##########
    def PostInit(s): 
        if s.dir0: s.fn_warn = '%s/%s' %(s.dir0, s.fn_warn)
        if s.dir0: s.fn_report = '%s/%s' %(s.dir0, s.fn_report)

        #if s.title == '':
        #    if s.dict['fn_save'] != 'munch.pdf': s.title = s.dict['fn_save'].replace('.pdf','')
        #    elif s.fn_table: s.title = s.fn_table[0].replace('.txt','').replace('.dat','')


        # Checks
        if s.ndim not in [1,2]: sys.exit("Need to specify coordinate vars with -coord <1 or 2 coords>")
        
        
        


    # ##################################################### CLASS METHODS
    # ##################################################### CLASS METHODS
    # ##################################################### CLASS METHODS
    # ##########
    def showHelp(s):
        print(" DESCRIPTION")
        print("   munch.py is a plottings script which is intended to simplify plotting with ROOT.")
        print("   The input is either")
        print("     a) a list of pickle files (each file containing a flat python dict which contains the coordinates and result variables of one point)")
        print("     b) or an ordered table with the header line describing the variable names;")
        print("        each line needs to contain the same number of variables.")
        print("          Ex: ")
        print("            mass  height  width  varN  var_M ")
        print("            2.1   4.2      3      4     4.2 ")
        print("            3     4.8      2.1    5.3   8   ")
        print("            8.22   4        34.2  9     9.11")
        print("            ... ")
        print() 
        print("          From this table an x-variable [a y-variable] and a list of result (z) variables are selected.")
        print("   (An arbitrary combination of the variables can also be used as variable to the plot, e.g. height+width or height*width or 2*height+width, ..)")
        print("     Ex: -coord height             : 1D-plot")
        print("         -coord height,width       : 2D-plot")
        print("         -coord height*width,mass  : 2D-plot with height*width as the x-variable")
        print("         -coord 2.1*height*width,mass  : 2D-plot with 2.1*height*width as the x-variable")
        print("         -coord height:H,width:W   : 2D-plot with 'height' and 'width' renamed as 'H' and 'W' in the plot")
        print("         -resvars  var_M           : var_M is the resvar") 
        print("         -resvars  var_M,mass      : var_M and mass are the result variables") 
        print("   For 2D-plots there will be as many plots as there are result variables.")
        print("   For 1D-plots all the result variables will be plotted in the same canvas.")
        print()
        print() 
        print() 
        print()
        print(" OPTIONS:") 
        print("    -coord <varlist>  : ex: -coord mu,M2              : selects these two columns in the input table")
        print("                              -coord mu:#mu,M2:M_{2}    : renames the axes in the plot")
        print("                              -coord 10*mu,M2+M1:a_sum  : coordinates can be (complex) combinations of the variables in the table")
        
        print("    -resvars <varlist> : ex: -resvars DM")

        print("    -title <new canvas title>")
        print("    -save <filename>     :  name of output file (D='munch.pdf'). If none of the allowed formats (%s)" %(s.dict['plot_formats']))
        print("                            is explicitly in the name, a '.pdf' is added.")  
        print("    -xrange xmin,xmax    :  uses a different (typically smaller) x-range than what is given in the table.")
        print("    -yrange ymin,ymax    :  uses a different (typically smaller) y-range than what is given in the table.")

        print("    -xyrange xymin,xymax : same as above, setting equal range for x and y")
        print("    -zrange zmin,zmax    : sets min and max for the result variables. Values below/above min/max are then set to min/max")
        print("    All the rangevars above can also be set onesided (min or max) with:")
        print("       -xrangemin / -yrangemin / -xyrangemin / -zrangemin / -xrangemax / -yrangemax / -xyrangemax / -zrangemax")
        print() 
        print("    --numbers          : [2D] adds numbers (not just markers). Format is intuitive.")
        print("      Format of numbers can be changed with:")
        print("      -numbers_angle <angle> : prints numbers at angle <angle> (D=0)")
        print("      -numbers_size  <size>  : size of numbers (D=0.030)")
        print("      -numbers_colour <col>  : colour of numbers (D=1=black)")
        print("      -numbers_asc <angle>,<size>,<col> : specify the three above in one go")
        print("      -textprecision <0/1/2/3>  : mode for precision in the numbers (D=1)")
        print("    --numbermarker     : [2D] adds a marker 'o' (instead of number) where the table has a value")
        print("    -numbermarker <x>  : [2D] adds a marker <x> (instead of number) where the table has a value")
        print()
        print("    [2D] Can plot contour in any of the variables in the table (on top of the result variable)")
        print("    -cont <var>:<varText>:<value>:<colour>:<style>:<width>")
        print("      Ex: -cont N1::100            : minimal, only give the var name and the value. Col/style/width default to 1,1,1")
        print("      Ex: -cont N1:mass_of_N1:100  : replace var name in the contour info put on the plot")
        print("      Ex: -cont N1:mass_of_N1:100:2:3:4  : also specify colour/style/width of the contour")
        print("    Can have several contours on the same plot (separated by ',,'):")
        print("    -cont <var1><var1Text>:<value>:<colour>:<style>:<width>,,<var2><var2Text>:<value>:<colour>:<style>:<width>")
        #print "    (Note: the contours are sometimes strange (that's a ROOT feature)"
        
        # Note: the numbermarker is simple asc-plotting, the built-in marker is not working with these contour plots (hm)
        #print "      Format of markers can be changed with:"  # 
        #print "      -marker_size <size>     : D=1.5"
        #print "      -marker_style <style>   : D=24=circle, "
        #print "      -marker_colour <colour> : D=1=black"
        #print "      -marker_ssc <style>,<size>,<colour> : all three in one go"
        print()
        print("    -text <some_text>,x,y[,size,col,font]  : Add text to canvas")
        print("        Need to specify <some_text>, typically inside '', then x and y coordinate (in range [0,1])")
        print("        Optional is to specify size (D=0.05), colour (D=1(=black)) and font")
        print("     Ex:  -text 'my testplot',0.8,0.2      # minimal")
        print("     Ex:  -text testplot,0.8,0.2,0.3       # size")
        print("     Ex:  -text testplot,0.8,0.2,0.3,2,42  # size,colour,font")
        print("     Ex:  -text testplot,0.8,0.2,0.3,2,42:'some more text',0.3,0.8    # another text, separated by :")
        
        
        print()
        print()
        print("  MORE OPTIONS:")
        print("    -h                : this help message")
        print("    --h               : brute info on implemented options")
        print("    --dict            : dump various default values (which can be changed)")
        print("    -scale <val>      : scale the result values by <val>")
        print("    -logx 0/1          , --logx   (D=%i)" %(s.dict['xLog']))
        print("    -logy 0/1          , --logy   (D=%i)" %(s.dict['yLog']))
        print("    -yLogMin <val>  (if on)  (D=%.3f)" %(s.dict['yLogMin']))
        print("    -zLogMin <val>  (if on)  (D=%.3f)" %(s.dict['zLogMin']))
        print("    -rightticky 0/1    , --rightticky   (D=%i)" %(s.dict['tickyright']))
        print("    -toptickx 0/1      , --toptickx     (D=%i)" %(s.dict['ticktop']))
        print("    -ticklengthx <val> ,  -ticklengthy <val>,  -ticklengthxy <val>   (D(x)=%.3f, D(y)=%.3f" %(s.dict['ticklengthx'], s.dict['ticklengthy']))
        print("    -gridx <val>       , --gridx")
        print("    -gridy <val>       , --gridy")
        print("    -gridxy <val>      , --gridxy")
        print("    -label_x1 <val>  (D=%.3f), -label_xwid <val>  (D=%.3f), -label_y2 <val>  (D=%.3f), -label_DY <val> (D=%.3f)" %( s.dict['label_x1'], s.dict['label_xwid'], s.dict['label_y2'], s.dict['label_DY']))
        print("    -label_maxchar <val>  (D=%i) , -label_maxcharoverflow <val> (D='%s')" %(s.dict['label_maxchar'], s.dict['label_maxcharoverflow']))
        #print "    -numbers"
        print("    -topmargin <val> , -bottommargin <val> , -leftmargin <val> , -rightmargin <val>    : change margins (D=0.10 for all).")
        print("        (Especially it may be relevant to increase the right margin to get the full colour-scale.)")
        print() 
        #print "    -cols <colnumbers/titles>"
        print()
        #print "  Some dict options: set with e.g.  -dict I,mcol,3:I,msize:2  [I=integer, F=float]"   # These do not work with contours
        #print "    I,mcol:  marker colour (1D)  D=linecolour"
        #print "    I,msize: marker size (1D), D=1"
        #print "    I,mstyle: marker style (1D), D=1"
        print() 
        #print "  - the x [and y] axis is the first [and second] column"
        #print "  - 1D: default is to plot all columns in same plot "
        #print
        #print "  - ROOT marker styles:  http://root.cern.ch/root/html520/MACRO_TAttMarker_3_c.gif"
        print()
        print(' EXAMPLES 2D') 
        #print '        %s  -dict test,txt1:I,testI,3:test2,txt2a,txt2b:F,testF,4.14   # for using autodict (NB: vars need to be defined in __init__)' %(s.myname)
        #print "    Ex: munch.py  -ndim 2  -f wide0.txt  -cols 2,1,4   # picks out cols for x,y,z" 

        #print "    Ex: munch.py -f fn_input_1d_2  -zrange 100,800  --autocolstyle   -titles '#mu [GeV]','masses [GeV]'  -label_DY 0.05  -label_xwid 0.25  -text 'M_{1}=M_{2}=100 GeV',0.05,0.02,0.04,2,42  -ctitle DGnoSL  -labelz '#tilde{#chi}_{3}^{0}'"
        #print "    Ex: munch.py -f xsec1D_MUeq3E3.txt,xsec1D_M2eq3E3.txt -titles 'M2/MU [GeV]','Total LO cross-section [pb]'  -text M_{1}='50GeV #mu/M2=3TeV',0.05,0.02,0.035,1,42 -ctitle DGnoSL  -fn_save xsec1D.pdf  --logy -vb 3  --grid -labelz  'bino+wino,bino+higgsino'  -label_xwid 0.25  -label_DY 0.06  # 2013-09-24"
        #print "    Ex: munch.py  -ndim 2  -f xsec2D_MU_M2_100_500.txt  -titles '#mu [GeV]','M2 [GeV]','LO Cross-section [pb]'  -text M_{1}='50GeV',0.05,0.02,0.035,1,42  -ctitle 'DGnoSL: Total LO Cross-section [pb]'  -fn_save xsec2D.pdf  -xyrange 100,500   -style2D  zcont4 --logz  -numbers 1  -vb 2  -right 0.12  -ticklengthxy -0.008   # 2013-09-24"
        
        #print "    ex: cat table.txt | munch.py  -ndim 1  -cols M2,3,4  -yrange 0,1  -ztitle 'Branching Ratio'  -title DGnoSL"
        #print "    ex: cat table.txt | munch.py  # defaults to 1D, takes all columns"
        #print 
        #print "  A typical command line could include the following"
        #print "    ls *slha | isa_simplereader_class.py [opts] | TableToTable.py [opts] | munch.py [opts]"

        print("   Ex: cat a_table | munch.py -coord mu,M2  -resvars N1          # plots mass of N1 in 2D mu-M2 plot  (table piped in: the output is munch.pdf)")
        print("   Ex: munch.py -coord mu,M2  -resvars N2  -fn_table a_table     # the table can also be input as a file. The output is now a_table.pdf") 
        print("   Ex: munch.py -coord mu,M2  -resvars N2another  -fn_table a_table,another_table     # several tables can be combined (in not all the required info is in the same table. Obviously the tables need to be ordered in the same way.")
        print("   Ex: ls *pickle | munch.py -coord mu,M2  -resvars N2   --pickle                         # if the info is in pickle files")
        print("   Ex: munch.py -coord mu,M2  -resvars N2   -ffn_pickle <file_with_list_of_pickles>       # (similar to the above)") 
        print()
        #print " EXAMPLES"
            
        print("    Ex: cat table1.txt | munch.py  -coord 'MU:#mu [GeV],M2:M_{2} [GeV]'   -resvars DMbestVSexp,DMwoCo") 
            
        print("    Ex: munch.py -coord mu:#mu[GeV],M2:M_{2}[GeV]   -resvars ~N2:N2  --mode2Da  -numbers 0  -vb 1   -cont ~C1::100,200:1:1:4,,~N3::300,400:2:2:2,,~N4::300,400:3:3:3  -fn_table table_elweak.txt") 
        print()
        print("  EXAMPLES 1D:")
        print("       cat table1D.txt | munch.py  -coord mu  -resvars h")
        print("       cat table1D.txt | munch.py  -coord mu  -resvars h  -zrange 110,130                        # To specify min/max in plot")
        print("       cat table1D.txt | munch.py  -coord mu:'#mu [GeV]'  -resvars h  -zrange 110,130            # Redefine text on x-axis")
        print("       cat table1D.txt | munch.py  -coord mu  -resvars h  -zrange 110,130  -ztitle 'm(h) [GeV]'  # Set title on z-axis")  
        print("       cat table1D.txt | munch.py  -coord mu  -resvars h:Higgs  -zrange 110,130                  # Redefine varname in legend")
        print("       cat table1D.txt | munch.py  -coord mu  -resvars h  -zrange 110,130  -title 'Higgs mass'   # Set title on plot")
        print("       cat table1D.txt | munch.py  -coord mu  -resvars h  -zrange 110,130  -msty 20                          # Put a marker on datapoints. Try  20-30, 1-7. See http://root.cern.ch/root/html534/gif/markers.gif")
        print("       cat table1D.txt | munch.py  -coord mu  -resvars h  -zrange 110,130  -msty 20  -msize 1.8              # Marker size")
        print("       cat table1D.txt | munch.py  -coord mu  -resvars h  -zrange 110,130  -msty 20  -msize 1.8  -mcol 3     # Marker colour")
        print("       cat table1D.txt | munch.py  -coord mu  -resvars h  -zrange 110,130  -marker_ssc1D 20,1.8,3            # Marker style,size,colour in one go")
        print("       cat table1D.txt | munch.py  -coord mu  -resvars h  -zrange 110,130  -lcol 5  -lsty 2  -wid 5          # Specify colour, style and width of line.")
        
        print("       cat table1D.txt | munch.py  -coord mu  -resvars h,~N1  -zrange 60,130       # Plot both h and ~N1 in same plot; will automatically get different colours ")
        print("       cat table1D.txt | munch.py  -coord mu  -resvars h,~N1  -zrange 60,130  -lcol 6,8  -lsty 2,9  -lwid 4,2  #  Specify colour/style/width of the two graphs")

        
        print()
        print("[Some/more documentation:  munch.py -h / --h / --dict]")
        # NOTE NOTE NOTE: the '#' above irritates the emacs indentation 

        
    # ##########
    def DumpWarnings(s):
        f = open(s.fn_warn,'w')
        for out in s.warn: f.write('%s\n' %(out))
        f.close()
        
    # ##########
    def Main(s):
        # if s.VB: print "INFO::%s  Main" %(s.myname)
        # for key in s.dict.keys(): print 'dict:  %-10s  %s' %(key, s.dict[key])

        s.PrepareXYZ()
        s.Plot()


    # ##########
    def Plot(s):

        # if s.dict['autoColorStyleUse']:
        #     s.dict['autoColorStyle'] = NextColourStyle(argtext=s.dict['autoColorStyleMode'], blackfirst=1)

        if s.ndim == 1:
            s.Plot1D() # make external loop here as well (probably)
            
        elif s.ndim == 2:
            #for iresvar in range(len(s.resvars)):
            #resvar = s.resvars[iresvar]
            s.Plot2D()
        
        


        

    # ##########
    def PrepareXYZ(s):

        # 1) Create the Dictlist (=tabledict)
        
        if not s.readfrompickle:
            # Case A:  Read table(s)
            # a) Read the table
            if not s.fn_table: 
                if s.VB > 0: print('INFO::munch::PrepareXYZ  Reading table names from stdin')
                lines = sys.stdin.readlines()  # This reads the fileNAME(s) of input tables (not the tables)
            else: 
                cmd = 'paste '
                for z in s.fn_table: cmd += '  ' + z.strip()
                # a) Get the table
                lines = os.popen(cmd).readlines()

            # b) Transfrom table to dict (then don't need the table any longer)
            s.tabledict = DictFromTable(lines)

        else:
            # Case B: Read list of dicts (from stdin only for now) 
            s.fns_pickle = sys.stdin.readlines()

            #print 'AUF fns_pickle: ', s.fns_pickle
            #print 'AUF vars2pick:  ', s.vars2pick
            
            s.tabledict = DictlistsFromListOfDicts(s.fns_pickle, first=1, vars2pick=s.vars2pick)

            for replacewiththis in s.moredicts:  # NOT YET IMPLEMENTED
                # usually one additional dictionary (or of course zero)
                # Create filename list of the additional pickles
                replacethis = s.moredicts[replacewiththis]
                
                fns_pickle = []

                # E.g. replace all dicts of type  DGnoSL_*_flatdict1.pickle -> DGnoSL_*_flatdict2.pickle by -moredict flatdict1:flatdict2
                for fn in s.fns_pickle:
                    fn_dict2 = fn.strip().replace(replacethis, replacewiththis)
                    if not os.path.exists(fn_dict2): sys.exit("Fatal::PrepareXY  Non-existent fn_dict2: %s" %(fn_dict2))
                    fns_pickle.append(fn_dict2)
                # Get the new (separate) Dictlist
                tabledict2 = DictlistsFromListOfDicts(fns_pickle, first=1, vars2pick=s.vars2pick)   # usage of 'first' is obsolete. should always be 1
                # Update the existing with the new
                s.tabledict.update(tabledict2)
                del tabledict2


        # result so far: s.tabledict


        # 2) Get the 2-3 arrays from the Dictlist
        # tja ... maybe just wait until plotting

        if s.VB>1: print("len(s.tabledict) = %i" %(len(s.tabledict)))
        if s.VB>2: print("s.tabledict: %s" %(s.tabledict))


        if s.VB>3: print('DEBUG:  s.xcoord: %s' %(s.xcoord))
        s.x = GetPlainArray(table=s.tabledict, var=s.xcoord, arraytype='f', VB=s.VB, protection=s.dict['operationprotection'])
        if s.ndim == 2:
            if s.VB>3: print('DEBUG:  s.xcoord: %s,  s.ycoord: %s' %(s.xcoord, s.ycoord))
            s.y = GetPlainArray(table=s.tabledict, var=s.ycoord, arraytype='f', protection=s.dict['operationprotection'])
        s.nvars = len(s.x)


        # Brute test if coordinates seem ok (no duplicates points)
        coordproblem = 0
        zzz = []
        zouts = []
        for i in range(s.nvars): 
            if s.ndim == 1: zz = s.x[i]
            else: zz = [s.x[i], s.y[i]]
            
            zout = '    coord  %2i  %s' %(i, zz)
            if zz in zzz: 
                zout += '    duplicate!' 
                coordproblem += 1
            else:
                zzz.append(zz)
            zouts.append(zout)

            #if coordproblem: print 'Coordinate problem: s.x: ', s.x 


        if coordproblem:
            print('WARNING  coordproblem = %i' %(coordproblem))
            print('  s.x (%i): %s' %(len(s.x), s.x))
            if s.ndim == 2: print('  s.y (%i): %s' %(len(s.y), s.y))
            for zout in zouts: print(zout)

        del zzz, zouts


        if min(s.x) == max(s.x):
            coordproblem += 1000
            print('Coordinate problem: min(s.x) == max(s.x) = %s' %(min(s.x)))
            print('  s.x = %s' %(s.x))
        if s.ndim == 2 and min(s.y) == max(s.y):
            coordproblem += 2000
            print('Coordinate problem: min(s.y) == max(s.y) = %s' %(min(s.y)))
            print('  s.y = %s' %(s.y))


        if coordproblem: 
            print('(IS NOW LIKELY TO CRASH OR BE ERRONEOUS)')

    def CalcSig(s):
        bkg = s.tabledict["nbkg"]
        sig = s.tabledict["nsig"]
        sign = []
        for bkg_i, sig_i in zip(bkg, sig):
            print(bkg_i, sig_i)
            sign.append(ROOT.RooStats.NumberCountingUtils.BinomialObsZ(int(bkg_i)+int(sig_i),int(bkg_i),0.2))
        s.tabledict["sign"] = sign

    # ##########
    def Plot2D(s):
        
        # 0) COMMON
        # Define canvas
        SetupROOT()
        if s.rootbatch: ROOT.gROOT.SetBatch(True)
        s.c1 = ROOT.TCanvas('sc1','sc1',s.dict['canv_xpos'],s.dict['canv_ypos'],s.dict['canv_xwidth'],s.dict['canv_ywidth'])
        s.c1.SetTopMargin(s.dict['topmargin'])
        s.c1.SetBottomMargin(s.dict['bottommargin'])
        s.c1.SetLeftMargin(s.dict['leftmargin'])
        s.c1.SetRightMargin(s.dict['rightmargin'])

        s.c1.SetTheta(s.dict['CanvasTheta'])
        s.c1.SetPhi(s.dict['CanvasPhi'])
        
        if s.xrangemin == s.undefined: s.xrangemin = min(s.x)
        if s.xrangemax == s.undefined: s.xrangemax = max(s.x)
        if s.yrangemin == s.undefined: s.yrangemin = min(s.y) # 2D
        if s.yrangemax == s.undefined: s.yrangemax = max(s.y) # 2D

        hDef = ROOT.TH2F('hDef', s.title, 1, s.xrangemin, s.xrangemax, 1, s.yrangemin, s.yrangemax)
        s.hDef = hDef
        # hDef.SetMinimum(s.zrangemin) # No effect here, need to do for hgr below (probably also other features)
        # hDef.SetMaximum(s.zrangemax)
        
        tax = hDef.GetXaxis()
        tay = hDef.GetYaxis()
        tax.SetTitle(r"\tilde \chi_{2}")
        tax.SetTitleSize(0.05)
        tax.SetLabelSize(0.035)
        tay.SetTitle(r"\tilde \chi_{1}")
        tay.SetTitleSize(0.05)
        tay.SetLabelSize(0.035)
        # print tax.GetTitleOffset()
        # print tay.GetTitleOffset()
        tax.SetTitleOffset(0.8)
        tay.SetTitleOffset(1)
        # print tax.GetTickLength()
        tax.SetTickLength(s.dict['ticklengthx'])  # no effect
        tay.SetTickLength(s.dict['ticklengthy'])


        # Contours
        for cont in s.conts:
            # Local vars: x,y,nvars (do this because they can be cut with cutrange
            #ndim = len(x)
            
            contT = s.contsT[cont]
            s.CalcSig()
            s.contval[cont] = GetPlainArray(table=s.tabledict, var=cont, arraytype='f', protection=s.dict['operationprotection'])  # should we allow this to also be scaled?

            # ------ hard xy-range (remove from table) [1. contours]
            x = array.array('f', list(s.x))
            y = array.array('f', list(s.y))
            if s.dict['cutrange']: 
                for i in range(s.nvars-1,-1,-1):  # iterate backwards
                    if not ( s.xrangemin <= x[i] <= s.xrangemax  and  s.yrangemin <= y[i] <= s.yrangemax ): 
                        x.pop(i); y.pop(i)
                        s.contval[cont].pop(i)
            # -----
            gr = ROOT.TGraph2D(len(x), x, y, s.contval[cont])  # needs 'f'
            gr.SetName(contT)
            s.gr_cont.append(gr)  # needed?
            hgr_cont = gr.GetHistogram()
            #a = s.contsDef[cont]['vals']
            #print type(list(a)), list(a)
            hgr_cont.SetContour(len(s.contsDef[cont]['vals']), s.contsDef[cont]['vals'])
            hgr_cont.SetLineColor(s.contsDef[cont]['col'])
            hgr_cont.SetLineStyle(s.contsDef[cont]['sty'])
            hgr_cont.SetLineWidth(s.contsDef[cont]['wid'])
            s.hgr_cont[cont] = hgr_cont
            

        for igr in range(len(s.resvars)):
            s.c1.Clear()
            
            # Get the z-values and define the TGraph2D
            resvar  = s.resvars[igr]
            resvarT = s.resvarsT[resvar]
            s.z = GetPlainArray(table=s.tabledict, var=resvar, arraytype='f', scale=s.dict['scale'], protection=s.dict['operationprotection'])
            # (could have used z instead of s.z)

            # ------ hard xy-range (remove from table) [2. graphs]
            x = array.array('f', list(s.x))
            y = array.array('f', list(s.y))
            if s.dict['cutrange']: 
                for i in range(s.nvars-1,-1,-1):  # iterate backwards
                    if not ( s.xrangemin <= x[i] <= s.xrangemax  and  s.yrangemin <= y[i] <= s.yrangemax ): 
                        x.pop(i); y.pop(i)
                        s.z.pop(i)
                    #else: print '%7.2f  %7.2f   %14.8f' %(x[i], y[i], s.z[i])  # Debugging 
            # -----

            #if s.dict['tableverbose_resvar']:   # 2014-02-20
            outs_table = ['  %4s  %10s  %10s  %14s' %('igr', s.xcoord, s.ycoord, resvar)]
            for i in range(len(x)):
                outs_table.append('  gr:%i  %10.3f  %10.3f  %14.5f' %(igr, x[i], y[i], s.z[i]))
            #for out in outs: print out

            gr = ROOT.TGraph2D(len(x), x, y, s.z)
            gr.SetName(resvarT)
            
            
            #grClone = ROOT.TGraph2D(len(x), x, y, s.z) ; grClone.SetName('grClone')  # DEBUGGING
            
            s.gr.append(gr)   # needed(?) to preserve?

            fn = s.dict['fn_save']
            if fn == '': fn = 'munch.pdf' 
            if fn == '': fn = 'munch_%s.pdf' %(resvarT)
            if len(s.resvars) > 1: fn = fn.replace('.pdf','_%s.pdf' %(resvarT))

            if s.zrangemin == s.undefined: s.zrangemin = min(s.z)
            if s.zrangemax == s.undefined: s.zrangemax = max(s.z)

            if s.dict['xLog'] and s.xrangemin <= 0: s.xrangemin = s.dict['xLogMin']
            if s.dict['yLog'] and s.yrangemin <= 0: s.yrangemin = s.dict['yLogMin']
            if s.dict['zLog'] and s.zrangemin <= 0: s.zrangemin = s.dict['zLogMin']


            # 2D
            # Plot
            hDef.SetTitle(resvarT)
            if s.title != "": hDef.SetTitle(s.title)
            hDef.Draw(s.dict['drawstyle_hdef'])
            if s.dict['tickyright']: ROOT.gPad.SetTicky(s.dict['tickyright'])  # 2014-02-06 uncommented (why was it commented?)
            if s.dict['ticktop']:   ROOT.gPad.SetTickx(s.dict['ticktop'])      # ditto
            

            if s.VB>3: print('DEBUG: len(s.nvars): %i, %i' %(s.nvars, len(s.gr)))
            #for igr in range(s.nvars):
            #gr = s.gr[igr]

            
            hgr = gr.GetHistogram()   # using the histogram (40x40) has no effect/gain
            s.hgr.append(hgr)

            # Very inelegant having to do this anew ...  (do I?)
            #hgr.SetTitle(s.title)
            #hgr.SetTitle(resvarT)
            tax2 = hgr.GetXaxis()
            tay2 = hgr.GetYaxis()
            tax2.SetTickLength(s.dict['ticklengthx'])
            tay2.SetTickLength(s.dict['ticklengthy'])
            tax2.SetTitle(s.xtitle)
            tay2.SetTitle(s.ytitle)
            tax2.CenterTitle()
            tay2.CenterTitle()
            tax2.SetTitleOffset(s.dict['xtitleoffset'])
            tay2.SetTitleOffset(s.dict['ytitleoffset'])
            tax2.SetTickLength(s.dict['ticklengthx'])  # no effect
            tay2.SetTickLength(s.dict['ticklengthy'])

            hgr.SetMinimum(s.zrangemin)
            hgr.SetMaximum(s.zrangemax)
        
            hgr.Draw(s.drawstyle_gr2D)  # 
            # Contour plotting
            
            tlegcont = ROOT.TLegend(0.5,0.5, 0.8,0.8)
            for cont in s.conts:
                hgr_cont = s.hgr_cont[cont]
                hgr_cont.Draw(s.dict['drawstyle_cont'])
                
            # Numbers on plot (manual)
            # ----------
            if s.dict['marker']: # 2014-02-07  BUT DOES NOT WORK ON TOP OF SURF.. rather plot a dot with TLatex??
                marker = ROOT.TMarker()
                marker.SetMarkerSize(s.dict['marker_size'])
                marker.SetMarkerColor(s.dict['marker_colour'])
                marker.SetMarkerStyle(s.dict['marker_style'])
                
            if s.dict['numbersonplot']:
                s.plotnumber = ROOT.TLatex()
                s.plotnumber.SetNDC()
                s.plotnumber.SetTextSize(s.dict['numbers_size'])
                s.plotnumber.SetTextAlign(22)
                s.plotnumber.SetTextAngle(s.dict['numbers_angle'])
                s.plotnumber.SetTextFont(42)
                s.plotnumber.SetTextColor(s.dict['numbers_colour'])
                
                xm = s.xrangemin
                xM = s.xrangemax
                ym = s.yrangemin
                yM = s.yrangemax
                L = s.c1.GetLeftMargin()
                R = s.c1.GetRightMargin()
                T = s.c1.GetTopMargin()
                B = s.c1.GetBottomMargin()
                W = 1.-L-R
                H = 1.-T-B
                for i in range(len(x)):
                    Tx = x[i]  # T is just for T[ext], to distinguish it from x
                    Ty = y[i]
                    xNDC = L + W*(Tx-xm)/(xM-xm)
                    yNDC = B + H*(Ty-ym)/(yM-ym)


                    if s.numbermarker:  # hack to plot a textmarker instead of the number (since TMarker does not work over surf)
                        s.plotnumber.DrawLatex(xNDC, yNDC, s.numbermarker)
                        continue

                    # Hack to change alignment at the four borders (to not write over the axis)
                    xal = 20; yal = 2
                    if 1: 
                        if Tx == s.xrangemin: xal = 10
                        if Tx == s.xrangemax: xal = 30
                        if Ty == s.yrangemin: yal = 1
                        if Ty == s.yrangemax: yal = 3
                    s.plotnumber.SetTextAlign(xal+yal)

                    Tz = s.z[i]
                    if s.dict['textprecision'] == 0: 
                        ztxt = '%.0f' %(Tz)
                        if Tz<1: ztxt = '%.1f' %(Tz)
                        if Tz<0.1: ztxt = '%.0f' %(Tz)  # interesting
                        if Tz==0: ztxt = '%.0f' %(Tz)
                        
                    elif s.dict['textprecision'] == 1: 
                        ztxt = '%.0f' %(Tz)
                        if Tz<1: ztxt = '%.1f' %(Tz)
                        if Tz<0.1: ztxt = '%.2f' %(Tz)
                        if Tz==0: ztxt = '%.0f' %(Tz)  # Note
                    elif s.dict['textprecision'] == 2: 
                        ztxt = '%.0f' %(Tz)
                        if Tz<10: ztxt = '%.1f' %(Tz)
                        if Tz<1: ztxt = '%.2f' %(Tz)
                        if Tz<0.1: ztxt = '%.3f' %(Tz)
                        if Tz==0: ztxt = '%.0f' %(Tz)  # Note
                    elif s.dict['textprecision'] == 3: 
                        ztxt = '%.1f' %(Tz)
                        if Tz<10: ztxt = '%.1f' %(Tz)
                        if Tz<1: ztxt = '%.2f' %(Tz)
                        if Tz<0.1: ztxt = '%.3f' %(Tz)
                        if Tz==0: ztxt = '%.0f' %(Tz)  # Note
                    else:
                        sys.exit("Error: Non-allowed value dict['textprecision'] = %s" %(s.dict['textprecision']))

                    s.plotnumber.DrawLatex(xNDC, yNDC, ztxt)
                    #print 'Plotting (%.1f,%.1f) %s' %(xNDC,yNDC,ztxt)

                    if s.dict['marker']:  # 2014-02-07  DOES NOT WORK ON SURF ; rather add -numbermarker <textmarker, e.g. o> [or just --numbermarker for using 'o']
                        marker.DrawMarker(Tx,Ty)
                    
            # ----------
                    
            

            # Labels/Legends
            # ... 

            # Logarithmic?, grid?
            if s.dict['xLog']: s.c1.SetLogx(s.dict['xLog'])
            if s.dict['yLog']: s.c1.SetLogy(s.dict['yLog'])
            if s.dict['zLog']: s.c1.SetLogz(s.dict['zLog'])
            s.c1.SetGrid(s.dict['gridx'],s.dict['gridy'])
            #print 'gridding: ', s.dict['gridx'], s.dict['gridy']  # doesn't work with SURF

            if s.conts:
                #tlegcont.Draw('L')
                # WORK-IN-PROGRESS
                pass

            # Additional text on the plot
            s.plottexts = TextsFromRaw(s.plottexts_raw, ROOT, VB=s.VB) 
            for plottext in s.plottexts: plottext.Draw()

            
            # Save
            s.c1.SaveAs(fn)

            if s.dict['saveCtoo']:
                fnC = fn.replace('.pdf','').replace('.eps','') + '.C'
                s.c1.SaveAs(fnC)
        
            if s.dict['tableverbose_resvar']:
                fn_table = fn.replace('.pdf','').replace('.eps','') + '_table.txt'
                WriteToFile(fn=fn_table, outs=outs_table, VB=s.VB-1)
                
        
    # ##########
    def Plot1D(s):

        # Define canvas, hDef
        SetupROOT()
        if s.rootbatch: ROOT.gROOT.SetBatch(True)
        s.c1 = ROOT.TCanvas('sc1','sc1',s.dict['canv_xpos'],s.dict['canv_ypos'],s.dict['canv_xwidth'],s.dict['canv_ywidth'])
        s.c1.SetTopMargin(s.dict['topmargin'])
        s.c1.SetBottomMargin(s.dict['bottommargin'])
        s.c1.SetLeftMargin(s.dict['leftmargin'])
        s.c1.SetRightMargin(s.dict['rightmargin'])

        s.c1.SetTheta(s.dict['CanvasTheta'])
        s.c1.SetPhi(s.dict['CanvasPhi'])

        if s.xrangemin == s.undefined: s.xrangemin = min(s.x)
        if s.xrangemax == s.undefined: s.xrangemax = max(s.x)
        if s.xtitle == '': s.xtitle = s.labelx

        hDef = ROOT.TH1F('hDef', s.title, 1, s.xrangemin, s.xrangemax)
        if s.title != "": hDef.SetTitle(s.title)  # not effective

        tax = hDef.GetXaxis()
        tay = hDef.GetYaxis()
        #tax.SetTitle(s.xtitle)
        tax.SetTitle(s.xcoordT)
        tay.SetTitle(s.ztitle)
        
        tax.CenterTitle()
        tay.CenterTitle()
        tax.SetTitleOffset(s.dict['xtitleoffset'])
        tay.SetTitleOffset(s.dict['ytitleoffset'])




        # Define graphs (not yet plotting)
        zrangemin =  1e9  # will be overwritten
        zrangemax = -1e9  # will be overwritten

        resvarsT = []
        gr_reorder = []   # order on size (optional)
        for igr in range(len(s.resvars)):

            resvar  = s.resvars[igr]
            resvarT = s.resvarsT[resvar]
            resvarsT.append(resvarT)
            s.z = GetPlainArray(table=s.tabledict, var=resvar, arraytype='f', scale=s.dict['scale'], protection=s.dict['operationprotection'])

            # ------ if logz: find non-zero minimum and replace zeros and negative with 'practicalZERO' (and give warning)
            nonzeromin = 9999999
            if s.dict['zLog']: 
                for i,zdum in enumerate(s.z): 
                    if zdum <= 0: 
                        if s.VB>0: print('Warning: var %i, %s : replacing value %i/%i:   %8.2e -> %8.2e' %(igr, resvarT, i, len(s.z), s.z[i], s.dict['practicalZERO']))
                        s.z[i] = s.dict['practicalZERO']
                    else: 
                        nonzeromin = min(nonzeromin, s.z[i])
            if nonzeromin == 9999999: nonzeromin = s.dict['practicalZERO']  # no nonzero minimum found ... 

            gr_reorder.append([sum(s.z), igr])   # (using the sum as ordering parameter)
 
            # ------ min/max
            if s.zrangemin == s.undefined: 
                if s.dict['zLog']: zrangemin = min(nonzeromin, zrangemin)  # hack for log-plots
                else: zrangemin = min(min(s.z), zrangemin)
            else: zrangemin = s.zrangemin
            if s.zrangemax == s.undefined: zrangemax = max(max(s.z), zrangemax)
            else: zrangemax = s.zrangemax



            # ------ hard xy-range (remove from table) [2. graphs]
            x = array.array('f', list(s.x))
            #y = array.array('f', list(s.y))
            if s.dict['cutrange']: 
                for i in range(s.nvars-1,-1,-1):  # iterate backwards
                    if not ( s.xrangemin <= x[i] <= s.xrangemax ): #  and  s.yrangemin <= y[i] <= s.yrangemax ): 
                        x.pop(i)  #; y.pop(i)
                        s.z.pop(i)
                    #else: print '%7.2f  %7.2f   %14.8f' %(x[i], y[i], s.z[i])  # Debugging 


            # hack to sort according to xaxis ..
            s.z = array.array('f', list(zip(*sorted(zip(x,s.z))))[1])   # wow...
            x = array.array('f', sorted(x))

             

            #if s.VB>2: print '1D  igr: %i   len(s.x): %i   len(s.z): %i' %(igr, len(s.x), len(s.y))
            if s.VB>3: print('  s.x: %s' %(s.x))
            if s.VB>3: print('  s.z: %s' %(s.z))

            if s.dict['zLog'] and min(s.z) <= 0: 
                print("Warning log y-axis and negative values: %i  %s : %s" %(igr, resvarT, s.z))


            gr = ROOT.TGraph(len(x), x, s.z)
            gr.SetName(resvarT)

            #gr.SetFillStyle(0)  # Not sure why I need to specify this  # needed to because I used default sty ('lpf') (f=fill)
            s.gr.append(gr)

        # done looping over res vars: fixing arrays, creating graphs, nicify them



        # Optional: reorder on size (largest first)
        if s.dict['1D_reorder'] ==  1: gr_reorder.sort(reverse=True) # largest first
        if s.dict['1D_reorder'] == -1: gr_reorder.sort() # smallest first
        # else: no reordering, keep the input order
        gr_reordered = []
        resvarsT_reordered = []
        for zdum,igr in gr_reorder: 
            gr_reordered.append(s.gr[igr])
            resvarsT_reordered.append(resvarsT[igr])



        # ----
        # Configure graph (colour, style, width, marker)
        for igr,gr in enumerate(gr_reordered):
 
            if s.dict['autoColorStyle']:
                [zcol,zsty] = next(s.dict['autoColorStyle'])
                if len(s.dict['lcol']) >= igr+1: zcol = s.dict['lcol'][igr]
                if len(s.dict['lsty']) >= igr+1: zsty = s.dict['lsty'][igr]
                
                gr.SetLineColor(zcol)
                gr.SetLineStyle(zsty)
                gr.SetMarkerColor(zcol)
                if s.dict['mcol']  > -1: gr.SetMarkerColor(s.dict['mcol'])
                if s.dict['msize'] > -1: gr.SetMarkerSize(s.dict['msize'])
                if s.dict['msty']  > -1: gr.SetMarkerStyle(s.dict['msty'])
                del zcol, zsty
                if s.VB>1: print("INFO::Plot1D  autoColorStyle : col:%i, size:%.3f, style:%i" %(gr.GetMarkerColor(), gr.GetMarkerSize(), gr.GetMarkerStyle()))

            # Overrides autoColorStyle (if present)
            #if len(s.dict['lcol']) >= igr+1: gr.SetLineColor(s.dict['lcol'][igr])
            #if len(s.dict['lsty']) >= igr+1: gr.SetLineStyle(s.dict['lsty'][igr])


            if len(s.dict['lcol']) >= igr+1: zcol = s.dict['lcol'][igr]  # need to redo in case of not s.dict['autoColorStyle'] (inelegant)
            if len(s.dict['lsty']) >= igr+1: zsty = s.dict['lsty'][igr]
            gr.SetLineWidth(s.dict['width_gr1D'])  # 1D
            if len(s.dict['lwid']) >= igr+1: gr.SetLineWidth(s.dict['lwid'][igr])
            # ----


            

        # 
        if s.VB>1: print(' DEBUG  range min / max = ', zrangemin, '   ', zrangemax)

        
        # Make small margin above/below max/min value
        if not s.dict['zLog']: 
            if s.zrangemin == s.undefined and s.dict['zmin_air']: zrangemin -= (zrangemax-zrangemin) * s.dict['zmin_air']
            if s.zrangemax == s.undefined and s.dict['zmax_air']: zrangemax += (zrangemax-zrangemin) * s.dict['zmax_air']
        else: 
            if s.zrangemin == s.undefined and s.dict['zmin_air_log']: zrangemin /= s.dict['zmin_air_log']
            if s.zrangemax == s.undefined and s.dict['zmax_air_log']: zrangemax *= s.dict['zmax_air_log']


        # safety if the above fails (if there are negative values) [however I think these are already caught, so might never be realised]
        if s.dict['zLog'] and zrangemin <= 0: 
            s.zrangemin = s.dict['zLogMin']  
            print("Warning: logarithmic y-axis and negative values ... ")



        #s.legend = ROOT.TLegend(s.dict['legX1'],s.dict['legY1'], s.dict['legX2'],s.dict['legY2'])


        # Plot
        hDef.SetMinimum(zrangemin)
        hDef.SetMaximum(zrangemax)
        hDef.Draw(s.dict['drawstyle_hdef'])
        if s.dict['tickyright']: ROOT.gPad.SetTicky(s.dict['tickyright'])
        if s.dict['ticktop']: ROOT.gPad.SetTickx(s.dict['ticktop'])

        if s.VB>1: print('Minimum value: ', hDef.GetMinimum())
        if s.VB>1: print('Maximum value: ', hDef.GetMaximum())

        for igr,gr in enumerate(gr_reordered): 
            gr.Draw(s.dict['drawstyle_gr1D'])  # 1D
            


        # Labels/Legends
        optD = {'sty':'lp', 'y2':s.dict['label_y2'], 'DY':s.dict['label_DY'], 'x1':s.dict['label_x1'], 'xwid':s.dict['label_xwid'], 'maxchar':s.dict['label_maxchar'], 'maxcharoverflow':s.dict['label_maxcharoverflow']}
        if s.dict['legend2right']  >= 0: optD['right']  = s.dict['legend2right']
        if s.dict['legend2bottom'] >= 0: optD['bottom'] = s.dict['legend2bottom']
        if s.VB>3: 
            print('s.gr: ', s.gr)
            print('s.resvarsT: ', s.resvarsT) #labelz
            print('optD: ', optD)

        res_legend = BLegends(ROOT.TLegend, objarr=gr_reordered,  titarr=resvarsT_reordered, optD=optD)
        # 'xwid':0.12, 'DY':.., 'y2':..
        
        if s.VB>2: print('res_legend:', res_legend)
        for legend in res_legend['legends']: legend.Draw('L')  # could be several

        # Logarithmic?, grid?
        if s.dict['xLog']: s.c1.SetLogy(s.dict['xLog'])
        if s.dict['zLog']: s.c1.SetLogy(s.dict['zLog'])
        s.c1.SetGrid(s.dict['gridx'],s.dict['gridy'])


        # Text on the plot
        s.plottexts = TextsFromRaw(s.plottexts_raw, ROOT, VB=s.VB) 
        for plottext in s.plottexts: plottext.Draw()

        # Save
        s.c1.SaveAs(s.dict['fn_save'])
        #s.c1.SaveAs(s.dict['fn_save'].replace('.pdf','.C')
        
        if s.dict['commonpdf'] and s.dict['fn_save'].endswith('.pdf'): 
            os.system('cp -p %s ~/munch.pdf' %(s.dict['fn_save'])) # keeps a "central" copy of the last plot (nice for plotting over ssh)

 
    # ##########
    def ReadArg(s):  
 
        # ################################### ARGUMENT READING
        Arg = bkgjelstenArgReader.ArgReader(s.argv, VB=0)

        '''
        if Arg.hasget('-alist'):  print 'a string list: ',Arg.list()
        if Arg.hasget('-alisti'): print 'an integer list: ',Arg.listI()
        if Arg.hasget('-alistf'): print 'a float list: ',Arg.listF()
        if Arg.hasget('-x'):  print 'a string: ',Arg.val()
        if Arg.hasget('-xI'): print 'an integer: ',Arg.valI()
        if Arg.hasget('-xF'): print 'a float: ',Arg.valF()
        '''

        # ====


        if Arg.hasget(['-coord']):  # Ex: -coord MU,M_2:M2
            zz = Arg.list()
            for iz in range(len(zz)):
                z = zz[iz]
                w = z.split(':')
                thevar = w[0]
                if len(w)>1: thevarT = w[1]
                else: thevarT = w[0]
                if iz == 0: s.xcoord,s.xcoordT = thevar,thevarT
                if iz == 1: s.ycoord,s.ycoordT = thevar,thevarT
            s.ndim = len(zz)


        if Arg.hasget(['-delim_resvar_comma']):  # a bit long, but ... 
            s.dict['delim_resvar_comma'] = Arg.val()
        if Arg.hasget(['-delim_resvar_colon']):
            s.dict['delim_resvar_colon'] = Arg.val()
            
        if Arg.has(['--delim_resvar_double']):  # a bit long, but ... 
            s.dict['delim_resvar_comma'] = ',,'
            s.dict['delim_resvar_colon'] = '::'

        if Arg.hasget(['-resvars']):  # Ex: -resvars N1,N2:N20,N3,N2-N1:N2mN1,N2/N1
            zz = Arg.list(s.dict['delim_resvar_comma'])
            for z in zz:
                w = z.split(s.dict['delim_resvar_colon'])
                zvar = w[0]
                if len(w)>1: zvarT = w[1]
                else: zvarT = w[0]
                s.resvars.append(zvar)
                s.resvarsT[zvar] = zvarT


        if Arg.hasget(['-cont','-conts']):  # Ex: -cont var:[varT]:val1,val2,val3:col:sty:wid,,
            zarg = Arg.val()
            zz = zarg.split(',,')
            for z in zz:
                w = z.split(':')  # this should at least give 3-6 values
                if len(w) not in [3,4,5,6]: sys.exit("Fatal::munch  non-allowed argument format: -conts %s" %(zarg))
                zvar  = w[0]
                zvarT = w[1]
                if zvarT == '': zvarT = zvar    # will happen with e.g. -conts N1::50,100:2
                zvals = array.array('d')  # need 'd' because SetContour takes 'd', not 'f'
                for z in w[2].split(','): zvals.append(float(z))
                if len(zvals) == 0: sys.exit("Fatal::munch  no contour values specified in  -conts %s" %(zarg))

                # then optional
                zcol, zsty, zwid = 1,1,1
                if len(w) >= 4: zcol = int(w[3])
                if len(w) >= 5: zsty = int(w[4])
                if len(w) >= 6: zwid = int(w[5])
                # print(zvals,zcol)
                s.conts.append(zvar)
                s.contsT[zvar] = zvarT
                s.contsDef[zvar] = {'vals': zvals, 'col':zcol, 'sty':zsty, 'wid':zwid}   # might be an idea to make a class

        if Arg.hasget(['-moredict','-moredicts']):
            thiswiththiss = Arg.list()  # ex -moredict dict1:dict0,dict2:dict0   # adds two more dicts to check (per scenario)
            for thiswiththis in thiswiththiss: 
                w = thiswiththis.split(':')
                replacewiththis = w[0]
                replacethis = w[1]
                s.moredicts[replacewiththis] = replacethis

        # ====


        if Arg.has(['--dict','--showdict']):
            print('DUMPING DEFAULT VALUES IN THE VARIABLE DICTIONARY, s.dict')
            print('   These can be changed by e.g. -dict I,var1,val1:var2,val2:F,var3,val3:var4,val4')
            print('   where I/F denote integer and float. In the example var2 inherits I from var1,')
            print('   while var4 is neither given as I or F and therefore is taken as string.')
            print('   You may need to inspect the code, though, to understand how they are used.')
            print() 
            for key in s.dict:
                print('%-20s  %s' %(key, s.dict[key]))
            sys.exit()

        if Arg.has(['--h','---help','--automan']):
            os.system("cat %s | grep 'if Arg\.has'" %(s.argv[0]))
            print("[ The above shows all implemented options. For full details inspect the code, %s  (Try also '-h') ]" %(s.argv[0]))
            sys.exit()
            
        if Arg.has(['help','-h','--help','-help']):
            s.showHelp()
            sys.exit()

        if Arg.hasget('-vb'):
            s.VB = Arg.valI()
            if s.VB: print('Verbosity level: %i' %(s.VB))

        if Arg.hasget(['-f','-fn_table']):
            s.fn_table = Arg.list()  # allows a list of inputs (in various modes)   # USAGE??
            s.dict['fn_save'] = s.fn_table[0].replace('.txt','.pdf').replace('.dat','.pdf')
            if s.dict['fn_save'] == s.fn_table[0]: s.dict['fn_save'] += '.pdf'

        if Arg.hasget(['-ffn_pickle']):
            s.ffn_pickle = Arg.val()
            s.readfrompickle = 1

        if Arg.has(['--pickle','--readfrompickle']):   # then will be piped in (if s.ffn_pickle == '')
            s.readfrompickle = 1  

        if Arg.has(['--mode2Da','--mode2D','--mode2d']):
            s.ndim = 2
            s.dict['ticklengthx'] = s.dict['ticklengthy'] = -0.008
            s.dict['rightmargin'] = 0.12
            if s.fn_table: s.title = s.fn_table[0].replace('.txt','').replace('.dat','')
            s.dict['numbersonplot'] = 1

        if Arg.hasget('-labelz'):
            s.labelz_byhand = Arg.list() 
        if Arg.hasget('-labelxy'):
            s.labelx,s.labely = Arg.list() 
        if Arg.hasget('-labelx'):
            s.labelx = Arg.val()

        #if Arg.hasget('-ndim'):
        #    s.ndim = Arg.valI()

        if Arg.hasget('-line1is'):
            s.dict['line1is'] = Arg.val()

        if Arg.hasget('-titles'):
            if s.ndim == 1: s.xtitle,s.ztitle = Arg.list()
            if s.ndim == 2: s.xtitle,s.ytitle,s.ztitle = Arg.list()
            
        if Arg.hasget(['-titlex','-xtitle']): s.xtitle = Arg.val()
        if Arg.hasget(['-titley','-ytitle']):
            if s.ndim == 1: s.ztitle = Arg.val()
            else: s.ytitle = Arg.val()
        if Arg.hasget(['-titlez','-ztitle']): s.ztitle = Arg.val()

        # Range: min and max
        if Arg.hasget('-xyrange'):
            s.xrangemin,s.xrangemax = Arg.listF()
            s.yrangemin = s.xrangemin
            s.yrangemax = s.xrangemax
        if Arg.hasget('-xrange'): s.xrangemin,s.xrangemax = Arg.listF()
        if Arg.hasget('-yrange'):
            if s.ndim == 1: s.zrangemin,s.zrangemax = Arg.listF()
            else: s.yrangemin,s.yrangemax = Arg.listF()
        if Arg.hasget('-zrange'): s.zrangemin,s.zrangemax = Arg.listF()

        # Range: one-sided values: min
        if Arg.hasget('-xyrangemin'):
            s.xrangemin = s.yrangemin = Arg.valF()
        if Arg.hasget('-xrangemin'): s.xrangemin = Arg.valF()
        if Arg.hasget('-yrangemin'):
            if s.ndim == 1: s.zrangemin = Arg.valF()
            else: s.yrangemin = Arg.valF()
        if Arg.hasget('-zrangemin'): s.zrangemin = Arg.valF()
        
        # Range: one-sided values: max
        if Arg.hasget('-xyrangemax'):
            s.xrangemax = s.yrangemax = Arg.valF()
        if Arg.hasget('-xrangemax'): s.xrangemax = Arg.valF()
        if Arg.hasget('-yrangemax'):
            if s.ndim == 1: s.zrangemax = Arg.valF()
            else: s.yrangemax = Arg.valF()
        if Arg.hasget('-zrangemax'): s.zrangemax = Arg.valF()
        



        if Arg.hasget(['-fn_save','-save']):
            s.dict['fn_save'] = Arg.val()
            zend = s.dict['fn_save'].split('.').pop()
            #print zend
            if zend not in s.dict['plot_formats']: 
                s.dict['fn_save'] += '.pdf'
            #if not (s.dict['fn_save'].endswith('.pdf') or s.dict['fn_save'].endswith('.eps')): s.dict['fn_save'] += '.pdf'  # Fragile? 

        if Arg.hasget('-autocolstylemode'):
            s.dict['autoColorStyleMode'] = Arg.val()
            s.dict['autoColorStyleUse'] = 1
        if Arg.hasget('-autocolstyle'):
            s.dict['autoColorStyleUse'] = Arg.valI()
        if Arg.has('--autocolstyle'):
            s.dict['autoColorStyleUse'] = 1

        if Arg.hasget('-lcol'): s.dict['lcol'] = Arg.listI()
        if Arg.hasget('-lsty'): s.dict['lsty'] = Arg.listI()
        if Arg.hasget('-lwid'): s.dict['lwid'] = Arg.listI()


        if Arg.hasget(['-righttick','-rightticky']):
            s.dict['tickyright'] = Arg.valI()
        if Arg.hasget(['-toptick','-toptick']): 
            s.dict['ticktop'] = Arg.valI()
        if Arg.has(['--righttick','--rightticky']):
            s.dict['tickyright'] = 1
        if Arg.has(['--toptick','--toptick']): 
            s.dict['ticktop'] = 1
        if Arg.has(['--rightlabel']): 
            s.dict['tickyright'] = 2  # will show both tick and label on the right axis too

        if Arg.hasget('-ticklengthx'): s.dict['ticklengthx'] = Arg.valF()
        if Arg.hasget('-ticklengthy'): s.dict['ticklengthy'] = Arg.valF()
        if Arg.hasget('-ticklengthxy'): s.dict['ticklengthx'] = s.dict['ticklengthy'] = Arg.valF()
    
        if Arg.hasget('-numbers_angle'): s.dict['numbers_angle'] = Arg.valF()
        if Arg.hasget('-numbers_size'): s.dict['numbers_size'] = Arg.valF()    # default = 0.030
        if Arg.hasget('-numbers_colour'): s.dict['numbers_colour'] = Arg.valI()
        if Arg.hasget('-numbers_asc'):
            z = Arg.listF()
            s.dict['numbers_angle']  = z[0]
            s.dict['numbers_size']   = z[1]
            s.dict['numbers_colour'] = int(z[2])

        # Note: these are also obtainable with -dict construction
        if Arg.hasget('-textprecision'):
            s.dict['textprecision'] = Arg.valF()   # can also set with -dict structure
            s.dict['numbersonplot'] = 1
        # the -marker* do not work on the (default) contour plots
        if Arg.hasget('-marker_size'): s.dict['marker_size'] = Arg.valF()
        if Arg.hasget('-marker_style'): s.dict['marker_style'] = Arg.valI()
        if Arg.hasget('-marker_colour'): s.dict['marker_colour'] = Arg.valI()
        if Arg.hasget('-marker_ssc'):
            z = Arg.listF()
            s.dict['marker_style']  = int(z[0])
            s.dict['marker_size']   = z[1]
            s.dict['marker_colour'] = int(z[2])

        if Arg.hasget('-msty'):  s.dict['msty']  = Arg.valI()    # 1D
        if Arg.hasget('-msize'): s.dict['msize'] = Arg.valF()  # 1D
        if Arg.hasget('-mcol'):  s.dict['mcol']  = Arg.valI()    # 1D

        if Arg.hasget('-marker_ssc1D'):  # 1D
            z = Arg.listF()
            s.dict['msty']  = int(z[0])
            s.dict['msize'] = z[1]
            s.dict['mcol']  = int(z[2])

        if Arg.hasget(['-xLog','-logx']): 
            s.dict['xLog'] = Arg.valI()
        if Arg.has(['--xLog','--logx']): 
            s.dict['xLog'] = 1
            
        #if Arg.hasget('-xLogMin'): 
        #    s.dict['xLogMin'] = Arg.valF()
        #    s.dict['xLog'] = 1


        if Arg.hasget(['-yLog','-logy']): 
            s.dict['yLog'] = Arg.valI()
        if Arg.has(['--yLog','--logy']): 
            s.dict['yLog'] = 1
            
        if Arg.hasget('-yLogMin'): 
            s.dict['yLogMin'] = Arg.valF()
            s.dict['yLog'] = 1


        if Arg.hasget(['-zLog','-logz']): 
            s.dict['zLog'] = Arg.valI()
        if Arg.has(['--zLog','--logz']): 
            s.dict['zLog'] = 1
            
        if Arg.hasget('-zLogMin'): 
            s.dict['zLogMin'] = Arg.valF()
            s.dict['zLog'] = 1

        if Arg.hasget(['-setgrid','-grid','-gridxy']):
            s.dict['gridx'],s.dict['gridy'] = Arg.listI()
        if Arg.hasget(['-setgridx','-gridx']):
            s.dict['gridx'] = Arg.valI()
        if Arg.hasget(['-setgridy','-gridy']):
            s.dict['gridy'] = Arg.valI()
        if Arg.has(['--setgrid','--grid','--gridxy']):
            s.dict['gridx'] = s.dict['gridy'] = 1


        if Arg.hasget('-label_x1'): s.dict['label_x1'] = Arg.valF()
        if Arg.hasget('-label_xwid'): s.dict['label_xwid'] = Arg.valF()
        if Arg.hasget('-label_y2'): s.dict['label_y2'] = Arg.valF()
        if Arg.hasget('-label_DY'): s.dict['label_DY'] = Arg.valF()
        if Arg.hasget('-label_maxchar'): s.dict['label_maxchar'] = Arg.valI()
        if Arg.hasget('-label_maxcharoverflow'): s.dict['label_maxcharoverflow'] = Arg.valI()


        if Arg.hasget(['-title','-tit','-ctitle','-ctit']): s.title = Arg.val()

        if Arg.hasget('-text'):  # <text>,x,y[,size[D=],color[D=1],font[D=]]:<text>,...
            s.plottexts_raw = Arg.val().split(':')
            #print s.plottexts_raw

        if Arg.hasget(['-style2D']):
            s.drawstyle_gr2D = Arg.val()
        if Arg.hasget(['-numbers','-textnumbers','-numbersonplot']):
            s.dict['numbersonplot'] = Arg.valI()
        if Arg.has(['--numbers','--textnumbers','--numbersonplot']):
            s.dict['numbersonplot'] = 1
            
        if Arg.hasget(['-topmargin','-top']): s.dict['topmargin'] = Arg.valF()
        if Arg.hasget(['-bottommargin','-bottom']): s.dict['bottommargin'] = Arg.valF()
        if Arg.hasget(['-leftmargin','-left']): s.dict['leftmargin'] = Arg.valF()
        if Arg.hasget(['-rightmargin','-right']): s.dict['rightmargin'] = Arg.valF()

        if Arg.has('--legend2bottom'): 
            s.dict['legend2bottom'] = s.dict['bottommargin'] + 0.02
        if Arg.has('--legend2right'): 
            s.dict['legend2right'] = 1. - s.dict['rightmargin'] - 0.02
        if Arg.has('--legend2bottomright'): 
            s.dict['legend2bottom'] = s.dict['bottommargin'] + 0.02
            s.dict['legend2right'] = 1. - s.dict['rightmargin'] - 0.02

        if Arg.hasget('-legend2bottom'): 
            s.dict['legend2bottom'] = Arg.valF()
        if Arg.hasget('-legend2right'): 
            s.dict['legend2right'] = Arg.valF()

        if Arg.has('--largestfirst'):  s.dict['1D_reorder'] = 1
        if Arg.has('--smallestfirst'): s.dict['1D_reorder'] = -1
        if Arg.has(['--noreorder','--keeporder']): s.dict['1D_reorder'] = 0
        if Arg.hasget('-1D_reorder'): s.dict['1D_reorder'] = Arg.valI()


        #if Arg.hasget('-cols'): s.cols = Arg.listI()
        #if Arg.hasget('-cols'): s.cols = Arg.list()  # can take colnumbers or titles ..  # 2014-08-21: is this deprecated / replaced to -coords <..>  -resvars <resvars> ?? Think so. Commenting it out. Changing use of -cols in help
        
        if Arg.hasget('-scale'): s.dict['scale'] = Arg.valF()

        if Arg.hasget(['-protection','-operationprotection']): s.dict['operationprotection'] = Arg.valI()  # 1 if use e.g. [+] to sum in resvalues

        if Arg.has(['--nobatch','--interactive']): 
            s.rootbatch = 0

        if Arg.has(['--c']):
            s.dict['saveCtoo'] = 1

        if Arg.has(['--marker']):
            s.dict['marker'] = 1
        
        if Arg.hasget(['-marker']):
            s.dict['marker'] = 1
            z = Arg.listF()
            if z: s.dict['marker_size']  = z.pop(0)
            if z: s.dict['marker_colour'] = int(z.pop(0))
            if z: s.dict['marker_style'] = int(z.pop(0))
            
        if Arg.hasget(['-numbermarker', '-numbers_marker']):
            s.numbermarker = Arg.val()
            s.dict['numbersonplot'] = 1
        if Arg.has(['--numbermarker', '--numbers_marker']):
            s.numbermarker = 'o' # default
            s.dict['numbersonplot'] = 1


        if Arg.has(['--maketable']): s.dict['tableverbose_resvar'] = 1
        
        # ----- The new general procedure for var input (should this be put into the ArgReader?)
        if Arg.hasget('-dict'):
            zs = Arg.list(':')
            # print zs
            for z in zs:
                zw = z.split(',')
                print(zw)
                
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
                    


        errors = Arg.ErrorMessages()
        if errors: 
            print('Problems...: showing help text')
            s.showHelp()
            print(80*'%','\n',80*'%','\n',80*'%', "FATAL:  ENDING DUE TO PROBLEMS OF ARGUMENTS:")
            for error in errors: print(error)
            sys.exit()
    
        # ################################### POST-INIT

        

############################## EXECUTE IF RUN AS SCRIPT (NOT JUST IMPORTED)
#if __name__ == '__main__':
t = munch(cmd=['ReadArg','PostInit','PrepareXYZ','Main'])
############################## 


