 #!/usr/bin/env python

import sys
import argparse

from array import array

import ROOT
from ROOT import *
from ROOT import TGraph
ROOT.gROOT.SetBatch(1) # Don't show graphics
ROOT.gROOT.LoadMacro( '~/atlasstyle/AtlasStyle.C' )
ROOT.SetAtlasStyle()
ROOT.gROOT.LoadMacro( '~/atlasstyle/AtlasLabels.C' )

import contourPlotter


# ===== Command line arguments ===== #

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--typee',              type=str,                        help='Type: original, new, both',                                            required=True)
parser.add_argument('-w', '--which',              type=str, default='both',        help='Which contours to plot. Default: both. Alternatives: expected, observed', required=False)
parser.add_argument('-b', '--backend',            type=str, default='jax',         help='Backend. Default: jax. Alternatives: numpy, pytorch, tensorflow',      required=False)
parser.add_argument('-o', '--optimizer',          type=str, default='scipy',       help='Optimizer. Default: scipy. Alternative: minuit',                       required=False)
parser.add_argument('-v', '--variation',          type=str, default='Nominal',     help='SigXSec variation. Default: Nominal. Alternatives: Up, Down',          required=False)
parser.add_argument('-a', '--addValues',          type=str, default='none,-',      help='Print CLs values in plot. Default: none"-. Alternative: clsobs/clsexp,grid, where grid can be 3L-onshell, allhad, 1Lbb, SS3L, 1L or 2tau', required=False)
parser.add_argument('-l', '--label',              type=str, default='none',        help='Optional filename label',                                              required=False)
parser.add_argument('-cl', '--confidenceLevel',   type=str, default='95',          help='Contour Confidence Level (CL). Default 95%%CL',                        required=False)
parser.add_argument('-xmax', '--xMax',            type=str, default='1250',        help='X axis maximum. Default: 1250',                                        required=False)
parser.add_argument('-O','--onlyOriginalGrid', action='store_true', help='Use only original workspaces for 3L-onshell, allhad and 1Lbb', required=False)


args = parser.parse_args()

if args.typee in ['original', 'new', 'both']:
    typee = args.typee
else:
    print(f'Not recognizing type {args.typee}')
    sys.exit(1)

if args.which in ['both', 'expected', 'observed']:
    which = args.which
else:
    print(f'Not recognizing "which" {args.which}')
    sys.exit(1)
    
if args.backend in ['pytorch', 'numpy', 'jax', 'tensorflow']:
    backend = args.backend
else:
    print(f'Not recognizing backend {args.backend}')
    sys.exit(1)
    
if args.optimizer in ['scipy', 'minuit']:
    optimizer = args.optimizer
else:
    print(f'Not recognizing optimizer {args.optimizer}')
    sys.exit(1)

if args.variation in ['Nominal', 'Up', 'Down']:
    variation = args.variation
else:
    print(f'Not recognizing variation {args.variation}')
    sys.exit(1)

addValues_list = [item for item in args.addValues.split(',')]
if addValues_list[0] in ['none', 'clsobs', 'clsexp'] and addValues_list[1] in ['-', '3L-onshell', 'allhad', '1Lbb', '1L', 'SS3L', '2tau']:
    addValues = addValues_list
else:
    print(f'Not recognizing addValues {args.addValues}')
    sys.exit(1)

label = args.label

if args.confidenceLevel in ['85','90','95']:
    level = args.confidenceLevel
else:
    print(f'Not recognizing CL {args.confidenceLevel}%%')
    sys.exit(1)

xmax_in = args.xMax

if args.onlyOriginalGrid:
    onlyOriginalGrid = True
else:
    onlyOriginalGrid = False
    
print('\n======================================')
print('         Plotting contours            ')
print('======================================')
print(f'  type               = {typee}')
print(f'  which              = {which}')
print(f'  backend            = {backend}')
print(f'  optimizer          = {optimizer}')
print(f'  SigXSec variation  = {variation}')
print(f'  addValues          = {addValues}')
print(f'  Confidence Level   = {level}')
print(f'  xmax               = {xmax_in}')
print(f'  filename label     = {label}')
print(f'  only original workspaces for 3L, allhad and 1Lbb = {onlyOriginalGrid}')
print('======================================\n')

if (which == 'observed' and addValues[0] == 'clsexp') or (which == 'expected' and addValues[0] == 'clsobs'):
    print(f'Does not make sense to print CLs values for {addValues[0]} when only plotting {which} contours')
    sys.exit(1)


# ===== Files and objects  ===== #

main_folder = '/home/elibry/comb/winobino/results'
folder_3Lonshell = f'{main_folder}/3L-onshell'
folder_allhad    = f'{main_folder}/allhad'
folder_1Lbb      = f'{main_folder}/1Lbb'
folder_1L        = f'{main_folder}/1L'
folder_SS3L      = f'{main_folder}/SS3L'
folder_2tau      = f'{main_folder}/2tau'

fname_1L             = f'{folder_1L}/1L_{backend}_{optimizer}_fixSigXSec{variation}'
fname_SS3L           = f'{folder_SS3L}/SS3L_{backend}_{optimizer}_fixSigXSec{variation}'
fname_2tau           = f'{folder_2tau}/2tau_{backend}_{optimizer}_fixSigXSec{variation}'
if onlyOriginalGrid:
      fname_3Lonshell      = f'{folder_3Lonshell}/3L-onshell_OriginalGridOnly_{backend}_{optimizer}_fixSigXSec{variation}'
      fname_allhad         = f'{folder_allhad}/allhad_OriginalGridOnly_{backend}_{optimizer}_fixSigXSec{variation}'
      fname_1Lbb           = f'{folder_1Lbb}/1Lbb_OriginalGridOnly_{backend}_{optimizer}_fixSigXSec{variation}'
else:
      fname_3Lonshell      = f'{folder_3Lonshell}/3L-onshell_{backend}_{optimizer}_fixSigXSec{variation}'
      fname_allhad         = f'{folder_allhad}/allhad_{backend}_{optimizer}_fixSigXSec{variation}'
      fname_1Lbb           = f'{folder_1Lbb}/1Lbb_{backend}_{optimizer}_fixSigXSec{variation}'

fname_orig_3Lonshell = f'{folder_3Lonshell}/original/3L-onshell_original_fixSigXSec{variation}'
fname_orig_allhad    = f'{folder_allhad}/original/allhad_original_fixSigXSec{variation}'
#fname_orig_1Lbb      = f'{folder_1Lbb}/original/1Lbb_original_fixSigXSec{variation}'
fname_orig_1Lbb      = f'{folder_1Lbb}/original/1Lbb_graphs.root' # Directly from analysis team
fname_orig_1L        = f'{folder_1L}/original/1L_original_fixSigXSec{variation}'
fname_orig_SS3L      = f'{folder_SS3L}/original/SS3L_original_fixSigXSec{variation}'
fname_orig_2tau      = f'{folder_2tau}/original/2tau_original_fixSigXSec{variation}'

if level != '95':
    fname_ending = f'_CL{level}_contours.root'
else:
    fname_ending = '_contours.root'
fname_ending_orig = '_contours.root'
    
f_3Lonshell      = TFile(fname_3Lonshell      + fname_ending)
f_orig_3Lonshell = TFile(fname_orig_3Lonshell + fname_ending_orig)
f_allhad         = TFile(fname_allhad         + fname_ending)
f_orig_allhad    = TFile(fname_orig_allhad    + fname_ending_orig)
f_1Lbb           = TFile(fname_1Lbb           + fname_ending)
f_orig_1Lbb      = TFile(fname_orig_1Lbb) #      + fname_ending) # Directly from analysis team, different file name
# f_1L             = TFile(fname_1L             + fname_ending)
# f_orig_1L        = TFile(fname_orig_1L        + fname_ending)
f_SS3L           = TFile(fname_SS3L           + fname_ending) 
f_orig_SS3L      = TFile(fname_orig_SS3L      + fname_ending_orig)
f_2tau           = TFile(fname_2tau           + fname_ending)
f_orig_2tau      = TFile(fname_orig_2tau      + fname_ending_orig)

### Contours
# Expected
exp_3Lonshell      = f_3Lonshell.Get('Exp_0')
exp_orig_3Lonshell = f_orig_3Lonshell.Get('Exp_0')
exp_allhad         = f_allhad.Get('Exp_0')
exp_orig_allhad    = f_orig_allhad.Get('Exp_0')
exp_1Lbb           = f_1Lbb.Get('Exp_0')
exp_orig_1Lbb      = f_orig_1Lbb.Get('Exp_0')
# exp_1L             = f_1L.Get('Exp_0')
# exp_orig_1L        = f_orig_1L.Get('Exp_0')
exp_SS3L           = f_SS3L.Get('Exp_0')
exp_orig_SS3L      = f_orig_SS3L.Get('Exp_0')
exp_2tau           = f_2tau.Get('Exp_0')
exp_orig_2tau      = f_orig_2tau.Get('Exp_0')

# Observed
obs_3Lonshell      = f_3Lonshell.Get('Obs_0')
obs_orig_3Lonshell = f_orig_3Lonshell.Get('Obs_0')
obs_allhad         = f_allhad.Get('Obs_0')
obs_orig_allhad    = f_orig_allhad.Get('Obs_0')
obs_1Lbb           = f_1Lbb.Get('Obs_0')
obs_orig_1Lbb      = f_orig_1Lbb.Get('Obs_0')
# obs_1L             = f_1L.Get('Obs_0')
# obs_orig_1L        = f_orig_1L.Get('Obs_0')
#obs_SS3L           = f_SS3L.Get('Obs_0') 
obs_SS3L           = f_SS3L.Get('Exp_0') # Use exp for obs when blinded
#obs_orig_SS3L      = f_orig_SS3L.Get('Obs_0')
obs_orig_SS3L      = f_orig_SS3L.Get('Exp_0') # Use exp for obs when blinded
obs_2tau           = f_2tau.Get('Obs_0')
obs_orig_2tau      = f_orig_2tau.Get('Obs_0')

### CLs values
if addValues[0] == 'clsobs':
    cls_3Lonshell      = f_3Lonshell.Get('CLs_gr')
    cls_allhad         = f_allhad.Get('CLs_gr')
    cls_1Lbb           = f_1Lbb.Get('CLs_gr')
    # cls_1L             = f_1L.Get('CLs_gr')
    cls_SS3L           = f_SS3L.Get('CLs_gr')
    cls_2tau           = f_2tau.Get('CLs_gr')
    cls_orig_3Lonshell = f_orig_3Lonshell.Get('CLs_gr')
    cls_orig_allhad    = f_orig_allhad.Get('CLs_gr')
    cls_orig_1Lbb      = f_orig_1Lbb.Get('CLs_gr')
    # cls_orig_1L        = f_orig_1L.Get('CLs_gr')
    cls_orig_SS3L      = f_orig_SS3L.Get('CLs_gr')
    cls_orig_2tau      = f_orig_2tau.Get('CLs_gr')
elif addValues[0] == 'clsexp':
    cls_3Lonshell      = f_3Lonshell.Get('CLsexp_gr')
    cls_allhad         = f_allhad.Get('CLsexp_gr')
    cls_1Lbb           = f_1Lbb.Get('CLsexp_gr')
    # cls_1L             = f_1L.Get('CLsexp_gr')
    cls_SS3L           = f_SS3L.Get('CLsexp_gr')
    cls_2tau           = f_2tau.Get('CLsexp_gr')
    cls_orig_3Lonshell = f_orig_3Lonshell.Get('CLsexp_gr')
    cls_orig_allhad    = f_orig_allhad.Get('CLsexp_gr')
    cls_orig_1Lbb      = f_orig_1Lbb.Get('CLsexp_gr')
    # cls_orig_1L        = f_orig_1L.Get('CLsexp_gr')
    cls_orig_SS3L      = f_orig_SS3L.Get('CLsexp_gr')
    cls_orig_2tau      = f_orig_2tau.Get('CLsexp_gr')
elif addValues[0] == 'none':
    pass


# ===== Plot settings ===== #

# Output filename
if typee == 'original':
    file_name = f'contour_{typee}_{variation}'
else:
    file_name = f'contour_{typee}_{backend}_{optimizer}_{variation}'
if which != 'both':
    file_name += f'_{which}Only'
if addValues[0] != 'none':
    file_name += f'_{addValues[0]}{addValues[1]}'
if label != 'none':
    file_name += f'_{label}'
if level != '95':
    file_name += f'_CL{level}'
if onlyOriginalGrid:
    file_name += '_onlyOriginalGridFor3LAH1Lbb'

# Process
process = '#tilde{#chi}^{0}_{2}#tilde{#chi}^{#pm}_{1}#rightarrowWh#tilde{#chi}^{0}_{1}#tilde{#chi}^{0}_{1}'

# Lumi label
lumilabel = '#sqrt{s}=13 TeV, 139 fb^{-1}'

# Axis labels
xlabel = 'm(#tilde{#chi}^{#pm}_{1}/#tilde{#chi}^{0}_{2}) [GeV]'
ylabel = 'm(#tilde{#chi}^{0}_{1}) [GeV]'

# Analysis legends
leg_1Lbb   = '1Lbb'
leg_3L     = '3L'
leg_allhad = 'all-hadronic'
leg_SS3L   = 'SS/3L'
leg_2tau   = '2tau'
leg_1L     = '1L'

# Axis ranges
xmin = 100 
ymin = 0
if addValues[0] != 'none':
    xmax = int(xmax_in)
    ymax = 550
else:
    xmax = int(xmax_in)
    ymax = 550

# Colors
blue       = TColor.GetColor('#68D3FB')
orange     = TColor.GetColor('#FEB144')
purple     = TColor.GetColor('#CC99C9')
yellow     = TColor.GetColor('#FEE227')
red        = TColor.GetColor('#FF6663')
green      = TColor.GetColor('#9EE09E')
pink       = TColor.GetColor('#FF33B8')
darkblue   = TColor.GetColor('#3368FF')
darkgreen  = TColor.GetColor('#00cc44') 
darkred    = TColor.GetColor('#ff3333')
darkpurple = TColor.GetColor('#a64dff') 
darkorange = TColor.GetColor('#cc7a00') 

orig_3Lonshell_color = darkblue
orig_allhad_color    = darkorange
orig_1Lbb_color      = darkpurple
#orig_1L_color        = pink
orig_SS3L_color      = darkgreen
orig_2tau_color      = darkred

new_3Lonshell_color = blue
new_allhad_color    = orange
new_1Lbb_color      = purple
#new_1L_color        = yellow
new_SS3L_color      = green
new_2tau_color      = red

# CLs values title
if addValues[0] == 'clsobs':
    numbers_title = 'Numbers represent observed CLs value'
elif addValues[0] == 'clsexp':
    numbers_title = 'Numbers represent expected CLs value'
elif addValues[0] == 'none':
    numbers_title = ''

# CLs values offset
cls_xoffset = 6
cls_yoffset = 0.1
    
# Text sizes
legTextSize = 0.02
valuesTextSize = 0.011 #0.013

# Line widths
linew_exp = 2

# ===== Plot ===== #

plots_folder = '/home/elibry/comb/winobino/scripts/plots/singleContours'

# General
plot = contourPlotter.contourPlotter(plots_folder+'/'+file_name, 800, 600)
plot.processLabel = process
plot.lumiLabel = lumilabel
plot.drawAxes([xmin, ymin, xmax, ymax])
plot.decorateCanvas()
plot.setXAxisLabel(xlabel)
plot.setYAxisLabel(ylabel)

# 'Expected' and 'Observed' legends
h1 = TH1F('observed','',1,0,1)
h2 = TH1F('expected','',1,0,1)
h2.SetLineStyle(2)
h1.SetLineColor(kBlack); h2.SetLineColor(kBlack)
h1.SetLineWidth(3); h2.SetLineWidth(1)
leg2 = TLegend(*(0.17,0.85,0.32,0.91)) # xmin, ymin, xmax, ymax
leg2.AddEntry(h1,'Observed limits','l')
leg2.AddEntry(h2,'Expected limits','l')
leg2.SetTextFont(42)
leg2.SetBorderSize(0)
#leg2.SetFillStyle(0)
leg2.SetTextSize(legTextSize)
leg2.SetMargin(0.5)

# 'Limits at X$ CL' text
l = TLatex()
l.SetTextAlign(12)
l.SetTextSize(legTextSize)
l.SetNDC()
l.SetTextFont(42)
l.SetTextColor(kBlack)
l.DrawLatex(0.245, 0.836, f'All limits at {level}% CL') # x, y, text 

# ATLAS label
ROOT.ATLASLabel( 0.4, 0.945, "Internal" ) # Defined in atlasstyle macro

# Legend shapes
#leg_shape = (0.17, 0.70, 0.47, 0.77) # 2 lines
#leg_shape = (0.17, 0.67, 0.47, 0.77) # 3 lines
#leg_shape = (0.17, 0.64, 0.47, 0.77) # 4 lines
leg_shape_5 = (0.17, 0.61, 0.47, 0.77) # 5 lines
leg_shape_10 = (0.17, 0.47, 0.47, 0.80) #10 lines

# Draw contours and legends
if typee == 'original':

    # In order to get legends also when plotting only expected contours
    if which == 'expected':
        plot.drawExpected(exp_orig_3Lonshell, color=orig_3Lonshell_color, title=f'{leg_3L} original (HF)',     legendOrder=0, linew=linew_exp)
        plot.drawExpected(exp_orig_allhad,    color=orig_allhad_color,    title=f'{leg_allhad} original (HF)', legendOrder=1, linew=linew_exp)
        plot.drawExpected(exp_orig_1Lbb,      color=orig_1Lbb_color,      title=f'{leg_1Lbb} original (HF)',   legendOrder=2, linew=linew_exp)
        plot.drawExpected(exp_orig_SS3L,      color=orig_SS3L_color,      title=f'{leg_SS3L} original (HF)',   legendOrder=3, linew=linew_exp) 
        #plot.drawExpected(exp_orig_1L,        color=orig_1L_color,         title=f'{leg_1L} original (HF)',     legendOrder=4, linew=linew_exp)
        plot.drawExpected(exp_orig_2tau,      color=orig_2tau_color,      title=f'{leg_2tau} original (HF)',   legendOrder=5, linew=linew_exp)

    if which == 'both': 
        plot.drawExpected(exp_orig_3Lonshell, color=orig_3Lonshell_color)
        plot.drawExpected(exp_orig_allhad,    color=orig_allhad_color)
        plot.drawExpected(exp_orig_1Lbb,      color=orig_1Lbb_color)
        plot.drawExpected(exp_orig_SS3L,      color=orig_SS3L_color) 
        #plot.drawExpected(exp_orig_1L,        color=orig_1L_color)
        plot.drawExpected(exp_orig_2tau,      color=orig_2tau_color)
        
    if which in ['both', 'observed']:
        plot.drawObserved(obs_orig_3Lonshell, color=orig_3Lonshell_color, title=f'{leg_3L} original (HF)',         legendOrder=0)
        plot.drawObserved(obs_orig_allhad,    color=orig_allhad_color,    title=f'{leg_allhad} original (HF)',     legendOrder=1)
        plot.drawObserved(obs_orig_1Lbb,      color=orig_1Lbb_color,      title=f'{leg_1Lbb} original (HF)',       legendOrder=2)
        plot.drawObserved(obs_orig_SS3L,      color=orig_SS3L_color,      title=f'{leg_SS3L} original (HF)',       legendOrder=3) 
        #plot.drawObserved(obs_orig_1L,        color=orig_1L_color,        title=f'{leg_1L} original (HF)',        legendOrder=4)
        plot.drawObserved(obs_orig_2tau,      color=orig_2tau_color,      title=f'{leg_2tau} original (HF)',      legendOrder=5)

    leg = plot.createLegend(shape=leg_shape_5) 
    leg.SetTextSize(legTextSize)
    leg.Draw('same')

    # Add legends for expected and observed
    leg2.Draw('same')
    
    # Add Cls values to plot 
    if not addValues[0] == 'none': 

        ### 3L-onshell
        if addValues[1] == '3L-onshell':
            # Values
            plot.drawTextFromTGraph2D(cls_orig_3Lonshell, size=valuesTextSize, title=numbers_title, color=orig_3Lonshell_color, xoffset=cls_xoffset, yoffset=cls_yoffset, borders=[xmin, ymin, xmax, ymax])       
            # Markers
            x_mark = []; y_mark = []
            x,y,z  = cls_orig_3Lonshell.GetX(), cls_orig_3Lonshell.GetY(), cls_orig_3Lonshell.GetZ() 
            for i in range(0,cls_orig_3Lonshell.GetN()):
                if(z[i]<0): continue
                x_mark.append(x[i])
                y_mark.append(y[i])
            bullets1 = TGraph(len(x_mark),array('d',x_mark),array('d',y_mark))
            bullets1.SetMarkerStyle(8)
            bullets1.SetMarkerSize(0.4)
            bullets1.Draw("PSAME")
        
            ### allhad
        elif addValues[1] == 'allhad':    
            # Values
            plot.drawTextFromTGraph2D(cls_orig_allhad, size=valuesTextSize, title=numbers_title, color=orig_allhad_color, xoffset=cls_xoffset, yoffset=cls_yoffset, borders=[xmin, ymin, xmax, ymax])     
            # Markers
            x_mark = []; y_mark = []
            x,y,z  = cls_orig_allhad.GetX(), cls_orig_allhad.GetY(), cls_orig_allhad.GetZ() 
            for i in range(0,cls_orig_allhad.GetN()):
                if(z[i]<0): continue
                x_mark.append(x[i])
                y_mark.append(y[i])
            bullets2 = TGraph(len(x_mark),array('d',x_mark),array('d',y_mark))
            bullets2.SetMarkerStyle(8)
            bullets2.SetMarkerSize(0.4)
            bullets2.Draw("PSAME")

        ### 1Lbb
        elif addValues[1] == '1Lbb':    
            # Values
            plot.drawTextFromTGraph2D(cls_orig_1Lbb, size=valuesTextSize, title=numbers_title, color=orig_1Lbb_color, xoffset=cls_xoffset, yoffset=cls_yoffset, borders=[xmin, ymin, xmax, ymax])
            # Markers
            x_mark = []; y_mark = []
            x,y,z  = cls_orig_1Lbb.GetX(), cls_orig_1Lbb.GetY(), cls_orig_1Lbb.GetZ() 
            for i in range(0,cls_orig_1Lbb.GetN()):
                if(z[i]<0): continue
                x_mark.append(x[i])
                y_mark.append(y[i])
            bullets3 = TGraph(len(x_mark),array('d',x_mark),array('d',y_mark))
            bullets3.SetMarkerStyle(8)
            bullets3.SetMarkerSize(0.4)
            bullets3.Draw("PSAME")

        ### SS3L
        elif addValues[1] == 'SS3L':
            # Values
            plot.drawTextFromTGraph2D(cls_orig_SS3L, size=valuesTextSize, title=numbers_title, color=orig_SS3L_color, xoffset=cls_xoffset, yoffset=cls_yoffset, borders=[xmin, ymin, xmax, ymax])
            # Markers
            x_mark = []; y_mark = []
            x,y,z  = cls_orig_SS3L.GetX(), cls_orig_SS3L.GetY(), cls_orig_SS3L.GetZ() 
            for i in range(0,cls_orig_SS3L.GetN()):
                if(z[i]<0): continue
                x_mark.append(x[i])
                y_mark.append(y[i])
            bullets4 = TGraph(len(x_mark),array('d',x_mark),array('d',y_mark))
            bullets4.SetMarkerStyle(8)
            bullets4.SetMarkerSize(0.4)
            bullets4.Draw("PSAME")

        ### 1L
        elif addValues[1] == '1L':    
            # Values
            plot.drawTextFromTGraph2D(cls_orig_1L, size=valuesTextSize, title=numbers_title, color=orig_1L_color, xoffset=cls_xoffset, yoffset=cls_yoffset, borders=[xmin, ymin, xmax, ymax])
            # Markers
            x_mark = []; y_mark = []
            x,y,z  = cls_orig_1L.GetX(), cls_orig_1L.GetY(), cls_orig_1L.GetZ() 
            for i in range(0,cls_orig_1L.GetN()):
                if(z[i]<0): continue
                x_mark.append(x[i])
                y_mark.append(y[i])
            bullets5 = TGraph(len(x_mark),array('d',x_mark),array('d',y_mark))
            bullets5.SetMarkerStyle(8)
            bullets5.SetMarkerSize(0.4)
            bullets5.Draw("PSAME")

        ### 2tau
        elif addValues[1] == '2tau':
            # Values
            plot.drawTextFromTGraph2D(cls_orig_2tau, size=valuesTextSize, title=numbers_title, color=orig_2tau_color, xoffset=cls_xoffset, yoffset=cls_yoffset, borders=[xmin, ymin, xmax, ymax])
            # Markers
            x_mark = []; y_mark = []
            x,y,z  = cls_orig_2tau.GetX(), cls_orig_2tau.GetY(), cls_orig_2tau.GetZ() 
            for i in range(0,cls_orig_2tau.GetN()):
                if(z[i]<0): continue
                x_mark.append(x[i])
                y_mark.append(y[i])
            bullets6 = TGraph(len(x_mark),array('d',x_mark),array('d',y_mark))
            bullets6.SetMarkerStyle(8)
            bullets6.SetMarkerSize(0.4)
            bullets6.Draw("PSAME")
                        
elif typee == 'new':

    if which == 'expected': 
        plot.drawExpected(exp_3Lonshell, color=new_3Lonshell_color, title=f'{leg_3L} (pyhf)',     legendOrder = 0, linew=linew_exp)
        plot.drawExpected(exp_allhad,    color=new_allhad_color,    title=f'{leg_allhad} (pyhf)', legendOrder = 1, linew=linew_exp)
        plot.drawExpected(exp_1Lbb,      color=new_1Lbb_color,      title=f'{leg_1Lbb} (pyhf)',   legendOrder = 2, linew=linew_exp)
        plot.drawExpected(exp_SS3L,      color=new_SS3L_color,      title=f'{leg_SS3L} (pyhf)',   legendOrder = 3, linew=linew_exp)
        #plot.drawExpected(exp_1L,        color=new_1L_color,        title=f'{leg_1L} (pyhf)',     legendOrder = 4, linew=linew_exp)
        plot.drawExpected(exp_2tau,      color=new_2tau_color,      title=f'{leg_2tau} (pyhf)',   legendOrder = 5, linew=linew_exp)

    if which == 'both': 
        plot.drawExpected(exp_3Lonshell, color=new_3Lonshell_color)
        plot.drawExpected(exp_allhad,    color=new_allhad_color)
        plot.drawExpected(exp_1Lbb,      color=new_1Lbb_color)
        plot.drawExpected(exp_SS3L,      color=new_SS3L_color)
        #plot.drawExpected(exp_1L,        color=new_1L_color)
        plot.drawExpected(exp_2tau,      color=new_2tau_color)

    if which in ['both', 'observed']:
        plot.drawObserved(obs_3Lonshell, color=new_3Lonshell_color, title=f'{leg_3L} (pyhf)', legendOrder=0) 
        plot.drawObserved(obs_allhad,    color=new_allhad_color,    title=f'{leg_allhad} (pyhf)',       legendOrder=1)
        plot.drawObserved(obs_1Lbb,      color=new_1Lbb_color,      title=f'{leg_1Lbb} (pyhf)',       legendOrder=2) 
        plot.drawObserved(obs_SS3L,      color=new_SS3L_color,      title=f'{leg_SS3L} (pyhf)',       legendOrder=3) 
        #plot.drawObserved(obs_1L,        color=new_1L_color,        title=f'{leg_1L} (pyhf)',         legendOrder=4)
        plot.drawObserved(obs_2tau,      color=new_2tau_color,      title=f'{leg_2tau} (pyhf)',       legendOrder=5) 

    leg = plot.createLegend(shape=leg_shape_5) 
    leg.SetTextSize(legTextSize)
    leg.Draw('same')

    # Add legends for expected and observed
    leg2.Draw('same')

    # Add values to plot
    if not addValues[0] == 'none': 

        ### 3L-onshell
        if addValues[1] == '3L-onshell':
            # Values
            plot.drawTextFromTGraph2D(cls_3Lonshell, size=valuesTextSize, title=numbers_title, color=new_3Lonshell_color, xoffset=cls_xoffset, yoffset=cls_yoffset, borders=[xmin, ymin, xmax, ymax])
            # Markers
            x_mark = []; y_mark = []
            x,y,z  = cls_3Lonshell.GetX(), cls_3Lonshell.GetY(), cls_3Lonshell.GetZ() 
            for i in range(0,cls_3Lonshell.GetN()):
                if(z[i]<0): continue
                x_mark.append(x[i])
                y_mark.append(y[i])
            bullets1 = TGraph(len(x_mark),array('d',x_mark),array('d',y_mark))
            bullets1.SetMarkerStyle(8)
            bullets1.SetMarkerSize(0.4)
            bullets1.Draw("PSAME")
        
        ### allhad
        elif addValues[1] == 'allhad':
            # Values
            plot.drawTextFromTGraph2D(cls_allhad, size=valuesTextSize, title=numbers_title, color=new_allhad_color, xoffset=cls_xoffset, yoffset=cls_yoffset, borders=[xmin, ymin, xmax, ymax])
            # Markers
            x_mark = []; y_mark = []
            x,y,z  = cls_allhad.GetX(), cls_allhad.GetY(), cls_allhad.GetZ() 
            for i in range(0,cls_allhad.GetN()):
                if(z[i]<0): continue
                x_mark.append(x[i])
                y_mark.append(y[i])
            bullets2 = TGraph(len(x_mark),array('d',x_mark),array('d',y_mark))
            bullets2.SetMarkerStyle(8)
            bullets2.SetMarkerSize(0.4)
            bullets2.Draw("PSAME")

        ### 1Lbb
        elif addValues[1] == '1Lbb':
            # Values
            plot.drawTextFromTGraph2D(cls_1Lbb, size=valuesTextSize, title=numbers_title, color=new_1Lbb_color, xoffset=cls_xoffset, yoffset=cls_yoffset, borders=[xmin, ymin, xmax, ymax])
            # Markers
            x_mark = []; y_mark = []
            x,y,z  = cls_1Lbb.GetX(), cls_1Lbb.GetY(), cls_1Lbb.GetZ() 
            for i in range(0,cls_1Lbb.GetN()):
                if(z[i]<0): continue
                x_mark.append(x[i])
                y_mark.append(y[i])
            bullets3 = TGraph(len(x_mark),array('d',x_mark),array('d',y_mark))
            bullets3.SetMarkerStyle(8)
            bullets3.SetMarkerSize(0.4)
            bullets3.Draw("PSAME")
        
        ### SS3L
        elif addValues[1] == 'SS3L':
            # Values
            plot.drawTextFromTGraph2D(cls_SS3L, size=valuesTextSize, title=numbers_title, color=new_SS3L_color, xoffset=cls_xoffset, yoffset=cls_yoffset, borders=[xmin, ymin, xmax, ymax])
            # Markers
            x_mark = []; y_mark = []
            x,y,z  = cls_SS3L.GetX(), cls_SS3L.GetY(), cls_SS3L.GetZ() 
            for i in range(0,cls_SS3L.GetN()):
                if(z[i]<0): continue
                x_mark.append(x[i])
                y_mark.append(y[i])
            bullets4 = TGraph(len(x_mark),array('d',x_mark),array('d',y_mark))
            bullets4.SetMarkerStyle(8)
            bullets4.SetMarkerSize(0.4)
            bullets4.Draw("PSAME")

        ### 1L
        elif addValues[1] == '1L':
            # Values
            plot.drawTextFromTGraph2D(cls_1L, size=valuesTextSize, title=numbers_title, color=new_1L_color, xoffset=cls_xoffset, yoffset=cls_yoffset, borders=[xmin, ymin, xmax, ymax])
            # Markers
            x_mark = []; y_mark = []
            x,y,z  = cls_1L.GetX(), cls_1L.GetY(), cls_1L.GetZ() 
            for i in range(0,cls_1L.GetN()):
                if(z[i]<0): continue
                x_mark.append(x[i])
                y_mark.append(y[i])
            bullets5 = TGraph(len(x_mark),array('d',x_mark),array('d',y_mark))
            bullets5.SetMarkerStyle(8)
            bullets5.SetMarkerSize(0.4)
            bullets5.Draw("PSAME")
        
        ### 2tau
        elif addValues[1] == '2tau':
            # Values
            plot.drawTextFromTGraph2D(cls_2tau, size=valuesTextSize, title=numbers_title, color=new_2tau_color, xoffset=cls_xoffset, yoffset=cls_yoffset, borders=[xmin, ymin, xmax, ymax])
            # Markers
            x_mark = []; y_mark = []
            x,y,z  = cls_2tau.GetX(), cls_2tau.GetY(), cls_2tau.GetZ() 
            for i in range(0,cls_2tau.GetN()):
                if(z[i]<0): continue
                x_mark.append(x[i])
                y_mark.append(y[i])
            bullets6 = TGraph(len(x_mark),array('d',x_mark),array('d',y_mark))
            bullets6.SetMarkerStyle(8)
            bullets6.SetMarkerSize(0.4)
            bullets6.Draw("PSAME")
                
elif typee == 'both':

    if which == 'expected':
        plot.drawExpected(exp_orig_3Lonshell, color=orig_3Lonshell_color, title=f'{leg_3L} original (HF)',     legendOrder=0, linew=linew_exp)
        plot.drawExpected(exp_3Lonshell,      color=new_3Lonshell_color,  title=f'{leg_3L} (pyhf)',             legendOrder=1, linew=linew_exp)
        plot.drawExpected(exp_orig_allhad,    color=orig_allhad_color,    title=f'{leg_allhad} original (HF)', legendOrder=2, linew=linew_exp)
        plot.drawExpected(exp_allhad,         color=new_allhad_color,     title=f'{leg_allhad} (pyhf)',         legendOrder=3, linew=linew_exp)
        plot.drawExpected(exp_orig_1Lbb,      color=orig_1Lbb_color,      title=f'{leg_1Lbb} original (HF)',   legendOrder=4, linew=linew_exp)
        plot.drawExpected(exp_1Lbb,           color=new_1Lbb_color,       title=f'{leg_1Lbb} (pyhf)',           legendOrder=5, linew=linew_exp)
        plot.drawExpected(exp_orig_SS3L,      color=orig_SS3L_color,      title=f'{leg_SS3L} original (HF)',   legendOrder=6, linew=linew_exp)
        plot.drawExpected(exp_SS3L,           color=new_SS3L_color,       title=f'{leg_SS3L} (pyhf)',           legendOrder=7, linew=linew_exp)
        #plot.drawExpected(exp_orig_1L,        color=orig_1L_color,        title=f'{leg_1L} original (HF)',     legendOrder=8, linew=linew_exp)
        #plot.drawExpected(exp_1L,             color=new_1L_color,         title=f'{leg_1L} (pyhf)',             legendOrder=9, linew=linew_exp)
        plot.drawExpected(exp_orig_2tau,      color=orig_2tau_color,      title=f'{leg_2tau} original (HF)',   legendOrder=10, linew=linew_exp)
        plot.drawExpected(exp_2tau,           color=new_2tau_color,       title=f'{leg_2tau} (pyhf)',           legendOrder=11, linew=linew_exp)

    if which == 'both':
        plot.drawExpected(exp_orig_3Lonshell, color=orig_3Lonshell_color)
        plot.drawExpected(exp_3Lonshell,      color=new_3Lonshell_color)
        plot.drawExpected(exp_orig_allhad,    color=orig_allhad_color)
        plot.drawExpected(exp_allhad,         color=new_allhad_color)
        plot.drawExpected(exp_orig_1Lbb,      color=orig_1Lbb_color)
        plot.drawExpected(exp_1Lbb,           color=new_1Lbb_color)
        plot.drawExpected(exp_orig_SS3L,      color=orig_SS3L_color)
        plot.drawExpected(exp_SS3L,           color=new_SS3L_color)
        #plot.drawExpected(exp_orig_1L,        color=orig_1L_color)
        #plot.drawExpected(exp_1L,             color=new_1L_color)
        plot.drawExpected(exp_orig_2tau,      color=orig_2tau_color)
        plot.drawExpected(exp_2tau,           color=new_2tau_color)

    if which in ['both', 'observed']:
        plot.drawObserved(obs_orig_3Lonshell, color=orig_3Lonshell_color,  title=f'{leg_3L} original (HF)',   legendOrder=0) 
        plot.drawObserved(obs_3Lonshell,      color=new_3Lonshell_color,   title=f'{leg_3L} (pyhf)',           legendOrder=1) 
        plot.drawObserved(obs_orig_allhad,    color=orig_allhad_color,     title=f'{leg_allhad} original (HF)',       legendOrder=2) 
        plot.drawObserved(obs_allhad,         color=new_allhad_color,      title=f'{leg_allhad} (pyhf)',               legendOrder=3)
        plot.drawObserved(obs_orig_1Lbb,      color=orig_1Lbb_color,       title=f'{leg_1Lbb} original (HF)',         legendOrder=4) 
        plot.drawObserved(obs_1Lbb,           color=new_1Lbb_color,        title=f'{leg_1Lbb} (pyhf)',                 legendOrder=5) 
        plot.drawObserved(obs_orig_SS3L,      color=orig_SS3L_color,       title=f'{leg_SS3L} original (HF)',         legendOrder=6) 
        plot.drawObserved(obs_SS3L,           color=new_SS3L_color,        title=f'{leg_SS3L} (pyhf)',                 legendOrder=7) 
        #plot.drawObserved(obs_orig_1L,        color=orig_1L_color,         title=f'{leg_1L} original (HF)',           legendOrder=8) 
        #plot.drawObserved(obs_1L,             color=new_1L_color,          title=f'{leg_1L} (pyhf)',                   legendOrder=9) 
        plot.drawObserved(obs_orig_2tau,      color=orig_2tau_color,       title=f'{leg_2tau} original (HF)',         legendOrder=10) 
        plot.drawObserved(obs_2tau,           color=new_2tau_color,        title=f'{leg_2tau} (pyhf)',                 legendOrder=11) 

    leg = plot.createLegend(shape=leg_shape_10)
    leg.SetTextSize(legTextSize)
    leg.Draw('same')    

    # Add legends for expected and observed
    leg2.Draw('same')
    
# Save plot
plot.writePlot()




### Manual 
# canvas = TCanvas('','',1000,1000)
# #canvas.SetLeftMargin(0.14)
# #canvas.SetRightMargin(0.07)
# #canvas.SetTopMargin(0.07)
# #canvas.SetBottomMargin(0.14)
# canvas.cd()

# exp.SetLineColorAlpha(kBlack,0.9)
# exp.SetLineStyle(1)
# exp.SetLineWidth(3)
# exp.Draw('AL')

# obs.SetLineColorAlpha(kBlack,0.9)
# obs.SetLineStyle(1)
# obs.SetLineWidth(3)
# obs.Draw('L same')

# #xax = exp.GetXaxis()
# #xax.SetTitle('hellooo')
# #print(xax)

# canvas.Update()
# canvas.SaveAs('./plots/test.pdf')
