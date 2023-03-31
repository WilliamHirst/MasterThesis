# SUSYPheno (Eirik Gramstad)
A program for doing scans of pMSSM models taking into account constraints from various measurements. 
The core code is written by Børge Kile Gjelsten and modified/extended by Eirik Gramstad.

## First time set up do
```
git clone https://github.com/eirikgr/SUSYPheno.git
cd SUSYPheno
source install.sh
```
## Every time you start a new terminal
It is important to add the environment variable *SUSYPHENO_PATH* pointing to the location of the clone of the SUSYPheno directory above. To do so issue the following two commands in the *SUSYPheno* directory: 
```
export SUSYPHENO_PATH=$PWD
export PATH=$PATH:$SUSYPHENO_PATH/bin
```
If you want to enable plotting you would also need to have ROOT installed ([ROOT Download Page](https://root.cern.ch/downloading-root))
## How to run the program
In order to do the plotting you will need ROOT (see [ROOT Download Page](https://root.cern.ch/downloading-root) on how to install root). You can continue to do the parameter scans, but the plotting will not work if ROOT is not installed.

## Make a first test run
Create a work directory, e.g. `susyscans`, some where on you computer
```
cd susyscans
mkdir scantest1
cd scantest1
```
Then calculate 1 MSSM point (with initial parameters set in the command).
```
signalgrid_loop.py  -fnadd _testgrid1  -mssmloop M1=100:M2=250:mu=400:TB=10:mA=2000 
```

The output is mostly given in so-called tlt-format: the first line shows the quantity name, the second line shows the value (for the given scenario). The output files are as follows (where <pointdesc> is susy_testgrid1_M1_100_M2_250_mu_400_TB_10_mA_2000 for the above example):

```
<pointdesc>_scanvars_all.tlt        : just the variables specified on input line (the ones not taking default values)
<pointdesc>_scanvars_multi.tlt      : similar, but only vars which on input line has more than one value
<pointdesc>_susyhit_slha_mass.tlt   : input vars + all SUSY masses
<pointdesc>_susyhit_slha_extpar.tlt : input vars + ~all Lagrangian parameters
<pointdesc>_susyhit_slha_elweak.tlt : input vars + parameters relevant for elweak sector
<pointdesc>_susyhit_slha_rest.tlt   : input vars + some additional vars (the rest)
<pointdesc>_susyhit_slha_mixing.tlt : input vars + mixing parameters in squark/slepton/neutralino/chargino sector
<pointdesc>_susyhit_errlowtune.tlt  : input vars + errors + finetuning parameters + experimentally constrained vars 
<pointdesc>_susyhit_MO.tlt          : input vars + MicrOmegas output (DarkMatter++, experimentally constrained vars)
<pointdesc>_susyhit_ds.tlt          : input vars + DarkSUSY output (DarkMatter, experimental bounds)
```

To see the above description and some more details on options etc. 
```
signalgrid_loop.py -h 
signalgrid_loop.py --showpmssmpars   : shows the notation for the 24 pmssm parameters
signalgrid_loop.py --showpmssmvals   : shows the default value of the 24 pmssm parameters
```

## Make a first scan
```
cd susyscans
mkdir scantest2
cd scantest2
```
Scan 7x7 in M2 and mu; Set M1,TB (=tan(beta)), mA ; the rest takes default (high values) and are then mostly not relevant at LHC
```
signalgrid_loop.py  -fnadd _testgrid1  -mssmloop M1=100:M2=100,150,200,250,300,350,400:mu=100,150,200,250,300,350,400:TB=10:mA=2000 
```
This takes a few minutes to run and produces for each of 7x7 pMSSM setups some ~16 files. These files contain the relevant output. Ex: Dark matter values are found in the files *_ds.tlt. 
To see the filenames do
```
ls -l *_ds.tlt  
```
To see the file content do 
```
cat *_ds.tlt
```
Make things a bit nicer by lining up and keeping header only on the first line:
```
cat *_ds.tlt | lineup_removeduplicates.sh
```
### Example: Check Dark Matter relic density
Look at the *DMbestVSexp* and see if predictions is
1. above 1: SUSY scenario predicts more DM than observed (i.e model is EXCLUDED). 
2. exactly 1 (within uncertainties): this SUSY scenario will account for all of the Dark Matter we have.  
3. below 1: this SUSY scenario is fine (in this variable), but does not account for all the Dark Matter. (There needs to be additional physics out there to give the rest of the Dark Matter.)

## Plot result of scan (**need ROOT**)

Make nice table of your tlt-files (e.g. from the scan above)
```
cat *_ds.tlt | lineup_removeduplicates.sh > table_ds.txt
```
And check if they look reasonable by inspecting it
```
less table_ds.txt
```
Then make a plot of the SUSY DM relative to experimental value (*DMwCoVSexp*) as a function of the two free parameters (*mu* and *M2*)
```
cat table_ds.txt | munch.py -coord mu,M2  -resvars DMwCoVSexp 
```
This produces a pdf-file (*much.pdf*) with your plot.

### More on plotting
Below are some useful options to refine the plot: 
Rename axes
```
cat table_ds.txt | munch.py -coord mu:#mu,M2:M_{2}  -resvars DMwCoVSexp  
```
Give new title
```
cat table_ds.txt | munch.py -coord mu:#mu,M2:M_{2}  -resvars DMwCoVSexp  -title 'Dark Matter with coannihilation VS experiment'  
```
Add numbers onto canvas
```
cat table_ds.txt | munch.py -coord mu:#mu,M2:M_{2}  -resvars DMwCoVSexp  --numbers  
```
Add a countour at 10, as well as changed the colour (to red, 2), the style (to 7) and the width (to 5)
```
cat table_ds.txt | munch.py -coord mu:#mu,M2:M_{2}  -resvars DMwCoVSexp  -cont DMwCoVSexp::1,10:2:7:5  
```
Specify a meaningful name on the output file
```
cat table_ds.txt | munch.py -coord mu:#mu,M2:M_{2}  -resvars DMwCoVSexp  -save DarkMatter1.pdf  # specify a meaningful filename
```
To see some of the options in munch.py, type 
```
munch.py -h 
```

## Look at other variables
The table  in the example above only have input from darksusy (ds). You might want to study other variables. In this case it is practical to make some more tables: 
```
cat *_elweak.tlt | lineup_removeduplicates.sh > table_elweak.txt
cat *_masses.tlt | lineup_removeduplicates.sh > table_masses.txt
cat *_errlowtune.tlt | lineup_removeduplicates.sh > table_errlowtune.txt
cat *_MO_mini.tlt | lineup_removeduplicates.sh > table_MO.txt
```
Tables can be combined (horizontally): 
```
paste table_elweak.txt table_ds.txt > table_elweak_ds.txt
```
The above is useful if the variables you want to use are not in the same table. With the combined table you can for instance plot the same as above, but with the chargino mass in contour at 100 GeV
```
cat table_elweak_ds.txt | munch.py -coord mu:#mu,M2:M_{2}  -resvars DMwCoVSexp  -cont ~C1::100
```
You can also add two contrours (chargino mass in contour at 100 GeV (in red) and *DMwCoVSexp* equals 1 (in green)
```
cat table_elweak_ds.txt | munch.py -coord mu:#mu,M2:M_{2}  -resvars DMwCoVSexp  -cont DMwCoVSexp::1:2:7:5,,~C1::100:3:4:4
```


## pMSSM Investigations

With the above setup you can scan the full pMSSM (24 input parameters) and get 
info about things like
 - sparticle masses [can compare against limits we have already set]  (table_masses.txt)
 - Dark Matter values [can compare against experiment]   (table_ds.txt, table_MO.txt)
 - Additional variables like g_mu-2, b->s+gamma [can compare against experiment]  (table_MO.txt)
 - Finetuning parameters  [these give some measure of the *likeliness* of a given scenario. Values much above 10 are on the unlikely side] (table_errlowtune.txt)

Many of these variables can be meaningfully investigated: compare against experimental value, against theoretical likeliness. 
It is therefore possible to get an idea about what input parameters give SUSY models 
which are already excluded (ex: if too much Dark Matter, if lightest Higgs is wrong, if too light sleptons, if too light squarks, ...) 
etc. 

One of the intentions about doing these pMSSM scans, is to learn a bit about 
how the input parameters affect measurable quantities. 

Since we have 24 parameters to vary, we need to first know where to start, we need default values. 
*signalgrid_loop.py* has a default setup where: 
 - all squark and slepton values are set to a *heavy* mass of 3 TeV (they would then not be produced at any noticeable rate at the LHC)
 - the gluino (M3) is also set to heavy (3 TeV)  (will also not be produced)
 - the trilinear couplings (At,Ab,AT) are set to 3400 TeV
 - mA = 2000 GeV  (controls the heavy Higgses)
 - tan(beta) = 10  (TB) 
 - (M1,M2,mu) = (100,200,500) GeV  [these, together with tan(beta), control the neutralino/chargino sector)
This is just a default, not motivated by any underlying physics. 
But it gives a simple scenario type since all squarks and sleptons are pushed high (and out of the way). 
With this default setup you can start investigating how output depends on input variables.  

### Tips for further studies

#### STUDY 1: 

Some relations are very straightforward. E.g. input eR controls the output mass of the right-handed selectron.

1.1 CHECK IF THIS IS THE CASE:
```
signalgrid_loop.py  -fnadd _ex1  -mssmloop eR=3000,300
cat susy_ex1_eR_300*masses.tlt | lineup_removeduplicates.sh   # you don't need to plot, just study the table
```
1.2 Whatever is true for the first generation of sleptons/squarks, should also be true for the second. CHECK IF THAT IS THE CASE (CHECK mR).

1.3 For lefthanded sleptons it is slightly more complicated. CHECK WHAT MASSES ARE AFFECTED WHEN YOU CHANGE L1.

1.4 THEN TRY L3.

1.5 THEN ALSO TRY TR (TOGETHER WITH L3).

1.6 WHICH SPARTICLES ARE CONTROLLED BY M3?

1.7 TRY CHANGING At AND CHECK THE EFFECT ON THE MASS OF THE LIGHTEST HIGGS, mh

1.8 HOW DID THE CHANGES ABOVE AFFECT THE DARK MATTER? (Check OmegaVSexp (what is this variable?) in susy_ex1*susyhit_MO_mini.tlt)

1.9 TRY CHANGING M1 AND SEE IF THAT AFFECTS DARK MATTER. (AND HOW DOES IT AFFECT THE MASS OF N1?)

1.10 TRY CHANGING M2 AND mu AND SEE HOW THAT AFFECTS THE CHARGINOS. ARE ANY OTHER QUANTITIES AFFECTED? 

#### STUDY 2: Neutralino sector
The neutralino is governed mainly by 3 parameters, M1, M2, mu 
(TB is also relevant, but can be neglected in a first approximation). 
These three input parameters control four netralinos. 
**M1**: bino-part, **M2**: wino-part, **mu**: higgsino-part. 
The neutralinos are in general linear combinations of the bino, wino and higgsino. 
In scenarios with considerable difference between some of M1, M2, mu we will have 
cases where a given neutralino can be e.g. very nearly pure bino. 
In the default setup of *signalgrid_loop.py*, this is roughly the case.

2.1 WHICH NEUTRALINOS ARE (MOSTLY) DETERMINED BY WHICH OF M1, M2 and mu?

2.2 DO VARIOUS USEFUL SCANS IN THE THREE VARIABLES AND TRY TO UNDERSTAND & DESCRIBE HOW THE NEUTRALINO SECTOR WORKS (HOW THE NEUTRALINO MASSES ARE GIVEN BY M1,M2,mu)

2.3 STUDY DARK MATTER AS FUNCTION OF M1,M2,mu  

