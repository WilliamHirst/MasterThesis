# FYS5555 Code Example for 13 TeV ATLAS Open Data

Prepared by Even Håland [edited by Eirik Gramstad January 2022]

This is a complete code example for use with the ATLAS 13 TeV Open Data. Use it to
get started quickly with your analysis. Try first to run the code as it is and check
that everything works, then you can move on to edit the relevant code files to
perform your analysis of interest.

## Running the code

Before running the code you need all the required software to be installed, e.g. by following the instructions here: https://github.uio.no/zpath/software#installation-instructions.

Then get hold of the 2-lepton skim (or any other skim of interest) of ATLAS Open Data from here: http://opendata.cern.ch/record/15003
(download the .zip file and extract it to a convenient directory).

The code is run by 
  
- `python GenericRunSelector.py <options>`

Do `python GenericRunSelector.py --help` to see the available options.

*Note* : one needs to change the `rootpath` on L89 in `GenericRunSelector.py` if one is not running on the hepp02 node.

Run the plotting for the electron channel:
- `python MakePlots.py ee`

Run the plotting for the muon channel:
- `python MakePlots.py uu`

If all commands ran without errors and you have some nice-looking plots in the Plots/ directory, you are
ready to move on to modify the analysis according to your needs.

## Modifying the code

To perform the analysis you are interested in, you will likely need to modify the event selection and
add additional histograms. Note that if your final state of interest does not contain (at least) two
leptons, you will need a different dataset (different skim), see here:
http://opendata.cern.ch/search?page=1&size=20&experiment=ATLAS&collision_energy=13TeV

### Changing the event selection (cuts)

The example analysis selects events with exactly two leptons of the same generation (electrons or muons)
with opposite charge. There are also additional cuts, mostly related to the quality of the two
leptons. The cuts are implemented in the file MySelector.C, for example:
```
// Require same flavour (2 electrons or 2 muons)
if(lep_type[0] != lep_type[1]){ return kTRUE; } 
```
This cut selects events where the lepton type is the same for the two leptons, i.e. they are both
electrons or both muons. If the lepton type is not the same, we hit the "return" statement, and
all the below code is not executed, including the filling of the histograms:
```
  // Fill histograms

  h_pt1[channel]->Fill(l1.Pt(), wgt); 
  h_pt2[channel]->Fill(l2.Pt(), wgt); 
  h_eta1[channel]->Fill(l1.Eta(), wgt); 
  h_eta2[channel]->Fill(l2.Eta(), wgt); 
  h_phi1[channel]->Fill(l1.Phi(), wgt); 
  h_phi2[channel]->Fill(l2.Phi(), wgt); 
  h_met[channel]->Fill(*met_et/1000.0, wgt); 
  h_mll[channel]->Fill(dileptons.M(), wgt);

```
Hence the histograms are only filled for events that pass all selection cuts. You can try adding
a cut requiring the missing transverse energy (MET) to be above 100 GeV like this:
```
//MET cut
if(*met_et < 100.e3){return kTRUE;}
```
(Note that variables with dimension of energy is stored in unit MeV, hence 100 GeV corresponds to 100.e3.)
Try rerunning the code and plotting. Which background is now dominating your invariant mass plot?

### Adding histograms

You will need to add histograms if you
- want to look at additional variables
- want to consider variations of the event selection without rerunning the code

As an example, let's consider the second use case. Say we want to look at the invariant mass distribution
both for the inclusive dilepton selection and the special case MET > 100 GeV considered above, but only
running the analysis code once. To add a new histogram, let's follow the example of the invariant mass
histogram already implemented. You can search for "h_mll" in the file MySelector.C, and you find the
following lines:

The creation of the histogram:
```
h_mll[chn] = new TH1D("h_"+chn+"_mll", chn+"_mll", 50, 0, 3000); 
```
(Note that here the histogram is defined to have 50 bins from 0 to 3000.)
The filling of the histogram:
```
h_mll[channel]->Fill(dileptons.M(), wgt);
```
The saving of the histogram to the output file:
```
h_mll[chn]->Write();
```
The resetting of the histogram when transitioning to a new MC sample:
```
h_mll[chn]->Reset();
```

Let's add these lines (in the different appropriate places) to add an invariant
mass histogram for events with MET > 100 GeV:
```
h_mll_highMET[chn] = new TH1D("h_"+chn+"_mll_highMET", chn+"_mll_highMET", 50, 0, 3000);
if(*met_et > 100.e3){h_mll_highMET[channel]->Fill(dileptons.M(), wgt);}
h_mll_highMET[chn]->Write();
h_mll_highMET[chn]->Reset();
```
Note the if statement, that makes sure this histogram is only filled with events that have
missing transverse energy above 100 GeV. We will also need to declare the histogram in the
file MySelector.h:
```
map<TString, TH1*> h_mll_highMET; //!
```
Note that this is actually a map of histograms, with individual histograms for the electron and
muon channels.

You can now run the modified code over data and MC and run the plotting script again, but you
may notice that the plotting script does not make a plot corresponding to the new histogram.
To get your new histogram plotted, you just need to modify a couple of lines near the top of
MakePlots.py, which after modification may look like this:
```
variables = ['pt1', 'pt2', 'eta1', 'eta2', 'phi1', 'phi2', 'mll', 'met', 'mll_highMET']

xtitles = {'pt1':'Leading lepton p_{T} (GeV)', 'pt2':'Subleading lepton p_{T} (GeV)', 'eta1':'Leading lepton #eta', 'eta2':'Subleading lepton #eta', 'phi1':'Leading lepton #phi', 'phi2':'Subleading lepton #phi', 'mll':'m_{ll} (GeV)', 'met':'E_{T}^{miss} (GeV)', 'mll_highMET':'m_{ll} (GeV)'}
```
Run the plotting and look at your new histogram.

To check that you have understood how this works: Add a new histogram showing the distribution of the
angle between the leptons in the transverse plane. You may calculate this angle as
```
l1.DeltaPhi(l2)
```
