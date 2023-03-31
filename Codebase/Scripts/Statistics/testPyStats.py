
import pyStats;
from ROOT import TMath;

countexp = pyStats.countingExperiment();
countexp.addChannel('myChannel', bkg = 100.0, bkgUnc = 10.0, Nobs = 120);
print(countexp);
countexp.setPriors(bkg='gauss');
print( "Significance = ",countexp.getSignificance() );
print( "Significance ( s/sqrt(b + sigma_b^2) ) = ",(120-100.)/TMath.Sqrt(100. + 10.**2) );

countexp = pyStats.countingExperiment();
countexp.addChannel('myChannel', bkg = 1.0, bkgUnc = 0.1, Nobs = 0);
print(countexp);
print( "Observed limit = ",countexp.getBayesianLimit() );

#-2sigma,-1sigma,expected,+1sigma,+2sigma
m2sig,m1sig,expected,p1sig,p2sig = countexp.getBayesianExpectedLimit();

print( "\nExpected limit and bands:" );
print( "  -2sigma                 -1sigma                 median                +1sigma                +2sigma" );
print( m2sig, "   ", m1sig, "   ", expected, "   ", p1sig, "   ", p2sig );


countexp = pyStats.countingExperiment(name = '300 GeV', intLum = 20.28, intLumUnc = 0.5678);
countexp.addChannel('electron channel', bkg = 12901.4, bkgUnc = 816.767129814, Nobs = 12717, eff = 0.22758, effUnc = 0.00930046154813);
countexp.addChannel('muon channel', bkg = 11266.4, bkgUnc = 773.744616417, Nobs = 10927, eff = 0.183817, effUnc = 0.00706668179581);
print(countexp);
print( "Observed limit = ",countexp.getBayesianLimit() );
m2sig,m1sig,expected,p1sig,p2sig = countexp.getBayesianExpectedLimit(Npseudos=200);
print( "\nExpected limit and bands:" );
print( "  -2sigma                 -1sigma                 median                +1sigma                +2sigma" );
print( m2sig, "   ", m1sig, "   ", expected, "   ", p1sig, "   ", p2sig );



import limitPlot;

limPlot = limitPlot.limitPlot('Electron and muon channels combination');
limPlotEl = limitPlot.limitPlot('Electron channel');
limPlotMu = limitPlot.limitPlot('Muon channel');

inputFile = open('inputs.txt','r');
for l in inputFile.readlines():

    exec(l); #Input file consists of valid Python statements such as mass=300 etc.

    if 'mass=' in l:
        countexp = pyStats.countingExperiment(name = str(mass)+' GeV', intLum = intLum, intLumUnc = intLumUncertainty);
        limPlot.addPoint(mass, countexp, theoryCrossSection);

        countexpEl = pyStats.countingExperiment(name = str(mass)+' GeV', intLum = intLum, intLumUnc = intLumUncertainty);
        limPlotEl.addPoint(mass, countexpEl, theoryCrossSection);

        countexpMu = pyStats.countingExperiment(name = str(mass)+' GeV', intLum = intLum, intLumUnc = intLumUncertainty);
        limPlotMu.addPoint(mass, countexpMu, theoryCrossSection);

    if 'channel=' in l:
        countexp.addChannel(name = channel, bkg = background, bkgUnc = backgroundUncertainty, Nobs = Nobs, eff = efficiency, effUnc = efficiencyUncertainty);

        if channel == 'electron':
            countexpEl.addChannel(name = channel, bkg = background, bkgUnc = backgroundUncertainty, Nobs = Nobs, eff = efficiency, effUnc = efficiencyUncertainty);
        if channel == 'muon':
            countexpMu.addChannel(name = channel, bkg = background, bkgUnc = backgroundUncertainty, Nobs = Nobs, eff = efficiency, effUnc = efficiencyUncertainty);

print(limPlotMu);

limPlotMu.calculate();
limPlotMu.dumpToFile();

print(limPlotMu);

limPlotMu.drawPlot(xtitle="W' mass [GeV]", ytitle='Cross section [fb]', yrange=[5.0e-2,5.0e3]);

from time import sleep;
sleep(2000);
