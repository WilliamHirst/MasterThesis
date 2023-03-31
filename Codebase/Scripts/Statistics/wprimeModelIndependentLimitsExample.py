
import pyStats;
import limitPlot;

#create the limit plot object
limPlot = limitPlot.limitPlot('Muon channel model independent limits');

#for each mT-threshold, we create a counting experiment object and add a point to our limit plot
countexp = pyStats.countingExperiment('mTmin = 1882 GeV', intLum = 139e3, intLumUnc = 2.4e3);
countexp.addChannel('Muon channel', bkg = 27.1, bkgUnc = 5.8, Nobs = 30);
limPlot.addPoint(1882e-3, countexp); 

countexp = pyStats.countingExperiment('mTmin = 2046 GeV', intLum = 139e3, intLumUnc = 2.4e3);
countexp.addChannel('Muon channel', bkg = 17.6, bkgUnc = 4.3, Nobs = 21);
limPlot.addPoint(2046e-3, countexp); 

countexp = pyStats.countingExperiment('mTmin = 2224 GeV', intLum = 139e3, intLumUnc = 2.4e3);
countexp.addChannel('Muon channel', bkg = 11.4, bkgUnc = 3.1, Nobs = 16);
limPlot.addPoint(2224e-3, countexp); 

countexp = pyStats.countingExperiment('mTmin = 2418 GeV', intLum = 139e3, intLumUnc = 2.4e3);
countexp.addChannel('Muon channel', bkg = 7.4, bkgUnc = 2.2, Nobs = 8);
limPlot.addPoint(2418e-3, countexp); 

countexp = pyStats.countingExperiment('mTmin = 2628 GeV', intLum = 139e3, intLumUnc = 2.4e3);
countexp.addChannel('Muon channel', bkg = 4.7, bkgUnc = 1.6, Nobs = 5);
limPlot.addPoint(2628e-3, countexp); 

countexp = pyStats.countingExperiment('mTmin = 2857 GeV', intLum = 139e3, intLumUnc = 2.4e3);
countexp.addChannel('Muon channel', bkg = 3.0, bkgUnc = 1.1, Nobs = 3);
limPlot.addPoint(2857e-3, countexp); 

countexp = pyStats.countingExperiment('mTmin = 3106 GeV', intLum = 139e3, intLumUnc = 2.4e3);
countexp.addChannel('Muon channel', bkg = 1.91, bkgUnc = 0.78, Nobs = 2);
limPlot.addPoint(3106e-3, countexp); 

countexp = pyStats.countingExperiment('mTmin = 3377 GeV', intLum = 139e3, intLumUnc = 2.4e3);
countexp.addChannel('Muon channel', bkg = 1.20, bkgUnc = 0.54, Nobs = 2);
limPlot.addPoint(3377e-3, countexp); 

countexp = pyStats.countingExperiment('mTmin = 3671 GeV', intLum = 139e3, intLumUnc = 2.4e3);
countexp.addChannel('Muon channel', bkg = 0.75, bkgUnc = 0.37, Nobs = 1);
limPlot.addPoint(3671e-3, countexp); 

countexp = pyStats.countingExperiment('mTmin = 3990 GeV', intLum = 139e3, intLumUnc = 2.4e3);
countexp.addChannel('Muon channel', bkg = 0.47, bkgUnc = 0.25, Nobs = 1);
limPlot.addPoint(3990e-3, countexp); 

countexp = pyStats.countingExperiment('mTmin = 4338 GeV', intLum = 139e3, intLumUnc = 2.4e3);
countexp.addChannel('Muon channel', bkg = 0.29, bkgUnc = 0.16, Nobs = 1);
limPlot.addPoint(4338e-3, countexp); 

countexp = pyStats.countingExperiment('mTmin = 4716 GeV', intLum = 139e3, intLumUnc = 2.4e3);
countexp.addChannel('Muon channel', bkg = 0.18, bkgUnc = 0.11, Nobs = 1);
limPlot.addPoint(4716e-3, countexp);

countexp = pyStats.countingExperiment('mTmin = 5127 GeV', intLum = 139e3, intLumUnc = 2.4e3);
countexp.addChannel('Muon channel', bkg = 0.11, bkgUnc = 0.07, Nobs = 0);
limPlot.addPoint(5127e-3, countexp); 

#then we calculate the limits ...
limPlot.calculate();

#... print the results to screen ...
print(limPlot);

#... and plot them
limPlot.drawPlot(xtitle = 'm_{T}^{min} [TeV]', ytitle = '#sigma_{vis} [pb]', yrange = [1e-5,1e-2]);

#these lines prevent the script from finishing such that the plot remains on screen
from time import sleep;
sleep(2000);
