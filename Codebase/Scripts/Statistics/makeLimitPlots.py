
import pyStats;
import limitPlot;

limitPlots = {'combined': limitPlot.limitPlot('Combined') };

#peek in the input file to check which channels are there
inputFile = open('inputs.txt','r');
lines = inputFile.readlines();
for l in lines:

    exec(l); #Input file consists of valid Python statements such as mass=300 etc.
    
    if 'channel=' in l:
        if not channel in limitPlots.keys():
            limitPlots[channel] = limitPlot.limitPlot(channel);

print( 'Here is what we have so far:' );
print(limitPlots);

countexp = {};
xtitle='Mass';
ytitle='Cross section';
yrange=[-999.0,-999.0];

#read the inputs and fill in the respective limit plots with background levels, observed counts, etc.
for l in lines:

    exec(l); #Input file consists of valid Python statements such as mass=300 etc.

    if 'mass=' in l: #Such a line defines a new point in the limit plot
        for channel in limitPlots.keys():
            countexp[channel] = pyStats.countingExperiment(name = 'mass = '+str(mass), intLum = intLum, intLumUnc = intLumUncertainty);
            limitPlots[channel].addPoint(mass, countexp[channel], theoryCrossSection);

    if 'channel=' in l: #Such a line gives the inputs for a given channel
        countexp['combined'].addChannel(name = channel, bkg = background, bkgUnc = backgroundUncertainty, Nobs = Nobs, eff = efficiency, effUnc = efficiencyUncertainty);
        countexp[channel].addChannel(name = channel, bkg = background, bkgUnc = backgroundUncertainty, Nobs = Nobs, eff = efficiency, effUnc = efficiencyUncertainty);

#do the actual calculations and plotting
for key in limitPlots.keys():

    #if you are testing things on your laptop, maybe you want to drop
    #the combined limits  because they take a lot longer to run
    #if key == 'combined':
    #    continue;

    print( limitPlots[key] );
    print( 'Running limits for ', limitPlots[key].__repr__() );

    limitPlots[key].calculate();
    limitPlots[key].dumpToFile('limits_'+key+'.txt');
    limitPlots[key].drawPlot(xtitle, ytitle, yrange, filename='limitPlot_'+key+'.png');
    

#for interactive running, these lines prevent the program from closing
#such that the plots remain on screen and can be manipulated
#from time import sleep;
#sleep(2000);
