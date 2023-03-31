import os;
from array import array;
import ROOT;
import multiprocessing as mp;
import warnings;
warnings.filterwarnings( action='ignore', category=RuntimeWarning, message='creating converter.*' );

#The C++ class that does all the calculations
ROOT.gROOT.ProcessLine(".L MC_Bayes_withBkgCorr.cpp+");

"""
Class for a channel in a counting experiment. A channel is characterized by
- background level and its uncertainty
- number of observed events
- signal efficiency and its uncertainty
Additionally, the channel can be given a name for bookkeeping/introspection
purposes. When correlations are treated between the background levels in the
different channels, the background uncertainty is split into a correlated
and an uncorrelated component.
"""
class channel:

    def __init__(self,name,bkg,bkgUnc,Nobs,eff,effUnc,bkgUncUncorr=-999.0,bkgUncCorr=-999.0):

        self.bkg = bkg;
        self.bkgUnc = bkgUnc;
        self.Nobs = Nobs;
        self.eff = eff;
        self.effUnc = effUnc;
        self.name = name;

        if bkgUncUncorr != -999.0:
            self.bkgUncUncorr = bkgUncUncorr;
            self.bkgUncCorr = bkgUncCorr;
            self.bkgUnc = ROOT.TMath.Sqrt( bkgUncUncorr*bkgUncUncorr + bkgUncCorr*bkgUncCorr );
            

"""
Class for a counting experiment for the combination of an arbitrary number
of channels without correlation between the background levels in the channels.
When the background levels are correlated (correlated systematic uncertainty),
the derived class "countingExperimentWithBkgCorr" defined further below should
be used.
"""
class countingExperiment:

    def __repr__(self):
        return 'Counting experiment "'+self.name+'"';

    def __str__(self):

        retVal = '';
        retVal += '---------------------------\n';
        retVal += 'Counting experiment "'+self.name+'"\n';
        retVal += '---------------------------\n';
        retVal += "Int. luminosity = "+str(self.intLum)+" +/- "+str(self.intLumUnc);
        retVal += '\n---------------------------\n';
        for chan in self.channels:
            retVal += 'Channel "'+chan.name+'":\n';
            retVal += self.bkgString(chan);
            retVal += "   Observed events = "+str(chan.Nobs)+'\n';
            retVal += "   Signal efficiency = "+str(chan.eff)+" +/- "+str(chan.effUnc)+'\n';
        return retVal;

    def bkgString(self, chan):
        return "   Background = "+str(chan.bkg)+" +/- "+str(chan.bkgUnc)+"\n";

    #Constructor
    def __init__(self, name='myCountingExperiment', intLum=1.0, intLumUnc=0.0):

        self.intLum = intLum;
        self.intLumUnc = intLumUnc;
        self.name = name;

        self.channels = [];

        self.background = array('d',[]);
        self.backgroundRelUncUncorr = array('d',[]);
        self.backgroundRelUncCorr = array('d',[]);

        self.Nobs = array('i',[]);
        self.efficiency = array('d',[]);
        self.efficiencyUnc = array('d',[]);
        
        self.Nmc = 5000;
        self.stepLength = -999.0;

        self.numCores = int(os.getenv('SLURM_JOB_CPUS_PER_NODE', os.cpu_count()));

        self.doLogNormalBkg = True;
        self.doLogNormalEff = True;
        self.doLogNormalLumi = True;

        self.posterior = 0;

    #Function to add a channel with given background level and uncertainty,
    #number of observed events, and signal efficiency and uncertainty
    def addChannel(self, name, bkg, bkgUnc, Nobs, eff=1.0, effUnc=0.0):

        self.channels.append(channel(name,bkg,bkgUnc,Nobs,eff,effUnc));

    #Choose between Gaussian or log-normal priors for the background level,
    #signal efficiency, and integrated luminosity
    def setPriors(self, bkg='lognormal', eff='lognormal', lumi='lognormal'):

        if bkg == 'lognormal':
            self.doLogNormalBkg = True;
        elif bkg == 'gauss':
            self.doLogNormalBkg = False;
        else:
            print( 'Please choose background prior/constraint term as "lognormal" or "gauss", defaulting to "lognormal"' );
            self.doLogNormalBkg = True;

        if eff == 'lognormal':
            self.doLogNormalEff = True;
        elif eff == 'gauss':
            self.doLogNormalEff = False;
        else:
            print( 'Please choose signal efficiency prior/constraint term as "lognormal" or "gauss", defaulting to "lognormal"' );
            self.doLogNormalEff = True;

        if lumi == 'lognormal':
            self.doLogNormalLumi = True;
        elif lumi == 'gauss':
            self.doLogNormalLumi = False;
        else:
            print( 'Please choose integrated luminosity prior/constraint term as "lognormal" or "gauss", defaulting to "lognormal"' );
            self.doLogNormalLumi = True;

    #To set a reasonable step length automatically, a very simple analytical calculation is
    #performed to obtain the correct order of magnitude of the limit based on the most sensitive channel.
    def getApproxLimit(self):

        approxLimit = 1.e100;
        for chan in self.channels:
            approxLimitEvts = 1.6 * ROOT.TMath.Sqrt(chan.bkgUnc**2 + chan.bkg + (chan.bkg*self.intLumUnc/self.intLum)**2);
            if approxLimitEvts < 3.0:
                approxLimitEvts = 3.0;
            thisApproxLimit = approxLimitEvts/(chan.eff*self.intLum) if chan.eff != 0.0 else 1.0e100;
            if thisApproxLimit < approxLimit:
                approxLimit = thisApproxLimit;
        return approxLimit;

    #Gives the step length that has been set by the user or auto-adjusted. Gives an initial guess based
    #on the analytical calculation if no step length has been set.
    def getStepLength(self):

        if self.stepLength == -999.0:
            return self.getApproxLimit()/500.0;
        else:
            return self.stepLength;

    #Counting experiments with large backgrounds are treated differently in the parameter optimization. In
    #this case, one can get outliers in the limit calculation due to the impact of systematics, and it's
    #better to use fewer pseudo-experiments with higher precision.
    def isHighStat(self):
        highStat = False;
        for chan in self.channels:
            if chan.bkg > 2500.:
                highStat = True;
        return highStat;

    #Counting experiments with low backgrounds are treated differently in the parameter optimization. In
    #this case the discreteness is obvious, and the caching of already calculated limits can lead to a
    #saturation of running time, so that one can afford a very large number of pseudo-experiments.
    def isLowStat(self):
        lowStat = True;
        for chan in self.channels:
            if chan.bkg > 5.:
                lowStat = False;
        return lowStat;

    #Optimized parameters for pseudo-experiment calculations. Low and high statistics are treated differently
    #as explained above.
    def getTunedParameters(self):

        Npseudos = 1000;
        step = self.getStepLength()*2.0;
        Nmc = 1500;

        if self.isHighStat():
            Npseudos = int(Npseudos/2);
            step *= 2;
            Nmc = int(Nmc*5.0);

        if self.isLowStat():
            Npseudos = 10000;

        return step,Nmc,Npseudos;

    #Prepares the arrays of background levels and uncertainties, number of observed events, and efficiencies
    #and uncertainties. These arrays are the inputs to the C++ code that does the actual calculations. Note
    #that the correlated component of the background uncertainty is set to zero. This function is overridden
    #in the derived class for counting experiments with background correlation, where both the correlated and
    #uncorrelated components are set.
    def prepareArrays(self):

        self.background = array('d',[]);
        self.backgroundRelUncUncorr = array('d',[]);
        self.backgroundRelUncCorr = array('d',[]);
        self.Nobs = array('i',[]);
        self.efficiency = array('d',[]);
        self.efficiencyUnc = array('d',[]);

        for chan in self.channels:
            self.background.append(chan.bkg/self.intLum);
            self.backgroundRelUncUncorr.append(chan.bkgUnc/chan.bkg);
            self.backgroundRelUncCorr.append(0.0);
            self.Nobs.append(chan.Nobs);
            self.efficiency.append(chan.eff);
            self.efficiencyUnc.append(chan.effUnc);

    #Sets up the C++ object that does the actual calculations
    def prepareBayesCalculator(self, Nmc = -1):

        if len(self.channels) == 0:
            raise RuntimeError('No channels have been added for counting experiment "'+self.name+'"');

        if Nmc == -1:
            Nmc = self.Nmc;
        Nchannels = len(self.channels);
        self.prepareArrays();

        return ROOT.MC_Bayes_withBkgCorr(Nchannels, self.intLum, self.intLumUnc, self.background, self.backgroundRelUncUncorr,  \
                                         self.Nobs, Nmc, self.efficiency, self.efficiencyUnc, self.backgroundRelUncCorr, False, \
                                         self.doLogNormalEff, self.doLogNormalLumi, self.doLogNormalBkg);

    #Calculates the significance of an excess using the asymptotic approximation, i.e. sqrt(q0).
    #The maximum value of the signal parameter used in the Minuit fit can be specified. In case the
    #Minuit fit (for signal+background) fails, you can try to tweak the maxSignal parameter.
    def getSignificance(self, verbose = False, maxSignal = -999.0):

        if maxSignal == -999.0:
            maxSignal = self.getApproxLimit()*10.0;
        
        return ROOT.TMath.Sqrt(self.prepareBayesCalculator().q0TestStatistic(maxSignal,verbose));

    #Calculates the observed Bayesian limit at the given confidence (or credibility) level (default: 95% CL)
    def getBayesianLimit(self, CL = 0.95):

        mcb = self.prepareBayesCalculator();
        mcb.confidenceLevel = CL;

        step = self.getStepLength();
        limit = mcb.excludedSignal(step, False);

        self.posterior = mcb.getPosteriorDistribution();

        #In case the initial guess for the step length was not good, we can adjust it now
        if limit/step < 200.0 or limit/step > 2000.0:
            print( "INFO: Adjusting step length and retrying limit calculation" );
            if limit < 0.0:
                print( "WARNING: Negative limit encountered, you may want to check the inputs..." );
                self.stepLength = self.getStepLength() / 10.0;
            else:
                self.stepLength = limit/500.0;
            return self.getBayesianLimit(CL);
            
        return limit;

    #Returns the marginalized posterior as a TGraph
    def getPosterior(self):

        if self.posterior == 0:
            self.getBayesianLimit();

        return self.posterior;

    #Calculates the expected Bayesian limit as well as the 1-sigma and 2-sigma fluctuations
    #(i.e. green and yellow bands). The number of pseudo-experiments and confidence (credibility)
    #level can be specified.
    def getBayesianExpectedLimit(self, Npseudos = -1, CL = 0.95):

        #In case of only one channel, we can run many pseudo-experiments very fast
        if len(self.channels) == 1:
            Npseudos = 100000;
            mcb = self.prepareBayesCalculator();
            mcb.confidenceLevel = CL;
            exp = mcb.expectedExclusion(self.getStepLength(),Npseudos,False);
            return exp[0], exp[1], exp[2], exp[3], exp[4];

        #If the user did not specify Npseudos, use pre-tuned accuracy parameters
        if Npseudos == -1:
            step,Nmc,Npseudos = self.getTunedParameters();
        else:
            step = self.getStepLength()*2.0;
            Nmc = self.Nmc;

        #Parallel calculation if numCores > 1
        if self.numCores == 1:
            mcb = self.prepareBayesCalculator(Nmc);
            mcb.confidenceLevel = CL;
            exp = mcb.expectedExclusion(step,Npseudos,False);
            return exp[0], exp[1], exp[2], exp[3], exp[4];
        else:
            return self.getBayesianExpectedLimitParallel(Npseudos, CL, step, Nmc);

    #Calculates the expected Bayesian limit as well as the 1-sigma and 2-sigma fluctuations
    #(i.e. green and yellow bands) using parallel processing.
    def getBayesianExpectedLimitParallel(self, Npseudos, CL, step, Nmc):

        #Number of pseudo-experiments to run in each process
        NpseudosLocal = max(1,int(Npseudos/self.numCores));

        #Each process runs a certain number of pseudo-experiments and returns all the limits in "outputDict"
        procs = [];
        outputDict = mp.Manager().dict();
        for i in range(self.numCores):
            p = mp.Process(target=self.doBayesianExpectedLimitPartial, args=(self.prepareBayesCalculator(Nmc), NpseudosLocal, CL, step, outputDict, i));
            procs.append(p);
            p.start();

        for p in procs:
            p.join();

        #All processes are done, collect all the limits in one list
        allExclusionVals = [];
        for i in range(NpseudosLocal):
            for j in range(self.numCores):
                allExclusionVals.append(outputDict[j][i]);
        allExclusionVals.sort();

        oneSigmaProb = 0.5*(1 + ROOT.TMath.Erf(-1.0/ROOT.TMath.Sqrt(2)));
        twoSigmaProb = 0.5*(1 + ROOT.TMath.Erf(-2.0/ROOT.TMath.Sqrt(2)));

        #Extract median and quantiles from combined list
        m2sig = allExclusionVals[int(twoSigmaProb*len(allExclusionVals))];
        m1sig = allExclusionVals[int(oneSigmaProb*len(allExclusionVals))];
        expected = allExclusionVals[int(0.5*len(allExclusionVals))];
        p1sig = allExclusionVals[int((1.0-oneSigmaProb)*len(allExclusionVals))];
        p2sig = allExclusionVals[int((1.0-twoSigmaProb)*len(allExclusionVals))];

        return m2sig,m1sig,expected,p1sig,p2sig;

    #Function for running the pseudo-experiments for one of the processes
    def doBayesianExpectedLimitPartial(self, mcb, Npseudos, CL, step, outputDict, procNum):

        mcb.confidenceLevel = CL;
        pseudos = mcb.expectedExclusion(step,Npseudos,False,False,True);
        
        pseudosList = [];
        for i in range(Npseudos):
            pseudosList.append(pseudos[i]);
        outputDict[procNum] = pseudosList;



"""
Class for a counting experiment for the combination of an arbitrary number
of channels with correlation between the background levels in the channels.
For each channel, the background uncertainty is split into a correlated and
an uncorrelated component. Only a few functions from the "countingExperiment"
base class are reimplemented, notably "prepareArrays" that sets the values of 
both the correlated and uncorrelated background uncertainties.
"""
class countingExperimentWithBkgCorr(countingExperiment):

    def bkgString(self, chan):
        return "   Background = "+str(chan.bkg)+" +/- "+str(chan.bkgUncUncorr)+"(uncorr) +/- "+str(chan.bkgUncCorr)+"(corr)\n";

    def addChannel(self, name, bkg, bkgUncUncorr, bkgUncCorr, Nobs, eff=1.0, effUnc=0.0):

        self.channels.append(channel(name,bkg,-999.0,Nobs,eff,effUnc,bkgUncUncorr,bkgUncCorr));

    def prepareArrays(self):

        self.background = array('d',[]);
        self.backgroundRelUncUncorr = array('d',[]);
        self.backgroundRelUncCorr = array('d',[]);
        self.Nobs = array('i',[]);
        self.efficiency = array('d',[]);
        self.efficiencyUnc = array('d',[]);

        for chan in self.channels:
            self.background.append(chan.bkg/self.intLum);
            self.backgroundRelUncUncorr.append(chan.bkgUncUncorr/chan.bkg);
            self.backgroundRelUncCorr.append(chan.bkgUncCorr/chan.bkg);
            self.Nobs.append(chan.Nobs);
            self.efficiency.append(chan.eff);
            self.efficiencyUnc.append(chan.effUnc);

