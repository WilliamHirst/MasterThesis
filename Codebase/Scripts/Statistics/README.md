# FYS5555 Statistics

Prepared by Magnar Bugge

This section allows the computation of "cut-and-count" ("single-bin") cross section limits on a new signal
process and production of standard limit plots including bands for the 1- and 2-sigma variations of the
limit around the expected (median) value. The significance of an excess can be calculated in the asymptotic
approximation.

The frontend classes that should be used by the user are defined in "pyStats.py", while the actual calculations
happen in the C++ backend code "MC_Bayes_withBkgCorr.cpp". There are two frontend classes, called "countingExperiment"
and "countingExperimentWithBkgCorr", where the latter allows for the treatment of a correlated systematic uncertainty
between all the channels. Examples of usage can be found in the Python code "testPyStats.py" and in the Jupyter
notebook "software/Notebooks/Statistics/statisticsNotebook.ipynb". For example, the combined limit for the 300 GeV mass point of the 8 TeV
ATLAS W' search ([arXiv:1407.7494](https://arxiv.org/abs/1407.7494)) can be reproduced as follows:

```
import pyStats;

countexp = pyStats.countingExperiment(name = '300 GeV', intLum = 20.28, intLumUnc = 0.5678);
countexp.addChannel('electron channel', bkg = 12901.4, bkgUnc = 816.767129814, Nobs = 12717, eff = 0.22758, effUnc = 0.00930046154813);
countexp.addChannel('muon channel', bkg = 11266.4, bkgUnc = 773.744616417, Nobs = 10927, eff = 0.183817, effUnc = 0.00706668179581);

print(countexp);

print( "Observed limit = ",countexp.getBayesianLimit() );

m2sig,m1sig,expected,p1sig,p2sig = countexp.getBayesianExpectedLimit();
print( "\nExpected limit and bands:" );
print( "  -2sigma                 -1sigma                 median                +1sigma                +2sigma" );
print( m2sig, "   ", m1sig, "   ", expected, "   ", p1sig, "   ", p2sig );
```

Keyword arguments are used here to make the code easy to read and understand. In the constructor of the counting experiment,
the integrated luminosity and its uncertainty are specified. For each channel, the following information is specified:
- Expected number of background events and its uncertainty
- Number of observed events
- Signal efficiency and its uncertainty

The line "print(countexp)" gives a clear overview of the information that has been entered. The observed limit is obtained from
the function "getBayesianLimit()" of the counting experiment object, and the expected limit with its 1- and 2-sigma expected
fluctuations is obtained from the function "getBayesianExpectedLimit()".

The significance of an excess in data is obtained from the function "getSignificance()", as in the following example

```
import pyStats;

countexp = pyStats.countingExperiment();
countexp.addChannel('myChannel', bkg = 100.0, bkgUnc = 10.0, Nobs = 120);
print(countexp);
print( "Significance = ",countexp.getSignificance() );
```

The significance is evaluated using the "asymptotic approximation", where the profile likelihood ratio test statistic is
assumed to follow a chi-squared distribution. This is generally a very good approximation for the so-called "q_0" test
statistic used for the significance evaluation, and the results are typically accurate also for low background expectations,
where simple expressions such as s/sqrt(b) are inappropriate.

A separate class is dedicated to the creation of standard limit plots, and the script "makeLimitPlots.py" can be used
to create a set of limit plots based on an input file "inputs.txt" giving the background levels, signal efficiencies etc.
for each point along the x-axis. The input file provided in the git repository is an example that contains the inputs
to reproduce the limit plots from the 8 TeV ATLAS W' search ([arXiv:1407.7494](https://arxiv.org/abs/1407.7494)).
Running "makeLimitPlots.py" "out of the box" will reproduce the plots from Fig. 2 (left) of that publication.
(The combined limit plot will not look identical to the published result, because we do not take into account
the partial correlation of the systematic uncertainties between the electron and muon channels.)
Note that the calculation of the combined limit can take some time depending on the available CPU resources, and the
script contains some by-default-commented lines to skip the combined limit calculation for quick tests.

The file "inputs.txt" consists of Python variable definitions, and should be reasonably self-explanatory. The first
line defines the total number of mass points in the plot, the number of channels, and the integrated luminosity and
its uncertainty. Then, for each mass point, there is
- one line defining the mass (i.e. x-axis value) and the theory cross section
- one line per channel giving the counting experiment inputs for that channel

Note that in the line containing the mass and theory cross section, we additionally kept track of the mT threshold
that was used to define the single-bin signal region used for that mass, although this information is not needed
to produce the limit plot. Finally, the last line of the inputs file specifies the x- and y-axis labels for the
plotting, along with the y-axis range. (When running interactively, you can anyway adjust the axis ranges with
the mouse directly.)

The basic entry of the limitPlot object is a "point", defined in terms of
- An x-axis value (the W' mass in the example inputs)
- A counting experiment
- A theoretical signal prediction (optional)

The limit plot object manages the calculation of the limits for all mass points (using parallel processing if
multiple cores are available), and takes care of the actual plotting of the results.

A basic example of how to use the limitPlot class can be found in the script "wprimeModelIndependentLimitsExample.py",
which reproduces the high-mass part of Fig. 4 (bottom) from the 13 TeV ATLAS W' search
([arXiv:1906.05609](https://arxiv.org/abs/1906.05609)). In this example, instead of using an input file, we manually
enter the inputs from Table 5 of the publication to make the code easy to understand, which is why only the high-mass
part of the plot is included to limit the number of points. Since this example concerns model independent limits,
no signal efficiency is specified. The results are thus quoted in terms of the number of expected signal events in
the signal region (with mT > mTmin) per unit of integrated luminosity, which is referred to as the "visible cross section".
To reproduce the results in terms of number of expected signal events in the full 139/fb, also given in Table 5,
you can set the integrated luminosity to 1 and the uncertainty to 0.

## Advanced usage

In the examples above, the user specifies only a minimum of the required information to define the counting
experiments. Behind the scenes, the code then auto-adjusts or sets default values for certain parameters related to
e.g. the numerical calculations. The auto-adjusted and default values are chosen in such a way that things should
"just work" in most cases. However, the user may encounter special cases where the default parameters are not
appropriate, in which case manual tuning may be necessary.

The Bayesian limit calculation performs the marginalization integral on a regular grid of points along the signal
cross section axis. The corresponding "step length" can be adjusted, and should as a rule of thumb be a few hundred
times smaller than the cross section limit. Additionally, the number of MC cycles performed for the integration
at each cross section point can be adjusted. These parameters are adjusted as in this example:

```
countexp.Nmc = 5000;
countexp.stepLength = 1.0;
```

For the calculation of an expected limit, the number of pseudo-experiments can be specified as an argument:

```
countexp.getBayesianExpectedLimit(Npseudos = 200);
```

When the number of pseudo-experiments is specified like this, the code will also use the user-set values
for "Nmc" and "stepLength", although multiplying the latter by 2, as a bit less precision and faster
calculation is generally a good idea for the pseudo-experiments. When the number of pseudo-experiments
is not specified, "pre-tuned" values are used for this as well as for "Nmc" and "stepLength". Note that
this is only relevant when combining multiple channels, since the calculation for a single channel
uses a "trick" so that the calculation runs very fast.

For the significance calculation, an upper bound on the signal cross section is used in the signal+background fit.
In case some tuning is needed, the upper bound can be specified as an argument to the significance function:

```
countexp.getSignificance(maxSignal = 30.0)
```

When running expected limits for a combination of multiple channels, the pseudo-experiments are by default
parallelized over all available CPUs. You can manually specify the number of jobs that run in parallel by
setting the "numCores" variable:

```
countexp.numCores = 4;
```

When calculating limits via a limitPlot object, the calculation is parallelized both by running multiple
mass points in parallel and by splitting over pseudo-experiments. By default, all mass points are run
in parallel, with additional splitting over pseudo-experiments if there are enough CPUs available.
Alternatively, the number of parallel mass points and number of pseudo-experiment jobs per mass point
can be specified as arguments to the "calculate" function:

```
limPlot.calculate(massPointsBatchSize = 4, NsplitPseudos = 2);
```

All the above configurations are of a purely technical nature. In addition there are some user options
which are more related to the statistical analysis. For example, the user can choose between log-normal (default)
or Gaussian priors / constraints for the background level, signal efficiency, and integrated luminosity:

```
countexp.setPriors(bkg = 'gauss', eff = 'lognormal', lumi = 'lognormal');
```

Note that setting Gaussian priors for the signal efficiency and/or integrated luminosity can lead to problems
in the Bayesian limit calculation. Finally, the confidence level (or "credibility level") can be specified
in the limit calculations, for example, to run with 90% CL:

```
countexp.getBayesianLimit(CL = 0.90)
countexp.getBayesianExpectedLimit(CL = 0.90)
```
