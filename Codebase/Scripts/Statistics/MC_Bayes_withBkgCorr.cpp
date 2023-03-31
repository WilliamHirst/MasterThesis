#include <TRandom3.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <TF1.h>
#include <TMinuit.h>


void* globalMCBayesPointer;
void MinusFrequentistLogLikelihoodWrapper(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);


class MC_Bayes_withBkgCorr {

public:
 
  vector<Int_t> exclusionsI;
  vector<Double_t> exclusionVals;

  vector<Int_t> observationsI;
  vector<Int_t> observationVals;
  Int_t medianObservation;

  vector<Double_t> likelihoods;
  vector<Double_t> signals;

  Double_t intLum;
  Double_t intLumSigma;
  Double_t* background;
  Double_t* backgroundRelSigmaUncorr;
  Double_t* backgroundRelSigmaCorr;
  Double_t limit;
  Int_t* Nobs;
  Int_t Nchan;
  Double_t* signalEff;
  Double_t* signalEffSigma;

  Double_t signalEffGaussTrunc;
  Double_t signalEffUpperTrunc;

  Double_t twoSigmaProb;
  Double_t oneSigmaProb;

  Double_t confidenceLevel;

  Bool_t lowstats;
  Bool_t doRefPrior;
  Bool_t doLogNormalEff;
  Bool_t doLogNormalLumi;
  Bool_t doLogNormalBkg;

  Double_t* logNormalMuEff;
  Double_t* logNormalSigmaEff;
  Double_t logNormalMuLumi;
  Double_t logNormalSigmaLumi;
  Double_t* logNormalMuBkg;
  Double_t* logNormalSigmaBkgUncorr;
  Double_t* logNormalSigmaBkgCorr;

  Bool_t doGammaEff;
  TF1** gammaEff;

  Bool_t doModGaussEff;
  TF1** modGaussEff;
  Double_t modGaussPower;

  Double_t backgroundRelSigmaCorrMax;

  ofstream* log;

  Int_t Nmc;

  TRandom* rnd;

  MC_Bayes_withBkgCorr(Int_t NchanA, Double_t intLumA, Double_t intLumSigmaA, Double_t* backgroundA, Double_t* backgroundRelSigmaUncorrA, 
		       Int_t* NobsA, Int_t NmcA, Double_t* signalEffA, Double_t* signalEffSigmaA, Double_t* backgroundRelSigmaCorrA, 
		       Bool_t doRefPriorA, Bool_t doLogNormalEffA, Bool_t doLogNormalLumiA, Bool_t doLogNormalBkgA,
		       Int_t logNormalMode = 2, Bool_t doGammaEffA = false, Int_t gammaMode = 1, Double_t signalEffGaussTruncA = 0.0,
		       Bool_t doModGaussEffA = false, Double_t modGaussPowerA = 1.0, Double_t signalEffUpperTruncA = 1.0) {

    intLum = intLumA;
    intLumSigma = intLumSigmaA;
    background = backgroundA;
    backgroundRelSigmaUncorr = backgroundRelSigmaUncorrA;
    backgroundRelSigmaCorr = backgroundRelSigmaCorrA;
    Nobs = NobsA;
    Nmc = NmcA;
    Nchan = NchanA;
    signalEff = signalEffA;
    signalEffSigma = signalEffSigmaA;
    doGammaEff = doGammaEffA;
    signalEffGaussTrunc = signalEffGaussTruncA;
    signalEffUpperTrunc = signalEffUpperTruncA;

    doRefPrior = doRefPriorA;
    doLogNormalEff = doLogNormalEffA;
    doLogNormalLumi = doLogNormalLumiA;
    doLogNormalBkg = doLogNormalBkgA;

    doModGaussEff = doModGaussEffA;
    modGaussPower = modGaussPowerA;

    logNormalMuEff = new Double_t[Nchan];
    logNormalSigmaEff = new Double_t[Nchan];
    logNormalMuBkg = new Double_t[Nchan];
    logNormalSigmaBkgUncorr = new Double_t[Nchan];
    logNormalSigmaBkgCorr = new Double_t[Nchan];

    if (logNormalMode == 1) //mean
      logNormalMuLumi = TMath::Log(intLum/TMath::Sqrt(1.0 + TMath::Power(intLumSigma/intLum,2.0)));
    else if (logNormalMode == 2) //median
      logNormalMuLumi = TMath::Log(intLum);
    else if (logNormalMode == 3) //mode
      logNormalMuLumi = TMath::Log(intLum*(1.0 + TMath::Power(intLumSigma/intLum,2.0)));
    else {
      cout << "Unrecognized log-normal mode\n";
      exit(1);
    }

    logNormalSigmaLumi = TMath::Sqrt(TMath::Log(1.0 + TMath::Power(intLumSigma/intLum,2.0)));

    for (Int_t chan = 0; chan < Nchan; chan++) {

      if (logNormalMode == 1) //mean
	logNormalMuEff[chan] = TMath::Log(signalEff[chan]/TMath::Sqrt(1.0 + TMath::Power(signalEffSigma[chan]/signalEff[chan],2.0)));
      else if (logNormalMode == 2) //median
	logNormalMuEff[chan] = TMath::Log(signalEff[chan]);
      else if (logNormalMode == 3) //mode
	logNormalMuEff[chan] = TMath::Log(signalEff[chan]*(1.0 + TMath::Power(signalEffSigma[chan]/signalEff[chan],2.0)));
      else {
	cout << "Unrecognized log-normal mode\n";
	exit(1);
      }

      logNormalSigmaEff[chan] = TMath::Sqrt(TMath::Log(1.0 + TMath::Power(signalEffSigma[chan]/signalEff[chan],2.0)));

      Double_t backgroundRelSigmaTotSq = backgroundRelSigmaCorr[chan]*backgroundRelSigmaCorr[chan] + backgroundRelSigmaUncorr[chan]*backgroundRelSigmaUncorr[chan];
      
      if (logNormalMode == 1) //mean
	logNormalMuBkg[chan] = TMath::Log(background[chan]/TMath::Sqrt(1.0 + backgroundRelSigmaTotSq));
      else if (logNormalMode == 2) //median
	logNormalMuBkg[chan] = TMath::Log(background[chan]);
      else if (logNormalMode == 3) //mode
	logNormalMuBkg[chan] = TMath::Log(background[chan]*(1.0 + backgroundRelSigmaTotSq));
      else {
	cout << "Unrecognized log-normal mode\n";
	exit(1);
      }

      logNormalSigmaBkgUncorr[chan] = TMath::Sqrt(TMath::Log(1.0 + TMath::Power(backgroundRelSigmaUncorr[chan],2.0)));
      logNormalSigmaBkgCorr[chan] = TMath::Sqrt(TMath::Log(1.0 + TMath::Power(backgroundRelSigmaCorr[chan],2.0)));
    }


    gammaEff = new TF1*[Nchan];

    if (doGammaEff + doLogNormalEff + doModGaussEff > 1) {

      cout << "Inconsistent initialization (more than one prior selected for the efficiency)\n";
      exit(1);
    }


    for (Int_t chan = 0; chan < Nchan; chan++) {

      Double_t gammaEffK = TMath::Power(signalEff[chan]/signalEffSigma[chan],2.0);
      Double_t gammaEffTheta;
      if (gammaMode == 1)  //mean
	gammaEffTheta = signalEff[chan]/gammaEffK;
      else if (gammaMode == 3)  //mode
	gammaEffTheta = signalEff[chan]/(gammaEffK - 1.0);
      else {

	cout << "Unrecognized gamma mode\n";
	exit(1);
      }

      TString name;
      if (chan == 0)
	name = "gammaEff1";
      else if (chan == 1)
	name = "gammaEff2";
      else
	name = "gammaEff";
      
      gammaEff[chan] = new TF1(name,"TMath::Exp(([0]-1)*TMath::Log(x) - x/[1] - [0]*TMath::Log([1]) - TMath::LnGamma([0]))",0.0,1.0);
      gammaEff[chan]->SetNpx(10000);
      gammaEff[chan]->SetParameter(0,gammaEffK);
      gammaEff[chan]->SetParameter(1,gammaEffTheta);
    }


    modGaussEff = new TF1*[Nchan];


    for (Int_t chan = 0; chan < Nchan; chan++) {

      TString name;
      if (chan == 0)
	name = "modGaussEff1";
      else if (chan == 1)
	name = "modGaussEff2";
      else
	name = "modGaussEff";
      
      modGaussEff[chan] = new TF1(name,"TMath::Power(x*(2.0*[0]-x),[1])*TMath::Exp(-(x-[0])*(x-[0])/(2.0*[2]*[2]))",0.0,2.0*signalEff[chan]);
      modGaussEff[chan]->SetNpx(10000);
      modGaussEff[chan]->SetParameter(0,signalEff[chan]);
      modGaussEff[chan]->SetParameter(1,modGaussPower);
      modGaussEff[chan]->SetParameter(2,signalEffSigma[chan]);
    }


    limit = -999.0;

    confidenceLevel = 0.95;
    
    rnd = new TRandom3(0);

    oneSigmaProb = 0.5*(1 + TMath::Erf(-1.0/TMath::Sqrt(2)));
    twoSigmaProb = 0.5*(1 + TMath::Erf(-2.0/TMath::Sqrt(2)));

    log = new ofstream("MC_Bayes_logfile.txt", ios::app);

    //We drop the factorial in the Poisson likelihood when statistics is low,
    //since we normalize afterwards in any case
    lowstats = true;
    for (Int_t chan = 0; chan < Nchan; chan++)
      if (Nobs[chan] > 10)
	lowstats = false;


    backgroundRelSigmaCorrMax = 0.0;
    for (Int_t chan = 0; chan < Nchan; chan++)
      if (backgroundRelSigmaCorr[chan] > backgroundRelSigmaCorrMax)
	backgroundRelSigmaCorrMax = backgroundRelSigmaCorr[chan];
  }


  ~MC_Bayes_withBkgCorr() {
    
    delete rnd;

    delete logNormalMuEff;
    delete logNormalSigmaEff;
    delete logNormalMuBkg;
    delete logNormalSigmaBkgUncorr;
    delete logNormalSigmaBkgCorr;

    log->close();
    delete log;

    for (Int_t chan = 0; chan < Nchan; chan++)
      delete gammaEff[chan];
    delete gammaEff;

    for (Int_t chan = 0; chan < Nchan; chan++)
      delete modGaussEff[chan];
    delete modGaussEff;
  }



  Bool_t notsortedExcl() {

    for (UInt_t i = 0; i < exclusionVals.size()-1; i++)
      if (exclusionVals[i+1] < exclusionVals[i])
	return true;

    return false;
  }


  void sortExcl() {

    Double_t temp;
    Int_t temp2;

    while (notsortedExcl())
      for (UInt_t i = 0; i < exclusionVals.size()-1; i++) {
	
	if (exclusionVals[i+1] < exclusionVals[i]) {
	  temp = exclusionVals[i+1];
	  exclusionVals[i+1] = exclusionVals[i];
	  exclusionVals[i] = temp;
	  
	  temp2 = exclusionsI[i+1];
	  exclusionsI[i+1] = exclusionsI[i];
	  exclusionsI[i] = temp2;

	}
      }
  }




  Int_t getIndexExcl(Double_t exclusion) {

    for (UInt_t N = 0; N < exclusionVals.size(); N++)
      if (exclusionVals[N] == exclusion)
	return N;

    exclusionVals.push_back(exclusion);
    exclusionsI.push_back(0);

    return exclusionVals.size() - 1;
  }



  Bool_t notsortedObs() {

    for (UInt_t i = 0; i < observationVals.size()-1; i++)
      if (observationVals[i+1] < observationVals[i])
	return true;
    
    return false;
  }


  void sortObservations() {

    Int_t temp;
    Int_t temp2;
    
    while (notsortedObs())
      for (UInt_t i = 0; i < observationVals.size()-1; i++) {
	
	if (observationVals[i+1] < observationVals[i]) {
	  temp = observationVals[i+1];
	  observationVals[i+1] = observationVals[i];
	  observationVals[i] = temp;
	  
	  temp2 = observationsI[i+1];
	  observationsI[i+1] = observationsI[i];
	  observationsI[i] = temp2;
	  
	}
      }
  }


  Int_t getIndexObs(Int_t observation) {
    
    for (UInt_t N = 0; N < observationVals.size(); N++)
      if (observationVals[N] == observation)
	return N;
    
    observationVals.push_back(observation);
    observationsI.push_back(0);
    
    return observationVals.size() - 1;
  }

  
  Double_t q0TestStatistic(Double_t maxSignal, Bool_t verbose = true) {

    TMinuit *minuit = new TMinuit(3 + 2*Nchan);
    minuit->SetPrintLevel(-1);
    
    globalMCBayesPointer = this;
    minuit->SetFCN(MinusFrequentistLogLikelihoodWrapper);

    Double_t arglist[2];
    Int_t ierflg = 0;
    vector<string> parNames;
    parNames.reserve(3 + 2*Nchan);
    double dummy;
    int dummyI;


    minuit->mnparm(0, "signal", 0.0, maxSignal/100000.0, 0.0, maxSignal, ierflg);
    parNames.push_back("signal");

    //fix signal to 0 for background-only maximization
    minuit->FixParameter(0);

    minuit->mnparm(1, "thetaLumi", 0.0, 0.001, -5.0, 5.0, ierflg);
    parNames.push_back("thetaLumi");
    minuit->mnparm(2, "thetaBkgCorr", 0.0, 0.001, -5.0, 5.0, ierflg);
    parNames.push_back("thetaBkgCorr");

    Int_t idx = 3;
    for (Int_t chan = 0; chan < Nchan; chan++) {

      minuit->mnparm(idx++, ("thetaBkgUncorr"+to_string(chan)).c_str(), 0.0, 0.001, -5.0, 5.0, ierflg);
      parNames.push_back("thetaBkgUncorr_channel"+to_string(chan));
    }

    for (Int_t chan = 0; chan < Nchan; chan++) {

      minuit->mnparm(idx++, ("thetaEff"+to_string(chan)).c_str(), 0.0, 0.001, -5.0, 5.0, ierflg);
      parNames.push_back("thetaEff_channel"+to_string(chan));
    }

    arglist[0] = 10000; //max number of calls to user function
    arglist[1] = .0001; //tolerance(?)

    minuit->mnexcm("MINIMIZE", arglist,2,ierflg); //background-only fit

    if (ierflg != 0) {
      cout << "WRONG RETURN VALUE from minimization (background-only)!\n";
      exit(1);
    }

    if (verbose) {
      cout << "--------------------------------------\nResults of background-only fit: \n------------------------------------\n";
      for (unsigned int i = 1; i < parNames.size(); i++) {
	
	double value,error;
	minuit->GetParameter(i,value,error);
	
	cout << parNames[i] << "  " << value << endl;
      }
    }

    double maxLogLikelihoodBonly;
    minuit->mnstat(maxLogLikelihoodBonly,dummy,dummy,dummyI,dummyI,dummyI);
    maxLogLikelihoodBonly = -maxLogLikelihoodBonly;
    
    if (verbose) cout << "Maximum log-likelihood (background-only) = " << maxLogLikelihoodBonly << endl;

    //Release signal parameter for signal+background fit
    minuit->Release(0);
    minuit->mnexcm("MINIMIZE", arglist,2,ierflg); //signal+background fit

    if (ierflg != 0) {
      cout << "WRONG RETURN VALUE from minimization (signal+background)!\n";
      exit(1);
    }

    double signal,sigErr;
    minuit->GetParameter(0,signal,sigErr);

    if (verbose) {
      cout << "-------------------------------------\nResults of signal+background fit: \n------------------------------------\n";
      for (unsigned int i = 0; i < parNames.size(); i++) {
	
	double value,error;
	minuit->GetParameter(i,value,error);
	
	cout << parNames[i] << "  " << value << endl;
      }
    }

    if (signal > 0.9*maxSignal) {
      cout << "ERROR: Too low maxSignal\n";
      exit(1);
    }

    double maxLogLikelihoodSplusB;
    minuit->mnstat(maxLogLikelihoodSplusB,dummy,dummy,dummyI,dummyI,dummyI);
    maxLogLikelihoodSplusB = -maxLogLikelihoodSplusB;

    if (verbose) cout << "Maximum log-likelihood (signal+background) = " << maxLogLikelihoodSplusB << endl;

    //Test statistic
    double q0 = -2.0 * (maxLogLikelihoodBonly - maxLogLikelihoodSplusB);

    //For downwards fluctuations, q0 is sometimes numerically negative - set to 0
    if (q0 < 0.0)
      q0 = 0.0;

    delete minuit;
    
    return q0;
  }

  
  Double_t FrequentistLogLikelihood(Double_t signal, Double_t thetaLumi, Double_t thetaBkgCorr, Double_t* thetaBkg, Double_t* thetaEff) {

    Double_t thisLogLikelihood;
    Double_t thisIntLum,thisBackground,s,b,lambda,thisSignalEff,bkgCorrFact,bkgUncorrFact;

    if (doLogNormalLumi)
      thisIntLum = TMath::Exp(logNormalMuLumi + thetaLumi*logNormalSigmaLumi);
    else
      thisIntLum = intLum + thetaLumi*intLumSigma;

    if (thisIntLum <= 0.0)
      thisIntLum = 1.0e-100;

    bkgCorrFact = thetaBkgCorr;

    thisLogLikelihood = 0.0;
    
    //Loop over channels
    for (Int_t chan = 0; chan < Nchan; chan++) {

      if (doLogNormalBkg) {
	thisBackground = TMath::Exp(logNormalMuBkg[chan] + bkgCorrFact*logNormalSigmaBkgCorr[chan] + thetaBkg[chan]*logNormalSigmaBkgUncorr[chan]);
      }
      else {
	thisBackground = background[chan]*(1.0 + bkgCorrFact*backgroundRelSigmaCorr[chan] + thetaBkg[chan]*backgroundRelSigmaUncorr[chan]);
      }

      if (thisBackground <= 0.0)
	thisBackground = 1.0e-100;

      if (doLogNormalEff) {
	thisSignalEff = TMath::Exp(logNormalMuEff[chan] + thetaEff[chan]*logNormalSigmaEff[chan]);
	//The log-normal expression is not well defined for 0 efficiency
	if (signalEff[chan] == 0.0 && signalEffSigma[chan] == 0.0)
	  thisSignalEff = 0.0;
      }
      else if (doGammaEff) {
	cout << "ERROR: Gamma efficiency prior not implemented for frequentist calculation\n";
	exit(1);
      }
      else if (doModGaussEff) {
	cout << "ERROR: Modified Gaussian efficiency prior not implemented for frequentist calculation\n";
	exit(1);
      }
      else {
	thisSignalEff = signalEff[chan] + thetaEff[chan]*signalEffSigma[chan];
      }

      if (thisSignalEff <= 0.0)
	thisSignalEff = 1.0e-100;


      //expected signal and background
      s = thisIntLum*signal*thisSignalEff;
      b = thisIntLum*thisBackground;
      lambda = s+b;
      
      if (!lowstats)
	thisLogLikelihood += Nobs[chan]*TMath::Log(lambda) - lambda - TMath::LnGamma(Nobs[chan]+1.0);
      else
	thisLogLikelihood += Nobs[chan]*TMath::Log(lambda) - lambda;
    }


    //Constraint terms for systematics
    thisLogLikelihood -= TMath::Power(thetaLumi,2) / 2.0;
    thisLogLikelihood -= TMath::Power(thetaBkgCorr,2) / 2.0;
    for (Int_t chan = 0; chan < Nchan; chan++) {

      thisLogLikelihood -= TMath::Power(thetaEff[chan],2) / 2.0;
      thisLogLikelihood -= TMath::Power(thetaBkg[chan],2) / 2.0;
    }

    return thisLogLikelihood;
  }

  Double_t CalcLikelihood(Double_t signal) {

    
    Double_t thisLogLikelihood;
    Double_t thisIntLum,thisBackground,s,b,lambda,thisSignalEff,bkgCorrFact,bkgUncorrFact;
    Double_t sum = 0.0;


    //Monte Carlo loop
    for (Int_t mc = 0; mc < Nmc; mc++) {
      
      if (intLumSigma == 0.0)
	thisIntLum = intLum;
      else {
	if (doLogNormalLumi)
	  thisIntLum = TMath::Exp(rnd->Gaus(logNormalMuLumi,logNormalSigmaLumi));
	else {
	  do
	    thisIntLum = rnd->Gaus(intLum,intLumSigma);
	  while (thisIntLum < 0.0);
	}
      }

      if (backgroundRelSigmaCorrMax != 0.0) {
	if (doLogNormalBkg) {
	  bkgCorrFact = rnd->Gaus(0.0,1.0);
	}
	else {
	  do 
	    bkgCorrFact = rnd->Gaus(0.0,1.0);
	  while (bkgCorrFact < -1.0/backgroundRelSigmaCorrMax);
	}
      }
      else
	bkgCorrFact = 0.0;

      thisLogLikelihood = 0.0;

      //Loop over channels
      for (Int_t chan = 0; chan < Nchan; chan++) {

	if (backgroundRelSigmaUncorr[chan] == 0.0) {
	  if (doLogNormalBkg)
	    thisBackground = TMath::Exp(logNormalMuBkg[chan] + bkgCorrFact*logNormalSigmaBkgCorr[chan]);
	  else
	    thisBackground = background[chan]*(1.0 + bkgCorrFact*backgroundRelSigmaCorr[chan]);
	}
	else {
	  if (doLogNormalBkg) {
	    bkgUncorrFact = rnd->Gaus(0.0,1.0);
	    thisBackground = TMath::Exp(logNormalMuBkg[chan] + bkgCorrFact*logNormalSigmaBkgCorr[chan] + bkgUncorrFact*logNormalSigmaBkgUncorr[chan]);
	  }
	  else {
	    do {
	      bkgUncorrFact = rnd->Gaus(0.0,1.0);
	      thisBackground = background[chan]*(1.0 + bkgCorrFact*backgroundRelSigmaCorr[chan] + bkgUncorrFact*backgroundRelSigmaUncorr[chan]);
	    }
	    while (thisBackground < 0.0);
	  }
	}
	
	if (signalEffSigma[chan] == 0.0)
	  thisSignalEff = signalEff[chan];
	else {
	  if (doLogNormalEff) {
	    do 
	      thisSignalEff = TMath::Exp(rnd->Gaus(logNormalMuEff[chan],logNormalSigmaEff[chan]));
	    while (thisSignalEff > signalEffUpperTrunc);
	  }
	  else if (doGammaEff) {
	    thisSignalEff = gammaEff[chan]->GetRandom(0.0,1.0);
	  }
	  else if (doModGaussEff) {
	    thisSignalEff = modGaussEff[chan]->GetRandom(0.0,2.0*signalEff[chan]);
	  }
	  else {
	    do
	      thisSignalEff = rnd->Gaus(signalEff[chan],signalEffSigma[chan]);
	    while (thisSignalEff < signalEffGaussTrunc || thisSignalEff > signalEffUpperTrunc);
	  }
	}
	
	//expected signal and background
	s = thisIntLum*signal*thisSignalEff;
	b = thisIntLum*thisBackground;
	lambda = s+b;
	
	if (!lowstats)
	  thisLogLikelihood += Nobs[chan]*TMath::Log(lambda) - lambda - TMath::LnGamma(Nobs[chan]+1.0);
	else
	  thisLogLikelihood += Nobs[chan]*TMath::Log(lambda) - lambda;

      }

      
      sum += TMath::Exp(thisLogLikelihood);
    }

    return sum/Nmc;
  }



  Double_t excludedSignal(Double_t h = 0.001, Bool_t verbose = true) {

    likelihoods.clear();
    signals.clear();

    Double_t signal = 0.0;
    Double_t sum = 0.0;
    Double_t likelihood;
    Double_t prior;

    //Calculate the posterior probability distribution
    do {

      signal += h;
      signals.push_back(signal);

      if (doRefPrior) {
	prior = 0.0;
	
	for (Int_t chan = 0; chan < Nchan; chan++)
	  prior += signalEff[chan]*signalEff[chan]*intLum / (signalEff[chan]*signal + background[chan]);
	
	prior = TMath::Sqrt(prior);
      }
      else
	prior = 1.0;

      likelihood = CalcLikelihood(signal) * prior;  //"likelihood" is really the posterior
      likelihoods.push_back(likelihood);
      sum += likelihood;

    } while (likelihood/(sum - likelihood) > 1.0e-8);

    //Normalize the distribution
    sum = 0.0;
    for (UInt_t i = 0; i < likelihoods.size(); i++) {

      sum += likelihoods[i]*h;
    }

    for (UInt_t i = 0; i < likelihoods.size(); i++)
      likelihoods[i] /= sum;

    //Find exclusion limit
    Double_t tail, prevTail = -999.0;
    for (UInt_t i = 0; i < likelihoods.size(); i++) {
      
      tail = 0.0;
      for (UInt_t j = i; j < likelihoods.size(); j++)
	tail += likelihoods[j]*h;

      if (tail <= (1.0 - confidenceLevel)) {
	if (TMath::Abs(tail - (1.0 - confidenceLevel)) < TMath::Abs(prevTail - (1.0 - confidenceLevel)))
	  limit = signals[i];
	else
	  limit = signals[i-1];
	break;
      }
      
      prevTail = tail;
    }

    if (verbose)
      cout << "Calculated limit = " << limit << endl;
    
    if (limit*2.0 > signals[signals.size()-1])
      *log << "Warning: posterior is cut off below limit*2, limit = " << limit << ", cutoff = " 
	   << signals[signals.size()-1] << ", limit/cutoff = " << limit/signals[signals.size()-1] << endl;

    return limit;
  }



  Double_t* expectedExclusionSingleChan(Double_t h = 0.001, Int_t Nmc2 = 5000, Bool_t verbose = true) {

    observationsI.clear();
    observationVals.clear();

    Int_t keepNobs = Nobs[0];
    Double_t keepLimit = limit;

    //Monte Carlo loop
    Int_t thisNobs;
    Double_t thisIntLum,thisBackground,bkgCorrFact,bkgUncorrFact;
    for (Int_t mc = 0; mc < Nmc2; mc++) {
      
      if (intLumSigma == 0.0)
	thisIntLum = intLum;
      else {
	if (doLogNormalLumi)
	  thisIntLum = TMath::Exp(rnd->Gaus(logNormalMuLumi,logNormalSigmaLumi));
	else {
	  do
	    thisIntLum = rnd->Gaus(intLum,intLumSigma);
	  while (thisIntLum < 0.0);
	}
      }
      

      if (backgroundRelSigmaCorrMax != 0.0) {
	if (doLogNormalBkg) {
	  bkgCorrFact = rnd->Gaus(0.0,1.0);
	}
	else {
	  do 
	    bkgCorrFact = rnd->Gaus(0.0,1.0);
	  while (bkgCorrFact < -1.0/backgroundRelSigmaCorrMax);
	}
      }
      else
	bkgCorrFact = 0.0;


      //No loop over channels
      Int_t chan = 0;
      
      if (backgroundRelSigmaUncorr[chan] == 0.0) {
	if (doLogNormalBkg)
	  thisBackground = TMath::Exp(logNormalMuBkg[chan] + bkgCorrFact*logNormalSigmaBkgCorr[chan]);
	else
	  thisBackground = background[chan]*(1.0 + bkgCorrFact*backgroundRelSigmaCorr[chan]);
      }
      else {
	if (doLogNormalBkg) {
	  bkgUncorrFact = rnd->Gaus(0.0,1.0);
	  thisBackground = TMath::Exp(logNormalMuBkg[chan] + bkgCorrFact*logNormalSigmaBkgCorr[chan] + bkgUncorrFact*logNormalSigmaBkgUncorr[chan]);
	}
	else {
	  do {
	    bkgUncorrFact = rnd->Gaus(0.0,1.0);
	    thisBackground = background[chan]*(1.0 + bkgCorrFact*backgroundRelSigmaCorr[chan] + bkgUncorrFact*backgroundRelSigmaUncorr[chan]);
	  }
	  while (thisBackground < 0.0);
	}
      }
	
      thisNobs = rnd->Poisson(thisIntLum*thisBackground);

      observationsI[getIndexObs(thisNobs)]++; 
    }
    
    sortObservations(); 


    vector<Int_t> allObservationValsSorted;
    allObservationValsSorted.clear();

    for (UInt_t i = 0; i < observationsI.size(); i++) {
      for (Int_t j = 0; j < observationsI[i]; j++)
	allObservationValsSorted.push_back(observationVals[i]);
    }

    Int_t* relevantObservations = new Int_t[5];
 
    relevantObservations[0] = allObservationValsSorted[(int)(twoSigmaProb*allObservationValsSorted.size())];
    relevantObservations[1] = allObservationValsSorted[(int)(oneSigmaProb*allObservationValsSorted.size())];
    relevantObservations[2] = allObservationValsSorted[(int)(0.5*allObservationValsSorted.size())];
    relevantObservations[3] = allObservationValsSorted[(int)((1.0-oneSigmaProb)*allObservationValsSorted.size())];
    relevantObservations[4] = allObservationValsSorted[(int)((1.0-twoSigmaProb)*allObservationValsSorted.size())];

    medianObservation = relevantObservations[2];

    Double_t* returnValue = new Double_t[7];

    for (Int_t i = 0; i < 5; i++) {
      Nobs[0] = relevantObservations[i];
      
      if (Nobs[0] == keepNobs && keepLimit > 0.0)
	returnValue[i] = keepLimit;
      else if (i > 0 && Nobs[0] == relevantObservations[i-1])
	returnValue[i] = returnValue[i-1];
      else
	returnValue[i] = excludedSignal(h, false);
    }

    //Cannot get average and sample std this way
    returnValue[5] = -999.0;
    returnValue[6] = -999.0;

    Nobs[0] = keepNobs;

    delete[] relevantObservations;

    if (verbose)
      std::cout << "Expected limit (xsec): " << returnValue[2] << endl;

    //Print distribution if less than 10 points
    if (verbose && observationVals.size() < 10) {
      cout << "Observations and their probabilities" << endl;
      for (UInt_t i = 0; i < observationVals.size(); i++)
	cout << observationVals[i] <<  " " << observationsI[i]*1.0/Nmc2 << endl;
    }
    
    if (verbose)
      for (Int_t i = 0; i < 5; i++)
	std::cout << "returnValue[" << i << "] = " << returnValue[i] << std::endl;

    return returnValue;
  }



  
  Double_t* expectedExclusion(Double_t h = 0.001, Int_t Nmc2 = 5000, Bool_t verbose = true, Bool_t needAvgAndStd = false, Bool_t returnPEresults = false) {

    if (Nchan == 1 && !needAvgAndStd)
      return expectedExclusionSingleChan(h,Nmc2,verbose);

    
    vector<Int_t*> NobsCalced;
    vector<Double_t> exclusionCalced;
    
    Int_t* NobsV;
    NobsCalced.clear();
    exclusionCalced.clear();
    exclusionsI.clear();
    exclusionVals.clear();

    Int_t* keepNobs = Nobs;

    Bool_t calced;
    Double_t thisIntLum,thisBackground,exclusion,bkgCorrFact,bkgUncorrFact;

    Double_t sum = 0.0;
    Double_t sum2 = 0.0;

    Int_t precycle = 0;

    //Include first already calculated observation
    if (limit > 0.0) {
      NobsV = new Int_t[Nchan];
      for (Int_t i = 0; i < Nchan; i++)
	NobsV[i] = Nobs[i];
    
      NobsCalced.push_back(NobsV);
      exclusionCalced.push_back(limit);

      sum += limit;
      sum2 += limit*limit;

      exclusionsI[getIndexExcl(limit)]++;

      precycle = 1;
    }

    Nobs = new Int_t[Nchan];

    //Monte Carlo loop
    for (Int_t mc = 0; mc < Nmc2-precycle; mc++) {
      
      if (intLumSigma == 0.0)
	thisIntLum = intLum;
      else {
	if (doLogNormalLumi)
	  thisIntLum = TMath::Exp(rnd->Gaus(logNormalMuLumi,logNormalSigmaLumi));
	else {
	  do
	    thisIntLum = rnd->Gaus(intLum,intLumSigma);
	  while (thisIntLum < 0.0);
	}
      }


      if (backgroundRelSigmaCorrMax != 0.0) {
	if (doLogNormalBkg) {
	  bkgCorrFact = rnd->Gaus(0.0,1.0);
	}
	else {
	  do 
	    bkgCorrFact = rnd->Gaus(0.0,1.0);
	  while (bkgCorrFact < -1.0/backgroundRelSigmaCorrMax);
	}
      }
      else
	bkgCorrFact = 0.0;


      //Loop over channels
      for (Int_t chan = 0; chan < Nchan; chan++) {

	if (backgroundRelSigmaUncorr[chan] == 0.0) {
	  if (doLogNormalBkg)
	    thisBackground = TMath::Exp(logNormalMuBkg[chan] + bkgCorrFact*logNormalSigmaBkgCorr[chan]);
	  else
	    thisBackground = background[chan]*(1.0 + bkgCorrFact*backgroundRelSigmaCorr[chan]);
	}
	else {
	  if (doLogNormalBkg) {
	    bkgUncorrFact = rnd->Gaus(0.0,1.0);
	    thisBackground = TMath::Exp(logNormalMuBkg[chan] + bkgCorrFact*logNormalSigmaBkgCorr[chan] + bkgUncorrFact*logNormalSigmaBkgUncorr[chan]);
	  }
	  else {
	    do {
	      bkgUncorrFact = rnd->Gaus(0.0,1.0);
	      thisBackground = background[chan]*(1.0 + bkgCorrFact*backgroundRelSigmaCorr[chan] + bkgUncorrFact*backgroundRelSigmaUncorr[chan]);
	    }
	    while (thisBackground < 0.0);
	  }
	}
	
	Nobs[chan] = rnd->Poisson(thisIntLum*thisBackground);
      }


      exclusion = -999.0;
      
      for (UInt_t i = 0; i < NobsCalced.size(); i++) {
	
	calced = true;
	for (Int_t chan = 0; chan < Nchan; chan++)
	  if (NobsCalced[i][chan] != Nobs[chan])
	    calced = false;

	if (calced)
	  exclusion = exclusionCalced[i];
      }
      
      if (exclusion < 0.0) {
	exclusion = excludedSignal(h, false);
	
	NobsV = new Int_t[Nchan];
	for (Int_t i = 0; i < Nchan; i++)
	  NobsV[i] = Nobs[i];
	
	NobsCalced.push_back(NobsV);
	exclusionCalced.push_back(exclusion);
      }
      
      sum += exclusion;
      sum2 += exclusion*exclusion;

      exclusionsI[getIndexExcl(exclusion)]++; 
    }
    
    //sort exclusion values
    sortExcl(); 

    //free memory
    delete[] Nobs;
    for (UInt_t i = 0; i < NobsCalced.size(); i++)
      delete[] NobsCalced[i];

    //reset original observation
    Nobs = keepNobs;

    //Sample standard deviation
    Double_t sigma = TMath::Sqrt(sum2/Nmc2 - (sum/Nmc2)*(sum/Nmc2));

    Double_t* returnValue = new Double_t[7];

    //Could somtimes be nice to have also the average limit (not only median)
    returnValue[6] = sum/Nmc2;

    vector<Double_t> allExclusionValsSorted;
    allExclusionValsSorted.clear();

    for (UInt_t i = 0; i < exclusionsI.size(); i++) {
      for (Int_t j = 0; j < exclusionsI[i]; j++)
	allExclusionValsSorted.push_back(exclusionVals[i]);
    }

 
    returnValue[0] = allExclusionValsSorted[(int)(twoSigmaProb*allExclusionValsSorted.size())];
    returnValue[1] = allExclusionValsSorted[(int)(oneSigmaProb*allExclusionValsSorted.size())];
    returnValue[2] = allExclusionValsSorted[(int)(0.5*allExclusionValsSorted.size())];
    returnValue[3] = allExclusionValsSorted[(int)((1.0-oneSigmaProb)*allExclusionValsSorted.size())];
    returnValue[4] = allExclusionValsSorted[(int)((1.0-twoSigmaProb)*allExclusionValsSorted.size())];


    if (verbose)
      std::cout << "Expected limit (xsec): " << returnValue[2] << endl; //sum/Nmc2 << " +/- " << sigma << std::endl;

    //Give also sample standard deviation
    returnValue[5] = sigma;
    
    //Print distribution if less than 10 points
    if (verbose && exclusionVals.size() < 10) {
      cout << "Exclusions and their probabilities" << endl;
      for (UInt_t i = 0; i < exclusionVals.size(); i++)
	cout << exclusionVals[i] <<  " " << exclusionsI[i]*1.0/Nmc2 << endl;
    }

    if (verbose)
      for (Int_t i = 0; i < 5; i++)
	std::cout << "returnValue[" << i << "] = " << returnValue[i] << std::endl;

    if (returnPEresults) {
      delete[] returnValue;
      returnValue = new Double_t[allExclusionValsSorted.size()];
      for (UInt_t i = 0; i < allExclusionValsSorted.size(); i++)
 	returnValue[i] = allExclusionValsSorted[i];
      return returnValue;
    }
    else
      return returnValue;
  }
  
  

  TGraph* getPosteriorDistribution() {

    if (limit == -999.0)
      std::cout << "Error: No calculated posterior distribution" << std::endl;

    Double_t* x = new Double_t[likelihoods.size()];
    Double_t* y = new Double_t[likelihoods.size()];

    for (UInt_t i = 0; i < likelihoods.size(); i++) {
      
      x[i] = signals[i]; 
      y[i] = likelihoods[i];
    }

    TGraph* posterior = new TGraph(likelihoods.size(),x,y);

    delete[] x;
    delete[] y;

    return posterior;
  }

};




void MinusFrequentistLogLikelihoodWrapper(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {

  int Nchan = ((MC_Bayes_withBkgCorr*)globalMCBayesPointer)->Nchan;
  
  f = - ((MC_Bayes_withBkgCorr*)globalMCBayesPointer)->FrequentistLogLikelihood(par[0], par[1], par[2], &(par[3]), &(par[3+Nchan]));
}
