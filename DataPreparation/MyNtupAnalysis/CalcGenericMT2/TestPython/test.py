from ROOT import *

gSystem.Load("../CalcGenericMT2/libBinnedLik.so")

Visap = TLorentzVector( -18.1222 , -14.4356 , 0 , 158.653)
Visbp = TLorentzVector( 48.2681 , 38.449 , 0 , 62.513)
METp = TLorentzVector( -30.1459 , -24.0134 , 0 , 74.8509)

mycalc = ComputeMT2(Visap,Visbp,METp,0.,80.)
print "The MT2 value for this event is:", mycalc.Compute()
#Should print: The MT2 value for this event is: 156.952

