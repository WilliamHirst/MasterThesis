#include <iostream>
#include "TLorentzVector.h"
#include "../CalcGenericMT2/MT2_ROOT.h"

void test() {

  myversion();

  TLorentzVector visa = TLorentzVector( -18.1222 , -14.4356 , 0 , 158.653);
  TLorentzVector visb = TLorentzVector( 48.2681 , 38.449 , 0 , 62.513);
  TLorentzVector met = TLorentzVector( -30.1459 , -24.0134 , 0 , 74.8509);

  ComputeMT2 mycalc = ComputeMT2(visa,visb,met,0.,80.);
  std::cout << "The MT2 value for this event is: " <<  mycalc.Compute() << std::endl;
  //Should print: The MT2 value for this event is: 156.952

}
