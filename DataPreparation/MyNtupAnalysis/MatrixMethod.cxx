/********
 * Class for doing fake estimations using Matrix Method      *
 * Contact: Eirik Gramstad (egramsta@cern.ch)                *
 * 06 Oct: Matrix inversion for 2L channel (4x4 matrix)      *
 ********/
#include "MatrixMethod.h"
#include <iostream>
#include <math.h>

/* 4-D inversion functions
 * The functions                         
 *      void N4_xx_LL()            
 * returns the estimate of RR,RF,FR or FF events in the loose region
 * Input data:
 * float r1,float fE,float r2,float f2 - real and fake eff. for the two leptons
 * float NTT, float NTl, float NlT, float Nll - number of TT,Tl,lT,ll events in the estimation region
 */
float MatrixMethod::N4_RR_LL(float r1,float f1,float r2,float f2, 
			     float NTT, float NTl, float NlT, float Nll){
  float den, num, num2;
  den = (r1 - f1)*(r2 - f2);
  num  = ((f1 - 1)*(f2 - 1)) * NTT +  (f2*(f1 - 1)) * NTl +  (f1*(f2 - 1)) * NlT +  (f1*f2) * Nll;
  num2 = ((1 - f1)*(1 - f2)) * NTT -  (f2*(1 - f1)) * NTl -  (f1*(1 - f2)) * NlT +  (f1*f2) * Nll;
  if (den==0) return -1; 
  else return num/den; 
}
float MatrixMethod::N4_RF_LL(float r1,float f1,float r2,float f2, 
			     float NTT, float NTl, float NlT, float Nll){
  float den, num, num2;
  den = (r1 - f1)*(r2 - f2);
  num =  -((f1 - 1)*(r2 - 1)) * NTT -  (r2*(f1 - 1)) * NTl -  (f1*(r2 - 1)) * NlT -  (f1*r2) * Nll;
  num2 = -((1 - f1)*(1 - r2)) * NTT +  (r2*(1 - f1)) * NTl +  (f1*(1 - r2)) * NlT -  (f1*r2) * Nll;
  if (den==0) return -1;
  else return num/den;
}

float MatrixMethod::N4_FR_LL(float r1,float f1,float r2,float f2, 
			     float NTT, float NTl, float NlT, float Nll){
  float den, num, num2;
  den = (r1 - f1)*(r2 - f2);
  num  = -((f2 - 1)*(r1 - 1)) * NTT -  (f2*(r1 - 1)) * NTl -  (r1*(f2 - 1)) * NlT -  (f2*r1) * Nll;
  num2 = -((1 - f2)*(1 - r1)) * NTT +  (f2*(1 - r1)) * NTl +  (r1*(1 - f2)) * NlT -  (f2*r1) * Nll;
  if (den==0) return -1;
  else return num/den; 
}

float MatrixMethod::N4_FF_LL(float r1,float f1,float r2,float f2, 
			     float NTT, float NTl, float NlT, float Nll){
  float den, num, num2;
  den = (r1 - f1)*(r2 - f2);
  num  =  ((r1 - 1)*(r2 - 1)) * NTT +  (r2*(r1 - 1)) * NTl +  (r1*(r2 - 1)) * NlT +  (r1*r2) * Nll;
  num2 =  ((1 - r1)*(1 - r2)) * NTT -  (r2*(1 - r1)) * NTl -  (r1*(1 - r2)) * NlT +  (r1*r2) * Nll;
  if (den==0) return -1; 
  else return num/den; 
}

/* 4-D inversion functions
 * The functions                         
 *      void N4_xx_TT()            
 * returns the estimate of RR,RF,FR or FF events in a tight region
 * Input data:
 * float r1,float f1,float r2,float f2 - real and fake eff. for the two leptons
 * float NTT, float NTl, float NlT, float Nll - number of TT,Tl,lT,ll events in the estimation region
 */
float MatrixMethod::N4_RR_TT(float r1,float f1,float r2,float f2, 
			     float NTT, float NTl, float NlT, float Nll){
  //cout<<NTT<<","<<NTl<<","<<NlT<<","<<Nll<<endl;
  return r1 * r2 * N4_RR_LL(r1,f1,r2,f2, NTT, NTl, NlT, Nll);
}

float MatrixMethod::N4_RF_TT(float r1,float f1,float r2,float f2, 
			     float NTT, float NTl, float NlT, float Nll){
  return r1 * f2 * N4_RF_LL(r1,f1,r2,f2, NTT, NTl, NlT, Nll);
}

float MatrixMethod::N4_FR_TT(float r1,float f1,float r2,float f2, 
			     float NTT, float NTl, float NlT, float Nll){
  return f1 * r2 * N4_FR_LL(r1,f1,r2,f2, NTT, NTl, NlT, Nll);
}

float MatrixMethod::N4_FF_TT(float r1,float f1,float r2,float f2, 
			     float NTT, float NTl, float NlT, float Nll){
  return f1 * f2 * N4_FF_LL(r1,f1,r2,f2, NTT, NTl, NlT, Nll);
}

/* 4-D error functions
 * The functions                         
 *      void N4_xx_TTerr()            
 * returns the error on the estimates RR,RF,FR,FF
 * Input data:
 * float r1,float f1,float r2,float f2 - real and fake eff. for the two leptons
 * float NTT, float NTl, float NlT, float Nll - number of TT,Tl,lT,ll events in the estimation region
 * float r1err,float f1err,float r2err,float f2err - errors on the real and fake efficiencies
 * formulas are extracted by doing standard extrapolation of the expressions after matrix inversion
 */
float MatrixMethod::N4_RR_TTerr(float r1,float f1,float r2,float f2,  float nTT,float nTl,float nlT,float nll,  float r1err,float f1err,float r2err,float f2err){
  float err_RR;

  float nTTerr, nTlerr, nlTerr,nllerr;
  nTTerr = sqrt(nTT);
  nTlerr = sqrt(nTl);
  nlTerr = sqrt(nlT);
  nllerr = sqrt(nll);

  float err = r1err*r1err*(r2*((nTT*(f1 - 1)*(f2 - 1))/((f1 - r1)*(f2 - r2)) + (f1*f2*nll)/((f1 - r1)*(f2 - r2)) + (f2*nTl*(f1 - 1))/((f1 - r1)*(f2 - r2)) + (f1*nlT*(f2 - 1))/((f1 - r1)*(f2 - r2))) + r1*r2*((nTT*(f1 - 1)*(f2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (f1*f2*nll)/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (f2*nTl*(f1 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (f1*nlT*(f2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2))))*(r2*((nTT*(f1 - 1)*(f2 - 1))/((f1 - r1)*(f2 - r2)) + (f1*f2*nll)/((f1 - r1)*(f2 - r2)) + (f2*nTl*(f1 - 1))/((f1 - r1)*(f2 - r2)) + (f1*nlT*(f2 - 1))/((f1 - r1)*(f2 - r2))) + r1*r2*((nTT*(f1 - 1)*(f2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (f1*f2*nll)/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (f2*nTl*(f1 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (f1*nlT*(f2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)))) + r2err*r2err*(r1*((nTT*(f1 - 1)*(f2 - 1))/((f1 - r1)*(f2 - r2)) + (f1*f2*nll)/((f1 - r1)*(f2 - r2)) + (f2*nTl*(f1 - 1))/((f1 - r1)*(f2 - r2)) + (f1*nlT*(f2 - 1))/((f1 - r1)*(f2 - r2))) + r1*r2*((nTT*(f1 - 1)*(f2 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (f1*f2*nll)/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (f2*nTl*(f1 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (f1*nlT*(f2 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2))))*(r1*((nTT*(f1 - 1)*(f2 - 1))/((f1 - r1)*(f2 - r2)) + (f1*f2*nll)/((f1 - r1)*(f2 - r2)) + (f2*nTl*(f1 - 1))/((f1 - r1)*(f2 - r2)) + (f1*nlT*(f2 - 1))/((f1 - r1)*(f2 - r2))) + r1*r2*((nTT*(f1 - 1)*(f2 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (f1*f2*nll)/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (f2*nTl*(f1 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (f1*nlT*(f2 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)))) + f2err*f2err*r1*r1*r2*r2*((f1*nlT)/((f1 - r1)*(f2 - r2)) + (f1*nll)/((f1 - r1)*(f2 - r2)) + (nTT*(f1 - 1))/((f1 - r1)*(f2 - r2)) + (nTl*(f1 - 1))/((f1 - r1)*(f2 - r2)) - (nTT*(f1 - 1)*(f2 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)) - (f1*f2*nll)/((f1 - r1)*(f2 - r2)*(f2 - r2)) - (f2*nTl*(f1 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)) - (f1*nlT*(f2 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)))*((f1*nlT)/((f1 - r1)*(f2 - r2)) + (f1*nll)/((f1 - r1)*(f2 - r2)) + (nTT*(f1 - 1))/((f1 - r1)*(f2 - r2)) + (nTl*(f1 - 1))/((f1 - r1)*(f2 - r2)) - (nTT*(f1 - 1)*(f2 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)) - (f1*f2*nll)/((f1 - r1)*(f2 - r2)*(f2 - r2)) - (f2*nTl*(f1 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)) - (f1*nlT*(f2 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2))) + f1err*f1err*r1*r1*r2*r2*((f2*nTl)/((f1 - r1)*(f2 - r2)) + (f2*nll)/((f1 - r1)*(f2 - r2)) + (nTT*(f2 - 1))/((f1 - r1)*(f2 - r2)) + (nlT*(f2 - 1))/((f1 - r1)*(f2 - r2)) - (nTT*(f1 - 1)*(f2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)) - (f1*f2*nll)/((f1 - r1)*(f1 - r1)*(f2 - r2)) - (f2*nTl*(f1 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)) - (f1*nlT*(f2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)))*((f2*nTl)/((f1 - r1)*(f2 - r2)) + (f2*nll)/((f1 - r1)*(f2 - r2)) + (nTT*(f2 - 1))/((f1 - r1)*(f2 - r2)) + (nlT*(f2 - 1))/((f1 - r1)*(f2 - r2)) - (nTT*(f1 - 1)*(f2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)) - (f1*f2*nll)/((f1 - r1)*(f1 - r1)*(f2 - r2)) - (f2*nTl*(f1 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)) - (f1*nlT*(f2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2))) + (nTTerr*nTTerr*r1*r1*r2*r2*(f1 - 1)*(f1 - 1)*(f2 - 1)*(f2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)*(f2 - r2)) + (f1*f1*f2*f2*nllerr*nllerr*r1*r1*r2*r2)/((f1 - r1)*(f1 - r1)*(f2 - r2)*(f2 - r2)) + (f2*f2*nTlerr*nTlerr*r1*r1*r2*r2*(f1 - 1)*(f1 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)*(f2 - r2)) + (f1*f1*nlTerr*nlTerr*r1*r1*r2*r2*(f2 - 1)*(f2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)*(f2 - r2));
  if(err<0){ // superfluous test
    printf("N4_FF_TTerr: imaginary !!!  f^2 = %f   (returning -sqrt(-f^2)",err);
    err_RR = -sqrt(-err);
  }else{ err_RR = sqrt(err);}

  return err_RR;
}


float MatrixMethod::N4_RF_TTerr(float r1,float f1,float r2,float f2,  float nTT,float nTl,float nlT,float nll,  float r1err,float f1err,float r2err,float f2err){
  float err_RF;

  float nTTerr, nTlerr, nlTerr,nllerr;
  nTTerr = sqrt(nTT);
  nTlerr = sqrt(nTl);
  nlTerr = sqrt(nlT);
  nllerr = sqrt(nll);

  float err = r1err*r1err*(f2*((nTT*(f1 - 1)*(r2 - 1))/((f1 - r1)*(f2 - r2)) + (f1*nll*r2)/((f1 - r1)*(f2 - r2)) + (nTl*r2*(f1 - 1))/((f1 - r1)*(f2 - r2)) + (f1*nlT*(r2 - 1))/((f1 - r1)*(f2 - r2))) + f2*r1*((nTT*(f1 - 1)*(r2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (f1*nll*r2)/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (nTl*r2*(f1 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (f1*nlT*(r2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2))))*(f2*((nTT*(f1 - 1)*(r2 - 1))/((f1 - r1)*(f2 - r2)) + (f1*nll*r2)/((f1 - r1)*(f2 - r2)) + (nTl*r2*(f1 - 1))/((f1 - r1)*(f2 - r2)) + (f1*nlT*(r2 - 1))/((f1 - r1)*(f2 - r2))) + f2*r1*((nTT*(f1 - 1)*(r2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (f1*nll*r2)/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (nTl*r2*(f1 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (f1*nlT*(r2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)))) + f2err*f2err*(r1*((nTT*(f1 - 1)*(r2 - 1))/((f1 - r1)*(f2 - r2)) + (f1*nll*r2)/((f1 - r1)*(f2 - r2)) + (nTl*r2*(f1 - 1))/((f1 - r1)*(f2 - r2)) + (f1*nlT*(r2 - 1))/((f1 - r1)*(f2 - r2))) - f2*r1*((nTT*(f1 - 1)*(r2 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (f1*nll*r2)/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (nTl*r2*(f1 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (f1*nlT*(r2 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2))))*(r1*((nTT*(f1 - 1)*(r2 - 1))/((f1 - r1)*(f2 - r2)) + (f1*nll*r2)/((f1 - r1)*(f2 - r2)) + (nTl*r2*(f1 - 1))/((f1 - r1)*(f2 - r2)) + (f1*nlT*(r2 - 1))/((f1 - r1)*(f2 - r2))) - f2*r1*((nTT*(f1 - 1)*(r2 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (f1*nll*r2)/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (nTl*r2*(f1 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (f1*nlT*(r2 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)))) + f2*f2*r1*r1*r2err*r2err*((f1*nlT)/((f1 - r1)*(f2 - r2)) + (f1*nll)/((f1 - r1)*(f2 - r2)) + (nTT*(f1 - 1))/((f1 - r1)*(f2 - r2)) + (nTl*(f1 - 1))/((f1 - r1)*(f2 - r2)) + (nTT*(f1 - 1)*(r2 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (f1*nll*r2)/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (nTl*r2*(f1 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (f1*nlT*(r2 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)))*((f1*nlT)/((f1 - r1)*(f2 - r2)) + (f1*nll)/((f1 - r1)*(f2 - r2)) + (nTT*(f1 - 1))/((f1 - r1)*(f2 - r2)) + (nTl*(f1 - 1))/((f1 - r1)*(f2 - r2)) + (nTT*(f1 - 1)*(r2 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (f1*nll*r2)/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (nTl*r2*(f1 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (f1*nlT*(r2 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2))) + f1err*f1err*f2*f2*r1*r1*((nTl*r2)/((f1 - r1)*(f2 - r2)) + (nll*r2)/((f1 - r1)*(f2 - r2)) + (nTT*(r2 - 1))/((f1 - r1)*(f2 - r2)) + (nlT*(r2 - 1))/((f1 - r1)*(f2 - r2)) - (nTT*(f1 - 1)*(r2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)) - (f1*nll*r2)/((f1 - r1)*(f1 - r1)*(f2 - r2)) - (nTl*r2*(f1 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)) - (f1*nlT*(r2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)))*((nTl*r2)/((f1 - r1)*(f2 - r2)) + (nll*r2)/((f1 - r1)*(f2 - r2)) + (nTT*(r2 - 1))/((f1 - r1)*(f2 - r2)) + (nlT*(r2 - 1))/((f1 - r1)*(f2 - r2)) - (nTT*(f1 - 1)*(r2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)) - (f1*nll*r2)/((f1 - r1)*(f1 - r1)*(f2 - r2)) - (nTl*r2*(f1 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)) - (f1*nlT*(r2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2))) + (f2*f2*nTTerr*nTTerr*r1*r1*(f1 - 1)*(f1 - 1)*(r2 - 1)*(r2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)*(f2 - r2)) + (f1*f1*f2*f2*nllerr*nllerr*r1*r1*r2*r2)/((f1 - r1)*(f1 - r1)*(f2 - r2)*(f2 - r2)) + (f2*f2*nTlerr*nTlerr*r1*r1*r2*r2*(f1 - 1)*(f1 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)*(f2 - r2)) + (f1*f1*f2*f2*nlTerr*nlTerr*r1*r1*(r2 - 1)*(r2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)*(f2 - r2));
  if(err<0){ // superfluous test
    printf("N4_FF_TTerr: imaginary !!!  f^2 = %f   (returning -sqrt(-f^2)",err);
    err_RF = -sqrt(-err);
  }else{ err_RF = sqrt(err);}

  return err_RF;
}

float MatrixMethod::N4_FR_TTerr(float r1,float f1,float r2,float f2,  float nTT,float nTl,float nlT,float nll,  float r1err,float f1err,float r2err,float f2err){
  float err_FR;

  float nTTerr, nTlerr, nlTerr,nllerr;
  nTTerr = sqrt(nTT);
  nTlerr = sqrt(nTl);
  nlTerr = sqrt(nlT);
  nllerr = sqrt(nll);

  float err = r2err*r2err*(f1*((nTT*(f2 - 1)*(r1 - 1))/((f1 - r1)*(f2 - r2)) + (f2*nll*r1)/((f1 - r1)*(f2 - r2)) + (f2*nTl*(r1 - 1))/((f1 - r1)*(f2 - r2)) + (nlT*r1*(f2 - 1))/((f1 - r1)*(f2 - r2))) + f1*r2*((nTT*(f2 - 1)*(r1 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (f2*nll*r1)/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (f2*nTl*(r1 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (nlT*r1*(f2 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2))))*(f1*((nTT*(f2 - 1)*(r1 - 1))/((f1 - r1)*(f2 - r2)) + (f2*nll*r1)/((f1 - r1)*(f2 - r2)) + (f2*nTl*(r1 - 1))/((f1 - r1)*(f2 - r2)) + (nlT*r1*(f2 - 1))/((f1 - r1)*(f2 - r2))) + f1*r2*((nTT*(f2 - 1)*(r1 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (f2*nll*r1)/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (f2*nTl*(r1 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (nlT*r1*(f2 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)))) + f1err*f1err*(r2*((nTT*(f2 - 1)*(r1 - 1))/((f1 - r1)*(f2 - r2)) + (f2*nll*r1)/((f1 - r1)*(f2 - r2)) + (f2*nTl*(r1 - 1))/((f1 - r1)*(f2 - r2)) + (nlT*r1*(f2 - 1))/((f1 - r1)*(f2 - r2))) - f1*r2*((nTT*(f2 - 1)*(r1 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (f2*nll*r1)/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (f2*nTl*(r1 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (nlT*r1*(f2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2))))*(r2*((nTT*(f2 - 1)*(r1 - 1))/((f1 - r1)*(f2 - r2)) + (f2*nll*r1)/((f1 - r1)*(f2 - r2)) + (f2*nTl*(r1 - 1))/((f1 - r1)*(f2 - r2)) + (nlT*r1*(f2 - 1))/((f1 - r1)*(f2 - r2))) - f1*r2*((nTT*(f2 - 1)*(r1 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (f2*nll*r1)/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (f2*nTl*(r1 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (nlT*r1*(f2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)))) + f1*f1*r1err*r1err*r2*r2*((f2*nTl)/((f1 - r1)*(f2 - r2)) + (f2*nll)/((f1 - r1)*(f2 - r2)) + (nTT*(f2 - 1))/((f1 - r1)*(f2 - r2)) + (nlT*(f2 - 1))/((f1 - r1)*(f2 - r2)) + (nTT*(f2 - 1)*(r1 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (f2*nll*r1)/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (f2*nTl*(r1 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (nlT*r1*(f2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)))*((f2*nTl)/((f1 - r1)*(f2 - r2)) + (f2*nll)/((f1 - r1)*(f2 - r2)) + (nTT*(f2 - 1))/((f1 - r1)*(f2 - r2)) + (nlT*(f2 - 1))/((f1 - r1)*(f2 - r2)) + (nTT*(f2 - 1)*(r1 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (f2*nll*r1)/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (f2*nTl*(r1 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (nlT*r1*(f2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2))) + f1*f1*f2err*f2err*r2*r2*((nlT*r1)/((f1 - r1)*(f2 - r2)) + (nll*r1)/((f1 - r1)*(f2 - r2)) + (nTT*(r1 - 1))/((f1 - r1)*(f2 - r2)) + (nTl*(r1 - 1))/((f1 - r1)*(f2 - r2)) - (nTT*(f2 - 1)*(r1 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)) - (f2*nll*r1)/((f1 - r1)*(f2 - r2)*(f2 - r2)) - (f2*nTl*(r1 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)) - (nlT*r1*(f2 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)))*((nlT*r1)/((f1 - r1)*(f2 - r2)) + (nll*r1)/((f1 - r1)*(f2 - r2)) + (nTT*(r1 - 1))/((f1 - r1)*(f2 - r2)) + (nTl*(r1 - 1))/((f1 - r1)*(f2 - r2)) - (nTT*(f2 - 1)*(r1 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)) - (f2*nll*r1)/((f1 - r1)*(f2 - r2)*(f2 - r2)) - (f2*nTl*(r1 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)) - (nlT*r1*(f2 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2))) + (f1*f1*nTTerr*nTTerr*r2*r2*(f2 - 1)*(f2 - 1)*(r1 - 1)*(r1 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)*(f2 - r2)) + (f1*f1*f2*f2*nllerr*nllerr*r1*r1*r2*r2)/((f1 - r1)*(f1 - r1)*(f2 - r2)*(f2 - r2)) + (f1*f1*f2*f2*nTlerr*nTlerr*r2*r2*(r1 - 1)*(r1 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)*(f2 - r2)) + (f1*f1*nlTerr*nlTerr*r1*r1*r2*r2*(f2 - 1)*(f2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)*(f2 - r2));
  if(err<0){ // superfluous test
    printf("N4_FF_TTerr: imaginary !!!  f^2 = %f   (returning -sqrt(-f^2)",err);
    err_FR = -sqrt(-err);
  }else{ err_FR = sqrt(err);}

  return err_FR;
}


float MatrixMethod::N4_FF_TTerr(float r1,float f1,float r2,float f2,  float nTT,float nTl,float nlT,float nll,  float r1err,float f1err,float r2err,float f2err){
  float err_FF;

  float nTTerr, nTlerr, nlTerr,nllerr;
  nTTerr = sqrt(nTT);
  nTlerr = sqrt(nTl);
  nlTerr = sqrt(nlT);
  nllerr = sqrt(nll);
  
  float err = f1err*f1err*(f2*((nTT*(r1 - 1)*(r2 - 1))/((f1 - r1)*(f2 - r2)) + (nll*r1*r2)/((f1 - r1)*(f2 - r2)) + (nTl*r2*(r1 - 1))/((f1 - r1)*(f2 - r2)) + (nlT*r1*(r2 - 1))/((f1 - r1)*(f2 - r2))) - f1*f2*((nTT*(r1 - 1)*(r2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (nll*r1*r2)/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (nTl*r2*(r1 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (nlT*r1*(r2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2))))*(f2*((nTT*(r1 - 1)*(r2 - 1))/((f1 - r1)*(f2 - r2)) + (nll*r1*r2)/((f1 - r1)*(f2 - r2)) + (nTl*r2*(r1 - 1))/((f1 - r1)*(f2 - r2)) + (nlT*r1*(r2 - 1))/((f1 - r1)*(f2 - r2))) - f1*f2*((nTT*(r1 - 1)*(r2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (nll*r1*r2)/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (nTl*r2*(r1 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (nlT*r1*(r2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)))) + f2err*f2err*(f1*((nTT*(r1 - 1)*(r2 - 1))/((f1 - r1)*(f2 - r2)) + (nll*r1*r2)/((f1 - r1)*(f2 - r2)) + (nTl*r2*(r1 - 1))/((f1 - r1)*(f2 - r2)) + (nlT*r1*(r2 - 1))/((f1 - r1)*(f2 - r2))) - f1*f2*((nTT*(r1 - 1)*(r2 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (nll*r1*r2)/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (nTl*r2*(r1 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (nlT*r1*(r2 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2))))*(f1*((nTT*(r1 - 1)*(r2 - 1))/((f1 - r1)*(f2 - r2)) + (nll*r1*r2)/((f1 - r1)*(f2 - r2)) + (nTl*r2*(r1 - 1))/((f1 - r1)*(f2 - r2)) + (nlT*r1*(r2 - 1))/((f1 - r1)*(f2 - r2))) - f1*f2*((nTT*(r1 - 1)*(r2 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (nll*r1*r2)/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (nTl*r2*(r1 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (nlT*r1*(r2 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)))) + f1*f1*f2*f2*r2err*r2err*((nlT*r1)/((f1 - r1)*(f2 - r2)) + (nll*r1)/((f1 - r1)*(f2 - r2)) + (nTT*(r1 - 1))/((f1 - r1)*(f2 - r2)) + (nTl*(r1 - 1))/((f1 - r1)*(f2 - r2)) + (nTT*(r1 - 1)*(r2 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (nll*r1*r2)/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (nTl*r2*(r1 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (nlT*r1*(r2 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)))*((nlT*r1)/((f1 - r1)*(f2 - r2)) + (nll*r1)/((f1 - r1)*(f2 - r2)) + (nTT*(r1 - 1))/((f1 - r1)*(f2 - r2)) + (nTl*(r1 - 1))/((f1 - r1)*(f2 - r2)) + (nTT*(r1 - 1)*(r2 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (nll*r1*r2)/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (nTl*r2*(r1 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2)) + (nlT*r1*(r2 - 1))/((f1 - r1)*(f2 - r2)*(f2 - r2))) + f1*f1*f2*f2*r1err*r1err*((nTl*r2)/((f1 - r1)*(f2 - r2)) + (nll*r2)/((f1 - r1)*(f2 - r2)) + (nTT*(r2 - 1))/((f1 - r1)*(f2 - r2)) + (nlT*(r2 - 1))/((f1 - r1)*(f2 - r2)) + (nTT*(r1 - 1)*(r2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (nll*r1*r2)/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (nTl*r2*(r1 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (nlT*r1*(r2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)))*((nTl*r2)/((f1 - r1)*(f2 - r2)) + (nll*r2)/((f1 - r1)*(f2 - r2)) + (nTT*(r2 - 1))/((f1 - r1)*(f2 - r2)) + (nlT*(r2 - 1))/((f1 - r1)*(f2 - r2)) + (nTT*(r1 - 1)*(r2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (nll*r1*r2)/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (nTl*r2*(r1 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)) + (nlT*r1*(r2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2))) + (f1*f1*f2*f2*nTTerr*nTTerr*(r1 - 1)*(r1 - 1)*(r2 - 1)*(r2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)*(f2 - r2)) + (f1*f1*f2*f2*nllerr*nllerr*r1*r1*r2*r2)/((f1 - r1)*(f1 - r1)*(f2 - r2)*(f2 - r2)) + (f1*f1*f2*f2*nTlerr*nTlerr*r2*r2*(r1 - 1)*(r1 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)*(f2 - r2)) + (f1*f1*f2*f2*nlTerr*nlTerr*r1*r1*(r2 - 1)*(r2 - 1))/((f1 - r1)*(f1 - r1)*(f2 - r2)*(f2 - r2));
  if(err<0){ // superfluous test
    printf("N4_FF_TTerr: imaginary !!!  f^2 = %f   (returning -sqrt(-f^2)",err);
    err_FF = -sqrt(-err);
  }else{ err_FF = sqrt(err);}

  return err_FF;
}


/* 3-D inversion functions
 * The functions                         
 *      void N3_xx_LL()            
 * returns the estimate of RR,RF or FF events in a loose region. In this case 
 * we do not differ between Tl,lT and RF,FR events. May be used in EE and MM channels.
 * Input data:
 * float r,float f - real and fake eff. for the two leptons (assumed to be the same)
 * float NTT, float NTl, float Nll - number of TT,Tl,ll events in the estimation region
 */
float MatrixMethod::N3_RR_LL(float r, float f, float NTT,float NTl,float Nll){
  float fac = (r-f)*(r-f);
  float ret = NTT*(1-f)*(1-f) - NTl*(f*(1-f)) + Nll*f*f;
  if(fac==0) return -1;
  
  return ret/fac;
  
}
float MatrixMethod::N3_RF_LL(float r, float f, float NTT,float NTl,float Nll){
  float fac = (r-f)*(r-f);
  float ret = NTT*(-2*(1-r)*(1-f)) + NTl*(r+f-2*r*f) + Nll*(-2*r*f);
  if(fac==0) return -1;
  
  return ret/fac;    
}
float MatrixMethod::N3_FF_LL(float r, float f, float NTT,float NTl,float Nll){
  
  float fac = (r-f)*(r-f);
  float ret = NTT*(1-r)*(1-r) - NTl*(r*(1-r)) + Nll*r*r;
  if(fac==0) return -1;

  return ret/fac;
}

/* 3-D inversion functions
 * The functions                         
 *      void N3_xx_TT()            
 * returns the estimate of RR,RF or FF events in a tight region. In this case 
 * we do not differ between Tl,lT and RF,FR events. May be used in EE and MM channels.
 * Input data:
 * float r,float f - real and fake eff. for the two leptons (assumed to be the same)
 * float NTT, float NTl, float Nll - number of TT,Tl,ll events in the estimation region
 */
float MatrixMethod::N3_RR_TT(float r, float f, float NTT,float NTl,float Nll){
  return r*r*N3_RR_LL(r,f,NTT,NTl,Nll);
}
float MatrixMethod::N3_RF_TT(float r, float f, float NTT,float NTl,float Nll){
  return r*f*N3_RF_LL(r,f,NTT,NTl,Nll);
}
float MatrixMethod::N3_FF_TT(float r, float f, float NTT,float NTl,float Nll){
  return f*f*N3_FF_LL(r,f,NTT,NTl,Nll);
}

/* 3-D error functions
 * The functions                         
 *      void N3_xx_TTerr()            
 * returns the error on the estimates RR,RF,FF
 * Input data:
 * float r,float f - real and fake eff. for the leptons
 * float NTT, float NTl, float Nll - number of TT,Tl (i.e. Tl+lT),ll events in the estimation region
 * float rerr,float ferr - errors on the real and fake efficiencies
 * formulas are extracted by doing standard extrapolation of the expressions after matrix inversion
 */
float MatrixMethod::N3_RR_TTerr(float r, float f, float NTT,float NTl,float Nll, float rerr, float ferr){

  float A = (float)NTT;
  float B = (float)NTl;
  float D = (float)Nll;

  float sigma_A = sqrt(A);
  float sigma_B = sqrt(B);
  float sigma_D = sqrt(D);

  float drrdA = ((1-f)*(1-f)*r*r)/((r-f)*(r-f));
  float drrdf = (((-2*A*(1-f)+2*D*f)*r*r)/((r-f)*(r-f))+((2*(A*(1-f)*(1-f)-B*f*(1-f)+D*f*f))*r*r)/((r-f)*(r-f)*(r-f)));
  float drrdr = ((-(2*(A*(1-f)*(1-f)-B*f*(1-f)+D*f*f))*r*r)/((r-f)*(r-f)*(r-f))+((2*(A*(1-f)*(1-f)-B*f*(1-f)+D*f*f))*r)/((r-f)*(r-f)));
  float drrdB = -(f*(1-f)*r*r)/((r-f)*(r-f));
  float drrdD = (f*f*r*r)/((r-f)*(r-f));

  float error_RR = sqrt(pow((drrdA*sigma_A),2) + pow((drrdB*sigma_B),2) + pow((drrdD*sigma_D),2) + pow((drrdr*rerr),2) + pow((drrdf*ferr),2));
  
  return error_RR;

}

float MatrixMethod::N3_RF_TTerr(float r, float f, float NTT,float NTl,float Nll, float rerr, float ferr){

  float A = (float)NTT;
  float B = (float)NTl;
  float D = (float)Nll;

  float sigma_A = sqrt(A);
  float sigma_B = sqrt(B);
  float sigma_D = sqrt(D);

  float drfdA = (-(2-2*r)*(1-f)*r*f)/((r-f)*(r-f));
  float drfdf = (((A*(2-2*r)+B*(1-2*r)-2*D*r)*r*f)/((r-f)*(r-f))+((-A*(2-2*r)*(1-f)+B*(r+f-2*r*f)-2*D*r*f)*r)/((r-f)*(r-f))+((2*(-A*(2-2*r)*(1-f)+B*(r+f-2*r*f)-2*D*r*f))*r*f)/((r-f)*(r-f)*(r-f)));
  float drfdr = (((2*A*(1-f)+B*(1-2*f)-2*D*f)*r*f)/((r-f)*(r-f))+((-A*(2-2*r)*(1-f)+B*(r+f-2*r*f)-2*D*r*f)*f)/((r-f)*(r-f))-((2*(-A*(2-2*r)*(1-f)+B*(r+f-2*r*f)-2*D*r*f))*r*f)/((r-f)*(r-f)*(r-f)));
  float drfdB = ((r+f-2*r*f)*r*f)/((r-f)*(r-f));
  float drfdD = -(2*f*f*r*r)/((r-f)*(r-f));

  float error_RF = sqrt(pow((drfdA*sigma_A),2) + pow((drfdB*sigma_B),2) + pow((drfdD*sigma_D),2) + pow((drfdr*rerr),2) + pow((drfdf*ferr),2));
  
  return error_RF;

}

float MatrixMethod::N3_FF_TTerr(float r, float f, float NTT,float NTl,float Nll, float rerr, float ferr){

  float A = (float)NTT;
  float B = (float)NTl;
  float D = (float)Nll;

  float sigma_A = sqrt(A);
  float sigma_B = sqrt(B);
  float sigma_D = sqrt(D);

  float dffdA = ((1-r)*(1-r)*f*f)/((f-r)*(f-r));
  float dffdf = ((-(2*(A*(1-r)*(1-r)-B*r*(1-r)+D*r*r))*f*f)/((f-r)*(f-r)*(f-r))+((2*(A*(1-r)*(1-r)-B*r*(1-r)+D*r*r))*f)/((f-r)*(f-r)));
  float dffdr = (((-2*A*(1-r)+2*D*r)*f*f)/((f-r)*(f-r))+((2*(A*(1-r)*(1-r)-B*r*(1-r)+D*r*r))*f*f)/((f-r)*(f-r)*(f-r)));
  float dffdB = -(r*(1-r)*f*f)/((f-r)*(f-r));
  float dffdD = (r*r*f*f)/((f-r)*(f-r));

  float error_FF = sqrt(pow((dffdA*sigma_A),2) + pow((dffdB*sigma_B),2) + pow((dffdD*sigma_D),2) + pow((dffdr*rerr),2) + pow((dffdf*ferr),2));

  
  return error_FF;

}

 
/* 8-D inversion functions
 * The functions                         
 *      void N8_xx_LLL()            
 * returns the estimate of RRR,RRF,RFR,FRR,RFF,FRF,FFR,FFF events in a loose region.
 * Input data:
 * float r1,float f1,float r2,float f2,float r3,float f3 - real and fake eff. for the three leptons
 * int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll - number of events in the various categories for the estimation region
 */
float MatrixMethod::N8_RRR_LLL(float r1,float f1,float r2,float f2,float r3,float f3, int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll){
  return ((lTT*f1 - TTT + TlT*f2 + TTT*f1 + TTT*f2 + TTl*f3 + TTT*f3 - llT*f1*f2 - lTT*f1*f2 - lTl*f1*f3 - lTT*f1*f3 - TlT*f1*f2 - Tll*f2*f3 - TlT*f2*f3 - TTT*f1*f2 - TTl*f1*f3 - TTT*f1*f3 - TTl*f2*f3 - TTT*f2*f3 + lll*f1*f2*f3 + llT*f1*f2*f3 + lTl*f1*f2*f3 + lTT*f1*f2*f3 + Tll*f1*f2*f3 + TlT*f1*f2*f3 + TTl*f1*f2*f3 + TTT*f1*f2*f3))/((f1 - r1)*(f2 - r2)*(f3 - r3));
}
float MatrixMethod::N8_RRF_LLL(float r1,float f1,float r2,float f2,float r3,float f3, int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll){
  return -((lTT*f1 - TTT + TlT*f2 + TTT*f1 + TTT*f2 + TTl*r3 + TTT*r3 - llT*f1*f2 - lTT*f1*f2 - TlT*f1*f2 - TTT*f1*f2 - lTl*f1*r3 - lTT*f1*r3 - Tll*f2*r3 - TlT*f2*r3 - TTl*f1*r3 - TTT*f1*r3 - TTl*f2*r3 - TTT*f2*r3 + lll*f1*f2*r3 + llT*f1*f2*r3 + lTl*f1*f2*r3 + lTT*f1*f2*r3 + Tll*f1*f2*r3 + TlT*f1*f2*r3 + TTl*f1*f2*r3 + TTT*f1*f2*r3))/((f1 - r1)*(f2 - r2)*(f3 - r3));
    
}
float MatrixMethod::N8_RFR_LLL(float r1,float f1,float r2,float f2,float r3,float f3, int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll){
  return -((lTT*f1 - TTT + TTT*f1 + TTl*f3 + TTT*f3 + TlT*r2 + TTT*r2 - lTl*f1*f3 - lTT*f1*f3 - TTl*f1*f3 - TTT*f1*f3 - llT*f1*r2 - lTT*f1*r2 - TlT*f1*r2 - Tll*f3*r2 - TlT*f3*r2 - TTT*f1*r2 - TTl*f3*r2 - TTT*f3*r2 + lll*f1*f3*r2 + llT*f1*f3*r2 + lTl*f1*f3*r2 + lTT*f1*f3*r2 + Tll*f1*f3*r2 + TlT*f1*f3*r2 + TTl*f1*f3*r2 + TTT*f1*f3*r2))/((f1 - r1)*(f2 - r2)*(f3 - r3));
}
float MatrixMethod::N8_FRR_LLL(float r1,float f1,float r2,float f2,float r3,float f3, int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll){
  return -((TlT*f2 - TTT + TTT*f2 + TTl*f3 + TTT*f3 + lTT*r1 + TTT*r1 - Tll*f2*f3 - TlT*f2*f3 - TTl*f2*f3 - TTT*f2*f3 - llT*f2*r1 - lTT*f2*r1 - lTl*f3*r1 - lTT*f3*r1 - TlT*f2*r1 - TTT*f2*r1 - TTl*f3*r1 - TTT*f3*r1 + lll*f2*f3*r1 + llT*f2*f3*r1 + lTl*f2*f3*r1 + lTT*f2*f3*r1 + Tll*f2*f3*r1 + TlT*f2*f3*r1 + TTl*f2*f3*r1 + TTT*f2*f3*r1))/((f1 - r1)*(f2 - r2)*(f3 - r3));
  
}
float MatrixMethod::N8_RFF_LLL(float r1,float f1,float r2,float f2,float r3,float f3, int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll){
  return  ((lTT*f1 - TTT + TTT*f1 + TlT*r2 + TTT*r2 + TTl*r3 + TTT*r3 - llT*f1*r2 - lTT*f1*r2 - lTl*f1*r3 - lTT*f1*r3 - TlT*f1*r2 - TTT*f1*r2 - TTl*f1*r3 - TTT*f1*r3 - Tll*r2*r3 - TlT*r2*r3 - TTl*r2*r3 - TTT*r2*r3 + lll*f1*r2*r3 + llT*f1*r2*r3 + lTl*f1*r2*r3 + lTT*f1*r2*r3 + Tll*f1*r2*r3 + TlT*f1*r2*r3 + TTl*f1*r2*r3 + TTT*f1*r2*r3))/((f1 - r1)*(f2 - r2)*(f3 - r3));
    
}
float MatrixMethod::N8_FRF_LLL(float r1,float f1,float r2,float f2,float r3,float f3, int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll){
  return  ((TlT*f2 - TTT + TTT*f2 + lTT*r1 + TTT*r1 + TTl*r3 + TTT*r3 - llT*f2*r1 - lTT*f2*r1 - TlT*f2*r1 - Tll*f2*r3 - TlT*f2*r3 - TTT*f2*r1 - TTl*f2*r3 - TTT*f2*r3 - lTl*r1*r3 - lTT*r1*r3 - TTl*r1*r3 - TTT*r1*r3 + lll*f2*r1*r3 + llT*f2*r1*r3 + lTl*f2*r1*r3 + lTT*f2*r1*r3 + Tll*f2*r1*r3 + TlT*f2*r1*r3 + TTl*f2*r1*r3 + TTT*f2*r1*r3))/((f1 - r1)*(f2 - r2)*(f3 - r3));
    
}
float MatrixMethod::N8_FFR_LLL(float r1,float f1,float r2,float f2,float r3,float f3, int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll){
  return  ((TTl*f3 - TTT + TTT*f3 + lTT*r1 + TlT*r2 + TTT*r1 + TTT*r2 - lTl*f3*r1 - lTT*f3*r1 - Tll*f3*r2 - TlT*f3*r2 - TTl*f3*r1 - TTT*f3*r1 - TTl*f3*r2 - TTT*f3*r2 - llT*r1*r2 - lTT*r1*r2 - TlT*r1*r2 - TTT*r1*r2 + lll*f3*r1*r2 + llT*f3*r1*r2 + lTl*f3*r1*r2 + lTT*f3*r1*r2 + Tll*f3*r1*r2 + TlT*f3*r1*r2 + TTl*f3*r1*r2 + TTT*f3*r1*r2))/((f1 - r1)*(f2 - r2)*(f3 - r3));
      
}
float MatrixMethod::N8_FFF_LLL(float r1,float f1,float r2,float f2,float r3,float f3, int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll){
  return -((lTT*r1 - TTT + TlT*r2 + TTT*r1 + TTT*r2 + TTl*r3 + TTT*r3 - llT*r1*r2 - lTT*r1*r2 - lTl*r1*r3 - lTT*r1*r3 - TlT*r1*r2 - Tll*r2*r3 - TlT*r2*r3 - TTT*r1*r2 - TTl*r1*r3 - TTT*r1*r3 - TTl*r2*r3 - TTT*r2*r3 + lll*r1*r2*r3 + llT*r1*r2*r3 + lTl*r1*r2*r3 + lTT*r1*r2*r3 + Tll*r1*r2*r3 + TlT*r1*r2*r3 + TTl*r1*r2*r3 + TTT*r1*r2*r3))/((f1 - r1)*(f2 - r2)*(f3 - r3));
    
}

/* 8-D inversion functions
 * The functions                         
 *      void N8_xx_TTT()            
 * returns the estimate of RRR,RRF,RFR,FRR,RFF,FRF,FFR,FFF events in a loose region.
 * Input data:
 * float r1,float f1,float r2,float f2,float r3,float f3 - real and fake eff. for the three leptons
 * int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll - number of events in the various categories for the estimation region
 */
float MatrixMethod::N8_RRR_TTT(float r1,float f1,float r2,float f2,float r3,float f3, int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll){
  return r1*r2*r3*N8_RRR_LLL(r1,f1,r2,f2,r3,f3,TTT,llT,lTl,Tll,lTT,TlT,TTl,lll);
}
float MatrixMethod::N8_RRF_TTT(float r1,float f1,float r2,float f2,float r3,float f3, int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll){
  return r1*r2*f3*N8_RRF_LLL(r1,f1,r2,f2,r3,f3,TTT,llT,lTl,Tll,lTT,TlT,TTl,lll);
}
float MatrixMethod::N8_RFR_TTT(float r1,float f1,float r2,float f2,float r3,float f3, int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll){
  return r1*f2*r3*N8_RFR_LLL(r1,f1,r2,f2,r3,f3,TTT,llT,lTl,Tll,lTT,TlT,TTl,lll);
}
float MatrixMethod::N8_FRR_TTT(float r1,float f1,float r2,float f2,float r3,float f3, int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll){
  return f1*r2*r3*N8_FRR_LLL(r1,f1,r2,f2,r3,f3,TTT,llT,lTl,Tll,lTT,TlT,TTl,lll);
}
float MatrixMethod::N8_RFF_TTT(float r1,float f1,float r2,float f2,float r3,float f3, int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll){
  return r1*f2*f3*N8_RFF_LLL(r1,f1,r2,f2,r3,f3,TTT,llT,lTl,Tll,lTT,TlT,TTl,lll);
}
float MatrixMethod::N8_FRF_TTT(float r1,float f1,float r2,float f2,float r3,float f3, int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll){
  return f1*r2*f3*N8_FRF_LLL(r1,f1,r2,f2,r3,f3,TTT,llT,lTl,Tll,lTT,TlT,TTl,lll);
}
float MatrixMethod::N8_FFR_TTT(float r1,float f1,float r2,float f2,float r3,float f3, int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll){
  return f1*f2*r3*N8_FFR_LLL(r1,f1,r2,f2,r3,f3,TTT,llT,lTl,Tll,lTT,TlT,TTl,lll);
}
float MatrixMethod::N8_FFF_TTT(float r1,float f1,float r2,float f2,float r3,float f3, int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll){
  return f1*f2*f3*N8_FFF_LLL(r1,f1,r2,f2,r3,f3,TTT,llT,lTl,Tll,lTT,TlT,TTl,lll);
}

/* 8-D error functions
 * The functions                         
 *      void N8_xx_TTTerr()            
 * returns the error on the estimates RRR,RRF etc.
 * Input data:
 * float r1,float f1,float r2,float f2,float r3,float f3 - real and fake eff. for the leptons
 * int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll - number of TT,Tl (i.e. Tl+lT),ll events in the estimation region
 * float r1err,float f1err,float r2err,float f2err,float r3err,float f3err - errors on the real and fake efficiencies
 * formulas are extracted by doing standard extrapolation of the expressions after matrix inversion
 */
float MatrixMethod::N8_RRR_TTTerr(float r1,float f1,float r2,float f2,float r3,float f3, int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll, float r1err,float f1err,float r2err,float f2err,float r3err,float f3err){ 

float TTTerr = sqrt(TTT);
float llTerr = sqrt(llT);
float lTlerr = sqrt(lTl);
float Tllerr = sqrt(Tll);
float lTTerr = sqrt(lTT);
float TlTerr = sqrt(TlT);
float TTlerr = sqrt(TTl);
float lllerr = sqrt(lll);
  
    float res = pow(f1err,2)*pow(((r1*r2*r3*(lTT + TTT - llT*f2 - lTT*f2 - lTl*f3 - lTT*f3 - TlT*f2 - TTT*f2 - TTl*f3 - TTT*f3 + lll*f2*f3 + llT*f2*f3 + lTl*f2*f3 + lTT*f2*f3 + Tll*f2*f3 + TlT*f2*f3 + TTl*f2*f3 + TTT*f2*f3))/((f1 - r1)*(f2 - r2)*(f3 - r3)) - (r1*r2*r3*(lTT*f1 - TTT + TlT*f2 + TTT*f1 + TTT*f2 + TTl*f3 + TTT*f3 - llT*f1*f2 - lTT*f1*f2 - lTl*f1*f3 - lTT*f1*f3 - TlT*f1*f2 - Tll*f2*f3 - TlT*f2*f3 - TTT*f1*f2 - TTl*f1*f3 - TTT*f1*f3 - TTl*f2*f3 - TTT*f2*f3 + lll*f1*f2*f3 + llT*f1*f2*f3 + lTl*f1*f2*f3 + lTT*f1*f2*f3 + Tll*f1*f2*f3 + TlT*f1*f2*f3 + TTl*f1*f2*f3 + TTT*f1*f2*f3))/(pow((f1 - r1),2)*(f2 - r2)*(f3 - r3))),2) + pow(f2err,2)*pow(((r1*r2*r3*(lTT*f1 - TTT + TlT*f2 + TTT*f1 + TTT*f2 + TTl*f3 + TTT*f3 - llT*f1*f2 - lTT*f1*f2 - lTl*f1*f3 - lTT*f1*f3 - TlT*f1*f2 - Tll*f2*f3 - TlT*f2*f3 - TTT*f1*f2 - TTl*f1*f3 - TTT*f1*f3 - TTl*f2*f3 - TTT*f2*f3 + lll*f1*f2*f3 + llT*f1*f2*f3 + lTl*f1*f2*f3 + lTT*f1*f2*f3 + Tll*f1*f2*f3 + TlT*f1*f2*f3 + TTl*f1*f2*f3 + TTT*f1*f2*f3))/((f1 - r1)*pow((f2 - r2),2)*(f3 - r3)) - (r1*r2*r3*(TlT + TTT - llT*f1 - lTT*f1 - TlT*f1 - Tll*f3 - TlT*f3 - TTT*f1 - TTl*f3 - TTT*f3 + lll*f1*f3 + llT*f1*f3 + lTl*f1*f3 + lTT*f1*f3 + Tll*f1*f3 + TlT*f1*f3 + TTl*f1*f3 + TTT*f1*f3))/((f1 - r1)*(f2 - r2)*(f3 - r3))),2) + pow(f3err,2)*pow(((r1*r2*r3*(lTT*f1 - TTT + TlT*f2 + TTT*f1 + TTT*f2 + TTl*f3 + TTT*f3 - llT*f1*f2 - lTT*f1*f2 - lTl*f1*f3 - lTT*f1*f3 - TlT*f1*f2 - Tll*f2*f3 - TlT*f2*f3 - TTT*f1*f2 - TTl*f1*f3 - TTT*f1*f3 - TTl*f2*f3 - TTT*f2*f3 + lll*f1*f2*f3 + llT*f1*f2*f3 + lTl*f1*f2*f3 + lTT*f1*f2*f3 + Tll*f1*f2*f3 + TlT*f1*f2*f3 + TTl*f1*f2*f3 + TTT*f1*f2*f3))/((f1 - r1)*(f2 - r2)*pow((f3 - r3),2)) - (r1*r2*r3*(TTl + TTT - lTl*f1 - lTT*f1 - Tll*f2 - TlT*f2 - TTl*f1 - TTT*f1 - TTl*f2 - TTT*f2 + lll*f1*f2 + llT*f1*f2 + lTl*f1*f2 + lTT*f1*f2 + Tll*f1*f2 + TlT*f1*f2 + TTl*f1*f2 + TTT*f1*f2))/((f1 - r1)*(f2 - r2)*(f3 - r3))),2) + pow(r1err,2)*pow(((r2*r3*(lTT*f1 - TTT + TlT*f2 + TTT*f1 + TTT*f2 + TTl*f3 + TTT*f3 - llT*f1*f2 - lTT*f1*f2 - lTl*f1*f3 - lTT*f1*f3 - TlT*f1*f2 - Tll*f2*f3 - TlT*f2*f3 - TTT*f1*f2 - TTl*f1*f3 - TTT*f1*f3 - TTl*f2*f3 - TTT*f2*f3 + lll*f1*f2*f3 + llT*f1*f2*f3 + lTl*f1*f2*f3 + lTT*f1*f2*f3 + Tll*f1*f2*f3 + TlT*f1*f2*f3 + TTl*f1*f2*f3 + TTT*f1*f2*f3))/((f1 - r1)*(f2 - r2)*(f3 - r3)) + (r1*r2*r3*(lTT*f1 - TTT + TlT*f2 + TTT*f1 + TTT*f2 + TTl*f3 + TTT*f3 - llT*f1*f2 - lTT*f1*f2 - lTl*f1*f3 - lTT*f1*f3 - TlT*f1*f2 - Tll*f2*f3 - TlT*f2*f3 - TTT*f1*f2 - TTl*f1*f3 - TTT*f1*f3 - TTl*f2*f3 - TTT*f2*f3 + lll*f1*f2*f3 + llT*f1*f2*f3 + lTl*f1*f2*f3 + lTT*f1*f2*f3 + Tll*f1*f2*f3 + TlT*f1*f2*f3 + TTl*f1*f2*f3 + TTT*f1*f2*f3))/(pow((f1 - r1),2)*(f2 - r2)*(f3 - r3))),2) + pow(r2err,2)*pow(((r1*r3*(lTT*f1 - TTT + TlT*f2 + TTT*f1 + TTT*f2 + TTl*f3 + TTT*f3 - llT*f1*f2 - lTT*f1*f2 - lTl*f1*f3 - lTT*f1*f3 - TlT*f1*f2 - Tll*f2*f3 - TlT*f2*f3 - TTT*f1*f2 - TTl*f1*f3 - TTT*f1*f3 - TTl*f2*f3 - TTT*f2*f3 + lll*f1*f2*f3 + llT*f1*f2*f3 + lTl*f1*f2*f3 + lTT*f1*f2*f3 + Tll*f1*f2*f3 + TlT*f1*f2*f3 + TTl*f1*f2*f3 + TTT*f1*f2*f3))/((f1 - r1)*(f2 - r2)*(f3 - r3)) + (r1*r2*r3*(lTT*f1 - TTT + TlT*f2 + TTT*f1 + TTT*f2 + TTl*f3 + TTT*f3 - llT*f1*f2 - lTT*f1*f2 - lTl*f1*f3 - lTT*f1*f3 - TlT*f1*f2 - Tll*f2*f3 - TlT*f2*f3 - TTT*f1*f2 - TTl*f1*f3 - TTT*f1*f3 - TTl*f2*f3 - TTT*f2*f3 + lll*f1*f2*f3 + llT*f1*f2*f3 + lTl*f1*f2*f3 + lTT*f1*f2*f3 + Tll*f1*f2*f3 + TlT*f1*f2*f3 + TTl*f1*f2*f3 + TTT*f1*f2*f3))/((f1 - r1)*pow((f2 - r2),2)*(f3 - r3))),2) + pow(r3err,2)*pow(((r1*r2*(lTT*f1 - TTT + TlT*f2 + TTT*f1 + TTT*f2 + TTl*f3 + TTT*f3 - llT*f1*f2 - lTT*f1*f2 - lTl*f1*f3 - lTT*f1*f3 - TlT*f1*f2 - Tll*f2*f3 - TlT*f2*f3 - TTT*f1*f2 - TTl*f1*f3 - TTT*f1*f3 - TTl*f2*f3 - TTT*f2*f3 + lll*f1*f2*f3 + llT*f1*f2*f3 + lTl*f1*f2*f3 + lTT*f1*f2*f3 + Tll*f1*f2*f3 + TlT*f1*f2*f3 + TTl*f1*f2*f3 + TTT*f1*f2*f3))/((f1 - r1)*(f2 - r2)*(f3 - r3)) + (r1*r2*r3*(lTT*f1 - TTT + TlT*f2 + TTT*f1 + TTT*f2 + TTl*f3 + TTT*f3 - llT*f1*f2 - lTT*f1*f2 - lTl*f1*f3 - lTT*f1*f3 - TlT*f1*f2 - Tll*f2*f3 - TlT*f2*f3 - TTT*f1*f2 - TTl*f1*f3 - TTT*f1*f3 - TTl*f2*f3 - TTT*f2*f3 + lll*f1*f2*f3 + llT*f1*f2*f3 + lTl*f1*f2*f3 + lTT*f1*f2*f3 + Tll*f1*f2*f3 + TlT*f1*f2*f3 + TTl*f1*f2*f3 + TTT*f1*f2*f3))/((f1 - r1)*(f2 - r2)*pow((f3 - r3),2))),2) + (pow(lTTerr,2)*pow(r1,2)*pow(r2,2)*pow(r3,2)*pow((f1 - f1*f2 - f1*f3 + f1*f2*f3),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(TTTerr,2)*pow(r1,2)*pow(r2,2)*pow(r3,2)*pow((f1 + f2 + f3 - f1*f2 - f1*f3 - f2*f3 + f1*f2*f3 - 1),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(TlTerr,2)*pow(r1,2)*pow(r2,2)*pow(r3,2)*pow((f2 - f1*f2 - f2*f3 + f1*f2*f3),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(TTlerr,2)*pow(r1,2)*pow(r2,2)*pow(r3,2)*pow((f3 - f1*f3 - f2*f3 + f1*f2*f3),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(llTerr,2)*pow(r1,2)*pow(r2,2)*pow(r3,2)*pow((f1*f2 - f1*f2*f3),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(lTlerr,2)*pow(r1,2)*pow(r2,2)*pow(r3,2)*pow((f1*f3 - f1*f2*f3),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(Tllerr,2)*pow(r1,2)*pow(r2,2)*pow(r3,2)*pow((f2*f3 - f1*f2*f3),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(lllerr,2)*pow(f1,2)*pow(f2,2)*pow(f3,2)*pow(r1,2)*pow(r2,2)*pow(r3,2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2));
    return sqrt(res);
}

float MatrixMethod::N8_RRF_TTTerr(float r1,float f1,float r2,float f2,float r3,float f3, int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll, float r1err,float f1err,float r2err,float f2err,float r3err,float f3err){ 

float TTTerr = sqrt(TTT);
float llTerr = sqrt(llT);
float lTlerr = sqrt(lTl);
float Tllerr = sqrt(Tll);
float lTTerr = sqrt(lTT);
float TlTerr = sqrt(TlT);
float TTlerr = sqrt(TTl);
float lllerr = sqrt(lll);

  float res = pow(r3err,2)*pow(((f3*r1*r2*(lTT*f1 - TTT + TlT*f2 + TTT*f1 + TTT*f2 + TTl*r3 + TTT*r3 - llT*f1*f2 - lTT*f1*f2 - TlT*f1*f2 - TTT*f1*f2 - lTl*f1*r3 - lTT*f1*r3 - Tll*f2*r3 - TlT*f2*r3 - TTl*f1*r3 - TTT*f1*r3 - TTl*f2*r3 - TTT*f2*r3 + lll*f1*f2*r3 + llT*f1*f2*r3 + lTl*f1*f2*r3 + lTT*f1*f2*r3 + Tll*f1*f2*r3 + TlT*f1*f2*r3 + TTl*f1*f2*r3 + TTT*f1*f2*r3))/((f1 - r1)*(f2 - r2)*pow((f3 - r3),2)) + (f3*r1*r2*(TTl + TTT - lTl*f1 - lTT*f1 - Tll*f2 - TlT*f2 - TTl*f1 - TTT*f1 - TTl*f2 - TTT*f2 + lll*f1*f2 + llT*f1*f2 + lTl*f1*f2 + lTT*f1*f2 + Tll*f1*f2 + TlT*f1*f2 + TTl*f1*f2 + TTT*f1*f2))/((f1 - r1)*(f2 - r2)*(f3 - r3))),2) + pow(r1err,2)*pow(((f3*r2*(lTT*f1 - TTT + TlT*f2 + TTT*f1 + TTT*f2 + TTl*r3 + TTT*r3 - llT*f1*f2 - lTT*f1*f2 - TlT*f1*f2 - TTT*f1*f2 - lTl*f1*r3 - lTT*f1*r3 - Tll*f2*r3 - TlT*f2*r3 - TTl*f1*r3 - TTT*f1*r3 - TTl*f2*r3 - TTT*f2*r3 + lll*f1*f2*r3 + llT*f1*f2*r3 + lTl*f1*f2*r3 + lTT*f1*f2*r3 + Tll*f1*f2*r3 + TlT*f1*f2*r3 + TTl*f1*f2*r3 + TTT*f1*f2*r3))/((f1 - r1)*(f2 - r2)*(f3 - r3)) + (f3*r1*r2*(lTT*f1 - TTT + TlT*f2 + TTT*f1 + TTT*f2 + TTl*r3 + TTT*r3 - llT*f1*f2 - lTT*f1*f2 - TlT*f1*f2 - TTT*f1*f2 - lTl*f1*r3 - lTT*f1*r3 - Tll*f2*r3 - TlT*f2*r3 - TTl*f1*r3 - TTT*f1*r3 - TTl*f2*r3 - TTT*f2*r3 + lll*f1*f2*r3 + llT*f1*f2*r3 + lTl*f1*f2*r3 + lTT*f1*f2*r3 + Tll*f1*f2*r3 + TlT*f1*f2*r3 + TTl*f1*f2*r3 + TTT*f1*f2*r3))/(pow((f1 - r1),2)*(f2 - r2)*(f3 - r3))),2) + pow(r2err,2)*pow(((f3*r1*(lTT*f1 - TTT + TlT*f2 + TTT*f1 + TTT*f2 + TTl*r3 + TTT*r3 - llT*f1*f2 - lTT*f1*f2 - TlT*f1*f2 - TTT*f1*f2 - lTl*f1*r3 - lTT*f1*r3 - Tll*f2*r3 - TlT*f2*r3 - TTl*f1*r3 - TTT*f1*r3 - TTl*f2*r3 - TTT*f2*r3 + lll*f1*f2*r3 + llT*f1*f2*r3 + lTl*f1*f2*r3 + lTT*f1*f2*r3 + Tll*f1*f2*r3 + TlT*f1*f2*r3 + TTl*f1*f2*r3 + TTT*f1*f2*r3))/((f1 - r1)*(f2 - r2)*(f3 - r3)) + (f3*r1*r2*(lTT*f1 - TTT + TlT*f2 + TTT*f1 + TTT*f2 + TTl*r3 + TTT*r3 - llT*f1*f2 - lTT*f1*f2 - TlT*f1*f2 - TTT*f1*f2 - lTl*f1*r3 - lTT*f1*r3 - Tll*f2*r3 - TlT*f2*r3 - TTl*f1*r3 - TTT*f1*r3 - TTl*f2*r3 - TTT*f2*r3 + lll*f1*f2*r3 + llT*f1*f2*r3 + lTl*f1*f2*r3 + lTT*f1*f2*r3 + Tll*f1*f2*r3 + TlT*f1*f2*r3 + TTl*f1*f2*r3 + TTT*f1*f2*r3))/((f1 - r1)*pow((f2 - r2),2)*(f3 - r3))),2) + pow(f3err,2)*pow(((r1*r2*(lTT*f1 - TTT + TlT*f2 + TTT*f1 + TTT*f2 + TTl*r3 + TTT*r3 - llT*f1*f2 - lTT*f1*f2 - TlT*f1*f2 - TTT*f1*f2 - lTl*f1*r3 - lTT*f1*r3 - Tll*f2*r3 - TlT*f2*r3 - TTl*f1*r3 - TTT*f1*r3 - TTl*f2*r3 - TTT*f2*r3 + lll*f1*f2*r3 + llT*f1*f2*r3 + lTl*f1*f2*r3 + lTT*f1*f2*r3 + Tll*f1*f2*r3 + TlT*f1*f2*r3 + TTl*f1*f2*r3 + TTT*f1*f2*r3))/((f1 - r1)*(f2 - r2)*(f3 - r3)) - (f3*r1*r2*(lTT*f1 - TTT + TlT*f2 + TTT*f1 + TTT*f2 + TTl*r3 + TTT*r3 - llT*f1*f2 - lTT*f1*f2 - TlT*f1*f2 - TTT*f1*f2 - lTl*f1*r3 - lTT*f1*r3 - Tll*f2*r3 - TlT*f2*r3 - TTl*f1*r3 - TTT*f1*r3 - TTl*f2*r3 - TTT*f2*r3 + lll*f1*f2*r3 + llT*f1*f2*r3 + lTl*f1*f2*r3 + lTT*f1*f2*r3 + Tll*f1*f2*r3 + TlT*f1*f2*r3 + TTl*f1*f2*r3 + TTT*f1*f2*r3))/((f1 - r1)*(f2 - r2)*pow((f3 - r3),2))),2) + pow(f1err,2)*pow(((f3*r1*r2*(lTT*f1 - TTT + TlT*f2 + TTT*f1 + TTT*f2 + TTl*r3 + TTT*r3 - llT*f1*f2 - lTT*f1*f2 - TlT*f1*f2 - TTT*f1*f2 - lTl*f1*r3 - lTT*f1*r3 - Tll*f2*r3 - TlT*f2*r3 - TTl*f1*r3 - TTT*f1*r3 - TTl*f2*r3 - TTT*f2*r3 + lll*f1*f2*r3 + llT*f1*f2*r3 + lTl*f1*f2*r3 + lTT*f1*f2*r3 + Tll*f1*f2*r3 + TlT*f1*f2*r3 + TTl*f1*f2*r3 + TTT*f1*f2*r3))/(pow((f1 - r1),2)*(f2 - r2)*(f3 - r3)) - (f3*r1*r2*(lTT + TTT - llT*f2 - lTT*f2 - TlT*f2 - TTT*f2 - lTl*r3 - lTT*r3 - TTl*r3 - TTT*r3 + lll*f2*r3 + llT*f2*r3 + lTl*f2*r3 + lTT*f2*r3 + Tll*f2*r3 + TlT*f2*r3 + TTl*f2*r3 + TTT*f2*r3))/((f1 - r1)*(f2 - r2)*(f3 - r3))),2) + pow(f2err,2)*pow(((f3*r1*r2*(lTT*f1 - TTT + TlT*f2 + TTT*f1 + TTT*f2 + TTl*r3 + TTT*r3 - llT*f1*f2 - lTT*f1*f2 - TlT*f1*f2 - TTT*f1*f2 - lTl*f1*r3 - lTT*f1*r3 - Tll*f2*r3 - TlT*f2*r3 - TTl*f1*r3 - TTT*f1*r3 - TTl*f2*r3 - TTT*f2*r3 + lll*f1*f2*r3 + llT*f1*f2*r3 + lTl*f1*f2*r3 + lTT*f1*f2*r3 + Tll*f1*f2*r3 + TlT*f1*f2*r3 + TTl*f1*f2*r3 + TTT*f1*f2*r3))/((f1 - r1)*pow((f2 - r2),2)*(f3 - r3)) - (f3*r1*r2*(TlT + TTT - llT*f1 - lTT*f1 - TlT*f1 - TTT*f1 - Tll*r3 - TlT*r3 - TTl*r3 - TTT*r3 + lll*f1*r3 + llT*f1*r3 + lTl*f1*r3 + lTT*f1*r3 + Tll*f1*r3 + TlT*f1*r3 + TTl*f1*r3 + TTT*f1*r3))/((f1 - r1)*(f2 - r2)*(f3 - r3))),2) + (pow(lTTerr,2)*pow(f3,2)*pow(r1,2)*pow(r2,2)*pow((f1 - f1*f2 - f1*r3 + f1*f2*r3),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(TlTerr,2)*pow(f3,2)*pow(r1,2)*pow(r2,2)*pow((f2 - f1*f2 - f2*r3 + f1*f2*r3),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(TTTerr,2)*pow(f3,2)*pow(r1,2)*pow(r2,2)*pow((f1 + f2 + r3 - f1*f2 - f1*r3 - f2*r3 + f1*f2*r3 - 1),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(TTlerr,2)*pow(f3,2)*pow(r1,2)*pow(r2,2)*pow((r3 - f1*r3 - f2*r3 + f1*f2*r3),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(llTerr,2)*pow(f3,2)*pow(r1,2)*pow(r2,2)*pow((f1*f2 - f1*f2*r3),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(lTlerr,2)*pow(f3,2)*pow(r1,2)*pow(r2,2)*pow((f1*r3 - f1*f2*r3),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(Tllerr,2)*pow(f3,2)*pow(r1,2)*pow(r2,2)*pow((f2*r3 - f1*f2*r3),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(lllerr,2)*pow(f1,2)*pow(f2,2)*pow(f3,2)*pow(r1,2)*pow(r2,2)*pow(r3,2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2));
  return sqrt(res);
}
 
float MatrixMethod::N8_RFR_TTTerr(float r1,float f1,float r2,float f2,float r3,float f3, int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll, float r1err,float f1err,float r2err,float f2err,float r3err,float f3err){ 

float TTTerr = sqrt(TTT);
float llTerr = sqrt(llT);
float lTlerr = sqrt(lTl);
float Tllerr = sqrt(Tll);
float lTTerr = sqrt(lTT);
float TlTerr = sqrt(TlT);
float TTlerr = sqrt(TTl);
float lllerr = sqrt(lll);

  float res = pow(r2err,2)*pow(((f2*r1*r3*(lTT*f1 - TTT + TTT*f1 + TTl*f3 + TTT*f3 + TlT*r2 + TTT*r2 - lTl*f1*f3 - lTT*f1*f3 - TTl*f1*f3 - TTT*f1*f3 - llT*f1*r2 - lTT*f1*r2 - TlT*f1*r2 - Tll*f3*r2 - TlT*f3*r2 - TTT*f1*r2 - TTl*f3*r2 - TTT*f3*r2 + lll*f1*f3*r2 + llT*f1*f3*r2 + lTl*f1*f3*r2 + lTT*f1*f3*r2 + Tll*f1*f3*r2 + TlT*f1*f3*r2 + TTl*f1*f3*r2 + TTT*f1*f3*r2))/((f1 - r1)*pow((f2 - r2),2)*(f3 - r3)) + (f2*r1*r3*(TlT + TTT - llT*f1 - lTT*f1 - TlT*f1 - Tll*f3 - TlT*f3 - TTT*f1 - TTl*f3 - TTT*f3 + lll*f1*f3 + llT*f1*f3 + lTl*f1*f3 + lTT*f1*f3 + Tll*f1*f3 + TlT*f1*f3 + TTl*f1*f3 + TTT*f1*f3))/((f1 - r1)*(f2 - r2)*(f3 - r3))),2) + pow(r1err,2)*pow(((f2*r3*(lTT*f1 - TTT + TTT*f1 + TTl*f3 + TTT*f3 + TlT*r2 + TTT*r2 - lTl*f1*f3 - lTT*f1*f3 - TTl*f1*f3 - TTT*f1*f3 - llT*f1*r2 - lTT*f1*r2 - TlT*f1*r2 - Tll*f3*r2 - TlT*f3*r2 - TTT*f1*r2 - TTl*f3*r2 - TTT*f3*r2 + lll*f1*f3*r2 + llT*f1*f3*r2 + lTl*f1*f3*r2 + lTT*f1*f3*r2 + Tll*f1*f3*r2 + TlT*f1*f3*r2 + TTl*f1*f3*r2 + TTT*f1*f3*r2))/((f1 - r1)*(f2 - r2)*(f3 - r3)) + (f2*r1*r3*(lTT*f1 - TTT + TTT*f1 + TTl*f3 + TTT*f3 + TlT*r2 + TTT*r2 - lTl*f1*f3 - lTT*f1*f3 - TTl*f1*f3 - TTT*f1*f3 - llT*f1*r2 - lTT*f1*r2 - TlT*f1*r2 - Tll*f3*r2 - TlT*f3*r2 - TTT*f1*r2 - TTl*f3*r2 - TTT*f3*r2 + lll*f1*f3*r2 + llT*f1*f3*r2 + lTl*f1*f3*r2 + lTT*f1*f3*r2 + Tll*f1*f3*r2 + TlT*f1*f3*r2 + TTl*f1*f3*r2 + TTT*f1*f3*r2))/(pow((f1 - r1),2)*(f2 - r2)*(f3 - r3))),2) + pow(r3err,2)*pow(((f2*r1*(lTT*f1 - TTT + TTT*f1 + TTl*f3 + TTT*f3 + TlT*r2 + TTT*r2 - lTl*f1*f3 - lTT*f1*f3 - TTl*f1*f3 - TTT*f1*f3 - llT*f1*r2 - lTT*f1*r2 - TlT*f1*r2 - Tll*f3*r2 - TlT*f3*r2 - TTT*f1*r2 - TTl*f3*r2 - TTT*f3*r2 + lll*f1*f3*r2 + llT*f1*f3*r2 + lTl*f1*f3*r2 + lTT*f1*f3*r2 + Tll*f1*f3*r2 + TlT*f1*f3*r2 + TTl*f1*f3*r2 + TTT*f1*f3*r2))/((f1 - r1)*(f2 - r2)*(f3 - r3)) + (f2*r1*r3*(lTT*f1 - TTT + TTT*f1 + TTl*f3 + TTT*f3 + TlT*r2 + TTT*r2 - lTl*f1*f3 - lTT*f1*f3 - TTl*f1*f3 - TTT*f1*f3 - llT*f1*r2 - lTT*f1*r2 - TlT*f1*r2 - Tll*f3*r2 - TlT*f3*r2 - TTT*f1*r2 - TTl*f3*r2 - TTT*f3*r2 + lll*f1*f3*r2 + llT*f1*f3*r2 + lTl*f1*f3*r2 + lTT*f1*f3*r2 + Tll*f1*f3*r2 + TlT*f1*f3*r2 + TTl*f1*f3*r2 + TTT*f1*f3*r2))/((f1 - r1)*(f2 - r2)*pow((f3 - r3),2))),2) + pow(f2err,2)*pow(((r1*r3*(lTT*f1 - TTT + TTT*f1 + TTl*f3 + TTT*f3 + TlT*r2 + TTT*r2 - lTl*f1*f3 - lTT*f1*f3 - TTl*f1*f3 - TTT*f1*f3 - llT*f1*r2 - lTT*f1*r2 - TlT*f1*r2 - Tll*f3*r2 - TlT*f3*r2 - TTT*f1*r2 - TTl*f3*r2 - TTT*f3*r2 + lll*f1*f3*r2 + llT*f1*f3*r2 + lTl*f1*f3*r2 + lTT*f1*f3*r2 + Tll*f1*f3*r2 + TlT*f1*f3*r2 + TTl*f1*f3*r2 + TTT*f1*f3*r2))/((f1 - r1)*(f2 - r2)*(f3 - r3)) - (f2*r1*r3*(lTT*f1 - TTT + TTT*f1 + TTl*f3 + TTT*f3 + TlT*r2 + TTT*r2 - lTl*f1*f3 - lTT*f1*f3 - TTl*f1*f3 - TTT*f1*f3 - llT*f1*r2 - lTT*f1*r2 - TlT*f1*r2 - Tll*f3*r2 - TlT*f3*r2 - TTT*f1*r2 - TTl*f3*r2 - TTT*f3*r2 + lll*f1*f3*r2 + llT*f1*f3*r2 + lTl*f1*f3*r2 + lTT*f1*f3*r2 + Tll*f1*f3*r2 + TlT*f1*f3*r2 + TTl*f1*f3*r2 + TTT*f1*f3*r2))/((f1 - r1)*pow((f2 - r2),2)*(f3 - r3))),2) + pow(f1err,2)*pow(((f2*r1*r3*(lTT*f1 - TTT + TTT*f1 + TTl*f3 + TTT*f3 + TlT*r2 + TTT*r2 - lTl*f1*f3 - lTT*f1*f3 - TTl*f1*f3 - TTT*f1*f3 - llT*f1*r2 - lTT*f1*r2 - TlT*f1*r2 - Tll*f3*r2 - TlT*f3*r2 - TTT*f1*r2 - TTl*f3*r2 - TTT*f3*r2 + lll*f1*f3*r2 + llT*f1*f3*r2 + lTl*f1*f3*r2 + lTT*f1*f3*r2 + Tll*f1*f3*r2 + TlT*f1*f3*r2 + TTl*f1*f3*r2 + TTT*f1*f3*r2))/(pow((f1 - r1),2)*(f2 - r2)*(f3 - r3)) - (f2*r1*r3*(lTT + TTT - lTl*f3 - lTT*f3 - TTl*f3 - TTT*f3 - llT*r2 - lTT*r2 - TlT*r2 - TTT*r2 + lll*f3*r2 + llT*f3*r2 + lTl*f3*r2 + lTT*f3*r2 + Tll*f3*r2 + TlT*f3*r2 + TTl*f3*r2 + TTT*f3*r2))/((f1 - r1)*(f2 - r2)*(f3 - r3))),2) + pow(f3err,2)*pow(((f2*r1*r3*(lTT*f1 - TTT + TTT*f1 + TTl*f3 + TTT*f3 + TlT*r2 + TTT*r2 - lTl*f1*f3 - lTT*f1*f3 - TTl*f1*f3 - TTT*f1*f3 - llT*f1*r2 - lTT*f1*r2 - TlT*f1*r2 - Tll*f3*r2 - TlT*f3*r2 - TTT*f1*r2 - TTl*f3*r2 - TTT*f3*r2 + lll*f1*f3*r2 + llT*f1*f3*r2 + lTl*f1*f3*r2 + lTT*f1*f3*r2 + Tll*f1*f3*r2 + TlT*f1*f3*r2 + TTl*f1*f3*r2 + TTT*f1*f3*r2))/((f1 - r1)*(f2 - r2)*pow((f3 - r3),2)) - (f2*r1*r3*(TTl + TTT - lTl*f1 - lTT*f1 - TTl*f1 - TTT*f1 - Tll*r2 - TlT*r2 - TTl*r2 - TTT*r2 + lll*f1*r2 + llT*f1*r2 + lTl*f1*r2 + lTT*f1*r2 + Tll*f1*r2 + TlT*f1*r2 + TTl*f1*r2 + TTT*f1*r2))/((f1 - r1)*(f2 - r2)*(f3 - r3))),2) + (pow(lTTerr,2)*pow(f2,2)*pow(r1,2)*pow(r3,2)*pow((f1 - f1*f3 - f1*r2 + f1*f3*r2),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(TTlerr,2)*pow(f2,2)*pow(r1,2)*pow(r3,2)*pow((f3 - f1*f3 - f3*r2 + f1*f3*r2),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(TTTerr,2)*pow(f2,2)*pow(r1,2)*pow(r3,2)*pow((f1 + f3 + r2 - f1*f3 - f1*r2 - f3*r2 + f1*f3*r2 - 1),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(TlTerr,2)*pow(f2,2)*pow(r1,2)*pow(r3,2)*pow((r2 - f1*r2 - f3*r2 + f1*f3*r2),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(lTlerr,2)*pow(f2,2)*pow(r1,2)*pow(r3,2)*pow((f1*f3 - f1*f3*r2),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(llTerr,2)*pow(f2,2)*pow(r1,2)*pow(r3,2)*pow((f1*r2 - f1*f3*r2),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(Tllerr,2)*pow(f2,2)*pow(r1,2)*pow(r3,2)*pow((f3*r2 - f1*f3*r2),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(lllerr,2)*pow(f1,2)*pow(f2,2)*pow(f3,2)*pow(r1,2)*pow(r2,2)*pow(r3,2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2));
  return sqrt(res);
}
 
float MatrixMethod::N8_FRR_TTTerr(float r1,float f1,float r2,float f2,float r3,float f3, int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll, float r1err,float f1err,float r2err,float f2err,float r3err,float f3err){ 

float TTTerr = sqrt(TTT);
float llTerr = sqrt(llT);
float lTlerr = sqrt(lTl);
float Tllerr = sqrt(Tll);
float lTTerr = sqrt(lTT);
float TlTerr = sqrt(TlT);
float TTlerr = sqrt(TTl);
float lllerr = sqrt(lll);

  float res = pow(r1err,2)*pow(((f1*r2*r3*(TlT*f2 - TTT + TTT*f2 + TTl*f3 + TTT*f3 + lTT*r1 + TTT*r1 - Tll*f2*f3 - TlT*f2*f3 - TTl*f2*f3 - TTT*f2*f3 - llT*f2*r1 - lTT*f2*r1 - lTl*f3*r1 - lTT*f3*r1 - TlT*f2*r1 - TTT*f2*r1 - TTl*f3*r1 - TTT*f3*r1 + lll*f2*f3*r1 + llT*f2*f3*r1 + lTl*f2*f3*r1 + lTT*f2*f3*r1 + Tll*f2*f3*r1 + TlT*f2*f3*r1 + TTl*f2*f3*r1 + TTT*f2*f3*r1))/(pow((f1 - r1),2)*(f2 - r2)*(f3 - r3)) + (f1*r2*r3*(lTT + TTT - llT*f2 - lTT*f2 - lTl*f3 - lTT*f3 - TlT*f2 - TTT*f2 - TTl*f3 - TTT*f3 + lll*f2*f3 + llT*f2*f3 + lTl*f2*f3 + lTT*f2*f3 + Tll*f2*f3 + TlT*f2*f3 + TTl*f2*f3 + TTT*f2*f3))/((f1 - r1)*(f2 - r2)*(f3 - r3))),2) + pow(r2err,2)*pow(((f1*r3*(TlT*f2 - TTT + TTT*f2 + TTl*f3 + TTT*f3 + lTT*r1 + TTT*r1 - Tll*f2*f3 - TlT*f2*f3 - TTl*f2*f3 - TTT*f2*f3 - llT*f2*r1 - lTT*f2*r1 - lTl*f3*r1 - lTT*f3*r1 - TlT*f2*r1 - TTT*f2*r1 - TTl*f3*r1 - TTT*f3*r1 + lll*f2*f3*r1 + llT*f2*f3*r1 + lTl*f2*f3*r1 + lTT*f2*f3*r1 + Tll*f2*f3*r1 + TlT*f2*f3*r1 + TTl*f2*f3*r1 + TTT*f2*f3*r1))/((f1 - r1)*(f2 - r2)*(f3 - r3)) + (f1*r2*r3*(TlT*f2 - TTT + TTT*f2 + TTl*f3 + TTT*f3 + lTT*r1 + TTT*r1 - Tll*f2*f3 - TlT*f2*f3 - TTl*f2*f3 - TTT*f2*f3 - llT*f2*r1 - lTT*f2*r1 - lTl*f3*r1 - lTT*f3*r1 - TlT*f2*r1 - TTT*f2*r1 - TTl*f3*r1 - TTT*f3*r1 + lll*f2*f3*r1 + llT*f2*f3*r1 + lTl*f2*f3*r1 + lTT*f2*f3*r1 + Tll*f2*f3*r1 + TlT*f2*f3*r1 + TTl*f2*f3*r1 + TTT*f2*f3*r1))/((f1 - r1)*pow((f2 - r2),2)*(f3 - r3))),2) + pow(r3err,2)*pow(((f1*r2*(TlT*f2 - TTT + TTT*f2 + TTl*f3 + TTT*f3 + lTT*r1 + TTT*r1 - Tll*f2*f3 - TlT*f2*f3 - TTl*f2*f3 - TTT*f2*f3 - llT*f2*r1 - lTT*f2*r1 - lTl*f3*r1 - lTT*f3*r1 - TlT*f2*r1 - TTT*f2*r1 - TTl*f3*r1 - TTT*f3*r1 + lll*f2*f3*r1 + llT*f2*f3*r1 + lTl*f2*f3*r1 + lTT*f2*f3*r1 + Tll*f2*f3*r1 + TlT*f2*f3*r1 + TTl*f2*f3*r1 + TTT*f2*f3*r1))/((f1 - r1)*(f2 - r2)*(f3 - r3)) + (f1*r2*r3*(TlT*f2 - TTT + TTT*f2 + TTl*f3 + TTT*f3 + lTT*r1 + TTT*r1 - Tll*f2*f3 - TlT*f2*f3 - TTl*f2*f3 - TTT*f2*f3 - llT*f2*r1 - lTT*f2*r1 - lTl*f3*r1 - lTT*f3*r1 - TlT*f2*r1 - TTT*f2*r1 - TTl*f3*r1 - TTT*f3*r1 + lll*f2*f3*r1 + llT*f2*f3*r1 + lTl*f2*f3*r1 + lTT*f2*f3*r1 + Tll*f2*f3*r1 + TlT*f2*f3*r1 + TTl*f2*f3*r1 + TTT*f2*f3*r1))/((f1 - r1)*(f2 - r2)*pow((f3 - r3),2))),2) + pow(f1err,2)*pow(((r2*r3*(TlT*f2 - TTT + TTT*f2 + TTl*f3 + TTT*f3 + lTT*r1 + TTT*r1 - Tll*f2*f3 - TlT*f2*f3 - TTl*f2*f3 - TTT*f2*f3 - llT*f2*r1 - lTT*f2*r1 - lTl*f3*r1 - lTT*f3*r1 - TlT*f2*r1 - TTT*f2*r1 - TTl*f3*r1 - TTT*f3*r1 + lll*f2*f3*r1 + llT*f2*f3*r1 + lTl*f2*f3*r1 + lTT*f2*f3*r1 + Tll*f2*f3*r1 + TlT*f2*f3*r1 + TTl*f2*f3*r1 + TTT*f2*f3*r1))/((f1 - r1)*(f2 - r2)*(f3 - r3)) - (f1*r2*r3*(TlT*f2 - TTT + TTT*f2 + TTl*f3 + TTT*f3 + lTT*r1 + TTT*r1 - Tll*f2*f3 - TlT*f2*f3 - TTl*f2*f3 - TTT*f2*f3 - llT*f2*r1 - lTT*f2*r1 - lTl*f3*r1 - lTT*f3*r1 - TlT*f2*r1 - TTT*f2*r1 - TTl*f3*r1 - TTT*f3*r1 + lll*f2*f3*r1 + llT*f2*f3*r1 + lTl*f2*f3*r1 + lTT*f2*f3*r1 + Tll*f2*f3*r1 + TlT*f2*f3*r1 + TTl*f2*f3*r1 + TTT*f2*f3*r1))/(pow((f1 - r1),2)*(f2 - r2)*(f3 - r3))),2) + pow(f2err,2)*pow(((f1*r2*r3*(TlT*f2 - TTT + TTT*f2 + TTl*f3 + TTT*f3 + lTT*r1 + TTT*r1 - Tll*f2*f3 - TlT*f2*f3 - TTl*f2*f3 - TTT*f2*f3 - llT*f2*r1 - lTT*f2*r1 - lTl*f3*r1 - lTT*f3*r1 - TlT*f2*r1 - TTT*f2*r1 - TTl*f3*r1 - TTT*f3*r1 + lll*f2*f3*r1 + llT*f2*f3*r1 + lTl*f2*f3*r1 + lTT*f2*f3*r1 + Tll*f2*f3*r1 + TlT*f2*f3*r1 + TTl*f2*f3*r1 + TTT*f2*f3*r1))/((f1 - r1)*pow((f2 - r2),2)*(f3 - r3)) - (f1*r2*r3*(TlT + TTT - Tll*f3 - TlT*f3 - TTl*f3 - TTT*f3 - llT*r1 - lTT*r1 - TlT*r1 - TTT*r1 + lll*f3*r1 + llT*f3*r1 + lTl*f3*r1 + lTT*f3*r1 + Tll*f3*r1 + TlT*f3*r1 + TTl*f3*r1 + TTT*f3*r1))/((f1 - r1)*(f2 - r2)*(f3 - r3))),2) + pow(f3err,2)*pow(((f1*r2*r3*(TlT*f2 - TTT + TTT*f2 + TTl*f3 + TTT*f3 + lTT*r1 + TTT*r1 - Tll*f2*f3 - TlT*f2*f3 - TTl*f2*f3 - TTT*f2*f3 - llT*f2*r1 - lTT*f2*r1 - lTl*f3*r1 - lTT*f3*r1 - TlT*f2*r1 - TTT*f2*r1 - TTl*f3*r1 - TTT*f3*r1 + lll*f2*f3*r1 + llT*f2*f3*r1 + lTl*f2*f3*r1 + lTT*f2*f3*r1 + Tll*f2*f3*r1 + TlT*f2*f3*r1 + TTl*f2*f3*r1 + TTT*f2*f3*r1))/((f1 - r1)*(f2 - r2)*pow((f3 - r3),2)) - (f1*r2*r3*(TTl + TTT - Tll*f2 - TlT*f2 - TTl*f2 - TTT*f2 - lTl*r1 - lTT*r1 - TTl*r1 - TTT*r1 + lll*f2*r1 + llT*f2*r1 + lTl*f2*r1 + lTT*f2*r1 + Tll*f2*r1 + TlT*f2*r1 + TTl*f2*r1 + TTT*f2*r1))/((f1 - r1)*(f2 - r2)*(f3 - r3))),2) + (pow(TlTerr,2)*pow(f1,2)*pow(r2,2)*pow(r3,2)*pow((f2 - f2*f3 - f2*r1 + f2*f3*r1),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(TTlerr,2)*pow(f1,2)*pow(r2,2)*pow(r3,2)*pow((f3 - f2*f3 - f3*r1 + f2*f3*r1),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(lTTerr,2)*pow(f1,2)*pow(r2,2)*pow(r3,2)*pow((r1 - f2*r1 - f3*r1 + f2*f3*r1),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(TTTerr,2)*pow(f1,2)*pow(r2,2)*pow(r3,2)*pow((f2 + f3 + r1 - f2*f3 - f2*r1 - f3*r1 + f2*f3*r1 - 1),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(Tllerr,2)*pow(f1,2)*pow(r2,2)*pow(r3,2)*pow((f2*f3 - f2*f3*r1),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(llTerr,2)*pow(f1,2)*pow(r2,2)*pow(r3,2)*pow((f2*r1 - f2*f3*r1),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(lTlerr,2)*pow(f1,2)*pow(r2,2)*pow(r3,2)*pow((f3*r1 - f2*f3*r1),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(lllerr,2)*pow(f1,2)*pow(f2,2)*pow(f3,2)*pow(r1,2)*pow(r2,2)*pow(r3,2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2));
  return sqrt(res);
}
 
float MatrixMethod::N8_RFF_TTTerr(float r1,float f1,float r2,float f2,float r3,float f3, int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll, float r1err,float f1err,float r2err,float f2err,float r3err,float f3err){ 

float TTTerr = sqrt(TTT);
float llTerr = sqrt(llT);
float lTlerr = sqrt(lTl);
float Tllerr = sqrt(Tll);
float lTTerr = sqrt(lTT);
float TlTerr = sqrt(TlT);
float TTlerr = sqrt(TTl);
float lllerr = sqrt(lll);

  float res = pow(r1err,2)*pow(((f2*f3*(lTT*f1 - TTT + TTT*f1 + TlT*r2 + TTT*r2 + TTl*r3 + TTT*r3 - llT*f1*r2 - lTT*f1*r2 - lTl*f1*r3 - lTT*f1*r3 - TlT*f1*r2 - TTT*f1*r2 - TTl*f1*r3 - TTT*f1*r3 - Tll*r2*r3 - TlT*r2*r3 - TTl*r2*r3 - TTT*r2*r3 + lll*f1*r2*r3 + llT*f1*r2*r3 + lTl*f1*r2*r3 + lTT*f1*r2*r3 + Tll*f1*r2*r3 + TlT*f1*r2*r3 + TTl*f1*r2*r3 + TTT*f1*r2*r3))/((f1 - r1)*(f2 - r2)*(f3 - r3)) + (f2*f3*r1*(lTT*f1 - TTT + TTT*f1 + TlT*r2 + TTT*r2 + TTl*r3 + TTT*r3 - llT*f1*r2 - lTT*f1*r2 - lTl*f1*r3 - lTT*f1*r3 - TlT*f1*r2 - TTT*f1*r2 - TTl*f1*r3 - TTT*f1*r3 - Tll*r2*r3 - TlT*r2*r3 - TTl*r2*r3 - TTT*r2*r3 + lll*f1*r2*r3 + llT*f1*r2*r3 + lTl*f1*r2*r3 + lTT*f1*r2*r3 + Tll*f1*r2*r3 + TlT*f1*r2*r3 + TTl*f1*r2*r3 + TTT*f1*r2*r3))/(pow((f1 - r1),2)*(f2 - r2)*(f3 - r3))),2) + pow(f2err,2)*pow(((f3*r1*(lTT*f1 - TTT + TTT*f1 + TlT*r2 + TTT*r2 + TTl*r3 + TTT*r3 - llT*f1*r2 - lTT*f1*r2 - lTl*f1*r3 - lTT*f1*r3 - TlT*f1*r2 - TTT*f1*r2 - TTl*f1*r3 - TTT*f1*r3 - Tll*r2*r3 - TlT*r2*r3 - TTl*r2*r3 - TTT*r2*r3 + lll*f1*r2*r3 + llT*f1*r2*r3 + lTl*f1*r2*r3 + lTT*f1*r2*r3 + Tll*f1*r2*r3 + TlT*f1*r2*r3 + TTl*f1*r2*r3 + TTT*f1*r2*r3))/((f1 - r1)*(f2 - r2)*(f3 - r3)) - (f2*f3*r1*(lTT*f1 - TTT + TTT*f1 + TlT*r2 + TTT*r2 + TTl*r3 + TTT*r3 - llT*f1*r2 - lTT*f1*r2 - lTl*f1*r3 - lTT*f1*r3 - TlT*f1*r2 - TTT*f1*r2 - TTl*f1*r3 - TTT*f1*r3 - Tll*r2*r3 - TlT*r2*r3 - TTl*r2*r3 - TTT*r2*r3 + lll*f1*r2*r3 + llT*f1*r2*r3 + lTl*f1*r2*r3 + lTT*f1*r2*r3 + Tll*f1*r2*r3 + TlT*f1*r2*r3 + TTl*f1*r2*r3 + TTT*f1*r2*r3))/((f1 - r1)*pow((f2 - r2),2)*(f3 - r3))),2) + pow(f3err,2)*pow(((f2*r1*(lTT*f1 - TTT + TTT*f1 + TlT*r2 + TTT*r2 + TTl*r3 + TTT*r3 - llT*f1*r2 - lTT*f1*r2 - lTl*f1*r3 - lTT*f1*r3 - TlT*f1*r2 - TTT*f1*r2 - TTl*f1*r3 - TTT*f1*r3 - Tll*r2*r3 - TlT*r2*r3 - TTl*r2*r3 - TTT*r2*r3 + lll*f1*r2*r3 + llT*f1*r2*r3 + lTl*f1*r2*r3 + lTT*f1*r2*r3 + Tll*f1*r2*r3 + TlT*f1*r2*r3 + TTl*f1*r2*r3 + TTT*f1*r2*r3))/((f1 - r1)*(f2 - r2)*(f3 - r3)) - (f2*f3*r1*(lTT*f1 - TTT + TTT*f1 + TlT*r2 + TTT*r2 + TTl*r3 + TTT*r3 - llT*f1*r2 - lTT*f1*r2 - lTl*f1*r3 - lTT*f1*r3 - TlT*f1*r2 - TTT*f1*r2 - TTl*f1*r3 - TTT*f1*r3 - Tll*r2*r3 - TlT*r2*r3 - TTl*r2*r3 - TTT*r2*r3 + lll*f1*r2*r3 + llT*f1*r2*r3 + lTl*f1*r2*r3 + lTT*f1*r2*r3 + Tll*f1*r2*r3 + TlT*f1*r2*r3 + TTl*f1*r2*r3 + TTT*f1*r2*r3))/((f1 - r1)*(f2 - r2)*pow((f3 - r3),2))),2) + pow(r2err,2)*pow(((f2*f3*r1*(TlT + TTT - llT*f1 - lTT*f1 - TlT*f1 - TTT*f1 - Tll*r3 - TlT*r3 - TTl*r3 - TTT*r3 + lll*f1*r3 + llT*f1*r3 + lTl*f1*r3 + lTT*f1*r3 + Tll*f1*r3 + TlT*f1*r3 + TTl*f1*r3 + TTT*f1*r3))/((f1 - r1)*(f2 - r2)*(f3 - r3)) + (f2*f3*r1*(lTT*f1 - TTT + TTT*f1 + TlT*r2 + TTT*r2 + TTl*r3 + TTT*r3 - llT*f1*r2 - lTT*f1*r2 - lTl*f1*r3 - lTT*f1*r3 - TlT*f1*r2 - TTT*f1*r2 - TTl*f1*r3 - TTT*f1*r3 - Tll*r2*r3 - TlT*r2*r3 - TTl*r2*r3 - TTT*r2*r3 + lll*f1*r2*r3 + llT*f1*r2*r3 + lTl*f1*r2*r3 + lTT*f1*r2*r3 + Tll*f1*r2*r3 + TlT*f1*r2*r3 + TTl*f1*r2*r3 + TTT*f1*r2*r3))/((f1 - r1)*pow((f2 - r2),2)*(f3 - r3))),2) + pow(r3err,2)*pow(((f2*f3*r1*(TTl + TTT - lTl*f1 - lTT*f1 - TTl*f1 - TTT*f1 - Tll*r2 - TlT*r2 - TTl*r2 - TTT*r2 + lll*f1*r2 + llT*f1*r2 + lTl*f1*r2 + lTT*f1*r2 + Tll*f1*r2 + TlT*f1*r2 + TTl*f1*r2 + TTT*f1*r2))/((f1 - r1)*(f2 - r2)*(f3 - r3)) + (f2*f3*r1*(lTT*f1 - TTT + TTT*f1 + TlT*r2 + TTT*r2 + TTl*r3 + TTT*r3 - llT*f1*r2 - lTT*f1*r2 - lTl*f1*r3 - lTT*f1*r3 - TlT*f1*r2 - TTT*f1*r2 - TTl*f1*r3 - TTT*f1*r3 - Tll*r2*r3 - TlT*r2*r3 - TTl*r2*r3 - TTT*r2*r3 + lll*f1*r2*r3 + llT*f1*r2*r3 + lTl*f1*r2*r3 + lTT*f1*r2*r3 + Tll*f1*r2*r3 + TlT*f1*r2*r3 + TTl*f1*r2*r3 + TTT*f1*r2*r3))/((f1 - r1)*(f2 - r2)*pow((f3 - r3),2))),2) + pow(f1err,2)*pow(((f2*f3*r1*(lTT + TTT - llT*r2 - lTT*r2 - lTl*r3 - lTT*r3 - TlT*r2 - TTT*r2 - TTl*r3 - TTT*r3 + lll*r2*r3 + llT*r2*r3 + lTl*r2*r3 + lTT*r2*r3 + Tll*r2*r3 + TlT*r2*r3 + TTl*r2*r3 + TTT*r2*r3))/((f1 - r1)*(f2 - r2)*(f3 - r3)) - (f2*f3*r1*(lTT*f1 - TTT + TTT*f1 + TlT*r2 + TTT*r2 + TTl*r3 + TTT*r3 - llT*f1*r2 - lTT*f1*r2 - lTl*f1*r3 - lTT*f1*r3 - TlT*f1*r2 - TTT*f1*r2 - TTl*f1*r3 - TTT*f1*r3 - Tll*r2*r3 - TlT*r2*r3 - TTl*r2*r3 - TTT*r2*r3 + lll*f1*r2*r3 + llT*f1*r2*r3 + lTl*f1*r2*r3 + lTT*f1*r2*r3 + Tll*f1*r2*r3 + TlT*f1*r2*r3 + TTl*f1*r2*r3 + TTT*f1*r2*r3))/(pow((f1 - r1),2)*(f2 - r2)*(f3 - r3))),2) + (pow(Tllerr,2)*pow(f2,2)*pow(f3,2)*pow(r1,2)*pow((r2*r3 - f1*r2*r3),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(lTTerr,2)*pow(f2,2)*pow(f3,2)*pow(r1,2)*pow((f1 - f1*r2 - f1*r3 + f1*r2*r3),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(TlTerr,2)*pow(f2,2)*pow(f3,2)*pow(r1,2)*pow((r2 - f1*r2 - r2*r3 + f1*r2*r3),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(TTlerr,2)*pow(f2,2)*pow(f3,2)*pow(r1,2)*pow((r3 - f1*r3 - r2*r3 + f1*r2*r3),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(TTTerr,2)*pow(f2,2)*pow(f3,2)*pow(r1,2)*pow((f1 + r2 + r3 - f1*r2 - f1*r3 - r2*r3 + f1*r2*r3 - 1),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(llTerr,2)*pow(f2,2)*pow(f3,2)*pow(r1,2)*pow((f1*r2 - f1*r2*r3),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(lTlerr,2)*pow(f2,2)*pow(f3,2)*pow(r1,2)*pow((f1*r3 - f1*r2*r3),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(lllerr,2)*pow(f1,2)*pow(f2,2)*pow(f3,2)*pow(r1,2)*pow(r2,2)*pow(r3,2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2));
  return sqrt(res);
}
 
float MatrixMethod::N8_FRF_TTTerr(float r1,float f1,float r2,float f2,float r3,float f3, int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll, float r1err,float f1err,float r2err,float f2err,float r3err,float f3err){ 

float TTTerr = sqrt(TTT);
float llTerr = sqrt(llT);
float lTlerr = sqrt(lTl);
float Tllerr = sqrt(Tll);
float lTTerr = sqrt(lTT);
float TlTerr = sqrt(TlT);
float TTlerr = sqrt(TTl);
float lllerr = sqrt(lll);

  float res = pow(r2err,2)*pow(((f1*f3*(TlT*f2 - TTT + TTT*f2 + lTT*r1 + TTT*r1 + TTl*r3 + TTT*r3 - llT*f2*r1 - lTT*f2*r1 - TlT*f2*r1 - Tll*f2*r3 - TlT*f2*r3 - TTT*f2*r1 - TTl*f2*r3 - TTT*f2*r3 - lTl*r1*r3 - lTT*r1*r3 - TTl*r1*r3 - TTT*r1*r3 + lll*f2*r1*r3 + llT*f2*r1*r3 + lTl*f2*r1*r3 + lTT*f2*r1*r3 + Tll*f2*r1*r3 + TlT*f2*r1*r3 + TTl*f2*r1*r3 + TTT*f2*r1*r3))/((f1 - r1)*(f2 - r2)*(f3 - r3)) + (f1*f3*r2*(TlT*f2 - TTT + TTT*f2 + lTT*r1 + TTT*r1 + TTl*r3 + TTT*r3 - llT*f2*r1 - lTT*f2*r1 - TlT*f2*r1 - Tll*f2*r3 - TlT*f2*r3 - TTT*f2*r1 - TTl*f2*r3 - TTT*f2*r3 - lTl*r1*r3 - lTT*r1*r3 - TTl*r1*r3 - TTT*r1*r3 + lll*f2*r1*r3 + llT*f2*r1*r3 + lTl*f2*r1*r3 + lTT*f2*r1*r3 + Tll*f2*r1*r3 + TlT*f2*r1*r3 + TTl*f2*r1*r3 + TTT*f2*r1*r3))/((f1 - r1)*pow((f2 - r2),2)*(f3 - r3))),2) + pow(f1err,2)*pow(((f3*r2*(TlT*f2 - TTT + TTT*f2 + lTT*r1 + TTT*r1 + TTl*r3 + TTT*r3 - llT*f2*r1 - lTT*f2*r1 - TlT*f2*r1 - Tll*f2*r3 - TlT*f2*r3 - TTT*f2*r1 - TTl*f2*r3 - TTT*f2*r3 - lTl*r1*r3 - lTT*r1*r3 - TTl*r1*r3 - TTT*r1*r3 + lll*f2*r1*r3 + llT*f2*r1*r3 + lTl*f2*r1*r3 + lTT*f2*r1*r3 + Tll*f2*r1*r3 + TlT*f2*r1*r3 + TTl*f2*r1*r3 + TTT*f2*r1*r3))/((f1 - r1)*(f2 - r2)*(f3 - r3)) - (f1*f3*r2*(TlT*f2 - TTT + TTT*f2 + lTT*r1 + TTT*r1 + TTl*r3 + TTT*r3 - llT*f2*r1 - lTT*f2*r1 - TlT*f2*r1 - Tll*f2*r3 - TlT*f2*r3 - TTT*f2*r1 - TTl*f2*r3 - TTT*f2*r3 - lTl*r1*r3 - lTT*r1*r3 - TTl*r1*r3 - TTT*r1*r3 + lll*f2*r1*r3 + llT*f2*r1*r3 + lTl*f2*r1*r3 + lTT*f2*r1*r3 + Tll*f2*r1*r3 + TlT*f2*r1*r3 + TTl*f2*r1*r3 + TTT*f2*r1*r3))/(pow((f1 - r1),2)*(f2 - r2)*(f3 - r3))),2) + pow(f3err,2)*pow(((f1*r2*(TlT*f2 - TTT + TTT*f2 + lTT*r1 + TTT*r1 + TTl*r3 + TTT*r3 - llT*f2*r1 - lTT*f2*r1 - TlT*f2*r1 - Tll*f2*r3 - TlT*f2*r3 - TTT*f2*r1 - TTl*f2*r3 - TTT*f2*r3 - lTl*r1*r3 - lTT*r1*r3 - TTl*r1*r3 - TTT*r1*r3 + lll*f2*r1*r3 + llT*f2*r1*r3 + lTl*f2*r1*r3 + lTT*f2*r1*r3 + Tll*f2*r1*r3 + TlT*f2*r1*r3 + TTl*f2*r1*r3 + TTT*f2*r1*r3))/((f1 - r1)*(f2 - r2)*(f3 - r3)) - (f1*f3*r2*(TlT*f2 - TTT + TTT*f2 + lTT*r1 + TTT*r1 + TTl*r3 + TTT*r3 - llT*f2*r1 - lTT*f2*r1 - TlT*f2*r1 - Tll*f2*r3 - TlT*f2*r3 - TTT*f2*r1 - TTl*f2*r3 - TTT*f2*r3 - lTl*r1*r3 - lTT*r1*r3 - TTl*r1*r3 - TTT*r1*r3 + lll*f2*r1*r3 + llT*f2*r1*r3 + lTl*f2*r1*r3 + lTT*f2*r1*r3 + Tll*f2*r1*r3 + TlT*f2*r1*r3 + TTl*f2*r1*r3 + TTT*f2*r1*r3))/((f1 - r1)*(f2 - r2)*pow((f3 - r3),2))),2) + pow(r1err,2)*pow(((f1*f3*r2*(lTT + TTT - llT*f2 - lTT*f2 - TlT*f2 - TTT*f2 - lTl*r3 - lTT*r3 - TTl*r3 - TTT*r3 + lll*f2*r3 + llT*f2*r3 + lTl*f2*r3 + lTT*f2*r3 + Tll*f2*r3 + TlT*f2*r3 + TTl*f2*r3 + TTT*f2*r3))/((f1 - r1)*(f2 - r2)*(f3 - r3)) + (f1*f3*r2*(TlT*f2 - TTT + TTT*f2 + lTT*r1 + TTT*r1 + TTl*r3 + TTT*r3 - llT*f2*r1 - lTT*f2*r1 - TlT*f2*r1 - Tll*f2*r3 - TlT*f2*r3 - TTT*f2*r1 - TTl*f2*r3 - TTT*f2*r3 - lTl*r1*r3 - lTT*r1*r3 - TTl*r1*r3 - TTT*r1*r3 + lll*f2*r1*r3 + llT*f2*r1*r3 + lTl*f2*r1*r3 + lTT*f2*r1*r3 + Tll*f2*r1*r3 + TlT*f2*r1*r3 + TTl*f2*r1*r3 + TTT*f2*r1*r3))/(pow((f1 - r1),2)*(f2 - r2)*(f3 - r3))),2) + pow(r3err,2)*pow(((f1*f3*r2*(TTl + TTT - Tll*f2 - TlT*f2 - TTl*f2 - TTT*f2 - lTl*r1 - lTT*r1 - TTl*r1 - TTT*r1 + lll*f2*r1 + llT*f2*r1 + lTl*f2*r1 + lTT*f2*r1 + Tll*f2*r1 + TlT*f2*r1 + TTl*f2*r1 + TTT*f2*r1))/((f1 - r1)*(f2 - r2)*(f3 - r3)) + (f1*f3*r2*(TlT*f2 - TTT + TTT*f2 + lTT*r1 + TTT*r1 + TTl*r3 + TTT*r3 - llT*f2*r1 - lTT*f2*r1 - TlT*f2*r1 - Tll*f2*r3 - TlT*f2*r3 - TTT*f2*r1 - TTl*f2*r3 - TTT*f2*r3 - lTl*r1*r3 - lTT*r1*r3 - TTl*r1*r3 - TTT*r1*r3 + lll*f2*r1*r3 + llT*f2*r1*r3 + lTl*f2*r1*r3 + lTT*f2*r1*r3 + Tll*f2*r1*r3 + TlT*f2*r1*r3 + TTl*f2*r1*r3 + TTT*f2*r1*r3))/((f1 - r1)*(f2 - r2)*pow((f3 - r3),2))),2) + pow(f2err,2)*pow(((f1*f3*r2*(TlT + TTT - llT*r1 - lTT*r1 - TlT*r1 - Tll*r3 - TlT*r3 - TTT*r1 - TTl*r3 - TTT*r3 + lll*r1*r3 + llT*r1*r3 + lTl*r1*r3 + lTT*r1*r3 + Tll*r1*r3 + TlT*r1*r3 + TTl*r1*r3 + TTT*r1*r3))/((f1 - r1)*(f2 - r2)*(f3 - r3)) - (f1*f3*r2*(TlT*f2 - TTT + TTT*f2 + lTT*r1 + TTT*r1 + TTl*r3 + TTT*r3 - llT*f2*r1 - lTT*f2*r1 - TlT*f2*r1 - Tll*f2*r3 - TlT*f2*r3 - TTT*f2*r1 - TTl*f2*r3 - TTT*f2*r3 - lTl*r1*r3 - lTT*r1*r3 - TTl*r1*r3 - TTT*r1*r3 + lll*f2*r1*r3 + llT*f2*r1*r3 + lTl*f2*r1*r3 + lTT*f2*r1*r3 + Tll*f2*r1*r3 + TlT*f2*r1*r3 + TTl*f2*r1*r3 + TTT*f2*r1*r3))/((f1 - r1)*pow((f2 - r2),2)*(f3 - r3))),2) + (pow(TlTerr,2)*pow(f1,2)*pow(f3,2)*pow(r2,2)*pow((f2 - f2*r1 - f2*r3 + f2*r1*r3),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(lTTerr,2)*pow(f1,2)*pow(f3,2)*pow(r2,2)*pow((r1 - f2*r1 - r1*r3 + f2*r1*r3),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(TTlerr,2)*pow(f1,2)*pow(f3,2)*pow(r2,2)*pow((r3 - f2*r3 - r1*r3 + f2*r1*r3),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(TTTerr,2)*pow(f1,2)*pow(f3,2)*pow(r2,2)*pow((f2 + r1 + r3 - f2*r1 - f2*r3 - r1*r3 + f2*r1*r3 - 1),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(llTerr,2)*pow(f1,2)*pow(f3,2)*pow(r2,2)*pow((f2*r1 - f2*r1*r3),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(Tllerr,2)*pow(f1,2)*pow(f3,2)*pow(r2,2)*pow((f2*r3 - f2*r1*r3),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(lTlerr,2)*pow(f1,2)*pow(f3,2)*pow(r2,2)*pow((r1*r3 - f2*r1*r3),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(lllerr,2)*pow(f1,2)*pow(f2,2)*pow(f3,2)*pow(r1,2)*pow(r2,2)*pow(r3,2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2));
  return sqrt(res);
}
 
float MatrixMethod::N8_FFR_TTTerr(float r1,float f1,float r2,float f2,float r3,float f3, int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll, float r1err,float f1err,float r2err,float f2err,float r3err,float f3err){ 

float TTTerr = sqrt(TTT);
float llTerr = sqrt(llT);
float lTlerr = sqrt(lTl);
float Tllerr = sqrt(Tll);
float lTTerr = sqrt(lTT);
float TlTerr = sqrt(TlT);
float TTlerr = sqrt(TTl);
float lllerr = sqrt(lll);

  float res = pow(r3err,2)*pow(((f1*f2*(TTl*f3 - TTT + TTT*f3 + lTT*r1 + TlT*r2 + TTT*r1 + TTT*r2 - lTl*f3*r1 - lTT*f3*r1 - Tll*f3*r2 - TlT*f3*r2 - TTl*f3*r1 - TTT*f3*r1 - TTl*f3*r2 - TTT*f3*r2 - llT*r1*r2 - lTT*r1*r2 - TlT*r1*r2 - TTT*r1*r2 + lll*f3*r1*r2 + llT*f3*r1*r2 + lTl*f3*r1*r2 + lTT*f3*r1*r2 + Tll*f3*r1*r2 + TlT*f3*r1*r2 + TTl*f3*r1*r2 + TTT*f3*r1*r2))/((f1 - r1)*(f2 - r2)*(f3 - r3)) + (f1*f2*r3*(TTl*f3 - TTT + TTT*f3 + lTT*r1 + TlT*r2 + TTT*r1 + TTT*r2 - lTl*f3*r1 - lTT*f3*r1 - Tll*f3*r2 - TlT*f3*r2 - TTl*f3*r1 - TTT*f3*r1 - TTl*f3*r2 - TTT*f3*r2 - llT*r1*r2 - lTT*r1*r2 - TlT*r1*r2 - TTT*r1*r2 + lll*f3*r1*r2 + llT*f3*r1*r2 + lTl*f3*r1*r2 + lTT*f3*r1*r2 + Tll*f3*r1*r2 + TlT*f3*r1*r2 + TTl*f3*r1*r2 + TTT*f3*r1*r2))/((f1 - r1)*(f2 - r2)*pow((f3 - r3),2))),2) + pow(f1err,2)*pow(((f2*r3*(TTl*f3 - TTT + TTT*f3 + lTT*r1 + TlT*r2 + TTT*r1 + TTT*r2 - lTl*f3*r1 - lTT*f3*r1 - Tll*f3*r2 - TlT*f3*r2 - TTl*f3*r1 - TTT*f3*r1 - TTl*f3*r2 - TTT*f3*r2 - llT*r1*r2 - lTT*r1*r2 - TlT*r1*r2 - TTT*r1*r2 + lll*f3*r1*r2 + llT*f3*r1*r2 + lTl*f3*r1*r2 + lTT*f3*r1*r2 + Tll*f3*r1*r2 + TlT*f3*r1*r2 + TTl*f3*r1*r2 + TTT*f3*r1*r2))/((f1 - r1)*(f2 - r2)*(f3 - r3)) - (f1*f2*r3*(TTl*f3 - TTT + TTT*f3 + lTT*r1 + TlT*r2 + TTT*r1 + TTT*r2 - lTl*f3*r1 - lTT*f3*r1 - Tll*f3*r2 - TlT*f3*r2 - TTl*f3*r1 - TTT*f3*r1 - TTl*f3*r2 - TTT*f3*r2 - llT*r1*r2 - lTT*r1*r2 - TlT*r1*r2 - TTT*r1*r2 + lll*f3*r1*r2 + llT*f3*r1*r2 + lTl*f3*r1*r2 + lTT*f3*r1*r2 + Tll*f3*r1*r2 + TlT*f3*r1*r2 + TTl*f3*r1*r2 + TTT*f3*r1*r2))/(pow((f1 - r1),2)*(f2 - r2)*(f3 - r3))),2) + pow(f2err,2)*pow(((f1*r3*(TTl*f3 - TTT + TTT*f3 + lTT*r1 + TlT*r2 + TTT*r1 + TTT*r2 - lTl*f3*r1 - lTT*f3*r1 - Tll*f3*r2 - TlT*f3*r2 - TTl*f3*r1 - TTT*f3*r1 - TTl*f3*r2 - TTT*f3*r2 - llT*r1*r2 - lTT*r1*r2 - TlT*r1*r2 - TTT*r1*r2 + lll*f3*r1*r2 + llT*f3*r1*r2 + lTl*f3*r1*r2 + lTT*f3*r1*r2 + Tll*f3*r1*r2 + TlT*f3*r1*r2 + TTl*f3*r1*r2 + TTT*f3*r1*r2))/((f1 - r1)*(f2 - r2)*(f3 - r3)) - (f1*f2*r3*(TTl*f3 - TTT + TTT*f3 + lTT*r1 + TlT*r2 + TTT*r1 + TTT*r2 - lTl*f3*r1 - lTT*f3*r1 - Tll*f3*r2 - TlT*f3*r2 - TTl*f3*r1 - TTT*f3*r1 - TTl*f3*r2 - TTT*f3*r2 - llT*r1*r2 - lTT*r1*r2 - TlT*r1*r2 - TTT*r1*r2 + lll*f3*r1*r2 + llT*f3*r1*r2 + lTl*f3*r1*r2 + lTT*f3*r1*r2 + Tll*f3*r1*r2 + TlT*f3*r1*r2 + TTl*f3*r1*r2 + TTT*f3*r1*r2))/((f1 - r1)*pow((f2 - r2),2)*(f3 - r3))),2) + pow(r1err,2)*pow(((f1*f2*r3*(lTT + TTT - lTl*f3 - lTT*f3 - TTl*f3 - TTT*f3 - llT*r2 - lTT*r2 - TlT*r2 - TTT*r2 + lll*f3*r2 + llT*f3*r2 + lTl*f3*r2 + lTT*f3*r2 + Tll*f3*r2 + TlT*f3*r2 + TTl*f3*r2 + TTT*f3*r2))/((f1 - r1)*(f2 - r2)*(f3 - r3)) + (f1*f2*r3*(TTl*f3 - TTT + TTT*f3 + lTT*r1 + TlT*r2 + TTT*r1 + TTT*r2 - lTl*f3*r1 - lTT*f3*r1 - Tll*f3*r2 - TlT*f3*r2 - TTl*f3*r1 - TTT*f3*r1 - TTl*f3*r2 - TTT*f3*r2 - llT*r1*r2 - lTT*r1*r2 - TlT*r1*r2 - TTT*r1*r2 + lll*f3*r1*r2 + llT*f3*r1*r2 + lTl*f3*r1*r2 + lTT*f3*r1*r2 + Tll*f3*r1*r2 + TlT*f3*r1*r2 + TTl*f3*r1*r2 + TTT*f3*r1*r2))/(pow((f1 - r1),2)*(f2 - r2)*(f3 - r3))),2) + pow(r2err,2)*pow(((f1*f2*r3*(TlT + TTT - Tll*f3 - TlT*f3 - TTl*f3 - TTT*f3 - llT*r1 - lTT*r1 - TlT*r1 - TTT*r1 + lll*f3*r1 + llT*f3*r1 + lTl*f3*r1 + lTT*f3*r1 + Tll*f3*r1 + TlT*f3*r1 + TTl*f3*r1 + TTT*f3*r1))/((f1 - r1)*(f2 - r2)*(f3 - r3)) + (f1*f2*r3*(TTl*f3 - TTT + TTT*f3 + lTT*r1 + TlT*r2 + TTT*r1 + TTT*r2 - lTl*f3*r1 - lTT*f3*r1 - Tll*f3*r2 - TlT*f3*r2 - TTl*f3*r1 - TTT*f3*r1 - TTl*f3*r2 - TTT*f3*r2 - llT*r1*r2 - lTT*r1*r2 - TlT*r1*r2 - TTT*r1*r2 + lll*f3*r1*r2 + llT*f3*r1*r2 + lTl*f3*r1*r2 + lTT*f3*r1*r2 + Tll*f3*r1*r2 + TlT*f3*r1*r2 + TTl*f3*r1*r2 + TTT*f3*r1*r2))/((f1 - r1)*pow((f2 - r2),2)*(f3 - r3))),2) + pow(f3err,2)*pow(((f1*f2*r3*(TTl + TTT - lTl*r1 - lTT*r1 - Tll*r2 - TlT*r2 - TTl*r1 - TTT*r1 - TTl*r2 - TTT*r2 + lll*r1*r2 + llT*r1*r2 + lTl*r1*r2 + lTT*r1*r2 + Tll*r1*r2 + TlT*r1*r2 + TTl*r1*r2 + TTT*r1*r2))/((f1 - r1)*(f2 - r2)*(f3 - r3)) - (f1*f2*r3*(TTl*f3 - TTT + TTT*f3 + lTT*r1 + TlT*r2 + TTT*r1 + TTT*r2 - lTl*f3*r1 - lTT*f3*r1 - Tll*f3*r2 - TlT*f3*r2 - TTl*f3*r1 - TTT*f3*r1 - TTl*f3*r2 - TTT*f3*r2 - llT*r1*r2 - lTT*r1*r2 - TlT*r1*r2 - TTT*r1*r2 + lll*f3*r1*r2 + llT*f3*r1*r2 + lTl*f3*r1*r2 + lTT*f3*r1*r2 + Tll*f3*r1*r2 + TlT*f3*r1*r2 + TTl*f3*r1*r2 + TTT*f3*r1*r2))/((f1 - r1)*(f2 - r2)*pow((f3 - r3),2))),2) + (pow(TTlerr,2)*pow(f1,2)*pow(f2,2)*pow(r3,2)*pow((f3 - f3*r1 - f3*r2 + f3*r1*r2),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(lTTerr,2)*pow(f1,2)*pow(f2,2)*pow(r3,2)*pow((r1 - f3*r1 - r1*r2 + f3*r1*r2),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(TlTerr,2)*pow(f1,2)*pow(f2,2)*pow(r3,2)*pow((r2 - f3*r2 - r1*r2 + f3*r1*r2),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(TTTerr,2)*pow(f1,2)*pow(f2,2)*pow(r3,2)*pow((f3 + r1 + r2 - f3*r1 - f3*r2 - r1*r2 + f3*r1*r2 - 1),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(lTlerr,2)*pow(f1,2)*pow(f2,2)*pow(r3,2)*pow((f3*r1 - f3*r1*r2),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(Tllerr,2)*pow(f1,2)*pow(f2,2)*pow(r3,2)*pow((f3*r2 - f3*r1*r2),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(llTerr,2)*pow(f1,2)*pow(f2,2)*pow(r3,2)*pow((r1*r2 - f3*r1*r2),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(lllerr,2)*pow(f1,2)*pow(f2,2)*pow(f3,2)*pow(r1,2)*pow(r2,2)*pow(r3,2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2));
  return sqrt(res);
}
 
float MatrixMethod::N8_FFF_TTTerr(float r1,float f1,float r2,float f2,float r3,float f3, int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll, float r1err,float f1err,float r2err,float f2err,float r3err,float f3err){ 

float TTTerr = sqrt(TTT);
float llTerr = sqrt(llT);
float lTlerr = sqrt(lTl);
float Tllerr = sqrt(Tll);
float lTTerr = sqrt(lTT);
float TlTerr = sqrt(TlT);
float TTlerr = sqrt(TTl);
float lllerr = sqrt(lll);

  float res = pow(r1err,2)*pow(((f1*f2*f3*(lTT + TTT - llT*r2 - lTT*r2 - lTl*r3 - lTT*r3 - TlT*r2 - TTT*r2 - TTl*r3 - TTT*r3 + lll*r2*r3 + llT*r2*r3 + lTl*r2*r3 + lTT*r2*r3 + Tll*r2*r3 + TlT*r2*r3 + TTl*r2*r3 + TTT*r2*r3))/((f1 - r1)*(f2 - r2)*(f3 - r3)) + (f1*f2*f3*(lTT*r1 - TTT + TlT*r2 + TTT*r1 + TTT*r2 + TTl*r3 + TTT*r3 - llT*r1*r2 - lTT*r1*r2 - lTl*r1*r3 - lTT*r1*r3 - TlT*r1*r2 - Tll*r2*r3 - TlT*r2*r3 - TTT*r1*r2 - TTl*r1*r3 - TTT*r1*r3 - TTl*r2*r3 - TTT*r2*r3 + lll*r1*r2*r3 + llT*r1*r2*r3 + lTl*r1*r2*r3 + lTT*r1*r2*r3 + Tll*r1*r2*r3 + TlT*r1*r2*r3 + TTl*r1*r2*r3 + TTT*r1*r2*r3))/(pow((f1 - r1),2)*(f2 - r2)*(f3 - r3))),2) + pow(f1err,2)*pow(((f2*f3*(lTT*r1 - TTT + TlT*r2 + TTT*r1 + TTT*r2 + TTl*r3 + TTT*r3 - llT*r1*r2 - lTT*r1*r2 - lTl*r1*r3 - lTT*r1*r3 - TlT*r1*r2 - Tll*r2*r3 - TlT*r2*r3 - TTT*r1*r2 - TTl*r1*r3 - TTT*r1*r3 - TTl*r2*r3 - TTT*r2*r3 + lll*r1*r2*r3 + llT*r1*r2*r3 + lTl*r1*r2*r3 + lTT*r1*r2*r3 + Tll*r1*r2*r3 + TlT*r1*r2*r3 + TTl*r1*r2*r3 + TTT*r1*r2*r3))/((f1 - r1)*(f2 - r2)*(f3 - r3)) - (f1*f2*f3*(lTT*r1 - TTT + TlT*r2 + TTT*r1 + TTT*r2 + TTl*r3 + TTT*r3 - llT*r1*r2 - lTT*r1*r2 - lTl*r1*r3 - lTT*r1*r3 - TlT*r1*r2 - Tll*r2*r3 - TlT*r2*r3 - TTT*r1*r2 - TTl*r1*r3 - TTT*r1*r3 - TTl*r2*r3 - TTT*r2*r3 + lll*r1*r2*r3 + llT*r1*r2*r3 + lTl*r1*r2*r3 + lTT*r1*r2*r3 + Tll*r1*r2*r3 + TlT*r1*r2*r3 + TTl*r1*r2*r3 + TTT*r1*r2*r3))/(pow((f1 - r1),2)*(f2 - r2)*(f3 - r3))),2) + pow(f2err,2)*pow(((f1*f3*(lTT*r1 - TTT + TlT*r2 + TTT*r1 + TTT*r2 + TTl*r3 + TTT*r3 - llT*r1*r2 - lTT*r1*r2 - lTl*r1*r3 - lTT*r1*r3 - TlT*r1*r2 - Tll*r2*r3 - TlT*r2*r3 - TTT*r1*r2 - TTl*r1*r3 - TTT*r1*r3 - TTl*r2*r3 - TTT*r2*r3 + lll*r1*r2*r3 + llT*r1*r2*r3 + lTl*r1*r2*r3 + lTT*r1*r2*r3 + Tll*r1*r2*r3 + TlT*r1*r2*r3 + TTl*r1*r2*r3 + TTT*r1*r2*r3))/((f1 - r1)*(f2 - r2)*(f3 - r3)) - (f1*f2*f3*(lTT*r1 - TTT + TlT*r2 + TTT*r1 + TTT*r2 + TTl*r3 + TTT*r3 - llT*r1*r2 - lTT*r1*r2 - lTl*r1*r3 - lTT*r1*r3 - TlT*r1*r2 - Tll*r2*r3 - TlT*r2*r3 - TTT*r1*r2 - TTl*r1*r3 - TTT*r1*r3 - TTl*r2*r3 - TTT*r2*r3 + lll*r1*r2*r3 + llT*r1*r2*r3 + lTl*r1*r2*r3 + lTT*r1*r2*r3 + Tll*r1*r2*r3 + TlT*r1*r2*r3 + TTl*r1*r2*r3 + TTT*r1*r2*r3))/((f1 - r1)*pow((f2 - r2),2)*(f3 - r3))),2) + pow(f3err,2)*pow(((f1*f2*(lTT*r1 - TTT + TlT*r2 + TTT*r1 + TTT*r2 + TTl*r3 + TTT*r3 - llT*r1*r2 - lTT*r1*r2 - lTl*r1*r3 - lTT*r1*r3 - TlT*r1*r2 - Tll*r2*r3 - TlT*r2*r3 - TTT*r1*r2 - TTl*r1*r3 - TTT*r1*r3 - TTl*r2*r3 - TTT*r2*r3 + lll*r1*r2*r3 + llT*r1*r2*r3 + lTl*r1*r2*r3 + lTT*r1*r2*r3 + Tll*r1*r2*r3 + TlT*r1*r2*r3 + TTl*r1*r2*r3 + TTT*r1*r2*r3))/((f1 - r1)*(f2 - r2)*(f3 - r3)) - (f1*f2*f3*(lTT*r1 - TTT + TlT*r2 + TTT*r1 + TTT*r2 + TTl*r3 + TTT*r3 - llT*r1*r2 - lTT*r1*r2 - lTl*r1*r3 - lTT*r1*r3 - TlT*r1*r2 - Tll*r2*r3 - TlT*r2*r3 - TTT*r1*r2 - TTl*r1*r3 - TTT*r1*r3 - TTl*r2*r3 - TTT*r2*r3 + lll*r1*r2*r3 + llT*r1*r2*r3 + lTl*r1*r2*r3 + lTT*r1*r2*r3 + Tll*r1*r2*r3 + TlT*r1*r2*r3 + TTl*r1*r2*r3 + TTT*r1*r2*r3))/((f1 - r1)*(f2 - r2)*pow((f3 - r3),2))),2) + pow(r2err,2)*pow(((f1*f2*f3*(TlT + TTT - llT*r1 - lTT*r1 - TlT*r1 - Tll*r3 - TlT*r3 - TTT*r1 - TTl*r3 - TTT*r3 + lll*r1*r3 + llT*r1*r3 + lTl*r1*r3 + lTT*r1*r3 + Tll*r1*r3 + TlT*r1*r3 + TTl*r1*r3 + TTT*r1*r3))/((f1 - r1)*(f2 - r2)*(f3 - r3)) + (f1*f2*f3*(lTT*r1 - TTT + TlT*r2 + TTT*r1 + TTT*r2 + TTl*r3 + TTT*r3 - llT*r1*r2 - lTT*r1*r2 - lTl*r1*r3 - lTT*r1*r3 - TlT*r1*r2 - Tll*r2*r3 - TlT*r2*r3 - TTT*r1*r2 - TTl*r1*r3 - TTT*r1*r3 - TTl*r2*r3 - TTT*r2*r3 + lll*r1*r2*r3 + llT*r1*r2*r3 + lTl*r1*r2*r3 + lTT*r1*r2*r3 + Tll*r1*r2*r3 + TlT*r1*r2*r3 + TTl*r1*r2*r3 + TTT*r1*r2*r3))/((f1 - r1)*pow((f2 - r2),2)*(f3 - r3))),2) + pow(r3err,2)*pow(((f1*f2*f3*(TTl + TTT - lTl*r1 - lTT*r1 - Tll*r2 - TlT*r2 - TTl*r1 - TTT*r1 - TTl*r2 - TTT*r2 + lll*r1*r2 + llT*r1*r2 + lTl*r1*r2 + lTT*r1*r2 + Tll*r1*r2 + TlT*r1*r2 + TTl*r1*r2 + TTT*r1*r2))/((f1 - r1)*(f2 - r2)*(f3 - r3)) + (f1*f2*f3*(lTT*r1 - TTT + TlT*r2 + TTT*r1 + TTT*r2 + TTl*r3 + TTT*r3 - llT*r1*r2 - lTT*r1*r2 - lTl*r1*r3 - lTT*r1*r3 - TlT*r1*r2 - Tll*r2*r3 - TlT*r2*r3 - TTT*r1*r2 - TTl*r1*r3 - TTT*r1*r3 - TTl*r2*r3 - TTT*r2*r3 + lll*r1*r2*r3 + llT*r1*r2*r3 + lTl*r1*r2*r3 + lTT*r1*r2*r3 + Tll*r1*r2*r3 + TlT*r1*r2*r3 + TTl*r1*r2*r3 + TTT*r1*r2*r3))/((f1 - r1)*(f2 - r2)*pow((f3 - r3),2))),2) + (pow(Tllerr,2)*pow(f1,2)*pow(f2,2)*pow(f3,2)*pow((r2*r3 - r1*r2*r3),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(lTTerr,2)*pow(f1,2)*pow(f2,2)*pow(f3,2)*pow((r1 - r1*r2 - r1*r3 + r1*r2*r3),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(TlTerr,2)*pow(f1,2)*pow(f2,2)*pow(f3,2)*pow((r2 - r1*r2 - r2*r3 + r1*r2*r3),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(TTlerr,2)*pow(f1,2)*pow(f2,2)*pow(f3,2)*pow((r3 - r1*r3 - r2*r3 + r1*r2*r3),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(TTTerr,2)*pow(f1,2)*pow(f2,2)*pow(f3,2)*pow((r1 + r2 + r3 - r1*r2 - r1*r3 - r2*r3 + r1*r2*r3 - 1),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(llTerr,2)*pow(f1,2)*pow(f2,2)*pow(f3,2)*pow((r1*r2 - r1*r2*r3),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(lTlerr,2)*pow(f1,2)*pow(f2,2)*pow(f3,2)*pow((r1*r3 - r1*r2*r3),2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2)) + (pow(lllerr,2)*pow(f1,2)*pow(f2,2)*pow(f3,2)*pow(r1,2)*pow(r2,2)*pow(r3,2))/(pow((f1 - r1),2)*pow((f2 - r2),2)*pow((f3 - r3),2));
  return sqrt(res);
}


/* 
 * Function for doing some simple testing
 */
void MatrixMethod::TestMatrixInversion(){

  float estimate;

  float TTTerr=3.; 
  float TTlerr=2.; 
  float TlTerr=3.; 
  float lTTerr=2.; 
  float Tllerr=2.2; 
  float lTlerr=1.5; 
  float llTerr=2.1; 
  float lllerr=1.1;  
  
  float r1err=0.02; 
  float r2err=0.02; 
  float r3err=0.02; 
  float f1err=0.04; 
  float f2err=0.03; 
  float f3err=0.05;
    
  float r1=0.90; 
  float r2=0.93; 
  float r3=0.86; 
  float f1=0.25; 
  float f2=0.22; 
  float f3=0.30;

  int TT=10; 
  int Tl=8; 
  int lT=7; 
  int ll=5; 

  int TTT=10; 
  int TTl=8; 
  int TlT=7; 
  int lTT=6; 
  int Tll=5; 
  int lTl=7; 
  int llT=6; 
  int lll=3;

  printf("A quick look at the 2l-3D values:\n");
  estimate = N3_RR_TT(r1,f1,TT,Tl,ll)+N3_RF_TT(r1,f1,TT,Tl,ll)+N3_FF_TT(r1,f1,TT,Tl,ll);
  printf("MINIMAl TEST: for TT=%.2i we find RR+RF+FF = %.2f \n",TT,estimate);
  if(fabs(TT-estimate)<0.1)printf("  Ok, at least all 3 weights sum up to the original TT value. \n");
  else printf("  Not good, all 3 weights do NOT sum up to the original TT value. \n");

  printf("A quick look at the 2l-4D values:\n");
  estimate = N4_RR_TT(r1,f1,r2,f2,TT,Tl,lT,ll)+N4_RF_TT(r1,f1,r2,f2,TT,Tl,lT,ll)+N4_FR_TT(r1,f1,r2,f2,TT,Tl,lT,ll)+N4_FF_TT(r1,f1,r2,f2,TT,Tl,lT,ll);
  printf("MINIMAl TEST: for TT=%.2i we find RR+RF+FR+FF = %.2f \n",TT,estimate);
  if(fabs(TT-estimate)<0.1)printf("  Ok, at least all 4 weights sum up to the original TT value. \n");
  else printf("  Not good, all 4 weights do NOT sum up to the original TT value. \n");
  
  printf("A quick look at the 3l-8D values:\n");
  estimate = N8_RRR_TTT(r1,f1,r2,f2,r3,f3,TTT,llT,lTl,Tll,lTT,TlT,TTl,lll)+N8_RRF_TTT(r1,f1,r2,f2,r3,f3,TTT,llT,lTl,Tll,lTT,TlT,TTl,lll)+N8_RFR_TTT(r1,f1,r2,f2,r3,f3,TTT,llT,lTl,Tll,lTT,TlT,TTl,lll)+N8_FRR_TTT(r1,f1,r2,f2,r3,f3,TTT,llT,lTl,Tll,lTT,TlT,TTl,lll)+N8_RFF_TTT(r1,f1,r2,f2,r3,f3,TTT,llT,lTl,Tll,lTT,TlT,TTl,lll)+N8_FRF_TTT(r1,f1,r2,f2,r3,f3,TTT,llT,lTl,Tll,lTT,TlT,TTl,lll)+N8_FFR_TTT(r1,f1,r2,f2,r3,f3,TTT,llT,lTl,Tll,lTT,TlT,TTl,lll)+N8_FFF_TTT(r1,f1,r2,f2,r3,f3,TTT,llT,lTl,Tll,lTT,TlT,TTl,lll);
  printf("MINIMAL TEST: for TTT=%.2i we find RRR+RRF+..+FFF = %.2f \n",TTT,estimate);
  if(fabs(TTT-estimate)<0.1)printf("  Ok, at least all 8 weights sum up to the original TTT value. \n");
  else printf("  Not good, all 8 weights do NOT sum up to the original TT value. \n");

  printf("\nA quick look at the errors:\n");
  printf(" err RRR = %.2f\n",(N8_RRR_TTTerr( r1, f1, r2, f2, r3, f3,  TTT,  llT,  lTl,  Tll,  lTT,  TlT,  TTl,  lll,  r1err, f1err, r2err, f2err, r3err, f3err)));
  printf(" err RRF = %.2f\n",(N8_RRF_TTTerr( r1, f1, r2, f2, r3, f3,  TTT,  llT,  lTl,  Tll,  lTT,  TlT,  TTl,  lll,  r1err, f1err, r2err, f2err, r3err, f3err)));
  printf(" err RFF = %.2f\n",(N8_RFF_TTTerr( r1, f1, r2, f2, r3, f3,  TTT,  llT,  lTl,  Tll,  lTT,  TlT,  TTl,  lll,  r1err, f1err, r2err, f2err, r3err, f3err)));
  printf(" err FRF = %.2f\n",(N8_FRF_TTTerr( r1, f1, r2, f2, r3, f3,  TTT,  llT,  lTl,  Tll,  lTT,  TlT,  TTl,  lll,  r1err, f1err, r2err, f2err, r3err, f3err)));
  printf(" err FFF = %.2f\n",(N8_FFF_TTTerr( r1, f1, r2, f2, r3, f3,  TTT,  llT,  lTl,  Tll,  lTT,  TlT,  TTl,  lll,  r1err, f1err, r2err, f2err, r3err, f3err)));
  
}
