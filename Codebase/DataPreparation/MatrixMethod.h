#ifndef MatrixMethod_h
#define MatrixMethod_h


class MatrixMethod
{

 public:

  ~MatrixMethod() {}
  MatrixMethod() {}

  void TestMatrixInversion();

  // 4D inversion functions
  // LL estimate
  float N4_RR_LL(float r1,float f1,float r2,float f2, 
		 float NTT, float NTl, float NlT, float Nll);
  float N4_RF_LL(float r1,float f1,float r2,float f2, 
		 float NTT, float NTl, float NlT, float Nll);
  float N4_FR_LL(float r1,float f1,float r2,float f2, 
		 float NTT, float NTl, float NlT, float Nll);
  float N4_FF_LL(float r1,float f1,float r2,float f2, 
		 float NTT, float NTl, float NlT, float Nll);
  //TT estimate
  float N4_RR_TT(float r1,float f1,float r2,float f2, 
		 float NTT, float NTl, float NlT, float Nll);
  float N4_RF_TT(float r1,float f1,float r2,float f2, 
		 float NTT, float NTl, float NlT, float Nll);
  float N4_FR_TT(float r1,float f1,float r2,float f2, 
		 float NTT, float NTl, float NlT, float Nll);
  float N4_FF_TT(float r1,float f1,float r2,float f2, 
		 float NTT, float NTl, float NlT, float Nll);
  // errors
  float N4_RF_TTerr(float r1,float f1,float r2,float f2,  
		    float nTT,float nTl,float nlT,float nll,  
		    float r1err,float f1err,float r2err,float f2err);
  float N4_RR_TTerr(float r1,float f1,float r2,float f2,  
		    float nTT,float nTl,float nlT,float nll,  
		    float r1err,float f1err,float r2err,float f2err);
  float N4_FR_TTerr(float r1,float f1,float r2,float f2,  
		    float nTT,float nTl,float nlT,float nll,  
		    float r1err,float f1err,float r2err,float f2err);
  float N4_FF_TTerr(float r1,float f1,float r2,float f2,  
		    float nTT,float nTl,float nlT,float nll,  
		    float r1err,float f1err,float r2err,float f2err);


  // 3D inversion functions
  // LL estimate
  float N3_RR_LL(float r, float f, float NTT,float NTl,float Nll);
  float N3_RF_LL(float r, float f, float NTT,float NTl,float Nll);
  float N3_FF_LL(float r, float f, float NTT,float NTl,float Nll);
  //TT estimate
  float N3_RR_TT(float r, float f, float NTT,float NTl,float Nll);
  float N3_RF_TT(float r, float f, float NTT,float NTl,float Nll);
  float N3_FF_TT(float r, float f, float NTT,float NTl,float Nll);
  // errors
  float N3_RR_TTerr(float r, float f, float NTT,float NTl,float Nll, float rerr, float ferr);
  float N3_RF_TTerr(float r, float f, float NTT,float NTl,float Nll, float rerr, float ferr);
  float N3_FF_TTerr(float r, float f, float NTT,float NTl,float Nll, float rerr, float ferr);


  // 8D inversion functions
  // LLL estimate
  float  N8_RRR_LLL(float r1,float f1,float r2,float f2,float r3,float f3, 
		    int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll);
  float  N8_RRF_LLL(float r1,float f1,float r2,float f2,float r3,float f3, 
		    int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll);
  float  N8_RFR_LLL(float r1,float f1,float r2,float f2,float r3,float f3, 
		    int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll);
  float  N8_FRR_LLL(float r1,float f1,float r2,float f2,float r3,float f3, 
		    int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll);
  float  N8_RFF_LLL(float r1,float f1,float r2,float f2,float r3,float f3, 
		    int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll);
  float  N8_FRF_LLL(float r1,float f1,float r2,float f2,float r3,float f3, 
		    int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll);
  float  N8_FFR_LLL(float r1,float f1,float r2,float f2,float r3,float f3, 
		    int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll);
  float  N8_FFF_LLL(float r1,float f1,float r2,float f2,float r3,float f3, 
		    int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll);
  // TTT estimate
  float  N8_RRR_TTT(float r1,float f1,float r2,float f2,float r3,float f3, 
		    int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll);
  float  N8_RRF_TTT(float r1,float f1,float r2,float f2,float r3,float f3, 
		    int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll);
  float  N8_RFR_TTT(float r1,float f1,float r2,float f2,float r3,float f3, 
		    int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll);
  float  N8_FRR_TTT(float r1,float f1,float r2,float f2,float r3,float f3, 
		    int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll);
  float  N8_RFF_TTT(float r1,float f1,float r2,float f2,float r3,float f3, 
		    int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll);
  float  N8_FRF_TTT(float r1,float f1,float r2,float f2,float r3,float f3, 
		    int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll);
  float  N8_FFR_TTT(float r1,float f1,float r2,float f2,float r3,float f3, 
		    int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll);
  float  N8_FFF_TTT(float r1,float f1,float r2,float f2,float r3,float f3, 
		    int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll);
  // errors
  float N8_FFF_TTTerr(float r1,float f1,float r2,float f2,float r3,float f3, 
		      int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll,
		      float r1err,float f1err,float r2err,float f2err,float r3err,float f3err); 
  float N8_FFR_TTTerr(float r1,float f1,float r2,float f2,float r3,float f3, 
		      int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll,
		      float r1err,float f1err,float r2err,float f2err,float r3err,float f3err); 
  float N8_FRF_TTTerr(float r1,float f1,float r2,float f2,float r3,float f3, 
		      int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll,
		      float r1err,float f1err,float r2err,float f2err,float r3err,float f3err); 
  float N8_RFF_TTTerr(float r1,float f1,float r2,float f2,float r3,float f3, 
		      int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll,
		      float r1err,float f1err,float r2err,float f2err,float r3err,float f3err); 
  float N8_FRR_TTTerr(float r1,float f1,float r2,float f2,float r3,float f3, 
		      int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll,
		      float r1err,float f1err,float r2err,float f2err,float r3err,float f3err); 
  float N8_RFR_TTTerr(float r1,float f1,float r2,float f2,float r3,float f3, 
		      int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll,
		      float r1err,float f1err,float r2err,float f2err,float r3err,float f3err); 
  float N8_RRF_TTTerr(float r1,float f1,float r2,float f2,float r3,float f3, 
		      int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll,
		      float r1err,float f1err,float r2err,float f2err,float r3err,float f3err); 
  float N8_RRR_TTTerr(float r1,float f1,float r2,float f2,float r3,float f3, 
		      int TTT, int llT, int lTl, int Tll, int lTT, int TlT, int TTl, int lll,
		      float r1err,float f1err,float r2err,float f2err,float r3err,float f3err); 


};


#endif
