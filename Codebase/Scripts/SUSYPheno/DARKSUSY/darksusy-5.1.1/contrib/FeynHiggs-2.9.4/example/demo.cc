// demo.cc
// demonstration program for calling FeynHiggs from C++
// this file is part of FeynHiggs
// last modified 13 May 11 th


#include <stdlib.h>
#include <stdio.h>
#include "CFeynHiggs.h"
#include "CSLHA.h"


/**********************************************************************/

static void setFlags()
{
  enum {
    mssmpart = 4,
    fieldren = 0,
    tanbren = 0,
    higgsmix = 2,
    p2approx = 0,
    looplevel = 2,
    runningMT = 1,
    botResum = 1,
    tlCplxApprox = 3
  };

  int error;

  FHSetFlags(&error, mssmpart, fieldren, tanbren, higgsmix,
    p2approx, looplevel, runningMT, botResum, tlCplxApprox);
  if( error ) exit(error);
}

/**********************************************************************/

static void setPara()
{
  const double invAlfa = -1;
  const double AlfasMZ = -1;
  const double GF = -1;

  const double ME = -1;
  const double MU = -1;
  const double MD = -1;
  const double MM = -1;
  const double MC = -1;
  const double MS = -1;
  const double ML = -1;
  const double MB = -1;
  const double MW = -1;
  const double MZ = -1;

  const double CKMlambda = -1;
  const double CKMA = -1;
  const double CKMrhobar = -1;
  const double CKMetabar = -1;

  const double MT = 172;
  const double TB = 5;
  const double MA0 = 250;
  const double MHp = -1;

  const double MSusy = 1000;
  const double M3SL = MSusy, M2SL = M3SL, M1SL = M2SL;
  const double M3SE = MSusy, M2SE = M3SE, M1SE = M2SE;
  const double M3SQ = MSusy, M2SQ = M3SQ, M1SQ = M2SQ;
  const double M3SU = MSusy, M2SU = M3SU, M1SU = M2SU;
  const double M3SD = MSusy, M2SD = M3SD, M1SD = M2SD;

  const double_complex Af = 2000;
  const double_complex At = Af, Ac = At, Au = Ac;
  const double_complex Ab = Af, As = Ab, Ad = As;
  const double_complex Atau = Af, Amu = Atau, Ae = Amu;

  const double_complex MUE = 200;
  const double_complex M_1 = 0;
  const double_complex M_2 = 500;
  const double_complex M_3 = 800;

  const double Qtau = 0;
  const double Qt = 0;
  const double Qb = 0;

  const double scalefactor = 1;

  int error;

  FHSetSMPara(&error,
    invAlfa, AlfasMZ, GF,
    ME, MU, MD, MM, MC, MS, ML, MB,
    MW, MZ,
    CKMlambda, CKMA, CKMrhobar, CKMetabar);
  if( error ) exit(error);

  FHSetPara(&error, scalefactor,
    MT, TB, MA0, MHp,
    M3SL, M3SE, M3SQ, M3SU, M3SD,
    M2SL, M2SE, M2SQ, M2SU, M2SD,
    M1SL, M1SE, M1SQ, M1SU, M1SD,
    MUE,
    Atau, At, Ab,
    Amu, Ac, As,
    Ae, Au, Ad,
    M_1, M_2, M_3,
    Qtau, Qt, Qb);
  if( error ) exit(error);
}

/**********************************************************************/

static void setSLHA(const char *filename)
{
  int error;
  double_complex slhadata[nslhadata];

  SLHARead(&error, slhadata, filename, 1);
  if( error ) exit(error);

  FHSetSLHA(&error, slhadata);
  if( error ) exit(error);
}

/**********************************************************************/

static void getPara()
{
  int error, nmfv;
  double MASf[4][6], MCha[2], MNeu[4];
  double_complex UASf[4][6][6];
  double_complex UCha[2][2], VCha[2][2], ZNeu[4][4];
  double_complex DeltaMB;
  double MGl;
  double MHtree[4], SAtree;

  FHGetPara(&error, &nmfv, MASf, UASf,
    MCha, UCha, VCha, MNeu, ZNeu, &DeltaMB, &MGl,
    MHtree, &SAtree);
  if( error ) exit(error);

// print some sample output:
  printf("MCha = %g %g\n"
         "MNeu = %g %g %g %g\n", 
    MCha[0], MCha[1],
    MNeu[0], MNeu[1], MNeu[2], MNeu[3]);
}

/**********************************************************************/

static void higgsCorr()
{
  int error;
  double MHiggs[4];
  double_complex SAeff, UHiggs[3][3], ZHiggs[3][3];

  FHHiggsCorr(&error, MHiggs, &SAeff, UHiggs, ZHiggs);
  if( error ) exit(error);

// print some sample output:
  printf("Mh1 = %g\n"
         "Mh2 = %g\n"
         "Mh3 = %g\n"
         "MHp = %g\n"
         "sin alpha_eff = %g %g\n",
    MHiggs[0], MHiggs[1], MHiggs[2], MHiggs[3],
    Re(SAeff), Im(SAeff));
}

/**********************************************************************/

int main()
{
  setFlags();

// either use setPara to set the parameters directly
// or use setSLHA to read them from an SLHA file
  setPara();
//  setSLHA("myfile.SLHA");

  getPara();

  higgsCorr();
}

