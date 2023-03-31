#include <string.h>
#include <math.h>
#include <stdio.h>

#include"../CalcHEP_src/c_source/ntools/include/vegas.h"
#include"micromegas.h"

#define q 1.6e-19
#define c 3.e8
#define ep0 8.85e-12
#define mu0 1.25e-6
#define meGeV 0.511e-3
#define GeV2Kg 1.78e-27
#define me meGeV*GeV2Kg

// input parameters 
static double Mhp, chi, sigmav;
static double mX;
static double freq;
static double theta;
static double B01; // normalisation of 2nd Bmag component
static double *SpE_stat;
static double alpha;
// micrOMEGAs parameters
//                    diffusion model
#define delta Delta_dif
#define L     L_dif
#define K0    K_dif
#define Rg    Rdisk

//                   DM profile
#define Rho0  rhoDM
#define rho0  Rsun

//step function definition
#define StepF(x) ((x)>0? 1:0)

//magnetic field profile
static double Bmag(double rho, double z) {
  /*
this is the function defining the magnetic field, cylindrical coord.
  */
#define B0 6. //microGauss
#define r2nd 2. // radius of 2nd component
 
  double Bm = 0., aux = 0., Bcore = 0., aux00 = 0;
 
  aux = exp(-(rho-rho0)/(delta*Rg))*exp(-fabs(z)/(delta*L));

  aux00 = exp(rho0/(delta*Rg));
  Bcore = B01/aux00 - B0; 
  Bm = B0*aux + StepF(r2nd-sqrt(rho*rho + z*z))*Bcore*aux;

//  Bm = B0*aux + StepF(r2nd-sqrt(rho*rho+z*z))*B01*aux;      

  return Bm;
}

// DM density profile
static double hProfileQ(double rho, double z)  
{ 
  double r=sqrt( rho*rho + z*z), prQ;
  if(r<1.E-10) r=1.E-10; 
  prQ=hProfile_(r); prQ*=prQ+rhoClumpEff_(r)/rhoDM; 
  return prQ;
}

static double spectrum(double En, double mx) { if(En>mx) return 0; else return  SpectdNdE(En,SpE_stat); }

// synchrotron energy loss rate
static double Psync(double En, double Bm) {

  double Psynch = 0.;
  
  Psynch = pow(q,4)*pow(En*Bm,2)/(6.*M_PI*ep0*pow(me,4));
  Psynch = Psynch*pow(10,-20)*GeV2Kg*pow(c,-3);
  Psynch = 2./3.*Psynch; // pitch angle average

  return Psynch;
}

// ICS energy loss rate
static double Pics(double En, double Bm) {
#define Urad 8. // eV/cm3
  double ics = 0.;
  ics = 2./3.*Urad*pow(Bm,-2)*2.*mu0;
  ics = ics*pow(10,17)*GeV2Kg*pow(c,2)*Psync(En,Bm);

  return ics;
  }

// total energy loss rate
static double Eloss(double En, double Bm) 
{
  double Et = Psync(En,Bm)+Pics(En,Bm);
  return Et;
}

// halo function, Green formalism
static double Halo(double lambda, double Ld, 
            double rho, double phi, double z, 
            double rhoS, double phiS, double zS)
{
  double sum = 0., expo = 0.;
  int n;
if(!lambda) return 0;  
  for (n=-10; n<11; n++)
    {
      expo = exp( - pow(2.*n*L + pow(-1,n)*zS - z,2)*pow(lambda,-2) );
      sum = sum + pow(-1,n)*expo;
    }
  double expoA = pow(rho,2)+pow(rhoS,2)-2.*rho*rhoS*cos(phi-phiS);
         expoA = exp(-expoA*pow(lambda,-2));
  double fac = 1/lambda/sqrt(M_PI);
  double Green = pow(fac,2)*expoA*fac*sum;
  double Halo = Green* hProfileQ(rhoS,zS);  
  return Halo;
}

// electron energy
// 1-to-1 corresp. with synch. frequency
static double Ene(double nu, double Bm) {

  double Ep = 0.;
  Ep = 4.*M_PI*pow(me,3)*nu/3./q/Bm * pow(10,19)*pow(GeV2Kg,-2);
  Ep = sqrt(Ep);
  Ep = Ep/sqrt(0.29);

  return Ep;
}

// diffusion length
static double lambda(double en, double es, double Bm) {

 double tt = 0.;
 tt = pow(Eloss(en,Bm)*pow(en,-2),-1);
 tt = -3.17*pow(10,-14)*tt*pow(en,delta-1) / (delta -1); 
 double lamb = 0.;
 lamb = tt;
 tt = pow(Eloss(es,Bm)*pow(es,-2),-1);
 tt = -3.17*pow(10,-14)*tt*pow(es,delta-1) / (delta -1); 

if(4.*K0*(lamb - tt)<0){ printf("es=%E en=%E\n", es,en);  return 0;}

 lamb = sqrt(4.*K0*(lamb - tt));

 return lamb;
}

// integrand of flux
static double primflux(double los, double en, 
                double es, double rhoS, double phiS, double zS) {
  //  is integrand of flux per unit solid angle

  double pF = 0.;
  double sv = 0.;

if(en>es) {printf("en=%E es=%E\n",en,es); exit(0);}
#define Dc 8.5 
#define sv sigmav
#define svth (3.*pow(10,-26))
#define prefact 1.21*pow(10,8)*0.5*0.5

  double rho = 0.; double z = 0.;
 rho = pow(los*sin(theta)*sin(alpha),2) + pow(Rsun - los*sin(theta)*cos(alpha),2);
 rho = sqrt(rho);
 z = los*cos(theta);



 double B = Bmag(rho,z);
 double Ee = Ene(freq,B);
 
 pF = Halo(lambda(Ee,es,B),L,rho,0.,z,rhoS,phiS,zS);
 pF = spectrum(es,mX)*pF; 
 pF = Psync(Ee,B)/Eloss(Ee,B)*pF/sqrt(B);
 
 pF = sv/svth*pow(Rho0,2)*pow(mX,-2)*pF;
 pF = pF/sqrt(freq); 
 pF = pF/sqrt(0.29);//  from dEe/dnu_peak
 pF = prefact*1./4./M_PI*pF;
 pF = 2.*pF;// from electrons and positrons 

 
 return pF;
}

// maximum line-of-sight
// [ consistency with mx > Ecrit = sqrt(nu/B[los])  ]
static double losMax(double mx, double nu, double thet)
{

  double losDEF = 0.;

  double logaux = log(0.1708709*pow(mx,2)/nu);
  double cT = cos(thet), sT = sin(thet);
  double funcA = 0.;

  funcA = pow(68.*sT-560.*cT*logaux,2);
  funcA = funcA - 4.*(784.*pow(logaux,2)-289.)*(100.*pow(cT,2)-4.*pow(sT,2));
  funcA = sqrt(funcA);
  funcA = 0.5*(560.*cT*logaux - funcA - 68.*sT)/(100.*pow(cT,2)-4.*pow(sT,2));
  
  double funcB = 0.;
  funcB = L/cT;

  if(funcA < funcB) { losDEF = funcA;}
  else {losDEF = funcB;}

  return losDEF;
}

// integrand to VEGAS
                        
static int Integrand(const int *ndim, const double xx[],
                     const int *ncomp, double ff[], void *userdata)
{

  double func = 0.;
  
#define sphit xx[0] /* azimuth integration */
#define szt   xx[1] /* z integration       */
#define srhot xx[2] /* distance to source  */ 
#define sEt   xx[3] /* source energy       */
#define lofst xx[4] /* line of sight integration */
#define f     ff[0]

#define srhoMin 1.e-50
#define srhoMax 20.
#define losMin 1.e-50
#define losmax (losMax(mX,freq,theta))

//  #define lofs ( losMin*exp(lofst*log(losmax/losMin)) )
//  #define Jaclofs ( losMin*log(losmax/losMin)*exp(lofst*log(losmax/losMin)) )

#define lofs (lofst*losmax)
#define Jaclofs (losmax)
 


  double R= lofs*sin(theta);
  double Z= lofs*cos(theta);
  R=sqrt( Rsun*Rsun + R*R -2*cos(alpha)*R*Rsun);
  
#define sEMin ( Ene(freq,Bmag(R,Z)) )
#define sEMax mX
  #define sE ( sEMin*exp(sEt*log(sEMax/sEMin)) )
  #define JacsE ( sEMin*log(sEMax/sEMin)*exp(sEt*log(sEMax/sEMin)) )

//  #define sphi ( sphiMin*exp(sphit*log(sphiMax/sphiMin)) )
//  #define Jacsphi ( sphiMin*log(sphiMax/sphiMin)*exp(sphit*log(sphiMax/sphiMin)) )
#define sphi (M_PI*sphit)
#define Jacsphi 2*M_PI


#define sz (2.*L*szt - L )
#define Jacsz (2.*L)

//  #define srho ( srhoMin*exp(srhot*log(srhoMax/srhoMin)) )
//  #define Jacsrho ( srho*srhoMin*log(srhoMax/srhoMin)*exp(srhot*log(srhoMax/srhoMin)) )

  #define srho  sqrt(srhot)*srhoMax
  #define Jacsrho (srhoMax*srhoMax/2) 


{ double B = Bmag(R,Z);
  double ecrit = Ene(freq,B);
//printf(" ecrit=%E %E\n", ecrit,   Ene(freq,Bmag(Rsun-lofs*sin(theta)  ,Z)));  
  double l=lambda(ecrit, sE,B);
  double xmin=exp(-srhoMax*srhoMax/l/l);
  
  double xr= xmin+srhot*(1-xmin);
  double srho_= l*sqrt(-log(xr));
  double Jacsrho_=l*l/xr/2.*(1-xmin);
 
   f = primflux(lofs,ecrit,sE,srho,sphi,sz)*Jaclofs*JacsE*Jacsrho*Jacsphi*Jacsz;
}
  return 0;

}

static double vInt(double*x,double wgt)
{ int ndim=5,ncomp=1;
  double f;
  Integrand(&ndim,x,&ncomp , &f,NULL);
  return f;
}  


double  Tsynch(double freq_,double theta_, double alpha_, double B01_,double sigmav_,double*SpE)
{
  double result;
  
  SpE_stat=SpE;
  
  freq  = freq_  ;
  theta = M_PI/2-theta_*M_PI/180 ;
  B01   = B01_   ;
  mX    = Mcdm   ;
  sigmav= sigmav_;
  alpha=alpha_*M_PI/180;

{ vegasGrid *vegPtr=NULL;
  double ti,dti;
  
  vegPtr=vegas_init(5,50);
  
  vegas_int(vegPtr, 10000 , 1.5, vInt, &ti, &dti);
  printf("ti=%E dti=%E\n",ti, dti);
  vegas_int(vegPtr, 100000 , 1.5, vInt, &ti, &dti);
  printf("ti=%E dti=%E\n",ti, dti);
  vegas_int(vegPtr, 100000 , 1.5, vInt, &ti, &dti);
  printf("ti=%E dti=%E\n",ti, dti);
  vegas_int(vegPtr, 100000 , 1.5, vInt, &ti, &dti);
  printf("ti=%E dti=%E\n",ti, dti);
  vegas_int(vegPtr, 100000 , 1.5, vInt, &ti, &dti);
  printf("ti=%E dti=%E\n",ti, dti);
  vegas_int(vegPtr, 100000 , 1.5, vInt, &ti, &dti);
  printf("ti=%E dti=%E\n",ti, dti);
  vegas_int(vegPtr, 100000 , 1.5, vInt, &ti, &dti);
  printf("ti=%E dti=%E\n",ti, dti);
  vegas_int(vegPtr, 100000 , 1.5, vInt, &ti, &dti);
  printf("ti=%E dti=%E\n",ti, dti);
  result=ti;
  vegas_finish(vegPtr);
}  

#define h 4.1e-15
#define ccm 3.e10
#define kB 8.6e-5
#define conv (1.e23 * 0.16e-11)
  double Temp = 0.; double fluxT = 0.;
#define nuHz (freq*1.e9)
#define nuMHz (freq*1.e3)

  Temp = 1. + 2.*conv*h*pow(nuHz,3)*pow(ccm,-2)/result;
  Temp = h*nuHz/kB/log(Temp);
  fluxT = Temp*pow(nuMHz,2.5);

  return fluxT;
}
