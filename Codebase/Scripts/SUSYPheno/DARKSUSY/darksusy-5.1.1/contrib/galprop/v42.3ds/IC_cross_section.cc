
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * IC_cross_section.cc *                         galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include<iostream.h>
#include<math.h>
//                       Compton-Scatter Cross section in cm^2 MeV-1
//         
//           Egamma in    MeV
//           Eelectron in MeV
//           Etarget in    eV           
//           See GFY thesis page 70
//           C++ version from fortran dsde1.f90
 
 

double IC_cross_section(double Egamma,double Eelectron,double Etarget){


float WQ;
float ro,mc2,Pi,x,G,k,F,q;
int Code;

ro = 2.817e-13;
mc2= 0.511     ;                                         
Pi = 3.14159    ;                                             
                                                
WQ = 0.0 ;                                                 

x = Egamma / Eelectron  ;                                     
G = Eelectron / mc2     ;                                    
k = 4.0 * Etarget * 1.0e-6 * G / mc2   ;
Code = 0       ;                 

if ( x > 0.999999 )           Code = 1;
         

F = 0.0;
if (Code == 0) {                                                         
       q = 1.0/k * x / (1.0-x)    ;                                   
       F=2.0*q*log (q)+(1.0-q)*(1.0+2.0*q)+0.5*pow(q*k, 2)/(1.0+k*q)*(1.0-q);
}        

if (F <= 0.0) {                                               
                 F=0.0    ;                                 
                 Code = 2 ;
		   }                                                             

if (Code == 0)  WQ = 2.0*Pi*ro*ro /(G*G) / Etarget * 1.0e6 * F ;   

// cout<<" IC_cross_section Egamma Eelectron Etarget cross section "<<Egamma <<" " <<Eelectron <<" " <<Etarget <<"     " <<WQ <<endl;

return WQ;
}
