
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * IC_anisotropy_factor.cc *                     galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include<iostream.h>
#include"galprop_classes.h"
#include"galprop.h"
#include"constants.h"

double IC_anisotropy_factor (double rho,double xi, double z, int isrf_comp){
  //cout<<">>>>IC_anisotropy_factor"<<endl;
  //cout<<"IC_anisotropy_factor:isrf_comp="<<isrf_comp<<" rho="<<rho<<" xi="<<xi<<" z="<<z<<endl;


  if(isrf_comp==2){return 1.0;}// this is the third component which is microwave background

extern int isrf_energy_density_i_comp;


 double nu       =1.e13;// typical target photon frequency
 double Eelectron=1000.;// typical electron energy in MeV

int kbg=1; //transparent disk
double E0=h_planck*nu*erg_to_eV/1.e6/m_electron;// target photon energy in mc^2 units
double E=Eelectron/m_electron;                  // electron energy in mc^2 units
double PLindex=3;                              // power law index
double RG=15.;                                 // radius of Galacti disc

double RS=8.5;                                 // distance of Sun from Galactic centre
 
 
 isrf_energy_density_i_comp=isrf_comp; // as required by cfactor_cc


 double factor=cfactor_cc(kbg,E0,E,PLindex,RG,rho,xi,z,RS);

 //cout<<"isrf_comp="<<isrf_comp<<" rho="<<rho<<" xi="<<xi<<" z="<<z<<" cfactor="<<factor<<endl;

 //cout<<"<<<<IC_anisotropy_factor"<<endl;
return factor;
}
//////////////////////////////////////////////
/*
      real*8 function CFactor(kbg,E0,E,PLindex,RG,rho,xi,z,RS)
c***********************************************************************
c                            *** I.Moskalenko, version of 31.03.1998 ***
c calculation of the correction factor for
c INVERSE COMPTON scattering in an ANISOTROPIC photon field
c    INPUT:
c kbg =-1 for isotropic background photon field (constant); 
c     = 0 for opaque disk photon field (~cos(theta));
c     = 1 for transparent disk photon field (~1/cos(theta)); 
c E0 is the energy of a background photon (in electron rest mass units mc^2);
c E  is the energy of scattered photon (in electron rest mass units mc^2);
c PLindex is the power-law index of the electron spectrum (>0);
c RG is the Galactic radius (kpc);
c (rho,xi,z) are the cylindrical coordinates of the electron position
c   in respect to the Galactic center (rho,z in kpc, 0< xi (radians) <2Pi 
c   is polar angle between the electron position and observer direction);
c RS is the Solar distance from the Galactic center (kpc).
c    internal OUTPUT:
c SPEC - the differential energy spectrum of gamma-rays (phot/ccm/sec/energy)
c    per one electron in ccm divided by Pi*r0^2*c;
c DENS - number density of background photons at the electron position
c    (the emissivity of unit area of the Galactic plane is equal to 1).
c    REFERENCES:
c I.V.Moskalenko & A.W.Strong 1999, ApJ submitted
c***********************************************************************
*/
