
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * global.cc *                                   galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

// global structures to be accessed via extern

#include "galprop_classes.h"

double global_var;


Spectrum S;
Distribution Dist;



Particle P;

int n_species;
int isrf_energy_density_i_comp; // required for cfactor

///////////////////////////////////
Configure configure;
Galdef galdef;

Particle *gcr; // all species
Galaxy    galaxy; 






