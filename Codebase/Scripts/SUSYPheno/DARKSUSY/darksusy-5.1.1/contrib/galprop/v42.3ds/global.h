
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * global.h *                                    galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#ifndef global_h
#define global_h

// global structures to be accessed via extern

extern double global_var;
extern Spectrum S;
extern Distribution Dist;
extern Particle P;
extern int n_species;
extern int isrf_energy_density_i_comp; // required for cfactor

/////////////////////////////////////////////
extern Configure configure;
extern Galdef   galdef;
extern Particle *gcr;
extern Galaxy   galaxy;

#endif



