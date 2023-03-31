This is the dark matter (DM) version of GALPROP code (imos 9/14/2005)
---------------------------------------------------------------------
The list of modified files:

create_galaxy.cc
create_gcr.cc
Galaxy.h
Galdef.cc
Galdef.h
galprop.cc
galprop.h
gen_bremss_emiss.cc
gen_DM_skymap.cc
gen_DM_source.cc
gen_IC_emiss.cc
gen_secondary_source.cc
gen_synch_emiss.cc
nuclei_normalize.cc
nuc_package.cc
store_DM_emiss.cc
store_DM_skymap.cc
pp_meson.f

Restrictions:
The 3D version is not fully consistent. The calculation of the
gamma-ray emission and skymap from the DM annihilation is calculated
assuming 2D cylindrical symmetry.
The integration over the line of sight is restricted by R=r_max,
and z_min<=Z<=z_max (see galdef-file for the definitions).

---------------------------------------------------------------------