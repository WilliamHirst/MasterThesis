### THIS IS THE NEW STANDARD FOR EXPERIMENTAL VALUES   # 2014-01-26

import math

exp = {}

# ### PARTICLE DATA GROUP
exp['Ztoee'] = 3.363e-2
exp['Ztomm'] = 3.366e-2   # ee is smaller by 0.09% = 1e-4 ; probably within errors : So Z asymmetry is absent
exp['ZtoTT'] = 3.370e-2
exp['Ztovv'] = 20.00e-2
#exp['ZtoQQ = 69.91e-2
exp['ZtoQQ'] = 1 - (exp['Ztoee'] + exp['Ztomm'] + exp['ZtoTT'] + exp['Ztovv'])   #=69.901e-2  #Sep20
exp['Ztobb'] = 15.12e-2

exp['Ztoll'] = exp['Ztoee'] + exp['Ztomm']
exp['Ztoqq'] = exp['ZtoQQ'] - exp['Ztobb']

exp['Wtoev'] = 10.75e-2
exp['Wtomv'] = 10.57e-2  # ev is 1.7% larger than mv   So W asymmetry is noticeable .. 
exp['WtoTv'] = 11.25e-2
#exp['WtoQQ'] = 67.60e-2
#exp['Wtoqq'] = 67.60e-2  # the b-contribution should be negligible
exp['Wtoqq'] = 67.43e-2  # the b-contribution should be negligible
# Sep20:had to manually take Wtoqq down to acheive Sum(W) = 1.0


exp['Wtolv'] = exp['Wtoev'] + exp['Wtomv']

exp['Ttoevv'] = 17.85e-2 #+ 1.75e-2  #2011-01-10: small bug discovered. The small numbers are included in the large (but I added them)
exp['Ttomvv'] = 17.36e-2 #+ 0.36e-2     # evv is 2.8% larger than mvv ; noticeable asymmetry
exp['Ttolvv'] = exp['Ttoevv'] + exp['Ttomvv']
exp['Ttovh']  = 1. - exp['Ttolvv']

exp['ttobW'] = 1.


# ##########

exp['amu'] = 0.00116592089  # amu'] = (gmu-2)/2  http://en.wikipedia.org/wiki/Anomalous_magnetic_dipole_moment
exp['amuErrStat'] = 5.4e-10
exp['amuErrSyst'] = 3.3e-10
exp['amuErr'] = math.sqrt(exp['amuErrStat']**2+exp['amuErrSyst']**2)  # well..  6.3e-10 .. this is what I must compare the output with

#exp['bsgam'] = 3.60e-4     # 2001, http://www.slac.stanford.edu/cgi-wrap/getdoc/slac-r-670.pdf, p.4 ; I've seen 3.55e-4 too
#exp['bsgamErr'] = 0.30e-4  

exp['bsgam'] = 3.55e-4  # http://arxiv.org/pdf/1207.2520v2.pdf
exp['bsgamErrStat'] = 0.24e-4
exp['bsgamErrSyst'] = 0.09e-4
exp['bsgamErr'] = math.sqrt(exp['bsgamErrStat']**2 + exp['bsgamErrSyst']**2)  # 0.265e-4

exp['rho'] = 1.0004  # http://pdg.lbl.gov/2013/reviews/rpp2013-rev-standard-model.pdf [p.37] (http://ilcphys.kek.jp/meeting/physics/archives/2013-05-18/Yokoya_HWW_130518.pdf)
exp['rhoErrPlus']  = 0.0003
exp['rhoErrMinus'] = 0.0004

#exp['rho'] = 1.0008  
#exp['rhoErrPlus '] = 0.0017
#exp['rhoErrMinus'] = 0.0007


# http://en.wikipedia.org/wiki/Lambda-CDM_model 
#exp['cdm'] = 0.1123         # Omega_c h^2  (Omega_c = 0.227 +- 0.014)
#exp['cdm_error'] = 0.0035

# http://pdg.lbl.gov/2013/reviews/rpp2013-rev-dark-matter.pdf
exp['cdm'] = 0.1198
exp['cdm_error'] = 0.0026

