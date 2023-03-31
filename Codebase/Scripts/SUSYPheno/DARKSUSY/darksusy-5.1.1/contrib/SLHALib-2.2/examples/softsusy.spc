# SUSY Les Houches Accord - MSSM spectrum
# SPS point 1a
# Compute with SOFTSUSY1.8
# B.C. Allanach, Comput. Phys. Commun. 143 (2002) 305-331, hep-ph/0104145
Block MODSEL                 # Select model
    1    1                   # sugra
Block MINPAR                 # Input parameters
    1   1.000000000e+02      # m0
    2   2.500000000e+02      # m12
    3   1.000000000e+01      # tanb
    4   1.000000000e+00      # sign(mu)
    5  -1.000000000e+02      # A0
Block SMINPUTS               # Standard Model inputs
    1   1.279340000e+02      # alpha^(-1) SM MSbar(MZ)
    2   1.166370000e-05      # G_Fermi
    3   1.172000000e-01      # alpha_s(MZ) SM MSbar
    4   9.118760000e+01      # MZ(pole)
    5   4.250000000e+00      # mb(mb) SM MSbar
    6   1.743000000e+02      # mtop(pole)
    7   1.777000000e+00      # mtau(pole)
#
Block SPINFO                 # Program information
    1    SOFTSUSY            # spectrum calculator
    2    1.8                 # version number
#
Block MASS  # Mass spectrum
#PDG code      mass              particle
          6    1.743000000e+02  # top
         25    1.096471686e+02  # h0
         35    3.905646065e+02  # H0
         36    3.849267509e+02  # A0
         37    3.963987424e+02  # H+
    1000001    5.537379281e+02  # ~d_L
    1000002    5.480648005e+02  # ~u_L
    1000003    5.536689385e+02  # ~s_L
    1000004    5.479950083e+02  # ~c_L
    1000005    4.990864878e+02  # ~b_1
    1000006    3.866681125e+02  # ~t_1
    1000011    2.005077001e+02  # ~e_L
    1000012    1.844822029e+02  # ~snue_L
    1000013    2.005050044e+02  # ~mu_L
    1000014    1.844792730e+02  # ~snumu_L
    1000015    1.339969762e+02  # ~stau_1
    1000016    1.836242253e+02  # ~snu_tau_L
    1000021    5.934756712e+02  # ~g
    1000022    9.701573617e+01  # ~neutralino(1)
    1000023    1.788864799e+02  # ~neutralino(2)
    1000024    1.782649096e+02  # ~chargino(1)
    1000025   -3.536102287e+02  # ~neutralino(3)
    1000035    3.733417082e+02  # ~neutralino(4)
    1000037    3.736128390e+02  # ~chargino(2)
    2000001    5.269676664e+02  # ~d_R
    2000002    5.311251030e+02  # ~u_R
    2000003    5.269652151e+02  # ~s_R
    2000004    5.309795680e+02  # ~c_R
    2000005    5.257115262e+02  # ~b_2
    2000006    5.704560875e+02  # ~t_2
    2000011    1.430886701e+02  # ~e_R
    2000013    1.430810123e+02  # ~mu_R
    2000015    2.043832731e+02  # ~stau_2
# Higgs mixing
Block alpha # Higgs mixing parameters
    -1.146864127e-01  # alpha
Block hmix Q= 4.520624648e+02  # Higgs mixing parameters
    1    3.439934743e+02  # mu
Block stopmix  # stop mixing matrix
  1  1    5.443784304e-01  # O_{11}
  1  2    8.388397490e-01  # O_{12}
  2  1    8.388397490e-01  # O_{21}
  2  2   -5.443784304e-01  # O_{22}
Block sbotmix  # sbottom mixing matrix
  1  1    9.355024721e-01  # O_{11}
  1  2    3.533201449e-01  # O_{12}
  2  1   -3.533201449e-01  # O_{21}
  2  2    9.355024721e-01  # O_{22}
Block staumix  # stau mixing matrix
  1  1    2.810947184e-01  # O_{11}
  1  2    9.596800297e-01  # O_{12}
  2  1    9.596800297e-01  # O_{21}
  2  2   -2.810947184e-01  # O_{22}
Block nmix  # neutralino mixing matrix
  1  1    9.849417415e-01  # N_{1,1}
  1  2   -5.795970738e-02  # N_{1,2}
  1  3    1.526931274e-01  # N_{1,3}
  1  4   -5.670314904e-02  # N_{1,4}
  2  1    1.090115410e-01  # N_{2,1}
  2  2    9.374300545e-01  # N_{2,2}
  2  3   -2.852021039e-01  # N_{2,3}
  2  4    1.673354023e-01  # N_{2,4}
  3  1   -6.143190096e-02  # N_{3,1}
  3  2    9.173963120e-02  # N_{3,2}
  3  3    6.949466769e-01  # N_{3,3}
  3  4    7.105343608e-01  # N_{3,4}
  4  1   -1.192995029e-01  # N_{4,1}
  4  2    3.308313851e-01  # N_{4,2}
  4  3    6.421788575e-01  # N_{4,3}
  4  4   -6.811200615e-01  # N_{4,4}
Block Umix  # chargino U mixing matrix 
  1  1    9.084497528e-01  # U_{1,1}
  1  2   -4.179940750e-01  # U_{1,2}
  2  1    4.179940750e-01  # U_{2,1}
  2  2    9.084497528e-01  # U_{2,2}
Block Vmix  # chargino V mixing matrix 
  1  1    9.691507874e-01  # V_{1,1}
  1  2   -2.464685606e-01  # V_{1,2}
  2  1    2.464685606e-01  # V_{2,1}
  2  2    9.691507874e-01  # V_{2,2}
Block gauge Q= 4.520624648e+02  # (SUSY scale)
    1    3.611138488e-01  # g'(Q)MSSM DRbar
    2    6.459589284e-01  # g(Q)MSSM DRbar
    3    1.076546223e+00  # g3(Q)MSSM DRbar
Block yu Q= 4.520624648e+02  # (SUSY scale)
  3  3   8.944749998e-01  # Yt(Q)MSSM DRbar
Block yd Q= 4.520624648e+02  # (SUSY scale)
  3  3   1.454916965e-01  # Yb(Q)MSSM DRbar
Block ye Q= 4.520624648e+02  # (SUSY scale)
  3  3   9.977348866e-02  # Ytau(Q)MSSM DRbar
Block au Q= 4.520624648e+02  # (SUSY scale)
  3  3  -4.713090005e+02  # At(Q)MSSM DRbar
Block ad Q= 4.520624648e+02  # (SUSY scale)
  3  3  -6.039636586e+02  # Ab(Q)MSSM DRbar
Block ae Q= 4.520624648e+02  # (SUSY scale)
  3  3  -2.513543943e+02  # Atau(Q)MSSM DRbar

