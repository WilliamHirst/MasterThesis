'''
PDG = {}
PDG['TtoL']    = 0.1736 + 0.1784
PDG['ZtoLL']   = 0.0363 + 0.0366
PDG['ZtoTT']   = 0.0337
PDG['Ztoqq']   = 1 - ZtoLL - ZtoTT
PDG['WtoL']    = 0.1075 + 0.1057
PDG['WtoT']    = 0.1125
PDG['WtoLtot'] = WtoL + WtoT*TtoL
PDG['Wtoq']    = 1 - WtoL - WtoT
PDG['Wtoqtot'] = Wtoq + WtoT*(1-TtoL)
'''

from math import sqrt, log
import array


TtoL  = 0.1736 + 0.1784
ZtoLL = 0.0363 + 0.0366
ZtoTT = 0.0337
ZtoLLtot = 0.0363 + 0.0366 + ZtoTT * TtoL * TtoL
Ztoqq = 1 - ZtoLL - ZtoTT
WtoL = 0.1075 + 0.1057
WtoT = 0.1125
WtoLtot = WtoL + WtoT*TtoL
Wtoq = 1 - WtoL - WtoT
Wtoqq = 1 - WtoL - WtoT
Wtoqtot = Wtoq + WtoT*(1-TtoL)



# http://en.wikipedia.org/wiki/Lambda-CDM_model: 
cdm = 0.1123
cdm_error = 0.0035


def Z_LLR(S,B):
    return sqrt( 2*( (S+B)*log(1+1.*S/B) - S) )

    


