# Reference:
#     Wallace P T, Capitaine N. Precession-nutation procedures consistent with 
#     IAU 2006 resolutions[J]. Astronomy & Astrophysics, 2006, 459(3): 981-985.

import numpy  as  np
import PyMsOfa  as  sf

utc1, utc2 = sf.pymCal2jd(2006, 1, 15)
D = sf.pymTf2d('+', 21, 24, 37.5)
utc2 +=D

dut1=0.3341

date1, date2 = sf.pymUtcut1(utc1, utc2, dut1)

print(date1, date2)
# =============================================================================
# 2400000.5 53750.892104561346
# =============================================================================

dpsi, deps, eosa, rb, rp, rbp, rn, rbpn = sf.pymPn06a(date1, date2)

print(rbpn)
# =============================================================================
# [[ 9.99998923e-01 -1.34606929e-03 -5.84803118e-04]
#  [ 1.34604476e-03  9.99999093e-01 -4.23222297e-05]
#  [ 5.84859557e-04  4.15350129e-05  9.99999828e-01]]
# =============================================================================

x, y, s = sf.pymXys06a(date1, date2)

CIO = sf.pymC2ixys(x, y, s)

print(CIO)
# =============================================================================
# [[ 9.99999829e-01  3.23202243e-10 -5.84859557e-04]
#  [-2.46153515e-08  9.99999999e-01 -4.15350056e-05]
#  [ 5.84859557e-04  4.15350129e-05  9.99999828e-01]]
# =============================================================================

ERA = sf.pymEra00(date1, date2)

R = sf.pymRz(ERA, CIO)

sign, ERA = sf.pymA2tf(5, sf.pymEra00(date1, date2))

print(sign, ERA)
# =============================================================================
# b'+' (5, 5, 3, 70345)
# =============================================================================

print(R)
# =============================================================================
# [[ 2.37424215e-01  9.71406048e-01 -1.79207215e-04]
#  [-9.71405888e-01  2.37424279e-01  5.58274693e-04]
#  [ 5.84859557e-04  4.15350129e-05  9.99999828e-01]]
# =============================================================================