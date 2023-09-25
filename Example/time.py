#Example-time

import PyMsOfa  as  sf

#Day->Hour, minute, second
ndp=0
days=0.34
sign, ihmsf=sf.pymD2tf(ndp, days)
print(sign, ihmsf)
#b'+' (8, 9, 36, 0)

#Hour, minute, second->Radian
s = b'-'
ihour, imin, sec = 15, 25, 12.5
rad = sf.pymTf2a(s, ihour, imin, sec)
print(rad)
#-4.036982920888967

#Calendar->Julian day
iyear, imon, iday = 2023 , 9, 1 
dj1,dj2 = sf.pymCal2jd(iyear, imon, iday)
print(dj1,dj2)
#2400000.5 60188.0

#Julian epoch->Julian day
epj = 2023.1
djm0, djm = sf.pymEpj2jd(epj)
print(djm0, djm)
#2400000.5 59981.774999999965

#UTC->TAI
utc1 = 2400000.5
utc2 = 57388.5
tai1, tai2 = sf.pymUtctai(utc1, utc2)
print(tai1, tai2)
#2400000.5 57388.50041666667

#TAI->TT
tt1, tt2 = sf.pymTaitt(tai1, tai2)
print(tt1, tt2)
#2400000.5 57388.50078916667
