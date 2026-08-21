#Example-coordinate

import PyMsOfa  as  sf
import numpy  as  np

#eq->ec
date1, date2 = 2400000.5, 57388.5
dr, dd = 2.8, -0.6
dl, db = sf.pymEqec06(date1, date2, dr, dd)
print(dl, db)
#3.108110717008802 -0.679015560765235

#ec->eq
dr, dd = sf.pymEceq06(date1, date2, dl, db)
print(dr, dd)
#2.8 -0.6000000000000002

#Star proper motion(Fig.2 in paper)
ra1,dec1=0.08788796644747092,-1.132188403551008
pmr1,pmd1=8.27454225597756e-06,5.647882757242327e-06
px1=116.18257854098105e-3
rv1=9.283571
ep1a,ep1b=2400000.5,57388.450625000056
ep2a,ep2b=2400000.5,60991.538933193
ra2,dec2,pmr2,pmd2,px2,rv2\
    =sf.pymPmsafe(ra1, dec1, pmr1, pmd1, 
                       px1, rv1, ep1a, ep1b, ep2a, ep2b)
print(ra2,dec2,pmr2,pmd2,px2,rv2)
#0.08796958189264041 -1.1321326881050646 8.272396933840485e-06 
#5.648019447684745e-06 0.11618131477951303 9.287244274111199

#Projection theorem(Fig.3 in paper)
ra1,dec1=0.09025370018535234,-1.1308656514899067
pmr1,pmd1=-1.4191556529435842e-08,4.093790696930718e-08
px1=1.6330602766416882e-3
rv1=-2.694851
ep1a,ep1b=2400000.5,57388.450625000056
ep2a,ep2b=2400000.5,60991.538933193

ra21,dec21,pmr2,pmd2,px2,rv2\
    =sf.pymPmsafe(ra1, dec1, pmr1, pmd1,
                  px1, rv1, ep1a, ep1b, ep2a, ep2b)
a,b=ra21,dec21
a0,b0=ra2,dec2
xi,eta=sf.pymTpxes(a, b, a0, b0)
print(xi,eta)
#0.0009726944659440729 0.001266436096867058

#ICRS->GCRS
date1, date2 = 2400000.5, 60188.0
pv   = np.array([[-6241497.16 ,  401346.896,  -1251136.04],    
                 [-29.264597,  -455.021831,  0.0266151194]])  
astrom = sf.pymApcs13(date1, date2, pv)
rc,dc=0.08788796644747092,-1.132188403551008
pr,pd=8.27454225597756e-06,5.647882757242327e-06
px=116.18257854098105e-3
rv=9.283571
ri, di = sf.pymAtciq(rc, dc, pr, pd, px, rv, astrom)
print(ri, di)
#0.0882705549440727 -1.1320010534992886

#mean longitude of Mars
t=2460188.5
lamda = sf.pymFama03(t)
print(lamda)
#1.9096191364710933