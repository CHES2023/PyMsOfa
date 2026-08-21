# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 02:44:08 2023
Done    on Mon Jun  5 02:11:37 2023
@author: Dr. Jianghui JI  (jijh@pmo.ac.cn)
"""
#jit  can speed up python
#from  numba  import jit
import ctypes as  ct
from   ctypes import *
import numpy  as  np
import numpy.ctypeslib  as  nt
import PyMsOfa  as  sf


def   line():
        print('-'*60)

def   line_c():
        print('/'*60)   


#dll = ctypes.cdll.LoadLibrary
#lib = CDLL('./libsofa_c.so')

line()
print("Testing pymPas...")
      
al  = c_double(np.pi/3.0)
ap  = c_double(np.pi/6.0)
bl  = c_double(np.pi/2.0)
bp  = c_double(np.pi/6.0)
ret_val  = sf.pymPas(al, ap, bl, bp)
print(ret_val)

line()
print("Testing pymS2c...")


theta = np.pi/4.0
phi   = np.pi/4.0

c = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C'))
sf.pymS2c_A(theta, phi, c)
print(c)

c = sf.pymS2c(theta, phi)
print(c)
 

print("Testing pymCal2jd...")
iyear = 2020
imon  = 12
iday  = 2
print(iyear, imon, iday)
 

djm0 = c_double()
djm  = c_double()   
sf.pymCal2jd_A(iyear, imon, iday, djm0, djm) 
print(djm0, djm)

jm0,jm = sf.pymCal2jd(iyear, imon, iday) 
print(jm0, jm)

print("Testing A2af angle in radians to degr, arcmin, arcsec, fraction...")
ndp = 3
#ndp = -2
angle = np.pi/2 - 0.01
sign = c_char()
sign, v = sf.pymA2af(ndp, angle)
print(sign,v)
line()

print("Testing pymJd2Cal...")
dj1 = 2400000.5E0
dj2 =   59185.0E0
iy, im, iday, fday = sf.pymJd2cal(dj1, dj2)
print(iy, im, iday, fday)
line()
print("Testing pymJdCalf...")
t1 = sf.pymJdcalf(ndp, dj1, dj2)
print(t1)  
 

print("Testing pymD2tf...")
ndp = 6
days = 0.961940171285E0
sign, v = sf.pymD2tf(ndp, days)
print(sign,v)
line()

scale = c_wchar_p('UTC')
ndp   = 6
dj1 = 2400000.5E0
dj2 =   59185.961940171285E0
#iy, im, iday, v = sf.pymD2dtf(scale, ndp, dj1, dj2) 
v = sf.pymD2dtf(scale, ndp, dj1, dj2)
print(v) 

print("Testing pymDAT...")
iy, im, iday, fday = 2020, 12, 2, 0.961940171285E0
delat = sf.pymDat(iy, im, iday, fday)
 
print(iy, im, iday, delat)
line()

date1, date2, ut, elong, u, v = 2400000.5E0, 59185.961, 0.0, 30.0, 6000.0, 5000.0
delt  = sf.pymDtdb(date1, date2, ut, elong, u, v)
print(delt)
line()

print("Testing pymDtf2d...")
iy, im, id, ihr, imn, sec = 2020, 12, 2, 23, 5, 11.6 
d1,d2 = sf.pymDtf2d(scale,  iy,  im,  id, ihr, imn, sec) 
scale = c_wchar_p('UTC')
print(scale.value, d1, d2)
line()

print("Testing pymEpb...")
epb = sf.pymEpb(d1, d2)
epc = sf.pymEpb_A(d1, d2)
print(epb, epc)
line()

print("Testing pymEpb2jd...")
d1, d2 = sf.pymEpb2jd(epb)
print(d1,d2)
d1, d2 = sf.pymEpb2jd_A(epb)
print(d1,d2)
line()
print("Testing pymEpj...")
epb = sf.pymEpj(d1, d2)
epc = sf.pymEpj_A(d1, d2)
print(epb, epc)
line()


print("Testing pymEpj2jd...")
d1, d2 = sf.pymEpj2jd(epb)
print(d1,d2)
d1, d2 = sf.pymEpj2jd_A(epc)
print(d1,d2)
line()

print("Testing pymTaitt...")
tai1, tai2 = 2400000.5E0, 59185.961E0
d1, d2 = sf.pymTaitt(tai1, tai2)
print(d1,d2)
line()

print("Testing pymTaiut1...")
tai1, tai2, dta = 2459185.5E0, 0.961E0, 0.5E0
ut11, ut12= sf.pymTaiut1(tai1, tai2, dta)
print(ut11,ut12)
line()


print("Testing pymTaiutc...")
tai1, tai2 = 2459185.5E0, 0.961E0
utc1, utc2 = sf.pymTaiutc(tai1, tai2)
print(utc1,utc2)

line()
print("Testing pymTcbtdb...")
tcb1, tcb2 = 2459185.5E0, 0.961E0
tdb1, tdb2 = sf.pymTcbtdb(tcb1, tcb2)
print(tdb1, tdb2)


line()
print("Testing pymTcgtt...")
tcg1, tcg2 = 2459185.5E0, 0.961E0
tt1,  tt2 = sf.pymTcgtt(tcg1, tcg2)
print(tt1, tt2)

line()
print("Testing pymTdbtcb...")
tdb1, tdb2 = 2459185.5E0, 0.960751265471194E0
tcb1, tcb2 = sf.pymTdbtcb(tdb1, tdb2)
print(tcb1, tcb2)


line()
print("Testing pymTdbtt...")
tdb1, tdb2, dtr = 2459185.5E0, 0.960751265471194E0, 0.01E0
tt1,  tt2 = sf.pymTdbtt(tdb1, tdb2, dtr)
print(tt1, tt2)

line()
print("Testing pymTf2d...")
s  = c_wchar('-')
ihour, imin, sec =  23, 5, 11.6 
days = sf.pymTf2d(s, ihour, imin, sec) 
print(days)
 


line()
print("Testing pymTttai...")
tt1, tt2=2400000.5E0, 59185.9613725E0
tai1, tai2 = sf.pymTttai(tt1, tt2) 
print(tai1, tai2)  
 
line()
print("Testing pymTcgtt...")
tt1,  tt2  = 2459185.5E0, 0.9609888198922069E0
tcg1, tcg2 = sf.pymTttcg(tt1,  tt2)
print(tcg1, tcg2)

line()
print("Testing pymTttdb...")
tt1,  tt2, dtr  = 2459185.5E0, 0.9607511497304533, 0.01E0
tdb1, tdb2 = sf.pymTttdb(tt1, tt2, dtr) 
print(tdb1, tdb2)

line()
print("Testing pymTtut1...")
tt1,  tt2, dtr  = 2459185.5E0, 0.9607511497304533, 0.01E0
ut11, ut12= sf.pymTtut1(tt1, tt2, dtr)
print(ut11, ut12)

line()
print("Testing pymUt1tai...")
ut11, ut12, dta  = 2459185.5E0, 0.961005787037037, 0.5E0
tai1, tai2 = sf.pymUt1tai(ut11, ut12, dta)
print(tai1, tai2) 

line()
print("Testing pymUt1tt...")
ut11, ut12, dt  = 2459185.5E0, 0.9607510339897125, 0.01E0
tt1,  tt2 = sf.pymUt1tt(ut11, ut12, dt)
print(tt1,  tt2)


line()
print("Testing pymUt1utc...")
ut11, ut12, dut1  = 2459185.5E0, 0.9607510339897125, 0.01E0
utc1, utc2  = sf.pymUt1utc(ut11, ut12, dut1)
print(utc1, utc2)
 
line()
print("Testing Utctai...")
utc1, utc2  = 2459185.5E0, 0.9605717592592592
tai1, tai2  = sf.pymUtctai(utc1, utc2)
print(tai1, tai2)


line()
print("Testing Utcut1...")
utc1, utc2, dut1  = 2459185.5E0, 0.9607509182489717, 0.01E0
ut11, ut12  = sf.pymUtcut1(utc1, utc2, dut1)
print(ut11, ut12)

line()
print("Testing Utc2TT...")
utc1, utc2 = 2455347E0, 0.090322
tai1, tai2  = sf.pymUtctai(utc1, utc2)
print('UTC:', utc1, utc2, utc1 + utc2)
print('TAI:', tai1, tai2, tai1 + tai2)
tt1, tt2 = sf.pymTaitt(tai1, tai2)
print('TT: ',  tt1,  tt2, tt1 + tt2)



line()
print("Testing Ab...")
c = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C'))

pnat = np.array([-0.76321968546737951, -0.60869453983060384, -0.21676408580639883]).reshape(1,3)
v    = np.array([2.1044018893653786e-5, -8.9108923304429319e-5, -3.8633714797716569e-5]).reshape(1,3)
s, bm1 = 0.99980921395708788,0.99999999506209258
c = sf.pymAb(pnat, v, s, bm1)
print(c)

line()
print("Testing Ae2hd...")
az, el, phi = 5.5, 1.1, 0.7
ha, dec = sf.pymAe2hd(az, el, phi)
print(ha,dec)

line()
print("Testing Apcg...")
date1, date2 = 2456165.5E0, 0.901310875
ebpv = np.array([[0.901310875, -0.417402664, -0.180982288],
                 [0.00742727954, 0.0140507459, 0.00609045792]]) 
ehp  = np.array([0.903358544, -0.415395237, -0.180084014]).reshape(1,3)
#astrom = sf.c_star_p()
#print(ehp)
astrom = sf.pymApcg(date1, date2, ebpv, ehp) 
#print('pmt', astrom.pmt)
#print('bm1', astrom.bm1)
sf.pymASTROM.print(astrom)

line()
print("Testing Apcg13...")
 
date1, date2 = 2456165.5E0, 0.401182685 
astrom = sf.pymApcg13(date1, date2) 
sf.pymASTROM.print(astrom)
#print('bm1', astrom.bm1)
sf.pymASTROM.print(astrom)

line()
print("Testing Apci...")
 
date1, date2 = 2456165.5E0, 0.401182685 
x, y, s =  0.0013122272,  -2.92808623e-5, 3.05749468e-8
astrom = sf.pymApci(date1, date2, ebpv, ehp, x, y, s) 

sf.pymASTROM.print(astrom) 

line()
print("Testing Apci13...")
date1, date2 = 2456165.5E0, 0.401182685

astrom, eo = sf.pymApci13(date1, date2)
sf.pymASTROM.print(astrom) 
print(eo)

line()
print("Testing Apco...")
date1, date2 = 2456384.5, 0.970031644
ebpv = np.array([[-0.974170438,  -0.211520082, -0.0917583024],
                 [0.00364365824, -0.0154287319, -0.00668922024]]) 
ehp  = np.array([-0.973458265, -0.209215307, -0.0906996477]).reshape(1,3)

x, y, s, theta = 0.0013122272, -2.92808623e-5,  3.05749468e-8,  3.14540971               
elong, phi, hm = -0.527800806, -1.2345856, 2738.0      
xp, yp, sp = 2.47230737e-7, 1.82640464e-6, -3.01974337e-11         
refa,  refb = 0.000201418779, -2.36140831e-7    
astrom=sf.pymApco(date1, date2, ebpv, ehp, x, y, s, theta, elong, phi, hm,  
                  xp, yp, sp, refa, refb)
sf.pymASTROM.print(astrom) 

line()
print("Testing Apco13...")
utc1, utc2, dut1  = 2456384.5, 0.969254051,  0.1550675     
elong, phi, hm = -0.527800806, -1.2345856, 2738.0             
xp, yp  = 2.47230737e-7, 1.82640464e-6    
phpa, tc = 731.0, 12.8          
rh, wl = 0.59, 0.55           

astrom, eo = sf.pymApco13(utc1, utc2, dut1, elong, phi, hm, xp, yp, 
                          phpa, tc, rh,  wl)
sf.pymASTROM.print(astrom) 
print(eo)


line()
print("Testing Apci...")
date1, date2 = 2456384.5, 0.970031644                     
pv   = np.array([[-1836024.09,   1056607.72,   -5998795.26],           
                 [-77.0361767,  -133.310856,  0.0971855934]])   

ebpv = np.array([[-0.974170438,  -0.211520082, -0.0917583024],           
                 [0.00364365824, -0.0154287319, -0.00668922024]])        
ehp  = np.array([-0.973458265, -0.209215307, -0.0906996477]).reshape(1,3)

astrom = sf.pymApcs(date1, date2, pv, ebpv, ehp)
sf.pymASTROM.print(astrom) 

line()
print("Testing Apci13...")
date1, date2 = 2456165.5, 0.401182685
pv   = np.array([[-6241497.16 ,  401346.896,  -1251136.04],    
                 [-29.264597,  -455.021831,  0.0266151194]])   

astrom = sf.pymApcs13(date1, date2, pv)
sf.pymASTROM.print(astrom) 

line()
print("Testing Aper...")
astrom  = sf.pymASTROM()
astrom.along = 1.234E0
#print(astrom.along)
theta = 5.678E0 
sf.pymAper(theta, astrom)
sf.pymASTROM.print(astrom) 
#print(astrom.eral)

line()
print("Testing Aper13...")
astrom  = sf.pymASTROM()
astrom.along = 1.234E0
ut11, ut12 = 2456165.5, 0.401182685
 
sf.pymAper13(ut11, ut12, astrom)
sf.pymASTROM.print(astrom) 


line()
print("Testing Apio...")
astrom  = sf.pymASTROM()

sp = -3.01974337e-11 
theta = 3.14540971   
elong = -0.527800806 
phi = -1.2345856     
hm = 2738.0          
xp = 2.47230737e-7   
yp = 1.82640464e-6   
refa = 0.000201418779
refb = -2.36140831e-7
 
astrom = sf.pymApio(sp, theta, elong, phi,  hm,  xp,  yp,  refa,  refb)
sf.pymASTROM.print(astrom) 


line()
print("Testing Apio13...")
utc1, utc2, dut1  = 2456384.5, 0.969254051,  0.1550675     
elong, phi, hm = -0.527800806, -1.2345856, 2738.0             
xp, yp  = 2.47230737e-7, 1.82640464e-6    
phpa, tc = 731.0, 12.8          
rh, wl = 0.59, 0.55           

astrom = sf.pymApio13(utc1, utc2, dut1, elong, phi, hm, xp, yp, 
                          phpa, tc, rh,  wl)
sf.pymASTROM.print(astrom) 
 

line()
print("Testing Atcc13...")
rc = 2.71            
dc = 0.174           
pr = 1e-5            
pd = 5e-6            
px = 0.1             
rv = 55.0            
date1, date2 = 2456165.5,  0.401182685  
ra, da = sf.pymAtcc13(rc,  dc,  pr,  pd,  px,  rv, date1, date2)
print(ra,da)

line()
print("Testing Atccq...")
#eo = c_double()
#astrom = sf.pymASTROM()
astrom, eo = sf.pymApci13(date1, date2)
ra, da = sf.pymAtccq(rc,  dc,  pr,  pd,  px,  rv,  astrom)
print(ra,da)
#sf.pymASTROM.print(astrom) 


line()
print("Testing Atci13...")
ra, da, eo = sf.pymAtci13(rc,  dc,  pr,  pd,  px,  rv,  date1, date2)
print(ra, da, eo)

line()
print("Testing Atciq...")
ra, da =  sf.pymAtciq(rc,  dc,  pr,  pd,  px,  rv,  astrom)
print(ra,da)
 

line()
print("Testing Atciqz...")
ra, da =  sf.pymAtciqz(rc,  dc,  astrom) 
print(ra,da)
                                  

line()
print("Testing Atciqz...")
ra, da =  sf.pymAtciqz(rc,  dc,  astrom) 
print(ra,da)

#end = time.time() 
#print("Elapsed time: %s" % (end - start))   

ri, di  = 2.710121572969038991 , 0.1729371367218230438 
date1, date2  = 2456165.5,   0.401182685 

line()
print("Testing Atic13...")
rc, dc, eo = sf.pymAtic13(ri, di, date1, date2) 
print(rc, dc, eo)

line()
print("Testing Aticq...")
rc, dc = sf.pymAticq(ri,  di,  astrom)
print(rc, dc)


line()
print("Testing Atio13...")

ri = 2.710121572969038991 
di = 0.1729371367218230438 

utc1, utc2, dut1  = 2456384.5, 0.969254051,  0.1550675     
elong, phi, hm = -0.527800806, -1.2345856, 2738.0             
xp, yp  = 2.47230737e-7, 1.82640464e-6    
phpa, tc = 731.0, 12.8          
rh, wl = 0.59, 0.55    

aob, zob, hob, dob, rob=sf.pymAtio13(ri, di, utc1, utc2, dut1, elong, phi, hm,
                                     xp, yp, phpa, tc, rh, wl) 
print(aob, zob, hob, dob, rob)


line()
print("Testing Atioq...")

astrom = sf.pymApio13(utc1, utc2, dut1, elong, phi, hm, xp, yp, 
                      phpa, tc, rh,  wl)

aob, zob, hob, dob, rob=sf.pymAtioq(ri, di, astrom)
print(aob, zob, hob, dob, rob)

line()
print("Testing Atoc13...")

ob1 = 2.710085107986886201
ob2 = 0.1717653435758265198
 
s  = c_char_p(b"R")
rc, dc = sf.pymAtoc13(s,   ob1,   ob2,
              utc1,  utc2,  dut1, 
              elong,  phi,    hm,  xp,  yp, 
              phpa,    tc,    rh,  wl) 
print(rc, dc)

ob1 = -0.09247619879782006106 
ob2 = 0.1717653435758265198 
s  = c_char_p(b"H")
rc, dc = sf.pymAtoc13(s,   ob1,   ob2,
              utc1,  utc2,  dut1, 
              elong,  phi,    hm,  xp,  yp, 
              phpa,    tc,    rh,  wl) 
print(rc, dc)

ob1 = 0.09233952224794989993 
ob2 = 1.407758704513722461 
s  = c_char_p(b"A")
rc, dc = sf.pymAtoc13(s,   ob1,   ob2,
              utc1,  utc2,  dut1, 
              elong,  phi,    hm,  xp,  yp, 
              phpa,    tc,    rh,  wl) 
print(rc, dc)


line()
print("Testing Atoi13...")

ob1 = 2.710085107986886201
ob2 = 0.1717653435758265198
 
s  = c_char_p(b"R")
ri, di = sf.pymAtoi13(s,   ob1,   ob2,
              utc1,  utc2,  dut1, 
              elong,  phi,    hm,  xp,  yp, 
              phpa,    tc,    rh,  wl) 
print(ri, di)

ob1 = -0.09247619879782006106 
ob2 = 0.1717653435758265198 
s  = c_char_p(b"H")
ri, di = sf.pymAtoi13(s,   ob1,   ob2,
              utc1,  utc2,  dut1, 
              elong,  phi,    hm,  xp,  yp, 
              phpa,    tc,    rh,  wl) 
print(ri, di)

ob1 = 0.09233952224794989993 
ob2 = 1.407758704513722461 
s  = c_char_p(b"A")
ri, di = sf.pymAtoi13(s,   ob1,   ob2,
              utc1,  utc2,  dut1, 
              elong,  phi,    hm,  xp,  yp, 
              phpa,    tc,    rh,  wl) 
print(ri, di)

line()
print("Testing Atoiq...")

ob1 = 2.710085107986886201
ob2 = 0.1717653435758265198
 
s  = c_char_p(b"R")
ri, di = sf.pymAtoiq(s,   ob1,   ob2,  astrom)
print(ri, di)

ob1 = -0.09247619879782006106 
ob2 = 0.1717653435758265198

s  = c_char_p(b"H")
ri, di = sf.pymAtoiq(s,   ob1,   ob2,  astrom)
print(ri, di)

ob1 = 0.09233952224794989993 
ob2 = 1.407758704513722461 

s  = c_char_p(b"A")
ri, di = sf.pymAtoiq(s,   ob1,   ob2,  astrom)
print(ri, di)

###
line()
print("Testing Hd2ae...")
h = 1.1 
d = 1.2 
p = 0.3 
az, el = sf.pymHd2ae(h, d, p)
print(az, el)

line()
print("Testing Hd2pa...")

q = sf.pymHd2pa(h, d, p)
print(q)

print("Testing Ld...")
p  = np.array([-0.763276255, -0.608633767, -0.216735543]).reshape(1,3)
q  = np.array([-0.763276255, -0.608633767, -0.216735543]).reshape(1,3)
e  = np.array([ 0.76700421,   0.605629598,  0.211937094]).reshape(1,3)

bm = 0.00028574 
em = 8.91276983 
dlim = 3e-10

p1 = sf.pymLd(bm, p, q, e, em, dlim)

#print(p1) 
print(p1[0,0], p1[0,1], p1[0,2]) 

print("Testing Ldsun...")
p  = np.array([-0.763276255, -0.608633767, -0.216735543]).reshape(1,3)
e  = np.array([-0.973644023, -0.20925523,  -0.0907169552]).reshape(1,3) 
em = 0.999809214

p1 = sf.pymLdsun(p, e, em)
#print each element of matrix p1[i,j]
print(p1[0,0], p1[0,1], p1[0,2]) 

line()
print("Testing Ldn...")
#using 2-dim array to replace iauLDBody
n   = 3 
ob  = np.array([-0.974170437, -0.2115201,   -0.0917583114]).reshape(1,3)
sc  = np.array([-0.763276255, -0.608633767, -0.216735543]).reshape(1,3)
b   = np.array([[0.00028574, 3e-10, 
                -7.81014427, -5.60956681, -1.98079819,
                 0.0030723249, -0.00406995477, -0.00181335842],
               [0.00095435, 3e-9,
                0.738098796, 4.63658692,  1.9693136,
                -0.00755816922, 0.00126913722, 0.000727999001],
               [1.0, 6e-6,
                -0.000712174377, -0.00230478303, -0.00105865966,
                6.29235213e-6, -3.30888387e-7, -2.96486623e-7
                   ]]).reshape(3,8)   
sn = sf.pymLdn(n, b, ob, sc)
print(sn[0,0], sn[0,1], sn[0,2]) 

line()
print("Testing Atciqn...") 
date1 = 2456165.5 
date2 = 0.401182685 
astrom, eo = sf.pymApci13(date1, date2)
rc = 2.71
dc = 0.174
pr = 1e-5 
pd = 5e-6 
px = 0.1 
rv = 55.0 
ri, di = sf.pymAtciqn(rc,  dc,  pr,  pd,  px,  rv,  astrom,  n,  b)
print(ri, di)

line()
print("Testing Aticqn...") 
#astrom, eo = sf.pymApci13(date1, date2)
ri = 2.709994899247599271 
di = 0.1728740720983623469 
rc, dc = sf.pymAticqn(ri,  di,  astrom,  n,  b)
print(rc, dc)

line()
print("Testing pymPmpx...") 
rc = 1.234
dc = 0.789
pr = 1e-5 
pd = -2e-5 
px = 1e-2 
rv = 10.0 
pmt = 8.75 

pob = np.array([0.9, 0.4 , 0.1]).reshape(1,3)
pco = sf.pymPmpx(rc, dc, pr, pd, px, rv, pmt, pob)
#print(poc)
print(pco[0,0], pco[0,1], pco[0,2]) 

line()
print("Testing pymPmsafe...") 

ra1 = 1.234 
dec1 = 0.789 
pmr1 = 1e-5 
pmd1 = -2e-5 
px1 = 1e-2 
rv1 = 10.0 
ep1a = 2400000.5 
ep1b = 48348.5625 
ep2a = 2400000.5 
ep2b = 51544.5 

ra2, dec2, pmr2, pmd2, px2, rv2 = sf.pymPmsafe(ra1,  dec1, pmr1, pmd1,
              px1,  rv1, ep1a, ep1b, ep2a, ep2b)
print(ra2, dec2, pmr2, pmd2, px2, rv2) 

line()
print("Testing pymPvtob...") 
elong = 2.0
phi = 0.5
hm = 3000.0
xp = 1e-6
yp = -0.5e-6
sp = 1e-8
theta = 5.0
pv = sf.pymPvtob(elong,  phi,  hm,  xp,  yp,  sp,  theta) 
print(pv[0,0], pv[0,1], pv[0,2])
print(pv[1,0], pv[1,1], pv[1,2])

line()
print("Testing pymRefco...") 
phpa = 800.0
tc = 10.0
rh = 0.9
wl = 0.4
refa, refb = sf.pymRefco(phpa,   tc,   rh,   wl) 
print(refa, refb)

line()
print("Testing pymTpors...") 
xi = -0.03   
eta = 0.07
a = 1.3  
b = 1.5 
a01, b01, a02, b02 = sf.pymTpors(xi,   eta,   a,   b) 
print(a01, b01, a02, b02)


line()
print("Testing pymTporv...") 
xi = -0.03   
eta = 0.07
ra = 1.3  
dec = 1.5 
 
v = sf.pymS2c(ra, dec)

v01, v02 = sf.pymTporv(xi,   eta,   v)
#print(v01, v02)
print(v01[0,0], v01[0,1], v01[0,2])
print(v02[0,0], v02[0,1], v02[0,2]) 

line()
print("Testing pymTpsts...")
xi = -0.03  
eta = 0.07  
a0 = 2.3    
b0 = 1.5    

a, b = sf.pymTpsts(xi,   eta,   a0,   b0)
print(a, b)

line()
print("Testing pymTpstv...")
xi = -0.03   
eta = 0.07
ra = 2.3  
dec = 1.5 
 
v0 = sf.pymS2c(ra, dec)
v  = sf.pymTpstv(xi,   eta,   v0)
print(v[0,0], v[0,1], v[0,2])

line()
print("Testing pymTpxes...")
ra = 1.3  
dec = 1.55
raz = 2.3 
decz = 1.5

xi, eta = sf.pymTpxes(ra,  dec,   raz,  decz) 
print(xi, eta)  

line()
print("Testing pymTpxev...")

v  = sf.pymS2c(ra,  dec)
v0 = sf.pymS2c(raz, decz)

xi, eta = sf.pymTpxev(v,   v0) 
print(xi, eta)  

#2023-05-03   Sofa Tools for Earth Attitude 
line()
print("Testing pymAnp...")
a = -0.1
b = sf.pymAnp(a)
print(b)

line()
print("Testing pymBi00...")
dpsibi, depsbi, dra = sf.pymBi00() 
print(dpsibi, depsbi, dra)

line()
print("Testing pymBpn2xy...")

rbpn = np.array([[9.999962358680738e-1, -2.516417057665452e-3, -1.093569785342370e-3],
                [2.516462370370876e-3,  9.999968329010883e-1,  4.006159587358310e-5],
                [1.093465510215479e-3, -4.281337229063151e-5,  9.999994012499173e-1]]) 
x, y = sf.pymBpn2xy(rbpn)
#x, y = sf.pymBpn2xy_A(rbpn)
print(x, y)

line()
print("Testing pymC2i00a...")
date1, date2 = 2400000.5, 53736.0
 
rc2i = sf.pymC2i00a(date1, date2)
print(rc2i)
#print(rc2i[0,0], rc2i[0,1], rc2i[0,2])
#print(rc2i[1,0], rc2i[1,1], rc2i[1,2])
#print(rc2i[2,0], rc2i[2,1], rc2i[2,2]) 

line()
print("Testing pymC2i00b...")
date1, date2 = 2400000.5, 53736.0
 
rc2i = sf.pymC2i00b(date1, date2)
print(rc2i)

line()
print("Testing pymC2ibpn...")
date1, date2 = 2400000.5, 50123.9999

rc2i = sf.pymC2ibpn(date1, date2, rbpn)
print(rc2i)

line()
print("Testing pymC2t00a...")
tta,  ttb = 2400000.5, 53736.0
uta,  utb = 2400000.5, 53736.0
xp,   yp  = 2.55060238e-7, 1.860359247e-6 
rc2t = sf.pymC2t00a(tta,  ttb,  uta,  utb,  xp,  yp)
print(rc2t)

line()
print("Testing pymC2t00b...")
 
rc2t = sf.pymC2t00b(tta,  ttb,  uta,  utb,  xp,  yp)
print(rc2t)

line()
print("Testing pymC2teqx...")
gst = 1.754166138040730516
rbpn = np.array([[0.9999989440476103608,  -0.1332881761240011518e-2, -0.5790767434730085097e-3],
                [0.1332858254308954453e-2, 0.9999991109044505944,    -0.4097782710401555759e-4],
                [0.5791308472168153320e-3, 0.4020595661593994396e-4,  0.9999998314954572365]]) 

rpom = np.array([[0.9999999999999674705,  -0.1367174580728847031e-10, 0.2550602379999972723e-6],
                [0.1414624947957029721e-10, 0.9999999999982694954,   -0.1860359246998866338e-5],
                [-0.2550602379741215275e-6, 0.1860359247002413923e-5, 0.9999999999982369658]])  


rc2t = sf.pymC2teqx(rbpn,  gst,  rpom) 
print(rc2t)

line()
print("Testing pymC2tpe...")
deps, dpsi = 0.4090789763356509900, -0.9630909107115582393e-5 
rc2t = sf.pymC2tpe(tta,  ttb,  uta,  utb,  dpsi,  deps,  xp,  yp)
print(rc2t)
#print(rc2t[0,0], rc2t[0,1], rc2t[0,2])
#print(rc2t[1,0], rc2t[1,1], rc2t[1,2])
#print(rc2t[2,0], rc2t[2,1], rc2t[2,2]) 

line()
print("Testing pymC2txy...")
 
x, y = 0.5791308486706011000e-3,  0.4020579816732961219e-4
rc2t = sf.pymC2txy(tta,  ttb,  uta,  utb,  x,  y,  xp,  yp)
print(rc2t)

line()
print("Testing pymEe00...")
date1, date2 = 2400000.5, 53736.0
epsa =  0.4090789763356509900 
dpsi = -0.9630909107115582393e-5 
d = sf.pymEe00(date1,  date2,   epsa,   dpsi)
print(d)


line()
print("Testing pymEe00a...")
date1, date2 = 2400000.5, 53736.0
d = sf.pymEe00a(date1,  date2)
print(d)


line()
print("Testing pymEe00b...")
date1, date2 = 2400000.5, 53736.0
d = sf.pymEe00b(date1,  date2)
print(d)

line()
print("Testing pymEe06a...")
date1, date2 = 2400000.5, 53736.0
d = sf.pymEe06a(date1,  date2)
print(d)

line()
print("Testing pymEect00...")
date1, date2 = 2400000.5, 53736.0
d = sf.pymEect00(date1,  date2)
print(d)

line()
print("Testing pymEo06a...")
date1, date2 = 2400000.5, 53736.0
d = sf.pymEo06a(date1,  date2)
print(d)

line()
print("Testing pymEors...")

rnpb = rbpn
s  = -0.1220040848472271978e-7 
eo = sf.pymEors(rnpb,  s)
print(eo)


line()
print("Testing pymEqeq94...")
date1, date2 = 2400000.5, 41234.0
d = sf.pymEqeq94(date1,  date2)
print(d)

line()
print("Testing pymEra00...")
date1, date2 = 2400000.5, 54388.0
d = sf.pymEra00(date1,  date2)
print(d)

line()
print("Testing pymFad03...")
t = 0.8
d = sf.pymFad03(t)
print(d)

line()
print("Testing pymFae03...")
 
d = sf.pymFae03(t)
print(d)

line()
print("Testing pymFaf03...")
 
d = sf.pymFaf03(t)
print(d)

line()
print("Testing pymFaju03...")
d = sf.pymFaju03(t)
print(d)


line()
print("Testing pymFal03...")
 
d = sf.pymFal03(t)
print(d)

line()
print("Testing pymFalp03...")
 
d = sf.pymFalp03(t)
print(d)

line()
print("Testing pymFama03...")
 
d = sf.pymFama03(t)
print(d)

line()
print("Testing pymFame03...")
 
d = sf.pymFame03(t)
print(d)

line()
print("Testing pymFane03...")
 
d = sf.pymFane03(t)
print(d)

line()
print("Testing pymFaom03...")
 
d = sf.pymFaom03(t)
print(d)

line()
print("Testing pymFapa03...")
 
d = sf.pymFapa03(t)
print(d)


line()
print("Testing pymFasa03...")
 
d = sf.pymFasa03(t)
print(d)

line()
print("Testing pymFaur03...")
 
d = sf.pymFaur03(t)
print(d)

line()
print("Testing pymFave03...")
 
d = sf.pymFave03(t)
print(d)


line()
print("Testing pymFw2m...")

gamb = -0.2243387670997992368e-5 
phib =  0.4091014602391312982 
psi  = -0.9501954178013015092e-3 
eps  =  0.4091014316587367472 
 
r= sf.pymFw2m(gamb,  phib,  psi,  eps)
print(r)

line()
print("Testing pymFw2xy...")
x, y = sf.pymFw2xy(gamb,   phib,  psi,  eps) 
print(x, y)

line()
print("Testing pymGmst00...")
t = sf.pymGmst00(uta,  utb,  tta,  ttb)
print(t)

line()
print("Testing pymGmst06...")
t = sf.pymGmst06(uta,  utb,  tta,  ttb)
print(t)

line()
print("Testing pymGmst82...")
t = sf.pymGmst82(uta,  utb)
print(t)

line()
print("Testing pymGst00a...")
t = sf.pymGst00a(uta,  utb,  tta,  ttb)
print(t)

line()
print("Testing pymGst00b...")
t = sf.pymGst00b(uta,  utb)
print(t)

line()
print("Testing pymGst06...")
t = sf.pymGst06(uta,  utb, tta,  ttb,  rnpb)
print(t)

line()
print("Testing pymGst06a...")
t = sf.pymGst06a(uta,  utb, tta,  ttb)
print(t)

line()
print("Testing pymGst94...")
t = sf.pymGst94(uta,  utb)
print(t)

line()
print("Testing pymNum00a...")
rmatn = sf.pymNum00a(uta,  utb)
print(rmatn)

line()
print("Testing pymNum00b...")
rmatn = sf.pymNum00b(uta,  utb)
print(rmatn)

line()
print("Testing pymNum06a...")
rmatn = sf.pymNum06a(uta,  utb)
print(rmatn)


line()
print("Testing pymNumat...")
epsa =  0.4090789763356509900 
dpsi = -0.9630909107115582393e-5 
deps =  0.4063239174001678826e-4 
rmatn = sf.pymNumat(epsa,   dpsi,   deps)
print(rmatn)

line()
print("Testing pymNut00a...")
dpsi, deps = sf.pymNut00a(uta,  utb)
print(dpsi, deps)


line()
print("Testing pymNut00b...")
dpsi, deps = sf.pymNut00b(uta,  utb)
print(dpsi, deps)

line()
print("Testing pymNut06a...")
dpsi, deps = sf.pymNut06a(uta,  utb)
print(dpsi, deps)


line()
print("Testing pymNut80...")
dpsi, deps = sf.pymNut80(uta,  utb)
print(dpsi, deps)


line()
print("Testing pymNutm80...")
rmatn = sf.pymNutm80(uta,  utb)
print(rmatn)

line()
print("Testing pymObl06...")
uta,  utb = 2400000.5, 54388.0
eps = sf.pymObl06(uta,  utb)
print(eps)

line()
print("Testing pymObl80...")
uta,  utb = 2400000.5, 54388.0
eps = sf.pymObl80(uta,  utb)
print(eps)


line()
print("Testing pymP06e...")
uta,  utb = 2400000.5, 52541.0
eps0, psia, oma, bpa, bqa, pia, bpia,  epsa, chia, za, zetaa, thetaa, pa, gam, phi, psi = sf.pymP06e(uta,  utb)
print(eps0, psia, oma, bpa, bqa, pia, bpia,  epsa, chia, za, zetaa, thetaa, pa, gam, phi, psi )


line()
print("Testing pymPb06...")
uta,  utb = 2400000.5, 50123.9999
bzeta,  bz,  btheta  = sf.pymPb06(uta,  utb)
print(bzeta,  bz,  btheta)

line()
print("Testing pymPn06a...")
uta,  utb = 2400000.5, 53736.0 
dpsi,   deps,   epsa,  rb,  rp,  rbp,  rn,  rbpn  = sf.pymPn06a(uta,  utb)
print(dpsi,   deps,   epsa)
print(rb)  
print(rp)  
print(rbp)
print(rn)
print(rbpn)

line()
print("Testing pymPnm00a...")
uta,  utb = 2400000.5, 50123.9999
rbpn = sf.pymPnm00a(uta,  utb)
print(rbpn)

line()
print("Testing pymPnm00b...")
uta,  utb = 2400000.5, 50123.9999
rbpn = sf.pymPnm00b(uta,  utb)
print(rbpn)

line()
print("Testing pymPnm06a...")
uta,  utb = 2400000.5, 50123.9999
rbpn = sf.pymPnm06a(uta,  utb)
print(rbpn)

line()
print("Testing pymPnm80...")
#uta,  utb = 2400000.5, 50123.9999
rbpn = sf.pymPnm80(uta,  utb)
print(rbpn)

line()
print("Testing pymPom00...")
xp,   yp,   sp =  2.55060238e-7, 1.860359247e-6,  -0.1367174580728891460e-10
rpom = sf.pymPom00(xp,   yp,   sp)
print( rpom)

line()
print("Testing pymPr00...")
uta,  utb = 2400000.5, 53736.0
dpsipr, depspr = sf.pymPr00(uta,  utb)
print(dpsipr, depspr)


line()
print("Testing pymPrec76...")
uta01,  uta02 = 2400000.5, 33282.0
uta11,  uta12 = 2400000.5, 51544.0
zeta, z, theta = sf.pymPrec76(uta01,  uta02,  uta11,  uta12)
print(zeta, z, theta)


line()
print("Testing pymS00...")
uta,  utb = 2400000.5, 53736.0
x, y  = 0.5791308486706011000e-3, 0.4020579816732961219e-4
s     = sf.pymS00(uta,  utb,   x,   y)
print(s)

line()
print("Testing pymS00a...")
uta,  utb = 2400000.5, 52541.0
s     = sf.pymS00a(uta,  utb)
print(s)

line()
print("Testing pymS00b...")
uta,  utb = 2400000.5, 52541.0
s     = sf.pymS00b(uta,  utb)
print(s)

line()
print("Testing pymS06...")
uta,  utb = 2400000.5, 53736.0
x, y  = 0.5791308486706011000e-3, 0.4020579816732961219e-4
s     = sf.pymS06(uta,  utb,   x,   y)
print(s)

line()
print("Testing pymS06a...")
uta,  utb = 2400000.5, 52541.0
s     = sf.pymS06a(uta,  utb)
print(s)


line()
print("Testing pymSp00...")
uta,  utb = 2400000.5, 52541.0
s     = sf.pymSp00(uta,  utb)
print(s)

line()
print("Testing pymXy06...")
uta,  utb = 2400000.5, 53736.0
x, y      = sf.pymXy06(uta,  utb)
print(x, y)

line()
print("Testing pymXys00a...")
uta,  utb = 2400000.5, 53736.0
x, y, s   = sf.pymXys00a(uta,  utb)
print(x, y, s)

line()
print("Testing pymXys00b...")
uta,  utb = 2400000.5, 53736.0
x, y, s   = sf.pymXys00b(uta,  utb)
print(x, y, s)

line()
print("Testing pymXys06a...")
uta,  utb = 2400000.5, 53736.0
x, y, s   = sf.pymXys06a(uta,  utb)
print(x, y, s)

line()
print("Testing pymA2tf...")
ndp = 4
angle = -3.01234
sign = c_char()
sign, v = sf.pymA2tf(ndp, angle)
print(sign,v)

line()
print("Testing pymAf2a...")
#single char using b... 
#s = c_char(b'-')
s =  b'-' 
ideg, iamin, asec = 45, 13, 27.2  
rad = sf.pymAf2a(s,   ideg,   iamin,   asec)
print(rad)

line()
print("Testing pymAnpm...")
#a = 2.283185307179586477
a  = -4.0
npi = sf.pymAnpm(a)
print(npi)

line()
print("Testing pymC2s...")
#a = 2.283185307179586477
#p  =  np.array([[100.0, -50.0, 25.0]])
p  =  np.array([100.0, -50.0, 25.0]).reshape(1,3)
theta, phi = sf.pymC2s(p)
print(theta, phi)


line()
print("Testing pymCp...")
 
p  =  np.array([0.3, 1.2, -2.5]).reshape(1,3)
c  = sf.pymCp(p)
#c1 = sf.pymCp_A(p)
print(c)
#print(c1)

line()
print("Testing pymCpv...")
 
pv  =  np.array([[0.3, 1.2, -2.5], [-0.5, 3.1, 0.9]]).reshape(2,3)
c   = sf.pymCpv(pv)
print(c)
#c   = sf.pymCpv_A(pv)
#print(c)

line()
print("Testing pymCr...")
r  =  np.array([[2.0, 3.0, 2.0], [3.0, 2.0, 3.0], [3.0, 4.0, 5.0]]).reshape(3,3)
c   = sf.pymCr(r)
print(c)
#c   = sf.pymCr_A(r)
#print(c)

line()
print("Testing pymIr...")
#r   = sf.pymIr().tolist() 
r   = sf.pymIr()
print(r)

line()
print("Testing pymP2pv...")
p  =  np.array([0.25, 1.2, 3.0]).reshape(1,3)
pv   = sf.pymP2pv(p)
print(pv)

line()
print("Testing pymP2s...")
p   =  np.array([100.0, -50.0, 25.0]).reshape(1,3)
theta, phi, r  = sf.pymP2s(p)
print(theta, phi, r)

line()
print("Testing pymPap...")
a   =  np.array([1.0,   0.1,  0.2]).reshape(1,3)
b   =  np.array([-3.0, 1e-3,  0.2]).reshape(1,3)
theta   = sf.pymPap(a, b)
print(theta)


line()
print("Testing pymPdp...")
a   =  np.array([2.0,   2.0,  3.0]).reshape(1,3)
b   =  np.array([1.0,   3.0,  4.0]).reshape(1,3)
adb   = sf.pymPdp(a, b)
print(adb)
a   =  np.array([2.0,   2.0,  3.0]) 
b   =  np.array([1.0,   3.0,  4.0]) 
adb   = sf.pymPdp_A(a, b)
print(adb)


line()
print("Testing pymPm...")
a   =  np.array([0.3,   1.2,  -2.5]).reshape(1,3)
an  =  sf.pymPm(a)
print(an)
a   =  np.array([0.3,   1.2,  -2.5]).reshape(1,3)
an  =  sf.pymPm_A(a)
print(an)


line()
print("Testing pymPmp...")
a   =  np.array([2.0,   2.0,  3.0]).reshape(1,3)
b   =  np.array([1.0,   3.0,  4.0]).reshape(1,3)
amb   = sf.pymPmp(a, b)
print(amb)

#amb   = sf.pymPmp_A(a, b)
#print(amb)


line()
print("Testing pymPn...")
p     =  np.array([0.3,   1.2,  -2.5]).reshape(1,3)
r, u  = sf.pymPn(p)
print(r, u)
p     =  np.array([0.3,   1.2,  -2.5]).reshape(1,3)
r, u  =  sf.pymPn_A(p)
print(r, u)

line()
print("Testing pymPpp...")
a   =  np.array([2.0,   2.0,  3.0]).reshape(1,3)
b   =  np.array([1.0,   3.0,  4.0]).reshape(1,3)
amb   = sf.pymPpp(a, b)
print(amb)
#amb   = sf.pymPpp_A(a, b)
#print(amb)

line()
print("Testing pymPpsp...")
 
s   = 5.0
apsb   = sf.pymPpsp(a, s, b)
print(apsb)
#apsb   = sf.pymPpsp_A(a, s, b)
#print(apsb)

line()
print("Testing pymPv2p...")
pv  =  np.array([[0.3, 1.2, -2.5], [-0.5, 3.1, 0.9]]).reshape(2,3)
p    = sf.pymPv2p(pv)
print(p)


line()
print("Testing pymPv2s...")
pv = np.array([[-0.4514964673880165, 0.03093394277342585, 0.05594668105108779], 
               [1.292270850663260e-5, 2.652814182060692e-6,2.568431853930293e-6]]).reshape(2,3)
theta, phi, r,  td, pd, rd = sf.pymPv2s(pv)
print(theta, phi, r,  td, pd, rd)

line()
print("Testing pymPvdpv...")
 
a   =  np.array([[2.0,   2.0,  3.0], [6.0,  0.0,  4.0]]).reshape(2,3)
b   =  np.array([[1.0,   3.0,  4.0], [0.0,  2.0,  8.0]]).reshape(2,3)
adb = sf.pymPvdpv(a, b)
print(adb)

line()
print("Testing pymPvm...")
pv   =  np.array([[0.3, 1.2, -2.5], [0.45, -0.25, 1.1]]).reshape(2,3)
r, s =  sf.pymPvm(pv)
print(r, s)

line()
print("Testing pymPvmpv...")
a   =  np.array([[2.0,   2.0,  3.0], [5.0,   6.0,  3.0]]).reshape(2,3)
b   =  np.array([[1.0,   3.0,  4.0], [3.0,   2.0,  1.0]]).reshape(2,3)
amb =  sf.pymPvmpv(a, b)
print(amb)
#amb =  sf.pymPvmpv_A(a, b)
#print(amb)


line()
print("Testing pymPvppv...")
a   =  np.array([[2.0,   2.0,  3.0], [5.0,   6.0,  3.0]]).reshape(2,3)
b   =  np.array([[1.0,   3.0,  4.0], [3.0,   2.0,  1.0]]).reshape(2,3)
apb =  sf.pymPvppv(a, b)
print(apb)
#apb =  sf.pymPvppv_A(a, b)
#print(apb)

line()
print("Testing pymPvu...")
pv  =  np.array([[126668.5912743160734, 2136.792716839935565, -245251.2339876830229],
                 [-0.4051854035740713039e-2, -0.6253919754866175788e-2, 0.1189353719774107615e-1]]).reshape(2,3)
dt  = 2920.0 
upv =  sf.pymPvu(dt, pv)
print(upv)


line()
print("Testing pymPvup...")
p =  sf.pymPvup(dt, pv)
print(p)


line()
print("Testing pymPvxpv...")
a   =  np.array([[2.0,   2.0,  3.0], [6.0,   0.0,  4.0]]).reshape(2,3)
b   =  np.array([[1.0,   3.0,  4.0], [0.0,   2.0,  8.0]]).reshape(2,3)   
axb =  sf.pymPvxpv(a, b)
print(axb)
 
#axb =  sf.pymPvxpv_A(a, b)
#print(axb)

line()
print("Testing pymPxp...")
a   =  np.array([2.0,   2.0,  3.0]).reshape(1,3)
b   =  np.array([1.0,   3.0,  4.0]).reshape(1,3)
axb   = sf.pymPxp(a, b)
print(axb)
axb   = sf.pymPxp_A(a, b)
print(axb)


line()
print("Testing pymRm2v...")
r   = np.array([[0.0, -0.8, -0.6], [0.8, -0.36, 0.48], [0.6, 0.48, -0.64]]).reshape(3,3)
w   = sf.pymRm2v(r)
print(w)


line()
print("Testing pymRv2m...")
w   = np.array([0.0, 1.41371669, -1.88495559]).reshape(1,3)
r   = sf.pymRv2m(w)
print(r)

line()
print("Testing pymRx...")
phi = 0.3456789
r   = np.array([[2.0, 3.0, 2.0], [3.0, 2.0, 3.0], [3.0, 4.0, 5.0]]).reshape(3,3)
rx  = sf.pymRx(phi, r)
print(rx)

line()
print("Testing pymRxp...")
r  = np.array([[2.0, 3.0, 2.0], [3.0, 2.0, 3.0], [3.0, 4.0, 5.0]]).reshape(3,3)
p  = np.array([0.2,  1.5, 0.1]).reshape(1,3)
rp = sf.pymRxp(r, p)
print(rp)
#p  = np.array([0.2,  1.5, 0.1]) 
#rp = sf.pymRxp_A(r, p)
#print(rp)


line()
print("Testing pymRxpv...")
r   = np.array([[2.0, 3.0, 2.0], [3.0, 2.0, 3.0], [3.0, 4.0, 5.0]]).reshape(3,3)
pv  = np.array([[0.2,  1.5,  0.1], [1.5,  0.2,  0.1]]).reshape(2, 3)
rpv = sf.pymRxpv(r, pv)
print(rpv)
 
pv  = pv.T 
rpv = sf.pymRxpv_A(r, pv)
print(rpv)


line()
print("Testing pymRxr...")
a   = np.array([[2.0, 3.0, 2.0], [3.0, 2.0, 3.0], [3.0, 4.0, 5.0]]).reshape(3,3)
b   = np.array([[1.0, 2.0, 2.0], [4.0, 1.0, 1.0], [3.0, 0.0, 1.0]]).reshape(3,3)
 
atb = sf.pymRxr(a, b)
print(atb)
atb = sf.pymRxr_A(a, b)
print(atb)
 
line()
print("Testing pymRy...")
theta = 0.3456789
r   = np.array([[2.0, 3.0, 2.0], [3.0, 2.0, 3.0], [3.0, 4.0, 5.0]]).reshape(3,3)
ry  = sf.pymRy(theta, r)
print(ry)
#r   = np.array([[2.0, 3.0, 2.0], [3.0, 2.0, 3.0], [3.0, 4.0, 5.0]]).reshape(3,3)
#ry  = sf.pymRy_A(theta, r)
#print(ry)

line()
print("Testing pymRz...")
psi = 0.3456789
r   = np.array([[2.0, 3.0, 2.0], [3.0, 2.0, 3.0], [3.0, 4.0, 5.0]]).reshape(3,3)
ry  = sf.pymRz(psi, r)
print(ry)

#r   = np.array([[2.0, 3.0, 2.0], [3.0, 2.0, 3.0], [3.0, 4.0, 5.0]]).reshape(3,3)
#ry  = sf.pymRz_A(psi, r)
#print(ry)

line()
print("Testing pymS2p...")
theta, phi, r = -3.21, 0.123, 0.456
p  = sf.pymS2p(theta, phi, r)
print(p)

line()
print("Testing pymS2pv...")
theta, phi, r, td, pd, rd = -3.21, 0.123, 0.456, -7.8e-6, 9.01e-6, -1.23e-5
pv  = sf.pymS2pv(theta, phi, r, td, pd, rd)
print(pv)

line()
print("Testing pymS2xpv...")
s1, s2 = 2.0, 3.0
pv  = np.array([[0.3,  1.2,  -2.5], [0.5,  2.3,  -0.4]]).reshape(2,3)
spv = sf.pymS2xpv(s1, s2, pv)
print(spv)

#spv = sf.pymS2xpv_A(s1, s2, pv)
#print(spv)


line()
print("Testing pymSepp...")

a   = np.array([ 1.0, 0.1,  0.2]).reshape(1,3)
b   = np.array([-3.0, 1e-3, 0.2]).reshape(1,3)   
s   = sf.pymSepp(a, b)
print(s)

#a   = np.array([ 1.0, 0.1,  0.2]) 
#b   = np.array([-3.0, 1e-3, 0.2])   
#s   = sf.pymSepp_A(a, b)
#print(s)

line()
print("Testing pymSeps...")
 
al, ap, bl, bp = 1.0, 0.1, 0.2, -3.0 
s   = sf.pymSeps(al,   ap,   bl,   bp)
print(s)

line()
print("Testing pymSxp...")
s1  = 2.0 

p  = np.array([0.3,  1.2,  -2.5]).reshape(1,3)
sp = sf.pymSxp(s1, p)
print(sp)

#p  = np.array([0.3,  1.2,  -2.5]) 
#sp = sf.pymSxp_A(s1, p)
#print(sp)

line()
print("Testing pymSxpv...")
s1 = 2.0 

pv  = np.array([[0.3,  1.2,  -2.5], [0.5,  3.2,  -0.7]]).reshape(2,3)
spv = sf.pymSxpv(s1, pv)
print(spv)


line()
print("Testing pymTf2a...")

s  = c_char(b'+')
ihour, imin, sec =  4, 58, 20.2
rad = sf.pymTf2a(s, ihour, imin, sec) 
print(rad)
s  = c_wchar('+')
rad = sf.pymTf2a_A(s, ihour, imin, sec) 
print(rad)

line()
print("Testing pymTr...")
r  =  np.array([[2.0, 3.0, 2.0], [3.0, 2.0, 3.0], [3.0, 4.0, 5.0]]).reshape(3,3)
rt = sf.pymTr(r)
print(rt)

rt = sf.pymTr_A(r)
print(rt)

line()
print("Testing pymTrxp...")
r  =  np.array([[2.0, 3.0, 2.0], [3.0, 2.0, 3.0], [3.0, 4.0, 5.0]]).reshape(3,3)
p  =  np.array([0.2,  1.5, 0.1]).reshape(1,3)
trp = sf.pymTrxp(r, p)
print(trp)

#p  =  np.array([0.2,  1.5, 0.1]) 
#trp = sf.pymTrxp_A(r, p)
#print(trp)

line()
print("Testing pymTrxpv...")
r  =  np.array([[2.0, 3.0, 2.0], [3.0, 2.0, 3.0], [3.0, 4.0, 5.0]]).reshape(3,3)
pv  = np.array([[0.2,  1.5, 0.1], [1.5, 0.2, 0.1]]).reshape(2,3)
trpv = sf.pymTrxpv(r, pv)
print(trpv)
 
#trpv = sf.pymTrxpv_A(r, pv)
#print(trpv)

line()
print("Testing pymZp...")
#print(p) 
p = sf.pymZp(p)
print(p)

#p = sf.pymZp_A(p)
#print(p)

line()
print("Testing pymZpv...")
#print(pv) 
pv = sf.pymZpv(pv)
print(pv)

#pv = sf.pymZpv_A(pv)
#print(pv)


line()
print("Testing pymZr...")
print(r) 
r = sf.pymZr(r)
print(r)

r = sf.pymZr_A(r)
print(r)

###########
line()
print("Testing pymBp00...")
date1 = 2400000.5
date2 = 50123.9999
rb, rp, rbp = sf.pymBp00(date1, date2)
print(rb,rp,rbp)

line()
print("Testing pymBp06...")
date1 = 2400000.5
date2 = 50123.9999
rb, rp, rbp = sf.pymBp06(date1, date2)
print(rb,rp,rbp)

line()
print("Testing pymC2i06a...")
date1 = 2400000.5
date2 = 53736.0
rc2i = sf.pymC2i06a(date1, date2)
print(rc2i)

line()
print("Testing pymC2ixy...")
date1 = 2400000.5
date2 = 53736.0
x = 0.5791308486706011000e-3
y = 0.4020579816732961219e-4
rc2i = sf.pymC2ixy(date1, date2, x, y)
print(rc2i)

line()
print("Testing pymC2ixys...")
x = 0.5791308486706011000e-3
y = 0.4020579816732961219e-4
s = -0.1220040848472271978e-7
rc2i = sf.pymC2ixys(x, y, s)
print(rc2i)

line()
print("Testing pymC2t06a...")
tta = 2400000.5
uta = 2400000.5
ttb = 53736.0
utb = 53736.0
xp = 2.55060238e-7
yp = 1.860359247e-6
rc2t = sf.pymC2t06a(tta, ttb, uta, utb, xp, yp)
print(rc2t)

line()
print("Testing pymC2tcio...")
rc2i = np.array([[0.9999998323037164738, 0.5581526271714303683e-9, -0.5791308477073443903e-3],
                 [-0.2384266227524722273e-7, 0.9999999991917404296, -0.4020594955030704125e-4],
                 [0.5791308472168153320e-3, 0.4020595661593994396e-4, 0.9999998314954572365]]).reshape(3,3)
era = 1.75283325530307
rpom = np.array([[0.9999999999999674705, -0.1367174580728847031e-10, 0.2550602379999972723e-6],
                 [0.1414624947957029721e-10, 0.9999999999982694954, -0.1860359246998866338e-5],
                 [-0.2550602379741215275e-6, 0.1860359247002413923e-5, 0.9999999999982369658]]).reshape(3,3)
rc2t = sf.pymC2tcio(rc2i, era, rpom)
print(rc2t)

line()
print("Testing pymEceq06...")
date1 = 2456165.5
date2 = 0.401182685
dl = 5.1
db = -0.9
dr, dd = sf.pymEceq06(date1, date2, dl, db)
print(dr, dd) 

line()
print("Testing pymEcm06...")
date1 = 2456165.5
date2 = 0.401182685
rm = sf.pymEcm06(date1, date2)
print(rm)

line()
print("Testing pymEform...")
n=1
a, f = sf.pymEform(n)
print(a, f)

line()
print("Testing pymEpv00...")
date1 = 2400000.5
date2 = 53411.52501161
pvh, pvb = sf.pymEpv00(date1, date2)
print(pvh, pvb)

line()
print("Testing pymEqec06...")
date1 = 1234.5
date2 = 2440000.5
dr = 1.234
dd = 0.987
dl, db = sf.pymEqec06(date1, date2, dr, dd)
print(dl, db)

line()
print("Testing pymFk5hip...")
r5h, s5h = sf.pymFk5hip()
print(r5h, s5h)

line()
print("Testing pymFk5hz...")
r5 =  1.76779433
d5 = -0.2917517103
date1 = 2400000.5
date2 = 54479.0
rh, dh = sf.pymFk5hz(r5, d5, date1, date2)
print(rh, dh)

line()
print("Testing pymFk45z...")
r1950 = 0.01602284975382960982
d1950 = -0.1164347929099906024
bepoch = 1954.677617625256806
r2000, d2000 = sf.pymFk45z(r1950, d1950, bepoch)
print(r2000, d2000)

line()
print("Testing pymFk52h...")
r5  =  1.76779433
d5  = -0.2917517103
dr5 = -1.91851572e-7
dd5 = -5.8468475e-6
px5 =  0.379210
rv5 = -7.6
rh, dh, drh, ddh, pxh, rvh = sf.pymFk52h(r5, d5, dr5, dd5, px5, rv5)
print(rh, dh, drh, ddh, pxh, rvh)

line()
print("Testing pymFk54z...")
r2000 = 0.02719026625066316119
d2000 = -0.1115815170738754813
bepoch = 1954.677308160316374
r1950, d1950, dr1950, dd1950 = sf.pymFk54z(r2000, d2000, bepoch)
print(r1950, d1950, dr1950, dd1950)

line()
print("Testing pymFk425...")
r1950 = 0.07626899753879587532
d1950 = -1.137405378399605780
dr1950 = 0.1973749217849087460e-4
dd1950 = 0.5659714913272723189e-5
p1950 = 0.134
v1950 = 8.7
r2000, d2000, dr2000, dd2000, p2000, v2000= sf.pymFk425(r1950, d1950, dr1950, dd1950, p1950, v1950)
print(r2000, d2000, dr2000, dd2000, p2000, v2000)

line()
print("Testing pymFk524...")
r2000 = 0.8723503576487275595
d2000 = -0.7517076365138887672
dr2000 = 0.2019447755430472323e-4
dd2000 = 0.3541563940505160433e-5
p2000 = 0.1559
v2000 = 86.87
r1950, d1950, dr1950, dd1950, p1950, v1950 = sf.pymFk524(r2000, d2000, dr2000, dd2000, p2000, v2000)
print(r1950, d1950, dr1950, dd1950, p1950, v1950)

line()
print("Testing pymG2icrs...")
dl =  5.5850536063818546461558105
db = -0.7853981633974483096156608
dr, dd = sf.pymG2icrs(dl, db)
print(dr, dd)

line()
print("Testing pymGc2gd...")
n=1
xyz=np.array([[2e6, 3e6, 5.244e6]]).reshape(1,3)
elong, phi, height = sf.pymGc2gd(n, xyz)
print(elong, phi, height)

line()
print("Testing pymGc2gde...")
a = 6378136.0
f = 0.0033528
xyz=np.array([[2e6, 3e6, 5.244e6]]).reshape(1,3)
elong, phi, height = sf.pymGc2gde(a, f, xyz)
print(elong, phi, height)

line()
print("Testing pymGd2gc...")
n=1
elong = 3.1
phi = -0.5
height = 2500.0
xyz = sf.pymGd2gc(n, elong, phi, height)
print(xyz)

line()
print("Testing pymGd2gce...")
a = 6378136.0
f = 0.0033528
elong = 3.1
phi = -0.5
height = 2500.0
xyz = sf.pymGd2gce(a, f, elong, phi, height)
print(xyz)

line()
print("Testing pymH2fk5...")
rh  =  1.767794352
dh  = -0.2917512594
drh = -2.76413026e-6
ddh = -5.92994449e-6
pxh =  0.379210
rvh = -7.6
r5, d5, dr5, dd5, px5, rv5= sf.pymH2fk5(rh, dh, drh, ddh, pxh, rvh)
print(r5, d5, dr5, dd5, px5, rv5)

line()
print("Testing pymHfk5z...")
rh =  1.767794352
dh = -0.2917512594
date1 = 2400000.5
date2 = 54479.0
r5, d5, dr5, dd5 = sf.pymHfk5z(rh, dh, date1, date2)
print(r5, d5, dr5, dd5)

line()
print("Testing pymIcrs2g...")
dr =  5.9338074302227188048671087
dd = -1.1784870613579944551540570
dl, db = sf.pymIcrs2g(dr, dd)
print(dl, db)

line()
print("Testing pymLteceq...")
epj = 2500.0
dl = 1.5
db = 0.6
dr, dd = sf.pymLteceq(epj, dl, db)
print(dr, dd)

line()
print("Testing pymLtecm...")
epj = -3000.0
rm = sf.pymLtecm(epj)
print(rm)

line()
print("Testing pymLteqec...")
epj = -1500.0
dr = 1.234
dd = 0.987
dl, db = sf.pymLteqec(epj, dr, dd)
print(dl ,db)

line()
print("Testing pymLtp...")
epj = 1666.666
rp = sf.pymLtp(epj)
print(rp)

line()
print("Testing pymLtpb...")
epj = 1666.666
rp = sf.pymLtpb(epj)
print(rp)

line()
print("Testing pymLtpecl...")
epj = -1500.0
vec = sf.pymLtpecl(epj)
print(vec)

line()
print("Testing pymLtpequ...")
epj = -2500.0
veq = sf.pymLtpequ(epj)
print(veq)

line()
print("Testing pymMoon98...")
date1 = 2400000.5
date2 = 43999.9
pv = sf.pymMoon98(date1, date2)
print(pv)

line()
print("Testing pymPfw06...")
date1 = 2400000.5
date2 = 50123.9999
gamb, phib, psib, epsa = sf.pymPfw06(date1, date2)
print(gamb, phib, psib, epsa)

line()
print("Testing pymPlan94...")
date1 = 2400000.5
date2 = 43999.9
npp = 1
pv = sf.pymPlan94(date1, date2, npp)
print(pv)

line()
print("Testing pymPmat00...")
date1 = 2400000.5
date2 = 50123.9999
rbp = sf.pymPmat00(date1, date2)
print(rbp)

line()
print("Testing pymPmat06...")
date1 = 2400000.5
date2 = 50123.9999
rbp = sf.pymPmat06(date1, date2)
print(rbp)

line()
print("Testing pymPmat76...")
date1 = 2400000.5
date2 = 50123.9999
rbp = sf.pymPmat76(date1, date2)
print(rbp)

line()
print("Testing pymPn00...")
date1 = 2400000.5
date2 = 53736.0
dpsi = -0.9632552291149335877e-5
deps =  0.4063197106621141414e-4
epsa, rb, rp, rbp, rn, rbpn = sf.pymPn00(date1, date2, dpsi, deps)
print(epsa, rb, rp, rbp, rn, rbpn)

line()
print("Testing pymPn00a...")
date1 = 2400000.5
date2 = 53736.0
dpsi, deps, epsa, rb, rp, rbp, rn, rbpn = sf.pymPn00a(date1, date2)
print(epsa, rb, rp, rbp, rn, rbpn)

line()
print("Testing pymPn00b...")
date1 = 2400000.5
date2 = 53736.0
dpsi, deps, epsa, rb, rp, rbp, rn, rbpn = sf.pymPn00b(date1, date2)
print(epsa, rb, rp, rbp, rn, rbpn)

line()
print("Testing pymPn06...")
date1 = 2400000.5
date2 = 53736.0
dpsi = -0.9632552291149335877e-5
deps =  0.4063197106621141414e-4
epsa, rb, rp, rbp, rn, rbpn = sf.pymPn06(date1, date2, dpsi, deps)
print(epsa, rb, rp, rbp, rn, rbpn)

line()
print("Testing pymPvstar...")
pv = np.array([[126668.5912743160601, 2136.792716839935195, -245251.2339876830091],
               [-0.4051854035740712739e-2, -0.6253919754866173866e-2, 0.1189353719774107189e-1]]).reshape(2,3)
ra, dec, pmr, pmd, px, rv = sf.pymPvstar(pv)
print(ra, dec, pmr, pmd, px, rv)

line()
print("Testing pymStarpm...")
ra1 =   0.01686756
dec1 = -1.093989828
pmr1 = -1.78323516e-5
pmd1 =  2.336024047e-6
px1 =   0.74723
rv1 = -21.6
ep1a = 2400000.5
ep1b = 50083.0
ep2a = 2400000.5
ep2b = 53736.0
ra2, dec2, pmr2, pmd2, px2, rv2 = sf.pymStarpm(ra1, dec1, pmr1, pmd1, px1, rv1, ep1a, ep1b, ep2a, ep2b)
print(ra2, dec2, pmr2, pmd2, px2, rv2)

line()
print("Testing pymStarpv...")
ra =   0.01686756
dec = -1.093989828
pmr = -1.78323516e-5
pmd =  2.336024047e-6
px =   0.74723
rv = -21.6
pv = sf.pymStarpv(ra, dec, pmr, pmd, px, rv)
print(pv)





