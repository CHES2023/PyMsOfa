import numpy  as  np
import numpy.ctypeslib  as  nt
import PyMsOfa  as  sf
import time

def   line():
        print('-'*60)

line()
print("Testing pymA2af...")
NDP = 4
ANGLE = 2.345
S, IDMSF = sf.pymA2af(NDP, ANGLE)
print(S, IDMSF)

line()
print("Testing pymA2tf...")
NDP=4
ANGLE=3.01234
S, IHMSF = sf.pymA2tf(NDP, ANGLE)
print(S, IHMSF)

line()
print("Testing pymAb...")
PNAT = [-0.76321968546737951, -0.60869453983060384, -0.21676408580639883]
V = [2.1044018893653786E-5, -8.9108923304429319E-5, -3.8633714797716569E-5]
S = 0.99980921395708788
BM1 = 0.99999999506209258
PPR = sf.pymAb(PNAT, V, S, BM1)
print(PPR)

line()
print("Testing pymAe2hd...")
a = 5.5
e = 1.1
p = 0.7
h, d = sf.pymAe2hd(a, e, p)
print(h, d)

line()
print("Testing pymAf2a...")
S = '-'
IDEG = 45
IAMIN = 13
ASEC = 27.2
RAD, J = sf.pymAf2a(S, IDEG, IAMIN, ASEC)
print(RAD, J)

line()
print("Testing pymAnp...")
A=-0.1
print(sf.pymAnp(A))

line()
print("Testing pymAnpm...")
A=-4.0
print(sf.pymAnpm(A))

line()
print("Testing pymApcg...")
ASTROM=[0 for i in range(30)]
DATE1 = 2456165.5
DATE2 = 0.401182685
EBPV = [[0.901310875, -0.417402664, -0.180982288],
        [0.00742727954, 0.0140507459, 0.00609045792]]
EHP = [0.903358544, -0.415395237, -0.180084014]
print(sf.pymApcg(DATE1, DATE2, EBPV, EHP, ASTROM))

line()
print("Testing pymApcg13...")
ASTROM=[0 for i in range(30)]
DATE1 = 2456165.5
DATE2 = 0.401182685
print(sf.pymApcg13(DATE1, DATE2, ASTROM))

line()
print("Testing pymApci...")
ASTROM=[0 for i in range(30)]
DATE1 = 2456165.5
DATE2 = 0.401182685
EBPV = [[0.901310875, -0.417402664, -0.180982288],
        [0.00742727954, 0.0140507459, 0.00609045792]]
EHP = [0.903358544, -0.415395237, -0.180084014]
X =  0.0013122272
Y = -2.92808623e-5
S =  3.05749468e-8
print(sf.pymApci(DATE1, DATE2, EBPV, EHP, X, Y, S, ASTROM))

line()
print("Testing pymApci13...")
ASTROM=[0 for i in range(30)]
DATE1 = 2456165.5
DATE2 = 0.401182685
print(sf.pymApci13(DATE1, DATE2, ASTROM))

line()
print("Testing pymApco...")
ASTROM=[0 for i in range(30)]
DATE1 = 2456384.5
DATE2 = 0.970031644
EBPV = [[-0.974170438, -0.211520082, -0.0917583024],
        [0.00364365824, -0.0154287319, -0.00668922024]]
EHP = [-0.973458265, -0.209215307, -0.0906996477]
X =  0.0013122272
Y = -2.92808623e-5
S =  3.05749468e-8
THETA = 3.14540971
ELONG = -0.527800806
PHI = -1.2345856
HM = 2738.0
XP = 2.47230737e-7
YP = 1.82640464e-6
SP = -3.01974337e-11
REFA = 0.000201418779
REFB = -2.36140831e-7
print(sf.pymApco(DATE1, DATE2, EBPV, EHP, X, Y, S, THETA, ELONG, 
                  PHI, HM, XP, YP, SP, REFA, REFB))

line()
print("Testing pymApco13...")
ASTROM=[0 for i in range(30)]
UTC1 = 2456384.5
UTC2 = 0.969254051
DUT1 = 0.1550675
ELONG = -0.527800806
PHI = -1.2345856
HM = 2738.0
XP = 2.47230737e-7
YP = 1.82640464e-6
PHPA = 731.0
TC = 12.8
RH = 0.59
WL = 0.55
print(sf.pymApco13(UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA,
                    TC, RH, WL))

line()
print("Testing pymApcs...")
ASTROM=[0 for i in range(30)]
DATE1 = 2456384.5
DATE2 = 0.970031644
PV = [[-1836024.09, 1056607.72, -5998795.26],
      [-77.0361767, -133.310856, 0.0971855934]]
EBPV = [[-0.974170438, -0.211520082, -0.0917583024],
        [0.00364365824, -0.0154287319, -0.00668922024]]
EHP = [-0.973458265, -0.209215307, -0.0906996477]
print(sf.pymApcs(DATE1, DATE2, PV, EBPV, EHP, ASTROM))

line()
print("Testing pymApcs13...")
ASTROM=[0 for i in range(30)]
DATE1 = 2456165.5
DATE2 = 0.401182685
PV = [[-6241497.16, 401346.896, -1251136.04],
      [-29.264597, -455.021831, 0.0266151194]]
print(sf.pymApcs13(DATE1, DATE2, PV, ASTROM))

line()
print("Testing pymAper...")
ASTROM=[0 for i in range(30)]
ASTROM[21] = 1.234
THETA = 5.678
print(sf.pymAper(THETA, ASTROM))

line()
print("Testing pymAper13...")
ASTROM=[0 for i in range(30)]
ASTROM[21] = 1.234
UT11 = 2456165.5
UT12 = 0.401182685
print(sf.pymAper13(UT11, UT12, ASTROM))

line()
print("Testing pymApio...")
ASTROM=[0 for i in range(30)]
SP = -3.01974337e-11
THETA = 3.14540971
ELONG = -0.527800806
PHI = -1.2345856
HM = 2738.0
XP = 2.47230737e-7
YP = 1.82640464e-6
REFA = 0.000201418779
REFB = -2.36140831e-7
print(sf.pymApio(SP, THETA, ELONG, PHI, HM, XP, YP, REFA, REFB, ASTROM))

line()
print("Testing pymApio13...")
ASTROM=[0 for i in range(30)]
UTC1 = 2456384.5
UTC2 = 0.969254051
DUT1 = 0.1550675
ELONG = -0.527800806
PHI = -1.2345856
HM = 2738.0
XP = 2.47230737e-7
YP = 1.82640464e-6
PHPA = 731.0
RH = 0.59
WL = 0.55
print(sf.pymApio13(UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL, ASTROM))

line()
print("Testing pymAtcc13...")
ASTROM=[0 for i in range(30)]
RC = 2.71
DC = 0.174
PR = 1e-5
PD = 5e-6
PX = 0.1
RV = 55.0
DATE1 = 2456165.5
DATE2 = 0.401182685
print(sf.pymAtcc13(RC, DC, PR, PD, PX, RV, DATE1, DATE2))

line()
print("Testing pymAtccq...")
ASTROM=[0 for i in range(30)]
RC = 2.71
DC = 0.174
PR = 1e-5
PD = 5e-6
PX = 0.1
RV = 55.0
DATE1 = 2456165.5
DATE2 = 0.401182685
ASTROM,EO=sf.pymApci13(DATE1, DATE2, ASTROM)
print(sf.pymAtccq(RC, DC, PR, PD, PX, RV, ASTROM))

line()
print("Testing pymAtci13...")
ASTROM=[0 for i in range(30)]
RC = 2.71
DC = 0.174
PR = 1e-5
PD = 5e-6
PX = 0.1
RV = 55.0
DATE1 = 2456165.5
DATE2 = 0.401182685
print(sf.pymAtci13(RC, DC, PR, PD, PX, RV, DATE1, DATE2))

line()
print("Testing pymAtciq...")
ASTROM=[0 for i in range(30)]
RC = 2.71
DC = 0.174
PR = 1e-5
PD = 5e-6
PX = 0.1
RV = 55.0
DATE1 = 2456165.5
DATE2 = 0.401182685
ASTROM, EO = sf.pymApci13(DATE1, DATE2, ASTROM)
print(sf.pymAtciq(RC, DC, PR, PD, PX, RV, ASTROM))

line()
print("Testing pymAtciqn...")
ASTROM=[0 for i in range(30)]
RC = 2.71
DC = 0.174
PR = 1e-5
PD = 5e-6
PX = 0.1
RV = 55.0
DATE1 = 2456165.5
DATE2 = 0.401182685
ASTROM, EO = sf.pymApci13(DATE1, DATE2, ASTROM)
N = 3
B = [[0.00028574, 3e-10, -7.81014427, -5.60956681,
      -1.98079819, 0.0030723249, -0.00406995477, -0.00181335842],
     [0.00095435, 3e-9, 0.738098796, 4.63658692, 
      1.9693136, -0.00755816922, 0.00126913722, 0.000727999001],
     [1.0, 6e-6, -0.000712174377, -0.00230478303, 
      -0.00105865966, 6.29235213e-6, -3.30888387e-7, -2.96486623e-7]]
print(sf.pymAtciqn(RC, DC, PR, PD, PX, RV, ASTROM, N, B))

line()
print("Testing pymAtciqz...")
ASTROM=[0 for i in range(30)]
RC = 2.71
DC = 0.174
DATE1 = 2456165.5
DATE2 = 0.401182685
ASTROM, EO = sf.pymApci13(DATE1, DATE2, ASTROM)
print(sf.pymAtciqz(RC, DC, ASTROM))

line()
print("Testing pymAtco13...")
ASTROM=[0 for i in range(30)]
RC = 2.71
DC = 0.174
PR = 1e-5
PD = 5e-6
PX = 0.1
RV = 55.0
UTC1 = 2456384.5
UTC2 = 0.969254051
DUT1 = 0.1550675
ELONG = -0.527800806
PHI = -1.2345856
HM = 2738.0
XP = 2.47230737e-7
YP = 1.82640464e-6
PHPA = 731.0
TC = 12.8
RH = 0.59
WL = 0.55
print(sf.pymAtco13(RC, DC, PR, PD, PX, RV, UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL))

line()
print("Testing pymAtic13...")
ASTROM=[0 for i in range(30)]
RI = 2.710121572969038991
DI = 0.1729371367218230438
DATE1 = 2456165.5
DATE2 = 0.401182685
print(sf.pymAtic13(RI, DI, DATE1, DATE2))

line()
print("Testing pymAticq...")
ASTROM=[0 for i in range(30)]
RI = 2.710121572969038991
DI = 0.1729371367218230438
DATE1 = 2456165.5
DATE2 = 0.401182685
ASTROM, EO = sf.pymApci13(DATE1, DATE2, ASTROM)
print(sf.pymAticq(RI, DI, ASTROM))

line()
print("Testing pymAticqn...")
ASTROM=[0 for i in range(30)]
RI = 2.709994899247599271
DI = 0.1728740720983623469
DATE1 = 2456165.5
DATE2 = 0.401182685
ASTROM, EO = sf.pymApci13(DATE1, DATE2, ASTROM)
B = [[0.00028574, 3e-10, -7.81014427, -5.60956681,
      -1.98079819, 0.0030723249, -0.00406995477, -0.00181335842],
     [0.00095435, 3e-9, 0.738098796, 4.63658692, 
      1.9693136, -0.00755816922, 0.00126913722, 0.000727999001],
     [1.0, 6e-6, -0.000712174377, -0.00230478303, 
      -0.00105865966, 6.29235213e-6, -3.30888387e-7, -2.96486623e-7]]
print(sf.pymAticqn(RI, DI, ASTROM, N, B))

line()
print("Testing pymAtio13...")
ASTROM=[0 for i in range(30)]
RI = 2.710121572969038991
DI = 0.1729371367218230438
UTC1 = 2456384.5
UTC2 = 0.969254051
DUT1 = 0.1550675
ELONG = -0.527800806
PHI = -1.2345856
HM = 2738.0
XP = 2.47230737e-7
YP = 1.82640464e-6
PHPA = 731.0
TC = 12.8
RH = 0.59
WL = 0.55
print(sf.pymAtio13(RI, DI, UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL))

line()
print("Testing pymAtioq...")
ASTROM=[0 for i in range(30)]
RI = 2.710121572969038991
DI = 0.1729371367218230438
UTC1 = 2456384.5
UTC2 = 0.969254051
DUT1 = 0.1550675
ELONG = -0.527800806
PHI = -1.2345856
HM = 2738.0
XP = 2.47230737e-7
YP = 1.82640464e-6
PHPA = 731.0
TC = 12.8
RH = 0.59
WL = 0.55
ASTROM, J = sf.pymApio13(UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL, ASTROM)
print(sf.pymAtioq(RI, DI, ASTROM))

line()
print("Testing pymAtoc13...")
ASTROM=[0 for i in range(30)]
UTC1 = 2456384.5
UTC2 = 0.969254051
DUT1 = 0.1550675
ELONG = -0.527800806
PHI = -1.2345856
HM = 2738.0
XP = 2.47230737e-7
YP = 1.82640464e-6
PHPA = 731.0
TC = 12.8
RH = 0.59
WL = 0.55
OB1 = 2.710085107986886201
OB2 = 0.1717653435758265198
TYPE = 'R'
print(sf.pymAtoc13(TYPE, OB1, OB2, UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL))
OB1 = -0.09247619879782006106
OB2 = 0.1717653435758265198
TYPE = 'H'
print(sf.pymAtoc13(TYPE, OB1, OB2, UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL))
OB1 = 0.09233952224794989993
OB2 = 1.407758704513722461
TYPE = 'A'
print(sf.pymAtoc13(TYPE, OB1, OB2, UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL))

line()
print("Testing pymAtoi13...")
ASTROM=[0 for i in range(30)]
UTC1 = 2456384.5
UTC2 = 0.969254051
DUT1 = 0.1550675
ELONG = -0.527800806
PHI = -1.2345856
HM = 2738.0
XP = 2.47230737e-7
YP = 1.82640464e-6
PHPA = 731.0
TC = 12.8
RH = 0.59
WL = 0.55
OB1 = 2.710085107986886201
OB2 = 0.1717653435758265198
TYPE = 'R'
print(sf.pymAtoi13(TYPE, OB1, OB2, UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL))
OB1 = -0.09247619879782006106
OB2 = 0.1717653435758265198
TYPE = 'H'
print(sf.pymAtoi13(TYPE, OB1, OB2, UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL))
OB1 = 0.09233952224794989993
OB2 = 1.407758704513722461
TYPE = 'A'
print(sf.pymAtoi13(TYPE, OB1, OB2, UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL))

line()
print("Testing pymAtoiq...")
ASTROM=[0 for i in range(30)]
UTC1 = 2456384.5
UTC2 = 0.969254051
DUT1 = 0.1550675
ELONG = -0.527800806
PHI = -1.2345856
HM = 2738.0
XP = 2.47230737e-7
YP = 1.82640464e-6
PHPA = 731.0
TC = 12.8
RH = 0.59
WL = 0.55
ASTROM, J = sf.pymApio13(UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL, ASTROM)
OB1 = 2.710085107986886201
OB2 = 0.1717653435758265198
TYPE = 'R'
print(sf.pymAtoiq(TYPE, OB1, OB2, ASTROM))
OB1 = -0.09247619879782006106
OB2 = 0.1717653435758265198
TYPE = 'H'
print(sf.pymAtoiq(TYPE, OB1, OB2, ASTROM))
OB1 = 0.09233952224794989993
OB2 = 1.407758704513722461
TYPE = 'A'
print(sf.pymAtoiq(TYPE, OB1, OB2, ASTROM))

line()
print("Testing pymBi00...")
print(sf.pymBi00())

line()
print("Testing pymBp00...")
DATE1 = 2400000.5
DATE2 = 50123.9999
print(sf.pymBp00(DATE1, DATE2))

line()
print("Testing pymBp06...")
DATE1 = 2400000.5
DATE2 = 50123.9999
print(sf.pymBp06(DATE1, DATE2))

line()
print("Testing pymBpn2xy...")
RBPN = [[9.999962358680738e-1, -2.516417057665452e-3, -1.093569785342370e-3],
        [2.516462370370876e-3, 9.999968329010883e-1, 4.006159587358310e-5],
        [1.093465510215479e-3, -4.281337229063151e-5, 9.999994012499173e-1]]
print(sf.pymBpn2xy(RBPN))

line()
print("Testing pymC2i00a...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymC2i00a(DATE1, DATE2))

line()
print("Testing pymC2i00b...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymC2i00b(DATE1, DATE2))

line()
print("Testing pymC2i06a...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymC2i06a(DATE1, DATE2))

line()
print("Testing pymC2ibpn...")
DATE1 = 2400000.5
DATE2 = 50123.9999
RBPN = [[9.999962358680738e-1, -2.516417057665452e-3, -1.093569785342370e-3],
        [2.516462370370876e-3, 9.999968329010883e-1, 4.006159587358310e-5],
        [1.093465510215479e-3, -4.281337229063151e-5, 9.999994012499173e-1]]
print(sf.pymC2ibpn(DATE1, DATE2, RBPN))

line()
print("Testing pymC2ixy...")
DATE1 = 2400000.5
DATE2 = 53736.0
X = 0.5791308486706011000e-3
Y = 0.4020579816732961219e-4
print(sf.pymC2ixy(DATE1, DATE2, X, Y))

line()
print("Testing pymC2ixys...")
X = 0.5791308486706011000e-3
Y = 0.4020579816732961219e-4
S = -0.1220040848472271978e-7
print(sf.pymC2ixys(X, Y, S))

line()
print("Testing pymC2s...")
P = [100.0, -50.0, 25.0]
print(sf.pymC2s(P))

line()
print("Testing pymC2t00a...")
TTA = 2400000.5
UTA = 2400000.5
TTB = 53736.0
UTB = 53736.0
XP = 2.55060238e-7
YP = 1.860359247e-6
print(sf.pymC2t00a(TTA, TTB, UTA, UTB, XP, YP))

line()
print("Testing pymC2t00b...")
TTA = 2400000.5
UTA = 2400000.5
TTB = 53736.0
UTB = 53736.0
XP = 2.55060238e-7
YP = 1.860359247e-6
print(sf.pymC2t00b(TTA, TTB, UTA, UTB, XP, YP))

line()
print("Testing pymC2t06a...")
TTA = 2400000.5
UTA = 2400000.5
TTB = 53736.0
UTB = 53736.0
XP = 2.55060238e-7
YP = 1.860359247e-6
print(sf.pymC2t06a(TTA, TTB, UTA, UTB, XP, YP))

line()
print("Testing pymC2tcio...")
RC2I = [[0.9999998323037164738, 0.5581526271714303683e-9, -0.5791308477073443903e-3],
        [-0.2384266227524722273e-7, 0.9999999991917404296, -0.4020594955030704125e-4],
        [0.5791308472168153320e-3, 0.4020595661593994396e-4, 0.9999998314954572365]]
ERA = 1.75283325530307
RPOM = [[0.9999999999999674705, -0.1367174580728847031e-10, 0.2550602379999972723e-6],
        [0.1414624947957029721e-10, 0.9999999999982694954, -0.1860359246998866338e-5],
        [-0.2550602379741215275e-6, 0.1860359247002413923e-5, 0.9999999999982369658]]
print(sf.pymC2tcio(RC2I, ERA, RPOM))

line()
print("Testing pymC2teqx...")
RBPN = [[0.9999989440476103608, -0.1332881761240011518e-2, -0.5790767434730085097e-3],
        [0.1332858254308954453e-2, 0.9999991109044505944, -0.4097782710401555759e-4],
        [0.5791308472168153320e-3, 0.4020595661593994396e-4, 0.9999998314954572365]]
GST = 1.754166138040730516
RPOM = [[0.9999999999999674705, -0.1367174580728847031e-10, 0.2550602379999972723e-6],
        [0.1414624947957029721e-10, 0.9999999999982694954, -0.1860359246998866338e-5],
        [-0.2550602379741215275e-6, 0.1860359247002413923e-5, 0.9999999999982369658]]
print(sf.pymC2teqx(RBPN, GST, RPOM))

line()
print("Testing pymC2tpe...")
TTA = 2400000.5
UTA = 2400000.5
TTB = 53736.0
UTB = 53736.0
DEPS = 0.4090789763356509900
DPSI = -0.9630909107115582393e-5
XP = 2.55060238e-7
YP = 1.860359247e-6
print(sf.pymC2tpe(TTA, TTB, UTA, UTB, DPSI, DEPS, XP, YP))

line()
print("Testing pymC2txy...")
TTA = 2400000.5
UTA = 2400000.5
TTB = 53736.0
UTB = 53736.0
X = 0.5791308486706011000e-3
Y = 0.4020579816732961219e-4
XP = 2.55060238e-7
YP = 1.860359247e-6
print(sf.pymC2txy(TTA, TTB, UTA, UTB, X, Y, XP, YP))

line()
print("Testing pymCal2jd...")
IY = 2003
IM  = 6
ID = 1
print(sf.pymCal2jd(IY, IM, ID))

line()
print("Testing pymCp...")
P = [0.3, 1.2, -2.5]
print(sf.pymCp(P))

line()
print("Testing pymCpv...")
PV = [[0.3, 1.2, -2.5],
      [-0.5, 3.1, 0.9]]
print(sf.pymCpv(PV))

line()
print("Testing pymCr...")
R = [[2.0, 3.0, 2.0],
      [3.0, 2.0, 3.0],
      [3.0, 4.0, 5.0]]
print(sf.pymCr(R))

line()
print("Testing pymD2dtf...")
SCALE = 'UTC'
NDP = 5
D1 = 2400000.5
D2 = 49533.99999
print(sf.pymD2dtf(SCALE, NDP, D1, D2))

line()
print("Testing pymD2tf...")
NDP = 4
DAYS = -0.987654321
print(sf.pymD2tf(NDP, DAYS))

line()
print("Testing pymDat...")
IY = 2017
IM  = 9
ID = 1
FD = 0.0
print(sf.pymDat(IY, IM, ID, FD))

line()
print("Testing pymDtdb...")
DATE1 = 2448939.5
DATE2 = 0.123
UT = 0.76543
ELONG = 5.0123
U = 5525.242
V = 3190.0
print(sf.pymDtdb(DATE1, DATE2, UT, ELONG, U, V))

line()
print("Testing pymDtf2d...")
SCALE = 'UTC'
IY = 1994
IM  = 6
ID = 30
IHR = 23
IMN = 59
SEC = 60.13599
print(sf.pymDtf2d(SCALE, IY, IM, ID, IHR, IMN, SEC))

line()
print("Testing pymEceq06...")
DATE1 = 2456165.5
DATE2 = 0.401182685
DL = 5.1
DB = -0.9
print(sf.pymEceq06(DATE1, DATE2, DL, DB))

line()
print("Testing pymEcm06...")
DATE1 = 2456165.5
DATE2 = 0.401182685
print(sf.pymEcm06(DATE1, DATE2))

line()
print("Testing pymEe00...")
DATE1 = 2400000.5
DATE2 = 53736.0
EPSA = 0.4090789763356509900
DPSI = -0.9630909107115582393e-5
print(sf.pymEe00(DATE1, DATE2, EPSA, DPSI))

line()
print("Testing pymEe00a...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymEe00a(DATE1, DATE2))

line()
print("Testing pymEe00b...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymEe00b(DATE1, DATE2))

line()
print("Testing pymEe06a...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymEe06a(DATE1, DATE2))

line()
print("Testing pymEect00...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymEect00(DATE1, DATE2))

line()
print("Testing pymEform...")
N = 1
print(sf.pymEform(N))

line()
print("Testing pymEo06a...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymEo06a(DATE1, DATE2))

line()
print("Testing pymEors...")
RNPB = [[0.9999989440476103608, -0.1332881761240011518e-2, -0.5790767434730085097e-3],
        [0.1332858254308954453e-2, 0.9999991109044505944, -0.4097782710401555759e-4],
        [0.5791308472168153320e-3, 0.4020595661593994396e-4, 0.9999998314954572365]]
S = -0.1220040848472271978e-7
print(sf.pymEors(RNPB, S))

line()
print("Testing pymEpb...")
DJ1 = 2415019.8135
DJ2 = 30103.18648
print(sf.pymEpb(DJ1, DJ2))

line()
print("Testing pymEpb2jd...")
EPB = 1957.3
print(sf.pymEpb2jd(EPB))

line()
print("Testing pymEpj...")
DJ1 = 2451545
DJ2 = -7392.5
print(sf.pymEpj(DJ1, DJ2))

line()
print("Testing pymEpj2jd...")
EPJ = 1996.8
print(sf.pymEpj2jd(EPJ))

line()
print("Testing pymEpv00...")
DATE1 = 2400000.5
DATE2 = 53411.52501161
print(sf.pymEpv00(DATE1, DATE2))

line()
print("Testing pymEqec06...")
DATE1 = 1234.5
DATE2 = 2440000.5
DR = 1.234
DD = 0.987
print(sf.pymEqec06(DATE1, DATE2, DR, DD))

line()
print("Testing pymEqeq94...")
DATE1 = 2400000.5
DATE2 = 41234.0
print(sf.pymEqeq94(DATE1, DATE2))

line()
print("Testing pymEra00...")
DJ1 = 2400000.5
DJ2 = 54388.0
print(sf.pymEra00(DJ1, DJ2))

line()
print("Testing pymFad03...")
T = 0.80
print(sf.pymFad03(T))

line()
print("Testing pymFae03...")
T = 0.80
print(sf.pymFae03(T))

line()
print("Testing pymFaf03...")
T = 0.80
print(sf.pymFaf03(T))

line()
print("Testing pymFaju03...")
T = 0.80
print(sf.pymFaju03(T))

line()
print("Testing pymFal03...")
T = 0.80
print(sf.pymFal03(T))

line()
print("Testing pymFalp03...")
T = 0.80
print(sf.pymFalp03(T))

line()
print("Testing pymFama03...")
T = 0.80
print(sf.pymFama03(T))

line()
print("Testing pymFame03...")
T = 0.80
print(sf.pymFame03(T))

line()
print("Testing pymFane03...")
T = 0.80
print(sf.pymFane03(T))

line()
print("Testing pymFaom03...")
T = 0.80
print(sf.pymFaom03(T))

line()
print("Testing pymFapa03...")
T = 0.80
print(sf.pymFapa03(T))

line()
print("Testing pymFasa03...")
T = 0.80
print(sf.pymFasa03(T))

line()
print("Testing pymFaur03...")
T = 0.80
print(sf.pymFaur03(T))

line()
print("Testing pymFave03...")
T = 0.80
print(sf.pymFave03(T))

line()
print("Testing pymFk425...")
R1950 = 0.07626899753879587532
D1950 = -1.137405378399605780
DR1950 = 0.1973749217849087460e-4
DD1950 = 0.5659714913272723189e-5
P1950 = 0.134
V1950 = 8.7
print(sf.pymFk425(R1950, D1950, DR1950, DD1950, P1950, V1950))

line()
print("Testing pymFk45z...")
R2000 = 0.01602284975382960982
D2000 = -0.1164347929099906024
BEPOCH = 1954.677617625256806
print(sf.pymFk45z(R2000, D2000, BEPOCH))

line()
print("Testing pymFk524...")
R2000 = 0.8723503576487275595
D2000 = -0.7517076365138887672
DR2000 = 0.2019447755430472323e-4
DD2000 = 0.3541563940505160433e-5
P2000 = 0.1559
V2000 = 86.87
print(sf.pymFk524(R2000, D2000, DR2000, DD2000, P2000, V2000))

line()
print("Testing pymFk52h...")
R5  =  1.76779433
D5  = -0.2917517103
DR5 = -1.91851572e-7
DD5 = -5.8468475e-6
PX5 =  0.379210
RV5 = -7.6
print(sf.pymFk52h(R5, D5, DR5, DD5, PX5, RV5))

line()
print("Testing pymFk54z...")
R2000 = 0.02719026625066316119
D2000 = -0.1115815170738754813
BEPOCH = 1954.677308160316374
print(sf.pymFk54z(R2000, D2000, BEPOCH))

line()
print("Testing pymFk5hip...")
print(sf.pymFk5hip())

line()
print("Testing pymFk5hz...")
R5 =  1.76779433
D5 = -0.2917517103
DATE1 = 2400000.5
DATE2 = 54479.0
print(sf.pymFk5hz(R5, D5, DATE1, DATE2))

line()
print("Testing pymFw2m...")
GAMB = -0.2243387670997992368e-5
PHIB =  0.4091014602391312982
PSI  = -0.9501954178013015092e-3
EPS  =  0.4091014316587367472
print(sf.pymFw2m(GAMB, PHIB, PSI, EPS))

line()
print("Testing pymFw2xy...")
GAMB = -0.2243387670997992368e-5
PHIB =  0.4091014602391312982
PSI  = -0.9501954178013015092e-3
EPS  =  0.4091014316587367472
print(sf.pymFw2xy(GAMB, PHIB, PSI, EPS))

line()
print("Testing pymG2icrs...")
DL =  5.5850536063818546461558105
DB = -0.7853981633974483096156608
print(sf.pymG2icrs(DL, DB))

line()
print("Testing pymGC2GD...")
N = 2
XYZ = [2e6, 3e6, 5.244e6]
print(sf.pymGC2GD(N, XYZ))

line()
print("Testing pymGc2gde...")
A = 6378136.0
F = 0.0033528
XYZ = [2e6, 3e6, 5.244e6]
print(sf.pymGc2gde(A, F, XYZ))

line()
print("Testing pymGd2gc...")
N = 1
ELONG = 3.1
PHI = -0.5
HEIGHT = 2500.0
print(sf.pymGd2gc(N, ELONG, PHI, HEIGHT))

line()
print("Testing pymGd2gce...")
A = 6378136.0
F = 0.0033528
ELONG = 3.1
PHI = -0.5
HEIGHT = 2500.0
print(sf.pymGd2gce(A, F, ELONG, PHI, HEIGHT))

line()
print("Testing pymGmst00...")
UTA = 2400000.5
UTB = 53736.0
TTA = 2400000.5
TTB = 53736.0
print(sf.pymGmst00(UTA, UTB, TTA, TTB))

line()
print("Testing pymGmst06...")
UTA = 2400000.5
UTB = 53736.0
TTA = 2400000.5
TTB = 53736.0
print(sf.pymGmst06(UTA, UTB, TTA, TTB))

line()
print("Testing pymGmst82...")
DJ1 = 2400000.5
DJ2 = 53736.0
print(sf.pymGmst82(DJ1, DJ2))

line()
print("Testing pymGst00a...")
TTA = 2400000.5
UTA = 2400000.5
TTB = 53736.0
UTB = 53736.0
print(sf.pymGst00a(UTA, UTB, TTA, TTB))

line()
print("Testing pymGst00b...")
UTA = 2400000.5
UTB = 53736.0
print(sf.pymGst00b(UTA, UTB))

line()
print("Testing pymGst06...")
TTA = 2400000.5
UTA = 2400000.5
TTB = 53736.0
UTB = 53736.0
RNPB = [[0.9999989440476103608, -0.1332881761240011518e-2, -0.5790767434730085097e-3],
        [0.1332858254308954453e-2, 0.9999991109044505944, -0.4097782710401555759e-4],
        [0.5791308472168153320e-3, 0.4020595661593994396e-4, 0.9999998314954572365]]
print(sf.pymGst06(UTA, UTB, TTA, TTB, RNPB))

line()
print("Testing pymGst06a...")
TTA = 2400000.5
UTA = 2400000.5
TTB = 53736.0
UTB = 53736.0
print(sf.pymGst06a(UTA, UTB, TTA, TTB))

line()
print("Testing pymGst94...")
UTA = 2400000.5
UTB = 53736.0
print(sf.pymGst94(UTA, UTB))

line()
print("Testing pymIcrs2g...")
DR =  5.9338074302227188048671087
DD = -1.1784870613579944551540570
print(sf.pymIcrs2g(DR, DD))

line()
print("Testing pymH2fk5...")
RH  =  1.767794352
DH  = -0.2917512594
DRH = -2.76413026e-6
DDH = -5.92994449e-6
PXH =  0.379210
RVH = -7.6
print(sf.pymH2fk5(RH, DH, DRH, DDH, PXH, RVH))

line()
print("Testing pymHd2ae...")
HA = 1.1
DEC = 1.2
PHI = 0.3
print(sf.pymHd2ae(HA, DEC, PHI))

line()
print("Testing pymHd2pa...")
HA = 1.1
DEC = 1.2
PHI = 0.3
print(sf.pymHd2pa(HA, DEC, PHI))

line()
print("Testing pymHfk5z...")
RH = 1.767794352
DH = -0.2917512594
DATE1 = 2400000.5
DATE2 = 54479.0
print(sf.pymHfk5z(RH, DH, DATE1, DATE2))

line()
print("Testing pymIr...")
print(sf.pymIr())

line()
print("Testing pymJd2cal...")
DJ1 = 2400000.5
DJ2 = 50123.9999
print(sf.pymJd2cal(DJ1, DJ2))

line()
print("Testing pymJdcalf...")
N=4
DJ1 = 2400000.5
DJ2 = 50123.9999
print(sf.pymJdcalf(NDP, DJ1, DJ2))

line()
print("Testing pymLd...")
BM = 0.00028574
P = [-0.763276255, -0.608633767, -0.216735543]
Q = [-0.763276255, -0.608633767, -0.216735543]
E = [0.76700421, 0.605629598, 0.211937094]
EM = 8.91276983
DLIM = 3e-10
print(sf.pymLd(BM, P, Q, E, EM, DLIM))

line()
print("Testing pymLdn...")
N = 3
B = [[0.00028574, 3e-10, -7.81014427, -5.60956681,
      -1.98079819, 0.0030723249, -0.00406995477, -0.00181335842],
     [0.00095435, 3e-9, 0.738098796, 4.63658692, 
      1.9693136, -0.00755816922, 0.00126913722, 0.000727999001],
     [1.0, 6e-6, -0.000712174377, -0.00230478303, 
      -0.00105865966, 6.29235213e-6, -3.30888387e-7, -2.96486623e-7]]
OB = [-0.974170437, -0.2115201, -0.0917583114]
SC = [-0.763276255, -0.608633767, -0.216735543]
print(sf.pymLdn(N, B, OB, SC))

line()
print("Testing pymLdsun...")
P = [-0.763276255, -0.608633767, -0.216735543]
E = [-0.973644023, -0.20925523, -0.0907169552]
EM = 0.999809214
print(sf.pymLdsun(P, E, EM))
#####################################
line()
print("Testing pymLteceq...")
EPJ = 2500.0
DL = 1.5
DB = 0.6
print(sf.pymLteceq(EPJ, DL, DB))

line()
print("Testing pymLtecm...")
EPJ = -3000.0
print(sf.pymLtecm(EPJ))

line()
print("Testing pymLteqec...")
EPJ = -1500
DR = 1.234
DD = 0.987
print(sf.pymLteqec(EPJ, DR, DD))

line()
print("Testing pymLtp...")
EPJ = 1666.666
print(sf.pymLtp(EPJ))

line()
print("Testing pymLtpb...")
EPJ = 1666.666
print(sf.pymLtpb(EPJ))

line()
print("Testing pymLtpecl...")
EPJ = -1500
print(sf.pymLtpecl(EPJ))

line()
print("Testing pymLtpequ...")
EPJ = -2500
print(sf.pymLtpequ(EPJ))

line()
print("Testing pymMoon98...")
DATE1 = 2400000.5
DATE2 = 43999.9
print(sf.pymMoon98(DATE1, DATE2))

line()
print("Testing pymNum00a...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymNum00a(DATE1, DATE2))

line()
print("Testing pymNum00b...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymNum00b(DATE1, DATE2))

line()
print("Testing pymNum06a...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymNum06a(DATE1, DATE2))

line()
print("Testing pymNumat...")
EPSA =  0.4090789763356509900
DPSI = -0.9630909107115582393e-5
DEPS =  0.4063239174001678826e-4
print(sf.pymNumat(EPSA, DPSI, DEPS))

line()
print("Testing pymNut00a...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymNut00a(DATE1, DATE2))

line()
print("Testing pymNut00b...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymNut00b(DATE1, DATE2))

line()
print("Testing pymNut06a...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymNut06a(DATE1, DATE2))

line()
print("Testing pymNut80...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymNut80(DATE1, DATE2))

line()
print("Testing pymNutm80...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymNutm80(DATE1, DATE2))

line()
print("Testing pymObl06...")
DATE1 = 2400000.5
DATE2 = 54388.0
print(sf.pymObl06(DATE1, DATE2))

line()
print("Testing pymObl80...")
DATE1 = 2400000.5
DATE2 = 54388.0
print(sf.pymObl80(DATE1, DATE2))

line()
print("Testing pymP06e...")
DATE1 = 2400000.5
DATE2 = 52541.0
print(sf.pymP06e(DATE1, DATE2))

line()
print("Testing pymP2pv...")
P = [0.25, 1.2, 3.0]
print(sf.pymP2pv(P))

line()
print("Testing pymP2s...")
P = [100.0, -50.0, 25.0]
print(sf.pymP2s(P))

line()
print("Testing pymPap...")
A = [1.0, 0.1, 0.2]
B = [-3.0, 1E-3, 0.2]
print(sf.pymPap(A, B))

line()
print("Testing pymPas...")
AL = 1.0
AP = 0.1
BL = 0.2
BP = -1.0
print(sf.pymPas(AL, AP, BL, BP))

line()
print("Testing pymPb06...")
DATE1 = 2400000.5
DATE2 = 50123.9999
print(sf.pymPb06(DATE1, DATE2))

line()
print("Testing pymPdp...")
A = [2.0, 2.0, 3.0]
B = [1.0, 3.0, 4.0]
print(sf.pymPdp(A, B))

line()
print("Testing pymPfw06...")
DATE1 = 2400000.5
DATE2 = 50123.9999
print(sf.pymPfw06(DATE1, DATE2))

line()
print("Testing pymPlan94...")
DATE1 = 2400000.5
DATE2 = 1E6
NP = 0
print(sf.pymPlan94(DATE1, DATE2, NP))
DATE1 = 2400000.5
DATE2 = 1E6
NP = 10
print(sf.pymPlan94(DATE1, DATE2, NP))
DATE1 = 2400000.5
DATE2 = -320000
NP = 3
print(sf.pymPlan94(DATE1, DATE2, NP))
DATE1 = 2400000.5
DATE2 = 43999.9
NP = 1
print(sf.pymPlan94(DATE1, DATE2, NP))

line()
print("Testing pymPmat00...")
DATE1 = 2400000.5
DATE2 = 50123.9999
print(sf.pymPmat00(DATE1, DATE2))

line()
print("Testing pymPmat06...")
DATE1 = 2400000.5
DATE2 = 50123.9999
print(sf.pymPmat06(DATE1, DATE2))

line()
print("Testing pymPmat76...")
DATE1 = 2400000.5
DATE2 = 50123.9999
print(sf.pymPmat76(DATE1, DATE2))

line()
print("Testing pymPm...")
P = [0.3, 1.2, -2.5]
print(sf.pymPm(P))

line()
print("Testing pymPmp...")
A = [2.0, 2.0, 3.0]
B = [1.0, 3.0, 4.0]
print(sf.pymPmp(A, B))

line()
print("Testing pymPmpx...")
RC = 1.234
DC = 0.789
PR = 1e-5
PD = -2e-5
PX = 1e-2
RV = 10.0
PMT = 8.75
POB = [0.9, 0.4, 0.1]
print(sf.pymPmpx(RC, DC, PR, PD, PX, RV, PMT, POB))

line()
print("Testing pymPmsafe...")
RA1 = 1.234
DEC1 = 0.789
PMR1 = 1e-5
PMD1 = -2e-5
PX1 = 1e-2
RV1 = 10.0
EP1A = 2400000.5
EP1B = 48348.5625
EP2A = 2400000.5
EP2B = 51544.5
print(sf.pymPmsafe(RA1, DEC1, PMR1, PMD1, PX1, RV1, EP1A, EP1B, EP2A, EP2B))

line()
print("Testing pymPn...")
P = [0.3, 1.2, -2.5]
print(sf.pymPn(P))

line()
print("Testing pymPn00...")
DATE1 = 2400000.5
DATE2 = 53736.0
DPSI = -0.9632552291149335877e-5
DEPS =  0.4063197106621141414e-4
print(sf.pymPn00(DATE1, DATE2, DPSI, DEPS))

line()
print("Testing pymPn00a...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymPn00a(DATE1, DATE2))

line()
print("Testing pymPn00b...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymPn00b(DATE1, DATE2))

line()
print("Testing pymPn06a...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymPn06a(DATE1, DATE2))

line()
print("Testing pymPn06...")
DATE1 = 2400000.5
DATE2 = 53736.0
DPSI = -0.9632552291149335877e-5
DEPS =  0.4063197106621141414e-4
print(sf.pymPn06(DATE1, DATE2, DPSI, DEPS))

line()
print("Testing pymPnm00a...")
DATE1 = 2400000.5
DATE2 = 50123.9999
print(sf.pymPnm00a(DATE1, DATE2))

line()
print("Testing pymPnm00b...")
DATE1 = 2400000.5
DATE2 = 50123.9999
print(sf.pymPnm00b(DATE1, DATE2))

line()
print("Testing pymPnm06a...")
DATE1 = 2400000.5
DATE2 = 50123.9999
print(sf.pymPnm06a(DATE1, DATE2))

line()
print("Testing pymPnm80...")
DATE1 = 2400000.5
DATE2 = 50123.9999
print(sf.pymPnm80(DATE1, DATE2))

line()
print("Testing pymPom00...")
XP = 2.55060238e-7
YP = 1.860359247e-6
SP = -0.1367174580728891460e-10
print(sf.pymPom00(XP, YP, SP))

line()
print("Testing pymPpp...")
A = [2.0, 2.0, 3.0]
B = [1.0, 3.0, 4.0]
print(sf.pymPpp(A, B))

line()
print("Testing pymPpsp...")
A = [2.0, 2.0, 3.0]
S = 5.0
B = [1.0, 3.0, 4.0]
print(sf.pymPpsp(A, S, B))

line()
print("Testing pymPr00...")
DATE1 = 2400000.5
DATE2 = 53736
print(sf.pymPr00(DATE1, DATE2))

line()
print("Testing pymPrec76...")
DATE01 = 2400000.5
DATE02 = 33282.0
DATE11 = 2400000.5
DATE12 = 51544.0
print(sf.pymPrec76(DATE01, DATE02, DATE11, DATE12))

line()
print("Testing pymPv2p...")
PV = [[0.3, 1.2, -2.5],
      [-0.5, 3.1, 0.9]]
print(sf.pymPv2p(PV))

line()
print("Testing pymPv2s...")
PV = [[-0.4514964673880165, 0.03093394277342585, 0.05594668105108779],
      [1.292270850663260e-5, 2.652814182060692e-6, 2.568431853930293e-6]]
print(sf.pymPv2s(PV))

line()
print("Testing pymPvdpv...")
A = [[2.0, 2.0, 3.0],
      [6.0, 0.0, 4.0]]
B = [[1.0, 3.0, 4.0],
      [0.0, 2.0, 8.0]]
print(sf.pymPvdpv(A, B))

line()
print("Testing pymPvm...")
PV = [[0.3, 1.2, -2.5],
      [0.45, -0.25, 1.1]]
print(sf.pymPvm(PV))

line()
print("Testing pymPvmpv...")
A = [[2.0, 2.0, 3.0],
      [5.0, 6.0, 3.0]]
B = [[1.0, 3.0, 4.0],
      [3.0, 2.0, 1.0]]
print(sf.pymPvmpv(A, B))

line()
print("Testing pymPvppv...")
A = [[2.0, 2.0, 3.0],
      [5.0, 6.0, 3.0]]
B = [[1.0, 3.0, 4.0],
      [3.0, 2.0, 1.0]]
print(sf.pymPvppv(A, B))

line()
print("Testing pymPvstar...")
PV = [[126668.5912743160601, 2136.792716839935195, -245251.2339876830091],
      [-0.4051854035740712739e-2, -0.6253919754866173866e-2, 0.1189353719774107189e-1]]
print(sf.pymPvstar(PV))

line()
print("Testing pymPvtob...")
ELONG = 2.0
PHI = 0.5
HM = 3000.0
XP = 1e-6
YP = -0.5e-6
SP = 1e-8
THETA = 5.0
print(sf.pymPvtob(ELONG, PHI, HM, XP, YP, SP, THETA))

line()
print("Testing pymPvu...")
DT = 2920.0
PV = [[126668.5912743160601, 2136.792716839935195, -245251.2339876830091],
      [-0.4051854035740712739e-2, -0.6253919754866173866e-2, 0.1189353719774107189e-1]]
print(sf.pymPvu(DT, PV))

line()
print("Testing pymPvup...")
DT = 2920.0
PV = [[126668.5912743160601, 2136.792716839935195, -245251.2339876830091],
      [-0.4051854035740712739e-2, -0.6253919754866173866e-2, 0.1189353719774107189e-1]]
print(sf.pymPvup(DT, PV))

line()
print("Testing pymPvxpv...")
A = [[2.0, 2.0, 3.0],
      [6.0, 0.0, 4.0]]
B = [[1.0, 3.0, 4.0],
      [0.0, 2.0, 8.0]]
print(sf.pymPvxpv(A, B))

line()
print("Testing pymPxp...")
A = [2.0, 2.0, 3.0]
B = [1.0, 3.0, 4.0]
print(sf.pymPxp(A, B))

line()
print("Testing pymRefco...")
PHPA = 800.0
TC = 10.0
RH = 0.9
WL = 0.4
print(sf.pymRefco(PHPA, TC, RH, WL))

line()
print("Testing pymRm2v...")
R = [[0.00, -0.80, -0.60],
     [0.80, -0.36, 0.48],
     [0.60, 0.48, -0.64]]
print(sf.pymRm2v(R))

line()
print("Testing pymRv2m...")
W = [0.0, 1.41371669, -1.88495559]
print(sf.pymRv2m(W))

line()
print("Testing pymRx...")
PHI = 0.3456789
R = [[2.0, 3.0, 2.0],
      [3.0, 2.0, 3.0],
      [3.0, 4.0, 5.0]]
print(sf.pymRx(PHI, R))

line()
print("Testing pymRxp...")
R = [[2.0, 3.0, 2.0],
      [3.0, 2.0, 3.0],
      [3.0, 4.0, 5.0]]
P = [0.2, 1.5, 0.1]
print(sf.pymRxp(R, P))

line()
print("Testing pymRxpv...")
R = [[2.0, 3.0, 2.0],
      [3.0, 2.0, 3.0],
      [3.0, 4.0, 5.0]]
PV = [[0.2, 1.5, 0.1],
      [1.5, 0.2, 0.1]]
print(sf.pymRxpv(R, PV))

line()
print("Testing pymRxr...")
A = [[2.0, 3.0, 2.0],
      [3.0, 2.0, 3.0],
      [3.0, 4.0, 5.0]]
B = [[1.0, 2.0, 2.0],
      [4.0, 1.0, 1.0],
      [3.0, 0.0, 1.0]]
print(sf.pymRxr(A, B))

line()
print("Testing pymRy...")
THETA = 0.3456789
R = [[2.0, 3.0, 2.0],
      [3.0, 2.0, 3.0],
      [3.0, 4.0, 5.0]]
print(sf.pymRy(THETA, R))

line()
print("Testing pymRz...")
PSI = 0.3456789
R = [[2.0, 3.0, 2.0],
      [3.0, 2.0, 3.0],
      [3.0, 4.0, 5.0]]
print(sf.pymRz(PSI, R))

line()
print("Testing pymS00a...")
DATE1 = 2400000.5
DATE2 = 52541.0
print(sf.pymS00a(DATE1, DATE2))

line()
print("Testing pymS00b...")
DATE1 = 2400000.5
DATE2 = 52541.0
print(sf.pymS00b(DATE1, DATE2))

line()
print("Testing pymS00...")
DATE1 = 2400000.5
DATE2 = 53736.0
X = 0.5791308486706011000e-3
Y = 0.4020579816732961219e-4
print(sf.pymS00(DATE1, DATE2, X, Y))

line()
print("Testing pymS06a...")
DATE1 = 2400000.5
DATE2 = 52541.0
print(sf.pymS06a(DATE1, DATE2))

line()
print("Testing pymS06...")
DATE1 = 2400000.5
DATE2 = 53736.0
X = 0.5791308486706011000e-3
Y = 0.4020579816732961219e-4
print(sf.pymS06(DATE1, DATE2, X, Y))

line()
print("Testing pymS2c...")
THETA = 3.0123
PHI = -0.999
print(sf.pymS2c(THETA, PHI))

line()
print("Testing pymS2p...")
THETA = -3.21
PHI = 0.123
R = 0.456
print(sf.pymS2p(THETA, PHI, R))

line()
print("Testing pymS2pv...")
THETA = -3.21
PHI = 0.123
R = 0.456
TD = -7.8e-6
PD = 9.01e-6
RD = -1.23e-5
print(sf.pymS2pv(THETA, PHI, R, TD, PD, RD))

line()
print("Testing pymS2xpv...")
S1 = 2.0
S2 = 3.0
PV = [[0.3, 1.2, -2.5],
      [0.5, 2.3, -0.4]]
print(sf.pymS2xpv(S1, S2, PV))

line()
print("Testing pymSepp...")
A = [1.0, 0.1, 0.2]
B = [-3.0, 1E-3, 0.2]
print(sf.pymSepp(A, B))

line()
print("Testing pymSeps...")
AL =  1.0
AP =  0.1
BL =  0.2
BP = -3.0
print(sf.pymSeps(AL, AP, BL, BP))

line()
print("Testing pymSp00...")
DATE1 = 2400000.5
DATE2 = 52541.0
print(sf.pymSp00(DATE1, DATE2))

line()
print("Testing pymStarpm...")
RA1 =   0.01686756
DEC1 = -1.093989828
PMR1 = -1.78323516e-5
PMD1 =  2.336024047e-6
PX1 =   0.74723
RV1 = -21.6
EP1A = 2400000.5
EP1B = 50083.0
EP2A = 2400000.5
EP2B = 53736.0
print(sf.pymStarpm(RA1, DEC1, PMR1, PMD1, PX1, RV1, EP1A, EP1B, EP2A, EP2B))

line()
print("Testing pymStarpv...")
RA =   0.01686756
DEC = -1.093989828
PMR = -1.78323516e-5
PMD =  2.336024047e-6
PX =   0.74723
RV = -21.6
EP1A = 2400000.5
EP1B = 50083.0
EP2A = 2400000.5
EP2B = 53736.0
print(sf.pymStarpv(RA, DEC, PMR, PMD, PX, RV))

line()
print("Testing pymSxp...")
S = 2.0
P = [0.3, 1.2, -2.5]
print(sf.pymSxp(S, P))

line()
print("Testing pymSxpv...")
S = 2.0
PV = [[0.3, 1.2, -2.5],
      [0.5, 3.2, -0.7]]
print(sf.pymSxpv(S, PV))

line()
print("Testing pymTaitt...")
TAI1 = 2453750.5
TAI2 = 0.892482639
print(sf.pymTaitt(TAI1, TAI2))

line()
print("Testing pymTaiut1...")
TAI1 = 2453750.5
TAI2 = 0.892482639
DTA = -32.6659
print(sf.pymTaiut1(TAI1, TAI2, DTA))

line()
print("Testing pymTaiutc...")
TAI1 = 2453750.5
TAI2 = 0.892482639
print(sf.pymTaiutc(TAI1, TAI2))

line()
print("Testing pymTcbtdb...")
TCB1 = 2453750.5
TCB2 = 0.893019599
print(sf.pymTcbtdb(TCB1, TCB2))

line()
print("Testing pymTcgtt...")
TCG1 = 2453750.5
TCG2 = 0.892862531
print(sf.pymTcgtt(TCG1, TCG2))

line()
print("Testing pymTdbtcb...")
TDB1 = 2453750.5
TDB2 = 0.892855137
print(sf.pymTdbtcb(TDB1, TDB2))

line()
print("Testing pymTdbtt...")
TDB1 = 2453750.5
TDB2 = 0.892855137
DTR = -0.000201
print(sf.pymTdbtt(TDB1, TDB2, DTR))

line()
print("Testing pymTf2a...")
S = '+'
IHOUR = 4
IMIN = 58
SEC = 20.2
print(sf.pymTf2a(S, IHOUR, IMIN, SEC))

line()
print("Testing pymTf2d...")
S = ' '
IHOUR = 23
IMIN = 55
SEC = 10.9
print(sf.pymTf2d(S, IHOUR, IMIN, SEC))

line()
print("Testing pymTpors...")
XI = -0.03
ETA = 0.07
A = 1.3
B = 1.5
print(sf.pymTpors(XI, ETA, A, B))

line()
print("Testing pymTporv...")
XI = -0.03
ETA = 0.07
V = sf.pymS2c(1.3, 1.5)
print(sf.pymTporv(XI, ETA, V))

line()
print("Testing pymTpsts...")
XI = -0.03
ETA = 0.07
A0 = 2.3
B0 = 1.5
print(sf.pymTpsts(XI, ETA, A0, B0))

line()
print("Testing pymTpstv...")
XI = -0.03
ETA = 0.07
A0 = 2.3
B0 = 1.5
V0 = sf.pymS2c(A0, B0)
print(sf.pymTpstv(XI, ETA, V0))

line()
print("Testing pymTpxes...")
A = 1.3
B = 1.55
A0 = 2.3
B0 = 1.5
print(sf.pymTpxes(A, B, A0, B0))

line()
print("Testing pymTpxev...")
A = 1.3
B = 1.55
A0 = 2.3
B0 = 1.5
V = sf.pymS2c(A, B)
V0 = sf.pymS2c(A0, B0)
print(sf.pymTpxev(V, V0))

line()
print("Testing pymTr...")
R = [[2.0, 3.0, 2.0],
      [3.0, 2.0, 3.0],
      [3.0, 4.0, 5.0]]
print(sf.pymTr(R))

line()
print("Testing pymTrxp...")
R = [[2.0, 3.0, 2.0],
      [3.0, 2.0, 3.0],
      [3.0, 4.0, 5.0]]
P = [0.2, 1.5, 0.1]
print(sf.pymTrxp(R, P))

line()
print("Testing pymTrxpv...")
R = [[2.0, 3.0, 2.0],
      [3.0, 2.0, 3.0],
      [3.0, 4.0, 5.0]]
PV = [[0.2, 1.5, 0.1],
      [1.5, 0.2, 0.1]]
print(sf.pymTrxpv(R, PV))

line()
print("Testing pymTttai...")
TT1 = 2453750.5
TT2 = 0.892482639
print(sf.pymTttai(TT1, TT2))

line()
print("Testing pymTttcg...")
TT1 = 2453750.5
TT2 = 0.892482639
print(sf.pymTttcg(TT1, TT2))

line()
print("Testing pymTttdb...")
TT1 = 2453750.5
TT2 = 0.892855139
DTR = -0.000201
print(sf.pymTttdb(TT1, TT2, DTR))

line()
print("Testing pymTtut1...")
TT1 = 2453750.5
TT2 = 0.892855139
DT = 64.8499
print(sf.pymTtut1(TT1, TT2, DT))

line()
print("Testing pymUt1tai...")
UT11 = 2453750.5
UT12 = 0.892104561
DTA = -32.6659
print(sf.pymUt1tai(UT11, UT12, DTA))

line()
print("Testing pymUt1tt...")
UT11 = 2453750.5
UT12 = 0.892104561
DTA = 64.8499
print(sf.pymUt1tt(UT11, UT12, DT))

line()
print("Testing pymUt1utc...")
UT11 = 2453750.5
UT12 = 0.892104561
DUT1 = 0.3341
print(sf.pymUt1utc(UT11, UT12, DUT1))

line()
print("Testing pymUtctai...")
UTC1 = 2453750.5
UTC2 = 0.892100694
print(sf.pymUtctai(UTC1, UTC2))

line()
print("Testing pymUtcut1...")
UTC1 = 2453750.5
UTC2 = 0.892100694
DUT1 = 0.3341
print(sf.pymUtcut1(UTC1, UTC2, DUT1))

line()
print("Testing pymXy06...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymXy06(DATE1, DATE2))

line()
print("Testing pymXys00a...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymXys00a(DATE1, DATE2))

line()
print("Testing pymXys00b...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymXys00b(DATE1, DATE2))

line()
print("Testing pymXys06a...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymXys06a(DATE1, DATE2))

line()
print("Testing pymZp..")
print(sf.pymZp())

line()
print("Testing pymZpv..")
print(sf.pymZpv())

line()
print("Testing pymZr..")
print(sf.pymZr())
