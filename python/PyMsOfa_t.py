import numpy  as  np
import numpy.ctypeslib  as  nt
import PyMsOfa  as  sf
import time

def   line():
        print('-'*60)

line()
print("Testing pymA2AF...")
NDP = 4
ANGLE = 2.345
S, IDMSF = sf.pymA2AF(NDP, ANGLE)
print(S, IDMSF)

line()
print("Testing pymA2TF...")
NDP=4
ANGLE=3.01234
S, IHMSF = sf.pymA2TF(NDP, ANGLE)
print(S, IHMSF)

line()
print("Testing pymAB...")
PNAT = [-0.76321968546737951, -0.60869453983060384, -0.21676408580639883]
V = [2.1044018893653786E-5, -8.9108923304429319E-5, -3.8633714797716569E-5]
S = 0.99980921395708788
BM1 = 0.99999999506209258
PPR = sf.pymAB(PNAT, V, S, BM1)
print(PPR)

line()
print("Testing pymAE2HD...")
a = 5.5
e = 1.1
p = 0.7
h, d = sf.pymAE2HD(a, e, p)
print(h, d)

line()
print("Testing pymAF2A...")
S = '-'
IDEG = 45
IAMIN = 13
ASEC = 27.2
RAD, J = sf.pymAF2A(S, IDEG, IAMIN, ASEC)
print(RAD, J)

line()
print("Testing pymANP...")
A=-0.1
print(sf.pymANP(A))

line()
print("Testing pymANPM...")
A=-4.0
print(sf.pymANPM(A))

line()
print("Testing pymAPCG...")
ASTROM=[0 for i in range(30)]
DATE1 = 2456165.5
DATE2 = 0.401182685
EBPV = [[0.901310875, -0.417402664, -0.180982288],
        [0.00742727954, 0.0140507459, 0.00609045792]]
EHP = [0.903358544, -0.415395237, -0.180084014]
print(sf.pymAPCG(DATE1, DATE2, EBPV, EHP, ASTROM))

line()
print("Testing pymAPCG13...")
ASTROM=[0 for i in range(30)]
DATE1 = 2456165.5
DATE2 = 0.401182685
print(sf.pymAPCG13(DATE1, DATE2, ASTROM))

line()
print("Testing pymAPCI...")
ASTROM=[0 for i in range(30)]
DATE1 = 2456165.5
DATE2 = 0.401182685
EBPV = [[0.901310875, -0.417402664, -0.180982288],
        [0.00742727954, 0.0140507459, 0.00609045792]]
EHP = [0.903358544, -0.415395237, -0.180084014]
X =  0.0013122272
Y = -2.92808623e-5
S =  3.05749468e-8
print(sf.pymAPCI(DATE1, DATE2, EBPV, EHP, X, Y, S, ASTROM))

line()
print("Testing pymAPCI13...")
ASTROM=[0 for i in range(30)]
DATE1 = 2456165.5
DATE2 = 0.401182685
print(sf.pymAPCI13(DATE1, DATE2, ASTROM))

line()
print("Testing pymAPCO...")
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
print(sf.pymAPCO(DATE1, DATE2, EBPV, EHP, X, Y, S, THETA, ELONG, 
                  PHI, HM, XP, YP, SP, REFA, REFB))

line()
print("Testing pymAPCO13...")
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
print(sf.pymAPCO13(UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA,
                    TC, RH, WL))

line()
print("Testing pymAPCS...")
ASTROM=[0 for i in range(30)]
DATE1 = 2456384.5
DATE2 = 0.970031644
PV = [[-1836024.09, 1056607.72, -5998795.26],
      [-77.0361767, -133.310856, 0.0971855934]]
EBPV = [[-0.974170438, -0.211520082, -0.0917583024],
        [0.00364365824, -0.0154287319, -0.00668922024]]
EHP = [-0.973458265, -0.209215307, -0.0906996477]
print(sf.pymAPCS(DATE1, DATE2, PV, EBPV, EHP, ASTROM))

line()
print("Testing pymAPCS13...")
ASTROM=[0 for i in range(30)]
DATE1 = 2456165.5
DATE2 = 0.401182685
PV = [[-6241497.16, 401346.896, -1251136.04],
      [-29.264597, -455.021831, 0.0266151194]]
print(sf.pymAPCS13(DATE1, DATE2, PV, ASTROM))

line()
print("Testing pymAPER...")
ASTROM=[0 for i in range(30)]
ASTROM[21] = 1.234
THETA = 5.678
print(sf.pymAPER(THETA, ASTROM))

line()
print("Testing pymAPER13...")
ASTROM=[0 for i in range(30)]
ASTROM[21] = 1.234
UT11 = 2456165.5
UT12 = 0.401182685
print(sf.pymAPER13(UT11, UT12, ASTROM))

line()
print("Testing pymAPIO...")
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
print(sf.pymAPIO(SP, THETA, ELONG, PHI, HM, XP, YP, REFA, REFB, ASTROM))

line()
print("Testing pymAPIO13...")
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
print(sf.pymAPIO13(UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL, ASTROM))

line()
print("Testing pymATCC13...")
ASTROM=[0 for i in range(30)]
RC = 2.71
DC = 0.174
PR = 1e-5
PD = 5e-6
PX = 0.1
RV = 55.0
DATE1 = 2456165.5
DATE2 = 0.401182685
print(sf.pymATCC13(RC, DC, PR, PD, PX, RV, DATE1, DATE2))

line()
print("Testing pymATCCQ...")
ASTROM=[0 for i in range(30)]
RC = 2.71
DC = 0.174
PR = 1e-5
PD = 5e-6
PX = 0.1
RV = 55.0
DATE1 = 2456165.5
DATE2 = 0.401182685
ASTROM,EO=sf.pymAPCI13(DATE1, DATE2, ASTROM)
print(sf.pymATCCQ(RC, DC, PR, PD, PX, RV, ASTROM))

line()
print("Testing pymATCI13...")
ASTROM=[0 for i in range(30)]
RC = 2.71
DC = 0.174
PR = 1e-5
PD = 5e-6
PX = 0.1
RV = 55.0
DATE1 = 2456165.5
DATE2 = 0.401182685
print(sf.pymATCI13(RC, DC, PR, PD, PX, RV, DATE1, DATE2))

line()
print("Testing pymATCIQ...")
ASTROM=[0 for i in range(30)]
RC = 2.71
DC = 0.174
PR = 1e-5
PD = 5e-6
PX = 0.1
RV = 55.0
DATE1 = 2456165.5
DATE2 = 0.401182685
ASTROM, EO = sf.pymAPCI13(DATE1, DATE2, ASTROM)
print(sf.pymATCIQ(RC, DC, PR, PD, PX, RV, ASTROM))

line()
print("Testing pymATCIQN...")
ASTROM=[0 for i in range(30)]
RC = 2.71
DC = 0.174
PR = 1e-5
PD = 5e-6
PX = 0.1
RV = 55.0
DATE1 = 2456165.5
DATE2 = 0.401182685
ASTROM, EO = sf.pymAPCI13(DATE1, DATE2, ASTROM)
N = 3
B = [[0.00028574, 3e-10, -7.81014427, -5.60956681,
      -1.98079819, 0.0030723249, -0.00406995477, -0.00181335842],
     [0.00095435, 3e-9, 0.738098796, 4.63658692, 
      1.9693136, -0.00755816922, 0.00126913722, 0.000727999001],
     [1.0, 6e-6, -0.000712174377, -0.00230478303, 
      -0.00105865966, 6.29235213e-6, -3.30888387e-7, -2.96486623e-7]]
print(sf.pymATCIQN(RC, DC, PR, PD, PX, RV, ASTROM, N, B))

line()
print("Testing pymATCIQZ...")
ASTROM=[0 for i in range(30)]
RC = 2.71
DC = 0.174
DATE1 = 2456165.5
DATE2 = 0.401182685
ASTROM, EO = sf.pymAPCI13(DATE1, DATE2, ASTROM)
print(sf.pymATCIQZ(RC, DC, ASTROM))

line()
print("Testing pymATCO13...")
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
print(sf.pymATCO13(RC, DC, PR, PD, PX, RV, UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL))

line()
print("Testing pymATIC13...")
ASTROM=[0 for i in range(30)]
RI = 2.710121572969038991
DI = 0.1729371367218230438
DATE1 = 2456165.5
DATE2 = 0.401182685
print(sf.pymATIC13(RI, DI, DATE1, DATE2))

line()
print("Testing pymATICQ...")
ASTROM=[0 for i in range(30)]
RI = 2.710121572969038991
DI = 0.1729371367218230438
DATE1 = 2456165.5
DATE2 = 0.401182685
ASTROM, EO = sf.pymAPCI13(DATE1, DATE2, ASTROM)
print(sf.pymATICQ(RI, DI, ASTROM))

line()
print("Testing pymATICQN...")
ASTROM=[0 for i in range(30)]
RI = 2.709994899247599271
DI = 0.1728740720983623469
DATE1 = 2456165.5
DATE2 = 0.401182685
ASTROM, EO = sf.pymAPCI13(DATE1, DATE2, ASTROM)
B = [[0.00028574, 3e-10, -7.81014427, -5.60956681,
      -1.98079819, 0.0030723249, -0.00406995477, -0.00181335842],
     [0.00095435, 3e-9, 0.738098796, 4.63658692, 
      1.9693136, -0.00755816922, 0.00126913722, 0.000727999001],
     [1.0, 6e-6, -0.000712174377, -0.00230478303, 
      -0.00105865966, 6.29235213e-6, -3.30888387e-7, -2.96486623e-7]]
print(sf.pymATICQN(RI, DI, ASTROM, N, B))

line()
print("Testing pymATIO13...")
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
print(sf.pymATIO13(RI, DI, UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL))

line()
print("Testing pymATIOQ...")
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
ASTROM, J = sf.pymAPIO13(UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL, ASTROM)
print(sf.pymATIOQ(RI, DI, ASTROM))

line()
print("Testing pymATOC13...")
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
print(sf.pymATOC13(TYPE, OB1, OB2, UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL))
OB1 = -0.09247619879782006106
OB2 = 0.1717653435758265198
TYPE = 'H'
print(sf.pymATOC13(TYPE, OB1, OB2, UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL))
OB1 = 0.09233952224794989993
OB2 = 1.407758704513722461
TYPE = 'A'
print(sf.pymATOC13(TYPE, OB1, OB2, UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL))

line()
print("Testing pymATOI13...")
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
print(sf.pymATOI13(TYPE, OB1, OB2, UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL))
OB1 = -0.09247619879782006106
OB2 = 0.1717653435758265198
TYPE = 'H'
print(sf.pymATOI13(TYPE, OB1, OB2, UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL))
OB1 = 0.09233952224794989993
OB2 = 1.407758704513722461
TYPE = 'A'
print(sf.pymATOI13(TYPE, OB1, OB2, UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL))

line()
print("Testing pymATOIQ...")
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
ASTROM, J = sf.pymAPIO13(UTC1, UTC2, DUT1, ELONG, PHI, HM, XP, YP, PHPA, TC, RH, WL, ASTROM)
OB1 = 2.710085107986886201
OB2 = 0.1717653435758265198
TYPE = 'R'
print(sf.pymATOIQ(TYPE, OB1, OB2, ASTROM))
OB1 = -0.09247619879782006106
OB2 = 0.1717653435758265198
TYPE = 'H'
print(sf.pymATOIQ(TYPE, OB1, OB2, ASTROM))
OB1 = 0.09233952224794989993
OB2 = 1.407758704513722461
TYPE = 'A'
print(sf.pymATOIQ(TYPE, OB1, OB2, ASTROM))

line()
print("Testing pymBI00...")
print(sf.pymBI00())

line()
print("Testing pymBP00...")
DATE1 = 2400000.5
DATE2 = 50123.9999
print(sf.pymBP00(DATE1, DATE2))

line()
print("Testing pymBP06...")
DATE1 = 2400000.5
DATE2 = 50123.9999
print(sf.pymBP06(DATE1, DATE2))

line()
print("Testing pymBPN2XY...")
RBPN = [[9.999962358680738e-1, -2.516417057665452e-3, -1.093569785342370e-3],
        [2.516462370370876e-3, 9.999968329010883e-1, 4.006159587358310e-5],
        [1.093465510215479e-3, -4.281337229063151e-5, 9.999994012499173e-1]]
print(sf.pymBPN2XY(RBPN))

line()
print("Testing pymC2I00A...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymC2I00A(DATE1, DATE2))

line()
print("Testing pymC2I00B...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymC2I00B(DATE1, DATE2))

line()
print("Testing pymC2I06A...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymC2I06A(DATE1, DATE2))

line()
print("Testing pymC2IBPN...")
DATE1 = 2400000.5
DATE2 = 50123.9999
RBPN = [[9.999962358680738e-1, -2.516417057665452e-3, -1.093569785342370e-3],
        [2.516462370370876e-3, 9.999968329010883e-1, 4.006159587358310e-5],
        [1.093465510215479e-3, -4.281337229063151e-5, 9.999994012499173e-1]]
print(sf.pymC2IBPN(DATE1, DATE2, RBPN))

line()
print("Testing pymC2IXY...")
DATE1 = 2400000.5
DATE2 = 53736.0
X = 0.5791308486706011000e-3
Y = 0.4020579816732961219e-4
print(sf.pymC2IXY(DATE1, DATE2, X, Y))

line()
print("Testing pymC2IXYS...")
X = 0.5791308486706011000e-3
Y = 0.4020579816732961219e-4
S = -0.1220040848472271978e-7
print(sf.pymC2IXYS(X, Y, S))

line()
print("Testing pymC2S...")
P = [100.0, -50.0, 25.0]
print(sf.pymC2S(P))

line()
print("Testing pymC2T00A...")
TTA = 2400000.5
UTA = 2400000.5
TTB = 53736.0
UTB = 53736.0
XP = 2.55060238e-7
YP = 1.860359247e-6
print(sf.pymC2T00A(TTA, TTB, UTA, UTB, XP, YP))

line()
print("Testing pymC2T00B...")
TTA = 2400000.5
UTA = 2400000.5
TTB = 53736.0
UTB = 53736.0
XP = 2.55060238e-7
YP = 1.860359247e-6
print(sf.pymC2T00B(TTA, TTB, UTA, UTB, XP, YP))

line()
print("Testing pymC2T06A...")
TTA = 2400000.5
UTA = 2400000.5
TTB = 53736.0
UTB = 53736.0
XP = 2.55060238e-7
YP = 1.860359247e-6
print(sf.pymC2T06A(TTA, TTB, UTA, UTB, XP, YP))

line()
print("Testing pymC2TCIO...")
RC2I = [[0.9999998323037164738, 0.5581526271714303683e-9, -0.5791308477073443903e-3],
        [-0.2384266227524722273e-7, 0.9999999991917404296, -0.4020594955030704125e-4],
        [0.5791308472168153320e-3, 0.4020595661593994396e-4, 0.9999998314954572365]]
ERA = 1.75283325530307
RPOM = [[0.9999999999999674705, -0.1367174580728847031e-10, 0.2550602379999972723e-6],
        [0.1414624947957029721e-10, 0.9999999999982694954, -0.1860359246998866338e-5],
        [-0.2550602379741215275e-6, 0.1860359247002413923e-5, 0.9999999999982369658]]
print(sf.pymC2TCIO(RC2I, ERA, RPOM))

line()
print("Testing pymC2TEQX...")
RBPN = [[0.9999989440476103608, -0.1332881761240011518e-2, -0.5790767434730085097e-3],
        [0.1332858254308954453e-2, 0.9999991109044505944, -0.4097782710401555759e-4],
        [0.5791308472168153320e-3, 0.4020595661593994396e-4, 0.9999998314954572365]]
GST = 1.754166138040730516
RPOM = [[0.9999999999999674705, -0.1367174580728847031e-10, 0.2550602379999972723e-6],
        [0.1414624947957029721e-10, 0.9999999999982694954, -0.1860359246998866338e-5],
        [-0.2550602379741215275e-6, 0.1860359247002413923e-5, 0.9999999999982369658]]
print(sf.pymC2TEQX(RBPN, GST, RPOM))

line()
print("Testing pymC2TPE...")
TTA = 2400000.5
UTA = 2400000.5
TTB = 53736.0
UTB = 53736.0
DEPS = 0.4090789763356509900
DPSI = -0.9630909107115582393e-5
XP = 2.55060238e-7
YP = 1.860359247e-6
print(sf.pymC2TPE(TTA, TTB, UTA, UTB, DPSI, DEPS, XP, YP))

line()
print("Testing pymC2TXY...")
TTA = 2400000.5
UTA = 2400000.5
TTB = 53736.0
UTB = 53736.0
X = 0.5791308486706011000e-3
Y = 0.4020579816732961219e-4
XP = 2.55060238e-7
YP = 1.860359247e-6
print(sf.pymC2TXY(TTA, TTB, UTA, UTB, X, Y, XP, YP))

line()
print("Testing pymCAL2JD...")
IY = 2003
IM  = 6
ID = 1
print(sf.pymCAL2JD(IY, IM, ID))

line()
print("Testing pymCP...")
P = [0.3, 1.2, -2.5]
print(sf.pymCP(P))

line()
print("Testing pymCPV...")
PV = [[0.3, 1.2, -2.5],
      [-0.5, 3.1, 0.9]]
print(sf.pymCPV(PV))

line()
print("Testing pymCR...")
R = [[2.0, 3.0, 2.0],
      [3.0, 2.0, 3.0],
      [3.0, 4.0, 5.0]]
print(sf.pymCR(R))

line()
print("Testing pymD2DTF...")
SCALE = 'UTC'
NDP = 5
D1 = 2400000.5
D2 = 49533.99999
print(sf.pymD2DTF(SCALE, NDP, D1, D2))

line()
print("Testing pymD2TF...")
NDP = 4
DAYS = -0.987654321
print(sf.pymD2TF(NDP, DAYS))

line()
print("Testing pymDAT...")
IY = 2017
IM  = 9
ID = 1
FD = 0.0
print(sf.pymDAT(IY, IM, ID, FD))

line()
print("Testing pymDTDB...")
DATE1 = 2448939.5
DATE2 = 0.123
UT = 0.76543
ELONG = 5.0123
U = 5525.242
V = 3190.0
print(sf.pymDTDB(DATE1, DATE2, UT, ELONG, U, V))

line()
print("Testing pymDTF2D...")
SCALE = 'UTC'
IY = 1994
IM  = 6
ID = 30
IHR = 23
IMN = 59
SEC = 60.13599
print(sf.pymDTF2D(SCALE, IY, IM, ID, IHR, IMN, SEC))

line()
print("Testing pymECEQ06...")
DATE1 = 2456165.5
DATE2 = 0.401182685
DL = 5.1
DB = -0.9
print(sf.pymECEQ06(DATE1, DATE2, DL, DB))

line()
print("Testing pymECM06...")
DATE1 = 2456165.5
DATE2 = 0.401182685
print(sf.pymECM06(DATE1, DATE2))

line()
print("Testing pymEE00...")
DATE1 = 2400000.5
DATE2 = 53736.0
EPSA = 0.4090789763356509900
DPSI = -0.9630909107115582393e-5
print(sf.pymEE00(DATE1, DATE2, EPSA, DPSI))

line()
print("Testing pymEE00A...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymEE00A(DATE1, DATE2))

line()
print("Testing pymEE00B...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymEE00B(DATE1, DATE2))

line()
print("Testing pymEE06A...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymEE06A(DATE1, DATE2))

line()
print("Testing pymEECT00...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymEECT00(DATE1, DATE2))

line()
print("Testing pymEFORM...")
N = 1
print(sf.pymEFORM(N))

line()
print("Testing pymEO06A...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymEO06A(DATE1, DATE2))

line()
print("Testing pymEORS...")
RNPB = [[0.9999989440476103608, -0.1332881761240011518e-2, -0.5790767434730085097e-3],
        [0.1332858254308954453e-2, 0.9999991109044505944, -0.4097782710401555759e-4],
        [0.5791308472168153320e-3, 0.4020595661593994396e-4, 0.9999998314954572365]]
S = -0.1220040848472271978e-7
print(sf.pymEORS(RNPB, S))

line()
print("Testing pymEPB...")
DJ1 = 2415019.8135
DJ2 = 30103.18648
print(sf.pymEPB(DJ1, DJ2))

line()
print("Testing pymEPB2JD...")
EPB = 1957.3
print(sf.pymEPB2JD(EPB))

line()
print("Testing pymEPJ...")
DJ1 = 2451545
DJ2 = -7392.5
print(sf.pymEPJ(DJ1, DJ2))

line()
print("Testing pymEPJ2JD...")
EPJ = 1996.8
print(sf.pymEPJ2JD(EPJ))

line()
print("Testing pymEPV00...")
DATE1 = 2400000.5
DATE2 = 53411.52501161
print(sf.pymEPV00(DATE1, DATE2))

line()
print("Testing pymEQEC06...")
DATE1 = 1234.5
DATE2 = 2440000.5
DR = 1.234
DD = 0.987
print(sf.pymEQEC06(DATE1, DATE2, DR, DD))

line()
print("Testing pymEQEQ94...")
DATE1 = 2400000.5
DATE2 = 41234.0
print(sf.pymEQEQ94(DATE1, DATE2))

line()
print("Testing pymERA00...")
DJ1 = 2400000.5
DJ2 = 54388.0
print(sf.pymERA00(DJ1, DJ2))

line()
print("Testing pymFAD03...")
T = 0.80
print(sf.pymFAD03(T))

line()
print("Testing pymFAE03...")
T = 0.80
print(sf.pymFAE03(T))

line()
print("Testing pymFAF03...")
T = 0.80
print(sf.pymFAF03(T))

line()
print("Testing pymFAJU03...")
T = 0.80
print(sf.pymFAJU03(T))

line()
print("Testing pymFAL03...")
T = 0.80
print(sf.pymFAL03(T))

line()
print("Testing pymFALP03...")
T = 0.80
print(sf.pymFALP03(T))

line()
print("Testing pymFAMA03...")
T = 0.80
print(sf.pymFAMA03(T))

line()
print("Testing pymFAME03...")
T = 0.80
print(sf.pymFAME03(T))

line()
print("Testing pymFANE03...")
T = 0.80
print(sf.pymFANE03(T))

line()
print("Testing pymFAOM03...")
T = 0.80
print(sf.pymFAOM03(T))

line()
print("Testing pymFAPA03...")
T = 0.80
print(sf.pymFAPA03(T))

line()
print("Testing pymFASA03...")
T = 0.80
print(sf.pymFASA03(T))

line()
print("Testing pymFAUR03...")
T = 0.80
print(sf.pymFAUR03(T))

line()
print("Testing pymFAVE03...")
T = 0.80
print(sf.pymFAVE03(T))

line()
print("Testing pymFK425...")
R1950 = 0.07626899753879587532
D1950 = -1.137405378399605780
DR1950 = 0.1973749217849087460e-4
DD1950 = 0.5659714913272723189e-5
P1950 = 0.134
V1950 = 8.7
print(sf.pymFK425(R1950, D1950, DR1950, DD1950, P1950, V1950))

line()
print("Testing pymFK45Z...")
R2000 = 0.01602284975382960982
D2000 = -0.1164347929099906024
BEPOCH = 1954.677617625256806
print(sf.pymFK45Z(R2000, D2000, BEPOCH))

line()
print("Testing pymFK524...")
R2000 = 0.8723503576487275595
D2000 = -0.7517076365138887672
DR2000 = 0.2019447755430472323e-4
DD2000 = 0.3541563940505160433e-5
P2000 = 0.1559
V2000 = 86.87
print(sf.pymFK524(R2000, D2000, DR2000, DD2000, P2000, V2000))

line()
print("Testing pymFK52H...")
R5  =  1.76779433
D5  = -0.2917517103
DR5 = -1.91851572e-7
DD5 = -5.8468475e-6
PX5 =  0.379210
RV5 = -7.6
print(sf.pymFK52H(R5, D5, DR5, DD5, PX5, RV5))

line()
print("Testing pymFK54Z...")
R2000 = 0.02719026625066316119
D2000 = -0.1115815170738754813
BEPOCH = 1954.677308160316374
print(sf.pymFK54Z(R2000, D2000, BEPOCH))

line()
print("Testing pymFK5HIP...")
print(sf.pymFK5HIP())

line()
print("Testing pymFK5HZ...")
R5 =  1.76779433
D5 = -0.2917517103
DATE1 = 2400000.5
DATE2 = 54479.0
print(sf.pymFK5HZ(R5, D5, DATE1, DATE2))

line()
print("Testing pymFW2M...")
GAMB = -0.2243387670997992368e-5
PHIB =  0.4091014602391312982
PSI  = -0.9501954178013015092e-3
EPS  =  0.4091014316587367472
print(sf.pymFW2M(GAMB, PHIB, PSI, EPS))

line()
print("Testing pymFW2XY...")
GAMB = -0.2243387670997992368e-5
PHIB =  0.4091014602391312982
PSI  = -0.9501954178013015092e-3
EPS  =  0.4091014316587367472
print(sf.pymFW2XY(GAMB, PHIB, PSI, EPS))

line()
print("Testing pymG2ICRS...")
DL =  5.5850536063818546461558105
DB = -0.7853981633974483096156608
print(sf.pymG2ICRS(DL, DB))

line()
print("Testing pymGC2GD...")
N = 2
XYZ = [2e6, 3e6, 5.244e6]
print(sf.pymGC2GD(N, XYZ))

line()
print("Testing pymGC2GDE...")
A = 6378136.0
F = 0.0033528
XYZ = [2e6, 3e6, 5.244e6]
print(sf.pymGC2GDE(A, F, XYZ))

line()
print("Testing pymGD2GC...")
N = 1
ELONG = 3.1
PHI = -0.5
HEIGHT = 2500.0
print(sf.pymGD2GC(N, ELONG, PHI, HEIGHT))

line()
print("Testing pymGD2GCE...")
A = 6378136.0
F = 0.0033528
ELONG = 3.1
PHI = -0.5
HEIGHT = 2500.0
print(sf.pymGD2GCE(A, F, ELONG, PHI, HEIGHT))

line()
print("Testing pymGMST00...")
UTA = 2400000.5
UTB = 53736.0
TTA = 2400000.5
TTB = 53736.0
print(sf.pymGMST00(UTA, UTB, TTA, TTB))

line()
print("Testing pymGMST06...")
UTA = 2400000.5
UTB = 53736.0
TTA = 2400000.5
TTB = 53736.0
print(sf.pymGMST06(UTA, UTB, TTA, TTB))

line()
print("Testing pymGMST82...")
DJ1 = 2400000.5
DJ2 = 53736.0
print(sf.pymGMST82(DJ1, DJ2))

line()
print("Testing pymGST00A...")
TTA = 2400000.5
UTA = 2400000.5
TTB = 53736.0
UTB = 53736.0
print(sf.pymGST00A(UTA, UTB, TTA, TTB))

line()
print("Testing pymGST00B...")
UTA = 2400000.5
UTB = 53736.0
print(sf.pymGST00B(UTA, UTB))

line()
print("Testing pymGST06...")
TTA = 2400000.5
UTA = 2400000.5
TTB = 53736.0
UTB = 53736.0
RNPB = [[0.9999989440476103608, -0.1332881761240011518e-2, -0.5790767434730085097e-3],
        [0.1332858254308954453e-2, 0.9999991109044505944, -0.4097782710401555759e-4],
        [0.5791308472168153320e-3, 0.4020595661593994396e-4, 0.9999998314954572365]]
print(sf.pymGST06(UTA, UTB, TTA, TTB, RNPB))

line()
print("Testing pymGST06A...")
TTA = 2400000.5
UTA = 2400000.5
TTB = 53736.0
UTB = 53736.0
print(sf.pymGST06A(UTA, UTB, TTA, TTB))

line()
print("Testing pymGST94...")
UTA = 2400000.5
UTB = 53736.0
print(sf.pymGST94(UTA, UTB))

line()
print("Testing pymICRS2G...")
DR =  5.9338074302227188048671087
DD = -1.1784870613579944551540570
print(sf.pymICRS2G(DR, DD))

line()
print("Testing pymH2FK5...")
RH  =  1.767794352
DH  = -0.2917512594
DRH = -2.76413026e-6
DDH = -5.92994449e-6
PXH =  0.379210
RVH = -7.6
print(sf.pymH2FK5(RH, DH, DRH, DDH, PXH, RVH))

line()
print("Testing pymHD2AE...")
HA = 1.1
DEC = 1.2
PHI = 0.3
print(sf.pymHD2AE(HA, DEC, PHI))

line()
print("Testing pymHD2PA...")
HA = 1.1
DEC = 1.2
PHI = 0.3
print(sf.pymHD2PA(HA, DEC, PHI))

line()
print("Testing pymHFK5Z...")
RH = 1.767794352
DH = -0.2917512594
DATE1 = 2400000.5
DATE2 = 54479.0
print(sf.pymHFK5Z(RH, DH, DATE1, DATE2))

line()
print("Testing pymIR...")
print(sf.pymIR())

line()
print("Testing pymJD2CAL...")
DJ1 = 2400000.5
DJ2 = 50123.9999
print(sf.pymJD2CAL(DJ1, DJ2))

line()
print("Testing pymJDCALF...")
N=4
DJ1 = 2400000.5
DJ2 = 50123.9999
print(sf.pymJDCALF(NDP, DJ1, DJ2))

line()
print("Testing pymLD...")
BM = 0.00028574
P = [-0.763276255, -0.608633767, -0.216735543]
Q = [-0.763276255, -0.608633767, -0.216735543]
E = [0.76700421, 0.605629598, 0.211937094]
EM = 8.91276983
DLIM = 3e-10
print(sf.pymLD(BM, P, Q, E, EM, DLIM))

line()
print("Testing pymLDN...")
N = 3
B = [[0.00028574, 3e-10, -7.81014427, -5.60956681,
      -1.98079819, 0.0030723249, -0.00406995477, -0.00181335842],
     [0.00095435, 3e-9, 0.738098796, 4.63658692, 
      1.9693136, -0.00755816922, 0.00126913722, 0.000727999001],
     [1.0, 6e-6, -0.000712174377, -0.00230478303, 
      -0.00105865966, 6.29235213e-6, -3.30888387e-7, -2.96486623e-7]]
OB = [-0.974170437, -0.2115201, -0.0917583114]
SC = [-0.763276255, -0.608633767, -0.216735543]
print(sf.pymLDN(N, B, OB, SC))

line()
print("Testing pymLDSUN...")
P = [-0.763276255, -0.608633767, -0.216735543]
E = [-0.973644023, -0.20925523, -0.0907169552]
EM = 0.999809214
print(sf.pymLDSUN(P, E, EM))
#####################################
line()
print("Testing pymLTECEQ...")
EPJ = 2500.0
DL = 1.5
DB = 0.6
print(sf.pymLTECEQ(EPJ, DL, DB))

line()
print("Testing pymLTECM...")
EPJ = -3000.0
print(sf.pymLTECM(EPJ))

line()
print("Testing pymLTEQEC...")
EPJ = -1500
DR = 1.234
DD = 0.987
print(sf.pymLTEQEC(EPJ, DR, DD))

line()
print("Testing pymLTP...")
EPJ = 1666.666
print(sf.pymLTP(EPJ))

line()
print("Testing pymLTPB...")
EPJ = 1666.666
print(sf.pymLTPB(EPJ))

line()
print("Testing pymLTPECL...")
EPJ = -1500
print(sf.pymLTPECL(EPJ))

line()
print("Testing pymLTPEQU...")
EPJ = -2500
print(sf.pymLTPEQU(EPJ))

line()
print("Testing pymMOON98...")
DATE1 = 2400000.5
DATE2 = 43999.9
print(sf.pymMOON98(DATE1, DATE2))

line()
print("Testing pymNUM00A...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymNUM00A(DATE1, DATE2))

line()
print("Testing pymNUM00B...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymNUM00B(DATE1, DATE2))

line()
print("Testing pymNUM06A...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymNUM06A(DATE1, DATE2))

line()
print("Testing pymNUMAT...")
EPSA =  0.4090789763356509900
DPSI = -0.9630909107115582393e-5
DEPS =  0.4063239174001678826e-4
print(sf.pymNUMAT(EPSA, DPSI, DEPS))

line()
print("Testing pymNUT00A...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymNUT00A(DATE1, DATE2))

line()
print("Testing pymNUT00B...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymNUT00B(DATE1, DATE2))

line()
print("Testing pymNUT06A...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymNUT06A(DATE1, DATE2))

line()
print("Testing pymNUT80...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymNUT80(DATE1, DATE2))

line()
print("Testing pymNUTM80...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymNUTM80(DATE1, DATE2))

line()
print("Testing pymOBL06...")
DATE1 = 2400000.5
DATE2 = 54388.0
print(sf.pymOBL06(DATE1, DATE2))

line()
print("Testing pymOBL80...")
DATE1 = 2400000.5
DATE2 = 54388.0
print(sf.pymOBL80(DATE1, DATE2))

line()
print("Testing pymP06E...")
DATE1 = 2400000.5
DATE2 = 52541.0
print(sf.pymP06E(DATE1, DATE2))

line()
print("Testing pymP2PV...")
P = [0.25, 1.2, 3.0]
print(sf.pymP2PV(P))

line()
print("Testing pymP2S...")
P = [100.0, -50.0, 25.0]
print(sf.pymP2S(P))

line()
print("Testing pymPAP...")
A = [1.0, 0.1, 0.2]
B = [-3.0, 1E-3, 0.2]
print(sf.pymPAP(A, B))

line()
print("Testing pymPAS...")
AL = 1.0
AP = 0.1
BL = 0.2
BP = -1.0
print(sf.pymPAS(AL, AP, BL, BP))

line()
print("Testing pymPB06...")
DATE1 = 2400000.5
DATE2 = 50123.9999
print(sf.pymPB06(DATE1, DATE2))

line()
print("Testing pymPDP...")
A = [2.0, 2.0, 3.0]
B = [1.0, 3.0, 4.0]
print(sf.pymPDP(A, B))

line()
print("Testing pymPFW06...")
DATE1 = 2400000.5
DATE2 = 50123.9999
print(sf.pymPFW06(DATE1, DATE2))

line()
print("Testing pymPLAN94...")
DATE1 = 2400000.5
DATE2 = 1E6
NP = 0
print(sf.pymPLAN94(DATE1, DATE2, NP))
DATE1 = 2400000.5
DATE2 = 1E6
NP = 10
print(sf.pymPLAN94(DATE1, DATE2, NP))
DATE1 = 2400000.5
DATE2 = -320000
NP = 3
print(sf.pymPLAN94(DATE1, DATE2, NP))
DATE1 = 2400000.5
DATE2 = 43999.9
NP = 1
print(sf.pymPLAN94(DATE1, DATE2, NP))

line()
print("Testing pymPMAT00...")
DATE1 = 2400000.5
DATE2 = 50123.9999
print(sf.pymPMAT00(DATE1, DATE2))

line()
print("Testing pymPMAT06...")
DATE1 = 2400000.5
DATE2 = 50123.9999
print(sf.pymPMAT06(DATE1, DATE2))

line()
print("Testing pymPMAT76...")
DATE1 = 2400000.5
DATE2 = 50123.9999
print(sf.pymPMAT76(DATE1, DATE2))

line()
print("Testing pymPM...")
P = [0.3, 1.2, -2.5]
print(sf.pymPM(P))

line()
print("Testing pymPMP...")
A = [2.0, 2.0, 3.0]
B = [1.0, 3.0, 4.0]
print(sf.pymPMP(A, B))

line()
print("Testing pymPMPX...")
RC = 1.234
DC = 0.789
PR = 1e-5
PD = -2e-5
PX = 1e-2
RV = 10.0
PMT = 8.75
POB = [0.9, 0.4, 0.1]
print(sf.pymPMPX(RC, DC, PR, PD, PX, RV, PMT, POB))

line()
print("Testing pymPMSAFE...")
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
print(sf.pymPMSAFE(RA1, DEC1, PMR1, PMD1, PX1, RV1, EP1A, EP1B, EP2A, EP2B))

line()
print("Testing pymPN...")
P = [0.3, 1.2, -2.5]
print(sf.pymPN(P))

line()
print("Testing pymPN00...")
DATE1 = 2400000.5
DATE2 = 53736.0
DPSI = -0.9632552291149335877e-5
DEPS =  0.4063197106621141414e-4
print(sf.pymPN00(DATE1, DATE2, DPSI, DEPS))

line()
print("Testing pymPN00A...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymPN00A(DATE1, DATE2))

line()
print("Testing pymPN00B...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymPN00B(DATE1, DATE2))

line()
print("Testing pymPN06A...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymPN06A(DATE1, DATE2))

line()
print("Testing pymPN06...")
DATE1 = 2400000.5
DATE2 = 53736.0
DPSI = -0.9632552291149335877e-5
DEPS =  0.4063197106621141414e-4
print(sf.pymPN06(DATE1, DATE2, DPSI, DEPS))

line()
print("Testing pymPNM00A...")
DATE1 = 2400000.5
DATE2 = 50123.9999
print(sf.pymPNM00A(DATE1, DATE2))

line()
print("Testing pymPNM00B...")
DATE1 = 2400000.5
DATE2 = 50123.9999
print(sf.pymPNM00B(DATE1, DATE2))

line()
print("Testing pymPNM06A...")
DATE1 = 2400000.5
DATE2 = 50123.9999
print(sf.pymPNM06A(DATE1, DATE2))

line()
print("Testing pymPNM80...")
DATE1 = 2400000.5
DATE2 = 50123.9999
print(sf.pymPNM80(DATE1, DATE2))

line()
print("Testing pymPOM00...")
XP = 2.55060238e-7
YP = 1.860359247e-6
SP = -0.1367174580728891460e-10
print(sf.pymPOM00(XP, YP, SP))

line()
print("Testing pymPPP...")
A = [2.0, 2.0, 3.0]
B = [1.0, 3.0, 4.0]
print(sf.pymPPP(A, B))

line()
print("Testing pymPPSP...")
A = [2.0, 2.0, 3.0]
S = 5.0
B = [1.0, 3.0, 4.0]
print(sf.pymPPSP(A, S, B))

line()
print("Testing pymPR00...")
DATE1 = 2400000.5
DATE2 = 53736
print(sf.pymPR00(DATE1, DATE2))

line()
print("Testing pymPREC76...")
DATE01 = 2400000.5
DATE02 = 33282.0
DATE11 = 2400000.5
DATE12 = 51544.0
print(sf.pymPREC76(DATE01, DATE02, DATE11, DATE12))

line()
print("Testing pymPV2P...")
PV = [[0.3, 1.2, -2.5],
      [-0.5, 3.1, 0.9]]
print(sf.pymPV2P(PV))

line()
print("Testing pymPV2S...")
PV = [[-0.4514964673880165, 0.03093394277342585, 0.05594668105108779],
      [1.292270850663260e-5, 2.652814182060692e-6, 2.568431853930293e-6]]
print(sf.pymPV2S(PV))

line()
print("Testing pymPVDPV...")
A = [[2.0, 2.0, 3.0],
      [6.0, 0.0, 4.0]]
B = [[1.0, 3.0, 4.0],
      [0.0, 2.0, 8.0]]
print(sf.pymPVDPV(A, B))

line()
print("Testing pymPVM...")
PV = [[0.3, 1.2, -2.5],
      [0.45, -0.25, 1.1]]
print(sf.pymPVM(PV))

line()
print("Testing pymPVMPV...")
A = [[2.0, 2.0, 3.0],
      [5.0, 6.0, 3.0]]
B = [[1.0, 3.0, 4.0],
      [3.0, 2.0, 1.0]]
print(sf.pymPVMPV(A, B))

line()
print("Testing pymPVPPV...")
A = [[2.0, 2.0, 3.0],
      [5.0, 6.0, 3.0]]
B = [[1.0, 3.0, 4.0],
      [3.0, 2.0, 1.0]]
print(sf.pymPVPPV(A, B))

line()
print("Testing pymPVSTAR...")
PV = [[126668.5912743160601, 2136.792716839935195, -245251.2339876830091],
      [-0.4051854035740712739e-2, -0.6253919754866173866e-2, 0.1189353719774107189e-1]]
print(sf.pymPVSTAR(PV))

line()
print("Testing pymPVTOB...")
ELONG = 2.0
PHI = 0.5
HM = 3000.0
XP = 1e-6
YP = -0.5e-6
SP = 1e-8
THETA = 5.0
print(sf.pymPVTOB(ELONG, PHI, HM, XP, YP, SP, THETA))

line()
print("Testing pymPVU...")
DT = 2920.0
PV = [[126668.5912743160601, 2136.792716839935195, -245251.2339876830091],
      [-0.4051854035740712739e-2, -0.6253919754866173866e-2, 0.1189353719774107189e-1]]
print(sf.pymPVU(DT, PV))

line()
print("Testing pymPVUP...")
DT = 2920.0
PV = [[126668.5912743160601, 2136.792716839935195, -245251.2339876830091],
      [-0.4051854035740712739e-2, -0.6253919754866173866e-2, 0.1189353719774107189e-1]]
print(sf.pymPVUP(DT, PV))

line()
print("Testing pymPVXPV...")
A = [[2.0, 2.0, 3.0],
      [6.0, 0.0, 4.0]]
B = [[1.0, 3.0, 4.0],
      [0.0, 2.0, 8.0]]
print(sf.pymPVXPV(A, B))

line()
print("Testing pymPXP...")
A = [2.0, 2.0, 3.0]
B = [1.0, 3.0, 4.0]
print(sf.pymPXP(A, B))

line()
print("Testing pymREFCO...")
PHPA = 800.0
TC = 10.0
RH = 0.9
WL = 0.4
print(sf.pymREFCO(PHPA, TC, RH, WL))

line()
print("Testing pymRM2V...")
R = [[0.00, -0.80, -0.60],
     [0.80, -0.36, 0.48],
     [0.60, 0.48, -0.64]]
print(sf.pymRM2V(R))

line()
print("Testing pymRV2M...")
W = [0.0, 1.41371669, -1.88495559]
print(sf.pymRV2M(W))

line()
print("Testing pymRX...")
PHI = 0.3456789
R = [[2.0, 3.0, 2.0],
      [3.0, 2.0, 3.0],
      [3.0, 4.0, 5.0]]
print(sf.pymRX(PHI, R))

line()
print("Testing pymRXP...")
R = [[2.0, 3.0, 2.0],
      [3.0, 2.0, 3.0],
      [3.0, 4.0, 5.0]]
P = [0.2, 1.5, 0.1]
print(sf.pymRXP(R, P))

line()
print("Testing pymRXPV...")
R = [[2.0, 3.0, 2.0],
      [3.0, 2.0, 3.0],
      [3.0, 4.0, 5.0]]
PV = [[0.2, 1.5, 0.1],
      [1.5, 0.2, 0.1]]
print(sf.pymRXPV(R, PV))

line()
print("Testing pymRXR...")
A = [[2.0, 3.0, 2.0],
      [3.0, 2.0, 3.0],
      [3.0, 4.0, 5.0]]
B = [[1.0, 2.0, 2.0],
      [4.0, 1.0, 1.0],
      [3.0, 0.0, 1.0]]
print(sf.pymRXR(A, B))

line()
print("Testing pymRY...")
THETA = 0.3456789
R = [[2.0, 3.0, 2.0],
      [3.0, 2.0, 3.0],
      [3.0, 4.0, 5.0]]
print(sf.pymRY(THETA, R))

line()
print("Testing pymRZ...")
PSI = 0.3456789
R = [[2.0, 3.0, 2.0],
      [3.0, 2.0, 3.0],
      [3.0, 4.0, 5.0]]
print(sf.pymRZ(PSI, R))

line()
print("Testing pymS00A...")
DATE1 = 2400000.5
DATE2 = 52541.0
print(sf.pymS00A(DATE1, DATE2))

line()
print("Testing pymS00B...")
DATE1 = 2400000.5
DATE2 = 52541.0
print(sf.pymS00B(DATE1, DATE2))

line()
print("Testing pymS00...")
DATE1 = 2400000.5
DATE2 = 53736.0
X = 0.5791308486706011000e-3
Y = 0.4020579816732961219e-4
print(sf.pymS00(DATE1, DATE2, X, Y))

line()
print("Testing pymS06A...")
DATE1 = 2400000.5
DATE2 = 52541.0
print(sf.pymS06A(DATE1, DATE2))

line()
print("Testing pymS06...")
DATE1 = 2400000.5
DATE2 = 53736.0
X = 0.5791308486706011000e-3
Y = 0.4020579816732961219e-4
print(sf.pymS06(DATE1, DATE2, X, Y))

line()
print("Testing pymS2C...")
THETA = 3.0123
PHI = -0.999
print(sf.pymS2C(THETA, PHI))

line()
print("Testing pymS2P...")
THETA = -3.21
PHI = 0.123
R = 0.456
print(sf.pymS2P(THETA, PHI, R))

line()
print("Testing pymS2PV...")
THETA = -3.21
PHI = 0.123
R = 0.456
TD = -7.8e-6
PD = 9.01e-6
RD = -1.23e-5
print(sf.pymS2PV(THETA, PHI, R, TD, PD, RD))

line()
print("Testing pymS2XPV...")
S1 = 2.0
S2 = 3.0
PV = [[0.3, 1.2, -2.5],
      [0.5, 2.3, -0.4]]
print(sf.pymS2XPV(S1, S2, PV))

line()
print("Testing pymSEPP...")
A = [1.0, 0.1, 0.2]
B = [-3.0, 1E-3, 0.2]
print(sf.pymSEPP(A, B))

line()
print("Testing pymSEPS...")
AL =  1.0
AP =  0.1
BL =  0.2
BP = -3.0
print(sf.pymSEPS(AL, AP, BL, BP))

line()
print("Testing pymSP00...")
DATE1 = 2400000.5
DATE2 = 52541.0
print(sf.pymSP00(DATE1, DATE2))

line()
print("Testing pymSTARPM...")
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
print(sf.pymSTARPM(RA1, DEC1, PMR1, PMD1, PX1, RV1, EP1A, EP1B, EP2A, EP2B))

line()
print("Testing pymSTARPV...")
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
print(sf.pymSTARPV(RA, DEC, PMR, PMD, PX, RV))

line()
print("Testing pymSXP...")
S = 2.0
P = [0.3, 1.2, -2.5]
print(sf.pymSXP(S, P))

line()
print("Testing pymSXPV...")
S = 2.0
PV = [[0.3, 1.2, -2.5],
      [0.5, 3.2, -0.7]]
print(sf.pymSXPV(S, PV))

line()
print("Testing pymTAITT...")
TAI1 = 2453750.5
TAI2 = 0.892482639
print(sf.pymTAITT(TAI1, TAI2))

line()
print("Testing pymTAIUT1...")
TAI1 = 2453750.5
TAI2 = 0.892482639
DTA = -32.6659
print(sf.pymTAIUT1(TAI1, TAI2, DTA))

line()
print("Testing pymTAIUTC...")
TAI1 = 2453750.5
TAI2 = 0.892482639
print(sf.pymTAIUTC(TAI1, TAI2))

line()
print("Testing pymTCBTDB...")
TCB1 = 2453750.5
TCB2 = 0.893019599
print(sf.pymTCBTDB(TCB1, TCB2))

line()
print("Testing pymTCGTT...")
TCG1 = 2453750.5
TCG2 = 0.892862531
print(sf.pymTCGTT(TCG1, TCG2))

line()
print("Testing pymTDBTCB...")
TDB1 = 2453750.5
TDB2 = 0.892855137
print(sf.pymTDBTCB(TDB1, TDB2))

line()
print("Testing pymTDBTT...")
TDB1 = 2453750.5
TDB2 = 0.892855137
DTR = -0.000201
print(sf.pymTDBTT(TDB1, TDB2, DTR))

line()
print("Testing pymTF2A...")
S = '+'
IHOUR = 4
IMIN = 58
SEC = 20.2
print(sf.pymTF2A(S, IHOUR, IMIN, SEC))

line()
print("Testing pymTF2D...")
S = ' '
IHOUR = 23
IMIN = 55
SEC = 10.9
print(sf.pymTF2D(S, IHOUR, IMIN, SEC))

line()
print("Testing pymTPORS...")
XI = -0.03
ETA = 0.07
A = 1.3
B = 1.5
print(sf.pymTPORS(XI, ETA, A, B))

line()
print("Testing pymTPORV...")
XI = -0.03
ETA = 0.07
V = sf.pymS2C(1.3, 1.5)
print(sf.pymTPORV(XI, ETA, V))

line()
print("Testing pymTPSTS...")
XI = -0.03
ETA = 0.07
A0 = 2.3
B0 = 1.5
print(sf.pymTPSTS(XI, ETA, A0, B0))

line()
print("Testing pymTPSTV...")
XI = -0.03
ETA = 0.07
A0 = 2.3
B0 = 1.5
V0 = sf.pymS2C(A0, B0)
print(sf.pymTPSTV(XI, ETA, V0))

line()
print("Testing pymTPXES...")
A = 1.3
B = 1.55
A0 = 2.3
B0 = 1.5
print(sf.pymTPXES(A, B, A0, B0))

line()
print("Testing pymTPXEV...")
A = 1.3
B = 1.55
A0 = 2.3
B0 = 1.5
V = sf.pymS2C(A, B)
V0 = sf.pymS2C(A0, B0)
print(sf.pymTPXEV(V, V0))

line()
print("Testing pymTR...")
R = [[2.0, 3.0, 2.0],
      [3.0, 2.0, 3.0],
      [3.0, 4.0, 5.0]]
print(sf.pymTR(R))

line()
print("Testing pymTRXP...")
R = [[2.0, 3.0, 2.0],
      [3.0, 2.0, 3.0],
      [3.0, 4.0, 5.0]]
P = [0.2, 1.5, 0.1]
print(sf.pymTRXP(R, P))

line()
print("Testing pymTRXPV...")
R = [[2.0, 3.0, 2.0],
      [3.0, 2.0, 3.0],
      [3.0, 4.0, 5.0]]
PV = [[0.2, 1.5, 0.1],
      [1.5, 0.2, 0.1]]
print(sf.pymTRXPV(R, PV))

line()
print("Testing pymTTTAI...")
TT1 = 2453750.5
TT2 = 0.892482639
print(sf.pymTTTAI(TT1, TT2))

line()
print("Testing pymTTTCG...")
TT1 = 2453750.5
TT2 = 0.892482639
print(sf.pymTTTCG(TT1, TT2))

line()
print("Testing pymTTTDB...")
TT1 = 2453750.5
TT2 = 0.892855139
DTR = -0.000201
print(sf.pymTTTDB(TT1, TT2, DTR))

line()
print("Testing pymTTUT1...")
TT1 = 2453750.5
TT2 = 0.892855139
DT = 64.8499
print(sf.pymTTUT1(TT1, TT2, DT))

line()
print("Testing pymUT1TAI...")
UT11 = 2453750.5
UT12 = 0.892104561
DTA = -32.6659
print(sf.pymUT1TAI(UT11, UT12, DTA))

line()
print("Testing pymUT1TT...")
UT11 = 2453750.5
UT12 = 0.892104561
DTA = 64.8499
print(sf.pymUT1TT(UT11, UT12, DT))

line()
print("Testing pymUT1UTC...")
UT11 = 2453750.5
UT12 = 0.892104561
DUT1 = 0.3341
print(sf.pymUT1UTC(UT11, UT12, DUT1))

line()
print("Testing pymUTCTAI...")
UTC1 = 2453750.5
UTC2 = 0.892100694
print(sf.pymUTCTAI(UTC1, UTC2))

line()
print("Testing pymUTCUT1...")
UTC1 = 2453750.5
UTC2 = 0.892100694
DUT1 = 0.3341
print(sf.pymUTCUT1(UTC1, UTC2, DUT1))

line()
print("Testing pymXY06...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymXY06(DATE1, DATE2))

line()
print("Testing pymXYS00A...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymXYS00A(DATE1, DATE2))

line()
print("Testing pymXYS00B...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymXYS00B(DATE1, DATE2))

line()
print("Testing pymXYS06A...")
DATE1 = 2400000.5
DATE2 = 53736.0
print(sf.pymXYS06A(DATE1, DATE2))

line()
print("Testing pymZP..")
print(sf.pymZP())

line()
print("Testing pymZPV..")
print(sf.pymZPV())

line()
print("Testing pymZR..")
print(sf.pymZR())
