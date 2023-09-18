#An example of the use of pymsofa in CHES. 
#Used to simulate the image of a target star in the tangent plane 
#when it is observed at point L2

import PyMsOfa  as  sf
import numpy  as  np
import xlrd
import math as ma
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

xl=xlrd.open_workbook(r'.\input.xls')
xltar=xl.sheets()[0]
xlref=xl.sheets()[1]
xlref15=xl.sheets()[2]

#input the information of the target star
tarinfo=[xltar.cell(2,i).value for i in range(10)]
tarinfo[1]=ma.radians(tarinfo[1])
tarinfo[2]=ma.radians(tarinfo[2])
tarinfo[4]=ma.radians(tarinfo[4]/1000/3600)
tarinfo[5]=ma.radians(tarinfo[5]/1000/3600)
tarinfo[6]=ma.radians(tarinfo[6]/1000/3600)
tarinfo[7]=ma.radians(tarinfo[7]/1000/3600)


#input the information of the reference star (Mag<13)
refinfo=[xlref.cell(1,0).value,int(xlref.cell(1,1).value)]
refinfo.append([[0 for k in range(9)] for i in range(refinfo[1])])
for i in range(refinfo[1]):
    for k in range(9):
        refinfo[2][i][k]=float(xlref.cell(2+i,k+1).value)
        if (k==0)|(k==1):
            refinfo[2][i][k]=ma.radians(refinfo[2][i][k])
        elif (k==3)|(k==4)|(k==5)|(k==6):
            refinfo[2][i][k]=ma.radians(refinfo[2][i][k]/1000/3600)

#input the information of the reference star (13<Mag<15)
refinfo15=[xlref15.cell(1,0).value,int(xlref15.cell(1,1).value)]
refinfo15.append([[0 for k in range(9)] for i in range(refinfo15[1])])

for i in range(refinfo15[1]):
    for k in range(9):
        refinfo15[2][i][k]=float(xlref15.cell(2+i,k+1).value)
        if (k==0)|(k==1):
            refinfo15[2][i][k]=ma.radians(refinfo15[2][i][k])
        elif (k==3)|(k==4)|(k==5)|(k==6):
            refinfo15[2][i][k]=ma.radians(refinfo15[2][i][k]/1000/3600)


#The estimated observation time(year, month, day, hour, minute, second)
T=[2025, 11, 12, 12, 56, 15]
D = sf.pymTf2d('+', 12, 56, 15)
JD1, JD2 = sf.pymCal2jd(2025, 11, 12)
JD2 +=D

#The epoch of the Gaia catalog
JD01, JD02 = sf.pymEpj2jd(2016.0)

#ICRS coordinates at the time of observation
tarcoord=[0,0]
refcoord=[[0,0] for i in range(refinfo[1])]
refcoord15=[[0,0] for i in range(refinfo15[1])]

tarcoord[0], tarcoord[1], pmr2,pmd2,px2,rv2\
    = sf.pymPmsafe(tarinfo[1], tarinfo[2], tarinfo[4], tarinfo[6],
                   tarinfo[3], tarinfo[9], JD01, JD02, JD1, JD2)

for i in range(refinfo[1]):
    refcoord[i][0], refcoord[i][1], pmr2,pmd2,px2,rv2\
        = sf.pymPmsafe(refinfo[2][i][0], refinfo[2][i][1], 
                       refinfo[2][i][3], refinfo[2][i][5], 
                       refinfo[2][i][2], refinfo[2][i][8], 
                       JD01, JD02, JD1, JD2)

for i in range(refinfo15[1]):
    refcoord15[i][0], refcoord15[i][1], pmr2,pmd2,px2,rv2\
        = sf.pymPmsafe(refinfo15[2][i][0], refinfo15[2][i][1], 
                       refinfo15[2][i][3], refinfo15[2][i][5],
                       refinfo15[2][i][2], refinfo15[2][i][8], 
                       JD01, JD02, JD1, JD2)

#Transform the coordinates of the star to 
#the station-centred coordinate system of the L2 point

#The longitude of the sun relative to the satellite at the time of observation
JD10, JD20 = sf.pymCal2jd(2025, 3, 21)
lamda_sun, beta_sun = (JD2-JD20)/365.25*2*ma.pi, 0

#Convert to equatorial coordinates
TAI1,TAI2 = sf.pymUtctai(JD1, JD2)
TT1, TT2 = sf .pymTaitt(TAI1, TAI2)

alpha_sun, delta_sun = sf.pymEceq06(TT1,TT2, lamda_sun, beta_sun)
suncoord= [alpha_sun, delta_sun]

#Position of the star in L2 station-centred coordinates
tartopo=[0,0]

cosAR=ma.cos(ma.pi/2-tarcoord[1])*ma.cos(ma.pi/2-suncoord[1])\
    +ma.sin(ma.pi/2-tarcoord[1])*ma.sin(ma.pi/2-suncoord[1])*\
                             ma.cos(np.abs(tarcoord[0]-suncoord[0]))
sinAR=ma.sqrt(1-cosAR**2)

sinp=sinAR*(149597870+1500000)/(1000/tarinfo[3]*30856775814671.915808)
p=ma.asin(sinp)

k=-p/sinAR

tartopo[0]=tarcoord[0]+k/ma.cos(tarcoord[1])*\
    ma.cos(suncoord[1])*ma.sin(tarcoord[0]-suncoord[0])

tartopo[1]=tarcoord[1]+k*(ma.sin(tarcoord[1])*\
    ma.cos(suncoord[1])*ma.cos(tarcoord[0]-suncoord[0])\
        -ma.cos(tarcoord[1])*ma.sin(suncoord[1]))

reftopo=[[0,0] for k in range(refinfo[1])]

for k in range(refinfo[1]):
    cosAR=ma.cos(ma.pi/2-refcoord[k][1])*ma.cos(ma.pi/2-suncoord[1])\
        +ma.sin(ma.pi/2-refcoord[k][1])*ma.sin(ma.pi/2-suncoord[1])*\
                                 ma.cos(np.abs(refcoord[k][0]-suncoord[0]))
    sinAR=ma.sqrt(1-cosAR**2)

    sinp=sinAR*(149597870+1500000)/(1000/refinfo[2][k][2]*30856775814671.915808)
    p=ma.asin(sinp)

    kk=-p/sinAR

    reftopo[k][0]=refcoord[k][0]+kk/ma.cos(refcoord[k][1])*\
        ma.cos(suncoord[1])*ma.sin(refcoord[k][0]-suncoord[0])

    reftopo[k][1]=refcoord[k][1]+kk*(ma.sin(refcoord[k][1])*\
        ma.cos(suncoord[1])*ma.cos(refcoord[k][0]-suncoord[0])\
            -ma.cos(refcoord[k][1])*ma.sin(suncoord[1]))

reftopo15=[[0,0] for k in range(refinfo15[1])]

for k in range(refinfo15[1]):
    cosAR=ma.cos(ma.pi/2-refcoord15[k][1])*ma.cos(ma.pi/2-suncoord[1])\
        +ma.sin(ma.pi/2-refcoord15[k][1])*ma.sin(ma.pi/2-suncoord[1])*\
                                 ma.cos(np.abs(refcoord15[k][0]-suncoord[0]))
    sinAR=ma.sqrt(1-cosAR**2)

    sinp=sinAR*(149597870+1500000)/(1000/refinfo15[2][k][2]*30856775814671.915808)
    p=ma.asin(sinp)

    kk=-p/sinAR

    reftopo15[k][0]=refcoord15[k][0]+kk/ma.cos(-refcoord15[k][1])*\
        ma.cos(-suncoord[1])*ma.sin(refcoord15[k][0]-suncoord[0])

    reftopo15[k][1]=refcoord15[k][1]+kk*(-ma.sin(refcoord15[k][1])*\
        ma.cos(-suncoord[1])*ma.cos(refcoord15[k][0]-suncoord[0])\
            -ma.cos(refcoord15[k][1])*ma.sin(-suncoord[1]))

#According to the projection theorem, using the target star as the tangent point, 
#the position of the reference star on the projection plane
refplt = [[0,0] for k in range(refinfo[1])]

for i in range(refinfo[1]):
    refplt[i][0], refplt[i][1] = sf.pymTpxes(reftopo[i][0],reftopo[i][1],
                                             tartopo[0], tartopo[1])

refplt15 = [[0,0] for k in range(refinfo15[1])]

for i in range(refinfo15[1]):
    refplt15[i][0], refplt15[i][1] = sf.pymTpxes(reftopo15[i][0],reftopo15[i][1],
                                                 tartopo[0], tartopo[1])

#diagram based on the coordinates
fig, ax = plt.subplots(figsize=(10,10))
color1=mcolors.to_rgba('white',alpha=1)
color2=mcolors.to_rgba('white',alpha=0.9)
color3=mcolors.to_rgba('white',alpha=0.4)

plt.scatter(0,0,c=color1, marker='.',s=1500)
for k in range(refinfo[1]):
    plt.scatter(-refplt[k][0],refplt[k][1],c=color2, marker='.',s=100)
    if k<9:
        if k==5:
            plt.text(-refplt[k][0]+0.00013,refplt[k][1],'R 0%d'%(k+1),c='w')
        else:
            plt.text(-refplt[k][0]-0.00013,refplt[k][1]+0.00008,'R 0%d'%(k+1),c='w')
    else:
        plt.text(-refplt[k][0]-0.00013,refplt[k][1]+0.00008,'R %d'%(k+1),c='w')
for k in range(refinfo15[1]):
    plt.scatter(-refplt15[k][0],refplt15[k][1],c=color3, marker='.',s=50)
plt.scatter(0,0,c='r',marker='+',s=400)
plt.scatter(0,0,c='w',marker=',',s=10)
plt.xticks([])
plt.yticks([])
plt.rcParams['axes.facecolor']='black'

plt.show()      