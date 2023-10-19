import math as ma
import numpy as np


def pymD2tf(NDP,DAYS):
    '''
    Decompose days to hours, minutes, seconds, fraction.

    Parameters
    ----------
    ndp : int
        resolution    
    days : float
        interval in days

    Returns
    -------
    sign : str
        '+' or '-'    
    ihmsf : list(4)    
        hours, minutes, seconds, fraction

         NDP         resolution
          :      ...0000 00 00
         -7         1000 00 00
         -6          100 00 00
         -5           10 00 00
         -4            1 00 00
         -3            0 10 00
         -2            0 01 00
         -1            0 00 10
          0            0 00 01
          1            0 00 00.1
          2            0 00 00.01
          3            0 00 00.001
          :            0 00 00.000...

    '''
    if DAYS>=0:
        SIGN='+'
    else:
        SIGN='-'
    
    #一天的秒数
    D2S=86400.0
    A=D2S*np.abs(DAYS)
    
    #根据精度要求处理
    if NDP<0:
        NRS=1.0
        N=1
        while N<=-NDP:
           if (N==2)|(N==4):
               NRS=NRS*6
           else:
               NRS=NRS*10
           N+=1
        RS=float(NRS)
        A=RS*float(int(A/RS+0.5))
    
    NRS=1.0
    N=1
    while N<=NDP:
        NRS=NRS*10
        N+=1
    RS=float(NRS)
    RM=RS*60.0
    RH=RM*60.0
    
    A=float(int(RS*A+0.5))
    
    AH=float(int(A/RH))
    A=A-AH*RH
    AM=float(int(A/RM))
    A=A-AM*RM
    AS=float(int(A/RS))
    AF=A-AS*RS
    
    IHMSF=[int(AH+0.5),int(AM+0.5),int(AS+0.5),int(AF+0.5)]
    return(SIGN,IHMSF)


def pymA2af(NDP,ANGLE):
    '''
    Decompose radians into degrees, arcminutes, arcseconds, fraction.

    Parameters
    ----------
    ndp : int
        resolution      
    angle : float
        angle in radians

    Returns
    -------
    sign : str
        '+' or '-'    
    idmsf : list(4)
        degrees, arcminutes, arcseconds, fraction 

         NDP         resolution
          :      ...0000 00 00
         -7         1000 00 00
         -6          100 00 00
         -5           10 00 00
         -4            1 00 00
         -3            0 10 00
         -2            0 01 00
         -1            0 00 10
          0            0 00 01
          1            0 00 00.1
          2            0 00 00.01
          3            0 00 00.001
          :            0 00 00.000...
    '''
    #调用D2TF中，天转换为时分秒的函数，通过F将列表中第一项小时转换为角度
    F=15.0/(2*ma.pi)
    SIGN,IDMSF=pymD2tf(NDP,ANGLE*F)

    return(SIGN,IDMSF)


def pymA2tf (NDP,ANGLE):
    '''
    Decompose radians into hours, minutes, seconds, fraction.

    Parameters
    ----------
    ndp : int
        resolution    
    angle : float
        angle in radians

    Returns
    -------
    sign : str 
        '+' or '-'    
    ihmsf : list(4)
        
    NDP         resolution
     :      ...0000 00 00
    -7         1000 00 00
    -6          100 00 00
    -5           10 00 00
    -4            1 00 00
    -3            0 10 00
    -2            0 01 00
    -1            0 00 10
     0            0 00 01
     1            0 00 00.1
     2            0 00 00.01
     3            0 00 00.001
     :            0 00 00.000...
    '''
    #直接调用天到时分秒的函数
    SIGN,IHMSF=pymD2tf(NDP,ANGLE/(2*ma.pi))
    return(SIGN,IHMSF)


def pymTf2d(S,IHOUR,IMIN,SEC):
    '''
    Convert hours, minutes, seconds to days

    Parameters
    ----------
    s : bytes
        sign:  '-' = negative, otherwise positive
    ihour : int
        hours
    imin : int
        minutes
    sec : float
        seconds
    
    Returns
    -------
    days : float
        interval in days
    J : ValueError
        1:'ihour outside range 0-23',
        2:'imin outside range 0-59',
        3:'sec outside range 0-59.999...'
    '''

    #一天的秒数
    D2S=86400.0
    
    #判断输入数据是否有误
    J=0
    if (SEC<0)|(SEC>=60):
        J=3
    if (IMIN<0)|(IMIN>59):
        J=2
    if (IHOUR<0)|(IHOUR>23):
        J=1
        
    W=(60.0*(60.0*float(np.abs(IHOUR))+float(np.abs(IMIN)))+np.abs(SEC))/D2S
    
    if (S=='-'):
        W=-W
    
    DAYS=W

    return(DAYS,J)


def pymCal2jd(IY,IM,ID):
    '''
    Gregorian Calendar to Julian Date.

    Parameters
    ----------
    iyear : int
        year in Gregorian calendar
    imon : int
        month in Gregorian calendar
    iday : int
        day in Gregorian calendar

    Returns
    -------
    djm0 : float
        MJD zero-point: always 2400000.5
    djm : float
        Modified Julian Date for 0 hrs
    J : ValueError
        -1: 'bad year,  the year is simply valid from -4800 March 1',
        -2: 'bad month, the month is not from 1 to 12',
        -3: 'bad day,   the day is not related to the month'
    '''

    #函数允许的最早日期
    IYMIN=-4799
    #每个月的天长
    MTAB=[31,28,31,30,31,30,31,31,30,31,30,31]
    
    J=0
    
    if (IY<IYMIN):
        J=-1
    else:
        if (IM>=1)&(IM<=12):
            #当月的天数
            NDAYS=MTAB[int(IM-1)]
            #判断是否闰年
            if (IM==2):
                if (IY%4==0):
                    NDAYS=29
                if (IY%100==0)&(IY%400!=0):
                    NDAYS=28
            if (ID<1)|(ID>NDAYS):
                J=-3
                
            MY=int((IM-14)/12)
            IYPMY=int(IY+MY)
            DJM0=2400000.5
            #儒略日自-4800年3月1日起
            DJM=float(int((1461*(IYPMY+4800))/4)+int((367*(IM-2-12*MY))/12)-\
                      int(3*(int((IYPMY+4900)/100))/4)+ID-2432076)
        else:
            J=-2
    return(DJM0,DJM,J)


def pymJd2cal(DJ1,DJ2):
    '''
    Julian Date to Gregorian year, month, day, and fraction of a day.

    Parameters
    ----------
    dj1 : float
        Julian Date    
    dj2 : float
        Julian Date

    Returns
    -------
    iy : int
        yaer    
    im : int
        month
    id : int
        day
    fd : float
        fraction of day
    J : ValueError
        -1: 'The valid date is -68569.5 (-4900 March 1) up to 1e9'
    '''
    #极小值，使1.0+EPS!=1.0
    EPS=2.2204460492503131e-16
    #最小，最大允许的输入儒略日
    DJMIN=-68569.5
    DJMAX=1e9
    
    DJ=DJ1+DJ2
    if (DJ<DJMIN)|(DJ>DJMAX):
        J=-1
    else:
        J=0
    
    #将天数和小数分离开(其中小数部分绝对值小于0.5).
        D=float(int(DJ1+0.5))
        F1=DJ1-D
        JD=int(D)
        D=float(int(DJ2+0.5))
        F2=DJ2-D
        JD=JD+int(D)
    
    #使用补偿求和计算F1+F2+0.5(Klein 2006)
        S=0.5
        CS=0.0
        V=[F1,F2]
        for i in range(2):
            X=V[i]
            T=S+X
            if (np.abs(S)>=np.abs(X)):
                C=(S-T)+X
            else:
                C=(X-T)+S
            CS=CS+C
            S=T
            if (S>=1.0):
                JD+=1
                S=S-1.0
        F=S+CS
        CS=F-S
            
    #当F为负数时.
        if (F<0.0):
        #补偿求和，假设 |S| <= 1.
            F=S+1.0
            CS=CS+((1.0-F)+S)
            S=F
            F=S+CS
            CS=F-S
            JD=JD-1
        
    #当F=1.0或者更大时.
        if ((F-1.0)>=(-EPS/4.0)):
        #补偿求和，假设 |S| <= 1. */
            T=S-1.0
            CS=CS+((S-T)-1.0)
            S=T
            F=S+CS
            if ((-EPS/2.0)<F):
                JD=JD+1
                F=max(F,0.0)
            
    #在公历下表示日期.
        L=JD+68569
        N=int((4*L)/146097)
        L=int(L-int((146097*N+3)/4))
        I=int((4000*(L+1))/1461001)
        L=int(L-int((1461*I)/4)+31)
        K=int((80*L)/2447)
        ID=int(L-int((2447*K)/80))
        L=int(K/11)
        IM=int(K+2-12*L)
        IY=int(100*(N-49)+I+L)
        FD=F      
        
    return(IY,IM,ID,FD,J)



def pymJdcalf(NDP,DJ1,DJ2):
    '''
    Julian Date to Gregorian Calendar, expressed in a form convenient 
    for formatting messages:  rounded to a specified precision.

    Parameters
    ----------
    ndp : int
        number of decimal places of days in fraction    
    dj1 : float
        Julian Date    
    dj2 : float
        Julian Date    

    Returns
    -------
    iymdf : list(4)
        year, month, day, fraction in Gregorian calendar

    '''
    #小数的位数 (e.g. NDP取2时，DENOM取100)
    if (NDP>=0)&(NDP<=9):
        J=0
        DENOM=float(10**NDP)
    else:
        J=1
        DENOM=1.0
    
    #D1为大值，D2为小值
    if (np.abs(DJ1)>=np.abs(DJ2)):
        D1=DJ1
        D2=DJ2
    else:
        D1=DJ2
        D2=DJ1            
    
    #调整到午夜(没有舍入误差)。
    D1=D1-0.5
    
    #尽可能精确地分离开日和小数.
    D=float(int(D1+0.5))
    F1=D1-D
    DJD=D
    D=float(int(D2+0.5))
    F2=D2-D
    DJD=DJD+D;
    D=float(int(F1+F2+0.5))
    F=(F1-D)+F2
    if (F<0.0):
        F=F+1.0
        D=D-1.0
    DJD+=D
    
    #在指定的精度要求下四舍五入.
    RF=float(int(F*DENOM+0.5))/DENOM
    
    #重新调整到中午
    DJD+=0.5
    
    #转换到公历.
    IYMDF=[0,0,0,0]
    IYMDF[0],IYMDF[1],IYMDF[2],F,JS=pymJd2cal(DJD, RF)
    if (JS==0):
        IYMDF[3]=int(F*DENOM+0.5)
    else:
        J=JS    
    
    return(IYMDF, J)


def pymEpb(DJ1,DJ2):
    '''
    Julian Date to Besselian Epoch.

    Parameters
    ----------
    dj1 : float
        Julian Date    
    dj2 : float
        Julian Date

    Returns
    -------
    function value : float
        Besselian Epoch

    '''

    #J2000.0
    DJ00=2451545.0
    #J2000.0 减去 B1900.0 (2415019.81352)的天数
    D1900=36524.68648
    #B1900的回归年长度(天)
    TY=365.242198781

    EPB=1900.0+((DJ1-DJ00)+(DJ2+D1900))/TY

    return(EPB)


def pymEpb2jd(EPB):
    '''
    Besselian Epoch to Julian Date

    Parameters
    ----------
    epb : float
        Besselian Epoch

    Returns
    -------
    djm0 : float
        MJD zero-point: always 2400000.5    
    djm : float
        Modified Julian Date

    '''

    #B1900的回归年长度(天)
    TY=365.242198781
    
    DJM0=2400000.5
    DJM=15019.81352+(EPB-1900.0)*TY
        
    return(DJM0,DJM)


def pymEpj(DJ1,DJ2):
    '''
    Julian Date to Julian Epoch

    Parameters
    ----------
    dj1 : float
        Julian Date    
    dj2 : float
        Julian Date

    Returns
    -------
    function value : float
        Julian Epoch

    '''

    #参考日期(J2000.0), JD
    DJ00=2451545.0
    #儒略年的天数
    DJY=365.25
    
    EPJ=2000.0+((DJ1-DJ00)+DJ2)/DJY
    
    return(EPJ)


def pymEpj2jd (EPJ):
    '''
    Julian Epoch to Julian Date

    Parameters
    ----------
    epj : float
        Julian Epoch (e.g. 1996.8)

    Returns
    -------
    djm0 : float
        MJD zero-point: always 2400000.5     
    djm : float
        Modified Julian Date

    '''
    DJM0=2400000.5
    DJM=51544.5+(EPJ-2000.0)*365.25
    
    return(DJM0,DJM)


def pymDat(IY,IM,ID,FD):
    
    '''
    For a given UTC date, calculate Delta(AT) = TAI-UTC.

    Parameters
    ----------
    iy : int
        year    
    im : int
        month    
    iday : int
        day    
    fday : float    
        fraction of day

    Returns
    -------
    deltat : float
        TAI minus UTC, seconds
    J : ValueError
        1: 'dubious year', 
       -1: 'bad year,  the year is simply valid from -4800 March 1',
       -2: 'bad month, the month is not from 1 to 12',
       -3: 'bad day,   the day is not related to the month',
       -4: 'bad fraction of day',
       -5: 'internal error', 
    '''

    #此函数版本发布日期
    IYV=2023
    '闰秒需要调整的参数，修改时年份'
    
    #闰秒的次数（当有新的闰秒出现时，+1）
    NDAT = 42
    '闰秒需要调整的参数，+1'
    
    #在闰秒被引入之前的时间差变化数目
    NERA1 = 14
    '为1972年前的次数，无需调整'
    IDAT=[[1960,1],[1961,1],[1961,8],[1962,1],[1963,11],[1964,1],[1964,4],
          [1964,9],[1965,1],[1965,3],[1965,7],[1965,9],[1966,1],[1968,2],
          [1972,1],[1972,7],[1973,1],[1974,1],[1975,1],[1976,1],[1977,1],
          [1978,1],[1979,1],[1980,1],[1981,7],[1982,7],[1983,7],[1985,7],
          [1988,1],[1990,1],[1991,1],[1992,7],[1993,7],[1994,7],[1996,1],
          [1997,7],[1999,1],[2006,1],[2009,1],[2012,7],[2015,7],[2017,1]]
    '闰秒需要调整的参数，最后增加一项闰秒时间'
    DATS=[1.4178180,1.4228180,1.3728180,1.8458580,1.9458580,3.2401300,
          3.3401300,3.4401300,3.5401300,3.6401300,3.7401300,3.8401300,
          4.3131700,4.2131700,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,
          18.0,19.0,20.0,21.0,22.0,23.0,24.0,25.0,26.0,27.0,28.0,29.0,
          30.0,31.0,32.0,33.0,34.0,35.0,36.0,37.0]
    '闰秒需要调整的参数，最后增加一项总闰秒时间，通常+-1'
    DRIFT=[[37300, 0.001296],[37300, 0.001296],[37300, 0.001296],
           [37665, 0.0011232],[37665, 0.0011232],[38761, 0.001296],
           [38761, 0.001296],[38761, 0.001296],[38761, 0.001296],
           [38761, 0.001296],[38761, 0.001296],[38761, 0.001296],
           [39126, 0.002592],[39126, 0.002592]]
    '为1972年前数据，无需调整'
    
    #初始化
    DA=0.0
    JS=0
    
    #用循环替代goto函数
    i=1
    while i<2:

        if (FD<0.0)|(FD>1.0):
            JS=-4
            print('ERROR1',JS)
            break
        
        #将日期转换为儒略日.
        DJM0,DJM,JS=pymCal2jd(IY,IM,ID)
        if (JS<0):
            print('ERROR2',JS)
            break
        
        if (IY<IDAT[0][0]):
            #早于1960年，报错
            JS=1
            print('ERROR3',JS)
            break
        #晚于函数版本5年后，报错
        if (IY>(IYV+5)):
            JS=1
        
        #年月结合
        M=12*IY+IM
        
        #找到最接近的项目.
        IS=0
        MORE=True
        N=NDAT
        while N>0:
            if MORE:
                IS=N
                MORE=M<(12*IDAT[N-1][0]+IDAT[N-1][1])    
            N-=1
        
        if (IS<1):
            JS=-5
            print('ERROR4',JS)
            break
        
        DA=DATS[IS-1]
        
        #早于1972年时，需要一个额外的调整
        if (IS<=NERA1):
            DA=DA+(DJM+FD-DRIFT[IS-1][0])*DRIFT[IS-1][1]
        
        i+=1
    DELTAT=DA
    J=JS
    
    return(DELTAT,J)


def pymDtf2d(SCALE,IY,IM,ID,IHR,IMN,SEC):
    '''
    Encode date and time fields into 2-part Julian Date (or in the case
    of UTC a quasi-JD form that includes special provision for leap seconds).

    Parameters
    ----------
    scale : str
        time scale ID(Only the value "UTC" is significant)
    iy : int
        year in Gregorian calendar
    im : int
        month in Gregorian calendar
    id : int
        day in Gregorian calendar
    ihr : int
        hour
    imn : int
        minute
    sec : float
        seconds

    Returns
    -------
    d1,d2 : float
        2-part Julian Date
    J : ValueError
        3: 'both of next two',
        2: 'time is after end of day',
        1: 'dubious year', 
       -1: 'bad year,  the year is simply valid from -4800 March 1',
       -2: 'bad month, the month is not from 1 to 12',
       -3: 'bad day,   the day is not related to the month',
       -4: 'bad hour',
       -5: 'bad minute',
       -6: 'bad second', 
    '''
    #一天的秒数
    D2S = 86400.0
    D1,D2=0,0
    i=1
    while i<2:
    
        #当日的儒略日数
        DJ,W,JS=pymCal2jd(IY,IM,ID)
    
        if (JS!=0):
            print('ERROR',JS)
            break
        DJ=DJ+W
        
        #一天的长度以及暂时认为的最后一分钟的分长（用于闰秒情况）
        DAY=D2S
        SECLIM=60.0
        
        #处理在UTC中出现闰秒时的情况
        if (SCALE=='UTC'):
            
            #当天0时的TAI-UTC时间差.
            DAT0,JS=pymDat(IY,IM,ID,0.0)
            if (JS<0):
                print('ERROR',JS)
                break
            
            #当天12h的TAI-UTC时间差（判断是否存在drift——1972年之前情况）
            DAT12,JS=pymDat(IY,IM,ID,0.5)
            if (JS<0):
                print('ERROR',JS)
                break
           
            #第二天0h的TAI-UTC时间差（判断是否存在跳秒）
            IY2,IM2,ID2,W,JS=pymJd2cal(DJ,1.5)
            if (JS!=0):
                print('ERROR',JS)
                break
            DAT24,JS=pymDat(IY2,IM2,ID2,0.0)
            if (JS<0):
                print('ERROR',JS)
                break
            
            #当天和第二天之间TAI-UTC是否有突然变化。
            DLEAP=DAT24-(2.0*DAT12-DAT0)
            
            #如果是闰秒日，改正日常以及当天最后一分钟的分长
            DAY=DAY+DLEAP
            if (IHR==23)&(IMN==59):
                SECLIM+=DLEAP
        
        #验证时间.
        if (IHR>=0)&(IHR<=23):
            if (IMN>=0)&(IMN<=59):
                if (SEC>=0.0):
                    if (SEC>=SECLIM):
                        JS=JS+2
                else:
                    JS=-6
            else:
                JS=-5
        else:
            JS=-4
        
        if (JS<0):
            print('ERROR',JS)
            break
        
        #以天为单位的时间.
        TIME=(60.0*float(60*IHR+IMN)+SEC)/DAY
        
        D1 = DJ
        D2 = TIME
  
        i+=1
    J=JS
    
    return(D1,D2,J)


def pymD2dtf(SCALE,NDP,D1,D2):
    '''
    Format for output a 2-part Julian Date (or in the case of UTC a 
    quasi-JD form that includes special provision for leap seconds).

    Parameters
    ----------
    scale : str
        time scale ID(Only the value "UTC" is significant)       
    ndp : int
        resolution    
    d1 : float
        time as a 2-part Julian Date    
    d2 : float
        time as a 2-part Julian Date

    Returns
    -------
    iy,im,id : int
        year, month, day in Gregorian calendar    
    ihmsf : list(4)
        hours, minutes, seconds, fraction
    J : ValueError
        1: 'dubious year',
       -1: 'unacceptable date',
       
         NDP         resolution
          :      ...0000 00 00
         -7         1000 00 00
         -6          100 00 00
         -5           10 00 00
         -4            1 00 00
         -3            0 10 00
         -2            0 01 00
         -1            0 00 10
          0            0 00 01
          1            0 00 00.1
          2            0 00 00.01
          3            0 00 00.001
          :            0 00 00.000...
    '''
    
    #一天的秒长
    D2S = 86400.0
    IY,IM,ID=0,0,0
    IHMSF=[0,0,0,0]
    i=1
    while i<2:
    
        A1 = D1
        B1 = D2
        
        #临时日历日期.
        IY1,IM1,ID1,FD,JS=pymJd2cal(A1,B1)
        if (JS!=0):
            print('ERROR',JS)
            break
        
        #判断是否是闰秒当天
        LEAP=False
        if (SCALE=='UTC'):
            
            #当天0时的TAI-UTC时间差.
            DAT0,JS=pymDat(IY1,IM1,ID1,0.0)
            if (JS<0):
                print('ERROR',JS)
                break
            
            #当天12h的TAI-UTC时间差（判断是否存在drift——1972年之前情况）
            DAT12,JS=pymDat(IY1,IM1,ID1,0.5)
            if (JS<0):
                print('ERROR',JS)
                break
            
            #第二天0h的TAI-UTC时间差（判断是否存在跳秒）
            IY2,IM2,ID2,W,JS=pymJd2cal(A1+1.5,B1-FD)
            if (JS!=0):
                print('ERROR',JS)
                break
            DAT24,JS=pymDat(IY2,IM2,ID2,0.0)
            if (JS<0):
                print('ERROR',JS)
                break
    
            #当天和第二天之间TAI-UTC是否有突然变化。
            DLEAP=DAT24-(2.0*DAT12-DAT0)
            
            #如果是闰秒日，则修改小数部分.
            LEAP=np.abs(DLEAP)>0.5
            if (LEAP):
                FD=FD+FD*DLEAP/D2S
        
        #临时时间.
        S,IHMSF1=pymD2tf(NDP,FD)
            
        #判断四舍五入的时间是否超过24点
        if (IHMSF1[0]>23):
            
            #是，则日期为第二天
            IY2,IM2,ID2,W,JS=pymJd2cal(A1+1.5,B1-FD)
            if (JS<0):
                print('ERROR',JS)
                break
            
            #判断是否当天有闰秒
            if (not LEAP):
                
                #否，则用第二天的零时
                IY1=IY2
                IM1=IM2
                ID1=ID2
                IHMSF1[0]=0
                IHMSF1[1]=0
                IHMSF1[2]=0
            
            else:
                
                #是，则判断是否已经过了闰秒的时间
                if (IHMSF1[2]>0):
                    
                    #是，则第二天要考虑到闰秒。
                    IY1=IY2
                    IM1=IM2
                    ID1=ID2
                    IHMSF1[0]=0
                    IHMSF1[1]=0
                    IHMSF1[2]=0
                
                else:
                    
                    #否，则当天采用23：59：60
                    IHMSF1[0]=23
                    IHMSF1[1]=59
                    IHMSF1[2]=60
            
                #如果四舍五入到10s，且大致上是要到第二天
                if (NDP<0)&(IHMSF1[2]==60):
                    IY1=IY2
                    IM1=IM2
                    ID1=ID2
                    IHMSF1[0]=0
                    IHMSF1[1]=0
                    IHMSF1[2]=0
        
        IY=IY1
        IM=IM1
        ID=ID1
        for i in range(4):
           IHMSF[i]=IHMSF1[i]
        
        i+=1
    J=JS
    
    return(IY,IM,ID,IHMSF,J)


def pymUtctai(UTC1,UTC2):
    '''
    Time scale transformation:  Coordinated Universal Time, UTC, to
    International Atomic Time, TAI.

    Parameters
    ----------
    utc1 : float
        UTC as a 2-part quasi Julian Date
    utc2 : float
        UTC as a 2-part quasi Julian Date

    Returns
    -------
    tai1 : float
        TAI as a 2-part Julian Date
    tai2 : float
        TAI as a 2-part Julian Date
    J : ValueError
        1: 'dubious year',
       -1: 'unacceptable date',
    '''
    TAI1,TAI2=0,0
    
    #一天的秒长
    D2S=86400.0
    
    #将两参数日期按先大后小排列
    BIG1=(np.abs(UTC1)>=np.abs(UTC2))
    if (BIG1):
        U1=UTC1
        U2=UTC2
    else:
        U1=UTC2
        U2=UTC1
    
    i=1
    while i<2:
        
        #当天0时的TAI-UTC时间差.
        IY,IM,ID,FD,JS=pymJd2cal(U1,U2)
        if (JS!=0):
            print('ERROR1',JS)
            break
        DAT0,JS=pymDat(IY,IM,ID,0.0)
        if (JS<0):
            print('ERROR2',JS)
            break
    
        #当天12h的TAI-UTC时间差（判断是否存在drift——1972年之前情况）
        DAT12,JS=pymDat(IY,IM,ID,0.5)
        if (JS<0):
            print('ERROR3',JS)
            break
    
        #第二天0h的TAI-UTC时间差（判断是否存在跳秒）
        IYT,IMT,IDT,W,JS=pymJd2cal(U1+1.5,U2-FD)
        if (JS!=0):
            print('ERROR4',JS)
            break
        DAT24,JS=pymDat(IYT,IMT,IDT,0.0)
        if (JS<0):
            print('ERROR5',JS)
            break

        #将TAI和UTC之间的时间差分解到两个部分，前一天的时差和跳秒的
        DLOD=2.0*(DAT12-DAT0)
        DLEAP=DAT24-(DAT0+DLOD)
        
        #在前一天中去除掉跳秒
        FD=FD*(D2S+DLEAP)/D2S
        
        #从(1972年以前)UTC秒到SI秒
        FD=FD*(D2S+DLOD)/D2S

        #将当天的日期转换到儒略日.
        Z1,Z2,JS=pymCal2jd(IY,IM,ID)
        if (JS!=0):
            print('ERROR6',JS)
            break

        #构建TAI的时间表达
        A2=Z1-U1
        A2=(A2+Z2)+(FD+DAT0/D2S)
        if (BIG1):
            TAI1=U1
            TAI2=A2
        else:
            TAI1=A2
            TAI2=U1
       
        i+=1
    J=JS
    
    return(TAI1,TAI2,J)


def pymTaiutc(TAI1,TAI2):
    '''
    Time scale transformation:  International Atomic Time, TAI, to
    Coordinated Universal Time, UTC.

    Parameters
    ----------
    tai1 : float
        TAI as a 2-part Julian Date
    tai2 : float
        TAI as a 2-part Julian Date

    Returns
    -------
    utc1 : float
        UTC as a 2-part quasi Julian Date
    utc2 : float
        UTC as a 2-part quasi Julian Date
    J : ValueError
        1: 'dubious year',
       -1: 'unacceptable date',
    '''

    UTC1,UTC2=0,0
    
    #将两参数日期按先大后小排列
    BIG1=np.abs(TAI1)>=np.abs(TAI2)
    if (BIG1):
        A1=TAI1
        A2=TAI2
    else:
        A1=TAI2
        A2=TAI1
    
    #初始假设日期
    U1=A1
    U2=A2
    
    i=1
    while i<2:
        
        #迭代
        for j in range(3):
            
            #假设UTC
            G1,G2,JS=pymUtctai(U1, U2)
            if (JS<0):
                print('ERROR',JS)
                break
        
            #调整UTC
            U2=U2+(A1-G1)
            U2=U2+(A2-G2)
    
        #得到UTC的结果
        if (BIG1):
            UTC1=U1
            UTC2=U2
        else:
            UTC1=U2
            UTC2=U1
    
        i+=1
    J=JS
    return(UTC1,UTC2,J)


def pymAf2a(S,IDEG,IAMIN,ASEC):
    '''
    Convert degrees, arcminutes, arcseconds to radians.

    Parameters
    ----------
    s : str
        sign:  '-' = negative, otherwise positive    
    ideg : int
        degrees    
    iamin : int
        arcminutes    
    asec : float
        arcseconds    

    Returns
    -------
    rad : float
        angle in radians
    J : ValueError
        1: 'ideg outside range 0-359',
        2: 'iamin outside range 0-59',
        3: 'asec outside range 0-59.999...',
    '''

    #1角秒对应的弧度
    DAS2R=4.848136811095359935899141e-6

    #重置J
    J=0
    
    #判断输入合法性
    if (ASEC<0.0)|(ASEC>=60.0):
        J=3
    if (IAMIN<0)|(IAMIN>59):
        J=2
    if (IDEG<0)|(IDEG>359):
        J=1
    
    #计算角度
    W=(60.0*(60.0*float(np.abs(IDEG))+float(np.abs(IAMIN)))+np.abs(ASEC))*DAS2R

    #应用符号.
    if (S=='-'):
        W=-W
    RAD=W
    
    return(RAD,J)

def pymEform(N):
    '''
    Earth reference ellipsoids.

    Parameters
    ----------
    n : int
        ellipsoid identifier
        
    n    ellipsoid    
    1     WGS84    
    2     GRS80    
    3     WGS72    
    
    Returns
    -------
    a : float
        equatorial radius (meters)    
    f : float
        flattening
    J : ValueError
        0: 'correct',
       -1: 'n error',
    '''

    J=0
    
    #根据编号查找对应的参数
    if (N==1):
        
        #WGS84.
        A=6378137.0
        F=1.0/298.257223563
    
    elif (N==2):
        
        #GRS80.
        A=6378137.0
        F=1.0/298.257222101
    
    elif (N==3):
        
        #WGS72.
        A=6378135.0
        F=1.0/298.26
    
    else:
        
        #编号错误情况
        A=0.0
        F=0.0
        J=-1
    
    return(A,F,J)


def pymGd2gce(A,F,ELONG,PHI,HEIGHT):
    '''
    Transform geodetic coordinates to geocentric for a reference
    ellipsoid of specified form.

    Parameters
    ----------
    a : float
        equatorial radius    
    f : float
        flattening     
    elong : float
        longitude (radians, east +ve)    
    phi : float
        latitude (geodetic, radians)    
    height : float
        height above ellipsoid (geodetic)

    Returns
    -------
    xyz : list(3)
        geocentric vector
    J : ValueError
        0: 'correct',
       -1: 'illegal case',
    '''

    XYZ=[0,0,0]
    
    #计算地心纬度.
    SP=ma.sin(PHI)
    CP=ma.cos(PHI)
    W=1.0-F
    W=W*W
    D=CP*CP+W*SP*SP
    if (D>0.0):
        AC=A/ma.sqrt(D)
        AS=W*AC
    
        #地心矢量.
        R=(AC+HEIGHT)*CP
        XYZ[0]=R*ma.cos(ELONG)
        XYZ[1]=R*ma.sin(ELONG)
        XYZ[2]=(AS+HEIGHT)*SP
        
        J=0
    else:
        J=-1
    
    return(XYZ,J)


def pymZp():
    '''
    Zero a p-vector

    Parameters
    ----------

    Returns
    -------
    p : list(3)
        zero p-vector
    '''

    P=[0,0,0]
    
    return(P)


def pymGd2gc (N,ELONG,PHI,HEIGHT):
    '''
    Transform geodetic coordinates to geocentric using the specified
    reference ellipsoid.

    Parameters
    ----------
    n : int
        ellipsoid identifier
    elong : float
        longitude (radians, east +ve)    
    phi : float
        latitude (geodetic, radians)    
    height : float
        height above ellipsoid (geodetic)

    n    ellipsoid    
    1     WGS84    
    2     GRS80    
    3     WGS72    
    
    Returns
    -------
    xyz : list(3)
        geocentric vector
    J : ValueError
        0 : OK
       -1 : illegal identifier 
       -2 : illegal case 
    '''
    XYZ=[0,0,0]
    
    #获得参考的椭球参数
    A,F,J=pymEform(N)
    
    #如果没有报错，调用pymGd2gce函数进行处理
    if (J==0):
        XYZ,J=pymGd2gce(A,F,ELONG,PHI,HEIGHT)
        if (J!=0):
            J=-2
    
    if (J!=0):
        XYZ=pymZp()
        
    return(XYZ,J)


def pymTaiut1(TAI1,TAI2,DTA):
    '''
    Time scale transformation:  International Atomic Time, TAI, to
    Universal Time, UT1.

    Parameters
    ----------
    tai1 : float
        TAI as a 2-part Julian Date    
    tai2 : float
        TAI as a 2-part Julian Date    
    dta : float
        UT1-TAI in seconds

    Returns
    -------
    ut11,ut12 : float
        UT1 as a 2-part Julian Date
    J : ValueError
        0 : OK
    '''
    #一天的秒长
    D2S=86400.0
    
    #维持精度
    DTAD=DTA/D2S
    if (np.abs(TAI1)>np.abs(TAI2)):
        UT11=TAI1
        UT12=TAI2+DTAD
    else:
        UT11=TAI1+DTAD
        UT12=TAI2
    
    J=0
    return(UT11,UT12,J)


def pymUtcut1(UTC1,UTC2,DUT1):
    '''
    Time scale transformation:  Coordinated Universal Time, UTC, to
    Universal Time, UT1.

    Parameters
    ----------
    utc1 : float
        UTC as a 2-part quasi Julian Date
    utc2 : float
        UTC as a 2-part quasi Julian Date
    dut1 : float
        Delta UT1 = UT1-UTC in seconds

    Returns
    -------
    ut11 : float
        UT1 as a 2-part Julian Date
    ut12 : float
        UT1 as a 2-part Julian Date
    J : ValueError
        1: 'dubious year',
       -1: 'unacceptable date',    
    '''   
    i=1
    while i<2:
            
        #调用TAI-UTC之间的关系.
        IY,IM,ID,W,JS=pymJd2cal(UTC1,UTC2)
        if (JS!=0):
            print('ERROR',JS)
            break
        DAT,JS=pymDat(IY,IM,ID,0.0)
        if (JS<0):
            print('ERROR',JS)
            break
    
        #构建UT1-TAI
        DTA=DUT1-DAT
        
        #由UTC到TAI，再到UT1
        TAI1,TAI2,JW=pymUtctai(UTC1,UTC2)
        if (JW<0):
            JS=JW
            print('ERROR',JS)
            break
        UT11,UT12,JW=pymTaiut1(TAI1,TAI2,DTA)
        
        i+=1
    J=JS
    return(UT11,UT12,J)


def pymTaitt(TAI1,TAI2):
    '''
    Time scale transformation:  International Atomic Time, TAI, to
    Terrestrial Time, TT.

    Parameters
    ----------
    tai1 : float
        TAI as a 2-part Julian Date    
    tai2 : float
        TAI as a 2-part Julian Date

    Returns
    -------
    tt1,tt2 : float
        TT as a 2-part Julian Date
    J : ValueError
        0 : OK
    '''
    #地球时与国际原子时之间的时差（单位：天）
    DTAT=32.184/86400.0
    
    #尽可能保留精度的结果
    if (np.abs(TAI1)>np.abs(TAI2)):
        TT1=TAI1
        TT2=TAI2+DTAT
    else:
        TT1=TAI1+DTAT
        TT2=TAI2
    
    J=0
    return(TT1,TT2,J)


def pymTttcg(TT1,TT2):
    '''
    Time scale transformation:  Terrestrial Time, TT, to Geocentric
    Coordinate Time, TCG.

    Parameters
    ----------
    tt1 : float
        TT as a 2-part Julian Date
    tt2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    tcg1 : float
        TCG as a 2-part Julian Date
    tcg2 : float
        TCG as a 2-part Julian Date
    J : ValueError
        0 : OK
    '''    
    #约简儒略日零点
    DJM0=2400000.5
    
    #1977 Jan 1 00:00:32.184 TT, 约简儒略日下的表达
    T77T=43144.0003725
    
    #L_G = 1 - dTT/dTCG，TT和TCG之间关系的定义
    ELG=6.969290134e-10
    
    #TT 与 TCG 之间的关系
    ELGG=ELG/(1.0-ELG)
    
    #尽可能保留精度的结果
    if (np.abs(TT1)>np.abs(TT2)):
        TCG1=TT1
        TCG2=TT2+((TT1-DJM0)+(TT2-T77T))*ELGG
    else:
        TCG1=TT1+((TT2-DJM0)+(TT1-T77T))*ELGG
        TCG2=TT2
    
    J=0    
    return(TCG1,TCG2,J)


def pymDtdb(DATE1,DATE2,UT,ELONG,U,V):
    '''
    An approximation to TDB-TT, the difference between barycentric
    dynamical time and terrestrial time, for an observer on the Earth.

    Parameters
    ----------
    date1 : float
        date, TDB    
    date2 : float
        date, TDB    
    ut : float
        universal time (UT1, fraction of one day)    
    elong : float
        longitude (east positive, radians)     
    u : float
        distance from Earth spin axis (km)    
    v : float
        distance north of equatorial plane (km)

    Returns
    -------
    DTDB : float
        TDB-TT (seconds)

    '''

    D2PI=2*ma.pi
    #1度对应的弧度
    DD2R=1.745329251994329576923691e-2
    #参考时期 (J2000.0), 儒略日期
    DJ00=2451545.0
    #每千儒略年的天数
    DJM=365250.0
# =============================================================================
#     *  =====================
#     *  Fairhead et al. model
#     *  =====================
#     *
#     *  787 sets of three coefficients.
#     *
#     *  Each set is amplitude (microseconds)
#     *              frequency (radians per Julian millennium since J2000.0),
#     *              phase (radians).
#     *
#     *  Sets   1-474 are the T**0 terms,
#     *   "   475-679  "   "  T**1   "
#     *   "   680-764  "   "  T**2   "
#     *   "   765-784  "   "  T**3   "
#     *   "   785-787  "   "  T**4   "  .
# =============================================================================
    #模型中的参数
    #FAIRHD(787,3)
    FAIRHD=[[1656.674564e-6,6283.075849991,6.240054195],
            [22.417471e-6,5753.384884897,4.296977442],
            [13.839792e-6,12566.151699983,6.196904410],
            [4.770086e-6,529.690965095,0.444401603],
            [4.676740e-6,6069.776754553,4.021195093],
            [2.256707e-6,213.299095438,5.543113262],
            [1.694205e-6,-3.523118349,5.025132748],
            [1.554905e-6,77713.771467920,5.198467090],
            [1.276839e-6,7860.419392439,5.988822341],
            [1.193379e-6,5223.693919802,3.649823730],
            [1.115322e-6,3930.209696220,1.422745069],
            [0.794185e-6,11506.769769794,2.322313077],
            [0.447061e-6,26.298319800,3.615796498],
            [0.435206e-6,-398.149003408,4.349338347],
            [0.600309e-6,1577.343542448,2.678271909],
            [0.496817e-6,6208.294251424,5.696701824],
            [0.486306e-6,5884.926846583,0.520007179],
            [0.432392e-6,74.781598567,2.435898309],
            [0.468597e-6,6244.942814354,5.866398759],
            [0.375510e-6,5507.553238667,4.103476804],
            [0.243085e-6,-775.522611324,3.651837925],
            [0.173435e-6,18849.227549974,6.153743485],
            [0.230685e-6,5856.477659115,4.773852582],
            [0.203747e-6,12036.460734888,4.333987818],
            [0.143935e-6,-796.298006816,5.957517795],
            [0.159080e-6,10977.078804699,1.890075226],
            [0.119979e-6,38.133035638,4.551585768],
            [0.118971e-6,5486.777843175,1.914547226],
            [0.116120e-6,1059.381930189,0.873504123],
            [0.137927e-6,11790.629088659,1.135934669],
            [0.098358e-6,2544.314419883,0.092793886],
            [0.101868e-6,-5573.142801634,5.984503847],
            [0.080164e-6,206.185548437,2.095377709],
            [0.079645e-6,4694.002954708,2.949233637],
            [0.062617e-6,20.775395492,2.654394814],
            [0.075019e-6,2942.463423292,4.980931759],
            [0.064397e-6,5746.271337896,1.280308748],
            [0.063814e-6,5760.498431898,4.167901731],
            [0.048042e-6,2146.165416475,1.495846011],
            [0.048373e-6,155.420399434,2.251573730],
            [0.058844e-6,426.598190876,4.839650148],
            [0.046551e-6,-0.980321068,0.921573539],
            [0.054139e-6,17260.154654690,3.411091093],
            [0.042411e-6,6275.962302991,2.869567043],
            [0.040184e-6,-7.113547001,3.565975565],
            [0.036564e-6,5088.628839767,3.324679049],
            [0.040759e-6,12352.852604545,3.981496998],
            [0.036507e-6,801.820931124,6.248866009],
            [0.036955e-6,3154.687084896,5.071801441],
            [0.042732e-6,632.783739313,5.720622217],
            [0.042560e-6,161000.685737473,1.270837679],
            [0.040480e-6,15720.838784878,2.546610123],
            [0.028244e-6,-6286.598968340,5.069663519],
            [0.033477e-6,6062.663207553,4.144987272],
            [0.034867e-6,522.577418094,5.210064075],
            [0.032438e-6,6076.890301554,0.749317412],
            [0.030215e-6,7084.896781115,3.389610345],
            [0.029247e-6,-71430.695617928,4.183178762],
            [0.033529e-6,9437.762934887,2.404714239],
            [0.032423e-6,8827.390269875,5.541473556],
            [0.027567e-6,6279.552731642,5.040846034],
            [0.029862e-6,12139.553509107,1.770181024],
            [0.022509e-6,10447.387839604,1.460726241],
            [0.020937e-6,8429.241266467,0.652303414],
            [0.020322e-6,419.484643875,3.735430632],
            [0.024816e-6,-1194.447010225,1.087136918],
            [0.025196e-6,1748.016413067,2.901883301],
            [0.021691e-6,14143.495242431,5.952658009],
            [0.017673e-6,6812.766815086,3.186129845],
            [0.022567e-6,6133.512652857,3.307984806],
            [0.016155e-6,10213.285546211,1.331103168],
            [0.014751e-6,1349.867409659,4.308933301],
            [0.015949e-6,-220.412642439,4.005298270],
            [0.015974e-6,-2352.866153772,6.145309371],
            [0.014223e-6,17789.845619785,2.104551349],
            [0.017806e-6,73.297125859,3.475975097],
            [0.013671e-6,-536.804512095,5.971672571],
            [0.011942e-6,8031.092263058,2.053414715],
            [0.014318e-6,16730.463689596,3.016058075],
            [0.012462e-6,103.092774219,1.737438797],
            [0.010962e-6,3.590428652,2.196567739],
            [0.015078e-6,19651.048481098,3.969480770],
            [0.010396e-6,951.718406251,5.717799605],
            [0.011707e-6,-4705.732307544,2.654125618],
            [0.010453e-6,5863.591206116,1.913704550],
            [0.012420e-6,4690.479836359,4.734090399],
            [0.011847e-6,5643.178563677,5.489005403],
            [0.008610e-6,3340.612426700,3.661698944],
            [0.011622e-6,5120.601145584,4.863931876],
            [0.010825e-6,553.569402842,0.842715011],
            [0.008666e-6,-135.065080035,3.293406547],
            [0.009963e-6,149.563197135,4.870690598],
            [0.009858e-6,6309.374169791,1.061816410],
            [0.007959e-6,316.391869657,2.465042647],
            [0.010099e-6,283.859318865,1.942176992],
            [0.007147e-6,-242.728603974,3.661486981],
            [0.007505e-6,5230.807466803,4.920937029],
            [0.008323e-6,11769.853693166,1.229392026],
            [0.007490e-6,-6256.777530192,3.658444681],
            [0.009370e-6,149854.400134205,0.673880395],
            [0.007117e-6,38.027672636,5.294249518],
            [0.007857e-6,12168.002696575,0.525733528],
            [0.007019e-6,6206.809778716,0.837688810],
            [0.006056e-6,955.599741609,4.194535082],
            [0.008107e-6,13367.972631107,3.793235253],
            [0.006731e-6,5650.292110678,5.639906583],
            [0.007332e-6,36.648562930,0.114858677],
            [0.006366e-6,4164.311989613,2.262081818],
            [0.006858e-6,5216.580372801,0.642063318],
            [0.006919e-6,6681.224853400,6.018501522],
            [0.006826e-6,7632.943259650,3.458654112],
            [0.005308e-6,-1592.596013633,2.500382359],
            [0.005096e-6,11371.704689758,2.547107806],
            [0.004841e-6,5333.900241022,0.437078094],
            [0.005582e-6,5966.683980335,2.246174308],
            [0.006304e-6,11926.254413669,2.512929171],
            [0.006603e-6,23581.258177318,5.393136889],
            [0.005123e-6,-1.484472708,2.999641028],
            [0.004648e-6,1589.072895284,1.275847090],
            [0.005119e-6,6438.496249426,1.486539246],
            [0.004521e-6,4292.330832950,6.140635794],
            [0.005680e-6,23013.539539587,4.557814849],
            [0.005488e-6,-3.455808046,0.090675389],
            [0.004193e-6,7234.794256242,4.869091389],
            [0.003742e-6,7238.675591600,4.691976180],
            [0.004148e-6,-110.206321219,3.016173439],
            [0.004553e-6,11499.656222793,5.554998314],
            [0.004892e-6,5436.993015240,1.475415597],
            [0.004044e-6,4732.030627343,1.398784824],
            [0.004164e-6,12491.370101415,5.650931916],
            [0.004349e-6,11513.883316794,2.181745369],
            [0.003919e-6,12528.018664345,5.823319737],
            [0.003129e-6,6836.645252834,0.003844094],
            [0.004080e-6,-7058.598461315,3.690360123],
            [0.003270e-6,76.266071276,1.517189902],
            [0.002954e-6,6283.143160294,4.447203799],
            [0.002872e-6,28.449187468,1.158692983],
            [0.002881e-6,735.876513532,0.349250250],
            [0.003279e-6,5849.364112115,4.893384368],
            [0.003625e-6,6209.778724132,1.473760578],
            [0.003074e-6,949.175608970,5.185878737],
            [0.002775e-6,9917.696874510,1.030026325],
            [0.002646e-6,10973.555686350,3.918259169],
            [0.002575e-6,25132.303399966,6.109659023],
            [0.003500e-6,263.083923373,1.892100742],
            [0.002740e-6,18319.536584880,4.320519510],
            [0.002464e-6,202.253395174,4.698203059],
            [0.002409e-6,2.542797281,5.325009315],
            [0.003354e-6,-90955.551694697,1.942656623],
            [0.002296e-6,6496.374945429,5.061810696],
            [0.003002e-6,6172.869528772,2.797822767],
            [0.003202e-6,27511.467873537,0.531673101],
            [0.002954e-6,-6283.008539689,4.533471191],
            [0.002353e-6,639.897286314,3.734548088],
            [0.002401e-6,16200.772724501,2.605547070],
            [0.003053e-6,233141.314403759,3.029030662],
            [0.003024e-6,83286.914269554,2.355556099],
            [0.002863e-6,17298.182327326,5.240963796],
            [0.002103e-6,-7079.373856808,5.756641637],
            [0.002303e-6,83996.847317911,2.013686814],
            [0.002303e-6,18073.704938650,1.089100410],
            [0.002381e-6,63.735898303,0.759188178],
            [0.002493e-6,6386.168624210,0.645026535],
            [0.002366e-6,3.932153263,6.215885448],
            [0.002169e-6,11015.106477335,4.845297676],
            [0.002397e-6,6243.458341645,3.809290043],
            [0.002183e-6,1162.474704408,6.179611691],
            [0.002353e-6,6246.427287062,4.781719760],
            [0.002199e-6,-245.831646229,5.956152284],
            [0.001729e-6,3894.181829542,1.264976635],
            [0.001896e-6,-3128.388765096,4.914231596],
            [0.002085e-6,35.164090221,1.405158503],
            [0.002024e-6,14712.317116458,2.752035928],
            [0.001737e-6,6290.189396992,5.280820144],
            [0.002229e-6,491.557929457,1.571007057],
            [0.001602e-6,14314.168113050,4.203664806],
            [0.002186e-6,454.909366527,1.402101526],
            [0.001897e-6,22483.848574493,4.167932508],
            [0.001825e-6,-3738.761430108,0.545828785],
            [0.001894e-6,1052.268383188,5.817167450],
            [0.001421e-6,20.355319399,2.419886601],
            [0.001408e-6,10984.192351700,2.732084787],
            [0.001847e-6,10873.986030480,2.903477885],
            [0.001391e-6,-8635.942003763,0.593891500],
            [0.001388e-6,-7.046236698,1.166145902],
            [0.001810e-6,-88860.057071188,0.487355242],
            [0.001288e-6,-1990.745017041,3.913022880],
            [0.001297e-6,23543.230504682,3.063805171],
            [0.001335e-6,-266.607041722,3.995764039],
            [0.001376e-6,10969.965257698,5.152914309],
            [0.001745e-6,244287.600007027,3.626395673],
            [0.001649e-6,31441.677569757,1.952049260],
            [0.001416e-6,9225.539273283,4.996408389],
            [0.001238e-6,4804.209275927,5.503379738],
            [0.001472e-6,4590.910180489,4.164913291],
            [0.001169e-6,6040.347246017,5.841719038],
            [0.001039e-6,5540.085789459,2.769753519],
            [0.001004e-6,-170.672870619,0.755008103],
            [0.001284e-6,10575.406682942,5.306538209],
            [0.001278e-6,71.812653151,4.713486491],
            [0.001321e-6,18209.330263660,2.624866359],
            [0.001297e-6,21228.392023546,0.382603541],
            [0.000954e-6,6282.095528923,0.882213514],
            [0.001145e-6,6058.731054289,1.169483931],
            [0.000979e-6,5547.199336460,5.448375984],
            [0.000987e-6,-6262.300454499,2.656486959],
            [0.001070e-6,-154717.609887482,1.827624012],
            [0.000991e-6,4701.116501708,4.387001801],
            [0.001155e-6,-14.227094002,3.042700750],
            [0.001176e-6,277.034993741,3.335519004],
            [0.000890e-6,13916.019109642,5.601498297],
            [0.000884e-6,-1551.045222648,1.088831705],
            [0.000876e-6,5017.508371365,3.969902609],
            [0.000806e-6,15110.466119866,5.142876744],
            [0.000773e-6,-4136.910433516,0.022067765],
            [0.001077e-6,175.166059800,1.844913056],
            [0.000954e-6,-6284.056171060,0.968480906],
            [0.000737e-6,5326.786694021,4.923831588],
            [0.000845e-6,-433.711737877,4.749245231],
            [0.000819e-6,8662.240323563,5.991247817],
            [0.000852e-6,199.072001436,2.189604979],
            [0.000723e-6,17256.631536341,6.068719637],
            [0.000940e-6,6037.244203762,6.197428148],
            [0.000885e-6,11712.955318231,3.280414875],
            [0.000706e-6,12559.038152982,2.824848947],
            [0.000732e-6,2379.164473572,2.501813417],
            [0.000764e-6,-6127.655450557,2.236346329],
            [0.000908e-6,131.541961686,2.521257490],
            [0.000907e-6,35371.887265976,3.370195967],
            [0.000673e-6,1066.495477190,3.876512374],
            [0.000814e-6,17654.780539750,4.627122566],
            [0.000630e-6,36.027866677,0.156368499],
            [0.000798e-6,515.463871093,5.151962502],
            [0.000798e-6,148.078724426,5.909225055],
            [0.000806e-6,309.278322656,6.054064447],
            [0.000607e-6,-39.617508346,2.839021623],
            [0.000601e-6,412.371096874,3.984225404],
            [0.000646e-6,11403.676995575,3.852959484],
            [0.000704e-6,13521.751441591,2.300991267],
            [0.000603e-6,-65147.619767937,4.140083146],
            [0.000609e-6,10177.257679534,0.437122327],
            [0.000631e-6,5767.611978898,4.026532329],
            [0.000576e-6,11087.285125918,4.760293101],
            [0.000674e-6,14945.316173554,6.270510511],
            [0.000726e-6,5429.879468239,6.039606892],
            [0.000710e-6,28766.924424484,5.672617711],
            [0.000647e-6,11856.218651625,3.397132627],
            [0.000678e-6,-5481.254918868,6.249666675],
            [0.000618e-6,22003.914634870,2.466427018],
            [0.000738e-6,6134.997125565,2.242668890],
            [0.000660e-6,625.670192312,5.864091907],
            [0.000694e-6,3496.032826134,2.668309141],
            [0.000531e-6,6489.261398429,1.681888780],
            [0.000611e-6,-143571.324284214,2.424978312],
            [0.000575e-6,12043.574281889,4.216492400],
            [0.000553e-6,12416.588502848,4.772158039],
            [0.000689e-6,4686.889407707,6.224271088],
            [0.000495e-6,7342.457780181,3.817285811],
            [0.000567e-6,3634.621024518,1.649264690],
            [0.000515e-6,18635.928454536,3.945345892],
            [0.000486e-6,-323.505416657,4.061673868],
            [0.000662e-6,25158.601719765,1.794058369],
            [0.000509e-6,846.082834751,3.053874588],
            [0.000472e-6,-12569.674818332,5.112133338],
            [0.000461e-6,6179.983075773,0.513669325],
            [0.000641e-6,83467.156352816,3.210727723],
            [0.000520e-6,10344.295065386,2.445597761],
            [0.000493e-6,18422.629359098,1.676939306],
            [0.000478e-6,1265.567478626,5.487314569],
            [0.000472e-6,-18.159247265,1.999707589],
            [0.000559e-6,11190.377900137,5.783236356],
            [0.000494e-6,9623.688276691,3.022645053],
            [0.000463e-6,5739.157790895,1.411223013],
            [0.000432e-6,16858.482532933,1.179256434],
            [0.000574e-6,72140.628666286,1.758191830],
            [0.000484e-6,17267.268201691,3.290589143],
            [0.000550e-6,4907.302050146,0.864024298],
            [0.000399e-6,14.977853527,2.094441910],
            [0.000491e-6,224.344795702,0.878372791],
            [0.000432e-6,20426.571092422,6.003829241],
            [0.000481e-6,5749.452731634,4.309591964],
            [0.000480e-6,5757.317038160,1.142348571],
            [0.000485e-6,6702.560493867,0.210580917],
            [0.000426e-6,6055.549660552,4.274476529],
            [0.000480e-6,5959.570433334,5.031351030],
            [0.000466e-6,12562.628581634,4.959581597],
            [0.000520e-6,39302.096962196,4.788002889],
            [0.000458e-6,12132.439962106,1.880103788],
            [0.000470e-6,12029.347187887,1.405611197],
            [0.000416e-6,-7477.522860216,1.082356330],
            [0.000449e-6,11609.862544012,4.179989585],
            [0.000465e-6,17253.041107690,0.353496295],
            [0.000362e-6,-4535.059436924,1.583849576],
            [0.000383e-6,21954.157609398,3.747376371],
            [0.000389e-6,17.252277143,1.395753179],
            [0.000331e-6,18052.929543158,0.566790582],
            [0.000430e-6,13517.870106233,0.685827538],
            [0.000368e-6,-5756.908003246,0.731374317],
            [0.000330e-6,10557.594160824,3.710043680],
            [0.000332e-6,20199.094959633,1.652901407],
            [0.000384e-6,11933.367960670,5.827781531],
            [0.000387e-6,10454.501386605,2.541182564],
            [0.000325e-6,15671.081759407,2.178850542],
            [0.000318e-6,138.517496871,2.253253037],
            [0.000305e-6,9388.005909415,0.578340206],
            [0.000352e-6,5749.861766548,3.000297967],
            [0.000311e-6,6915.859589305,1.693574249],
            [0.000297e-6,24072.921469776,1.997249392],
            [0.000363e-6,-640.877607382,5.071820966],
            [0.000323e-6,12592.450019783,1.072262823],
            [0.000341e-6,12146.667056108,4.700657997],
            [0.000290e-6,9779.108676125,1.812320441],
            [0.000342e-6,6132.028180148,4.322238614],
            [0.000329e-6,6268.848755990,3.033827743],
            [0.000374e-6,17996.031168222,3.388716544],
            [0.000285e-6,-533.214083444,4.687313233],
            [0.000338e-6,6065.844601290,0.877776108],
            [0.000276e-6,24.298513841,0.770299429],
            [0.000336e-6,-2388.894020449,5.353796034],
            [0.000290e-6,3097.883822726,4.075291557],
            [0.000318e-6,709.933048357,5.941207518],
            [0.000271e-6,13095.842665077,3.208912203],
            [0.000331e-6,6073.708907816,4.007881169],
            [0.000292e-6,742.990060533,2.714333592],
            [0.000362e-6,29088.811415985,3.215977013],
            [0.000280e-6,12359.966151546,0.710872502],
            [0.000267e-6,10440.274292604,4.730108488],
            [0.000262e-6,838.969287750,1.327720272],
            [0.000250e-6,16496.361396202,0.898769761],
            [0.000325e-6,20597.243963041,0.180044365],
            [0.000268e-6,6148.010769956,5.152666276],
            [0.000284e-6,5636.065016677,5.655385808],
            [0.000301e-6,6080.822454817,2.135396205],
            [0.000294e-6,-377.373607916,3.708784168],
            [0.000236e-6,2118.763860378,1.733578756],
            [0.000234e-6,5867.523359379,5.575209112],
            [0.000268e-6,-226858.238553767,0.069432392],
            [0.000265e-6,167283.761587465,4.369302826],
            [0.000280e-6,28237.233459389,5.304829118],
            [0.000292e-6,12345.739057544,4.096094132],
            [0.000223e-6,19800.945956225,3.069327406],
            [0.000301e-6,43232.306658416,6.205311188],
            [0.000264e-6,18875.525869774,1.417263408],
            [0.000304e-6,-1823.175188677,3.409035232],
            [0.000301e-6,109.945688789,0.510922054],
            [0.000260e-6,813.550283960,2.389438934],
            [0.000299e-6,316428.228673312,5.384595078],
            [0.000211e-6,5756.566278634,3.789392838],
            [0.000209e-6,5750.203491159,1.661943545],
            [0.000240e-6,12489.885628707,5.684549045],
            [0.000216e-6,6303.851245484,3.862942261],
            [0.000203e-6,1581.959348283,5.549853589],
            [0.000200e-6,5642.198242609,1.016115785],
            [0.000197e-6,-70.849445304,4.690702525],
            [0.000227e-6,6287.008003254,2.911891613],
            [0.000197e-6,533.623118358,1.048982898],
            [0.000205e-6,-6279.485421340,1.829362730],
            [0.000209e-6,-10988.808157535,2.636140084],
            [0.000208e-6,-227.526189440,4.127883842],
            [0.000191e-6,415.552490612,4.401165650],
            [0.000190e-6,29296.615389579,4.175658539],
            [0.000264e-6,66567.485864652,4.601102551],
            [0.000256e-6,-3646.350377354,0.506364778],
            [0.000188e-6,13119.721102825,2.032195842],
            [0.000185e-6,-209.366942175,4.694756586],
            [0.000198e-6,25934.124331089,3.832703118],
            [0.000195e-6,4061.219215394,3.308463427],
            [0.000234e-6,5113.487598583,1.716090661],
            [0.000188e-6,1478.866574064,5.686865780],
            [0.000222e-6,11823.161639450,1.942386641],
            [0.000181e-6,10770.893256262,1.999482059],
            [0.000171e-6,6546.159773364,1.182807992],
            [0.000206e-6,70.328180442,5.934076062],
            [0.000169e-6,20995.392966449,2.169080622],
            [0.000191e-6,10660.686935042,5.405515999],
            [0.000228e-6,33019.021112205,4.656985514],
            [0.000184e-6,-4933.208440333,3.327476868],
            [0.000220e-6,-135.625325010,1.765430262],
            [0.000166e-6,23141.558382925,3.454132746],
            [0.000191e-6,6144.558353121,5.020393445],
            [0.000180e-6,6084.003848555,0.602182191],
            [0.000163e-6,17782.732072784,4.960593133],
            [0.000225e-6,16460.333529525,2.596451817],
            [0.000222e-6,5905.702242076,3.731990323],
            [0.000204e-6,227.476132789,5.636192701],
            [0.000159e-6,16737.577236597,3.600691544],
            [0.000200e-6,6805.653268085,0.868220961],
            [0.000187e-6,11919.140866668,2.629456641],
            [0.000161e-6,127.471796607,2.862574720],
            [0.000205e-6,6286.666278643,1.742882331],
            [0.000189e-6,153.778810485,4.812372643],
            [0.000168e-6,16723.350142595,0.027860588],
            [0.000149e-6,11720.068865232,0.659721876],
            [0.000189e-6,5237.921013804,5.245313000],
            [0.000143e-6,6709.674040867,4.317625647],
            [0.000146e-6,4487.817406270,4.815297007],
            [0.000144e-6,-664.756045130,5.381366880],
            [0.000175e-6,5127.714692584,4.728443327],
            [0.000162e-6,6254.626662524,1.435132069],
            [0.000187e-6,47162.516354635,1.354371923],
            [0.000146e-6,11080.171578918,3.369695406],
            [0.000180e-6,-348.924420448,2.490902145],
            [0.000148e-6,151.047669843,3.799109588],
            [0.000157e-6,6197.248551160,1.284375887],
            [0.000167e-6,146.594251718,0.759969109],
            [0.000133e-6,-5331.357443741,5.409701889],
            [0.000154e-6,95.979227218,3.366890614],
            [0.000148e-6,-6418.140930027,3.384104996],
            [0.000128e-6,-6525.804453965,3.803419985],
            [0.000130e-6,11293.470674356,0.939039445],
            [0.000152e-6,-5729.506447149,0.734117523],
            [0.000138e-6,210.117701700,2.564216078],
            [0.000123e-6,6066.595360816,4.517099537],
            [0.000140e-6,18451.078546566,0.642049130],
            [0.000126e-6,11300.584221356,3.485280663],
            [0.000119e-6,10027.903195729,3.217431161],
            [0.000151e-6,4274.518310832,4.404359108],
            [0.000117e-6,6072.958148291,0.366324650],
            [0.000165e-6,-7668.637425143,4.298212528],
            [0.000117e-6,-6245.048177356,5.379518958],
            [0.000130e-6,-5888.449964932,4.527681115],
            [0.000121e-6,-543.918059096,6.109429504],
            [0.000162e-6,9683.594581116,5.720092446],
            [0.000141e-6,6219.339951688,0.679068671],
            [0.000118e-6,22743.409379516,4.881123092],
            [0.000129e-6,1692.165669502,0.351407289],
            [0.000126e-6,5657.405657679,5.146592349],
            [0.000114e-6,728.762966531,0.520791814],
            [0.000120e-6,52.596639600,0.948516300],
            [0.000115e-6,65.220371012,3.504914846],
            [0.000126e-6,5881.403728234,5.577502482],
            [0.000158e-6,163096.180360983,2.957128968],
            [0.000134e-6,12341.806904281,2.598576764],
            [0.000151e-6,16627.370915377,3.985702050],
            [0.000109e-6,1368.660252845,0.014730471],
            [0.000131e-6,6211.263196841,0.085077024],
            [0.000146e-6,5792.741760812,0.708426604],
            [0.000146e-6,-77.750543984,3.121576600],
            [0.000107e-6,5341.013788022,0.288231904],
            [0.000138e-6,6281.591377283,2.797450317],
            [0.000113e-6,-6277.552925684,2.788904128],
            [0.000115e-6,-525.758811831,5.895222200],
            [0.000138e-6,6016.468808270,6.096188999],
            [0.000139e-6,23539.707386333,2.028195445],
            [0.000146e-6,-4176.041342449,4.660008502],
            [0.000107e-6,16062.184526117,4.066520001],
            [0.000142e-6,83783.548222473,2.936315115],
            [0.000128e-6,9380.959672717,3.223844306],
            [0.000135e-6,6205.325306007,1.638054048],
            [0.000101e-6,2699.734819318,5.481603249],
            [0.000104e-6,-568.821874027,2.205734493],
            [0.000103e-6,6321.103522627,2.440421099],
            [0.000119e-6,6321.208885629,2.547496264],
            [0.000138e-6,1975.492545856,2.314608466],
            [0.000121e-6,137.033024162,4.539108237],
            [0.000123e-6,19402.796952817,4.538074405],
            [0.000119e-6,22805.735565994,2.869040566],
            [0.000133e-6,64471.991241142,6.056405489],
            [0.000129e-6,-85.827298831,2.540635083],
            [0.000131e-6,13613.804277336,4.005732868],
            [0.000104e-6,9814.604100291,1.959967212],
            [0.000112e-6,16097.679950283,3.589026260],
            [0.000123e-6,2107.034507542,1.728627253],
            [0.000121e-6,36949.230808424,6.072332087],
            [0.000108e-6,-12539.853380183,3.716133846],
            [0.000113e-6,-7875.671863624,2.725771122],
            [0.000109e-6,4171.425536614,4.033338079],
            [0.000101e-6,6247.911759770,3.441347021],
            [0.000113e-6,7330.728427345,0.656372122],
            [0.000113e-6,51092.726050855,2.791483066],
            [0.000106e-6,5621.842923210,1.815323326],
            [0.000101e-6,111.430161497,5.711033677],
            [0.000103e-6,909.818733055,2.812745443],
            [0.000101e-6,1790.642637886,1.965746028],
            #T
            [102.156724e-6,6283.075849991,4.249032005],
            [1.706807e-6,12566.151699983,4.205904248],
            [0.269668e-6,213.299095438,3.400290479],
            [0.265919e-6,529.690965095,5.836047367],
            [0.210568e-6,-3.523118349,6.262738348],
            [0.077996e-6,5223.693919802,4.670344204],
            [0.054764e-6,1577.343542448,4.534800170],
            [0.059146e-6,26.298319800,1.083044735],
            [0.034420e-6,-398.149003408,5.980077351],
            [0.032088e-6,18849.227549974,4.162913471],
            [0.033595e-6,5507.553238667,5.980162321],
            [0.029198e-6,5856.477659115,0.623811863],
            [0.027764e-6,155.420399434,3.745318113],
            [0.025190e-6,5746.271337896,2.980330535],
            [0.022997e-6,-796.298006816,1.174411803],
            [0.024976e-6,5760.498431898,2.467913690],
            [0.021774e-6,206.185548437,3.854787540],
            [0.017925e-6,-775.522611324,1.092065955],
            [0.013794e-6,426.598190876,2.699831988],
            [0.013276e-6,6062.663207553,5.845801920],
            [0.011774e-6,12036.460734888,2.292832062],
            [0.012869e-6,6076.890301554,5.333425680],
            [0.012152e-6,1059.381930189,6.222874454],
            [0.011081e-6,-7.113547001,5.154724984],
            [0.010143e-6,4694.002954708,4.044013795],
            [0.009357e-6,5486.777843175,3.416081409],
            [0.010084e-6,522.577418094,0.749320262],
            [0.008587e-6,10977.078804699,2.777152598],
            [0.008628e-6,6275.962302991,4.562060226],
            [0.008158e-6,-220.412642439,5.806891533],
            [0.007746e-6,2544.314419883,1.603197066],
            [0.007670e-6,2146.165416475,3.000200440],
            [0.007098e-6,74.781598567,0.443725817],
            [0.006180e-6,-536.804512095,1.302642751],
            [0.005818e-6,5088.628839767,4.827723531],
            [0.004945e-6,-6286.598968340,0.268305170],
            [0.004774e-6,1349.867409659,5.808636673],
            [0.004687e-6,-242.728603974,5.154890570],
            [0.006089e-6,1748.016413067,4.403765209],
            [0.005975e-6,-1194.447010225,2.583472591],
            [0.004229e-6,951.718406251,0.931172179],
            [0.005264e-6,553.569402842,2.336107252],
            [0.003049e-6,5643.178563677,1.362634430],
            [0.002974e-6,6812.766815086,1.583012668],
            [0.003403e-6,-2352.866153772,2.552189886],
            [0.003030e-6,419.484643875,5.286473844],
            [0.003210e-6,-7.046236698,1.863796539],
            [0.003058e-6,9437.762934887,4.226420633],
            [0.002589e-6,12352.852604545,1.991935820],
            [0.002927e-6,5216.580372801,2.319951253],
            [0.002425e-6,5230.807466803,3.084752833],
            [0.002656e-6,3154.687084896,2.487447866],
            [0.002445e-6,10447.387839604,2.347139160],
            [0.002990e-6,4690.479836359,6.235872050],
            [0.002890e-6,5863.591206116,0.095197563],
            [0.002498e-6,6438.496249426,2.994779800],
            [0.001889e-6,8031.092263058,3.569003717],
            [0.002567e-6,801.820931124,3.425611498],
            [0.001803e-6,-71430.695617928,2.192295512],
            [0.001782e-6,3.932153263,5.180433689],
            [0.001694e-6,-4705.732307544,4.641779174],
            [0.001704e-6,-1592.596013633,3.997097652],
            [0.001735e-6,5849.364112115,0.417558428],
            [0.001643e-6,8429.241266467,2.180619584],
            [0.001680e-6,38.133035638,4.164529426],
            [0.002045e-6,7084.896781115,0.526323854],
            [0.001458e-6,4292.330832950,1.356098141],
            [0.001437e-6,20.355319399,3.895439360],
            [0.001738e-6,6279.552731642,0.087484036],
            [0.001367e-6,14143.495242431,3.987576591],
            [0.001344e-6,7234.794256242,0.090454338],
            [0.001438e-6,11499.656222793,0.974387904],
            [0.001257e-6,6836.645252834,1.509069366],
            [0.001358e-6,11513.883316794,0.495572260],
            [0.001628e-6,7632.943259650,4.968445721],
            [0.001169e-6,103.092774219,2.838496795],
            [0.001162e-6,4164.311989613,3.408387778],
            [0.001092e-6,6069.776754553,3.617942651],
            [0.001008e-6,17789.845619785,0.286350174],
            [0.001008e-6,639.897286314,1.610762073],
            [0.000918e-6,10213.285546211,5.532798067],
            [0.001011e-6,-6256.777530192,0.661826484],
            [0.000753e-6,16730.463689596,3.905030235],
            [0.000737e-6,11926.254413669,4.641956361],
            [0.000694e-6,3340.612426700,2.111120332],
            [0.000701e-6,3894.181829542,2.760823491],
            [0.000689e-6,-135.065080035,4.768800780],
            [0.000700e-6,13367.972631107,5.760439898],
            [0.000664e-6,6040.347246017,1.051215840],
            [0.000654e-6,5650.292110678,4.911332503],
            [0.000788e-6,6681.224853400,4.699648011],
            [0.000628e-6,5333.900241022,5.024608847],
            [0.000755e-6,-110.206321219,4.370971253],
            [0.000628e-6,6290.189396992,3.660478857],
            [0.000635e-6,25132.303399966,4.121051532],
            [0.000534e-6,5966.683980335,1.173284524],
            [0.000543e-6,-433.711737877,0.345585464],
            [0.000517e-6,-1990.745017041,5.414571768],
            [0.000504e-6,5767.611978898,2.328281115],
            [0.000485e-6,5753.384884897,1.685874771],
            [0.000463e-6,7860.419392439,5.297703006],
            [0.000604e-6,515.463871093,0.591998446],
            [0.000443e-6,12168.002696575,4.830881244],
            [0.000570e-6,199.072001436,3.899190272],
            [0.000465e-6,10969.965257698,0.476681802],
            [0.000424e-6,-7079.373856808,1.112242763],
            [0.000427e-6,735.876513532,1.994214480],
            [0.000478e-6,-6127.655450557,3.778025483],
            [0.000414e-6,10973.555686350,5.441088327],
            [0.000512e-6,1589.072895284,0.107123853],
            [0.000378e-6,10984.192351700,0.915087231],
            [0.000402e-6,11371.704689758,4.107281715],
            [0.000453e-6,9917.696874510,1.917490952],
            [0.000395e-6,149.563197135,2.763124165],
            [0.000371e-6,5739.157790895,3.112111866],
            [0.000350e-6,11790.629088659,0.440639857],
            [0.000356e-6,6133.512652857,5.444568842],
            [0.000344e-6,412.371096874,5.676832684],
            [0.000383e-6,955.599741609,5.559734846],
            [0.000333e-6,6496.374945429,0.261537984],
            [0.000340e-6,6055.549660552,5.975534987],
            [0.000334e-6,1066.495477190,2.335063907],
            [0.000399e-6,11506.769769794,5.321230910],
            [0.000314e-6,18319.536584880,2.313312404],
            [0.000424e-6,1052.268383188,1.211961766],
            [0.000307e-6,63.735898303,3.169551388],
            [0.000329e-6,29.821438149,6.106912080],
            [0.000357e-6,6309.374169791,4.223760346],
            [0.000312e-6,-3738.761430108,2.180556645],
            [0.000301e-6,309.278322656,1.499984572],
            [0.000268e-6,12043.574281889,2.447520648],
            [0.000257e-6,12491.370101415,3.662331761],
            [0.000290e-6,625.670192312,1.272834584],
            [0.000256e-6,5429.879468239,1.913426912],
            [0.000339e-6,3496.032826134,4.165930011],
            [0.000283e-6,3930.209696220,4.325565754],
            [0.000241e-6,12528.018664345,3.832324536],
            [0.000304e-6,4686.889407707,1.612348468],
            [0.000259e-6,16200.772724501,3.470173146],
            [0.000238e-6,12139.553509107,1.147977842],
            [0.000236e-6,6172.869528772,3.776271728],
            [0.000296e-6,-7058.598461315,0.460368852],
            [0.000306e-6,10575.406682942,0.554749016],
            [0.000251e-6,17298.182327326,0.834332510],
            [0.000290e-6,4732.030627343,4.759564091],
            [0.000261e-6,5884.926846583,0.298259862],
            [0.000249e-6,5547.199336460,3.749366406],
            [0.000213e-6,11712.955318231,5.415666119],
            [0.000223e-6,4701.116501708,2.703203558],
            [0.000268e-6,-640.877607382,0.283670793],
            [0.000209e-6,5636.065016677,1.238477199],
            [0.000193e-6,10177.257679534,1.943251340],
            [0.000182e-6,6283.143160294,2.456157599],
            [0.000184e-6,-227.526189440,5.888038582],
            [0.000182e-6,-6283.008539689,0.241332086],
            [0.000228e-6,-6284.056171060,2.657323816],
            [0.000166e-6,7238.675591600,5.930629110],
            [0.000167e-6,3097.883822726,5.570955333],
            [0.000159e-6,-323.505416657,5.786670700],
            [0.000154e-6,-4136.910433516,1.517805532],
            [0.000176e-6,12029.347187887,3.139266834],
            [0.000167e-6,12132.439962106,3.556352289],
            [0.000153e-6,202.253395174,1.463313961],
            [0.000157e-6,17267.268201691,1.586837396],
            [0.000142e-6,83996.847317911,0.022670115],
            [0.000152e-6,17260.154654690,0.708528947],
            [0.000144e-6,6084.003848555,5.187075177],
            [0.000135e-6,5756.566278634,1.993229262],
            [0.000134e-6,5750.203491159,3.457197134],
            [0.000144e-6,5326.786694021,6.066193291],
            [0.000160e-6,11015.106477335,1.710431974],
            [0.000133e-6,3634.621024518,2.836451652],
            [0.000134e-6,18073.704938650,5.453106665],
            [0.000134e-6,1162.474704408,5.326898811],
            [0.000128e-6,5642.198242609,2.511652591],
            [0.000160e-6,632.783739313,5.628785365],
            [0.000132e-6,13916.019109642,0.819294053],
            [0.000122e-6,14314.168113050,5.677408071],
            [0.000125e-6,12359.966151546,5.251984735],
            [0.000121e-6,5749.452731634,2.210924603],
            [0.000136e-6,-245.831646229,1.646502367],
            [0.000120e-6,5757.317038160,3.240883049],
            [0.000134e-6,12146.667056108,3.059480037],
            [0.000137e-6,6206.809778716,1.867105418],
            [0.000141e-6,17253.041107690,2.069217456],
            [0.000129e-6,-7477.522860216,2.781469314],
            [0.000116e-6,5540.085789459,4.281176991],
            [0.000116e-6,9779.108676125,3.320925381],
            [0.000129e-6,5237.921013804,3.497704076],
            [0.000113e-6,5959.570433334,0.983210840],
            [0.000122e-6,6282.095528923,2.674938860],
            [0.000140e-6,-11.045700264,4.957936982],
            [0.000108e-6,23543.230504682,1.390113589],
            [0.000106e-6,-12569.674818332,0.429631317],
            [0.000110e-6,-266.607041722,5.501340197],
            [0.000115e-6,12559.038152982,4.691456618],
            [0.000134e-6,-2388.894020449,0.577313584],
            [0.000109e-6,10440.274292604,6.218148717],
            [0.000102e-6,-543.918059096,1.477842615],
            [0.000108e-6,21228.392023546,2.237753948],
            [0.000101e-6,-4535.059436924,3.100492232],
            [0.000103e-6,76.266071276,5.594294322],
            [0.000104e-6,949.175608970,5.674287810],
            [0.000101e-6,13517.870106233,2.196632348],
            [0.000100e-6,11933.367960670,4.056084160],
            #T^2
            [4.322990e-6,6283.075849991,2.642893748],
            [0.406495e-6,0.000000000,4.712388980],
            [0.122605e-6,12566.151699983,2.438140634],
            [0.019476e-6,213.299095438,1.642186981],
            [0.016916e-6,529.690965095,4.510959344],
            [0.013374e-6,-3.523118349,1.502210314],
            [0.008042e-6,26.298319800,0.478549024],
            [0.007824e-6,155.420399434,5.254710405],
            [0.004894e-6,5746.271337896,4.683210850],
            [0.004875e-6,5760.498431898,0.759507698],
            [0.004416e-6,5223.693919802,6.028853166],
            [0.004088e-6,-7.113547001,0.060926389],
            [0.004433e-6,77713.771467920,3.627734103],
            [0.003277e-6,18849.227549974,2.327912542],
            [0.002703e-6,6062.663207553,1.271941729],
            [0.003435e-6,-775.522611324,0.747446224],
            [0.002618e-6,6076.890301554,3.633715689],
            [0.003146e-6,206.185548437,5.647874613],
            [0.002544e-6,1577.343542448,6.232904270],
            [0.002218e-6,-220.412642439,1.309509946],
            [0.002197e-6,5856.477659115,2.407212349],
            [0.002897e-6,5753.384884897,5.863842246],
            [0.001766e-6,426.598190876,0.754113147],
            [0.001738e-6,-796.298006816,2.714942671],
            [0.001695e-6,522.577418094,2.629369842],
            [0.001584e-6,5507.553238667,1.341138229],
            [0.001503e-6,-242.728603974,0.377699736],
            [0.001552e-6,-536.804512095,2.904684667],
            [0.001370e-6,-398.149003408,1.265599125],
            [0.001889e-6,-5573.142801634,4.413514859],
            [0.001722e-6,6069.776754553,2.445966339],
            [0.001124e-6,1059.381930189,5.041799657],
            [0.001258e-6,553.569402842,3.849557278],
            [0.000831e-6,951.718406251,2.471094709],
            [0.000767e-6,4694.002954708,5.363125422],
            [0.000756e-6,1349.867409659,1.046195744],
            [0.000775e-6,-11.045700264,0.245548001],
            [0.000597e-6,2146.165416475,4.543268798],
            [0.000568e-6,5216.580372801,4.178853144],
            [0.000711e-6,1748.016413067,5.934271972],
            [0.000499e-6,12036.460734888,0.624434410],
            [0.000671e-6,-1194.447010225,4.136047594],
            [0.000488e-6,5849.364112115,2.209679987],
            [0.000621e-6,6438.496249426,4.518860804],
            [0.000495e-6,-6286.598968340,1.868201275],
            [0.000456e-6,5230.807466803,1.271231591],
            [0.000451e-6,5088.628839767,0.084060889],
            [0.000435e-6,5643.178563677,3.324456609],
            [0.000387e-6,10977.078804699,4.052488477],
            [0.000547e-6,161000.685737473,2.841633844],
            [0.000522e-6,3154.687084896,2.171979966],
            [0.000375e-6,5486.777843175,4.983027306],
            [0.000421e-6,5863.591206116,4.546432249],
            [0.000439e-6,7084.896781115,0.522967921],
            [0.000309e-6,2544.314419883,3.172606705],
            [0.000347e-6,4690.479836359,1.479586566],
            [0.000317e-6,801.820931124,3.553088096],
            [0.000262e-6,419.484643875,0.606635550],
            [0.000248e-6,6836.645252834,3.014082064],
            [0.000245e-6,-1592.596013633,5.519526220],
            [0.000225e-6,4292.330832950,2.877956536],
            [0.000214e-6,7234.794256242,1.605227587],
            [0.000205e-6,5767.611978898,0.625804796],
            [0.000180e-6,10447.387839604,3.499954526],
            [0.000229e-6,199.072001436,5.632304604],
            [0.000214e-6,639.897286314,5.960227667],
            [0.000175e-6,-433.711737877,2.162417992],
            [0.000209e-6,515.463871093,2.322150893],
            [0.000173e-6,6040.347246017,2.556183691],
            [0.000184e-6,6309.374169791,4.732296790],
            [0.000227e-6,149854.400134205,5.385812217],
            [0.000154e-6,8031.092263058,5.120720920],
            [0.000151e-6,5739.157790895,4.815000443],
            [0.000197e-6,7632.943259650,0.222827271],
            [0.000197e-6,74.781598567,3.910456770],
            [0.000138e-6,6055.549660552,1.397484253],
            [0.000149e-6,-6127.655450557,5.333727496],
            [0.000137e-6,3894.181829542,4.281749907],
            [0.000135e-6,9437.762934887,5.979971885],
            [0.000139e-6,-2352.866153772,4.715630782],
            [0.000142e-6,6812.766815086,0.513330157],
            [0.000120e-6,-4705.732307544,0.194160689],
            [0.000131e-6,-71430.695617928,0.000379226],
            [0.000124e-6,6279.552731642,2.122264908],
            [0.000108e-6,-6256.777530192,0.883445696],
            #T^3
            [0.143388e-6,6283.075849991,1.131453581],
            [0.006671e-6,12566.151699983,0.775148887],
            [0.001480e-6,155.420399434,0.480016880],
            [0.000934e-6,213.299095438,6.144453084],
            [0.000795e-6,529.690965095,2.941595619],
            [0.000673e-6,5746.271337896,0.120415406],
            [0.000672e-6,5760.498431898,5.317009738],
            [0.000389e-6,-220.412642439,3.090323467],
            [0.000373e-6,6062.663207553,3.003551964],
            [0.000360e-6,6076.890301554,1.918913041],
            [0.000316e-6,-21.340641002,5.545798121],
            [0.000315e-6,-242.728603974,1.884932563],
            [0.000278e-6,206.185548437,1.266254859],
            [0.000238e-6,-536.804512095,4.532664830],
            [0.000185e-6,522.577418094,4.578313856],
            [0.000245e-6,18849.227549974,0.587467082],
            [0.000180e-6,426.598190876,5.151178553],
            [0.000200e-6,553.569402842,5.355983739],
            [0.000141e-6,5223.693919802,1.336556009],
            [0.000104e-6,5856.477659115,4.239842759],
            #T^4
            [0.003826e-6,6283.075849991,5.705257275],
            [0.000303e-6,12566.151699983,5.407132842],
            [0.000209e-6,155.420399434,1.989815753]]

    #自J200.0经过的时间，单位：千儒略年
    T=((DATE1-DJ00)+DATE2)/DJM
    
    #=================
    #Topocentric terms（以地面上的某点为中心的）
    #=================

    #将UT转换为以弧度为单位的当地太阳时间
    TSOL=(UT%1.0)*D2PI+ELONG
    
    #基于以下理论:  Simon et al. 1994.

    #将时间与度和角秒之间的关系联系起来.
    W=T/3600.0

    #太阳的平赤经.
    ELSUN=(280.46645683+1296027711.03429*W)%360.0*DD2R

    #太阳的平近点角.
    EMSUN=(357.52910918+1295965810.481*W)%360.0*DD2R

    #日月之间平均距离
    D=(297.85019547+16029616012.090*W)%360.0*DD2R

    #木星的平赤经.
    ELJ=(34.35151874+109306899.89453*W)%360.0*DD2R

    #土星的平赤经.
    ELS=(50.07744430+44046398.47038*W)%360.0*DD2R
    
    #TOPOCENTRIC TERMS:  Moyer 1981 and Murray 1983.
    WT=+0.00029e-10*U*ma.sin(TSOL+ELSUN-ELS)\
        +0.00100e-10*U*ma.sin(TSOL-2.0*EMSUN)\
        +0.00133e-10*U*ma.sin(TSOL-D)\
        +0.00133e-10*U*ma.sin(TSOL+ELSUN-ELJ)\
        -0.00229e-10*U*ma.sin(TSOL+2.0*ELSUN+EMSUN)\
        -0.02200e-10*V*ma.cos(ELSUN+EMSUN)\
        +0.05312e-10*U*ma.sin(TSOL-EMSUN)\
        -0.13677e-10*U*ma.sin(TSOL + 2.0*ELSUN)\
        -1.31840e-10*V*ma.cos(ELSUN)\
        +3.17679e-10*U*ma.sin(TSOL)
    
    #=====================
    #Fairhead et al. model
    #=====================
    
    #T**0
    W0=0.0
    for J in range(473,-1,-1):
        W0=W0+FAIRHD[J][0]*ma.sin(FAIRHD[J][1]*T + FAIRHD[J][2])
    
    #T**1
    W1=0.0
    for J in range(678,473,-1):
        W1=W1+FAIRHD[J][0]*ma.sin(FAIRHD[J][1]*T + FAIRHD[J][2])
    
    #T**2
    W2=0.0
    for J in range(763,678,-1):
        W2=W2+FAIRHD[J][0]*ma.sin(FAIRHD[J][1]*T + FAIRHD[J][2])
    
    #T**3
    W3=0.0
    for J in range(783,763,-1):
        W3=W3+FAIRHD[J][0]*ma.sin(FAIRHD[J][1]*T + FAIRHD[J][2])
    
    #T**4
    W4=0.0
    for J in range(786,783,-1):
        W4=W4+FAIRHD[J][0]*ma.sin(FAIRHD[J][1]*T + FAIRHD[J][2])
    
    #T的不同次项结合
    WF=T*(T*(T*(T*W4+W3)+W2)+W1)+W0
    
    #调整使用JPL的行星质量而不是IAU.
    WJ=0.00065e-6*ma.sin(6069.776754*T+4.021194)+\
        0.00033e-6*ma.sin(213.299095*T+5.543132)+\
            (-0.00196e-6*ma.sin(6208.294251*T+5.696701))+\
                (-0.00173e-6*ma.sin(74.781599*T+2.435900))+\
                    0.03638e-6*T*T
    #以秒为单位的TDB-TT时间差
    DTDB=WT+WF+WJ
    
    return(DTDB)


def pymTttdb(TT1,TT2,DTR):
    '''
    Time scale transformation:  Terrestrial Time, TT, to Barycentric
    Dynamical Time, TDB.

    Parameters
    ----------
    tt1 : float
        TT as a 2-part Julian Date
    tt2 : float
        TT as a 2-part Julian Date
    dtr : float
        TDB-TT in seconds

    Returns
    -------
    tdb1 : float
        TDB as a 2-part Julian Date
    tdb2 : float
        TDB as a 2-part Julian Date
    J : ValueError
        0 : OK
    '''
    #一天的秒长
    D2S=86400.0
    
    DTRD=DTR/D2S
    if (np.abs(TT1)>np.abs(TT2)):
        TDB1=TT1
        TDB2=TT2+DTRD
    else:
        TDB1=TT1+DTRD
        TDB2=TT2
    
    J=0   
    return(TDB1,TDB2,J)


def pymTttai(TT1,TT2):
    '''
    Time scale transformation:  Terrestrial Time, TT, to International
    Atomic Time, TAI.

    Parameters
    ----------
    tt1 : float
        TT as a 2-part Julian Date
    tt2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    tai1 : float
        TAI as a 2-part Julian Date
    tai2 : float
        TAI as a 2-part Julian Date
    J : ValueError
        0 : OK
    '''
    #TT-TAI的时差，单位：天.
    DTAT = 32.184/86400.0
    
    if (np.abs(TT1)>np.abs(TT2)):
        TAI1=TT1
        TAI2=TT2-DTAT
    else:
        TAI1=TT1-DTAT
        TAI2=TT2
    
    J=0
    return(TAI1,TAI2,J)


def pymTtut1(TT1,TT2,DT):
    '''
    Time scale transformation:  Terrestrial Time, TT, to Universal Time,
    UT1.

    Parameters
    ----------
    tt1 : float
        TT as a 2-part Julian Date
    tt2 : float
        TT as a 2-part Julian Date
    dtr : float
        TT-UT1 in seconds

    Returns
    -------
    ut11 : float
        UT1 as a 2-part Julian Date
    ut12 : float
        UT1 as a 2-part Julian Date
    J : ValueError
        0 : OK
    '''
    #一天的秒长
    D2S=86400.0
    DTD=DT/D2S
    if (np.abs(TT1)>np.abs(TT2)):
        UT11=TT1
        UT12=TT2-DTD
    else:
        UT11=TT1-DTD
        UT12=TT2
    
    J=0    
    return(UT11,UT12,J)


def pymTf2a(S,IHOUR,IMIN,SEC):
    '''
    Convert hours, minutes, seconds to radians.

    Parameters
    ----------
    s : str  
        sign:  '-' = negative, otherwise positive      
    ihour : int    
        hours    
    imin : int    
        minutes    
    sec : float    
        seconds    
        
    Returns
    -------
    rad : float
        angle in radians
    J : ValueError
        1:'ihour outside range 0-23',
        2:'imin outside range 0-59',
        3:'sec outside range 0-59.999...'    
    '''
    #1秒对应的弧度
    DS2R=7.272205216643039903848712e-5
    
    J=0
    if (SEC<0.0)|(SEC>=60.0):
        J=3
    if (IMIN<0)|(IMIN>59):
        J=2
    if (IHOUR<0)|(IHOUR>23):
        J=1
        
    W=(60.0*(60.0*float(np.abs(IHOUR))+float(np.abs(IMIN)))+np.abs(SEC))*DS2R
    
    if (S=='-'):
        W=-W
        
    RAD=W
    
    return(RAD, J)


def pymTcgtt(TCG1,TCG2):
    '''
    Time scale transformation:  Geocentric Coordinate Time, TCG, to
    Terrestrial Time, TT.

    Parameters
    ----------
    tcg1 : float
        TCG as a 2-part Julian Date
    tcg2 : float
        TCG as a 2-part Julian Date

    Returns
    -------
    tt1 : float
        TT as a 2-part Julian Date
    tt2 : float
        TT as a 2-part Julian Date
    J : ValueError
        0 : OK
    '''  
    #约简儒略日零点
    DJM0=2400000.5
    
    #1977 Jan 1 00:00:32.184 TT, 儒略日期
    T77T=43144.0003725
    
    #L_G = 1 - dTT/dTCG，TT和TCG之间的关系
    ELG=6.969290134e-10
    
    #尽可能保留精度
    if (np.abs(TCG1)>(np.abs(TCG2))):
        TT1=TCG1
        TT2=TCG2-((TCG1-DJM0)+(TCG2-T77T))*ELG
    else:
        TT1=TCG1-((TCG2-DJM0)+(TCG1-T77T))*ELG
        TT2=TCG2
    
    J=0
    return(TT1,TT2,J)


def pymTdbtt(TDB1,TDB2,DTR):
    '''
    Time scale transformation:  Barycentric Dynamical Time, TDB, to
    Terrestrial Time, TT.

    Parameters
    ----------
    tdb1 : float
        TDB as a 2-part Julian Date    
    tdb2 : float
        TDB as a 2-part Julian Date    
    dtr : float
        TDB-TT in seconds

    Returns
    -------
    tt1 : float
        TT as a 2-part Julian Date    
    tt2 : float
        TT as a 2-part Julian Date    
    J : ValueError
        0 : OK
    '''
    #一天的秒长
    D2S=86400.0
    
    DTRD=DTR/D2S
    if (np.abs(TDB1)>np.abs(TDB2)):
        TT1=TDB1
        TT2=TDB2-DTRD
    else:
        TT1=TDB1-DTRD
        TT2=TDB2
    
    J=0  
    return(TT1,TT2,J)


def pymTdbtcb(TDB1,TDB2):
    '''
    Time scale transformation:  Barycentric Dynamical Time, TDB, to
    Barycentric Coordinate Time, TCB.

    Parameters
    ----------
    tdb1 : float
        TDB as a 2-part Julian Date
    tdb2 : float
        TDB as a 2-part Julian Date

    Returns
    -------
    tcb1 : float
        TCB as a 2-part Julian Date
    tcb2 : float
        TCB as a 2-part Julian Date
    J : ValueError
        0 : OK
    '''    
    #1977 Jan 1.0 TAI = 1977/1/1 00:00:32.184 TCB, 作为两参数的儒略日期
    T77TD=2443144.0
    T77TF=0.5003725
    
    #L_B， TDB0 (d)，TDB-TCB之间的关系
    ELB=1.550519768e-8
    TDB0=-6.55e-5/86400.0
    ELBB=ELB/(1.0-ELB)
    
    if (np.abs(TDB1)>np.abs(TDB2)):
        D=T77TD-TDB1
        F=TDB2-TDB0
        TCB1=TDB1
        TCB2=F-(D-(F-T77TF))*ELBB
    else:
        D=T77TD-TDB2
        F=TDB1-TDB0
        TCB1=F-(D-(F-T77TF))*ELBB
        TCB2=TDB2
    
    J=0
    return(TCB1,TCB2,J)


def pymTcbtdb(TCB1,TCB2):
    '''
    Time scale transformation:  Barycentric Coordinate Time, TCB, to
    Barycentric Dynamical Time, TDB.

    Parameters
    ----------
    tcb1 : float
        TCB as a 2-part Julian Date
    tcb2 : float
        TCB as a 2-part Julian Date

    Returns
    -------
    tdb1 : float
        TDB as a 2-part Julian Date
    tdb2 :float
        TDB as a 2-part Julian Date
    J : ValueError
        0 : OK
    '''
    #1977 Jan 1.0 TAI = 1977/1/1 00:00:32.184 TCB, 两参数的儒略日
    T77TD=2443144.0
    T77TF=0.5003725
    
    #L_B, and TDB0 (d)，TCB-TDB之间的关系
    ELB=1.550519768e-8
    TDB0=-6.55e-5/86400.0
    
    if (np.abs(TCB1)>np.abs(TCB2)):
        D=TCB1-T77TD
        TDB1=TCB1
        TDB2=TCB2+TDB0-(D+(TCB2-T77TF))*ELB
    else:
        D=TCB2-T77TD
        TDB1=TCB1+TDB0-(D+(TCB1-T77TF))*ELB
        TDB2=TCB2
    
    J=0
    return(TDB1,TDB2,J)


def pymUt1tai(UT11,UT12,DTA):
    '''
    Time scale transformation:  Universal Time, UT1, to International
    Atomic Time, TAI.

    Parameters
    ----------
    ut11 : float
        UT1 as a 2-part Julian Date
    ut12 : float
        UT1 as a 2-part Julian Date
    dta : float
        UT1-TAI in seconds

    Returns
    -------
    tai1 : float
        TAI as a 2-part Julian Date
    tai2 : float
        TAI as a 2-part Julian Date
    J : ValueError
        0 : OK
    '''
    #一天的秒长
    D2S=86400.0
       
    DTAD=DTA/D2S
    if (np.abs(UT11)>np.abs(UT12)):
        TAI1=UT11
        TAI2=UT12-DTAD
    else:
        TAI1=UT11-DTAD
        TAI2=UT12
    
    J=0
    return(TAI1,TAI2,J)


def pymUt1tt(UT11,UT12,DT):
    '''
    Time scale transformation:  Universal Time, UT1, to Terrestrial
    Time, TT.

    Parameters
    ----------
    ut11 : float
        UT1 as a 2-part Julian Date
    ut12 : float
        UT1 as a 2-part Julian Date
    dt : float
        TT-UT1 in seconds

    Returns
    -------
    tt1 : float
        TT as a 2-part Julian Date
    tt2 : float
        TT as a 2-part Julian Date
    J : ValueError
        0 : OK
    '''
    #一天的秒长
    D2S = 86400.0
    
    DTD=DT/D2S
    if (np.abs(UT11)>np.abs(UT12)):
        TT1=UT11
        TT2=UT12+DTD
    else:
        TT1=UT11+DTD
        TT2=UT12
    
    J=0
    return(TT1,TT2,J)


def pymUt1utc(UT11,UT12,DUT1):
    '''
    Time scale transformation:  Universal Time, UT1, to Coordinated
    Universal Time, UTC.

    Parameters
    ----------
    ut11 : float
        UT1 as a 2-part Julian Date
    ut12 : float
        UT1 as a 2-part Julian Date
    dut1 : float
        Delta UT1: UT1-UTC in seconds

    Returns
    -------
    utc1 : float
        UTC as a 2-part quasi Julian Date
    utc2 : float
        UTC as a 2-part quasi Julian Date
    J : ValueError
        1: 'dubious year',
       -1: 'unacceptable date',    
    '''
    #一天的秒长
    D2S=86400.0
    
    #UT1-UTC的时差，单位：秒
    DUTS=DUT1
    #大的在前.
    BIG1=np.abs(UT11)>=np.abs(UT12)
    if (BIG1):
        U1=UT11
        U2=UT12
    else:
        U1=UT12
        U2=UT11
        
    #判断是否存在闰秒   
    D1=U1
    DATS1=0.0
    
    for j in range(-1,4,1):
        D2=U2+float(j)
        IY,IM,ID,FD,JS=pymJd2cal(D1, D2)
        if (JS!=0):
            break
        DATS2,JS=pymDat(IY,IM,ID,0.0)
        if (JS<0):
            break
        if (j==-1):
            DATS1=DATS2
        DDATS=DATS2-DATS1
        if (np.abs(DDATS)>=0.5):
            
            #是的，闰秒附近:确保UT1-UTC是“前”值
            if ((DDATS*DUTS)>=0.0):
                DUTS=DUTS-DDATS
            
            #将UT1作为一个UTC日的起始，并将闰秒作为结尾
            D1,D2,JS=pymCal2jd(IY,IM,ID)
            US1=D1
            US2=D2-1.0+DUTS/D2S
            
            #UT1是在这一个时间点之后么？
            DU=U1-US1
            DU=DU+(U2-US2)
            
            if (DU>0.0):
                
                #是，则这个时间点是过去的
                FD=DU*D2S/(D2S+DDATS)
                
                #在UT1-UTC中引入SOFA's JD(UTC)公约.
                DUTS=DUTS+DDATS*min(FD,1.0)
            
            break
        DATS1=DATS2
        
    #从UT1中减去UT1-UTC(可能已调整)得到UTC.
    U2 = U2 - DUTS/D2S
    
    #尽可能保留结果的精度
    if (BIG1):
        UTC1=U1
        UTC2=U2
    else:
        UTC1=U2
        UTC2=U1
 
    J=JS 
    return(UTC1,UTC2,J)


def pymGc2gde(A,F,XYZ):
    '''
    Transform geocentric coordinates to geodetic for a reference
    ellipsoid of specified form.

    Parameters
    ----------
    a : float
        equatorial radius    
    f : float
        flattening    
    xyz : list(3)
        geocentric vector

    Returns
    -------
    elong : float
        longitude (radians, east +ve)    
    phi : float
        latitude (geodetic, radians)    
    height : float
        height above ellipsoid (geodetic)
    J : ValueError
        0 = OK
       -1 = illegal F
       -2 = illegal A
    '''

    ELONG,PHI,HEIGHT=0,0,0
    
    i=1
    while i<2:
    
        #检验椭球参数
        if (F<0.0)|(F>=1.0):
            J=-1
            print('ERROR',J)
            break
        elif (A<=0.0):
            J=-2
            print('ERROR',J)
            break
        
        #椭圆参数的功能(进一步验证F)
        AEPS2=A*A*1e-32
        E2=(2.0-F)*F
        E4T=E2*E2*1.5
        EC2=1.0-E2
        if (EC2<=0.0):
            J=-1
            print('ERROR',J)
            break
        EC=ma.sqrt(EC2)
        B=A*EC
        
        #笛卡尔坐标系分量
        X=XYZ[0]
        Y=XYZ[1]
        Z=XYZ[2]
        
        #离极轴的距离.
        P2=X*X+Y*Y
        
        #经度.
        if (P2>0.0):
            ELONG=ma.atan2(Y,X)
        else:
            ELONG=0.0
        
        #无符号z坐标.
        ABSZ=np.abs(Z)
        
        #判断是否在极点.
        if (P2>AEPS2):
            
            #离极轴的距离
            P=ma.sqrt(P2)
            
            #归一化.
            S0=ABSZ/A
            PN=P/A
            ZC=EC*S0
        
            #准备牛顿修正因子
            C0=EC*PN
            C02=C0*C0
            C03=C02*C0
            S02=S0*S0
            S03=S02*S0
            A02=C02+S02
            A0=ma.sqrt(A02)
            A03=A02*A0
            D0=ZC*A03+E2*S03
            F0=PN*A03-E2*C03
            
            #准备哈雷修正系数
            B0=E4T*S02*C02*PN*(A0-EC)
            S1=D0*F0-B0*S0
            CC=EC*(F0*F0-B0*C0)
            
            #评估纬度和高度.
            PHI=ma.atan(S1/CC)
            S12=S1*S1
            CC2=CC*CC
            HEIGHT=(P*CC+ABSZ*S1-A*ma.sqrt(EC2*S12+CC2))/ma.sqrt(S12+CC2)
            
        else:
            
            #例外情况：在极点.
            PHI=ma.pi/2.0
            HEIGHT=ABSZ-B
            
        #重新赋予z的正负
        if (Z<0.0):
            PHI=-PHI
        
        J=0
        
        i+=1
            
    return(ELONG,PHI,HEIGHT,J)


def pymGC2GD(N,XYZ):
    '''
    Transform geocentric coordinates to geodetic using the specified
    reference ellipsoid.

    Parameters
    ----------
    n : int
        ellipsoid identifier    
    xyz : list(3)
        geocentric vector
        
    n    ellipsoid    
    1     WGS84    
    2     GRS80    
    3     WGS72    
    
    Returns
    -------
    elong : float
        longitude (radians, east +ve)    
    phi : float
        latitude (geodetic, radians)    
    height : float
        height above ellipsoid (geodetic)
    J : ValueError
        0 = OK
       -1 = illegal identifier 
       -2 = internal error 
    '''
    ELONG,PHI,HEIGHT=0,0,0
    
    #获得参考椭球的参数
    A,F,J=pymEform(N)
    
    #无误则调用相关函数获得所需要参数
    if (J==0):
        ELONG,PHI,HEIGHT,J=pymGc2gde(A,F,XYZ)
        if (J<0):
            J=-2
    
    #出现错误时的处理
    if (J<0):
        ELONG=-1e9
        PHI=-1e9
        HEIGHT=-1e9        
    
    return(ELONG,PHI,HEIGHT,J)


def pymPdp(A,B):
    '''
    p-vector inner (=scalar=dot) product.

    Parameters
    ----------
    a : list(3)
        first p-vector    
    b : list(3)
        second p-vector

    Returns
    -------
    function value : float
        a . b
    '''
    W=0.0
    for i in range(3):
        W=W+A[i]*B[i]
    
    ADB=W
    return(ADB)


def pymAb(PNAT,V,S,BM1):
    '''
    Apply aberration to transform natural direction into proper direction.

    Parameters
    ----------
    pnat : list(3)
        natural direction to the source (unit vector)
    v : list(3)
        observer barycentric velocity in units of c
    s : float
        distance between the Sun and the observer (au)
    bm1 : float
        sqrt(1-|v|^2): reciprocal of Lorenz factor

    Returns
    -------
    ppr : list(3)
        proper direction to source (unit vector)

    '''
    #太阳的施瓦西半径，单位：au
    #= 2 * 1.32712440041 e20 / (2.99792458 e8)^2 / 1.49597870700 e11
    SRS=1.97412574336e-08
    
    PPR=[0,0,0]
    P=[0,0,0]
    PDV=pymPdp(PNAT,V)
    W1=1.0+PDV/(1.0+BM1)
    W2=SRS/S
    R2=0.0
    
    for i in range(3):
        W=PNAT[i]*BM1+W1*V[i]+W2*(V[i]-PDV*PNAT[i])
        P[i]=W
        R2=R2+W*W
    
    R=ma.sqrt(R2)
    for i in range(3):
        PPR[i]=P[i]/R
    
    return(PPR)


def pymAe2hd(AZ,EL,PHI):
    '''
    Horizon to equatorial coordinates:  transform azimuth and altitude
    to hour angle and declination.

    Parameters
    ----------
    az : float
        azimuth
    el : float
        altitude (informally, elevation)
    phi : float
        site latitude

    Returns
    -------
    ha : float
        hour angle (local)
    dec : float
        declination

    '''   
    SA=ma.sin(AZ)
    CA=ma.cos(AZ)
    SE=ma.sin(EL)
    CE=ma.cos(EL)
    SP=ma.sin(PHI)
    CP=ma.cos(PHI)
    
    #单位向量
    X=-CA*CE*SP+SE*CP
    Y=-SA*CE
    Z=CA*CE*CP+SE*SP
    
    #转换到球面
    R=ma.sqrt(X*X+Y*Y)
    if (R==0):
        HA=0.0
    else:
        HA=ma.atan2(Y,X)
    DEC=ma.atan2(Z,R)
    
    return(HA,DEC)


def pymAnp(A):
    '''
    Normalize angle into the range 0 <= a < 2pi.

    Parameters
    ----------
    a : float
        angle (radians)

    Returns
    -------
    function value : float
        angle in range 0-2pi

    '''
    W=np.abs(A)%(2*ma.pi)
    if A<0:
        W=-W
    if (W<0.0):
        W=W+2*ma.pi
    
    return(W)
    

def pymAnpm(A):
    '''
    Normalize angle into the range -pi <= a < +pi.

    Parameters
    ----------
    a : float
        angle (radians)

    Returns
    -------
    function value : float
        angle in range +/-pi

    '''
    W=np.abs(A)%(2*ma.pi)
    if A<0:
        W=-W
    if (np.abs(W)>=ma.pi):
        if A>=0:
            W=W-2*ma.pi
        else:
            W=W+2*ma.pi
        
    return(W)


def pymSxp(S,P):
    '''
    Multiply a p-vector by a scalar.

    Parameters
    ----------
    s : float
        scalar    
    p : list(3)
        p-vector

    Returns
    -------
    sp : list(3)
        s * p

    '''
    SP=[0,0,0]
    
    for i in range(3):
        SP[i]=S*P[i]
    
    return(SP)


def pymPm(P):
    '''
    Modulus of p-vector.

    Parameters
    ----------
    p : list(3)
        p-vector

    Returns
    -------
    function value : float
        modulus

    '''
    W=0.0
    for i in range(3):
        C=P[i]
        W+=C**2
    R=ma.sqrt(W)
        
    return(R)


def pymPn(P):
    '''
    Convert a p-vector into modulus and unit vector.

    Parameters
    ----------
    p : list(3)
        p-vector

    Returns
    -------
    r : float
        modulus    
    u : list(3)
        unit vector

    '''
    U=[0,0,0]
    
    #调用pymPm函数获得向量的模，并判断是否为零向量
    W=pymPm(P)
    if (W==0.0):
        
        #零向量
        U=pymZp()
    else:
        
        #单位向量
        U=pymSxp(1.0/W, P)
    
    R=W
    return(R,U)


def pymIr():
    '''
    Initialize an r-matrix to the identity matrix.

    Returns
    -------
    r : list(3,3)
        r-matrix

    '''

    R=[[1,0,0],[0,1,0],[0,0,1]]
    
    return(R)


def pymCp(P):
    '''
    Copy a p-vector.

    Parameters
    ----------
    p : list(3)
        p-vector to be copied

    Returns
    -------
    c : list(3)
        copy

    '''
    C=[P[i] for i in range(3)]
  
    return(C)


def pymApcs(DATE1,DATE2,PV,EBPV,EHP,ASTROM):
    '''
    For an observer whose geocentric position and velocity are known,
    prepare star-independent astrometry parameters for transformations
    between ICRS and GCRS.  The Earth ephemeris is supplied by the
    caller.
    
    The parameters produced by this function are required in the space
    motion, parallax, light deflection and aberration parts of the
    astrometric transformation chain.

    Parameters
    ----------
    date1 : float
        TDB as a 2-part Julian Date    
    date2 : float
        TDB as a 2-part Julian Date    
    pv : list(2,3)
        observer's geocentric pos/vel (m, m/s)    
    ebpv : list(2,3)
        Earth barycentric position/velocity (au, au/day)     
    ehp : list(1,3)
        Earth heliocentric position (au)    

    Returns
    -------
    astrom : list(30)
        star-independent astrometry parameters     
        [0]>pmt : PM time interval (SSB, Julian years)
        [1-3]>eb : SSB to observer (vector, au)
        [4-6]>eh : Sun to observer (unit vector)
        [7]>em : distance from Sun to observer (au)
        [8-10]>v : barycentric observer velocity (vector, c)
        [11]>bm1 : sqrt(1-|v|^2): reciprocal of Lorenz factor
        [12-20]>bpn : bias-precession-nutation matrix
        [21]>along : unchanged
        [22]>xpl : unchanged
        [23]>ypl : unchanged
        [24]>sphi : unchanged
        [25]>cphi : unchanged
        [26]>diurab : unchanged
        [27]>eral : unchanged
        [28]>refa : unchanged 
        [29]>refb : unchanged

    '''
    PB,VB,PH=[0,0,0],[0,0,0],[0,0,0]
    #参考历元
    DJ00=2451545.0
    
    #一个儒略年的天数
    DJY=365.25
    
    #每天的秒数
    DAYSEC=86400.0
    
    #光速(m/s)
    CMPS=299792458.0
    
    #天文单位AU(m, IAU 2012)
    AUM=149597870.7e3
    
    #光走过1AU所需要的时间(s)
    AULT=AUM/CMPS
    
    #AU/天到m/s的单位转换
    AUDMS=AUM/DAYSEC
    
    #光走过1AU所需要的时间(day)
    CR=AULT/DAYSEC
    
    #相对于参考历元经过的儒略年数
    ASTROM[0]=((DATE1-DJ00)+DATE2)/DJY
    
    #根据地球历元调整观测者的位置、速度参数.
    for i in range(3):
        DP=PV[0][i]/AUM
        DV=PV[1][i]/AUDMS
        PB[i]=EBPV[0][i]+DP
        VB[i]=EBPV[1][i]+DV
        PH[i]=EHP[i]+DP
    
    #观察者相对于太阳系质心的位置(au).
    A=pymCp(PB)
    ASTROM[1]=A[0]
    ASTROM[2]=A[1]
    ASTROM[3]=A[2]
    
    #观察者相对于太阳的方向和距离(单位向量|单位： au).
    ASTROM[7],B=pymPn(PH)
    ASTROM[4]=B[0]
    ASTROM[5]=B[1]
    ASTROM[6]=B[2]
    
    #观察者的质心速度，单位：c，以及洛伦兹因子的倒数
    V2=0.0
    for i in range(3):
        W=VB[i]*CR
        ASTROM[8+i]=W
        V2=V2+W*W
    ASTROM[11]=ma.sqrt(1.0-V2)
    
    #将NPB矩阵设置为单位矩阵
    for i in range(3):
        ASTROM[12+i*4]=1
    
    return(ASTROM)


def pymZpv():
    '''
    Zero a pv-vector

    Parameters
    ----------
    pv : list(2,3)
        pv-vector
        
    Returns
    -------
    pv : list(2,3)
        zero pv-vector

    '''

    PV=[[0,0,0],[0,0,0]]
    
    return(PV)


def pymApcg(DATE1,DATE2,EBPV,EHP,ASTROM):
    '''
    For a geocentric observer, prepare star-independent astrometry
    parameters for transformations between ICRS and GCRS coordinates.
    The Earth ephemeris is supplied by the caller.
    
    The parameters produced by this function are required in the
    parallax, light deflection and aberration parts of the astrometric
    transformation chain.

    Parameters
    ----------
    date1 : float
        TDB as a 2-part Julian Date    
    date2 : float
        TDB as a 2-part Julian Date    
    ebpv : list(2,3)
        Earth barycentric pos/vel (au, au/day)     
    ehp : list(3)
        Earth heliocentric position (au)

    Returns
    -------
    astrom : list(30)
        star-independent astrometry parameters    
        [0]>pmt : PM time interval (SSB, Julian years)
        [1-3]>eb : SSB to observer (vector, au)
        [4-6]>eh : Sun to observer (unit vector)
        [7]>em : distance from Sun to observer (au)
        [8-10]>v : barycentric observer velocity (vector, c)
        [11]>bm1 : sqrt(1-|v|^2): reciprocal of Lorenz factor
        [12-20]>bpn : bias-precession-nutation matrix
        [21]>along : unchanged
        [22]>xpl : unchanged
        [23]>ypl : unchanged
        [24]>sphi : unchanged
        [25]>cphi : unchanged
        [26]>diurab : unchanged
        [27]>eral : unchanged
        [28]>refa : unchanged 
        [29]>refb : unchanged
    '''
    #观测者位于地心，其相对于地心的位置、速度矢量都为0
    PV=pymZpv()
    
    #调用pymApcs函数，将为0的PV向量带入
    ASTROM=pymApcs(DATE1,DATE2,PV,EBPV,EHP,ASTROM)
    
    return(ASTROM)


def pymEpv00(DATE1,DATE2):
    '''
    Earth position and velocity, heliocentric and barycentric, with
    respect to the Barycentric Celestial Reference System.

    Parameters
    ----------
    date1 : float
        TDB date    
    date2 : float
        TDB date

    Returns
    -------
    pvh : list(2,3)
        heliocentric Earth position/velocity    
    pvb : list(2,3)
        barycentric Earth position/velocity
    J : ValueError
        0 : OK
       +1 : date out of range(1900-2100)
    '''
    
    PVB=[[0 for i in range(3)] for j in range(2)]
    PVH=[[0 for i in range(3)] for j in range(2)]    
    JSTAT=0
    
    PH,VH,PB,VB=[0,0,0],[0,0,0],[0,0,0],[0,0,0]
    
    #1儒略日年的天数
    DJY=365.25
    
    #参考历元(J2000.0), JD
    DJ00=2451545.0
    
    #用于将分析模型定向到DE405/MICRF的矩阵参数
    #其对应的欧拉角（经验公式）
    #                        d  '  "
    #    1st rotation    -  23 26 21.4091 about the x-axis  (obliquity)
    #    2nd rotation    +         0.0475 about the z-axis  (RA offset)
    #
    #  These were obtained empirically, by comparisons with DE405 over
    #  1900-2100.
    
    AM12=+0.000000211284
    AM13=-0.000000091603
    AM21=-0.000000230286
    AM22=+0.917482137087
    AM23=-0.397776982902
    AM32=+0.397776982902
    AM33=+0.917482137087
    
    #星表系数
    #Xi[j]
    #X;E：地球相对于太阳，S：太阳相对于太阳系质心
    #i;0:T^0项，1：T^1项，2：T^2项
    #j;0:X组份，1：Y组份，2:Z组份
    #其中3列数据分别表示：振幅、相位、频率
    
    NE0=[0,0,0]
    NE1=[0,0,0]
    NE2=[0,0,0]
    NS0=[0,0,0]
    NS1=[0,0,0]
    NS2=[0,0,0]
    NE0X = 501
    NE0Y = 501
    NE0Z = 137
    ME0 = NE0X
    NE1X =  79
    NE1Y =  80
    NE1Z =  12
    ME1 = NE1Y
    NE2X =   5
    NE2Y =   5
    NE2Z =   3
    ME2 = NE2X
    NS0X = 212
    NS0Y = 213
    NS0Z =  69
    MS0 = NS0Y
    NS1X =  50
    NS1Y =  50
    NS1Z =  14
    MS1 = NS1X
    NS2X =   9
    NS2Y =   9
    NS2Z =   2
    MS2 = NS2X 
    
    E0=[0,0,0]
    E1=[0,0,0]
    E2=[0,0,0]
    S0=[0,0,0]
    S1=[0,0,0]
    S2=[0,0,0]
    
    NE0=[NE0X, NE0Y, NE0Z]
    NE1=[NE1X, NE1Y, NE1Z]
    NE2=[NE2X, NE2Y, NE2Z]
    NS0=[NS0X, NS0Y, NS0Z]
    NS1=[NS1X, NS1Y, NS1Z]
    NS2=[NS2X, NS2Y, NS2Z]
    
    #Sun-to-Earth, T^0, X
    E0[0]=[[0.9998292878132e+00,0.1753485171504e+01,0.6283075850446e+01],
             [0.8352579567414e-02,0.1710344404582e+01,0.1256615170089e+02],
             [0.5611445335148e-02,0.0000000000000e+00,0.0000000000000e+00],
             [0.1046664295572e-03,0.1667225416770e+01,0.1884922755134e+02],
             [0.3110842534677e-04,0.6687513390251e+00,0.8399684731857e+02],
             [0.2552413503550e-04,0.5830637358413e+00,0.5296909721118e+00],
             [0.2137207845781e-04,0.1092330954011e+01,0.1577343543434e+01],
             [0.1680240182951e-04,0.4955366134987e+00,0.6279552690824e+01],
             [0.1679012370795e-04,0.6153014091901e+01,0.6286599010068e+01],
             [0.1445526946777e-04,0.3472744100492e+01,0.2352866153506e+01],
             [0.1091038246184e-04,0.3689845786119e+01,0.5223693906222e+01],
             [0.9344399733932e-05,0.6073934645672e+01,0.1203646072878e+02],
             [0.8993182910652e-05,0.3175705249069e+01,0.1021328554739e+02],
             [0.5665546034116e-05,0.2152484672246e+01,0.1059381944224e+01],
             [0.6844146703035e-05,0.1306964099750e+01,0.5753384878334e+01],
             [0.7346610905565e-05,0.4354980070466e+01,0.3981490189893e+00],
             [0.6815396474414e-05,0.2218229211267e+01,0.4705732307012e+01],
             [0.6112787253053e-05,0.5384788425458e+01,0.6812766822558e+01],
             [0.4518120711239e-05,0.6087604012291e+01,0.5884926831456e+01],
             [0.4521963430706e-05,0.1279424524906e+01,0.6256777527156e+01],
             [0.4497426764085e-05,0.5369129144266e+01,0.6309374173736e+01],
             [0.4062190566959e-05,0.5436473303367e+00,0.6681224869435e+01],
             [0.5412193480192e-05,0.7867838528395e+00,0.7755226100720e+00],
             [0.5469839049386e-05,0.1461440311134e+01,0.1414349524433e+02],
             [0.5205264083477e-05,0.4432944696116e+01,0.7860419393880e+01],
             [0.2149759935455e-05,0.4502237496846e+01,0.1150676975667e+02],
             [0.2279109618501e-05,0.1239441308815e+01,0.7058598460518e+01],
             [0.2259282939683e-05,0.3272430985331e+01,0.4694002934110e+01],
             [0.2558950271319e-05,0.2265471086404e+01,0.1216800268190e+02],
             [0.2561581447555e-05,0.1454740653245e+01,0.7099330490126e+00],
             [0.1781441115440e-05,0.2962068630206e+01,0.7962980379786e+00],
             [0.1612005874644e-05,0.1473255041006e+01,0.5486777812467e+01],
             [0.1818630667105e-05,0.3743903293447e+00,0.6283008715021e+01],
             [0.1818601377529e-05,0.6274174354554e+01,0.6283142985870e+01],
             [0.1554475925257e-05,0.1624110906816e+01,0.2513230340178e+02],
             [0.2090948029241e-05,0.5852052276256e+01,0.1179062909082e+02],
             [0.2000176345460e-05,0.4072093298513e+01,0.1778984560711e+02],
             [0.1289535917759e-05,0.5217019331069e+01,0.7079373888424e+01],
             [0.1281135307881e-05,0.4802054538934e+01,0.3738761453707e+01],
             [0.1518229005692e-05,0.8691914742502e+00,0.2132990797783e+00],
             [0.9450128579027e-06,0.4601859529950e+01,0.1097707878456e+02],
             [0.7781119494996e-06,0.1844352816694e+01,0.8827390247185e+01],
             [0.7733407759912e-06,0.3582790154750e+01,0.5507553240374e+01],
             [0.7350644318120e-06,0.2695277788230e+01,0.1589072916335e+01],
             [0.6535928827023e-06,0.3651327986142e+01,0.1176985366291e+02],
             [0.6324624183656e-06,0.2241302375862e+01,0.6262300422539e+01],
             [0.6298565300557e-06,0.4407122406081e+01,0.6303851278352e+01],
             [0.8587037089179e-06,0.3024307223119e+01,0.1672837615881e+03],
             [0.8299954491035e-06,0.6192539428237e+01,0.3340612434717e+01],
             [0.6311263503401e-06,0.2014758795416e+01,0.7113454667900e-02],
             [0.6005646745452e-06,0.3399500503397e+01,0.4136910472696e+01],
             [0.7917715109929e-06,0.2493386877837e+01,0.6069776770667e+01],
             [0.7556958099685e-06,0.4159491740143e+01,0.6496374930224e+01],
             [0.6773228244949e-06,0.4034162934230e+01,0.9437762937313e+01],
             [0.5370708577847e-06,0.1562219163734e+01,0.1194447056968e+01],
             [0.5710804266203e-06,0.2662730803386e+01,0.6282095334605e+01],
             [0.5709824583726e-06,0.3985828430833e+01,0.6284056366286e+01],
             [0.5143950896447e-06,0.1308144688689e+01,0.6290189305114e+01],
             [0.5088010604546e-06,0.5352817214804e+01,0.6275962395778e+01],
             [0.4960369085172e-06,0.2644267922349e+01,0.6127655567643e+01],
             [0.4803137891183e-06,0.4008844192080e+01,0.6438496133249e+01],
             [0.5731747768225e-06,0.3794550174597e+01,0.3154687086868e+01],
             [0.4735947960579e-06,0.6107118308982e+01,0.3128388763578e+01],
             [0.4808348796625e-06,0.4771458618163e+01,0.8018209333619e+00],
             [0.4115073743137e-06,0.3327111335159e+01,0.8429241228195e+01],
             [0.5230575889287e-06,0.5305708551694e+01,0.1336797263425e+02],
             [0.5133977889215e-06,0.5784230738814e+01,0.1235285262111e+02],
             [0.5065815825327e-06,0.2052064793679e+01,0.1185621865188e+02],
             [0.4339831593868e-06,0.3644994195830e+01,0.1726015463500e+02],
             [0.3952928638953e-06,0.4930376436758e+01,0.5481254917084e+01],
             [0.4898498111942e-06,0.4542084219731e+00,0.9225539266174e+01],
             [0.4757490209328e-06,0.3161126388878e+01,0.5856477690889e+01],
             [0.4727701669749e-06,0.6214993845446e+00,0.2544314396739e+01],
             [0.3800966681863e-06,0.3040132339297e+01,0.4265981595566e+00],
             [0.3257301077939e-06,0.8064977360087e+00,0.3930209696940e+01],
             [0.3255810528674e-06,0.1974147981034e+01,0.2146165377750e+01],
             [0.3252029748187e-06,0.2845924913135e+01,0.4164311961999e+01],
             [0.3255505635308e-06,0.3017900824120e+01,0.5088628793478e+01],
             [0.2801345211990e-06,0.6109717793179e+01,0.1256967486051e+02],
             [0.3688987740970e-06,0.2911550235289e+01,0.1807370494127e+02],
             [0.2475153429458e-06,0.2179146025856e+01,0.2629832328990e-01],
             [0.3033457749150e-06,0.1994161050744e+01,0.4535059491685e+01],
             [0.2186743763110e-06,0.5125687237936e+01,0.1137170464392e+02],
             [0.2764777032774e-06,0.4822646860252e+00,0.1256262854127e+02],
             [0.2199028768592e-06,0.4637633293831e+01,0.1255903824622e+02],
             [0.2046482824760e-06,0.1467038733093e+01,0.7084896783808e+01],
             [0.2611209147507e-06,0.3044718783485e+00,0.7143069561767e+02],
             [0.2286079656818e-06,0.4764220356805e+01,0.8031092209206e+01],
             [0.1855071202587e-06,0.3383637774428e+01,0.1748016358760e+01],
             [0.2324669506784e-06,0.6189088449251e+01,0.1831953657923e+02],
             [0.1709528015688e-06,0.5874966729774e+00,0.4933208510675e+01],
             [0.2168156875828e-06,0.4302994009132e+01,0.1044738781244e+02],
             [0.2106675556535e-06,0.3800475419891e+01,0.7477522907414e+01],
             [0.1430213830465e-06,0.1294660846502e+01,0.2942463415728e+01],
             [0.1388396901944e-06,0.4594797202114e+01,0.8635942003952e+01],
             [0.1922258844190e-06,0.4943044543591e+00,0.1729818233119e+02],
             [0.1888460058292e-06,0.2426943912028e+01,0.1561374759853e+03],
             [0.1789449386107e-06,0.1582973303499e+00,0.1592596075957e+01],
             [0.1360803685374e-06,0.5197240440504e+01,0.1309584267300e+02],
             [0.1504038014709e-06,0.3120360916217e+01,0.1649636139783e+02],
             [0.1382769533389e-06,0.6164702888205e+01,0.7632943190217e+01],
             [0.1438059769079e-06,0.1437423770979e+01,0.2042657109477e+02],
             [0.1326303260037e-06,0.3609688799679e+01,0.1213955354133e+02],
             [0.1159244950540e-06,0.5463018167225e+01,0.5331357529664e+01],
             [0.1433118149136e-06,0.6028909912097e+01,0.7342457794669e+01],
             [0.1234623148594e-06,0.3109645574997e+01,0.6279485555400e+01],
             [0.1233949875344e-06,0.3539359332866e+01,0.6286666145492e+01],
             [0.9927196061299e-07,0.1259321569772e+01,0.7234794171227e+01],
             [0.1242302191316e-06,0.1065949392609e+01,0.1511046609763e+02],
             [0.1098402195201e-06,0.2192508743837e+01,0.1098880815746e+02],
             [0.1158191395315e-06,0.4054411278650e+01,0.5729506548653e+01],
             [0.9048475596241e-07,0.5429764748518e+01,0.9623688285163e+01],
             [0.8889853269023e-07,0.5046586206575e+01,0.6148010737701e+01],
             [0.1048694242164e-06,0.2628858030806e+01,0.6836645152238e+01],
             [0.1112308378646e-06,0.4177292719907e+01,0.1572083878776e+02],
             [0.8631729709901e-07,0.1601345232557e+01,0.6418140963190e+01],
             [0.8527816951664e-07,0.2463888997513e+01,0.1471231707864e+02],
             [0.7892139456991e-07,0.3154022088718e+01,0.2118763888447e+01],
             [0.1051782905236e-06,0.4795035816088e+01,0.1349867339771e+01],
             [0.1048219943164e-06,0.2952983395230e+01,0.5999216516294e+01],
             [0.7435760775143e-07,0.5420547991464e+01,0.6040347114260e+01],
             [0.9869574106949e-07,0.3695646753667e+01,0.6566935184597e+01],
             [0.9156886364226e-07,0.3922675306609e+01,0.5643178611111e+01],
             [0.7006834356188e-07,0.1233968624861e+01,0.6525804586632e+01],
             [0.9806170182601e-07,0.1919542280684e+01,0.2122839202813e+02],
             [0.9052289673607e-07,0.4615902724369e+01,0.4690479774488e+01],
             [0.7554200867893e-07,0.1236863719072e+01,0.1253985337760e+02],
             [0.8215741286498e-07,0.3286800101559e+00,0.1097355562493e+02],
             [0.7185178575397e-07,0.5880942158367e+01,0.6245048154254e+01],
             [0.7130726476180e-07,0.7674871987661e+00,0.6321103546637e+01],
             [0.6650894461162e-07,0.6987129150116e+00,0.5327476111629e+01],
             [0.7396888823688e-07,0.3576824794443e+01,0.5368044267797e+00],
             [0.7420588884775e-07,0.5033615245369e+01,0.2354323048545e+02],
             [0.6141181642908e-07,0.9449927045673e+00,0.1296430071988e+02],
             [0.6373557924058e-07,0.6206342280341e+01,0.9517183207817e+00],
             [0.6359474329261e-07,0.5036079095757e+01,0.1990745094947e+01],
             [0.5740173582646e-07,0.6105106371350e+01,0.9555997388169e+00],
             [0.7019864084602e-07,0.7237747359018e+00,0.5225775174439e+00],
             [0.6398054487042e-07,0.3976367969666e+01,0.2407292145756e+02],
             [0.7797092650498e-07,0.4305423910623e+01,0.2200391463820e+02],
             [0.6466760000900e-07,0.3500136825200e+01,0.5230807360890e+01],
             [0.7529417043890e-07,0.3514779246100e+01,0.1842262939178e+02],
             [0.6924571140892e-07,0.2743457928679e+01,0.1554202828031e+00],
             [0.6220798650222e-07,0.2242598118209e+01,0.1845107853235e+02],
             [0.5870209391853e-07,0.2332832707527e+01,0.6398972393349e+00],
             [0.6263953473888e-07,0.2191105358956e+01,0.6277552955062e+01],
             [0.6257781390012e-07,0.4457559396698e+01,0.6288598745829e+01],
             [0.5697304945123e-07,0.3499234761404e+01,0.1551045220144e+01],
             [0.6335438746791e-07,0.6441691079251e+00,0.5216580451554e+01],
             [0.6377258441152e-07,0.2252599151092e+01,0.5650292065779e+01],
             [0.6484841818165e-07,0.1992812417646e+01,0.1030928125552e+00],
             [0.4735551485250e-07,0.3744672082942e+01,0.1431416805965e+02],
             [0.4628595996170e-07,0.1334226211745e+01,0.5535693017924e+00],
             [0.6258152336933e-07,0.4395836159154e+01,0.2608790314060e+02],
             [0.6196171366594e-07,0.2587043007997e+01,0.8467247584405e+02],
             [0.6159556952126e-07,0.4782499769128e+01,0.2394243902548e+03],
             [0.4987741172394e-07,0.7312257619924e+00,0.7771377146812e+02],
             [0.5459280703142e-07,0.3001376372532e+01,0.6179983037890e+01],
             [0.4863461189999e-07,0.3767222128541e+01,0.9027992316901e+02],
             [0.5349912093158e-07,0.3663594450273e+01,0.6386168663001e+01],
             [0.5673725607806e-07,0.4331187919049e+01,0.6915859635113e+01],
             [0.4745485060512e-07,0.5816195745518e+01,0.6282970628506e+01],
             [0.4745379005326e-07,0.8323672435672e+00,0.6283181072386e+01],
             [0.4049002796321e-07,0.3785023976293e+01,0.6254626709878e+01],
             [0.4247084014515e-07,0.2378220728783e+01,0.7875671926403e+01],
             [0.4026912363055e-07,0.2864103423269e+01,0.6311524991013e+01],
             [0.4062935011774e-07,0.2415408595975e+01,0.3634620989887e+01],
             [0.5347771048509e-07,0.3343479309801e+01,0.2515860172507e+02],
             [0.4829494136505e-07,0.2821742398262e+01,0.5760498333002e+01],
             [0.4342554404599e-07,0.5624662458712e+01,0.7238675589263e+01],
             [0.4021599184361e-07,0.5557250275009e+00,0.1101510648075e+02],
             [0.4104900474558e-07,0.3296691780005e+01,0.6709674010002e+01],
             [0.4376532905131e-07,0.3814443999443e+01,0.6805653367890e+01],
             [0.3314590480650e-07,0.3560229189250e+01,0.1259245002418e+02],
             [0.3232421839643e-07,0.5185389180568e+01,0.1066495398892e+01],
             [0.3541176318876e-07,0.3921381909679e+01,0.9917696840332e+01],
             [0.3689831242681e-07,0.4190658955386e+01,0.1192625446156e+02],
             [0.3890605376774e-07,0.5546023371097e+01,0.7478166569050e-01],
             [0.3038559339780e-07,0.6231032794494e+01,0.1256621883632e+02],
             [0.3137083969782e-07,0.6207063419190e+01,0.4292330755499e+01],
             [0.4024004081854e-07,0.1195257375713e+01,0.1334167431096e+02],
             [0.3300234879283e-07,0.1804694240998e+01,0.1057540660594e+02],
             [0.3635399155575e-07,0.5597811343500e+01,0.6208294184755e+01],
             [0.3032668691356e-07,0.3191059366530e+01,0.1805292951336e+02],
             [0.2809652069058e-07,0.4094348032570e+01,0.3523159621801e-02],
             [0.3696955383823e-07,0.5219282738794e+01,0.5966683958112e+01],
             [0.3562894142503e-07,0.1037247544554e+01,0.6357857516136e+01],
             [0.3510598524148e-07,0.1430020816116e+01,0.6599467742779e+01],
             [0.3617736142953e-07,0.3002911403677e+01,0.6019991944201e+01],
             [0.2624524910730e-07,0.2437046757292e+01,0.6702560555334e+01],
             [0.2535824204490e-07,0.1581594689647e+01,0.3141537925223e+02],
             [0.3519787226257e-07,0.5379863121521e+01,0.2505706758577e+03],
             [0.2578406709982e-07,0.4904222639329e+01,0.1673046366289e+02],
             [0.3423887981473e-07,0.3646448997315e+01,0.6546159756691e+01],
             [0.2776083886467e-07,0.3307829300144e+01,0.1272157198369e+02],
             [0.3379592818379e-07,0.1747541251125e+01,0.1494531617769e+02],
             [0.3050255426284e-07,0.1784689432607e-01,0.4732030630302e+01],
             [0.2652378350236e-07,0.4420055276260e+01,0.5863591145557e+01],
             [0.2374498173768e-07,0.3629773929208e+01,0.2388894113936e+01],
             [0.2716451255140e-07,0.3079623706780e+01,0.1202934727411e+02],
             [0.3038583699229e-07,0.3312487903507e+00,0.1256608456547e+02],
             [0.2220681228760e-07,0.5265520401774e+01,0.1336244973887e+02],
             [0.3044156540912e-07,0.4766664081250e+01,0.2908881142201e+02],
             [0.2731859923561e-07,0.5069146530691e+01,0.1391601904066e+02],
             [0.2285603018171e-07,0.5954935112271e+01,0.6076890225335e+01],
             [0.2025006454555e-07,0.4061789589267e+01,0.4701116388778e+01],
             [0.2012597519804e-07,0.2485047705241e+01,0.6262720680387e+01],
             [0.2003406962258e-07,0.4163779209320e+01,0.6303431020504e+01],
             [0.2207863441371e-07,0.6923839133828e+00,0.6489261475556e+01],
             [0.2481374305624e-07,0.5944173595676e+01,0.1204357418345e+02],
             [0.2130923288870e-07,0.4641013671967e+01,0.5746271423666e+01],
             [0.2446370543391e-07,0.6125796518757e+01,0.1495633313810e+00],
             [0.1932492759052e-07,0.2234572324504e+00,0.1352175143971e+02],
             [0.2600122568049e-07,0.4281012405440e+01,0.4590910121555e+01],
             [0.2431754047488e-07,0.1429943874870e+00,0.1162474756779e+01],
             [0.1875902869209e-07,0.9781803816948e+00,0.6279194432410e+01],
             [0.1874381139426e-07,0.5670368130173e+01,0.6286957268481e+01],
             [0.2156696047173e-07,0.2008985006833e+01,0.1813929450232e+02],
             [0.1965076182484e-07,0.2566186202453e+00,0.4686889479442e+01],
             [0.2334816372359e-07,0.4408121891493e+01,0.1002183730415e+02],
             [0.1869937408802e-07,0.5272745038656e+01,0.2427287361862e+00],
             [0.2436236460883e-07,0.4407720479029e+01,0.9514313292143e+02],
             [0.1761365216611e-07,0.1943892315074e+00,0.1351787002167e+02],
             [0.2156289480503e-07,0.1418570924545e+01,0.6037244212485e+01],
             [0.2164748979255e-07,0.4724603439430e+01,0.2301353951334e+02],
             [0.2222286670853e-07,0.2400266874598e+01,0.1266924451345e+02],
             [0.2070901414929e-07,0.5230348028732e+01,0.6528907488406e+01],
             [0.1792745177020e-07,0.2099190328945e+01,0.6819880277225e+01],
             [0.1841802068445e-07,0.3467527844848e+00,0.6514761976723e+02],
             [0.1578401631718e-07,0.7098642356340e+00,0.2077542790660e-01],
             [0.1561690152531e-07,0.5943349620372e+01,0.6272439236156e+01],
             [0.1558591045463e-07,0.7040653478980e+00,0.6293712464735e+01],
             [0.1737356469576e-07,0.4487064760345e+01,0.1765478049437e+02],
             [0.1434755619991e-07,0.2993391570995e+01,0.1102062672231e+00],
             [0.1482187806654e-07,0.2278049198251e+01,0.1052268489556e+01],
             [0.1424812827089e-07,0.1682114725827e+01,0.1311972100268e+02],
             [0.1380282448623e-07,0.3262668602579e+01,0.1017725758696e+02],
             [0.1811481244566e-07,0.3187771221777e+01,0.1887552587463e+02],
             [0.1504446185696e-07,0.5650162308647e+01,0.7626583626240e-01],
             [0.1740776154137e-07,0.5487068607507e+01,0.1965104848470e+02],
             [0.1374339536251e-07,0.5745688172201e+01,0.6016468784579e+01],
             [0.1761377477704e-07,0.5748060203659e+01,0.2593412433514e+02],
             [0.1535138225795e-07,0.6226848505790e+01,0.9411464614024e+01],
             [0.1788140543676e-07,0.6189318878563e+01,0.3301902111895e+02],
             [0.1375002807996e-07,0.5371812884394e+01,0.6327837846670e+00],
             [0.1242115758632e-07,0.1471687569712e+01,0.3894181736510e+01],
             [0.1450977333938e-07,0.4143836662127e+01,0.1277945078067e+02],
             [0.1297579575023e-07,0.9003477661957e+00,0.6549682916313e+01],
             [0.1462667934821e-07,0.5760505536428e+01,0.1863592847156e+02],
             [0.1381774374799e-07,0.1085471729463e+01,0.2379164476796e+01],
             [0.1682333169307e-07,0.5409870870133e+01,0.1620077269078e+02],
             [0.1190812918837e-07,0.1397205174601e+01,0.1149965630200e+02],
             [0.1221434762106e-07,0.9001804809095e+00,0.1257326515556e+02],
             [0.1549934644860e-07,0.4262528275544e+01,0.1820933031200e+02],
             [0.1252138953050e-07,0.1411642012027e+01,0.6993008899458e+01],
             [0.1237078905387e-07,0.2844472403615e+01,0.2435678079171e+02],
             [0.1446953389615e-07,0.5295835522223e+01,0.3813291813120e-01],
             [0.1388446457170e-07,0.4969428135497e+01,0.2458316379602e+00],
             [0.1019339179228e-07,0.2491369561806e+01,0.6112403035119e+01],
             [0.1258880815343e-07,0.4679426248976e+01,0.5429879531333e+01],
             [0.1297768238261e-07,0.1074509953328e+01,0.1249137003520e+02],
             [0.9913505718094e-08,0.4735097918224e+01,0.6247047890016e+01],
             [0.9830453155969e-08,0.4158649187338e+01,0.6453748665772e+01],
             [0.1192615865309e-07,0.3438208613699e+01,0.6290122169689e+01],
             [0.9835874798277e-08,0.1913300781229e+01,0.6319103810876e+01],
             [0.9639087569277e-08,0.9487683644125e+00,0.8273820945392e+01],
             [0.1175716107001e-07,0.3228141664287e+01,0.6276029531202e+01],
             [0.1018926508678e-07,0.2216607854300e+01,0.1254537627298e+02],
             [0.9500087869225e-08,0.2625116459733e+01,0.1256517118505e+02],
             [0.9664192916575e-08,0.5860562449214e+01,0.6259197520765e+01],
             [0.9612858712203e-08,0.7885682917381e+00,0.6306954180126e+01],
             [0.1117645675413e-07,0.3932148831189e+01,0.1779695906178e+02],
             [0.1158864052160e-07,0.9995605521691e+00,0.1778273215245e+02],
             [0.9021043467028e-08,0.5263769742673e+01,0.6172869583223e+01],
             [0.8836134773563e-08,0.1496843220365e+01,0.1692165728891e+01],
             [0.1045872200691e-07,0.7009039517214e+00,0.2204125344462e+00],
             [0.1211463487798e-07,0.4041544938511e+01,0.8257698122054e+02],
             [0.8541990804094e-08,0.1447586692316e+01,0.6393282117669e+01],
             [0.1038720703636e-07,0.4594249718112e+00,0.1550861511662e+02],
             [0.1126722351445e-07,0.3925550579036e+01,0.2061856251104e+00],
             [0.8697373859631e-08,0.4411341856037e+01,0.9491756770005e+00],
             [0.8869380028441e-08,0.2402659724813e+01,0.3903911373650e+01],
             [0.9247014693258e-08,0.1401579743423e+01,0.6267823317922e+01],
             [0.9205062930950e-08,0.5245978000814e+01,0.6298328382969e+01],
             [0.8000745038049e-08,0.3590803356945e+01,0.2648454860559e+01],
             [0.9168973650819e-08,0.2470150501679e+01,0.1498544001348e+03],
             [0.1075444949238e-07,0.1328606161230e+01,0.3694923081589e+02],
             [0.7817298525817e-08,0.6162256225998e+01,0.4804209201333e+01],
             [0.9541469226356e-08,0.3942568967039e+01,0.1256713221673e+02],
             [0.9821910122027e-08,0.2360246287233e+00,0.1140367694411e+02],
             [0.9897822023777e-08,0.4619805634280e+01,0.2280573557157e+02],
             [0.7737289283765e-08,0.3784727847451e+01,0.7834121070590e+01],
             [0.9260204034710e-08,0.2223352487601e+01,0.2787043132925e+01],
             [0.7320252888486e-08,0.1288694636874e+01,0.6282655592598e+01],
             [0.7319785780946e-08,0.5359869567774e+01,0.6283496108294e+01],
             [0.7147219933778e-08,0.5516616675856e+01,0.1725663147538e+02],
             [0.7946502829878e-08,0.2630459984567e+01,0.1241073141809e+02],
             [0.9001711808932e-08,0.2849815827227e+01,0.6281591679874e+01],
             [0.8994041507257e-08,0.3795244450750e+01,0.6284560021018e+01],
             [0.8298582787358e-08,0.5236413127363e+00,0.1241658836951e+02],
             [0.8526596520710e-08,0.4794605424426e+01,0.1098419223922e+02],
             [0.8209822103197e-08,0.1578752370328e+01,0.1096996532989e+02],
             [0.6357049861094e-08,0.5708926113761e+01,0.1596186371003e+01],
             [0.7370473179049e-08,0.3842402530241e+01,0.4061219149443e+01],
             [0.7232154664726e-08,0.3067548981535e+01,0.1610006857377e+03],
             [0.6328765494903e-08,0.1313930030069e+01,0.1193336791622e+02],
             [0.8030064908595e-08,0.3488500408886e+01,0.8460828644453e+00],
             [0.6275464259232e-08,0.1532061626198e+01,0.8531963191132e+00],
             [0.7051897446325e-08,0.3285859929993e+01,0.5849364236221e+01],
             [0.6161593705428e-08,0.1477341999464e+01,0.5573142801433e+01],
             [0.7754683957278e-08,0.1586118663096e+01,0.8662240327241e+01],
             [0.5889928990701e-08,0.1304887868803e+01,0.1232342296471e+02],
             [0.5705756047075e-08,0.4555333589350e+01,0.1258692712880e+02],
             [0.5964178808332e-08,0.3001762842062e+01,0.5333900173445e+01],
             [0.6712446027467e-08,0.4886780007595e+01,0.1171295538178e+02],
             [0.5941809275464e-08,0.4701509603824e+01,0.9779108567966e+01],
             [0.5466993627395e-08,0.4588357817278e+01,0.1884211409667e+02],
             [0.6340512090980e-08,0.1164543038893e+01,0.5217580628120e+02],
             [0.6325505710045e-08,0.3919171259645e+01,0.1041998632314e+02],
             [0.6164789509685e-08,0.2143828253542e+01,0.6151533897323e+01],
             [0.5263330812430e-08,0.6066564434241e+01,0.1885275071096e+02],
             [0.5597087780221e-08,0.2926316429472e+01,0.4337116142245e+00],
             [0.5396556236817e-08,0.3244303591505e+01,0.6286362197481e+01],
             [0.5396615148223e-08,0.3404304703662e+01,0.6279789503410e+01],
             [0.7091832443341e-08,0.8532377803192e+00,0.4907302013889e+01],
             [0.6572352589782e-08,0.4901966774419e+01,0.1176433076753e+02],
             [0.5960236060795e-08,0.1874672315797e+01,0.1422690933580e-01],
             [0.5125480043511e-08,0.3735726064334e+01,0.1245594543367e+02],
             [0.5928241866410e-08,0.4502033899935e+01,0.6414617803568e+01],
             [0.5249600357424e-08,0.4372334799878e+01,0.1151388321134e+02],
             [0.6059171276087e-08,0.2581617302908e+01,0.6062663316000e+01],
             [0.5295235081662e-08,0.2974811513158e+01,0.3496032717521e+01],
             [0.5820561875933e-08,0.1796073748244e+00,0.2838593341516e+00],
             [0.4754696606440e-08,0.1981998136973e+01,0.3104930017775e+01],
             [0.6385053548955e-08,0.2559174171605e+00,0.6133512519065e+01],
             [0.6589828273941e-08,0.2750967106776e+01,0.4087944051283e+02],
             [0.5383376567189e-08,0.6325947523578e+00,0.2248384854122e+02],
             [0.5928941683538e-08,0.1672304519067e+01,0.1581959461667e+01],
             [0.4816060709794e-08,0.3512566172575e+01,0.9388005868221e+01],
             [0.6003381586512e-08,0.5610932219189e+01,0.5326786718777e+01],
             [0.5504225393105e-08,0.4037501131256e+01,0.6503488384892e+01],
             [0.5353772620129e-08,0.6122774968240e+01,0.1735668374386e+03],
             [0.5786253768544e-08,0.5527984999515e+01,0.1350651127443e+00],
             [0.5065706702002e-08,0.9980765573624e+00,0.1248988586463e+02],
             [0.5972838885276e-08,0.6044489493203e+01,0.2673594526851e+02],
             [0.5323585877961e-08,0.3924265998147e+01,0.4171425416666e+01],
             [0.5210772682858e-08,0.6220111376901e+01,0.2460261242967e+02],
             [0.4726549040535e-08,0.3716043206862e+01,0.7232251527446e+01],
             [0.6029425105059e-08,0.8548704071116e+00,0.3227113045244e+03],
             [0.4481542826513e-08,0.1426925072829e+01,0.5547199253223e+01],
             [0.5836024505068e-08,0.7135651752625e-01,0.7285056171570e+02],
             [0.4137046613272e-08,0.5330767643283e+01,0.1087398597200e+02],
             [0.5171977473924e-08,0.4494262335353e+00,0.1884570439172e+02],
             [0.5694429833732e-08,0.2952369582215e+01,0.9723862754494e+02],
             [0.4009158925298e-08,0.3500003416535e+01,0.6244942932314e+01],
             [0.4784939596873e-08,0.6196709413181e+01,0.2929661536378e+02],
             [0.3983725022610e-08,0.5103690031897e+01,0.4274518229222e+01],
             [0.3870535232462e-08,0.3187569587401e+01,0.6321208768577e+01],
             [0.5140501213951e-08,0.1668924357457e+01,0.1232032006293e+02],
             [0.3849034819355e-08,0.4445722510309e+01,0.1726726808967e+02],
             [0.4002383075060e-08,0.5226224152423e+01,0.7018952447668e+01],
             [0.3890719543549e-08,0.4371166550274e+01,0.1491901785440e+02],
             [0.4887084607881e-08,0.5973556689693e+01,0.1478866649112e+01],
             [0.3739939287592e-08,0.2089084714600e+01,0.6922973089781e+01],
             [0.5031925918209e-08,0.4658371936827e+01,0.1715706182245e+02],
             [0.4387748764954e-08,0.4825580552819e+01,0.2331413144044e+03],
             [0.4147398098865e-08,0.3739003524998e+01,0.1376059875786e+02],
             [0.3719089993586e-08,0.1148941386536e+01,0.6297302759782e+01],
             [0.3934238461056e-08,0.1559893008343e+01,0.7872148766781e+01],
             [0.3672471375622e-08,0.5516145383612e+01,0.6268848941110e+01],
             [0.3768911277583e-08,0.6116053700563e+01,0.4157198507331e+01],
             [0.4033388417295e-08,0.5076821746017e+01,0.1567108171867e+02],
             [0.3764194617832e-08,0.8164676232075e+00,0.3185192151914e+01],
             [0.4840628226284e-08,0.1360479453671e+01,0.1252801878276e+02],
             [0.4949443923785e-08,0.2725622229926e+01,0.1617106187867e+03],
             [0.4117393089971e-08,0.6054459628492e+00,0.5642198095270e+01],
             [0.3925754020428e-08,0.8570462135210e+00,0.2139354194808e+02],
             [0.3630551757923e-08,0.3552067338279e+01,0.6294805223347e+01],
             [0.3627274802357e-08,0.3096565085313e+01,0.6271346477544e+01],
             [0.3806143885093e-08,0.6367751709777e+00,0.1725304118033e+02],
             [0.4433254641565e-08,0.4848461503937e+01,0.7445550607224e+01],
             [0.3712319846576e-08,0.1331950643655e+01,0.4194847048887e+00],
             [0.3849847534783e-08,0.4958368297746e+00,0.9562891316684e+00],
             [0.3483955430165e-08,0.2237215515707e+01,0.1161697602389e+02],
             [0.3961912730982e-08,0.3332402188575e+01,0.2277943724828e+02],
             [0.3419978244481e-08,0.5785600576016e+01,0.1362553364512e+02],
             [0.3329417758177e-08,0.9812676559709e-01,0.1685848245639e+02],
             [0.4207206893193e-08,0.9494780468236e+00,0.2986433403208e+02],
             [0.3268548976410e-08,0.1739332095686e+00,0.5749861718712e+01],
             [0.3321880082685e-08,0.1423354800666e+01,0.6279143387820e+01],
             [0.4503173010852e-08,0.2314972675293e+00,0.1385561574497e+01],
             [0.4316599090954e-08,0.1012646782616e+00,0.4176041334900e+01],
             [0.3283493323850e-08,0.5233306881265e+01,0.6287008313071e+01],
             [0.3164033542343e-08,0.4005597257511e+01,0.2099539292909e+02],
             [0.4159720956725e-08,0.5365676242020e+01,0.5905702259363e+01],
             [0.3565176892217e-08,0.4284440620612e+01,0.3932462625300e-02],
             [0.3514440950221e-08,0.4270562636575e+01,0.7335344340001e+01],
             [0.3540596871909e-08,0.5953553201060e+01,0.1234573916645e+02],
             [0.2960769905118e-08,0.1115180417718e+01,0.2670964694522e+02],
             [0.2962213739684e-08,0.3863811918186e+01,0.6408777551755e+00],
             [0.3883556700251e-08,0.1268617928302e+01,0.6660449441528e+01],
             [0.2919225516346e-08,0.4908605223265e+01,0.1375773836557e+01],
             [0.3115158863370e-08,0.3744519976885e+01,0.3802769619140e-01],
             [0.4099438144212e-08,0.4173244670532e+01,0.4480965020977e+02],
             [0.2899531858964e-08,0.5910601428850e+01,0.2059724391010e+02],
             [0.3289733429855e-08,0.2488050078239e+01,0.1081813534213e+02],
             [0.3933075612875e-08,0.1122363652883e+01,0.3773735910827e+00],
             [0.3021403764467e-08,0.4951973724904e+01,0.2982630633589e+02],
             [0.2798598949757e-08,0.5117057845513e+01,0.1937891852345e+02],
             [0.3397421302707e-08,0.6104159180476e+01,0.6923953605621e+01],
             [0.3720398002179e-08,0.1184933429829e+01,0.3066615496545e+02],
             [0.3598484186267e-08,0.3505282086105e+01,0.6147450479709e+01],
             [0.3694594027310e-08,0.2286651088141e+01,0.2636725487657e+01],
             [0.2680444152969e-08,0.1871816775482e+00,0.6816289982179e+01],
             [0.3497574865641e-08,0.3143251755431e+01,0.6418701221183e+01],
             [0.3130274129494e-08,0.2462167316018e+01,0.1235996607578e+02],
             [0.3241119069551e-08,0.4256374004686e+01,0.1652265972112e+02],
             [0.2601960842061e-08,0.4970362941425e+01,0.1045450126711e+02],
             [0.2690601527504e-08,0.2372657824898e+01,0.3163918923335e+00],
             [0.2908688152664e-08,0.4232652627721e+01,0.2828699048865e+02],
             [0.3120456131875e-08,0.3925747001137e+00,0.2195415756911e+02],
             [0.3148855423384e-08,0.3093478330445e+01,0.1172006883645e+02],
             [0.3051044261017e-08,0.5560948248212e+01,0.6055599646783e+01],
             [0.2826006876660e-08,0.5072790310072e+01,0.5120601093667e+01],
             [0.3100034191711e-08,0.4998530231096e+01,0.1799603123222e+02],
             [0.2398771640101e-08,0.2561739802176e+01,0.6255674361143e+01],
             [0.2384002842728e-08,0.4087420284111e+01,0.6310477339748e+01],
             [0.2842146517568e-08,0.2515048217955e+01,0.5469525544182e+01],
             [0.2847674371340e-08,0.5235326497443e+01,0.1034429499989e+02],
             [0.2903722140764e-08,0.1088200795797e+01,0.6510552054109e+01],
             [0.3187610710605e-08,0.4710624424816e+01,0.1693792562116e+03],
             [0.3048869992813e-08,0.2857975896445e+00,0.8390110365991e+01],
             [0.2860216950984e-08,0.2241619020815e+01,0.2243449970715e+00],
             [0.2701117683113e-08,0.6651573305272e-01,0.6129297044991e+01],
             [0.2509891590152e-08,0.1285135324585e+01,0.1044027435778e+02],
             [0.2623200252223e-08,0.2981229834530e+00,0.6436854655901e+01],
             [0.2622541669202e-08,0.6122470726189e+01,0.9380959548977e+01],
             [0.2818435667099e-08,0.4251087148947e+01,0.5934151399930e+01],
             [0.2365196797465e-08,0.3465070460790e+01,0.2470570524223e+02],
             [0.2358704646143e-08,0.5791603815350e+01,0.8671969964381e+01],
             [0.2388299481390e-08,0.4142483772941e+01,0.7096626156709e+01],
             [0.1996041217224e-08,0.2101901889496e+01,0.1727188400790e+02],
             [0.2687593060336e-08,0.1526689456959e+01,0.7075506709219e+02],
             [0.2618913670810e-08,0.2397684236095e+01,0.6632000300961e+01],
             [0.2571523050364e-08,0.5751929456787e+00,0.6206810014183e+01],
             [0.2582135006946e-08,0.5595464352926e+01,0.4873985990671e+02],
             [0.2372530190361e-08,0.5092689490655e+01,0.1590676413561e+02],
             [0.2357178484712e-08,0.4444363527851e+01,0.3097883698531e+01],
             [0.2451590394723e-08,0.3108251687661e+01,0.6612329252343e+00],
             [0.2370045949608e-08,0.2608133861079e+01,0.3459636466239e+02],
             [0.2268997267358e-08,0.3639717753384e+01,0.2844914056730e-01],
             [0.1731432137906e-08,0.1741898445707e+00,0.2019909489111e+02],
             [0.1629869741622e-08,0.3902225646724e+01,0.3035599730800e+02],
             [0.2206215801974e-08,0.4971131250731e+01,0.6281667977667e+01],
             [0.2205469554680e-08,0.1677462357110e+01,0.6284483723224e+01],
             [0.2148792362509e-08,0.4236259604006e+01,0.1980482729015e+02],
             [0.1873733657847e-08,0.5926814998687e+01,0.2876692439167e+02],
             [0.2026573758959e-08,0.4349643351962e+01,0.2449240616245e+02],
             [0.1807770325110e-08,0.5700940482701e+01,0.2045286941806e+02],
             [0.1881174408581e-08,0.6601286363430e+00,0.2358125818164e+02],
             [0.1368023671690e-08,0.2211098592752e+01,0.2473415438279e+02],
             [0.1720017916280e-08,0.4942488551129e+01,0.1679593901136e+03],
             [0.1702427665131e-08,0.1452233856386e+01,0.3338575901272e+03],
             [0.1414032510054e-08,0.5525357721439e+01,0.1624205518357e+03],
             [0.1652626045364e-08,0.4108794283624e+01,0.8956999012000e+02],
             [0.1642957769686e-08,0.7344335209984e+00,0.5267006960365e+02],
             [0.1614952403624e-08,0.3541213951363e+01,0.3332657872986e+02],
             [0.1535988291188e-08,0.4031094072151e+01,0.3852657435933e+02],
             [0.1593193738177e-08,0.4185136203609e+01,0.2282781046519e+03],
             [0.1074569126382e-08,0.1720485636868e+01,0.8397383534231e+02],
             [0.1074408214509e-08,0.2758613420318e+01,0.8401985929482e+02],
             [0.9700199670465e-09,0.4216686842097e+01,0.7826370942180e+02],
             [0.1258433517061e-08,0.2575068876639e+00,0.3115650189215e+03],
             [0.1240303229539e-08,0.4800844956756e+00,0.1784300471910e+03],
             [0.9018345948127e-09,0.3896756361552e+00,0.5886454391678e+02],
             [0.1135301432805e-08,0.3700805023550e+00,0.7842370451713e+02],
             [0.9215887951370e-09,0.4364579276638e+01,0.1014262087719e+03],
             [0.1055401054147e-08,0.2156564222111e+01,0.5660027930059e+02],
             [0.1008725979831e-08,0.5454015785234e+01,0.4245678405627e+02],
             [0.7217398104321e-09,0.1597772562175e+01,0.2457074661053e+03],
             [0.6912033134447e-09,0.5824090621461e+01,0.1679936946371e+03],
             [0.6833881523549e-09,0.3578778482835e+01,0.6053048899753e+02],
             [0.4887304205142e-09,0.3724362812423e+01,0.9656299901946e+02],
             [0.5173709754788e-09,0.5422427507933e+01,0.2442876000072e+03],
             [0.4671353097145e-09,0.2396106924439e+01,0.1435713242844e+03],
             [0.5652608439480e-09,0.2804028838685e+01,0.8365903305582e+02],
             [0.5604061331253e-09,0.1638816006247e+01,0.8433466158131e+02],
             [0.4712723365400e-09,0.8979003224474e+00,0.3164282286739e+03],
             [0.4909967465112e-09,0.3210426725516e+01,0.4059982187939e+03],
             [0.4771358267658e-09,0.5308027211629e+01,0.1805255418145e+03],
             [0.3943451445989e-09,0.2195145341074e+01,0.2568537517081e+03],
             [0.3952109120244e-09,0.5081189491586e+01,0.2449975330562e+03],
             [0.3788134594789e-09,0.4345171264441e+01,0.1568131045107e+03],
             [0.3738330190479e-09,0.2613062847997e+01,0.3948519331910e+03],
             [0.3099866678136e-09,0.2846760817689e+01,0.1547176098872e+03],
             [0.2002962716768e-09,0.4921360989412e+01,0.2268582385539e+03],
             [0.2198291338754e-09,0.1130360117454e+00,0.1658638954901e+03],
             [0.1491958330784e-09,0.4228195232278e+01,0.2219950288015e+03],
             [0.1475384076173e-09,0.3005721811604e+00,0.3052819430710e+03],
             [0.1661626624624e-09,0.7830125621203e+00,0.2526661704812e+03],
             [0.9015823460025e-10,0.3807792942715e+01,0.4171445043968e+03]]
    
    #Sun-to-Earth, T^1, X
    E1[0]=[[0.1234046326004e-05,0.0000000000000e+00,0.0000000000000e+00],
             [0.5150068824701e-06,0.6002664557501e+01,0.1256615170089e+02],
             [0.1290743923245e-07,0.5959437664199e+01,0.1884922755134e+02],
             [0.1068615564952e-07,0.2015529654209e+01,0.6283075850446e+01],
             [0.2079619142538e-08,0.1732960531432e+01,0.6279552690824e+01],
             [0.2078009243969e-08,0.4915604476996e+01,0.6286599010068e+01],
             [0.6206330058856e-09,0.3616457953824e+00,0.4705732307012e+01],
             [0.5989335313746e-09,0.3802607304474e+01,0.6256777527156e+01],
             [0.5958495663840e-09,0.2845866560031e+01,0.6309374173736e+01],
             [0.4866923261539e-09,0.5213203771824e+01,0.7755226100720e+00],
             [0.4267785823142e-09,0.4368189727818e+00,0.1059381944224e+01],
             [0.4610675141648e-09,0.1837249181372e-01,0.7860419393880e+01],
             [0.3626989993973e-09,0.2161590545326e+01,0.5753384878334e+01],
             [0.3563071194389e-09,0.1452631954746e+01,0.5884926831456e+01],
             [0.3557015642807e-09,0.4470593393054e+01,0.6812766822558e+01],
             [0.3210412089122e-09,0.5195926078314e+01,0.6681224869435e+01],
             [0.2875473577986e-09,0.5916256610193e+01,0.2513230340178e+02],
             [0.2842913681629e-09,0.1149902426047e+01,0.6127655567643e+01],
             [0.2751248215916e-09,0.5502088574662e+01,0.6438496133249e+01],
             [0.2481432881127e-09,0.2921989846637e+01,0.5486777812467e+01],
             [0.2059885976560e-09,0.3718070376585e+01,0.7079373888424e+01],
             [0.2015522342591e-09,0.5979395259740e+01,0.6290189305114e+01],
             [0.1995364084253e-09,0.6772087985494e+00,0.6275962395778e+01],
             [0.1957436436943e-09,0.2899210654665e+01,0.5507553240374e+01],
             [0.1651609818948e-09,0.6228206482192e+01,0.1150676975667e+02],
             [0.1822980550699e-09,0.1469348746179e+01,0.1179062909082e+02],
             [0.1675223159760e-09,0.3813910555688e+01,0.7058598460518e+01],
             [0.1706491764745e-09,0.3004380506684e+00,0.7113454667900e-02],
             [0.1392952362615e-09,0.1440393973406e+01,0.7962980379786e+00],
             [0.1209868266342e-09,0.4150425791727e+01,0.4694002934110e+01],
             [0.1009827202611e-09,0.3290040429843e+01,0.3738761453707e+01],
             [0.1047261388602e-09,0.4229590090227e+01,0.6282095334605e+01],
             [0.1047006652004e-09,0.2418967680575e+01,0.6284056366286e+01],
             [0.9609993143095e-10,0.4627943659201e+01,0.6069776770667e+01],
             [0.9590900593873e-10,0.1894393939924e+01,0.4136910472696e+01],
             [0.9146249188071e-10,0.2010647519562e+01,0.6496374930224e+01],
             [0.8545274480290e-10,0.5529846956226e-01,0.1194447056968e+01],
             [0.8224377881194e-10,0.1254304102174e+01,0.1589072916335e+01],
             [0.6183529510410e-10,0.3360862168815e+01,0.8827390247185e+01],
             [0.6259255147141e-10,0.4755628243179e+01,0.8429241228195e+01],
             [0.5539291694151e-10,0.5371746955142e+01,0.4933208510675e+01],
             [0.7328259466314e-10,0.4927699613906e+00,0.4535059491685e+01],
             [0.6017835843560e-10,0.5776682001734e-01,0.1255903824622e+02],
             [0.7079827775243e-10,0.4395059432251e+01,0.5088628793478e+01],
             [0.5170358878213e-10,0.5154062619954e+01,0.1176985366291e+02],
             [0.4872301838682e-10,0.6289611648973e+00,0.6040347114260e+01],
             [0.5249869411058e-10,0.5617272046949e+01,0.3154687086868e+01],
             [0.4716172354411e-10,0.3965901800877e+01,0.5331357529664e+01],
             [0.4871214940964e-10,0.4627507050093e+01,0.1256967486051e+02],
             [0.4598076850751e-10,0.6023631226459e+01,0.6525804586632e+01],
             [0.4562196089485e-10,0.4138562084068e+01,0.3930209696940e+01],
             [0.4325493872224e-10,0.1330845906564e+01,0.7632943190217e+01],
             [0.5673781176748e-10,0.2558752615657e+01,0.5729506548653e+01],
             [0.3961436642503e-10,0.2728071734630e+01,0.7234794171227e+01],
             [0.5101868209058e-10,0.4113444965144e+01,0.6836645152238e+01],
             [0.5257043167676e-10,0.6195089830590e+01,0.8031092209206e+01],
             [0.5076613989393e-10,0.2305124132918e+01,0.7477522907414e+01],
             [0.3342169352778e-10,0.5415998155071e+01,0.1097707878456e+02],
             [0.3545881983591e-10,0.3727160564574e+01,0.4164311961999e+01],
             [0.3364063738599e-10,0.2901121049204e+00,0.1137170464392e+02],
             [0.3357039670776e-10,0.1652229354331e+01,0.5223693906222e+01],
             [0.4307412268687e-10,0.4938909587445e+01,0.1592596075957e+01],
             [0.3405769115435e-10,0.2408890766511e+01,0.3128388763578e+01],
             [0.3001926198480e-10,0.4862239006386e+01,0.1748016358760e+01],
             [0.2778264787325e-10,0.5241168661353e+01,0.7342457794669e+01],
             [0.2676159480666e-10,0.3423593942199e+01,0.2146165377750e+01],
             [0.2954273399939e-10,0.1881721265406e+01,0.5368044267797e+00],
             [0.3309362888795e-10,0.1931525677349e+01,0.8018209333619e+00],
             [0.2810283608438e-10,0.2414659495050e+01,0.5225775174439e+00],
             [0.3378045637764e-10,0.4238019163430e+01,0.1554202828031e+00],
             [0.2558134979840e-10,0.1828225235805e+01,0.5230807360890e+01],
             [0.2273755578447e-10,0.5858184283998e+01,0.7084896783808e+01],
             [0.2294176037690e-10,0.4514589779057e+01,0.1726015463500e+02],
             [0.2533506099435e-10,0.2355717851551e+01,0.5216580451554e+01],
             [0.2716685375812e-10,0.2221003625100e+01,0.8635942003952e+01],
             [0.2419043435198e-10,0.5955704951635e+01,0.4690479774488e+01],
             [0.2521232544812e-10,0.1395676848521e+01,0.5481254917084e+01],
             [0.2630195021491e-10,0.5727468918743e+01,0.2629832328990e-01],
             [0.2548395840944e-10,0.2628351859400e-03,0.1349867339771e+01]]
    
    #Sun-to-Earth, T^2, X
    E2[0]=[[-0.4143818297913e-10,0.0000000000000e+00,0.0000000000000e+00],
             [0.2171497694435e-10,0.4398225628264e+01,0.1256615170089e+02],
             [0.9845398442516e-11,0.2079720838384e+00,0.6283075850446e+01],
             [0.9256833552682e-12,0.4191264694361e+01,0.1884922755134e+02],
             [0.1022049384115e-12,0.5381133195658e+01,0.8399684731857e+02]]
    
    #Sun-to-Earth, T^0, Y
    E0[1]=[[0.9998921098898e+00,0.1826583913846e+00,0.6283075850446e+01],
             [-0.2442700893735e-01,0.0000000000000e+00,0.0000000000000e+00],
             [0.8352929742915e-02,0.1395277998680e+00,0.1256615170089e+02],
             [0.1046697300177e-03,0.9641423109763e-01,0.1884922755134e+02],
             [0.3110841876663e-04,0.5381140401712e+01,0.8399684731857e+02],
             [0.2570269094593e-04,0.5301016407128e+01,0.5296909721118e+00],
             [0.2147389623610e-04,0.2662510869850e+01,0.1577343543434e+01],
             [0.1680344384050e-04,0.5207904119704e+01,0.6279552690824e+01],
             [0.1679117312193e-04,0.4582187486968e+01,0.6286599010068e+01],
             [0.1440512068440e-04,0.1900688517726e+01,0.2352866153506e+01],
             [0.1135139664999e-04,0.5273108538556e+01,0.5223693906222e+01],
             [0.9345482571018e-05,0.4503047687738e+01,0.1203646072878e+02],
             [0.9007418719568e-05,0.1605621059637e+01,0.1021328554739e+02],
             [0.5671536712314e-05,0.5812849070861e+00,0.1059381944224e+01],
             [0.7451401861666e-05,0.2807346794836e+01,0.3981490189893e+00],
             [0.6393470057114e-05,0.6029224133855e+01,0.5753384878334e+01],
             [0.6814275881697e-05,0.6472990145974e+00,0.4705732307012e+01],
             [0.6113705628887e-05,0.3813843419700e+01,0.6812766822558e+01],
             [0.4503851367273e-05,0.4527804370996e+01,0.5884926831456e+01],
             [0.4522249141926e-05,0.5991783029224e+01,0.6256777527156e+01],
             [0.4501794307018e-05,0.3798703844397e+01,0.6309374173736e+01],
             [0.5514927480180e-05,0.3961257833388e+01,0.5507553240374e+01],
             [0.4062862799995e-05,0.5256247296369e+01,0.6681224869435e+01],
             [0.5414900429712e-05,0.5499032014097e+01,0.7755226100720e+00],
             [0.5463153987424e-05,0.6173092454097e+01,0.1414349524433e+02],
             [0.5071611859329e-05,0.2870244247651e+01,0.7860419393880e+01],
             [0.2195112094455e-05,0.2952338617201e+01,0.1150676975667e+02],
             [0.2279139233919e-05,0.5951775132933e+01,0.7058598460518e+01],
             [0.2278386100876e-05,0.4845456398785e+01,0.4694002934110e+01],
             [0.2559088003308e-05,0.6945321117311e+00,0.1216800268190e+02],
             [0.2561079286856e-05,0.6167224608301e+01,0.7099330490126e+00],
             [0.1792755796387e-05,0.1400122509632e+01,0.7962980379786e+00],
             [0.1818715656502e-05,0.4703347611830e+01,0.6283142985870e+01],
             [0.1818744924791e-05,0.5086748900237e+01,0.6283008715021e+01],
             [0.1554518791390e-05,0.5331008042713e-01,0.2513230340178e+02],
             [0.2063265737239e-05,0.4283680484178e+01,0.1179062909082e+02],
             [0.1497613520041e-05,0.6074207826073e+01,0.5486777812467e+01],
             [0.2000617940427e-05,0.2501426281450e+01,0.1778984560711e+02],
             [0.1289731195580e-05,0.3646340599536e+01,0.7079373888424e+01],
             [0.1282657998934e-05,0.3232864804902e+01,0.3738761453707e+01],
             [0.1528915968658e-05,0.5581433416669e+01,0.2132990797783e+00],
             [0.1187304098432e-05,0.5453576453694e+01,0.9437762937313e+01],
             [0.7842782928118e-06,0.2823953922273e+00,0.8827390247185e+01],
             [0.7352892280868e-06,0.1124369580175e+01,0.1589072916335e+01],
             [0.6570189360797e-06,0.2089154042840e+01,0.1176985366291e+02],
             [0.6324967590410e-06,0.6704855581230e+00,0.6262300422539e+01],
             [0.6298289872283e-06,0.2836414855840e+01,0.6303851278352e+01],
             [0.6476686465855e-06,0.4852433866467e+00,0.7113454667900e-02],
             [0.8587034651234e-06,0.1453511005668e+01,0.1672837615881e+03],
             [0.8068948788113e-06,0.9224087798609e+00,0.6069776770667e+01],
             [0.8353786011661e-06,0.4631707184895e+01,0.3340612434717e+01],
             [0.6009324532132e-06,0.1829498827726e+01,0.4136910472696e+01],
             [0.7558158559566e-06,0.2588596800317e+01,0.6496374930224e+01],
             [0.5809279504503e-06,0.5516818853476e+00,0.1097707878456e+02],
             [0.5374131950254e-06,0.6275674734960e+01,0.1194447056968e+01],
             [0.5711160507326e-06,0.1091905956872e+01,0.6282095334605e+01],
             [0.5710183170746e-06,0.2415001635090e+01,0.6284056366286e+01],
             [0.5144373590610e-06,0.6020336443438e+01,0.6290189305114e+01],
             [0.5103108927267e-06,0.3775634564605e+01,0.6275962395778e+01],
             [0.4960654697891e-06,0.1073450946756e+01,0.6127655567643e+01],
             [0.4786385689280e-06,0.2431178012310e+01,0.6438496133249e+01],
             [0.6109911263665e-06,0.5343356157914e+01,0.3154687086868e+01],
             [0.4839898944024e-06,0.5830833594047e-01,0.8018209333619e+00],
             [0.4734822623919e-06,0.4536080134821e+01,0.3128388763578e+01],
             [0.4834741473290e-06,0.2585090489754e+00,0.7084896783808e+01],
             [0.5134858581156e-06,0.4213317172603e+01,0.1235285262111e+02],
             [0.5064004264978e-06,0.4814418806478e+00,0.1185621865188e+02],
             [0.3753476772761e-06,0.1599953399788e+01,0.8429241228195e+01],
             [0.4935264014283e-06,0.2157417556873e+01,0.2544314396739e+01],
             [0.3950929600897e-06,0.3359394184254e+01,0.5481254917084e+01],
             [0.4895849789777e-06,0.5165704376558e+01,0.9225539266174e+01],
             [0.4215241688886e-06,0.2065368800993e+01,0.1726015463500e+02],
             [0.3796773731132e-06,0.1468606346612e+01,0.4265981595566e+00],
             [0.3114178142515e-06,0.3615638079474e+01,0.2146165377750e+01],
             [0.3260664220838e-06,0.4417134922435e+01,0.4164311961999e+01],
             [0.3976996123008e-06,0.4700866883004e+01,0.5856477690889e+01],
             [0.2801459672924e-06,0.4538902060922e+01,0.1256967486051e+02],
             [0.3638931868861e-06,0.1334197991475e+01,0.1807370494127e+02],
             [0.2487013269476e-06,0.3749275558275e+01,0.2629832328990e-01],
             [0.3034165481994e-06,0.4236622030873e+00,0.4535059491685e+01],
             [0.2676278825586e-06,0.5970848007811e+01,0.3930209696940e+01],
             [0.2764903818918e-06,0.5194636754501e+01,0.1256262854127e+02],
             [0.2485149930507e-06,0.1002434207846e+01,0.5088628793478e+01],
             [0.2199305540941e-06,0.3066773098403e+01,0.1255903824622e+02],
             [0.2571106500435e-06,0.7588312459063e+00,0.1336797263425e+02],
             [0.2049751817158e-06,0.3444977434856e+01,0.1137170464392e+02],
             [0.2599707296297e-06,0.1873128542205e+01,0.7143069561767e+02],
             [0.1785018072217e-06,0.5015891306615e+01,0.1748016358760e+01],
             [0.2324833891115e-06,0.4618271239730e+01,0.1831953657923e+02],
             [0.1709711119545e-06,0.5300003455669e+01,0.4933208510675e+01],
             [0.2107159351716e-06,0.2229819815115e+01,0.7477522907414e+01],
             [0.1750333080295e-06,0.6161485880008e+01,0.1044738781244e+02],
             [0.2000598210339e-06,0.2967357299999e+01,0.8031092209206e+01],
             [0.1380920248681e-06,0.3027007923917e+01,0.8635942003952e+01],
             [0.1412460470299e-06,0.6037597163798e+01,0.2942463415728e+01],
             [0.1888459803001e-06,0.8561476243374e+00,0.1561374759853e+03],
             [0.1788370542585e-06,0.4869736290209e+01,0.1592596075957e+01],
             [0.1360893296167e-06,0.3626411886436e+01,0.1309584267300e+02],
             [0.1506846530160e-06,0.1550975377427e+01,0.1649636139783e+02],
             [0.1800913376176e-06,0.2075826033190e+01,0.1729818233119e+02],
             [0.1436261390649e-06,0.6148876420255e+01,0.2042657109477e+02],
             [0.1220227114151e-06,0.4382583879906e+01,0.7632943190217e+01],
             [0.1337883603592e-06,0.2036644327361e+01,0.1213955354133e+02],
             [0.1159326650738e-06,0.3892276994687e+01,0.5331357529664e+01],
             [0.1352853128569e-06,0.1447950649744e+01,0.1673046366289e+02],
             [0.1433408296083e-06,0.4457854692961e+01,0.7342457794669e+01],
             [0.1234701666518e-06,0.1538818147151e+01,0.6279485555400e+01],
             [0.1234027192007e-06,0.1968523220760e+01,0.6286666145492e+01],
             [0.1244024091797e-06,0.5779803499985e+01,0.1511046609763e+02],
             [0.1097934945516e-06,0.6210975221388e+00,0.1098880815746e+02],
             [0.1254611329856e-06,0.2591963807998e+01,0.1572083878776e+02],
             [0.1158247286784e-06,0.2483612812670e+01,0.5729506548653e+01],
             [0.9039078252960e-07,0.3857554579796e+01,0.9623688285163e+01],
             [0.9108024978836e-07,0.5826368512984e+01,0.7234794171227e+01],
             [0.8887068108436e-07,0.3475694573987e+01,0.6148010737701e+01],
             [0.8632374035438e-07,0.3059070488983e-01,0.6418140963190e+01],
             [0.7893186992967e-07,0.1583194837728e+01,0.2118763888447e+01],
             [0.8297650201172e-07,0.8519770534637e+00,0.1471231707864e+02],
             [0.1019759578988e-06,0.1319598738732e+00,0.1349867339771e+01],
             [0.1010037696236e-06,0.9937860115618e+00,0.6836645152238e+01],
             [0.1047727548266e-06,0.1382138405399e+01,0.5999216516294e+01],
             [0.7351993881086e-07,0.3833397851735e+01,0.6040347114260e+01],
             [0.9868771092341e-07,0.2124913814390e+01,0.6566935184597e+01],
             [0.7007321959390e-07,0.5946305343763e+01,0.6525804586632e+01],
             [0.6861411679709e-07,0.4574654977089e+01,0.7238675589263e+01],
             [0.7554519809614e-07,0.5949232686844e+01,0.1253985337760e+02],
             [0.9541880448335e-07,0.3495242990564e+01,0.2122839202813e+02],
             [0.7185606722155e-07,0.4310113471661e+01,0.6245048154254e+01],
             [0.7131360871710e-07,0.5480309323650e+01,0.6321103546637e+01],
             [0.6651142021039e-07,0.5411097713654e+01,0.5327476111629e+01],
             [0.8538618213667e-07,0.1827849973951e+01,0.1101510648075e+02],
             [0.8634954288044e-07,0.5443584943349e+01,0.5643178611111e+01],
             [0.7449415051484e-07,0.2011535459060e+01,0.5368044267797e+00],
             [0.7421047599169e-07,0.3464562529249e+01,0.2354323048545e+02],
             [0.6140694354424e-07,0.5657556228815e+01,0.1296430071988e+02],
             [0.6353525143033e-07,0.3463816593821e+01,0.1990745094947e+01],
             [0.6221964013447e-07,0.1532259498697e+01,0.9517183207817e+00],
             [0.5852480257244e-07,0.1375396598875e+01,0.9555997388169e+00],
             [0.6398637498911e-07,0.2405645801972e+01,0.2407292145756e+02],
             [0.7039744069878e-07,0.5397541799027e+01,0.5225775174439e+00],
             [0.6977997694382e-07,0.4762347105419e+01,0.1097355562493e+02],
             [0.7460629558396e-07,0.2711944692164e+01,0.2200391463820e+02],
             [0.5376577536101e-07,0.2352980430239e+01,0.1431416805965e+02],
             [0.7530607893556e-07,0.1943940180699e+01,0.1842262939178e+02],
             [0.6822928971605e-07,0.4337651846959e+01,0.1554202828031e+00],
             [0.6220772380094e-07,0.6716871369278e+00,0.1845107853235e+02],
             [0.6586950799043e-07,0.2229714460505e+01,0.5216580451554e+01],
             [0.5873800565771e-07,0.7627013920580e+00,0.6398972393349e+00],
             [0.6264346929745e-07,0.6202785478961e+00,0.6277552955062e+01],
             [0.6257929115669e-07,0.2886775596668e+01,0.6288598745829e+01],
             [0.5343536033409e-07,0.1977241012051e+01,0.4690479774488e+01],
             [0.5587849781714e-07,0.1922923484825e+01,0.1551045220144e+01],
             [0.6905100845603e-07,0.3570757164631e+01,0.1030928125552e+00],
             [0.6178957066649e-07,0.5197558947765e+01,0.5230807360890e+01],
             [0.6187270224331e-07,0.8193497368922e+00,0.5650292065779e+01],
             [0.5385664291426e-07,0.5406336665586e+01,0.7771377146812e+02],
             [0.6329363917926e-07,0.2837760654536e+01,0.2608790314060e+02],
             [0.4546018761604e-07,0.2933580297050e+01,0.5535693017924e+00],
             [0.6196091049375e-07,0.4157871494377e+01,0.8467247584405e+02],
             [0.6159555108218e-07,0.3211703561703e+01,0.2394243902548e+03],
             [0.4995340539317e-07,0.1459098102922e+01,0.4732030630302e+01],
             [0.5457031243572e-07,0.1430457676136e+01,0.6179983037890e+01],
             [0.4863461418397e-07,0.2196425916730e+01,0.9027992316901e+02],
             [0.5342947626870e-07,0.2086612890268e+01,0.6386168663001e+01],
             [0.5674296648439e-07,0.2760204966535e+01,0.6915859635113e+01],
             [0.4745783120161e-07,0.4245368971862e+01,0.6282970628506e+01],
             [0.4745676961198e-07,0.5544725787016e+01,0.6283181072386e+01],
             [0.4049796869973e-07,0.2213984363586e+01,0.6254626709878e+01],
             [0.4248333596940e-07,0.8075781952896e+00,0.7875671926403e+01],
             [0.4027178070205e-07,0.1293268540378e+01,0.6311524991013e+01],
             [0.4066543943476e-07,0.3986141175804e+01,0.3634620989887e+01],
             [0.4858863787880e-07,0.1276112738231e+01,0.5760498333002e+01],
             [0.5277398263530e-07,0.4916111741527e+01,0.2515860172507e+02],
             [0.4105635656559e-07,0.1725805864426e+01,0.6709674010002e+01],
             [0.4376781925772e-07,0.2243642442106e+01,0.6805653367890e+01],
             [0.3235827894693e-07,0.3614135118271e+01,0.1066495398892e+01],
             [0.3073244740308e-07,0.2460873393460e+01,0.5863591145557e+01],
             [0.3088609271373e-07,0.5678431771790e+01,0.9917696840332e+01],
             [0.3393022279836e-07,0.3814017477291e+01,0.1391601904066e+02],
             [0.3038686508802e-07,0.4660216229171e+01,0.1256621883632e+02],
             [0.4019677752497e-07,0.5906906243735e+01,0.1334167431096e+02],
             [0.3288834998232e-07,0.9536146445882e+00,0.1620077269078e+02],
             [0.3889973794631e-07,0.3942205097644e+01,0.7478166569050e-01],
             [0.3050438987141e-07,0.1624810271286e+01,0.1805292951336e+02],
             [0.3601142564638e-07,0.4030467142575e+01,0.6208294184755e+01],
             [0.3689015557141e-07,0.3648878818694e+01,0.5966683958112e+01],
             [0.3563471893565e-07,0.5749584017096e+01,0.6357857516136e+01],
             [0.2776183170667e-07,0.2630124187070e+01,0.3523159621801e-02],
             [0.2922350530341e-07,0.1790346403629e+01,0.1272157198369e+02],
             [0.3511076917302e-07,0.6142198301611e+01,0.6599467742779e+01],
             [0.3619351007632e-07,0.1432421386492e+01,0.6019991944201e+01],
             [0.2561254711098e-07,0.2302822475792e+01,0.1259245002418e+02],
             [0.2626903942920e-07,0.8660470994571e+00,0.6702560555334e+01],
             [0.2550187397083e-07,0.6069721995383e+01,0.1057540660594e+02],
             [0.2535873526138e-07,0.1079020331795e-01,0.3141537925223e+02],
             [0.3519786153847e-07,0.3809066902283e+01,0.2505706758577e+03],
             [0.3424651492873e-07,0.2075435114417e+01,0.6546159756691e+01],
             [0.2372676630861e-07,0.2057803120154e+01,0.2388894113936e+01],
             [0.2710980779541e-07,0.1510068488010e+01,0.1202934727411e+02],
             [0.3038710889704e-07,0.5043617528901e+01,0.1256608456547e+02],
             [0.2220364130585e-07,0.3694793218205e+01,0.1336244973887e+02],
             [0.3025880825460e-07,0.5450618999049e-01,0.2908881142201e+02],
             [0.2784493486864e-07,0.3381164084502e+01,0.1494531617769e+02],
             [0.2294414142438e-07,0.4382309025210e+01,0.6076890225335e+01],
             [0.2012723294724e-07,0.9142212256518e+00,0.6262720680387e+01],
             [0.2036357831958e-07,0.5676172293154e+01,0.4701116388778e+01],
             [0.2003474823288e-07,0.2592767977625e+01,0.6303431020504e+01],
             [0.2207144900109e-07,0.5404976271180e+01,0.6489261475556e+01],
             [0.2481664905135e-07,0.4373284587027e+01,0.1204357418345e+02],
             [0.2674949182295e-07,0.5859182188482e+01,0.4590910121555e+01],
             [0.2450554720322e-07,0.4555381557451e+01,0.1495633313810e+00],
             [0.2601975986457e-07,0.3933165584959e+01,0.1965104848470e+02],
             [0.2199860022848e-07,0.5227977189087e+01,0.1351787002167e+02],
             [0.2448121172316e-07,0.4858060353949e+01,0.1162474756779e+01],
             [0.1876014864049e-07,0.5690546553605e+01,0.6279194432410e+01],
             [0.1874513219396e-07,0.4099539297446e+01,0.6286957268481e+01],
             [0.2156380842559e-07,0.4382594769913e+00,0.1813929450232e+02],
             [0.1981691240061e-07,0.1829784152444e+01,0.4686889479442e+01],
             [0.2329992648539e-07,0.2836254278973e+01,0.1002183730415e+02],
             [0.1765184135302e-07,0.2803494925833e+01,0.4292330755499e+01],
             [0.2436368366085e-07,0.2836897959677e+01,0.9514313292143e+02],
             [0.2164089203889e-07,0.6127522446024e+01,0.6037244212485e+01],
             [0.1847755034221e-07,0.3683163635008e+01,0.2427287361862e+00],
             [0.1674798769966e-07,0.3316993867246e+00,0.1311972100268e+02],
             [0.2222542124356e-07,0.8294097805480e+00,0.1266924451345e+02],
             [0.2071074505925e-07,0.3659492220261e+01,0.6528907488406e+01],
             [0.1608224471835e-07,0.4774492067182e+01,0.1352175143971e+02],
             [0.1857583439071e-07,0.2873120597682e+01,0.8662240327241e+01],
             [0.1793018836159e-07,0.5282441177929e+00,0.6819880277225e+01],
             [0.1575391221692e-07,0.1320789654258e+01,0.1102062672231e+00],
             [0.1840132009557e-07,0.1917110916256e+01,0.6514761976723e+02],
             [0.1760917288281e-07,0.2972635937132e+01,0.5746271423666e+01],
             [0.1561779518516e-07,0.4372569261981e+01,0.6272439236156e+01],
             [0.1558687885205e-07,0.5416424926425e+01,0.6293712464735e+01],
             [0.1951359382579e-07,0.3094448898752e+01,0.2301353951334e+02],
             [0.1569144275614e-07,0.2802103689808e+01,0.1765478049437e+02],
             [0.1479130389462e-07,0.2136435020467e+01,0.2077542790660e-01],
             [0.1467828510764e-07,0.7072627435674e+00,0.1052268489556e+01],
             [0.1627627337440e-07,0.3947607143237e+01,0.6327837846670e+00],
             [0.1503498479758e-07,0.4079248909190e+01,0.7626583626240e-01],
             [0.1297967708237e-07,0.6269637122840e+01,0.1149965630200e+02],
             [0.1374416896634e-07,0.4175657970702e+01,0.6016468784579e+01],
             [0.1783812325219e-07,0.1476540547560e+01,0.3301902111895e+02],
             [0.1525884228756e-07,0.4653477715241e+01,0.9411464614024e+01],
             [0.1451067396763e-07,0.2573001128225e+01,0.1277945078067e+02],
             [0.1297713111950e-07,0.5612799618771e+01,0.6549682916313e+01],
             [0.1462784012820e-07,0.4189661623870e+01,0.1863592847156e+02],
             [0.1384185980007e-07,0.2656915472196e+01,0.2379164476796e+01],
             [0.1221497599801e-07,0.5612515760138e+01,0.1257326515556e+02],
             [0.1560574525896e-07,0.4783414317919e+01,0.1887552587463e+02],
             [0.1544598372036e-07,0.2694431138063e+01,0.1820933031200e+02],
             [0.1531678928696e-07,0.4105103489666e+01,0.2593412433514e+02],
             [0.1349321503795e-07,0.3082437194015e+00,0.5120601093667e+01],
             [0.1252030290917e-07,0.6124072334087e+01,0.6993008899458e+01],
             [0.1459243816687e-07,0.3733103981697e+01,0.3813291813120e-01],
             [0.1226103625262e-07,0.1267127706817e+01,0.2435678079171e+02],
             [0.1019449641504e-07,0.4367790112269e+01,0.1725663147538e+02],
             [0.1380789433607e-07,0.3387201768700e+01,0.2458316379602e+00],
             [0.1019453421658e-07,0.9204143073737e+00,0.6112403035119e+01],
             [0.1297929434405e-07,0.5786874896426e+01,0.1249137003520e+02],
             [0.9912677786097e-08,0.3164232870746e+01,0.6247047890016e+01],
             [0.9829386098599e-08,0.2586762413351e+01,0.6453748665772e+01],
             [0.1226807746104e-07,0.6239068436607e+01,0.5429879531333e+01],
             [0.1192691755997e-07,0.1867380051424e+01,0.6290122169689e+01],
             [0.9836499227081e-08,0.3424716293727e+00,0.6319103810876e+01],
             [0.9642862564285e-08,0.5661372990657e+01,0.8273820945392e+01],
             [0.1165184404862e-07,0.5768367239093e+01,0.1778273215245e+02],
             [0.1175794418818e-07,0.1657351222943e+01,0.6276029531202e+01],
             [0.1018948635601e-07,0.6458292350865e+00,0.1254537627298e+02],
             [0.9500383606676e-08,0.1054306140741e+01,0.1256517118505e+02],
             [0.1227512202906e-07,0.2505278379114e+01,0.2248384854122e+02],
             [0.9664792009993e-08,0.4289737277000e+01,0.6259197520765e+01],
             [0.9613285666331e-08,0.5500597673141e+01,0.6306954180126e+01],
             [0.1117906736211e-07,0.2361405953468e+01,0.1779695906178e+02],
             [0.9611378640782e-08,0.2851310576269e+01,0.2061856251104e+00],
             [0.8845354852370e-08,0.6208777705343e+01,0.1692165728891e+01],
             [0.1054046966600e-07,0.5413091423934e+01,0.2204125344462e+00],
             [0.1215539124483e-07,0.5613969479755e+01,0.8257698122054e+02],
             [0.9932460955209e-08,0.1106124877015e+01,0.1017725758696e+02],
             [0.8785804715043e-08,0.2869224476477e+01,0.9491756770005e+00],
             [0.8538084097562e-08,0.6159640899344e+01,0.6393282117669e+01],
             [0.8648994369529e-08,0.1374901198784e+01,0.4804209201333e+01],
             [0.1039063219067e-07,0.5171080641327e+01,0.1550861511662e+02],
             [0.8867983926439e-08,0.8317320304902e+00,0.3903911373650e+01],
             [0.8327495955244e-08,0.3605591969180e+01,0.6172869583223e+01],
             [0.9243088356133e-08,0.6114299196843e+01,0.6267823317922e+01],
             [0.9205657357835e-08,0.3675153683737e+01,0.6298328382969e+01],
             [0.1033269714606e-07,0.3313328813024e+01,0.5573142801433e+01],
             [0.8001706275552e-08,0.2019980960053e+01,0.2648454860559e+01],
             [0.9171858254191e-08,0.8992015524177e+00,0.1498544001348e+03],
             [0.1075327150242e-07,0.2898669963648e+01,0.3694923081589e+02],
             [0.9884866689828e-08,0.4946715904478e+01,0.1140367694411e+02],
             [0.9541835576677e-08,0.2371787888469e+01,0.1256713221673e+02],
             [0.7739903376237e-08,0.2213775190612e+01,0.7834121070590e+01],
             [0.7311962684106e-08,0.3429378787739e+01,0.1192625446156e+02],
             [0.9724904869624e-08,0.6195878564404e+01,0.2280573557157e+02],
             [0.9251628983612e-08,0.6511509527390e+00,0.2787043132925e+01],
             [0.7320763787842e-08,0.6001083639421e+01,0.6282655592598e+01],
             [0.7320296650962e-08,0.3789073265087e+01,0.6283496108294e+01],
             [0.7947032271039e-08,0.1059659582204e+01,0.1241073141809e+02],
             [0.9005277053115e-08,0.1280315624361e+01,0.6281591679874e+01],
             [0.8995601652048e-08,0.2224439106766e+01,0.6284560021018e+01],
             [0.8288040568796e-08,0.5234914433867e+01,0.1241658836951e+02],
             [0.6359381347255e-08,0.4137989441490e+01,0.1596186371003e+01],
             [0.8699572228626e-08,0.1758411009497e+01,0.6133512519065e+01],
             [0.6456797542736e-08,0.5919285089994e+01,0.1685848245639e+02],
             [0.7424573475452e-08,0.5414616938827e+01,0.4061219149443e+01],
             [0.7235671196168e-08,0.1496516557134e+01,0.1610006857377e+03],
             [0.8104015182733e-08,0.1919918242764e+01,0.8460828644453e+00],
             [0.8098576535937e-08,0.3819615855458e+01,0.3894181736510e+01],
             [0.6275292346625e-08,0.6244264115141e+01,0.8531963191132e+00],
             [0.6052432989112e-08,0.5037731872610e+00,0.1567108171867e+02],
             [0.5705651535817e-08,0.2984557271995e+01,0.1258692712880e+02],
             [0.5789650115138e-08,0.6087038140697e+01,0.1193336791622e+02],
             [0.5512132153377e-08,0.5855668994076e+01,0.1232342296471e+02],
             [0.7388890819102e-08,0.2443128574740e+01,0.4907302013889e+01],
             [0.5467593991798e-08,0.3017561234194e+01,0.1884211409667e+02],
             [0.6388519802999e-08,0.5887386712935e+01,0.5217580628120e+02],
             [0.6106777149944e-08,0.3483461059895e+00,0.1422690933580e-01],
             [0.7383420275489e-08,0.5417387056707e+01,0.2358125818164e+02],
             [0.5505208141738e-08,0.2848193644783e+01,0.1151388321134e+02],
             [0.6310757462877e-08,0.2349882520828e+01,0.1041998632314e+02],
             [0.6166904929691e-08,0.5728575944077e+00,0.6151533897323e+01],
             [0.5263442042754e-08,0.4495796125937e+01,0.1885275071096e+02],
             [0.5591828082629e-08,0.1355441967677e+01,0.4337116142245e+00],
             [0.5397051680497e-08,0.1673422864307e+01,0.6286362197481e+01],
             [0.5396992745159e-08,0.1833502206373e+01,0.6279789503410e+01],
             [0.6572913000726e-08,0.3331122065824e+01,0.1176433076753e+02],
             [0.5123421866413e-08,0.2165327142679e+01,0.1245594543367e+02],
             [0.5930495725999e-08,0.2931146089284e+01,0.6414617803568e+01],
             [0.6431797403933e-08,0.4134407994088e+01,0.1350651127443e+00],
             [0.5003182207604e-08,0.3805420303749e+01,0.1096996532989e+02],
             [0.5587731032504e-08,0.1082469260599e+01,0.6062663316000e+01],
             [0.5935263407816e-08,0.8384333678401e+00,0.5326786718777e+01],
             [0.4756019827760e-08,0.3552588749309e+01,0.3104930017775e+01],
             [0.6599951172637e-08,0.4320826409528e+01,0.4087944051283e+02],
             [0.5902606868464e-08,0.4811879454445e+01,0.5849364236221e+01],
             [0.5921147809031e-08,0.9942628922396e-01,0.1581959461667e+01],
             [0.5505382581266e-08,0.2466557607764e+01,0.6503488384892e+01],
             [0.5353771071862e-08,0.4551978748683e+01,0.1735668374386e+03],
             [0.5063282210946e-08,0.5710812312425e+01,0.1248988586463e+02],
             [0.5926120403383e-08,0.1333998428358e+01,0.2673594526851e+02],
             [0.5211016176149e-08,0.4649315360760e+01,0.2460261242967e+02],
             [0.5347075084894e-08,0.5512754081205e+01,0.4171425416666e+01],
             [0.4872609773574e-08,0.1308025299938e+01,0.5333900173445e+01],
             [0.4727711321420e-08,0.2144908368062e+01,0.7232251527446e+01],
             [0.6029426018652e-08,0.5567259412084e+01,0.3227113045244e+03],
             [0.4321485284369e-08,0.5230667156451e+01,0.9388005868221e+01],
             [0.4476406760553e-08,0.6134081115303e+01,0.5547199253223e+01],
             [0.5835268277420e-08,0.4783808492071e+01,0.7285056171570e+02],
             [0.5172183602748e-08,0.5161817911099e+01,0.1884570439172e+02],
             [0.5693571465184e-08,0.1381646203111e+01,0.9723862754494e+02],
             [0.4060634965349e-08,0.3876705259495e+00,0.4274518229222e+01],
             [0.3967398770473e-08,0.5029491776223e+01,0.3496032717521e+01],
             [0.3943754005255e-08,0.1923162955490e+01,0.6244942932314e+01],
             [0.4781323427824e-08,0.4633332586423e+01,0.2929661536378e+02],
             [0.3871483781204e-08,0.1616650009743e+01,0.6321208768577e+01],
             [0.5141741733997e-08,0.9817316704659e-01,0.1232032006293e+02],
             [0.4002385978497e-08,0.3656161212139e+01,0.7018952447668e+01],
             [0.4901092604097e-08,0.4404098713092e+01,0.1478866649112e+01],
             [0.3740932630345e-08,0.5181188732639e+00,0.6922973089781e+01],
             [0.4387283718538e-08,0.3254859566869e+01,0.2331413144044e+03],
             [0.5019197802033e-08,0.3086773224677e+01,0.1715706182245e+02],
             [0.3834931695175e-08,0.2797882673542e+01,0.1491901785440e+02],
             [0.3760413942497e-08,0.2892676280217e+01,0.1726726808967e+02],
             [0.3719717204628e-08,0.5861046025739e+01,0.6297302759782e+01],
             [0.4145623530149e-08,0.2168239627033e+01,0.1376059875786e+02],
             [0.3932788425380e-08,0.6271811124181e+01,0.7872148766781e+01],
             [0.3686377476857e-08,0.3936853151404e+01,0.6268848941110e+01],
             [0.3779077950339e-08,0.1404148734043e+01,0.4157198507331e+01],
             [0.4091334550598e-08,0.2452436180854e+01,0.9779108567966e+01],
             [0.3926694536146e-08,0.6102292739040e+01,0.1098419223922e+02],
             [0.4841000253289e-08,0.6072760457276e+01,0.1252801878276e+02],
             [0.4949340130240e-08,0.1154832815171e+01,0.1617106187867e+03],
             [0.3761557737360e-08,0.5527545321897e+01,0.3185192151914e+01],
             [0.3647396268188e-08,0.1525035688629e+01,0.6271346477544e+01],
             [0.3932405074189e-08,0.5570681040569e+01,0.2139354194808e+02],
             [0.3631322501141e-08,0.1981240601160e+01,0.6294805223347e+01],
             [0.4130007425139e-08,0.2050060880201e+01,0.2195415756911e+02],
             [0.4433905965176e-08,0.3277477970321e+01,0.7445550607224e+01],
             [0.3851814176947e-08,0.5210690074886e+01,0.9562891316684e+00],
             [0.3485807052785e-08,0.6653274904611e+00,0.1161697602389e+02],
             [0.3979772816991e-08,0.1767941436148e+01,0.2277943724828e+02],
             [0.3402607460500e-08,0.3421746306465e+01,0.1087398597200e+02],
             [0.4049993000926e-08,0.1127144787547e+01,0.3163918923335e+00],
             [0.3420511182382e-08,0.4214794779161e+01,0.1362553364512e+02],
             [0.3640772365012e-08,0.5324905497687e+01,0.1725304118033e+02],
             [0.3323037987501e-08,0.6135761838271e+01,0.6279143387820e+01],
             [0.4503141663637e-08,0.1802305450666e+01,0.1385561574497e+01],
             [0.4314560055588e-08,0.4812299731574e+01,0.4176041334900e+01],
             [0.3294226949110e-08,0.3657547059723e+01,0.6287008313071e+01],
             [0.3215657197281e-08,0.4866676894425e+01,0.5749861718712e+01],
             [0.4129362656266e-08,0.3809342558906e+01,0.5905702259363e+01],
             [0.3137762976388e-08,0.2494635174443e+01,0.2099539292909e+02],
             [0.3514010952384e-08,0.2699961831678e+01,0.7335344340001e+01],
             [0.3327607571530e-08,0.3318457714816e+01,0.5436992986000e+01],
             [0.3541066946675e-08,0.4382703582466e+01,0.1234573916645e+02],
             [0.3216179847052e-08,0.5271066317054e+01,0.3802769619140e-01],
             [0.2959045059570e-08,0.5819591585302e+01,0.2670964694522e+02],
             [0.3884040326665e-08,0.5980934960428e+01,0.6660449441528e+01],
             [0.2922027539886e-08,0.3337290282483e+01,0.1375773836557e+01],
             [0.4110846382042e-08,0.5742978187327e+01,0.4480965020977e+02],
             [0.2934508411032e-08,0.2278075804200e+01,0.6408777551755e+00],
             [0.3966896193000e-08,0.5835747858477e+01,0.3773735910827e+00],
             [0.3286695827610e-08,0.5838898193902e+01,0.3932462625300e-02],
             [0.3720643094196e-08,0.1122212337858e+01,0.1646033343740e+02],
             [0.3285508906174e-08,0.9182250996416e+00,0.1081813534213e+02],
             [0.3753880575973e-08,0.5174761973266e+01,0.5642198095270e+01],
             [0.3022129385587e-08,0.3381611020639e+01,0.2982630633589e+02],
             [0.2798569205621e-08,0.3546193723922e+01,0.1937891852345e+02],
             [0.3397872070505e-08,0.4533203197934e+01,0.6923953605621e+01],
             [0.3708099772977e-08,0.2756168198616e+01,0.3066615496545e+02],
             [0.3599283541510e-08,0.1934395469918e+01,0.6147450479709e+01],
             [0.3688702753059e-08,0.7149920971109e+00,0.2636725487657e+01],
             [0.2681084724003e-08,0.4899819493154e+01,0.6816289982179e+01],
             [0.3495993460759e-08,0.1572418915115e+01,0.6418701221183e+01],
             [0.3130770324995e-08,0.8912190180489e+00,0.1235996607578e+02],
             [0.2744353821941e-08,0.3800821940055e+01,0.2059724391010e+02],
             [0.2842732906341e-08,0.2644717440029e+01,0.2828699048865e+02],
             [0.3046882682154e-08,0.3987793020179e+01,0.6055599646783e+01],
             [0.2399072455143e-08,0.9908826440764e+00,0.6255674361143e+01],
             [0.2384306274204e-08,0.2516149752220e+01,0.6310477339748e+01],
             [0.2977324500559e-08,0.5849195642118e+01,0.1652265972112e+02],
             [0.3062835258972e-08,0.1681660100162e+01,0.1172006883645e+02],
             [0.3109682589231e-08,0.5804143987737e+00,0.2751146787858e+02],
             [0.2903920355299e-08,0.5800768280123e+01,0.6510552054109e+01],
             [0.2823221989212e-08,0.9241118370216e+00,0.5469525544182e+01],
             [0.3187949696649e-08,0.3139776445735e+01,0.1693792562116e+03],
             [0.2922559771655e-08,0.3549440782984e+01,0.2630839062450e+00],
             [0.2436302066603e-08,0.4735540696319e+01,0.3946258593675e+00],
             [0.3049473043606e-08,0.4998289124561e+01,0.8390110365991e+01],
             [0.2863682575784e-08,0.6709515671102e+00,0.2243449970715e+00],
             [0.2641750517966e-08,0.5410978257284e+01,0.2986433403208e+02],
             [0.2704093466243e-08,0.4778317207821e+01,0.6129297044991e+01],
             [0.2445522177011e-08,0.6009020662222e+01,0.1171295538178e+02],
             [0.2623608810230e-08,0.5010449777147e+01,0.6436854655901e+01],
             [0.2079259704053e-08,0.5980943768809e+01,0.2019909489111e+02],
             [0.2820225596771e-08,0.2679965110468e+01,0.5934151399930e+01],
             [0.2365221950927e-08,0.1894231148810e+01,0.2470570524223e+02],
             [0.2359682077149e-08,0.4220752950780e+01,0.8671969964381e+01],
             [0.2387577137206e-08,0.2571783940617e+01,0.7096626156709e+01],
             [0.1982102089816e-08,0.5169765997119e+00,0.1727188400790e+02],
             [0.2687502389925e-08,0.6239078264579e+01,0.7075506709219e+02],
             [0.2207751669135e-08,0.2031184412677e+01,0.4377611041777e+01],
             [0.2618370214274e-08,0.8266079985979e+00,0.6632000300961e+01],
             [0.2591951887361e-08,0.8819350522008e+00,0.4873985990671e+02],
             [0.2375055656248e-08,0.3520944177789e+01,0.1590676413561e+02],
             [0.2472019978911e-08,0.1551431908671e+01,0.6612329252343e+00],
             [0.2368157127199e-08,0.4178610147412e+01,0.3459636466239e+02],
             [0.1764846605693e-08,0.1506764000157e+01,0.1980094587212e+02],
             [0.2291769608798e-08,0.2118250611782e+01,0.2844914056730e-01],
             [0.2209997316943e-08,0.3363255261678e+01,0.2666070658668e+00],
             [0.2292699097923e-08,0.4200423956460e+00,0.1484170571900e-02],
             [0.1629683015329e-08,0.2331362582487e+01,0.3035599730800e+02],
             [0.2206492862426e-08,0.3400274026992e+01,0.6281667977667e+01],
             [0.2205746568257e-08,0.1066051230724e+00,0.6284483723224e+01],
             [0.2026310767991e-08,0.2779066487979e+01,0.2449240616245e+02],
             [0.1762977622163e-08,0.9951450691840e+00,0.2045286941806e+02],
             [0.1368535049606e-08,0.6402447365817e+00,0.2473415438279e+02],
             [0.1720598775450e-08,0.2303524214705e+00,0.1679593901136e+03],
             [0.1702429015449e-08,0.6164622655048e+01,0.3338575901272e+03],
             [0.1414033197685e-08,0.3954561185580e+01,0.1624205518357e+03],
             [0.1573768958043e-08,0.2028286308984e+01,0.3144167757552e+02],
             [0.1650705184447e-08,0.2304040666128e+01,0.5267006960365e+02],
             [0.1651087618855e-08,0.2538461057280e+01,0.8956999012000e+02],
             [0.1616409518983e-08,0.5111054348152e+01,0.3332657872986e+02],
             [0.1537175173581e-08,0.5601130666603e+01,0.3852657435933e+02],
             [0.1593191980553e-08,0.2614340453411e+01,0.2282781046519e+03],
             [0.1499480170643e-08,0.3624721577264e+01,0.2823723341956e+02],
             [0.1493807843235e-08,0.4214569879008e+01,0.2876692439167e+02],
             [0.1074571199328e-08,0.1496911744704e+00,0.8397383534231e+02],
             [0.1074406983417e-08,0.1187817671922e+01,0.8401985929482e+02],
             [0.9757576855851e-09,0.2655703035858e+01,0.7826370942180e+02],
             [0.1258432887565e-08,0.4969896184844e+01,0.3115650189215e+03],
             [0.1240336343282e-08,0.5192460776926e+01,0.1784300471910e+03],
             [0.9016107005164e-09,0.1960356923057e+01,0.5886454391678e+02],
             [0.1135392360918e-08,0.5082427809068e+01,0.7842370451713e+02],
             [0.9216046089565e-09,0.2793775037273e+01,0.1014262087719e+03],
             [0.1061276615030e-08,0.3726144311409e+01,0.5660027930059e+02],
             [0.1010110596263e-08,0.7404080708937e+00,0.4245678405627e+02],
             [0.7217424756199e-09,0.2697449980577e-01,0.2457074661053e+03],
             [0.6912003846756e-09,0.4253296276335e+01,0.1679936946371e+03],
             [0.6871814664847e-09,0.5148072412354e+01,0.6053048899753e+02],
             [0.4887158016343e-09,0.2153581148294e+01,0.9656299901946e+02],
             [0.5161802866314e-09,0.3852750634351e+01,0.2442876000072e+03],
             [0.5652599559057e-09,0.1233233356270e+01,0.8365903305582e+02],
             [0.4710812608586e-09,0.5610486976767e+01,0.3164282286739e+03],
             [0.4909977500324e-09,0.1639629524123e+01,0.4059982187939e+03],
             [0.4772641839378e-09,0.3737100368583e+01,0.1805255418145e+03],
             [0.4487562567153e-09,0.1158417054478e+00,0.8433466158131e+02],
             [0.3943441230497e-09,0.6243502862796e+00,0.2568537517081e+03],
             [0.3952236913598e-09,0.3510377382385e+01,0.2449975330562e+03],
             [0.3788898363417e-09,0.5916128302299e+01,0.1568131045107e+03],
             [0.3738329328831e-09,0.1042266763456e+01,0.3948519331910e+03],
             [0.2451199165151e-09,0.1166788435700e+01,0.1435713242844e+03],
             [0.2436734402904e-09,0.3254726114901e+01,0.2268582385539e+03],
             [0.2213605274325e-09,0.1687210598530e+01,0.1658638954901e+03],
             [0.1491521204829e-09,0.2657541786794e+01,0.2219950288015e+03],
             [0.1474995329744e-09,0.5013089805819e+01,0.3052819430710e+03],
             [0.1661939475656e-09,0.5495315428418e+01,0.2526661704812e+03],
             [0.9015946748003e-10,0.2236989966505e+01,0.4171445043968e+03]]
    
    #Sun-to-Earth, T^1, Y
    E1[1]=[[0.9304690546528e-06,0.0000000000000e+00,0.0000000000000e+00],
             [0.5150715570663e-06,0.4431807116294e+01,0.1256615170089e+02],
             [0.1290825411056e-07,0.4388610039678e+01,0.1884922755134e+02],
             [0.4645466665386e-08,0.5827263376034e+01,0.6283075850446e+01],
             [0.2079625310718e-08,0.1621698662282e+00,0.6279552690824e+01],
             [0.2078189850907e-08,0.3344713435140e+01,0.6286599010068e+01],
             [0.6207190138027e-09,0.5074049319576e+01,0.4705732307012e+01],
             [0.5989826532569e-09,0.2231842216620e+01,0.6256777527156e+01],
             [0.5961360812618e-09,0.1274975769045e+01,0.6309374173736e+01],
             [0.4874165471016e-09,0.3642277426779e+01,0.7755226100720e+00],
             [0.4283834034360e-09,0.5148765510106e+01,0.1059381944224e+01],
             [0.4652389287529e-09,0.4715794792175e+01,0.7860419393880e+01],
             [0.3751707476401e-09,0.6617207370325e+00,0.5753384878334e+01],
             [0.3559998806198e-09,0.6155548875404e+01,0.5884926831456e+01],
             [0.3558447558857e-09,0.2898827297664e+01,0.6812766822558e+01],
             [0.3211116927106e-09,0.3625813502509e+01,0.6681224869435e+01],
             [0.2875609914672e-09,0.4345435813134e+01,0.2513230340178e+02],
             [0.2843109704069e-09,0.5862263940038e+01,0.6127655567643e+01],
             [0.2744676468427e-09,0.3926419475089e+01,0.6438496133249e+01],
             [0.2481285237789e-09,0.1351976572828e+01,0.5486777812467e+01],
             [0.2060338481033e-09,0.2147556998591e+01,0.7079373888424e+01],
             [0.2015822358331e-09,0.4408358972216e+01,0.6290189305114e+01],
             [0.2001195944195e-09,0.5385829822531e+01,0.6275962395778e+01],
             [0.1953667642377e-09,0.1304933746120e+01,0.5507553240374e+01],
             [0.1839744078713e-09,0.6173567228835e+01,0.1179062909082e+02],
             [0.1643334294845e-09,0.4635942997523e+01,0.1150676975667e+02],
             [0.1768051018652e-09,0.5086283558874e+01,0.7113454667900e-02],
             [0.1674874205489e-09,0.2243332137241e+01,0.7058598460518e+01],
             [0.1421445397609e-09,0.6186899771515e+01,0.7962980379786e+00],
             [0.1255163958267e-09,0.5730238465658e+01,0.4694002934110e+01],
             [0.1013945281961e-09,0.1726055228402e+01,0.3738761453707e+01],
             [0.1047294335852e-09,0.2658801228129e+01,0.6282095334605e+01],
             [0.1047103879392e-09,0.8481047835035e+00,0.6284056366286e+01],
             [0.9530343962826e-10,0.3079267149859e+01,0.6069776770667e+01],
             [0.9604637611690e-10,0.3258679792918e+00,0.4136910472696e+01],
             [0.9153518537177e-10,0.4398599886584e+00,0.6496374930224e+01],
             [0.8562458214922e-10,0.4772686794145e+01,0.1194447056968e+01],
             [0.8232525360654e-10,0.5966220721679e+01,0.1589072916335e+01],
             [0.6150223411438e-10,0.1780985591923e+01,0.8827390247185e+01],
             [0.6272087858000e-10,0.3184305429012e+01,0.8429241228195e+01],
             [0.5540476311040e-10,0.3801260595433e+01,0.4933208510675e+01],
             [0.7331901699361e-10,0.5205948591865e+01,0.4535059491685e+01],
             [0.6018528702791e-10,0.4770139083623e+01,0.1255903824622e+02],
             [0.5150530724804e-10,0.3574796899585e+01,0.1176985366291e+02],
             [0.6471933741811e-10,0.2679787266521e+01,0.5088628793478e+01],
             [0.5317460644174e-10,0.9528763345494e+00,0.3154687086868e+01],
             [0.4832187748783e-10,0.5329322498232e+01,0.6040347114260e+01],
             [0.4716763555110e-10,0.2395235316466e+01,0.5331357529664e+01],
             [0.4871509139861e-10,0.3056663648823e+01,0.1256967486051e+02],
             [0.4598417696768e-10,0.4452762609019e+01,0.6525804586632e+01],
             [0.5674189533175e-10,0.9879680872193e+00,0.5729506548653e+01],
             [0.4073560328195e-10,0.5939127696986e+01,0.7632943190217e+01],
             [0.5040994945359e-10,0.4549875824510e+01,0.8031092209206e+01],
             [0.5078185134679e-10,0.7346659893982e+00,0.7477522907414e+01],
             [0.3769343537061e-10,0.1071317188367e+01,0.7234794171227e+01],
             [0.4980331365299e-10,0.2500345341784e+01,0.6836645152238e+01],
             [0.3458236594757e-10,0.3825159450711e+01,0.1097707878456e+02],
             [0.3578859493602e-10,0.5299664791549e+01,0.4164311961999e+01],
             [0.3370504646419e-10,0.5002316301593e+01,0.1137170464392e+02],
             [0.3299873338428e-10,0.2526123275282e+01,0.3930209696940e+01],
             [0.4304917318409e-10,0.3368078557132e+01,0.1592596075957e+01],
             [0.3402418753455e-10,0.8385495425800e+00,0.3128388763578e+01],
             [0.2778460572146e-10,0.3669905203240e+01,0.7342457794669e+01],
             [0.2782710128902e-10,0.2691664812170e+00,0.1748016358760e+01],
             [0.2711725179646e-10,0.4707487217718e+01,0.5296909721118e+00],
             [0.2981760946340e-10,0.3190260867816e+00,0.5368044267797e+00],
             [0.2811672977772e-10,0.3196532315372e+01,0.7084896783808e+01],
             [0.2863454474467e-10,0.2263240324780e+00,0.5223693906222e+01],
             [0.3333464634051e-10,0.3498451685065e+01,0.8018209333619e+00],
             [0.3312991747609e-10,0.5839154477412e+01,0.1554202828031e+00],
             [0.2813255564006e-10,0.8268044346621e+00,0.5225775174439e+00],
             [0.2665098083966e-10,0.3934021725360e+01,0.5216580451554e+01],
             [0.2349795705216e-10,0.5197620913779e+01,0.2146165377750e+01],
             [0.2330352293961e-10,0.2984999231807e+01,0.1726015463500e+02],
             [0.2728001683419e-10,0.6521679638544e+00,0.8635942003952e+01],
             [0.2484061007669e-10,0.3468955561097e+01,0.5230807360890e+01],
             [0.2646328768427e-10,0.1013724533516e+01,0.2629832328990e-01],
             [0.2518630264831e-10,0.6108081057122e+01,0.5481254917084e+01],
             [0.2421901455384e-10,0.1651097776260e+01,0.1349867339771e+01],
             [0.6348533267831e-11,0.3220226560321e+01,0.8433466158131e+02]]
    
    #Sun-to-Earth, T^2, Y
    E2[1]=[[0.5063375872532e-10,0.0000000000000e+00,0.0000000000000e+00],
             [0.2173815785980e-10,0.2827805833053e+01,0.1256615170089e+02],
             [0.1010231999920e-10,0.4634612377133e+01,0.6283075850446e+01],
             [0.9259745317636e-12,0.2620612076189e+01,0.1884922755134e+02],
             [0.1022202095812e-12,0.3809562326066e+01,0.8399684731857e+02]]    
    
    #Sun-to-Earth, T^0, Z
    E0[2]=[[0.2796207639075e-05,0.3198701560209e+01,0.8433466158131e+02],
             [0.1016042198142e-05,0.5422360395913e+01,0.5507553240374e+01],
             [0.8044305033647e-06,0.3880222866652e+01,0.5223693906222e+01],
             [0.4385347909274e-06,0.3704369937468e+01,0.2352866153506e+01],
             [0.3186156414906e-06,0.3999639363235e+01,0.1577343543434e+01],
             [0.2272412285792e-06,0.3984738315952e+01,0.1047747311755e+01],
             [0.1645620103007e-06,0.3565412516841e+01,0.5856477690889e+01],
             [0.1815836921166e-06,0.4984507059020e+01,0.6283075850446e+01],
             [0.1447461676364e-06,0.3702753570108e+01,0.9437762937313e+01],
             [0.1430760876382e-06,0.3409658712357e+01,0.1021328554739e+02],
             [0.1120445753226e-06,0.4829561570246e+01,0.1414349524433e+02],
             [0.1090232840797e-06,0.2080729178066e+01,0.6812766822558e+01],
             [0.9715727346551e-07,0.3476295881948e+01,0.4694002934110e+01],
             [0.1036267136217e-06,0.4056639536648e+01,0.7109288135493e+02],
             [0.8752665271340e-07,0.4448159519911e+01,0.5753384878334e+01],
             [0.8331864956004e-07,0.4991704044208e+01,0.7084896783808e+01],
             [0.6901658670245e-07,0.4325358994219e+01,0.6275962395778e+01],
             [0.9144536848998e-07,0.1141826375363e+01,0.6620890113188e+01],
             [0.7205085037435e-07,0.3624344170143e+01,0.5296909721118e+00],
             [0.7697874654176e-07,0.5554257458998e+01,0.1676215758509e+03],
             [0.5197545738384e-07,0.6251760961735e+01,0.1807370494127e+02],
             [0.5031345378608e-07,0.2497341091913e+01,0.4705732307012e+01],
             [0.4527110205840e-07,0.2335079920992e+01,0.6309374173736e+01],
             [0.4753355798089e-07,0.7094148987474e+00,0.5884926831456e+01],
             [0.4296951977516e-07,0.1101916352091e+01,0.6681224869435e+01],
             [0.3855341568387e-07,0.1825495405486e+01,0.5486777812467e+01],
             [0.5253930970990e-07,0.4424740687208e+01,0.7860419393880e+01],
             [0.4024630496471e-07,0.5120498157053e+01,0.1336797263425e+02],
             [0.4061069791453e-07,0.6029771435451e+01,0.3930209696940e+01],
             [0.3797883804205e-07,0.4435193600836e+00,0.3154687086868e+01],
             [0.2933033225587e-07,0.5124157356507e+01,0.1059381944224e+01],
             [0.3503000930426e-07,0.5421830162065e+01,0.6069776770667e+01],
             [0.3670096214050e-07,0.4582101667297e+01,0.1219403291462e+02],
             [0.2905609437008e-07,0.1926566420072e+01,0.1097707878456e+02],
             [0.2466827821713e-07,0.6090174539834e+00,0.6496374930224e+01],
             [0.2691647295332e-07,0.1393432595077e+01,0.2200391463820e+02],
             [0.2150554667946e-07,0.4308671715951e+01,0.5643178611111e+01],
             [0.2237481922680e-07,0.8133968269414e+00,0.8635942003952e+01],
             [0.1817741038157e-07,0.3755205127454e+01,0.3340612434717e+01],
             [0.2227820762132e-07,0.2759558596664e+01,0.1203646072878e+02],
             [0.1944713772307e-07,0.5699645869121e+01,0.1179062909082e+02],
             [0.1527340520662e-07,0.1986749091746e+01,0.3981490189893e+00],
             [0.1577282574914e-07,0.3205017217983e+01,0.5088628793478e+01],
             [0.1424738825424e-07,0.6256747903666e+01,0.2544314396739e+01],
             [0.1616563121701e-07,0.2601671259394e+00,0.1729818233119e+02],
             [0.1401210391692e-07,0.4686939173506e+01,0.7058598460518e+01],
             [0.1488726974214e-07,0.2815862451372e+01,0.2593412433514e+02],
             [0.1692626442388e-07,0.4956894109797e+01,0.1564752902480e+03],
             [0.1123571582910e-07,0.2381192697696e+01,0.3738761453707e+01],
             [0.9903308606317e-08,0.4294851657684e+01,0.9225539266174e+01],
             [0.9174533187191e-08,0.3075171510642e+01,0.4164311961999e+01],
             [0.8645985631457e-08,0.5477534821633e+00,0.8429241228195e+01],
             [-0.1085876492688e-07,0.0000000000000e+00,0.0000000000000e+00],
             [0.9264309077815e-08,0.5968571670097e+01,0.7079373888424e+01],
             [0.8243116984954e-08,0.1489098777643e+01,0.1044738781244e+02],
             [0.8268102113708e-08,0.3512977691983e+01,0.1150676975667e+02],
             [0.9043613988227e-08,0.1290704408221e+00,0.1101510648075e+02],
             [0.7432912038789e-08,0.1991086893337e+01,0.2608790314060e+02],
             [0.8586233727285e-08,0.4238357924414e+01,0.2986433403208e+02],
             [0.7612230060131e-08,0.2911090150166e+01,0.4732030630302e+01],
             [0.7097787751408e-08,0.1908938392390e+01,0.8031092209206e+01],
             [0.7640237040175e-08,0.6129219000168e+00,0.7962980379786e+00],
             [0.7070445688081e-08,0.1380417036651e+01,0.2146165377750e+01],
             [0.7690770957702e-08,0.1680504249084e+01,0.2122839202813e+02],
             [0.8051292542594e-08,0.5127423484511e+01,0.2942463415728e+01],
             [0.5902709104515e-08,0.2020274190917e+01,0.7755226100720e+00],
             [0.5134567496462e-08,0.2606778676418e+01,0.1256615170089e+02],
             [0.5525802046102e-08,0.1613011769663e+01,0.8018209333619e+00],
             [0.5880724784221e-08,0.4604483417236e+01,0.4690479774488e+01],
             [0.5211699081370e-08,0.5718964114193e+01,0.8827390247185e+01],
             [0.4891849573562e-08,0.3689658932196e+01,0.2132990797783e+00],
             [0.5150246069997e-08,0.4099769855122e+01,0.6480980550449e+02],
             [0.5102434319633e-08,0.5660834602509e+01,0.3379454372902e+02],
             [0.5083405254252e-08,0.9842221218974e+00,0.4136910472696e+01],
             [0.4206562585682e-08,0.1341363634163e+00,0.3128388763578e+01],
             [0.4663249683579e-08,0.8130132735866e+00,0.5216580451554e+01],
             [0.4099474416530e-08,0.5791497770644e+01,0.4265981595566e+00],
             [0.4628251220767e-08,0.1249802769331e+01,0.1572083878776e+02],
             [0.5024068728142e-08,0.4795684802743e+01,0.6290189305114e+01],
             [0.5120234327758e-08,0.3810420387208e+01,0.5230807360890e+01],
             [0.5524029815280e-08,0.1029264714351e+01,0.2397622045175e+03],
             [0.4757415718860e-08,0.3528044781779e+01,0.1649636139783e+02],
             [0.3915786131127e-08,0.5593889282646e+01,0.1589072916335e+01],
             [0.4869053149991e-08,0.3299636454433e+01,0.7632943190217e+01],
             [0.3649365703729e-08,0.1286049002584e+01,0.6206810014183e+01],
             [0.3992493949002e-08,0.3100307589464e+01,0.2515860172507e+02],
             [0.3320247477418e-08,0.6212683940807e+01,0.1216800268190e+02],
             [0.3287123739696e-08,0.4699118445928e+01,0.7234794171227e+01],
             [0.3472776811103e-08,0.2630507142004e+01,0.7342457794669e+01],
             [0.3423253294767e-08,0.2946432844305e+01,0.9623688285163e+01],
             [0.3896173898244e-08,0.1224834179264e+01,0.6438496133249e+01],
             [0.3388455337924e-08,0.1543807616351e+01,0.1494531617769e+02],
             [0.3062704716523e-08,0.1191777572310e+01,0.8662240327241e+01],
             [0.3270075600400e-08,0.5483498767737e+01,0.1194447056968e+01],
             [0.3101209215259e-08,0.8000833804348e+00,0.3772475342596e+02],
             [0.2780883347311e-08,0.4077980721888e+00,0.5863591145557e+01],
             [0.2903605931824e-08,0.2617490302147e+01,0.1965104848470e+02],
             [0.2682014743119e-08,0.2634703158290e+01,0.7238675589263e+01],
             [0.2534360108492e-08,0.6102446114873e+01,0.6836645152238e+01],
             [0.2392564882509e-08,0.3681820208691e+01,0.5849364236221e+01],
             [0.2656667254856e-08,0.6216045388886e+01,0.6133512519065e+01],
             [0.2331242096773e-08,0.5864949777744e+01,0.4535059491685e+01],
             [0.2287898363668e-08,0.4566628532802e+01,0.7477522907414e+01],
             [0.2336944521306e-08,0.2442722126930e+01,0.1137170464392e+02],
             [0.3156632236269e-08,0.1626628050682e+01,0.2509084901204e+03],
             [0.2982612402766e-08,0.2803604512609e+01,0.1748016358760e+01],
             [0.2774031674807e-08,0.4654002897158e+01,0.8223916695780e+02],
             [0.2295236548638e-08,0.4326518333253e+01,0.3378142627421e+00],
             [0.2190714699873e-08,0.4519614578328e+01,0.2908881142201e+02],
             [0.2191495845045e-08,0.3012626912549e+01,0.1673046366289e+02],
             [0.2492901628386e-08,0.1290101424052e+00,0.1543797956245e+03],
             [0.1993778064319e-08,0.3864046799414e+01,0.1778984560711e+02],
             [0.1898146479022e-08,0.5053777235891e+01,0.2042657109477e+02],
             [0.1918280127634e-08,0.2222470192548e+01,0.4165496312290e+02],
             [0.1916351061607e-08,0.8719067257774e+00,0.7737595720538e+02],
             [0.1834720181466e-08,0.4031491098040e+01,0.2358125818164e+02],
             [0.1249201523806e-08,0.5938379466835e+01,0.3301902111895e+02],
             [0.1477304050539e-08,0.6544722606797e+00,0.9548094718417e+02],
             [0.1264316431249e-08,0.2059072853236e+01,0.8399684731857e+02],
             [0.1203526495039e-08,0.3644813532605e+01,0.4558517281984e+02],
             [0.9221681059831e-09,0.3241815055602e+01,0.7805158573086e+02],
             [0.7849278367646e-09,0.5043812342457e+01,0.5217580628120e+02],
             [0.7983392077387e-09,0.5000024502753e+01,0.1501922143975e+03],
             [0.7925395431654e-09,0.1398734871821e-01,0.9061773743175e+02],
             [0.7640473285886e-09,0.5067111723130e+01,0.4951538251678e+02],
             [0.5398937754482e-09,0.5597382200075e+01,0.1613385000004e+03],
             [0.5626247550193e-09,0.2601338209422e+01,0.7318837597844e+02],
             [0.5525197197855e-09,0.5814832109256e+01,0.1432335100216e+03],
             [0.5407629837898e-09,0.3384820609076e+01,0.3230491187871e+03],
             [0.3856739119801e-09,0.1072391840473e+01,0.2334791286671e+03],
             [0.3856425239987e-09,0.2369540393327e+01,0.1739046517013e+03],
             [0.4350867755983e-09,0.5255575751082e+01,0.1620484330494e+03],
             [0.3844113924996e-09,0.5482356246182e+01,0.9757644180768e+02],
             [0.2854869155431e-09,0.9573634763143e+00,0.1697170704744e+03],
             [0.1719227671416e-09,0.1887203025202e+01,0.2265204242912e+03],
             [0.1527846879755e-09,0.3982183931157e+01,0.3341954043900e+03],
             [0.1128229264847e-09,0.2787457156298e+01,0.3119028331842e+03]]
    
    #Sun-to-Earth, T^1, Z
    E1[2]=[[0.2278290449966e-05,0.3413716033863e+01,0.6283075850446e+01],
             [0.5429458209830e-07,0.0000000000000e+00,0.0000000000000e+00],
             [0.1903240492525e-07,0.3370592358297e+01,0.1256615170089e+02],
             [0.2385409276743e-09,0.3327914718416e+01,0.1884922755134e+02],
             [0.8676928342573e-10,0.1824006811264e+01,0.5223693906222e+01],
             [0.7765442593544e-10,0.3888564279247e+01,0.5507553240374e+01],
             [0.7066158332715e-10,0.5194267231944e+01,0.2352866153506e+01],
             [0.7092175288657e-10,0.2333246960021e+01,0.8399684731857e+02],
             [0.5357582213535e-10,0.2224031176619e+01,0.5296909721118e+00],
             [0.3828035865021e-10,0.2156710933584e+01,0.6279552690824e+01],
             [0.3824857220427e-10,0.1529755219915e+01,0.6286599010068e+01],
             [0.3286995181628e-10,0.4879512900483e+01,0.1021328554739e+02]]
    
    #Sun-to-Earth, T^2, Z
    E2[2]=[[0.9722666114891e-10,0.5152219582658e+01,0.6283075850446e+01],
             [-0.3494819171909e-11,0.0000000000000e+00,0.0000000000000e+00],
             [0.6713034376076e-12,0.6440188750495e+00,0.1256615170089e+02]]

    #SSB-to-Sun, T^0, X
    S0[0]=[[0.4956757536410e-02,0.3741073751789e+01,0.5296909721118e+00],
             [0.2718490072522e-02,0.4016011511425e+01,0.2132990797783e+00],
             [0.1546493974344e-02,0.2170528330642e+01,0.3813291813120e-01],
             [0.8366855276341e-03,0.2339614075294e+01,0.7478166569050e-01],
             [0.2936777942117e-03,0.0000000000000e+00,0.0000000000000e+00],
             [0.1201317439469e-03,0.4090736353305e+01,0.1059381944224e+01],
             [0.7578550887230e-04,0.3241518088140e+01,0.4265981595566e+00],
             [0.1941787367773e-04,0.1012202064330e+01,0.2061856251104e+00],
             [0.1889227765991e-04,0.3892520416440e+01,0.2204125344462e+00],
             [0.1937896968613e-04,0.4797779441161e+01,0.1495633313810e+00],
             [0.1434506110873e-04,0.3868960697933e+01,0.5225775174439e+00],
             [0.1406659911580e-04,0.4759766557397e+00,0.5368044267797e+00],
             [0.1179022300202e-04,0.7774961520598e+00,0.7626583626240e-01],
             [0.8085864460959e-05,0.3254654471465e+01,0.3664874755930e-01],
             [0.7622752967615e-05,0.4227633103489e+01,0.3961708870310e-01],
             [0.6209171139066e-05,0.2791828325711e+00,0.7329749511860e-01],
             [0.4366435633970e-05,0.4440454875925e+01,0.1589072916335e+01],
             [0.3792124889348e-05,0.5156393842356e+01,0.7113454667900e-02],
             [0.3154548963402e-05,0.6157005730093e+01,0.4194847048887e+00],
             [0.3088359882942e-05,0.2494567553163e+01,0.6398972393349e+00],
             [0.2788440902136e-05,0.4934318747989e+01,0.1102062672231e+00],
             [0.3039928456376e-05,0.4895077702640e+01,0.6283075850446e+01],
             [0.2272258457679e-05,0.5278394064764e+01,0.1030928125552e+00],
             [0.2162007057957e-05,0.5802978019099e+01,0.3163918923335e+00],
             [0.1767632855737e-05,0.3415346595193e-01,0.1021328554739e+02],
             [0.1349413459362e-05,0.2001643230755e+01,0.1484170571900e-02],
             [0.1170141900476e-05,0.2424750491620e+01,0.6327837846670e+00],
             [0.1054355266820e-05,0.3123311487576e+01,0.4337116142245e+00],
             [0.9800822461610e-06,0.3026258088130e+01,0.1052268489556e+01],
             [0.1091203749931e-05,0.3157811670347e+01,0.1162474756779e+01],
             [0.6960236715913e-06,0.8219570542313e+00,0.1066495398892e+01],
             [0.5689257296909e-06,0.1323052375236e+01,0.9491756770005e+00],
             [0.6613172135802e-06,0.2765348881598e+00,0.8460828644453e+00],
             [0.6277702517571e-06,0.5794064466382e+01,0.1480791608091e+00],
             [0.6304884066699e-06,0.7323555380787e+00,0.2243449970715e+00],
             [0.4897850467382e-06,0.3062464235399e+01,0.3340612434717e+01],
             [0.3759148598786e-06,0.4588290469664e+01,0.3516457698740e-01],
             [0.3110520548195e-06,0.1374299536572e+01,0.6373574839730e-01],
             [0.3064708359780e-06,0.4222267485047e+01,0.1104591729320e-01],
             [0.2856347168241e-06,0.3714202944973e+01,0.1510475019529e+00],
             [0.2840945514288e-06,0.2847972875882e+01,0.4110125927500e-01],
             [0.2378951599405e-06,0.3762072563388e+01,0.2275259891141e+00],
             [0.2714229481417e-06,0.1036049980031e+01,0.2535050500000e-01],
             [0.2323551717307e-06,0.4682388599076e+00,0.8582758298370e-01],
             [0.1881790512219e-06,0.4790565425418e+01,0.2118763888447e+01],
             [0.2261353968371e-06,0.1669144912212e+01,0.7181332454670e-01],
             [0.2214546389848e-06,0.3937717281614e+01,0.2968341143800e-02],
             [0.2184915594933e-06,0.1129169845099e+00,0.7775000683430e-01],
             [0.2000164937936e-06,0.4030009638488e+01,0.2093666171530e+00],
             [0.1966105136719e-06,0.8745955786834e+00,0.2172315424036e+00],
             [0.1904742332624e-06,0.5919743598964e+01,0.2022531624851e+00],
             [0.1657399705031e-06,0.2549141484884e+01,0.7358765972222e+00],
             [0.1574070533987e-06,0.5277533020230e+01,0.7429900518901e+00],
             [0.1832261651039e-06,0.3064688127777e+01,0.3235053470014e+00],
             [0.1733615346569e-06,0.3011432799094e+01,0.1385174140878e+00],
             [0.1549124014496e-06,0.4005569132359e+01,0.5154640627760e+00],
             [0.1637044713838e-06,0.1831375966632e+01,0.8531963191132e+00],
             [0.1123420082383e-06,0.1180270407578e+01,0.1990721704425e+00],
             [0.1083754165740e-06,0.3414101320863e+00,0.5439178814476e+00],
             [0.1156638012655e-06,0.6130479452594e+00,0.5257585094865e+00],
             [0.1142548785134e-06,0.3724761948846e+01,0.5336234347371e+00],
             [0.7921463895965e-07,0.2435425589361e+01,0.1478866649112e+01],
             [0.7428600285231e-07,0.3542144398753e+01,0.2164800718209e+00],
             [0.8323211246747e-07,0.3525058072354e+01,0.1692165728891e+01],
             [0.7257595116312e-07,0.1364299431982e+01,0.2101180877357e+00],
             [0.7111185833236e-07,0.2460478875808e+01,0.4155522422634e+00],
             [0.6868090383716e-07,0.4397327670704e+01,0.1173197218910e+00],
             [0.7226419974175e-07,0.4042647308905e+01,0.1265567569334e+01],
             [0.6955642383177e-07,0.2865047906085e+01,0.9562891316684e+00],
             [0.7492139296331e-07,0.5014278994215e+01,0.1422690933580e-01],
             [0.6598363128857e-07,0.2376730020492e+01,0.6470106940028e+00],
             [0.7381147293385e-07,0.3272990384244e+01,0.1581959461667e+01],
             [0.6402909624032e-07,0.5302290955138e+01,0.9597935788730e-01],
             [0.6237454263857e-07,0.5444144425332e+01,0.7084920306520e-01],
             [0.5241198544016e-07,0.4215359579205e+01,0.5265099800692e+00],
             [0.5144463853918e-07,0.1218916689916e+00,0.5328719641544e+00],
             [0.5868164772299e-07,0.2369402002213e+01,0.7871412831580e-01],
             [0.6233195669151e-07,0.1254922242403e+01,0.2608790314060e+02],
             [0.6068463791422e-07,0.5679713760431e+01,0.1114304132498e+00],
             [0.4359361135065e-07,0.6097219641646e+00,0.1375773836557e+01],
             [0.4686510366826e-07,0.4786231041431e+01,0.1143987543936e+00],
             [0.3758977287225e-07,0.1167368068139e+01,0.1596186371003e+01],
             [0.4282051974778e-07,0.1519471064319e+01,0.2770348281756e+00],
             [0.5153765386113e-07,0.1860532322984e+01,0.2228608264996e+00],
             [0.4575129387188e-07,0.7632857887158e+00,0.1465949902372e+00],
             [0.3326844933286e-07,0.1298219485285e+01,0.5070101000000e-01],
             [0.3748617450984e-07,0.1046510321062e+01,0.4903339079539e+00],
             [0.2816756661499e-07,0.3434522346190e+01,0.2991266627620e+00],
             [0.3412750405039e-07,0.2523766270318e+01,0.3518164938661e+00],
             [0.2655796761776e-07,0.2904422260194e+01,0.6256703299991e+00],
             [0.2963597929458e-07,0.5923900431149e+00,0.1099462426779e+00],
             [0.2539523734781e-07,0.4851947722567e+01,0.1256615170089e+02],
             [0.2283087914139e-07,0.3400498595496e+01,0.6681224869435e+01],
             [0.2321309799331e-07,0.5789099148673e+01,0.3368040641550e-01],
             [0.2549657649750e-07,0.3991856479792e-01,0.1169588211447e+01],
             [0.2290462303977e-07,0.2788567577052e+01,0.1045155034888e+01],
             [0.1945398522914e-07,0.3290896998176e+01,0.1155361302111e+01],
             [0.1849171512638e-07,0.2698060129367e+01,0.4452511715700e-02],
             [0.1647199834254e-07,0.3016735644085e+01,0.4408250688924e+00],
             [0.1529530765273e-07,0.5573043116178e+01,0.6521991896920e-01],
             [0.1433199339978e-07,0.1481192356147e+01,0.9420622223326e+00],
             [0.1729134193602e-07,0.1422817538933e+01,0.2108507877249e+00],
             [0.1716463931346e-07,0.3469468901855e+01,0.2157473718317e+00],
             [0.1391206061378e-07,0.6122436220547e+01,0.4123712502208e+00],
             [0.1404746661924e-07,0.1647765641936e+01,0.4258542984690e-01],
             [0.1410452399455e-07,0.5989729161964e+01,0.2258291676434e+00],
             [0.1089828772168e-07,0.2833705509371e+01,0.4226656969313e+00],
             [0.1047374564948e-07,0.5090690007331e+00,0.3092784376656e+00],
             [0.1358279126532e-07,0.5128990262836e+01,0.7923417740620e-01],
             [0.1020456476148e-07,0.9632772880808e+00,0.1456308687557e+00],
             [0.1033428735328e-07,0.3223779318418e+01,0.1795258541446e+01],
             [0.1412435841540e-07,0.2410271572721e+01,0.1525316725248e+00],
             [0.9722759371574e-08,0.2333531395690e+01,0.8434341241180e-01],
             [0.9657334084704e-08,0.6199270974168e+01,0.1272681024002e+01],
             [0.1083641148690e-07,0.2864222292929e+01,0.7032915397480e-01],
             [0.1067318403838e-07,0.5833458866568e+00,0.2123349582968e+00],
             [0.1062366201976e-07,0.4307753989494e+01,0.2142632012598e+00],
             [0.1236364149266e-07,0.2873917870593e+01,0.1847279083684e+00],
             [0.1092759489593e-07,0.2959887266733e+01,0.1370332435159e+00],
             [0.8912069362899e-08,0.5141213702562e+01,0.2648454860559e+01],
             [0.9656467707970e-08,0.4532182462323e+01,0.4376440768498e+00],
             [0.8098386150135e-08,0.2268906338379e+01,0.2880807454688e+00],
             [0.7857714675000e-08,0.4055544260745e+01,0.2037373330570e+00],
             [0.7288455940646e-08,0.5357901655142e+01,0.1129145838217e+00],
             [0.9450595950552e-08,0.4264926963939e+01,0.5272426800584e+00],
             [0.9381718247537e-08,0.7489366976576e-01,0.5321392641652e+00],
             [0.7079052646038e-08,0.1923311052874e+01,0.6288513220417e+00],
             [0.9259004415344e-08,0.2970256853438e+01,0.1606092486742e+00],
             [0.8259801499742e-08,0.3327056314697e+01,0.8389694097774e+00],
             [0.6476334355779e-08,0.2954925505727e+01,0.2008557621224e+01],
             [0.5984021492007e-08,0.9138753105829e+00,0.2042657109477e+02],
             [0.5989546863181e-08,0.3244464082031e+01,0.2111650433779e+01],
             [0.6233108606023e-08,0.4995232638403e+00,0.4305306221819e+00],
             [0.6877299149965e-08,0.2834987233449e+01,0.9561746721300e-02],
             [0.8311234227190e-08,0.2202951835758e+01,0.3801276407308e+00],
             [0.6599472832414e-08,0.4478581462618e+01,0.1063314406849e+01],
             [0.6160491096549e-08,0.5145858696411e+01,0.1368660381889e+01],
             [0.6164772043891e-08,0.3762976697911e+00,0.4234171675140e+00],
             [0.6363248684450e-08,0.3162246718685e+01,0.1253008786510e-01],
             [0.6448587520999e-08,0.3442693302119e+01,0.5287268506303e+00],
             [0.6431662283977e-08,0.8977549136606e+00,0.5306550935933e+00],
             [0.6351223158474e-08,0.4306447410369e+01,0.5217580628120e+02],
             [0.5476721393451e-08,0.3888529177855e+01,0.2221856701002e+01],
             [0.5341772572619e-08,0.2655560662512e+01,0.7466759693650e-01],
             [0.5337055758302e-08,0.5164990735946e+01,0.7489573444450e-01],
             [0.5373120816787e-08,0.6041214553456e+01,0.1274714967946e+00],
             [0.5392351705426e-08,0.9177763485932e+00,0.1055449481598e+01],
             [0.6688495850205e-08,0.3089608126937e+01,0.2213766559277e+00],
             [0.5072003660362e-08,0.4311316541553e+01,0.2132517061319e+00],
             [0.5070726650455e-08,0.5790675464444e+00,0.2133464534247e+00],
             [0.5658012950032e-08,0.2703945510675e+01,0.7287631425543e+00],
             [0.4835509924854e-08,0.2975422976065e+01,0.7160067364790e-01],
             [0.6479821978012e-08,0.1324168733114e+01,0.2209183458640e-01],
             [0.6230636494980e-08,0.2860103632836e+01,0.3306188016693e+00],
             [0.4649239516213e-08,0.4832259763403e+01,0.7796265773310e-01],
             [0.6487325792700e-08,0.2726165825042e+01,0.3884652414254e+00],
             [0.4682823682770e-08,0.6966602455408e+00,0.1073608853559e+01],
             [0.5704230804976e-08,0.5669634104606e+01,0.8731175355560e-01],
             [0.6125413585489e-08,0.1513386538915e+01,0.7605151500000e-01],
             [0.6035825038187e-08,0.1983509168227e+01,0.9846002785331e+00],
             [0.4331123462303e-08,0.2782892992807e+01,0.4297791515992e+00],
             [0.4681107685143e-08,0.5337232886836e+01,0.2127790306879e+00],
             [0.4669105829655e-08,0.5837133792160e+01,0.2138191288687e+00],
             [0.5138823602365e-08,0.3080560200507e+01,0.7233337363710e-01],
             [0.4615856664534e-08,0.1661747897471e+01,0.8603097737811e+00],
             [0.4496916702197e-08,0.2112508027068e+01,0.7381754420900e-01],
             [0.4278479042945e-08,0.5716528462627e+01,0.7574578717200e-01],
             [0.3840525503932e-08,0.6424172726492e+00,0.3407705765729e+00],
             [0.4866636509685e-08,0.4919244697715e+01,0.7722995774390e-01],
             [0.3526100639296e-08,0.2550821052734e+01,0.6225157782540e-01],
             [0.3939558488075e-08,0.3939331491710e+01,0.5268983110410e-01],
             [0.4041268772576e-08,0.2275337571218e+01,0.3503323232942e+00],
             [0.3948761842853e-08,0.1999324200790e+01,0.1451108196653e+00],
             [0.3258394550029e-08,0.9121001378200e+00,0.5296435984654e+00],
             [0.3257897048761e-08,0.3428428660869e+01,0.5297383457582e+00],
             [0.3842559031298e-08,0.6132927720035e+01,0.9098186128426e+00],
             [0.3109920095448e-08,0.7693650193003e+00,0.3932462625300e-02],
             [0.3132237775119e-08,0.3621293854908e+01,0.2346394437820e+00],
             [0.3942189421510e-08,0.4841863659733e+01,0.3180992042600e-02],
             [0.3796972285340e-08,0.1814174994268e+01,0.1862120789403e+00],
             [0.3995640233688e-08,0.1386990406091e+01,0.4549093064213e+00],
             [0.2875013727414e-08,0.9178318587177e+00,0.1905464808669e+01],
             [0.3073719932844e-08,0.2688923811835e+01,0.3628624111593e+00],
             [0.2731016580075e-08,0.1188259127584e+01,0.2131850110243e+00],
             [0.2729549896546e-08,0.3702160634273e+01,0.2134131485323e+00],
             [0.3339372892449e-08,0.7199163960331e+00,0.2007689919132e+00],
             [0.2898833764204e-08,0.1916709364999e+01,0.5291709230214e+00],
             [0.2894536549362e-08,0.2424043195547e+01,0.5302110212022e+00],
             [0.3096872473843e-08,0.4445894977497e+01,0.2976424921901e+00],
             [0.2635672326810e-08,0.3814366984117e+01,0.1485980103780e+01],
             [0.3649302697001e-08,0.2924200596084e+01,0.6044726378023e+00],
             [0.3127954585895e-08,0.1842251648327e+01,0.1084620721060e+00],
             [0.2616040173947e-08,0.4155841921984e+01,0.1258454114666e+01],
             [0.2597395859860e-08,0.1158045978874e+00,0.2103781122809e+00],
             [0.2593286172210e-08,0.4771850408691e+01,0.2162200472757e+00],
             [0.2481823585747e-08,0.4608842558889e+00,0.1062562936266e+01],
             [0.2742219550725e-08,0.1538781127028e+01,0.5651155736444e+00],
             [0.3199558469610e-08,0.3226647822878e+00,0.7036329877322e+00],
             [0.2666088542957e-08,0.1967991731219e+00,0.1400015846597e+00],
             [0.2397067430580e-08,0.3707036669873e+01,0.2125476091956e+00],
             [0.2376570772738e-08,0.1182086628042e+01,0.2140505503610e+00],
             [0.2547228007887e-08,0.4906256820629e+01,0.1534957940063e+00],
             [0.2265575594114e-08,0.3414949866857e+01,0.2235935264888e+00],
             [0.2464381430585e-08,0.4599122275378e+01,0.2091065926078e+00],
             [0.2433408527044e-08,0.2830751145445e+00,0.2174915669488e+00],
             [0.2443605509076e-08,0.4212046432538e+01,0.1739420156204e+00],
             [0.2319779262465e-08,0.9881978408630e+00,0.7530171478090e-01],
             [0.2284622835465e-08,0.5565347331588e+00,0.7426161660010e-01],
             [0.2467268750783e-08,0.5655708150766e+00,0.2526561439362e+00],
             [0.2808513492782e-08,0.1418405053408e+01,0.5636314030725e+00],
             [0.2329528932532e-08,0.4069557545675e+01,0.1056200952181e+01],
             [0.9698639532817e-09,0.1074134313634e+01,0.7826370942180e+02]]

    #SSB-to-Sun, T^1, X        
    S1[0]=[[-0.1296310361520e-07,0.0000000000000e+00,0.0000000000000e+00],
             [0.8975769009438e-08,0.1128891609250e+01,0.4265981595566e+00],
             [0.7771113441307e-08,0.2706039877077e+01,0.2061856251104e+00],
             [0.7538303866642e-08,0.2191281289498e+01,0.2204125344462e+00],
             [0.6061384579336e-08,0.3248167319958e+01,0.1059381944224e+01],
             [0.5726994235594e-08,0.5569981398610e+01,0.5225775174439e+00],
             [0.5616492836424e-08,0.5057386614909e+01,0.5368044267797e+00],
             [0.1010881584769e-08,0.3473577116095e+01,0.7113454667900e-02],
             [0.7259606157626e-09,0.3651858593665e+00,0.6398972393349e+00],
             [0.8755095026935e-09,0.1662835408338e+01,0.4194847048887e+00],
             [0.5370491182812e-09,0.1327673878077e+01,0.4337116142245e+00],
             [0.5743773887665e-09,0.4250200846687e+01,0.2132990797783e+00],
             [0.4408103140300e-09,0.3598752574277e+01,0.1589072916335e+01],
             [0.3101892374445e-09,0.4887822983319e+01,0.1052268489556e+01],
             [0.3209453713578e-09,0.9702272295114e+00,0.5296909721118e+00],
             [0.3017228286064e-09,0.5484462275949e+01,0.1066495398892e+01],
             [0.3200700038601e-09,0.2846613338643e+01,0.1495633313810e+00],
             [0.2137637279911e-09,0.5692163292729e+00,0.3163918923335e+00],
             [0.1899686386727e-09,0.2061077157189e+01,0.2275259891141e+00],
             [0.1401994545308e-09,0.4177771136967e+01,0.1102062672231e+00],
             [0.1578057810499e-09,0.5782460597335e+01,0.7626583626240e-01],
             [0.1237713253351e-09,0.5705900866881e+01,0.5154640627760e+00],
             [0.1313076837395e-09,0.5163438179576e+01,0.3664874755930e-01],
             [0.1184963304860e-09,0.3054804427242e+01,0.6327837846670e+00],
             [0.1238130878565e-09,0.2317292575962e+01,0.3961708870310e-01],
             [0.1015959527736e-09,0.2194643645526e+01,0.7329749511860e-01],
             [0.9017954423714e-10,0.2868603545435e+01,0.1990721704425e+00],
             [0.8668024955603e-10,0.4923849675082e+01,0.5439178814476e+00],
             [0.7756083930103e-10,0.3014334135200e+01,0.9491756770005e+00],
             [0.7536503401741e-10,0.2704886279769e+01,0.1030928125552e+00],
             [0.5483308679332e-10,0.6010983673799e+01,0.8531963191132e+00],
             [0.5184339620428e-10,0.1952704573291e+01,0.2093666171530e+00],
             [0.5108658712030e-10,0.2958575786649e+01,0.2172315424036e+00],
             [0.5019424524650e-10,0.1736317621318e+01,0.2164800718209e+00],
             [0.4909312625978e-10,0.3167216416257e+01,0.2101180877357e+00],
             [0.4456638901107e-10,0.7697579923471e+00,0.3235053470014e+00],
             [0.4227030350925e-10,0.3490910137928e+01,0.6373574839730e-01],
             [0.4095456040093e-10,0.5178888984491e+00,0.6470106940028e+00],
             [0.4990537041422e-10,0.3323887668974e+01,0.1422690933580e-01],
             [0.4321170010845e-10,0.4288484987118e+01,0.7358765972222e+00],
             [0.3544072091802e-10,0.6021051579251e+01,0.5265099800692e+00],
             [0.3480198638687e-10,0.4600027054714e+01,0.5328719641544e+00],
             [0.3440287244435e-10,0.4349525970742e+01,0.8582758298370e-01],
             [0.3330628322713e-10,0.2347391505082e+01,0.1104591729320e-01],
             [0.2973060707184e-10,0.4789409286400e+01,0.5257585094865e+00],
             [0.2932606766089e-10,0.5831693799927e+01,0.5336234347371e+00],
             [0.2876972310953e-10,0.2692638514771e+01,0.1173197218910e+00],
             [0.2827488278556e-10,0.2056052487960e+01,0.2022531624851e+00],
             [0.2515028239756e-10,0.7411863262449e+00,0.9597935788730e-01],
             [0.2853033744415e-10,0.3948481024894e+01,0.2118763888447e+01]]
    
    #SSB-to-Sun, T^2, X
    S2[0]=[[0.1603551636587e-11,0.4404109410481e+01,0.2061856251104e+00],
             [0.1556935889384e-11,0.4818040873603e+00,0.2204125344462e+00],
             [0.1182594414915e-11,0.9935762734472e+00,0.5225775174439e+00],
             [0.1158794583180e-11,0.3353180966450e+01,0.5368044267797e+00],
             [0.9597358943932e-12,0.5567045358298e+01,0.2132990797783e+00],
             [0.6511516579605e-12,0.5630872420788e+01,0.4265981595566e+00],
             [0.7419792747688e-12,0.2156188581957e+01,0.5296909721118e+00],
             [0.3951972655848e-12,0.1981022541805e+01,0.1059381944224e+01],
             [0.4478223877045e-12,0.0000000000000e+00,0.0000000000000e+00]]

    #SSB-to-Sun, T^0, Y
    S0[1]=[[0.4955392320126e-02,0.2170467313679e+01,0.5296909721118e+00],
             [0.2722325167392e-02,0.2444433682196e+01,0.2132990797783e+00],
             [0.1546579925346e-02,0.5992779281546e+00,0.3813291813120e-01],
             [0.8363140252966e-03,0.7687356310801e+00,0.7478166569050e-01],
             [0.3385792683603e-03,0.0000000000000e+00,0.0000000000000e+00],
             [0.1201192221613e-03,0.2520035601514e+01,0.1059381944224e+01],
             [0.7587125720554e-04,0.1669954006449e+01,0.4265981595566e+00],
             [0.1964155361250e-04,0.5707743963343e+01,0.2061856251104e+00],
             [0.1891900364909e-04,0.2320960679937e+01,0.2204125344462e+00],
             [0.1937373433356e-04,0.3226940689555e+01,0.1495633313810e+00],
             [0.1437139941351e-04,0.2301626908096e+01,0.5225775174439e+00],
             [0.1406267683099e-04,0.5188579265542e+01,0.5368044267797e+00],
             [0.1178703080346e-04,0.5489483248476e+01,0.7626583626240e-01],
             [0.8079835186041e-05,0.1683751835264e+01,0.3664874755930e-01],
             [0.7623253594652e-05,0.2656400462961e+01,0.3961708870310e-01],
             [0.6248667483971e-05,0.4992775362055e+01,0.7329749511860e-01],
             [0.4366353695038e-05,0.2869706279678e+01,0.1589072916335e+01],
             [0.3829101568895e-05,0.3572131359950e+01,0.7113454667900e-02],
             [0.3175733773908e-05,0.4535372530045e+01,0.4194847048887e+00],
             [0.3092437902159e-05,0.9230153317909e+00,0.6398972393349e+00],
             [0.2874168812154e-05,0.3363143761101e+01,0.1102062672231e+00],
             [0.3040119321826e-05,0.3324250895675e+01,0.6283075850446e+01],
             [0.2699723308006e-05,0.2917882441928e+00,0.1030928125552e+00],
             [0.2134832683534e-05,0.4220997202487e+01,0.3163918923335e+00],
             [0.1770412139433e-05,0.4747318496462e+01,0.1021328554739e+02],
             [0.1377264209373e-05,0.4305058462401e+00,0.1484170571900e-02],
             [0.1127814538960e-05,0.8538177240740e+00,0.6327837846670e+00],
             [0.1055608090130e-05,0.1551800742580e+01,0.4337116142245e+00],
             [0.9802673861420e-06,0.1459646735377e+01,0.1052268489556e+01],
             [0.1090329461951e-05,0.1587351228711e+01,0.1162474756779e+01],
             [0.6959590025090e-06,0.5534442628766e+01,0.1066495398892e+01],
             [0.5664914529542e-06,0.6030673003297e+01,0.9491756770005e+00],
             [0.6607787763599e-06,0.4989507233927e+01,0.8460828644453e+00],
             [0.6269725742838e-06,0.4222951804572e+01,0.1480791608091e+00],
             [0.6301889697863e-06,0.5444316669126e+01,0.2243449970715e+00],
             [0.4891042662861e-06,0.1490552839784e+01,0.3340612434717e+01],
             [0.3457083123290e-06,0.3030475486049e+01,0.3516457698740e-01],
             [0.3032559967314e-06,0.2652038793632e+01,0.1104591729320e-01],
             [0.2841133988903e-06,0.1276744786829e+01,0.4110125927500e-01],
             [0.2855564444432e-06,0.2143368674733e+01,0.1510475019529e+00],
             [0.2765157135038e-06,0.5444186109077e+01,0.6373574839730e-01],
             [0.2382312465034e-06,0.2190521137593e+01,0.2275259891141e+00],
             [0.2808060365077e-06,0.5735195064841e+01,0.2535050500000e-01],
             [0.2332175234405e-06,0.9481985524859e-01,0.7181332454670e-01],
             [0.2322488199659e-06,0.5180499361533e+01,0.8582758298370e-01],
             [0.1881850258423e-06,0.3219788273885e+01,0.2118763888447e+01],
             [0.2196111392808e-06,0.2366941159761e+01,0.2968341143800e-02],
             [0.2183810335519e-06,0.4825445110915e+01,0.7775000683430e-01],
             [0.2002733093326e-06,0.2457148995307e+01,0.2093666171530e+00],
             [0.1967111767229e-06,0.5586291545459e+01,0.2172315424036e+00],
             [0.1568473250543e-06,0.3708003123320e+01,0.7429900518901e+00],
             [0.1852528314300e-06,0.4310638151560e+01,0.2022531624851e+00],
             [0.1832111226447e-06,0.1494665322656e+01,0.3235053470014e+00],
             [0.1746805502310e-06,0.1451378500784e+01,0.1385174140878e+00],
             [0.1555730966650e-06,0.1068040418198e+01,0.7358765972222e+00],
             [0.1554883462559e-06,0.2442579035461e+01,0.5154640627760e+00],
             [0.1638380568746e-06,0.2597913420625e+00,0.8531963191132e+00],
             [0.1159938593640e-06,0.5834512021280e+01,0.1990721704425e+00],
             [0.1083427965695e-06,0.5054033177950e+01,0.5439178814476e+00],
             [0.1156480369431e-06,0.5325677432457e+01,0.5257585094865e+00],
             [0.1141308860095e-06,0.2153403923857e+01,0.5336234347371e+00],
             [0.7913146470946e-07,0.8642846847027e+00,0.1478866649112e+01],
             [0.7439752463733e-07,0.1970628496213e+01,0.2164800718209e+00],
             [0.7280277104079e-07,0.6073307250609e+01,0.2101180877357e+00],
             [0.8319567719136e-07,0.1954371928334e+01,0.1692165728891e+01],
             [0.7137705549290e-07,0.8904989440909e+00,0.4155522422634e+00],
             [0.6900825396225e-07,0.2825717714977e+01,0.1173197218910e+00],
             [0.7245757216635e-07,0.2481677513331e+01,0.1265567569334e+01],
             [0.6961165696255e-07,0.1292955312978e+01,0.9562891316684e+00],
             [0.7571804456890e-07,0.3427517575069e+01,0.1422690933580e-01],
             [0.6605425721904e-07,0.8052192701492e+00,0.6470106940028e+00],
             [0.7375477357248e-07,0.1705076390088e+01,0.1581959461667e+01],
             [0.7041664951470e-07,0.4848356967891e+00,0.9597935788730e-01],
             [0.6322199535763e-07,0.3878069473909e+01,0.7084920306520e-01],
             [0.5244380279191e-07,0.2645560544125e+01,0.5265099800692e+00],
             [0.5143125704988e-07,0.4834486101370e+01,0.5328719641544e+00],
             [0.5871866319373e-07,0.7981472548900e+00,0.7871412831580e-01],
             [0.6300822573871e-07,0.5979398788281e+01,0.2608790314060e+02],
             [0.6062154271548e-07,0.4108655402756e+01,0.1114304132498e+00],
             [0.4361912339976e-07,0.5322624319280e+01,0.1375773836557e+01],
             [0.4417005920067e-07,0.6240817359284e+01,0.2770348281756e+00],
             [0.4686806749936e-07,0.3214977301156e+01,0.1143987543936e+00],
             [0.3758892132305e-07,0.5879809634765e+01,0.1596186371003e+01],
             [0.5151351332319e-07,0.2893377688007e+00,0.2228608264996e+00],
             [0.4554683578572e-07,0.5475427144122e+01,0.1465949902372e+00],
             [0.3442381385338e-07,0.5992034796640e+01,0.5070101000000e-01],
             [0.2831093954933e-07,0.5367350273914e+01,0.3092784376656e+00],
             [0.3756267090084e-07,0.5758171285420e+01,0.4903339079539e+00],
             [0.2816374679892e-07,0.1863718700923e+01,0.2991266627620e+00],
             [0.3419307025569e-07,0.9524347534130e+00,0.3518164938661e+00],
             [0.2904250494239e-07,0.5304471615602e+01,0.1099462426779e+00],
             [0.2471734511206e-07,0.1297069793530e+01,0.6256703299991e+00],
             [0.2539620831872e-07,0.3281126083375e+01,0.1256615170089e+02],
             [0.2281017868007e-07,0.1829122133165e+01,0.6681224869435e+01],
             [0.2275319473335e-07,0.5797198160181e+01,0.3932462625300e-02],
             [0.2547755368442e-07,0.4752697708330e+01,0.1169588211447e+01],
             [0.2285979669317e-07,0.1223205292886e+01,0.1045155034888e+01],
             [0.1913386560994e-07,0.1757532993389e+01,0.1155361302111e+01],
             [0.1809020525147e-07,0.4246116108791e+01,0.3368040641550e-01],
             [0.1649213300201e-07,0.1445162890627e+01,0.4408250688924e+00],
             [0.1834972793932e-07,0.1126917567225e+01,0.4452511715700e-02],
             [0.1439550648138e-07,0.6160756834764e+01,0.9420622223326e+00],
             [0.1487645457041e-07,0.4358761931792e+01,0.4123712502208e+00],
             [0.1731729516660e-07,0.6134456753344e+01,0.2108507877249e+00],
             [0.1717747163567e-07,0.1898186084455e+01,0.2157473718317e+00],
             [0.1418190430374e-07,0.4180286741266e+01,0.6521991896920e-01],
             [0.1404844134873e-07,0.7654053565412e-01,0.4258542984690e-01],
             [0.1409842846538e-07,0.4418612420312e+01,0.2258291676434e+00],
             [0.1090948346291e-07,0.1260615686131e+01,0.4226656969313e+00],
             [0.1357577323612e-07,0.3558248818690e+01,0.7923417740620e-01],
             [0.1018154061960e-07,0.5676087241256e+01,0.1456308687557e+00],
             [0.1412073972109e-07,0.8394392632422e+00,0.1525316725248e+00],
             [0.1030938326496e-07,0.1653593274064e+01,0.1795258541446e+01],
             [0.1180081567104e-07,0.1285802592036e+01,0.7032915397480e-01],
             [0.9708510575650e-08,0.7631889488106e+00,0.8434341241180e-01],
             [0.9637689663447e-08,0.4630642649176e+01,0.1272681024002e+01],
             [0.1068910429389e-07,0.5294934032165e+01,0.2123349582968e+00],
             [0.1063716179336e-07,0.2736266800832e+01,0.2142632012598e+00],
             [0.1234858713814e-07,0.1302891146570e+01,0.1847279083684e+00],
             [0.8912631189738e-08,0.3570415993621e+01,0.2648454860559e+01],
             [0.1036378285534e-07,0.4236693440949e+01,0.1370332435159e+00],
             [0.9667798501561e-08,0.2960768892398e+01,0.4376440768498e+00],
             [0.8108314201902e-08,0.6987781646841e+00,0.2880807454688e+00],
             [0.7648364324628e-08,0.2499017863863e+01,0.2037373330570e+00],
             [0.7286136828406e-08,0.3787426951665e+01,0.1129145838217e+00],
             [0.9448237743913e-08,0.2694354332983e+01,0.5272426800584e+00],
             [0.9374276106428e-08,0.4787121277064e+01,0.5321392641652e+00],
             [0.7100226287462e-08,0.3530238792101e+00,0.6288513220417e+00],
             [0.9253056659571e-08,0.1399478925664e+01,0.1606092486742e+00],
             [0.6636432145504e-08,0.3479575438447e+01,0.1368660381889e+01],
             [0.6469975312932e-08,0.1383669964800e+01,0.2008557621224e+01],
             [0.7335849729765e-08,0.1243698166898e+01,0.9561746721300e-02],
             [0.8743421205855e-08,0.3776164289301e+01,0.3801276407308e+00],
             [0.5993635744494e-08,0.5627122113596e+01,0.2042657109477e+02],
             [0.5981008479693e-08,0.1674336636752e+01,0.2111650433779e+01],
             [0.6188535145838e-08,0.5214925208672e+01,0.4305306221819e+00],
             [0.6596074017566e-08,0.2907653268124e+01,0.1063314406849e+01],
             [0.6630815126226e-08,0.2127643669658e+01,0.8389694097774e+00],
             [0.6156772830040e-08,0.5082160803295e+01,0.4234171675140e+00],
             [0.6446960563014e-08,0.1872100916905e+01,0.5287268506303e+00],
             [0.6429324424668e-08,0.5610276103577e+01,0.5306550935933e+00],
             [0.6302232396465e-08,0.1592152049607e+01,0.1253008786510e-01],
             [0.6399244436159e-08,0.2746214421532e+01,0.5217580628120e+02],
             [0.5474965172558e-08,0.2317666374383e+01,0.2221856701002e+01],
             [0.5339293190692e-08,0.1084724961156e+01,0.7466759693650e-01],
             [0.5334733683389e-08,0.3594106067745e+01,0.7489573444450e-01],
             [0.5392665782110e-08,0.5630254365606e+01,0.1055449481598e+01],
             [0.6682075673789e-08,0.1518480041732e+01,0.2213766559277e+00],
             [0.5079130495960e-08,0.2739765115711e+01,0.2132517061319e+00],
             [0.5077759793261e-08,0.5290711290094e+01,0.2133464534247e+00],
             [0.4832037368310e-08,0.1404473217200e+01,0.7160067364790e-01],
             [0.6463279674802e-08,0.6038381695210e+01,0.2209183458640e-01],
             [0.6240592771560e-08,0.1290170653666e+01,0.3306188016693e+00],
             [0.4672013521493e-08,0.3261895939677e+01,0.7796265773310e-01],
             [0.6500650750348e-08,0.1154522312095e+01,0.3884652414254e+00],
             [0.6344161389053e-08,0.6206111545062e+01,0.7605151500000e-01],
             [0.4682518370646e-08,0.5409118796685e+01,0.1073608853559e+01],
             [0.5329460015591e-08,0.1202985784864e+01,0.7287631425543e+00],
             [0.5701588675898e-08,0.4098715257064e+01,0.8731175355560e-01],
             [0.6030690867211e-08,0.4132033218460e+00,0.9846002785331e+00],
             [0.4336256312655e-08,0.1211415991827e+01,0.4297791515992e+00],
             [0.4688498808975e-08,0.3765479072409e+01,0.2127790306879e+00],
             [0.4675578609335e-08,0.4265540037226e+01,0.2138191288687e+00],
             [0.4225578112158e-08,0.5237566010676e+01,0.3407705765729e+00],
             [0.5139422230028e-08,0.1507173079513e+01,0.7233337363710e-01],
             [0.4619995093571e-08,0.9023957449848e-01,0.8603097737811e+00],
             [0.4494776255461e-08,0.5414930552139e+00,0.7381754420900e-01],
             [0.4274026276788e-08,0.4145735303659e+01,0.7574578717200e-01],
             [0.5018141789353e-08,0.3344408829055e+01,0.3180992042600e-02],
             [0.4866163952181e-08,0.3348534657607e+01,0.7722995774390e-01],
             [0.4111986020501e-08,0.4198823597220e+00,0.1451108196653e+00],
             [0.3356142784950e-08,0.5609144747180e+01,0.1274714967946e+00],
             [0.4070575554551e-08,0.7028411059224e+00,0.3503323232942e+00],
             [0.3257451857278e-08,0.5624697983086e+01,0.5296435984654e+00],
             [0.3256973703026e-08,0.1857842076707e+01,0.5297383457582e+00],
             [0.3830771508640e-08,0.4562887279931e+01,0.9098186128426e+00],
             [0.3725024005962e-08,0.2358058692652e+00,0.1084620721060e+00],
             [0.3136763921756e-08,0.2049731526845e+01,0.2346394437820e+00],
             [0.3795147256194e-08,0.2432356296933e+00,0.1862120789403e+00],
             [0.2877342229911e-08,0.5631101279387e+01,0.1905464808669e+01],
             [0.3076931798805e-08,0.1117615737392e+01,0.3628624111593e+00],
             [0.2734765945273e-08,0.5899826516955e+01,0.2131850110243e+00],
             [0.2733405296885e-08,0.2130562964070e+01,0.2134131485323e+00],
             [0.2898552353410e-08,0.3462387048225e+00,0.5291709230214e+00],
             [0.2893736103681e-08,0.8534352781543e+00,0.5302110212022e+00],
             [0.3095717734137e-08,0.2875061429041e+01,0.2976424921901e+00],
             [0.2636190425832e-08,0.2242512846659e+01,0.1485980103780e+01],
             [0.3645512095537e-08,0.1354016903958e+01,0.6044726378023e+00],
             [0.2808173547723e-08,0.6705114365631e-01,0.6225157782540e-01],
             [0.2625012866888e-08,0.4775705748482e+01,0.5268983110410e-01],
             [0.2572233995651e-08,0.2638924216139e+01,0.1258454114666e+01],
             [0.2604238824792e-08,0.4826358927373e+01,0.2103781122809e+00],
             [0.2596886385239e-08,0.3200388483118e+01,0.2162200472757e+00],
             [0.3228057304264e-08,0.5384848409563e+01,0.2007689919132e+00],
             [0.2481601798252e-08,0.5173373487744e+01,0.1062562936266e+01],
             [0.2745977498864e-08,0.6250966149853e+01,0.5651155736444e+00],
             [0.2669878833811e-08,0.4906001352499e+01,0.1400015846597e+00],
             [0.3203986611711e-08,0.5034333010005e+01,0.7036329877322e+00],
             [0.3354961227212e-08,0.6108262423137e+01,0.4549093064213e+00],
             [0.2400407324558e-08,0.2135399294955e+01,0.2125476091956e+00],
             [0.2379905859802e-08,0.5893721933961e+01,0.2140505503610e+00],
             [0.2550844302187e-08,0.3331940762063e+01,0.1534957940063e+00],
             [0.2268824211001e-08,0.1843418461035e+01,0.2235935264888e+00],
             [0.2464700891204e-08,0.3029548547230e+01,0.2091065926078e+00],
             [0.2436814726024e-08,0.4994717970364e+01,0.2174915669488e+00],
             [0.2443623894745e-08,0.2645102591375e+01,0.1739420156204e+00],
             [0.2318701783838e-08,0.5700547397897e+01,0.7530171478090e-01],
             [0.2284448700256e-08,0.5268898905872e+01,0.7426161660010e-01],
             [0.2468848123510e-08,0.5276280575078e+01,0.2526561439362e+00],
             [0.2814052350303e-08,0.6130168623475e+01,0.5636314030725e+00],
             [0.2243662755220e-08,0.6631692457995e+00,0.8886590321940e-01],
             [0.2330795855941e-08,0.2499435487702e+01,0.1056200952181e+01],
             [0.9757679038404e-09,0.5796846023126e+01,0.7826370942180e+02]]

    #SSB-to-Sun, T^1, Y
    S1[1]=[[0.8989047573576e-08,0.5840593672122e+01,0.4265981595566e+00],
             [0.7815938401048e-08,0.1129664707133e+01,0.2061856251104e+00],
             [0.7550926713280e-08,0.6196589104845e+00,0.2204125344462e+00],
             [0.6056556925895e-08,0.1677494667846e+01,0.1059381944224e+01],
             [0.5734142698204e-08,0.4000920852962e+01,0.5225775174439e+00],
             [0.5614341822459e-08,0.3486722577328e+01,0.5368044267797e+00],
             [0.1028678147656e-08,0.1877141024787e+01,0.7113454667900e-02],
             [0.7270792075266e-09,0.5077167301739e+01,0.6398972393349e+00],
             [0.8734141726040e-09,0.9069550282609e-01,0.4194847048887e+00],
             [0.5377371402113e-09,0.6039381844671e+01,0.4337116142245e+00],
             [0.4729719431571e-09,0.2153086311760e+01,0.2132990797783e+00],
             [0.4458052820973e-09,0.5059830025565e+01,0.5296909721118e+00],
             [0.4406855467908e-09,0.2027971692630e+01,0.1589072916335e+01],
             [0.3101659310977e-09,0.3317677981860e+01,0.1052268489556e+01],
             [0.3016749232545e-09,0.3913703482532e+01,0.1066495398892e+01],
             [0.3198541352656e-09,0.1275513098525e+01,0.1495633313810e+00],
             [0.2142065389871e-09,0.5301351614597e+01,0.3163918923335e+00],
             [0.1902615247592e-09,0.4894943352736e+00,0.2275259891141e+00],
             [0.1613410990871e-09,0.2449891130437e+01,0.1102062672231e+00],
             [0.1576992165097e-09,0.4211421447633e+01,0.7626583626240e-01],
             [0.1241637259894e-09,0.4140803368133e+01,0.5154640627760e+00],
             [0.1313974830355e-09,0.3591920305503e+01,0.3664874755930e-01],
             [0.1181697118258e-09,0.1506314382788e+01,0.6327837846670e+00],
             [0.1238239742779e-09,0.7461405378404e+00,0.3961708870310e-01],
             [0.1010107068241e-09,0.6271010795475e+00,0.7329749511860e-01],
             [0.9226316616509e-10,0.1259158839583e+01,0.1990721704425e+00],
             [0.8664946419555e-10,0.3353244696934e+01,0.5439178814476e+00],
             [0.7757230468978e-10,0.1447677295196e+01,0.9491756770005e+00],
             [0.7693168628139e-10,0.1120509896721e+01,0.1030928125552e+00],
             [0.5487897454612e-10,0.4439380426795e+01,0.8531963191132e+00],
             [0.5196118677218e-10,0.3788856619137e+00,0.2093666171530e+00],
             [0.5110853339935e-10,0.1386879372016e+01,0.2172315424036e+00],
             [0.5027804534813e-10,0.1647881805466e+00,0.2164800718209e+00],
             [0.4922485922674e-10,0.1594315079862e+01,0.2101180877357e+00],
             [0.6155599524400e-10,0.0000000000000e+00,0.0000000000000e+00],
             [0.4447147832161e-10,0.5480720918976e+01,0.3235053470014e+00],
             [0.4144691276422e-10,0.1931371033660e+01,0.6373574839730e-01],
             [0.4099950625452e-10,0.5229611294335e+01,0.6470106940028e+00],
             [0.5060541682953e-10,0.1731112486298e+01,0.1422690933580e-01],
             [0.4293615946300e-10,0.2714571038925e+01,0.7358765972222e+00],
             [0.3545659845763e-10,0.4451041444634e+01,0.5265099800692e+00],
             [0.3479112041196e-10,0.3029385448081e+01,0.5328719641544e+00],
             [0.3438516493570e-10,0.2778507143731e+01,0.8582758298370e-01],
             [0.3297341285033e-10,0.7898709807584e+00,0.1104591729320e-01],
             [0.2972585818015e-10,0.3218785316973e+01,0.5257585094865e+00],
             [0.2931707295017e-10,0.4260731012098e+01,0.5336234347371e+00],
             [0.2897198149403e-10,0.1120753978101e+01,0.1173197218910e+00],
             [0.2832293240878e-10,0.4597682717827e+00,0.2022531624851e+00],
             [0.2864348326612e-10,0.2169939928448e+01,0.9597935788730e-01],
             [0.2852714675471e-10,0.2377659870578e+01,0.2118763888447e+01]]
        
    #SSB-to-Sun, T^2, Y
    S2[1]=[[0.1609114495091e-11,0.2831096993481e+01,0.2061856251104e+00],
             [0.1560330784946e-11,0.5193058213906e+01,0.2204125344462e+00],
             [0.1183535479202e-11,0.5707003443890e+01,0.5225775174439e+00],
             [0.1158183066182e-11,0.1782400404928e+01,0.5368044267797e+00],
             [0.1032868027407e-11,0.4036925452011e+01,0.2132990797783e+00],
             [0.6540142847741e-12,0.4058241056717e+01,0.4265981595566e+00],
             [0.7305236491596e-12,0.6175401942957e+00,0.5296909721118e+00],
             [-0.5580725052968e-12,0.0000000000000e+00,0.0000000000000e+00],
             [0.3946122651015e-12,0.4108265279171e+00,0.1059381944224e+01]]
    
    #SSB-to-Sun, T^0, Z
    S0[2]=[[0.1181255122986e-03,0.4607918989164e+00,0.2132990797783e+00],
             [0.1127777651095e-03,0.4169146331296e+00,0.5296909721118e+00],
             [0.4777754401806e-04,0.4582657007130e+01,0.3813291813120e-01],
             [0.1129354285772e-04,0.5758735142480e+01,0.7478166569050e-01],
             [-0.1149543637123e-04,0.0000000000000e+00,0.0000000000000e+00],
             [0.3298730512306e-05,0.5978801994625e+01,0.4265981595566e+00],
             [0.2733376706079e-05,0.7665413691040e+00,0.1059381944224e+01],
             [0.9426389657270e-06,0.3710201265838e+01,0.2061856251104e+00],
             [0.8187517749552e-06,0.3390675605802e+00,0.2204125344462e+00],
             [0.4080447871819e-06,0.4552296640088e+00,0.5225775174439e+00],
             [0.3169973017028e-06,0.3445455899321e+01,0.5368044267797e+00],
             [0.2438098615549e-06,0.5664675150648e+01,0.3664874755930e-01],
             [0.2601897517235e-06,0.1931894095697e+01,0.1495633313810e+00],
             [0.2314558080079e-06,0.3666319115574e+00,0.3961708870310e-01],
             [0.1962549548002e-06,0.3167411699020e+01,0.7626583626240e-01],
             [0.2180518287925e-06,0.1544420746580e+01,0.7113454667900e-02],
             [0.1451382442868e-06,0.1583756740070e+01,0.1102062672231e+00],
             [0.1358439007389e-06,0.5239941758280e+01,0.6398972393349e+00],
             [0.1050585898028e-06,0.2266958352859e+01,0.3163918923335e+00],
             [0.1050029870186e-06,0.2711495250354e+01,0.4194847048887e+00],
             [0.9934920679800e-07,0.1116208151396e+01,0.1589072916335e+01],
             [0.1048395331560e-06,0.3408619600206e+01,0.1021328554739e+02],
             [0.8370147196668e-07,0.3810459401087e+01,0.2535050500000e-01],
             [0.7989856510998e-07,0.3769910473647e+01,0.7329749511860e-01],
             [0.5441221655233e-07,0.2416994903374e+01,0.1030928125552e+00],
             [0.4610812906784e-07,0.5858503336994e+01,0.4337116142245e+00],
             [0.3923022803444e-07,0.3354170010125e+00,0.1484170571900e-02],
             [0.2610725582128e-07,0.5410600646324e+01,0.6327837846670e+00],
             [0.2455279767721e-07,0.6120216681403e+01,0.1162474756779e+01],
             [0.2375530706525e-07,0.6055443426143e+01,0.1052268489556e+01],
             [0.1782967577553e-07,0.3146108708004e+01,0.8460828644453e+00],
             [0.1581687095238e-07,0.6255496089819e+00,0.3340612434717e+01],
             [0.1594657672461e-07,0.3782604300261e+01,0.1066495398892e+01],
             [0.1563448615040e-07,0.1997775733196e+01,0.2022531624851e+00],
             [0.1463624258525e-07,0.1736316792088e+00,0.3516457698740e-01],
             [0.1331585056673e-07,0.4331941830747e+01,0.9491756770005e+00],
             [0.1130634557637e-07,0.6152017751825e+01,0.2968341143800e-02],
             [0.1028949607145e-07,0.2101792614637e+00,0.2275259891141e+00],
             [0.1024074971618e-07,0.4071833211074e+01,0.5070101000000e-01],
             [0.8826956060303e-08,0.4861633688145e+00,0.2093666171530e+00],
             [0.8572230171541e-08,0.5268190724302e+01,0.4110125927500e-01],
             [0.7649332643544e-08,0.5134543417106e+01,0.2608790314060e+02],
             [0.8581673291033e-08,0.2920218146681e+01,0.1480791608091e+00],
             [0.8430589300938e-08,0.3604576619108e+01,0.2172315424036e+00],
             [0.7776165501012e-08,0.3772942249792e+01,0.6373574839730e-01],
             [0.8311070234408e-08,0.6200412329888e+01,0.3235053470014e+00],
             [0.6927365212582e-08,0.4543353113437e+01,0.8531963191132e+00],
             [0.6791574208598e-08,0.2882188406238e+01,0.7181332454670e-01],
             [0.5593100811839e-08,0.1776646892780e+01,0.7429900518901e+00],
             [0.4553381853021e-08,0.3949617611240e+01,0.7775000683430e-01],
             [0.5758000450068e-08,0.3859251775075e+01,0.1990721704425e+00],
             [0.4281283457133e-08,0.1466294631206e+01,0.2118763888447e+01],
             [0.4206935661097e-08,0.5421776011706e+01,0.1104591729320e-01],
             [0.4213751641837e-08,0.3412048993322e+01,0.2243449970715e+00],
             [0.5310506239878e-08,0.5421641370995e+00,0.5154640627760e+00],
             [0.3827450341320e-08,0.8887314524995e+00,0.1510475019529e+00],
             [0.4292435241187e-08,0.1405043757194e+01,0.1422690933580e-01],
             [0.3189780702289e-08,0.1060049293445e+01,0.1173197218910e+00],
             [0.3226611928069e-08,0.6270858897442e+01,0.2164800718209e+00],
             [0.2893897608830e-08,0.5117563223301e+01,0.6470106940028e+00],
             [0.3239852024578e-08,0.4079092237983e+01,0.2101180877357e+00],
             [0.2956892222200e-08,0.1594917021704e+01,0.3092784376656e+00],
             [0.2980177912437e-08,0.5258787667564e+01,0.4155522422634e+00],
             [0.3163725690776e-08,0.3854589225479e+01,0.8582758298370e-01],
             [0.2662262399118e-08,0.3561326430187e+01,0.5257585094865e+00],
             [0.2766689135729e-08,0.3180732086830e+00,0.1385174140878e+00],
             [0.2411600278464e-08,0.3324798335058e+01,0.5439178814476e+00],
             [0.2483527695131e-08,0.4169069291947e+00,0.5336234347371e+00],
             [0.7788777276590e-09,0.1900569908215e+01,0.5217580628120e+02]]

    #SSB-to-Sun, T^1, Z
    S1[2]=[[0.5444220475678e-08,0.1803825509310e+01,0.2132990797783e+00],
             [0.3883412695596e-08,0.4668616389392e+01,0.5296909721118e+00],
             [0.1334341434551e-08,0.0000000000000e+00,0.0000000000000e+00],
             [0.3730001266883e-09,0.5401405918943e+01,0.2061856251104e+00],
             [0.2894929197956e-09,0.4932415609852e+01,0.2204125344462e+00],
             [0.2857950357701e-09,0.3154625362131e+01,0.7478166569050e-01],
             [0.2499226432292e-09,0.3657486128988e+01,0.4265981595566e+00],
             [0.1937705443593e-09,0.5740434679002e+01,0.1059381944224e+01],
             [0.1374894396320e-09,0.1712857366891e+01,0.5368044267797e+00],
             [0.1217248678408e-09,0.2312090870932e+01,0.5225775174439e+00],
             [0.7961052740870e-10,0.5283368554163e+01,0.3813291813120e-01],
             [0.4979225949689e-10,0.4298290471860e+01,0.4194847048887e+00],
             [0.4388552286597e-10,0.6145515047406e+01,0.7113454667900e-02],
             [0.2586835212560e-10,0.3019448001809e+01,0.6398972393349e+00]]
    
    #SSB-to-Sun, T^2, Z
    S2[2]=[[0.3749920358054e-12,0.3230285558668e+01,0.2132990797783e+00],
             [0.2735037220939e-12,0.6154322683046e+01,0.5296909721118e+00]]

    #相对于参考历元的时间间隔，单位：年
    T=((DATE1-DJ00)+DATE2)/DJY
    T2=T*T
    
    #判断日期是否有误.
    if (np.abs(T)<=100.0):
        JSTAT=0
    else:
        JSTAT=1
        
    #对三个纬度的数据进行处理（X,Y,Z）
    for K in range(3):
        
    #初始化位置和速度参数
        XYZ=0.0
        XYZD=0.0
    
    # ------------------------------------------------
    #获得日心相对于地球的黄道向量
    # ------------------------------------------------

    # Sun to Earth, T^0 terms.
        for J in range(NE0[K]):
            A=E0[K][J][0]
            B=E0[K][J][1]
            C=E0[K][J][2]
            P=B+C*T
            XYZ=XYZ+A*ma.cos(P)
            XYZD=XYZD-A*C*ma.sin(P)
            
    #Sun to Earth, T^1 terms.
        for J in range(NE1[K]):
            A=E1[K][J][0]
            B=E1[K][J][1]
            C=E1[K][J][2]
            CT=C*T
            P=B+CT
            CP=ma.cos(P)
            XYZ=XYZ+A*T*CP
            XYZD=XYZD+A*(CP-CT*ma.sin(P))
        
    #Sun to Earth, T^2 terms. 
        for J in range(NE2[K]):
            A=E2[K][J][0]
            B=E2[K][J][1]
            C=E2[K][J][2]   
            CT=C*T
            P=B+CT
            CP=ma.cos(P)
            XYZ=XYZ+A*T2*CP
            XYZD=XYZD+A*T*(2.0*CP-CT*ma.sin(P))
            
    #日心相对于地球的位置速度参数
        PH[K]=XYZ
        VH[K]=XYZD/DJY
        
    # ------------------------------------------------
    # 获得质心相对与地球的黄道向量
    # ------------------------------------------------
    
    #SSB to Sun, T^0 terms.
        for J in range(NS0[K]):
            A=S0[K][J][0]
            B=S0[K][J][1]
            C=S0[K][J][2]
            P=B+C*T
            XYZ=XYZ+A*ma.cos(P)
            XYZD=XYZD-A*C*ma.sin(P)
            
    #SSB to Sun, T^1 terms.
        for J in range(NS1[K]):    
            A=S1[K][J][0]
            B=S1[K][J][1]
            C=S1[K][J][2]
            CT=C*T
            P=B+CT
            CP=ma.cos(P)
            XYZ=XYZ+A*T*CP
            XYZD=XYZD+A*(CP-CT*ma.sin(P))
    
    #SSB to Sun, T^2 terms.
        for J in range(NS2[K]):
            A=S2[K][J][0]
            B=S2[K][J][1]
            C=S2[K][J][2]   
            CT=C*T
            P=B+CT
            CP=ma.cos(P)
            XYZ=XYZ+A*T2*CP
            XYZD=XYZD+A*T*(2.0*CP-CT*ma.sin(P))
    
    #质心相对与地球的位置速度参数
        PB[K]=XYZ
        VB[K]=XYZD/DJY
    
    #通过旋转，从黄道变换到ICRF坐标 
    X=PH[0]
    Y=PH[1]
    Z=PH[2]
    PVH[0][0]=X+AM12*Y+AM13*Z
    PVH[0][1]=AM21*X+AM22*Y+AM23*Z
    PVH[0][2]=AM32*Y+AM33*Z
    X=VH[0]
    Y=VH[1]
    Z=VH[2]
    PVH[1][0]=X+AM12*Y+AM13*Z
    PVH[1][1]=AM21*X+AM22*Y+AM23*Z
    PVH[1][2]=AM32*Y+AM33*Z
    X=PB[0]
    Y=PB[1]
    Z=PB[2]
    PVB[0][0]=X+AM12*Y+AM13*Z
    PVB[0][1]=AM21*X+AM22*Y+AM23*Z
    PVB[0][2]=AM32*Y+AM33*Z
    X=VB[0]
    Y=VB[1]
    Z=VB[2]
    PVB[1][0]=X+AM12*Y+AM13*Z
    PVB[1][1]=AM21*X+AM22*Y+AM23*Z
    PVB[1][2]=AM32*Y+AM33*Z
        
    return(PVH,PVB,JSTAT)


def pymApcg13(DATE1,DATE2,ASTROM):
    '''
    For a geocentric observer, prepare star-independent astrometry
    parameters for transformations between ICRS and GCRS coordinates.
    The caller supplies the date, and SOFA models are used to predict
    the Earth ephemeris.
    
    The parameters produced by this function are required in the
    parallax, light deflection and aberration parts of the astrometric
    transformation chain.

    Parameters
    ----------
    date1 : float
        TDB as a 2-part Julian Date    
    date2 : float
        TDB as a 2-part Julian Date    

    Returns
    -------
    astrom : list(30)
        star-independent astrometry parameters    
        [0]>pmt : PM time interval (SSB, Julian years)
        [1-3]>eb : SSB to observer (vector, au)
        [4-6]>eh : Sun to observer (unit vector)
        [7]>em : distance from Sun to observer (au)
        [8-10]>v : barycentric observer velocity (vector, c)
        [11]>bm1 : sqrt(1-|v|^2): reciprocal of Lorenz factor
        [12-20]>bpn : bias-precession-nutation matrix
        [21]>along : unchanged
        [22]>xpl : unchanged
        [23]>ypl : unchanged
        [24]>sphi : unchanged
        [25]>cphi : unchanged
        [26]>diurab : unchanged
        [27]>eral : unchanged
        [28]>refa : unchanged 
        [29]>refb : unchanged

    '''
    #地球相对于太阳系质心以及日心的位置、速度信息 (au, au/天).
    EHPV,EBPV,J=pymEpv00(DATE1,DATE2)
    EHP = EHPV[0]
    #调用pymApcg函数获得ASTROM
    ASTROM=pymApcg(DATE1,DATE2,EBPV,EHP,ASTROM)
    
    return(ASTROM)


def pymRx(PHI,R):
    '''
    Rotate an r-matrix about the x-axis.

    Parameters
    ----------
    phi : float    
        angle (radians)    
    r : list(3,3)    
        r-matrix    

    Returns
    -------
    r : list(3,3)    
        r-matrix, rotated

    '''
    S=ma.sin(PHI)
    C=ma.cos(PHI)
    
    A21=C*R[1][0]+S*R[2][0]
    A22=C*R[1][1]+S*R[2][1]
    A23=C*R[1][2]+S*R[2][2]
    A31=-S*R[1][0]+C*R[2][0]
    A32=-S*R[1][1]+C*R[2][1]
    A33=-S*R[1][2]+C*R[2][2]

    R[1][0]=A21
    R[1][1]=A22
    R[1][2]=A23
    R[2][0]=A31
    R[2][1]=A32
    R[2][2]=A33
    
    return(R)


def pymRy(THETA,R):
    '''
    Rotate an r-matrix about the y-axis.

    Parameters
    ----------
    theta : float
        angle (radians)    
    r : list(3,3)
        r-matrix

    Returns
    -------
    r : list(3,3)
        r-matrix, rotated

    '''
    S=ma.sin(THETA)
    C=ma.cos(THETA)

    A11=C*R[0][0]-S*R[2][0]
    A12=C*R[0][1]-S*R[2][1]
    A13=C*R[0][2]-S*R[2][2]
    A31=S*R[0][0]+C*R[2][0]
    A32=S*R[0][1]+C*R[2][1]
    A33=S*R[0][2]+C*R[2][2]

    R[0][0]=A11
    R[0][1]=A12
    R[0][2]=A13
    R[2][0]=A31
    R[2][1]=A32
    R[2][2]=A33
    
    return(R)


def pymRz(PSI,R):
    '''
    Rotate an r-matrix about the z-axis.

    Parameters
    ----------
    psi : float
        angle (radians)    
    r : list(3,3)
        r-matrix

    Returns
    -------
    r : list(3,3)
        r-matrix, rotated

    '''
    S=ma.sin(PSI)
    C=ma.cos(PSI)

    A11=C*R[0][0]+S*R[1][0]
    A12=C*R[0][1]+S*R[1][1]
    A13=C*R[0][2]+S*R[1][2]
    A21=-S*R[0][0]+C*R[1][0]
    A22=-S*R[0][1]+C*R[1][1]
    A23=-S*R[0][2]+C*R[1][2]

    R[0][0]=A11
    R[0][1]=A12
    R[0][2]=A13
    R[1][0]=A21
    R[1][1]=A22
    R[1][2]=A23
 
    return(R)


def pymC2ixys(X,Y,S):
    '''
    Form the celestial to intermediate-frame-of-date matrix given the CIP
    X,Y and the CIO locator s.

    Parameters
    ----------
    x : float
        Celestial Intermediate Pole    
    y : float
        Celestial Intermediate Pole    
    s : float
        the CIO locator s 

    Returns
    -------
    rc2i : list(3,3)
        celestial-to-intermediate matrix

    '''
    #得到球面角E和d
    R2=X**2+Y**2
    if (R2>0.0):
        E=ma.atan2(Y,X)
    else:
        E=0.0
    D=ma.atan(ma.sqrt(R2/(1.0-R2)))
    
    #构建矩阵
    RC2I=pymIr()
    RC2I=pymRz(E,RC2I)
    RC2I=pymRy(D,RC2I)
    RC2I=pymRz(-(E+S),RC2I)    
    
    return(RC2I)


def pymApci(DATE1,DATE2,EBPV,EHP,X,Y,S,ASTROM):
    '''
    For a terrestrial observer, prepare star-independent astrometry
    parameters for transformations between ICRS and geocentric CIRS
    coordinates.  The Earth ephemeris and CIP/CIO are supplied by the
    caller.
   
    The parameters produced by this function are required in the
    parallax, light deflection, aberration, and bias-precession-nutation
    parts of the astrometric transformation chain.

    Parameters
    ----------
    date1 : float
        TDB as a 2-part Julian Date    
    date2 : float
        TDB as a 2-part Julian Date     
    ebpv : list(2,3)
        Earth barycentric position/velocity (au, au/day)     
    ehp : list(3)
        Earth heliocentric position (au)    
    x : flaot
        CIP X,Y (components of unit vector)    
    y : float
        CIP X,Y (components of unit vector)    
    s : float
        the CIO locator s (radians)

    Returns
    -------
    astrom : list(30)
        star-independent astrometry parameters    
        [0]>pmt : PM time interval (SSB, Julian years)    
        [1-3]>eb : SSB to observer (vector, au)    
        [4-6]>eh : Sun to observer (unit vector)     
        [7]>em : distance from Sun to observer (au)    
        [8-10]>v : barycentric observer velocity (vector, c)    
        [11]>bm1 : sqrt(1-|v|^2): reciprocal of Lorenz factor    
        [12-20]>bpn : bias-precession-nutation matrix    
        [21]>along : unchanged    
        [22]>xpl : unchanged    
        [23]>ypl : unchanged    
        [24]>sphi : unchanged    
        [25]>cphi : unchanged    
        [26]>diurab : unchanged    
        [27]>eral : unchanged    
        [28]>refa : unchanged     
        [29]>refb : unchanged    

    '''
    #在地心的与恒星无关的天体测量参数.
    ASTROM=pymApcg(DATE1,DATE2,EBPV,EHP,ASTROM)
    
    #基于CIO的BPN矩阵.
    A=pymC2ixys(X,Y,S)
    ASTROM[12]=A[0][0]
    ASTROM[13]=A[0][1]
    ASTROM[14]=A[0][2]
    ASTROM[15]=A[1][0]
    ASTROM[16]=A[1][1]
    ASTROM[17]=A[1][2]
    ASTROM[18]=A[2][0]
    ASTROM[19]=A[2][1]
    ASTROM[20]=A[2][2]
    
    return(ASTROM)


def pymFw2m(GAMB,PHIB,PSI,EPS):
    '''
    Form rotation matrix given the Fukushima-Williams angles.

    Parameters
    ----------
    gamb : float
        F-W angle gamma_bar (radians)    
    phib : float
        F-W angle phi_bar (radians)    
    psi : float
        F-W angle psi (radians)    
    eps : float
        F-W angle epsilon (radians)

    Returns
    -------
    r : list(3,3)
        rotation matrix

    '''
    R=pymIr()
    R=pymRz(GAMB, R)
    R=pymRx(PHIB, R)
    R=pymRz(-PSI, R)
    R=pymRx(-EPS, R)
   
    return(R)


def pymFad03(T):
    '''
    Fundamental argument, IERS Conventions (2003):
    mean elongation of the Moon from the Sun.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0
 
    Returns
    -------
    function value : float
        D, radians

    '''
    #1角秒对应的弧度
    DAS2R=4.848136811095359935899141e-6
    
    #一个圆内的总角秒数
    TURNAS=1296000.0
    
    #月亮到太阳的平均距角 
    A=1072260.703692+T*(1602961601.2090+T*(-6.3706+\
            T*(0.006593+T*(-0.00003169))))
    FAD=(A%TURNAS)*DAS2R
    
    return(FAD)


def pymFae03(T):
    '''
    Fundamental argument, IERS Conventions (2003):
    mean longitude of Earth.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0
 
    Returns
    -------
    function value : float
        mean longitude of Earth, radians

    '''
    FAE=(1.753470314+628.3075849991*T)%(2*ma.pi)    
    
    return(FAE)


def pymFaf03(T):
    '''
    Fundamental argument, IERS Conventions (2003):
    mean longitude of the Moon minus mean longitude of the ascending
    node.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0
 
    Returns
    -------
    function value : float
        F, radians

    '''
    #1角秒对应的弧度
    DAS2R=4.848136811095359935899141e-6
    
    #一个圆内的总角秒数
    TURNAS=1296000.0
    
    #月球的平黄经减去其上升点的平黄经 
    A=335779.526232+T*(1739527262.8478+T*(-12.7512+\
        T*(-0.001037+T*(0.00000417))))
    FAF=(A%TURNAS)*DAS2R
    
    return(FAF)


def pymFaju03(T):
    '''
    Fundamental argument, IERS Conventions (2003):
    mean longitude of Jupiter.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0
 
    Returns
    -------
    function value : float
        mean longitude of Jupiter, radians

    '''
    FAJU=(0.599546497+52.9690962641*T)%(2*ma.pi)

    return(FAJU)


def pymFal03(T):
    '''
    Fundamental argument, IERS Conventions (2003):
    mean anomaly of the Moon.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0
 
    Returns
    -------
    function value : float
        l, radians

    '''
    #1角秒对应的弧度
    DAS2R=4.848136811095359935899141e-6
    
    #一个圆内的总角秒数
    TURNAS=1296000.0
    
    #月亮的平近点角
    A=485868.249036+T*(1717915923.2178+T*(31.8792+\
        T*(0.051635+T*(-0.00024470))))
    FAL=(A%TURNAS)*DAS2R
    
    return(FAL)


def pymFalp03(T):
    '''
    Fundamental argument, IERS Conventions (2003):
    mean anomaly of the Sun.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0
 
    Returns
    -------
    function value : float
        l', radians

    '''
    #1角秒对应的弧度
    DAS2R=4.848136811095359935899141e-6
    
    #一个圆内的总角秒数
    TURNAS=1296000.0
    
    #太阳的平近点角
    A=1287104.793048+T*(129596581.0481+T*(-0.5532+\
        T*(0.000136+T*(-0.00001149))))
    FALP=(A%TURNAS)*DAS2R

    return(FALP)


def pymFama03(T):
    '''
    Fundamental argument, IERS Conventions (2003):
    mean longitude of Mars.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0
 
    Returns
    -------
    function value : float
        mean longitude of Mars, radians

    '''
    FAMA=(6.203480913+334.0612426700*T)%(2*ma.pi)
    
    return(FAMA)


def pymFame03(T):
    '''
    Fundamental argument, IERS Conventions (2003):
    mean longitude of Mercury.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0
 
    Returns
    -------
    function value : float
        mean longitude of Mercury, radians

    '''
    
    FAME=(4.402608842+2608.7903141574*T)%(2*ma.pi)
    
    return(FAME)


def pymFane03(T):
    '''
    Fundamental argument, IERS Conventions (2003):
    mean longitude of Neptune.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0

    Returns
    -------
    function value : float
        mean longitude of Neptune, radians

    '''    
    FANE=(5.311886287+3.8133035638*T)%(2*ma.pi)
    
    return(FANE)


def pymFaom03(T):
    '''
    Fundamental argument, IERS Conventions (2003):
    mean longitude of the Moon's ascending node.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0

    Returns
    -------
    function value : float
        Omega, radians 

    '''
    #1角秒对应的弧度
    DAS2R=4.848136811095359935899141e-6
    
    #一个圆内的总角秒数
    TURNAS=1296000.0
    
    #月亮升交点的平黄经
    A=450160.398036+T*(-6962890.5431+T*(7.4722+\
        T*(0.007702+T*(-0.00005939))))
    FAOM=(A%TURNAS)*DAS2R
    
    return(FAOM)


def pymFapa03(T):
    '''
    Fundamental argument, IERS Conventions (2003):
    general accumulated precession in longitude.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0

    Returns
    -------
    function value : float
        general precession in longitude, radians

    '''
    FAPA=(0.024381750+0.00000538691*T)*T
    
    return(FAPA)


def pymFasa03(T):
    '''
    Fundamental argument, IERS Conventions (2003):
    mean longitude of Saturn.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0

    Returns
    -------
    function value : float
        mean longitude of Saturn, radians

    '''
    FASA=(0.874016757+21.3299104960*T)%(2*ma.pi)
    
    return(FASA)


def pymFaur03(T):
    '''
    Fundamental argument, IERS Conventions (2003):
    mean longitude of Uranus.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0

    Returns
    -------
    function value : float
        mean longitude of Uranus, radians

    '''
    FAUR=(5.481293872+7.4781598567*T)%(2*ma.pi)
    
    return(FAUR)


def pymFave03(T):
    '''
    Fundamental argument, IERS Conventions (2003):
    mean longitude of Venus.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0

    Returns
    -------
    function value : float
        mean longitude of Venus, radians

    '''
    FAVE=(3.176146697+1021.3285546211*T)%(2*ma.pi)
    
    return(FAVE)


def pymNut00a(DATE1,DATE2):
    '''
    Nutation, IAU 2000A model (MHB2000 luni-solar and planetary nutation
    with free core nutation omitted).

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    dpsi : float
        nutation, luni-solar + planetary    
    deps : float
        nutation, luni-solar + planetary

    '''
    #1角秒对应的弧度
    DAS2R=4.848136811095359935899141e-6
    
    #一个圆内的总角秒数
    TURNAS=1296000.0
    
    #2Pi
    D2PI = 6.283185307179586476925287
    
    #0.1微角秒对应的弧度
    U2R=DAS2R/1e7
    
    #参考历元(J2000.0), JD
    DJ00=2451545.0
    
    #1儒略世纪对应的天数
    DJC=36525.0
    
    #日月章动模型
    #模型中的项数
    NLS=678
    
    #模型中的项数
    #行星章动模型
    NPL=687
    
    #日月章动的参数和项系数的表格
    #  L     L'    F     D     Om
    NALS=[[0,0,0,0,1],
            [0,0,2,-2,2],
            [0,0,2,0,2],
            [0,0,0,0,2],
            [0,1,0,0,0],
            [0,1,2,-2,2],
            [1,0,0,0,0],
            [0,0,2,0,1],
            [1,0,2,0,2],
            [0,-1,2,-2,2],
            [0,0,2,-2,1],
            [-1,0,2,0,2],
            [-1,0,0,2,0],
            [1,0,0,0,1],
            [-1,0,0,0,1],
            [-1,0,2,2,2],
            [1,0,2,0,1],
            [-2,0,2,0,1],
            [0,0,0,2,0],
            [0,0,2,2,2],
            [0,-2,2,-2,2],
            [-2,0,0,2,0],
            [2,0,2,0,2],
            [1,0,2,-2,2],
            [-1,0,2,0,1],
            [2,0,0,0,0],
            [0,0,2,0,0],
            [0,1,0,0,1],
            [-1,0,0,2,1],
            [0,2,2,-2,2],
            [0,0,-2,2,0],
            [1,0,0,-2,1],
            [0,-1,0,0,1],
            [-1,0,2,2,1],
            [0,2,0,0,0],
            [1,0,2,2,2],
            [-2,0,2,0,0],
            [0,1,2,0,2],
            [0,0,2,2,1],
            [0,-1,2,0,2],
            [0,0,0,2,1],
            [1,0,2,-2,1],
            [2,0,2,-2,2],
            [-2,0,0,2,1],
            [2,0,2,0,1],
            [0,-1,2,-2,1],
            [0,0,0,-2,1],
            [-1,-1,0,2,0],
            [2,0,0,-2,1],
            [1,0,0,2,0],
            [0,1,2,-2,1],
            [1,-1,0,0,0],
            [-2,0,2,0,2],
            [3,0,2,0,2],
            [0,-1,0,2,0],
            [1,-1,2,0,2],
            [0,0,0,1,0],
            [-1,-1,2,2,2],
            [-1,0,2,0,0],
            [0,-1,2,2,2],
            [-2,0,0,0,1],
            [1,1,2,0,2],
            [2,0,0,0,1],
            [-1,1,0,1,0],
            [1,1,0,0,0],
            [1,0,2,0,0],
            [-1,0,2,-2,1],
            [1,0,0,0,2],
            [-1,0,0,1,0],
            [0,0,2,1,2],
            [-1,0,2,4,2],
            [-1,1,0,1,1],
            [0,-2,2,-2,1],
            [1,0,2,2,1],
            [-2,0,2,2,2],
            [-1,0,0,0,2],
            [1,1,2,-2,2],
            [-2,0,2,4,2],
            [-1,0,4,0,2],
            [2,0,2,-2,1],
            [2,0,2,2,2],
            [1,0,0,2,1],
            [3,0,0,0,0],
            [3,0,2,-2,2],
            [0,0,4,-2,2],
            [0,1,2,0,1],
            [0,0,-2,2,1],
            [0,0,2,-2,3],
            [-1,0,0,4,0],
            [2,0,-2,0,1],
            [-2,0,0,4,0],
            [-1,-1,0,2,1],
            [-1,0,0,1,1],
            [0,1,0,0,2],
            [0,0,-2,0,1],
            [0,-1,2,0,1],
            [0,0,2,-1,2],
            [0,0,2,4,2],
            [-2,-1,0,2,0],
            [1,1,0,-2,1],
            [-1,1,0,2,0],
            [-1,1,0,1,2],
            [1,-1,0,0,1],
            [1,-1,2,2,2],
            [-1,1,2,2,2],
            [3,0,2,0,1],
            [0,1,-2,2,0],
            [-1,0,0,-2,1],
            [0,1,2,2,2],
            [-1,-1,2,2,1],
            [0,-1,0,0,2],
            [1,0,2,-4,1],
            [-1,0,-2,2,0],
            [0,-1,2,2,1],
            [2,-1,2,0,2],
            [0,0,0,2,2],
            [1,-1,2,0,1],
            [-1,1,2,0,2],
            [0,1,0,2,0],
            [0,-1,-2,2,0],
            [0,3,2,-2,2],
            [0,0,0,1,1],
            [-1,0,2,2,0],
            [2,1,2,0,2],
            [1,1,0,0,1],
            [1,1,2,0,1],
            [2,0,0,2,0],
            [1,0,-2,2,0],
            [-1,0,0,2,2],
            [0,1,0,1,0],
            [0,1,0,-2,1],
            [-1,0,2,-2,2],
            [0,0,0,-1,1],
            [-1,1,0,0,1],
            [1,0,2,-1,2],
            [1,-1,0,2,0],
            [0,0,0,4,0],
            [1,0,2,1,2],
            [0,0,2,1,1],
            [1,0,0,-2,2],
            [-1,0,2,4,1],
            [1,0,-2,0,1],
            [1,1,2,-2,1],
            [0,0,2,2,0],
            [-1,0,2,-1,1],
            [-2,0,2,2,1],
            [4,0,2,0,2],
            [2,-1,0,0,0],
            [2,1,2,-2,2],
            [0,1,2,1,2],
            [1,0,4,-2,2],
            [-1,-1,0,0,1],
            [0,1,0,2,1],
            [-2,0,2,4,1],
            [2,0,2,0,0],
            [1,0,0,1,0],
            [-1,0,0,4,1],
            [-1,0,4,0,1],
            [2,0,2,2,1],
            [0,0,2,-3,2],
            [-1,-2,0,2,0],
            [2,1,0,0,0],
            [0,0,4,0,2],
            [0,0,0,0,3],
            [0,3,0,0,0],
            [0,0,2,-4,1],
            [0,-1,0,2,1],
            [0,0,0,4,1],
            [-1,-1,2,4,2],
            [1,0,2,4,2],
            [-2,2,0,2,0],
            [-2,-1,2,0,1],
            [-2,0,0,2,2],
            [-1,-1,2,0,2],
            [0,0,4,-2,1],
            [3,0,2,-2,1],
            [-2,-1,0,2,1],
            [1,0,0,-1,1],
            [0,-2,0,2,0],
            [-2,0,0,4,1],
            [-3,0,0,0,1],
            [1,1,2,2,2],
            [0,0,2,4,1],
            [3,0,2,2,2],
            [-1,1,2,-2,1],
            [2,0,0,-4,1],
            [0,0,0,-2,2],
            [2,0,2,-4,1],
            [-1,1,0,2,1],
            [0,0,2,-1,1],
            [0,-2,2,2,2],
            [2,0,0,2,1],
            [4,0,2,-2,2],
            [2,0,0,-2,2],
            [0,2,0,0,1],
            [1,0,0,-4,1],
            [0,2,2,-2,1],
            [-3,0,0,4,0],
            [-1,1,2,0,1],
            [-1,-1,0,4,0],
            [-1,-2,2,2,2],
            [-2,-1,2,4,2],
            [1,-1,2,2,1],
            [-2,1,0,2,0],
            [-2,1,2,0,1],
            [2,1,0,-2,1],
            [-3,0,2,0,1],
            [-2,0,2,-2,1],
            [-1,1,0,2,2],
            [0,-1,2,-1,2],
            [-1,0,4,-2,2],
            [0,-2,2,0,2],
            [-1,0,2,1,2],
            [2,0,0,0,2],
            [0,0,2,0,3],
            [-2,0,4,0,2],
            [-1,0,-2,0,1],
            [-1,1,2,2,1],
            [3,0,0,0,1],
            [-1,0,2,3,2],
            [2,-1,2,0,1],
            [0,1,2,2,1],
            [0,-1,2,4,2],
            [2,-1,2,2,2],
            [0,2,-2,2,0],
            [-1,-1,2,-1,1],
            [0,-2,0,0,1],
            [1,0,2,-4,2],
            [1,-1,0,-2,1],
            [-1,-1,2,0,1],
            [1,-1,2,-2,2],
            [-2,-1,0,4,0],
            [-1,0,0,3,0],
            [-2,-1,2,2,2],
            [0,2,2,0,2],
            [1,1,0,2,0],
            [2,0,2,-1,2],
            [1,0,2,1,1],
            [4,0,0,0,0],
            [2,1,2,0,1],
            [3,-1,2,0,2],
            [-2,2,0,2,1],
            [1,0,2,-3,1],
            [1,1,2,-4,1],
            [-1,-1,2,-2,1],
            [0,-1,0,-1,1],
            [0,-1,0,-2,1],
            [-2,0,0,0,2],
            [-2,0,-2,2,0],
            [-1,0,-2,4,0],
            [1,-2,0,0,0],
            [0,1,0,1,1],
            [-1,2,0,2,0],
            [1,-1,2,-2,1],
            [1,2,2,-2,2],
            [2,-1,2,-2,2],
            [1,0,2,-1,1],
            [2,1,2,-2,1],
            [-2,0,0,-2,1],
            [1,-2,2,0,2],
            [0,1,2,1,1],
            [1,0,4,-2,1],
            [-2,0,4,2,2],
            [1,1,2,1,2],
            [1,0,0,4,0],
            [1,0,2,2,0],
            [2,0,2,1,2],
            [3,1,2,0,2],
            [4,0,2,0,1],
            [-2,-1,2,0,0],
            [0,1,-2,2,1],
            [1,0,-2,1,0],
            [0,-1,-2,2,1],
            [2,-1,0,-2,1],
            [-1,0,2,-1,2],
            [1,0,2,-3,2],
            [0,1,2,-2,3],
            [0,0,2,-3,1],
            [-1,0,-2,2,1],
            [0,0,2,-4,2],
            [-2,1,0,0,1],
            [-1,0,0,-1,1],
            [2,0,2,-4,2],
            [0,0,4,-4,4],
            [0,0,4,-4,2],
            [-1,-2,0,2,1],
            [-2,0,0,3,0],
            [1,0,-2,2,1],
            [-3,0,2,2,2],
            [-3,0,2,2,1],
            [-2,0,2,2,0],
            [2,-1,0,0,1],
            [-2,1,2,2,2],
            [1,1,0,1,0],
            [0,1,4,-2,2],
            [-1,1,0,-2,1],
            [0,0,0,-4,1],
            [1,-1,0,2,1],
            [1,1,0,2,1],
            [-1,2,2,2,2],
            [3,1,2,-2,2],
            [0,-1,0,4,0],
            [2,-1,0,2,0],
            [0,0,4,0,1],
            [2,0,4,-2,2],
            [-1,-1,2,4,1],
            [1,0,0,4,1],
            [1,-2,2,2,2],
            [0,0,2,3,2],
            [-1,1,2,4,2],
            [3,0,0,2,0],
            [-1,0,4,2,2],
            [1,1,2,2,1],
            [-2,0,2,6,2],
            [2,1,2,2,2],
            [-1,0,2,6,2],
            [1,0,2,4,1],
            [2,0,2,4,2],
            [1,1,-2,1,0],
            [-3,1,2,1,2],
            [2,0,-2,0,2],
            [-1,0,0,1,2],
            [-4,0,2,2,1],
            [-1,-1,0,1,0],
            [0,0,-2,2,2],
            [1,0,0,-1,2],
            [0,-1,2,-2,3],
            [-2,1,2,0,0],
            [0,0,2,-2,4],
            [-2,-2,0,2,0],
            [-2,0,-2,4,0],
            [0,-2,-2,2,0],
            [1,2,0,-2,1],
            [3,0,0,-4,1],
            [-1,1,2,-2,2],
            [1,-1,2,-4,1],
            [1,1,0,-2,2],
            [-3,0,2,0,0],
            [-3,0,2,0,2],
            [-2,0,0,1,0],
            [0,0,-2,1,0],
            [-3,0,0,2,1],
            [-1,-1,-2,2,0],
            [0,1,2,-4,1],
            [2,1,0,-4,1],
            [0,2,0,-2,1],
            [1,0,0,-3,1],
            [-2,0,2,-2,2],
            [-2,-1,0,0,1],
            [-4,0,0,2,0],
            [1,1,0,-4,1],
            [-1,0,2,-4,1],
            [0,0,4,-4,1],
            [0,3,2,-2,2],
            [-3,-1,0,4,0],
            [-3,0,0,4,1],
            [1,-1,-2,2,0],
            [-1,-1,0,2,2],
            [1,-2,0,0,1],
            [1,-1,0,0,2],
            [0,0,0,1,2],
            [-1,-1,2,0,0],
            [1,-2,2,-2,2],
            [0,-1,2,-1,1],
            [-1,0,2,0,3],
            [1,1,0,0,2],
            [-1,1,2,0,0],
            [1,2,0,0,0],
            [-1,2,2,0,2],
            [-1,0,4,-2,1],
            [3,0,2,-4,2],
            [1,2,2,-2,1],
            [1,0,4,-4,2],
            [-2,-1,0,4,1],
            [0,-1,0,2,2],
            [-2,1,0,4,0],
            [-2,-1,2,2,1],
            [2,0,-2,2,0],
            [1,0,0,1,1],
            [0,1,0,2,2],
            [1,-1,2,-1,2],
            [-2,0,4,0,1],
            [2,1,0,0,1],
            [0,1,2,0,0],
            [0,-1,4,-2,2],
            [0,0,4,-2,4],
            [0,2,2,0,1],
            [-3,0,0,6,0],
            [-1,-1,0,4,1],
            [1,-2,0,2,0],
            [-1,0,0,4,2],
            [-1,-2,2,2,1],
            [-1,0,0,-2,2],
            [1,0,-2,-2,1],
            [0,0,-2,-2,1],
            [-2,0,-2,0,1],
            [0,0,0,3,1],
            [0,0,0,3,0],
            [-1,1,0,4,0],
            [-1,-1,2,2,0],
            [-2,0,2,3,2],
            [1,0,0,2,2],
            [0,-1,2,1,2],
            [3,-1,0,0,0],
            [2,0,0,1,0],
            [1,-1,2,0,0],
            [0,0,2,1,0],
            [1,0,2,0,3],
            [3,1,0,0,0],
            [3,-1,2,-2,2],
            [2,0,2,-1,1],
            [1,1,2,0,0],
            [0,0,4,-1,2],
            [1,2,2,0,2],
            [-2,0,0,6,0],
            [0,-1,0,4,1],
            [-2,-1,2,4,1],
            [0,-2,2,2,1],
            [0,-1,2,2,0],
            [-1,0,2,3,1],
            [-2,1,2,4,2],
            [2,0,0,2,2],
            [2,-2,2,0,2],
            [-1,1,2,3,2],
            [3,0,2,-1,2],
            [4,0,2,-2,1],
            [-1,0,0,6,0],
            [-1,-2,2,4,2],
            [-3,0,2,6,2],
            [-1,0,2,4,0],
            [3,0,0,2,1],
            [3,-1,2,0,1],
            [3,0,2,0,0],
            [1,0,4,0,2],
            [5,0,2,-2,2],
            [0,-1,2,4,1],
            [2,-1,2,2,1],
            [0,1,2,4,2],
            [1,-1,2,4,2],
            [3,-1,2,2,2],
            [3,0,2,2,1],
            [5,0,2,0,2],
            [0,0,2,6,2],
            [4,0,2,2,2],
            [0,-1,1,-1,1],
            [-1,0,1,0,3],
            [0,-2,2,-2,3],
            [1,0,-1,0,1],
            [2,-2,0,-2,1],
            [-1,0,1,0,2],
            [-1,0,1,0,1],
            [-1,-1,2,-1,2],
            [-2,2,0,2,2],
            [-1,0,1,0,0],
            [-4,1,2,2,2],
            [-3,0,2,1,1],
            [-2,-1,2,0,2],
            [1,0,-2,1,1],
            [2,-1,-2,0,1],
            [-4,0,2,2,0],
            [-3,1,0,3,0],
            [-1,0,-1,2,0],
            [0,-2,0,0,2],
            [0,-2,0,0,2],
            [-3,0,0,3,0],
            [-2,-1,0,2,2],
            [-1,0,-2,3,0],
            [-4,0,0,4,0],
            [2,1,-2,0,1],
            [2,-1,0,-2,2],
            [0,0,1,-1,0],
            [-1,2,0,1,0],
            [-2,1,2,0,2],
            [1,1,0,-1,1],
            [1,0,1,-2,1],
            [0,2,0,0,2],
            [1,-1,2,-3,1],
            [-1,1,2,-1,1],
            [-2,0,4,-2,2],
            [-2,0,4,-2,1],
            [-2,-2,0,2,1],
            [-2,0,-2,4,0],
            [1,2,2,-4,1],
            [1,1,2,-4,2],
            [-1,2,2,-2,1],
            [2,0,0,-3,1],
            [-1,2,0,0,1],
            [0,0,0,-2,0],
            [-1,-1,2,-2,2],
            [-1,1,0,0,2],
            [0,0,0,-1,2],
            [-2,1,0,1,0],
            [1,-2,0,-2,1],
            [1,0,-2,0,2],
            [-3,1,0,2,0],
            [-1,1,-2,2,0],
            [-1,-1,0,0,2],
            [-3,0,0,2,0],
            [-3,-1,0,2,0],
            [2,0,2,-6,1],
            [0,1,2,-4,2],
            [2,0,0,-4,2],
            [-2,1,2,-2,1],
            [0,-1,2,-4,1],
            [0,1,0,-2,2],
            [-1,0,0,-2,0],
            [2,0,-2,-2,1],
            [-4,0,2,0,1],
            [-1,-1,0,-1,1],
            [0,0,-2,0,2],
            [-3,0,0,1,0],
            [-1,0,-2,1,0],
            [-2,0,-2,2,1],
            [0,0,-4,2,0],
            [-2,-1,-2,2,0],
            [1,0,2,-6,1],
            [-1,0,2,-4,2],
            [1,0,0,-4,2],
            [2,1,2,-4,2],
            [2,1,2,-4,1],
            [0,1,4,-4,4],
            [0,1,4,-4,2],
            [-1,-1,-2,4,0],
            [-1,-3,0,2,0],
            [-1,0,-2,4,1],
            [-2,-1,0,3,0],
            [0,0,-2,3,0],
            [-2,0,0,3,1],
            [0,-1,0,1,0],
            [-3,0,2,2,0],
            [1,1,-2,2,0],
            [-1,1,0,2,2],
            [1,-2,2,-2,1],
            [0,0,1,0,2],
            [0,0,1,0,1],
            [0,0,1,0,0],
            [-1,2,0,2,1],
            [0,0,2,0,2],
            [-2,0,2,0,2],
            [2,0,0,-1,1],
            [3,0,0,-2,1],
            [1,0,2,-2,3],
            [1,2,0,0,1],
            [2,0,2,-3,2],
            [-1,1,4,-2,2],
            [-2,-2,0,4,0],
            [0,-3,0,2,0],
            [0,0,-2,4,0],
            [-1,-1,0,3,0],
            [-2,0,0,4,2],
            [-1,0,0,3,1],
            [2,-2,0,0,0],
            [1,-1,0,1,0],
            [-1,0,0,2,0],
            [0,-2,2,0,1],
            [-1,0,1,2,1],
            [-1,1,0,3,0],
            [-1,-1,2,1,2],
            [0,-1,2,0,0],
            [-2,1,2,2,1],
            [2,-2,2,-2,2],
            [1,1,0,1,1],
            [1,0,1,0,1],
            [1,0,1,0,0],
            [0,2,0,2,0],
            [2,-1,2,-2,1],
            [0,-1,4,-2,1],
            [0,0,4,-2,3],
            [0,1,4,-2,1],
            [4,0,2,-4,2],
            [2,2,2,-2,2],
            [2,0,4,-4,2],
            [-1,-2,0,4,0],
            [-1,-3,2,2,2],
            [-3,0,2,4,2],
            [-3,0,2,-2,1],
            [-1,-1,0,-2,1],
            [-3,0,0,0,2],
            [-3,0,-2,2,0],
            [0,1,0,-4,1],
            [-2,1,0,-2,1],
            [-4,0,0,0,1],
            [-1,0,0,-4,1],
            [-3,0,0,-2,1],
            [0,0,0,3,2],
            [-1,1,0,4,1],
            [1,-2,2,0,1],
            [0,1,0,3,0],
            [-1,0,2,2,3],
            [0,0,2,2,2],
            [-2,0,2,2,2],
            [-1,1,2,2,0],
            [3,0,0,0,2],
            [2,1,0,1,0],
            [2,-1,2,-1,2],
            [0,0,2,0,1],
            [0,0,3,0,3],
            [0,0,3,0,2],
            [-1,2,2,2,1],
            [-1,0,4,0,0],
            [1,2,2,0,1],
            [3,1,2,-2,1],
            [1,1,4,-2,2],
            [-2,-1,0,6,0],
            [0,-2,0,4,0],
            [-2,0,0,6,1],
            [-2,-2,2,4,2],
            [0,-3,2,2,2],
            [0,0,0,4,2],
            [-1,-1,2,3,2],
            [-2,0,2,4,0],
            [2,-1,0,2,1],
            [1,0,0,3,0],
            [0,1,0,4,1],
            [0,1,0,4,0],
            [1,-1,2,1,2],
            [0,0,2,2,3],
            [1,0,2,2,2],
            [-1,0,2,2,2],
            [-2,0,4,2,1],
            [2,1,0,2,1],
            [2,1,0,2,0],
            [2,-1,2,0,0],
            [1,0,2,1,0],
            [0,1,2,2,0],
            [2,0,2,0,3],
            [3,0,2,0,2],
            [1,0,2,0,2],
            [1,0,3,0,3],
            [1,1,2,1,1],
            [0,2,2,2,2],
            [2,1,2,0,0],
            [2,0,4,-2,1],
            [4,1,2,-2,2],
            [-1,-1,0,6,0],
            [-3,-1,2,6,2],
            [-1,0,0,6,1],
            [-3,0,2,6,1],
            [1,-1,0,4,1],
            [1,-1,0,4,0],
            [-2,0,2,5,2],
            [1,-2,2,2,1],
            [3,-1,0,2,0],
            [1,-1,2,2,0],
            [0,0,2,3,1],
            [-1,1,2,4,1],
            [0,1,2,3,2],
            [-1,0,4,2,1],
            [2,0,2,1,1],
            [5,0,0,0,0],
            [2,1,2,1,2],
            [1,0,4,0,1],
            [3,1,2,0,1],
            [3,0,4,-2,2],
            [-2,-1,2,6,2],
            [0,0,0,6,0],
            [0,-2,2,4,2],
            [-2,0,2,6,1],
            [2,0,0,4,1],
            [2,0,0,4,0],
            [2,-2,2,2,2],
            [0,0,2,4,0],
            [1,0,2,3,2],
            [4,0,0,2,0],
            [2,0,2,2,0],
            [0,0,4,2,2],
            [4,-1,2,0,2],
            [3,0,2,1,2],
            [2,1,2,2,1],
            [4,1,2,0,2],
            [-1,-1,2,6,2],
            [-1,0,2,6,1],
            [1,-1,2,4,1],
            [1,1,2,4,2],
            [3,1,2,2,2],
            [5,0,2,0,1],
            [2,-1,2,4,2],
            [2,0,2,4,1]]
    
    #日月章动系数，单位：0.1微角秒
    #*  longitude (sin, t*sin, cos), obliquity (cos, t*cos, sin)
    CLS=[[-172064161.0,-174666.0,33386.0,92052331.0,9086.0,15377.0],
            [-13170906.0,-1675.0,-13696.0,5730336.0,-3015.0,-4587.0],
            [-2276413.0,-234.0,2796.0,978459.0,-485.0,1374.0],
            [2074554.0,207.0,-698.0,-897492.0,470.0,-291.0],
            [1475877.0,-3633.0,11817.0,73871.0,-184.0,-1924.0],
            [-516821.0,1226.0,-524.0,224386.0,-677.0,-174.0],
            [711159.0,73.0,-872.0,-6750.0,0.0,358.0],
            [-387298.0,-367.0,380.0,200728.0,18.0,318.0],
            [-301461.0,-36.0,816.0,129025.0,-63.0,367.0],
            [215829.0,-494.0,111.0,-95929.0,299.0,132.0],
            [128227.0,137.0,181.0,-68982.0,-9.0,39.0],
            [123457.0,11.0,19.0,-53311.0,32.0,-4.0],
            [156994.0,10.0,-168.0,-1235.0,0.0,82.0],
            [63110.0,63.0,27.0,-33228.0,0.0,-9.0],
            [-57976.0,-63.0,-189.0,31429.0,0.0,-75.0],
            [-59641.0,-11.0,149.0,25543.0,-11.0,66.0],
            [-51613.0,-42.0,129.0,26366.0,0.0,78.0],
            [45893.0,50.0,31.0,-24236.0,-10.0,20.0],
            [63384.0,11.0,-150.0,-1220.0,0.0,29.0],
            [-38571.0,-1.0,158.0,16452.0,-11.0,68.0],
            [32481.0,0.0,0.0,-13870.0,0.0,0.0],
            [-47722.0,0.0,-18.0,477.0,0.0,-25.0],
            [-31046.0,-1.0,131.0,13238.0,-11.0,59.0],
            [28593.0,0.0,-1.0,-12338.0,10.0,-3.0],
            [20441.0,21.0,10.0,-10758.0,0.0,-3.0],
            [29243.0,0.0,-74.0,-609.0,0.0,13.0],
            [25887.0,0.0,-66.0,-550.0,0.0,11.0],
            [-14053.0,-25.0,79.0,8551.0,-2.0,-45.0],
            [15164.0,10.0,11.0,-8001.0,0.0,-1.0],
            [-15794.0,72.0,-16.0,6850.0,-42.0,-5.0],
            [21783.0,0.0,13.0,-167.0,0.0,13.0],
            [-12873.0,-10.0,-37.0,6953.0,0.0,-14.0],
            [-12654.0,11.0,63.0,6415.0,0.0,26.0],
            [-10204.0,0.0,25.0,5222.0,0.0,15.0],
            [16707.0,-85.0,-10.0,168.0,-1.0,10.0],
            [-7691.0,0.0,44.0,3268.0,0.0,19.0],
            [-11024.0,0.0,-14.0,104.0,0.0,2.0],
            [7566.0,-21.0,-11.0,-3250.0,0.0,-5.0],
            [-6637.0,-11.0,25.0,3353.0,0.0,14.0],
            [-7141.0,21.0,8.0,3070.0,0.0,4.0],
            [-6302.0,-11.0,2.0,3272.0,0.0,4.0],
            [5800.0,10.0,2.0,-3045.0,0.0,-1.0],
            [6443.0,0.0,-7.0,-2768.0,0.0,-4.0],
            [-5774.0,-11.0,-15.0,3041.0,0.0,-5.0],
            [-5350.0,0.0,21.0,2695.0,0.0,12.0],
            [-4752.0,-11.0,-3.0,2719.0,0.0,-3.0],
            [-4940.0,-11.0,-21.0,2720.0,0.0,-9.0],
            [7350.0,0.0,-8.0,-51.0,0.0,4.0],
            [4065.0,0.0,6.0,-2206.0,0.0,1.0],
            [6579.0,0.0,-24.0,-199.0,0.0,2.0],
            [3579.0,0.0,5.0,-1900.0,0.0,1.0],
            [4725.0,0.0,-6.0,-41.0,0.0,3.0],
            [-3075.0,0.0,-2.0,1313.0,0.0,-1.0],
            [-2904.0,0.0,15.0,1233.0,0.0,7.0],
            [4348.0,0.0,-10.0,-81.0,0.0,2.0],
            [-2878.0,0.0,8.0,1232.0,0.0,4.0],
            [-4230.0,0.0,5.0,-20.0,0.0,-2.0],
            [-2819.0,0.0,7.0,1207.0,0.0,3.0],
            [-4056.0,0.0,5.0,40.0,0.0,-2.0],
            [-2647.0,0.0,11.0,1129.0,0.0,5.0],
            [-2294.0,0.0,-10.0,1266.0,0.0,-4.0],
            [2481.0,0.0,-7.0,-1062.0,0.0,-3.0],
            [2179.0,0.0,-2.0,-1129.0,0.0,-2.0],
            [3276.0,0.0,1.0,-9.0,0.0,0.0],
            [-3389.0,0.0,5.0,35.0,0.0,-2.0],
            [3339.0,0.0,-13.0,-107.0,0.0,1.0],
            [-1987.0,0.0,-6.0,1073.0,0.0,-2.0],
            [-1981.0,0.0,0.0,854.0,0.0,0.0],
            [4026.0,0.0,-353.0,-553.0,0.0,-139.0],
            [1660.0,0.0,-5.0,-710.0,0.0,-2.0],
            [-1521.0,0.0,9.0,647.0,0.0,4.0],
            [1314.0,0.0,0.0,-700.0,0.0,0.0],
            [-1283.0,0.0,0.0,672.0,0.0,0.0],
            [-1331.0,0.0,8.0,663.0,0.0,4.0],
            [1383.0,0.0,-2.0,-594.0,0.0,-2.0],
            [1405.0,0.0,4.0,-610.0,0.0,2.0],
            [1290.0,0.0,0.0,-556.0,0.0,0.0],
            [-1214.0,0.0,5.0,518.0,0.0,2.0],
            [1146.0,0.0,-3.0,-490.0,0.0,-1.0],
            [1019.0,0.0,-1.0,-527.0,0.0,-1.0],
            [-1100.0,0.0,9.0,465.0,0.0,4.0],
            [-970.0,0.0,2.0,496.0,0.0,1.0],
            [1575.0,0.0,-6.0,-50.0,0.0,0.0],
            [934.0,0.0,-3.0,-399.0,0.0,-1.0],
            [922.0,0.0,-1.0,-395.0,0.0,-1.0],
            [815.0,0.0,-1.0,-422.0,0.0,-1.0],
            [834.0,0.0,2.0,-440.0,0.0,1.0],
            [1248.0,0.0,0.0,-170.0,0.0,1.0],
            [1338.0,0.0,-5.0,-39.0,0.0,0.0],
            [716.0,0.0,-2.0,-389.0,0.0,-1.0],
            [1282.0,0.0,-3.0,-23.0,0.0,1.0],
            [742.0,0.0,1.0,-391.0,0.0,0.0],
            [1020.0,0.0,-25.0,-495.0,0.0,-10.0],
            [715.0,0.0,-4.0,-326.0,0.0,2.0],
            [-666.0,0.0,-3.0,369.0,0.0,-1.0],
            [-667.0,0.0,1.0,346.0,0.0,1.0],
            [-704.0,0.0,0.0,304.0,0.0,0.0],
            [-694.0,0.0,5.0,294.0,0.0,2.0],
            [-1014.0,0.0,-1.0,4.0,0.0,-1.0],
            [-585.0,0.0,-2.0,316.0,0.0,-1.0],
            [-949.0,0.0,1.0,8.0,0.0,-1.0],
            [-595.0,0.0,0.0,258.0,0.0,0.0],
            [528.0,0.0,0.0,-279.0,0.0,0.0],
            [-590.0,0.0,4.0,252.0,0.0,2.0],
            [570.0,0.0,-2.0,-244.0,0.0,-1.0],
            [-502.0,0.0,3.0,250.0,0.0,2.0],
            [-875.0,0.0,1.0,29.0,0.0,0.0],
            [-492.0,0.0,-3.0,275.0,0.0,-1.0],
            [535.0,0.0,-2.0,-228.0,0.0,-1.0],
            [-467.0,0.0,1.0,240.0,0.0,1.0],
            [591.0,0.0,0.0,-253.0,0.0,0.0],
            [-453.0,0.0,-1.0,244.0,0.0,-1.0],
            [766.0,0.0,1.0,9.0,0.0,0.0],
            [-446.0,0.0,2.0,225.0,0.0,1.0],
            [-488.0,0.0,2.0,207.0,0.0,1.0],
            [-468.0,0.0,0.0,201.0,0.0,0.0],
            [-421.0,0.0,1.0,216.0,0.0,1.0],
            [463.0,0.0,0.0,-200.0,0.0,0.0],
            [-673.0,0.0,2.0,14.0,0.0,0.0],
            [658.0,0.0,0.0,-2.0,0.0,0.0],
            [-438.0,0.0,0.0,188.0,0.0,0.0],
            [-390.0,0.0,0.0,205.0,0.0,0.0],
            [639.0,-11.0,-2.0,-19.0,0.0,0.0],
            [412.0,0.0,-2.0,-176.0,0.0,-1.0],
            [-361.0,0.0,0.0,189.0,0.0,0.0],
            [360.0,0.0,-1.0,-185.0,0.0,-1.0],
            [588.0,0.0,-3.0,-24.0,0.0,0.0],
            [-578.0,0.0,1.0,5.0,0.0,0.0],
            [-396.0,0.0,0.0,171.0,0.0,0.0],
            [565.0,0.0,-1.0,-6.0,0.0,0.0],
            [-335.0,0.0,-1.0,184.0,0.0,-1.0],
            [357.0,0.0,1.0,-154.0,0.0,0.0],
            [321.0,0.0,1.0,-174.0,0.0,0.0],
            [-301.0,0.0,-1.0,162.0,0.0,0.0],
            [-334.0,0.0,0.0,144.0,0.0,0.0],
            [493.0,0.0,-2.0,-15.0,0.0,0.0],
            [494.0,0.0,-2.0,-19.0,0.0,0.0],
            [337.0,0.0,-1.0,-143.0,0.0,-1.0],
            [280.0,0.0,-1.0,-144.0,0.0,0.0],
            [309.0,0.0,1.0,-134.0,0.0,0.0],
            [-263.0,0.0,2.0,131.0,0.0,1.0],
            [253.0,0.0,1.0,-138.0,0.0,0.0],
            [245.0,0.0,0.0,-128.0,0.0,0.0],
            [416.0,0.0,-2.0,-17.0,0.0,0.0],
            [-229.0,0.0,0.0,128.0,0.0,0.0],
            [231.0,0.0,0.0,-120.0,0.0,0.0],
            [-259.0,0.0,2.0,109.0,0.0,1.0],
            [375.0,0.0,-1.0,-8.0,0.0,0.0],
            [252.0,0.0,0.0,-108.0,0.0,0.0],
            [-245.0,0.0,1.0,104.0,0.0,0.0],
            [243.0,0.0,-1.0,-104.0,0.0,0.0],
            [208.0,0.0,1.0,-112.0,0.0,0.0],
            [199.0,0.0,0.0,-102.0,0.0,0.0],
            [-208.0,0.0,1.0,105.0,0.0,0.0],
            [335.0,0.0,-2.0,-14.0,0.0,0.0],
            [-325.0,0.0,1.0,7.0,0.0,0.0],
            [-187.0,0.0,0.0,96.0,0.0,0.0],
            [197.0,0.0,-1.0,-100.0,0.0,0.0],
            [-192.0,0.0,2.0,94.0,0.0,1.0],
            [-188.0,0.0,0.0,83.0,0.0,0.0],
            [276.0,0.0,0.0,-2.0,0.0,0.0],
            [-286.0,0.0,1.0,6.0,0.0,0.0],
            [186.0,0.0,-1.0,-79.0,0.0,0.0],
            [-219.0,0.0,0.0,43.0,0.0,0.0],
            [276.0,0.0,0.0,2.0,0.0,0.0],
            [-153.0,0.0,-1.0,84.0,0.0,0.0],
            [-156.0,0.0,0.0,81.0,0.0,0.0],
            [-154.0,0.0,1.0,78.0,0.0,0.0],
            [-174.0,0.0,1.0,75.0,0.0,0.0],
            [-163.0,0.0,2.0,69.0,0.0,1.0],
            [-228.0,0.0,0.0,1.0,0.0,0.0],
            [91.0,0.0,-4.0,-54.0,0.0,-2.0],
            [175.0,0.0,0.0,-75.0,0.0,0.0],
            [-159.0,0.0,0.0,69.0,0.0,0.0],
            [141.0,0.0,0.0,-72.0,0.0,0.0],
            [147.0,0.0,0.0,-75.0,0.0,0.0],
            [-132.0,0.0,0.0,69.0,0.0,0.0],
            [159.0,0.0,-28.0,-54.0,0.0,11.0],
            [213.0,0.0,0.0,-4.0,0.0,0.0],
            [123.0,0.0,0.0,-64.0,0.0,0.0],
            [-118.0,0.0,-1.0,66.0,0.0,0.0],
            [144.0,0.0,-1.0,-61.0,0.0,0.0],
            [-121.0,0.0,1.0,60.0,0.0,0.0],
            [-134.0,0.0,1.0,56.0,0.0,1.0],
            [-105.0,0.0,0.0,57.0,0.0,0.0],
            [-102.0,0.0,0.0,56.0,0.0,0.0],
            [120.0,0.0,0.0,-52.0,0.0,0.0],
            [101.0,0.0,0.0,-54.0,0.0,0.0],
            [-113.0,0.0,0.0,59.0,0.0,0.0],
            [-106.0,0.0,0.0,61.0,0.0,0.0],
            [-129.0,0.0,1.0,55.0,0.0,0.0],
            [-114.0,0.0,0.0,57.0,0.0,0.0],
            [113.0,0.0,-1.0,-49.0,0.0,0.0],
            [-102.0,0.0,0.0,44.0,0.0,0.0],
            [-94.0,0.0,0.0,51.0,0.0,0.0],
            [-100.0,0.0,-1.0,56.0,0.0,0.0],
            [87.0,0.0,0.0,-47.0,0.0,0.0],
            [161.0,0.0,0.0,-1.0,0.0,0.0],
            [96.0,0.0,0.0,-50.0,0.0,0.0],
            [151.0,0.0,-1.0,-5.0,0.0,0.0],
            [-104.0,0.0,0.0,44.0,0.0,0.0],
            [-110.0,0.0,0.0,48.0,0.0,0.0],
            [-100.0,0.0,1.0,50.0,0.0,0.0],
            [92.0,0.0,-5.0,12.0,0.0,-2.0],
            [82.0,0.0,0.0,-45.0,0.0,0.0],
            [82.0,0.0,0.0,-45.0,0.0,0.0],
            [-78.0,0.0,0.0,41.0,0.0,0.0],
            [-77.0,0.0,0.0,43.0,0.0,0.0],
            [2.0,0.0,0.0,54.0,0.0,0.0],
            [94.0,0.0,0.0,-40.0,0.0,0.0],
            [-93.0,0.0,0.0,40.0,0.0,0.0],
            [-83.0,0.0,10.0,40.0,0.0,-2.0],
            [83.0,0.0,0.0,-36.0,0.0,0.0],
            [-91.0,0.0,0.0,39.0,0.0,0.0],
            [128.0,0.0,0.0,-1.0,0.0,0.0],
            [-79.0,0.0,0.0,34.0,0.0,0.0],
            [-83.0,0.0,0.0,47.0,0.0,0.0],
            [84.0,0.0,0.0,-44.0,0.0,0.0],
            [83.0,0.0,0.0,-43.0,0.0,0.0],
            [91.0,0.0,0.0,-39.0,0.0,0.0],
            [-77.0,0.0,0.0,39.0,0.0,0.0],
            [84.0,0.0,0.0,-43.0,0.0,0.0],
            [-92.0,0.0,1.0,39.0,0.0,0.0],
            [-92.0,0.0,1.0,39.0,0.0,0.0],
            [-94.0,0.0,0.0,0.0,0.0,0.0],
            [68.0,0.0,0.0,-36.0,0.0,0.0],
            [-61.0,0.0,0.0,32.0,0.0,0.0],
            [71.0,0.0,0.0,-31.0,0.0,0.0],
            [62.0,0.0,0.0,-34.0,0.0,0.0],
            [-63.0,0.0,0.0,33.0,0.0,0.0],
            [-73.0,0.0,0.0,32.0,0.0,0.0],
            [115.0,0.0,0.0,-2.0,0.0,0.0],
            [-103.0,0.0,0.0,2.0,0.0,0.0],
            [63.0,0.0,0.0,-28.0,0.0,0.0],
            [74.0,0.0,0.0,-32.0,0.0,0.0],
            [-103.0,0.0,-3.0,3.0,0.0,-1.0],
            [-69.0,0.0,0.0,30.0,0.0,0.0],
            [57.0,0.0,0.0,-29.0,0.0,0.0],
            [94.0,0.0,0.0,-4.0,0.0,0.0],
            [64.0,0.0,0.0,-33.0,0.0,0.0],
            [-63.0,0.0,0.0,26.0,0.0,0.0],
            [-38.0,0.0,0.0,20.0,0.0,0.0],
            [-43.0,0.0,0.0,24.0,0.0,0.0],
            [-45.0,0.0,0.0,23.0,0.0,0.0],
            [47.0,0.0,0.0,-24.0,0.0,0.0],
            [-48.0,0.0,0.0,25.0,0.0,0.0],
            [45.0,0.0,0.0,-26.0,0.0,0.0],
            [56.0,0.0,0.0,-25.0,0.0,0.0],
            [88.0,0.0,0.0,2.0,0.0,0.0],
            [-75.0,0.0,0.0,0.0,0.0,0.0],
            [85.0,0.0,0.0,0.0,0.0,0.0],
            [49.0,0.0,0.0,-26.0,0.0,0.0],
            [-74.0,0.0,-3.0,-1.0,0.0,-1.0],
            [-39.0,0.0,0.0,21.0,0.0,0.0],
            [45.0,0.0,0.0,-20.0,0.0,0.0],
            [51.0,0.0,0.0,-22.0,0.0,0.0],
            [-40.0,0.0,0.0,21.0,0.0,0.0],
            [41.0,0.0,0.0,-21.0,0.0,0.0],
            [-42.0,0.0,0.0,24.0,0.0,0.0],
            [-51.0,0.0,0.0,22.0,0.0,0.0],
            [-42.0,0.0,0.0,22.0,0.0,0.0],
            [39.0,0.0,0.0,-21.0,0.0,0.0],
            [46.0,0.0,0.0,-18.0,0.0,0.0],
            [-53.0,0.0,0.0,22.0,0.0,0.0],
            [82.0,0.0,0.0,-4.0,0.0,0.0],
            [81.0,0.0,-1.0,-4.0,0.0,0.0],
            [47.0,0.0,0.0,-19.0,0.0,0.0],
            [53.0,0.0,0.0,-23.0,0.0,0.0],
            [-45.0,0.0,0.0,22.0,0.0,0.0],
            [-44.0,0.0,0.0,-2.0,0.0,0.0],
            [-33.0,0.0,0.0,16.0,0.0,0.0],
            [-61.0,0.0,0.0,1.0,0.0,0.0],
            [28.0,0.0,0.0,-15.0,0.0,0.0],
            [-38.0,0.0,0.0,19.0,0.0,0.0],
            [-33.0,0.0,0.0,21.0,0.0,0.0],
            [-60.0,0.0,0.0,0.0,0.0,0.0],
            [48.0,0.0,0.0,-10.0,0.0,0.0],
            [27.0,0.0,0.0,-14.0,0.0,0.0],
            [38.0,0.0,0.0,-20.0,0.0,0.0],
            [31.0,0.0,0.0,-13.0,0.0,0.0],
            [-29.0,0.0,0.0,15.0,0.0,0.0],
            [28.0,0.0,0.0,-15.0,0.0,0.0],
            [-32.0,0.0,0.0,15.0,0.0,0.0],
            [45.0,0.0,0.0,-8.0,0.0,0.0],
            [-44.0,0.0,0.0,19.0,0.0,0.0],
            [28.0,0.0,0.0,-15.0,0.0,0.0],
            [-51.0,0.0,0.0,0.0,0.0,0.0],
            [-36.0,0.0,0.0,20.0,0.0,0.0],
            [44.0,0.0,0.0,-19.0,0.0,0.0],
            [26.0,0.0,0.0,-14.0,0.0,0.0],
            [-60.0,0.0,0.0,2.0,0.0,0.0],
            [35.0,0.0,0.0,-18.0,0.0,0.0],
            [-27.0,0.0,0.0,11.0,0.0,0.0],
            [47.0,0.0,0.0,-1.0,0.0,0.0],
            [36.0,0.0,0.0,-15.0,0.0,0.0],
            [-36.0,0.0,0.0,20.0,0.0,0.0],
            [-35.0,0.0,0.0,19.0,0.0,0.0],
            [-37.0,0.0,0.0,19.0,0.0,0.0],
            [32.0,0.0,0.0,-16.0,0.0,0.0],
            [35.0,0.0,0.0,-14.0,0.0,0.0],
            [32.0,0.0,0.0,-13.0,0.0,0.0],
            [65.0,0.0,0.0,-2.0,0.0,0.0],
            [47.0,0.0,0.0,-1.0,0.0,0.0],
            [32.0,0.0,0.0,-16.0,0.0,0.0],
            [37.0,0.0,0.0,-16.0,0.0,0.0],
            [-30.0,0.0,0.0,15.0,0.0,0.0],
            [-32.0,0.0,0.0,16.0,0.0,0.0],
            [-31.0,0.0,0.0,13.0,0.0,0.0],
            [37.0,0.0,0.0,-16.0,0.0,0.0],
            [31.0,0.0,0.0,-13.0,0.0,0.0],
            [49.0,0.0,0.0,-2.0,0.0,0.0],
            [32.0,0.0,0.0,-13.0,0.0,0.0],
            [23.0,0.0,0.0,-12.0,0.0,0.0],
            [-43.0,0.0,0.0,18.0,0.0,0.0],
            [26.0,0.0,0.0,-11.0,0.0,0.0],
            [-32.0,0.0,0.0,14.0,0.0,0.0],
            [-29.0,0.0,0.0,14.0,0.0,0.0],
            [-27.0,0.0,0.0,12.0,0.0,0.0],
            [30.0,0.0,0.0,0.0,0.0,0.0],
            [-11.0,0.0,0.0,5.0,0.0,0.0],
            [-21.0,0.0,0.0,10.0,0.0,0.0],
            [-34.0,0.0,0.0,15.0,0.0,0.0],
            [-10.0,0.0,0.0,6.0,0.0,0.0],
            [-36.0,0.0,0.0,0.0,0.0,0.0],
            [-9.0,0.0,0.0,4.0,0.0,0.0],
            [-12.0,0.0,0.0,5.0,0.0,0.0],
            [-21.0,0.0,0.0,5.0,0.0,0.0],
            [-29.0,0.0,0.0,-1.0,0.0,0.0],
            [-15.0,0.0,0.0,3.0,0.0,0.0],
            [-20.0,0.0,0.0,0.0,0.0,0.0],
            [28.0,0.0,0.0,0.0,0.0,-2.0],
            [17.0,0.0,0.0,0.0,0.0,0.0],
            [-22.0,0.0,0.0,12.0,0.0,0.0],
            [-14.0,0.0,0.0,7.0,0.0,0.0],
            [24.0,0.0,0.0,-11.0,0.0,0.0],
            [11.0,0.0,0.0,-6.0,0.0,0.0],
            [14.0,0.0,0.0,-6.0,0.0,0.0],
            [24.0,0.0,0.0,0.0,0.0,0.0],
            [18.0,0.0,0.0,-8.0,0.0,0.0],
            [-38.0,0.0,0.0,0.0,0.0,0.0],
            [-31.0,0.0,0.0,0.0,0.0,0.0],
            [-16.0,0.0,0.0,8.0,0.0,0.0],
            [29.0,0.0,0.0,0.0,0.0,0.0],
            [-18.0,0.0,0.0,10.0,0.0,0.0],
            [-10.0,0.0,0.0,5.0,0.0,0.0],
            [-17.0,0.0,0.0,10.0,0.0,0.0],
            [9.0,0.0,0.0,-4.0,0.0,0.0],
            [16.0,0.0,0.0,-6.0,0.0,0.0],
            [22.0,0.0,0.0,-12.0,0.0,0.0],
            [20.0,0.0,0.0,0.0,0.0,0.0],
            [-13.0,0.0,0.0,6.0,0.0,0.0],
            [-17.0,0.0,0.0,9.0,0.0,0.0],
            [-14.0,0.0,0.0,8.0,0.0,0.0],
            [0.0,0.0,0.0,-7.0,0.0,0.0],
            [14.0,0.0,0.0,0.0,0.0,0.0],
            [19.0,0.0,0.0,-10.0,0.0,0.0],
            [-34.0,0.0,0.0,0.0,0.0,0.0],
            [-20.0,0.0,0.0,8.0,0.0,0.0],
            [9.0,0.0,0.0,-5.0,0.0,0.0],
            [-18.0,0.0,0.0,7.0,0.0,0.0],
            [13.0,0.0,0.0,-6.0,0.0,0.0],
            [17.0,0.0,0.0,0.0,0.0,0.0],
            [-12.0,0.0,0.0,5.0,0.0,0.0],
            [15.0,0.0,0.0,-8.0,0.0,0.0],
            [-11.0,0.0,0.0,3.0,0.0,0.0],
            [13.0,0.0,0.0,-5.0,0.0,0.0],
            [-18.0,0.0,0.0,0.0,0.0,0.0],
            [-35.0,0.0,0.0,0.0,0.0,0.0],
            [9.0,0.0,0.0,-4.0,0.0,0.0],
            [-19.0,0.0,0.0,10.0,0.0,0.0],
            [-26.0,0.0,0.0,11.0,0.0,0.0],
            [8.0,0.0,0.0,-4.0,0.0,0.0],
            [-10.0,0.0,0.0,4.0,0.0,0.0],
            [10.0,0.0,0.0,-6.0,0.0,0.0],
            [-21.0,0.0,0.0,9.0,0.0,0.0],
            [-15.0,0.0,0.0,0.0,0.0,0.0],
            [9.0,0.0,0.0,-5.0,0.0,0.0],
            [-29.0,0.0,0.0,0.0,0.0,0.0],
            [-19.0,0.0,0.0,10.0,0.0,0.0],
            [12.0,0.0,0.0,-5.0,0.0,0.0],
            [22.0,0.0,0.0,-9.0,0.0,0.0],
            [-10.0,0.0,0.0,5.0,0.0,0.0],
            [-20.0,0.0,0.0,11.0,0.0,0.0],
            [-20.0,0.0,0.0,0.0,0.0,0.0],
            [-17.0,0.0,0.0,7.0,0.0,0.0],
            [15.0,0.0,0.0,-3.0,0.0,0.0],
            [8.0,0.0,0.0,-4.0,0.0,0.0],
            [14.0,0.0,0.0,0.0,0.0,0.0],
            [-12.0,0.0,0.0,6.0,0.0,0.0],
            [25.0,0.0,0.0,0.0,0.0,0.0],
            [-13.0,0.0,0.0,6.0,0.0,0.0],
            [-14.0,0.0,0.0,8.0,0.0,0.0],
            [13.0,0.0,0.0,-5.0,0.0,0.0],
            [-17.0,0.0,0.0,9.0,0.0,0.0],
            [-12.0,0.0,0.0,6.0,0.0,0.0],
            [-10.0,0.0,0.0,5.0,0.0,0.0],
            [10.0,0.0,0.0,-6.0,0.0,0.0],
            [-15.0,0.0,0.0,0.0,0.0,0.0],
            [-22.0,0.0,0.0,0.0,0.0,0.0],
            [28.0,0.0,0.0,-1.0,0.0,0.0],
            [15.0,0.0,0.0,-7.0,0.0,0.0],
            [23.0,0.0,0.0,-10.0,0.0,0.0],
            [12.0,0.0,0.0,-5.0,0.0,0.0],
            [29.0,0.0,0.0,-1.0,0.0,0.0],
            [-25.0,0.0,0.0,1.0,0.0,0.0],
            [22.0,0.0,0.0,0.0,0.0,0.0],
            [-18.0,0.0,0.0,0.0,0.0,0.0],
            [15.0,0.0,0.0,3.0,0.0,0.0],
            [-23.0,0.0,0.0,0.0,0.0,0.0],
            [12.0,0.0,0.0,-5.0,0.0,0.0],
            [-8.0,0.0,0.0,4.0,0.0,0.0],
            [-19.0,0.0,0.0,0.0,0.0,0.0],
            [-10.0,0.0,0.0,4.0,0.0,0.0],
            [21.0,0.0,0.0,-9.0,0.0,0.0],
            [23.0,0.0,0.0,-1.0,0.0,0.0],
            [-16.0,0.0,0.0,8.0,0.0,0.0],
            [-19.0,0.0,0.0,9.0,0.0,0.0],
            [-22.0,0.0,0.0,10.0,0.0,0.0],
            [27.0,0.0,0.0,-1.0,0.0,0.0],
            [16.0,0.0,0.0,-8.0,0.0,0.0],
            [19.0,0.0,0.0,-8.0,0.0,0.0],
            [9.0,0.0,0.0,-4.0,0.0,0.0],
            [-9.0,0.0,0.0,4.0,0.0,0.0],
            [-9.0,0.0,0.0,4.0,0.0,0.0],
            [-8.0,0.0,0.0,4.0,0.0,0.0],
            [18.0,0.0,0.0,-9.0,0.0,0.0],
            [16.0,0.0,0.0,-1.0,0.0,0.0],
            [-10.0,0.0,0.0,4.0,0.0,0.0],
            [-23.0,0.0,0.0,9.0,0.0,0.0],
            [16.0,0.0,0.0,-1.0,0.0,0.0],
            [-12.0,0.0,0.0,6.0,0.0,0.0],
            [-8.0,0.0,0.0,4.0,0.0,0.0],
            [30.0,0.0,0.0,-2.0,0.0,0.0],
            [24.0,0.0,0.0,-10.0,0.0,0.0],
            [10.0,0.0,0.0,-4.0,0.0,0.0],
            [-16.0,0.0,0.0,7.0,0.0,0.0],
            [-16.0,0.0,0.0,7.0,0.0,0.0],
            [17.0,0.0,0.0,-7.0,0.0,0.0],
            [-24.0,0.0,0.0,10.0,0.0,0.0],
            [-12.0,0.0,0.0,5.0,0.0,0.0],
            [-24.0,0.0,0.0,11.0,0.0,0.0],
            [-23.0,0.0,0.0,9.0,0.0,0.0],
            [-13.0,0.0,0.0,5.0,0.0,0.0],
            [-15.0,0.0,0.0,7.0,0.0,0.0],
            [0.0,0.0,-1988.0,0.0,0.0,-1679.0],
            [0.0,0.0,-63.0,0.0,0.0,-27.0],
            [-4.0,0.0,0.0,0.0,0.0,0.0],
            [0.0,0.0,5.0,0.0,0.0,4.0],
            [5.0,0.0,0.0,-3.0,0.0,0.0],
            [0.0,0.0,364.0,0.0,0.0,176.0],
            [0.0,0.0,-1044.0,0.0,0.0,-891.0],
            [-3.0,0.0,0.0,1.0,0.0,0.0],
            [4.0,0.0,0.0,-2.0,0.0,0.0],
            [0.0,0.0,330.0,0.0,0.0,0.0],
            [5.0,0.0,0.0,-2.0,0.0,0.0],
            [3.0,0.0,0.0,-2.0,0.0,0.0],
            [-3.0,0.0,0.0,1.0,0.0,0.0],
            [-5.0,0.0,0.0,2.0,0.0,0.0],
            [3.0,0.0,0.0,-1.0,0.0,0.0],
            [3.0,0.0,0.0,0.0,0.0,0.0],
            [3.0,0.0,0.0,0.0,0.0,0.0],
            [0.0,0.0,5.0,0.0,0.0,0.0],
            [0.0,0.0,0.0,1.0,0.0,0.0],
            [4.0,0.0,0.0,-2.0,0.0,0.0],
            [6.0,0.0,0.0,0.0,0.0,0.0],
            [5.0,0.0,0.0,-2.0,0.0,0.0],
            [-7.0,0.0,0.0,0.0,0.0,0.0],
            [-12.0,0.0,0.0,0.0,0.0,0.0],
            [5.0,0.0,0.0,-3.0,0.0,0.0],
            [3.0,0.0,0.0,-1.0,0.0,0.0],
            [-5.0,0.0,0.0,0.0,0.0,0.0],
            [3.0,0.0,0.0,0.0,0.0,0.0],
            [-7.0,0.0,0.0,3.0,0.0,0.0],
            [7.0,0.0,0.0,-4.0,0.0,0.0],
            [0.0,0.0,-12.0,0.0,0.0,-10.0],
            [4.0,0.0,0.0,-2.0,0.0,0.0],
            [3.0,0.0,0.0,-2.0,0.0,0.0],
            [-3.0,0.0,0.0,2.0,0.0,0.0],
            [-7.0,0.0,0.0,3.0,0.0,0.0],
            [-4.0,0.0,0.0,2.0,0.0,0.0],
            [-3.0,0.0,0.0,1.0,0.0,0.0],
            [0.0,0.0,0.0,0.0,0.0,0.0],
            [-3.0,0.0,0.0,1.0,0.0,0.0],
            [7.0,0.0,0.0,-3.0,0.0,0.0],
            [-4.0,0.0,0.0,2.0,0.0,0.0],
            [4.0,0.0,0.0,-2.0,0.0,0.0],
            [-5.0,0.0,0.0,3.0,0.0,0.0],
            [5.0,0.0,0.0,0.0,0.0,0.0],
            [-5.0,0.0,0.0,2.0,0.0,0.0],
            [5.0,0.0,0.0,-2.0,0.0,0.0],
            [-8.0,0.0,0.0,3.0,0.0,0.0],
            [9.0,0.0,0.0,0.0,0.0,0.0],
            [6.0,0.0,0.0,-3.0,0.0,0.0],
            [-5.0,0.0,0.0,2.0,0.0,0.0],
            [3.0,0.0,0.0,0.0,0.0,0.0],
            [-7.0,0.0,0.0,0.0,0.0,0.0],
            [-3.0,0.0,0.0,1.0,0.0,0.0],
            [5.0,0.0,0.0,0.0,0.0,0.0],
            [3.0,0.0,0.0,0.0,0.0,0.0],
            [-3.0,0.0,0.0,2.0,0.0,0.0],
            [4.0,0.0,0.0,-2.0,0.0,0.0],
            [3.0,0.0,0.0,-1.0,0.0,0.0],
            [-5.0,0.0,0.0,2.0,0.0,0.0],
            [4.0,0.0,0.0,-2.0,0.0,0.0],
            [9.0,0.0,0.0,-3.0,0.0,0.0],
            [4.0,0.0,0.0,0.0,0.0,0.0],
            [4.0,0.0,0.0,-2.0,0.0,0.0],
            [-3.0,0.0,0.0,2.0,0.0,0.0],
            [-4.0,0.0,0.0,2.0,0.0,0.0],
            [9.0,0.0,0.0,-3.0,0.0,0.0],
            [-4.0,0.0,0.0,0.0,0.0,0.0],
            [-4.0,0.0,0.0,0.0,0.0,0.0],
            [3.0,0.0,0.0,-2.0,0.0,0.0],
            [8.0,0.0,0.0,0.0,0.0,0.0],
            [3.0,0.0,0.0,0.0,0.0,0.0],
            [-3.0,0.0,0.0,2.0,0.0,0.0],
            [3.0,0.0,0.0,-1.0,0.0,0.0],
            [3.0,0.0,0.0,-1.0,0.0,0.0],
            [-3.0,0.0,0.0,1.0,0.0,0.0],
            [6.0,0.0,0.0,-3.0,0.0,0.0],
            [3.0,0.0,0.0,0.0,0.0,0.0],
            [-3.0,0.0,0.0,1.0,0.0,0.0],
            [-7.0,0.0,0.0,0.0,0.0,0.0],
            [9.0,0.0,0.0,0.0,0.0,0.0],
            [-3.0,0.0,0.0,2.0,0.0,0.0],
            [-3.0,0.0,0.0,0.0,0.0,0.0],
            [-4.0,0.0,0.0,0.0,0.0,0.0],
            [-5.0,0.0,0.0,3.0,0.0,0.0],
            [-13.0,0.0,0.0,0.0,0.0,0.0],
            [-7.0,0.0,0.0,0.0,0.0,0.0],
            [10.0,0.0,0.0,0.0,0.0,0.0],
            [3.0,0.0,0.0,-1.0,0.0,0.0],
            [10.0,0.0,13.0,6.0,0.0,-5.0],
            [0.0,0.0,30.0,0.0,0.0,14.0],
            [0.0,0.0,-162.0,0.0,0.0,-138.0],
            [0.0,0.0,75.0,0.0,0.0,0.0],
            [-7.0,0.0,0.0,4.0,0.0,0.0],
            [-4.0,0.0,0.0,2.0,0.0,0.0],
            [4.0,0.0,0.0,-2.0,0.0,0.0],
            [5.0,0.0,0.0,-2.0,0.0,0.0],
            [5.0,0.0,0.0,-3.0,0.0,0.0],
            [-3.0,0.0,0.0,0.0,0.0,0.0],
            [-3.0,0.0,0.0,2.0,0.0,0.0],
            [-4.0,0.0,0.0,2.0,0.0,0.0],
            [-5.0,0.0,0.0,2.0,0.0,0.0],
            [6.0,0.0,0.0,0.0,0.0,0.0],
            [9.0,0.0,0.0,0.0,0.0,0.0],
            [5.0,0.0,0.0,0.0,0.0,0.0],
            [-7.0,0.0,0.0,0.0,0.0,0.0],
            [-3.0,0.0,0.0,1.0,0.0,0.0],
            [-4.0,0.0,0.0,2.0,0.0,0.0],
            [7.0,0.0,0.0,0.0,0.0,0.0],
            [-4.0,0.0,0.0,0.0,0.0,0.0],
            [4.0,0.0,0.0,0.0,0.0,0.0],
            [-6.0,0.0,-3.0,3.0,0.0,1.0],
            [0.0,0.0,-3.0,0.0,0.0,-2.0],
            [11.0,0.0,0.0,0.0,0.0,0.0],
            [3.0,0.0,0.0,-1.0,0.0,0.0],
            [11.0,0.0,0.0,0.0,0.0,0.0],
            [-3.0,0.0,0.0,2.0,0.0,0.0],
            [-1.0,0.0,3.0,3.0,0.0,-1.0],
            [4.0,0.0,0.0,-2.0,0.0,0.0],
            [0.0,0.0,-13.0,0.0,0.0,-11.0],
            [3.0,0.0,6.0,0.0,0.0,0.0],
            [-7.0,0.0,0.0,0.0,0.0,0.0],
            [5.0,0.0,0.0,-3.0,0.0,0.0],
            [-3.0,0.0,0.0,1.0,0.0,0.0],
            [3.0,0.0,0.0,0.0,0.0,0.0],
            [5.0,0.0,0.0,-3.0,0.0,0.0],
            [-7.0,0.0,0.0,3.0,0.0,0.0],
            [8.0,0.0,0.0,-3.0,0.0,0.0],
            [-4.0,0.0,0.0,2.0,0.0,0.0],
            [11.0,0.0,0.0,0.0,0.0,0.0],
            [-3.0,0.0,0.0,1.0,0.0,0.0],
            [3.0,0.0,0.0,-1.0,0.0,0.0],
            [-4.0,0.0,0.0,2.0,0.0,0.0],
            [8.0,0.0,0.0,-4.0,0.0,0.0],
            [3.0,0.0,0.0,-1.0,0.0,0.0],
            [11.0,0.0,0.0,0.0,0.0,0.0],
            [-6.0,0.0,0.0,3.0,0.0,0.0],
            [-4.0,0.0,0.0,2.0,0.0,0.0],
            [-8.0,0.0,0.0,4.0,0.0,0.0],
            [-7.0,0.0,0.0,3.0,0.0,0.0],
            [-4.0,0.0,0.0,2.0,0.0,0.0],
            [3.0,0.0,0.0,-1.0,0.0,0.0],
            [6.0,0.0,0.0,-3.0,0.0,0.0],
            [-6.0,0.0,0.0,3.0,0.0,0.0],
            [6.0,0.0,0.0,0.0,0.0,0.0],
            [6.0,0.0,0.0,-1.0,0.0,0.0],
            [5.0,0.0,0.0,-2.0,0.0,0.0],
            [-5.0,0.0,0.0,2.0,0.0,0.0],
            [-4.0,0.0,0.0,0.0,0.0,0.0],
            [-4.0,0.0,0.0,2.0,0.0,0.0],
            [4.0,0.0,0.0,0.0,0.0,0.0],
            [6.0,0.0,0.0,-3.0,0.0,0.0],
            [-4.0,0.0,0.0,2.0,0.0,0.0],
            [0.0,0.0,-26.0,0.0,0.0,-11.0],
            [0.0,0.0,-10.0,0.0,0.0,-5.0],
            [5.0,0.0,0.0,-3.0,0.0,0.0],
            [-13.0,0.0,0.0,0.0,0.0,0.0],
            [3.0,0.0,0.0,-2.0,0.0,0.0],
            [4.0,0.0,0.0,-2.0,0.0,0.0],
            [7.0,0.0,0.0,-3.0,0.0,0.0],
            [4.0,0.0,0.0,0.0,0.0,0.0],
            [5.0,0.0,0.0,0.0,0.0,0.0],
            [-3.0,0.0,0.0,2.0,0.0,0.0],
            [-6.0,0.0,0.0,2.0,0.0,0.0],
            [-5.0,0.0,0.0,2.0,0.0,0.0],
            [-7.0,0.0,0.0,3.0,0.0,0.0],
            [5.0,0.0,0.0,-2.0,0.0,0.0],
            [13.0,0.0,0.0,0.0,0.0,0.0],
            [-4.0,0.0,0.0,2.0,0.0,0.0],
            [-3.0,0.0,0.0,0.0,0.0,0.0],
            [5.0,0.0,0.0,-2.0,0.0,0.0],
            [-11.0,0.0,0.0,0.0,0.0,0.0],
            [5.0,0.0,0.0,-2.0,0.0,0.0],
            [4.0,0.0,0.0,0.0,0.0,0.0],
            [4.0,0.0,0.0,-2.0,0.0,0.0],
            [-4.0,0.0,0.0,2.0,0.0,0.0],
            [6.0,0.0,0.0,-3.0,0.0,0.0],
            [3.0,0.0,0.0,-2.0,0.0,0.0],
            [-12.0,0.0,0.0,0.0,0.0,0.0],
            [4.0,0.0,0.0,0.0,0.0,0.0],
            [-3.0,0.0,0.0,0.0,0.0,0.0],
            [-4.0,0.0,0.0,0.0,0.0,0.0],
            [3.0,0.0,0.0,0.0,0.0,0.0],
            [3.0,0.0,0.0,-1.0,0.0,0.0],
            [-3.0,0.0,0.0,1.0,0.0,0.0],
            [0.0,0.0,-5.0,0.0,0.0,-2.0],
            [-7.0,0.0,0.0,4.0,0.0,0.0],
            [6.0,0.0,0.0,-3.0,0.0,0.0],
            [-3.0,0.0,0.0,0.0,0.0,0.0],
            [5.0,0.0,0.0,-3.0,0.0,0.0],
            [3.0,0.0,0.0,-1.0,0.0,0.0],
            [3.0,0.0,0.0,0.0,0.0,0.0],
            [-3.0,0.0,0.0,1.0,0.0,0.0],
            [-5.0,0.0,0.0,3.0,0.0,0.0],
            [-3.0,0.0,0.0,2.0,0.0,0.0],
            [-3.0,0.0,0.0,2.0,0.0,0.0],
            [12.0,0.0,0.0,0.0,0.0,0.0],
            [3.0,0.0,0.0,-1.0,0.0,0.0],
            [-4.0,0.0,0.0,2.0,0.0,0.0],
            [4.0,0.0,0.0,0.0,0.0,0.0],
            [6.0,0.0,0.0,0.0,0.0,0.0],
            [5.0,0.0,0.0,-3.0,0.0,0.0],
            [4.0,0.0,0.0,-2.0,0.0,0.0],
            [-6.0,0.0,0.0,3.0,0.0,0.0],
            [4.0,0.0,0.0,-2.0,0.0,0.0],
            [6.0,0.0,0.0,-3.0,0.0,0.0],
            [6.0,0.0,0.0,0.0,0.0,0.0],
            [-6.0,0.0,0.0,3.0,0.0,0.0],
            [3.0,0.0,0.0,-2.0,0.0,0.0],
            [7.0,0.0,0.0,-4.0,0.0,0.0],
            [4.0,0.0,0.0,-2.0,0.0,0.0],
            [-5.0,0.0,0.0,2.0,0.0,0.0],
            [5.0,0.0,0.0,0.0,0.0,0.0],
            [-6.0,0.0,0.0,3.0,0.0,0.0],
            [-6.0,0.0,0.0,3.0,0.0,0.0],
            [-4.0,0.0,0.0,2.0,0.0,0.0],
            [10.0,0.0,0.0,0.0,0.0,0.0],
            [-4.0,0.0,0.0,2.0,0.0,0.0],
            [7.0,0.0,0.0,0.0,0.0,0.0],
            [7.0,0.0,0.0,-3.0,0.0,0.0],
            [4.0,0.0,0.0,0.0,0.0,0.0],
            [11.0,0.0,0.0,0.0,0.0,0.0],
            [5.0,0.0,0.0,-2.0,0.0,0.0],
            [-6.0,0.0,0.0,2.0,0.0,0.0],
            [4.0,0.0,0.0,-2.0,0.0,0.0],
            [3.0,0.0,0.0,-2.0,0.0,0.0],
            [5.0,0.0,0.0,-2.0,0.0,0.0],
            [-4.0,0.0,0.0,2.0,0.0,0.0],
            [-4.0,0.0,0.0,2.0,0.0,0.0],
            [-3.0,0.0,0.0,2.0,0.0,0.0],
            [4.0,0.0,0.0,-2.0,0.0,0.0],
            [3.0,0.0,0.0,-1.0,0.0,0.0],
            [-3.0,0.0,0.0,1.0,0.0,0.0],
            [-3.0,0.0,0.0,1.0,0.0,0.0],
            [-3.0,0.0,0.0,2.0,0.0,0.0]]
    
    #行星章动参数
    #*   L   L'  F   D   Om  Me  Ve  E  Ma  Ju  Sa  Ur  Ne  pre
    NAPL=[[0,0,0,0,0,0,0,8,-16,4,5,0,0,0],
            [0,0,0,0,0,0,0,-8,16,-4,-5,0,0,2],
            [0,0,0,0,0,0,0,8,-16,4,5,0,0,2],
            [0,0,0,0,0,0,0,0,0,0,0,-1,2,2],
            [0,0,0,0,0,0,0,-4,8,-1,-5,0,0,2],
            [0,0,0,0,0,0,0,4,-8,3,0,0,0,1],
            [0,0,1,-1,1,0,0,3,-8,3,0,0,0,0],
            [-1,0,0,0,0,0,10,-3,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,-2,6,-3,0,2],
            [0,0,0,0,0,0,0,4,-8,3,0,0,0,0],
            [0,0,1,-1,1,0,0,-5,8,-3,0,0,0,0],
            [0,0,0,0,0,0,0,-4,8,-3,0,0,0,1],
            [0,0,0,0,0,0,0,4,-8,1,5,0,0,2],
            [0,0,0,0,0,0,-5,6,4,0,0,0,0,2],
            [0,0,0,0,0,0,0,0,0,2,-5,0,0,2],
            [0,0,0,0,0,0,0,0,0,2,-5,0,0,1],
            [0,0,1,-1,1,0,0,-1,0,2,-5,0,0,0],
            [0,0,0,0,0,0,0,0,0,2,-5,0,0,0],
            [0,0,1,-1,1,0,0,-1,0,-2,5,0,0,0],
            [0,0,0,0,0,0,0,0,0,-2,5,0,0,1],
            [0,0,0,0,0,0,0,0,0,-2,5,0,0,2],
            [2,0,-1,-1,0,0,0,3,-7,0,0,0,0,0],
            [1,0,0,-2,0,0,19,-21,3,0,0,0,0,0],
            [0,0,1,-1,1,0,2,-4,0,-3,0,0,0,0],
            [1,0,0,-1,1,0,0,-1,0,2,0,0,0,0],
            [0,0,1,-1,1,0,0,-1,0,-4,10,0,0,0],
            [-2,0,0,2,1,0,0,2,0,0,-5,0,0,0],
            [0,0,0,0,0,0,3,-7,4,0,0,0,0,0],
            [0,0,-1,1,0,0,0,1,0,1,-1,0,0,0],
            [-2,0,0,2,1,0,0,2,0,-2,0,0,0,0],
            [-1,0,0,0,0,0,18,-16,0,0,0,0,0,0],
            [-2,0,1,1,2,0,0,1,0,-2,0,0,0,0],
            [-1,0,1,-1,1,0,18,-17,0,0,0,0,0,0],
            [-1,0,0,1,1,0,0,2,-2,0,0,0,0,0],
            [0,0,0,0,0,0,-8,13,0,0,0,0,0,2],
            [0,0,2,-2,2,0,-8,11,0,0,0,0,0,0],
            [0,0,0,0,0,0,-8,13,0,0,0,0,0,1],
            [0,0,1,-1,1,0,-8,12,0,0,0,0,0,0],
            [0,0,0,0,0,0,8,-13,0,0,0,0,0,0],
            [0,0,1,-1,1,0,8,-14,0,0,0,0,0,0],
            [0,0,0,0,0,0,8,-13,0,0,0,0,0,1],
            [-2,0,0,2,1,0,0,2,0,-4,5,0,0,0],
            [-2,0,0,2,2,0,3,-3,0,0,0,0,0,0],
            [-2,0,0,2,0,0,0,2,0,-3,1,0,0,0],
            [0,0,0,0,1,0,3,-5,0,2,0,0,0,0],
            [-2,0,0,2,0,0,0,2,0,-4,3,0,0,0],
            [0,0,-1,1,0,0,0,0,2,0,0,0,0,0],
            [0,0,0,0,1,0,0,-1,2,0,0,0,0,0],
            [0,0,1,-1,2,0,0,-2,2,0,0,0,0,0],
            [-1,0,1,0,1,0,3,-5,0,0,0,0,0,0],
            [-1,0,0,1,0,0,3,-4,0,0,0,0,0,0],
            [-2,0,0,2,0,0,0,2,0,-2,-2,0,0,0],
            [-2,0,2,0,2,0,0,-5,9,0,0,0,0,0],
            [0,0,1,-1,1,0,0,-1,0,0,0,-1,0,0],
            [0,0,0,0,0,0,0,0,0,0,0,1,0,0],
            [0,0,1,-1,1,0,0,-1,0,0,0,0,2,0],
            [0,0,0,0,0,0,0,0,0,0,0,0,2,1],
            [0,0,0,0,0,0,0,0,0,0,0,0,2,2],
            [-1,0,0,1,0,0,0,3,-4,0,0,0,0,0],
            [0,0,-1,1,0,0,0,1,0,0,2,0,0,0],
            [0,0,1,-1,2,0,0,-1,0,0,2,0,0,0],
            [0,0,0,0,1,0,0,-9,17,0,0,0,0,0],
            [0,0,0,0,2,0,-3,5,0,0,0,0,0,0],
            [0,0,1,-1,1,0,0,-1,0,-1,2,0,0,0],
            [0,0,0,0,0,0,0,0,0,1,-2,0,0,0],
            [1,0,0,-2,0,0,17,-16,0,-2,0,0,0,0],
            [0,0,1,-1,1,0,0,-1,0,1,-3,0,0,0],
            [-2,0,0,2,1,0,0,5,-6,0,0,0,0,0],
            [0,0,-2,2,0,0,0,9,-13,0,0,0,0,0],
            [0,0,1,-1,2,0,0,-1,0,0,1,0,0,0],
            [0,0,0,0,1,0,0,0,0,0,1,0,0,0],
            [0,0,-1,1,0,0,0,1,0,0,1,0,0,0],
            [0,0,-2,2,0,0,5,-6,0,0,0,0,0,0],
            [0,0,-1,1,1,0,5,-7,0,0,0,0,0,0],
            [-2,0,0,2,0,0,6,-8,0,0,0,0,0,0],
            [2,0,1,-3,1,0,-6,7,0,0,0,0,0,0],
            [0,0,0,0,2,0,0,0,0,1,0,0,0,0],
            [0,0,-1,1,1,0,0,1,0,1,0,0,0,0],
            [0,0,1,-1,1,0,0,-1,0,0,0,2,0,0],
            [0,0,0,0,0,0,0,0,0,0,0,2,0,1],
            [0,0,0,0,0,0,0,0,0,0,0,2,0,2],
            [0,0,0,0,0,0,0,-8,15,0,0,0,0,2],
            [0,0,0,0,0,0,0,-8,15,0,0,0,0,1],
            [0,0,1,-1,1,0,0,-9,15,0,0,0,0,0],
            [0,0,0,0,0,0,0,8,-15,0,0,0,0,0],
            [1,0,-1,-1,0,0,0,8,-15,0,0,0,0,0],
            [2,0,0,-2,0,0,2,-5,0,0,0,0,0,0],
            [-2,0,0,2,0,0,0,2,0,-5,5,0,0,0],
            [2,0,0,-2,1,0,0,-6,8,0,0,0,0,0],
            [2,0,0,-2,1,0,0,-2,0,3,0,0,0,0],
            [-2,0,1,1,0,0,0,1,0,-3,0,0,0,0],
            [-2,0,1,1,1,0,0,1,0,-3,0,0,0,0],
            [-2,0,0,2,0,0,0,2,0,-3,0,0,0,0],
            [-2,0,0,2,0,0,0,6,-8,0,0,0,0,0],
            [-2,0,0,2,0,0,0,2,0,-1,-5,0,0,0],
            [-1,0,0,1,0,0,0,1,0,-1,0,0,0,0],
            [-1,0,1,1,1,0,-20,20,0,0,0,0,0,0],
            [1,0,0,-2,0,0,20,-21,0,0,0,0,0,0],
            [0,0,0,0,1,0,0,8,-15,0,0,0,0,0],
            [0,0,2,-2,1,0,0,-10,15,0,0,0,0,0],
            [0,0,-1,1,0,0,0,1,0,1,0,0,0,0],
            [0,0,0,0,1,0,0,0,0,1,0,0,0,0],
            [0,0,1,-1,2,0,0,-1,0,1,0,0,0,0],
            [0,0,1,-1,1,0,0,-1,0,-2,4,0,0,0],
            [2,0,0,-2,1,0,-6,8,0,0,0,0,0,0],
            [0,0,-2,2,1,0,5,-6,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,0,-1,0,0,1],
            [0,0,1,-1,1,0,0,-1,0,0,-1,0,0,0],
            [0,0,0,0,0,0,0,0,0,0,1,0,0,0],
            [0,0,1,-1,1,0,0,-1,0,0,1,0,0,0],
            [0,0,0,0,0,0,0,0,0,0,1,0,0,1],
            [0,0,0,0,0,0,0,0,0,0,1,0,0,2],
            [0,0,2,-2,1,0,0,-9,13,0,0,0,0,0],
            [0,0,0,0,1,0,0,7,-13,0,0,0,0,0],
            [-2,0,0,2,0,0,0,5,-6,0,0,0,0,0],
            [0,0,0,0,0,0,0,9,-17,0,0,0,0,0],
            [0,0,0,0,0,0,0,-9,17,0,0,0,0,2],
            [1,0,0,-1,1,0,0,-3,4,0,0,0,0,0],
            [1,0,0,-1,1,0,-3,4,0,0,0,0,0,0],
            [0,0,0,0,2,0,0,-1,2,0,0,0,0,0],
            [0,0,-1,1,1,0,0,0,2,0,0,0,0,0],
            [0,0,-2,2,0,1,0,-2,0,0,0,0,0,0],
            [0,0,0,0,0,0,3,-5,0,2,0,0,0,0],
            [-2,0,0,2,1,0,0,2,0,-3,1,0,0,0],
            [-2,0,0,2,1,0,3,-3,0,0,0,0,0,0],
            [0,0,0,0,1,0,8,-13,0,0,0,0,0,0],
            [0,0,-1,1,0,0,8,-12,0,0,0,0,0,0],
            [0,0,2,-2,1,0,-8,11,0,0,0,0,0,0],
            [-1,0,0,1,0,0,0,2,-2,0,0,0,0,0],
            [-1,0,0,0,1,0,18,-16,0,0,0,0,0,0],
            [0,0,1,-1,1,0,0,-1,0,-1,1,0,0,0],
            [0,0,0,0,1,0,3,-7,4,0,0,0,0,0],
            [-2,0,1,1,1,0,0,-3,7,0,0,0,0,0],
            [0,0,1,-1,2,0,0,-1,0,-2,5,0,0,0],
            [0,0,0,0,1,0,0,0,0,-2,5,0,0,0],
            [0,0,0,0,1,0,0,-4,8,-3,0,0,0,0],
            [1,0,0,0,1,0,-10,3,0,0,0,0,0,0],
            [0,0,2,-2,1,0,0,-2,0,0,0,0,0,0],
            [-1,0,0,0,1,0,10,-3,0,0,0,0,0,0],
            [0,0,0,0,1,0,0,4,-8,3,0,0,0,0],
            [0,0,0,0,1,0,0,0,0,2,-5,0,0,0],
            [0,0,-1,1,0,0,0,1,0,2,-5,0,0,0],
            [2,0,-1,-1,1,0,0,3,-7,0,0,0,0,0],
            [-2,0,0,2,0,0,0,2,0,0,-5,0,0,0],
            [0,0,0,0,1,0,-3,7,-4,0,0,0,0,0],
            [-2,0,0,2,0,0,0,2,0,-2,0,0,0,0],
            [1,0,0,0,1,0,-18,16,0,0,0,0,0,0],
            [-2,0,1,1,1,0,0,1,0,-2,0,0,0,0],
            [0,0,1,-1,2,0,-8,12,0,0,0,0,0,0],
            [0,0,0,0,1,0,-8,13,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,1,-2,0,0,0,0,1],
            [0,0,1,-1,1,0,0,0,-2,0,0,0,0,0],
            [0,0,0,0,0,0,0,1,-2,0,0,0,0,0],
            [0,0,1,-1,1,0,0,-2,2,0,0,0,0,0],
            [0,0,0,0,0,0,0,-1,2,0,0,0,0,1],
            [-1,0,0,1,1,0,3,-4,0,0,0,0,0,0],
            [-1,0,0,1,1,0,0,3,-4,0,0,0,0,0],
            [0,0,1,-1,1,0,0,-1,0,0,-2,0,0,0],
            [0,0,1,-1,1,0,0,-1,0,0,2,0,0,0],
            [0,0,0,0,0,0,0,0,0,0,2,0,0,1],
            [0,0,0,0,0,0,0,0,0,0,2,0,0,2],
            [0,0,1,-1,0,0,3,-6,0,0,0,0,0,0],
            [0,0,0,0,1,0,-3,5,0,0,0,0,0,0],
            [0,0,1,-1,2,0,-3,4,0,0,0,0,0,0],
            [0,0,0,0,1,0,0,-2,4,0,0,0,0,0],
            [0,0,2,-2,1,0,-5,6,0,0,0,0,0,0],
            [0,0,-1,1,0,0,5,-7,0,0,0,0,0,0],
            [0,0,0,0,1,0,5,-8,0,0,0,0,0,0],
            [-2,0,0,2,1,0,6,-8,0,0,0,0,0,0],
            [0,0,0,0,1,0,0,-8,15,0,0,0,0,0],
            [-2,0,0,2,1,0,0,2,0,-3,0,0,0,0],
            [-2,0,0,2,1,0,0,6,-8,0,0,0,0,0],
            [1,0,0,-1,1,0,0,-1,0,1,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,3,-5,0,0,0],
            [0,0,1,-1,1,0,0,-1,0,-1,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,-1,0,0,0,1],
            [0,0,0,0,0,0,0,0,0,1,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,1,0,0,0,1],
            [0,0,1,-1,1,0,0,-1,0,1,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,1,0,0,0,1],
            [0,0,0,0,0,0,0,0,0,1,0,0,0,2],
            [0,0,1,-1,2,0,0,-1,0,0,-1,0,0,0],
            [0,0,0,0,1,0,0,0,0,0,-1,0,0,0],
            [0,0,-1,1,0,0,0,1,0,0,-1,0,0,0],
            [0,0,0,0,0,0,0,-7,13,0,0,0,0,2],
            [0,0,0,0,0,0,0,7,-13,0,0,0,0,0],
            [2,0,0,-2,1,0,0,-5,6,0,0,0,0,0],
            [0,0,2,-2,1,0,0,-8,11,0,0,0,0,0],
            [0,0,2,-2,1,-1,0,2,0,0,0,0,0,0],
            [-2,0,0,2,0,0,0,4,-4,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,2,-2,0,0,0],
            [0,0,1,-1,1,0,0,-1,0,0,3,0,0,0],
            [0,0,0,0,0,0,0,0,0,0,3,0,0,1],
            [0,0,0,0,0,0,0,0,0,0,3,0,0,2],
            [-2,0,0,2,0,0,3,-3,0,0,0,0,0,0],
            [0,0,0,0,2,0,0,-4,8,-3,0,0,0,0],
            [0,0,0,0,2,0,0,4,-8,3,0,0,0,0],
            [2,0,0,-2,1,0,0,-2,0,2,0,0,0,0],
            [0,0,1,-1,2,0,0,-1,0,2,0,0,0,0],
            [0,0,1,-1,2,0,0,0,-2,0,0,0,0,0],
            [0,0,0,0,1,0,0,1,-2,0,0,0,0,0],
            [0,0,-1,1,0,0,0,2,-2,0,0,0,0,0],
            [0,0,-1,1,0,0,0,1,0,0,-2,0,0,0],
            [0,0,2,-2,1,0,0,-2,0,0,2,0,0,0],
            [0,0,1,-1,1,0,3,-6,0,0,0,0,0,0],
            [0,0,0,0,0,0,3,-5,0,0,0,0,0,1],
            [0,0,0,0,0,0,3,-5,0,0,0,0,0,0],
            [0,0,1,-1,1,0,-3,4,0,0,0,0,0,0],
            [0,0,0,0,0,0,-3,5,0,0,0,0,0,1],
            [0,0,0,0,0,0,-3,5,0,0,0,0,0,2],
            [0,0,2,-2,2,0,-3,3,0,0,0,0,0,0],
            [0,0,0,0,0,0,-3,5,0,0,0,0,0,2],
            [0,0,0,0,0,0,0,2,-4,0,0,0,0,1],
            [0,0,1,-1,1,0,0,1,-4,0,0,0,0,0],
            [0,0,0,0,0,0,0,2,-4,0,0,0,0,0],
            [0,0,0,0,0,0,0,-2,4,0,0,0,0,1],
            [0,0,1,-1,1,0,0,-3,4,0,0,0,0,0],
            [0,0,0,0,0,0,0,-2,4,0,0,0,0,1],
            [0,0,0,0,0,0,0,-2,4,0,0,0,0,2],
            [0,0,0,0,0,0,-5,8,0,0,0,0,0,2],
            [0,0,2,-2,2,0,-5,6,0,0,0,0,0,0],
            [0,0,0,0,0,0,-5,8,0,0,0,0,0,2],
            [0,0,0,0,0,0,-5,8,0,0,0,0,0,1],
            [0,0,1,-1,1,0,-5,7,0,0,0,0,0,0],
            [0,0,0,0,0,0,-5,8,0,0,0,0,0,1],
            [0,0,0,0,0,0,5,-8,0,0,0,0,0,0],
            [0,0,1,-1,2,0,0,-1,0,-1,0,0,0,0],
            [0,0,0,0,1,0,0,0,0,-1,0,0,0,0],
            [0,0,-1,1,0,0,0,1,0,-1,0,0,0,0],
            [0,0,2,-2,1,0,0,-2,0,1,0,0,0,0],
            [0,0,0,0,0,0,0,-6,11,0,0,0,0,2],
            [0,0,0,0,0,0,0,6,-11,0,0,0,0,0],
            [0,0,0,0,0,-1,0,4,0,0,0,0,0,2],
            [0,0,0,0,0,1,0,-4,0,0,0,0,0,0],
            [2,0,0,-2,1,0,-3,3,0,0,0,0,0,0],
            [-2,0,0,2,0,0,0,2,0,0,-2,0,0,0],
            [0,0,2,-2,1,0,0,-7,9,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,4,-5,0,0,2],
            [0,0,0,0,0,0,0,0,0,2,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,2,0,0,0,1],
            [0,0,1,-1,1,0,0,-1,0,2,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,2,0,0,0,1],
            [0,0,0,0,0,0,0,0,0,2,0,0,0,2],
            [0,0,2,-2,2,0,0,-2,0,2,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,0,5,0,0,2],
            [0,0,0,0,1,0,3,-5,0,0,0,0,0,0],
            [0,0,-1,1,0,0,3,-4,0,0,0,0,0,0],
            [0,0,2,-2,1,0,-3,3,0,0,0,0,0,0],
            [0,0,0,0,1,0,0,2,-4,0,0,0,0,0],
            [0,0,2,-2,1,0,0,-4,4,0,0,0,0,0],
            [0,0,1,-1,2,0,-5,7,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,3,-6,0,0,0,0,0],
            [0,0,0,0,0,0,0,-3,6,0,0,0,0,1],
            [0,0,1,-1,1,0,0,-4,6,0,0,0,0,0],
            [0,0,0,0,0,0,0,-3,6,0,0,0,0,1],
            [0,0,0,0,0,0,0,-3,6,0,0,0,0,2],
            [0,0,-1,1,0,0,2,-2,0,0,0,0,0,0],
            [0,0,0,0,1,0,2,-3,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,-5,9,0,0,0,0,2],
            [0,0,0,0,0,0,0,-5,9,0,0,0,0,1],
            [0,0,0,0,0,0,0,5,-9,0,0,0,0,0],
            [0,0,-1,1,0,0,0,1,0,-2,0,0,0,0],
            [0,0,2,-2,1,0,0,-2,0,2,0,0,0,0],
            [-2,0,1,1,1,0,0,1,0,0,0,0,0,0],
            [0,0,-2,2,0,0,3,-3,0,0,0,0,0,0],
            [0,0,0,0,0,0,-6,10,0,0,0,0,0,1],
            [0,0,0,0,0,0,-6,10,0,0,0,0,0,2],
            [0,0,0,0,0,0,-2,3,0,0,0,0,0,2],
            [0,0,0,0,0,0,-2,3,0,0,0,0,0,1],
            [0,0,1,-1,1,0,-2,2,0,0,0,0,0,0],
            [0,0,0,0,0,0,2,-3,0,0,0,0,0,0],
            [0,0,0,0,0,0,2,-3,0,0,0,0,0,1],
            [0,0,0,0,0,0,0,0,0,3,0,0,0,1],
            [0,0,1,-1,1,0,0,-1,0,3,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,3,0,0,0,1],
            [0,0,0,0,0,0,0,0,0,3,0,0,0,2],
            [0,0,0,0,0,0,0,4,-8,0,0,0,0,0],
            [0,0,0,0,0,0,0,-4,8,0,0,0,0,2],
            [0,0,-2,2,0,0,0,2,0,-2,0,0,0,0],
            [0,0,0,0,0,0,0,-4,7,0,0,0,0,2],
            [0,0,0,0,0,0,0,-4,7,0,0,0,0,1],
            [0,0,0,0,0,0,0,4,-7,0,0,0,0,0],
            [0,0,0,0,1,0,-2,3,0,0,0,0,0,0],
            [0,0,2,-2,1,0,0,-2,0,3,0,0,0,0],
            [0,0,0,0,0,0,0,-5,10,0,0,0,0,2],
            [0,0,0,0,1,0,-1,2,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,4,0,0,0,2],
            [0,0,0,0,0,0,0,-3,5,0,0,0,0,2],
            [0,0,0,0,0,0,0,-3,5,0,0,0,0,1],
            [0,0,0,0,0,0,0,3,-5,0,0,0,0,0],
            [0,0,0,0,0,0,1,-2,0,0,0,0,0,1],
            [0,0,1,-1,1,0,1,-3,0,0,0,0,0,0],
            [0,0,0,0,0,0,1,-2,0,0,0,0,0,0],
            [0,0,0,0,0,0,-1,2,0,0,0,0,0,1],
            [0,0,0,0,0,0,-1,2,0,0,0,0,0,2],
            [0,0,0,0,0,0,-7,11,0,0,0,0,0,2],
            [0,0,0,0,0,0,-7,11,0,0,0,0,0,1],
            [0,0,-2,2,0,0,4,-4,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,2,-3,0,0,0,0,0],
            [0,0,2,-2,1,0,-4,4,0,0,0,0,0,0],
            [0,0,-1,1,0,0,4,-5,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,1,-1,0,0,0,0,0],
            [0,0,0,0,0,0,-4,7,0,0,0,0,0,1],
            [0,0,1,-1,1,0,-4,6,0,0,0,0,0,0],
            [0,0,0,0,0,0,-4,7,0,0,0,0,0,2],
            [0,0,0,0,0,0,-4,6,0,0,0,0,0,2],
            [0,0,0,0,0,0,-4,6,0,0,0,0,0,1],
            [0,0,1,-1,1,0,-4,5,0,0,0,0,0,0],
            [0,0,0,0,0,0,-4,6,0,0,0,0,0,1],
            [0,0,0,0,0,0,4,-6,0,0,0,0,0,0],
            [-2,0,0,2,0,0,2,-2,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,1,0,0,0,0,0],
            [0,0,-1,1,0,0,1,0,0,0,0,0,0,0],
            [0,0,0,0,1,0,1,-1,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,-1,0,5,0,0,0,2],
            [0,0,0,0,0,0,0,1,-3,0,0,0,0,0],
            [0,0,0,0,0,0,0,-1,3,0,0,0,0,2],
            [0,0,0,0,0,0,0,-7,12,0,0,0,0,2],
            [0,0,0,0,0,0,-1,1,0,0,0,0,0,2],
            [0,0,0,0,0,0,-1,1,0,0,0,0,0,1],
            [0,0,1,-1,1,0,-1,0,0,0,0,0,0,0],
            [0,0,0,0,0,0,1,-1,0,0,0,0,0,0],
            [0,0,0,0,0,0,1,-1,0,0,0,0,0,1],
            [0,0,1,-1,1,0,1,-2,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,-2,5,0,0,0,0,2],
            [0,0,0,0,0,0,0,-1,0,4,0,0,0,2],
            [0,0,0,0,0,0,0,1,0,-4,0,0,0,0],
            [0,0,0,0,1,0,-1,1,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,-6,10,0,0,0,0,2],
            [0,0,0,0,0,0,0,-6,10,0,0,0,0,0],
            [0,0,2,-2,1,0,0,-3,0,3,0,0,0,0],
            [0,0,0,0,0,0,0,-3,7,0,0,0,0,2],
            [-2,0,0,2,0,0,4,-4,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,-5,8,0,0,0,0,2],
            [0,0,0,0,0,0,0,5,-8,0,0,0,0,0],
            [0,0,0,0,0,0,0,-1,0,3,0,0,0,2],
            [0,0,0,0,0,0,0,-1,0,3,0,0,0,1],
            [0,0,0,0,0,0,0,1,0,-3,0,0,0,0],
            [0,0,0,0,0,0,2,-4,0,0,0,0,0,0],
            [0,0,0,0,0,0,-2,4,0,0,0,0,0,1],
            [0,0,1,-1,1,0,-2,3,0,0,0,0,0,0],
            [0,0,0,0,0,0,-2,4,0,0,0,0,0,2],
            [0,0,0,0,0,0,-6,9,0,0,0,0,0,2],
            [0,0,0,0,0,0,-6,9,0,0,0,0,0,1],
            [0,0,0,0,0,0,6,-9,0,0,0,0,0,0],
            [0,0,0,0,1,0,0,1,0,-2,0,0,0,0],
            [0,0,2,-2,1,0,-2,2,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,-4,6,0,0,0,0,2],
            [0,0,0,0,0,0,0,4,-6,0,0,0,0,0],
            [0,0,0,0,1,0,3,-4,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,-1,0,2,0,0,0,2],
            [0,0,0,0,0,0,0,1,0,-2,0,0,0,0],
            [0,0,0,0,1,0,0,1,0,-1,0,0,0,0],
            [0,0,0,0,0,0,-5,9,0,0,0,0,0,2],
            [0,0,0,0,0,0,0,3,-4,0,0,0,0,0],
            [0,0,0,0,0,0,-3,4,0,0,0,0,0,2],
            [0,0,0,0,0,0,-3,4,0,0,0,0,0,1],
            [0,0,0,0,0,0,3,-4,0,0,0,0,0,0],
            [0,0,0,0,0,0,3,-4,0,0,0,0,0,1],
            [0,0,0,0,1,0,0,2,-2,0,0,0,0,0],
            [0,0,0,0,1,0,0,-1,0,2,0,0,0,0],
            [0,0,0,0,0,0,0,1,0,0,-3,0,0,0],
            [0,0,0,0,0,0,0,1,0,1,-5,0,0,0],
            [0,0,0,0,0,0,0,-1,0,1,0,0,0,1],
            [0,0,0,0,0,0,0,1,0,-1,0,0,0,0],
            [0,0,0,0,0,0,0,1,0,-1,0,0,0,1],
            [0,0,0,0,0,0,0,1,0,-3,5,0,0,0],
            [0,0,0,0,1,0,-3,4,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,1,0,0,-2,0,0,0],
            [0,0,0,0,0,0,0,2,-2,0,0,0,0,0],
            [0,0,0,0,0,0,0,1,0,0,-1,0,0,0],
            [0,0,0,0,1,0,0,-1,0,1,0,0,0,0],
            [0,0,0,0,1,0,0,-2,2,0,0,0,0,0],
            [0,0,0,0,0,0,-8,14,0,0,0,0,0,2],
            [0,0,0,0,0,0,0,1,0,2,-5,0,0,0],
            [0,0,0,0,0,0,0,5,-8,3,0,0,0,0],
            [0,0,0,0,0,0,0,5,-8,3,0,0,0,2],
            [0,0,0,0,0,0,0,-1,0,0,0,0,0,1],
            [0,0,0,0,0,0,0,1,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,3,-8,3,0,0,0,0],
            [0,0,0,0,0,0,0,-3,8,-3,0,0,0,2],
            [0,0,0,0,0,0,0,1,0,-2,5,0,0,2],
            [0,0,0,0,0,0,-8,12,0,0,0,0,0,2],
            [0,0,0,0,0,0,-8,12,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,1,0,1,-2,0,0,0],
            [0,0,0,0,0,0,0,1,0,0,1,0,0,2],
            [0,0,0,0,0,0,0,0,2,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,2,0,0,0,0,2],
            [0,0,0,0,0,0,0,1,0,0,2,0,0,2],
            [0,0,2,-2,1,0,-5,5,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,1,0,1,0,0,0,0],
            [0,0,0,0,0,0,0,1,0,1,0,0,0,1],
            [0,0,0,0,0,0,0,1,0,1,0,0,0,2],
            [0,0,0,0,0,0,3,-6,0,0,0,0,0,0],
            [0,0,0,0,0,0,-3,6,0,0,0,0,0,1],
            [0,0,0,0,0,0,-3,6,0,0,0,0,0,2],
            [0,0,0,0,0,0,0,-1,4,0,0,0,0,2],
            [0,0,0,0,0,0,-5,7,0,0,0,0,0,2],
            [0,0,0,0,0,0,-5,7,0,0,0,0,0,1],
            [0,0,1,-1,1,0,-5,6,0,0,0,0,0,0],
            [0,0,0,0,0,0,5,-7,0,0,0,0,0,0],
            [0,0,2,-2,1,0,0,-1,0,1,0,0,0,0],
            [0,0,0,0,0,0,0,-1,0,1,0,0,0,0],
            [0,0,0,0,0,-1,0,3,0,0,0,0,0,2],
            [0,0,0,0,0,0,0,1,0,2,0,0,0,2],
            [0,0,0,0,0,0,0,-2,6,0,0,0,0,2],
            [0,0,0,0,1,0,2,-2,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,-6,9,0,0,0,0,2],
            [0,0,0,0,0,0,0,6,-9,0,0,0,0,0],
            [0,0,0,0,0,0,-2,2,0,0,0,0,0,1],
            [0,0,1,-1,1,0,-2,1,0,0,0,0,0,0],
            [0,0,0,0,0,0,2,-2,0,0,0,0,0,0],
            [0,0,0,0,0,0,2,-2,0,0,0,0,0,1],
            [0,0,0,0,0,0,0,1,0,3,0,0,0,2],
            [0,0,0,0,0,0,0,-5,7,0,0,0,0,2],
            [0,0,0,0,0,0,0,5,-7,0,0,0,0,0],
            [0,0,0,0,1,0,-2,2,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,4,-5,0,0,0,0,0],
            [0,0,0,0,0,0,1,-3,0,0,0,0,0,0],
            [0,0,0,0,0,0,-1,3,0,0,0,0,0,1],
            [0,0,1,-1,1,0,-1,2,0,0,0,0,0,0],
            [0,0,0,0,0,0,-1,3,0,0,0,0,0,2],
            [0,0,0,0,0,0,-7,10,0,0,0,0,0,2],
            [0,0,0,0,0,0,-7,10,0,0,0,0,0,1],
            [0,0,0,0,0,0,0,3,-3,0,0,0,0,0],
            [0,0,0,0,0,0,-4,8,0,0,0,0,0,2],
            [0,0,0,0,0,0,-4,5,0,0,0,0,0,2],
            [0,0,0,0,0,0,-4,5,0,0,0,0,0,1],
            [0,0,0,0,0,0,4,-5,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,1,1,0,0,0,0,2],
            [0,0,0,0,0,0,0,-2,0,5,0,0,0,2],
            [0,0,0,0,0,0,0,0,3,0,0,0,0,2],
            [0,0,0,0,0,0,1,0,0,0,0,0,0,0],
            [0,0,0,0,0,0,1,0,0,0,0,0,0,2],
            [0,0,0,0,0,0,-9,13,0,0,0,0,0,2],
            [0,0,0,0,0,0,0,-1,5,0,0,0,0,2],
            [0,0,0,0,0,0,0,-2,0,4,0,0,0,2],
            [0,0,0,0,0,0,0,2,0,-4,0,0,0,0],
            [0,0,0,0,0,0,0,-2,7,0,0,0,0,2],
            [0,0,0,0,0,0,0,2,0,-3,0,0,0,0],
            [0,0,0,0,0,0,-2,5,0,0,0,0,0,1],
            [0,0,0,0,0,0,-2,5,0,0,0,0,0,2],
            [0,0,0,0,0,0,-6,8,0,0,0,0,0,2],
            [0,0,0,0,0,0,-6,8,0,0,0,0,0,1],
            [0,0,0,0,0,0,6,-8,0,0,0,0,0,0],
            [0,0,0,0,1,0,0,2,0,-2,0,0,0,0],
            [0,0,0,0,0,0,0,-3,9,0,0,0,0,2],
            [0,0,0,0,0,0,0,5,-6,0,0,0,0,0],
            [0,0,0,0,0,0,0,5,-6,0,0,0,0,2],
            [0,0,0,0,0,0,0,2,0,-2,0,0,0,0],
            [0,0,0,0,0,0,0,2,0,-2,0,0,0,1],
            [0,0,0,0,0,0,0,2,0,-2,0,0,0,2],
            [0,0,0,0,0,0,-5,10,0,0,0,0,0,2],
            [0,0,0,0,0,0,0,4,-4,0,0,0,0,0],
            [0,0,0,0,0,0,0,4,-4,0,0,0,0,2],
            [0,0,0,0,0,0,-3,3,0,0,0,0,0,1],
            [0,0,0,0,0,0,3,-3,0,0,0,0,0,0],
            [0,0,0,0,0,0,3,-3,0,0,0,0,0,1],
            [0,0,0,0,0,0,3,-3,0,0,0,0,0,2],
            [0,0,0,0,0,0,0,2,0,0,-3,0,0,0],
            [0,0,0,0,0,0,0,-5,13,0,0,0,0,2],
            [0,0,0,0,0,0,0,2,0,-1,0,0,0,0],
            [0,0,0,0,0,0,0,2,0,-1,0,0,0,2],
            [0,0,0,0,0,0,0,2,0,0,-2,0,0,0],
            [0,0,0,0,0,0,0,2,0,0,-2,0,0,1],
            [0,0,0,0,0,0,0,3,-2,0,0,0,0,0],
            [0,0,0,0,0,0,0,3,-2,0,0,0,0,2],
            [0,0,0,0,0,0,0,2,0,0,-1,0,0,2],
            [0,0,0,0,0,0,0,-6,15,0,0,0,0,2],
            [0,0,0,0,0,0,-8,15,0,0,0,0,0,2],
            [0,0,0,0,0,0,-3,9,-4,0,0,0,0,2],
            [0,0,0,0,0,0,0,2,0,2,-5,0,0,2],
            [0,0,0,0,0,0,0,-2,8,-1,-5,0,0,2],
            [0,0,0,0,0,0,0,6,-8,3,0,0,0,2],
            [0,0,0,0,0,0,0,2,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,2,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,2,0,0,0,0,0,1],
            [0,0,1,-1,1,0,0,1,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,2,0,0,0,0,0,1],
            [0,0,0,0,0,0,0,2,0,0,0,0,0,2],
            [0,0,0,0,0,0,0,-6,16,-4,-5,0,0,2],
            [0,0,0,0,0,0,0,-2,8,-3,0,0,0,2],
            [0,0,0,0,0,0,0,-2,8,-3,0,0,0,2],
            [0,0,0,0,0,0,0,6,-8,1,5,0,0,2],
            [0,0,0,0,0,0,0,2,0,-2,5,0,0,2],
            [0,0,0,0,0,0,3,-5,4,0,0,0,0,2],
            [0,0,0,0,0,0,-8,11,0,0,0,0,0,2],
            [0,0,0,0,0,0,-8,11,0,0,0,0,0,1],
            [0,0,0,0,0,0,-8,11,0,0,0,0,0,2],
            [0,0,0,0,0,0,0,11,0,0,0,0,0,2],
            [0,0,0,0,0,0,0,2,0,0,1,0,0,2],
            [0,0,0,0,0,0,3,-3,0,2,0,0,0,2],
            [0,0,2,-2,1,0,0,4,-8,3,0,0,0,0],
            [0,0,1,-1,0,0,0,1,0,0,0,0,0,0],
            [0,0,2,-2,1,0,0,-4,8,-3,0,0,0,0],
            [0,0,0,0,0,0,0,1,2,0,0,0,0,2],
            [0,0,0,0,0,0,0,2,0,1,0,0,0,2],
            [0,0,0,0,0,0,-3,7,0,0,0,0,0,2],
            [0,0,0,0,0,0,0,0,4,0,0,0,0,2],
            [0,0,0,0,0,0,-5,6,0,0,0,0,0,2],
            [0,0,0,0,0,0,-5,6,0,0,0,0,0,1],
            [0,0,0,0,0,0,5,-6,0,0,0,0,0,0],
            [0,0,0,0,0,0,5,-6,0,0,0,0,0,2],
            [0,0,0,0,0,0,0,2,0,2,0,0,0,2],
            [0,0,0,0,0,0,0,-1,6,0,0,0,0,2],
            [0,0,0,0,0,0,0,7,-9,0,0,0,0,2],
            [0,0,0,0,0,0,2,-1,0,0,0,0,0,0],
            [0,0,0,0,0,0,2,-1,0,0,0,0,0,2],
            [0,0,0,0,0,0,0,6,-7,0,0,0,0,2],
            [0,0,0,0,0,0,0,5,-5,0,0,0,0,2],
            [0,0,0,0,0,0,-1,4,0,0,0,0,0,1],
            [0,0,0,0,0,0,-1,4,0,0,0,0,0,2],
            [0,0,0,0,0,0,-7,9,0,0,0,0,0,2],
            [0,0,0,0,0,0,-7,9,0,0,0,0,0,1],
            [0,0,0,0,0,0,0,4,-3,0,0,0,0,2],
            [0,0,0,0,0,0,0,3,-1,0,0,0,0,2],
            [0,0,0,0,0,0,-4,4,0,0,0,0,0,1],
            [0,0,0,0,0,0,4,-4,0,0,0,0,0,0],
            [0,0,0,0,0,0,4,-4,0,0,0,0,0,1],
            [0,0,0,0,0,0,4,-4,0,0,0,0,0,2],
            [0,0,0,0,0,0,0,2,1,0,0,0,0,2],
            [0,0,0,0,0,0,0,-3,0,5,0,0,0,2],
            [0,0,0,0,0,0,1,1,0,0,0,0,0,0],
            [0,0,0,0,0,0,1,1,0,0,0,0,0,1],
            [0,0,0,0,0,0,1,1,0,0,0,0,0,2],
            [0,0,0,0,0,0,-9,12,0,0,0,0,0,2],
            [0,0,0,0,0,0,0,3,0,-4,0,0,0,0],
            [0,0,2,-2,1,0,1,-1,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,7,-8,0,0,0,0,2],
            [0,0,0,0,0,0,0,3,0,-3,0,0,0,0],
            [0,0,0,0,0,0,0,3,0,-3,0,0,0,2],
            [0,0,0,0,0,0,-2,6,0,0,0,0,0,2],
            [0,0,0,0,0,0,-6,7,0,0,0,0,0,1],
            [0,0,0,0,0,0,6,-7,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,6,-6,0,0,0,0,2],
            [0,0,0,0,0,0,0,3,0,-2,0,0,0,0],
            [0,0,0,0,0,0,0,3,0,-2,0,0,0,2],
            [0,0,0,0,0,0,0,5,-4,0,0,0,0,2],
            [0,0,0,0,0,0,3,-2,0,0,0,0,0,0],
            [0,0,0,0,0,0,3,-2,0,0,0,0,0,2],
            [0,0,0,0,0,0,0,3,0,-1,0,0,0,2],
            [0,0,0,0,0,0,0,3,0,-1,0,0,0,2],
            [0,0,0,0,0,0,0,3,0,0,-2,0,0,2],
            [0,0,0,0,0,0,0,4,-2,0,0,0,0,2],
            [0,0,0,0,0,0,0,3,0,0,-1,0,0,2],
            [0,0,2,-2,1,0,0,1,0,-1,0,0,0,0],
            [0,0,0,0,0,0,-8,16,0,0,0,0,0,2],
            [0,0,0,0,0,0,0,3,0,2,-5,0,0,2],
            [0,0,0,0,0,0,0,7,-8,3,0,0,0,2],
            [0,0,0,0,0,0,0,-5,16,-4,-5,0,0,2],
            [0,0,0,0,0,0,0,3,0,0,0,0,0,2],
            [0,0,0,0,0,0,0,-1,8,-3,0,0,0,2],
            [0,0,0,0,0,0,-8,10,0,0,0,0,0,2],
            [0,0,0,0,0,0,-8,10,0,0,0,0,0,1],
            [0,0,0,0,0,0,-8,10,0,0,0,0,0,2],
            [0,0,0,0,0,0,0,2,2,0,0,0,0,2],
            [0,0,0,0,0,0,0,3,0,1,0,0,0,2],
            [0,0,0,0,0,0,-3,8,0,0,0,0,0,2],
            [0,0,0,0,0,0,-5,5,0,0,0,0,0,1],
            [0,0,0,0,0,0,5,-5,0,0,0,0,0,0],
            [0,0,0,0,0,0,5,-5,0,0,0,0,0,1],
            [0,0,0,0,0,0,5,-5,0,0,0,0,0,2],
            [0,0,0,0,0,0,2,0,0,0,0,0,0,0],
            [0,0,0,0,0,0,2,0,0,0,0,0,0,1],
            [0,0,0,0,0,0,2,0,0,0,0,0,0,2],
            [0,0,0,0,0,0,0,7,-7,0,0,0,0,2],
            [0,0,0,0,0,0,0,7,-7,0,0,0,0,2],
            [0,0,0,0,0,0,0,6,-5,0,0,0,0,2],
            [0,0,0,0,0,0,7,-8,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,5,-3,0,0,0,0,2],
            [0,0,0,0,0,0,4,-3,0,0,0,0,0,2],
            [0,0,0,0,0,0,1,2,0,0,0,0,0,2],
            [0,0,0,0,0,0,-9,11,0,0,0,0,0,2],
            [0,0,0,0,0,0,-9,11,0,0,0,0,0,1],
            [0,0,0,0,0,0,0,4,0,-4,0,0,0,2],
            [0,0,0,0,0,0,0,4,0,-3,0,0,0,2],
            [0,0,0,0,0,0,-6,6,0,0,0,0,0,1],
            [0,0,0,0,0,0,6,-6,0,0,0,0,0,0],
            [0,0,0,0,0,0,6,-6,0,0,0,0,0,1],
            [0,0,0,0,0,0,0,4,0,-2,0,0,0,2],
            [0,0,0,0,0,0,0,6,-4,0,0,0,0,2],
            [0,0,0,0,0,0,3,-1,0,0,0,0,0,0],
            [0,0,0,0,0,0,3,-1,0,0,0,0,0,1],
            [0,0,0,0,0,0,3,-1,0,0,0,0,0,2],
            [0,0,0,0,0,0,0,4,0,-1,0,0,0,2],
            [0,0,0,0,0,0,0,4,0,0,-2,0,0,2],
            [0,0,0,0,0,0,0,5,-2,0,0,0,0,2],
            [0,0,0,0,0,0,0,4,0,0,0,0,0,0],
            [0,0,0,0,0,0,8,-9,0,0,0,0,0,0],
            [0,0,0,0,0,0,5,-4,0,0,0,0,0,2],
            [0,0,0,0,0,0,2,1,0,0,0,0,0,2],
            [0,0,0,0,0,0,2,1,0,0,0,0,0,1],
            [0,0,0,0,0,0,2,1,0,0,0,0,0,1],
            [0,0,0,0,0,0,-7,7,0,0,0,0,0,1],
            [0,0,0,0,0,0,7,-7,0,0,0,0,0,0],
            [0,0,0,0,0,0,4,-2,0,0,0,0,0,1],
            [0,0,0,0,0,0,4,-2,0,0,0,0,0,2],
            [0,0,0,0,0,0,4,-2,0,0,0,0,0,0],
            [0,0,0,0,0,0,4,-2,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,5,0,-4,0,0,0,2],
            [0,0,0,0,0,0,0,5,0,-3,0,0,0,2],
            [0,0,0,0,0,0,0,5,0,-2,0,0,0,2],
            [0,0,0,0,0,0,3,0,0,0,0,0,0,2],
            [0,0,0,0,0,0,-8,8,0,0,0,0,0,1],
            [0,0,0,0,0,0,8,-8,0,0,0,0,0,0],
            [0,0,0,0,0,0,5,-3,0,0,0,0,0,1],
            [0,0,0,0,0,0,5,-3,0,0,0,0,0,2],
            [0,0,0,0,0,0,-9,9,0,0,0,0,0,1],
            [0,0,0,0,0,0,-9,9,0,0,0,0,0,1],
            [0,0,0,0,0,0,-9,9,0,0,0,0,0,1],
            [0,0,0,0,0,0,9,-9,0,0,0,0,0,0],
            [0,0,0,0,0,0,6,-4,0,0,0,0,0,1],
            [0,0,0,0,0,0,0,6,0,0,0,0,0,2],
            [0,0,0,0,0,0,0,6,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,6,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,6,0,0,0,0,0,1],
            [0,0,0,0,0,0,0,6,0,0,0,0,0,2],
            [0,0,0,0,0,0,0,6,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,6,0,0,0,0,0,1],
            [0,0,0,0,0,0,0,6,0,0,0,0,0,2],
            [0,0,0,0,0,0,0,0,0,0,0,0,0,2],
            [1,0,0,-2,0,0,0,2,0,-2,0,0,0,0],
            [1,0,0,-2,0,0,2,-2,0,0,0,0,0,0],
            [1,0,0,-2,0,0,0,1,0,-1,0,0,0,0],
            [1,0,0,-2,0,0,1,-1,0,0,0,0,0,0],
            [-1,0,0,0,0,0,3,-3,0,0,0,0,0,0],
            [-1,0,0,0,0,0,0,2,0,-2,0,0,0,0],
            [-1,0,0,2,0,0,0,4,-8,3,0,0,0,0],
            [1,0,0,-2,0,0,0,4,-8,3,0,0,0,0],
            [-2,0,0,2,0,0,0,4,-8,3,0,0,0,0],
            [-1,0,0,0,0,0,0,2,0,-3,0,0,0,0],
            [-1,0,0,0,0,0,0,1,0,-1,0,0,0,0],
            [-1,0,0,0,0,0,1,-1,0,0,0,0,0,0],
            [-1,0,0,2,0,0,2,-2,0,0,0,0,0,0],
            [1,0,-1,1,0,0,0,1,0,0,0,0,0,0],
            [-1,0,0,2,0,0,0,2,0,-3,0,0,0,0],
            [-2,0,0,0,0,0,0,2,0,-3,0,0,0,0],
            [1,0,0,0,0,0,0,4,-8,3,0,0,0,0],
            [-1,0,1,-1,1,0,0,-1,0,0,0,0,0,0],
            [1,0,1,-1,1,0,0,-1,0,0,0,0,0,0],
            [-1,0,0,0,0,0,0,4,-8,3,0,0,0,0],
            [-1,0,0,2,1,0,0,2,0,-2,0,0,0,0],
            [0,0,0,0,0,0,0,2,0,-2,0,0,0,0],
            [-1,0,0,2,0,0,0,2,0,-2,0,0,0,0],
            [-1,0,0,2,0,0,3,-3,0,0,0,0,0,0],
            [1,0,0,-2,1,0,0,-2,0,2,0,0,0,0],
            [1,0,2,-2,2,0,-3,3,0,0,0,0,0,0],
            [1,0,2,-2,2,0,0,-2,0,2,0,0,0,0],
            [1,0,0,0,0,0,1,-1,0,0,0,0,0,0],
            [1,0,0,0,0,0,0,1,0,-1,0,0,0,0],
            [0,0,0,-2,0,0,2,-2,0,0,0,0,0,0],
            [0,0,0,-2,0,0,0,1,0,-1,0,0,0,0],
            [0,0,2,0,2,0,-2,2,0,0,0,0,0,0],
            [0,0,2,0,2,0,0,-1,0,1,0,0,0,0],
            [0,0,2,0,2,0,-1,1,0,0,0,0,0,0],
            [0,0,2,0,2,0,-2,3,0,0,0,0,0,0],
            [0,0,0,2,0,0,0,2,0,-2,0,0,0,0],
            [0,0,1,1,2,0,0,1,0,0,0,0,0,0],
            [1,0,2,0,2,0,0,1,0,0,0,0,0,0],
            [-1,0,2,0,2,0,10,-3,0,0,0,0,0,0],
            [0,0,1,1,1,0,0,1,0,0,0,0,0,0],
            [1,0,2,0,2,0,0,1,0,0,0,0,0,0],
            [0,0,2,0,2,0,0,4,-8,3,0,0,0,0],
            [0,0,2,0,2,0,0,-4,8,-3,0,0,0,0],
            [-1,0,2,0,2,0,0,-4,8,-3,0,0,0,0],
            [2,0,2,-2,2,0,0,-2,0,3,0,0,0,0],
            [1,0,2,0,1,0,0,-2,0,3,0,0,0,0],
            [0,0,1,1,0,0,0,1,0,0,0,0,0,0],
            [-1,0,2,0,1,0,0,1,0,0,0,0,0,0],
            [-2,0,2,2,2,0,0,2,0,-2,0,0,0,0],
            [0,0,2,0,2,0,2,-3,0,0,0,0,0,0],
            [0,0,2,0,2,0,1,-1,0,0,0,0,0,0],
            [0,0,2,0,2,0,0,1,0,-1,0,0,0,0],
            [0,0,2,0,2,0,2,-2,0,0,0,0,0,0],
            [-1,0,2,2,2,0,0,-1,0,1,0,0,0,0],
            [1,0,2,0,2,0,-1,1,0,0,0,0,0,0],
            [-1,0,2,2,2,0,0,2,0,-3,0,0,0,0],
            [2,0,2,0,2,0,0,2,0,-3,0,0,0,0],
            [1,0,2,0,2,0,0,-4,8,-3,0,0,0,0],
            [1,0,2,0,2,0,0,4,-8,3,0,0,0,0],
            [1,0,1,1,1,0,0,1,0,0,0,0,0,0],
            [0,0,2,0,2,0,0,1,0,0,0,0,0,0],
            [2,0,2,0,1,0,0,1,0,0,0,0,0,0],
            [-1,0,2,2,2,0,0,2,0,-2,0,0,0,0],
            [-1,0,2,2,2,0,3,-3,0,0,0,0,0,0],
            [1,0,2,0,2,0,1,-1,0,0,0,0,0,0],
            [0,0,2,2,2,0,0,2,0,-2,0,0,0,0]]
    
    #行星章动系数，单位：0.1微角秒
    #*  longitude (sin, cos), obliquity (sin, cos)
    ICPL=[[1440,0,0,0],
            [56,-117,-42,-40],
            [125,-43,0,-54],
            [0,5,0,0],
            [3,-7,-3,0],
            [3,0,0,-2],
            [-114,0,0,61],
            [-219,89,0,0],
            [-3,0,0,0],
            [-462,1604,0,0],
            [99,0,0,-53],
            [-3,0,0,2],
            [0,6,2,0],
            [3,0,0,0],
            [-12,0,0,0],
            [14,-218,117,8],
            [31,-481,-257,-17],
            [-491,128,0,0],
            [-3084,5123,2735,1647],
            [-1444,2409,-1286,-771],
            [11,-24,-11,-9],
            [26,-9,0,0],
            [103,-60,0,0],
            [0,-13,-7,0],
            [-26,-29,-16,14],
            [9,-27,-14,-5],
            [12,0,0,-6],
            [-7,0,0,0],
            [0,24,0,0],
            [284,0,0,-151],
            [226,101,0,0],
            [0,-8,-2,0],
            [0,-6,-3,0],
            [5,0,0,-3],
            [-41,175,76,17],
            [0,15,6,0],
            [425,212,-133,269],
            [1200,598,319,-641],
            [235,334,0,0],
            [11,-12,-7,-6],
            [5,-6,3,3],
            [-5,0,0,3],
            [6,0,0,-3],
            [15,0,0,0],
            [13,0,0,-7],
            [-6,-9,0,0],
            [266,-78,0,0],
            [-460,-435,-232,246],
            [0,15,7,0],
            [-3,0,0,2],
            [0,131,0,0],
            [4,0,0,0],
            [0,3,0,0],
            [0,4,2,0],
            [0,3,0,0],
            [-17,-19,-10,9],
            [-9,-11,6,-5],
            [-6,0,0,3],
            [-16,8,0,0],
            [0,3,0,0],
            [11,24,11,-5],
            [-3,-4,-2,1],
            [3,0,0,-1],
            [0,-8,-4,0],
            [0,3,0,0],
            [0,5,0,0],
            [0,3,2,0],
            [-6,4,2,3],
            [-3,-5,0,0],
            [-5,0,0,2],
            [4,24,13,-2],
            [-42,20,0,0],
            [-10,233,0,0],
            [-3,0,0,1],
            [78,-18,0,0],
            [0,3,1,0],
            [0,-3,-1,0],
            [0,-4,-2,1],
            [0,-8,-4,-1],
            [0,-5,3,0],
            [-7,0,0,3],
            [-14,8,3,6],
            [0,8,-4,0],
            [0,19,10,0],
            [45,-22,0,0],
            [-3,0,0,0],
            [0,-3,0,0],
            [0,3,0,0],
            [3,5,3,-2],
            [89,-16,-9,-48],
            [0,3,0,0],
            [-3,7,4,2],
            [-349,-62,0,0],
            [-15,22,0,0],
            [-3,0,0,0],
            [-53,0,0,0],
            [5,0,0,-3],
            [0,-8,0,0],
            [15,-7,-4,-8],
            [-3,0,0,1],
            [-21,-78,0,0],
            [20,-70,-37,-11],
            [0,6,3,0],
            [5,3,2,-2],
            [-17,-4,-2,9],
            [0,6,3,0],
            [32,15,-8,17],
            [174,84,45,-93],
            [11,56,0,0],
            [-66,-12,-6,35],
            [47,8,4,-25],
            [0,8,4,0],
            [10,-22,-12,-5],
            [-3,0,0,2],
            [-24,12,0,0],
            [5,-6,0,0],
            [3,0,0,-2],
            [4,3,1,-2],
            [0,29,15,0],
            [-5,-4,-2,2],
            [8,-3,-1,-5],
            [0,-3,0,0],
            [10,0,0,0],
            [3,0,0,-2],
            [-5,0,0,3],
            [46,66,35,-25],
            [-14,7,0,0],
            [0,3,2,0],
            [-5,0,0,0],
            [-68,-34,-18,36],
            [0,14,7,0],
            [10,-6,-3,-5],
            [-5,-4,-2,3],
            [-3,5,2,1],
            [76,17,9,-41],
            [84,298,159,-45],
            [3,0,0,-1],
            [-3,0,0,2],
            [-3,0,0,1],
            [-82,292,156,44],
            [-73,17,9,39],
            [-9,-16,0,0],
            [3,0,-1,-2],
            [-3,0,0,0],
            [-9,-5,-3,5],
            [-439,0,0,0],
            [57,-28,-15,-30],
            [0,-6,-3,0],
            [-4,0,0,2],
            [-40,57,30,21],
            [23,7,3,-13],
            [273,80,43,-146],
            [-449,430,0,0],
            [-8,-47,-25,4],
            [6,47,25,-3],
            [0,23,13,0],
            [-3,0,0,2],
            [3,-4,-2,-2],
            [-48,-110,-59,26],
            [51,114,61,-27],
            [-133,0,0,57],
            [0,4,0,0],
            [-21,-6,-3,11],
            [0,-3,-1,0],
            [-11,-21,-11,6],
            [-18,-436,-233,9],
            [35,-7,0,0],
            [0,5,3,0],
            [11,-3,-1,-6],
            [-5,-3,-1,3],
            [-53,-9,-5,28],
            [0,3,2,1],
            [4,0,0,-2],
            [0,-4,0,0],
            [-50,194,103,27],
            [-13,52,28,7],
            [-91,248,0,0],
            [6,49,26,-3],
            [-6,-47,-25,3],
            [0,5,3,0],
            [52,23,10,-23],
            [-3,0,0,1],
            [0,5,3,0],
            [-4,0,0,0],
            [-4,8,3,2],
            [10,0,0,0],
            [3,0,0,-2],
            [0,8,4,0],
            [0,8,4,1],
            [-4,0,0,0],
            [-4,0,0,0],
            [-8,4,2,4],
            [8,-4,-2,-4],
            [0,15,7,0],
            [-138,0,0,0],
            [0,-7,-3,0],
            [0,-7,-3,0],
            [54,0,0,-29],
            [0,10,4,0],
            [-7,0,0,3],
            [-37,35,19,20],
            [0,4,0,0],
            [-4,9,0,0],
            [8,0,0,-4],
            [-9,-14,-8,5],
            [-3,-9,-5,3],
            [-145,47,0,0],
            [-10,40,21,5],
            [11,-49,-26,-7],
            [-2150,0,0,932],
            [-12,0,0,5],
            [85,0,0,-37],
            [4,0,0,-2],
            [3,0,0,-2],
            [-86,153,0,0],
            [-6,9,5,3],
            [9,-13,-7,-5],
            [-8,12,6,4],
            [-51,0,0,22],
            [-11,-268,-116,5],
            [0,12,5,0],
            [0,7,3,0],
            [31,6,3,-17],
            [140,27,14,-75],
            [57,11,6,-30],
            [-14,-39,0,0],
            [0,-6,-2,0],
            [4,15,8,-2],
            [0,4,0,0],
            [-3,0,0,1],
            [0,11,5,0],
            [9,6,0,0],
            [-4,10,4,2],
            [5,3,0,0],
            [16,0,0,-9],
            [-3,0,0,0],
            [0,3,2,-1],
            [7,0,0,-3],
            [-25,22,0,0],
            [42,223,119,-22],
            [-27,-143,-77,14],
            [9,49,26,-5],
            [-1166,0,0,505],
            [-5,0,0,2],
            [-6,0,0,3],
            [-8,0,1,4],
            [0,-4,0,0],
            [117,0,0,-63],
            [-4,8,4,2],
            [3,0,0,-2],
            [-5,0,0,2],
            [0,31,0,0],
            [-5,0,1,3],
            [4,0,0,-2],
            [-4,0,0,2],
            [-24,-13,-6,10],
            [3,0,0,0],
            [0,-32,-17,0],
            [8,12,5,-3],
            [3,0,0,-1],
            [7,13,0,0],
            [-3,16,0,0],
            [50,0,0,-27],
            [0,-5,-3,0],
            [13,0,0,0],
            [0,5,3,1],
            [24,5,2,-11],
            [5,-11,-5,-2],
            [30,-3,-2,-16],
            [18,0,0,-9],
            [8,614,0,0],
            [3,-3,-1,-2],
            [6,17,9,-3],
            [-3,-9,-5,2],
            [0,6,3,-1],
            [-127,21,9,55],
            [3,5,0,0],
            [-6,-10,-4,3],
            [5,0,0,0],
            [16,9,4,-7],
            [3,0,0,-2],
            [0,22,0,0],
            [0,19,10,0],
            [7,0,0,-4],
            [0,-5,-2,0],
            [0,3,1,0],
            [-9,3,1,4],
            [17,0,0,-7],
            [0,-3,-2,-1],
            [-20,34,0,0],
            [-10,0,1,5],
            [-4,0,0,2],
            [22,-87,0,0],
            [-4,0,0,2],
            [-3,-6,-2,1],
            [-16,-3,-1,7],
            [0,-3,-2,0],
            [4,0,0,0],
            [-68,39,0,0],
            [27,0,0,-14],
            [0,-4,0,0],
            [-25,0,0,0],
            [-12,-3,-2,6],
            [3,0,0,-1],
            [3,66,29,-1],
            [490,0,0,-213],
            [-22,93,49,12],
            [-7,28,15,4],
            [-3,13,7,2],
            [-46,14,0,0],
            [-5,0,0,0],
            [2,1,0,0],
            [0,-3,0,0],
            [-28,0,0,15],
            [5,0,0,-2],
            [0,3,0,0],
            [-11,0,0,5],
            [0,3,1,0],
            [-3,0,0,1],
            [25,106,57,-13],
            [5,21,11,-3],
            [1485,0,0,0],
            [-7,-32,-17,4],
            [0,5,3,0],
            [-6,-3,-2,3],
            [30,-6,-2,-13],
            [-4,4,0,0],
            [-19,0,0,10],
            [0,4,2,-1],
            [0,3,0,0],
            [4,0,0,-2],
            [0,-3,-1,0],
            [-3,0,0,0],
            [5,3,1,-2],
            [0,11,0,0],
            [118,0,0,-52],
            [0,-5,-3,0],
            [-28,36,0,0],
            [5,-5,0,0],
            [14,-59,-31,-8],
            [0,9,5,1],
            [-458,0,0,198],
            [0,-45,-20,0],
            [9,0,0,-5],
            [0,-3,0,0],
            [0,-4,-2,-1],
            [11,0,0,-6],
            [6,0,0,-2],
            [-16,23,0,0],
            [0,-4,-2,0],
            [-5,0,0,2],
            [-166,269,0,0],
            [15,0,0,-8],
            [10,0,0,-4],
            [-78,45,0,0],
            [0,-5,-2,0],
            [7,0,0,-4],
            [-5,328,0,0],
            [3,0,0,-2],
            [5,0,0,-2],
            [0,3,1,0],
            [-3,0,0,0],
            [-3,0,0,0],
            [0,-4,-2,0],
            [-1223,-26,0,0],
            [0,7,3,0],
            [3,0,0,0],
            [0,3,2,0],
            [-6,20,0,0],
            [-368,0,0,0],
            [-75,0,0,0],
            [11,0,0,-6],
            [3,0,0,-2],
            [-3,0,0,1],
            [-13,-30,0,0],
            [21,3,0,0],
            [-3,0,0,1],
            [-4,0,0,2],
            [8,-27,0,0],
            [-19,-11,0,0],
            [-4,0,0,2],
            [0,5,2,0],
            [-6,0,0,2],
            [-8,0,0,0],
            [-1,0,0,0],
            [-14,0,0,6],
            [6,0,0,0],
            [-74,0,0,32],
            [0,-3,-1,0],
            [4,0,0,-2],
            [8,11,0,0],
            [0,3,2,0],
            [-262,0,0,114],
            [0,-4,0,0],
            [-7,0,0,4],
            [0,-27,-12,0],
            [-19,-8,-4,8],
            [202,0,0,-87],
            [-8,35,19,5],
            [0,4,2,0],
            [16,-5,0,0],
            [5,0,0,-3],
            [0,-3,0,0],
            [1,0,0,0],
            [-35,-48,-21,15],
            [-3,-5,-2,1],
            [6,0,0,-3],
            [3,0,0,-1],
            [0,-5,0,0],
            [12,55,29,-6],
            [0,5,3,0],
            [-598,0,0,0],
            [-3,-13,-7,1],
            [-5,-7,-3,2],
            [3,0,0,-1],
            [5,-7,0,0],
            [4,0,0,-2],
            [16,-6,0,0],
            [8,-3,0,0],
            [8,-31,-16,-4],
            [0,3,1,0],
            [113,0,0,-49],
            [0,-24,-10,0],
            [4,0,0,-2],
            [27,0,0,0],
            [-3,0,0,1],
            [0,-4,-2,0],
            [5,0,0,-2],
            [0,-3,0,0],
            [-13,0,0,6],
            [5,0,0,-2],
            [-18,-10,-4,8],
            [-4,-28,0,0],
            [-5,6,3,2],
            [-3,0,0,1],
            [-5,-9,-4,2],
            [17,0,0,-7],
            [11,4,0,0],
            [0,-6,-2,0],
            [83,15,0,0],
            [-4,0,0,2],
            [0,-114,-49,0],
            [117,0,0,-51],
            [-5,19,10,2],
            [-3,0,0,0],
            [-3,0,0,2],
            [0,-3,-1,0],
            [3,0,0,0],
            [0,-6,-2,0],
            [393,3,0,0],
            [-4,21,11,2],
            [-6,0,-1,3],
            [-3,8,4,1],
            [8,0,0,0],
            [18,-29,-13,-8],
            [8,34,18,-4],
            [89,0,0,0],
            [3,12,6,-1],
            [54,-15,-7,-24],
            [0,3,0,0],
            [3,0,0,-1],
            [0,35,0,0],
            [-154,-30,-13,67],
            [15,0,0,0],
            [0,4,2,0],
            [0,9,0,0],
            [80,-71,-31,-35],
            [0,-20,-9,0],
            [11,5,2,-5],
            [61,-96,-42,-27],
            [14,9,4,-6],
            [-11,-6,-3,5],
            [0,-3,-1,0],
            [123,-415,-180,-53],
            [0,0,0,-35],
            [-5,0,0,0],
            [7,-32,-17,-4],
            [0,-9,-5,0],
            [0,-4,2,0],
            [-89,0,0,38],
            [0,-86,-19,-6],
            [0,0,-19,6],
            [-123,-416,-180,53],
            [0,-3,-1,0],
            [12,-6,-3,-5],
            [-13,9,4,6],
            [0,-15,-7,0],
            [3,0,0,-1],
            [-62,-97,-42,27],
            [-11,5,2,5],
            [0,-19,-8,0],
            [-3,0,0,1],
            [0,4,2,0],
            [0,3,0,0],
            [0,4,2,0],
            [-85,-70,-31,37],
            [163,-12,-5,-72],
            [-63,-16,-7,28],
            [-21,-32,-14,9],
            [0,-3,-1,0],
            [3,0,0,-2],
            [0,8,0,0],
            [3,10,4,-1],
            [3,0,0,-1],
            [0,-7,-3,0],
            [0,-4,-2,0],
            [6,19,0,0],
            [5,-173,-75,-2],
            [0,-7,-3,0],
            [7,-12,-5,-3],
            [-3,0,0,2],
            [3,-4,-2,-1],
            [74,0,0,-32],
            [-3,12,6,2],
            [26,-14,-6,-11],
            [19,0,0,-8],
            [6,24,13,-3],
            [83,0,0,0],
            [0,-10,-5,0],
            [11,-3,-1,-5],
            [3,0,1,-1],
            [3,0,0,-1],
            [-4,0,0,0],
            [5,-23,-12,-3],
            [-339,0,0,147],
            [0,-10,-5,0],
            [5,0,0,0],
            [3,0,0,-1],
            [0,-4,-2,0],
            [18,-3,0,0],
            [9,-11,-5,-4],
            [-8,0,0,4],
            [3,0,0,-1],
            [0,9,0,0],
            [6,-9,-4,-2],
            [-4,-12,0,0],
            [67,-91,-39,-29],
            [30,-18,-8,-13],
            [0,0,0,0],
            [0,-114,-50,0],
            [0,0,0,23],
            [517,16,7,-224],
            [0,-7,-3,0],
            [143,-3,-1,-62],
            [29,0,0,-13],
            [-4,0,0,2],
            [-6,0,0,3],
            [5,12,5,-2],
            [-25,0,0,11],
            [-3,0,0,1],
            [0,4,2,0],
            [-22,12,5,10],
            [50,0,0,-22],
            [0,7,4,0],
            [0,3,1,0],
            [-4,4,2,2],
            [-5,-11,-5,2],
            [0,4,2,0],
            [4,17,9,-2],
            [59,0,0,0],
            [0,-4,-2,0],
            [-8,0,0,4],
            [-3,0,0,0],
            [4,-15,-8,-2],
            [370,-8,0,-160],
            [0,0,-3,0],
            [0,3,1,0],
            [-6,3,1,3],
            [0,6,0,0],
            [-10,0,0,4],
            [0,9,4,0],
            [4,17,7,-2],
            [34,0,0,-15],
            [0,5,3,0],
            [-5,0,0,2],
            [-37,-7,-3,16],
            [3,13,7,-2],
            [40,0,0,0],
            [0,-3,-2,0],
            [-184,-3,-1,80],
            [-3,0,0,1],
            [-3,0,0,0],
            [0,-10,-6,-1],
            [31,-6,0,-13],
            [-3,-32,-14,1],
            [-7,0,0,3],
            [0,-8,-4,0],
            [3,-4,0,0],
            [0,4,0,0],
            [0,3,1,0],
            [19,-23,-10,2],
            [0,0,0,-10],
            [0,3,2,0],
            [0,9,5,-1],
            [28,0,0,0],
            [0,-7,-4,0],
            [8,-4,0,-4],
            [0,0,-2,0],
            [0,3,0,0],
            [-3,0,0,1],
            [-9,0,1,4],
            [3,12,5,-1],
            [17,-3,-1,0],
            [0,7,4,0],
            [19,0,0,0],
            [0,-5,-3,0],
            [14,-3,0,-1],
            [0,0,-1,0],
            [0,0,0,-5],
            [0,5,3,0],
            [13,0,0,0],
            [0,-3,-2,0],
            [2,9,4,3],
            [0,0,0,-4],
            [8,0,0,0],
            [0,4,2,0],
            [6,0,0,-3],
            [6,0,0,0],
            [0,3,1,0],
            [5,0,0,-2],
            [3,0,0,-1],
            [-3,0,0,0],
            [6,0,0,0],
            [7,0,0,0],
            [-4,0,0,0],
            [4,0,0,0],
            [6,0,0,0],
            [0,-4,0,0],
            [0,-4,0,0],
            [5,0,0,0],
            [-3,0,0,0],
            [4,0,0,0],
            [-5,0,0,0],
            [4,0,0,0],
            [0,3,0,0],
            [13,0,0,0],
            [21,11,0,0],
            [0,-5,0,0],
            [0,-5,-2,0],
            [0,5,3,0],
            [0,-5,0,0],
            [-3,0,0,2],
            [20,10,0,0],
            [-34,0,0,0],
            [-19,0,0,0],
            [3,0,0,-2],
            [-3,0,0,1],
            [-6,0,0,3],
            [-4,0,0,0],
            [3,0,0,0],
            [3,0,0,0],
            [4,0,0,0],
            [3,0,0,-1],
            [6,0,0,-3],
            [-8,0,0,3],
            [0,3,1,0],
            [-3,0,0,0],
            [0,-3,-2,0],
            [126,-63,-27,-55],
            [-5,0,1,2],
            [-3,28,15,2],
            [5,0,1,-2],
            [0,9,4,1],
            [0,9,4,-1],
            [-126,-63,-27,55],
            [3,0,0,-1],
            [21,-11,-6,-11],
            [0,-4,0,0],
            [-21,-11,-6,11],
            [-3,0,0,1],
            [0,3,1,0],
            [8,0,0,-4],
            [-6,0,0,3],
            [-3,0,0,1],
            [3,0,0,-1],
            [-3,0,0,1],
            [-5,0,0,2],
            [24,-12,-5,-11],
            [0,3,1,0],
            [0,3,1,0],
            [0,3,2,0],
            [-24,-12,-5,10],
            [4,0,-1,-2],
            [13,0,0,-6],
            [7,0,0,-3],
            [3,0,0,-1],
            [3,0,0,-1]]
    
    #给定日期与参考日期之间的时间间隔，单位：儒略世纪数
    T=((DATE1-DJ00)+DATE2)/DJC
    
    #日月章动
    #基本参数（Delaunay）

    #月亮的平均近点角
    EL=pymFal03(T)
    
    #太阳的平均近点角
    A=1287104.79305+T*(129596581.0481+T*(-0.5532+\
        T*(0.000136+T*(-0.00001149))))
    ELP=(A%TURNAS)*DAS2R
    
    #月亮的平黄经减去升交点的平黄经
    F=pymFaf03(T)
    
    #月亮相对于太阳的平均距角
    B=1072260.70369+T*(1602961601.2090+T*(-6.3706+\
        T*(0.006593+T*(-0.00003169))))
    D=(B%TURNAS)*DAS2R
    
    #月亮升交点的平黄经
    OM=pymFaom03(T)
    
    #初始化章动量
    DP=0.0
    DE=0.0
    
    #从后到前对日月章动的各项求和
    for I in range(NLS-1,-1,-1):
        C=float(NALS[I][0])*EL+float(NALS[I][1])*ELP+\
            float(NALS[I][2])*F+float(NALS[I][3])*D+float(NALS[I][4])*OM
        ARG=C%D2PI
        SARG=ma.sin(ARG)
        CARG=ma.cos(ARG)
    
        #Term.
        DP=DP+(CLS[I][0]+CLS[I][1]*T)*SARG+CLS[I][2]*CARG
        DE=DE+(CLS[I][3]+CLS[I][4]*T)*CARG+CLS[I][5]*SARG
    
    #从0.1微弧秒单位转换为弧度。
    DPSILS=DP*U2R
    DEPSLS=DE*U2R
    
    #行星章动

    #月亮的平近点角(MHB2000).
    AL=(2.35555598+8328.6914269554*T)%D2PI
    
    #太阳的平近点角(MHB2000).
    ALSU=(6.24006013+628.301955*T)%D2PI
    
    #月亮的平黄经减去升交点的平黄经(MHB2000).
    AF=(1.627905234+8433.466158131*T)%D2PI

    #月亮相对于太阳的平均距角(MHB2000).
    AD=(5.198466741+7771.3771468121*T)%D2PI

    #月亮升交点的平黄经(MHB2000).
    AOM=(2.18243920-33.757045*T)%D2PI

    #黄经的一般累计进动(IERS 2003).
    APA=pymFapa03(T)
    
    #行星的平黄经，水星到天王星(IERS 2003).
    ALME=pymFame03(T)
    ALVE=pymFave03(T)
    ALEA=pymFae03(T)
    ALMA=pymFama03(T)
    ALJU=pymFaju03(T)
    ALSA=pymFasa03(T)
    ALUR=pymFaur03(T)
    
    #海王星的平黄经(MHB2000).
    ALNE=(5.321159000+3.8127774000*T)%D2PI
    
    #初始化章动量
    DP=0.0
    DE=0.0
    
    #从后到前对行星章动的各项求和
    for I in range(NPL-1,-1,-1):
        D=float(NAPL[I][0])*AL+float(NAPL[I][1])*ALSU+float(NAPL[I][2])*AF+\
            float(NAPL[I][3])*AD+float(NAPL[I][4])*AOM+float(NAPL[I][5])*ALME+\
            float(NAPL[I][6])*ALVE+float(NAPL[I][7])*ALEA+\
            float(NAPL[I][8])*ALMA+float(NAPL[I][9])*ALJU+\
            float(NAPL[I][10])*ALSA+float(NAPL[I][11])*ALUR+\
            float(NAPL[I][12])*ALNE+float(NAPL[I][13])*APA
        ARG=D%D2PI
        SARG=ma.sin(ARG)
        CARG=ma.cos(ARG)
    
        #Term.
        DP=DP+float(ICPL[I][0])*SARG+float(ICPL[I][1])*CARG
        DE=DE+float(ICPL[I][2])*SARG+float(ICPL[I][3])*CARG
    
    #从0.1微弧秒单位转换为弧度。
    DPSIPL=DP*U2R
    DEPSPL=DE*U2R

    #日月章动和行星章动相加
    DPSI=DPSILS+DPSIPL
    DEPS=DEPSLS+DEPSPL
    
    return(DPSI,DEPS)


def pymNut06a(DATE1,DATE2):
    '''
    IAU 2000A nutation with adjustments to match the IAU 2006
    precession.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    dpsi : float
        nutation, luni-solar + planetary    
    deps : float
        nutation, luni-solar + planetary

    '''
    #参考历元(J2000.0), JD
    DJ00=2451545.0
    
    #1儒略世纪的天数
    DJC=36525.0
    
    #相对于参考历元的时间间隔，单位：儒略世纪数.
    T=((DATE1-DJ00)+DATE2)/DJC
    
    #根据IAU 2006，对于J2的长期变化进行改正
    FJ2=-2.7774e-6*T
    
    #调用pymNut00a 获得IAU 2000A 的章动.
    DP,DE=pymNut00a(DATE1, DATE2)
    
    #应用改正项(Wallace & Capitaine, 2006, Eqs.5).
    DPSI=DP+DP*(0.4697e-6+FJ2)
    DEPS=DE+DE*FJ2
    
    return(DPSI,DEPS)


def pymObl06(DATE1,DATE2):
    '''
    Mean obliquity of the ecliptic, IAU 2006 precession model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    function value : float
        obliquity of the ecliptic (radians)

    '''
    #1角秒对应的弧度
    DAS2R=4.848136811095359935899141e-6
    
    #参考历元(J2000.0), JD
    DJ00=2451545.0
    
    #1儒略世纪的天数
    DJC=36525.0
    
    #相对于参考历元的时间间隔，单位：儒略世纪数.
    T=((DATE1-DJ00)+DATE2)/DJC
    
    #黄道平均倾角.
    OBL=(84381.406+(-46.836769+(-0.0001831+(0.00200340+(-0.000000576+\
            (-0.0000000434)*T)*T)*T)*T)*T)*DAS2R

    return(OBL)


def pymPfw06(DATE1,DATE2):
    '''
    Precession angles, IAU 2006 (Fukushima-Williams 4-angle formulation).

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    gamb : flaot
        F-W angle gamma_bar (radians)    
    phib : float
        F-W angle phi_bar (radians)    
    psib : float
        F-W angle psi_bar (radians)    
    epsa : float
        F-W angle epsilon_A (radians)

    '''
    #1角秒对应的弧度
    DAS2R=4.848136811095359935899141e-6
    
    #参考历元(J2000.0), JD
    DJ00=2451545.0
    
    #1儒略世纪的天数
    DJC=36525.0
    
    #相对于参考历元的时间间隔，单位：儒略世纪数.
    T=((DATE1-DJ00)+DATE2)/DJC
    
    #P03极移+进动角.
    GAMB=(-0.052928+(10.556378+(0.4932044+(-0.00031238+\
           (-0.000002788+(0.0000000260)*T)*T)*T)*T)*T)*DAS2R
    PHIB=(84381.412819+(-46.811016+(0.0511268+(0.00053289+\
           (-0.000000440+(-0.0000000176)*T)*T)*T)*T)*T)*DAS2R
    PSIB=(-0.041775+(5038.481484+(1.5584175+(-0.00018522+\
           (-0.000026452+(-0.0000000148)*T)*T)*T)*T)*T)*DAS2R      
    EPSA=pymObl06(DATE1,DATE2)
    
    return(GAMB,PHIB,PSIB,EPSA)


def pymPnm06a(DATE1,DATE2):
    '''
    Form the matrix of precession-nutation for a given date (including
    frame bias), equinox based, IAU 2006 precession and IAU 2000A
    nutation models.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    rbpn : list(3,3)
        bias-precession-nutation matrix

    '''
    #调用pymPfw06获得四个进动角
    GAMB,PHIB,PSIB,EPSA=pymPfw06(DATE1, DATE2)
    
    #调用pymNut06a获得章动参数
    DP,DE=pymNut06a(DATE1,DATE2)
    
    #结合进动角，章动参数，获得用于计算的极移进动章动矩阵
    RBPN=pymFw2m(GAMB,PHIB,PSIB+DP,EPSA+DE)
    
    return(RBPN)


def pymRefco(PHPA,TC,RH,WL):
    '''
    Determine the constants A and B in the atmospheric refraction model
    dZ = A tan Z + B tan^3 Z.

    Parameters
    ----------
    phpa : float
        pressure at the observer (hPa = millibar)    
    tc : float
        ambient temperature at the observer (deg C)     
    rh : float
        relative humidity at the observer (range 0-1)    
    wl : float
        wavelength (micrometers)

    Returns
    -------
    refa : float
        tan Z coefficient (radians)     
    refb : float
        tan^3 Z coefficient (radians)    

    '''
    #判断观测的波长是属于光学近红外还是射电，判断依据为100微米
    OPTIC=(WL<=100.0)
    
    #将输入参数保持在合理范围内
    T=min(max(TC,-150),200)
    P=min(max(PHPA,0),10000)
    R=min(max(RH,0),1)
    W=min(max(WL,0.1),1e6)
    
    #观察者处的水汽压力
    if (P>0.0):
        A=(0.7859+0.03477*T)/(1.0+0.00412*T)
        PS=10.0**A*(1.0+P*(4.5e-6+6e-10*T*T))
        PW=R*PS/(1.0-(1.0-R)*PS/P)
    else:
        PW=0.0
        
    #gamma为公式中的n0-1
    TK=T+273.15
    if (OPTIC):
        WLSQ=W*W
        GAMMA=((77.53484e-6+(4.39108e-7+3.666e-9/WLSQ)/WLSQ)*P\
               -11.2684e-6*PW)/TK
    else:
        GAMMA=(77.6890e-6*P-(6.3938e-6-0.375463/TK)*PW)/TK
    
    #对于BETA的公式是改编自(Stone, Ronald C., P.A.S.P. 108, 1051-1058, 1996.)
    #并且其中根据经验进行了调整
    BETA=4.4474e-6*TK
    if (not OPTIC):
        BETA=BETA-0.0074*PW*BETA
    
    #根据 Green 给出折射常数.
    REFA=GAMMA*(1.0-BETA)
    REFB=-GAMMA*(BETA-GAMMA/2.0)    
    
    return(REFA,REFB)


def pymEra00(DJ1,DJ2):
    '''
    Earth rotation angle (IAU 2000 model).

    Parameters
    ----------
    dj1 : float
        UT1 as a 2-part Julian Date    
    dj2 : float
        UT1 as a 2-part Julian Date

    Returns
    -------
    function value : float
        Earth rotation angle (radians), range 0-2pi

    '''
    #参考历元(J2000.0), JD
    DJ00=2451545.0
    
    #自参考历元经过的天数
    if (DJ1<=DJ2):
        D1=DJ1
        D2=DJ2
    else:
        D1=DJ2
        D2=DJ1
    T=D1+(D2-DJ00)

    #时间的小数部分(days).
    F=D1%1.0+D2%1.0
    
    #在这个UT1时刻的地球自转角.
    A=2*ma.pi*(F+0.7790572732640+0.00273781191135448*T)
    ERA=pymAnp(A)
    
    return(ERA)


def pymSp00(DATE1,DATE2):
    '''
    The TIO locator s', positioning the Terrestrial Intermediate Origin
    on the equator of the Celestial Intermediate Pole.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    function value : float
        the TIO locator s' in radians

    '''
    #1角秒对应的弧度
    DAS2R=4.848136811095359935899141e-6
    
    #参考历元(J2000.0), JD
    DJ00=2451545.0
    
    #1儒略世纪的天数
    DJC=36525.0
    
    #相对于J2000.0的儒略世纪数
    T=((DATE1-DJ00)+DATE2)/DJC
    
    #近似TIO定位角
    SP=-47e-6*T*DAS2R
    
    return(SP)


def pymCr(R):
    '''
    Copy an r-matrix.

    Parameters
    ----------
    r : list(3,3)
        r-matrix to be copied

    Returns
    -------
    c : list(3,3)
        r-matrix to be copied

    '''
    C=[0,0,0]
    
    for k in range(3):
        C[k]=pymCp(R[k])
    
    return(C)


def pymTr(R):
    '''
    Transpose an r-matrix.

    Parameters
    ----------
    r : list(3,3)
        r-matrix

    Returns
    -------
    rt : list(3,3)
        transpose

    '''
    WM=[[0 for i in range(3)] for j in range(3)]
    
    for i in range(3):
        for j in range(3):
            WM[i][j]=R[j][i]
    
    RT=pymCr(WM)    
    
    return(RT)


def pymRxp(R,P):
    '''
    Multiply a p-vector by an r-matrix.

    Parameters
    ----------
    r : list(3,3)
        r-matrix    
    p : list(3)
        p-vector    

    Returns
    -------
    rp : list(3,3)
        r * p

    '''
    WRP=[0,0,0]

    for j in range(3):
        W=0.0
        for i in range(3):
            W=W+R[j][i]*P[i]
        WRP[j]=W
    RP=pymCp(WRP)   
    
    return(RP)


def pymTrxp(R,P):
    '''
    Multiply a p-vector by the transpose of an r-matrix.

    Parameters
    ----------
    r : list(3,3)
        r-matrix    
    p : list(3)
        p-vector

    Returns
    -------
    trp : list(3)
        r^T * p

    '''
    RI=pymTr(R)

    TRP=pymRxp(RI,P)
    
    return(TRP)


def pymPom00(XP,YP,SP):
    '''
    Form the matrix of polar motion for a given date, IAU 2000.

    Parameters
    ----------
    xp : float
        coordinates of the pole (radians)    
    yp : float
        coordinates of the pole (radians)    
    sp : float
        the TIO locator s' (radians)

    Returns
    -------
    rpom : list(3,3)
        polar-motion matrix

    '''
    RPOM=pymIr()
    RPOM=pymRz(SP,RPOM)
    RPOM=pymRy(-XP,RPOM)
    RPOM=pymRx(-YP,RPOM)
    
    return(RPOM)


def pymPvtob(ELONG,PHI,HM,XP,YP,SP,THETA):
    '''
    Position and velocity of a terrestrial observing station.

    Parameters
    ----------
    elong : float
        longitude (radians, east +ve)    
    phi : float
        latitude (geodetic, radians)    
    hm : float
        height above ref. ellipsoid (geodetic, m)    
    xp : float
        coordinates of the pole (radians)    
    yp : float
        coordinates of the pole (radians)     
    sp : float
        the TIO locator s' (radians)    
    theta : float
        Earth rotation angle (radians)

    Returns
    -------
    pv : list(2,3)
        position/velocity vector (m, m/s, CIRS)

    '''
    PV=[[0 for i in range(3)] for j in range(2)]
    
    #一天的秒长
    D2S=86400.0
    
    #2Pi
    D2PI=6.283185307179586476925287
    
    #地球在没UT1的1秒内转动的速率，弧度制
    OM=1.00273781191135448*D2PI/D2S
    
    #大地测量学坐标到地心坐标的转换(WGS84).
    XYZM,J=pymGd2gc(1,ELONG,PHI,HM)
    
    #极移以及TIO定位角的位置
    RPM=pymPom00(XP,YP,SP)
    XYZ=pymTrxp(RPM,XYZM)
    X=XYZ[0]
    Y=XYZ[1]
    Z=XYZ[2]

    #地球自转角的旋转.
    S=ma.sin(THETA)
    C=ma.cos(THETA)
    
    #位置向量
    PV[0][0]=C*X-S*Y
    PV[0][1]=S*X+C*Y
    PV[0][2]=Z
    
    #速度向量
    PV[1][0]=OM*(-S*X-C*Y)
    PV[1][1]=OM*(C*X-S*Y)
    PV[1][2]=0.0
    
    return(PV)


def pymApio(SP,THETA,ELONG,PHI,HM,XP,YP,REFA,REFB,ASTROM):
    '''
    For a terrestrial observer, prepare star-independent astrometry
    parameters for transformations between CIRS and observed
    coordinates.  The caller supplies the Earth orientation information
    and the refraction constants as well as the site coordinates.

    Parameters
    ----------
    sp : float
        the TIO locator s' (radians)    
    theta : float
        Earth rotation angle (radians)    
    elong : float
        longitude (radians, east +ve)    
    phi : float
        geodetic latitude (radians)    
    hm : float
        height above ellipsoid (m, geodetic)    
    xp : float
        polar motion coordinates (radians)    
    yp : float
        polar motion coordinates (radians)    
    refa : float
        refraction constant A (radians)    
    refb : float
        refraction constant B (radians)

    Returns
    -------
    astrom : list(30)
        star-independent astrometry parameters     
        [0]>pmt : unchanged
        [1-3]>eb : unchanged
        [4-6]>eh : unchanged
        [7]>em : unchanged
        [8-10]>v : unchanged
        [11]>bm1 : unchanged
        [12-20]>bpn : unchanged
        [21]>along : adjusted longitude (radians)
        [22]>xpl : polar motion xp wrt local meridian (radians)
        [23]>ypl : polar motion yp wrt local meridian (radians)
        [24]>sphi : sine of geodetic latitude
        [25]>cphi : cosine of geodetic latitude
        [26]>diurab : magnitude of diurnal aberration vector
        [27]>eral : "local" Earth rotation angle (radians)
        [28]>refa : refraction constant A (radians)
        [29]>refb : refraction constant B (radians)

    '''
    #光速(m/s)
    CMPS=299792458.0
    
    #构建从天球中间参考架（CIRS）到地面（地平坐标）的旋转矩阵
    R=pymIr()
    R=pymRz(THETA+SP,R)
    R=pymRy(-XP,R)
    R=pymRx(-YP,R)
    R=pymRz(ELONG,R)
    
    #获得‘当地’地球自转角度
    A=R[0][0]
    B=R[1][0]
    if (A!=0.0)|(B!=0.0):
        ERAL=ma.atan2(B,A)
    else:
        ERAL=0.0
    ASTROM[27]=ERAL
    
    #获得相对于当地子午线的极运动坐标[X,Y]
    A=R[0][0]
    C=R[2][0]
    ASTROM[22]=ma.atan2(C,ma.sqrt(A*A+B*B))
    A=R[2][1]
    B=R[2][2]
    if (A!=0.0)|(B!=0.0):
        ASTROM[23]=-ma.atan2(A,B)
    else:
        ASTROM[23]=0.0
    
    #调整经度
    ASTROM[21]=pymAnpm(ERAL-THETA)
    
    #纬度的函数
    ASTROM[24]=ma.sin(PHI)
    ASTROM[25]=ma.cos(PHI)
    
    #观察者的地心位置和速度(m, m/s, CIRS).
    PV=pymPvtob(ELONG,PHI,HM,XP,YP,SP,THETA)
    
    #周日光行差矢量的模
    ASTROM[26]=ma.sqrt(PV[1][0]*PV[1][0]+PV[1][1]*PV[1][1])/CMPS
    
    #折射常数
    ASTROM[28]=REFA
    ASTROM[29]=REFB
    
    return(ASTROM)


def pymApio13(UTC1,UTC2,DUT1,ELONG,PHI,HM,XP,YP,PHPA,TC,RH,WL,ASTROM):
    '''
    For a terrestrial observer, prepare star-independent astrometry
    parameters for transformations between CIRS and observed
    coordinates.  The caller supplies UTC, site coordinates, ambient air
    conditions and observing wavelength.

    Parameters
    ----------
    utc1 : float
        UTC as a 2-part Julian Date    
    utc2 : float
        UTC as a 2-part Julian Date    
    dut1 : float
        UT1-UTC (seconds)    
    elong : float
        longitude (radians, east +ve)    
    phi : float
        geodetic latitude (radians)    
    hm : float
        height above ellipsoid (m, geodetic)    
    xp : float
        polar motion coordinates (radians)    
    yp : float
        polar motion coordinates (radians)    
    phpa : float
        pressure at the observer (hPa = mB)    
    tc : float
        ambient temperature at the observer (deg C)    
    rh : float
        relative humidity at the observer (range 0-1)    
    wl : float
        wavelength (micrometers)    

    Raises
    ------
    ValueError
        1: 'dubious year',
       -1: 'unacceptable date',

    Returns
    -------
    astrom : list(30)
        star-independent astrometry parameters     
        [0]>pmt : unchanged
        [1-3]>eb : unchanged
        [4-6]>eh : unchanged
        [7]>em : unchanged
        [8-10]>v : unchanged
        [11]>bm1 : unchanged
        [12-20]>bpn : unchanged
        [21]>along : adjusted longitude (radians)
        [22]>xpl : polar motion xp wrt local meridian (radians)
        [23]>ypl : polar motion yp wrt local meridian (radians)
        [24]>sphi : sine of geodetic latitude
        [25]>cphi : cosine of geodetic latitude
        [26]>diurab : magnitude of diurnal aberration vector
        [27]>eral : "local" Earth rotation angle (radians)
        [28]>refa : refraction constant A (radians)
        [29]>refb : refraction constant B (radians)

    '''
    #从UTC转换到其他时制
    TAI1,TAI2,JS=pymUtctai(UTC1,UTC2)
    TT1,TT2,JS=pymTaitt(TAI1,TAI2)
    UT11,UT12,JS=pymUtcut1(UTC1,UTC2,DUT1)
    
    i=1
    while i<2:
        if (JS<0):
            print('ERROR',JS)
            break
        
        #获得TIO定位角
        SP=pymSp00(TT1,TT2)
        
        #获得地球自转角
        THETA=pymEra00(UT11,UT12)
        
        #获得折射常数
        REFA,REFB=pymRefco(PHPA,TC,RH,WL)
        
        #调用pymApio计算天体测量参数
        ASTROM=pymApio(SP,THETA,ELONG,PHI,HM,XP,YP,REFA,REFB,ASTROM)
        
        i+=1
    J=JS

    return(ASTROM,J)


def pymAper(THETA,ASTROM):
    '''
    In the star-independent astrometry parameters, update only the
    Earth rotation angle, supplied by the caller explicitly.

    Parameters
    ----------
    theta : float
        Earth rotation angle (radians)    
    astrom : list(30)
        star-independent astrometry parameters     
        [0]>pmt : not used
        [1-3]>eb : not used
        [4-6]>eh : not used
        [7]>em : not used
        [8-10]>v : not used
        [11]>bm1 : not used
        [12-20]>bpn : not used
        [21]>along : longitude + s' (radians)
        [22]>xpl : not used
        [23]>ypl : not used
        [24]>sphi : not used
        [25]>cphi : not used
        [26]>diurab : not used
        [27]>eral : not used
        [28]>refa : not used
        [29]>refb : not used

    Returns
    -------
    astrom : list(30)
        star-independent astrometry parameters     
        [0]>pmt : unchanged
        [1-3]>eb : unchanged
        [4-6]>eh : unchanged
        [7]>em : unchanged
        [8-10]>v : unchanged
        [11]>bm1 : unchanged
        [12-20]>bpn : unchanged
        [21]>along : unchanged
        [22]>xpl : unchanged
        [23]>ypl : unchanged
        [24]>sphi : unchanged
        [25]>cphi : unchanged
        [26]>diurab : unchanged
        [27]>eral : "local" Earth rotation angle (radians)
        [28]>refa : unchanged 
        [29]>refb : unchanged

    '''
    ASTROM[27]=THETA+ASTROM[21]
    
    return(ASTROM)


def pymAper13(UT11,UT12,ASTROM):
    '''
    In the star-independent astrometry parameters, update only the
    Earth rotation angle.  The caller provides UT1, (n.b. not UTC).

    Parameters
    ----------
    ut11 : float
        UT1 as a 2-part Julian Date    
    ut12 : float
        UT1 as a 2-part Julian Date    
    astrom : list(30)
        star-independent astrometry parameters     
        [0]>pmt : not used
        [1-3]>eb : not used
        [4-6]>eh : not used
        [7]>em : not used
        [8-10]>v : not used
        [11]>bm1 : not used
        [12-20]>bpn : not used
        [21]>along : longitude + s' (radians)
        [22]>xpl : not used
        [23]>ypl : not used
        [24]>sphi : not used
        [25]>cphi : not used
        [26]>diurab : not used
        [27]>eral : not used
        [28]>refa : not used
        [29]>refb : not used

    Returns
    -------
    astrom : list(30)
        star-independent astrometry parameters     
        [0]>pmt : unchanged
        [1-3]>eb : unchanged
        [4-6]>eh : unchanged
        [7]>em : unchanged
        [8-10]>v : unchanged
        [11]>bm1 : unchanged
        [12-20]>bpn : unchanged
        [21]>along : unchanged
        [22]>xpl : unchanged
        [23]>ypl : unchanged
        [24]>sphi : unchanged
        [25]>cphi : unchanged
        [26]>diurab : unchanged
        [27]>eral : "local" Earth rotation angle (radians)
        [28]>refa : unchanged 
        [29]>refb : unchanged

    '''
    A=pymEra00(UT11,UT12)
    ASTROM=pymAper(A,ASTROM)
    
    return(ASTROM)


def pymApcs13(DATE1,DATE2,PV,ASTROM):
    '''
    For an observer whose geocentric position and velocity are known,
    prepare star-independent astrometry parameters for transformations
    between ICRS and GCRS.  The Earth ephemeris is from SOFA models.
    
    The parameters produced by this function are required in the space
    motion, parallax, light deflection and aberration parts of the
    astrometric transformation chain.

    Parameters
    ----------
    date1 : float
        TDB as a 2-part Julian Date    
    date2 : float
        TDB as a 2-part Julian Date    
    pv : list(2,3)
        observer's geocentric pos/vel (m, m/s)    

    Returns
    -------
    astrom : list(30)
        star-independent astrometry parameters     
        [0]>pmt : PM time interval (SSB, Julian years)
        [1-3]>eb : SSB to observer (vector, au)
        [4-6]>eh : Sun to observer (unit vector)
        [7]>em : distance from Sun to observer (au)
        [8-10]>v : barycentric observer velocity (vector, c)
        [11]>bm1 : sqrt(1-|v|^2): reciprocal of Lorenz factor
        [12-20]>bpn : bias-precession-nutation matrix
        [21]>along : unchanged
        [22]>xpl : unchanged
        [23]>ypl : unchanged
        [24]>sphi : unchanged
        [25]>cphi : unchanged
        [26]>diurab : unchanged
        [27]>eral : unchanged
        [28]>refa : unchanged 
        [29]>refb : unchanged

    '''

    #地球相对于太阳系质心和日心的位置和速度 (au, au/d).
    EHPV,EBPV,J=pymEpv00(DATE1,DATE2)
    EHP = EHPV[0]
    #计算不依赖于恒星的天体测量参数
    ASTROM=pymApcs(DATE1,DATE2,PV,EBPV,EHP,ASTROM)   
    
    return(ASTROM)


def pymS06(DATE1,DATE2,X,Y):
    '''
    The CIO locator s, positioning the Celestial Intermediate Origin on
    the equator of the Celestial Intermediate Pole, given the CIP's X,Y
    coordinates.  Compatible with IAU 2006/2000A precession-nutation.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date    
    x : float
        CIP coordinates    
    y : float
        CIP coordinates

    Returns
    -------
    function value : float
        the CIO locator s in radians

    '''
    #1角秒对应的弧度
    DAS2R=4.848136811095359935899141e-6
    
    #参考历元(J2000.0), JD
    DJ00=2451545.0
    
    #1儒略世纪对应的天数
    DJC=36525.0
    
    #自J2000.0起经过的儒略世纪数，此处为定义
    T=0
    
    #基本的参数
    FA=[0 for i in range(8)]
    
    # *  ---------------------
    # *  The series for s+XY/2
    # *  ---------------------
    
    NSP=6
    NS0=33
    NS1=3
    NS2=25
    NS3=4
    NS4=1
    
    #多项式系数
    SP=[94e-6,3808.65e-6,-122.68e-6,-72574.11e-6,27.98e-6,15.62e-6]
    
    #Coefficients of l,l',F,D,Om,LVe,LE,pA
    #Argument coefficients for t^0
    KS0=[[0,0,0,0,1,0,0,0],
            [0,0,0,0,2,0,0,0],
            [0,0,2,-2,3,0,0,0],
            [0,0,2,-2,1,0,0,0],
            [0,0,2,-2,2,0,0,0],
            [0,0,2,0,3,0,0,0],
            [0,0,2,0,1,0,0,0],
            [0,0,0,0,3,0,0,0],
            [0,1,0,0,1,0,0,0],
            [0,1,0,0,-1,0,0,0],
            [1,0,0,0,-1,0,0,0],
            [1,0,0,0,1,0,0,0],
            [0,1,2,-2,3,0,0,0],
            [0,1,2,-2,1,0,0,0],
            [0,0,4,-4,4,0,0,0],
            [0,0,1,-1,1,-8,12,0],
            [0,0,2,0,0,0,0,0],
            [0,0,2,0,2,0,0,0],
            [1,0,2,0,3,0,0,0],
            [1,0,2,0,1,0,0,0],
            [0,0,2,-2,0,0,0,0],
            [0,1,-2,2,-3,0,0,0],
            [0,1,-2,2,-1,0,0,0],
            [0,0,0,0,0,8,-13,-1],
            [0,0,0,2,0,0,0,0],
            [2,0,-2,0,-1,0,0,0],
            [0,1,2,-2,2,0,0,0],
            [1,0,0,-2,1,0,0,0],
            [1,0,0,-2,-1,0,0,0],
            [0,0,4,-2,4,0,0,0],
            [0,0,2,-2,4,0,0,0],
            [1,0,-2,0,-3,0,0,0],
            [1,0,-2,0,-1,0,0,0]]
    
    #Argument coefficients for t^1
    KS1=[[0,0,0,0,2,0,0,0],
            [0,0,0,0,1,0,0,0],
            [0,0,2,-2,3,0,0,0]]
    
    #Argument coefficients for t^2
    KS2=[[0,0,0,0,1,0,0,0],
            [0,0,2,-2,2,0,0,0],
            [0,0,2,0,2,0,0,0],
            [0,0,0,0,2,0,0,0],
            [0,1,0,0,0,0,0,0],
            [1,0,0,0,0,0,0,0],
            [0,1,2,-2,2,0,0,0],
            [0,0,2,0,1,0,0,0],
            [1,0,2,0,2,0,0,0],
            [0,1,-2,2,-2,0,0,0],
            [1,0,0,-2,0,0,0,0],
            [0,0,2,-2,1,0,0,0],
            [1,0,-2,0,-2,0,0,0],
            [0,0,0,2,0,0,0,0],
            [1,0,0,0,1,0,0,0],
            [1,0,-2,-2,-2,0,0,0],
            [1,0,0,0,-1,0,0,0],
            [1,0,2,0,1,0,0,0],
            [2,0,0,-2,0,0,0,0],
            [2,0,-2,0,-1,0,0,0],
            [0,0,2,2,2,0,0,0],
            [2,0,2,0,2,0,0,0],
            [2,0,0,0,0,0,0,0],
            [1,0,2,-2,2,0,0,0],
            [0,0,2,0,0,0,0,0]]
    
    #Argument coefficients for t^3
    KS3=[[0,0,0,0,1,0,0,0],
            [0,0,2,-2,2,0,0,0],
            [0,0,2,0,2,0,0,0],
            [0,0,0,0,2,0,0,0]]
    
    #Argument coefficients for t^4
    KS4=[[0,0,0,0,1,0,0,0]]
 
    #Sine and cosine coefficients
    #Sine and cosine coefficients for t^0
    SS0=[[-2640.73e-6,+0.39e-6],
            [-63.53e-6,+0.02e-6],
            [-11.75e-6,-0.01e-6],
            [-11.21e-6,-0.01e-6],
            [+4.57e-6,0.00e-6],
            [-2.02e-6,0.00e-6],
            [-1.98e-6,0.00e-6],
            [+1.72e-6,0.00e-6],
            [+1.41e-6,+0.01e-6],
            [+1.26e-6,+0.01e-6],
            [+0.63e-6,0.00e-6],
            [+0.63e-6,0.00e-6],
            [-0.46e-6,0.00e-6],
            [-0.45e-6,0.00e-6],
            [-0.36e-6,0.00e-6],
            [+0.24e-6,+0.12e-6],
            [-0.32e-6,0.00e-6],
            [-0.28e-6,0.00e-6],
            [-0.27e-6,0.00e-6],
            [-0.26e-6,0.00e-6],
            [+0.21e-6,0.00e-6],
            [-0.19e-6,0.00e-6],
            [-0.18e-6,0.00e-6],
            [+0.10e-6,-0.05e-6],
            [-0.15e-6,0.00e-6],
            [+0.14e-6,0.00e-6],
            [+0.14e-6,0.00e-6],
            [-0.14e-6,0.00e-6],
            [-0.14e-6,0.00e-6],
            [-0.13e-6,0.00e-6],
            [+0.11e-6,0.00e-6],
            [-0.11e-6,0.00e-6],
            [-0.11e-6,0.00e-6]]
    
    #Sine and cosine coefficients for t^1
    SS1=[[-0.07e-6,+3.57e-6],
            [+1.73e-6,-0.03e-6],
            [0.00e-6,+0.48e-6]]
    
    #Sine and cosine coefficients for t^2
    SS2=[[+743.52e-6,-0.17e-6],
            [+56.91e-6,+0.06e-6],
            [+9.84e-6,-0.01e-6],
            [-8.85e-6,+0.01e-6],
            [-6.38e-6,-0.05e-6],
            [-3.07e-6,0.00e-6],
            [+2.23e-6,0.00e-6],
            [+1.67e-6,0.00e-6],
            [+1.30e-6,0.00e-6],
            [+0.93e-6,0.00e-6],
            [+0.68e-6,0.00e-6],
            [-0.55e-6,0.00e-6],
            [+0.53e-6,0.00e-6],
            [-0.27e-6,0.00e-6],
            [-0.27e-6,0.00e-6],
            [-0.26e-6,0.00e-6],
            [-0.25e-6,0.00e-6],
            [+0.22e-6,0.00e-6],
            [-0.21e-6,0.00e-6],
            [+0.20e-6,0.00e-6],
            [+0.17e-6,0.00e-6],
            [+0.13e-6,0.00e-6],
            [-0.13e-6,0.00e-6],
            [-0.12e-6,0.00e-6],
            [-0.11e-6,0.00e-6]]
    
    #Sine and cosine coefficients for t^3
    SS3=[[+0.30e-6,-23.42e-6],
            [-0.03e-6,-1.46e-6],
            [-0.01e-6,-0.25e-6],
            [0.00e-6,+0.23e-6]]
    
    #Sine and cosine coefficients for t^4
    SS4=[[-0.26e-6,-0.01e-6]]
    
    #参考历元J2000.0和当前日期之间的间隔，儒略世纪数
    T=((DATE1-DJ00)+DATE2)/DJC

    #Fundamental Arguments (from IERS Conventions 2003)

    #月亮的平近点角
    FA[0]=pymFal03(T)
    
    #太阳的平近点角
    FA[1]=pymFalp03(T)

    #月亮的平黄经减去升交点黄经
    FA[2]=pymFaf03(T)

    #月亮到太阳的平均距角
    FA[3]=pymFad03(T)

    #月亮的升交点平黄经
    FA[4]=pymFaom03(T)

    #金星的平黄经
    FA[5]=pymFave03(T)

    #地球的平黄经
    FA[6]=pymFae03(T)

    #黄经上的累计进动
    FA[7]=pymFapa03(T)
    
    #估值 s.
    S0 = SP[0]
    S1 = SP[1]
    S2 = SP[2]
    S3 = SP[3]
    S4 = SP[4]
    S5 = SP[5]
    
    for I in range(NS0-1,-1,-1):
        A=0.0
        for J in range(8):
            A=A+float(KS0[I][J])*FA[J]
        S0=S0+(SS0[I][0]*ma.sin(A)+SS0[I][1]*ma.cos(A))
    
    for I in range(NS1-1,-1,-1):
        A=0.0    
        for J in range(8):
            A=A+float(KS1[I][J])*FA[J]
        S1=S1+(SS1[I][0]*ma.sin(A)+SS1[I][1]*ma.cos(A))
    
    for I in range(NS2-1,-1,-1):
        A=0.0
        for J in range(8):
            A=A+float(KS2[I][J])*FA[J]
        S2=S2+(SS2[I][0]*ma.sin(A)+SS2[I][1]*ma.cos(A))
    
    for I in range(NS3-1,-1,-1):
        A=0.0    
        for J in range(8):
            A=A+float(KS3[I][J])*FA[J]
        S3=S3+(SS3[I][0]*ma.sin(A)+SS3[I][1]*ma.cos(A))
    
    for I in range(NS4-1,-1,-1):
        A=0.0    
        for J in range(8):
            A=A+float(KS4[I][J])*FA[J]
        S4=S4+(SS4[I][0]*ma.sin(A)+SS4[I][1]*ma.cos(A))
        
    S06=(S0+(S1+(S2+(S3+(S4+S5*T)*T)*T)*T)*T)*DAS2R-X*Y/2.0
        
    return(S06)


def pymEors(RNPB,S):
    '''
    Equation of the origins, given the classical NPB matrix and the
    quantity s.

    Parameters
    ----------
    rnpb : list(3,3) 
        classical nutation x precession x bias matrix    
    s : float
        the quantity s (the CIO locator) in radians

    Returns
    -------
    function value : float
        the equation of the origins in radians

    '''
    #由 Wallace & Capitaine (2006) expression (16)估计.
    X=RNPB[2][0]
    AX=X/(1.0+RNPB[2][2])
    XS=1.0-AX*X
    YS=-AX*RNPB[2][1]
    ZS=-X
    P=RNPB[0][0]*XS+RNPB[0][1]*YS+RNPB[0][2]*ZS
    Q=RNPB[1][0]*XS+RNPB[1][1]*YS+RNPB[1][2]*ZS
    
    if (P!=0.0)|(Q!=0.0):
        EORS=S-ma.atan2(Q,P)
    else:
        EORS=S
    
    return(EORS)


def pymBpn2xy(RBPN):
    '''
    Extract from the bias-precession-nutation matrix the X,Y coordinates
    of the Celestial Intermediate Pole.

    Parameters
    ----------
    rbpn : list(3,3)
        celestial-to-true matrix

    Returns
    -------
    x : float
        Celestial Intermediate Pole    
    y : float
        Celestial Intermediate Pole

    '''
    #分离出X,Y坐标
    X=RBPN[2][0]
    Y=RBPN[2][1]

    return(X,Y)


def pymApci13(DATE1,DATE2,ASTROM):
    '''
    For a terrestrial observer, prepare star-independent astrometry
    parameters for transformations between ICRS and geocentric CIRS
    coordinates.  The caller supplies the date, and SOFA models are used
    to predict the Earth ephemeris and CIP/CIO.
    
    The parameters produced by this function are required in the
    parallax, light deflection, aberration, and bias-precession-nutation
    parts of the astrometric transformation chain.

    Parameters
    ----------
    date1 : float
        TDB as a 2-part Julian Date    
    date2 : float
        TDB as a 2-part Julian Date     

    Returns
    -------
    astrom : list(30)
        star-independent astrometry parameters     
    eo : float
        equation of the origins (ERA-GST)
        astrom:
        [0]>pmt : PM time interval (SSB, Julian years)
        [1-3]>eb : SSB to observer (vector, au)
        [4-6]>eh : Sun to observer (unit vector)
        [7]>em : distance from Sun to observer (au)
        [8-10]>v : barycentric observer velocity (vector, c)
        [11]>bm1 : sqrt(1-|v|^2): reciprocal of Lorenz factor
        [12-20]>bpn : bias-precession-nutation matrix
        [21]>along : unchanged
        [22]>xpl : unchanged
        [23]>ypl : unchanged
        [24]>sphi : unchanged
        [25]>cphi : unchanged
        [26]>diurab : unchanged
        [27]>eral : unchanged
        [28]>refa : unchanged 
        [29]>refb : unchanged
    
    '''
    #调用pymEpv00获得地球相对于质心和日心的位置速度信息(au, au/d).
    PVH,PVB,J=pymEpv00(DATE1,DATE2)
    PV = PVH[0]
    #调用pymPnm06a获得极移进动章动矩阵 IAU 2006/2000A.
    R=pymPnm06a(DATE1,DATE2)
    
    #从极移进洞章动矩阵中获得 CIP X,Y.
    X,Y=pymBpn2xy(R)
    
    #调用pymS06获得CIO定位角.
    S=pymS06(DATE1,DATE2,X,Y)
    
    #计算与恒星无关的天体测量参数.
    ASTROM=pymApci(DATE1,DATE2,PVB,PV,X,Y,S,ASTROM)
    
    #调用pymEors获得真春分点与天球中间原点之间的距离
    EO=pymEors(R,S)
    
    return(ASTROM,EO)


def pymRxpv(R,PV):
    '''
    Multiply a pv-vector by an r-matrix.

    Parameters
    ----------
    r : list(3,3)
        r-matrix    
    pv : list(2,3)
        pv-vector

    Returns
    -------
    rpv : list(2,3)
        r * pv

    '''
    RPV=[0,0]
    
    RPV[0]=pymRxp(R,PV[0])
    RPV[1]=pymRxp(R,PV[1])    
    
    return(RPV)


def pymTrxpv(R,PV):
    '''
    Multiply a pv-vector by the transpose of an r-matrix.

    Parameters
    ----------
    r : list(3,3)
        r-matrix    
    pv : list(2,3)
        pv-vector

    Returns
    -------
    trpv : list(2,3)
        r^T * pv

    '''
    RI=pymTr(R)
    
    TRPV=pymRxpv(RI,PV)    
    
    return(TRPV)


def pymApco(DATE1,DATE2,EBPV,EHP,X,Y,S,THETA,
             ELONG,PHI,HM,XP,YP,SP,REFA,REFB):
    '''
    For a terrestrial observer, prepare star-independent astrometry
    parameters for transformations between ICRS and observed
    coordinates.  The caller supplies the Earth ephemeris, the Earth
    rotation information and the refraction constants as well as the
    site coordinates.

    Parameters
    ----------
    date1 : float
        TDB as a 2-part Julian Date    
    date2 : float
        TDB as a 2-part Julian Date     
    ebpv : list(2,3)
        Earth barycentric position/velocity (au, au/day)     
    ehp : list(3)
        Earth heliocentric position (au)    
    x : flaot
        CIP X,Y (components of unit vector)    
    y : float
        CIP X,Y (components of unit vector)    
    s : float
        the CIO locator s (radians)    
    theta : float
        Earth rotation angle (radians)
    elong : float
        longitude (radians, east +ve)
    phi : float
        latitude (geodetic, radians)
    hm : float
        height above ellipsoid (m, geodetic)
    xp : float
        polar motion coordinates (radians)
    yp : float
        polar motion coordinates (radians)
    sp : float
        the TIO locator s' (radians)
    refa : float
        refraction constant A (radians)
    refb : float
        refraction constant B (radians)

    Returns
    -------
    astrom : list(30)
        star-independent astrometry parameters    
        [0]>pmt : PM time interval (SSB, Julian years)
        [1-3]>eb : SSB to observer (vector, au)
        [4-6]>eh : Sun to observer (unit vector)
        [7]>em : distance from Sun to observer (au)
        [8-10]>v : barycentric observer velocity (vector, c)
        [11]>bm1 : sqrt(1-|v|^2): reciprocal of Lorenz factor
        [12-20]>bpn : bias-precession-nutation matrix
        [21]>along : adjusted longitude (radians)
        [22]>xpl : polar motion xp wrt local meridian (radians)
        [23]>ypl : polar motion yp wrt local meridian (radians)
        [24]>sphi : sine of geodetic latitude
        [25]>cphi : cosine of geodetic latitude
        [26]>diurab : magnitude of diurnal aberration vector
        [27]>eral : "local" Earth rotation angle (radians)
        [28]>refa : refraction constant A (radians)
        [29]>refb : refraction constant B (radians)

    '''
    #构建从天球中间参考系到时角坐标的旋转矩阵.
    R=pymIr()
    R=pymRz(THETA+SP,R)
    R=pymRy(-XP,R)
    R=pymRx(-YP,R)
    R=pymRz(ELONG,R)
    ASTROM=[0 for i in range(30)]
    
    #求‘当地’地球自转角度
    A=R[0][0]
    B=R[1][0]
    if (A!=0.0)|(B!=0.0):
        ERAL=ma.atan2(B,A)
    else:
        ERAL=0.0
    ASTROM[27]=ERAL
    
    #求关于局部子午线的极运动[X,Y]
    A=R[0][0]
    C=R[2][0]
    ASTROM[22]=ma.atan2(C,ma.sqrt(A*A+B*B))
    A=R[2][1]
    B=R[2][2]
    if (A!=0.0)|(B!=0.0):
        ASTROM[23]=-ma.atan2(A,B)
    else:
        ASTROM[23]=0.0
    
    #调整的经度，引入地球自转角
    D=ERAL-THETA
    ASTROM[21]=pymAnpm(D)
    
    #纬度相关函数   
    ASTROM[24]=ma.sin(PHI)
    ASTROM[25]=ma.cos(PHI)
    
    #折射率参数
    ASTROM[28]=REFA
    ASTROM[29]=REFB
    
    #忽略周日光行差
    ASTROM[26]=0.0
    
    #基于CIO的极移进动章动矩阵
    R=pymC2ixys(X,Y,S)
    
    #调用pymPvtob获得地面观测站的位置速度信息 (m, m/s, CIRS).
    PVC=pymPvtob(ELONG,PHI,HM,XP,YP,SP,THETA)
    
    #将位置速度向量旋转到地心天球参考系GCRS.
    PV=pymTrxpv(R,PVC)
    
    #调用pymApcs获得相关的天体测量参数（ICRS <-> GCRS）
    ASTROM=pymApcs(DATE1,DATE2,PV,EBPV,EHP,ASTROM)
    
    #保存基于极移进动章动矩阵（BPN）获得的CIO（天球中间零点）
    k=0
    for i in range(3):
        for j in range(3):
            ASTROM[12+k]=R[i][j]
            k+=1
            
    return(ASTROM)


def pymApco13(UTC1,UTC2,DUT1,ELONG,PHI,HM,XP,YP,PHPA,TC,RH,WL):
    '''
    For a terrestrial observer, prepare star-independent astrometry
    parameters for transformations between ICRS and observed
    coordinates.  The caller supplies UTC, site coordinates, ambient air
    conditions and observing wavelength, and SOFA models are used to
    obtain the Earth ephemeris, CIP/CIO and refraction constants.
   
    The parameters produced by this function are required in the
    parallax, light deflection, aberration, and bias-precession-nutation
    parts of the ICRS/CIRS transformations.

    Parameters
    ----------
    utc1 : float
        UTC as a 2-part quasi Julian Date    
    utc2 : float
        UTC as a 2-part quasi Julian Date    
    dut1 : float
        UT1-UTC (seconds)    
    elong : float
        longitude (radians, east +ve)    
    phi : float
        latitude (geodetic, radians)    
    hm : float
        height above ellipsoid (m, geodetic)    
    xp : float
        polar motion coordinates (radians)    
    yp : float 
        polar motion coordinates (radians)    
    phpa : float
        pressure at the observer (hPa = mB)    
    tc : float
        ambient temperature at the observer (deg C)    
    rh : float
        relative humidity at the observer (range 0-1)    
    wl : float
        wavelength (micrometers)    
        
    Raises
    ------
    ValueError
        1: 'dubious year',
       -1: 'unacceptable date',
    
    Returns
    -------
    astrom : list(30)
        star-independent astrometry parameters     
    eo : float
        equation of the origins (ERA-GST)
    J : ValueError
        1: 'dubious year',
       -1: 'unacceptable date',
        astrom:
        [0]>pmt : PM time interval (SSB, Julian years)
        [1-3]>eb : SSB to observer (vector, au)
        [4-6]>eh : Sun to observer (unit vector)
        [7]>em : distance from Sun to observer (au)
        [8-10]>v : barycentric observer velocity (vector, c)
        [11]>bm1 : sqrt(1-|v|^2): reciprocal of Lorenz factor
        [12-20]>bpn : bias-precession-nutation matrix
        [21]>along : adjusted longitude (radians)
        [22]>xpl : polar motion xp wrt local meridian (radians)
        [23]>ypl : polar motion yp wrt local meridian (radians)
        [24]>sphi : sine of geodetic latitude
        [25]>cphi : cosine of geodetic latitude
        [26]>diurab : magnitude of diurnal aberration vector
        [27]>eral : "local" Earth rotation angle (radians)
        [28]>refa : refraction constant A (radians)
        [29]>refb : refraction constant B (radians)

    '''
    
    i=1
    while i<2:

        #从协调世界时UTC转换到其他时间
        TAI1,TAI2,JS=pymUtctai(UTC1,UTC2)
        if (JS<0):
            print('ERROR1',JS)
            break
        TT1,TT2,JS=pymTaitt(TAI1,TAI2)
        UT11,UT12,JS=pymUtcut1(UTC1,UTC2,DUT1)
        if (JS<0):
            print('ERROR2',JS)
            break
        
        #调用pymEpv00获得地球相对与质心与日心的位置速度向量.
        EHPV,EBPV,JW=pymEpv00(TT1,TT2)
        EHP = EHPV[0]
        #构建基于春分的极移进动章动（BPN）矩阵, IAU 2006/2000A.
        R=pymPnm06a(TT1,TT2)
    
        #从极移进动章动矩阵中分理处 CIP X,Y.
        X,Y=pymBpn2xy(R)
        
        #调用pymS06获得CIO定位角s.
        S=pymS06(TT1,TT2,X,Y)
        
        #地球的自转角
        THETA=pymEra00(UT11,UT12)
        
        #调用pymSp00获得TIO定位角s'.
        SP=pymSp00(TT1,TT2)
        
        #调用pymRefco，通过环境参数获得折射参数
        REFA,REFB=pymRefco(PHPA,TC,RH,WL)
        
        #通过pymApco计算与恒星无关的天体测量参数
        ASTROM=pymApco(TT1,TT2,EBPV,EHP,X,Y,S,THETA,ELONG,
                        PHI,HM,XP,YP,SP,REFA,REFB)
        
        #获得真春分点与天球中间原点之间的距离
        EO=pymEors(R,S)
        
        i+=1
    
    J=JS
        
    return(ASTROM, EO, J)


def pymPmpx(RC,DC,PR,PD,PX,RV,PMT,POB):
    '''
    Proper motion and parallax.

    Parameters
    ----------
    rc : float
        ICRS RA,Dec at catalog epoch (radians)    
    dc : float
        ICRS RA,Dec at catalog epoch (radians)    
    pr : float
        RA proper motion (radians/year)    
    pd : float
        Dec proper motion (radians/year)    
    px : float
        parallax (arcsec)    
    rv : float
        radial velocity (km/s, +ve if receding)    
    pmt : float
        proper motion time interval (SSB, Julian years)    
    pob : list(3)
        SSB to observer vector (au)    

    Returns
    -------
    pco : list(3)
        coordinate direction (BCRS unit vector)

    '''
    #一天的秒长
    D2S=86400.0
    
    #1儒略年的天数
    DJY=365.25
    
    #1千儒略年的天数
    DJM=365250.0
    
    #光速(m/s)
    CMPS=299792458.0
    
    #天文单位(m, IAU 2012)
    AUM=149597870.7e3
    
    #光走过1AU经过的时间，单位：儒略年
    AULTY=AUM/CMPS/D2S/DJY
    
    #Km/s到au/year的单位换算
    VF=D2S*DJM/AUM
    
    #1角秒对应的弧度
    DAS2R=4.848136811095359935899141e-6
    
    P=[0,0,0]
    PM=[0,0,0]
    
    #将在球面上的星表坐标转换为单位向量
    SR=ma.sin(RC)
    CR=ma.cos(RC)
    SD=ma.sin(DC)
    CD=ma.cos(DC)
    X=CR*CD
    Y=SR*CD
    Z=SD
    P[0]=X
    P[1]=Y
    P[2]=Z    
    
    #观测矢量在恒星方向上的分量
    PDB=pymPdp(P,POB)
    
    #自行的时间间隔(y)，包括罗默效应（关于光速的累计效应？）
    DT=PMT+PDB*AULTY
    
    #空间运动(radians per year).
    PXR=PX*DAS2R
    W=VF*RV*PXR
    PDZ=PD*Z
    PM[0]=-PR*Y-PDZ*CR+W*X
    PM[1]=PR*X-PDZ*SR+W*Y
    PM[2]=PD*CD+W*Z
    
    #相对于太阳的坐标方向(unit vector, BCRS).
    for i in range(3):
        P[i]=P[i]+DT*PM[i]-PXR*POB[i]
    W,PCO=pymPn(P)
    
    return(PCO)


def pymC2s(P):
    '''
    P-vector to spherical coordinates.

    Parameters
    ----------
    p : list(3)
        p-vector

    Returns
    -------
    theta : float
        longitude angle (radians)    
    phi : float
        latitude angle (radians)

    '''
    X=P[0]
    Y=P[1]
    Z=P[2]
    D2=X*X+Y*Y
    
    if (D2==0.0):
        THETA=0.0
    else:
        THETA=ma.atan2(Y,X)
    
    if (Z==0.0):
        PHI=0.0
    else:
        PHI=ma.atan2(Z,ma.sqrt(D2))
        
    return(THETA,PHI)


def pymAtccq(RC,DC,PR,PD,PX,RV,ASTROM):
    '''
    Quick transformation of a star's ICRS catalog entry (epoch J2000.0)
    into ICRS astrometric place, given precomputed star-independent
    astrometry parameters.
    
    Use of this function is appropriate when efficiency is important and
    where many star positions are to be transformed for one date.  The
    star-independent parameters can be obtained by calling one of the
    functions iauApci[13], iauApcg[13], iauApco[13] or iauApcs[13].
    
    If the parallax and proper motions are zero the transformation has
    no effect.

    Parameters
    ----------
    rc : float
        ICRS right ascension at J2000.0 (radians)    
    dc : float
        ICRS declination at J2000.0 (radians)    
    pr : float
        RA proper motion (radians/year)       
    pd : float
        Dec proper motion (radians/year)    
    px : float
        parallax (arcsec)    
    rv : float
        radial velocity (km/s, +ve if receding)     
    astrom : list(30)
        star-independent astrometry parameters     
        [0]>pmt : PM time interval (SSB, Julian years)
        [1-3]>eb : SSB to observer (vector, au)
        [4-6]>eh : Sun to observer (unit vector)
        [7]>em : distance from Sun to observer (au)
        [8-10]>v : barycentric observer velocity (vector, c)
        [11]>bm1 : sqrt(1-|v|^2): reciprocal of Lorenz factor
        [12-20]>bpn : bias-precession-nutation matrix
        [21]>along : adjusted longitude (radians)
        [22]>xpl : polar motion xp wrt local meridian (radians)
        [23]>ypl : polar motion yp wrt local meridian (radians)
        [24]>sphi : sine of geodetic latitude
        [25]>cphi : cosine of geodetic latitude
        [26]>diurab : magnitude of diurnal aberration vector
        [27]>eral : "local" Earth rotation angle (radians)
        [28]>refa : refraction constant A (radians)
        [29]>refb : refraction constant B (radians)

    Returns
    -------
    ra : float
        ICRS astrometric RA (radians)    
    da : float
        ICRS astrometric Dec (radians)

    '''
    #调用pymPmpx获得恒星BCRS的坐标方向
    B=[ASTROM[1],ASTROM[2],ASTROM[3]]
    P=pymPmpx(RC,DC,PR,PD,PX,RV,ASTROM[0],B)
    
    #转换为ICRS赤经赤纬坐标
    W,DA=pymC2s(P)
    RA=pymAnp(W)
    
    return(RA,DA)


def pymAtcc13(RC,DC,PR,PD,PX,RV,DATE1,DATE2):
    '''
    Transform a star's ICRS catalog entry (epoch J2000.0) into ICRS
    astrometric place.

    Parameters
    ----------
    rc : float
        ICRS right ascension at J2000.0 (radians)    
    dc : float
        ICRS declination at J2000.0 (radians)    
    pr : float
        RA proper motion (radians/year)       
    pd : float
        Dec proper motion (radians/year)    
    px : float
        parallax (arcsec)    
    rv : float
        radial velocity (km/s, +ve if receding)     
    date1 : float
        TDB as a 2-part Julian Date    
    date2 : float
        TDB as a 2-part Julian Date    

    Returns
    -------
    ra : float
        ICRS astrometric RA (radians)    
    da : float
        ICRS astrometric Dec (radians)

    '''
    ASTROM=[0 for i in range(30)]
    
    #调用pymApci13获得与恒星无关的天体测量参数ASTROM
    ASTROM,W=pymApci13(DATE1,DATE2,ASTROM)
    
    RA,DA=pymAtccq(RC,DC,PR,PD,PX,RV,ASTROM)
    
    return(RA,DA)


def pymPxp (A,B):
    ''''
    p-vector outer (=vector=cross) product.

    Parameters
    ----------
    a : list(3) 
        first p-vector    
    b : list(3) 
        second p-vector

    Returns
    -------
    axb : list(3) 
        a x b

    '''

    AXB=[0,0,0]
    
    XA=A[0]
    YA=A[1]
    ZA=A[2]
    XB=B[0]
    YB=B[1]
    ZB=B[2]
    AXB[0]=YA*ZB-ZA*YB
    AXB[1]=ZA*XB-XA*ZB
    AXB[2]=XA*YB-YA*XB
       
    return(AXB)


def pymLd(BM,P,Q,E,EM,DLIM):
    '''
    Apply light deflection by a solar-system body, as part of
    transforming coordinate direction into natural direction.

    Parameters
    ----------
    bm : float
        mass of the gravitating body (solar masses)    
    p : list(3)
        direction from observer to source (unit vector)    
    q : list(3)
        direction from body to source (unit vector)    
    e : list(3)
        direction from body to observer (unit vector)    
    em : flaot
        distance from body to observer (au)    
    dlim : flaot
        deflection limiter

    Returns
    -------
    p1 : list(3)
        observer to deflected source (unit vector)

    '''
    #  太阳的施瓦西半径 (au)
    #  = 2 * 1.32712440041 D20 / (2.99792458 D8)^2 / 1.49597870700 D11
    SRS=1.97412574336e-08
    
    QPE=[0,0,0]
    P1=[0,0,0]
    
    #Q . (Q + E).
    for i in range(3):
        QPE[i]=Q[i]+E[i]
    QDQPE=pymPdp(Q, QPE)
    
    #2 x G x BM / ( EM x c^2 x ( Q . (Q + E) ) ).
    W=BM*SRS/EM/max(QDQPE,DLIM)
    
    #P x (E x Q).
    EQ=pymPxp(E,Q)
    PEQ=pymPxp(P,EQ)
    
    #应用偏转
    for i in range(3):
        P1[i]=P[i]+W*PEQ[i]
    
    return(P1)


def pymPmp(A,B):
    '''
    P-vector subtraction.

    Parameters
    ----------
    a : list(3)
        first p-vector    
    b : list(3)
        second p-vector

    Returns
    -------
    amb : list(3)
        a - b

    '''
    AMB=[0,0,0]
    for i in range(3):
        AMB[i]=A[i]-B[i]
    
    return(AMB)


def pymPpsp (A,S,B):
    '''
    P-vector plus scaled p-vector.

    Parameters
    ----------
    a : list(3)
        first p-vector    
    s : float     
        scalar (multiplier for b)
    b : list(3)
        second p-vector

    Returns
    -------
    apsb : list(3)
        a + s*b

    '''
    APSB=[0,0,0]
    
    for i in range(3):
        APSB[i]=A[i]+S*B[i]
    
    return(APSB)


def pymLdn (N,B,OB,SC):
    '''
    For a star, apply light deflection by multiple solar-system bodies,
    as part of transforming coordinate direction into natural direction.

    Parameters
    ----------
    n : int
        number of bodies    
    b: list(n,8)
        data for each of the n bodies    
        [i][0]>bm : mass of the body (solar masses)    
        [i][1]>dl : deflection limiter
        [i][2-7]>pv : barycentric PV of the body (au, au/day)
    ob : list(3)
        barycentric position of the observer (au)    
    sc : list(3)
        observer to star coord direction (unit vector)

    Returns
    -------
    sn : list(3)
        observer to deflected star (unit vector)

    '''          
    #一天的秒长
    DAYSEC=86400.0
    
    #光速(m/s)
    CMPS=299792458.0
    
    #天文单位 (m, IAU 2012)
    AUM=149597870.7e3
    
    #光经过1AU所需要的时间(day)
    CR=AUM/CMPS/DAYSEC
    
    #信号源在偏转之前的方向
    S=pymCp(SC)
    
    #对每个天体依序处理
    for i in range(N):
        
        #天体到观察者的位置向量，单位：AU
        C=[B[i][2],B[i][3],B[i][4]]
        V=pymPmp(OB,C)
        
        #减去光线经过天体的时间(days)
        D=pymPdp(S,V)
        DT=D*CR
        
        #修正处理信号源在观察者背后的情况
        DT=min(DT,0.0)
        
        #从天体回溯到光线经过天体的时间。
        H=[B[i][5],B[i][6],B[i][7]]
        EV=pymPpsp(V,-DT,H)
        
        #将天体到观察者矢量分离成大小和方向。
        EM,E=pymPn(EV)
        
        #对这个天体应用光线偏转
        V=pymLd(B[i][0],S,S,E,EM,B[i][1])
        
        #更新信号源的方向
        S=pymCp(V)
    
    #给出信号源最终的方向
    SN=pymCp(S)
    
    return(SN)


def pymLdsun (P,E,EM):
    '''
    Deflection of starlight by the Sun.

    Parameters
    ----------
    p : list(3)
        direction from observer to star (unit vector)    
    e : list(3)
        direction from Sun to observer (unit vector)    
    em : float
        distance from Sun to observer (au)

    Returns
    -------
    p1 : list(3)
        observer to deflected star (unit vector)

    '''
    #偏移限制(对遥远的观察者来说更小)
    DLIM=1e-6/max(EM*EM,1.0)
    #应用偏转.
    P1=pymLd(1.0,P,P,E,EM,DLIM)

    return(P1)


def pymAtciq(RC,DC,PR,PD,PX,RV,ASTROM):
    '''
    Quick ICRS, epoch J2000.0, to CIRS transformation, given precomputed
    star-independent astrometry parameters.
    
    Use of this function is appropriate when efficiency is important and
    where many star positions are to be transformed for one date.  The
    star-independent parameters can be obtained by calling one of the
    functions iauApci[13], iauApcg[13], iauApco[13] or iauApcs[13].
    
    If the parallax and proper motions are zero the iauAtciqz function
    can be used instead.
    
    Parameters
    ----------
    rc : float
        ICRS right ascension at J2000.0 (radians)    
    dc : float
        ICRS declination at J2000.0 (radians)    
    pr : float
        RA proper motion (radians/year)       
    pd : float
        Dec proper motion (radians/year)    
    px : float
        parallax (arcsec)    
    rv : float
        radial velocity (km/s, +ve if receding)     
    astrom : list(30)
        star-independent astrometry parameters     
        [0]>pmt : PM time interval (SSB, Julian years)
        [1-3]>eb : SSB to observer (vector, au)
        [4-6]>eh : Sun to observer (unit vector)
        [7]>em : distance from Sun to observer (au)
        [8-10]>v : barycentric observer velocity (vector, c)
        [11]>bm1 : sqrt(1-|v|^2): reciprocal of Lorenz factor
        [12-20]>bpn : bias-precession-nutation matrix
        [21]>along : adjusted longitude (radians)
        [22]>xpl : polar motion xp wrt local meridian (radians)
        [23]>ypl : polar motion yp wrt local meridian (radians)
        [24]>sphi : sine of geodetic latitude
        [25]>cphi : cosine of geodetic latitude
        [26]>diurab : magnitude of diurnal aberration vector
        [27]>eral : "local" Earth rotation angle (radians)
        [28]>refa : refraction constant A (radians)
        [29]>refb : refraction constant B (radians)

    Returns
    -------
    ri : float
        CIRS geocentric RA,Dec (radians)    
    di : float
        CIRS geocentric RA,Dec (radians)    

    '''
    #调用pymPmpx获得恒星BCRS的坐标方向
    B=[ASTROM[1],ASTROM[2],ASTROM[3]]
    PCO=pymPmpx(RC,DC,PR,PD,PX,RV, ASTROM[0], B)
    
    #调用pymLdsun，处理太阳对光线的偏转，给出BCRS真实方向
    C=[ASTROM[4],ASTROM[5],ASTROM[6]]
    PNAT=pymLdsun(PCO,C,ASTROM[7])
    
    #调用pymAb，处理光行差，给出GCRS方向.
    D=[ASTROM[8],ASTROM[9],ASTROM[10]]
    PPR=pymAb(PNAT,D,ASTROM[7],ASTROM[11])
    
    #调用pymRxp，处理极移进动章动，给出CIRS方向
    F=[[ASTROM[12],ASTROM[13],ASTROM[14]],
       [ASTROM[15],ASTROM[16],ASTROM[17]],
       [ASTROM[18],ASTROM[19],ASTROM[20]]]
    PI=pymRxp(F,PPR)
    
    #CIRS RA,Dec.
    W,DI=pymC2s(PI)
    RI=pymAnp(W)
    
    return(RI,DI)


def pymAtci13(RC,DC,PR,PD,PX,RV,DATE1,DATE2):
    '''
    Transform ICRS star data, epoch J2000.0, to CIRS.

    Parameters
    ----------
    rc : float
        ICRS right ascension at J2000.0 (radians)    
    dc : float
        ICRS declination at J2000.0 (radians)    
    pr : float
        RA proper motion (radians/year)       
    pd : float
        Dec proper motion (radians/year)    
    px : float
        parallax (arcsec)    
    rv : float
        radial velocity (km/s, +ve if receding)     
    date1 : float
        TDB as a 2-part Julian Date    
    date2 : float
        TDB as a 2-part Julian Date    

    Returns
    -------
    ri : float
        CIRS geocentric RA,Dec (radians)    
    di : float
        CIRS geocentric RA,Dec (radians)    
    eo : float
        equation of the origins (ERA-GST)

    '''
    ASTROM=[0 for i in range(30)]
    
    #调用pymApci13获得与恒星无关的天体测量参数ASTROM
    ASTROM,EO=pymApci13(DATE1,DATE2,ASTROM)
    
    #ICRS (epoch J2000.0) 到 CIRS的位置变换.
    RI,DI=pymAtciq(RC,DC,PR,PD,PX,RV,ASTROM)
    
    return(RI,DI,EO)


def pymAtciqn (RC,DC,PR,PD,PX,RV,ASTROM,N,B):
    '''
    Quick ICRS, epoch J2000.0, to CIRS transformation, given precomputed
    star-independent astrometry parameters plus a list of light-
    deflecting bodies.
    
    Use of this function is appropriate when efficiency is important and
    where many star positions are to be transformed for one date.  The
    star-independent parameters can be obtained by calling one of the
    functions iauApci[13], iauApcg[13], iauApco[13] or iauApcs[13].
    
    If the only light-deflecting body to be taken into account is the
    Sun, the iauAtciq function can be used instead.  If in addition the
    parallax and proper motions are zero, the iauAtciqz function can be
    used.

    Parameters
    ----------
    rc : float
        ICRS right ascension at J2000.0 (radians)    
    dc : float
        ICRS declination at J2000.0 (radians)    
    pr : float
        RA proper motion (radians/year)       
    pd : float
        Dec proper motion (radians/year)    
    px : float
        parallax (arcsec)    
    rv : float
        radial velocity (km/s, +ve if receding)     
    astrom : list(30)
        star-independent astrometry parameters     
        [0]>pmt : PM time interval (SSB, Julian years)
        [1-3]>eb : SSB to observer (vector, au)
        [4-6]>eh : Sun to observer (unit vector)
        [7]>em : distance from Sun to observer (au)
        [8-10]>v : barycentric observer velocity (vector, c)
        [11]>bm1 : sqrt(1-|v|^2): reciprocal of Lorenz factor
        [12-20]>bpn : bias-precession-nutation matrix
        [21]>along : adjusted longitude (radians)
        [22]>xpl : polar motion xp wrt local meridian (radians)
        [23]>ypl : polar motion yp wrt local meridian (radians)
        [24]>sphi : sine of geodetic latitude
        [25]>cphi : cosine of geodetic latitude
        [26]>diurab : magnitude of diurnal aberration vector
        [27]>eral : "local" Earth rotation angle (radians)
        [28]>refa : refraction constant A (radians)
        [29]>refb : refraction constant B (radians)
    n : int
        number of bodies    
    b: list(n,8)
        data for each of the n bodies    
        [i][0]>bm :  : mass of the body (solar masses)    
        [i][1]>dl : deflection limiter
        [i][2-7]>pv : barycentric PV of the body (au, au/day)

    Returns
    -------
    ri : float
        CIRS geocentric RA,Dec (radians)    
    di : float
        CIRS geocentric RA,Dec (radians)  

    '''
    #调用pymPmpx,给出在质心天球参考系(BCRS)下的坐标方向
    G=[ASTROM[1],ASTROM[2],ASTROM[3]]
    PCO=pymPmpx(RC,DC,PR,PD,PX,RV,ASTROM[0], G)
    
    #调用pymLdn,处理多个天体对光线的偏转，给出BCRS的真实方向.
    PNAT=pymLdn(N,B,G,PCO)
    
    #调用pymAb，处理光行差，给出GCRS方向.
    D=[ASTROM[8],ASTROM[9],ASTROM[10]]
    PPR=pymAb(PNAT,D,ASTROM[7],ASTROM[11])
    
    #调用pymRxp，处理极移进动章动，给出CIRS方向
    F=[[ASTROM[12],ASTROM[13],ASTROM[14]],
       [ASTROM[15],ASTROM[16],ASTROM[17]],
       [ASTROM[18],ASTROM[19],ASTROM[20]]]
    PI=pymRxp(F,PPR)
    
    #CIRS RA,Dec.
    W,DI=pymC2s(PI)
    RI=pymAnp(W)
    
    return(RI,DI)


def pymS2c(THETA,PHI):
    '''
    Convert spherical coordinates to Cartesian.

    Parameters
    ----------
    theta : float
        longitude angle (radians)    
    phi : float
        latitude angle (radians)

    Returns
    -------
    c : numpy.matrix
        direction cosines

    '''
    C=[0,0,0]
    
    CP=ma.cos(PHI)
    C[0]=ma.cos(THETA)*CP
    C[1]=ma.sin(THETA)*CP
    C[2]=ma.sin(PHI)
    
    return(C)


def pymAtciqz (RC,DC,ASTROM):
    '''
    Quick ICRS to CIRS transformation, given precomputed star-
    independent astrometry parameters, and assuming zero parallax and
    proper motion.
    
    Use of this function is appropriate when efficiency is important and
    where many star positions are to be transformed for one date.  The
    star-independent parameters can be obtained by calling one of the
    functions iauApci[13], iauApcg[13], iauApco[13] or iauApcs[13].
    
    The corresponding function for the case of non-zero parallax and
    proper motion is iauAtciq.

    Parameters
    ----------
    rc : float
        ICRS right ascension at J2000.0 (radians)    
    dc : float
        ICRS declination at J2000.0 (radians)    
    astrom : list(30)
        star-independent astrometry parameters     
        [0]>pmt : PM time interval (SSB, Julian years)
        [1-3]>eb : SSB to observer (vector, au)
        [4-6]>eh : Sun to observer (unit vector)
        [7]>em : distance from Sun to observer (au)
        [8-10]>v : barycentric observer velocity (vector, c)
        [11]>bm1 : sqrt(1-|v|^2): reciprocal of Lorenz factor
        [12-20]>bpn : bias-precession-nutation matrix
        [21]>along : adjusted longitude (radians)
        [22]>xpl : polar motion xp wrt local meridian (radians)
        [23]>ypl : polar motion yp wrt local meridian (radians)
        [24]>sphi : sine of geodetic latitude
        [25]>cphi : cosine of geodetic latitude
        [26]>diurab : magnitude of diurnal aberration vector
        [27]>eral : "local" Earth rotation angle (radians)
        [28]>refa : refraction constant A (radians)
        [29]>refb : refraction constant B (radians)

    Returns
    -------
    ri : float
        CIRS geocentric RA,Dec (radians)    
    di : float
        CIRS geocentric RA,Dec (radians)  

    '''
    #转换为BCRS坐标方向，单位向量
    PCO=pymS2c(RC,DC)
    
    #调用pymLdsun，处理太阳对光线的偏转，给出BCRS真实方向
    C=[ASTROM[4],ASTROM[5],ASTROM[6]]
    PNAT=pymLdsun(PCO,C,ASTROM[7])
    
    #调用pymAb，处理光行差，给出GCRS方向.
    D=[ASTROM[8],ASTROM[9],ASTROM[10]]
    PPR=pymAb(PNAT,D,ASTROM[7],ASTROM[11])
    
    #调用pymRxp，处理极移进动章动，给出CIRS方向
    F=[[ASTROM[12],ASTROM[13],ASTROM[14]],
       [ASTROM[15],ASTROM[16],ASTROM[17]],
       [ASTROM[18],ASTROM[19],ASTROM[20]]]
    PI=pymRxp(F,PPR)
    
    #CIRS RA,Dec.
    W,DI=pymC2s(PI)
    RI=pymAnp(W)
    
    return(RI,DI)


def pymAtioq(RI,DI,ASTROM):
    '''
    Quick CIRS to observed place transformation.
    
    Use of this function is appropriate when efficiency is important and
    where many star positions are all to be transformed for one date.
    The star-independent astrometry parameters can be obtained by
    calling iauApio[13] or iauApco[13].

    Parameters
    ----------
    ri : float
        CIRS geocentric RA,Dec (radians)    
    di : float
        CIRS geocentric RA,Dec (radians)     
    astrom : list(30)
        star-independent astrometry parameters     
        [0]>pmt : PM time interval (SSB, Julian years)
        [1-3]>eb : SSB to observer (vector, au)
        [4-6]>eh : Sun to observer (unit vector)
        [7]>em : distance from Sun to observer (au)
        [8-10]>v : barycentric observer velocity (vector, c)
        [11]>bm1 : sqrt(1-|v|^2): reciprocal of Lorenz factor
        [12-20]>bpn : bias-precession-nutation matrix
        [21]>along : adjusted longitude (radians)
        [22]>xpl : polar motion xp wrt local meridian (radians)
        [23]>ypl : polar motion yp wrt local meridian (radians)
        [24]>sphi : sine of geodetic latitude
        [25]>cphi : cosine of geodetic latitude
        [26]>diurab : magnitude of diurnal aberration vector
        [27]>eral : "local" Earth rotation angle (radians)
        [28]>refa : refraction constant A (radians)
        [29]>refb : refraction constant B (radians)

    Returns
    -------
    aob : float
        observed azimuth (radians: N=0,E=90)     
    zob : float
        observed zenith distance (radians)    
    hob : float
        observed hour angle (radians)    
    dob : float
        observed declination (radians)    
    rob : float
        observed right ascension (CIO-based, radians)   

    '''
    #用于折射的最小正弦和余弦高度
    SELMIN=0.05
    CELMIN=1e-6

    #CIRS RA,Dec 到笛卡尔-HA,Dec的转换.
    V=pymS2c(RI-ASTROM[27],DI)    
    X=V[0]
    Y=V[1]
    Z=V[2]
    
    #极移
    SX=ma.sin(ASTROM[22])
    CX=ma.cos(ASTROM[22])
    SY=ma.sin(ASTROM[23])
    CY=ma.cos(ASTROM[23])
    XHD=CX*X+SX*Z
    YHD=SX*SY*X+CY*Y-CX*SY*Z
    ZHD=-SX*CY*X+SY*Y+CX*CY*Z

    #周日光行差
    F=(1.0-ASTROM[26]*YHD )
    XHDT=F*XHD
    YHDT=F*(YHD+ASTROM[26])
    ZHDT=F*ZHD
    
    #由笛卡尔-HA,Dec到笛卡尔Az,El的转换 (S=0,E=90).
    XAET=ASTROM[24]*XHDT-ASTROM[25]*ZHDT
    YAET=YHDT
    ZAET=ASTROM[25]*XHDT+ASTROM[24]*ZHDT
    
    #方位角 (N=0,E=90).
    if (XAET!=0.0)|(YAET!=0.0):
        AZOBS=ma.atan2(YAET,-XAET)
    else:
        AZOBS=0.0
    
    # 折射
    
    #在避免出错的情况下，给出正弦和余弦高度
    R=max(ma.sqrt(XAET*XAET+YAET*YAET),CELMIN)
    Z=max(ZAET,SELMIN)
    
    #A*tan(z)+B*tan^3(z) 折射模型.
    TZ=R/Z
    W=ASTROM[29]*TZ*TZ
    DEL=(ASTROM[28]+W)*TZ/(1.0+(ASTROM[28]+3.0*W)/(Z*Z))

    #应用变化，给出观察向量。
    COSDEL=1.0-DEL*DEL/2.0
    F=COSDEL-DEL*Z/R
    XAEO=XAET*F
    YAEO=YAET*F
    ZAEO=COSDEL*ZAET+DEL*R
    
    # 获得天顶距
    ZDOBS=ma.atan2(ma.sqrt(XAEO*XAEO+YAEO*YAEO),ZAEO)
    
    #Az/El 向量到HA,Dec向量的转换 (皆为右手系).
    V[0]=ASTROM[24]*XAEO+ASTROM[25]*ZAEO
    V[1]=YAEO
    V[2]=-ASTROM[25]*XAEO+ASTROM[24]*ZAEO
    
    #从笛卡尔到球面 -HA,Dec.
    HMOBS,DCOBS=pymC2s(V)

    #赤经(根据CIO).
    RAOBS=ASTROM[27]+HMOBS

    #最终结果
    AOB=pymAnp(AZOBS)
    ZOB=ZDOBS
    HOB=-HMOBS
    DOB=DCOBS
    ROB=pymAnp(RAOBS)
    
    return(AOB,ZOB,HOB,DOB,ROB)


def pymAtco13(RC,DC,PR,PD,PX,RV,UTC1,UTC2,DUT1,ELONG,
                PHI,HM,XP,YP,PHPA,TC,RH,WL):
    '''
    ICRS RA,Dec to observed place.  The caller supplies UTC, site
    coordinates, ambient air conditions and observing wavelength.

    Parameters
    ----------
    rc : float
        ICRS right ascension at J2000.0 (radians)    
    dc : float
        ICRS declination at J2000.0 (radians)    
    pr : float
        RA proper motion (radians/year)       
    pd : float
        Dec proper motion (radians/year)    
    px : float
        parallax (arcsec)    
    rv : float
        radial velocity (km/s, +ve if receding)     
    utc1 : float
        UTC as a 2-part Julian Date    
    utc2 : float
        UTC as a 2-part Julian Date    
    dut1 : float
        UT1-UTC (seconds)    
    elong : float
        longitude (radians, east +ve)    
    phi : float
        geodetic latitude (radians)    
    hm : float
        height above ellipsoid (m, geodetic)    
    xp : float
        polar motion coordinates (radians)    
    yp : float
        polar motion coordinates (radians)    
    phpa : float
        pressure at the observer (hPa = mB)    
    tc : float
        ambient temperature at the observer (deg C)    
    rh : float
        relative humidity at the observer (range 0-1)    
    wl : float
        wavelength (micrometers)    

    Returns
    -------
    aob : float
        observed azimuth (radians: N=0,E=90)     
    zob : float
        observed zenith distance (radians)    
    hob : float
        observed hour angle (radians)    
    dob : float
        observed declination (radians)    
    rob : float
        observed right ascension (CIO-based, radians)    
    eo : float
        equation of the origins (ERA-GST)
    J : ValueError
        1: 'dubious year',
       -1: 'unacceptable date',
    '''
    JS=0
    i=1    
    while i<2:

        #调用pymApco13获得与恒星无关的天体测量参数
        ASTROM,EO,JS=pymApco13(UTC1,UTC2,DUT1,ELONG,
                                PHI,HM,XP,YP,PHPA,TC,RH,WL)
        if (JS<0):
            print('ERROR',JS)
            break
        
        #从ICRS变换到CIRS.
        RI,DI=pymAtciq(RC,DC,PR,PD,PX,RV,ASTROM)
        
        #从CIRS变换到observed.
        AOB,ZOB,HOB,DOB,ROB=pymAtioq(RI,DI,ASTROM)
    
        i+=1
    J=JS

    return(AOB,ZOB,HOB,DOB,ROB,EO,J)


def pymAticq(RI,DI,ASTROM):
    '''
    Quick CIRS RA,Dec to ICRS astrometric place, given the star-
    independent astrometry parameters.

    Parameters
    ----------
    ri : float
        CIRS geocentric RA,Dec (radians)    
    di : float
        CIRS geocentric RA,Dec (radians)  
    astrom : list(30)
        star-independent astrometry parameters     
        [0]>pmt : PM time interval (SSB, Julian years)
        [1-3]>eb : SSB to observer (vector, au)
        [4-6]>eh : Sun to observer (unit vector)
        [7]>em : distance from Sun to observer (au)
        [8-10]>v : barycentric observer velocity (vector, c)
        [11]>bm1 : sqrt(1-|v|^2): reciprocal of Lorenz factor
        [12-20]>bpn : bias-precession-nutation matrix
        [21]>along : adjusted longitude (radians)
        [22]>xpl : polar motion xp wrt local meridian (radians)
        [23]>ypl : polar motion yp wrt local meridian (radians)
        [24]>sphi : sine of geodetic latitude
        [25]>cphi : cosine of geodetic latitude
        [26]>diurab : magnitude of diurnal aberration vector
        [27]>eral : "local" Earth rotation angle (radians)
        [28]>refa : refraction constant A (radians)
        [29]>refb : refraction constant B (radians)

    Returns
    -------
    rc : float
        ICRS right ascension at J2000.0 (radians)    
    dc : float
        ICRS declination at J2000.0 (radians)    

    '''
    BEFORE,PNAT,PCO=[0,0,0],[0,0,0],[0,0,0]
    #CIRS RA,Dec转换到笛卡尔坐标.
    PI=pymS2c(RI,DI)
    
    #调用pymRxp，处理极移进动章动，给出GCRS方向
    F=[[ASTROM[12],ASTROM[13],ASTROM[14]],
       [ASTROM[15],ASTROM[16],ASTROM[17]],
       [ASTROM[18],ASTROM[19],ASTROM[20]]]
    PPR=pymTrxp(F,PI)
    
    #处理光行差，给出GCRS方向
    D=pymZp()
    E=[ASTROM[8],ASTROM[9],ASTROM[10]]
    for j in range(2):
        R2=0.0
        for i in range(3):
            W=PPR[i]-D[i]
            BEFORE[i]=W
            R2=R2+W*W
        R=ma.sqrt(R2)
        for i in range(3):
            BEFORE[i]=BEFORE[i]/R
        AFTER=pymAb(BEFORE,E,ASTROM[7],ASTROM[11])
        R2=0.0
        for i in range(3):
            D[i]=AFTER[i]-BEFORE[i]
            W=PPR[i]-D[i]
            PNAT[i]=W
            R2=R2+W*W
        R=ma.sqrt(R2)
        for i in range(3):
            PNAT[i]=PNAT[i]/R
    
    #调用pymLdsun，处理太阳对光线的偏转，给出BCRS真实方向
    D=pymZp()
    C=[ASTROM[4],ASTROM[5],ASTROM[6]]
    for j in range(5):
        R2=0.0
        for i in range(3):
            W=PNAT[i]-D[i]
            BEFORE[i]=W
            R2=R2+W*W
        R=ma.sqrt(R2)
        for i in range(3):
            BEFORE[i]=BEFORE[i]/R
        AFTER=pymLdsun(BEFORE,C,ASTROM[7])
        R2=0.0
        for i in range(3):
            D[i]=AFTER[i]-BEFORE[i]
            W=PNAT[i]-D[i]
            PCO[i]=W
            R2=R2+W*W
        R=ma.sqrt(R2)
        for i in range(3):
            PCO[i]=PCO[i]/R
    
    #ICRS的天体测量RA,Dec.
    W,DC=pymC2s(PCO)
    RC=pymAnp(W)  
        
    return(RC,DC)


def pymAtic13(RI,DI,DATE1,DATE2):
    '''
    Transform star RA,Dec from geocentric CIRS to ICRS astrometric.

    Parameters
    ----------
    ri : float
        CIRS geocentric RA,Dec (radians)    
    di : float
        CIRS geocentric RA,Dec (radians)  
    date1 : float
        TDB as a 2-part Julian Date    
    date2 : float
        TDB as a 2-part Julian Date    

    Returns
    -------
    rc : float
        ICRS right ascension at J2000.0 (radians)    
    dc : float
        ICRS declination at J2000.0 (radians)    
    eo : float
        equation of the origins (ERA-GST)

    '''
    ASTROM=[0 for i in range(30)]
    #调用pymApci13获得与恒星无关的天体测量参数
    ASTROM, EO=pymApci13(DATE1,DATE2,ASTROM)
    
    #从CIRS到ICRS转换
    RC,DC=pymAticq(RI,DI,ASTROM)
    
    return(RC,DC,EO)


def pymAticqn(RI,DI,ASTROM,N,B):
    '''
    Quick CIRS to ICRS astrometric place transformation, given the star-
    independent astrometry parameters plus a list of light-deflecting
    bodies.
    
    Use of this function is appropriate when efficiency is important and
    where many star positions are all to be transformed for one date.
    The star-independent astrometry parameters can be obtained by
    calling one of the functions iauApci[13], iauApcg[13], iauApco[13]
    or iauApcs[13].
    
    If the only light-deflecting body to be taken into account is the
    Sun, the iauAticq function can be used instead.

    Parameters
    ----------
    ri : float
        CIRS geocentric RA,Dec (radians)    
    di : float
        CIRS geocentric RA,Dec (radians)    
    astrom : list(30)
        star-independent astrometry parameters     
        [0]>pmt : PM time interval (SSB, Julian years)
        [1-3]>eb : SSB to observer (vector, au)
        [4-6]>eh : Sun to observer (unit vector)
        [7]>em : distance from Sun to observer (au)
        [8-10]>v : barycentric observer velocity (vector, c)
        [11]>bm1 : sqrt(1-|v|^2): reciprocal of Lorenz factor
        [12-20]>bpn : bias-precession-nutation matrix
        [21]>along : adjusted longitude (radians)
        [22]>xpl : polar motion xp wrt local meridian (radians)
        [23]>ypl : polar motion yp wrt local meridian (radians)
        [24]>sphi : sine of geodetic latitude
        [25]>cphi : cosine of geodetic latitude
        [26]>diurab : magnitude of diurnal aberration vector
        [27]>eral : "local" Earth rotation angle (radians)
        [28]>refa : refraction constant A (radians)
        [29]>refb : refraction constant B (radians)
    n : int
        number of bodies    
    b: list(n,8)
        data for each of the n bodies    
        [i][0]>bm :  : mass of the body (solar masses)    
        [i][1]>dl : deflection limiter
        [i][2-7]>pv : barycentric PV of the body (au, au/day)

    Returns
    -------
    rc : float
        ICRS right ascension at J2000.0 (radians)    
    dc : float
        ICRS declination at J2000.0 (radians)    

    '''
    
    BEFORE,PNAT,PCO=[0,0,0],[0,0,0],[0,0,0]
    #CIRS RA,Dec转换到笛卡尔坐标.
    PI=pymS2c(RI,DI)
    
    #调用pymRxp，处理极移进动章动，给出GCRS方向
    F=[[ASTROM[12],ASTROM[13],ASTROM[14]],
       [ASTROM[15],ASTROM[16],ASTROM[17]],
       [ASTROM[18],ASTROM[19],ASTROM[20]]]
    PPR=pymTrxp(F,PI)
    
    #处理光行差，给出GCRS方向
    D=pymZp()
    E=[ASTROM[8],ASTROM[9],ASTROM[10]]
    for j in range(2):
        R2=0.0
        for i in range(3):
            W=PPR[i]-D[i]
            BEFORE[i]=W
            R2=R2+W*W
        R=ma.sqrt(R2)
        for i in range(3):
            BEFORE[i]=BEFORE[i]/R
        AFTER=pymAb(BEFORE,E,ASTROM[7],ASTROM[11])
        R2=0.0
        for i in range(3):
            D[i]=AFTER[i]-BEFORE[i]
            W=PPR[i]-D[i]
            PNAT[i]=W
            R2=R2+W*W
        R=ma.sqrt(R2)
        for i in range(3):
            PNAT[i]=PNAT[i]/R
       
    #调用pymLdn，处理多个天体对光线的偏转，给出BCRS真实方向
    D=pymZp()
    G=[ASTROM[1],ASTROM[2],ASTROM[3]]
    for j in range(5):
        R2=0.0
        for i in range(3):
            W=PNAT[i]-D[i]
            BEFORE[i]=W
            R2=R2+W*W
        R=ma.sqrt(R2)
        for i in range(3):
            BEFORE[i]=BEFORE[i]/R
        AFTER=pymLdn(N,B,G,BEFORE)
        R2=0.0
        for i in range(3):
            D[i]=AFTER[i]-BEFORE[i]
            W=PNAT[i]-D[i]
            PCO[i]=W
            R2=R2+W*W
        R=ma.sqrt(R2)
        for i in range(3):
            PCO[i]=PCO[i]/R
    
    #ICRS astrometric RA,Dec.
    W,DC=pymC2s(PCO)
    RC=pymAnp(W)    
    
    return(RC,DC)


def pymAtio13(RI,DI,UTC1,UTC2,DUT1,ELONG,PHI,HM,XP,YP,PHPA,TC,RH,WL):
    '''
    CIRS RA,Dec to observed place.  The caller supplies UTC, site
    coordinates, ambient air conditions and observing wavelength.

    Parameters
    ----------
    ri : float
        CIRS geocentric RA,Dec (radians)    
    di : float
        CIRS geocentric RA,Dec (radians)     
    utc1 : float
        UTC as a 2-part Julian Date    
    utc2 : float
        UTC as a 2-part Julian Date    
    dut1 : float
        UT1-UTC (seconds)    
    elong : float
        longitude (radians, east +ve)    
    phi : float
        geodetic latitude (radians)    
    hm : float
        height above ellipsoid (m, geodetic)    
    xp : float
        polar motion coordinates (radians)    
    yp : float
        polar motion coordinates (radians)    
    phpa : float
        pressure at the observer (hPa = mB)    
    tc : float
        ambient temperature at the observer (deg C)    
    rh : float
        relative humidity at the observer (range 0-1)    
    wl : float
        wavelength (micrometers)    

    Returns
    -------
    aob : float
        observed azimuth (radians: N=0,E=90)     
    zob : float
        observed zenith distance (radians)    
    hob : float
        observed hour angle (radians)    
    dob : float
        observed declination (radians)    
    rob : float
        observed right ascension (CIO-based, radians)   
    J : ValueError
        1: 'dubious year',
       -1: 'unacceptable date',
    '''
    ASTROM=[0 for i in range(30)]
    
    i=1
    while i<2:
        #获得从CIRS转换到观察者的与恒星无关的天体测量参数
        ASTROM,JS=pymApio13(UTC1,UTC2,DUT1,ELONG,PHI,HM, 
                              XP,YP,PHPA,TC,RH,WL,ASTROM)
        
        #当UTC错误时终止函数.
        if (JS<0):
            print('ERROR',JS)
            break
        
        #从CIRS转换到观察者
        AOB,ZOB,HOB,DOB,ROB=pymAtioq(RI,DI,ASTROM) 

        i+=1
    
    J=JS
    
    return(AOB,ZOB,HOB,DOB,ROB,J)


def pymAtoiq(TYPE,OB1,OB2,ASTROM):
    '''
    Quick observed place to CIRS, given the star-independent astrometry
    parameters.
    
    Use of this function is appropriate when efficiency is important and
    where many star positions are all to be transformed for one date.
    The star-independent astrometry parameters can be obtained by
    calling iauApio[13] or iauApco[13].

    Parameters
    ----------
    stype : str
        type of coordinates - "R", "H" or "A"     
    ob1 : float
        observed Az, HA or RA (radians; Az is N=0,E=90)    
    ob2 : float
        observed ZD or Dec (radians)    
    astrom : list(30)
        star-independent astrometry parameters     
        [0]>pmt : PM time interval (SSB, Julian years)
        [1-3]>eb : SSB to observer (vector, au)
        [4-6]>eh : Sun to observer (unit vector)
        [7]>em : distance from Sun to observer (au)
        [8-10]>v : barycentric observer velocity (vector, c)
        [11]>bm1 : sqrt(1-|v|^2): reciprocal of Lorenz factor
        [12-20]>bpn : bias-precession-nutation matrix
        [21]>along : adjusted longitude (radians)
        [22]>xpl : polar motion xp wrt local meridian (radians)
        [23]>ypl : polar motion yp wrt local meridian (radians)
        [24]>sphi : sine of geodetic latitude
        [25]>cphi : cosine of geodetic latitude
        [26]>diurab : magnitude of diurnal aberration vector
        [27]>eral : "local" Earth rotation angle (radians)
        [28]>refa : refraction constant A (radians)
        [29]>refb : refraction constant B (radians)

    Returns
    -------
    ri : float
        CIRS geocentric RA,Dec (radians)    
    di : float
        CIRS geocentric RA,Dec (radians)     

    '''
    V=[0,0,0]
    #用于折射的最小正弦和余弦高度
    SELMIN=0.05
    
    #参考架类型
    C=TYPE
    
    #坐标
    C1=OB1
    C2=OB2
    
    #纬度的正、余弦值
    SPHI=ASTROM[24]
    CPHI=ASTROM[25]
    
    #标准化坐标系类型
    if (C=='r')|(C=='R'):
        C='R'
    elif (C=='h')|(C=='H'):
        C='H'
    else:
        C='A'
        
    #如果是地平坐标系，转换为笛卡尔坐标 (S=0,E=90).
    if (C=='A'):
        CE=ma.sin(C2)
        XAEO=-ma.cos(C1)*CE
        YAEO=ma.sin(C1)*CE
        ZAEO=ma.cos(C2)
    
    else:
        #如果是赤道坐标系，利用ASTROM[27]，地球自转角，转换为时角
        if (C=='R'):
            C1=ASTROM[27]-C1
        #转换为笛卡尔时角 -HA,DeC.
        V=pymS2c(-C1,C2)
        XMHDO=V[0]
        YMHDO=V[1]
        ZMHDO=V[2]
        
        #转换为笛卡尔地平 Az,El (S=0,E=90).
        XAEO=SPHI*XMHDO-CPHI*ZMHDO
        YAEO=YMHDO
        ZAEO=CPHI*XMHDO+SPHI*ZMHDO
        
    #方位角 (S=0,E=90).
    if (XAEO!=0.0)|(YAEO!=0.0):
        AZ=ma.atan2(YAEO,XAEO)
    else:
        AZ=0.0
    
    #观测者天顶距的正弦值，以及天顶距
    SZ=ma.sqrt(XAEO*XAEO+YAEO*YAEO)
    ZDO=ma.atan2(SZ,ZAEO)
    
    #折射

    # 采用双常数模型的快速算法。
    REFA=ASTROM[28]
    REFB=ASTROM[29]
    ZAEO=max(ZAEO,SELMIN)
    TZ=SZ/ZAEO
    DREF=(REFA+REFB*TZ*TZ)*TZ
    ZDT=ZDO+DREF
    
    #转换为笛卡尔地平 Az,ZD.
    CE=ma.sin(ZDT)
    XAET=ma.cos(AZ)*CE
    YAET=ma.sin(AZ)*CE
    ZAET=ma.cos(ZDT)
    
    #笛卡尔地平Az,ZD到笛卡尔时角-HA,DeC.
    XMHDA=SPHI*XAET+CPHI*ZAET
    YMHDA=YAET
    ZMHDA=-CPHI*XAET+SPHI*ZAET

    #周日光行差
    F=(1.0+ASTROM[26]*YMHDA)
    XHD=F*XMHDA
    YHD=F*(YMHDA-ASTROM[26])
    ZHD=F*ZMHDA

    #极移
    SX=ma.sin(ASTROM[22])
    CX=ma.cos(ASTROM[22])
    SY=ma.sin(ASTROM[23])
    CY=ma.cos(ASTROM[23])
    V[0]=CX*XHD+SX*SY*YHD-SX*CY*ZHD
    V[1]=CY*YHD+SY*ZHD
    V[2]=SX*XHD-CX*SY*YHD+CX*CY*ZHD
    
    #转换为球面赤道 -HA,Dec.
    HMA,DI=pymC2s(V)
    
    #赤经
    RI=pymAnp(ASTROM[27]+HMA)
    
    return(RI,DI)


def pymAtoi13(TYPE,OB1,OB2,UTC1,UTC2,DUT1,ELONG,PHI,HM,XP,YP,PHPA,TC,RH,WL):
    '''
    Observed place to CIRS.  The caller supplies UTC, site coordinates,
    ambient air conditions and observing wavelength.

    Parameters
    ----------
    stype : str
        type of coordinates - "R", "H" or "A"     
    ob1 : float
        observed Az, HA or RA (radians; Az is N=0,E=90)    
    ob2 : float
        observed ZD or Dec (radians)    
    utc1 : float
        UTC as a 2-part Julian Date    
    utc2 : float
        UTC as a 2-part Julian Date    
    dut1 : float
        UT1-UTC (seconds)    
    elong : float
        longitude (radians, east +ve)    
    phi : float
        geodetic latitude (radians)    
    hm : float
        height above ellipsoid (m, geodetic)    
    xp : float
        polar motion coordinates (radians)    
    yp : float
        polar motion coordinates (radians)    
    phpa : float
        pressure at the observer (hPa = mB)    
    tc : float
        ambient temperature at the observer (deg C)    
    rh : float
        relative humidity at the observer (range 0-1)    
    wl : float
        wavelength (micrometers)  
    
    Returns
    -------
    ri : float
        CIRS geocentric RA,Dec (radians)    
    di : float
        CIRS geocentric RA,Dec (radians)     
    J : ValueError
        1: 'dubious year',
       -1: 'unacceptable date',
    '''
    ASTROM=[0 for i in range(30)]
    i=1
    while i<2:
        
        #获得从CIRS转换到观察者的与恒星无关的天体测量参数
        ASTROM,JS=pymApio13(UTC1,UTC2,DUT1,ELONG,PHI,HM,XP,\
                             YP,PHPA,TC,RH,WL,ASTROM)
    
        #当UTC错误时终止函数.
        if (JS<0):
            print('ERROR',JS)
            break
        
        #从观察者转换到CIRS
        RI,DI=pymAtoiq(TYPE,OB1,OB2,ASTROM)
        
        i+=1
    
    J=JS
    
    return(RI,DI,J)


def pymAtoc13(TYPE,OB1,OB2,UTC1,UTC2,DUT1,ELONG,PHI,HM,XP,YP,PHPA,TC,RH,WL,):
    '''
    Observed place at a groundbased site to to ICRS astrometric RA,Dec.
    The caller supplies UTC, site coordinates, ambient air conditions
    and observing wavelength.

    Parameters
    ----------
    stype : str
        type of coordinates - "R", "H" or "A"     
    ob1 : float
        observed Az, HA or RA (radians; Az is N=0,E=90)    
    ob2 : float
        observed ZD or Dec (radians)    
    utc1 : float
        UTC as a 2-part Julian Date    
    utc2 : float
        UTC as a 2-part Julian Date    
    dut1 : float
        UT1-UTC (seconds)    
    elong : float
        longitude (radians, east +ve)    
    phi : float
        geodetic latitude (radians)    
    hm : float
        height above ellipsoid (m, geodetic)    
    xp : float
        polar motion coordinates (radians)    
    yp : float
        polar motion coordinates (radians)    
    phpa : float
        pressure at the observer (hPa = mB)    
    tc : float
        ambient temperature at the observer (deg C)    
    rh : float
        relative humidity at the observer (range 0-1)    
    wl : float
        wavelength (micrometers)  

    Returns
    -------
    rc : float
        ICRS right ascension at J2000.0 (radians)    
    dc : float
        ICRS declination at J2000.0 (radians)    
    J : ValueError
        1: 'dubious year',
       -1: 'unacceptable date',
    '''
    i=1
    while i<2:
        #获得从ICRS转换到观察者的与恒星无关的天体测量参数
        ASTROM,EO,JS=pymApco13(UTC1,UTC2,DUT1,ELONG,PHI,HM,XP,YP,PHPA,TC,RH,WL)
        
        #当UTC错误时终止函数.
        if (JS<0):
            print('ERROR',JS)
            break
    
        #从观察者转换到CIRS
        RI,DI=pymAtoiq(TYPE,OB1,OB2,ASTROM)
        
        #从CIRS转换到ICRS
        RC,DC=pymAticq(RI,DI,ASTROM)
    
        i+=1
    J=JS
    
    return(RC,DC,J)


def pymBi00():
    '''
    Frame bias components of IAU 2000 precession-nutation models;  part
    of the Mathews-Herring-Buffett (MHB2000) nutation series, with
    additions.

    Returns
    -------
    dpsibi : float
        longitude and obliquity corrections    
    depsbi : float
        longitude and obliquity corrections    
    dra : flaot
        the ICRS RA of the J2000.0 mean equinox

    '''
    #1角秒对应的弧度
    DAS2R=4.848136811095359935899141e-6
    
    #在经度和倾角上的框架极移改正
    DPBIAS=-0.041775*DAS2R
    DEBIAS=-0.0068192*DAS2R
    
    #J2000.0平春分点的ICRS赤经(Chapront et al., 2002)
    DRA0=-0.0146*DAS2R
    
    #输出结果.
    DPSIBI=DPBIAS
    DEPSBI=DEBIAS
    DRA=DRA0
    
    return(DPSIBI,DEPSBI,DRA)


def pymPr00(DATE1,DATE2):
    '''
    Precession-rate part of the IAU 2000 precession-nutation models
    (part of MHB2000).

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    dpsipr : float 
        precession corrections    
    depspr : float
        precession corrections

    '''
    #1角秒对应的弧度
    DAS2R=4.848136811095359935899141e-6
    
    #参考历元(J2000.0), JD
    DJ00=2451545

    #1儒略世纪的天数
    DJC=36525.0
    
    # 进动和倾角的改正，单位：弧度/世纪
    PRECOR=-0.29965*DAS2R
    OBLCOR=-0.02524*DAS2R
    
    #与参考历元之间的时间间隔，单位：儒略世纪数.
    T=((DATE1-DJ00)+DATE2)/DJC
    
    #进动速率
    DPSIPR=PRECOR*T
    DEPSPR=OBLCOR*T
    
    return(DPSIPR,DEPSPR)


def pymRxr(A,B):
    '''
    Multiply two r-matrices.

    Parameters
    ----------
    a : list(3,3)
        first r-matrix    
    b : list(3,3)
        second r-matrix

    Returns
    -------
    atb : list(3,3)
        a * b

    '''
    WM=[[0 for i in range(3)] for j in range(3)]
    
    for i in range(3):
        for j in range(3):
            W=0.0
            for k in range(3):
                W=W+A[i][k]*B[k][j]
            WM[i][j]=W
    ATB=pymCr(WM)
    
    return(ATB)


def pymBp00(DATE1,DATE2):
    '''
    Frame bias and precession, IAU 2000.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    rb : list(3,3)
        frame bias matrix     
    rp : list(3,3)   
        precession matrix    
    rbp : list(3,3)
        bias-precession matrix
    '''
    #1角秒对应的弧度
    DAS2R=4.848136811095359935899141e-6
    
    #参考历元(J2000.0), JD
    DJ00=2451545.0

    #1儒略世纪对应的天数
    DJC=36525.0

    #J2000.0 倾角(Lieske et al. 1977)
    EPS0=84381.448e0*DAS2R
    
    #与参考历元之间的时间间隔，单位：儒略世纪数.
    T=((DATE1-DJ00)+DATE2)/DJC
    
    #参考架偏差
    DPSIBI,DEPSBI,DRA0=pymBi00()

    #岁差角(Lieske et al. 1977)
    PSIA77=(5038.7784+(-1.07259+(-0.001147)*T)*T)*T*DAS2R
    OMA77=EPS0+((0.05127+(-0.007726)*T)*T)*T*DAS2R
    CHIA=(10.5526+(-2.38064+(-0.001125)*T)*T)*T*DAS2R
    
    #应用IAU 2000岁差改正.
    DPSIPR, DEPSPR=pymPr00(DATE1, DATE2)
    PSIA=PSIA77+DPSIPR
    OMA=OMA77+DEPSPR
    
    #参考架偏差矩阵: GCRS to J2000.0.
    RBW=pymIr()
    RBW=pymRz(DRA0,RBW)
    RBW=pymRy(DPSIBI*ma.sin(EPS0),RBW)
    RBW=pymRx(-DEPSBI,RBW)
    RB=pymCr(RBW)
    
    #岁差矩阵: J2000.0 to mean of date.
    RP=pymIr()
    RP=pymRx(EPS0,RP)
    RP=pymRz(-PSIA,RP)
    RP=pymRx(-OMA,RP)
    RP=pymRz(CHIA,RP)
    
    #偏差-岁差矩阵: GCRS to mean of date.
    RBP=pymRxr(RP,RBW)
    
    return(RB,RP,RBP)


def pymPmat06 (DATE1,DATE2):
    '''
    Precession matrix (including frame bias) from GCRS to a specified
    date, IAU 2006 model.

    Parameters
    ----------
    date1 : flaot
        TT as a 2-part Julian Date     
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    rbp : list(3,3)
        bias-precession matrix

    '''

    #Fukushima-Williams 偏差-岁差角
    GAMB,PHIB,PSIB,EPSA=pymPfw06(DATE1,DATE2)
    
    #构建矩阵
    RBP=pymFw2m(GAMB, PHIB, PSIB, EPSA)
    
    return(RBP)


def pymBp06(DATE1,DATE2):
    '''
    Frame bias and precession, IAU 2006.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    rb : list(3,3)
        frame bias matrix    
    rp : list(3,3)
        precession matrix    
    rbp : list(3,3)
        bias-precession matrix

    '''
    #约简儒略日零点
    DJM0=2400000.5
    
    #参考历元(J2000.0), 约简儒略日
    DJM00 = 51544.5
    
    #B matrix.
    GAMB,PHIB,PSIB,EPSA=pymPfw06(DJM0,DJM00)
    RB=pymFw2m(GAMB,PHIB,PSIB,EPSA)
    
    #PxB matrix (temporary).
    RBPW=pymPmat06(DATE1, DATE2)
    
    #P matrix.
    RBT=pymTr(RB)
    RP=pymRxr(RBPW,RBT)
    
    #PxB matrix.
    RBP=pymCr(RBPW)
    
    return(RB,RP,RBP)


def pymObl80(DATE1,DATE2):
    '''
    Mean obliquity of the ecliptic, IAU 1980 model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    function value : float
        obliquity of the ecliptic (radians)

    '''
    #1角秒对应的弧度
    DAS2R=4.848136811095359935899141e-6
    
    #参考历元(J2000.0), JD
    DJ00=2451545.0

    #1儒略世纪对应的天数
    DJC=36525.0
    
    #与参考历元之间的时间间隔(JC).
    T=((DATE1-DJ00)+DATE2)/DJC
    
    #平均倾角
    OBL=DAS2R*(84381.448+(-46.8150+(-0.00059+0.001813*T)*T)*T)
    
    return(OBL)


def pymNumat(EPSA,DPSI,DEPS):
    '''
    Form the matrix of nutation.

    Parameters
    ----------
    epsa : float 
        mean obliquity of date    
    dpsi : float
        nutation    
    deps : float
        nutation

    Returns
    -------
    rmatn : list(3,3)
        nutation matrix

    '''
    #构建旋转矩阵
    RMATN=pymIr()
    RMATN=pymRx(EPSA,RMATN)
    RMATN=pymRz(-DPSI,RMATN)
    RMATN=pymRx(-(EPSA+DEPS),RMATN)    
    
    return(RMATN)


def pymPn00(DATE1,DATE2,DPSI,DEPS):
    '''
    Precession-nutation, IAU 2000 model:  a multi-purpose function,
    supporting classical (equinox-based) use directly and CIO-based
    use indirectly.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date    
    dpsi : float
        nutation    
    deps : float
        nutation

    Returns
    -------
    epsa : flaot
        mean obliquity    
    rb : list(3,3)
        frame bias matrix    
    rp : list(3,3)
        precession matrix    
    rbp : list(3,3)
        bias-precession matrix    
    rn : list(3,3)
        nutation matrix    
    rbpn : list(3,3)
        GCRS-to-true matrix
        
    '''  
    #IAU 2000 岁差速率修正.
    DPSIPR,DEPSPR=pymPr00(DATE1,DATE2)
    
    #平均倾角，与IAU2000的岁差章动模型一致   
    EPSA=pymObl80(DATE1,DATE2)+DEPSPR
    
    #参考架偏差，岁差，以及两者的乘积
    RB,RP,RBP=pymBp00(DATE1,DATE2)
    
    #章动矩阵
    RN=pymNumat(EPSA,DPSI,DEPS)
    
    #参考架偏差-岁差-章动矩阵(经典).
    RBPN=pymRxr(RN,RBP)
    
    return(EPSA,RB,RP,RBP,RN,RBPN)


def pymPn00a(DATE1,DATE2):
    '''
    Precession-nutation, IAU 2000A model:  a multi-purpose function,
    supporting classical (equinox-based) use directly and CIO-based
    use indirectly.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date    

    Returns
    -------
    dpsi : float
        nutation    
    deps : float
        nutation    
    epsa : flaot
        mean obliquity    
    rb : list(3,3)
        frame bias matrix    
    rp : list(3,3)
        precession matrix    
    rbp : list(3,3)
        bias-precession matrix    
    rn : list(3,3)
        nutation matrix    
    rbpn : list(3,3)
        GCRS-to-true matrix

    '''
    #章动.
    DPSI,DEPS=pymNut00a(DATE1,DATE2)
    
    #调用pymPn00
    EPSA,RB,RP,RBP,RN,RBPN=pymPn00(DATE1,DATE2,DPSI,DEPS)
    
    return(DPSI,DEPS,EPSA,RB,RP,RBP,RN,RBPN)


def pymPnm00a(DATE1,DATE2):
    '''
    Form the matrix of precession-nutation for a given date (including
    frame bias), equinox based, IAU 2000A model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    rbpn : list(3,3)
        bias-precession-nutation matrix

    '''
    #调用pymPn00a获得矩阵，舍去其他参数
    DPSI,DEPS,EPSA,RB,RP,RBP,RN,RBPN=pymPn00a(DATE1, DATE2)
    
    return(RBPN)


def pymS00(DATE1,DATE2,X,Y):
    '''
    The CIO locator s, positioning the Celestial Intermediate Origin on
    the equator of the Celestial Intermediate Pole, given the CIP's X,Y
    coordinates.  Compatible with IAU 2000A precession-nutation.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date    
    x : float
        CIP coordinates    
    y : float
        CIP coordinates

    Returns
    -------
    Tfunction value : float
        the CIO locator s in radians

    '''
    #1角秒对应的弧度
    DAS2R=4.848136811095359935899141e-6

    #参考历元(J2000.0), JD
    DJ00=2451545.0

    #1儒略世纪对应的天数
    DJC=36525.0
    
    #与参考历元之间的时间间隔，占位，单位：儒略世纪数
    T=0
    FA=[0 for i in range(8)]
    
    # s+XY/2的系数
    #各系数的项数
    NSP=6
    NS0=33
    NS1=3
    NS2=25
    NS3=4
    NS4=1
    
    #多项式系数
    SP=[94e-6,3808.35e-6,-119.94e-6,-72574.09e-6,27.70e-6,15.61e-6]
    
    #参数系数  t^0
    KS0=[[0,0,0,0,1,0,0,0],
            [0,0,0,0,2,0,0,0],
            [0,0,2,-2,3,0,0,0],
            [0,0,2,-2,1,0,0,0],
            [0,0,2,-2,2,0,0,0],
            [0,0,2,0,3,0,0,0],
            [0,0,2,0,1,0,0,0],
            [0,0,0,0,3,0,0,0],
            [0,1,0,0,1,0,0,0],
            [0,1,0,0,-1,0,0,0],
            [1,0,0,0,-1,0,0,0],
            [1,0,0,0,1,0,0,0],
            [0,1,2,-2,3,0,0,0],
            [0,1,2,-2,1,0,0,0],
            [0,0,4,-4,4,0,0,0],
            [0,0,1,-1,1,-8,12,0],
            [0,0,2,0,0,0,0,0],
            [0,0,2,0,2,0,0,0],
            [1,0,2,0,3,0,0,0],
            [1,0,2,0,1,0,0,0],
            [0,0,2,-2,0,0,0,0],
            [0,1,-2,2,-3,0,0,0],
            [0,1,-2,2,-1,0,0,0],
            [0,0,0,0,0,8,-13,-1],
            [0,0,0,2,0,0,0,0],
            [2,0,-2,0,-1,0,0,0],
            [0,1,2,-2,2,0,0,0],
            [1,0,0,-2,1,0,0,0],
            [1,0,0,-2,-1,0,0,0],
            [0,0,4,-2,4,0,0,0],
            [0,0,2,-2,4,0,0,0],
            [1,0,-2,0,-3,0,0,0],
            [1,0,-2,0,-1,0,0,0]]
    
    #参数系数 t^1
    KS1=[[0,0,0,0,2,0,0,0],
             [0,0,0,0,1,0,0,0],
             [0,0,2,-2,3,0,0,0]]
    
    #参数系数 t^2
    KS2=[[0,0,0,0,1,0,0,0],
            [0,0,2,-2,2,0,0,0],
            [0,0,2,0,2,0,0,0],
            [0,0,0,0,2,0,0,0],
            [0,1,0,0,0,0,0,0],
            [1,0,0,0,0,0,0,0],
            [0,1,2,-2,2,0,0,0],
            [0,0,2,0,1,0,0,0],
            [1,0,2,0,2,0,0,0],
            [0,1,-2,2,-2,0,0,0],
            [1,0,0,-2,0,0,0,0],
            [0,0,2,-2,1,0,0,0],
            [1,0,-2,0,-2,0,0,0],
            [0,0,0,2,0,0,0,0],
            [1,0,0,0,1,0,0,0],
            [1,0,-2,-2,-2,0,0,0],
            [1,0,0,0,-1,0,0,0],
            [1,0,2,0,1,0,0,0],
            [2,0,0,-2,0,0,0,0],
            [2,0,-2,0,-1,0,0,0],
            [0,0,2,2,2,0,0,0],
            [2,0,2,0,2,0,0,0],
            [2,0,0,0,0,0,0,0],
            [1,0,2,-2,2,0,0,0],
            [0,0,2,0,0,0,0,0]]
    
    #参数系数 t^3
    KS3=[[0,0,0,0,1,0,0,0],
            [0,0,2,-2,2,0,0,0],
            [0,0,2,0,2,0,0,0],
            [0,0,0,0,2,0,0,0]]
    
    #参数系数 t^4
    KS4=[[0,0,0,0,1,0,0,0]]

    #正余弦系数 t^0
    SS0=[[-2640.73e-6,+0.39e-6],
            [-63.53e-6,+0.02e-6],
            [-11.75e-6,-0.01e-6],
            [-11.21e-6,-0.01e-6],
            [+4.57e-6,0.00e-6],
            [-2.02e-6,0.00e-6],
            [-1.98e-6,0.00e-6],
            [+1.72e-6,0.00e-6],
            [+1.41e-6,+0.01e-6],
            [+1.26e-6,+0.01e-6],
            [+0.63e-6,0.00e-6],
            [+0.63e-6,0.00e-6],
            [-0.46e-6,0.00e-6],
            [-0.45e-6,0.00e-6],
            [-0.36e-6,0.00e-6],
            [+0.24e-6,+0.12e-6],
            [-0.32e-6,0.00e-6],
            [-0.28e-6,0.00e-6],
            [-0.27e-6,0.00e-6],
            [-0.26e-6,0.00e-6],
            [+0.21e-6,0.00e-6],
            [-0.19e-6,0.00e-6],
            [-0.18e-6,0.00e-6],
            [+0.10e-6,-0.05e-6],
            [-0.15e-6,0.00e-6],
            [+0.14e-6,0.00e-6],
            [+0.14e-6,0.00e-6],
            [-0.14e-6,0.00e-6],
            [-0.14e-6,0.00e-6],
            [-0.13e-6,0.00e-6],
            [+0.11e-6,0.00e-6],
            [-0.11e-6,0.00e-6],
            [-0.11e-6,0.00e-6]]
    
    #正余弦系数 t^1
    SS1=[[-0.07e-6,+3.57e-6],
            [+1.71e-6,-0.03e-6],
            [0.00e-6,+0.48e-6]]
    
    #正余弦系数 t^2
    SS2=[[+743.53e-6,-0.17e-6],
            [+56.91e-6,+0.06e-6],
            [+9.84e-6,-0.01e-6],
            [-8.85e-6,+0.01e-6],
            [-6.38e-6,-0.05e-6],
            [-3.07e-6,0.00e-6],
            [+2.23e-6,0.00e-6],
            [+1.67e-6,0.00e-6],
            [+1.30e-6,0.00e-6],
            [+0.93e-6,0.00e-6],
            [+0.68e-6,0.00e-6],
            [-0.55e-6,0.00e-6],
            [+0.53e-6,0.00e-6],
            [-0.27e-6,0.00e-6],
            [-0.27e-6,0.00e-6],
            [-0.26e-6,0.00e-6],
            [-0.25e-6,0.00e-6],
            [+0.22e-6,0.00e-6],
            [-0.21e-6,0.00e-6],
            [+0.20e-6,0.00e-6],
            [+0.17e-6,0.00e-6],
            [+0.13e-6,0.00e-6],
            [-0.13e-6,0.00e-6],
            [-0.12e-6,0.00e-6],
            [-0.11e-6,0.00e-6]]
    
    #正余弦系数 t^3
    SS3=[[+0.30e-6,-23.51e-6],
            [-0.03e-6,-1.39e-6],
            [-0.01e-6,-0.24e-6],
            [0.00e-6,+0.22e-6]]
    
    #正余弦系数 t^4
    SS4=[[-0.26e-6,-0.01e-6]]
    
    #当前日期和参考历元之间的时间间隔，儒略世纪数.
    T=((DATE1-DJ00)+DATE2)/DJC
    
    #基本参数(from IERS Conventions 2003)
    #月亮的平近点角
    FA[0]=pymFal03(T)
    
    #太阳的的平近点角
    FA[1]=pymFalp03(T)

    #月球的平黄经减去其上升点的平黄经
    FA[2]=pymFaf03(T)

    #月亮到太阳的平均距角
    FA[3]=pymFad03(T)

    #月亮升交点的平黄经
    FA[4]=pymFaom03(T)

    #金星的平黄经
    FA[5]=pymFave03(T)

    #地球的平黄经
    FA[6]=pymFae03(T)

    #黄经的一般累计岁差
    FA[7]=pymFapa03(T)

    #估计 s.
    S0=SP[0]
    S1=SP[1]
    S2=SP[2]
    S3=SP[3]
    S4=SP[4]
    S5=SP[5]

    for i in range(NS0-1,-1,-1):
        A=0.0
        for j in range(8):
            A=A+float(KS0[i][j])*FA[j]
        S0=S0+(SS0[i][0]*ma.sin(A)+SS0[i][1]*ma.cos(A))
    
    for i in range(NS1-1,-1,-1):
        A=0.0
        for j in range(8):
            A=A+float(KS1[i][j])*FA[j]
        S1=S1+(SS1[i][0]*ma.sin(A)+SS1[i][1]*ma.cos(A))

    for i in range(NS2-1,-1,-1):
        A=0.0
        for j in range(8):
            A=A+float(KS2[i][j])*FA[j]
        S2=S2+(SS2[i][0]*ma.sin(A)+SS2[i][1]*ma.cos(A))
        
    for i in range(NS3-1,-1,-1):
        A=0.0
        for j in range(8):
            A=A+float(KS3[i][j])*FA[j]
        S3=S3+(SS3[i][0]*ma.sin(A)+SS3[i][1]*ma.cos(A))
        
    for i in range(NS4-1,-1,-1):
        A=0.0
        for j in range(8):
            A=A+float(KS4[i][j])*FA[j]
        S4=S4+(SS4[i][0]*ma.sin(A)+SS4[i][1]*ma.cos(A))

    S00=(S0+(S1+(S2+(S3+(S4+S5*T)*T)*T)*T)*T)*DAS2R-X*Y/2.0
    
    return(S00)


def pymC2ixy(DATE1,DATE2,X,Y):
    '''
    Form the celestial to intermediate-frame-of-date matrix for a given
    date when the CIP X,Y coordinates are known.  IAU 2000.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date    
    x : float
        Celestial Intermediate Pole    
    y : float
        Celestial Intermediate Pole

    Returns
    -------
    rc2i : list(3,3)
        celestial-to-intermediate matrix

    '''
    #计算s，并获得矩阵
    A=pymS00(DATE1,DATE2,X,Y)
    RC2I=pymC2ixys(X,Y,A)
    
    return(RC2I)


def pymC2ibpn(DATE1,DATE2,RBPN):
    '''
    Form the celestial-to-intermediate matrix for a given date given
    the bias-precession-nutation matrix.  IAU 2000.
    
    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date     
    date2 : float
        TT as a 2-part Julian Date    
    rbpn : list(3,3)
        celestial-to-true matrix
    
    Returns
    -------
    rc2i : list(3,3)
        celestial-to-intermediate matrix
    
    '''
    #分离出CIP的X,Y分量
    X,Y=pymBpn2xy(RBPN)
    
    #调用pymC2ixy，获得天球到中间参考系的矩阵  
    RC2I=pymC2ixy(DATE1,DATE2,X,Y)
    
    return(RC2I)


def pymC2i00a(DATE1,DATE2):
    '''
    Form the celestial-to-intermediate matrix for a given date using the
    IAU 2000A precession-nutation model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date     
    date2 : float
        TT as a 2-part Julian Date   

    Returns
    -------
    rc2i : list(3,3)
        celestial-to-intermediate matrix

    '''
    #调用pymPnm00a获得参考架偏差-岁差-章动矩阵(IAU 2000A).
    RBPN=pymPnm00a(DATE1,DATE2)
    
    #构建天球到中间参考系的矩阵
    RC2I=pymC2ibpn(DATE1,DATE2,RBPN)
    
    return(RC2I)


def pymNut00b(DATE1,DATE2):
    '''
    Nutation, IAU 2000B model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    dpsi : float
        nutation, luni-solar + planetary    
    deps : float
        nutation, luni-solar + planetary

    '''
    #1角秒对应的弧度
    DAS2R=4.848136811095359935899141e-6

    #1毫角秒对应的弧度
    DMAS2R=DAS2R/1e3
    
    #360度对应的角秒
    TURNAS=1296000.0
    
    #2Pi
    D2PI=6.283185307179586476925287
    
    #0.1微角秒对应的弧度
    U2R=DAS2R/1e7
    
    #参考历元(J2000.0), JD
    DJ00=2451545.0

    #1儒略世纪对应的天数
    DJC=36525.0
    
    #日月章动模型
    #模型中的项数
    NLS=77
    
    #基本参数的系数
    NALS=[]
    
    #经度和倾角系数
    CLS=[]
    
    #修正的偏移代替行星项，弧度制
    DPPLAN=-0.135*DMAS2R
    DEPLAN=+0.388*DMAS2R
    
    #各项参数的系数表
    
    # *  Luni-Solar argument multipliers:
    # *
    # *               L     L'    F     D     Om
    NALS=[[0,0,0,0,1],
            [0,0,2,-2,2],
            [0,0,2,0,2],
            [0,0,0,0,2],
            [0,1,0,0,0],
            [0,1,2,-2,2],
            [1,0,0,0,0],
            [0,0,2,0,1],
            [1,0,2,0,2],
            [0,-1,2,-2,2],
            [0,0,2,-2,1],
            [-1,0,2,0,2],
            [-1,0,0,2,0],
            [1,0,0,0,1],
            [-1,0,0,0,1],
            [-1,0,2,2,2],
            [1,0,2,0,1],
            [-2,0,2,0,1],
            [0,0,0,2,0],
            [0,0,2,2,2],
            [0,-2,2,-2,2],
            [-2,0,0,2,0],
            [2,0,2,0,2],
            [1,0,2,-2,2],
            [-1,0,2,0,1],
            [2,0,0,0,0],
            [0,0,2,0,0],
            [0,1,0,0,1],
            [-1,0,0,2,1],
            [0,2,2,-2,2],
            [0,0,-2,2,0],
            [1,0,0,-2,1],
            [0,-1,0,0,1],
            [-1,0,2,2,1],
            [0,2,0,0,0],
            [1,0,2,2,2],
            [-2,0,2,0,0],
            [0,1,2,0,2],
            [0,0,2,2,1],
            [0,-1,2,0,2],
            [0,0,0,2,1],
            [1,0,2,-2,1],
            [2,0,2,-2,2],
            [-2,0,0,2,1],
            [2,0,2,0,1],
            [0,-1,2,-2,1],
            [0,0,0,-2,1],
            [-1,-1,0,2,0],
            [2,0,0,-2,1],
            [1,0,0,2,0],
            [0,1,2,-2,1],
            [1,-1,0,0,0],
            [-2,0,2,0,2],
            [3,0,2,0,2],
            [0,-1,0,2,0],
            [1,-1,2,0,2],
            [0,0,0,1,0],
            [-1,-1,2,2,2],
            [-1,0,2,0,0],
            [0,-1,2,2,2],
            [-2,0,0,0,1],
            [1,1,2,0,2],
            [2,0,0,0,1],
            [-1,1,0,1,0],
            [1,1,0,0,0],
            [1,0,2,0,0],
            [-1,0,2,-2,1],
            [1,0,0,0,2],
            [-1,0,0,1,0],
            [0,0,2,1,2],
            [-1,0,2,4,2],
            [-1,1,0,1,1],
            [0,-2,2,-2,1],
            [1,0,2,2,1],
            [-2,0,2,2,2],
            [-1,0,0,0,2],
            [1,1,2,-2,2]]
    
    # *  Luni-Solar nutation coefficients, unit 1e-7 arcsec:
    # *  longitude (sin, t*sin, cos), obliquity (cos, t*cos, sin)
    CLS=[[-172064161e0,-174666e0,33386e0,92052331e0,9086e0,15377e0],
            [-13170906e0,-1675e0,-13696e0,5730336e0,-3015e0,-4587e0],
            [-2276413e0,-234e0,2796e0,978459e0,-485e0,1374e0],
            [2074554e0,207e0,-698e0,-897492e0,470e0,-291e0],
            [1475877e0,-3633e0,11817e0,73871e0,-184e0,-1924e0],
            [-516821e0,1226e0,-524e0,224386e0,-677e0,-174e0],
            [711159e0,73e0,-872e0,-6750e0,0e0,358e0],
            [-387298e0,-367e0,380e0,200728e0,18e0,318e0],
            [-301461e0,-36e0,816e0,129025e0,-63e0,367e0],
            [215829e0,-494e0,111e0,-95929e0,299e0,132e0],
            [128227e0,137e0,181e0,-68982e0,-9e0,39e0],
            [123457e0,11e0,19e0,-53311e0,32e0,-4e0],
            [156994e0,10e0,-168e0,-1235e0,0e0,82e0],
            [63110e0,63e0,27e0,-33228e0,0e0,-9e0],
            [-57976e0,-63e0,-189e0,31429e0,0e0,-75e0],
            [-59641e0,-11e0,149e0,25543e0,-11e0,66e0],
            [-51613e0,-42e0,129e0,26366e0,0e0,78e0],
            [45893e0,50e0,31e0,-24236e0,-10e0,20e0],
            [63384e0,11e0,-150e0,-1220e0,0e0,29e0],
            [-38571e0,-1e0,158e0,16452e0,-11e0,68e0],
            [32481e0,0e0,0e0,-13870e0,0e0,0e0],
            [-47722e0,0e0,-18e0,477e0,0e0,-25e0],
            [-31046e0,-1e0,131e0,13238e0,-11e0,59e0],
            [28593e0,0e0,-1e0,-12338e0,10e0,-3e0],
            [20441e0,21e0,10e0,-10758e0,0e0,-3e0],
            [29243e0,0e0,-74e0,-609e0,0e0,13e0],
            [25887e0,0e0,-66e0,-550e0,0e0,11e0],
            [-14053e0,-25e0,79e0,8551e0,-2e0,-45e0],
            [15164e0,10e0,11e0,-8001e0,0e0,-1e0],
            [-15794e0,72e0,-16e0,6850e0,-42e0,-5e0],
            [21783e0,0e0,13e0,-167e0,0e0,13e0],
            [-12873e0,-10e0,-37e0,6953e0,0e0,-14e0],
            [-12654e0,11e0,63e0,6415e0,0e0,26e0],
            [-10204e0,0e0,25e0,5222e0,0e0,15e0],
            [16707e0,-85e0,-10e0,168e0,-1e0,10e0],
            [-7691e0,0e0,44e0,3268e0,0e0,19e0],
            [-11024e0,0e0,-14e0,104e0,0e0,2e0],
            [7566e0,-21e0,-11e0,-3250e0,0e0,-5e0],
            [-6637e0,-11e0,25e0,3353e0,0e0,14e0],
            [-7141e0,21e0,8e0,3070e0,0e0,4e0],
            [-6302e0,-11e0,2e0,3272e0,0e0,4e0],
            [5800e0,10e0,2e0,-3045e0,0e0,-1e0],
            [6443e0,0e0,-7e0,-2768e0,0e0,-4e0],
            [-5774e0,-11e0,-15e0,3041e0,0e0,-5e0],
            [-5350e0,0e0,21e0,2695e0,0e0,12e0],
            [-4752e0,-11e0,-3e0,2719e0,0e0,-3e0],
            [-4940e0,-11e0,-21e0,2720e0,0e0,-9e0],
            [7350e0,0e0,-8e0,-51e0,0e0,4e0],
            [4065e0,0e0,6e0,-2206e0,0e0,1e0],
            [6579e0,0e0,-24e0,-199e0,0e0,2e0],
            [3579e0,0e0,5e0,-1900e0,0e0,1e0],
            [4725e0,0e0,-6e0,-41e0,0e0,3e0],
            [-3075e0,0e0,-2e0,1313e0,0e0,-1e0],
            [-2904e0,0e0,15e0,1233e0,0e0,7e0],
            [4348e0,0e0,-10e0,-81e0,0e0,2e0],
            [-2878e0,0e0,8e0,1232e0,0e0,4e0],
            [-4230e0,0e0,5e0,-20e0,0e0,-2e0],
            [-2819e0,0e0,7e0,1207e0,0e0,3e0],
            [-4056e0,0e0,5e0,40e0,0e0,-2e0],
            [-2647e0,0e0,11e0,1129e0,0e0,5e0],
            [-2294e0,0e0,-10e0,1266e0,0e0,-4e0],
            [2481e0,0e0,-7e0,-1062e0,0e0,-3e0],
            [2179e0,0e0,-2e0,-1129e0,0e0,-2e0],
            [3276e0,0e0,1e0,-9e0,0e0,0e0],
            [-3389e0,0e0,5e0,35e0,0e0,-2e0],
            [3339e0,0e0,-13e0,-107e0,0e0,1e0],
            [-1987e0,0e0,-6e0,1073e0,0e0,-2e0],
            [-1981e0,0e0,0e0,854e0,0e0,0e0],
            [4026e0,0e0,-353e0,-553e0,0e0,-139e0],
            [1660e0,0e0,-5e0,-710e0,0e0,-2e0],
            [-1521e0,0e0,9e0,647e0,0e0,4e0],
            [1314e0,0e0,0e0,-700e0,0e0,0e0],
            [-1283e0,0e0,0e0,672e0,0e0,0e0],
            [-1331e0,0e0,8e0,663e0,0e0,4e0],
            [1383e0,0e0,-2e0,-594e0,0e0,-2e0],
            [1405e0,0e0,4e0,-610e0,0e0,2e0],
            [1290e0,0e0,0e0,-556e0,0e0,0e0]]
    
    #给定日期与参考历元之间的时间间隔，儒略世纪数.
    T=((DATE1-DJ00)+DATE2)/DJC
    
    #日月章动
    
    #Fundamental (Delaunay) arguments from Simon et al. (1994)
    
    #月亮的平近点角
    A=(485868.249036+(+1717915923.2178)*T)%TURNAS
    EL=A*DAS2R

    #太阳的平近点角
    A=(1287104.79305+(+129596581.0481)*T)%TURNAS
    ELP=A*DAS2R

    #月亮的平均升交角距
    A=(335779.526232+(+1739527262.8478)*T)%TURNAS
    F=A*DAS2R

    #月亮到太阳的平均角距
    A=(1072260.70369+(+1602961601.2090)*T)%TURNAS
    D=A*DAS2R

    #月亮的升交点平黄经
    A=(450160.398036+(-6962890.5431)*T)%TURNAS
    OM=A*DAS2R
    
    #初始化章动量
    DP=0.0
    DE=0.0
    
    #日月章动级数求和(倒序)
    for i in range(NLS-1,-1,-1):
        B=float(NALS[i][0])*EL+float(NALS[i][1])*ELP+\
           float(NALS[i][2])*F+float(NALS[i][3])*D+\
           float(NALS[i][4])*OM
        ARG=B%D2PI
        SARG=ma.sin(ARG)
        CARG=ma.cos(ARG)
    
        #Term.
        DP=DP+(CLS[i][0]+CLS[i][1]*T)*SARG+CLS[i][2]*CARG
        DE=DE+(CLS[i][3]+CLS[i][4]*T)*CARG+CLS[i][5]*SARG
    
    #从0.1微弧秒单位转换为弧度。
    DPSILS=DP*U2R
    DEPSLS=DE*U2R
    
    #代替行星章动项
    
    #修正偏移量，以纠正截断系数中缺失的项
    DPSIPL=DPPLAN
    DEPSPL=DEPLAN
    
    #日月和行星章动相加
    DPSI=DPSILS+DPSIPL
    DEPS=DEPSLS+DEPSPL
    
    return(DPSI,DEPS)


def pymPn00b(DATE1,DATE2):
    '''
    Precession-nutation, IAU 2000B model:  a multi-purpose function,
    supporting classical (equinox-based) use directly and CIO-based
    use indirectly.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date    

    Returns
    -------
    dpsi : float
        nutation    
    deps : float
        nutation    
    epsa : flaot
        mean obliquity    
    rb : list(3,3)
        frame bias matrix    
    rp : list(3,3)
        precession matrix    
    rbp : list(3,3)
        bias-precession matrix    
    rn : list(3,3)
        nutation matrix    
    rbpn : list(3,3)
        GCRS-to-true matrix

    '''  
    
    #章动，调用pymNut00b
    DPSI,DEPS=pymNut00b(DATE1,DATE2)
    
    #Remaining results.
    EPSA,RB,RP,RBP,RN,RBPN=pymPn00(DATE1,DATE2,DPSI,DEPS)
    
    return(DPSI,DEPS,EPSA,RB,RP,RBP,RN,RBPN)


def pymPnm00b(DATE1,DATE2):
    '''
    Form the matrix of precession-nutation for a given date (including
    frame bias), equinox-based, IAU 2000B model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    rbpn : list(3,3)
        bias-precession-nutation matrix

    '''
    #获得需要的矩阵，舍去其他无关参数
    DPSI,DEPS,EPSA,RB,RP,RBP,RN,RBPN=pymPn00b(DATE1,DATE2)
    
    return(RBPN)


def pymC2i00b(DATE1,DATE2):
    '''
    Form the celestial-to-intermediate matrix for a given date using the
    IAU 2000B precession-nutation model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date     
    date2 : float
        TT as a 2-part Julian Date   

    Returns
    -------
    rc2i : list(3,3)
        celestial-to-intermediate matrix

    '''
    #获得天球到真赤道和春分点的矩阵 (IAU 2000B).
    RBPN=pymPnm00b(DATE1,DATE2)
    
    #构建天球到中间参考系的矩阵
    RC2I=pymC2ibpn(DATE1,DATE2,RBPN)
    
    return(RC2I)


def pymC2i06a(DATE1,DATE2):
    '''
    Form the celestial-to-intermediate matrix for a given date using the
    IAU 2006 precession and IAU 2000A nutation models.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    rc2i : list(3,3)
        celestial-to-intermediate matrix

    '''
    #获得天球到真赤道和春分点的矩阵(IAU 2006/2000A).
    RBPN=pymPnm06a(DATE1,DATE2)
    
    #分离出CIP的X,Y分量
    X,Y=pymBpn2xy(RBPN)
    
    #获得CIO定位角s
    S=pymS06(DATE1,DATE2,X,Y)

    #构建天球到中间参考系的矩阵
    RC2I=pymC2ixys(X,Y,S)
    
    return(RC2I)


def pymC2tcio(RC2I,ERA,RPOM):
    '''
    Assemble the celestial to terrestrial matrix from CIO-based
    components (the celestial-to-intermediate matrix, the Earth Rotation
    Angle and the polar motion matrix).

    Parameters
    ----------
    rc2i : list(3,3)
        celestial-to-intermediate matrix     
    era : float
        Earth rotation angle (radians)     
    rpom : list(3,3)
        polar-motion matrix

    Returns
    -------
    rc2t : list(3,3)
        celestial-to-terrestrial matrix

    '''    
    #构建矩阵
    R=pymCr(RC2I)
    R=pymRz(ERA,R)
    RC2T=pymRxr(RPOM,R)
    
    return(RC2T)


def pymC2t00a(TTA,TTB,UTA,UTB,XP,YP):
    '''
    Form the celestial to terrestrial matrix given the date, the UT1 and
    the polar motion, using the IAU 2000A precession-nutation model.

    Parameters
    ----------
    tta : float
        TT as a 2-part Julian Date    
    ttb : float
        TT as a 2-part Julian Date    
    uta : float
        UT1 as a 2-part Julian Date    
    utb : float
        UT1 as a 2-part Julian Date    
    xp : float
        coordinates of the pole    
    yp : float
        coordinates of the pole  

    Returns
    -------
    rc2t : list(3,3)
        celestial-to-terrestrial matrix

    '''
    #构建在TT时刻的天球到中间参考系的矩阵(IAU 2000A).
    RC2I=pymC2i00a(TTA,TTB)
    
    #调用pymEra00预测在UT1时刻的地球自转角
    ERA=pymEra00(UTA,UTB)
    
    #预测TIO定位角 s'.
    SP=pymSp00(TTA,TTB)
    
    #构建极移矩阵
    RPOM=pymPom00(XP,YP,SP)
    
    #组合构建天球到地球的矩阵
    RC2T=pymC2tcio(RC2I,ERA,RPOM)
    
    return(RC2T)


def pymC2t00b(TTA,TTB,UTA,UTB,XP,YP):
    '''
    Form the celestial to terrestrial matrix given the date, the UT1 and
    the polar motion, using the IAU 2000B precession-nutation model.

    Parameters
    ----------
    tta : float
        TT as a 2-part Julian Date    
    ttb : float
        TT as a 2-part Julian Date    
    uta : float
        UT1 as a 2-part Julian Date    
    utb : float
        UT1 as a 2-part Julian Date    
    xp : float
        coordinates of the pole    
    yp : float
        coordinates of the pole  

    Returns
    -------
    rc2t : list(3,3)
        celestial-to-terrestrial matrix
        
    '''
    #构建在TT时刻的天球到中间参考系的矩阵(IAU 2000B).
    RC2I=pymC2i00b(TTA,TTB)
    
    #调用pymEra00预测在UT1时刻的地球自转角
    ERA=pymEra00(UTA,UTB)
    
    #构建极移矩阵，忽略TIO定位角s'
    RPOM=pymPom00(XP,YP,0.0)
    
    #组合构建天球到地球的矩阵
    RC2T=pymC2tcio(RC2I,ERA,RPOM)
    
    return(RC2T)


def pymC2t06a(TTA,TTB,UTA,UTB,XP,YP):
    '''
    Form the celestial to terrestrial matrix given the date, the UT1 and
    the polar motion, using the IAU 2006/2000A precession-nutation
    model.

    Parameters
    ----------
    tta : float
        TT as a 2-part Julian Date     
    ttb : float
        TT as a 2-part Julian Date    
    uta : float
        UT1 as a 2-part Julian Date    
    utb : float
        UT1 as a 2-part Julian Date     
    xp : float
        coordinates of the pole (radians)    
    yp : float
        coordinates of the pole (radians)

    Returns
    -------
    rc2t : list(3,3)
        celestial-to-terrestrial matrix

    '''
    #构建在TT时刻的天球到中间参考系的矩阵(IAU 2006/2000A).
    RC2I=pymC2i06a(TTA,TTB)
    
    #调用pymEra00预测在UT1时刻的地球自转角
    ERA=pymEra00(UTA,UTB)
    
    #预测TIO定位角 s'.
    SP=pymSp00(TTA,TTB)
    
    #构建极移矩阵
    RPOM=pymPom00(XP,YP,SP)
    
    #组合构建天球到地球的矩阵
    RC2T=pymC2tcio(RC2I,ERA,RPOM)
    
    return(RC2T)


def pymC2teqx(RBPN,GST,RPOM):
    '''
    Assemble the celestial to terrestrial matrix from equinox-based
    components (the celestial-to-true matrix, the Greenwich Apparent
    Sidereal Time and the polar motion matrix).

    Parameters
    ----------
    rbpn : list(3,3)
        celestial-to-true matrix
    gst : float
        Greenwich (apparent) Sidereal Time (radians)
    rpom : list(3,3)
        polar-motion matrix

    Returns
    -------
    rc2t : list(3,3)
        celestial-to-terrestrial matrix

    '''
    #构建矩阵
    R=pymCr(RBPN)
    R=pymRz(GST,R)
    RC2T=pymRxr(RPOM,R)
    
    return(RC2T)


def pymGmst00(UTA,UTB,TTA,TTB):
    '''
    Greenwich mean sidereal time (model consistent with IAU 2000
    resolutions).

    Parameters
    ----------
    uta : float
        UT1 as a 2-part Julian Date    
    utb : float
        UT1 as a 2-part Julian Date    
    tta : float
        TT as a 2-part Julian Date    
    ttb : float
        TT as a 2-part Julian Date

    Returns
    -------
    function value : float
        Greenwich mean sidereal time (radians)

    '''
    #1角秒对应的弧度
    DAS2R=4.848136811095359935899141e-6

    #参考历元(J2000.0), JD
    DJ00=2451545.0

    #1儒略世纪对应的天数
    DJC=36525.0
    
    #自参考历元起的世界时间隔
    T=((TTA-DJ00)+TTB)/DJC
    
    #格林尼治平恒星时, IAU 2000.
    A=pymEra00(UTA,UTB)
    B=A+(0.014506+(4612.15739966+(+1.39667721+(-0.00009344+(+0.00001882)
                                               *T)*T)*T)*T)*DAS2R
    GMST=pymAnp(B)  
    
    return(GMST)


def pymEect00(DATE1,DATE2):
    '''
    Equation of the equinoxes complementary terms, consistent with
    IAU 2000 resolutions.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    function value : float
        complementary terms

    '''
    #1角秒对应的弧度
    DAS2R=4.848136811095359935899141e-6

    #参考历元(J2000.0), JD
    DJ00=2451545.0

    #1儒略世纪的天数
    DJC=36525.0
    
    #基本参数
    FA=[0 for i in range(14)]
    
    # EE补充项部分
    
    #项数
    NE0=33
    NE1=1
    
    #Argument coefficients for t^0
    KE0=[[0,0,0,0,1,0,0,0],
            [0,0,0,0,2,0,0,0],
            [0,0,2,-2,3,0,0,0],
            [0,0,2,-2,1,0,0,0],
            [0,0,2,-2,2,0,0,0],
            [0,0,2,0,3,0,0,0],
            [0,0,2,0,1,0,0,0],
            [0,0,0,0,3,0,0,0],
            [0,1,0,0,1,0,0,0],
            [0,1,0,0,-1,0,0,0],
            [1,0,0,0,-1,0,0,0],
            [1,0,0,0,1,0,0,0],
            [0,1,2,-2,3,0,0,0],
            [0,1,2,-2,1,0,0,0],
            [0,0,4,-4,4,0,0,0],
            [0,0,1,-1,1,-8,12,0],
            [0,0,2,0,0,0,0,0],
            [0,0,2,0,2,0,0,0],
            [1,0,2,0,3,0,0,0],
            [1,0,2,0,1,0,0,0],
            [0,0,2,-2,0,0,0,0],
            [0,1,-2,2,-3,0,0,0],
            [0,1,-2,2,-1,0,0,0],
            [0,0,0,0,0,8,-13,-1],
            [0,0,0,2,0,0,0,0],
            [2,0,-2,0,-1,0,0,0],
            [1,0,0,-2,1,0,0,0],
            [0,1,2,-2,2,0,0,0],
            [1,0,0,-2,-1,0,0,0],
            [0,0,4,-2,4,0,0,0],
            [0,0,2,-2,4,0,0,0],
            [1,0,-2,0,-3,0,0,0],
            [1,0,-2,0,-1,0,0,0]]
    
    #Argument coefficients for t^1
    KE1=[[0,0,0,0,1,0,0,0]]
    
    #Sine and cosine coefficients for t^0
    SE0=[[+2640.96e-6,-0.39e-6],
            [+63.52e-6,-0.02e-6],
            [+11.75e-6,+0.01e-6],
            [+11.21e-6,+0.01e-6],
            [-4.55e-6,+0.00e-6],
            [+2.02e-6,+0.00e-6],
            [+1.98e-6,+0.00e-6],
            [-1.72e-6,+0.00e-6],
            [-1.41e-6,-0.01e-6],
            [-1.26e-6,-0.01e-6],
            [-0.63e-6,+0.00e-6],
            [-0.63e-6,+0.00e-6],
            [+0.46e-6,+0.00e-6],
            [+0.45e-6,+0.00e-6],
            [+0.36e-6,+0.00e-6],
            [-0.24e-6,-0.12e-6],
            [+0.32e-6,+0.00e-6],
            [+0.28e-6,+0.00e-6],
            [+0.27e-6,+0.00e-6],
            [+0.26e-6,+0.00e-6],
            [-0.21e-6,+0.00e-6],
            [+0.19e-6,+0.00e-6],
            [+0.18e-6,+0.00e-6],
            [-0.10e-6,+0.05e-6],
            [+0.15e-6,+0.00e-6],
            [-0.14e-6,+0.00e-6],
            [+0.14e-6,+0.00e-6],
            [-0.14e-6,+0.00e-6],
            [+0.14e-6,+0.00e-6],
            [+0.13e-6,+0.00e-6],
            [-0.11e-6,+0.00e-6],
            [+0.11e-6,+0.00e-6],
            [+0.11e-6,+0.00e-6]]
    
    #Sine and cosine coefficients for t^1
    SE1=[[-0.87e-6,+0.00e-6]]
    
    #给定日期和参考历元之间的时间间隔儒略世纪数.
    T=((DATE1-DJ00)+DATE2)/DJC
    
    #基本参数(依据IERS Conventions 2003)
    #月亮的平近点角
    FA[0]=pymFal03(T)
    
    #太阳的平近点角
    FA[1]=pymFalp03(T)

    #月亮的平黄经减去升交点黄经
    FA[2]=pymFaf03(T)

    #月亮到太阳的平均距角
    FA[3]=pymFad03(T)

    #月亮的升交点平黄经
    FA[4]=pymFaom03(T)

    #金星的平黄经
    FA[5]=pymFave03(T)

    #地球的平黄经
    FA[6]=pymFae03(T)

    #黄经上的累计岁差
    FA[7]=pymFapa03(T)
    
    #估计EE补充项
    S0=0.0
    S1=0.0
    
    for i in range(NE0-1,-1,-1):
        A=0.0
        for j in range(8):
            A=A+float(KE0[i][j])*FA[j]
        S0=S0+(SE0[i][0]*ma.sin(A)+SE0[i][1]*ma.cos(A))
    
    for i in range(NE1-1,-1,-1):
        A=0.0
        for j in range(8):
            A=A+float(KE1[i][j])*FA[j]
        S1=S1+(SE1[i][0]*ma.sin(A)+SE1[i][1]*ma.cos(A))
    
    EECT=(S0+S1*T)*DAS2R    
    
    return(EECT)


def pymEe00(DATE1,DATE2,EPSA,DPSI):
    '''
    The equation of the equinoxes, compatible with IAU 2000 resolutions,
    given the nutation in longitude and the mean obliquity.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date    
    epsa : float
        mean obliquity    
    dpsi : float
        nutation in longitude

    Returns
    -------
    function value : float
        equation of the equinoxes

    '''
    #赤经章动
    EE=DPSI*ma.cos(EPSA)+pymEect00(DATE1,DATE2)
    
    return(EE)


def pymC2tpe(TTA,TTB,UTA,UTB,DPSI,DEPS,XP,YP):
    '''
    Assemble the celestial to terrestrial matrix from equinox-based
    components (the celestial-to-true matrix, the Greenwich Apparent
    Sidereal Time and the polar motion matrix).

    Parameters
    ----------
    rbpn : list(3,3)
        celestial-to-true matrix
    gst : float
        Greenwich (apparent) Sidereal Time (radians)
    rpom : list(3,3)
        polar-motion matrix

    Returns
    -------
    rc2t : list(3,3)
        celestial-to-terrestrial matrix

    '''
    #构建给定时间的从天球到真赤道的矩阵
    EPSA,RB,RP,RBP,RN,RBPN=pymPn00(TTA,TTB,DPSI,DEPS)
    
    #预测给定时间的格林尼治平恒星时.
    GMST=pymGmst00(UTA,UTB,TTA,TTB)
    
    #基于给定的时间和章动量，预测赤经章动的补充项.
    EE=pymEe00(TTA,TTB,EPSA,DPSI)
    
    #估计TIO定位角s'.
    SP=pymSp00(TTA,TTB)
    
    #构建极移矩阵
    RPOM=pymPom00(XP,YP,SP)
    
    #组合构建天球到地球的矩阵
    RC2T=pymC2teqx(RBPN,GMST+EE,RPOM)
    
    return(RC2T)


def pymC2txy(TTA,TTB,UTA,UTB,X,Y,XP,YP):
    '''
    Form the celestial to terrestrial matrix given the date, the UT1,
    the CIP coordinates and the polar motion.  IAU 2000.

    Parameters
    ----------
    tta : float
        TT as a 2-part Julian Date    
    ttb : float
        TT as a 2-part Julian Date    
    uta : float
        UT1 as a 2-part Julian Date    
    utb : float
        UT1 as a 2-part Julian Date    
    x : float
        Celestial Intermediate Pole    
    y : float
        Celestial Intermediate Pole    
    xp : float
        coordinates of the pole    
    yp : float
        coordinates of the pole    

    Returns
    -------
    rc2t : list(3,3)
        celestial-to-terrestrial matrix

    '''
    #构建给定时间的从天球到中间参考架的矩阵.
    RC2I=pymC2ixy(TTA,TTB,X,Y)
    
    #预测在UT1时刻的地球自转角
    ERA=pymEra00(UTA,UTB)
    
    #估计 s'.
    SP=pymSp00(TTA,TTB)
 
    #构建极移矩阵
    RPOM=pymPom00(XP,YP,SP)
    
    #组合构建天球到地球的矩阵
    RC2T=pymC2tcio(RC2I,ERA,RPOM)
    
    return(RC2T)


def pymCpv(PV):
    '''
    Copy a position/velocity vector.

    Parameters
    ----------
    pv : list(2,3)
        position/velocity vector to be copied

    Returns
    -------
    c : list(2,3)
        copy

    '''
    C=[0,0]
    
    C[0]=pymCp(PV[0])
    C[1]=pymCp(PV[1])
    
    return(C)


def pymEcm06(DATE1,DATE2):
    '''
    ICRS equatorial to ecliptic rotation matrix, IAU 2006.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date     

    Returns
    -------
    rm : list(3,3)
        ICRS to ecliptic rotation matrix

    '''
    #倾角, IAU 2006.
    OB=pymObl06(DATE1,DATE2)
    
    #参考架偏差-岁差矩阵, IAU 2006.
    BP=pymPmat06(DATE1,DATE2)
    
    #给定日期赤道到黄道的旋转矩阵
    E=pymIr()
    E=pymRx(OB,E)
    
    #ICRS到黄道的旋转矩阵
    RM=pymRxr(E,BP)
    
    return(RM)


def pymEceq06(DATE1,DATE2,DL,DB):
    '''
    Transformation from ecliptic coordinates (mean equinox and ecliptic
    of date) to ICRS RA,Dec, using the IAU 2006 precession model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date     
    dl : float
        ecliptic longitude (radians)
    db : float
        ecliptic latitude (radians)

    Returns
    -------
    dr : float
        ICRS right ascension (radians)    
    dd : float
        ICRS declination (radians)

    '''
    #球面到笛卡尔坐标
    V1=pymS2c(DL,DB)
    
    #调用pymEcm06获得从ICRS赤道到黄道的旋转矩阵
    RM=pymEcm06(DATE1,DATE2)
    
    #通过转置旋转矩阵，从黄道转换到ICRS赤道
    V2=pymTrxp(RM,V1)
    
    #笛卡尔系到球面
    A,B=pymC2s(V2)
    
    #规范化表示
    DR=pymAnp(A)
    DD=pymAnpm(B)    
    
    return(DR,DD)


def pymEe00a(DATE1,DATE2):
    '''
    Equation of the equinoxes, compatible with IAU 2000 resolutions.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    function value : float
        equation of the equinoxes

    '''
    #IAU 2000岁差速率改正.
    DPSIPR,DEPSPR=pymPr00(DATE1,DATE2)
    
    #平均倾角，与IAU 2000岁差章动模型一致.
    EPSA=pymObl80(DATE1,DATE2)+DEPSPR
    
    #经度上的章动
    DPSI,DEPS=pymNut00a(DATE1,DATE2)
    
    #赤经章动
    EE=pymEe00(DATE1,DATE2,EPSA,DPSI)
    
    return(EE)


def pymEe00b(DATE1,DATE2):
    '''
    Equation of the equinoxes, compatible with IAU 2000 resolutions but
    using the truncated nutation model IAU 2000B.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    function value : float
        equation of the equinoxes

    '''
    #IAU 2000岁差速率改正.
    DPSIPR,DEPSPR=pymPr00(DATE1,DATE2)
    
    #平均倾角，与IAU 2000岁差章动模型一致.
    EPSA=pymObl80(DATE1,DATE2)+DEPSPR
    
    #经度上的章动，采用删减后的章动模型IAU 2000B  
    DPSI,DEPS=pymNut00b(DATE1,DATE2)
    
    #赤经章动
    EE=pymEe00(DATE1,DATE2,EPSA,DPSI)
    
    return(EE)


def pymGmst06(UTA,UTB,TTA,TTB):
    '''
    Greenwich mean sidereal time (consistent with IAU 2006 precession).

    Parameters
    ----------
    uta : float
        UT1 as a 2-part Julian Date    
    utb : float
        UT1 as a 2-part Julian Date    
    tta : float
        TT as a 2-part Julian Date    
    ttb : float
        TT as a 2-part Julian Date

    Returns
    -------
    function value : float
        Greenwich mean sidereal time (radians)

    '''
    #1角秒对应的弧度
    DAS2R=4.848136811095359935899141e-6

    #参考历元(J2000.0), JD
    DJ00=2451545.0

    #1儒略世纪对应的天数
    DJC=36525.0
    
    #自参考历元起的世界时间隔
    T=((TTA-DJ00)+TTB)/DJC
    
    #格林尼治平恒星时, IAU 2006.
    A=pymEra00(UTA,UTB)
    GMST=pymAnp(A+(0.014506+(4612.156534+(1.3915817+(-0.00000044+
                    (-0.000029956+(-0.0000000368)*T)*T)*T)*T)*T)*DAS2R)
    
    return(GMST)


def pymGst06(UTA,UTB,TTA,TTB,RNPB):
    '''
    Greenwich apparent sidereal time, IAU 2006, given the NPB matrix.

    Parameters
    ----------
    uta : float
        UT1 as a 2-part Julian Date    
    utb : float
        UT1 as a 2-part Julian Date    
    tta : float
        TT as a 2-part Julian Date     
    ttb : float
        TT as a 2-part Julian Date    
    rnpb : float
        nutation x precession x bias matrix

    Returns
    -------
    function value : float
        Greenwich apparent sidereal time (radians)

    '''
    #从CIP坐标中分离出X,Y
    X,Y=pymBpn2xy(RNPB)
    
    #CIO定位角, s.
    S=pymS06(TTA,TTB,X,Y)
    
    #格林尼治视恒星时
    GST=pymAnp(pymEra00(UTA,UTB)-pymEors(RNPB,S))
    
    return(GST)


def pymGst06a(UTA,UTB,TTA,TTB):
    '''
    Greenwich apparent sidereal time (consistent with IAU 2000 and 2006
    resolutions).

    Parameters
    ----------
    uta : float
        UT1 as a 2-part Julian Date    
    utb : float
        UT1 as a 2-part Julian Date     
    tta : float
        TT as a 2-part Julian Date    
    ttb : float
        TT as a 2-part Julian Date

    Returns
    -------
    function value : float
        Greenwich apparent sidereal time (radians)

    '''
    #调用pymPnm06a获得NPB（参考架偏差-岁差-章动）矩阵, IAU 2000A/2006.
    RNPB=pymPnm06a(TTA,TTB)
    
    #格林尼治视恒星时
    GST=pymGst06(UTA,UTB,TTA,TTB,RNPB)
    
    return(GST)


def pymEe06a(DATE1,DATE2):
    '''
    Equation of the equinoxes, compatible with IAU 2000 resolutions and
    IAU 2006/2000A precession-nutation.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    function value : float
        equation of the equinoxes

    '''
    #赤经章动
    A=pymGst06a(0.0,0.0,DATE1,DATE2)-pymGmst06(0.0,0.0,DATE1,DATE2)
    EE=pymAnpm(A)
    
    return(EE)


def pymEo06a (DATE1,DATE2):
    '''
    Equation of the origins, IAU 2006 precession and IAU 2000A nutation.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    function value : float
        the equation of the origins in radians

    '''
    #NPB（参考架偏差-岁差-章动）矩阵
    R=pymPnm06a(DATE1,DATE2)
    
    #从CIP中分离出X,Y
    X,Y=pymBpn2xy(R)
    
    #CIO定位角, s.
    S=pymS06(DATE1,DATE2,X,Y)
    
    #求解EO
    EO=pymEors(R,S)
    
    return(EO)


def pymEqec06(DATE1,DATE2,DR,DD):
    '''
    Transformation from ICRS equatorial coordinates to ecliptic
    coordinates (mean equinox and ecliptic of date) using IAU 2006
    precession model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian date    
    date2 : float
        TT as a 2-part Julian date    
    dr : float
        ICRS right ascension (radians)    
    dd : float
        ICRS declination (radians)

    Returns
    -------
    dl : flaot
        ecliptic longitude (radians)    
    db : float
        ecliptic latitude (radians)

    '''
    #球面到笛卡尔坐标
    V1=pymS2c(DR,DD)
    
    #旋转矩阵，从ICRS 赤道到黄道.
    RM=pymEcm06(DATE1,DATE2)
    
    #将坐标从ICRS变换到黄道.
    V2=pymRxp(RM,V1)
    
    #笛卡尔坐标到球面
    A,B=pymC2s(V2)
    
    #规范化表示
    DL=pymAnp(A)
    DB=pymAnpm(B)
     
    return(DL,DB)


def pymNut80(DATE1,DATE2):
    '''
    Nutation, IAU 1980 model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    dpsi : float
        nutation in longitude (radians)    
    deps : float
        nutation in obliquity (radians)

    '''
    #1角秒对应的弧度
    DAS2R=4.848136811095359935899141e-6

    #2Pi
    D2PI=6.283185307179586476925287
    
    #0.1微角秒对应的弧度
    U2R=DAS2R/1e4
    
    #参考历元(J2000.0), JD
    DJ00=2451545.0

    #1儒略世纪对应的天数
    DJC=36525.0
    
    # 参数和多项式系数表
    #系数单位为0.1微角秒，变化速率为角秒/千儒略年
    # *                Multiple of            Longitude        Obliquity
    # *           L    L'   F    D  Omega   coeff. of sin    coeff. of cos
    # *                                         1       t        1     t
    
    X=[[0.,0.,0.,0.,1.,-171996.,-1742.,92025.,89.],
            [0.,0.,0.,0.,2.,2062.,2.,-895.,5.],
            [-2.,0.,2.,0.,1.,46.,0.,-24.,0.],
            [2.,0.,-2.,0.,0.,11.,0.,0.,0.],
            [-2.,0.,2.,0.,2.,-3.,0.,1.,0.],
            [1.,-1.,0.,-1.,0.,-3.,0.,0.,0.],
            [0.,-2.,2.,-2.,1.,-2.,0.,1.,0.],
            [2.,0.,-2.,0.,1.,1.,0.,0.,0.],
            [0.,0.,2.,-2.,2.,-13187.,-16.,5736.,-31.],
            [0.,1.,0.,0.,0.,1426.,-34.,54.,-1.],
            [0.,1.,2.,-2.,2.,-517.,12.,224.,-6.],
            [0.,-1.,2.,-2.,2.,217.,-5.,-95.,3.],
            [0.,0.,2.,-2.,1.,129.,1.,-70.,0.],
            [2.,0.,0.,-2.,0.,48.,0.,1.,0.],
            [0.,0.,2.,-2.,0.,-22.,0.,0.,0.],
            [0.,2.,0.,0.,0.,17.,-1.,0.,0.],
            [0.,1.,0.,0.,1.,-15.,0.,9.,0.],
            [0.,2.,2.,-2.,2.,-16.,1.,7.,0.],
            [0.,-1.,0.,0.,1.,-12.,0.,6.,0.],
            [-2.,0.,0.,2.,1.,-6.,0.,3.,0.],
            [0.,-1.,2.,-2.,1.,-5.,0.,3.,0.],
            [2.,0.,0.,-2.,1.,4.,0.,-2.,0.],
            [0.,1.,2.,-2.,1.,4.,0.,-2.,0.],
            [1.,0.,0.,-1.,0.,-4.,0.,0.,0.],
            [2.,1.,0.,-2.,0.,1.,0.,0.,0.],
            [0.,0.,-2.,2.,1.,1.,0.,0.,0.],
            [0.,1.,-2.,2.,0.,-1.,0.,0.,0.],
            [0.,1.,0.,0.,2.,1.,0.,0.,0.],
            [-1.,0.,0.,1.,1.,1.,0.,0.,0.],
            [0.,1.,2.,-2.,0.,-1.,0.,0.,0.],
            [0.,0.,2.,0.,2.,-2274.,-2.,977.,-5.],
            [1.,0.,0.,0.,0.,712.,1.,-7.,0.],
            [0.,0.,2.,0.,1.,-386.,-4.,200.,0.],
            [1.,0.,2.,0.,2.,-301.,0.,129.,-1.],
            [1.,0.,0.,-2.,0.,-158.,0.,-1.,0.],
            [-1.,0.,2.,0.,2.,123.,0.,-53.,0.],
            [0.,0.,0.,2.,0.,63.,0.,-2.,0.],
            [1.,0.,0.,0.,1.,63.,1.,-33.,0.],
            [-1.,0.,0.,0.,1.,-58.,-1.,32.,0.],
            [-1.,0.,2.,2.,2.,-59.,0.,26.,0.],
            [1.,0.,2.,0.,1.,-51.,0.,27.,0.],
            [0.,0.,2.,2.,2.,-38.,0.,16.,0.],
            [2.,0.,0.,0.,0.,29.,0.,-1.,0.],
            [1.,0.,2.,-2.,2.,29.,0.,-12.,0.],
            [2.,0.,2.,0.,2.,-31.,0.,13.,0.],
            [0.,0.,2.,0.,0.,26.,0.,-1.,0.],
            [-1.,0.,2.,0.,1.,21.,0.,-10.,0.],
            [-1.,0.,0.,2.,1.,16.,0.,-8.,0.],
            [1.,0.,0.,-2.,1.,-13.,0.,7.,0.],
            [-1.,0.,2.,2.,1.,-10.,0.,5.,0.],
            [1.,1.,0.,-2.,0.,-7.,0.,0.,0.],
            [0.,1.,2.,0.,2.,7.,0.,-3.,0.],
            [0.,-1.,2.,0.,2.,-7.,0.,3.,0.],
            [1.,0.,2.,2.,2.,-8.,0.,3.,0.],
            [1.,0.,0.,2.,0.,6.,0.,0.,0.],
            [2.,0.,2.,-2.,2.,6.,0.,-3.,0.],
            [0.,0.,0.,2.,1.,-6.,0.,3.,0.],
            [0.,0.,2.,2.,1.,-7.,0.,3.,0.],
            [1.,0.,2.,-2.,1.,6.,0.,-3.,0.],
            [0.,0.,0.,-2.,1.,-5.,0.,3.,0.],
            [1.,-1.,0.,0.,0.,5.,0.,0.,0.],
            [2.,0.,2.,0.,1.,-5.,0.,3.,0.],
            [0.,1.,0.,-2.,0.,-4.,0.,0.,0.],
            [1.,0.,-2.,0.,0.,4.,0.,0.,0.],
            [0.,0.,0.,1.,0.,-4.,0.,0.,0.],
            [1.,1.,0.,0.,0.,-3.,0.,0.,0.],
            [1.,0.,2.,0.,0.,3.,0.,0.,0.],
            [1.,-1.,2.,0.,2.,-3.,0.,1.,0.],
            [-1.,-1.,2.,2.,2.,-3.,0.,1.,0.],
            [-2.,0.,0.,0.,1.,-2.,0.,1.,0.],
            [3.,0.,2.,0.,2.,-3.,0.,1.,0.],
            [0.,-1.,2.,2.,2.,-3.,0.,1.,0.],
            [1.,1.,2.,0.,2.,2.,0.,-1.,0.],
            [-1.,0.,2.,-2.,1.,-2.,0.,1.,0.],
            [2.,0.,0.,0.,1.,2.,0.,-1.,0.],
            [1.,0.,0.,0.,2.,-2.,0.,1.,0.],
            [3.,0.,0.,0.,0.,2.,0.,0.,0.],
            [0.,0.,2.,1.,2.,2.,0.,-1.,0.],
            [-1.,0.,0.,0.,2.,1.,0.,-1.,0.],
            [1.,0.,0.,-4.,0.,-1.,0.,0.,0.],
            [-2.,0.,2.,2.,2.,1.,0.,-1.,0.],
            [-1.,0.,2.,4.,2.,-2.,0.,1.,0.],
            [2.,0.,0.,-4.,0.,-1.,0.,0.,0.],
            [1.,1.,2.,-2.,2.,1.,0.,-1.,0.],
            [1.,0.,2.,2.,1.,-1.,0.,1.,0.],
            [-2.,0.,2.,4.,2.,-1.,0.,1.,0.],
            [-1.,0.,4.,0.,2.,1.,0.,0.,0.],
            [1.,-1.,0.,-2.,0.,1.,0.,0.,0.],
            [2.,0.,2.,-2.,1.,1.,0.,-1.,0.],
            [2.,0.,2.,2.,2.,-1.,0.,0.,0.],
            [1.,0.,0.,2.,1.,-1.,0.,0.,0.],
            [0.,0.,4.,-2.,2.,1.,0.,0.,0.],
            [3.,0.,2.,-2.,2.,1.,0.,0.,0.],
            [1.,0.,2.,-2.,0.,-1.,0.,0.,0.],
            [0.,1.,2.,0.,1.,1.,0.,0.,0.],
            [-1.,-1.,0.,2.,1.,1.,0.,0.,0.],
            [0.,0.,-2.,0.,1.,-1.,0.,0.,0.],
            [0.,0.,2.,-1.,2.,-1.,0.,0.,0.],
            [0.,1.,0.,2.,0.,-1.,0.,0.,0.],
            [1.,0.,-2.,-2.,0.,-1.,0.,0.,0.],
            [0.,-1.,2.,0.,1.,-1.,0.,0.,0.],
            [1.,1.,0.,-2.,1.,-1.,0.,0.,0.],
            [1.,0.,-2.,2.,0.,-1.,0.,0.,0.],
            [2.,0.,0.,2.,0.,1.,0.,0.,0.],
            [0.,0.,2.,4.,2.,-1.,0.,0.,0.],
            [0.,1.,0.,1.,0.,1.,0.,0.,0.]]

    #给定日期与参考历元之间的时间间隔，儒略世纪数.
    T=((DATE1-DJ00)+DATE2)/DJC
    
    #FK5参考系中的基本参数
    #月球的平黄经减去月球近地点的平黄经。
    C=np.abs(1325.0*T)%1.0
    if (1325.0*T)<0:
        B=-C
    else:
        B=C
    A=(485866.733+(715922.633+(31.310+0.064*T)*T)*T)*DAS2R+B*D2PI
    EL=pymAnpm(A)
    
    #太阳的平黄经减去近地点的平黄经。
    C=np.abs(99.0*T)%1.0
    if (99.0*T)<0:
        B=-C
    else:
        B=C
    A=(1287099.804+(1292581.224+(-0.577-0.012*T)*T)*T)*DAS2R+B*D2PI
    ELP=pymAnpm(A)
    
    #月球的平黄经减去月球交点的平黄经。
    C=np.abs(1342.0*T)%1.0
    if (1342.0*T)<0:
        B=-C
    else:
        B=C
    A=(335778.877+(295263.137+(-13.257+0.011*T)*T)*T)*DAS2R+B*D2PI
    F=pymAnpm(A)
    
    #月亮与太阳之间的平角距
    C=np.abs(1236.0*T)%1.0
    if (1236.0*T)<0:
        B=-C
    else:
        B=C
    A=(1072261.307+(1105601.328+(-6.891+0.019*T)*T)*T)*DAS2R+B*D2PI
    D=pymAnpm(A)
    
    #月球轨道在黄道上的升交点平黄经，起量点为春分点
    C=np.abs(-5.0*T)%1.0
    if (-5.0*T)<0:
        B=-C
    else:
        B=C
    A=(450160.280+(-482890.539+(7.455+0.008*T)*T)*T)*DAS2R+B*D2PI
    OM=pymAnpm(A)
    
    # 章动部分
    
    #将时间单位由世纪转换为千年
    T=T/10.0
    
    #初始化章动参数
    DP=0.0
    DE=0.0
    
    #对章动项求和，从小到大
    for j in range(105,-1,-1):
        
        #构建当前项的参数
        ARG=float(X[j][0])*EL+float(X[j][1])*ELP+float(X[j][2])*F+\
            float(X[j][3])*D+float(X[j][4])*OM
            
        #累加当前章动项。
        S=float(X[j][5])+float(X[j][6])*T
        C=float(X[j][7])+float(X[j][8])*T
        if (S!=0.0):
            DP=DP+S*ma.sin(ARG)
        if (C!=0.0):
            DE=DE+C*ma.cos(ARG)
        
    #转换单位到弧度制
    DPSI=DP*U2R
    DEPS=DE*U2R
    
    return(DPSI,DEPS)
    

def pymEqeq94(DATE1,DATE2):
    '''
    Equation of the equinoxes, IAU 1994 model.

    Parameters
    ----------
    date1 : float
        TDB date    
    date2 : float
        TDB date

    Returns
    -------
    function value : float
        equation of the equinoxes

    '''
    #1角秒对应的弧度
    DAS2R=4.848136811095359935899141e-6

    #2Pi
    D2PI=6.283185307179586476925287
    
    #参考历元(J2000.0), JD
    DJ00=2451545.0

    #1儒略世纪的天数
    DJC=36525.0
    
    #给定日期和参考历元之间的时间间隔，儒略世纪数.
    T=((DATE1-DJ00)+DATE2)/DJC
    
    #月球轨道在黄道上的升交点平黄经，起量点为春分点
    C=np.abs(-5.0*T)%1.0
    if (-5.0*T)<0:
        B=-C
    else:
        B=C
    A=(450160.280+(-482890.539+(7.455+0.008*T)*T)*T)*DAS2R+B*D2PI
    OM=pymAnpm(A)
    
    #章动参数和平均倾角
    DPSI,DEPS=pymNut80(DATE1,DATE2)
    EPS0=pymObl80(DATE1,DATE2)
    
    #赤经章动
    EQEQ=DPSI*ma.cos(EPS0)+DAS2R*(0.00264*ma.sin(OM)+0.000063*ma.sin(OM+OM))
    
    return(EQEQ)


def pymRv2m(W):
    '''
    Form the r-matrix corresponding to a given r-vector.

    Parameters
    ----------
    w : list(3)     
        rotation vector     

    Returns
    -------
    r : list(3,3)     
        rotation matrix    

    '''
    R=[[0 for i in range(3)] for j in range(3)]
    
    #欧拉角（即为旋转向量的模）
    X=W[0]
    Y=W[1]
    Z=W[2]
    PHI=ma.sqrt(X*X+Y*Y+Z*Z)
    S=ma.sin(PHI)
    C=ma.cos(PHI)
    F=1.0-C 
    
    #欧拉轴(旋转向量的方向), 可能为空值.
    if (PHI>0.0):
        X=X/PHI
        Y=Y/PHI
        Z=Z/PHI
    
    #构建旋转矩阵
    R[0][0]=X*X*F+C
    R[0][1]=X*Y*F+Z*S
    R[0][2]=X*Z*F-Y*S
    R[1][0]=Y*X*F-Z*S
    R[1][1]=Y*Y*F+C
    R[1][2]=Y*Z*F+X*S
    R[2][0]=Z*X*F+Y*S
    R[2][1]=Z*Y*F-X*S
    R[2][2]=Z*Z*F+C
    
    return(R)


def pymRm2v(R):
    '''
    Express an r-matrix as an r-vector.

    Parameters
    ----------
    r : list(3,3)  
        rotation matrix

    Returns
    -------
    w : list(3)  
        rotation vector

    '''
    W=[0,0,0]
    
    X=R[1][2]-R[2][1]
    Y=R[2][0]-R[0][2]
    Z=R[0][1]-R[1][0]
    S2=ma.sqrt(X*X+Y*Y+Z*Z)
    if (S2>0.0):
        C2=R[0][0]+R[1][1]+R[2][2]-1.0
        PHI=ma.atan2(S2,C2)
        F=PHI/S2
        W[0]=X*F
        W[1]=Y*F
        W[2]=Z*F
    else:
        W[0]=0.0
        W[1]=0.0
        W[2]=0.0
    
    return(W)


def pymFk5hip():
    '''
    FK5 to Hipparcos rotation and spin.

    Returns
    -------
    r5h : list(3,3)
        r-matrix: FK5 rotation wrt Hipparcos    
    s5h : list(3)
        r-vector: FK5 spin wrt Hipparcos

    '''
    #1角秒对应的弧度
    DAS2R=4.848136811095359935899141e-6
    
    #从FK5到依巴谷的指向变化和绕轴旋转(弧度,弧度/年)
    EPX=-19.9e-3*DAS2R
    EPY=-9.1e-3*DAS2R
    EPZ=+22.9e-3*DAS2R
    
    OMX=-0.30e-3*DAS2R
    OMY=+0.60e-3*DAS2R
    OMZ=+0.70e-3*DAS2R
    
    #以旋转向量表示指向变化
    V=[0,0,0]
    V[0]=EPX
    V[1]=EPY
    V[2]=EPZ
    
    #将旋转向量转换为旋转矩阵
    R5H=pymRv2m(V)
    
    #从FK5到依巴谷的绕轴旋转矩阵
    S5H=[0,0,0]
    S5H[0]=OMX
    S5H[1]=OMY
    S5H[2]=OMZ
    
    return(R5H,S5H)


def pymFk5hz(R5,D5,DATE1,DATE2):
    '''
    Transform an FK5 (J2000.0) star position into the system of the
    Hipparcos catalog, assuming zero Hipparcos proper motion.

    Parameters
    ----------
    r5 : flaot
        FK5 RA (radians), equinox J2000.0, at date    
    d5 : flaot
        FK5 Dec (radians), equinox J2000.0, at date     
    date1 : flaot
        TDB date    
    date2 : flaot
        TDB date    

    Returns
    -------
    rh : float
        Hipparcos RA (radians)    
    dh : float
        Hipparcos Dec (radians)

    '''
    #参考历元(J2000.0), JD
    DJ00=2451545.0

    #1儒略世纪对应的天数
    DJY=36525.0
    
    #给定日期和参考历元之间的时间间隔，儒略年.
    T=-((DATE1-DJ00)+DATE2)/DJY
    
    #FK5中的质心位置向量
    P5E=pymS2c(R5,D5)
    
    #从FK5到依巴谷的旋转和绕轴旋转矩阵
    R5H,S5H=pymFk5hip()
    
    #在这个时间间隔上累计的依巴谷相对于FK5的绕轴旋转量
    VST=pymSxp(T,S5H)     

    #将累计的绕轴旋转量转换为旋转矩阵
    RST=pymRv2m(VST)    

    #通过旋转矩阵改正质心位置向量.
    P5=pymTrxp(RST,P5E)     

    #将这个向量旋转到依巴谷系中
    PH=pymRxp(R5H, P5)     

    #依巴谷中，向量到球面
    W,DH=pymC2s(PH)
    RH=pymAnp(W)
    
    return(RH,DH)


def pymPvu(DT,PV):
    '''
    Update a pv-vector.

    Parameters
    ----------
    dt : float
        time interval    
    pv : list(2,3) 
        pv-vector

    Returns
    -------
    upv : list(2,3) 
        p updated, v unchanged

    '''    
    UPV=[0,0]
    UPV[0]=pymPpsp(PV[0],DT,PV[1])
    UPV[1]=PV[1]
    
    return(UPV)


def pymFk45z (R1950,D1950,BEPOCH):
    '''
    Convert a B1950.0 FK4 star position to J2000.0 FK5, assuming zero
    proper motion in the FK5 system.

    Parameters
    ----------
    r1950 : float
        B1950.0 FK4 RA at epoch (rad)    
    d1950 : float
        B1950.0 FK4 Dec at epoch (rad)    
    bepoch : float
        Besselian epoch (e.g. 1979.3)

    Returns
    -------
    r2000 : float
        J2000.0 FK5 RA (rad)     
    d2000 : float
        J2000.0 FK5 Dec (rad)

    '''
    #2pi
    D2PI=6.283185307179586476925287
    
    #每百年弧度对应的角秒（1年为2pi）
    PMF=100.0*60.0*60.0*360.0/D2PI
    
    PV1=[[0 for i in range(3)] for j in range(2)]
    PV2=[[0 for i in range(3)] for j in range(2)]
    
    #规范的常量
    
    #向量 A 和 Adot (Seidelmann 3.591-2)
    A=[-1.62557e-6,-0.31919e-6,-0.13843e-6]
    AD=[-1.62557e-6,-0.31919e-6,-0.13843e-6]
    
    #位置-速度矩阵 (cf. Seidelmann 3.591-4, matrix M)
    EM=[0,0]
    EM[0]=[[+0.9999256782,-0.0111820611,-0.0048579477],
           [+0.0111820610,+0.9999374784,-0.0000271765],
           [+0.0048579479,-0.0000271474,+0.9999881997]]
    EM[1]=[[-0.000551,-0.238565,+0.435739],
           [+0.238514,-0.002667,-0.008541],
           [-0.435623,+0.012254,+0.002117]]
    
    #球面到笛卡尔坐标
    R0=pymS2c(R1950,D1950)
    
    #调整位置向量到0自行的FK5中
    W=(BEPOCH-1950.0)/PMF
    A1=pymPpsp(A,W,AD)
    
    #移除E项
    W=pymPdp(R0,A1)
    P1=pymPpsp(A1,-W,R0)
    P2=pymPmp(R0,P1)
    
    #转换到FK中的pv向量 (cf. Seidelmann 3.591-3).
    for k in range(2):
        for j in range(3):
            W=0.0
            for i in range(3):
                W=W+EM[k][j][i]*P2[i]
            PV1[k][j]=W
    
    #允许虚拟一个自行
    DJM0,DJM=pymEpb2jd(BEPOCH)
    W=(pymEpj(DJM0,DJM)-2000.0)/PMF
    PV2=pymPvu(W,PV1)
    P = PV2[0]
    #转换回球面
    W,D2000=pymC2s(P)
    R2000=pymAnp(W)
    
    return(R2000,D2000)


def pymPpp(A,B):
    '''
    P-vector addition.

    Parameters
    ----------
    a : list(3)
        first p-vector    
    b : list(3)
        second p-vector

    Returns
    -------
    apb : list(3)
        a + b

    '''
    APB=[A[0]+B[0],A[1]+B[1],A[2]+B[2]]
    
    return(APB)


def pymPv2s(PV):
    '''
    Convert position/velocity from Cartesian to spherical coordinates.

    Parameters
    ----------
    pv : list(2,3)
        pv-vector

    Returns
    -------
    theta : float
        longitude angle (radians)    
    phi : float
        latitude angle (radians)    
    r : float
        radial distance    
    td : float
        rate of change of theta    
    pd : float
        rate of change of phi    
    rd : float
        rate of change of r

    '''
    #位置-速度向量的各个分量
    X=PV[0][0]
    Y=PV[0][1]
    Z=PV[0][2]
    XD=PV[1][0]
    YD=PV[1][1]
    ZD=PV[1][2]

    #XY平面的平方和
    RXY2=X*X+Y*Y
    
    #模的平方
    R2=RXY2+Z*Z
    
    #模
    RTRUE=ma.sqrt(R2)
    
    #如果时空向量，沿着运动的方向移动原点
    RW=RTRUE
    if (RTRUE==0.0):
        X=XD
        Y=YD
        Z=ZD
        RXY2=X*X+Y*Y
        R2=RXY2+Z*Z
        RW=ma.sqrt(R2)
    
    #在球坐标中的位置和速度向量
    RXY=ma.sqrt(RXY2)
    XYP=X*XD+Y*YD
    if (RXY2!=0.0):
        THETA=ma.atan2(Y,X)
        PHI=ma.atan2(Z,RXY)
        TD=(X*YD-Y*XD)/RXY2
        PD=(ZD*RXY2-Z*XYP)/(R2*RXY)
    else:
        THETA=0.0
        if (Z!=0.0):
            PHI=ma.atan2(Z,RXY)
        else:
            PHI=0.0
        TD=0.0
        PD=0.0
    
    R=RTRUE
    if (RW!=0.0):
        RD=(XYP+Z*ZD)/RW
    else:
        RD=0.0
    
    return(THETA,PHI,R,TD,PD,RD)


def pymPvstar(PV):
    '''
    Convert star position+velocity vector to catalog coordinates.

    Parameters
    ----------
    pv : list(2,3)
        pv-vector (au, au/day)

    Returns
    -------
    ra : float
        right ascension (radians)    
    dec : float
        declination (radians)    
    pmr : float
        RA proper motion (radians/year)     
    pmd : float
        Dec proper motion (radians/year)    
    px : float
        parallax (arcsec)    
    rv : float
        radial velocity (km/s, positive = receding)
    J : ValueError
        0: 'OK',
       -1: 'superluminal speed',
       -2: 'null position vector' 
    '''
    RA,DEC,PMR,PMD,PX,RV=0,0,0,0,0,0
    
    #1儒略年的天数
    Y2D=365.25
    
    #1弧度对应的角秒
    DR2AS=206264.8062470963551564734
    
    #1天的秒数
    D2S=86400.0
    
    #光速(m/s)
    CMPS = 299792458.0
    
    #天文单位(m, IAU 2012)
    AUM=149597870.7e3
    
    #光速 (au per day)
    C=D2S*CMPS/AUM
 
    #分离出速度的径向分量 (au/day,惯性系).
    R,X=pymPn(PV[0])
    VR=pymPdp(X,PV[1])
    UR=pymSxp(VR,X)
    
    #分离出速度的横向分量 (au/day, 惯性系).
    UT=pymPmp(PV[1],UR)
    VT=pymPm(UT)
    
    #狭义相对论无量纲参数
    BETT=VT/C
    BETR=VR/C
    
    I=1
    while I<2:
    
        #惯性-观测修正项
        D=1.0+BETR
        W=BETR*BETR+BETT*BETT
        if (D==0.0)|(W>=1.0):
            J=-1
            print('ERROR',J)
            break
        DEL=-W/(ma.sqrt(1.0-W)+1.0)
        
        UST=pymSxp(1.0/D,UT)
        
        B=C*(BETR-DEL)/D
        USR=pymSxp(B, X)
        
        #将两者结合得到观测到的速度矢量(au/day).
        PV[1]=pymPpp(USR,UST)
    
        #笛卡尔坐标到球面
        A,DEC,R,RAD,DECD,RD=pymPv2s(PV)
        if (R==0.0):
            J=-2
            print('ERROR',J)
            break
        
        #将赤经规范化到0-2pi
        RA=pymAnp(A)
    
        #将自行单位转换为弧度每年
        PMR=RAD*Y2D
        PMD=DECD*Y2D
        
        #将视差单位转换为角秒
        PX=DR2AS/R
        
        #将视向速度单位转换为km/s.
        RV=1e-3*RD*AUM/D2S
        
        J=0
        
        I+=1
    
    return(RA,DEC,PMR,PMD,PX,RV,J)


def pymS2pv(THETA,PHI,R,TD,PD,RD):
    '''
    Convert position/velocity from spherical to Cartesian coordinates.

    Parameters
    ----------
    theta : float
        longitude angle (radians)    
    phi : float
        latitude angle (radians)    
    r : float
        radial distance    
    td : float
        rate of change of theta    
    pd : float
        rate of change of phi    
    rd : float
        rate of change of r

    Returns
    -------
    pv : list(2,3)
        pv-vector

    '''
    PV=[[0 for i in range(3)] for j in range(2)]
    
    ST=ma.sin(THETA)
    CT=ma.cos(THETA)
    SP=ma.sin(PHI)
    CP=ma.cos(PHI)
    RCP=R*CP
    X=RCP*CT
    Y=RCP*ST
    RPD=R*PD
    W=RPD*SP-CP*RD

    PV[0][0]=X
    PV[0][1]=Y
    PV[0][2]=R*SP
    PV[1][0]=-Y*TD-W*CT
    PV[1][1]=X*TD-W*ST
    PV[1][2]=RPD*CP+SP*RD

    return(PV)


def pymStarpv(RA,DEC,PMR,PMD,PX,RV):
    '''
    Convert star catalog coordinates to position+velocity vector.

    Parameters
    ----------
    ra : float
        right ascension (radians)    
    dec : float
        declination (radians)    
    pmr : float
        RA proper motion (radians/year)     
    pmd : float
        Dec proper motion (radians/year)    
    px : float
        parallax (arcsec)    
    rv : float
        radial velocity (km/s, positive = receding)

    Returns
    -------
    pv : list(2,3)
        pv-vector (au, au/day)
    J : ValueError
        0: 'no warnings',
        1: 'distance overridden',
        2: 'excessive speed',
        4: 'solution did not converge '    
    '''
    #允许的最小视差
    PXMIN=1e-7

    #允许的最大速度，单位：c
    VMAX=0.5

    #1儒略年的天数
    Y2D=365.25

    #弧度到角秒
    DR2AS=206264.8062470963551564734

    #1天的秒数
    D2S=86400.0

    #光速(m/s)
    CMPS=299792458.0

    #天文单位 (m, IAU 2012)
    AUM=149597870.7e3

    #光速(au per day)
    C=D2S*CMPS/AUM

    #相对论解的最大迭代次数
    IMAX=100
    
    #距离(au).
    if (PX>=PXMIN):
        W=PX
        IWARN=0
    else:
        W=IWARN
        IWARN=1
    R=DR2AS/W
    
    #视向速度(au/day).
    RD=D2S*RV*1e3/AUM

    #自行(radian/day).
    RAD=PMR/Y2D
    DECD=PMD/Y2D

    #转换到位置-速度向量(au,au/day).
    PV=pymS2pv(RA,DEC,R,RAD,DECD,RD)
    
    #当速度过大时，设为零
    V=pymPm(PV[1])
    A=V/C
    if (A>VMAX):
        PV[1]=pymZp()
        IWARN=IWARN+2
    
    #分离速度的径向分量(au/day).
    W,X=pymPn(PV[0])
    VSR=pymPdp(X,PV[1])
    USR=pymSxp(VSR,X)

    #分离速度的横向分量(au/day)
    UST=pymPmp(PV[1],USR)
    VST=pymPm(UST)
    
    #狭义无量纲参数
    BETSR=VSR/C
    BETST=VST/C
    
    #确定惯性-观测相对论修正项.
    OD=0.0
    ODEL=0.0
    ODD=0.0
    ODDEL=0.0
    BETT=BETST
    BETR=BETSR
    for i in range(IMAX):
        D=1.0+BETR
        W=BETR*BETR+BETT*BETT
        DEL=-W/(ma.sqrt(1.0-W)+1.0)
        BETR=D*BETSR+DEL
        BETT=D*BETST
        if (i>1):
            DD=np.abs(D-OD)
            DDEL=np.abs(DEL-ODEL)
            if ((i>2)&(DD>=ODD)&(DDEL>=ODDEL)):
                break
            if (i>=IMAX):
                IWARN=IWARN+4
            ODD=DD
            ODDEL=DDEL
        OD=D
        ODEL=DEL
    
    UT=pymSxp(D,UST)
    
    B=C*(D*BETSR+DEL)
    UR=pymSxp(B,X)
     
    #结合两者获得惯性空间速度.
    PV[1]=pymPpp(UR,UT)
        
    J=IWARN
    
    return(PV,J)


def pymFk52h(R5,D5,DR5,DD5,PX5,RV5):
    '''
    Transform FK5 (J2000.0) star data into the Hipparcos system.

    Parameters
    ----------
    all FK5, equinox J2000.0, epoch J2000.0     
    r5 : float
        RA (radians)    
    d5 : float
        Dec (radians)    
    dr5 : float
        proper motion in RA (dRA/dt, rad/Jyear)    
    dd5 : float
        proper motion in Dec (dDec/dt, rad/Jyear)    
    px5 : float
        parallax (arcsec)    
    rv5 : float
        radial velocity (km/s, positive = receding)

    Returns
    -------
    all Hipparcos, epoch J2000.0    
    rh : float
        RA (radians)    
    dh : float
        Dec (radians)    
    drh : float
        proper motion in RA (dRA/dt, rad/Jyear)    
    ddh : float
        proper motion in Dec (dDec/dt, rad/Jyear)    
    pxh : float
        parallax (arcsec)    
    rvh : float
        radial velocity (km/s, positive = receding)

    '''
    PVH=[[0 for i in range(3)] for j in range(2)]
    #FK5中的质心位置向量
    PV5,J=pymStarpv(R5,D5,DR5,DD5,PX5,RV5)
    
    #从FK5到依巴谷的旋转和绕轴旋转矩阵
    R5H,S5H=pymFk5hip()

    #将绕轴旋转的单位从每天变换到每年
    for i in range(3):
        S5H[i]=S5H[i]/365.25
        
    #将FK5的位置转换到依巴谷中
    PVH[0]=pymRxp(R5H,PV5[0])
 
    #应用旋转的位置，给予额外的空间运动参数
    WXP=pymPxp(PV5[0],S5H)
    
    #将此参数添加到FK5空间运动中
    VV=pymPpp(WXP,PV5[1])

    #将FK5的空间运动转换到依巴谷中
    PVH[1]=pymRxp(R5H,VV)
    
    #将依巴谷的位置-速度向量转换到球面
    RH,DH,DRH,DDH,PXH,RVH,J=pymPvstar(PVH)
    
    return(RH,DH,DRH,DDH,PXH,RVH)


def pymFk524(R2000,D2000,DR2000,DD2000,P2000,V2000):
    '''
    Convert J2000.0 FK5 star catalog data to B1950.0 FK4.

    Parameters
    ----------
    all J2000.0, FK5    
    r2000 : float
        J2000.0 RA (rad)    
    d2000 : float
        J2000.0 Dec (rad)    
    dr2000 : float
        J2000.0 proper motions (rad/Jul.yr)     
    dd2000 : float
        J2000.0 proper motions (rad/Jul.yr)     
    p2000 : float
        parallax (arcsec)    
    v2000 : float
        radial velocity (km/s, +ve = moving away)
        
    Returns
    -------
    all B1950.0, FK4    
    r1950 : float
        B1950.0 RA (rad)    
    d1950 : float
        B1950.0 Dec (rad)    
    dr1950 : float
        B1950.0 proper motions (rad/trop.yr)     
    dd1950 : float
        B1950.0 proper motions (rad/trop.yr)     
    p1950 : float
        parallax (arcsec)    
    v1950 : float
        radial velocity (km/s, +ve = moving away)

    '''
    R1=[[0 for i in range(3)] for j in range(2)]
    PV=[[0 for i in range(3)] for j in range(2)]
    #2Pi
    D2PI=6.283185307179586476925287
    
    #每百年弧度对应的角秒（1年为2pi）
    PMF=100.0*60.0*60.0*360.0/D2PI
    
    #最小值，以免出错
    TINY = 1e-30
    
    #规范的常量  (Seidelmann 1992)
    #公里每秒到AU每回归世纪的变换= 86400 * 36524.2198782 / 149597870
    VF=21.095
    
    #位置速度向量常数(cf. Seidelmann 3.591-2, vectors A and Adot)
    A=[[-1.62557e-6,-0.31919e-6,-0.13843e-6],
       [+1.245e-3,-1.580e-3,-0.659e-3]]
    
    #位置-速度向量(cf. Seidelmann 3.592-1, matrix M^-1)
    EMI=[[[[ +0.9999256795,     +0.0111814828,     +0.0048590039,    ],
           [ -0.00000242389840, -0.00000002710544, -0.00000001177742 ]],
          [[ -0.0111814828,     +0.9999374849,     -0.0000271771,    ],
           [ +0.00000002710544, -0.00000242392702, +0.00000000006585 ]],
          [[ -0.0048590040,     -0.0000271557,     +0.9999881946,    ],
           [ +0.00000001177742, +0.00000000006585, -0.00000242404995 ]]],
         [[[ -0.000551,         +0.238509,         -0.435614,        ],
           [ +0.99990432,       +0.01118145,       +0.00485852       ]],
          [[ -0.238560,         -0.002667,         +0.012254,        ],
           [ -0.01118145,       +0.99991613,       -0.00002717       ]],
          [[ +0.435730,         -0.008541,         +0.002117,        ],
           [ -0.00485852,       -0.00002716,       +0.99996684       ]]]]
    
    #FK5数据，单位：弧度，角秒/儒略世纪
    R=R2000
    D=D2000
    UR=DR2000*PMF
    UD=DD2000*PMF
    PX=P2000
    RV=V2000
    
    #表示成pv向量
    PXVF=PX*VF
    W=RV*PXVF
    R0=pymS2pv(R,D,1.0,UR,UD,W)
    
    #将pv-vector转换为 Bessel-Newcomb system (cf. Seidelmann 3.592-1).
    for L in range(2):
        for K in range(3):
            W=0.0
            for J in range(2):
                for I in range(3):
                    W=W+EMI[L][K][J][I]*R0[J][I]
            R1[L][K]=W
    
    #应用E项(equivalent to Seidelmann 3.592-3, one iteration).
    #方向.
    WR=pymPm(R1[0])
    W=pymPdp(R1[0],A[0])
    P1=pymSxp(W,R1[0])
    P2=pymSxp(WR,A[0])
    P1=pymPmp(P2,P1)
    P1=pymPpp(R1[0],P1)
    
    #再计算长度
    WR=pymPm(P1)
    
    #方向.
    W=pymPdp(R1[0],A[0])
    P1=pymSxp(W,R1[0])
    P2=pymSxp(WR,A[0])
    P1=pymPmp(P2,P1)
    PV[0]=pymPpp(R1[0],P1)
    
    #微商
    W=pymPdp(R1[0],A[1])
    P1=pymSxp(W,PV[0])
    P2=pymSxp(WR,A[1])
    P1=pymPmp(P2,P1)
    PV[1]=pymPpp(R1[1],P1)
    
    #返回星表形式
    R,D,W,UR,UD,RD=pymPv2s(PV)
    if (PX>TINY):
        RV=RD/PXVF
        PX=PX/W
    
    #给出结果.
    R1950=pymAnp(R)
    D1950=D
    DR1950=UR/PMF
    DD1950=UD/PMF
    P1950=PX
    V1950=RV
    
    return(R1950,D1950,DR1950,DD1950,P1950,V1950)


def pymPvmpv(A,B):
    '''
    Subtract one pv-vector from another.

    Parameters
    ----------
    a : list(2,3)
        first pv-vector    
    b : list(2,3)
        second pv-vector

    Returns
    -------
    amb : list(2,3)
        a - b

    '''
    AMB=[0,0]
    for i in range(2):
        AMB[i]=pymPmp(A[i], B[i])
    
    return(AMB)


def pymPvppv(A,B):
    '''
    Add one pv-vector to another.

    Parameters
    ----------
    a : list(2,3) 
        first pv-vector    
    b : list(2,3) 
        second pv-vector

    Returns
    -------
    apb : list(2,3) 
        a + b

    '''
    APB=[0,0]
    for i in range(2):
        APB[i]=pymPpp(A[i], B[i])
    
    return(APB)


def pymFk425(R1950,D1950,DR1950,DD1950,P1950,V1950):
    '''
    Convert B1950.0 FK4 star catalog data to J2000.0 FK5.

    Parameters
    ----------
    all B1950.0, FK4    
    r1950 : float
        B1950.0 RA (rad)    
    d1950 : float
        B1950.0 Dec (rad)    
    dr1950 : float
        B1950.0 proper motions (rad/trop.yr)     
    dd1950 : float
        B1950.0 proper motions (rad/trop.yr)     
    p1950 : float
        parallax (arcsec)    
    v1950 : float
        radial velocity (km/s, +ve = moving away)

    Returns
    -------
    all J2000.0, FK5    
    r2000 : float
        J2000.0 RA (rad)    
    d2000 : float
        J2000.0 Dec (rad)    
    dr2000 : float
        J2000.0 proper motions (rad/Jul.yr)     
    dd2000 : float
        J2000.0 proper motions (rad/Jul.yr)     
    p2000 : float
        parallax (arcsec)    
    v2000 : float
        radial velocity (km/s, +ve = moving away)

    '''
    R0=[[0 for i in range(3)] for j in range(2)]
    PV1=[[0 for i in range(3)] for j in range(2)]
    PV2=[[0 for i in range(3)] for j in range(2)]
    
    #2Pi
    D2PI=6.283185307179586476925287
    
    #每百年弧度对应的角秒（1年为2pi）
    PMF=100.0*60.0*60.0*360.0/D2PI
    
    #最小值，以免出错
    TINY = 1e-30
    
    #规范的常量  (Seidelmann 1992)
    #公里每秒到AU每回归世纪的变换= 86400 * 36524.2198782 / 149597870
    VF=21.095
    
    #位置速度向量常数(cf. Seidelmann 3.591-2, vectors A and Adot)
    A=[[-1.62557e-6, -0.31919e-6, -0.13843e-6],
       [ +1.245e-3,   -1.580e-3,   -0.659e-3 ]]
    
    #位置-速度向量(cf. Seidelmann 3.592-1, matrix M^-1)    
    EM=[[[[ +0.9999256782,     -0.0111820611,     -0.0048579477     ],
          [ +0.00000242395018, -0.00000002710663, -0.00000001177656 ]],
         [[ +0.0111820610,     +0.9999374784,     -0.0000271765     ],
          [ +0.00000002710663, +0.00000242397878, -0.00000000006587 ]],
         [[ +0.0048579479,     -0.0000271474,     +0.9999881997,    ],
          [ +0.00000001177656, -0.00000000006582, +0.00000242410173 ]]],
        [[[ -0.000551,         -0.238565,         +0.435739        ],
          [ +0.99994704,       -0.01118251,       -0.00485767       ]],
         [[ +0.238514,         -0.002667,         -0.008541        ],
          [ +0.01118251,       +0.99995883,       -0.00002718       ]],
         [[ -0.435623,         +0.012254,         +0.002117         ],
          [ +0.00485767,       -0.00002714,       +1.00000956       ]]]]

    
    #FK4数据，单位：弧度，角秒/儒略世纪
    R=R1950
    D=D1950
    UR=DR1950*PMF
    UD=DD1950*PMF
    PX=P1950
    RV=V1950
    
    #表示成pv向量
    PXVF=PX*VF
    W=RV*PXVF
    R0=pymS2pv(R,D,1.0,UR,UD,W)
    
    #E项 (cf. Seidelmann 3.591-2).
    PV1=pymPvmpv(R0,A)
    W=pymPdp(R0[0],A[0])
    PV2[0]=pymSxp(W,R0[0])
    W=pymPdp(R0[0],A[1])
    PV2[1]=pymSxp(W,R0[0])
    PV1=pymPvppv(PV1,PV2)
    
    #将位置-速度向量转换到Fricke system (cf. Seidelmann 3.591-3).
    for I in range(2):
        for J in range(3):
            W=0.0
            for K in range(2):
                for L in range(3):
                    W=W+EM[I][J][K][L]*PV1[K][L]
            PV2[I][J]=W
    
    #返回到星表形式
    R,D,W,UR,UD,RD=pymPv2s(PV2)
    if (PX>TINY):
        RV=RD/PXVF
        PX=PX/W
    
    #给出结果
    R2000=pymAnp(R)
    D2000=D
    DR2000=UR/PMF
    DD2000=UD/PMF
    V2000=RV
    P2000=PX
    
    return(R2000,D2000,DR2000,DD2000,P2000,V2000)


def pymFk54z(R2000,D2000,BEPOCH):
    '''
    Convert a J2000.0 FK5 star position to B1950.0 FK4, assuming zero
    proper motion in FK5 and parallax.

    Parameters
    ----------
    r2000 : float
        J2000.0 FK5 RA (rad)    
    d2000 : float
        J2000.0 FK5 Dec (rad)    
    bepoch : float
        Besselian epoch (e.g. 1950.0)

    Returns
    -------
    r1950 : float
        B1950.0 FK4 RA (rad) at epoch BEPOCH    
    d1950 : float
        B1950.0 FK4 Dec (rad) at epoch BEPOCH    
    dr1950 : float
        B1950.0 FK4 proper motions (rad/trop.yr)    
    dd1950 : float
        B1950.0 FK4 proper motions (rad/trop.yr)    

    '''    
    V=[0,0,0]
    #FK5 J2000.0春分点到FK4 B1950.0春分点.
    R,D,PR,PD,PX,RV=pymFk524(R2000,D2000,0.0,0.0,0.0,0.0)
    
    #球面到笛卡尔坐标
    P=pymS2c(R,D)
    
    #自行
    V[0]=-PR*P[1]-PD*ma.cos(R)*ma.sin(D)
    V[1]=PR*P[0]-PD*ma.sin(R)*ma.sin(D)
    V[2]=PD*ma.cos(D)
    
    #应用自行
    W =BEPOCH-1950.0
    for i in range(3):
        P[i]=P[i]+W*V[i]

    #笛卡尔到球面
    W,D1950=pymC2s(P)
    R1950=pymAnp(W)
    
    #自行
    DR1950=PR
    DD1950=PD
    
    return(R1950,D1950,DR1950,DD1950)


def pymFw2xy(GAMB,PHIB,PSI,EPS):
    '''
    CIP X,Y given Fukushima-Williams bias-precession-nutation angles.

    Parameters
    ----------
    gamb : float
        F-W angle gamma_bar (radians)    
    phib : float
        F-W angle phi_bar (radians)    
    psi : float
        F-W angle psi (radians)    
    eps : float
        F-W angle epsilon (radians)

    Returns
    -------
    x : float
        CIP unit vector X,Y     
    y : float
        CIP unit vector X,Y

    '''
    #构建NxPxB矩阵.
    R=pymFw2m(GAMB,PHIB,PSI,EPS)
    
    #从CIP中分解出X,Y
    X,Y=pymBpn2xy(R)
    
    return(X,Y)


def pymG2icrs(DL,DB):
    '''
    Transformation from Galactic Coordinates to ICRS.

    Parameters
    ----------
    dl : float
        galactic longitude (radians)    
    db : float
        galactic latitude (radians)

    Returns
    -------
    dr : float
        ICRS right ascension (radians)    
    dd : float
        ICRS declination (radians)

    '''
    # *  L2,B2 system of galactic coordinates in the form presented in the
    # *  Hipparcos catalog.  In degrees:
    # *
    # *  P = 192.85948    right ascension of the Galactic north pole in ICRS
    # *  Q =  27.12825    declination of the Galactic north pole in ICRS
    # *  R =  32.93192    Galactic longitude of the ascending node of
    # *                   the Galactic equator on the ICRS equator
    # *
    # *  ICRS to galactic rotation matrix, obtained by computing
    # *  R_3(-R) R_1(pi/2-Q) R_3(pi/2+P) to the full precision shown:
    R=[[-0.054875560416215368492398900454,
        -0.873437090234885048760383168409,
        -0.483835015548713226831774175116],
       [+0.494109427875583673525222371358,
        -0.444829629960011178146614061616,
        +0.746982244497218890527388004556],
       [-0.867666149019004701181616534570,
        -0.198076373431201528180486091412,
        +0.455983776175066922272100478348]]
    
    #球面到笛卡尔
    V1=pymS2c(DL,DB)
    
    #银道到ICRS
    V2=pymTrxp(R,V1)
    
    #笛卡尔到球面
    DR,DD=pymC2s(V2)      
   
    #用常规范围表示角度
    DR=pymAnp(DR)
    DD=pymAnpm(DD)
    
    return(DR,DD)


def pymGmst82(DJ1,DJ2):
    '''
    Universal Time to Greenwich mean sidereal time (IAU 1982 model).

    Parameters
    ----------
    dj1 : float
        UT1 Julian Date    
    dj2 : float
        UT1 Julian Date

    Returns
    -------
    function value : float
        Greenwich mean sidereal time (radians)

    '''
    #1秒对应的弧度
    DS2R=7.272205216643039903848712e-5
    
    #参考历元(J2000.0), JD
    DJ00=2451545.0
    
    #1天的秒数，1世纪的天数
    D2S=86400.0
    DJC=36525.0
    
    #IAU 1982 GMST-UT1模型的系数
    A=24110.54841-D2S/2.0
    B=8640184.812866
    C=0.093104
    D=-6.2e-6
    #UT1的一天的起始点为当天正午，因此第一个系数A要调整12个小时 
    #给定时间与参考历元之间的时间间隔，儒略世纪数
    if (DJ1<DJ2):
        D1=DJ1
        D2=DJ2
    else:
        D1=DJ2
        D2=DJ1
    T=(D1+(D2-DJ00))/DJC
    
    #儒略日的小数部分(UT1), 单位：秒.
    C1=np.abs(D1)%1.0
    C2=np.abs(D2)%1.0
    if D1<0:
        B1=-C1
    else:
        B1=C1
    if D2<0:
        B2=-C2
    else:
        B2=C2
    F=D2S*(B1+B2)
    
    #给定时刻的格林尼治平恒星时
    GMST=pymAnp(DS2R*((A+(B+(C+D*T)*T)*T)+F))
    
    return(GMST)


def pymGst00a(UTA,UTB,TTA,TTB):
    '''
    Greenwich apparent sidereal time (consistent with IAU 2000
    resolutions).

    Parameters
    ----------
    uta : float
        UT1 as a 2-part Julian Date    
    utb : float
        UT1 as a 2-part Julian Date    
    tta : float
        TT as a 2-part Julian Date    
    ttb : float
        TT as a 2-part Julian Date

    Returns
    -------
    GST : float
        Greenwich apparent sidereal time (radians)

    '''
    A=pymGmst00(UTA,UTB,TTA,TTB)
    B=pymEe00a(TTA,TTB)
    GST=pymAnp(A+B)
    
    return(GST)


def pymGst00b(UTA,UTB):
    '''
    Greenwich apparent sidereal time (consistent with IAU 2000
    resolutions but using the truncated nutation model IAU 2000B).

    Parameters
    ----------
    uta : float
        UT1 as a 2-part Julian Date    
    utb : float
        UT1 as a 2-part Julian Date

    Returns
    -------
    GST : float
        Greenwich apparent sidereal time (radians)

    '''
    A=pymGmst00(UTA,UTB,UTA,UTB)
    B=pymEe00b(UTA,UTB)
    GST=pymAnp(A+B)

    return(GST)


def pymGst94(UTA,UTB):
    '''
    Greenwich apparent sidereal time (consistent with IAU 1982/94
    resolutions).

    Parameters
    ----------
    uta : float
        UT1 as a 2-part Julian Date    
    utb : float
        UT1 as a 2-part Julian Date

    Returns
    -------
    GST : float
        Greenwich apparent sidereal time (radians)

    '''
    A=pymGmst82(UTA,UTB)
    B=pymEqeq94(UTA,UTB)
    GST=pymAnp(A+B)

    return(GST)
    

def pymH2fk5 (RH,DH,DRH,DDH,PXH,RVH):
    '''
    Transform Hipparcos star data into the FK5 (J2000.0) system.

    Parameters
    ----------
    all Hipparcos, epoch J2000.0    
    rh : float
        RA (radians)    
    dh : float
        Dec (radians)    
    drh : float
        proper motion in RA (dRA/dt, rad/Jyear)    
    ddh : float
        proper motion in Dec (dDec/dt, rad/Jyear)    
    pxh : float
        parallax (arcsec)    
    rvh : float
        radial velocity (km/s, positive = receding)
        
    Returns
    -------
    all FK5, equinox J2000.0, epoch J2000.0     
    r5 : float
        RA (radians)    
    d5 : float
        Dec (radians)    
    dr5 : float
        proper motion in RA (dRA/dt, rad/Jyear)    
    dd5 : float
        proper motion in Dec (dDec/dt, rad/Jyear)    
    px5 : float
        parallax (arcsec)    
    rv5 : float
        radial velocity (km/s, positive = receding)

    '''
    PV5=[0,0]
    #依巴谷中的质心位置-速度向量
    PVH,J=pymStarpv(RH,DH,DRH,DDH,PXH,RVH)     

    #FK5到依巴谷的旋转矩阵和绕轴旋转向量
    R5H,S5H=pymFk5hip()

    #将绕轴旋转的单位由每年转换为每天
    for i in range(3):
        S5H[i]=S5H[i]/365.25

    #将绕轴旋转定向到依巴谷星表中
    SH=pymRxp(R5H,S5H)
    
    #重新定向依巴谷中的位置到FK5中
    PV5[0]=pymTrxp(R5H,PVH[0])    

    #将绕轴旋转应用到位置向量上，给出额外的空间运动参数
    WXP=pymPxp(PVH[0],SH)     

    #从依巴谷空间运动中减去这个分量
    VV=pymPmp(PVH[1],WXP)     

    #重新定向依巴谷中的空间运动到FK5中
    PV5[1]=pymTrxp(R5H,VV)

    #FK5位置速度向量转换到球面
    R5,D5,DR5,DD5,PX5,RV5,J=pymPvstar(PV5)
    
    return(R5,D5,DR5,DD5,PX5,RV5)


def pymHd2ae(HA,DEC,PHI):
    '''
    Equatorial to horizon coordinates:  transform hour angle and
    declination to azimuth and altitude.

    Parameters
    ----------
    ha : float
        hour angle (local)    
    dec : float
        declination    
    phi : float
        site latitude

    Returns
    -------
    az : float
        azimuth    
    el : float
        altitude (informally, elevation)

    '''
    D2PI=6.283185307179586476925287
    
    SH=ma.sin(HA)
    CH=ma.cos(HA)
    SD=ma.sin(DEC)
    CD=ma.cos(DEC)
    SP=ma.sin(PHI)
    CP=ma.cos(PHI)
    
    #单位向量
    X=-CH*CD*SP+SD*CP
    Y=-SH*CD
    Z=CH*CD*CP+SD*SP
    
    #到球面
    R=ma.sqrt(X*X+Y*Y)
    if (R==0.0):
        A=0.0
    else:
        A=ma.atan2(Y,X)
    
    if (A<0.0):
        A=A+D2PI
    AZ=A
    EL=ma.atan2(Z,R)
    
    return(AZ,EL)


def pymHd2pa(HA,DEC,PHI):
    '''
    Parallactic angle for a given hour angle and declination.

    Parameters
    ----------
    ha : float
        hour angle    
    dec : float
        declination    
    phi : float
        site latitude

    Returns
    -------
    function value : flaot
        parallactic angle

    '''
    CP=ma.cos(PHI)
    SQSZ=CP*ma.sin(HA)
    CQSZ=ma.sin(PHI)*ma.cos(DEC)-CP*ma.sin(DEC)*ma.cos(HA)
    if (SQSZ==0.0)&(CQSZ==0.0):
        CQSZ=1.0
    HD2PA=ma.atan2(SQSZ,CQSZ)
    
    return(HD2PA)


def pymHfk5z(RH,DH,DATE1,DATE2):
    '''
    Transform a Hipparcos star position into FK5 J2000.0, assuming
    =zero Hipparcos proper motion.

    Parameters
    ----------
    rh : float
        Hipparcos RA (radians)    
    dh : float
        Hipparcos Dec (radians)    
    date1 : float
        TDB date    
    date2 : float
        TDB date

    Returns
    -------
    all FK5, equinox J2000.0, date date1+date2      
    r5 : float
        RA (radians)    
    d5 : float
        Dec (radians)    
    dr5 : float
        proper motion in RA (dRA/dt, rad/Jyear)    
    dd5 : float
        proper motion in Dec (dDec/dt, rad/Jyear)    

    '''
    PV5E=[0,0]
    #参考历元(J2000.0), JD
    DJ00=2451545.0
    
    #1儒略年的天数
    DJY=365.25
    
    #给定日期到参考历元之间的时间间隔，儒略年.
    T=((DATE1-DJ00)+DATE2)/DJY

    #依巴谷中的质心位置向量
    PH=pymS2c(RH,DH)

    #FK5到依巴谷的旋转矩阵和绕轴旋转向量
    R5H,S5H=pymFk5hip()

    #将绕轴旋转添加到依巴谷系统中
    SH=pymRxp(R5H,S5H)

    #累计依巴谷相对于FK5的绕轴旋转时间间隔。
    VST=pymSxp(T,S5H)

    #将累积的绕轴旋转表示为旋转矩阵.
    RST=pymRv2m(VST)

    #旋转矩阵:累积的绕轴旋转，从FK5到依巴谷。
    R5HT=pymRxr(R5H,RST)

    #重新定向和旋转，将依巴谷中的位置转换到 FK5 J2000.0.
    PV5E[0]=pymTrxp(R5HT,PH)

    #通过绕轴旋转给出一个空间运动
    VV=pymPxp(SH,PH)

    #重新定向和旋转，将依巴谷中的速度转换到 FK5 J2000.0.
    PV5E[1]=pymTrxp(R5HT,VV)

    #FK5 位置-速度向量转换到球面
    W,D5,R,DR5,DD5,V=pymPv2s(PV5E)
    R5=pymAnp(W)
    
    return(R5,D5,DR5,DD5)


def pymIcrs2g(DR,DD):
    '''
    Transformation from ICRS to Galactic Coordinates.

    Parameters
    ----------
    dr : float
        ICRS right ascension (radians)    
    dd : float
        ICRS declination (radians)

    Returns
    -------
    dl : float
        galactic longitude (radians)    
    db : float
        galactic latitude (radians)

    '''
    # *  L2,B2 system of galactic coordinates in the form presented in the
    # *  Hipparcos catalog.  In degrees:
    # *
    # *  P = 192.85948    right ascension of the Galactic north pole in ICRS
    # *  Q =  27.12825    declination of the Galactic north pole in ICRS
    # *  R =  32.93192    Galactic longitude of the ascending node of
    # *                   the Galactic equator on the ICRS equator
    # *
    # *  ICRS to galactic rotation matrix, obtained by computing
    # *  R_3(-R) R_1(pi/2-Q) R_3(pi/2+P) to the full precision shown:
    R=[[-0.054875560416215368492398900454,
        -0.873437090234885048760383168409,
        -0.483835015548713226831774175116],
       [+0.494109427875583673525222371358,
        -0.444829629960011178146614061616,
        +0.746982244497218890527388004556],
       [-0.867666149019004701181616534570,
        -0.198076373431201528180486091412,
        +0.455983776175066922272100478348]]
    
    #球面到笛卡尔
    V1=pymS2c(DR,DD)

    #ICRS到银道.
    V2=pymRxp(R,V1)

    #笛卡尔到球面
    DL,DB=pymC2s(V2)

    #用常规范围表示角度
    DL=pymAnp(DL)
    DB=pymAnpm(DB)
    
    return(DL,DB)


def pymLtpequ(EPJ):
    '''
    Long-term precession of the equator.

    Parameters
    ----------
    epj : float
        Julian epoch (TT)   

    Returns
    -------
    veq : list(3)
        equator pole unit vector

    '''
    VEQ=[0,0,0]
    
    #1角秒对应的弧度
    AS2R=4.848136811095359935899141e-6
    
    #2Pi
    D2PI=6.283185307179586476925287
    
    #多项式数目
    NPOL=4
    
    #周期项数目
    NPER=14
    
    #多项式
    XYPOL=[[+5453.282155,+0.4252841,-0.00037173,-0.000000152],
           [-73750.930350,-0.7675452,-0.00018725,+0.000000231]]
    
    #周期
    XYPER=[[256.75,-819.940624,75004.344875,81491.287984,1558.515853],
           [708.15,-8444.676815,624.033993,787.163481,7774.939698],
           [274.20,2600.009459,1251.136893,1251.296102,-2219.534038],
           [241.45,2755.175630,-1102.212834,-1257.950837,-2523.969396],
           [2309.00,-167.659835,-2660.664980,-2966.799730,247.850422],
           [492.20,871.855056,699.291817,639.744522,-846.485643],
           [396.10,44.769698,153.167220,131.600209,-1393.124055],
           [288.90,-512.313065,-950.865637,-445.040117,368.526116],
           [231.10,-819.415595,499.754645,584.522874,749.045012],
           [1610.00,-538.071099,-145.188210,-89.756563,444.704518],
           [620.00,-189.793622,558.116553,524.429630,235.934465],
           [157.87,-402.922932,-23.923029,-13.549067,374.049623],
           [220.30,179.516345,-165.405086,-210.157124,-171.330180],
           [1200.00,-9.814756,9.344131,-44.919798,-22.899655]]
    
    #自J2000的世纪数.
    T=(EPJ-2000.0)/100.0
    
    #初始化X,Y
    X=0.0
    Y=0.0
    
    #周期项
    W=D2PI*T
    for i in range(NPER):
        A=W/XYPER[i][0]
        S=ma.sin(A)
        C=ma.cos(A)
        X=X+C*XYPER[i][1]+S*XYPER[i][3]
        Y=Y+C*XYPER[i][2]+S*XYPER[i][4]
    
    #多项式
    W=1.0
    for i in range(NPOL):
        X=X+XYPOL[0][i]*W
        Y=Y+XYPOL[1][i]*W
        W=W*T
    
    #X,Y（方向余弦）
    X=X*AS2R
    Y=Y*AS2R

    #构建赤道极向量
    VEQ[0]=X
    VEQ[1]=Y
    VEQ[2]=ma.sqrt(max(1.0-X*X-Y*Y,0.0))
    
    return(VEQ)


def pymLtpecl(EPJ):
    '''
    Long-term precession of the ecliptic.

    Parameters
    ----------
    epj : float
        Julian epoch (TT)    
        
    Returns
    -------
    vec : list(3)
        ecliptic pole unit vector    

    '''
    VEC=[0,0,0]
    
    #1角秒对应的弧度
    DAS2R=4.848136811095359935899141e-6
    
    #2Pi
    D2PI=6.283185307179586476925287
    
    #J2000.0 时的倾角，弧度制.
    EPS0=84381.406*DAS2R
    
    #多项式数目
    NPOL=4
    
    #周期项数目
    NPER=8
    
    #多项式
    PQPOL=[[+5851.607687,-0.1189000,-0.00028913,+0.000000101],
           [-1600.886300,+1.1689818,-0.00000020,-0.000000437]]
    
    #周期
    PQPER=[[708.15,-5486.751211,-684.661560,667.666730,-5523.863691],
           [2309.00,-17.127623,2446.283880,-2354.886252,-549.747450],
           [1620.00,-617.517403,399.671049,-428.152441,-310.998056],
           [492.20,413.442940,-356.652376,376.202861,421.535876],
           [1183.00,78.614193,-186.387003,184.778874,-36.776172],
           [622.00,-180.732815,-316.800070,335.321713,-145.278396],
           [882.00,-87.676083,198.296701,-185.138669,-34.744450],
           [547.00,46.140315,101.135679,-120.972830,22.885731]]
    
    #自J2000的世纪数.
    T=(EPJ-2000.0)/100.0
    
    #初始化P,Q
    P=0.0
    Q=0.0
    
    #周期项
    W=D2PI*T
    for i in range(NPER):
        A=W/PQPER[i][0]
        S=ma.sin(A)
        C=ma.cos(A)
        P=P+C*PQPER[i][1]+S*PQPER[i][3]
        Q=Q+C*PQPER[i][2]+S*PQPER[i][4]
    
    #多项式
    W = 1.0
    for i in range(NPOL):
        P=P+PQPOL[0][i]*W
        Q=Q+PQPOL[1][i]*W
        W=W*T
    
    #P,Q，弧度制.
    P=P*DAS2R
    Q=Q*DAS2R
    
    #构建黄道极向量
    W=ma.sqrt(max(1.0-P*P-Q*Q,0.0))
    S=ma.sin(EPS0)
    C=ma.cos(EPS0)
    VEC[0]=P
    VEC[1]=-Q*C-W*S
    VEC[2]=-Q*S+W*C
    
    return(VEC)


def pymLtecm(EPJ):
    '''
    ICRS equatorial to ecliptic rotation matrix, long-term.

    Parameters
    ----------
    epj : float
        Julian epoch (TT)

    Returns
    -------
    rm : list(3,3)
        ICRS to ecliptic rotation matrix

    '''
    RM=[[0 for i in range(3)] for j in range(3)]
    
    #1角秒对应的弧度
    AS2R=4.848136811095359935899141e-6
    
    #参考架偏差(IERS Conventions 2010, Eqs. 5.21 and 5.33)
    DX=-0.016617*AS2R
    DE=-0.0068192*AS2R
    DR=-0.0146*AS2R
    
    #赤道极
    P=pymLtpequ(EPJ)
    
    #黄道极
    Z=pymLtpecl(EPJ)
    
    #春分点(top row of matrix).
    W=pymPxp(P,Z)
    S,X=pymPn(W)
    
    #中间一行
    Y=pymPxp(Z,X)
    
    #结合参考架偏差
    RM[0][0]=X[0]-X[1]*DR+X[2]*DX
    RM[0][1]=X[0]*DR+X[1]+X[2]*DE
    RM[0][2]=-X[0]*DX-X[1]*DE+X[2]
    RM[1][0]=Y[0]-Y[1]*DR+Y[2]*DX
    RM[1][1]=Y[0]*DR+Y[1]+Y[2]*DE
    RM[1][2]=-Y[0]*DX-Y[1]*DE+Y[2]
    RM[2][0]=Z[0]-Z[1]*DR+Z[2]*DX
    RM[2][1]=Z[0]*DR+Z[1]+Z[2]*DE
    RM[2][2]=-Z[0]*DX-Z[1]*DE+Z[2]

    return(RM)


def pymLteceq(EPJ,DL,DB):
    '''
    Transformation from ecliptic coordinates (mean equinox and ecliptic
    of date) to ICRS RA,Dec, using a long-term precession model.

    Parameters
    ----------
    epj : float
        Julian epoch (TT)
    dl : float
        ecliptic longitude (radians)    
    db : float
        ecliptic latitude (radians)

    Returns
    -------
    dr : float
        ICRS right ascension (radians)     
    dd : float
        ICRS right declination (radians)

    '''
    #球面到笛卡尔
    V1=pymS2c(DL, DB)

    #旋转矩阵, ICRS 赤道到黄道
    RM=pymLtecm(EPJ)

    #从黄道转换到 ICRS赤道.
    V2=pymTrxp(RM,V1)

    #笛卡尔到球面
    A,B=pymC2s(V2)

    #规范化表达
    DR=pymAnp(A)
    DD=pymAnpm(B)
    
    return(DR,DD)


def pymLteqec(EPJ,DR,DD):
    '''
    Transformation from ICRS equatorial coordinates to ecliptic
    coordinates (mean equinox and ecliptic of date) using a long-term
    precession model.

    Parameters
    ----------
    epj : float
        Julian epoch (TT)    
    dr : float
        ICRS right ascension (radians)    
    dd : float
        ICRS declination (radians)

    Returns
    -------
    dl : float
        ecliptic longitude (radians)    
    db : float
        ecliptic latitude (radians)

    '''
    #球面到笛卡尔
    V1=pymS2c(DR,DD)

    #旋转矩阵, ICRS 赤道到黄道
    RM=pymLtecm(EPJ)

    #从ICRS赤道转换到黄道.
    V2=pymRxp(RM,V1)

    #笛卡尔到球面
    A,B=pymC2s(V2)

    #规范化表达
    DL=pymAnp(A)
    DB=pymAnpm(B)
    
    return(DL,DB)


def pymLtp(EPJ):
    '''
    Long-term precession matrix.

    Parameters
    ----------
    epj : float
        Julian epoch (TT)     

    Returns
    -------
    rp : list(3,3)
        precession matrix, J2000.0 to date

    '''
    RP=[[0 for i in range(3)] for j in range(3)]
    
    #赤道极
    PEQR=pymLtpequ(EPJ)
    
    #黄道极
    PECL=pymLtpecl(EPJ)
    
    #春分点
    V=pymPxp(PEQR,PECL)
    W,EQX=pymPn(V)
    
    #矩阵中间行
    V=pymPxp(PEQR,EQX)
    
    #组建矩阵
    for i in range(3):
        RP[0][i]=EQX[i]
        RP[1][i]=V[i]
        RP[2][i]=PEQR[i]
    
    return(RP)


def pymLtpb(EPJ):
    '''
    Long-term precession matrix, including ICRS frame bias.

    Parameters
    ----------
    epj : float
        Julian epoch (TT)    

    Returns
    -------
    rpb : list(3,3)
        precession-bias matrix, J2000.0 to date

    '''
    RPB=[[0 for i in range(3)] for j in range(3)]
    #1角秒对应的弧度
    DAS2R=4.848136811095359935899141e-6
    
    #参考架偏差(IERS Conventions 2010, Eqs. 5.21 and 5.33)
    DX=-0.016617*DAS2R
    DE=-0.0068192*DAS2R
    DR=-0.0146*DAS2R
    
    #岁差矩阵
    RP=pymLtp(EPJ)

    #应用参考架偏差
    for i in range(3):
        RPB[i][0] =  RP[i][0]    - RP[i][1]*DR + RP[i][2]*DX
        RPB[i][1] =  RP[i][0]*DR + RP[i][1]    + RP[i][2]*DE
        RPB[i][2] = -RP[i][0]*DX - RP[i][1]*DE + RP[i][2]
    
    return(RPB)


def pymMoon98(DATE1,DATE2):
    '''
    Approximate geocentric position and velocity of the Moon.

    Parameters
    ----------
    date1 : float
        TT date part A    
    date2 : float
        TT date part A

    Returns
    -------
    pv : list(2,3)
        Moon p,v, GCRS (AU, AU/d)

    '''
    #天文单位(m)
    DAU=149597870.7e3
    
    #参考历元
    DJ00=2451545.0
    
    #1度对应的弧度
    DD2R=1.745329251994329576923691e-2
    
    #1儒略世纪对应的天数
    DJC=36525.0
    
    #基本参数的系数
    #时间在儒略世纪数上的次方，单位：度
    
    #月亮的平黄经(基于给定日期的平春分点和倾角)
    ELP0=218.31665436
    ELP1=481267.88123421
    ELP2=-0.0015786
    ELP3=1.0/538841.0
    ELP4=-1.0/65194000.0
    
    #月亮的平近点角
    EMP0=134.9633964
    EMP1=477198.8675055
    EMP2=0.0087414
    EMP3=1.0/69699.0
    EMP4=-1.0/14712000.0
    
    #月亮的平距角
    D0=297.8501921
    D1=445267.1114034
    D2=-0.0018819
    D3=1.0/545868
    D4=1.0/113065000.0
    
    #太阳的平近点角
    EM0=357.5291092
    EM1=35999.0502909
    EM2=-0.0001536
    EM3=1/2449000
    EM4=0.0
    
    #月亮到其升交点的平均距离
    F0=93.2720950
    F1=483202.0175233
    F2=-0.0036539
    F3=1.0/3526000.0
    F4=1.0/863310000 
    
    #其他参数
    
    #Meeus A_1, 金星相关 (deg)
    A10=119.75
    A11=131.849
    
    #Meeus A_2, 木星相关 (deg)
    A20=53.09
    A21=479264.290
    
    #Meeus A_3, due to sidereal motion of the Moon in longitude (deg)
    A30=313.45
    A31=481266.484
    
    #Coefficients for Meeus "additive terms" (deg)
    AL1=0.003958
    AL2=0.001962
    AL3=0.000318
    
    AB1=-0.002235
    AB2=0.000382
    AB3=0.000175
    AB4=0.000175
    AB5=0.000127
    AB6=-0.000115
    
    #Fixed term in distance (m)
    R0=385000560.0
    
    #Coefficients for (dimensionless) E factor
    E1=-0.002516
    E2=-0.0000074
    
    # *  Coefficients for Moon position series (L,B,R)
    # *
    # *   TLR(1,N)      =  coefficient of L sine term (deg)
    # *   TLR(2,N)      =  coefficient of R cosine term (m)
    # *   TB(N)         =  coefficient B sine term (deg)
    # *   ITx(1-4,N)    =  coefficients of D, M, M', F in argument
    NLR=60
    NB=60
    
    #Longitude and distance series
    TLR=[[6.288774e0,-20905355e0],
            [1.274027e0,-3699111e0],
            [0.658314e0,-2955968e0],
            [0.213618e0,-569925e0],
            [-0.185116e0,48888e0],
            [-0.114332e0,-3149e0],
            [0.058793e0,246158e0],
            [0.057066e0,-152138e0],
            [0.053322e0,-170733e0],
            [0.045758e0,-204586e0],
            [-0.040923e0,-129620e0],
            [-0.034720e0,108743e0],
            [-0.030383e0,104755e0],
            [0.015327e0,10321e0],
            [-0.012528e0,0e0],
            [0.010980e0,79661e0],
            [0.010675e0,-34782e0],
            [0.010034e0,-23210e0],
            [0.008548e0,-21636e0],
            [-0.007888e0,24208e0],
            [-0.006766e0,30824e0],
            [-0.005163e0,-8379e0],
            [0.004987e0,-16675e0],
            [0.004036e0,-12831e0],
            [0.003994e0,-10445e0],
            [0.003861e0,-11650e0],
            [0.003665e0,14403e0],
            [-0.002689e0,-7003e0],
            [-0.002602e0,0e0],
            [0.002390e0,10056e0],
            [-0.002348e0,6322e0],
            [0.002236e0,-9884e0],
            [-0.002120e0,5751e0],
            [-0.002069e0,0e0],
            [0.002048e0,-4950e0],
            [-0.001773e0,4130e0],
            [-0.001595e0,0e0],
            [0.001215e0,-3958e0],
            [-0.001110e0,0e0],
            [-0.000892e0,3258e0],
            [-0.000810e0,2616e0],
            [0.000759e0,-1897e0],
            [-0.000713e0,-2117e0],
            [-0.000700e0,2354e0],
            [0.000691e0,0e0],
            [0.000596e0,0e0],
            [0.000549e0,-1423e0],
            [0.000537e0,-1117e0],
            [0.000520e0,-1571e0],
            [-0.000487e0,-1739e0],
            [-0.000399e0,0e0],
            [-0.000381e0,-4421e0],
            [0.000351e0,0e0],
            [-0.000340e0,0e0],
            [0.000330e0,0e0],
            [0.000327e0,0e0],
            [-0.000323e0,1165e0],
            [0.000299e0,0e0],
            [0.000294e0,0e0],
            [0.000000e0,8752e0]]
    
    #D   M   M'  F
    ITLR=[[0,0,1,0],
            [2,0,-1,0],
            [2,0,0,0],
            [0,0,2,0],
            [0,1,0,0],
            [0,0,0,2],
            [2,0,-2,0],
            [2,-1,-1,0],
            [2,0,1,0],
            [2,-1,0,0],
            [0,1,-1,0],
            [1,0,0,0],
            [0,1,1,0],
            [2,0,0,-2],
            [0,0,1,2],
            [0,0,1,-2],
            [4,0,-1,0],
            [0,0,3,0],
            [4,0,-2,0],
            [2,1,-1,0],
            [2,1,0,0],
            [1,0,-1,0],
            [1,1,0,0],
            [2,-1,1,0],
            [2,0,2,0],
            [4,0,0,0],
            [2,0,-3,0],
            [0,1,-2,0],
            [2,0,-1,2],
            [2,-1,-2,0],
            [1,0,1,0],
            [2,-2,0,0],
            [0,1,2,0],
            [0,2,0,0],
            [2,-2,-1,0],
            [2,0,1,-2],
            [2,0,0,2],
            [4,-1,-1,0],
            [0,0,2,2],
            [3,0,-1,0],
            [2,1,1,0],
            [4,-1,-2,0],
            [0,2,-1,0],
            [2,2,-1,0],
            [2,1,-2,0],
            [2,-1,0,-2],
            [4,0,1,0],
            [0,0,4,0],
            [4,-1,0,0],
            [1,0,-2,0],
            [2,1,0,-2],
            [0,0,2,-2],
            [1,1,1,0],
            [3,0,-2,0],
            [4,0,-3,0],
            [2,-1,2,0],
            [0,2,1,0],
            [1,1,-1,0],
            [2,0,3,0],
            [2,0,-1,-2]]
    
    #Latitude series
    TB=[5.128122e0,
            0.280602e0,
            0.277693e0,
            0.173237e0,
            0.055413e0,
            0.046271e0,
            0.032573e0,
            0.017198e0,
            0.009266e0,
            0.008822e0,
            0.008216e0,
            0.004324e0,
            0.004200e0,
            -0.003359e0,
            0.002463e0,
            0.002211e0,
            0.002065e0,
            -0.001870e0,
            0.001828e0,
            -0.001794e0,
            -0.001749e0,
            -0.001565e0,
            -0.001491e0,
            -0.001475e0,
            -0.001410e0,
            -0.001344e0,
            -0.001335e0,
            0.001107e0,
            0.001021e0,
            0.000833e0,
            0.000777e0,
            0.000671e0,
            0.000607e0,
            0.000596e0,
            0.000491e0,
            -0.000451e0,
            0.000439e0,
            0.000422e0,
            0.000421e0,
            -0.000366e0,
            -0.000351e0,
            0.000331e0,
            0.000315e0,
            0.000302e0,
            -0.000283e0,
            -0.000229e0,
            0.000223e0,
            0.000223e0,
            -0.000220e0,
            -0.000220e0,
            -0.000185e0,
            0.000181e0,
            -0.000177e0,
            0.000176e0,
            0.000166e0,
            -0.000164e0,
            0.000132e0,
            -0.000119e0,
            0.000115e0,
            0.000107e0]
    
    #D   M   M'  F
    ITB=[[0,0,0,1],
            [0,0,1,1],
            [0,0,1,-1],
            [2,0,0,-1],
            [2,0,-1,1],
            [2,0,-1,-1],
            [2,0,0,1],
            [0,0,2,1],
            [2,0,1,-1],
            [0,0,2,-1],
            [2,-1,0,-1],
            [2,0,-2,-1],
            [2,0,1,1],
            [2,1,0,-1],
            [2,-1,-1,1],
            [2,-1,0,1],
            [2,-1,-1,-1],
            [0,1,-1,-1],
            [4,0,-1,-1],
            [0,1,0,1],
            [0,0,0,3],
            [0,1,-1,1],
            [1,0,0,1],
            [0,1,1,1],
            [0,1,1,-1],
            [0,1,0,-1],
            [1,0,0,-1],
            [0,0,3,1],
            [4,0,0,-1],
            [4,0,-1,1],
            [0,0,1,-3],
            [4,0,-2,1],
            [2,0,0,-3],
            [2,0,2,-1],
            [2,-1,1,-1],
            [2,0,-2,1],
            [0,0,3,-1],
            [2,0,2,1],
            [2,0,-3,-1],
            [2,1,-1,1],
            [2,1,0,1],
            [4,0,0,1],
            [2,-1,1,1],
            [2,-2,0,-1],
            [0,0,1,3],
            [2,1,1,-1],
            [1,1,0,-1],
            [1,1,0,1],
            [0,1,-2,-1],
            [2,1,-1,-1],
            [1,0,1,1],
            [2,-1,-2,-1],
            [0,1,2,1],
            [4,0,-2,-1],
            [4,-1,-1,-1],
            [1,0,1,-1],
            [4,0,1,-1],
            [1,0,-1,-1],
            [4,-1,0,-1],
            [2,-2,0,1]]
    
    #自J2000其的儒略世纪数.
    T=((DATE1-DJ00)+DATE2)/DJC
    
    #基本参数
    
    #给定日期的参数(弧度)和导数(每儒略世纪的弧度)

    #月亮的平黄经
    A=ELP0+(ELP1+(ELP2+(ELP3+ELP4*T)*T)*T)*T
    ELP=DD2R*(A%360.0)
    DELP=DD2R*(ELP1+(ELP2*2.0+(ELP3*3.0+ELP4*4.0*T)*T)*T)
    
    #月亮的平距角
    A=D0+(D1+(D2+(D3+D4*T)*T)*T)*T
    D=DD2R*(A%360.0)
    DD=DD2R*(D1+(D2*2.0+(D3*3.0+D4*4.0*T)*T)*T)

    #太阳的平近点角
    A=EM0+(EM1+(EM2+(EM3+EM4*T)*T)*T)*T
    EM=DD2R*(A%360.0)
    DEM=DD2R*(EM1+(EM2*2.0+(EM3*3.0+EM4*4.0*T)*T)*T)

    #月亮的平近点角.
    A=EMP0+(EMP1+(EMP2+(EMP3+EMP4*T)*T)*T)*T
    EMP=DD2R*(A%360.0)
    DEMP=DD2R*(EMP1+(EMP2*2.0+(EMP3*3.0+EMP4*4.0*T)*T)*T)
    
    #月亮到其升交点的平均距离
    A=F0+(F1+(F2+(F3+F4*T)*T)*T)*T
    F=DD2R*(A%360.0)
    DF=DD2R*(F1+(F2*2.0+(F3*3.0+F4*4.0*T)*T)*T)

    #Meeus更进一步的参数.
    A1=DD2R*(A10+A11*T)
    DA1=DD2R*AL1
    A2=DD2R*(A20+A21*T)
    DA2=DD2R*A21
    A3=DD2R*(A30+A31*T)
    DA3=DD2R*A31

    #E-因子和平方项
    E=1.0+(E1+E2*T)*T
    DE=E1+2.0*E2*T
    ESQ=E*E
    DESQ=2.0*E*DE

    #Use the Meeus additive terms (deg) to start off the summations.
    ELPMF=ELP-F
    DELPMF=DELP-DF;
    VEL=AL1*ma.sin(A1)+AL2*ma.sin(ELPMF)+AL3*ma.sin(A2)
    VDEL=AL1*ma.cos(A1)*DA1+AL2*ma.cos(ELPMF)*DELPMF+AL3*ma.cos(A2)*DA2

    VR=0.0
    VDR=0.0
    
    A1MF=A1-F
    DA1MF=DA1-DF
    A1PF=A1+F
    DA1PF=DA1+DF
    DLPMP=ELP-EMP
    SLPMP=ELP+EMP
    VB=AB1*ma.sin(ELP)+AB2*ma.sin(A3)+AB3*ma.sin(A1MF)+AB4*ma.sin(A1PF)\
        +AB5*ma.sin(DLPMP)+AB6*ma.sin(SLPMP)
    VDB=AB1*ma.cos(ELP)*DELP+AB2*ma.cos(A3)*DA3+AB3*ma.cos(A1MF)*DA1MF\
        +AB4*ma.cos(A1PF)*DA1PF+AB5*ma.cos(DLPMP)*(DELP-DEMP)\
            +AB6*ma.cos(SLPMP)*(DELP+DEMP)

    #级数展开
    
    #Longitude and distance plus derivatives.
    for N in range(NLR-1,-1,-1):
        DN=float(ITLR[N][0])
        I=ITLR[N][1]
        EMN=float(I)
        EMPN=float(ITLR[N][2])
        FN=float(ITLR[N][3])
        I=np.abs(I)
        if (I==1):
            EN=E
            DEN=DE
        elif (I==2):
            EN=ESQ
            DEN=DESQ
        else:
            EN=1.0
            DEN=0.0
        ARG=DN*D+EMN*EM+EMPN*EMP+FN*F
        DARG=DN*DD+EMN*DEM+EMPN*DEMP+FN*DF
        FARG=ma.sin(ARG)
        V=FARG*EN
        DV=ma.cos(ARG)*DARG*EN+FARG*DEN
        COEFF=TLR[N][0]
        VEL=VEL+COEFF*V
        VDEL=VDEL+COEFF*DV
        FARG=ma.cos(ARG)
        V=FARG*EN
        DV=-ma.sin(ARG)*DARG*EN+FARG*DEN
        COEFF=TLR[N][1]
        VR=VR+COEFF*V
        VDR=VDR+COEFF*DV
    EL=ELP+DD2R*VEL
    DEL=(DELP+DD2R*VDEL)/DJC
    R=(VR+R0)/DAU
    DR=VDR/DAU/DJC

    #Latitude plus derivative.
    for N in range(NB-1,-1,-1):
        DN=float(ITB[N][0])
        I=ITB[N][1]
        EMN=float(I)
        EMPN=float(ITB[N][2])
        FN=float(ITB[N][3])
        I=np.abs(I)
        if (I==1):
            EN=E
            DEN=DE
        elif (I==2):
            EN=ESQ
            DEN=DESQ
        else:
            EN=1.0
            DEN=0.0
        ARG=DN*D+EMN*EM+EMPN*EMP+FN*F
        DARG=DN*DD+EMN*DEM+EMPN*DEMP+FN*DF
        FARG=ma.sin(ARG)
        V=FARG*EN
        DV=ma.cos(ARG)*DARG*EN+FARG*DEN
        COEFF=TB[N]
        VB=VB+COEFF*V
        VDB=VDB+COEFF*DV
    B=VB*DD2R
    DB=VDB*DD2R/DJC
    
    #转换到最终形式
    
    #Longitude, latitude to x, y, z (AU).
    PV=pymS2pv(EL,B,R,DEL,DB,DR)
    
    #IAU 2006 Fukushima-Williams bias+precession angles.
    GAMB,PHIB,PSIB,EPSA=pymPfw06(DATE1,DATE2)
    
    #Mean ecliptic coordinates to GCRS rotation matrix.
    RM=pymIr()
    RM=pymRz(PSIB,RM)
    RM=pymRx(-PHIB,RM)
    RM=pymRz(-GAMB,RM)

    #Rotate the Moon position and velocity into GCRS (Note 7).
    PV=pymRxpv(RM,PV)
    
    return(PV)


def pymNum00a(DATE1,DATE2):
    '''
    Form the matrix of nutation for a given date, IAU 2000A model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    rmatn : list(3,3)
        nutation matrix

    '''
    #调用au_PN00A获得所需的矩阵，忽略其他项
    DPSI,DEPS,EPSA,RB,RP,RBP,RMATN,RBPN=pymPn00a(DATE1, DATE2)
    
    return(RMATN)


def pymNum00b (DATE1,DATE2):
    '''
    Form the matrix of nutation for a given date, IAU 2000B model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    rmatn : list(3,3)
        nutation matrix

    '''
    #调用au_PN00B获得所需的矩阵，忽略其他项
    DPSI,DEPS,EPSA,RB,RP,RBP,RMATN,RBPN=pymPn00b(DATE1,DATE2)
    
    return(RMATN)


def pymNum06a(DATE1,DATE2):
    '''
    Form the matrix of nutation for a given date, IAU 2006/2000A model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    rmatn : list(3,3)
        nutation matrix

    '''
    #平均倾角
    EPS=pymObl06(DATE1,DATE2)

    #章动量
    DP,DE=pymNut06a(DATE1,DATE2)

    #章动矩阵
    RMATN=pymNumat(EPS,DP,DE)
    
    return(RMATN)


def pymNutm80(DATE1,DATE2):
    '''
    Form the matrix of nutation for a given date, IAU 1980 model.

    Parameters
    ----------
    date1 : float
        TDB date      
    date2 : float
        TDB date

    Returns
    -------
    rmatn : list(3,3)
        nutation matrix

    '''
    #章动量以及平均倾角
    DPSI,DEPS=pymNut80(DATE1,DATE2)
    EPSA=pymObl80(DATE1, DATE2)
    
    #构建旋转矩阵
    RMATN=pymNumat(EPSA,DPSI,DEPS)
    
    return(RMATN)


def pymP2pv(P):
    '''
    Extend a p-vector to a pv-vector by appending a zero velocity.

    Parameters
    ----------
    p : list(3)
        p-vector

    Returns
    -------
    pv : list(2,3)
        pv-vector

    ''' 
    PV=[0,0]
    PV[0]=pymCp(P)
    PV[1]=pymZp()
    
    return(PV)


def pymP2s(P):
    '''
    P-vector to spherical polar coordinates.

    Parameters
    ----------
    p : list(3)
        p-vector

    Returns
    -------
    theta : float
        longitude angle (radians)    
    phi : float
        latitude angle (radians)    
    r : float
        radial distance

    '''
    THETA,PHI=pymC2s(P)
    R=pymPm(P)
    
    return(THETA,PHI,R)


def pymP06e(DATE1,DATE2):
    '''
    Precession angles, IAU 2006, equinox based.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    eps0 : float
        epsilon_0    
    psia : float
        psi_A    
    oma : float
        omega_A    
    bpa : float
        P_A    
    bqa : float
        Q_A    
    pia : float
        pi_A    
    bpia : float
        Pi_A    
    epsa : float
        obliquity epsilon_A    
    chia : float
        chi_A     
    za : float
        z_A     
    zetaa : float
        zeta_A    
    thetaa : float 
        theta_A    
    pa : float
        p_A    
    gma : float
        F-W angle gamma_J2000    
    phi : float
        F-W angle phi_J2000    
    psi : float
        F-W angle psi_J2000    

    '''
    #1角秒对应的弧度
    DAS2R=4.848136811095359935899141e-6
    
    #参考历元(J2000.0), JD
    DJ00=2451545.0
    
    #1儒略世纪对应的天数
    DJC=36525.0
    
    #给定日期到参考历元的时间间隔，儒略世纪数
    T=((DATE1-DJ00)+DATE2)/DJC
    
    #倾角J2000.0.
    EPS0=84381.406*DAS2R
    
    #日月岁差
    PSIA=(5038.481507+(-1.0790069+(-0.00114045+(0.000132851+(-0.0000000951)\
                                                 *T)*T)*T)*T)*T*DAS2R
    
    #平赤道相对于J2000.0黄道的倾角。
    OMA=EPS0+(-0.025754+(0.0512623+(-0.00772503+(-0.000000467+(0.0000003337)\
                                                 *T)*T)*T)*T)*T*DAS2R

    #黄极 x, J2000.0 
    BPA=(4.199094+(0.1939873+(-0.00022466+(-0.000000912+(0.0000000120)\
                                           *T)*T)*T)*T)*T*DAS2R

    #黄极 -y, J2000.0 ecliptic triad.
    BQA=(-46.811015+(0.0510283+(0.00052413+(-0.000000646+(-0.0000000172)\
                                            *T)*T)*T)*T)*T*DAS2R

    #变化的黄道与J2000.0黄道之间的夹角
    PIA=(46.998973+(-0.0334926+(-0.00012559+(0.000000113+(-0.0000000022)\
                                             *T)*T)*T)*T)*T*DAS2R

    #变化的黄道的升交点黄经.
    BPIA=(629546.7936+(-867.95758+(0.157992+(-0.0005371+(-0.00004797\
                                        +(0.000000072)*T)*T)*T)*T)*T)*DAS2R

    #黄道的平均倾角
    EPSA=pymObl06(DATE1,DATE2)

    #行星岁差.
    CHIA=(10.556403+(-2.3814292+(-0.00121197+(0.000170663+(-0.0000000560)\
                                              *T)*T)*T)*T)*T*DAS2R

    #赤道岁差: 负的第三个 323欧拉角.
    ZA=(-2.650545+(2306.077181+(1.0927348+(0.01826837+(-0.000028596\
                                       +(-0.0000002904)*T)*T)*T)*T)*T)*DAS2R

    #赤道岁差: 负的第一个 323欧拉角.
    ZETAA=(2.650545+(2306.083227+(0.2988499+(0.01801828+(-0.000005971\
                                     +(-0.0000003173)*T)*T)*T)*T)*T)*DAS2R

    #赤道岁差: 第二个 323欧拉角.
    THETAA=(2004.191903+(-0.4294934+(-0.04182264+(-0.000007089\
                                      +(-0.0000001274)*T)*T)*T)*T)*T*DAS2R

    #累计岁差
    PA=(5028.796195+(1.1054348+(0.00007964+(-0.000023857+(-0.0000000383)\
                                            *T)*T)*T)*T)*T*DAS2R
          

    #Fukushima-Williams岁差角
    GAM=(10.556403+(0.4932044+(-0.00031238+(-0.000002788+(0.0000000260)\
                                           *T)*T)*T)*T)*T*DAS2R

    PHI=EPS0+(-46.811015+(0.0511269+(0.00053289+(-0.000000440\
                                     +(-0.0000000176)*T)*T)*T)*T)*T*DAS2R

    PSI=(5038.481507+(1.5584176+(-0.00018522+(-0.000026452+(-0.0000000148)\
                                              *T)*T)*T)*T)*T*DAS2R

    return(EPS0,PSIA,OMA,BPA,BQA,PIA,BPIA,EPSA,CHIA,
           ZA,ZETAA,THETAA,PA,GAM,PHI,PSI)


def pymPap(A,B):
    '''
    Position-angle from two p-vectors.

    Parameters
    ----------
    a : list(3)
        direction of reference point    
    b : list(3)
        direction of point whose PA is required

    Returns
    -------
    function value : float
        position angle of b with respect to a (radians)

    '''
    ETA=[0,0,0]
    #A向量的模和方向
    AM,AU=pymPn(A)
    
    #B向量的模
    BM=pymPm(B)
    
    #对空向量进行处理
    if (AM==0.0)|(BM==0.0):
        ST=0.0
        CT=1.0
    else:
        
        #自A点到“北”轴切线(任意长度)
        XA=A[0]
        YA=A[1]
        ZA=A[2]
        ETA[0]=-XA*ZA
        ETA[1]=-YA*ZA
        ETA[2]=XA*XA+YA*YA
    
        #自A点到“东”轴切线(任意长度)。
        XI=pymPxp(ETA,AU)
    
        #从A到B的向量
        A2B=pymPmp(B,A)
        
        #沿着北轴和东轴分解参数
        ST=pymPdp(A2B,XI)
        CT=pymPdp(A2B,ETA)
        
        #处理有误情况
        if (ST==0.0)&(CT==0.0):
            CT=1.0
    
    #方位角
    THETA=ma.atan2(ST,CT)
        
    return(THETA)


def pymPas(AL,AP,BL,BP):
    '''
    Position-angle from spherical coordinates. 

    Parameters
    ----------
    al : float
        longitude of point A (e.g. RA) in radians     
    ap : float
        latitude of point A (e.g. Dec) in radians    
    bl : float
        longitude of point B     
    bp : float
        latitude of point B     

    Returns
    -------
    THETA : float
        position angle of B with respect to A

    '''
    DL=BL-AL
    Y=ma.sin(DL)*ma.cos(BP)
    X=ma.sin(BP)*ma.cos(AP)-ma.cos(BP)*ma.sin(AP)*ma.cos(DL)
    if (X!=0.0)|(Y!=0.0):
        THETA=ma.atan2(Y,X)
    else:
        THETA=0.0
    
    return(THETA)


def pymPb06(DATE1,DATE2):
    '''
    This function forms three Euler angles which implement general
    precession from epoch J2000.0, using the IAU 2006 model.  Frame
    bias (the offset between ICRS and mean J2000.0) is included.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    bzeta : float
        1st rotation: radians cw around z    
    ca : float
        3rd rotation: radians cw around z     
    btheta : float
        2nd rotation: radians ccw around y

    '''
    #岁差矩阵基于 Fukushima-Williams angles.
    R=pymPmat06(DATE1,DATE2)
    
    #求解z，选取+/-
    Y=R[1][2]
    X=-R[0][2]
    if (X<0.0):
        Y=-Y
        X=-X
    if (X!=0.0)|(Y!=0.0):
        BZ=-ma.atan2(Y,X)
    else:
        BZ=0.0
        
    #把它从矩阵中划去
    R=pymRz(BZ,R)
    
    #解出剩余两个角
    Y=R[0][2]
    X=R[2][2]
    if (X!=0.0)|(Y!=0.0):
        BTHETA=-ma.atan2(Y,X)
    else:
        BTHETA=0.0
    
    Y=-R[1][0]
    X=R[1][1]
    if (X!=0.0)|(Y!=0.0):
        BZETA=-ma.atan2(Y,X)
    else:
        BZETA=0.0
    
    return(BZETA,BZ,BTHETA)


def pymPlan94(DATE1,DATE2,NP):
    '''
    Approximate heliocentric position and velocity of a nominated major
    planet:  Mercury, Venus, EMB, Mars, Jupiter, Saturn, Uranus or
    Neptune (but not the Earth itself).

    Parameters
    ----------
    date1 : float
        TDB date part A    
    date2 : float
        TDB date part B    
    npp : int
        planet (1=Mercury, 2=Venus, 3=EMB, 4=Mars,
        5=Jupiter, 6=Saturn, 7=Uranus, 8=Neptune)

    Returns
    -------
    pv : list(2,3)
        planet p,v (heliocentric, J2000.0, au,au/d)
    J : ValueError
       -1 = illegal NP (outside 1-8)
        0 = OK
       +1 = warning: date outside 1000-3000 CE
       +2 = warning: solution failed to converge
    '''  
    PV=[[0 for i in range(3)] for j in range(2)]
    
    #解开普勒方程的最大迭代数
    KMAX=10
    
    #2Pi
    D2PI=6.283185307179586476925287
    
    #1角秒对应的弧度
    DAS2R=4.848136811095359935899141e-6

    #参考历元(J2000.0), JD
    DJ00=2451545.0
    
    #1千年的儒略日数
    DJM=365250.0
    
    #J2000.0平倾角的正余弦 (IAU 1976)
    SINEPS=0.3977771559319137
    COSEPS=0.9174820620691818
    
    #高斯常数
    GK=0.017202098950
    
    #行星质量倒数
    AMAS=[6023600,408523.5,328900.5,3098710,1047.355,3498.5,22869.0,19314.0]
    
    # Tables giving the mean Keplerian elements, limited to T**2 terms:
    
    #          A       semi-major axis (au)
    #          DLM     mean longitude (degree and arcsecond)
    #          E       eccentricity
    #          PI      longitude of the perihelion (degree and arcsecond)
    #          DINC    inclination (degree and arcsecond)
    #          OMEGA   longitude of the ascending node (degree and arcsecond)
    
    A=[[0.3870983098,0.0,0.0],
        [0.7233298200,0.0,0.0],
        [1.0000010178,0.0,0.0],
        [1.5236793419,3e-10,0.0],
        [5.2026032092,19132e-10,-39e-10],
        [9.5549091915,-0.0000213896,444e-10],
        [19.2184460618,-3716e-10,979e-10],
        [30.1103868694,-16635e-10,686e-10]]

    DLM=[[252.25090552e0,5381016286.88982e0,-1.92789e0],
        [181.97980085e0,2106641364.33548e0,0.59381e0],
        [100.46645683e0,1295977422.83429e0,-2.04411e0],
        [355.43299958e0,689050774.93988e0,0.94264e0],
        [34.35151874e0,109256603.77991e0,-30.60378e0],
        [50.07744430e0,43996098.55732e0,75.61614e0],
        [314.05500511e0,15424811.93933e0,-1.75083e0],
        [304.34866548e0,7865503.20744e0,0.21103e0]]
            
    E=[[0.2056317526e0,0.0002040653e0,-28349e-10],
        [0.0067719164e0,-0.0004776521e0,98127e-10],
        [0.0167086342e0,-0.0004203654e0,-0.0000126734e0],
        [0.0934006477e0,0.0009048438e0,-80641e-10],
        [0.0484979255e0,0.0016322542e0,-0.0000471366e0],
        [0.0555481426e0,-0.0034664062e0,-0.0000643639e0],
        [0.0463812221e0,-0.0002729293e0,0.0000078913e0],
        [0.0094557470e0,0.0000603263e0,0e0]]
    
    PI=[[77.45611904e0,5719.11590e0,-4.83016e0],
        [131.56370300e0,175.48640e0,-498.48184e0],
        [102.93734808e0,11612.35290e0,53.27577e0],
        [336.06023395e0,15980.45908e0,-62.32800e0],
        [14.33120687e0,7758.75163e0,259.95938e0],
        [93.05723748e0,20395.49439e0,190.25952e0],
        [173.00529106e0,3215.56238e0,-34.09288e0],
        [48.12027554e0,1050.71912e0,27.39717e0]]
    
    DINC=[[7.00498625e0,-214.25629e0,0.28977e0],
        [3.39466189e0,-30.84437e0,-11.67836e0],
        [0e0,469.97289e0,-3.35053e0],
        [1.84972648e0,-293.31722e0,-8.11830e0],
        [1.30326698e0,-71.55890e0,11.95297e0],
        [2.48887878e0,91.85195e0,-17.66225e0],
        [0.77319689e0,-60.72723e0,1.25759e0],
        [1.76995259e0,8.12333e0,0.08135e0]]
    
    OMEGA=[[48.33089304e0,-4515.21727e0,-31.79892e0],
        [76.67992019e0,-10008.48154e0,-51.32614e0],
        [174.87317577e0,-8679.27034e0,15.34191e0],
        [49.55809321e0,-10620.90088e0,-230.57416e0],
        [100.46440702e0,6362.03561e0,326.52178e0],
        [113.66550252e0,-9240.19942e0,-66.23743e0],
        [74.00595701e0,2669.15033e0,145.93964e0],
        [131.78405702e0,-221.94322e0,-0.78728e0]]
    
    # Tables for trigonometric terms to be added to the mean elements
    # of the semi-major axes.
    
    KP=[[69613,75645,88306,59899,15746,71087,142173,3086,0],
        [21863,32794,26934,10931,26250,43725,53867,28939,0],
        [16002,21863,32004,10931,14529,16368,15318,32794,0],
        [6345,7818,15636,7077,8184,14163,1107,4872,0],
        [1760,1454,1167,880,287,2640,19,2047,1454],
        [574,0,880,287,19,1760,1167,306,574],
        [204,0,177,1265,4,385,200,208,204],
        [0,102,106,4,98,1367,487,204,0]]
    
    CA=[[4,-13,11,-9,-9,-3,-1,4,0],
        [-156,59,-42,6,19,-20,-10,-12,0],
        [64,-152,62,-8,32,-41,19,-11,0],
        [124,621,-145,208,54,-57,30,15,0],
        [-23437,-2634,6601,6259,-1507,-1821,2620,-2115,-1489],
        [62911,-119919,79336,17814,-24241,12068,8306,-4893,8902],
        [389061,-262125,-44088,8387,-22976,-2093,-615,-9720,6633],
        [-412235,-157046,-31430,37817,-9740,-13,-7449,9644,0]]
    
    SA=[[-29,-1,9,6,-6,5,4,0,0],
        [-48,-125,-26,-37,18,-13,-20,-2,0],
        [-150,-46,68,54,14,24,-28,22,0],
        [-621,532,-694,-20,192,-94,71,-73,0],
        [-14614,-19828,-5869,1881,-4372,-2255,782,930,913],
        [139737,0,24667,51123,-5102,7429,-4095,-1976,-9566],
        [-138081,0,37205,-49039,-41901,-33872,-27037,-12474,18797],
        [0,28492,133236,69654,52322,-49577,-26430,-3593,0]]
    
    # Tables giving the trigonometric terms to be added to the mean
    # elements of the mean longitudes.
    
    KQ=[[3086,15746,69613,59899,75645,88306,12661,2658,0,0],
        [21863,32794,10931,73,4387,26934,1473,2157,0,0],
        [10,16002,21863,10931,1473,32004,4387,73,0,0],
        [10,6345,7818,1107,15636,7077,8184,532,10,0],
        [19,1760,1454,287,1167,880,574,2640,19,1454],
        [19,574,287,306,1760,12,31,38,19,574],
        [4,204,177,8,31,200,1265,102,4,204],
        [4,102,106,8,98,1367,487,204,4,102]]
    
    CL=[[21,-95,-157,41,-5,42,23,30,0,0],
        [-160,-313,-235,60,-74,-76,-27,34,0,0],
        [-325,-322,-79,232,-52,97,55,-41,0,0],
        [2268,-979,802,602,-668,-33,345,201,-55,0],
        [7610,-4997,-7689,-5841,-2617,1115,-748,-607,6074,354],
        [-18549,30125,20012,-730,824,23,1289,-352,-14767,-2062],
        [-135245,-14594,4197,-4030,-5630,-2898,2540,-306,2939,1986],
        [89948,2103,8963,2695,3682,1648,866,-154,-1963,-283]]
    
    SL=[[-342,136,-23,62,66,-52,-33,17,0,0],
        [524,-149,-35,117,151,122,-71,-62,0,0],
        [-105,-137,258,35,-116,-88,-112,-80,0,0],
        [854,-205,-936,-240,140,-341,-97,-232,536,0],
        [-56980,8016,1012,1448,-3024,-3710,318,503,3767,577],
        [138606,-13478,-4964,1441,-1319,-1482,427,1236,-9167,-1918],
        [71234,-41116,5334,-4935,-1848,66,434,-1748,3780,-701],
        [-47645,11647,2166,3194,679,0,-244,-419,-2531,48]]
    
    #验证行星编号
    if (NP<1)|(NP>8):
        JSTAT=-1
        
        #验证行星编号
        for k in range(2):
            for i in range(3):
                PV[k][i]=0.0
    
    else:
        
        #时间：自 J2000.0起，千儒略年.
        T=((DATE1-DJ00)+DATE2)/DJM
        
        #判断日期是否超出限制.
        if (np.abs(T)<=1.0):
            JSTAT=0
        else:
            JSTAT=1
        
        #计算平均参数
        DA=A[NP-1][0]+(A[NP-1][1]+A[NP-1][2]*T)*T
        DL=(3600.0*DLM[NP-1][0]+(DLM[NP-1][1]+DLM[NP-1][2]*T)*T)*DAS2R
        DE=E[NP-1][0]+(E[NP-1][1]+E[NP-1][2]*T)*T
        DP=pymAnpm((3600.0*PI[NP-1][0]+(PI[NP-1][1]+PI[NP-1][2]*T)*T)*DAS2R)
        DI=(3600.0*DINC[NP-1][0]+(DINC[NP-1][1]+DINC[NP-1][2]*T)*T)*DAS2R
        DOM=pymAnpm((3600.0*OMEGA[NP-1][0]+(OMEGA[NP-1][1]+OMEGA[NP-1][2]\
                                             *T)*T)*DAS2R)
        
        #应用三角函数项
        DMU=0.35953620*T
        for K in range(8):
            ARGA=KP[NP-1][K]*DMU
            ARGL=KQ[NP-1][K]*DMU
            DA=DA+(CA[NP-1][K]*ma.cos(ARGA)+SA[NP-1][K]*ma.sin(ARGA))*1e-7
            DL=DL+(CL[NP-1][K]*ma.cos(ARGL)+SL[NP-1][K]*ma.sin(ARGL))*1e-7
            
        ARGA=KP[NP-1][8]*DMU
        DA=DA+T*(CA[NP-1][8]*ma.cos(ARGA)+SA[NP-1][8]*ma.sin(ARGA))*1e-7
        
        for K in range(8,10,1):
            ARGL=KQ[NP-1][K]*DMU
            DL=DL+T*(CL[NP-1][K]*ma.cos(ARGL)+SL[NP-1][K]*ma.sin(ARGL))*1e-7
        C=np.abs(DL)%D2PI
        if DL<0:
            B=-C
        else:
            B=C
        DL=B
        
        #迭代求解开普勒方程获得偏近点角
        AM=DL-DP
        AE=AM+DE*ma.sin(AM)
        K=0
        DAE=1
        while (K<KMAX)&(np.abs(DAE)>1e-12):
            DAE=(AM-AE+DE*ma.sin(AE))/(1.0-DE*ma.cos(AE))
            AE=AE+DAE
            K=K+1
            if (K>=KMAX):
                JSTAT=2
            
        #真实近点角.
        AE2=AE/2.0
        AT=2.0*ma.atan2(ma.sqrt((1.0+DE)/(1.0-DE))*ma.sin(AE2),ma.cos(AE2))
        
        #距离(au) 和速度(弧度/day).
        R=DA*(1.0-DE*ma.cos(AE))
        V=GK*ma.sqrt((1.0+1.0/AMAS[NP])/(DA*DA*DA))
        
        SI2=ma.sin(DI/2.0)
        XQ=SI2*ma.cos(DOM)
        XP=SI2*ma.sin(DOM)
        TL=AT+DP
        XSW=ma.sin(TL)
        XCW=ma.cos(TL)
        XM2=2.0*(XP*XCW-XQ*XSW)
        XF=DA/ma.sqrt(1.0-DE*DE)
        CI2=ma.cos(DI/2.0)
        XMS=(DE*ma.sin(DP)+XSW)*XF
        XMC=(DE*ma.cos(DP)+XCW)*XF
        XPXQ2=2.0*XP*XQ

        #位置(J2000.0 ecliptic x,y,z in au).
        X=R*(XCW-XM2*XP)
        Y=R*(XSW+XM2*XQ)
        Z=R*(-XM2*CI2)
        
        #旋转到赤道
        PV[0][0]=X
        PV[0][1]=Y*COSEPS-Z*SINEPS
        PV[0][2]=Y*SINEPS+Z*COSEPS
        
        #速度 (J2000.0 ecliptic xdot,ydot,zdot in au/d).
        X=V*((-1.0+2.0*XP*XP)*XMS+XPXQ2*XMC)
        Y=V*((1.0-2.0*XQ*XQ)*XMC-XPXQ2*XMS)
        Z=V*(2.0*CI2*(XP*XMS+XQ*XMC))

        #旋转到赤道
        PV[1][0]=X
        PV[1][1]=Y*COSEPS-Z*SINEPS
        PV[1][2]=Y*SINEPS+Z*COSEPS
        
    J=JSTAT

    return(PV,J)


def pymXy06(DATE1,DATE2):
    '''
    X,Y coordinates of celestial intermediate pole from series based
    on IAU 2006 precession and IAU 2000A nutation.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date     
    date2 : float
        TT as a 2-part Julian Date     

    Returns
    -------
    x : float
        CIP X,Y coordinates    
    y : float
        CIP X,Y coordinates

    '''
    #1角秒对应的弧度
    DAS2R=4.848136811095359935899141e-6

    #参考历元(J2000.0), JD
    DJ00=2451545.0

    #1儒略世纪对应的天数
    DJC=36525.0
    
    #T在X和Y的多项式中的最大次幂
    MAXPT=5
    
    #频率数:日月、行星、总计
    NFLS=653
    NFPL=656
    NF=NFLS+NFPL
    
    #振幅系数数
    NA=4755
    
    PT=[0 for i in range(MAXPT+1)]
    FA=[0 for i in range(14)]
    XYPR=[0,0]
    XYPL=[0,0]
    XYLS=[0,0]
    SC=[0,0]
    
    #多项式系数(角秒).
    XYP=[[-0.016617,+2004.191898,-0.4297829,-0.19861834,
          +0.000007578,+0.0000059285],
         [-0.006951,-0.025896,-22.4072747,+0.00190059,
          +0.001112526,+0.0000001358]]
    
    #基本参数:日月项.
    MFALS=[[0,0,0,0,1],
            [0,0,2,-2,2],
            [0,0,2,0,2],
            [0,0,0,0,2],
            [0,1,0,0,0],
            [0,1,2,-2,2],
            [1,0,0,0,0],
            [0,0,2,0,1],
            [1,0,2,0,2],
            [0,1,-2,2,-2],
            [0,0,2,-2,1],
            [1,0,-2,0,-2],
            [1,0,0,-2,0],
            [1,0,0,0,1],
            [1,0,0,0,-1],
            [1,0,-2,-2,-2],
            [1,0,2,0,1],
            [2,0,-2,0,-1],
            [0,0,0,2,0],
            [0,0,2,2,2],
            [2,0,0,-2,0],
            [0,2,-2,2,-2],
            [2,0,2,0,2],
            [1,0,2,-2,2],
            [1,0,-2,0,-1],
            [2,0,0,0,0],
            [0,0,2,0,0],
            [0,1,0,0,1],
            [1,0,0,-2,-1],
            [0,2,2,-2,2],
            [0,0,2,-2,0],
            [1,0,0,-2,1],
            [0,1,0,0,-1],
            [0,2,0,0,0],
            [1,0,-2,-2,-1],
            [1,0,2,2,2],
            [0,1,2,0,2],
            [2,0,-2,0,0],
            [0,0,2,2,1],
            [0,1,-2,0,-2],
            [0,0,0,2,1],
            [1,0,2,-2,1],
            [2,0,0,-2,-1],
            [2,0,2,-2,2],
            [2,0,2,0,1],
            [0,0,0,2,-1],
            [0,1,-2,2,-1],
            [1,1,0,-2,0],
            [2,0,0,-2,1],
            [1,0,0,2,0],
            [0,1,2,-2,1],
            [1,-1,0,0,0],
            [0,1,-1,1,-1],
            [2,0,-2,0,-2],
            [0,1,0,-2,0],
            [1,0,0,-1,0],
            [3,0,2,0,2],
            [0,0,0,1,0],
            [1,-1,2,0,2],
            [1,1,-2,-2,-2],
            [1,0,-2,0,0],
            [2,0,0,0,-1],
            [0,1,-2,-2,-2],
            [1,1,2,0,2],
            [2,0,0,0,1],
            [1,1,0,0,0],
            [1,0,-2,2,-1],
            [1,0,2,0,0],
            [1,-1,0,-1,0],
            [1,0,0,0,2],
            [1,0,-1,0,-1],
            [0,0,2,1,2],
            [1,0,-2,-4,-2],
            [1,-1,0,-1,-1],
            [1,0,2,2,1],
            [0,2,-2,2,-1],
            [1,0,0,0,-2],
            [2,0,-2,-2,-2],
            [1,1,2,-2,2],
            [2,0,-2,-4,-2],
            [1,0,-4,0,-2],
            [2,0,2,-2,1],
            [1,0,0,-1,-1],
            [2,0,2,2,2],
            [3,0,0,0,0],
            [1,0,0,2,1],
            [0,0,2,-2,-1],
            [3,0,2,-2,2],
            [0,0,4,-2,2],
            [1,0,0,-4,0],
            [0,1,2,0,1],
            [2,0,0,-4,0],
            [1,1,0,-2,-1],
            [2,0,-2,0,1],
            [0,0,2,0,-1],
            [0,1,-2,0,-1],
            [0,1,0,0,2],
            [0,0,2,-1,2],
            [0,0,2,4,2],
            [2,1,0,-2,0],
            [1,1,0,-2,1],
            [1,-1,0,-2,0],
            [1,-1,0,-1,-2],
            [1,-1,0,0,1],
            [0,1,-2,2,0],
            [0,1,0,0,-2],
            [1,-1,2,2,2],
            [1,0,0,2,-1],
            [1,-1,-2,-2,-2],
            [3,0,2,0,1],
            [0,1,2,2,2],
            [1,0,2,-2,0],
            [1,1,-2,-2,-1],
            [1,0,2,-4,1],
            [0,1,-2,-2,-1],
            [2,-1,2,0,2],
            [0,0,0,2,2],
            [1,-1,2,0,1],
            [1,-1,-2,0,-2],
            [0,1,0,2,0],
            [0,1,2,-2,0],
            [0,0,0,1,1],
            [1,0,-2,-2,0],
            [0,3,2,-2,2],
            [2,1,2,0,2],
            [1,1,0,0,1],
            [2,0,0,2,0],
            [1,1,2,0,1],
            [1,0,0,-2,-2],
            [1,0,-2,2,0],
            [1,0,-1,0,-2],
            [0,1,0,-2,1],
            [0,1,0,1,0],
            [0,0,0,1,-1],
            [1,0,-2,2,-2],
            [1,-1,0,0,-1],
            [0,0,0,4,0],
            [1,-1,0,2,0],
            [1,0,2,1,2],
            [1,0,2,-1,2],
            [0,0,2,1,1],
            [1,0,0,-2,2],
            [1,0,-2,0,1],
            [1,0,-2,-4,-1],
            [0,0,2,2,0],
            [1,1,2,-2,1],
            [1,0,-2,1,-1],
            [0,0,1,0,1],
            [2,0,-2,-2,-1],
            [4,0,2,0,2],
            [2,-1,0,0,0],
            [2,1,2,-2,2],
            [0,1,2,1,2],
            [1,0,4,-2,2],
            [1,1,0,0,-1],
            [2,0,2,0,0],
            [2,0,-2,-4,-1],
            [1,0,-1,0,0],
            [1,0,0,1,0],
            [0,1,0,2,1],
            [1,0,-4,0,-1],
            [1,0,0,-4,-1],
            [2,0,2,2,1],
            [2,1,0,0,0],
            [0,0,2,-3,2],
            [1,2,0,-2,0],
            [0,3,0,0,0],
            [0,0,4,0,2],
            [0,0,2,-4,1],
            [2,0,0,-2,-2],
            [1,1,-2,-4,-2],
            [0,1,0,-2,-1],
            [0,0,0,4,1],
            [3,0,2,-2,1],
            [1,0,2,4,2],
            [1,1,-2,0,-2],
            [0,0,4,-2,1],
            [2,-2,0,-2,0],
            [2,1,0,-2,-1],
            [0,2,0,-2,0],
            [1,0,0,-1,1],
            [1,1,2,2,2],
            [3,0,0,0,-1],
            [2,0,0,-4,-1],
            [3,0,2,2,2],
            [0,0,2,4,1],
            [0,2,-2,-2,-2],
            [1,-1,0,-2,-1],
            [0,0,2,-1,1],
            [2,0,0,2,1],
            [1,-1,-2,2,-1],
            [0,0,0,2,-2],
            [2,0,0,-4,1],
            [1,0,0,-4,1],
            [2,0,2,-4,1],
            [4,0,2,-2,2],
            [2,1,-2,0,-1],
            [2,1,-2,-4,-2],
            [3,0,0,-4,0],
            [1,-1,2,2,1],
            [1,-1,-2,0,-1],
            [0,2,0,0,1],
            [1,2,-2,-2,-2],
            [1,1,0,-4,0],
            [2,0,0,-2,2],
            [0,2,2,-2,1],
            [1,0,2,0,-1],
            [2,1,0,-2,1],
            [2,-1,-2,0,-1],
            [1,-1,-2,-2,-1],
            [0,1,-2,1,-2],
            [1,0,-4,2,-2],
            [0,1,2,2,1],
            [3,0,0,0,1],
            [2,-1,2,2,2],
            [0,1,-2,-4,-2],
            [1,0,-2,-3,-2],
            [2,0,0,0,2],
            [1,-1,0,-2,-2],
            [2,0,-2,2,-1],
            [0,2,-2,0,-2],
            [3,0,-2,0,-1],
            [2,-1,2,0,1],
            [1,0,-2,-1,-2],
            [0,0,2,0,3],
            [2,0,-4,0,-2],
            [2,1,0,-4,0],
            [1,1,-2,1,-1],
            [0,2,2,0,2],
            [1,-1,2,-2,2],
            [1,-1,0,-2,1],
            [2,1,2,0,1],
            [1,0,2,-4,2],
            [1,1,-2,0,-1],
            [1,1,0,2,0],
            [1,0,0,-3,0],
            [2,0,2,-1,2],
            [0,2,0,0,-1],
            [2,-1,0,-2,0],
            [4,0,0,0,0],
            [2,1,-2,-2,-2],
            [0,2,-2,2,0],
            [1,0,2,1,1],
            [1,0,-1,0,-3],
            [3,-1,2,0,2],
            [2,0,2,-2,0],
            [1,-2,0,0,0],
            [2,0,0,0,-2],
            [1,0,0,4,0],
            [0,1,0,1,1],
            [1,0,2,2,0],
            [0,1,0,2,-1],
            [0,1,0,1,-1],
            [0,0,2,-2,3],
            [3,1,2,0,2],
            [1,1,2,1,2],
            [1,1,-2,2,-1],
            [2,-1,2,-2,2],
            [1,-2,2,0,2],
            [1,0,2,-4,0],
            [0,0,1,0,0],
            [1,0,2,-3,1],
            [1,-2,0,-2,0],
            [2,0,0,2,-1],
            [1,1,2,-4,1],
            [4,0,2,0,1],
            [0,1,2,1,1],
            [1,2,2,-2,2],
            [2,0,2,1,2],
            [2,1,2,-2,1],
            [1,0,2,-1,1],
            [1,0,4,-2,1],
            [1,-1,2,-2,1],
            [0,1,0,-4,0],
            [3,0,-2,-2,-2],
            [0,0,4,-4,2],
            [2,0,-4,-2,-2],
            [2,-2,0,-2,-1],
            [1,0,2,-2,-1],
            [2,0,-2,-6,-2],
            [1,0,-2,1,-2],
            [1,0,-2,2,1],
            [1,-1,0,2,-1],
            [1,0,-2,1,0],
            [2,-1,0,-2,1],
            [1,-1,0,2,1],
            [2,0,-2,-2,0],
            [1,0,2,-3,2],
            [0,0,0,4,-1],
            [2,-1,0,0,1],
            [2,0,4,-2,2],
            [0,0,2,3,2],
            [0,1,4,-2,2],
            [0,1,-2,2,1],
            [1,1,0,2,1],
            [1,0,0,4,1],
            [0,0,4,0,1],
            [2,0,0,-3,0],
            [1,0,0,-1,-2],
            [1,-2,-2,-2,-2],
            [3,0,0,2,0],
            [2,0,2,-4,2],
            [1,1,-2,-4,-1],
            [1,0,-2,-6,-2],
            [2,-1,0,0,-1],
            [2,-1,0,2,0],
            [0,1,2,-2,-1],
            [1,1,0,1,0],
            [1,2,0,-2,-1],
            [1,0,0,1,-1],
            [0,0,1,0,2],
            [3,1,2,-2,2],
            [1,0,-4,-2,-2],
            [1,0,2,4,1],
            [1,-2,2,2,2],
            [1,-1,-2,-4,-2],
            [0,0,2,-4,2],
            [0,0,2,-3,1],
            [2,1,-2,0,0],
            [3,0,-2,-2,-1],
            [2,0,2,4,2],
            [0,0,0,0,3],
            [2,-1,-2,-2,-2],
            [2,0,0,-1,0],
            [3,0,2,-4,2],
            [2,1,2,2,2],
            [0,0,3,0,3],
            [1,1,2,2,1],
            [2,1,0,0,-1],
            [1,2,0,-2,1],
            [3,0,2,2,1],
            [1,-1,-2,2,-2],
            [1,1,0,-1,0],
            [1,2,0,0,0],
            [1,0,4,0,2],
            [1,-1,2,4,2],
            [2,1,0,0,1],
            [1,0,0,2,2],
            [1,-1,-2,2,0],
            [0,2,-2,-2,-1],
            [2,0,-2,0,2],
            [5,0,2,0,2],
            [3,0,-2,-6,-2],
            [1,-1,2,-1,2],
            [3,0,0,-4,-1],
            [1,0,0,1,1],
            [1,0,-4,2,-1],
            [0,1,2,-4,1],
            [1,2,2,0,2],
            [0,1,0,-2,-2],
            [0,0,2,-1,0],
            [1,0,1,0,1],
            [0,2,0,-2,1],
            [3,0,2,0,0],
            [1,1,-2,1,0],
            [2,1,-2,-4,-1],
            [3,-1,0,0,0],
            [2,-1,-2,0,0],
            [4,0,2,-2,1],
            [2,0,-2,2,0],
            [1,1,2,-2,0],
            [1,0,-2,4,-1],
            [1,0,-2,-2,1],
            [2,0,2,-4,0],
            [1,1,0,-2,-2],
            [1,1,-2,-2,0],
            [1,0,1,-2,1],
            [2,-1,-2,-4,-2],
            [3,0,-2,0,-2],
            [0,1,-2,-2,0],
            [3,0,0,-2,-1],
            [1,0,-2,-3,-1],
            [0,1,0,-4,-1],
            [1,-2,2,-2,1],
            [0,1,-2,1,-1],
            [1,-1,0,0,2],
            [2,0,0,1,0],
            [1,-2,0,2,0],
            [1,2,-2,-2,-1],
            [0,0,4,-4,1],
            [0,1,2,4,2],
            [0,1,-4,2,-2],
            [3,0,-2,0,0],
            [2,-1,2,2,1],
            [0,1,-2,-4,-1],
            [4,0,2,2,2],
            [2,0,-2,-3,-2],
            [2,0,0,-6,0],
            [1,0,2,0,3],
            [3,1,0,0,0],
            [3,0,0,-4,1],
            [1,-1,2,0,0],
            [1,-1,0,-4,0],
            [2,0,-2,2,-2],
            [1,1,0,-2,2],
            [4,0,0,-2,0],
            [2,2,0,-2,0],
            [0,1,2,0,0],
            [1,1,0,-4,1],
            [1,0,0,-4,-2],
            [0,0,0,1,2],
            [3,0,0,2,1],
            [1,1,0,-4,-1],
            [0,0,2,2,-1],
            [1,1,2,0,0],
            [1,-1,2,-4,1],
            [1,1,0,0,2],
            [0,0,2,6,2],
            [4,0,-2,-2,-1],
            [2,1,0,-4,-1],
            [0,0,0,3,1],
            [1,-1,-2,0,0],
            [0,0,2,1,0],
            [1,0,0,2,-2],
            [3,-1,2,2,2],
            [3,-1,2,-2,2],
            [1,0,0,-1,2],
            [1,-2,2,-2,2],
            [0,1,0,2,2],
            [0,1,-2,-1,-2],
            [1,1,-2,0,0],
            [0,2,2,-2,0],
            [3,-1,-2,-1,-2],
            [1,0,0,-6,0],
            [1,0,-2,-4,0],
            [2,1,0,-4,1],
            [2,0,2,0,-1],
            [2,0,-4,0,-1],
            [0,0,3,0,2],
            [2,1,-2,-2,-1],
            [1,-2,0,0,1],
            [2,-1,0,-4,0],
            [0,0,0,3,0],
            [5,0,2,-2,2],
            [1,2,-2,-4,-2],
            [1,0,4,-4,2],
            [0,0,4,-1,2],
            [3,1,0,-4,0],
            [3,0,0,-6,0],
            [2,0,0,2,2],
            [2,-2,2,0,2],
            [1,0,0,-3,1],
            [1,-2,-2,0,-2],
            [1,-1,-2,-3,-2],
            [0,0,2,-2,-2],
            [2,0,-2,-4,0],
            [1,0,-4,0,0],
            [0,1,0,-1,0],
            [4,0,0,0,-1],
            [3,0,2,-1,2],
            [3,-1,2,0,1],
            [2,0,2,-1,1],
            [1,2,2,-2,1],
            [1,1,0,2,-1],
            [0,2,2,0,1],
            [3,1,2,0,1],
            [1,1,2,1,1],
            [1,1,0,-1,1],
            [1,-2,0,-2,-1],
            [4,0,0,-4,0],
            [2,1,0,2,0],
            [1,-1,0,4,0],
            [0,1,0,-2,2],
            [0,0,2,0,-2],
            [1,0,-1,0,1],
            [3,0,2,-2,0],
            [2,0,2,2,0],
            [1,2,0,-4,0],
            [1,-1,0,-3,0],
            [0,1,0,4,0],
            [0,1,-2,0,0],
            [2,2,2,-2,2],
            [0,0,0,1,-2],
            [0,2,-2,0,-1],
            [4,0,2,-4,2],
            [2,0,-4,2,-2],
            [2,-1,-2,0,-2],
            [1,1,4,-2,2],
            [1,1,2,-4,2],
            [1,0,2,3,2],
            [1,0,0,4,-1],
            [0,0,0,4,2],
            [2,0,0,4,0],
            [1,1,-2,2,0],
            [2,1,2,1,2],
            [2,1,2,-4,1],
            [2,0,2,1,1],
            [2,0,-4,-2,-1],
            [2,0,-2,-6,-1],
            [2,-1,2,-1,2],
            [1,-2,2,0,1],
            [1,-2,0,-2,1],
            [1,-1,0,-4,-1],
            [0,2,2,2,2],
            [0,2,-2,-4,-2],
            [0,1,2,3,2],
            [0,1,0,-4,1],
            [3,0,0,-2,1],
            [2,1,-2,0,1],
            [2,0,4,-2,1],
            [2,0,0,-3,-1],
            [2,-2,0,-2,1],
            [2,-1,2,-2,1],
            [1,0,0,-6,-1],
            [1,-2,0,0,-1],
            [1,-2,-2,-2,-1],
            [0,1,4,-2,1],
            [0,0,2,3,1],
            [2,-1,0,-1,0],
            [1,3,0,-2,0],
            [0,3,0,-2,0],
            [2,-2,2,-2,2],
            [0,0,4,-2,0],
            [4,-1,2,0,2],
            [2,2,-2,-4,-2],
            [4,1,2,0,2],
            [4,-1,-2,-2,-2],
            [2,1,0,-2,-2],
            [2,1,-2,-6,-2],
            [2,0,0,-1,1],
            [2,-1,-2,2,-1],
            [1,1,-2,2,-2],
            [1,1,-2,-3,-2],
            [1,0,3,0,3],
            [1,0,-2,1,1],
            [1,0,-2,0,2],
            [1,-1,2,1,2],
            [1,-1,0,0,-2],
            [1,-1,-4,2,-2],
            [0,3,-2,-2,-2],
            [0,1,0,4,1],
            [0,0,4,2,2],
            [3,0,-2,-2,0],
            [2,-2,0,0,0],
            [1,1,2,-4,0],
            [1,1,0,-3,0],
            [1,0,2,-3,0],
            [1,-1,2,-2,0],
            [0,2,0,2,0],
            [0,0,2,4,0],
            [1,0,1,0,0],
            [3,1,2,-2,1],
            [3,0,4,-2,2],
            [3,0,2,1,2],
            [3,0,0,2,-1],
            [3,0,0,0,2],
            [3,0,-2,2,-1],
            [2,0,4,-4,2],
            [2,0,2,-3,2],
            [2,0,0,4,1],
            [2,0,0,-3,1],
            [2,0,-4,2,-1],
            [2,0,-2,-2,1],
            [2,-2,2,2,2],
            [2,-2,0,-2,-2],
            [2,-1,0,2,1],
            [2,-1,0,2,-1],
            [1,1,2,4,2],
            [1,1,0,1,1],
            [1,1,0,1,-1],
            [1,1,-2,-6,-2],
            [1,0,0,-3,-1],
            [1,0,-4,-2,-1],
            [1,0,-2,-6,-1],
            [1,-2,2,2,1],
            [1,-2,-2,2,-1],
            [1,-1,-2,-4,-1],
            [0,2,0,0,2],
            [0,1,2,-4,2],
            [0,1,-2,4,-1],
            [5,0,0,0,0],
            [3,0,0,-3,0],
            [2,2,0,-4,0],
            [1,-1,2,2,0],
            [0,1,0,3,0],
            [4,0,-2,0,-1],
            [3,0,-2,-6,-1],
            [3,0,-2,-1,-1],
            [2,1,2,2,1],
            [2,1,0,2,1],
            [2,0,2,4,1],
            [2,0,2,-6,1],
            [2,0,2,-2,-1],
            [2,0,0,-6,-1],
            [2,-1,-2,-2,-1],
            [1,2,2,0,1],
            [1,2,0,0,1],
            [1,0,4,0,1],
            [1,0,2,-6,1],
            [1,0,2,-4,-1],
            [1,0,-1,-2,-1],
            [1,-1,2,4,1],
            [1,-1,2,-3,1],
            [1,-1,0,4,1],
            [1,-1,-2,1,-1],
            [0,1,2,-2,3],
            [3,0,0,-2,0],
            [1,0,1,-2,0],
            [0,2,0,-4,0],
            [0,0,2,-4,0],
            [0,0,1,-1,0],
            [0,0,0,6,0],
            [0,2,0,0,-2],
            [0,1,-2,2,-3],
            [4,0,0,2,0],
            [3,0,0,-1,0],
            [3,-1,0,2,0],
            [2,1,0,1,0],
            [2,1,0,-6,0],
            [2,-1,2,0,0],
            [1,0,2,-1,0],
            [1,-1,0,1,0],
            [1,-1,-2,-2,0],
            [0,1,2,2,0],
            [0,0,2,-3,0],
            [2,2,0,-2,-1],
            [2,-1,-2,0,1],
            [1,2,2,-4,1],
            [0,1,4,-4,2],
            [0,0,0,3,2],
            [5,0,2,0,1],
            [4,1,2,-2,2],
            [4,0,-2,-2,0],
            [3,1,2,2,2],
            [3,1,0,-2,0],
            [3,1,-2,-6,-2],
            [3,0,0,0,-2],
            [3,0,-2,-4,-2],
            [3,-1,0,-3,0],
            [3,-1,0,-2,0],
            [2,1,2,0,0],
            [2,1,2,-4,2],
            [2,1,2,-2,0],
            [2,1,0,-3,0],
            [2,1,-2,0,-2],
            [2,0,0,-4,2],
            [2,0,0,-4,-2],
            [2,0,-2,-5,-2],
            [2,-1,2,4,2],
            [2,-1,0,-2,2],
            [1,3,-2,-2,-2],
            [1,1,0,0,-2],
            [1,1,0,-6,0],
            [1,1,-2,1,-2],
            [1,1,-2,-1,-2],
            [1,0,2,1,0],
            [1,0,0,3,0],
            [1,0,0,-4,2],
            [1,0,-2,4,-2],
            [1,-2,0,-1,0],
            [0,1,-4,2,-1],
            [1,0,-2,0,-3],
            [0,0,4,-4,4]]
    
    #基本参数:行星项.
    MFAPL=[[0,0,1,-1,1,0,0,-1,0,-2,5,0,0,0],
            [0,0,0,0,0,0,0,0,0,2,-5,0,0,-1],
            [0,0,0,0,0,0,3,-5,0,0,0,0,0,-2],
            [0,0,1,-1,1,0,-8,12,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,2,0,0,0,2],
            [0,0,0,0,0,0,0,4,-8,3,0,0,0,0],
            [0,0,0,0,0,0,1,-1,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,8,-16,4,5,0,0,0],
            [0,0,0,0,0,0,0,1,0,-1,0,0,0,0],
            [0,0,0,0,1,0,0,-1,2,0,0,0,0,0],
            [0,0,0,0,0,0,8,-13,0,0,0,0,0,-1],
            [0,0,1,-1,1,0,0,-1,0,2,-5,0,0,0],
            [0,0,2,-2,1,0,-5,6,0,0,0,0,0,0],
            [0,0,0,0,0,0,4,-6,0,0,0,0,0,-2],
            [0,0,0,0,0,0,0,3,0,-1,0,0,0,2],
            [0,0,0,0,0,0,0,2,-8,3,0,0,0,-2],
            [0,0,0,0,0,0,2,-4,0,0,0,0,0,-2],
            [0,0,0,0,0,0,0,6,-8,3,0,0,0,2],
            [0,0,0,0,0,0,0,1,-2,0,0,0,0,0],
            [0,0,0,0,0,0,2,-3,0,0,0,0,0,0],
            [0,0,0,0,0,0,2,-2,0,0,0,0,0,0],
            [0,0,0,0,0,0,2,0,0,0,0,0,0,2],
            [0,0,0,0,1,0,0,-4,8,-3,0,0,0,0],
            [0,0,0,0,1,0,0,4,-8,3,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,2,-5,0,0,0],
            [0,0,0,0,0,0,1,1,0,0,0,0,0,2],
            [0,0,1,-1,1,0,0,0,-2,0,0,0,0,0],
            [2,0,0,-2,-1,0,0,-2,0,2,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,2,0,0,0,1],
            [2,0,0,-2,0,0,0,-2,0,2,0,0,0,0],
            [0,0,0,0,0,0,0,2,0,-2,0,0,0,0],
            [0,0,0,0,0,0,8,-13,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,1,0,1,0,0,0,2],
            [0,0,0,0,0,0,5,-8,0,0,0,0,0,-2],
            [0,0,0,0,0,0,0,2,-2,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,2,-5,0,0,1],
            [2,0,0,-2,0,0,0,-2,0,3,0,0,0,0],
            [0,0,1,-1,1,0,0,-1,0,-1,0,0,0,0],
            [0,0,0,0,0,0,3,-4,0,0,0,0,0,0],
            [0,0,1,-1,1,0,0,-1,0,0,-1,0,0,0],
            [0,0,0,0,0,0,0,1,0,-2,0,0,0,0],
            [0,0,0,0,0,0,5,-7,0,0,0,0,0,-2],
            [0,0,1,-1,0,0,0,0,-2,0,0,0,0,0],
            [0,0,0,0,0,0,0,4,0,-2,0,0,0,2],
            [0,0,0,0,0,0,8,-13,0,0,0,0,0,-2],
            [0,0,0,0,0,0,0,0,0,1,0,0,0,0],
            [0,0,0,0,0,0,2,-1,0,0,0,0,0,2],
            [1,0,0,0,0,0,-18,16,0,0,0,0,0,0],
            [0,0,1,-1,1,0,0,-1,0,2,0,0,0,0],
            [0,0,0,0,0,0,0,2,0,1,0,0,0,2],
            [0,0,1,-1,1,0,-5,7,0,0,0,0,0,0],
            [1,0,0,0,0,0,-10,3,0,0,0,0,0,0],
            [0,0,2,-2,0,0,-5,6,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,2,0,-1,0,0,0,2],
            [1,0,2,0,2,0,0,1,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,4,-2,0,0,0,0,2],
            [0,0,0,0,0,0,0,0,0,0,2,0,0,1],
            [1,0,-2,0,-2,0,0,4,-8,3,0,0,0,0],
            [0,0,1,-1,1,0,0,-1,0,0,2,0,0,0],
            [0,0,2,-2,1,0,-3,3,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,0,2,0,0,2],
            [0,0,0,0,0,0,0,8,-16,4,5,0,0,-2],
            [0,0,1,-1,1,0,0,3,-8,3,0,0,0,0],
            [0,0,0,0,0,0,8,-11,0,0,0,0,0,-2],
            [0,0,0,0,0,0,0,0,0,3,0,0,0,2],
            [0,0,0,0,0,0,0,8,-16,4,5,0,0,2],
            [0,0,0,0,0,0,1,-1,0,0,0,0,0,-1],
            [0,0,0,0,0,0,4,-6,0,0,0,0,0,-1],
            [0,0,0,0,0,0,0,1,0,-3,0,0,0,-2],
            [0,0,0,0,0,0,0,2,-4,0,0,0,0,0],
            [0,0,0,0,0,0,6,-8,0,0,0,0,0,-2],
            [0,0,0,0,0,0,3,-2,0,0,0,0,0,2],
            [0,0,0,0,0,0,8,-15,0,0,0,0,0,-2],
            [0,0,0,0,0,0,2,-5,0,0,0,0,0,-2],
            [0,0,0,0,0,0,1,-3,0,0,0,0,0,-2],
            [0,0,0,0,0,0,0,3,0,-2,0,0,0,2],
            [0,0,1,-1,1,0,0,-5,8,-3,0,0,0,0],
            [0,0,0,0,0,0,0,1,2,0,0,0,0,2],
            [0,0,0,0,0,0,0,3,-2,0,0,0,0,2],
            [0,0,0,0,0,0,3,-5,0,0,0,0,0,0],
            [2,0,0,-2,1,0,0,-2,0,3,0,0,0,0],
            [0,0,0,0,0,0,5,-8,0,0,0,0,0,-1],
            [2,0,0,-2,0,0,-3,3,0,0,0,0,0,0],
            [0,0,0,0,1,0,8,-13,0,0,0,0,0,0],
            [0,0,0,0,1,0,0,0,0,-2,5,0,0,0],
            [1,0,0,-1,0,0,-3,4,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,2,0,0,0,0,0,2],
            [1,0,0,0,-1,0,-18,16,0,0,0,0,0,0],
            [0,0,0,0,1,0,0,0,0,2,-5,0,0,0],
            [0,0,0,0,1,0,0,0,0,1,0,0,0,0],
            [1,0,0,-2,0,0,19,-21,3,0,0,0,0,0],
            [0,0,0,0,1,0,-8,13,0,0,0,0,0,0],
            [0,0,1,-1,1,0,0,-1,0,0,1,0,0,0],
            [0,0,0,0,0,0,7,-9,0,0,0,0,0,-2],
            [0,0,0,0,0,0,0,0,2,0,0,0,0,2],
            [1,0,0,0,1,0,-18,16,0,0,0,0,0,0],
            [0,0,0,0,0,0,2,-4,0,0,0,0,0,-1],
            [0,0,0,0,0,0,0,6,-16,4,5,0,0,-2],
            [0,0,0,0,0,0,4,-7,0,0,0,0,0,-2],
            [0,0,0,0,0,0,3,-7,0,0,0,0,0,-2],
            [0,0,0,0,0,0,2,-2,0,0,0,0,0,-1],
            [0,0,0,0,0,0,0,0,0,1,0,0,0,1],
            [2,0,0,-2,1,0,0,-2,0,2,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,1,0,0,0,-1],
            [0,0,0,0,0,0,0,3,-4,0,0,0,0,0],
            [0,0,0,0,0,0,1,-2,0,0,0,0,0,0],
            [2,0,0,-2,-1,0,0,-2,0,3,0,0,0,0],
            [0,0,0,0,0,0,3,-3,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,2,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,1,0,2,0,0,0,2],
            [0,0,0,0,1,0,0,1,-2,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,1,0,0,0,2],
            [0,0,2,-2,1,0,0,-2,0,2,0,0,0,0],
            [0,0,0,0,0,0,0,2,0,-3,0,0,0,0],
            [0,0,0,0,0,0,3,-5,0,0,0,0,0,-1],
            [0,0,0,0,0,0,3,-3,0,0,0,0,0,2],
            [0,0,0,0,0,0,4,-4,0,0,0,0,0,0],
            [0,0,1,-1,0,0,0,-1,0,-1,0,0,0,0],
            [2,0,0,-2,0,0,-6,8,0,0,0,0,0,0],
            [0,0,1,-1,1,0,0,-2,2,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,0,1,0,0,1],
            [0,0,1,-1,1,0,0,-1,0,1,0,0,0,0],
            [0,0,0,0,0,0,0,1,-2,0,0,0,0,-1],
            [0,0,0,0,0,0,0,2,-3,0,0,0,0,0],
            [0,0,0,0,0,0,0,2,-4,0,0,0,0,-2],
            [0,0,0,0,0,0,0,1,0,0,-1,0,0,0],
            [0,0,0,0,0,0,8,-10,0,0,0,0,0,-2],
            [0,0,1,-1,1,0,-3,4,0,0,0,0,0,0],
            [0,0,0,0,0,0,6,-9,0,0,0,0,0,-2],
            [1,0,0,-1,1,0,0,-1,0,2,0,0,0,0],
            [0,0,0,0,0,0,5,-7,0,0,0,0,0,-1],
            [0,0,0,0,0,0,5,-5,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,0,1,0,0,-1],
            [0,0,0,0,0,0,3,-3,0,0,0,0,0,-1],
            [0,0,0,0,0,0,0,0,0,0,1,0,0,0],
            [0,0,0,0,0,0,0,0,4,0,0,0,0,2],
            [0,0,0,0,0,0,0,4,0,-3,0,0,0,2],
            [0,0,0,0,0,0,1,-1,0,0,0,0,0,1],
            [0,0,0,0,0,0,0,2,0,0,0,0,0,1],
            [0,0,0,0,1,0,2,-3,0,0,0,0,0,0],
            [1,0,0,-1,0,0,0,-1,0,1,0,0,0,0],
            [0,0,0,0,0,0,1,-3,0,0,0,0,0,-1],
            [0,0,0,0,0,0,0,5,-4,0,0,0,0,2],
            [0,0,0,0,0,0,0,4,-4,0,0,0,0,2],
            [0,0,0,0,0,0,9,-11,0,0,0,0,0,-2],
            [0,0,0,0,0,0,2,-3,0,0,0,0,0,-1],
            [0,0,0,0,0,0,0,8,-15,0,0,0,0,0],
            [0,0,1,-1,1,0,-4,5,0,0,0,0,0,0],
            [0,0,0,0,0,0,4,-6,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,4,0,-1,0,0,0,2],
            [1,0,0,-1,1,0,-3,4,0,0,0,0,0,0],
            [0,0,1,1,1,0,0,1,0,0,0,0,0,0],
            [0,0,1,-1,1,0,0,-1,0,-4,10,0,0,0],
            [0,0,0,0,1,0,1,-1,0,0,0,0,0,0],
            [0,0,1,-1,0,0,0,-1,0,0,-1,0,0,0],
            [0,0,0,0,0,0,0,1,0,-3,0,0,0,0],
            [0,0,0,0,0,0,3,-1,0,0,0,0,0,2],
            [0,0,0,0,0,0,0,1,0,-4,0,0,0,-2],
            [0,0,0,0,0,0,0,0,0,2,-5,0,0,-2],
            [0,0,2,-2,1,0,-4,4,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,3,0,0,-1,0,0,2],
            [0,0,0,0,0,0,0,4,-3,0,0,0,0,2],
            [0,0,1,-1,1,0,0,-1,0,0,0,0,2,0],
            [0,0,0,0,0,0,4,-4,0,0,0,0,0,-1],
            [0,0,0,0,0,0,0,2,-4,0,0,0,0,-1],
            [0,0,0,0,0,0,5,-8,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,1,-2,0,0,0,0,1],
            [0,0,0,0,1,0,0,0,0,0,1,0,0,0],
            [0,0,2,-2,1,0,0,-9,13,0,0,0,0,0],
            [2,0,2,0,2,0,0,2,0,-3,0,0,0,0],
            [0,0,0,0,0,0,3,-6,0,0,0,0,0,-2],
            [0,0,1,-1,2,0,0,-1,0,0,2,0,0,0],
            [1,0,0,-1,-1,0,-3,4,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,3,-6,0,0,0,0,-2],
            [0,0,0,0,0,0,6,-6,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,3,0,0,0,1],
            [1,0,2,0,1,0,0,-2,0,3,0,0,0,0],
            [1,0,-2,0,-1,0,0,-1,0,0,0,0,0,0],
            [0,0,0,0,1,0,0,-2,4,0,0,0,0,0],
            [0,0,0,0,0,0,0,3,-5,0,0,0,0,0],
            [0,0,0,0,0,0,2,1,0,0,0,0,0,2],
            [0,0,0,0,0,0,1,1,0,0,0,0,0,1],
            [0,0,2,0,2,0,0,1,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,1,-8,3,0,0,0,-2],
            [0,0,0,0,0,0,6,-10,0,0,0,0,0,-2],
            [0,0,0,0,0,0,0,7,-8,3,0,0,0,2],
            [0,0,0,0,1,0,-3,5,0,0,0,0,0,0],
            [0,0,1,-1,1,0,-1,0,0,0,0,0,0,0],
            [0,0,1,-1,0,0,-5,7,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,2,0,-2,0,0,0,1],
            [0,0,0,0,0,0,0,2,0,-1,0,0,0,0],
            [0,0,0,0,0,0,7,-10,0,0,0,0,0,-2],
            [1,0,0,-2,0,0,0,-2,0,2,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,2,0,0,0,0],
            [0,0,0,0,0,0,0,1,0,2,-5,0,0,0],
            [0,0,0,0,0,0,6,-8,0,0,0,0,0,-1],
            [0,0,1,-1,1,0,0,-9,15,0,0,0,0,0],
            [0,0,0,0,1,0,-2,3,0,0,0,0,0,0],
            [0,0,0,0,1,0,-1,1,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,3,-6,0,0,0,0,0],
            [0,0,0,0,0,0,0,1,-4,0,0,0,0,-2],
            [0,0,0,0,0,0,0,0,3,0,0,0,0,2],
            [0,0,0,0,0,0,0,2,0,0,-1,0,0,2],
            [2,0,0,-2,1,0,-6,8,0,0,0,0,0,0],
            [0,0,0,0,0,0,5,-5,0,0,0,0,0,-1],
            [0,0,1,-1,1,0,3,-6,0,0,0,0,0,0],
            [0,0,1,-1,1,0,-2,2,0,0,0,0,0,0],
            [0,0,1,-1,1,0,8,-14,0,0,0,0,0,0],
            [0,0,0,0,0,0,1,0,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,1,0,0,0,0,0,0],
            [0,0,0,0,1,0,0,8,-15,0,0,0,0,0],
            [0,0,0,0,0,0,0,4,-6,0,0,0,0,0],
            [0,0,0,0,0,0,7,-7,0,0,0,0,0,0],
            [2,0,0,-2,1,0,-3,3,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,3,-1,0,0,0,0,2],
            [0,0,0,0,0,0,0,2,0,0,1,0,0,2],
            [2,0,-1,-1,0,0,0,3,-7,0,0,0,0,0],
            [0,0,0,0,0,0,0,4,-7,0,0,0,0,-2],
            [0,0,0,0,0,0,0,3,-3,0,0,0,0,0],
            [0,0,1,-1,1,0,0,-3,4,0,0,0,0,0],
            [2,0,0,-2,0,0,0,-6,8,0,0,0,0,0],
            [2,0,0,-2,0,0,0,-5,6,0,0,0,0,0],
            [0,0,0,0,1,0,0,0,0,-1,0,0,0,0],
            [0,0,0,0,0,0,2,0,0,0,0,0,0,1],
            [0,0,0,0,0,0,2,1,0,0,0,0,0,1],
            [0,0,0,0,0,0,1,2,0,0,0,0,0,2],
            [0,0,0,0,1,0,0,1,0,-1,0,0,0,0],
            [0,0,0,0,0,0,0,1,-1,0,0,0,0,0],
            [0,0,0,0,0,0,3,-9,4,0,0,0,0,-2],
            [0,0,0,0,0,0,0,3,-5,0,0,0,0,-2],
            [0,0,0,0,0,0,0,2,0,-4,0,0,0,-2],
            [0,0,0,0,0,0,0,0,0,0,0,0,2,1],
            [0,0,0,0,0,0,7,-11,0,0,0,0,0,-2],
            [0,0,0,0,0,0,3,-5,4,0,0,0,0,2],
            [0,0,1,-1,0,0,0,-1,0,-1,1,0,0,0],
            [2,0,0,0,0,0,0,-2,0,3,0,0,0,0],
            [0,0,0,0,0,0,0,8,-15,0,0,0,0,-2],
            [0,0,1,-1,2,0,0,-2,2,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,0,3,0,0,2],
            [0,0,0,0,0,0,6,-6,0,0,0,0,0,-1],
            [0,0,1,-1,1,0,0,-1,0,-1,1,0,0,0],
            [0,0,0,0,0,0,2,-2,0,0,0,0,0,1],
            [0,0,0,0,0,0,0,4,-7,0,0,0,0,0],
            [0,0,0,0,0,0,0,3,-8,3,0,0,0,0],
            [0,0,1,-1,1,0,2,-4,0,-3,0,0,0,0],
            [0,0,0,0,1,0,3,-5,0,2,0,0,0,0],
            [0,0,0,0,0,0,0,3,0,-3,0,0,0,2],
            [0,0,2,-2,2,0,-8,11,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,5,-8,3,0,0,0,0],
            [0,0,0,0,0,0,0,1,0,0,-2,0,0,0],
            [0,0,0,0,0,0,0,1,0,0,1,0,0,2],
            [0,0,0,0,0,0,0,5,-9,0,0,0,0,-2],
            [0,0,0,0,0,0,0,5,-5,0,0,0,0,2],
            [0,0,0,0,0,0,7,-9,0,0,0,0,0,-1],
            [0,0,0,0,0,0,4,-7,0,0,0,0,0,-1],
            [0,0,0,0,0,0,2,-1,0,0,0,0,0,0],
            [1,0,-2,-2,-2,0,0,-2,0,2,0,0,0,0],
            [0,0,0,0,0,0,0,1,1,0,0,0,0,2],
            [0,0,0,0,0,0,0,2,0,-2,5,0,0,2],
            [0,0,0,0,0,0,3,-3,0,0,0,0,0,1],
            [0,0,0,0,0,0,0,6,0,0,0,0,0,2],
            [0,0,0,0,0,0,0,2,0,2,-5,0,0,2],
            [2,0,0,-2,-1,0,0,-2,0,0,5,0,0,0],
            [2,0,0,-2,-1,0,-6,8,0,0,0,0,0,0],
            [1,0,0,-2,0,0,-3,3,0,0,0,0,0,0],
            [0,0,0,0,0,0,8,-8,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,3,0,2,-5,0,0,2],
            [0,0,0,0,1,0,3,-7,4,0,0,0,0,0],
            [0,0,2,-2,1,0,-2,2,0,0,0,0,0,0],
            [0,0,0,0,1,0,0,-1,0,1,0,0,0,0],
            [0,0,1,-1,0,0,0,-1,0,-2,5,0,0,0],
            [0,0,0,0,0,0,0,3,0,-3,0,0,0,0],
            [0,0,0,0,0,0,3,-1,0,0,0,0,0,1],
            [0,0,0,0,0,0,2,-3,0,0,0,0,0,-2],
            [0,0,0,0,0,0,0,11,0,0,0,0,0,2],
            [0,0,0,0,0,0,0,6,-15,0,0,0,0,-2],
            [0,0,0,0,0,0,0,3,0,1,0,0,0,2],
            [1,0,0,-1,0,0,0,-3,4,0,0,0,0,0],
            [0,0,0,0,1,0,-3,7,-4,0,0,0,0,0],
            [0,0,0,0,0,0,0,5,0,-2,0,0,0,2],
            [0,0,0,0,0,0,3,-5,0,0,0,0,0,1],
            [0,0,2,-2,2,0,-5,6,0,0,0,0,0,0],
            [0,0,2,-2,2,0,-3,3,0,0,0,0,0,0],
            [0,0,0,0,0,0,3,0,0,0,0,0,0,2],
            [0,0,0,0,0,0,0,6,0,0,0,0,0,0],
            [0,0,0,0,0,0,4,-4,0,0,0,0,0,2],
            [0,0,0,0,0,0,0,4,-8,0,0,0,0,-2],
            [0,0,0,0,0,0,0,4,-5,0,0,0,0,0],
            [0,0,0,0,0,0,5,-7,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,6,-11,0,0,0,0,-2],
            [0,0,0,0,0,0,0,1,-3,0,0,0,0,-2],
            [0,0,1,-1,1,0,0,-1,0,3,0,0,0,0],
            [0,0,1,-1,0,0,0,-1,0,2,0,0,0,0],
            [0,0,0,0,0,0,1,-2,0,0,0,0,0,1],
            [0,0,0,0,0,0,9,-12,0,0,0,0,0,-2],
            [0,0,0,0,0,0,4,-4,0,0,0,0,0,1],
            [0,0,1,-1,0,0,-8,12,0,0,0,0,0,0],
            [0,0,1,-1,1,0,-2,3,0,0,0,0,0,0],
            [0,0,0,0,0,0,7,-7,0,0,0,0,0,-1],
            [0,0,0,0,0,0,0,3,-6,0,0,0,0,-1],
            [0,0,0,0,0,0,0,6,-6,0,0,0,0,2],
            [0,0,0,0,0,1,0,-4,0,0,0,0,0,-2],
            [0,0,1,-1,1,0,0,1,0,0,0,0,0,0],
            [0,0,0,0,0,0,6,-9,0,0,0,0,0,-1],
            [0,0,1,-1,-1,0,0,0,-2,0,0,0,0,0],
            [0,0,0,0,0,0,0,1,-5,0,0,0,0,-2],
            [2,0,0,-2,0,0,0,-2,0,3,-1,0,0,0],
            [0,0,0,0,0,0,0,2,0,0,-2,0,0,0],
            [0,0,0,0,0,0,0,5,-9,0,0,0,0,0],
            [0,0,0,0,0,0,5,-6,0,0,0,0,0,2],
            [0,0,0,0,0,0,9,-9,0,0,0,0,0,-1],
            [0,0,1,-1,1,0,0,-1,0,0,3,0,0,0],
            [0,0,0,0,1,0,0,2,-4,0,0,0,0,0],
            [0,0,0,0,0,0,5,-3,0,0,0,0,0,2],
            [0,0,0,0,0,0,0,0,0,0,3,0,0,1],
            [0,0,1,-1,2,0,0,-1,0,2,0,0,0,0],
            [0,0,0,0,0,0,5,-9,0,0,0,0,0,-2],
            [0,0,0,0,0,0,0,5,-3,0,0,0,0,2],
            [0,0,0,0,0,0,0,0,0,4,0,0,0,2],
            [0,0,2,0,2,0,0,4,-8,3,0,0,0,0],
            [0,0,2,0,2,0,0,-4,8,-3,0,0,0,0],
            [0,0,0,0,0,0,0,5,0,-3,0,0,0,2],
            [0,0,0,0,0,0,0,1,0,1,0,0,0,0],
            [2,0,-1,-1,-1,0,0,-1,0,3,0,0,0,0],
            [0,0,0,0,0,0,4,-3,0,0,0,0,0,2],
            [0,0,0,0,0,0,4,-2,0,0,0,0,0,2],
            [0,0,0,0,0,0,5,-10,0,0,0,0,0,-2],
            [0,0,0,0,0,0,8,-13,0,0,0,0,0,1],
            [0,0,2,-2,1,-1,0,2,0,0,0,0,0,0],
            [0,0,1,-1,1,0,0,-1,0,0,0,2,0,0],
            [0,0,0,0,1,0,3,-5,0,0,0,0,0,0],
            [1,0,0,-2,0,0,0,-2,0,3,0,0,0,0],
            [0,0,2,-2,0,0,-3,3,0,0,0,0,0,0],
            [0,0,0,0,0,0,9,-9,0,0,0,0,0,0],
            [0,0,2,0,2,0,1,-1,0,0,0,0,0,0],
            [0,0,2,-2,1,0,0,-8,11,0,0,0,0,0],
            [0,0,2,-2,1,0,0,-2,0,0,2,0,0,0],
            [0,0,1,-1,1,0,0,-1,0,-1,2,0,0,0],
            [0,0,0,0,0,0,5,-5,0,0,0,0,0,2],
            [0,0,0,0,0,0,2,-6,0,0,0,0,0,-2],
            [0,0,0,0,0,0,0,8,-15,0,0,0,0,-1],
            [0,0,0,0,0,0,0,5,-2,0,0,0,0,2],
            [0,0,0,0,0,0,0,0,0,0,1,0,0,2],
            [0,0,0,0,0,0,0,7,-13,0,0,0,0,-2],
            [0,0,0,0,0,0,0,3,0,-2,0,0,0,0],
            [0,0,0,0,0,0,0,1,0,3,0,0,0,2],
            [0,0,2,-2,1,0,0,-2,0,3,0,0,0,0],
            [0,0,0,0,0,0,8,-8,0,0,0,0,0,-1],
            [0,0,0,0,0,0,8,-10,0,0,0,0,0,-1],
            [0,0,0,0,0,0,4,-2,0,0,0,0,0,1],
            [0,0,0,0,0,0,3,-6,0,0,0,0,0,-1],
            [0,0,0,0,0,0,3,-4,0,0,0,0,0,-1],
            [0,0,0,0,0,0,0,0,0,2,-5,0,0,2],
            [0,0,0,0,0,0,1,0,0,0,0,0,0,2],
            [0,0,0,0,0,0,0,2,0,-4,0,0,0,0],
            [2,0,0,-2,-1,0,0,-5,6,0,0,0,0,0],
            [0,0,0,0,0,0,0,2,-5,0,0,0,0,-2],
            [2,0,-1,-1,-1,0,0,3,-7,0,0,0,0,0],
            [0,0,0,0,0,0,0,5,-8,0,0,0,0,0],
            [0,0,2,0,2,0,-1,1,0,0,0,0,0,0],
            [2,0,0,-2,0,0,0,-2,0,4,-3,0,0,0],
            [0,0,0,0,0,0,0,6,-11,0,0,0,0,0],
            [2,0,0,-2,1,0,0,-6,8,0,0,0,0,0],
            [0,0,0,0,0,0,0,4,-8,1,5,0,0,-2],
            [0,0,0,0,0,0,0,6,-5,0,0,0,0,2],
            [1,0,-2,-2,-2,0,-3,3,0,0,0,0,0,0],
            [0,0,1,-1,2,0,0,0,-2,0,0,0,0,0],
            [0,0,0,0,2,0,0,4,-8,3,0,0,0,0],
            [0,0,0,0,2,0,0,-4,8,-3,0,0,0,0],
            [0,0,0,0,0,0,0,6,0,0,0,0,0,1],
            [0,0,0,0,0,0,0,6,-7,0,0,0,0,2],
            [0,0,0,0,0,0,0,4,0,0,-2,0,0,2],
            [0,0,0,0,0,0,0,3,0,0,-2,0,0,2],
            [0,0,0,0,0,0,0,1,0,-1,0,0,0,1],
            [0,0,0,0,0,0,0,1,-6,0,0,0,0,-2],
            [0,0,0,0,0,0,0,0,0,4,-5,0,0,2],
            [0,0,0,0,0,0,0,0,0,0,0,2,0,2],
            [0,0,0,0,0,0,3,-5,0,2,0,0,0,0],
            [0,0,0,0,0,0,0,7,-13,0,0,0,0,0],
            [0,0,0,0,0,0,0,2,0,-2,0,0,0,2],
            [0,0,1,-1,0,0,0,-1,0,0,2,0,0,0],
            [0,0,0,0,1,0,0,-8,15,0,0,0,0,0],
            [2,0,0,-2,-2,0,-3,3,0,0,0,0,0,0],
            [2,0,-1,-1,-1,0,0,-1,0,2,0,0,0,0],
            [1,0,2,-2,2,0,0,-2,0,2,0,0,0,0],
            [1,0,-1,1,-1,0,-18,17,0,0,0,0,0,0],
            [0,0,2,0,2,0,0,1,0,-1,0,0,0,0],
            [0,0,2,0,2,0,0,-1,0,1,0,0,0,0],
            [0,0,2,-2,-1,0,-5,6,0,0,0,0,0,0],
            [0,0,1,-1,2,0,0,-1,0,1,0,0,0,0],
            [0,0,0,0,1,0,2,-2,0,0,0,0,0,0],
            [0,0,0,0,0,0,8,-16,0,0,0,0,0,-2],
            [0,0,0,0,0,0,0,0,0,0,5,0,0,2],
            [0,0,0,0,0,0,0,0,0,0,0,0,2,2],
            [0,0,0,0,2,0,0,-1,2,0,0,0,0,0],
            [2,0,-1,-1,-2,0,0,-1,0,2,0,0,0,0],
            [0,0,0,0,0,0,6,-10,0,0,0,0,0,-1],
            [0,0,1,-1,1,0,0,-1,0,-2,4,0,0,0],
            [0,0,0,0,0,0,0,2,2,0,0,0,0,2],
            [2,0,0,-2,-1,0,0,-2,0,4,-5,0,0,0],
            [2,0,0,-2,-1,0,-3,3,0,0,0,0,0,0],
            [2,0,-1,-1,-1,0,0,-1,0,0,0,0,0,0],
            [1,0,1,-1,1,0,0,-1,0,0,0,0,0,0],
            [1,0,0,-1,-1,0,0,-2,2,0,0,0,0,0],
            [1,0,-1,-1,-1,0,20,-20,0,0,0,0,0,0],
            [0,0,2,-2,1,0,0,-1,0,1,0,0,0,0],
            [0,0,1,-1,1,0,1,-2,0,0,0,0,0,0],
            [0,0,1,-1,1,0,-2,1,0,0,0,0,0,0],
            [0,0,0,0,1,0,5,-8,0,0,0,0,0,0],
            [0,0,0,0,1,0,0,0,0,0,-1,0,0,0],
            [0,0,0,0,0,0,9,-11,0,0,0,0,0,-1],
            [0,0,0,0,0,0,5,-3,0,0,0,0,0,1],
            [0,0,0,0,0,0,0,1,0,-3,0,0,0,-1],
            [0,0,0,0,0,0,0,0,0,0,0,2,0,1],
            [0,0,0,0,0,0,6,-7,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,3,-2,0,0,0,0,0],
            [0,0,0,0,0,0,1,-2,0,0,0,0,0,-2],
            [0,0,1,-1,1,0,0,-1,0,0,-2,0,0,0],
            [0,0,1,-1,2,0,0,-1,0,-2,5,0,0,0],
            [0,0,0,0,0,0,0,5,-7,0,0,0,0,0],
            [0,0,0,0,0,0,1,-3,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,5,-8,0,0,0,0,-2],
            [0,0,0,0,0,0,0,2,-6,0,0,0,0,-2],
            [1,0,0,-2,0,0,20,-21,0,0,0,0,0,0],
            [0,0,0,0,0,0,8,-12,0,0,0,0,0,0],
            [0,0,0,0,0,0,5,-6,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,4,-4,0,0,0,0,0],
            [0,0,1,-1,2,0,0,-1,0,-1,0,0,0,0],
            [0,0,0,0,0,0,8,-12,0,0,0,0,0,-2],
            [0,0,0,0,0,0,0,9,-17,0,0,0,0,0],
            [0,0,0,0,0,0,0,5,-6,0,0,0,0,2],
            [0,0,0,0,0,0,0,4,-8,1,5,0,0,2],
            [0,0,0,0,0,0,0,4,-6,0,0,0,0,-2],
            [0,0,0,0,0,0,0,2,-7,0,0,0,0,-2],
            [1,0,0,-1,1,0,0,-3,4,0,0,0,0,0],
            [1,0,-2,0,-2,0,-10,3,0,0,0,0,0,0],
            [0,0,0,0,1,0,0,-9,17,0,0,0,0,0],
            [0,0,0,0,0,0,1,-4,0,0,0,0,0,-2],
            [1,0,-2,-2,-2,0,0,-2,0,3,0,0,0,0],
            [1,0,-1,1,-1,0,0,1,0,0,0,0,0,0],
            [0,0,2,-2,2,0,0,-2,0,2,0,0,0,0],
            [0,0,1,-1,2,0,0,-1,0,0,1,0,0,0],
            [0,0,1,-1,2,0,-5,7,0,0,0,0,0,0],
            [0,0,0,0,1,0,0,2,-2,0,0,0,0,0],
            [0,0,0,0,0,0,4,-5,0,0,0,0,0,-1],
            [0,0,0,0,0,0,3,-4,0,0,0,0,0,-2],
            [0,0,0,0,0,0,2,-4,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,5,-10,0,0,0,0,-2],
            [0,0,0,0,0,0,0,4,0,-4,0,0,0,2],
            [0,0,0,0,0,0,0,2,0,-5,0,0,0,-2],
            [0,0,0,0,0,0,0,1,0,-5,0,0,0,-2],
            [0,0,0,0,0,0,0,1,0,-2,5,0,0,2],
            [0,0,0,0,0,0,0,1,0,-2,0,0,0,-2],
            [0,0,0,0,0,0,2,-3,0,0,0,0,0,1],
            [1,0,0,-2,0,0,0,1,0,-1,0,0,0,0],
            [0,0,0,0,0,0,3,-7,4,0,0,0,0,0],
            [2,0,2,0,1,0,0,1,0,0,0,0,0,0],
            [0,0,1,-1,-1,0,0,-1,0,-1,0,0,0,0],
            [0,0,0,0,1,0,0,1,0,-2,0,0,0,0],
            [0,0,0,0,0,0,0,6,-10,0,0,0,0,-2],
            [1,0,0,-1,1,0,0,-1,0,1,0,0,0,0],
            [0,0,2,-2,1,0,0,4,-8,3,0,0,0,0],
            [0,0,2,-2,1,0,0,1,0,-1,0,0,0,0],
            [0,0,2,-2,1,0,0,-4,8,-3,0,0,0,0],
            [0,0,2,-2,1,0,0,-3,0,3,0,0,0,0],
            [0,0,2,-2,1,0,-5,5,0,0,0,0,0,0],
            [0,0,1,-1,1,0,1,-3,0,0,0,0,0,0],
            [0,0,1,-1,1,0,0,-4,6,0,0,0,0,0],
            [0,0,1,-1,1,0,0,-1,0,0,0,-1,0,0],
            [0,0,1,-1,1,0,-5,6,0,0,0,0,0,0],
            [0,0,0,0,1,0,3,-4,0,0,0,0,0,0],
            [0,0,0,0,1,0,-2,2,0,0,0,0,0,0],
            [0,0,0,0,0,0,7,-10,0,0,0,0,0,-1],
            [0,0,0,0,0,0,5,-5,0,0,0,0,0,1],
            [0,0,0,0,0,0,4,-5,0,0,0,0,0,-2],
            [0,0,0,0,0,0,3,-8,0,0,0,0,0,-2],
            [0,0,0,0,0,0,2,-5,0,0,0,0,0,-1],
            [0,0,0,0,0,0,1,-2,0,0,0,0,0,-1],
            [0,0,0,0,0,0,0,7,-9,0,0,0,0,2],
            [0,0,0,0,0,0,0,7,-8,0,0,0,0,2],
            [0,0,0,0,0,0,0,3,0,0,0,0,0,2],
            [0,0,0,0,0,0,0,3,-8,3,0,0,0,-2],
            [0,0,0,0,0,0,0,2,0,0,-2,0,0,1],
            [0,0,0,0,0,0,0,2,-4,0,0,0,0,1],
            [0,0,0,0,0,0,0,1,0,0,0,0,0,-1],
            [0,0,0,0,0,0,0,1,0,-1,0,0,0,-1],
            [2,0,0,-2,-1,0,0,-6,8,0,0,0,0,0],
            [2,0,-1,-1,1,0,0,3,-7,0,0,0,0,0],
            [0,0,2,-2,1,0,0,-7,9,0,0,0,0,0],
            [0,0,0,0,0,0,0,3,-5,0,0,0,0,-1],
            [0,0,1,-1,2,0,-8,12,0,0,0,0,0,0],
            [1,0,0,0,0,0,0,-2,0,2,0,0,0,0],
            [1,0,0,-2,0,0,2,-2,0,0,0,0,0,0],
            [0,0,0,0,0,0,7,-8,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,2,0,0,0,0,0],
            [2,0,0,-2,1,0,0,-5,6,0,0,0,0,0],
            [2,0,0,-2,-1,0,0,-2,0,3,-1,0,0,0],
            [1,0,1,1,1,0,0,1,0,0,0,0,0,0],
            [1,0,0,-2,1,0,0,-2,0,2,0,0,0,0],
            [1,0,0,-2,-1,0,0,-2,0,2,0,0,0,0],
            [1,0,0,-1,-1,0,0,-3,4,0,0,0,0,0],
            [1,0,-1,0,-1,0,-3,5,0,0,0,0,0,0],
            [0,0,2,-2,1,0,0,-4,4,0,0,0,0,0],
            [0,0,2,-2,1,0,0,-2,0,0,0,0,0,0],
            [0,0,2,-2,1,0,-8,11,0,0,0,0,0,0],
            [0,0,2,-2,0,0,0,-9,13,0,0,0,0,0],
            [0,0,1,1,2,0,0,1,0,0,0,0,0,0],
            [0,0,1,-1,1,0,0,1,-4,0,0,0,0,0],
            [0,0,1,-1,1,0,0,-1,0,1,-3,0,0,0],
            [0,0,0,0,1,0,0,7,-13,0,0,0,0,0],
            [0,0,0,0,1,0,0,2,0,-2,0,0,0,0],
            [0,0,0,0,1,0,0,-2,2,0,0,0,0,0],
            [0,0,0,0,1,0,-3,4,0,0,0,0,0,0],
            [0,0,0,0,0,1,0,-4,0,0,0,0,0,0],
            [0,0,0,0,0,0,7,-11,0,0,0,0,0,-1],
            [0,0,0,0,0,0,6,-6,0,0,0,0,0,1],
            [0,0,0,0,0,0,6,-4,0,0,0,0,0,1],
            [0,0,0,0,0,0,5,-6,0,0,0,0,0,-1],
            [0,0,0,0,0,0,4,-2,0,0,0,0,0,0],
            [0,0,0,0,0,0,3,-4,0,0,0,0,0,1],
            [0,0,0,0,0,0,1,-4,0,0,0,0,0,-1],
            [0,0,0,0,0,0,0,9,-17,0,0,0,0,-2],
            [0,0,0,0,0,0,0,7,-7,0,0,0,0,2],
            [0,0,0,0,0,0,0,4,-8,3,0,0,0,1],
            [0,0,0,0,0,0,0,4,-8,3,0,0,0,-1],
            [0,0,0,0,0,0,0,4,-8,0,0,0,0,0],
            [0,0,0,0,0,0,0,4,-7,0,0,0,0,-1],
            [0,0,0,0,0,0,0,1,0,1,0,0,0,1],
            [0,0,0,0,0,0,0,1,0,-4,0,0,0,0],
            [2,0,0,-2,0,0,0,-4,8,-3,0,0,0,0],
            [2,0,0,-2,0,0,-2,2,0,0,0,0,0,0],
            [1,0,0,0,0,0,0,4,-8,3,0,0,0,0],
            [1,0,0,0,0,0,0,-4,8,-3,0,0,0,0],
            [1,0,0,0,0,0,-1,1,0,0,0,0,0,0],
            [1,0,0,-2,0,0,17,-16,0,-2,0,0,0,0],
            [1,0,0,-1,0,0,0,-2,2,0,0,0,0,0],
            [0,0,2,-2,0,0,0,-2,0,2,0,0,0,0],
            [0,0,0,0,0,0,0,6,-9,0,0,0,0,0],
            [0,0,0,0,0,0,0,4,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,3,0,-4,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,0,0,1,-2,-2],
            [0,0,0,0,0,0,0,2,1,0,0,0,0,2],
            [2,0,0,-2,0,0,0,-4,4,0,0,0,0,0],
            [2,0,0,-2,0,0,0,-2,0,2,2,0,0,0],
            [1,0,0,0,0,0,1,-1,0,0,0,0,0,0],
            [1,0,0,0,0,0,0,-1,0,1,0,0,0,0],
            [1,0,0,0,0,0,-3,3,0,0,0,0,0,0],
            [1,0,0,-2,0,0,1,-1,0,0,0,0,0,0],
            [1,0,0,-2,0,0,0,4,-8,3,0,0,0,0],
            [1,0,0,-2,0,0,0,-4,8,-3,0,0,0,0],
            [1,0,0,-2,0,0,-2,2,0,0,0,0,0,0],
            [0,0,2,-2,0,0,-4,4,0,0,0,0,0,0],
            [0,0,1,1,0,0,0,1,0,0,0,0,0,0],
            [0,0,1,-1,0,0,3,-6,0,0,0,0,0,0],
            [0,0,1,-1,0,0,0,-2,2,0,0,0,0,0],
            [0,0,1,-1,0,0,0,-1,0,1,0,0,0,0],
            [0,0,1,-1,0,0,0,-1,0,0,1,0,0,0],
            [0,0,1,-1,0,0,-4,5,0,0,0,0,0,0],
            [0,0,1,-1,0,0,-3,4,0,0,0,0,0,0],
            [0,0,0,2,0,0,0,-1,0,1,0,0,0,0],
            [0,0,0,0,0,0,8,-9,0,0,0,0,0,0],
            [0,0,0,0,0,0,3,-6,0,0,0,0,0,0],
            [0,0,0,0,0,0,1,1,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,3,-5,0,0,0],
            [0,0,0,0,0,0,0,0,0,2,-2,0,0,0],
            [2,0,-2,-2,-2,0,0,-2,0,2,0,0,0,0],
            [1,0,0,0,1,0,-10,3,0,0,0,0,0,0],
            [1,0,0,0,-1,0,-10,3,0,0,0,0,0,0],
            [0,0,2,0,2,0,2,-3,0,0,0,0,0,0],
            [0,0,2,0,2,0,2,-2,0,0,0,0,0,0],
            [0,0,2,0,2,0,-2,3,0,0,0,0,0,0],
            [0,0,2,0,2,0,-2,2,0,0,0,0,0,0],
            [0,0,0,0,2,0,0,0,0,1,0,0,0,0],
            [0,0,0,0,1,0,0,-1,0,2,0,0,0,0],
            [2,0,2,-2,2,0,0,-2,0,3,0,0,0,0],
            [2,0,1,-3,1,0,-6,7,0,0,0,0,0,0],
            [2,0,0,-2,0,0,2,-5,0,0,0,0,0,0],
            [2,0,0,-2,0,0,0,-2,0,5,-5,0,0,0],
            [2,0,0,-2,0,0,0,-2,0,1,5,0,0,0],
            [2,0,0,-2,0,0,0,-2,0,0,5,0,0,0],
            [2,0,0,-2,0,0,0,-2,0,0,2,0,0,0],
            [2,0,0,-2,0,0,-4,4,0,0,0,0,0,0],
            [2,0,-2,0,-2,0,0,5,-9,0,0,0,0,0],
            [2,0,-1,-1,0,0,0,-1,0,3,0,0,0,0],
            [1,0,2,0,2,0,1,-1,0,0,0,0,0,0],
            [1,0,2,0,2,0,0,4,-8,3,0,0,0,0],
            [1,0,2,0,2,0,0,-4,8,-3,0,0,0,0],
            [1,0,2,0,2,0,-1,1,0,0,0,0,0,0],
            [1,0,2,-2,2,0,-3,3,0,0,0,0,0,0],
            [1,0,0,0,0,0,0,1,0,-1,0,0,0,0],
            [1,0,0,0,0,0,0,-2,0,3,0,0,0,0],
            [1,0,0,-2,0,0,0,2,0,-2,0,0,0,0],
            [1,0,-2,-2,-2,0,0,1,0,-1,0,0,0,0],
            [1,0,-1,1,0,0,0,1,0,0,0,0,0,0],
            [1,0,-1,-1,0,0,0,8,-15,0,0,0,0,0],
            [0,0,2,2,2,0,0,2,0,-2,0,0,0,0],
            [0,0,2,-2,1,0,1,-1,0,0,0,0,0,0],
            [0,0,2,-2,1,0,0,-2,0,1,0,0,0,0],
            [0,0,2,-2,1,0,0,-10,15,0,0,0,0,0],
            [0,0,2,-2,0,-1,0,2,0,0,0,0,0,0],
            [0,0,1,-1,2,0,0,-1,0,0,-1,0,0,0],
            [0,0,1,-1,2,0,-3,4,0,0,0,0,0,0],
            [0,0,1,-1,1,0,-4,6,0,0,0,0,0,0],
            [0,0,1,-1,1,0,-1,2,0,0,0,0,0,0],
            [0,0,1,-1,0,0,0,1,0,0,0,0,0,0],
            [0,0,1,-1,0,0,0,-1,0,0,-2,0,0,0],
            [0,0,1,-1,0,0,-2,2,0,0,0,0,0,0],
            [0,0,1,-1,0,0,-1,0,0,0,0,0,0,0],
            [0,0,1,-1,-1,0,-5,7,0,0,0,0,0,0],
            [0,0,0,2,0,0,0,2,0,-2,0,0,0,0],
            [0,0,0,2,0,0,-2,2,0,0,0,0,0,0],
            [0,0,0,0,2,0,-3,5,0,0,0,0,0,0],
            [0,0,0,0,1,0,-1,2,0,0,0,0,0,0],
            [0,0,0,0,0,0,9,-13,0,0,0,0,0,-2],
            [0,0,0,0,0,0,8,-14,0,0,0,0,0,-2],
            [0,0,0,0,0,0,8,-11,0,0,0,0,0,-1],
            [0,0,0,0,0,0,6,-9,0,0,0,0,0,0],
            [0,0,0,0,0,0,6,-8,0,0,0,0,0,0],
            [0,0,0,0,0,0,6,-7,0,0,0,0,0,-1],
            [0,0,0,0,0,0,5,-6,0,0,0,0,0,-2],
            [0,0,0,0,0,0,5,-6,-4,0,0,0,0,-2],
            [0,0,0,0,0,0,5,-4,0,0,0,0,0,2],
            [0,0,0,0,0,0,4,-8,0,0,0,0,0,-2],
            [0,0,0,0,0,0,4,-5,0,0,0,0,0,0],
            [0,0,0,0,0,0,3,-3,0,2,0,0,0,2],
            [0,0,0,0,0,0,3,-1,0,0,0,0,0,0],
            [0,0,0,0,0,0,2,0,0,0,0,0,0,0],
            [0,0,0,0,0,0,1,-1,0,0,0,0,0,-2],
            [0,0,0,0,0,0,0,7,-12,0,0,0,0,-2],
            [0,0,0,0,0,0,0,6,-9,0,0,0,0,-2],
            [0,0,0,0,0,0,0,6,-8,1,5,0,0,2],
            [0,0,0,0,0,0,0,6,-4,0,0,0,0,2],
            [0,0,0,0,0,0,0,6,-10,0,0,0,0,0],
            [0,0,0,0,0,0,0,5,0,-4,0,0,0,2],
            [0,0,0,0,0,0,0,5,-9,0,0,0,0,-1],
            [0,0,0,0,0,0,0,5,-8,3,0,0,0,2],
            [0,0,0,0,0,0,0,5,-7,0,0,0,0,-2],
            [0,0,0,0,0,0,0,5,-6,0,0,0,0,0],
            [0,0,0,0,0,0,0,5,-16,4,5,0,0,-2],
            [0,0,0,0,0,0,0,5,-13,0,0,0,0,-2],
            [0,0,0,0,0,0,0,3,0,-5,0,0,0,-2],
            [0,0,0,0,0,0,0,3,-9,0,0,0,0,-2],
            [0,0,0,0,0,0,0,3,-7,0,0,0,0,-2],
            [0,0,0,0,0,0,0,2,0,2,0,0,0,2],
            [0,0,0,0,0,0,0,2,0,0,-3,0,0,0],
            [0,0,0,0,0,0,0,2,-8,1,5,0,0,-2],
            [0,0,0,0,0,0,0,1,0,1,-5,0,0,0],
            [0,0,0,0,0,0,0,1,0,0,2,0,0,2],
            [0,0,0,0,0,0,0,1,0,0,-3,0,0,0],
            [0,0,0,0,0,0,0,1,0,-3,5,0,0,0],
            [0,0,0,0,0,0,0,1,-3,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,2,-6,3,0,-2],
            [0,0,0,0,0,0,0,0,0,1,-2,0,0,0],
            [0,0,0,0,0,0,0,0,0,0,0,1,0,0],
            [0,0,0,0,0,0,0,0,0,0,0,0,0,2],
            [0,0,0,0,0,0,0,0,1,0,0,0,0,0]]
    
    #指向振幅数组的指针，每个频率一个指针
    NC=[1,21,37,51,65,79,91,103,115,127,
        139,151,163,172,184,196,207,219,231,240,
        252,261,273,285,297,309,318,327,339,351,
        363,372,384,396,405,415,423,435,444,452,
        460,467,474,482,490,498,506,513,521,528,
        536,543,551,559,566,574,582,590,597,605,
        613,620,628,636,644,651,658,666,674,680,
        687,695,702,710,717,725,732,739,746,753,
        760,767,774,782,790,798,805,812,819,826,
        833,840,846,853,860,867,874,881,888,895,
        901,908,914,921,928,934,941,948,955,962,
        969,976,982,989,996,1003,1010,1017,1024,1031,
        1037,1043,1050,1057,1064,1071,1078,1084,1091,1098,
        1104,1112,1118,1124,1131,1138,1145,1151,1157,1164,
        1171,1178,1185,1192,1199,1205,1212,1218,1226,1232,
        1239,1245,1252,1259,1266,1272,1278,1284,1292,1298,
        1304,1310,1316,1323,1329,1335,1341,1347,1353,1359,
        1365,1371,1377,1383,1389,1396,1402,1408,1414,1420,
        1426,1434,1440,1446,1452,1459,1465,1471,1477,1482,
        1488,1493,1499,1504,1509,1514,1520,1527,1532,1538,
        1543,1548,1553,1558,1564,1569,1574,1579,1584,1589,
        1594,1596,1598,1600,1602,1605,1608,1610,1612,1617,
        1619,1623,1625,1627,1629,1632,1634,1640,1642,1644,
        1646,1648,1650,1652,1654,1658,1660,1662,1664,1668,
        1670,1672,1673,1675,1679,1681,1683,1684,1686,1688,
        1690,1693,1695,1697,1701,1703,1705,1707,1709,1711,
        1712,1715,1717,1721,1723,1725,1727,1729,1731,1733,
        1735,1737,1739,1741,1743,1745,1747,1749,1751,1753,
        1755,1757,1759,1761,1762,1764,1766,1768,1769,1771,
        1773,1775,1777,1779,1781,1783,1785,1787,1788,1790,
        1792,1794,1796,1798,1800,1802,1804,1806,1807,1809,
        1811,1815,1817,1819,1821,1823,1825,1827,1829,1831,
        1833,1835,1837,1839,1840,1842,1844,1848,1850,1852,
        1854,1856,1858,1859,1860,1862,1864,1866,1868,1869,
        1871,1873,1875,1877,1879,1881,1883,1885,1887,1889,
        1891,1892,1896,1898,1900,1901,1903,1905,1907,1909,
        1910,1911,1913,1915,1919,1921,1923,1927,1929,1931,
        1933,1935,1937,1939,1943,1945,1947,1948,1949,1951,
        1953,1955,1957,1958,1960,1962,1964,1966,1968,1970,
        1971,1973,1974,1975,1977,1979,1980,1981,1982,1984,
        1986,1988,1990,1992,1994,1995,1997,1999,2001,2003,
        2005,2007,2008,2009,2011,2013,2015,2017,2019,2021,
        2023,2024,2025,2027,2029,2031,2033,2035,2037,2041,
        2043,2045,2046,2047,2049,2051,2053,2055,2056,2057,
        2059,2061,2063,2065,2067,2069,2070,2071,2072,2074,
        2076,2078,2080,2082,2084,2086,2088,2090,2092,2094,
        2095,2096,2097,2099,2101,2105,2106,2107,2108,2109,
        2110,2111,2113,2115,2119,2121,2123,2125,2127,2129,
        2131,2133,2135,2136,2137,2139,2141,2143,2145,2147,
        2149,2151,2153,2155,2157,2159,2161,2163,2165,2167,
        2169,2171,2173,2175,2177,2179,2181,2183,2185,2186,
        2187,2188,2192,2193,2195,2197,2199,2201,2203,2205,
        2207,2209,2211,2213,2217,2219,2221,2223,2225,2227,
        2229,2231,2233,2234,2235,2236,2237,2238,2239,2240,
        2241,2244,2246,2248,2250,2252,2254,2256,2258,2260,
        2262,2264,2266,2268,2270,2272,2274,2276,2278,2280,
        2282,2284,2286,2288,2290,2292,2294,2296,2298,2300,
        2302,2303,2304,2305,2306,2307,2309,2311,2313,2315,
        2317,2319,2321,2323,2325,2327,2329,2331,2333,2335,
        2337,2341,2343,2345,2347,2349,2351,2352,2355,2356,
        2357,2358,2359,2361,2363,2364,2365,2366,2367,2368,
        2369,2370,2371,2372,2373,2374,2376,2378,2380,2382,
        2384,2385,2386,2387,2388,2389,2390,2391,2392,2393,
        2394,2395,2396,2397,2398,2399,2400,2401,2402,2403,
        2404,2405,2406,2407,2408,2409,2410,2411,2412,2413,
        2414,2415,2417,2418,2430,2438,2445,2453,2460,2468,
        2474,2480,2488,2496,2504,2512,2520,2527,2535,2543,
        2550,2558,2566,2574,2580,2588,2596,2604,2612,2619,
        2627,2634,2642,2648,2656,2664,2671,2679,2685,2693,
        2701,2709,2717,2725,2733,2739,2747,2753,2761,2769,
        2777,2785,2793,2801,2809,2817,2825,2833,2841,2848,
        2856,2864,2872,2878,2884,2892,2898,2906,2914,2922,
        2930,2938,2944,2952,2958,2966,2974,2982,2988,2996,
        3001,3009,3017,3025,3032,3039,3045,3052,3059,3067,
        3069,3076,3083,3090,3098,3105,3109,3111,3113,3120,
        3124,3128,3132,3136,3140,3144,3146,3150,3158,3161,
        3165,3166,3168,3172,3176,3180,3182,3185,3189,3193,
        3194,3197,3200,3204,3208,3212,3216,3219,3221,3222,
        3226,3230,3234,3238,3242,3243,3247,3251,3254,3258,
        3262,3266,3270,3274,3275,3279,3283,3287,3289,3293,
        3296,3300,3303,3307,3311,3315,3319,3321,3324,3327,
        3330,3334,3338,3340,3342,3346,3350,3354,3358,3361,
        3365,3369,3373,3377,3381,3385,3389,3393,3394,3398,
        3402,3406,3410,3413,3417,3421,3425,3429,3433,3435,
        3439,3443,3446,3450,3453,3457,3458,3461,3464,3468,
        3472,3476,3478,3481,3485,3489,3493,3497,3501,3505,
        3507,3511,3514,3517,3521,3524,3525,3527,3529,3533,
        3536,3540,3541,3545,3548,3551,3555,3559,3563,3567,
        3569,3570,3574,3576,3578,3582,3586,3590,3593,3596,
        3600,3604,3608,3612,3616,3620,3623,3626,3630,3632,
        3636,3640,3643,3646,3648,3652,3656,3660,3664,3667,
        3669,3671,3675,3679,3683,3687,3689,3693,3694,3695,
        3699,3703,3705,3707,3710,3713,3717,3721,3725,3729,
        3733,3736,3740,3744,3748,3752,3754,3757,3759,3763,
        3767,3770,3773,3777,3779,3783,3786,3790,3794,3798,
        3801,3805,3809,3813,3817,3821,3825,3827,3831,3835,
        3836,3837,3840,3844,3848,3852,3856,3859,3863,3867,
        3869,3871,3875,3879,3883,3887,3890,3894,3898,3901,
        3905,3909,3913,3917,3921,3922,3923,3924,3926,3930,
        3932,3936,3938,3940,3944,3948,3952,3956,3959,3963,
        3965,3969,3973,3977,3979,3981,3982,3986,3989,3993,
        3997,4001,4004,4006,4009,4012,4016,4020,4024,4026,
        4028,4032,4036,4040,4044,4046,4050,4054,4058,4060,
        4062,4063,4064,4068,4071,4075,4077,4081,4083,4087,
        4089,4091,4095,4099,4101,4103,4105,4107,4111,4115,
        4119,4123,4127,4129,4131,4135,4139,4141,4143,4145,
        4149,4153,4157,4161,4165,4169,4173,4177,4180,4183,
        4187,4191,4195,4198,4201,4205,4209,4212,4213,4216,
        4217,4221,4223,4226,4230,4234,4236,4240,4244,4248,
        4252,4256,4258,4262,4264,4266,4268,4270,4272,4276,
        4279,4283,4285,4287,4289,4293,4295,4299,4300,4301,
        4305,4309,4313,4317,4319,4323,4325,4329,4331,4333,
        4335,4337,4341,4345,4349,4351,4353,4357,4361,4365,
        4367,4369,4373,4377,4381,4383,4387,4389,4391,4395,
        4399,4403,4407,4411,4413,4414,4415,4418,4419,4421,
        4423,4427,4429,4431,4433,4435,4437,4439,4443,4446,
        4450,4452,4456,4458,4460,4462,4466,4469,4473,4477,
        4481,4483,4487,4489,4491,4493,4497,4499,4501,4504,
        4506,4510,4513,4514,4515,4518,4521,4522,4525,4526,
        4527,4530,4533,4534,4537,4541,4542,4543,4544,4545,
        4546,4547,4550,4553,4554,4555,4558,4561,4564,4567,
        4568,4571,4574,4575,4578,4581,4582,4585,4586,4588,
        4590,4592,4596,4598,4602,4604,4608,4612,4613,4616,
        4619,4622,4623,4624,4625,4626,4629,4632,4633,4636,
        4639,4640,4641,4642,4643,4644,4645,4648,4649,4650,
        4651,4652,4653,4656,4657,4660,4661,4664,4667,4670,
        4671,4674,4675,4676,4677,4678,4681,4682,4683,4684,
        4687,4688,4689,4692,4693,4696,4697,4700,4701,4702,
        4703,4704,4707,4708,4711,4712,4715,4716,4717,4718,
        4719,4720,4721,4722,4723,4726,4729,4730,4733,4736,
        4737,4740,4741,4742,4745,4746,4749,4752,4753]
    
    #振幅系数(微角秒);使用NC数组进行索引。
    A= [-6844318.44e0,9205236.26e0,1328.67e0,1538.18e0,205833.11e0,
        153041.79e0,-3309.73e0,853.32e0,2037.98e0,-2301.27e0,81.46e0,
        120.56e0,-20.39e0,-15.22e0,1.73e0,-1.61e0,-0.1e0,0.11e0,-0.02e0,
        -0.02e0,-523908.04e0,573033.42e0,-544.75e0,-458.66e0,12814.01e0,
        11714.49e0,198.97e0,-290.91e0,155.74e0,-143.27e0,-2.75e0,-1.03e0,
        -1.27e0,-1.16e0,0e0,-0.01e0,-90552.22e0,97846.69e0,111.23e0,
        137.41e0,2187.91e0,2024.68e0,41.44e0,-51.26e0,26.92e0,-24.46e0,
        -0.46e0,-0.28e0,-0.22e0,-0.2e0,82168.76e0,-89618.24e0,-27.64e0,
        -29.05e0,-2004.36e0,-1837.32e0,-36.07e0,48e0,-24.43e0,22.41e0,
        0.47e0,0.24e0,0.2e0,0.18e0,58707.02e0,7387.02e0,470.05e0,
        -192.4e0,164.33e0,-1312.21e0,-179.73e0,-28.93e0,-17.36e0,-1.83e0,
        -0.5e0,3.57e0,0e0,0.13e0,-20557.78e0,22438.42e0,-20.84e0,
        -17.4e0,501.82e0,459.68e0,59.2e0,-67.3e0,6.08e0,-5.61e0,-1.36e0,
        -1.19e0,28288.28e0,-674.99e0,-34.69e0,35.8e0,-15.07e0,-632.54e0,
        -11.19e0,0.78e0,-8.41e0,0.17e0,0.01e0,0.07e0,-15406.85e0,
        20069.5e0,15.12e0,31.8e0,448.76e0,344.5e0,-5.77e0,1.41e0,4.59e0,
        -5.02e0,0.17e0,0.24e0,-11991.74e0,12902.66e0,32.46e0,36.7e0,
        288.49e0,268.14e0,5.7e0,-7.06e0,3.57e0,-3.23e0,-0.06e0,-0.04e0,
        -8584.95e0,-9592.72e0,4.42e0,-13.2e0,-214.5e0,192.06e0,23.87e0,
        29.83e0,2.54e0,2.4e0,0.6e0,-0.48e0,5095.5e0,-6918.22e0,
        7.19e0,3.92e0,-154.91e0,-113.94e0,2.86e0,-1.04e0,-1.52e0,
        1.73e0,-0.07e0,-0.1e0,-4910.93e0,-5331.13e0,0.76e0,0.4e0,
        -119.21e0,109.81e0,2.16e0,3.2e0,1.46e0,1.33e0,0.04e0,-0.02e0,
        -6245.02e0,-123.48e0,-6.68e0,-8.2e0,-2.76e0,139.64e0,2.71e0,
        0.15e0,1.86e0,2511.85e0,-3323.89e0,1.07e0,-0.9e0,-74.33e0,
        -56.17e0,1.16e0,-0.01e0,-0.75e0,0.83e0,-0.02e0,-0.04e0,
        2307.58e0,3143.98e0,-7.52e0,7.5e0,70.31e0,-51.6e0,1.46e0,0.16e0,
        -0.69e0,-0.79e0,0.02e0,-0.05e0,2372.58e0,2554.51e0,5.93e0,
        -6.6e0,57.12e0,-53.05e0,-0.96e0,-1.24e0,-0.71e0,-0.64e0,-0.01e0,
        -2053.16e0,2636.13e0,5.13e0,7.8e0,58.94e0,45.91e0,-0.42e0,
        -0.12e0,0.61e0,-0.66e0,0.02e0,0.03e0,-1825.49e0,-2423.59e0,
        1.23e0,-2e0,-54.19e0,40.82e0,-1.07e0,-1.02e0,0.54e0,0.61e0,
        -0.04e0,0.04e0,2521.07e0,-122.28e0,-5.97e0,2.9e0,-2.73e0,
        -56.37e0,-0.82e0,0.13e0,-0.75e0,-1534.09e0,1645.01e0,6.29e0,
        6.8e0,36.78e0,34.3e0,0.92e0,-1.25e0,0.46e0,-0.41e0,-0.02e0,
        -0.01e0,1898.27e0,47.7e0,-0.72e0,2.5e0,1.07e0,-42.45e0,-0.94e0,
        0.02e0,-0.56e0,-1292.02e0,-1387e0,0e0,0e0,-31.01e0,28.89e0,
        0.68e0,0e0,0.38e0,0.35e0,-0.01e0,-0.01e0,-1234.96e0,
        1323.81e0,5.21e0,5.9e0,29.6e0,27.61e0,0.74e0,-1.22e0,0.37e0,
        -0.33e0,-0.02e0,-0.01e0,1137.48e0,-1233.89e0,-0.04e0,-0.3e0,
        -27.59e0,-25.43e0,-0.61e0,1e0,-0.34e0,0.31e0,0.01e0,0.01e0,
        -813.13e0,-1075.6e0,0.4e0,0.3e0,-24.05e0,18.18e0,-0.4e0,-0.01e0,
        0.24e0,0.27e0,-0.01e0,0.01e0,1163.22e0,-60.9e0,-2.94e0,1.3e0,
        -1.36e0,-26.01e0,-0.58e0,0.07e0,-0.35e0,1029.7e0,-55.55e0,
        -2.63e0,1.1e0,-1.25e0,-23.02e0,-0.52e0,0.06e0,-0.31e0,
        -556.26e0,852.85e0,3.16e0,-4.48e0,19.06e0,12.44e0,-0.81e0,
        -0.27e0,0.17e0,-0.21e0,0e0,0.02e0,-603.52e0,-800.34e0,
        0.44e0,0.1e0,-17.9e0,13.49e0,-0.08e0,-0.01e0,0.18e0,0.2e0,
        -0.01e0,0.01e0,-628.24e0,684.99e0,-0.64e0,-0.5e0,15.32e0,
        14.05e0,3.18e0,-4.19e0,0.19e0,-0.17e0,-0.09e0,-0.07e0,
        -866.48e0,-16.26e0,0.52e0,-1.3e0,-0.36e0,19.37e0,0.43e0,-0.01e0,
        0.26e0,-512.37e0,695.54e0,-1.47e0,-1.4e0,15.55e0,11.46e0,
        -0.16e0,0.03e0,0.15e0,-0.17e0,0.01e0,0.01e0,506.65e0,
        643.75e0,2.54e0,-2.62e0,14.4e0,-11.33e0,-0.77e0,-0.06e0,-0.15e0,
        -0.16e0,0e0,0.01e0,664.57e0,16.81e0,-0.4e0,1e0,0.38e0,
        -14.86e0,-3.71e0,-0.09e0,-0.2e0,405.91e0,522.11e0,0.99e0,-1.5e0,
        11.67e0,-9.08e0,-0.25e0,-0.02e0,-0.12e0,-0.13e0,-305.78e0,
        326.6e0,1.75e0,1.9e0,7.3e0,6.84e0,0.2e0,-0.04e0,300.99e0,
        -325.03e0,-0.44e0,-0.5e0,-7.27e0,-6.73e0,-1.01e0,0.01e0,0e0,
        0.08e0,0e0,0.02e0,438.51e0,10.47e0,-0.56e0,-0.2e0,0.24e0,
        -9.81e0,-0.24e0,0.01e0,-0.13e0,-264.02e0,335.24e0,0.99e0,1.4e0,
        7.49e0,5.9e0,-0.27e0,-0.02e0,284.09e0,307.03e0,0.32e0,-0.4e0,
        6.87e0,-6.35e0,-0.99e0,-0.01e0,-250.54e0,327.11e0,0.08e0,0.4e0,
        7.31e0,5.6e0,-0.3e0,230.72e0,-304.46e0,0.08e0,-0.1e0,-6.81e0,
        -5.16e0,0.27e0,229.78e0,304.17e0,-0.6e0,0.5e0,6.8e0,-5.14e0,
        0.33e0,0.01e0,256.3e0,-276.81e0,-0.28e0,-0.4e0,-6.19e0,-5.73e0,
        -0.14e0,0.01e0,-212.82e0,269.45e0,0.84e0,1.2e0,6.02e0,4.76e0,
        0.14e0,-0.02e0,196.64e0,272.05e0,-0.84e0,0.9e0,6.08e0,-4.4e0,
        0.35e0,0.02e0,188.95e0,272.22e0,-0.12e0,0.3e0,6.09e0,-4.22e0,
        0.34e0,-292.37e0,-5.1e0,-0.32e0,-0.4e0,-0.11e0,6.54e0,0.14e0,
        0.01e0,161.79e0,-220.67e0,0.24e0,0.1e0,-4.93e0,-3.62e0,-0.08e0,
        261.54e0,-19.94e0,-0.95e0,0.2e0,-0.45e0,-5.85e0,-0.13e0,0.02e0,
        142.16e0,-190.79e0,0.2e0,0.1e0,-4.27e0,-3.18e0,-0.07e0,187.95e0,
        -4.11e0,-0.24e0,0.3e0,-0.09e0,-4.2e0,-0.09e0,0.01e0,0e0,
        0e0,-79.08e0,167.9e0,0.04e0,0e0,3.75e0,1.77e0,121.98e0,
        131.04e0,-0.08e0,0.1e0,2.93e0,-2.73e0,-0.06e0,-172.95e0,
        -8.11e0,-0.4e0,-0.2e0,-0.18e0,3.87e0,0.09e0,0.01e0,
        -160.15e0,-55.3e0,-14.04e0,13.9e0,-1.23e0,3.58e0,0.4e0,0.31e0,
        -115.4e0,123.2e0,0.6e0,0.7e0,2.75e0,2.58e0,0.08e0,-0.01e0,
        -168.26e0,-2e0,0.2e0,-0.2e0,-0.04e0,3.76e0,0.08e0,
        -114.49e0,123.2e0,0.32e0,0.4e0,2.75e0,2.56e0,0.07e0,-0.01e0,
        112.14e0,120.7e0,0.28e0,-0.3e0,2.7e0,-2.51e0,-0.07e0,-0.01e0,
        161.34e0,4.03e0,0.2e0,0.2e0,0.09e0,-3.61e0,-0.08e0,91.31e0,
        126.64e0,-0.4e0,0.4e0,2.83e0,-2.04e0,-0.04e0,0.01e0,105.29e0,
        112.9e0,0.44e0,-0.5e0,2.52e0,-2.35e0,-0.07e0,-0.01e0,98.69e0,
        -106.2e0,-0.28e0,-0.3e0,-2.37e0,-2.21e0,-0.06e0,0.01e0,86.74e0,
        -112.94e0,-0.08e0,-0.2e0,-2.53e0,-1.94e0,-0.05e0,-134.81e0,
        3.51e0,0.2e0,-0.2e0,0.08e0,3.01e0,0.07e0,79.03e0,107.31e0,
        -0.24e0,0.2e0,2.4e0,-1.77e0,-0.04e0,0.01e0,132.81e0,
        -10.77e0,-0.52e0,0.1e0,-0.24e0,-2.97e0,-0.07e0,0.01e0,
        -130.31e0,-0.9e0,0.04e0,0e0,0e0,2.91e0,-78.56e0,85.32e0,
        0e0,0e0,1.91e0,1.76e0,0.04e0,0e0,0e0,-41.53e0,
        89.1e0,0.02e0,0e0,1.99e0,0.93e0,66.03e0,-71e0,-0.2e0,
        -0.2e0,-1.59e0,-1.48e0,-0.04e0,60.5e0,64.7e0,0.36e0,-0.4e0,
        1.45e0,-1.35e0,-0.04e0,-0.01e0,-52.27e0,-70.01e0,0e0,0e0,
        -1.57e0,1.17e0,0.03e0,-52.95e0,66.29e0,0.32e0,0.4e0,1.48e0,
        1.18e0,0.04e0,-0.01e0,51.02e0,67.25e0,0e0,0e0,1.5e0,
        -1.14e0,-0.03e0,-55.66e0,-60.92e0,0.16e0,-0.2e0,-1.36e0,1.24e0,
        0.03e0,-54.81e0,-59.2e0,-0.08e0,0.2e0,-1.32e0,1.23e0,0.03e0,
        51.32e0,-55.6e0,0e0,0e0,-1.24e0,-1.15e0,-0.03e0,48.29e0,
        51.8e0,0.2e0,-0.2e0,1.16e0,-1.08e0,-0.03e0,-45.59e0,-49e0,
        -0.12e0,0.1e0,-1.1e0,1.02e0,0.03e0,40.54e0,-52.69e0,-0.04e0,
        -0.1e0,-1.18e0,-0.91e0,-0.02e0,-40.58e0,-49.51e0,-1e0,1e0,
        -1.11e0,0.91e0,0.04e0,0.02e0,-43.76e0,46.5e0,0.36e0,0.4e0,
        1.04e0,0.98e0,0.03e0,-0.01e0,62.65e0,-5e0,-0.24e0,0e0,
        -0.11e0,-1.4e0,-0.03e0,0.01e0,-38.57e0,49.59e0,0.08e0,0.1e0,
        1.11e0,0.86e0,0.02e0,-33.22e0,-44.04e0,0.08e0,-0.1e0,-0.98e0,
        0.74e0,0.02e0,37.15e0,-39.9e0,-0.12e0,-0.1e0,-0.89e0,-0.83e0,
        -0.02e0,36.68e0,-39.5e0,-0.04e0,-0.1e0,-0.88e0,-0.82e0,-0.02e0,
        -53.22e0,-3.91e0,-0.2e0,0e0,-0.09e0,1.19e0,0.03e0,32.43e0,
        -42.19e0,-0.04e0,-0.1e0,-0.94e0,-0.73e0,-0.02e0,-51e0,-2.3e0,
        -0.12e0,-0.1e0,0e0,1.14e0,-29.53e0,-39.11e0,0.04e0,0e0,
        -0.87e0,0.66e0,0.02e0,28.5e0,-38.92e0,-0.08e0,-0.1e0,-0.87e0,
        -0.64e0,-0.02e0,26.54e0,36.95e0,-0.12e0,0.1e0,0.83e0,-0.59e0,
        -0.01e0,26.54e0,34.59e0,0.04e0,-0.1e0,0.77e0,-0.59e0,-0.02e0,
        28.35e0,-32.55e0,-0.16e0,0.2e0,-0.73e0,-0.63e0,-0.01e0,-28e0,
        30.4e0,0e0,0e0,0.68e0,0.63e0,0.01e0,-27.61e0,29.4e0,
        0.2e0,0.2e0,0.66e0,0.62e0,0.02e0,40.33e0,0.4e0,-0.04e0,
        0.1e0,0e0,-0.9e0,-23.28e0,31.61e0,-0.08e0,-0.1e0,0.71e0,
        0.52e0,0.01e0,37.75e0,0.8e0,0.04e0,0.1e0,0e0,-0.84e0,
        23.66e0,25.8e0,0e0,0e0,0.58e0,-0.53e0,-0.01e0,21.01e0,
        -27.91e0,0e0,0e0,-0.62e0,-0.47e0,-0.01e0,-34.81e0,2.89e0,
        0.04e0,0e0,0e0,0.78e0,-23.49e0,-25.31e0,0e0,0e0,
        -0.57e0,0.53e0,0.01e0,-23.47e0,25.2e0,0.16e0,0.2e0,0.56e0,
        0.52e0,0.02e0,19.58e0,27.5e0,-0.12e0,0.1e0,0.62e0,-0.44e0,
        -0.01e0,-22.67e0,-24.4e0,-0.08e0,0.1e0,-0.55e0,0.51e0,0.01e0,
        -19.97e0,25e0,0.12e0,0.2e0,0.56e0,0.45e0,0.01e0,21.28e0,
        -22.8e0,-0.08e0,-0.1e0,-0.51e0,-0.48e0,-0.01e0,-30.47e0,0.91e0,
        0.04e0,0e0,0e0,0.68e0,18.58e0,24e0,0.04e0,-0.1e0,
        0.54e0,-0.42e0,-0.01e0,-18.02e0,24.4e0,-0.04e0,-0.1e0,0.55e0,
        0.4e0,0.01e0,17.74e0,22.5e0,0.08e0,-0.1e0,0.5e0,-0.4e0,
        -0.01e0,-19.41e0,20.7e0,0.08e0,0.1e0,0.46e0,0.43e0,0.01e0,
        -18.64e0,20.11e0,0e0,0e0,0.45e0,0.42e0,0.01e0,-16.75e0,
        21.6e0,0.04e0,0.1e0,0.48e0,0.37e0,0.01e0,-18.42e0,-20e0,
        0e0,0e0,-0.45e0,0.41e0,0.01e0,-26.77e0,1.41e0,0.08e0,
        0e0,0e0,0.6e0,-26.17e0,-0.19e0,0e0,0e0,0e0,
        0.59e0,-15.52e0,20.51e0,0e0,0e0,0.46e0,0.35e0,0.01e0,
        -25.42e0,-1.91e0,-0.08e0,0e0,-0.04e0,0.57e0,0.45e0,-17.42e0,
        18.1e0,0e0,0e0,0.4e0,0.39e0,0.01e0,16.39e0,-17.6e0,
        -0.08e0,-0.1e0,-0.39e0,-0.37e0,-0.01e0,-14.37e0,18.91e0,0e0,
        0e0,0.42e0,0.32e0,0.01e0,23.39e0,-2.4e0,-0.12e0,0e0,
        0e0,-0.52e0,14.32e0,-18.5e0,-0.04e0,-0.1e0,-0.41e0,-0.32e0,
        -0.01e0,15.69e0,17.08e0,0e0,0e0,0.38e0,-0.35e0,-0.01e0,
        -22.99e0,0.5e0,0.04e0,0e0,0e0,0.51e0,0e0,0e0,
        14.47e0,-17.6e0,-0.01e0,0e0,-0.39e0,-0.32e0,-13.33e0,18.4e0,
        -0.04e0,-0.1e0,0.41e0,0.3e0,22.47e0,-0.6e0,-0.04e0,0e0,
        0e0,-0.5e0,-12.78e0,-17.41e0,0.04e0,0e0,-0.39e0,0.29e0,
        0.01e0,-14.1e0,-15.31e0,0.04e0,0e0,-0.34e0,0.32e0,0.01e0,
        11.98e0,16.21e0,-0.04e0,0e0,0.36e0,-0.27e0,-0.01e0,19.65e0,
        -1.9e0,-0.08e0,0e0,0e0,-0.44e0,19.61e0,-1.5e0,-0.08e0,
        0e0,0e0,-0.44e0,13.41e0,-14.3e0,-0.04e0,-0.1e0,-0.32e0,
        -0.3e0,-0.01e0,-13.29e0,14.4e0,0e0,0e0,0.32e0,0.3e0,
        0.01e0,11.14e0,-14.4e0,-0.04e0,0e0,-0.32e0,-0.25e0,-0.01e0,
        12.24e0,-13.38e0,0.04e0,0e0,-0.3e0,-0.27e0,-0.01e0,10.07e0,
        -13.81e0,0.04e0,0e0,-0.31e0,-0.23e0,-0.01e0,10.46e0,13.1e0,
        0.08e0,-0.1e0,0.29e0,-0.23e0,-0.01e0,16.55e0,-1.71e0,-0.08e0,
        0e0,0e0,-0.37e0,9.75e0,-12.8e0,0e0,0e0,-0.29e0,
        -0.22e0,-0.01e0,9.11e0,12.8e0,0e0,0e0,0.29e0,-0.2e0,
        0e0,0e0,-6.44e0,-13.8e0,0e0,0e0,-0.31e0,0.14e0,
        -9.19e0,-12e0,0e0,0e0,-0.27e0,0.21e0,-10.3e0,10.9e0,
        0.08e0,0.1e0,0.24e0,0.23e0,0.01e0,14.92e0,-0.8e0,-0.04e0,
        0e0,0e0,-0.33e0,10.02e0,-10.8e0,0e0,0e0,-0.24e0,
        -0.22e0,-0.01e0,-9.75e0,10.4e0,0.04e0,0e0,0.23e0,0.22e0,
        0.01e0,9.67e0,-10.4e0,-0.04e0,0e0,-0.23e0,-0.22e0,-0.01e0,
        -8.28e0,-11.2e0,0.04e0,0e0,-0.25e0,0.19e0,13.32e0,-1.41e0,
        -0.08e0,0e0,0e0,-0.3e0,8.27e0,10.5e0,0.04e0,0e0,
        0.23e0,-0.19e0,0e0,0e0,13.13e0,0e0,0e0,0e0,
        0e0,-0.29e0,-12.93e0,0.7e0,0.04e0,0e0,0e0,0.29e0,
        7.91e0,-10.2e0,0e0,0e0,-0.23e0,-0.18e0,-7.84e0,-10e0,
        -0.04e0,0e0,-0.22e0,0.18e0,7.44e0,9.6e0,0e0,0e0,
        0.21e0,-0.17e0,-7.64e0,9.4e0,0.08e0,0.1e0,0.21e0,0.17e0,
        0.01e0,-11.38e0,0.6e0,0.04e0,0e0,0e0,0.25e0,-7.48e0,
        8.3e0,0e0,0e0,0.19e0,0.17e0,-10.98e0,-0.2e0,0e0,
        0e0,0e0,0.25e0,10.98e0,0.2e0,0e0,0e0,0e0,
        -0.25e0,7.4e0,-7.9e0,-0.04e0,0e0,-0.18e0,-0.17e0,-6.09e0,
        8.4e0,-0.04e0,0e0,0.19e0,0.14e0,-6.94e0,-7.49e0,0e0,
        0e0,-0.17e0,0.16e0,6.92e0,7.5e0,0.04e0,0e0,0.17e0,
        -0.15e0,6.2e0,8.09e0,0e0,0e0,0.18e0,-0.14e0,-6.12e0,
        7.8e0,0.04e0,0e0,0.17e0,0.14e0,5.85e0,-7.5e0,0e0,
        0e0,-0.17e0,-0.13e0,-6.48e0,6.9e0,0.08e0,0.1e0,0.15e0,
        0.14e0,0.01e0,6.32e0,6.9e0,0e0,0e0,0.15e0,-0.14e0,
        5.61e0,-7.2e0,0e0,0e0,-0.16e0,-0.13e0,9.07e0,0e0,
        0e0,0e0,0e0,-0.2e0,5.25e0,6.9e0,0e0,0e0,
        0.15e0,-0.12e0,-8.47e0,-0.4e0,0e0,0e0,0e0,0.19e0,
        6.32e0,-5.39e0,-1.11e0,1.1e0,-0.12e0,-0.14e0,0.02e0,0.02e0,
        5.73e0,-6.1e0,-0.04e0,0e0,-0.14e0,-0.13e0,4.7e0,6.6e0,
        -0.04e0,0e0,0.15e0,-0.11e0,-4.9e0,-6.4e0,0e0,0e0,
        -0.14e0,0.11e0,-5.33e0,5.6e0,0.04e0,0.1e0,0.13e0,0.12e0,
        0.01e0,-4.81e0,6e0,0.04e0,0e0,0.13e0,0.11e0,5.13e0,
        5.5e0,0.04e0,0e0,0.12e0,-0.11e0,4.5e0,5.9e0,0e0,
        0e0,0.13e0,-0.1e0,-4.22e0,6.1e0,0e0,0e0,0.14e0,
        -4.53e0,5.7e0,0e0,0e0,0.13e0,0.1e0,4.18e0,5.7e0,
        0e0,0e0,0.13e0,-4.75e0,-5.19e0,0e0,0e0,-0.12e0,
        0.11e0,-4.06e0,5.6e0,0e0,0e0,0.13e0,-3.98e0,5.6e0,
        -0.04e0,0e0,0.13e0,4.02e0,-5.4e0,0e0,0e0,-0.12e0,
        4.49e0,-4.9e0,-0.04e0,0e0,-0.11e0,-0.1e0,-3.62e0,-5.4e0,
        -0.16e0,0.2e0,-0.12e0,0e0,0.01e0,4.38e0,4.8e0,0e0,
        0e0,0.11e0,-6.4e0,-0.1e0,0e0,0e0,0e0,0.14e0,
        -3.98e0,5e0,0.04e0,0e0,0.11e0,-3.82e0,-5e0,0e0,
        0e0,-0.11e0,-3.71e0,5.07e0,0e0,0e0,0.11e0,4.14e0,
        4.4e0,0e0,0e0,0.1e0,-6.01e0,-0.5e0,-0.04e0,0e0,
        0e0,0.13e0,-4.04e0,4.39e0,0e0,0e0,0.1e0,3.45e0,
        -4.72e0,0e0,0e0,-0.11e0,3.31e0,4.71e0,0e0,0e0,
        0.11e0,3.26e0,-4.5e0,0e0,0e0,-0.1e0,-3.26e0,-4.5e0,
        0e0,0e0,-0.1e0,-3.34e0,-4.4e0,0e0,0e0,-0.1e0,
        -3.74e0,-4e0,3.7e0,4e0,3.34e0,-4.3e0,3.3e0,-4.3e0,
        -3.66e0,3.9e0,0.04e0,3.66e0,3.9e0,0.04e0,-3.62e0,-3.9e0,
        -3.61e0,3.9e0,-0.2e0,5.3e0,0e0,0e0,0.12e0,3.06e0,
        4.3e0,3.3e0,4e0,0.4e0,0.2e0,3.1e0,4.1e0,-3.06e0,
        3.9e0,-3.3e0,-3.6e0,-3.3e0,3.36e0,0.01e0,3.14e0,3.4e0,
        -4.57e0,-0.2e0,0e0,0e0,0e0,0.1e0,-2.7e0,-3.6e0,
        2.94e0,-3.2e0,-2.9e0,3.2e0,2.47e0,-3.4e0,2.55e0,-3.3e0,
        2.8e0,-3.08e0,2.51e0,3.3e0,-4.1e0,0.3e0,-0.12e0,-0.1e0,
        4.1e0,0.2e0,-2.74e0,3e0,2.46e0,3.23e0,-3.66e0,1.2e0,
        -0.2e0,0.2e0,3.74e0,-0.4e0,-2.51e0,-2.8e0,-3.74e0,2.27e0,
        -2.9e0,0e0,0e0,-2.5e0,2.7e0,-2.51e0,2.6e0,-3.5e0,
        0.2e0,3.38e0,-2.22e0,-2.5e0,3.26e0,-0.4e0,1.95e0,-2.6e0,
        3.22e0,-0.4e0,-0.04e0,-1.79e0,-2.6e0,1.91e0,2.5e0,0.74e0,
        3.05e0,-0.04e0,0.08e0,2.11e0,-2.3e0,-2.11e0,2.2e0,-1.87e0,
        -2.4e0,2.03e0,-2.2e0,-2.03e0,2.2e0,2.98e0,0e0,0e0,
        2.98e0,-1.71e0,2.4e0,2.94e0,-0.1e0,-0.12e0,0.1e0,1.67e0,
        2.4e0,-1.79e0,2.3e0,-1.79e0,2.2e0,-1.67e0,2.2e0,1.79e0,
        -2e0,1.87e0,-1.9e0,1.63e0,-2.1e0,-1.59e0,2.1e0,1.55e0,
        -2.1e0,-1.55e0,2.1e0,-2.59e0,-0.2e0,-1.75e0,-1.9e0,-1.75e0,
        1.9e0,-1.83e0,-1.8e0,1.51e0,2e0,-1.51e0,-2e0,1.71e0,
        1.8e0,1.31e0,2.1e0,-1.43e0,2e0,1.43e0,2e0,-2.43e0,
        -1.51e0,1.9e0,-1.47e0,1.9e0,2.39e0,0.2e0,-2.39e0,1.39e0,
        1.9e0,1.39e0,-1.8e0,1.47e0,-1.6e0,1.47e0,-1.6e0,1.43e0,
        -1.5e0,-1.31e0,1.6e0,1.27e0,-1.6e0,-1.27e0,1.6e0,1.27e0,
        -1.6e0,2.03e0,1.35e0,1.5e0,-1.39e0,-1.4e0,1.95e0,-0.2e0,
        -1.27e0,1.49e0,1.19e0,1.5e0,1.27e0,1.4e0,1.15e0,1.5e0,
        1.87e0,-0.1e0,-1.12e0,-1.5e0,1.87e0,-1.11e0,-1.5e0,-1.11e0,
        -1.5e0,0e0,0e0,1.19e0,1.4e0,1.27e0,-1.3e0,-1.27e0,
        -1.3e0,-1.15e0,1.4e0,-1.23e0,1.3e0,-1.23e0,-1.3e0,1.22e0,
        -1.29e0,1.07e0,-1.4e0,1.75e0,-0.2e0,-1.03e0,-1.4e0,-1.07e0,
        1.2e0,-1.03e0,1.15e0,1.07e0,1.1e0,1.51e0,-1.03e0,1.1e0,
        1.03e0,-1.1e0,0e0,0e0,-1.03e0,-1.1e0,0.91e0,-1.2e0,
        -0.88e0,-1.2e0,-0.88e0,1.2e0,-0.95e0,1.1e0,-0.95e0,-1.1e0,
        1.43e0,-1.39e0,0.95e0,-1e0,-0.95e0,1e0,-0.8e0,1.1e0,
        0.91e0,-1e0,-1.35e0,0.88e0,1e0,-0.83e0,1e0,-0.91e0,
        0.9e0,0.91e0,0.9e0,0.88e0,-0.9e0,-0.76e0,-1e0,-0.76e0,
        1e0,0.76e0,1e0,-0.72e0,1e0,0.84e0,-0.9e0,0.84e0,
        0.9e0,1.23e0,0e0,0e0,-0.52e0,-1.1e0,-0.68e0,1e0,
        1.19e0,-0.2e0,1.19e0,0.76e0,0.9e0,1.15e0,-0.1e0,1.15e0,
        -0.1e0,0.72e0,-0.9e0,-1.15e0,-1.15e0,0.68e0,0.9e0,-0.68e0,
        0.9e0,-1.11e0,0e0,0e0,0.2e0,0.79e0,0.8e0,-1.11e0,
        -0.1e0,0e0,0e0,-0.48e0,-1e0,-0.76e0,-0.8e0,-0.72e0,
        -0.8e0,-1.07e0,-0.1e0,0.64e0,0.8e0,-0.64e0,-0.8e0,0.64e0,
        0.8e0,0.4e0,0.6e0,0.52e0,-0.5e0,-0.6e0,-0.8e0,-0.71e0,
        0.7e0,-0.99e0,0.99e0,0.56e0,0.8e0,-0.56e0,0.8e0,0.68e0,
        -0.7e0,0.68e0,0.7e0,-0.95e0,-0.64e0,0.7e0,0.64e0,0.7e0,
        -0.6e0,0.7e0,-0.6e0,-0.7e0,-0.91e0,-0.1e0,-0.51e0,0.76e0,
        -0.91e0,-0.56e0,0.7e0,0.88e0,0.88e0,-0.63e0,-0.6e0,0.55e0,
        -0.6e0,-0.8e0,0.8e0,-0.8e0,-0.52e0,0.6e0,0.52e0,0.6e0,
        0.52e0,-0.6e0,-0.48e0,0.6e0,0.48e0,0.6e0,0.48e0,0.6e0,
        -0.76e0,0.44e0,-0.6e0,0.52e0,-0.5e0,-0.52e0,0.5e0,0.4e0,
        0.6e0,-0.4e0,-0.6e0,0.4e0,-0.6e0,0.72e0,-0.72e0,-0.51e0,
        -0.5e0,-0.48e0,0.5e0,0.48e0,-0.5e0,-0.48e0,0.5e0,-0.48e0,
        0.5e0,0.48e0,-0.5e0,-0.48e0,-0.5e0,-0.68e0,-0.68e0,0.44e0,
        0.5e0,-0.64e0,-0.1e0,-0.64e0,-0.1e0,-0.4e0,0.5e0,0.4e0,
        0.5e0,0.4e0,0.5e0,0e0,0e0,-0.4e0,-0.5e0,-0.36e0,
        -0.5e0,0.36e0,-0.5e0,0.6e0,-0.6e0,0.4e0,-0.4e0,0.4e0,
        0.4e0,-0.4e0,0.4e0,-0.4e0,0.4e0,-0.56e0,-0.56e0,0.36e0,
        -0.4e0,-0.36e0,0.4e0,0.36e0,-0.4e0,-0.36e0,-0.4e0,0.36e0,
        0.4e0,0.36e0,0.4e0,-0.52e0,0.52e0,0.52e0,0.32e0,0.4e0,
        -0.32e0,0.4e0,-0.32e0,0.4e0,-0.32e0,0.4e0,0.32e0,-0.4e0,
        -0.32e0,-0.4e0,0.32e0,-0.4e0,0.28e0,-0.4e0,-0.28e0,0.4e0,
        0.28e0,-0.4e0,0.28e0,0.4e0,0.48e0,-0.48e0,0.48e0,0.36e0,
        -0.3e0,-0.36e0,-0.3e0,0e0,0e0,0.2e0,0.4e0,-0.44e0,
        0.44e0,-0.44e0,-0.44e0,-0.44e0,-0.44e0,0.32e0,-0.3e0,0.32e0,
        0.3e0,0.24e0,0.3e0,-0.12e0,-0.1e0,-0.28e0,0.3e0,0.28e0,
        0.3e0,0.28e0,0.3e0,0.28e0,-0.3e0,0.28e0,-0.3e0,0.28e0,
        -0.3e0,0.28e0,0.3e0,-0.28e0,0.3e0,0.4e0,0.4e0,-0.24e0,
        0.3e0,0.24e0,-0.3e0,0.24e0,-0.3e0,-0.24e0,-0.3e0,0.24e0,
        0.3e0,0.24e0,-0.3e0,-0.24e0,0.3e0,0.24e0,-0.3e0,-0.24e0,
        -0.3e0,0.24e0,-0.3e0,0.24e0,0.3e0,-0.24e0,0.3e0,-0.24e0,
        0.3e0,0.2e0,-0.3e0,0.2e0,-0.3e0,0.2e0,-0.3e0,0.2e0,
        0.3e0,0.2e0,-0.3e0,0.2e0,-0.3e0,0.2e0,0.3e0,0.2e0,
        0.3e0,-0.2e0,-0.3e0,0.2e0,-0.3e0,0.2e0,-0.3e0,-0.36e0,
        -0.36e0,-0.36e0,-0.04e0,0.3e0,0.12e0,-0.1e0,-0.32e0,-0.24e0,
        0.2e0,0.24e0,0.2e0,0.2e0,-0.2e0,-0.2e0,-0.2e0,-0.2e0,
        -0.2e0,0.2e0,0.2e0,0.2e0,-0.2e0,0.2e0,0.2e0,0.2e0,
        0.2e0,-0.2e0,-0.2e0,0e0,0e0,-0.2e0,-0.2e0,-0.2e0,
        0.2e0,-0.2e0,0.2e0,0.2e0,-0.2e0,-0.2e0,-0.2e0,0.2e0,
        0.2e0,0.2e0,0.2e0,0.2e0,-0.2e0,0.2e0,-0.2e0,0.28e0,
        0.28e0,0.28e0,0.28e0,0.28e0,0.28e0,-0.28e0,0.28e0,0.12e0,
        0e0,0.24e0,0.16e0,-0.2e0,0.16e0,-0.2e0,0.16e0,-0.2e0,
        0.16e0,0.2e0,-0.16e0,0.2e0,0.16e0,0.2e0,-0.16e0,0.2e0,
        -0.16e0,0.2e0,-0.16e0,0.2e0,0.16e0,-0.2e0,0.16e0,0.2e0,
        0.16e0,-0.2e0,-0.16e0,0.2e0,-0.16e0,-0.2e0,-0.16e0,0.2e0,
        0.16e0,0.2e0,0.16e0,-0.2e0,0.16e0,-0.2e0,0.16e0,0.2e0,
        0.16e0,0.2e0,0.16e0,0.2e0,-0.16e0,-0.2e0,0.16e0,0.2e0,
        -0.16e0,0.2e0,0.16e0,0.2e0,-0.16e0,-0.2e0,0.16e0,-0.2e0,
        0.16e0,-0.2e0,-0.16e0,-0.2e0,0.24e0,-0.24e0,-0.24e0,0.24e0,
        0.24e0,0.12e0,0.2e0,0.12e0,0.2e0,-0.12e0,-0.2e0,0.12e0,
        -0.2e0,0.12e0,-0.2e0,-0.12e0,0.2e0,-0.12e0,0.2e0,-0.12e0,
        -0.2e0,0.12e0,0.2e0,0.12e0,0.2e0,0.12e0,-0.2e0,-0.12e0,
        0.2e0,0.12e0,-0.2e0,-0.12e0,0.2e0,0.12e0,0.2e0,0e0,
        0e0,-0.12e0,0.2e0,-0.12e0,0.2e0,0.12e0,-0.2e0,-0.12e0,
        0.2e0,0.12e0,0.2e0,0e0,-0.21e0,-0.2e0,0e0,0e0,
        0.2e0,-0.2e0,-0.2e0,-0.2e0,0.2e0,-0.16e0,-0.1e0,0e0,
        0.17e0,0.16e0,0.16e0,0.16e0,0.16e0,-0.16e0,0.16e0,0.16e0,
        -0.16e0,0.16e0,-0.16e0,0.16e0,0.12e0,0.1e0,0.12e0,-0.1e0,
        -0.12e0,0.1e0,-0.12e0,0.1e0,0.12e0,-0.1e0,-0.12e0,0.12e0,
        -0.12e0,0.12e0,-0.12e0,0.12e0,-0.12e0,-0.12e0,-0.12e0,-0.12e0,
        -0.12e0,-0.12e0,-0.12e0,0.12e0,0.12e0,0.12e0,0.12e0,-0.12e0,
        -0.12e0,0.12e0,0.12e0,0.12e0,-0.12e0,0.12e0,-0.12e0,-0.12e0,
        -0.12e0,0.12e0,-0.12e0,-0.12e0,0.12e0,0e0,0.11e0,0.11e0,
        -122.67e0,164.7e0,203.78e0,273.5e0,3.58e0,2.74e0,6.18e0,-4.56e0,
        0e0,-0.04e0,0e0,-0.07e0,57.44e0,-77.1e0,95.82e0,128.6e0,
        -1.77e0,-1.28e0,2.85e0,-2.14e0,82.14e0,89.5e0,0e0,0e0,
        2e0,-1.84e0,-0.04e0,47.73e0,-64.1e0,23.79e0,31.9e0,-1.45e0,
        -1.07e0,0.69e0,-0.53e0,-46.38e0,50.5e0,0e0,0e0,1.13e0,
        1.04e0,0.02e0,-18.38e0,0e0,63.8e0,0e0,0e0,0.41e0,
        0e0,-1.43e0,59.07e0,0e0,0e0,0e0,0e0,-1.32e0,
        57.28e0,0e0,0e0,0e0,0e0,-1.28e0,-48.65e0,0e0,
        -1.15e0,0e0,0e0,1.09e0,0e0,0.03e0,-18.3e0,24.6e0,
        -17.3e0,-23.2e0,0.56e0,0.41e0,-0.51e0,0.39e0,-16.91e0,26.9e0,
        8.43e0,13.3e0,0.6e0,0.38e0,0.31e0,-0.19e0,1.23e0,-1.7e0,
        -19.13e0,-25.7e0,-0.03e0,-0.03e0,-0.58e0,0.43e0,-0.72e0,0.9e0,
        -17.34e0,-23.3e0,0.03e0,0.02e0,-0.52e0,0.39e0,-19.49e0,-21.3e0,
        0e0,0e0,-0.48e0,0.44e0,0.01e0,20.57e0,-20.1e0,0.64e0,
        0.7e0,-0.45e0,-0.46e0,0e0,-0.01e0,4.89e0,5.9e0,-16.55e0,
        19.9e0,0.14e0,-0.11e0,0.44e0,0.37e0,18.22e0,19.8e0,0e0,
        0e0,0.44e0,-0.41e0,-0.01e0,4.89e0,-5.3e0,-16.51e0,-18e0,
        -0.11e0,-0.11e0,-0.41e0,0.37e0,-17.86e0,0e0,17.1e0,0e0,
        0e0,0.4e0,0e0,-0.38e0,0.32e0,0e0,24.42e0,0e0,
        0e0,-0.01e0,0e0,-0.55e0,-23.79e0,0e0,0e0,0e0,
        0e0,0.53e0,14.72e0,-16e0,-0.32e0,0e0,-0.36e0,-0.33e0,
        -0.01e0,0.01e0,3.34e0,-4.5e0,11.86e0,15.9e0,-0.11e0,-0.07e0,
        0.35e0,-0.27e0,-3.26e0,4.4e0,11.62e0,15.6e0,0.09e0,0.07e0,
        0.35e0,-0.26e0,-19.53e0,0e0,5.09e0,0e0,0e0,0.44e0,
        0e0,-0.11e0,-13.48e0,14.7e0,0e0,0e0,0.33e0,0.3e0,
        0.01e0,10.86e0,-14.6e0,3.18e0,4.3e0,-0.33e0,-0.24e0,0.09e0,
        -0.07e0,-11.3e0,-15.1e0,0e0,0e0,-0.34e0,0.25e0,0.01e0,
        2.03e0,-2.7e0,10.82e0,14.5e0,-0.07e0,-0.05e0,0.32e0,-0.24e0,
        17.46e0,0e0,0e0,0e0,0e0,-0.39e0,16.43e0,0e0,
        0.52e0,0e0,0e0,-0.37e0,0e0,-0.01e0,9.35e0,0e0,
        13.29e0,0e0,0e0,-0.21e0,0e0,-0.3e0,-10.42e0,11.4e0,
        0e0,0e0,0.25e0,0.23e0,0.01e0,0.44e0,0.5e0,-10.38e0,
        11.3e0,0.02e0,-0.01e0,0.25e0,0.23e0,-14.64e0,0e0,0e0,
        0e0,0e0,0.33e0,0.56e0,0.8e0,-8.67e0,11.7e0,0.02e0,
        -0.01e0,0.26e0,0.19e0,13.88e0,0e0,-2.47e0,0e0,0e0,
        -0.31e0,0e0,0.06e0,-1.99e0,2.7e0,7.72e0,10.3e0,0.06e0,
        0.04e0,0.23e0,-0.17e0,-0.2e0,0e0,13.05e0,0e0,0e0,
        0e0,0e0,-0.29e0,6.92e0,-9.3e0,3.34e0,4.5e0,-0.21e0,
        -0.15e0,0.1e0,-0.07e0,-6.6e0,0e0,10.7e0,0e0,0e0,
        0.15e0,0e0,-0.24e0,-8.04e0,-8.7e0,0e0,0e0,-0.19e0,
        0.18e0,-10.58e0,0e0,-3.1e0,0e0,0e0,0.24e0,0e0,
        0.07e0,-7.32e0,8e0,-0.12e0,-0.1e0,0.18e0,0.16e0,1.63e0,
        1.7e0,6.96e0,-7.6e0,0.03e0,-0.04e0,-0.17e0,-0.16e0,-3.62e0,
        0e0,9.86e0,0e0,0e0,0.08e0,0e0,-0.22e0,0.2e0,
        -0.2e0,-6.88e0,-7.5e0,0e0,0e0,-0.17e0,0.15e0,-8.99e0,
        0e0,4.02e0,0e0,0e0,0.2e0,0e0,-0.09e0,-1.07e0,
        1.4e0,-5.69e0,-7.7e0,0.03e0,0.02e0,-0.17e0,0.13e0,6.48e0,
        -7.2e0,-0.48e0,-0.5e0,-0.16e0,-0.14e0,-0.01e0,0.01e0,5.57e0,
        -7.5e0,1.07e0,1.4e0,-0.17e0,-0.12e0,0.03e0,-0.02e0,8.71e0,
        0e0,3.54e0,0e0,0e0,-0.19e0,0e0,-0.08e0,0.4e0,
        0e0,9.27e0,0e0,0e0,-0.01e0,0e0,-0.21e0,-6.13e0,
        6.7e0,-1.19e0,-1.3e0,0.15e0,0.14e0,-0.03e0,0.03e0,5.21e0,
        -5.7e0,-2.51e0,-2.6e0,-0.13e0,-0.12e0,-0.06e0,0.06e0,5.69e0,
        -6.2e0,-0.12e0,-0.1e0,-0.14e0,-0.13e0,-0.01e0,2.03e0,-2.7e0,
        4.53e0,6.1e0,-0.06e0,-0.05e0,0.14e0,-0.1e0,5.01e0,5.5e0,
        -2.51e0,2.7e0,0.12e0,-0.11e0,0.06e0,0.06e0,-1.91e0,2.6e0,
        -4.38e0,-5.9e0,0.06e0,0.04e0,-0.13e0,0.1e0,4.65e0,-6.3e0,
        0e0,0e0,-0.14e0,-0.1e0,-5.29e0,5.7e0,0e0,0e0,
        0.13e0,0.12e0,-2.23e0,-4e0,-4.65e0,4.2e0,-0.09e0,0.05e0,
        0.1e0,0.1e0,-4.53e0,6.1e0,0e0,0e0,0.14e0,0.1e0,
        2.47e0,2.7e0,-4.46e0,4.9e0,0.06e0,-0.06e0,0.11e0,0.1e0,
        -5.05e0,5.5e0,0.84e0,0.9e0,0.12e0,0.11e0,0.02e0,-0.02e0,
        4.97e0,-5.4e0,-1.71e0,0e0,-0.12e0,-0.11e0,0e0,0.04e0,
        -0.99e0,-1.3e0,4.22e0,-5.7e0,-0.03e0,0.02e0,-0.13e0,-0.09e0,
        0.99e0,1.4e0,4.22e0,-5.6e0,0.03e0,-0.02e0,-0.13e0,-0.09e0,
        -4.69e0,-5.2e0,0e0,0e0,-0.12e0,0.1e0,-3.42e0,0e0,
        6.09e0,0e0,0e0,0.08e0,0e0,-0.14e0,-4.65e0,-5.1e0,
        0e0,0e0,-0.11e0,0.1e0,0e0,0e0,-4.53e0,-5e0,
        0e0,0e0,-0.11e0,0.1e0,-2.43e0,-2.7e0,-3.82e0,4.2e0,
        -0.06e0,0.05e0,0.1e0,0.09e0,0e0,0e0,-4.53e0,4.9e0,
        0e0,0e0,0.11e0,0.1e0,-4.49e0,-4.9e0,0e0,0e0,
        -0.11e0,0.1e0,2.67e0,-2.9e0,-3.62e0,-3.9e0,-0.06e0,-0.06e0,
        -0.09e0,0.08e0,3.94e0,-5.3e0,0e0,0e0,-0.12e0,-3.38e0,
        3.7e0,-2.78e0,-3.1e0,0.08e0,0.08e0,-0.07e0,0.06e0,3.18e0,
        -3.5e0,-2.82e0,-3.1e0,-0.08e0,-0.07e0,-0.07e0,0.06e0,-5.77e0,
        0e0,1.87e0,0e0,0e0,0.13e0,0e0,-0.04e0,3.54e0,
        -4.8e0,-0.64e0,-0.9e0,-0.11e0,0e0,-0.02e0,-3.5e0,-4.7e0,
        0.68e0,-0.9e0,-0.11e0,0e0,-0.02e0,5.49e0,0e0,0e0,
        0e0,0e0,-0.12e0,1.83e0,-2.5e0,2.63e0,3.5e0,-0.06e0,
        0e0,0.08e0,3.02e0,-4.1e0,0.68e0,0.9e0,-0.09e0,0e0,
        0.02e0,0e0,0e0,5.21e0,0e0,0e0,0e0,0e0,
        -0.12e0,-3.54e0,3.8e0,2.7e0,3.6e0,-1.35e0,1.8e0,0.08e0,
        0e0,0.04e0,-2.9e0,3.9e0,0.68e0,0.9e0,0.09e0,0e0,
        0.02e0,0.8e0,-1.1e0,-2.78e0,-3.7e0,-0.02e0,0e0,-0.08e0,
        4.1e0,0e0,-2.39e0,0e0,0e0,-0.09e0,0e0,0.05e0,
        -1.59e0,2.1e0,2.27e0,3e0,0.05e0,0e0,0.07e0,-2.63e0,
        3.5e0,-0.48e0,-0.6e0,-2.94e0,-3.2e0,-2.94e0,3.2e0,2.27e0,
        -3e0,-1.11e0,-1.5e0,-0.07e0,0e0,-0.03e0,-0.56e0,-0.8e0,
        -2.35e0,3.1e0,0e0,-0.6e0,-3.42e0,1.9e0,-0.12e0,-0.1e0,
        2.63e0,-2.9e0,2.51e0,2.8e0,-0.64e0,0.7e0,-0.48e0,-0.6e0,
        2.19e0,-2.9e0,0.24e0,-0.3e0,2.15e0,2.9e0,2.15e0,-2.9e0,
        0.52e0,0.7e0,2.07e0,-2.8e0,-3.1e0,0e0,1.79e0,0e0,
        0e0,0.07e0,0e0,-0.04e0,0.88e0,0e0,-3.46e0,2.11e0,
        2.8e0,-0.36e0,0.5e0,3.54e0,-0.2e0,-3.5e0,-1.39e0,1.5e0,
        -1.91e0,-2.1e0,-1.47e0,2e0,1.39e0,1.9e0,2.07e0,-2.3e0,
        0.91e0,1e0,1.99e0,-2.7e0,3.3e0,0e0,0.6e0,-0.44e0,
        -0.7e0,-1.95e0,2.6e0,2.15e0,-2.4e0,-0.6e0,-0.7e0,3.3e0,
        0.84e0,0e0,-3.1e0,-3.1e0,0e0,-0.72e0,-0.32e0,0.4e0,
        -1.87e0,-2.5e0,1.87e0,-2.5e0,0.32e0,0.4e0,-0.24e0,0.3e0,
        -1.87e0,-2.5e0,-0.24e0,-0.3e0,1.87e0,-2.5e0,-2.7e0,0e0,
        1.55e0,2.03e0,2.2e0,-2.98e0,-1.99e0,-2.2e0,0.12e0,-0.1e0,
        -0.4e0,0.5e0,1.59e0,2.1e0,0e0,0e0,-1.79e0,2e0,
        -1.03e0,1.4e0,-1.15e0,-1.6e0,0.32e0,0.5e0,1.39e0,-1.9e0,
        2.35e0,-1.27e0,1.7e0,0.6e0,0.8e0,-0.32e0,-0.4e0,1.35e0,
        -1.8e0,0.44e0,0e0,2.23e0,-0.84e0,0.9e0,-1.27e0,-1.4e0,
        -1.47e0,1.6e0,-0.28e0,-0.3e0,-0.28e0,0.4e0,-1.27e0,-1.7e0,
        0.28e0,-0.4e0,-1.43e0,-1.5e0,0e0,0e0,-1.27e0,-1.7e0,
        2.11e0,-0.32e0,-0.4e0,-1.23e0,1.6e0,1.19e0,-1.3e0,-0.72e0,
        -0.8e0,0.72e0,-0.8e0,-1.15e0,-1.3e0,-1.35e0,-1.5e0,-1.19e0,
        -1.6e0,-0.12e0,0.2e0,1.79e0,0e0,-0.88e0,-0.28e0,0.4e0,
        1.11e0,1.5e0,-1.83e0,0e0,0.56e0,-0.12e0,0.1e0,-1.27e0,
        -1.4e0,0e0,0e0,1.15e0,1.5e0,-0.12e0,0.2e0,1.11e0,
        1.5e0,0.36e0,-0.5e0,-1.07e0,-1.4e0,-1.11e0,1.5e0,1.67e0,
        0e0,0.8e0,-1.11e0,0e0,1.43e0,1.23e0,-1.3e0,-0.24e0,
        -1.19e0,-1.3e0,-0.24e0,0.2e0,-0.44e0,-0.9e0,-0.95e0,1.1e0,
        1.07e0,-1.4e0,1.15e0,-1.3e0,1.03e0,-1.1e0,-0.56e0,-0.6e0,
        -0.68e0,0.9e0,-0.76e0,-1e0,-0.24e0,-0.3e0,0.95e0,-1.3e0,
        0.56e0,0.7e0,0.84e0,-1.1e0,-0.56e0,0e0,-1.55e0,0.91e0,
        -1.3e0,0.28e0,0.3e0,0.16e0,-0.2e0,0.95e0,1.3e0,0.4e0,
        -0.5e0,-0.88e0,-1.2e0,0.95e0,-1.1e0,-0.48e0,-0.5e0,0e0,
        0e0,-1.07e0,1.2e0,0.44e0,-0.5e0,0.95e0,1.1e0,0e0,
        0e0,0.92e0,-1.3e0,0.95e0,1e0,-0.52e0,0.6e0,1.59e0,
        0.24e0,-0.4e0,0.91e0,1.2e0,0.84e0,-1.1e0,-0.44e0,-0.6e0,
        0.84e0,1.1e0,-0.44e0,0.6e0,-0.44e0,0.6e0,-0.84e0,-1.1e0,
        -0.8e0,0e0,1.35e0,0.76e0,0.2e0,-0.91e0,-1e0,0.2e0,
        -0.3e0,-0.91e0,-1.2e0,-0.95e0,1e0,-0.48e0,-0.5e0,0.88e0,
        1e0,0.48e0,-0.5e0,-0.95e0,-1.1e0,0.2e0,-0.2e0,-0.99e0,
        1.1e0,-0.84e0,1.1e0,-0.24e0,-0.3e0,0.2e0,-0.3e0,0.84e0,
        1.1e0,-1.39e0,0e0,-0.28e0,-0.16e0,0.2e0,0.84e0,1.1e0,
        0e0,0e0,1.39e0,0e0,0e0,-0.95e0,1e0,1.35e0,
        -0.99e0,0e0,0.88e0,-0.52e0,0e0,-1.19e0,0.2e0,0.2e0,
        0.76e0,-1e0,0e0,0e0,0.76e0,1e0,0e0,0e0,
        0.76e0,1e0,-0.76e0,1e0,0e0,0e0,1.23e0,0.76e0,
        0.8e0,-0.32e0,0.4e0,-0.72e0,0.8e0,-0.4e0,-0.4e0,0e0,
        0e0,-0.8e0,-0.9e0,-0.68e0,0.9e0,-0.16e0,-0.2e0,-0.16e0,
        -0.2e0,0.68e0,-0.9e0,-0.36e0,0.5e0,-0.56e0,-0.8e0,0.72e0,
        -0.9e0,0.44e0,-0.6e0,-0.48e0,-0.7e0,-0.16e0,0e0,-1.11e0,
        0.32e0,0e0,-1.07e0,0.6e0,-0.8e0,-0.28e0,-0.4e0,-0.64e0,
        0e0,0.91e0,1.11e0,0.64e0,-0.9e0,0.76e0,-0.8e0,0e0,
        0e0,-0.76e0,-0.8e0,1.03e0,0e0,-0.36e0,-0.64e0,-0.7e0,
        0.36e0,-0.4e0,1.07e0,0.36e0,-0.5e0,-0.52e0,-0.7e0,0.6e0,
        0e0,0.88e0,0.95e0,0e0,0.48e0,0.16e0,-0.2e0,0.6e0,
        0.8e0,0.16e0,-0.2e0,-0.6e0,-0.8e0,0e0,-1e0,0.12e0,
        0.2e0,0.16e0,-0.2e0,0.68e0,0.7e0,0.59e0,-0.8e0,-0.99e0,
        -0.56e0,-0.6e0,0.36e0,-0.4e0,-0.68e0,-0.7e0,-0.68e0,-0.7e0,
        -0.36e0,-0.5e0,-0.44e0,0.6e0,0.64e0,0.7e0,-0.12e0,0.1e0,
        -0.52e0,0.6e0,0.36e0,0.4e0,0e0,0e0,0.95e0,-0.84e0,
        0e0,0.44e0,0.56e0,0.6e0,0.32e0,-0.3e0,0e0,0e0,
        0.6e0,0.7e0,0e0,0e0,0.6e0,0.7e0,-0.12e0,-0.2e0,
        0.52e0,-0.7e0,0e0,0e0,0.56e0,0.7e0,-0.12e0,0.1e0,
        -0.52e0,-0.7e0,0e0,0e0,0.88e0,-0.76e0,0e0,-0.44e0,
        0e0,0e0,-0.52e0,-0.7e0,0.52e0,-0.7e0,0.36e0,-0.4e0,
        -0.44e0,-0.5e0,0e0,0e0,0.6e0,0.6e0,0.84e0,0e0,
        0.12e0,-0.24e0,0e0,0.8e0,-0.56e0,0.6e0,-0.32e0,-0.3e0,
        0.48e0,-0.5e0,0.28e0,-0.3e0,-0.48e0,-0.5e0,0.12e0,0.2e0,
        0.48e0,-0.6e0,0.48e0,0.6e0,-0.12e0,0.2e0,0.24e0,0e0,
        0.76e0,-0.52e0,-0.6e0,-0.52e0,0.6e0,0.48e0,-0.5e0,-0.24e0,
        -0.3e0,0.12e0,-0.1e0,0.48e0,0.6e0,0.52e0,-0.2e0,0.36e0,
        0.4e0,-0.44e0,0.5e0,-0.24e0,-0.3e0,-0.48e0,-0.6e0,-0.44e0,
        -0.6e0,-0.12e0,0.1e0,0.76e0,0.76e0,0.2e0,-0.2e0,0.48e0,
        0.5e0,0.4e0,-0.5e0,-0.24e0,-0.3e0,0.44e0,-0.6e0,0.44e0,
        -0.6e0,0.36e0,0e0,-0.64e0,0.72e0,0e0,-0.12e0,0e0,
        -0.1e0,-0.4e0,-0.6e0,-0.2e0,-0.2e0,-0.44e0,0.5e0,-0.44e0,
        0.5e0,0.2e0,0.2e0,-0.44e0,-0.5e0,0.2e0,-0.2e0,-0.2e0,
        0.2e0,-0.44e0,-0.5e0,0.64e0,0e0,0.32e0,-0.36e0,0.5e0,
        -0.2e0,-0.3e0,0.12e0,-0.1e0,0.48e0,0.5e0,-0.12e0,0.3e0,
        -0.36e0,-0.5e0,0e0,0e0,0.48e0,0.5e0,-0.48e0,0.5e0,
        0.68e0,0e0,-0.12e0,0.56e0,-0.4e0,0.44e0,-0.5e0,-0.12e0,
        -0.1e0,0.24e0,0.3e0,-0.4e0,0.4e0,0.64e0,0e0,-0.24e0,
        0.64e0,0e0,-0.2e0,0e0,0e0,0.44e0,-0.5e0,0.44e0,
        0.5e0,-0.12e0,0.2e0,-0.36e0,-0.5e0,0.12e0,0e0,0.64e0,
        -0.4e0,0.5e0,0e0,0.1e0,0e0,0e0,-0.4e0,0.5e0,
        0e0,0e0,-0.4e0,-0.5e0,0.56e0,0e0,0.28e0,0e0,
        0.1e0,0.36e0,0.5e0,0e0,-0.1e0,0.36e0,-0.5e0,0.36e0,
        0.5e0,0e0,-0.1e0,0.24e0,-0.2e0,-0.36e0,-0.4e0,0.16e0,
        0.2e0,0.4e0,-0.4e0,0e0,0e0,-0.36e0,-0.5e0,-0.36e0,
        -0.5e0,-0.32e0,-0.5e0,-0.12e0,0.1e0,0.2e0,0.2e0,-0.36e0,
        0.4e0,-0.6e0,0.6e0,0.28e0,0e0,0.52e0,0.12e0,-0.1e0,
        0.4e0,0.4e0,0e0,-0.5e0,0.2e0,-0.2e0,-0.32e0,0.4e0,
        0.16e0,0.2e0,-0.16e0,0.2e0,0.32e0,0.4e0,0.56e0,0e0,
        -0.12e0,0.32e0,-0.4e0,-0.16e0,-0.2e0,0e0,0e0,0.4e0,
        0.4e0,-0.4e0,-0.4e0,-0.4e0,0.4e0,-0.36e0,0.4e0,0.12e0,
        0.1e0,0e0,0.1e0,0.36e0,0.4e0,0e0,-0.1e0,0.36e0,
        0.4e0,-0.36e0,0.4e0,0e0,0.1e0,0.32e0,0e0,0.44e0,
        0.12e0,0.2e0,0.28e0,-0.4e0,0e0,0e0,0.36e0,0.4e0,
        0.32e0,-0.4e0,-0.16e0,0.12e0,0.1e0,0.32e0,-0.4e0,0.2e0,
        0.3e0,-0.24e0,0.3e0,0e0,0.1e0,0.32e0,0.4e0,0e0,
        -0.1e0,-0.32e0,-0.4e0,-0.32e0,0.4e0,0e0,0.1e0,-0.52e0,
        -0.52e0,0.52e0,0.32e0,-0.4e0,0e0,0e0,0.32e0,0.4e0,
        0.32e0,-0.4e0,0e0,0e0,-0.32e0,-0.4e0,-0.32e0,0.4e0,
        0.32e0,0.4e0,0e0,0e0,0.32e0,0.4e0,0e0,0e0,
        -0.32e0,-0.4e0,0e0,0e0,0.32e0,0.4e0,0.16e0,0.2e0,
        0.32e0,-0.3e0,-0.16e0,0e0,-0.48e0,-0.2e0,0.2e0,-0.28e0,
        -0.3e0,0.28e0,-0.4e0,0e0,0e0,0.28e0,-0.4e0,0e0,
        0e0,0.28e0,-0.4e0,0e0,0e0,-0.28e0,-0.4e0,0.28e0,
        0.4e0,-0.28e0,-0.4e0,-0.48e0,-0.2e0,0.2e0,0.24e0,0.3e0,
        0.44e0,0e0,0.16e0,0.24e0,0.3e0,0.16e0,-0.2e0,0.24e0,
        0.3e0,-0.12e0,0.2e0,0.2e0,0.3e0,-0.16e0,0.2e0,0e0,
        0e0,0.44e0,-0.32e0,0.3e0,0.24e0,0e0,-0.36e0,0.36e0,
        0e0,0.24e0,0.12e0,-0.2e0,0.2e0,0.3e0,-0.12e0,0e0,
        -0.28e0,0.3e0,-0.24e0,0.3e0,0.12e0,0.1e0,-0.28e0,-0.3e0,
        -0.28e0,0.3e0,0e0,0e0,-0.28e0,-0.3e0,0e0,0e0,
        -0.28e0,-0.3e0,0e0,0e0,0.28e0,0.3e0,0e0,0e0,
        -0.28e0,-0.3e0,-0.28e0,0.3e0,0e0,0e0,-0.28e0,-0.3e0,
        0e0,0e0,0.28e0,0.3e0,0e0,0e0,-0.28e0,0.3e0,
        0.28e0,-0.3e0,-0.28e0,0.3e0,0.4e0,0.4e0,-0.24e0,0.3e0,
        0e0,-0.1e0,0.16e0,0e0,0.36e0,-0.2e0,0.3e0,-0.12e0,
        -0.1e0,-0.24e0,-0.3e0,0e0,0e0,-0.24e0,0.3e0,-0.24e0,
        0.3e0,0e0,0e0,-0.24e0,0.3e0,-0.24e0,0.3e0,0.24e0,
        -0.3e0,0e0,0e0,0.24e0,-0.3e0,0e0,0e0,0.24e0,
        0.3e0,0.24e0,-0.3e0,0.24e0,0.3e0,-0.24e0,0.3e0,-0.24e0,
        0.3e0,-0.2e0,0.2e0,-0.16e0,-0.2e0,0e0,0e0,-0.32e0,
        0.2e0,0e0,0.1e0,0.2e0,-0.3e0,0.2e0,-0.2e0,0.12e0,
        0.2e0,-0.16e0,0.2e0,0.16e0,0.2e0,0.2e0,0.3e0,0.2e0,
        0.3e0,0e0,0e0,-0.2e0,0.3e0,0e0,0e0,0.2e0,
        0.3e0,-0.2e0,-0.3e0,-0.2e0,-0.3e0,0.2e0,-0.3e0,0e0,
        0e0,0.2e0,0.3e0,0e0,0e0,0.2e0,0.3e0,0e0,
        0e0,0.2e0,0.3e0,0e0,0e0,0.2e0,0.3e0,0e0,
        0e0,0.2e0,-0.3e0,0e0,0e0,-0.2e0,-0.3e0,0e0,
        0e0,-0.2e0,0.3e0,0e0,0e0,-0.2e0,0.3e0,0e0,
        0e0,0.36e0,0e0,0e0,0.36e0,0.12e0,0.1e0,-0.24e0,
        0.2e0,0.12e0,-0.2e0,-0.16e0,-0.2e0,-0.13e0,0.1e0,0.22e0,
        0.21e0,0.2e0,0e0,-0.28e0,0.32e0,0e0,-0.12e0,-0.2e0,
        -0.2e0,0.12e0,-0.1e0,0.12e0,0.1e0,-0.2e0,0.2e0,0e0,
        0e0,-0.32e0,0.32e0,0e0,0e0,0.32e0,0.32e0,0e0,
        0e0,-0.24e0,-0.2e0,0.24e0,0.2e0,0.2e0,0e0,-0.24e0,
        0e0,0e0,-0.24e0,-0.2e0,0e0,0e0,0.24e0,0.2e0,
        -0.24e0,-0.2e0,0e0,0e0,-0.24e0,0.2e0,0.16e0,-0.2e0,
        0.12e0,0.1e0,0.2e0,0.2e0,0e0,-0.1e0,-0.12e0,0.1e0,
        -0.16e0,-0.2e0,-0.12e0,-0.1e0,-0.16e0,0.2e0,0.2e0,0.2e0,
        0e0,0e0,-0.2e0,0.2e0,-0.2e0,0.2e0,-0.2e0,0.2e0,
        -0.2e0,0.2e0,0.2e0,-0.2e0,-0.2e0,-0.2e0,0e0,0e0,
        -0.2e0,0.2e0,0.2e0,0e0,-0.2e0,0e0,0e0,-0.2e0,
        0.2e0,-0.2e0,0.2e0,-0.2e0,-0.2e0,-0.2e0,-0.2e0,0e0,
        0e0,0.2e0,0.2e0,0.2e0,0.2e0,0.12e0,-0.2e0,-0.12e0,
        -0.1e0,0.28e0,-0.28e0,0.16e0,-0.2e0,0e0,-0.1e0,0e0,
        0.1e0,-0.16e0,0.2e0,0e0,-0.1e0,-0.16e0,-0.2e0,0e0,
        -0.1e0,0.16e0,-0.2e0,0.16e0,-0.2e0,0e0,0e0,0.16e0,
        0.2e0,-0.16e0,0.2e0,0e0,0e0,0.16e0,0.2e0,0.16e0,
        -0.2e0,0.16e0,-0.2e0,-0.16e0,0.2e0,0.16e0,-0.2e0,0e0,
        0e0,0.16e0,0.2e0,0e0,0e0,0.16e0,0.2e0,0e0,
        0e0,-0.16e0,-0.2e0,0.16e0,-0.2e0,-0.16e0,-0.2e0,0e0,
        0e0,-0.16e0,-0.2e0,0e0,0e0,-0.16e0,0.2e0,0e0,
        0e0,0.16e0,-0.2e0,0.16e0,0.2e0,0.16e0,0.2e0,0e0,
        0e0,-0.16e0,-0.2e0,0e0,0e0,-0.16e0,-0.2e0,0e0,
        0e0,0.16e0,0.2e0,0.16e0,0.2e0,0e0,0e0,0.16e0,
        0.2e0,0.16e0,-0.2e0,0.16e0,0.2e0,0e0,0e0,-0.16e0,
        0.2e0,0e0,0.1e0,0.12e0,-0.2e0,0.12e0,-0.2e0,0e0,
        -0.1e0,0e0,-0.1e0,0.12e0,0.2e0,0e0,-0.1e0,-0.12e0,
        0.2e0,-0.15e0,0.2e0,-0.24e0,0.24e0,0e0,0e0,0.24e0,
        0.24e0,0.12e0,-0.2e0,-0.12e0,-0.2e0,0e0,0e0,0.12e0,
        0.2e0,0.12e0,-0.2e0,0.12e0,0.2e0,0.12e0,0.2e0,0.12e0,
        0.2e0,0.12e0,-0.2e0,-0.12e0,0.2e0,0e0,0e0,0.12e0,
        0.2e0,0.12e0,0e0,-0.2e0,0e0,0e0,-0.12e0,-0.2e0,
        0.12e0,-0.2e0,0e0,0e0,0.12e0,0.2e0,-0.12e0,0.2e0,
        -0.12e0,0.2e0,0.12e0,-0.2e0,0e0,0e0,0.12e0,0.2e0,
        0.2e0,0e0,0.12e0,0e0,0e0,-0.12e0,0.2e0,0e0,
        0e0,-0.12e0,-0.2e0,0e0,0e0,-0.12e0,-0.2e0,-0.12e0,
        -0.2e0,0e0,0e0,0.12e0,-0.2e0,0.12e0,-0.2e0,0.12e0,
        0.2e0,-0.12e0,-0.2e0,0e0,0e0,0.12e0,-0.2e0,0.12e0,
        -0.2e0,0.12e0,0.2e0,0.12e0,0e0,0.2e0,-0.12e0,-0.2e0,
        0e0,0e0,0.12e0,0.2e0,-0.16e0,0e0,0.16e0,-0.2e0,
        0.2e0,0e0,0e0,-0.2e0,0e0,0e0,-0.2e0,0.2e0,
        0e0,0e0,0.2e0,0.2e0,-0.2e0,0e0,0e0,-0.2e0,
        0.12e0,0e0,-0.16e0,0.2e0,0e0,0e0,0.2e0,0.12e0,
        -0.1e0,0e0,0.1e0,0.16e0,-0.16e0,-0.16e0,-0.16e0,-0.16e0,
        -0.16e0,0e0,0e0,-0.16e0,0e0,0e0,-0.16e0,-0.16e0,
        -0.16e0,0e0,0e0,-0.16e0,0e0,0e0,0.16e0,0e0,
        0e0,0.16e0,0e0,0e0,0.16e0,0.16e0,0e0,0e0,
        -0.16e0,0e0,0e0,-0.16e0,-0.16e0,0e0,0e0,0.16e0,
        0e0,0e0,-0.16e0,-0.16e0,0e0,0e0,-0.16e0,-0.16e0,
        0.12e0,0.1e0,0.12e0,-0.1e0,0.12e0,0.1e0,0e0,0e0,
        0.12e0,0.1e0,-0.12e0,0.1e0,0e0,0e0,0.12e0,0.1e0,
        0.12e0,-0.1e0,0e0,0e0,-0.12e0,-0.1e0,0e0,0e0,
        0.12e0,0.1e0,0.12e0,0e0,0e0,0.12e0,0e0,0e0,
        -0.12e0,0e0,0e0,0.12e0,0.12e0,0.12e0,0.12e0,0.12e0,
        0e0,0e0,0.12e0,0e0,0e0,0.12e0,0.12e0,0e0,
        0e0,0.12e0,0e0,0e0,0.12e0,-0.12e0,-0.12e0,0.12e0,
        0.12e0,-0.12e0,-0.12e0,0e0,0e0,0.12e0,-0.12e0,0.12e0,
        0.12e0,-0.12e0,-0.12e0,0e0,0e0,-0.12e0,-0.12e0,0e0,
        0e0,-0.12e0,0.12e0,0e0,0e0,0.12e0,0e0,0e0,
        0.12e0,0e0,0e0,0.12e0,-0.12e0,0e0,0e0,-0.12e0,
        0.12e0,-0.12e0,-0.12e0,0.12e0,0e0,0e0,0.12e0,0.12e0,
        0.12e0,-0.12e0,0e0,0e0,-0.12e0,-0.12e0,-0.12e0,0e0,
        0e0,-0.12e0,-0.12e0,0e0,0e0,0.12e0,0.12e0,0e0,
        0e0,-0.12e0,-0.12e0,-0.12e0,-0.12e0,0.12e0,0e0,0e0,
        0.12e0,-0.12e0,0e0,0e0,-0.12e0,-0.12e0,0e0,0e0,
        0.12e0,-0.12e0,-0.12e0,-0.12e0,-0.12e0,0.12e0,0.12e0,-0.12e0,
        -0.12e0,0e0,0e0,-0.12e0,0e0,0e0,-0.12e0,0.12e0,
        0e0,0e0,0.12e0,0e0,0e0,-0.12e0,-0.12e0,0e0,
        0e0,-0.12e0,-0.12e0,0.12e0,0e0,0e0,0.12e0,0.12e0,
        0e0,0e0,0.12e0,0e0,0e0,0.12e0,0.12e0,0.08e0,
        0e0,0.04e0]
    
    #振幅使用:X或Y, sin或cos, T的幂。
    JAXY=[0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1]
    JASC=[0,1,1,0,1,0,0,1,0,1,1,0,1,0,0,1,0,1,1,0]
    JAPT=[0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4]
    
    #给定日期与参考历元之间的时间间隔，儒略世纪数.
    T=((DATE1-DJ00)+DATE2)/DJC
    
    #T的幂.
    W=1.0
    for JPT in range(MAXPT+1):
        PT[JPT]=W
        W=W*T
    
    #日月基本(Delaunay)参数(IERS 2003)
    
    #月亮的平近点角.
    FA[0]=pymFal03(T)

    #太阳的平近点角.
    FA[1]=pymFalp03(T)

    #月球纬度的平均辐角.
    FA[2]=pymFaf03(T)

    #月亮到太阳的平均距角
    FA[3]=pymFad03(T)

    #月亮升交点的平黄经
    FA[4]=pymFaom03(T)

    #行星的平黄经，水星到海王星
    FA[5]=pymFame03(T)
    FA[6]=pymFave03(T)
    FA[7]=pymFae03(T)
    FA[8]=pymFama03(T)
    FA[9]=pymFaju03(T)
    FA[10]=pymFasa03(T)
    FA[11]=pymFaur03(T)
    FA[12]=pymFane03(T)
    
    #黄经上的累计岁差
    FA[13]=pymFapa03(T)
    
    #岁差章动的多项式部分
    for JXY in range(2):
        XYPR[JXY]=0.0
        for J in range(MAXPT,-1,-1):
            XYPR[JXY]=XYPR[JXY]+XYP[JXY][J]*PT[J]
            
    #章动的周期项，行星
    
    #初始化X,Y.
    for JXY in range(2):
        XYPL[JXY]=0.0
        
    #每个频率的系数从后往前处理
    IALAST=NA
    for IFREQ in range(NFPL-1,-1,-1):
        
        #获得系数方程
        ARG=0.0
        for I in range(14):
            M=MFAPL[IFREQ][I]
            if (M!=0):
                ARG=ARG+float(M)*FA[I]
        SC[0]=ma.sin(ARG)
        SC[1]=ma.cos(ARG)
        
        #对这个频率的幅度从后往前处理
        IA=NC[IFREQ+NFLS]
        for I in range(IALAST,IA-1,-1):
            
            #系数数目 (0 = 1st).
            J=I-IA
            
            #X/Y.
            JXY=JAXY[J]
    
            #Sin/cos.
            JSC=JASC[J]

            # T的幂次.
            JPT=JAPT[J]
    
            #总计该分量
            XYPL[JXY]=XYPL[JXY]+A[I-1]*SC[JSC]*PT[JPT]
        IALAST=IA-1
            
    #章动的周期项，日月
    
    #初始化X,Y.
    for JXY in range(2):
        XYLS[JXY]=0.0
    
    #每个频率的系数从后往前处理
    for IFREQ in range(NFLS-1,-1,-1):
        
        #获得系数方程
        ARG=0.0
        for I in range(5):
            M=MFALS[IFREQ][I]
            if (M!=0):
                ARG=ARG+float(M)*FA[I]
        SC[0]=ma.sin(ARG)
        SC[1]=ma.cos(ARG)
        
        #对这个频率的幅度从后往前处理
        IA=NC[IFREQ]
        for I in range(IALAST,IA-1,-1):
           
            #系数数目 (0 = 1st).
            J=I-IA
            #X/ Y.
            JXY=JAXY[J]
    
            #Sin/cos.
            JSC=JASC[J]

            #T的幂次.
            JPT=JAPT[J]
    
            #总计该分量
            XYLS[JXY]=XYLS[JXY]+A[I-1]*SC[JSC]*PT[JPT]
        IALAST=IA-1
                
    #CIP的单位向量分量
    X=DAS2R*(XYPR[0]+(XYLS[0]+XYPL[0])/1e6)
    Y=DAS2R*(XYPR[1]+(XYLS[1]+XYPL[1])/1e6)

    return(X,Y)


def pymPmat00(DATE1,DATE2):
    '''
    Precession matrix (including frame bias) from GCRS to a specified
    date, IAU 2000 model.

    Parameters
    ----------
    date1 : flaot
        TT as a 2-part Julian Date     
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    rbp : list(3,3)
        bias-precession matrix

    '''
    #调用pymBp00获得需要的矩阵，舍去无关项
    RB,RP,RBP=pymBp00(DATE1,DATE2)
    
    return(RBP)


def pymPrec76 (DATE01,DATE02,DATE11,DATE12):
    '''
    IAU 1976 precession model.
    
    This function forms the three Euler angles which implement general
    precession between two dates, using the IAU 1976 model (as for the
    FK5 catalog).

    Parameters
    ----------
    date01 : float
        TDB starting date    
    date02 : float
        TDB starting date    
    date11 : float
        TDB ending date    
    date12 : flaot
        TDB ending date    

    Returns
    -------
    zeta : float
        1st rotation: radians cw around z    
    z : float
        3rd rotation: radians cw around z    
    theta : float
        2nd rotation: radians ccw around y

    '''
    #1角秒对应的弧度
    DAS2R=4.848136811095359935899141e-6

    #参考历元(J2000.0), JD
    DJ00=2451545.0

    #1儒略世纪对应的天数
    DJC=36525.0
    
    #给定的开始时间与参考历元之间的时间间隔儒略世纪数.
    T0=((DATE01-DJ00)+DATE02)/DJC
    
    #计算的岁差的持续时间（结束时间-开始时间），儒略世纪数.
    T=((DATE11-DATE01)+(DATE12-DATE02))/DJC
    
    #欧拉角
    TAS2R=T*DAS2R
    W=2306.2181+(1.39656-0.000139*T0)*T0

    ZETA=(W+((0.30188-0.000344*T0)+0.017998*T)*T)*TAS2R

    Z=(W+((1.09468+0.000066*T0)+0.018203*T)*T)*TAS2R

    THETA=((2004.3109+(-0.85330-0.000217*T0)*T0)+\
           ((-0.42665-0.000217*T0)-0.041833*T)*T)*TAS2R 
    
    return(ZETA,Z,THETA)


def pymPmat76(DATE1,DATE2):
    '''
    Precession matrix from J2000.0 to a specified date, IAU 1976 model.

    Parameters
    ----------
    date1 : float
        ending date, TT    
    date2 : float
        ending date, TT

    Returns
    -------
    rmatp : list(3,3)
        precession matrix, J2000.0 -> date1+date2

    ''' 
    #参考历元(J2000.0), JD
    DJ00=2451545.0
    
    #岁差欧拉角, J2000.0到给定日期.
    ZETA,Z,THETA=pymPrec76(DJ00,0.0,DATE1,DATE2)
    
    #构建旋转矩阵
    WMAT=pymIr()
    WMAT=pymRz(-ZETA,WMAT)
    WMAT=pymRy(THETA,WMAT)
    WMAT=pymRz(-Z,WMAT)
    RMATP=pymCr(WMAT)
    
    return(RMATP)


def pymSepp (A,B):
    '''
    Angular separation between two p-vectors.

    Parameters
    ----------
    a : list(3)
        first p-vector (not necessarily unit length)    
    b : list(3)
        second p-vector (not necessarily unit length)

    Returns
    -------
    function value : float
        angular separation (radians, always positive)

    '''
    #向量夹角的正弦，乘以两个模
    AXB=pymPxp(A,B)
    SS=pymPm(AXB)
    
    #角的余弦，乘以两个模。
    CS=pymPdp(A,B)
    
    #角度
    if (SS!=0.0)|(CS!=0.0):
        S=ma.atan2(SS,CS)
    else:
        S=0.0
    
    return(S)


def pymSeps(AL,AP,BL,BP):
    '''
    Angular separation between two sets of spherical coordinates.

    Parameters
    ----------
    al : float
        first longitude (radians)    
    ap : float
        first latitude (radians)    
    bl : float
        second longitude (radians)    
    bp : float
        second latitude (radians)

    Returns
    -------
    function value : float
        angular separation (radians)

    '''
    #球面到笛卡尔
    AC=pymS2c(AL,AP)
    BC=pymS2c(BL,BP)
    
    #两个向量之间的角度差
    S=pymSepp(AC,BC)
    
    return(S)


def pymStarpm (RA1,DEC1,PMR1,PMD1,PX1,RV1,EP1A,EP1B,EP2A,EP2B):
    '''
    Star proper motion:  update star catalog data for space motion.

    Parameters
    ----------
    ra1 : float
        right ascension (radians), before    
    dec1 : float
        declination (radians), before    
    pmr1 : float
        RA proper motion (radians/year), before    
    pmd1 : float
        Dec proper motion (radians/year), before    
    px1 : float
        Dec proper motion (radians/year), before    
    rv1 : float
        radial velocity (km/s, +ve = receding), before    
    ep1a : float
        "before" epoch, part A    
    ep1b : float
        "before" epoch, part B    
    ep2a : float
        "after" epoch, part A    
    ep2b : float
        "after" epoch, part B

    Returns
    -------
    ra2 : float
        right ascension (radians), after    
    dec2 : float
        declination (radians), after    
    pmr2 : float
        RA proper motion (radians/year), after    
    pmd2 : float
        Dec proper motion (radians/year), after    
    px2 : float
        parallax (arcseconds), after    
    rv2 : flaot
        radial velocity (km/s, +ve = receding), after    
    J : ValueError
        4: 'solution not converge',
        2: 'excessive velocity',
        1: 'distance overridden',
        0: 'no warnings or errors',
       -1: 'system error',
    '''
    #1天的秒数
    D2S=86400.0
    
    #光速 (m/s)
    CMPS=299792458.0
    
    #天文单位(m, IAU 2012)
    AUM=149597870.7e3
    
    #光速(au/ day)
    C=D2S*CMPS/AUM
  
    #之前的坐标、自行转换为位置-速度矩阵.
    PV1,J1=pymStarpv(RA1,DEC1,PMR1,PMD1,PX1,RV1)
    
    #观测到时经历的时间(days).
    R=pymPm(PV1[0])
    TL1=R/C
    
    #之前和之后之间的时间间隔(days).
    DT=(EP2A-EP1A)+(EP2B-EP1B)
    
    #将星沿轨道从之前观测位置移动到之后的几何位置。
    PV=pymPvu(DT+TL1, PV1)
    
    I=1
    while I<2:
    
        #从之后的位置，减少观测到时的经历的时间 (days)
        R2=pymPdp(PV[0],PV[0])
        RDV=pymPdp(PV[0],PV[1])
        V2=pymPdp(PV[1],PV[1])    
        C2MV2=C*C-V2
        if (C2MV2<=0.0):
            J=-1
            print('ERROR',J)
            break
        TL2=(-RDV+ma.sqrt(RDV*RDV+C2MV2*R2))/C2MV2
        
        #沿着轨迹的方向更新恒星之后的位置速度向量
        PV2=pymPvu(DT+(TL1-TL2),PV1)
        
        #之后的空间运动的位置速度向量到赤经、赤纬等
        RA2,DEC2,PMR2,PMD2,PX2,RV2,J2=pymPvstar(PV2)
        
        if (J2!=0):
            J1=-1
        J=J1
    
        I+=1
    
    return(RA2,DEC2,PMR2,PMD2,PX2,RV2,J)


def pymPmsafe(RA1,DEC1,PMR1,PMD1,PX1,RV1,EP1A,EP1B,EP2A,EP2B):
    '''
    Star proper motion:  update star catalog data for space motion, with
    special handling to handle the zero parallax case.

    Parameters
    ----------
    ra1 : float
        right ascension (radians), before    
    dec1 : float
        declination (radians), before    
    pmr1 : float
        RA proper motion (radians/year), before     
    pmd1 : float
        Dec proper motion (radians/year), before     
    px1 : float
        parallax (arcseconds), before    
    rv1 : float
        radial velocity (km/s, +ve = receding), before    
    ep1a : float
        "before" epoch, part A    
    ep1b : float
        "before" epoch, part B    
    ep2a : float
        "after" epoch, part A     
    ep2b : float
        "after" epoch, part B    

    Raises
    ------
    ValueError
        4: 'solution not converge',
        2: 'excessive velocity',
        1: 'distance overridden',
        0: 'no warnings or errors',
       -1: 'system error',
        
    Returns
    -------
    ra2 : float
        right ascension (radians), after    
    dec2 : float
        declination (radians), after    
    pmr2 : float
        RA proper motion (radians/year), after      
    pmd2 : float
        Dec proper motion (radians/year), after       
    px2 : float
        parallax (arcseconds), after    
    rv2 : float
        radial velocity (km/s, +ve = receding), after    
    J : ValueError
        4: 'solution not converge',
        2: 'excessive velocity',
        1: 'distance overridden',
        0: 'no warnings or errors',
       -1: 'system error',
    '''
    #允许的最小视差(arcsec)
    PXMIN=5e-7
    
    #给出最大的允许横向运动的因子，大约为1%的光速
    F=326.0
    
    #1年的自行(radians).
    PM=pymSeps(RA1,DEC1,RA1+PMR1,DEC1+PMD1)
    
    #对视差大小进行判断，避免出错
    JPX=0
    PX1A=PX1
    PM=PM*F
    if (PX1A<PM):
        JPX=1
        PX1A=PM
    if (PX1A<PXMIN):
        JPX=1
        PX1A=PXMIN
    
    #使用改进的视差进行计算
    RA2,DEC2,PMR2,PMD2,PX2,RV2,J=pymStarpm(RA1,DEC1,PMR1,PMD1,
                                            PX1A,RV1,EP1A,EP1B,EP2A,EP2B)
   
    C=np.abs(J)%2
    if J<0:
        B=-C
    else:
        B=C
    if (B==0):
        J=J+JPX
        
    return(RA2,DEC2,PMR2,PMD2,PX2,RV2,J)


def pymPn06 (DATE1,DATE2,DPSI,DEPS):
    '''
    Precession-nutation, IAU 2006 model:  a multi-purpose function,
    supporting classical (equinox-based) use directly and CIO-based use
    indirectly.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date       
    dpsi : float
        nutation    
    deps : float
        nutation    

    Returns
    -------
    epsa : flaot
        mean obliquity    
    rb : list(3,3)
        frame bias matrix    
    rp : list(3,3)
        precession matrix    
    rbp : list(3,3)
        bias-precession matrix    
    rn : list(3,3)
        nutation matrix    
    rbpn : list(3,3)
        GCRS-to-true matrix

    '''
    #约简的儒略日零点
    DJM0=2400000.5
    
    #参考历元(J2000.0), MJD
    DJM00=51544.5
    
    #参考架偏差-岁差角Fukushima-Williams J2000.0.
    GAMB,PHIB,PSIB,EPS=pymPfw06(DJM0,DJM00)
    
    #参考架偏差矩阵.
    RB=pymFw2m(GAMB,PHIB,PSIB,EPS)
    
    #参考架偏差-岁差 Fukushima-Williams 角.
    GAMB,PHIB,PSIB,EPS=pymPfw06(DATE1,DATE2)
    
    #参考架偏差-岁差矩阵
    RBP=pymFw2m(GAMB,PHIB,PSIB,EPS)
    
    #解出岁差矩阵
    RT=pymTr(RB)
    RP=pymRxr(RBP,RT)
    
    #基于春分点的参考架偏差-岁差-章动矩阵
    RBPN=pymFw2m(GAMB,PHIB,PSIB+DPSI,EPS+DEPS)
    
    #解出章动矩阵
    RT=pymTr(RBP)
    RN=pymRxr(RBPN,RT)
    
    #倾角，给定日期的平均值
    EPSA=EPS
    
    return(EPSA,RB,RP,RBP,RN,RBPN)

    
def pymPn06a(DATE1,DATE2):
    '''
    Precession-nutation, IAU 2006/2000A models:  a multi-purpose function,
    supporting classical (equinox-based) use directly and CIO-based use
    indirectly.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    dpsi : float
        nutation    
    deps : flaot
        nutation     
    epsa : flaot
        mean obliquity    
    rb : list(3,3)
        frame bias matrix    
    rp : list(3,3)
        precession matrix     
    rbp : list(3,3)
        bias-precession matrix    
    rn : list(3,3)
        nutation matrix    
    rbpn : list(3,3)
        GCRS-to-true matrix

    '''    
    #章动
    DPSI,DEPS=pymNut06a(DATE1,DATE2)
    
    #调用pymPn06获得结果
    EPSA,RB,RP,RBP,RN,RBPN=pymPn06(DATE1,DATE2,DPSI,DEPS)
    
    return(DPSI,DEPS,EPSA,RB,RP,RBP,RN,RBPN)


def pymPnm80(DATE1,DATE2):
    '''
    Form the matrix of precession/nutation for a given date, IAU 1976
    precession model, IAU 1980 nutation model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : flaot
        TT as a 2-part Julian Date

    Returns
    -------
    rmatpn : list(3,3)
        combined precession/nutation matrix

    '''
    #岁差矩阵 J2000.0到给定日期.
    RMATP=pymPmat76(DATE1,DATE2)
    
    #章动矩阵
    RMATN=pymNutm80(DATE1,DATE2)
    
    #两个矩阵结合  PN = N x P.
    RMATPN=pymRxr(RMATN,RMATP)
    
    return(RMATPN)


def pymPv2p(PV):
    '''
    Discard velocity component of a pv-vector.

    Parameters
    ----------
    pv : list(2,3)
        pv-vector

    Returns
    -------
    p : list(3)
        p-vector

    '''
    P=pymCp(PV[0])
    
    return(P)


def pymPvdpv (A,B):
    '''
    Inner (=scalar=dot) product of two pv-vectors.

    Parameters
    ----------
    a : list(2,3)
        first pv-vector   
    b : list(2,3)
        second pv-vector   

    Returns
    -------
    adb : list(2)
        A . B 

    '''
    ADB=[0,0]
    #位置向量内积
    ADB[0]=pymPdp(A[0],B[0])
    
    #A位置向量·B速度向量
    ADBD=pymPdp(A[0],B[1])
    
    #A速度向量·B位置向量
    ADDB=pymPdp(A[1],B[0])
    
    #速度项
    ADB[1]=ADBD+ADDB
    
    return(ADB)


def pymPvm(PV):
    '''
    Modulus of pv-vector.

    Parameters
    ----------
    pv : list(2,3)
        pv-vector

    Returns
    -------
    r : float
        modulus of position component
    s : float
        modulus of velocity component

    '''
    #距离
    R=pymPm(PV[0])
    
    #速度
    S=pymPm(PV[1])
    
    return(R,S)


def pymPvup(DT,PV):
    '''
    Update a pv-vector, discarding the velocity component.

    Parameters
    ----------
    dt : float
        time interval    
    pv : list(2,3) 
        pv-vector

    Returns
    -------
    p : list(3) 
        p-vector

    '''
    P=[0,0,0]
    
    for i in range(3):
        P[i]=PV[0][i]+PV[1][i]*DT
    
    return(P)
    

def pymPvxpv(A,B):
    '''
    Outer (=vector=cross) product of two pv-vectors.

    Parameters
    ----------
    a : list(2,3) 
        first pv-vector    
    b : list(2,3) 
        second pv-vector

    Returns
    -------
    axb : list(2,3) 
        a x b

    '''
    AXB=[0,0]
    #输入向量的副本
    WA=pymCpv(A)
    WB=pymCpv(B)
    
    #A x B的位置结果
    AXB[0]=pymPxp(WA[0],WB[0])

    #A x Bdot + Adot x B  的速度结果
    AXBD=pymPxp(WA[0],WB[1])
    ADXB=pymPxp(WA[1],WB[0])
    AXB[1]=pymPpp(AXBD,ADXB)    
    
    return(AXB)


def pymS00a(DATE1,DATE2):
    '''
    The CIO locator s, positioning the Celestial Intermediate Origin on
    the equator of the Celestial Intermediate Pole, using the IAU 2000A
    precession-nutation model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    function value : float
        the CIO locator s in radians

    '''
    #参考架偏差-岁差-章动矩阵 IAU 2000A.
    RBPN=pymPnm00a(DATE1,DATE2)
    
    #分离出CIP的X,Y坐标
    X,Y=pymBpn2xy(RBPN)
    
    #基于X,Y坐标，计算CIO定位角s
    S00=pymS00(DATE1,DATE2,X,Y)
    
    return(S00)


def pymS00b(DATE1,DATE2):
    '''
    The CIO locator s, positioning the Celestial Intermediate Origin on
    the equator of the Celestial Intermediate Pole, using the IAU 2000B
    precession-nutation model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    function value : float
        the CIO locator s in radians

    '''
    #参考架偏差-岁差-章动矩阵 IAU 2000B.
    RBPN=pymPnm00b(DATE1,DATE2)
    
    #分离出CIP的X,Y坐标
    X,Y=pymBpn2xy(RBPN)
    
    #基于X,Y坐标，计算CIO定位角s
    S00=pymS00(DATE1,DATE2,X,Y)              
    
    return(S00)


def pymS2p(THETA,PHI,R):
    '''
    Convert spherical polar coordinates to p-vector.

    Parameters
    ----------
    theta : float
        longitude angle (radians)    
    phi : float
        latitude angle (radians)    
    r : float
        radial distance

    Returns
    -------
    p : list(3)
        Cartesian coordinates

    '''
    U=pymS2c(THETA,PHI)
    P=pymSxp(R,U)
    
    return(P)


def pymS2xpv(S1,S2,PV):
    '''
    Multiply a pv-vector by two scalars.

    Parameters
    ----------
    s1 : float
        scalar to multiply position component by    
    s2 : float
        scalar to multiply velocity component by    
    pv : list(2,3)
        pv-vector

    Returns
    -------
    spv : list(2,3)
        pv-vector: p scaled by s1, v scaled by s2

    '''
    SPV=[0,0]
    SPV[0]=pymSxp(S1,PV[0])
    SPV[1]=pymSxp(S2,PV[1])
    
    return(SPV)


def pymS06a(DATE1,DATE2):
    '''
    The CIO locator s, positioning the Celestial Intermediate Origin on
    the equator of the Celestial Intermediate Pole, using the IAU 2006
    precession and IAU 2000A nutation models.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    function value : float
        the CIO locator s in radians

    '''
    #参考架偏差-岁差-章动矩阵, IAU 20006/2000A.
    RBPN=pymPnm06a(DATE1,DATE2)
    
    #从天球中间极（CIP）中分解出X，Y坐标
    X,Y=pymBpn2xy(RBPN)
    
    #基于X,Y坐标，计算CIO定位角s
    S06=pymS06(DATE1,DATE2,X,Y)  
    
    return(S06)


def pymSxpv(S,PV):
    '''
    Multiply a pv-vector by a scalar.

    Parameters
    ----------
    s : float
        scalar    
    pv : list(2,3)
        pv-vector

    Returns
    -------
    spv : list(2,3)
        s * pv

    '''
    SPV=pymS2xpv(S,S,PV)
    
    return(SPV)


def pymTpors( XI,ETA,A,B):
    '''
    In the tangent plane projection, given the rectangular coordinates
    of a star and its spherical coordinates, determine the spherical
    coordinates of the tangent point.

    Parameters
    ----------
    xi : float
        rectangular coordinates of star image    
    eta : float
        rectangular coordinates of star image    
    a : float
        star's spherical coordinates    
    b : float
        star's spherical coordinates

    Returns
    -------
    a01 : float
        tangent point's spherical coordinates, Soln. 1    
    b01 : float
        tangent point's spherical coordinates, Soln. 1    
    a02 : float
        tangent point's spherical coordinates, Soln. 2    
    b02 : float
        tangent point's spherical coordinates, Soln. 2

    '''
    XI2=XI*XI
    R=ma.sqrt(1.0+XI2+ETA*ETA)
    SB=ma.sin(B)
    CB=ma.cos(B)
    RSB=R*SB
    RCB=R*CB
    W2=RCB*RCB-XI2
    if (W2>=0.0):
        W=ma.sqrt(W2)
        S=RSB-ETA*W
        C=RSB*ETA+W
        if (XI==0.0)&(W==0.0):
            W=1.0
        A01=pymAnp(A-ma.atan2(XI,W))
        B01=ma.atan2(S,C)
        W=-W
        S=RSB-ETA*W
        C=RSB*ETA+W
        A02=pymAnp(A-ma.atan2(XI,W))
        B02=ma.atan2(S,C)
        if (np.abs(RSB)<1.0):
            N=1
        else:
            N=2
    else:
        N=0
    
    return(A01,B01,A02,B02,N)


def pymTporv(XI,ETA,V):
    '''
    In the tangent plane projection, given the rectangular coordinates
    of a star and its direction cosines, determine the direction
    cosines of the tangent point.

    Parameters
    ----------
    xi : float
        rectangular coordinates of star image    
    eta : float
        rectangular coordinates of star image     
    v : list(3,3)
        star's direction cosines

    Returns
    -------
    v01 : list(3,3)
        tangent point's direction cosines, Solution 1    
    v02 : list(3,3)
        tangent point's direction cosines, Solution 2    
    N : number of solutions:
        0: 'no solutions returned',
        1: 'only the first solution is useful',
        2: 'both solutions are useful',
    '''
    V01=[0,0,0]
    V02=[0,0,0]
    
    X=V[0]
    Y=V[1]
    Z=V[2]
    RXY2=X*X+Y*Y
    XI2=XI*XI
    ETA2P1=ETA*ETA+1.0
    R=ma.sqrt(XI2+ETA2P1)
    RSB=R*Z
    RCB=R*ma.sqrt(X*X+Y*Y)
    W2=RCB*RCB-XI2
    if (W2>0.0):
        W=ma.sqrt(W2)
        C=(RSB*ETA+W)/(ETA2P1*ma.sqrt(RXY2*(W2+XI2)))
        V01[0]=C*(X*W+Y*XI)
        V01[1]=C*(Y*W-X*XI)
        V01[2]=(RSB-ETA*W)/ETA2P1
        W=-W
        C=(RSB*ETA+W)/(ETA2P1*ma.sqrt(RXY2*(W2+XI2)))
        V02[0]=C*(X*W+Y*XI)
        V02[1]=C*(Y*W-X*XI)
        V02[2]=(RSB-ETA*W)/ETA2P1
        if (np.abs(RSB)<1.0):
            N=1
        else:
            N=2
    else:
        N=0
 
    return(V01,V02,N)


def pymTpsts(XI,ETA,A0,B0):
    '''
    In the tangent plane projection, given the star's rectangular
    coordinates and the spherical coordinates of the tangent point,
    solve for the spherical coordinates of the star.

    Parameters
    ----------
    xi : float
        rectangular coordinates of star image    
    eta : float
        rectangular coordinates of star image        
    a0 : float
        tangent point's spherical coordinates    
    b0 : float
        tangent point's spherical coordinates

    Returns
    -------
    a : float
        star's spherical coordinates    
    b : float
        star's spherical coordinates

    '''
    SB0=ma.sin(B0)
    CB0=ma.cos(B0)
    D=CB0-ETA*SB0
    A=pymAnp(ma.atan2(XI,D)+A0)
    B=ma.atan2(SB0+ETA*CB0,ma.sqrt(XI*XI+D*D))
    
    return(A,B)


def pymTpstv(XI,ETA,V0):
    '''
    In the tangent plane projection, given the star's rectangular
    coordinates and the direction cosines of the tangent point, solve
    for the direction cosines of the star.

    Parameters
    ----------
    xi : float
        rectangular coordinates of star image    
    eta : float
        rectangular coordinates of star image     
    v0 : list(3,3)
        tangent point's direction cosines

    Returns
    -------
    v : list(3,3)
        star's direction cosines

    '''
    V=[0,0,0]
    #T切点
    X=V0[0]
    Y=V0[1]
    Z=V0[2]

    #对在极点的情况处理
    R=ma.sqrt(X*X+Y*Y)
    if (R==0.0):
        R=1e-20
        X=R

    #恒星到切面的向量长度
    F=ma.sqrt(1.0+XI*XI+ETA*ETA)

    #应用转换并归一化
    V[0]=(X-(XI*Y+ETA*X*Z)/R)/F
    V[1]=(Y+(XI*X-ETA*Y*Z)/R)/F
    V[2]=(Z+ETA*R)/F
    
    return(V)


def pymTpxes(A,B,A0,B0):
    '''
    In the tangent plane projection, given celestial spherical
    coordinates for a star and the tangent point, solve for the star's
    rectangular coordinates in the tangent plane.

    Parameters
    ----------
    a : float
        star's spherical coordinates    
    b : float
        star's spherical coordinates    
    a0 : float
        tangent point's spherical coordinates    
    b0 : float
        tangent point's spherical coordinates 

    Returns
    -------
    xi : float
        rectangular coordinates of star image    
    eta : float
        rectangular coordinates of star image
    J : ValueError
        0: 'OK',
        1: 'star too far from axis',
        2: 'antistar on tangent plane',
        3: 'antistar too far from axis',
    '''
    TINY=1e-6
    
    #球面坐标
    SB0=ma.sin(B0)
    SB=ma.sin(B)
    CB0=ma.cos(B0)
    CB=ma.cos(B)
    DA=A-A0
    SDA=ma.sin(DA)
    CDA=ma.cos(DA)

    #恒星向量到切平面的长度的倒数
    D=SB*SB0+CB*CB0*CDA

    #检查是否存在错误
    if (D>TINY):
        J=0    
    elif (D>=0.0):
        J=1
        D=TINY
    elif (D>=-TINY):
        J=2
        D=-TINY
    else:
        J=3

    #返回在切平面上的坐标
    XI=CB*SDA/D
    ETA=(SB*CB0-CB*SB0*CDA)/D
       
    return(XI,ETA,J)


def pymTpxev(V,V0):
    '''
    In the tangent plane projection, given celestial direction cosines
    for a star and the tangent point, solve for the star's rectangular
    coordinates in the tangent plane.

    Parameters
    ----------
    v : list(3,3)
        direction cosines of star    
    v0 : list(3,3)
        direction cosines of tangent point

    Returns
    -------
    xi : float
        tangent plane coordinates of star    
    eta : float
        tangent plane coordinates of star
    J : ValueError
        0: 'OK',
        1: 'star too far from axis',
        2: 'antistar on tangent plane',
        3: 'antistar too far from axis',
    '''
    TINY=1e-6
    
    #恒星以及切点
    X=V[0]
    Y=V[1]
    Z=V[2]
    X0=V0[0]
    Y0=V0[1]
    Z0=V0[2]

    #对在极点的情况进行处理
    R2=X0*X0+Y0*Y0
    R=ma.sqrt(R2)
    if (R==0.0):
        R=1e-20
        X0=R

    #恒星向量到切平面的长度的倒数
    W=X*X0+Y*Y0
    D=W+Z*Z0

    #检查是否存在错误
    if (D>TINY):
        J=0    
    elif (D>=0.0):
        J=1
        D=TINY
    elif (D>=-TINY):
        J=2
        D=-TINY
    else:
        J=3

    #返回在切平面上的坐标
    D=R*D
    XI=(Y*X0-X*Y0)/D
    ETA=(Z*R2-Z0*W)/D

    return(XI,ETA,J)


def pymXys00a(DATE1,DATE2):
    '''
    For a given TT date, compute the X,Y coordinates of the Celestial
    Intermediate Pole and the CIO locator s, using the IAU 2000A
    precession-nutation model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    x : float
        Celestial Intermediate Pole    
    y : float
        Celestial Intermediate Pole    
    s : float
        the CIO locator s
        
    '''
    #调用pymPmN00A给出参考架偏差-岁差-章动矩阵，IAU 2000A.
    RBPN=pymPnm00a(DATE1,DATE2)

    #分离出X,Y坐标
    X,Y=pymBpn2xy(RBPN)

    #获得 s.
    S=pymS00(DATE1,DATE2,X,Y)
    
    return(X,Y,S)


def pymXys00b(DATE1,DATE2):
    '''
    For a given TT date, compute the X,Y coordinates of the Celestial
    Intermediate Pole and the CIO locator s, using the IAU 2000B
    precession-nutation model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    x : float
        Celestial Intermediate Pole    
    y : float
        Celestial Intermediate Pole    
    s : float
        the CIO locator s

    '''
    #调用pymPmN00B给出参考架偏差-岁差-章动矩阵，IAU 2000B.
    RBPN=pymPnm00b(DATE1,DATE2)

    #分离出X,Y坐标
    X,Y=pymBpn2xy(RBPN)

    #获得 s.
    S=pymS00(DATE1,DATE2,X,Y)
    
    return(X,Y,S)


def pymXys06a(DATE1,DATE2):
    '''
    For a given TT date, compute the X,Y coordinates of the Celestial
    Intermediate Pole and the CIO locator s, using the IAU 2006
    precession and IAU 2000A nutation models.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date    
    date2 : float
        TT as a 2-part Julian Date

    Returns
    -------
    x : float
        Celestial Intermediate Pole    
    y : float
        Celestial Intermediate Pole    
    s : float
        the CIO locator s

    '''
    #调用pymPmN06A给出参考架偏差-岁差-章动矩阵，IAU 2006/2000A..
    RBPN=pymPnm06a(DATE1,DATE2)

    #分离出X,Y坐标
    X,Y=pymBpn2xy(RBPN)

    #获得 s.
    S=pymS06(DATE1,DATE2,X,Y)
    
    return(X,Y,S)


def pymZr():
    '''
    Initialize an r-matrix to the null matrix.

    Parameters
    ----------
    r : list(3,3)
        r-matrix

    Returns
    -------
    r : list(3,3)
        null matrix

    '''
    R=[[0 for i in range(3)] for j in range(3)]
    
    return(R)
