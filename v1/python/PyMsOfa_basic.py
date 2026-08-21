import math as ma
import numpy as np

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
    c : list(3)
        direction cosines

    '''
    C=[0,0,0]
    
    CP=ma.cos(PHI)
    C[0]=ma.cos(THETA)*CP
    C[1]=ma.sin(THETA)*CP
    C[2]=ma.sin(PHI)
    
    return(C)
    
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

def pymPxp(A,B):
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
