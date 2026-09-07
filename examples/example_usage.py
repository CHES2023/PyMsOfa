# -*- coding: utf-8 -*-
"""
PyMsOfa 2.x usage examples, with the equivalent 1.x (1.1.6) idiom alongside.

Run it after installing PyMsOfa:

    pip install PyMsOfa
    python examples/example_usage.py
"""

import numpy as np

# ===========================================================================
# 第一版 (1.1.6) 的导入方式：
#     from PyMsOfa import python as sf      # 纯 Python 后端
#     from PyMsOfa import ctypes as sf      # ctypes 后端（Windows 上从 pypi 装完不可用）
#     from PyMsOfa import cffi   as sf
#
# 第二版 (2.0.0)：只有一种，全部平铺在顶层
# ===========================================================================
import PyMsOfa as sf

print("PyMsOfa", sf.__version__, "| SOFA release", sf.SOFA_RELEASE)
print("-" * 70)


# --- 1. 时间系统 ------------------------------------------------------------
print("\n[1] 时间系统")

# 1.x:  djm0, djm, j = sf.pymCal2jd(2003, 6, 1);  if j < 0: ...
# 2.x:  少了状态码 j，出错直接抛 ValueError
djm0, djm = sf.pymCal2jd(2003, 6, 1)
print("  pymCal2jd(2003,6,1)      ->", djm0, djm)

try:
    sf.pymCal2jd(2003, 13, 1)          # 月份非法
except ValueError as err:
    print("  非法月份直接抛异常        ->", err)

# UTC -> TAI -> TT，链式转换。1.x 每一步都多返回一个状态码
utc1, utc2 = 2453750.5, 0.892100694
tai1, tai2 = sf.pymUtctai(utc1, utc2)
tt1, tt2 = sf.pymTaitt(tai1, tai2)
print("  UTC -> TAI               ->", tai1, tai2)
print("  TAI -> TT                ->", tt1, tt2)

# 角度 / 时间格式化
print("  pymD2tf(3, 0.5)          ->", sf.pymD2tf(3, 0.5))
print("  pymA2af(4, 2.345)        ->", sf.pymA2af(4, 2.345))
print("  pymJd2cal(2400000.5,50123.9999) ->", sf.pymJd2cal(2400000.5, 50123.9999))


# --- 2. 向量 / 球面：新版支持数组广播 --------------------------------------
print("\n[2] 向量与球面（2.x 新增：可以直接喂数组）")

print("  标量  pymS2c(0.3, 0.4)         ->", sf.pymS2c(0.3, 0.4))

theta = np.array([0.0, 0.3, 0.6, 0.9])
phi = np.array([0.0, 0.4, 0.8, 1.2])
v = sf.pymS2c(theta, phi)              # 1.x 这里必须写 for 循环
print("  数组  pymS2c(theta, phi).shape ->", v.shape)

th, ph = sf.pymC2s(v)                  # 反变换，同样支持数组
print("  回代  pymC2s(v) 最大误差       ->",
      max(np.abs(sf.pymAnpm(th - theta)).max(), np.abs(ph - phi).max()))

print("  pymAnp(-0.1)                   ->", sf.pymAnp(-0.1))
print("  pymAnp([-0.1, 7.0])            ->", sf.pymAnp(np.array([-0.1, 7.0])))


# --- 3. 岁差 / 章动 / 地球姿态 ---------------------------------------------
print("\n[3] 岁差、章动、地球自转")

tta, ttb = 2400000.5, 53736.0
print("  pymObl06                 ->", sf.pymObl06(tta, ttb))
print("  pymEra00                 ->", sf.pymEra00(2400000.5, 54388.0))
print("  pymGmst06                ->", sf.pymGmst06(2400000.5, 53736.0, 2400000.5, 53736.0))

rnpb = sf.pymPnm06a(tta, ttb)          # 偏置-岁差-章动矩阵
print("  pymPnm06a shape          ->", np.asarray(rnpb).shape)
x, y = sf.pymBpn2xy(rnpb)
print("  pymBpn2xy                ->", x, y)


# --- 4. 天体测量 -----------------------------------------------------------
print("\n[4] 天体测量")

# ICRS -> CIRS。1.x 要先自己造一个 astrom 容器再传进去：
#     astrom = sf.pymASTROM()
#     astrom, eo = sf.pymApci13(date1, date2, astrom)
# 2.x 直接返回：
astrom, eo = sf.pymApci13(2456165.5, 0.401182685)
print("  pymApci13 -> astrom.pmt  ->", astrom.pmt, "| eo =", eo)

ri, di, eo = sf.pymAtci13(2.71, 0.174,
                          -354.45e-3, 595.35e-3, 164.99e-3, 0.0,
                          2456165.5, 0.401182685)
print("  pymAtci13                ->", ri, di, eo)

# 光行差、光线偏折
pnat = [-0.76321968546737951, -0.60869453983060384, -0.21676408580639883]
vel = [2.1044018893653786e-5, -8.9108923304429319e-5, -3.8633714797716569e-5]
print("  pymAb                    ->", sf.pymAb(pnat, vel, 0.99980921395708788,
                                                0.99999999506209258))

# 多天体光线偏折：LDBODY 在 2.x 里是普通 Python 类
body = sf.pymLDBODY(0.00028574, 3e-10,
                    [[-7.81014427, -5.60956681, -1.98079819],
                     [0.0030723249, -0.00406995477, -0.00181335842]])
print("  pymLdn                   ->", sf.pymLdn(1, [body],
                                                 [-0.974170437, -0.2115201, -0.0917583114],
                                                 [-0.763276255, -0.608633767, -0.216735543]))


# --- 5. 大地坐标 -----------------------------------------------------------
print("\n[5] 大地坐标（注意 1.x 纯 Python 后端里这个函数叫 pymGC2GD，2.x 统一成 pymGc2gd）")

elong, phi_, height = sf.pymGc2gd(sf.WGS84, [2e6, 3e6, 5.244e6])
print("  pymGc2gd(WGS84, xyz)     ->", elong, phi_, height)

# 2.x 起：大地测量这一族也统一成"出错抛 ValueError"，不再返回状态码/哨兵值
xyz = sf.pymGd2gc(sf.WGS84, 3.1, -0.5, 2500.0)
print("  pymGd2gc -> xyz           ->", xyz)
a, f = sf.pymEform(sf.WGS84)
print("  pymEform -> a, f          ->", a, f)
try:
    sf.pymGc2gd(4, [2e6, 3e6, 5.244e6])
except ValueError as err:
    print("  非法椭球号 pymGc2gd(4,..) -> ValueError:", str(err)[:45], "...")


# --- 6. 星历 ---------------------------------------------------------------
print("\n[6] 星历")

pvh, pvb = sf.pymEpv00(2400000.5, 53411.52501161)
print("  pymEpv00 日心位置        ->", np.asarray(pvh)[0])

# pymPlan94 也统一成抛异常；非法行星号直接 ValueError
pv = sf.pymPlan94(2400000.5, 43999.9, 1)      # 1 = 水星
print("  pymPlan94(水星) 位置     ->", np.asarray(pv)[0])
try:
    sf.pymPlan94(2400000.5, 43999.9, 12)
except ValueError as err:
    print("  pymPlan94(12)            -> ValueError:", str(err)[:45], "...")

pv = sf.pymMoon98(2400000.5, 43999.9)
print("  pymMoon98 月球位置       ->", np.asarray(pv)[0])


# --- 7. 常量 ---------------------------------------------------------------
print("\n[7] SOFA 常量（原来在 sofam.h 里）")
for name in ["DPI", "D2PI", "DAS2R", "DR2AS", "DJ00", "DJC", "DAYSEC",
             "DAU", "CMPS", "SRS", "WGS84"]:
    print("  sf.%-8s = %r" % (name, getattr(sf, name)))


# --- 8. 想按主题分模块用也可以 ---------------------------------------------
print("\n[8] 分主题导入")
from PyMsOfa import PyMsOfa_time as t
print("  PyMsOfa_time 模块导入成功:", t.pymCal2jd(2024, 2, 27))

print("\n" + "-" * 70)
print("done")
