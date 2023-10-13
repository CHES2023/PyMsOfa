# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 23:01:06 2023
Done    on Mon Aug  6 02:11:37 2023

@author: Dr. Jianghui JI  (jijh@pmo.ac.cn)

@author: PMO
"""
import numpy  as  np
from   cffi  import  FFI
import platform as pf
import sys, os
import warnings  as  ws 

ffi = FFI()


ffi.cdef(
"""
/*sofa.h                                  */

/* Star-independent astrometry parameters */
typedef struct {
   double pmt;        /* PM time interval (SSB, Julian years) */
   double eb[3];      /* SSB to observer (vector, au) */
   double eh[3];      /* Sun to observer (unit vector) */
   double em;         /* distance from Sun to observer (au) */
   double v[3];       /* barycentric observer velocity (vector, c) */
   double bm1;        /* sqrt(1-|v|^2): reciprocal of Lorenz factor */
   double bpn[3][3];  /* bias-precession-nutation matrix */
   double along;      /* longitude + s' + dERA(DUT) (radians) */
   double phi;        /* geodetic latitude (radians) */
   double xpl;        /* polar motion xp wrt local meridian (radians) */
   double ypl;        /* polar motion yp wrt local meridian (radians) */
   double sphi;       /* sine of geodetic latitude */
   double cphi;       /* cosine of geodetic latitude */
   double diurab;     /* magnitude of diurnal aberration vector */
   double eral;       /* "local" Earth rotation angle (radians) */
   double refa;       /* refraction constant A (radians) */
   double refb;       /* refraction constant B (radians) */
} iauASTROM;
/* (Vectors eb, eh, em and v are all with respect to BCRS axes.) */

/* Body parameters for light deflection */
typedef struct {
   double bm;         /* mass of the body (solar masses) */
   double dl;         /* deflection limiter (radians^2/2) */
   double pv[2][3];   /* barycentric PV of the body (au, au/day) */
} iauLDBODY;

/* Astronomy/Calendars */
int    iauCal2jd(int iy, int im, int id, double *djm0, double *djm);
double iauEpb(double dj1, double dj2);
void   iauEpb2jd(double epb, double *djm0, double *djm);
double iauEpj(double dj1, double dj2);
void   iauEpj2jd(double epj, double *djm0, double *djm);
int    iauJd2cal(double dj1, double dj2,
                     int *iy, int *im, int *id, double *fd);
int    iauJdcalf(int ndp, double dj1, double dj2, int iymdf[4]);

/* Astronomy/Astrometry */
void iauAb(double pnat[3], double v[3], double s, double bm1,
           double ppr[3]);
void iauApcg(double date1, double date2,
             double ebpv[2][3], double ehp[3],
             iauASTROM *astrom);
void iauApcg13(double date1, double date2, iauASTROM *astrom);
void iauApci(double date1, double date2,
             double ebpv[2][3], double ehp[3],
             double x, double y, double s,
             iauASTROM *astrom);
void iauApci13(double date1, double date2,
               iauASTROM *astrom, double *eo);
void iauApco(double date1, double date2,
             double ebpv[2][3], double ehp[3],
             double x, double y, double s, double theta,
             double elong, double phi, double hm,
             double xp, double yp, double sp,
             double refa, double refb,
             iauASTROM *astrom);
int iauApco13(double utc1, double utc2, double dut1,
              double elong, double phi, double hm, double xp, double yp,
              double phpa, double tc, double rh, double wl,
              iauASTROM *astrom, double *eo);
void iauApcs(double date1, double date2, double pv[2][3],
             double ebpv[2][3], double ehp[3],
             iauASTROM *astrom);
void iauApcs13(double date1, double date2, double pv[2][3],
               iauASTROM *astrom);
void iauAper(double theta, iauASTROM *astrom);
void iauAper13(double ut11, double ut12, iauASTROM *astrom);
void iauApio(double sp, double theta,
             double elong, double phi, double hm, double xp, double yp,
             double refa, double refb,
             iauASTROM *astrom);
int iauApio13(double utc1, double utc2, double dut1,
              double elong, double phi, double hm, double xp, double yp,
              double phpa, double tc, double rh, double wl,
              iauASTROM *astrom);
void iauAtcc13(double rc, double dc,
               double pr, double pd, double px, double rv,
               double date1, double date2,
               double *ra, double *da);
void iauAtccq(double rc, double dc,
              double pr, double pd, double px, double rv,
              iauASTROM *astrom, double *ra, double *da);
void iauAtci13(double rc, double dc,
               double pr, double pd, double px, double rv,
               double date1, double date2,
               double *ri, double *di, double *eo);
void iauAtciq(double rc, double dc, double pr, double pd,
              double px, double rv, iauASTROM *astrom,
              double *ri, double *di);
void iauAtciqn(double rc, double dc, double pr, double pd,
               double px, double rv, iauASTROM *astrom,
               int n, iauLDBODY b[], double *ri, double *di);
void iauAtciqz(double rc, double dc, iauASTROM *astrom,
               double *ri, double *di);
int iauAtco13(double rc, double dc,
              double pr, double pd, double px, double rv,
              double utc1, double utc2, double dut1,
              double elong, double phi, double hm, double xp, double yp,
              double phpa, double tc, double rh, double wl,
              double *aob, double *zob, double *hob,
              double *dob, double *rob, double *eo);
void iauAtic13(double ri, double di,
               double date1, double date2,
               double *rc, double *dc, double *eo);
void iauAticq(double ri, double di, iauASTROM *astrom,
              double *rc, double *dc);
void iauAticqn(double ri, double di, iauASTROM *astrom,
               int n, iauLDBODY b[], double *rc, double *dc);
int iauAtio13(double ri, double di,
              double utc1, double utc2, double dut1,
              double elong, double phi, double hm, double xp, double yp,
              double phpa, double tc, double rh, double wl,
              double *aob, double *zob, double *hob,
              double *dob, double *rob);
void iauAtioq(double ri, double di, iauASTROM *astrom,
              double *aob, double *zob,
              double *hob, double *dob, double *rob);
int iauAtoc13(const char *type, double ob1, double ob2,
              double utc1, double utc2, double dut1,
              double elong, double phi, double hm, double xp, double yp,
              double phpa, double tc, double rh, double wl,
              double *rc, double *dc);
int iauAtoi13(const char *type, double ob1, double ob2,
              double utc1, double utc2, double dut1,
              double elong, double phi, double hm, double xp, double yp,
              double phpa, double tc, double rh, double wl,
              double *ri, double *di);
void iauAtoiq(const char *type,
              double ob1, double ob2, iauASTROM *astrom,
              double *ri, double *di);
void iauLd(double bm, double p[3], double q[3], double e[3],
           double em, double dlim, double p1[3]);
void iauLdn(int n, iauLDBODY b[], double ob[3], double sc[3],
            double sn[3]);
void iauLdsun(double p[3], double e[3], double em, double p1[3]);
void iauPmpx(double rc, double dc, double pr, double pd,
             double px, double rv, double pmt, double pob[3],
             double pco[3]);
int iauPmsafe(double ra1, double dec1, double pmr1, double pmd1,
              double px1, double rv1,
              double ep1a, double ep1b, double ep2a, double ep2b,
              double *ra2, double *dec2, double *pmr2, double *pmd2,
              double *px2, double *rv2);
void iauPvtob(double elong, double phi, double height, double xp,
              double yp, double sp, double theta, double pv[2][3]);
void iauRefco(double phpa, double tc, double rh, double wl,
              double *refa, double *refb);

/* Astronomy/Ephemerides */
int iauEpv00(double date1, double date2,
             double pvh[2][3], double pvb[2][3]);
void iauMoon98(double date1, double date2, double pv[2][3]);
int iauPlan94(double date1, double date2, int np, double pv[2][3]);

/* Astronomy/FundamentalArgs */
double iauFad03(double t);
double iauFae03(double t);
double iauFaf03(double t);
double iauFaju03(double t);
double iauFal03(double t);
double iauFalp03(double t);
double iauFama03(double t);
double iauFame03(double t);
double iauFane03(double t);
double iauFaom03(double t);
double iauFapa03(double t);
double iauFasa03(double t);
double iauFaur03(double t);
double iauFave03(double t);

/* Astronomy/PrecNutPolar */
void iauBi00(double *dpsibi, double *depsbi, double *dra);
void iauBp00(double date1, double date2,
             double rb[3][3], double rp[3][3], double rbp[3][3]);
void iauBp06(double date1, double date2,
             double rb[3][3], double rp[3][3], double rbp[3][3]);
void iauBpn2xy(double rbpn[3][3], double *x, double *y);
void iauC2i00a(double date1, double date2, double rc2i[3][3]);
void iauC2i00b(double date1, double date2, double rc2i[3][3]);
void iauC2i06a(double date1, double date2, double rc2i[3][3]);
void iauC2ibpn(double date1, double date2, double rbpn[3][3],
               double rc2i[3][3]);
void iauC2ixy(double date1, double date2, double x, double y,
              double rc2i[3][3]);
void iauC2ixys(double x, double y, double s, double rc2i[3][3]);
void iauC2t00a(double tta, double ttb, double uta, double utb,
               double xp, double yp, double rc2t[3][3]);
void iauC2t00b(double tta, double ttb, double uta, double utb,
               double xp, double yp, double rc2t[3][3]);
void iauC2t06a(double tta, double ttb, double uta, double utb,
               double xp, double yp, double rc2t[3][3]);
void iauC2tcio(double rc2i[3][3], double era, double rpom[3][3],
               double rc2t[3][3]);
void iauC2teqx(double rbpn[3][3], double gst, double rpom[3][3],
               double rc2t[3][3]);
void iauC2tpe(double tta, double ttb, double uta, double utb,
              double dpsi, double deps, double xp, double yp,
              double rc2t[3][3]);
void iauC2txy(double tta, double ttb, double uta, double utb,
              double x, double y, double xp, double yp,
              double rc2t[3][3]);
double iauEo06a(double date1, double date2);
double iauEors(double rnpb[3][3], double s);
void iauFw2m(double gamb, double phib, double psi, double eps,
             double r[3][3]);
void iauFw2xy(double gamb, double phib, double psi, double eps,
              double *x, double *y);
void iauLtp(double epj, double rp[3][3]);
void iauLtpb(double epj, double rpb[3][3]);
void iauLtpecl(double epj, double vec[3]);
void iauLtpequ(double epj, double veq[3]);
void iauNum00a(double date1, double date2, double rmatn[3][3]);
void iauNum00b(double date1, double date2, double rmatn[3][3]);
void iauNum06a(double date1, double date2, double rmatn[3][3]);
void iauNumat(double epsa, double dpsi, double deps, double rmatn[3][3]);
void iauNut00a(double date1, double date2, double *dpsi, double *deps);
void iauNut00b(double date1, double date2, double *dpsi, double *deps);
void iauNut06a(double date1, double date2, double *dpsi, double *deps);
void iauNut80(double date1, double date2, double *dpsi, double *deps);
void iauNutm80(double date1, double date2, double rmatn[3][3]);
double iauObl06(double date1, double date2);
double iauObl80(double date1, double date2);
void iauP06e(double date1, double date2,
             double *eps0, double *psia, double *oma, double *bpa,
             double *bqa, double *pia, double *bpia,
             double *epsa, double *chia, double *za, double *zetaa,
             double *thetaa, double *pa,
             double *gam, double *phi, double *psi);
void iauPb06(double date1, double date2,
             double *bzeta, double *bz, double *btheta);
void iauPfw06(double date1, double date2,
              double *gamb, double *phib, double *psib, double *epsa);
void iauPmat00(double date1, double date2, double rbp[3][3]);
void iauPmat06(double date1, double date2, double rbp[3][3]);
void iauPmat76(double date1, double date2, double rmatp[3][3]);
void iauPn00(double date1, double date2, double dpsi, double deps,
             double *epsa,
             double rb[3][3], double rp[3][3], double rbp[3][3],
             double rn[3][3], double rbpn[3][3]);
void iauPn00a(double date1, double date2,
              double *dpsi, double *deps, double *epsa,
              double rb[3][3], double rp[3][3], double rbp[3][3],
              double rn[3][3], double rbpn[3][3]);
void iauPn00b(double date1, double date2,
              double *dpsi, double *deps, double *epsa,
              double rb[3][3], double rp[3][3], double rbp[3][3],
              double rn[3][3], double rbpn[3][3]);
void iauPn06(double date1, double date2, double dpsi, double deps,
             double *epsa,
             double rb[3][3], double rp[3][3], double rbp[3][3],
             double rn[3][3], double rbpn[3][3]);
void iauPn06a(double date1, double date2,
              double *dpsi, double *deps, double *epsa,
              double rb[3][3], double rp[3][3], double rbp[3][3],
              double rn[3][3], double rbpn[3][3]);
void iauPnm00a(double date1, double date2, double rbpn[3][3]);
void iauPnm00b(double date1, double date2, double rbpn[3][3]);
void iauPnm06a(double date1, double date2, double rnpb[3][3]);
void iauPnm80(double date1, double date2, double rmatpn[3][3]);
void iauPom00(double xp, double yp, double sp, double rpom[3][3]);
void iauPr00(double date1, double date2,
             double *dpsipr, double *depspr);
void iauPrec76(double date01, double date02,
               double date11, double date12,
               double *zeta, double *z, double *theta);
double iauS00(double date1, double date2, double x, double y);
double iauS00a(double date1, double date2);
double iauS00b(double date1, double date2);
double iauS06(double date1, double date2, double x, double y);
double iauS06a(double date1, double date2);
double iauSp00(double date1, double date2);
void iauXy06(double date1, double date2, double *x, double *y);
void iauXys00a(double date1, double date2,
               double *x, double *y, double *s);
void iauXys00b(double date1, double date2,
               double *x, double *y, double *s);
void iauXys06a(double date1, double date2,
               double *x, double *y, double *s);

/* Astronomy/Rotation And Time */
double iauEe00(double date1, double date2, double epsa, double dpsi);
double iauEe00a(double date1, double date2);
double iauEe00b(double date1, double date2);
double iauEe06a(double date1, double date2);
double iauEect00(double date1, double date2);
double iauEqeq94(double date1, double date2);
double iauEra00(double dj1, double dj2);
double iauGmst00(double uta, double utb, double tta, double ttb);
double iauGmst06(double uta, double utb, double tta, double ttb);
double iauGmst82(double dj1, double dj2);
double iauGst00a(double uta, double utb, double tta, double ttb);
double iauGst00b(double uta, double utb);
double iauGst06(double uta, double utb, double tta, double ttb,
                double rnpb[3][3]);
double iauGst06a(double uta, double utb, double tta, double ttb);
double iauGst94(double uta, double utb);

/* Astronomy/SpaceMotion */
int iauPvstar(double pv[2][3], double *ra, double *dec,
              double *pmr, double *pmd, double *px, double *rv);
int iauStarpv(double ra, double dec,
              double pmr, double pmd, double px, double rv,
              double pv[2][3]);

/* Astronomy/StarCatalogs */

void iauFk425(double r1950, double d1950,
              double dr1950, double dd1950,
              double p1950, double v1950,
              double *r2000, double *d2000,
              double *dr2000, double *dd2000,
              double *p2000, double *v2000);
void iauFk45z(double r1950, double d1950, double bepoch,
              double *r2000, double *d2000);
void iauFk524(double r2000, double d2000,
              double dr2000, double dd2000,
              double p2000, double v2000,
              double *r1950, double *d1950,
              double *dr1950, double *dd1950,
              double *p1950, double *v1950);
void iauFk52h(double r5, double d5,
              double dr5, double dd5, double px5, double rv5,
              double *rh, double *dh,
              double *drh, double *ddh, double *pxh, double *rvh);
void iauFk54z(double r2000, double d2000, double bepoch,
              double *r1950, double *d1950,
              double *dr1950, double *dd1950);
void iauFk5hip(double r5h[3][3], double s5h[3]);
void iauFk5hz(double r5, double d5, double date1, double date2,
              double *rh, double *dh);
void iauH2fk5(double rh, double dh,
              double drh, double ddh, double pxh, double rvh,
              double *r5, double *d5,
              double *dr5, double *dd5, double *px5, double *rv5);
void iauHfk5z(double rh, double dh, double date1, double date2,
              double *r5, double *d5, double *dr5, double *dd5);
int iauStarpm(double ra1, double dec1,
              double pmr1, double pmd1, double px1, double rv1,
              double ep1a, double ep1b, double ep2a, double ep2b,
              double *ra2, double *dec2,
              double *pmr2, double *pmd2, double *px2, double *rv2);

/* Astronomy/EclipticCoordinates */
void iauEceq06(double date1, double date2, double dl, double db,
               double *dr, double *dd);
void iauEcm06(double date1, double date2, double rm[3][3]);
void iauEqec06(double date1, double date2, double dr, double dd,
               double *dl, double *db);
void iauLteceq(double epj, double dl, double db, double *dr, double *dd);
void iauLtecm(double epj, double rm[3][3]);
void iauLteqec(double epj, double dr, double dd, double *dl, double *db);

/* Astronomy/GalacticCoordinates */
void iauG2icrs(double dl, double db, double *dr, double *dd);
void iauIcrs2g(double dr, double dd, double *dl, double *db);

/* Astronomy/GeodeticGeocentric */
int iauEform(int n, double *a, double *f);
int iauGc2gd(int n, double xyz[3],
             double *elong, double *phi, double *height);
int iauGc2gde(double a, double f, double xyz[3],
              double *elong, double *phi, double *height);
int iauGd2gc(int n, double elong, double phi, double height,
             double xyz[3]);
int iauGd2gce(double a, double f,
              double elong, double phi, double height, double xyz[3]);

/* Astronomy/Timescales */
int iauD2dtf(const char *scale, int ndp, double d1, double d2,
             int *iy, int *im, int *id, int ihmsf[4]);
int iauDat(int iy, int im, int id, double fd, double *deltat);
double iauDtdb(double date1, double date2,
               double ut, double elong, double u, double v);
int iauDtf2d(const char *scale, int iy, int im, int id,
             int ihr, int imn, double sec, double *d1, double *d2);
int iauTaitt(double tai1, double tai2, double *tt1, double *tt2);
int iauTaiut1(double tai1, double tai2, double dta,
              double *ut11, double *ut12);
int iauTaiutc(double tai1, double tai2, double *utc1, double *utc2);
int iauTcbtdb(double tcb1, double tcb2, double *tdb1, double *tdb2);
int iauTcgtt(double tcg1, double tcg2, double *tt1, double *tt2);
int iauTdbtcb(double tdb1, double tdb2, double *tcb1, double *tcb2);
int iauTdbtt(double tdb1, double tdb2, double dtr,
             double *tt1, double *tt2);
int iauTttai(double tt1, double tt2, double *tai1, double *tai2);
int iauTttcg(double tt1, double tt2, double *tcg1, double *tcg2);
int iauTttdb(double tt1, double tt2, double dtr,
             double *tdb1, double *tdb2);
int iauTtut1(double tt1, double tt2, double dt,
             double *ut11, double *ut12);
int iauUt1tai(double ut11, double ut12, double dta,
              double *tai1, double *tai2);
int iauUt1tt(double ut11, double ut12, double dt,
             double *tt1, double *tt2);
int iauUt1utc(double ut11, double ut12, double dut1,
              double *utc1, double *utc2);
int iauUtctai(double utc1, double utc2, double *tai1, double *tai2);
int iauUtcut1(double utc1, double utc2, double dut1,
              double *ut11, double *ut12);

/* Astronomy/HorizonEquatorial */
void iauAe2hd(double az, double el, double phi,
              double *ha, double *dec);
void iauHd2ae(double ha, double dec, double phi,
              double *az, double *el);
double iauHd2pa(double ha, double dec, double phi);

/* Astronomy/Gnomonic */
int iauTpors(double xi, double eta, double a, double b,
             double *a01, double *b01, double *a02, double *b02);
int iauTporv(double xi, double eta, double v[3],
             double v01[3], double v02[3]);
void iauTpsts(double xi, double eta, double a0, double b0,
              double *a, double *b);
void iauTpstv(double xi, double eta, double v0[3], double v[3]);
int iauTpxes(double a, double b, double a0, double b0,
             double *xi, double *eta);
int iauTpxev(double v[3], double v0[3], double *xi, double *eta);

/* VectorMatrix/AngleOps */
void iauA2af(int ndp, double angle, char *sign, int idmsf[4]);
void iauA2tf(int ndp, double angle, char *sign, int ihmsf[4]);
int iauAf2a(char s, int ideg, int iamin, double asec, double *rad);
double iauAnp(double a);
double iauAnpm(double a);
void iauD2tf(int ndp, double days, char *sign, int ihmsf[4]);
int iauTf2a(char s, int ihour, int imin, double sec, double *rad);
int iauTf2d(char s, int ihour, int imin, double sec, double *days);

/* VectorMatrix/BuildRotations */
void iauRx(double phi, double r[3][3]);
void iauRy(double theta, double r[3][3]);
void iauRz(double psi, double r[3][3]);

/* VectorMatrix/CopyExtendExtract */
void iauCp(double p[3], double c[3]);
void iauCpv(double pv[2][3], double c[2][3]);
void iauCr(double r[3][3], double c[3][3]);
void iauP2pv(double p[3], double pv[2][3]);
void iauPv2p(double pv[2][3], double p[3]);

/* VectorMatrix/Initialization */
void iauIr(double r[3][3]);
void iauZp(double p[3]);
void iauZpv(double pv[2][3]);
void iauZr(double r[3][3]);

/* VectorMatrix/MatrixOps */
void iauRxr(double a[3][3], double b[3][3], double atb[3][3]);
void iauTr(double r[3][3], double rt[3][3]);

/* VectorMatrix/MatrixVectorProducts */
void iauRxp(double r[3][3], double p[3], double rp[3]);
void iauRxpv(double r[3][3], double pv[2][3], double rpv[2][3]);
void iauTrxp(double r[3][3], double p[3], double trp[3]);
void iauTrxpv(double r[3][3], double pv[2][3], double trpv[2][3]);

/* VectorMatrix/RotationVectors */
void iauRm2v(double r[3][3], double w[3]);
void iauRv2m(double w[3], double r[3][3]);

/* VectorMatrix/SeparationAndAngle */
double iauPap(double a[3], double b[3]);
double iauPas(double al, double ap, double bl, double bp);
double iauSepp(double a[3], double b[3]);
double iauSeps(double al, double ap, double bl, double bp);

/* VectorMatrix/SphericalCartesian */
void iauC2s(double p[3], double *theta, double *phi);
void iauP2s(double p[3], double *theta, double *phi, double *r);
void iauPv2s(double pv[2][3],
             double *theta, double *phi, double *r,
             double *td, double *pd, double *rd);
void iauS2c(double theta, double phi, double c[3]);
void iauS2p(double theta, double phi, double r, double p[3]);
void iauS2pv(double theta, double phi, double r,
             double td, double pd, double rd,
             double pv[2][3]);

/* VectorMatrix/VectorOps */
double iauPdp(double a[3], double b[3]);
double iauPm(double p[3]);
void iauPmp(double a[3], double b[3], double amb[3]);
void iauPn(double p[3], double *r, double u[3]);
void iauPpp(double a[3], double b[3], double apb[3]);
void iauPpsp(double a[3], double s, double b[3], double apsb[3]);
void iauPvdpv(double a[2][3], double b[2][3], double adb[2]);
void iauPvm(double pv[2][3], double *r, double *s);
void iauPvmpv(double a[2][3], double b[2][3], double amb[2][3]);
void iauPvppv(double a[2][3], double b[2][3], double apb[2][3]);
void iauPvu(double dt, double pv[2][3], double upv[2][3]);
void iauPvup(double dt, double pv[2][3], double p[3]);
void iauPvxpv(double a[2][3], double b[2][3], double axb[2][3]);
void iauPxp(double a[3], double b[3], double axb[3]);
void iauS2xpv(double s1, double s2, double pv[2][3], double spv[2][3]);
void iauSxp(double s, double p[3], double sp[3]);
void iauSxpv(double s, double pv[2][3], double spv[2][3]);

""")

'''
   Descr: return file name/path information
#      Input: 
          filepath: full path of a file

#      Return values:
          dirname: directory name of a file
          fullname: full name of a file
          filename: file name of a file
          extension: file extension
'''
def getFileName(filepath):

    __file__  = filepath
    
    dirname   = os.path.dirname(__file__)
    fullname  = os.path.split(__file__)[-1]
    filename  = os.path.split(__file__)[-1].split('.')[0]
    extension = os.path.split(__file__)[-1].split('.')[1]
    
    return dirname, fullname, filename, extension

     
def getFilePath(filepath):

    __file__  = filepath
    
    dirname   = os.path.dirname(__file__)
    fullname  = os.path.split(__file__)[-1]
    filename  = os.path.split(__file__)[-1].split('.')[0]
    extension = os.path.split(__file__)[-1].split('.')[1]
    
    return dirname 

#place libsofa_c.so in the same directory with PyMsOfa.py
cur_dir    = sys.argv[0]
#dirname, fullname, filename, extension = getFileName(cur_dir)
dirname    = getFilePath(cur_dir)
libs_file  = 'libsofa_c.so'
libs_path  = os.path.join(dirname, libs_file) 


#'''
#for windows
if   pf.system().lower() == 'windows':
           
          #user_path  = './'
          #libs_file  = 'libsofa_c.so'
          #libs_path  = user_path + libs_file
          #print(libs_path)

          lib    = ffi.dlopen(libs_path)
           
#for linux 
elif pf.system().lower() == 'linux':  
          dirname = os.getcwd()
          lib     = ffi.dlopen(os.path.join(dirname, libs_file)) 
 
#for macOS
elif pf.system().lower() == 'darwin':  

          lib  = ffi.dlopen(libs_path)



#2023-08-02
def double_ptr():
     
    return ffi.new("double *")

#ctypes, double_ptr similar to c_double in ctypes
def c_double():
 
    return ffi.new("double *")

def c_int():
    
    return ffi.new("int *")

def c_char():
     
    return ffi.new("char *")

'''
def c_double():
    dbl_ptr   = ffi.new("double *")
    return dbl_ptr
    
def double_ptr():
    dbl_ptr   = ffi.new("double *")
    return dbl_ptr

def int_ptr():
    int_ptr   = ffi.new("int *")
    return int_ptr 

def char_ptr():
    char_ptr   = ffi.new("char *")
    return char_ptr 
    
'''

def int_ptr():
    #int_ptr   = ffi.new("int *")
    return ffi.new("int *")

def char_ptr():
    #char_ptr   = ffi.new("char *")
    return ffi.new("char *")

def int_array_4():
    
    idmsf   = np.array(np.zeros(shape=(1,4), dtype=int, order='C'))
    idmsf_p = ffi.cast("int [4]",  ffi.from_buffer(idmsf))
    
    return idmsf, idmsf_p 


def pymASTROM():  
     
    return ffi.new("iauASTROM *")

'''
typedef struct {
   double pmt;        /* PM time interval (SSB, Julian years) */
   double eb[3];      /* SSB to observer (vector, au) */
   double eh[3];      /* Sun to observer (unit vector) */
   double em;         /* distance from Sun to observer (au) */
   double v[3];       /* barycentric observer velocity (vector, c) */
   double bm1;        /* sqrt(1-|v|^2): reciprocal of Lorenz factor */
   double bpn[3][3];  /* bias-precession-nutation matrix */
   double along;      /* longitude + s' + dERA(DUT) (radians) */
   double phi;        /* geodetic latitude (radians) */
   double xpl;        /* polar motion xp wrt local meridian (radians) */
   double ypl;        /* polar motion yp wrt local meridian (radians) */
   double sphi;       /* sine of geodetic latitude */
   double cphi;       /* cosine of geodetic latitude */
   double diurab;     /* magnitude of diurnal aberration vector */
   double eral;       /* "local" Earth rotation angle (radians) */
   double refa;       /* refraction constant A (radians) */
   double refb;       /* refraction constant B (radians) */
'''

def pymASTROM_print(a):
    
        eb_c = tuple([x for x in a.eb])
        eh_c = tuple([x for x in a.eh])
        v_c  = tuple([x for x in a.v])
       
        #bpn_c= tuple([x for x in a.bpn[0]])
        #bpn_c= tuple([x for x in a.bpn[1]])
        #bpn_c= tuple([x for x in a.bpn[2]])   
        
        #eb_c    = np.array([x for x in a.eb])
        #eh_c    = np.array([x for x in a.eh])
        #v_c     = np.array([x for x in a.v])
        
        bpn_c0  = np.array([x for x in a.bpn[0]] ).tolist()
        bpn_c1  = np.array([x for x in a.bpn[1]] ).tolist()
        bpn_c2  = np.array([x for x in a.bpn[2]] ).tolist()       
       
        bpn_c  = [bpn_c0, bpn_c1, bpn_c2]
        print("pmt:",    a.pmt)
        print("eb:",     eb_c)
        print("eh:",     eh_c)
        print("em:",     a.em) 
        print("v:",      v_c)
        print("bm1:",    a.bm1) 
        print("bpn:",    bpn_c)
        print("along:",  a.along)   
        print("phi:",    a.phi)    
        print("xpl:",    a.xpl)     
        print("ypl:",    a.ypl)     
        print("sphi:",   a.sphi)    
        print("cphi:",   a.cphi)    
        print("diurab:", a.diurab) 
        print("eral:",   a.eral)  
        print("refa:",   a.refa)     
        print("refb:",   a.refb)
        return 
 
def pymLDBODY(c):

    c_vec = ffi.cast("iauLDBODY *",  ffi.from_buffer(c))
    
    return c_vec     

    
#for input arrays passing values
def vector_double(c):

    c_vec = ffi.cast("double [3]",  ffi.from_buffer(c))
    
    return c_vec 

def vector_double2(c):

    c_vec = ffi.cast("double [2][3]",  ffi.from_buffer(c))
    
    return c_vec 

def vector_double3(c):

    c_vec = ffi.cast("double [3][3]",  ffi.from_buffer(c))
    
    return c_vec 


#for output arrays returning values
#vector_double   = nt.ndpointer(shape=(1,3), dtype=np.double, flags='C') 
def vector_double_ptr():
    
    c     = np.array(np.zeros(shape=(1,3), dtype=float, order='C'))
    c_vec = ffi.cast("double [3]",  ffi.from_buffer(c))
    
    return c, c_vec 

def vector_double2_ptr():
    
    c     = np.array(np.zeros(shape=(2,3), dtype=float, order='C'))
    c_vec = ffi.cast("double [2][3]",  ffi.from_buffer(c))
    
    return c, c_vec 

def vector_double3_ptr():
    
    c     = np.array(np.zeros(shape=(3,3), dtype=float, order='C'))
    c_vec = ffi.cast("double [3][3]",  ffi.from_buffer(c))
    
    return c, c_vec 

def array_double_ptr():
    
    c     = np.array(np.zeros(shape=(1,2), dtype=float, order='C'))
    c_vec = ffi.cast("double [2]",  ffi.from_buffer(c))
    
    return c, c_vec 

#####

def pymPas(al, ap, bl, bp):
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
    function value : float
        position angle of B with respect to A

    '''
    return lib.iauPas(al, ap, bl, bp)


def pymS2c(theta, phi):
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
    #c     = np.array(np.zeros(shape=(1,3), dtype=float, order='C'))
    #c_vec = ffi.cast("double [3]",  ffi.from_buffer(c))
 
    c, c_vec = vector_double_ptr()
    
    lib.iauS2c(theta, phi, c_vec)
    
    #c  = list(np.array(c).flatten())
    return c 



#int iauCal2jd(int iy, int im, int id, double *djm0, double *djm)
iau_cal2jd_msg = {
                 -1: 'bad year,  the year is simply valid from -4800 March 1',
                 -2: 'bad month, the month is not from 1 to 12',
                 -3: 'bad day,   the day is not related to the month'
                 } 

def pymCal2jd(iyear, imon, iday):
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

    Raises
    ------
    ValueError
        -1: 'bad year,  the year is simply valid from -4800 March 1',
        -2: 'bad month, the month is not from 1 to 12',
        -3: 'bad day,   the day is not related to the month'

    Returns
    -------
    djm0 : float
        MJD zero-point: always 2400000.5
    djm : float
        Modified Julian Date for 0 hrs

    '''
    djm0_p  = double_ptr()
    djm_p   = double_ptr()
   
    j       = lib.iauCal2jd(iyear, imon, iday, djm0_p, djm_p)
   
    if j != 0:
        raise ValueError(iau_cal2jd_msg[j]) 
        
    djm0   =  djm0_p[0]
    djm    =   djm_p[0]
        
    return  djm0, djm

#void iauA2af(int ndp, double angle, char *sign, int idmsf[4])

def pymA2af(ndp, angle):
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
    sign : bytes
        '+' or '-'    
    idmsf : tuple
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
   
    sign_p  = char_ptr()
    
    #idmsf   = np.array(np.zeros(shape=(1,4), dtype=int, order='C'))
    #idmsf_p = ffi.cast("int [4]",  ffi.from_buffer(idmsf))
   
    idmsf, idmsf_p = int_array_4()
    
    lib.iauA2af(ndp, angle, sign_p, idmsf_p)
    
    sign   = sign_p[0]
    idmsf  = list(np.array(idmsf).flatten())
    
    return sign, idmsf 


def pymA2af_a(ndp, angle):
   
    sign_p  = char_ptr()
    
    idmsf   = np.array(np.zeros(shape=(1,4), dtype=int, order='C'))
    idmsf_p = ffi.cast("int [4]",  ffi.from_buffer(idmsf))
    
    lib.iauA2af(ndp, angle, sign_p, idmsf_p)
    
    sign   = sign_p[0]
    idmsf  = list(np.array(idmsf).flatten())
    
    return sign, idmsf 


#int iauJd2cal(double dj1, double dj2,
#                     int *iy, int *im, int *id, double *fd);

iau_jd2cal_msg ={
                -1: 'The valid date is -68569.5 (-4900 March 1) up to 1e9'
                }  

def pymJd2cal(dj1, dj2):
    '''
    Julian Date to Gregorian year, month, day, and fraction of a day.

    Parameters
    ----------
    dj1 : float
        Julian Date    
    dj2 : float
        Julian Date

    Raises
    ------
    ValueError
        -1: 'The valid date is -68569.5 (-4900 March 1) up to 1e9'

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

    '''
    iy_p   = int_ptr()
    im_p   = int_ptr()
    iday_p = int_ptr()
    fday_p = double_ptr()
    
    j = lib.iauJd2cal(dj1, dj2, iy_p, im_p, iday_p, fday_p)

    if j != 0:
        raise ValueError(iau_jd2cal_msg[j])

    iy   = iy_p[0]
    im   = im_p[0]
    iday = iday_p[0]
    fday = fday_p[0]
    
    return iy, im, iday, fday 

#int iauJdcalf(int ndp, double dj1, double dj2, int iymdf[4])


def pymJdcalf(ndp, dj1, dj2):
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
    iymdf : tuple
        year, month, day, fraction in Gregorian calendar

    '''
    #iymdf   = np.array(np.zeros(shape=(1,4), dtype=int, order='C'))
    #iymdf_p = ffi.cast("int [4]",  ffi.from_buffer(iymdf))
   
    iymdf, iymdf_p = int_array_4()
    
    lib.iauJdcalf(ndp, dj1,  dj2, iymdf_p)
    
    iymdf  = list(np.array(iymdf).flatten())
    
    return iymdf


#void iauD2tf(int ndp, double days, char *sign, int ihmsf[4])
def pymD2tf(ndp, days):
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
    sign : 'bytes'
        '+' or '-'    
    ihmsf : tuple    
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
    sign_p = char_ptr()
    
    #ihmsf   = np.array(np.zeros(shape=(1,4), dtype=int, order='C'))
    #ihmsf_p = ffi.cast("int [4]",  ffi.from_buffer(ihmsf))
    
    ihmsf, ihmsf_p = int_array_4()
    lib.iauD2tf(ndp, days, sign_p, ihmsf_p)
    
    sign   = sign_p[0]
    ihmsf  = list(np.array(ihmsf).flatten())
    
    return sign, ihmsf


#2023-08-03
#int iauD2dtf(const char *scale, int ndp, double d1, double d2,
#             int *iy, int *im, int *id, int ihmsf[4])
iau_d2dtf_msg = {
                 1: 'dubious year',
                -1: 'unacceptable date',
                }

def pymD2dtf(scale, ndp, d1, d2):
    '''
    Format for output a 2-part Julian Date (or in the case of UTC a 
    quasi-JD form that includes special provision for leap seconds).

    Parameters
    ----------
    scale : ctypes.c_char_p
        time scale ID(Only the value "UTC" is significant)       
    ndp : int
        resolution    
    d1 : float
        time as a 2-part Julian Date    
    d2 : float
        time as a 2-part Julian Date

    Raises
    ------
    ValueError
        1: 'dubious year',
       -1: 'unacceptable date',

    Returns
    -------
    iy,im,id : int
        year, month, day in Gregorian calendar    
    ihmsf : tuple
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
    iy_p   = int_ptr()
    im_p   = int_ptr()
    iday_p = int_ptr()
    
    ihmsf, ihmsf_p = int_array_4()
    
    j=lib.iauD2dtf(scale, ndp, d1, d2, iy_p, im_p, iday_p, ihmsf_p)
    
    if j < 0:
         raise ValueError(iau_d2dtf_msg[j])
    elif j > 0:
         ws.warn(iau_d2dtf_msg[j], UserWarning, 2)
         
    iy   = iy_p[0]
    im   = im_p[0]
    iday = iday_p[0]     
    
    ihmsf  = list(np.array(ihmsf).flatten())
    
    return (iy, im, iday) + tuple([x for x in ihmsf])


#int iauDat(int iy, int im, int id, double fd, double *deltat)
#Note: current leap second is 2021, this should update and recompile iauDat for new date    
iau_dat_msg =   {
                  1: 'dubious year', 
                 -1: 'bad year,  the year is simply valid from -4800 March 1',
                 -2: 'bad month, the month is not from 1 to 12',
                 -3: 'bad day,   the day is not related to the month',
                 -4: 'bad fraction of day',
                 -5: 'internal error', 
                 } 

def pymDat(iy, im, iday, fday):
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

    Raises
    ------
    ValueError
        1: 'dubious year', 
       -1: 'bad year,  the year is simply valid from -4800 March 1',
       -2: 'bad month, the month is not from 1 to 12',
       -3: 'bad day,   the day is not related to the month',
       -4: 'bad fraction of day',
       -5: 'internal error', 

    Returns
    -------
    deltat : float
        TAI minus UTC, seconds

    '''
    deltat_p = double_ptr()
    
    j = lib.iauDat(iy, im, iday, fday, deltat_p)
    
    if   j < 0:
          raise ValueError(iau_dat_msg[j])
    elif j > 0:
          ws.warn(iau_dat_msg[j], UserWarning, 2)    
          
    deltat = deltat_p[0]
    
    return deltat  

#double iauDtdb(double date1, double date2,
#               double ut, double elong, double u, double v)   
def pymDtdb(date1, date2, ut, elong, u, v):
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
    function value : float
        TDB-TT (seconds)

    '''    
    return lib.iauDtdb(date1, date2, ut, elong, u, v)


#int iauDtf2d(const char *scale, int iy, int im, int id,           
#             int ihr, int imn, double sec, double *d1, double *d2)
iau_dtf2d_msg = {
                  3: 'both of next two',
                  2: 'time is after end of day',
                  1: 'dubious year', 
                 -1: 'bad year,  the year is simply valid from -4800 March 1',
                 -2: 'bad month, the month is not from 1 to 12',
                 -3: 'bad day,   the day is not related to the month',
                 -4: 'bad hour',
                 -5: 'bad minute',
                 -6: 'bad second', 
                } 

def pymDtf2d(scale,  iy,  im,  id,  ihr, imn,  sec):
    '''
    Encode date and time fields into 2-part Julian Date (or in the case
    of UTC a quasi-JD form that includes special provision for leap seconds).

    Parameters
    ----------
    scale : ctypes.c_char_p
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

    Raises
    ------
    ValueError
        3: 'both of next two',
        2: 'time is after end of day',
        1: 'dubious year', 
       -1: 'bad year,  the year is simply valid from -4800 March 1',
       -2: 'bad month, the month is not from 1 to 12',
       -3: 'bad day,   the day is not related to the month',
       -4: 'bad hour',
       -5: 'bad minute',
       -6: 'bad second', 

    Returns
    -------
    d1,d2 : float
        2-part Julian Date

    '''
    d1_p   = double_ptr()
    d2_p   = double_ptr()
   
    j      = lib.iauDtf2d(scale, iy, im, id, ihr, imn, sec, d1_p, d2_p)
    
    if j < 0:
         raise ValueError(iau_dtf2d_msg[j])
    elif j > 0:
         ws.warn(iau_dtf2d_msg[j], UserWarning, 2)
         
    d1 = d1_p[0]
    d2 = d2_p[0]
    
    return d1, d2

#double iauEpb(double dj1, double dj2)
def pymEpb(dj1, dj2):
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
    return lib.iauEpb(dj1, dj2)


#void iauEpb2jd(double epb, double *djm0, double *djm)
def pymEpb2jd(epb):
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
    djm0_p  = double_ptr()
    djm_p   = double_ptr()
    
    lib.iauEpb2jd(epb, djm0_p, djm_p)
    
    djm0   =  djm0_p[0]
    djm    =   djm_p[0]
    
    return djm0, djm

#double iauEpj(double dj1, double dj2)    
def pymEpj(dj1, dj2):
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
    return lib.iauEpj(dj1, dj2)  

#void iauEpj2jd(double epj, double *djm0, double *djm)
def pymEpj2jd(epj):
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
    djm0_p  = double_ptr()
    djm_p   = double_ptr()
    
    lib.iauEpj2jd(epj, djm0_p, djm_p)
    
    djm0   =  djm0_p[0]
    djm    =   djm_p[0]
    
    return djm0, djm

#int iauTaitt(double tai1, double tai2, double *tt1, double *tt2)
def pymTaitt(tai1, tai2):
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

    '''
    tt1_p  = double_ptr()
    tt2_p  = double_ptr()
    
    j = lib.iauTaitt(tai1, tai2, tt1_p, tt2_p)
    
    tt1 = tt1_p[0]  
    tt2 = tt2_p[0]
    
    return tt1, tt2


#int iauTaiut1(double tai1, double tai2, double dta,
#              double *ut11, double *ut12)
    
def pymTaiut1(tai1, tai2, dta):  
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

    '''
    ut11_p  = double_ptr()
    ut12_p  = double_ptr()
    
    j = lib.iauTaiut1(tai1, tai2, dta, ut11_p, ut12_p)
    
    ut11  = ut11_p[0]
    ut12  = ut12_p[0]
    
    return ut11, ut12


#int iauTaiutc(double tai1, double tai2, double *utc1, double *utc2)
iau_taiutc_msg = {
                  1: 'dubious year',
                 -1: 'unacceptable date',
                  }    
def pymTaiutc(tai1, tai2):  
    '''
    Time scale transformation:  International Atomic Time, TAI, to
    Coordinated Universal Time, UTC.

    Parameters
    ----------
    tai1 : float
        TAI as a 2-part Julian Date
    tai2 : float
        TAI as a 2-part Julian Date

    Raises
    ------
    ValueError
        1: 'dubious year',
       -1: 'unacceptable date',

    Returns
    -------
    utc1 : float
        UTC as a 2-part quasi Julian Date
    utc2 : float
        UTC as a 2-part quasi Julian Date

    '''
    utc1_p  = double_ptr()
    utc2_p  = double_ptr()
    
    j = lib.iauTaiutc(tai1, tai2, utc1_p, utc2_p)
    if   j < 0:
        raise ValueError(iau_taiutc_msg[j])
    elif j > 0:
        ws.warn(iau_taiutc_msg[j], UserWarning, 2)
    
    utc1 = utc1_p [0]
    utc2 = utc2_p [0]
    
    return utc1, utc2

#int iauTcbtdb(double tcb1, double tcb2, double *tdb1, double *tdb2)
def pymTcbtdb(tcb1, tcb2):  
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

    '''
    tdb1_p  = double_ptr()
    tdb2_p  = double_ptr()
    j       = lib.iauTcbtdb(tcb1, tcb2, tdb1_p, tdb2_p)
    
    tdb1    = tdb1_p[0]
    tdb2    = tdb2_p[0]
    
    return tdb1, tdb2 


#int iauTcgtt(double tcg1, double tcg2, double *tt1, double *tt2) 
def pymTcgtt(tcg1, tcg2):  
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

    '''    
    tt1_p  = double_ptr()
    tt2_p  = double_ptr()
    j      = lib.iauTcgtt(tcg1, tcg2,  tt1_p,  tt2_p)
    
    tt1    =  tt1_p[0]
    tt2    =  tt2_p[0]

    return tt1, tt2   


#int iauTdbtcb(double tdb1, double tdb2, double *tcb1, double *tcb2)
def pymTdbtcb(tdb1, tdb2):  
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

    '''
    tcb1_p  = double_ptr()
    tcb2_p  = double_ptr()
    
    j       = lib.iauTdbtcb(tdb1, tdb2, tcb1_p, tcb2_p)
   
    tcb1    = tcb1_p[0]
    tcb2    = tcb2_p[0]
    
    return tcb1, tcb2 

#2023-08-04
#int iauTdbtt(double tdb1, double tdb2, double dtr,
#             double *tt1, double *tt2 )      
def pymTdbtt(tdb1, tdb2, dtr):  
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

    '''
    tt1_p  = double_ptr()
    tt2_p  = double_ptr()
    j = lib.iauTdbtt(tdb1, tdb2, dtr, tt1_p, tt2_p)
    
    tt1  = tt1_p[0]
    tt2  = tt2_p[0]
 
    return tt1, tt2


#int iauTf2d(char s, int ihour, int imin, double sec, double *days)
iau_tf2d_msg = {
                1:'ihour outside range 0-23',
                2:'imin outside range 0-59',
                3:'sec outside range 0-59.999...'
                }    
def pymTf2d(s, ihour, imin, sec): 
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
    
    Raises
    ------
    ValueError
        1:'ihour outside range 0-23',
        2:'imin outside range 0-59',
        3:'sec outside range 0-59.999...'

    Returns
    -------
    days : float
        interval in days

    '''
    days_p  = double_ptr()
     
    j = lib.iauTf2d(s, ihour, imin, sec, days_p)
    
    if j > 0:
        ws.warn(iau_tf2d_msg[j], UserWarning, 2)
    days = days_p[0]    
    return days 

#int iauTttai(double tt1, double tt2, double *tai1, double *tai2)
def pymTttai(tt1, tt2):
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

    '''
    
    tai1_p  = double_ptr()
    tai2_p  = double_ptr()
    j       = lib.iauTttai(tt1, tt2, tai1_p, tai2_p)
    
    tai1    = tai1_p[0] 
    tai2    = tai2_p[0]
    
    return tai1, tai2

#int iauTttcg(double tt1, double tt2, double *tcg1, double *tcg2)
def pymTttcg(tt1, tt2):  
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

    '''
    tcg1_p  = double_ptr()
    tcg2_p  = double_ptr()
    
    j       = lib.iauTttcg(tt1, tt2, tcg1_p, tcg2_p)
    
    tcg1    = tcg1_p[0]
    tcg2    = tcg2_p[0]
   
    return tcg1, tcg2 


#int iauTttdb(double tt1, double tt2, double dtr,
#             double *tdb1, double *tdb2)
def pymTttdb(tt1, tt2, dtr):  
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

    '''

    tdb1_p  = double_ptr()
    tdb2_p  = double_ptr()
    j       = lib.iauTttdb(tt1, tt2, dtr, tdb1_p, tdb2_p)
   
    tdb1    = tdb1_p[0]
    tdb2    = tdb2_p[0]
    
    return tdb1, tdb2  

#int iauTtut1(double tt1, double tt2, double dt,
#             double *ut11, double *ut12)         
def pymTtut1(tt1, tt2, dtr):  
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

    '''
    ut11_p  = double_ptr()
    ut12_p  = double_ptr()
    j       = lib.iauTtut1(tt1, tt2, dtr, ut11_p, ut12_p)
   
    ut11    = ut11_p[0]
    ut12    = ut12_p[0]
    return ut11, ut12


#int iauUt1tai(double ut11, double ut12, double dta,
#              double *tai1, double *tai2)   
#The argument dta, i.e. UT1-TAI, is an observed quantity, and is
#     available from IERS tabulations.   
def pymUt1tai(ut11, ut12, dta):  
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

    '''
    tai1_p  = double_ptr()
    tai2_p  = double_ptr()
    j       = lib.iauUt1tai(ut11, ut12, dta, tai1_p, tai2_p)
    tai1    = tai1_p[0]
    tai2    = tai2_p[0]
    return tai1, tai2 


#int iauUt1tt(double ut11, double ut12, double dt,
#             double *tt1, double *tt2)
def pymUt1tt(ut11, ut12, dt):  
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

    '''
  
    tt1_p  = double_ptr()
    tt2_p  = double_ptr()
    j      = lib.iauUt1tt(ut11, ut12, dt,  tt1_p, tt2_p)
   
    tt1    = tt1_p[0]
    tt2    = tt2_p[0]
    
    return tt1, tt2  


#int iauUt1utc(double ut11, double ut12, double dut1,
#              double *utc1, double *utc2)     
iau_ut1utc_msg = {
                  1: 'dubious year',
                 -1: 'unacceptable date',
                  }       
def pymUt1utc(ut11, ut12, dut1):  
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

    Raises
    ------
    ValueError
        1: 'dubious year',
       -1: 'unacceptable date',

    Returns
    -------
    utc1 : float
        UTC as a 2-part quasi Julian Date
    utc2 : float
        UTC as a 2-part quasi Julian Date

    '''
    utc1_p  = double_ptr()
    utc2_p  = double_ptr()
    
    j   = lib.iauUt1utc(ut11, ut12, dut1, utc1_p, utc2_p)
    if   j < 0:
        raise ValueError(iau_ut1utc_msg[j])
    elif j > 0:
        ws.warn(iau_ut1utc_msg[j], UserWarning, 2)
        
    utc1 = utc1_p [0]
    utc2 = utc2_p [0]
    
    return utc1, utc2   


#int iauUtctai(double utc1, double utc2, double *tai1, double *tai2)
iau_utctai_msg = {
                  1: 'dubious year',
                 -1: 'unacceptable date',
                  }    
def pymUtctai(utc1, utc2):  
    '''
    Time scale transformation:  Coordinated Universal Time, UTC, to
    International Atomic Time, TAI.

    Parameters
    ----------
    utc1 : float
        UTC as a 2-part quasi Julian Date
    utc2 : float
        UTC as a 2-part quasi Julian Date

    Raises
    ------
    ValueError
        1: 'dubious year',
       -1: 'unacceptable date',

    Returns
    -------
    tai1 : float
        TAI as a 2-part Julian Date
    tai2 : float
        TAI as a 2-part Julian Date

    '''
    tai1_p  = double_ptr()
    tai2_p  = double_ptr()
    
    j = lib.iauUtctai(utc1, utc2,  tai1_p, tai2_p)
    if   j < 0:
        raise ValueError(iau_utctai_msg[j])
    elif j > 0:
        ws.warn(iau_utctai_msg[j], UserWarning, 2)
    
    tai1    = tai1_p[0]
    tai2    = tai2_p[0] 
    
    return tai1, tai2  


#Hereafter use c_double() to replace double_ptr()
    
#int iauUtcut1(double utc1, double utc2, double dut1,
#              double *ut11, double *ut12)
iau_utcut1_msg = iau_ut1utc_msg
def pymUtcut1(utc1, utc2, dut1):  
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

    Raises
    ------
    ValueError
        1: 'dubious year',
       -1: 'unacceptable date',

    Returns
    -------
    ut11 : float
        UT1 as a 2-part Julian Date
    ut12 : float
        UT1 as a 2-part Julian Date

    '''
    ut11 = c_double()
    ut12 = c_double()
    
    j    = lib.iauUtcut1(utc1, utc2, dut1, ut11, ut12)
    
    if   j < 0:
        raise ValueError(iau_utcut1_msg[j])
    elif j > 0:
        ws.warn(iau_utcut1_msg[j], UserWarning, 2)
        
    #return ut11, ut12
    return ut11[0], ut12[0]   


#/* Astronomy/Astrometry */

#void iauAb(double pnat[3], double v[3], double s, double bm1,
#           double ppr[3]);
def pymAb(pnat, v, s, bm1):
    '''
    Apply aberration to transform natural direction into proper direction.

    Parameters
    ----------
    pnat : numpy.matrix(1,3)
        natural direction to the source (unit vector)
    v : numpy.matrix(1,3)
        observer barycentric velocity in units of c
    s : float
        distance between the Sun and the observer (au)
    bm1 : float
        sqrt(1-|v|^2): reciprocal of Lorenz factor

    Returns
    -------
    ppr : numpy.matrix(1,3)
        proper direction to source (unit vector)

    '''
    #pnat_vec = ffi.cast("double [3]",  ffi.from_buffer(pnat))
    #v_vec    = ffi.cast("double [3]",  ffi.from_buffer(v))
    
    pnat_vec       = vector_double(pnat)
    v_vec          = vector_double(v)
    ppr,  ppr_vec  = vector_double_ptr()
    
    lib.iauAb(pnat_vec, v_vec, s, bm1, ppr_vec)
    
    return ppr

#2023-08-05
#void iauAe2hd (double az, double el, double phi,
#               double *ha, double *dec)
def pymAe2hd(az, el, phi):
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
    ha  = c_double()
    dec = c_double()
    
    lib.iauAe2hd(az, el, phi, ha, dec) 
    
    return ha[0], dec[0] 



#void iauApcg(double date1, double date2,
#             double ebpv[2][3], double ehp[3],
#             iauASTROM *astrom)    

def pymApcg(date1, date2, ebpv,  ehp):
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
    ebpv : numpy.matrix(2,3)
        Earth barycentric pos/vel (au, au/day)     
    ehp : numpy.matrix(1,3)
        Earth heliocentric position (au)

    Returns
    -------
    astrom : pymASTROM class
        star-independent astrometry parameters    
        >pmt : float : PM time interval (SSB, Julian years)
        >eb : numpy.matrix(1,3) : SSB to observer (vector, au)
        >eh : numpy.matrix(1,3) : Sun to observer (unit vector)
        >em : float : distance from Sun to observer (au)
        >v : numpy.matrix(1,3) : barycentric observer velocity (vector, c)
        >bm1 : float : sqrt(1-|v|^2): reciprocal of Lorenz factor
        >bpn : numpy.matrix(3,3) : bias-precession-nutation matrix
        >along : float : unchanged
        >xpl : float : unchanged
        >ypl : float : unchanged
        >sphi : float : unchanged
        >cphi : float : unchanged
        >diurab : float : unchanged
        >eral : float : unchanged
        >refa : float : unchanged 
        >refb : float : unchanged
    '''
    astrom  = pymASTROM()
    
    ebpv_v  = vector_double2(ebpv)
    ehp_v   = vector_double (ehp)
   
    lib.iauApcg(date1, date2,  ebpv_v,  ehp_v,  astrom)
    
    return  astrom


#void iauApcg13(double date1, double date2, iauASTROM *astrom)
def pymApcg13(date1, date2):
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
    astrom : pymASTROM class
        star-independent astrometry parameters    
        >pmt : float : PM time interval (SSB, Julian years)
        >eb : numpy.matrix(1,3) : SSB to observer (vector, au)
        >eh : numpy.matrix(1,3) : Sun to observer (unit vector)
        >em : float : distance from Sun to observer (au)
        >v : numpy.matrix(1,3) : barycentric observer velocity (vector, c)
        >bm1 : float : sqrt(1-|v|^2): reciprocal of Lorenz factor
        >bpn : numpy.matrix(3,3) : bias-precession-nutation matrix
        >along : float : unchanged
        >xpl : float : unchanged
        >ypl : float : unchanged
        >sphi : float : unchanged
        >cphi : float : unchanged
        >diurab : float : unchanged
        >eral : float : unchanged
        >refa : float : unchanged 
        >refb : float : unchanged

    '''
    astrom  = pymASTROM()
   
    lib.iauApcg13(date1, date2, astrom)
    
    return  astrom


#void iauApci(double date1, double date2,
#             double ebpv[2][3], double ehp[3],
#             double x, double y, double s,
#             iauASTROM *astrom)
def pymApci(date1, date2, ebpv,  ehp,  x, y, s): 
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
    ebpv : numpy.matrix(2,3)
        Earth barycentric position/velocity (au, au/day)     
    ehp : numpy.matrix(1,3)
        Earth heliocentric position (au)    
    x : flaot
        CIP X,Y (components of unit vector)    
    y : float
        CIP X,Y (components of unit vector)    
    s : float
        the CIO locator s (radians)

    Returns
    -------
    astrom : pymASTROM class
        star-independent astrometry parameters    
        >pmt : float : PM time interval (SSB, Julian years)
        >eb : numpy.matrix(1,3) : SSB to observer (vector, au)
        >eh : numpy.matrix(1,3) : Sun to observer (unit vector)
        >em : float : distance from Sun to observer (au)
        >v : numpy.matrix(1,3) : barycentric observer velocity (vector, c)
        >bm1 : float : sqrt(1-|v|^2): reciprocal of Lorenz factor
        >bpn : numpy.matrix(3,3) : bias-precession-nutation matrix
        >along : float : unchanged
        >xpl : float : unchanged
        >ypl : float : unchanged
        >sphi : float : unchanged
        >cphi : float : unchanged
        >diurab : float : unchanged
        >eral : float : unchanged
        >refa : float : unchanged 
        >refb : float : unchanged

    '''
    ebpv_v  = vector_double2(ebpv)
    ehp_v   = vector_double (ehp)
    
    astrom  = pymASTROM()
    
    lib.iauApci(date1, date2, ebpv_v,  ehp_v,  x,  y,  s,  astrom)
    
    return  astrom

#void iauApci13(double date1, double date2,
#               iauASTROM *astrom, double *eo)    

def pymApci13(date1, date2): 
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
    astrom : pymASTROM class
        star-independent astrometry parameters     
    eo : float
        equation of the origins (ERA-GST)
        astrom:
        >pmt : float : PM time interval (SSB, Julian years)
        >eb : numpy.matrix(1,3) : SSB to observer (vector, au)
        >eh : numpy.matrix(1,3) : Sun to observer (unit vector)
        >em : float : distance from Sun to observer (au)
        >v : numpy.matrix(1,3) : barycentric observer velocity (vector, c)
        >bm1 : float : sqrt(1-|v|^2): reciprocal of Lorenz factor
        >bpn : numpy.matrix(3,3) : bias-precession-nutation matrix
        >along : float : unchanged
        >xpl : float : unchanged
        >ypl : float : unchanged
        >sphi : float : unchanged
        >cphi : float : unchanged
        >diurab : float : unchanged
        >eral : float : unchanged
        >refa : float : unchanged 
        >refb : float : unchanged
    
    '''
    astrom  = pymASTROM()
    eo      = c_double()
    lib.iauApci13(date1, date2, astrom, eo) 
    
    return  astrom, eo[0]


#void iauApco(double date1, double date2,
#             double ebpv[2][3], double ehp[3],
#             double x, double y, double s, double theta,
#             double elong, double phi, double hm,
#             double xp, double yp, double sp,
#             double refa, double refb,
#             iauASTROM *astrom)   
def pymApco(date1, date2, ebpv, ehp, x, y, s, theta, elong, phi, hm, 
            xp, yp, sp, refa, refb): 
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
    ebpv : numpy.matrix(2,3)
        Earth barycentric position/velocity (au, au/day)     
    ehp : numpy.matrix(1,3)
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
    astrom : pymASTROM class
        star-independent astrometry parameters    
        >pmt : float : PM time interval (SSB, Julian years)
        >eb : numpy.matrix(1,3) : SSB to observer (vector, au)
        >eh : numpy.matrix(1,3) : Sun to observer (unit vector)
        >em : float : distance from Sun to observer (au)
        >v : numpy.matrix(1,3) : barycentric observer velocity (vector, c)
        >bm1 : float : sqrt(1-|v|^2): reciprocal of Lorenz factor
        >bpn : numpy.matrix(3,3) : bias-precession-nutation matrix
        >along : float : adjusted longitude (radians)
        >xpl : float : polar motion xp wrt local meridian (radians)
        >ypl : float : polar motion yp wrt local meridian (radians)
        >sphi : float : sine of geodetic latitude
        >cphi : float : cosine of geodetic latitude
        >diurab : float : magnitude of diurnal aberration vector
        >eral : float : "local" Earth rotation angle (radians)
        >refa : float : refraction constant A (radians)
        >refb : float : refraction constant B (radians)

    '''
    ebpv_v  = vector_double2(ebpv)
    ehp_v   = vector_double (ehp)
    
    astrom  = pymASTROM()
    
    lib.iauApco(date1, date2, ebpv_v, ehp_v, x, y, s, theta, elong, phi, hm, 
            xp, yp, sp, refa, refb, astrom)
    
    return  astrom


#int iauApco13(double utc1, double utc2, double dut1,
#              double elong, double phi, double hm, double xp, double yp,
#              double phpa, double tc, double rh, double wl,
#              iauASTROM *astrom, double *eo)
def pymApco13(utc1, utc2, dut1, elong, phi, hm, xp, yp, 
              phpa, tc, rh,  wl): 
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
    astrom : pymASTROM class
        star-independent astrometry parameters     
    eo : float
        equation of the origins (ERA-GST)
        astrom:
        >pmt : float : PM time interval (SSB, Julian years)
        >eb : numpy.matrix(1,3) : SSB to observer (vector, au)
        >eh : numpy.matrix(1,3) : Sun to observer (unit vector)
        >em : float : distance from Sun to observer (au)
        >v : numpy.matrix(1,3) : barycentric observer velocity (vector, c)
        >bm1 : float : sqrt(1-|v|^2): reciprocal of Lorenz factor
        >bpn : numpy.matrix(3,3) : bias-precession-nutation matrix
        >along : float : adjusted longitude (radians)
        >xpl : float : polar motion xp wrt local meridian (radians)
        >ypl : float : polar motion yp wrt local meridian (radians)
        >sphi : float : sine of geodetic latitude
        >cphi : float : cosine of geodetic latitude
        >diurab : float : magnitude of diurnal aberration vector
        >eral : float : "local" Earth rotation angle (radians)
        >refa : float : refraction constant A (radians)
        >refb : float : refraction constant B (radians)

    '''
    astrom  = pymASTROM()
    eo      = c_double()
    
    j = lib.iauApco13(utc1, utc2, dut1, elong, phi, hm, xp, yp, 
              phpa, tc, rh,  wl, astrom, eo)
#debug message from previous procedures iauTaiutc   
    if   j < 0:
        raise ValueError(iau_taiutc_msg[j])
    elif j > 0:
        ws.warn(iau_taiutc_msg[j], UserWarning, 2)
    return  astrom, eo[0]

#void iauApcs(double date1, double date2, double pv[2][3],
#             double ebpv[2][3], double ehp[3],
#             iauASTROM *astrom)

def pymApcs(date1, date2, pv, ebpv, ehp): 
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
    pv : numpy.matrix(2,3)
        observer's geocentric pos/vel (m, m/s)    
    ebpv : numpy.matrix(2,3)
        Earth barycentric position/velocity (au, au/day)     
    ehp : numpy.matrix(1,3)
        Earth heliocentric position (au)    

    Returns
    -------
    astrom : pymASTROM class
        star-independent astrometry parameters     
        >pmt : float : PM time interval (SSB, Julian years)
        >eb : numpy.matrix(1,3) : SSB to observer (vector, au)
        >eh : numpy.matrix(1,3) : Sun to observer (unit vector)
        >em : float : distance from Sun to observer (au)
        >v : numpy.matrix(1,3) : barycentric observer velocity (vector, c)
        >bm1 : float : sqrt(1-|v|^2): reciprocal of Lorenz factor
        >bpn : numpy.matrix(3,3) : bias-precession-nutation matrix
        >along : float : unchanged
        >xpl : float : unchanged
        >ypl : float : unchanged
        >sphi : float : unchanged
        >cphi : float : unchanged
        >diurab : float : unchanged
        >eral : float : unchanged
        >refa : float : unchanged 
        >refb : float : unchanged

    '''
    pv_v    = vector_double2(pv)
    ebpv_v  = vector_double2(ebpv)
    ehp_v   = vector_double (ehp)
    
    astrom  = pymASTROM()
    
    lib.iauApcs(date1, date2, pv_v, ebpv_v, ehp_v, astrom)
    
    return  astrom  



#void iauApcs13(double date1, double date2, double pv[2][3],
#               iauASTROM *astrom)  
def pymApcs13(date1, date2, pv): 
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
    pv : numpy.matrix(2,3)
        observer's geocentric pos/vel (m, m/s)    

    Returns
    -------
    astrom : pymASTROM class
        star-independent astrometry parameters     
        >pmt : float : PM time interval (SSB, Julian years)
        >eb : numpy.matrix(1,3) : SSB to observer (vector, au)
        >eh : numpy.matrix(1,3) : Sun to observer (unit vector)
        >em : float : distance from Sun to observer (au)
        >v : numpy.matrix(1,3) : barycentric observer velocity (vector, c)
        >bm1 : float : sqrt(1-|v|^2): reciprocal of Lorenz factor
        >bpn : numpy.matrix(3,3) : bias-precession-nutation matrix
        >along : float : unchanged
        >xpl : float : unchanged
        >ypl : float : unchanged
        >sphi : float : unchanged
        >cphi : float : unchanged
        >diurab : float : unchanged
        >eral : float : unchanged
        >refa : float : unchanged 
        >refb : float : unchanged

    '''
    pv_v    = vector_double2(pv)
    astrom  = pymASTROM()
    
    lib.iauApcs13(date1, date2, pv_v, astrom)
    
    return  astrom 


#void iauAper(double theta, iauASTROM *astrom)  
def pymAper(theta, astrom):
    '''
    In the star-independent astrometry parameters, update only the
    Earth rotation angle, supplied by the caller explicitly.

    Parameters
    ----------
    theta : float
        Earth rotation angle (radians)    
    astrom : pymASTROM class
        star-independent astrometry parameters     
        >pmt : float : not used
        >eb : numpy.matrix(1,3) : not used
        >eh : numpy.matrix(1,3) : not used
        >em : float : not used
        >v : numpy.matrix(1,3) : not used
        >bm1 : float : not used
        >bpn : numpy.matrix(3,3) : not used
        >along : float : longitude + s' (radians)
        >xpl : float : not used
        >ypl : float : not used
        >sphi : float : not used
        >cphi : float : not used
        >diurab : float : not used
        >eral : float : not used
        >refa : float : not used
        >refb : float : not used

    Returns
    -------
    astrom : pymASTROM class
        star-independent astrometry parameters     
        >pmt : float : unchanged
        >eb : numpy.matrix(1,3) : unchanged
        >eh : numpy.matrix(1,3) : unchanged
        >em : float : unchanged
        >v : numpy.matrix(1,3) : unchanged
        >bm1 : float : unchanged
        >bpn : numpy.matrix(3,3) : unchanged
        >along : float : unchanged
        >xpl : float : unchanged
        >ypl : float : unchanged
        >sphi : float : unchanged
        >cphi : float : unchanged
        >diurab : float : unchanged
        >eral : float : "local" Earth rotation angle (radians)
        >refa : float : unchanged 
        >refb : float : unchanged

    '''
    lib.iauAper(theta, astrom)
    
    return   


#void iauAper13(double ut11, double ut12, iauASTROM *astrom)  
def pymAper13(ut11, ut12, astrom):
    '''
    In the star-independent astrometry parameters, update only the
    Earth rotation angle.  The caller provides UT1, (n.b. not UTC).

    Parameters
    ----------
    ut11 : float
        UT1 as a 2-part Julian Date    
    ut12 : float
        UT1 as a 2-part Julian Date    
    astrom : pymASTROM class
        star-independent astrometry parameters     
        >pmt : float : not used
        >eb : numpy.matrix(1,3) : not used
        >eh : numpy.matrix(1,3) : not used
        >em : float : not used
        >v : numpy.matrix(1,3) : not used
        >bm1 : float : not used
        >bpn : numpy.matrix(3,3) : not used
        >along : float : longitude + s' (radians)
        >xpl : float : not used
        >ypl : float : not used
        >sphi : float : not used
        >cphi : float : not used
        >diurab : float : not used
        >eral : float : not used
        >refa : float : not used
        >refb : float : not used

    Returns
    -------
    astrom : pymASTROM class
        star-independent astrometry parameters     
        >pmt : float : unchanged
        >eb : numpy.matrix(1,3) : unchanged
        >eh : numpy.matrix(1,3) : unchanged
        >em : float : unchanged
        >v : numpy.matrix(1,3) : unchanged
        >bm1 : float : unchanged
        >bpn : numpy.matrix(3,3) : unchanged
        >along : float : unchanged
        >xpl : float : unchanged
        >ypl : float : unchanged
        >sphi : float : unchanged
        >cphi : float : unchanged
        >diurab : float : unchanged
        >eral : float : "local" Earth rotation angle (radians)
        >refa : float : unchanged 
        >refb : float : unchanged

    '''
    lib.iauAper13(ut11, ut12, astrom)
    
    return  astrom


#void iauApio(double sp, double theta,
#             double elong, double phi, double hm, double xp, double yp,
#             double refa, double refb,
#             iauASTROM *astrom)    

def pymApio(sp, theta, elong, phi,  hm,  xp,  yp,  refa,  refb):
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
    astrom : pymASTROM class
        star-independent astrometry parameters     
        >pmt : float : unchanged
        >eb : numpy.matrix(1,3) : unchanged
        >eh : numpy.matrix(1,3) : unchanged
        >em : float : unchanged
        >v : numpy.matrix(1,3) : unchanged
        >bm1 : float : unchanged
        >bpn : numpy.matrix(3,3) : unchanged
        >along : float : adjusted longitude (radians)
        >xpl : float : polar motion xp wrt local meridian (radians)
        >ypl : float : polar motion yp wrt local meridian (radians)
        >sphi : float : sine of geodetic latitude
        >cphi : float : cosine of geodetic latitude
        >diurab : float : magnitude of diurnal aberration vector
        >eral : float : "local" Earth rotation angle (radians)
        >refa : float : refraction constant A (radians)
        >refb : float : refraction constant B (radians)

    '''
    astrom  = pymASTROM()
    
    lib.iauApio(sp, theta, elong, phi, hm, xp,  yp,  refa,  refb, astrom)
    
    return   astrom 



#int iauApio13(double utc1, double utc2, double dut1,
#              double elong, double phi, double hm, double xp, double yp,
#              double phpa, double tc, double rh, double wl,
#              iauASTROM *astrom)
    
def pymApio13(utc1, utc2, dut1, elong, phi, hm, xp, yp, 
              phpa, tc, rh,  wl): 
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
    astrom : pymASTROM class
        star-independent astrometry parameters     
        >pmt : float : unchanged
        >eb : numpy.matrix(1,3) : unchanged
        >eh : numpy.matrix(1,3) : unchanged
        >em : float : unchanged
        >v : numpy.matrix(1,3) : unchanged
        >bm1 : float : unchanged
        >bpn : numpy.matrix(3,3) : unchanged
        >along : float : adjusted longitude (radians)
        >xpl : float : polar motion xp wrt local meridian (radians)
        >ypl : float : polar motion yp wrt local meridian (radians)
        >sphi : float : sine of geodetic latitude
        >cphi : float : cosine of geodetic latitude
        >diurab : float : magnitude of diurnal aberration vector
        >eral : float : "local" Earth rotation angle (radians)
        >refa : float : refraction constant A (radians)
        >refb : float : refraction constant B (radians)

    '''
    astrom  = pymASTROM()
      
    j = lib.iauApio13(utc1, utc2, dut1, elong, phi, hm, xp, yp, 
                     phpa, tc, rh,  wl, astrom)
#debug message from previous procedures iauTaiutc   
   
    if   j < 0:
        raise ValueError(iau_taiutc_msg[j])
    elif j > 0:
        ws.warn(iau_taiutc_msg[j], UserWarning, 2)
    
    return  astrom


#void iauAtcc13(double rc, double dc,
#               double pr, double pd, double px, double rv,
#               double date1, double date2,
#               double *ra, double *da)

def pymAtcc13(rc,  dc,  pr,  pd,  px,  rv, date1, date2):
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
    ra  = c_double()
    da  = c_double()
    
    lib.iauAtcc13(rc,  dc,  pr,  pd,  px,  rv, date1, date2, 
                  ra, da)
    
    return  ra[0], da[0]



#void iauAtccq(double rc, double dc,
#              double pr, double pd, double px, double rv,
#              iauASTROM *astrom, double *ra, double *da) 
def pymAtccq(rc,  dc,  pr,  pd,  px,  rv,  astrom): 
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
    astrom : pymASTROM class
        star-independent astrometry parameters     
        >pmt : float : PM time interval (SSB, Julian years)
        >eb : numpy.matrix(1,3) : SSB to observer (vector, au)
        >eh : numpy.matrix(1,3) : Sun to observer (unit vector)
        >em : float : distance from Sun to observer (au)
        >v : numpy.matrix(1,3) : barycentric observer velocity (vector, c)
        >bm1 : float : sqrt(1-|v|^2): reciprocal of Lorenz factor
        >bpn : numpy.matrix(3,3) : bias-precession-nutation matrix
        >along : float : adjusted longitude (radians)
        >xpl : float : polar motion xp wrt local meridian (radians)
        >ypl : float : polar motion yp wrt local meridian (radians)
        >sphi : float : sine of geodetic latitude
        >cphi : float : cosine of geodetic latitude
        >diurab : float : magnitude of diurnal aberration vector
        >eral : float : "local" Earth rotation angle (radians)
        >refa : float : refraction constant A (radians)
        >refb : float : refraction constant B (radians)

    Returns
    -------
    ra : float
        ICRS astrometric RA (radians)    
    da : float
        ICRS astrometric Dec (radians)

    '''
    ra  = c_double()
    da  = c_double()
      
    lib.iauAtccq(rc,  dc,  pr,  pd,  px,  rv,  astrom, 
                     ra, da)
    return  ra[0], da[0] 


#void iauAtci13(double rc, double dc,
#               double pr, double pd, double px, double rv,
#               double date1, double date2,
#               double *ri, double *di, double *eo)   
def pymAtci13(rc,  dc,  pr,  pd,  px,  rv,  date1, date2):
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
    ri  = c_double()
    di  = c_double()
    eo  = c_double()
    
    lib.iauAtci13(rc,  dc,  pr,  pd,  px,  rv,  date1, date2, 
                  ri,  di,  eo)
    
    return  ri[0], di[0], eo[0]


#void iauAtciq(double rc, double dc,
#              double pr, double pd, double px, double rv,
#              iauASTROM *astrom, double *ri, double *di)
def pymAtciq(rc,  dc,  pr,  pd,  px,  rv,  astrom):
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
    astrom : pymASTROM class
        star-independent astrometry parameters     
        >pmt : float : PM time interval (SSB, Julian years)
        >eb : numpy.matrix(1,3) : SSB to observer (vector, au)
        >eh : numpy.matrix(1,3) : Sun to observer (unit vector)
        >em : float : distance from Sun to observer (au)
        >v : numpy.matrix(1,3) : barycentric observer velocity (vector, c)
        >bm1 : float : sqrt(1-|v|^2): reciprocal of Lorenz factor
        >bpn : numpy.matrix(3,3) : bias-precession-nutation matrix
        >along : float : adjusted longitude (radians)
        >xpl : float : polar motion xp wrt local meridian (radians)
        >ypl : float : polar motion yp wrt local meridian (radians)
        >sphi : float : sine of geodetic latitude
        >cphi : float : cosine of geodetic latitude
        >diurab : float : magnitude of diurnal aberration vector
        >eral : float : "local" Earth rotation angle (radians)
        >refa : float : refraction constant A (radians)
        >refb : float : refraction constant B (radians)
     
    Returns
    -------
    ri : float
        CIRS geocentric RA,Dec (radians)    
    di : float
        CIRS geocentric RA,Dec (radians)    
     
    '''
    
    ri  = c_double()
    di  = c_double()
     
    lib.iauAtciq(rc,  dc,  pr,  pd,  px,  rv,  astrom, 
                 ri,  di) 
                 
    return  ri[0], di[0]

 
#void iauAtciqn(double rc, double dc, double pr, double pd,
#               double px, double rv, iauASTROM *astrom,
#               int n, iauLDBODY b[], double *ri, double *di)

#to examine

def pymAtciqn(rc,  dc,  pr,  pd,  px,  rv,  astrom,  n,  b):
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
    astrom : pymASTROM class
        star-independent astrometry parameters     
        >pmt : float : PM time interval (SSB, Julian years)
        >eb : numpy.matrix(1,3) : SSB to observer (vector, au)
        >eh : numpy.matrix(1,3) : Sun to observer (unit vector)
        >em : float : distance from Sun to observer (au)
        >v : numpy.matrix(1,3) : barycentric observer velocity (vector, c)
        >bm1 : float : sqrt(1-|v|^2): reciprocal of Lorenz factor
        >bpn : numpy.matrix(3,3) : bias-precession-nutation matrix
        >along : float : adjusted longitude (radians)
        >xpl : float : polar motion xp wrt local meridian (radians)
        >ypl : float : polar motion yp wrt local meridian (radians)
        >sphi : float : sine of geodetic latitude
        >cphi : float : cosine of geodetic latitude
        >diurab : float : magnitude of diurnal aberration vector
        >eral : float : "local" Earth rotation angle (radians)
        >refa : float : refraction constant A (radians)
        >refb : float : refraction constant B (radians)
    n : int
        number of bodies    
    b : pymLDBODY class
        data for each of the n bodies    
        >bm : float : mass of the body (solar masses)    
        >dl : float : deflection limiter
        >pv : numpy.matrix(2,3) : barycentric PV of the body (au, au/day)

    Returns
    -------
    ri : float
        CIRS geocentric RA,Dec (radians)    
    di : float
        CIRS geocentric RA,Dec (radians)  

    '''
    ri     = c_double()
    di     = c_double()
    b_v    = pymLDBODY(b)
         
    lib.iauAtciqn(rc,  dc,  pr,  pd,  px,  rv,  astrom,  n,  b_v,
                  ri,  di) 
    
    return  ri[0], di[0]


#void iauAtciqz(double rc, double dc, iauASTROM *astrom,
#               double *ri, double *di)

def pymAtciqz(rc,  dc,  astrom): 
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
    astrom : pymASTROM class
        star-independent astrometry parameters     
        >pmt : float : PM time interval (SSB, Julian years)
        >eb : numpy.matrix(1,3) : SSB to observer (vector, au)
        >eh : numpy.matrix(1,3) : Sun to observer (unit vector)
        >em : float : distance from Sun to observer (au)
        >v : numpy.matrix(1,3) : barycentric observer velocity (vector, c)
        >bm1 : float : sqrt(1-|v|^2): reciprocal of Lorenz factor
        >bpn : numpy.matrix(3,3) : bias-precession-nutation matrix
        >along : float : adjusted longitude (radians)
        >xpl : float : polar motion xp wrt local meridian (radians)
        >ypl : float : polar motion yp wrt local meridian (radians)
        >sphi : float : sine of geodetic latitude
        >cphi : float : cosine of geodetic latitude
        >diurab : float : magnitude of diurnal aberration vector
        >eral : float : "local" Earth rotation angle (radians)
        >refa : float : refraction constant A (radians)
        >refb : float : refraction constant B (radians)

    Returns
    -------
    ri : float
        CIRS geocentric RA,Dec (radians)    
    di : float
        CIRS geocentric RA,Dec (radians)  

    '''
    ri  = c_double()
    di  = c_double()
      
    lib.iauAtciqz(rc,  dc,  astrom,  ri,  di) 
    
    return  ri[0], di[0]

#2023-08-06
#int iauAtco13(double rc, double dc,
#              double pr, double pd, double px, double rv,
#              double utc1, double utc2, double dut1,
#              double elong, double phi, double hm, double xp, double yp,
#              double phpa, double tc, double rh, double wl,
#              double *aob, double *zob, double *hob,
#              double *dob, double *rob, double *eo) 

def pymAtco13(rc,  dc,  pr,  pd,  px,  rv,
              utc1, utc2, dut1, elong, phi, hm, xp, yp, 
              phpa, tc, rh,  wl): 
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

    Raises
    ------
    ValueError
        1: 'dubious year',
       -1: 'unacceptable date',

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

    '''
    aob     = c_double()
    zob     = c_double()
    hob     = c_double()
    dob     = c_double()
    rob     = c_double()
    eo      = c_double()
   
    j = lib.iauAtco13(rc,  dc,  pr,  pd,  px,  rv,
                      utc1, utc2, dut1, elong, phi, hm, xp, yp, 
                      phpa, tc, rh,  wl,  
                      aob, zob, hob, 
                      dob, rob, eo)
#debug message from previous procedures iauTaiutc   
    if   j < 0:
        raise ValueError(iau_taiutc_msg[j])
    elif j > 0:
        ws.warn(iau_taiutc_msg[j], UserWarning, 2)
        
    return  aob[0], zob[0], hob[0], dob[0], rob[0], eo[0]   
    

#void iauAtic13(double ri, double di, double date1, double date2,
#               double *rc, double *dc, double *eo)  
def pymAtic13(ri,  di,  date1,  date2):
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
    rc  = c_double()
    dc  = c_double()
    eo  = c_double()
    lib.iauAtic13(ri,  di,  date1,  date2,  rc,  dc, eo)
    
    return  rc[0], dc[0], eo[0] 


#void iauAticq(double ri, double di, iauASTROM *astrom,
#              double *rc, double *dc)      
def pymAticq(ri,  di,  astrom): 
    '''
    Quick CIRS RA,Dec to ICRS astrometric place, given the star-
    independent astrometry parameters.

    Parameters
    ----------
    ri : float
        CIRS geocentric RA,Dec (radians)    
    di : float
        CIRS geocentric RA,Dec (radians)  
    astrom : pymASTROM class
        star-independent astrometry parameters     
        >pmt : float : PM time interval (SSB, Julian years)
        >eb : numpy.matrix(1,3) : SSB to observer (vector, au)
        >eh : numpy.matrix(1,3) : Sun to observer (unit vector)
        >em : float : distance from Sun to observer (au)
        >v : numpy.matrix(1,3) : barycentric observer velocity (vector, c)
        >bm1 : float : sqrt(1-|v|^2): reciprocal of Lorenz factor
        >bpn : numpy.matrix(3,3) : bias-precession-nutation matrix
        >along : float : adjusted longitude (radians)
        >xpl : float : polar motion xp wrt local meridian (radians)
        >ypl : float : polar motion yp wrt local meridian (radians)
        >sphi : float : sine of geodetic latitude
        >cphi : float : cosine of geodetic latitude
        >diurab : float : magnitude of diurnal aberration vector
        >eral : float : "local" Earth rotation angle (radians)
        >refa : float : refraction constant A (radians)
        >refb : float : refraction constant B (radians)

    Returns
    -------
    rc : float
        ICRS right ascension at J2000.0 (radians)    
    dc : float
        ICRS declination at J2000.0 (radians)    

    '''
    rc  = c_double()
    dc  = c_double()
      
    lib.iauAticq(ri,  di,  astrom,  rc,  dc) 
    return  rc[0], dc[0] 


'''    
void iauAticqn(double ri, double di, iauASTROM *astrom,
               int n, iauLDBODY b[], double *rc, double *dc)  

'''
def pymAticqn(ri,  di,  astrom,  n,  b):
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
    astrom : pymASTROM class
        star-independent astrometry parameters     
        >pmt : float : PM time interval (SSB, Julian years)
        >eb : numpy.matrix(1,3) : SSB to observer (vector, au)
        >eh : numpy.matrix(1,3) : Sun to observer (unit vector)
        >em : float : distance from Sun to observer (au)
        >v : numpy.matrix(1,3) : barycentric observer velocity (vector, c)
        >bm1 : float : sqrt(1-|v|^2): reciprocal of Lorenz factor
        >bpn : numpy.matrix(3,3) : bias-precession-nutation matrix
        >along : float : adjusted longitude (radians)
        >xpl : float : polar motion xp wrt local meridian (radians)
        >ypl : float : polar motion yp wrt local meridian (radians)
        >sphi : float : sine of geodetic latitude
        >cphi : float : cosine of geodetic latitude
        >diurab : float : magnitude of diurnal aberration vector
        >eral : float : "local" Earth rotation angle (radians)
        >refa : float : refraction constant A (radians)
        >refb : float : refraction constant B (radians)
    n : int
        number of bodies    
    b : pymLDBODY class
        data for each of the n bodies    
        >bm : float : mass of the body (solar masses)    
        >dl : float : deflection limiter
        >pv : numpy.matrix(2,3) : barycentric PV of the body (au, au/day)

    Returns
    -------
    rc : float
        ICRS right ascension at J2000.0 (radians)    
    dc : float
        ICRS declination at J2000.0 (radians)    

    '''
    
    b_v    = pymLDBODY(b)
    rc     = c_double()
    dc     = c_double()
     
    lib.iauAticqn(ri,  di,  astrom,  n,  b_v,
                  rc,  dc)
    
    return  rc[0], dc[0]  


'''
int iauAtio13(double ri, double di,
              double utc1, double utc2, double dut1,
              double elong, double phi, double hm, double xp, double yp,
              double phpa, double tc, double rh, double wl,
              double *aob, double *zob, double *hob,
              double *dob, double *rob)
'''
def pymAtio13(ri,  di, 
              utc1, utc2, dut1, 
              elong, phi, hm, xp, yp, 
              phpa, tc, rh,  wl): 
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

    Raises
    ------
    ValueError
        1: 'dubious year',
       -1: 'unacceptable date',

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
    aob     = c_double()
    zob     = c_double()
    hob     = c_double()
    dob     = c_double()
    rob     = c_double()
   
    j = lib.iauAtio13(ri,  di,   
                      utc1, utc2, dut1, 
                      elong, phi, hm, xp, yp, 
                      phpa, tc, rh,  wl,  
                      aob, zob, hob, 
                      dob, rob)
    
#debug message from previous procedures iauTaiutc   
    if   j < 0:
        raise ValueError(iau_taiutc_msg[j])
    elif j > 0:
        ws.warn(iau_taiutc_msg[j], UserWarning, 2)
        
    return  aob[0], zob[0], hob[0], dob[0], rob[0]  
    

'''
void iauAtioq(double ri, double di, iauASTROM *astrom,
              double *aob, double *zob,
              double *hob, double *dob, double *rob)
'''

def pymAtioq(ri,  di,  astrom):
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
    astrom : pymASTROM class
        star-independent astrometry parameters     
        >pmt : float : PM time interval (SSB, Julian years)
        >eb : numpy.matrix(1,3) : SSB to observer (vector, au)
        >eh : numpy.matrix(1,3) : Sun to observer (unit vector)
        >em : float : distance from Sun to observer (au)
        >v : numpy.matrix(1,3) : barycentric observer velocity (vector, c)
        >bm1 : float : sqrt(1-|v|^2): reciprocal of Lorenz factor
        >bpn : numpy.matrix(3,3) : bias-precession-nutation matrix
        >along : float : adjusted longitude (radians)
        >xpl : float : polar motion xp wrt local meridian (radians)
        >ypl : float : polar motion yp wrt local meridian (radians)
        >sphi : float : sine of geodetic latitude
        >cphi : float : cosine of geodetic latitude
        >diurab : float : magnitude of diurnal aberration vector
        >eral : float : "local" Earth rotation angle (radians)
        >refa : float : refraction constant A (radians)
        >refb : float : refraction constant B (radians)

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
    aob     = c_double()
    zob     = c_double()
    hob     = c_double()
    dob     = c_double()
    rob     = c_double()
   
    j = lib.iauAtioq( ri,  di,   astrom,
                      aob, zob, hob, 
                      dob, rob)   
    
    return  aob[0], zob[0], hob[0], dob[0], rob[0]


'''
int iauAtoc13(const char *type, double ob1, double ob2,
              double utc1, double utc2, double dut1,
              double elong, double phi, double hm, double xp, double yp,
              double phpa, double tc, double rh, double wl,
              double *rc, double *dc)
'''
def pymAtoc13(stype,   ob1,   ob2,
              utc1,  utc2,  dut1, 
              elong,  phi,    hm,  xp,  yp, 
              phpa,    tc,    rh,  wl): 
    '''
    Observed place at a groundbased site to to ICRS astrometric RA,Dec.
    The caller supplies UTC, site coordinates, ambient air conditions
    and observing wavelength.

    Parameters
    ----------
    stype : ctypes.c_char_p
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

    Raises
    ------
    ValueError
        1: 'dubious year',
       -1: 'unacceptable date',

    Returns
    -------
    rc : float
        ICRS right ascension at J2000.0 (radians)    
    dc : float
        ICRS declination at J2000.0 (radians)    

    '''
    rc  = c_double()
    dc  = c_double()
   
    j = lib.iauAtoc13(stype,   ob1,   ob2,
                      utc1,  utc2,   dut1, 
                      elong, phi, hm, xp, yp, 
                      phpa, tc, rh,  wl,  
                      rc,  dc)
#debug message from previous procedures iauTaiutc   
    if   j < 0:
        raise ValueError(iau_taiutc_msg[j])
    elif j > 0:
        ws.warn(iau_taiutc_msg[j], UserWarning, 2)
        
    return  rc[0], dc[0]  



'''
int iauAtoi13(const char *type, double ob1, double ob2,
              double utc1, double utc2, double dut1,
              double elong, double phi, double hm, double xp, double yp,
              double phpa, double tc, double rh, double wl,
              double *ri, double *di)
'''
def pymAtoi13(stype,   ob1,   ob2,
              utc1,  utc2,  dut1, 
              elong,  phi,    hm,  xp,  yp, 
              phpa,    tc,    rh,  wl): 
    '''
    Observed place to CIRS.  The caller supplies UTC, site coordinates,
    ambient air conditions and observing wavelength.

    Parameters
    ----------
    stype : ctypes.c_char_p
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
    
    Raises
    ------
    ValueError
        1: 'dubious year',
       -1: 'unacceptable date',
    
    Returns
    -------
    ri : float
        CIRS geocentric RA,Dec (radians)    
    di : float
        CIRS geocentric RA,Dec (radians)     

    '''
   
    ri  = c_double()
    di  = c_double()
   
    j = lib.iauAtoi13(stype,   ob1,   ob2,
                      utc1,  utc2,  dut1, 
                      elong, phi, hm, xp, yp, 
                      phpa, tc, rh,  wl,  
                      ri,  di)
#debug message from previous procedures iauTaiutc   
    if   j < 0:
        raise ValueError(iau_taiutc_msg[j])
    elif j > 0:
        ws.warn(iau_taiutc_msg[j], UserWarning, 2)
        
    return  ri[0], di[0] 


'''
void iauAtoiq(const char *type,
              double ob1, double ob2, iauASTROM *astrom,
              double *ri, double *di)
'''
def pymAtoiq(stype,   ob1,   ob2,  astrom):
    '''
    Quick observed place to CIRS, given the star-independent astrometry
    parameters.
    
    Use of this function is appropriate when efficiency is important and
    where many star positions are all to be transformed for one date.
    The star-independent astrometry parameters can be obtained by
    calling iauApio[13] or iauApco[13].

    Parameters
    ----------
    stype : ctypes.c_char_p
        type of coordinates - "R", "H" or "A"     
    ob1 : float
        observed Az, HA or RA (radians; Az is N=0,E=90)    
    ob2 : float
        observed ZD or Dec (radians)    
    astrom : pymASTROM class
        star-independent astrometry parameters     
        >pmt : float : PM time interval (SSB, Julian years)
        >eb : numpy.matrix(1,3) : SSB to observer (vector, au)
        >eh : numpy.matrix(1,3) : Sun to observer (unit vector)
        >em : float : distance from Sun to observer (au)
        >v : numpy.matrix(1,3) : barycentric observer velocity (vector, c)
        >bm1 : float : sqrt(1-|v|^2): reciprocal of Lorenz factor
        >bpn : numpy.matrix(3,3) : bias-precession-nutation matrix
        >along : float : adjusted longitude (radians)
        >xpl : float : polar motion xp wrt local meridian (radians)
        >ypl : float : polar motion yp wrt local meridian (radians)
        >sphi : float : sine of geodetic latitude
        >cphi : float : cosine of geodetic latitude
        >diurab : float : magnitude of diurnal aberration vector
        >eral : float : "local" Earth rotation angle (radians)
        >refa : float : refraction constant A (radians)
        >refb : float : refraction constant B (radians)

    Returns
    -------
    ri : float
        CIRS geocentric RA,Dec (radians)    
    di : float
        CIRS geocentric RA,Dec (radians)     

    '''
    ri  = c_double()
    di  = c_double()
   
    lib.iauAtoiq(stype,   ob1,   ob2,  astrom,
                 ri,  di)
        
    return  ri[0], di[0] 


'''
void iauHd2ae (double ha, double dec, double phi,
               double *az, double *el)
'''    
def pymHd2ae(ha, dec, phi):
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
    az  = c_double()
    el  = c_double()
   
    lib.iauHd2ae(ha, dec, phi, 
                 az, el) 
        
    return  az[0], el[0]  


'''
double iauHd2pa (double ha, double dec, double phi)

'''
def pymHd2pa(ha, dec, phi):
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
    return  lib.iauHd2pa(ha, dec, phi)



'''
void iauLd(double bm, double p[3], double q[3], double e[3],
           double em, double dlim, double p1[3])
'''
def pymLd(bm, p, q, e, em, dlim):
    '''
    Apply light deflection by a solar-system body, as part of
    transforming coordinate direction into natural direction.

    Parameters
    ----------
    bm : float
        mass of the gravitating body (solar masses)    
    p : numpy.matrix(1,3)
        direction from observer to source (unit vector)    
    q : numpy.matrix(1,3)
        direction from body to source (unit vector)    
    e : numpy.matrix(1,3)
        direction from body to observer (unit vector)    
    em : flaot
        distance from body to observer (au)    
    dlim : flaot
        deflection limiter

    Returns
    -------
    p1 : numpy.matrix(1,3)
        observer to deflected source (unit vector)

    '''
    p_v   = vector_double (p)
    q_v   = vector_double (q)
    e_v   = vector_double (e)
    
    p1, p1_v = vector_double_ptr()
    
    lib.iauLd(bm, p_v, q_v, e_v, em, dlim, p1_v) 

    return  p1



'''
void iauLdsun(double p[3], double e[3], double em, double p1[3])
'''
def pymLdsun(p, e, em):
    '''
    Deflection of starlight by the Sun.

    Parameters
    ----------
    p : numpy.matrix(1,3)
        direction from observer to star (unit vector)    
    e : numpy.matrix(1,3)
        direction from Sun to observer (unit vector)    
    em : float
        distance from Sun to observer (au)

    Returns
    -------
    p1 : numpy.matrix(1,3)
        observer to deflected star (unit vector)

    '''
    p_v   = vector_double (p)
    e_v   = vector_double (e)
    
    p1, p1_v = vector_double_ptr()
    
    lib.iauLdsun(p_v, e_v, em, p1_v) 

    return  p1



'''
void iauLdn(int n, iauLDBODY b[], double ob[3], double sc[3],
            double sn[3])
'''

def pymLdn(n, b, ob, sc):
    '''
    For a star, apply light deflection by multiple solar-system bodies,
    as part of transforming coordinate direction into natural direction.

    Parameters
    ----------
    n : int
        number of bodies    
    b : pymLDBODY class
        data for each of the n bodies    
        >bm : float : mass of the body (solar masses)    
        >dl : float : deflection limiter
        >pv : numpy.matrix(2,3) : barycentric PV of the body (au, au/day)
    ob : numpy.matrix(1,3)
        barycentric position of the observer (au)    
    sc : numpy.matrix(1,3)
        observer to star coord direction (unit vector)

    Returns
    -------
    sn : numpy.matrix(1,3)
        observer to deflected star (unit vector)

    '''                                   
    b_v    = pymLDBODY(b)
    ob_v   = vector_double (ob)
    sc_v   = vector_double (sc)
    
    sn, sn_v = vector_double_ptr()
    
    lib.iauLdn(n, b_v, ob_v, sc_v, sn_v) 

    return  sn



'''
void iauPmpx(double rc, double dc, double pr, double pd,
             double px, double rv, double pmt, double pob[3],
             double pco[3])
'''
def pymPmpx(rc, dc, pr, pd, px, rv, pmt, pob):
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
    pob : numpy.matrix(1,3)
        SSB to observer vector (au)    

    Returns
    -------
    pco : numpy.matrix(1,3)
        coordinate direction (BCRS unit vector)

    '''
   
    pob_v      = vector_double (pob)
    
    pco, pco_v = vector_double_ptr()
    
    lib.iauPmpx(rc, dc, pr, pd, px, rv, pmt, pob_v, pco_v) 

    return  pco



'''
int iauPmsafe(double ra1, double dec1, double pmr1, double pmd1,
              double px1, double rv1,
              double ep1a, double ep1b, double ep2a, double ep2b,
              double *ra2, double *dec2, double *pmr2, double *pmd2,
              double *px2, double *rv2)
'''

iau_pmsafe_msg = {
                  4: 'solution not converge',
                  2: 'excessive velocity',
                  1: 'distance overridden',
                  0: 'no warnings or errors',
                 -1: 'system error',
                  }  
 
def pymPmsafe(ra1,  dec1, pmr1, pmd1,
              px1,  rv1,
              ep1a, ep1b, ep2a, ep2b):
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

    '''
 
    ra2   = c_double()
    dec2  = c_double()
    pmr2  = c_double()
    pmd2  = c_double()
    
    px2   = c_double()
    rv2   = c_double()
   
    j = lib.iauPmsafe(ra1,  dec1, pmr1, pmd1,
                      px1,  rv1,
                      ep1a, ep1b,  ep2a, ep2b, 
                      (ra2),  (dec2),  (pmr2),  (pmd2),
                      (px2),  (rv2))

#debug message should check 
    if   j < 0:
         raise ValueError(iau_pmsafe_msg[j])
    elif j > 0:
         ws.warn(iau_pmsafe_msg[j], UserWarning, 2)
        
    return ra2[0], dec2[0], pmr2[0], pmd2[0], px2[0], rv2[0]  




'''
void iauPvtob(double elong, double phi, double hm,
              double xp, double yp, double sp, double theta,
              double pv[2][3])
'''
def pymPvtob(elong,  phi,  hm,  xp,  yp,  sp,  theta):
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
    pv : numpy.matrix(2,3)
        position/velocity vector (m, m/s, CIRS)

    '''
    pv, pv_v = vector_double2_ptr()
    
    lib.iauPvtob(elong,  phi,  hm,
                 xp,  yp,  sp,  theta,  pv_v) 

    return  pv



'''
void iauRefco(double phpa, double tc, double rh, double wl,
              double *refa, double *refb) 
''' 
def pymRefco(phpa,   tc,   rh,   wl):
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
    refa   = c_double()
    refb   = c_double()
    
    lib.iauRefco(phpa,   tc,   rh,   wl,
                (refa),   (refb)) 

    return  refa[0], refb[0]


#/* Astronomy/Gnomonic */
'''
int iauTpors(double xi, double eta, double a, double b,
             double *a01, double *b01, double *a02, double *b02)
'''
def pymTpors(xi,   eta,   a,   b):
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
    a01   = c_double()
    b01   = c_double()
    a02   = c_double()
    b02   = c_double()
    
    lib.iauTpors(xi,   eta,   a,   b,
                (a01),  (b01), (a02), (b02)) 

    return  a01[0], b01[0], a02[0], b02[0] 


'''
int iauTporv(double xi, double eta, double v[3],
             double v01[3], double v02[3])
'''

iau_tporv_msg = {
                  0: 'no solutions returned',
                  1: 'only the first solution is useful',
                  2: 'both solutions are useful',
                 }  
def pymTporv(xi,   eta,   v):
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
    v : numpy.matrix(3,3)
        star's direction cosines

    Raises
    ------
    ValueError
        0: 'no solutions returned',
        1: 'only the first solution is useful',
        2: 'both solutions are useful',

    Returns
    -------
    v01 : numpy.matrix(3,3)
        tangent point's direction cosines, Solution 1    
    v02 : numpy.matrix(3,3)
        tangent point's direction cosines, Solution 2    

    '''
    v_v      = vector_double (v)
    
    v01, v01_v = vector_double_ptr()
    v02, v02_v = vector_double_ptr()
    
    j = lib.iauTporv(xi,   eta,   v_v,  v01_v,  v02_v)

    return  v01, v02



'''
void iauTpsts(double xi, double eta, double a0, double b0,
              double *a, double *b)
'''
def pymTpsts(xi,   eta,   a0,   b0):
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
    a   = c_double()
    b   = c_double()
    
    lib.iauTpsts(xi,   eta,   a0,   b0,
                (a),  (b)) 

    return  a[0], b[0]



'''
void iauTpstv(double xi, double eta, double v0[3], double v[3])
'''
def pymTpstv(xi,   eta,   v0):
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
    v0 : numpy.matrix(3,3)
        tangent point's direction cosines

    Returns
    -------
    v : numpy.matrix(3,3)
        star's direction cosines

    '''
 
    v0_v     = vector_double (v0)
    
    v, v_v   = vector_double_ptr()
    
    lib.iauTpstv(xi,   eta,   v0_v,   v_v)
     
    return  v



'''
int iauTpxes(double a, double b, double a0, double b0,
             double *xi, double *eta)
'''


iau_tpxes_msg = {
                  0: 'OK',
                  1: 'star too far from axis',
                  2: 'antistar on tangent plane',
                  3: 'antistar too far from axis',
                 }  

def pymTpxes(a,   b,   a0,   b0):
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

    Raises
    ------
    ValueError
        0: 'OK',
        1: 'star too far from axis',
        2: 'antistar on tangent plane',
        3: 'antistar too far from axis',

    Returns
    -------
    xi : float
        rectangular coordinates of star image    
    eta : float
        rectangular coordinates of star image

    '''
    
    xi   = c_double()
    eta  = c_double()
    
    j = lib.iauTpxes(a,   b,   a0,   b0, 
                    (xi),   (eta))
    
    if  j > 0:
        ws.warn(iau_tpxes_msg[j], UserWarning, 2) 
        
    return  xi[0], eta[0]


'''
int iauTpxev(double v[3], double v0[3], double *xi, double *eta)
'''
def pymTpxev(v,   v0):
    '''
    In the tangent plane projection, given celestial direction cosines
    for a star and the tangent point, solve for the star's rectangular
    coordinates in the tangent plane.

    Parameters
    ----------
    v : numpy.matrix(3,3)
        direction cosines of star    
    v0 : numpy.matrix(3,3)
        direction cosines of tangent point

    Returns
    -------
    xi : float
        tangent plane coordinates of star    
    eta : float
        tangent plane coordinates of star

    '''
    v_v      = vector_double(v)
    v0_v     = vector_double(v0)
    
    xi   = c_double()
    eta  = c_double()
    
    j = lib.iauTpxev(v_v,  v0_v,   (xi),   (eta))
    
    if  j > 0:
        ws.warn(iau_tpxes_msg[j], UserWarning, 2)   
     
    return  xi[0], eta[0]


#2023-08-06   Sofa Tools for Earth Attitude 
###sofa_pn_f.pdf###
'''    
double iauAnp(double a)
'''
def pymAnp(a):
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
    return lib.iauAnp(a)


'''
void iauBi00(double *dpsibi, double *depsbi, double *dra)
'''
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

    dpsibi   = c_double()
    depsbi   = c_double()
    dra      = c_double()
    
    lib.iauBi00((dpsibi),  (depsbi),  (dra))
    
    return  dpsibi[0], depsbi[0], dra[0]


'''
void iauBpn2xy(double rbpn[3][3], double *x, double *y)
'''
def pymBpn2xy(rbpn):
    '''
    Extract from the bias-precession-nutation matrix the X,Y coordinates
    of the Celestial Intermediate Pole.

    Parameters
    ----------
    rbpn : numpy.matrix(3,3)
        celestial-to-true matrix

    Returns
    -------
    x : float
        Celestial Intermediate Pole    
    y : float
        Celestial Intermediate Pole

    '''
    rbpn_v = vector_double3(rbpn)
    
    x   = c_double()
    y   = c_double()
    
    lib.iauBpn2xy(rbpn_v,  (x),  (y))
    
    return  x[0], y[0]  



'''
void iauC2i00a(double date1, double date2, double rc2i[3][3])
'''
def pymC2i00a(date1, date2):
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
    rc2i : numpy.matrix(3,3)
        celestial-to-intermediate matrix

    '''
    rc2i, rc2i_v = vector_double3_ptr()
    lib.iauC2i00a(date1, date2, rc2i_v)
    
    return  rc2i 


'''
void iauC2i00b(double date1, double date2, double rc2i[3][3])
'''    
def pymC2i00b(date1, date2):
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
    rc2i : numpy.matrix(3,3)
        celestial-to-intermediate matrix

    '''
    rc2i, rc2i_v = vector_double3_ptr()
    lib.iauC2i00b(date1, date2, rc2i_v)
    
    return  rc2i 


'''
void iauC2ibpn(double date1, double date2, double rbpn[3][3],
               double rc2i[3][3])
'''
def pymC2ibpn(date1, date2, rbpn):
    '''
    Form the celestial-to-intermediate matrix for a given date given
    the bias-precession-nutation matrix.  IAU 2000.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date     
    date2 : float
        TT as a 2-part Julian Date    
    rbpn : numpy.matrix(3,3)
        celestial-to-true matrix

    Returns
    -------
    rc2i : numpy.matrix(3,3)
        celestial-to-intermediate matrix

    '''
    rbpn_v       = vector_double3(rbpn) 
    
    rc2i, rc2i_v = vector_double3_ptr()
     
    lib.iauC2ibpn(date1, date2, rbpn_v, rc2i_v)
    
    return  rc2i 


'''
void iauC2t00a(double tta, double ttb, double uta, double utb,
               double xp, double yp, double rc2t[3][3])
'''
def pymC2t00a(tta,  ttb,  uta,  utb,  xp,  yp):
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
    rc2t : numpy.matrix(3,3)
        celestial-to-terrestrial matrix

    '''
    rc2t, rc2t_v = vector_double3_ptr()
    
    lib.iauC2t00a(tta,  ttb,  uta,  utb,  xp,  yp, rc2t_v)
    
    return  rc2t


'''
void iauC2t00b(double tta, double ttb, double uta, double utb,
               double xp, double yp, double rc2t[3][3])
'''
def pymC2t00b(tta,  ttb,  uta,  utb,  xp,  yp):
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
    rc2t : numpy.matrix(3,3)
        celestial-to-terrestrial matrix
        
    '''
    rc2t, rc2t_v = vector_double3_ptr()
    
    lib.iauC2t00b(tta,  ttb,  uta,  utb,  xp,  yp, rc2t_v)
    
    return  rc2t



'''
void iauC2teqx(double rbpn[3][3], double gst, double rpom[3][3],
               double rc2t[3][3])
'''
def pymC2teqx(rbpn,  gst,  rpom):
    '''
    Assemble the celestial to terrestrial matrix from equinox-based
    components (the celestial-to-true matrix, the Greenwich Apparent
    Sidereal Time and the polar motion matrix).

    Parameters
    ----------
    rbpn : numpy.matrix(3,3)
        celestial-to-true matrix
    gst : float
        Greenwich (apparent) Sidereal Time (radians)
    rpom : numpy.matrix(3,3)
        polar-motion matrix

    Returns
    -------
    rc2t : numpy.matrix(3,3)
        celestial-to-terrestrial matrix

    '''
    rbpn_v       = vector_double3(rbpn) 
    rpom_v       = vector_double3(rpom) 
    
    rc2t, rc2t_v = vector_double3_ptr()
    
    lib.iauC2teqx(rbpn_v,  gst,  rpom_v,  rc2t_v)
    
    return  rc2t



'''
void iauC2tpe(double tta, double ttb, double uta, double utb,
              double dpsi, double deps, double xp, double yp,
              double rc2t[3][3])
'''
def pymC2tpe(tta,  ttb,  uta,  utb,  dpsi,  deps,  xp,  yp):
    '''
    Form the celestial to terrestrial matrix given the date, the UT1,
    the nutation and the polar motion.  IAU 2000.

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
    dpsi : float
        nutation     
    deps : float
        nutation         
    xp : float
        coordinates of the pole    
    yp : float
        coordinates of the pole  

    Returns
    -------
    rc2t : numpy.matrix(3,3)
        celestial-to-terrestrial matrix

    '''
    rc2t, rc2t_v = vector_double3_ptr()
    
    lib.iauC2tpe(tta,  ttb,  uta,  utb, dpsi,  deps,  xp,  yp, rc2t_v)
    
    return  rc2t


'''
void iauC2txy(double tta, double ttb, double uta, double utb,
              double x, double y, double xp, double yp,
              double rc2t[3][3])
'''
def pymC2txy(tta,  ttb,  uta,  utb,  x,  y,  xp,  yp):
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
    rc2t : numpy.matrix(3,3)
        celestial-to-terrestrial matrix

    '''
    rc2t, rc2t_v = vector_double3_ptr()
    
    lib.iauC2txy(tta,  ttb,  uta,  utb, x,  y,  xp,  yp, rc2t_v)
    
    return  rc2t



'''
double iauEe00(double date1, double date2, double epsa, double dpsi)
'''
def pymEe00(date1,  date2,   epsa,   dpsi):
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
    return   lib.iauEe00(date1,  date2,   epsa,  dpsi)


'''
double iauEe00a(double date1, double date2)
''' 
def pymEe00a(date1,  date2):
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
    return   lib.iauEe00a(date1,  date2)


'''
double iauEe00b(double date1, double date2)
''' 
def pymEe00b(date1,  date2):
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
    return   lib.iauEe00b(date1,  date2)


'''
double iauEe06a(double date1, double date2)
''' 
def pymEe06a(date1,  date2):
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
    return   lib.iauEe06a(date1,  date2)


'''
double iauEect00(double date1, double date2)
'''
def pymEect00(date1,  date2):
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
    return   lib.iauEect00(date1,  date2)


'''
double iauEo06a(double date1, double date2)
'''
def pymEo06a(date1,  date2):
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
    return   lib.iauEo06a(date1,  date2)


'''
double iauEors(double rnpb[3][3], double s)
'''
def pymEors(rnpb,  s):
    '''
    Equation of the origins, given the classical NPB matrix and the
    quantity s.

    Parameters
    ----------
    rnpb : numpy.matrix(3,3) 
        classical nutation x precession x bias matrix    
    s : float
        the quantity s (the CIO locator) in radians

    Returns
    -------
    function value : float
        the equation of the origins in radians

    '''
    rnpb_v       = vector_double3(rnpb) 
    
    return   lib.iauEors(rnpb_v,  s)


'''
double iauEqeq94(double date1, double date2)
'''
def pymEqeq94(date1,  date2):
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
    return   lib.iauEqeq94(date1,  date2)


'''
double iauEra00(double dj1, double dj2)
'''
def pymEra00(dj1,  dj2):
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
    return   lib.iauEra00(dj1,  dj2)


#/* Astronomy/FundamentalArgs */
'''
double iauFad03(double t)
'''
def pymFad03(t):
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
    return   lib.iauFad03(t)


'''
double iauFae03(double t)
'''
def pymFae03(t):
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
    return   lib.iauFae03(t)


'''
double iauFaf03(double t)
'''
def pymFaf03(t):
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
    return   lib.iauFaf03(t)

'''
double iauFaju03(double t) #Jupiter,  mean longitude
'''

def pymFaju03(t):
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
    return   lib.iauFaju03(t)

'''
double iauFal03(double t)
'''
def pymFal03(t):
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
    return   lib.iauFal03(t)

'''
double iauFalp03(double t)
'''
def pymFalp03(t):
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
    return   lib.iauFalp03(t)


'''
double iauFama03(double t)  #Mars,  mean longitude
'''
def pymFama03(t):
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
    return   lib.iauFama03(t)


'''
double iauFame03(double t)   #Mercucy,  mean longitude
'''
def pymFame03(t):
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
    return   lib.iauFame03(t)

'''
double iauFane03(double t)  #Neptune,  mean longitude
'''
def pymFane03(t):
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
    return   lib.iauFane03(t)

'''
double iauFaom03(double t)  #mean long. of Moon’s asc. node,
'''
def pymFaom03(t):
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
    return   lib.iauFaom03(t)

'''
double iauFapa03(double t)  # 
'''
def pymFapa03(t):
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
    return   lib.iauFapa03(t)

'''
double iauFasa03(double t)  #Saturn,  mean longitude 
'''
def pymFasa03(t):
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
    return   lib.iauFasa03(t)


'''
double iauFaur03(double t)  #Venus,  mean longitude  
'''
def pymFaur03(t):
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
    return   lib.iauFaur03(t)

'''
double iauFave03(double t)  # 
'''
def pymFave03(t):
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
    return   lib.iauFave03(t)


'''
void iauFw2m(double gamb, double phib, double psi, double eps,
             double r[3][3])
'''

def pymFw2m(gamb,  phib,  psi,  eps):
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
    r : numpy.matrix(3,3)
        rotation matrix

    '''
    r, r_v = vector_double3_ptr()
    
    lib.iauFw2m(gamb,  phib,  psi,  eps,  r_v)
    
    return  r 


'''
void iauFw2xy(double gamb, double phib, double psi, double eps,
              double *x, double *y)
'''    
def pymFw2xy(gamb,   phib,  psi,  eps):
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
    x    = c_double()
    y    = c_double()
    
    lib.iauFw2xy(gamb,   phib,  psi,  eps, 
                 (x),   (y))
    
    return  x[0], y[0] 


'''
double iauGmst00(double uta, double utb, double tta, double ttb)
'''
def pymGmst00(uta,  utb,  tta,  ttb):
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
    return  lib.iauGmst00(uta,  utb,  tta,  ttb)


'''
double iauGmst06(double uta, double utb, double tta, double ttb)
'''
def pymGmst06(uta,  utb,  tta,  ttb):
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
    return  lib.iauGmst06(uta,  utb,  tta,  ttb)

'''
double iauGmst82(double dj1, double dj2)
'''
def pymGmst82(dj1,  dj2):
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
    return  lib.iauGmst82(dj1,  dj2)


'''
double iauGst00a(double uta, double utb, double tta, double ttb)
'''
def pymGst00a(uta,  utb,  tta,  ttb):
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
    function value : float
        Greenwich apparent sidereal time (radians)

    '''
    return  lib.iauGst00a(uta,  utb,  tta,  ttb)

'''
double iauGst00b(double uta, double utb)
'''
def pymGst00b(uta,  utb):
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
    function value : float
        Greenwich apparent sidereal time (radians)

    '''
    return  lib.iauGst00b(uta,  utb)


'''
double iauGst06(double uta, double utb, double tta, double ttb,
                double rnpb[3][3])
'''
def pymGst06(uta,  utb,  tta,  ttb,  rnpb):
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
 
    rnpb_v       = vector_double3(rnpb) 
    
    return  lib.iauGst06(uta,  utb,  tta,  ttb,  rnpb_v)

'''
double iauGst06a(double uta, double utb, double tta, double ttb)
'''
def pymGst06a(uta,  utb, tta,  ttb):
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
    return  lib.iauGst06a(uta,  utb, tta,  ttb)


'''
double iauGst94(double uta, double utb)
'''
def pymGst94(dj1,  dj2):
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
    function value : float
        Greenwich apparent sidereal time (radians)

    '''
    return  lib.iauGst94(dj1,  dj2)



'''
void iauNum00a(double date1, double date2, double rmatn[3][3])
'''
def pymNum00a(date1, date2):
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
    rmatn : numpy.matrix(3,3)
        nutation matrix

    '''
    rmatn, rmatn_v = vector_double3_ptr()
    
    lib.iauNum00a(date1, date2, rmatn_v)
    
    return  rmatn


'''
void iauNum00b(double date1, double date2, double rmatn[3][3])
'''
def pymNum00b(date1, date2):
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
    rmatn : numpy.matrix(3,3)
        nutation matrix

    '''
    rmatn, rmatn_v = vector_double3_ptr()
    
    lib.iauNum00b(date1, date2, rmatn_v)
    
    return  rmatn


'''
void iauNum06a(double date1, double date2, double rmatn[3][3])
'''
def pymNum06a(date1, date2):
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
    rmatn : numpy.matrix(3,3)
        nutation matrix

    '''
    rmatn, rmatn_v = vector_double3_ptr()
    
    lib.iauNum06a(date1, date2, rmatn_v)
    
    return  rmatn


'''
void iauNumat(double epsa, double dpsi, double deps, double rmatn[3][3])
'''
def pymNumat(epsa,   dpsi,   deps):
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
    rmatn : numpy.matrix(3,3)
        nutation matrix

    '''
    rmatn, rmatn_v = vector_double3_ptr()
    
    lib.iauNumat(epsa,   dpsi,   deps,  rmatn_v)
    
    return  rmatn


'''
void iauNut00a(double date1, double date2, double *dpsi, double *deps)
'''
def pymNut00a(date1,  date2):
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
    dpsi = c_double()
    deps = c_double()
    
    lib.iauNut00a(date1,  date2,  (dpsi),   (deps))
    
    return  dpsi[0], deps[0] 


'''
void iauNut00b(double date1, double date2, double *dpsi, double *deps)
'''
def pymNut00b(date1,  date2):
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
    dpsi = c_double()
    deps = c_double()
    
    lib.iauNut00b(date1,  date2,  (dpsi),   (deps))
    
    return  dpsi[0], deps[0] 


'''
void iauNut06a(double date1, double date2, double *dpsi, double *deps)
'''
def pymNut06a(date1,  date2):
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
    dpsi = c_double()
    deps = c_double()
    
    lib.iauNut06a(date1,  date2,  (dpsi),   (deps))
    
    return  dpsi[0], deps[0] 



'''
void iauNut80(double date1, double date2, double *dpsi, double *deps)
'''
def pymNut80(date1,  date2):
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
    dpsi = c_double()
    deps = c_double()
    
    lib.iauNut80(date1,  date2,  (dpsi),   (deps))
    
    return  dpsi[0], deps[0] 

'''
void iauNutm80(double date1, double date2, double rmatn[3][3])
'''
def pymNutm80(date1, date2):
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
    rmatn : numpy.matrix(3,3)
        nutation matrix

    '''
    rmatn, rmatn_v = vector_double3_ptr() 
    
    lib.iauNutm80(date1, date2, rmatn_v)
    
    return  rmatn


'''
double iauObl06(double date1, double date2)
'''    
def pymObl06(date1, date2):
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
    return  lib.iauObl06(date1, date2)

'''
double iauObl80(double date1, double date2)
'''
def pymObl80(date1, date2):
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
    return  lib.iauObl80(date1, date2)



'''
void iauP06e(double date1, double date2,
             double *eps0, double *psia, double *oma, double *bpa,
             double *bqa,  double *pia, double *bpia,
             double *epsa, double *chia, double *za, double *zetaa,
             double *thetaa, double *pa,
             double *gam, double *phi, double *psi)
'''
def pymP06e(date1, date2):
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
    eps0  = c_double()
    psia  = c_double()
    oma   = c_double()
    bpa   = c_double()
  
    bqa   = c_double() 
    pia   = c_double() 
    bpia  = c_double()  
       
    epsa  = c_double()
    chia  = c_double()
    za    = c_double()
    zetaa = c_double()
    
    thetaa = c_double() 
    pa     = c_double()
    
    gam   = c_double() 
    phi   = c_double() 
    psi   = c_double()  
    
    lib.iauP06e(date1, date2, 
                (eps0), (psia), (oma), (bpa),
                (bqa),  (pia),  (bpia),
                (epsa), (chia), (za), (zetaa),
                (thetaa), (pa),
                (gam),  (phi),  (psi))
    
    return   eps0[0],  psia[0],   oma[0],  bpa[0],  \
             bqa[0],   pia[0],  bpia[0],            \
             epsa[0],  chia[0],    za[0], zetaa[0], \
             thetaa[0],   pa[0],                    \
             gam[0],   phi[0],  psi[0]
             
             
'''
void iauPb06(double date1, double date2,
             double *bzeta, double *bz, double *btheta)
'''             
def pymPb06(date1,  date2):
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
    bzeta  = c_double()
    bz     = c_double()
    btheta = c_double()
    
    lib.iauPb06(date1,  date2,  (bzeta),  (bz),  (btheta))
    
    return  bzeta[0],  bz[0],  btheta[0]             
             


'''
void iauPn06a(double date1, double date2,
              double *dpsi, double *deps, double *epsa,
              double rb[3][3], double rp[3][3], double rbp[3][3],
              double rn[3][3], double rbpn[3][3])
'''  
def pymPn06a(date1,  date2):
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
    rb : numpy.matrix(3,3)
        frame bias matrix    
    rp : numpy.matrix(3,3)
        precession matrix     
    rbp : numpy.matrix(3,3)
        bias-precession matrix    
    rn : numpy.matrix(3,3)
        nutation matrix    
    rbpn : numpy.matrix(3,3)
        GCRS-to-true matrix

    '''
    dpsi    = c_double()
    deps    = c_double()
    epsa    = c_double()
  
    rb,  rb_v   = vector_double3_ptr()
    rp,  rp_v   = vector_double3_ptr()
    rbp, rbp_v  = vector_double3_ptr()
    rn,  rn_v   = vector_double3_ptr()
    rbpn,rbpn_v = vector_double3_ptr()
    
    lib.iauPn06a(date1,  date2,  (dpsi),  (deps),  (epsa),  
                 rb_v,  rp_v,  rbp_v,  rn_v,  rbpn_v)
    
    return  dpsi[0],   deps[0],   epsa[0],  rb,  rp,  rbp,  rn,  rbpn 


'''
void iauPnm00a(double date1, double date2, double rbpn[3][3])
'''            
def pymPnm00a(date1, date2):
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
    rbpn : numpy.matrix(3,3)
        bias-precession-nutation matrix

    '''
    rbpn,rbpn_v = vector_double3_ptr()
    
    lib.iauPnm00a(date1, date2, rbpn_v)
    
    return  rbpn

'''
void iauPnm00b(double date1, double date2, double rbpn[3][3])
'''
def pymPnm00b(date1, date2):
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
    rbpn : numpy.matrix(3,3)
        bias-precession-nutation matrix

    '''
    rbpn,rbpn_v = vector_double3_ptr()
    
    lib.iauPnm00b(date1, date2, rbpn_v)
    
    return  rbpn



'''
void iauPnm06a(double date1, double date2, double rbpn[3][3])
'''
def pymPnm06a(date1, date2):
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
    rbpn : numpy.matrix(3,3)
        bias-precession-nutation matrix

    '''
    rbpn,rbpn_v = vector_double3_ptr()
   
    lib.iauPnm06a(date1, date2, rbpn_v)
    
    return  rbpn

'''
void iauPnm80(double date1, double date2, double rmatpn[3][3])
'''
def pymPnm80(date1, date2):
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
    rmatpn : numpy.matrix(3,3)
        combined precession/nutation matrix

    '''
    rbpn,rbpn_v = vector_double3_ptr()
 
    lib.iauPnm80(date1, date2, rbpn_v)
    
    return  rbpn


'''
void iauPom00(double xp, double yp, double sp, double rpom[3][3])
'''
def pymPom00(xp,   yp,   sp):
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
    rpom : numpy.matrix(3,3)
        polar-motion matrix

    '''
    rpom, rpom_v = vector_double3_ptr()
     
    lib. iauPom00(xp,   yp,   sp,  rpom_v)
    
    return  rpom


'''
double iauS00(double date1, double date2, double x, double y)
'''
def pymS00(date1,   date2,   x,   y):
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
    
    return  lib.iauS00(date1,   date2,   x,   y)


'''
void iauPr00(double date1, double date2, double *dpsipr, double *depspr)
'''
def pymPr00(date1, date2):
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
    dpsipr   = c_double()
    depspr   = c_double()
    
    lib.iauPr00(date1,  date2,  dpsipr,  depspr)
    
    return   dpsipr[0], depspr[0]


'''
void iauPrec76(double date01, double date02, double date11, double date12,
               double *zeta, double *z, double *theta)
'''
def pymPrec76(date01, date02, date11, date12):
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
    zeta   = c_double()
    z      = c_double()
    theta  = c_double()

    lib.iauPrec76(date01,  date02,  date11,  date12, 
                  (zeta), (z), (theta))
    
    return   zeta[0], z[0], theta[0]

'''
double iauS00a(double date1, double date2)
'''
def pymS00a(date1,   date2):
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
    return  lib.iauS00a(date1,   date2)

'''
double iauS00b(double date1, double date2)
'''
def pymS00b(date1,   date2):
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
    return  lib.iauS00b(date1,   date2)

'''
double iauS06(double date1, double date2, double x, double y)
'''
def pymS06(date1,   date2,   x,   y):
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
    return  lib.iauS06(date1,   date2,   x,   y)

'''
double iauS06a(double date1, double date2)
'''
def pymS06a(date1,   date2):
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
    return  lib.iauS06a(date1,   date2)

'''
double iauSp00(double date1, double date2)
'''
def pymSp00(date1,   date2):
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
    return  lib.iauSp00(date1,   date2)



'''
void iauXy06(double date1, double date2, double *x, double *y)
'''
def pymXy06(date1,   date2):
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
    x      = c_double()
    y      = c_double()
    
    lib.iauXy06(date1,   date2,  (x),  (y))
    
    return  x[0],  y[0]

'''
void iauXys00a(double date1, double date2,
               double *x, double *y, double *s)
'''
def pymXys00a(date1,   date2):
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
    x      = c_double()
    y      = c_double()
    s      = c_double()
    
    lib.iauXys00a(date1,   date2,  (x),  (y),  (s))
    
    return  x[0],  y[0],  s[0] 


'''
void iauXys00b(double date1, double date2,
               double *x, double *y, double *s)
'''
def pymXys00b(date1,   date2):
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
    x      = c_double()
    y      = c_double()
    s      = c_double()
    
    lib.iauXys00b(date1,   date2,  (x),  (y),  (s))
    
    return  x[0],  y[0],  s[0] 


'''
void iauXys06a(double date1, double date2,
               double *x, double *y, double *s)
'''
def pymXys06a(date1,   date2):
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
    x      = c_double()
    y      = c_double()
    s      = c_double()
    
    lib.iauXys06a(date1,   date2,  (x),  (y),  (s))
    
    return  x[0],  y[0],  s[0] 


#sofa vector_matrix model 
 
'''
void iauA2tf(int ndp, double angle, char *sign, int ihmsf[4])
'''
def pymA2tf(ndp, angle):
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
    sign : bytes  
        '+' or '-'    
    ihmsf : tuple
        
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
    sign   = char_ptr()
    ihmsf, ihmsf_p = int_array_4()
    
    lib.iauA2tf(ndp, angle, (sign), ihmsf_p)
    
    ihmsf  = list(np.array(ihmsf).flatten())
    
    return sign[0], tuple([x for x in ihmsf]) 


'''
int iauAf2a(char s, int ideg, int iamin, double asec, double *rad)
'''

iau_af2a_msg = {
                  1: 'ideg outside range 0-359',
                  2: 'iamin outside range 0-59',
                  3: 'asec outside range 0-59.999...',
                  }    

def pymAf2a(s,   ideg,   iamin,   asec):
    '''
    Convert degrees, arcminutes, arcseconds to radians.

    Parameters
    ----------
    s : bytes
        sign:  '-' = negative, otherwise positive    
    ideg : int
        degrees    
    iamin : int
        arcminutes    
    asec : float
        arcseconds    

    Raises
    ------
    ValueError
        1: 'ideg outside range 0-359',
        2: 'iamin outside range 0-59',
        3: 'asec outside range 0-59.999...',
        
    Returns
    -------
    rad : float
        angle in radians

    '''
    rad = c_double()
    
    j   = lib.iauAf2a(s,   ideg,   iamin,   asec,  (rad))
    
    if j > 0:
        ws.warn(iau_af2a_msg[j], UserWarning, 2)
 
    return rad[0]


'''
double iauAnpm(double a)
'''
def pymAnpm(a):
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
    return lib.iauAnpm(a)


'''
void iauC2s(double p[3], double *theta, double *phi)
'''
def pymC2s(p):
    '''
    P-vector to spherical coordinates.

    Parameters
    ----------
    p : numpy.matrix(1,3)
        p-vector

    Returns
    -------
    theta : float
        longitude angle (radians)    
    phi : float
        latitude angle (radians)

    '''
    p_v   = vector_double(p)
    
    theta = c_double()
    phi   = c_double()

    lib.iauC2s(p_v, (theta),  (phi))
    
    return theta[0], phi[0]


'''
void iauCp(double p[3], double c[3])
'''
def pymCp(p):
    '''
    Copy a p-vector.
     
    Parameters
    ----------
    p : numpy.matrix(1,3)
        p-vector to be copied
     
    Returns
    -------
    c : numpy.matrix(1,3)
        copy
     
    '''
    p_v    = vector_double(p)
    
    c, c_v = vector_double_ptr()
    
    lib.iauCp(p_v, c_v)
    
    return c

def pymCp_A(p):
    
    return p 


'''
void iauCpv(double pv[2][3], double c[2][3])
'''
def pymCpv(pv):
    '''
    Copy a position/velocity vector.

    Parameters
    ----------
    pv : numpy.matrix(2,3)
        position/velocity vector to be copied

    Returns
    -------
    c : numpy.matrix(2,3)
        copy

    '''
    pv_v    = vector_double2(pv)
    
    c, c_v  = vector_double2_ptr()
    
    lib.iauCpv(pv_v, c_v)
    
    return c

def pymCpv_A(pv):

    return pv


'''
void iauCr(double r[3][3], double c[3][3])
'''
def pymCr(r):
    '''
    Copy an r-matrix.

    Parameters
    ----------
    r : numpy.matrix(3,3)
        r-matrix to be copied

    Returns
    -------
    c : numpy.matrix(3,3)
        r-matrix to be copied

    '''
    r_v    = vector_double3(r)
    
    c, c_v  = vector_double3_ptr()

    lib.iauCr(r_v, c_v)
    
    return c

def pymCr_A(r):
    
    return r


'''
void iauIr(double r[3][3])
'''
def pymIr():
    '''
    Initialize an r-matrix to the identity matrix.

    Returns
    -------
    r : numpy.matrix(3,3)
        r-matrix

    '''
    r, r_v  = vector_double3_ptr()
    
    lib.iauIr(r_v)
    
    return r


'''
void iauP2pv(double p[3], double pv[2][3])
'''
def pymP2pv(p):
    '''
    Extend a p-vector to a pv-vector by appending a zero velocity.

    Parameters
    ----------
    p : numpy.matrix(1,3)
        p-vector

    Returns
    -------
    pv : numpy.matrix(2,3)
        pv-vector

    '''
    p_v     = vector_double(p)
    
    pv, pv_v  = vector_double2_ptr()
    
    lib.iauP2pv(p_v, pv_v)
    
    return pv



'''
void iauP2s(double p[3], double *theta, double *phi, double *r)
'''
def pymP2s(p):
    '''
    P-vector to spherical polar coordinates.

    Parameters
    ----------
    p : numpy.matrix(1,3)
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
    p_v   = vector_double(p)
    
    theta = c_double()
    phi   = c_double()
    r     = c_double()

    lib.iauP2s(p_v, (theta), (phi), (r))
    
    return theta[0], phi[0], r[0] 


'''
double iauPap(double a[3], double b[3])
'''
def pymPap(a, b):
    '''
    Position-angle from two p-vectors.

    Parameters
    ----------
    a : numpy.matrix(1,3)
        direction of reference point    
    b : numpy.matrix(1,3)
        direction of point whose PA is required

    Returns
    -------
    function value : float
        position angle of b with respect to a (radians)

    '''
    a_v   = vector_double(a)
    b_v   = vector_double(b)
    
    return lib.iauPap(a_v, b_v)


'''
double iauPdp(double a[3], double b[3])
'''
def pymPdp(a, b):
    '''
    p-vector inner (=scalar=dot) product.

    Parameters
    ----------
    a : numpy.matrix(1,3)
        first p-vector    
    b : numpy.matrix(1,3)
        second p-vector

    Returns
    -------
    function value : float
        a . b

    '''
    a_v   = vector_double(a)
    b_v   = vector_double(b)
    
    return lib.iauPdp(a_v, b_v)

#Two vectors dot to have scalar product,  using np.dot
def pymPdp_A(a, b):
 
    return  np.dot(a, b)


'''
double iauPm(double p[3])
'''
def pymPm(p):
    '''
    Modulus of p-vector.

    Parameters
    ----------
    p : numpy.matrix(1,3)
        p-vector

    Returns
    -------
    function value : float
        modulus

    '''
    p_v   = vector_double(p)
    
    return lib.iauPm(p_v)

def pymPm_A(p):
        
    return np.linalg.norm(p, ord=2)


'''
void iauPmp(double a[3], double b[3], double amb[3])
'''
def pymPmp(a, b):
    '''
    P-vector subtraction.

    Parameters
    ----------
    a : numpy.matrix(1,3)
        first p-vector    
    b : numpy.matrix(1,3)
        second p-vector

    Returns
    -------
    amb : numpy.matrix(1,3)
        a - b

    '''
    a_v   = vector_double(a)
    b_v   = vector_double(b)
    amb, amb_v  = vector_double_ptr()
    
    lib.iauPmp(a_v, b_v, amb_v)
    
    return  amb

def pymPmp_A(a, b):

    return  a - b


'''
void iauPn(double p[3], double *r, double u[3])
'''
def pymPn(p):
    '''
    Convert a p-vector into modulus and unit vector.

    Parameters
    ----------
    p : numpy.matrix(1,3)
        p-vector

    Returns
    -------
    r : float
        modulus    
    u : numpy.matrix(1,3)
        unit vector

    '''
    p_v   = vector_double(p) 
    r     = c_double()
     
    u, u_v  = vector_double_ptr()
    
    lib.iauPn(p_v, (r), u_v)
    
    return r[0], u

#Python 
def pymPn_A(p):
    
    r = np.linalg.norm(p, ord=2) 
    u = p/r 
    
    return r, u



'''
void iauPpp(double a[3], double b[3], double apb[3])
'''
def pymPpp(a, b):
    '''
    P-vector addition.

    Parameters
    ----------
    a : numpy.matrix(1,3)
        first p-vector    
    b : numpy.matrix(1,3)
        second p-vector

    Returns
    -------
    apb : numpy.matrix(1,3)
        a + b

    '''
    a_v         = vector_double(a)
    b_v         = vector_double(b)
    apb, apb_v  = vector_double_ptr()
    
    lib.iauPpp(a_v, b_v, apb_v)
    
    return  apb

def pymPpp_A(a, b):

    return  a + b


'''
void iauPpsp(double a[3], double s, double b[3], double apsb[3])
'''
def pymPpsp(a, s, b):
    '''
    P-vector plus scaled p-vector.

    Parameters
    ----------
    a : numpy.matrix(1,3)
        first p-vector    
    s : float     
        scalar (multiplier for b)
    b : numpy.matrix(1,3)
        second p-vector

    Returns
    -------
    apsb : numpy.matrix(1,3)
        a + s*b

    '''
    a_v           = vector_double(a)
    b_v           = vector_double(b)
    apsb, apsb_v  = vector_double_ptr() 
    
    lib.iauPpsp(a_v, s, b_v, apsb_v)
    
    return  apsb

def pymPpsp_A(a, s, b):
    
    return a + s*b


'''
void iauPv2p(double pv[2][3], double p[3])
'''
def pymPv2p(pv):
    '''
    Discard velocity component of a pv-vector.

    Parameters
    ----------
    pv : numpy.matrix(2,3)
        pv-vector

    Returns
    -------
    p : numpy.matrix(1,3)
        p-vector

    '''
 
    pv_v    = vector_double2(pv) 
    
    p, p_v  = vector_double_ptr() 
    
    lib.iauPv2p(pv_v, p_v)
    
    return p 


'''
void iauPv2s(double pv[2][3],
             double *theta, double *phi, double *r,
             double *td, double *pd, double *rd)
'''
def pymPv2s(pv):
    '''
    Convert position/velocity from Cartesian to spherical coordinates.

    Parameters
    ----------
    pv : numpy.matrix(2,3)
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
    pv_v  = vector_double2(pv) 
 
    theta = c_double()
    phi   = c_double()
    r     = c_double()
    
    td    = c_double()
    pd    = c_double()
    rd    = c_double()
    
    lib.iauPv2s(pv_v, (theta), (phi), (r),
                (td),    (pd),  (rd))
    
    return theta[0],  phi[0],  r[0],  td[0], pd[0], rd[0]



'''
void iauPvdpv(double a[2][3], double b[2][3], double adb[2])
'''
def pymPvdpv(a, b):
    '''
    Inner (=scalar=dot) product of two pv-vectors.

    Parameters
    ----------
    a : numpy.matrix(2,3)
        first pv-vector   
    b : numpy.matrix(2,3)
        second pv-vector   

    Returns
    -------
    adb : numpy.matrix(1,2)
        DESCRIPTION.

    '''
    a_v         = vector_double2(a) 
    b_v         = vector_double2(b) 
    adb, adb_v  = array_double_ptr()
    
    lib.iauPvdpv(a_v, b_v, adb_v)
    
    return adb 


'''
void iauPvm(double pv[2][3], double *r, double *s)
'''
def pymPvm(pv):
    '''
    Modulus of pv-vector.

    Parameters
    ----------
    pv : numpy.matrix(2,3)
        pv-vector

    Returns
    -------
    r : float
        modulus of position component
    s : float
        modulus of velocity component

    '''
    pv_v  = vector_double2(pv)  
    r     = c_double()
    s     = c_double()
   
    lib.iauPvm(pv_v, (r), (s))
    
    return r[0], s[0]



'''
void iauPvmpv(double a[2][3], double b[2][3], double amb[2][3])
'''
def pymPvmpv(a, b):
    '''
    Subtract one pv-vector from another.

    Parameters
    ----------
    a : numpy.matrix(2,3)
        first pv-vector    
    b : numpy.matrix(2,3)
        second pv-vector

    Returns
    -------
    amb : numpy.matrix(2,3)
        a - b

    '''
    a_v           = vector_double2(a)
    b_v           = vector_double2(b)
    amb,   amb_v  = vector_double2_ptr() 
    
    lib.iauPvmpv(a_v, b_v, amb_v)
    
    return amb

def pymPvmpv_A(a, b):
    
    return a - b


'''
void iauPvppv(double a[2][3], double b[2][3], double apb[2][3])
'''
def pymPvppv(a, b):
    '''
    Add one pv-vector to another.

    Parameters
    ----------
    a : numpy.matrix(2,3) 
        first pv-vector    
    b : numpy.matrix(2,3) 
        second pv-vector

    Returns
    -------
    apb : numpy.matrix(2,3) 
        a + b

    '''
    a_v           = vector_double2(a)
    b_v           = vector_double2(b)
    apb,   apb_v  = vector_double2_ptr() 
    
    lib.iauPvppv(a_v, b_v, apb_v)
    
    return apb

def pymPvppv_A(a, b):
    
    return a + b


'''
void iauPvu(double dt, double pv[2][3], double upv[2][3])
''' 
def pymPvu(dt, pv):
    '''
    Update a pv-vector.

    Parameters
    ----------
    dt : float
        time interval    
    pv : numpy.matrix(2,3) 
        pv-vector

    Returns
    -------
    upv : numpy.matrix(2,3) 
        p updated, v unchanged

    '''
    pv_v          = vector_double2(pv)
    
    upv,   upv_v  = vector_double2_ptr()  
    
    lib.iauPvu(dt, pv_v, upv_v)
    
    return upv


'''
void iauPvup(double dt, double pv[2][3], double p[3])
'''
def pymPvup(dt, pv):
    '''
    Update a pv-vector, discarding the velocity component.

    Parameters
    ----------
    dt : float
        time interval    
    pv : numpy.matrix(2,3) 
        pv-vector

    Returns
    -------
    p : numpy.matrix(1,3) 
        p-vector

    '''
    pv_v      = vector_double2(pv)
    
    p,   p_v  = vector_double_ptr()  
    
    lib.iauPvup(dt, pv_v, p_v)
    
    return p

'''
void iauPvxpv(double a[2][3], double b[2][3], double axb[2][3])
'''
def pymPvxpv(a, b):
    '''
    Outer (=vector=cross) product of two pv-vectors.

    Parameters
    ----------
    a : numpy.matrix(2,3) 
        first pv-vector    
    b : numpy.matrix(2,3) 
        second pv-vector

    Returns
    -------
    axb : numpy.matrix(2,3) 
        a x b

    '''
    a_v          = vector_double2(a)
    b_v          = vector_double2(b)
    
    axb,  axb_v  = vector_double2_ptr() 
    
    lib.iauPvxpv(a_v, b_v, axb_v)
    
    return axb


'''
void iauPxp(double a[3], double b[3], double axb[3])
'''
def pymPxp(a, b):
    '''
    p-vector outer (=vector=cross) product.

    Parameters
    ----------
    a : numpy.matrix(1,3) 
        first p-vector    
    b : numpy.matrix(1,3) 
        second p-vector

    Returns
    -------
    axb : numpy.matrix(1,3) 
        a x b

    '''
    a_v          = vector_double(a)
    b_v          = vector_double(b)
    
    axb,  axb_v  = vector_double_ptr() 
    
    lib.iauPxp(a_v, b_v, axb_v)
    
    return axb

def pymPxp_A(a, b):
    
    return np.cross(a, b)


'''
void iauRm2v(double r[3][3], double w[3])
'''
def pymRm2v(r):
    '''
    Express an r-matrix as an r-vector.

    Parameters
    ----------
    r : numpy.matrix(3,3)  
        rotation matrix

    Returns
    -------
    w : numpy.matrix(1,3)  
        rotation vector

    '''
    r_v        = vector_double3(r)
    
    w,  w_v    = vector_double_ptr() 
    
    lib.iauRm2v(r_v, w_v)
    
    return w


'''
void iauRv2m(double w[3], double r[3][3])
'''
def pymRv2m(w):
    '''
    Form the r-matrix corresponding to a given r-vector.

    Parameters
    ----------
    w : numpy.matrix(1,3)     
        rotation vector     

    Returns
    -------
    r : numpy.matrix(3,3)     
        rotation matrix    

    '''
    w_v        = vector_double(w)
    
    r,  r_v    = vector_double3_ptr() 
    
    lib.iauRv2m(w_v, r_v)
    
    return r


'''
void iauRx(double phi, double r[3][3])
'''
def pymRx(phi, r):
    '''
    Rotate an r-matrix about the x-axis.

    Parameters
    ----------
    phi : float    
        angle (radians)    
    r : numpy.matrix(3,3)    
        r-matrix    

    Returns
    -------
    r : numpy.matrix(3,3)    
        r-matrix, rotated

    '''
    r_v    = vector_double3(r) 
    
    lib.iauRx(phi, r_v)
    
    return r

#Python codes for iauRx
def pymRx_A(phi, r):
    
    cos_phi = np.math.cos(phi)
    sin_phi = np.math.sin(phi)

# to be in accordance with SOFA definition for Rx
    mat_rx  = np.array(
                  [ [1,         0,          0],
                    [0,    cos_phi,   sin_phi],
                    [0,   -sin_phi,   cos_phi]
                  ]).reshape(3,3)
    
    return np.dot(mat_rx, r)


'''
void iauRxp(double r[3][3], double p[3], double rp[3])
'''
def pymRxp(r, p):
    '''
    Multiply a p-vector by an r-matrix.

    Parameters
    ----------
    r : numpy.matrix(3,3)
        r-matrix    
    p : numpy.matrix(1,3)
        p-vector    

    Returns
    -------
    rp : numpy.matrix(3,3)
        r * p

    '''
    r_v        = vector_double3(r)
    p_v        = vector_double (p)
    
    rp,  rp_v  = vector_double_ptr() 
    
    lib.iauRxp(r_v, p_v, rp_v)
    
    return rp 

def pymRxp_A(r, p):
    
    return np.dot(r, p) 

'''
void iauRxpv(double r[3][3], double pv[2][3], double rpv[2][3])
''' 
def pymRxpv(r, pv):
    '''
    Multiply a pv-vector by an r-matrix.

    Parameters
    ----------
    r : numpy.matrix(3,3)
        r-matrix    
    pv : numpy.matrix(2,3)
        pv-vector

    Returns
    -------
    rpv : numpy.matrix(2,3)
        r * pv

    '''
    r_v        = vector_double3(r)
    pv_v       = vector_double2(pv)
    
    rpv, rpv_v = vector_double2_ptr()  
    
    lib.iauRxpv(r_v, pv_v, rpv_v)
    
    return rpv 

def pymRxpv_A(r, pv):

#transpose 3,2 => 2, 3  
    
    return np.dot(r, pv).T 


'''
void iauRxr(double a[3][3], double b[3][3], double atb[3][3])
'''
def pymRxr(a, b):
    '''
    Multiply two r-matrices.

    Parameters
    ----------
    a : numpy.matrix(3,3)
        first r-matrix    
    b : numpy.matrix(3,3)
        second r-matrix

    Returns
    -------
    atb : numpy.matrix(3,3)
        a * b

    '''
    a_v          = vector_double3(a)
    b_v          = vector_double3(b)
    
    atb,  atb_v  = vector_double3_ptr() 
    
    lib.iauRxr(a_v, b_v, atb_v)
    
    return atb


def pymRxr_A(a, b):
    
    return np.dot(a, b)



'''
void iauRy(double theta, double r[3][3])
'''
def pymRy(theta, r):
    '''
    Rotate an r-matrix about the y-axis.

    Parameters
    ----------
    theta : float
        angle (radians)    
    r : numpy.matrix(3,3)
        r-matrix

    Returns
    -------
    r : numpy.matrix(3,3)
        r-matrix, rotated

    '''
    r_v    = vector_double3(r) 
    
    lib.iauRy(theta, r_v)
    
    return r

#Python codes for iauRy
def pymRy_A(theta, r):
    
    cos_theta = np.math.cos(theta)
    sin_theta = np.math.sin(theta)

# to be in accordance with SOFA definition for Ry
    mat_ry  = np.array(
                  [ [cos_theta,    0,      -sin_theta],
                    [0,            1,              0 ],
                    [sin_theta,    0,       cos_theta]
                  ]).reshape(3,3)
    
    return np.dot(mat_ry, r)


'''
void iauRz(double psi, double r[3][3])
'''
def pymRz(psi, r):
    '''
    Rotate an r-matrix about the z-axis.

    Parameters
    ----------
    psi : float
        angle (radians)    
    r : numpy.matrix(3,3)
        r-matrix

    Returns
    -------
    r : numpy.matrix(3,3)
        r-matrix, rotated

    '''
    r_v    = vector_double3(r) 
    
    lib.iauRz(psi, r_v)
    
    return r

#Python codes for iauRz
def pymRz_A(psi, r):
    
    cos_psi = np.math.cos(psi)
    sin_psi = np.math.sin(psi)

# to be in accordance with SOFA definition for Ry
    mat_rz  = np.array(
                  [ [ cos_psi,     sin_psi,     0],
                    [-sin_psi,     cos_psi,     0],
                    [0,            0,           1]
                  ]).reshape(3,3)
    
    return np.dot(mat_rz, r)


'''
void iauS2p(double theta, double phi, double r, double p[3])
'''    
def pymS2p(theta, phi, r):
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
    p : numpy.matrix(1,3)
        Cartesian coordinates

    '''
    p, p_v     = vector_double_ptr()
   
    lib.iauS2p(theta, phi, r, p_v)
    
    return p


'''
void iauS2pv(double theta, double phi, double r,
             double td, double pd, double rd,
             double pv[2][3])
'''
def pymS2pv(theta, phi, r, td, pd, rd):
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
    pv : numpy.matrix(2,3)
        pv-vector

    '''
    pv, pv_v     = vector_double2_ptr()
    
    lib.iauS2pv(theta, phi, r, td, pd, rd, pv_v)
    
    return pv



'''
void iauS2xpv(double s1, double s2, double pv[2][3], double spv[2][3])
'''
def pymS2xpv(s1, s2, pv):
    '''
    Multiply a pv-vector by two scalars.

    Parameters
    ----------
    s1 : float
        scalar to multiply position component by    
    s2 : float
        scalar to multiply velocity component by    
    pv : numpy.matrix(2,3)
        pv-vector

    Returns
    -------
    spv : numpy.matrix(2,3)
        pv-vector: p scaled by s1, v scaled by s2

    '''
    pv_v         = vector_double2(pv)
    spv, spv_v   = vector_double2_ptr()
    
    lib.iauS2xpv(s1, s2, pv_v, spv_v)
    
    return spv

def pymS2xpv_A(s1, s2, pv):
   
    spv = np.asmatrix(np.zeros(shape=(2,3), dtype=float, order='C')) 
    
    spv[0] = s1 * pv[0]
    spv[1] = s2 * pv[1]
    
    return   spv


'''
double iauSepp(double a[3], double b[3])
'''
def pymSepp(a, b):
    '''
    Angular separation between two p-vectors.

    Parameters
    ----------
    a : numpy.matrix(1,3)
        first p-vector (not necessarily unit length)    
    b : numpy.matrix(1,3)
        second p-vector (not necessarily unit length)

    Returns
    -------
    function value : float
        angular separation (radians, always positive)

    '''
    a_v      = vector_double(a)
    b_v      = vector_double(b)
    
    return lib.iauSepp(a_v, b_v)

#Python codes for iauSepp
def pymSepp_A(a, b):
     
    axb     = np.cross(a, b)
    sin_phi = np.linalg.norm(axb, ord=2)
    cos_phi = np.dot(a, b)
        
    return  np.math.atan2(sin_phi, cos_phi)


'''
double iauSeps(double al, double ap, double bl, double bp)
'''
def pymSeps(al,   ap,   bl,   bp):
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
    return lib.iauSeps(al,   ap,   bl,   bp)



'''
void iauSxp(double s, double p[3], double sp[3])
'''
def pymSxp(s, p):
    '''
    Multiply a p-vector by a scalar.

    Parameters
    ----------
    s : float
        scalar    
    p : numpy.matrix(1,3)
        p-vector

    Returns
    -------
    sp : numpy.matrix(1,3)
        s * p

    '''
    p_v         = vector_double(p)
    
    sp, sp_v   = vector_double_ptr()
    
    lib.iauSxp(s, p_v, sp_v)
    
    return sp 

 
def pymSxp_A(s, p):

    return   np.asmatrix(s*p)


'''
void iauSxpv(double s, double pv[2][3], double spv[2][3])
'''    
def pymSxpv(s, pv):
    '''
    Multiply a pv-vector by a scalar.

    Parameters
    ----------
    s : float
        scalar    
    pv : numpy.matrix(2,3)
        pv-vector

    Returns
    -------
    spv : numpy.matrix(2,3)
        s * pv

    '''
    pv_v         = vector_double2(pv)
    
    spv, spv_v   = vector_double2_ptr()
    
    lib.iauSxpv(s, pv_v, spv_v)
    
    return spv 


'''
int iauTf2a(char s, int ihour, int imin, double sec, double *rad)
'''

iau_tf2a_msg = {
                1:'ihour outside range 0-23',
                2:'imin outside range 0-59',
                3:'sec outside range 0-59.999...'
                }   
 
def pymTf2a(s, ihour, imin, sec): 
    '''
    Convert hours, minutes, seconds to radians.

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
        
    Raises
    ------
    ValueError
        1:'ihour outside range 0-23',
        2:'imin outside range 0-59',
        3:'sec outside range 0-59.999...'

    Returns
    -------
    rad : float
        angle in radians

    '''
    rad = c_double()
     
    j = lib.iauTf2a(s, ihour, imin, sec, (rad))
    if j > 0:
        ws.warn(iau_tf2a_msg[j], UserWarning, 2)
        
    return rad[0]

'''
void iauTr(double r[3][3], double rt[3][3])
'''
def pymTr(r):
    '''
    Transpose an r-matrix.

    Parameters
    ----------
    r : numpy.matrix(3,3)
        r-matrix

    Returns
    -------
    rt : numpy.matrix(3,3)
        transpose

    '''
    r_v        = vector_double3(r)  
    
    rt, rt_v   = vector_double3_ptr()
    
    lib.iauTr(r_v, rt_v)
    
    return rt

def pymTr_A(r):
    
    return r.T


'''
void iauTrxp(double r[3][3], double p[3], double trp[3])
'''
def pymTrxp(r, p):
    '''
    Multiply a p-vector by the transpose of an r-matrix.

    Parameters
    ----------
    r : numpy.matrix(3,3)
        r-matrix    
    p : numpy.matrix(1,3)
        p-vector

    Returns
    -------
    trp : numpy.matrix(1,3)
        r^T * p

    '''
    r_v        = vector_double3(r)  
    p_v        = vector_double(p)  
    
    trp, trp_v   = vector_double_ptr()
    
    lib.iauTrxp(r_v, p_v, trp_v)
    
    return trp 

def pymTrxp_A(r, p):
   
#   trp  = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C'))
    trp  = np.dot(r.T, p)
    
    return np.asmatrix(trp)


'''
void iauTrxpv(double r[3][3], double pv[2][3], double trpv[2][3])
'''
def pymTrxpv(r, pv):
    '''
    Multiply a pv-vector by the transpose of an r-matrix.

    Parameters
    ----------
    r : numpy.matrix(3,3)
        r-matrix    
    pv : numpy.matrix(2,3)
        pv-vector

    Returns
    -------
    trpv : numpy.matrix(2,3)
        r^T * pv

    '''
    r_v           = vector_double3(r)  
    pv_v          = vector_double2(pv)  
    
    trpv, trpv_v  = vector_double2_ptr()
    
    lib.iauTrxpv(r_v, pv_v, trpv_v)
    
    return trpv


def pymTrxpv_A(r, pv):
    
    return np.dot(r.T, pv.T).T 



'''
void iauZp(double p[3])
'''
def pymZp(p):
    '''
    Zero a p-vector

    Parameters
    ----------
    p : numpy.matrix(1,3)
        p-vector

    Returns
    -------
    p : numpy.matrix(1,3)
        zero p-vector
    '''
    p, p_v  = vector_double_ptr()
    
    lib.iauZp(p_v)
    
    return p

def pymZp_A(p):
   
    pt  = np.asmatrix(np.zeros(shape=(1,3), dtype=float, order='C'))
     
    return pt


'''
void iauZpv(double pv[2][3])
'''
def pymZpv(pv):
    '''
    Zero a pv-vector

    Parameters
    ----------
    pv : numpy.matrix(2,3)
        pv-vector
        
    Returns
    -------
    pv : numpy.matrix(2,3)
        zero pv-vector

    '''
    pv, pv_v  = vector_double2_ptr()
    
    lib.iauZpv(pv_v)
    
    return pv

def pymZpv_A(pv):
   
    pvt  = np.asmatrix(np.zeros(shape=(2,3), dtype=float, order='C'))
     
    return pvt


'''
void iauZr(double r[3][3])
'''
def pymZr(r):
    '''
    Initialize an r-matrix to the null matrix.

    Parameters
    ----------
    r : numpy.matrix(3,3)
        r-matrix

    Returns
    -------
    r : numpy.matrix(3,3)
        null matrix

    '''
    r, r_v  = vector_double3_ptr()
 
    lib.iauZr(r_v)
    
    return r


def pymZr_A(r):
   
    rt  = np.asmatrix(np.zeros(shape=(3,3), dtype=float, order='C'))
     
    return rt

#2023-09-28
#void iauBp00(double date1, double date2,
#             double rb[3][3], double rp[3][3], double rbp[3][3]);

'''
    r_v        = vector_double3(r)  
    p_v        = vector_double(p)  
    
    trp, trp_v   = vector_double_ptr()
    
    lib.iauTrxp(r_v, p_v, trp_v)
'''

def pymBp00(date1,  date2):
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
    rb : numpy.matrix(3,3)
        frame bias matrix     
    rp : numpy.matrix(3,3)   
        precession matrix    
    rbp : numpy.matrix(3,3)
        bias-precession matrix
    '''
    rb,  rb_v  = vector_double3_ptr()
    rp,  rp_v  = vector_double3_ptr()
    rbp, rbp_v = vector_double3_ptr()
 
    lib.iauBp00(date1,  date2,  rb_v,  rp_v,  rbp_v)
    
    return   rb,  rp,  rbp 
    

#void iauBp06(double date1, double date2,
#             double rb[3][3], double rp[3][3], double rbp[3][3]);   
def pymBp06(date1,  date2):
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
    rb : numpy.matrix(3,3)
        frame bias matrix    
    rp : numpy.matrix(3,3)
        precession matrix    
    rbp : numpy.matrix(3,3)
        bias-precession matrix

    '''
    rb,  rb_v  = vector_double3_ptr()
    rp,  rp_v  = vector_double3_ptr()
    rbp, rbp_v = vector_double3_ptr()
 
    lib.iauBp06(date1,  date2,  rb_v,  rp_v,  rbp_v)
    
    return   rb,  rp,  rbp 

#2023-10-01
#void iauC2i06a(double date1, double date2, double rc2i[3][3])   
def pymC2i06a(date1, date2):
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
    rc2i : numpy.matrix(3,3)
        celestial-to-intermediate matrix

    '''
    rc2i, rc2i_v = vector_double3_ptr()
     
    lib.iauC2i06a(date1, date2, rc2i_v)
    
    return  rc2i 


#void iauC2ixy(double date1, double date2, double x, double y,
#              double rc2i[3][3])
def pymC2ixy(date1, date2, x, y):
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
    rc2i : numpy.matrix(3,3)
        celestial-to-intermediate matrix

    '''
    rc2i, rc2i_v = vector_double3_ptr()
    
    lib.iauC2ixy(date1, date2, x, y, rc2i_v)
    
    return  rc2i 


#void iauC2ixys(double x, double y, double s, double rc2i[3][3])
def pymC2ixys(x, y, s):
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
    rc2i : numpy.matrix(3,3)
        celestial-to-intermediate matrix

    '''
    rc2i, rc2i_v = vector_double3_ptr()
    
    lib.iauC2ixys(x, y, s, rc2i_v)
    
    return  rc2i 


#void iauC2t06a(double tta, double ttb, double uta, double utb,
#               double xp, double yp, double rc2t[3][3])
def pymC2t06a(tta,  ttb,  uta,  utb,  xp,  yp):
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
    rc2t : numpy.matrix(3,3)
        celestial-to-terrestrial matrix

    '''
    
    rc2t , rc2t_v = vector_double3_ptr()
    
    lib.iauC2t06a(tta,  ttb,  uta,  utb,  xp,  yp, rc2t_v)
    
    return  rc2t


#2023-10-02
#void iauC2tcio(double rc2i[3][3], double era, double rpom[3][3],
#               double rc2t[3][3])
def pymC2tcio(rc2i, era, rpom):
    '''
    Assemble the celestial to terrestrial matrix from CIO-based
    components (the celestial-to-intermediate matrix, the Earth Rotation
    Angle and the polar motion matrix).

    Parameters
    ----------
    rc2i : numpy.matrix(3,3)
        celestial-to-intermediate matrix     
    era : float
        Earth rotation angle (radians)     
    rpom : numpy.matrix(3,3)
        polar-motion matrix

    Returns
    -------
    rc2t : numpy.matrix(3,3)
        celestial-to-terrestrial matrix

    '''    
    rc2i_v        = vector_double3(rc2i)
    rpom_v        = vector_double3(rpom)
    
    rc2t , rc2t_v = vector_double3_ptr()
    
    lib.iauC2tcio(rc2i_v, era, rpom_v, rc2t_v)
    
    return  rc2t  


#void iauEceq06(double date1, double date2, double dl, double db,
#               double *dr, double *dd)  
def pymEceq06(date1,   date2,  dl,  db):
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
    dr     = c_double()
    dd     = c_double()
    
    lib.iauEceq06(date1,   date2,   dl,   db,   (dr),   (dd))
    
    return  dr[0],  dd[0] 


#void iauEcm06(double date1, double date2, double rm[3][3])
def pymEcm06(date1, date2):
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
    rm : numpy.matrix(3,3)
        ICRS to ecliptic rotation matrix

    '''
    rm , rm_v = vector_double3_ptr()
      
    lib.iauEcm06(date1, date2, rm_v)
    
    return  rm 


#int iauEform ( int n, double *a, double *f )
def pymEform(n):
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

    '''
    a     = c_double()
    f     = c_double()
    
    lib.iauEform(n,   (a),   (f))
    
    return  a[0],  f[0]


#2023-10-03
#int iauEpv00(double date1, double date2,
#             double pvh[2][3], double pvb[2][3])     
def pymEpv00(date1, date2): 
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
    pvh : numpy.matrix(2,3)
        heliocentric Earth position/velocity    
    pvb : numpy.matrix(2,3)
        barycentric Earth position/velocity

    '''
    pvh, pvh_v = vector_double2_ptr()
    pvb, pvb_v = vector_double2_ptr()
    
    lib.iauEpv00(date1, date2, pvh_v, pvb_v)
    
    return  pvh, pvb


#void iauEqec06(double date1, double date2, double dr, double dd,
#               double *dl, double *db)

def pymEqec06(date1,  date2,  dr,  dd):
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
   
    dl     = c_double()
    db     = c_double()
  
    
    lib.iauEqec06(date1,   date2,  dr,  dd,   (dl),  (db))
    
    return  dl[0],  db[0] 


#void iauFk5hip(double r5h[3][3], double s5h[3])    
def pymFk5hip():
    '''
    FK5 to Hipparcos rotation and spin.

    Returns
    -------
    r5h : numpy.matrix(3,3)
        r-matrix: FK5 rotation wrt Hipparcos    
    s5h : numpy.matrix(1,3)
        r-vector: FK5 spin wrt Hipparcos

    '''
    r5h, r5h_v = vector_double3_ptr()
    s5h, s5h_v = vector_double_ptr()
    
    lib.iauFk5hip(r5h_v, s5h_v)
    
    return  r5h, s5h 

#void iauFk5hz(double r5, double d5, double date1, double date2,
#              double *rh, double *dh)
def pymFk5hz(r5,   d5,   date1,   date2):
    '''
    Transform an FK5 (J2000.0) star position into the system of the
    Hipparcos catalogue, assuming zero Hipparcos proper motion.

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
    rh     = c_double()
    dh     = c_double()
    
    lib.iauFk5hz(r5,   d5,   date1,   date2,  (rh),   (dh))
    
    return  rh[0],  dh[0]


#void iauFk45z(double r1950, double d1950, double bepoch,
#              double *r2000, double *d2000)
def pymFk45z(r1950,  d1950,  bepoch):
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
    r2000    = c_double()
    d2000    = c_double()
    
    lib.iauFk45z(r1950,  d1950,  bepoch,  (r2000),  (d2000))
    
    return  r2000[0],  d2000[0]   


#void iauFk52h(double r5, double d5,
#              double dr5, double dd5, double px5, double rv5,
#              double *rh, double *dh,
#              double *drh, double *ddh, double *pxh, double *rvh)       
def pymFk52h(r5,  d5,  dr5,  dd5,  px5,  rv5):
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
    rh     = c_double()
    dh     = c_double()
    drh    = c_double()
    ddh    = c_double()
    pxh    = c_double()
    rvh    = c_double()
    
    lib.iauFk52h(r5,   d5,   dr5,  dd5,  px5,  rv5,  
                 (rh),  (dh),
                 (drh),  (ddh), (pxh),  (rvh))
    
    return  rh[0],  dh[0],  drh[0],  ddh[0],  pxh[0],  rvh[0] 



#void iauFk54z(double r2000, double d2000, double bepoch,
#              double *r1950, double *d1950,
#              double *dr1950, double *dd1950)    
def pymFk54z(r2000,  d2000,  bepoch):
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
    r1950    = c_double()
    d1950    = c_double()
    dr1950   = c_double()
    dd1950   = c_double()
    
    lib.iauFk54z(r2000,  d2000,  bepoch,  (r1950),  (d1950),
                 (dr1950),  (dd1950))
    
    return  r1950[0],  d1950[0],  dr1950[0],  dd1950[0] 


#void iauFk425(double r1950,  double d1950,
#              double dr1950, double dd1950,
#              double p1950,  double v1950,
#              double *r2000, double *d2000,
#              double *dr2000,double *dd2000,
#              double *p2000, double *v2000)   
    
def pymFk425(r1950,  d1950,  dr1950,  dd1950,  p1950,  v1950):
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
    r2000    = c_double()
    d2000    = c_double()
    
    dr2000   = c_double()
    dd2000   = c_double()
    
    p2000    = c_double()
    v2000    = c_double()
   
    lib.iauFk425(r1950,  d1950,  dr1950,  dd1950,  p1950,  v1950,
                 (r2000),   (d2000),
                 (dr2000),  (dd2000),
                 (p2000),   (v2000))
    
    return  r2000[0], d2000[0], dr2000[0], dd2000[0], p2000[0], v2000[0] 



#void iauFk524(double r2000, double d2000,
#              double dr2000, double dd2000,
#              double p2000, double v2000,
#              double *r1950, double *d1950,
#              double *dr1950, double *dd1950,
#              double *p1950, double *v1950) 
    
def pymFk524(r2000, d2000, dr2000, dd2000, p2000, v2000):
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
    r1950    = c_double()
    d1950    = c_double()
    
    dr1950   = c_double()
    dd1950   = c_double()
    
    p1950    = c_double()
    v1950    = c_double()
    
    
    lib.iauFk524(r2000, d2000, dr2000, dd2000, p2000, v2000,
                 (r1950),   (d1950),
                 (dr1950),  (dd1950),
                 (p1950),   (v1950))
    
    return  r1950[0], d1950[0], dr1950[0], dd1950[0], p1950[0], v1950[0]


#void iauG2icrs (double dl, double db, double *dr, double *dd )

def pymG2icrs(dl, db):
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
 
    dr     = c_double()
    dd     = c_double()
    
    lib.iauG2icrs(dl, db,  (dr),  (dd))
    
    return  dr[0],  dd[0]


#int iauGc2gd ( int n, double xyz[3],
#               double *elong, double *phi, double *height )
    
def pymGc2gd(n, xyz):
    '''
    Transform geocentric coordinates to geodetic using the specified
    reference ellipsoid.

    Parameters
    ----------
    n : int
        ellipsoid identifier    
    xyz : numpy.matrix(1,3)
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

    '''
    xyz_v     = vector_double(xyz)
    
    elong     = c_double()
    phi       = c_double()
    height    = c_double()
    
    lib.iauGc2gd(n, xyz_v,  (elong),  (phi), (height))
    
    return  elong[0],  phi[0],  height[0]


#int iauGc2gde ( double a, double f, double xyz[3],
#                double *elong, double *phi, double *height )
def pymGc2gde(a, f, xyz):
    '''
    Transform geocentric coordinates to geodetic for a reference
    ellipsoid of specified form.

    Parameters
    ----------
    a : float
        equatorial radius    
    f : float
        flattening    
    xyz : numpy.matrix(1,3)
        geocentric vector

    Returns
    -------
    elong : float
        longitude (radians, east +ve)    
    phi : float
        latitude (geodetic, radians)    
    height : float
        height above ellipsoid (geodetic)

    '''
    xyz_v     = vector_double(xyz)
    
    elong     = c_double()
    phi       = c_double()
    height    = c_double()
    
    lib.iauGc2gde(a, f, xyz_v,  (elong),  (phi),  (height))
    
    return  elong[0],  phi[0],  height[0] 


#int iauGd2gc ( int n, double elong, double phi, double height,
#               double xyz[3] )     
def pymGd2gc(n, elong,  phi,  height):
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
    xyz : numpy.matrix(1,3)
        geocentric vector

    '''
    xyz, xyz_v  =  vector_double_ptr()
    
    lib.iauGd2gc(n,  elong,  phi,  height,  xyz_v)
    
    return  xyz


#int iauGd2gce ( double a, double f, double elong, double phi,
#                double height, double xyz[3] )
def pymGd2gce(a,  f,  elong,  phi,  height):
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
    xyz : numpy.matrix(1,3)
        geocentric vector

    '''
    xyz, xyz_v  =  vector_double_ptr()
    
    lib.iauGd2gce(a,  f,   elong,  phi,  height,  xyz_v)
    
    return  xyz  



#void iauH2fk5(double rh, double dh,
#              double drh, double ddh, double pxh, double rvh,
#              double *r5, double *d5,
#              double *dr5, double *dd5, double *px5, double *rv5) 
def pymH2fk5(rh,  dh,  drh,  ddh,  pxh,  rvh):
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
    r5     = c_double()
    d5     = c_double()
    dr5    = c_double()
    dd5    = c_double()
    px5    = c_double()
    rv5    = c_double()
    
    
    lib.iauH2fk5(rh,  dh,  drh,  ddh,  pxh,  rvh,  
                 (r5),  (d5),
                 (dr5),  (dd5), (px5), (rv5))
    
    return  r5[0],  d5[0],  dr5[0],  dd5[0],  px5[0],  rv5[0]



#void iauHfk5z(double rh, double dh, double date1, double date2,
#              double *r5, double *d5, double *dr5, double *dd5)
def pymHfk5z(rh,   dh,   date1,   date2):
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
    r5     = c_double()
    d5     = c_double()
    dr5    = c_double()
    dd5    = c_double()
    
    lib.iauHfk5z(rh,   dh,   date1,   date2,  
                 (r5),  (d5), (dr5),  (dd5))
    
    return  r5[0],  d5[0], dr5[0],  dd5[0] 


#void iauIcrs2g ( double dr, double dd, double *dl, double *db ) 
def pymIcrs2g(dr, dd):
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
    dl    = c_double()
    db    = c_double()
    
    lib.iauIcrs2g(dr,  dd,  (dl),  (db))
    
    return  dl[0],  db[0]   
 


#void iauLteceq(double epj, double dl, double db, double *dr, double *dd)
def pymLteceq(epj, dl, db):
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
    dr     = c_double()
    dd     = c_double()
    
    lib.iauLteceq(epj, dl, db,  (dr),  (dd))
    
    return  dr[0],  dd[0] 


#void iauLtecm(double epj, double rm[3][3])
def pymLtecm(epj):
    '''
    ICRS equatorial to ecliptic rotation matrix, long-term.

    Parameters
    ----------
    epj : float
        Julian epoch (TT)

    Returns
    -------
    rm : numpy.matrix(3,3)
        ICRS to ecliptic rotation matrix

    '''
    rm, rm_v = vector_double3_ptr()
       
    lib.iauLtecm(epj, rm_v)
    
    return rm  


#void iauLteqec(double epj, double dr, double dd, double *dl, double *db)
def pymLteqec(epj, dr, dd):
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
    dl    = c_double()
    db    = c_double()
    
    lib.iauLteqec(epj, dr,  dd,  (dl),  (db))
    
    return  dl[0],  db[0]  


#void iauLtp(double epj, double rp[3][3])
def pymLtp(epj):
    '''
    Long-term precession matrix.

    Parameters
    ----------
    epj : float
        Julian epoch (TT)     

    Returns
    -------
    rp : numpy.matrix(3,3)
        precession matrix, J2000.0 to date

    '''
    rp, rp_v = vector_double3_ptr()
    
    lib.iauLtp(epj, rp_v)
    
    return rp 


#void iauLtpb(double epj, double rpb[3][3])
def pymLtpb(epj):
    '''
    Long-term precession matrix, including ICRS frame bias.

    Parameters
    ----------
    epj : float
        Julian epoch (TT)    

    Returns
    -------
    rpb : numpy.matrix(3,3)
        precession-bias matrix, J2000.0 to date

    '''
    rpb, rpb_v = vector_double3_ptr()
    
    lib.iauLtpb(epj, rpb_v)
    
    return rpb


#void iauLtpecl(double epj, double vec[3])
def pymLtpecl(epj):
    '''
    Long-term precession of the ecliptic.

    Parameters
    ----------
    epj : float
        Julian epoch (TT)    
        
    Returns
    -------
    vec : numpy.matrix(1,3)
        ecliptic pole unit vector    

    '''
    vec, vec_v  =  vector_double_ptr()
    
    lib.iauLtpecl(epj, vec_v)
    
    return vec


#void iauLtpequ(double epj, double veq[3])
def pymLtpequ(epj):
    '''
    Long-term precession of the equator.

    Parameters
    ----------
    epj : float
        Julian epoch (TT)   

    Returns
    -------
    veq : numpy.matrix(1,3)
        equator pole unit vector

    '''
    veq, veq_v  =  vector_double_ptr()
    
    lib.iauLtpequ(epj, veq_v)
    
    return veq  


#void iauMoon98 ( double date1, double date2, double pv[2][3] )
def pymMoon98(date1, date2): 
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
    pv : numpy.matrix(2,3)
        Moon p,v, GCRS (AU, AU/d)

    '''
    pv, pv_v  =  vector_double2_ptr()
    
    lib.iauMoon98(date1, date2, pv_v)
    
    return  pv  


#void iauPfw06(double date1, double date2,
#              double *gamb, double *phib, double *psib, double *epsa)                                 
def pymPfw06(date1,   date2):
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
    gamb     = c_double()
    phib     = c_double()
    psib     = c_double()
    epsa     = c_double()
 
    
    lib.iauPfw06(date1,   date2,   
                 (gamb),  (phib), (psib),  (epsa))
    
    return  gamb[0],  phib[0],  psib[0],  epsa[0] 


#int iauPlan94(double date1, double date2, int np, double pv[2][3])
def pymPlan94(date1, date2, n): 
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
    pv : numpy.matrix(2,3)
        planet p,v (heliocentric, J2000.0, au,au/d)

    '''
    pv, pv_v  =  vector_double2_ptr()
    
    lib.iauPlan94(date1, date2, n, pv_v)
    
    return  pv 


#void iauPmat00(double date1, double date2, double rbp[3][3])
def pymPmat00(date1, date2): 
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
    rbp : numpy.matrix(3,3)
        bias-precession matrix

    '''
    rbp, rbp_v = vector_double3_ptr()
    
    lib.iauPmat00(date1, date2,  rbp_v)
    
    return  rbp  


#void iauPmat06(double date1, double date2, double rbp[3][3]) 
def pymPmat06(date1, date2): 
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
    rbp : numpy.matrix(3,3)
        bias-precession matrix

    '''
    rbp, rbp_v = vector_double3_ptr()
    
    lib.iauPmat06(date1, date2,  rbp_v)
    
    return  rbp  


#void iauPmat76(double date1, double date2, double rmatp[3][3])
def pymPmat76(date1, date2): 
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
    rmatp : numpy.matrix(3,3)
        precession matrix, J2000.0 -> date1+date2

    '''    
    rmatp, rmatp_v = vector_double3_ptr()
    
    lib.iauPmat76(date1, date2,  rmatp_v)
    
    return  rmatp 



#void iauPn00(double date1, double date2, double dpsi, double deps,
#             double *epsa,
#             double rb[3][3], double rp[3][3], double rbp[3][3],
#             double rn[3][3], double rbpn[3][3])

def pymPn00(date1,  date2,  dpsi,  deps):
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
    rb : numpy.matrix(3,3)
        frame bias matrix    
    rp : numpy.matrix(3,3)
        precession matrix    
    rbp : numpy.matrix(3,3)
        bias-precession matrix    
    rn : numpy.matrix(3,3)
        nutation matrix    
    rbpn : numpy.matrix(3,3)
        GCRS-to-true matrix
        
    '''   
    epsa    = c_double()
    
    rb,   rb_v   = vector_double3_ptr()
    rp,   rp_v   = vector_double3_ptr()
    rbp,  rbp_v  = vector_double3_ptr()
    rn,   rn_v   = vector_double3_ptr()
    rbpn, rbpn_v = vector_double3_ptr()
 
    lib.iauPn00(date1,  date2,   dpsi,  deps,  (epsa),  
                 rb_v,  rp_v,  rbp_v,  rn_v,  rbpn_v)
    
    return   epsa[0],  rb,  rp,  rbp,  rn,  rbpn   



#void iauPn00a(double date1, double date2,
#              double *dpsi, double *deps, double *epsa,
#              double rb[3][3], double rp[3][3], double rbp[3][3],
#              double rn[3][3], double rbpn[3][3]) 

def pymPn00a(date1,  date2):
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
    rb : numpy.matrix(3,3)
        frame bias matrix    
    rp : numpy.matrix(3,3)
        precession matrix    
    rbp : numpy.matrix(3,3)
        bias-precession matrix    
    rn : numpy.matrix(3,3)
        nutation matrix    
    rbpn : numpy.matrix(3,3)
        GCRS-to-true matrix

    '''
    dpsi    = c_double()
    deps    = c_double()
    epsa    = c_double()
    
    rb,   rb_v   = vector_double3_ptr()
    rp,   rp_v   = vector_double3_ptr()
    rbp,  rbp_v  = vector_double3_ptr()
    rn,   rn_v   = vector_double3_ptr()
    rbpn, rbpn_v = vector_double3_ptr()
    
    lib.iauPn00a(date1,  date2,  (dpsi),  (deps),  (epsa),  
                 rb_v,  rp_v,  rbp_v,  rn_v,  rbpn_v)
    
    return  dpsi[0],   deps[0],   epsa[0],  rb,  rp,  rbp,  rn,  rbpn    
 

#void iauPn00b(double date1, double date2,
#              double *dpsi, double *deps, double *epsa,
#              double rb[3][3], double rp[3][3], double rbp[3][3],
#              double rn[3][3], double rbpn[3][3])
def pymPn00b(date1,  date2):
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
    dpsi : flaot
        nutation    
    deps : float
        nutation    
    epsa : flaot
        mean obliquity    
    rb : numpy.matrix(3,3)
        frame bias matrix    
    rp : numpy.matrix(3,3)
        precession matrix    
    rbp : numpy.matrix(3,3)
        bias-precession matrix    
    rn : numpy.matrix(3,3)
        nutation matrix    
    rbpn : numpy.matrix(3,3)
        GCRS-to-true matrix

    '''
    dpsi    = c_double()
    deps    = c_double()
    epsa    = c_double()
    
    rb,   rb_v   = vector_double3_ptr()
    rp,   rp_v   = vector_double3_ptr()
    rbp,  rbp_v  = vector_double3_ptr()
    rn,   rn_v   = vector_double3_ptr()
    rbpn, rbpn_v = vector_double3_ptr()
    
    lib.iauPn00b(date1,  date2,  (dpsi),  (deps),  (epsa),  
                 rb_v,  rp_v,  rbp_v,  rn_v,  rbpn_v)
    
    return  dpsi[0],   deps[0],   epsa[0],  rb,  rp,  rbp,  rn,  rbpn   



#void iauPn06(double date1, double date2, double dpsi, double deps,
#             double *epsa,
#             double rb[3][3], double rp[3][3], double rbp[3][3],
#             double rn[3][3], double rbpn[3][3])
def pymPn06(date1,  date2,  dpsi,  deps):
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
    rb : numpy.matrix(3,3)
        frame bias matrix    
    rp : numpy.matrix(3,3)
        precession matrix    
    rbp : numpy.matrix(3,3)
        bias-precession matrix    
    rn : numpy.matrix(3,3)
        nutation matrix    
    rbpn : numpy.matrix(3,3)
        GCRS-to-true matrix

    '''
    epsa    = c_double()
    
    rb,   rb_v   = vector_double3_ptr()
    rp,   rp_v   = vector_double3_ptr()
    rbp,  rbp_v  = vector_double3_ptr()
    rn,   rn_v   = vector_double3_ptr()
    rbpn, rbpn_v = vector_double3_ptr()
 
    lib.iauPn06(date1,  date2,   dpsi,  deps,  (epsa),  
                 rb_v,  rp_v,  rbp_v,  rn_v,  rbpn_v)
    
    return   epsa[0],  rb,  rp,  rbp,  rn,  rbpn   


iau_pvstar_msg = {
                  0: 'OK',
                 -1: 'superluminal speed',
                 -2: 'null position vector' 
                  }  
#int iauPvstar(double pv[2][3], double *ra, double *dec,
#              double *pmr, double *pmd, double *px, double *rv)  
def pymPvstar(pv):
    '''
    Convert star position+velocity vector to catalog coordinates.

    Parameters
    ----------
    pv : numpy.matrix(2,3)
        pv-vector (au, au/day)

    Raises
    ------
    ValueError
        0: 'OK',
       -1: 'superluminal speed',
       -2: 'null position vector' 

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

    '''
    pv_v = vector_double2(pv)
    
    ra   = c_double()
    dec  = c_double()
    pmr  = c_double()
    pmd  = c_double()
    
    px   = c_double()
    rv   = c_double()
   
    j = lib.iauPvstar(pv_v,
                      (ra),  (dec),  (pmr),  (pmd),
                      (px),  (rv))

    if   j < 0:
        raise ValueError(iau_pvstar_msg[j])

        
    return ra[0], dec[0], pmr[0], pmd[0], px[0], rv[0]



#int iauStarpm(double ra1, double dec1,
#              double pmr1, double pmd1, double px1, double rv1,
#              double ep1a, double ep1b, double ep2a, double ep2b,
#              double *ra2, double *dec2,
#              double *pmr2, double *pmd2, double *px2, double *rv2)     
 
def pymStarpm(ra1,  dec1, pmr1, pmd1,
              px1,  rv1,
              ep1a, ep1b, ep2a, ep2b):
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
    rv2 : flaot
        radial velocity (km/s, +ve = receding), after    

    '''
    ra2   = c_double()
    dec2  = c_double()
    pmr2  = c_double()
    pmd2  = c_double()
    
    px2   = c_double()
    rv2   = c_double()
   
    j = lib.iauStarpm(ra1,  dec1, pmr1, pmd1,
                      px1,  rv1,
                      ep1a, ep1b,  ep2a, ep2b, 
                      (ra2),  (dec2),  (pmr2),  (pmd2),
                      (px2),  (rv2))
  

    if   j < 0:
         raise ValueError(iau_pmsafe_msg[j])
    elif j > 0:
         ws.warn(iau_pmsafe_msg[j], UserWarning, 2)
        
    return ra2[0], dec2[0], pmr2[0], pmd2[0], px2[0], rv2[0]



#int iauStarpv(double ra, double dec,
#              double pmr, double pmd, double px, double rv,
#              double pv[2][3])  

iau_starpv_msg = {
                  0: 'no warnings',
                  1: 'distance overridden',
                  2: 'excessive speed',
                  4: 'solution did not converge '
                  }  
 
def pymStarpv(ra, dec, pmr, pmd, px, rv):    
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

    Raises
    ------
    ValueError
        0: 'no warnings',
        1: 'distance overridden',
        2: 'excessive speed',
        4: 'solution did not converge '

    Returns
    -------
    pv : numpy.matrix(2,3)
        pv-vector (au, au/day)

    '''
   
    pv, pv_v  =  vector_double2_ptr()
    
    j = lib.iauStarpv(ra,  dec,  pmr,  pmd,   px,  rv, 
                      pv_v)

    if   j < 0:
         raise ValueError(iau_starpv_msg[j])
    elif j > 0:
         ws.warn(iau_starpv_msg[j], UserWarning, 2)
        
    return pv
