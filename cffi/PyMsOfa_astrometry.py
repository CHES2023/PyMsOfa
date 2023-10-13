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
iau_taiutc_msg = {
                  1: 'dubious year',
                 -1: 'unacceptable date',
                  }  

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


