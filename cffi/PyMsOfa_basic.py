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
