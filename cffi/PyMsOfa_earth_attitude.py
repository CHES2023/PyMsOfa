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
