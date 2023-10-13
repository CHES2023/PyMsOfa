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

