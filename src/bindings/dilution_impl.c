/******************************************************************************
 *                      Code generated with sympy 0.7.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                        This file is part of 'P2VV'                         *
 ******************************************************************************/
#include "dilution.h"
#include <math.h>

double dD2_ddms(double st, double dms, double sf1, double f, double sf2) {

   return (f*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2) + (-f + 1)*exp(-pow(dms, 2)*pow(sf1, 2)*pow(st, 2)/2))*(-2*dms*f*pow(sf2, 2)*pow(st, 2)*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2) - 2*dms*pow(sf1, 2)*pow(st, 2)*(-f + 1)*exp(-pow(dms, 2)*pow(sf1, 2)*pow(st, 2)/2));

}

double dD2_dsf1(double st, double dms, double sf1, double f, double sf2) {

   return -2*pow(dms, 2)*sf1*pow(st, 2)*(-f + 1)*(f*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2) + (-f + 1)*exp(-pow(dms, 2)*pow(sf1, 2)*pow(st, 2)/2))*exp(-pow(dms, 2)*pow(sf1, 2)*pow(st, 2)/2);

}

double dD2_dsf2(double st, double dms, double sf1, double f, double sf2) {

   return -2*pow(dms, 2)*f*sf2*pow(st, 2)*(f*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2) + (-f + 1)*exp(-pow(dms, 2)*pow(sf1, 2)*pow(st, 2)/2))*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2);

}

double dD2_df(double st, double dms, double sf1, double f, double sf2) {

   return (f*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2) + (-f + 1)*exp(-pow(dms, 2)*pow(sf1, 2)*pow(st, 2)/2))*(2*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2) - 2*exp(-pow(dms, 2)*pow(sf1, 2)*pow(st, 2)/2));

}

double dDc2_df(double st, double dms, double sfc, double f, double sf2) {

   return (f*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2) + (-f + 1)*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*sf2 + sfc, 2)/(2*pow(-f + 1, 2))))*(2*(-f + 1)*(pow(dms, 2)*sf2*pow(st, 2)*(-f*sf2 + sfc)/pow(-f + 1, 2) - pow(dms, 2)*pow(st, 2)*pow(-f*sf2 + sfc, 2)/pow(-f + 1, 3))*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*sf2 + sfc, 2)/(2*pow(-f + 1, 2))) - 2*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*sf2 + sfc, 2)/(2*pow(-f + 1, 2))) + 2*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2));

}

double dDc2_dsf2(double st, double dms, double sfc, double f, double sf2) {

   return (f*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2) + (-f + 1)*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*sf2 + sfc, 2)/(2*pow(-f + 1, 2))))*(-2*pow(dms, 2)*f*sf2*pow(st, 2)*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2) + 2*pow(dms, 2)*f*pow(st, 2)*(-f*sf2 + sfc)*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*sf2 + sfc, 2)/(2*pow(-f + 1, 2)))/(-f + 1));

}

double dDc2_dsfc(double st, double dms, double sfc, double f, double sf2) {

   return -pow(dms, 2)*pow(st, 2)*(-2*f*sf2 + 2*sfc)*(f*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2) + (-f + 1)*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*sf2 + sfc, 2)/(2*pow(-f + 1, 2))))*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*sf2 + sfc, 2)/(2*pow(-f + 1, 2)))/(-f + 1);

}

double dDc2_ddms(double st, double dms, double sfc, double f, double sf2) {

   return (f*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2) + (-f + 1)*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*sf2 + sfc, 2)/(2*pow(-f + 1, 2))))*(-2*dms*f*pow(sf2, 2)*pow(st, 2)*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2) - 2*dms*pow(st, 2)*pow(-f*sf2 + sfc, 2)*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*sf2 + sfc, 2)/(2*pow(-f + 1, 2)))/(-f + 1));

}

double dDs2_dsfs(double st, double dms, double sfc, double f, double sfs) {

   return (f*exp(-pow(dms, 2)*pow(st, 2)*pow(sfc + sfs*sqrt((-f + 1)/f), 2)/2) + (-f + 1)*exp(-pow(dms, 2)*pow(st, 2)*pow(sfc - sfs*sqrt(f/(-f + 1)), 2)/2))*(-2*pow(dms, 2)*f*pow(st, 2)*sqrt((-f + 1)/f)*(sfc + sfs*sqrt((-f + 1)/f))*exp(-pow(dms, 2)*pow(st, 2)*pow(sfc + sfs*sqrt((-f + 1)/f), 2)/2) + 2*pow(dms, 2)*pow(st, 2)*sqrt(f/(-f + 1))*(-f + 1)*(sfc - sfs*sqrt(f/(-f + 1)))*exp(-pow(dms, 2)*pow(st, 2)*pow(sfc - sfs*sqrt(f/(-f + 1)), 2)/2));

}

double dDs2_df(double st, double dms, double sfc, double f, double sfs) {

   return (f*exp(-pow(dms, 2)*pow(st, 2)*pow(sfc + sfs*sqrt((-f + 1)/f), 2)/2) + (-f + 1)*exp(-pow(dms, 2)*pow(st, 2)*pow(sfc - sfs*sqrt(f/(-f + 1)), 2)/2))*(-2*pow(dms, 2)*pow(f, 2)*sfs*pow(st, 2)*sqrt((-f + 1)/f)*(-1/(2*f) - (-f + 1)/(2*pow(f, 2)))*(sfc + sfs*sqrt((-f + 1)/f))*exp(-pow(dms, 2)*pow(st, 2)*pow(sfc + sfs*sqrt((-f + 1)/f), 2)/2)/(-f + 1) + 2*pow(dms, 2)*sfs*pow(st, 2)*sqrt(f/(-f + 1))*pow(-f + 1, 2)*(sfc - sfs*sqrt(f/(-f + 1)))*(f/(2*pow(-f + 1, 2)) + 1/(2*(-f + 1)))*exp(-pow(dms, 2)*pow(st, 2)*pow(sfc - sfs*sqrt(f/(-f + 1)), 2)/2)/f - 2*exp(-pow(dms, 2)*pow(st, 2)*pow(sfc - sfs*sqrt(f/(-f + 1)), 2)/2) + 2*exp(-pow(dms, 2)*pow(st, 2)*pow(sfc + sfs*sqrt((-f + 1)/f), 2)/2));

}

double dDs2_dsfc(double st, double dms, double sfc, double f, double sfs) {

   return (f*exp(-pow(dms, 2)*pow(st, 2)*pow(sfc + sfs*sqrt((-f + 1)/f), 2)/2) + (-f + 1)*exp(-pow(dms, 2)*pow(st, 2)*pow(sfc - sfs*sqrt(f/(-f + 1)), 2)/2))*(-pow(dms, 2)*f*pow(st, 2)*(2*sfc + 2*sfs*sqrt((-f + 1)/f))*exp(-pow(dms, 2)*pow(st, 2)*pow(sfc + sfs*sqrt((-f + 1)/f), 2)/2) - pow(dms, 2)*pow(st, 2)*(-f + 1)*(2*sfc - 2*sfs*sqrt(f/(-f + 1)))*exp(-pow(dms, 2)*pow(st, 2)*pow(sfc - sfs*sqrt(f/(-f + 1)), 2)/2));

}

double dDs2_ddms(double st, double dms, double sfc, double f, double sfs) {

   return (f*exp(-pow(dms, 2)*pow(st, 2)*pow(sfc + sfs*sqrt((-f + 1)/f), 2)/2) + (-f + 1)*exp(-pow(dms, 2)*pow(st, 2)*pow(sfc - sfs*sqrt(f/(-f + 1)), 2)/2))*(-2*dms*f*pow(st, 2)*pow(sfc + sfs*sqrt((-f + 1)/f), 2)*exp(-pow(dms, 2)*pow(st, 2)*pow(sfc + sfs*sqrt((-f + 1)/f), 2)/2) - 2*dms*pow(st, 2)*(-f + 1)*pow(sfc - sfs*sqrt(f/(-f + 1)), 2)*exp(-pow(dms, 2)*pow(st, 2)*pow(sfc - sfs*sqrt(f/(-f + 1)), 2)/2));

}

double dDcc2_dsf2o(double st, double stm, double dms, double sfco, double sfcs, double f, double sf2o, double sf2s) {

   return (f*exp(-pow(dms, 2)*pow(st, 2)*pow(sf2o + sf2s*(st - stm), 2)/2) + (-f + 1)*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*(sf2o + sf2s*(st - stm)) + sfco + sfcs*(st - stm), 2)/(2*pow(-f + 1, 2))))*(-pow(dms, 2)*f*pow(st, 2)*(2*sf2o + 2*sf2s*(st - stm))*exp(-pow(dms, 2)*pow(st, 2)*pow(sf2o + sf2s*(st - stm), 2)/2) + 2*pow(dms, 2)*f*pow(st, 2)*(-f*(sf2o + sf2s*(st - stm)) + sfco + sfcs*(st - stm))*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*(sf2o + sf2s*(st - stm)) + sfco + sfcs*(st - stm), 2)/(2*pow(-f + 1, 2)))/(-f + 1));

}

double dDcc2_ddms(double st, double stm, double dms, double sfco, double sfcs, double f, double sf2o, double sf2s) {

   return (f*exp(-pow(dms, 2)*pow(st, 2)*pow(sf2o + sf2s*(st - stm), 2)/2) + (-f + 1)*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*(sf2o + sf2s*(st - stm)) + sfco + sfcs*(st - stm), 2)/(2*pow(-f + 1, 2))))*(-2*dms*f*pow(st, 2)*pow(sf2o + sf2s*(st - stm), 2)*exp(-pow(dms, 2)*pow(st, 2)*pow(sf2o + sf2s*(st - stm), 2)/2) - 2*dms*pow(st, 2)*pow(-f*(sf2o + sf2s*(st - stm)) + sfco + sfcs*(st - stm), 2)*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*(sf2o + sf2s*(st - stm)) + sfco + sfcs*(st - stm), 2)/(2*pow(-f + 1, 2)))/(-f + 1));

}

double dDcc2_dsfcs(double st, double stm, double dms, double sfco, double sfcs, double f, double sf2o, double sf2s) {

   return -pow(dms, 2)*pow(st, 2)*(2*st - 2*stm)*(f*exp(-pow(dms, 2)*pow(st, 2)*pow(sf2o + sf2s*(st - stm), 2)/2) + (-f + 1)*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*(sf2o + sf2s*(st - stm)) + sfco + sfcs*(st - stm), 2)/(2*pow(-f + 1, 2))))*(-f*(sf2o + sf2s*(st - stm)) + sfco + sfcs*(st - stm))*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*(sf2o + sf2s*(st - stm)) + sfco + sfcs*(st - stm), 2)/(2*pow(-f + 1, 2)))/(-f + 1);

}

double dDcc2_dsfco(double st, double stm, double dms, double sfco, double sfcs, double f, double sf2o, double sf2s) {

   return -pow(dms, 2)*pow(st, 2)*(f*exp(-pow(dms, 2)*pow(st, 2)*pow(sf2o + sf2s*(st - stm), 2)/2) + (-f + 1)*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*(sf2o + sf2s*(st - stm)) + sfco + sfcs*(st - stm), 2)/(2*pow(-f + 1, 2))))*(-2*f*(sf2o + sf2s*(st - stm)) + 2*sfco + 2*sfcs*(st - stm))*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*(sf2o + sf2s*(st - stm)) + sfco + sfcs*(st - stm), 2)/(2*pow(-f + 1, 2)))/(-f + 1);

}

double dDcc2_df(double st, double stm, double dms, double sfco, double sfcs, double f, double sf2o, double sf2s) {

   return (f*exp(-pow(dms, 2)*pow(st, 2)*pow(sf2o + sf2s*(st - stm), 2)/2) + (-f + 1)*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*(sf2o + sf2s*(st - stm)) + sfco + sfcs*(st - stm), 2)/(2*pow(-f + 1, 2))))*(2*(-f + 1)*(-pow(dms, 2)*pow(st, 2)*(-2*sf2o - 2*sf2s*(st - stm))*(-f*(sf2o + sf2s*(st - stm)) + sfco + sfcs*(st - stm))/(2*pow(-f + 1, 2)) - pow(dms, 2)*pow(st, 2)*pow(-f*(sf2o + sf2s*(st - stm)) + sfco + sfcs*(st - stm), 2)/pow(-f + 1, 3))*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*(sf2o + sf2s*(st - stm)) + sfco + sfcs*(st - stm), 2)/(2*pow(-f + 1, 2))) - 2*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*(sf2o + sf2s*(st - stm)) + sfco + sfcs*(st - stm), 2)/(2*pow(-f + 1, 2))) + 2*exp(-pow(dms, 2)*pow(st, 2)*pow(sf2o + sf2s*(st - stm), 2)/2));

}

double dDcc2_dsf2s(double st, double stm, double dms, double sfco, double sfcs, double f, double sf2o, double sf2s) {

   return (f*exp(-pow(dms, 2)*pow(st, 2)*pow(sf2o + sf2s*(st - stm), 2)/2) + (-f + 1)*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*(sf2o + sf2s*(st - stm)) + sfco + sfcs*(st - stm), 2)/(2*pow(-f + 1, 2))))*(-pow(dms, 2)*f*pow(st, 2)*(sf2o + sf2s*(st - stm))*(2*st - 2*stm)*exp(-pow(dms, 2)*pow(st, 2)*pow(sf2o + sf2s*(st - stm), 2)/2) + 2*pow(dms, 2)*f*pow(st, 2)*(st - stm)*(-f*(sf2o + sf2s*(st - stm)) + sfco + sfcs*(st - stm))*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*(sf2o + sf2s*(st - stm)) + sfco + sfcs*(st - stm), 2)/(2*pow(-f + 1, 2)))/(-f + 1));

}

double dDsc2_dsfco(double st, double stm, double dms, double sfco, double sfcs, double f, double sfso, double sfss) {

   return (f*exp(-pow(dms, 2)*pow(st, 2)*pow(sfco + sfcs*(st - stm) + sqrt((-f + 1)/f)*(sfso + sfss*(st - stm)), 2)/2) + (-f + 1)*exp(-pow(dms, 2)*pow(st, 2)*pow(sfco + sfcs*(st - stm) - sqrt(f/(-f + 1))*(sfso + sfss*(st - stm)), 2)/2))*(-pow(dms, 2)*f*pow(st, 2)*(2*sfco + 2*sfcs*(st - stm) + 2*sqrt((-f + 1)/f)*(sfso + sfss*(st - stm)))*exp(-pow(dms, 2)*pow(st, 2)*pow(sfco + sfcs*(st - stm) + sqrt((-f + 1)/f)*(sfso + sfss*(st - stm)), 2)/2) - pow(dms, 2)*pow(st, 2)*(-f + 1)*(2*sfco + 2*sfcs*(st - stm) - 2*sqrt(f/(-f + 1))*(sfso + sfss*(st - stm)))*exp(-pow(dms, 2)*pow(st, 2)*pow(sfco + sfcs*(st - stm) - sqrt(f/(-f + 1))*(sfso + sfss*(st - stm)), 2)/2));

}

double dDsc2_df(double st, double stm, double dms, double sfco, double sfcs, double f, double sfso, double sfss) {

   return (f*exp(-pow(dms, 2)*pow(st, 2)*pow(sfco + sfcs*(st - stm) + sqrt((-f + 1)/f)*(sfso + sfss*(st - stm)), 2)/2) + (-f + 1)*exp(-pow(dms, 2)*pow(st, 2)*pow(sfco + sfcs*(st - stm) - sqrt(f/(-f + 1))*(sfso + sfss*(st - stm)), 2)/2))*(-2*pow(dms, 2)*pow(f, 2)*pow(st, 2)*sqrt((-f + 1)/f)*(-1/(2*f) - (-f + 1)/(2*pow(f, 2)))*(sfso + sfss*(st - stm))*(sfco + sfcs*(st - stm) + sqrt((-f + 1)/f)*(sfso + sfss*(st - stm)))*exp(-pow(dms, 2)*pow(st, 2)*pow(sfco + sfcs*(st - stm) + sqrt((-f + 1)/f)*(sfso + sfss*(st - stm)), 2)/2)/(-f + 1) + 2*pow(dms, 2)*pow(st, 2)*sqrt(f/(-f + 1))*pow(-f + 1, 2)*(sfso + sfss*(st - stm))*(f/(2*pow(-f + 1, 2)) + 1/(2*(-f + 1)))*(sfco + sfcs*(st - stm) - sqrt(f/(-f + 1))*(sfso + sfss*(st - stm)))*exp(-pow(dms, 2)*pow(st, 2)*pow(sfco + sfcs*(st - stm) - sqrt(f/(-f + 1))*(sfso + sfss*(st - stm)), 2)/2)/f - 2*exp(-pow(dms, 2)*pow(st, 2)*pow(sfco + sfcs*(st - stm) - sqrt(f/(-f + 1))*(sfso + sfss*(st - stm)), 2)/2) + 2*exp(-pow(dms, 2)*pow(st, 2)*pow(sfco + sfcs*(st - stm) + sqrt((-f + 1)/f)*(sfso + sfss*(st - stm)), 2)/2));

}

double dDsc2_dsfss(double st, double stm, double dms, double sfco, double sfcs, double f, double sfso, double sfss) {

   return (f*exp(-pow(dms, 2)*pow(st, 2)*pow(sfco + sfcs*(st - stm) + sqrt((-f + 1)/f)*(sfso + sfss*(st - stm)), 2)/2) + (-f + 1)*exp(-pow(dms, 2)*pow(st, 2)*pow(sfco + sfcs*(st - stm) - sqrt(f/(-f + 1))*(sfso + sfss*(st - stm)), 2)/2))*(-2*pow(dms, 2)*f*pow(st, 2)*sqrt((-f + 1)/f)*(st - stm)*(sfco + sfcs*(st - stm) + sqrt((-f + 1)/f)*(sfso + sfss*(st - stm)))*exp(-pow(dms, 2)*pow(st, 2)*pow(sfco + sfcs*(st - stm) + sqrt((-f + 1)/f)*(sfso + sfss*(st - stm)), 2)/2) + 2*pow(dms, 2)*pow(st, 2)*sqrt(f/(-f + 1))*(-f + 1)*(st - stm)*(sfco + sfcs*(st - stm) - sqrt(f/(-f + 1))*(sfso + sfss*(st - stm)))*exp(-pow(dms, 2)*pow(st, 2)*pow(sfco + sfcs*(st - stm) - sqrt(f/(-f + 1))*(sfso + sfss*(st - stm)), 2)/2));

}

double dDsc2_dsfso(double st, double stm, double dms, double sfco, double sfcs, double f, double sfso, double sfss) {

   return (f*exp(-pow(dms, 2)*pow(st, 2)*pow(sfco + sfcs*(st - stm) + sqrt((-f + 1)/f)*(sfso + sfss*(st - stm)), 2)/2) + (-f + 1)*exp(-pow(dms, 2)*pow(st, 2)*pow(sfco + sfcs*(st - stm) - sqrt(f/(-f + 1))*(sfso + sfss*(st - stm)), 2)/2))*(-2*pow(dms, 2)*f*pow(st, 2)*sqrt((-f + 1)/f)*(sfco + sfcs*(st - stm) + sqrt((-f + 1)/f)*(sfso + sfss*(st - stm)))*exp(-pow(dms, 2)*pow(st, 2)*pow(sfco + sfcs*(st - stm) + sqrt((-f + 1)/f)*(sfso + sfss*(st - stm)), 2)/2) + 2*pow(dms, 2)*pow(st, 2)*sqrt(f/(-f + 1))*(-f + 1)*(sfco + sfcs*(st - stm) - sqrt(f/(-f + 1))*(sfso + sfss*(st - stm)))*exp(-pow(dms, 2)*pow(st, 2)*pow(sfco + sfcs*(st - stm) - sqrt(f/(-f + 1))*(sfso + sfss*(st - stm)), 2)/2));

}

double dDsc2_ddms(double st, double stm, double dms, double sfco, double sfcs, double f, double sfso, double sfss) {

   return (f*exp(-pow(dms, 2)*pow(st, 2)*pow(sfco + sfcs*(st - stm) + sqrt((-f + 1)/f)*(sfso + sfss*(st - stm)), 2)/2) + (-f + 1)*exp(-pow(dms, 2)*pow(st, 2)*pow(sfco + sfcs*(st - stm) - sqrt(f/(-f + 1))*(sfso + sfss*(st - stm)), 2)/2))*(-2*dms*f*pow(st, 2)*pow(sfco + sfcs*(st - stm) + sqrt((-f + 1)/f)*(sfso + sfss*(st - stm)), 2)*exp(-pow(dms, 2)*pow(st, 2)*pow(sfco + sfcs*(st - stm) + sqrt((-f + 1)/f)*(sfso + sfss*(st - stm)), 2)/2) - 2*dms*pow(st, 2)*(-f + 1)*pow(sfco + sfcs*(st - stm) - sqrt(f/(-f + 1))*(sfso + sfss*(st - stm)), 2)*exp(-pow(dms, 2)*pow(st, 2)*pow(sfco + sfcs*(st - stm) - sqrt(f/(-f + 1))*(sfso + sfss*(st - stm)), 2)/2));

}

double dDsc2_dsfcs(double st, double stm, double dms, double sfco, double sfcs, double f, double sfso, double sfss) {

   return (f*exp(-pow(dms, 2)*pow(st, 2)*pow(sfco + sfcs*(st - stm) + sqrt((-f + 1)/f)*(sfso + sfss*(st - stm)), 2)/2) + (-f + 1)*exp(-pow(dms, 2)*pow(st, 2)*pow(sfco + sfcs*(st - stm) - sqrt(f/(-f + 1))*(sfso + sfss*(st - stm)), 2)/2))*(-pow(dms, 2)*f*pow(st, 2)*(2*st - 2*stm)*(sfco + sfcs*(st - stm) + sqrt((-f + 1)/f)*(sfso + sfss*(st - stm)))*exp(-pow(dms, 2)*pow(st, 2)*pow(sfco + sfcs*(st - stm) + sqrt((-f + 1)/f)*(sfso + sfss*(st - stm)), 2)/2) - pow(dms, 2)*pow(st, 2)*(-f + 1)*(2*st - 2*stm)*(sfco + sfcs*(st - stm) - sqrt(f/(-f + 1))*(sfso + sfss*(st - stm)))*exp(-pow(dms, 2)*pow(st, 2)*pow(sfco + sfcs*(st - stm) - sqrt(f/(-f + 1))*(sfso + sfss*(st - stm)), 2)/2));

}
