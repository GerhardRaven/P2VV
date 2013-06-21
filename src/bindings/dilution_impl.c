/******************************************************************************
 *                      Code generated with sympy 0.7.2                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                        This file is part of 'P2VV'                         *
 ******************************************************************************/
#include "dilution.h"
#include <math.h>

double dD1_ddms(double st, double dms, double sf1, double f, double sf2) {

   return -dms*f*pow(sf2, 2)*pow(st, 2)*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2) - dms*pow(sf1, 2)*pow(st, 2)*(-f + 1)*exp(-pow(dms, 2)*pow(sf1, 2)*pow(st, 2)/2);

}

double dD1_df(double st, double dms, double sf1, double f, double sf2) {

   return exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2) - exp(-pow(dms, 2)*pow(sf1, 2)*pow(st, 2)/2);

}

double dD1_dst(double st, double dms, double sf1, double f, double sf2) {

   return -pow(dms, 2)*f*pow(sf2, 2)*st*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2) - pow(dms, 2)*pow(sf1, 2)*st*(-f + 1)*exp(-pow(dms, 2)*pow(sf1, 2)*pow(st, 2)/2);

}

double dD1_dsf1(double st, double dms, double sf1, double f, double sf2) {

   return -pow(dms, 2)*sf1*pow(st, 2)*(-f + 1)*exp(-pow(dms, 2)*pow(sf1, 2)*pow(st, 2)/2);

}

double dD1_dsf2(double st, double dms, double sf1, double f, double sf2) {

   return -pow(dms, 2)*f*sf2*pow(st, 2)*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2);

}

double dDc_dst(double st, double dms, double sfc, double f, double sf2) {

   return -pow(dms, 2)*f*pow(sf2, 2)*st*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2) - pow(dms, 2)*st*pow(-f*sf2 + sfc, 2)*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*sf2 + sfc, 2)/(2*pow(-f + 1, 2)))/(-f + 1);

}

double dDc_df(double st, double dms, double sfc, double f, double sf2) {

   return (-f + 1)*(pow(dms, 2)*sf2*pow(st, 2)*(-f*sf2 + sfc)/pow(-f + 1, 2) - pow(dms, 2)*pow(st, 2)*pow(-f*sf2 + sfc, 2)/pow(-f + 1, 3))*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*sf2 + sfc, 2)/(2*pow(-f + 1, 2))) - exp(-pow(dms, 2)*pow(st, 2)*pow(-f*sf2 + sfc, 2)/(2*pow(-f + 1, 2))) + exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2);

}

double dDc_dsfc(double st, double dms, double sfc, double f, double sf2) {

   return -pow(dms, 2)*pow(st, 2)*(-2*f*sf2 + 2*sfc)*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*sf2 + sfc, 2)/(2*pow(-f + 1, 2)))/(2*(-f + 1));

}

double dDc_dsf2(double st, double dms, double sfc, double f, double sf2) {

   return -pow(dms, 2)*f*sf2*pow(st, 2)*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2) + pow(dms, 2)*f*pow(st, 2)*(-f*sf2 + sfc)*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*sf2 + sfc, 2)/(2*pow(-f + 1, 2)))/(-f + 1);

}

double dDc_ddms(double st, double dms, double sfc, double f, double sf2) {

   return -dms*f*pow(sf2, 2)*pow(st, 2)*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2) - dms*pow(st, 2)*pow(-f*sf2 + sfc, 2)*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*sf2 + sfc, 2)/(2*pow(-f + 1, 2)))/(-f + 1);

}
