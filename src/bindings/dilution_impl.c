/******************************************************************************
 *                      Code generated with sympy 0.7.2                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                        This file is part of 'P2VV'                         *
 ******************************************************************************/
#include "dilution.h"
#include <math.h>

double dDs2_df(double st, double dms, double sf1, double f, double sf2) {

   return (f*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2) + (-f + 1)*exp(-pow(dms, 2)*pow(sf1, 2)*pow(st, 2)/2))*(2*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2) - 2*exp(-pow(dms, 2)*pow(sf1, 2)*pow(st, 2)/2));

}

double dDs2_ddms(double st, double dms, double sf1, double f, double sf2) {

   return (f*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2) + (-f + 1)*exp(-pow(dms, 2)*pow(sf1, 2)*pow(st, 2)/2))*(-2*dms*f*pow(sf2, 2)*pow(st, 2)*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2) - 2*dms*pow(sf1, 2)*pow(st, 2)*(-f + 1)*exp(-pow(dms, 2)*pow(sf1, 2)*pow(st, 2)/2));

}

double dDs2_dsf2(double st, double dms, double sf1, double f, double sf2) {

   return -2*pow(dms, 2)*f*sf2*pow(st, 2)*(f*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2) + (-f + 1)*exp(-pow(dms, 2)*pow(sf1, 2)*pow(st, 2)/2))*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2);

}

double dDs2_dsf1(double st, double dms, double sf1, double f, double sf2) {

   return -2*pow(dms, 2)*sf1*pow(st, 2)*(-f + 1)*(f*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2) + (-f + 1)*exp(-pow(dms, 2)*pow(sf1, 2)*pow(st, 2)/2))*exp(-pow(dms, 2)*pow(sf1, 2)*pow(st, 2)/2);

}

double dDs2_dst(double st, double dms, double sf1, double f, double sf2) {

   return (f*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2) + (-f + 1)*exp(-pow(dms, 2)*pow(sf1, 2)*pow(st, 2)/2))*(-2*pow(dms, 2)*f*pow(sf2, 2)*st*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2) - 2*pow(dms, 2)*pow(sf1, 2)*st*(-f + 1)*exp(-pow(dms, 2)*pow(sf1, 2)*pow(st, 2)/2));

}

double dDc2_df(double st, double dms, double sfc, double f, double sf2) {

   return (f*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2) + (-f + 1)*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*sf2 + sfc, 2)/(2*pow(-f + 1, 2))))*(2*(-f + 1)*(pow(dms, 2)*sf2*pow(st, 2)*(-f*sf2 + sfc)/pow(-f + 1, 2) - pow(dms, 2)*pow(st, 2)*pow(-f*sf2 + sfc, 2)/pow(-f + 1, 3))*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*sf2 + sfc, 2)/(2*pow(-f + 1, 2))) - 2*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*sf2 + sfc, 2)/(2*pow(-f + 1, 2))) + 2*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2));

}

double dDc2_ddms(double st, double dms, double sfc, double f, double sf2) {

   return (f*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2) + (-f + 1)*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*sf2 + sfc, 2)/(2*pow(-f + 1, 2))))*(-2*dms*f*pow(sf2, 2)*pow(st, 2)*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2) - 2*dms*pow(st, 2)*pow(-f*sf2 + sfc, 2)*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*sf2 + sfc, 2)/(2*pow(-f + 1, 2)))/(-f + 1));

}

double dDc2_dsf2(double st, double dms, double sfc, double f, double sf2) {

   return (f*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2) + (-f + 1)*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*sf2 + sfc, 2)/(2*pow(-f + 1, 2))))*(-2*pow(dms, 2)*f*sf2*pow(st, 2)*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2) + 2*pow(dms, 2)*f*pow(st, 2)*(-f*sf2 + sfc)*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*sf2 + sfc, 2)/(2*pow(-f + 1, 2)))/(-f + 1));

}

double dDc2_dsfc(double st, double dms, double sfc, double f, double sf2) {

   return -pow(dms, 2)*pow(st, 2)*(-2*f*sf2 + 2*sfc)*(f*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2) + (-f + 1)*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*sf2 + sfc, 2)/(2*pow(-f + 1, 2))))*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*sf2 + sfc, 2)/(2*pow(-f + 1, 2)))/(-f + 1);

}

double dDc2_dst(double st, double dms, double sfc, double f, double sf2) {

   return (f*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2) + (-f + 1)*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*sf2 + sfc, 2)/(2*pow(-f + 1, 2))))*(-2*pow(dms, 2)*f*pow(sf2, 2)*st*exp(-pow(dms, 2)*pow(sf2, 2)*pow(st, 2)/2) - 2*pow(dms, 2)*st*pow(-f*sf2 + sfc, 2)*exp(-pow(dms, 2)*pow(st, 2)*pow(-f*sf2 + sfc, 2)/(2*pow(-f + 1, 2)))/(-f + 1));

}
