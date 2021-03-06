/******************************************************************************
 *                      Code generated with sympy 0.7.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                        This file is part of 'P2VV'                         *
 ******************************************************************************/


#ifndef P2VV__DILUTION__H
#define P2VV__DILUTION__H

double dD2_ddms(double st, double dms, double sf1, double f, double sf2);
double dD2_dsf1(double st, double dms, double sf1, double f, double sf2);
double dD2_dsf2(double st, double dms, double sf1, double f, double sf2);
double dD2_df(double st, double dms, double sf1, double f, double sf2);
double dDc2_df(double st, double dms, double sfc, double f, double sf2);
double dDc2_dsf2(double st, double dms, double sfc, double f, double sf2);
double dDc2_dsfc(double st, double dms, double sfc, double f, double sf2);
double dDc2_ddms(double st, double dms, double sfc, double f, double sf2);
double dDs2_dsfs(double st, double dms, double sfc, double f, double sfs);
double dDs2_df(double st, double dms, double sfc, double f, double sfs);
double dDs2_dsfc(double st, double dms, double sfc, double f, double sfs);
double dDs2_ddms(double st, double dms, double sfc, double f, double sfs);
double dDcc2_dsf2o(double st, double stm, double dms, double sfco, double sfcs, double f, double sf2o, double sf2s);
double dDcc2_ddms(double st, double stm, double dms, double sfco, double sfcs, double f, double sf2o, double sf2s);
double dDcc2_dsfcs(double st, double stm, double dms, double sfco, double sfcs, double f, double sf2o, double sf2s);
double dDcc2_dsfco(double st, double stm, double dms, double sfco, double sfcs, double f, double sf2o, double sf2s);
double dDcc2_df(double st, double stm, double dms, double sfco, double sfcs, double f, double sf2o, double sf2s);
double dDcc2_dsf2s(double st, double stm, double dms, double sfco, double sfcs, double f, double sf2o, double sf2s);
double dDsc2_dsfco(double st, double stm, double dms, double sfco, double sfcs, double f, double sfso, double sfss);
double dDsc2_df(double st, double stm, double dms, double sfco, double sfcs, double f, double sfso, double sfss);
double dDsc2_dsfss(double st, double stm, double dms, double sfco, double sfcs, double f, double sfso, double sfss);
double dDsc2_dsfso(double st, double stm, double dms, double sfco, double sfcs, double f, double sfso, double sfss);
double dDsc2_ddms(double st, double stm, double dms, double sfco, double sfcs, double f, double sfso, double sfss);
double dDsc2_dsfcs(double st, double stm, double dms, double sfco, double sfcs, double f, double sfso, double sfss);

#endif

