cdef extern from "dilution.h":
   double dDs2_df(double st, double dms, double sf1, double f, double sf2)

def dDs2_df_c(double st, double dms, double sf1, double f, double sf2):
   return dDs2_df(st, dms, sf1, f, sf2)

cdef extern from "dilution.h":
   double dDs2_dsf1(double st, double dms, double sf1, double f, double sf2)

def dDs2_dsf1_c(double st, double dms, double sf1, double f, double sf2):
   return dDs2_dsf1(st, dms, sf1, f, sf2)

cdef extern from "dilution.h":
   double dDs2_dsf2(double st, double dms, double sf1, double f, double sf2)

def dDs2_dsf2_c(double st, double dms, double sf1, double f, double sf2):
   return dDs2_dsf2(st, dms, sf1, f, sf2)

cdef extern from "dilution.h":
   double dDs2_ddms(double st, double dms, double sf1, double f, double sf2)

def dDs2_ddms_c(double st, double dms, double sf1, double f, double sf2):
   return dDs2_ddms(st, dms, sf1, f, sf2)

cdef extern from "dilution.h":
   double dDc2_df(double st, double dms, double sfc, double f, double sf2)

def dDc2_df_c(double st, double dms, double sfc, double f, double sf2):
   return dDc2_df(st, dms, sfc, f, sf2)

cdef extern from "dilution.h":
   double dDc2_dsf2(double st, double dms, double sfc, double f, double sf2)

def dDc2_dsf2_c(double st, double dms, double sfc, double f, double sf2):
   return dDc2_dsf2(st, dms, sfc, f, sf2)

cdef extern from "dilution.h":
   double dDc2_dsfc(double st, double dms, double sfc, double f, double sf2)

def dDc2_dsfc_c(double st, double dms, double sfc, double f, double sf2):
   return dDc2_dsfc(st, dms, sfc, f, sf2)

cdef extern from "dilution.h":
   double dDc2_ddms(double st, double dms, double sfc, double f, double sf2)

def dDc2_ddms_c(double st, double dms, double sfc, double f, double sf2):
   return dDc2_ddms(st, dms, sfc, f, sf2)

cdef extern from "dilution.h":
   double dDcc2_dsf2o(double st, double stm, double dms, double sfco, double sfcs, double f, double sf2o, double sf2s)

def dDcc2_dsf2o_c(double st, double stm, double dms, double sfco, double sfcs, double f, double sf2o, double sf2s):
   return dDcc2_dsf2o(st, stm, dms, sfco, sfcs, f, sf2o, sf2s)

cdef extern from "dilution.h":
   double dDcc2_ddms(double st, double stm, double dms, double sfco, double sfcs, double f, double sf2o, double sf2s)

def dDcc2_ddms_c(double st, double stm, double dms, double sfco, double sfcs, double f, double sf2o, double sf2s):
   return dDcc2_ddms(st, stm, dms, sfco, sfcs, f, sf2o, sf2s)

cdef extern from "dilution.h":
   double dDcc2_dsfcs(double st, double stm, double dms, double sfco, double sfcs, double f, double sf2o, double sf2s)

def dDcc2_dsfcs_c(double st, double stm, double dms, double sfco, double sfcs, double f, double sf2o, double sf2s):
   return dDcc2_dsfcs(st, stm, dms, sfco, sfcs, f, sf2o, sf2s)

cdef extern from "dilution.h":
   double dDcc2_dsfco(double st, double stm, double dms, double sfco, double sfcs, double f, double sf2o, double sf2s)

def dDcc2_dsfco_c(double st, double stm, double dms, double sfco, double sfcs, double f, double sf2o, double sf2s):
   return dDcc2_dsfco(st, stm, dms, sfco, sfcs, f, sf2o, sf2s)

cdef extern from "dilution.h":
   double dDcc2_df(double st, double stm, double dms, double sfco, double sfcs, double f, double sf2o, double sf2s)

def dDcc2_df_c(double st, double stm, double dms, double sfco, double sfcs, double f, double sf2o, double sf2s):
   return dDcc2_df(st, stm, dms, sfco, sfcs, f, sf2o, sf2s)

cdef extern from "dilution.h":
   double dDcc2_dsf2s(double st, double stm, double dms, double sfco, double sfcs, double f, double sf2o, double sf2s)

def dDcc2_dsf2s_c(double st, double stm, double dms, double sfco, double sfcs, double f, double sf2o, double sf2s):
   return dDcc2_dsf2s(st, stm, dms, sfco, sfcs, f, sf2o, sf2s)

