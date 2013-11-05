cdef extern from "dilution.h":
   double dD2_ddms(double st, double dms, double sf1, double f, double sf2)

def dD2_ddms_c(double st, double dms, double sf1, double f, double sf2):
   return dD2_ddms(st, dms, sf1, f, sf2)

cdef extern from "dilution.h":
   double dD2_dsf1(double st, double dms, double sf1, double f, double sf2)

def dD2_dsf1_c(double st, double dms, double sf1, double f, double sf2):
   return dD2_dsf1(st, dms, sf1, f, sf2)

cdef extern from "dilution.h":
   double dD2_dsf2(double st, double dms, double sf1, double f, double sf2)

def dD2_dsf2_c(double st, double dms, double sf1, double f, double sf2):
   return dD2_dsf2(st, dms, sf1, f, sf2)

cdef extern from "dilution.h":
   double dD2_df(double st, double dms, double sf1, double f, double sf2)

def dD2_df_c(double st, double dms, double sf1, double f, double sf2):
   return dD2_df(st, dms, sf1, f, sf2)

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
   double dDs2_dsfs(double st, double dms, double sfc, double f, double sfs)

def dDs2_dsfs_c(double st, double dms, double sfc, double f, double sfs):
   return dDs2_dsfs(st, dms, sfc, f, sfs)

cdef extern from "dilution.h":
   double dDs2_df(double st, double dms, double sfc, double f, double sfs)

def dDs2_df_c(double st, double dms, double sfc, double f, double sfs):
   return dDs2_df(st, dms, sfc, f, sfs)

cdef extern from "dilution.h":
   double dDs2_dsfc(double st, double dms, double sfc, double f, double sfs)

def dDs2_dsfc_c(double st, double dms, double sfc, double f, double sfs):
   return dDs2_dsfc(st, dms, sfc, f, sfs)

cdef extern from "dilution.h":
   double dDs2_ddms(double st, double dms, double sfc, double f, double sfs)

def dDs2_ddms_c(double st, double dms, double sfc, double f, double sfs):
   return dDs2_ddms(st, dms, sfc, f, sfs)

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

cdef extern from "dilution.h":
   double dDsc2_dsfco(double st, double stm, double dms, double sfco, double sfcs, double f, double sfso, double sfss)

def dDsc2_dsfco_c(double st, double stm, double dms, double sfco, double sfcs, double f, double sfso, double sfss):
   return dDsc2_dsfco(st, stm, dms, sfco, sfcs, f, sfso, sfss)

cdef extern from "dilution.h":
   double dDsc2_df(double st, double stm, double dms, double sfco, double sfcs, double f, double sfso, double sfss)

def dDsc2_df_c(double st, double stm, double dms, double sfco, double sfcs, double f, double sfso, double sfss):
   return dDsc2_df(st, stm, dms, sfco, sfcs, f, sfso, sfss)

cdef extern from "dilution.h":
   double dDsc2_dsfss(double st, double stm, double dms, double sfco, double sfcs, double f, double sfso, double sfss)

def dDsc2_dsfss_c(double st, double stm, double dms, double sfco, double sfcs, double f, double sfso, double sfss):
   return dDsc2_dsfss(st, stm, dms, sfco, sfcs, f, sfso, sfss)

cdef extern from "dilution.h":
   double dDsc2_dsfso(double st, double stm, double dms, double sfco, double sfcs, double f, double sfso, double sfss)

def dDsc2_dsfso_c(double st, double stm, double dms, double sfco, double sfcs, double f, double sfso, double sfss):
   return dDsc2_dsfso(st, stm, dms, sfco, sfcs, f, sfso, sfss)

cdef extern from "dilution.h":
   double dDsc2_ddms(double st, double stm, double dms, double sfco, double sfcs, double f, double sfso, double sfss)

def dDsc2_ddms_c(double st, double stm, double dms, double sfco, double sfcs, double f, double sfso, double sfss):
   return dDsc2_ddms(st, stm, dms, sfco, sfcs, f, sfso, sfss)

cdef extern from "dilution.h":
   double dDsc2_dsfcs(double st, double stm, double dms, double sfco, double sfcs, double f, double sfso, double sfss)

def dDsc2_dsfcs_c(double st, double stm, double dms, double sfco, double sfcs, double f, double sfso, double sfss):
   return dDsc2_dsfcs(st, stm, dms, sfco, sfcs, f, sfso, sfss)

