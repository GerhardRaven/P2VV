cdef extern from "dilution.h":
   double dD1_ddms(double st, double dms, double sf1, double sf2, double f)

def dD1_ddms_c(double st, double dms, double sf1, double sf2, double f):
   return dD1_ddms(st, dms, sf1, sf2, f)

cdef extern from "dilution.h":
   double dD1_df(double st, double dms, double sf1, double sf2, double f)

def dD1_df_c(double st, double dms, double sf1, double sf2, double f):
   return dD1_df(st, dms, sf1, sf2, f)

cdef extern from "dilution.h":
   double dD1_dst(double st, double dms, double sf1, double sf2, double f)

def dD1_dst_c(double st, double dms, double sf1, double sf2, double f):
   return dD1_dst(st, dms, sf1, sf2, f)

cdef extern from "dilution.h":
   double dD1_dsf1(double st, double dms, double sf1, double sf2, double f)

def dD1_dsf1_c(double st, double dms, double sf1, double sf2, double f):
   return dD1_dsf1(st, dms, sf1, sf2, f)

cdef extern from "dilution.h":
   double dD1_dsf2(double st, double dms, double sf1, double sf2, double f)

def dD1_dsf2_c(double st, double dms, double sf1, double sf2, double f):
   return dD1_dsf2(st, dms, sf1, sf2, f)

cdef extern from "dilution.h":
   double dDc_dst(double st, double dms, double sfc, double sf2, double f)

def dDc_dst_c(double st, double dms, double sfc, double sf2, double f):
   return dDc_dst(st, dms, sfc, sf2, f)

cdef extern from "dilution.h":
   double dDc_df(double st, double dms, double sfc, double sf2, double f)

def dDc_df_c(double st, double dms, double sfc, double sf2, double f):
   return dDc_df(st, dms, sfc, sf2, f)

cdef extern from "dilution.h":
   double dDc_dsfc(double st, double dms, double sfc, double sf2, double f)

def dDc_dsfc_c(double st, double dms, double sfc, double sf2, double f):
   return dDc_dsfc(st, dms, sfc, sf2, f)

cdef extern from "dilution.h":
   double dDc_dsf2(double st, double dms, double sfc, double sf2, double f)

def dDc_dsf2_c(double st, double dms, double sfc, double sf2, double f):
   return dDc_dsf2(st, dms, sfc, sf2, f)

cdef extern from "dilution.h":
   double dDc_ddms(double st, double dms, double sfc, double sf2, double f)

def dDc_ddms_c(double st, double dms, double sfc, double sf2, double f):
   return dDc_ddms(st, dms, sfc, sf2, f)

