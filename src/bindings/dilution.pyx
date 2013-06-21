cdef extern from "dilution.h":
   double dD1_ddms(double st, double dms, double sf1, double f, double sf2)

def dD1_ddms_c(double st, double dms, double sf1, double f, double sf2):
   return dD1_ddms(st, dms, sf1, f, sf2)

cdef extern from "dilution.h":
   double dD1_df(double st, double dms, double sf1, double f, double sf2)

def dD1_df_c(double st, double dms, double sf1, double f, double sf2):
   return dD1_df(st, dms, sf1, f, sf2)

cdef extern from "dilution.h":
   double dD1_dst(double st, double dms, double sf1, double f, double sf2)

def dD1_dst_c(double st, double dms, double sf1, double f, double sf2):
   return dD1_dst(st, dms, sf1, f, sf2)

cdef extern from "dilution.h":
   double dD1_dsf1(double st, double dms, double sf1, double f, double sf2)

def dD1_dsf1_c(double st, double dms, double sf1, double f, double sf2):
   return dD1_dsf1(st, dms, sf1, f, sf2)

cdef extern from "dilution.h":
   double dD1_dsf2(double st, double dms, double sf1, double f, double sf2)

def dD1_dsf2_c(double st, double dms, double sf1, double f, double sf2):
   return dD1_dsf2(st, dms, sf1, f, sf2)

cdef extern from "dilution.h":
   double dDc_dst(double st, double dms, double sfc, double f, double sf2)

def dDc_dst_c(double st, double dms, double sfc, double f, double sf2):
   return dDc_dst(st, dms, sfc, f, sf2)

cdef extern from "dilution.h":
   double dDc_df(double st, double dms, double sfc, double f, double sf2)

def dDc_df_c(double st, double dms, double sfc, double f, double sf2):
   return dDc_df(st, dms, sfc, f, sf2)

cdef extern from "dilution.h":
   double dDc_dsfc(double st, double dms, double sfc, double f, double sf2)

def dDc_dsfc_c(double st, double dms, double sfc, double f, double sf2):
   return dDc_dsfc(st, dms, sfc, f, sf2)

cdef extern from "dilution.h":
   double dDc_dsf2(double st, double dms, double sfc, double f, double sf2)

def dDc_dsf2_c(double st, double dms, double sfc, double f, double sf2):
   return dDc_dsf2(st, dms, sfc, f, sf2)

cdef extern from "dilution.h":
   double dDc_ddms(double st, double dms, double sfc, double f, double sf2)

def dDc_ddms_c(double st, double dms, double sfc, double f, double sf2):
   return dDc_ddms(st, dms, sfc, f, sf2)

