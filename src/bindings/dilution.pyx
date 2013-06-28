cdef extern from "dilution.h":
   double dDs2_df(double st, double dms, double sf1, double f, double sf2)

def dDs2_df_c(double st, double dms, double sf1, double f, double sf2):
   return dDs2_df(st, dms, sf1, f, sf2)

cdef extern from "dilution.h":
   double dDs2_ddms(double st, double dms, double sf1, double f, double sf2)

def dDs2_ddms_c(double st, double dms, double sf1, double f, double sf2):
   return dDs2_ddms(st, dms, sf1, f, sf2)

cdef extern from "dilution.h":
   double dDs2_dsf2(double st, double dms, double sf1, double f, double sf2)

def dDs2_dsf2_c(double st, double dms, double sf1, double f, double sf2):
   return dDs2_dsf2(st, dms, sf1, f, sf2)

cdef extern from "dilution.h":
   double dDs2_dsf1(double st, double dms, double sf1, double f, double sf2)

def dDs2_dsf1_c(double st, double dms, double sf1, double f, double sf2):
   return dDs2_dsf1(st, dms, sf1, f, sf2)

cdef extern from "dilution.h":
   double dDs2_dst(double st, double dms, double sf1, double f, double sf2)

def dDs2_dst_c(double st, double dms, double sf1, double f, double sf2):
   return dDs2_dst(st, dms, sf1, f, sf2)

cdef extern from "dilution.h":
   double dDc2_df(double st, double dms, double sfc, double f, double sf2)

def dDc2_df_c(double st, double dms, double sfc, double f, double sf2):
   return dDc2_df(st, dms, sfc, f, sf2)

cdef extern from "dilution.h":
   double dDc2_ddms(double st, double dms, double sfc, double f, double sf2)

def dDc2_ddms_c(double st, double dms, double sfc, double f, double sf2):
   return dDc2_ddms(st, dms, sfc, f, sf2)

cdef extern from "dilution.h":
   double dDc2_dsf2(double st, double dms, double sfc, double f, double sf2)

def dDc2_dsf2_c(double st, double dms, double sfc, double f, double sf2):
   return dDc2_dsf2(st, dms, sfc, f, sf2)

cdef extern from "dilution.h":
   double dDc2_dsfc(double st, double dms, double sfc, double f, double sf2)

def dDc2_dsfc_c(double st, double dms, double sfc, double f, double sf2):
   return dDc2_dsfc(st, dms, sfc, f, sf2)

cdef extern from "dilution.h":
   double dDc2_dst(double st, double dms, double sfc, double f, double sf2)

def dDc2_dst_c(double st, double dms, double sfc, double f, double sf2):
   return dDc2_dst(st, dms, sfc, f, sf2)

