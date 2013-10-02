#!/usr/bin/env python
import optparse
import sys
import os

parser = optparse.OptionParser(usage = '%prog output_dir')

(options, args) = parser.parse_args()

if len(args) != 1:
    print parser.usage
    sys.exit(-1)

install_dir = os.path.realpath(args[0])
if not os.path.exists(install_dir):
    os.makedirs(install_dir)

# Import sympy and create the derivatives
from sympy import symbols
from sympy import exp
from sympy import sqrt
from sympy import diff

st, dms, sf1, sf2, f = symbols('st dms sf1 sf2 f')
D = (1 - f) * exp(- dms ** 2 * sf1 ** 2 * st ** 2 / 2) + f * exp(- dms ** 2 * sf2 ** 2 * st ** 2 / 2)
derivs = {}
for s in symbols('dms sf1 sf2 f'):
    derivs['dD2_d' + s.name] = diff(D ** 2, s)

sfc = symbols('sfc')
D_sfc = D.subs(sf1, (sfc - f * sf2) / (1 - f))
derivs_sfc = {}
for s in symbols('dms sf2 sfc f'):
    derivs_sfc['dDc2_d' + s.name] = diff(D_sfc ** 2, s)

sfs = symbols('sfs')
D_sfs = D.subs(sf1, - sqrt(f / (1 - f)) * sfs + sfc)
D_sfs = D_sfs.subs(sf2, sqrt((1 - f) / f) * sfs + sfc)
derivs_sfs = {}
for s in symbols('dms sfc sfs f'):
    derivs_sfs['dDs2_d' + s.name] = diff(D_sfs ** 2, s)

stm, sfco, sfcs, sf2o, sf2s = symbols('stm sfco sfcs sf2o sf2s')
D_sfc_calib = D_sfc.subs([(sfc, sfco + sfcs * (st - stm)), (sf2, sf2o + sf2s * (st - stm))])
derivs_sfc_calib = {}
for s in symbols('dms sfco sfcs sf2o sf2s f'):
    derivs_sfc_calib['dDcc2_d' + s.name] = diff(D_sfc_calib ** 2, s)

stm, sfco, sfcs, sfso, sfss,  = symbols('stm sfco sfcs sfso sfss')
D_sfs_calib = D_sfs.subs([(sfc, sfco + sfcs * (st - stm)), (sfs, sfso + sfss * (st - stm))])
derivs_sfs_calib = {}
for s in symbols('dms sfco sfcs sfso sfss f'):
    derivs_sfs_calib['dDsc2_d' + s.name] = diff(D_sfs_calib ** 2, s)

# Use codegen and autowrap to write c code and wrappers
from sympy.utilities.autowrap import CythonCodeWrapper
from sympy.utilities.codegen import CCodeGen
from sympy.utilities.codegen import Routine

routines = []
for ds, args in [(derivs, symbols('st dms sf1 f sf2')), (derivs_sfc, symbols('st dms sfc f sf2')),
                 (derivs_sfs, symbols('st dms sfc f sfs')),
                 (derivs_sfc_calib, symbols('st stm dms sfco sfcs f sf2o sf2s')),
                 (derivs_sfs_calib, symbols('st stm dms sfco sfcs f sfso sfss'))]:
    for name, expr in ds.iteritems():
        routines.append(Routine(name, expr, argument_sequence = args))

cgen = CCodeGen(project = 'P2VV')
with open(os.path.join(install_dir, "dilution_impl.c"), "w") as c_file:
    cgen.dump_c(routines, c_file, 'dilution')

with open(os.path.join(install_dir, "dilution.h"), "w") as h_file:
    cgen.dump_h(routines, h_file, 'dilution')

wrapper = CythonCodeWrapper(cgen)
with open(os.path.join(install_dir, "dilution.pyx"), "w") as wrapper_file:
    wrapper.dump_pyx(routines, wrapper_file, 'dilution')

# Call cython to create the python wrapper c file.
import subprocess
pwd = os.path.realpath('.')
wd = os.chdir(install_dir)
subprocess.call(['cython', 'dilution.pyx'])
os.chdir(pwd)
