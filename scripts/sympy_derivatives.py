#!/usr/bin/env python
import optparse
import sys
import os
from math import sqrt

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
from sympy import diff

st, dms, sf1, sf2, sfc, f = symbols('st dms sf1 sf2 sfc f')
D = (1 - f) * exp(- dms ** 2 * sf1 ** 2 * st ** 2 / 2) + f * exp(- dms ** 2 * sf2 ** 2 * st ** 2 / 2)
derivs = {}
for s in symbols('st dms sf1 sf2 f'):
    derivs['dD1_d' + s.name] = diff(D, s)

D_sfc = D.subs(sf1, (sfc - f * sf2) / (1 - f))
derivs_sfc = {}
for s in symbols('st dms sf2 sfc f'):
    derivs_sfc['dDc_d' + s.name] = diff(D_sfc, s)

# Use codegen and autowrap to write c code and wrappers
from sympy.utilities.autowrap import CythonCodeWrapper
from sympy.utilities.codegen import CCodeGen
from sympy.utilities.codegen import Routine

routines = []
for ds, args in [(derivs, symbols('st dms sf1 sf2 f')), (derivs_sfc, symbols('st dms sfc sf2 f'))]:
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
