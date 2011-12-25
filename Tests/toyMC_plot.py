from optparse import OptionParser
import sys
parser = OptionParser(usage = '%prog input_filename')

parser.add_option("-p", "--param_file", dest = "param_file", default = 'params.root',
                  type = 'string', help = 'file containing generated parameter values')

(options, args) = parser.parse_args()
if len(args) != 1:
    print parser.usage
    sys.exit(-1)
input_filename = args[0]

from ROOT import (RooPullVar, RooConstVar,
                  RooFormulaVar, TFile,
                  TCanvas, RooArgList,
                  RooErrorVar, RooFit,
                  RooArgSet)

param_file = TFile.Open(options.param_file)
if not param_file:
    print 'Error: could not open file with generator parameters: %s' % options.param_file
    sys.exit(-2)
gen_params = param_file.Get('gen_params')

root_file = TFile.Open(input_filename)
if not root_file:
    print 'Error: could not open input file: %s' % input_filename
    sys.exit(-2)

data = root_file.Get('result_data')
result_params = data.get()
error_params = {}
pull_vars = RooArgSet()
it = gen_params.createIterator()
gp = it.Next()
while gp:
    rp = result_params.find(gp.GetName())
    if not gp:
        print 'Warning: cannot find result parameter for %s' % gp.GetName()
        continue
    name = gp.GetName() + '_error'
    _error = RooErrorVar(name, name, rp)
    error = data.addColumn(_error)
    error_params[name] = error
    name = gp.GetName() + '_pull'
    _pull = RooFormulaVar(name, name, '(@0 - @1) / @2', RooArgList(rp, gp, error))
    pull = data.addColumn(_pull)
    pull_vars.add(pull)
    gp = it.Next()
    
from ROOT import TCanvas
hsize = 300
vsize = 300
ncolumns = 5
if pull_vars.getSize() > 15:
    ncolumns = pull_vars.getSize() / 3 if pull_vars.getSize() % 3 == 0 else pull_vars.getSize() / 3 + 1
    hsize = int(1500 / ncolumns)
    vsize = hsize
canvas = TCanvas('canvas', 'canvas', ncolumns * hsize, 3 * vsize)
canvas.Divide(ncolumns, 3)

it = pull_vars.createIterator()
pv = it.Next()
i = 1
while pv:
    canvas.cd(i)
    frame = pv.frame(RooFit.Range(-3, 3))
    data.plotOn(frame)
    frame.Draw()
    pv = it.Next()
    i += 1
