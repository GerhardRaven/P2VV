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
gen_params = param_file.get('gen_params')

root_file = TFile.Open(input_filename)
if not root_file:
    print 'Error: could not open input file: %s' % input_filename
    sys.exit(-2)

data = root_file.Get('data')
result_params = data.get()
error_params = RooArgSet()
_pull_vars = RooArgList()
for rp in result_params:
    name = rp.GetName() + '_error'
    error = RooErrorVar(name, name, rp)
    name = rp.GetName() + '_pull'
    gp = gen_params.find(rp.GetName())
    if not gp:
        print 'Error: cannot find generator parameter for %s' % rp.GetName()
    _pull = RooFormulaVar(name, name, '(@0 - @1) / @2', RooArgList(rp, error, gp))
    _pull_vars.add(_pull)
pull_vars = data.addColumns(_pull_vars)
    
from ROOT import TCanvas
hsize = 300
vsize = 300
ncolumns = 5
if pull_vars.getSize() > 15:
    ncolumns = pull_vars.getSize / 3 if pull_vars.getSize % 3 == 0 else pull_vars.getSize() / 3 + 1
    hsize = int(1500 / ncolumns)
    vsize = hsize
canvas = TCanvas('canvas', 'canvas', ncolumns * hsize, 3 * vsize)
canvas.Divide(ncolumns, 3)

for p in pull_vars:
    frame = p.frame(RooFit.Range(-3, 3))
    data.plotOn(frame)
    frame.Draw()
