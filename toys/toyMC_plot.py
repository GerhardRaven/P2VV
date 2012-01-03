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
pull_vars = []
it = gen_params.createIterator()
gp = it.Next()
while gp:
    rp = result_params.find(gp.GetName())
    if not gp:
        print 'Warning: cannot find result parameter for %s' % gp.GetName()
        continue
    # Make error var
    name = gp.GetName() + '_error'
    _error = RooErrorVar(name, name, rp)
    error = data.addColumn(_error)

    # Make pull var
    name = gp.GetName() + '_pull'
    _pull = RooFormulaVar(name, name, '(@0 - %f) / @1' % gp.getVal(), RooArgList(rp, error))
    pull = data.addColumn(_pull)
    pull_vars.append(pull)
    gp = it.Next()

plot_vars = []
for pv in pull_vars:
    mean = data.mean(pv)
    plot = True
    if mean < -5. or mean > 5.:
        plot = False
    sigma = data.sigma(pv)
    if sigma < 0.1 or sigma > 5.:
        plot = False
    if plot:
        print "Plotting %s" % pv.GetName()
        plot_vars.append(pv)
    else:
        print "Not plotting %s" % pv.GetName()
    
from ROOT import TCanvas
hsize = 300
vsize = 300
from math import sqrt
l = len(plot_vars)
nrows = int(sqrt(l))
ncolumns = 1
if nrows == 1:
    ncolumns = l
elif l % nrows == 0:
    ncolumns = l / nrows
else:
    ncolumns = l / nrows + 1
hsize = 300 * ncolumns
vsize = 300 * nrows
if hsize > 1500:
    hsize = int(1500 / ncolumns)
    vsize = hsize
elif vsize > 900:
    vsize = 900 / nrows
    hsize = vsize
canvas = TCanvas('canvas', 'canvas', hsize, vsize)
canvas.Divide(ncolumns, nrows)

from ROOT import (RooRealVar, RooGaussian,
                  SetOwnership, RooFit)
def make_gaus(pv):
    name = pv.GetName()
    mean = RooRealVar('mean', '#mu', 0., -5, 5)
    sigma = RooRealVar('sigma', '#sigma', 1., 0.1, 3)
    pull_pdf = RooGaussian(name + '_gaus', name + '_gaus', pv, mean, sigma)
    SetOwnership(mean, False)
    SetOwnership(sigma, False)
    SetOwnership(pull_pdf, False)
    return pull_pdf

from ROOT import gStyle
gStyle.SetOptTitle(0)

frames = {}
for i, pv in enumerate(plot_vars):
    canvas.cd(i + 1)
    pull_pdf = make_gaus(pv)
    r = pull_pdf.fitTo(data, RooFit.Minimizer('Minuit2'), RooFit.Save(True))
    frame = pv.frame(RooFit.Range(-5, 5))
    data.plotOn(frame)
    if r.status() == 0:
        pull_pdf.plotOn(frame)
        pull_pdf.paramOn(frame)
    frames[pv.GetName()] = frame
    frame.Draw()
canvas.Update()