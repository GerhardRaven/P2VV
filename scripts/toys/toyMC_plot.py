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
                  RooArgSet, RooCategory,
                  RooAddPdf)
from P2VV.RooFitDecorators import *

param_file = TFile.Open(options.param_file)
if not param_file:
    print 'Error: could not open file with generator parameters: %s' % options.param_file
    sys.exit(-2)
gen_params = param_file.Get('gen_params')

root_file = TFile.Open(input_filename)
if not root_file:
    print 'Error: could not open input file: %s' % input_filename
    sys.exit(-2)

values = {}

data = root_file.Get('result_data')
result_params = data.get()
pull_vars = []
for gp in gen_params:
    rp = result_params.find(gp.GetName())
    if not rp:
        print 'Warning: cannot find result parameter for %s' % gp.GetName()
        continue
    if isinstance(gp, RooCategory):
        continue
    # Make error var
    name = gp.GetName() + '_error'
    print name
    _error = RooErrorVar(name, name, rp)
    error = data.addColumn(_error)

    # Make pull var
    name = gp.GetName() + '_pull'

    _pull = RooFormulaVar(name, name, '(@0 - %f) / @1' % gp.getVal(), RooArgList(rp, error))
    pull = data.addColumn(_pull)
    pull_vars.append(pull)

plot_vars = []
for gp, pv in zip([p for p in gen_params if result_params.find(p.GetName())], pull_vars):
    plot = True
    rp = result_params.find(gp.GetName())
    if rp.isConstant():
        plot = False
    ## if mean < -5. or mean > 5.:
    ##     plot = False
    ## sigma = data.sigma(pv)
    ## if sigma < 0.1 or sigma > 5.:
    ##     plot = False
    if plot:
        print "Plotting %s" % pv.GetName()
        plot_vars.append(pv)
    else:
        print "Not plotting %s" % pv.GetName()
    
from ROOT import TCanvas
from math import ceil
n_canvas = int(ceil(float(len(plot_vars)) / 9.))
canvases = {'pulls' : [TCanvas('canvas_%i' % i, 'canvas_%i' % i, 900, 900) for i in range(n_canvas)]}
for c in canvases['pulls']:
    c.Divide(3, 3)

from ROOT import (RooRealVar, RooGaussian,
                  SetOwnership, RooFit)
def make_gaus(pv):
    name = pv.GetName()
    mean = RooRealVar('mean', '#mu', 0., -5, 5)
    sigma1 = RooRealVar('sigma1', '#sigma_{1}', 1., 0.1, 3)
    sigma2 = RooRealVar('sigma2', '#sigma_{2}', 4., 0.1, 20)
    frac = RooRealVar('frac_g2', '#frac_g2', 0.1, 0.001, 0.999)
    
    values[pv.GetName()] = (mean, sigma1, sigma2)
    g1 = RooGaussian(name + '_gaus_1', name + '_gaus_1', pv, mean, sigma1)
    g2 = RooGaussian(name + '_gaus_2', name + '_gaus_2', pv, mean, sigma2)
    pull_pdf = RooAddPdf(name + '_pull_pdf', name + '_pull_pdf', g2, g1, frac)
    for ro in [mean, sigma1, sigma2, frac, g1, g2, pull_pdf]:
        SetOwnership(ro, False)
    return pull_pdf

from ROOT import gStyle
gStyle.SetOptTitle(0)

def chunks(l, n):
    """ Yield successive n-sized chunks from l.
    """
    for i in range(0, len(l), n):
        yield l[i:i+n]

frames = {}
results = {}
for canvas, pvs in zip(canvases['pulls'], chunks(plot_vars, 9)):
    for i, pv in enumerate(pvs):
        canvas.cd(i + 1)
        pull_pdf = make_gaus(pv)
        r = pull_pdf.fitTo(data, RooFit.Minimizer('Minuit2'), RooFit.Save(True))
        frame = pv.frame(RooFit.Range(-5, 5))
        data.plotOn(frame)
        if r.status() == 0:
            pull_pdf.plotOn(frame)
            pull_pdf.paramOn(frame)
        frames[pv.GetName()] = frame
        results[pv.GetName()] = r
        frame.Draw()
    canvas.Update()
 
from collections import defaultdict
means = defaultdict(dict)
sigmas = defaultdict(dict)

for k, (m, s1, s2) in values.iteritems():
    if k.find('hlt') == -1:
        continue
    means[k[:4]][k[-12:-5]] = m
    sigmas[k[:4]][k[-12:-5]] = s1

from array import array
from operator import itemgetter
from ROOT import TGraphErrors

__graphs = []
for name, cont in {'means' : means, 'sigmas' : sigmas}.iteritems():
    canvas = TCanvas('canvas_%s' % name, 'canvas_%s' % name, 600, 300)
    canvas.Divide(2, 1)
    canvases[name] = canvas
    for i, (l, bins) in enumerate(cont.iteritems()):
        y = array('d', [v.getVal() for k, v in sorted(bins.iteritems(), key = itemgetter(0))])
        ey = array('d', [v.getError() for k, v in sorted(bins.iteritems(), key = itemgetter(0))])
        x = array('d', [j for j in range(len(y))])
        ex = array('d', [0. for j in range(len(x))])
        graph = TGraphErrors(len(x), x, y, ex, ey)
        graph.SetName('%s_%s' % (l, name))
        graph.SetTitle('%s_%s' % (l, name))
        canvas.cd(i + 1)
        graph.Draw("AP")
        __graphs.append(graph)
