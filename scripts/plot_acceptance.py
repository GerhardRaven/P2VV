from P2VV.Load import P2VVLibrary
from P2VV.Load import LHCbStyle
from P2VV.Load import RooFitOutput
from P2VV.RooFitWrappers import *

obj = RooObject( workspace = 'w')
w = obj.ws()

t = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax=(0.3, 14))
run_period = Category('runPeriod', States = {'p2011' : 2011, 'p2012' : 2012}, Observable = True)

from ROOT import TFile
input_file = TFile("acceptance_fit_result_uniform_UB.root")
result = input_file.Get("result")
fpf = dict([(p.GetName(), p) for p in result.floatParsFinal()])
for p in fpf.itervalues():
    w.put(p, Silence = True)

hist_file = TFile.Open('/project/bfys/raaij/p2vv/data/Bs_HltPropertimeAcceptance_Data_2011_40bins.root')
hist = hist_file.Get('Bs_HltPropertimeAcceptance_Data_2011_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted')
from array import array
bname = 'acc_binning'
bounds = array('d', (hist.GetBinLowEdge(1+i) for i in range(hist.GetNbinsX() + 1)))
from ROOT import RooBinning
binning = RooBinning(len(bounds) - 1, bounds, bname)
t.setBinning(binning, bname)

__frames = []
def plot_shape(p, o, shape, errorOpts = {}, pdfOpts = {}):
    from operator import itemgetter
    p.cd()
    p.SetLogx(True)
    if type(shape) == TGraphErrors:
        shape.Draw("AP")
        frame = shape
    else:
        frame = o.frame()
        if errorOpts:
            r = errorOpts.pop('result')
            errorPlots = dict([(x, c) for x, c in errorOpts.iteritems() if type(x) == int])
            for x in errorPlots.keys():
                errorOpts.pop(x)
            entries = sorted(errorPlots.iteritems(), key = itemgetter(0))
            entries.reverse()
            for x, colour in entries:
                shape.plotOn(frame, VisualizeError = (r, x), FillColor = colour, **errorOpts)
        shape.plotOn(frame, **pdfOpts)
    frame.GetXaxis().SetTitle('t [ps]')
    n = shape.GetName()
    pos = n.find('hlt')
    title = '#varepsilon_{%s}' % n[pos : pos + 4].upper()
    if '1' in title:
        title += ' [a.u.]'
    frame.GetYaxis().SetTitle(title)
    frame.GetYaxis().SetTitleOffset(1)
    frame.GetXaxis().SetTitleOffset(1.05)
    if type(shape) != TGraphErrors:
        frame.Draw()
    __frames.append(frame)

x = array('d', [(bounds[i] + bounds[i + 1]) / 2 for i in range(len(bounds) - 1)])
ex = array('d', [(bounds[i + 1] - bounds[i]) / 2 for i in range(len(bounds) - 1)])

from ROOT import TGraphErrors
shape_patterns = ['runPeriod_%s_' + p + '_%03d' for p in ['hlt1_excl_biased_exclB_bin', 'hlt2_biased_B_bin']]
shapes = dict([(s.GetName(), []) for s in run_period])
from itertools import product
for s, pat in product(shapes.iterkeys(), shape_patterns):
    l = [fpf[pat % (s, i)] for i in range(1, 41)]
    pos = pat.find('hlt')
    n = '%s_%s' % (s, pat[pos : pos + 4])
    y = array('d', [p.getVal() for p in l])
    ey = array('d', [p.getError() for p in l])
    g = TGraphErrors(len(l), x, y, ex, ey)
    g.SetName(n)
    shapes[s].append(g)
    ## shapes[s].append(BinnedFun(Name = '_'.join((n, 'shape')), ObsVar = t,
    ##                  Binning = binning, Coefficients = l))

eff_canvases = {}
from ROOT import TCanvas, kYellow, kOrange
for period, s in shapes.iteritems():
    n = 'eff_canvas_' + period
    eff_canvas = TCanvas(n, n, 1200, 400)
    eff_canvases[period] = eff_canvas
    for p, shape in zip(eff_canvas.pads(2, 1), sorted(s, key = lambda e: e.GetName())):
        plot_shape(p, t, shape, errorOpts = {'result' : result, 2 : kYellow, 1 : kOrange})
