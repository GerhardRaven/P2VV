import ROOT
from ROOT import (RooDataHist, RooWorkspace,
                  RooBrentRootFinder, RooRealBinding,
                  RooArgSet, RooArgList, RooHistFunc,
                  RooFit)
from ROOT import (TH1F, TFile, TCanvas, Double)

w = RooWorkspace('w')

w.factory('t[-3., 14.]')
w.factory("expr::eff_shape('(t > 0.) ? (1 / (1 + (a * t) ** (-c))) : 0.0001', t, a[1.45], c[2.37])")
eff = w.function('eff_shape')
t = w.var('t')

# Bind the expression to a function and create the root finder object
binding = RooRealBinding(eff, RooArgSet(t))
root_finder = RooBrentRootFinder(binding)

# Some parameter values
t_min = t.getMin()
t_max = t.getMax()
bin_max = 1.
dy = 0.025

# The bin boundaries
from array import array
boundaries = array('f')

boundaries.append(t_min)

t.setVal(t_max)
eff_max = eff.getVal()

# Loop in the vertical direction and use the root finder to create bin boundaries.
dt = boundaries[0]
i = 1
while dt < t_max:
    low = boundaries[i - 1]
    high = low + bin_max
    t.setVal(low)
    val = eff.getVal() + dy
    result = Double(0.)
    r = root_finder.findRoot(result, low, high, val)
    if r:
        boundaries.append(result)
    else:
        boundaries.append(high)
    dt = boundaries[i]
    #print '{:>5d} {: 7.5f} {: 7.5f} {: 7.5f} {: b} {: 7.5f} {: 7.5f}'.format(i, low, high, val, r, result, dt)
    i += 1
boundaries[-1] = t_max

# Create a histogram from the bin boundaries and use the function to set its content
_hist = TH1F('eff_hist', 'eff_hist', len(boundaries) - 1, boundaries)
for i, low in enumerate(boundaries):
    high = boundaries[i + 1]
    mid = (low + high) / 2.
    t.setVal(mid)
    val = eff.getVal()
    _hist.SetBinContent(i + 1, val)
    if i == len(boundaries) - 2:
        break
_hist.SetEntries(len(boundaries) - 1)

# Create a RooHistFunc
eff_hist = RooDataHist('eff_hist', 'eff_hist', RooArgList(t), RooFit.Import(_hist, False))
eff_func = RooHistFunc('acceptance', 'time acceptance', RooArgSet(t), eff_hist)
frame = t.frame()
canvas = TCanvas('canvas', 'canvas', 500, 500)
canvas.cd(1)
eff.plotOn(frame, RooFit.LineColor(ROOT.kBlack))
eff_func.plotOn(frame)
frame.Draw()

_import = getattr(w, 'import')
_import(eff_hist)

root_file = TFile.Open('LP_binning.root', 'recreate')
root_file.WriteTObject(w, 'w')
root_file.Close()
