import sys
import os
from ToyMCUtils import Toy

toy = Toy()
parser = toy.parser()
(options, args) = toy.configure()
if len(args) != 1:
    print parser.usage
    sys.exit(-2)

dy = None
try:
    dy = float(args[0])
except ValueError:
    print 'Bad dy size given, must be convertible to float'
    sys.exit(-2)
from ROOT import (RooDataHist, RooWorkspace,
                  RooBrentRootFinder, RooRealBinding,
                  RooArgSet, RooArgList, RooHistFunc,
                  RooFit)
from ROOT import (TH1F, TFile, TCanvas, Double)
from P2VVLoad import P2VVLibrary
from RooFitWrappers import *

# Create a HistFunc
obj = RooObject(workspace = 'w')
w = obj.ws()

t = RealVar('t', Title = 'decay time', Unit='ps', Observable = True,  MinMax=(-3,14))
a = RealVar('a', Title = 'a', Value = 1.45, MinMax = (1, 2))
c = RealVar('c', Title = 'c', Value = -2.37, MinMax = (-3, 2))

eff = FormulaVar('eff_shape', "(t > 0.) ? (1 / (1 + (a * t) ** (c))) : 0.0001", [t, a, c])

# Bind the expression to a function and create the root finder object
binding = RooRealBinding(eff._target_(), RooArgSet(t._target_()))
root_finder = RooBrentRootFinder(binding)

# Some parameter values
t_min = t.getMin()
t_max = t.getMax()
bin_max = 1.

# The bin boundaries
from array import array
from ROOT import vector
boundaries = array('d')

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
    print '{:>5d} {: 7.5f} {: 7.5f} {: 7.5f} {: b} {: 7.5f} {: 7.5f}'.format(i, low, high, val, r, result, dt)
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
    print '{:>5d} {: 7.5f} {: 7.5f} {: 7.5f} {: 7.5f}'.format(i, low, mid, high, val)
    if i == len(boundaries) - 2:
        break
_hist.SetEntries(len(boundaries) - 1)

from ROOT import RooBinning
binning = RooBinning(len(boundaries) - 1, boundaries, "binning")

## Use EffHistProd to generate events
t_flat = UniformPdf('t_flat', Arguments = [t])
eff_func = HistFunc('t_acceptance', Histogram = _hist, Observables = [t])
pdf = eff_func * t_flat

from ROOT import RooBinnedPdf
binned_pdf = RooBinnedPdf('binned_pdf', 'binned_pdf', t._target_(), "binning", eff._target_())

toy.run(Observables = w.argSet('t'), Pdf = binned_pdf, GenPdf = pdf)

toy.write_output()
