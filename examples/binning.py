#!/usr/bin/env python
import optparse
import sys
import os

parser = optparse.OptionParser(usage = '%prog dy_size')
parser.add_option("-p", '--prefix', dest = "prefix", default = 'LP',
                  action = 'store', help = 'which prefix to use for the output filename')
parser.add_option('--plot', dest = "plot", default = False,
                  action = 'store_true', help = 'plot function and histo')

(options, args) = parser.parse_args()

if len(args) != 1:
    print parser.usage
    sys.exit(-2)

dy = None
try:
    dy = float(args[0])
except ValueError:
    print 'Bad dy size given, must be convertible to float'
    sys.exit(-2)

import ROOT
from ROOT import (RooDataHist, RooWorkspace,
                  RooBrentRootFinder, RooRealBinding,
                  RooArgSet, RooArgList, RooHistFunc,
                  RooFit)
from ROOT import (TH1F, TFile, TCanvas, Double)
from RooFitWrappers import *

# Create a HistFunc
obj = RooObject(workspace = 'w')
w = obj.ws()

t = RealVar('t', Title = 'decay time', Unit='ps', Observable = True,  MinMax=(-3,14)  )
a = ConstVar(Name = 'a', Value = 1.45)
c = ConstVar(Name = 'c', Value = 2.37)
eff = FormulaVar('eff_shape', "(t > 0.) ? (1 / (1 + (%f * t) ** (-%f))) : 0.0001" % (1.45, 2.37), [t])

# Bind the expression to a function and create the root finder object
binding = RooRealBinding(eff._target_(), RooArgSet(t._target_()))
root_finder = RooBrentRootFinder(binding)

# Some parameter values
t_min = t.getMin()
t_max = t.getMax()
bin_max = 1.

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

eff_func = HistFunc('t_acceptance', Histogram = _hist, Observables = [t])

if options.plot:
    frame = t.frame()
    canvas = TCanvas('canvas', 'canvas', 1000, 500)
    canvas.Divide(2, 1)
    canvas.cd(1)
    eff.plotOn(frame, RooFit.LineColor(ROOT.kBlack))
    eff_func.plotOn(frame)
    frame.Draw()
    canvas.cd(2)
    frame = t.frame()    
    frame.Draw()

root_file = TFile.Open('%s.root' % '_'.join((options.prefix, str(dy).replace('.', '_'))), 'recreate')
root_file.WriteTObject(_hist, 'eff_hist')
root_file.Close()
