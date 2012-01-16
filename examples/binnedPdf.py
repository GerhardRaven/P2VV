#!/usr/bin/env python
import optparse
import sys
import os

parser = optparse.OptionParser(usage = '%prog dy_size')
parser.add_option("-p", '--prefix', dest = "prefix", default = 'LP',
                  action = 'store', help = 'which prefix to use for the output filename')

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

t = RealVar('t', Title = 'decay time', Unit='ps', Observable = True,  MinMax=(-3,14))
a = RealVar('a', Title = 'a', Value = 1.45, MinMax = (1, 2))
c = RealVar('c', Title = 'c', Value = -2.37, MinMax = (-3, -2))

eff = FormulaVar('eff_shape', "(t > 0.) ? (1 / (1 + (a * t) ** (c))) : 0.0001", [t, a, c])

from P2VVBinningBuilders import build1DVerticalBinning
binning, eff_func = build1DVerticalBinning('binning', eff, t, dy, 1.)

## Use EffHistProd to generate events
t_flat = UniformPdf('t_flat', Arguments = [t])
acceptance = BinnedPdf(Name = 'time_acceptance', Observables = [t], Function = eff,
                       Binning = binning)
pdf = acceptance * t_flat
pdf.Print('t')

data = pdf.generate([t], 100000)

result = pdf.fitTo(data, Save = True, Minimizer = 'Minuit2')

frame = t.frame()
canvas = TCanvas('canvas', 'canvas', 1000, 500)
canvas.Divide(2, 1)
canvas.cd(1)
data.plotOn(frame)
pdf.plotOn(frame)
frame.Draw()

