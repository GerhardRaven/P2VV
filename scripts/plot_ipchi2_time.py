#!/usr/bin/env python
import gc
import optparse
import sys
import os
from math import sqrt

from P2VV.Utilities.Resolution import prefix
input_data = {'2011_Reco14' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_2011_Reco14_ntupleAB_20140430.root'),
              '2012' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_2012_ntupleAB_20140430.root')}

parser = optparse.OptionParser(usage = 'usage: %prog year')
(options, args) = parser.parse_args()

if len(args) != 1:
    print parser.print_usage()
    sys.exit(-2)

if args[0] not in input_data.keys():
    print parser.print_usage()
    print "Possible samples are: %s" % ' '.join(input_data.keys())
    sys.exit(-2)

from P2VV.Utilities.Resolution import prefix
from P2VV.Load import LHCbStyle
from P2VV.ROOTDecorators import *
from ROOT import TFile

input_file = TFile(input_data[args[0]])
tree = input_file.Get("DecayTree")

from ROOT import TH2D
from P2VV.Utilities.General import make_exp_binning
bins = make_exp_binning(10, 0.3, 14)
hist = TH2D("hist", "hist", len(bins) - 1, bins, 10, 0, 25)

from ROOT import TCanvas
canvas = TCanvas("canvas", "canvas", 600, 600)
canvas.SetLogx()
tree.Draw("B_s0_IPCHI2_OWNPV:time >> hist", "sweight * (sel == 1)", "box")
hists = [hist.ProjectionX("%s_%02d" % (hist.GetName(), i), i, i + 1, "oe") for i in range(1, len(bins))]
for h in hists:
    h.Scale(1 / h.Integral())
