#!/usr/bin/env python
from array import array
import optparse
import sys
import os
from math import sqrt

parser = optparse.OptionParser(usage = '%prog data_type')
(options, args) = parser.parse_args()

prefix = '/stuff/PhD' if os.path.exists('/stuff') else '/project/bfys/raaij'
input_data = {'2011_Reco14' : {'data' : os.path.join(prefix, 'p2vv/data/P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_HLT2B_20140309.root'),
                               'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_2011_Reco14_Prescaled_cache.root'),
                               'dataset' : '1bin_15500.00fs_simple/2027465761870101697/sdata/sig_sdata'},
              '2012' : {'data' : os.path.join(prefix, 'p2vv/data/P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_HLT2B_20140309.root'),
                               'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_2012_Prescaled_cache.root'),
                               'dataset' : '1bin_15500.00fs_simple/2027465761870101697/sdata/sig_sdata'},
              'MC2011_Sim08a' : {'sig_cache'      : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_MC2011_Sim08a_Prescaled_cache.root'),
                                 'sig_dataset'    : '1bin_19000.00fs_simple/m4074732224475561764/sdata/sig_sdata',
                                 'prompt_cache'   : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_MC2011_Sim08a_incl_Jpsi_Prescaled_cache.root'),
                                 'prompt_dataset' : '1bin_15500.00fs_simple/2027465761870101697/sdata/sig_sdata'},
              'MC2012'        : {'sig_cache'      : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_MC2012_Prescaled_cache.root'),
                                 'sig_dataset'    : '1bin_19000.00fs_simple/m4074732224475561764/sdata/sig_sdata',
                                 'prompt_cache'   : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_MC2012_incl_Jpsi_Prescaled_cache.root'),
                                 'prompt_dataset' : '1bin_15500.00fs_simple/2027465761870101697/sdata/sig_sdata'}
              }

if len(args) != 1:
    print parser.print_usage()
    sys.exit(-2)
elif args[0] not in input_data:
    print parser.print_usage()
    sys.exit(-2)

from P2VV.Load import LHCbStyle
from P2VV.RooFitWrappers import *

from ROOT import RooRealVar
from ROOT import RooDataSet
from ROOT import RooArgSet
from ROOT import RooCategory
input_data = input_data[args[0]]

weight = RooRealVar('sWeights_ipatia', 'sWeights_ipatia', -1e3, 1e3)
momentum = RooRealVar('B_P', 'B_P', 0, 1e6, 'MeV')
runPeriod = RooCategory('runPeriod', 'runPeriod')
runPeriod.defineType('p2011', 2011)
runPeriod.defineType('p2012', 2012)

from ROOT import TFile
if args[0].startswith('MC'):
    from ROOT import TFile
    sig_file = TFile(input_data['sig_cache'])
    sig_data = sig_file.Get(input_data['sig_dataset'])
    prompt_file = TFile(input_data['prompt_cache'])
    prompt_data = prompt_file.Get(input_data['prompt_dataset'])
else:
    sig_file = TFile("/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/fitNTuple_peakBkg_2011_2012_Reco14_TOS_HLT2B_20140415.root")
    sig_tree = sig_file.Get("DecayTree")
    period = 'p' + args[0][:4]
    sig_data = RooDataSet('sig_data', 'sig_data', RooArgSet(momentum, weight, runPeriod),
                          RooFit.Import(sig_tree), RooFit.WeightVar(weight),
                          RooFit.Cut('runPeriod == runPeriod::%s' % period))

    from ROOT import TFile
    prompt_file = TFile(input_data['cache'])
    prompt_data = prompt_file.Get(input_data['dataset'])

    
from ROOT import TCanvas
canvas = TCanvas('momentum_canvas', 'momentum_canvas', 600, 400)
canvas.SetLogy()
frame = momentum.frame(Range = (0, 500000))
from ROOT import kBlue, kGreen
sig_data.plotOn(frame, Rescale = (1 / sig_data.sumEntries()), MarkerColor = kBlue)
prompt_data.plotOn(frame, Rescale = (1 / prompt_data.sumEntries()), MarkerColor = kGreen)
frame.GetXaxis().SetTitle("P_{B} [MeV]")
frame.GetYaxis().SetTitle("scaled candidates")
frame.GetYaxis().SetRangeUser(0.00005, 0.1)
frame.Draw()
