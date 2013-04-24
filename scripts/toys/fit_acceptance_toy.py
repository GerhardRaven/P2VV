import sys
import os
from P2VV.ToyMCUtils import Toy

toy = Toy()
parser = toy.parser()
(options, args) = toy.configure()

from itertools import product
from P2VV.RooFitWrappers import *
from P2VV.Load import P2VVLibrary
from ROOT import RooCBShape as CrystalBall
from P2VV.Parameterizations.GeneralUtils import valid_combinations

obj = RooObject( workspace = 'w')
w = obj.ws()

from math import pi
t = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax=(0.3, 14))
m = RealVar('mass', Title = 'B mass', Unit = 'MeV', Observable = True, MinMax = (5250, 5550))
nPV = RealVar('nPV', Title = 'nPV', Observable = True, MinMax = (0, 15))
mpsi = RealVar('mdau1', Title = 'J/psi mass', Unit = 'MeV', Observable = True, MinMax = (3030, 3150))
st = RealVar('sigmat',Title = '#sigma(t)', Unit = 'ps', Observable = True, MinMax = (0.0001, 0.12))

# Categories
hlt1_biased = Category('hlt1_biased', States = {'biased' : 1, 'not_biased' : 0}, Observable = True)
hlt1_unbiased = Category('hlt1_unbiased', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
hlt1_excl_biased = Category('hlt1_excl_biased', States = {'excl_biased' : 1, 'unbiased' : 0}, Observable = True)
hlt2_biased = Category('hlt2_biased', States = {'biased' : 1, 'not_biased' : 0}, Observable = True)
hlt2_unbiased = Category('hlt2_unbiased', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)

categories = [hlt1_biased, hlt1_unbiased, hlt1_excl_biased, hlt2_biased, hlt2_unbiased]
categories = dict([(c.GetName(), c) for c in categories])

## project_vars = [hlt1_biased, hlt1_unbiased, hlt2_biased, hlt2_unbiased, st]
project_vars = [hlt1_excl_biased, hlt1_unbiased, hlt2_biased, hlt2_unbiased, st]

selected = Category('sel', States = {'Selected' : 1, 'NotSelected' : 0})

observables = [t, m, mpsi, st, hlt1_biased, hlt1_unbiased, hlt1_excl_biased,
               hlt2_biased, hlt2_unbiased, selected, nPV]

# now build the actual signal PDF...
from ROOT import RooGaussian as Gaussian
from ROOT import RooExponential as Exponential
from ROOT import RooDecay as Decay

signal_tau = RealVar('signal_tau', Title = 'mean lifetime', Unit = 'ps', Value =  1.5,
                     MinMax = (1., 2.5))

# Time resolution model
from P2VV.Parameterizations.TimeResolution import Moriond2012_TimeResolution
tres = Moriond2012_TimeResolution(time = t, timeResSFConstraint = True, sigmat = st,
                                  timeResSF =  dict(Value = 1.46, MinMax = ( 0.5, 5. )))
## from P2VV.Parameterizations.TimeResolution import LP2011_TimeResolution
## tres = LP2011_TimeResolution(time = t)
## from P2VV.Parameterizations.TimeResolution import Gaussian_TimeResolution as TimeResolution
## tres = TimeResolution(time = t)

# Signal time pdf
from P2VV.Parameterizations.TimePDFs import Single_Exponent_Time
sig_t = Single_Exponent_Time(Name = 'sig_t', time = t, resolutionModel = tres.model())

# Build the acceptance using the histogram as starting values
input_file = 'start_values_20bins.root'
hlt1_histogram = 'hlt1_shape'
hlt2_histogram = 'hlt2_shape'

from ROOT import TFile
acceptance_file = TFile.Open(input_file)
if not acceptance_file:
    raise ValueError, "Cannot open histogram file %s" % input_file
hlt1_histogram = acceptance_file.Get(hlt1_histogram)
if not hlt1_histogram:
    raise ValueError, 'Cannot get acceptance historgram %s from file' % hlt1_histogram
xaxis = hlt1_histogram.GetXaxis()
hlt2_histogram = acceptance_file.Get(hlt2_histogram)
if not hlt2_histogram:
    raise ValueError, 'Cannot get acceptance historgram %s from file' % hlt1_histogram

from array import array
biased_bins = array('d', (xaxis.GetBinLowEdge(i) for i in range(1, hlt1_histogram.GetNbinsX() + 2)))
unbiased_bins = array('d', [0.3, 14])

hlt1_biased_heights = [hlt1_histogram.GetBinContent(i) for i in range(1, hlt1_histogram.GetNbinsX() + 1)]
hlt1_unbiased_heights = [0.5]

hlt2_biased_heights = [hlt2_histogram.GetBinContent(i) for i in range(1, hlt2_histogram.GetNbinsX() + 1)]
hlt2_unbiased_heights = [0.5]

## valid_definition = [[(hlt1_biased, 'biased'), (hlt1_unbiased, 'unbiased')], [(hlt2_biased, 'biased'), (hlt2_unbiased, 'unbiased')]]
valid_definition = [[(hlt1_excl_biased, 'excl_biased'), (hlt1_excl_biased, 'unbiased')], [(hlt2_biased, 'biased'), (hlt2_unbiased, 'unbiased')]]
valid = valid_combinations(valid_definition)

bin_spec = {hlt1_excl_biased : {'excl_biased' : {'bins'    : biased_bins,
                                                 'heights' : hlt1_biased_heights,
                                                 'average' : (6.285e-01, 1.633e-02)},
                                'unbiased'    : {'bins'    : unbiased_bins,
                                                 'heights' : hlt1_unbiased_heights}
                                },
            hlt2_biased      : {'biased'      : {'bins'    : biased_bins,
                                                 'heights' : hlt2_biased_heights,
                                                 'average' : (6.329e-01, 1.3e-02)}
                                },
            hlt2_unbiased    : {'unbiased'    : {'bins'    : unbiased_bins,
                                                 'heights' : hlt2_unbiased_heights}
                                }
            }
        
## Fit options
rel_spec = {(('hlt1_excl_biased', 'excl_biased'), ('hlt2_biased', 'biased'), ('hlt2_unbiased', 'not_unbiased')) : 0.078,
            (('hlt1_excl_biased', 'unbiased'), ('hlt2_biased', 'not_biased'), ('hlt2_unbiased', 'unbiased')) : 0.027,
            (('hlt1_excl_biased', 'unbiased'), ('hlt2_biased', 'biased'), ('hlt2_unbiased', 'unbiased')) : 0.383,
            (('hlt1_excl_biased', 'excl_biased'), ('hlt2_biased', 'not_biased'), ('hlt2_unbiased', 'unbiased')) : 0.01,
            (('hlt1_excl_biased', 'unbiased'), ('hlt2_biased', 'biased'), ('hlt2_unbiased', 'not_unbiased')) : 0.433,
            (('hlt1_excl_biased', 'excl_biased'), ('hlt2_biased', 'biased'), ('hlt2_unbiased', 'unbiased')) : None}
rel_spec = dict([(tuple((categories[c], l) for c, l in k), {'Constant' : True, 'Value' : v} if v else None) for k, v in rel_spec.iteritems()])

res_model = MultiHistEfficiencyModel(Name = "RMHE", Original = sig_t.pdf(), Observable = t,
                                     ConditionalCategories = True, UseSingleBinConstraint = False,
                                     ResolutionModel = tres.model(), Bins = bin_spec,
                                     Relative = rel_spec)

pdf = Single_Exponent_Time(Name = 'pdf', time = t, resolutionModel = res_model)
pdf = pdf.pdf()

gen_observables = [t, hlt1_excl_biased, hlt2_unbiased, hlt2_biased]

## Get proto data
proto_data = None
if options.protodata:
    proto_file = TFile.Open(options.protodata)
    if not proto_file.IsOpen():
        raise OSError
    for key in proto_file.GetListOfKeys():
        if key.GetClassName() == "RooDataSet":
            proto_data = key.ReadObj()
            break
    else:
        raise RunTimeError
    if not proto_data:
        raise RunTimeError

## run the toy
toy.set_fit_opts(**dict(Verbose = False))
toy.run(Observables = gen_observables, Pdf = pdf, GenPdf = pdf, ProtoData = proto_data)

toy.write_output()
