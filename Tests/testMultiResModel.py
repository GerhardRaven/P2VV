import sys
import os
from ToyMCUtils import Toy

toy = Toy()
parser = toy.parser()
(options, args) = toy.configure()

from itertools import product
from RooFitWrappers import *
from P2VVLoad import P2VVLibrary
from ROOT import RooCBShape as CrystalBall
from P2VVParameterizations.GeneralUtils import valid_combinations

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
## from P2VVParameterizations.TimeResolution import Moriond2012_TimeResolution
## tres = Moriond2012_TimeResolution(time = t, timeResSFConstraint = True, sigmat = st,
##                                   timeResSF =  dict(Value = 1.46, MinMax = ( 0.5, 5. ),
##                                                     Constant = True))
## from P2VVParameterizations.TimeResolution import LP2011_TimeResolution
## tres = LP2011_TimeResolution(time = t)
from P2VVParameterizations.TimeResolution import Gaussian_TimeResolution as TimeResolution
tres = TimeResolution(time = t)

# Build the acceptance using the histogram as starting values
input_file = 'start_values.root'
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

from P2VVParameterizations.TimePDFs import Single_Exponent_Time
sig_t = Single_Exponent_Time(Name = 'sig_t_orig', time = t, resolutionModel = tres.model())

res_model = MultiHistEfficiencyModel(Name = "RMHE", Original = sig_t.pdf(), Observable = t,
                                     ConditionalCategories = True, UseSingleBinConstraint = False,
                                     Relative = rel_spec, Bins = bin_spec,
                                     ResolutionModel = tres.model())

# Signal time pdf
pdf = Single_Exponent_Time(Name = 'sig_t', time = t, resolutionModel = res_model)
pdf = pdf.pdf()

gen_observables = [t, hlt1_excl_biased, hlt2_unbiased, hlt2_biased]

data = pdf.generate(gen_observables, 20000)

states_signal = set([(state, label) for d in valid_definition for state, label in d])
def sort_combination(combination):
    valid_def = valid_definition[:]
    valid_def.reverse()
    level_left = 0
    n = 0
    c = set(combination)
    for level in valid_def:
        for j, state in enumerate(level):
            n |= int(state in c) << (level_left + j)
        level_left += len(level)
    return n - 1

def make_title(combination):
    title = []
    for level in valid_definition:
        l = level[0][0].GetName()[ : 4]
        level_states = set(level)
        s = [c for c in combination if c in level_states and c in states_signal]
        if len(s) == 1:
            title.append('%s_only_%s' % (l, s[0][1]))
        elif len(s) == 2:
            title.append('%s_both' % l)
    return '_X_'.join(title)

from ROOT import TCanvas, RooBinning
from ROOT import kDashed, kRed, kGreen, kBlue, kBlack

canv = TCanvas('canvas', 'canvas', 900, 700)
obs = [t]
for states, (p, o) in zip(sorted(rel_spec.keys(), key = sort_combination),
                          (i for i in product(canv.pads(3, 2), obs))):
    name = '__'.join(['%s_%s' % (state.GetName(), label) for state, label in states])
    title = make_title(states)
    cuts = ' && '.join(['{0} == {0}::{1}'.format(state.GetName(), label) for state, label in states])
    cat_data = data.reduce(cuts)
    project_set = RooArgSet(*project_vars)
    pdfOpts = dict(ProjWData = (project_set, cat_data, True))
    from P2VVGeneralUtils import plot
    binning = RooBinning(len(biased_bins) - 1, biased_bins)
    plot( p, o, cat_data, pdf, components = {'sig*' : dict(LineColor = kGreen, LineStyle = kDashed)}
          , plotResidHist = True
          , dataOpts = dict(MarkerSize = 0.8, MarkerColor = kBlack, Binning = binning)
          , frameOpts = {'Title' : title}
          , pdfOpts  = dict(LineWidth = 2, **pdfOpts)
          , logy = False
          , logx = True
          )
    
## run the toy
## toy.set_fit_opts(**dict(Verbose = False))
## toy.run(Observables = gen_observables, Pdf = pdf, GenPdf = pdf)

## toy.write_output()