import os
import sys
if sys.argv[1] not in ['MC', 'real_data', 'generate']:
    print 'usage: fit_acceptance_multi.py [real_data|MC|generate]'
    sys.exit(-1)

real_data = sys.argv[1] == 'real_data'
MC = sys.argv[1] == 'MC'

from itertools import product
from RooFitWrappers import *
from P2VVLoad import P2VVLibrary
from ROOT import RooCBShape as CrystalBall
from P2VVParameterizations.GeneralUtils import valid_combinations

## from ROOT import RooMsgService
## RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.Eval))
## RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.Integration))

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
## project_vars = [hlt1_biased, hlt1_unbiased, hlt2_biased, hlt2_unbiased, st]
categories = [hlt1_biased, hlt1_unbiased, hlt1_excl_biased, hlt2_biased, hlt2_unbiased]
categories = dict([(c.GetName(), c) for c in categories])

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
from P2VVParameterizations.TimeResolution import Moriond2012_TimeResolution
tres = Moriond2012_TimeResolution(time = t, timeResSFConstraint = True, sigmat = st,
                                  timeResSF =  dict(Value = 1.46, MinMax = ( 0.5, 5. ),
                                                    Constant = True))
## from P2VVParameterizations.TimeResolution import LP2011_TimeResolution
## tres = LP2011_TimeResolution(time = t)
## from P2VVParameterizations.TimeResolution import Gaussian_TimeResolution as TimeResolution
## tres = TimeResolution(time = t)

# Signal time pdf
from P2VVParameterizations.TimePDFs import Single_Exponent_Time
sig_t = Single_Exponent_Time(Name = 'sig_t', time = t, resolutionModel = tres.model())

# B mass pdf
from P2VVParameterizations.MassPDFs import LP2011_Signal_Mass as Signal_BMass, LP2011_Background_Mass as Background_BMass
sig_m = Signal_BMass(     Name = 'sig_m', mass = m, m_sig_mean = dict( Value = 5365, MinMax = (5363,5372) ) )

# J/psi mass pdf
mpsi_mean  = RealVar('mpsi_mean',   Unit = 'MeV', Value = 3097, MinMax = (3070, 3110))
mpsi_sigma = RealVar('mpsi_sigma',  Unit = 'MeV', Value = 10, MinMax = (5, 20))
mpsi_alpha = RealVar('mpsi_alpha',  Unit = '', Value = 1.8, MinMax = (0.5, 3), Constant = True)
mpsi_n = RealVar('mpsi_n',  Unit = '', Value = 2, MinMax = (0.1, 4), Constant = True)
psi_m  = Pdf(Name = 'psi_m', Type = CrystalBall, Parameters = [mpsi, mpsi_mean, mpsi_sigma, mpsi_alpha, mpsi_n])

# J/psi background
psi_c = RealVar( 'psi_c',  Unit = '1/MeV', Value = -0.0004, MinMax = (-0.1, -0.0000001))
bkg_mpsi = Pdf(Name = 'bkg_mpsi',  Type = Exponential, Parameters = [mpsi, psi_c])

# Create psi background component
from P2VVParameterizations.TimePDFs import LP2011_Background_Time as Background_Time
psi_t = Background_Time( Name = 'psi_t', time = t, resolutionModel = tres.model()
                         , psi_t_fll = dict( Name = 'psi_t_fll',    Value = 0.2 )
                         , psi_t_ll_tau = dict( Name = 'psi_t_ll_tau', Value = 1.25, MinMax = (0.5,2.5) )
                         , psi_t_ml_tau = dict( Name = 'psi_t_ml_tau', Value = 0.16, MinMax = (0.01,0.5) )
                         , ExternalConstraints = tres.model().ExternalConstraints())

# Create combinatorical background component
bkg_m = Background_BMass( Name = 'bkg_m', mass = m, m_bkg_exp  = dict( Name = 'm_bkg_exp' ) )

bkg_tau = RealVar('bkg_tau', Title = 'comb background lifetime', Unit = 'ps', Value = 1, MinMax = (0.0001, 5))

bkg_t = Background_Time( Name = 'bkg_t', time = t, resolutionModel = tres.model()
                       , t_bkg_fll    = dict( Name = 't_bkg_fll',    Value = 0.3 )
                       , t_bkg_ll_tau = dict( Name = 't_bkg_ll_tau', Value = 1.92, MinMax = (0.5,2.5) )
                       , t_bkg_ml_tau = dict( Name = 't_bkg_ml_tau', Value = 0.21, MinMax = (0.01,0.5) ) )

# Create components
signal = Component('signal', (sig_m.pdf(), psi_m, sig_t.pdf()), Yield = (30000,100,100000))
psi_background = Component('psi_background', (bkg_m.pdf(), psi_m, psi_t.pdf()), Yield= (100000,500,200000) )
background = Component('background', (bkg_m.pdf(), bkg_mpsi, bkg_t.pdf()), Yield = (100000,100,300000) )

## Build mass PDF
mass_pdf = buildPdf(Components = (signal, background), Observables = (m, ), Name='mass_pdf')
mass_pdf.Print("t")

## base_location = '/home/raaij'
base_location = '/stuff/PhD/p2vv'

# Build the acceptance using the histogram as starting values
input_file = os.path.join(base_location, 'data/start_values.root')
hlt1_histogram = 'hlt1_shape'
hlt2_histogram = 'hlt2_shape'

## if MC:
    ## hlt1_histogram += '_MC'
    ## hlt2_histogram += '_MC'

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

spec = {'Bins' : {hlt1_excl_biased : {'excl_biased' : {'bins'    : biased_bins,
                                                       'heights' : hlt1_biased_heights,
                                                       'average' : (7.285e-01, 1.633e-02)},
                                      'unbiased'    : {'bins'    : unbiased_bins,
                                                       'heights' : hlt1_unbiased_heights}
                                      },
                  hlt2_biased      : {'biased'      : {'bins'    : biased_bins,
                                                       'heights' : hlt2_biased_heights,
                                                       'average' : (7.330e-01, 1.402e-02)}
                                      },
                  hlt2_unbiased    : {'unbiased'    : {'bins'    : unbiased_bins,
                                                       'heights' : hlt2_unbiased_heights}
                                      }
                  }
        }
        
# Read input data
from P2VVGeneralUtils import readData
tree_name = 'DecayTree'
## input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhiPrescaled_ntupleB_for_fitting_20120110.root'

## Fit options
fitOpts = dict(NumCPU = 4, Timer = 1, Save = True, Verbose = True, Optimize = 1, Minimizer = 'Minuit2')

data = None
if real_data:
    input_file = os.path.join(base_location, 'data/Bs2JpsiPhi_2011_biased_unbiased.root')
    data = readData(input_file, tree_name, cuts = 'sel == 1 && (hlt1_biased == 1 || hlt1_unbiased == 1) && (hlt2_biased == 1 || hlt2_unbiased == 1)',
                    NTuple = True, observables = observables)
    total = data.sumEntries()

    rel_spec = {}
    for comb in valid:
        cuts = ' && '.join(['{0} == {0}::{1}'.format(state.GetName(), label) for state, label in comb])
        rel_spec[comb] = {'Value' : data.sumEntries(cuts) / total, "Constant" : True}

    spec['Relative'] = rel_spec
    pdf = MultiHistEfficiency(Name = "RMHE", Original = sig_t.pdf(), Observable = t,
                              ConditionalCategories = True, UseSingleBinConstraint = False,
                              **spec)
    pdf.Print('v')

    mass_pdf.fitTo(data, **fitOpts)
    # Plot mass pdf
    from ROOT import kDashed, kRed, kGreen, kBlue, kBlack
    from ROOT import TCanvas
    canvas = TCanvas('mass_canvas', 'mass_canvas', 1000, 500)
    obs = [m, mpsi]
    for (p,o) in zip(canvas.pads(len(obs)), obs):
        from P2VVGeneralUtils import plot
        pdfOpts  = dict()
        plot(p, o, pdf = mass_pdf, data = data
             , dataOpts = dict(MarkerSize = 0.8, MarkerColor = kBlack)
             , pdfOpts  = dict(LineWidth = 2, **pdfOpts)
             , plotResidHist = True
             , components = { 'bkg_*'     : dict( LineColor = kRed,   LineStyle = kDashed ),
                              'psi_*'  : dict( LineColor = kGreen, LineStyle = kDashed ),
                              'sig_*'     : dict( LineColor = kBlue,  LineStyle = kDashed )
                              }
             )
    # Do the sWeights
    # make sweighted dataset. TODO: use mumu mass as well...
    from P2VVGeneralUtils import SData, splot

    for p in mass_pdf.Parameters() : p.setConstant( not p.getAttribute('Yield') )
    splot = SData(Pdf = mass_pdf, Data = data, Name = 'MassSplot')
    data = splot.data('signal')
    ## psi_sdata = splot.data('psi_background')
    bkg_sdata = splot.data('background')
elif MC:
    input_file = os.path.join(base_location, 'data/Bs2JpsiPhi_MC11a_biased_unbiased.root')
    data = readData(input_file, tree_name, cuts = 'sel == 1 && (hlt1_biased == 1 || hlt1_unbiased == 1) && hlt2_unbiased == 1',
                    NTuple = True, observables = observables)

    import random
    # generate some efficiency as a function of t
    # make a NEW dataset with hit-miss on the efficienty, add 1/eff as weight 
    new_data = RooDataSet("new_data", "new_data", data.get())
    for i, obs in enumerate(data):
        b2 = obs.find('hlt2_biased')
        ub2 = obs.find('hlt2_unbiased')
        if b2.getIndex() == 0:
            pass
        elif random.random() < 0.5:
            ub2.setIndex(0)
        new_data.add(obs)
        if i >= 50000:
            break

    new_data.table(RooArgSet(hlt1_excl_biased, hlt2_biased, hlt2_unbiased)).Print('v')
    del data
    data = new_data
    total = data.sumEntries()

    rel_spec = {}
    for comb in valid:
        cuts = ' && '.join(['{0} == {0}::{1}'.format(state.GetName(), label) for state, label in comb])
        ## rel_spec[comb] = {'Value' : data.sumEntries(cuts) / total, "Constant" : True}
        rel_spec[comb] = {'Value' : 0.5, "Constant" : True}

    spec['Relative'] = rel_spec
    pdf = MultiHistEfficiency(Name = "RMHE", Original = sig_t.pdf(), Observable = t,
                              ConditionalCategories = True, UseSingleBinConstraint = False,
                              **spec)
    pdf.Print('v')
else:
    rel_spec = {(('hlt1_excl_biased', 'excl_biased'), ('hlt2_biased', 'biased'), ('hlt2_unbiased', 'not_unbiased')) : 0.078,
                (('hlt1_excl_biased', 'unbiased'), ('hlt2_biased', 'not_biased'), ('hlt2_unbiased', 'unbiased')) : 0.027,
                (('hlt1_excl_biased', 'unbiased'), ('hlt2_biased', 'biased'), ('hlt2_unbiased', 'unbiased')) : 0.383,
                (('hlt1_excl_biased', 'excl_biased'), ('hlt2_biased', 'not_biased'), ('hlt2_unbiased', 'unbiased')) : 0.01,
                (('hlt1_excl_biased', 'unbiased'), ('hlt2_biased', 'biased'), ('hlt2_unbiased', 'not_unbiased')) : 0.433,
                (('hlt1_excl_biased', 'excl_biased'), ('hlt2_biased', 'biased'), ('hlt2_unbiased', 'unbiased')) : None}
    spec['Relative'] = dict([(tuple((categories[c], l) for c, l in k), {'Constant' : True, 'Value' : v} if v else None) for k, v in rel_spec.iteritems()])
    for comb in valid:
        cuts = ' && '.join(['{0} == {0}::{1}'.format(state.GetName(), label) for state, label in comb])
        rel_spec[comb] = {'Value' : 1. / len(valid), "Constant" : True},
    pdf = MultiHistEfficiency(Name = "RMHE", Original = sig_t.pdf(), Observable = t,
                              ConditionalCategories = True, UseSingleBinConstraint = False,
                              **spec)
    pdf.Print('v')        
    data = pdf.generate([t, hlt1_excl_biased, hlt2_unbiased, hlt2_biased], 30000)

## Fit
print 'fitting data'
## from profiler import profiler_start, profiler_stop
## profiler_start("acceptance.log")
result = pdf.fitTo(data, **fitOpts)
## profiler_stop()

from ROOT import kDashed, kRed, kGreen, kBlue, kBlack
from ROOT import TCanvas, RooBinning
canvas = {}
print 'plotting'

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
    
# Plot the lifetime shapes
canv = TCanvas('canvas', 'canvas', 900, 700)
obs = [t]
for states, (p, o) in zip(sorted(spec['Relative'].keys(), key = sort_combination),
                          (i for i in product(canv.pads(3, 2), obs))):
    name = '__'.join(['%s_%s' % (state.GetName(), label) for state, label in states])
    title = make_title(states)
    cuts = ' && '.join(['{0} == {0}::{1}'.format(state.GetName(), label) for state, label in states])
    cat_data = data.reduce(cuts)
    project_set = RooArgSet(*project_vars)
    pdfOpts = dict(ProjWData = (project_set, cat_data))
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
    p.SetLogx(1)
    p.Update()
    
# plot the efficiency shapes
def plot_shape(p, o, shape, errorOpts = {}, pdfOpts = {}):
    from operator import itemgetter
    i = shape.createIntegral(RooArgSet(o))
    n = i.getVal()
    p.cd()
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
    frame.Draw()

shapes = pdf.shapes()
eff_canvas = TCanvas('eff_canvas', 'eff_canvas', 1000, 500)
from ROOT import kYellow, kOrange
for p, shape in zip(eff_canvas.pads(len(shapes), 1), shapes):
    plot_shape(p, t, shape, errorOpts = {'result' : result, 3 : kYellow, 1 : kOrange})

output = {'hlt1_shape' : 'hlt1_excl_biased_excl_biased_bin_%03d',
          'hlt2_shape' : 'hlt2_biased_biased_bin_%03d'}
output_file = TFile.Open('efficiencies.root', 'recreate')
from ROOT import TH1D
for name, pat in output.iteritems():
    n = len(biased_bins)
    heights = [obj.ws().var(pat % i) for i in range(1, n)]
    v = [(h.getVal(), h.getError()) for h in heights]
    hist = TH1D(name, name, n - 1, biased_bins)
    for i in range(1, n):
        hist.SetBinContent(i, v[i - 1][0])
        hist.SetBinError(i, v[i - 1][1])
    output_file.WriteTObject(hist)
output_file.Close()
