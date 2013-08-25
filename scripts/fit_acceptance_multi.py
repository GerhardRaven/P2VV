import os
import sys
if sys.argv[1] not in ['MC', 'data', 'generate']:
    print 'usage: fit_acceptance_multi.py [data|MC|generate]'
    sys.exit(-1)

real_data = sys.argv[1] == 'data'
MC = sys.argv[1] == 'MC'

from itertools import product
from P2VV.RooFitWrappers import *
from ROOT import RooCBShape as CrystalBall
from P2VV.Parameterizations.GeneralUtils import valid_combinations

from ROOT import RooMsgService
## RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.ObjectHandling))
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
hlt1_unbiased = Category('hlt1_unbiased_dec', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
hlt1_excl_biased_dec = Category('hlt1_excl_biased_dec', States = {'excl_biased' : 1, 'unbiased' : 0}, Observable = True)
hlt2_biased = Category('hlt2_biased', States = {'biased' : 1, 'not_biased' : 0}, Observable = True)
hlt2_unbiased = Category('hlt2_unbiased', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
hlt2_excl_biased = Category('hlt2_excl_biased', States = {'excl_biased' : 1, 'unbiased' : 0}, Observable = True)

## project_vars = [hlt1_biased, hlt1_unbiased, hlt2_biased, hlt2_unbiased, st]
categories = [hlt1_biased, hlt1_unbiased, hlt1_excl_biased_dec,
              hlt2_biased, hlt2_unbiased, hlt2_excl_biased]
categories = dict([(c.GetName(), c) for c in categories])

project_vars = [hlt1_excl_biased_dec, hlt1_unbiased, hlt2_biased, hlt2_unbiased, st]

selected = Category('sel', States = {'Selected' : 1, 'NotSelected' : 0})

observables = [t, m, mpsi, st, hlt1_biased, hlt1_unbiased, hlt1_excl_biased_dec,
               hlt2_biased, hlt2_unbiased, hlt2_excl_biased, selected, nPV]

# now build the actual signal PDF...
from ROOT import RooGaussian as Gaussian
from ROOT import RooExponential as Exponential
from ROOT import RooDecay as Decay

signal_tau = RealVar('signal_tau', Title = 'mean lifetime', Unit = 'ps', Value =  1.5,
                     MinMax = (1., 2.5))

# Time resolution model
from P2VV.Parameterizations.TimeResolution import Paper2012_TimeResolution as TimeResolution
tres = TimeResolution(time = t, timeResSigma = st, Cache = False,
                      timeResComb = dict(Value = 1.4918, Error = 4.08e-03, MinMax = ( 0.1, 5. ), Constant = True),
                      timeResSigmaSF2 = dict(Value = 6.0074, Error = 1.89e-01, MinMax = (1, 10), Constant = True),
                      timeResSigmaFrac2 = dict(Value = 1.5818e-02, Error = 1.07e-03, MinMax = (0.001, 0.999), Constant = True),
                      Covariance = {('timeResComb', 'timeResComb'): 1.663e-05,
                                    ('timeResComb', 'timeResSigmaFrac2'): 1.322e-06,
                                    ('timeResComb', 'timeResSigmaSF2'): 0.0001297,
                                    ('timeResSigmaFrac2', 'timeResSigmaFrac2'): 1.146e-06,
                                    ('timeResSigmaFrac2', 'timeResSigmaSF2'): -0.0001486,
                                    ('timeResSigmaSF2', 'timeResSigmaSF2'): 0.03556},
                      nGauss = 2, timeResMean = dict(Value = -4.0735e-03, Error = 1.33e-04),
                      timeResMeanConstraint = 'fixed')

## from P2VV.Parameterizations.TimeResolution import LP2011_TimeResolution
## tres = LP2011_TimeResolution(time = t)
## from P2VV.Parameterizations.TimeResolution import Gaussian_TimeResolution as TimeResolution
## tres = TimeResolution(time = t)

# Signal time pdf
from P2VV.Parameterizations.TimePDFs import Single_Exponent_Time
sig_t = Single_Exponent_Time(Name = 'sig_t', time = t, resolutionModel = tres.model())

# B mass pdf
from P2VV.Parameterizations.MassPDFs import LP2011_Signal_Mass as Signal_BMass, LP2011_Background_Mass as Background_BMass
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

# Create combinatorical background component
bkg_m = Background_BMass( Name = 'bkg_m', mass = m, m_bkg_exp  = dict( Name = 'm_bkg_exp' ) )

# Create components
signal_mass = Component('signal', (sig_m.pdf(), psi_m), Yield = (30000,100,100000))
psi_background_mass = Component('psi_background', (bkg_m.pdf(), psi_m), Yield= (100000,500,200000) )
background_mass = Component('background', (bkg_m.pdf(), bkg_mpsi), Yield = (100000,100,300000) )

## Build mass PDF
mass_pdf = buildPdf(Components = (signal_mass, background_mass), Observables = (m, ), Name='mass_pdf')
mass_pdf.Print("t")

## base_location = '/home/raaij'
base_location = '/stuff/PhD/p2vv'

# Build the acceptance using the histogram as starting values
input_file = os.path.join(base_location, 'data/start_values.root')
hlt1_histogram = 'hlt1_shape'
hlt2_histogram = 'hlt2_shape'

## input_file = os.path.join(base_location, 'data/Bs_HltPropertimeAcceptance_Data-20120816.root')
## hlt1_histogram = 'hlt1_shape'
## hlt2_histogram = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_Data_20bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'

## if MC:
    ## hlt1_histogram += '_MC'
    ## hlt2_histogram += '_MC'

from ROOT import TFile
acceptance_file = TFile.Open(input_file)
if not acceptance_file:
    raise ValueError, "Cannot open histogram file %s" % input_file
hlt2_histogram = acceptance_file.Get(hlt2_histogram)
if not hlt2_histogram:
    raise ValueError, 'Cannot get acceptance historgram %s from file' % hlt2_histogram
xaxis = hlt2_histogram.GetXaxis()

hlt1_histogram = acceptance_file.Get(hlt1_histogram)
if not hlt1_histogram:
    raise ValueError, 'Cannot get acceptance historgram %s from file' % hlt1_histogram


from array import array
biased_bins = array('d', (xaxis.GetBinLowEdge(i) for i in range(1, hlt2_histogram.GetNbinsX() + 2)))
unbiased_bins = array('d', [0.3, 14])

hlt1_biased_heights = [hlt1_histogram.GetBinContent(i) for i in range(1, hlt1_histogram.GetNbinsX() + 1)]
hlt1_biased_average = (6.285e-01, 1.633e-02)
## hlt1_biased_heights = [hlt1_biased_average[0] for i in range(1, hlt2_histogram.GetNbinsX() + 1)]

hlt1_unbiased_heights = [0.5]

hlt2_biased_average = (6.33e-01, 1.65e-02)
hlt2_biased_heights = [hlt2_histogram.GetBinContent(i) for i in range(1, hlt2_histogram.GetNbinsX() + 1)]
av = sum(hlt2_biased_heights) / float(len(hlt2_biased_heights))
scale = av / hlt2_biased_average[0]
hlt2_biased_heights = [h / scale for h in hlt2_biased_heights]

hlt2_unbiased_heights = [0.5]

## valid_definition = [[(hlt1_biased, 'biased'), (hlt1_unbiased, 'unbiased')], [(hlt2_biased, 'biased'), (hlt2_unbiased, 'unbiased')]]
## valid_definition = [[(hlt2_biased, 'biased'), (hlt2_unbiased, 'unbiased')]]
valid_definition = [[(hlt1_excl_biased_dec, 'excl_biased'), (hlt1_excl_biased_dec, 'unbiased')], [(hlt2_biased, 'biased'), (hlt2_unbiased, 'unbiased')]]
## valid_definition = [[(hlt1_excl_biased_dec, 'excl_biased'), (hlt1_excl_biased_dec, 'unbiased')], [(hlt2_excl_biased, 'excl_biased'), (hlt2_excl_biased, 'unbiased')]]
valid = valid_combinations(valid_definition)

spec = {'Bins' : {hlt1_excl_biased_dec : {'excl_biased' : {'bins'    : biased_bins,
                                                       'heights' : hlt1_biased_heights,
                                                       'average' : hlt1_biased_average},
                                      'unbiased'    : {'bins'    : unbiased_bins,
                                                       'heights' : hlt1_unbiased_heights}
                                      },
                  hlt2_biased      : {'biased'      : {'bins'    : biased_bins,
                                                       'heights' : hlt2_biased_heights,
                                                       'average' : hlt2_biased_average}
                                      },
                  hlt2_unbiased    : {'unbiased'    : {'bins'    : unbiased_bins,
                                                       'heights' : hlt2_unbiased_heights}
                                      }
                  }
        }

## spec = {'Bins' : {hlt2_biased      : {'biased'      : {'bins'    : biased_bins,
##                                                        'heights' : hlt2_biased_heights,
##                                                        'average' : (6.330e-01, 1.65e-02)}
##                                       },
##                   hlt2_unbiased    : {'unbiased'    : {'bins'    : unbiased_bins,
##                                                        'heights' : hlt2_unbiased_heights}
##                                       }
##                   }
##         }

## spec = {'Bins' : {hlt1_excl_biased_dec : {'excl_biased' : {'bins'    : biased_bins,
##                                                        'heights' : hlt1_biased_heights,
##                                                        'average' : (6.285e-01, 1.633e-02)},
##                                       'unbiased'    : {'bins'    : unbiased_bins,
##                                                        'heights' : hlt1_unbiased_heights}
##                                       },
##                   hlt2_excl_biased : {'excl_biased' : {'bins'    : biased_bins,
##                                                        'heights' : hlt2_biased_heights,
##                                                        'average' : (6.330e-01, 1.65e-02)},
##                                       'unbiased'    : {'bins'    : unbiased_bins,
##                                                        'heights' : hlt2_unbiased_heights}
##                                       }
##                   }
##         }

# Read input data
from P2VV.Utilities.DataHandling import readData
tree_name = 'DecayTree'
## input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhiPrescaled_ntupleB_for_fitting_20120110.root'

## Fit options
fitOpts = dict(NumCPU = 4, Timer = 1, Save = True, Verbose = True, Optimize = 1,
               Strategy = 2, Minimizer = 'Minuit2')

data = None
if real_data:
    ## input_file = os.path.join(base_location, 'data/Bs2JpsiPhi_2011_biased_unbiased.root')
    input_file = os.path.join(base_location, 'data/Bs2JpsiPhi_ntupleB_for_fitting_20121012_MagDownMagUp.root')
    data = readData(input_file, tree_name, cuts = 'sel == 1 && (hlt1_biased == 1 || hlt1_unbiased_dec == 1) && (hlt2_biased == 1 || hlt2_unbiased == 1)',
                    NTuple = True, observables = observables)
    ## data = readData(input_file, tree_name, cuts = 'sel == 1 && hlt1_unbiased == 1 && (hlt2_biased == 1 || hlt2_unbiased == 1)',
    ##                 NTuple = True, observables = observables)
    total = data.sumEntries()

    rel_spec = {}
    for comb in valid:
        cuts = ' && '.join(['{0} == {0}::{1}'.format(state.GetName(), label) for state, label in comb])
        rel_spec[comb] = {'Value' : data.sumEntries(cuts) / total, "Constant" : True}

    spec['Relative'] = rel_spec
    res_model = MultiHistEfficiencyModel(Name = "RMHE", Original = sig_t.pdf(), Observable = t,
                                         ConditionalCategories = True, UseSingleBinConstraint = False,
                                         ResolutionModel = tres.model(), Spline = True,
                                         SmoothSpline = 1, **spec)
    pdf = Single_Exponent_Time(Name = 'pdf', time = t, resolutionModel = res_model)
    pdf = pdf.pdf()
    pdf.Print('v')

    mass_pdf.fitTo(data, **fitOpts)
    # Plot mass pdf
    from ROOT import kDashed, kRed, kGreen, kBlue, kBlack
    from ROOT import TCanvas
    canvas = TCanvas('mass_canvas', 'mass_canvas', 500, 500)
    obs = [m]
    for (p,o) in zip(canvas.pads(len(obs)), obs):
        from P2VV.Utilities.Plotting import plot
        pdfOpts  = dict()
        plot(p, o, pdf = mass_pdf, data = data
             , dataOpts = dict(MarkerSize = 0.8, MarkerColor = kBlack)
             , pdfOpts  = dict(LineWidth = 2, **pdfOpts)
             , plotResidHist = True
             , components = { 'bkg_*'     : dict( LineColor = kRed,   LineStyle = kDashed ),
                              ## 'psi_*'  : dict( LineColor = kGreen, LineStyle = kDashed ),
                              'sig_*'     : dict( LineColor = kBlue,  LineStyle = kDashed )
                              }
             )
    # Do the sWeights
    # make sweighted dataset. TODO: use mumu mass as well...
    from P2VV.Utilities.SWeights import SData, splot

    for p in mass_pdf.Parameters() : p.setConstant( not p.getAttribute('Yield') )
    splot = SData(Pdf = mass_pdf, Data = data, Name = 'MassSplot')
    data = splot.data('signal')
    ## psi_sdata = splot.data('psi_background')
    bkg_sdata = splot.data('background')
elif MC:
    input_file = os.path.join(base_location, 'data/Bs2JpsiPhi_MC11a_biased_unbiased.root')
    data = readData(input_file, tree_name, cuts = 'sel == 1 && (hlt1_biased == 1 || hlt1_unbiased == 1) && (hlt2_biased == 1 || hlt2_unbiased == 1)',
                    NTuple = True, observables = observables)
    import random
    # generate some efficiency as a function of t
    # make a NEW dataset with hit-miss on the efficienty, add 1/eff as weight 
    new_data = RooDataSet("new_data", "new_data", data.get())
    for i, obs in enumerate(data):
        b2 = obs.find('hlt2_biased')
        ub2 = obs.find('hlt2_unbiased')
        eb2 = obs.find('hlt2_excl_biased')
        if b2.getIndex() == 0:
            pass
        elif random.random() < 0.5:
            ub2.setIndex(0)
            eb2.setIndex(1)
        new_data.add(obs)
        if i >= 50000:
            break

    new_data.table(RooArgSet(hlt1_excl_biased_dec, hlt2_unbiased, hlt2_biased)).Print('v')
    del data
    data = new_data
    total = data.sumEntries()

    rel_spec = {}
    for comb in valid:
        cuts = ' && '.join(['{0} == {0}::{1}'.format(state.GetName(), label) for state, label in comb])
        rel_spec[comb] = {'Value' : data.sumEntries(cuts) / total, "Constant" : True}
        ## rel_spec[comb] = {'Value' : 0.5, "Constant" : True}

    spec['Relative'] = rel_spec
    res_model = MultiHistEfficiencyModel(Name = "RMHE", Original = sig_t.pdf(), Observable = t,
                                         ConditionalCategories = True, UseSingleBinConstraint = False,
                                         ResolutionModel = tres.model(), **spec)
    pdf = Single_Exponent_Time(Name = 'pdf', time = t, resolutionModel = res_model)
    pdf = pdf.pdf()
    
    pdf.Print('v')
else:
    rel_spec = {(('hlt1_excl_biased_dec', 'excl_biased'), ('hlt2_biased', 'biased'), ('hlt2_unbiased', 'not_unbiased')) : 0.078,
                (('hlt1_excl_biased_dec', 'unbiased'), ('hlt2_biased', 'not_biased'), ('hlt2_unbiased', 'unbiased')) : 0.027,
                (('hlt1_excl_biased_dec', 'unbiased'), ('hlt2_biased', 'biased'), ('hlt2_unbiased', 'unbiased')) : 0.383,
                (('hlt1_excl_biased_dec', 'excl_biased'), ('hlt2_biased', 'not_biased'), ('hlt2_unbiased', 'unbiased')) : 0.01,
                (('hlt1_excl_biased_dec', 'unbiased'), ('hlt2_biased', 'biased'), ('hlt2_unbiased', 'not_unbiased')) : 0.433,
                (('hlt1_excl_biased_dec', 'excl_biased'), ('hlt2_biased', 'biased'), ('hlt2_unbiased', 'unbiased')) : None}
    spec['Relative'] = dict([(tuple((categories[c], l) for c, l in k), {'Constant' : True, 'Value' : v} if v else None) for k, v in rel_spec.iteritems()])
    for comb in valid:
        cuts = ' && '.join(['{0} == {0}::{1}'.format(state.GetName(), label) for state, label in comb])
        rel_spec[comb] = {'Value' : 1. / len(valid), "Constant" : True},

    from P2VV.Parameterizations.TimeResolution import Gaussian_TimeResolution as TimeResolution
    tres = TimeResolution(time = t, timeResSigma = dict(Value = 0.5, Constant = True),
                          timeResMu = dict(Value = 0, Constant = True),
                          BiasScaleFactor = False)

    res_model = MultiHistEfficiencyModel(Name = "RMHE", Original = sig_t.pdf(), Observable = t,
                                         ConditionalCategories = True, UseSingleBinConstraint = False,
                                         ResolutionModel = tres.model(), **spec)
    pdf = Single_Exponent_Time(Name = 'pdf', time = t, resolutionModel = res_model)
    pdf = pdf.pdf()

    pdf.Print('v')        
    data = pdf.generate([t, hlt1_excl_biased_dec, hlt2_unbiased, hlt2_biased], 30000)

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
    pdfOpts = dict(ProjWData = (project_set, cat_data, True))
    from P2VV.Utilities.Plotting import plot
    binning = RooBinning(len(biased_bins) - 1, biased_bins)
    plot( p, o, cat_data, pdf, components = {'sig*' : dict(LineColor = kGreen, LineStyle = kDashed)}
          , plotResidHist = True
          , dataOpts = dict(MarkerSize = 0.8, MarkerColor = kBlack, Binning = binning)
          , frameOpts = {'Title' : title}
          , pdfOpts  = dict(LineWidth = 2, **pdfOpts)
          , logy = False
          , logx = True
          )
    
# plot the efficiency shapes
__frames = []
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
    __frames.append(frame)
    
shapes = res_model.shapes()
eff_canvas = TCanvas('eff_canvas', 'eff_canvas', 1000, 500)
from ROOT import kYellow, kOrange
for p, shape in zip(eff_canvas.pads(len(shapes), 1), shapes):
    plot_shape(p, t, shape, errorOpts = {'result' : result, 3 : kYellow, 1 : kOrange})

output = {'hlt1_shape' : 'hlt1_excl_biased_dec_excl_biased_bin',
          'hlt2_shape' : 'hlt2_biased_biased_bin'}
output_file = TFile.Open('efficiencies.root', 'recreate')
allVars = w.allVars()
from ROOT import TH1D
for name, pat in output.iteritems():
    n = len(biased_bins)
    heights = [v for v in allVars if v.GetName().find(pat) != -1]
    heights = sorted(heights, key = lambda v: int(v.GetName().split('_', 1)[-1]))
    v = [(h.getVal(), h.getError()) for h in heights]
    hist = TH1D(name, name, n - 1, biased_bins)
    for i in range(1, n):
        hist.SetBinContent(i, v[i - 1][0])
        hist.SetBinError(i, v[i - 1][1])
    output_file.WriteTObject(hist)
output_file.Close()
