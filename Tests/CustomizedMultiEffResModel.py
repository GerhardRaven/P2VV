import os
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
hlt1_unbiased = Category('hlt1_unbiased_dec', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
hlt1_excl_biased = Category('hlt1_excl_biased', States = {'excl_biased' : 1, 'unbiased' : 0}, Observable = True)
hlt2_biased = Category('hlt2_biased', States = {'biased' : 1, 'not_biased' : 0}, Observable = True)
hlt2_unbiased = Category('hlt2_unbiased', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
hlt2_excl_biased = Category('hlt2_excl_biased', States = {'excl_biased' : 1, 'unbiased' : 0}, Observable = True)

## project_vars = [hlt1_biased, hlt1_unbiased, hlt2_biased, hlt2_unbiased, st]
categories = [hlt1_biased, hlt1_unbiased, hlt1_excl_biased,
              hlt2_biased, hlt2_unbiased, hlt2_excl_biased]
categories = dict([(c.GetName(), c) for c in categories])

project_vars = [hlt1_excl_biased, hlt1_unbiased, hlt2_biased, hlt2_unbiased, st]

selected = Category('sel', States = {'Selected' : 1, 'NotSelected' : 0})

observables = [t, m, mpsi, st, hlt1_biased, hlt1_unbiased, hlt1_excl_biased,
               hlt2_biased, hlt2_unbiased, hlt2_excl_biased, selected, nPV]

# now build the actual signal PDF...
from ROOT import RooGaussian as Gaussian
from ROOT import RooExponential as Exponential
from ROOT import RooDecay as Decay

signal_tau = RealVar('signal_tau', Title = 'mean lifetime', Unit = 'ps', Value =  1.5,
                     MinMax = (1., 2.5))

# Time resolution model
from P2VVParameterizations.TimeResolution import Gaussian_TimeResolution as TimeResolution
tres = TimeResolution(Name = 'tres', time = t, sigmat = st, PerEventError = True,
                      BiasScaleFactor = False, Cache = True,
                      bias = dict(Value = 0, MinMax = (-1, 1), Constant = True),
                      sigmaSF  = dict(Value = 1.46, MinMax = (0.1, 5), Constant = True))

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
## input_file = os.path.join(base_location, 'data/start_values.root')
## hlt1_histogram = 'hlt1_shape'
## hlt2_histogram = 'hlt2_shape'

input_file = os.path.join(base_location, 'data/Bs_HltPropertimeAcceptance_Data-20120816.root')
##hlt1_histogram = 'hlt1_shape'
hlt2_histogram = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_Data_20bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'

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

## hlt1_histogram = acceptance_file.Get(hlt1_histogram)
## if not hlt1_histogram:
##     raise ValueError, 'Cannot get acceptance historgram %s from file' % hlt1_histogram


from array import array
biased_bins = array('d', (xaxis.GetBinLowEdge(i) for i in range(1, hlt2_histogram.GetNbinsX() + 2, 2)))
unbiased_bins = array('d', [0.3, 14])

## hlt1_biased_heights = [hlt1_histogram.GetBinContent(i) for i in range(1, hlt1_histogram.GetNbinsX() + 1)]
hlt1_biased_average = (6.285e-01 - 0.25, 1.633e-02)
hlt1_biased_heights = [hlt1_biased_average[0] + i / 2 * 0.05 for i in range(1, hlt2_histogram.GetNbinsX() + 1, 2)]

hlt1_unbiased_heights = [0.5]

hlt2_biased_average = (6.33e-01, 1.65e-02)
hlt2_biased_heights = [hlt2_histogram.GetBinContent(i) for i in range(1, hlt2_histogram.GetNbinsX() + 1, 2)]
av = sum(hlt2_biased_heights) / float(len(hlt2_biased_heights))
scale = av / hlt2_biased_average[0]
hlt2_biased_heights = [h / scale for h in hlt2_biased_heights]

hlt2_unbiased_heights = [0.5]

## valid_definition = [[(hlt1_biased, 'biased'), (hlt1_unbiased, 'unbiased')], [(hlt2_biased, 'biased'), (hlt2_unbiased, 'unbiased')]]
## valid_definition = [[(hlt2_biased, 'biased'), (hlt2_unbiased, 'unbiased')]]
valid_definition = [[(hlt1_excl_biased, 'excl_biased'), (hlt1_excl_biased, 'unbiased')], [(hlt2_biased, 'biased'), (hlt2_unbiased, 'unbiased')]]
## valid_definition = [[(hlt1_excl_biased, 'excl_biased'), (hlt1_excl_biased, 'unbiased')], [(hlt2_excl_biased, 'excl_biased'), (hlt2_excl_biased, 'unbiased')]]
valid = valid_combinations(valid_definition)

spec = {'Bins' : {hlt1_excl_biased : {'excl_biased' : {'bins'    : biased_bins,
                                                       'heights' : hlt1_biased_heights,
                                                       'average' : hlt1_biased_average},
                                      'unbiased'    : {'bins'    : unbiased_bins,
                                                       'heights' : hlt1_unbiased_heights}
                                      }
                  }
        }

## spec = {'Bins' : {hlt1_excl_biased : {'excl_biased' : {'bins'    : biased_bins,
##                                                        'heights' : hlt1_biased_heights,
##                                                        'average' : hlt1_biased_average},
##                                       'unbiased'    : {'bins'    : unbiased_bins,
##                                                        'heights' : hlt1_unbiased_heights}
##                                       },
##                   hlt2_biased      : {'biased'      : {'bins'    : biased_bins,
##                                                        'heights' : hlt2_biased_heights,
##                                                        'average' : hlt2_biased_average}
##                                       },
##                   hlt2_unbiased    : {'unbiased'    : {'bins'    : unbiased_bins,
##                                                        'heights' : hlt2_unbiased_heights}
##                                       }
##                   }
##         }

## spec = {'Bins' : {hlt2_biased      : {'biased'      : {'bins'    : biased_bins,
##                                                        'heights' : hlt2_biased_heights,
##                                                        'average' : (6.330e-01, 1.65e-02)}
##                                       },
##                   hlt2_unbiased    : {'unbiased'    : {'bins'    : unbiased_bins,
##                                                        'heights' : hlt2_unbiased_heights}
##                                       }
##                   }
##         }

## spec = {'Bins' : {hlt1_excl_biased : {'excl_biased' : {'bins'    : biased_bins,
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
from P2VVGeneralUtils import readData
tree_name = 'DecayTree'
## input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhiPrescaled_ntupleB_for_fitting_20120110.root'

## Fit options
fitOpts = dict(NumCPU = 4, Timer = 1, Save = True, Verbose = True, Optimize = 1,
               Strategy = 2, Minimizer = 'Minuit2')

rel_spec = {(('hlt1_excl_biased', 'excl_biased'),) : 0.5,
            (('hlt1_excl_biased', 'unbiased'),)    : None}
spec['Relative'] = dict([(tuple((categories[c], l) for c, l in k), {'Constant' : True, 'Value' : v} if v else None) for k, v in rel_spec.iteritems()])
for comb in valid:
    cuts = ' && '.join(['{0} == {0}::{1}'.format(state.GetName(), label) for state, label in comb])
    rel_spec[comb] = {'Value' : 1. / len(valid), "Constant" : True},

res_model = MultiHistEfficiencyModel(Name = "RMHE", Original = sig_t.pdf(), Observable = t,
                                     ConditionalCategories = True, UseSingleBinConstraint = False,
                                     ResolutionModel = tres.model(), Fit = False, **spec)
pdf = Single_Exponent_Time(Name = 'pdf', time = t, resolutionModel = res_model)
pdf = pdf.pdf()

## data = pdf.generate([t, hlt1_excl_biased, hlt2_unbiased, hlt2_biased], 30000)

from RooFitWrappers import RealVar, Customizer
newTime = RealVar('time_new', Title = 'Decay time', Unit = 'ps', Observable = True, Value = 0.5, MinMax = (0.3, 14.))

custPdf = Customizer(Pdf = pdf, OriginalArgs = [t], SubstituteArgs = [newTime], NameSuffix = 'newTime')

from ROOT import RooArgSet
intSet      = RooArgSet()
normSetOrig = RooArgSet(t)
normSetCust = RooArgSet(newTime)

print '\noriginal PDF:'
pdf.Print('t')

print '\ncustomized PDF:'
custPdf.Print('t')

print '\noriginal PDF integral:'
pdfInt  = pdf.createIntegral(normSetOrig)
## pdfInt.Print('t')
print pdfInt.getVal()

print '\ncustomized PDF integral:'
custInt = custPdf.createIntegral(normSetCust)
## custInt.Print('t')
print custInt.getVal()
