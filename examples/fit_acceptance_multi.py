from itertools import product
def valid_combinations(states):
    all_states = []
    for level in states:
        all_states.extend(level.keys())
    labels = [[(state, label.GetName()) for label in state] for state in all_states]
    all_combinations = list(product(*labels))
    valid = []
    def good(combination):
        s = set(combination)
        for level in states:
            level_good = False
            for entry in level.iteritems():
                if entry in s:
                    level_good = True
                    break
            if not level_good:
                return level_good
        return True
    return filter(good, all_combinations)

from RooFitWrappers import *
from P2VVLoad import P2VVLibrary
from ROOT import RooCBShape as CrystalBall

from ROOT import RooMsgService

# RooMsgService.instance().addStream(RooFit.INFO,RooFit.Topic(RooFit.Optimization))

obj = RooObject( workspace = 'w')
w = obj.ws()

from math import pi
t = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax=(0.3, 14))
m = RealVar('mass', Title = 'B mass', Unit = 'MeV', Observable = True, MinMax = (5250, 5550))
mpsi = RealVar('mdau1', Title = 'J/psi mass', Unit = 'MeV', Observable = True, MinMax = (3030, 3150))
st = RealVar('sigmat',Title = '#sigma(t)', Unit = 'ps', Observable = True, MinMax = (0.0001, 0.12))

# Categories
hlt1_biased = Category('hlt1_biased', States = {'biased' : 1, 'not_biased' : 0}, Observable = True)
hlt1_unbiased = Category('hlt1_unbiased', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
hlt2_biased = Category('hlt2_biased', States = {'biased' : 1, 'not_biased' : 0}, Observable = True)
hlt2_unbiased = Category('hlt2_unbiased', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
project_vars = [hlt1_biased, hlt1_unbiased, hlt2_biased, hlt2_unbiased, st]

selected = Category('sel', States = {'Selected' : 1, 'NotSelected' : 0})

observables = [t, m, mpsi, st, hlt1_biased, hlt1_unbiased, hlt2_biased, hlt2_unbiased, selected]

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
                                                    Constant = False))
## from P2VVParameterizations.TimeResolution import Gaussian_TimeResolution as TimeResolution
## tres = TimeResolution(time = t).model()

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
mass_pdf = buildPdf(Components = (signal, psi_background, background), Observables = (m, mpsi), Name='mass_pdf')
mass_pdf.Print("t")

# Build the acceptance using the histogram as starting values
input_file = '/home/raaij/data/start_values.root'
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


valid = valid_combinations([{hlt1_biased : 'biased', hlt1_unbiased : 'unbiased'}, {hlt2_biased : 'biased', hlt2_unbiased : 'unbiased'}])

# Spec to build efficiency shapes
spec = {"Bins" : {hlt1_biased : {'state'   : 'biased',
                                 'bounds'  : biased_bins,
                                 'heights' : hlt1_biased_heights},
                  hlt1_unbiased : {'state' : 'unbiased',
                                   'bounds' : unbiased_bins,
                                   'heights' : hlt1_unbiased_heights},
                  hlt2_biased : {'state'   : 'biased',
                                 'bounds'  : biased_bins,
                                 'heights' : hlt2_biased_heights},
                  hlt2_unbiased : {'state' : 'unbiased',
                                   'bounds' : unbiased_bins,
                                   'heights' : hlt2_unbiased_heights}
                  }
        }
# Read input data
from P2VVGeneralUtils import readData
tree_name = 'DecayTree'
## input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhiPrescaled_ntupleB_for_fitting_20120110.root'
input_file = '/home/raaij/data/Bs2JpsiPhi_2011_biased_unbiased.root'

## Fit options
fitOpts = dict(NumCPU = 12, Timer = 1, Save = True, Verbose = True, Optimize = 1, Minimizer = 'Minuit2')

data = None
real_data = True
if real_data:
    data = readData(input_file, tree_name, cuts = 'sel == 1 && (hlt1_biased == 1 || hlt1_unbiased == 1) && (hlt2_biased == 1 || hlt2_unbiased == 1)',
                    NTuple = True, observables = observables)
    total = data.sumEntries()

    rel_spec = {}
    for comb in valid:
        cuts = ' && '.join(['{0} == {0}::{1}'.format(state.GetName(), label) for state, label in comb])
        rel_spec[comb] = {'Value' : data.sumEntries(cuts) / total, "Constant" : True}

    spec['Relative'] = rel_spec
    pdf = MultiHistEfficiency(Name = "RMHE", Original = sig_t.pdf(), Observable = t,
                              ConditionalCategories = True, **spec)
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
else:
    rel_spec = {}
    for comb in valid:
        cuts = ' && '.join(['{0} == {0}::{1}'.format(state.GetName(), label) for state, label in comb])
        rel_spec[comb] = {'Value' : 1. / len(valid), "Constant" : True},
    pdf = MultiHistEfficiency(Name = "RMHE", Original = sig_t.pdf(), Observable = t,
                              ConditionalCategories = True, **spec)
    pdf.Print('v')        
    data = pdf.generate([t, hlt1_biased, hlt1_unbiased, hlt2_unbiased, hlt2_biased], 25000)

# Get the SuperCategory from the MultiHistEfficiency and add it to the data.
entries = pdf.getEntries()
super_cat = pdf.getSuper()

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

# Plot the lifetime shapes
canv = TCanvas('canvas', 'canvas', 900, 900)
obs = [t]
for states, (p, o) in zip(spec['Relative'].keys(), (i for i in product(canv.pads(3, 3), obs))):
    name = '__'.join(['%s_%s' % (state.GetName(), label) for state, label in states])
    cuts = ' && '.join(['{0} == {0}::{1}'.format(state.GetName(), label) for state, label in states])
    cat_data = data.reduce(cuts)
    pdfOpts = dict(ProjWData = (RooArgSet(*project_vars), cat_data))
    from P2VVGeneralUtils import plot
    plot( p, o, cat_data, pdf, components = { 'sig*' : dict(LineColor = kGreen, LineStyle = kDashed)
                                              ## , 'bkg*' : dict(LineColor = kBlue,  LineStyle = kDashed)
                                              }
          , dataOpts = dict( MarkerSize = 0.8, MarkerColor = kBlack,
                             Binning = RooBinning(len(biased_bins) - 1, biased_bins))
          , pdfOpts  = dict( LineWidth = 2, **pdfOpts )
          , logy = False
          )

# plot the efficiency shapes
