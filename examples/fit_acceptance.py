from RooFitWrappers import *
from P2VVLoad import P2VVLibrary

from ROOT import RooMsgService

# RooMsgService.instance().addStream(RooFit.INFO,RooFit.Topic(RooFit.Optimization))

obj = RooObject( workspace = 'w')
w = obj.ws()

from math import pi
t = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax=(0.2, 14))
m = RealVar('mass', Title = 'B mass', Unit = 'MeV', Observable = True, MinMax = (5250, 5550))

# Categories
biased = Category('triggerDecision', States = {'biased' : 1, 'not_biased' : 0})
unbiased = Category('triggerDecisionUnbiased', States = {'unbiased' : 1, 'not_unbiased' : 0})
selected = Category('sel', States = {'Selected' : 1, 'NotSelected' : 0})

observables = [t, m, biased, unbiased, selected]

# Build acceptance
a   = RealVar('a', Title = 'a', Value = 1.45, MinMax = (1, 2))
c   = RealVar('c', Title = 'c', Value = -2.37, MinMax = (-3, -1))
pre = RealVar('eff_pre', Title = 'effective prescale', Value = 0.85, MinMax = (0.5, 0.999), Constant = True)
ub_eff = FormulaVar('pre', '@0', [pre, t])
b_eff  = FormulaVar('det', "(@0 > 0.) ? (1 / (1 + (@1 * @0) ** (@2))) : 0.0001", [t, a, c])

biased_eff   = FormulaVar('biased_eff',   "@0 * (1 - @1)", [b_eff, ub_eff])
unbiased_eff = FormulaVar('unbiased_eff', "@0 * (1 - @1)", [ub_eff, b_eff])
both_eff     = FormulaVar('both_eff',     "@0 * @1",       [ub_eff, b_eff])

from P2VVBinningBuilders import build1DVerticalBinning
binning, eff_func = build1DVerticalBinning('time_binning', b_eff, t, 0.05, 1.)

# now build the actual signal PDF...
from ROOT import RooGaussian as Gaussian
from ROOT import RooExponential as Exponential
from ROOT import RooDecay as Decay

signal_tau = RealVar('signal_tau', Title = 'mean lifetime', Unit = 'ps', Value =  1.5,
                     MinMax = (1., 2.5))

# Time resolution model
from P2VVParameterizations.TimeResolution import LP2011_TimeResolution
tres = LP2011_TimeResolution(time = t)['model']

# Signal time pdf
sig_t = Pdf(Name = 'sig_t', Type = Decay,  Parameters = [t, signal_tau, tres, 'SingleSided'])

# B mass pdf
from P2VVParameterizations.MassPDFs import LP2011_Signal_Mass as Signal_BMass, LP2011_Background_Mass as Background_BMass
sig_m = Signal_BMass(     Name = 'sig_m', mass = m, m_sig_mean = dict( Value = 5365, MinMax = (5363,5372) ) )

# Apply acceptance
from P2VVGeneralUtils import readData
tree_name = 'DecayTree'
## input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhiPrescaled_ntupleB_for_fitting_20120110.root'
input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhi_ntupleB_for_fitting_20120118.root'
data = readData(input_file, tree_name, cuts = '(sel == 1)',
                NTuple = True, observables = observables)

# Build the acceptance using the histogram as starting values
input_file = '/stuff/PhD/p2vv/data/efficiencies.root'
histogram = 'signal_efficiency_histo_20bins'

from ROOT import TFile
acceptance_file = TFile.Open(input_file)
if not acceptance_file:
    raise ValueError, "Cannot open histogram file %s" % input_file
histogram = acceptance_file.Get(histogram)
if not histogram:
    raise ValueError, 'Cannot get acceptance historgram %s from file' % histogram
xaxis = histogram.GetXaxis()

from array import array
biased_bins = array('d', (xaxis.GetBinLowEdge(i) for i in range(1, histogram.GetNbinsX() + 2)))
biased_heights = [histogram.GetBinContent(i) for i in range(1, histogram.GetNbinsX() + 1)]

unbiased_bins = array('d', [0, 10])
unbiased_heights = [0.2]

# Spec to build efficiency shapes
spec = {"Bins" : {biased : {'state'   : 'biased',
                            'bounds'  : biased_bins,
                            'heights' : biased_heights},
                  unbiased : {'state' : 'unbiased',
                              'bounds' : unbiased_bins,
                              'heights' : unbiased_heights}
                  },
        "Relative" : {((biased, "biased"),     (unbiased, "unbiased")) : {'Value' : 0.2, 'MinMax' : (0.1, 0.45)},
                      ((biased, "not_biased"), (unbiased, "unbiased")) : {'Value' : 0.3, 'MinMax' : (0.1, 0.45)},
                      ((biased, "biased"),     (unbiased, "not_unbiased")) : None}
        }
mhe = MultiHistEfficiency(Name = "RMHE", Original = sig_t, Observable = t,
                          ConditionalCategories = True, **spec)

entries = mhe.getEntries()
superCat = mhe.getSuper()
data.addColumn(superCat)

states = {}
for entry in entries:
    label = superCat.lookupType(entry.first).GetName()
    states[label] = entry.second.effProd()

sig_t = SimultaneousPdf(Name = 'signal_time', States = states, SplitCategory = superCat)

# Create signal component
signal = Component('signal', (sig_m.pdf(), sig_t), Yield = (30000,100,100000))

# Create combinatorical background component
bkg_m = Background_BMass( Name = 'bkg_m', mass = m, m_bkg_exp  = dict( Name = 'm_bkg_exp' ) )

bkg_tau = RealVar('bkg_tau', Title = 'comb background lifetime', Unit = 'ps', Value = 1, MinMax = (0.0001, 5))

from P2VVParameterizations.TimePDFs import LP2011_Background_Time as Background_Time
bkg_t = Background_Time( Name = 'bkg_t', time = t, resolutionModel = tres
                       , t_bkg_fll    = dict( Name = 't_bkg_fll',    Value = 0.3 )
                       , t_bkg_ll_tau = dict( Name = 't_bkg_ll_tau', Value = 1.92, MinMax = (0.5,2.5) )
                       , t_bkg_ml_tau = dict( Name = 't_bkg_ml_tau', Value = 0.21, MinMax = (0.01,0.5) ) )
background = Component('background', (bkg_m.pdf(), bkg_t.pdf()), Yield = (50000,100,300000) )


pdf = buildPdf(Components = (signal, background), Observables = (m, t), Name='pdf')
pdf.Print("t")

## from Helpers import Mapping
## mapping = Mapping({m : 'm', t : 'tau'}, data)

## Fit
print 'fitting data'
## from profiler import profiler_start, profiler_stop
## profiler_start("acceptance.log")
pdf.fitTo(data, NumCPU = 4 , Timer = 1, Minimizer = 'Minuit2', Verbose = True, Optimize = 1)
## profiler_stop()

from ROOT import kDashed, kRed, kGreen, kBlue, kBlack
from ROOT import TCanvas
canvas = {}
print 'plotting'
for cat in split_cat:
    idx, label = cat.getVal(), cat.GetName()
    if label == 'NotMapped':
        continue
    name = split_cat.GetName() + '_' + label
    canv = canvas[name] = TCanvas(name, name, 1000, 500)
    obs =  [o for o in pdf.Observables() if hasattr(o,'frame')]
    for (p,o) in zip(canv.pads(len(obs)), obs):
        dataOpts = dict(Cut = '{0} == {0}::{1}'.format(split_cat['Name'], label) )
        pdfOpts  = dict(Slice = (split_cat, label), ProjWData = (RooArgSet(split_cat), data))
        from P2VVGeneralUtils import plot
        plot( p, o, data, pdf, components = { 'sig*' : dict(LineColor = kGreen, LineStyle = kDashed)
                                            , 'bkg*' : dict(LineColor = kBlue,  LineStyle = kDashed)
                                              }
              , dataOpts = dict( MarkerSize = 0.8, MarkerColor = kBlack, **dataOpts )
              , pdfOpts  = dict( LineWidth = 2, **pdfOpts )
              , logy = ( o == t )
              )
